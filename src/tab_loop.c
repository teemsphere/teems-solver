#include <teems_solver.h>

/* Loops in TAB files (manual 11.18) and their LHS-mapping use (11.9.8.1).

   tab_loop_transform rewrites the preprocessed TAB (one statement per
   line). Each loop becomes a pair of marker statements,
     loop (begin) <id> <index> <set> <name> <fused>
     loop (end) <id>
   and its body statements are made ordinary statements over the loop:

   - general loop: every body Formula and Assertion gets the quantifier
     (all,<index>,lp@<id>) of every enclosing loop, outermost first.
     lp@<id> is a one-element subset of the loop set, created after the
     set read (tab_loop_sets_link) and pointed at the current element
     each iteration (loop_set_point), so the body sees the loop index
     through the ordinary subset routing. statements_execute runs the
     body in file order once per element.
   - fused loop: a body of exactly one Formula with no quantifier of
     its own (the 11.9.8.1 shape, "Loop (all,g,GRD); Formula
     WGT2(M1(g),M2(g)) = WGT2(M1(g),M2(g)) + WGT(g); Loop (end)") runs
     as that Formula over (all,<index>,<set>) with the (loopserial)
     qualifier: serial in loop order, an element assigned more than
     once allowed (the order makes it unambiguous), no loop driver.

   BREAK and CYCLE become probe assertions,
     assertion (loopctl) (lc=<b|c>,<e|a>,<target id>) <pins> <rest>
   which assertions_execute evaluates only in probe mode and counts
   rather than reports: EVERY fires when the condition holds for every
   tuple, ANY when it holds for at least one. */

teems_loop_syn *teems_loop_syns = NULL;
int teems_loop_nsyn = 0;
static int loop_next_id = 1;

typedef struct {
  int id, fused, depth, parent, begin, end;
  char idx[NAMESIZE], set[NAMESIZE], name[NAMESIZE];
} lt_loop;

static int lt_kw(const char *s, char *kw, int cap) {
  int k = 0;
  while (*s == ' ') s++;
  while (*s != '\0' && *s != ' ' && *s != '(' && *s != ';' && k < cap - 1) kw[k++] = *s++;
  kw[k] = '\0';
  return k;
}

/* end of a balanced group starting at s ('(' at s[0]) */
static const char *lt_group_end(const char *s) {
  int d = 0;
  for (; *s != '\0'; s++) {
    if (*s == '(') d++;
    else if (*s == ')' && --d == 0) return s + 1;
  }
  return NULL;
}

/* offset in a statement after its keyword, qualifiers and # label # --
   where the quantifiers start */
static size_t lt_quant_at(const char *s) {
  const char *p = s;
  char kw[NAMESIZE];
  while (*p == ' ') p++;
  p += lt_kw(p, kw, sizeof(kw));
  for (;;) {
    while (*p == ' ') p++;
    if (*p == '(' && strncmp(p, "(all,", 5) != 0) {
      const char *e = lt_group_end(p);
      if (e == NULL) break;
      p = e;
      continue;
    }
    if (*p == '#') {
      const char *e = strchr(p + 1, '#');
      if (e == NULL) break;
      p = e + 1;
      continue;
    }
    break;
  }
  return (size_t)(p - s);
}

static int lt_insert(char *line, size_t cap, size_t at, const char *ins) {
  size_t n = strlen(line), m = strlen(ins);
  if (n + m + 2 >= cap) return -1;
  memmove(line + at + m, line + at, n - at + 1);
  memcpy(line + at, ins, m);
  return 0;
}

static void lt_trim(char *s) {
  size_t n = strlen(s);
  while (n > 0 && (s[n - 1] == '\n' || s[n - 1] == '\r' || s[n - 1] == ' ' || s[n - 1] == ';')) s[--n] = '\0';
}

/* "(begin, name=x)" / "(end)" / "(every, name=x)" -> flags */
static int lt_qual(const char *g, int *begin, int *end, int *any, int *every, char *name) {
  char buf[TABREADLINE], *t, *sv = NULL;
  size_t n = strlen(g);
  if (n < 2 || n >= sizeof(buf)) return -1;
  memcpy(buf, g + 1, n - 2);
  buf[n - 2] = '\0';
  for (t = strtok_r(buf, ",", &sv); t != NULL; t = strtok_r(NULL, ",", &sv)) {
    while (*t == ' ') t++;
    lt_trim(t);
    if (strcmp(t, "begin") == 0) *begin = 1;
    else if (strcmp(t, "end") == 0) *end = 1;
    else if (strcmp(t, "any") == 0) *any = 1;
    else if (strcmp(t, "every") == 0) *every = 1;
    else if (strncmp(t, "name=", 5) == 0 || strncmp(t, "name =", 6) == 0) {
      const char *v = strchr(t, '=') + 1;
      while (*v == ' ') v++;
      if (strlen(v) == 0 || strlen(v) >= NAMESIZE) return -1;
      strcpy(name, v);
    }
    else return -1;
  }
  return 0;
}

/* does the statement use the index as a name (not inside quotes)? A
   body statement that does not needs no quantifier over the loop: the
   driver runs it once per element anyway */
static int lt_uses(const char *s, const char *idx) {
  size_t n = strlen(idx);
  int inq = 0;
  const char *p;
  for (p = s; *p != '\0'; p++) {
    if (*p == '"') { inq = !inq; continue; }
    if (inq) continue;
    if (strncmp(p, idx, n) == 0) {
      char a = (p == s) ? ' ' : p[-1], b = p[n];
      if (!(isalnum((unsigned char)a) || a == '_' || a == '@') && !(isalnum((unsigned char)b) || b == '_' || b == '@')) return 1;
    }
  }
  return 0;
}

static const char *const lt_allowed[] = {"formula", "assertion", "break", "cycle", "loop", "zerodivide", NULL};

int tab_loop_transform(char *fname) {
  FILE *f, *fo;
  char (*ln)[TABREADLINE] = NULL;
  int nl = 0, cap = 0, i, j, nloop = 0, top = -1, rc = 0;
  int *stack = NULL;
  lt_loop *lp = NULL;
  char tmpname[TABREADLINE], kw[NAMESIZE];
  f = teems_fopen(fname, "r");
  if (f == NULL) return 0;
  {
    char buf[TABREADLINE];
    int any = 0;
    while (fgets(buf, sizeof(buf), f)) {
      if (nl == cap) {
        cap = cap ? 2 * cap : 256;
        ln = realloc(ln, (size_t)cap * TABREADLINE);
      }
      strcpy(ln[nl], buf);
      lt_kw(buf, kw, sizeof(kw));
      if (strcmp(kw, "loop") == 0 || strcmp(kw, "break") == 0 || strcmp(kw, "cycle") == 0) any = 1;
      nl++;
    }
    fclose(f);
    if (!any) {
      free(ln);
      return 0;
    }
  }
  lp = calloc(nl + 1, sizeof(lt_loop));
  stack = calloc(nl + 1, sizeof(int));
  /* structure: begin/end pairs, the loop quantifier, names */
  for (i = 0; i < nl; i++) {
    char *s = ln[i], *p;
    int b = 0, e = 0, an = 0, ev = 0;
    char nm[NAMESIZE] = "";
    lt_kw(s, kw, sizeof(kw));
    if (strcmp(kw, "loop") != 0) {
      if (top < 0 && (strcmp(kw, "break") == 0 || strcmp(kw, "cycle") == 0)) {
        lt_trim(s);
        errmsg("Error: %s outside any loop (manual 11.18): %.120s\n", kw[0] == 'b' ? "BREAK" : "CYCLE", s);
        rc = -1;
        goto done;
      }
      continue;
    }
    p = s;
    while (*p == ' ') p++;
    p += 4;
    while (*p == ' ') p++;
    if (*p != '(' || lt_group_end(p) == NULL) {
      lt_trim(s);
      errmsg("Error: a Loop statement needs (BEGIN) or (END) (manual 11.18): %.120s\n", s);
      rc = -1;
      goto done;
    }
    {
      char g[TABREADLINE];
      const char *ge = lt_group_end(p);
      size_t gl = (size_t)(ge - p);
      if (gl >= sizeof(g)) gl = sizeof(g) - 1;
      memcpy(g, p, gl);
      g[gl] = '\0';
      if (lt_qual(g, &b, &e, &an, &ev, nm) < 0 || an || ev || b == e) {
        lt_trim(s);
        errmsg("Error: malformed Loop qualifier %s; Loop (BEGIN[, name=LoopName]) (all,index,set) or Loop (END) (manual 11.18): %.120s\n", g, s);
        rc = -1;
        goto done;
      }
      p = (char *)ge;
    }
    if (b) {
      lt_loop *L = &lp[nloop];
      char q[TABREADLINE], *a, *c;
      while (*p == ' ') p++;
      strcpy(q, p);
      lt_trim(q);
      if (strncmp(q, "(all,", 5) != 0 || lt_group_end(q) == NULL || *lt_group_end(q) != '\0') {
        errmsg("Error: Loop (BEGIN) takes exactly one quantifier (all,index,set) (manual 11.18): %.120s\n", s);
        rc = -1;
        goto done;
      }
      a = q + 5;
      c = strchr(a, ',');
      if (c == NULL || strchr(a, ':') != NULL) {
        errmsg("Error: Loop (BEGIN) takes a plain quantifier (all,index,set), without a condition (manual 11.18): %.120s\n", s);
        rc = -1;
        goto done;
      }
      *c = '\0';
      c[strlen(c + 1)] = '\0';
      if (strlen(a) == 0 || strlen(a) >= NAMESIZE || strlen(c + 1) == 0 || strlen(c + 1) >= NAMESIZE) {
        errmsg("Error: malformed Loop quantifier: %.120s\n", s);
        rc = -1;
        goto done;
      }
      strcpy(L->idx, a);
      strcpy(L->set, c + 1);
      {
        size_t k = strlen(L->set);
        if (k > 0 && L->set[k - 1] == ')') L->set[k - 1] = '\0';
      }
      strcpy(L->name, nm[0] ? nm : "-");
      L->id = loop_next_id++;
      L->begin = i;
      L->end = -1;
      L->parent = top >= 0 ? stack[top] : -1;
      L->depth = top + 1;
      for (j = 0; j <= top; j++) {
        if (strcmp(lp[stack[j]].idx, L->idx) == 0) {
          errmsg("Error: nested loops reuse the index %s (manual 11.18): %.120s\n", L->idx, s);
          rc = -1;
          goto done;
        }
        if (nm[0] && strcmp(lp[stack[j]].name, nm) == 0) {
          errmsg("Error: nested loops share the name %s (manual 11.18)\n", nm);
          rc = -1;
          goto done;
        }
      }
      stack[++top] = nloop++;
    }
    else {
      if (top < 0) {
        errmsg("Error: Loop (END) without a matching Loop (BEGIN) (manual 11.18)\n");
        rc = -1;
        goto done;
      }
      lp[stack[top--]].end = i;
    }
  }
  if (top >= 0) {
    errmsg("Error: Loop (BEGIN) over %s has no matching Loop (END) (manual 11.18)\n", lp[stack[top]].set);
    rc = -1;
    goto done;
  }
  /* fused: one Formula, no quantifier of its own, nothing else */
  for (j = 0; j < nloop; j++) {
    lt_loop *L = &lp[j];
    if (L->end == L->begin + 2) {
      char *s = ln[L->begin + 1];
      lt_kw(s, kw, sizeof(kw));
      if (strcmp(kw, "formula") == 0 && strncmp(s + lt_quant_at(s), "(all,", 5) != 0 && lt_uses(s, L->idx)) L->fused = 1;
    }
  }
  /* rewrite */
  top = -1;
  for (i = 0; i < nl && rc == 0; i++) {
    char *s = ln[i];
    int k;
    lt_kw(s, kw, sizeof(kw));
    if (strcmp(kw, "loop") == 0) {
      for (j = 0; j < nloop; j++) if (lp[j].begin == i || lp[j].end == i) break;
      if (lp[j].begin == i) {
        snprintf(s, TABREADLINE, "loop (begin) %d %s %s %s %d;\n", lp[j].id, lp[j].idx, lp[j].set, lp[j].name, lp[j].fused);
        stack[++top] = j;
        if (!lp[j].fused) {
          teems_loop_syns = realloc(teems_loop_syns, (size_t)(teems_loop_nsyn + 1) * sizeof(teems_loop_syn));
          teems_loop_syns[teems_loop_nsyn].id = lp[j].id;
          strcpy(teems_loop_syns[teems_loop_nsyn].parent, lp[j].set);
          teems_loop_nsyn++;
        }
      }
      else {
        snprintf(s, TABREADLINE, "loop (end) %d;\n", lp[j].id);
        top--;
      }
      continue;
    }
    if (top < 0) continue;
    for (k = 0; lt_allowed[k] != NULL; k++) if (strcmp(kw, lt_allowed[k]) == 0) break;
    if (lt_allowed[k] == NULL) {
      lt_trim(s);
      if (strcmp(kw, "write") == 0)
        errmsg("Error: Write inside a loop is not supported (TEEMS writes once, after the pass); move it after Loop (END) (manual 11.18): %.120s\n", s);
      else
        errmsg("Error: a %s statement is not allowed inside a loop; only Formula, Assertion, Break, Cycle and ZeroDivide are (manual 11.18): %.120s\n", kw, s);
      rc = -1;
      break;
    }
    if (strcmp(kw, "zerodivide") == 0) continue;
    {
      char pins[TABREADLINE] = "";
      int fusedhere = lp[stack[top]].fused;
      for (j = 0; j <= top; j++) {
        lt_loop *L = &lp[stack[j]];
        char q[3 * NAMESIZE];
        if (!lt_uses(s, L->idx)) continue;
        if (L->fused) snprintf(q, sizeof(q), "(loopserial) (all,%s,%s) ", L->idx, L->set);
        else snprintf(q, sizeof(q), "(all,%s,lp@%d) ", L->idx, L->id);
        if (strlen(pins) + strlen(q) + 1 >= sizeof(pins)) {
          errmsg("Error: loops nested too deeply\n");
          rc = -1;
          break;
        }
        strcat(pins, q);
      }
      if (rc) break;
      if (fusedhere && strcmp(kw, "formula") != 0) {
        rc = -1;
        break;
      }
      if (strcmp(kw, "break") == 0 || strcmp(kw, "cycle") == 0) {
        char rest[TABREADLINE], nm[NAMESIZE] = "", ctl[3 * NAMESIZE];
        char *p = s;
        int b = 0, e = 0, an = 0, ev = 0, tgt = stack[top];
        while (*p == ' ') p++;
        p += strlen(kw);
        while (*p == ' ') p++;
        if (*p == '(' && strncmp(p, "(all,", 5) != 0) {
          char g[TABREADLINE];
          const char *ge = lt_group_end(p);
          size_t gl;
          if (ge == NULL) {
            lt_trim(s);
            errmsg("Error: malformed %s qualifier: %.120s\n", kw, s);
            rc = -1;
            break;
          }
          gl = (size_t)(ge - p);
          if (gl >= sizeof(g)) gl = sizeof(g) - 1;
          memcpy(g, p, gl);
          g[gl] = '\0';
          if (lt_qual(g, &b, &e, &an, &ev, nm) < 0 || b || e || (an && ev)) {
            lt_trim(s);
            errmsg("Error: malformed %s qualifier %s; (EVERY|ANY[, name=LoopName]) (manual 11.18): %.120s\n", kw, g, s);
            rc = -1;
            break;
          }
          p = (char *)ge;
          while (*p == ' ') p++;
        }
        if (nm[0]) {
          for (j = top; j >= 0; j--) if (strcmp(lp[stack[j]].name, nm) == 0) break;
          if (j < 0) {
            lt_trim(s);
            errmsg("Error: %s names loop %s, which does not enclose it (manual 11.18): %.120s\n", kw, nm, s);
            rc = -1;
            break;
          }
          tgt = stack[j];
        }
        strcpy(rest, p);
        lt_trim(rest);
        if (rest[0] == '\0') {
          lt_trim(s);
          errmsg("Error: %s needs a condition (manual 11.18): %.120s\n", kw, s);
          rc = -1;
          break;
        }
        snprintf(ctl, sizeof(ctl), "assertion (loopctl) (lc=%c,%c,%d) ", kw[0], an ? 'a' : 'e', lp[tgt].id);
        if (strlen(ctl) + strlen(pins) + strlen(rest) + 3 >= TABREADLINE) {
          errmsg("Error: %s statement too long after the loop rewrite\n", kw);
          rc = -1;
          break;
        }
        strcpy(s, ctl);
        strcat(s, pins);
        strcat(s, rest);
        strcat(s, ";\n");
        continue;
      }
      if (lt_insert(s, TABREADLINE, lt_quant_at(s), pins) < 0) {
        errmsg("Error: statement too long after the loop rewrite: %.120s\n", s);
        rc = -1;
        break;
      }
    }
  }
  if (rc == 0) {
    strcpy(tmpname, fname);
    strcat(tmpname, "_lp");
    fo = teems_fopen(tmpname, "w");
    if (fo == NULL) {
      errmsg("Error: cannot write %s\n", tmpname);
      rc = -1;
      goto done;
    }
    for (i = 0; i < nl; i++) fputs(ln[i], fo);
    fclose(fo);
    if (rename(tmpname, fname) != 0) {
      errmsg("Error: cannot rename %s\n", tmpname);
      rc = -1;
    }
    else rc = nloop;
  }
done:
  free(ln);
  free(lp);
  free(stack);
  return rc;
}

/* after the set read: point each loop set at its parent's first
   element and register it as a subset of the parent and of every
   superset of the parent */
int tab_loop_sets_link(set_element *se, set_def *sets, dim_t nset) {
  int k;
  dim_t s, par, t;
  for (k = 0; k < teems_loop_nsyn; k++) {
    char nm[NAMESIZE];
    snprintf(nm, sizeof(nm), "lp@%d", teems_loop_syns[k].id);
    for (s = 0; s < nset; s++) if (strcmp(sets[s].setname, nm) == 0) break;
    for (par = 0; par < nset; par++) if (strcmp(sets[par].setname, teems_loop_syns[k].parent) == 0) break;
    if (s == nset) {
      errmsg("Error: loop set %s missing (internal)\n", nm);
      return -1;
    }
    if (par == nset) {
      errmsg("Error: Loop over %s, which is not a declared set (manual 11.18)\n", teems_loop_syns[k].parent);
      return -1;
    }
    if (sets[par].size == 0) sets[s].size = 0;
    teems_loop_syns[k].setid = s;
    teems_loop_syns[k].parentid = par;
    sets[s].subsetid[1] = par;
    for (t = 1; t < MAXSUPSET - 1; t++) {
      if (sets[par].subsetid[t] == -1) break;
      sets[s].subsetid[t + 1] = sets[par].subsetid[t];
    }
    sets[s].intertemp = sets[par].intertemp;
    if (sets[s].size > 0) loop_set_point(se, sets, k, 0);
  }
  return 0;
}

/* point loop set k at element e of its parent */
void loop_set_point(set_element *se, set_def *sets, int k, dim_t e) {
  dim_t s = teems_loop_syns[k].setid, par = teems_loop_syns[k].parentid, t, u;
  set_element *x = &se[sets[s].offset], *pe = &se[sets[par].offset + e];
  strcpy(x->setele, pe->setele);
  x->superset_pos[0] = 0;
  for (t = 1; t < MAXSUPSET; t++) {
    dim_t sup = sets[s].subsetid[t];
    if (sup == -1) break;
    if (sup == par) { x->superset_pos[t] = e; continue; }
    for (u = 1; u < MAXSUPSET; u++) if (sets[par].subsetid[u] == sup) break;
    x->superset_pos[t] = (u < MAXSUPSET) ? pe->superset_pos[u] : -1;
  }
}

/* parent set of an index's set when that set is a loop set, else -1 */
dim_t loop_set_parent(dim_t s) {
  int k;
  for (k = 0; k < teems_loop_nsyn; k++) if (teems_loop_syns[k].setid == s) return teems_loop_syns[k].parentid;
  return -1;
}
