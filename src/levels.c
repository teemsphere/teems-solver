#include <teems_solver.h>

/* C0: minimal levels-equation linearization (GEMPACK manual 9.2, 18.1,
   18.2.1, 10.9.1; design: docs/mapping_complementarity_design.md §5).

   tab_levels_transform runs on the preprocessed one-statement-per-line
   stream (after tab_postsim_split, before the readers) and rewrites:

   - each levels Variable declaration X into the GEMPACK associated
     triple (9.2.2): the declaration itself passes through VERBATIM
     (variables_read re-derives level_par/change_real from the same
     qualifiers and positional defaults) and the pass appends
       coefficient (non_parameter) <quants> X<args> ;
       update [(change)] <quants> X<args> = p_X<args> ;
     The pair shares one name ("auto-pair"): bare tokens bind the
     coefficient (value reference), p_/c_ tokens bind the variable
     (names_validate exemption + binder ordering in formula.c).

   - Formula&Equation into Formula (initial) + Equation (levels)
     (10.9.1), the latter then linearized like any levels equation.

   - Equation (levels) NAME <quants> LHS = RHS ; into an ordinary
     linear equation by CHANGE differentiation (18.1; TABLO's ACD
     mode -- every corpus levels equation meets a 9.2.6 change-diff
     trigger, and ACD is a sanctioned global TABLO option):
       d(A+B)=dA+dB   d(A-B)=dA-dB   d(A*B)=B*dA+A*dB
       d(A/B)=(1/B)*dA-(A/B^2)*dB
       d(A^B)=B*A^(B-1)*dA+loge(A)*A^B*dB
       d(sum(j,S,A))=sum(j,S,dA)     d(parameter)=d(constant)=0
     terminals: dX -> p_X (change) or (X/100)*p_X (percent).

   Levels equations may reference ONLY levels variables, parameters
   and constants (11.4.8) -- non-parameter coefficients, linear
   variables and unknown names are named fatals here. Functions are
   deferred (SQRT/EXP/LOGE/LOG10 are the only legal ones, zero corpus
   uses).

   Levels variable names starting p_/c_ (GMig2 P_L; design doc
   section 7, C1a): the pair coefficient cannot share the declared
   name -- a bare value token would be indistinguishable from the
   p_-prefixed column text the equation scanners key on -- so it gets
   a generated name gen_lvN (no p_/c_ substring) and every bare
   value-reference token of the declared name is rewritten to gen_lvN
   across the statement stream (whole-token match, double-quoted
   spans excluded) and the _ps companion. The Variable declaration
   keeps the declared name (closures, shocks, solution, compose all
   see the user's name); explicit p_<name>/c_<name> column tokens are
   distinct tokens, untouched, and resolve via linvar_resolve
   (section 6). names_validate's auto-pair exemption excludes
   p_/c_-leading names, so a user coefficient sharing such a name is
   a genuine clash fatal there.

   Zero-cost contract: a TAB with no levels statements is returned
   untouched (the goldens' bit-identity gate rides on this). */

#define LV_MAXLV 512
#define LV_MAXNAMES 4096
#define LV_MAXNODES 2048
#define LV_MAXTERMS 512
#define LV_MAXSUM 8
#define LV_TERMBUF 4096

typedef struct {
  char name[NAMESIZE];      /* declared name (variable keeps it) */
  char valname[NAMESIZE];   /* pair-coefficient name: == name unless
                               p_/c_-leading, then gen_lvN (C1a) */
  bool change;
} lv_pair;

typedef struct {
  char name[NAMESIZE];
  bool param;
} lv_coefrec;

typedef enum { LN_NUM, LN_REF, LN_ADD, LN_SUB, LN_MUL, LN_DIV, LN_POW, LN_NEG, LN_SUM } lv_kind;
/* REF classification */
typedef enum { LR_CHANGE, LR_PERCENT, LR_PARAM } lv_refkind;

typedef struct {
  lv_kind kind;
  int a, b;                 /* child node ids (-1 none) */
  const char *s; int slen;  /* source span of the whole subexpression */
  const char *name; int namelen;   /* REF name / SUM index */
  const char *args; int argslen;   /* REF index group incl. parens ("" if scalar) */
  const char *set; int setlen;     /* SUM set name */
  lv_refkind refkind;
  int pairidx;                     /* lv[] index for levels-var REFs */
} lv_node;

typedef struct {
  int sign;                     /* +1 / -1 */
  char coef[LV_TERMBUF];        /* '*'-joined factor text, "" if none */
  int var;                      /* node id of the dX terminal */
  int nsum;
  const char *sumidx[LV_MAXSUM]; int sumidxlen[LV_MAXSUM];
  const char *sumset[LV_MAXSUM]; int sumsetlen[LV_MAXSUM];
} lv_term;

typedef struct {
  lv_pair lv[LV_MAXLV]; int nlv;
  int nrenamed;                 /* gen_lvN counter (C1a) */
  lv_coefrec cf[LV_MAXNAMES]; int ncf;
  char linvar[LV_MAXNAMES][NAMESIZE]; int nlinvar; /* non-levels variables */
  /* per-statement expression workspace */
  lv_node nodes[LV_MAXNODES]; int nnodes;
  lv_term terms[LV_MAXTERMS]; int nterms;
  const char *stmt;             /* current statement text for messages */
} lv_ctx;

static int lv_strcat_b(char *dst, const char *src, size_t cap) {
  size_t dl = strlen(dst), sl = strlen(src);
  if (dl + sl + 1 > cap) return -1;
  memcpy(dst + dl, src, sl + 1);
  return 0;
}
static int lv_ncat_b(char *dst, const char *src, int n, size_t cap) {
  size_t dl = strlen(dst);
  if (dl + (size_t)n + 1 > cap) return -1;
  memcpy(dst + dl, src, n);
  dst[dl + n] = '\0';
  return 0;
}

static int lv_err(lv_ctx *c, const char *msg) {
  printf("Error: %s in levels statement: %s\n", msg, c->stmt ? c->stmt : "");
  return -1;
}

static bool lv_isnamec(char ch) {
  return (ch >= 'a' && ch <= 'z') || (ch >= '0' && ch <= '9') || ch == '_' || ch == '@' || (ch >= 'A' && ch <= 'Z');
}

/* expression tokens arrive post-rename: match the VALUE name */
static int lv_find_lv(lv_ctx *c, const char *name, int len) {
  int i;
  for (i = 0; i < c->nlv; i++)
    if ((int)strlen(c->lv[i].valname) == len && strncmp(c->lv[i].valname, name, len) == 0) return i;
  return -1;
}
/* declaration lines keep the DECLARED name */
static int lv_find_lv_decl(lv_ctx *c, const char *name, int len) {
  int i;
  for (i = 0; i < c->nlv; i++)
    if ((int)strlen(c->lv[i].name) == len && strncmp(c->lv[i].name, name, len) == 0) return i;
  return -1;
}
static int lv_find_cf(lv_ctx *c, const char *name, int len) {
  int i;
  for (i = 0; i < c->ncf; i++)
    if ((int)strlen(c->cf[i].name) == len && strncmp(c->cf[i].name, name, len) == 0) return i;
  return -1;
}
static int lv_find_linvar(lv_ctx *c, const char *name, int len) {
  int i;
  for (i = 0; i < c->nlinvar; i++)
    if ((int)strlen(c->linvar[i]) == len && strncmp(c->linvar[i], name, len) == 0) return i;
  return -1;
}

/* -------- statement head parsing (shared by scan and rewrite) -------- */

/* after the keyword: skip leading qualifier groups (a '(' group whose
   first token is not "all"), collecting the tokens; then skip
   quantifier groups "(all,...)" (capturing their verbatim span); then
   the declared name and optional index group. Returns 0 on success. */
typedef struct {
  bool has_levels, has_linear, has_change, has_percent;
  const char *quants; int quantslen;   /* verbatim quantifier span */
  const char *name; int namelen;
  const char *args; int argslen;       /* incl. parens, "" if none */
} lv_head;

static int lv_group_end(const char *p) {
  int depth = 0, i = 0;
  if (p[0] != '(') return -1;
  for (i = 0; p[i] != '\0'; i++) {
    if (p[i] == '(') depth++;
    else if (p[i] == ')') { depth--; if (depth == 0) return i; }
  }
  return -1;
}

static int lv_parse_head(const char *after_kw, lv_head *h) {
  const char *p = after_kw;
  int ge;
  memset(h, 0, sizeof(*h));
  h->args = ""; h->argslen = 0;
  while (*p == ' ') p++;
  /* qualifier groups */
  while (*p == '(') {
    if (strncmp(p, "(all,", 5) == 0 || strncmp(p, "(all ", 5) == 0) break;
    ge = lv_group_end(p);
    if (ge < 0) return -1;
    {
      /* tokenize group body on ',' with blanks stripped */
      char tok[NAMESIZE]; int ti = 0, i;
      for (i = 1; i <= ge; i++) {
        char ch = p[i];
        if (ch == ',' || ch == ')') {
          tok[ti] = '\0';
          if (strcmp(tok, "levels") == 0) h->has_levels = true;
          else if (strcmp(tok, "linear") == 0) h->has_linear = true;
          else if (strcmp(tok, "change") == 0) h->has_change = true;
          else if (strcmp(tok, "percent_change") == 0) h->has_percent = true;
          /* everything else (orig_level=, vpqtype=, bounds...) is
             validated by tab_qualifiers_parse downstream */
          ti = 0;
        } else if (ch != ' ' && ti < NAMESIZE - 1) tok[ti++] = ch;
      }
    }
    p += ge + 1;
    while (*p == ' ') p++;
  }
  /* quantifier groups */
  h->quants = p;
  while (*p == '(') {
    if (strncmp(p, "(all,", 5) != 0 && strncmp(p, "(all ", 5) != 0) break;
    ge = lv_group_end(p);
    if (ge < 0) return -1;
    p += ge + 1;
    while (*p == ' ') p++;
  }
  h->quantslen = (int)(p - h->quants);
  /* trim trailing blanks out of the quant span */
  while (h->quantslen > 0 && h->quants[h->quantslen - 1] == ' ') h->quantslen--;
  /* name */
  h->name = p;
  while (lv_isnamec(*p)) p++;
  h->namelen = (int)(p - h->name);
  if (h->namelen == 0) return -1;
  while (*p == ' ') p++;
  if (*p == '(') {
    ge = lv_group_end(p);
    if (ge < 0) return -1;
    h->args = p; h->argslen = ge + 1;
  }
  return 0;
}

/* separate parameter-ish token scan for coefficient statements (the
   head parser above only tracks variable tokens) */
static bool lv_coef_is_param_stmt(const char *after_kw, bool sticky_param) {
  const char *p = after_kw;
  bool par = sticky_param;
  int ge, i, ti;
  char tok[NAMESIZE];
  while (*p == ' ') p++;
  while (*p == '(') {
    if (strncmp(p, "(all,", 5) == 0 || strncmp(p, "(all ", 5) == 0) break;
    ge = lv_group_end(p);
    if (ge < 0) break;
    ti = 0;
    for (i = 1; i <= ge; i++) {
      char ch = p[i];
      if (ch == ',' || ch == ')') {
        tok[ti] = '\0';
        if (strcmp(tok, "parameter") == 0) par = true;
        else if (strcmp(tok, "non_parameter") == 0) par = false;
        else if (strcmp(tok, "integer") == 0) par = true; /* 10.3: INTEGER default PARAMETER */
        ti = 0;
      } else if (ch != ' ' && ti < NAMESIZE - 1) tok[ti++] = ch;
    }
    p += ge + 1;
    while (*p == ' ') p++;
  }
  return par;
}

/* -------- expression parser (recursive descent, spaces stripped) ---- */

typedef struct { const char *p; lv_ctx *c; int fail; } lv_cur;

static int lv_parse_expr(lv_cur *cur);

static int lv_newnode(lv_ctx *c, lv_kind k) {
  if (c->nnodes >= LV_MAXNODES) return -1;
  memset(&c->nodes[c->nnodes], 0, sizeof(lv_node));
  c->nodes[c->nnodes].kind = k;
  c->nodes[c->nnodes].a = -1;
  c->nodes[c->nnodes].b = -1;
  c->nodes[c->nnodes].args = ""; c->nodes[c->nnodes].argslen = 0;
  return c->nnodes++;
}

static const char *lv_funcs[] = { "sqrt", "exp", "loge", "log10", "abs", "max", "min", "if",
                                  "id01", "id0v", "round", "trunc0", "truncb", "prod", NULL };

static int lv_parse_primary(lv_cur *cur) {
  lv_ctx *c = cur->c;
  const char *start = cur->p;
  int n, ge, f;
  if (*cur->p == '(') {
    ge = lv_group_end(cur->p);
    if (ge < 0) { cur->fail = 1; lv_err(c, "unbalanced parentheses"); return -1; }
    {
      lv_cur inner = { 0 };
      char save = 0;
      /* parse the inside by bounding with a temporary terminator */
      inner.c = c;
      inner.p = cur->p + 1;
      save = *(cur->p + ge);
      *(char *)(cur->p + ge) = '\0';
      n = lv_parse_expr(&inner);
      *(char *)(cur->p + ge) = save;
      if (n < 0 || inner.fail) { cur->fail = 1; return -1; }
      /* the terminator was restored above: end-of-span is a pointer
         test, not a NUL test */
      if (inner.p != cur->p + ge) { cur->fail = 1; lv_err(c, "trailing text inside parentheses"); return -1; }
      cur->p += ge + 1;
      /* keep the parenthesized span on the node for verbatim reuse */
      c->nodes[n].s = start; c->nodes[n].slen = (int)(cur->p - start);
      return n;
    }
  }
  if ((*cur->p >= '0' && *cur->p <= '9') || *cur->p == '.') {
    n = lv_newnode(c, LN_NUM);
    if (n < 0) { cur->fail = 1; lv_err(c, "expression too large"); return -1; }
    while ((*cur->p >= '0' && *cur->p <= '9') || *cur->p == '.') cur->p++;
    c->nodes[n].s = start; c->nodes[n].slen = (int)(cur->p - start);
    return n;
  }
  if (lv_isnamec(*cur->p)) {
    const char *nm = cur->p;
    int nmlen;
    while (lv_isnamec(*cur->p)) cur->p++;
    nmlen = (int)(cur->p - nm);
    if (*cur->p == '(') {
      ge = lv_group_end(cur->p);
      if (ge < 0) { cur->fail = 1; lv_err(c, "unbalanced parentheses"); return -1; }
      if (nmlen == 3 && strncmp(nm, "sum", 3) == 0) {
        /* sum(idx,SET,body) */
        const char *g = cur->p + 1, *gend = cur->p + ge;
        const char *c1 = NULL, *c2 = NULL, *q;
        int depth = 0;
        for (q = g; q < gend; q++) {
          if (*q == '(') depth++;
          else if (*q == ')') depth--;
          else if (*q == ',' && depth == 0) { if (c1 == NULL) c1 = q; else { c2 = q; break; } }
          else if (*q == ':' && depth == 0) {
            /* the loop stops at the second top-level comma, so this
               covers the index and set parts (sum(r,REG: cond, ...)) */
            cur->fail = 1; lv_err(c, "conditional sums are not supported in levels equations yet"); return -1;
          }
        }
        if (c1 == NULL || c2 == NULL) { cur->fail = 1; lv_err(c, "malformed sum"); return -1; }
        n = lv_newnode(c, LN_SUM);
        if (n < 0) { cur->fail = 1; lv_err(c, "expression too large"); return -1; }
        c->nodes[n].name = g; c->nodes[n].namelen = (int)(c1 - g);
        c->nodes[n].set = c1 + 1; c->nodes[n].setlen = (int)(c2 - c1 - 1);
        {
          lv_cur inner = { 0 };
          char save;
          int b;
          inner.c = c; inner.p = c2 + 1;
          save = *gend; *(char *)gend = '\0';
          b = lv_parse_expr(&inner);
          *(char *)gend = save;
          if (b < 0 || inner.fail) { cur->fail = 1; return -1; }
          /* terminator restored above: end-of-span is a pointer test */
          if (inner.p != gend) { cur->fail = 1; lv_err(c, "trailing text inside sum"); return -1; }
          c->nodes[n].b = b;
        }
        cur->p += ge + 1;
        c->nodes[n].s = start; c->nodes[n].slen = (int)(cur->p - start);
        return n;
      }
      for (f = 0; lv_funcs[f] != NULL; f++) {
        if ((int)strlen(lv_funcs[f]) == nmlen && strncmp(nm, lv_funcs[f], nmlen) == 0) {
          cur->fail = 1;
          lv_err(c, "functions are not supported in levels equations yet (manual 11.4.10 allows only SQRT/EXP/LOGE/LOG10; linearize by hand or wait for the follow-on)");
          return -1;
        }
      }
      /* array reference */
      n = lv_newnode(c, LN_REF);
      if (n < 0) { cur->fail = 1; lv_err(c, "expression too large"); return -1; }
      c->nodes[n].name = nm; c->nodes[n].namelen = nmlen;
      c->nodes[n].args = cur->p; c->nodes[n].argslen = ge + 1;
      cur->p += ge + 1;
    } else {
      n = lv_newnode(c, LN_REF);
      if (n < 0) { cur->fail = 1; lv_err(c, "expression too large"); return -1; }
      c->nodes[n].name = nm; c->nodes[n].namelen = nmlen;
    }
    /* classify the reference (11.4.8) */
    {
      int k = lv_find_lv(c, c->nodes[n].name, c->nodes[n].namelen);
      if (k >= 0) { c->nodes[n].refkind = c->lv[k].change ? LR_CHANGE : LR_PERCENT; c->nodes[n].pairidx = k; }
      else {
        k = lv_find_cf(c, c->nodes[n].name, c->nodes[n].namelen);
        if (k >= 0) {
          if (!c->cf[k].param) {
            cur->fail = 1;
            lv_err(c, "non-parameter coefficient in a levels equation (manual 11.4.8: only levels variables, parameters and constants)");
            return -1;
          }
          c->nodes[n].refkind = LR_PARAM;
        } else if (lv_find_linvar(c, c->nodes[n].name, c->nodes[n].namelen) >= 0) {
          cur->fail = 1;
          lv_err(c, "linear variable in a levels equation (manual 11.4.8)");
          return -1;
        } else {
          cur->fail = 1;
          lv_err(c, "unknown name in a levels equation (not a levels variable, parameter or constant)");
          return -1;
        }
      }
    }
    c->nodes[n].s = start; c->nodes[n].slen = (int)(cur->p - start);
    return n;
  }
  cur->fail = 1;
  lv_err(c, "unexpected character in expression");
  return -1;
}

static int lv_parse_unary(lv_cur *cur) {
  const char *start = cur->p;
  int n, a;
  if (*cur->p == '-') {
    cur->p++;
    a = lv_parse_unary(cur);
    if (a < 0) return -1;
    n = lv_newnode(cur->c, LN_NEG);
    if (n < 0) { cur->fail = 1; lv_err(cur->c, "expression too large"); return -1; }
    cur->c->nodes[n].a = a;
    cur->c->nodes[n].s = start; cur->c->nodes[n].slen = (int)(cur->p - start);
    return n;
  }
  if (*cur->p == '+') { cur->p++; return lv_parse_unary(cur); }
  return lv_parse_primary(cur);
}

static int lv_parse_pow(lv_cur *cur) {
  const char *start = cur->p;
  int a = lv_parse_unary(cur), b, n;
  if (a < 0) return -1;
  while (*cur->p == '^') {
    cur->p++;
    b = lv_parse_unary(cur);
    if (b < 0) return -1;
    n = lv_newnode(cur->c, LN_POW);
    if (n < 0) { cur->fail = 1; lv_err(cur->c, "expression too large"); return -1; }
    cur->c->nodes[n].a = a; cur->c->nodes[n].b = b;
    cur->c->nodes[n].s = start; cur->c->nodes[n].slen = (int)(cur->p - start);
    a = n;
  }
  return a;
}

static int lv_parse_muldiv(lv_cur *cur) {
  const char *start = cur->p;
  int a = lv_parse_pow(cur), b, n;
  if (a < 0) return -1;
  while (*cur->p == '*' || *cur->p == '/') {
    char op = *cur->p;
    cur->p++;
    b = lv_parse_pow(cur);
    if (b < 0) return -1;
    n = lv_newnode(cur->c, op == '*' ? LN_MUL : LN_DIV);
    if (n < 0) { cur->fail = 1; lv_err(cur->c, "expression too large"); return -1; }
    cur->c->nodes[n].a = a; cur->c->nodes[n].b = b;
    cur->c->nodes[n].s = start; cur->c->nodes[n].slen = (int)(cur->p - start);
    a = n;
  }
  return a;
}

static int lv_parse_expr(lv_cur *cur) {
  const char *start = cur->p;
  int a = lv_parse_muldiv(cur), b, n;
  if (a < 0) return -1;
  while (*cur->p == '+' || *cur->p == '-') {
    char op = *cur->p;
    cur->p++;
    b = lv_parse_muldiv(cur);
    if (b < 0) return -1;
    n = lv_newnode(cur->c, op == '+' ? LN_ADD : LN_SUB);
    if (n < 0) { cur->fail = 1; lv_err(cur->c, "expression too large"); return -1; }
    cur->c->nodes[n].a = a; cur->c->nodes[n].b = b;
    cur->c->nodes[n].s = start; cur->c->nodes[n].slen = (int)(cur->p - start);
    a = n;
  }
  return a;
}

/* -------- differentiation (change mode, 18.1) ----------------------- */

/* append "(span)" of node id to a factor buffer */
static int lv_factor_span(lv_ctx *c, char *buf, int id) {
  if (lv_strcat_b(buf, "(", LV_TERMBUF) < 0) return -1;
  if (lv_ncat_b(buf, c->nodes[id].s, c->nodes[id].slen, LV_TERMBUF) < 0) return -1;
  if (lv_strcat_b(buf, ")", LV_TERMBUF) < 0) return -1;
  return 0;
}

typedef struct {
  int nsum;
  const char *sumidx[LV_MAXSUM]; int sumidxlen[LV_MAXSUM];
  const char *sumset[LV_MAXSUM]; int sumsetlen[LV_MAXSUM];
} lv_sumstack;

static int lv_diff(lv_ctx *c, int id, int sign, const char *factors, const lv_sumstack *ss) {
  lv_node *nd = &c->nodes[id];
  char fbuf[LV_TERMBUF];
  switch (nd->kind) {
  case LN_NUM:
    return 0;
  case LN_REF:
    if (nd->refkind == LR_PARAM) return 0;
    if (c->nterms >= LV_MAXTERMS) return lv_err(c, "too many terms after linearization");
    {
      lv_term *t = &c->terms[c->nterms++];
      int i;
      t->sign = sign;
      strcpy(t->coef, factors);
      t->var = id;
      t->nsum = ss->nsum;
      for (i = 0; i < ss->nsum; i++) {
        t->sumidx[i] = ss->sumidx[i]; t->sumidxlen[i] = ss->sumidxlen[i];
        t->sumset[i] = ss->sumset[i]; t->sumsetlen[i] = ss->sumsetlen[i];
      }
    }
    return 0;
  case LN_NEG:
    return lv_diff(c, nd->a, -sign, factors, ss);
  case LN_ADD:
    if (lv_diff(c, nd->a, sign, factors, ss) < 0) return -1;
    return lv_diff(c, nd->b, sign, factors, ss);
  case LN_SUB:
    if (lv_diff(c, nd->a, sign, factors, ss) < 0) return -1;
    return lv_diff(c, nd->b, -sign, factors, ss);
  case LN_MUL:
    /* d(A*B) = B*dA + A*dB */
    strcpy(fbuf, factors);
    if (fbuf[0] != '\0' && lv_strcat_b(fbuf, "*", LV_TERMBUF) < 0) return lv_err(c, "term too large");
    if (lv_factor_span(c, fbuf, nd->b) < 0) return lv_err(c, "term too large");
    if (lv_diff(c, nd->a, sign, fbuf, ss) < 0) return -1;
    strcpy(fbuf, factors);
    if (fbuf[0] != '\0' && lv_strcat_b(fbuf, "*", LV_TERMBUF) < 0) return lv_err(c, "term too large");
    if (lv_factor_span(c, fbuf, nd->a) < 0) return lv_err(c, "term too large");
    return lv_diff(c, nd->b, sign, fbuf, ss);
  case LN_DIV:
    /* d(A/B) = (1/B)*dA - (A/B^2)*dB */
    strcpy(fbuf, factors);
    if (fbuf[0] != '\0' && lv_strcat_b(fbuf, "*", LV_TERMBUF) < 0) return lv_err(c, "term too large");
    if (lv_strcat_b(fbuf, "(1/", LV_TERMBUF) < 0) return lv_err(c, "term too large");
    if (lv_factor_span(c, fbuf, nd->b) < 0) return lv_err(c, "term too large");
    if (lv_strcat_b(fbuf, ")", LV_TERMBUF) < 0) return lv_err(c, "term too large");
    if (lv_diff(c, nd->a, sign, fbuf, ss) < 0) return -1;
    strcpy(fbuf, factors);
    if (fbuf[0] != '\0' && lv_strcat_b(fbuf, "*", LV_TERMBUF) < 0) return lv_err(c, "term too large");
    if (lv_strcat_b(fbuf, "(", LV_TERMBUF) < 0) return lv_err(c, "term too large");
    if (lv_factor_span(c, fbuf, nd->a) < 0) return lv_err(c, "term too large");
    if (lv_strcat_b(fbuf, "/", LV_TERMBUF) < 0) return lv_err(c, "term too large");
    if (lv_factor_span(c, fbuf, nd->b) < 0) return lv_err(c, "term too large");
    if (lv_strcat_b(fbuf, "^2)", LV_TERMBUF) < 0) return lv_err(c, "term too large");
    return lv_diff(c, nd->b, -sign, fbuf, ss);
  case LN_POW:
    /* d(A^B) = B*A^(B-1)*dA + loge(A)*A^B*dB */
    strcpy(fbuf, factors);
    if (fbuf[0] != '\0' && lv_strcat_b(fbuf, "*", LV_TERMBUF) < 0) return lv_err(c, "term too large");
    if (lv_strcat_b(fbuf, "(", LV_TERMBUF) < 0) return lv_err(c, "term too large");
    if (lv_factor_span(c, fbuf, nd->b) < 0) return lv_err(c, "term too large");
    if (lv_strcat_b(fbuf, "*", LV_TERMBUF) < 0) return lv_err(c, "term too large");
    if (lv_factor_span(c, fbuf, nd->a) < 0) return lv_err(c, "term too large");
    if (lv_strcat_b(fbuf, "^(", LV_TERMBUF) < 0) return lv_err(c, "term too large");
    if (lv_ncat_b(fbuf, c->nodes[nd->b].s, c->nodes[nd->b].slen, LV_TERMBUF) < 0) return lv_err(c, "term too large");
    if (lv_strcat_b(fbuf, "-1))", LV_TERMBUF) < 0) return lv_err(c, "term too large");
    if (lv_diff(c, nd->a, sign, fbuf, ss) < 0) return -1;
    strcpy(fbuf, factors);
    if (fbuf[0] != '\0' && lv_strcat_b(fbuf, "*", LV_TERMBUF) < 0) return lv_err(c, "term too large");
    if (lv_strcat_b(fbuf, "(loge", LV_TERMBUF) < 0) return lv_err(c, "term too large");
    if (lv_factor_span(c, fbuf, nd->a) < 0) return lv_err(c, "term too large");
    if (lv_strcat_b(fbuf, "*", LV_TERMBUF) < 0) return lv_err(c, "term too large");
    if (lv_factor_span(c, fbuf, nd->a) < 0) return lv_err(c, "term too large");
    if (lv_strcat_b(fbuf, "^", LV_TERMBUF) < 0) return lv_err(c, "term too large");
    if (lv_factor_span(c, fbuf, nd->b) < 0) return lv_err(c, "term too large");
    if (lv_strcat_b(fbuf, ")", LV_TERMBUF) < 0) return lv_err(c, "term too large");
    return lv_diff(c, nd->b, sign, fbuf, ss);
  case LN_SUM: {
    lv_sumstack ss2 = *ss;
    if (ss2.nsum >= LV_MAXSUM) return lv_err(c, "sums nested too deeply");
    ss2.sumidx[ss2.nsum] = nd->name; ss2.sumidxlen[ss2.nsum] = nd->namelen;
    ss2.sumset[ss2.nsum] = nd->set; ss2.sumsetlen[ss2.nsum] = nd->setlen;
    ss2.nsum++;
    return lv_diff(c, nd->b, sign, factors, &ss2);
  }
  }
  return lv_err(c, "internal: unhandled node kind");
}

/* differentiate one side into out (starts with "0", then +/- terms) */
static int lv_diff_side(lv_ctx *c, const char *expr, char *out, size_t outcap) {
  lv_cur cur = { 0 };
  lv_sumstack ss = { 0 };
  char stripped[DATREADLINE];
  int i, j, root, t0 = c->nterms;
  /* strip blanks: spans stay valid inside this stable local copy */
  j = 0;
  for (i = 0; expr[i] != '\0'; i++) {
    if (expr[i] == ' ') continue;
    if (j >= DATREADLINE - 1) return lv_err(c, "levels expression too long");
    stripped[j++] = expr[i];
  }
  stripped[j] = '\0';
  if (j == 0 || (j == 1 && stripped[0] == '0')) {
    if (lv_strcat_b(out, "0", outcap) < 0) return lv_err(c, "output too large");
    return 0;
  }
  cur.c = c; cur.p = stripped;
  root = lv_parse_expr(&cur);
  if (root < 0 || cur.fail) return -1;
  if (*cur.p != '\0') return lv_err(c, "trailing text after expression");
  if (lv_diff(c, root, +1, "", &ss) < 0) return -1;
  if (lv_strcat_b(out, "0", outcap) < 0) return lv_err(c, "output too large");
  for (i = t0; i < c->nterms; i++) {
    lv_term *t = &c->terms[i];
    lv_node *v = &c->nodes[t->var];
    if (lv_strcat_b(out, t->sign > 0 ? "+" : "-", outcap) < 0) return lv_err(c, "output too large");
    for (j = 0; j < t->nsum; j++) {
      if (lv_strcat_b(out, "sum(", outcap) < 0) return lv_err(c, "output too large");
      if (lv_ncat_b(out, t->sumidx[j], t->sumidxlen[j], outcap) < 0) return lv_err(c, "output too large");
      if (lv_strcat_b(out, ",", outcap) < 0) return lv_err(c, "output too large");
      if (lv_ncat_b(out, t->sumset[j], t->sumsetlen[j], outcap) < 0) return lv_err(c, "output too large");
      if (lv_strcat_b(out, ",", outcap) < 0) return lv_err(c, "output too large");
    }
    if (t->coef[0] != '\0') {
      if (lv_strcat_b(out, t->coef, outcap) < 0) return lv_err(c, "output too large");
      if (lv_strcat_b(out, "*", outcap) < 0) return lv_err(c, "output too large");
    }
    if (v->refkind == LR_PERCENT) {
      /* (X/100)*p_X */
      if (lv_strcat_b(out, "(", outcap) < 0) return lv_err(c, "output too large");
      if (lv_ncat_b(out, v->name, v->namelen, outcap) < 0) return lv_err(c, "output too large");
      if (lv_ncat_b(out, v->args, v->argslen, outcap) < 0) return lv_err(c, "output too large");
      if (lv_strcat_b(out, "/100)*", outcap) < 0) return lv_err(c, "output too large");
    }
    /* the column text is p_ + DECLARED name (differs from the value
       token for gen_lv-renamed pairs); linvar_resolve resolves
       p_<p_-leading name> tokens (section 6) */
    if (lv_strcat_b(out, "p_", outcap) < 0) return lv_err(c, "output too large");
    if (lv_strcat_b(out, c->lv[v->pairidx].name, outcap) < 0) return lv_err(c, "output too large");
    if (lv_ncat_b(out, v->args, v->argslen, outcap) < 0) return lv_err(c, "output too large");
    for (j = 0; j < t->nsum; j++)
      if (lv_strcat_b(out, ")", outcap) < 0) return lv_err(c, "output too large");
  }
  return 0;
}

/* -------- scan pass -------------------------------------------------- */

static int lv_scan(lv_ctx *c, char *fname, bool *any) {
  FILE *f;
  char line[TABREADLINE], defval[NAMESIZE];
  bool sticky_level = false, sticky_change = false, sticky_param = false;
  lv_head h;
  *any = false;
  f = fopen(fname, "r");
  if (f == NULL) { printf("Error: cannot open %s\n", fname); return -1; }
  while (fgets(line, TABREADLINE, f)) {
    if (strncmp(line, "variable", 8) == 0) {
      if (tab_default_value(line, defval)) {
        if (strcmp(defval, "levels") == 0) sticky_level = true;
        else if (strcmp(defval, "linear") == 0) sticky_level = false;
        else if (strcmp(defval, "change") == 0) sticky_change = true;
        else if (strcmp(defval, "percent_change") == 0) sticky_change = false;
        continue;
      }
      if (lv_parse_head(line + 8, &h) < 0) continue; /* readers report malformed decls */
      {
        bool lv = sticky_level, ch = sticky_change;
        if (h.has_levels) lv = true;
        if (h.has_linear) lv = false;
        if (h.has_change) ch = true;
        if (h.has_percent) ch = false;
        if (lv) {
          if (c->nlv >= LV_MAXLV) { printf("Error: too many levels variables\n"); fclose(f); return -1; }
          {
            int l = h.namelen < NAMESIZE - 1 ? h.namelen : NAMESIZE - 1;
            strncpy(c->lv[c->nlv].name, h.name, l);
            c->lv[c->nlv].name[l] = '\0';
            c->lv[c->nlv].change = ch;
            /* p_-leading declared name: the pair coefficient gets a
               generated name and value references are rewritten to it
               (C1a; design doc section 7). c_-leading LEVELS names
               stay fatal: the preprocess c_->p_ rewrite folds their
               value references on equation/update lines into p_<tail>
               BEFORE this pass, indistinguishable from a change-
               reference column of the same variable (zero corpus
               uses; GMig2's c_* class is linear, section 6). */
            if (l >= 2 && c->lv[c->nlv].name[0] == 'c' && c->lv[c->nlv].name[1] == '_') {
              printf("Error: levels variable %s: c_-leading names are indistinguishable from change-reference columns after the c_->p_ rewrite; rename the variable\n", c->lv[c->nlv].name);
              fclose(f);
              return -1;
            }
            if (l >= 2 && c->lv[c->nlv].name[0] == 'p' && c->lv[c->nlv].name[1] == '_')
              sprintf(c->lv[c->nlv].valname, "gen_lv%d", c->nrenamed++);
            else strcpy(c->lv[c->nlv].valname, c->lv[c->nlv].name);
            c->nlv++;
          }
          *any = true;
        } else {
          if (c->nlinvar < LV_MAXNAMES) {
            int l = h.namelen < NAMESIZE - 1 ? h.namelen : NAMESIZE - 1;
            strncpy(c->linvar[c->nlinvar], h.name, l);
            c->linvar[c->nlinvar][l] = '\0';
            c->nlinvar++;
          }
        }
      }
    } else if (strncmp(line, "coefficient", 11) == 0) {
      if (tab_default_value(line, defval)) {
        if (strcmp(defval, "parameter") == 0) sticky_param = true;
        else if (strcmp(defval, "non_parameter") == 0) sticky_param = false;
        continue;
      }
      if (lv_parse_head(line + 11, &h) < 0) continue;
      if (c->ncf < LV_MAXNAMES) {
        int l = h.namelen < NAMESIZE - 1 ? h.namelen : NAMESIZE - 1;
        strncpy(c->cf[c->ncf].name, h.name, l);
        c->cf[c->ncf].name[l] = '\0';
        c->cf[c->ncf].param = lv_coef_is_param_stmt(line + 11, sticky_param);
        c->ncf++;
      }
    } else if (strncmp(line, "formula", 7) == 0) {
      if (strstr(line, "&equation") != NULL || strstr(line, "& equation") != NULL) *any = true;
    } else if (strncmp(line, "equation", 8) == 0) {
      const char *p = line + 8;
      while (*p == ' ') p++;
      if (strncmp(p, "(levels", 7) == 0) *any = true;
    }
  }
  fclose(f);
  return 0;
}

/* -------- rewrite pass ----------------------------------------------- */

/* C1a: whole-token rewrite of renamed pair value references (declared
   p_/c_-leading name -> gen_lvN). Double-quoted spans (set elements)
   are copied verbatim; p_<name>/c_<name> column references are longer
   tokens and never match. Returns 0 (also when nothing changed), -1
   on overflow. */
static int lv_rename_line(lv_ctx *c, char *line, size_t cap) {
  char buf[TABREADLINE];
  int i = 0, k, tl, changed = 0;
  size_t bi = 0, vl;
  if (c->nrenamed == 0) return 0;
  while (line[i] != '\0') {
    if (line[i] == '"') {
      if (bi + 1 >= sizeof(buf)) return -1;
      buf[bi++] = line[i++];
      while (line[i] != '\0' && line[i] != '"') {
        if (bi + 1 >= sizeof(buf)) return -1;
        buf[bi++] = line[i++];
      }
      if (line[i] == '"') {
        if (bi + 1 >= sizeof(buf)) return -1;
        buf[bi++] = line[i++];
      }
      continue;
    }
    if (lv_isnamec(line[i]) && (i == 0 || !lv_isnamec(line[i - 1]))) {
      tl = 0;
      while (lv_isnamec(line[i + tl])) tl++;
      for (k = 0; k < c->nlv; k++)
        if (strcmp(c->lv[k].name, c->lv[k].valname) != 0 &&
            (int)strlen(c->lv[k].name) == tl &&
            strncmp(c->lv[k].name, line + i, tl) == 0) break;
      if (k < c->nlv) {
        vl = strlen(c->lv[k].valname);
        if (bi + vl + 1 >= sizeof(buf)) return -1;
        memcpy(buf + bi, c->lv[k].valname, vl);
        bi += vl;
        i += tl;
        changed = 1;
        continue;
      }
      if (bi + (size_t)tl + 1 >= sizeof(buf)) return -1;
      memcpy(buf + bi, line + i, tl);
      bi += tl;
      i += tl;
      continue;
    }
    if (bi + 1 >= sizeof(buf)) return -1;
    buf[bi++] = line[i++];
  }
  buf[bi] = '\0';
  if (changed) {
    if (bi + 1 > cap) return -1;
    memcpy(line, buf, bi + 1);
  }
  return 0;
}

/* linearize one "equation (levels) ..." line (keyword already
   consumed by the caller: p points after "equation") into fout */
static int lv_emit_linearized(lv_ctx *c, const char *p, FILE *fout) {
  char out[DATREADLINE];
  const char *name, *label = NULL, *quants, *lhs;
  int namelen, labellen = 0, quantslen, ge;
  char lhsbuf[DATREADLINE], rhsbuf[DATREADLINE];
  const char *eq, *semi;
  c->nnodes = 0; c->nterms = 0;
  while (*p == ' ') p++;
  if (strncmp(p, "(levels)", 8) != 0) {
    if (strncmp(p, "(levels", 7) == 0)
      return lv_err(c, "unsupported qualifier combination on Equation (levels)");
    return lv_err(c, "internal: expected (levels)");
  }
  p += 8;
  while (*p == ' ') p++;
  name = p;
  while (lv_isnamec(*p)) p++;
  namelen = (int)(p - name);
  if (namelen == 0) return lv_err(c, "missing equation name");
  while (*p == ' ') p++;
  if (*p == '#') {
    label = p;
    p++;
    while (*p != '\0' && *p != '#') p++;
    if (*p != '#') return lv_err(c, "unterminated # label #");
    p++;
    labellen = (int)(p - label);
    while (*p == ' ') p++;
  }
  quants = p;
  while (*p == '(') {
    if (strncmp(p, "(all,", 5) != 0 && strncmp(p, "(all ", 5) != 0) break;
    ge = lv_group_end(p);
    if (ge < 0) return lv_err(c, "unbalanced quantifier");
    {
      const char *q;
      for (q = p; q <= p + ge; q++)
        if (*q == ':') return lv_err(c, "conditions on Equation quantifiers are not supported (manual 11.4.11)");
    }
    p += ge + 1;
    while (*p == ' ') p++;
  }
  quantslen = (int)(p - quants);
  while (quantslen > 0 && quants[quantslen - 1] == ' ') quantslen--;
  /* LHS = RHS ; */
  eq = strchr(p, '=');
  if (eq == NULL) return lv_err(c, "missing = in levels equation");
  semi = strchr(eq, ';');
  if (semi == NULL) return lv_err(c, "missing ; in levels equation");
  {
    int ll = (int)(eq - p), rl = (int)(semi - eq - 1);
    if (ll >= DATREADLINE || rl >= DATREADLINE) return lv_err(c, "levels equation too long");
    strncpy(lhsbuf, p, ll); lhsbuf[ll] = '\0';
    strncpy(rhsbuf, eq + 1, rl); rhsbuf[rl] = '\0';
  }
  lhs = lhsbuf;
  out[0] = '\0';
  if (lv_strcat_b(out, "equation ", sizeof(out)) < 0) return lv_err(c, "output too large");
  if (lv_ncat_b(out, name, namelen, sizeof(out)) < 0) return lv_err(c, "output too large");
  if (lv_strcat_b(out, " ", sizeof(out)) < 0) return lv_err(c, "output too large");
  if (labellen > 0) {
    if (lv_ncat_b(out, label, labellen, sizeof(out)) < 0) return lv_err(c, "output too large");
    if (lv_strcat_b(out, " ", sizeof(out)) < 0) return lv_err(c, "output too large");
  }
  if (quantslen > 0) {
    if (lv_ncat_b(out, quants, quantslen, sizeof(out)) < 0) return lv_err(c, "output too large");
    if (lv_strcat_b(out, " ", sizeof(out)) < 0) return lv_err(c, "output too large");
  }
  if (lv_diff_side(c, lhs, out, sizeof(out)) < 0) return -1;
  if (lv_strcat_b(out, " = ", sizeof(out)) < 0) return lv_err(c, "output too large");
  if (lv_diff_side(c, rhsbuf, out, sizeof(out)) < 0) return -1;
  if (lv_strcat_b(out, " ;", sizeof(out)) < 0) return lv_err(c, "output too large");
  /* sibling generated sums can share an index name inside one
     statement; the preprocess dedup pass has already run, so apply it
     to the generated line here */
  sum_dedup_indices(out);
  if (strlen(out) >= TABREADLINE - 2) {
    printf("Error: linearized form of levels equation %.*s exceeds the statement buffer; split the equation\n", namelen, name);
    return -1;
  }
  fprintf(fout, "%s\n", out);
  return 0;
}

/* ==================================================================
   C1: Complementarity parse/validation + derived statements (manual
   10.17/11.14/51.7.2; design doc section 7). Runs BEFORE
   tab_levels_transform: the derived statements are levels statements
   the C0 machinery then auto-pairs and linearizes.

     complementarity (variable = x[,lower_bound = b][,upper_bound = b])
         name (all,i,S)... expr ;

   Derived per statement (GEMPACK's $-prefixed names use '@' suffixes
   instead -- '$' does not survive the tokenizers, '@' does and is
   illegal in user names, so the namespace is collision-free):
     variable (change,levels) <quants> name@e(idxs) ;
     formula (initial) <quants> name@e(idxs) = expr ;
     equation (levels) e_name@e <quants> name@e(idxs) = expr ;
     variable (change) <quants> name@d(idxs) ;           [dummy]
     variable (change) del_comp@ ;                       [once per TAB]
   plus, when a bound is a levels variable (51.7.2 (iv)/(v)):
     variable (change,levels) <quants> name@l(idxs) ;
     formula (initial)/equation (levels): name@l = x(idxs) - l(idxs)
   The statement itself is consumed. Set matching (11.14 points 2-3)
   is validated in complementarities_validate once set elements
   exist; closure integration (auto-exogenize the dummy and del_comp@,
   the 11.14.1 backsolve guard and the C2 state-machinery guard) in
   comp_closure_check. The equation row E_$comp that lets the
   complementarity BIND arrives with C2; at C1 a complementarity is
   solvable only in the inert closure (x fully exogenous). */

/* per-arg-position quantifier sets of a declaration head: arg k's
   index name is looked up in the head's quantifier span */
static int cp_argsets(lv_ctx *c, const lv_head *h, char argsets[][NAMESIZE], dim_t *nargs) {
  const char *p;
  char idx[NAMESIZE], find[NAMESIZE + 8];
  int ti, k = 0;
  *nargs = 0;
  if (h->argslen == 0) return 0;
  p = h->args + 1;                      /* inside '(' */
  while (*p != '\0' && *p != ')' && p < h->args + h->argslen) {
    ti = 0;
    while (*p != ',' && *p != ')' && *p != '\0') {
      if (*p != ' ' && ti < NAMESIZE - 1) idx[ti++] = *p;
      p++;
    }
    idx[ti] = '\0';
    if (*p == ',') p++;
    if (ti == 0) continue;
    if (k >= MAXVARDIM) return lv_err(c, "too many arguments in declaration");
    /* find "(all,<idx>," in the quantifier span */
    snprintf(find, sizeof(find), "(all,%s,", idx);
    {
      const char *q = h->quants, *hit = NULL, *e;
      int fl = (int)strlen(find);
      while (q + fl <= h->quants + h->quantslen) {
        if (strncmp(q, find, fl) == 0) { hit = q + fl; break; }
        q++;
      }
      if (hit == NULL) return lv_err(c, "declaration argument has no matching (all,...) quantifier");
      e = hit;
      while (*e != ')' && *e != '\0') e++;
      if ((int)(e - hit) >= NAMESIZE) return lv_err(c, "set name too long");
      strncpy(argsets[k], hit, e - hit);
      argsets[k][e - hit] = '\0';
    }
    k++;
  }
  *nargs = k;
  return 0;
}

/* parse one complementarity statement (post-preprocess: lowercase,
   one line, comments stripped) into cp; expression text into expr */
static int cp_parse_stmt(lv_ctx *c, char *line, comp_def *cp, char *expr, size_t exprcap) {
  const char *p = line + 15;            /* past "complementarity" */
  int ge, i;
  c->stmt = line;
  memset(cp, 0, sizeof(*cp));
  while (*p == ' ') p++;
  if (*p != '(') return lv_err(c, "Complementarity needs a (variable = ..., lower_bound/upper_bound = ...) qualifier (manual 10.17)");
  ge = lv_group_end(p);
  if (ge < 0) return lv_err(c, "unbalanced Complementarity qualifier");
  {
    /* tokenize the group body on commas, strip blanks, split at '=' */
    char part[NAMESIZE * 2];
    int ti = 0;
    for (i = 1; i <= ge; i++) {
      char ch = p[i];
      if (ch == ',' || ch == ')') {
        part[ti] = '\0';
        if (ti > 0) {
          char *eq = strchr(part, '=');
          if (eq == NULL) return lv_err(c, "malformed Complementarity qualifier entry (expect variable/lower_bound/upper_bound = value)");
          *eq = '\0';
          if (strcmp(part, "variable") == 0) {
            if (strlen(eq + 1) == 0 || strlen(eq + 1) >= NAMESIZE) return lv_err(c, "malformed variable = entry");
            strcpy(cp->varname, eq + 1);
          } else if (strcmp(part, "lower_bound") == 0 || strcmp(part, "upper_bound") == 0) {
            int lower = (part[0] == 'l');
            char *val = eq + 1, *end;
            double d;
            if (strlen(val) == 0) return lv_err(c, "empty bound in Complementarity qualifier");
            d = strtod(val, &end);
            if (end != val && *end == '\0') {
              if (lower) { cp->lower_kind = 1; cp->lower_const = d; }
              else { cp->upper_kind = 1; cp->upper_const = d; }
            } else {
              if (strlen(val) >= NAMESIZE) return lv_err(c, "bound name too long");
              if (lv_find_lv(c, val, (int)strlen(val)) >= 0 || lv_find_lv_decl(c, val, (int)strlen(val)) >= 0) {
                if (lower) { cp->lower_kind = 2; strcpy(cp->lower_name, val); }
                else { cp->upper_kind = 2; strcpy(cp->upper_name, val); }
              } else if (lv_find_cf(c, val, (int)strlen(val)) >= 0) {
                if (!c->cf[lv_find_cf(c, val, (int)strlen(val))].param)
                  return lv_err(c, "a Complementarity bound must be a levels variable, a Coefficient(parameter) or a real constant (manual 10.17); non-parameter coefficient");
                if (lower) { cp->lower_kind = 3; strcpy(cp->lower_name, val); }
                else { cp->upper_kind = 3; strcpy(cp->upper_name, val); }
              } else if (lv_find_linvar(c, val, (int)strlen(val)) >= 0) {
                return lv_err(c, "a Complementarity bound must be a levels variable, a Coefficient(parameter) or a real constant (manual 10.17); linear variable");
              } else {
                return lv_err(c, "unknown name as Complementarity bound (must be a levels variable, a Coefficient(parameter) or a real constant)");
              }
            }
          } else {
            return lv_err(c, "unknown Complementarity qualifier keyword (expect variable, lower_bound, upper_bound)");
          }
        }
        ti = 0;
      } else if (ch != ' ' && ti < (int)sizeof(part) - 1) part[ti++] = ch;
    }
  }
  p += ge + 1;
  if (cp->varname[0] == '\0') return lv_err(c, "Complementarity needs a VARIABLE qualifier (manual 11.14)");
  if (cp->lower_kind == 0 && cp->upper_kind == 0) return lv_err(c, "Complementarity needs at least one of LOWER_BOUND/UPPER_BOUND (manual 10.17)");
  /* X must be a LEVELS variable */
  if (lv_find_lv(c, cp->varname, (int)strlen(cp->varname)) < 0 &&
      lv_find_lv_decl(c, cp->varname, (int)strlen(cp->varname)) < 0) {
    if (lv_find_linvar(c, cp->varname, (int)strlen(cp->varname)) >= 0)
      return lv_err(c, "the Complementarity variable must be a LEVELS variable (manual 11.14); linear variable");
    if (lv_find_cf(c, cp->varname, (int)strlen(cp->varname)) >= 0)
      return lv_err(c, "the Complementarity variable must be a LEVELS variable (manual 11.14); coefficient");
    return lv_err(c, "the Complementarity variable is not declared (must be a LEVELS variable, manual 11.14)");
  }
  /* name */
  while (*p == ' ') p++;
  {
    const char *nm = p;
    int nmlen;
    while (lv_isnamec(*p)) p++;
    nmlen = (int)(p - nm);
    if (nmlen == 0) return lv_err(c, "missing Complementarity name (manual 10.17)");
    if (nmlen > 10) return lv_err(c, "Complementarity name is limited to 10 characters (manual 11.14/11.2.1)");
    strncpy(cp->name, nm, nmlen);
    cp->name[nmlen] = '\0';
  }
  /* quantifiers */
  while (*p == ' ') p++;
  while (*p == '(') {
    if (strncmp(p, "(all,", 5) != 0) break;
    ge = lv_group_end(p);
    if (ge < 0) return lv_err(c, "unbalanced Complementarity quantifier");
    {
      char part[NAMESIZE];
      int ti = 0, fld = 0;
      for (i = 1; i <= ge; i++) {
        char ch = p[i];
        if (ch == ':') return lv_err(c, "conditions on Complementarity quantifiers are not supported");
        if (ch == ',' || ch == ')') {
          part[ti] = '\0';
          if (fld == 1) {
            if (cp->nquant >= MAXVARDIM) return lv_err(c, "too many Complementarity quantifiers");
            strcpy(cp->qidx[cp->nquant], part);
          } else if (fld == 2) strcpy(cp->qset[cp->nquant], part);
          fld++;
          ti = 0;
        } else if (ch != ' ' && ti < NAMESIZE - 1) part[ti++] = ch;
      }
      if (fld != 3) return lv_err(c, "malformed Complementarity quantifier (expect (all,index,set))");
      cp->nquant++;
    }
    p += ge + 1;
    while (*p == ' ') p++;
  }
  /* expression */
  {
    const char *semi = strrchr(p, ';');
    int el;
    if (semi == NULL) return lv_err(c, "missing ; in Complementarity statement");
    el = (int)(semi - p);
    while (el > 0 && p[el - 1] == ' ') el--;
    if (el <= 0) return lv_err(c, "missing Complementarity expression (manual 10.17)");
    if ((size_t)el >= exprcap) return lv_err(c, "Complementarity expression too long");
    strncpy(expr, p, el);
    expr[el] = '\0';
  }
  return 0;
}

/* second declaration scan: per-arg-position set names for X and the
   variable/coefficient bounds, plus the argument-count checks (11.14
   points 2-3) */
static int cp_fill_decl_sets(lv_ctx *c, char *fname) {
  FILE *f;
  char line[TABREADLINE], defval[NAMESIZE];
  lv_head h;
  dim_t k, nargs;
  char argsets[MAXVARDIM][NAMESIZE];
  f = fopen(fname, "r");
  if (f == NULL) { printf("Error: cannot open %s\n", fname); return -1; }
  while (fgets(line, TABREADLINE, f)) {
    int iscoef = (strncmp(line, "coefficient", 11) == 0);
    int isvar = (strncmp(line, "variable", 8) == 0);
    size_t sl = strlen(line);
    while (sl > 0 && (line[sl - 1] == '\n' || line[sl - 1] == '\r')) line[--sl] = '\0';
    if (!iscoef && !isvar) continue;
    if (tab_default_value(line, defval)) continue;
    if (lv_parse_head(line + (iscoef ? 11 : 8), &h) < 0) continue;
    for (k = 0; k < teems_ncomp; k++) {
      comp_def *cp = &teems_comps[k];
      c->stmt = line;
      if (isvar && (int)strlen(cp->varname) == h.namelen && strncmp(cp->varname, h.name, h.namelen) == 0) {
        if (cp_argsets(c, &h, argsets, &nargs) < 0) { fclose(f); return -1; }
        if (nargs != cp->nquant) {
          printf("Error: Complementarity %s has %d quantifiers but its variable %s has %d arguments (manual 11.14)\n", cp->name, (int)cp->nquant, cp->varname, (int)nargs);
          fclose(f);
          return -1;
        }
        memcpy(cp->xset, argsets, sizeof(argsets));
      }
      if (cp->lower_kind >= 2 && ((isvar && cp->lower_kind == 2) || (iscoef && cp->lower_kind == 3)) &&
          (int)strlen(cp->lower_name) == h.namelen && strncmp(cp->lower_name, h.name, h.namelen) == 0) {
        if (cp_argsets(c, &h, argsets, &nargs) < 0) { fclose(f); return -1; }
        if (nargs != cp->nquant) {
          printf("Error: Complementarity %s has %d quantifiers but its lower bound %s has %d arguments (manual 11.14)\n", cp->name, (int)cp->nquant, cp->lower_name, (int)nargs);
          fclose(f);
          return -1;
        }
        memcpy(cp->lset, argsets, sizeof(argsets));
      }
      if (cp->upper_kind >= 2 && ((isvar && cp->upper_kind == 2) || (iscoef && cp->upper_kind == 3)) &&
          (int)strlen(cp->upper_name) == h.namelen && strncmp(cp->upper_name, h.name, h.namelen) == 0) {
        if (cp_argsets(c, &h, argsets, &nargs) < 0) { fclose(f); return -1; }
        if (nargs != cp->nquant) {
          printf("Error: Complementarity %s has %d quantifiers but its upper bound %s has %d arguments (manual 11.14)\n", cp->name, (int)cp->nquant, cp->upper_name, (int)nargs);
          fclose(f);
          return -1;
        }
        memcpy(cp->uset, argsets, sizeof(argsets));
      }
    }
  }
  fclose(f);
  return 0;
}

/* emit the derived statements for one complementarity */
static int cp_emit_derived(FILE *fout, comp_def *cp, const char *expr) {
  char quants[TABREADLINE], idxs[TABREADLINE];
  dim_t k;
  quants[0] = '\0';
  idxs[0] = '\0';
  for (k = 0; k < cp->nquant; k++) {
    char one[NAMESIZE * 2 + 16];
    snprintf(one, sizeof(one), "(all,%s,%s)", cp->qidx[k], cp->qset[k]);
    if (lv_strcat_b(quants, one, sizeof(quants)) < 0) return -1;
    if (lv_strcat_b(idxs, k == 0 ? "(" : ",", sizeof(idxs)) < 0) return -1;
    if (lv_strcat_b(idxs, cp->qidx[k], sizeof(idxs)) < 0) return -1;
  }
  if (cp->nquant > 0 && lv_strcat_b(idxs, ")", sizeof(idxs)) < 0) return -1;
  /* emitted lines must survive the downstream TABREADLINE readers */
  if (strlen(quants) + 2 * strlen(idxs) + strlen(expr) + strlen(cp->name) + 64 >= TABREADLINE - 2) return -1;
  fprintf(fout, "variable (change,levels) %s %s@e%s ;\n", quants, cp->name, idxs);
  fprintf(fout, "formula (initial) %s %s@e%s = %s ;\n", quants, cp->name, idxs, expr);
  fprintf(fout, "equation (levels) e_%s@e %s %s@e%s = %s ;\n", cp->name, quants, cp->name, idxs, expr);
  fprintf(fout, "variable (change) %s %s@d%s ;\n", quants, cp->name, idxs);
  if (cp->lower_kind == 2) {
    fprintf(fout, "variable (change,levels) %s %s@l%s ;\n", quants, cp->name, idxs);
    fprintf(fout, "formula (initial) %s %s@l%s = %s%s - %s%s ;\n", quants, cp->name, idxs, cp->varname, idxs, cp->lower_name, idxs);
    fprintf(fout, "equation (levels) e_%s@l %s %s@l%s = %s%s - %s%s ;\n", cp->name, quants, cp->name, idxs, cp->varname, idxs, cp->lower_name, idxs);
  }
  if (cp->upper_kind == 2) {
    fprintf(fout, "variable (change,levels) %s %s@u%s ;\n", quants, cp->name, idxs);
    fprintf(fout, "formula (initial) %s %s@u%s = %s%s - %s%s ;\n", quants, cp->name, idxs, cp->varname, idxs, cp->upper_name, idxs);
    fprintf(fout, "equation (levels) e_%s@u %s %s@u%s = %s%s - %s%s ;\n", cp->name, quants, cp->name, idxs, cp->varname, idxs, cp->upper_name, idxs);
  }
  return 0;
}

int tab_complementarity_transform(char *fname) {
  FILE *f, *fout;
  char line[TABREADLINE], stmt[TABREADLINE], tmpname[TABREADLINE];
  char (*exprs)[TABREADLINE] = NULL;
  lv_ctx *c;
  bool anylv = false;
  int rc = 0;
  dim_t ncomp = 0, ci = 0;
  f = fopen(fname, "r");
  if (f == NULL) { printf("Error: cannot open %s\n", fname); return -1; }
  while (fgets(line, TABREADLINE, f))
    if (strncmp(line, "complementarity", 15) == 0) ncomp++;
  fclose(f);
  if (ncomp == 0) return 0;
  c = (lv_ctx *)calloc(1, sizeof(lv_ctx));
  if (c == NULL) { printf("Error: out of memory in tab_complementarity_transform\n"); return -1; }
  if (lv_scan(c, fname, &anylv) < 0) { free(c); return -1; }
  teems_comps = (comp_def *)calloc(ncomp, sizeof(comp_def));
  exprs = calloc(ncomp, sizeof(*exprs));
  if (teems_comps == NULL || exprs == NULL) { printf("Error: out of memory in tab_complementarity_transform\n"); free(c); free(exprs); return -1; }
  teems_ncomp = 0;
  f = fopen(fname, "r");
  if (f == NULL) { printf("Error: cannot open %s\n", fname); free(c); free(exprs); return -1; }
  while (fgets(line, TABREADLINE, f)) {
    size_t sl = strlen(line);
    while (sl > 0 && (line[sl - 1] == '\n' || line[sl - 1] == '\r')) line[--sl] = '\0';
    if (strncmp(line, "complementarity", 15) != 0) continue;
    if (cp_parse_stmt(c, line, &teems_comps[teems_ncomp], exprs[teems_ncomp], TABREADLINE) < 0) { rc = -1; break; }
    teems_ncomp++;
  }
  fclose(f);
  if (rc == 0) rc = cp_fill_decl_sets(c, fname);
  if (rc == 0) {
    /* rewrite: consume complementarity statements, emit the derived
       statements in their place (del_comp@ once, before the first) */
    f = fopen(fname, "r");
    strcpy(tmpname, fname);
    strcat(tmpname, "_cp");
    fout = f == NULL ? NULL : fopen(tmpname, "w");
    if (f == NULL || fout == NULL) {
      printf("Error: cannot open %s\n", f == NULL ? fname : tmpname);
      if (f != NULL) fclose(f);
      rc = -1;
    } else {
      while (fgets(line, TABREADLINE, f)) {
        if (strncmp(line, "complementarity", 15) == 0) {
          if (ci == 0) fprintf(fout, "variable (change) del_comp@ ;\n");
          strcpy(stmt, line);
          c->stmt = stmt;
          if (cp_emit_derived(fout, &teems_comps[ci], exprs[ci]) < 0) { rc = lv_err(c, "derived statement too long"); break; }
          ci++;
          continue;
        }
        fputs(line, fout);
      }
      fclose(f);
      fclose(fout);
      if (rc == 0) {
        if (rename(tmpname, fname) != 0) {
          printf("Error: cannot rename %s\n", tmpname);
          rc = -1;
        }
      } else remove(tmpname);
    }
  }
  free(exprs);
  free(c);
  if (rc != 0) { free(teems_comps); teems_comps = NULL; teems_ncomp = 0; }
  return rc;
}

/* 11.14 points 2-3: each statement quantifier set must be equal to or
   an ordered subset of the corresponding set of X (and of each
   variable/coefficient bound). Runs once set elements exist. */
static int cp_ordered_subset(const char *small, const char *big, set_def *sets, dim_t nset, set_element *se, const char *compname, const char *ofwhat) {
  dim_t si = nset, bi = nset, i, j;
  for (i = 0; i < nset; i++) {
    if (strcmp(sets[i].setname, small) == 0) si = i;
    if (strcmp(sets[i].setname, big) == 0) bi = i;
  }
  if (si == nset || bi == nset) {
    printf("Error: Complementarity %s references undeclared set %s\n", compname, si == nset ? small : big);
    return -1;
  }
  if (si == bi) return 0;
  j = 0;
  for (i = 0; i < (dim_t)sets[si].size; i++) {
    while (j < (dim_t)sets[bi].size && strcmp(se[sets[si].offset + i].setele, se[sets[bi].offset + j].setele) != 0) j++;
    if (j == (dim_t)sets[bi].size) {
      printf("Error: Complementarity %s: quantifier set %s is not an equal or same-ordered subset of %s (%s) (manual 11.14)\n", compname, small, big, ofwhat);
      return -1;
    }
    j++;
  }
  return 0;
}

int complementarities_validate(set_def *sets, dim_t nset, set_element *set_elems) {
  dim_t k, q;
  for (k = 0; k < teems_ncomp; k++) {
    comp_def *cp = &teems_comps[k];
    for (q = 0; q < cp->nquant; q++) {
      if (cp_ordered_subset(cp->qset[q], cp->xset[q], sets, nset, set_elems, cp->name, cp->varname) == -1) return -1;
      if (cp->lower_kind >= 2 &&
          cp_ordered_subset(cp->qset[q], cp->lset[q], sets, nset, set_elems, cp->name, cp->lower_name) == -1) return -1;
      if (cp->upper_kind >= 2 &&
          cp_ordered_subset(cp->qset[q], cp->uset[q], sets, nset, set_elems, cp->name, cp->upper_name) == -1) return -1;
    }
  }
  return 0;
}

/* closure integration (51.7.2 (b)-(e) + 11.14.1, C1 subset): the
   dummy variables and del_comp@ are auto-exogenized (the user cannot
   mention them -- closure_var_find fatals on '@' names); comp@e/@l/@u
   stay endogenous by omission. The complementarity variable must not
   be backsolved (11.14.1) and -- until the C2 state machinery -- must
   be fully exogenous (the E_$comp equation row does not exist yet). */
int comp_closure_check(closure_entry *closure_vals, array_def *vars, offset_t nvar, offset_t *nexo) {
  offset_t i, j;
  dim_t k;
  if (teems_ncomp == 0) return 0;
  for (i = 0; i < nvar; i++) {
    size_t ln = strlen(vars[i].cofname);
    offset_t nele = vars[i].nelem == 0 ? 1 : vars[i].nelem; /* scalars: nelem 0, one slot */
    bool dummy;
    if (strchr(vars[i].cofname, '@') == NULL) continue;
    dummy = (ln >= 2 && vars[i].cofname[ln - 2] == '@' && vars[i].cofname[ln - 1] == 'd') ||
            strcmp(vars[i].cofname, "del_comp@") == 0;
    for (j = 0; j < nele; j++) {
      if (closure_vals[vars[i].offset + j].is_backsolved) {
        printf("Error: derived complementarity variable %s cannot be backsolved\n", vars[i].cofname);
        return -1;
      }
      if (dummy && !closure_vals[vars[i].offset + j].is_exogenous) {
        closure_vals[vars[i].offset + j].is_exogenous = true;
        (*nexo)++;
      }
    }
  }
  for (k = 0; k < teems_ncomp; k++) {
    comp_def *cp = &teems_comps[k];
    offset_t nele;
    for (i = 0; i < nvar; i++) if (strcmp(vars[i].cofname, cp->varname) == 0) break;
    if (i == nvar) {
      printf("Error: Complementarity %s: variable %s not found after reading declarations\n", cp->name, cp->varname);
      return -1;
    }
    nele = vars[i].nelem == 0 ? 1 : vars[i].nelem;
    for (j = 0; j < nele; j++) {
      if (closure_vals[vars[i].offset + j].is_backsolved) {
        printf("Error: the Complementarity variable %s must not be backsolved (manual 11.14.1)\n", cp->varname);
        return -1;
      }
      if (!closure_vals[vars[i].offset + j].is_exogenous) {
        printf("Error: Complementarity %s: the approximate-run state machinery is not implemented yet (C2); a complementarity is currently solvable only with its variable %s fully exogenous\n", cp->name, cp->varname);
        return -1;
      }
    }
  }
  return 0;
}

int tab_levels_transform(char *fname) {
  lv_ctx *c;
  FILE *f, *fout;
  char line[TABREADLINE], stmt[TABREADLINE], tmpname[TABREADLINE], defval[NAMESIZE];
  bool any = false;
  lv_head h;
  int rc = 0;
  c = (lv_ctx *)calloc(1, sizeof(lv_ctx));
  if (c == NULL) { printf("Error: out of memory in tab_levels_transform\n"); return -1; }
  if (lv_scan(c, fname, &any) < 0) { free(c); return -1; }
  if (!any) { free(c); return 0; }
  f = fopen(fname, "r");
  if (f == NULL) { printf("Error: cannot open %s\n", fname); free(c); return -1; }
  strcpy(tmpname, fname);
  strcat(tmpname, "_lv");
  fout = fopen(tmpname, "w");
  if (fout == NULL) { printf("Error: cannot open %s\n", tmpname); fclose(f); free(c); return -1; }
  while (fgets(line, TABREADLINE, f)) {
    /* keep a pristine copy for messages and verbatim passthrough:
       the head/emit parsers only read, but tab_default_value and the
       expression cursor bound spans inside `stmt` */
    strcpy(stmt, line);
    {
      size_t sl = strlen(stmt);
      while (sl > 0 && (stmt[sl - 1] == '\n' || stmt[sl - 1] == '\r')) stmt[--sl] = '\0';
    }
    c->stmt = stmt;
    /* C1a: rewrite bare value references of renamed pairs before any
       further handling -- value-reference contexts only (declaration
       statements keep declared names for variables_read) */
    if (c->nrenamed > 0 &&
        (strncmp(stmt, "formula", 7) == 0 || strncmp(stmt, "equation", 8) == 0 ||
         strncmp(stmt, "update", 6) == 0 || strncmp(stmt, "read", 4) == 0 ||
         strncmp(stmt, "write", 5) == 0 || strncmp(stmt, "display", 7) == 0 ||
         strncmp(stmt, "assertion", 9) == 0)) {
      if (lv_rename_line(c, stmt, sizeof(stmt)) < 0) { rc = lv_err(c, "statement too long after levels value rename"); break; }
    }
    if (strncmp(stmt, "variable", 8) == 0 && !tab_default_value(stmt, defval)) {
      if (lv_parse_head(stmt + 8, &h) == 0 && lv_find_lv_decl(c, h.name, h.namelen) >= 0) {
        int k = lv_find_lv_decl(c, h.name, h.namelen);
        /* declaration passes through verbatim (variables_read
           re-derives level_par/change_real from the same statement);
           append the value coefficient and its update (9.2.2). The
           name+args ride verbatim so the argument order stays exactly
           as declared (it may differ from quantifier order); the
           coefficient carries the pair's VALUE name (== declared
           unless gen_lv-renamed, C1a), the update RHS the declared
           column reference. */
        fprintf(fout, "%s\n", stmt);
        fprintf(fout, "coefficient (non_parameter) %.*s %s%.*s ;\n",
                h.quantslen, h.quants, c->lv[k].valname, h.argslen, h.args);
        fprintf(fout, "update %s%.*s %s%.*s = p_%.*s%.*s ;\n",
                c->lv[k].change ? "(change) " : "",
                h.quantslen, h.quants, c->lv[k].valname, h.argslen, h.args,
                h.namelen, h.name, h.argslen, h.args);
        continue;
      }
      fprintf(fout, "%s\n", stmt);
      continue;
    }
    if (strncmp(stmt, "formula", 7) == 0) {
      char *amp = strstr(stmt, "&equation");
      int ampskip = 9;
      if (amp == NULL) { amp = strstr(stmt, "& equation"); ampskip = 10; }
      if (amp != NULL) {
        /* formula [ (initial) ] & equation [ (levels) ] NAME [#lbl#]
           quants LHS = RHS ; (10.9.1) */
        const char *q = stmt + 7;
        const char *tail;
        char eqline[TABREADLINE];
        while (*q == ' ') q++;
        if (strncmp(q, "(initial)", 9) == 0) { q += 9; while (*q == ' ') q++; }
        if (q != amp) { rc = lv_err(c, "malformed Formula&Equation head"); break; }
        tail = amp + ampskip;
        while (*tail == ' ') tail++;
        if (strncmp(tail, "(levels)", 8) == 0) { tail += 8; while (*tail == ' ') tail++; }
        /* tail = NAME [#lbl#] quants LHS=RHS; -- the formula drops the
           name/label, the equation keeps them */
        {
          const char *nm = tail, *body;
          int nmlen;
          while (lv_isnamec(*tail)) tail++;
          nmlen = (int)(tail - nm);
          if (nmlen == 0) { rc = lv_err(c, "missing Formula&Equation name"); break; }
          body = tail;
          while (*body == ' ') body++;
          if (*body == '#') {
            body++;
            while (*body != '\0' && *body != '#') body++;
            if (*body != '#') { rc = lv_err(c, "unterminated # label #"); break; }
            body++;
            while (*body == ' ') body++;
          }
          fprintf(fout, "formula (initial) %s\n", body);
          snprintf(eqline, sizeof(eqline), "equation (levels) %.*s %s", nmlen, nm, body);
          c->stmt = eqline;
          if (lv_emit_linearized(c, eqline + 8, fout) < 0) { rc = -1; break; }
        }
        continue;
      }
      fprintf(fout, "%s\n", stmt);
      continue;
    }
    if (strncmp(stmt, "equation", 8) == 0) {
      const char *q = stmt + 8;
      while (*q == ' ') q++;
      if (strncmp(q, "(levels", 7) == 0) {
        if (lv_emit_linearized(c, stmt + 8, fout) < 0) { rc = -1; break; }
        continue;
      }
      fprintf(fout, "%s\n", stmt);
      continue;
    }
    fprintf(fout, "%s\n", stmt);
  }
  fclose(f);
  fclose(fout);
  if (rc == 0) {
    if (rename(tmpname, fname) != 0) {
      printf("Error: cannot rename %s\n", tmpname);
      rc = -1;
    }
  } else remove(tmpname);
  /* C1a: the _ps companion references pair coefficients by bare name
     (PostSim binds levels values); apply the same value rename there */
  if (rc == 0 && c->nrenamed > 0) {
    char psname[TABREADLINE];
    strcpy(psname, fname);
    strcat(psname, "_ps");
    f = fopen(psname, "r");
    if (f != NULL) {
      strcpy(tmpname, psname);
      strcat(tmpname, "_lv");
      fout = fopen(tmpname, "w");
      if (fout == NULL) { printf("Error: cannot open %s\n", tmpname); fclose(f); free(c); return -1; }
      while (fgets(line, TABREADLINE, f)) {
        c->stmt = line;
        if (lv_rename_line(c, line, sizeof(line)) < 0) { rc = lv_err(c, "statement too long after levels value rename"); break; }
        fputs(line, fout);
      }
      fclose(f);
      fclose(fout);
      if (rc == 0) {
        if (rename(tmpname, psname) != 0) {
          printf("Error: cannot rename %s\n", tmpname);
          rc = -1;
        }
      } else remove(tmpname);
    }
  }
  free(c);
  return rc;
}
