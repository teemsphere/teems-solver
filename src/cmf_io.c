#include <teems_solver.h>

/* bounded append/copy for tab_preprocess's growable line buffers: the
   auto-index transform accumulates one subset-declaration block per
   quoted set element into readline/readline1 and rebuilds the statement
   in line2, which can grow past the fixed buffers on a pathological
   statement (many quoted elements). These return -1 (caller aborts)
   instead of overflowing; a no-op for valid models, which stay under
   the buffer. */
static int cmf_strcat_bounded(char *dst, const char *src, size_t cap) {
  size_t dl = strlen(dst), sl = strlen(src);
  if (dl + sl + 1 > cap) return -1;
  memcpy(dst + dl, src, sl + 1);
  return 0;
}
static int cmf_strcpy_bounded(char *dst, const char *src, size_t cap) {
  size_t sl = strlen(src);
  if (sl + 1 > cap) return -1;
  memcpy(dst, src, sl + 1);
  return 0;
}

/* Normalize the `:`-condition segments of a statement line (manual
   11.4.11) ahead of the singleton-subset transform: the GEMPACK word
   comparisons EQ/NE/GT/LT/GE/LE become = <> > < >= <= (while blanks
   still delimit them: space-stripped, "SC GT 1" is indistinguishable
   from a name), and a quoted element RHS
   immediately after `=` is unquoted and lowercased (quoted strings
   elsewhere stay verbatim for the transform). A segment runs from a
   ':' inside brackets to the next ',' or closing bracket at the same
   depth; brace and paren sum forms both count. In-place (the result
   never grows). */
static int cond_is_namec(char ch) {
  return (ch>='a'&&ch<='z')||(ch>='A'&&ch<='Z')||(ch>='0'&&ch<='9')||ch=='_'||ch=='@'||ch==MAPMARK;
}
/* Formula qualifier lists (manual 10.8, 11.11): "(initial, always,
   by_elements, ...)" and the "write updated value to file F header "H""
   form are split into one group per qualifier, which is the form every
   later reader recognises; the updated-value write request is dropped
   with a warning (its value is in the coefficient dump, sol.cof/.cbin,
   like every coefficient's). A compound group used to stay in place: its
   '(' counted as a quantifier and formulas_execute read past the
   statement (SEGV), and its quoted header was lowered as an element
   literal (potential-models segfaults.md #4, G-S4). */
static void formula_qualifiers_normalize(char *line) {
  char out[TABREADLINE],item[TABREADLINE];
  char *p=line;
  size_t o=0;
  int changed=0;
  while (*p==' ') p++;
  if (str_find_ci(p,"formula")!=0) return;
  p+=7;
  if (*p!=' '&&*p!='(') return;
  memcpy(out,line,p-line);
  o=p-line;
  for (;;) {
    char *g=p,*e;
    int d=0,inq=0;
    while (*g==' ') g++;
    if (*g!='('||str_find_ci(g,"(all,")==0) break;
    for (e=g; *e!='\0'; e++) {
      if (*e=='\"') inq=!inq;
      else if (!inq&&*e=='(') d++;
      else if (!inq&&*e==')') { d--; if (d==0) break; }
    }
    if (*e!=')') return;
    {
      char *c=g+1,*st=g+1;
      int dq=0,ni=0;
      inq=0;
      for (;; c++) {
        if (*c=='\"') inq=!inq;
        else if (!inq&&*c=='(') dq++;
        else if (!inq&&*c==')'&&dq>0) dq--;
        if (c==e||(!inq&&dq==0&&*c==',')) {
          size_t n=c-st;
          while (n>0&&*st==' ') { st++; n--; }
          while (n>0&&st[n-1]==' ') n--;
          if (n>=sizeof(item)) return;
          memcpy(item,st,n);
          item[n]='\0';
          ni++;
          if (str_find_ci(item,"write")==0) {
            printf("Warning: Formula qualifier \"%s\" is accepted but no separate updated-data file is written (the value is in the coefficient dump)\n",item);
            changed=1;
          } else if (n>0) {
            if (o+n+3>=sizeof(out)) return;
            out[o++]=' ';
            out[o++]='(';
            memcpy(out+o,item,n);
            o+=n;
            out[o++]=')';
          }
          if (c==e) break;
          st=c+1;
        }
      }
      if (ni>1) changed=1;
    }
    p=e+1;
  }
  if (!changed) return;
  if (o+strlen(p)+1>=sizeof(out)) return;
  if (*p!=' ') out[o++]=' ';
  strcpy(out+o,p);
  strcpy(line,out);
}

static void cond_segment_normalize(char *line) {
  char buf[TABREADLINE];
  int depth=0,incond=0,conddepth=0,changed=0;
  char lastns='\0';
  size_t bi=0,i;
  for (i=0; line[i]!='\0'&&bi+1<sizeof(buf); i++) {
    char ch=line[i];
    if (!incond) {
      if (ch=='('||ch=='{'||ch=='[') depth++;
      else if (ch==')'||ch=='}'||ch==']') depth--;
      else if (ch==':'&&depth>0) {
        incond=1;
        conddepth=depth;
        lastns='\0';
      }
      buf[bi++]=ch;
      continue;
    }
    if (ch=='('||ch=='{'||ch=='[') depth++;
    else if (ch==')'||ch=='}'||ch==']') {
      if (depth==conddepth) incond=0;
      depth--;
    }
    else if (ch==','&&depth==conddepth) incond=0;
    if (!incond) {
      buf[bi++]=ch;
      continue;
    }
    if (i>0&&!cond_is_namec(line[i-1])&&line[i+1]!='\0'&&!cond_is_namec(line[i+2])&&
        line[i+2]!='('&&line[i+2]!='['&&line[i+2]!='{') {
      static const char *wops[6][2]={{"eq","="},{"ne","<>"},{"gt",">"},{"lt","<"},{"ge",">="},{"le","<="}};
      int w;
      for (w=0; w<6; w++)
        if (tolower((int)ch)==wops[w][0][0]&&tolower((int)line[i+1])==wops[w][0][1]) break;
      if (w<6) {
        const char *s=wops[w][1];
        while (*s!='\0') buf[bi++]=*s++;
        lastns=wops[w][1][strlen(wops[w][1])-1];
        i++;
        changed=1;
        continue;
      }
    }
    if (ch=='\"'&&lastns=='=') {
      for (i++; line[i]!='\0'&&line[i]!='\"'&&bi+1<sizeof(buf); i++)
        buf[bi++]=(char)tolower((int)line[i]);
      if (line[i]!='\"') { changed=0; break; } /* unterminated: leave the line (and the fatal) to the readers */
      changed=1;
      lastns='\0';
      continue;
    }
    if (ch=='\"') {
      /* any other quoted element is copied as written: "NE" is an
         element, not the word operator */
      buf[bi++]=ch;
      for (i++; line[i]!='\0'&&line[i]!='\"'&&bi+2<sizeof(buf); i++) buf[bi++]=line[i];
      if (line[i]!='\"') { changed=0; break; }
      buf[bi++]=line[i];
      lastns='\"';
      continue;
    }
    buf[bi++]=ch;
    if (ch!=' '&&ch!='\t') lastns=ch;
  }
  if (changed&&line[i]=='\0') {
    buf[bi]='\0';
    strcpy(line,buf);
  }
}

int cmf_count_files(char *fname,char *comsyntax) {
  FILE * filehandle;
  char line[TABREADLINE]="\0";
  int j=0;
  filehandle = teems_fopen(fname,"r");
  if(filehandle==NULL){
    errmsg("Error: cannot open %s\n",fname);
    return -1;
  }
  while (tab_next_statement(comsyntax,filehandle,line,TABREADLINE)) {
    logmsg(2,"Com %s file %s\n",comsyntax,line);
    j++;
  }
  fclose(filehandle);
  return j;
}

int datafile_read_header_info(char *varname, char *filename,dim_t *vsize, char *longname,dim_t *d1) {
  FILE * filehandle;
  char line[TABREADLINE+1],linecopy[TABREADLINE+1];
  dim_t nlength=0,nlength1=0,vsizein=0,din1=0;
  int succ=0;
  char *readitem=NULL;
  while (varname[nlength] != '\0') nlength++;
  filehandle = teems_fopen(filename,"r");
  if(filehandle==NULL){
    errmsg("Error: cannot open %s\n",filename);
    return -1;
  }

  while (fgets(line,TABREADLINE,filehandle)) {
    strcpy(linecopy,line);
    readitem = strtok(line,"\"");
    readitem = strtok(NULL,"\"");
    if (readitem != NULL) {
      nlength1=0;
      while (readitem[nlength1] != '\0'&&readitem[nlength1] != ' ') nlength1++;
      if (nlength1==nlength&&str_ncmp_ci(readitem,varname,nlength) == 0) {
        succ=1;
        readitem = strtok(line," ");
        if (readitem==NULL) { fclose(filehandle); return -1; }
        din1=atoi(readitem);//strtol(readitem,NULL,10);//sscanf(readitem, "%d", &din1);
        readitem = strtok(NULL," ");
        readitem = strtok(NULL," ");
        readitem = strtok(NULL," ");
        if (readitem==NULL) { fclose(filehandle); return -1; }
        vsizein=atoi(readitem);//sscanf(readitem, "%d", &vsizein);
        readitem = strtok(linecopy,"\"");
        readitem = strtok(NULL,"\"");
        readitem = strtok(NULL,"\"");
        readitem = strtok(NULL,"\"");
        if (readitem==NULL) { fclose(filehandle); return -1; }
        strcpy(longname,readitem);
      }
    }
  }
  fclose(filehandle);
  *d1=din1;
  *vsize=vsizein;
  if(succ==0)errmsg("Error: header \"%s\" not found in %s\n",varname,filename);
  return succ;
}

int datafile_read_labels(char *varname, char *filename,dim_t d1, datafile_labels *record) {
  FILE * dfile;
  char line[DATREADLINE],header[NAMESIZE],varnamecpy[NAMESIZE];
  dim_t nlength=0,nhead=0,reccount = 0,count1=0;
  char *readitem=NULL;
  /* callers pass a TAB token here (tab_parse.c mapping reader), not only
     a short set header, so this copy needs the same bound as the rest */
  if (cmf_strcpy_bounded(varnamecpy,varname,sizeof(varnamecpy))) {
    errmsg("Error: header name %s exceeds %d characters\n",varname,NAMESIZE-1);
    return -1;
  }
  str_delete_char(varnamecpy,' ');
  while (varnamecpy[nlength] != '\0') nlength++;
  /* nothing to fill: `record` is calloc'd from d1, so d1 0 leaves a
     zero-size allocation the label loop below wrote into (fuzz batch 14) */
  if (d1<=0) return 0;
  dfile = teems_fopen(filename,"r");
  if (dfile==NULL) {
    errmsg("Error: cannot open %s\n",filename);
    return -1;
  }

  while (fgets(line,DATREADLINE,dfile)) {
    readitem = strtok(line,"\"");
    readitem = strtok(NULL,"\"");
    if (readitem != NULL) {
      /* readitem is a quoted token out of a DATREADLINE line; header is
         NAMESIZE, and a long one overran it into the neighbouring
         buffers (fuzz batch 14). A header that cannot fit cannot match
         the name being looked up, so skipping the line is correct. */
      if (cmf_strcpy_bounded(header,readitem,sizeof(header))) continue;
      str_delete_char(header,' ');
      nhead=0;
      while (header[nhead] != '\0') nhead++;
      if(nhead<nlength)nhead=nlength;
      if (str_ncmp_ci(readitem,varnamecpy,nhead) == 0) {
        while (fgets(line,DATREADLINE,dfile)) {
          str_delete_char(line,'\r');
          while (str_replace_all(line,"  ", " "));
          if (line[0]=='\n') count1=1;
          if (count1!=1) {
            readitem = strtok(line,"\n");
            if (readitem != NULL) {
              /* the bound was checked only AFTER the write below, so a
                 second block matching the same header re-entered with
                 reccount already at d1 (fuzz batch 14) */
              if (reccount>=d1) break;
              record[reccount].dim1=reccount;
              /* element labels are NAMESIZE; a long data line overran
                 the record array (fuzz batch 14) */
              if (cmf_strcpy_bounded(record[reccount].ch,readitem,sizeof(record[reccount].ch))) {
                errmsg("Error: element label for header %s exceeds %d characters in %s\n",varname,NAMESIZE-1,filename);
                fclose(dfile);
                MPI_Abort(PETSC_COMM_WORLD,1);
                return -1;
              }
              reccount++;
              if (reccount>=d1) break;
            }
          }

        }
      }

    }
  }
  fclose(dfile);
  return 0;
}

/* the manifest's optional subtotals statement (Tier B3): subtotals
   "<path>"; names the file of shock groups. Returns 1 and the path in
   out when present, 0 when absent; a second statement is fatal. */
int cmf_subtotals_file(char *filename, char *out) {
  FILE *filehandle;
  char line[TABREADLINE],commsyntax[NAMESIZE]="subtotals",*readitem;
  int n=0;
  out[0]='\0';
  filehandle = teems_fopen(filename,"r");
  if(filehandle==NULL)return 0;
  while (tab_next_statement(commsyntax,filehandle,line,TABREADLINE)) {
    readitem = strtok(line,"\"");
    readitem = strtok(NULL,"\"");
    if(readitem==NULL) {
      errmsg("Error: the subtotals statement of the manifest (.cmf) file names no quoted file\n");
      fclose(filehandle);
      MPI_Abort(PETSC_COMM_WORLD,1);
    }
    if(++n>1) {
      errmsg("Error: the manifest (.cmf) file has more than one subtotals statement; put every subtotal in one file\n");
      fclose(filehandle);
      MPI_Abort(PETSC_COMM_WORLD,1);
    }
    snprintf(out,TABREADLINE,"%s",readitem);
  }
  fclose(filehandle);
  return n;
}

int cmf_read(char *filename, int niodata, cmf_file_entry *iodata, char *tabfile, char *closure, char *shock) {
  FILE * filehandle;
  char line[TABREADLINE],*readitem,commsyntax[NAMESIZE];
  int j=0,k;
  strcpy(commsyntax,"iodata");
  filehandle = teems_fopen(filename,"r");
  if(filehandle==NULL){
    errmsg("Error: cannot open %s\n",filename);
    return -1;
  }
  while (tab_next_statement(commsyntax,filehandle,line,TABREADLINE)) {
    readitem = strtok(line,"\"");
    readitem = strtok(NULL,"\"");
    k=0;
    while (readitem[k]!= '\0') {
      readitem[k]=tolower((int)readitem[k]);
      k++;
    }
    strcpy(iodata[j].logname,readitem);
    k=0;
    readitem = strtok(NULL,"\"");
    readitem = strtok(NULL,"\"");
    strcpy(iodata[j].filname,readitem);
    j++;
  }
  fclose(filehandle);

  strcpy(commsyntax,"outdata");
  filehandle = teems_fopen(filename,"r");
  while (tab_next_statement(commsyntax,filehandle,line,TABREADLINE)) {
    readitem = strtok(line,"\"");
    readitem = strtok(NULL,"\"");
    k=0;
    while (readitem[k]!= '\0') {
      readitem[k]=tolower((int)readitem[k]);
      k++;
    }
    strcpy(iodata[j].logname,readitem);
    k=0;
    readitem = strtok(NULL,"\"");
    readitem = strtok(NULL,"\"");
    strcpy(iodata[j].filname,readitem);
    j++;
  }
  fclose(filehandle);

  strcpy(commsyntax,"soldata");
  filehandle = teems_fopen(filename,"r");
  while (tab_next_statement(commsyntax,filehandle,line,TABREADLINE)) {
    readitem = strtok(line,"\"");
    readitem = strtok(NULL,"\"");
    k=0;
    while (readitem[k]!= '\0') {
      readitem[k]=tolower((int)readitem[k]);
      k++;
    }
    strcpy(iodata[j].logname,readitem);
    k=0;
    readitem = strtok(NULL,"\"");
    readitem = strtok(NULL,"\"");
    strcpy(iodata[j].filname,readitem);
    j++;
  }
  fclose(filehandle);

  filehandle = teems_fopen(filename,"r");
  strcpy(commsyntax,"tabfile");
  while (tab_next_statement(commsyntax,filehandle,line,TABREADLINE)) {
    readitem = strtok(line,"\"");
    readitem = strtok(NULL,"\"");
    strcpy(tabfile,readitem);
  }
  fclose(filehandle);
  filehandle = teems_fopen(filename,"r");
  strcpy(commsyntax,"closure");
  while (tab_next_statement(commsyntax,filehandle,line,TABREADLINE)) {
    readitem = strtok(line,"\"");
    readitem = strtok(NULL,"\"");
    strcpy(closure,readitem);
  }
  fclose(filehandle);
  filehandle = teems_fopen(filename,"r");
  strcpy(commsyntax,"shock");
  while (tab_next_statement(commsyntax,filehandle,line,TABREADLINE)) {
    readitem = strtok(line,"\"");
    readitem = strtok(NULL,"\"");
    strcpy(shock,readitem);
  }
  fclose(filehandle);
  return 1;
}

/* token-boundary finds throughout: "wsum(b," must not scan as a
   repeated-index sum over b (the *sum-name class the bordered-map
   kit's condsum leg caught -- the rename here rewrote a variable's
   argument list) */
/* One rename pass of sum_dedup_indices: replace every occurrence of `a`
   in line1 (which points into a TABREADLINE buffer with `cap` bytes left)
   by `b`.  Forward-scanning: the old restart-from-the-start loop spun
   forever when the replacement re-created the pattern at its boundary
   (an empty sum index makes the pattern two spaces and b ends in one;
   fuzz batch 13 hang class), and could grow the line past its buffer. */
static void dedup_subst(char *line1, const char *a, const char *b, size_t cap) {
  if (a[0]=='\0') {
    errmsg("Error: cannot rename a sum index: malformed sum(<index>,<set>,...) in a formula or equation\n");
    MPI_Abort(PETSC_COMM_WORLD,1);
  }
  if (str_subst_all_bounded(line1,a,b,cap)!=0) {
    errmsg("Error: renaming sum index %s does not fit the statement buffer (%d characters)\n",a,(int)TABREADLINE);
    MPI_Abort(PETSC_COMM_WORLD,1);
  }
}

/* Write the renamed sum span `b` back over its original text `a` (first
   occurrence) in the whole statement `line` (capacity cap).  The rename
   loops above bound the span against the buffer that holds it, but the
   span was cut out of the statement by sum_extract, so the statement's
   tail after the span is not in that bound: the old unbounded
   str_replace_all write-back overflowed `line` by the tail's length
   (fuzz batch 13 stack overflow, tab_parse.c str_replace_all). */
static void dedup_subst_back(char *line, const char *a, const char *b, size_t cap) {
  size_t alen=strlen(a),blen=strlen(b),len=strlen(line);
  char *p=strstr(line,a);
  if (p==NULL || alen==0) return;
  if (len-alen+blen>=cap) {
    errmsg("Error: renaming a repeated sum index does not fit the statement buffer (%d characters): %s\n",(int)TABREADLINE,a);
    MPI_Abort(PETSC_COMM_WORLD,1);
  }
  memmove(p+blen,p+alen,len-(size_t)(p-line)-alen+1);
  memcpy(p,b,blen);
}

/* a quoted element inside a "# label #" is label text, not an element
   argument (manual 11.2.3): literal lowering scanned it and synthesized
   a subset from it (TERM "# Check HOUPUR_H = PUR_S("Hou") #" became
   S(hou) subset of COM). Mask those quotes with the \001 marker the
   lowering already restores on output. */
static void label_quotes_mask(char *s) {
  int in_label=0;
  for (; *s!='\0'; s++) {
    if (*s=='#') in_label=!in_label;
    else if (in_label&&*s=='\"') *s='\001';
  }
}

/* the delimiter pairs above miss a bare index beside an operator, which
   is how a sum condition reads it (sum(r,reg: r<>d, ...)): the index
   was renamed everywhere but there, so the condition tested an unbound
   name (TERM E_pdelimp). Rename every remaining whole-identifier
   occurrence in the sum span, case-insensitively, outside quoted
   element literals. */
static int ident_char(char c) {
  return isalnum((unsigned char)c)||c=='_'||c=='@'||c==MAPMARK||c=='?';
}

static void dedup_ident_rename(char *span, const char *a, const char *b, size_t cap) {
  size_t alen=strlen(a),blen=strlen(b);
  char *p=span;
  int in_quote=0;
  if (alen==0) return;
  while (*p!='\0') {
    if (*p=='\"') {
      in_quote=!in_quote;
      p++;
      continue;
    }
    if (!in_quote&&strncasecmp(p,a,alen)==0&&(p==span||!ident_char(p[-1]))&&!ident_char(p[alen])) {
      size_t len=strlen(span);
      if (len-alen+blen>=cap) {
        errmsg("Error: renaming sum index %s does not fit the statement buffer (%d characters)\n",a,(int)TABREADLINE);
        MPI_Abort(PETSC_COMM_WORLD,1);
      }
      memmove(p+blen,p+alen,strlen(p+alen)+1);
      memcpy(p,b,blen);
      p+=blen;
      continue;
    }
    p++;
  }
}

int sum_dedup_indices(char *formulain) {
  char line[TABREADLINE],finditem[TABREADLINE],replitem[TABREADLINE],newreplitem[TABREADLINE],temp[TABREADLINE],replitem1[TABREADLINE],newreplitem1[TABREADLINE];
  char*readitem,*line1;
  char syntax[]="sum(";
  int nsum,i,j,l,k,k1,k2;
  nsum=sum_count(formulain,syntax);
  strcpy(line,formulain);
  readitem=line;
  l=0;
  for(i=0;i<nsum-1;i++){
    strcpy(line,formulain);
    k=str_find_token_ci(line,readitem,syntax);
    readitem=readitem+k;
    k1=str_find_ci(readitem,",");
    for(j=0;j<k1+1;j++)finditem[j]=readitem[j];
    finditem[j]='\0';
    readitem=readitem+k1;
    k2=str_find_token_ci(line,readitem,finditem);
    while(k2!=-1) {
      line1=readitem+k2;
      sprintf(temp, "%d", l);
      strcat(temp,"?");
      l++;
      /* the later sum was found case-insensitively (str_find_token_ci),
         but the rename below is a literal replace: take the index as it
         is spelled at the hit, not as the first sum spelled it, or a
         sum(S,...) found from sum(s,...) is never renamed and the
         re-find below loops forever (fuzz batch 13 hang class) */
      for(j=4;j<k1;j++){
        replitem[j-4]=line1[j];
        newreplitem[j-4]=line1[j];
      }
      replitem[j-4]='\0';
      newreplitem[j-4]='\0';
      strcat(newreplitem,temp);
      sum_extract(line1);
      strcpy(temp,line1);
      replitem1[0]='(';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,")");
      newreplitem1[0]='(';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,")");
      dedup_subst(line1,replitem1,newreplitem1,(size_t)TABREADLINE-(size_t)(line1-line));

      replitem1[0]='(';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,"]");
      newreplitem1[0]='(';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,"]");
      dedup_subst(line1,replitem1,newreplitem1,(size_t)TABREADLINE-(size_t)(line1-line));

      replitem1[0]='(';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,"}");
      newreplitem1[0]='(';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,"}");
      dedup_subst(line1,replitem1,newreplitem1,(size_t)TABREADLINE-(size_t)(line1-line));

      replitem1[0]='(';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,",");
      newreplitem1[0]='(';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,",");
      dedup_subst(line1,replitem1,newreplitem1,(size_t)TABREADLINE-(size_t)(line1-line));

      replitem1[0]='(';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1," ");
      newreplitem1[0]='(';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1," ");
      dedup_subst(line1,replitem1,newreplitem1,(size_t)TABREADLINE-(size_t)(line1-line));


      replitem1[0]='[';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,")");
      newreplitem1[0]='[';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,")");
      dedup_subst(line1,replitem1,newreplitem1,(size_t)TABREADLINE-(size_t)(line1-line));

      replitem1[0]='[';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,"]");
      newreplitem1[0]='[';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,"]");
      dedup_subst(line1,replitem1,newreplitem1,(size_t)TABREADLINE-(size_t)(line1-line));

      replitem1[0]='[';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,"}");
      newreplitem1[0]='[';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,"}");
      dedup_subst(line1,replitem1,newreplitem1,(size_t)TABREADLINE-(size_t)(line1-line));

      replitem1[0]='[';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,",");
      newreplitem1[0]='[';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,",");
      dedup_subst(line1,replitem1,newreplitem1,(size_t)TABREADLINE-(size_t)(line1-line));

      replitem1[0]='[';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1," ");
      newreplitem1[0]='[';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1," ");
      dedup_subst(line1,replitem1,newreplitem1,(size_t)TABREADLINE-(size_t)(line1-line));

      replitem1[0]='{';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,")");
      newreplitem1[0]='{';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,")");
      dedup_subst(line1,replitem1,newreplitem1,(size_t)TABREADLINE-(size_t)(line1-line));

      replitem1[0]='{';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,"]");
      newreplitem1[0]='{';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,"]");
      dedup_subst(line1,replitem1,newreplitem1,(size_t)TABREADLINE-(size_t)(line1-line));

      replitem1[0]='{';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,"}");
      newreplitem1[0]='{';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,"}");
      dedup_subst(line1,replitem1,newreplitem1,(size_t)TABREADLINE-(size_t)(line1-line));

      replitem1[0]='{';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,",");
      newreplitem1[0]='{';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,",");
      dedup_subst(line1,replitem1,newreplitem1,(size_t)TABREADLINE-(size_t)(line1-line));

      replitem1[0]='{';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1," ");
      newreplitem1[0]='{';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1," ");
      dedup_subst(line1,replitem1,newreplitem1,(size_t)TABREADLINE-(size_t)(line1-line));

      replitem1[0]=' ';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,")");
      newreplitem1[0]=' ';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,")");
      dedup_subst(line1,replitem1,newreplitem1,(size_t)TABREADLINE-(size_t)(line1-line));

      replitem1[0]=' ';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,"]");
      newreplitem1[0]=' ';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,"]");
      dedup_subst(line1,replitem1,newreplitem1,(size_t)TABREADLINE-(size_t)(line1-line));

      replitem1[0]=' ';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,"}");
      newreplitem1[0]=' ';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,"}");
      dedup_subst(line1,replitem1,newreplitem1,(size_t)TABREADLINE-(size_t)(line1-line));

      replitem1[0]=' ';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,",");
      newreplitem1[0]=' ';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,",");
      dedup_subst(line1,replitem1,newreplitem1,(size_t)TABREADLINE-(size_t)(line1-line));

      replitem1[0]=' ';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1," ");
      newreplitem1[0]=' ';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1," ");
      dedup_subst(line1,replitem1,newreplitem1,(size_t)TABREADLINE-(size_t)(line1-line));

      replitem1[0]=',';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,")");
      newreplitem1[0]=',';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,")");
      dedup_subst(line1,replitem1,newreplitem1,(size_t)TABREADLINE-(size_t)(line1-line));

      replitem1[0]=',';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,"]");
      newreplitem1[0]=',';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,"]");
      dedup_subst(line1,replitem1,newreplitem1,(size_t)TABREADLINE-(size_t)(line1-line));

      replitem1[0]=',';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,"}");
      newreplitem1[0]=',';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,"}");
      dedup_subst(line1,replitem1,newreplitem1,(size_t)TABREADLINE-(size_t)(line1-line));

      replitem1[0]=',';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,",");
      newreplitem1[0]=',';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,",");
      dedup_subst(line1,replitem1,newreplitem1,(size_t)TABREADLINE-(size_t)(line1-line));

      replitem1[0]=',';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1," ");
      newreplitem1[0]=',';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1," ");
      dedup_subst(line1,replitem1,newreplitem1,(size_t)TABREADLINE-(size_t)(line1-line));

      dedup_ident_rename(line1,replitem,newreplitem,(size_t)TABREADLINE-(size_t)(line1-line));

      if (strcmp(temp,line1)==0) {
        errmsg("Error: cannot rename the repeated sum index %s in: %s\n",replitem,temp);
        MPI_Abort(PETSC_COMM_WORLD,1);
      }
      dedup_subst_back(formulain,temp,line1,(size_t)TABREADLINE);
      strcpy(line,formulain);
      k2=str_find_token_ci(line,readitem,finditem);
    }
  }
  return 0;
}


static int tab_preprocess_run(char *filename, char *newtabfile);

/* the declaration index tab_read_set_name builds lives for one
   preprocess run, whichever way the run ends */
/* is the identifier starting at tok the name of a declared coefficient
   (declaration index of the preprocessed TAB)? */
static void ps_decl_name(char *line, int kwlen, char *out);
static int decl_index_build(char *filename);
static char **decl_cof;
static int n_decl_cof;
static int decl_coef_declared(char *tabfile, const char *tok) {
  char nm[NAMESIZE],dn[NAMESIZE];
  int k=0,i;
  while (cond_is_namec(tok[k])&&k<NAMESIZE-1) { nm[k]=(char)tolower((int)tok[k]); k++; }
  nm[k]='\0';
  if (k==0||decl_index_build(tabfile)<0) return 0;
  for (i=0; i<n_decl_cof; i++) {
    ps_decl_name(decl_cof[i],11,dn);
    if (str_cmp_ci(dn,nm)==0) return 1;
  }
  return 0;
}

/* p_/c_ names (manual 9.2.2): a levels variable X brings one linear
   variable, p_X (percent change) or c_X (change), or the LINEAR_NAME=
   it is given; a linear variable is referred to by its own name. The
   variables' kinds are read once from the declarations, with Variable
   (default=...) statements applied in order. */
typedef struct { char name[NAMESIZE]; int kind; char lin[NAMESIZE]; } pc_var; /* kind 1 linear, 2 levels percent, 3 levels change, 4 levels with LINEAR_NAME/VAR (its linear variable in lin) */
static pc_var *pc_vars=NULL;
static int n_pc_vars=0;

static int pc_kinds_build(const char *tabfile) {
  FILE *f=teems_fopen((char *)tabfile,"r");
  char line[TABREADLINE+1],defval[NAMESIZE];
  int lev=0,chg=0,cap=0;
  n_pc_vars=0;
  if (f==NULL) return -1;
  while (fgets(line,TABREADLINE,f)) {
    char *p,name[NAMESIZE],linname[NAMESIZE]="",linvarnm[NAMESIZE]="";
    int l=0,c=chg,linvar=0,k;
    if (strncmp(line,"variable",8)!=0) continue;
    if (tab_default_value(line,defval)) {
      if (strcmp(defval,"levels")==0) lev=1;
      else if (strcmp(defval,"linear")==0) lev=0;
      else if (strcmp(defval,"change")==0) chg=1;
      else if (strcmp(defval,"percent_change")==0) chg=0;
      continue;
    }
    l=lev;
    for (p=line+8; *p==' '; p++) {}
    while (*p=='('&&strncmp(p,"(all,",5)!=0&&strncmp(p,"(all ",5)!=0) {
      char tok[NAMESIZE];
      int t=0;
      for (p++; *p!='\0'&&*p!=')'; p++) {
        if (*p==','||*p==' ') {
          tok[t]='\0';
          if (strcmp(tok,"levels")==0) l=1;
          else if (strcmp(tok,"linear")==0) l=0;
          else if (strcmp(tok,"change")==0) c=1;
          else if (strcmp(tok,"percent_change")==0) c=0;
          else if (strncmp(tok,"linear_name=",12)==0) strcpy(linname,tok+12);
          else if (strncmp(tok,"linear_var=",11)==0) { linvar=1; strcpy(linvarnm,tok+11); }
          t=0;
          continue;
        }
        if (t<NAMESIZE-1) tok[t++]=*p;
      }
      tok[t]='\0';
      if (strcmp(tok,"levels")==0) l=1;
      else if (strcmp(tok,"linear")==0) l=0;
      else if (strcmp(tok,"change")==0) c=1;
      else if (strcmp(tok,"percent_change")==0) c=0;
      else if (strncmp(tok,"linear_name=",12)==0) strcpy(linname,tok+12);
      else if (strncmp(tok,"linear_var=",11)==0) { linvar=1; strcpy(linvarnm,tok+11); }
      if (*p==')') p++;
      while (*p==' ') p++;
    }
    while (strncmp(p,"(all",4)==0) {
      int d=0;
      for (; *p!='\0'; p++) { if (*p=='(') d++; else if (*p==')') { d--; if (d==0) { p++; break; } } }
      while (*p==' ') p++;
    }
    for (k=0; cond_is_namec(*p)&&k<NAMESIZE-1; p++) name[k++]=*p;
    name[k]='\0';
    if (k==0) continue;
    if (n_pc_vars+2>cap) {
      pc_var *g;
      cap=cap?2*cap:128;
      g=(pc_var *)realloc(pc_vars,cap*sizeof(pc_var));
      if (g==NULL) { fclose(f); return -1; }
      pc_vars=g;
    }
    /* LINEAR_NAME/LINEAR_VAR: the declared name is a coefficient, the
       linear variable is the named one */
    if (l&&(linname[0]!='\0'||linvar)) {
      if (linname[0]!='\0') {
        strcpy(pc_vars[n_pc_vars].name,linname);
        pc_vars[n_pc_vars++].kind=1;
      }
      strcpy(pc_vars[n_pc_vars].name,name);
      strcpy(pc_vars[n_pc_vars].lin,linname[0]!='\0'?linname:linvarnm);
      pc_vars[n_pc_vars++].kind=4;
      continue;
    }
    strcpy(pc_vars[n_pc_vars].name,name);
    pc_vars[n_pc_vars++].kind=l?(c?3:2):1;
  }
  fclose(f);
  return 0;
}

static int pc_find(const char *nm, size_t len) {
  int i;
  for (i=0; i<n_pc_vars; i++) if (strlen(pc_vars[i].name)==len&&strncmp(pc_vars[i].name,nm,len)==0) return i;
  return -1;
}

static int pc_kind(const char *nm, size_t len) {
  int i=pc_find(nm,len);
  return (i<0||pc_vars[i].kind==4)?0:pc_vars[i].kind;
}

/* check (and, in equations and updates, rewrite c_X to the p_X column
   token) every p_/c_ reference of one statement line; -1 on a named
   error */
static int pc_tokens_line(char *line, const char *tabfile, int rewrite) {
  char *n;
  int q=0;
  for (n=line; *n!='\0'; n++) {
    char tok[NAMESIZE];
    int k=0,kind;
    if (*n=='"') { q=!q; continue; }
    if (q||(n[0]!='p'&&n[0]!='c')||n[1]!='_'||(n>line&&cond_is_namec(n[-1]))) continue;
    while (cond_is_namec(n[k])&&k<NAMESIZE-1) { tok[k]=n[k]; k++; }
    tok[k]='\0';
    if (pc_kind(tok,strlen(tok))>0||decl_coef_declared((char *)tabfile,n)) { n+=k-1; continue; }
    kind=pc_kind(tok+2,strlen(tok+2));
    {
      int w=pc_find(tok+2,strlen(tok+2));
      if (w>=0&&pc_vars[w].kind==4) {
        errmsg("Error: %s names no variable: the linear variable of levels variable %s is %s (LINEAR_NAME/LINEAR_VAR, manual 9.2.2): %.200s\n",tok,tok+2,pc_vars[w].lin,line);
        return -1;
      }
    }
    if (kind==1) {
      errmsg("Error: %s refers to linear variable %s; a linear variable is named by itself (p_/c_ names belong to levels variables, manual 9.2.2): %.200s\n",tok,tok+2,line);
      return -1;
    }
    if (kind==2&&tok[0]=='c') {
      errmsg("Error: %s: levels variable %s is a percentage-change variable; its linear variable is p_%s (manual 9.2.2): %.200s\n",tok,tok+2,tok+2,line);
      return -1;
    }
    if (kind==3&&tok[0]=='p') {
      errmsg("Error: %s: levels variable %s is a change variable; its linear variable is c_%s (manual 9.2.2): %.200s\n",tok,tok+2,tok+2,line);
      return -1;
    }
    if (kind==3&&rewrite) n[0]='p';
    n+=k-1;
  }
  return 0;
}

/* PROD(i,S,e), MAXS(...) and MINS(...) (manual 11.4.4) become
   sum(<mark>i,S,e) in place -- "prod(" and "sum(<mark>" have the same
   length -- and fold as products, maxima or minima in the SUM
   machinery (sum_parse moves the mark into sum_def.fold) */
static void reduce_marks(char *line) {
  static const char *kw[3]={"prod","maxs","mins"};
  static const char mk[3]={SUM_MARK_PROD,SUM_MARK_MAXS,SUM_MARK_MINS};
  char *c;
  int q=0,k;
  for (c=line; *c!='\0'; c++) {
    if (*c=='"') { q=!q; continue; }
    if (q||(c>line&&cond_is_namec(c[-1]))) continue;
    for (k=0; k<3; k++) if (strncmp(c,kw[k],4)==0) {
        char *b=c+4;
        while (*b==' ') b++;
        if (*b!='('&&*b!='{'&&*b!='[') continue;
        memcpy(c,"sum",3);
        memmove(c+3,c+4,(size_t)(b-(c+4))+1);
        b[0]=mk[k];
        break;
      }
  }
}

/* a PROD/MAXS/MINS body takes no linear variable (manual 11.4.4):
   a declared linear variable, or the p_/c_ linear name of a levels
   variable, inside a marked reduction stops here, before the equation
   scans read the reduction as a sum */
static int reduce_linvar_check(const char *line) {
  const char *m;
  for (m=line; *m!='\0'; m++) {
    const char *e,*t;
    int d=0;
    if (sum_mark_fold(*m)==SUM_FOLD_SUM) continue;
    for (e=m-1; *e!='\0'; e++) {
      if (*e=='('||*e=='{'||*e=='[') d++;
      else if (*e==')'||*e=='}'||*e==']') { d--; if (d==0) break; }
    }
    for (t=m+1; t<e; ) {
      const char *s0=t;
      size_t len;
      int kd,pre;
      if (!cond_is_namec(*t)||(t>line&&cond_is_namec(t[-1]))) { t++; continue; }
      while (t<e&&cond_is_namec(*t)) t++;
      len=(size_t)(t-s0);
      kd=pc_kind(s0,len);
      pre=(len>2&&(s0[0]=='p'||s0[0]=='c')&&s0[1]=='_'&&pc_kind(s0+2,len-2)>=2);
      if (kd==1||pre) {
        errmsg("Error: linear variable %.*s inside PROD, MAXS or MINS; linear variables are not permitted there (manual 11.4.4)\n",(int)len,s0);
        return -1;
      }
    }
  }
  return 0;
}

/* A coefficient named p_NAME (legal beside a linear variable NAME,
   manual 9.2.2) reads like the p_ column token every linear variable
   reference becomes, so formulas bound it to the variable and equations
   took it for a column. Its tokens become p@NAME throughout -- '@'
   cannot begin a declared name -- and the coefficient dump reports the
   declared name. A linear variable p_NAME beside a change levels
   variable NAME, whose c_NAME references become p_NAME column tokens,
   is renamed the same way; the structure files, closures and shocks
   use the declared name. Runs before the c_ -> p_ rewrite. */
static int follows_file_kw(const char *line, const char *c) {
  const char *b=c;
  while (b>line&&b[-1]==' ') b--;
  while (b>line&&b[-1]==')') {
    int d=0;
    do { b--; if (*b==')') d++; else if (*b=='(') d--; } while (b>line&&d>0);
    while (b>line&&b[-1]==' ') b--;
  }
  return b-line>=4&&strncmp(b-4,"file",4)==0&&(b-4==line||!cond_is_namec(b[-5]));
}

static int pcoef_rename(const char *tabfile) {
  FILE *f=teems_fopen((char *)tabfile,"r"),*fo;
  char line[TABREADLINE+1],tmp[TABREADLINE],out[TABREADLINE+64];
  char (*nm)[NAMESIZE]=NULL,(*both)[NAMESIZE]=NULL,*bothvar=NULL;
  int n=0,cap=0,i,nb=0,capb=0;
  if (f==NULL) return -1;
  while (fgets(line,TABREADLINE,f)) {
    char *p,name[NAMESIZE];
    int k=0,isvar=strncmp(line,"variable",8)==0;
    if (strncmp(line,"coefficient",11)!=0&&!isvar) continue;
    for (p=line+(isvar?8:11); *p==' '; p++) {}
    while (*p=='(') {
      int d=0;
      for (; *p!='\0'; p++) { if (*p=='(') d++; else if (*p==')') { d--; if (d==0) { p++; break; } } }
      while (*p==' ') p++;
    }
    while (cond_is_namec(*p)&&k<NAMESIZE-1) name[k++]=*p++;
    name[k]='\0';
    if (k<3||name[0]!='p'||name[1]!='_') continue;
    /* every p_-leading declaration, to leave a coefficient and a
       variable of one name to the name check (11.2.1) */
    if (nb==capb) {
      char (*g)[NAMESIZE],*gv;
      capb=capb?2*capb:16;
      g=realloc(both,capb*sizeof(*both));
      if (g!=NULL) both=g;
      gv=realloc(bothvar,capb);
      if (gv!=NULL) bothvar=gv;
      if (g==NULL||gv==NULL) { free(nm); free(both); free(bothvar); fclose(f); return -1; }
    }
    strcpy(both[nb],name);
    bothvar[nb++]=(char)isvar;
    if (isvar) {
      if (n_pc_vars==0&&pc_vars==NULL&&pc_kinds_build(tabfile)<0) { free(nm); free(both); free(bothvar); fclose(f); return -1; }
      if (pc_kind(name+2,strlen(name+2))!=3) continue;
    }
    if (n==cap) {
      char (*g)[NAMESIZE];
      cap=cap?2*cap:16;
      g=realloc(nm,cap*sizeof(*nm));
      if (g==NULL) { free(nm); free(both); free(bothvar); fclose(f); return -1; }
      nm=g;
    }
    strcpy(nm[n++],name);
  }
  fclose(f);
  {
    int a,b,w=0;
    for (a=0; a<n; a++) {
      int cf=0,vr=0;
      for (b=0; b<nb; b++) if (strcmp(both[b],nm[a])==0) { if (bothvar[b]) vr=1; else cf=1; }
      if (!(cf&&vr)) { if (w!=a) memcpy(nm[w],nm[a],NAMESIZE); w++; }
    }
    n=w;
  }
  free(both);
  free(bothvar);
  if (n==0) { free(nm); return 0; }
  snprintf(tmp,sizeof(tmp),"%s_pc",tabfile);
  f=teems_fopen((char *)tabfile,"r");
  fo=teems_fopen(tmp,"w");
  if (f==NULL||fo==NULL) { if (f) fclose(f); if (fo) fclose(fo); free(nm); return -1; }
  while (fgets(line,TABREADLINE,f)) {
    size_t o=0;
    int q=0,isfile=strncmp(line,"file",4)==0&&!cond_is_namec(line[4]);
    char *c;
    for (c=line; *c!='\0'; ) {
      if (*c=='"') q=!q;
      /* a logical file name is not a coefficient: File P_X and
         "Write P_X to file P_X" keep it, as the CMF binds it */
      if (!q&&!isfile&&c[0]=='p'&&c[1]=='_'&&(c==line||!cond_is_namec(c[-1]))&&!follows_file_kw(line,c)) {
        int k=0;
        while (cond_is_namec(c[k])) k++;
        for (i=0; i<n; i++) if ((int)strlen(nm[i])==k&&strncmp(nm[i],c,k)==0) break;
        if (i<n) {
          if (o+k+2>=sizeof(out)) { errmsg("Error: statement too long after renaming coefficient %s\n",nm[i]); fclose(f); fclose(fo); free(nm); return -1; }
          out[o++]='p';
          out[o++]='@';
          memcpy(out+o,c+2,k-2);
          o+=k-2;
          c+=k;
          continue;
        }
      }
      if (o+2>=sizeof(out)) { fclose(f); fclose(fo); free(nm); errmsg("Error: statement too long after renaming p_ coefficients\n"); return -1; }
      out[o++]=*c++;
    }
    out[o]='\0';
    fputs(out,fo);
  }
  fclose(f);
  fclose(fo);
  free(nm);
  if (rename(tmp,tabfile)!=0) { errmsg("Error: cannot rename %s\n",tmp); return -1; }
  return 0;
}

int tab_preprocess(char *filename, char *newtabfile) {
  int r;
  free(pc_vars);
  pc_vars=NULL;
  n_pc_vars=0;
  r=tab_preprocess_run(filename,newtabfile);
  tab_decl_index_free();
  free(pc_vars);
  pc_vars=NULL;
  n_pc_vars=0;
  return r;
}

/* a Formula whose target is a scalar coefficient: no quantifier list and
   no argument list on the left of '=' (qualifier groups aside) */
static int formula_lhs_scalar(const char *line) {
  const char *eq=strchr(line,'='),*p=line;
  if (eq==NULL) return 0;
  if (str_find_ci((char *)line,"formula")==0) p+=7;
  while (p<eq) {
    if (*p=='(') {
      if (strncmp(p,"(initial)",9)==0||strncmp(p,"(always)",8)==0) {
        p=strchr(p,')');
        if (p==NULL) return 0;
        p++;
        continue;
      }
      return 0;
    }
    p++;
  }
  return 1;
}

/* strong comments (manual 11.1.5): the markers nest -- one !]]! cancels
   one ![[! -- and everything between the outermost pair is comment. The
   depth rides across lines; the removed text becomes one blank. The two
   markers used to map to one character stripped in pairs by parity, so
   the text between an inner opener and its closer came back to life.
   Returns -1 on a closer with nothing open. */
static int strong_comment_strip(char *s, int *depth) {
  char *r=s,*w=s;
  while (*r!='\0') {
    if (strncmp(r,"![[!",4)==0) {
      if (*depth==0) *w++=' ';
      (*depth)++;
      r+=4;
      continue;
    }
    if (strncmp(r,"!]]!",4)==0) {
      if (*depth==0) return -1;
      (*depth)--;
      r+=4;
      continue;
    }
    if (*depth>0) {
      if (*r=='\n') *w++='\n';
      r++;
      continue;
    }
    *w++=*r++;
  }
  *w='\0';
  return 0;
}

/* first word of a statement, lowercased (letters and '_' only) */
static void stmt_first_word(const char *s, char *out, int cap) {
  int k=0;
  while (*s==' '||*s=='\t') s++;
  while (((*s>='a'&&*s<='z')||(*s>='A'&&*s<='Z')||*s=='_')&&k<cap-1) out[k++]=(char)tolower((int)*s++);
  out[k]='\0';
  if (cond_is_namec(*s)) out[0]='\0';
}

/* numeric constants in exponent notation (1e-5, 2.5e+3, 1.e-3; manual
   11.4.9 discourages them, R-generated text can carry them) are
   rewritten as plain decimals before any tokenizer runs: the expression
   splitters cut at the sign, so 1e-5 used to read as 1 - 5. The digits
   move, nothing is rounded, so strtod reads the same double either way.
   Quoted text is left alone; a token glued to a name (c1e5) is a name. */
static int num_exp_expand(char *s, size_t cap) {
  char out[TABREADLINE];
  size_t o=0,i=0,n=strlen(s);
  int inq=0;
  while (i<n) {
    char c=s[i];
    if (c=='"') inq=!inq;
    if (!inq&&(isdigit((unsigned char)c)||(c=='.'&&isdigit((unsigned char)s[i+1])))&&(i==0||(!cond_is_namec(s[i-1])&&s[i-1]!='.'))) {
      size_t j=i,ms,me,dp=0,nd=0,ee;
      char dig[TABREADLINE];
      int neg=0,k=0,any=0;
      while (isdigit((unsigned char)s[j])) { dig[nd++]=s[j++]; dp++; }
      if (s[j]=='.') { j++; while (isdigit((unsigned char)s[j])) dig[nd++]=s[j++]; }
      ms=i; me=j;
      ee=j;
      if (s[ee]=='e'||s[ee]=='E') {
        ee++;
        if (s[ee]=='+'||s[ee]=='-') { neg=(s[ee]=='-'); ee++; }
        while (isdigit((unsigned char)s[ee])) { if (k<1000) k=k*10+(s[ee]-'0'); ee++; any=1; }
      }
      if (any&&!cond_is_namec(s[ee])&&s[ee]!='.') {
        char num[TABREADLINE];
        long p,q,m=0;
        if (k>99||nd+(size_t)k+3>=sizeof(num)) {
          errmsg("Error: numeric constant %.*s is out of the supported range (exponent notation, manual 11.4.9)\n",(int)(ee-ms)>200?200:(int)(ee-ms),s+ms);
          return -1;
        }
        p=(long)dp+(neg?-k:k);
        if (p<=0) {
          num[m++]='0'; num[m++]='.';
          for (q=0; q<-p; q++) num[m++]='0';
          for (q=0; q<(long)nd; q++) num[m++]=dig[q];
        } else {
          for (q=0; q<p; q++) num[m++]=(q<(long)nd)?dig[q]:'0';
          if (p<(long)nd) {
            num[m++]='.';
            for (q=p; q<(long)nd; q++) num[m++]=dig[q];
          }
        }
        num[m]='\0';
        if (o+m>=sizeof(out)) return -1;
        memcpy(out+o,num,m);
        o+=m;
        i=ee;
        continue;
      }
      if (o+(me-ms)>=sizeof(out)) return -1;
      memcpy(out+o,s+ms,me-ms);
      o+=me-ms;
      i=me;
      continue;
    }
    if (o+1>=sizeof(out)) return -1;
    out[o++]=c;
    i++;
  }
  out[o]='\0';
  if (o>=cap) return -1;
  strcpy(s,out);
  return 0;
}

/* end of the operand starting at s[i] (a number, a name with an
   optional argument group, or a bracket group), -1 if none */
static long unary_operand_end(const char *s, long i) {
  long j=i;
  if (s[j]=='(') {
    int d=0;
    for (; s[j]!='\0'; j++) {
      if (s[j]=='(') d++;
      else if (s[j]==')') { d--; if (d==0) return j+1; }
    }
    return -1;
  }
  if (isdigit((unsigned char)s[j])||s[j]=='.') {
    while (isdigit((unsigned char)s[j])||s[j]=='.') j++;
    return j;
  }
  if (!cond_is_namec(s[j])&&s[j]!='$') return -1;
  if (s[j]=='$') j++;
  while (cond_is_namec(s[j])) j++;
  {
    long k=j;
    while (s[k]==' ') k++;
    if (s[k]=='(') {
      long e=unary_operand_end(s,k);
      if (e>0) return e;
    }
  }
  return j;
}

/* unary + and - bind before ^ (manual 11.4.1: -A+B^C/D is
   ((-A)+((B^C)/D))). levels.c parses that way but the formula engine
   raised first and negated after, so -X^2 was -(X^2) there and +(X)^2
   in the levels half of the same model. A signed operand next to a ^
   -- -X^2, 2*-3^2, X^-1 -- is bracketed here, before either engine
   reads the text, so both see (-X)^2 and X^(-1). A signed operand
   after * or / (A*-B, A/-B, 11.4.1) is bracketed the same way; the
   engines took the sign as an empty operand there. */
static int unary_pow_bracket(char *s, size_t cap) {
  long i,n;
  int inq=0;
  for (i=0; s[i]!='\0'; i++) {
    long b,e,k;
    char pc;
    if (s[i]=='"') { inq=!inq; continue; }
    if (inq||(s[i]!='-'&&s[i]!='+')) continue;
    for (b=i-1; b>=0&&s[b]==' '; b--) {}
    pc=(b>=0)?s[b]:'=';
    if (strchr("=(,*/^<>:+-",pc)==NULL) continue;
    for (k=i+1; s[k]==' '; k++) {}
    e=unary_operand_end(s,k);
    if (e<0) continue;
    {
      long f=e;
      while (s[f]==' ') f++;
      if (pc!='^'&&pc!='*'&&pc!='/'&&s[f]!='^') continue;
    }
    n=(long)strlen(s);
    if ((size_t)(n+3)>=cap) return -1;
    memmove(s+e+1,s+e,n-e+1);
    s[e]=')';
    memmove(s+i+1,s+i,e+1-i+(n+1-e));
    s[i]='(';
    i++;
  }
  return 0;
}

/* A statement with no keyword continues the previous statement's kind
   (manual 11.1.1), so a misspelt keyword was read as that kind:
   "Coeficient (all,r,REG) K(r);" after a Coefficient declared a
   coefficient named coeficient and exited 0. Where the continued kind
   cannot open with a name followed by a quantifier or by another word
   (declarations, formulas, updates; assertions for the quantifier), the
   first word is an unknown keyword. Returns 1 and the word in `word`. */
static int kwless_unknown(const char *s,const char *sticky,char *word,int cap) {
  int k=0,sp=0,decl;
  const char *p;
  decl=strcmp(sticky,"coefficient")==0||strcmp(sticky,"variable")==0||strcmp(sticky,"formula")==0||strcmp(sticky,"update")==0;
  if (!decl&&strcmp(sticky,"assertion")!=0) return 0;
  while (*s==' '||*s=='\t') s++;
  if (!isalpha((unsigned char)*s)) return 0;
  while ((isalnum((unsigned char)*s)||*s=='_')&&k<cap-1) word[k++]=(char)tolower((unsigned char)*s++);
  word[k]='\0';
  if (isalnum((unsigned char)*s)||*s=='_') return 0;
  while (*s==' '||*s=='\t') { s++; sp=1; }
  if (*s=='(') {
    for (p=s+1; *p==' '; p++);
    if (strncmp(p,"all,",4)==0||strncmp(p,"all ,",5)==0) return 1;
    return 0;
  }
  return decl&&sp&&isalpha((unsigned char)*s);
}

/* IF(i IN S, v) (manual 11.4.7) becomes SUM(j, S intersect T: j = i, v)
   with j renamed for i in v: T is the set i ranges over where the IF
   stands (rule 1), j ranges over a subset of S inside v as rule 2
   requires, and the one-element sum is v where i is in S and 0
   elsewhere, in any position of a formula, update, assertion or
   equation. The intersection set and its subset statements go to
   `pre`, written ahead of the statement. Returns -1 on a named error. */
static int ifin_namec(char c) {
  return isalnum((unsigned char)c)||c=='_'||c=='@';
}

static int ifin_word_at(const char *s, long k, const char *w) {
  size_t n=strlen(w);
  return strncmp(s+k,w,n)==0&&(k==0||!ifin_namec(s[k-1]))&&!ifin_namec(s[k+n]);
}

static int if_in_lower(char *s, size_t cap, char *pre, size_t precap, int *seq) {
  long k;
  int inq=0;
  for (k=0; s[k]!='\0'; k++) {
    long o,c,cm,b,d,j;
    char idx[NAMESIZE],set[NAMESIZE],rng[NAMESIZE],nidx[NAMESIZE],jidx[NAMESIZE+16],nset[NAMESIZE],cond[TABREADLINE];
    int n;
    if (s[k]=='"') { inq=!inq; continue; }
    if (inq||!ifin_word_at(s,k,"if")) continue;
    for (o=k+2; s[o]==' '; o++) {}
    if (s[o]!='(') continue;
    for (c=o,d=0,cm=-1; s[c]!='\0'; c++) {
      if (s[c]=='(') d++;
      else if (s[c]==')') { d--; if (d==0) break; }
      else if (s[c]==','&&d==1&&cm<0) cm=c;
    }
    if (s[c]=='\0'||cm<0||cm-o-1>=(long)sizeof(cond)) continue;
    memcpy(cond,s+o+1,cm-o-1);
    cond[cm-o-1]='\0';
    if (sscanf(cond," %255[a-z0-9_@] in %255[a-z0-9_@] %n",idx,set,&n)<2||cond[n]!='\0') {
      if (strstr(cond," in ")!=NULL&&strpbrk(cond,"&|`")!=NULL) {
        errmsg("Error: a condition \"index IN set\" cannot be combined with AND, OR or NOT (manual 11.4.7 rule 5): %.200s\n",cond);
        return -1;
      }
      continue;
    }
    /* the set the index ranges over: its nearest quantifier or sum
       before the IF */
    rng[0]='\0';
    for (b=k-1; b>=0&&rng[0]=='\0'; b--) {
      long p;
      if (strncmp(s+b,"(all",4)==0) p=b+4;
      else if (ifin_word_at(s,b,"sum")) {
        for (p=b+3; s[p]==' '; p++) {}
        if (s[p]!='(') continue;
      } else continue;
      p++;
      while (s[p]==' '||(s[p]==','&&s[b]=='(')) p++;
      if (sscanf(s+p,"%255[a-z0-9_@]",nidx)!=1||strcmp(nidx,idx)!=0) continue;
      p+=(long)strlen(nidx);
      while (s[p]==' ') p++;
      if (s[p]!=',') continue;
      for (p++; s[p]==' '; p++) {}
      if (sscanf(s+p,"%255[a-z0-9_@]",rng)!=1) rng[0]='\0';
    }
    if (rng[0]=='\0') {
      errmsg("Error: index %s in \"%s in %s\" is not active where the IF stands; it must come from an ALL quantifier or an enclosing SUM (manual 11.4.7 rule 1)\n",idx,idx,set);
      return -1;
    }
    (*seq)++;
    snprintf(nset,sizeof(nset),"if_in%d",*seq);
    snprintf(jidx,sizeof(jidx),"%s@in%d",idx,*seq);
    {
      char val[TABREADLINE],out[TABREADLINE],tail[TABREADLINE];
      size_t vo=0;
      long v;
      int vq=0;
      for (v=cm+1; v<c&&vo+NAMESIZE<sizeof(val); v++) {
        if (s[v]=='"') vq=!vq;
        if (!vq&&ifin_word_at(s,v,idx)) {
          vo+=snprintf(val+vo,sizeof(val)-vo,"%s",jidx);
          v+=(long)strlen(idx)-1;
          continue;
        }
        val[vo++]=s[v];
      }
      val[vo]='\0';
      strcpy(tail,s+c+1);
      s[k]='\0';
      if ((size_t)snprintf(out,sizeof(out),"%ssum(%s,%s: %s = %s,%s)%s",s,jidx,nset,jidx,idx,val,tail)>=cap) {
        errmsg("Error: statement too long after rewriting \"%s in %s\" (manual 11.4.7)\n",idx,set);
        return -1;
      }
      strcpy(s,out);
    }
    j=(long)strlen(pre);
    if ((size_t)snprintf(pre+j,precap-j,"set %s = %s intersect %s;\nsubset %s is subset of %s;\nsubset %s is subset of %s;\n",nset,set,rng,nset,set,nset,rng)>=precap-j) {
      errmsg("Error: too many \"index IN set\" conditions in one statement\n");
      return -1;
    }
    k=-1;
    inq=0;
  }
  return 0;
}

/* In an equation, a sum condition naming the index of an enclosing sum
   (sum(s,S, sum(j,T: j = s, w(j))) -- the IN rewrite's shape, and any
   index condition of that form, manual 11.4.11) cannot gate the columns
   of the inner sum's variables: the enclosing index is summed inside
   their coefficient, not looped with the column. The condition moves
   into the coefficient, sum(j,T, if(cond,1)*(body)), where both indices
   are in scope; other sums keep their text. */
static int eq_sumcond_outer(char *s, size_t cap) {
  long k;
  int changed=1,guard=0;
  while (changed&&guard++<256) {
    changed=0;
    for (k=0; s[k]!='\0'&&!changed; k++) {
      long o,c,d,colon=-1,c1=-1,c2=-1,b;
      char cond[TABREADLINE],outer[64][NAMESIZE];
      int nout=0,q=0,hit=0;
      if (!ifin_word_at(s,k,"sum")) continue;
      for (o=k+3; s[o]==' '; o++) {}
      if (s[o]!='(') continue;
      for (c=o,d=0; s[c]!='\0'; c++) {
        if (s[c]=='"') { q=!q; continue; }
        if (q) continue;
        if (s[c]=='(') d++;
        else if (s[c]==')') { d--; if (d==0) break; }
        else if (d==1&&s[c]==':'&&colon<0&&c1>=0&&c2<0) colon=c;
        else if (d==1&&s[c]==',') { if (c1<0) c1=c; else if (c2<0) c2=c; }
      }
      if (s[c]=='\0'||colon<0||c2<0||c2-colon-1>=(long)sizeof(cond)) continue;
      memcpy(cond,s+colon+1,c2-colon-1);
      cond[c2-colon-1]='\0';
      /* the indices of the sums enclosing this one */
      for (b=0; b<k&&nout<64; b++) {
        long bo,bc,bd=0;
        int bq=0;
        if (!ifin_word_at(s,b,"sum")) continue;
        for (bo=b+3; s[bo]==' '; bo++) {}
        if (s[bo]!='(') continue;
        for (bc=bo; s[bc]!='\0'; bc++) {
          if (s[bc]=='"') { bq=!bq; continue; }
          if (bq) continue;
          if (s[bc]=='(') bd++;
          else if (s[bc]==')') { bd--; if (bd==0) break; }
        }
        if (bc<c) continue;
        if (sscanf(s+bo+1," %255[a-z0-9_@]",outer[nout])==1) nout++;
      }
      {
        int i;
        for (i=0; i<nout&&!hit; i++) {
          long v;
          int vq=0;
          for (v=0; cond[v]!='\0'; v++) {
            if (cond[v]=='"') { vq=!vq; continue; }
            if (!vq&&ifin_word_at(cond,v,outer[i])) { hit=1; break; }
          }
        }
      }
      if (!hit) continue;
      {
        char out[TABREADLINE];
        long n;
        n=snprintf(out,sizeof(out),"%.*s, if(%s,1)*(%.*s))%s",(int)colon,s,cond,(int)(c-c2-1),s+c2+1,s+c+1);
        if (n<0||(size_t)n>=cap) {
          errmsg("Error: equation too long after moving a sum condition into its coefficient: %.120s\n",s);
          return -1;
        }
        strcpy(s,out);
        changed=1;
      }
    }
  }
  return 0;
}

#define MAXLITDEDUP 64
static int tab_preprocess_run(char *filename, char *newtabfile) {
  FILE * filehandle,*fout;
  char dedup_set[MAXLITDEDUP][NAMESIZE],dedup_ele[MAXLITDEDUP][NAMESIZE],dedup_idx[MAXLITDEDUP][NAMESIZE];
  int ndedup=0;
  char line[TABREADLINE]="\0",line1[TABREADLINE],line2[TABREADLINE],indx[NAMESIZE],indx1[NAMESIZE],indx2[NAMESIZE],*readitem,*readitem1,commsyntax[NAMESIZE],readline[TABREADLINE],readline1[TABREADLINE],*n,newtabfile1[TABREADLINE];
  char setname[NAMESIZE],newset[NAMESIZE],varname[NAMESIZE],*n1,setelement[TABREADLINE];//,*ne,*np;//,*n2;
  char assertmsg[TABREADLINE],*am1,*am2;
  char rawline[TABLINESIZE],*rawpos;
  char ifin_pre[TABREADLINE]="";
  int ifin_seq=0;
  filehandle = teems_fopen(filename,"r");
  if(filehandle==NULL){
    errmsg("Error: cannot open %s\n",filename);
    return -1;
  }
  int check,i1,i,setindx,l1,l2,l3,l4,k1,k2,l5;//,necheck,npcheck;//,j;,check1
  strcpy(newtabfile1,newtabfile);
  str_replace_all(newtabfile1,".","1.");
  fout = teems_fopen(newtabfile1,"w");
  readline[0]='\0';
  commsyntax[0]='\0';
  /* GEMPACK allows several statements on one physical line: each raw
     line is fed to the accumulator piecewise, cut after every ';' that
     lies outside a ! comment, a ![[! !]]! block comment and a # label #
     (parity taken over the pending statement + the piece so far) and
     is followed by further text, so every iteration below sees at most
     one statement terminator */
  rawpos=NULL;
  int sdepth=0;
  long rawno=0;
  while (rawpos!=NULL||fgets(rawline,TABLINESIZE,filehandle)) {
    if (rawpos==NULL) {
      rawno++;
      /* a DOS end-of-file mark (Ctrl-Z) after the last statement read
         as a keyword-less statement continuing the previous kind; tabs
         and form feeds are blanks (manual 11.1.2) */
      {
        char *z;
        for (z=rawline; *z!='\0'; z++) if (*z=='\032'||*z=='\t'||*z=='\f'||*z=='\v') *z=' ';
      }
      if (strong_comment_strip(rawline,&sdepth)<0) {
        errmsg("Error: '!]]!' at line %ld of the TAB file closes a strong comment that was never opened (manual 11.1.5)\n",rawno);
        fclose(filehandle);
        fclose(fout);
        return -1;
      }
      rawpos=rawline;
    }
    {
      int pb=str_count_char(readline,'!')%2,ph=str_count_char(readline,'#')%2;
      char *c,*cut=NULL;
      for (c=rawpos; *c!='\0'; c++) {
        if (*c=='!') pb^=1;
        else if (*c=='#') ph^=1;
        else if (*c==';'&&!pb&&!ph) {
          char *d=c+1;
          while (*d==' '||*d=='\t') d++;
          if (*d!='\0'&&*d!='\n'&&*d!='\r') cut=c+1;
          break;
        }
      }
      if (cut!=NULL) {
        size_t len=(size_t)(cut-rawpos);
        memcpy(line,rawpos,len);
        line[len]='\0';
        rawpos=cut;
      } else {
        strcpy(line,rawpos);
        rawpos=NULL;
      }
    }
    if (strlen(readline)+strlen(line)>=sizeof(readline)) {
      errmsg("Error: TAB statement too long (exceeds %d chars; a Backsolve/Substitute expansion in teems-R can grow a statement past this): %.120s ...\n",TABREADLINE,readline);
      return -1;
    }
    strcat(readline,line);
    n=strrchr(line,';');
    i1=str_count_char(readline,'!');
    while (i1>1) {
      while (str_strip_comment(readline,"!"));
      i1-=2;
    }
    if (n!=NULL&&i1==0) {
      check=0;
      /* capture the first # label # BEFORE the strip and the lowercase
         pass: assertion messages must survive preprocessing verbatim
         (manual 10.14; batch-1 residual (a)) -- reinserted at write-out
         for assertion statements only, dropped for everything else */
      assertmsg[0]='\0';
      am1=strchr(readline,'#');
      if (am1!=NULL) {
        am2=strchr(am1+1,'#');
        if (am2!=NULL&&am2-am1>1&&(size_t)(am2-am1)<sizeof(assertmsg)) {
          strncpy(assertmsg,am1+1,am2-am1-1);
          assertmsg[am2-am1-1]='\0';
          for (am1=assertmsg; *am1!='\0'; am1++)
            if (*am1=='\n'||*am1=='\r') *am1=' ';
        }
      }
      while (str_strip_comment(readline,"#"));
      while (str_replace_all(readline,"\n", " "));
      while (str_replace_all(readline,"\r", " "));
      while (str_replace_all(readline,"  ", " "));
      while (str_replace_all(readline,", ", ","));
      while (str_replace_all(readline," ,", ","));
      while (str_replace_all(readline,"( ", "("));
      while (str_replace_all(readline," )", ")"));
      /* the keyword/set spacing below inserts at most 14 characters
         ("set(" and "(" once each, ten "keyword(" forms) with the
         unbounded str_replace_* helpers: a statement within that of
         the buffer overflowed it (fuzz batch 13) */
      if (strlen(readline)+16>=sizeof(readline)) {
        errmsg("Error: TAB statement too long (exceeds %d chars; a Backsolve/Substitute expansion in teems-R can grow a statement past this): %.120s ...\n",TABREADLINE,readline);
        fclose(filehandle);
        fclose(fout);
        return -1;
      }
      k1=0;
      k2=0;
      while (readline[k1]!= '\0') {
        if(readline[k1]=='\"') {
          if(k2==0) k2=1;
          else k2=0;
        } else {
          if(k2==0)readline[k1]=tolower((int)readline[k1]);
        }
        k1++;
      }
      str_replace_first(readline,"set(", "set (");
      str_replace_first(readline,"set[", "set [");
      str_replace_first(readline,"set{", "set {");
      if(str_find_ci(readline,"set ")==1||str_find_ci(readline,"set ")==0) {
        strcpy(commsyntax,"set");
        str_replace_all(readline,"(", " (");
        while (str_replace_all(readline,"  ", " "));
        check=1;
      }
      if(check==0){
          while (str_replace_all(readline,"[", "("));
          while (str_replace_all(readline,"]", ")"));
          while (str_replace_all(readline,"{", "("));
          while (str_replace_all(readline,"}", ")"));
      }
      str_replace_first(readline,"equation(", "equation (");
      str_replace_first(readline,"formula(", "formula (");
      str_replace_first(readline,"coefficient(", "coefficient (");
      str_replace_first(readline,"variable(", "variable (");
      str_replace_first(readline,"update(", "update (");
      str_replace_first(readline,"read(", "read (");
      str_replace_first(readline,"write(", "write (");
      str_replace_first(readline,"zerodivide(", "zerodivide (");
      str_replace_first(readline,"mapping(", "mapping (");
      str_replace_first(readline,"complementarity(", "complementarity (");
      str_replace_first(readline,"subset(", "subset (");
      if(str_find_ci(readline,"subset ")==1||str_find_ci(readline,"subset ")==0) {
        strcpy(commsyntax,"subset");
        check=1;
        /* Subset (BY_ELEMENTS) is the default reading (manual 10.2);
           (BY_NUMBERS) is obsolete */
        if(str_find_ci(readline,"(by_numbers)")>-1) {
          errmsg("Error: Subset (by_numbers) is obsolete and not supported; list the subset's elements by name (manual 10.2)\n");
          fclose(filehandle);
          fclose(fout);
          return -1;
        }
        str_replace_first(readline,"(by_elements)","");
      }
      if(str_find_ci(readline,"file ")==1||str_find_ci(readline,"file ")==0) {
        strcpy(commsyntax,"file");
        check=1;
      }
      /* MAPPING declarations (manual 11.9); without recognition the
         sticky-keyword continuation would hand them to the previous
         statement's scanner */
      if(str_find_ci(readline,"mapping ")==1||str_find_ci(readline,"mapping ")==0) {
        strcpy(commsyntax,"mapping");
        check=1;
      }
      /* COMPLEMENTARITY statements (manual 10.17; design doc section
         7); the embedded "variable =" qualifier cannot false-match the
         variable check below because those tests anchor at position
         0/1 */
      if(str_find_ci(readline,"complementarity ")==1||str_find_ci(readline,"complementarity ")==0) {
        strcpy(commsyntax,"complementarity");
        check=1;
      }
      if (str_find_ci(readline,"coefficient ")==1||str_find_ci(readline,"coefficient ")==0) {//if (ha_cgefind(readline,"coefficient ")>-1) {
        strcpy(commsyntax,"coefficient");
        check=1;
      }
      if (str_find_ci(readline,"variable ")==1||str_find_ci(readline,"variable ")==0) {//if (ha_cgefind(readline,"variable ")>-1) {
        strcpy(commsyntax,"variable");
        check=1;
      }
      if((str_find_ci(readline,"read ")==1||str_find_ci(readline,"read ")==0)&&check==0&&str_find_ci(readline,"read elements")==-1) {
        strcpy(commsyntax,"read");
        check=1;
      }
      if (str_find_ci(readline,"formula ")==1||str_find_ci(readline,"formula ")==0
          ||str_find_ci(readline,"formula&")==1||str_find_ci(readline,"formula&")==0) {
        /* the no-space Formula&Equation spelling (10.9.1) must not
           fall through to the sticky-keyword prepend */
        strcpy(commsyntax,"formula");
        check=1;
      }
      if (str_find_ci(readline,"equation ")==1||str_find_ci(readline,"equation ")==0) {//if (ha_cgefind(readline,"equation ")>-1) {
        strcpy(commsyntax,"equation");
        check=1;
      }

      if(str_find_ci(readline,"update ")==1||str_find_ci(readline,"update ")==0) {//if(ha_cgefind(readline,"update ")>-1) {
        strcpy(commsyntax,"update");
        check=1;
      }
      if(str_find_ci(readline,"zerodivide ")==1||str_find_ci(readline,"zerodivide ")==0) {//if(ha_cgefind(readline,"zerodivide ")>-1) {
        strcpy(commsyntax,"zerodivide");
        check=1;
      }

      if(str_find_ci(readline,"write ")==1||str_find_ci(readline,"write ")==0) {//if(ha_cgefind(readline,"write ")>-1) {
        strcpy(commsyntax,"write");
        check=1;
      }

      if(str_find_ci(readline,"assertion ")==1||str_find_ci(readline,"assertion ")==0) {//if(ha_cgefind(readline,"assertion ")>-1) {
        strcpy(commsyntax,"assertion");
        check=1;
      }
      /* condensation statements (GEMPACK manual 10.16): backsolve is
         honored downstream; omit/substitute must be resolved by teems-R
         before deployment and abort in backsolve_read */
      if(str_find_ci(readline,"backsolve ")==1||str_find_ci(readline,"backsolve ")==0) {
        strcpy(commsyntax,"backsolve");
        check=1;
      }
      if(str_find_ci(readline,"omit ")==1||str_find_ci(readline,"omit ")==0) {
        strcpy(commsyntax,"omit");
        check=1;
      }
      if(str_find_ci(readline,"substitute ")==1||str_find_ci(readline,"substitute ")==0) {
        strcpy(commsyntax,"substitute");
        check=1;
      }
      /* POSTSIM section markers (manual 10.18): recognized so the
         sticky-keyword mechanism cannot mangle them; the post-
         preprocess splitter routes the section contents */
      if(str_find_ci(readline,"postsim ")==1||str_find_ci(readline,"postsim ")==0||str_find_ci(readline,"postsim(")==1||str_find_ci(readline,"postsim(")==0) {
        strcpy(commsyntax,"postsim");
        check=1;
      }
      /* keywords the solver does not carry (manual 11.1.1 lists them):
         a keyword-less statement inherits the previous keyword, so
         "Loop (begin) ..." used to become "formula loop ..." and
         "Display X" a Write or a bogus declaration. TAB-file loops are
         fatal; DISPLAY and TRANSFER write files TEEMS does not produce
         (every coefficient is in the coefficient dump) and are dropped
         with a warning, their keyword-less continuations with them */
      if (check==0) {
        char kw[16];
        stmt_first_word(readline,kw,sizeof(kw));
        if (strcmp(kw,"loop")==0||strcmp(kw,"break")==0||strcmp(kw,"cycle")==0) {
          char up[16];
          for (k1=0; kw[k1]!='\0'; k1++) up[k1]=(char)toupper((int)kw[k1]);
          up[k1]='\0';
          errmsg("Error: %s statements are not supported (loops in TAB files, manual 11.18): %.120s\n",up,readline);
          fclose(filehandle);
          fclose(fout);
          return -1;
        }
        if (strcmp(kw,"display")==0||strcmp(kw,"transfer")==0) {
          strcpy(commsyntax,kw);
          check=1;
        }
        else if (commsyntax[0]!='\0') {
          char uw[NAMESIZE];
          if (kwless_unknown(readline,commsyntax,uw,sizeof(uw))) {
            char *t=readline;
            while (*t==' ') t++;
            errmsg("Error: unknown statement keyword '%s' (a statement without a keyword continues the previous %s statement, manual 11.1.1, and this one cannot): %.120s\n",uw,commsyntax,t);
            fclose(filehandle);
            fclose(fout);
            return -1;
          }
        }
      }
      if (strchr(readline,';')!=NULL&&(strcmp(commsyntax,"display")==0||strcmp(commsyntax,"transfer")==0)) {
        char *t=readline;
        while (*t==' ') t++;
        if (commsyntax[0]=='d') printf("Warning: DISPLAY statement ignored (TEEMS writes every coefficient to the coefficient dump; manual 10.12): %.120s\n",t);
        else printf("Warning: TRANSFER statement ignored (TEEMS writes no Header Array output files; manual 10.15): %.120s\n",t);
        readline[0]='\0';
        continue;
      }
      /* word operators (manual 11.4.5): comparisons to symbols and
         AND/OR/NOT to single characters while the blanks still
         delimit them; every reader strips blanks */
      if (strcmp(commsyntax,"formula")==0||strcmp(commsyntax,"equation")==0||strcmp(commsyntax,"update")==0||strcmp(commsyntax,"assertion")==0
          ||(strcmp(commsyntax,"set")==0&&strchr(readline,':')!=NULL)) {
        tab_wordops_normalize(readline);
        if (tab_logicops_normalize(readline,sizeof(readline))<0) {
          errmsg("Error: TAB statement too long after normalizing AND/OR/NOT: %.120s ...\n",readline);
          fclose(filehandle);
          fclose(fout);
          return -1;
        }
        if (strcmp(commsyntax,"set")!=0&&strchr(readline,';')!=NULL&&strstr(readline," in ")!=NULL
            &&if_in_lower(readline,sizeof(readline),ifin_pre,sizeof(ifin_pre),&ifin_seq)<0) {
          fclose(filehandle);
          fclose(fout);
          return -1;
        }
        if (strcmp(commsyntax,"equation")==0&&strchr(readline,';')!=NULL&&strchr(readline,':')!=NULL
            &&eq_sumcond_outer(readline,sizeof(readline))<0) {
          fclose(filehandle);
          fclose(fout);
          return -1;
        }
      }
      /* exponent constants and signed operands of ^ (manual 11.4.9,
         11.4.1) are rewritten once here, for every later reader */
      if (strcmp(commsyntax,"formula")==0||strcmp(commsyntax,"equation")==0||strcmp(commsyntax,"update")==0||strcmp(commsyntax,"assertion")==0
          ||strcmp(commsyntax,"coefficient")==0||strcmp(commsyntax,"variable")==0||strcmp(commsyntax,"zerodivide")==0||strcmp(commsyntax,"complementarity")==0) {
        if (strpbrk(readline,"eE")!=NULL&&num_exp_expand(readline,sizeof(readline))<0) {
          errmsg("Error: TAB statement too long after expanding exponent-notation constants: %.120s ...\n",readline);
          fclose(filehandle);
          fclose(fout);
          return -1;
        }
        if (strcmp(commsyntax,"coefficient")!=0&&strcmp(commsyntax,"variable")!=0&&strcmp(commsyntax,"zerodivide")!=0&&strcmp(commsyntax,"complementarity")!=0
            &&strpbrk(readline,"^*/")!=NULL&&unary_pow_bracket(readline,sizeof(readline))<0) {
          errmsg("Error: TAB statement too long after bracketing signed operands: %.120s ...\n",readline);
          fclose(filehandle);
          fclose(fout);
          return -1;
        }
      }
      /* reinsert the captured message ahead of the terminator; only
         assertions keep their label (sticky commsyntax covers the
         keyword-less continuation form) */
      if (assertmsg[0]!='\0'&&strcmp(commsyntax,"assertion")==0&&
          strlen(readline)+strlen(assertmsg)+8<sizeof(readline)&&
          (am2=strrchr(readline,';'))!=NULL) {
        sprintf(am2," #%s# ;",assertmsg);
      }
      if (strchr(readline,';')!=NULL) {
        /* every later pass reads this file a line at a time into
           TABREADLINE buffers: a statement whose keyword prefix pushes
           it past that was split mid-statement and the pieces
           re-joined past the reader's buffer (fuzz batch 13) */
        if (strlen(commsyntax)+strlen(readline)+3>TABREADLINE) {
          errmsg("Error: TAB statement too long (exceeds %d chars; a Backsolve/Substitute expansion in teems-R can grow a statement past this): %.120s ...\n",TABREADLINE,readline);
          fclose(filehandle);
          fclose(fout);
          return -1;
        }
        if (ifin_pre[0]!='\0') {
          fprintf(fout,"%s",ifin_pre);
          ifin_pre[0]='\0';
        }
        if (check==1) {
          if (readline[0]==' ') fprintf(fout,"%s\n",readline+1);
          else fprintf(fout,"%s\n",readline);
          readline[0]='\0';
        } else {
          if (readline[0]==' ') fprintf(fout,"%s%s\n",commsyntax,readline);
          else fprintf(fout,"%s %s\n",commsyntax,readline);
          readline[0]='\0';
        }
      }
    }
  }
  if (sdepth>0) {
    errmsg("Error: a strong comment opened with '![[!' in the TAB file is never closed by '!]]!' (%d still open at the end of the file; manual 11.1.5)\n",sdepth);
    fclose(filehandle);
    fclose(fout);
    return -1;
  }
  fclose(filehandle);
  fclose(fout);
  if (pcoef_rename(newtabfile1)<0) return -1;
  free(pc_vars);
  pc_vars=NULL;
  n_pc_vars=0;
  filehandle = teems_fopen(newtabfile1,"r");
  if (filehandle==NULL) return -1;
  fout = teems_fopen(newtabfile,"w");
  i=0;
  while (fgets(line,TABREADLINE,filehandle)) {
    formula_qualifiers_normalize(line);
    sum_dedup_indices(line);
    l1=str_find_ci(line,"equation");
    l2=str_find_ci(line,"formula");
    l3=str_find_ci(line,"read");
    l4=str_find_ci(line,"update");
    /* assertions too (manual 10.14): a quoted element argument used to
       fall through to the executor, which read element 0 silently */
    l5=str_find_ci(line,"assertion");
    if (l3==0&&str_count_ci(line,"\"")<3) l3=-1;
    if (l1==0||l2==0||l3==0||l4==0||l5==0) {
      /* p_X and c_X name the linear variable of a levels variable X
         only (manual 9.2.2); in equations and updates c_X becomes the
         p_X column token. A declared coefficient or variable of that
         name is itself (vetting dynamic G8 / small S11). */
      if (l1==0||l2==0||l4==0||l5==0) {
        reduce_marks(line);
        if ((l1==0||l4==0)&&(strchr(line,SUM_MARK_PROD)!=NULL||strchr(line,SUM_MARK_MAXS)!=NULL||strchr(line,SUM_MARK_MINS)!=NULL)) {
          if (n_pc_vars==0&&pc_vars==NULL&&pc_kinds_build(newtabfile1)<0) { fclose(filehandle); fclose(fout); return -1; }
          if (reduce_linvar_check(line)<0) { fclose(filehandle); fclose(fout); return -1; }
        }
      }
      if ((l1==0||l2==0||l4==0||l5==0)&&(strstr(line,"p_")!=NULL||strstr(line,"c_")!=NULL)) {
        if (n_pc_vars==0&&pc_vars==NULL&&pc_kinds_build(newtabfile1)<0) { fclose(filehandle); fclose(fout); return -1; }
        if (pc_tokens_line(line,newtabfile1,l1==0||l4==0)<0) { fclose(filehandle); fclose(fout); return -1; }
      }
      /* normalize `:`-condition segments (manual 11.4.11) BEFORE the
         singleton-subset transform below: the GEMPACK `EQ` spelling
         becomes `=`, and a quoted element RHS (`= "ele"`) is
         unquoted+lowercased -- otherwise the transform treats the
         element as an index and rewrites it into a singleton subset
         (the former M3 deferral) */
      cond_segment_normalize(line);
      label_quotes_mask(line);
      strcpy(line1,line);
      n=strchr(line1,'\"');
      if (n!=NULL) {
        readline[0]='\0';
        readline1[0]='\0';
        ndedup=0;
        while (n!=NULL) {
          strncpy(line2,line1,n-line1);
          line2[n-line1]='\0';
          /* the owner of the quoted argument is the call whose '(' is
             still open at the quote, and the argument position is the
             comma count at that depth -- strrchr took the last '(' and
             mis-attributed a literal following a mapping call,
             X(MAP(i),"el") (manual 11.9.5) */
          {
            char *sc;
            int sdepth=0;
            n1=NULL;
            setindx=0;
            for (sc=line2; *sc!='\0'; sc++) {
              if (*sc=='(') { if (sdepth==0) { n1=sc; setindx=0; } sdepth++; }
              else if (*sc==')') { if (sdepth>0) sdepth--; if (sdepth==0) n1=NULL; }
              else if (*sc==','&&sdepth==1) setindx++;
            }
            /* nested: the innermost still-open '(' */
            if (n1!=NULL&&sdepth>1) {
              int d=0;
              setindx=0;
              for (sc=line2; *sc!='\0'; sc++) {
                if (*sc=='(') { d++; if (d==sdepth) { n1=sc; setindx=0; } }
                else if (*sc==')') d--;
                else if (*sc==','&&d==sdepth) setindx++;
              }
            }
          }
          /* any Formula, not only one spelling "(by_elements)" exactly:
             a qualifier list or a (default=...) by_elements formula is
             the same statement (vetting S5); a stray literal that is no
             mapping value fails in the executor with its own message */
          if (n1==NULL&&l2==0) {
            /* a quoted element outside every call: the codomain element
               a Formula (by_elements) assigns to a mapping (manual
               10.13.1), `M(u) = "ele"`. It has no owner to take a set
               from -- masked for the executor like the literals below,
               where it used to stop at the guard that follows */
            char *q2=strchr(n+1,'\"');
            if (q2==NULL) {
              errmsg("Error: unterminated element literal in TAB file: %s\n",line);
              return -1;
            }
            *n='\001';
            *q2='\001';
            strcpy(line,line1);
            strcpy(line2,line1);
            n=strchr(line1,'\"');
            continue;
          }
          if (n1==NULL) {
            errmsg("Error: malformed indexed expression in TAB file: %s\n",line);
            return -1;
          }
          /* the owner is the identifier run immediately before that
             '(' -- scanning back to the last blank or operator missed
             an owner after a bracket or glued to a quantifier,
             "*(p(c,"dom")" and "(all,c,mar)sales(c,"dom")", and the
             literal passed through unlowered (potential-models vetting
             SW1, 2026-09-25). A '$' prefix stays with the name so a
             $POS call never resolves to a declaration. */
          {
            int oe=(int)(n1-line2),os=oe;
            while (os>0&&cond_is_namec(line2[os-1])) os--;
            if (os>0&&line2[os-1]=='$') os--;
            /* no identifier before the '(' (a qualifier group such as
               "(initial, write ... header "COST")"): not an owner --
               ")(" used to match the first declaration with two
               quantifier groups (segfaults.md #4a) */
            if (oe-os+2>(int)sizeof(varname)) {
              errmsg("Error: TAB statement too complex in tab_preprocess: %s\n",line);
              return -1;
            }
            memcpy(varname,&line2[os],oe-os+1);
            varname[oe-os+1]='\0';
          }
          n1=varname;
          if (!cond_is_namec(varname[0])||tab_read_set_name(newtabfile1,n1,setindx,setname)==-1) {
            /* the quoted element is not an argument of a declared
               coefficient or variable -- a $POS("el",S) literal or a
               mapping's domain/codomain element in a Formula (manual
               11.5.6, 10.13.1): leave it for the statement's executor.
               The pair is masked so the scan moves on and restored on
               output. */
            char *q2=strchr(n+1,'\"');
            if (q2==NULL) {
              errmsg("Error: unterminated element literal in TAB file: %s\n",line);
              return -1;
            }
            *n='\001';
            *q2='\001';
            strcpy(line,line1);
            strcpy(line2,line1);
            n=strchr(line1,'\"');
            if (l3==0&&str_count_ci(line1,"\"")<3) n=NULL;
            continue;
          }
          strcpy(indx,"i");
          sprintf(indx1, "%d",i);
          strcat(indx,indx1);
          strcpy(newset,"sub_");
          strcat(newset,indx);
          i1=n-line1;
          readitem=strtok(n,"\"");
          if (readitem==NULL||cmf_strcpy_bounded(setelement,readitem,sizeof(setelement))) {
            errmsg("Error: malformed indexed expression in TAB file: %s\n",line);
            return -1;
          }
          k1=0;
          while(setelement[k1]!='\0') {
            setelement[k1]=tolower((int)setelement[k1]);
            k1++;
          }
          if (cmf_strcpy_bounded(indx1,"\"",sizeof(indx1)) ||
              cmf_strcat_bounded(indx1,readitem,sizeof(indx1)) ||
              cmf_strcat_bounded(indx1,"\"",sizeof(indx1))) {
            errmsg("Error: TAB statement too complex in tab_preprocess: %s\n",line);
            return -1;
          }
          /* one synthesized set per (owning set, element) PER STATEMENT:
             a literal repeated in a statement used to get a fresh set
             and quantifier each time, and E_regx0_A of ORANI-G (ten
             "dom" arguments, thirty-six after condensation) overran the
             statement frame so a sum index bound to the wrong set
             (2026-09-24). A repeat reuses the earlier index name. */
          {
            int dk,dhit=-1;
            for (dk=0; dk<ndedup; dk++)
              if (strcmp(dedup_set[dk],setname)==0&&strcmp(dedup_ele[dk],setelement)==0) { dhit=dk; break; }
            if (dhit>=0) {
              if (str_subst_first_bounded(line,indx1,dedup_idx[dhit],sizeof(line))) {
                errmsg("Error: TAB statement too complex in tab_preprocess: %s\n",line);
                return -1;
              }
              strcpy(line2,line);
              strcpy(line1,line);
              n=strchr(line1,'\"');
              if (l3==0&&str_count_ci(line1,"\"")<3) n=NULL;
              continue;
            }
            if (ndedup<MAXLITDEDUP) {
              strcpy(dedup_set[ndedup],setname);
              strcpy(dedup_ele[ndedup],setelement);
              strcpy(dedup_idx[ndedup],indx);
              ndedup++;
            }
          }
          if (cmf_strcat_bounded(readline,"set ",sizeof(readline)) ||
              cmf_strcat_bounded(readline,newset,sizeof(readline)) ||
              cmf_strcat_bounded(readline," (",sizeof(readline)) ||
              cmf_strcat_bounded(readline,setelement,sizeof(readline)) ||
              cmf_strcat_bounded(readline,");\n",sizeof(readline)) ||
              cmf_strcat_bounded(readline1,"subset sub_",sizeof(readline1)) ||
              cmf_strcat_bounded(readline1,indx,sizeof(readline1)) ||
              cmf_strcat_bounded(readline1," is subset of ",sizeof(readline1)) ||
              cmf_strcat_bounded(readline1,setname,sizeof(readline1)) ||
              cmf_strcat_bounded(readline1," ;\n",sizeof(readline1))) {
            errmsg("Error: TAB statement too complex in tab_preprocess: %s\n",line);
            return -1;
          }
          /* the index name replacing the quoted element can be longer
             than it (a short element with a multi-digit index), so this
             substitution grows the statement: bounded, or a near-full
             statement overruns `line` (fuzz batch 13) */
          if (str_subst_first_bounded(line,indx1,indx,sizeof(line))) {
            errmsg("Error: TAB statement too complex in tab_preprocess: %s\n",line);
            return -1;
          }
          /* a scalar-target Formula, "X = sum{r,REG, V("ele",r)}": a
             leading quantifier over the singleton set gave the formula a
             frame its scalar target does not have and the executor bound
             the literal's index to the wrong set. The literal is bound
             by wrapping the right-hand side in a sum over the singleton
             instead (a one-element sum is the term itself), 2026-09-24 */
          if (l2==0&&str_find_ci(line,"(all,")==-1&&formula_lhs_scalar(line)) {
            char *eq=strchr(line,'='),*sc=strrchr(line,';');
            if (eq!=NULL&&sc!=NULL&&sc>eq) {
              *eq='\0';
              *sc='\0';
              if (snprintf(line2,sizeof(line2),"%s= sum(%s,%s, %s);\n",line,indx,newset,eq+1)>=(int)sizeof(line2)) {
                errmsg("Error: TAB statement too complex in tab_preprocess: %s\n",line);
                return -1;
              }
              strcpy(line,line2);
              strcpy(line1,line);
              n=strchr(line1,'\"');
              if (l3==0&&str_count_ci(line1,"\"")<3) n=NULL;
              i++;
              continue;
            }
          }
          /* a qualifier glued to the quantifiers, "formula
             (initial)(all,c,comm)(all,r,reg: ...)", was one token: the
             synthesized quantifier landed inside the condition */
          {
            char *g=line;
            while ((g=strstr(g,")(all,"))!=NULL) {
              char *o=g;
              int d=0;
              for (; o>=line; o--) {
                if (*o==')') d++;
                else if (*o=='(') { d--; if (d==0) break; }
              }
              if (o>=line&&strncmp(o,"(all,",5)!=0&&strlen(line)+2<TABREADLINE) {
                memmove(g+2,g+1,strlen(g+1)+1);
                g[1]=' ';
              }
              g+=2;
            }
          }
          readitem1=strtok(line," ");
          if (readitem1==NULL||
              cmf_strcpy_bounded(line2,readitem1,sizeof(line2)) ||
              cmf_strcat_bounded(line2," ",sizeof(line2))) {
            errmsg("Error: malformed indexed expression in TAB file: %s\n",line);
            return -1;
          }
          readitem1=strtok(NULL," ");
          if (readitem1==NULL) {
            errmsg("Error: malformed indexed expression in TAB file: %s\n",line);
            return -1;
          }
          strcpy(indx2,indx);
          strcat(indx2,",");
          if(strstr(readitem1,indx2)==NULL) {
            /* a qualifier group after the keyword, "equation (levels)
               e_x ...", "formula (initial) ...", rides with the keyword:
               the synthesized quantifier used to land between it and
               the name, and the levels transform found no name (GTAP-W
               E_TAXBAS*, 2026-09-24) */
            int fe=(l2==0&&(strstr(line1,"&equation")!=NULL||strstr(line1,"& equation")!=NULL));
            if(l1==0||l2==0||l4==0) {
              while ((readitem1[0]=='('&&strncmp(readitem1,"(all,",5)!=0)||
                     (fe&&(strcmp(readitem1,"&")==0||strcmp(readitem1,"&equation")==0||strcmp(readitem1,"equation")==0))) {
                if (cmf_strcat_bounded(line2,readitem1,sizeof(line2)) ||
                    cmf_strcat_bounded(line2," ",sizeof(line2))) {
                  errmsg("Error: TAB statement too complex in tab_preprocess: %s\n",line);
                  return -1;
                }
                readitem1=strtok(NULL," ");
                if (readitem1==NULL) {
                  errmsg("Error: malformed indexed expression in TAB file: %s\n",line);
                  return -1;
                }
              }
            }
            if(l1==0||fe){//if(l1==0&&strchr(line2,'(')!=NULL) {
              if (cmf_strcat_bounded(line2,readitem1,sizeof(line2)) ||
                  cmf_strcat_bounded(line2," ",sizeof(line2))) {
                errmsg("Error: TAB statement too complex in tab_preprocess: %s\n",line);
                return -1;
              }
              readitem1=strtok(NULL," ");
              if (readitem1==NULL) {
                errmsg("Error: malformed indexed expression in TAB file: %s\n",line);
                return -1;
              }
            }
            /* the quantifier used to be glued to the body, so a
               body-leading variable "(all,i5,sub_i5)contbot" was
               never prefixed by tab_write_variables (ORANI-G,
               2026-09-24): a space separates them */
            if (cmf_strcat_bounded(line2," (all,",sizeof(line2)) ||
                cmf_strcat_bounded(line2,indx,sizeof(line2)) ||
                cmf_strcat_bounded(line2,",",sizeof(line2)) ||
                cmf_strcat_bounded(line2,newset,sizeof(line2)) ||
                cmf_strcat_bounded(line2,") ",sizeof(line2)) ||
                cmf_strcat_bounded(line2,readitem1,sizeof(line2)) ||
                cmf_strcat_bounded(line2," ",sizeof(line2))) {
              errmsg("Error: TAB statement too complex in tab_preprocess: %s\n",line);
              return -1;
            }
          } else {
            if (cmf_strcat_bounded(line2," (all,",sizeof(line2)) ||
                cmf_strcat_bounded(line2,indx,sizeof(line2)) ||
                cmf_strcat_bounded(line2,",",sizeof(line2)) ||
                cmf_strcat_bounded(line2,newset,sizeof(line2)) ||
                cmf_strcat_bounded(line2,")",sizeof(line2)) ||
                cmf_strcat_bounded(line2," ",sizeof(line2)) ||
                cmf_strcat_bounded(line2,readitem1,sizeof(line2))) {
              errmsg("Error: TAB statement too complex in tab_preprocess: %s\n",line);
              return -1;
            }
          }
          readitem1=strtok(NULL,"\n");
          if(readitem1!=NULL&&cmf_strcat_bounded(line2,readitem1,sizeof(line2))) {
            errmsg("Error: TAB statement too complex in tab_preprocess: %s\n",line);
            return -1;
          }
          if (cmf_strcat_bounded(line2,"\n",sizeof(line2))) {
            errmsg("Error: TAB statement too complex in tab_preprocess: %s\n",line);
            return -1;
          }
          strcpy(line,line2);
          strcpy(line1,line);
          n=strchr(line1,'\"');
          if (l3==0&&str_count_ci(line1,"\"")<3) n=NULL;
          i++;
        }
        str_replace_char_all(readline,'\001','\"');
        str_replace_char_all(readline1,'\001','\"');
        str_replace_char_all(line2,'\001','\"');
        fprintf(fout,"%s%s%s",readline,readline1,line2);
      } else {
        str_replace_char_all(line,'\001','\"');
        fprintf(fout,"%s",line);
      }
    } else fprintf(fout,"%s",line);
  }
  fclose(filehandle);
  fclose(fout);
  remove(newtabfile1);
  return 1;
}

int outputs_write_csv(char *filename, char *newdatlogname, char *newdatfile,set_def *sets,dim_t nset, set_element *set_elems,array_def *coefs,offset_t ncof,offset_t ncofele,array_def *vars,offset_t nvar,offset_t nvarele, elem_value *elem_vals) {
  FILE * filehandle,*fout;
  char line[TABREADLINE]="\0",*readline,comsyntax[TABREADLINE],longname[TABREADLINE],varname[NAMESIZE],*vname1,header[NAMESIZE],setsize[DATREADLINE],tempname[NAMESIZE];
  filehandle = teems_fopen(filename,"r");
  long int i,n,j,j1,innerloop,outerloop,l,indx;
  long int setindx[MAXVARDIM],antidim[MAXVARDIM];
  strcpy(comsyntax,"to file ");
  strcat(comsyntax,newdatlogname);
  strcat(comsyntax," header \"");
  n=strlen(comsyntax);
  fout = fopen(newdatfile,"w");
  /* named failure instead of writing through a NULL handle; the caller
     goes on to the remaining outputs and the errmsg count makes the
     exit status nonzero */
  if(filehandle==NULL) {
    errmsg("Error: cannot open %s for reading\n",filename);
    if(fout!=NULL)fclose(fout);
    return -1;
  }
  if(fout==NULL) {
    errmsg("Error: cannot open %s for writing: %s (the solver runs as uid %d)\n",newdatfile,strerror(errno),(int)getuid());
    fclose(filehandle);
    return -1;
  }
  i=0;
  while (fgets(line,TABREADLINE,filehandle)) {
    if(strncmp(line,"write",5)==0&&strstr(line,comsyntax)!=NULL) {
      if(strstr(line,"(set)")!=NULL) {
        i=str_find_ci(line," ");
        readline=line+i+1;
        i=str_find_ci(readline," ");
        readline=readline+i+1;
        i=str_find_ci(readline," ");
        strncpy(varname,readline,i);
        varname[i]='\0';
        i=str_find_ci(line,comsyntax);
        readline=line+i+n;
        i=str_find_ci(readline,"\"");
        strncpy(header,readline,i);
        header[i]='\0';

        i=str_find_ci(line,"longname \"");
        if(i<0) longname[0]='\0';
        else {
          readline=line+i+10;
          i=str_find_ci(readline,"\"");
          if(i<0) i=0;
          strncpy(longname,readline,i);
          longname[i]='\0';
        }
        setsize[0]='\0';
        for (i=0; i<nset; i++) {
          if (strcmp(sets[i].setname,varname)==0) {
            sprintf(setsize, "%d", sets[i].size);
            strcat(setsize," Strings Length 12 Header \"");
            strcat(setsize,header);
            strcat(setsize,"\" LongName \"");
            strcat(setsize,longname);
            strcat(setsize,"\";\n");
            fprintf(fout,"%s",setsize);
            for (j=0; j<sets[i].size; j++)fprintf(fout,"%s\n",set_elems[sets[i].offset+j].setele);
            fprintf(fout,"\n");
          }
        }
      } else {
        int wbyele=0;
        i=str_find_ci(line," ");
        readline=line+i+1;
        if(strncmp(readline,"(by_elements) ",14)==0) { wbyele=1; readline+=14; }
        i=str_find_ci(readline," ");
        strncpy(varname,readline,i);
        varname[i]='\0';
        i=str_find_ci(line,comsyntax);
        if(i==-1)continue;
        readline=line+i+n;
        i=str_find_ci(readline,"\"");
        strncpy(header,readline,i);
        header[i]='\0';

        i=str_find_ci(line,"longname \"");
        if(i<0) longname[0]='\0';
        else {
          readline=line+i+10;
          i=str_find_ci(readline,"\"");
          if(i<0) i=0;
          strncpy(longname,readline,i);
          longname[i]='\0';
        }
        setsize[0]='\0';
        for (i=0; i<ncof; i++) {
          vname1= strtok(coefs[i].cofname,"(");

          if (strcmp(vname1,varname)==0) {
            if(coefs[i].size==0) {
              innerloop=1;
              outerloop=0;
            }
            else if(coefs[i].size==1) {
              innerloop=sets[coefs[i].setid[0]].size;
              outerloop=0;
            }
            else if(coefs[i].size==2) {
              innerloop=sets[coefs[i].setid[0]].size*sets[coefs[i].setid[1]].size;
              outerloop=1;
            }
            else {
              innerloop=sets[coefs[i].setid[0]].size*sets[coefs[i].setid[1]].size;
              outerloop=1;
              for(j=2; j<coefs[i].size; j++)outerloop*=sets[coefs[i].setid[j]].size;
            }
            for(j=0; j<coefs[i].size; j++) {
              sprintf(tempname, "%d", sets[coefs[i].setid[j]].size);
              strcat(setsize,tempname);
              strcat(setsize," ");
            }
            if(coefs[i].size==0){
              sprintf(tempname, "%d", 1);
              strcat(setsize,tempname);
              strcat(setsize," ");
            }
            strcat(setsize,"Real SpreadSheet Header \"");
            strcat(setsize,header);
            strcat(setsize,"\" LongName \"");
            strcat(setsize,longname);
            strcat(setsize,"\";\n");
            fprintf(fout,"%s",setsize);
            indx=0;
            if(coefs[i].size<2) {
              for(j=0; j<innerloop; j++) {
                fprintf(fout,"%f\n",elem_vals[coefs[i].offset+j].value);
              }
              fprintf(fout,"\n");
            } else {
              antidim[2]=1;
              for (l=3; l<coefs[i].size; l++){
                antidim[l]=antidim[l-1]*sets[coefs[i].setid[l-1]].size;
              }
              for(j1=0; j1<outerloop; j1++) {
                indx=j1;
                for (l=coefs[i].size-1; l>1; l--) {
                  setindx[l]=indx/antidim[l];
                  indx-=setindx[l]*antidim[l];
                }
                for(j=0; j<innerloop; j++) {
                  setindx[0]=j/sets[coefs[i].setid[1]].size;
                  setindx[1]=j-sets[coefs[i].setid[1]].size*setindx[0];
                  indx=0;
                  for (l=0; l<coefs[i].size; l++)indx+=coefs[i].strides[l]*setindx[l];
                  if(setindx[1]==sets[coefs[i].setid[1]].size-1){
                    fprintf(fout,"%f\n",elem_vals[coefs[i].offset+indx].value);
                  }else{
                    fprintf(fout,"%f,",elem_vals[coefs[i].offset+indx].value);
                  }
                }
                fprintf(fout,"\n");
              }
            }
            break;
          }
        }
        /* a mapping (manual 11.9.10): its element numbers, or with
           (by_elements) its element names, one per domain element */
        {
          dim_t mm,e;
          for (mm=0; mm<teems_nmap; mm++) if (strcmp(teems_maps[mm].mapname,varname)==0) break;
          if (mm<teems_nmap) {
            map_def *md=&teems_maps[mm];
            dim_t n1=sets[md->fromset].size;
            if (!mapping_ready(mm)) {
              errmsg("Error: mapping %s is written before all of its values are assigned (manual 11.9.10)\n",md->mapname);
              fclose(filehandle);
              fclose(fout);
              return -1;
            }
            if (wbyele) {
              size_t len=1;
              for (e=0; e<n1; e++) {
                size_t l1=strlen(set_elems[sets[md->toset].offset+md->values[e]].setele);
                if (l1>len) len=l1;
              }
              fprintf(fout,"%ld Strings Length %ld Header \"%s\" LongName \"%s\";\n",(long)n1,(long)len,header,longname);
              for (e=0; e<n1; e++) fprintf(fout,"%s\n",set_elems[sets[md->toset].offset+md->values[e]].setele);
            } else {
              fprintf(fout,"%ld Integer SpreadSheet Header \"%s\" LongName \"%s\";\n",(long)n1,header,longname);
              for (e=0; e<n1; e++) fprintf(fout,"%ld\n",(long)md->values[e]+1);
            }
            fprintf(fout,"\n");
          }
        }
      }
    }
  }
  fclose(filehandle);
  fclose(fout);
  return 1;
}

/* an IF condition may read coefficients and levels variables, not
   linear variables (manual 11.4.5): in an equation or update the
   condition would be evaluated on the substituted columns. The linear
   references carry their p_ prefix by now. */
static int if_cond_linear_var(const char *line, array_def *vars, offset_t nvar, char *bad) {
  const char *p=line;
  while ((p=strstr(p,"if"))!=NULL) {
    const char *o=p+2,*c;
    int d=0,q=0;
    if ((p>line&&ifin_namec(p[-1]))) { p+=2; continue; }
    while (*o==' ') o++;
    if (*o!='('&&*o!='['&&*o!='{') { p+=2; continue; }
    for (c=o; *c!='\0'; c++) {
      if (*c=='"') { q=!q; continue; }
      if (q) continue;
      if (*c=='('||*c=='['||*c=='{') d++;
      else if (*c==')'||*c==']'||*c=='}') { d--; if (d==0) break; }
      else if (*c==','&&d==1) break;
      if (strncmp(c,"p_",2)==0&&(c==o||!ifin_namec(c[-1]))) {
        char nm[NAMESIZE];
        int k=0;
        offset_t v;
        const char *t=c+2;
        while (ifin_namec(*t)&&k<NAMESIZE-1) nm[k++]=*t++;
        nm[k]='\0';
        for (v=0; v<nvar; v++) if (strcmp(nm,vars[v].cofname)==0) { strcpy(bad,nm); return 1; }
      }
    }
    p+=2;
  }
  return 0;
}

static int pfx_depth(const char *nm) {
  int d=0;
  while (nm[0]=='p'&&nm[1]=='_') { d++; nm+=2; }
  return d;
}

int tab_write_variables(char *filename, char *newtabfile,array_def *vars,offset_t nvar) {
  FILE * filehandle,*fout;
  char line[TABREADLINE+1]="\0",*p;//,nvarname[nvar][NAMESIZE+2],*p;//,line1[DATREADLINE];//,*ne,*np;//,*n2;
  filehandle = teems_fopen(filename,"r");
  offset_t i,n,j,l,l1,linelght;
  int lvar;
  if (filehandle==NULL) return -1;
  /* every linear variable gets the p_ column marker, p_-leading names
     too (a linear zz and a linear p_zz are distinct, manual 9.2.2);
     the deepest p_ nesting goes first so a freshly marked zz is never
     taken for a declared p_zz */
  offset_t *ord=(offset_t *)malloc((nvar>0?nvar:1)*sizeof(offset_t)),oi;
  int dmax=0,dd;
  if (ord==NULL) { fclose(filehandle); return -1; }
  for (i=0; i<nvar; i++) if ((dd=pfx_depth(vars[i].cofname))>dmax) dmax=dd;
  for (oi=0,dd=dmax; dd>=0; dd--) for (i=0; i<nvar; i++) if (pfx_depth(vars[i].cofname)==dd) ord[oi++]=i;
  fout = teems_fopen(newtabfile,"w");
  while (fgets(line,TABREADLINE,filehandle)) {
    int eqpos=str_find_ci(line,"equation ");
    int updpos=str_find_ci(line,"update ");
    /* conditional Update quantifiers are evaluated by updates_apply
       (vetting S9; this used to be a fatal) */
    if(eqpos>-1||updpos>-1) {
      linelght=strlen(line);
      for (oi=0; oi<nvar; oi++) {
        i=ord[oi];
        p=strchr(line,';');
        if(p==NULL) break;
        line[p-line+1]='\n';
        line[p-line+2]='\0';
        linelght=strlen(line);
        n=str_count_ci(line,vars[i].cofname);
        lvar=strlen(vars[i].cofname);
        l=0;
        for (j=0; j<n; j++) {
          l1=str_find_ci(&line[l],vars[i].cofname);
          if(l1<0) break;
          l=l+l1;
          /* a scalar variable closing a group, "... - pxwwld)", was not prefixed: the term was bound as a value (its column dropped) -- ')' and '}' close the follower set (2026-09-06) */
          if(vars[i].level_par==false) if(line[l+lvar]==')'||line[l+lvar]=='}'||line[l+lvar]==' '||line[l+lvar]=='('||line[l+lvar]=='+'||line[l+lvar]=='-'||line[l+lvar]=='*'||line[l+lvar]=='/'||line[l+lvar]=='^'||line[l+lvar]==']'||line[l+lvar]==','||line[l+lvar]==';'||line[l+lvar]=='=')if(l==0||line[l-1]==' '||line[l-1]=='+'||line[l-1]=='-'||line[l-1]=='*'||line[l-1]=='/'||line[l-1]=='^'||line[l-1]=='['||line[l-1]=='('||line[l-1]==')'||line[l-1]==','||line[l-1]=='=') { /* ')' : a quantifier written flush against the body, "(all,c,com)p(c)" (2026-09-24) */
                /* each p_ prefix grows the statement by 2 and nothing
                   bounded the growth (fuzz batch 13 stack overflow);
                   the move carries the terminator along */
                if(linelght+3>TABREADLINE) {
                  errmsg("Error: equation too long after prefixing its variables (exceeds %d chars)\n",TABREADLINE);
                  fclose(filehandle);
                  fclose(fout);
                  MPI_Abort(PETSC_COMM_WORLD,1);
                  return -1;
                }
                memmove(&line[l+2],&line[l],linelght-l+1);
                line[l]='p';
                line[l+1]='_';
                l=l+2;
                linelght+=2;
              }
          l=l+strlen(vars[i].cofname);
        }
      }
      /* one-shot lowering of mapping calls in Equation statements
         (design doc M2b): every downstream equation consumer -- the
         ordering scans, preallocation, the statement builder, the
         backsolve validator -- reads this rewritten file, so the
         nested MAP(i) index form is rewritten to the flat map~i token
         exactly once, here.  Updates keep their named fatal
         (mapping_use_guards). */
      if((eqpos==0||eqpos==1)&&teems_nmap>0) mapping_lower_calls(line);
      /* equation-level quantifier conditions prune ROWS -- that
         changes VecSize/eq_addr/closure squareness and is not
         supported; before M3 the ':' made the set lookup miss and the
         equation silently expanded over sets[0] */
      if(eqpos==0||eqpos==1) {
        offset_t qk=0,qf;
        while((qf=str_find_ci(&line[qk],"(all,"))>-1) {
          for(qk=qk+qf+5; line[qk]!='\0'&&line[qk]!=')'; qk++) {
            if(line[qk]==':') {
              errmsg("Error: conditions on Equation quantifiers are not supported (row pruning); put the condition on a sum inside the equation (manual 11.4.11)\n");
              fclose(filehandle);
              fclose(fout);
              MPI_Abort(PETSC_COMM_WORLD,1);
              return -1;
            }
          }
        }
      }
    }
    /* compound and index IF conditions (manual 11.4.5, 11.4.11) become
       one numeric test here, ahead of every reader */
    int fopos=str_find_ci(line,"formula"),aspos=str_find_ci(line,"assertion");
    if ((eqpos==0||eqpos==1||updpos==0||updpos==1)&&strstr(line,"if")!=NULL) {
      char badv[NAMESIZE];
      if (if_cond_linear_var(line,vars,nvar,badv)) {
        errmsg("Error: an IF condition reads variable %s; conditions take coefficients and levels variables, not linear variables (manual 11.4.5): %.200s\n",badv,line);
        fclose(filehandle);
        fclose(fout);
        MPI_Abort(PETSC_COMM_WORLD,1);
        return -1;
      }
    }
    if (strstr(line,"if")!=NULL&&(eqpos==0||eqpos==1||updpos==0||updpos==1||fopos==0||fopos==1||aspos==0||aspos==1)) {
      if (!cond_text_lower(line,TABREADLINE,NULL,0,"an IF condition")) {
        fclose(filehandle);
        fclose(fout);
        MPI_Abort(PETSC_COMM_WORLD,1);
        return -1;
      }
    }
    fprintf(fout,"%s",line);
  }
  fclose(filehandle);
  fclose(fout);
  free(ord);
  return 1;
}


/* Declaration index for tab_read_set_name: the preprocessed TAB's
   variable and coefficient statements, read once per file (spaces
   already stripped) and searched in memory. Before it, every quoted
   element in the TAB reopened the file and rescanned every statement
   through 20 KB line buffers cleared per statement -- quadratic in
   the statement count (fuzz batch 13: 15 s on a 100 KB input of tiny
   statements). The lookup logic below is the original's, run over the
   arrays instead of the file. */
static char *decl_index_file=NULL;
static char **decl_var=NULL,**decl_cof=NULL;
static int n_decl_var=0,n_decl_cof=0;

static int decl_index_read(char *filename, char *commsyntax, char ***out, int *nout) {
  FILE *filehandle=teems_fopen(filename,"r");
  char line[TABREADLINE+1];
  offset_t lsize=TABREADLINE+1;
  int n=0,cap=0;
  char **arr=NULL;
  if (filehandle==NULL) return -1;
  while (tab_next_statement(commsyntax,filehandle,line,lsize)) {
    str_delete_char(line,' ');
    if (n==cap) {
      cap=cap?2*cap:64;
      char **grown=(char **)realloc(arr,cap*sizeof(char *));
      if (grown==NULL) { fclose(filehandle); for (int i=0; i<n; i++) free(arr[i]); free(arr); return -1; }
      arr=grown;
    }
    arr[n]=strdup(line);
    if (arr[n]==NULL) { fclose(filehandle); for (int i=0; i<n; i++) free(arr[i]); free(arr); return -1; }
    n++;
  }
  fclose(filehandle);
  *out=arr;
  *nout=n;
  return 0;
}

void tab_decl_index_free(void) {
  int i;
  for (i=0; i<n_decl_var; i++) free(decl_var[i]);
  for (i=0; i<n_decl_cof; i++) free(decl_cof[i]);
  free(decl_var); free(decl_cof); free(decl_index_file);
  decl_var=decl_cof=NULL; decl_index_file=NULL;
  n_decl_var=n_decl_cof=0;
}

static int decl_index_build(char *filename) {
  if (decl_index_file!=NULL&&strcmp(decl_index_file,filename)==0) return 0;
  tab_decl_index_free();
  decl_index_file=strdup(filename);
  if (decl_index_file==NULL) return -1;
  if (decl_index_read(filename,"variable",&decl_var,&n_decl_var)<0) { tab_decl_index_free(); return -1; }
  if (decl_index_read(filename,"coefficient",&decl_cof,&n_decl_cof)<0) { tab_decl_index_free(); return -1; }
  return 0;
}

/* one pass of the lookup over an indexed statement kind for the
   declared name varname1 (")NAME("). Returns 1 found, -1 malformed,
   0 not in this kind. */
static int tab_read_set_name_pass(char **stmts, int nstmt, char *varname1, int indx, char *setname) {
  int n,i,k;
  char indxname[NAMESIZE],line[TABREADLINE+1],line1[TABREADLINE+1],*p,tmp[TABREADLINE+1];
  for (k=0; k<nstmt; k++) {
    strcpy(line,stmts[k]);
    strcpy(line1,line);
    n=str_find_ci(line,varname1);
    if (n>-1) {
      p=strtok(line+n,"(");
      p=strtok(NULL,")");
      if (p==NULL||strlen(p)+1>=sizeof(tmp)) return -1;
      strcpy(tmp,p);
      strcat(tmp,",");
      p=NULL;
      for (i=0; i<indx+1; i++) {
        if(i==0) p=strtok(tmp,",");
        else p=strtok(NULL,",");
        if (p==NULL) return -1;
      }
      if (strlen(p)+2>=sizeof(indxname)) return -1;
      /* ",idx," -- the bare ",idx" prefix matched inside qualifier
         groups such as (linear,change) for an index named c */
      strcpy(indxname,",");
      strcat(indxname,p);
      strcat(indxname,",");
      strcpy(line,line1);
      n=str_find_ci(line,indxname);
      if (n<0) return -1;
      p=strtok(line+n,",");
      p=strtok(NULL,")");
      if (p==NULL||strlen(p)>=NAMESIZE) return -1;
      strcpy(setname,p);
      return 1;
    }
  }
  return 0;
}

int tab_read_set_name(char *filename, char *varname, int indx, char *setname) {
  int r;
  char varname1[NAMESIZE+2];
  /* an owner is a declared name: "(" alone matched ")(" in the first
     declaration with two quantifier groups */
  if (!cond_is_namec(varname[0])) return -1;
  strcpy(varname1,")");
  strcat(varname1,varname);
  if (decl_index_build(filename)<0) return -1;
  r=tab_read_set_name_pass(decl_var,n_decl_var,varname1,indx,setname);
  if (r!=0) return r;
  r=tab_read_set_name_pass(decl_cof,n_decl_cof,varname1,indx,setname);
  if (r!=0) return r;
  /* p_X / c_X with no declaration of its own is the linear variable of
     a levels X (manual 9.2.2); a declared p_X or c_X coefficient or
     variable resolves above and its literals are checked like any
     other's */
  if ((tolower((unsigned char)varname[0])=='p'||tolower((unsigned char)varname[0])=='c')&&varname[1]=='_') {
    strcpy(varname1,")");
    strcat(varname1,varname+2);
    r=tab_read_set_name_pass(decl_var,n_decl_var,varname1,indx,setname);
    if (r!=0) return r;
  }
  return -1;
}




/* Validate every Default statement up front (manual 10.19; audit A6),
   one scan of the preprocessed TAB (one statement per line,
   lowercased). Supported -- applied positionally by the readers:
   Coefficient parameter/non_parameter, Variable linear/levels/change/
   percent_change, Formula initial/always, Equation linear/levels (the
   levels transform qualifies unqualified equations positionally).
   Equation add_homotopy[=name]/not_add_homotopy are applied by the
   levels transform (manual 26.7.5). Fatal: Coefficient
   lower_bound/upper_bound defaults (single bound slot, audit A9), and
   any unknown keyword or value. Returns 0 ok, -1 fatal. */
int tab_defaults_validate(char *fname) {
  FILE *f;
  char line[TABREADLINE],val[NAMESIZE];
  int bad=0;
  f=teems_fopen(fname,"r");
  if(f==NULL)return 0;
  while(fgets(line,TABREADLINE,f)) {
    if(strstr(line,"(default")==NULL)continue;
    tab_default_value(line,val);
    if(strncmp(line,"coefficient",11)==0) {
      if(strcmp(val,"parameter")==0||strcmp(val,"non_parameter")==0)continue;
      if(strncmp(val,"lower_bound",11)==0||strncmp(val,"upper_bound",11)==0)
        errmsg("Error: Coefficient (default=%s): bound defaults are not supported\n",val);
      else errmsg("Error: unknown Coefficient default '%s'\n",val);
    } else if(strncmp(line,"variable",8)==0) {
      if(strcmp(val,"linear")==0||strcmp(val,"levels")==0||strcmp(val,"change")==0||strcmp(val,"percent_change")==0)continue;
      errmsg("Error: unknown Variable default '%s'\n",val);
    } else if(strncmp(line,"formula",7)==0) {
      if(strcmp(val,"initial")==0||strcmp(val,"always")==0)continue;
      errmsg("Error: unknown Formula default '%s'\n",val);
    } else if(strncmp(line,"equation",8)==0) {
      if(strcmp(val,"linear")==0||strcmp(val,"levels")==0||strcmp(val,"not_add_homotopy")==0)continue;
      if(strcmp(val,"add_homotopy")==0||(strncmp(val,"add_homotopy=",13)==0&&val[13]!='\0'))continue;
      errmsg("Error: unknown Equation default '%s'\n",val);
    } else {
      errmsg("Error: Default statements apply only to Coefficient/Variable/Formula/Equation: %s",line);
    }
    bad=1;
  }
  fclose(f);
  return bad?-1:0;
}

/* --- PostSim scope helpers (Tier 0 residuals; manual 12.2.1-12.2.3) --- */

/* declared name of a preprocessed declaration line: skip the keyword,
   any (qualifier)/(quantifier) groups and blanks; the name is the next
   identifier run (works for set/subset/coefficient/file forms) */
static void ps_decl_name(char *line, int kwlen, char *out) {
  char *p=line+kwlen;
  int k=0;
  out[0]='\0';
  for(;;) {
    while(*p==' ')p++;
    if(*p!='(')break;
    while(*p!=')'&&*p!='\0')p++;
    if(*p==')')p++;
  }
  while((isalnum((int)*p)||*p=='_'||*p=='@'||*p==MAPMARK)&&k<NAMESIZE-1)out[k++]=*p++;
  out[k]='\0';
}

/* word-boundary identifier search (lines and names are lowercase) */
static int line_has_ident(char *line, char *name) {
  char *p=line;
  size_t n=strlen(name);
  if(n==0)return 0;
  while((p=strstr(p,name))!=NULL) {
    int lb=(p==line)?0:(isalnum((int)p[-1])||p[-1]=='_');
    int rb=(isalnum((int)p[n])||p[n]=='_');
    if(!lb&&!rb)return 1;
    p++;
  }
  return 0;
}

/* logical file name of a read statement: the token after " file "
   (terminal/text forms without one yield "") */
static void ps_read_logname(char *line, char *out) {
  char *p=strstr(line," file ");
  int k=0;
  out[0]='\0';
  if(p==NULL)return;
  p+=6;
  while(*p==' ')p++;
  while(isalnum((int)*p)||*p=='_') {
    if(k<NAMESIZE-1)out[k++]=*p;
    p++;
  }
  out[k]='\0';
}

/* Split the preprocessed TAB around POSTSIM (BEGIN)/(END) sections
   (manual 10.18, 12.2.1). Declarations (set/subset/coefficient/file
   and "read elements") stay in the ordinary file -- a single
   namespace; what separates PostSim is EXECUTION order -- while
   executables (formula/assertion/zerodivide/read) move to the _ps
   companion consumed once after the solve. Write/Display in sections
   are dropped (outputs ride the write-all dump); the forbidden
   statements are fatal. Scope rules enforced here: ordinary
   statements may not reference PostSim-declared names (12.2.1), and
   no file serves both normal and PostSim Reads (12.2.3); PostSim
   coefficient names are recorded for the 12.2.2/12.2.3 LHS and
   Read-target rules. Returns the PostSim executable count (0 = no
   sections), -1 on error. */
int tab_postsim_split(char *newtabfile, char *psfile) {
  FILE *fin,*fmain,*fps;
  char line[TABREADLINE],tmpname[TABREADLINE],nm[NAMESIZE];
  int inps=0,nps=0,found=0;
  char (*psnames)[NAMESIZE]=NULL,(*ordlogs)[NAMESIZE]=NULL,(*pslogs)[NAMESIZE]=NULL;
  int npsn=0,nordlog=0,npslog=0,k;
  fin=teems_fopen(newtabfile,"r");
  if(fin==NULL)return 0;
  while (fgets(line,TABREADLINE,fin)) {
    if(strncmp(line,"postsim (begin)",15)==0||strncmp(line,"postsim(begin)",14)==0) {
      found=1;
      break;
    }
  }
  if(!found) {
    fclose(fin);
    return 0;
  }
  /* pass A: collect the names declared inside sections (scope checks
     below) and record the PostSim coefficients for the readers */
  rewind(fin);
  while (fgets(line,TABREADLINE,fin)) {
    if(strncmp(line,"postsim (begin)",15)==0||strncmp(line,"postsim(begin)",14)==0) {
      inps=1;
      continue;
    }
    if(strncmp(line,"postsim (end)",13)==0||strncmp(line,"postsim(end)",12)==0) {
      inps=0;
      continue;
    }
    if(!inps)continue;
    nm[0]='\0';
    if(strncmp(line,"set ",4)==0)ps_decl_name(line,4,nm);
    else if(strncmp(line,"subset ",7)==0)ps_decl_name(line,7,nm);
    else if(strncmp(line,"file ",5)==0)ps_decl_name(line,5,nm);
    else if(strncmp(line,"mapping ",8)==0)ps_decl_name(line,8,nm);
    else if(strncmp(line,"coefficient ",12)==0) {
      ps_decl_name(line,12,nm);
      if(nm[0]!='\0') {
        teems_ps_coefnames=realloc(teems_ps_coefnames,(teems_ps_ncoefs+1)*sizeof(*teems_ps_coefnames));
        strcpy(teems_ps_coefnames[teems_ps_ncoefs++],nm);
      }
    }
    if(nm[0]!='\0') {
      psnames=realloc(psnames,(npsn+1)*sizeof(*psnames));
      strcpy(psnames[npsn++],nm);
    }
  }
  inps=0;
  rewind(fin);
  strcpy(tmpname,newtabfile);
  strcat(tmpname,"_o");
  fmain=teems_fopen(tmpname,"w");
  fps=teems_fopen(psfile,"w");
  if(fmain==NULL||fps==NULL) {
    errmsg("Error: cannot open PostSim split scratch files\n");
    fclose(fin);
    if(fmain!=NULL)fclose(fmain);
    if(fps!=NULL)fclose(fps);
    return -1;
  }
  while (fgets(line,TABREADLINE,fin)) {
    if(strncmp(line,"postsim (begin)",15)==0||strncmp(line,"postsim(begin)",14)==0) {
      inps=1;
      continue;
    }
    if(strncmp(line,"postsim (end)",13)==0||strncmp(line,"postsim(end)",12)==0) {
      inps=0;
      continue;
    }
    if(!inps) {
      /* scope isolation (12.2.1): PostSim names are PostSim-only.
         Set element lists and quoted element literals are not names
         (a PostSim coefficient INVESTMENT vs the element investment of
         an ordinary set): scan a copy with those masked out */
      char scan[TABREADLINE];
      {
        int si=0,inq=0,inlist=0,isset=(strncmp(line,"set ",4)==0),k3;
        /* a write's longname is free text (it may nest quotes) */
        char *lnm=(strncmp(line,"write",5)==0)?strstr(line," longname "):NULL;
        for(k3=0; line[k3]!='\0'&&si<TABREADLINE-1; k3++) {
          char ch=line[k3];
          if(lnm!=NULL&&line+k3>=lnm) break;
          if(inq) { if(ch=='"') inq=0; continue; }
          if(ch=='"') { inq=1; continue; }
          if(isset&&ch=='(') { inlist=1; continue; }
          if(isset&&inlist) { if(ch==')') inlist=0; continue; }
          scan[si++]=ch;
        }
        scan[si]='\0';
      }
      for(k=0;k<npsn;k++)if(line_has_ident(scan,psnames[k])) {
        errmsg("Error: ordinary statement references PostSim-declared name %s (manual 12.2.1): %s",psnames[k],line);
        fclose(fin);
        fclose(fmain);
        fclose(fps);
        free(psnames);
        free(ordlogs);
        free(pslogs);
        return -1;
      }
      if(strncmp(line,"read ",5)==0&&strncmp(line,"read elements",13)!=0) {
        ps_read_logname(line,nm);
        if(nm[0]!='\0') {
          ordlogs=realloc(ordlogs,(nordlog+1)*sizeof(*ordlogs));
          strcpy(ordlogs[nordlog++],nm);
        }
      }
      fputs(line,fmain);
      continue;
    }
    if(strncmp(line,"variable ",9)==0||strncmp(line,"equation ",9)==0||strncmp(line,"update ",7)==0||strncmp(line,"transfer ",9)==0||strncmp(line,"omit ",5)==0||strncmp(line,"substitute ",11)==0||strncmp(line,"backsolve ",10)==0||strncmp(line,"complementarity",15)==0||strstr(line,"(default")!=NULL) {
      errmsg("Error: statement not allowed in a PostSim section (manual 12.2.1): %s",line);
      fclose(fin);
      fclose(fmain);
      fclose(fps);
      free(psnames);
      free(ordlogs);
      free(pslogs);
      return -1;
    }
    if(strncmp(line,"read elements",13)==0) {
      /* set-definition machinery: rides the declarations */
      fputs(line,fmain);
      continue;
    }
    if(strncmp(line,"read ",5)==0) {
      ps_read_logname(line,nm);
      if(nm[0]!='\0') {
        pslogs=realloc(pslogs,(npslog+1)*sizeof(*pslogs));
        strcpy(pslogs[npslog++],nm);
      }
      fputs(line,fps);
      nps++;
      continue;
    }
    if(strncmp(line,"set ",4)==0||strncmp(line,"subset ",7)==0||strncmp(line,"coefficient ",12)==0||strncmp(line,"file ",5)==0||strncmp(line,"mapping ",8)==0) {
      /* declarations (mappings included: their formula-assigned values
         arrive in the PostSim pass, manual 10.13.1) ride the main stream */
      fputs(line,fmain);
      continue;
    }
    if(strncmp(line,"write ",6)==0||strncmp(line,"display ",8)==0) {
      /* outputs ride the write-all coefficient dump */
      continue;
    }
    if(strncmp(line,"formula ",8)==0||strncmp(line,"assertion ",10)==0||strncmp(line,"zerodivide",10)==0) {
      fputs(line,fps);
      nps++;
      continue;
    }
    errmsg("Error: unrecognized statement in a PostSim section: %s",line);
    fclose(fin);
    fclose(fmain);
    fclose(fps);
    free(psnames);
    free(ordlogs);
    free(pslogs);
    return -1;
  }
  fclose(fin);
  fclose(fmain);
  fclose(fps);
  /* 12.2.3: no file serves both normal and PostSim Reads */
  for(k=0;k<npslog;k++) {
    int k2;
    for(k2=0;k2<nordlog;k2++)if(strcmp(pslogs[k],ordlogs[k2])==0) {
      errmsg("Error: file %s is read in both the ordinary and PostSim parts (manual 12.2.3); split the data across two files\n",pslogs[k]);
      free(psnames);
      free(ordlogs);
      free(pslogs);
      return -1;
    }
  }
  free(psnames);
  free(ordlogs);
  free(pslogs);
  if(rename(tmpname,newtabfile)!=0) {
    errmsg("Error: cannot finalize the PostSim split\n");
    return -1;
  }
  return nps;
}


/* ---------- conditional set builders (manual 10.1.2; survey
   2026-08-06) --------------------------------------------------------
   `Set NAME = (all,i,SRC: <cond>);` where <cond> is one of the three
   corpus shapes:
     COEF(i)               <op> <const>          [GDYN/gtapep SLUG]
     COEF(i,"ele")/COEF("ele",i)... <op> <const> [v7 ENDOWFLAG]
     sum{j,S2: MAP(j) = i, COEF2(j)} <op> <const> [GTAP-E/-EP/-AEZ
                                                   UNITD* flags]
   The condition is DATA-dependent, so it is evaluated here -- before
   set resolution -- straight from the input files (the condition
   coefficient and any mapping must be file-Read; formula-computed
   operands are a named fatal), and the statement is rewritten into an
   explicit element list plus the subset relation:
     set NAME (e1,...);  subset NAME is subset of SRC ;
   Everything downstream (sets_read, subset_map_build, superset_pos)
   then works untouched. Elements keep SRC order; an empty selection
   is a named fatal. Zero-cost when the TAB has no builder statements. */

static char *sb_iodata_path(cmf_file_entry *iodata, int nio, const char *logname) {
  int i;
  size_t j;
  for (i=0; i<nio; i++) {
    /* logical file names are case-insensitive (TAB text is lowercased
       by the preprocess; CMF entries keep the author's case) */
    for (j=0; ; j++) {
      if (tolower((int)iodata[i].logname[j])!=tolower((int)logname[j])) break;
      if (logname[j]=='\0') return iodata[i].filname;
    }
  }
  return NULL;
}

/* elements of a set declared with an explicit list or a
   read-elements statement; returns count or -1 */
#define SB_MAXELE 4096
/* nesting allowed while resolving a synthesized narrowing set: the IF
   rewrite emits at most one level over a declared set, the limit only
   stops a malformed tab looping */
#define SB_MAXDEPTH 8

static int sb_elements_d(char *tabfile, cmf_file_entry *iodata, int nio, const char *setname, char (*ele)[NAMESIZE], int depth);

static int sb_eqi(const char *a, const char *b) {
  for (; *a!='\0'&&*b!='\0'; a++,b++) if (tolower((int)*a)!=tolower((int)*b)) return 0;
  return *a=='\0'&&*b=='\0';
}

/* derived set expression for a set builder's source: named sets and
   quoted elements joined by + / union / - / intersect (no brackets) --
   anything else stays with the original forms */
static int sb_expr_form(const char *p) {
  while (*p==' ') p++;
  if (*p=='\"'||*p=='(') return 0;
  for (; *p!='\0'&&*p!=';'&&*p!='\n'&&*p!='#'; p++) if (*p=='('||*p=='\\') return 0;
  return 1;
}

static int sb_expr_elements(char *tabfile, cmf_file_entry *iodata, int nio, const char *p, char (*ele)[NAMESIZE], int depth) {
  char (*tmp)[NAMESIZE];
  int n=0,first=1,op='+';
  if (depth>=SB_MAXDEPTH) return -1;
  tmp=calloc(SB_MAXELE,NAMESIZE);
  if (tmp==NULL) return -1;
  for (;;) {
    char tok[NAMESIZE];
    int tl=0,nt,i,j;
    while (*p==' ') p++;
    if (*p=='\0'||*p==';'||*p=='\n'||*p=='\r'||*p=='#') break;
    if (*p=='+'||*p=='-'||*p=='&') { op=*p; p++; continue; }
    if (*p=='\"') {
      p++;
      while (*p!='\0'&&*p!='\"'&&tl<NAMESIZE-1) tok[tl++]=(char)tolower((int)*p++);
      tok[tl]='\0';
      if (*p=='\"') p++;
      strcpy(tmp[0],tok);
      nt=1;
    } else {
      while (*p!='\0'&&*p!=' '&&*p!='+'&&*p!='-'&&*p!='&'&&*p!=';'&&*p!='\n'&&*p!='\r'&&tl<NAMESIZE-1) tok[tl++]=*p++;
      tok[tl]='\0';
      if (sb_eqi(tok,"union")) { op='+'; continue; }
      if (sb_eqi(tok,"intersect")) { op='&'; continue; }
      nt=sb_elements_d(tabfile,iodata,nio,tok,tmp,depth+1);
      if (nt<0) { free(tmp); return -1; }
    }
    if (first) { for (i=0; i<nt&&n<SB_MAXELE; i++) strcpy(ele[n++],tmp[i]); first=0; continue; }
    if (op=='+') {
      for (i=0; i<nt; i++) {
        for (j=0; j<n; j++) if (sb_eqi(ele[j],tmp[i])) break;
        if (j==n&&n<SB_MAXELE) strcpy(ele[n++],tmp[i]);
      }
    } else {
      int m=0;
      for (j=0; j<n; j++) {
        int hit=0;
        for (i=0; i<nt; i++) if (sb_eqi(ele[j],tmp[i])) { hit=1; break; }
        if ((op=='&')==hit) { if (m!=j) strcpy(ele[m],ele[j]); m++; }
      }
      n=m;
    }
  }
  free(tmp);
  return first?-1:n;
}

static int sb_elements(char *tabfile, cmf_file_entry *iodata, int nio, const char *setname, char (*ele)[NAMESIZE]) {
  return sb_elements_d(tabfile,iodata,nio,setname,ele,0);
}

static int sb_elements_d(char *tabfile, cmf_file_entry *iodata, int nio, const char *setname, char (*ele)[NAMESIZE], int depth) {
  FILE *f;
  char line[TABREADLINE];
  int n=-1;
  size_t snlen=strlen(setname);
  f=teems_fopen(tabfile,"r");
  if (f==NULL) return -1;
  while (fgets(line,TABREADLINE,f)) {
    char *p=line;
    if (strncmp(p,"set",3)!=0) continue;
    p+=3;
    while (*p==' ') p++;
    if (strncmp(p,setname,snlen)!=0||(p[snlen]!=' '&&p[snlen]!='('&&p[snlen]!='#')) continue;
    p+=snlen;
    while (*p==' ') p++;
    if (*p=='#') { p=strchr(p+1,'#'); if (p==NULL) break; p++; while (*p==' ') p++; }
    if (strncmp(p,"(intertemporal)",15)==0) { p+=15; while (*p==' ') p++; }
    else if (strncmp(p,"(non_intertemporal)",19)==0) { p+=19; while (*p==' ') p++; }
    if (*p=='(') {
      /* explicit list, element ranges (c1 - c5) expanded as sets_read
         does (manual 11.2.2) */
      char lst[TABREADLINE],*e=strchr(p,')'),*q;
      if (e==NULL||(size_t)(e-p-1)>=sizeof(lst)) break;
      memcpy(lst,p+1,e-p-1);
      lst[e-p-1]='\0';
      if (elem_list_expand(lst,line,sizeof(line),setname)<0) break;
      n=0;
      for (q=line; *q!='\0'&&n<SB_MAXELE; ) {
        int tl=0;
        while (*q==' '||*q==',') q++;
        while (*q!='\0'&&*q!=','&&*q!=' '&&tl<NAMESIZE-1) ele[n][tl++]=*q++;
        ele[n][tl]='\0';
        if (tl>0) n++;
      }
      break;
    }
    if (strstr(p,"read")!=NULL) {
      char *hd=strstr(p,"header");
      char *fl=strstr(p,"file");
      char logname[NAMESIZE],header[NAMESIZE],*path;
      int tl=0;
      if (hd==NULL||fl==NULL) break;
      fl+=4;
      while (*fl==' ') fl++;
      while (*fl!='\0'&&*fl!=' '&&tl<NAMESIZE-1) logname[tl++]=*fl++;
      logname[tl]='\0';
      hd=strchr(hd,'\"');
      if (hd==NULL) break;
      hd++;
      tl=0;
      while (*hd!='\0'&&*hd!='\"'&&tl<NAMESIZE-1) header[tl++]=*hd++;
      header[tl]='\0';
      path=sb_iodata_path(iodata,nio,logname);
      if (path==NULL) break;
      {
        datafile_labels *lab=(datafile_labels *)calloc(SB_MAXELE,sizeof(datafile_labels));
        int i;
        if (lab==NULL) break;
        datafile_read_labels(header,path,SB_MAXELE,lab);
        n=0;
        for (i=0; i<SB_MAXELE&&lab[i].ch[0]!='\0'; i++) {
          int j;
          /* data-file labels enter the synthesized set in lowercase,
             the case sets_read gives every other element */
          for (j=0; j<NAMESIZE-1&&lab[i].ch[j]!='\0'; j++) ele[n][j]=(char)tolower((int)lab[i].ch[j]);
          ele[n][j]='\0';
          n++;
        }
        free(lab);
      }
      break;
    }
    if (*p=='='&&sb_expr_form(p+1)) {
      /* a derived source set, A + B, A union B, A - B, A intersect B,
         or a plain copy (evaluated left to right over the resolvable
         sets; vetting G-S6: REG = DST intersect ORG) */
      n=sb_expr_elements(tabfile,iodata,nio,p+1,ele,depth);
      break;
    }
    if (*p=='=') {
      /* a set the IF rewrite synthesized for an element-equality
         condition: "ele" & RANGE. Set builders are evaluated ahead of
         set resolution, so the one narrowing form teems-R emits is
         resolved here rather than left to the resolver. */
      char lit[NAMESIZE],rng[NAMESIZE];
      char (*rele)[NAMESIZE];
      int tl=0,nr,i;
      p++;
      while (*p==' ') p++;
      if (*p!='\"') break;
      p++;
      while (*p!='\0'&&*p!='\"'&&tl<NAMESIZE-1) lit[tl++]=*p++;
      lit[tl]='\0';
      if (*p!='\"') break;
      p++;
      while (*p==' ') p++;
      if (*p!='&') break;
      p++;
      while (*p==' ') p++;
      tl=0;
      while (*p!='\0'&&*p!=' '&&*p!=';'&&*p!='\n'&&*p!='\r'&&tl<NAMESIZE-1) rng[tl++]=*p++;
      rng[tl]='\0';
      if (rng[0]=='\0'||depth>=SB_MAXDEPTH) break;
      rele=calloc(SB_MAXELE,NAMESIZE);
      if (rele==NULL) break;
      nr=sb_elements_d(tabfile,iodata,nio,rng,rele,depth+1);
      n=(nr<0)?-1:0;
      /* case-insensitive: the literal comes from the lowercased TAB,
         the range from an explicit list (lowercased) or a data file
         (author's case); the range's spelling is kept so the caller's
         element lookups against that range match */
      for (i=0; i<nr; i++) {
        if (sb_eqi(rele[i],lit)) {
          strncpy(ele[0],rele[i],NAMESIZE-1);
          ele[0][NAMESIZE-1]='\0';
          n=1;
          break;
        }
      }
      free(rele);
      break;
    }
    break;
  }
  fclose(f);
  return n;
}

/* numeric text header: dims from the header line, values in the
   text-file layout the main reader (tab_parse.c Read) uses -- 2-D
   slices of dim-1 rows x dim-2 columns, slices over dims 3+ with dim 3
   fastest, blank lines between slices (the teems-R writer's layout;
   1-D/2-D degenerate to plain row-major). Values are returned
   re-ordered to C order (last dimension fastest) so the caller's
   stride arithmetic is layout-agnostic. Returns total count or -1 */
static int sb_read_reals(const char *path, const char *header, double **vals) {
  FILE *f;
  char line[DATREADLINE];
  int total=0,got=0;
  double *seq=NULL;
  *vals=NULL;
  f=teems_fopen((char *)path,"r");
  if (f==NULL) return -1;
  while (fgets(line,DATREADLINE,f)) {
    char *q=strchr(line,'\"');
    char hdr[NAMESIZE];
    int tl=0,dims[8],nd=0;
    if (q==NULL||strstr(line,"Header")==NULL) continue;
    q++;
    while (*q!='\0'&&*q!='\"'&&tl<NAMESIZE-1) hdr[tl++]=*q++;
    hdr[tl]='\0';
    if (str_cmp_ci(hdr,header)!=0) continue;
    {
      char *p=line;
      while (*p==' ') p++;
      while (*p>='0'&&*p<='9') {
        if (nd<8) dims[nd++]=atoi(p);
        while (*p>='0'&&*p<='9') p++;
        while (*p==' ') p++;
      }
    }
    if (nd==0) { fclose(f); return -1; }
    total=1;
    { int i; for (i=0; i<nd; i++) total*=dims[i]; }
    seq=(double *)malloc((size_t)total*sizeof(double));
    *vals=(double *)malloc((size_t)total*sizeof(double));
    if (seq==NULL||*vals==NULL) { free(seq); free(*vals); *vals=NULL; fclose(f); return -1; }
    while (got<total&&fgets(line,DATREADLINE,f)) {
      char *p=line;
      while (*p==' '||*p=='\t') p++;
      if (*p=='\n'||*p=='\r'||*p=='\0') continue; /* slice separator */
      if (strchr(line,'\"')!=NULL) break;          /* next header */
      while (*p!='\0'&&got<total) {
        while (*p==' '||*p==','||*p=='\t') p++;
        if (*p=='\0'||*p=='\n'||*p=='\r') break;
        seq[got++]=strtod(p,&p);
      }
    }
    fclose(f);
    if (got!=total) { free(seq); free(*vals); *vals=NULL; return -1; }
    {
      /* sequential position -> per-dimension index (antidim scheme of
         the main reader) -> C-order position */
      long antidim[8],cstride[8];
      int i,d;
      if (nd==1) { antidim[0]=1; }
      else {
        antidim[1]=1;
        antidim[0]=dims[1];
        for (d=2; d<nd; d++) antidim[d]=(d==2)?(long)dims[0]*dims[1]:antidim[d-1]*dims[d-1];
      }
      cstride[nd-1]=1;
      for (d=nd-2; d>=0; d--) cstride[d]=cstride[d+1]*dims[d+1];
      for (i=0; i<total; i++) {
        long l1=i,cpos=0;
        int idx[8];
        if (nd==1) idx[0]=(int)l1;
        else {
          for (d=nd-1; d>1; d--) { idx[d]=(int)(l1/antidim[d]); l1-=antidim[d]*idx[d]; }
          idx[0]=(int)(l1/antidim[0]); l1-=antidim[0]*idx[0];
          idx[1]=(int)(l1/antidim[1]);
        }
        for (d=0; d<nd; d++) cpos+=idx[d]*cstride[d];
        (*vals)[cpos]=seq[i];
      }
    }
    free(seq);
    return total;
  }
  fclose(f);
  return -1;
}

/* the coefficient's Read statement: logical file + header */
static int sb_coef_read_stmt(char *tabfile, const char *coef, char *logname, char *header) {
  FILE *f;
  char line[TABREADLINE];
  size_t cl=strlen(coef);
  int found=0;
  f=teems_fopen(tabfile,"r");
  if (f==NULL) return 0;
  while (fgets(line,TABREADLINE,f)) {
    char *p=line,*hd,*fl;
    int tl;
    if (strncmp(p,"read",4)!=0) continue;
    p+=4;
    while (*p==' ') p++;
    if (*p=='(') { p=strchr(p,')'); if (p==NULL) continue; p++; while (*p==' ') p++; }
    if (strncmp(p,coef,cl)!=0||(p[cl]!=' '&&p[cl]!='f')) continue;
    fl=strstr(p,"file");
    hd=strstr(p,"header");
    if (fl==NULL||hd==NULL) continue;
    fl+=4;
    while (*fl==' ') fl++;
    tl=0;
    while (*fl!='\0'&&*fl!=' '&&tl<NAMESIZE-1) logname[tl++]=*fl++;
    logname[tl]='\0';
    hd=strchr(hd,'\"');
    if (hd==NULL) continue;
    hd++;
    tl=0;
    while (*hd!='\0'&&*hd!='\"'&&tl<NAMESIZE-1) header[tl++]=*hd++;
    header[tl]='\0';
    found=1;
    break;
  }
  fclose(f);
  return found;
}

/* the coefficient declaration's per-dimension set names */
static int sb_coef_dims(char *tabfile, const char *coef, char dimset[][NAMESIZE]) {
  FILE *f;
  char line[TABREADLINE];
  int nd=-1;
  size_t cl=strlen(coef);
  f=teems_fopen(tabfile,"r");
  if (f==NULL) return -1;
  while (fgets(line,TABREADLINE,f)) {
    char *p=line,*nm;
    if (strncmp(p,"coefficient",11)!=0) continue;
    /* find " <coef>(" or " <coef> " after the quantifiers */
    nm=strstr(line,coef);
    while (nm!=NULL) {
      char before=(nm==line)?' ':nm[-1];
      char after=nm[cl];
      if ((before==' '||before==')')&&(after=='('||after==' '||after=='#'||after==';')) break;
      nm=strstr(nm+1,coef);
    }
    if (nm==NULL) continue;
    /* count and read the (all,idx,SET) quantifiers before the name */
    nd=0;
    p=line;
    while ((p=strstr(p,"(all,"))!=NULL&&p<nm) {
      char *c1=strchr(p+5,','),*c2;
      int tl=0;
      if (c1==NULL) { nd=-1; break; }
      c2=c1+1;
      while (*c2!='\0'&&*c2!=')'&&tl<NAMESIZE-1) dimset[nd][tl++]=*c2++;
      dimset[nd][tl]='\0';
      { /* trim spaces */
        char *s=dimset[nd],*d=dimset[nd];
        for (; *s!='\0'; s++) if (*s!=' ') *d++=*s;
        *d='\0';
      }
      nd++;
      if (nd>=MAXVARDIM) break;
      p=c2;
    }
    break;
  }
  fclose(f);
  return nd;
}

static int sb_ele_find(char (*ele)[NAMESIZE], int n, const char *name) {
  int i;
  for (i=0; i<n; i++) if (strcmp(ele[i],name)==0) return i;
  return -1;
}

static int sb_op_test(double v, const char *op, double c) {
  if (strcmp(op,"ne")==0||strcmp(op,"<>")==0) return v!=c;
  if (strcmp(op,"eq")==0||strcmp(op,"=")==0) return v==c;
  if (strcmp(op,"gt")==0||strcmp(op,">")==0) return v>c;
  if (strcmp(op,"lt")==0||strcmp(op,"<")==0) return v<c;
  if (strcmp(op,"ge")==0||strcmp(op,">=")==0) return v>=c;
  if (strcmp(op,"le")==0||strcmp(op,"<=")==0) return v<=c;
  return -1;
}

/* Indicator-formula operand (GTAP-AEZ UNITD* flags): a coefficient
   that is not file-Read but assigned only constants -- `C(i) = c`
   over a set and `C(i) = C(i) + [c]` / `- [c]` over a (sub)set, the
   shape the front end's IF rewrite emits for `C(i) = 0 + IF[i in S,
   1]`.  Applied in file order to the `nele` elements `ele` (elements
   of a formula set outside `ele` are skipped).  Returns 1 with *out
   allocated when every formula targeting `coef` has that shape and
   at least one exists, else 0 (the caller keeps its named fatal). */
static int sb_indicator_eval(char *tabfile, cmf_file_entry *iodata, int niodata, const char *coef, char (*ele)[NAMESIZE], int nele, double **out) {
  FILE *f;
  char line[TABREADLINE];
  char (*sele)[NAMESIZE]=NULL;
  double *v=NULL;
  int nf=0,ok=1;
  size_t cl=strlen(coef);
  f=teems_fopen(tabfile,"r");
  if (f==NULL) return 0;
  v=calloc(nele>0?nele:1,sizeof(double));
  sele=calloc(SB_MAXELE,NAMESIZE);
  if (v==NULL||sele==NULL) { free(v); free(sele); fclose(f); return 0; }
  while (ok&&fgets(line,TABREADLINE,f)) {
    char *p=line,*q,idx[NAMESIZE],set[NAMESIZE],rhs[TABREADLINE];
    int nq=0,tl,mode,k,ns,e;
    double c;
    if (strncmp(p,"formula",7)!=0) continue;
    p+=7;
    /* qualifier groups and the one quantifier */
    for (;;) {
      while (*p==' ') p++;
      if (*p!='(') break;
      if (strncmp(p,"(all,",5)==0) {
        p+=5;
        tl=0; while (*p!='\0'&&*p!=','&&tl<NAMESIZE-1) { if (*p!=' ') idx[tl++]=*p; p++; }
        idx[tl]='\0'; if (*p==',') p++;
        tl=0; while (*p!='\0'&&*p!=')'&&*p!=':'&&tl<NAMESIZE-1) { if (*p!=' ') set[tl++]=*p; p++; }
        set[tl]='\0';
        if (*p==':') nq=99; /* conditional quantifier: not an indicator shape */
        while (*p!='\0'&&*p!=')') p++;
        if (*p==')') p++;
        nq++;
      } else {
        while (*p!='\0'&&*p!=')') p++;
        if (*p==')') p++;
      }
    }
    /* the target: coef(idx) = */
    if (strncmp(p,coef,cl)!=0) continue;
    q=p+cl;
    while (*q==' ') q++;
    if (*q!='(') continue;
    nf++;
    if (nq!=1) { ok=0; break; }
    q++;
    tl=0; while (*q!='\0'&&*q!=')'&&tl<NAMESIZE-1) { if (*q!=' ') rhs[tl++]=*q; q++; }
    rhs[tl]='\0';
    if (strcmp(rhs,idx)!=0||*q!=')') { ok=0; break; }
    q++;
    while (*q==' ') q++;
    if (*q!='=') { ok=0; break; }
    q++;
    /* rhs without blanks, brackets or the terminator */
    tl=0;
    for (; *q!='\0'&&*q!=';'&&*q!='\n'&&*q!='\r'; q++) if (*q!=' '&&*q!='['&&*q!=']'&&*q!='('&&*q!=')'&&tl<TABREADLINE-1) rhs[tl++]=*q;
    rhs[tl]='\0';
    /* C idx + c  |  C idx - c  |  c   (the brackets are gone: "coef" "idx") */
    {
      char self[NAMESIZE*2];
      size_t sl;
      strcpy(self,coef); strcat(self,idx);
      sl=strlen(self);
      if (strncmp(rhs,self,sl)==0&&(rhs[sl]=='+'||rhs[sl]=='-')) {
        char *endp=NULL;
        mode=1;
        c=strtod(rhs+sl+1,&endp);
        if (endp==rhs+sl+1||*endp!='\0') { ok=0; break; }
        if (rhs[sl]=='-') c=-c;
      } else {
        char *endp=NULL;
        mode=0;
        c=strtod(rhs,&endp);
        if (endp==rhs||*endp!='\0') { ok=0; break; }
      }
    }
    ns=sb_elements(tabfile,iodata,niodata,set,sele);
    if (ns<0) { ok=0; break; }
    /* an empty quantifier set (a narrowing set whose element is not
       in the range) makes the formula a no-op, as TABLO has it */
    if (ns==0) continue;
    for (k=0; k<ns; k++) {
      e=sb_ele_find(ele,nele,sele[k]);
      if (e<0) continue;
      if (mode==0) v[e]=c; else v[e]+=c;
    }
  }
  fclose(f);
  free(sele);
  if (!ok||nf==0) { free(v); return 0; }
  *out=v;
  return 1;
}

/* one comparison of a set builder's condition: keep[k] for each source
   element. An index compared with a quoted element (r <> "usa", manual
   11.4.11) selects by name; the other forms are listed above. */
static int sb_leaf_keep(char *fname, cmf_file_entry *iodata, int niodata, const char *name, const char *idx, const char *src, char (*srcele)[NAMESIZE], int nsrc, char *cond, char *keep) {
  char op[8],*p,*q;
  int tl;
  double cval;
  memset(keep,0,SB_MAXELE);
  {
    char a[TABREADLINE],c[TABREADLINE],*o;
    int k=0,ci,cq=0,ol;
    for (ci=0; cond[ci]!='\0'&&k<TABREADLINE-1; ci++) { if (cond[ci]=='"') cq=!cq; if (cond[ci]!=' '||cq) a[k++]=cond[ci]; }
    a[k]='\0';
    o=strpbrk(a,"<>=");
    if (o!=NULL&&o>a) {
      ol=(o[0]=='<'&&o[1]=='>')?2:1;
      if ((o[0]=='='||ol==2)&&o[ol]!='\0') {
        int eq=(o[0]=='=');
        strcpy(c,o+ol);
        *o='\0';
        if ((sb_eqi(a,idx)&&c[0]=='"')||(sb_eqi(c,idx)&&a[0]=='"')) {
          char el[NAMESIZE],*e=(a[0]=='"')?a+1:c+1;
          int n=0;
          while (*e!='\0'&&*e!='"'&&n<NAMESIZE-1) el[n++]=*e++;
          el[n]='\0';
          for (k=0; k<nsrc; k++) keep[k]=(char)(sb_eqi(srcele[k],el)==eq);
          return 0;
        }
      }
    }
  }
  {
      /* split <operand> <op> <const>: find the comparison at depth 0,
         word ops need surrounding blanks stripped later */
      char opnd[TABREADLINE];
      int depth=0,i,oi=-1,olen=0;
      for (i=0; cond[i]!='\0'; i++) {
        char c=cond[i];
        if (c=='('||c=='{'||c=='[') depth++;
        else if (c==')'||c=='}'||c==']') depth--;
        else if (depth==0) {
          if (c=='<'||c=='>') { oi=i; olen=(cond[i+1]=='='||cond[i+1]=='>')?2:1; break; }
          if (c=='='&&(i==0||(cond[i-1]!='<'&&cond[i-1]!='>'))) { oi=i; olen=1; break; }
          if ((c==' ')&&((strncmp(cond+i+1,"ne ",3)==0)||(strncmp(cond+i+1,"eq ",3)==0)||
                         (strncmp(cond+i+1,"gt ",3)==0)||(strncmp(cond+i+1,"lt ",3)==0)||
                         (strncmp(cond+i+1,"ge ",3)==0)||(strncmp(cond+i+1,"le ",3)==0))) { oi=i+1; olen=2; break; }
        }
      }
      if (oi<0) {
        errmsg("Error: set builder %s: unsupported condition '%s' (supported: COEF(...) <op> const, or a mapping-conditional sum <op> const; manual 10.1.2)\n",name,cond);
        return -1;
      }
      strncpy(op,cond+oi,olen);
      op[olen]='\0';
      {
        char *endp=NULL;
        cval=strtod(cond+oi+olen,&endp);
        while (endp!=NULL&&*endp==' ') endp++;
        if (endp==cond+oi+olen||endp==NULL||*endp!='\0') {
          errmsg("Error: set builder %s: unsupported condition '%s' (a single comparison against a numeric constant; compound conditions are not supported; manual 10.1.2)\n",name,cond);
          return -1;
        }
      }
      strncpy(opnd,cond,oi);
      opnd[oi]='\0';
      { /* trim operand blanks at both ends */
        char *s=opnd;
        int e=(int)strlen(opnd);
        while (*s==' ') s++;
        while (e>0&&(opnd[e-1]==' ')) opnd[--e]='\0';
        memmove(opnd,s,strlen(s)+1);
      }
      {
        int k;
        if (str_find_ci(opnd,"$pos")==0) {
          /* $POS(i) or $POS(i,S) (manual 11.5.6): the position of each
             source element in the source set, or in S -- which must hold
             every element of the index's set (round 3: TERM $pos(c,COM)) */
          char a1[NAMESIZE],a2[NAMESIZE];
          char (*sele)[NAMESIZE]=NULL;
          int ns=0,ok=1;
          a1[0]=a2[0]='\0';
          p=opnd+4;
          while (*p==' '||*p=='('||*p=='{'||*p=='[') p++;
          tl=0; while (*p!='\0'&&*p!=','&&*p!=')'&&*p!='}'&&*p!=']'&&tl<NAMESIZE-1) { if (*p!=' ') a1[tl++]=*p; p++; }
          a1[tl]='\0';
          if (*p==',') {
            p++;
            tl=0; while (*p!='\0'&&*p!=')'&&*p!='}'&&*p!=']'&&tl<NAMESIZE-1) { if (*p!=' ') a2[tl++]=*p; p++; }
            a2[tl]='\0';
          }
          if (!sb_eqi(a1,idx)) {
            errmsg("Error: set builder %s: $POS(%s) must take the builder's index %s (manual 11.5.6)\n",name,a1,idx);
            ok=0;
          }
          if (ok&&a2[0]!='\0') {
            sele=calloc(SB_MAXELE,NAMESIZE);
            ns=(sele==NULL)?-1:sb_elements(fname,iodata,niodata,a2,sele);
            if (ns<=0) {
              errmsg("Error: set builder %s: cannot resolve the elements of set %s in $POS(%s,%s)\n",name,a2,a1,a2);
              ok=0;
            }
          }
          for (k=0; ok&&k<nsrc; k++) {
            double pos=k+1;
            if (sele!=NULL) {
              int e;
              for (e=0; e<ns; e++) if (sb_eqi(sele[e],srcele[k])) break;
              if (e==ns) {
                errmsg("Error: set builder %s: element %s of %s is not in %s; $POS(%s,%s) needs %s to range over a subset of %s (manual 11.5.6)\n",name,srcele[k],src,a2,a1,a2,a1,a2);
                ok=0;
                break;
              }
              pos=e+1;
            }
            keep[k]=(char)sb_op_test(pos,op,cval);
          }
          free(sele);
          if (!ok) return -1;
        }
        else if (strncmp(opnd,"sum",3)==0) {
          /* sum{j,S2: MAP(j) = idx, COEF2(j)} */
          char sset[NAMESIZE],mapname[NAMESIZE],c2[NAMESIZE];
          char logname[NAMESIZE],header[NAMESIZE],*path;
          char (*s2ele)[NAMESIZE]=NULL;
          datafile_labels *mlab=NULL;
          double *v2=NULL;
          int ns2,j,ok=1;
          p=opnd+3;
          while (*p=='('||*p=='{'||*p=='['||*p==' ') p++;
          tl=0; while (*p!='\0'&&*p!=','&&tl<NAMESIZE-1) { if (*p!=' ') tl++; p++; }
          if (*p==',') p++;
          tl=0; while (*p!='\0'&&*p!=':'&&tl<NAMESIZE-1) { if (*p!=' ') sset[tl++]=*p; p++; }
          sset[tl]='\0'; if (*p==':') p++;
          while (*p==' ') p++;
          tl=0; while (*p!='\0'&&*p!='('&&*p!='{'&&tl<NAMESIZE-1) { if (*p!=' ') mapname[tl++]=*p; p++; }
          mapname[tl]='\0';
          /* require MAP(sidx) = idx , COEF2(sidx) */
          q=strchr(p,',');
          if (q==NULL) ok=0;
          else {
            char inner[NAMESIZE*2];
            tl=0; q++;
            while (*q==' ') q++;
            while (*q!='\0'&&*q!='('&&*q!='{'&&tl<(int)sizeof(inner)-1) { if (*q!=' ') inner[tl++]=*q; q++; }
            inner[tl]='\0';
            strncpy(c2,inner,NAMESIZE-1);
            c2[NAMESIZE-1]='\0';
          }
          if (ok) {
            s2ele=calloc(SB_MAXELE,NAMESIZE);
            mlab=calloc(SB_MAXELE,sizeof(datafile_labels));
            ns2=(s2ele!=NULL&&mlab!=NULL)?sb_elements(fname,iodata,niodata,sset,s2ele):-1;
            if (ns2<=0) ok=0;
            else {
              /* the mapping's by_elements read */
              if (!sb_coef_read_stmt(fname,mapname,logname,header)) ok=0;
              else {
                path=sb_iodata_path(iodata,niodata,logname);
                if (path==NULL) ok=0;
                else datafile_read_labels(header,path,ns2,mlab);
              }
            }
            if (ok) {
              if (!sb_coef_read_stmt(fname,c2,logname,header)) {
                /* not file-Read: an indicator-formula operand (GTAP-AEZ) */
                if (!sb_indicator_eval(fname,iodata,niodata,c2,s2ele,ns2,&v2)) ok=0;
              }
              else {
                path=sb_iodata_path(iodata,niodata,logname);
                if (path==NULL||sb_read_reals(path,header,&v2)!=ns2) ok=0;
              }
            }
            if (ok) {
              for (k=0; k<nsrc; k++) {
                double acc=0;
                for (j=0; j<ns2; j++) if (strcmp(mlab[j].ch,srcele[k])==0) acc+=v2[j];
                keep[k]=(char)sb_op_test(acc,op,cval);
              }
            }
          }
          free(s2ele);
          free(mlab);
          free(v2);
          if (!ok) {
            errmsg("Error: set builder %s: cannot evaluate the mapping-conditional sum '%s' (the mapping must be file-Read and the summed coefficient file-Read or an indicator assigned only constants; manual 10.1.2)\n",name,opnd);
            return -1;
          }
        }
        else {
          /* COEF(args) with the loop index and optional quoted elements */
          char coef[NAMESIZE],args[MAXVARDIM][NAMESIZE];
          char dimset[MAXVARDIM][NAMESIZE];
          char logname[NAMESIZE],header[NAMESIZE],*path;
          double *cv=NULL;
          int nargs=0,nd,k,d,ok=1,loopdim=-1;
          long stride[MAXVARDIM],fixoff;
          int dsz[MAXVARDIM];
          p=opnd;
          tl=0;
          while (*p!='\0'&&*p!='('&&*p!='{'&&tl<NAMESIZE-1) { if (*p!=' ') coef[tl++]=*p; p++; }
          coef[tl]='\0';
          if (*p!='\0') p++;
          while (*p!='\0'&&*p!=')'&&*p!='}'&&nargs<MAXVARDIM) {
            tl=0;
            while (*p==' ') p++;
            if (*p=='\"') { p++; while (*p!='\0'&&*p!='\"'&&tl<NAMESIZE-1) args[nargs][tl++]=*p++; if (*p=='\"') p++; }
            else while (*p!='\0'&&*p!=','&&*p!=')'&&*p!='}'&&tl<NAMESIZE-1) { if (*p!=' ') args[nargs][tl++]=*p; p++; }
            args[nargs][tl]='\0';
            nargs++;
            if (*p==',') p++;
          }
          nd=sb_coef_dims(fname,coef,dimset);
          if (nd!=nargs||nd<=0) ok=0;
          if (ok&&!sb_coef_read_stmt(fname,coef,logname,header)) {
            /* a 1-D indicator assigned only constants evaluates from
               its formulas (GTAP-AEZ); anything else stays fatal */
            if (nd==1&&strcmp(args[0],idx)==0&&sb_indicator_eval(fname,iodata,niodata,coef,srcele,nsrc,&cv)) {
              for (k=0; k<nsrc; k++) keep[k]=(char)sb_op_test(cv[k],op,cval);
              free(cv);
              cv=NULL;
              ok=2;
            } else {
              errmsg("Error: set builder %s: condition coefficient %s must be Read from an input file or be an indicator assigned only constants (formula-computed operands cannot drive set resolution; manual 10.1.2)\n",name,coef);
              return -1;
            }
          }
          if (ok==1) {
            /* per-dim sizes + fixed-element offsets */
            char (*dele)[NAMESIZE]=calloc(SB_MAXELE,NAMESIZE);
            char (*ldele)[NAMESIZE]=calloc(SB_MAXELE,NAMESIZE);
            int nde,*lpos=NULL;
            fixoff=0;
            if (dele==NULL||ldele==NULL) ok=0;
            for (d=0; ok&&d<nd; d++) {
              nde=sb_elements(fname,iodata,niodata,dimset[d],dele);
              if (nde<=0) { ok=0; break; }
              dsz[d]=nde;
              if (strcmp(args[d],idx)==0) { loopdim=d; memcpy(ldele,dele,(size_t)nde*NAMESIZE); }
              else {
                int e=sb_ele_find(dele,nde,args[d]);
                if (e<0) { ok=0; break; }
                stride[d]=e; /* park the element position; strides resolved below */
              }
            }
            free(dele);
            if (ok&&loopdim>=0) {
              long st=1;
              for (d=nd-1; d>=0; d--) { long tmpst=st; st*=dsz[d]; stride[d]=(d==loopdim)?tmpst:stride[d]*tmpst; }
              for (d=0; d<nd; d++) if (d!=loopdim) fixoff+=stride[d];
            }
            else ok=0;
            if (ok) {
              path=sb_iodata_path(iodata,niodata,logname);
              if (path==NULL) ok=0;
              else {
                long total=1;
                for (d=0; d<nd; d++) total*=dsz[d];
                if (sb_read_reals(path,header,&cv)!=total) ok=0;
              }
            }
            /* the source set may be the dimension set or a subset of it
               (a derived source, REG = DST intersect ORG; vetting G-S6):
               each source element reads its position in the dimension */
            if (ok) {
              lpos=malloc((nsrc>0?nsrc:1)*sizeof(int));
              if (lpos==NULL) ok=0;
              for (k=0; ok&&k<nsrc; k++) {
                lpos[k]=(dsz[loopdim]==nsrc)?k:sb_ele_find(ldele,dsz[loopdim],srcele[k]);
                if (lpos[k]<0) {
                  errmsg("Error: set builder %s: element %s of source set %s is not in %s's dimension set %s\n",name,srcele[k],src,coef,dimset[loopdim]);
                  ok=0;
                }
              }
            }
            if (ok) for (k=0; k<nsrc; k++) keep[k]=(char)sb_op_test(cv[fixoff+(long)lpos[k]*stride[loopdim]],op,cval);
            free(lpos);
            free(ldele);
            free(cv);
          }
          if (!ok) {
            errmsg("Error: set builder %s: cannot evaluate condition '%s' (declaration/read/dimension resolution failed; manual 10.1.2)\n",name,cond);
            return -1;
          }
        }
      }
  }
  return 0;
}

int tab_setbuilder_transform(char *fname, cmf_file_entry *iodata, int niodata) {
  FILE *f,*fout;
  char line[TABREADLINE],tmpname[TABREADLINE];
  int any=0,rc=0;
  f=teems_fopen(fname,"r");
  if (f==NULL) { errmsg("Error: cannot open %s\n",fname); return -1; }
  while (fgets(line,TABREADLINE,f)) {
    if (strncmp(line,"set",3)==0&&strstr(line,"(all,")!=NULL&&strchr(line,':')!=NULL&&strchr(line,'=')!=NULL) any=1;
  }
  fclose(f);
  if (!any) return 0;
  f=teems_fopen(fname,"r");
  strcpy(tmpname,fname);
  strcat(tmpname,"_sb");
  fout=f==NULL?NULL:teems_fopen(tmpname,"w");
  if (f==NULL||fout==NULL) {
    if (f!=NULL) fclose(f);
    errmsg("Error: cannot open set-builder scratch file\n");
    return -1;
  }
  while (rc==0&&fgets(line,TABREADLINE,f)) {
    char name[NAMESIZE],idx[NAMESIZE],src[NAMESIZE],cond[TABREADLINE];
    char *p,*q;
    int tl;
    if (!(strncmp(line,"set",3)==0&&(p=strstr(line,"="))!=NULL&&(q=strstr(line,"(all,"))!=NULL&&q>p&&strchr(q,':')!=NULL)) {
      fputs(line,fout);
      continue;
    }
    /* set NAME [# label #] = (all,idx,SRC: cond) ; */
    p=line+3;
    while (*p==' ') p++;
    tl=0;
    while (*p!='\0'&&*p!=' '&&*p!='='&&*p!='#'&&tl<NAMESIZE-1) name[tl++]=*p++;
    name[tl]='\0';
    q=strstr(line,"(all,")+5;
    tl=0;
    while (*q!='\0'&&*q!=','&&tl<NAMESIZE-1) { if (*q!=' ') idx[tl++]=*q; q++; }
    idx[tl]='\0';
    if (*q==',') q++;
    tl=0;
    while (*q!='\0'&&*q!=':'&&tl<NAMESIZE-1) { if (*q!=' ') src[tl++]=*q; q++; }
    src[tl]='\0';
    if (*q!=':') { errmsg("Error: malformed set-builder statement: %s",line); rc=-1; break; }
    q++;
    /* condition text up to the builder's closing ')' (depth-aware) */
    {
      int depth=0;
      tl=0;
      for (; *q!='\0'; q++) {
        if (*q=='('||*q=='{'||*q=='[') depth++;
        else if (*q==')'||*q=='}'||*q==']') { if (depth==0) break; depth--; }
        if (tl<TABREADLINE-1) cond[tl++]=*q;
      }
      cond[tl]='\0';
    }
    cond_unwrap(cond);
    {
      /* the condition's boolean structure (manual 11.4.5): one keep
         vector per comparison, combined by AND/OR/NOT */
      char (*srcele)[NAMESIZE]=calloc(SB_MAXELE,NAMESIZE);
      char *keep=calloc(SB_MAXELE,1),*lk=NULL,cs[TABREADLINE];
      int nsrc,k,nkept=0,ci,cq=0,nleaf=0,nprog=0,lb[ICOND_MAXLEAF],le[ICOND_MAXLEAF];
      signed char prog[ICOND_MAXPROG];
      if (srcele==NULL||keep==NULL) { free(srcele); free(keep); rc=-1; break; }
      nsrc=sb_elements(fname,iodata,niodata,src,srcele);
      if (nsrc<=0) {
        errmsg("Error: set builder %s: cannot resolve the elements of source set %s (explicit list or read-elements declarations only)\n",name,src);
        free(srcele);
        free(keep);
        rc=-1;
        break;
      }
      for (ci=0,k=0; cond[ci]!='\0'&&k<TABREADLINE-1; ci++) { if (cond[ci]=='"') cq=!cq; if (cond[ci]!=' '||cq) cs[k++]=cond[ci]; }
      cs[k]='\0';
      if (strchr(cs,COND_AND)==NULL&&strchr(cs,COND_OR)==NULL&&strchr(cs,COND_NOT)==NULL) rc=sb_leaf_keep(fname,iodata,niodata,name,idx,src,srcele,nsrc,cond,keep);
      else if (cond_tree_rpn(cs,&nleaf,lb,le,prog,&nprog,name)<0||(lk=calloc((size_t)nleaf*SB_MAXELE,1))==NULL) rc=-1;
      else {
        int lf,st[ICOND_MAXPROG];
        for (lf=0; rc==0&&lf<nleaf; lf++) {
          char leaf[TABREADLINE];
          memcpy(leaf,cs+lb[lf],le[lf]-lb[lf]);
          leaf[le[lf]-lb[lf]]='\0';
          rc=sb_leaf_keep(fname,iodata,niodata,name,idx,src,srcele,nsrc,leaf,lk+(size_t)lf*SB_MAXELE);
        }
        for (k=0; rc==0&&k<nsrc; k++) {
          int n=0,j;
          st[0]=0;
          for (j=0; j<nprog; j++) {
            if (prog[j]>=0) st[n++]=lk[(size_t)prog[j]*SB_MAXELE+k];
            else if (prog[j]==-3) st[n-1]=!st[n-1];
            else { n--; st[n-1]=(prog[j]==-1)?(st[n-1]&&st[n]):(st[n-1]||st[n]); }
          }
          keep[k]=(char)st[0];
        }
      }
      free(lk);
      if (rc!=0) { free(srcele); free(keep); break; }
        for (k=0; k<nsrc; k++) if (keep[k]) nkept++;
        /* an empty selection is a legal empty set (manual 11.7.9):
           statements over it have no tuples, sums over it are zero */
        fprintf(fout,"set %s (",name);
        {
          int first=1,k2;
          for (k2=0; k2<nsrc; k2++) if (keep[k2]) { fprintf(fout,"%s%s",first?"":",",srcele[k2]); first=0; }
        }
        fprintf(fout,");\n");
        fprintf(fout,"subset %s is subset of %s ;\n",name,src);
        printf("set builder %s: %d of %d elements of %s selected\n",name,nkept,nsrc,src);
        free(srcele);
        free(keep);
    }
  }
  fclose(f);
  fclose(fout);
  if (rc==0) {
    if (rename(tmpname,fname)!=0) { errmsg("Error: cannot rename %s\n",tmpname); rc=-1; }
  }
  else remove(tmpname);
  return rc;
}
