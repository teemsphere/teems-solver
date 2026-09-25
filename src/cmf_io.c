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
  filehandle = fopen(fname,"r");
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
  filehandle = fopen(filename,"r");
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
  dfile = fopen(filename,"r");
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

int cmf_read(char *filename, int niodata, cmf_file_entry *iodata, char *tabfile, char *closure, char *shock) {
  FILE * filehandle;
  char line[TABREADLINE],*readitem,commsyntax[NAMESIZE];
  int j=0,k;
  strcpy(commsyntax,"iodata");
  filehandle = fopen(filename,"r");
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
  filehandle = fopen(filename,"r");
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
  filehandle = fopen(filename,"r");
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

  filehandle = fopen(filename,"r");
  strcpy(commsyntax,"tabfile");
  while (tab_next_statement(commsyntax,filehandle,line,TABREADLINE)) {
    readitem = strtok(line,"\"");
    readitem = strtok(NULL,"\"");
    strcpy(tabfile,readitem);
  }
  fclose(filehandle);
  filehandle = fopen(filename,"r");
  strcpy(commsyntax,"closure");
  while (tab_next_statement(commsyntax,filehandle,line,TABREADLINE)) {
    readitem = strtok(line,"\"");
    readitem = strtok(NULL,"\"");
    strcpy(closure,readitem);
  }
  fclose(filehandle);
  filehandle = fopen(filename,"r");
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
int tab_preprocess(char *filename, char *newtabfile) {
  int r=tab_preprocess_run(filename,newtabfile);
  tab_decl_index_free();
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

#define MAXLITDEDUP 64
static int tab_preprocess_run(char *filename, char *newtabfile) {
  FILE * filehandle,*fout;
  char dedup_set[MAXLITDEDUP][NAMESIZE],dedup_ele[MAXLITDEDUP][NAMESIZE],dedup_idx[MAXLITDEDUP][NAMESIZE];
  int ndedup=0;
  char line[TABREADLINE]="\0",line1[TABREADLINE],line2[TABREADLINE],indx[NAMESIZE],indx1[NAMESIZE],indx2[NAMESIZE],*readitem,*readitem1,commsyntax[NAMESIZE],readline[TABREADLINE],readline1[TABREADLINE],*n,newtabfile1[TABREADLINE];
  char setname[NAMESIZE],newset[NAMESIZE],varname[NAMESIZE],*n1,setelement[TABREADLINE];//,*ne,*np;//,*n2;
  char assertmsg[TABREADLINE],*am1,*am2;
  char rawline[TABLINESIZE],*rawpos;
  filehandle = fopen(filename,"r");
  if(filehandle==NULL){
    errmsg("Error: cannot open %s\n",filename);
    return -1;
  }
  int check,i1,i2,i,setindx,l1,l2,l3,l4,k1,k2,l5;//,necheck,npcheck;//,j;,check1
  strcpy(newtabfile1,newtabfile);
  str_replace_all(newtabfile1,".","1.");
  fout = fopen(newtabfile1,"w");
  readline[0]='\0';
  commsyntax[0]='\0';
  /* GEMPACK allows several statements on one physical line: each raw
     line is fed to the accumulator piecewise, cut after every ';' that
     lies outside a ! comment, a ![[! !]]! block comment and a # label #
     (parity taken over the pending statement + the piece so far) and
     is followed by further text, so every iteration below sees at most
     one statement terminator */
  rawpos=NULL;
  while (rawpos!=NULL||fgets(rawline,TABLINESIZE,filehandle)) {
    if (rawpos==NULL) {
      str_replace_all(rawline,"\v"," ");
      str_replace_all(rawline,"![[!","\v");
      str_replace_all(rawline,"!]]!","\v");
      rawpos=rawline;
    }
    {
      int pb=str_count_char(readline,'!')%2,pv=str_count_char(readline,'\v')%2,ph=str_count_char(readline,'#')%2;
      char *c,*cut=NULL;
      for (c=rawpos; *c!='\0'; c++) {
        if (*c=='!') pb^=1;
        else if (*c=='\v') pv^=1;
        else if (*c=='#') ph^=1;
        else if (*c==';'&&!pb&&!pv&&!ph) {
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
      errmsg("Error: TAB statement too long (exceeds %d chars)\n",TABREADLINE);
      return -1;
    }
    strcat(readline,line);
    n=strrchr(line,';');
    i1=str_count_char(readline,'!');
    while (i1>1) {
      while (str_strip_comment(readline,"!"));
      i1-=2;
    }
    i2=str_count_char(readline,'\v');
    while (i2>1) {
      while (str_strip_comment(readline,"\v"));
      i2-=2;
    }
    if (n!=NULL&&i1==0&&i2==0) {
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
            if (*am1=='\n'||*am1=='\r'||*am1=='\v') *am1=' ';
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
        errmsg("Error: TAB statement too long (exceeds %d chars)\n",TABREADLINE);
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
      if(str_find_ci(readline,"subset ")==1||str_find_ci(readline,"subset ")==0) {
        strcpy(commsyntax,"subset");
        check=1;
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
          errmsg("Error: TAB statement too long (exceeds %d chars)\n",TABREADLINE);
          fclose(filehandle);
          fclose(fout);
          return -1;
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
  fclose(filehandle);
  fclose(fout);
  filehandle = fopen(newtabfile1,"r");
  if (filehandle==NULL) return -1;
  fout = fopen(newtabfile,"w");
  i=0;
  while (fgets(line,TABREADLINE,filehandle)) {
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
      if (l1==0||l4==0) {
        n=strstr(line,"c_");
        while (n!=NULL) {
          if (line[n-line-1]==' '||line[n-line-1]=='+'||line[n-line-1]=='-'||line[n-line-1]=='*'||line[n-line-1]=='/'||line[n-line-1]=='^'||line[n-line-1]==','||line[n-line-1]=='=') {
            line[n-line]='p';
            n=strstr(line,"c_");
          } else n=strstr(n+1,"c_");
        }
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
          if (n1==NULL&&l2==0&&strstr(line,"(by_elements)")!=NULL) {
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
            if (oe-os+2>(int)sizeof(varname)) {
              errmsg("Error: TAB statement too complex in tab_preprocess: %s\n",line);
              return -1;
            }
            memcpy(varname,&line2[os],oe-os+1);
            varname[oe-os+1]='\0';
          }
          n1=varname;
          if (tab_read_set_name(newtabfile1,n1,setindx,setname)==-1) {
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
  filehandle = fopen(filename,"r");
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
    errmsg("Error: cannot open %s for writing\n",newdatfile);
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
        readline=line+i+10;
        i=str_find_ci(readline,"\"");
        strncpy(longname,readline,i);
        longname[i]='\0';
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
        i=str_find_ci(line," ");
        readline=line+i+1;
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
        readline=line+i+10;
        i=str_find_ci(readline,"\"");
        strncpy(longname,readline,i);
        longname[i]='\0';
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
      }
    }
  }
  fclose(filehandle);
  fclose(fout);
  return 1;
}

int tab_write_variables(char *filename, char *newtabfile,array_def *vars,offset_t nvar) {
  FILE * filehandle,*fout;
  char line[TABREADLINE+1]="\0",*p;//,nvarname[nvar][NAMESIZE+2],*p;//,line1[DATREADLINE];//,*ne,*np;//,*n2;
  filehandle = fopen(filename,"r");
  offset_t i,n,j,l,l1,linelght;
  int lvar;
  if (filehandle==NULL) return -1;
  fout = fopen(newtabfile,"w");
  while (fgets(line,TABREADLINE,filehandle)) {
    int eqpos=str_find_ci(line,"equation ");
    int updpos=str_find_ci(line,"update ");
    /* update statements have no condition machinery: a ':' used to
       make the set lookup miss and expand over sets[0] in silence;
       the per-step update pass may never run (a broken system dies in
       the solver first), so fail here (M3) */
    if((updpos==0||updpos==1)&&strchr(line,':')!=NULL) {
      errmsg("Error: conditions in Update statements are not supported\n");
      fclose(filehandle);
      fclose(fout);
      MPI_Abort(PETSC_COMM_WORLD,1);
      return -1;
    }
    if(eqpos>-1||updpos>-1) {
      linelght=strlen(line);
      for (i=0; i<nvar; i++) {
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
          if(strncmp(vars[i].cofname,"p_",2)!=0&&vars[i].level_par==false) if(line[l+lvar]==')'||line[l+lvar]=='}'||line[l+lvar]==' '||line[l+lvar]=='('||line[l+lvar]=='+'||line[l+lvar]=='-'||line[l+lvar]=='*'||line[l+lvar]=='/'||line[l+lvar]=='^'||line[l+lvar]==']'||line[l+lvar]==','||line[l+lvar]==';'||line[l+lvar]=='=')if(l==0||line[l-1]==' '||line[l-1]=='+'||line[l-1]=='-'||line[l-1]=='*'||line[l-1]=='/'||line[l-1]=='^'||line[l-1]=='['||line[l-1]=='('||line[l-1]==')'||line[l-1]==','||line[l-1]=='=') { /* ')' : a quantifier written flush against the body, "(all,c,com)p(c)" (2026-09-24) */
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
    fprintf(fout,"%s",line);
  }
  fclose(filehandle);
  fclose(fout);
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
  FILE *filehandle=fopen(filename,"r");
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

/* one pass of the lookup over an indexed statement kind; varname1 is
   shared between the passes on purpose (the p_/c_ fallback rewrites
   it once and the coefficient pass then sees the rewritten name, as
   the file-scanning original did). Returns 1 found, -1 malformed,
   0 not in this kind. */
static int tab_read_set_name_pass(char **stmts, int nstmt, char *varname, char *varname1, int indx, char *setname) {
  int n,i,k;
  char indxname[NAMESIZE],line[TABREADLINE+1],line1[TABREADLINE+1],*p,tmp[TABREADLINE+1];
  for (k=0; k<nstmt; k++) {
    strcpy(line,stmts[k]);
    strcpy(line1,line);
    n=str_find_ci(line,varname1);
    if (n==-1&&((varname[0]=='p'&&varname[1]=='_')||(varname[0]=='c'&&varname[1]=='_'))) {
      strcpy(varname1,")");
      strcat(varname1,varname+2);
      n=str_find_ci(line,varname1);
    }
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
  strcpy(varname1,")");
  strcat(varname1,varname);
  if (decl_index_build(filename)<0) return -1;
  r=tab_read_set_name_pass(decl_var,n_decl_var,varname,varname1,indx,setname);
  if (r!=0) return r;
  r=tab_read_set_name_pass(decl_cof,n_decl_cof,varname,varname1,indx,setname);
  if (r!=0) return r;
  return -1;
}




/* Validate every Default statement up front (manual 10.19; audit A6),
   one scan of the preprocessed TAB (one statement per line,
   lowercased). Supported -- applied positionally by the readers:
   Coefficient parameter/non_parameter, Variable linear/levels/change/
   percent_change, Formula initial/always, Equation linear/levels (the
   levels transform qualifies unqualified equations positionally).
   Accepted no-op: Equation not_add_homotopy. Fatal: Equation
   add_homotopy (unsupported semantics), Coefficient
   lower_bound/upper_bound defaults (single bound slot, audit A9), and
   any unknown keyword or value. Returns 0 ok, -1 fatal. */
int tab_defaults_validate(char *fname) {
  FILE *f;
  char line[TABREADLINE],val[NAMESIZE];
  int bad=0;
  f=fopen(fname,"r");
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
      if(strncmp(val,"add_homotopy",12)==0)
        errmsg("Error: Equation (default=add_homotopy) is not supported\n");
      else errmsg("Error: unknown Equation default '%s'\n",val);
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
  fin=fopen(newtabfile,"r");
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
  fmain=fopen(tmpname,"w");
  fps=fopen(psfile,"w");
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

static int sb_elements(char *tabfile, cmf_file_entry *iodata, int nio, const char *setname, char (*ele)[NAMESIZE]) {
  return sb_elements_d(tabfile,iodata,nio,setname,ele,0);
}

static int sb_elements_d(char *tabfile, cmf_file_entry *iodata, int nio, const char *setname, char (*ele)[NAMESIZE], int depth) {
  FILE *f;
  char line[TABREADLINE];
  int n=-1;
  size_t snlen=strlen(setname);
  f=fopen(tabfile,"r");
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
    if (*p=='(') {
      /* explicit list */
      n=0;
      p++;
      while (*p!='\0'&&*p!=')') {
        int tl=0;
        while (*p==' '||*p==',') p++;
        while (*p!='\0'&&*p!=','&&*p!=')'&&*p!=' '&&tl<NAMESIZE-1) ele[n][tl++]=*p++;
        ele[n][tl]='\0';
        if (tl>0) n++;
        if (n>=SB_MAXELE) break;
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
  f=fopen((char *)path,"r");
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
  f=fopen(tabfile,"r");
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
  f=fopen(tabfile,"r");
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
  f=fopen(tabfile,"r");
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

int tab_setbuilder_transform(char *fname, cmf_file_entry *iodata, int niodata) {
  FILE *f,*fout;
  char line[TABREADLINE],tmpname[TABREADLINE];
  int any=0,rc=0;
  f=fopen(fname,"r");
  if (f==NULL) { errmsg("Error: cannot open %s\n",fname); return -1; }
  while (fgets(line,TABREADLINE,f)) {
    if (strncmp(line,"set",3)==0&&strstr(line,"(all,")!=NULL&&strchr(line,':')!=NULL&&strchr(line,'=')!=NULL) any=1;
  }
  fclose(f);
  if (!any) return 0;
  f=fopen(fname,"r");
  strcpy(tmpname,fname);
  strcat(tmpname,"_sb");
  fout=f==NULL?NULL:fopen(tmpname,"w");
  if (f==NULL||fout==NULL) {
    if (f!=NULL) fclose(f);
    errmsg("Error: cannot open set-builder scratch file\n");
    return -1;
  }
  while (rc==0&&fgets(line,TABREADLINE,f)) {
    char name[NAMESIZE],idx[NAMESIZE],src[NAMESIZE],op[8],cond[TABREADLINE];
    char *p,*q;
    int tl;
    double cval;
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
        rc=-1;
        break;
      }
      strncpy(op,cond+oi,olen);
      op[olen]='\0';
      {
        char *endp=NULL;
        cval=strtod(cond+oi+olen,&endp);
        while (endp!=NULL&&*endp==' ') endp++;
        if (endp==cond+oi+olen||endp==NULL||*endp!='\0') {
          errmsg("Error: set builder %s: unsupported condition '%s' (a single comparison against a numeric constant; compound conditions are not supported; manual 10.1.2)\n",name,cond);
          rc=-1;
          break;
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
      /* SRC elements */
      {
        char (*srcele)[NAMESIZE]=calloc(SB_MAXELE,NAMESIZE);
        char keep[SB_MAXELE];
        int nsrc,k,nkept=0;
        if (srcele==NULL) { rc=-1; break; }
        nsrc=sb_elements(fname,iodata,niodata,src,srcele);
        if (nsrc<=0) {
          errmsg("Error: set builder %s: cannot resolve the elements of source set %s (explicit list or read-elements declarations only)\n",name,src);
          free(srcele);
          rc=-1;
          break;
        }
        memset(keep,0,sizeof(keep));
        if (strncmp(opnd,"sum",3)==0) {
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
            free(srcele);
            rc=-1;
            break;
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
              free(srcele);
              rc=-1;
              break;
            }
          }
          if (ok==1) {
            /* per-dim sizes + fixed-element offsets */
            char (*dele)[NAMESIZE]=calloc(SB_MAXELE,NAMESIZE);
            int nde;
            fixoff=0;
            if (dele==NULL) ok=0;
            for (d=0; ok&&d<nd; d++) {
              nde=sb_elements(fname,iodata,niodata,dimset[d],dele);
              if (nde<=0) { ok=0; break; }
              dsz[d]=nde;
              if (strcmp(args[d],idx)==0) loopdim=d;
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
            if (ok) {
              if (dsz[loopdim]!=nsrc) {
                errmsg("Error: set builder %s: the loop index %s must range over %s's dimension set %s exactly\n",name,idx,coef,dimset[loopdim]);
                ok=0;
              }
            }
            if (ok) for (k=0; k<nsrc; k++) keep[k]=(char)sb_op_test(cv[fixoff+(long)k*stride[loopdim]],op,cval);
            free(cv);
          }
          if (!ok&&rc==0) {
            errmsg("Error: set builder %s: cannot evaluate condition '%s' (declaration/read/dimension resolution failed; manual 10.1.2)\n",name,cond);
            free(srcele);
            rc=-1;
            break;
          }
        }
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
      }
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
