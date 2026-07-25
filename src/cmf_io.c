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

int cmf_count_files(char *fname,char *comsyntax) {
  FILE * filehandle;
  char line[TABREADLINE]="\0";
  int j=0;
  filehandle = fopen(fname,"r");
  if(filehandle==NULL){
    printf("Error: cannot open %s\n",fname);
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
    printf("Error: cannot open %s\n",filename);
    return -1;
  }

  while (fgets(line,TABREADLINE,filehandle)) {
    strcpy(linecopy,line);
    readitem = strtok(line,"\"");
    readitem = strtok(NULL,"\"");
    if (readitem != NULL) {
      nlength1=0;
      while (readitem[nlength1] != '\0'&&readitem[nlength1] != ' ') nlength1++;
      if (nlength1==nlength&&strncmp(readitem,varname,nlength) == 0) {
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
  if(succ==0)printf("Error: header \"%s\" not found in %s\n",varname,filename);
  return succ;
}

int datafile_read_labels(char *varname, char *filename,dim_t d1, datafile_labels *record) {
  FILE * dfile;
  char line[DATREADLINE],header[NAMESIZE],varnamecpy[NAMESIZE];
  dim_t nlength=0,nhead=0,reccount = 0, count1=0,i;
  char *readitem=NULL;
  strcpy(varnamecpy,varname);
  while (str_replace_all(varnamecpy," ", ""));
  while (varnamecpy[nlength] != '\0') nlength++;
  dfile = fopen(filename,"r");

  while (fgets(line,DATREADLINE,dfile)) {
    readitem = strtok(line,"\"");
    readitem = strtok(NULL,"\"");
    if (readitem != NULL) {
      strcpy(header,readitem);
      while (str_replace_all(header," ", ""));
      nhead=0;
      while (header[nhead] != '\0') nhead++;
      if(nhead<nlength)nhead=nlength;
      if (strncmp(readitem,varnamecpy,nhead) == 0) {
        while (fgets(line,DATREADLINE,dfile)) {
          while (str_replace_all(line,"\r", ""));
          while (str_replace_all(line,"  ", " "));
          if (line[0]=='\n') count1=1;
          if (count1!=1) {
            readitem = strtok(line,"\n");
            if (readitem != NULL) {
              record[reccount].dim1=reccount;
              strcpy(record[reccount].ch,readitem);
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
    printf("Error: cannot open %s\n",filename);
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
    k=str_find_ci(readitem,syntax);
    readitem=readitem+k;
    k1=str_find_ci(readitem,",");
    for(j=0;j<k1+1;j++)finditem[j]=readitem[j];
    finditem[j]='\0';
    readitem=readitem+k1;
    k2=str_find_ci(readitem,finditem);
    while(k2!=-1) {
      line1=readitem+k2;
      sprintf(temp, "%d", l);
      strcat(temp,"?");
      l++;
      for(j=4;j<k1;j++){
        replitem[j-4]=finditem[j];
        newreplitem[j-4]=finditem[j];
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
      while(str_replace_all(line1,replitem1,newreplitem1));

      replitem1[0]='(';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,"]");
      newreplitem1[0]='(';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,"]");
      while(str_replace_all(line1,replitem1,newreplitem1));

      replitem1[0]='(';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,"}");
      newreplitem1[0]='(';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,"}");
      while(str_replace_all(line1,replitem1,newreplitem1));

      replitem1[0]='(';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,",");
      newreplitem1[0]='(';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,",");
      while(str_replace_all(line1,replitem1,newreplitem1));

      replitem1[0]='(';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1," ");
      newreplitem1[0]='(';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1," ");
      while(str_replace_all(line1,replitem1,newreplitem1));


      replitem1[0]='[';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,")");
      newreplitem1[0]='[';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,")");
      while(str_replace_all(line1,replitem1,newreplitem1));

      replitem1[0]='[';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,"]");
      newreplitem1[0]='[';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,"]");
      while(str_replace_all(line1,replitem1,newreplitem1));

      replitem1[0]='[';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,"}");
      newreplitem1[0]='[';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,"}");
      while(str_replace_all(line1,replitem1,newreplitem1));

      replitem1[0]='[';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,",");
      newreplitem1[0]='[';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,",");
      while(str_replace_all(line1,replitem1,newreplitem1));

      replitem1[0]='[';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1," ");
      newreplitem1[0]='[';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1," ");
      while(str_replace_all(line1,replitem1,newreplitem1));

      replitem1[0]='{';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,")");
      newreplitem1[0]='{';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,")");
      while(str_replace_all(line1,replitem1,newreplitem1));

      replitem1[0]='{';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,"]");
      newreplitem1[0]='{';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,"]");
      while(str_replace_all(line1,replitem1,newreplitem1));

      replitem1[0]='{';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,"}");
      newreplitem1[0]='{';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,"}");
      while(str_replace_all(line1,replitem1,newreplitem1));

      replitem1[0]='{';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,",");
      newreplitem1[0]='{';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,",");
      while(str_replace_all(line1,replitem1,newreplitem1));

      replitem1[0]='{';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1," ");
      newreplitem1[0]='{';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1," ");
      while(str_replace_all(line1,replitem1,newreplitem1));

      replitem1[0]=' ';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,")");
      newreplitem1[0]=' ';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,")");
      while(str_replace_all(line1,replitem1,newreplitem1));

      replitem1[0]=' ';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,"]");
      newreplitem1[0]=' ';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,"]");
      while(str_replace_all(line1,replitem1,newreplitem1));

      replitem1[0]=' ';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,"}");
      newreplitem1[0]=' ';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,"}");
      while(str_replace_all(line1,replitem1,newreplitem1));

      replitem1[0]=' ';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,",");
      newreplitem1[0]=' ';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,",");
      while(str_replace_all(line1,replitem1,newreplitem1));

      replitem1[0]=' ';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1," ");
      newreplitem1[0]=' ';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1," ");
      while(str_replace_all(line1,replitem1,newreplitem1));

      replitem1[0]=',';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,")");
      newreplitem1[0]=',';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,")");
      while(str_replace_all(line1,replitem1,newreplitem1));

      replitem1[0]=',';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,"]");
      newreplitem1[0]=',';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,"]");
      while(str_replace_all(line1,replitem1,newreplitem1));

      replitem1[0]=',';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,"}");
      newreplitem1[0]=',';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,"}");
      while(str_replace_all(line1,replitem1,newreplitem1));

      replitem1[0]=',';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1,",");
      newreplitem1[0]=',';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1,",");
      while(str_replace_all(line1,replitem1,newreplitem1));

      replitem1[0]=',';
      replitem1[1]='\0';
      strcat(replitem1,replitem);
      strcat(replitem1," ");
      newreplitem1[0]=',';
      newreplitem1[1]='\0';
      strcat(newreplitem1,newreplitem);
      strcat(newreplitem1," ");
      while(str_replace_all(line1,replitem1,newreplitem1));

      str_replace_all(formulain,temp,line1);
      strcpy(line,formulain);
      k2=str_find_ci(readitem,finditem);
    }
  }
  return 0;
}


int tab_preprocess(char *filename, char *newtabfile) {
  FILE * filehandle,*fout;
  char line[TABREADLINE]="\0",line1[TABREADLINE],line2[TABREADLINE],indx[NAMESIZE],indx1[NAMESIZE],indx2[NAMESIZE],*readitem,*readitem1,commsyntax[NAMESIZE],readline[TABREADLINE],readline1[TABREADLINE],*n,newtabfile1[TABREADLINE];
  char setname[NAMESIZE],newset[NAMESIZE],varname[NAMESIZE],*n1,setelement[TABREADLINE];//,*ne,*np;//,*n2;
  char msetele[MAXVARDIM][NAMESIZE],msetsymb[MAXVARDIM][NAMESIZE],mset[MAXVARDIM][NAMESIZE];
  char assertmsg[TABREADLINE],*am1,*am2;
  filehandle = fopen(filename,"r");
  if(filehandle==NULL){
    printf("Error: cannot open %s\n",filename);
    return -1;
  }
  int check,i1,i2,i,setindx,varindx,l1,l2,l3,l4,k1,k2,j,j1,j2;//,necheck,npcheck;//,j;,check1
  strcpy(newtabfile1,newtabfile);
  str_replace_all(newtabfile1,".","1.");
  fout = fopen(newtabfile1,"w");
  readline[0]='\0';
  commsyntax[0]='\0';
  while (fgets(line,TABLINESIZE,filehandle)) {
    str_replace_all(line,"\v"," ");
    str_replace_all(line,"![[!","\v");
    str_replace_all(line,"!]]!","\v");
    if (strlen(readline)+strlen(line)>=sizeof(readline)) {
      printf("Error: TAB statement too long (exceeds %d chars)\n",TABREADLINE);
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
      if(str_find_ci(readline,"subset ")==1||str_find_ci(readline,"subset ")==0) {
        strcpy(commsyntax,"subset");
        check=1;
      }
      if(str_find_ci(readline,"file ")==1||str_find_ci(readline,"file ")==0) {
        strcpy(commsyntax,"file");
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
      if (str_find_ci(readline,"formula ")==1||str_find_ci(readline,"formula ")==0) {
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
    if (l3==0&&str_count_ci(line,"\"")<3) l3=-1;
    if (l1==0||l2==0||l3==0||l4==0) {
      if (l1==0||l4==0) {
        n=strstr(line,"c_");
        while (n!=NULL) {
          if (line[n-line-1]==' '||line[n-line-1]=='+'||line[n-line-1]=='-'||line[n-line-1]=='*'||line[n-line-1]=='/'||line[n-line-1]=='^'||line[n-line-1]==','||line[n-line-1]=='=') {
            line[n-line]='p';
            n=strstr(line,"c_");
          } else n=strstr(n+1,"c_");
        }
      }
      strcpy(line1,line);
      n=strchr(line1,'\"');
      if (n!=NULL) {
        readline[0]='\0';
        readline1[0]='\0';
        while (n!=NULL) {
          strncpy(line2,line1,n-line1);
          line2[n-line1]='\0';
          n1=strrchr(line2,'(');
          if (n1==NULL) {
            printf("Error: malformed indexed expression in TAB file: %s\n",line);
            return -1;
          }
          setindx=str_count_ci(n1,",");
          varindx=n1-line2+1;
          line2[varindx]='\0';
          n1=strrchr(line2,' ');
          if (n1!=NULL) varindx=n1-line2;
          n1=strrchr(line2,',');
          if (n1!=NULL&&(n1-line2)>varindx) varindx=n1-line2;
          n1=strrchr(line2,'*');
          if (n1!=NULL&&(n1-line2)>varindx) varindx=n1-line2;
          n1=strrchr(line2,'/');
          if (n1!=NULL&&(n1-line2)>varindx) varindx=n1-line2;
          n1=strrchr(line2,'+');
          if (n1!=NULL&&(n1-line2)>varindx) varindx=n1-line2;
          n1=strrchr(line2,'-');
          if (n1!=NULL&&(n1-line2)>varindx) varindx=n1-line2;
          n1=strrchr(line2,'^');
          if (n1!=NULL&&(n1-line2)>varindx) varindx=n1-line2;
          n1=strrchr(line2,'=');
          if (n1!=NULL&&(n1-line2)>varindx) varindx=n1-line2;
          if (cmf_strcpy_bounded(varname,&line2[varindx+1],sizeof(varname))) {
            printf("Error: TAB statement too complex in tab_preprocess: %s\n",line);
            return -1;
          }
          n1=strrchr(varname,'(');
          if (n1==NULL) {
            printf("Error: malformed indexed expression in TAB file: %s\n",line);
            return -1;
          }
          varname[n1-varname+1]='\0';
           j1=0;
           j2=0;
           while (varname[j1]!= '\0'){
             if(varname[j1]!= '(')j2++;
             j1++;
           }
           n1=varname;
           if(j2>1){
             j2=0;
             while (j1>-1){
             if(varname[j1]== '(')j2++;
             if(j2==2){
               n1=varname+j1+1;
               break;
             }
             j1--;
             }
           }
          tab_read_set_name(newtabfile1,n1,setindx,setname);
          strcpy(indx,"i");
          sprintf(indx1, "%d",i);
          strcat(indx,indx1);
          strcpy(newset,"sub_");
          strcat(newset,indx);
          i1=n-line1;
          readitem=strtok(n,"\"");
          if (readitem==NULL||cmf_strcpy_bounded(setelement,readitem,sizeof(setelement))) {
            printf("Error: malformed indexed expression in TAB file: %s\n",line);
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
            printf("Error: TAB statement too complex in tab_preprocess: %s\n",line);
            return -1;
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
            printf("Error: TAB statement too complex in tab_preprocess: %s\n",line);
            return -1;
          }
          str_replace_all(line,indx1,indx);
          readitem1=strtok(line," ");
          if (readitem1==NULL||
              cmf_strcpy_bounded(line2,readitem1,sizeof(line2)) ||
              cmf_strcat_bounded(line2," ",sizeof(line2))) {
            printf("Error: malformed indexed expression in TAB file: %s\n",line);
            return -1;
          }
          readitem1=strtok(NULL," ");
          if (readitem1==NULL) {
            printf("Error: malformed indexed expression in TAB file: %s\n",line);
            return -1;
          }
          strcpy(indx2,indx);
          strcat(indx2,",");
          if(strstr(readitem1,indx2)==NULL) {
            if(l1==0){//if(l1==0&&strchr(line2,'(')!=NULL) {
              if (cmf_strcat_bounded(line2,readitem1,sizeof(line2)) ||
                  cmf_strcat_bounded(line2," ",sizeof(line2))) {
                printf("Error: TAB statement too complex in tab_preprocess: %s\n",line);
                return -1;
              }
              readitem1=strtok(NULL," ");
              if (readitem1==NULL) {
                printf("Error: malformed indexed expression in TAB file: %s\n",line);
                return -1;
              }
            }
            if (cmf_strcat_bounded(line2," (all,",sizeof(line2)) ||
                cmf_strcat_bounded(line2,indx,sizeof(line2)) ||
                cmf_strcat_bounded(line2,",",sizeof(line2)) ||
                cmf_strcat_bounded(line2,newset,sizeof(line2)) ||
                cmf_strcat_bounded(line2,")",sizeof(line2)) ||
                cmf_strcat_bounded(line2,readitem1,sizeof(line2)) ||
                cmf_strcat_bounded(line2," ",sizeof(line2))) {
              printf("Error: TAB statement too complex in tab_preprocess: %s\n",line);
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
              printf("Error: TAB statement too complex in tab_preprocess: %s\n",line);
              return -1;
            }
          }
          readitem1=strtok(NULL,"\n");
          if(readitem1!=NULL&&cmf_strcat_bounded(line2,readitem1,sizeof(line2))) {
            printf("Error: TAB statement too complex in tab_preprocess: %s\n",line);
            return -1;
          }
          if (cmf_strcat_bounded(line2,"\n",sizeof(line2))) {
            printf("Error: TAB statement too complex in tab_preprocess: %s\n",line);
            return -1;
          }
          strcpy(line,line2);
          strcpy(line1,line);
          n=strchr(line1,'\"');
          if (l3==0&&str_count_ci(line1,"\"")<3) n=NULL;
          i++;
        }
        fprintf(fout,"%s%s%s",readline,readline1,line2);
      } else fprintf(fout,"%s",line);
    } else fprintf(fout,"%s",line);
  }
  fclose(filehandle);
  fclose(fout);
  remove(newtabfile1);
  return 1;
}

int outputs_write_csv(char *filename, char *newdatlogname, char *newdatfile,set_def *sets,dim_t nset, set_element *set_elems,array_def *coefs,offset_t ncof,offset_t ncofele,array_def *vars,offset_t nvar,offset_t nvarele, elem_value *elem_vals) {
  FILE * filehandle,*fout;
  char line[TABREADLINE]="\0",*readline,comsyntax[TABREADLINE],longname[TABREADLINE],datline[DATREADLINE],varname[NAMESIZE],*vname1,header[NAMESIZE],setsize[DATREADLINE],tempname[NAMESIZE];
  filehandle = fopen(filename,"r");
  long int i,n,j,j1,j2,innerloop,outerloop,l,indx;
  long int setindx[MAXVARDIM],antidim[MAXVARDIM];
  strcpy(comsyntax,"to file ");
  strcat(comsyntax,newdatlogname);
  strcat(comsyntax," header \"");
  n=strlen(comsyntax);
  fout = fopen(newdatfile,"w");
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
            if(coefs[i].size==1) {
              innerloop=sets[coefs[i].setid[0]].size;
              outerloop=0;
            }
            if(coefs[i].size==2) {
              innerloop=sets[coefs[i].setid[0]].size*sets[coefs[i].setid[1]].size;
              outerloop=1;
            }
            if(coefs[i].size>2) {
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
    if(str_find_ci(line,"equation ")>-1||str_find_ci(line,"update ")>-1) {
      linelght=strlen(line);
      for (i=0; i<nvar; i++) {
        p=strchr(line,';');
        if(p==NULL) break;
        line[p-line+1]='\n';
        line[p-line+2]='\0';
        n=str_count_ci(line,vars[i].cofname);
        lvar=strlen(vars[i].cofname);
        l=0;
        for (j=0; j<n; j++) {
          l1=str_find_ci(&line[l],vars[i].cofname);
          if(l1<0) break;
          l=l+l1;
          if(strncmp(vars[i].cofname,"p_",2)!=0&&vars[i].level_par==false) if(line[l+lvar]==' '||line[l+lvar]=='('||line[l+lvar]=='+'||line[l+lvar]=='-'||line[l+lvar]=='*'||line[l+lvar]=='/'||line[l+lvar]=='^'||line[l+lvar]==']'||line[l+lvar]==','||line[l+lvar]==';'||line[l+lvar]=='=')if(l==0||line[l-1]==' '||line[l-1]=='+'||line[l-1]=='-'||line[l-1]=='*'||line[l-1]=='/'||line[l-1]=='^'||line[l-1]=='['||line[l-1]=='('||line[l-1]==','||line[l-1]=='=') {
                memmove(&line[l+2],&line[l],linelght-l);
                line[l]='p';
                line[l+1]='_';
                l=l+2;
                linelght+=2;
              }
          l=l+strlen(vars[i].cofname);
        }
      }

    }
    fprintf(fout,"%s",line);
  }
  fclose(filehandle);
  fclose(fout);
  return 1;
}


int tab_read_set_name(char *filename, char *varname, int indx, char *setname) {
  FILE * filehandle;
  int n,i;
  offset_t lsize;
  lsize=TABREADLINE+1;
  char commsyntax[NAMESIZE],varname1[NAMESIZE+2],indxname[NAMESIZE],line[TABREADLINE+1],line1[TABREADLINE+1],*p,tmp[TABREADLINE+1];//,varname2[NAMESIZE+2],varname3[NAMESIZE+2]
  strcpy(varname1,")");
  strcat(varname1,varname);
  filehandle = fopen(filename,"r");
  if (filehandle==NULL) return -1;
  strcpy(commsyntax,"variable");
  while (tab_next_statement(commsyntax,filehandle,line,lsize)) {
    while (str_replace_all(line," ", ""));
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
      if (p==NULL||strlen(p)+1>=sizeof(tmp)) { fclose(filehandle); return -1; }
      strcpy(tmp,p);
      strcat(tmp,",");
      p=NULL;
      for (i=0; i<indx+1; i++) {
        if(i==0) p=strtok(tmp,",");
        else p=strtok(NULL,",");
        if (p==NULL) { fclose(filehandle); return -1; }
      }
      if (strlen(p)+1>=sizeof(indxname)) { fclose(filehandle); return -1; }
      strcpy(indxname,",");
      strcat(indxname,p);
      strcpy(line,line1);
      n=str_find_ci(line,indxname);
      if (n<0) { fclose(filehandle); return -1; }
      p=strtok(line+n,",");
      p=strtok(NULL,")");
      if (p==NULL||strlen(p)>=NAMESIZE) { fclose(filehandle); return -1; }
      strcpy(setname,p);
      return 1;
    }
  }
  fclose(filehandle);
  filehandle = fopen(filename,"r");
  if (filehandle==NULL) return -1;
  strcpy(commsyntax,"coefficient");
  while (tab_next_statement(commsyntax,filehandle,line,lsize)) {
    while (str_replace_all(line," ", ""));
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
      if (p==NULL||strlen(p)+1>=sizeof(tmp)) { fclose(filehandle); return -1; }
      strcpy(tmp,p);
      strcat(tmp,",");
      p=NULL;
      for (i=0; i<indx+1; i++) {
        if(i==0) p=strtok(tmp,",");
        else p=strtok(NULL,",");
        if (p==NULL) { fclose(filehandle); return -1; }
      }
      if (strlen(p)+2>=sizeof(indxname)) { fclose(filehandle); return -1; }
      strcpy(indxname,",");
      strcat(indxname,p);
      strcat(indxname,",");
      strcpy(line,line1);
      n=str_find_ci(line,indxname);
      if (n<0) { fclose(filehandle); return -1; }
      p=strtok(line+n,",");
      p=strtok(NULL,")");
      if (p==NULL||strlen(p)>=NAMESIZE) { fclose(filehandle); return -1; }
      strcpy(setname,p);
      return 1;
    }
  }
  fclose(filehandle);
  return -1;
}


/* CMF "Assertions = yes|no|warn ;" (GEMPACK manual 25.3).
   Returns 0/1/2 (no/warn/yes); default yes. */
int cmf_assertions_mode(char *filename) {
  FILE *f;
  char l[TABREADLINE];
  char *p;
  int k,mode=2;
  f=fopen(filename,"r");
  if(f==NULL)return 2;
  while(fgets(l,TABREADLINE,f)!=NULL) {
    k=str_find_ci(l,"assertions");
    if(k<0)continue;
    p=strchr(l+k,'=');
    if(p==NULL)continue;
    p++;
    while(*p==' '||*p=='\t')p++;
    if(str_find_ci(p,"no")==0)mode=0;
    else if(str_find_ci(p,"warn")==0)mode=1;
    else if(str_find_ci(p,"yes")==0)mode=2;
  }
  fclose(f);
  return mode;
}

/* CMF "range test initial values = yes|no|warn ;" and
   "range test updated values = updated|extrapolated|both|no|warn ;"
   (GEMPACK manual 25.4.4). Modes 0/1/2 = no/warn/fatal; the manual's
   default outside automatic accuracy is WARN for both legs. The
   updated-leg selectors updated/extrapolated/both all map to fatal at
   our update-time checks (extrapolation is not separately
   instrumented). */
void cmf_range_test_modes(char *filename, int *ini, int *upd) {
  FILE *f;
  char l[TABREADLINE];
  char *p;
  int k;
  *ini=1;
  *upd=1;
  f=fopen(filename,"r");
  if(f==NULL)return;
  while(fgets(l,TABREADLINE,f)!=NULL) {
    k=str_find_ci(l,"range test initial values");
    if(k>=0) {
      p=strchr(l+k,'=');
      if(p==NULL)continue;
      p++;
      while(*p==' '||*p=='\t')p++;
      if(str_find_ci(p,"no")==0)*ini=0;
      else if(str_find_ci(p,"warn")==0)*ini=1;
      else if(str_find_ci(p,"yes")==0)*ini=2;
      continue;
    }
    k=str_find_ci(l,"range test updated values");
    if(k<0)continue;
    p=strchr(l+k,'=');
    if(p==NULL)continue;
    p++;
    while(*p==' '||*p=='\t')p++;
    if(str_find_ci(p,"no")==0)*upd=0;
    else if(str_find_ci(p,"warn")==0)*upd=1;
    else if(str_find_ci(p,"updated")==0||str_find_ci(p,"extrapolated")==0||str_find_ci(p,"both")==0)*upd=2;
  }
  fclose(f);
}

/* Validate every Default statement up front (manual 10.19; audit A6),
   one scan of the preprocessed TAB (one statement per line,
   lowercased). Supported -- applied positionally by the readers:
   Coefficient parameter/non_parameter, Variable linear/levels/change/
   percent_change, Formula initial/always. Accepted no-ops: Equation
   linear/not_add_homotopy (the solver's only equation kind). Fatal:
   Equation levels/add_homotopy (unsupported semantics), Coefficient
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
        printf("Error: Coefficient (default=%s): bound defaults are not supported (single bound slot, audit A9)\n",val);
      else printf("Error: unknown Coefficient default '%s'\n",val);
    } else if(strncmp(line,"variable",8)==0) {
      if(strcmp(val,"linear")==0||strcmp(val,"levels")==0||strcmp(val,"change")==0||strcmp(val,"percent_change")==0)continue;
      printf("Error: unknown Variable default '%s'\n",val);
    } else if(strncmp(line,"formula",7)==0) {
      if(strcmp(val,"initial")==0||strcmp(val,"always")==0)continue;
      printf("Error: unknown Formula default '%s'\n",val);
    } else if(strncmp(line,"equation",8)==0) {
      if(strcmp(val,"linear")==0||strcmp(val,"not_add_homotopy")==0)continue;
      if(strcmp(val,"levels")==0)
        printf("Error: Equation (default=levels) is not supported -- the solver handles linearized equations only\n");
      else if(strncmp(val,"add_homotopy",12)==0)
        printf("Error: Equation (default=add_homotopy) is not supported\n");
      else printf("Error: unknown Equation default '%s'\n",val);
    } else {
      printf("Error: Default statements apply only to Coefficient/Variable/Formula/Equation: %s",line);
    }
    bad=1;
  }
  fclose(f);
  return bad?-1:0;
}

/* Split the preprocessed TAB around POSTSIM (BEGIN)/(END) sections
   (manual 10.18, 12.2.1). Declarations (set/subset/coefficient/file)
   stay in the ordinary file -- a single namespace; what separates
   PostSim is EXECUTION order -- while executables (formula/assertion/
   zerodivide) move to the _ps companion consumed once after the solve.
   Write/Display in sections are dropped (outputs ride the write-all
   dump); PostSim Read and the forbidden statements are fatal.
   Returns the PostSim executable count (0 = no sections), -1 on
   error. */
int tab_postsim_split(char *newtabfile, char *psfile) {
  FILE *fin,*fmain,*fps;
  char line[TABREADLINE],tmpname[TABREADLINE];
  int inps=0,nps=0,found=0;
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
  rewind(fin);
  strcpy(tmpname,newtabfile);
  strcat(tmpname,"_o");
  fmain=fopen(tmpname,"w");
  fps=fopen(psfile,"w");
  if(fmain==NULL||fps==NULL) {
    printf("Error: cannot open PostSim split scratch files\n");
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
      fputs(line,fmain);
      continue;
    }
    if(strncmp(line,"variable ",9)==0||strncmp(line,"equation ",9)==0||strncmp(line,"update ",7)==0||strncmp(line,"transfer ",9)==0||strncmp(line,"omit ",5)==0||strncmp(line,"substitute ",11)==0||strncmp(line,"backsolve ",10)==0||strncmp(line,"complementarity",15)==0||strstr(line,"(default")!=NULL) {
      printf("Error: statement not allowed in a PostSim section (manual 12.2.1): %s",line);
      fclose(fin);
      fclose(fmain);
      fclose(fps);
      return -1;
    }
    if(strncmp(line,"read ",5)==0) {
      printf("Error: PostSim Read is not supported yet; read the data in the ordinary part into a (parameter) coefficient instead\n");
      fclose(fin);
      fclose(fmain);
      fclose(fps);
      return -1;
    }
    if(strncmp(line,"set ",4)==0||strncmp(line,"subset ",7)==0||strncmp(line,"coefficient ",12)==0||strncmp(line,"file ",5)==0) {
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
    printf("Error: unrecognized statement in a PostSim section: %s",line);
    fclose(fin);
    fclose(fmain);
    fclose(fps);
    return -1;
  }
  fclose(fin);
  fclose(fmain);
  fclose(fps);
  if(rename(tmpname,newtabfile)!=0) {
    printf("Error: cannot finalize the PostSim split\n");
    return -1;
  }
  return nps;
}

/* CMF "PostSim = yes|no ;" run-time switch (manual ch.12); default yes */
int cmf_postsim_on(char *filename) {
  FILE *f;
  char l[TABREADLINE];
  char *p;
  int k,on=1;
  f=fopen(filename,"r");
  if(f==NULL)return 1;
  while(fgets(l,TABREADLINE,f)!=NULL) {
    k=str_find_ci(l,"postsim");
    if(k<0)continue;
    if(str_find_ci(l,"xpostsim")>-1)continue;
    p=strchr(l+k,'=');
    if(p==NULL)continue;
    p++;
    while(*p==' '||*p=='\t')p++;
    if(str_find_ci(p,"no")==0)on=0;
    else if(str_find_ci(p,"yes")==0)on=1;
  }
  fclose(f);
  return on;
}
