#include <teems_solver.h>

/* bounded in-place replace-all used by the declaration parsers' set-symbol
   substitution: see the definition below str_replace_all. */
static int str_subst_all_bounded(char *line, const char *finditem, const char *replitem, size_t linesz);

int formula_normalize(char *fomulain) {
  int index,i,i1,i2,j;
  char fpart1[TABREADLINE],*p=NULL,*p1=NULL;
  while (str_replace_all(fomulain, " ", ""));
  p=strchr(fomulain,')');
  if (p==NULL) {
    while (str_replace_char(fomulain, '[', '('));
    while (str_replace_char(fomulain, ']', ')'));
    return 0;
  }
  strncpy(fpart1, fomulain, p-fomulain);
  fpart1[p-fomulain] = '\0';
  p1=strrchr(fpart1,'(');
  if (p1==NULL) {
    printf("Error: unbalanced parentheses in formula: %s\n",fomulain);
    return -1;
  }
  index=p1-fpart1;//ha_cgerevfind(fpart1,"(");
  fpart1[index]='\0';
  for(i=index-1; i>-1; i--) {
    if(fpart1[i]=='+'||fpart1[i]=='-'||fpart1[i]=='*'||fpart1[i]=='/'||fpart1[i]=='['||fpart1[i]=='('||fpart1[i]==','||fpart1[i]=='^') {
      break;
    }
  }
  if (i>-1) {
    i=index-i;
    switch (i) {
    case 1:
      fomulain[p-fomulain]=']';
      fomulain[index]='[';
      break;
    case 3:
      j=index-i;
      i=0;
      if(fpart1[j+1]=='i') {
        i++;
      }
      if(fpart1[j+2]=='f') {
        i++;
      }
      if (i==2) {
        fomulain[p-fomulain]=']';
        fomulain[index]='[';
      } else {
        fomulain[p-fomulain]='}';
        fomulain[index]='{';
      }
      break;
    case 4:
      j=index-i;
      i=0;
      if(fpart1[j+1]=='s') {
        i++;
      }
      if(fpart1[j+2]=='u') {
        i++;
      }
      if(fpart1[j+3]=='m') {
        i++;
      }
      i1=0;
      if(fpart1[j+1]=='a') {
        i1++;
      }
      if(fpart1[j+2]=='b') {
        i1++;
      }
      if(fpart1[j+3]=='s') {
        i1++;
      }
      /* exp( -- intrinsic (parity plan 3.1) */
      i2=0;
      if(fpart1[j+1]=='e') {
        i2++;
      }
      if(fpart1[j+2]=='x') {
        i2++;
      }
      if(fpart1[j+3]=='p') {
        i2++;
      }
      /* max( / min( -- multi-arg intrinsics (parity plan 3.1) */
      if (i==3||i1==3||i2==3||strncmp(fpart1+j+1,"max",3)==0||strncmp(fpart1+j+1,"min",3)==0) {
        fomulain[p-fomulain]=']';
        fomulain[index]='[';
      } else {
        fomulain[p-fomulain]='}';
        fomulain[index]='{';
      }
      break;
    case 5:
      j=index-i;
      i=0;
      if(fpart1[j+1]=='i') {
        i++;
      }
      if(fpart1[j+2]=='d') {
        i++;
      }
      if(fpart1[j+3]=='0') {
        i++;
      }
      if(fpart1[j+4]=='1') {
        i++;
      }
      i1=0;
      if(fpart1[j+1]=='l') {
        i1++;
      }
      if(fpart1[j+2]=='o') {
        i1++;
      }
      if(fpart1[j+3]=='g') {
        i1++;
      }
      if(fpart1[j+4]=='e') {
        i1++;
      }
      /* sqrt( -- intrinsic (parity plan 3.1) */
      i2=0;
      if(fpart1[j+1]=='s') {
        i2++;
      }
      if(fpart1[j+2]=='q') {
        i2++;
      }
      if(fpart1[j+3]=='r') {
        i2++;
      }
      if(fpart1[j+4]=='t') {
        i2++;
      }
      /* id0v( -- multi-arg intrinsic (parity plan 3.1) */
      if (i==4||i1==4||i2==4||strncmp(fpart1+j+1,"id0v",4)==0) {
        fomulain[p-fomulain]=']';
        fomulain[index]='[';
      } else {
        fomulain[p-fomulain]='}';
        fomulain[index]='{';
      }
      break;
    case 6:
      /* log10( / round( -- intrinsics (parity plan 3.1) */
      j=index-i;
      i=0;
      if(fpart1[j+1]=='l'&&fpart1[j+2]=='o'&&fpart1[j+3]=='g'&&fpart1[j+4]=='1'&&fpart1[j+5]=='0')i=5;
      i1=0;
      if(fpart1[j+1]=='r'&&fpart1[j+2]=='o'&&fpart1[j+3]=='u'&&fpart1[j+4]=='n'&&fpart1[j+5]=='d')i1=5;
      if (i==5||i1==5) {
        fomulain[p-fomulain]=']';
        fomulain[index]='[';
      } else {
        fomulain[p-fomulain]='}';
        fomulain[index]='{';
      }
      break;
    case 7:
      /* trunc0( / truncb( -- intrinsics (parity plan 3.1) */
      j=index-i;
      if (strncmp(fpart1+j+1,"trunc0",6)==0||strncmp(fpart1+j+1,"truncb",6)==0) {
        fomulain[p-fomulain]=']';
        fomulain[index]='[';
      } else {
        fomulain[p-fomulain]='}';
        fomulain[index]='{';
      }
      break;
    default:
      fomulain[p-fomulain]='}';
      fomulain[index]='{';
      break;
    }
  } else {
    i=index-i;
    switch (i) {
    case 1:
      fomulain[p-fomulain]=']';
      fomulain[index]='[';
      break;
    case 3:
      i=0;
      if(fpart1[0]=='i') {
        i++;
      }
      if(fpart1[1]=='f') {
        i++;
      }
      if (i==2) {
        fomulain[p-fomulain]=']';
        fomulain[index]='[';
      } else {
        fomulain[p-fomulain]='}';
        fomulain[index]='{';
      }
      break;
    case 4:
      i=0;
      if(fpart1[0]=='s') {
        i++;
      }
      if(fpart1[1]=='u') {
        i++;
      }
      if(fpart1[2]=='m') {
        i++;
      }
      i1=0;
      if(fpart1[0]=='a') {
        i1++;
      }
      if(fpart1[1]=='b') {
        i1++;
      }
      if(fpart1[2]=='s') {
        i1++;
      }
      /* exp( -- intrinsic (parity plan 3.1) */
      i2=0;
      if(fpart1[0]=='e') {
        i2++;
      }
      if(fpart1[1]=='x') {
        i2++;
      }
      if(fpart1[2]=='p') {
        i2++;
      }
      /* max( / min( -- multi-arg intrinsics (parity plan 3.1) */
      if (i==3||i1==3||i2==3||strncmp(fpart1,"max",3)==0||strncmp(fpart1,"min",3)==0) {
        fomulain[p-fomulain]=']';
        fomulain[index]='[';
      } else {
        fomulain[p-fomulain]='}';
        fomulain[index]='{';
      }
      break;
    case 5:
      i=0;
      if(fpart1[0]=='i') {
        i++;
      }
      if(fpart1[1]=='d') {
        i++;
      }
      if(fpart1[2]=='0') {
        i++;
      }
      if(fpart1[3]=='1') {
        i++;
      }
      /* loge( at expression start was missing (asymmetry with the
         operator-preceded switch); sqrt( is new (parity plan 3.1) */
      i1=0;
      if(fpart1[0]=='l'&&fpart1[1]=='o'&&fpart1[2]=='g'&&fpart1[3]=='e')i1=4;
      i2=0;
      if(fpart1[0]=='s'&&fpart1[1]=='q'&&fpart1[2]=='r'&&fpart1[3]=='t')i2=4;
      /* id0v( -- multi-arg intrinsic (parity plan 3.1) */
      if (i==4||i1==4||i2==4||strncmp(fpart1,"id0v",4)==0) {
        fomulain[p-fomulain]=']';
        fomulain[index]='[';
      } else {
        fomulain[p-fomulain]='}';
        fomulain[index]='{';
      }
      break;
    case 6:
      /* log10( / round( -- intrinsics (parity plan 3.1) */
      i=0;
      if(fpart1[0]=='l'&&fpart1[1]=='o'&&fpart1[2]=='g'&&fpart1[3]=='1'&&fpart1[4]=='0')i=5;
      i1=0;
      if(fpart1[0]=='r'&&fpart1[1]=='o'&&fpart1[2]=='u'&&fpart1[3]=='n'&&fpart1[4]=='d')i1=5;
      if (i==5||i1==5) {
        fomulain[p-fomulain]=']';
        fomulain[index]='[';
      } else {
        fomulain[p-fomulain]='}';
        fomulain[index]='{';
      }
      break;
    case 7:
      /* trunc0( / truncb( -- intrinsics (parity plan 3.1) */
      if (strncmp(fpart1,"trunc0",6)==0||strncmp(fpart1,"truncb",6)==0) {
        fomulain[p-fomulain]=']';
        fomulain[index]='[';
      } else {
        fomulain[p-fomulain]='}';
        fomulain[index]='{';
      }
      break;
    default:
      fomulain[p-fomulain]='}';
      fomulain[index]='{';
      break;
    }
  }
  return 1;
}

int sum_extract(char *formula) {
  char *t1=")",*t2="(";
  int i,l,count=0;
  l=strchr(formula,'(')-formula;//ha_cgefind(formula,"(");

  for (i=l; formula[i]; i++) {
    if (formula[i]==*t1) {
      count=count-1;
    }
    if (formula[i]==*t2) {
      count=count+1;
    }
    if (count==0) {
      break;
    }
  }
  formula[i+1]='\0';
  return 1;
}

int eq_replace_linvar(char *formulain,int linindx) {
  char *p,*p1,*p2,leftline[TABREADLINE];//,*p3
  int np,i,l,pcheck;//,pl,mn;//d,
  np=str_count_ci(formulain,"p_");
  p=&formulain[0];
  p2=&formulain[0];
  for (i=0; i<np; i++) {
    p=p+str_find_ci(p,"p_");
    if(p==&formulain[0]) {
      pcheck=0;
      p1=p;
      while (pcheck==0) {
        p1=strpbrk(p1,"+-");
        if (p1==NULL) {
          strcpy(leftline,p2);
        } else {
          strncpy(leftline,p2,p1-p2);
          leftline[p1-p2]='\0';
        }
        if(str_count_ci(leftline,"(")==str_count_ci(leftline,")")) {
          pcheck=1;
        }
        p1++;
      }
      pcheck=0;
      p1=&leftline[0];
      p2=p1;
      l=0;
      while (p1!=NULL) {
        p1=strpbrk(p1,"+-");
        if (p1==NULL) {
          break;
        }
        if (str_count_ci(leftline,"(")==0&&str_count_ci(leftline,")")==0) {
          p2=p1;
          p1++;
          l++;
        } else {
          if (str_count_ci(p1,"(")!=str_count_ci(p1,")")) {
            l++;
            p1++;
          } else {
            p2=p1;
            p1++;
            l++;
          }
        }
      }
      break;
    } else {
      if(*(p-1)=='*'||*(p-1)=='+'||*(p-1)=='-'||*(p-1)=='('||*(p-1)==',') {
        pcheck=0;
        p1=p;
        while (pcheck==0) {
          p1=strpbrk(p1,"+-");
          if (p1==NULL) {
            strcpy(leftline,p2);
          } else {
            strncpy(leftline,p2,p1-p2);
            leftline[p1-p2]='\0';
          }
          if(str_count_ci(leftline,"(")==str_count_ci(leftline,")")) {
            pcheck=1;
          }
          p1++;
        }
        pcheck=0;
        p1=&leftline[0];
        p2=p1;
        l=0;
        while (p1!=NULL) {
          p1=strpbrk(p1,"+-");
          if (p1==NULL) {
            break;
          }
          if (str_count_ci(leftline,"(")==0&&str_count_ci(leftline,")")==0) {
            p2=p1;
            p1++;
            l++;
          } else {
            if (str_count_ci(p1,"(")!=str_count_ci(p1,")")) {
              l++;
              p1++;
            } else {
              p2=p1;
              p1++;
              l++;
            }
          }
        }
        break;
      }
    }
    p++;
  }
  strcpy(formulain,p2);
  return 1;
}

int eq_zero_linvar(char *formulain,int linindx) {
  char *p,*p1=NULL;//,line1[TABREADLINE];
  int np,i,l=linindx;//,l,d;
  np=str_count_ci(formulain,"p_");
  p=&formulain[0];
  for (i=0; i<np; i++) {
    p=p+str_find_ci(p,"p_");

    if(p==&formulain[0]) {
      p1=strpbrk(p,"}+*-/^)");
      if (p1==NULL) {
        p1=strchr(p,'\0');
      } else if (*p1=='}') {
        p1++;
      }
      if (i==l) {
        p=p1;
      } else {
        if (*p1=='\0') {
          *p='0';
          p++;
          *p='\0';
        } else {
          memmove(p+1,p1,strlen(p1)+1);  //formulain[l-d+1]='\0';
          *p='0';
        }
      }
    } else {
      if(*(p-1)=='*'||*(p-1)=='+'||*(p-1)=='-'||*(p-1)=='('||*(p-1)==',') {
        p1=strpbrk(p,"}+*-/^)");
        if (p1==NULL) {
          p1=strchr(p,'\0');
        } else if (*p1=='}') {
          p1++;
        }
        if (i==l) {
          p=p1;
        } else {
          if (*p1=='\0') {
            *p++='0';
            p++;
            *p='\0';
          } else {
            memmove(p+1,p1,strlen(p1)+1);  //formulain[l-d+1]='\0';
            *p++='0';
          }
        }
      } else {
        p++;
        l++;
      }
    }
  }
  return 1;
}

offset_t data_read_files(char *fname, int niodata, cmf_file_entry *iodata, char *commsyntax,set_def *sets,dim_t nset, set_element *set_elems,array_def *coefs,offset_t ncof, elem_store *coef_store,offset_t ncofele,array_def *vars,offset_t nvar, elem_store *var_store,offset_t nvarele) {
  FILE * filehandle, * filehandle1;
  char line[DATREADLINE]="\0",linecopy[DATREADLINE],line1[TABREADLINE],*p,*p1;//,line1[DATREADLINE]
  char vname[NAMESIZE],setindx[NAMESIZE],header[NAMESIZE],*vname1=NULL,argu[TABREADLINE],varset[MAXVARDIM][NAMESIZE];//,varindx[MAXVARDIM][NAMESIZE],varindx1[MAXVARDIM][NAMESIZE];
  offset_t j=0,i,l1,count1,recount,recount1,n,n1,n2,dims,antidim[MAXVARDIM],dim[MAXVARDIM],offset[MAXVARDIM],supsetid[MAXVARDIM],index[MAXVARDIM];
  int count2;//0 no finding;1 find var; 2 find header cof; 4 find header var;
  char *readitem=NULL,*copyvar;
  solve_real val;
  int k0,k1;

  filehandle = fopen(fname,"r");
  if(filehandle==NULL){
    printf("Error: cannot open %s\n",fname);
    return -1;
  }

  while (tab_next_statement(commsyntax,filehandle,line,DATREADLINE)) {
    strcpy(linecopy,line);
    /* Read (by_elements) statements carry CHARACTER data for set
       mappings (manual 11.9.1a) and are consumed by
       mapping_values_read; the numeric reader must not see them
       (the "(" would route them into the indexed-read branch) */
    if (str_find_ci(line,"by_elements")>-1) continue;
    /* Read (IfHeaderExists) qualifier (manual 10.6/11.11.8; plan 3.9):
       the read runs only when the header is on the file; an absent
       header is NOT an error -- values set up by earlier statements
       stand. Header-form reads only (the headerless forms error
       below regardless). */
    int ifhdr=0;
    if (str_find_ci(line,"(ifheaderexists)")>-1) {
      ifhdr=1;
      str_replace_first(line,"(ifheaderexists)","");
      while (str_replace_all(line,"  "," "));
    }
    if (str_find_ci(line,"from terminal")>-1) {
      printf("Error: Read ... from terminal is not supported (read from a file instead)\n");
      return -1;
    }
    k0=str_find_ci(line,"from file ");
    if(k0>-1) {
      k1=str_find_ci(line+k0+10," ");
      strncpy(line1,line+k0+10,k1);
      line1[k1]='\0';
      for (k0=0; k0<niodata; k0++) if (strcmp(line1,iodata[k0].logname)==0) {
          break;
        }
      if(k0==niodata){
        printf("Error: cannot open file %s named in the CMF file\n",line1);
        return -1;
      }
    }
    if (strchr(line,'(')==NULL) {
      readitem = strtok(line," ");
      readitem = strtok(NULL," ");
      if (readitem==NULL) {
        printf("Error: malformed Read statement: %s\n",linecopy);
        return -1;
      }
      strcpy(vname,readitem);
      readitem = strtok(NULL,"\"");
      readitem = strtok(NULL,"\"");
      if (readitem==NULL) {
        printf("Error: Read without a header is not supported (use 'Read X from file <log> header \"H\"'): %s\n",linecopy);
        return -1;
      }
      strcpy(header,readitem);
      n1=0;
      while (header[n1]!='\0'){
        if(header[n1]==' '){
          header[n1]='\0';
          break;
        }
        n1++;
      }
      count2=0;
      for (i=0; i<nvar; i++) {
        vname1= strtok(vars[i].cofname,"(");
        if (strcmp(vname1,vname)==0) { //,vsize ha_var[i].cofname
          vars[i].suplval=true;
          dims=1;
          for (n1=0; n1<vars[i].size; n1++) {
            dims=dims*sets[vars[i].setid[n1]].size;
          }
          antidim[0]=sets[vars[i].setid[1]].size;
          antidim[1]=1;
          antidim[2]=antidim[0]*sets[vars[i].setid[0]].size;
          if (vars[i].size>3) {
            for (n1=3; n1<vars[i].size; n1++) {
              antidim[n1]=antidim[n1-1]*sets[vars[i].setid[n1-1]].size;
            }
          }
          strcpy(vname,readitem);
          n1=0;
          while (vname[n1]!='\0'){
            if(vname[n1]==' '){
              vname[n1]='\0';
              break;
            }
            n1++;
          }
          filehandle1 = fopen(iodata[k0].filname,"r");
          if(filehandle1==NULL){
            printf("Error: cannot open data file %s\n",iodata[k0].filname);
            return -1;
          }
          while (fgets(line,DATREADLINE,filehandle1)) {
            readitem = strtok(line,"\"");
            readitem = strtok(NULL,"\"");
            if (readitem != NULL) {
              if (strcmp(header,vname) == 0) {
                count2=4;
                recount=0;
                while (fgets(line,DATREADLINE,filehandle1)) {
                  count1=0;
                  copyvar=line;
                  while (*copyvar==' ' || *copyvar=='\t') {
                    copyvar++;  // skip whitespaces
                  }
                  if (*copyvar=='\r') {
                    copyvar++;
                  }
                  if (*copyvar=='\n') {
                    count1=1;
                  }
                  if (count1!=1) {
                    n=str_count_char(line,',');
                    if (n==0) {
                      readitem = strtok(line,"\n");
                      val=atof(readitem);//sscanf(readitem, "%lf", &val);
                      l1=recount;
                      if(vars[i].size>1) {
                        for (n1=vars[i].size-1; n1>1; n1--) {
                          index[n1]=(offset_t) l1/antidim[n1];
                          l1=l1-antidim[n1]*index[n1];
                        }
                        index[0]=(offset_t) l1/antidim[0];
                        l1=l1-antidim[0]*index[0];
                        index[1]=(offset_t) l1/antidim[1];
                        l1=l1-antidim[1]*index[1];
                      } else {
                        index[0]=recount;
                      }
                      recount1=0;
                      for (n1=0; n1<vars[i].size; n1++) {
                        recount1=recount1+index[n1]*vars[i].strides[n1];
                      }
                      var_store[vars[i].offset+recount1].value=val;
                      recount++;
                      if (recount>=dims) {
                        break;
                      }
                    } else {
                      for (n2=0; n2<n; n2++) {
                        if (n2==0) {
                          readitem = strtok(line,",");
                        } else {
                          readitem = strtok(NULL,",");
                        }
                        val=atof(readitem);//sscanf(readitem, "%lf", &val);
                        l1=recount;
                        if(vars[i].size>1) {
                          for (n1=vars[i].size-1; n1>1; n1--) {
                            index[n1]=(offset_t) l1/antidim[n1];
                            l1=l1-antidim[n1]*index[n1];
                          }
                          index[0]=(offset_t) l1/antidim[0];
                          l1=l1-antidim[0]*index[0];
                          index[1]=(offset_t) l1/antidim[1];
                          l1=l1-antidim[1]*index[1];
                        } else {
                          index[0]=recount;
                        }
                        recount1=0;
                        for (n1=0; n1<vars[i].size; n1++) {
                          recount1=recount1+index[n1]*vars[i].strides[n1];
                        }
                        var_store[vars[i].offset+recount1].value=val;
                        recount++;
                      }
                      readitem = strtok(NULL,"\n");
                      val=atof(readitem);//sscanf(readitem, "%lf", &val);
                      l1=recount;
                      if(vars[i].size>1) {
                        for (n1=vars[i].size-1; n1>1; n1--) {
                          index[n1]=(offset_t) l1/antidim[n1];
                          l1=l1-antidim[n1]*index[n1];
                        }
                        index[0]=(offset_t) l1/antidim[0];
                        l1=l1-antidim[0]*index[0];
                        index[1]=(offset_t) l1/antidim[1];
                        l1=l1-antidim[1]*index[1];
                      } else {
                        index[0]=recount;
                      }
                      recount1=0;
                      for (n1=0; n1<vars[i].size; n1++) {
                        recount1=recount1+index[n1]*vars[i].strides[n1];
                      }
                      var_store[vars[i].offset+recount1].value=val;
                      recount++;
                      if (recount>=dims) {
                        break;
                      }
                    }
                  }

                }
                break;
              }

            }
            if(feof(filehandle1)) {
              break;
            }
          }
          if(count2!=4){
            if(ifhdr) {
              /* absent header: skip the read silently (11.11.8) */
              fclose(filehandle1);
              break;
            }
            printf("Error: header \"%s\" not found in the data file\n",header);
            return -1;
          }
          fclose(filehandle1);
          break;
        }
      }
      if(i<nvar)count2=1;
      for (i=0; i<ncof; i++) {
        vname1= strtok(coefs[i].cofname,"(");
        if (strcmp(vname1,vname)==0) { //,vsize ha_cof[i].cofname
          coefs[i].suplval=true;
          dims=1;
          for (n1=0; n1<coefs[i].size; n1++) {
            dims=dims*sets[coefs[i].setid[n1]].size;
          }
          antidim[0]=sets[coefs[i].setid[1]].size;
          antidim[1]=1;
          antidim[2]=antidim[0]*sets[coefs[i].setid[0]].size;
          if (coefs[i].size>3) {
            for (n1=3; n1<coefs[i].size; n1++) {
              antidim[n1]=antidim[n1-1]*sets[coefs[i].setid[n1-1]].size;
            }
          }
          filehandle1 = fopen(iodata[k0].filname,"r");
          if(filehandle1==NULL){
            printf("Error: cannot open data file %s\n",iodata[k0].filname);
            return -1;
          }
          while (fgets(line,DATREADLINE,filehandle1)) {
            readitem = strtok(line,"\"");
            readitem = strtok(NULL,"\"");
            if (readitem != NULL) {
              strcpy(vname,readitem);
              n1=0;
              while (vname[n1]!='\0'){
                if(vname[n1]==' '){
                  vname[n1]='\0';
                  break;
                }
                n1++;
              }
              if (strncmp(vname,header,strlen(header))==0&&strlen(header)==strlen(vname)) {
                count2=2;
                recount=0;
                while (fgets(line,DATREADLINE,filehandle1)) {
                  count1=0;
                  copyvar=line;
                  while (*copyvar==' ' || *copyvar=='\t') {
                    copyvar++;  // skip whitespaces
                  }
                  if (*copyvar=='\r') {
                    copyvar++;
                  }
                  if (*copyvar=='\n') {
                    count1=1;
                  }
                  if (count1!=1) {
                    n=str_count_char(line,',');
                    if (n==0) {
                      readitem = strtok(line,"\n");
                      val=atof(readitem);//sscanf(readitem, "%lf", &val);
                      l1=recount;
                      if(coefs[i].size>1) {
                        for (n1=coefs[i].size-1; n1>1; n1--) {
                          index[n1]=(offset_t) l1/antidim[n1];
                          l1=l1-antidim[n1]*index[n1];
                        }
                        index[0]=(offset_t) l1/antidim[0];
                        l1=l1-antidim[0]*index[0];
                        index[1]=(offset_t) l1/antidim[1];
                        l1=l1-antidim[1]*index[1];
                      } else {
                        index[0]=recount;
                      }
                      recount1=0;
                      for (n1=0; n1<coefs[i].size; n1++) {
                        recount1=recount1+index[n1]*coefs[i].strides[n1];
                      }
                      coef_store[coefs[i].offset+recount1].value=val;
                      recount++;
                      if (recount>=dims) {
                        break;
                      }
                    } else {
                      for (n2=0; n2<n; n2++) {
                        if (n2==0) {
                          readitem = strtok(line,",");
                        } else {
                          readitem = strtok(NULL,",");
                        }
                        val=atof(readitem);//sscanf(readitem, "%lf", &val);
                        l1=recount;
                        if(coefs[i].size>1) {
                          for (n1=coefs[i].size-1; n1>1; n1--) {
                            index[n1]=(offset_t) l1/antidim[n1];
                            l1=l1-antidim[n1]*index[n1];
                          }
                          index[0]=(offset_t) l1/antidim[0];
                          l1=l1-antidim[0]*index[0];
                          index[1]=(offset_t) l1/antidim[1];
                          l1=l1-antidim[1]*index[1];
                        } else {
                          index[0]=recount;
                        }
                        recount1=0;
                        for (n1=0; n1<coefs[i].size; n1++) {
                          recount1=recount1+index[n1]*coefs[i].strides[n1];
                        }
                        coef_store[coefs[i].offset+recount1].value=val;
                        recount++;
                      }
                      readitem = strtok(NULL,"\n");
                      val=atof(readitem);//sscanf(readitem, "%lf", &val);
                      l1=recount;
                      if(coefs[i].size>1) {
                        for (n1=coefs[i].size-1; n1>1; n1--) {
                          index[n1]=(offset_t) l1/antidim[n1];
                          l1=l1-antidim[n1]*index[n1];
                        }
                        index[0]=(offset_t) l1/antidim[0];
                        l1=l1-antidim[0]*index[0];
                        index[1]=(offset_t) l1/antidim[1];
                        l1=l1-antidim[1]*index[1];
                      } else {
                        index[0]=recount;
                      }
                      recount1=0;
                      for (n1=0; n1<coefs[i].size; n1++) {
                        recount1=recount1+index[n1]*coefs[i].strides[n1];
                      }
                      coef_store[coefs[i].offset+recount1].value=val;
                      recount++;
                      if (recount>=dims) {
                        break;
                      }
                    }
                  }

                }
                break;
              }

            }
            if(feof(filehandle1)) {
              break;
            }
          }
          if(count2!=2){
            if(ifhdr) {
              /* absent header: skip the read silently (11.11.8) */
              fclose(filehandle1);
              break;
            }
            printf("Error: header \"%s\" not found in the data file\n",header);
            return -1;
          }
          /* satisfied conditional read: formulas assigning this
             coefficient are superseded (11.11.8 idiom, plan 3.9) */
          if(ifhdr&&teems_coef_ifhdr!=NULL&&i<ncof)teems_coef_ifhdr[i]=true;
          fclose(filehandle1);
          break;
        }
      }
      if(count2==0&&i==ncof){
        printf("Error: %s is not a declared variable, coefficient, or parameter\n",vname);
        return -1;
      }
    } else {
      n=str_count_char(line,'(');
      readitem = strtok(line,")");
      for (n1=2; n1<n; n1++) {
        readitem = strtok(NULL,")");
      }
      readitem = strtok(NULL,"(");
      strcpy(vname,readitem);
      str_replace_all(vname," ","");
      readitem = strtok(NULL,")");
      if (readitem==NULL) {
        printf("Error: malformed partial Read statement in TAB file\n");
        return -1;
      }
      strcpy(argu,readitem);
      strcat(argu,",");
      readitem = strtok(NULL,"\"");
      readitem = strtok(NULL,"\"");
      if (readitem==NULL) {
        printf("Error: Read without a header is not supported (use 'Read X(...) from file <log> header \"H\"')\n");
        return -1;
      }
      strcpy(header,readitem);
      n1=0;
      while (header[n1]!='\0'){
        if(header[n1]==' '){
          header[n1]='\0';
          break;
        }
        n1++;
      }
      strcpy(line,linecopy);
      n=str_count_char(argu,',');
      for (n1=0; n1<n; n1++) {
        if (n1==0) {
          readitem = strtok(argu,",");
        } else {
          readitem = strtok(NULL,",");
        }
        strcpy(setindx,",");
        strcat(setindx,readitem);
        strcat(setindx,",");
        p=strstr(line,setindx);
        p1=strchr(p+1,',');
        p=strchr(p1,')');
        strncpy(varset[n1],p1+1,p-p1-1);
        varset[n1][p-p1-1]='\0';
      }
      count2=0;
      for (i=0; i<ncof; i++) {
        vname1= strtok(coefs[i].cofname,"(");
        if (strcmp(vname1,vname)==0) { //,vsize ha_var[i].cofname
          dims=1;
          for (n1=0; n1<coefs[i].size; n1++) {
            if (strcmp(varset[n1],sets[coefs[i].setid[n1]].setname)==0) {
              dims=dims*sets[coefs[i].setid[n1]].size;
              dim[n1]=sets[coefs[i].setid[n1]].size;
              offset[n1]=sets[coefs[i].setid[n1]].offset;
              supsetid[n1]=0;
            } else {
              for (n=0; n<nset; n++) {
                if (strcmp(varset[n1],sets[n].setname)==0) {
                  dims=dims*sets[n].size;
                  dim[n1]=sets[n].size;
                  offset[n1]=sets[n].offset;
                  for(l1=0; l1<MAXSUPSET; l1++)if(strcmp(sets[sets[n].subsetid[l1]].setname,sets[coefs[i].setid[n1]].setname)==0) {
                      supsetid[n1]=l1;
                      break;
                    }
                  break;
                }
              }
            }
          }

          antidim[0]=dim[1];
          antidim[1]=1;
          antidim[2]=antidim[0]*dim[0];
          if (coefs[i].size>3) {
            for (n1=3; n1<coefs[i].size; n1++) {
              antidim[n1]=antidim[n1-1]*dim[n1-1];
            }
          }
          filehandle1 = fopen(iodata[k0].filname,"r");
          if(filehandle1==NULL){
            printf("Error: cannot open data file %s\n",iodata[k0].filname);
            return -1;
          }
          while (fgets(line,DATREADLINE,filehandle1)) {
            readitem = strtok(line,"\"");
            readitem = strtok(NULL,"\"");
            if (readitem != NULL) {
              strcpy(vname1,readitem);              //printf("read %s\n",readitem);
              n1=0;
              while (vname1[n1]!='\0'){
                if(vname1[n1]==' '){
                  vname1[n1]='\0';
                  break;
                }
                n1++;
              }
              if (strncmp(vname1,header,strlen(header)) == 0&&strlen(header)==strlen(vname1)) {
                logmsg(2,"dim %s\n",readitem);
                count2=2;
                recount=0;
                while (fgets(line,DATREADLINE,filehandle1)) {
                  count1=0;
                  copyvar=line;
                  while (*copyvar==' ' || *copyvar=='\t') {
                    copyvar++;  // skip whitespaces
                  }
                  if (*copyvar=='\r') {
                    copyvar++;
                  }
                  if (*copyvar=='\n') {
                    count1=1;
                  }
                  if (count1!=1) {
                    n=str_count_char(line,',');
                    if (n==0) {
                      readitem = strtok(line,"\n");
                      val=atof(readitem);//sscanf(readitem, "%lf", &val);
                      l1=recount;
                      if (coefs[i].size>1) {
                        for (n1=coefs[i].size-1; n1>1; n1--) {
                          index[n1]=(offset_t) l1/antidim[n1];
                          l1=l1-antidim[n1]*index[n1];
                        }
                        index[0]=(offset_t) l1/antidim[0];
                        l1=l1-antidim[0]*index[0];
                        index[1]=(offset_t) l1/antidim[1];
                        l1=l1-antidim[1]*index[1];
                      } else {
                        index[0]=recount;
                      }
                      recount1=0;
                      for (n1=0; n1<coefs[i].size; n1++) {
                        recount1=recount1+set_elems[offset[n1]+index[n1]].superset_pos[supsetid[n1]]*coefs[i].strides[n1];
                      }
                      coef_store[coefs[i].offset+recount1].value=val;
                      recount++;
                      if (recount>=dims) {
                        break;
                      }
                    } else {
                      for (n2=0; n2<n; n2++) {
                        if (n2==0) {
                          readitem = strtok(line,",");
                        } else {
                          readitem = strtok(NULL,",");
                        }
                        val=atof(readitem);//sscanf(readitem, "%lf", &val);
                        l1=recount;
                        if (coefs[i].size>1) {
                          for (n1=coefs[i].size-1; n1>1; n1--) {
                            index[n1]=(offset_t) l1/antidim[n1];
                            l1=l1-antidim[n1]*index[n1];
                          }
                          index[0]=(offset_t) l1/antidim[0];
                          l1=l1-antidim[0]*index[0];
                          index[1]=(offset_t) l1/antidim[1];
                          l1=l1-antidim[1]*index[1];
                        } else {
                          index[0]=recount;
                        }
                        recount1=0;
                        for (n1=0; n1<coefs[i].size; n1++) {
                          recount1=recount1+set_elems[offset[n1]+index[n1]].superset_pos[supsetid[n1]]*coefs[i].strides[n1];
                        }
                        coef_store[coefs[i].offset+recount1].value=val;
                        recount++;
                      }
                      readitem = strtok(NULL,"\n");
                      val=atof(readitem);//sscanf(readitem, "%lf", &val);
                      l1=recount;
                      if (coefs[i].size>1) {
                        for (n1=coefs[i].size-1; n1>1; n1--) {
                          index[n1]=(offset_t) l1/antidim[n1];
                          l1=l1-antidim[n1]*index[n1];
                        }
                        index[0]=(offset_t) l1/antidim[0];
                        l1=l1-antidim[0]*index[0];
                        index[1]=(offset_t) l1/antidim[1];
                        l1=l1-antidim[1]*index[1];
                      } else {
                        index[0]=recount;
                      }
                      recount1=0;
                      for (n1=0; n1<coefs[i].size; n1++) {
                        recount1=recount1+set_elems[offset[n1]+index[n1]].superset_pos[supsetid[n1]]*coefs[i].strides[n1];
                      }
                      coef_store[coefs[i].offset+recount1].value=val;
                      recount++;
                      if (recount>=dims) {
                        break;
                      }
                    }
                  }

                }
              }

            }
            if(feof(filehandle1)) {
              break;
            }
          }
          if(count2!=2){
            if(ifhdr) {
              /* absent header: skip the read silently (11.11.8) */
              fclose(filehandle1);
              break;
            }
            printf("Error: header \"%s\" not found in the data file\n",header);
            return -1;
          }
          /* satisfied conditional read: formulas assigning this
             coefficient are superseded (11.11.8 idiom, plan 3.9) */
          if(ifhdr&&teems_coef_ifhdr!=NULL&&i<ncof)teems_coef_ifhdr[i]=true;
          fclose(filehandle1);
          break;
        }
      }
      if(i<ncof)count2=1;
      for (i=0; i<nvar; i++) {
        vname1= strtok(vars[i].cofname,"(");
        if (strcmp(vname1,vname)==0) { //,vsize ha_var[i].cofname
          dims=1;
          for (n1=0; n1<vars[i].size; n1++) {
            if (strcmp(varset[n1],sets[vars[i].setid[n1]].setname)==0) {
              dims=dims*sets[vars[i].setid[n1]].size;
              dim[n1]=sets[vars[i].setid[n1]].size;
              offset[n1]=sets[vars[i].setid[n1]].offset;
              supsetid[n1]=0;
            } else {
              for (n=0; n<nset; n++) {
                if (strcmp(varset[n1],sets[n].setname)==0) {
                  dims=dims*sets[n].size;
                  dim[n1]=sets[n].size;
                  offset[n1]=sets[n].offset;
                  for(l1=0; l1<MAXSUPSET; l1++)if(strcmp(sets[sets[n].subsetid[l1]].setname,sets[vars[i].setid[n1]].setname)==0) {
                      supsetid[n1]=l1;
                      break;
                    }
                  break;
                }
              }
            }
          }

          antidim[0]=dim[1];
          antidim[1]=1;
          antidim[2]=antidim[0]*dim[0];
          if (vars[i].size>3) {
            for (n1=3; n1<vars[i].size; n1++) {
              antidim[n1]=antidim[n1-1]*dim[n1-1];
            }
          }
          filehandle1 = fopen(iodata[k0].filname,"r");
          if(filehandle1==NULL){
            printf("Error: cannot open data file %s\n",iodata[k0].filname);
            return -1;
          }
          while (fgets(line,DATREADLINE,filehandle1)) {
            readitem = strtok(line,"\"");
            readitem = strtok(NULL,"\"");
            if (readitem != NULL) {
              strcpy(vname1,readitem);              //printf("read %s\n",readitem);
              n1=0;
              while (vname1[n1]!='\0'){
                if(vname1[n1]==' '){
                  vname1[n1]='\0';
                  break;
                }
                n1++;
              }
              if (strncmp(vname1,header,strlen(header)) == 0&&strlen(header)==strlen(vname1)) {
                count2=4;
                recount=0;
                while (fgets(line,DATREADLINE,filehandle1)) {
                  count1=0;
                  copyvar=line;
                  while (*copyvar==' ' || *copyvar=='\t') {
                    copyvar++;  // skip whitespaces
                  }
                  if (*copyvar=='\r') {
                    copyvar++;
                  }
                  if (*copyvar=='\n') {
                    count1=1;
                  }
                  if (count1!=1) {
                    n=str_count_char(line,',');
                    if (n==0) {
                      readitem = strtok(line,"\n");
                      val=atof(readitem);//sscanf(readitem, "%lf", &val);
                      l1=recount;
                      if (vars[i].size>1) {
                        for (n1=vars[i].size-1; n1>1; n1--) {
                          index[n1]=(offset_t) l1/antidim[n1];
                          l1=l1-antidim[n1]*index[n1];
                        }
                        index[0]=(offset_t) l1/antidim[0];
                        l1=l1-antidim[0]*index[0];
                        index[1]=(offset_t) l1/antidim[1];
                        l1=l1-antidim[1]*index[1];
                      } else {
                        index[0]=recount;
                      }
                      recount1=0;
                      for (n1=0; n1<vars[i].size; n1++) {
                        recount1=recount1+set_elems[offset[n1]+index[n1]].superset_pos[supsetid[n1]]*vars[i].strides[n1];
                      }
                      var_store[vars[i].offset+recount1].value=val;
                      recount++;
                      if (recount>=dims) {
                        break;
                      }
                    } else {
                      for (n2=0; n2<n; n2++) {
                        if (n2==0) {
                          readitem = strtok(line,",");
                        } else {
                          readitem = strtok(NULL,",");
                        }
                        val=atof(readitem);//sscanf(readitem, "%lf", &val);
                        l1=recount;
                        if (vars[i].size>1) {
                          for (n1=vars[i].size-1; n1>1; n1--) {
                            index[n1]=(offset_t) l1/antidim[n1];
                            l1=l1-antidim[n1]*index[n1];
                          }
                          index[0]=(offset_t) l1/antidim[0];
                          l1=l1-antidim[0]*index[0];
                          index[1]=(offset_t) l1/antidim[1];
                          l1=l1-antidim[1]*index[1];
                        } else {
                          index[0]=recount;
                        }
                        recount1=0;
                        for (n1=0; n1<vars[i].size; n1++) {
                          recount1=recount1+set_elems[offset[n1]+index[n1]].superset_pos[supsetid[n1]]*vars[i].strides[n1];
                        }
                        var_store[vars[i].offset+recount1].value=val;
                        recount++;
                      }
                      readitem = strtok(NULL,"\n");
                      val=atof(readitem);//sscanf(readitem, "%lf", &val);
                      l1=recount;
                      if (vars[i].size>1) {
                        for (n1=vars[i].size-1; n1>1; n1--) {
                          index[n1]=(offset_t) l1/antidim[n1];
                          l1=l1-antidim[n1]*index[n1];
                        }
                        index[0]=(offset_t) l1/antidim[0];
                        l1=l1-antidim[0]*index[0];
                        index[1]=(offset_t) l1/antidim[1];
                        l1=l1-antidim[1]*index[1];
                      } else {
                        index[0]=recount;
                      }
                      recount1=0;
                      for (n1=0; n1<vars[i].size; n1++) {
                        recount1=recount1+set_elems[offset[n1]+index[n1]].superset_pos[supsetid[n1]]*vars[i].strides[n1];
                      }
                      var_store[vars[i].offset+recount1].value=val;
                      recount++;
                      if (recount>=dims) {
                        break;
                      }
                    }
                  }

                }
              }

            }
            if(feof(filehandle1)) {
              break;
            }
          }
          if(count2!=4){
            if(ifhdr) {
              /* absent header: skip the read silently (11.11.8) */
              fclose(filehandle1);
              break;
            }
            printf("Error: header \"%s\" not found in the data file\n",header);
            return -1;
          }
          fclose(filehandle1);
          break;
        }
      }
      if(count2==0&&i==ncof){
        printf("Error: %s is not a declared variable, coefficient, or parameter\n",vname);
        return -1;
      }
    }
  }
  fclose(filehandle);
  return j;
}


offset_t sum_parse(char *formulain, char *commsyntax, sum_def *sum_cof,quantifier *arSet,set_def *sets,dim_t nset,dim_t fdim,int j) {
  char *readitem,*p,*p1,*p2,interchar2[TABREADLINE],argu[TABREADLINE];//,line5[TABREADLINE]
  char interchar[TABREADLINE],interchar1[TABREADLINE],line[TABREADLINE],line1[TABREADLINE],line2[TABREADLINE],line3[TABREADLINE],line4[TABREADLINE],tempname[NAMESIZE];
  offset_t i=0,k=0,k1=0,length,ncur=0,ncuri,l,l1,l2,l3,l4,l5,l6,l7,sup;
  length=strlen(formulain);
  readitem=formulain;
  while (i<length) {
    k=str_find_ci(readitem,commsyntax);
    if (k==-1) {
      return 0;
    }
    if (k==0) {
      readitem=formulain+i+k;
      strcpy(line,readitem);
      strcpy(line1,readitem);
      sum_extract(line);
      k1=str_find_ci(line+4,commsyntax);
      if (k1!=-1) {
        i=i+k+4;
        readitem=formulain+i;
      } else {
        strcpy(line2,line);
        strcpy(line4,line);
        sprintf(interchar1, "%d", j);
        interchar[0]='\0';
        if (j<10) {
          strcat(interchar,"gen_sum000");
        }
        if (9<j&&j<100) {
          strcat(interchar,"gen_sum00");
        }
        if (99<j&&j<1000) {
          strcat(interchar,"gen_sum0");
        }
        if (999<j&&j<10000) {
          strcat(interchar,"gen_sum");
        }
        if (i>10000) {
          strcat(interchar,"gen_sum");
          printf("Error: too many sum() terms in one statement\n");
        }
        strcat(interchar,interchar1);
        strcpy(sum_cof[j].sumname,interchar);
        strcat(interchar,"{");
        strcpy(interchar2,interchar);
        strcpy(line,line1);
        p = strtok(line,"(");
        p = strtok(NULL,",");
        strcpy(sum_cof[j].sumindx,p);
        p = strtok(NULL,",");
        sum_cond_parse(p,sum_cof[j].sumindx,&sum_cof[j].cond_mapid,sum_cof[j].cond_rhs);
        for (l7=0; l7<nset; l7++) if(strcmp(p,sets[l7].setname)==0) {
            sum_cof[j].sumsetid=l7;
            break;
          }
        sum_cond_domain_check(&sum_cof[j],sets);

        ncur=str_count_ci(line2, "{");
        l2=0;
        l3=0;
        for (ncuri=0; ncuri<ncur; ncuri++) {
          strcpy(line2,line4);
          p = &line2[l2];
          p = strtok(p,"{");
          p = strtok(NULL,"}");
          l2=p-line2;
          strcpy(argu,p);
          strcat(argu,",");
          l=str_count_ci(argu, ",");
          if (l<2) {
            if(strcmp(sum_dim_identity(p),sum_cof[j].sumindx)!=0) {
              strcat(interchar,sum_cof[j].sumindx);
              for (l4=0; l4<l3; l4++) {
                if(strcmp(sum_dim_identity(p),sum_cof[j].dimnames[l4])==0) {
                  break;
                }
              }
              if (l4==l3||l3==0) {
                strcpy(sum_cof[j].dimnames[l3],sum_dim_identity(p));
                l6=0;
                for (l5=0; l5<fdim-1; l5++) if(strcmp(sum_dim_identity(p),arSet[l5].index_name)==0) {
                    sum_cof[j].setid[l3]=arSet[l5].setid;
                    l6++;
                  }
                if (l6==0) {
                  interchar1[0]='\0';
                  strcat(interchar1,"sum(");
                  strcat(interchar1,sum_dim_identity(p));
                  strcpy(line3,formulain);
                  line3[readitem-formulain]='\0';
                  l7=str_rfind_ci(line3,interchar1);
                  p1=&line3[l7+2];
                  p1 = strtok(p1,",");
                  for (l7=0; l7<nset; l7++) if(strcmp(p1,sets[l7].setname)==0) {
                      sum_cof[j].setid[l3]=l7;
                      break;
                    }
                }
                l3++;
              }
            }
          } else {
            for (l1=0; l1<l; l1++) {
              if(l1==0) {
                p = strtok(argu,",");
              } else {
                p = strtok(NULL,",");
              }
              if(strcmp(sum_dim_identity(p),sum_cof[j].sumindx)!=0) {
                for (l4=0; l4<l3; l4++) {
                  if(strcmp(sum_dim_identity(p),sum_cof[j].dimnames[l4])==0) {
                    break;
                  }
                }
                if (l4==l3||l3==0) {
                  strcpy(sum_cof[j].dimnames[l3],sum_dim_identity(p));
                  strcat(interchar,sum_cof[j].dimnames[l3]);
                  l6=0;
                  for (l5=0; l5<fdim-1; l5++) if(strcmp(sum_dim_identity(p),arSet[l5].index_name)==0) {
                      sum_cof[j].setid[l3]=arSet[l5].setid;
                      l6++;
                      break;
                    }
                  if (l6==0) {
                    interchar1[0]='\0';
                    strcat(interchar1,"sum(");
                    strcat(interchar1,sum_dim_identity(p));
                    strcpy(line3,formulain);
                    line3[readitem-formulain]='\0';
                    l7=str_rfind_ci(line3,interchar1);
                    p1=&line3[l7+2];
                    p2=strchr(p1,',');
                    strncpy(tempname,p1,p2-p1);
                    tempname[p2-p1]='\0';
                    for (l7=0; l7<nset; l7++) if(strcmp(tempname,sets[l7].setname)==0) {
                        sum_cof[j].setid[l3]=l7;
                        break;
                      }
                  }
                  l3++;
                }
              }
            }
          }
        }
        l3=sum_cond_carry_rhs(&sum_cof[j],arSet,fdim,l3,interchar);
        if (interchar[strlen(interchar)-1]==',') {
          interchar[strlen(interchar)-1]='}';
        } else {
          if (interchar[strlen(interchar)-1]=='{') {
            interchar[strlen(interchar)-1]='\0';
          } else {
            strcat(interchar,"}");
          }
        }
        while(str_replace_all(formulain,line4,interchar));
        sum_cof[j].size=l3;
        return 1;
      }
    } else if (formulain[i+k-1]=='+'||formulain[i+k-1]=='-'||formulain[i+k-1]=='*'||formulain[i+k-1]=='/'||formulain[i+k-1]=='^'||formulain[i+k-1]=='('||formulain[i+k-1]==',') {
      readitem=formulain+i+k;
      strcpy(line,readitem);
      strcpy(line1,readitem);
      sum_extract(line);
      k1=str_find_ci(line+4,commsyntax);
      if (k1!=-1) {
        i=i+k+4;
        readitem=formulain+i;
      } else {
        strcpy(line2,line);
        strcpy(line4,line);
        sprintf(interchar1, "%d", j);
        interchar[0]='\0';
        if (j<10) {
          strcat(interchar,"gen_sum000");
        }
        if (9<j&&j<100) {
          strcat(interchar,"gen_sum00");
        }
        if (99<j&&j<1000) {
          strcat(interchar,"gen_sum0");
        }
        if (999<j&&j<10000) {
          strcat(interchar,"gen_sum");
        }
        if (i>10000) {
          strcat(interchar,"gen_sum");
          printf("Error: too many sum() terms in one statement\n");
        }
        strcat(interchar,interchar1);
        strcpy(sum_cof[j].sumname,interchar);
        strcat(interchar,"{");
        strcpy(interchar2,interchar);
        strcpy(line,line1);
        p = strtok(line,"(");
        p = strtok(NULL,",");
        strcpy(sum_cof[j].sumindx,p);
        p = strtok(NULL,",");
        sum_cond_parse(p,sum_cof[j].sumindx,&sum_cof[j].cond_mapid,sum_cof[j].cond_rhs);
        for (l7=0; l7<nset; l7++) if(strcmp(p,sets[l7].setname)==0) {
            sum_cof[j].sumsetid=l7;
            break;
          }
        sum_cond_domain_check(&sum_cof[j],sets);

        ncur=str_count_ci(line2, "{");
        l2=0;
        l3=0;
        for (ncuri=0; ncuri<ncur; ncuri++) {
          strcpy(line2,line4);
          p = &line2[l2];
          p = strtok(p,"{");
          p = strtok(NULL,"}");
          l2=p-line2;
          strcpy(argu,p);
          strcat(argu,",");
          l=str_count_ci(argu, ",");
          if (l<2) {
            if(strcmp(sum_dim_identity(p),sum_cof[j].sumindx)!=0) {
              for (l4=0; l4<l3; l4++) {
                if(strcmp(sum_dim_identity(p),sum_cof[j].dimnames[l4])==0) {
                  break;
                }
              }
              if (l4==l3||l3==0) {
                strcpy(sum_cof[j].dimnames[l3],sum_dim_identity(p));
                strcat(interchar,sum_dim_identity(p));
                l6=0;
                for (l5=0; l5<fdim-1; l5++) if(strcmp(sum_dim_identity(p),arSet[l5].index_name)==0) {
                    sum_cof[j].setid[l3]=arSet[l5].setid;
                    l6++;
                  }
                if (l6==0) {
                  interchar1[0]='\0';
                  strcat(interchar1,"sum(");
                  strcat(interchar1,sum_dim_identity(p));
                  strcpy(line3,formulain);
                  line3[readitem-formulain]='\0';
                  l7=str_rfind_ci(line3,interchar1);
                  p1=&line3[l7+2];
                  p1 = strtok(p1,",");
                  for (l7=0; l7<nset; l7++) if(strcmp(p1,sets[l7].setname)==0) {
                      sum_cof[j].setid[l3]=l7;
                      break;
                    }
                }
                l3++;
              }
            }
          } else {
            for (l1=0; l1<l; l1++) {
              if(l1==0) {
                p = strtok(argu,",");
              } else {
                p = strtok(NULL,",");
              }
              if(strcmp(sum_dim_identity(p),sum_cof[j].sumindx)!=0) {
                for (l4=0; l4<l3; l4++) {
                  if(strcmp(sum_dim_identity(p),sum_cof[j].dimnames[l4])==0) {
                    break;
                  }
                }
                if (l4==l3||l3==0) {
                  strcpy(sum_cof[j].dimnames[l3],sum_dim_identity(p));
                  strcat(interchar,sum_cof[j].dimnames[l3]);
                  strcat(interchar,",");
                  l6=0;
                  for (l5=0; l5<fdim-1; l5++) if(strcmp(sum_dim_identity(p),arSet[l5].index_name)==0) {
                      sum_cof[j].setid[l3]=arSet[l5].setid;
                      l6++;
                    }
                  if (l6==0) {
                    interchar1[0]='\0';
                    strcat(interchar1,"sum(");
                    strcat(interchar1,sum_dim_identity(p));
                    strcpy(line3,formulain);
                    line3[readitem-formulain]='\0';
                    l7=str_rfind_ci(line3,interchar1);
                    p1=&line3[l7+2];
                    p2=strchr(p1,',');
                    strncpy(tempname,p1,p2-p1);
                    tempname[p2-p1]='\0';
                    for (l7=0; l7<nset; l7++) if(strcmp(tempname,sets[l7].setname)==0) {
                        sum_cof[j].setid[l3]=l7;
                        break;
                      }
                  }
                  l3++;
                }
              }
            }
          }
        }
        l3=sum_cond_carry_rhs(&sum_cof[j],arSet,fdim,l3,interchar);
        if (interchar[strlen(interchar)-1]==',') {
          interchar[strlen(interchar)-1]='}';
        } else {
          if (interchar[strlen(interchar)-1]=='{') {
            interchar[strlen(interchar)-1]='\0';
          } else {
            strcat(interchar,"}");
          }
        }
        while(str_replace_all(formulain,line4,interchar));
        sum_cof[j].size=l3;
        return 1;
      }
    } else {
      i=i+k+4;
      readitem=formulain+i;
    }
  }
  return 0;
}

int sum_count(char *formulain, char *commsyntax) {
  char *readitem;
  int j=0,i=0,k=0,length;
  length=strlen(formulain);
  readitem=formulain;
  while (i<length) {
    k=str_find_ci(readitem,commsyntax);
    if (k==-1) {
      break;
    }
    if (k==0) {
      j++;
      i=i+k+4;
      readitem=formulain+i;
    } else if (formulain[i+k-1]=='+'||formulain[i+k-1]=='-'||formulain[i+k-1]=='*'||formulain[i+k-1]=='/'||formulain[i+k-1]=='^'||formulain[i+k-1]=='('||formulain[i+k-1]==','||formulain[i+k-1]=='='||formulain[i+k-1]=='['||formulain[i+k-1]=='{'||formulain[i+k-1]==' ') {
      j++;
      i=i+k+4;
      readitem=formulain+i;
    } else {
      i=i+k+4;
      readitem=formulain+i;
    }
  }
  return j;
}

/* Name-uniqueness validation (manual 11.2.1; closes the 12.2.2 name-
   resolution spec pass): "duplication of a name for two different
   purposes is not allowed" -- names are case-insensitive, so
   coefficient POP and variable pop are the SAME name. The binder
   searches coefficients, variables and sets from several contexts
   (a PostSim formula RHS admits both coefficients and variables), so
   a collision silently binds whichever list is searched first -- the
   class behind the POP-vs-pop PostSim finding (solution zeros instead
   of data). Enforced across and within coefficient/variable/set.
   Logical FILE and EQUATION names are exempt: the corpus write-all
   pattern names output Files after the coefficients they dump, and a
   real corpus TAB shares a coefficient/equation name; neither
   namespace cross-binds here. Reserved function words are rejected as
   identifiers (the expression machinery would misparse them).
   Returns 0 ok, -1 fatal. */
int names_validate(set_def *sets, dim_t nset, array_def *coefs, offset_t ncof, array_def *vars, offset_t nvar, map_def *maps, dim_t nmap) {
  static const char *reserved[]={"sum","if","abs","max","min","sqrt","exp","loge","log10","id01","id0v","round","trunc0","truncb","prod","maxs","mins","random","normal","cumnormal","lognormal","cumlognormal","gperf","gperfc","ras_matrix","all",NULL};
  offset_t i,j;
  dim_t k;
  int r;
  for(i=0; i<ncof; i++) {
    for(j=0; j<nvar; j++)if(strcmp(coefs[i].cofname,vars[j].cofname)==0) {
        printf("Error: name %s is declared as both a coefficient and a variable; names are case-insensitive and must be unique (manual 11.2.1)\n",coefs[i].cofname);
        return -1;
      }
    for(k=0; k<nset; k++)if(strcmp(coefs[i].cofname,sets[k].setname)==0) {
        printf("Error: name %s is declared as both a coefficient and a set (manual 11.2.1)\n",coefs[i].cofname);
        return -1;
      }
    for(j=i+1; j<ncof; j++)if(strcmp(coefs[i].cofname,coefs[j].cofname)==0) {
        printf("Error: coefficient %s is declared more than once (manual 11.2.1)\n",coefs[i].cofname);
        return -1;
      }
    for(r=0; reserved[r]!=NULL; r++)if(strcmp(coefs[i].cofname,reserved[r])==0) {
        printf("Error: coefficient name %s is a reserved word (manual 11.2.1)\n",coefs[i].cofname);
        return -1;
      }
  }
  for(i=0; i<nvar; i++) {
    for(k=0; k<nset; k++)if(strcmp(vars[i].cofname,sets[k].setname)==0) {
        printf("Error: name %s is declared as both a variable and a set (manual 11.2.1)\n",vars[i].cofname);
        return -1;
      }
    for(j=i+1; j<nvar; j++)if(strcmp(vars[i].cofname,vars[j].cofname)==0) {
        printf("Error: variable %s is declared more than once (manual 11.2.1)\n",vars[i].cofname);
        return -1;
      }
    for(r=0; reserved[r]!=NULL; r++)if(strcmp(vars[i].cofname,reserved[r])==0) {
        printf("Error: variable name %s is a reserved word (manual 11.2.1)\n",vars[i].cofname);
        return -1;
      }
  }
  for(k=0; k<nset; k++)for(r=0; reserved[r]!=NULL; r++)if(strcmp(sets[k].setname,reserved[r])==0) {
      printf("Error: set name %s is a reserved word (manual 11.2.1)\n",sets[k].setname);
      return -1;
    }
  for(k=0; k<nmap; k++) {
    for(i=0; i<ncof; i++)if(strcmp(maps[k].mapname,coefs[i].cofname)==0) {
        printf("Error: name %s is declared as both a mapping and a coefficient (manual 11.2.1)\n",maps[k].mapname);
        return -1;
      }
    for(i=0; i<nvar; i++)if(strcmp(maps[k].mapname,vars[i].cofname)==0) {
        printf("Error: name %s is declared as both a mapping and a variable (manual 11.2.1)\n",maps[k].mapname);
        return -1;
      }
    for(j=0; j<nset; j++)if(strcmp(maps[k].mapname,sets[j].setname)==0) {
        printf("Error: name %s is declared as both a mapping and a set (manual 11.2.1)\n",maps[k].mapname);
        return -1;
      }
    for(j=k+1; j<nmap; j++)if(strcmp(maps[k].mapname,maps[j].mapname)==0) {
        printf("Error: mapping %s is declared more than once (manual 11.2.1)\n",maps[k].mapname);
        return -1;
      }
    for(r=0; reserved[r]!=NULL; r++)if(strcmp(maps[k].mapname,reserved[r])==0) {
        printf("Error: mapping name %s is a reserved word (manual 11.2.1)\n",maps[k].mapname);
        return -1;
      }
  }
  return 0;
}

/* mark coefficients declared inside POSTSIM sections (names recorded
   by the split) -- drives the 12.2.2 LHS rule and the 12.2.3
   Read-target rule; teems_coef_is_ps stays NULL when the TAB has no
   sections (zero cost) */
void postsim_mark_coefs(array_def *coefs, offset_t ncof) {
  offset_t i;
  int k;
  if(teems_ps_ncoefs==0)return;
  free(teems_coef_is_ps);
  teems_coef_is_ps= (bool *) calloc (ncof+1,sizeof(bool));
  for(i=0; i<ncof; i++)
    for(k=0; k<teems_ps_ncoefs; k++)
      if(strcmp(coefs[i].cofname,teems_ps_coefnames[k])==0) {
        teems_coef_is_ps[i]=true;
        break;
      }
}

/* PostSim Read execution (manual 12.2.1/12.2.3; Tier 0 residual):
   runs the _ps companion's Read statements once after the solve.
   Targets must be PostSim Coefficients -- never Variables, never
   ordinary Coefficients. Reuses data_read_files against scratch
   stores, then copies only the read coefficients' element ranges into
   the live elem_vals, preserving every post-simulation value.
   Returns the read count (0 = zero cost), -1 on error. */
offset_t postsim_reads_execute(char *psname, int niodata, cmf_file_entry *iodata, set_def *sets, dim_t nset, set_element *set_elems, array_def *coefs, offset_t ncof, offset_t ncofele, array_def *vars, offset_t nvar, offset_t nvarele, elem_value *elem_vals) {
  FILE *f;
  char line[TABREADLINE],name[NAMESIZE],*p;
  offset_t i,j,l,nreads;
  elem_store *cstore,*vstore;
  bool *docopy;
  nreads=tab_count_statements(psname,"read");
  if(nreads==0)return 0;
  docopy= (bool *) calloc (ncof+1,sizeof(bool));
  f=fopen(psname,"r");
  if(f==NULL) {
    free(docopy);
    return -1;
  }
  while (tab_next_statement("read",f,line,TABREADLINE)) {
    /* target = the first token after "read" ("read elements" stays
       with the declarations and never reaches the _ps file) */
    p=line+4;
    while(*p==' ')p++;
    l=0;
    while(*p!=' '&&*p!='\0'&&*p!=';'&&l<NAMESIZE-1)name[l++]=*p++;
    name[l]='\0';
    for(i=0; i<ncof; i++)if(strcmp(coefs[i].cofname,name)==0)break;
    if(i<ncof) {
      if(teems_coef_is_ps==NULL||!teems_coef_is_ps[i]) {
        printf("Error: PostSim Read into ordinary coefficient %s; targets must be PostSim Coefficients (manual 12.2.3)\n",name);
        fclose(f);
        free(docopy);
        return -1;
      }
      docopy[i]=true;
      continue;
    }
    for(j=0; j<nvar; j++)if(strcmp(vars[j].cofname,name)==0)break;
    if(j<nvar) {
      printf("Error: PostSim Read into variable %s; simulation results cannot be changed (manual 12.2.3)\n",name);
      fclose(f);
      free(docopy);
      return -1;
    }
    printf("Error: PostSim Read target %s is not a declared coefficient\n",name);
    fclose(f);
    free(docopy);
    return -1;
  }
  fclose(f);
  cstore= (elem_store *) calloc (ncofele+1,sizeof(elem_store));
  vstore= (elem_store *) calloc (nvarele+1,sizeof(elem_store));
  if(data_read_files(psname,niodata,iodata,"read",sets,nset,set_elems,coefs,ncof,cstore,ncofele,vars,nvar,vstore,nvarele)==-1) {
    free(cstore);
    free(vstore);
    free(docopy);
    return -1;
  }
  for(i=0; i<ncof; i++)if(docopy[i])
      for(l=0; l<coefs[i].nelem; l++)elem_vals[coefs[i].offset+l].value=cstore[coefs[i].offset+l].value;
  free(cstore);
  free(vstore);
  free(docopy);
  return nreads;
}

offset_t tab_count_statements(char *fname, char *commsyntax) {
  FILE * filehandle;
  char line[TABREADLINE]="\0";
  offset_t j=0;
  filehandle = fopen(fname,"r");
  if (filehandle==NULL) return 0;
  while (tab_next_statement(commsyntax,filehandle,line,TABREADLINE)) {
    if (strstr(line,"(default")==NULL) {
      j++;
    }
  }
  fclose(filehandle);
  return j;
}

/* MAPPING declarations (manual 11.9):
     Mapping [(onto)] <name> from <set1> to <set2> ;
   Both sets must already be declared; the (onto) qualifier arms the
   codomain-coverage check in mappings_validate (11.9.3). */
int mappings_read(char *fname, map_def *maps, dim_t nmap, set_def *sets, dim_t nset) {
  FILE *filehandle;
  char line[TABREADLINE]="\0",linecopy[TABREADLINE];
  char *commsyntax="mapping";
  char *readitem=NULL;
  dim_t j=0,i;
  filehandle = fopen(fname,"r");
  while (tab_next_statement(commsyntax,filehandle,line,TABREADLINE)) {
    if (j>=nmap) break;
    strcpy(linecopy,line);
    readitem = strtok(line," ");   /* "mapping" */
    readitem = strtok(NULL," ");
    if (readitem!=NULL&&strcmp(readitem,"(onto)")==0) {
      maps[j].onto=true;
      readitem = strtok(NULL," ");
    }
    if (readitem==NULL||strlen(readitem)>=NAMESIZE) {
      printf("Error: malformed Mapping statement in TAB file: %s\n",linecopy);
      fclose(filehandle);
      MPI_Abort(PETSC_COMM_WORLD,1);
      return -1;
    }
    strcpy(maps[j].mapname,readitem);
    readitem = strtok(NULL," ");
    if (readitem==NULL||strcmp(readitem,"from")!=0) {
      printf("Error: malformed Mapping statement in TAB file: %s\n",linecopy);
      fclose(filehandle);
      MPI_Abort(PETSC_COMM_WORLD,1);
      return -1;
    }
    readitem = strtok(NULL," ");
    for (i=0; i<nset; i++) if (readitem!=NULL&&strcmp(readitem,sets[i].setname)==0) break;
    if (readitem==NULL||i==nset) {
      printf("Error: set %s in Mapping %s is not declared\n",readitem==NULL?"":readitem,maps[j].mapname);
      fclose(filehandle);
      MPI_Abort(PETSC_COMM_WORLD,1);
      return -1;
    }
    maps[j].fromset=i;
    readitem = strtok(NULL," ");
    if (readitem==NULL||strcmp(readitem,"to")!=0) {
      printf("Error: malformed Mapping statement in TAB file: %s\n",linecopy);
      fclose(filehandle);
      MPI_Abort(PETSC_COMM_WORLD,1);
      return -1;
    }
    readitem = strtok(NULL," ;");
    for (i=0; i<nset; i++) if (readitem!=NULL&&strcmp(readitem,sets[i].setname)==0) break;
    if (readitem==NULL||i==nset) {
      printf("Error: set %s in Mapping %s is not declared\n",readitem==NULL?"":readitem,maps[j].mapname);
      fclose(filehandle);
      MPI_Abort(PETSC_COMM_WORLD,1);
      return -1;
    }
    maps[j].toset=i;
    maps[j].values= (dim_t *) calloc (sets[maps[j].fromset].size>0?sets[maps[j].fromset].size:1,sizeof(dim_t));
    j++;
  }
  fclose(filehandle);
  return 0;
}

/* Mapping values via Read (by_elements) (manual 11.9.1a): CHARACTER
   data, one codomain element NAME per domain element, resolved to
   codomain positions here. Also the guard sweep for the deferred
   mapping forms: integer reads of a mapping (11.9.1c), formulas
   assigning a mapping (11.9.1b/d/e), and writes of a mapping
   (11.9.10) are named fatals, not silent mis-binds. */
int mapping_values_read(char *fname, int niodata, cmf_file_entry *iodata, map_def *maps, dim_t nmap, set_def *sets, set_element *set_elems) {
  FILE *filehandle;
  char line[TABREADLINE]="\0",linecopy[TABREADLINE],line1[TABREADLINE],longname[TABREADLINE];
  char *readitem=NULL;
  dim_t j,i,vsize,dim1;
  int k0,k1,nj;
  filehandle = fopen(fname,"r");
  while (tab_next_statement("read",filehandle,line,TABREADLINE)) {
    strcpy(linecopy,line);
    int byele=(str_find_ci(line,"by_elements")>-1);
    if (byele) {
      str_replace_first(line,"(by_elements)","");
      while (str_replace_all(line,"  "," "));
    }
    readitem = strtok(line," ");   /* "read" */
    readitem = strtok(NULL," ");
    if (readitem==NULL) continue;
    for (j=0; j<nmap; j++) if (strcmp(readitem,maps[j].mapname)==0) break;
    if (j==nmap) {
      if (byele) {
        printf("Error: Read (by_elements) target %s is not a declared mapping (manual 11.9.1)\n",readitem);
        fclose(filehandle);
        MPI_Abort(PETSC_COMM_WORLD,1);
        return -1;
      }
      continue;
    }
    if (!byele) {
      printf("Error: integer Read of mapping %s is not supported; use Read (by_elements) (manual 11.9.1)\n",maps[j].mapname);
      fclose(filehandle);
      MPI_Abort(PETSC_COMM_WORLD,1);
      return -1;
    }
    k0=str_find_ci(linecopy,"from file ");
    if (k0<0) {
      printf("Error: Read (by_elements) for mapping %s has no file clause\n",maps[j].mapname);
      fclose(filehandle);
      MPI_Abort(PETSC_COMM_WORLD,1);
      return -1;
    }
    k1=str_find_ci(linecopy+k0+10," ");
    if (k1<0||k1>=(int)sizeof(line1)) {
      printf("Error: malformed Read (by_elements) statement for mapping %s\n",maps[j].mapname);
      fclose(filehandle);
      MPI_Abort(PETSC_COMM_WORLD,1);
      return -1;
    }
    strncpy(line1,linecopy+k0+10,k1);
    line1[k1]='\0';
    for (k0=0; k0<niodata; k0++) if (strcmp(line1,iodata[k0].logname)==0) break;
    if (k0==niodata) {
      printf("Error: cannot open file %s named in the CMF file\n",line1);
      fclose(filehandle);
      MPI_Abort(PETSC_COMM_WORLD,1);
      return -1;
    }
    readitem = strtok(NULL,"\"");
    readitem = strtok(NULL,"\"");
    if (readitem==NULL||strlen(readitem)>=HEADERSIZE) {
      printf("Error: Read (by_elements) for mapping %s needs a header (manual 11.9.1)\n",maps[j].mapname);
      fclose(filehandle);
      MPI_Abort(PETSC_COMM_WORLD,1);
      return -1;
    }
    datafile_read_header_info(readitem,iodata[k0].filname,&vsize,longname,&dim1);
    if (dim1!=sets[maps[j].fromset].size) {
      printf("Error: Read (by_elements) for mapping %s supplies %ld values for the %ld elements of set %s\n",maps[j].mapname,(long)dim1,(long)sets[maps[j].fromset].size,sets[maps[j].fromset].setname);
      fclose(filehandle);
      MPI_Abort(PETSC_COMM_WORLD,1);
      return -1;
    }
    datafile_labels *matvar1= (datafile_labels *) calloc (dim1,sizeof(datafile_labels));
    datafile_read_labels(readitem,iodata[k0].filname,dim1,matvar1);
    for (i=0; i<dim1; i++) {
      nj=0;
      while (matvar1[i].ch[nj]!='\0') {
        matvar1[i].ch[nj]=tolower((int)matvar1[i].ch[nj]);
        nj++;
      }
      dim_t k;
      for (k=0; k<sets[maps[j].toset].size; k++)
        if (strcmp(matvar1[i].ch,set_elems[sets[maps[j].toset].offset+k].setele)==0) break;
      if (k==sets[maps[j].toset].size) {
        printf("Error: %s in the data for mapping %s is not an element of set %s (manual 11.9.2)\n",matvar1[i].ch,maps[j].mapname,sets[maps[j].toset].setname);
        free(matvar1);
        fclose(filehandle);
        MPI_Abort(PETSC_COMM_WORLD,1);
        return -1;
      }
      maps[j].values[i]=k;
    }
    free(matvar1);
    maps[j].has_values=true;
  }
  fclose(filehandle);
  return 0;
}

/* Deferred-form guards (run AFTER names_validate so a name clash gets
   its own message first): a formula whose LHS is a mapping
   (11.9.1b/d/e) or a write of a mapping (11.9.10) would otherwise
   mis-bind silently */
int mapping_use_guards(char *fname, map_def *maps, dim_t nmap) {
  FILE *filehandle;
  char line[TABREADLINE]="\0",linecopy[TABREADLINE];
  char *readitem=NULL;
  dim_t j;
  filehandle = fopen(fname,"r");
  while (tab_next_statement("formula",filehandle,line,TABREADLINE)) {
    strcpy(linecopy,line);
    /* a mapping equality in a Formula QUANTIFIER condition would have
       its RHS read as a numeric constant (0) and filter silently
       wrong; checked here, before the LHS walk below mistakes the
       condition's mapping call for a formula-assigned mapping (M3) */
    {
      int qk=0,qf,qbeg,depth;
      while ((qf=str_find_ci(line+qk,"(all,"))>-1) {
        int hascolon=0;
        char qsave;
        qk+=qf+5;
        qbeg=qk;
        depth=1;
        for (; line[qk]!='\0'&&depth>0; qk++) {
          if (line[qk]=='(') depth++;
          else if (line[qk]==')') depth--;
          else if (line[qk]==':') hascolon=1;
        }
        if (hascolon) {
          dim_t qm;
          char qfind[NAMESIZE+2];
          qsave=line[qk];
          line[qk]='\0';
          for (qm=0; qm<nmap; qm++) {
            sprintf(qfind,"%s(",maps[qm].mapname);
            if (str_find_ci(line+qbeg,qfind)>-1) {
              printf("Error: mapping equalities in Formula quantifier conditions are not supported; move the condition into a sum (manual 11.4.11)\n");
              fclose(filehandle);
              MPI_Abort(PETSC_COMM_WORLD,1);
              return -1;
            }
          }
          line[qk]=qsave;
        }
      }
      strcpy(line,linecopy);
    }
    readitem = strtok(line," ");
    readitem = strtok(NULL," (");
    /* skip quantifier groups ("all,c,com)" -- the "(" is a delimiter)
       and Formula qualifiers ("by_elements)", "initial)", ...) */
    while (readitem!=NULL&&(strncmp(readitem,"all,",4)==0||strchr(readitem,')')!=NULL)) {
      readitem = strtok(NULL," (");
    }
    if (readitem==NULL) continue;
    strtok(readitem,"(=");
    for (j=0; j<nmap; j++) if (strcmp(readitem,maps[j].mapname)==0) {
        printf("Error: formula-assigned mappings are not supported; give mapping %s its values with Read (by_elements) (manual 11.9.1)\n",maps[j].mapname);
        fclose(filehandle);
        MPI_Abort(PETSC_COMM_WORLD,1);
        return -1;
      }
  }
  fclose(filehandle);
  /* updates cannot carry mapping calls yet (corpus: none do); equations
     are lowered in the tab_write_variables rewrite and compiled through
     the mapped linear-variable machinery (design doc M2b) */
  filehandle = fopen(fname,"r");
  while (tab_next_statement("update",filehandle,line,TABREADLINE)) {
    dim_t gm;
    int gpos;
    char gfind[NAMESIZE+2];
    for (gm=0; gm<nmap; gm++) {
      sprintf(gfind,"%s(",maps[gm].mapname);
      gpos=str_find_ci(line,gfind);
      if (gpos==0||(gpos>0&&!isalnum((int)line[gpos-1])&&line[gpos-1]!='_')) {
        printf("Error: mapping-valued indices in Update statements are not supported yet (mapping %s)\n",maps[gm].mapname);
        fclose(filehandle);
        MPI_Abort(PETSC_COMM_WORLD,1);
        return -1;
      }
    }
  }
  fclose(filehandle);
  filehandle = fopen(fname,"r");
  while (tab_next_statement("write",filehandle,line,TABREADLINE)) {
    readitem = strtok(line," ");
    readitem = strtok(NULL," ");
    if (readitem!=NULL&&readitem[0]=='(') readitem = strtok(NULL," ");
    if (readitem==NULL) continue;
    for (j=0; j<nmap; j++) if (strcmp(readitem,maps[j].mapname)==0) {
        printf("Error: writing mapping %s is not supported (manual 11.9.10)\n",maps[j].mapname);
        fclose(filehandle);
        MPI_Abort(PETSC_COMM_WORLD,1);
        return -1;
      }
  }
  fclose(filehandle);
  return 0;
}

/* pre-lower mapping calls "map(i)" to the flat token "map@i" so the
   brace tokenizers never meet a nested index list (design doc M2);
   composition map1(map2(i)) stays a named fatal (11.9.6 deferred) */
void mapping_lower_calls(char *line) {
  dim_t m;
  int k,k2,from;
  char find[NAMESIZE+2];
  for (m=0; m<teems_nmap; m++) {
    sprintf(find,"%s(",teems_maps[m].mapname);
    from=0;
    while ((k=str_find_ci(line+from,find))>-1) {
      k+=from;
      if (k>0&&(isalnum((int)line[k-1])||line[k-1]=='_')) {
        from=k+1;
        continue;
      }
      int op=k+(int)strlen(find)-1;
      for (k2=op+1; line[k2]!='\0'&&line[k2]!=')'&&line[k2]!='('; k2++);
      /* a '(' before the ')' is an unlowered inner call; an '@' in the
         span is an already-lowered one -- both are composition */
      {
        int k3;
        int inner=(line[k2]=='(');
        for (k3=op+1; k3<k2; k3++) if (line[k3]=='@') inner=1;
        if (inner) {
          printf("Error: composition of set mappings is not supported (mapping %s; manual 11.9.6)\n",teems_maps[m].mapname);
          MPI_Abort(PETSC_COMM_WORLD,1);
        }
      }
      if (line[k2]=='\0') {
        from=k+1;
        continue;
      }
      line[op]='@';
      memmove(line+k2,line+k2+1,strlen(line+k2+1)+1);
      from=k+1;
    }
  }
}

/* Identity of an index token for sum carried-dim discovery (design
   doc M2c): a lowered mapping token map@idx is summed or carried by
   its domain index idx -- the mapping routes the OPERAND's dimension,
   not the loop dimension -- so the summed-index test, carried-dim
   dedupe, set resolution and the generated-sum replacement token all
   work on the part after the '@'. */
char *sum_dim_identity(char *p) {
  char *at=strchr(p,'@');
  return at==NULL?p:at+1;
}

/* Parse the ":<condition>" tail of a sum's set token (manual 11.4.11;
   design doc M3).  The one supported form is a mapping equality on
   the summed index -- MAP(idx) = rhs, arriving lowered as
   map@idx=rhs -- and every other condition is a named fatal: before
   M3 a ':' in the set token made the set lookup miss and the sum
   silently expanded over sets[0].  The set part is truncated in
   place for the caller to resolve; the RHS token is kept raw and
   resolves at build/eval time against the frame (outer quantifier
   index) or the mapping's codomain elements. */
void sum_cond_parse(char *settok, const char *sumindx, int *cond_mapid, char *cond_rhs) {
  char *colon,*eq,*lhs;
  int mp=0;
  *cond_mapid=0;
  cond_rhs[0]='\0';
  colon=strchr(settok,':');
  if (colon==NULL) return;
  *colon='\0';
  lhs=colon+1;
  eq=strchr(lhs,'=');
  if (eq==NULL) {
    printf("Error: unsupported sum condition '%s'; only a mapping equality MAPPING(index) = value is supported (manual 11.4.11)\n",lhs);
    MPI_Abort(PETSC_COMM_WORLD,1);
    return;
  }
  *eq='\0';
  lhs=mapping_token_split(lhs,&mp);
  if (mp==0) {
    printf("Error: unsupported sum condition on '%s'; only a mapping equality MAPPING(index) = value is supported (manual 11.4.11)\n",lhs);
    MPI_Abort(PETSC_COMM_WORLD,1);
    return;
  }
  if (strcmp(lhs,sumindx)!=0) {
    printf("Error: the sum condition on mapping %s must test the summed index %s, not %s\n",teems_maps[mp-1].mapname,sumindx,lhs);
    MPI_Abort(PETSC_COMM_WORLD,1);
    return;
  }
  if (eq[1]=='\0'||strlen(eq+1)>=NAMESIZE) {
    printf("Error: malformed sum condition RHS for mapping %s\n",teems_maps[mp-1].mapname);
    MPI_Abort(PETSC_COMM_WORLD,1);
    return;
  }
  *cond_mapid=mp;
  {
    /* keep the RHS raw but unquoted: a fixed codomain element may
       arrive as "ele" */
    char *src=eq+1,*dst=cond_rhs;
    for (; *src!='\0'; src++) if (*src!='\"') *dst++=*src;
    *dst='\0';
  }
}

/* Resolve a sum condition's RHS against an evaluation frame: an index
   name in the frame gives a per-tuple compare (the quantifier must
   range over the codomain exactly), a codomain element name a fixed
   position; anything else is a named fatal. */
void sum_cond_rhs_resolve(int cond_mapid, const char *cond_rhs, quantifier *frame, dim_t nframe, set_def *sets, set_element *set_elems, int *condpos, offset_t *condfix) {
  dim_t l;
  offset_t cs;
  *condpos=-1;
  *condfix=-1;
  if (cond_mapid==0) return;
  for (l=0; l<nframe; l++) if (strcmp(cond_rhs,frame[l].index_name)==0) {
      if (frame[l].setid!=(offset_t)teems_maps[cond_mapid-1].toset) {
        printf("Error: the RHS index %s of the condition on mapping %s does not range over the mapping's codomain set\n",cond_rhs,teems_maps[cond_mapid-1].mapname);
        MPI_Abort(PETSC_COMM_WORLD,1);
        return;
      }
      *condpos=(int)l;
      return;
    }
  cs=(offset_t)teems_maps[cond_mapid-1].toset;
  for (l=0; l<(dim_t)sets[cs].size; l++) if (strcmp(set_elems[sets[cs].offset+l].setele,cond_rhs)==0) {
      *condfix=(offset_t)l;
      return;
    }
  printf("Error: the sum condition RHS %s for mapping %s is neither a quantifier index nor an element of the codomain set %s\n",cond_rhs,teems_maps[cond_mapid-1].mapname,sets[cs].setname);
  MPI_Abort(PETSC_COMM_WORLD,1);
}

/* The condition's value table is indexed by domain position, so the
   summed set must be the mapping's domain exactly (the M2 exact-match
   contract). */
void sum_cond_domain_check(sum_def *sc, set_def *sets) {
  if (sc->cond_mapid==0) return;
  if ((offset_t)teems_maps[sc->cond_mapid-1].fromset!=sc->sumsetid) {
    printf("Error: the condition on mapping %s sums over set %s, not the mapping's domain set %s\n",teems_maps[sc->cond_mapid-1].mapname,sets[sc->sumsetid].setname,sets[teems_maps[sc->cond_mapid-1].fromset].setname);
    MPI_Abort(PETSC_COMM_WORLD,1);
  }
}

/* A mapping-equality condition whose RHS is an outer quantifier index
   makes the sum's value depend on that index even when the body never
   references it (the GTAP-E shape: tot(b) = sum over the domain of
   the r with MAP(r)=b): append it as a carried dim so the sum store
   holds one slot per RHS element.  The RHS quantifier must range over
   the codomain exactly.  Returns the updated carried-dim count. */
dim_t sum_cond_carry_rhs(sum_def *sc, quantifier *arSet, dim_t fdim, dim_t l3, char *interchar) {
  dim_t l4,l5;
  if (sc->cond_mapid==0) return l3;
  for (l5=0; l5<fdim-1; l5++) if (strcmp(sc->cond_rhs,arSet[l5].index_name)==0) break;
  if (l5>=fdim-1) return l3; /* not a quantifier (scalar statements have fdim 0): codomain element, resolved at eval */
  if (arSet[l5].setid!=(offset_t)teems_maps[sc->cond_mapid-1].toset) {
    printf("Error: the RHS quantifier %s of the condition on mapping %s does not range over the mapping's codomain set\n",sc->cond_rhs,teems_maps[sc->cond_mapid-1].mapname);
    MPI_Abort(PETSC_COMM_WORLD,1);
    return l3;
  }
  for (l4=0; l4<l3; l4++) if (strcmp(sc->cond_rhs,sc->dimnames[l4])==0) return l3;
  strcpy(sc->dimnames[l3],sc->cond_rhs);
  sc->setid[l3]=arSet[l5].setid;
  strcat(interchar,sc->cond_rhs);
  return l3+1;
}

/* named fatal for statement kinds that cannot carry mapping calls yet
   -- without this the reference mis-parses in silence */
void mapping_reject_in(char *line, const char *what) {
  dim_t m;
  int k;
  char find[NAMESIZE+2];
  if (teems_nmap==0) return;
  for (m=0; m<teems_nmap; m++) {
    sprintf(find,"%s(",teems_maps[m].mapname);
    k=str_find_ci(line,find);
    if (k==0||(k>0&&!isalnum((int)line[k-1])&&line[k-1]!='_')) {
      printf("Error: mapping-valued indices in %s statements are not supported yet (mapping %s)\n",what,teems_maps[m].mapname);
      MPI_Abort(PETSC_COMM_WORLD,1);
    }
  }
}

/* Interim guard for the bordered matrix methods (design doc M2b): a
   mapped variable reference can place a column in another diagonal
   block invisibly to the SBBD/DBBD/NDBBD border classification, so a
   tabfile whose equations carry lowered mapping tokens runs only under
   -matsol 0 (LU) until the border marking learns mapped references.
   Scans the post-lowering tabfile; caller gates on the matrix method. */
void mapping_eq_matsol_guard(char *fname) {
  FILE *filehandle;
  char line[TABREADLINE]="\0";
  if (teems_nmap==0) return;
  filehandle = fopen(fname,"r");
  if (filehandle==NULL) return;
  while (tab_next_statement("equation",filehandle,line,TABREADLINE)) {
    if (strchr(line,'@')!=NULL) {
      printf("Error: mapping-valued indices in Equation statements need -matsol 0 (LU); the bordered methods' block classification does not follow mapped references yet\n");
      fclose(filehandle);
      MPI_Abort(PETSC_COMM_WORLD,1);
    }
  }
  fclose(filehandle);
}

/* Pre-use mapping validation (manual 11.9.2/11.9.3): every mapping
   must have values, and an (onto) mapping must cover its codomain.
   Range validity is guaranteed by construction in
   mapping_values_read (names resolve to positions or abort). */
int mappings_validate(map_def *maps, dim_t nmap, set_def *sets, set_element *set_elems) {
  dim_t j,i;
  for (j=0; j<nmap; j++) {
    if (!maps[j].has_values) {
      printf("Error: mapping %s is never given values (manual 11.9.1)\n",maps[j].mapname);
      MPI_Abort(PETSC_COMM_WORLD,1);
      return -1;
    }
    if (maps[j].onto) {
      dim_t tosize=sets[maps[j].toset].size;
      bool *hit= (bool *) calloc (tosize>0?tosize:1,sizeof(bool));
      for (i=0; i<sets[maps[j].fromset].size; i++) hit[maps[j].values[i]]=true;
      for (i=0; i<tosize; i++) if (!hit[i]) {
          printf("Error: mapping %s is not onto: element %s of set %s is not mapped to (manual 11.9.3)\n",maps[j].mapname,set_elems[sets[maps[j].toset].offset+i].setele,sets[maps[j].toset].setname);
          free(hit);
          MPI_Abort(PETSC_COMM_WORLD,1);
          return -1;
        }
      free(hit);
    }
  }
  return 0;
}

offset_t closure_read(char *fname, char *commsyntax,closure_entry *closure_vals, array_def *vars,offset_t nvar,set_def *sets,dim_t nset, set_element *set_elems) {
  FILE * filehandle;
  char line[TABREADLINE]="\0",*readitem=NULL,*p=NULL,*p1=NULL,vname[TABREADLINE],argu[TABREADLINE];//,linecopy[TABREADLINE]
  offset_t i,j,n=0,l=0,l1,l2,n1,m,dims,dcount,exodims,supsetid[MAXSUPSET],sup,doublepr[MAXSUPSET];
  int k1,k2;
  bool check;
  filehandle = fopen(fname,"r");
  if (filehandle==NULL) {
    printf("Error: cannot open closure file %s\n",fname);
    MPI_Abort(PETSC_COMM_WORLD,1);
    return -1;
  }
  while (closure_next_statement(commsyntax,filehandle,line)) {
    while (str_strip_comment(line,"!"));
    str_replace_all(line,";", " ;");
    while (str_replace_all(line,"\n", " "));
    while (str_replace_all(line,"\r", " "));
    while (str_replace_first(line," p_", " "));
    while (str_replace_first(line," c_", " "));
    while (str_replace_all(line,"  ", " "));
    while (str_replace_all(line," ,", ","));
    while (str_replace_all(line,", ", ","));
    while (str_replace_all(line," )", ")"));
    while (str_replace_all(line," (", "("));
    while (str_replace_all(line," \"", "\""));
    while (str_replace_all(line,"\" ", "\""));
    k1=0;
    k2=0;
    while (line[k1]!= '\0') {
      if(line[k1]=='\"') {
        if(k2==0) {
          k2=1;
        } else {
          k2=0;
        }
      } else {
        if(k2==0) {
          line[k1]=tolower((int)line[k1]);
        }
      }
      k1++;
    }
    n1=str_count_char(line,' ');
    readitem = strchr(line,' ');
    if (readitem==NULL) continue;
    readitem++;
    for (i=0; i<n1-1; i++) {
      p = strchr(readitem,' ');
      if (p==NULL) break;
      strncpy(vname,readitem,p-readitem);
      vname[p-readitem]='\0';
      readitem=p;
      readitem++;
      check=true;
      if (strchr(vname,'(')==NULL) {
        for (j=0; j<nvar; j++) if (strcmp(vname,vars[j].cofname)==0) {
            dims=vars[j].nelem;
            if (dims==0) {
              closure_vals[vars[j].offset].is_exogenous=true;
              n=n+1;
            } else for (l=0; l<dims; l++) {
                closure_vals[vars[j].offset+l].is_exogenous=true;
                n=n+1;
              }
            check=false;
            break;
          }
          if(j==nvar&&vname[0]!=';') {
            printf("Error: variable %s is not declared (closure file)\n",vname);
            fclose(filehandle);
            MPI_Abort(PETSC_COMM_WORLD,1);
            return -1;
          }
      } else {
        dims=1;
        p = strtok(vname,"(");
        for(sup=0; sup<MAXSUPSET; sup++)supsetid[sup]=0;
        for (j=0; j<nvar; j++) if (strcmp(p,vars[j].cofname)==0) {
            quantifier *arSet= (quantifier *) calloc (vars[j].size,sizeof(quantifier));
            offset_t *exoantidim= (offset_t *) calloc (vars[j].size,sizeof(offset_t));
            switch (vars[j].size) {
            case 1:
              p = strtok(NULL,")");//p = strtok(vname,")");
              if (p==NULL) {
                printf("Error: malformed entry for variable %s (closure file)\n",vars[j].cofname);
                free(arSet);
                free(exoantidim);
                fclose(filehandle);
                MPI_Abort(PETSC_COMM_WORLD,1);
                return -1;
              }
              if (strchr(p,'"')!=NULL) {
                strcpy(argu,p);
                p=strchr(argu,'"');
                p+=1;
                p1=strchr(p,'"');
                *p1='\0';
                for (l1=0; l1<sets[vars[j].setid[0]].size; l1++)
                  if (strcmp(p,set_elems[sets[vars[j].setid[0]].offset+l1].setele)==0) {
                    closure_vals[vars[j].offset+l1].is_exogenous=true;
                    n=n+1;
                    check=false;
                    break;
                  }
                  if(l1==sets[vars[j].setid[0]].size) {
                    printf("Error: element %s is not in set %s (in %s)\n",p,sets[vars[j].setid[0]].setname,vars[j].cofname);
                    free(arSet);
                    free(exoantidim);
                    fclose(filehandle);
                    MPI_Abort(PETSC_COMM_WORLD,1);
                    return -1;
                  }
              } else {
                if(strcmp(p,sets[vars[j].setid[0]].setname)==0) {
                  arSet[0].setid=vars[j].setid[0];
                  supsetid[0]=0;
                  exoantidim[0]=1;//arSet[0].SetSize;
                  dims=sets[arSet[0].setid].size;
                } else {
                  for (l1=0; l1<nset; l1++) if (strcmp(p,sets[l1].setname)==0) {
                      arSet[0].setid=l1;
                      for(sup=1; sup<MAXSUPSET; sup++)if(sets[l1].subsetid[sup]==vars[j].setid[0])supsetid[0]=sup;
                      exoantidim[0]=1;//arSet[0].SetSize;
                      dims=sets[arSet[0].setid].size;
                      break;
                    }
                    if(l1==nset) {
                      printf("Error: set %s is not declared (in %s)\n",p,vars[j].cofname);
                      free(arSet);
                      free(exoantidim);
                      fclose(filehandle);
                      MPI_Abort(PETSC_COMM_WORLD,1);
                      return -1;
                    }
                }
              }
              break;
            default:
              p = strtok(NULL,")");
              if (p==NULL) {
                printf("Error: malformed entry for variable %s (closure file)\n",vars[j].cofname);
                free(arSet);
                free(exoantidim);
                fclose(filehandle);
                MPI_Abort(PETSC_COMM_WORLD,1);
                return -1;
              }
              strcat(p,",");
              for (l=0; l<vars[j].size; l++) {
                if (l==0) {
                  p = strtok(p,",");
                } else {
                  p = strtok(NULL,",");
                }
                if (p==NULL) {
                  printf("Error: wrong number of arguments for variable %s (closure file)\n",vars[j].cofname);
                  free(arSet);
                  free(exoantidim);
                  fclose(filehandle);
                  MPI_Abort(PETSC_COMM_WORLD,1);
                  return -1;
                }
                if (*p=='"') {
                  strcpy(argu,p+1);
                  str_replace_char(argu,'"','\0');
                  for(l1=0; l1<sets[vars[j].setid[l]].size; l1++) if (strcmp(argu,set_elems[sets[vars[j].setid[l]].offset+l1].setele)==0) {
                      doublepr[l]=l1;
                      supsetid[l]=-2;
                      break;
                    }
                    if(l1==sets[vars[j].setid[l]].size) {
                      printf("Error: element %s is not in set %s (in %s)\n",argu,sets[vars[j].setid[l]].setname,vars[j].cofname);
                      free(arSet);
                      free(exoantidim);
                      fclose(filehandle);
                      MPI_Abort(PETSC_COMM_WORLD,1);
                      return -1;
                    }
                } else {
                  if(strcmp(p,sets[vars[j].setid[l]].setname)==0) {
                    arSet[l].setid=vars[j].setid[l];
                    dims=dims*sets[arSet[l].setid].size;
                  } else {
                    for(l1=0; l1<nset; l1++) if (strcmp(p,sets[l1].setname)==0) {
                        arSet[l].setid=l1;
                        for(sup=1; sup<MAXSUPSET; sup++)if(sets[l1].subsetid[sup]==vars[j].setid[l])supsetid[l]=sup;
                        dims=dims*sets[arSet[l].setid].size;
                        break;
                      }
                      if(l1==nset) {
                        printf("Error: set %s is not declared (in %s)\n",p,vars[j].cofname);
                        free(arSet);
                        free(exoantidim);
                        fclose(filehandle);
                        MPI_Abort(PETSC_COMM_WORLD,1);
                        return -1;
                      }
                  }
                }
              }
              exoantidim[vars[j].size-1]=1;
              for (l1=vars[j].size-2; l1>-1; l1--) {
                if(supsetid[l1+1]>-2)exoantidim[l1]=exoantidim[l1+1]*sets[arSet[l1+1].setid].size;
                else exoantidim[l1]=exoantidim[l1+1];
              }
              break;
            }
            if(check) for (l=0; l<dims; l++) {
                m=l;
                exodims=0;
                for (dcount=0; dcount<vars[j].size; dcount++) {
                  l1=(offset_t) m/exoantidim[dcount];
                  m=m-l1*exoantidim[dcount];
                  if(supsetid[dcount]==-2) {
                    l2=doublepr[dcount];//ha_set[arSet[dcount].setid].begadd;  //Special case
                  } else {
                    l2=set_elems[sets[arSet[dcount].setid].offset+l1].superset_pos[supsetid[dcount]];
                  }
                  exodims=exodims+l2*vars[j].strides[dcount];
                }
                closure_vals[vars[j].offset+exodims].is_exogenous=true;
                n=n+1;
              }
            free(arSet);
            free(exoantidim);
            break;
          }
          if(j==nvar) {
            printf("Error: variable %s is not declared (closure file)\n",p);
            fclose(filehandle);
            MPI_Abort(PETSC_COMM_WORLD,1);
            return -1;
          }
      }
    }
  }
  fclose(filehandle);
  return n;
}

offset_t shocks_read(char *fname, char *commsyntax,closure_entry *closure_vals,offset_t nvarele, array_def *vars,offset_t nvar,set_def *sets,dim_t nset, set_element *set_elems,dim_t subints) {
  FILE * filehandle;
  char line[DATREADLINE],linecopy[DATREADLINE],argu[TABREADLINE];
  line[0]='\0';
  linecopy[0]='\0';
  argu[0]='\0';
  char *readitem=NULL,*p=NULL;
  dim_t varnset;
  dim_t k1,k2,dimindx[MAXVARDIM],dimbegadd[MAXVARDIM],dimssize[MAXVARDIM],arSet[MAXVARDIM],antidim[MAXVARDIM];
  offset_t j,l=0,dims,n1,l1,l2,dcount,supsetid[MAXSUPSET],sup;
  solve_real val;
      if ( (filehandle = fopen(fname,"r")) == NULL ) {
        printf("Error: cannot open shock file %s\n",fname);
        MPI_Abort(PETSC_COMM_WORLD,1);
        return -1;
      }
  while (tab_next_statement(commsyntax,filehandle,line,DATREADLINE)) {
    str_replace_char_all(line,'\r',' ');
    str_replace_char_all(line,'\n',' ');
    k1=0;
    k2=0;
    while (line[k1]!= '\0') {
      line[k1]=tolower((int)line[k1]);
      k1++;
    }
    if(k1>=DATREADLINE){
      printf("Error: shock statement exceeds the line buffer; increase DATREADLINE\n");
      fclose(filehandle);
      MPI_Abort(PETSC_COMM_WORLD,1);
      return -1;
    }
    str_replace_all_bounded(line,";", " ;",DATREADLINE);
    while (str_replace_all(line,"= ", "="));
    while (str_replace_all_bounded(line,"  ", " ",DATREADLINE));
    while (str_replace_first_bounded(line," p_", " ",DATREADLINE));
    while (str_replace_first_bounded(line," c_", " ",DATREADLINE));
    varnset=str_count_char(line,',');
    j=str_count_char(line,'(');
    for(sup=0; sup<MAXSUPSET; sup++)supsetid[sup]=-1;
    if (j>0) {
      varnset=varnset+1;
    }
    if(varnset==0){
      while(str_replace_all(line," =", "="));
    }

    if (line[0]==' ') {
      readitem = strtok(line+1," ");
    } else {
      readitem = strtok(line," ");
    }
    if (varnset==0) {
      readitem = strtok(NULL,"=");
      if (readitem==NULL) {
        printf("Error: malformed shock statement (shock file)\n");
        fclose(filehandle);
        MPI_Abort(PETSC_COMM_WORLD,1);
        return -1;
      }
      for (j=0; j<nvar; j++) if (strcmp(readitem,vars[j].cofname)==0) {
          readitem = strtok(NULL,";");
          if (readitem==NULL) {
            printf("Error: shock statement for variable %s has no value (shock file)\n",vars[j].cofname);
            fclose(filehandle);
            MPI_Abort(PETSC_COMM_WORLD,1);
            return -1;
          }
          k1=str_rfind_ci(readitem,"uniform");
          if(k1!=-1){
            closure_vals[vars[j].offset].shock_value=atof(readitem+k1+1)/subints;
          }else{
            closure_vals[vars[j].offset].shock_value=atof(readitem)/subints;
          }
          l=l+1;
          break;
        }
      /* previously fell through in silence: the shock never landed and
         the run solved an unshocked model */
      if (j==nvar&&readitem[0]!=';') {
        printf("Error: %s in the shock file is not a declared variable\n",readitem);
        fclose(filehandle);
        MPI_Abort(PETSC_COMM_WORLD,1);
        return -1;
      }
    } else {
      readitem = strtok(NULL,"(");
      for (j=0; j<nvar; j++) {
        if (strcmp(readitem,vars[j].cofname)==0) {
          dims=1;
          readitem = strtok(NULL,")");
          if (readitem==NULL) {
            printf("Error: malformed shock statement for variable %s (shock file)\n",vars[j].cofname);
            fclose(filehandle);
            MPI_Abort(PETSC_COMM_WORLD,1);
            return -1;
          }
          strcpy(argu,readitem);
          readitem = strtok(NULL,"=");
          readitem = strtok(NULL,";");
          if (readitem==NULL) {
            printf("Error: shock statement for variable %s has no value (shock file)\n",vars[j].cofname);
            fclose(filehandle);
            MPI_Abort(PETSC_COMM_WORLD,1);
            return -1;
          }
          strcpy(linecopy,readitem);
          strcat(linecopy," ");
          strcat(argu,",");
          while (str_replace_all(argu," ", ""));
          for (n1=0; n1<MAXVARDIM; n1++) {
            dimindx[n1]=-1;
            dimbegadd[n1]=-1;
            dimssize[n1]=0;
            arSet[n1]=0;
            antidim[n1]=0;
          }
          for (n1=0; n1<vars[j].size; n1++) {
            if(n1==0) {
              p=strtok(argu,",");
            } else {
              p=strtok(NULL,",");
            }
            if (p==NULL) {
              printf("Error: wrong number of arguments for variable %s (shock file)\n",vars[j].cofname);
              fclose(filehandle);
              MPI_Abort(PETSC_COMM_WORLD,1);
              return -1;
            }
            if(strchr(p,'\"')==NULL) {
              if((strcmp(sets[vars[j].setid[n1]].setname,p)==0)) {
                dimbegadd[n1]=sets[vars[j].setid[n1]].offset;//ha_set[k1].begadd;
                dimssize[n1]=sets[vars[j].setid[n1]].size;//ha_set[k1].size;
                dims=dims*dimssize[n1];
                supsetid[n1]=0;
              } else {
                for (k1=0; k1<nset; k1++)if (strcmp(sets[k1].setname,p)==0) {
                    dimbegadd[n1]=sets[k1].offset;
                    dimssize[n1]=sets[k1].size;
                    dims=dims*sets[k1].size;
                    for(sup=0; sup<MAXSUPSET; sup++)if(sets[k1].subsetid[sup]==vars[j].setid[n1]) {
                        supsetid[n1]=sup;
                        break;
                      }
                    if(sup==MAXSUPSET){
                      printf("Error: %s is not a valid set at that position of variable %s (shock file)\n",p,vars[j].cofname);
                      MPI_Abort(PETSC_COMM_WORLD,1);
                      return -1;
                    }
                    break;
                  }
                  if(k1==nset){
                    printf("Error: %s in variable %s (shock file) is not a set\n",p,vars[j].cofname);
                    MPI_Abort(PETSC_COMM_WORLD,1);
                    return -1;
                  }
              }
            } else {
              k2=0;
              while (p[k2] != '\0') {
                k2++;
              }
              for (k1=0; k1<sets[vars[j].setid[n1]].size; k1++)if (strncmp(set_elems[sets[vars[j].setid[n1]].offset+k1].setele,p+1,k2-2)==0&&set_elems[sets[vars[j].setid[n1]].offset+k1].setele[k2-2]=='\0') {
                  dimindx[n1]=set_elems[sets[vars[j].setid[n1]].offset+k1].superset_pos[0];
                  dimssize[n1]=1;
                  supsetid[n1]=-2;
                  break;
                }
              if(k1==sets[vars[j].setid[n1]].size){
                printf("Error: %s in variable %s (shock file) is not a set element\n",p,vars[j].cofname);
                MPI_Abort(PETSC_COMM_WORLD,1);
                return -1;
              }
            }
          }
          /* GEMPACK-standard value ordering: first subscript varies
             fastest (column-major), fixed dims collapsed to size 1 */
          if (vars[j].size>0)antidim[0]=1;
          for (n1=1; n1<vars[j].size; n1++) {
            antidim[n1]=antidim[n1-1]*dimssize[n1-1];
          }
          k1=str_rfind_ci(linecopy,"uniform");
          if(k1>-1) {
            val=atof(linecopy+k1+1);//sscanf(readitem, "%lf", &val);
            for (n1=0; n1<dims; n1++) {
              l2=n1;
              if(vars[j].size>1) {
                for (dcount=vars[j].size-1; dcount>=0; dcount--) {
                  if(dimindx[dcount]>-1) {
                    arSet[dcount]=dimindx[dcount];
                  } else {
                    l1=(offset_t) l2/antidim[dcount];
                    if(supsetid[dcount]>-2)arSet[dcount]=set_elems[dimbegadd[dcount]+l1].superset_pos[supsetid[dcount]];
                    else arSet[dcount]=set_elems[dimbegadd[dcount]+l1].superset_pos[0];
                    l2=l2-l1*antidim[dcount];
                  }
                }
              } else {
                arSet[0]=l2;
              }
              l2=0;
              for (dcount=0; dcount<vars[j].size; dcount++) {
                l2=l2+arSet[dcount]*vars[j].strides[dcount];
              }
              closure_vals[vars[j].offset+l2].shock_value=val/subints;
              l=l+1;
            }
          } else {
            for (n1=0; n1<dims; n1++) {
              if(n1==0) {
                if(linecopy[0]==' ') {
                  readitem=strtok(linecopy+1," ");
                } else {
                  readitem=strtok(linecopy," ");
                }
              } else {
                readitem=strtok(NULL," ");
              }
              if (readitem==NULL) {
                printf("Error: shock statement for variable %s supplies fewer values than elements (%ld expected) (shock file)\n",vars[j].cofname,(long)dims);
                fclose(filehandle);
                MPI_Abort(PETSC_COMM_WORLD,1);
                return -1;
              }
              val=atof(readitem);
              l2=n1;
              if(vars[j].size>1) {
                for (dcount=vars[j].size-1; dcount>=0; dcount--) {
                  if(dimindx[dcount]>-1) {
                    arSet[dcount]=dimindx[dcount];
                  } else {
                    l1=(offset_t) l2/antidim[dcount];
                    if(supsetid[dcount]>-2)arSet[dcount]=set_elems[dimbegadd[dcount]+l1].superset_pos[supsetid[dcount]];
                    else arSet[dcount]=set_elems[dimbegadd[dcount]+l1].superset_pos[0];
                    l2=l2-l1*antidim[dcount];
                  }
                }
              } else {
                if(supsetid[0]>-2)arSet[0]=set_elems[dimbegadd[0]+l2].superset_pos[supsetid[0]];
                else arSet[0]=dimindx[0];
              }
              l2=0;
              for (dcount=0; dcount<vars[j].size; dcount++) {
                l2=l2+arSet[dcount]*vars[j].strides[dcount];
              }
              closure_vals[vars[j].offset+l2].shock_value=val/subints;
              l=l+1;
            }
          }
          break;
        }
      }
      if(j==nvar){
        printf("Error: %s in the shock file is not a declared variable\n",readitem);
        MPI_Abort(PETSC_COMM_WORLD,1);
        return -1;
      }
    }
  }
  fclose(filehandle);
  return l;
}

/* Default statements (manual 10.19): "Keyword (DEFAULT = value) ;" --
   DEFAULT is a qualifier, not a keyword; the statement resets the
   named default for the declarations that FOLLOW it (positional).
   Detects one and extracts the space-squeezed value ("initial",
   "add_homotopy=h2", ...) into out (NAMESIZE); returns 1 for a
   default statement, 0 otherwise. Value validation is centralized in
   tab_defaults_validate (audit A6). */
int tab_default_value(char *line, char *out) {
  char *p=strstr(line,"(default");
  int k=0;
  out[0]='\0';
  if(p==NULL)return 0;
  p+=8;
  while(*p==' ')p++;
  if(*p=='=')p++;
  while(*p!=')'&&*p!=';'&&*p!='\0'&&k<NAMESIZE-1) {
    if(*p!=' ')out[k++]=*p;
    p++;
  }
  out[k]='\0';
  return 1;
}

/* Tokenized qualifier-list parsing for VARIABLE and COEFFICIENT
   declarations (manual 10.3/10.4; audit A4, closes batch-1 residual
   (b)). Strips every LEADING (qualifier,...) group from the statement
   -- the first (all,...) quantifier ends the qualifier region -- and
   applies the recognized semantics to the record. The old enumerated
   substring matching let any unlisted qualifier or combination fall
   through into the quantifier walk and silently corrupt the
   declaration; unknown tokens are fatal instead. Semantics-bearing
   qualifiers we do not implement (LINEAR_NAME=/LINEAR_VAR= renames the
   generated linear variable, NO_SPLIT changes multi-step shock
   treatment) are fatal until their features land; reporting/metadata
   qualifiers (ORIG_LEVEL=, VPQTYPE=) parse and are ignored. Caller
   gates Default statements out (legacy sticky handling, audit A6).
   is_variable selects the legal token set; explicit_param
   (coefficients only) reports an explicit PARAMETER/NON_PARAMETER
   token (1/0, -1 = none) and isint_out an INTEGER token -- the
   default resolution (10.3/10.19) is the caller's. A declaration may
   carry one lower (GE/GT) and one upper (LE/LT) bound (manual
   10.19.1; audit A9): the first goes to rec->gltype/glval, the
   opposite-direction second to gltype2/glval2; a same-direction
   duplicate is fatal. Returns 1 on success, -1 on error. */
static int tab_qualifiers_parse(char *line, offset_t ncommsyntax, array_def *rec, int *explicit_param, bool *isint_out, int *gltype2, store_real *glval2, int is_variable) {
  char group[TABREADLINE],tok[TABREADLINE],*p,*q,*t,*s;
  offset_t len;
  int k,nbtype;
  p=line+ncommsyntax;
  while(*p==' ')p++;
  while(*p=='(') {
    if(strncmp(p+1,"all,",4)==0||strncmp(p+1,"all ",4)==0)break;
    q=strchr(p,')');
    if(q==NULL) {
      printf("Error: unbalanced parentheses in %s qualifier list: %s\n",is_variable?"variable":"coefficient",line);
      return -1;
    }
    len=q-p-1;
    if(len>=(offset_t)sizeof(group))len=sizeof(group)-1;
    strncpy(group,p+1,len);
    group[len]='\0';
    /* excise the group and any blanks that followed it */
    t=q+1;
    while(*t==' ')t++;
    memmove(p,t,strlen(t)+1);
    for(t=strtok(group,","); t!=NULL; t=strtok(NULL,",")) {
      k=0;
      for(s=t; *s!='\0'; s++)if(*s!=' ')tok[k++]=*s;
      tok[k]='\0';
      if(k==0) {
        printf("Error: empty qualifier in %s declaration\n",is_variable?"variable":"coefficient");
        return -1;
      }
      if(is_variable) {
        if(strcmp(tok,"change")==0) {
          rec->change_real=true;
          continue;
        }
        if(strcmp(tok,"percent_change")==0) {
          rec->change_real=false;
          continue;
        }
        if(strcmp(tok,"linear")==0) {
          rec->level_par=false;
          continue;
        }
        if(strcmp(tok,"levels")==0) {
          rec->level_par=true;
          continue;
        }
        /* reporting-only: levels values of a linear variable (11.6.5) */
        if(strncmp(tok,"orig_level=",11)==0)continue;
        /* homogeneity-test metadata only (ch.57): no solution impact */
        if(strncmp(tok,"vpqtype=",8)==0)continue;
        if(strcmp(tok,"no_split")==0) {
          printf("Error: variable qualifier NO_SPLIT (full shock at every step) is not supported\n");
          return -1;
        }
        if(strncmp(tok,"linear_name=",12)==0||strncmp(tok,"linear_var=",11)==0) {
          printf("Error: variable qualifier LINEAR_NAME=/LINEAR_VAR= is not supported yet; use the default p_/c_ linear name\n");
          return -1;
        }
      } else {
        if(strcmp(tok,"real")==0)continue;
        if(strcmp(tok,"integer")==0) {
          *isint_out=true;
          continue;
        }
        if(strcmp(tok,"parameter")==0) {
          *explicit_param=1;
          continue;
        }
        if(strcmp(tok,"non_parameter")==0) {
          *explicit_param=0;
          continue;
        }
      }
      nbtype=0;
      if(tok[0]=='g'&&tok[1]=='e')nbtype=BT_GE;
      if(tok[0]=='g'&&tok[1]=='t')nbtype=BT_GT;
      if(tok[0]=='l'&&tok[1]=='e')nbtype=BT_LE;
      if(tok[0]=='l'&&tok[1]=='t')nbtype=BT_LT;
      if(nbtype>0&&(isdigit((int)tok[2])||tok[2]=='-'||tok[2]=='+'||tok[2]=='.')) {
        if(rec->gltype>0) {
          if((rec->gltype==BT_GE||rec->gltype==BT_GT)==(nbtype==BT_GE||nbtype==BT_GT)) {
            printf("Error: duplicate %s bound on a %s declaration (one lower GE/GT and one upper LE/LT allowed)\n",(nbtype==BT_GE||nbtype==BT_GT)?"lower":"upper",is_variable?"variable":"coefficient");
            return -1;
          }
          *gltype2=nbtype;
          *glval2=atof(tok+2);
        } else {
          rec->gltype=nbtype;
          rec->glval=atof(tok+2);
        }
        continue;
      }
      printf("Error: unknown %s qualifier '%s'\n",is_variable?"variable":"coefficient",tok);
      return -1;
    }
    while(*p==' ')p++;
  }
  return 1;
}

offset_t variables_read(char *fname, char *commsyntax, array_def *record, offset_t ncof, set_def *sets,dim_t nset) {
  FILE * filehandle;
  char line[TABREADLINE]="\0",linecopy[TABREADLINE],setname1[NAMESIZE],setname[NAMESIZE],setname2[NAMESIZE],setname3[NAMESIZE],finditem[NAMESIZE],finditem1[NAMESIZE],finditem2[NAMESIZE],finditem3[NAMESIZE],vname[NAMESIZE],defval[NAMESIZE];//,vnamecopy[NAMESIZE];
  int vgltype2=0;
  store_real vglval2=0;
  offset_t n,m,l,ncommsyntax=0,i,j=0,addi=0,add=0;
  dim_t dcount;
  char *readitem=NULL;
  bool IsLevel=false,IsChange=false;
  while (commsyntax[ncommsyntax] != '\0') {
    ncommsyntax++;
  }
  filehandle = fopen(fname,"r");

  while (tab_next_statement(commsyntax,filehandle,line,TABREADLINE)) {
    /* positional defaults (manual 10.19; audit A6): reset the LINEAR/
       LEVELS and PERCENT_CHANGE/CHANGE defaults for the declarations
       that follow; both directions and spaced forms. Values are
       validated up front by tab_defaults_validate. */
    if(tab_default_value(line,defval)) {
      if(strcmp(defval,"levels")==0)IsLevel=true;
      else if(strcmp(defval,"linear")==0)IsLevel=false;
      else if(strcmp(defval,"change")==0)IsChange=true;
      else if(strcmp(defval,"percent_change")==0)IsChange=false;
      continue;
    }
    if(IsLevel) {
      record[j].level_par=true;
    }
    if(IsChange) {
      record[j].change_real=true;
    }
    /* tokenized qualifier parsing (audit A4): explicit qualifiers
       override the defaults above */
    if (tab_qualifiers_parse(line,ncommsyntax,&record[j],NULL,NULL,&vgltype2,&vglval2,1)<0) {
      fclose(filehandle);
      return -1;
    }
    while (str_replace_all(line," ", ""));
    strcpy(linecopy,line);
    {
      n=str_count_char(line,')');
      if (n>1) {
        for (i=1; i<n; i++) {
          if (i==1) {
            readitem = strtok(line,",");
          } else {
            readitem = strtok(NULL,",");
          }
          readitem = strtok(NULL,",");
          if (readitem==NULL||strlen(readitem)+3>sizeof(finditem)) {
            printf("Error: malformed %s declaration in TAB file: %s\n",commsyntax,linecopy);
            return -1;
          }
          finditem[0]='(';
          finditem[1]='\0';
          finditem1[0]=',';
          finditem1[1]='\0';
          finditem2[0]=',';
          finditem2[1]='\0';
          strcat(finditem,readitem);
          strcat(finditem1,readitem);
          strcat(finditem2,readitem);
          strcpy(finditem3,finditem);
          strcat(finditem,",");
          strcat(finditem1,",");
          strcat(finditem2,")");
          strcat(finditem3,")");
          readitem = strtok(NULL,")");
          if (readitem==NULL||strlen(readitem)+3>sizeof(setname)) {
            printf("Error: malformed %s declaration in TAB file: %s\n",commsyntax,linecopy);
            return -1;
          }
          setname[0]='(';
          setname[1]='\0';
          setname1[0]=',';
          setname1[1]='\0';
          setname2[0]=',';
          setname2[1]='\0';
          strcat(setname,readitem);
          strcat(setname1,readitem);
          strcat(setname2,readitem);
          strcpy(setname3,setname);
          strcat(setname,",");
          strcat(setname1,",");
          strcat(setname2,")");
          strcat(setname3,")");
          if (str_subst_all_bounded(linecopy, finditem,  setname,  sizeof(linecopy)) ||
              str_subst_all_bounded(linecopy, finditem1, setname1, sizeof(linecopy)) ||
              str_subst_all_bounded(linecopy, finditem2, setname2, sizeof(linecopy)) ||
              str_subst_all_bounded(linecopy, finditem3, setname3, sizeof(linecopy))) {
            printf("Error: malformed %s declaration in TAB file: %s\n",commsyntax,linecopy);
            return -1;
          }
        }
        for (i=1; i<n; i++) {
          if (i==1) {
            readitem = strtok(linecopy,")");
          } else {
            readitem = strtok(NULL,")");
          }
        }
        readitem = strtok(NULL,";");
        if (readitem!=NULL) readitem = strtok(readitem,"(");
        if (readitem==NULL||strlen(readitem)>=sizeof(record[j].cofname)) {
          printf("Error: malformed %s declaration in TAB file\n",commsyntax);
          return -1;
        }
        strcpy(record[j].cofname,readitem);
        readitem = strtok(NULL,")");
        if (readitem==NULL||strlen(readitem)+2>sizeof(vname)) {
          printf("Error: malformed %s declaration in TAB file\n",commsyntax);
          return -1;
        }
        strcpy(vname,readitem);
        strcat(vname,",");
        dcount=str_count_char(vname,',');
        add=1;
        for (m=0; m<dcount; m++) {
          if(m==0) {
            readitem = strtok(vname,",");
          } else {
            readitem = strtok(NULL,",");
          }
          if (readitem==NULL) {
            printf("Error: malformed %s declaration in TAB file\n",commsyntax);
            return -1;
          }
          for (l=0; l<nset; l++) {
            if (strcmp(sets[l].setname,readitem)==0) {
              add=add*sets[l].size;
              record[j].setid[m]=l;
              break;
            }
          }
        }
        record[j].offset=addi;
        record[j].size=dcount;
        record[j].nelem=add;
        record[j].strides[dcount-1]=1;
        for (l=record[j].size-2; l>-1; l--) {
          record[j].strides[l]=record[j].strides[l+1]*sets[record[j].setid[l+1]].size;
        }
        addi=addi+add;
      } else {
        if (n==1) {
          printf("Error: unbalanced parentheses in statement\n");
        } else {
          strcpy(setname,line+ncommsyntax);
          str_replace_all(setname,";", "");
          str_replace_all(setname,"\n", "");
          if (strlen(setname)>=sizeof(record[j].cofname)) {
            printf("Error: malformed %s declaration in TAB file\n",commsyntax);
            return -1;
          }
          strcpy(record[j].cofname,setname);
          record[j].offset=addi;
          record[j].size=0;
          record[j].nelem=1;
          addi=addi+1;
        }
      }
      j++;
    }
  }
  fclose(filehandle);
  return addi;
}

offset_t set_find_alltime(set_def *sets,dim_t nset) {
  dim_t i;
  for(i=0; i<nset; i++) {
    if(sets[i].intertemp&&sets[i].subsetid[1]==-1) {
      break;
    }
  }
  if(i==nset) {
    return -1;
  } else {
    return i;
  }
}


offset_t coefficients_read(char *fname, char *commsyntax, array_def *record, offset_t ncof, set_def *sets,offset_t nset) {
  FILE * filehandle;
  char line[TABREADLINE]="\0",linecopy[TABREADLINE],setname1[TABREADLINE],setname[TABREADLINE],setname2[TABREADLINE],setname3[TABREADLINE],finditem[TABREADLINE],finditem1[TABREADLINE],finditem2[TABREADLINE],finditem3[TABREADLINE],vname[TABREADLINE];//,vnamecopy[NAMESIZE];
  offset_t n,m,l,ncommsyntax=0,i=0,j=0,addi=0,add=0;
  dim_t dcount;
  bool DefParam=false,isint=false;
  int expparam=-1,cgltype2=0;
  store_real cglval2=0;
  char defval[NAMESIZE];
  char *readitem=NULL;
  while (commsyntax[ncommsyntax] != '\0') {
    ncommsyntax++;
  }
  filehandle = fopen(fname,"r");
  /* PostSim foundation F2: (parameter) tracked in a parallel array --
     array_def is binary-locked to the R-side sol.var parser */
  free(teems_coef_is_param);
  teems_coef_is_param= (bool *) calloc (ncof+1,sizeof(bool));
  free(teems_coef_gltype2);
  free(teems_coef_glval2);
  teems_coef_gltype2= (int *) calloc (ncof+1,sizeof(int));
  teems_coef_glval2= (store_real *) calloc (ncof+1,sizeof(store_real));
  free(teems_coef_ifhdr);
  teems_coef_ifhdr= (bool *) calloc (ncof+1,sizeof(bool));

  while (tab_next_statement(commsyntax,filehandle,line,TABREADLINE)) {
    /* positional PARAMETER/NON_PARAMETER default for real
       coefficients (manual 10.19; audit A6); values validated up
       front by tab_defaults_validate */
    if(tab_default_value(line,defval)) {
      if(strcmp(defval,"parameter")==0)DefParam=true;
      else if(strcmp(defval,"non_parameter")==0)DefParam=false;
      continue;
    }
    /* tokenized qualifier parsing (audit A4; closes batch-1 residual
       (b)): an explicit PARAMETER/NON_PARAMETER token wins, INTEGER
       defaults to PARAMETER (10.3), real coefficients follow the
       positional default */
    expparam=-1;
    isint=false;
    cgltype2=0;
    cglval2=0;
    if (tab_qualifiers_parse(line,ncommsyntax,&record[j],&expparam,&isint,&cgltype2,&cglval2,0)<0) {
      fclose(filehandle);
      return -1;
    }
    if (j<ncof) {
      if (expparam>=0?expparam==1:(isint||DefParam))teems_coef_is_param[j]=true;
      teems_coef_gltype2[j]=cgltype2;
      teems_coef_glval2[j]=cglval2;
    }
    while (str_replace_all(line," ", ""));
    strcpy(linecopy,line);
    {
      n=str_count_char(line,')');
      if (n>1) {
        for (i=1; i<n; i++) {
          if (i==1) {
            readitem = strtok(line,",");
          } else {
            readitem = strtok(NULL,",");
          }
          readitem = strtok(NULL,",");
          if (readitem==NULL||strlen(readitem)+3>sizeof(finditem)) {
            printf("Error: malformed %s declaration in TAB file: %s\n",commsyntax,linecopy);
            return -1;
          }
          finditem[0]='(';
          finditem[1]='\0';
          finditem1[0]=',';
          finditem1[1]='\0';
          finditem2[0]=',';
          finditem2[1]='\0';
          strcat(finditem,readitem);
          strcat(finditem1,readitem);
          strcat(finditem2,readitem);
          strcpy(finditem3,finditem);
          strcat(finditem,",");
          strcat(finditem1,",");
          strcat(finditem2,")");
          strcat(finditem3,")");
          readitem = strtok(NULL,")");
          if (readitem==NULL||strlen(readitem)+3>sizeof(setname)) {
            printf("Error: malformed %s declaration in TAB file: %s\n",commsyntax,linecopy);
            return -1;
          }
          setname[0]='(';
          setname[1]='\0';
          setname1[0]=',';
          setname1[1]='\0';
          setname2[0]=',';
          setname2[1]='\0';
          strcat(setname,readitem);
          strcat(setname1,readitem);
          strcat(setname2,readitem);
          strcpy(setname3,setname);
          strcat(setname,",");
          strcat(setname1,",");
          strcat(setname2,")");
          strcat(setname3,")");
          if (str_subst_all_bounded(linecopy, finditem,  setname,  sizeof(linecopy)) ||
              str_subst_all_bounded(linecopy, finditem1, setname1, sizeof(linecopy)) ||
              str_subst_all_bounded(linecopy, finditem2, setname2, sizeof(linecopy)) ||
              str_subst_all_bounded(linecopy, finditem3, setname3, sizeof(linecopy))) {
            printf("Error: malformed %s declaration in TAB file: %s\n",commsyntax,linecopy);
            return -1;
          }
        }
        for (i=1; i<n; i++) {
          if (i==1) {
            readitem = strtok(linecopy,")");
          } else {
            readitem = strtok(NULL,")");
          }
        }
        readitem = strtok(NULL,";");
        if (readitem!=NULL) readitem = strtok(readitem,"(");
        if (readitem==NULL||strlen(readitem)>=sizeof(record[j].cofname)) {
          printf("Error: malformed %s declaration in TAB file\n",commsyntax);
          return -1;
        }
        strcpy(record[j].cofname,readitem);
        if(record[j].cofname[0]=='c'&&record[j].cofname[1]=='_'){
          printf("Error: the c_/C_ prefix is reserved for change variables; rename coefficient %s\n",record[j].cofname);
          return -1;
        }
        readitem = strtok(NULL,")");
        if (readitem==NULL||strlen(readitem)+2>sizeof(vname)) {
          printf("Error: malformed %s declaration in TAB file\n",commsyntax);
          return -1;
        }
        strcpy(vname,readitem);
        strcat(vname,",");
        dcount=str_count_char(vname,',');
        add=1;
        for (m=0; m<dcount; m++) {
          if(m==0) {
            readitem = strtok(vname,",");
          } else {
            readitem = strtok(NULL,",");
          }
          if (readitem==NULL) {
            printf("Error: malformed %s declaration in TAB file\n",commsyntax);
            return -1;
          }
          for (l=0; l<nset; l++) {
            if (strcmp(sets[l].setname,readitem)==0) {
              add=add*sets[l].size;
              record[j].setid[m]=l;
              break;
            }
          }
        }
        record[j].offset=addi;
        record[j].size=dcount;
        record[j].nelem=add;
        record[j].strides[dcount-1]=1;
        for (l=record[j].size-2; l>-1; l--) {
          record[j].strides[l]=record[j].strides[l+1]*sets[record[j].setid[l+1]].size;
        }
        addi=addi+add;
      } else {
        if (n==1) {
          printf("Error: unbalanced parentheses in statement: %s\n",line);
        } else {
          strcpy(setname,line+ncommsyntax);
          str_replace_all(setname,";", "");
          str_replace_all(setname,"\n", "");
          if (strlen(setname)>=sizeof(record[j].cofname)) {
            printf("Error: malformed %s declaration in TAB file\n",commsyntax);
            return -1;
          }
          strcpy(record[j].cofname,setname);
          record[j].offset=addi;
          record[j].size=0;
          record[j].nelem=1;
          addi=addi+1;
        }
      }
      j++;
    }
  }
  fclose(filehandle);

  return addi;
}


int coef_resolve_sets(array_def *coefs,offset_t ncof, set_def *sets,dim_t nset, elem_store *coef_store) {
  dim_t l,m,dcount;//,dcount1;//,l3=0,l4=0,ld2=0,ld3=0,ld4=0
  offset_t l1=0,i;//,sizeele,dcountdim1[MAXVARDIM],l2=0,dcountdim[MAXVARDIM];//,ld1=0
  for (i=0; i<ncof; i++) {
    for (dcount=0; dcount<coefs[i].size; dcount++) {
      for (l=0; l<nset; l++) {
        if (strcmp(sets[l].setname,sets[coefs[i].setid[dcount]].setname)==0) {
          sets[coefs[i].setid[dcount]].offset=sets[l].offset;
          break;
        }
      }
    }
  }
  return 1;
}

dim_t sets_count(char *fname) {
  FILE * filehandle;
  char line[TABREADLINE]="\0";
  char *commsyntax="set";
  dim_t j=0;
  filehandle = fopen(fname,"r");
  while (tab_next_statement(commsyntax,filehandle,line,TABREADLINE)) {
    j++;
  }
  fclose(filehandle);
  return j;
}


/* "maximum size <n>" on a Set declaration is silently ignored by
   modern GEMPACK (manual 11.7.2); excise it so legacy TABs parse
   (parity plan 3.3) */
static void set_maxsize_excise(char *line) {
  int k;
  char *p,*q;
  k=str_find_ci(line,"maximum size");
  if(k<0)return;
  p=line+k;
  q=p+12;
  while(*q==' '||*q=='\t')q++;
  while(*q>='0'&&*q<='9')q++;
  while(*q==' '||*q=='\t')q++;
  memmove(p,q,strlen(q)+1);
}

int sets_read_intertemporal(char *fname, int niodata, cmf_file_entry *iodata, set_def *record,dim_t nset) {
  FILE * filehandle;//, *fileout;
  char line[TABREADLINE]="\0",linecopy[TABREADLINE];//,line1[TABREADLINE];longname[TABREADLINE],
  char header[HEADERSIZE],floginame[NAMESIZE];
  char commsyntax[TABREADLINE],varname[NAMESIZE];
  dim_t j=0,i,intindx[4],sign[2],intvar[2];//,inttype[4]
  int k0,k1=0;
  char *readitem=NULL;//,*p;
  for (j=0; j<nset; j++) {
    if(str_find_ci(record[j].readele,"intertemporal")>-1) {
      intindx[0]=0;
      intindx[1]=0;
      intindx[2]=0;
      intindx[3]=0;
      sign[0]=0;
      sign[1]=0;
      strcpy(line,record[j].readele);
      while (str_replace_all(line," ", ""));
      strcpy(linecopy,line);
      k0=str_count_ci(line,",");
      if(k0==1) {
        readitem = strtok(line,",");
        readitem = strtok(NULL,"\0");
        if(readitem==NULL) {
          printf("Error: malformed intertemporal set declaration in TAB file\n");
          return -1;
        }
        if(strchr(readitem,'+')!=NULL||strchr(readitem,'-')!=NULL) {
          if (strchr(readitem,'+')!=NULL) {
            sign[0]=1;
          } else {
            sign[0]=2;
          }
          readitem = strtok(readitem,"-+");
          if(readitem==NULL) {
            printf("Error: malformed intertemporal set declaration in TAB file\n");
            return -1;
          }
          if(readitem[0]>='0'&&readitem[0]<='9') {
            intindx[0]=atoi(readitem);
          } else {
            if (strlen(readitem)>=sizeof(varname)) {
              printf("Error: malformed intertemporal set declaration in TAB file\n");
              return -1;
            }
            strcpy(varname,readitem);
            strcpy(commsyntax,"read ");
            strcat(commsyntax,readitem);
            filehandle = fopen(fname,"r");
            while (tab_next_statement(commsyntax,filehandle,line,TABREADLINE)) {
              readitem = strtok(line," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              if (readitem==NULL||strlen(readitem)>=sizeof(floginame)) {
                printf("Error: malformed intertemporal set declaration in TAB file\n");
                fclose(filehandle);
                return -1;
              }
              strcpy(floginame,readitem);
              readitem = strtok(NULL,"\"");
              readitem = strtok(NULL,"\"");
              if (readitem==NULL||strlen(readitem)>=sizeof(header)) {
                printf("Error: malformed intertemporal set declaration in TAB file\n");
                fclose(filehandle);
                return -1;
              }
              strcpy(header,readitem);
              for (k1=0; k1<niodata; k1++) if (strcmp(floginame,iodata[k1].logname)==0) {
                  break;
                }
              break;
            }
            fclose(filehandle);
            filehandle = fopen(iodata[k1].filname,"r");
            while (fgets(line,TABREADLINE,filehandle)) {
              if (str_count_char(line,'\"')>1) {
                readitem = strtok(line,"\"");
                readitem = strtok(NULL,"\"");
                if (strcmp(readitem,header) == 0) {
                  while (fgets(line,TABREADLINE,filehandle)) {
                    readitem = strtok(line,"\n");
                    intindx[0]=atoi(readitem);
                    break;
                  }
                  break;
                }
              }
            }
            fclose(filehandle);
          }
          strcpy(line,linecopy);
          readitem = strtok(line,",");
          readitem = strtok(NULL,"-+");
          readitem = strtok(NULL,"\0");
          if(readitem==NULL) {
            printf("Error: malformed intertemporal set declaration in TAB file\n");
            return -1;
          }
          if(readitem[0]>='0'&&readitem[0]<='9') {
            intindx[1]=atoi(readitem);
          } else {
            if (strlen(readitem)>=sizeof(varname)) {
              printf("Error: malformed intertemporal set declaration in TAB file\n");
              return -1;
            }
            strcpy(varname,readitem);
            strcpy(commsyntax,"read ");
            strcat(commsyntax,readitem);
            filehandle = fopen(fname,"r");
            while (tab_next_statement(commsyntax,filehandle,line,TABREADLINE)) {
              readitem = strtok(line," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              if (readitem==NULL||strlen(readitem)>=sizeof(floginame)) {
                printf("Error: malformed intertemporal set declaration in TAB file\n");
                fclose(filehandle);
                return -1;
              }
              strcpy(floginame,readitem);
              readitem = strtok(NULL,"\"");
              readitem = strtok(NULL,"\"");
              if (readitem==NULL||strlen(readitem)>=sizeof(header)) {
                printf("Error: malformed intertemporal set declaration in TAB file\n");
                fclose(filehandle);
                return -1;
              }
              strcpy(header,readitem);
              for (k1=0; k1<niodata; k1++) if (strcmp(floginame,iodata[k1].logname)==0) {
                  break;
                }
              break;
            }
            fclose(filehandle);
            filehandle = fopen(iodata[k1].filname,"r");
            while (fgets(line,TABREADLINE,filehandle)) {
              if (str_count_char(line,'\"')>1) {
                readitem = strtok(line,"\"");
                readitem = strtok(NULL,"\"");
                if (strcmp(readitem,header) == 0) {
                  while (fgets(line,TABREADLINE,filehandle)) {
                    readitem = strtok(line,"\n");
                    intindx[1]=atoi(readitem);
                    break;
                  }
                  break;
                }
              }
            }
            fclose(filehandle);
          }
        } else {
          if(readitem==NULL) {
            printf("Error: malformed intertemporal set declaration in TAB file\n");
            return -1;
          }
          if(readitem[0]>='0'&&readitem[0]<='9') {
            intindx[0]=atoi(readitem);
          } else {
            if (strlen(readitem)>=sizeof(varname)) {
              printf("Error: malformed intertemporal set declaration in TAB file\n");
              return -1;
            }
            strcpy(varname,readitem);
            strcpy(commsyntax,"read ");
            strcat(commsyntax,readitem);
            filehandle = fopen(fname,"r");
            while (tab_next_statement(commsyntax,filehandle,line,TABREADLINE)) {
              readitem = strtok(line," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              if (readitem==NULL||strlen(readitem)>=sizeof(floginame)) {
                printf("Error: malformed intertemporal set declaration in TAB file\n");
                fclose(filehandle);
                return -1;
              }
              strcpy(floginame,readitem);
              readitem = strtok(NULL,"\"");
              readitem = strtok(NULL,"\"");
              if (readitem==NULL||strlen(readitem)>=sizeof(header)) {
                printf("Error: malformed intertemporal set declaration in TAB file\n");
                fclose(filehandle);
                return -1;
              }
              strcpy(header,readitem);
              for (k1=0; k1<niodata; k1++) if (strcmp(floginame,iodata[k1].logname)==0) {
                  break;
                }
              break;
            }
            fclose(filehandle);
            filehandle = fopen(iodata[k1].filname,"r");
            while (fgets(line,TABREADLINE,filehandle)) {
              if (str_count_char(line,'\"')>1) {
                readitem = strtok(line,"\"");
                readitem = strtok(NULL,"\"");
                if (strcmp(readitem,header) == 0) {
                  while (fgets(line,TABREADLINE,filehandle)) {
                    readitem = strtok(line,"\n");
                    intindx[0]=atoi(readitem);
                    break;
                  }
                  break;
                }
              }
            }
            fclose(filehandle);
          }
        }
        if(sign[0]>0) if(sign[0]==1) {
            sprintf(line, "%d", intindx[0]+intindx[1]);
            strcpy(record[j].readele,line);
          } else {
            sprintf(line, "%d", intindx[0]-intindx[1]);
            strcpy(record[j].readele,line);
          }
        else {
          sprintf(line, "%d", intindx[0]);
          strcpy(record[j].readele,line);
        }
        record[j].size=1;
      } else {
        strcpy(line,linecopy);
        readitem = strtok(line,",");
        readitem = strtok(NULL,",");
        if(readitem==NULL) {
          printf("Error: malformed intertemporal set declaration in TAB file\n");
          return -1;
        }
        if(strchr(readitem,'+')!=NULL||strchr(readitem,'-')!=NULL) {
          if (strchr(readitem,'+')!=NULL) {
            sign[0]=1;
          } else {
            sign[0]=2;
          }
          readitem = strtok(readitem,"-+");
          if(readitem==NULL) {
            printf("Error: malformed intertemporal set declaration in TAB file\n");
            return -1;
          }
          if(readitem[0]>='0'&&readitem[0]<='9') {
            intindx[0]=atoi(readitem);
          } else {
            if (strlen(readitem)>=sizeof(varname)) {
              printf("Error: malformed intertemporal set declaration in TAB file\n");
              return -1;
            }
            strcpy(varname,readitem);
            strcpy(commsyntax,"read ");
            strcat(commsyntax,readitem);
            filehandle = fopen(fname,"r");
            while (tab_next_statement(commsyntax,filehandle,line,TABREADLINE)) {
              readitem = strtok(line," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              if (readitem==NULL||strlen(readitem)>=sizeof(floginame)) {
                printf("Error: malformed intertemporal set declaration in TAB file\n");
                fclose(filehandle);
                return -1;
              }
              strcpy(floginame,readitem);
              readitem = strtok(NULL,"\"");
              readitem = strtok(NULL,"\"");
              if (readitem==NULL||strlen(readitem)>=sizeof(header)) {
                printf("Error: malformed intertemporal set declaration in TAB file\n");
                fclose(filehandle);
                return -1;
              }
              strcpy(header,readitem);
              for (k1=0; k1<niodata; k1++) if (strcmp(floginame,iodata[k1].logname)==0) {
                  break;
                }
              break;
            }
            fclose(filehandle);
            filehandle = fopen(iodata[k1].filname,"r");
            while (fgets(line,TABREADLINE,filehandle)) {
              if (str_count_char(line,'\"')>1) {
                readitem = strtok(line,"\"");
                readitem = strtok(NULL,"\"");
                if (strcmp(readitem,header) == 0) {
                  while (fgets(line,TABREADLINE,filehandle)) {
                    readitem = strtok(line,"\n");
                    intindx[0]=atoi(readitem);
                    break;
                  }
                  break;
                }
              }
            }
            fclose(filehandle);
          }
          strcpy(line,linecopy);
          readitem = strtok(line,",");
          readitem = strtok(NULL,"-+");
          readitem = strtok(NULL,",");
          if(readitem==NULL) {
            printf("Error: malformed intertemporal set declaration in TAB file\n");
            return -1;
          }
          if(readitem[0]>='0'&&readitem[0]<='9') {
            intindx[1]=atoi(readitem);
          } else {
            if (strlen(readitem)>=sizeof(varname)) {
              printf("Error: malformed intertemporal set declaration in TAB file\n");
              return -1;
            }
            strcpy(varname,readitem);
            strcpy(commsyntax,"read ");
            strcat(commsyntax,readitem);
            filehandle = fopen(fname,"r");
            while (tab_next_statement(commsyntax,filehandle,line,TABREADLINE)) {
              readitem = strtok(line," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              if (readitem==NULL||strlen(readitem)>=sizeof(floginame)) {
                printf("Error: malformed intertemporal set declaration in TAB file\n");
                fclose(filehandle);
                return -1;
              }
              strcpy(floginame,readitem);
              readitem = strtok(NULL,"\"");
              readitem = strtok(NULL,"\"");
              if (readitem==NULL||strlen(readitem)>=sizeof(header)) {
                printf("Error: malformed intertemporal set declaration in TAB file\n");
                fclose(filehandle);
                return -1;
              }
              strcpy(header,readitem);
              for (k1=0; k1<niodata; k1++) if (strcmp(floginame,iodata[k1].logname)==0) {
                  break;
                }
              break;
            }
            fclose(filehandle);
            filehandle = fopen(iodata[k1].filname,"r");
            while (fgets(line,TABREADLINE,filehandle)) {
              if (str_count_char(line,'\"')>1) {
                readitem = strtok(line,"\"");
                readitem = strtok(NULL,"\"");
                if (strcmp(readitem,header) == 0) {
                  while (fgets(line,TABREADLINE,filehandle)) {
                    readitem = strtok(line,"\n");
                    intindx[1]=atoi(readitem);
                    break;
                  }
                  break;
                }
              }
            }
            fclose(filehandle);
          }
        } else {
          if(readitem==NULL) {
            printf("Error: malformed intertemporal set declaration in TAB file\n");
            return -1;
          }
          if(readitem[0]>='0'&&readitem[0]<='9') {
            intindx[0]=atoi(readitem);
          } else {
            if (strlen(readitem)>=sizeof(varname)) {
              printf("Error: malformed intertemporal set declaration in TAB file\n");
              return -1;
            }
            strcpy(varname,readitem);
            strcpy(commsyntax,"read ");
            strcat(commsyntax,readitem);
            filehandle = fopen(fname,"r");
            while (tab_next_statement(commsyntax,filehandle,line,TABREADLINE)) {
              readitem = strtok(line," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              if (readitem==NULL||strlen(readitem)>=sizeof(floginame)) {
                printf("Error: malformed intertemporal set declaration in TAB file\n");
                fclose(filehandle);
                return -1;
              }
              strcpy(floginame,readitem);
              readitem = strtok(NULL,"\"");
              readitem = strtok(NULL,"\"");
              if (readitem==NULL||strlen(readitem)>=sizeof(header)) {
                printf("Error: malformed intertemporal set declaration in TAB file\n");
                fclose(filehandle);
                return -1;
              }
              strcpy(header,readitem);
              for (k1=0; k1<niodata; k1++) if (strcmp(floginame,iodata[k1].logname)==0) {
                  break;
                }
              break;
            }
            fclose(filehandle);
            filehandle = fopen(iodata[k1].filname,"r");
            while (fgets(line,TABREADLINE,filehandle)) {
              if (str_count_char(line,'\"')>1) {
                readitem = strtok(line,"\"");
                readitem = strtok(NULL,"\"");
                if (strcmp(readitem,header) == 0) {
                  while (fgets(line,TABREADLINE,filehandle)) {
                    readitem = strtok(line,"\n");
                    intindx[0]=atoi(readitem);
                    break;
                  }
                  break;
                }
              }
            }
            fclose(filehandle);
          }
        }
        if(sign[0]>0) {
          if(sign[0]==1) {
            intvar[0]=intindx[0]+intindx[1];
          } else {
            intvar[0]=intindx[0]+intindx[1];
          }
        } else {
          intvar[0]=intindx[0];
        }
        strcpy(line,linecopy);
        readitem = strtok(line,",");
        readitem = strtok(NULL,",");
        readitem = strtok(NULL,"\0");
        if(readitem==NULL) {
          printf("Error: malformed intertemporal set declaration in TAB file\n");
          return -1;
        }
        if(strchr(readitem,'+')!=NULL||strchr(readitem,'-')!=NULL) {
          if (strchr(readitem,'+')!=NULL) {
            sign[1]=1;
          } else {
            sign[1]=2;
          }
          readitem = strtok(readitem,"-+");
          if(readitem==NULL) {
            printf("Error: malformed intertemporal set declaration in TAB file\n");
            return -1;
          }
          if(readitem[0]>='0'&&readitem[0]<='9') {
            intindx[2]=atoi(readitem);
          } else {
            if (strlen(readitem)>=sizeof(varname)) {
              printf("Error: malformed intertemporal set declaration in TAB file\n");
              return -1;
            }
            strcpy(varname,readitem);
            strcpy(commsyntax,"read ");
            strcat(commsyntax,readitem);
            filehandle = fopen(fname,"r");
            while (tab_next_statement(commsyntax,filehandle,line,TABREADLINE)) {
              readitem = strtok(line," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              if (readitem==NULL||strlen(readitem)>=sizeof(floginame)) {
                printf("Error: malformed intertemporal set declaration in TAB file\n");
                fclose(filehandle);
                return -1;
              }
              strcpy(floginame,readitem);
              readitem = strtok(NULL,"\"");
              readitem = strtok(NULL,"\"");
              if (readitem==NULL||strlen(readitem)>=sizeof(header)) {
                printf("Error: malformed intertemporal set declaration in TAB file\n");
                fclose(filehandle);
                return -1;
              }
              strcpy(header,readitem);
              for (k1=0; k1<niodata; k1++) if (strcmp(floginame,iodata[k1].logname)==0) {
                  break;
                }
              break;
            }
            fclose(filehandle);
            filehandle = fopen(iodata[k1].filname,"r");
            while (fgets(line,TABREADLINE,filehandle)) {
              if (str_count_char(line,'\"')>1) {
                readitem = strtok(line,"\"");
                readitem = strtok(NULL,"\"");
                if (strcmp(readitem,header) == 0) {
                  while (fgets(line,TABREADLINE,filehandle)) {
                    readitem = strtok(line,"\n");
                    intindx[2]=atoi(readitem);
                    break;
                  }
                  break;
                }
              }
            }
            fclose(filehandle);
          }
          strcpy(line,linecopy);
          readitem = strtok(line,",");
          readitem = strtok(NULL,",");
          readitem = strtok(NULL,"-+");
          readitem = strtok(NULL,"\0");
          if(readitem==NULL) {
            printf("Error: malformed intertemporal set declaration in TAB file\n");
            return -1;
          }
          if(readitem[0]>='0'&&readitem[0]<='9') {
            intindx[3]=atoi(readitem);
          } else {
            if (strlen(readitem)>=sizeof(varname)) {
              printf("Error: malformed intertemporal set declaration in TAB file\n");
              return -1;
            }
            strcpy(varname,readitem);
            strcpy(commsyntax,"read ");
            strcat(commsyntax,readitem);
            filehandle = fopen(fname,"r");
            while (tab_next_statement(commsyntax,filehandle,line,TABREADLINE)) {
              readitem = strtok(line," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              if (readitem==NULL||strlen(readitem)>=sizeof(floginame)) {
                printf("Error: malformed intertemporal set declaration in TAB file\n");
                fclose(filehandle);
                return -1;
              }
              strcpy(floginame,readitem);
              readitem = strtok(NULL,"\"");
              readitem = strtok(NULL,"\"");
              if (readitem==NULL||strlen(readitem)>=sizeof(header)) {
                printf("Error: malformed intertemporal set declaration in TAB file\n");
                fclose(filehandle);
                return -1;
              }
              strcpy(header,readitem);
              for (k1=0; k1<niodata; k1++) if (strcmp(floginame,iodata[k1].logname)==0) {
                  break;
                }
              break;
            }
            fclose(filehandle);
            filehandle = fopen(iodata[k1].filname,"r");
            while (fgets(line,TABREADLINE,filehandle)) {
              if (str_count_char(line,'\"')>1) {
                readitem = strtok(line,"\"");
                readitem = strtok(NULL,"\"");
                if (strcmp(readitem,header) == 0) {
                  while (fgets(line,TABREADLINE,filehandle)) {
                    readitem = strtok(line,"\n");
                    intindx[3]=atoi(readitem);
                    break;
                  }
                  break;
                }
              }
            }
            fclose(filehandle);
          }
        } else {
          if(readitem==NULL) {
            printf("Error: malformed intertemporal set declaration in TAB file\n");
            return -1;
          }
          if(readitem[0]>='0'&&readitem[0]<='9') {
            intindx[2]=atoi(readitem);
          } else {
            if (strlen(readitem)>=sizeof(varname)) {
              printf("Error: malformed intertemporal set declaration in TAB file\n");
              return -1;
            }
            strcpy(varname,readitem);
            strcpy(commsyntax,"read ");
            strcat(commsyntax,readitem);
            filehandle = fopen(fname,"r");
            while (tab_next_statement(commsyntax,filehandle,line,TABREADLINE)) {
              readitem = strtok(line," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              readitem = strtok(NULL," ");
              if (readitem==NULL||strlen(readitem)>=sizeof(floginame)) {
                printf("Error: malformed intertemporal set declaration in TAB file\n");
                fclose(filehandle);
                return -1;
              }
              strcpy(floginame,readitem);
              readitem = strtok(NULL,"\"");
              readitem = strtok(NULL,"\"");
              if (readitem==NULL||strlen(readitem)>=sizeof(header)) {
                printf("Error: malformed intertemporal set declaration in TAB file\n");
                fclose(filehandle);
                return -1;
              }
              strcpy(header,readitem);
              for (k1=0; k1<niodata; k1++) if (strcmp(floginame,iodata[k1].logname)==0) {
                  break;
                }
              break;
            }
            fclose(filehandle);
            filehandle = fopen(iodata[k1].filname,"r");
            while (fgets(line,TABREADLINE,filehandle)) {
              if (str_count_char(line,'\"')>1) {
                readitem = strtok(line,"\"");
                readitem = strtok(NULL,"\"");
                if (strcmp(readitem,header) == 0) {
                  while (fgets(line,TABREADLINE,filehandle)) {
                    readitem = strtok(line,"\n");
                    intindx[2]=atoi(readitem);
                    break;
                  }
                  break;
                }
              }
            }
            fclose(filehandle);
          }
        }
        if(sign[1]>0) if(sign[1]==1) {
            intvar[1]=intindx[2]+intindx[3];
          } else {
            intvar[1]=intindx[2]-intindx[3];
          }
        else {
          intvar[1]=intindx[2];
        }
        if (intvar[1]+1-intvar[0]<=0) {
          printf("Error: intertemporal set has an empty or inverted time range in TAB file\n");
          return -1;
        }
        record[j].readele[0]='\0';
        for (i=intvar[0]; i<intvar[1]+1; i++) {
          sprintf(line, "%d",i);
          if (strlen(record[j].readele)+strlen(line)+2>=sizeof(record[j].readele)) {
            printf("Error: intertemporal set time range is too large in TAB file\n");
            return -1;
          }
          strcat(record[j].readele,line);
          strcat(record[j].readele,",");
        }
        record[j].readele[strlen(record[j].readele)-1]='\0';
        record[j].size=intvar[1]+1-intvar[0];
      }
    }
  }
  return 0;
}

int sets_read(char *fname, int niodata, cmf_file_entry *iodata, set_def *record,dim_t nset) {
  FILE * filehandle;//, *fileout;
  char line[TABREADLINE]="\0",linecopy[TABREADLINE],line1[TABREADLINE],line2[TABREADLINE];
  char longname[TABREADLINE],tempvar[256];
  char *commsyntax="set";
  dim_t j=0,vsize,dim1,i;
  int k0,k1,k2;
  char *readitem=NULL;//,*p;

  filehandle = fopen(fname,"r");

  while (tab_next_statement(commsyntax,filehandle,line,TABREADLINE)) {
    set_maxsize_excise(line);
    /* audit A5: the explicit default qualifier must not substring-match
       the intertemporal route */
    str_replace_first(line,"(non_intertemporal)","");
    strcpy(linecopy,line);
    k2=str_find_ci(line,"intertemporal");
    if(k2>-1) {
      readitem = strtok(line+k2," ");
      readitem = strtok(NULL," ");
      if (readitem==NULL||strlen(readitem)>=sizeof(record[j].setname)) {
        printf("Error: malformed set declaration in TAB file\n");
        return -1;
      }
      strcpy(record[j].setname,readitem);
      strcpy(line1,"intertemporal,");
      readitem = strtok(NULL,"[");
      readitem = strtok(NULL,"]");
      if (readitem==NULL) {
        printf("Error: malformed set declaration in TAB file\n");
        return -1;
      }
      strcat(line1,readitem);
      readitem = strtok(NULL,";");
      if (readitem!=NULL&&strchr(readitem,'[')!=NULL) {
        strcat(line1,",");
        readitem = strtok(readitem,"[");
        readitem = strtok(NULL,"]");
        if (readitem==NULL) {
          printf("Error: malformed set declaration in TAB file\n");
          return -1;
        }
        strcat(line1,readitem);
      }
      strcpy(record[j].readele,line1);
      record[j].intertemp=true;
    } else {
      if (strchr(line,'=')!=NULL) {
        str_replace_all(line,"union", "^");
        /* tokenize line2, not line1: the strtok(NULL,...) tokens below
           stay live inside the buffer while line1 is rebuilt in place —
           tokenizing line1 made those strcat(line1,readitem) calls
           overlapping copies */
        readitem = &line[4];
        strcpy(line2,readitem);
        while (str_replace_all(line2," ", ""));
        readitem = strtok(line2,"=");
        strcpy(record[j].setname,readitem);
        /* GEMPACK set expressions (manual 10.1.1.1): UNION, INTERSECT,
           '+', '-' and '\' over already-declared sets and quoted
           single elements, with '(' ')' grouping, evaluated left to
           right. Normalized here (UNION -> '^', INTERSECT -> '&',
           '\' -> '-', spaces stripped) and stored as "@<expr>" for
           set_expr_build(). The single-operator unbracketed forms keep
           the legacy pairwise encodings below byte-for-byte. NB the
           keyword rewrite is substring-based (as the legacy path
           already was), so set names must not contain "union" or
           "intersect". */
        int expr_route=0;
        {
          char exprbuf[TABREADLINE];
          strcpy(exprbuf,linecopy);
          while (str_replace_first(exprbuf,"intersect","&"));
          while (str_replace_first(exprbuf,"union","^"));
          while (str_replace_all(exprbuf," ", ""));
          str_replace_char(exprbuf,'\\','-');
          char *crhs=strchr(exprbuf,'=');
          if (crhs!=NULL) {
            crhs++;
            char *p;
            int nops=0,special=0;
            for (p=crhs; *p!='\0'&&*p!=';'; p++) {
              if (*p=='+'||*p=='-'||*p=='^') nops++;
              if (*p=='&'||*p=='('||*p=='"') special=1;
              if (*p=='&') nops++;
            }
            if (nops>=2||special) {
              expr_route=1;
              /* size upper bound: sum of named-operand sizes plus one
                 per quoted element; corrected to the exact count when
                 the elements build (set_expr_build) */
              char nm[NAMESIZE];
              dim_t sz=0;
              int ni=0,inq=0;
              record[j].readele[0]='@';
              record[j].readele[1]='\0';
              for (p=crhs;; p++) {
                if (inq) {
                  if (*p=='"'||*p=='\0') inq=0;
                  if (*p=='\0') break;
                  continue;
                }
                if (*p=='+'||*p=='-'||*p=='^'||*p=='&'||*p=='('||*p==')'||*p=='"'||*p==';'||*p=='\0') {
                  if (ni>0) {
                    nm[ni]='\0';
                    ni=0;
                    for (i=0; i<nset; i++) if (strcmp(nm,record[i].setname)==0) break;
                    if (i==nset) {
                      printf("Error: set %s in the definition of %s is not declared before use\n",nm,record[j].setname);
                      return -1;
                    }
                    else sz+=record[i].size;
                  }
                  if (*p=='"') { sz++; inq=1; }
                  if (*p==';'||*p=='\0') break;
                } else if (ni<NAMESIZE-1) {
                  nm[ni++]=*p;
                }
              }
              record[j].size=sz;
              {
                size_t rl=strlen(record[j].readele);
                for (p=crhs; *p!='\0'&&*p!=';'&&rl<TABREADLINE-1; p++) record[j].readele[rl++]=*p;
                record[j].readele[rl]='\0';
              }
            }
          }
        }
        if (expr_route) {
          /* encoded above; element build happens in set_expr_build */
        } else if (strchr(line,'-')!=NULL) {
          line1[0]='-';
          line1[1]=',';
          line1[2]='\0';
          readitem = strtok(NULL,"-");
          if (readitem==NULL) { printf("Error: malformed set declaration in TAB file\n"); return -1; }
          for (i=0; i<nset; i++) {
            if (strcmp(readitem,record[i].setname)==0) {
              dim1=record[i].size;
              strcat(line1,readitem);
              strcat(line1,",");
              break;
            }
          }
          readitem = strtok(NULL,";");
          if (readitem==NULL) { printf("Error: malformed set declaration in TAB file\n"); return -1; }
          for (i=0; i<nset; i++) {
            if (strcmp(readitem,record[i].setname)==0) {
              record[j].size=dim1-record[i].size;
              strcat(line1,readitem);
              break;
            }
          }
          if (record[j].size<0) {
            printf("Error: set difference subtracts a larger set in TAB file\n");
            return -1;
          }
          strcpy(record[j].readele,line1);
        } else {
          if(strchr(line,'+')!=NULL) {
            line1[0]='+';
            line1[1]=',';
            line1[2]='\0';
            readitem = strtok(NULL,"+");
            if (readitem==NULL) { printf("Error: malformed set declaration in TAB file\n"); return -1; }
            for (i=0; i<nset; i++) {
              if (strcmp(readitem,record[i].setname)==0) {
                dim1=record[i].size;
                strcat(line1,readitem);
                strcat(line1,",");
                break;
              }
            }
            readitem = strtok(NULL,";");
            if (readitem==NULL) { printf("Error: malformed set declaration in TAB file\n"); return -1; }
            for (i=0; i<nset; i++) {
              if (strcmp(readitem,record[i].setname)==0) {
                record[j].size=dim1+record[i].size;
                strcat(line1,readitem);
                break;
              }
            }
            strcpy(record[j].readele,line1);
          } else {
            if(strchr(line,'^')!=NULL) {
              line1[0]='^';
              line1[1]=',';
              line1[2]='\0';
              readitem = strtok(NULL,"^");
              if (readitem==NULL) { printf("Error: malformed set declaration in TAB file\n"); return -1; }
              for (i=0; i<nset; i++) {
                if (strcmp(readitem,record[i].setname)==0) {
                  dim1=record[i].size;
                  strcat(line1,readitem);
                  strcat(line1,",");
                  break;
                }
              }
              readitem = strtok(NULL,";");
              if (readitem==NULL) { printf("Error: malformed set declaration in TAB file\n"); return -1; }
              for (i=0; i<nset; i++) {
                if (strcmp(readitem,record[i].setname)==0) {
                  record[j].size=dim1+record[i].size;
                  strcat(line1,readitem);
                  break;
                }
              }
              strcpy(record[j].readele,line1);
            } else {
              line1[0]='=';
              line1[1]='\0';
              readitem = strtok(NULL,";");
              if (readitem==NULL) { printf("Error: malformed set declaration in TAB file\n"); return -1; }
              for (i=0; i<nset; i++) {
                if (strcmp(readitem,record[i].setname)==0) {
                  if (i==j) {
                    printf("Error: set %s is defined as equal to itself in TAB file\n",record[j].setname);
                    return -1;
                  }
                  record[j].size=record[i].size;
                  /* audit A2: do NOT copy the header -- a copied header
                     sent main down the read-from-file path with a stale
                     fileid and skipped set_equality_build (elements +
                     the two implied SUBSETs); the '=' readele drives
                     the equality path for every source kind */
                  sprintf(tempvar, "%d",i);
                  strcat(line1,tempvar);
                  break;
                }
              }
              if (i>=nset) {
                printf("Error: the right-hand side of 'Set %s = %s' is not a declared set (set products 'A x B' and data-dependent sets are not supported)\n",record[j].setname,readitem);
                return -1;
              }
              strcpy(record[j].readele,line1);
            }
          }
        }
      } else {
        k0=str_find_ci(line,"read elements from file");
        if(k0>-1) {
          k1=str_find_ci(line+k0+24," ");
          if (k1<0||k1>=(int)sizeof(line1)) {
            printf("Error: malformed 'read elements from file' clause in set declaration in TAB file\n");
            return -1;
          }
          strncpy(line1,line+k0+24,k1);
          line1[k1]='\0';
          for (k0=0; k0<niodata; k0++) if (strcmp(line1,iodata[k0].logname)==0) {
              break;
            }
          record[j].fileid=k0;
        }
        readitem = strtok(line," ");
        readitem = strtok(NULL," ");
        if (readitem==NULL||strlen(readitem)>=sizeof(record[j].setname)) {
          printf("Error: malformed set declaration in TAB file\n");
          return -1;
        }
        strcpy(record[j].setname,readitem);
        readitem = strtok(NULL,"\"");
        readitem = strtok(NULL,"\"");
        if (readitem!=NULL) {
          if (strlen(readitem)>=sizeof(record[j].header)) {
            printf("Error: header too long in set declaration in TAB file: %s\n",readitem);
            return -1;
          }
          strcpy(record[j].header,readitem);
          record[j].readele[0]='\0';
          datafile_read_header_info(readitem,iodata[record[j].fileid].filname,&vsize,longname,&dim1);
          record[j].size=dim1;
        } else {
          record[j].header[0]='\0';
          readitem = strtok(linecopy,"(");
          if (readitem!=NULL) readitem = strtok(NULL,")");
          if (readitem==NULL) {
            printf("Error: malformed set declaration in TAB file\n");
            return -1;
          }
          if (strchr(readitem,'-')!=NULL) {
            printf("Error: element range abbreviation '(first - last)' in set %s is not supported; list the elements explicitly\n",record[j].setname);
            return -1;
          }
          strcpy(record[j].readele,readitem);
          dim1=str_count_char(readitem, ',');
          record[j].size=dim1+1;
        }
      }
    }
    j++;
  }
  fclose(filehandle);
  return 0;
}

dim_t set_union_named(set_element *set_elems, set_def *sets,dim_t nset,dim_t i) {
  dim_t j,l,n,m,dim1=0,dim2=0,j1,sup1=MAXSUPSET,sup2=MAXSUPSET;
  char line[TABREADLINE],*readitem;
  dim_t bound=sets[i].size; /* parse-time upper bound = allocated element slots */
  strcpy(line,sets[i].readele);
  readitem = strtok(line,",");
  readitem = strtok(NULL,",");
  if (readitem==NULL) return 0;
  for (l=0; l<nset; l++) {
    if (strcmp(readitem,sets[l].setname)==0) {
      dim1=sets[l].size;
      for (sup1=1; sup1<MAXSUPSET; sup1++)if(sets[l].subsetid[sup1]==-1) {
          sets[l].subsetid[sup1]=i;
          break;
        }
      break;
    }
  }
  if(sup1>=MAXSUPSET){printf("Error: superset size exceeded; increase MAXSUPSET in teems_solver.h\n");return 0;}
  if (dim1>0&&l==i) {printf("Error: set %s references itself in a set expression\n",sets[i].setname);return 0;}
  m=0;
  for (n=0; n<dim1; n++) {
    if (m>=bound) {printf("Error: set expression produces more elements than declared for set %s\n",sets[i].setname);return 0;}
    if(sets[i].offset+m!=sets[l].offset+n)strcpy(set_elems[sets[i].offset+m].setele,set_elems[sets[l].offset+n].setele);
    set_elems[sets[i].offset+m].superset_pos[0]=n;
    set_elems[sets[l].offset+m].superset_pos[sup1]=n;
    m++;
  }
  readitem = strtok(NULL,",");
  if (readitem==NULL) return 0;
  for (j=0; j<nset; j++) {
    if (strcmp(readitem,sets[j].setname)==0) {
      dim2=sets[j].size;
      for (sup2=1; sup2<MAXSUPSET; sup2++)if(sets[j].subsetid[sup2]==-1) {
          sets[j].subsetid[sup2]=i;
          break;
        }
      break;
    }
  }
  if(sup2>=MAXSUPSET){printf("Error: superset size exceeded; increase MAXSUPSET in teems_solver.h\n");return 0;}
  if (dim2>0&&j==i) {printf("Error: set %s references itself in a set expression\n",sets[i].setname);return 0;}
  for (n=0; n<dim2; n++) {
    for (j1=0; j1<dim1; j1++) if(strcmp(set_elems[sets[l].offset+j1].setele,set_elems[sets[j].offset+n].setele)==0) break;
    if(j1==dim1) {
      if (m>=bound) {printf("Error: set expression produces more elements than declared for set %s\n",sets[i].setname);return 0;}
      if(sets[i].offset+m!=sets[j].offset+n)strcpy(set_elems[sets[i].offset+m].setele,set_elems[sets[j].offset+n].setele);
      set_elems[sets[i].offset+m].superset_pos[0]=m;
      set_elems[sets[j].offset+n].superset_pos[sup2]=m;
      m++;
    } else {
      set_elems[sets[j].offset+n].superset_pos[sup2]=j1;
    }
  }
  sets[i].size=m;
  return m;
}
dim_t set_union_op(set_element *set_elems, set_def *sets,dim_t nset,dim_t i) {
  dim_t j,l,n,m,dim1=0,dim2=0,j1,sup1=MAXSUPSET,sup2=MAXSUPSET;
  char line[TABREADLINE],*readitem;
  dim_t bound=sets[i].size; /* parse-time upper bound = allocated element slots */
  strcpy(line,sets[i].readele);
  readitem = strtok(line,",");
  readitem = strtok(NULL,",");
  if (readitem==NULL) return 0;
  for (l=0; l<nset; l++) {
    if (strcmp(readitem,sets[l].setname)==0) {
      dim1=sets[l].size;
      for (sup1=1; sup1<MAXSUPSET; sup1++)if(sets[l].subsetid[sup1]==-1) {
          sets[l].subsetid[sup1]=i;
          break;
        }
      break;
    }
  }
  if(sup1>=MAXSUPSET){printf("Error: superset size exceeded; increase MAXSUPSET in teems_solver.h\n");return 0;}
  if (dim1>0&&l==i) {printf("Error: set %s references itself in a set expression\n",sets[i].setname);return 0;}
  m=0;
  for (n=0; n<dim1; n++) {
    if (m>=bound) {printf("Error: set expression produces more elements than declared for set %s\n",sets[i].setname);return 0;}
    if(sets[i].offset+m!=sets[l].offset+n)strcpy(set_elems[sets[i].offset+m].setele,set_elems[sets[l].offset+n].setele);
    set_elems[sets[i].offset+m].superset_pos[0]=n;
    set_elems[sets[l].offset+m].superset_pos[sup1]=n;
    m++;
  }
  readitem = strtok(NULL,",");
  if (readitem==NULL) return 0;
  for (j=0; j<nset; j++) {
    if (strcmp(readitem,sets[j].setname)==0) {
      dim2=sets[j].size;
      for (sup2=1; sup2<MAXSUPSET; sup2++)if(sets[j].subsetid[sup2]==-1) {
          sets[j].subsetid[sup2]=i;
          break;
        }
      break;
    }
  }
  if(sup2>=MAXSUPSET){printf("Error: superset size exceeded; increase MAXSUPSET in teems_solver.h\n");return 0;}
  if (dim2>0&&j==i) {printf("Error: set %s references itself in a set expression\n",sets[i].setname);return 0;}
  for (n=0; n<dim2; n++) {
    if (m>=bound) {printf("Error: set expression produces more elements than declared for set %s\n",sets[i].setname);return 0;}
    if(sets[i].offset+m!=sets[j].offset+n)strcpy(set_elems[sets[i].offset+m].setele,set_elems[sets[j].offset+n].setele);
    set_elems[sets[i].offset+m].superset_pos[0]=m;
    set_elems[sets[j].offset+n].superset_pos[sup2]=m;
    m++;
  }
  sets[i].size=m; /* correct a parse-time upper bound (expression operands) */
  return m;
}
dim_t set_difference(set_element *set_elems, set_def *sets,dim_t nset,dim_t i) {
  dim_t j,l,n,m,v,indi,dim1=0,dim2=0,sup1=MAXSUPSET;
  char line[TABREADLINE],*readitem;
  dim_t bound=sets[i].size; /* parse-time upper bound = allocated element slots */
  strcpy(line,sets[i].readele);
  readitem = strtok(line,",");
  readitem = strtok(NULL,",");
  if (readitem==NULL) return 0;
  for (l=0; l<nset; l++) {
    if (strcmp(readitem,sets[l].setname)==0) {
      dim1=sets[l].size;
      for (sup1=1; sup1<MAXSUPSET; sup1++)if(sets[i].subsetid[sup1]==-1) {
          sets[i].subsetid[sup1]=l;
          break;
        }
      break;
    }
  }
  if(dim1>0&&sup1>=MAXSUPSET){printf("Error: superset size exceeded; increase MAXSUPSET in teems_solver.h\n");return 0;}
  if (dim1>0&&l==i) {printf("Error: set %s references itself in a set expression\n",sets[i].setname);return 0;}
  readitem = strtok(NULL,",");
  if (readitem==NULL) return 0;
  for (j=0; j<nset; j++) {
    if (strcmp(readitem,sets[j].setname)==0) {
      dim2=sets[j].size;
      break;
    }
  }
  if (dim2>0&&j==i) {printf("Error: set %s references itself in a set expression\n",sets[i].setname);return 0;}
  m=0;
  for (n=0; n<dim1; n++) {
    indi=0;
    for (v=0; v<dim2; v++) {
      if (strcmp(set_elems[sets[l].offset+n].setele,set_elems[sets[j].offset+v].setele)==0) {
        indi=1;
        break;
      }
    }
    if(indi==0) {
      if (m>=bound) {printf("Error: set expression produces more elements than declared for set %s\n",sets[i].setname);return 0;}
      if(sets[i].offset+m!=sets[l].offset+n)strcpy(set_elems[sets[i].offset+m].setele,set_elems[sets[l].offset+n].setele);
      set_elems[sets[i].offset+m].superset_pos[sup1]=n;
      set_elems[sets[i].offset+m].superset_pos[0]=m;
      m++;
    }
  }
  sets[i].size=m; /* correct a parse-time upper bound (expression operands) */
  return m;
}
/* ---- GEMPACK set expressions (manual 10.1.1.1) -----------------------
   readele "@<expr>" holds the normalized RHS of
     SET <name> = <set-expression>;
   with UNION as '^', INTERSECT as '&', '\' folded into '-', spaces
   stripped. Terms are already-declared set names, quoted single
   elements, or parenthesized subexpressions; operators apply left to
   right. Validity (per the manual): '+' operands must be disjoint;
   '-' may only remove elements that are present. */

/* record sup as a superset of sub and fill element positions */
static void set_register_subset(set_element *se, set_def *sets, dim_t sub, dim_t sup) {
  dim_t s,n,j1;
  for (s=1; s<MAXSUPSET; s++) {
    if (sets[sub].subsetid[s]==sup) return; /* already registered */
    if (sets[sub].subsetid[s]==-1) break;
  }
  if (s==MAXSUPSET) {
    printf("Error: superset count exceeds MAXSUPSET; increase MAXSUPSET in teems_solver.h\n");
    return;
  }
  sets[sub].subsetid[s]=sup;
  for (n=0; n<sets[sub].size; n++) {
    for (j1=0; j1<sets[sup].size; j1++)
      if (strcmp(se[sets[sub].offset+n].setele,se[sets[sup].offset+j1].setele)==0) break;
    if (j1<sets[sup].size) se[sets[sub].offset+n].superset_pos[s]=j1;
  }
}

static dim_t set_expr_eval(char **pp, set_element *se, set_def *sets, dim_t nset,
                           char (*out)[NAMESIZE], dim_t cap, const char *owner);

/* one term into out; returns element count */
static dim_t set_expr_term(char **pp, set_element *se, set_def *sets, dim_t nset,
                           char (*out)[NAMESIZE], dim_t cap, const char *owner) {
  char *p=*pp;
  dim_t n=0,l,k;
  if (*p=='(') {
    p++;
    *pp=p;
    n=set_expr_eval(pp,se,sets,nset,out,cap,owner);
    if (**pp==')') (*pp)++;
    else printf("Error: unbalanced '(' in the definition of %s\n",owner);
    return n;
  }
  if (*p=='"') {
    /* quoted single element; stored lowercase like read elements */
    p++;
    k=0;
    while (*p!='"'&&*p!='\0'&&k<NAMESIZE-1) out[0][k++]=tolower((int)*p++);
    out[0][k]='\0';
    if (*p=='"') p++;
    else printf("Error: unterminated quote in the definition of %s\n",owner);
    *pp=p;
    return 1;
  }
  {
    char nm[NAMESIZE];
    k=0;
    while (*p!='\0'&&*p!='+'&&*p!='-'&&*p!='^'&&*p!='&'&&*p!='('&&*p!=')'&&*p!='"'&&k<NAMESIZE-1) nm[k++]=*p++;
    nm[k]='\0';
    *pp=p;
    for (l=0; l<nset; l++) if (strcmp(nm,sets[l].setname)==0) break;
    if (l==nset) {
      printf("Error: set %s in the definition of %s is not declared before use\n",nm,owner);
      return 0;
    }
    for (k=0; (dim_t)k<sets[l].size&&(dim_t)k<cap; k++) strcpy(out[k],se[sets[l].offset+k].setele);
    if ((dim_t)k<sets[l].size) printf("Error: set %s has more elements than the declared size of %s\n",sets[l].setname,owner);
    return k;
  }
}

static dim_t set_expr_eval(char **pp, set_element *se, set_def *sets, dim_t nset,
                           char (*out)[NAMESIZE], dim_t cap, const char *owner) {
  dim_t n,m,a,b,w;
  char op;
  char (*tmp)[NAMESIZE];
  n=set_expr_term(pp,se,sets,nset,out,cap,owner);
  while (**pp=='+'||**pp=='-'||**pp=='^'||**pp=='&') {
    op=**pp;
    (*pp)++;
    tmp=malloc((size_t)cap*NAMESIZE);
    m=set_expr_term(pp,se,sets,nset,tmp,cap,owner);
    if (op=='+') {
      for (b=0; b<m; b++) {
        for (a=0; a<n; a++) if (strcmp(out[a],tmp[b])==0) break;
        if (a<n) printf("Error: '+' operands in the definition of %s are not disjoint (element %s); use UNION for overlapping sets\n",owner,tmp[b]);
        else if (n<cap) strcpy(out[n++],tmp[b]);
      }
    } else if (op=='^') {
      for (b=0; b<m; b++) {
        for (a=0; a<n; a++) if (strcmp(out[a],tmp[b])==0) break;
        if (a==n&&n<cap) strcpy(out[n++],tmp[b]);
      }
    } else if (op=='-') {
      for (b=0; b<m; b++) {
        for (a=0; a<n; a++) if (strcmp(out[a],tmp[b])==0) break;
        if (a==n) printf("Error: '-' in the definition of %s removes element %s, which is not present\n",owner,tmp[b]);
        else {
          for (w=a; w<n-1; w++) strcpy(out[w],out[w+1]);
          n--;
        }
      }
    } else { /* '&' INTERSECT */
      w=0;
      for (a=0; a<n; a++) {
        for (b=0; b<m; b++) if (strcmp(out[a],tmp[b])==0) break;
        if (b<m) {
          if (w!=a) strcpy(out[w],out[a]);
          w++;
        }
      }
      n=w;
    }
    free(tmp);
  }
  return n;
}

/* Build elements for set i from its "@<expr>" encoding, then apply the
   manual's implied-SUBSET rules:
     - all operators UNION/'+' (at any depth): every RHS set is a
       subset of the result;
     - all operators INTERSECT: the result is a subset of every RHS set;
     - expression ends "... UNION <set>" at top level: that set is a
       subset of the result;
     - expression ends "... INTERSECT <set>" at top level: the result
       is a subset of that set.
   (The simple two-set complement keeps its legacy path and rule.)
   Anything beyond these must be declared with an explicit Subset
   statement, exactly as in GEMPACK. */
dim_t set_expr_build(set_element *se, set_def *sets, dim_t nset, dim_t i) {
  char expr[TABREADLINE],*p,*cursor;
  char (*out)[NAMESIZE];
  dim_t cap=sets[i].size,m,n,l;
  int depth=0,inq=0,allplusun=1,allint=1,anyop=0;
  char lastop=0;
  char lastterm[NAMESIZE];
  lastterm[0]='\0';
  strcpy(expr,sets[i].readele+1);
  if (cap<1) cap=1;
  out=malloc((size_t)cap*NAMESIZE);
  cursor=expr;
  m=set_expr_eval(&cursor,se,sets,nset,out,cap,i<nset?sets[i].setname:"?");
  if (*cursor!='\0') printf("Error: trailing characters in the definition of %s: %s\n",sets[i].setname,cursor);
  for (n=0; n<m; n++) {
    strcpy(se[sets[i].offset+n].setele,out[n]);
    se[sets[i].offset+n].superset_pos[0]=n;
  }
  free(out);
  sets[i].size=m;
  /* operator census + last top-level term for the implied rules */
  {
    int k=0;
    for (p=expr; *p!='\0'; p++) {
      if (inq) {
        if (*p=='"') inq=0;
        continue;
      }
      if (*p=='"') { inq=1; lastterm[0]='\0'; k=0; continue; }
      if (*p=='(') { depth++; lastterm[0]='\0'; k=0; continue; }
      if (*p==')') { depth--; k=0; continue; }
      if (*p=='+'||*p=='-'||*p=='^'||*p=='&') {
        anyop=1;
        if (*p=='-') { allplusun=0; allint=0; }
        else if (*p=='&') allplusun=0;
        else allint=0;
        if (depth==0) { lastop=*p; lastterm[0]='\0'; k=0; }
        continue;
      }
      if (depth==0) {
        if (k<NAMESIZE-1) { lastterm[k++]=*p; lastterm[k]='\0'; }
      }
    }
  }
  if (anyop&&(allplusun||allint)) {
    /* register against every named RHS set */
    char nm[NAMESIZE];
    int k=0;
    inq=0;
    for (p=expr;; p++) {
      if (inq) {
        if (*p=='"'||*p=='\0') inq=0;
        if (*p=='\0') break;
        continue;
      }
      if (*p=='+'||*p=='-'||*p=='^'||*p=='&'||*p=='('||*p==')'||*p=='"'||*p=='\0') {
        if (k>0) {
          nm[k]='\0';
          k=0;
          for (l=0; l<nset; l++) if (strcmp(nm,sets[l].setname)==0) break;
          if (l<nset) {
            if (allplusun) set_register_subset(se,sets,l,i);
            else set_register_subset(se,sets,i,l);
          }
        }
        if (*p=='"') inq=1;
        if (*p=='\0') break;
      } else if (k<NAMESIZE-1) {
        nm[k++]=*p;
      }
    }
  } else if (lastterm[0]!='\0'&&(lastop=='^'||lastop=='&')) {
    for (l=0; l<nset; l++) if (strcmp(lastterm,sets[l].setname)==0) break;
    if (l<nset) {
      if (lastop=='^') set_register_subset(se,sets,l,i);
      else set_register_subset(se,sets,i,l);
    }
  }
  return m;
}

/* SET <new> = <old>; (manual 10.1.2.1): readele "=<source idx>" from
   sets_read. Elements are copied at the source's CURRENT size
   (expression sources carry a parse-time upper bound) and the two
   implied SUBSET statements are registered in both directions. */
void set_equality_build(set_element *se, set_def *sets, dim_t i) {
  dim_t j,j1;
  dim_t bound=sets[i].size; /* parse-time upper bound = allocated element slots */
  j1=(dim_t)atoi(sets[i].readele+1);
  sets[i].size=sets[j1].size;
  for (j=0; j<sets[i].size; j++) {
    if (j>=bound) {
      printf("Error: set equality produces more elements than declared for set %s\n",sets[i].setname);
      sets[i].size=j;
      return;
    }
    strcpy(se[sets[i].offset+j].setele,se[sets[j1].offset+j].setele);
    se[sets[i].offset+j].superset_pos[0]=j;
  }
  set_register_subset(se,sets,i,j1);
  set_register_subset(se,sets,j1,i);
}

offset_t subsets_read(char *fname, set_element *set_elems, set_def *sets,dim_t nset) {
  FILE * filehandle;//, *fileout;
  char line[TABREADLINE]="\0";
  char set[NAMESIZE],subset[NAMESIZE];
  char *commsyntax="subset";
  dim_t i,setd,subsetd,sup1;//,nlength;
  offset_t jj,jjj,j=0,succ=0,ssize=0;
  char *readitem=NULL;

  filehandle = fopen(fname,"r");
  if (filehandle==NULL) return -1;

  while (tab_next_statement(commsyntax,filehandle,line,TABREADLINE)) {
    while (str_replace_all(line," ;", ";"));
    readitem = strtok(line," ");
    readitem = strtok(NULL," ");
    if (readitem==NULL||strlen(readitem)>=sizeof(subset)) {
      printf("Error: malformed %s statement in TAB file\n",commsyntax);
      fclose(filehandle);
      return -1;
    }
    strcpy(subset,readitem);
    readitem = strtok(NULL," ");
    readitem = strtok(NULL," ");
    readitem = strtok(NULL," ");
    readitem = strtok(NULL,";");
    if (readitem==NULL||strlen(readitem)>=sizeof(set)) {
      printf("Error: malformed %s statement in TAB file\n",commsyntax);
      fclose(filehandle);
      return -1;
    }
    strcpy(set,readitem);
    succ=0;
    for (i=0; i<nset; i++) {
      if (strcmp(sets[i].setname,subset)==0) { //,subsetd)==0) {
        for (j=0; j<nset; j++) {
          if (strcmp(sets[j].setname,set)==0) { //,setd)==0) {
            ssize=sets[j].size;
            for (sup1=1; sup1<MAXSUPSET; sup1++)if(sets[i].subsetid[sup1]==-1) {
                sets[i].subsetid[sup1]=j;
                break;
              }
            if(sup1==MAXSUPSET) {
              printf("Error: superset count exceeds MAXSUPSET; increase MAXSUPSET in teems_solver.h\n");
              fclose(filehandle);
              return -1;
            }
            for (jj=sets[j].offset; jj<sets[j].offset+sets[j].size; jj++) {
              for (jjj=sets[i].offset; jjj<sets[i].offset+sets[i].size; jjj++) {
                if (strcmp(set_elems[jj].setele,set_elems[jjj].setele)==0) {
                  set_elems[jjj].superset_pos[sup1]=set_elems[jj].superset_pos[0];
                  succ++;
                  break;
                }
              }
            }
            break;
          }
        }
        if(j==nset)printf("Error: set %s is not declared\n",set);
        break;
      }
      j++;
    }
    if(i==nset)printf("Error: set %s is not declared\n",subset);
    else if(succ-sets[i].size!=0)printf("Error: some elements of set %s are not in set %s\n",subset,set);
  }
  fclose(filehandle);
  return j;
}

dim_t subset_map_build(set_element *set_elems, set_def *sets,dim_t nset,offset_t* contin) {
  dim_t i,setd,subsetd,sup1,sup2,sup3,b;//,nlength;
  offset_t jj,jjj,j=0;
  *contin=0;
  for (i=0; i<nset; i++) {
    b=1;
    while(b==1) {
      for (sup1=1; sup1<MAXSUPSET; sup1++) {
        b=0;
        if(sets[i].subsetid[sup1]==-1) {
          b=2;
          break;
        }
        for (sup2=1; sup2<MAXSUPSET; sup2++) {
          b=0;
          if(sets[sets[i].subsetid[sup1]].subsetid[sup2]==-1) {
            b=2;
            break;
          }
          for (sup3=1; sup3<MAXSUPSET; sup3++){
            if(sets[sets[i].subsetid[sup1]].subsetid[sup2]==sets[i].subsetid[sup3]) {
              b=2;
              break;
            }
          }
          if(b==0)break;
        }
        if(b==0)break;
      }
      if(b==0) {
        *contin=1;
        j=sets[sets[i].subsetid[sup1]].subsetid[sup2];
        for (sup1=1; sup1<MAXSUPSET; sup1++){
          if(sets[i].subsetid[sup1]==-1)break;
        }
        if(sup1==MAXSUPSET){
          printf("Error: superset size exceeded; increase MAXSUPSET in teems_solver.h\n");
          *contin=0;
          return 0;
        }
        sets[i].subsetid[sup1]=j;
        for (jj=sets[j].offset; jj<sets[j].offset+sets[j].size; jj++) {
          for (jjj=sets[i].offset; jjj<sets[i].offset+sets[i].size; jjj++) {
            if (strcmp(set_elems[jj].setele,set_elems[jjj].setele)==0) {
              set_elems[jjj].superset_pos[sup1]=set_elems[jj].superset_pos[0];
              break;
            }
          }
        }
        b=1;
      }
    }
  }
  for (i=nset-1; i>-1; i--) {
    b=1;
    while(b==1) {
      for (sup1=1; sup1<MAXSUPSET; sup1++) {
        b=0;
        if(sets[i].subsetid[sup1]==-1) {
          b=2;
          break;
        }
        for (sup2=1; sup2<MAXSUPSET; sup2++) {
          b=0;
          if(sets[sets[i].subsetid[sup1]].subsetid[sup2]==-1) {
            b=2;
            break;
          }
          for (sup3=1; sup3<MAXSUPSET; sup3++){
            if(sets[sets[i].subsetid[sup1]].subsetid[sup2]==sets[i].subsetid[sup3]) {
              b=2;
              break;
            }
          }
          if(b==0)break;
        }
        if(b==0)break;
      }
      if(b==0) {
        *contin=1;
        j=sets[sets[i].subsetid[sup1]].subsetid[sup2];
        for (sup1=1; sup1<MAXSUPSET; sup1++){
          if(sets[i].subsetid[sup1]==-1)break;
        }
        if(sup1==MAXSUPSET){
          printf("Error: superset size exceeded; increase MAXSUPSET in teems_solver.h\n");
          *contin=0;
          return 0;
        }
        sets[i].subsetid[sup1]=j;
        for (jj=sets[j].offset; jj<sets[j].offset+sets[j].size; jj++) {
          for (jjj=sets[i].offset; jjj<sets[i].offset+sets[i].size; jjj++) {
            if (strcmp(set_elems[jj].setele,set_elems[jjj].setele)==0) {
              set_elems[jjj].superset_pos[sup1]=set_elems[jj].superset_pos[0];
              break;
            }
          }
        }
        b=1;
      }
    }
  }
  return 1;
}

char *closure_next_statement(char *commsyntax, FILE *filehandle, char *readline) {
  int check1=0,i,count1=0;
  while (commsyntax[count1] != '\0') {
    count1++;
  }
  char uppercomsyn[TABLINESIZE],*n,*p;
  char line[TABLINESIZE],linecomm[TABLINESIZE],upperlinecomm[TABLINESIZE],*finditem=";";//,linecopy[TABLINESIZE+2]
  for (i=0; i<count1; i++) {
    uppercomsyn[ i ] = toupper( (int) commsyntax[ i ] );
  }
  strcpy(readline,"\0");

  while (fgets(line,TABLINESIZE,filehandle)) {
    while (str_replace_all(line,"  ", " "));
    if (line[0]==' ') {
      memmove(line,&line[0]+1,strlen(line)-1);
    }
    strncpy(linecomm,line,count1);
    for ( i = 0; i<count1; i++) {
      upperlinecomm[ i ] = toupper( (int) linecomm[ i ] );
    }
    if (strncmp(upperlinecomm,uppercomsyn,count1) == 0&&check1==0) {
      p=strpbrk(line,"!");
      while(p!=NULL) {
        line[p-line]='\0';
        p=strpbrk(line,"!");
      }
      check1=1;
      n=strstr(line,finditem);//ha_cgefendofc
      if (n==NULL) {
        strcpy(readline,line);
      } else {
        strcpy(readline,line);
        return readline;
      }
      while (fgets(line,TABLINESIZE,filehandle)) {
        p=strpbrk(line,"!");
        while(p!=NULL) {
          line[p-line]='\0';
          p=strpbrk(line,"!");
        }
        n=strstr(line,finditem);//ha_cgefendofc
        if (n==NULL) {
          strcat(readline, line);
        } else {
          strcat(readline, line);
          return readline;
        }
      }

    }
  }
  return NULL;
}

/* Extract the equation name from a full "equation ..." statement into
   eqname (lowercased by tab_preprocess; an optional "(linear)" qualifier
   may precede the name).  Returns 0 for option statements such as
   "equation (default=...)" that carry no name. */
int tab_equation_name(char *stmt, char *eqname) {
  char *p=stmt+strlen("equation");
  int i=0;
  while (*p==' '||*p=='\n'||*p=='\r'||*p=='\t') p++;
  if (strncmp(p,"(linear)",8)==0) {
    p+=8;
    while (*p==' '||*p=='\n'||*p=='\r'||*p=='\t') p++;
  }
  if (*p=='('||*p=='\0') return 0;
  while (*p!='\0'&&*p!=' '&&*p!='('&&*p!='\n'&&*p!='\r'&&*p!='\t'&&*p!='#'&&*p!=';'&&i<NAMESIZE-1) {
    eqname[i]=tolower((int)*p);
    i++;
    p++;
  }
  eqname[i]='\0';
  return i>0;
}

/* Backsolve filter for "equation" scans: in SKIP mode the nominated
   defining equations are invisible to every consumer (count, ordering,
   preallocation, fill, structural detection), so they all see the same
   condensed system; ONLY mode inverts the filter for the recovery build.
   Returns 1 when stmt must be rejected by the current scan. */
static int backsolve_eq_reject(char *commsyntax, char *stmt) {
  char eqname[NAMESIZE];
  int i,hit=0;
  if (nbacksolve==0) return 0;
  if (strcmp(commsyntax,"equation")!=0) return 0;
  if (tab_equation_name(stmt,eqname)) {
    for (i=0; i<nbacksolve; i++) {
      if (strcmp(backsolves[i].eqname,eqname)==0) {
        hit=1;
        break;
      }
    }
  }
  if (backsolve_scan_mode==BS_SCAN_ONLY) return !hit;
  return hit;
}

static char *tab_next_statement_raw(char *commsyntax, FILE *filehandle, char *readline,offset_t rlinesize) {
  int check1=0,i,count1=0;
  while (commsyntax[count1] != '\0') {
    count1++;
  }
  char uppercomsyn[TABLINESIZE],*n;
  char line[TABLINESIZE],linecomm[TABLINESIZE],upperlinecomm[TABLINESIZE],*finditem=";";//,linecopy[TABLINESIZE+2]
  for (i=0; i<count1; i++) {
    uppercomsyn[ i ] = toupper( (int) commsyntax[ i ] );
  }
  strcpy(readline,"\0");

  while (fgets(line,TABLINESIZE,filehandle)) {
    while (str_replace_all(line,"  ", " "));
    if (line[0]==' ') {
      memmove(line,&line[0]+1,strlen(line)-1);
    }
    strncpy(linecomm,line,count1);
    for ( i = 0; i<count1; i++) {
      upperlinecomm[ i ] = toupper( (int) linecomm[ i ] );
    }
    if (strncmp(upperlinecomm,uppercomsyn,count1) == 0&&check1==0) {
      check1=1;
      n=strstr(line,finditem);//ha_cgefendofc
      if (n==NULL) {
        strcpy(readline,line);
      } else {
        strcpy(readline,line);
          if(strlen(readline)>=((size_t)(9*rlinesize/10)))printf("Warning: statement near the line-buffer limit (%ld of %ld bytes) in tab_next_statement; increase TABREADLINE/DATREADLINE\n",strlen(readline),rlinesize);
        return readline;
      }
      while (fgets(line,TABLINESIZE,filehandle)) {
        n=strstr(line,finditem);//ha_cgefendofc
        if (n==NULL) {
          strcat(readline, line);
        } else {
          strcat(readline, line);
          if(strlen(readline)>=((size_t)(9*rlinesize/10)))printf("Warning: statement near the line-buffer limit (%ld of %ld bytes) in tab_next_statement; increase TABREADLINE/DATREADLINE\n",strlen(readline),rlinesize);
          return readline;
        }
      }

    }
  }
  return NULL;
}
char *tab_next_statement(char *commsyntax, FILE *filehandle, char *readline,offset_t rlinesize) {
  char *r;
  while ((r=tab_next_statement_raw(commsyntax,filehandle,readline,rlinesize))!=NULL) {
    if (backsolve_eq_reject(commsyntax,r)) continue;
    return r;
  }
  return NULL;
}

static char *tab_next_statement_resolved_raw(char *commsyntax, FILE *filehandle, char *readline, elem_value *record, array_def *coefs,offset_t ncof,solve_real *zerodivide,offset_t rlinesize) {
  int check1=0,count1=0;
  char *zerosyntax="zerodivide default",*p,*zerosyntax1="zerodivide (",*zerosyntax2="zerodivide(";
  while (commsyntax[count1] != '\0') {
    count1++;
  }
  char *n,line[TABLINESIZE],*finditem=";";//,linecopy[TABLINESIZE+2]
  strcpy(readline,"\0");
  while (fgets(line,TABLINESIZE,filehandle)) {
    while (str_replace_all(line,"  ", " "));
    if (line[0]==' ') {
      memmove(line,&line[0]+1,strlen(line)-1);
    }
    if (strncmp(line,commsyntax,count1) == 0&&check1==0) {
      check1=1;
      n=strstr(line,finditem);//ha_cgefendofc
      if (n==NULL) {
        strcpy(readline,line);
      } else {
        strcpy(readline,line);
          if(strlen(readline)>=((size_t)(9*rlinesize/10)))printf("Warning: statement near the line-buffer limit (%ld of %ld bytes) in tab_next_statement_resolved; increase TABREADLINE/DATREADLINE\n",strlen(readline),rlinesize);
        return readline;
      }
      while (fgets(line,TABLINESIZE,filehandle)) {
        n=strstr(line,finditem);//ha_cgefendofc
        if (n==NULL) {
          strcat(readline, line);
        } else {
          strcat(readline, line);
          if(strlen(readline)>=((size_t)(9*rlinesize/10)))printf("Warning: statement near the line-buffer limit (%ld of %ld bytes) in tab_next_statement_resolved; increase TABREADLINE/DATREADLINE\n",strlen(readline),rlinesize);
          return readline;
        }
      }

    }
    if (strncmp(line,commsyntax,count1) != 0) {
      /* dual-class state (plan A1), parsed from a copy BEFORE the
         legacy strtok chain mangles the line; the legacy single
         default below stays byte-identical for -gpzerodivide 0 */
      if (strncmp(line,"zerodivide",10) == 0) {
        char zline[TABLINESIZE];
        char *zp;
        solve_real zval;
        int znbz;
        strcpy(zline,line);
        znbz=(strstr(zline,"nonzero_by_zero")!=NULL);
        /* DEFAULT checked first: OFF is only legal without a value, so
           " off" inside a scalar-coefficient default name cannot misroute */
        if ((zp=strstr(zline,"default"))!=NULL) {
          zp+=7;
          zp=strtok(zp,";");
          zval=formula_subst_scalar(zp,record,coefs,ncof);
          if(znbz) {
            teems_zdiv_scan.nbz_val=zval;
            teems_zdiv_scan.nbz_on=1;
          } else {
            teems_zdiv_scan.zbz_val=zval;
            teems_zdiv_scan.zbz_on=1;
          }
        } else if (strstr(zline," off")!=NULL||strstr(zline,")off")!=NULL) {
          if(znbz)teems_zdiv_scan.nbz_on=0;
          else teems_zdiv_scan.zbz_on=0;
        }
      }
      if (strncmp(line,zerosyntax,18) == 0) {
        p=strtok(&line[18],";");
        *zerodivide=formula_subst_scalar(p,record,coefs,ncof);//atof(p);
      }
      if (strncmp(line,zerosyntax1,12) == 0||strncmp(line,zerosyntax2,11) == 0) {
        p=strtok(line,")");
        p=strtok(NULL,"default");
        p=strtok(NULL," ");
        p=strtok(NULL,";");
        *zerodivide=formula_subst_scalar(p,record,coefs,ncof);//atof(p);
      }
    }
  }
  return NULL;
}

char *tab_next_statement_resolved(char *commsyntax, FILE *filehandle, char *readline, elem_value *record, array_def *coefs,offset_t ncof,solve_real *zerodivide,offset_t rlinesize) {
  char *r;
  while ((r=tab_next_statement_resolved_raw(commsyntax,filehandle,readline,record,coefs,ncof,zerodivide,rlinesize))!=NULL) {
    if (backsolve_eq_reject(commsyntax,r)) continue;
    return r;
  }
  return NULL;
}

/* Read the "backsolve <var> using <eq> ;" statements (GEMPACK manual
   10.16; lowercased and one-per-line after tab_preprocess).  Marks every
   element of each named variable in closure_vals and assigns its compact
   slot in the recovered-value array.  Returns the total number of
   backsolved elements, or -1 on an invalid statement. */
offset_t backsolve_read(char *fname, array_def *vars, offset_t nvar, closure_entry *closure_vals) {
  FILE *filehandle;
  char line[TABREADLINE]="\0",commsyntax[NAMESIZE],*p,*vname,*using,*eqname,*term;
  offset_t i,j,l;
  strcpy(commsyntax,"backsolve");
  filehandle=fopen(fname,"r");
  if (filehandle==NULL) {
    printf("Error: cannot open %s\n",fname);
    return -1;
  }
  while (tab_next_statement_raw(commsyntax,filehandle,line,TABREADLINE)) {
    while (str_replace_all(line,"\n"," "));
    while (str_replace_all(line,"\r"," "));
    while (str_replace_all(line,"  "," "));
    p=line+strlen("backsolve");
    vname=strtok(p," ");
    using=strtok(NULL," ");
    eqname=strtok(NULL," ;");
    term=strtok(NULL," ;");
    if (vname==NULL||using==NULL||eqname==NULL||strcmp(using,"using")!=0||(term!=NULL&&term[0]!='\0')) {
      printf("Error: malformed backsolve statement \"%s\"; expected \"backsolve <variable> using <equation> ;\"\n",line);
      fclose(filehandle);
      return -1;
    }
    j=-1;
    for (i=0; i<nvar; i++) {
      if (strcmp(vars[i].cofname,vname)==0) {
        j=i;
        break;
      }
    }
    /* GEMPACK 10.16: the linear name p_X/c_X may stand for the variable */
    if (j==-1&&(strncmp(vname,"p_",2)==0||strncmp(vname,"c_",2)==0)) {
      for (i=0; i<nvar; i++) {
        if (strcmp(vars[i].cofname,vname+2)==0) {
          j=i;
          break;
        }
      }
    }
    if (j==-1) {
      printf("Error: backsolve names variable %s but no such variable is declared\n",vname);
      fclose(filehandle);
      return -1;
    }
    for (i=0; i<nbacksolve; i++) {
      if (backsolves[i].varindx==j) {
        printf("Error: variable %s is backsolved more than once\n",vars[j].cofname);
        fclose(filehandle);
        return -1;
      }
      if (strcmp(backsolves[i].eqname,eqname)==0) {
        printf("Error: equation %s is nominated by more than one backsolve statement\n",eqname);
        fclose(filehandle);
        return -1;
      }
    }
    for (l=vars[j].offset; l<vars[j].offset+vars[j].nelem; l++) {
      if (closure_vals[l].is_exogenous) {
        printf("Error: backsolved variable %s is exogenous in the closure; a backsolved variable must be endogenous (GEMPACK manual 14.1.3)\n",vars[j].cofname);
        fclose(filehandle);
        return -1;
      }
      closure_vals[l].is_backsolved=true;
      closure_vals[l].exo_index=(exo_idx_t)(nbselems+(l-vars[j].offset));
    }
    backsolves=realloc(backsolves,(nbacksolve+1)*sizeof(backsolve_def));
    strncpy(backsolves[nbacksolve].eqname,eqname,NAMESIZE-1);
    backsolves[nbacksolve].eqname[NAMESIZE-1]='\0';
    backsolves[nbacksolve].varindx=j;
    backsolves[nbacksolve].elem_base=nbselems;
    nbacksolve++;
    nbselems+=vars[j].nelem;
  }
  fclose(filehandle);
  /* omit/substitute are symbolic condensation actions the solver cannot
     perform; teems-R resolves them during model preparation */
  strcpy(commsyntax,"omit");
  filehandle=fopen(fname,"r");
  if (tab_next_statement_raw(commsyntax,filehandle,line,TABREADLINE)!=NULL) {
    printf("Error: the TAB file contains an omit statement; omission is resolved during model preparation (ems_model(omit=)) and must not reach the solver\n");
    fclose(filehandle);
    return -1;
  }
  fclose(filehandle);
  strcpy(commsyntax,"substitute");
  filehandle=fopen(fname,"r");
  if (tab_next_statement_raw(commsyntax,filehandle,line,TABREADLINE)!=NULL) {
    printf("Error: the TAB file contains a substitute statement; substitution is resolved during model preparation (ems_model(backsolve=)) and must not reach the solver\n");
    fclose(filehandle);
    return -1;
  }
  fclose(filehandle);
  return nbselems;
}

/* Validate the condensed system's references once the backsolve pairs are
   known: every nominated defining equation must exist and reference its
   backsolved variable, and no retained equation may reference any
   backsolved variable (teems-R eliminates such references during
   condensation; a leftover reference means the TAB and the backsolve
   statements disagree).  Equation text references linear variables with
   the p_ prefix added by tab_write_variables.  Returns -1 on error. */
int backsolve_validate_refs(char *fname, array_def *vars) {
  FILE *filehandle;
  char line[TABREADLINE]="\0",commsyntax[NAMESIZE],eqname[NAMESIZE],ref[NAMESIZE+4];
  int i,k,hit;
  offset_t l;
  int *eqfound;
  if (nbacksolve==0) return 0;
  eqfound=calloc(nbacksolve,sizeof(int));
  strcpy(commsyntax,"equation");
  filehandle=fopen(fname,"r");
  if (filehandle==NULL) {
    printf("Error: cannot open %s\n",fname);
    free(eqfound);
    return -1;
  }
  while (tab_next_statement_raw(commsyntax,filehandle,line,TABREADLINE)) {
    if (strstr(line,"(default")!=NULL) continue;
    if (!tab_equation_name(line,eqname)) continue;
    hit=-1;
    for (i=0; i<nbacksolve; i++) {
      if (strcmp(backsolves[i].eqname,eqname)==0) {
        hit=i;
        break;
      }
    }
    if (hit>=0) eqfound[hit]=1;
    /* M2b: backsolve recovery inverts a row->element bijection of the
       defining equation; a mapped (many-to-one) reference breaks it by
       construction, so a nominated defining equation carrying lowered
       mapping tokens is rejected up front rather than failing the
       runtime one-to-one check with a generic message */
    if (hit>=0&&strchr(line,'@')!=NULL) {
      printf("Error: defining equation %s for backsolved variable %s uses a set mapping; backsolving through mapped references is not supported (the row-to-element recovery is no longer one-to-one)\n",eqname,vars[backsolves[hit].varindx].cofname);
      fclose(filehandle);
      free(eqfound);
      MPI_Abort(PETSC_COMM_WORLD,1);
      return -1;
    }
    for (i=0; i<nbacksolve; i++) {
      strcpy(ref,"p_");
      strcat(ref,vars[backsolves[i].varindx].cofname);
      /* a reference is p_<name> delimited on both sides (so p_pm inside
         p_pms or exp_pm does not count) */
      k=str_find_ci(line,ref);
      while (k>-1) {
        l=k+strlen(ref);
        if ((k==0||(!isalnum((int)line[k-1])&&line[k-1]!='_'))
            &&!isalnum((int)line[l])&&line[l]!='_') break;
        {
          int k2=str_find_ci(line+l,ref);
          k=(k2==-1)?-1:(int)(l+k2);
        }
      }
      if (k>-1) {
        if (hit==i) continue;                     /* its own defining equation */
        if (hit>=0) {
          printf("Error: defining equation %s (backsolves %s) also references backsolved variable %s; each retained defining equation may reference only surviving variables and its own backsolved variable\n",eqname,vars[backsolves[hit].varindx].cofname,vars[backsolves[i].varindx].cofname);
        }
        else {
          printf("Error: equation %s references backsolved variable %s; a backsolved variable must be eliminated from every retained equation (redeploy the model so the condensation rewrites this equation)\n",eqname,vars[backsolves[i].varindx].cofname);
        }
        fclose(filehandle);
        free(eqfound);
        return -1;
      }
      if (hit==i) {
        printf("Error: defining equation %s does not reference its backsolved variable %s\n",eqname,vars[backsolves[i].varindx].cofname);
        fclose(filehandle);
        free(eqfound);
        return -1;
      }
    }
  }
  fclose(filehandle);
  for (i=0; i<nbacksolve; i++) {
    if (!eqfound[i]) {
      printf("Error: backsolve for %s nominates equation %s but no such equation is declared\n",vars[backsolves[i].varindx].cofname,backsolves[i].eqname);
      free(eqfound);
      return -1;
    }
  }
  free(eqfound);
  return 0;
}

char *str_replace_first(char *line, char *finditem, char *replitem) {
  char buffer[DATREADLINE];
  int count2 = 0,index=0;
  while (finditem[count2] != '\0') {
    count2++;
  }
  index=str_find_ci(line, finditem);
  if (index==-1) {
    return NULL;
  }
  strncpy(buffer, line, index);
  buffer[index] = '\0';
  strcat(buffer, replitem);
  strcat(buffer, line+index+count2);
  strcpy(line,buffer);
  return line;
}

char *str_replace_first_bounded(char *line, char *finditem, char *replitem,dim_t nbuffer) {
  char *buffer= (char *) calloc (nbuffer,sizeof(char));
  int count2 = 0,index=0;
  while (finditem[count2] != '\0') {
    count2++;
  }
  index=str_find_ci(line, finditem);
  if (index==-1) {
    free(buffer);
    return NULL;
  }
  strncpy(buffer, line, index);
  buffer[index] = '\0';
  strcat(buffer, replitem);
  strcat(buffer, line+index+count2);
  strcpy(line,buffer);
  free(buffer);
  return line;
}

char *str_replace_all(char *line, char *finditem, char *replitem) {
  char buffer[DATREADLINE];
  char *p;
  unsigned short int count2 = 0,index;
  if (line==NULL) return NULL;
  while (finditem[count2] != '\0') {
    count2++;
  }

  if (!(p = strstr(line, finditem))) {
    return NULL;
  }
  index=p-line;
  strncpy(buffer, line, index);
  buffer[index] = '\0';
  strcat(buffer, replitem);
  strcat(buffer, (line+index+count2));
  strcpy(line,buffer);
  return line;
}

/* Replace every occurrence of finditem with replitem in `line` (buffer
   capacity `linesz`), in place. For the parsers' set-symbol substitution
   this is behaviourally identical to `while (str_replace_all(line,
   finditem, replitem));` — replitem never contains finditem, so the
   left-to-right restart terminates — but it refuses to grow the result
   past linesz-1 rather than overflowing `line` via str_replace_all's
   unbounded strcpy-back. Returns 0 on success, -1 if the result would
   not fit (valid models stay well under the buffer, so this fires only
   on the malformed/pathological input the old code overflowed on). */
static int str_subst_all_bounded(char *line, const char *finditem, const char *replitem, size_t linesz) {
  size_t flen = strlen(finditem), rlen = strlen(replitem);
  char *p, *start = line;
  if (flen == 0) return 0;
  while ((p = strstr(start, finditem)) != NULL) {
    size_t taillen = strlen(p + flen);
    if (strlen(line) - flen + rlen + 1 > linesz) return -1;
    memmove(p + rlen, p + flen, taillen + 1);
    memcpy(p, replitem, rlen);
    start = p + rlen; /* resume past the replacement so a replitem containing
                         finditem cannot re-match and loop forever */
  }
  return 0;
}

char *str_replace_all_bounded(char *line, char *finditem, char *replitem,dim_t nbuffer) {
  char *buffer= (char *) calloc (nbuffer,sizeof(char));
  char *p;
  long int count2 = 0,index;
  while (finditem[count2] != '\0') {
    count2++;
  }

  if (!(p = strstr(line, finditem))) {
    free(buffer);
    return NULL;
  }
  index=p-line;
  strncpy(buffer, line, index);
  buffer[index] = '\0';
  strcat(buffer, replitem);
  strcat(buffer, (line+index+count2));
  strcpy(line,buffer);
  free(buffer);
  return line;
}

char *str_replace_char(char *line, int finditem, int replitem) {
  char *p;
  if (!(p = strchr(line, finditem))) {
    return NULL;
  }
  *p=replitem;
  return line;
}

int str_find_ci(char *line, char *finditem) {
  int i;
  char *t1,*t2;
  for (i=0; line[i]; i++) {
    t1 = &line[i];
    t2=&finditem[0];
    while (*t2 && tolower((int)*t2)==tolower((int)*t1)) {
      t1++;
      t2++;
    }
    if (!*t2) {
      return i;
    }
  }
  return -1;
}//ref: http://www.java2s.com/Code/C/String/Findsubstringourownfunction.htm

char *str_strip_comment(char *line, char *token) {
  char buffer[TABREADLINE];
  char *p,*readitem=NULL;
  int count1 = 0;
  while (line[count1] != '\0') {
    count1++;
  }

  if (!(p = strstr(line, token))) {
    return NULL;
  }
  strncpy(buffer, line, p-line);
  buffer[p-line] = '\0';
  strcat(buffer," ");
  readitem = strtok(line,token);
  if ((p-line)!=0) {
    readitem = strtok(NULL,token);
  }
  readitem = strtok(NULL,"\0");
  if(readitem!=NULL) {
    strcat(buffer,readitem);
  }
  strcpy(line,buffer);
  return line;
}

