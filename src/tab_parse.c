#include <teems_solver.h>

int formula_normalize(char *fomulain) {
  int index,i,i1,j;
  char fpart1[TABREADLINE],*p=NULL;//,*p1=NULL;
  while (str_replace_all(fomulain, " ", ""));
  p=strchr(fomulain,')');
  if (p==NULL) {
    while (str_replace_char(fomulain, '[', '('));
    while (str_replace_char(fomulain, ']', ')'));
    return 0;
  }
  strncpy(fpart1, fomulain, p-fomulain);
  fpart1[p-fomulain] = '\0';
  index=strrchr(fpart1,'(')-fpart1;//ha_cgerevfind(fpart1,"(");
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
      if (i==3||i1==3) {
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
      if (i==4||i1==4) {
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
      if (i==3||i1==3) {
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
      if (i==4) {
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

offset_t data_read_files(char *fname, int niodata, cmf_file_entry *iodata, char *commsyntax,set_def *ha_set,dim_t nset, set_element *ha_setele,array_def *ha_cof,offset_t ncof, elem_store *ha_cofele,offset_t ncofele,array_def *ha_var,offset_t nvar, elem_store *ha_varele,offset_t nvarele) {
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
    printf("Error!!! No such %s file!\n",fname);
    return -1;
  }

  while (tab_next_statement(commsyntax,filehandle,line,DATREADLINE)) {
    strcpy(linecopy,line);
    k0=str_find_ci(line,"from file ");
    if(k0>-1) {
      k1=str_find_ci(line+k0+10," ");
      strncpy(line1,line+k0+10,k1);
      line1[k1]='\0';
      for (k0=0; k0<niodata; k0++) if (strcmp(line1,iodata[k0].logname)==0) {
          break;
        }
      if(k0==niodata){
        printf("Error!!!! No such file %s in cmf file!\n",line1);
        return -1;
      }
    }
    if (strchr(line,'(')==NULL) {
      readitem = strtok(line," ");
      readitem = strtok(NULL," ");
      strcpy(vname,readitem);
      readitem = strtok(NULL,"\"");
      readitem = strtok(NULL,"\"");
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
        vname1= strtok(ha_var[i].cofname,"(");
        if (strcmp(vname1,vname)==0) { //,vsize ha_var[i].cofname
          ha_var[i].suplval=true;
          dims=1;
          for (n1=0; n1<ha_var[i].size; n1++) {
            dims=dims*ha_set[ha_var[i].setid[n1]].size;
          }
          antidim[0]=ha_set[ha_var[i].setid[1]].size;
          antidim[1]=1;
          antidim[2]=antidim[0]*ha_set[ha_var[i].setid[0]].size;
          if (ha_var[i].size>3) {
            for (n1=3; n1<ha_var[i].size; n1++) {
              antidim[n1]=antidim[n1-1]*ha_set[ha_var[i].setid[n1-1]].size;
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
            printf("Error!!! No such %s file!\n",iodata[k0].filname);
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
                      if(ha_var[i].size>1) {
                        for (n1=ha_var[i].size-1; n1>1; n1--) {
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
                      for (n1=0; n1<ha_var[i].size; n1++) {
                        recount1=recount1+index[n1]*ha_var[i].strides[n1];
                      }
                      ha_varele[ha_var[i].offset+recount1].value=val;
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
                        if(ha_var[i].size>1) {
                          for (n1=ha_var[i].size-1; n1>1; n1--) {
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
                        for (n1=0; n1<ha_var[i].size; n1++) {
                          recount1=recount1+index[n1]*ha_var[i].strides[n1];
                        }
                        ha_varele[ha_var[i].offset+recount1].value=val;
                        recount++;
                      }
                      readitem = strtok(NULL,"\n");
                      val=atof(readitem);//sscanf(readitem, "%lf", &val);
                      l1=recount;
                      if(ha_var[i].size>1) {
                        for (n1=ha_var[i].size-1; n1>1; n1--) {
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
                      for (n1=0; n1<ha_var[i].size; n1++) {
                        recount1=recount1+index[n1]*ha_var[i].strides[n1];
                      }
                      ha_varele[ha_var[i].offset+recount1].value=val;
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
            printf("Error!!!! %s is not a header in data file!\n",header);
            return -1;
          }
          fclose(filehandle1);
          break;
        }
      }
      if(i<nvar)count2=1;
      for (i=0; i<ncof; i++) {
        vname1= strtok(ha_cof[i].cofname,"(");
        if (strcmp(vname1,vname)==0) { //,vsize ha_cof[i].cofname
          ha_cof[i].suplval=true;
          dims=1;
          for (n1=0; n1<ha_cof[i].size; n1++) {
            dims=dims*ha_set[ha_cof[i].setid[n1]].size;
          }
          antidim[0]=ha_set[ha_cof[i].setid[1]].size;
          antidim[1]=1;
          antidim[2]=antidim[0]*ha_set[ha_cof[i].setid[0]].size;
          if (ha_cof[i].size>3) {
            for (n1=3; n1<ha_cof[i].size; n1++) {
              antidim[n1]=antidim[n1-1]*ha_set[ha_cof[i].setid[n1-1]].size;
            }
          }
          filehandle1 = fopen(iodata[k0].filname,"r");
          if(filehandle1==NULL){
            printf("Error!!! No such %s file!\n",iodata[k0].filname);
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
                      if(ha_cof[i].size>1) {
                        for (n1=ha_cof[i].size-1; n1>1; n1--) {
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
                      for (n1=0; n1<ha_cof[i].size; n1++) {
                        recount1=recount1+index[n1]*ha_cof[i].strides[n1];
                      }
                      ha_cofele[ha_cof[i].offset+recount1].value=val;
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
                        if(ha_cof[i].size>1) {
                          for (n1=ha_cof[i].size-1; n1>1; n1--) {
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
                        for (n1=0; n1<ha_cof[i].size; n1++) {
                          recount1=recount1+index[n1]*ha_cof[i].strides[n1];
                        }
                        ha_cofele[ha_cof[i].offset+recount1].value=val;
                        recount++;
                      }
                      readitem = strtok(NULL,"\n");
                      val=atof(readitem);//sscanf(readitem, "%lf", &val);
                      l1=recount;
                      if(ha_cof[i].size>1) {
                        for (n1=ha_cof[i].size-1; n1>1; n1--) {
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
                      for (n1=0; n1<ha_cof[i].size; n1++) {
                        recount1=recount1+index[n1]*ha_cof[i].strides[n1];
                      }
                      ha_cofele[ha_cof[i].offset+recount1].value=val;
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
            printf("Error1!!!! %s is not a header in data file!\n",header);
            return -1;
          }
          fclose(filehandle1);
          break;
        }
      }
      if(count2==0&&i==ncof){
        printf("Error!!!! %s is not a variable coefficient or parameter!\n",vname);
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
      strcpy(argu,readitem);
      strcat(argu,",");
      readitem = strtok(NULL,"\"");
      readitem = strtok(NULL,"\"");
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
        vname1= strtok(ha_cof[i].cofname,"(");
        if (strcmp(vname1,vname)==0) { //,vsize ha_var[i].cofname
          dims=1;
          for (n1=0; n1<ha_cof[i].size; n1++) {
            if (strcmp(varset[n1],ha_set[ha_cof[i].setid[n1]].setname)==0) {
              dims=dims*ha_set[ha_cof[i].setid[n1]].size;
              dim[n1]=ha_set[ha_cof[i].setid[n1]].size;
              offset[n1]=ha_set[ha_cof[i].setid[n1]].offset;
              supsetid[n1]=0;
            } else {
              for (n=0; n<nset; n++) {
                if (strcmp(varset[n1],ha_set[n].setname)==0) {
                  dims=dims*ha_set[n].size;
                  dim[n1]=ha_set[n].size;
                  offset[n1]=ha_set[n].offset;
                  for(l1=0; l1<MAXSUPSET; l1++)if(strcmp(ha_set[ha_set[n].subsetid[l1]].setname,ha_set[ha_cof[i].setid[n1]].setname)==0) {
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
          if (ha_cof[i].size>3) {
            for (n1=3; n1<ha_cof[i].size; n1++) {
              antidim[n1]=antidim[n1-1]*dim[n1-1];
            }
          }
          filehandle1 = fopen(iodata[k0].filname,"r");
          if(filehandle1==NULL){
            printf("Error!!! No such %s file!\n",iodata[k0].filname);
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
                printf("dim %s\n",readitem);
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
                      if (ha_cof[i].size>1) {
                        for (n1=ha_cof[i].size-1; n1>1; n1--) {
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
                      for (n1=0; n1<ha_cof[i].size; n1++) {
                        recount1=recount1+ha_setele[offset[n1]+index[n1]].superset_pos[supsetid[n1]]*ha_cof[i].strides[n1];
                      }
                      ha_cofele[ha_cof[i].offset+recount1].value=val;
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
                        if (ha_cof[i].size>1) {
                          for (n1=ha_cof[i].size-1; n1>1; n1--) {
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
                        for (n1=0; n1<ha_cof[i].size; n1++) {
                          recount1=recount1+ha_setele[offset[n1]+index[n1]].superset_pos[supsetid[n1]]*ha_cof[i].strides[n1];
                        }
                        ha_cofele[ha_cof[i].offset+recount1].value=val;
                        recount++;
                      }
                      readitem = strtok(NULL,"\n");
                      val=atof(readitem);//sscanf(readitem, "%lf", &val);
                      l1=recount;
                      if (ha_cof[i].size>1) {
                        for (n1=ha_cof[i].size-1; n1>1; n1--) {
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
                      for (n1=0; n1<ha_cof[i].size; n1++) {
                        recount1=recount1+ha_setele[offset[n1]+index[n1]].superset_pos[supsetid[n1]]*ha_cof[i].strides[n1];
                      }
                      ha_cofele[ha_cof[i].offset+recount1].value=val;
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
            printf("Error!!!! %s is not a header in data file!\n",header);
            return -1;
          }
          fclose(filehandle1);
          break;
        }
      }
      if(i<ncof)count2=1;
      for (i=0; i<nvar; i++) {
        vname1= strtok(ha_var[i].cofname,"(");
        if (strcmp(vname1,vname)==0) { //,vsize ha_var[i].cofname
          dims=1;
          for (n1=0; n1<ha_var[i].size; n1++) {
            if (strcmp(varset[n1],ha_set[ha_var[i].setid[n1]].setname)==0) {
              dims=dims*ha_set[ha_var[i].setid[n1]].size;
              dim[n1]=ha_set[ha_var[i].setid[n1]].size;
              offset[n1]=ha_set[ha_var[i].setid[n1]].offset;
              supsetid[n1]=0;
            } else {
              for (n=0; n<nset; n++) {
                if (strcmp(varset[n1],ha_set[n].setname)==0) {
                  dims=dims*ha_set[n].size;
                  dim[n1]=ha_set[n].size;
                  offset[n1]=ha_set[n].offset;
                  for(l1=0; l1<MAXSUPSET; l1++)if(strcmp(ha_set[ha_set[n].subsetid[l1]].setname,ha_set[ha_var[i].setid[n1]].setname)==0) {
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
          if (ha_var[i].size>3) {
            for (n1=3; n1<ha_var[i].size; n1++) {
              antidim[n1]=antidim[n1-1]*dim[n1-1];
            }
          }
          filehandle1 = fopen(iodata[k0].filname,"r");
          if(filehandle1==NULL){
            printf("Error!!! No such %s file!\n",iodata[k0].filname);
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
                      if (ha_var[i].size>1) {
                        for (n1=ha_var[i].size-1; n1>1; n1--) {
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
                      for (n1=0; n1<ha_var[i].size; n1++) {
                        recount1=recount1+ha_setele[offset[n1]+index[n1]].superset_pos[supsetid[n1]]*ha_var[i].strides[n1];
                      }
                      ha_varele[ha_var[i].offset+recount1].value=val;
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
                        if (ha_var[i].size>1) {
                          for (n1=ha_var[i].size-1; n1>1; n1--) {
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
                        for (n1=0; n1<ha_var[i].size; n1++) {
                          recount1=recount1+ha_setele[offset[n1]+index[n1]].superset_pos[supsetid[n1]]*ha_var[i].strides[n1];
                        }
                        ha_varele[ha_var[i].offset+recount1].value=val;
                        recount++;
                      }
                      readitem = strtok(NULL,"\n");
                      val=atof(readitem);//sscanf(readitem, "%lf", &val);
                      l1=recount;
                      if (ha_var[i].size>1) {
                        for (n1=ha_var[i].size-1; n1>1; n1--) {
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
                      for (n1=0; n1<ha_var[i].size; n1++) {
                        recount1=recount1+ha_setele[offset[n1]+index[n1]].superset_pos[supsetid[n1]]*ha_var[i].strides[n1];
                      }
                      ha_varele[ha_var[i].offset+recount1].value=val;
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
            printf("Error!!!! %s is not a header in data file!\n",header);
            return -1;
          }
          fclose(filehandle1);
          break;
        }
      }
      if(count2==0&&i==ncof){
        printf("Error!!!! %s is not a variable coefficient or parameter!\n",vname);
        return -1;
      }
    }
  }
  fclose(filehandle);
  return j;
}


offset_t sum_parse(char *formulain, char *commsyntax, sum_def *sum_cof,quantifier *arSet,set_def *ha_set,dim_t nset,dim_t fdim,int j) {
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
          strcat(interchar,"ha_cgesum000");
        }
        if (9<j&&j<100) {
          strcat(interchar,"ha_cgesum00");
        }
        if (99<j&&j<1000) {
          strcat(interchar,"ha_cgesum0");
        }
        if (999<j&&j<10000) {
          strcat(interchar,"ha_cgesum");
        }
        if (i>10000) {
          strcat(interchar,"ha_cgesum");
          printf("Error: Too many sum\n");
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
        for (l7=0; l7<nset; l7++) if(strcmp(p,ha_set[l7].setname)==0) {
            sum_cof[j].sumsetid=l7;
            break;
          }

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
            if(strcmp(p,sum_cof[j].sumindx)!=0) {
              strcat(interchar,sum_cof[j].sumindx);
              for (l4=0; l4<l3; l4++) {
                if(strcmp(p,sum_cof[j].dimnames[l4])==0) {
                  break;
                }
              }
              if (l4==l3||l3==0) {
                strcpy(sum_cof[j].dimnames[l3],p);
                l6=0;
                for (l5=0; l5<fdim-1; l5++) if(strcmp(p,arSet[l5].index_name)==0) {
                    sum_cof[j].setid[l3]=arSet[l5].setid;
                    l6++;
                  }
                if (l6==0) {
                  interchar1[0]='\0';
                  strcat(interchar1,"sum(");
                  strcat(interchar1,p);
                  strcpy(line3,formulain);
                  line3[readitem-formulain]='\0';
                  l7=str_rfind_ci(line3,interchar1);
                  p1=&line3[l7+2];
                  p1 = strtok(p1,",");
                  for (l7=0; l7<nset; l7++) if(strcmp(p1,ha_set[l7].setname)==0) {
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
              if(strcmp(p,sum_cof[j].sumindx)!=0) {
                for (l4=0; l4<l3; l4++) {
                  if(strcmp(p,sum_cof[j].dimnames[l4])==0) {
                    break;
                  }
                }
                if (l4==l3||l3==0) {
                  strcpy(sum_cof[j].dimnames[l3],p);
                  strcat(interchar,sum_cof[j].dimnames[l3]);
                  l6=0;
                  for (l5=0; l5<fdim-1; l5++) if(strcmp(p,arSet[l5].index_name)==0) {
                      sum_cof[j].setid[l3]=arSet[l5].setid;
                      l6++;
                      break;
                    }
                  if (l6==0) {
                    interchar1[0]='\0';
                    strcat(interchar1,"sum(");
                    strcat(interchar1,p);
                    strcpy(line3,formulain);
                    line3[readitem-formulain]='\0';
                    l7=str_rfind_ci(line3,interchar1);
                    p1=&line3[l7+2];
                    p2=strchr(p1,',');
                    strncpy(tempname,p1,p2-p1);
                    tempname[p2-p1]='\0';
                    for (l7=0; l7<nset; l7++) if(strcmp(tempname,ha_set[l7].setname)==0) {
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
          strcat(interchar,"ha_cgesum000");
        }
        if (9<j&&j<100) {
          strcat(interchar,"ha_cgesum00");
        }
        if (99<j&&j<1000) {
          strcat(interchar,"ha_cgesum0");
        }
        if (999<j&&j<10000) {
          strcat(interchar,"ha_cgesum");
        }
        if (i>10000) {
          strcat(interchar,"ha_cgesum");
          printf("Error: Too many sum\n");
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
        for (l7=0; l7<nset; l7++) if(strcmp(p,ha_set[l7].setname)==0) {
            sum_cof[j].sumsetid=l7;
            break;
          }

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
            if(strcmp(p,sum_cof[j].sumindx)!=0) {
              for (l4=0; l4<l3; l4++) {
                if(strcmp(p,sum_cof[j].dimnames[l4])==0) {
                  break;
                }
              }
              if (l4==l3||l3==0) {
                strcpy(sum_cof[j].dimnames[l3],p);
                strcat(interchar,p);
                l6=0;
                for (l5=0; l5<fdim-1; l5++) if(strcmp(p,arSet[l5].index_name)==0) {
                    sum_cof[j].setid[l3]=arSet[l5].setid;
                    l6++;
                  }
                if (l6==0) {
                  interchar1[0]='\0';
                  strcat(interchar1,"sum(");
                  strcat(interchar1,p);
                  strcpy(line3,formulain);
                  line3[readitem-formulain]='\0';
                  l7=str_rfind_ci(line3,interchar1);
                  p1=&line3[l7+2];
                  p1 = strtok(p1,",");
                  for (l7=0; l7<nset; l7++) if(strcmp(p1,ha_set[l7].setname)==0) {
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
              if(strcmp(p,sum_cof[j].sumindx)!=0) {
                for (l4=0; l4<l3; l4++) {
                  if(strcmp(p,sum_cof[j].dimnames[l4])==0) {
                    break;
                  }
                }
                if (l4==l3||l3==0) {
                  strcpy(sum_cof[j].dimnames[l3],p);
                  strcat(interchar,sum_cof[j].dimnames[l3]);
                  strcat(interchar,",");
                  l6=0;
                  for (l5=0; l5<fdim-1; l5++) if(strcmp(p,arSet[l5].index_name)==0) {
                      sum_cof[j].setid[l3]=arSet[l5].setid;
                      l6++;
                    }
                  if (l6==0) {
                    interchar1[0]='\0';
                    strcat(interchar1,"sum(");
                    strcat(interchar1,p);
                    strcpy(line3,formulain);
                    line3[readitem-formulain]='\0';
                    l7=str_rfind_ci(line3,interchar1);
                    p1=&line3[l7+2];
                    p2=strchr(p1,',');
                    strncpy(tempname,p1,p2-p1);
                    tempname[p2-p1]='\0';
                    for (l7=0; l7<nset; l7++) if(strcmp(tempname,ha_set[l7].setname)==0) {
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

offset_t tab_count_statements(char *fname, char *commsyntax) {
  FILE * filehandle;
  char line[TABREADLINE]="\0";
  offset_t j=0;
  filehandle = fopen(fname,"r");
  while (tab_next_statement(commsyntax,filehandle,line,TABREADLINE)) {
    if (strstr(line,"(default")==NULL) {
      j++;
    }
  }
  fclose(filehandle);
  return j;
}

offset_t closure_read(char *fname, char *commsyntax,closure_entry *ha_cgeshock, array_def *ha_var,offset_t nvar,set_def *ha_set,dim_t nset, set_element *ha_setele) {
  FILE * filehandle;
  char line[TABREADLINE]="\0",*readitem=NULL,*p=NULL,*p1=NULL,vname[TABREADLINE],argu[TABREADLINE];//,linecopy[TABREADLINE]
  offset_t i,j,n=0,l=0,l1,l2,n1,m,dims,dcount,exodims,supsetid[MAXSUPSET],sup,doublepr[MAXSUPSET];
  int k1,k2;
  bool check;
  filehandle = fopen(fname,"r");
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
    readitem++;
    for (i=0; i<n1-1; i++) {
      p = strchr(readitem,' ');
      strncpy(vname,readitem,p-readitem);
      vname[p-readitem]='\0';
      readitem=p;
      readitem++;
      check=true;
      if (strchr(vname,'(')==NULL) {
        for (j=0; j<nvar; j++) if (strcmp(vname,ha_var[j].cofname)==0) {
            dims=ha_var[j].nelem;
            if (dims==0) {
              ha_cgeshock[ha_var[j].offset].is_exogenous=true;
              n=n+1;
            } else for (l=0; l<dims; l++) {
                ha_cgeshock[ha_var[j].offset+l].is_exogenous=true;
                n=n+1;
              }
            check=false;
            break;
          }
          if(j==nvar&&vname[0]!=';')printf("Error!!! There is no such var %s\n",vname);
      } else {
        dims=1;
        p = strtok(vname,"(");
        for(sup=0; sup<MAXSUPSET; sup++)supsetid[sup]=0;
        for (j=0; j<nvar; j++) if (strcmp(p,ha_var[j].cofname)==0) {
            quantifier *arSet= (quantifier *) calloc (ha_var[j].size,sizeof(quantifier));
            offset_t *exoantidim= (offset_t *) calloc (ha_var[j].size,sizeof(offset_t));
            switch (ha_var[j].size) {
            case 1:
              p = strtok(NULL,")");//p = strtok(vname,")");
              if (strchr(p,'"')!=NULL) {
                strcpy(argu,p);
                p=strchr(argu,'"');
                p+=1;
                p1=strchr(p,'"');
                *p1='\0';
                for (l1=0; l1<ha_set[ha_var[j].setid[0]].size; l1++)
                  if (strcmp(p,ha_setele[ha_set[ha_var[j].setid[0]].offset+l1].setele)==0) {
                    ha_cgeshock[ha_var[j].offset+l1].is_exogenous=true;
                    n=n+1;
                    check=false;
                    break;
                  }
                  if(l1==ha_set[ha_var[j].setid[0]].size)printf("Error!!! There is no element %s in the set %s in cof/var %s\n",p,ha_set[ha_var[j].setid[0]].setname,ha_var[j].cofname);
              } else {
                if(strcmp(p,ha_set[ha_var[j].setid[0]].setname)==0) {
                  arSet[0].setid=ha_var[j].setid[0];
                  supsetid[0]=0;
                  exoantidim[0]=1;//arSet[0].SetSize;
                  dims=ha_set[arSet[0].setid].size;
                } else {
                  for (l1=0; l1<nset; l1++) if (strcmp(p,ha_set[l1].setname)==0) {
                      arSet[0].setid=l1;
                      for(sup=1; sup<MAXSUPSET; sup++)if(ha_set[l1].subsetid[sup]==ha_var[j].setid[0])supsetid[0]=sup;
                      exoantidim[0]=1;//arSet[0].SetSize;
                      dims=ha_set[arSet[0].setid].size;
                      break;
                    }
                    if(l1==nset)printf("Error!!! There is no set %s in cof/var %s\n",p,ha_var[j].cofname);
                }
              }
              break;
            default:
              p = strtok(NULL,")");
              strcat(p,",");
              for (l=0; l<ha_var[j].size; l++) {
                if (l==0) {
                  p = strtok(p,",");
                } else {
                  p = strtok(NULL,",");
                }
                if (*p=='"') {
                  strcpy(argu,p+1);
                  str_replace_char(argu,'"','\0');
                  for(l1=0; l1<ha_set[ha_var[j].setid[l]].size; l1++) if (strcmp(argu,ha_setele[ha_set[ha_var[j].setid[l]].offset+l1].setele)==0) {
                      doublepr[l]=l1;
                      supsetid[l]=-2;
                      break;
                    }
                    if(l1==ha_set[ha_var[j].setid[l]].size)printf("Error!!! There is no element %s in the set %s in cof/var %s\n",argu,ha_set[ha_var[j].setid[l]].setname,ha_var[j].cofname);
                } else {
                  if(strcmp(p,ha_set[ha_var[j].setid[l]].setname)==0) {
                    arSet[l].setid=ha_var[j].setid[l];
                    dims=dims*ha_set[arSet[l].setid].size;
                  } else {
                    for(l1=0; l1<nset; l1++) if (strcmp(p,ha_set[l1].setname)==0) {
                        arSet[l].setid=l1;
                        for(sup=1; sup<MAXSUPSET; sup++)if(ha_set[l1].subsetid[sup]==ha_var[j].setid[l])supsetid[l]=sup;
                        dims=dims*ha_set[arSet[l].setid].size;
                        break;
                      }
                      if(l1==nset)printf("Error!!! There is no set %s in cof/var %s\n",p,ha_var[j].cofname);
                  }
                }
              }
              exoantidim[ha_var[j].size-1]=1;
              for (l1=ha_var[j].size-2; l1>-1; l1--) {
                if(supsetid[l1+1]>-2)exoantidim[l1]=exoantidim[l1+1]*ha_set[arSet[l1+1].setid].size;
                else exoantidim[l1]=exoantidim[l1+1];
              }
              break;
            }
            if(check) for (l=0; l<dims; l++) {
                m=l;
                exodims=0;
                for (dcount=0; dcount<ha_var[j].size; dcount++) {
                  l1=(offset_t) m/exoantidim[dcount];
                  m=m-l1*exoantidim[dcount];
                  if(supsetid[dcount]==-2) {
                    l2=doublepr[dcount];//ha_set[arSet[dcount].setid].begadd;  //Special case
                  } else {
                    l2=ha_setele[ha_set[arSet[dcount].setid].offset+l1].superset_pos[supsetid[dcount]];
                  }
                  exodims=exodims+l2*ha_var[j].strides[dcount];
                }
                ha_cgeshock[ha_var[j].offset+exodims].is_exogenous=true;
                n=n+1;
              }
            free(arSet);
            free(exoantidim);
            break;
          }
          if(j==nvar)printf("Error!!! There is no such var %s\n",p);
      }
    }
  }
  fclose(filehandle);
  return n;
}

offset_t shocks_read(char *fname, char *commsyntax,closure_entry *ha_cgeshock,offset_t nvarele, array_def *ha_var,offset_t nvar,set_def *ha_set,dim_t nset, set_element *ha_setele,dim_t subints) {
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
        printf("Error opening file %s\n",fname);
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
      printf("Error!!!! The shock file line (from Shock to ;) is too long!\n");
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
      for (j=0; j<nvar; j++) if (strcmp(readitem,ha_var[j].cofname)==0) {
          readitem = strtok(NULL,";");
          k1=str_rfind_ci(readitem,"uniform");
          if(k1!=-1){
            ha_cgeshock[ha_var[j].offset].shock_value=atof(readitem+k1+1)/subints;
          }else{
            ha_cgeshock[ha_var[j].offset].shock_value=atof(readitem)/subints;
          }
          l=l+1;
          break;
        }
    } else {
      readitem = strtok(NULL,"(");
      for (j=0; j<nvar; j++) {
        if (strcmp(readitem,ha_var[j].cofname)==0) {
          dims=1;
          readitem = strtok(NULL,")");
          strcpy(argu,readitem);
          readitem = strtok(NULL,"=");
          readitem = strtok(NULL,";");
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
          for (n1=0; n1<ha_var[j].size; n1++) {
            if(n1==0) {
              p=strtok(argu,",");
            } else {
              p=strtok(NULL,",");
            }
            if(strchr(p,'\"')==NULL) {
              if((strcmp(ha_set[ha_var[j].setid[n1]].setname,p)==0)) {
                dimbegadd[n1]=ha_set[ha_var[j].setid[n1]].offset;//ha_set[k1].begadd;
                dimssize[n1]=ha_set[ha_var[j].setid[n1]].size;//ha_set[k1].size;
                dims=dims*dimssize[n1];
                supsetid[n1]=0;
              } else {
                for (k1=0; k1<nset; k1++)if (strcmp(ha_set[k1].setname,p)==0) {
                    dimbegadd[n1]=ha_set[k1].offset;
                    dimssize[n1]=ha_set[k1].size;
                    dims=dims*ha_set[k1].size;
                    for(sup=0; sup<MAXSUPSET; sup++)if(ha_set[k1].subsetid[sup]==ha_var[j].setid[n1]) {
                        supsetid[n1]=sup;
                        break;
                      }
                    if(sup==MAXSUPSET){
                      printf("Error!!!! %s in variable %s in shock file is not a correct set at that position!\n",p,ha_var[j].cofname);
                      return -1;
                    }
                    break;
                  }
                  if(k1==nset){
                    printf("Error!!!! %s in variable %s in shock file is not a set!\n",p,ha_var[j].cofname);
                    return -1;
                  }
              }
            } else {
              k2=0;
              while (p[k2] != '\0') {
                k2++;
              }
              for (k1=0; k1<ha_set[ha_var[j].setid[n1]].size; k1++)if (strncmp(ha_setele[ha_set[ha_var[j].setid[n1]].offset+k1].setele,p+1,k2-2)==0&&ha_setele[ha_set[ha_var[j].setid[n1]].offset+k1].setele[k2-2]=='\0') {
                  dimindx[n1]=ha_setele[ha_set[ha_var[j].setid[n1]].offset+k1].superset_pos[0];
                  dimssize[n1]=1;
                  supsetid[n1]=-2;
                  break;
                }
              if(k1==ha_set[ha_var[j].setid[n1]].size){
                printf("Error!!!! %s in variable %s in shock file is not a set element!\n",p,ha_var[j].cofname);
                return -1;
              }
            }
          }
          if (ha_var[j].size>0)antidim[0]=dimssize[1];
          if (ha_var[j].size>1)antidim[1]=1;
          if (ha_var[j].size>2)antidim[2]=antidim[0]*dimssize[0];
          if (ha_var[j].size>3) {
            for (n1=3; n1<ha_var[j].size; n1++) {
              antidim[n1]=antidim[n1-1]*dimssize[n1-1];
            }
          }
          k1=str_rfind_ci(linecopy,"uniform");
          if(k1>-1) {
            val=atof(linecopy+k1+1);//sscanf(readitem, "%lf", &val);
            for (n1=0; n1<dims; n1++) {
              l2=n1;
              if(ha_var[j].size>1) {
                for (dcount=ha_var[j].size-1; dcount>1; dcount--) {
                  if(dimindx[dcount]>-1) {
                    arSet[dcount]=dimindx[dcount];//ha_setele[dimbegadd[dcount]+l1].setsh;
                  } else {
                    l1=(offset_t) l2/antidim[dcount];
                    if(supsetid[dcount]>-2)arSet[dcount]=ha_setele[dimbegadd[dcount]+l1].superset_pos[supsetid[dcount]];
                    else arSet[dcount]=dimindx[dcount];
                    l2=l2-l1*antidim[dcount];
                  }
                }
                if(dimindx[0]>-1) {
                  arSet[0]=dimindx[0];
                } else {
                  l1=(offset_t) l2/antidim[0];
                  if(supsetid[0]>-2)arSet[0]=ha_setele[dimbegadd[0]+l1].superset_pos[supsetid[0]];
                  else arSet[0]=ha_setele[dimbegadd[0]+l1].superset_pos[0];//dimindx[0];
                  l2=l2-l1*antidim[0];
                }
                if(dimindx[1]>-1) {
                  arSet[1]=dimindx[1];
                } else {
                  l1=(offset_t) l2/antidim[1];
                  if(supsetid[1]>-2)arSet[1]=ha_setele[dimbegadd[1]+l1].superset_pos[supsetid[1]];
                  else arSet[1]=ha_setele[dimbegadd[1]+l1].superset_pos[0];//l1;//dimindx[1];
                }
              } else {
                arSet[0]=l2;
              }
              l2=0;
              for (dcount=0; dcount<ha_var[j].size; dcount++) {
                l2=l2+arSet[dcount]*ha_var[j].strides[dcount];
              }
              ha_cgeshock[ha_var[j].offset+l2].shock_value=val/subints;
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
              val=atof(readitem);
              l2=n1;
              if(ha_var[j].size>1) {
                for (dcount=ha_var[j].size-1; dcount>1; dcount--) {
                  if(dimindx[dcount]>-1) {
                    arSet[dcount]=dimindx[dcount];//ha_setele[dimbegadd[dcount]+l1].setsh;
                  } else {
                    l1=(offset_t) l2/antidim[dcount];
                    if(supsetid[dcount]>-2)arSet[dcount]=ha_setele[dimbegadd[dcount]+l1].superset_pos[supsetid[dcount]];
                    else arSet[dcount]=dimindx[dcount];
                    l2=l2-l1*antidim[dcount];
                  }
                }
                if(ha_var[j].size==2){//compliance
                if(dimindx[1]>-1) {
                  arSet[1]=dimindx[1];
                } else {
                  l1=(offset_t) l2/dimssize[0];
                  if(supsetid[1]>-2)arSet[1]=ha_setele[dimbegadd[1]+l1].superset_pos[supsetid[1]];
                  else arSet[1]=ha_setele[dimbegadd[1]+l1].superset_pos[0];//dimindx[1];
                  l2=l2-l1*dimssize[0];
                }
                if(dimindx[0]>-1) {
                  arSet[0]=dimindx[0];
                } else {
                  l1=(offset_t) l2;
                  if(supsetid[0]>-2)arSet[0]=ha_setele[dimbegadd[0]+l1].superset_pos[supsetid[0]];
                  else arSet[0]=ha_setele[dimbegadd[0]+l1].superset_pos[0];//l1;//dimindx[0];
                }
                }else{
                if(dimindx[0]>-1) {
                  arSet[0]=dimindx[0];//ha_setele[dimbegadd[dcount]+l1].setsh;
                } else {
                  l1=(offset_t) l2/antidim[0];
                  if(supsetid[0]>-2)arSet[0]=ha_setele[dimbegadd[0]+l1].superset_pos[supsetid[0]];
                  else arSet[0]=ha_setele[dimbegadd[0]+l1].superset_pos[0];//dimindx[0];
                  l2=l2-l1*antidim[0];
                }
                if(dimindx[1]>-1) {
                  arSet[1]=dimindx[1];//ha_setele[dimbegadd[dcount]+l1].setsh;
                } else {
                  l1=(offset_t) l2/antidim[1];
                  if(supsetid[1]>-2)arSet[1]=ha_setele[dimbegadd[1]+l1].superset_pos[supsetid[1]];
                  else arSet[1]=ha_setele[dimbegadd[1]+l1].superset_pos[0];//dimindx[0];
                }
                }
              } else {
                if(supsetid[0]>-2)arSet[0]=ha_setele[dimbegadd[0]+l2].superset_pos[supsetid[0]];
                else arSet[0]=dimindx[0];
              }
              l2=0;
              for (dcount=0; dcount<ha_var[j].size; dcount++) {
                l2=l2+arSet[dcount]*ha_var[j].strides[dcount];
              }
              ha_cgeshock[ha_var[j].offset+l2].shock_value=val/subints;
              l=l+1;
            }
          }
          break;
        }
      }
      if(j==nvar){
        printf("Error!!!! %s in shock file is not a variable!\n",readitem);
        return -1;
      }
    }
  }
  fclose(filehandle);
  return l;
}

int variables_read_defaults(char *fname, array_def *record, offset_t ncof) {
  char line[TABREADLINE];
  FILE * filehandle;
  char commsyntax[]="variable";
  filehandle = fopen(fname,"r");
  offset_t i;
  while (tab_next_statement(commsyntax,filehandle,line,TABREADLINE)) {
    if(str_find_ci(line,"default")>0) {
      while (str_replace_all(line," ", ""));
      if(str_find_ci(line,"default=levels")>0) for (i=0; i<ncof; i++) {
          record[i].level_par=true;
        }
      if(str_find_ci(line,"default=change")>0) for (i=0; i<ncof; i++) {
          record[i].change_real=true;
        }
    }
  }
  fclose(filehandle);
  return 1;
}

offset_t variables_read(char *fname, char *commsyntax, array_def *record, offset_t ncof, set_def *ha_set,dim_t nset) {
  FILE * filehandle;
  char line[TABREADLINE]="\0",linecopy[TABREADLINE],setname1[NAMESIZE],setname[NAMESIZE],setname2[NAMESIZE],setname3[NAMESIZE],finditem[NAMESIZE],finditem1[NAMESIZE],finditem2[NAMESIZE],finditem3[NAMESIZE],vname[NAMESIZE];//,vnamecopy[NAMESIZE];
  char *tpnt=NULL;
  offset_t n,m,l,ncommsyntax=0,i,j=0,addi=0,add=0,orig;
  dim_t dcount;
  char *readitem=NULL;
  bool IsLevel=false,IsChange=false;
  while (commsyntax[ncommsyntax] != '\0') {
    ncommsyntax++;
  }
  filehandle = fopen(fname,"r");

  while (tab_next_statement(commsyntax,filehandle,line,TABREADLINE)) {
    if(strstr(line,"(default=levels)")!=NULL) {
      IsLevel=true;
    }
    if(strstr(line,"(default=change)")!=NULL) {
      IsChange=true;
    }
    if(IsLevel) {
      record[j].level_par=true;
    }
    if(IsChange) {
      record[j].change_real=true;
    }
    readitem=strstr(line, "(ge ");
    if(readitem!=NULL){
      record[j].gltype=1;
      strcpy(linecopy,line);
      tpnt=strstr(readitem, ")");
      tpnt+=1;
      memmove(readitem,tpnt,strlen(tpnt)+1);
      readitem=strstr(linecopy, "(ge ");
      readitem+=3;
      readitem = strtok(readitem,")");
      record[j].glval=atof(readitem);
    }
    readitem=strstr(line, "(gt ");
    if(readitem!=NULL){
      record[j].gltype=2;
      strcpy(linecopy,line);
      tpnt=strstr(readitem, ")");
      tpnt+=1;
      memmove(readitem,tpnt,strlen(tpnt)+1);
      readitem=strstr(linecopy, "(gt ");
      readitem+=3;
      readitem = strtok(readitem,")");
      record[j].glval=atof(readitem);
    }
    readitem=strstr(line, "(le ");
    if(readitem!=NULL){
      record[j].gltype=3;
      strcpy(linecopy,line);
      tpnt=strstr(readitem, ")");
      tpnt+=1;
      memmove(readitem,tpnt,strlen(tpnt)+1);
      readitem=strstr(linecopy, "(le ");
      readitem+=3;
      readitem = strtok(readitem,")");
      record[j].glval=atof(readitem);
    }
    readitem=strstr(line, "(lt ");
    if(readitem!=NULL){
      record[j].gltype=4;
      strcpy(linecopy,line);
      tpnt=strstr(readitem, ")");
      tpnt+=1;
      memmove(readitem,tpnt,strlen(tpnt)+1);
      readitem=strstr(linecopy, "(lt ");
      readitem+=3;
      readitem = strtok(readitem,")");
      record[j].glval=atof(readitem);
    }

    readitem=strstr(line, ",ge ");
    if(readitem!=NULL){
      record[j].gltype=1;
      strcpy(linecopy,line);
      tpnt=strstr(readitem, ")");
      memmove(readitem,tpnt,strlen(tpnt)+1);
      readitem=strstr(linecopy, ",ge ");
      readitem+=3;
      readitem = strtok(readitem,")");
      record[j].glval=atof(readitem);
      printf("line %s\n type %d val %lf read %s\n",linecopy,record[j].gltype,record[j].glval,readitem);
      printf("line %s\n type %d val %lf\n",line,record[j].gltype,record[j].glval);
    }
    readitem=strstr(line, ",gt ");
    if(readitem!=NULL){
      record[j].gltype=2;
      strcpy(linecopy,line);
      tpnt=strstr(readitem, ")");
      memmove(readitem,tpnt,strlen(tpnt)+1);
      readitem=strstr(linecopy, ",gt ");
      readitem+=3;
      readitem = strtok(readitem,")");
      record[j].glval=atof(readitem);
    }
    readitem=strstr(line, ",le ");
    if(readitem!=NULL){
      record[j].gltype=3;
      strcpy(linecopy,line);
      tpnt=strstr(readitem, ")");
      memmove(readitem,tpnt,strlen(tpnt)+1);
      readitem=strstr(linecopy, ",le ");
      readitem+=3;
      readitem = strtok(readitem,")");
      record[j].glval=atof(readitem);
    }
    readitem=strstr(line, ",lt ");
    if(readitem!=NULL){
      record[j].gltype=4;
      strcpy(linecopy,line);
      tpnt=strstr(readitem, ")");
      memmove(readitem,tpnt,strlen(tpnt)+1);
      readitem=strstr(linecopy, ",lt ");
      readitem+=3;
      readitem = strtok(readitem,")");
      record[j].glval=atof(readitem);
    }
    while (str_replace_all(line," ", ""));
    if(str_find_ci(line,"(change)")>0) {
      record[j].change_real=true;
      str_replace_first(line,"(change)", "");
    }
    if(str_find_ci(line,"(percent_change)")>0) {
      record[j].change_real=false;
      str_replace_first(line,"(percent_change)", "");
    }
    if(str_find_ci(line,"(levels)")>0) {
      record[j].level_par=true;
      str_replace_first(line,"(levels)", "");
    }
    if(str_find_ci(line,"(linear)")>0) {
      record[j].level_par=false;
      str_replace_first(line,"(linear)", "");
    }
    if(str_find_ci(line,"(linear,change)")>0) {
      record[j].change_real=true;
      record[j].level_par=false;
      str_replace_first(line,"(linear,change)", "");
    }
    if(str_find_ci(line,"(change,linear)")>0) {
      record[j].change_real=true;
      record[j].level_par=false;
      str_replace_first(line,"(change,linear)", "");
    }
    strcpy(linecopy,line);
    if (strstr(line,"(default")==NULL) {
      n=str_count_char(line,')');
      orig=0;
      if (strstr(line,"(orig_level")!=NULL) orig=1;
      if (n-orig>1) {
        for (i=1; i<n-orig; i++) {
          if (i==1) {
            readitem = strtok(line,",");
          } else {
            readitem = strtok(NULL,",");
          }
          readitem = strtok(NULL,",");
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
          while (str_replace_all(linecopy, finditem, setname));
          while (str_replace_all(linecopy, finditem1, setname1));
          while (str_replace_all(linecopy, finditem2, setname2));
          while (str_replace_all(linecopy, finditem3, setname3));
        }
        for (i=1; i<n; i++) {
          if (i==1) {
            readitem = strtok(linecopy,")");
          } else {
            readitem = strtok(NULL,")");
          }
        }
        readitem = strtok(NULL,";");
        readitem = strtok(readitem,"(");
        strcpy(record[j].cofname,readitem);
        readitem = strtok(NULL,")");
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
          for (l=0; l<nset; l++) {
            if (strcmp(ha_set[l].setname,readitem)==0) {
              add=add*ha_set[l].size;
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
          record[j].strides[l]=record[j].strides[l+1]*ha_set[record[j].setid[l+1]].size;
        }
        addi=addi+add;
      } else {
        if (n-orig==1) {
          printf("%s\n","Syntax error! Not enough parences");
        } else {
          if(orig==1) {
            strcpy(setname,line+ncommsyntax);
            str_replace_all(setname,";", "");
            str_replace_all(setname,"\n", "");
            readitem=strchr(setname,')');
            readitem++;
            strcpy(record[j].cofname,readitem);
            record[j].offset=addi;
            record[j].size=0;
            record[j].nelem=1;
            addi=addi+1;
          } else {
            strcpy(setname,line+ncommsyntax);
            str_replace_all(setname,";", "");
            str_replace_all(setname,"\n", "");
            strcpy(record[j].cofname,setname);
            record[j].offset=addi;
            record[j].size=0;
            record[j].nelem=1;
            addi=addi+1;
          }
        }
      }
      j++;
    }
  }
  fclose(filehandle);
  return addi;
}

offset_t set_find_alltime(set_def *ha_set,dim_t nset) {
  dim_t i;
  for(i=0; i<nset; i++) {
    if(ha_set[i].intertemp&&ha_set[i].subsetid[1]==-1) {
      break;
    }
  }
  if(i==nset) {
    return -1;
  } else {
    return i;
  }
}


offset_t coefficients_read(char *fname, char *commsyntax, array_def *record, offset_t ncof, set_def *ha_set,offset_t nset) {
  FILE * filehandle;
  char line[TABREADLINE]="\0",linecopy[TABREADLINE],setname1[TABREADLINE],setname[TABREADLINE],setname2[TABREADLINE],setname3[TABREADLINE],finditem[TABREADLINE],finditem1[TABREADLINE],finditem2[TABREADLINE],finditem3[TABREADLINE],vname[TABREADLINE];//,vnamecopy[NAMESIZE];
  offset_t n,m,l,ncommsyntax=0,i=0,j=0,addi=0,add=0;
  dim_t dcount;
  char *readitem=NULL,*tpnt=NULL;
  while (commsyntax[ncommsyntax] != '\0') {
    ncommsyntax++;
  }
  filehandle = fopen(fname,"r");

  while (tab_next_statement(commsyntax,filehandle,line,TABREADLINE)) {
    readitem=strstr(line, "(ge ");
    if(readitem!=NULL){
      record[j].gltype=1;
      strcpy(linecopy,line);
      tpnt=strstr(readitem, ")");
      tpnt+=1;
      memmove(readitem,tpnt,strlen(tpnt)+1);
      readitem=strstr(linecopy, "(ge ");
      readitem+=3;
      readitem = strtok(readitem,")");
      record[j].glval=atof(readitem);
    }
    readitem=strstr(line, "(gt ");
    if(readitem!=NULL){
      record[j].gltype=2;
      strcpy(linecopy,line);
      tpnt=strstr(readitem, ")");
      tpnt+=1;
      memmove(readitem,tpnt,strlen(tpnt)+1);
      readitem=strstr(linecopy, "(gt ");
      readitem+=3;
      readitem = strtok(readitem,")");
      record[j].glval=atof(readitem);
    }
    readitem=strstr(line, "(le ");
    if(readitem!=NULL){
      record[j].gltype=3;
      strcpy(linecopy,line);
      tpnt=strstr(readitem, ")");
      tpnt+=1;
      memmove(readitem,tpnt,strlen(tpnt)+1);
      readitem=strstr(linecopy, "(le ");
      readitem+=3;
      readitem = strtok(readitem,")");
      record[j].glval=atof(readitem);
    }
    readitem=strstr(line, "(lt ");
    if(readitem!=NULL){
      record[j].gltype=4;
      strcpy(linecopy,line);
      tpnt=strstr(readitem, ")");
      tpnt+=1;
      memmove(readitem,tpnt,strlen(tpnt)+1);
      readitem=strstr(linecopy, "(lt ");
      readitem+=3;
      readitem = strtok(readitem,")");
      record[j].glval=atof(readitem);
    }

    readitem=strstr(line, ",ge ");
    if(readitem!=NULL){
      record[j].gltype=1;
      strcpy(linecopy,line);
      tpnt=strstr(readitem, ")");
      memmove(readitem,tpnt,strlen(tpnt)+1);
      readitem=strstr(linecopy, ",ge ");
      readitem+=3;
      readitem = strtok(readitem,")");
      record[j].glval=atof(readitem);
    }
    readitem=strstr(line, ",gt ");
    if(readitem!=NULL){
      record[j].gltype=2;
      strcpy(linecopy,line);
      tpnt=strstr(readitem, ")");
      memmove(readitem,tpnt,strlen(tpnt)+1);
      readitem=strstr(linecopy, ",gt ");
      readitem+=3;
      readitem = strtok(readitem,")");
      record[j].glval=atof(readitem);
    }
    readitem=strstr(line, ",le ");
    if(readitem!=NULL){
      record[j].gltype=3;
      strcpy(linecopy,line);
      tpnt=strstr(readitem, ")");
      memmove(readitem,tpnt,strlen(tpnt)+1);
      readitem=strstr(linecopy, ",le ");
      readitem+=3;
      readitem = strtok(readitem,")");
      record[j].glval=atof(readitem);
    }
    readitem=strstr(line, ",lt ");
    if(readitem!=NULL){
      record[j].gltype=4;
      strcpy(linecopy,line);
      tpnt=strstr(readitem, ")");
      memmove(readitem,tpnt,strlen(tpnt)+1);
      readitem=strstr(linecopy, ",lt ");
      readitem+=3;
      readitem = strtok(readitem,")");
      record[j].glval=atof(readitem);
    }
    
    while (str_replace_all(line," ", ""));
    str_replace_first(line,"parameter", "");
    str_replace_first(line,"change", "");
    str_replace_first(line,"integer", "");
    str_replace_first(line,"()", "");
    str_replace_first(line,"(,)", "");
    strcpy(linecopy,line);
    if (strstr(line,"(default")==NULL) {
      n=str_count_char(line,')');
      if (n>1) {
        for (i=1; i<n; i++) {
          if (i==1) {
            readitem = strtok(line,",");
          } else {
            readitem = strtok(NULL,",");
          }
          readitem = strtok(NULL,",");
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
          while (str_replace_all(linecopy, finditem, setname));
          while (str_replace_all(linecopy, finditem1, setname1));
          while (str_replace_all(linecopy, finditem2, setname2));
          while (str_replace_all(linecopy, finditem3, setname3));
        }
        for (i=1; i<n; i++) {
          if (i==1) {
            readitem = strtok(linecopy,")");
          } else {
            readitem = strtok(NULL,")");
          }
        }
        readitem = strtok(NULL,";");
        readitem = strtok(readitem,"(");
        strcpy(record[j].cofname,readitem);
        if(record[j].cofname[0]=='c'&&record[j].cofname[1]=='_'){
          printf("Errors!!!! c_ and C_ are reserved for change variable!!!\nPlease change the name of coefficient %s\n",record[j].cofname);
          return -1;
        }
        readitem = strtok(NULL,")");
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
          for (l=0; l<nset; l++) {
            if (strcmp(ha_set[l].setname,readitem)==0) {
              add=add*ha_set[l].size;
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
          record[j].strides[l]=record[j].strides[l+1]*ha_set[record[j].setid[l+1]].size;
        }
        addi=addi+add;
      } else {
        if (n==1) {
          printf("%s\n","Syntax error! Not enough parences");
          printf("line %s\n",line);
        } else {
          strcpy(setname,line+ncommsyntax);
          str_replace_all(setname,";", "");
          str_replace_all(setname,"\n", "");
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


int coef_resolve_sets(array_def *ha_cof,offset_t ncof, set_def *ha_set,dim_t nset, elem_store *ha_cofele) {
  dim_t l,m,dcount;//,dcount1;//,l3=0,l4=0,ld2=0,ld3=0,ld4=0
  offset_t l1=0,i;//,sizeele,dcountdim1[MAXVARDIM],l2=0,dcountdim[MAXVARDIM];//,ld1=0
  for (i=0; i<ncof; i++) {
    for (dcount=0; dcount<ha_cof[i].size; dcount++) {
      for (l=0; l<nset; l++) {
        if (strcmp(ha_set[l].setname,ha_set[ha_cof[i].setid[dcount]].setname)==0) {
          ha_set[ha_cof[i].setid[dcount]].offset=ha_set[l].offset;
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
        if(strchr(readitem,'+')!=NULL||strchr(readitem,'-')!=NULL) {
          if (strchr(readitem,'+')!=NULL) {
            sign[0]=1;
          } else {
            sign[0]=2;
          }
          readitem = strtok(readitem,"-+");
          if(readitem[0]>='0'&&readitem[0]<='9') {
            intindx[0]=atoi(readitem);
          } else {
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
              strcpy(floginame,readitem);
              readitem = strtok(NULL,"\"");
              readitem = strtok(NULL,"\"");
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
          if(readitem[0]>='0'&&readitem[0]<='9') {
            intindx[1]=atoi(readitem);
          } else {
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
              strcpy(floginame,readitem);
              readitem = strtok(NULL,"\"");
              readitem = strtok(NULL,"\"");
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
          if(readitem[0]>='0'&&readitem[0]<='9') {
            intindx[0]=atoi(readitem);
          } else {
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
              strcpy(floginame,readitem);
              readitem = strtok(NULL,"\"");
              readitem = strtok(NULL,"\"");
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
        if(strchr(readitem,'+')!=NULL||strchr(readitem,'-')!=NULL) {
          if (strchr(readitem,'+')!=NULL) {
            sign[0]=1;
          } else {
            sign[0]=2;
          }
          readitem = strtok(readitem,"-+");
          if(readitem[0]>='0'&&readitem[0]<='9') {
            intindx[0]=atoi(readitem);
          } else {
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
              strcpy(floginame,readitem);
              readitem = strtok(NULL,"\"");
              readitem = strtok(NULL,"\"");
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
          if(readitem[0]>='0'&&readitem[0]<='9') {
            intindx[1]=atoi(readitem);
          } else {
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
              strcpy(floginame,readitem);
              readitem = strtok(NULL,"\"");
              readitem = strtok(NULL,"\"");
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
          if(readitem[0]>='0'&&readitem[0]<='9') {
            intindx[0]=atoi(readitem);
          } else {
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
              strcpy(floginame,readitem);
              readitem = strtok(NULL,"\"");
              readitem = strtok(NULL,"\"");
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
        if(strchr(readitem,'+')!=NULL||strchr(readitem,'-')!=NULL) {
          if (strchr(readitem,'+')!=NULL) {
            sign[1]=1;
          } else {
            sign[1]=2;
          }
          readitem = strtok(readitem,"-+");
          if(readitem[0]>='0'&&readitem[0]<='9') {
            intindx[2]=atoi(readitem);
          } else {
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
              strcpy(floginame,readitem);
              readitem = strtok(NULL,"\"");
              readitem = strtok(NULL,"\"");
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
          if(readitem[0]>='0'&&readitem[0]<='9') {
            intindx[3]=atoi(readitem);
          } else {
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
              strcpy(floginame,readitem);
              readitem = strtok(NULL,"\"");
              readitem = strtok(NULL,"\"");
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
          if(readitem[0]>='0'&&readitem[0]<='9') {
            intindx[2]=atoi(readitem);
          } else {
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
              strcpy(floginame,readitem);
              readitem = strtok(NULL,"\"");
              readitem = strtok(NULL,"\"");
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
        record[j].readele[0]='\0';
        for (i=intvar[0]; i<intvar[1]+1; i++) {
          sprintf(line, "%d",i);
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
  char line[TABREADLINE]="\0",linecopy[TABREADLINE],line1[TABREADLINE];
  char longname[TABREADLINE],tempvar[256];
  char *commsyntax="set";
  dim_t j=0,vsize,dim1,i;
  int k0,k1,k2;
  char *readitem=NULL;//,*p;

  filehandle = fopen(fname,"r");

  while (tab_next_statement(commsyntax,filehandle,line,TABREADLINE)) {
    strcpy(linecopy,line);
    k2=str_find_ci(line,"intertemporal");
    if(k2>-1) {
      readitem = strtok(line+k2," ");
      readitem = strtok(NULL," ");
      strcpy(record[j].setname,readitem);
      strcpy(line1,"intertemporal,");
      readitem = strtok(NULL,"[");
      readitem = strtok(NULL,"]");
      strcat(line1,readitem);
      readitem = strtok(NULL,";");
      if (strchr(readitem,'[')!=NULL) {
        strcat(line1,",");
        readitem = strtok(readitem,"[");
        readitem = strtok(NULL,"]");
        strcat(line1,readitem);
      }
      strcpy(record[j].readele,line1);
      record[j].intertemp=true;
    } else {
      if (strchr(line,'=')!=NULL) {
        str_replace_all(line,"union", "^");
        readitem = &line[4];
        strcpy(line1,readitem);
        while (str_replace_all(line1," ", ""));
        readitem = strtok(line1,"=");
        strcpy(record[j].setname,readitem);
        if (strchr(line,'-')!=NULL) {
          line1[0]='-';
          line1[1]=',';
          line1[2]='\0';
          readitem = strtok(NULL,"-");
          for (i=0; i<nset; i++) {
            if (strcmp(readitem,record[i].setname)==0) {
              dim1=record[i].size;
              strcat(line1,readitem);
              strcat(line1,",");
              break;
            }
          }
          readitem = strtok(NULL,";");
          for (i=0; i<nset; i++) {
            if (strcmp(readitem,record[i].setname)==0) {
              record[j].size=dim1-record[i].size;
              strcat(line1,readitem);
              break;
            }
          }
          strcpy(record[j].readele,line1);
        } else {
          if(strchr(line,'+')!=NULL) {
            line1[0]='+';
            line1[1]=',';
            line1[2]='\0';
            readitem = strtok(NULL,"+");
            for (i=0; i<nset; i++) {
              if (strcmp(readitem,record[i].setname)==0) {
                dim1=record[i].size;
                strcat(line1,readitem);
                strcat(line1,",");
                break;
              }
            }
            readitem = strtok(NULL,";");
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
              for (i=0; i<nset; i++) {
                if (strcmp(readitem,record[i].setname)==0) {
                  dim1=record[i].size;
                  strcat(line1,readitem);
                  strcat(line1,",");
                  break;
                }
              }
              readitem = strtok(NULL,";");
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
              for (i=0; i<nset; i++) {
                if (strcmp(readitem,record[i].setname)==0) {
                  record[j].size=record[i].size;
                  strcpy(record[j].header,record[i].header);
                  sprintf(tempvar, "%d",i);
                  strcat(line1,tempvar);
                  break;
                }
              }
              strcpy(record[j].readele,line1);
            }
          }
        }
      } else {
        k0=str_find_ci(line,"read elements from file");
        if(k0>-1) {
          k1=str_find_ci(line+k0+24," ");
          strncpy(line1,line+k0+24,k1);
          line1[k1]='\0';
          for (k0=0; k0<niodata; k0++) if (strcmp(line1,iodata[k0].logname)==0) {
              break;
            }
          record[j].fileid=k0;
        }
        readitem = strtok(line," ");
        readitem = strtok(NULL," ");
        strcpy(record[j].setname,readitem);
        readitem = strtok(NULL,"\"");
        readitem = strtok(NULL,"\"");
        if (readitem!=NULL) {
          strcpy(record[j].header,readitem);
          record[j].readele[0]='\0';
          datafile_read_header_info(readitem,iodata[record[j].fileid].filname,&vsize,longname,&dim1);
          record[j].size=dim1;
        } else {
          record[j].header[0]='\0';
          readitem = strtok(linecopy,"(");
          readitem = strtok(NULL,")");
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

dim_t set_union_named(set_element *ha_setele, set_def *ha_set,dim_t nset,dim_t i) {
  dim_t j,l,n,m,dim1=0,dim2=0,j1,sup1,sup2;
  char line[TABREADLINE],*readitem;
  strcpy(line,ha_set[i].readele);
  readitem = strtok(line,",");
  readitem = strtok(NULL,",");
  for (l=0; l<nset; l++) {
    if (strcmp(readitem,ha_set[l].setname)==0) {
      dim1=ha_set[l].size;
      for (sup1=1; sup1<MAXSUPSET; sup1++)if(ha_set[l].subsetid[sup1]==-1) {
          ha_set[l].subsetid[sup1]=i;
          break;
        }
      break;
    }
  }
  m=0;
  for (n=0; n<dim1; n++) {
    strcpy(ha_setele[ha_set[i].offset+m].setele,ha_setele[ha_set[l].offset+n].setele);
    ha_setele[ha_set[i].offset+m].superset_pos[0]=n;
    ha_setele[ha_set[l].offset+m].superset_pos[sup1]=n;
    m++;
  }
  readitem = strtok(NULL,",");
  for (j=0; j<nset; j++) {
    if (strcmp(readitem,ha_set[j].setname)==0) {
      dim2=ha_set[j].size;
      for (sup2=1; sup2<MAXSUPSET; sup2++)if(ha_set[j].subsetid[sup2]==-1) {
          ha_set[j].subsetid[sup2]=i;
          break;
        }
      break;
    }
  }
  for (n=0; n<dim2; n++) {
    for (j1=0; j1<dim1; j1++) if(strcmp(ha_setele[ha_set[l].offset+j1].setele,ha_setele[ha_set[j].offset+n].setele)==0) break;
    if(j1==dim1) {
      strcpy(ha_setele[ha_set[i].offset+m].setele,ha_setele[ha_set[j].offset+n].setele);
      ha_setele[ha_set[i].offset+m].superset_pos[0]=m;
      ha_setele[ha_set[j].offset+n].superset_pos[sup2]=m;
      m++;
    } else {
      ha_setele[ha_set[j].offset+n].superset_pos[sup2]=j1;
    }
  }
  ha_set[i].size=m;
  return m;
}
dim_t set_union_op(set_element *ha_setele, set_def *ha_set,dim_t nset,dim_t i) {
  dim_t j,l,n,m,dim1=0,dim2=0,j1,sup1,sup2;
  char line[TABREADLINE],*readitem;
  strcpy(line,ha_set[i].readele);
  readitem = strtok(line,",");
  readitem = strtok(NULL,",");
  for (l=0; l<nset; l++) {
    if (strcmp(readitem,ha_set[l].setname)==0) {
      dim1=ha_set[l].size;
      for (sup1=1; sup1<MAXSUPSET; sup1++)if(ha_set[l].subsetid[sup1]==-1) {
          ha_set[l].subsetid[sup1]=i;
          break;
        }
      break;
    }
  }
  m=0;
  for (n=0; n<dim1; n++) {
    strcpy(ha_setele[ha_set[i].offset+m].setele,ha_setele[ha_set[l].offset+n].setele);
    ha_setele[ha_set[i].offset+m].superset_pos[0]=n;
    ha_setele[ha_set[l].offset+m].superset_pos[sup1]=n;
    m++;
  }
  readitem = strtok(NULL,",");
  for (j=0; j<nset; j++) {
    if (strcmp(readitem,ha_set[j].setname)==0) {
      dim2=ha_set[j].size;
      for (sup2=1; sup2<MAXSUPSET; sup2++)if(ha_set[j].subsetid[sup2]==-1) {
          ha_set[j].subsetid[sup2]=i;
          break;
        }
      break;
    }
  }
  for (n=0; n<dim2; n++) {
    strcpy(ha_setele[ha_set[i].offset+m].setele,ha_setele[ha_set[j].offset+n].setele);
    ha_setele[ha_set[i].offset+m].superset_pos[0]=m;
    ha_setele[ha_set[j].offset+n].superset_pos[sup2]=m;
    m++;
  }
  return m;
}
dim_t set_difference(set_element *ha_setele, set_def *ha_set,dim_t nset,dim_t i) {
  dim_t j,l,n,m,v,indi,dim1=0,dim2=0,sup1;
  char line[TABREADLINE],*readitem;
  strcpy(line,ha_set[i].readele);
  readitem = strtok(line,",");
  readitem = strtok(NULL,",");
  for (l=0; l<nset; l++) {
    if (strcmp(readitem,ha_set[l].setname)==0) {
      dim1=ha_set[l].size;
      for (sup1=1; sup1<MAXSUPSET; sup1++)if(ha_set[i].subsetid[sup1]==-1) {
          ha_set[i].subsetid[sup1]=l;
          break;
        }
      break;
    }
  }
  readitem = strtok(NULL,",");
  for (j=0; j<nset; j++) {
    if (strcmp(readitem,ha_set[j].setname)==0) {
      dim2=ha_set[j].size;
      break;
    }
  }
  m=0;
  for (n=0; n<dim1; n++) {
    indi=0;
    for (v=0; v<dim2; v++) {
      if (strcmp(ha_setele[ha_set[l].offset+n].setele,ha_setele[ha_set[j].offset+v].setele)==0) {
        indi=1;
        break;
      }
    }
    if(indi==0) {
      strcpy(ha_setele[ha_set[i].offset+m].setele,ha_setele[ha_set[l].offset+n].setele);
      ha_setele[ha_set[i].offset+m].superset_pos[sup1]=n;
      ha_setele[ha_set[i].offset+m].superset_pos[0]=m;
      m++;
    }
  }
  return m;
}
offset_t subsets_read(char *fname, set_element *ha_setele, set_def *ha_set,dim_t nset) {
  FILE * filehandle;//, *fileout;
  char line[TABREADLINE]="\0";
  char set[NAMESIZE],subset[NAMESIZE];
  char *commsyntax="subset";
  dim_t i,setd,subsetd,sup1;//,nlength;
  offset_t jj,jjj,j=0,succ=0,ssize=0;
  char *readitem=NULL;

  filehandle = fopen(fname,"r");

  while (tab_next_statement(commsyntax,filehandle,line,TABREADLINE)) {
    while (str_replace_all(line," ;", ";"));
    readitem = strtok(line," ");
    readitem = strtok(NULL," ");
    strcpy(subset,readitem);
    readitem = strtok(NULL," ");
    readitem = strtok(NULL," ");
    readitem = strtok(NULL," ");
    readitem = strtok(NULL,";");
    strcpy(set,readitem);
    succ=0;
    for (i=0; i<nset; i++) {
      if (strcmp(ha_set[i].setname,subset)==0) { //,subsetd)==0) {
        for (j=0; j<nset; j++) {
          if (strcmp(ha_set[j].setname,set)==0) { //,setd)==0) {
            ssize=ha_set[j].size;
            for (sup1=1; sup1<MAXSUPSET; sup1++)if(ha_set[i].subsetid[sup1]==-1) {
                ha_set[i].subsetid[sup1]=j;
                break;
              }
            if(sup1==MAXSUPSET)printf("Error!!! Number of superset exceeded! Increase MAXSUPSET!\n");
            for (jj=ha_set[j].offset; jj<ha_set[j].offset+ha_set[j].size; jj++) {
              for (jjj=ha_set[i].offset; jjj<ha_set[i].offset+ha_set[i].size; jjj++) {
                if (strcmp(ha_setele[jj].setele,ha_setele[jjj].setele)==0) {
                  ha_setele[jjj].superset_pos[sup1]=ha_setele[jj].superset_pos[0];
                  succ++;
                  break;
                }
              }
            }
            break;
          }
        }
        if(j==nset)printf("Error!!! No such set %s\n",set);
        break;
      }
      j++;
    }
    if(i==nset)printf("Error!!! No such set %s\n",subset);
    if(succ-ha_set[i].size!=0)printf("Error!!! No set or some elements of set %s are not elements of set %s\n",subset,set);
  }
  fclose(filehandle);
  return j;
}

dim_t subset_map_build(set_element *ha_setele, set_def *ha_set,dim_t nset,offset_t* contin) {
  dim_t i,setd,subsetd,sup1,sup2,sup3,b;//,nlength;
  offset_t jj,jjj,j=0;
  *contin=0;
  for (i=0; i<nset; i++) {
    b=1;
    while(b==1) {
      for (sup1=1; sup1<MAXSUPSET; sup1++) {
        b=0;
        if(ha_set[i].subsetid[sup1]==-1) {
          b=2;
          break;
        }
        for (sup2=1; sup2<MAXSUPSET; sup2++) {
          b=0;
          if(ha_set[ha_set[i].subsetid[sup1]].subsetid[sup2]==-1) {
            b=2;
            break;
          }
          for (sup3=1; sup3<MAXSUPSET; sup3++){
            if(ha_set[ha_set[i].subsetid[sup1]].subsetid[sup2]==ha_set[i].subsetid[sup3]) {
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
        j=ha_set[ha_set[i].subsetid[sup1]].subsetid[sup2];
        for (sup1=1; sup1<MAXSUPSET; sup1++){
          if(ha_set[i].subsetid[sup1]==-1)break;
        }
        if(sup1==MAXSUPSET){
          printf("Error!!!! Superset size exceeded!!\n");
          *contin=0;
          return 0;
        }
        ha_set[i].subsetid[sup1]=j;
        for (jj=ha_set[j].offset; jj<ha_set[j].offset+ha_set[j].size; jj++) {
          for (jjj=ha_set[i].offset; jjj<ha_set[i].offset+ha_set[i].size; jjj++) {
            if (strcmp(ha_setele[jj].setele,ha_setele[jjj].setele)==0) {
              ha_setele[jjj].superset_pos[sup1]=ha_setele[jj].superset_pos[0];
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
        if(ha_set[i].subsetid[sup1]==-1) {
          b=2;
          break;
        }
        for (sup2=1; sup2<MAXSUPSET; sup2++) {
          b=0;
          if(ha_set[ha_set[i].subsetid[sup1]].subsetid[sup2]==-1) {
            b=2;
            break;
          }
          for (sup3=1; sup3<MAXSUPSET; sup3++){
            if(ha_set[ha_set[i].subsetid[sup1]].subsetid[sup2]==ha_set[i].subsetid[sup3]) {
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
        j=ha_set[ha_set[i].subsetid[sup1]].subsetid[sup2];
        for (sup1=1; sup1<MAXSUPSET; sup1++){
          if(ha_set[i].subsetid[sup1]==-1)break;
        }
        if(sup1==MAXSUPSET){
          printf("Error!!!! Superset size exceeded!!\n");
          *contin=0;
          return 0;
        }
        ha_set[i].subsetid[sup1]=j;
        for (jj=ha_set[j].offset; jj<ha_set[j].offset+ha_set[j].size; jj++) {
          for (jjj=ha_set[i].offset; jjj<ha_set[i].offset+ha_set[i].size; jjj++) {
            if (strcmp(ha_setele[jj].setele,ha_setele[jjj].setele)==0) {
              ha_setele[jjj].superset_pos[sup1]=ha_setele[jj].superset_pos[0];
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

char *tab_next_statement(char *commsyntax, FILE *filehandle, char *readline,offset_t rlinesize) {
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
          if(strlen(readline)>=((size_t)(9*rlinesize/10)))printf("Warning!!!! Overflow may happen in function ha_cgertabl. Increase TABREADLINE or DATREADLINE may help!%ld %ld\n",strlen(readline),rlinesize);
        return readline;
      }
      while (fgets(line,TABLINESIZE,filehandle)) {
        n=strstr(line,finditem);//ha_cgefendofc
        if (n==NULL) {
          strcat(readline, line);
        } else {
          strcat(readline, line);
          if(strlen(readline)>=((size_t)(9*rlinesize/10)))printf("Warning!!!! Overflow may happen in function ha_cgertabl. Increase TABREADLINE or DATREADLINE may help! %ld %ld\n",strlen(readline),rlinesize);
          return readline;
        }
      }

    }
  }
  return NULL;
}
char *tab_next_statement_resolved(char *commsyntax, FILE *filehandle, char *readline, elem_value *record, array_def *ha_cof,offset_t ncof,solve_real *zerodivide,offset_t rlinesize) {
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
          if(strlen(readline)>=((size_t)(9*rlinesize/10)))printf("Warning!!!! Overflow may happen in function ha_cgertabl1. Increase TABREADLINE or DATREADLINE may help! %ld %ld\n",strlen(readline),rlinesize);
        return readline;
      }
      while (fgets(line,TABLINESIZE,filehandle)) {
        n=strstr(line,finditem);//ha_cgefendofc
        if (n==NULL) {
          strcat(readline, line);
        } else {
          strcat(readline, line);
          if(strlen(readline)>=((size_t)(9*rlinesize/10)))printf("Warning!!!! Overflow may happen in function ha_cgertabl1. Increase TABREADLINE or DATREADLINE may help! %ld %ld\n",strlen(readline),rlinesize);
          return readline;
        }
      }

    }
    if (strncmp(line,commsyntax,count1) != 0) {
      if (strncmp(line,zerosyntax,18) == 0) {
        p=strtok(&line[18],";");
        *zerodivide=formula_subst_scalar(p,record,ha_cof,ncof);//atof(p);
      }
      if (strncmp(line,zerosyntax1,12) == 0||strncmp(line,zerosyntax2,11) == 0) {
        p=strtok(line,")");
        p=strtok(NULL,"default");
        p=strtok(NULL," ");
        p=strtok(NULL,";");
        *zerodivide=formula_subst_scalar(p,record,ha_cof,ncof);//atof(p);
      }
    }
  }
  return NULL;
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

