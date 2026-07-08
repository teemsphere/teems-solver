#include <teems_solver.h>

solve_real formula_subst_scalar(char *var2, elem_value *record, array_def *coefs,offset_t ncof) {
  offset_t index;
  while (str_replace_all(var2," ", ""));
  solve_real eval=0;
  if (var2[0]>='0'&&var2[0]<='9') {
    eval=atof(var2);
    return eval;
  }
  index=ncof-1;
  do {
    if (strcmp(coefs[index].cofname,var2)==0) {
      eval=record[coefs[index].offset].value;
      break;
    }
  } while (index--);
  return eval;
}


int formula_bind_operand(char *var2, set_def *sets,array_def *coefs,offset_t ncof, array_def *vars,offset_t nvar,offset_t ncofele,sum_def *sum_cof,int totalsum,formula_op *ops,int nops,quantifier *arSet,dim_t fdim,int varindex) {
  offset_t index;
  char *p=NULL;//,copyvar[TABREADLINE];//,*p1=NULL,*p2=NULL,*p3=NULL,*p4=NULL;
  dim_t l1,l,sup;//,svar2;//=0,i2=0,i3=0,i4=0,svar1,svar2,checkvar20=0,checkvar10=0,checkvar11=0,checkvar12=0,checkvar16=0,checkvar14=0,l;
  int leadlag;
  bool IsChange=false;
  p= strtok(var2,"{");
  if (p==NULL) {
    p=&var2[0];
  }
  if(var2[0]=='p'&&var2[1]=='_') {
    IsChange=true;
    p++;
    p++;
  }
  index=ncof-1;
  do {
    if (strcmp(coefs[index].cofname,p)==0) {
      if(!coefs[index].suplval)printf("Warning!!!! coefficient %s has not been supplied with values!\n",coefs[index].cofname);
      if(varindex==2) {
        ops[nops].Var2BegAdd=coefs[index].offset;
        ops[nops].Var2Type=OT_ARRAY;
      } else {
        ops[nops].Var1BegAdd=coefs[index].offset;
        ops[nops].Var1Type=OT_ARRAY;
      }
      switch(coefs[index].size) {
      case 0:
        for (l=0; l<fdim; l++) {
          if(varindex==2) {
            ops[nops].Var2Dims[l].ADims=0;
          } else {
            ops[nops].Var1Dims[l].ADims=0;
          }
        }
        break;
      case 1:
        p=strtok(NULL,"}");
        leadlag=0;
        parse_index_leadlag(p,&leadlag);
       for (l=0; l<fdim; l++) {
             if(varindex==2) {
              ops[nops].Var2Dims[l].ADims=0;
              ops[nops].Var2Dims[l].leadlag=0;
              ops[nops].Var2Dims[l].SupSet=0;
            } else {
              ops[nops].Var1Dims[l].ADims=0;
              ops[nops].Var1Dims[l].leadlag=0;
              ops[nops].Var1Dims[l].SupSet=0;
            }
          if (strcmp(p,arSet[l].index_name)==0) {
            if(varindex==2) {
              if (sets[coefs[index].setid[0]].size>sets[arSet[l].setid].size) {
                ops[nops].Var2Dims[l].SupSet=1;
                for(sup=1;sup<MAXSUPSET;sup++)if(sets[arSet[l].setid].subsetid[sup]==coefs[index].setid[0]){ops[nops].Var2Dims[l].SSIndx=sup;break;}
              }
              ops[nops].Var2Dims[l].ADims=coefs[index].strides[0];
              ops[nops].Var2Dims[l].leadlag=leadlag;
            } else {
              if (sets[coefs[index].setid[0]].size>sets[arSet[l].setid].size) {
                ops[nops].Var1Dims[l].SupSet=1;
                for(sup=1;sup<MAXSUPSET;sup++)if(sets[arSet[l].setid].subsetid[sup]==coefs[index].setid[0]){ops[nops].Var1Dims[l].SSIndx=sup;break;}
              }
              ops[nops].Var1Dims[l].ADims=coefs[index].strides[0];
              ops[nops].Var1Dims[l].leadlag=leadlag;
            }
          }
        }
        break;
      default:
        for (l1=0; l1<fdim; l1++) {
          if(varindex==2) {
            ops[nops].Var2Dims[l1].ADims=0;
            ops[nops].Var2Dims[l1].leadlag=0;
            ops[nops].Var2Dims[l1].SupSet=0;
          } else {
            ops[nops].Var1Dims[l1].ADims=0;
            ops[nops].Var1Dims[l1].leadlag=0;
            ops[nops].Var1Dims[l1].SupSet=0;
          }
        }
        for (l=0; l<coefs[index].size-1; l++) {
          p=strtok(NULL,",");
          leadlag=0;
          parse_index_leadlag(p,&leadlag);
          for (l1=0; l1<fdim; l1++) {
            if (strcmp(p,arSet[l1].index_name)==0) {
              if(varindex==2) {
                if (sets[coefs[index].setid[l]].size>sets[arSet[l1].setid].size) {
                  ops[nops].Var2Dims[l1].SupSet=1;
                  for(sup=1;sup<MAXSUPSET;sup++)if(sets[arSet[l1].setid].subsetid[sup]==coefs[index].setid[l]){ops[nops].Var2Dims[l1].SSIndx=sup;break;}
                }
                ops[nops].Var2Dims[l1].ADims=coefs[index].strides[l];
                ops[nops].Var2Dims[l1].leadlag=leadlag;
              } else {
                if (sets[coefs[index].setid[l]].size>sets[arSet[l1].setid].size) {
                  ops[nops].Var1Dims[l1].SupSet=1;
                  for(sup=1;sup<MAXSUPSET;sup++)if(sets[arSet[l1].setid].subsetid[sup]==coefs[index].setid[l]){ops[nops].Var1Dims[l1].SSIndx=sup;break;}
                }
                ops[nops].Var1Dims[l1].ADims=coefs[index].strides[l];
                ops[nops].Var1Dims[l1].leadlag=leadlag;
              }
              break;
            }
          }
        }
        p=strtok(NULL,"}");
        leadlag=0;
        parse_index_leadlag(p,&leadlag);
        for (l1=0; l1<fdim; l1++) {
          if (strcmp(p,arSet[l1].index_name)==0) {
            if(varindex==2) {
              if (sets[coefs[index].setid[l]].size>sets[arSet[l1].setid].size) {
                ops[nops].Var2Dims[l1].SupSet=1;
                for(sup=1;sup<MAXSUPSET;sup++)if(sets[arSet[l1].setid].subsetid[sup]==coefs[index].setid[l]){ops[nops].Var2Dims[l1].SSIndx=sup;break;}
              }
              ops[nops].Var2Dims[l1].ADims=coefs[index].strides[l];
              ops[nops].Var2Dims[l1].leadlag=leadlag;
            } else {
              if (sets[coefs[index].setid[l]].size>sets[arSet[l1].setid].size) {
                ops[nops].Var1Dims[l1].SupSet=1;
                for(sup=1;sup<MAXSUPSET;sup++)if(sets[arSet[l1].setid].subsetid[sup]==coefs[index].setid[l]){ops[nops].Var1Dims[l1].SSIndx=sup;break;}
              }
              ops[nops].Var1Dims[l1].ADims=coefs[index].strides[l];
              ops[nops].Var1Dims[l1].leadlag=leadlag;
            }
            break;
          }
        }
        break;
      }
      return 1;
    }
  } while (index--);

  index=nvar-1;
  do {
    if (strcmp(vars[index].cofname,p)==0) {
      if(varindex==2) {
        ops[nops].Var2BegAdd=vars[index].offset+ncofele;
        if(IsChange) ops[nops].Var2Type=OT_CHANGE;
        else ops[nops].Var2Type=OT_ARRAY;
      } else {
        ops[nops].Var1BegAdd=vars[index].offset+ncofele;
        if(IsChange) ops[nops].Var1Type=OT_CHANGE;
        else ops[nops].Var1Type=OT_ARRAY;
      }
      switch(vars[index].size) {
      case 0:
        for (l=0; l<fdim; l++) {
          if(varindex==2) {
            ops[nops].Var2Dims[l].ADims=0;
            ops[nops].Var2Dims[l].leadlag=0;
          } else {
            ops[nops].Var1Dims[l].ADims=0;
            ops[nops].Var1Dims[l].leadlag=0;
          }
        }
        break;
      case 1:
        p=strtok(NULL,"}");
        leadlag=0;
        parse_index_leadlag(p,&leadlag);
        for (l=0; l<fdim; l++) {
             if(varindex==2) {
              ops[nops].Var2Dims[l].ADims=0;
              ops[nops].Var2Dims[l].leadlag=0;
              ops[nops].Var2Dims[l].SupSet=0;
            } else {
              ops[nops].Var1Dims[l].ADims=0;
              ops[nops].Var1Dims[l].leadlag=0;
              ops[nops].Var1Dims[l].SupSet=0;
            }
          if (strcmp(p,arSet[l].index_name)==0) {
            if(varindex==2) {
              if (sets[vars[index].setid[0]].size>sets[arSet[l].setid].size) {
                ops[nops].Var2Dims[l].SupSet=1;
                for(sup=1;sup<MAXSUPSET;sup++)if(sets[arSet[l].setid].subsetid[sup]==vars[index].setid[0]){ops[nops].Var2Dims[l].SSIndx=sup;break;}
              }
              ops[nops].Var2Dims[l].ADims=vars[index].strides[0];
              ops[nops].Var2Dims[l].leadlag=leadlag;
            } else {
              if (sets[vars[index].setid[0]].size>sets[arSet[l].setid].size) {
                ops[nops].Var1Dims[l].SupSet=1;
                for(sup=1;sup<MAXSUPSET;sup++)if(sets[arSet[l].setid].subsetid[sup]==vars[index].setid[0]){ops[nops].Var1Dims[l].SSIndx=sup;break;}
              }
              ops[nops].Var1Dims[l].ADims=vars[index].strides[0];
              ops[nops].Var1Dims[l].leadlag=leadlag;
            }
          }
        }
        break;
      default:
        for (l1=0; l1<fdim; l1++) {
          if(varindex==2) {
            ops[nops].Var2Dims[l1].ADims=0;
            ops[nops].Var2Dims[l1].leadlag=0;
            ops[nops].Var2Dims[l1].SupSet=0;
          } else {
            ops[nops].Var1Dims[l1].ADims=0;
            ops[nops].Var1Dims[l1].leadlag=0;
            ops[nops].Var1Dims[l1].SupSet=0;
          }
        }
        for (l=0; l<vars[index].size-1; l++) {
          p=strtok(NULL,",");
          leadlag=0;
          parse_index_leadlag(p,&leadlag);
          for (l1=0; l1<fdim; l1++) {
            if (strcmp(p,arSet[l1].index_name)==0) {
              if(varindex==2) {
                if (sets[vars[index].setid[l]].size>sets[arSet[l1].setid].size) {
                  ops[nops].Var2Dims[l1].SupSet=1;
                  for(sup=1;sup<MAXSUPSET;sup++)if(sets[arSet[l1].setid].subsetid[sup]==vars[index].setid[l]){ops[nops].Var2Dims[l1].SSIndx=sup;break;}
                }
                ops[nops].Var2Dims[l1].ADims=vars[index].strides[l];
                ops[nops].Var2Dims[l1].leadlag=leadlag;
              } else {
                if (sets[vars[index].setid[l]].size>sets[arSet[l1].setid].size) {
                  ops[nops].Var1Dims[l1].SupSet=1;
                  for(sup=1;sup<MAXSUPSET;sup++)if(sets[arSet[l1].setid].subsetid[sup]==vars[index].setid[l]){ops[nops].Var1Dims[l1].SSIndx=sup;break;}
                }
                ops[nops].Var1Dims[l1].ADims=vars[index].strides[l];
                ops[nops].Var1Dims[l1].leadlag=leadlag;
              }
              break;
            }
          }
        }
        p=strtok(NULL,"}");
        leadlag=0;
        parse_index_leadlag(p,&leadlag);
        for (l1=0; l1<fdim; l1++) {
          if (strcmp(p,arSet[l1].index_name)==0) {
            if(varindex==2) {
              if (sets[vars[index].setid[l]].size>sets[arSet[l1].setid].size) {
                ops[nops].Var2Dims[l1].SupSet=1;
                for(sup=1;sup<MAXSUPSET;sup++)if(sets[arSet[l1].setid].subsetid[sup]==vars[index].setid[l]){ops[nops].Var2Dims[l1].SSIndx=sup;break;}
              }
              ops[nops].Var2Dims[l1].ADims=vars[index].strides[l];
              ops[nops].Var2Dims[l1].leadlag=leadlag;
            } else {
              if (sets[vars[index].setid[l]].size>sets[arSet[l1].setid].size) {
                ops[nops].Var1Dims[l1].SupSet=1;
                for(sup=1;sup<MAXSUPSET;sup++)if(sets[arSet[l1].setid].subsetid[sup]==vars[index].setid[l]){ops[nops].Var1Dims[l1].SSIndx=sup;break;}
              }
              ops[nops].Var1Dims[l1].ADims=vars[index].strides[l];
              ops[nops].Var1Dims[l1].leadlag=leadlag;
            }
            break;
          }
        }
        break;
      }
      return 1;
    }
  } while (index--);

  for (index=totalsum-1; index>-1; index--) {
    if (strcmp(sum_cof[index].sumname,p)==0) {
      if(varindex==2) {
        ops[nops].Var2BegAdd=sum_cof[index].offset;
        ops[nops].Var2Type=OT_SUM;
        for (l1=0; l1<sum_cof[index].size; l1++)ops[nops].Var2Dims[l1].SSIndx=0;
      } else {
        ops[nops].Var1BegAdd=sum_cof[index].offset;
        ops[nops].Var1Type=OT_SUM;
        for (l1=0; l1<sum_cof[index].size; l1++)ops[nops].Var1Dims[l1].SSIndx=0;
      }
      switch(sum_cof[index].size) {
      case 0:
        for (l1=0; l1<fdim; l1++) {
          if(varindex==2) ops[nops].Var2Dims[l1].ADims=0;
          else ops[nops].Var1Dims[l1].ADims=0;
        }
        break;
      case 1:
        p=strtok(NULL,"}");
        for (l=0; l<fdim; l++) {
          if (strcmp(p,arSet[l].index_name)==0) {
            if(varindex==2) {
              ops[nops].Var2Dims[l].ADims=sum_cof[index].strides[0];
            } else {
              ops[nops].Var1Dims[l].ADims=sum_cof[index].strides[0];
            }
          } else {
            if(varindex==2) {
              ops[nops].Var2Dims[l].ADims=0;
            } else {
              ops[nops].Var1Dims[l].ADims=0;
            }
          }
        }
        break;
      default:
        for (l1=0; l1<fdim; l1++) {
          if(varindex==2) ops[nops].Var2Dims[l1].ADims=0;
          else ops[nops].Var1Dims[l1].ADims=0;
        }
        for (l=0; l<sum_cof[index].size-1; l++) {
          p=strtok(NULL,",");
          for (l1=0; l1<fdim; l1++) {
            if (strcmp(p,arSet[l1].index_name)==0) {
              if(varindex==2) {
                ops[nops].Var2Dims[l1].ADims=sum_cof[index].strides[l];
              } else {
                ops[nops].Var1Dims[l1].ADims=sum_cof[index].strides[l];
              }
              break;
            }
          }
        }
        p=strtok(NULL,"}");
        for (l1=0; l1<fdim; l1++) {
          if (strcmp(p,arSet[l1].index_name)==0) {
            if(varindex==2) {
              ops[nops].Var2Dims[l1].ADims=sum_cof[index].strides[l];
            } else {
              ops[nops].Var1Dims[l1].ADims=sum_cof[index].strides[l];
            }
            break;
          }
        }
        break;
      }
      return 1;
    }
  }

  for (l1=0; l1<nops; l1++) {
    if (strcmp(var2,ops[l1].TmpVarName)==0) {
      if(varindex==2) {
        ops[nops].Var2Type=OT_TEMP;
        ops[nops].Var2BegAdd=l1;
      } else {
        ops[nops].Var1Type=OT_TEMP;
        ops[nops].Var1BegAdd=l1;
      }
      return 1;
    }
  }
  if (var2[0]>='0'&&var2[0]<='9') {
    if(varindex==2) {
      ops[nops].Var2Type=OT_CONST;
      ops[nops].Var2Val=atof(var2);
    } else {
      ops[nops].Var1Type=OT_CONST;
      ops[nops].Var1Val=atof(var2);
    }
    return 1;
  }
  return 0;
}

int leadlag_encode(char *line) {
  char *p1,*p2,*p3;
  p1=strchr(line,'{');
  while (p1!=NULL) {
    p3=strchr(p1,'}');
    p2=strchr(p1,'+');
    if(p2!=NULL) {
      while(p3-p2>0) {
        *p2='#';
        p2=strchr(p1,'+');
        if(p2==NULL)p2=p3;
      }
    }
    p2=strchr(p1,'-');
    if(p2!=NULL) {
      while(p3-p2>0) {
        *p2='!';
        p2=strchr(p1,'-');
        if(p2==NULL)p2=p3;
      }
    }
    p1=strchr(p3,'{');
  }
  return 1;
}

int parse_index_leadlag(char *p,int *leadlag) {
  char *plussign,*minsign;
  plussign=strchr(p,'#');
  minsign=strchr(p,'!');
  if(plussign!=NULL) {
    *leadlag=atoi(plussign+1);
    *plussign='\0';
  }
  if(minsign!=NULL) {
    *leadlag=-atoi(minsign+1);
    *minsign='\0';
  }
  return 1;
}

int formula_compile(char *fomulain, set_def *sets,array_def *coefs, offset_t ncof, array_def *vars,offset_t nvar,offset_t ncofele,sum_def *sum_cof,dim_t totalsum,formula_op *ops,dim_t *nops,quantifier *arSet,dim_t fdim) {
  int npar=0,npow=0,nmul=0,ndiv=0,nplu=0,nmin=0,j;
  *nops=0;
  npar=str_count_char(fomulain, ')');
  if (npar==0) {
    npow=str_count_char(fomulain, '^');
    if (npow>0) {
      formula_compile_pow(fomulain,sets,npow,0,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,nops,arSet,fdim);
    }
    nmul=str_count_char(fomulain, '*');
    ndiv=str_count_char(fomulain, '/');
    nmul=nmul+ndiv;
    if (nmul>0) {
      formula_compile_muldiv(fomulain,sets,nmul,0,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,nops,arSet,fdim);
    }
    nplu=str_count_char(fomulain, '+');
    nmin=str_count_char(fomulain, '-');
    nplu=nplu+nmin;
    if (nplu>0) {
      formula_compile_addsub(fomulain,sets,nplu,0,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,nops,arSet,fdim);
    }
      formula_bind_operand(fomulain,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,*nops,arSet,fdim,1);
      ops[*nops].Oper=OP_LOAD;
      ops[*nops].TmpVarName[0]='\0';
      *nops=*nops+1;
    return 1;
  }
  int i;
  char *p=NULL;
  char fpart1[TABREADLINE],fpart2[TABREADLINE],fpart3[TABREADLINE];
  char interchar[TABREADLINE],interchar1[TABREADLINE];

  for (i=1; i<npar+2; i++) {
    p=strchr(fomulain,')');
    if (p!=NULL) {
      strncpy(fpart1, fomulain, p-fomulain);
      fpart1[p-fomulain] = '\0';
      strcpy(fpart3,p+1);
      p=strrchr(fpart1,'(');
      if (p!=NULL) {
        strcpy(fpart2, p+1);
        strncpy(fpart1,fomulain,p-fpart1);
        fpart1[p-fpart1] = '\0';
      } else {
        printf("Error in Formula!\n");
        return 0;
      }
    } else {
      fpart1[0]='\0';
      strcpy(fpart2,fomulain);
      fpart3[0]='\0';
    }
    npow=str_count_char(fpart2, '^');
    if (npow>0) {
      formula_compile_pow(fpart2,sets,npow,i,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,nops,arSet,fdim);
    }
    nmul=str_count_char(fpart2, '*');
    ndiv=str_count_char(fpart2, '/');
    nmul=nmul+ndiv;
    if (nmul>0) {
      formula_compile_muldiv(fpart2,sets,nmul,i,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,nops,arSet,fdim);
    }
    nplu=str_count_char(fpart2, '+');
    nmin=str_count_char(fpart2, '-');
    nplu=nplu+nmin;
    if (nplu>0) {
      formula_compile_addsub(fpart2,sets,nplu,i,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,nops,arSet,fdim);
    }
      if(strpbrk(fpart2,"=<>")==NULL){
      formula_bind_operand(fpart2,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,*nops,arSet,fdim,1);
      ops[*nops].Oper=OP_LOAD;
      ops[*nops].TmpVarName[0]='\0';
      *nops=*nops+1;
      }
    sprintf(interchar1, "%d", i);
    interchar[0]='\0';
    if (i<10) {
      strcat(interchar,"gen_par000");
    }
    if (9<i&&i<100) {
      strcat(interchar,"gen_par00");
    }
    if (99<i&&i<1000) {
      strcat(interchar,"gen_par0");
    }
    if (999<i&&i<10000) {
      strcat(interchar,"gen_par");
    }
    if (i>10000) {
      strcat(interchar,"gen_par");
      printf("Warning: Too many parenthesises\n");
    }
    strcat(interchar,interchar1);
    strcpy(ops[*nops].TmpVarName,interchar);
    ops[*nops].Oper=OP_LOAD;
    ops[*nops].Var1BegAdd=*nops-1;
    j=strlen(fpart1);
    ops[*nops].Var1Type=OT_TEMP;
    if (j>3) if (fpart1[j-1]=='1'&&fpart1[j-2]=='0'&&tolower((int)fpart1[j-3])=='d'&&tolower((int)fpart1[j-4])=='i') {
        ops[*nops].Var1Type=OT_TEMP_ID01;
        fpart1[j-4]='\0';
      }
    if (j>3) if (fpart1[j-1]=='e'&&fpart1[j-2]=='g'&&fpart1[j-3]=='o'&&fpart1[j-4]=='l') {
        ops[*nops].Var1Type=OT_TEMP_LOG;
        fpart1[j-4]='\0';
      }
    if (j==3) if (fpart1[j-1]=='s'&&fpart1[j-2]=='b'&&fpart1[j-3]=='a') {
        ops[*nops].Var1Type=OT_TEMP_ABS;
        fpart1[j-3]='\0';
      }
    if (j>3) if (fpart1[j-1]=='s'&&fpart1[j-2]=='b'&&fpart1[j-3]=='a') if(j==3||fpart1[j-4]==' '||fpart1[j-4]=='('||fpart1[j-4]=='+'||fpart1[j-4]=='-'||fpart1[j-4]=='*'||fpart1[j-4]=='/'||fpart1[j-4]=='^'||fpart1[j-4]==',') {
          ops[*nops].Var1Type=OT_TEMP_ABS;
          fpart1[j-3]='\0';
        }
    if (j==2) if (fpart1[j-1]=='f'&&fpart1[j-2]=='i') {
        formula_compile_if(fpart2,sets,2,i,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,nops,arSet,fdim);
        fpart1[j-2]='\0';
      }
    if (j>2) if (fpart1[j-1]=='f'&&fpart1[j-2]=='i') if(j==2||fpart1[j-3]==' '||fpart1[j-3]=='('||fpart1[j-3]=='+'||fpart1[j-3]=='-'||fpart1[j-3]=='*'||fpart1[j-3]=='/'||fpart1[j-3]=='^'||fpart1[j-3]==',') {
          formula_compile_if(fpart2,sets,2,i,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,nops,arSet,fdim);
          fpart1[j-2]='\0';
        }
    *nops=*nops+1;
    strcat(fpart1, interchar);
    strcat(fpart1, fpart3);
    strcpy(fomulain,fpart1);
  }
  return 1;
}

solve_real formula_eval(elem_value *record,set_def *sets,set_element *set_elems,sum_value *sum_vals,formula_op *ops,int nops,quantifier *arSet,dim_t fdim, solve_real zerodivide) {
  int i;
  dim_t j;
  offset_t l=0,l1=0;
  solve_real eval1=0,eval2=0,eval3=0;
  for (i=0; i<nops; i++) {
    switch(ops[i].Oper) {
    case OP_LOAD:
      if (ops[i].Var1Type==OT_ARRAY) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ops[i].Var1Dims[j].SupSet==1) {
            l+=ops[i].Var1Dims[j].ADims*(set_elems[sets[arSet[j].setid].offset+arSet[j].indx].superset_pos[ops[i].Var1Dims[j].SSIndx]+ops[i].Var1Dims[j].leadlag);
          } else {
            l+=ops[i].Var1Dims[j].ADims*(arSet[j].indx+ops[i].Var1Dims[j].leadlag);
          }
        }
        ops[i].TmpVarVal=record[ops[i].Var1BegAdd+l].value;
        break;
      }
      if (ops[i].Var1Type==OT_SUM) {
        l=0;
        for (j=0; j<fdim; j++) {
          l+=ops[i].Var1Dims[j].ADims*arSet[j].indx;
        }
        ops[i].TmpVarVal=sum_vals[ops[i].Var1BegAdd+l].value;
        break;
      }
      if (ops[i].Var1Type==OT_TEMP) {
        ops[i].TmpVarVal=ops[ops[i].Var1BegAdd].TmpVarVal;
        break;
      }
      if (ops[i].Var1Type==OT_TEMP_ID01) {
        if(ops[ops[i].Var1BegAdd].TmpVarVal==0)ops[i].TmpVarVal=1;
        else ops[i].TmpVarVal=ops[ops[i].Var1BegAdd].TmpVarVal;
        break;
      }
      if (ops[i].Var1Type==OT_TEMP_ABS) {
        if(ops[ops[i].Var1BegAdd].TmpVarVal<0)ops[i].TmpVarVal=-ops[ops[i].Var1BegAdd].TmpVarVal;
        else ops[i].TmpVarVal=ops[ops[i].Var1BegAdd].TmpVarVal;
        break;
      }
      if (ops[i].Var1Type==OT_TEMP_LOG) {
        ops[i].TmpVarVal=log(ops[ops[i].Var1BegAdd].TmpVarVal);
        break;
      }
      if (ops[i].Var1Type==OT_CONST) {
        ops[i].TmpVarVal=ops[i].Var1Val;
        break;
      }
      if (ops[i].Var1Type==OT_CHANGE) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ops[i].Var1Dims[j].SupSet==1) {
            l+=ops[i].Var1Dims[j].ADims*(set_elems[sets[arSet[j].setid].offset+arSet[j].indx].superset_pos[ops[i].Var1Dims[j].SSIndx]+ops[i].Var1Dims[j].leadlag);
          } else {
            l+=ops[i].Var1Dims[j].ADims*(arSet[j].indx+ops[i].Var1Dims[j].leadlag);
          }
        }
        ops[i].TmpVarVal=record[ops[i].Var1BegAdd+l].substep_base;
        break;
      }
      break;
    case OP_MUL:
      if(ops[i].Var1Type<3||ops[i].Var1Type==OT_CHANGE) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ops[i].Var1Dims[j].SupSet==1) {
            l+=ops[i].Var1Dims[j].ADims*(set_elems[sets[arSet[j].setid].offset+arSet[j].indx].superset_pos[ops[i].Var1Dims[j].SSIndx]+ops[i].Var1Dims[j].leadlag);
          } else {
            l+=ops[i].Var1Dims[j].ADims*(arSet[j].indx+ops[i].Var1Dims[j].leadlag);
          }
        }
      }
      if(ops[i].Var2Type<3||ops[i].Var2Type==OT_CHANGE) {
        l1=0;
        for (j=0; j<fdim; j++) {
          if(ops[i].Var2Dims[j].SupSet==1) {
            l1+=ops[i].Var2Dims[j].ADims*(set_elems[sets[arSet[j].setid].offset+arSet[j].indx].superset_pos[ops[i].Var2Dims[j].SSIndx]+ops[i].Var2Dims[j].leadlag);
          } else {
            l1+=ops[i].Var2Dims[j].ADims*(arSet[j].indx+ops[i].Var2Dims[j].leadlag);
          }
        }
      }
      if(ops[i].Var1Type==OT_ARRAY) eval1=record[ops[i].Var1BegAdd+l].value;
      if(ops[i].Var1Type==OT_SUM) eval1=sum_vals[ops[i].Var1BegAdd+l].value;
      if(ops[i].Var1Type==OT_TEMP) eval1=ops[ops[i].Var1BegAdd].TmpVarVal;
      if(ops[i].Var1Type==OT_CONST) eval1=ops[i].Var1Val;
      if(ops[i].Var1Type==OT_CHANGE) eval1=record[ops[i].Var1BegAdd+l].substep_base;
      if(ops[i].Var2Type==OT_ARRAY) eval2=record[ops[i].Var2BegAdd+l1].value;
      if(ops[i].Var2Type==OT_SUM) eval2=sum_vals[ops[i].Var2BegAdd+l1].value;
      if(ops[i].Var2Type==OT_TEMP) eval2=ops[ops[i].Var2BegAdd].TmpVarVal;
      if(ops[i].Var2Type==OT_CONST) eval2=ops[i].Var2Val;
      if(ops[i].Var2Type==OT_CHANGE) eval2=record[ops[i].Var2BegAdd+l1].substep_base;
      ops[i].TmpVarVal=eval1*eval2;
      break;
    case OP_DIV:
      if(ops[i].Var1Type<3||ops[i].Var1Type==OT_CHANGE) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ops[i].Var1Dims[j].SupSet==1) {
            l+=ops[i].Var1Dims[j].ADims*(set_elems[sets[arSet[j].setid].offset+arSet[j].indx].superset_pos[ops[i].Var1Dims[j].SSIndx]+ops[i].Var1Dims[j].leadlag);
          } else {
            l+=ops[i].Var1Dims[j].ADims*(arSet[j].indx+ops[i].Var1Dims[j].leadlag);
          }
        }
      }
      if(ops[i].Var2Type<3||ops[i].Var2Type==OT_CHANGE) {
        l1=0;
        for (j=0; j<fdim; j++) {
          if(ops[i].Var2Dims[j].SupSet==1) {
            l1+=ops[i].Var2Dims[j].ADims*(set_elems[sets[arSet[j].setid].offset+arSet[j].indx].superset_pos[ops[i].Var2Dims[j].SSIndx]+ops[i].Var2Dims[j].leadlag);
          } else {
            l1+=ops[i].Var2Dims[j].ADims*(arSet[j].indx+ops[i].Var2Dims[j].leadlag);
          }
        }
      }
      if(ops[i].Var1Type==OT_ARRAY) eval1=record[ops[i].Var1BegAdd+l].value;
      if(ops[i].Var1Type==OT_SUM) eval1=sum_vals[ops[i].Var1BegAdd+l].value;
      if(ops[i].Var1Type==OT_TEMP) eval1=ops[ops[i].Var1BegAdd].TmpVarVal;
      if(ops[i].Var1Type==OT_CONST) eval1=ops[i].Var1Val;
      if(ops[i].Var1Type==OT_CHANGE) eval1=record[ops[i].Var1BegAdd+l].substep_base;
      if(ops[i].Var2Type==OT_ARRAY) eval2=record[ops[i].Var2BegAdd+l1].value;
      if(ops[i].Var2Type==OT_SUM) eval2=sum_vals[ops[i].Var2BegAdd+l1].value;
      if(ops[i].Var2Type==OT_TEMP) eval2=ops[ops[i].Var2BegAdd].TmpVarVal;
      if(ops[i].Var2Type==OT_CONST) eval2=ops[i].Var2Val;
      if(ops[i].Var2Type==OT_CHANGE) eval2=record[ops[i].Var2BegAdd+l1].substep_base;
      if(eval2==0) {
        ops[i].TmpVarVal=zerodivide;
      } else {
        ops[i].TmpVarVal=eval1/eval2;
      }
      break;
    case OP_ADD:
      if(ops[i].Var1Type==OT_ARRAY) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ops[i].Var1Dims[j].SupSet==1) {
            l+=ops[i].Var1Dims[j].ADims*(set_elems[sets[arSet[j].setid].offset+arSet[j].indx].superset_pos[ops[i].Var1Dims[j].SSIndx]+ops[i].Var1Dims[j].leadlag);
          } else {
            l+=ops[i].Var1Dims[j].ADims*(arSet[j].indx+ops[i].Var1Dims[j].leadlag);
          }
        }
        eval1=record[ops[i].Var1BegAdd+l].value;
      }
      if(ops[i].Var1Type==OT_SUM) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ops[i].Var1Dims[j].SupSet==1) {
            l+=ops[i].Var1Dims[j].ADims*set_elems[sets[arSet[j].setid].offset+arSet[j].indx].superset_pos[ops[i].Var1Dims[j].SSIndx];
          } else {
            l+=ops[i].Var1Dims[j].ADims*arSet[j].indx;
          }
        }
        eval1=sum_vals[ops[i].Var1BegAdd+l].value;
      }
      if (ops[i].Var1Type==OT_TEMP) {
        eval1=ops[ops[i].Var1BegAdd].TmpVarVal;
      }
      if (ops[i].Var1Type==OT_CONST) {
        eval1=ops[i].Var1Val;
      }
      if(ops[i].Var1Type==OT_CHANGE) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ops[i].Var1Dims[j].SupSet==1) {
            l+=ops[i].Var1Dims[j].ADims*(set_elems[sets[arSet[j].setid].offset+arSet[j].indx].superset_pos[ops[i].Var1Dims[j].SSIndx]+ops[i].Var1Dims[j].leadlag);
          } else {
            l+=ops[i].Var1Dims[j].ADims*(arSet[j].indx+ops[i].Var1Dims[j].leadlag);
          }
        }
        eval1=record[ops[i].Var1BegAdd+l].substep_base;
      }
      if(ops[i].Var2Type==OT_ARRAY) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ops[i].Var2Dims[j].SupSet==1) {
            l+=ops[i].Var2Dims[j].ADims*(set_elems[sets[arSet[j].setid].offset+arSet[j].indx].superset_pos[ops[i].Var2Dims[j].SSIndx]+ops[i].Var2Dims[j].leadlag);
          } else {
            l+=ops[i].Var2Dims[j].ADims*(arSet[j].indx+ops[i].Var2Dims[j].leadlag);
          }
        }
        eval2=record[ops[i].Var2BegAdd+l].value;
      }
      if(ops[i].Var2Type==OT_SUM) {
        l=0;
        for (j=0; j<fdim; j++) {
            l1=l;
          if(ops[i].Var2Dims[j].SupSet==1) {
            l+=ops[i].Var2Dims[j].ADims*set_elems[sets[arSet[j].setid].offset+arSet[j].indx].superset_pos[ops[i].Var2Dims[j].SSIndx];
          } else {
            l+=ops[i].Var2Dims[j].ADims*arSet[j].indx;
          }
        }
        eval2=sum_vals[ops[i].Var2BegAdd+l].value;
      }
      if (ops[i].Var2Type==OT_TEMP) {
        eval2=ops[ops[i].Var2BegAdd].TmpVarVal;
      }
      if (ops[i].Var2Type==OT_CONST) {
        eval2=ops[i].Var2Val;
      }
      if(ops[i].Var2Type==OT_CHANGE) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ops[i].Var2Dims[j].SupSet==1) {
            l+=ops[i].Var2Dims[j].ADims*(set_elems[sets[arSet[j].setid].offset+arSet[j].indx].superset_pos[ops[i].Var2Dims[j].SSIndx]+ops[i].Var2Dims[j].leadlag);
          } else {
            l+=ops[i].Var2Dims[j].ADims*(arSet[j].indx+ops[i].Var2Dims[j].leadlag);
          }
        }
        eval2=record[ops[i].Var2BegAdd+l].substep_base;
      }
      ops[i].TmpVarVal=eval1+eval2;
      break;
    case OP_SUB:
      if(ops[i].Var1Type==OT_ARRAY) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ops[i].Var1Dims[j].SupSet==1) {
            l+=ops[i].Var1Dims[j].ADims*(set_elems[sets[arSet[j].setid].offset+arSet[j].indx].superset_pos[ops[i].Var1Dims[j].SSIndx]+ops[i].Var1Dims[j].leadlag);
          } else {
            l+=ops[i].Var1Dims[j].ADims*(arSet[j].indx+ops[i].Var1Dims[j].leadlag);
          }
        }
        eval1=record[ops[i].Var1BegAdd+l].value;
      }
      if(ops[i].Var1Type==OT_SUM) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ops[i].Var1Dims[j].SupSet==1) {
            l+=ops[i].Var1Dims[j].ADims*set_elems[sets[arSet[j].setid].offset+arSet[j].indx].superset_pos[ops[i].Var1Dims[j].SSIndx];
          } else {
            l+=ops[i].Var1Dims[j].ADims*arSet[j].indx;
          }
        }
        eval1=sum_vals[ops[i].Var1BegAdd+l].value;
      }
      if (ops[i].Var1Type==OT_TEMP) {
        eval1=ops[ops[i].Var1BegAdd].TmpVarVal;
      }
      if (ops[i].Var1Type==OT_CONST) {
        eval1=ops[i].Var1Val;
      }
      if(ops[i].Var1Type==OT_CHANGE) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ops[i].Var1Dims[j].SupSet==1) {
            l+=ops[i].Var1Dims[j].ADims*(set_elems[sets[arSet[j].setid].offset+arSet[j].indx].superset_pos[ops[i].Var1Dims[j].SSIndx]+ops[i].Var1Dims[j].leadlag);
          } else {
            l+=ops[i].Var1Dims[j].ADims*(arSet[j].indx+ops[i].Var1Dims[j].leadlag);
          }
        }
        eval1=record[ops[i].Var1BegAdd+l].substep_base;
      }
      if(ops[i].Var2Type==OT_ARRAY) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ops[i].Var2Dims[j].SupSet==1) {
            l+=ops[i].Var2Dims[j].ADims*(set_elems[sets[arSet[j].setid].offset+arSet[j].indx].superset_pos[ops[i].Var2Dims[j].SSIndx]+ops[i].Var2Dims[j].leadlag);
          } else {
            l+=ops[i].Var2Dims[j].ADims*(arSet[j].indx+ops[i].Var2Dims[j].leadlag);
          }
        }
        eval2=record[ops[i].Var2BegAdd+l].value;
      }
      if(ops[i].Var2Type==OT_SUM) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ops[i].Var2Dims[j].SupSet==1) {
            l+=ops[i].Var2Dims[j].ADims*set_elems[sets[arSet[j].setid].offset+arSet[j].indx].superset_pos[ops[i].Var2Dims[j].SSIndx];
          } else {
            l+=ops[i].Var2Dims[j].ADims*arSet[j].indx;
          }
        }
        eval2=sum_vals[ops[i].Var2BegAdd+l].value;
      }
      if (ops[i].Var2Type==OT_TEMP) {
        eval2=ops[ops[i].Var2BegAdd].TmpVarVal;
      }
      if (ops[i].Var2Type==OT_CONST) {
        eval2=ops[i].Var2Val;
      }
      if(ops[i].Var2Type==OT_CHANGE) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ops[i].Var2Dims[j].SupSet==1) {
            l+=ops[i].Var2Dims[j].ADims*(set_elems[sets[arSet[j].setid].offset+arSet[j].indx].superset_pos[ops[i].Var2Dims[j].SSIndx]+ops[i].Var2Dims[j].leadlag);
          } else {
            l+=ops[i].Var2Dims[j].ADims*(arSet[j].indx+ops[i].Var2Dims[j].leadlag);
          }
        }
        eval2=record[ops[i].Var2BegAdd+l].substep_base;
      }
      ops[i].TmpVarVal=eval1-eval2;
      break;
    case OP_POW:
      if(ops[i].Var1Type<3||ops[i].Var1Type==OT_CHANGE) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ops[i].Var1Dims[j].SupSet==1) {
            l+=ops[i].Var1Dims[j].ADims*(set_elems[sets[arSet[j].setid].offset+arSet[j].indx].superset_pos[ops[i].Var1Dims[j].SSIndx]+ops[i].Var1Dims[j].leadlag);
          } else {
            l+=ops[i].Var1Dims[j].ADims*(arSet[j].indx+ops[i].Var1Dims[j].leadlag);
          }
        }
      }
      if(ops[i].Var2Type<3||ops[i].Var2Type==OT_CHANGE) {
        l1=0;
        for (j=0; j<fdim; j++) {
          if(ops[i].Var2Dims[j].SupSet==1) {
            l1+=ops[i].Var2Dims[j].ADims*(set_elems[sets[arSet[j].setid].offset+arSet[j].indx].superset_pos[ops[i].Var2Dims[j].SSIndx]+ops[i].Var2Dims[j].leadlag);
          } else {
            l1+=ops[i].Var2Dims[j].ADims*(arSet[j].indx+ops[i].Var2Dims[j].leadlag);
          }
        }
      }
      if(ops[i].Var1Type==OT_ARRAY) eval1=record[ops[i].Var1BegAdd+l].value;
      if(ops[i].Var1Type==OT_SUM) eval1=sum_vals[ops[i].Var1BegAdd+l].value;
      if(ops[i].Var1Type==OT_TEMP) eval1=ops[ops[i].Var1BegAdd].TmpVarVal;
      if(ops[i].Var1Type==OT_CONST) eval1=ops[i].Var1Val;
      if(ops[i].Var1Type==OT_CHANGE) eval1=record[ops[i].Var1BegAdd+l].substep_base;
      if(ops[i].Var2Type==OT_ARRAY) eval2=record[ops[i].Var2BegAdd+l1].value;
      if(ops[i].Var2Type==OT_SUM) eval2=sum_vals[ops[i].Var2BegAdd+l1].value;
      if(ops[i].Var2Type==OT_TEMP) eval2=ops[ops[i].Var2BegAdd].TmpVarVal;
      if(ops[i].Var2Type==OT_CONST) eval2=ops[i].Var2Val;
      if(ops[i].Var2Type==OT_CHANGE) eval2=record[ops[i].Var2BegAdd+l1].substep_base;
      if(eval1==0&&eval2<0) {
        ops[i].TmpVarVal=zerodivide;
      } else {
        if(eval1<0&&eval2-floor(eval2)!=0)printf("Serious errors: fraction power of negative number!!!!!!!");
        ops[i].TmpVarVal=pow(eval1,eval2);
      }
      break;
    default:
      if(ops[i].Var1Type==OT_ARRAY) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ops[i].Var1Dims[j].SupSet==1) {
            l+=ops[i].Var1Dims[j].ADims*(set_elems[sets[arSet[j].setid].offset+arSet[j].indx].superset_pos[ops[i].Var1Dims[j].SSIndx]+ops[i].Var1Dims[j].leadlag);
          } else {
            l+=ops[i].Var1Dims[j].ADims*(arSet[j].indx+ops[i].Var1Dims[j].leadlag);
          }
        }
        eval1=record[ops[i].Var1BegAdd+l].value;
      }
      if(ops[i].Var1Type==OT_SUM) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ops[i].Var1Dims[j].SupSet==1) {
            l+=ops[i].Var1Dims[j].ADims*set_elems[sets[arSet[j].setid].offset+arSet[j].indx].superset_pos[ops[i].Var1Dims[j].SSIndx];
          } else {
            l+=ops[i].Var1Dims[j].ADims*arSet[j].indx;
          }
        }
        eval1=sum_vals[ops[i].Var1BegAdd+l].value;
      }
      if (ops[i].Var1Type==OT_TEMP) {
        eval1=ops[ops[i].Var1BegAdd].TmpVarVal;
      }
      if (ops[i].Var1Type==OT_CONST) {
        eval1=ops[i].Var1Val;
      }
      if(ops[i].Var1Type==OT_CHANGE) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ops[i].Var1Dims[j].SupSet==1) {
            l+=ops[i].Var1Dims[j].ADims*(set_elems[sets[arSet[j].setid].offset+arSet[j].indx].superset_pos[ops[i].Var1Dims[j].SSIndx]+ops[i].Var1Dims[j].leadlag);
          } else {
            l+=ops[i].Var1Dims[j].ADims*(arSet[j].indx+ops[i].Var1Dims[j].leadlag);
          }
        }
        eval1=record[ops[i].Var1BegAdd+l].substep_base;
      }
      if(ops[i].Var2Type==OT_ARRAY) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ops[i].Var2Dims[j].SupSet==1) {
            l+=ops[i].Var2Dims[j].ADims*(set_elems[sets[arSet[j].setid].offset+arSet[j].indx].superset_pos[ops[i].Var2Dims[j].SSIndx]+ops[i].Var2Dims[j].leadlag);
          } else {
            l+=ops[i].Var2Dims[j].ADims*(arSet[j].indx+ops[i].Var2Dims[j].leadlag);
          }
        }
        eval2=record[ops[i].Var2BegAdd+l].value;
      }
      if(ops[i].Var2Type==OT_SUM) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ops[i].Var2Dims[j].SupSet==1) {
            l+=ops[i].Var2Dims[j].ADims*set_elems[sets[arSet[j].setid].offset+arSet[j].indx].superset_pos[ops[i].Var2Dims[j].SSIndx];
          } else {
            l+=ops[i].Var2Dims[j].ADims*arSet[j].indx;
          }
        }
        eval2=sum_vals[ops[i].Var2BegAdd+l].value;
      }
      if (ops[i].Var2Type==OT_TEMP) {
        eval2=ops[ops[i].Var2BegAdd].TmpVarVal;
      }
      if (ops[i].Var2Type==OT_CONST) {
        eval2=ops[i].Var2Val;
      }
      if(ops[i].Var2Type==OT_CHANGE) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ops[i].Var2Dims[j].SupSet==1) {
            l+=ops[i].Var2Dims[j].ADims*(set_elems[sets[arSet[j].setid].offset+arSet[j].indx].superset_pos[ops[i].Var2Dims[j].SSIndx]+ops[i].Var2Dims[j].leadlag);
          } else {
            l+=ops[i].Var2Dims[j].ADims*(arSet[j].indx+ops[i].Var2Dims[j].leadlag);
          }
        }
        eval2=record[ops[i].Var2BegAdd+l].substep_base;
      }

      if(ops[i].Var3Type==OT_ARRAY) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ops[i].Var3Dims[j].SupSet==1) {
            l+=ops[i].Var3Dims[j].ADims*(set_elems[sets[arSet[j].setid].offset+arSet[j].indx].superset_pos[ops[i].Var3Dims[j].SSIndx]+ops[i].Var3Dims[j].leadlag);
          } else {
            l+=ops[i].Var3Dims[j].ADims*(arSet[j].indx+ops[i].Var3Dims[j].leadlag);
          }
        }
        eval3=record[ops[i].Var3BegAdd+l].value;
      }
      if(ops[i].Var3Type==OT_SUM) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ops[i].Var3Dims[j].SupSet==1) {
            l+=ops[i].Var3Dims[j].ADims*set_elems[sets[arSet[j].setid].offset+arSet[j].indx].superset_pos[ops[i].Var3Dims[j].SSIndx];
          } else {
            l+=ops[i].Var3Dims[j].ADims*arSet[j].indx;
          }
        }
        eval3=sum_vals[ops[i].Var3BegAdd+l].value;
      }
      if (ops[i].Var3Type==OT_TEMP) {
        eval3=ops[ops[i].Var3BegAdd].TmpVarVal;
      }
      if (ops[i].Var3Type==OT_CONST) {
        eval3=ops[i].Var3Val;
      }
      if(ops[i].Var3Type==OT_CHANGE) {
        l=0;
        for (j=0; j<fdim; j++) {
          if(ops[i].Var3Dims[j].SupSet==1) {
            l+=ops[i].Var3Dims[j].ADims*(set_elems[sets[arSet[j].setid].offset+arSet[j].indx].superset_pos[ops[i].Var3Dims[j].SSIndx]+ops[i].Var3Dims[j].leadlag);
          } else {
            l+=ops[i].Var3Dims[j].ADims*(arSet[j].indx+ops[i].Var3Dims[j].leadlag);
          }
        }
        eval3=record[ops[i].Var3BegAdd+l].substep_base;
      }
      if(ops[i].Oper==OP_IF_EQ)if(eval1==eval2)ops[i].TmpVarVal=eval3;else ops[i].TmpVarVal=0;
      if(ops[i].Oper==OP_IF_GT)if(eval1>eval2)ops[i].TmpVarVal=eval3;else ops[i].TmpVarVal=0;
      if(ops[i].Oper==OP_IF_LT)if(eval1<eval2)ops[i].TmpVarVal=eval3;else ops[i].TmpVarVal=0;
      if(ops[i].Oper==OP_IF_NE)if(eval1!=eval2)ops[i].TmpVarVal=eval3;else ops[i].TmpVarVal=0;
      if(ops[i].Oper==OP_IF_LE)if(eval1<=eval2)ops[i].TmpVarVal=eval3;else ops[i].TmpVarVal=0;
      if(ops[i].Oper==OP_IF_GE)if(eval1>=eval2)ops[i].TmpVarVal=eval3;else ops[i].TmpVarVal=0;
      break;
    }
  }
  return ops[i-1].TmpVarVal;
}




int formula_compile_pow(char *fomulain, set_def *sets,int npow,int ipar,array_def *coefs,offset_t ncof, array_def *vars,offset_t nvar,offset_t ncofele,sum_def *sum_cof,int totalsum,formula_op *ops,int *nops,quantifier *arSet,dim_t fdim) {
  int i,i1,i2,i3,i4,i0,ibar=0,index,j,j1,i5,p1;
  char *p=NULL;//,*p1=NULL,*p2=NULL,*p3=NULL,*p4=NULL;
  char fpart1[TABREADLINE],fpart2[TABREADLINE],fpart3[TABREADLINE],var1[TABREADLINE],var2[TABREADLINE],interchar[TABREADLINE],interchar1[TABREADLINE];

  for (i=1; i<npow+1; i++) {
    p=strchr(fomulain,'^');
    index=p-fomulain;
    strncpy(fpart1, fomulain, index);
    fpart1[index] = '\0';
    strcpy(fpart3,fomulain+index+1);

    i1=-1;
    p=str_rfind_any(fpart1,"^*/+-=<>");
    if (p!=NULL) {
      i1=p-fpart1;
    }
    i5=-1;
    j=strlen(fpart1)-1;
    while(j>-1){
      if(fpart1[j]==','){
        p1=-1;
        for(j1=j-1;j1>-1;j1--){
          if(fpart1[j1]=='{')p1=j1;
          if(fpart1[j1]=='}')break;
        }
        if(p1==-1)break;
        else j=p1;
      }
      j--;
    }
    if(j>0)i5=j;
    if (i1==i5==-1) {//if (i0==i1==i2==i3==i4==i5==-1) {
      index=0;
    } else {
      index=i1;
      if (index<i5) {
        index=i5;
      }
    }
    strcpy(var1, fpart1+index+1);
    strncpy(fpart1,fomulain,index+1);
    fpart1[index+1] = '\0';

    ibar=0;
    while (fpart3[ibar] != '\0') {
      ibar++;
    }
    p=strpbrk(fpart3,"^*/+-=<>");
    if(p==NULL)index=ibar;
    else index=p-fpart3;
    strncpy(var2, fpart3, index);
    var2[index] = '\0';
    strcpy(fpart2,fpart3+index);

    formula_bind_operand(var1,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,*nops,arSet,fdim,1);
    formula_bind_operand(var2,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,*nops,arSet,fdim,2);

    sprintf(interchar1, "%d", ipar);
    interchar[0]='\0';
    if (i<10) {
      strcat(interchar,"gen_pow");
      strcat(interchar,interchar1);
      strcat(interchar,"000");
    }
    if (9<i&&i<100) {
      strcat(interchar,"gen_pow");
      strcat(interchar,interchar1);
      strcat(interchar,"00");
    }
    if (99<i&&i<1000) {
      strcat(interchar,"gen_pow");
      strcat(interchar,interchar1);
      strcat(interchar,"0");
    }
    if (999<i&&i<10000) {
      strcat(interchar,"gen_pow");
      strcat(interchar,interchar1);
    }
    sprintf(interchar1, "%d", i);
    strcat(interchar,interchar1);
    strcpy(ops[*nops].TmpVarName,interchar);
    ops[*nops].Oper=OP_POW;
    *nops=*nops+1;
    strcat(fpart1, interchar);
    strcat(fpart1, fpart2);
    strcpy(fomulain,fpart1);
  }
  return 1;
}
int formula_compile_muldiv(char *fomulain, set_def *sets,int nmul,int ipar,array_def *coefs,offset_t ncof, array_def *vars,offset_t nvar,offset_t ncofele,sum_def *sum_cof,int totalsum,formula_op *ops,int *nops,quantifier *arSet,dim_t fdim) {
  int i,i1,i2,i3,i4,i5,ibar=0,index,j,j1,p1;
  char *p=NULL;//,*p1=NULL,*p2=NULL,*p3=NULL,*p4=NULL;
  char fpart1[TABREADLINE],fpart2[TABREADLINE],fpart3[TABREADLINE],var1[TABREADLINE],var2[TABREADLINE],interchar[TABREADLINE],interchar1[TABREADLINE];

  for (i=1; i<nmul+1; i++) {
    index=0;
    p=strpbrk(fomulain,"*/");
    if (*p=='/') {
      ops[*nops].Oper=OP_DIV;
    } else {
      ops[*nops].Oper=OP_MUL;
    }
    index=p-fomulain;

    strncpy(fpart1, fomulain, index);
    fpart1[index] = '\0';
    strcpy(fpart3,fomulain+index+1);
    i1=-1;
    p=str_rfind_any(fpart1,"*/+-=<>");
    if (p!=NULL) {
      i1=p-fpart1;
    }
    i5=-1;
    j=strlen(fpart1)-1;
    while(j>-1){
      if(fpart1[j]==','){
        p1=-1;
        for(j1=j-1;j1>-1;j1--){
          if(fpart1[j1]=='{')p1=j1;
          if(fpart1[j1]=='}')break;
        }
        if(p1==-1)break;
        else j=p1;
      }
      j--;
    }
    if(j>0)i5=j;
    if (i1==i5==-1) {//if (i1==i2==i3==i4==i5==-1) {
      index=0;
    } else {
      index=i1;
      if (index<i5) {
        index=i5;
      }
    }
    strcpy(var1, fpart1+index+1);
    strncpy(fpart1,fomulain,index+1);
    fpart1[index+1] = '\0';

    ibar=0;
    while (fpart3[ibar] != '\0') {
      ibar++;
    }
    p=strpbrk(fpart3,"*/+-=<>");
    if(p==NULL)index=ibar;
    else index=p-fpart3;
    strncpy(var2, fpart3, index);
    var2[index] = '\0';
    strcpy(fpart2,fpart3+index);

    formula_bind_operand(var1,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,*nops,arSet,fdim,1);
    formula_bind_operand(var2,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,*nops,arSet,fdim,2);

    sprintf(interchar1, "%d", ipar);
    interchar[0]='\0';
    if (i<10) {
      strcat(interchar,"gen_mul");
      strcat(interchar,interchar1);
      strcat(interchar,"000");
    }
    if (9<i&&i<100) {
      strcat(interchar,"gen_mul");
      strcat(interchar,interchar1);
      strcat(interchar,"00");
    }
    if (99<i&&i<1000) {
      strcat(interchar,"gen_mul");
      strcat(interchar,interchar1);
      strcat(interchar,"0");
    }
    if (999<i&&i<10000) {
      strcat(interchar,"gen_mul");
      strcat(interchar,interchar1);
    }
    sprintf(interchar1, "%d", i);
    strcat(interchar,interchar1);
    strcpy(ops[*nops].TmpVarName,interchar);
    *nops=*nops+1;
    strcat(fpart1, interchar);
    strcat(fpart1, fpart2);
    strcpy(fomulain,fpart1);
  }
  return 1;
}

int formula_compile_addsub(char *fomulain, set_def *sets,int nplu,int ipar,array_def *coefs,offset_t ncof, array_def *vars,offset_t nvar,offset_t ncofele,sum_def *sum_cof,int totalsum,formula_op *ops,int *nops,quantifier *arSet,dim_t fdim) {
  int i,i1,i3,i4,i5,ibar=0,index,j,j1,p1;
  char *p=NULL;//,*p1=NULL,*p2=NULL,*p3=NULL,*p4=NULL;
  char fpart1[TABREADLINE],fpart2[TABREADLINE],fpart3[TABREADLINE],var1[TABREADLINE],var2[TABREADLINE],interchar[TABREADLINE],interchar1[TABREADLINE];

  for (i=1; i<nplu+1; i++) {
    index=0;
    p=strpbrk(fomulain,"+-");
    if (*p=='+') {
      ops[*nops].Oper=OP_ADD;
    } else {
      ops[*nops].Oper=OP_SUB;
    }
    index=p-fomulain;

    strncpy(fpart1, fomulain, index);
    fpart1[index] = '\0';
    strcpy(fpart3,fomulain+index+1);
    i1=-1;
    p=str_rfind_any(fpart1,"+-=<>");
    if (p!=NULL) {
      i1=p-fpart1;
    }
    i5=-1;
    j=strlen(fpart1)-1;
    while(j>-1){
      if(fpart1[j]==','){
        p1=-1;
        for(j1=j-1;j1>-1;j1--){
          if(fpart1[j1]=='{')p1=j1;
          if(fpart1[j1]=='}')break;
        }
        if(p1==-1)break;
        else j=p1;
      }
      j--;
    }
    if(j>0)i5=j;
    if (i1==i5==-1) {//if (i3==i4==i5==-1) {
      index=0;
    } else {
      index=i1;
      if (index<i5) {
        index=i5;
      }
    }
    strcpy(var1, fpart1+index+1);
    strncpy(fpart1,fomulain,index+1);
    fpart1[index+1] = '\0';

    ibar=0;
    while (fpart3[ibar] != '\0') {
      ibar++;
    }
    p=strpbrk(fpart3,"+-=<>");
    if(p==NULL)index=ibar;
    else index=p-fpart3;
    strncpy(var2, fpart3, index);
    var2[index] = '\0';
    strcpy(fpart2,fpart3+index);

    if(i==1&&var1[0]=='\0'){
        ops[*nops].Var1Type=OT_CONST;
        ops[*nops].Var1Val=0;
    }
    else formula_bind_operand(var1,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,*nops,arSet,fdim,1);
    if(i==1&&var2[0]=='\0'){
        ops[*nops].Var2Type=OT_CONST;
        ops[*nops].Var2Val=0;
    }
    else formula_bind_operand(var2,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,*nops,arSet,fdim,2);

    sprintf(interchar1, "%d", ipar);
    interchar[0]='\0';
    if (i<10) {
      strcat(interchar,"gen_add");
      strcat(interchar,interchar1);
      strcat(interchar,"000");
    }
    if (9<i&&i<100) {
      strcat(interchar,"gen_add");
      strcat(interchar,interchar1);
      strcat(interchar,"00");
    }
    if (99<i&&i<1000) {
      strcat(interchar,"gen_add");
      strcat(interchar,interchar1);
      strcat(interchar,"0");
    }
    if (999<i&&i<10000) {
      strcat(interchar,"gen_add");
      strcat(interchar,interchar1);
    }
    sprintf(interchar1, "%d", i);
    strcat(interchar,interchar1);
    strcpy(ops[*nops].TmpVarName,interchar);
    *nops=*nops+1;
    strcat(fpart1, interchar);
    strcat(fpart1, fpart2);
    strcpy(fomulain,fpart1);
  }
  return 1;
}


int formula_compile_if(char *fomulain, set_def *sets,int nif,int ipar,array_def *coefs,offset_t ncof, array_def *vars,offset_t nvar,offset_t ncofele,sum_def *sum_cof,int totalsum,formula_op *ops,int *nops,quantifier *arSet,dim_t fdim) {
  char *p=NULL,*p1,*p3,var1[NAMESIZE],var2[NAMESIZE],var3[NAMESIZE];
  int i,j1,j2,j3,l;//,varindex;
  p1=fomulain;
  p=strpbrk(p1,"=<>");
  strncpy(var1,p1,p-p1);
  var1[p-p1]='\0';
  p3=p+1;
  if(*p3=='='){
    if(*p=='<')ops[*nops].Oper=OP_IF_LE;
    if(*p=='>')ops[*nops].Oper=OP_IF_GE;
    p++;
    p++;
  }else if(*p3=='>'){
          ops[*nops].Oper=OP_IF_NE;
          p++;
          p++;
          }else{
          if(*p=='=')ops[*nops].Oper=OP_IF_EQ;
          if(*p=='>')ops[*nops].Oper=OP_IF_GT;
          if(*p=='<')ops[*nops].Oper=OP_IF_LT;
          p++;
        }
  l=strlen(p);
  j1=-1;j2=-1;j3=-1;
  for(i=0;i<l;i++){
    if(*(p+i)==',')j1=i;
    if(*(p+i)=='}')j2=i;
    if(*(p+i)=='{')j3=i;
    if(j3>-1)if(j1>-1&&j1>j2)break;
    else if(j1>-1)break;
  }
  strncpy(var2,p,j1);
  var2[j1]='\0';
  strcpy(var3,p+j1+1);
  formula_bind_operand(var1,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,*nops,arSet,fdim,1);
  formula_bind_operand(var3,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,*nops,arSet,fdim,2);
  ops[*nops].Var3Type=ops[*nops].Var2Type;
  ops[*nops].Var3BegAdd=ops[*nops].Var2BegAdd;
  for(i=0;i<fdim;i++){
  ops[*nops].Var3Dims[i].leadlag=ops[*nops].Var2Dims[i].leadlag;
  ops[*nops].Var3Dims[i].SupSet=ops[*nops].Var2Dims[i].SupSet;
  ops[*nops].Var3Dims[i].SSIndx=ops[*nops].Var2Dims[i].SSIndx;
  ops[*nops].Var3Dims[i].ADims=ops[*nops].Var2Dims[i].ADims;
  }
  ops[*nops].Var3Val=ops[*nops].Var2Val;
  formula_bind_operand(var2,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,*nops,arSet,fdim,2);
  return 1;
}

offset_t formulas_execute(char *fname, char *commsyntax,set_def *sets,dim_t nset, set_element *set_elems, array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar, elem_value *elem_vals,offset_t ncofvar,offset_t ncofele,bool IsIni) {
  FILE * filehandle;
  char line[TABREADLINE],line1[TABREADLINE],line2[TABREADLINE],linecopy[TABREADLINE],condvar[MAXVARDIM][NAMESIZE];
  char vname[NAMESIZE],sumsyntax[NAMESIZE],argu[NAMESIZE],tempset[NAMESIZE];
  char *readitem=NULL,*p=NULL,*p1=NULL;
  offset_t i,i1,i3,i4,l,l2=0,j=0,nsumele,dcountdim1[4*MAXVARDIM],ncond,nloops,logioper[MAXVARDIM],logi,logiantidim[MAXVARDIM][MAXVARDIM],logisup[MAXVARDIM][MAXVARDIM],logivarindx[MAXVARDIM],logivartype[MAXVARDIM],index;//m,
  dim_t fdim,dcount,neqsign=0,sup,varsupsetid[MAXVARDIM];
  int nops=0,totalsum,sumcount=1,npow,nmul,ndiv,nplu,nmin,npar,sumindx,b=0;
  offset_t varantidim[MAXVARDIM],varsubset[MAXVARDIM],vararset[MAXVARDIM];
  solve_real zerodivide=0,cond[MAXVARDIM],eval;
  bool IsFomIni=false,IsDefFomIni=false;
  quantifier *arSet1=NULL;
  formula_op *ops1= NULL;

  filehandle = fopen(fname,"r");
  while (tab_next_statement_resolved(commsyntax,filehandle,line,elem_vals,coefs,ncof,&zerodivide,TABREADLINE)) {
    if (strstr(line,"(default=initial)")!=NULL) IsDefFomIni=true;
    if (strstr(line,"(default")==NULL) {
      IsFomIni=IsDefFomIni;
      if(strstr(line, "(initial)")!=NULL) {
        str_replace_first(line, "(initial)", "");
        IsFomIni=true;
      }
      if(strstr(line, "(always)")!=NULL) {
        str_replace_first(line, "(always)", "");
        IsFomIni=false;
      }
      if(!(IsFomIni&&!IsIni)) {
        ncond=0;
        for (i=0; i<MAXVARDIM; i++)logioper[i]=0;
        str_replace_first(line, commsyntax, "");
        while (str_replace_all(line," ", ""));
        while (str_replace_char(line, '[', '('));
        while (str_replace_char(line, ']', ')'));
        while (str_replace_char(line, '{', '('));
        while (str_replace_char(line, '}', ')'));
        strcpy(linecopy,line);
        totalsum=str_count_ci(line, "sum(");
        neqsign=str_count_char(line, '=');
        readitem = strtok(line,"=");
        for(i=1; i<neqsign; i++)readitem = strtok(NULL,"=");
        readitem = strtok(NULL,";");
        npow=str_count_char(readitem, '^');
        nmul=str_count_char(readitem, '*');
        ndiv=str_count_char(readitem, '/');
        nmul=nmul+ndiv;
        nplu=str_count_char(readitem, '+');
        nmin=str_count_char(readitem, '-');
        nplu=nplu+nmin;

        strcpy(line,linecopy);
        readitem =str_rfind_toplevel(line,'=');//readitem =strrchr(line,'=');
        line[readitem-line]='\0';
        fdim=str_count_char(line, '(');
        if (fdim==1) {
          fdim=fdim+1;
        }
        quantifier *arSet= (quantifier *) calloc (fdim+1,sizeof(quantifier));

        nloops=1;
        if (fdim==0) {
          readitem = line;
        } else {
          for (i=0; i<fdim-1; i++) {
            if(i==0) {
              readitem = strtok(line,",");
            } else {
              readitem = strtok(NULL,",");
            }
            readitem = strtok(NULL,",");
            strcpy(arSet[i].index_name,readitem);
            readitem = strtok(NULL,")");
            if(strchr(readitem,':')==NULL)strcpy(tempset,readitem);
            else {
              p = strchr(readitem,':');
              strncpy(tempset,readitem,p-readitem);
              tempset[p-readitem]='\0';
              if(strchr(readitem,'(')!=NULL) {
                fdim--;
                p++;
                strcpy(condvar[i],p);
                strcat(condvar[i],")");
                readitem = strtok(NULL,")");
                if(strstr(readitem,"=")){
                  logioper[i]=1;
                  if(strstr(readitem,"<="))logioper[i]=5;
                  if(strstr(readitem,">="))logioper[i]=6;
                }else{
                  if(strstr(readitem,">"))logioper[i]=2;
                  if(strstr(readitem,"<"))logioper[i]=3;
                  if(strstr(readitem,"<>"))logioper[i]=4;
                }
                cond[i]=atof(readitem);
              } else {
                p++;
                p1=strstr(readitem,"=");
                if(p1!=NULL) {
                  logioper[i]=1;
                  strncpy(condvar[i],p,p1-p);
                  cond[i]=atof(readitem);
                p1=strstr(readitem,"<=");
                if(p1!=NULL) {
                  logioper[i]=5;
                  strncpy(condvar[i],p,p1-p);
                  cond[i]=atof(readitem);
                }
                p1=strstr(readitem,">=");
                if(p1!=NULL) {
                  logioper[i]=6;
                  strncpy(condvar[i],p,p1-p);
                  cond[i]=atof(readitem);
                }
                }else{
                p1=strstr(readitem,">");
                if(p1!=NULL) {
                  logioper[i]=2;
                  strncpy(condvar[i],p,p1-p);
                  cond[i]=atof(readitem);
                }
                p1=strstr(readitem,"<");
                if(p1!=NULL) {
                  logioper[i]=3;
                  strncpy(condvar[i],p,p1-p);
                  cond[i]=atof(readitem);
                }
                p1=strstr(readitem,"<>");
                if(p1!=NULL) {
                  logioper[i]=4;
                  strncpy(condvar[i],p,p1-p);
                  cond[i]=atof(readitem);
                }
                }
              }
              ncond++;
            }
            for (i4=0; i4<nset; i4++) if(strcmp(tempset,sets[i4].setname)==0) {
                arSet[i].setid=i4;
                break;
              }
            nloops=nloops*sets[arSet[i].setid].size;
          }
          readitem = strtok(NULL,"=");
          dcountdim1[fdim-2]=1;
          for (i=fdim-3; i>-1; i--) {
            dcountdim1[i]=sets[arSet[i+1].setid].size*dcountdim1[i+1];
          }
        }
        strcpy(vname,readitem);
        strcpy(line,linecopy);
        readitem =str_rfind_toplevel(line,'=');
        readitem++;
        readitem = strtok(readitem,";");
        while (formula_normalize(readitem)==1);
        leadlag_encode(readitem);
        npar=str_count_char(readitem, '(');
        strcpy(sumsyntax,"sum(");
        totalsum=sum_count(readitem,sumsyntax);
        sum_def *sum_cof= (sum_def *) calloc (totalsum,sizeof(sum_def));
        sumcount=0;
        strcpy(line1,readitem);
        strcpy(line2,line1);
        readitem=line2;
        while (sum_parse(readitem,sumsyntax,sum_cof,arSet,sets,nset,fdim,sumcount)==1) {
          sumcount++;
        }
        totalsum=sumcount;
        i3=0;
        for (i=0; i<totalsum; i++) {
          i1=1;
          for(j=0; j<sum_cof[i].size; j++) {
            i1=i1*sets[sum_cof[i].setid[j]].size;
          }
          sum_cof[i].offset=i3;
          sum_cof[i].summatsize=i1;
          i3=i3+i1;
        }
        nsumele=i3;
        for (i=0; i<totalsum; i++) {
          i1=1;
          sum_cof[i].strides[sum_cof[i].size-1]=1;
          for(j=sum_cof[i].size-2; j>-1; j--) {
            sum_cof[i].strides[j]=sum_cof[i].strides[j+1]*sets[sum_cof[i].setid[j+1]].size;
          }
        }
        formula_op *ops= (formula_op *) calloc (npow+nmul+nplu+2*(npar+2),sizeof(formula_op));
        sum_value *sum_vals= (sum_value *) calloc (nsumele,sizeof(sum_value));
        sumcount=0;
        strcpy(line2,line1);
        readitem=line2;
        sumindx=0;
        while (sum_eval(readitem,sumsyntax,sets,nset,set_elems,elem_vals,ncofvar,ncofele,coefs,ncof,vars,nvar,sum_cof,totalsum,sum_vals,nsumele,ops,arSet,fdim,&sumindx,sumcount,zerodivide)==1) {
          sumcount++;
        }
        strcpy(line1,readitem);
        offset_t index=ncof-1, offset=0;//,simpl=0;
        bool check10=true;
        offset_t varsize=0;
        p=strtok(vname,"(");
        do {
          if (strcmp(coefs[index].cofname,p)==0) {
            coefs[index].suplval=true;
            offset=coefs[index].offset;
            varsize=coefs[index].size;
            if(coefs[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
            check10=false;
            break;
          }
        } while (index--);
        if (check10) {
          index=nvar-1;
          do {
            if (strcmp(vars[index].cofname,p)==0) {
              offset=ncofele+vars[index].offset;
              varsize=vars[index].size;
              vars[index].suplval=true;
              if(vars[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
              break;
            }
          } while (index--);
        }
        for (l=0; l<MAXVARDIM; l++){varantidim[l]=0;varsubset[l]=0;varsupsetid[l]=0;}
        if (check10) {
            for (dcount=0; dcount<vars[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].index_name,p)==0) {
                  varantidim[dcount]=vars[index].strides[dcount];
                  vararset[dcount]=l+1;
                  if (sets[arSet[l].setid].size!=sets[vars[index].setid[dcount]].size){
                    varsubset[dcount]=1;
                    for(sup=1;sup<MAXSUPSET;sup++)if(sets[arSet[l].setid].subsetid[sup]==vars[index].setid[dcount]){varsupsetid[dcount]=sup;break;}
                  }
                  break;
                }
            }
        } else {
            for (dcount=0; dcount<coefs[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].index_name,p)==0) {
                  varantidim[dcount]=coefs[index].strides[dcount];
                  vararset[dcount]=l+1;
                  if (sets[arSet[l].setid].size!=sets[coefs[index].setid[dcount]].size){
                    varsubset[dcount]=1;
                    for(sup=1;sup<MAXSUPSET;sup++)if(sets[arSet[l].setid].subsetid[sup]==coefs[index].setid[dcount]){varsupsetid[dcount]=sup;break;}
                  }
                  break;
                }
            }
        }
        formula_compile(line1,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,&nops,arSet,fdim-1);
        if(ncond>0) {
          for(i=0; i<MAXVARDIM; i++)for(j=0; j<MAXVARDIM; j++){
            logiantidim[i][j]=0;
            logisup[i][j]=0;
          }
          for(i1=0; i1<fdim-1; i1++) if(logioper[i1]>0){
            index=ncof-1;
            p=strtok(condvar[i1],"(");
            b=0;
            do {
              if (strcmp(coefs[index].cofname,condvar[i1])==0) {
                if(!coefs[index].suplval)printf("Warning!!!! coefficient %s has not been supplied with values!\n",coefs[index].cofname);
                logivarindx[i1]=index;
                logivartype[i1]=0;
                b++;
                if(p!=NULL)p=strtok(NULL,")");
                strcpy(argu,p);
                strcat(argu,",");
                for(i=0; i<coefs[index].size; i++) {
                  if(i==0)p=strtok(argu,",");
                  else p=strtok(NULL,",");
                  for(j=0; j<fdim; j++) {
                    if(strcmp(arSet[j].index_name,p)==0) {
                      logiantidim[i1][j]=coefs[index].strides[i];
                      for(sup=1; sup<MAXSUPSET; sup++)if(sets[arSet[j].setid].subsetid[sup]==coefs[index].setid[i]){logisup[i1][j]=sup;break;}
                      break;
                    }
                  }
                }
                break;
              }
            } while (index--);
            if(b==0) {
              index=nvar-1;
              p=strtok(condvar[i1],"(");
              do {
                if (strcmp(vars[index].cofname,condvar[i1])==0) {
                if(!vars[index].suplval)printf("Warning!!!! coefficient %s has not been supplied with values!\n",vars[index].cofname);
                  logivarindx[i1]=index;
                  logivartype[i1]=1;
                  if(p!=NULL)p=strtok(NULL,")");
                  strcpy(argu,p);
                  strcat(argu,",");
                  for(i=0; i<vars[index].size; i++) {
                    if(i==0)p=strtok(argu,",");
                    else p=strtok(NULL,",");
                    for(j=0; j<fdim-1; j++) {
                      if(strcmp(arSet[j].index_name,p)==0) {
                        logiantidim[i1][j]=vars[index].strides[i];
                        for(sup=1; sup<MAXSUPSET; sup++)if(sets[arSet[j].setid].subsetid[sup]==vars[index].setid[i]){logisup[i1][j]=sup;break;}
                        break;
                      }
                    }
                  }
                  break;
                }
              } while (index--);
            }
          }
        #pragma omp parallel private(l,l2,i4,dcount,i3,i1,arSet1,logi,index,eval,ops1) shared(elem_vals,arSet)
        {
        if(omp_get_thread_num()!=0){
          arSet1=realloc(arSet1,(fdim+1)*sizeof(quantifier));
          memcpy (arSet1,arSet,(fdim+1)*sizeof(quantifier));
          ops1=realloc(ops1,(npow+nmul+nplu+2*(npar+2))*sizeof(formula_op));
          memcpy (ops1,ops,(npow+nmul+nplu+2*(npar+2))*sizeof(formula_op));
        }else{
          ops1=ops;
          arSet1=arSet;
        }
        #pragma omp for
          for (l=0; l<nloops; l++) {
            l2=0;
            i4=l;
            for (dcount=0; dcount<fdim-1; dcount++) {
              i3=(offset_t) i4/dcountdim1[dcount];
              arSet1[dcount].indx=i3;
              i4=i4-i3*dcountdim1[dcount];
              if(varsize<fdim-1||varsize==fdim-1) {
                for(i1=0; i1<varsize; i1++) {
                  if(vararset[i1]-1==dcount) {
                    if(varsubset[i1]==1) {
                      l2=l2+set_elems[sets[arSet1[dcount].setid].offset+i3].superset_pos[varsupsetid[i1]]*varantidim[i1];
                    } else {
                      l2=l2+i3*varantidim[i1];
                    }
                    break;
                  }
                }
              } else {
                for(i1=0; i1<varsize; i1++) {
                  if(vararset[i1]-1==dcount) {
                    if(varsubset[i1]==1) {
                      l2=l2+set_elems[sets[arSet1[dcount].setid].offset+i3].superset_pos[varsupsetid[i1]]*varantidim[i1];
                    } else {
                      l2=l2+i3*varantidim[i1];
                    }
                  }
                }
              }
            }
            logi=0;
            index=0;
            for(i1=0; i1<fdim-1; i1++) {
              if(logioper[i1]>0){
              for(i=0; i<fdim-1; i++){
                index+=set_elems[sets[arSet1[i].setid].offset+arSet1[i].indx].superset_pos[logisup[i1][i]]*logiantidim[i1][i];
              }
              if(logivartype[i1]==0)eval=elem_vals[coefs[logivarindx[i1]].offset+index].value;
              else eval=elem_vals[ncofele+vars[logivarindx[i1]].offset+index].value;
              if(logioper[i1]==1)if(eval==cond[i1])logi++;
              if(logioper[i1]==2)if(eval>cond[i1])logi++;
              if(logioper[i1]==3)if(eval<cond[i1])logi++;
              if(logioper[i1]==4)if(eval!=cond[i1])logi++;
              if(logioper[i1]==5)if(eval<=cond[i1])logi++;
              if(logioper[i1]==6)if(eval>=cond[i1])logi++;
              }
            }
            if(logi==ncond)elem_vals[offset+l2].value=formula_eval(elem_vals,sets,set_elems,sum_vals,ops1,nops,arSet1,fdim-1,zerodivide);
          }
        if(omp_get_thread_num()!=0){
          free(arSet1);
          arSet1=NULL;
          free(ops1);
          ops1=NULL;
        }else{
          ops1=NULL;
          arSet1=NULL;
        }
        }
        } else {
        #pragma omp parallel private(l,l2,i4,dcount,i3,i1,arSet1,ops1) shared(elem_vals,arSet)
        {
        if(omp_get_thread_num()!=0){
          arSet1=realloc(arSet1,(fdim+1)*sizeof(quantifier));
          memcpy (arSet1,arSet,(fdim+1)*sizeof(quantifier));
          ops1=realloc(ops1,(npow+nmul+nplu+2*(npar+2))*sizeof(formula_op));
          memcpy (ops1,ops,(npow+nmul+nplu+2*(npar+2))*sizeof(formula_op));
        }else{
          ops1=ops;
          arSet1=arSet;
        }
        #pragma omp for
          for (l=0; l<nloops; l++) {
            l2=0;
            i4=l;
            for (dcount=0; dcount<fdim-1; dcount++) {
              i3=(offset_t) i4/dcountdim1[dcount];
              arSet1[dcount].indx=i3;
              i4=i4-i3*dcountdim1[dcount];
              if(varsize<=fdim-1) {
                for(i1=0; i1<varsize; i1++) {
                  if(vararset[i1]-1==dcount) {
                    if(varsubset[i1]==1) {
                      l2=l2+set_elems[sets[arSet1[dcount].setid].offset+i3].superset_pos[varsupsetid[i1]]*varantidim[i1];
                    } else {
                      l2=l2+i3*varantidim[i1];
                    }
                    break;
                  }
                }
              } else {
                for(i1=0; i1<varsize; i1++) {
                  if(vararset[i1]-1==dcount) {
                    if(varsubset[i1]==1) {
                      l2=l2+set_elems[sets[arSet1[dcount].setid].offset+i3].superset_pos[varsupsetid[i1]]*varantidim[i1];
                    } else {
                      l2=l2+i3*varantidim[i1];
                    }
                  }
                }
              }
            }
            elem_vals[offset+l2].value=formula_eval(elem_vals,sets,set_elems,sum_vals,ops1,nops,arSet1,fdim-1,zerodivide);
          }
        if(omp_get_thread_num()!=0){
          free(arSet1);
          arSet1=NULL;
          free(ops1);
          ops1=NULL;
        }else{
          ops1=NULL;
          arSet1=NULL;
        }
        }
        }

        if(coefs[index].gltype>0){
        #pragma omp parallel private(l) shared(elem_vals,coefs,ncof,offset,index)
        {
          if(coefs[index].gltype==BT_GE){
          #pragma omp for
          for (l=0; l<varsize; l++) {
            if(elem_vals[offset+l].value<coefs[index].glval){
              printf("Error!!! Condition not met var %s type 1!\n",coefs[index].cofname);
              l=varsize;
            }
          }
          }
          if(coefs[index].gltype==BT_GT){
          #pragma omp for
          for (l=0; l<varsize; l++) {
            if(elem_vals[offset+l].value<=coefs[index].glval){
              printf("Error!!! Condition not met var %s type 2!\n",coefs[index].cofname);
              l=varsize;
            }
          }
          }
          if(coefs[index].gltype==BT_LE){
          #pragma omp for
          for (l=0; l<varsize; l++) {
            if(elem_vals[offset+l].value>coefs[index].glval){
              printf("Error!!! Condition not met var %s type 3!\n",coefs[index].cofname);
              l=varsize;
            }
          }
          }
          if(coefs[index].gltype==BT_LT){
          #pragma omp for
          for (l=0; l<varsize; l++) {
            if(elem_vals[offset+l].value>=coefs[index].glval){
              printf("Error!!! Condition not met var %s type 4!\n",coefs[index].cofname);
              l=varsize;
            }
          }
          }
        }
        }
        
        
        
        free(sum_cof);
        free(sum_vals);
        free(arSet);
        free(ops);
      }
    }
  }
  fclose(filehandle);
  return j;
}

/* Applies the TAB update statements to coefficient values. With
   midpoint!=0 the modified-midpoint correction is used: the value is
   advanced from the sub-step base (csolpupd) by twice the computed
   change, and csolpupd retains the pre-update value. */
offset_t updates_apply(char *fname,set_def *sets,dim_t nset, set_element *set_elems, array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar, elem_value *elem_vals,offset_t ncofvar,offset_t ncofele,int midpoint) {
  FILE * filehandle;
  char commsyntax[NAMESIZE],line[TABREADLINE],line1[TABREADLINE],line2[TABREADLINE],linecopy[TABREADLINE];
  char vname[NAMESIZE],sumsyntax[NAMESIZE],argu[NAMESIZE];
  char *readitem=NULL,*p=NULL;
  offset_t i,i1,i3,i4,l,l2=0,j=0,nsumele,dcountdim1[4*MAXVARDIM],nloops;//m,
  dim_t fdim,dcount,sup,varsupsetid[MAXVARDIM];
  int nops=0,totalsum,sumcount=1,npow,nmul,ndiv,nplu,nmin,npar,sumindx;
  bool IsChange=false,IsExplicit=false;
  solve_real zerodivide=0,temp1,temp2;
  offset_t varantidim[MAXVARDIM],varsubset[MAXVARDIM],vararset[MAXVARDIM];
  quantifier *arSet1=NULL;
  formula_op *ops1= NULL;
  strcpy(commsyntax,"update");
  filehandle = fopen(fname,"r");
  while (tab_next_statement_resolved(commsyntax,filehandle,line,elem_vals,coefs,ncof,&zerodivide,TABREADLINE)) {
    IsChange=false;
    IsExplicit=false;
    if(strstr(line, "(change)")!=NULL) {
      IsChange=true;
      str_replace_first(line, "(change)", "");
    }
    if(strstr(line, "(explicit)")!=NULL) {
      IsExplicit=true;
      str_replace_first(line, "(explicit)", "");
    }
    str_replace_first(line, commsyntax, "");
    while (str_replace_all(line," ", ""));
    while (str_replace_char(line, '[', '('));
    while (str_replace_char(line, ']', ')'));
    while (str_replace_char(line, '{', '('));
    while (str_replace_char(line, '}', ')'));
    strcpy(linecopy,line);
    readitem = strtok(line,"=");
    fdim=str_count_char(readitem, '(');
    if (fdim==1) {
      fdim=fdim+1;
    }
    quantifier *arSet= (quantifier *) calloc (fdim+1,sizeof(quantifier));

    nloops=1;
    if (fdim==0) {
      readitem = strtok(line,"=");
    } else {
      for (i=0; i<fdim-1; i++) {
        if(i==0) {
          readitem = strtok(line,",");
        } else {
          readitem = strtok(NULL,",");
        }
        readitem = strtok(NULL,",");
        strcpy(arSet[i].index_name,readitem);
        readitem = strtok(NULL,")");
        for (i4=0; i4<nset; i4++) if(strcmp(readitem,sets[i4].setname)==0) {
            arSet[i].setid=i4;
            break;
          }
        nloops=nloops*sets[arSet[i].setid].size;
      }
      readitem = strtok(NULL,"=");
      dcountdim1[fdim-2]=1;
      for (i=fdim-3; i>-1; i--) {
        dcountdim1[i]=sets[arSet[i+1].setid].size*dcountdim1[i+1];
      }
    }
    strcpy(vname,readitem);
    strcpy(line,linecopy);
    readitem = strtok(line,"=");
    readitem = strtok(NULL,";");
    if(IsChange==false&&IsExplicit==false) {
      while (str_replace_all(readitem,"*", "+"));
      strcpy(line1,vname);
      strcat(line1,"*(1+(");
      strcat(line1,readitem);
      strcat(line1,")/100)");
      readitem=line1;
    }
    if(IsChange==true) {
      strcpy(line1,vname);
      if(readitem[0]!='+'||readitem[0]!='-') strcat(line1,"+");
      strcat(line1,readitem);
      readitem=line1;
    }
    npow=str_count_char(readitem, '^');
    nmul=str_count_char(readitem, '*');
    ndiv=str_count_char(readitem, '/');
    nmul=nmul+ndiv;
    nplu=str_count_char(readitem, '+');
    nmin=str_count_char(readitem, '-');
    nplu=nplu+nmin;

    while (formula_normalize(readitem)==1);
    leadlag_encode(readitem);
    npar=str_count_char(readitem, '(');
    strcpy(sumsyntax,"sum(");
    totalsum=sum_count(readitem,sumsyntax);
    sum_def *sum_cof= (sum_def *) calloc (totalsum,sizeof(sum_def));
    sumcount=0;
    strcpy(line1,readitem);
    strcpy(line2,line1);
    readitem=line2;
    while (sum_parse(readitem,sumsyntax,sum_cof,arSet,sets,nset,fdim,sumcount)==1) {
      sumcount++;
    }
    totalsum=sumcount;
    i3=0;
    for (i=0; i<totalsum; i++) {
      i1=1;
      for(j=0; j<sum_cof[i].size; j++) {
        i1=i1*sets[sum_cof[i].setid[j]].size;
      }
      sum_cof[i].offset=i3;
      i3=i3+i1;
    }
    nsumele=i3;
    for (i=0; i<totalsum; i++) {
      i1=1;
      sum_cof[i].strides[sum_cof[i].size-1]=1;
      for(j=sum_cof[i].size-2; j>-1; j--) {
        sum_cof[i].strides[j]=sum_cof[i].strides[j+1]*sets[sum_cof[i].setid[j+1]].size;
      }
    }
    formula_op *ops= (formula_op *) calloc (npow+nmul+nplu+2*npar+2,sizeof(formula_op));
    sum_value *sum_vals= (sum_value *) calloc (nsumele,sizeof(sum_value));
    sumcount=0;
    strcpy(line2,line1);
    readitem=line2;
    sumindx=0;
    while (sum_eval(readitem,sumsyntax,sets,nset,set_elems,elem_vals,ncofvar,ncofele,coefs,ncof,vars,nvar,sum_cof,totalsum,sum_vals,nsumele,ops,arSet,fdim,&sumindx,sumcount,zerodivide)==1) {
      sumcount++;
    }
    strcpy(line1,readitem);
    offset_t index=ncof-1, offset=0;//,simpl=0;
    bool check10=true;
        offset_t varsize=0;
        p=strtok(vname,"(");
        do {
          if (strcmp(coefs[index].cofname,p)==0) {
            offset=coefs[index].offset;
            varsize=coefs[index].size;
            if(coefs[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
            check10=false;
            break;
          }
        } while (index--);
        if (check10) {
          index=nvar-1;
          do {
            if (strcmp(vars[index].cofname,p)==0) {
              offset=ncofele+vars[index].offset;
              varsize=vars[index].size;
              if(vars[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
              break;
            }
          } while (index--);
        }
        for (l=0; l<MAXVARDIM; l++){varantidim[l]=0;varsubset[l]=0;varsupsetid[l]=0;}
        if (check10) {
            for (dcount=0; dcount<vars[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].index_name,p)==0) {
                  varantidim[dcount]=vars[index].strides[dcount];
                  vararset[dcount]=l+1;
                  if (sets[arSet[l].setid].size!=sets[vars[index].setid[dcount]].size){
                    varsubset[dcount]=1;
                    for(sup=1;sup<MAXSUPSET;sup++)if(sets[arSet[l].setid].subsetid[sup]==vars[index].setid[dcount]){varsupsetid[dcount]=sup;break;}
                  }
                  break;
                }
            }
        } else {
            for (dcount=0; dcount<coefs[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].index_name,p)==0) {
                  varantidim[dcount]=coefs[index].strides[dcount];
                  vararset[dcount]=l+1;
                  if (sets[arSet[l].setid].size!=sets[coefs[index].setid[dcount]].size){
                    varsubset[dcount]=1;
                    for(sup=1;sup<MAXSUPSET;sup++)if(sets[arSet[l].setid].subsetid[sup]==coefs[index].setid[dcount]){varsupsetid[dcount]=sup;break;}
                  }
                  break;
                }
            }
        }
    formula_compile(line1,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,&nops,arSet,fdim-1);
        #pragma omp parallel private(l,l2,i4,dcount,i3,i1,temp1,temp2,arSet1,ops1) shared(elem_vals,arSet)
        {
        if(omp_get_thread_num()!=0){
          arSet1=realloc(arSet1,(fdim+1)*sizeof(quantifier));
          memcpy (arSet1,arSet,(fdim+1)*sizeof(quantifier));
          ops1=realloc(ops1,nops*sizeof(formula_op));
          memcpy (ops1,ops,nops*sizeof(formula_op));
        }else{
          ops1=ops;
          arSet1=arSet;
        }
        #pragma omp for
    for (l=0; l<nloops; l++) {
      l2=0;
      i4=l;
            for (dcount=0; dcount<fdim-1; dcount++) {
              i3=(offset_t) i4/dcountdim1[dcount];
              arSet1[dcount].indx=i3;
              i4=i4-i3*dcountdim1[dcount];
              if(varsize<=fdim-1) {
                for(i1=0; i1<varsize; i1++) {
                  if(vararset[i1]-1==dcount) {
                    if(varsubset[i1]==1) {
                      l2=l2+set_elems[sets[arSet1[dcount].setid].offset+i3].superset_pos[varsupsetid[i1]]*varantidim[i1];
                    } else {
                      l2=l2+i3*varantidim[i1];
                    }
                    break;
                  }
                }
              } else {
                for(i1=0; i1<varsize; i1++) {
                  if(vararset[i1]-1==dcount) {
                    if(varsubset[i1]==1) {
                      l2=l2+set_elems[sets[arSet1[dcount].setid].offset+i3].superset_pos[varsupsetid[i1]]*varantidim[i1];
                    } else {
                      l2=l2+i3*varantidim[i1];
                    }
                  }
                }
              }
            }
      if(midpoint){
        temp2=elem_vals[offset+l2].value;
        temp1=formula_eval(elem_vals,sets,set_elems,sum_vals,ops1,nops,arSet1,fdim-1,zerodivide);
        if(temp1-elem_vals[offset+l2].value>0.000000001||temp1-elem_vals[offset+l2].value<-0.000000001)elem_vals[offset+l2].value=elem_vals[offset+l2].substep_base+2*(temp1-elem_vals[offset+l2].value);
        elem_vals[offset+l2].substep_base=temp2;
      }else{
        elem_vals[offset+l2].substep_base=elem_vals[offset+l2].value;
        temp1=formula_eval(elem_vals,sets,set_elems,sum_vals,ops1,nops,arSet1,fdim-1,zerodivide);
        if(temp1-elem_vals[offset+l2].value>0.000000001||temp1-elem_vals[offset+l2].value<-0.000000001)elem_vals[offset+l2].value=temp1;
      }
    }
        if(omp_get_thread_num()!=0){
          free(arSet1);
          arSet1=NULL;
          free(ops1);
          ops1=NULL;
        }else{
          ops1=NULL;
          arSet1=NULL;
        }
        }
    free(sum_cof);
    free(sum_vals);
    free(arSet);
    free(ops);

        if(coefs[index].gltype>0){
        #pragma omp parallel private(l) shared(elem_vals,coefs,ncof,offset,index)
        {
          if(coefs[index].gltype==BT_GE){
          #pragma omp for
          for (l=0; l<varsize; l++) {
            if(elem_vals[offset+l].value<coefs[index].glval){
              printf("Error!!! Condition not met var %s type 1!\n",coefs[index].cofname);
              l=varsize;
            }
          }
          }
          if(coefs[index].gltype==BT_GT){
          #pragma omp for
          for (l=0; l<varsize; l++) {
            if(elem_vals[offset+l].value<=coefs[index].glval){
              printf("Error!!! Condition not met var %s type 2!\n",coefs[index].cofname);
              l=varsize;
            }
          }
          }
          if(coefs[index].gltype==BT_LE){
          #pragma omp for
          for (l=0; l<varsize; l++) {
            if(elem_vals[offset+l].value>coefs[index].glval){
              printf("Error!!! Condition not met var %s type 3!\n",coefs[index].cofname);
              l=varsize;
            }
          }
          }
          if(coefs[index].gltype==BT_LT){
          #pragma omp for
          for (l=0; l<varsize; l++) {
            if(elem_vals[offset+l].value>=coefs[index].glval){
              printf("Error!!! Condition not met var %s type 4!\n",coefs[index].cofname);
              l=varsize;
            }
          }
          }
        }
        }
    
  }
  fclose(filehandle);
  return j;
}



offset_t updates_apply_product(char *fname,set_def *sets,dim_t nset, set_element *set_elems, array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar, elem_value *elem_vals,offset_t ncofvar,offset_t ncofele) {
  FILE * filehandle;
  char commsyntax[NAMESIZE],line[TABREADLINE],line1[TABREADLINE],line2[TABREADLINE],linecopy[TABREADLINE];
  char vname[NAMESIZE],sumsyntax[NAMESIZE],argu[NAMESIZE];
  char *readitem=NULL,*p=NULL;
  offset_t i,i1,i3,i4,l,l2=0,j=0,nsumele,dcountdim1[4*MAXVARDIM],nloops;//m,
  dim_t fdim,dcount,sup,varsupsetid[MAXVARDIM];
  int nops=0,totalsum,sumcount=1,npow,nmul,ndiv,nplu,nmin,npar,sumindx;
  bool IsChange=false,IsExplicit=false;
  solve_real zerodivide=0,temp1;
  offset_t varantidim[MAXVARDIM],varsubset[MAXVARDIM],vararset[MAXVARDIM];
  quantifier *arSet1=NULL;
  formula_op *ops1= NULL;
  strcpy(commsyntax,"update");
  filehandle = fopen(fname,"r");
  while (tab_next_statement_resolved(commsyntax,filehandle,line,elem_vals,coefs,ncof,&zerodivide,TABREADLINE)) {
    IsChange=false;
    IsExplicit=false;
    if(strstr(line, "(change)")!=NULL) {
      IsChange=true;
      str_replace_first(line, "(change)", "");
    }
    if(strstr(line, "(explicit)")!=NULL) {
      IsExplicit=true;
      str_replace_first(line, "(explicit)", "");
    }
    str_replace_first(line, commsyntax, "");
    while (str_replace_all(line," ", ""));
    while (str_replace_char(line, '[', '('));
    while (str_replace_char(line, ']', ')'));
    while (str_replace_char(line, '{', '('));
    while (str_replace_char(line, '}', ')'));
    strcpy(linecopy,line);
    readitem = strtok(line,"=");
    fdim=str_count_char(readitem, '(');
    if (fdim==1) {
      fdim=fdim+1;
    }
    quantifier *arSet= (quantifier *) calloc (fdim+1,sizeof(quantifier));

    nloops=1;
    if (fdim==0) {
      readitem = strtok(line,"=");
    } else {
      for (i=0; i<fdim-1; i++) {
        if(i==0) {
          readitem = strtok(line,",");
        } else {
          readitem = strtok(NULL,",");
        }
        readitem = strtok(NULL,",");
        strcpy(arSet[i].index_name,readitem);
        readitem = strtok(NULL,")");
        for (i4=0; i4<nset; i4++) if(strcmp(readitem,sets[i4].setname)==0) {
            arSet[i].setid=i4;
            break;
          }
        nloops=nloops*sets[arSet[i].setid].size;
      }
      readitem = strtok(NULL,"=");
      dcountdim1[fdim-2]=1;
      for (i=fdim-3; i>-1; i--) {
        dcountdim1[i]=sets[arSet[i+1].setid].size*dcountdim1[i+1];
      }
    }
    strcpy(vname,readitem);
    strcpy(line,linecopy);
    readitem = strtok(line,"=");
    readitem = strtok(NULL,";");
    if(IsChange==false&&IsExplicit==false) {
      while (str_replace_all(readitem,"*", "/100)!(1+"));
      while (str_replace_all(readitem,"!", "*"));
      strcpy(line1,vname);
      strcat(line1,"*(1+");
      strcat(line1,readitem);
      strcat(line1,"/100)");
      readitem=line1;
    }
    if(IsChange==true) {
      strcpy(line1,vname);
      if(readitem[0]!='+'||readitem[0]!='-') strcat(line1,"+");
      strcat(line1,readitem);
      readitem=line1;
    }
    npow=str_count_char(readitem, '^');
    nmul=str_count_char(readitem, '*');
    ndiv=str_count_char(readitem, '/');
    nmul=nmul+ndiv;
    nplu=str_count_char(readitem, '+');
    nmin=str_count_char(readitem, '-');
    nplu=nplu+nmin;

    while (formula_normalize(readitem)==1);
    leadlag_encode(readitem);
    npar=str_count_char(readitem, '(');
    strcpy(sumsyntax,"sum(");
    totalsum=sum_count(readitem,sumsyntax);
    sum_def *sum_cof= (sum_def *) calloc (totalsum,sizeof(sum_def));
    sumcount=0;
    strcpy(line1,readitem);
    strcpy(line2,line1);
    readitem=line2;
    while (sum_parse(readitem,sumsyntax,sum_cof,arSet,sets,nset,fdim,sumcount)==1) {
      sumcount++;
    }
    totalsum=sumcount;
    i3=0;
    for (i=0; i<totalsum; i++) {
      i1=1;
      for(j=0; j<sum_cof[i].size; j++) {
        i1=i1*sets[sum_cof[i].setid[j]].size;
      }
      sum_cof[i].offset=i3;
      i3=i3+i1;
    }
    nsumele=i3;
    for (i=0; i<totalsum; i++) {
      i1=1;
      sum_cof[i].strides[sum_cof[i].size-1]=1;
      for(j=sum_cof[i].size-2; j>-1; j--) {
        sum_cof[i].strides[j]=sum_cof[i].strides[j+1]*sets[sum_cof[i].setid[j+1]].size;
      }
    }
    formula_op *ops= (formula_op *) calloc (npow+nmul+nplu+2*npar+2,sizeof(formula_op));
    sum_value *sum_vals= (sum_value *) calloc (nsumele,sizeof(sum_value));
    sumcount=0;
    strcpy(line2,line1);
    readitem=line2;
    sumindx=0;
    while (sum_eval(readitem,sumsyntax,sets,nset,set_elems,elem_vals,ncofvar,ncofele,coefs,ncof,vars,nvar,sum_cof,totalsum,sum_vals,nsumele,ops,arSet,fdim,&sumindx,sumcount,zerodivide)==1) {
      sumcount++;
    }
    strcpy(line1,readitem);
    offset_t index=ncof-1, offset=0;//,simpl=0;
    bool check10=true;
        offset_t varsize=0;
        p=strtok(vname,"(");
        do {
          if (strcmp(coefs[index].cofname,p)==0) {
            offset=coefs[index].offset;
            varsize=coefs[index].size;
            if(coefs[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
            check10=false;
            break;
          }
        } while (index--);
        if (check10) {
          index=nvar-1;
          do {
            if (strcmp(vars[index].cofname,p)==0) {
              offset=ncofele+vars[index].offset;
              varsize=vars[index].size;
              if(vars[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
              break;
            }
          } while (index--);
        }
        for (l=0; l<MAXVARDIM; l++){varantidim[l]=0;varsubset[l]=0;varsupsetid[l]=0;}
        if (check10) {
            for (dcount=0; dcount<vars[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].index_name,p)==0) {
                  varantidim[dcount]=vars[index].strides[dcount];
                  vararset[dcount]=l+1;
                  if (sets[arSet[l].setid].size!=sets[vars[index].setid[dcount]].size){
                    varsubset[dcount]=1;
                    for(sup=1;sup<MAXSUPSET;sup++)if(sets[arSet[l].setid].subsetid[sup]==vars[index].setid[dcount]){varsupsetid[dcount]=sup;break;}
                  }
                  break;
                }
            }
        } else {
            for (dcount=0; dcount<coefs[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].index_name,p)==0) {
                  varantidim[dcount]=coefs[index].strides[dcount];
                  vararset[dcount]=l+1;
                  if (sets[arSet[l].setid].size!=sets[coefs[index].setid[dcount]].size){
                    varsubset[dcount]=1;
                    for(sup=1;sup<MAXSUPSET;sup++)if(sets[arSet[l].setid].subsetid[sup]==coefs[index].setid[dcount]){varsupsetid[dcount]=sup;break;}
                  }
                  break;
                }
            }
        }
    formula_compile(line1,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,&nops,arSet,fdim-1);
        #pragma omp parallel private(l,l2,i4,dcount,i3,i1,temp1,arSet1,ops1) shared(elem_vals,arSet)
        {
        if(omp_get_thread_num()!=0){
          arSet1=realloc(arSet1,(fdim+1)*sizeof(quantifier));
          memcpy (arSet1,arSet,(fdim+1)*sizeof(quantifier));
          ops1=realloc(ops1,nops*sizeof(formula_op));
          memcpy (ops1,ops,nops*sizeof(formula_op));
        }else{
          ops1=ops;
          arSet1=arSet;
        }
        #pragma omp for
    for (l=0; l<nloops; l++) {
      l2=0;
      i4=l;
            for (dcount=0; dcount<fdim-1; dcount++) {
              i3=(offset_t) i4/dcountdim1[dcount];
              arSet1[dcount].indx=i3;
              i4=i4-i3*dcountdim1[dcount];
              if(varsize<=fdim-1) {
                for(i1=0; i1<varsize; i1++) {
                  if(vararset[i1]-1==dcount) {
                    if(varsubset[i1]==1) {
                      l2=l2+set_elems[sets[arSet1[dcount].setid].offset+i3].superset_pos[varsupsetid[i1]]*varantidim[i1];
                    } else {
                      l2=l2+i3*varantidim[i1];
                    }
                    break;
                  }
                }
              } else {
                for(i1=0; i1<varsize; i1++) {
                  if(vararset[i1]-1==dcount) {
                    if(varsubset[i1]==1) {
                      l2=l2+set_elems[sets[arSet1[dcount].setid].offset+i3].superset_pos[varsupsetid[i1]]*varantidim[i1];
                    } else {
                      l2=l2+i3*varantidim[i1];
                    }
                  }
                }
              }
            }
      temp1=formula_eval(elem_vals,sets,set_elems,sum_vals,ops1,nops,arSet1,fdim-1,zerodivide);
      if(temp1-elem_vals[offset+l2].value>0.000000001||temp1-elem_vals[offset+l2].value<-0.000000001)elem_vals[offset+l2].value=temp1;
    }
        if(omp_get_thread_num()!=0){
          free(arSet1);
          arSet1=NULL;
          free(ops1);
          ops1=NULL;
        }else{
          ops1=NULL;
          arSet1=NULL;
        }
        }
    free(sum_cof);
    free(sum_vals);
    free(arSet);
    free(ops);

        if(coefs[index].gltype>0){
        #pragma omp parallel private(l) shared(elem_vals,coefs,ncof,offset,index)
        {
          if(coefs[index].gltype==BT_GE){
          #pragma omp for
          for (l=0; l<varsize; l++) {
            if(elem_vals[offset+l].value<coefs[index].glval){
              printf("Error!!! Condition not met var %s type 1!\n",coefs[index].cofname);
              l=varsize;
            }
          }
          }
          if(coefs[index].gltype==BT_GT){
          #pragma omp for
          for (l=0; l<varsize; l++) {
            if(elem_vals[offset+l].value<=coefs[index].glval){
              printf("Error!!! Condition not met var %s type 2!\n",coefs[index].cofname);
              l=varsize;
            }
          }
          }
          if(coefs[index].gltype==BT_LE){
          #pragma omp for
          for (l=0; l<varsize; l++) {
            if(elem_vals[offset+l].value>coefs[index].glval){
              printf("Error!!! Condition not met var %s type 3!\n",coefs[index].cofname);
              l=varsize;
            }
          }
          }
          if(coefs[index].gltype==BT_LT){
          #pragma omp for
          for (l=0; l<varsize; l++) {
            if(elem_vals[offset+l].value>=coefs[index].glval){
              printf("Error!!! Condition not met var %s type 4!\n",coefs[index].cofname);
              l=varsize;
            }
          }
          }
        }
        }
    
  }
  fclose(filehandle);
  return j;
}

int sum_eval(char *formulain, char *commsyntax,set_def *sets,dim_t nset, set_element *set_elems,elem_value *elem_vals,offset_t ncofvar,offset_t ncofele, array_def *coefs,offset_t ncof, array_def *vars,offset_t nvar,sum_def *sum_cof,int totalsum,sum_value *sum_vals,offset_t nsumele,formula_op *ops,quantifier *arSet1,dim_t fdim,int *sumindx,int j, solve_real zerodivide) {
  char *readitem,*p;//,*p1,interchar2[NAMESIZE],line5[TABREADLINE];
  char interchar[NAMESIZE],line[TABREADLINE],line1[TABREADLINE],line2[TABREADLINE];//,line3[TABREADLINE],line4[TABREADLINE];//,interchar1[NAMESIZE]
  int nops,length,k=0,k1=0,i=0;
  dim_t fdimsumcof,superset_pos,dcount;
  offset_t l,l1,l2,nloops,dcountdim1[4*MAXVARDIM];
  solve_real vval;
  quantifier *arSet2=NULL;
  formula_op *ops1= NULL;
  offset_t arsetsize;
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
      sum_extract(line);
      k1=str_find_ci(line+4,commsyntax);
      if (k1!=-1) {
        i=i+k+4;
        readitem=formulain+i;
      } else {
        strcpy(line1,line);
        p=strtok(line,",");
        p=strtok(NULL,",");
        p=strtok(NULL,"\0");
        p[strlen(p)-1]='\0';
        strcpy(line2,p);
        arsetsize=sum_cof[j].size+1;
        quantifier *arSet= (quantifier *) calloc (arsetsize,sizeof(quantifier));
        for (l=0; l<sum_cof[j].size; l++) {
          arSet[l].setid=sum_cof[j].setid[l];
          strcpy(arSet[l].index_name,sum_cof[j].dimnames[l]);
        }
        nloops=1;
        for (l=0; l<sum_cof[j].size; l++) {
          nloops=nloops*sets[arSet[l].setid].size;
          dcount=sum_cof[j].size-l;
          if(dcount==sum_cof[j].size) {
            dcountdim1[dcount-1]=1;
          } else {
            dcountdim1[dcount-1]=dcountdim1[dcount]*sets[arSet[dcount].setid].size;
          }
        }
        arSet[sum_cof[j].size].setid=sum_cof[j].sumsetid;
        strcpy(arSet[sum_cof[j].size].index_name,sum_cof[j].sumindx);
        fdimsumcof=sum_cof[j].size+1;
        nops=0;
        formula_compile(p,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,&nops,arSet,fdimsumcof);
        #pragma omp parallel private(l,l1,l2,dcount,superset_pos,vval,arSet2,ops1) shared(elem_vals,arSet,sum_vals)
        {
        if(omp_get_thread_num()!=0){
          arSet2=realloc(arSet2,arsetsize*sizeof(quantifier));
          memcpy(arSet2,arSet,arsetsize*sizeof(quantifier));
          ops1=realloc(ops1,nops*sizeof(formula_op));
          memcpy(ops1,ops,nops*sizeof(formula_op));
        }else{
          ops1=ops;
          arSet2=arSet;
        }
        #pragma omp for
        for (l=0; l<nloops; l++) {
          l2=l;
          for (dcount=0; dcount<sum_cof[j].size; dcount++) {
            superset_pos=(dim_t) l2/dcountdim1[dcount];
            arSet2[dcount].indx=superset_pos;
            l2=l2-superset_pos*dcountdim1[dcount];
          }
          vval=0;
          for (l1=0; l1<sets[sum_cof[j].sumsetid].size; l1++) {
            arSet2[sum_cof[j].size].indx=l1;
            vval+=formula_eval(elem_vals,sets,set_elems,sum_vals,ops1,nops,arSet2,fdimsumcof,zerodivide);
          }
          sum_vals[*sumindx+l].value=vval;
        }
        if(omp_get_thread_num()!=0){
          free(arSet2);
          arSet2=NULL;
          free(ops1);
          ops1=NULL;
        }else{
          ops1=NULL;
          arSet2=NULL;
        }
        }
        *sumindx=*sumindx+nloops;
        strcpy(interchar,sum_cof[j].sumname);
        strcat(interchar,"{");
        for (l=0; l<sum_cof[j].size; l++) {
          strcat(interchar,sum_cof[j].dimnames[l]);
          strcat(interchar,",");
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
        while(str_replace_all(formulain,line1,interchar)!=NULL);
        free(arSet);
        return 1;
      }
    } else if (formulain[i+k-1]=='+'||formulain[i+k-1]=='-'||formulain[i+k-1]=='*'||formulain[i+k-1]=='/'||formulain[i+k-1]=='^'||formulain[i+k-1]=='('||formulain[i+k-1]==',') {
      readitem=formulain+i+k;
      strcpy(line,readitem);
      sum_extract(line);
      k1=str_find_ci(line+4,commsyntax);
      if (k1!=-1) {
        i=i+k+4;
        readitem=formulain+i;
      } else {
        strcpy(line1,line);
        p=strtok(line,",");
        p=strtok(NULL,",");
        p=strtok(NULL,"\0");
        p[strlen(p)-1]='\0';
        arsetsize=sum_cof[j].size+1;
        quantifier *arSet= (quantifier *) calloc (arsetsize,sizeof(quantifier));
        for (l=0; l<sum_cof[j].size; l++) {
          arSet[l].setid=sum_cof[j].setid[l];
          strcpy(arSet[l].index_name,sum_cof[j].dimnames[l]);
        }
        nloops=1;
        for (l=0; l<sum_cof[j].size; l++) {
          nloops=nloops*sets[arSet[l].setid].size;//sum_cof[j].dims[l];
          dcount=sum_cof[j].size-l;
          if(dcount==sum_cof[j].size) {
            dcountdim1[dcount-1]=1;
          } else {
            dcountdim1[dcount-1]=dcountdim1[dcount]*sets[arSet[dcount].setid].size;
          }
        }
        arSet[sum_cof[j].size].setid=sum_cof[j].sumsetid;
        strcpy(arSet[sum_cof[j].size].index_name,sum_cof[j].sumindx);
        fdimsumcof=sum_cof[j].size+1;
        nops=0;
        formula_compile(p,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,&nops,arSet,fdimsumcof);
        #pragma omp parallel private(l,l1,l2,dcount,superset_pos,vval,arSet2,ops1) shared(elem_vals,arSet,sum_vals)
        {
        if(omp_get_thread_num()!=0){
          arSet2=realloc(arSet2,arsetsize*sizeof(quantifier));
          memcpy(arSet2,arSet,arsetsize*sizeof(quantifier));
          ops1=realloc(ops1,nops*sizeof(formula_op));
          memcpy(ops1,ops,nops*sizeof(formula_op));
        }else{
          ops1=ops;
          arSet2=arSet;
        }
        #pragma omp for
        for (l=0; l<nloops; l++) {
          l2=l;
          for (dcount=0; dcount<sum_cof[j].size; dcount++) {
            superset_pos=(dim_t) l2/dcountdim1[dcount];
            arSet2[dcount].indx=superset_pos;
            l2=l2-superset_pos*dcountdim1[dcount];
          }
          vval=0;
          for (l1=0; l1<sets[sum_cof[j].sumsetid].size; l1++) {
            arSet2[sum_cof[j].size].indx=l1;
            vval+=formula_eval(elem_vals,sets,set_elems,sum_vals,ops1,nops,arSet2,fdimsumcof,zerodivide);
          }
          sum_vals[*sumindx+l].value=vval;//ha_sumele[*sumindx+l2].varval=vval;
        }
        if(omp_get_thread_num()!=0){
          free(arSet2);
          arSet2=NULL;
          free(ops1);
          ops1=NULL;
        }else{
          ops1=NULL;
          arSet2=NULL;
        }
        }
        *sumindx=*sumindx+nloops;
        strcpy(interchar,sum_cof[j].sumname);
        strcat(interchar,"{");
        for (l=0; l<sum_cof[j].size; l++) {
          strcat(interchar,sum_cof[j].dimnames[l]);
          strcat(interchar,",");
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
        while(str_replace_all(formulain,line1,interchar));
        free(arSet);
        return 1;
      }
    } else {
      i=i+k+4;
      readitem=formulain+i;
    }
  }
  return 0;
}

offset_t subinterval_update(PetscInt rank,char *fname,set_def *sets,dim_t nset, set_element *set_elems, array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar, elem_value *elem_vals,offset_t ncofvar,offset_t ncofele,closure_entry *closure_vals,offset_t nvarele,int laA,dim_t subints,bool IsIni,int IsSplint,int nsteps) {
  FILE * filehandle,*fout;
  char commsyntax[NAMESIZE],line[TABREADLINE],line1[TABREADLINE],line2[TABREADLINE],linecopy[TABREADLINE];
  char vname[NAMESIZE],sumsyntax[NAMESIZE],argu[NAMESIZE];
  char *readitem=NULL,*p=NULL;
  char filename[1024],j1name[1024];
  offset_t i,i1,i3,i4,l,l2=0,j=0,nsumele,dcountdim1[4*MAXVARDIM],nloops,nloops1,tsize,stosize,bsize,rsize,vvar,xvar,yvar,jvar,wvar,svar,vlmu,dlmu,wlmu,varsize,l2v;//m,
  dim_t fdim,fdim1,dcount,sup;
  bool IsChange=false;
  offset_t index=ncof-1, offset=0,sbegadd=0,address=0;//,simpl=0;
  bool check10=true;
  solve_real zerodivide=0,temp1,temp2,temp3,maxerr=0,temp4,temp5,sx0,sxn;
  offset_t varantidim[MAXVARDIM],varsubset[MAXVARDIM],vararset[MAXVARDIM],varsupsetid[MAXVARDIM];
  offset_t xantidim[MAXVARDIM],xsubset[MAXVARDIM],xarset[MAXVARDIM],xsupsetid[MAXVARDIM];
  offset_t yantidim[MAXVARDIM],ysubset[MAXVARDIM],yarset[MAXVARDIM],ysupsetid[MAXVARDIM];
  offset_t jantidim[MAXVARDIM],jsubset[MAXVARDIM],jarset[MAXVARDIM],jsupsetid[MAXVARDIM];
  offset_t wantidim[MAXVARDIM],wsubset[MAXVARDIM],warset[MAXVARDIM],wsupsetid[MAXVARDIM];
  offset_t santidim[MAXVARDIM],ssubset[MAXVARDIM],sarset[MAXVARDIM],ssupsetid[MAXVARDIM];
  offset_t mvantidim[MAXVARDIM],mvsubset[MAXVARDIM],mvarset[MAXVARDIM],mvsupsetid[MAXVARDIM];
  offset_t mdantidim[MAXVARDIM],mdsubset[MAXVARDIM],mdarset[MAXVARDIM],mdsupsetid[MAXVARDIM];
  offset_t mwantidim[MAXVARDIM],mwsubset[MAXVARDIM],mwarset[MAXVARDIM],mwsupsetid[MAXVARDIM];
  size_t freadresult;
  quantifier *arSet1=NULL;
  formula_op *ops1= NULL;
  strcpy(commsyntax,"splinter");
      if(rank<10)strcpy(j1name,"000");
      if(rank<100&&rank>9)strcpy(j1name,"00");
      if(rank<1000&&rank>99)strcpy(j1name,"0");
      if(rank>=1000)j1name[0]='\0';
      sprintf(filename, "%d",rank);
      strcat(j1name,filename);
      strcpy(filename,"_biupd");
      strcat(filename,j1name);
      strcat(filename,".bin");
  filehandle = fopen(fname,"r");
  while (tab_next_statement_resolved(commsyntax,filehandle,line,elem_vals,coefs,ncof,&zerodivide,TABREADLINE)) {
    IsChange=false;
    if(strstr(line, "(change)")!=NULL) {
      IsChange=true;
      str_replace_first(line, "(change)", "");
    }
    str_replace_first(line, commsyntax, "");
    while (str_replace_all(line," ", ""));
    while (str_replace_char(line, '[', '('));
    while (str_replace_char(line, ']', ')'));
    while (str_replace_char(line, '{', '('));
    while (str_replace_char(line, '}', ')'));
    strcpy(linecopy,line);
    readitem = strtok(line,"=");//Syntax (linear or change)(other sets):(inter set):(spline set)var_plus()=var():x():x_jump():weight();//
    fdim=str_count_char(readitem, '(');
    if (fdim==1) {
      fdim=fdim+1;
    }
    fdim1=fdim-1;
    quantifier *arSet= (quantifier *) calloc (fdim+1,sizeof(quantifier));
    if(fdim<3){printf("We should at least have an intertemporal variable and a spline set!\n");return -1;}
    nloops=1;
    nloops1=1;
    if (fdim==3) {
      readitem = strtok(line,"=");
    } else {
      for (i=0; i<fdim-1; i++) {
        if(i==0) {
          readitem = strtok(line,",");
        } else {
          readitem = strtok(NULL,",");
        }
        readitem = strtok(NULL,",");
        strcpy(arSet[i].index_name,readitem);
        readitem = strtok(NULL,")");
        for (i4=0; i4<nset; i4++) if(strcmp(readitem,sets[i4].setname)==0) {
            arSet[i].setid=i4;
            break;
          }
        if(i<fdim-3)nloops=nloops*sets[arSet[i].setid].size;
        if(i>=fdim-3)nloops1=nloops1*sets[arSet[i].setid].size;
      }
      readitem = strtok(NULL,"=");
      dcountdim1[fdim-4]=1;
      for (i=fdim-5; i>-1; i--) {
        dcountdim1[i]=sets[arSet[i+1].setid].size*dcountdim1[i+1];
      }
    }
    printf("nloop %ld nloop1 %ld\n",nloops,nloops1);
    strcpy(vname,readitem);
        index=ncof-1;
        p=strtok(vname,"(");
        do {
          if (strcmp(coefs[index].cofname,p)==0) {
            vvar=index;
            if(coefs[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
            break;
          }
        } while (index--);
        for (l=0; l<MAXVARDIM; l++){varantidim[l]=0;varsubset[l]=0;varsupsetid[l]=0;}
            for (dcount=0; dcount<coefs[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].index_name,p)==0) {
                  varantidim[dcount]=coefs[index].strides[dcount];
                  vararset[dcount]=l+1;
                  if (sets[arSet[l].setid].size!=sets[coefs[index].setid[dcount]].size){
                    varsubset[dcount]=1;
                    for(sup=1;sup<MAXSUPSET;sup++)if(sets[arSet[l].setid].subsetid[sup]==coefs[index].setid[dcount]){varsupsetid[dcount]=sup;break;}
                  }
                  break;
                }
            }
    
    strcpy(line,linecopy);
    readitem = strtok(line,"=");
    readitem = strtok(NULL,";");
    i=str_count_char(readitem, ':');
    if(i!=7){printf("Syntax Error in %s!\n",readitem);return -1;}
    strcpy(line,readitem);
    strcat(line,":");
    strcpy(line1,line);
    ///x
    readitem = strtok(line,":");
    strcpy(vname,readitem);
        index=ncof-1;
        p=strtok(vname,"(");
        do {
          if (strcmp(coefs[index].cofname,p)==0) {
            xvar=index;
            if(coefs[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
            break;
          }
        } while (index--);
        for (l=0; l<MAXVARDIM; l++){xantidim[l]=0;xsubset[l]=0;xsupsetid[l]=0;}
            for (dcount=0; dcount<coefs[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].index_name,p)==0) {
                  xantidim[dcount]=coefs[index].strides[dcount];
                  xarset[dcount]=l+1;
                  if (sets[arSet[l].setid].size!=sets[coefs[index].setid[dcount]].size){
                    xsubset[dcount]=1;
                    for(sup=1;sup<MAXSUPSET;sup++)if(sets[arSet[l].setid].subsetid[sup]==coefs[index].setid[dcount]){xsupsetid[dcount]=sup;break;}
                  }
                  break;
                }
            }
    ///y
    strcpy(line,line1);
    readitem = strtok(line,":");
    readitem=strtok(NULL,":");
    strcpy(vname,readitem);
        index=ncof-1;
        p=strtok(vname,"(");
        do {
          if (strcmp(coefs[index].cofname,p)==0) {
            yvar=index;
            if(coefs[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
            break;
          }
        } while (index--);
        for (l=0; l<MAXVARDIM; l++){yantidim[l]=0;ysubset[l]=0;ysupsetid[l]=0;}
            for (dcount=0; dcount<coefs[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].index_name,p)==0) {
                  yantidim[dcount]=coefs[index].strides[dcount];
                  yarset[dcount]=l+1;
                  if (sets[arSet[l].setid].size!=sets[coefs[index].setid[dcount]].size){
                    ysubset[dcount]=1;
                    for(sup=1;sup<MAXSUPSET;sup++)if(sets[arSet[l].setid].subsetid[sup]==coefs[index].setid[dcount]){ysupsetid[dcount]=sup;break;}
                  }
                  break;
                }
            }
    
    ///j
    strcpy(line,line1);
    readitem = strtok(line,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    strcpy(vname,readitem);
        index=ncof-1;
        p=strtok(vname,"(");
        do {
          if (strcmp(coefs[index].cofname,p)==0) {
            jvar=index;
            if(coefs[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
            break;
          }
        } while (index--);
        for (l=0; l<MAXVARDIM; l++){jantidim[l]=0;jsubset[l]=0;jsupsetid[l]=0;}
            for (dcount=0; dcount<coefs[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].index_name,p)==0) {
                  jantidim[dcount]=coefs[index].strides[dcount];
                  jarset[dcount]=l+1;
                  if (sets[arSet[l].setid].size!=sets[coefs[index].setid[dcount]].size){
                    jsubset[dcount]=1;
                    for(sup=1;sup<MAXSUPSET;sup++)if(sets[arSet[l].setid].subsetid[sup]==coefs[index].setid[dcount]){jsupsetid[dcount]=sup;break;}
                  }
                  break;
                }
            }
    if(coefs[vvar].size!=coefs[jvar].size){
      printf("Syntax Error!!! %s and %s must have the same size\n",coefs[vvar].cofname,coefs[jvar].cofname);
      return 0;
    }
    for(l=0;l<coefs[vvar].size;l++){
      if(vararset[l]!=jarset[l]){
        printf("Syntax Error!!! %s and %s must have the same arguments\n",coefs[vvar].cofname,coefs[jvar].cofname);
        return 0;
      }
    }
    

    ///w
    strcpy(line,line1);
    readitem = strtok(line,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    strcpy(vname,readitem);
        index=ncof-1;
        p=strtok(vname,"(");
        do {
          if (strcmp(coefs[index].cofname,p)==0) {
            wvar=index;
            if(coefs[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
            break;
          }
        } while (index--);
        for (l=0; l<MAXVARDIM; l++){wantidim[l]=0;wsubset[l]=0;wsupsetid[l]=0;}
            for (dcount=0; dcount<coefs[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].index_name,p)==0) {
                  wantidim[dcount]=coefs[index].strides[dcount];
                  warset[dcount]=l+1;
                  if (sets[arSet[l].setid].size!=sets[coefs[index].setid[dcount]].size){
                    wsubset[dcount]=1;
                    for(sup=1;sup<MAXSUPSET;sup++)if(sets[arSet[l].setid].subsetid[sup]==coefs[index].setid[dcount]){wsupsetid[dcount]=sup;break;}
                  }
                  break;
                }
            }
    

    ///s
    strcpy(line,line1);
    readitem = strtok(line,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    strcpy(vname,readitem);
        index=nvar-1;
        p=strtok(vname,"(");
        do {
          if (strcmp(vars[index].cofname,p)==0) {
            svar=index;
            if(vars[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
            break;
          }
        } while (index--);
        sbegadd=ncofele+vars[svar].offset;
    printf("read %s\n",argu);
        for (l=0; l<MAXVARDIM; l++){santidim[l]=0;ssubset[l]=0;ssupsetid[l]=0;}
            for (dcount=0; dcount<vars[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].index_name,p)==0) {
                  santidim[dcount]=vars[index].strides[dcount];
                  sarset[dcount]=l+1;
                  if (sets[arSet[l].setid].size!=sets[vars[index].setid[dcount]].size){
                    ssubset[dcount]=1;
                    for(sup=1;sup<MAXSUPSET;sup++)if(sets[arSet[l].setid].subsetid[sup]==vars[index].setid[dcount]){ssupsetid[dcount]=sup;break;}
                  }
                  break;
                }
            }
    if(coefs[vvar].size!=coefs[svar].size){
      printf("Syntax Error!!! %s and %s must have the same size\n",coefs[vvar].cofname,coefs[svar].cofname);
      return 0;
    }
    for(l=0;l<coefs[vvar].size;l++){
      if(vararset[l]!=sarset[l]){
        printf("Syntax Error!!! %s and %s must have the same arguments\n",coefs[vvar].cofname,coefs[svar].cofname);
        return 0;
      }
    }

    ///vlmu
    strcpy(line,line1);
    readitem = strtok(line,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    strcpy(vname,readitem);
    printf("var1 %s\n",vname);
        index=ncof-1;
        p=strtok(vname,"(");
        do {
          if (strcmp(coefs[index].cofname,p)==0) {
            vlmu=index;
            if(coefs[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
            break;
          }
        } while (index--);
        for (l=0; l<MAXVARDIM; l++){mvantidim[l]=0;mvsubset[l]=0;mvsupsetid[l]=0;}
            for (dcount=0; dcount<coefs[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].index_name,p)==0) {
                  mvantidim[dcount]=coefs[index].strides[dcount];
                  mvarset[dcount]=l+1;
                  if (sets[arSet[l].setid].size!=sets[coefs[index].setid[dcount]].size){
                    mvsubset[dcount]=1;
                    for(sup=1;sup<MAXSUPSET;sup++)if(sets[arSet[l].setid].subsetid[sup]==coefs[index].setid[dcount]){mvsupsetid[dcount]=sup;break;}
                  }
                  break;
                }
            }
    
    ///dmlu
    strcpy(line,line1);
    readitem = strtok(line,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    strcpy(vname,readitem);
        index=ncof-1;
        p=strtok(vname,"(");
        do {
          if (strcmp(coefs[index].cofname,p)==0) {
            dlmu=index;
            if(coefs[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
            break;
          }
        } while (index--);
        for (l=0; l<MAXVARDIM; l++){mdantidim[l]=0;mdsubset[l]=0;mdsupsetid[l]=0;}
            for (dcount=0; dcount<coefs[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].index_name,p)==0) {
                  mdantidim[dcount]=coefs[index].strides[dcount];
                  mdarset[dcount]=l+1;
                  if (sets[arSet[l].setid].size!=sets[coefs[index].setid[dcount]].size){
                    mdsubset[dcount]=1;
                    for(sup=1;sup<MAXSUPSET;sup++)if(sets[arSet[l].setid].subsetid[sup]==coefs[index].setid[dcount]){mdsupsetid[dcount]=sup;break;}
                  }
                  break;
                }
            }
        for(l=0;l<4;l++)printf("var2 %lf dmlu %ldbeg %ld\n",elem_vals[coefs[dlmu].offset+l].value,dlmu,coefs[dlmu].offset);
        for(l=10;l<14;l++)printf("var2 %lf\n",elem_vals[coefs[dlmu].offset+l].value);

    ///wlmu
    strcpy(line,line1);
    readitem = strtok(line,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    readitem=strtok(NULL,":");
    strcpy(vname,readitem);
        index=ncof-1;
        p=strtok(vname,"(");
        do {
          if (strcmp(coefs[index].cofname,p)==0) {
            wlmu=index;
            if(vars[index].size>0){strcpy(argu,strtok(NULL,")"));strcat(argu,",");}
            break;
          }
        } while (index--);
        for (l=0; l<MAXVARDIM; l++){mwantidim[l]=0;mwsubset[l]=0;mwsupsetid[l]=0;}
            for (dcount=0; dcount<coefs[index].size; dcount++) {
              if(dcount==0)p=strtok(argu,",");
              else p=strtok(NULL,",");
              for (l=0; l<fdim-1; l++) if (strcmp(arSet[l].index_name,p)==0) {
                  mwantidim[dcount]=coefs[index].strides[dcount];
                  mwarset[dcount]=l+1;
                  if (sets[arSet[l].setid].size!=sets[coefs[index].setid[dcount]].size){
                    mwsubset[dcount]=1;
                    for(sup=1;sup<MAXSUPSET;sup++)if(sets[arSet[l].setid].subsetid[sup]==coefs[index].setid[dcount]){mwsupsetid[dcount]=sup;break;}
                  }
                  break;
                }
            }
    
    rsize=sets[arSet[fdim1-5].setid].size;
    tsize=sets[arSet[fdim1-4].setid].size;
    stosize=sets[arSet[fdim1-2].setid].size;
    bsize=sets[arSet[fdim1-1].setid].size;
    printf("rsize %ld tsize %ld stosize %ld bsize %ld\n",rsize,tsize,stosize,bsize);
    if(bsize!=4)printf("Errors!!! Splinter set with wrong size or in wrong position!!!\n");
    solve_real ***array = (solve_real ***) malloc(sizeof(solve_real **)*rsize);
    for (i = 0; i < rsize; i++)array[i] = (solve_real**) malloc(sizeof(solve_real *)*stosize);
    solve_real vecy2;
    solve_real *matuw= (solve_real *) calloc (rsize*stosize*stosize*bsize,sizeof(solve_real));
    solve_real *matuwold= (solve_real *) calloc (rsize*stosize*stosize*bsize,sizeof(solve_real));
    solve_real *matuv= (solve_real *) calloc (rsize*stosize*stosize,sizeof(solve_real));
    solve_real *matuk= (solve_real *) calloc (rsize*stosize,sizeof(solve_real));
    solve_real *curk= (solve_real *) calloc (rsize,sizeof(solve_real));
    offset_t *curpos= (offset_t *) calloc (rsize,sizeof(offset_t));
    solve_real *vecy= (solve_real *) calloc (rsize*tsize,sizeof(solve_real));
    solve_real *vecx= (solve_real *) calloc (rsize*tsize,sizeof(solve_real));
    solve_real *curx= (solve_real *) calloc (stosize,sizeof(solve_real));
    solve_real *cury= (solve_real *) calloc (stosize,sizeof(solve_real));
    solve_real *curw= (solve_real *) calloc (stosize*bsize,sizeof(solve_real));
    solve_real *vecwf= (solve_real *) calloc (coefs[wvar].nelem,sizeof(solve_real));
    solve_real *vecw=NULL,*vecwold=NULL;
    solve_real *vecwoldf= (solve_real *) calloc (coefs[wvar].nelem,sizeof(solve_real));
    vecw=vecwf;
    vecwold=vecwoldf;
    if(IsChange==false){
          if((!IsIni)||(IsSplint==2)){
                fout = fopen(filename, "rb");
                if (fout==NULL)printf("Weight file opening error\n");
                freadresult=fread(matuwold,sizeof(solve_real),rsize*stosize*stosize,fout);
                fclose(fout);
          }
          for (l=0; l<rsize; l++){
            arSet[0].indx=l;
            for (i=0; i<stosize; i++){
              arSet[2].indx=i;
              arSet[3].indx=i;//special for uk, only needs 2 set of k
              //dmlu
              l2=0;
              for (dcount=0; dcount<fdim-1; dcount++){
                for(i1=0; i1<coefs[dlmu].size; i1++) {
                  if(mvarset[i1]-1==dcount) {
                    if(mvsubset[i1]==1) {
                      l2=l2+set_elems[sets[coefs[dlmu].setid[i1]].offset+arSet[dcount].indx].superset_pos[mdsupsetid[i1]]*mdantidim[i1];
                    } else {
                      l2=l2+arSet[dcount].indx*mdantidim[i1];
                    }
                    break;
                  }
                }
              }
              matuk[l*stosize+i]=elem_vals[coefs[dlmu].offset+l2].value;
              printf("k1 %lf v %lf l2 %ld\n",elem_vals[coefs[dlmu].offset+l2].value,elem_vals[coefs[vlmu].offset+l2].value,l2);
              for (j=0; j<stosize; j++){
                arSet[3].indx=j;
                //vlmu
              l2=0;
              for (dcount=0; dcount<fdim-1; dcount++){
                for(i1=0; i1<coefs[vlmu].size; i1++) {
                  if(mvarset[i1]-1==dcount) {
                    if(mvsubset[i1]==1) {
                      l2=l2+set_elems[sets[coefs[vlmu].setid[i1]].offset+arSet[dcount].indx].superset_pos[mvsupsetid[i1]]*mvantidim[i1];
                    } else {
                      l2=l2+arSet[dcount].indx*mvantidim[i1];
                    }
                    break;
                  }
                }
              }
              matuv[l*stosize*stosize+i*stosize+j]=elem_vals[coefs[vlmu].offset+l2].value;
              }
            }
            for(i4=0;i4<tsize;i4++){
              arSet[1].indx=i4;
              //x
              l2=0;
              for (dcount=0; dcount<fdim-1; dcount++){
                for(i1=0; i1<coefs[xvar].size; i1++) {
                  if(xarset[i1]-1==dcount) {
                    if(xsubset[i1]==1) {
                      l2=l2+set_elems[sets[coefs[xvar].setid[i1]].offset+arSet[dcount].indx].superset_pos[xsupsetid[i1]]*xantidim[i1];
                    } else {
                      l2=l2+arSet[dcount].indx*xantidim[i1];
                    }
                    break;
                  }
                }
              }
              vecx[l*tsize+i4]=elem_vals[coefs[xvar].offset+l2].value;
              //y
              l2=0;
              for (dcount=0; dcount<fdim-2; dcount++){
                for(i1=0; i1<coefs[yvar].size; i1++) {
                  if(yarset[i1]-1==dcount) {
                    if(ysubset[i1]==1) {
                      l2=l2+set_elems[sets[coefs[yvar].setid[i1]].offset+arSet[dcount].indx].superset_pos[ysupsetid[i1]]*yantidim[i1];
                    } else {
                      l2=l2+arSet[dcount].indx*yantidim[i1];
                    }
                    break;
                  }
                }
              }
              vecy[l*tsize+i4]=elem_vals[coefs[yvar].offset+l2].value;
              printf("var %s l %ld l2 %ld vecy %lf vecx %lf\n",coefs[yvar].cofname,l,l2,vecy[l*tsize+i4],vecx[l*tsize+i4]);
            }
          }
          for (l=0; l<fdim+1; l++)arSet[l].indx=0;
          for (l=0; l<rsize; l++){
            arSet[0].indx=l;
            if(IsSplint==1||IsSplint==2){
              for(i4=0;i4<stosize;i4++){
                l2=l*stosize*stosize+i4*stosize;
                for(i1=0;i1<stosize;i1++)printf("k %lf v %lf\n",matuk[stosize+i1],matuv[l2+i1]);
                cubic_spline(matuv+l2,matuk+stosize,0,0,stosize-1,matuw+l*stosize*stosize*bsize+i4*stosize*bsize,laA);
              }
                memcpy(matuwold,matuw,rsize*stosize*stosize*bsize*sizeof(solve_real));
            } else memcpy(matuw,matuwold,rsize*stosize*stosize*bsize*sizeof(solve_real));
            }
          for (l=0; l<fdim+1; l++)arSet[l].indx=0;
          for (l=0; l<rsize; l++){
            arSet[0].indx=l;
            for(i4=0;i4<tsize;i4++){
              arSet[1].indx=i4;
              l2=0;
              for (dcount=0; dcount<fdim-2; dcount++){
                for(i1=0; i1<coefs[jvar].size; i1++) {
                  if(jarset[i1]-1==dcount) {
                    if(jsubset[i1]==1) {
                      l2=l2+set_elems[sets[coefs[jvar].setid[i1]].offset+arSet[dcount].indx].superset_pos[jsupsetid[i1]]*jantidim[i1];
                    } else {
                      l2=l2+arSet[dcount].indx*jantidim[i1];
                    }
                    break;
                  }
                }
              }
              temp3=(1+elem_vals[coefs[jvar].offset+l2].value/100)*vecy[l*tsize+i4];
              printf("temp3 %lf\n",temp3);
              for (i1=0; i1<rsize; i1++)curk[i1]=vecy[i1*tsize+i4];
              curk[l]=temp3;
              for (i1=0; i1<rsize; i1++)printf("curk0 %lf curk1 %lf\n",curk[0],vecy[l*tsize+i4]);
              //New x position
              for (i1=0; i1<rsize; i1++){
              if(matuk[i1*stosize]<matuk[i1*stosize+stosize-1]){
              if(curk[i1]<=matuk[i1*stosize])curpos[i1]=0;//force beg point
              else{
              for(j=1;j<stosize;j++){
                if(j!=stosize-1){
                if(curk[i1]<matuk[i1*stosize+j]&&curk[i1]>=matuk[i1*stosize+j-1]){
                  curpos[i1]=j-1;
                  break;
                }
                }else{
                  curpos[i1]=j-1;
                }
              }
              }
              }else{
              if(curk[i1]>=matuk[i1*stosize])curpos[i1]=0;//force beg point
              else{
              for(j=1;j<stosize;j++){
                if(j!=stosize-1){
                if(curk[i1]>matuk[i1*stosize+j]&&curk[i1]<=matuk[i1*stosize+j-1]){
                  curpos[i1]=j-1;
                  break;
                }
                }else{
                  curpos[i1]=j-1;
                }
              }
              }
              }
              printf("i1 %ld curk %lf cupos %ld\n",i1,curk[i1],curpos[i1]);
              }
              temp3=curk[rsize-1];
              for(j=0;j<stosize;j++){
                l2=l*stosize*stosize*bsize+j*stosize*bsize+curpos[rsize-1];
                curx[j]=matuwold[l2]+matuwold[l2+stosize-1]*temp3+matuwold[l2+2*(stosize-1)]*temp3*temp3+matuwold[l2+3*(stosize-1)]*temp3*temp3*temp3;
              }
              for(j=0;j<stosize;j++)printf("l %ld curx %lf curk %lfcurpos %ld curk0 %lf\n",l,curx[j],matuk[j],curpos[0],curk[0]);
              cubic_spline(curx,matuk,0,0,stosize-1,curw,laA);
              temp1=curw[curpos[0]]+curw[curpos[0]+stosize-1]*curk[0]+curw[curpos[0]+2*(stosize-1)]*curk[0]*curk[0]+curw[curpos[0]+3*(stosize-1)]*curk[0]*curk[0]*curk[0];

              l2=0;
              for (dcount=0; dcount<fdim-1; dcount++){
                for(i1=0; i1<coefs[vvar].size; i1++) {
                  if(vararset[i1]-1==dcount) {
                    if(varsubset[i1]==1) {
                      l2=l2+set_elems[sets[coefs[vvar].setid[i1]].offset+arSet[dcount].indx].superset_pos[varsupsetid[i1]]*varantidim[i1];
                    } else {
                      l2=l2+arSet[dcount].indx*varantidim[i1];
                    }
                    break;
                  }
                }
              }
                address=coefs[vvar].offset+l2;
                temp2=temp1/elem_vals[address].value;
                temp2=temp2*100-100;
                closure_vals[vars[svar].offset+l2].shock_value=temp2/subints;//shock, linear only, mup must be exo
                printf("l2 %ld shock %lf var %s temp1 %lf vval %lf shock %lf vvar %lf\n",l2,closure_vals[vars[svar].offset+l2].shock_value,vars[svar].cofname,temp1,elem_vals[address].value,closure_vals[vars[svar].offset+l2].shock_value,elem_vals[coefs[vvar].offset+l2].value);
              
            }
          }
    }
    printf("size1 %ld max error in percent %lf\n",tsize,maxerr);
    
    for (i = 0; i < rsize; i++) free(array[i]);
    free(matuk);
    free(matuv);
    free(matuw);
    free(matuwold);
    free(curk);
    free(curpos);
    free(vecy);
    free(vecx);
    free(curx);
    free(cury);
    free(curw);
    free(vecwf);
    free(vecwoldf);
    free(arSet);
  }
  fclose(filehandle);
  return j;
}

