#include <teems_solver.h>
#include <hsl_kernels.h>

int jacobian_fill(char *fname, char *commsyntax,set_def *sets,offset_t nset, set_element *set_elems, array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar, elem_value *elem_vals,offset_t ncofvar,offset_t ncofele,closure_entry *closure_vals,offset_t ndblock,offset_t alltimeset,offset_t allregset,offset_t *eq_addr,offset_t *counteq,offset_t nintraeq,Mat A,Mat B) {
  FILE * filehandle;
  char tline[TABREADLINE],line[TABREADLINE],line1[TABREADLINE],leftline[TABREADLINE],linecopy[TABREADLINE];//,set1[NAMESIZE],set2[NAMESIZE];
  char vname[TABREADLINE],sumsyntax[NAMESIZE],lintmp[TABREADLINE];//,*p1=NULL;
  char *readitem=NULL,*p=NULL,*p1=NULL;//,*p2=NULL,*varpnts;
  solve_real vval;
  PetscScalar    vvalpetsc;
  PetscInt Iindx,Jindx=0,Istart1,Iend1,matrow,rowindx;//,jfix=0
  PetscErrorCode ierr;
  solve_real zerodivide=0;
  PetscMPIInt  mpisize1;
  bool isinproc;
  ierr = MatGetOwnershipRange(A,&Istart1,&Iend1);
  MPI_Comm_size(PETSC_COMM_WORLD,&mpisize1);
  CHKERRQ(ierr);
  dim_t fdim,np,dcount,fdimlin=0,i4,sup,supset[MAXSUPSET];
  int totalsum,sumcount=1,sumcount1=0,lvar,lvar1,lvar2,lvar3,lvar4;
  offset_t lj,l1,i1=0,sumbegadd,dcountdim1[4*MAXVARDIM],dcountdim2[4*MAXVARDIM],dcountdim3[4*MAXVARDIM],nloops,nloopslin,nloopsfac,li3,nsumele,nsumele1,l2,eqindx=0,ltime;//,sizelinvars,totlinvars,templinvars
  int sumindx,npow,npar,nmul,nplu,ndiv,nmin,nops=0,nlinvars,leadlag,varindx1,varindx2;
  offset_t j,sj,l,i3,i5,i,arsetdim=0,nops_alloc=0;
  quantifier *arSet1=NULL;
  formula_op *ops1= NULL;

  filehandle = fopen(fname,"r");
  matrow=0;//nintraeq;

  while (tab_next_statement_resolved(commsyntax,filehandle,line,elem_vals,coefs,ncof,&zerodivide,TABREADLINE)) {
    if (strstr(line,"(default")==NULL) {
      str_replace_first(line, commsyntax, "");
      str_replace_first(line, "(linear)", "");
      while (str_replace_all(line,"  ", " "));
      while (str_replace_char(line, '[', '('));
      while (str_replace_char(line, ']', ')'));
      while (str_replace_char(line, '{', '('));
      while (str_replace_char(line, '}', ')'));
      strcpy(linecopy,line);
      fdim=str_count_ci(line, "(all,");
      if (fdim==0) {
        readitem = strtok(line+1," ");
        readitem = strtok(NULL,"=");
        strcpy(vname,readitem);
        strcpy(line,linecopy);
        readitem = strtok(line,"=");
        readitem = strtok(NULL,";");
        strcat(readitem,"-");
        strcat(readitem,"(");
        strcat(readitem,vname);
        strcat(readitem,")");
      }
      else {
        i=str_rfind_ci(line, "(all,");
        readitem=line+i;
        readitem = strtok(readitem,")");
        readitem = strtok(NULL,"=");
        strcpy(vname,readitem);
        strcpy(line,linecopy);
        readitem = strtok(line,"=");
        readitem = strtok(NULL,";");
        strcat(readitem,"-");
        strcat(readitem,"(");
        strcat(readitem,vname);
        strcat(readitem,")");
      }
      while (str_replace_all(readitem," ", ""));
      while (formula_normalize(readitem)==1);
      leadlag_encode(readitem);
      strcpy(tline,readitem);
      strcpy(line1,readitem);
      npow=str_count_char(readitem, '^');
      nmul=str_count_char(readitem, '*');
      ndiv=str_count_char(readitem, '/');
      nmul=nmul+ndiv;
      nplu=str_count_char(readitem, '+');
      nmin=str_count_char(readitem, '-');
      nplu=nplu+nmin;
      npar=str_count_char(readitem, '(');

      strcpy(line,line1);
      readitem=line;
      np=str_count_ci(readitem,"p_");
      nops_alloc=2*(npow+nmul+nplu+npar+1);
      formula_op *ops= (formula_op *) calloc (nops_alloc,sizeof(formula_op));
      eq_var_ref *LinVars= (eq_var_ref *) calloc (np,sizeof(eq_var_ref));
      i3=0;
      lvar=0;
      for (i=0; i<np; i++) {
        varindx2=0;
        while(-1<0) {
          varindx1=str_find_ci(readitem+varindx2,"p_");
          if(varindx1==-1) break;
          varindx2=varindx2+varindx1;
          if(varindx2==0||readitem[varindx2-1]=='*'||readitem[varindx2-1]=='+'||readitem[varindx2-1]=='-'||readitem[varindx2-1]=='('||readitem[varindx2-1]==',') break;
          else varindx2++;
        }
        if(varindx1==-1) break;
        lvar=lvar+varindx2;
        readitem=readitem+varindx2+2;

        p=strpbrk(readitem,"{+*-/^)");
        if(p!=NULL) {
          strncpy(vname,readitem,p-readitem);
          vname[p-readitem]='\0';
          strcpy(LinVars[i3].LinVarName,vname);
          for (l=0; l<nvar; l++) {
            if (strcmp(vars[l].cofname,vname)==0) {
              LinVars[i3].LinVarIndx=l;
              break;
            }
          }
          p=strpbrk(readitem,"}+*-/^)");
          if (*p=='}') {
            p++;
          }
          strncpy(tline,readitem,p-readitem);
          tline[p-readitem]='\0';
          switch (vars[l].size) {
          case 0:
            break;
          case 1:
            p = strtok(tline,"{");
            p = strtok(NULL,"}");
            leadlag=0;
            parse_index_leadlag(p,&leadlag);
            strcpy(LinVars[i3].dimnames[0],p);
            LinVars[i3].dimleadlag[0]=leadlag;
            strcpy(lintmp,"(all,");
            strcat(lintmp,p);
            strcat(lintmp,",");
            l1=str_find_ci(linecopy,lintmp);
            if (l1>-1) {
              p1=&linecopy[0]+l1;
              strncpy(LinVars[i3].dimsetnames[0],p1+strlen(lintmp),strchr(p1,')')-p1-strlen(lintmp));
            }
            else {
              strcpy(lintmp,"sum(");
              strcat(lintmp,p);
              strcat(lintmp,",");
              lvar1=str_count_ci(linecopy,lintmp);
              lvar3=str_find_ci(linecopy,lintmp);
              if (lvar1>1) for(lvar2=0; lvar2<lvar1; lvar2++) {
                  lvar4=str_find_ci(&linecopy[lvar3+4],lintmp);
                  if (lvar4>-1&&lvar4<lvar) {
                    lvar3=lvar3+lvar4+4;
                  }
                  else {
                    break;
                  }
                }
              p1=&linecopy[0]+lvar3;
              p1=p1+strlen(lintmp);
              strncpy(LinVars[i3].dimsetnames[0],p1,strchr(p1,',')-p1);
            }
            break;
          default:
            p = strtok(tline,"{");
            for (i4=0; i4<vars[l].size-1; i4++) {
              p = strtok(NULL,",");
              leadlag=0;
              parse_index_leadlag(p,&leadlag);
              strcpy(LinVars[i3].dimnames[i4],p);
              LinVars[i3].dimleadlag[i4]=leadlag;
              strcpy(lintmp,"(all,");
              strcat(lintmp,p);
              strcat(lintmp,",");
              l1=str_find_ci(linecopy,lintmp);
              if (l1>-1) {
                p1=&linecopy[0]+l1;
                strncpy(LinVars[i3].dimsetnames[i4],p1+strlen(lintmp),strchr(p1,')')-p1-strlen(lintmp));
              }
              else {
                strcpy(lintmp,"sum(");
                strcat(lintmp,p);
                strcat(lintmp,",");
                lvar1=str_count_ci(linecopy,lintmp);
                lvar3=str_find_ci(linecopy,lintmp);
                if (lvar1>1) for(lvar2=0; lvar2<lvar1; lvar2++) {
                    lvar4=str_find_ci(&linecopy[lvar3+4],lintmp);
                    if (lvar4>-1&&lvar4<lvar) {
                      lvar3=lvar3+lvar4+4;
                    }
                    else {
                      break;
                    }
                  }
                p1=&linecopy[0]+lvar3;
                p1=p1+strlen(lintmp);
                strncpy(LinVars[i3].dimsetnames[i4],p1,strchr(p1,',')-p1);
              }
            }
            p = strtok(NULL,"}");
            leadlag=0;
            parse_index_leadlag(p,&leadlag);
            strcpy(LinVars[i3].dimnames[i4],p);
            LinVars[i3].dimleadlag[i4]=leadlag;
            strcpy(lintmp,"(all,");
            strcat(lintmp,p);
            strcat(lintmp,",");
            l1=str_find_ci(linecopy,lintmp);
            if (l1>-1) {
              p1=&linecopy[0]+l1;
              strncpy(LinVars[i3].dimsetnames[i4],p1+strlen(lintmp),strchr(p1,')')-p1-strlen(lintmp));
            }
            else {
              strcpy(lintmp,"sum(");
              strcat(lintmp,p);
              strcat(lintmp,",");
              lvar1=str_count_ci(linecopy,lintmp);
              lvar3=str_find_ci(linecopy,lintmp);
              if (lvar1>1) for(lvar2=0; lvar2<lvar1; lvar2++) {
                  lvar4=str_find_ci(&linecopy[lvar3+4],lintmp);
                  if (lvar4>-1&&lvar4<lvar) {
                    lvar3=lvar3+lvar4+4;
                  }
                  else {
                    break;
                  }
                }
              p1=&linecopy[0]+lvar3;
              p1=p1+strlen(lintmp);
              strncpy(LinVars[i3].dimsetnames[i4],p1,strchr(p1,',')-p1);
            }
            break;
          }
          i3++;
        }
        else {
          strcpy(vname,readitem);
          strcpy(LinVars[i3].LinVarName,vname);
          for (l=0; l<nvar; l++) {
            if (strcmp(vars[l].cofname,vname)==0) {
              LinVars[i3].LinVarIndx=l;
              break;
            }
          }
          i3++;
        }
      }
      nlinvars=i3;
      i3=0;
      for (i=0; i<nlinvars; i++) i3+=vars[LinVars[i].LinVarIndx].size;
      arsetdim=fdim+i3;
      quantifier *arSet= (quantifier *) calloc (arsetdim,sizeof(quantifier));

      strcpy(line,linecopy);
      nloops=1;
      if (fdim>0) {
        for (i=0; i<fdim; i++) {
          if(i==0) {
            readitem = strtok(line,",");
          }
          else {
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
        dcountdim1[fdim-1]=1;
        for (dcount=fdim-2; dcount>-1; dcount--) {
          dcountdim1[dcount]=sets[arSet[dcount+1].setid].size*dcountdim1[dcount+1];
        }
      }
      isinproc=false;
      if(mpisize1>1) {
        for (lj=0; lj<nloops; lj++) {
          Jindx=eq_addr[matrow+lj];
          if(Jindx>=Istart1&&Jindx<Iend1) {
            isinproc=true;
            break;
          }
        }
      }
      else isinproc=true;
      if(isinproc) {
        strcpy(line,line1);
        readitem=line;
        strcpy(sumsyntax,"sum(");
        totalsum=sum_count(readitem,sumsyntax);
        sum_def *sum_cof= (sum_def *) calloc (totalsum*nlinvars+1,sizeof(sum_def));
        sumcount=0;
        while (sum_parse(readitem,sumsyntax,sum_cof,arSet,sets,nset,fdim+1,sumcount)==1) {
          sumcount++;
        }
        totalsum=sumcount;
        li3=0;
        for (i4=0; i4<totalsum; i4++) {
          i1=1;
          for(sj=0; sj<sum_cof[i4].size; sj++) {
            i1=i1*sets[sum_cof[i4].setid[sj]].size;
          }
          sum_cof[i4].offset=li3;
          li3=li3+i1;
        }
        nsumele=li3;
        sum_value *sum_vals= (sum_value *) calloc (nsumele*nlinvars,sizeof(sum_value));

        strcpy(line,line1);
        readitem=line;
        sumcount=0;
        while (eq_sum_parse(readitem,sumsyntax,sum_cof,arSet,sets,nset,fdim+1,sumcount)==1) {
          sumcount++;
        }
        totalsum=sumcount;
        li3=0;
        for (i4=0; i4<totalsum; i4++) {
          i1=1;
          for(sj=0; sj<sum_cof[i4].size; sj++) {
            i1=i1*sets[sum_cof[i4].setid[sj]].size;
          }
          sum_cof[i4].offset=li3;
          li3=li3+i1;
        }
        nsumele=li3;
        for (i3=0; i3<totalsum; i3++) {
          i1=1;
          sum_cof[i3].strides[sum_cof[i3].size-1]=1;
          for(sj=sum_cof[i3].size-2; sj>-1; sj--) {
            sum_cof[i3].strides[sj]=sum_cof[i3].strides[sj+1]*sets[sum_cof[i3].setid[sj+1]].size;
          }
        }
        sumcount=0;
        strcpy(line,line1);
        readitem=line;
        sumindx=0;
        while (eq_sum_eval(readitem,sumsyntax,sets,nset,set_elems,elem_vals,ncofvar,ncofele,coefs,ncof,vars,nvar,sum_cof,totalsum,sum_vals,nsumele,ops,arSet,fdim+1,&sumindx,sumcount,zerodivide)==1) {
          sumcount++;
        }
        strcpy(line1,readitem);
        sumbegadd=nsumele;
        sumcount1=sumcount;
        for (i=0; i<nlinvars; i++) {
          Jindx=eq_addr[matrow];
          if(Jindx>=Iend1)continue;
          i3=0;
          nloopslin=nloops;
          if (fdim==0) {
            for (l2=0; l2<vars[LinVars[i].LinVarIndx].size; l2++) {
              strcpy(arSet[i3].index_name,LinVars[i].dimnames[l2]);
              for (i4=0; i4<nset; i4++) if(strcmp(sets[i4].setname,LinVars[i].dimsetnames[l2])==0) {
                  break;
                }
              arSet[i3].setid=i4;
              nloopslin=nloopslin*sets[arSet[i3].setid].size;
              i3++;
            }
          }
          else for (l2=0; l2<vars[LinVars[i].LinVarIndx].size; l2++) {
              for (sj=0; sj<fdim; sj++) {
                if(strcmp(arSet[sj].index_name,LinVars[i].dimnames[l2])==0) {
                  break;
                }
                else if (sj==fdim-1) {
                  strcpy(arSet[fdim+i3].index_name,LinVars[i].dimnames[l2]);
                  for (i4=0; i4<nset; i4++) if(strcmp(sets[i4].setname,LinVars[i].dimsetnames[l2])==0) {
                      break;
                    }
                  arSet[fdim+i3].setid=i4;
                  nloopslin=nloopslin*sets[arSet[fdim+i3].setid].size;
                  i3++;
                }
              }
            }
          nloopsfac=(offset_t)nloopslin/nloops;
          Jindx=eq_addr[matrow+(offset_t)(nloopslin-1)/nloopsfac];
          if(Jindx<Istart1)continue;
          
          fdimlin=fdim+i3;
          if (i3>0) {
            dcountdim2[fdimlin-1]=1;
            for (dcount=fdimlin-2; dcount>-1; dcount--) {
              dcountdim2[dcount]=sets[arSet[dcount+1].setid].size*dcountdim2[dcount+1];
            }
          }
          else for (dcount=0; dcount<fdim; dcount++) {
              dcountdim2[dcount]=dcountdim1[dcount];
            }

          strcpy(leftline,line1);
          eq_zero_linvar(leftline,i);
          eq_replace_linvar(leftline,0);
          strcpy(sumsyntax,"sum(");
          eq_sum_replace(leftline,sumsyntax,i,LinVars,vars);
          eq_linvar_read(leftline,LinVars,i,vars);
          strcpy(line,leftline);
          readitem=line;
          sumcount1=sumcount;
          while (sum_parse(readitem,sumsyntax,sum_cof,arSet,sets,nset,fdimlin+1,sumcount1)==1) {
            sumcount1++;
          }
          totalsum=sumcount1;
          for (i4=sumcount; i4<totalsum; i4++) {
            i1=1;
            for(sj=0; sj<sum_cof[i4].size; sj++) {
              i1=i1*sets[sum_cof[i4].setid[sj]].size;
            }
            sum_cof[i4].offset=sumbegadd;
            sumbegadd=sumbegadd+i1;
          }
          nsumele1=sumbegadd+nsumele;
          for (i3=sumcount; i3<totalsum; i3++) {
            i1=1;
            sum_cof[i3].strides[sum_cof[i3].size-1]=1;
            for(sj=sum_cof[i3].size-2; sj>-1; sj--) {
              sum_cof[i3].strides[sj]=sum_cof[i3].strides[sj+1]*sets[sum_cof[i3].setid[sj+1]].size;
            }
          }
          strcpy(line,leftline);
          readitem=line;
          while (sum_eval(readitem,sumsyntax,sets,nset,set_elems,elem_vals,ncofvar,ncofele,coefs,ncof,vars,nvar,sum_cof,totalsum,sum_vals,nsumele1,ops,arSet,fdimlin+1,&sumindx,sumcount,zerodivide)==1) {
            sumcount++;
          }
          nops=0;
          formula_compile(readitem,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,&nops,arSet,fdimlin);
          for (dcount=0; dcount<vars[LinVars[i].LinVarIndx].size; dcount++) {
            for (i4=0; i4<fdimlin; i4++) {
              if (strcmp(LinVars[i].dimnames[dcount],arSet[i4].index_name)==0) {
                dcountdim3[dcount]=i4;
                break;
              }
            }
          }
          for(dcount=0; dcount<MAXSUPSET; dcount++)supset[dcount]=0;
          for (dcount=0; dcount<vars[LinVars[i].LinVarIndx].size; dcount++) {
            if(sets[vars[LinVars[i].LinVarIndx].setid[dcount]].size!=sets[arSet[dcountdim3[dcount]].setid].size) {
              for(sup=1; sup<MAXSUPSET; sup++)if(vars[LinVars[i].LinVarIndx].setid[dcount]==sets[arSet[dcountdim3[dcount]].setid].subsetid[sup]) {
                  supset[dcount]=sup;
                  break;
                }
            }
            else supset[dcount]=0;
          }
        #pragma omp parallel private(lj,Jindx,i3,sj,i5,l2,dcount,l1,li3,Iindx,ierr,arSet1,ops1,vval) shared(elem_vals,arSet)
        {
        if(omp_get_thread_num()!=0){
          arSet1=realloc(arSet1,arsetdim*sizeof(quantifier));
          memcpy(arSet1,arSet,arsetdim*sizeof(quantifier));
          ops1=realloc(ops1,nops_alloc*sizeof(formula_op));
          memcpy(ops1,ops,nops_alloc*sizeof(formula_op));
        }else{
          ops1=ops;
          arSet1=arSet;
        }
          solve_real *value= (solve_real *) calloc (nloopsfac,sizeof(solve_real));
          PetscInt *jcn= (PetscInt *) calloc (nloopsfac,sizeof(PetscInt));
          solve_real *valueb= (solve_real *) calloc (nloopsfac,sizeof(solve_real));
          PetscInt *jcnb= (PetscInt *) calloc (nloopsfac,sizeof(PetscInt));
        #pragma omp for
          for (i5=0; i5<nloops; i5++) {
            Jindx=eq_addr[matrow+i5];//Jindx=ha_eqadd[matrow+(uvadd)lj/nloopsfac];
            if(Jindx>=Istart1&&Jindx<Iend1) {
              i3=0;
              sj=0;
              for (lj=i5*nloopsfac; lj<(i5+1)*nloopsfac; lj++) {
              l2=lj;
              for (dcount=0; dcount<fdimlin; dcount++) {
                l1=(offset_t) l2/dcountdim2[dcount];
                arSet1[dcount].indx=l1;
                l2=l2-l1*dcountdim2[dcount];
              }
              li3=0;
              for (dcount=0; dcount<vars[LinVars[i].LinVarIndx].size; dcount++) {
                if(supset[dcount]==0) {
                  li3=li3+(arSet1[dcountdim3[dcount]].indx+LinVars[i].dimleadlag[dcount])*vars[LinVars[i].LinVarIndx].strides[dcount];
                }
                else {
                  li3=li3+(set_elems[sets[arSet1[dcountdim3[dcount]].setid].offset+arSet1[dcountdim3[dcount]].indx].superset_pos[supset[dcount]]+LinVars[i].dimleadlag[dcount])*vars[LinVars[i].LinVarIndx].strides[dcount];
                }
              }
              vval=formula_eval(elem_vals,sets,set_elems,sum_vals,ops1,nops,arSet1,fdimlin,zerodivide);
              Iindx=closure_vals[vars[LinVars[i].LinVarIndx].offset+li3].exo_index;
              if (!closure_vals[vars[LinVars[i].LinVarIndx].offset+li3].is_exogenous&&vval!=0) {
                value[i3]=vval;
                jcn[i3]=Iindx;
                i3++;
              }
              if (closure_vals[vars[LinVars[i].LinVarIndx].offset+li3].is_exogenous&&vval!=0) {
                valueb[sj]=-vval;
                jcnb[sj]=Iindx;
                sj++;
              }
            }
            if(i3>0)MatSetValues(A,1,&Jindx,i3,jcn,value,ADD_VALUES);
            if(sj>0)MatSetValues(B,1,&Jindx,sj,jcnb,valueb,ADD_VALUES);
            }
          }
          free(value);
          free(jcn);
          free(valueb);
          free(jcnb);
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
        free(sum_cof);
        free(sum_vals);
      }
      matrow+=nloops;
      eqindx++;

      free(ops);
      free(LinVars);
      free(arSet);
    }
  }
  fclose(filehandle);
  return 1;
}

int eq_sum_parse(char *formulain, char *commsyntax, sum_def *sum_cof,quantifier *arSet,set_def *sets,dim_t nset,dim_t fdim,int j) {
  char *readitem,*p,*p1,*p2,interchar2[NAMESIZE],argu[TABREADLINE],tempname[NAMESIZE];//,line5[TABREADLINE]
  char interchar[NAMESIZE],interchar1[NAMESIZE],line[TABREADLINE],line1[TABREADLINE],line2[TABREADLINE],line3[TABREADLINE],line4[TABREADLINE];
  dim_t l,l1,l2,l3,l4,l5,l6,l7;
  int ncur=0,ncuri,i=0,k=0,k1=0,length;
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
      }
      else {
        if(str_find_ci(line,",p_")>-1||str_find_ci(line,"*p_")>-1||str_find_ci(line,"+p_")>-1||str_find_ci(line,"-p_")>-1||str_find_ci(line,"(p_")>-1) {
          i=i+k+4;
          readitem=formulain+i;
        }
        else {
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
          for (l7=0; l7<nset; l7++) if(strcmp(p,sets[l7].setname)==0) {
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
                      sum_cof[j].setid[l3]=l5;
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
                    for (l7=0; l7<nset; l7++) if(strcmp(p1,sets[l7].setname)==0) {
                        sum_cof[j].setid[l3]=l7;
                        break;
                      }
                  }
                  l3++;
                }
              }
            }
            else {
              for (l1=0; l1<l; l1++) {
                if(l1==0) {
                  p = strtok(argu,",");
                }
                else {
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
                        sum_cof[j].setid[l3]=l5;
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
          if (interchar[strlen(interchar)-1]==',') {
            interchar[strlen(interchar)-1]='}';
          }
          else {
            if (interchar[strlen(interchar)-1]=='{') {
              interchar[strlen(interchar)-1]='\0';
            }
            else {
              strcat(interchar,"}");
            }
          }
          while(str_replace_all(formulain,line4,interchar));
          sum_cof[j].size=l3;
          return 1;
        }
      }

    }
    else if (formulain[i+k-1]=='+'||formulain[i+k-1]=='-'||formulain[i+k-1]=='*'||formulain[i+k-1]=='/'||formulain[i+k-1]=='^'||formulain[i+k-1]=='('||formulain[i+k-1]==',') {
      readitem=formulain+i+k;
      strcpy(line,readitem);
      strcpy(line1,readitem);
      sum_extract(line);
      k1=str_find_ci(line+4,commsyntax);
      if (k1!=-1) {
        i=i+k+4;
        readitem=formulain+i;
      }
      else {
        if(str_find_ci(line,",p_")>-1||str_find_ci(line,"*p_")>-1||str_find_ci(line,"+p_")>-1||str_find_ci(line,"-p_")>-1||str_find_ci(line,"(p_")>-1) {
          i=i+k+4;
          readitem=formulain+i;
        }
        else {
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
          for (l7=0; l7<nset; l7++) if(strcmp(p,sets[l7].setname)==0) {
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
                  strcat(interchar,sum_cof[j].dimnames[l3]);
                  l6=0;
                  for (l5=0; l5<fdim-1; l5++) if(strcmp(p,arSet[l5].index_name)==0) {
                      sum_cof[j].setid[l3]=l5;
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
                    for (l7=0; l7<nset; l7++) if(strcmp(p1,sets[l7].setname)==0) {
                        sum_cof[j].setid[l3]=l7;
                        break;
                      }
                  }
                  l3++;
                }
              }
            }
            else {
              for (l1=0; l1<l; l1++) {
                if(l1==0) {
                  p = strtok(argu,",");
                }
                else {
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
                        sum_cof[j].setid[l3]=arSet[l5].setid;//l5;
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
          if (interchar[strlen(interchar)-1]==',') {
            interchar[strlen(interchar)-1]='}';
          }
          else {
            if (interchar[strlen(interchar)-1]=='{') {
              interchar[strlen(interchar)-1]='\0';
            }
            else {
              strcat(interchar,"}");
            }
          }
          while(str_replace_all(formulain,line4,interchar));
          sum_cof[j].size=l3;
          return 1;
        }
      }

    }
    else {
      i=i+k+4;
      readitem=formulain+i;
    }
  }
  return 0;
}

int eq_sum_eval(char *formulain, char *commsyntax,set_def *sets,dim_t nset, set_element *set_elems,elem_value *elem_vals,offset_t ncofvar,offset_t ncofele, array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar,sum_def *sum_cof,int totalsum,sum_value *sum_vals,offset_t nsumele,formula_op *ops,quantifier *arSet1,dim_t fdim,int *sumindx,int j, solve_real zerodivide) {
  char *readitem,*p;//,*p1,interchar2[NAMESIZE],line5[TABREADLINE];
  char interchar[NAMESIZE],line[TABREADLINE],line1[TABREADLINE],line2[TABREADLINE];//,line3[TABREADLINE],line4[TABREADLINE];//,interchar1[NAMESIZE]
  int i=0,k=0,k1=0,length;//,simpl=0;//,ncur=0,ncuri,l3,l4,l5,l6,l7
  dim_t dcount,superset_pos,fdimsumcof,l;
  offset_t dcountdim1[4*MAXVARDIM],nloops,l1,l2,l3;
  int nops;
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
      }
      else {
        if(str_find_ci(line,",p_")>-1||str_find_ci(line,"*p_")>-1||str_find_ci(line,"+p_")>-1||str_find_ci(line,"-p_")>-1||str_find_ci(line,"(p_")>-1) {
          i=i+k+4;
          readitem=formulain+i;
        }
        else {
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
            }
            else {
              dcountdim1[dcount-1]=dcountdim1[dcount]*sets[arSet[dcount].setid].size;
            }
          }
          arSet[sum_cof[j].size].setid=sum_cof[j].sumsetid;
          strcpy(arSet[sum_cof[j].size].index_name,sum_cof[j].sumindx);
          fdimsumcof=sum_cof[j].size+1;
          nops=0;
          formula_compile(p,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,&nops,arSet,fdimsumcof);
        #pragma omp parallel private(l3,l1,l2,dcount,superset_pos,vval,arSet2,ops1) shared(elem_vals,arSet,sum_vals)
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
          for (l3=0; l3<nloops; l3++) {
            l2=l3;
            for (dcount=0; dcount<sum_cof[j].size; dcount++) {
              superset_pos=(offset_t) l2/dcountdim1[dcount];
              arSet2[dcount].indx=superset_pos;
              l2=l2-superset_pos*dcountdim1[dcount];
            }
            vval=0;
            for (l1=0; l1<sets[sum_cof[j].sumsetid].size; l1++) {
              arSet2[sum_cof[j].size].indx=l1;
              vval+=formula_eval(elem_vals,sets,set_elems,sum_vals,ops1,nops,arSet2,fdimsumcof,zerodivide);
            }
            sum_vals[*sumindx+l3].value=vval;
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
          }
          else {
            if (interchar[strlen(interchar)-1]=='{') {
              interchar[strlen(interchar)-1]='\0';
            }
            else {
              strcat(interchar,"}");
            }
          }
          while(str_replace_all(formulain,line1,interchar)!=NULL);
          free(arSet);
          return 1;
        }
      }

    }
    else if (formulain[i+k-1]=='+'||formulain[i+k-1]=='-'||formulain[i+k-1]=='*'||formulain[i+k-1]=='/'||formulain[i+k-1]=='^'||formulain[i+k-1]=='('||formulain[i+k-1]==',') {
      readitem=formulain+i+k;
      strcpy(line,readitem);
      sum_extract(line);
      k1=str_find_ci(line+4,commsyntax);
      if (k1!=-1) {
        i=i+k+4;
        readitem=formulain+i;
      }
      else {
        if(str_find_ci(line,",p_")>-1||str_find_ci(line,"*p_")>-1||str_find_ci(line,"+p_")>-1||str_find_ci(line,"-p_")>-1||str_find_ci(line,"(p_")>-1) {
          i=i+k+4;
          readitem=formulain+i;
        }
        else {
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
            }
            else {
              dcountdim1[dcount-1]=dcountdim1[dcount]*sets[arSet[dcount].setid].size;
            }
          }
          arSet[sum_cof[j].size].setid=sum_cof[j].sumsetid;
          strcpy(arSet[sum_cof[j].size].index_name,sum_cof[j].sumindx);
          fdimsumcof=sum_cof[j].size+1;
          nops=0;
          formula_compile(p,sets,coefs,ncof,vars,nvar,ncofele,sum_cof,totalsum,ops,&nops,arSet,fdimsumcof);
        #pragma omp parallel private(l3,l1,l2,dcount,superset_pos,vval,arSet2,ops1) shared(elem_vals,arSet,sum_vals)
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
          for (l3=0; l3<nloops; l3++) {
            l2=l3;
            for (dcount=0; dcount<sum_cof[j].size; dcount++) {
              superset_pos=(offset_t) l2/dcountdim1[dcount];
              arSet2[dcount].indx=superset_pos;
              l2=l2-superset_pos*dcountdim1[dcount];
            }
            vval=0;
            for (l1=0; l1<sets[sum_cof[j].sumsetid].size; l1++) {
              arSet2[sum_cof[j].size].indx=l1;
              vval+=formula_eval(elem_vals,sets,set_elems,sum_vals,ops1,nops,arSet2,fdimsumcof,zerodivide);
            }
            sum_vals[(offset_t)*sumindx+l3].value=vval;//ha_sumele[*sumindx+l2].varval=vval;
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
          }
          else {
            if (interchar[strlen(interchar)-1]=='{') {
              interchar[strlen(interchar)-1]='\0';
            }
            else {
              strcat(interchar,"}");
            }
          }
          while(str_replace_all(formulain,line1,interchar));
          free(arSet);
          return 1;
        }
      }
    }
    else {
      i=i+k+4;
      readitem=formulain+i;
    }
  }
  return 0;
}

int eq_sum_replace(char *formulain, char *commsyntax,int LinIndx, eq_var_ref *LinVars,array_def *vars) {
  char *readitem,line[TABREADLINE],line1[TABREADLINE],line2[TABREADLINE],line3[TABREADLINE],argu[TABREADLINE],*p0,*p1;
  int i=0,k=0,k1=0,k2,k3,length,tsum,j,l;
  bool IsRemSum;
  length=strlen(formulain);
  readitem=formulain;
  while (i<length) {
    k=str_find_ci(formulain+i,commsyntax);
    if (k==-1) {
      return 0;
    }
    if (k==0) {
      readitem=formulain+i+k;
      strcpy(line,readitem);
      sum_extract(line);
      strcpy(line3,line);
      if(str_find_ci(line,",p_")>-1||str_find_ci(line,"*p_")>-1||str_find_ci(line,"+p_")>-1||str_find_ci(line,"-p_")>-1||str_find_ci(line,"(p_")>-1) {
        tsum=sum_count(line,commsyntax);
        IsRemSum=false;
        for (j=tsum; j>1; j--) {
          k2=1;
          k1=0;
          while(k2<j) {
            k3=k1;
            k1=str_find_ci(line+k1+3,commsyntax);
            k1=k3+k1+3;//printf("k1 %d line %s\n",k1,line+k1-1);
            if (line[k1-1]=='+'||line[k1-1]=='-'||line[k1-1]=='*'||line[k1-1]=='/'||line[k1-1]=='^'||line[k1-1]=='('||line[k1-1]==',') k2++;//printf("k1 %d line %s\n",k1,line+k1-1);}
          }
          strcpy(line1,&line[k1]);
          sum_extract(line1);
          p0=strchr(line1,'(');
          p0++;
          p1=strchr(line1,',');
          strncpy(argu,p0,p1-p0);
          argu[p1-p0]='\0';
          for(l=0; l<vars[LinVars[LinIndx].LinVarIndx].size; l++) {
            if(strcmp(argu,LinVars[LinIndx].dimnames[l])==0) {
              p1++;
              p1=strchr(p1,',');
              p1++;
              strncpy(line2,p1,strlen(line1)-(p1-line1)-1);
              line2[strlen(line1)-(p1-line1)-1]='\0';
              str_replace_all(line,line1,line2);
              IsRemSum=true;
              break;
            }
          }
        }
        p0=strchr(line,'(');
        p0++;
        p1=strchr(line,',');
        strncpy(argu,p0,p1-p0);
        argu[p1-p0]='\0';
        for(l=0; l<vars[LinVars[LinIndx].LinVarIndx].size; l++) {
          if(strcmp(argu,LinVars[LinIndx].dimnames[l])==0) {
            p1++;
            p1=strchr(p1,',');
            p1++;
            line2[0]='(';
            strncpy(&line2[1],p1,strlen(line)-(p1-line)-1);
            line2[strlen(line)-(p1-line)]=')';
            line2[strlen(line)-(p1-line)+1]='\0';
            strcpy(line,line2);
            IsRemSum=true;
            break;
          }
        }
        if(IsRemSum) {
          str_replace_all(formulain,line3,line);
          i=i+k+strlen(line);
        }
        else i=i+k+strlen(line3);
      }
      else {
        str_replace_all(formulain,line,"0");
        i=i+k+1;
      }
    }
    else if (formulain[i+k-1]=='+'||formulain[i+k-1]=='-'||formulain[i+k-1]=='*'||formulain[i+k-1]=='/'||formulain[i+k-1]=='^'||formulain[i+k-1]=='('||formulain[i+k-1]==',') {
      readitem=formulain+i+k;
      strcpy(line,readitem);
      sum_extract(line);
      strcpy(line3,line);
      if(str_find_ci(line,",p_")>-1||str_find_ci(line,"*p_")>-1||str_find_ci(line,"+p_")>-1||str_find_ci(line,"-p_")>-1||str_find_ci(line,"(p_")>-1) {
        tsum=sum_count(line,commsyntax);
        IsRemSum=false;
        for (j=tsum; j>1; j--) {
          k2=1;
          k1=0;
          while(k2<j) {
            k3=k1;
            k1=str_find_ci(line+k1+3,commsyntax);
            k1=k3+k1+3;//printf("k1 %d line %s\n",k1,line+k1-1);
            if (line[k1-1]=='+'||line[k1-1]=='-'||line[k1-1]=='*'||line[k1-1]=='/'||line[k1-1]=='^'||line[k1-1]=='('||line[k1-1]==',') k2++;//printf("k1 %d line %s\n",k1,line+k1-1);}
          }
          strcpy(line1,&line[k1]);
          sum_extract(line1);
          p0=strchr(line1,'(');
          p0++;
          p1=strchr(line1,',');
          strncpy(argu,p0,p1-p0);
          argu[p1-p0]='\0';
          for(l=0; l<vars[LinVars[LinIndx].LinVarIndx].size; l++) {
            if(strcmp(argu,LinVars[LinIndx].dimnames[l])==0) {
              p1++;
              p1=strchr(p1,',');
              p1++;
              strncpy(line2,p1,strlen(line1)-(p1-line1)-1);
              line2[strlen(line1)-(p1-line1)-1]='\0';
              str_replace_all(line,line1,line2);
              IsRemSum=true;
              break;
            }
          }
        }
        p0=strchr(line,'(');
        p0++;
        p1=strchr(line,',');
        strncpy(argu,p0,p1-p0);
        argu[p1-p0]='\0';
        for(l=0; l<vars[LinVars[LinIndx].LinVarIndx].size; l++) {
          if(strcmp(argu,LinVars[LinIndx].dimnames[l])==0) {
            p1++;
            p1=strchr(p1,',');
            p1++;
            line2[0]='(';
            strncpy(&line2[1],p1,strlen(line)-(p1-line)-1);
            line2[strlen(line)-(p1-line)]=')';
            line2[strlen(line)-(p1-line)+1]='\0';
            strcpy(line,line2);
            IsRemSum=true;
            break;
          }
        }
        if(IsRemSum) {
          str_replace_all(formulain,line3,line);
          i=i+k+strlen(line);
        }
        else i=i+k+strlen(line3);
      }
      else {
        str_replace_all(formulain,line,"0");
        i=i+k+1;
      }
    }
    else {
      i=i+k+4;
      readitem=formulain+i;
    }
  }
  return 0;
}

int eq_linvar_read(char *formulain,eq_var_ref *LinVars,int linindx,array_def *vars) {
  char line[TABREADLINE];
  int i1,i2,i3,l;//,l,d;
  line[0]='p';
  line[1]='_';
  line[2]='\0';
  l=strlen(formulain);
  strcat(line,LinVars[linindx].LinVarName);
  i3=0;
  i1=0;
  while (i3==0) {
    i1=str_find_ci(formulain+i1,line);
    if (i1==-1) break;
    if (i1==0) {
      if(vars[LinVars[linindx].LinVarIndx].size==0) {
        i2=l-1;
      }
      else
        i2=str_find_ci(formulain,"}");
      formulain[0]='1';
      memmove(formulain+1,formulain+i2+1,l-i2);
      i3=1;
    }
    else {
      if (formulain[i1-1]=='+'||formulain[i1-1]=='-'||formulain[i1-1]=='*'||formulain[i1-1]=='/'||formulain[i1-1]=='^'||formulain[i1-1]=='('||formulain[i1-1]==',') {
        if(vars[LinVars[linindx].LinVarIndx].size==0) {
          i2=strlen(line)-1;
        }
        else i2=str_find_ci(formulain+i1,"}");
        formulain[i1]='1';
        memmove(formulain+i1+1,formulain+i1+i2+1,l-i2-i1);
        i3=1;
      }
      else i1++;
    }
  }
  return 1;
}

int equation_order_read(char *fname, char *commsyntax,set_def *sets,dim_t nset,set_element *set_elems,array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar,elem_value *elem_vals,offset_t ncofvar,offset_t ncofele,closure_entry *closure_vals,bool *var_inter,array_def *eq_defs,bool *eq_intertemp,dim_t *eq_orderintra,dim_t *eq_orderreg,offset_t allregset,offset_t alltimeset,dim_t *orderintra,dim_t *orderreg) {
  FILE * filehandle;
  char tline[TABREADLINE],line[TABREADLINE],line1[TABREADLINE],linecopy[TABREADLINE];//,set1[NAMESIZE],set2[NAMESIZE];
  char vname[TABREADLINE],lintmp[TABREADLINE];//,*p1=NULL;
  char *readitem=NULL,*p=NULL,*p1=NULL;//,*p2=NULL,*varpnts;
  solve_real zerodivide=0;
  dim_t fdim,np,i4;
  offset_t rowindx,j,l,l1,lj,dcountdim1[4*MAXVARDIM],dcountdim2[4*MAXVARDIM],nloops,nloopslin,li3,l2,matrow,eqindx=0,nelem,l01,j01;//,sizelinvars,totlinvars,templinvars
  int sj,i,i3,nlinvars,lvar,lvar1,lvar2,lvar3,lvar4,dcount,fdimlin=0,leadlag,varindx1,varindx2,tempint;
  solve_real dimmat[4*MAXVARDIM];

  filehandle = fopen(fname,"r");
  matrow=0;

  while (tab_next_statement_resolved(commsyntax,filehandle,line,elem_vals,coefs,ncof,&zerodivide,TABREADLINE)) {
    if (strstr(line,"(default")==NULL) {
      str_replace_first(line, commsyntax, "");
      str_replace_first(line, "(linear)", "");
      while (str_replace_all(line,"  ", " "));
      while (str_replace_char(line, '[', '('));
      while (str_replace_char(line, ']', ')'));
      while (str_replace_char(line, '{', '('));
      while (str_replace_char(line, '}', ')'));
      strcpy(linecopy,line);
      fdim=str_count_ci(line, "(all,");
      if (fdim==0) {
        readitem = strtok(line+1," ");
        strcpy(eq_defs[eqindx].cofname,readitem);
        readitem = strtok(NULL,"=");
        strcpy(vname,readitem);
        strcpy(line,linecopy);
        readitem = strtok(line,"=");
        readitem = strtok(NULL,";");
        strcat(readitem,"-");
        strcat(readitem,"(");
        strcat(readitem,vname);
        strcat(readitem,")");
      }
      else {
        readitem = strtok(line+1,"(");
        strcpy(eq_defs[eqindx].cofname,readitem);
        strcpy(line,linecopy);
        i=str_rfind_ci(line, "(all,");
        readitem=line+i;
        readitem = strtok(readitem,")");
        readitem = strtok(NULL,"=");
        strcpy(vname,readitem);
        strcpy(line,linecopy);
        readitem = strtok(line,"=");
        readitem = strtok(NULL,";");
        strcat(readitem,"-");
        strcat(readitem,"(");
        strcat(readitem,vname);
        strcat(readitem,")");
      }
      while (str_replace_all(readitem," ", ""));
      while (formula_normalize(readitem)==1);
      leadlag_encode(readitem);
      strcpy(tline,readitem);
      strcpy(line1,readitem);
      np=str_count_ci(readitem,"p_");
      eq_var_ref *LinVars= (eq_var_ref *) calloc (np+1,sizeof(eq_var_ref));
      i3=0;
      lvar=0;
      for (i=0; i<np; i++) {
        varindx2=0;
        while(-1<0) {
          varindx1=str_find_ci(readitem+varindx2,"p_");
          if(varindx1==-1) break;
          varindx2=varindx2+varindx1;
          if(varindx2==0||readitem[varindx2-1]=='*'||readitem[varindx2-1]=='+'||readitem[varindx2-1]=='-'||readitem[varindx2-1]=='('||readitem[varindx2-1]==',') break;
          else varindx2++;
        }
        if(varindx1==-1) break;
        lvar=lvar+varindx2;
        readitem=readitem+varindx2+2;
        p=strpbrk(readitem,"{+*-/^)");
        if(p!=NULL) {
          strncpy(vname,readitem,p-readitem);
          vname[p-readitem]='\0';
          strcpy(LinVars[i3].LinVarName,vname);
          for (l=0; l<nvar; l++) {
            if (strcmp(vars[l].cofname,vname)==0) {
              LinVars[i3].LinVarIndx=l;
              break;
            }
          }
          p=strpbrk(readitem,"}+*-/^)");
          if (*p=='}') {
            p++;
          }
          strncpy(tline,readitem,p-readitem);
          tline[p-readitem]='\0';
          switch (vars[l].size) {
          case 0:
            break;
          case 1:
            p = strtok(tline,"{");
            p = strtok(NULL,"}");
            leadlag=0;
            parse_index_leadlag(p,&leadlag);
            if(leadlag!=0)var_inter[l]=true;//printf("var %s\n",ha_var[l].cofname);}
            strcpy(LinVars[i3].dimnames[0],p);
            strcpy(lintmp,"(all,");
            strcat(lintmp,p);
            strcat(lintmp,",");
            l1=str_find_ci(linecopy,lintmp);
            if (l1>-1) {
              p1=&linecopy[0]+l1;
              strncpy(LinVars[i3].dimsetnames[0],p1+strlen(lintmp),strchr(p1,')')-p1-strlen(lintmp));
            }
            else {
              strcpy(lintmp,"sum(");
              strcat(lintmp,p);
              strcat(lintmp,",");
              lvar1=str_count_ci(linecopy,lintmp);
              lvar3=str_find_ci(linecopy,lintmp);
              if (lvar1>1) for(lvar2=0; lvar2<lvar1; lvar2++) {
                  lvar4=str_find_ci(&linecopy[lvar3+4],lintmp);
                  if (lvar4>-1&&lvar4<lvar) {
                    lvar3=lvar3+lvar4+4;
                  }
                  else {
                    break;
                  }
                }
              p1=&linecopy[0]+lvar3;
              p1=p1+strlen(lintmp);
              strncpy(LinVars[i3].dimsetnames[0],p1,strchr(p1,',')-p1);
            }
            break;
          default:
            p = strtok(tline,"{");
            for (i4=0; i4<vars[l].size-1; i4++) {
              p = strtok(NULL,",");
              leadlag=0;
              parse_index_leadlag(p,&leadlag);
              if(leadlag!=0)var_inter[l]=true;//printf("var %s\n",ha_var[l].cofname);}
              strcpy(LinVars[i3].dimnames[i4],p);
              strcpy(lintmp,"(all,");
              strcat(lintmp,p);
              strcat(lintmp,",");
              l1=str_find_ci(linecopy,lintmp);
              if (l1>-1) {
                p1=&linecopy[0]+l1;
                strncpy(LinVars[i3].dimsetnames[i4],p1+strlen(lintmp),strchr(p1,')')-p1-strlen(lintmp));
              }
              else {
                strcpy(lintmp,"sum(");
                strcat(lintmp,p);
                strcat(lintmp,",");
                lvar1=str_count_ci(linecopy,lintmp);
                lvar3=str_find_ci(linecopy,lintmp);
                if (lvar1>1) for(lvar2=0; lvar2<lvar1; lvar2++) {
                    lvar4=str_find_ci(&linecopy[lvar3+4],lintmp);
                    if (lvar4>-1&&lvar4<lvar) {
                      lvar3=lvar3+lvar4+4;
                    }
                    else {
                      break;
                    }
                  }
                p1=&linecopy[0]+lvar3;
                p1=p1+strlen(lintmp);
                strncpy(LinVars[i3].dimsetnames[i4],p1,strchr(p1,',')-p1);
              }
            }
            p = strtok(NULL,"}");
            leadlag=0;
            parse_index_leadlag(p,&leadlag);
            if(leadlag!=0)var_inter[l]=true;//printf("var %s\n",ha_var[l].cofname);}
            strcpy(LinVars[i3].dimnames[i4],p);
            strcpy(lintmp,"(all,");
            strcat(lintmp,p);
            strcat(lintmp,",");
            l1=str_find_ci(linecopy,lintmp);
            if (l1>-1) {
              p1=&linecopy[0]+l1;
              strncpy(LinVars[i3].dimsetnames[i4],p1+strlen(lintmp),strchr(p1,')')-p1-strlen(lintmp));
            }
            else {
              strcpy(lintmp,"sum(");
              strcat(lintmp,p);
              strcat(lintmp,",");
              lvar1=str_count_ci(linecopy,lintmp);
              lvar3=str_find_ci(linecopy,lintmp);
              if (lvar1>1) for(lvar2=0; lvar2<lvar1; lvar2++) {
                  lvar4=str_find_ci(&linecopy[lvar3+4],lintmp);
                  if (lvar4>-1&&lvar4<lvar) {
                    lvar3=lvar3+lvar4+4;
                  }
                  else {
                    break;
                  }
                }
              p1=&linecopy[0]+lvar3;
              p1=p1+strlen(lintmp);
              strncpy(LinVars[i3].dimsetnames[i4],p1,strchr(p1,',')-p1);
            }
            break;
          }
          i3++;
        }
        else {
          strcpy(vname,readitem);
          strcpy(LinVars[i3].LinVarName,vname);
          for (l=0; l<nvar; l++) {
            if (strcmp(vars[l].cofname,vname)==0) {
              LinVars[i3].LinVarIndx=l;
              break;
            }
          }
          i3++;
        }
      }
      nlinvars=i3;
      for (i=0; i<nlinvars; i++) i3+=vars[LinVars[i].LinVarIndx].size;
      quantifier *arSet= (quantifier *) calloc (fdim+i3+1,sizeof(quantifier));
      bool *linvarrcount= (bool *) calloc (nlinvars,sizeof(bool));

      strcpy(line,linecopy);
      nloops=1;
      if (fdim>0) {
        for (i=0; i<fdim; i++) {
          if(i==0) {
            readitem = strtok(line,",");
          }
          else {
            readitem = strtok(NULL,",");
          }
          readitem = strtok(NULL,",");
          strcpy(arSet[i].index_name,readitem);
          readitem = strtok(NULL,")");
          for (i4=0; i4<nset; i4++) if(strcmp(readitem,sets[i4].setname)==0) {
              arSet[i].setid=i4;
              eq_defs[eqindx].setid[i]=i4;
              break;
            }
          nloops=nloops*sets[arSet[i].setid].size;
        }
      }
      eq_defs[eqindx].size=fdim;
      eq_defs[eqindx].nelem=nloops;
      if(alltimeset>=0) {
        if(eqindx==0) {
          for(i=0; i<nvar; i++) {
            if(var_inter[i]) continue;
            j=0;
            for(i4=0; i4<vars[i].size; i4++)if(sets[vars[i].setid[i4]].intertemp) {
                orderintra[i]=i4;
                j++;
                break;
              }
            if(j==0)var_inter[i]=true;
          }
        }
        j=0;
        for (i=0; i<fdim; i++)if(sets[eq_defs[eqindx].setid[i]].intertemp) {
            j++;
            eq_orderintra[eqindx]=i;
          }
        if(j==0) {
          if(allregset<0)for (i4=0; i4<nlinvars; i4++) var_inter[LinVars[i4].LinVarIndx]=true;
          eq_intertemp[eqindx]=true;
        }

      }
      if(allregset>=0) {
        if(eqindx==0) {
          for(i=0; i<nvar; i++) {
            j=0;
            for(i4=0; i4<vars[i].size; i4++)if(sets[vars[i].setid[i4]].regional) {
                orderreg[i]=i4;
                j++;
              }
            if(j==0)var_inter[i]=true;
          }
        }
        j=0;
        for(i4=0; i4<nlinvars; i4++)for(i=0; i<vars[LinVars[i4].LinVarIndx].size; i++) if(sets[vars[LinVars[i4].LinVarIndx].setid[i]].regional)j++;
        if(j==0)eq_intertemp[eqindx]=true;
        j=0;
        for (i=0; i<fdim; i++)if(sets[eq_defs[eqindx].setid[i]].regional) {
            j++;
            eq_orderreg[eqindx]=i;
          }
        if(j==0) {
          eq_intertemp[eqindx]=true;
        }
        if(j>=2) {
          tempint=eq_orderreg[eqindx];
          eq_orderreg[eqindx]=-1;
          for (i=0; i<fdim; i++)dimmat[i]=-1;
          for (i=0; i<fdim; i++) {
            memset(linvarrcount,0,nlinvars*sizeof(bool));
            for(i4=0; i4<nlinvars; i4++) {
              for(lvar4=0; lvar4<i4; lvar4++)if(strcmp(LinVars[i4].LinVarName,LinVars[lvar4].LinVarName)==0)break;//{printf("here!!! linvarrcount[lvar4] %d\n",linvarrcount[lvar4]);
              if(lvar4<i4&&linvarrcount[lvar4]) {
                linvarrcount[i4]=true;
                continue;
              }
              for(i3=0; i3<vars[LinVars[i4].LinVarIndx].size; i3++) {
                if(strcmp(LinVars[i4].dimnames[i3],arSet[i].index_name)==0&&sets[eq_defs[eqindx].setid[i]].regional&&orderreg[LinVars[i4].LinVarIndx]==i3) { //ha_set[ha_var[LinVars[i4].LinVarIndx].setid[i3]].regional) {
                  dimmat[i]+=vars[LinVars[i4].LinVarIndx].nelem;
                  l01=0;
                  for(j01=0; j01<vars[LinVars[i4].LinVarIndx].nelem; j01++)if(closure_vals[vars[LinVars[i4].LinVarIndx].offset+j01].is_exogenous)l01++;
                  dimmat[i]-=l01;
                  linvarrcount[i4]=true;
                  if(lvar4<i4)linvarrcount[lvar4]=true;
                }
              }
            }
          }
          nelem=0;
          for (i=0; i<fdim; i++)if(dimmat[i]>nelem)nelem=dimmat[i];
          for (i=0; i<fdim; i++)if(dimmat[i]==nelem)break;
          if(nelem==0) {
            eq_orderreg[eqindx]=tempint;
            eq_intertemp[eqindx]=true;
          }
          else {
            eq_orderreg[eqindx]=i;
          }

        }


        if(eq_orderreg[eqindx]>-1) {
          nelem=0;
          memset(linvarrcount,0,nlinvars*sizeof(bool));
          for(i4=0; i4<nlinvars; i4++) {
            for(lvar4=0; lvar4<i4; lvar4++)if(strcmp(LinVars[i4].LinVarName,LinVars[lvar4].LinVarName)==0)break;
            if(lvar4<i4&&linvarrcount[lvar4]) {
              linvarrcount[i4]=true;
              continue;
            }
            if(strcmp(LinVars[i4].dimnames[orderreg[LinVars[i4].LinVarIndx]],arSet[eq_orderreg[eqindx]].index_name)!=0&&!var_inter[LinVars[i4].LinVarIndx]) {
              nelem+=vars[LinVars[i4].LinVarIndx].nelem;
              l01=0;
              for(j01=0; j01<vars[LinVars[i4].LinVarIndx].nelem; j01++)if(closure_vals[vars[LinVars[i4].LinVarIndx].offset+j01].is_exogenous)l01++;
              nelem+=vars[LinVars[i4].LinVarIndx].nelem-l01;
              linvarrcount[i4]=true;
              if(lvar4<i4)linvarrcount[lvar4]=true;
            }
          }
          if(nelem>=eq_defs[eqindx].nelem)eq_intertemp[eqindx]=true;
          else {
            for(i4=0; i4<nlinvars; i4++)
              if(strcmp(LinVars[i4].dimnames[orderreg[LinVars[i4].LinVarIndx]],arSet[eq_orderreg[eqindx]].index_name)!=0&&!var_inter[LinVars[i4].LinVarIndx]) {
                var_inter[LinVars[i4].LinVarIndx]=true;
              }
          }
        }
      }
      free(arSet);
      free(linvarrcount);
      free(LinVars);
      eqindx++;
    }
  }
  fclose(filehandle);
  return 1;
}

int equation_order_read_nested(char *fname, char *commsyntax,set_def *sets,dim_t nset,set_element *set_elems,array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar,elem_value *elem_vals,offset_t ncofvar,offset_t ncofele,closure_entry *closure_vals,bool *var_inter,array_def *eq_defs,bool *eq_intertemp,dim_t *eq_orderintra,dim_t *eq_orderreg,offset_t allregset,offset_t alltimeset,dim_t *orderintra,dim_t *orderreg) {
  FILE * filehandle;
  char tline[TABREADLINE],line[TABREADLINE],line1[TABREADLINE],linecopy[TABREADLINE];//,set1[NAMESIZE],set2[NAMESIZE];
  char vname[TABREADLINE],lintmp[TABREADLINE];//,*p1=NULL;
  char *readitem=NULL,*p=NULL,*p1=NULL;//,*p2=NULL,*varpnts;
  solve_real zerodivide=0;
  dim_t fdim,np;
  offset_t rowindx,j,j01,l,l01,lj,dcountdim1[4*MAXVARDIM],dcountdim2[4*MAXVARDIM],nloops,nloopslin,li3,l2,matrow,eqindx=0,nelem;//,sizelinvars,totlinvars,templinvars
  offset_t sj,i,i3,i4,nlinvars,lvar,lvar1,lvar2,lvar3,lvar4,dcount,fdimlin=0,varindx1,varindx2,tempint;
  int leadlag;
  solve_real dimmat[4*MAXVARDIM];

  filehandle = fopen(fname,"r");
  matrow=0;

  while (tab_next_statement_resolved(commsyntax,filehandle,line,elem_vals,coefs,ncof,&zerodivide,TABREADLINE)) {
    if (strstr(line,"(default")==NULL) {
      str_replace_first(line, commsyntax, "");
      str_replace_first(line, "(linear)", "");
      while (str_replace_all(line,"  ", " "));
      while (str_replace_char(line, '[', '('));
      while (str_replace_char(line, ']', ')'));
      while (str_replace_char(line, '{', '('));
      while (str_replace_char(line, '}', ')'));
      strcpy(linecopy,line);
      fdim=str_count_ci(line, "(all,");
      if (fdim==0) {
        readitem = strtok(line+1," ");
        strcpy(eq_defs[eqindx].cofname,readitem);
        readitem = strtok(NULL,"=");
        strcpy(vname,readitem);
        strcpy(line,linecopy);
        readitem = strtok(line,"=");
        readitem = strtok(NULL,";");
        strcat(readitem,"-");
        strcat(readitem,"(");
        strcat(readitem,vname);
        strcat(readitem,")");
      }
      else {
        readitem = strtok(line+1,"(");
        strcpy(eq_defs[eqindx].cofname,readitem);
        strcpy(line,linecopy);
        i=str_rfind_ci(line, "(all,");
        readitem=line+i;
        readitem = strtok(readitem,")");
        readitem = strtok(NULL,"=");
        strcpy(vname,readitem);
        strcpy(line,linecopy);
        readitem = strtok(line,"=");
        readitem = strtok(NULL,";");
        strcat(readitem,"-");
        strcat(readitem,"(");
        strcat(readitem,vname);
        strcat(readitem,")");
      }
      while (str_replace_all(readitem," ", ""));
      while (formula_normalize(readitem)==1);
      leadlag_encode(readitem);
      strcpy(tline,readitem);
      strcpy(line1,readitem);
      np=str_count_ci(readitem,"p_");
      eq_var_ref *LinVars= (eq_var_ref *) calloc (np+1,sizeof(eq_var_ref));
      i3=0;
      lvar=0;
      for (i=0; i<np; i++) {
        varindx2=0;
        while(-1<0) {
          varindx1=str_find_ci(readitem+varindx2,"p_");
          if(varindx1==-1) break;
          varindx2=varindx2+varindx1;
          if(varindx2==0||readitem[varindx2-1]=='*'||readitem[varindx2-1]=='+'||readitem[varindx2-1]=='-'||readitem[varindx2-1]=='('||readitem[varindx2-1]==',') break;
          else varindx2++;
        }
        if(varindx1==-1) break;
        lvar=lvar+varindx2;
        readitem=readitem+varindx2+2;
        p=strpbrk(readitem,"{+*-/^)");
        if(p!=NULL) {
          strncpy(vname,readitem,p-readitem);
          vname[p-readitem]='\0';
          strcpy(LinVars[i3].LinVarName,vname);
          for (l=0; l<nvar; l++) {
            if (strcmp(vars[l].cofname,vname)==0) {
              LinVars[i3].LinVarIndx=l;
              break;
            }
          }
          p=strpbrk(readitem,"}+*-/^)");
          if (*p=='}') {
            p++;
          }
          strncpy(tline,readitem,p-readitem);
          tline[p-readitem]='\0';
          switch (vars[l].size) {
          case 0:
            break;
          case 1:
            p = strtok(tline,"{");
            p = strtok(NULL,"}");
            leadlag=0;
            parse_index_leadlag(p,&leadlag);
            if(leadlag!=0) {
              var_inter[l]=true;
              orderintra[l]=-1;
              orderreg[l]=-1;
            }//printf("var %s\n",ha_var[l].cofname);}
            strcpy(LinVars[i3].dimnames[0],p);
            strcpy(lintmp,"(all,");
            strcat(lintmp,p);
            strcat(lintmp,",");
            l01=str_find_ci(linecopy,lintmp);
            if (l01>-1) {
              p1=&linecopy[0]+l01;
              strncpy(LinVars[i3].dimsetnames[0],p1+strlen(lintmp),strchr(p1,')')-p1-strlen(lintmp));
            }
            else {
              strcpy(lintmp,"sum(");
              strcat(lintmp,p);
              strcat(lintmp,",");
              lvar1=str_count_ci(linecopy,lintmp);
              lvar3=str_find_ci(linecopy,lintmp);
              if (lvar1>1) for(lvar2=0; lvar2<lvar1; lvar2++) {
                  lvar4=str_find_ci(&linecopy[lvar3+4],lintmp);
                  if (lvar4>-1&&lvar4<lvar) {
                    lvar3=lvar3+lvar4+4;
                  }
                  else {
                    break;
                  }
                }
              p1=&linecopy[0]+lvar3;
              p1=p1+strlen(lintmp);
              strncpy(LinVars[i3].dimsetnames[0],p1,strchr(p1,',')-p1);
            }
            break;
          default:
            p = strtok(tline,"{");
            for (i4=0; i4<vars[l].size-1; i4++) {
              p = strtok(NULL,",");
              leadlag=0;
              parse_index_leadlag(p,&leadlag);
              if(leadlag!=0) {
                var_inter[l]=true;
                orderintra[l]=-1;
                orderreg[l]=-1;
              }//printf("var %s\n",ha_var[l].cofname);}
              strcpy(LinVars[i3].dimnames[i4],p);
              strcpy(lintmp,"(all,");
              strcat(lintmp,p);
              strcat(lintmp,",");
              l01=str_find_ci(linecopy,lintmp);
              if (l01>-1) {
                p1=&linecopy[0]+l01;
                strncpy(LinVars[i3].dimsetnames[i4],p1+strlen(lintmp),strchr(p1,')')-p1-strlen(lintmp));
              }
              else {
                strcpy(lintmp,"sum(");
                strcat(lintmp,p);
                strcat(lintmp,",");
                lvar1=str_count_ci(linecopy,lintmp);
                lvar3=str_find_ci(linecopy,lintmp);
                if (lvar1>1) for(lvar2=0; lvar2<lvar1; lvar2++) {
                    lvar4=str_find_ci(&linecopy[lvar3+4],lintmp);
                    if (lvar4>-1&&lvar4<lvar) {
                      lvar3=lvar3+lvar4+4;
                    }
                    else {
                      break;
                    }
                  }
                p1=&linecopy[0]+lvar3;
                p1=p1+strlen(lintmp);
                strncpy(LinVars[i3].dimsetnames[i4],p1,strchr(p1,',')-p1);
              }
            }
            p = strtok(NULL,"}");
            leadlag=0;
            parse_index_leadlag(p,&leadlag);
            if(leadlag!=0) {
              var_inter[l]=true;
              orderintra[l]=-1;
              orderreg[l]=-1;
            }//printf("var %s\n",ha_var[l].cofname);}
            strcpy(LinVars[i3].dimnames[i4],p);
            strcpy(lintmp,"(all,");
            strcat(lintmp,p);
            strcat(lintmp,",");
            l01=str_find_ci(linecopy,lintmp);
            if (l01>-1) {
              p1=&linecopy[0]+l01;
              strncpy(LinVars[i3].dimsetnames[i4],p1+strlen(lintmp),strchr(p1,')')-p1-strlen(lintmp));
            }
            else {
              strcpy(lintmp,"sum(");
              strcat(lintmp,p);
              strcat(lintmp,",");
              lvar1=str_count_ci(linecopy,lintmp);
              lvar3=str_find_ci(linecopy,lintmp);
              if (lvar1>1) for(lvar2=0; lvar2<lvar1; lvar2++) {
                  lvar4=str_find_ci(&linecopy[lvar3+4],lintmp);
                  if (lvar4>-1&&lvar4<lvar) {
                    lvar3=lvar3+lvar4+4;
                  }
                  else {
                    break;
                  }
                }
              p1=&linecopy[0]+lvar3;
              p1=p1+strlen(lintmp);
              strncpy(LinVars[i3].dimsetnames[i4],p1,strchr(p1,',')-p1);
            }
            break;
          }
          i3++;
        }
        else {
          strcpy(vname,readitem);
          strcpy(LinVars[i3].LinVarName,vname);
          for (l=0; l<nvar; l++) {
            if (strcmp(vars[l].cofname,vname)==0) {
              LinVars[i3].LinVarIndx=l;
              break;
            }
          }
          i3++;
        }
      }
      nlinvars=i3;
      for (i=0; i<nlinvars; i++) i3+=vars[LinVars[i].LinVarIndx].size;
      quantifier *arSet= (quantifier *) calloc (fdim+i3+1,sizeof(quantifier));
      bool *linvarrcount= (bool *) calloc (nlinvars,sizeof(bool));

      strcpy(line,linecopy);
      nloops=1;
      if (fdim>0) {
        for (i=0; i<fdim; i++) {
          if(i==0) {
            readitem = strtok(line,",");
          }
          else {
            readitem = strtok(NULL,",");
          }
          readitem = strtok(NULL,",");
          strcpy(arSet[i].index_name,readitem);
          readitem = strtok(NULL,")");
          for (i4=0; i4<nset; i4++) if(strcmp(readitem,sets[i4].setname)==0) {
              arSet[i].setid=i4;
              eq_defs[eqindx].setid[i]=i4;
              break;
            }
          nloops=nloops*sets[arSet[i].setid].size;
        }
      }
      eq_defs[eqindx].size=fdim;
      eq_defs[eqindx].nelem=nloops;


      if(eqindx==0) {
        for(i=0; i<nvar; i++) {
          if(var_inter[i]) continue;
          j=0;
          for(i4=0; i4<vars[i].size; i4++)if(sets[vars[i].setid[i4]].intertemp) {
              orderintra[i]=i4;
              j++;
              break;
            }
          if(j==0) {
            var_inter[i]=true;
            orderreg[i]=-1;
            orderintra[i]=-1;
          }
          else {
            j01=0;
            for(i4=0; i4<vars[i].size; i4++)if(sets[vars[i].setid[i4]].regional) {
                orderreg[i]=i4;
                j01++;
              }
          }
        }
      }


      j=0;
      for(i4=0; i4<nlinvars; i4++)for(i=0; i<vars[LinVars[i4].LinVarIndx].size; i++) if(sets[vars[LinVars[i4].LinVarIndx].setid[i]].intertemp)j++;
      if(j==0) {
        eq_intertemp[eqindx]=true;//eq_orderreg[eqindx]=-1;eq_orderintra[eqindx]=-1;}
      }
      j=0;
      for (i=0; i<fdim; i++)if(sets[eq_defs[eqindx].setid[i]].regional) {
          j++;
          eq_orderreg[eqindx]=i;
        }
      if(j>=2) {
        for (i=0; i<fdim; i++)dimmat[i]=-1;
        for (i=0; i<fdim; i++) {
          memset(linvarrcount,0,nlinvars*sizeof(bool));
          for(i4=0; i4<nlinvars; i4++) {
            for(lvar4=0; lvar4<i4; lvar4++)if(strcmp(LinVars[i4].LinVarName,LinVars[lvar4].LinVarName)==0)break;
            if(lvar4<i4&&linvarrcount[lvar4]) {
              linvarrcount[i4]=true;
              continue;
            }
            for(i3=0; i3<vars[LinVars[i4].LinVarIndx].size; i3++) {
              if(strcmp(LinVars[i4].dimnames[i3],arSet[i].index_name)==0&&sets[eq_defs[eqindx].setid[i]].regional&&orderreg[LinVars[i4].LinVarIndx]==i3) {
                dimmat[i]+=vars[LinVars[i4].LinVarIndx].nelem;
                linvarrcount[i4]=true;
                if(lvar4<i4)linvarrcount[lvar4]=true;
                l01=0;
                for(j01=0; j01<vars[LinVars[i4].LinVarIndx].nelem; j01++)if(closure_vals[vars[LinVars[i4].LinVarIndx].offset+j01].is_exogenous)l01++;
                dimmat[i]-=l01;
              }
            }
          }
        }
        nelem=0;
        for (i=0; i<fdim; i++)if(dimmat[i]>nelem)nelem=dimmat[i];
        for (i=0; i<fdim; i++)if(dimmat[i]==nelem)break;
        if(nelem==0) {
          eq_orderreg[eqindx]=-1;
        }
        else {
          eq_orderreg[eqindx]=i;
        }
      }
      if(eq_orderreg[eqindx]>-1) {
        nelem=0;
        memset(linvarrcount,0,nlinvars*sizeof(bool));
        for(i4=0; i4<nlinvars; i4++) {
          for(lvar4=0; lvar4<i4; lvar4++)if(strcmp(LinVars[i4].LinVarName,LinVars[lvar4].LinVarName)==0)break;
          if(lvar4<i4&&linvarrcount[lvar4]) {
            linvarrcount[i4]=true;
            continue;
          }
          if(strcmp(LinVars[i4].dimnames[orderreg[LinVars[i4].LinVarIndx]],arSet[eq_orderreg[eqindx]].index_name)!=0&&orderreg[LinVars[i4].LinVarIndx]!=-1) {
            l01=0;
            for(j01=0; j01<vars[LinVars[i4].LinVarIndx].nelem; j01++)if(closure_vals[vars[LinVars[i4].LinVarIndx].offset+j01].is_exogenous)l01++;
            nelem+=vars[LinVars[i4].LinVarIndx].nelem-l01;
            linvarrcount[i4]=true;
            if(lvar4<i4)linvarrcount[lvar4]=true;
          }
        }
        if(nelem>=eq_defs[eqindx].nelem)eq_orderreg[eqindx]=-1;//printf("OOOOOOOOOOOOOOOOO\n");}
        else {
          for(i4=0; i4<nlinvars; i4++)
            if(strcmp(LinVars[i4].dimnames[orderreg[LinVars[i4].LinVarIndx]],arSet[eq_orderreg[eqindx]].index_name)!=0&&orderreg[LinVars[i4].LinVarIndx]!=-1) {
              orderreg[LinVars[i4].LinVarIndx]=-1;//orderintra[LinVars[i4].LinVarIndx]=-1;
            }
        }
      }


      j=0;
      for (i=0; i<fdim; i++)if(sets[eq_defs[eqindx].setid[i]].intertemp) {
          j++;
          eq_orderintra[eqindx]=i;
        }
      if(j==0) {
        eq_intertemp[eqindx]=true;
      }

      if(eq_intertemp[eqindx]==true) {
        eq_orderreg[eqindx]=-1;
        eq_orderintra[eqindx]=-1;
      }

      free(arSet);
      free(LinVars);
      free(linvarrcount);
      eqindx++;
    }
  }
  fclose(filehandle);
  return 1;
}

int jacobian_preallocate(char *fname, char *commsyntax,set_def *sets,dim_t nset,set_element *set_elems,array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar,elem_value *elem_vals,offset_t ncofvar,offset_t ncofele, offset_t nexo,closure_entry *closure_vals,offset_t ndblock,offset_t alltimeset,offset_t allregset,bool *eq_intertemp,offset_t *eq_addr,dim_t *eq_time,dim_t *eq_reg,offset_t *counteq,offset_t nintraeq,bool *sbbd_overrid,PetscInt Istart,PetscInt Iend,PetscInt *dnz,PetscInt *dnnz,PetscInt *onz,PetscInt *onnz,PetscInt *dnzB,PetscInt *dnnzB,PetscInt *onzB,PetscInt *onnzB,int nesteddbbd) {
  FILE * filehandle;
  char tline[TABREADLINE],line[TABREADLINE],line1[TABREADLINE],linecopy[TABREADLINE];//,set1[NAMESIZE],set2[NAMESIZE];
  char vname[TABREADLINE],lintmp[TABREADLINE];//,*p1=NULL;
  char *readitem=NULL,*p=NULL,*p1=NULL;//,*p2=NULL,*varpnts;
  PetscInt Iindx=0,Jindx;
  solve_real zerodivide=0;
  dim_t fdim=0,np,i4,sup,supset[MAXSUPSET];
  offset_t rowindx,rowindxorg,l,l1,lj,dcountdim1[4*MAXVARDIM],dcountdim2[4*MAXVARDIM],dcountdim3[4*MAXVARDIM],dcountdim4[4*MAXVARDIM],dcountdim5[4*MAXVARDIM],nloops,nloopslin,nloopsfac,li3,l2,matrow,matroworg,ltime,lreg,leq=0,eqindx=0;//,sizelinvars,totlinvars,templinvars
  offset_t nreg=0,nint=0,sj,i,i3;
  if(allregset>-1)nreg=sets[allregset].size;
  if(alltimeset>-1)nint=sets[alltimeset].size;
  int nlinvars,lvar,lvar1,lvar2,lvar3,lvar4,dcount,fdimlin=0,leadlag,varindx1,varindx2;
  offset_t *counteq1= (offset_t *) calloc (ndblock,sizeof(offset_t));
  for(i=0; i<ndblock; i++)counteq1[i]=counteq[i];
  filehandle = fopen(fname,"r");
  matrow=nintraeq;
  matroworg=0;
  Jindx=0;

  while (tab_next_statement_resolved(commsyntax,filehandle,line,elem_vals,coefs,ncof,&zerodivide,TABREADLINE)) {
    if (strstr(line,"(default")==NULL) {
      str_replace_first(line, commsyntax, "");
      str_replace_first(line, "(linear)", "");
      while (str_replace_all(line,"  ", " "));
      while (str_replace_char(line, '[', '('));
      while (str_replace_char(line, ']', ')'));
      while (str_replace_char(line, '{', '('));
      while (str_replace_char(line, '}', ')'));
      strcpy(linecopy,line);
      fdim=str_count_ci(line, "(all,");
      if (fdim==0) {
        readitem = strtok(line+1," ");
        readitem = strtok(NULL,"=");
        strcpy(vname,readitem);
        strcpy(line,linecopy);
        readitem = strtok(line,"=");
        readitem = strtok(NULL,";");
        strcat(readitem,"-");
        strcat(readitem,"(");
        strcat(readitem,vname);
        strcat(readitem,")");
      }
      else {
        i=str_rfind_ci(line, "(all,");
        readitem=line+i;
        readitem = strtok(readitem,")");
        readitem = strtok(NULL,"=");
        strcpy(vname,readitem);
        strcpy(line,linecopy);
        readitem = strtok(line,"=");
        readitem = strtok(NULL,";");
        strcat(readitem,"-");
        strcat(readitem,"(");
        strcat(readitem,vname);
        strcat(readitem,")");
      }
      while (str_replace_all(readitem," ", ""));
      while (formula_normalize(readitem)==1);
      leadlag_encode(readitem);
      strcpy(tline,readitem);
      strcpy(line1,readitem);
      np=str_count_ci(readitem,"p_");
      eq_var_ref *LinVars= (eq_var_ref *) calloc (np+1,sizeof(eq_var_ref));
      i3=0;
      lvar=0;
      for (i=0; i<np; i++) {
        varindx2=0;
        while(-1<0) {
          varindx1=str_find_ci(readitem+varindx2,"p_");
          if(varindx1==-1) break;
          varindx2=varindx2+varindx1;
          if(varindx2==0||readitem[varindx2-1]=='*'||readitem[varindx2-1]=='+'||readitem[varindx2-1]=='-'||readitem[varindx2-1]=='('||readitem[varindx2-1]==',') break;
          else varindx2++;
        }
        if(varindx1==-1) break;
        lvar=lvar+varindx2;
        readitem=readitem+varindx2+2;
        p=strpbrk(readitem,"{+*-/^)");
        if(p!=NULL) {
          strncpy(vname,readitem,p-readitem);
          vname[p-readitem]='\0';
          strcpy(LinVars[i3].LinVarName,vname);
          for (l=0; l<nvar; l++) {
            if (strcmp(vars[l].cofname,vname)==0) {
              LinVars[i3].LinVarIndx=l;
              break;
            }
          }
          p=strpbrk(readitem,"}+*-/^)");
          if (*p=='}') {
            p++;
          }
          strncpy(tline,readitem,p-readitem);
          tline[p-readitem]='\0';
          switch (vars[l].size) {
          case 0:
            break;
          case 1:
            p = strtok(tline,"{");
            p = strtok(NULL,"}");
            leadlag=0;
            parse_index_leadlag(p,&leadlag);
            strcpy(LinVars[i3].dimnames[0],p);
            LinVars[i3].dimleadlag[0]=leadlag;
            strcpy(lintmp,"(all,");
            strcat(lintmp,p);
            strcat(lintmp,",");
            l1=str_find_ci(linecopy,lintmp);
            if (l1>-1) {
              p1=&linecopy[0]+l1;
              strncpy(LinVars[i3].dimsetnames[0],p1+strlen(lintmp),strchr(p1,')')-p1-strlen(lintmp));
            }
            else {
              strcpy(lintmp,"sum(");
              strcat(lintmp,p);
              strcat(lintmp,",");
              lvar1=str_count_ci(linecopy,lintmp);
              lvar3=str_find_ci(linecopy,lintmp);
              if (lvar1>1) for(lvar2=0; lvar2<lvar1; lvar2++) {
                  lvar4=str_find_ci(&linecopy[lvar3+4],lintmp);
                  if (lvar4>-1&&lvar4<lvar) {
                    lvar3=lvar3+lvar4+4;
                  }
                  else {
                    break;
                  }
                }
              p1=&linecopy[0]+lvar3;
              p1=p1+strlen(lintmp);
              strncpy(LinVars[i3].dimsetnames[0],p1,strchr(p1,',')-p1);
            }
            break;
          default:
            p = strtok(tline,"{");
            for (i4=0; i4<vars[l].size-1; i4++) {
              p = strtok(NULL,",");
              leadlag=0;
              parse_index_leadlag(p,&leadlag);
              strcpy(LinVars[i3].dimnames[i4],p);
              LinVars[i3].dimleadlag[i4]=leadlag;
              strcpy(lintmp,"(all,");
              strcat(lintmp,p);
              strcat(lintmp,",");
              l1=str_find_ci(linecopy,lintmp);
              if (l1>-1) {
                p1=&linecopy[0]+l1;
                strncpy(LinVars[i3].dimsetnames[i4],p1+strlen(lintmp),strchr(p1,')')-p1-strlen(lintmp));
              }
              else {
                strcpy(lintmp,"sum(");
                strcat(lintmp,p);
                strcat(lintmp,",");
                lvar1=str_count_ci(linecopy,lintmp);
                lvar3=str_find_ci(linecopy,lintmp);
                if (lvar1>1) for(lvar2=0; lvar2<lvar1; lvar2++) {
                    lvar4=str_find_ci(&linecopy[lvar3+4],lintmp);
                    if (lvar4>-1&&lvar4<lvar) {
                      lvar3=lvar3+lvar4+4;
                    }
                    else {
                      break;
                    }
                  }
                p1=&linecopy[0]+lvar3;
                p1=p1+strlen(lintmp);
                strncpy(LinVars[i3].dimsetnames[i4],p1,strchr(p1,',')-p1);
              }
            }
            p = strtok(NULL,"}");
            leadlag=0;
            parse_index_leadlag(p,&leadlag);
            strcpy(LinVars[i3].dimnames[i4],p);
            LinVars[i3].dimleadlag[i4]=leadlag;
            strcpy(lintmp,"(all,");
            strcat(lintmp,p);
            strcat(lintmp,",");
            l1=str_find_ci(linecopy,lintmp);
            if (l1>-1) {
              p1=&linecopy[0]+l1;
              strncpy(LinVars[i3].dimsetnames[i4],p1+strlen(lintmp),strchr(p1,')')-p1-strlen(lintmp));
            }
            else {
              strcpy(lintmp,"sum(");
              strcat(lintmp,p);
              strcat(lintmp,",");
              lvar1=str_count_ci(linecopy,lintmp);
              lvar3=str_find_ci(linecopy,lintmp);
              if (lvar1>1) for(lvar2=0; lvar2<lvar1; lvar2++) {
                  lvar4=str_find_ci(&linecopy[lvar3+4],lintmp);
                  if (lvar4>-1&&lvar4<lvar) {
                    lvar3=lvar3+lvar4+4;
                  }
                  else {
                    break;
                  }
                }
              p1=&linecopy[0]+lvar3;
              p1=p1+strlen(lintmp);
              strncpy(LinVars[i3].dimsetnames[i4],p1,strchr(p1,',')-p1);
            }
            break;
          }
          i3++;
        }
        else {
          strcpy(vname,readitem);
          strcpy(LinVars[i3].LinVarName,vname);
          for (l=0; l<nvar; l++) {
            if (strcmp(vars[l].cofname,vname)==0) {
              LinVars[i3].LinVarIndx=l;
              break;
            }
          }
          i3++;
        }
      }
      nlinvars=i3;
      for (i=0; i<nlinvars; i++) i3+=vars[LinVars[i].LinVarIndx].size;
      quantifier *arSet= (quantifier *) calloc (fdim+i3+1,sizeof(quantifier));

      strcpy(line,linecopy);
      nloops=1;
      if (fdim>0) {
        for (i=0; i<fdim; i++) {
          if(i==0) {
            readitem = strtok(line,",");
          }
          else {
            readitem = strtok(NULL,",");
          }
          readitem = strtok(NULL,",");
          strcpy(arSet[i].index_name,readitem);
          readitem = strtok(NULL,")");
          for (i4=0; i4<nset; i4++) if(strcmp(readitem,sets[i4].setname)==0) {
              arSet[i].setid=i4;
              if(sets[i4].intertemp) *sbbd_overrid=true;
              break;
            }
          nloops=nloops*sets[arSet[i].setid].size;
        }
        //*****Adjust ha_eqreg*******************//
        //*****End Adjust ha_eqreg*******************//
        dcountdim1[fdim-1]=1;
        for (dcount=fdim-2; dcount>-1; dcount--) {
          dcountdim1[dcount]=sets[arSet[dcount+1].setid].size*dcountdim1[dcount+1];
        }
        if(alltimeset>=0&&allregset<0) {
          for (dcount=eq_time[eqindx]-1; dcount>-1; dcount--) {
            dcountdim3[dcount]=dcountdim1[dcount]/sets[arSet[eq_time[eqindx]].setid].size;
          }
          for (dcount=fdim-1; dcount>eq_time[eqindx]-1; dcount--) {
            dcountdim3[dcount]=dcountdim1[dcount];
          }
        }
        if(alltimeset<0&&allregset>=0) {
          for (dcount=eq_reg[eqindx]-1; dcount>-1; dcount--) {
            dcountdim4[dcount]=dcountdim1[dcount]/sets[arSet[eq_reg[eqindx]].setid].size;
          }
          for (dcount=fdim-1; dcount>eq_reg[eqindx]-1; dcount--) {
            dcountdim4[dcount]=dcountdim1[dcount];
          }
        }
        if(alltimeset>=0&&allregset>=0) {
          i3=1;
          dcountdim4[fdim-1]=dcountdim1[fdim-1];
          for (dcount=fdim-2; dcount>-1; dcount--) {
            if(dcount+1==eq_reg[eqindx])i3*=sets[arSet[eq_reg[eqindx]].setid].size;
            if(dcount+1==eq_time[eqindx])i3*=sets[arSet[eq_time[eqindx]].setid].size;
            dcountdim4[dcount]=dcountdim1[dcount]/i3;
          }
        }
      }
      for (i=0; i<nlinvars; i++) {
        i3=0;
        nloopslin=nloops;
        if (fdim==0) {
          for (l2=0; l2<vars[LinVars[i].LinVarIndx].size; l2++) {
            strcpy(arSet[i3].index_name,LinVars[i].dimnames[l2]);
            for (i4=0; i4<nset; i4++) if(strcmp(sets[i4].setname,LinVars[i].dimsetnames[l2])==0) {
                break;
              }
            arSet[i3].setid=i4;
            nloopslin=nloopslin*sets[arSet[i3].setid].size;
            if(sets[i4].intertemp)*sbbd_overrid=true;
            i3++;
          }
        }
        else for (l2=0; l2<vars[LinVars[i].LinVarIndx].size; l2++) {
            for (sj=0; sj<fdim; sj++) {
              if(strcmp(arSet[sj].index_name,LinVars[i].dimnames[l2])==0) {
                break;
              }
              else if (sj==fdim-1) {
                strcpy(arSet[fdim+i3].index_name,LinVars[i].dimnames[l2]);
                for (i4=0; i4<nset; i4++) if(strcmp(sets[i4].setname,LinVars[i].dimsetnames[l2])==0) {
                    break;
                  }
                arSet[fdim+i3].setid=i4;
                nloopslin=nloopslin*sets[arSet[fdim+i3].setid].size;
                if(sets[i4].intertemp)*sbbd_overrid=true;
                i3++;
              }
            }
          }
        fdimlin=fdim+i3;
        if (i3>0) {
          dcountdim2[fdimlin-1]=1;
          for (dcount=fdimlin-2; dcount>-1; dcount--) {
            dcountdim2[dcount]=sets[arSet[dcount+1].setid].size*dcountdim2[dcount+1];
          }
        }
        else for (dcount=0; dcount<fdim; dcount++) {
            dcountdim2[dcount]=dcountdim1[dcount];
          }
        nloopsfac=(offset_t)nloopslin/nloops;
        for (dcount=0; dcount<vars[LinVars[i].LinVarIndx].size; dcount++) {
          for (i4=0; i4<fdimlin; i4++) {
            if (strcmp(LinVars[i].dimnames[dcount],arSet[i4].index_name)==0) {
              dcountdim5[dcount]=i4;
              break;
            }
          }
        }
        for(dcount=0; dcount<MAXSUPSET; dcount++)supset[dcount]=0;
        for (dcount=0; dcount<vars[LinVars[i].LinVarIndx].size; dcount++) {
          if(sets[vars[LinVars[i].LinVarIndx].setid[dcount]].size!=sets[arSet[dcountdim5[dcount]].setid].size) {
            for(sup=1; sup<MAXSUPSET; sup++)if(vars[LinVars[i].LinVarIndx].setid[dcount]==sets[arSet[dcountdim5[dcount]].setid].subsetid[sup]) {
                supset[dcount]=sup;
                break;
              }
          }
          else supset[dcount]=0;
        }
        for (lj=0; lj<nloopslin; lj++) {
          l2=lj;
          rowindx=0;
          rowindxorg=0;
          for (dcount=0; dcount<fdimlin; dcount++) {
            l1=(offset_t) l2/dcountdim2[dcount];
            arSet[dcount].indx=l1;
            if(eq_intertemp[eqindx]) {
              if(alltimeset>=0&&allregset<0) {
                if (dcount<fdim&&dcount!=eq_time[eqindx]) {
                  rowindx=rowindx+l1*dcountdim3[dcount];
                }
                if(dcount==eq_time[eqindx])ltime=l1;
              }
              if(alltimeset<0&&allregset>=0) {
                if (dcount<fdim&&dcount!=eq_reg[eqindx]) {
                  rowindx=rowindx+l1*dcountdim4[dcount];
                }
                if(dcount==eq_reg[eqindx])ltime=l1;
              }
              if(alltimeset>=0&&allregset>=0) {
                if(nesteddbbd==1) {
                  if(eq_reg[eqindx]>-1) {
                    if (dcount<fdim&&dcount!=eq_time[eqindx]&&dcount!=eq_reg[eqindx])rowindx=rowindx+l1*dcountdim4[dcount];
                  }
                  else {
                    if (dcount<fdim&&dcount!=eq_time[eqindx])rowindx=rowindx+l1*dcountdim4[dcount];
                  }
                  if(dcount==eq_time[eqindx])ltime=l1;
                  if(dcount==eq_reg[eqindx])lreg=l1;
                }
                else {
                  if (dcount<fdim&&dcount!=eq_time[eqindx]&&dcount!=eq_reg[eqindx]) {
                    rowindx=rowindx+l1*dcountdim4[dcount];
                  }
                  if(dcount==eq_time[eqindx])ltime=l1;
                  if(dcount==eq_reg[eqindx])lreg=l1;
                }
              }
            }
            else {
              if (dcount<fdim) {
                rowindx=rowindx+l1*dcountdim1[dcount];
              }
            }
            l2=l2-l1*dcountdim2[dcount];
          }
          li3=0;
          for (dcount=0; dcount<vars[LinVars[i].LinVarIndx].size; dcount++) {
            if(supset[dcount]==0) {
              li3=li3+(arSet[dcountdim5[dcount]].indx+LinVars[i].dimleadlag[dcount])*vars[LinVars[i].LinVarIndx].strides[dcount];
            }
            else {
              li3=li3+(set_elems[sets[arSet[dcountdim5[dcount]].setid].offset+arSet[dcountdim5[dcount]].indx].superset_pos[supset[dcount]]+LinVars[i].dimleadlag[dcount])*vars[LinVars[i].LinVarIndx].strides[dcount];
            }
          }

          if(eq_intertemp[eqindx]) {
            if(alltimeset>=0&&allregset>=0) {
              if(nesteddbbd==1) {
                if(eq_reg[eqindx]>-1)Jindx=counteq1[(set_elems[sets[arSet[eq_time[eqindx]].setid].offset+ltime].superset_pos[sets[arSet[eq_time[eqindx]].setid].intsup])*(nreg+1)+set_elems[sets[arSet[eq_reg[eqindx]].setid].offset+lreg].superset_pos[sets[arSet[eq_reg[eqindx]].setid].regsup]]+rowindx;
                else Jindx=counteq1[(set_elems[sets[arSet[eq_time[eqindx]].setid].offset+ltime].superset_pos[sets[arSet[eq_time[eqindx]].setid].intsup])*(nreg+1)+nreg]+rowindx;
              }
              else {
                if(arSet[eq_time[eqindx]].setid==alltimeset)
                  Jindx=counteq1[(set_elems[sets[arSet[eq_time[eqindx]].setid].offset+ltime].superset_pos[0])*sets[arSet[eq_reg[eqindx]].setid].size+set_elems[sets[arSet[eq_reg[eqindx]].setid].offset+lreg].superset_pos[0]]+rowindx;
                else {
                  for(i4=1; i4<MAXSUPSET; i4++)if(sets[arSet[eq_time[eqindx]].setid].subsetid[i4]=alltimeset)break;
                  Jindx=counteq1[(set_elems[sets[arSet[eq_time[eqindx]].setid].offset+ltime].superset_pos[i4])*sets[arSet[eq_reg[eqindx]].setid].size+set_elems[sets[arSet[eq_reg[eqindx]].setid].offset+lreg].superset_pos[0]]+rowindx;
                }
              }
            }
            if(alltimeset>=0&&allregset<0) {
              if(arSet[eq_time[eqindx]].setid==alltimeset)
                Jindx=counteq1[set_elems[sets[arSet[eq_time[eqindx]].setid].offset+ltime].superset_pos[0]]+rowindx;
              else {
                for(i4=1; i4<MAXSUPSET; i4++)if(sets[arSet[eq_time[eqindx]].setid].subsetid[i4]=alltimeset)break;
                Jindx=counteq1[set_elems[sets[arSet[eq_time[eqindx]].setid].offset+ltime].superset_pos[i4]]+rowindx;
              }
            }
            if(alltimeset<0&&allregset>=0) Jindx=counteq1[set_elems[sets[arSet[eq_reg[eqindx]].setid].offset+ltime].superset_pos[0]]+rowindx;
          }
          else Jindx=matrow+rowindx;
          leq=matroworg+(offset_t)lj/nloopsfac;//matroworg+rowindxorg;
          eq_addr[leq]=Jindx;
          Iindx=closure_vals[vars[LinVars[i].LinVarIndx].offset+li3].exo_index;
          if(Istart<=Jindx&&Jindx<Iend) {
            if (!closure_vals[vars[LinVars[i].LinVarIndx].offset+li3].is_exogenous) {
              if (Istart<=Iindx&&Iindx<Iend) {
                dnnz[Jindx-Istart]=dnnz[Jindx-Istart]+1;
              }
              else {
                onnz[Jindx-Istart]=onnz[Jindx-Istart]+1;
              }
            }
            if (closure_vals[vars[LinVars[i].LinVarIndx].offset+li3].is_exogenous) {
              if (Istart<=Iindx&&Iindx<Iend) { //&&Iindx<nexo
                dnnzB[Jindx-Istart]=dnnzB[Jindx-Istart]+1;
              }
              else {
                onnzB[Jindx-Istart]=onnzB[Jindx-Istart]+1;
              }
            }
          }
        }
      }
      if(eq_intertemp[eqindx]) {
        if(alltimeset>=0&&allregset>=0) {
          if(nesteddbbd==1) {
            if(eq_reg[eqindx]>-1)for(lj=0; lj<sets[arSet[eq_time[eqindx]].setid].size; lj++)for(l2=0; l2<sets[arSet[eq_reg[eqindx]].setid].size; l2++) {
                  counteq1[set_elems[sets[arSet[eq_time[eqindx]].setid].offset+lj].superset_pos[sets[arSet[eq_time[eqindx]].setid].intsup]*(nreg+1)+set_elems[sets[arSet[eq_reg[eqindx]].setid].offset+l2].superset_pos[sets[arSet[eq_reg[eqindx]].setid].regsup]]+=nloops/sets[arSet[eq_time[eqindx]].setid].size/sets[arSet[eq_reg[eqindx]].setid].size;
                  if(sets[arSet[eq_reg[eqindx]].setid].regsup!=0)printf("Error: subsets of the regional partition set are not supported in NDBBD ordering\n");
                }
            else for(lj=0; lj<sets[arSet[eq_time[eqindx]].setid].size; lj++)
                counteq1[set_elems[sets[arSet[eq_time[eqindx]].setid].offset+lj].superset_pos[sets[arSet[eq_time[eqindx]].setid].intsup]*(nreg+1)+nreg]+=nloops/sets[arSet[eq_time[eqindx]].setid].size;
          }
          else {
            if(arSet[eq_time[eqindx]].setid==alltimeset)
              for(lj=0; lj<sets[arSet[eq_time[eqindx]].setid].size; lj++)for(l2=0; l2<sets[arSet[eq_reg[eqindx]].setid].size; l2++)
                  counteq1[set_elems[sets[arSet[eq_time[eqindx]].setid].offset+lj].superset_pos[0]*sets[arSet[eq_reg[eqindx]].setid].size+set_elems[sets[arSet[eq_reg[eqindx]].setid].offset+l2].superset_pos[0]]+=nloops/sets[arSet[eq_time[eqindx]].setid].size/sets[arSet[eq_reg[eqindx]].setid].size;
            else {
              for(lj=0; lj<sets[arSet[eq_time[eqindx]].setid].size; lj++)for(l2=0; l2<sets[arSet[eq_reg[eqindx]].setid].size; l2++) {
                  for(i4=1; i4<MAXSUPSET; i4++)if(sets[arSet[eq_time[eqindx]].setid].subsetid[i4]=alltimeset)break;
                  counteq1[set_elems[sets[arSet[eq_time[eqindx]].setid].offset+lj].superset_pos[i4]*sets[arSet[eq_reg[eqindx]].setid].size+set_elems[sets[arSet[eq_reg[eqindx]].setid].offset+l2].superset_pos[0]]+=nloops/sets[arSet[eq_time[eqindx]].setid].size/sets[arSet[eq_reg[eqindx]].setid].size;
                }
            }
          }
        }
        if(alltimeset>=0&&allregset<0) {
          if(arSet[eq_time[eqindx]].setid==alltimeset)
            for(lj=0; lj<sets[arSet[eq_time[eqindx]].setid].size; lj++)counteq1[set_elems[sets[arSet[eq_time[eqindx]].setid].offset+lj].superset_pos[0]]+=nloops/sets[arSet[eq_time[eqindx]].setid].size;
          else {
            for(i4=1; i4<MAXSUPSET; i4++)if(sets[arSet[eq_time[eqindx]].setid].subsetid[i4]=alltimeset)break;
            for(lj=0; lj<sets[arSet[eq_time[eqindx]].setid].size; lj++)counteq1[set_elems[sets[arSet[eq_time[eqindx]].setid].offset+lj].superset_pos[i4]]+=nloops/sets[arSet[eq_time[eqindx]].setid].size;
          }
        }
        if(alltimeset<0&&allregset>=0)for(lj=0; lj<sets[arSet[eq_reg[eqindx]].setid].size; lj++)counteq1[set_elems[sets[arSet[eq_reg[eqindx]].setid].offset+lj].superset_pos[0]]+=nloops/sets[arSet[eq_reg[eqindx]].setid].size;
      }
      else matrow+=nloops;
      matroworg+=nloops;
      eqindx++;
      free(LinVars);
      free(arSet);
    }
  }
  free(counteq1);
  fclose(filehandle);
  return 1;
}

