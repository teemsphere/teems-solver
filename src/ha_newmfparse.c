#include <ha_cgeglobal.h>

extern void spec51_rank_(int *INSIZE,ha_cgetype *cntl6,int *IRN, int *JCN, ha_cgetype *VA);
extern void spec51m_rank_(int *INSIZE,ha_cgetype *cntl6,int *IRN, int *JCN, ha_cgetype *VA,int *IRNA, int *JCNA, int *KEEP,ha_cgetype *w51, int *iw51);
extern void spec48_ssol2la_(int *INSIZE,int *IRN, int *JCN, ha_cgetype *VA, ha_cgetype *B, ha_cgetype *X);
extern void spec48m_ssol2la_(int *INSIZE,int *IRN, int *JCN, ha_cgetype *VA, ha_cgetype *B, ha_cgetype *X);
extern void spec48_msol_(int *INSIZE,int *IRN, int *JCN, ha_cgetype *VA, ha_cgetype *B, ha_cgetype *X, PetscInt *IRNC, PetscInt *JCNC, PetscScalar *VAC,int *IRNB,int *JCNB,PetscScalar *VALUESB,ha_cgetype *VECBIVI,int *bivinzrow0,int *bivinzcol0);//, forint *IRNV, forint *JCNV, ha_cgetype *VAV
extern void spec48m_msol_(int *INSIZE,int *IRN, int *JCN, ha_cgetype *VA, ha_cgetype *B, ha_cgetype *X, PetscInt *IRNC, PetscInt *JCNC, PetscScalar *VAC,int *IRNB,int *JCNB,PetscScalar *VALUESB,ha_cgetype *VECBIVI,int *bivinzrow0,int *bivinzcol0);//, forint *IRNV, forint *JCNV, ha_cgetype *VAV
extern void spec48_esol_(int *INSIZE,int *IRN, ha_cgetype *VA,int *KEEP, ha_cgetype *B, ha_cgetype *X);
extern void spec48m_esol_(int *INSIZE,int *IRN, ha_cgetype *VA,int *KEEP, ha_cgetype *B, ha_cgetype *X);
extern void spec48_rpesol_(int *INSIZE,int *IRN, ha_cgetype *VA,int *KEEP, ha_cgetype *B, ha_cgetype *X,ha_cgetype *cntl,ha_cgetype *rinfo,ha_cgetype *error1,int *icntl,int *info,ha_cgetype *w,int *iw);
extern void spec48m_rpesol_(int *INSIZE,int *IRN, ha_cgetype *VA,int *KEEP, ha_cgetype *B, ha_cgetype *X,ha_cgetype *cntl,ha_cgetype *rinfo,ha_cgetype *error1,int *icntl,int *info,ha_cgetype *w,int *iw);
extern void spec48_single_(forint *indata,int *irn, int *jcn,ha_cgetype *b1, ha_cgetype *values,ha_cgetype *x1, int *neleperrow,int *ai1, MPI_Fint *fcomm);
extern void spec48_nomc66_(forint *indata, int *jcn,ha_cgetype *b1, ha_cgetype *values,ha_cgetype *x1, int *neleperrow, MPI_Fint *fcomm,forint *rowptrin, forint *colptrin);
extern void my_spar_add_(ha_cgetype *vecbivi, int *biviindx,int *nz1,ha_cgetype *vecbivi0,int *biviindx0,int *nz0,int *nz2);
extern void my_spar_addl_(ha_cgetype *vecbivi, long int *biviindx,int *nz1,ha_cgetype *vecbivi0,long int *biviindx0,int *nz0,int *nz2);
extern void my_spar_add3l_(ha_cgetype *vecbivi, long int *biviindx,long int *nz1,ha_cgetype *vecbivi0,long int *biviindx0,long int *nz0,long int *nz2);
extern void my_spar_add1_(ha_cgetype *vecbivi, int *biviindx,int *irn, int *jcn,int *nz1,ha_cgetype *vecbivi0,int *biviindx0,int *nz0,int *nz2,int *ncol);
extern void my_spar_add1l_(ha_cgetype *vecbivi, long int *biviindx,int *irn, int *jcn,int *nz1,ha_cgetype *vecbivi0,long int *biviindx0,int *nz0,int *nz2,int *ncol);
extern void my_spar_add4l_(ha_cgetype *vecbivi, long int *biviindx,int *irn, int *jcn,long int *nz1,ha_cgetype *vecbivi0,long int *biviindx0,long int *nz0,long int *nz2,int *ncol);
extern void my_spar_add2_(ha_cgetype *vecbivi, int *biviindx,int *irn, int *jcn,int *nz1,ha_cgetype *vecbivi0,int *biviindx0,int *nz0,int *nz2,int *ncol,ha_cgetype *vecbivi2,int *irn2, int *jcn2,int *j2,ha_cgetype *cntl3);
extern void my_spar_comp_(int *biviindx,int *nz1,int *biviindx0,int *nz0,int *nz2);
extern void my_spar_compl_(long int *biviindx,long int *nz1,long int *biviindx0,long int *nz0,long int *nz2);
extern void prep48_alu_(int *INSIZE,int *IRN,int* JCN,ha_cgetype *VA);
extern void prep48_alu1_(int *INSIZE,int *IRN,int* JCN,ha_cgetype *VA,ha_cgetype *W,int *IW,int *KEEP);
extern void prep48_msol_(int *INSIZE,int *IRN, int *JCN, ha_cgetype *VA, PetscInt *IRNC, PetscInt *JCNC, PetscScalar *VAC,int *IRNB,int *JCNB,PetscScalar *VALUESB,ha_cgetype *VECBIVI,long int *bivinzrow0,int *bivinzcol0);//, forint *IRNV, forint *JCNV, ha_cgetype *VAV
extern void prep48m_msol_(int *INSIZE,int *IRN, int *JCN, ha_cgetype *VA, PetscInt *IRNC, PetscInt *JCNC, PetscScalar *VAC,int *IRNB,int *JCNB,PetscScalar *VALUESB,ha_cgetype *VECBIVI,long int *bivinzrow0,int *bivinzcol0,int* jcnb1,ha_cgetype *sol,ha_cgetype *b,ha_cgetype *w,int* iw,int* keep);
extern void spar_mulmin_(ha_cgetype* sol,int* nrow,int* nz,int* irn,int* jcn,ha_cgetype* va,ha_cgetype* res);
extern void spar_muladd_(ha_cgetype* sol,int* nrow,int* nz,int* irn,int* jcn,ha_cgetype* va,ha_cgetype* res);
extern void spar_mulnoadd_(ha_cgetype* sol,int* nrow,int* nz,int* irn,int* jcn,ha_cgetype* va,ha_cgetype* res);
extern void spar_vbiviadd_(ha_cgetype* sol,int* bvcol,long int* bvrow,long int* bvsize,int* nrow,int *ncol,int* nz,int* irn,int* jcn,ha_cgetype* va,ha_cgetype* res);
extern void patio_mat_(int* insizeda,int* IRN,int* JCN,ha_cgetype* VBIVI,int* IRN1A,int* JCN1A);

int HaNewMatVal(char *fname, char *commsyntax,ha_cgeset *ha_set,uvadd nset, ha_cgesetele *ha_setele, hcge_cof *ha_cof,uvadd ncof,hcge_cof *ha_var,uvadd nvar, ha_cgevar *ha_cofvar,uvadd ncofvar,uvadd ncofele,ha_cgeexovar *ha_cgeshock,uvadd ndblock,uvadd alltimeset,uvadd allregset,uvadd *ha_eqadd,uvadd *counteq,uvadd nintraeq,Mat A,Mat B) {
  FILE * filehandle;
  char tline[TABREADLINE],line[TABREADLINE],line1[TABREADLINE],leftline[TABREADLINE],linecopy[TABREADLINE];//,set1[NAMESIZE],set2[NAMESIZE];
  char vname[TABREADLINE],sumsyntax[NAMESIZE],lintmp[TABREADLINE];//,*p1=NULL;
  char *readitem=NULL,*p=NULL,*p1=NULL;//,*p2=NULL,*varpnts;
  ha_cgetype vval;
  PetscScalar    vvalpetsc;
  PetscInt Iindx,Jindx=0,Istart1,Iend1,matrow,rowindx;//,jfix=0
  PetscErrorCode ierr;
  ha_cgetype zerodivide=0;
  PetscMPIInt  mpisize1;
  bool isinproc;
  ierr = MatGetOwnershipRange(A,&Istart1,&Iend1);
  MPI_Comm_size(PETSC_COMM_WORLD,&mpisize1);
  CHKERRQ(ierr);
  uvdim fdim,np,dcount,fdimlin=0,i4,sup,supset[MAXSUPSET];
  int totalsum,sumcount=1,sumcount1=0,lvar,lvar1,lvar2,lvar3,lvar4;
  uvadd lj,l1,i1=0,sumbegadd,dcountdim1[4*MAXVARDIM],dcountdim2[4*MAXVARDIM],dcountdim3[4*MAXVARDIM],nloops,nloopslin,nloopsfac,li3,nsumele,nsumele1,l2,eqindx=0,ltime;//,sizelinvars,totlinvars,templinvars
  int sumindx,npow,npar,nmul,nplu,ndiv,nmin,ha_calvarsize=0,nlinvars,leadlag,varindx1,varindx2;
  uvadd j,sj,l,i3,i5,i,arsetdim=0,ha_calvardim=0;
  ha_cgesetindx *arSet1=NULL;
  hcge_calvars *ha_calvar1= NULL;

  filehandle = fopen(fname,"r");
  matrow=0;//nintraeq;

  while (ha_cgertabl1(commsyntax,filehandle,line,ha_cofvar,ha_cof,ncof,&zerodivide,TABREADLINE)) {
    if (strstr(line,"(default")==NULL) {
      ha_cgefrstr1(line, commsyntax, "");
      ha_cgefrstr1(line, "(linear)", "");
      while (ha_cgefrstr(line,"  ", " "));
      while (ha_cgefrchr(line, '[', '('));
      while (ha_cgefrchr(line, ']', ')'));
      while (ha_cgefrchr(line, '{', '('));
      while (ha_cgefrchr(line, '}', ')'));
      strcpy(linecopy,line);
      fdim=ha_cgenfind(line, "(all,");
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
        i=ha_cgerevfind(line, "(all,");
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
      while (ha_cgefrstr(readitem," ", ""));
      while (ha_cgerecovar(readitem)==1);
      hnew_intrpl(readitem);
      strcpy(tline,readitem);
      strcpy(line1,readitem);
      npow=ha_cgenchf(readitem, '^');
      nmul=ha_cgenchf(readitem, '*');
      ndiv=ha_cgenchf(readitem, '/');
      nmul=nmul+ndiv;
      nplu=ha_cgenchf(readitem, '+');
      nmin=ha_cgenchf(readitem, '-');
      nplu=nplu+nmin;
      npar=ha_cgenchf(readitem, '(');

      strcpy(line,line1);
      readitem=line;
      np=ha_cgenfind(readitem,"p_");
      ha_calvardim=2*(npow+nmul+nplu+npar+1);
      hcge_calvars *ha_calvar= (hcge_calvars *) calloc (ha_calvardim,sizeof(hcge_calvars));
      hcge_linvars *LinVars= (hcge_linvars *) calloc (np,sizeof(hcge_linvars));
      i3=0;
      lvar=0;
      for (i=0; i<np; i++) {
        varindx2=0;
        while(-1<0) {
          varindx1=ha_cgefind(readitem+varindx2,"p_");
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
            if (strcmp(ha_var[l].cofname,vname)==0) {
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
          switch (ha_var[l].size) {
          case 0:
            break;
          case 1:
            p = strtok(tline,"{");
            p = strtok(NULL,"}");
            leadlag=0;
            hnew_arset(p,&leadlag);
            strcpy(LinVars[i3].dimnames[0],p);
            LinVars[i3].dimleadlag[0]=leadlag;
            strcpy(lintmp,"(all,");
            strcat(lintmp,p);
            strcat(lintmp,",");
            l1=ha_cgefind(linecopy,lintmp);
            if (l1>-1) {
              p1=&linecopy[0]+l1;
              strncpy(LinVars[i3].dimsetnames[0],p1+strlen(lintmp),strchr(p1,')')-p1-strlen(lintmp));
            }
            else {
              strcpy(lintmp,"sum(");
              strcat(lintmp,p);
              strcat(lintmp,",");
              lvar1=ha_cgenfind(linecopy,lintmp);
              lvar3=ha_cgefind(linecopy,lintmp);
              if (lvar1>1) for(lvar2=0; lvar2<lvar1; lvar2++) {
                  lvar4=ha_cgefind(&linecopy[lvar3+4],lintmp);
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
            for (i4=0; i4<ha_var[l].size-1; i4++) {
              p = strtok(NULL,",");
              leadlag=0;
              hnew_arset(p,&leadlag);
              strcpy(LinVars[i3].dimnames[i4],p);
              LinVars[i3].dimleadlag[i4]=leadlag;
              strcpy(lintmp,"(all,");
              strcat(lintmp,p);
              strcat(lintmp,",");
              l1=ha_cgefind(linecopy,lintmp);
              if (l1>-1) {
                p1=&linecopy[0]+l1;
                strncpy(LinVars[i3].dimsetnames[i4],p1+strlen(lintmp),strchr(p1,')')-p1-strlen(lintmp));
              }
              else {
                strcpy(lintmp,"sum(");
                strcat(lintmp,p);
                strcat(lintmp,",");
                lvar1=ha_cgenfind(linecopy,lintmp);
                lvar3=ha_cgefind(linecopy,lintmp);
                if (lvar1>1) for(lvar2=0; lvar2<lvar1; lvar2++) {
                    lvar4=ha_cgefind(&linecopy[lvar3+4],lintmp);
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
            hnew_arset(p,&leadlag);
            strcpy(LinVars[i3].dimnames[i4],p);
            LinVars[i3].dimleadlag[i4]=leadlag;
            strcpy(lintmp,"(all,");
            strcat(lintmp,p);
            strcat(lintmp,",");
            l1=ha_cgefind(linecopy,lintmp);
            if (l1>-1) {
              p1=&linecopy[0]+l1;
              strncpy(LinVars[i3].dimsetnames[i4],p1+strlen(lintmp),strchr(p1,')')-p1-strlen(lintmp));
            }
            else {
              strcpy(lintmp,"sum(");
              strcat(lintmp,p);
              strcat(lintmp,",");
              lvar1=ha_cgenfind(linecopy,lintmp);
              lvar3=ha_cgefind(linecopy,lintmp);
              if (lvar1>1) for(lvar2=0; lvar2<lvar1; lvar2++) {
                  lvar4=ha_cgefind(&linecopy[lvar3+4],lintmp);
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
            if (strcmp(ha_var[l].cofname,vname)==0) {
              LinVars[i3].LinVarIndx=l;
              break;
            }
          }
          i3++;
        }
      }
      nlinvars=i3;
      i3=0;
      for (i=0; i<nlinvars; i++) i3+=ha_var[LinVars[i].LinVarIndx].size;
      arsetdim=fdim+i3;
      ha_cgesetindx *arSet= (ha_cgesetindx *) calloc (arsetdim,sizeof(ha_cgesetindx));

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
          strcpy(arSet[i].arIndx,readitem);
          readitem = strtok(NULL,")");
          for (i4=0; i4<nset; i4++) if(strcmp(readitem,ha_set[i4].setname)==0) {
              arSet[i].setid=i4;
              break;
            }
          nloops=nloops*ha_set[arSet[i].setid].size;
        }
        dcountdim1[fdim-1]=1;
        for (dcount=fdim-2; dcount>-1; dcount--) {
          dcountdim1[dcount]=ha_set[arSet[dcount+1].setid].size*dcountdim1[dcount+1];
        }
      }
      isinproc=false;
      if(mpisize1>1) {
        for (lj=0; lj<nloops; lj++) {
          Jindx=ha_eqadd[matrow+lj];
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
        totalsum=hcge_nsum(readitem,sumsyntax);
        hcge_sumcof *sum_cof= (hcge_sumcof *) calloc (totalsum*nlinvars+1,sizeof(hcge_sumcof));
        sumcount=0;
        while (hcge_dsum(readitem,sumsyntax,sum_cof,arSet,ha_set,nset,fdim+1,sumcount)==1) {
          sumcount++;
        }
        totalsum=sumcount;
        li3=0;
        for (i4=0; i4<totalsum; i4++) {
          i1=1;
          for(sj=0; sj<sum_cof[i4].size; sj++) {
            i1=i1*ha_set[sum_cof[i4].setid[sj]].size;
          }
          sum_cof[i4].begadd=li3;
          li3=li3+i1;
        }
        nsumele=li3;
        ha_cgesumele *ha_sumele= (ha_cgesumele *) calloc (nsumele*nlinvars,sizeof(ha_cgesumele));

        strcpy(line,line1);
        readitem=line;
        sumcount=0;
        while (hlin_dsum(readitem,sumsyntax,sum_cof,arSet,ha_set,nset,fdim+1,sumcount)==1) {
          sumcount++;
        }
        totalsum=sumcount;
        li3=0;
        for (i4=0; i4<totalsum; i4++) {
          i1=1;
          for(sj=0; sj<sum_cof[i4].size; sj++) {
            i1=i1*ha_set[sum_cof[i4].setid[sj]].size;
          }
          sum_cof[i4].begadd=li3;
          li3=li3+i1;
        }
        nsumele=li3;
        for (i3=0; i3<totalsum; i3++) {
          i1=1;
          sum_cof[i3].antidims[sum_cof[i3].size-1]=1;
          for(sj=sum_cof[i3].size-2; sj>-1; sj--) {
            sum_cof[i3].antidims[sj]=sum_cof[i3].antidims[sj+1]*ha_set[sum_cof[i3].setid[sj+1]].size;
          }
        }
        sumcount=0;
        strcpy(line,line1);
        readitem=line;
        sumindx=0;
        while (hlin_calsum(readitem,sumsyntax,ha_set,nset,ha_setele,ha_cofvar,ncofvar,ncofele,ha_cof,ncof,ha_var,nvar,sum_cof,totalsum,ha_sumele,nsumele,ha_calvar,arSet,fdim+1,&sumindx,sumcount,zerodivide)==1) {
          sumcount++;
        }
        strcpy(line1,readitem);
        sumbegadd=nsumele;
        sumcount1=sumcount;
        for (i=0; i<nlinvars; i++) {
          Jindx=ha_eqadd[matrow];
          if(Jindx>=Iend1)continue;
          i3=0;
          nloopslin=nloops;
          if (fdim==0) {
            for (l2=0; l2<ha_var[LinVars[i].LinVarIndx].size; l2++) {
              strcpy(arSet[i3].arIndx,LinVars[i].dimnames[l2]);
              for (i4=0; i4<nset; i4++) if(strcmp(ha_set[i4].setname,LinVars[i].dimsetnames[l2])==0) {
                  break;
                }
              arSet[i3].setid=i4;
              nloopslin=nloopslin*ha_set[arSet[i3].setid].size;
              i3++;
            }
          }
          else for (l2=0; l2<ha_var[LinVars[i].LinVarIndx].size; l2++) {
              for (sj=0; sj<fdim; sj++) {
                if(strcmp(arSet[sj].arIndx,LinVars[i].dimnames[l2])==0) {
                  break;
                }
                else if (sj==fdim-1) {
                  strcpy(arSet[fdim+i3].arIndx,LinVars[i].dimnames[l2]);
                  for (i4=0; i4<nset; i4++) if(strcmp(ha_set[i4].setname,LinVars[i].dimsetnames[l2])==0) {
                      break;
                    }
                  arSet[fdim+i3].setid=i4;
                  nloopslin=nloopslin*ha_set[arSet[fdim+i3].setid].size;
                  i3++;
                }
              }
            }
          nloopsfac=(uvadd)nloopslin/nloops;
          Jindx=ha_eqadd[matrow+(uvadd)(nloopslin-1)/nloopsfac];
          if(Jindx<Istart1)continue;
          
          fdimlin=fdim+i3;
          if (i3>0) {
            dcountdim2[fdimlin-1]=1;
            for (dcount=fdimlin-2; dcount>-1; dcount--) {
              dcountdim2[dcount]=ha_set[arSet[dcount+1].setid].size*dcountdim2[dcount+1];
            }
          }
          else for (dcount=0; dcount<fdim; dcount++) {
              dcountdim2[dcount]=dcountdim1[dcount];
            }

          strcpy(leftline,line1);
          hcge_rlinzero(leftline,i);
          hcge_repllin(leftline,0);
          strcpy(sumsyntax,"sum(");
          hlin_replsum(leftline,sumsyntax,i,LinVars,ha_var);
          hlin_rlinone(leftline,LinVars,i,ha_var);
          strcpy(line,leftline);
          readitem=line;
          sumcount1=sumcount;
          while (hcge_dsum(readitem,sumsyntax,sum_cof,arSet,ha_set,nset,fdimlin+1,sumcount1)==1) {
            sumcount1++;
          }
          totalsum=sumcount1;
          for (i4=sumcount; i4<totalsum; i4++) {
            i1=1;
            for(sj=0; sj<sum_cof[i4].size; sj++) {
              i1=i1*ha_set[sum_cof[i4].setid[sj]].size;
            }
            sum_cof[i4].begadd=sumbegadd;
            sumbegadd=sumbegadd+i1;
          }
          nsumele1=sumbegadd+nsumele;
          for (i3=sumcount; i3<totalsum; i3++) {
            i1=1;
            sum_cof[i3].antidims[sum_cof[i3].size-1]=1;
            for(sj=sum_cof[i3].size-2; sj>-1; sj--) {
              sum_cof[i3].antidims[sj]=sum_cof[i3].antidims[sj+1]*ha_set[sum_cof[i3].setid[sj+1]].size;
            }
          }
          strcpy(line,leftline);
          readitem=line;
          while (hnew_calsum(readitem,sumsyntax,ha_set,nset,ha_setele,ha_cofvar,ncofvar,ncofele,ha_cof,ncof,ha_var,nvar,sum_cof,totalsum,ha_sumele,nsumele1,ha_calvar,arSet,fdimlin+1,&sumindx,sumcount,zerodivide)==1) {
            sumcount++;
          }
          ha_calvarsize=0;
          ha_newfparse(readitem,ha_set,ha_cof,ncof,ha_var,nvar,ncofele,sum_cof,totalsum,ha_calvar,&ha_calvarsize,arSet,fdimlin);
          for (dcount=0; dcount<ha_var[LinVars[i].LinVarIndx].size; dcount++) {
            for (i4=0; i4<fdimlin; i4++) {
              if (strcmp(LinVars[i].dimnames[dcount],arSet[i4].arIndx)==0) {
                dcountdim3[dcount]=i4;
                break;
              }
            }
          }
          for(dcount=0; dcount<MAXSUPSET; dcount++)supset[dcount]=0;
          for (dcount=0; dcount<ha_var[LinVars[i].LinVarIndx].size; dcount++) {
            if(ha_set[ha_var[LinVars[i].LinVarIndx].setid[dcount]].size!=ha_set[arSet[dcountdim3[dcount]].setid].size) {
              for(sup=1; sup<MAXSUPSET; sup++)if(ha_var[LinVars[i].LinVarIndx].setid[dcount]==ha_set[arSet[dcountdim3[dcount]].setid].subsetid[sup]) {
                  supset[dcount]=sup;
                  break;
                }
            }
            else supset[dcount]=0;
          }
        #pragma omp parallel private(lj,Jindx,i3,sj,i5,l2,dcount,l1,li3,Iindx,ierr,arSet1,ha_calvar1,vval) shared(ha_cofvar,arSet)
        {
        if(omp_get_thread_num()!=0){
          arSet1=realloc(arSet1,arsetdim*sizeof(ha_cgesetindx));
          memcpy(arSet1,arSet,arsetdim*sizeof(ha_cgesetindx));
          ha_calvar1=realloc(ha_calvar1,ha_calvardim*sizeof(hcge_calvars));
          memcpy(ha_calvar1,ha_calvar,ha_calvardim*sizeof(hcge_calvars));
        }else{
          ha_calvar1=ha_calvar;
          arSet1=arSet;
        }
          ha_cgetype *value= (ha_cgetype *) calloc (nloopsfac,sizeof(ha_cgetype));
          PetscInt *jcn= (PetscInt *) calloc (nloopsfac,sizeof(PetscInt));
          ha_cgetype *valueb= (ha_cgetype *) calloc (nloopsfac,sizeof(ha_cgetype));
          PetscInt *jcnb= (PetscInt *) calloc (nloopsfac,sizeof(PetscInt));
        #pragma omp for
          for (i5=0; i5<nloops; i5++) {
            Jindx=ha_eqadd[matrow+i5];//Jindx=ha_eqadd[matrow+(uvadd)lj/nloopsfac];
            if(Jindx>=Istart1&&Jindx<Iend1) {
              i3=0;
              sj=0;
              for (lj=i5*nloopsfac; lj<(i5+1)*nloopsfac; lj++) {
              l2=lj;
              for (dcount=0; dcount<fdimlin; dcount++) {
                l1=(uvadd) l2/dcountdim2[dcount];
                arSet1[dcount].indx=l1;
                l2=l2-l1*dcountdim2[dcount];
              }
              li3=0;
              for (dcount=0; dcount<ha_var[LinVars[i].LinVarIndx].size; dcount++) {
                if(supset[dcount]==0) {
                  li3=li3+(arSet1[dcountdim3[dcount]].indx+LinVars[i].dimleadlag[dcount])*ha_var[LinVars[i].LinVarIndx].antidims[dcount];
                }
                else {
                  li3=li3+(ha_setele[ha_set[arSet1[dcountdim3[dcount]].setid].begadd+arSet1[dcountdim3[dcount]].indx].setsh[supset[dcount]]+LinVars[i].dimleadlag[dcount])*ha_var[LinVars[i].LinVarIndx].antidims[dcount];
                }
              }
              vval=ha_newfpcal(ha_cofvar,ha_set,ha_setele,ha_sumele,ha_calvar1,ha_calvarsize,arSet1,fdimlin,zerodivide);
              Iindx=ha_cgeshock[ha_var[LinVars[i].LinVarIndx].begadd+li3].ExoIndx;
              if (!ha_cgeshock[ha_var[LinVars[i].LinVarIndx].begadd+li3].ShockId&&vval!=0) {
                value[i3]=vval;
                jcn[i3]=Iindx;
                i3++;
              }
              if (ha_cgeshock[ha_var[LinVars[i].LinVarIndx].begadd+li3].ShockId&&vval!=0) {
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
          free(ha_calvar1);
          ha_calvar1=NULL;
        }else{
          ha_calvar1=NULL;
          arSet1=NULL;
        }
        }
        }
        free(sum_cof);
        free(ha_sumele);
      }
      matrow+=nloops;
      eqindx++;

      free(ha_calvar);
      free(LinVars);
      free(arSet);
    }
  }
  fclose(filehandle);
  return 1;
}

int hlin_dsum(char *formulain, char *commsyntax, hcge_sumcof *sum_cof,ha_cgesetindx *arSet,ha_cgeset *ha_set,uvdim nset,uvdim fdim,int j) {
  char *readitem,*p,*p1,*p2,interchar2[NAMESIZE],argu[TABREADLINE],tempname[NAMESIZE];//,line5[TABREADLINE]
  char interchar[NAMESIZE],interchar1[NAMESIZE],line[TABREADLINE],line1[TABREADLINE],line2[TABREADLINE],line3[TABREADLINE],line4[TABREADLINE];
  uvdim l,l1,l2,l3,l4,l5,l6,l7;
  int ncur=0,ncuri,i=0,k=0,k1=0,length;
  length=strlen(formulain);
  readitem=formulain;
  while (i<length) {
    k=ha_cgefind(readitem,commsyntax);
    if (k==-1) {
      return 0;
    }
    if (k==0) {
      readitem=formulain+i+k;
      strcpy(line,readitem);
      strcpy(line1,readitem);
      ha_cgecutsum(line);
      k1=ha_cgefind(line+4,commsyntax);
      if (k1!=-1) {
        i=i+k+4;
        readitem=formulain+i;
      }
      else {
        if(ha_cgefind(line,",p_")>-1||ha_cgefind(line,"*p_")>-1||ha_cgefind(line,"+p_")>-1||ha_cgefind(line,"-p_")>-1||ha_cgefind(line,"(p_")>-1) {
          i=i+k+4;
          readitem=formulain+i;
        }
        else {
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

          ncur=ha_cgenfind(line2, "{");
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
            l=ha_cgenfind(argu, ",");
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
                  for (l5=0; l5<fdim-1; l5++) if(strcmp(p,arSet[l5].arIndx)==0) {
                      sum_cof[j].setid[l3]=l5;
                      l6++;
                    }
                  if (l6==0) {
                    interchar1[0]='\0';
                    strcat(interchar1,"sum(");
                    strcat(interchar1,p);
                    strcpy(line3,formulain);
                    line3[readitem-formulain]='\0';
                    l7=ha_cgerevfind(line3,interchar1);
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
                    for (l5=0; l5<fdim-1; l5++) if(strcmp(p,arSet[l5].arIndx)==0) {
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
                      l7=ha_cgerevfind(line3,interchar1);
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
          }
          else {
            if (interchar[strlen(interchar)-1]=='{') {
              interchar[strlen(interchar)-1]='\0';
            }
            else {
              strcat(interchar,"}");
            }
          }
          while(ha_cgefrstr(formulain,line4,interchar));
          sum_cof[j].size=l3;
          return 1;
        }
      }

    }
    else if (formulain[i+k-1]=='+'||formulain[i+k-1]=='-'||formulain[i+k-1]=='*'||formulain[i+k-1]=='/'||formulain[i+k-1]=='^'||formulain[i+k-1]=='('||formulain[i+k-1]==',') {
      readitem=formulain+i+k;
      strcpy(line,readitem);
      strcpy(line1,readitem);
      ha_cgecutsum(line);
      k1=ha_cgefind(line+4,commsyntax);
      if (k1!=-1) {
        i=i+k+4;
        readitem=formulain+i;
      }
      else {
        if(ha_cgefind(line,",p_")>-1||ha_cgefind(line,"*p_")>-1||ha_cgefind(line,"+p_")>-1||ha_cgefind(line,"-p_")>-1||ha_cgefind(line,"(p_")>-1) {
          i=i+k+4;
          readitem=formulain+i;
        }
        else {
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

          ncur=ha_cgenfind(line2, "{");
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
            l=ha_cgenfind(argu, ",");
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
                  for (l5=0; l5<fdim-1; l5++) if(strcmp(p,arSet[l5].arIndx)==0) {
                      sum_cof[j].setid[l3]=l5;
                      l6++;
                    }
                  if (l6==0) {
                    interchar1[0]='\0';
                    strcat(interchar1,"sum(");
                    strcat(interchar1,p);
                    strcpy(line3,formulain);
                    line3[readitem-formulain]='\0';
                    l7=ha_cgerevfind(line3,interchar1);
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
                    for (l5=0; l5<fdim-1; l5++) if(strcmp(p,arSet[l5].arIndx)==0) {
                        sum_cof[j].setid[l3]=arSet[l5].setid;//l5;
                        l6++;
                      }
                    if (l6==0) {
                      interchar1[0]='\0';
                      strcat(interchar1,"sum(");
                      strcat(interchar1,p);
                      strcpy(line3,formulain);
                      line3[readitem-formulain]='\0';
                      l7=ha_cgerevfind(line3,interchar1);
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
          }
          else {
            if (interchar[strlen(interchar)-1]=='{') {
              interchar[strlen(interchar)-1]='\0';
            }
            else {
              strcat(interchar,"}");
            }
          }
          while(ha_cgefrstr(formulain,line4,interchar));
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

int hlin_calsum(char *formulain, char *commsyntax,ha_cgeset *ha_set,uvdim nset, ha_cgesetele *ha_setele,ha_cgevar *ha_cofvar,uvadd ncofvar,uvadd ncofele, hcge_cof *ha_cof,uvadd ncof,hcge_cof *ha_var,uvadd nvar,hcge_sumcof *sum_cof,int totalsum,ha_cgesumele *ha_sumele,uvadd nsumele,hcge_calvars *ha_calvar,ha_cgesetindx *arSet1,uvdim fdim,int *sumindx,int j, ha_cgetype zerodivide) {
  char *readitem,*p;//,*p1,interchar2[NAMESIZE],line5[TABREADLINE];
  char interchar[NAMESIZE],line[TABREADLINE],line1[TABREADLINE],line2[TABREADLINE];//,line3[TABREADLINE],line4[TABREADLINE];//,interchar1[NAMESIZE]
  int i=0,k=0,k1=0,length;//,simpl=0;//,ncur=0,ncuri,l3,l4,l5,l6,l7
  uvdim dcount,setsh,fdimsumcof,l;
  uvadd dcountdim1[4*MAXVARDIM],nloops,l1,l2,l3;
  int ha_calvarsize;
  ha_cgetype vval;
  ha_cgesetindx *arSet2=NULL;
  hcge_calvars *ha_calvar1= NULL;
  uvadd arsetsize;
  length=strlen(formulain);
  readitem=formulain;
  while (i<length) {
    k=ha_cgefind(readitem,commsyntax);
    if (k==-1) {
      return 0;
    }
    if (k==0) {
      readitem=formulain+i+k;
      strcpy(line,readitem);
      ha_cgecutsum(line);
      k1=ha_cgefind(line+4,commsyntax);
      if (k1!=-1) {
        i=i+k+4;
        readitem=formulain+i;
      }
      else {
        if(ha_cgefind(line,",p_")>-1||ha_cgefind(line,"*p_")>-1||ha_cgefind(line,"+p_")>-1||ha_cgefind(line,"-p_")>-1||ha_cgefind(line,"(p_")>-1) {
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
          ha_cgesetindx *arSet= (ha_cgesetindx *) calloc (arsetsize,sizeof(ha_cgesetindx));
          for (l=0; l<sum_cof[j].size; l++) {
            arSet[l].setid=sum_cof[j].setid[l];
            strcpy(arSet[l].arIndx,sum_cof[j].dimnames[l]);
          }
          nloops=1;
          for (l=0; l<sum_cof[j].size; l++) {
            nloops=nloops*ha_set[arSet[l].setid].size;
            dcount=sum_cof[j].size-l;
            if(dcount==sum_cof[j].size) {
              dcountdim1[dcount-1]=1;
            }
            else {
              dcountdim1[dcount-1]=dcountdim1[dcount]*ha_set[arSet[dcount].setid].size;
            }
          }
          arSet[sum_cof[j].size].setid=sum_cof[j].sumsetid;
          strcpy(arSet[sum_cof[j].size].arIndx,sum_cof[j].sumindx);
          fdimsumcof=sum_cof[j].size+1;
          ha_calvarsize=0;
          ha_newfparse(p,ha_set,ha_cof,ncof,ha_var,nvar,ncofele,sum_cof,totalsum,ha_calvar,&ha_calvarsize,arSet,fdimsumcof);
        #pragma omp parallel private(l3,l1,l2,dcount,setsh,vval,arSet2,ha_calvar1) shared(ha_cofvar,arSet,ha_sumele)
        {
        if(omp_get_thread_num()!=0){
          arSet2=realloc(arSet2,arsetsize*sizeof(ha_cgesetindx));
          memcpy(arSet2,arSet,arsetsize*sizeof(ha_cgesetindx));
          ha_calvar1=realloc(ha_calvar1,ha_calvarsize*sizeof(hcge_calvars));
          memcpy(ha_calvar1,ha_calvar,ha_calvarsize*sizeof(hcge_calvars));
        }else{
          ha_calvar1=ha_calvar;
          arSet2=arSet;
        }
        #pragma omp for
          for (l3=0; l3<nloops; l3++) {
            l2=l3;
            for (dcount=0; dcount<sum_cof[j].size; dcount++) {
              setsh=(uvadd) l2/dcountdim1[dcount];
              arSet2[dcount].indx=setsh;
              l2=l2-setsh*dcountdim1[dcount];
            }
            vval=0;
            for (l1=0; l1<ha_set[sum_cof[j].sumsetid].size; l1++) {
              arSet2[sum_cof[j].size].indx=l1;
              vval+=ha_newfpcal(ha_cofvar,ha_set,ha_setele,ha_sumele,ha_calvar1,ha_calvarsize,arSet2,fdimsumcof,zerodivide);
            }
            ha_sumele[*sumindx+l3].varval=vval;
          }
        if(omp_get_thread_num()!=0){
          free(arSet2);
          arSet2=NULL;
          free(ha_calvar1);
          ha_calvar1=NULL;
        }else{
          ha_calvar1=NULL;
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
          while(ha_cgefrstr(formulain,line1,interchar)!=NULL);
          free(arSet);
          return 1;
        }
      }

    }
    else if (formulain[i+k-1]=='+'||formulain[i+k-1]=='-'||formulain[i+k-1]=='*'||formulain[i+k-1]=='/'||formulain[i+k-1]=='^'||formulain[i+k-1]=='('||formulain[i+k-1]==',') {
      readitem=formulain+i+k;
      strcpy(line,readitem);
      ha_cgecutsum(line);
      k1=ha_cgefind(line+4,commsyntax);
      if (k1!=-1) {
        i=i+k+4;
        readitem=formulain+i;
      }
      else {
        if(ha_cgefind(line,",p_")>-1||ha_cgefind(line,"*p_")>-1||ha_cgefind(line,"+p_")>-1||ha_cgefind(line,"-p_")>-1||ha_cgefind(line,"(p_")>-1) {
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
          ha_cgesetindx *arSet= (ha_cgesetindx *) calloc (arsetsize,sizeof(ha_cgesetindx));
          for (l=0; l<sum_cof[j].size; l++) {
            arSet[l].setid=sum_cof[j].setid[l];
            strcpy(arSet[l].arIndx,sum_cof[j].dimnames[l]);
          }
          nloops=1;
          for (l=0; l<sum_cof[j].size; l++) {
            nloops=nloops*ha_set[arSet[l].setid].size;//sum_cof[j].dims[l];
            dcount=sum_cof[j].size-l;
            if(dcount==sum_cof[j].size) {
              dcountdim1[dcount-1]=1;
            }
            else {
              dcountdim1[dcount-1]=dcountdim1[dcount]*ha_set[arSet[dcount].setid].size;
            }
          }
          arSet[sum_cof[j].size].setid=sum_cof[j].sumsetid;
          strcpy(arSet[sum_cof[j].size].arIndx,sum_cof[j].sumindx);
          fdimsumcof=sum_cof[j].size+1;
          ha_calvarsize=0;
          ha_newfparse(p,ha_set,ha_cof,ncof,ha_var,nvar,ncofele,sum_cof,totalsum,ha_calvar,&ha_calvarsize,arSet,fdimsumcof);
        #pragma omp parallel private(l3,l1,l2,dcount,setsh,vval,arSet2,ha_calvar1) shared(ha_cofvar,arSet,ha_sumele)
        {
        if(omp_get_thread_num()!=0){
          arSet2=realloc(arSet2,arsetsize*sizeof(ha_cgesetindx));
          memcpy(arSet2,arSet,arsetsize*sizeof(ha_cgesetindx));
          ha_calvar1=realloc(ha_calvar1,ha_calvarsize*sizeof(hcge_calvars));
          memcpy(ha_calvar1,ha_calvar,ha_calvarsize*sizeof(hcge_calvars));
        }else{
          ha_calvar1=ha_calvar;
          arSet2=arSet;
        }
        #pragma omp for
          for (l3=0; l3<nloops; l3++) {
            l2=l3;
            for (dcount=0; dcount<sum_cof[j].size; dcount++) {
              setsh=(uvadd) l2/dcountdim1[dcount];
              arSet2[dcount].indx=setsh;
              l2=l2-setsh*dcountdim1[dcount];
            }
            vval=0;
            for (l1=0; l1<ha_set[sum_cof[j].sumsetid].size; l1++) {
              arSet2[sum_cof[j].size].indx=l1;
              vval+=ha_newfpcal(ha_cofvar,ha_set,ha_setele,ha_sumele,ha_calvar1,ha_calvarsize,arSet2,fdimsumcof,zerodivide);
            }
            ha_sumele[(uvadd)*sumindx+l3].varval=vval;//ha_sumele[*sumindx+l2].varval=vval;
          }
        if(omp_get_thread_num()!=0){
          free(arSet2);
          arSet2=NULL;
          free(ha_calvar1);
          ha_calvar1=NULL;
        }else{
          ha_calvar1=NULL;
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
          while(ha_cgefrstr(formulain,line1,interchar));
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

int hlin_replsum(char *formulain, char *commsyntax,int LinIndx, hcge_linvars *LinVars,hcge_cof *ha_var) {
  char *readitem,line[TABREADLINE],line1[TABREADLINE],line2[TABREADLINE],line3[TABREADLINE],argu[TABREADLINE],*p0,*p1;
  int i=0,k=0,k1=0,k2,k3,length,tsum,j,l;
  bool IsRemSum;
  length=strlen(formulain);
  readitem=formulain;
  while (i<length) {
    k=ha_cgefind(formulain+i,commsyntax);
    if (k==-1) {
      return 0;
    }
    if (k==0) {
      readitem=formulain+i+k;
      strcpy(line,readitem);
      ha_cgecutsum(line);
      strcpy(line3,line);
      if(ha_cgefind(line,",p_")>-1||ha_cgefind(line,"*p_")>-1||ha_cgefind(line,"+p_")>-1||ha_cgefind(line,"-p_")>-1||ha_cgefind(line,"(p_")>-1) {
        tsum=hcge_nsum(line,commsyntax);
        IsRemSum=false;
        for (j=tsum; j>1; j--) {
          k2=1;
          k1=0;
          while(k2<j) {
            k3=k1;
            k1=ha_cgefind(line+k1+3,commsyntax);
            k1=k3+k1+3;//printf("k1 %d line %s\n",k1,line+k1-1);
            if (line[k1-1]=='+'||line[k1-1]=='-'||line[k1-1]=='*'||line[k1-1]=='/'||line[k1-1]=='^'||line[k1-1]=='('||line[k1-1]==',') k2++;//printf("k1 %d line %s\n",k1,line+k1-1);}
          }
          strcpy(line1,&line[k1]);
          ha_cgecutsum(line1);
          p0=strchr(line1,'(');
          p0++;
          p1=strchr(line1,',');
          strncpy(argu,p0,p1-p0);
          argu[p1-p0]='\0';
          for(l=0; l<ha_var[LinVars[LinIndx].LinVarIndx].size; l++) {
            if(strcmp(argu,LinVars[LinIndx].dimnames[l])==0) {
              p1++;
              p1=strchr(p1,',');
              p1++;
              strncpy(line2,p1,strlen(line1)-(p1-line1)-1);
              line2[strlen(line1)-(p1-line1)-1]='\0';
              ha_cgefrstr(line,line1,line2);
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
        for(l=0; l<ha_var[LinVars[LinIndx].LinVarIndx].size; l++) {
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
          ha_cgefrstr(formulain,line3,line);
          i=i+k+strlen(line);
        }
        else i=i+k+strlen(line3);
      }
      else {
        ha_cgefrstr(formulain,line,"0");
        i=i+k+1;
      }
    }
    else if (formulain[i+k-1]=='+'||formulain[i+k-1]=='-'||formulain[i+k-1]=='*'||formulain[i+k-1]=='/'||formulain[i+k-1]=='^'||formulain[i+k-1]=='('||formulain[i+k-1]==',') {
      readitem=formulain+i+k;
      strcpy(line,readitem);
      ha_cgecutsum(line);
      strcpy(line3,line);
      if(ha_cgefind(line,",p_")>-1||ha_cgefind(line,"*p_")>-1||ha_cgefind(line,"+p_")>-1||ha_cgefind(line,"-p_")>-1||ha_cgefind(line,"(p_")>-1) {
        tsum=hcge_nsum(line,commsyntax);
        IsRemSum=false;
        for (j=tsum; j>1; j--) {
          k2=1;
          k1=0;
          while(k2<j) {
            k3=k1;
            k1=ha_cgefind(line+k1+3,commsyntax);
            k1=k3+k1+3;//printf("k1 %d line %s\n",k1,line+k1-1);
            if (line[k1-1]=='+'||line[k1-1]=='-'||line[k1-1]=='*'||line[k1-1]=='/'||line[k1-1]=='^'||line[k1-1]=='('||line[k1-1]==',') k2++;//printf("k1 %d line %s\n",k1,line+k1-1);}
          }
          strcpy(line1,&line[k1]);
          ha_cgecutsum(line1);
          p0=strchr(line1,'(');
          p0++;
          p1=strchr(line1,',');
          strncpy(argu,p0,p1-p0);
          argu[p1-p0]='\0';
          for(l=0; l<ha_var[LinVars[LinIndx].LinVarIndx].size; l++) {
            if(strcmp(argu,LinVars[LinIndx].dimnames[l])==0) {
              p1++;
              p1=strchr(p1,',');
              p1++;
              strncpy(line2,p1,strlen(line1)-(p1-line1)-1);
              line2[strlen(line1)-(p1-line1)-1]='\0';
              ha_cgefrstr(line,line1,line2);
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
        for(l=0; l<ha_var[LinVars[LinIndx].LinVarIndx].size; l++) {
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
          ha_cgefrstr(formulain,line3,line);
          i=i+k+strlen(line);
        }
        else i=i+k+strlen(line3);
      }
      else {
        ha_cgefrstr(formulain,line,"0");
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

int hlin_rlinone(char *formulain,hcge_linvars *LinVars,int linindx,hcge_cof *ha_var) {
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
    i1=ha_cgefind(formulain+i1,line);
    if (i1==-1) break;
    if (i1==0) {
      if(ha_var[LinVars[linindx].LinVarIndx].size==0) {
        i2=l-1;
      }
      else
        i2=ha_cgefind(formulain,"}");
      formulain[0]='1';
      memmove(formulain+1,formulain+i2+1,l-i2);
      i3=1;
    }
    else {
      if (formulain[i1-1]=='+'||formulain[i1-1]=='-'||formulain[i1-1]=='*'||formulain[i1-1]=='/'||formulain[i1-1]=='^'||formulain[i1-1]=='('||formulain[i1-1]==',') {
        if(ha_var[LinVars[linindx].LinVarIndx].size==0) {
          i2=strlen(line)-1;
        }
        else i2=ha_cgefind(formulain+i1,"}");
        formulain[i1]='1';
        memmove(formulain+i1+1,formulain+i1+i2+1,l-i2-i1);
        i3=1;
      }
      else i1++;
    }
  }
  return 1;
}

int NewMatvarRead(char *fname, char *commsyntax,ha_cgeset *ha_set,uvdim nset,ha_cgesetele *ha_setele,hcge_cof *ha_cof,uvadd ncof,hcge_cof *ha_var,uvadd nvar,ha_cgevar *ha_cofvar,uvadd ncofvar,uvadd ncofele,ha_cgeexovar *ha_cgeshock,bool *var_inter,hcge_cof *ha_eq,bool *ha_eqint,uvdim *eq_orderintra,uvdim *eq_orderreg,uvadd allregset,uvadd alltimeset,uvdim *orderintra,uvdim *orderreg) {
  FILE * filehandle;
  char tline[TABREADLINE],line[TABREADLINE],line1[TABREADLINE],linecopy[TABREADLINE];//,set1[NAMESIZE],set2[NAMESIZE];
  char vname[TABREADLINE],lintmp[TABREADLINE];//,*p1=NULL;
  char *readitem=NULL,*p=NULL,*p1=NULL;//,*p2=NULL,*varpnts;
  ha_cgetype zerodivide=0;
  uvdim fdim,np,i4;
  uvadd rowindx,j,l,l1,lj,dcountdim1[4*MAXVARDIM],dcountdim2[4*MAXVARDIM],nloops,nloopslin,li3,l2,matrow,eqindx=0,matsize,l01,j01;//,sizelinvars,totlinvars,templinvars
  int sj,i,i3,nlinvars,lvar,lvar1,lvar2,lvar3,lvar4,dcount,fdimlin=0,leadlag,varindx1,varindx2,tempint;
  ha_cgetype dimmat[4*MAXVARDIM];

  filehandle = fopen(fname,"r");
  matrow=0;

  while (ha_cgertabl1(commsyntax,filehandle,line,ha_cofvar,ha_cof,ncof,&zerodivide,TABREADLINE)) {
    if (strstr(line,"(default")==NULL) {
      ha_cgefrstr1(line, commsyntax, "");
      ha_cgefrstr1(line, "(linear)", "");
      while (ha_cgefrstr(line,"  ", " "));
      while (ha_cgefrchr(line, '[', '('));
      while (ha_cgefrchr(line, ']', ')'));
      while (ha_cgefrchr(line, '{', '('));
      while (ha_cgefrchr(line, '}', ')'));
      strcpy(linecopy,line);
      fdim=ha_cgenfind(line, "(all,");
      if (fdim==0) {
        readitem = strtok(line+1," ");
        strcpy(ha_eq[eqindx].cofname,readitem);
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
        strcpy(ha_eq[eqindx].cofname,readitem);
        strcpy(line,linecopy);
        i=ha_cgerevfind(line, "(all,");
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
      while (ha_cgefrstr(readitem," ", ""));
      while (ha_cgerecovar(readitem)==1);
      hnew_intrpl(readitem);
      strcpy(tline,readitem);
      strcpy(line1,readitem);
      np=ha_cgenfind(readitem,"p_");
      hcge_linvars *LinVars= (hcge_linvars *) calloc (np+1,sizeof(hcge_linvars));
      i3=0;
      lvar=0;
      for (i=0; i<np; i++) {
        varindx2=0;
        while(-1<0) {
          varindx1=ha_cgefind(readitem+varindx2,"p_");
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
            if (strcmp(ha_var[l].cofname,vname)==0) {
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
          switch (ha_var[l].size) {
          case 0:
            break;
          case 1:
            p = strtok(tline,"{");
            p = strtok(NULL,"}");
            leadlag=0;
            hnew_arset(p,&leadlag);
            if(leadlag!=0)var_inter[l]=true;//printf("var %s\n",ha_var[l].cofname);}
            strcpy(LinVars[i3].dimnames[0],p);
            strcpy(lintmp,"(all,");
            strcat(lintmp,p);
            strcat(lintmp,",");
            l1=ha_cgefind(linecopy,lintmp);
            if (l1>-1) {
              p1=&linecopy[0]+l1;
              strncpy(LinVars[i3].dimsetnames[0],p1+strlen(lintmp),strchr(p1,')')-p1-strlen(lintmp));
            }
            else {
              strcpy(lintmp,"sum(");
              strcat(lintmp,p);
              strcat(lintmp,",");
              lvar1=ha_cgenfind(linecopy,lintmp);
              lvar3=ha_cgefind(linecopy,lintmp);
              if (lvar1>1) for(lvar2=0; lvar2<lvar1; lvar2++) {
                  lvar4=ha_cgefind(&linecopy[lvar3+4],lintmp);
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
            for (i4=0; i4<ha_var[l].size-1; i4++) {
              p = strtok(NULL,",");
              leadlag=0;
              hnew_arset(p,&leadlag);
              if(leadlag!=0)var_inter[l]=true;//printf("var %s\n",ha_var[l].cofname);}
              strcpy(LinVars[i3].dimnames[i4],p);
              strcpy(lintmp,"(all,");
              strcat(lintmp,p);
              strcat(lintmp,",");
              l1=ha_cgefind(linecopy,lintmp);
              if (l1>-1) {
                p1=&linecopy[0]+l1;
                strncpy(LinVars[i3].dimsetnames[i4],p1+strlen(lintmp),strchr(p1,')')-p1-strlen(lintmp));
              }
              else {
                strcpy(lintmp,"sum(");
                strcat(lintmp,p);
                strcat(lintmp,",");
                lvar1=ha_cgenfind(linecopy,lintmp);
                lvar3=ha_cgefind(linecopy,lintmp);
                if (lvar1>1) for(lvar2=0; lvar2<lvar1; lvar2++) {
                    lvar4=ha_cgefind(&linecopy[lvar3+4],lintmp);
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
            hnew_arset(p,&leadlag);
            if(leadlag!=0)var_inter[l]=true;//printf("var %s\n",ha_var[l].cofname);}
            strcpy(LinVars[i3].dimnames[i4],p);
            strcpy(lintmp,"(all,");
            strcat(lintmp,p);
            strcat(lintmp,",");
            l1=ha_cgefind(linecopy,lintmp);
            if (l1>-1) {
              p1=&linecopy[0]+l1;
              strncpy(LinVars[i3].dimsetnames[i4],p1+strlen(lintmp),strchr(p1,')')-p1-strlen(lintmp));
            }
            else {
              strcpy(lintmp,"sum(");
              strcat(lintmp,p);
              strcat(lintmp,",");
              lvar1=ha_cgenfind(linecopy,lintmp);
              lvar3=ha_cgefind(linecopy,lintmp);
              if (lvar1>1) for(lvar2=0; lvar2<lvar1; lvar2++) {
                  lvar4=ha_cgefind(&linecopy[lvar3+4],lintmp);
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
            if (strcmp(ha_var[l].cofname,vname)==0) {
              LinVars[i3].LinVarIndx=l;
              break;
            }
          }
          i3++;
        }
      }
      nlinvars=i3;
      for (i=0; i<nlinvars; i++) i3+=ha_var[LinVars[i].LinVarIndx].size;
      ha_cgesetindx *arSet= (ha_cgesetindx *) calloc (fdim+i3+1,sizeof(ha_cgesetindx));
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
          strcpy(arSet[i].arIndx,readitem);
          readitem = strtok(NULL,")");
          for (i4=0; i4<nset; i4++) if(strcmp(readitem,ha_set[i4].setname)==0) {
              arSet[i].setid=i4;
              ha_eq[eqindx].setid[i]=i4;
              break;
            }
          nloops=nloops*ha_set[arSet[i].setid].size;
        }
      }
      ha_eq[eqindx].size=fdim;
      ha_eq[eqindx].matsize=nloops;
      if(alltimeset>=0) {
        if(eqindx==0) {
          for(i=0; i<nvar; i++) {
            if(var_inter[i]) continue;
            j=0;
            for(i4=0; i4<ha_var[i].size; i4++)if(ha_set[ha_var[i].setid[i4]].intertemp) {
                orderintra[i]=i4;
                j++;
                break;
              }
            if(j==0)var_inter[i]=true;
          }
        }
        j=0;
        for (i=0; i<fdim; i++)if(ha_set[ha_eq[eqindx].setid[i]].intertemp) {
            j++;
            eq_orderintra[eqindx]=i;
          }
        if(j==0) {
          if(allregset<0)for (i4=0; i4<nlinvars; i4++) var_inter[LinVars[i4].LinVarIndx]=true;
          ha_eqint[eqindx]=true;
        }

      }
      if(allregset>=0) {
        if(eqindx==0) {
          for(i=0; i<nvar; i++) {
            j=0;
            for(i4=0; i4<ha_var[i].size; i4++)if(ha_set[ha_var[i].setid[i4]].regional) {
                orderreg[i]=i4;
                j++;
              }
            if(j==0)var_inter[i]=true;
          }
        }
        j=0;
        for(i4=0; i4<nlinvars; i4++)for(i=0; i<ha_var[LinVars[i4].LinVarIndx].size; i++) if(ha_set[ha_var[LinVars[i4].LinVarIndx].setid[i]].regional)j++;
        if(j==0)ha_eqint[eqindx]=true;
        j=0;
        for (i=0; i<fdim; i++)if(ha_set[ha_eq[eqindx].setid[i]].regional) {
            j++;
            eq_orderreg[eqindx]=i;
          }
        if(j==0) {
          ha_eqint[eqindx]=true;
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
              for(i3=0; i3<ha_var[LinVars[i4].LinVarIndx].size; i3++) {
                if(strcmp(LinVars[i4].dimnames[i3],arSet[i].arIndx)==0&&ha_set[ha_eq[eqindx].setid[i]].regional&&orderreg[LinVars[i4].LinVarIndx]==i3) { //ha_set[ha_var[LinVars[i4].LinVarIndx].setid[i3]].regional) {
                  dimmat[i]+=ha_var[LinVars[i4].LinVarIndx].matsize;
                  l01=0;
                  for(j01=0; j01<ha_var[LinVars[i4].LinVarIndx].matsize; j01++)if(ha_cgeshock[ha_var[LinVars[i4].LinVarIndx].begadd+j01].ShockId)l01++;
                  dimmat[i]-=l01;
                  linvarrcount[i4]=true;
                  if(lvar4<i4)linvarrcount[lvar4]=true;
                }
              }
            }
          }
          matsize=0;
          for (i=0; i<fdim; i++)if(dimmat[i]>matsize)matsize=dimmat[i];
          for (i=0; i<fdim; i++)if(dimmat[i]==matsize)break;
          if(matsize==0) {
            eq_orderreg[eqindx]=tempint;
            ha_eqint[eqindx]=true;
          }
          else {
            eq_orderreg[eqindx]=i;
          }

        }


        if(eq_orderreg[eqindx]>-1) {
          matsize=0;
          memset(linvarrcount,0,nlinvars*sizeof(bool));
          for(i4=0; i4<nlinvars; i4++) {
            for(lvar4=0; lvar4<i4; lvar4++)if(strcmp(LinVars[i4].LinVarName,LinVars[lvar4].LinVarName)==0)break;
            if(lvar4<i4&&linvarrcount[lvar4]) {
              linvarrcount[i4]=true;
              continue;
            }
            if(strcmp(LinVars[i4].dimnames[orderreg[LinVars[i4].LinVarIndx]],arSet[eq_orderreg[eqindx]].arIndx)!=0&&!var_inter[LinVars[i4].LinVarIndx]) {
              matsize+=ha_var[LinVars[i4].LinVarIndx].matsize;
              l01=0;
              for(j01=0; j01<ha_var[LinVars[i4].LinVarIndx].matsize; j01++)if(ha_cgeshock[ha_var[LinVars[i4].LinVarIndx].begadd+j01].ShockId)l01++;
              matsize+=ha_var[LinVars[i4].LinVarIndx].matsize-l01;
              linvarrcount[i4]=true;
              if(lvar4<i4)linvarrcount[lvar4]=true;
            }
          }
          if(matsize>=ha_eq[eqindx].matsize)ha_eqint[eqindx]=true;
          else {
            for(i4=0; i4<nlinvars; i4++)
              if(strcmp(LinVars[i4].dimnames[orderreg[LinVars[i4].LinVarIndx]],arSet[eq_orderreg[eqindx]].arIndx)!=0&&!var_inter[LinVars[i4].LinVarIndx]) {
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

int NestedMatvarRead(char *fname, char *commsyntax,ha_cgeset *ha_set,uvdim nset,ha_cgesetele *ha_setele,hcge_cof *ha_cof,uvadd ncof,hcge_cof *ha_var,uvadd nvar,ha_cgevar *ha_cofvar,uvadd ncofvar,uvadd ncofele,ha_cgeexovar *ha_cgeshock,bool *var_inter,hcge_cof *ha_eq,bool *ha_eqint,uvdim *eq_orderintra,uvdim *eq_orderreg,uvadd allregset,uvadd alltimeset,uvdim *orderintra,uvdim *orderreg) {
  FILE * filehandle;
  char tline[TABREADLINE],line[TABREADLINE],line1[TABREADLINE],linecopy[TABREADLINE];//,set1[NAMESIZE],set2[NAMESIZE];
  char vname[TABREADLINE],lintmp[TABREADLINE];//,*p1=NULL;
  char *readitem=NULL,*p=NULL,*p1=NULL;//,*p2=NULL,*varpnts;
  ha_cgetype zerodivide=0;
  uvdim fdim,np;
  uvadd rowindx,j,j01,l,l01,lj,dcountdim1[4*MAXVARDIM],dcountdim2[4*MAXVARDIM],nloops,nloopslin,li3,l2,matrow,eqindx=0,matsize;//,sizelinvars,totlinvars,templinvars
  uvadd sj,i,i3,i4,nlinvars,lvar,lvar1,lvar2,lvar3,lvar4,dcount,fdimlin=0,varindx1,varindx2,tempint;
  int leadlag;
  ha_cgetype dimmat[4*MAXVARDIM];

  filehandle = fopen(fname,"r");
  matrow=0;

  while (ha_cgertabl1(commsyntax,filehandle,line,ha_cofvar,ha_cof,ncof,&zerodivide,TABREADLINE)) {
    if (strstr(line,"(default")==NULL) {
      ha_cgefrstr1(line, commsyntax, "");
      ha_cgefrstr1(line, "(linear)", "");
      while (ha_cgefrstr(line,"  ", " "));
      while (ha_cgefrchr(line, '[', '('));
      while (ha_cgefrchr(line, ']', ')'));
      while (ha_cgefrchr(line, '{', '('));
      while (ha_cgefrchr(line, '}', ')'));
      strcpy(linecopy,line);
      fdim=ha_cgenfind(line, "(all,");
      if (fdim==0) {
        readitem = strtok(line+1," ");
        strcpy(ha_eq[eqindx].cofname,readitem);
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
        strcpy(ha_eq[eqindx].cofname,readitem);
        strcpy(line,linecopy);
        i=ha_cgerevfind(line, "(all,");
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
      while (ha_cgefrstr(readitem," ", ""));
      while (ha_cgerecovar(readitem)==1);
      hnew_intrpl(readitem);
      strcpy(tline,readitem);
      strcpy(line1,readitem);
      np=ha_cgenfind(readitem,"p_");
      hcge_linvars *LinVars= (hcge_linvars *) calloc (np+1,sizeof(hcge_linvars));
      i3=0;
      lvar=0;
      for (i=0; i<np; i++) {
        varindx2=0;
        while(-1<0) {
          varindx1=ha_cgefind(readitem+varindx2,"p_");
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
            if (strcmp(ha_var[l].cofname,vname)==0) {
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
          switch (ha_var[l].size) {
          case 0:
            break;
          case 1:
            p = strtok(tline,"{");
            p = strtok(NULL,"}");
            leadlag=0;
            hnew_arset(p,&leadlag);
            if(leadlag!=0) {
              var_inter[l]=true;
              orderintra[l]=-1;
              orderreg[l]=-1;
            }//printf("var %s\n",ha_var[l].cofname);}
            strcpy(LinVars[i3].dimnames[0],p);
            strcpy(lintmp,"(all,");
            strcat(lintmp,p);
            strcat(lintmp,",");
            l01=ha_cgefind(linecopy,lintmp);
            if (l01>-1) {
              p1=&linecopy[0]+l01;
              strncpy(LinVars[i3].dimsetnames[0],p1+strlen(lintmp),strchr(p1,')')-p1-strlen(lintmp));
            }
            else {
              strcpy(lintmp,"sum(");
              strcat(lintmp,p);
              strcat(lintmp,",");
              lvar1=ha_cgenfind(linecopy,lintmp);
              lvar3=ha_cgefind(linecopy,lintmp);
              if (lvar1>1) for(lvar2=0; lvar2<lvar1; lvar2++) {
                  lvar4=ha_cgefind(&linecopy[lvar3+4],lintmp);
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
            for (i4=0; i4<ha_var[l].size-1; i4++) {
              p = strtok(NULL,",");
              leadlag=0;
              hnew_arset(p,&leadlag);
              if(leadlag!=0) {
                var_inter[l]=true;
                orderintra[l]=-1;
                orderreg[l]=-1;
              }//printf("var %s\n",ha_var[l].cofname);}
              strcpy(LinVars[i3].dimnames[i4],p);
              strcpy(lintmp,"(all,");
              strcat(lintmp,p);
              strcat(lintmp,",");
              l01=ha_cgefind(linecopy,lintmp);
              if (l01>-1) {
                p1=&linecopy[0]+l01;
                strncpy(LinVars[i3].dimsetnames[i4],p1+strlen(lintmp),strchr(p1,')')-p1-strlen(lintmp));
              }
              else {
                strcpy(lintmp,"sum(");
                strcat(lintmp,p);
                strcat(lintmp,",");
                lvar1=ha_cgenfind(linecopy,lintmp);
                lvar3=ha_cgefind(linecopy,lintmp);
                if (lvar1>1) for(lvar2=0; lvar2<lvar1; lvar2++) {
                    lvar4=ha_cgefind(&linecopy[lvar3+4],lintmp);
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
            hnew_arset(p,&leadlag);
            if(leadlag!=0) {
              var_inter[l]=true;
              orderintra[l]=-1;
              orderreg[l]=-1;
            }//printf("var %s\n",ha_var[l].cofname);}
            strcpy(LinVars[i3].dimnames[i4],p);
            strcpy(lintmp,"(all,");
            strcat(lintmp,p);
            strcat(lintmp,",");
            l01=ha_cgefind(linecopy,lintmp);
            if (l01>-1) {
              p1=&linecopy[0]+l01;
              strncpy(LinVars[i3].dimsetnames[i4],p1+strlen(lintmp),strchr(p1,')')-p1-strlen(lintmp));
            }
            else {
              strcpy(lintmp,"sum(");
              strcat(lintmp,p);
              strcat(lintmp,",");
              lvar1=ha_cgenfind(linecopy,lintmp);
              lvar3=ha_cgefind(linecopy,lintmp);
              if (lvar1>1) for(lvar2=0; lvar2<lvar1; lvar2++) {
                  lvar4=ha_cgefind(&linecopy[lvar3+4],lintmp);
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
            if (strcmp(ha_var[l].cofname,vname)==0) {
              LinVars[i3].LinVarIndx=l;
              break;
            }
          }
          i3++;
        }
      }
      nlinvars=i3;
      for (i=0; i<nlinvars; i++) i3+=ha_var[LinVars[i].LinVarIndx].size;
      ha_cgesetindx *arSet= (ha_cgesetindx *) calloc (fdim+i3+1,sizeof(ha_cgesetindx));
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
          strcpy(arSet[i].arIndx,readitem);
          readitem = strtok(NULL,")");
          for (i4=0; i4<nset; i4++) if(strcmp(readitem,ha_set[i4].setname)==0) {
              arSet[i].setid=i4;
              ha_eq[eqindx].setid[i]=i4;
              break;
            }
          nloops=nloops*ha_set[arSet[i].setid].size;
        }
      }
      ha_eq[eqindx].size=fdim;
      ha_eq[eqindx].matsize=nloops;


      if(eqindx==0) {
        for(i=0; i<nvar; i++) {
          if(var_inter[i]) continue;
          j=0;
          for(i4=0; i4<ha_var[i].size; i4++)if(ha_set[ha_var[i].setid[i4]].intertemp) {
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
            for(i4=0; i4<ha_var[i].size; i4++)if(ha_set[ha_var[i].setid[i4]].regional) {
                orderreg[i]=i4;
                j01++;
              }
          }
        }
      }


      j=0;
      for(i4=0; i4<nlinvars; i4++)for(i=0; i<ha_var[LinVars[i4].LinVarIndx].size; i++) if(ha_set[ha_var[LinVars[i4].LinVarIndx].setid[i]].intertemp)j++;
      if(j==0) {
        ha_eqint[eqindx]=true;//eq_orderreg[eqindx]=-1;eq_orderintra[eqindx]=-1;}
      }
      j=0;
      for (i=0; i<fdim; i++)if(ha_set[ha_eq[eqindx].setid[i]].regional) {
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
            for(i3=0; i3<ha_var[LinVars[i4].LinVarIndx].size; i3++) {
              if(strcmp(LinVars[i4].dimnames[i3],arSet[i].arIndx)==0&&ha_set[ha_eq[eqindx].setid[i]].regional&&orderreg[LinVars[i4].LinVarIndx]==i3) {
                dimmat[i]+=ha_var[LinVars[i4].LinVarIndx].matsize;
                linvarrcount[i4]=true;
                if(lvar4<i4)linvarrcount[lvar4]=true;
                l01=0;
                for(j01=0; j01<ha_var[LinVars[i4].LinVarIndx].matsize; j01++)if(ha_cgeshock[ha_var[LinVars[i4].LinVarIndx].begadd+j01].ShockId)l01++;
                dimmat[i]-=l01;
              }
            }
          }
        }
        matsize=0;
        for (i=0; i<fdim; i++)if(dimmat[i]>matsize)matsize=dimmat[i];
        for (i=0; i<fdim; i++)if(dimmat[i]==matsize)break;
        if(matsize==0) {
          eq_orderreg[eqindx]=-1;
        }
        else {
          eq_orderreg[eqindx]=i;
        }
      }
      if(eq_orderreg[eqindx]>-1) {
        matsize=0;
        memset(linvarrcount,0,nlinvars*sizeof(bool));
        for(i4=0; i4<nlinvars; i4++) {
          for(lvar4=0; lvar4<i4; lvar4++)if(strcmp(LinVars[i4].LinVarName,LinVars[lvar4].LinVarName)==0)break;
          if(lvar4<i4&&linvarrcount[lvar4]) {
            linvarrcount[i4]=true;
            continue;
          }
          if(strcmp(LinVars[i4].dimnames[orderreg[LinVars[i4].LinVarIndx]],arSet[eq_orderreg[eqindx]].arIndx)!=0&&orderreg[LinVars[i4].LinVarIndx]!=-1) {
            l01=0;
            for(j01=0; j01<ha_var[LinVars[i4].LinVarIndx].matsize; j01++)if(ha_cgeshock[ha_var[LinVars[i4].LinVarIndx].begadd+j01].ShockId)l01++;
            matsize+=ha_var[LinVars[i4].LinVarIndx].matsize-l01;
            linvarrcount[i4]=true;
            if(lvar4<i4)linvarrcount[lvar4]=true;
          }
        }
        if(matsize>=ha_eq[eqindx].matsize)eq_orderreg[eqindx]=-1;//printf("OOOOOOOOOOOOOOOOO\n");}
        else {
          for(i4=0; i4<nlinvars; i4++)
            if(strcmp(LinVars[i4].dimnames[orderreg[LinVars[i4].LinVarIndx]],arSet[eq_orderreg[eqindx]].arIndx)!=0&&orderreg[LinVars[i4].LinVarIndx]!=-1) {
              orderreg[LinVars[i4].LinVarIndx]=-1;//orderintra[LinVars[i4].LinVarIndx]=-1;
            }
        }
      }


      j=0;
      for (i=0; i<fdim; i++)if(ha_set[ha_eq[eqindx].setid[i]].intertemp) {
          j++;
          eq_orderintra[eqindx]=i;
        }
      if(j==0) {
        ha_eqint[eqindx]=true;
      }

      if(ha_eqint[eqindx]==true) {
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

int NewMatreadele(char *fname, char *commsyntax,ha_cgeset *ha_set,uvdim nset,ha_cgesetele *ha_setele,hcge_cof *ha_cof,uvadd ncof,hcge_cof *ha_var,uvadd nvar,ha_cgevar *ha_cofvar,uvadd ncofvar,uvadd ncofele, uvadd nexo,ha_cgeexovar *ha_cgeshock,uvadd ndblock,uvadd alltimeset,uvadd allregset,bool *ha_eqint,uvadd *ha_eqadd,uvdim *ha_eqtime,uvdim *ha_eqreg,uvadd *counteq,uvadd nintraeq,bool *sbbd_overrid,PetscInt Istart,PetscInt Iend,PetscInt *dnz,PetscInt *dnnz,PetscInt *onz,PetscInt *onnz,PetscInt *dnzB,PetscInt *dnnzB,PetscInt *onzB,PetscInt *onnzB,int nesteddbbd) {
  FILE * filehandle;
  char tline[TABREADLINE],line[TABREADLINE],line1[TABREADLINE],linecopy[TABREADLINE];//,set1[NAMESIZE],set2[NAMESIZE];
  char vname[TABREADLINE],lintmp[TABREADLINE];//,*p1=NULL;
  char *readitem=NULL,*p=NULL,*p1=NULL;//,*p2=NULL,*varpnts;
  PetscInt Iindx=0,Jindx;
  ha_cgetype zerodivide=0;
  uvdim fdim=0,np,i4,sup,supset[MAXSUPSET];
  uvadd rowindx,rowindxorg,l,l1,lj,dcountdim1[4*MAXVARDIM],dcountdim2[4*MAXVARDIM],dcountdim3[4*MAXVARDIM],dcountdim4[4*MAXVARDIM],dcountdim5[4*MAXVARDIM],nloops,nloopslin,nloopsfac,li3,l2,matrow,matroworg,ltime,lreg,leq=0,eqindx=0;//,sizelinvars,totlinvars,templinvars
  uvadd nreg=0,nint=0,sj,i,i3;
  if(allregset>-1)nreg=ha_set[allregset].size;
  if(alltimeset>-1)nint=ha_set[alltimeset].size;
  int nlinvars,lvar,lvar1,lvar2,lvar3,lvar4,dcount,fdimlin=0,leadlag,varindx1,varindx2;
  uvadd *counteq1= (uvadd *) calloc (ndblock,sizeof(uvadd));
  for(i=0; i<ndblock; i++)counteq1[i]=counteq[i];
  filehandle = fopen(fname,"r");
  matrow=nintraeq;
  matroworg=0;
  Jindx=0;

  while (ha_cgertabl1(commsyntax,filehandle,line,ha_cofvar,ha_cof,ncof,&zerodivide,TABREADLINE)) {
    if (strstr(line,"(default")==NULL) {
      ha_cgefrstr1(line, commsyntax, "");
      ha_cgefrstr1(line, "(linear)", "");
      while (ha_cgefrstr(line,"  ", " "));
      while (ha_cgefrchr(line, '[', '('));
      while (ha_cgefrchr(line, ']', ')'));
      while (ha_cgefrchr(line, '{', '('));
      while (ha_cgefrchr(line, '}', ')'));
      strcpy(linecopy,line);
      fdim=ha_cgenfind(line, "(all,");
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
        i=ha_cgerevfind(line, "(all,");
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
      while (ha_cgefrstr(readitem," ", ""));
      while (ha_cgerecovar(readitem)==1);
      hnew_intrpl(readitem);
      strcpy(tline,readitem);
      strcpy(line1,readitem);
      np=ha_cgenfind(readitem,"p_");
      hcge_linvars *LinVars= (hcge_linvars *) calloc (np+1,sizeof(hcge_linvars));
      i3=0;
      lvar=0;
      for (i=0; i<np; i++) {
        varindx2=0;
        while(-1<0) {
          varindx1=ha_cgefind(readitem+varindx2,"p_");
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
            if (strcmp(ha_var[l].cofname,vname)==0) {
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
          switch (ha_var[l].size) {
          case 0:
            break;
          case 1:
            p = strtok(tline,"{");
            p = strtok(NULL,"}");
            leadlag=0;
            hnew_arset(p,&leadlag);
            strcpy(LinVars[i3].dimnames[0],p);
            LinVars[i3].dimleadlag[0]=leadlag;
            strcpy(lintmp,"(all,");
            strcat(lintmp,p);
            strcat(lintmp,",");
            l1=ha_cgefind(linecopy,lintmp);
            if (l1>-1) {
              p1=&linecopy[0]+l1;
              strncpy(LinVars[i3].dimsetnames[0],p1+strlen(lintmp),strchr(p1,')')-p1-strlen(lintmp));
            }
            else {
              strcpy(lintmp,"sum(");
              strcat(lintmp,p);
              strcat(lintmp,",");
              lvar1=ha_cgenfind(linecopy,lintmp);
              lvar3=ha_cgefind(linecopy,lintmp);
              if (lvar1>1) for(lvar2=0; lvar2<lvar1; lvar2++) {
                  lvar4=ha_cgefind(&linecopy[lvar3+4],lintmp);
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
            for (i4=0; i4<ha_var[l].size-1; i4++) {
              p = strtok(NULL,",");
              leadlag=0;
              hnew_arset(p,&leadlag);
              strcpy(LinVars[i3].dimnames[i4],p);
              LinVars[i3].dimleadlag[i4]=leadlag;
              strcpy(lintmp,"(all,");
              strcat(lintmp,p);
              strcat(lintmp,",");
              l1=ha_cgefind(linecopy,lintmp);
              if (l1>-1) {
                p1=&linecopy[0]+l1;
                strncpy(LinVars[i3].dimsetnames[i4],p1+strlen(lintmp),strchr(p1,')')-p1-strlen(lintmp));
              }
              else {
                strcpy(lintmp,"sum(");
                strcat(lintmp,p);
                strcat(lintmp,",");
                lvar1=ha_cgenfind(linecopy,lintmp);
                lvar3=ha_cgefind(linecopy,lintmp);
                if (lvar1>1) for(lvar2=0; lvar2<lvar1; lvar2++) {
                    lvar4=ha_cgefind(&linecopy[lvar3+4],lintmp);
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
            hnew_arset(p,&leadlag);
            strcpy(LinVars[i3].dimnames[i4],p);
            LinVars[i3].dimleadlag[i4]=leadlag;
            strcpy(lintmp,"(all,");
            strcat(lintmp,p);
            strcat(lintmp,",");
            l1=ha_cgefind(linecopy,lintmp);
            if (l1>-1) {
              p1=&linecopy[0]+l1;
              strncpy(LinVars[i3].dimsetnames[i4],p1+strlen(lintmp),strchr(p1,')')-p1-strlen(lintmp));
            }
            else {
              strcpy(lintmp,"sum(");
              strcat(lintmp,p);
              strcat(lintmp,",");
              lvar1=ha_cgenfind(linecopy,lintmp);
              lvar3=ha_cgefind(linecopy,lintmp);
              if (lvar1>1) for(lvar2=0; lvar2<lvar1; lvar2++) {
                  lvar4=ha_cgefind(&linecopy[lvar3+4],lintmp);
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
            if (strcmp(ha_var[l].cofname,vname)==0) {
              LinVars[i3].LinVarIndx=l;
              break;
            }
          }
          i3++;
        }
      }
      nlinvars=i3;
      for (i=0; i<nlinvars; i++) i3+=ha_var[LinVars[i].LinVarIndx].size;
      ha_cgesetindx *arSet= (ha_cgesetindx *) calloc (fdim+i3+1,sizeof(ha_cgesetindx));

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
          strcpy(arSet[i].arIndx,readitem);
          readitem = strtok(NULL,")");
          for (i4=0; i4<nset; i4++) if(strcmp(readitem,ha_set[i4].setname)==0) {
              arSet[i].setid=i4;
              if(ha_set[i4].intertemp) *sbbd_overrid=true;
              break;
            }
          nloops=nloops*ha_set[arSet[i].setid].size;
        }
        //*****Adjust ha_eqreg*******************//
        //*****End Adjust ha_eqreg*******************//
        dcountdim1[fdim-1]=1;
        for (dcount=fdim-2; dcount>-1; dcount--) {
          dcountdim1[dcount]=ha_set[arSet[dcount+1].setid].size*dcountdim1[dcount+1];
        }
        if(alltimeset>=0&&allregset<0) {
          for (dcount=ha_eqtime[eqindx]-1; dcount>-1; dcount--) {
            dcountdim3[dcount]=dcountdim1[dcount]/ha_set[arSet[ha_eqtime[eqindx]].setid].size;
          }
          for (dcount=fdim-1; dcount>ha_eqtime[eqindx]-1; dcount--) {
            dcountdim3[dcount]=dcountdim1[dcount];
          }
        }
        if(alltimeset<0&&allregset>=0) {
          for (dcount=ha_eqreg[eqindx]-1; dcount>-1; dcount--) {
            dcountdim4[dcount]=dcountdim1[dcount]/ha_set[arSet[ha_eqreg[eqindx]].setid].size;
          }
          for (dcount=fdim-1; dcount>ha_eqreg[eqindx]-1; dcount--) {
            dcountdim4[dcount]=dcountdim1[dcount];
          }
        }
        if(alltimeset>=0&&allregset>=0) {
          i3=1;
          dcountdim4[fdim-1]=dcountdim1[fdim-1];
          for (dcount=fdim-2; dcount>-1; dcount--) {
            if(dcount+1==ha_eqreg[eqindx])i3*=ha_set[arSet[ha_eqreg[eqindx]].setid].size;
            if(dcount+1==ha_eqtime[eqindx])i3*=ha_set[arSet[ha_eqtime[eqindx]].setid].size;
            dcountdim4[dcount]=dcountdim1[dcount]/i3;
          }
        }
      }
      for (i=0; i<nlinvars; i++) {
        i3=0;
        nloopslin=nloops;
        if (fdim==0) {
          for (l2=0; l2<ha_var[LinVars[i].LinVarIndx].size; l2++) {
            strcpy(arSet[i3].arIndx,LinVars[i].dimnames[l2]);
            for (i4=0; i4<nset; i4++) if(strcmp(ha_set[i4].setname,LinVars[i].dimsetnames[l2])==0) {
                break;
              }
            arSet[i3].setid=i4;
            nloopslin=nloopslin*ha_set[arSet[i3].setid].size;
            if(ha_set[i4].intertemp)*sbbd_overrid=true;
            i3++;
          }
        }
        else for (l2=0; l2<ha_var[LinVars[i].LinVarIndx].size; l2++) {
            for (sj=0; sj<fdim; sj++) {
              if(strcmp(arSet[sj].arIndx,LinVars[i].dimnames[l2])==0) {
                break;
              }
              else if (sj==fdim-1) {
                strcpy(arSet[fdim+i3].arIndx,LinVars[i].dimnames[l2]);
                for (i4=0; i4<nset; i4++) if(strcmp(ha_set[i4].setname,LinVars[i].dimsetnames[l2])==0) {
                    break;
                  }
                arSet[fdim+i3].setid=i4;
                nloopslin=nloopslin*ha_set[arSet[fdim+i3].setid].size;
                if(ha_set[i4].intertemp)*sbbd_overrid=true;
                i3++;
              }
            }
          }
        fdimlin=fdim+i3;
        if (i3>0) {
          dcountdim2[fdimlin-1]=1;
          for (dcount=fdimlin-2; dcount>-1; dcount--) {
            dcountdim2[dcount]=ha_set[arSet[dcount+1].setid].size*dcountdim2[dcount+1];
          }
        }
        else for (dcount=0; dcount<fdim; dcount++) {
            dcountdim2[dcount]=dcountdim1[dcount];
          }
        nloopsfac=(uvadd)nloopslin/nloops;
        for (dcount=0; dcount<ha_var[LinVars[i].LinVarIndx].size; dcount++) {
          for (i4=0; i4<fdimlin; i4++) {
            if (strcmp(LinVars[i].dimnames[dcount],arSet[i4].arIndx)==0) {
              dcountdim5[dcount]=i4;
              break;
            }
          }
        }
        for(dcount=0; dcount<MAXSUPSET; dcount++)supset[dcount]=0;
        for (dcount=0; dcount<ha_var[LinVars[i].LinVarIndx].size; dcount++) {
          if(ha_set[ha_var[LinVars[i].LinVarIndx].setid[dcount]].size!=ha_set[arSet[dcountdim5[dcount]].setid].size) {
            for(sup=1; sup<MAXSUPSET; sup++)if(ha_var[LinVars[i].LinVarIndx].setid[dcount]==ha_set[arSet[dcountdim5[dcount]].setid].subsetid[sup]) {
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
            l1=(uvadd) l2/dcountdim2[dcount];
            arSet[dcount].indx=l1;
            if(ha_eqint[eqindx]) {
              if(alltimeset>=0&&allregset<0) {
                if (dcount<fdim&&dcount!=ha_eqtime[eqindx]) {
                  rowindx=rowindx+l1*dcountdim3[dcount];
                }
                if(dcount==ha_eqtime[eqindx])ltime=l1;
              }
              if(alltimeset<0&&allregset>=0) {
                if (dcount<fdim&&dcount!=ha_eqreg[eqindx]) {
                  rowindx=rowindx+l1*dcountdim4[dcount];
                }
                if(dcount==ha_eqreg[eqindx])ltime=l1;
              }
              if(alltimeset>=0&&allregset>=0) {
                if(nesteddbbd==1) {
                  if(ha_eqreg[eqindx]>-1) {
                    if (dcount<fdim&&dcount!=ha_eqtime[eqindx]&&dcount!=ha_eqreg[eqindx])rowindx=rowindx+l1*dcountdim4[dcount];
                  }
                  else {
                    if (dcount<fdim&&dcount!=ha_eqtime[eqindx])rowindx=rowindx+l1*dcountdim4[dcount];
                  }
                  if(dcount==ha_eqtime[eqindx])ltime=l1;
                  if(dcount==ha_eqreg[eqindx])lreg=l1;
                }
                else {
                  if (dcount<fdim&&dcount!=ha_eqtime[eqindx]&&dcount!=ha_eqreg[eqindx]) {
                    rowindx=rowindx+l1*dcountdim4[dcount];
                  }
                  if(dcount==ha_eqtime[eqindx])ltime=l1;
                  if(dcount==ha_eqreg[eqindx])lreg=l1;
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
          for (dcount=0; dcount<ha_var[LinVars[i].LinVarIndx].size; dcount++) {
            if(supset[dcount]==0) {
              li3=li3+(arSet[dcountdim5[dcount]].indx+LinVars[i].dimleadlag[dcount])*ha_var[LinVars[i].LinVarIndx].antidims[dcount];
            }
            else {
              li3=li3+(ha_setele[ha_set[arSet[dcountdim5[dcount]].setid].begadd+arSet[dcountdim5[dcount]].indx].setsh[supset[dcount]]+LinVars[i].dimleadlag[dcount])*ha_var[LinVars[i].LinVarIndx].antidims[dcount];
            }
          }

          if(ha_eqint[eqindx]) {
            if(alltimeset>=0&&allregset>=0) {
              if(nesteddbbd==1) {
                if(ha_eqreg[eqindx]>-1)Jindx=counteq1[(ha_setele[ha_set[arSet[ha_eqtime[eqindx]].setid].begadd+ltime].setsh[ha_set[arSet[ha_eqtime[eqindx]].setid].intsup])*(nreg+1)+ha_setele[ha_set[arSet[ha_eqreg[eqindx]].setid].begadd+lreg].setsh[ha_set[arSet[ha_eqreg[eqindx]].setid].regsup]]+rowindx;
                else Jindx=counteq1[(ha_setele[ha_set[arSet[ha_eqtime[eqindx]].setid].begadd+ltime].setsh[ha_set[arSet[ha_eqtime[eqindx]].setid].intsup])*(nreg+1)+nreg]+rowindx;
              }
              else {
                if(arSet[ha_eqtime[eqindx]].setid==alltimeset)
                  Jindx=counteq1[(ha_setele[ha_set[arSet[ha_eqtime[eqindx]].setid].begadd+ltime].setsh[0])*ha_set[arSet[ha_eqreg[eqindx]].setid].size+ha_setele[ha_set[arSet[ha_eqreg[eqindx]].setid].begadd+lreg].setsh[0]]+rowindx;
                else {
                  for(i4=1; i4<MAXSUPSET; i4++)if(ha_set[arSet[ha_eqtime[eqindx]].setid].subsetid[i4]=alltimeset)break;
                  Jindx=counteq1[(ha_setele[ha_set[arSet[ha_eqtime[eqindx]].setid].begadd+ltime].setsh[i4])*ha_set[arSet[ha_eqreg[eqindx]].setid].size+ha_setele[ha_set[arSet[ha_eqreg[eqindx]].setid].begadd+lreg].setsh[0]]+rowindx;
                }
              }
            }
            if(alltimeset>=0&&allregset<0) {
              if(arSet[ha_eqtime[eqindx]].setid==alltimeset)
                Jindx=counteq1[ha_setele[ha_set[arSet[ha_eqtime[eqindx]].setid].begadd+ltime].setsh[0]]+rowindx;
              else {
                for(i4=1; i4<MAXSUPSET; i4++)if(ha_set[arSet[ha_eqtime[eqindx]].setid].subsetid[i4]=alltimeset)break;
                Jindx=counteq1[ha_setele[ha_set[arSet[ha_eqtime[eqindx]].setid].begadd+ltime].setsh[i4]]+rowindx;
              }
            }
            if(alltimeset<0&&allregset>=0) Jindx=counteq1[ha_setele[ha_set[arSet[ha_eqreg[eqindx]].setid].begadd+ltime].setsh[0]]+rowindx;
          }
          else Jindx=matrow+rowindx;
          leq=matroworg+(uvadd)lj/nloopsfac;//matroworg+rowindxorg;
          ha_eqadd[leq]=Jindx;
          Iindx=ha_cgeshock[ha_var[LinVars[i].LinVarIndx].begadd+li3].ExoIndx;
          if(Istart<=Jindx&&Jindx<Iend) {
            if (!ha_cgeshock[ha_var[LinVars[i].LinVarIndx].begadd+li3].ShockId) {
              if (Istart<=Iindx&&Iindx<Iend) {
                dnnz[Jindx-Istart]=dnnz[Jindx-Istart]+1;
              }
              else {
                onnz[Jindx-Istart]=onnz[Jindx-Istart]+1;
              }
            }
            if (ha_cgeshock[ha_var[LinVars[i].LinVarIndx].begadd+li3].ShockId) {
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
      if(ha_eqint[eqindx]) {
        if(alltimeset>=0&&allregset>=0) {
          if(nesteddbbd==1) {
            if(ha_eqreg[eqindx]>-1)for(lj=0; lj<ha_set[arSet[ha_eqtime[eqindx]].setid].size; lj++)for(l2=0; l2<ha_set[arSet[ha_eqreg[eqindx]].setid].size; l2++) {
                  counteq1[ha_setele[ha_set[arSet[ha_eqtime[eqindx]].setid].begadd+lj].setsh[ha_set[arSet[ha_eqtime[eqindx]].setid].intsup]*(nreg+1)+ha_setele[ha_set[arSet[ha_eqreg[eqindx]].setid].begadd+l2].setsh[ha_set[arSet[ha_eqreg[eqindx]].setid].regsup]]+=nloops/ha_set[arSet[ha_eqtime[eqindx]].setid].size/ha_set[arSet[ha_eqreg[eqindx]].setid].size;
                  if(ha_set[arSet[ha_eqreg[eqindx]].setid].regsup!=0)printf("Sub reg not supported in NDBBD!!!\n");
                }
            else for(lj=0; lj<ha_set[arSet[ha_eqtime[eqindx]].setid].size; lj++)
                counteq1[ha_setele[ha_set[arSet[ha_eqtime[eqindx]].setid].begadd+lj].setsh[ha_set[arSet[ha_eqtime[eqindx]].setid].intsup]*(nreg+1)+nreg]+=nloops/ha_set[arSet[ha_eqtime[eqindx]].setid].size;
          }
          else {
            if(arSet[ha_eqtime[eqindx]].setid==alltimeset)
              for(lj=0; lj<ha_set[arSet[ha_eqtime[eqindx]].setid].size; lj++)for(l2=0; l2<ha_set[arSet[ha_eqreg[eqindx]].setid].size; l2++)
                  counteq1[ha_setele[ha_set[arSet[ha_eqtime[eqindx]].setid].begadd+lj].setsh[0]*ha_set[arSet[ha_eqreg[eqindx]].setid].size+ha_setele[ha_set[arSet[ha_eqreg[eqindx]].setid].begadd+l2].setsh[0]]+=nloops/ha_set[arSet[ha_eqtime[eqindx]].setid].size/ha_set[arSet[ha_eqreg[eqindx]].setid].size;
            else {
              for(lj=0; lj<ha_set[arSet[ha_eqtime[eqindx]].setid].size; lj++)for(l2=0; l2<ha_set[arSet[ha_eqreg[eqindx]].setid].size; l2++) {
                  for(i4=1; i4<MAXSUPSET; i4++)if(ha_set[arSet[ha_eqtime[eqindx]].setid].subsetid[i4]=alltimeset)break;
                  counteq1[ha_setele[ha_set[arSet[ha_eqtime[eqindx]].setid].begadd+lj].setsh[i4]*ha_set[arSet[ha_eqreg[eqindx]].setid].size+ha_setele[ha_set[arSet[ha_eqreg[eqindx]].setid].begadd+l2].setsh[0]]+=nloops/ha_set[arSet[ha_eqtime[eqindx]].setid].size/ha_set[arSet[ha_eqreg[eqindx]].setid].size;
                }
            }
          }
        }
        if(alltimeset>=0&&allregset<0) {
          if(arSet[ha_eqtime[eqindx]].setid==alltimeset)
            for(lj=0; lj<ha_set[arSet[ha_eqtime[eqindx]].setid].size; lj++)counteq1[ha_setele[ha_set[arSet[ha_eqtime[eqindx]].setid].begadd+lj].setsh[0]]+=nloops/ha_set[arSet[ha_eqtime[eqindx]].setid].size;
          else {
            for(i4=1; i4<MAXSUPSET; i4++)if(ha_set[arSet[ha_eqtime[eqindx]].setid].subsetid[i4]=alltimeset)break;
            for(lj=0; lj<ha_set[arSet[ha_eqtime[eqindx]].setid].size; lj++)counteq1[ha_setele[ha_set[arSet[ha_eqtime[eqindx]].setid].begadd+lj].setsh[i4]]+=nloops/ha_set[arSet[ha_eqtime[eqindx]].setid].size;
          }
        }
        if(alltimeset<0&&allregset>=0)for(lj=0; lj<ha_set[arSet[ha_eqreg[eqindx]].setid].size; lj++)counteq1[ha_setele[ha_set[arSet[ha_eqreg[eqindx]].setid].begadd+lj].setsh[0]]+=nloops/ha_set[arSet[ha_eqreg[eqindx]].setid].size;
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

int HaDBBDMatOder(Mat A, uvadd VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend, uvadd nvarele, uvadd *ha_eqadd,int *ha_rows,int *ha_cols, uvadd ndblock,int *ha_ndblocks, uvadd *countvarintra1, uvadd *counteq, uvadd *counteqnoadd,uvdim laA,ha_cgetype cntl6) {
  IS *rowindices,*colindices;//,isrow,iscol;
  PetscInt bfirst,bend,sblockin,nmatin,nmatinplus,nrowcolin,sumrowcolin;
  Mat *submatA;
  PetscInt i,j,j0,j1,j2,j3,j4,nrow,ncol,nz,nz1,*ai,*aj,la;
  PetscScalar *vals;
  PetscErrorCode ierr;
  PetscViewer viewer;
  uvadd lasize;
  int *ha_rows1= (int *) calloc (VecSize,sizeof(int));
  int *ha_cols1= (int *) calloc (VecSize,sizeof(int));
  int *ha_ndblocks1= (int *) calloc (ndblock,sizeof(int));

  //*****New order test**********//
  nmatin=(uvadd)ndblock/mpisize;
  if(nmatin*mpisize<ndblock)nmatinplus=nmatin+1;
  else nmatinplus=nmatin;
  for(i=0; i<mpisize; i++)if(rank+1<=ndblock-mpisize*nmatin)nmatin++;
  int *begblock= (int *) calloc (mpisize,sizeof(int));
  ierr = PetscMalloc(nmatin*sizeof(IS **),&rowindices);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatin*sizeof(IS **),&colindices);
  CHKERRQ(ierr);
  printf("rank %d nmatin %d\n",rank,nmatin);
  begblock[rank]=nmatin;
  for(i=0; i<mpisize; i++) {
    j=rank;
    MPI_Bcast(&j, 1, MPI_INT,i, PETSC_COMM_WORLD);
    MPI_Bcast(&begblock[j], 1, MPI_INT,i, PETSC_COMM_WORLD);
  }
  j2=begblock[0];
  begblock[0]=0;
  for(i=1; i<mpisize; i++) {
    j1=begblock[i];
    begblock[i]=begblock[i-1]+j2;
    j2=j1;
  }
  for(i=0; i<nmatin; i++) {
    bfirst=counteq[i+begblock[rank]];
    bend=counteqnoadd[i+begblock[rank]]+counteq[i+begblock[rank]];
    ISCreateStride(PETSC_COMM_SELF,bend-bfirst,bfirst,1,rowindices+i);
    bfirst=countvarintra1[i+begblock[rank]];
    bend=countvarintra1[i+1+begblock[rank]];
    ISCreateStride(PETSC_COMM_SELF,bend-bfirst,bfirst,1,colindices+i);
  }
  printf("OK???? rank %d\n",rank);
  ierr = MatCreateSubMatrices(A,nmatin,rowindices,colindices,MAT_INITIAL_MATRIX,&submatA);
  CHKERRQ(ierr);
  for (i=0; i<nmatin; i++) {
    ierr = ISDestroy(&rowindices[i]);
    CHKERRQ(ierr);
  }
  for (i=0; i<nmatin; i++) {
    ierr = ISDestroy(&colindices[i]);
    CHKERRQ(ierr);
  }
  PetscFree(rowindices);
  PetscFree(colindices);
  MPI_Barrier(PETSC_COMM_WORLD);
  nrowcolin=0;
  for(j1=0; j1<nmatin; j1++) {
    Mat_SeqAIJ         *aa=(Mat_SeqAIJ*)submatA[j1]->data;//*aa=subA->data;
    ai= aa->i;
    aj= aa->j;
    vals=aa->a;
    nz=aa->nz;
    nrow=submatA[j1]->rmap->n;
    ncol=submatA[j1]->cmap->n;
    printf("rank %d nz %d\n",rank,nz);
    nz1=nz;
    if(nz1<nrow)nz1=nrow;
    if(nz1<ncol)nz1=ncol;
    lasize=ceil((laA/100.0)*nz1);
    int *irn=(int *) calloc (lasize,sizeof(int));
    int *jcn=(int *) calloc (lasize,sizeof(int));
    int *irn1=(int *)malloc(nrow*sizeof(int));
    int *jcn1=(int *)malloc(ncol*sizeof(int));
    ha_cgetype *values= (ha_cgetype *) calloc (lasize,sizeof(ha_cgetype));
    int *keep=(int *) malloc ((nrow+5*ncol+4*ncol+7)*sizeof(int));//ICNTL(6)=1
  ha_cgetype *w51;
  if(nrow>ncol) w51=(ha_cgetype *) malloc ((5*nrow)*sizeof(ha_cgetype));
  else w51=(ha_cgetype *) malloc ((5*ncol)*sizeof(ha_cgetype));
  int *iw51=(int *) malloc ((6*nrow+3*ncol)*sizeof(int));
    int *insize=(int *) calloc (6,sizeof(int));
    insize[0]=nrow;
    insize[1]=ncol;
    insize[2]=nz;
    insize[4]=laA;
    insize[5]=lasize;
    for(i=0; i<nrow-1; i++)for(j=ai[i]; j<ai[i+1]; j++) {
        irn[j]=i+1;
        jcn[j]=aj[j]+1;
        values[j]=vals[j];
      }
    for(j=ai[nrow-1]; j<nz; j++) {
      irn[j]=nrow;
      jcn[j]=aj[j]+1;
      values[j]=vals[j];
    }
    j=0;
    spec51m_rank_(insize,&cntl6,irn,jcn,values,irn1,jcn1,keep,w51,iw51);
    for(i=0; i<nrow; i++) {
      ha_rows1[i+counteq[j1+begblock[rank]]]=irn1[i]-1;
    }
    for(i=0; i<ncol; i++) {
      ha_cols1[i+countvarintra1[j1+begblock[rank]]]=jcn1[i]-1;
    }
    printf("rank %d j1 %d proc %d\n",insize[3],j1,rank);
    ha_ndblocks1[j1+begblock[rank]]=insize[3];
    free(insize);
    free(irn);
    free(jcn);
    free(irn1);
    free(jcn1);
    free(values);
    free(keep);
    free(w51);
    free(iw51);
    MatDestroy(&submatA[j1]);
    nrowcolin+=ncol;
  }
  PetscFree(submatA);
  MPI_Barrier(PETSC_COMM_WORLD);
  MPI_Allreduce(ha_cols1,ha_cols,VecSize, MPI_INT, MPI_SUM,PETSC_COMM_WORLD);
  free(ha_cols1);
  MPI_Allreduce(ha_rows1,ha_rows,VecSize, MPI_INT, MPI_SUM,PETSC_COMM_WORLD);
  free(ha_rows1);
  MPI_Allreduce(ha_ndblocks1,ha_ndblocks,ndblock, MPI_INT, MPI_SUM,PETSC_COMM_WORLD);
  free(ha_ndblocks1);
  free(begblock);
  MPI_Barrier(PETSC_COMM_WORLD);
  return 1;
}


int HaNDBBDMatOderPre(Mat A, uvadd VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int nreg, int ntime, uvadd nvarele, uvadd *ha_eqadd,int *ha_rows,int *ha_cols, uvadd ndblock,int *ha_ndblocks, uvadd *countvarintra1, uvadd *counteq, uvadd *counteqnoadd,uvdim laA,uvdim laDi,ha_cgetype cntl6,PetscInt* ndbbdrank,PetscBool presol) {
  FILE *presolfile;
  char j1name[1024],filename[1024],rankname[1024];
  size_t frd;
  IS *rowindices=NULL,*colindices=NULL,*rowindicesD=NULL,*colindicesD=NULL;//,isrow,iscol;
  PetscInt bfirst,bend,sblockin,nmatin,nmatint,nmatinplus,nrowcolin,sumrowcolin,rfirst,cfirst;
  Mat *submatA=NULL;//,*submatD=NULL;
  PetscInt i,j,j0,j1,j2,j3,j4,j5,j6,j7,j8,nrow,ncol,nz,nz1,*ai,*aj,la;
  uvadd lasize;
  ha_cgetype cntl6in;
  PetscScalar *vals;
  PetscErrorCode ierr;
  PetscViewer viewer;
  MatInfo           matinfo;
  MatGetInfo(A,MAT_LOCAL,&matinfo);
  printf("rank %d matinfo.nz_used %g\n",rank,matinfo.nz_used);
  int *ha_ndblocks1= (int *) calloc (ndblock,sizeof(int));
  uvadd *counteq2= (uvadd *) calloc (ndblock+1,sizeof(uvadd));
  uvadd *counteqnoadd1= (uvadd *) calloc (ndblock,sizeof(uvadd));
  uvadd *counteqnoadd2= (uvadd *) calloc (ndblock,sizeof(uvadd));
  uvadd *countvarintra2= (uvadd *) calloc (ndblock+1,sizeof(uvadd));
  uvadd *countvarintra4= (uvadd *) calloc (ndblock+1,sizeof(uvadd));
  uvadd *countvarintra6= (uvadd *) calloc (ndblock+1,sizeof(uvadd));
  if(rank<10)strcpy(rankname,"000");
  if(rank<100&&rank>9)strcpy(rankname,"00");
  if(rank<1000&&rank>99)strcpy(rankname,"0");
  if(rank>=1000)rankname[0]='\0';
  sprintf(filename, "%d",rank);
  strcat(rankname,filename);
  //*****New order test**********//
  nmatint=(uvadd)ntime/mpisize;
  for(i=0; i<mpisize; i++)if(rank+1<=ntime-mpisize*nmatint)nmatint++;
  nmatin=(nreg+1)*nmatint;
  int *begblock= (int *) calloc (mpisize,sizeof(int));
  begblock[rank]=nmatin;
  printf("rank %d nmatin %d nmatint %d ndblock %ld\n",rank,nmatin,nmatint,ndblock);
  for(i=0; i<mpisize; i++) {
    j=rank;
    j2=j;
    MPI_Bcast(&j2,1, MPI_INT,i, PETSC_COMM_WORLD);
    j=j2;
    MPI_Bcast(&begblock[j2],1, MPI_INT,i, PETSC_COMM_WORLD);
  }
  j2=begblock[0];
  begblock[0]=0;
  for(i=1; i<mpisize; i++) {
    j1=begblock[i];
    begblock[i]=begblock[i-1]+j2;
    j2=j1;
  }
  ierr = PetscMalloc(nmatint*nreg*sizeof(IS **),&rowindices);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatint*nreg*sizeof(IS **),&colindices);
  CHKERRQ(ierr);
  j2=0;
  for(i=0; i<nmatint; i++) {
    for(j=0; j<nreg; j++) {
      j1=i*(nreg+1)+j;
      bfirst=counteq[j1+begblock[rank]];
      bend=counteqnoadd[j1+begblock[rank]]+counteq[j1+begblock[rank]];
      ISCreateStride(PETSC_COMM_SELF,bend-bfirst,bfirst,1,rowindices+j2);
      bfirst=countvarintra1[j1+begblock[rank]];
      bend=countvarintra1[j1+1+begblock[rank]];
      ISCreateStride(PETSC_COMM_SELF,bend-bfirst,bfirst,1,colindices+j2);
      j2++;
    }
  }
  ierr = MatCreateSubMatrices(A,nmatint*nreg,rowindices,colindices,MAT_INITIAL_MATRIX,&submatA);
  CHKERRQ(ierr);
  j2=0;
  for (i=0; i<nmatint; i++) {
    for(j=0; j<nreg; j++) {
      ierr = ISDestroy(&rowindices[j2]);
      CHKERRQ(ierr);
      ierr = ISDestroy(&colindices[j2]);
      CHKERRQ(ierr);
      j2++;
    }
  }
  PetscFree(rowindices);
  PetscFree(colindices);
  MPI_Barrier(PETSC_COMM_WORLD);
  int *ha_rows2= (int *) calloc (VecSize,sizeof(int));
  int *ha_cols2= (int *) calloc (VecSize,sizeof(int));
  if(cntl6==0&&SORD==0)cntl6in=1e-5;
  else cntl6in=cntl6;
  #pragma omp parallel private(j1,j2,j3,j4,j5,j6,j7,j8,ai,aj,vals,nz,nrow,ncol,nz1,i,j,bfirst,lasize) shared(cntl6in,ha_rows2,ha_cols2,ha_ndblocks1,counteqnoadd1,countvarintra2,countvarintra1,submatA)
  {
  long int lasizemax=0,nrowmax=0,ncolmax=0;
  #pragma omp for schedule (static)
  for(j1=0; j1<nmatint; j1++) {
    j8=j1*nreg;
    for(j2=0; j2<nreg; j2++) {
      Mat_SeqAIJ         *aa=(Mat_SeqAIJ*)submatA[j8]->data;//*aa=subA->data;
      nz=aa->nz;
      nrow=submatA[j8]->rmap->n;
      if(nrowmax<nrow)nrowmax=nrow;
      ncol=submatA[j8]->cmap->n;
      if(ncolmax<ncol)ncolmax=ncol;
      nz1=nz;
      if(nz1<nrow)nz1=nrow;
      if(nz1<ncol)nz1=ncol;
      lasize=ceil((laA/100.0)*nz1);
      if(lasizemax<lasize)lasizemax=lasize;
      j8++;
    }
  }
  lasizemax+=10;
    int *irn=(int *) calloc (lasizemax,sizeof(int));
    int *jcn=(int *) calloc (lasizemax,sizeof(int));
    int*irn1=NULL,*jcn1=NULL;
  irn1=realloc(irn1,nrowmax*sizeof(int));
  jcn1=realloc(jcn1,ncolmax*sizeof(int));
    ha_cgetype *values= (ha_cgetype *) calloc (lasizemax,sizeof(ha_cgetype));
  int *keep=(int *) malloc ((nrowmax+5*ncolmax+4*ncolmax+7)*sizeof(int));//ICNTL(6)=1
  ha_cgetype *w51;
  if(nrowmax>ncolmax) w51=(ha_cgetype *) malloc ((5*nrowmax)*sizeof(ha_cgetype));
  else w51=(ha_cgetype *) malloc ((5*ncolmax)*sizeof(ha_cgetype));
  int *iw51=(int *) malloc ((6*nrowmax+3*ncolmax)*sizeof(int));
    int *insize=(int *) calloc (6,sizeof(int));
  #pragma omp for schedule (static)
  for(j1=0; j1<nmatint; j1++) {
    j4=0;
    j5=0;
    j6=0;
    j7=0;
    j8=j1*nreg;
    for(j2=0; j2<nreg; j2++) {
      j3=j1*(nreg+1)+j2;
      Mat_SeqAIJ         *aa=(Mat_SeqAIJ*)submatA[j8]->data;//*aa=subA->data;
      ai= aa->i;
      aj= aa->j;
      vals=aa->a;
      nz=aa->nz;
      nrow=submatA[j8]->rmap->n;
      ncol=submatA[j8]->cmap->n;
      nz1=nz;
      if(nz1<nrow)nz1=nrow;
      if(nz1<ncol)nz1=ncol;
      lasize=ceil((laA/100.0)*nz1);
      insize[0]=nrow;
      insize[1]=ncol;
      insize[2]=nz;
      insize[4]=laA;
      insize[5]=lasize;
      for(i=0; i<nrow-1; i++) {
        for(j=ai[i]; j<ai[i+1]; j++) {
          irn[j]=i+1;
          jcn[j]=aj[j]+1;
          values[j]=vals[j];
        }
      }
      for(j=ai[nrow-1]; j<nz; j++) {
        irn[j]=nrow;
        jcn[j]=aj[j]+1;
        values[j]=vals[j];
      }
      spec51m_rank_(insize,&cntl6in,irn,jcn,values,irn1,jcn1,keep,w51,iw51);
      bfirst=counteq[j3+begblock[rank]];
      for(i=0; i<nrow; i++) {
        ha_rows2[i+counteq[j3+begblock[rank]]]=irn1[i]-1+bfirst;
      }
      bfirst=countvarintra1[j3+begblock[rank]];
      for(i=0; i<ncol; i++) {
        ha_cols2[i+countvarintra1[j3+begblock[rank]]]=jcn1[i]-1+bfirst;
      }
      printf("rank %d mat rank %d nrow %d ncol %d j3 %d proc %d\n",rank,insize[3],nrow,ncol,j3,rank);
      ha_ndblocks1[j3+begblock[rank]]=insize[3];
      counteqnoadd1[j3+begblock[rank]]=insize[3];
      countvarintra2[j3+begblock[rank]]=insize[3];
      j4+=counteqnoadd1[j3+begblock[rank]];
      j5+=countvarintra2[j3+begblock[rank]];
      j6+=counteqnoadd[j3+begblock[rank]];
      j7+=(countvarintra1[j3+1+begblock[rank]]-countvarintra1[j3+begblock[rank]]);
      MatDestroy(&submatA[j8]);
      j8++;
    }
    counteqnoadd1[j1*(nreg+1)+nreg+begblock[rank]]=j6-j4+counteqnoadd[j1*(nreg+1)+nreg+begblock[rank]];
    countvarintra2[j1*(nreg+1)+nreg+begblock[rank]]=j7-j5+countvarintra1[j1*(nreg+1)+nreg+1+begblock[rank]]-countvarintra1[j1*(nreg+1)+nreg+begblock[rank]];
  }
      free(insize);
      free(irn);
      free(jcn);
      free(irn1);
      free(jcn1);
    free(keep);
    free(w51);
    free(iw51);
      free(values);
  }
  PetscFree(submatA);
  MPI_Allreduce(counteqnoadd1,counteqnoadd2,ntime*(nreg+1), MPI_LONG, MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(countvarintra2,countvarintra6,ntime*(nreg+1), MPI_LONG, MPI_SUM,PETSC_COMM_WORLD);
  counteq2[0]=0;
  countvarintra4[0]=0;
  for(i=1; i<ndblock; i++) {
    counteq2[i]=counteq2[i-1]+counteqnoadd2[i-1];
    countvarintra4[i]=countvarintra4[i-1]+countvarintra6[i-1];
  }
  counteq2[i]=VecSize;//Not change
  countvarintra4[i]=countvarintra4[i-1]+countvarintra6[i-1];

  int *ha_rows1= (int *) calloc (VecSize,sizeof(int));
  int *ha_cols1= (int *) calloc (VecSize,sizeof(int));
  for(j1=0; j1<nmatint; j1++) {
    for(j2=0; j2<nreg; j2++) {
      j3=j1*(nreg+1)+j2;
      bfirst=counteq2[j3+begblock[rank]];
      for(i=bfirst; i<counteq2[j3+begblock[rank]]+counteqnoadd2[j3+begblock[rank]]; i++) {
        ha_rows1[i]=ha_rows2[i-bfirst+counteq[j3+begblock[rank]]];
      }
      bfirst=countvarintra4[j3+begblock[rank]];
      for(i=bfirst; i<countvarintra4[j3+1+begblock[rank]]; i++) {
        ha_cols1[i]=ha_cols2[i-bfirst+countvarintra1[j3+begblock[rank]]];
      }
    }
    j3=j1*(nreg+1)+j2;
    j4=counteq2[j3+begblock[rank]];
    j5=countvarintra4[j3+begblock[rank]];
    for(j2=0; j2<nreg; j2++) {
      j6=j1*(nreg+1)+j2;
      for(i=counteqnoadd[j6+begblock[rank]]+counteq[j6+begblock[rank]]-1; i>counteqnoadd2[j6+begblock[rank]]+counteq[j6+begblock[rank]]-1; i--) {
        ha_rows1[j4]=ha_rows2[i];
        j4++;
      }
      for(i=countvarintra1[j6+1+begblock[rank]]-1; i>countvarintra1[j6+begblock[rank]]+countvarintra6[j6+begblock[rank]]-1; i--) {
        ha_cols1[j5]=ha_cols2[i];
        j5++;
      }
    }
    for(i=j4; i<j4+counteqnoadd[j3+begblock[rank]]; i++) {
      ha_rows1[i]=i-j4+counteq[j3+begblock[rank]];
    }
    for(i=j5; i<j5+countvarintra1[j3+1+begblock[rank]]-countvarintra1[j3+begblock[rank]]; i++) {
      ha_cols1[i]=i-j5+countvarintra1[j3+begblock[rank]];
    }
  }
  free(ha_cols2);
  free(ha_rows2);
  for(j1=0; j1<nmatint; j1++) {
    j3=j1*(nreg+1)+nreg;
    nrow=counteqnoadd2[j3+begblock[rank]];
    ncol=countvarintra6[j3+begblock[rank]];
    j4=j1+begblock[rank]/(nreg+1);

    if(nrow<ncol)ndbbdrank[j4]=nrow;
    else ndbbdrank[j4]=ncol;

    if(ndbbdrank[j4]>0) {
      ha_ndblocks1[j3+begblock[rank]]=ndbbdrank[j4];
    }

  }
  memcpy(counteq,counteq2,(ndblock+1)*sizeof(uvadd));
  memcpy(counteqnoadd,counteqnoadd2,(ndblock)*sizeof(uvadd));
  memcpy(countvarintra1,countvarintra4,(ndblock+1)*sizeof(uvadd));
  MPI_Barrier(PETSC_COMM_WORLD);
  MPI_Allreduce(ha_cols1,ha_cols,VecSize, MPI_INT, MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(ha_rows1,ha_rows,VecSize, MPI_INT, MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(ha_ndblocks1,ha_ndblocks,ndblock, MPI_INT, MPI_SUM,PETSC_COMM_WORLD);
  free(ha_cols1);
  free(ha_rows1);
  free(ha_ndblocks1);
  free(begblock);
  free(counteq2);
  free(counteqnoadd1);
  free(counteqnoadd2);
  free(countvarintra2);
  free(countvarintra4);
  free(countvarintra6);
  return 1;
}

int HaNDBBDMatOder(Mat A, uvadd VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int nreg, int ntime, uvadd nvarele, uvadd *ha_eqadd,int *ha_rows,int *ha_cols, uvadd ndblock,int *ha_ndblocks, uvadd *countvarintra1, uvadd *counteq, uvadd *counteqnoadd,uvdim laA,uvdim laDi,ha_cgetype cntl6,PetscInt* ndbbdrank,PetscBool presol) {
  FILE *presolfile;
  char j1name[1024],filename[1024],rankname[1024];
  size_t frd;
  PetscInt bfirst,bend,sblockin,nmatin,nmatint,nmatinplus,nrowcolin,sumrowcolin,rfirst,cfirst;
  PetscInt i,j,j0,j1,j2,j3,j4,j5,j6,j7,j8,nrow,ncol,nz,nz1,*ai,*aj,la;
  if(rank<10)strcpy(rankname,"000");
  if(rank<100&&rank>9)strcpy(rankname,"00");
  if(rank<1000&&rank>99)strcpy(rankname,"0");
  if(rank>=1000)rankname[0]='\0';
  sprintf(filename, "%d",rank);
  strcat(rankname,filename);
  //*****New order test**********//
  nmatint=(uvadd)ntime/mpisize;
  for(i=0; i<mpisize; i++)if(rank+1<=ntime-mpisize*nmatint)nmatint++;
  nmatin=(nreg+1)*nmatint;
  int *begblock= (int *) calloc (mpisize,sizeof(int));
  begblock[rank]=nmatin;
  printf("rank %d nmatin %d nmatint %d ndblock %ld\n",rank,nmatin,nmatint,ndblock);
  for(i=0; i<mpisize; i++) {
    j=rank;
    j2=j;
    MPI_Bcast(&j2,1, MPI_INT,i, PETSC_COMM_WORLD);
    j=j2;
    MPI_Bcast(&begblock[j2],1, MPI_INT,i, PETSC_COMM_WORLD);
  }
  j2=begblock[0];
  begblock[0]=0;
  for(i=1; i<mpisize; i++) {
    j1=begblock[i];
    begblock[i]=begblock[i-1]+j2;
    j2=j1;
  }

  int *ha_rows3= (int *) calloc (VecSize,sizeof(int));
  int *ha_cols3= (int *) calloc (VecSize,sizeof(int));
  for(i=0; i<nmatin; i++) {
    bfirst=counteq[i+begblock[rank]];
    bend=counteqnoadd[i+begblock[rank]]+counteq[i+begblock[rank]];
    for(j=bfirst; j<bend; j++)ha_rows3[j]=ha_rows[j];
    bfirst=countvarintra1[i+begblock[rank]];
    bend=countvarintra1[i+1+begblock[rank]];
    for(j=bfirst; j<bend; j++)ha_cols3[j]=ha_cols[j];
  }
  int *ha_ndblocks2= (int *) calloc (ndblock,sizeof(int));
  for(j1=0; j1<nmatin; j1++)ha_ndblocks2[j1+begblock[rank]]=ha_ndblocks[j1+begblock[rank]];

  for(j1=0; j1<nmatint; j1++) {
    j3=j1*(nreg+1)+nreg;
    j4=j1+begblock[rank]/(nreg+1);



    if(!presol) {
      int *insized=(int *) calloc (5,sizeof(int));
      if(j1<10)strcpy(j1name,"000");
      if(j1<100&&j1>9)strcpy(j1name,"00");
      if(j1<1000&&j1>99)strcpy(j1name,"0");
      if(j1>=1000)j1name[0]='\0';
      sprintf(filename, "%d",j1);
      strcat(j1name,filename);
      strcpy(filename,"_rank");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "r"))==NULL) {
        printf("Cannot open file. Please run presol by setting -presol 1!\n");
      }
      frd=fread(insized, sizeof(int), 5, presolfile);
      fclose(presolfile);
      if( frd== 0)printf("File read error. No presol, take min rank!\n");
      if(ndbbdrank[j4]>insized[3]) {
        int *irn1=(int *) calloc (insized[0],sizeof(int));
        int *jcn1=(int *) calloc (insized[1],sizeof(int));
        strcpy(filename,"_row");
        strcat(filename,rankname);
        strcat(filename,j1name);
        strcat(filename,".bin");
        if((presolfile=fopen(filename, "r"))==NULL) {
          printf("Cannot open file.\n");
        }
        frd=fread(irn1, sizeof(int), insized[0], presolfile);
        if(frd == 0)printf("File read error.");
        fclose(presolfile);
        strcpy(filename,"_col");
        strcat(filename,rankname);
        strcat(filename,j1name);
        strcat(filename,".bin");
        if((presolfile=fopen(filename, "r"))==NULL) {
          printf("Cannot open file.\n");
        }
        frd=fread(jcn1, sizeof(int), insized[0], presolfile);
        if( frd== 0) printf("File read error.");
        fclose(presolfile);
        int *indices= (int *) calloc (insized[0],sizeof(int));
        for(i=0; i<insized[0]; i++) {
          indices[i]=ha_rows3[counteq[j3+begblock[rank]]+irn1[i]-1];
        }
        memcpy(&ha_rows3[counteq[j3+begblock[rank]]],indices,insized[0]*sizeof(int));
        for(i=0; i<insized[0]; i++) {
          indices[i]=ha_cols3[countvarintra1[j3+begblock[rank]]+jcn1[i]-1];
        }
        memcpy(&ha_cols3[countvarintra1[j3+begblock[rank]]],indices,insized[0]*sizeof(int));
        ndbbdrank[j4]=insized[3];
        free(indices);
        free(irn1);
        free(jcn1);
      }
      free(insized);
    }

    if(ndbbdrank[j4]>0) {
      ha_ndblocks2[j3+begblock[rank]]=ndbbdrank[j4];
    }

  }
  MPI_Barrier(PETSC_COMM_WORLD);
  MPI_Allreduce(ha_cols3,ha_cols,VecSize, MPI_INT, MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(ha_rows3,ha_rows,VecSize, MPI_INT, MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(ha_ndblocks2,ha_ndblocks,ndblock, MPI_INT, MPI_SUM,PETSC_COMM_WORLD);
  printf("OK\n");
  free(ha_cols3);
  free(ha_rows3);
  free(ha_ndblocks2);
  free(begblock);
  return 1;
}


int HaDBBDParSol(Mat A, Vec b, ha_cgetype *x1, uvadd VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int *ha_rows,int *ha_cols, uvadd ndblock,int *ha_ndblocks, uvadd *countvarintra1, uvadd *counteq, uvadd *counteqnoadd,uvdim laA,uvdim laD,PetscReal cntl3) {//,bool iter
  IS *rowindices,*colindices,*Cindices,*Bindices,*BBindices,*BBiindices;
  const PetscInt *nindices;
  PetscInt bfirst,bend,sblockin,nmatin,nmatinplus,nrowcolin,sumrowcolin,i,i1,j,j0,j1,j2,j3,j4,j5,j6,l0,l1,l2,l3,l4,l5,rank1,proc1=0,nnzmax,j1nz,j1irnbs;
  Mat *submatA,*submatC,*submatB,*submatD,submatBT,submatCT,*submatBB;
  Vec vecxd;
  PetscInt *ai,*aj,*aic,*ajc;
  PetscInt nrow,ncol,nz,nrowc,ncolc,nrowb,ncolb,nzc,nzv;
  uvadd lasize,ldsize;
  MPI_Status   status;
  clock_t timestr,timeend;
  int j2int,la1;
  size_t freadresult;
  PetscScalar *vals,*valsc,vecval;
  ha_cgetype vval;
  PetscErrorCode ierr;
  PetscViewer viewer;
  FILE* fp1,*fp2,*fp3;
  int *begblock= (int *) calloc (mpisize,sizeof(int));
  char filename[1024],rankname[1024],j1name[1024],tempchar[1024];//,fn1[1024],fn2[1024],fn3[1024]
  if(rank<10)strcpy(rankname,"000");
  if(rank<100&&rank>9)strcpy(rankname,"00");
  if(rank<1000&&rank>99)strcpy(rankname,"0");
  if(rank>=1000)rankname[0]='\0';
  sprintf(filename, "%d",rank);
  strcat(rankname,filename);
  nmatin=(uvadd)ndblock/mpisize;
  int nmatinBB=1;
  if(nmatin*mpisize<ndblock)nmatinplus=nmatin+1;
  else nmatinplus=nmatin;
  for(i=0; i<mpisize; i++)if(rank+1<=ndblock-mpisize*nmatin)nmatin++;
  begblock[rank]=nmatin;
  for(i=0; i<mpisize; i++) {
    j=rank;
    j2=j;
    MPI_Bcast(&j2,1, MPI_INT,i, PETSC_COMM_WORLD);
    j=j2;
    MPI_Bcast(&begblock[j2],1, MPI_INT,i, PETSC_COMM_WORLD);
  }
  j2=begblock[0];
  begblock[0]=0;
  for(i=1; i<mpisize; i++) {
    j1=begblock[i];
    begblock[i]=begblock[i-1]+j2;
    j2=j1;
  }
  ierr = PetscMalloc(nmatin*sizeof(IS **),&rowindices);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatin*sizeof(IS **),&colindices);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatin*sizeof(IS **),&Cindices);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatin*sizeof(IS **),&Bindices);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatinBB*sizeof(IS **),&BBindices);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatin*sizeof(Mat *),&submatA);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatin*sizeof(Mat *),&submatC);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatin*sizeof(Mat *),&submatB);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatinBB*sizeof(Mat *),&submatBB);
  CHKERRQ(ierr);
  PetscScalar **yi1= (PetscScalar**)calloc(nmatin,sizeof(PetscScalar*));
  for (i=0; i<nmatin; i++) yi1[i] = (PetscScalar*)calloc(ha_ndblocks[i+begblock[rank]],sizeof(PetscScalar));
  sumrowcolin=0;
  for(j=0; j<ndblock; j++)sumrowcolin+=ha_ndblocks[j];
  PetscInt *indicesC= (PetscInt *) calloc (VecSize-sumrowcolin,sizeof(PetscInt));
  PetscInt *indicesB= (PetscInt *) calloc (VecSize-sumrowcolin,sizeof(PetscInt));
  PetscInt *offblock= (PetscInt *) calloc (ndblock+1,sizeof(PetscInt));
  PetscInt *offblockrow= (PetscInt *) calloc (ndblock+1,sizeof(PetscInt));
  ha_cgetype *vecbiui= (ha_cgetype *) calloc (VecSize-sumrowcolin,sizeof(ha_cgetype));
  offblock[0]=0;
  for(j=0; j<ndblock; j++)offblock[j+1]=offblock[j]+countvarintra1[j+1]-countvarintra1[j]-ha_ndblocks[j];
  offblockrow[0]=0;
  for(j=0; j<ndblock; j++)offblockrow[j+1]=offblockrow[j]+counteqnoadd[j]-ha_ndblocks[j];

  char **fn01= (char**)calloc(nmatinplus,sizeof(char*));
  for (i=0; i<nmatinplus; i++) fn01[i] = (char*)calloc(1024,sizeof(char));
  char **fn02= (char**)calloc(nmatinplus,sizeof(char*));
  for (i=0; i<nmatinplus; i++) fn02[i] = (char*)calloc(1024,sizeof(char));
  char **fn03= (char**)calloc(nmatinplus,sizeof(char*));
  for (i=0; i<nmatinplus; i++) fn03[i] = (char*)calloc(1024,sizeof(char));

  timestr=clock();
  for(j1=0; j1<nmatinplus; j1++) {
    if(j1<nmatin) {
      if(j1<10)strcpy(j1name,"000");
      if(j1<100&&j1>9)strcpy(j1name,"00");
      if(j1<1000&&j1>99)strcpy(j1name,"0");
      if(j1>=1000)j1name[0]='\0';
      sprintf(filename, "%d",j1);
      strcat(j1name,filename);
      strcpy(filename,"_irnv");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      strcpy(fn01[j1],filename);
      strcpy(filename,"_keep");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      strcpy(fn02[j1],filename);
      strcpy(filename,"_vav");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      strcpy(fn03[j1],filename);
    }
  }
  printf("OK012345 rank %d mpis %d\n",rank,mpisize);
  for(j=0; j<mpisize; j++) {
    printf("j %d\n",j);
    for(j3=0; j3<nmatinplus; j3++) {
      if(j3==nmatin)i=j3-1;
      else i=j3;
      bfirst=counteq[i+begblock[rank]];
      bend=ha_ndblocks[i+begblock[rank]];
      MPI_Bcast(&bfirst,1, MPI_INT,j, PETSC_COMM_WORLD);
      MPI_Bcast(&bend,1, MPI_INT,j, PETSC_COMM_WORLD);
      PetscInt *indices= (PetscInt *) calloc (bend,sizeof(PetscInt));
      ha_cgetype *yi0 = (ha_cgetype*)calloc(bend,sizeof(ha_cgetype));
      for(j1=0; j1<bend; j1++) { //ha_ndblocks[i+begblock[j]]
        j2=ha_rows[bfirst+j1]+bfirst;
        if(j2>=Istart&&j2<Iend)indices[j1]=j2;
        else indices[j1]=-1;
      }
      VecGetValues(b,bend,indices,yi0);

      HaReduce(yi0,(forint)bend,mpisize,rank,j);
      if(rank==j&&i<nmatin) {
        memcpy (yi1[i],yi0,bend*sizeof(ha_cgetype));
      }
      free(indices);
      free(yi0);
    }
  }
  printf("OK12345 rank %d\n",rank);
  j1=0;
  for(i=0; i<nmatin; i++) {
    bfirst=counteq[i+begblock[rank]];
    PetscInt *indices= (PetscInt *) calloc (ha_ndblocks[i+begblock[rank]],sizeof(PetscInt));
    for(j=0; j<ha_ndblocks[i+begblock[rank]]; j++) {
      indices[j]=ha_rows[bfirst+j]+bfirst;
    }

    bend=counteq[i+begblock[rank]]+counteqnoadd[i+begblock[rank]];
    for(j=offblockrow[i+begblock[rank]]; j<offblockrow[i+1+begblock[rank]]; j++) {
      j1=j-offblockrow[i+begblock[rank]];
      indicesB[j]=ha_rows[bend-j1-1]+bfirst;
    }
    ISCreateGeneral(PETSC_COMM_SELF,ha_ndblocks[i+begblock[rank]],indices,PETSC_COPY_VALUES,rowindices+i);
    bfirst=countvarintra1[i+begblock[rank]];
    for(j=0; j<ha_ndblocks[i+begblock[rank]]; j++) {
      indices[j]=ha_cols[bfirst+j]+bfirst;
    }
    bend=countvarintra1[i+1+begblock[rank]];
    for(j=offblock[i+begblock[rank]]; j<offblock[i+1+begblock[rank]]; j++) {
      j1=j-offblock[i+begblock[rank]];
      indicesC[j]=ha_cols[bend-j1-1]+bfirst;
    }

    ISCreateGeneral(PETSC_COMM_SELF,ha_ndblocks[i+begblock[rank]],indices,PETSC_COPY_VALUES,colindices+i);
    free(indices);
  }
  printf("0okkkk\n");
  j1=offblock[ndblock];
  for(i=countvarintra1[ndblock]; i<VecSize; i++) {
    indicesC[j1]=i;
    j1++;
  }
  printf("1okkkk\n");
  j1=offblockrow[ndblock];
  for(i=counteq[ndblock-1]+counteqnoadd[ndblock-1]; i<VecSize; i++) {
    indicesB[j1]=i;
    j1++;
  }
  printf("2okkkk\n");
  for(i=0; i<mpisize; i++) {
    j=offblockrow[begblock[rank]];
    j1=offblockrow[nmatin+begblock[rank]]-offblockrow[begblock[rank]];
    MPI_Bcast(&j,1, MPI_INT,i, PETSC_COMM_WORLD);
    MPI_Bcast(&j1,1, MPI_INT,i, PETSC_COMM_WORLD);
    MPI_Bcast(&indicesB[j],j1, MPI_INT,i, PETSC_COMM_WORLD);
  }
  printf("rank1 %d\n",rank);

  for(i=0; i<VecSize-sumrowcolin; i++)if(indicesB[i]>=Istart) {
      j=i;
      break;
    }
  for(j2=i; j2<VecSize-sumrowcolin; j2++)if(indicesB[j2]>=Iend) {
      j1=j2;
      break;
    }
  if(j2==VecSize-sumrowcolin)j1=VecSize-sumrowcolin;
  printf("rank %d j %d j1 %d istart %d iend %d\n",rank,j,j1,Istart,Iend);
  VecGetValues(b,j1-j,&indicesB[j],&vecbiui[j]);//implicite yd save mem
  ierr = VecDestroy(&b);
  CHKERRQ(ierr);
  printf("rank2 %d\n",rank);
  for(i=0; i<mpisize; i++) {
    j=offblock[begblock[rank]];
    j1=offblock[nmatin+begblock[rank]]-offblock[begblock[rank]];
    MPI_Bcast(&j,1, MPI_INT,i, PETSC_COMM_WORLD);
    MPI_Bcast(&j1,1, MPI_INT,i, PETSC_COMM_WORLD);
    MPI_Bcast(&indicesC[j],j1, MPI_INT,i, PETSC_COMM_WORLD);
  }
  printf("OK132\n");
  for(i=0; i<nmatin; i++) {
    ISCreateGeneral(PETSC_COMM_SELF,VecSize-sumrowcolin,indicesC,PETSC_COPY_VALUES,Cindices+i);
    ISCreateGeneral(PETSC_COMM_SELF,VecSize-sumrowcolin,indicesB,PETSC_COPY_VALUES,Bindices+i);
  }
  free(indicesB);
  free(offblock);
  free(offblockrow);
  PetscInt *indexBB= (PetscInt *) calloc (VecSize,sizeof(PetscInt));
  for(i=0; i<VecSize; i++)indexBB[i]=i;
  ISCreateGeneral(PETSC_COMM_SELF,VecSize,indexBB,PETSC_COPY_VALUES,BBindices);
  printf("Begin partitioning rank %d\n",rank);
  ierr = MatCreateSubMatrices(A,nmatin,rowindices,colindices,MAT_INITIAL_MATRIX,&submatA);
  printf("A\n");
  ierr = MatCreateSubMatrices(A,nmatin,rowindices,Cindices,MAT_INITIAL_MATRIX,&submatC);
  printf("C\n");
  ierr = MatCreateSubMatrices(A,nmatinBB,Bindices,BBindices,MAT_INITIAL_MATRIX,&submatBB);
  ierr = MatDestroy(&A);
  CHKERRQ(ierr);
  PetscInt BBrow,BBcol,CCrow,CCcol;
  MatGetSize(submatBB[0],&BBrow,&BBcol);
  MatGetSize(submatC[0],&CCrow,&CCcol);
  printf("bbrow %d bbcol %d ccrow %d cccol %d\n",BBrow,BBcol,BBcol,CCcol);
  PetscInt *indexBBi= (PetscInt *) calloc (BBrow,sizeof(PetscInt));
  for(i=0; i<BBrow; i++)indexBBi[i]=i;
  for(i=0; i<nmatin; i++) {
    ISCreateGeneral(PETSC_COMM_SELF,BBrow,indexBBi,PETSC_COPY_VALUES,Bindices+i);
  }
  ierr = MatCreateSubMatrices(submatBB[0],nmatin,Bindices,colindices,MAT_INITIAL_MATRIX,&submatB);
  printf("B\n");
  if(rank==mpisize-1)proc1=1;
  ierr = PetscMalloc(proc1*sizeof(Mat *),&submatD);
  CHKERRQ(ierr);
  printf("D\n");
  MatCreateSubMatrices(submatBB[0],proc1,Bindices,Cindices,MAT_INITIAL_MATRIX,&submatD);
  printf("End partitionning rank %d time %f\n",rank,((double)clock()-timestr)/CLOCKS_PER_SEC);
  for (i=0; i<nmatin; i++) {
    ierr = ISDestroy(&rowindices[i]);
    CHKERRQ(ierr);
    ierr = ISDestroy(&Cindices[i]);
    CHKERRQ(ierr);
    ierr = ISDestroy(&Bindices[i]);
    CHKERRQ(ierr);
  }
  ierr = ISDestroy(&BBindices[0]);
  PetscFree(rowindices);
  PetscFree(Cindices);
  PetscFree(Bindices);
  PetscFree(BBindices);
  PetscFree(indexBB);
  PetscFree(indexBBi);
  ierr = MatDestroy(&submatBB[0]);
  PetscFree(submatBB);
  CHKERRQ(ierr);
  printf("OK1 rank %d sumrowcolin %d\n",rank,sumrowcolin);
  ha_cgetype *xi1 = (ha_cgetype*)calloc(sumrowcolin,sizeof(ha_cgetype));
  int insizes=17;
  int *insize=(int *) calloc (insizes*nmatin,sizeof(int));
  int *bivinzrow=(int *) calloc (VecSize-sumrowcolin,sizeof(int));
  int *bivinzcol=(int *) calloc (VecSize-sumrowcolin,sizeof(int));
  PetscReal *ccolnorms=(PetscReal *) calloc (BBrow,sizeof(PetscReal));
  for(j1=0; j1<nmatinplus; j1++) {
    if(j1<nmatin) {
      MatGetColumnNorms(submatC[j1],NORM_1,ccolnorms);
      for(i=0; i<BBrow; i++)if(ccolnorms[i]>0)bivinzcol[i]++;
      Mat_SeqAIJ         *ab=(Mat_SeqAIJ*)submatB[j1]->data;//*aa=subA->data;
      ai= ab->i;
      nz=ab->nz;
      nrow=submatB[j1]->rmap->n;
      for(i=0; i<nrow-1; i++)if(ai[i]!=ai[i+1])bivinzrow[i]++;
      if(ai[i]<nz)bivinzrow[i]++;
    }
  }
  free(ccolnorms);
  int bivirowsize=1,bivicolsize=1;
  long int vecbivisize=0,li,lj;
  for(i=0; i<VecSize-sumrowcolin; i++) {
    if(bivinzrow[i]>0) {
      bivinzrow[i]=bivirowsize;
      bivirowsize++;
    }
    if(bivinzcol[i]>0) {
      bivinzcol[i]=bivicolsize;
      bivicolsize++;
    }
  }
  bivirowsize--;
  bivicolsize--;
  long int *bivinzrow1=(long int *) calloc (bivirowsize,sizeof(long int));
  long int *bivinzcol1=(long int *) calloc (bivicolsize,sizeof(long int));
  j=0;
  j2=0;
  for(li=0; li<VecSize-sumrowcolin; li++) {
    if(bivinzrow[li]>0) {
      bivinzrow1[j]=li*BBrow;
      j++;
    }
    if(bivinzcol[li]>0) {
      bivinzcol1[j2]=li;
      j2++;
    }
    bivinzrow[li]=(bivinzrow[li]-1)*bivicolsize;//faster MSOL
  }
  vecbivisize=((long int)bivirowsize)*((long int)bivicolsize);
  if(vecbivisize<0)printf("Overflow detected, the interface problem is too big!!!!!!\n");
  ha_cgetype *vecbivi= (ha_cgetype *) calloc (vecbivisize,sizeof(ha_cgetype));
  printf("vecbivisize %ld rank %d\n",vecbivisize,rank);

  ha_cgetype *xi1point;
  uvadd xi1indx=0;
  int jthrd,nthrd=1;
  nthrd=omp_get_max_threads();
  uvadd *nthrds= (uvadd *) calloc (nthrd+1,sizeof(uvadd));
  uvadd *nthrds1= (uvadd *) calloc (nthrd+1,sizeof(uvadd));
  #pragma omp parallel private(jthrd,j1,j,timestr,submatCT,ai,aj,vals,nz,nrow,ncol,lasize,aic,ajc,valsc,nzc,nrowc,ncolc,nrowb,ncolb,xi1point,xi1indx,i,timeend) shared(submatA,insize,nthrds,nthrds1)
  {
  jthrd=omp_get_thread_num();
  int windx=0,bindx,eindx;
  windx=90/nthrd;
  bindx=windx*jthrd;
  eindx=bindx+windx;
  windx=bindx;
  #pragma omp for
  for(j1=0; j1<nmatinplus; j1++){
    if(j1<nmatin) {
      nthrds[jthrd+1]+=ha_ndblocks[j1+begblock[rank]];
    }
    nthrds1[jthrd+1]+=1;
  }
  if(jthrd==0){
  for(j1=1; j1<nthrd+1; j1++){
    nthrds[j1]+=nthrds[j1-1];
    nthrds1[j1]+=nthrds1[j1-1];
  }
  }
  #pragma omp barrier
  xi1indx=nthrds[jthrd];
  for(j1=nthrds1[jthrd]; j1<nthrds1[jthrd+1]; j1++) {
    if(j1<nmatin) {
      time(&timestr);//=clock();
      Mat_SeqAIJ         *aa=(Mat_SeqAIJ*)submatA[j1]->data;//*aa=subA->data;
      MatTranspose(submatC[j1],MAT_INITIAL_MATRIX,&submatCT);
      Mat_SeqAIJ         *ac=(Mat_SeqAIJ*)submatCT->data;//*aa=subA->data;
      ai= aa->i;
      aj= aa->j;
      vals=aa->a;
      nz=aa->nz;
      nrow=submatA[j1]->rmap->n;
      ncol=submatA[j1]->cmap->n;
      lasize=ceil((laA/100.0)*nz);
      int *irn=(int *) calloc (lasize,sizeof(int));
      int *jcn=(int *) calloc (lasize,sizeof(int));
      ha_cgetype *values= (ha_cgetype *) calloc (lasize,sizeof(ha_cgetype));
      insize[j1*insizes+13]=bivirowsize;
      insize[j1*insizes+14]=bivicolsize;
      insize[j1*insizes]=nrow;
      insize[j1*insizes+1]=ncol;
      insize[j1*insizes+2]=nz;
      aic= ac->i;
      ajc= ac->j;
      valsc=ac->a;
      nzc=ac->nz;
      nrowc=submatCT->rmap->n;
      ncolc=submatCT->cmap->n;
      insize[j1*insizes+3]=nrowc;
      insize[j1*insizes+4]=ncolc;
      memcpy (irn+nz,ai,(nrow+1)*sizeof(PetscInt));
      memcpy (jcn,aj,nz*sizeof(PetscInt));
      memcpy (values,vals,nz*sizeof(ha_cgetype));
      insize[j1*insizes+2]=nz;
      insize[j1*insizes+5]=nzc;
      MatDestroy(&submatA[j1]);
      insize[j1*insizes+9]=laA;
      insize[j1*insizes+16]=lasize;
      insize[j1*insizes+10]=rank;
      insize[j1*insizes+11]=j1;
      nrowb=submatB[j1]->rmap->n;
      ncolb=submatB[j1]->cmap->n;
      Mat_SeqAIJ         *ab=(Mat_SeqAIJ*)submatB[j1]->data;//*aa=subA->data;
      ai= ab->i;//irnbs+j1irnbs;
      aj= ab->j;//jcnbs+j1nz;
      vals=ab->a;//valbs+j1nz;
      nz=ab->nz;//nzbs[j1];
      insize[j1*insizes+6]=nrowb;//ncolb;
      insize[j1*insizes+7]=ncolb;//nrowb;
      insize[j1*insizes+8]=nz;
      insize[j1*insizes+15]=7+windx;//7+rank+jthrd+windx;//(j4+jthrd)%90+7;
      windx++;
      if(windx==eindx)windx=bindx;
      xi1point=xi1+xi1indx;
      spec48m_msol_(insize+j1*insizes,irn,jcn,values,yi1[j1],xi1point,aic,ajc,valsc,ai,aj,vals,vecbivi,bivinzrow,bivinzcol);
      MatDestroy(&submatCT);
      free(yi1[j1]);
      free(irn);
      free(jcn);
      free(values);
      //Multiply Bi by ui:
      for(i=0; i<nrowb-1; i++) {
        for(j=ai[i]; j<ai[i+1]; j++) {
          #pragma omp atomic
          vecbiui[i]-=vals[j]*xi1point[aj[j]];//xi[i];
        }
      }
      for(j=ai[i]; j<nz; j++) {
        #pragma omp atomic
        vecbiui[i]-=vals[j]*xi1point[aj[j]];//xi[i];
      }
      MatDestroy(&submatB[j1]);//submatBT);
      xi1indx+=ha_ndblocks[j1+begblock[rank]];
      time(&timeend);
      printf("Submatrix %d rank %d thrd %d calculation time %f\n",j1,rank,jthrd,difftime(timeend,timestr));
    }
  }
  }
  free(nthrds);
  free(nthrds1);
  free(bivinzcol);
  free(bivinzrow);
  ierr = PetscFree(submatA);
  ierr = PetscFree(submatB);
  CHKERRQ(ierr);
  free(yi1);
  printf("Completed calculation of partitioned matrices! Rank %d\n",rank);
  timestr=clock();

  #pragma omp parallel private(li) reduction(+:lj)
  {
  lj=0;
  #pragma omp for
  for(li=0; li<vecbivisize; ++li){
    if(vecbivi[li]!=0)lj++;
  }
  }
  long int *obiviindx1,*biviindx1,*obiviindx0,*biviindx0;
  obiviindx1=realloc(obiviindx1,(lj+1)*sizeof(long int));
  obiviindx1[0]=-1;
  biviindx1=obiviindx1;
  biviindx1++;
  long int nz0,nz1,nz2,nz3=lj,halfj2;
  long int lj2=0;
  for(li=0; li<vecbivisize; ++li) {
    if(vecbivi[li]!=0) {
      vecbivi[lj2]=vecbivi[li];
      biviindx1[lj2]=bivinzcol1[li%bivicolsize]+bivinzrow1[li/bivicolsize];//[i-j1*bivicolsize]+bivinzrow1[j1];//*nrowb;
      ++lj2;
    }
  }
  free(bivinzrow1);
  free(bivinzcol1);
  vecbivi=realloc(vecbivi,lj2*sizeof(ha_cgetype));
  nz0=lj2;
  printf("jjjj1 %ld j2 %ld rank %d nz0 %ld nz3 %ld\n",lj,lj2,rank,nz0,nz3);
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len,name_len_max,name_beg,class_size,color,group_size,ha_id,group_size1,group_size11,ha_id1;
  MPI_Comm_rank( HA_COMM, &ha_id);
  MPI_Comm_size(HA_COMM,&group_size);
  if(ha_id==group_size-1)color=1;
  else color=0;
  MPI_Comm_rank( HA1_COMM, &ha_id1);
  MPI_Comm_size(HA1_COMM,&group_size1);
  group_size11=group_size1;
  MPI_Bcast(&group_size11,1, MPI_INT,mpisize-1, PETSC_COMM_WORLD);
  int sindx01,sindx02,sindx03;

  if(ha_id!=group_size-1) {
    lj2=nz3;
    MPI_Send(&lj2,1, MPI_LONG,group_size-1, 10, HA_COMM);
    sindx01=lj2/MAXSSIZE;
    sindx03=0;
    for(sindx02=0; sindx02<sindx01; sindx02++) {
      if(SORD==1)MPI_Send(vecbivi+sindx03,MAXSSIZE, MPI_DOUBLE,group_size-1, 11*sindx01+sindx02, HA_COMM);
      else MPI_Send(vecbivi+sindx03,MAXSSIZE, MPI_FLOAT,group_size-1, 11*sindx01+sindx02, HA_COMM);
      sindx03+=MAXSSIZE;
    }
    if(SORD==1)MPI_Send(vecbivi+sindx03,lj2-sindx01*MAXSSIZE, MPI_DOUBLE,group_size-1, 11, HA_COMM);
    else MPI_Send(vecbivi+sindx03,lj2-sindx01*MAXSSIZE, MPI_FLOAT,group_size-1, 11, HA_COMM);
    free(vecbivi);
    vecbivi=NULL;
    sindx03=0;
    for(sindx02=0; sindx02<sindx01; sindx02++) {
      MPI_Send(biviindx1+sindx03,MAXSSIZE, MPI_LONG,group_size-1, 12*sindx01+sindx02, HA_COMM);
      sindx03+=MAXSSIZE;
    }
    MPI_Send(biviindx1+sindx03,lj2-sindx01*MAXSSIZE, MPI_LONG,group_size-1, 12, HA_COMM);
    free(obiviindx1);
    obiviindx1=NULL;
  }
  if(ha_id==group_size-1) {
    for(j6=0; j6<group_size-1; j6++) {
      MPI_Recv(&lj2,1, MPI_LONG,  MPI_ANY_SOURCE,10, HA_COMM,&status);
      i=status.MPI_SOURCE;
      obiviindx0=realloc(obiviindx0,(lj2+1)*sizeof(long int));
      obiviindx0[0]=-1;
      biviindx0=obiviindx0;
      biviindx0++;
      ha_cgetype *vecbivi0= (ha_cgetype *) calloc (lj2,sizeof(ha_cgetype));
      sindx01=lj2/MAXSSIZE;
      sindx03=0;
      for(sindx02=0; sindx02<sindx01; sindx02++) {
        if(SORD==1)MPI_Recv(vecbivi0+sindx03,MAXSSIZE, MPI_DOUBLE,i,11*sindx01+sindx02, HA_COMM,&status);
        else MPI_Recv(vecbivi0+sindx03,MAXSSIZE, MPI_FLOAT,i,11*sindx01+sindx02, HA_COMM,&status);
        sindx03+=MAXSSIZE;
      }
      if(SORD==1)MPI_Recv(vecbivi0+sindx03,lj2-sindx01*MAXSSIZE, MPI_DOUBLE,i,11, HA_COMM,&status);
      else MPI_Recv(vecbivi0+sindx03,lj2-sindx01*MAXSSIZE, MPI_FLOAT,i,11, HA_COMM,&status);
      sindx03=0;
      for(sindx02=0; sindx02<sindx01; sindx02++) {
        MPI_Recv(biviindx0+sindx03,MAXSSIZE, MPI_LONG,i,12*sindx01+sindx02, HA_COMM,&status);
        sindx03+=MAXSSIZE;
      }
      MPI_Recv(biviindx0+sindx03,lj2-sindx01*MAXSSIZE, MPI_LONG,i,12, HA_COMM,&status);
      my_spar_compl_(biviindx1,&nz0,biviindx0,&lj2,&nz1);
      vecbivi=realloc(vecbivi,nz1*sizeof(ha_cgetype));
      obiviindx1=realloc(obiviindx1,(nz1+1)*sizeof(long int));
      obiviindx1[0]=-1;
      biviindx1=obiviindx1;
      biviindx1++;
      my_spar_add3l_(vecbivi,obiviindx1,&nz0,vecbivi0,obiviindx0,&lj2,&nz1);
      free(obiviindx0);
      obiviindx0=NULL;
      free(vecbivi0);
      vecbivi0=NULL;
      printf("!!!!!!!!!!!! rank11 %d nz0 %ld j2 %ld nz1 %ld\n",rank,nz0,lj2,nz1);
      nz0=nz1;
    }
  }

  if(group_size11>1&&color==1) {
    nz3=nz0;
    if(ha_id1!=group_size1-1) {
      lj2=nz3;
      MPI_Send(&lj2,1, MPI_LONG,group_size1-1, 10, HA1_COMM);
      sindx01=lj2/MAXSSIZE;
      sindx03=0;
      for(sindx02=0; sindx02<sindx01; sindx02++) {
        if(SORD==1) MPI_Send(vecbivi+sindx03,MAXSSIZE, MPI_DOUBLE,group_size1-1, 11*sindx01+sindx02, HA1_COMM);
        else MPI_Send(vecbivi+sindx03,MAXSSIZE, MPI_FLOAT,group_size1-1, 11*sindx01+sindx02, HA1_COMM);
        sindx03+=MAXSSIZE;
      }
      if(SORD==1) MPI_Send(vecbivi+sindx03,lj2-sindx01*MAXSSIZE, MPI_DOUBLE,group_size1-1, 11, HA1_COMM);
      else MPI_Send(vecbivi+sindx03,lj2-sindx01*MAXSSIZE, MPI_FLOAT,group_size1-1, 11, HA1_COMM);
      sindx03=0;
      for(sindx02=0; sindx02<sindx01; sindx02++) {
        MPI_Send(biviindx1+sindx03,MAXSSIZE, MPI_LONG,group_size1-1, 12*sindx01+sindx02, HA1_COMM);
        sindx03+=MAXSSIZE;
      }
      MPI_Send(biviindx1+sindx03,lj2-sindx01*MAXSSIZE, MPI_LONG,group_size1-1, 12, HA1_COMM);
      free(obiviindx1);
      obiviindx1=NULL;
      free(vecbivi);
      vecbivi=NULL;
    }
    if(ha_id1==group_size1-1) {
      for(j6=0; j6<group_size1-1; j6++) {
        MPI_Recv(&lj2,1, MPI_LONG,  MPI_ANY_SOURCE,10, HA1_COMM,&status);
        i=status.MPI_SOURCE;
        obiviindx0=realloc(obiviindx0,(lj2+1)*sizeof(long int));
        obiviindx0[0]=-1;
        biviindx0=obiviindx0;
        biviindx0++;
        ha_cgetype *vecbivi0= (ha_cgetype *) calloc (lj2,sizeof(ha_cgetype));
        sindx01=lj2/MAXSSIZE;
        sindx03=0;
        for(sindx02=0; sindx02<sindx01; sindx02++) {
          if(SORD==1) MPI_Recv(vecbivi0+sindx03,MAXSSIZE, MPI_DOUBLE,i,11*sindx01+sindx02, HA1_COMM,&status);
          else MPI_Recv(vecbivi0+sindx03,MAXSSIZE, MPI_FLOAT,i,11*sindx01+sindx02, HA1_COMM,&status);
          sindx03+=MAXSSIZE;
        }
        if(SORD==1) MPI_Recv(vecbivi0+sindx03,lj2-sindx01*MAXSSIZE, MPI_DOUBLE,i,11, HA1_COMM,&status);
        else MPI_Recv(vecbivi0+sindx03,lj2-sindx01*MAXSSIZE, MPI_FLOAT,i,11, HA1_COMM,&status);
        sindx03=0;
        for(sindx02=0; sindx02<sindx01; sindx02++) {
          MPI_Recv(biviindx0+sindx03,MAXSSIZE, MPI_LONG,i,12*sindx01+sindx02, HA1_COMM,&status);
          sindx03+=MAXSSIZE;
        }
        MPI_Recv(biviindx0+sindx03,lj2-sindx01*MAXSSIZE, MPI_LONG,i,12, HA1_COMM,&status);
        my_spar_compl_(biviindx1,&nz0,biviindx0,&lj2,&nz1);
        vecbivi=realloc(vecbivi,nz1*sizeof(ha_cgetype));
        obiviindx1=realloc(obiviindx1,(nz1+1)*sizeof(long int));
        obiviindx1[0]=-1;
        biviindx1=obiviindx1;
        biviindx1++;
        my_spar_add3l_(vecbivi,obiviindx1,&nz0,vecbivi0,obiviindx0,&lj2,&nz1);
        free(obiviindx0);
        obiviindx0=NULL;
        free(vecbivi0);
        printf("!!!!!!!!!!!! rank22 %d nz0 %ld j2 %ld nz1 %ld\n",rank,nz0,lj2,nz1);
        nz0=nz1;
      }
    }
  }

  forint vecbiuisize=0;
  vecbiuisize=VecSize-sumrowcolin;
  printf("after bivi rank %d\n",rank);
  HaReduceNoComp(vecbiui,vecbiuisize,mpisize,rank,mpisize-1);
  printf("Completed MPI_Reduce Operation! Rank %d time %f\n",rank,((double)clock()-timestr)/CLOCKS_PER_SEC);
  if(rank!=mpisize-1){
    free(vecbiui);
    vecbiui=NULL;
  }
  timestr=clock();
  ha_cgetype *xd;//= (ha_cgetype*)calloc(vecbiuisize,sizeof(ha_cgetype));
  ha_cgetype *x0;//,*vecbivi0;//= (ha_cgetype*)calloc(VecSize,sizeof(ha_cgetype));
  long int j7,lnz;
  if(rank==mpisize-1) {
    Mat_SeqAIJ         *aa=(Mat_SeqAIJ*)submatD[0]->data;//*aa=subA->data;
    ai= aa->i;
    aj= aa->j;
    vals=aa->a;
    lnz=aa->nz;
    nrow=submatD[0]->rmap->n;
    ncol=submatD[0]->cmap->n;
    printf("D ncol %d nrow %d nz %ld\n",ncol,nrow,lnz);
    obiviindx0=realloc(obiviindx0,(lnz+1)*sizeof(long int));
    obiviindx0[0]=-1;
    biviindx0=obiviindx0;
    biviindx0++;
    lj2=0;
    for(li=0; li<nrow-1; li++) {
      j7=li*ncol;
      for(j=ai[li]; j<ai[li+1]; j++) {
        biviindx0[lj2]=aj[j]+j7;
        lj2++;
      }
    }
    j7=li*ncol;
    for(j=ai[li]; j<lnz; j++) {
      biviindx0[lj2]=aj[j]+j7;
      lj2++;
    }
    my_spar_compl_(biviindx1,&nz0,biviindx0,&lnz,&nz1);
    printf("nz0 %ld nz %ld nz1 %ld\n",nz0,lnz,nz1);



    printf("nz1 %ld\n",nz1);
    ldsize=ceil((laD/100.0)*nz1);
    vecbivi=realloc(vecbivi,ldsize*sizeof(ha_cgetype));
    int *irn1=(int *) calloc (nz1,sizeof(int));
    int *jcn=(int *) calloc (ldsize,sizeof(int));
    printf("indx0 %ld indx1 %ld\n",obiviindx1[0],obiviindx0[0]);
    my_spar_add4l_(vecbivi,obiviindx1,irn1,jcn,&nz0,vals,obiviindx0,&lnz,&nz1,&ncol);//my_spar_add1l_



    MatDestroy(&submatD[0]);
    lnz=nz1;
    free(obiviindx1);
    obiviindx1=NULL;
    free(obiviindx0);
    obiviindx0=NULL;
    irn1=realloc(irn1,ldsize*sizeof(int));
    int *insizeD=(int *) calloc (6,sizeof(int));
    insizeD[0]=nrow;
    insizeD[1]=ncol;
    insizeD[2]=lnz;
    insizeD[3]=laD;
    insizeD[5]=ldsize;
    printf("ok\n");
    xd=(ha_cgetype *) calloc (vecbiuisize,sizeof(ha_cgetype));//realloc (xd,vecbiuisize*sizeof(ha_cgetype));
    printf("ok!!!!!!!! vecbiuisize %ld nrow %d\n",vecbiuisize,nrow);
    spec48_ssol2la_(insizeD,irn1,jcn,vecbivi,vecbiui,xd);
    printf("ok!!!!!!!!\n");



    free(vecbivi);
    vecbivi=NULL;
    free(vecbiui);
    vecbiui=NULL;
    free(irn1);
    free(jcn);
    free(insizeD);
    x0=(ha_cgetype *) calloc (VecSize,sizeof(ha_cgetype));//realloc (x0,*sizeof(ha_cgetype));
    for(i=0; i<vecbiuisize; i++) {
      x0[indicesC[i]]=xd[i];
    }
  }
  free(vecbiui);
    obiviindx1=realloc(obiviindx1,sizeof(long int));
    free(obiviindx1);
    obiviindx1=NULL;
    obiviindx0=realloc(obiviindx0,sizeof(long int));
    free(obiviindx0);
    obiviindx0=NULL;
  free(vecbivi);
  free(indicesC);
  if(rank!=mpisize-1) {
    xd=(ha_cgetype *) calloc (vecbiuisize,sizeof(ha_cgetype));//realloc (xd,vecbiuisize*sizeof(ha_cgetype));
    x0=(ha_cgetype *) calloc (VecSize,sizeof(ha_cgetype));//realloc (x0,*sizeof(ha_cgetype));
  }
  PetscFree(submatD);
  MPI_Barrier(PETSC_COMM_WORLD);
  if(SORD==1)MPI_Bcast(xd, vecbiuisize, MPI_DOUBLE,mpisize-1, PETSC_COMM_WORLD);
  else MPI_Bcast(xd, vecbiuisize, MPI_FLOAT,mpisize-1, PETSC_COMM_WORLD);
  printf("Interface problem solution Rank %d time %f\n",rank,((double)clock()-timestr)/CLOCKS_PER_SEC);
  timestr=clock();
  xi1indx=0;
  for(j1=0; j1<nmatinplus; j1++) {
    if(j1<nmatin) {
      ha_cgetype *biui0= (ha_cgetype *) calloc (ha_ndblocks[j1+begblock[rank]],sizeof(ha_cgetype));

      fp1 = fopen(fn01[j1], "rb");
      if (fp1==NULL)printf("File opening error\n");
      fp2 = fopen(fn02[j1], "rb");
      if (fp2==NULL)printf("File opening error\n");
      fp3 = fopen(fn03[j1], "rb");
      if (fp3==NULL)printf("File opening error\n");
      la1=ceil((insize[j1*insizes+9]/100.0)*insize[j1*insizes+2]);
      int *irne = (int*)calloc(la1,sizeof(int));
      int *keep = (int*)calloc(insize[j1*insizes+12],sizeof(int));
      ha_cgetype *vale = (ha_cgetype*)calloc(la1,sizeof(ha_cgetype));
      freadresult=fread(irne,sizeof(int),la1,fp1);
      freadresult=fread(keep,sizeof(int),insize[j1*insizes+12],fp2);
      freadresult=fread(vale,sizeof(ha_cgetype),la1,fp3);
      fclose(fp1);
      fclose(fp2);
      fclose(fp3);
      remove(fn01[j1]);
      remove(fn02[j1]);
      remove(fn03[j1]);
      Mat_SeqAIJ         *ac=(Mat_SeqAIJ*)submatC[j1]->data;//*aa=subA->data;
      ai= ac->i;
      aj= ac->j;
      vals=ac->a;
      nz=ac->nz;
      nrow=submatC[j1]->rmap->n;
      ncol=submatC[j1]->cmap->n;
      ha_cgetype *be0 = (ha_cgetype*)calloc(nrow,sizeof(ha_cgetype));
      for(i=0; i<nrow-1; i++) {
        for(j=ai[i]; j<ai[i+1]; j++) {
          be0[i]+=vals[j]*xd[aj[j]];
        }
      }
      for(j=ai[i]; j<nz; j++) {
        be0[i]+=vals[j]*xd[aj[j]];
      }
      MatDestroy(&submatC[j1]);
      if(insize[j1*insizes+16]!=la1)insize[j1*insizes+16]=la1;
      spec48m_esol_(insize+j1*insizes,irne,vale,keep,be0,biui0);
      free(irne);
      free(keep);
      free(vale);
      free(be0);
      ISGetIndices(colindices[j1],&nindices);
      xi1point=xi1+xi1indx;
      for(i=0; i<ha_ndblocks[j1+begblock[rank]]; i++) {
        x0[nindices[i]]+=xi1point[i]-biui0[i];//vals[i]-biui0[i];//vecval;
      }
      ISRestoreIndices(colindices[j1],&nindices);
      free(biui0);
      xi1indx+=ha_ndblocks[j1+begblock[rank]];
    }
  }
  free(xi1);
  free(insize);
  ierr = PetscFree(submatC);
  free(begblock);
  printf("Solution calculation Rank %d time %f\n",rank,((double)clock()-timestr)/CLOCKS_PER_SEC);
  MPI_Barrier(PETSC_COMM_WORLD);
  free(xd);
  timestr=clock();
  if(SORD==1)MPI_Allreduce(x0,x1,VecSize, MPI_DOUBLE, MPI_SUM,PETSC_COMM_WORLD);
  else MPI_Allreduce(x0,x1,VecSize, MPI_FLOAT, MPI_SUM,PETSC_COMM_WORLD);
  printf("Reduce solution Rank %d time %f\n",rank,((double)clock()-timestr)/CLOCKS_PER_SEC);
  free(x0);
  for (i=0; i<nmatinplus; i++){
    free(fn01[i]);
    free(fn02[i]);
    free(fn03[i]);
  }
  free(fn01);
  free(fn02);
  free(fn03);
  for (i=0; i<nmatin; i++) {
    ierr = ISDestroy(&colindices[i]);
    CHKERRQ(ierr);
  }
  PetscFree(colindices);
  printf("OK rui rank %d\n",rank);
  return 0;
}

int HaNDBBDParPre(Mat A, Vec b, ha_cgetype *x1, uvadd VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int *ha_rows,int *ha_cols, uvadd ndblock,uvadd nreg,uvadd ntime,int *ha_ndblocks, uvadd *countvarintra1, uvadd *counteq, uvadd *counteqnoadd,uvdim laA,uvdim laDi,uvdim laD,PetscReal cntl3,PetscReal cntl6,PetscBool presol) {//,bool iter
  IS *rowindices=NULL,*colindices=NULL,*colindicesbc1=NULL,*rowBBij=NULL,*colBBij=NULL;//,*colindicesbcpm,*colindicesbcpm1
  PetscInt bfirst,bend,sblockin,nmatin,nmatinplus,nmatint,nmatinplust,nmatminust,nrowcolin,i,i1,j,j0,j1,j2,j3,j4,j5,j6,j7,l0,l1,l2,l3,l4,l5,rank1,proc1=0,nnzmax,j1nz,j1irnbs;//,sumrowcolin
  Mat *submatAij=NULL,*submatBBij=NULL;//,*submatCij,*submatBij;,*submatB
  Vec vecxd;
  PetscInt *ai,*aj,*a1i,*a1j,*a2i,*a2j,*aic,*ajc;
  PetscReal cntl6in;
  uvadd lasize,ldsize;
  size_t freadresult,fwrt;
  PetscInt nrow,ncol,nz,nrowc,ncolc,nrowb,ncolb,nzc,nzv;
  MPI_Status   status;
  int j2int,la1,fd1,fd2,fd3,unequal=0;
  bool ifremove=false;
  PetscScalar *vals,*val1s,*val2s,*valsc,vecval;
  ha_cgetype vval;
  PetscErrorCode ierr;
  PetscViewer viewer;
  FILE* presolfile;
  int *begblock= (int *) calloc (mpisize,sizeof(int));
  int *ndblockinrank= (int *) calloc (mpisize,sizeof(int));
  char filename[1024],rankname[1024],j1name[1024],tempchar[1024];//,fn1[1024],fn2[1024],fn3[1024]
  if(rank<10)strcpy(rankname,"000");
  if(rank<100&&rank>9)strcpy(rankname,"00");
  if(rank<1000&&rank>99)strcpy(rankname,"0");
  if(rank>=1000)rankname[0]='\0';
  sprintf(filename, "%d",rank);
  strcat(rankname,filename);
  int nmatinBB=1;
  nmatint=(uvadd)ntime/mpisize;
  if(nmatint*mpisize<ntime)unequal=1;
  nmatminust=nmatint;
  if(nmatint*mpisize<ntime)nmatinplust=nmatint+1;
  else nmatinplust=nmatint;
  for(i=0; i<mpisize; i++)if(rank+1<=ntime-mpisize*nmatint)nmatint++;
  nmatin=(nreg+1)*nmatint;
  nmatinplus=(nreg+1)*nmatinplust;
  begblock[rank]=nmatin;
  printf("rank %d nmatin %d nmatint %d nmplus %d\n",rank,nmatin,nmatint,nmatinplus);
  for(i=0; i<mpisize; i++) {
    j=rank;
    j2=j;
    MPI_Bcast(&j2,1, MPI_INT,i, PETSC_COMM_WORLD);
    j=j2;
    MPI_Bcast(&begblock[j2],1, MPI_INT,i, PETSC_COMM_WORLD);
  }
  j2=begblock[0];
  begblock[0]=0;
  for(i=1; i<mpisize; i++) {
    j1=begblock[i];
    begblock[i]=begblock[i-1]+j2;
    j2=j1;
  }
  j2=nmatin;
  MPI_Bcast(&j2,1, MPI_INT,0, PETSC_COMM_WORLD);
  ndblockinrank[0]=j2;
  for(i=1; i<mpisize; i++) {
    j2=nmatin;
    MPI_Bcast(&j2,1, MPI_INT,i, PETSC_COMM_WORLD);
    ndblockinrank[i]=ndblockinrank[i-1]+j2;
  }
  for(i=0; i<ndblock; i++) {
    j2=i-ndblockinrank[rank];
    if(j2<0&&j2>=-nmatin) {
      j1=i/(nreg+1);
      j=i-(nreg+1)*j1;
      if(rank==0)j3=0;
      else j3=ndblockinrank[rank-1];
    }
  }
  j1=0;
  for(j=0; j<mpisize; j++)if(j<rank)j1+=ndblockinrank[j];
  free(ndblockinrank);
  ierr = PetscMalloc(nmatin*sizeof(IS **),&rowindices);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatin*sizeof(IS **),&colindices);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatint*sizeof(IS **),&colindicesbc1);
  CHKERRQ(ierr);
  Mat** submatBij = (Mat**)calloc(nmatin,sizeof(Mat*));
  Mat** submatCij = (Mat**)calloc(nmatin,sizeof(Mat*));


  char **fn01= (char**)calloc(nmatinplus,sizeof(char*));
  for (i=0; i<nmatinplus; i++) fn01[i] = (char*)calloc(1024,sizeof(char));
  char **fn02= (char**)calloc(nmatinplus,sizeof(char*));
  for (i=0; i<nmatinplus; i++) fn02[i] = (char*)calloc(1024,sizeof(char));
  char **fn03= (char**)calloc(nmatinplus,sizeof(char*));
  for (i=0; i<nmatinplus; i++) fn03[i] = (char*)calloc(1024,sizeof(char));

  for(j1=0; j1<nmatinplus; j1++) {
    if(j1<nmatin) {
      if(j1<10)strcpy(j1name,"000");
      if(j1<100&&j1>9)strcpy(j1name,"00");
      if(j1<1000&&j1>99)strcpy(j1name,"0");
      if(j1>=1000)j1name[0]='\0';
      sprintf(filename, "%d",j1);
      strcat(j1name,filename);
      strcpy(filename,"_irnv");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      strcpy(fn01[j1],filename);
      strcpy(filename,"_keep");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      strcpy(fn02[j1],filename);
      strcpy(filename,"_vav");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      strcpy(fn03[j1],filename);
    }
  }
  for(j1=0; j1<nmatinplus; j1++) {
    free(fn01[j1]);//1
    free(fn02[j1]);//1
    free(fn03[j1]);//1
  }
  free(fn01);//1
  free(fn02);//1
  free(fn03);//1
  j1=0;
  for(i=0; i<nmatin; i++) {
    bfirst=counteq[i+begblock[rank]];
    PetscInt *indices= (PetscInt *) calloc (ha_ndblocks[i+begblock[rank]],sizeof(PetscInt));
    for(j=0; j<ha_ndblocks[i+begblock[rank]]; j++) {
      indices[j]=ha_rows[bfirst+j];//+bfirst;
    }

    ISCreateGeneral(PETSC_COMM_SELF,ha_ndblocks[i+begblock[rank]],indices,PETSC_COPY_VALUES,rowindices+i);
    bfirst=countvarintra1[i+begblock[rank]];
    for(j=0; j<ha_ndblocks[i+begblock[rank]]; j++) {
      indices[j]=ha_cols[bfirst+j];//+bfirst;
    }
    ISCreateGeneral(PETSC_COMM_SELF,ha_ndblocks[i+begblock[rank]],indices,PETSC_COPY_VALUES,colindices+i);
    free(indices);//1
  }
  for(i=0; i<nmatint; i++) {
    j2=0;
    for(j1=0; j1<nreg+1; j1++)j2+=ha_ndblocks[i*(nreg+1)+j1+begblock[rank]];
    PetscInt *indices= (PetscInt *) calloc (j2,sizeof(PetscInt));
    j2=0;
    for(j1=0; j1<nreg; j1++) {
      bfirst=countvarintra1[i*(nreg+1)+j1+begblock[rank]];
      for(j=0; j<ha_ndblocks[i*(nreg+1)+j1+begblock[rank]]; j++) {
        indices[j2]=ha_cols[bfirst+j];//+bfirst;
        j2++;
      }
    }
    ISCreateGeneral(PETSC_COMM_SELF,j2,indices,PETSC_COPY_VALUES,colindicesbc1+i);
    free(indices);
  }
  printf("Begin preparation rank %d\n",rank);
  ierr = MatCreateSubMatrices(A,nmatin,rowindices,colindices,MAT_INITIAL_MATRIX,&submatAij);
  int insizes=17;
  int *insize=(int *) calloc (insizes*nmatin,sizeof(int));
  printf("rank %d A\n",rank);
  PetscInt *indicesbbij=NULL;//(int *) calloc (1,sizeof(PetscInt));
  ierr = PetscMalloc(1*sizeof(IS **),&rowBBij);
  CHKERRQ(ierr);
  ierr = PetscMalloc(1*sizeof(IS **),&colBBij);
  CHKERRQ(ierr);

  for(i=0; i<nmatminust; i++) {
    j2=i*(nreg+1)+nreg;
    ierr = MatCreateSubMatrices(A,1,rowindices+j2,colindicesbc1+i,MAT_INITIAL_MATRIX,&submatBBij);
    nrow=submatBBij[0]->rmap->n;
    ncol=submatBBij[0]->cmap->n;
    if(nrow>ncol) {
      indicesbbij=realloc(indicesbbij,nrow*sizeof(PetscInt));
      for(j=0; j<nrow; j++)indicesbbij[j]=j;
    }
    else {
      indicesbbij=realloc(indicesbbij,ncol*sizeof(PetscInt));
      for(j=0; j<ncol; j++)indicesbbij[j]=j;
    }
    ISCreateGeneral(PETSC_COMM_SELF,nrow,indicesbbij,PETSC_COPY_VALUES,rowBBij);
    j3=0;
    for(j=0; j<nreg; j++) {
      j1=i*(nreg+1)+j;
      ierr = MatCreateSubMatrices(A,1,rowindices+j1,colindices+j2,MAT_INITIAL_MATRIX,&submatCij[j1]);
      ISCreateGeneral(PETSC_COMM_SELF,ha_ndblocks[j1+begblock[rank]],indicesbbij+j3,PETSC_COPY_VALUES,colBBij);
      ierr = MatCreateSubMatrices(submatBBij[0],1,rowBBij,colBBij,MAT_INITIAL_MATRIX,&submatBij[j1]);
      ierr = ISDestroy(&colBBij[0]);
      j3+=ha_ndblocks[j1+begblock[rank]];
    }
    ierr = MatDestroy(&submatBBij[0]);
    ierr = PetscFree(submatBBij);//1
    submatBBij=NULL;
    ierr = ISDestroy(&rowBBij[0]);
  }

  if(unequal==1){
  if(nmatint<nmatinplust) for(i=0; i<nreg; i++) {//to create again! MatCreateSubMatrices must be the same across
      MatDestroy(&submatCij[(nmatint-1)*(nreg+1)+i][0]);
      ierr = PetscFree(submatCij[(nmatint-1)*(nreg+1)+i]);
      MatDestroy(&submatBij[(nmatint-1)*(nreg+1)+i][0]);
      ierr = PetscFree(submatBij[(nmatint-1)*(nreg+1)+i]);
    }
  j2=(nmatint-1)*(nreg+1)+nreg;
  ierr = MatCreateSubMatrices(A,1,rowindices+j2,colindicesbc1+(nmatint-1),MAT_INITIAL_MATRIX,&submatBBij);
  nrow=submatBBij[0]->rmap->n;
  ncol=submatBBij[0]->cmap->n;
  if(nrow>ncol) {
    indicesbbij=realloc(indicesbbij,nrow*sizeof(PetscInt));
    for(j=0; j<nrow; j++)indicesbbij[j]=j;
  }
  else {
    indicesbbij=realloc(indicesbbij,ncol*sizeof(PetscInt));
    for(j=0; j<ncol; j++)indicesbbij[j]=j;
  }
  ISCreateGeneral(PETSC_COMM_SELF,nrow,indicesbbij,PETSC_COPY_VALUES,rowBBij);
  j3=0;
  for(j=0; j<nreg; j++) {
    j1=(nmatint-1)*(nreg+1)+j;
    ierr = MatCreateSubMatrices(A,1,rowindices+j1,colindices+j2,MAT_INITIAL_MATRIX,&submatCij[j1]);
    ISCreateGeneral(PETSC_COMM_SELF,ha_ndblocks[j1+begblock[rank]],indicesbbij+j3,PETSC_COPY_VALUES,colBBij);
    ierr = MatCreateSubMatrices(submatBBij[0],1,rowBBij,colBBij,MAT_INITIAL_MATRIX,&submatBij[j1]);
    ierr = ISDestroy(&colBBij[0]);
    j3+=ha_ndblocks[j1+begblock[rank]];
  }
  ierr = MatDestroy(&submatBBij[0]);
  ierr = PetscFree(submatBBij);//1
  submatBBij=NULL;
  ierr = ISDestroy(&rowBBij[0]);
  }
  free(begblock);
  free(indicesbbij);
  PetscFree(rowBBij);
  PetscFree(colBBij);

  for (i=0; i<nmatin; i++) {
    ierr = ISDestroy(&rowindices[i]);//1
    CHKERRQ(ierr);
    ierr = ISDestroy(&colindices[i]);//1
    CHKERRQ(ierr);
  }
  for (i=0; i<nmatint; i++) {
    ierr = ISDestroy(&colindicesbc1[i]);//1
    CHKERRQ(ierr);
  }
  PetscFree(colindices);//1
  PetscFree(rowindices);//1
  PetscFree(colindicesbc1);//1


  int bivirowsize=1,bivicolsize=1,bbrowij,ddrowi;
  long int vecbivisize=0,li,lj;//,halfvec;
  long int nz0,nz1,nz2,nz3=j,halfj2;
  long int lj2;
  int jthrd,nthrd;
  #pragma omp parallel private(jthrd,nthrd,j3,j4,bivirowsize,bivicolsize,bbrowij,ai,nz,nrow,i,j,j1,j2,li,lj,ddrowi,vecbivisize,aic,ajc,valsc,nzc,nrowc,ncolc,ncolb,nrowb,aj,vals,lj2,nz0,j1name,filename,presolfile,fwrt,fd1,nz1,cntl6in,ncol,lasize,ldsize) shared(insize,submatAij,submatBij,submatCij)
  {
  long int *bivinzrow=NULL,irnmems=0;//(long int *) calloc (1,sizeof(long int));
  PetscInt *bivinzcol=NULL;//(PetscInt *) calloc (1,sizeof(PetscInt));
  PetscReal *ccolnorms=NULL;//(PetscReal *) calloc (1,sizeof(PetscReal));
  long int *bivinzrow1=NULL;//(long int *) calloc (1,sizeof(long int));
  long int *bivinzcol1=NULL;//(long int *) calloc (1,sizeof(long int));
  int *irn=NULL,*jcn=NULL,*irn1=NULL,*jcn1=NULL;
  int windx=0,bindx,eindx;
  jthrd=omp_get_thread_num();
  nthrd=omp_get_max_threads( );
  nthrd=90/nthrd;
  bindx=nthrd*jthrd;
  eindx=bindx+nthrd;
  windx=bindx;
  Mat_SeqAIJ         *aa,*ac,*ab;
  Mat submatCT;
  ha_cgetype *vecbivi=NULL;
  long int bbrowijmax=0,bivirowsizemax=0,bivicolsizemax=0,vecbivisizemax=0,lasizemax=0,nz1max=0,nrowmax=0,ncolmax=0,nrowamax=0,ncolamax=0,nebmax=0;
  #pragma omp for schedule (static)
  for(j3=0; j3<nmatint; j3++) {
    bbrowij=submatBij[j3*(nreg+1)][0]->rmap->n;
    if(bbrowij>bbrowijmax)bbrowijmax=bbrowij;
  }
  ccolnorms=realloc(ccolnorms,bbrowijmax*sizeof(PetscReal));
  bivinzrow=realloc(bivinzrow,bbrowijmax*sizeof(long int));
  bivinzcol=realloc(bivinzcol,bbrowijmax*sizeof(PetscInt));
  #pragma omp for schedule (static)
  for(j3=0; j3<nmatint; j3++) {
    bivirowsize=1;
    bivicolsize=1;
    bbrowij=submatBij[j3*(nreg+1)][0]->rmap->n;
    memset(bivinzrow,0,bbrowij*sizeof(long int));
    memset(bivinzcol,0,bbrowij*sizeof(PetscInt));
    for(j1=0; j1<nreg; j1++) {
      MatGetColumnNorms(submatCij[j1+j3*(nreg+1)][0],NORM_1,ccolnorms);
      for(i=0; i<bbrowij; i++)if(ccolnorms[i]>0)bivinzcol[i]++;
      ab=(Mat_SeqAIJ*)submatBij[j1+j3*(nreg+1)][0]->data;//*aa=subA->data;
      ai= ab->i;
      nz=ab->nz;
      if(nebmax<nz)nebmax=nz;
      nrow=submatBij[j1+j3*(nreg+1)][0]->rmap->n;
      for(i=0; i<nrow-1; i++)if(ai[i]!=ai[i+1])bivinzrow[i]++;
      if(ai[i]<nz)bivinzrow[i]++;
    }
    ddrowi=submatBij[j3*(nreg+1)][0]->rmap->n;
    for(i=0; i<ddrowi; i++) {
      if(bivinzrow[i]>0) {
        bivinzrow[i]=bivirowsize;
        bivirowsize++;
      }
      if(bivinzcol[i]>0) {
        bivinzcol[i]=bivicolsize;
        bivicolsize++;
      }
    }
    bivirowsize--;
    bivicolsize--;
    if(bivirowsize>bivirowsizemax)bivirowsizemax=bivirowsize;
    if(bivicolsize>bivicolsizemax)bivicolsizemax=bivicolsize;
    vecbivisize=((long int)bivirowsize)*((long int)bivicolsize);
    if(vecbivisize>vecbivisizemax)vecbivisizemax=vecbivisize;
    for(j1=0; j1<nreg; j1++) {
      j4=j1+j3*(nreg+1);
      aa=(Mat_SeqAIJ*)submatAij[j4]->data;//*aa=subA->data;
      nrow=submatAij[j4]->rmap->n;
      if(nrowamax<nrow)nrowamax=nrow;
      ncol=submatAij[j4]->cmap->n;
      if(ncolamax<ncol)ncolamax=ncol;
      nz=aa->nz;
      lasize=ceil((laA/100.0)*nz);
      if(lasize>lasizemax)lasizemax=lasize;
    }
  }
  lasizemax+=10;

      int *insized=(int *) calloc (5+nreg*insizes,sizeof(int));
    bivinzrow1=realloc (bivinzrow1,bivirowsizemax*sizeof(long int));
    if(bivinzrow1==NULL)printf("Error!!!! No memory left for bivinzrow1!\n");
    bivinzcol1=realloc (bivinzcol1,bivicolsizemax*sizeof(long int));
    if(bivinzcol1==NULL)printf("Error!!!! No memory left for bivinzcol1!\n");
      vecbivi= (ha_cgetype *) malloc (vecbivisizemax*sizeof(ha_cgetype));
    if(vecbivi==NULL)printf("Error!!!! No memory left for vecbivi!\n");
      irn=(int *) malloc (lasizemax*sizeof(int));
    if(irn==NULL)printf("Error!!!! No memory left for irn!\n");
      jcn=(int *) malloc (lasizemax*sizeof(int));
    if(jcn==NULL)printf("Error!!!! No memory left for jcn!\n");
      ha_cgetype *values= (ha_cgetype *) calloc (lasizemax,sizeof(ha_cgetype));
    if(values==NULL)printf("Error!!!! No memory left for values!\n");
      int* jcnb1=(int *) malloc (nebmax*sizeof(int));
    if(jcnb1==NULL)printf("Error!!!! No memory left for jcnb1!\n");
      ha_cgetype *sol48= (ha_cgetype *) malloc (nrowamax*sizeof(ha_cgetype));
    if(sol48==NULL)printf("Error!!!! No memory left for sol48!\n");
      ha_cgetype *b48= (ha_cgetype *) malloc (nrowamax*sizeof(ha_cgetype));
    if(b48==NULL)printf("Error!!!! No memory left for b48!\n");
      int* keep=NULL,*iw51=NULL;
  keep=(int *) malloc ((nrowamax+5*ncolamax+4*ncolamax+7)*sizeof(int));//ICNTL(6)=1
    if(keep==NULL)printf("Error!!!! No memory left for keep!\n");
  ha_cgetype *w51=NULL;
  if(nrowamax>ncolamax) w51=(ha_cgetype *) malloc ((5*nrowamax)*sizeof(ha_cgetype));
  else w51=(ha_cgetype *) malloc ((5*ncolamax)*sizeof(ha_cgetype));
    if(w51==NULL)printf("Error!!!! No memory left for w51!\n");
  iw51=(int *) malloc ((6*nrowamax+3*ncolamax)*sizeof(int));
    if(iw51==NULL)printf("Error!!!! No memory left for iw51!\n");
  #pragma omp for schedule (static)
  for(j3=0; j3<nmatint; j3++) {
    bivirowsize=1;
    bivicolsize=1;
    bbrowij=submatBij[j3*(nreg+1)][0]->rmap->n;
    printf("j3 %d bbrowij %djth %d windx %d\n",j3,bbrowij,jthrd,windx);
    memset(bivinzrow,0,bbrowij*sizeof(long int));
    memset(bivinzcol,0,bbrowij*sizeof(PetscInt));
    for(j1=0; j1<nreg; j1++) {
      MatGetColumnNorms(submatCij[j1+j3*(nreg+1)][0],NORM_1,ccolnorms);
      for(i=0; i<bbrowij; i++)if(ccolnorms[i]>0)bivinzcol[i]++;
      ab=(Mat_SeqAIJ*)submatBij[j1+j3*(nreg+1)][0]->data;//*aa=subA->data;
      ai= ab->i;
      nz=ab->nz;
      nrow=submatBij[j1+j3*(nreg+1)][0]->rmap->n;
      for(i=0; i<nrow-1; i++)if(ai[i]!=ai[i+1])bivinzrow[i]++;
      if(ai[i]<nz)bivinzrow[i]++;
    }
    ddrowi=submatBij[j3*(nreg+1)][0]->rmap->n;
    for(i=0; i<ddrowi; i++) {
      if(bivinzrow[i]>0) {
        bivinzrow[i]=bivirowsize;
        bivirowsize++;
      }
      if(bivinzcol[i]>0) {
        bivinzcol[i]=bivicolsize;
        bivicolsize++;
      }
    }
    bivirowsize--;
    bivicolsize--;
    j=0;
    j2=0;
    for(li=0; li<ddrowi; li++) {
      if(bivinzrow[li]>0) {
        bivinzrow1[j]=li*bbrowij;
        j++;
      }
      if(bivinzcol[li]>0) {
        bivinzcol1[j2]=li;
        j2++;
      }
      bivinzrow[li]=(bivinzrow[li]-1)*bivicolsize;//faster MSOL
    }
    vecbivisize=((long int)bivirowsize)*((long int)bivicolsize);
    printf("rank %d vecbivisize %ld bivirowsize %d bivicolsize %d\n",rank,vecbivisize,bivirowsize,bivicolsize);
    if(vecbivisize<0)printf("Overflow detected, the interface problem is too big!!!!!!\n");
    memset(vecbivi,0,vecbivisize*sizeof(ha_cgetype));
    for(j1=0; j1<nreg; j1++) {
      j4=j1+j3*(nreg+1);
      aa=(Mat_SeqAIJ*)submatAij[j4]->data;//*aa=subA->data;
      MatTranspose(submatCij[j4][0],MAT_INITIAL_MATRIX,&submatCT);
      MatDestroy(&submatCij[j4][0]);//1
      ac=(Mat_SeqAIJ*)submatCT->data;//*aa=subA->data;
      ai= aa->i;
      aj= aa->j;
      vals=aa->a;
      nz=aa->nz;
      nrow=submatAij[j4]->rmap->n;
      ncol=submatAij[j4]->cmap->n;
      lasize=ceil((laA/100.0)*nz);
      insize[j4*insizes+13]=bivirowsize;
      insize[j4*insizes+14]=bivicolsize;
      insize[j4*insizes]=nrow;
      insize[j4*insizes+1]=ncol;
      insize[j4*insizes+2]=nz;
      insize[j4*insizes+15]=7+windx;//7+rank+jthrd+windx;//(j4+jthrd)%90+7;
      windx++;
      if(windx==eindx)windx=bindx;
      aic= ac->i;
      ajc= ac->j;
      valsc=ac->a;
      nzc=ac->nz;
      nrowc=submatCT->rmap->n;
      ncolc=submatCT->cmap->n;
      insize[j4*insizes+3]=nrowc;
      insize[j4*insizes+4]=ncolc;
      memcpy (irn+nz,ai,(nrow+1)*sizeof(PetscInt));
      memcpy (jcn,aj,nz*sizeof(PetscInt));
      memcpy (values,vals,nz*sizeof(ha_cgetype));
      insize[j4*insizes+2]=nz;
      insize[j4*insizes+5]=nzc;
      MatDestroy(&submatAij[j4]);//1
      insize[j4*insizes+9]=laA;
      insize[j4*insizes+16]=lasize;
      insize[j4*insizes+10]=rank;
      insize[j4*insizes+11]=j4;
      MatGetSize(submatBij[j4][0],&nrowb,&ncolb);
      ab=(Mat_SeqAIJ*)submatBij[j4][0]->data;//*aa=subA->data;
      ai= ab->i;//irnbs+j1irnbs;
      aj= ab->j;//jcnbs+j1nz;
      vals=ab->a;//valbs+j1nz;
      nz=ab->nz;//nzbs[j1];
      insize[j4*insizes+6]=nrowb;//ncolb;
      insize[j4*insizes+7]=ncolb;//nrowb;
      insize[j4*insizes+8]=nz;
      prep48m_msol_(insize+j4*insizes,irn,jcn,values,aic,ajc,valsc,ai,aj,vals,vecbivi,bivinzrow,bivinzcol,jcnb1,sol48,b48,w51,iw51,keep);
      insize[j4*insizes+15]=0;
      MatDestroy(&submatCT);
      MatDestroy(&submatBij[j4][0]);//1
    }

    j4=nreg+j3*(nreg+1);
    aa=(Mat_SeqAIJ*)submatAij[j4]->data;//*aa=subA->data;
    ai= aa->i;
    aj= aa->j;
    vals=aa->a;
    nz=aa->nz;

    lj=0;
    for(li=vecbivisize-1; li!=-1; --li)if(vecbivi[li]!=0)lj++;
    irnmems=nz+lj;
    #pragma omp critical
    {
    irn1=realloc(irn1,irnmems*sizeof(int));
    if(irn1==NULL)printf("Error!!!! No memory left for irn1!\n");
    jcn1=realloc(jcn1,irnmems*sizeof(int));
    if(jcn1==NULL)printf("Error!!!! No memory left for jcn1!\n");

    lj2=0;
    for(i=0; i<bivirowsize; i++) {
      bivinzrow1[i]/=bbrowij;
      bivinzrow1[i]++;
    }
    for(i=0; i<bivicolsize; i++) {
      bivinzcol1[i]++;
    }
    for(li=0; li<vecbivisize; ++li) {
      if(vecbivi[li]!=0) {
        vecbivi[lj2]=vecbivi[li];
        irn1[lj2]=bivinzrow1[li/bivicolsize];
        jcn1[lj2]=bivinzcol1[li%bivicolsize];
        lj2++;
      }
    }
    nz0=lj2;
    nrow=submatAij[j4]->rmap->n;
    if(nrow>nrowmax)nrowmax=nrow;
    ncol=submatAij[j4]->cmap->n;
    if(ncol>ncolmax)ncolmax=ncol;
    printf("D ncol %d nrow %d nz %d\n",ncol,nrow,nz);
    for(i=1; i<nrow; i++) {
      for(j=ai[i-1]; j<ai[i]; j++) {
        vecbivi[lj2]=vals[j];
        irn1[lj2]=i;
        jcn1[lj2]=aj[j]+1;
        lj2++;
      }
    }
    for(j=ai[i-1]; j<nz; j++) {
      vecbivi[lj2]=vals[j];
      irn1[lj2]=i;
      jcn1[lj2]=aj[j]+1;
      lj2++;
    }
    nz1=lj2;
    if(nz1>nz1max)nz1max=nz1;
    MatDestroy(&submatAij[j4]);//1
    //begin save for sol
    if(j3<10)strcpy(j1name,"000");
    if(j3<100&&j3>9)strcpy(j1name,"00");
    if(j3<1000&&j3>99)strcpy(j1name,"0");
    if(j3>=1000)j1name[0]='\0';
    sprintf(filename, "%d",j3);
    strcat(j1name,filename);
    if(isLinux==0) {
      strcpy(filename,"_bivi");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "wb"))==NULL) {
        printf("Cannot open file.\n");
      }
      fwrt=fwrite(vecbivi, sizeof(ha_cgetype),nz1, presolfile);
      if(fwrt==0)printf("File write error.");
      fclose(presolfile);
      strcpy(filename,"_rbvi");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "wb"))==NULL) {
        printf("Cannot open file.\n");
      }
      fwrt=fwrite(irn1, sizeof(int), nz1, presolfile);
      if(fwrt== 0) printf("File write error.");
      fclose(presolfile);
      strcpy(filename,"_cbvi");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "wb"))==NULL) {
        printf("Cannot open file.\n");
      }
      fwrt=fwrite(jcn1, sizeof(int), nz1, presolfile);
      if(fwrt== 0) printf("File write error.");
      fclose(presolfile);
    }
    else {
      strcpy(filename,"_bivi");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      fd1 = open(filename,O_WRONLY|O_CREAT|O_TRUNC, 0644);
      fwrt = write(fd1,vecbivi, sizeof(ha_cgetype)*nz1);
      close(fd1);
      strcpy(filename,"_rbvi");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      fd1 = open(filename,O_WRONLY|O_CREAT|O_TRUNC, 0644);
      fwrt = write(fd1,irn1, sizeof(int)*nz1);
      close(fd1);
      strcpy(filename,"_cbvi");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      fd1 = open(filename,O_WRONLY|O_CREAT|O_TRUNC, 0644);
      fwrt = write(fd1,jcn1, sizeof(int)*nz1);
      close(fd1);
    }
    free(irn1);//1
    irn1=NULL;
    free(jcn1);//1
    jcn1=NULL;
    }
    //end save for sol

    nz=nz1;
    insize[j4*insizes]=nrow;
    insize[j4*insizes+1]=ncol;
    insize[j4*insizes+2]=nz;
    insize[j4*insizes+10]=rank;
    insize[j4*insizes+11]=j4;
    insize[j4*insizes+9]=laDi;
    insize[j4*insizes+15]=7+windx;//7+rank+jthrd+windx;//(j4+jthrd)%90+7;
      insized[0]=nrow;
      insized[1]=ncol;
      insized[2]=nz;
      insized[4]=laDi;
      for(i=0; i<nreg; i++) {
        for(j=0; j<insizes; j++) {
          insized[i*insizes+5+j]=insize[(i+j3*(nreg+1))*insizes+j];
        }
      }
      if(j3<10)strcpy(j1name,"000");
      if(j3<100&&j3>9)strcpy(j1name,"00");
      if(j3<1000&&j3>99)strcpy(j1name,"0");
      if(j3>=1000)j1name[0]='\0';
      sprintf(filename, "%d",j3);
      strcat(j1name,filename);
      strcpy(filename,"_rank");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "wb"))==NULL) {
        printf("Cannot open file.\n");
      }
      printf("Rank %d RRRRR! ncol %d nrow %d rank %d nz %d cntl6 %lf\n",rank,nrow,ncol,insized[3],nz,cntl6);
      fwrt=fwrite(insized, sizeof(int), 5+nreg*insizes, presolfile);
      if(fwrt==0)printf("File write error.");
      fclose(presolfile);
    insize[j4*insizes+15]=0;
  }
      free(irn);//1
      irn=NULL;
      free(jcn);//1
      jcn=NULL;
      free(values);//1
    free(jcnb1);
    jcnb1=NULL;
    free(sol48);
    sol48=NULL;
    free(b48);
    b48=NULL;
    free(keep);
    keep=NULL;
    free(w51);
    w51=NULL;
    free(iw51);
    iw51=NULL;
    free(insized);
    free(ccolnorms);//1
    ccolnorms=NULL;
    free(bivinzcol);//1
    bivinzcol=NULL;
    free(bivinzrow);//1
    bivinzrow=NULL;
    free(bivinzrow1);//1
    bivinzrow1=NULL;
    free(bivinzcol1);//1
    bivinzcol1=NULL;
    free(vecbivi);//1
    vecbivi=NULL;
  ldsize=ceil((laDi/100.0)*nz1max);
  ldsize+=10;
  irn1=realloc(irn1,ldsize*sizeof(int));
  if(irn1==NULL)printf("Error!!!! No memory left for irn1!\n");
  jcn1=realloc(jcn1,ldsize*sizeof(int));
  if(jcn1==NULL)printf("Error!!!! No memory left for jcn1!\n");
  irn=realloc(irn,nrowmax*sizeof(int));
  if(irn==NULL)printf("Error!!!! No memory left for irn!\n");
  jcn=realloc(jcn,ncolmax*sizeof(int));
  if(jcn==NULL)printf("Error!!!! No memory left for jcn!\n");
  keep=(int *) malloc ((nrowmax+5*ncolmax+4*ncolmax+7)*sizeof(int));//ICNTL(6)=1
  if(keep==NULL)printf("Error!!!! No memory left for keep!\n");
  if(nrowmax>ncolmax) w51=(ha_cgetype *) malloc ((5*nrowmax)*sizeof(ha_cgetype));
  else w51=(ha_cgetype *) malloc ((5*ncolmax)*sizeof(ha_cgetype));
  if(w51==NULL)printf("Error!!!! No memory left for w51!\n");
  iw51=(int *) malloc ((6*nrowmax+3*ncolmax)*sizeof(int));
  if(iw51==NULL)printf("Error!!!! No memory left for iw51!\n");
  vecbivi=realloc(vecbivi,ldsize*sizeof(ha_cgetype));
  if(vecbivi==NULL)printf("Error!!!! No memory left for vecbivi!\n");
  int *insizeda=(int *) calloc (5+nreg*insizes,sizeof(int));
  #pragma omp for schedule (static)
  for(j3=0; j3<nmatint; j3++) {
    j4=nreg+j3*(nreg+1);
    if(j3<10)strcpy(j1name,"000");
    if(j3<100&&j3>9)strcpy(j1name,"00");
    if(j3<1000&&j3>99)strcpy(j1name,"0");
    if(j3>=1000)j1name[0]='\0';
    sprintf(filename, "%d",j3);
    strcat(j1name,filename);
    strcpy(filename,"_rank");
    strcat(filename,rankname);
    strcat(filename,j1name);
    strcat(filename,".bin");
    if((presolfile=fopen(filename, "r"))==NULL) {
      printf("Cannot open file. Please run presol by setting -presol 1!\n");
    }
    fwrt=fread(insizeda, sizeof(int), 5+nreg*insizes, presolfile);
    fclose(presolfile);
    remove(filename);
    if( fwrt== 0)printf("File read error. No presol!\n");
    nrow=insizeda[0];
    ncol=insizeda[1];
    nz1=insizeda[2];
    nz=nz1;
      strcpy(filename,"_bivi");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "r"))==NULL) {
        printf("Cannot open file.\n");
      }
      fwrt=fread(vecbivi, sizeof(ha_cgetype), insizeda[2], presolfile);
      if(fwrt== 0) printf("File read error %ld %d.\n",fwrt,insizeda[2]);
      fclose(presolfile);
      strcpy(filename,"_rbvi");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "r"))==NULL) {
        printf("Cannot open file.\n");
      }
      fwrt=fread(irn1, sizeof(int), insizeda[2], presolfile);
      if(fwrt== 0) printf("File read error.\n");
      fclose(presolfile);
      strcpy(filename,"_cbvi");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "r"))==NULL) {
        printf("Cannot open file.\n");
      }
      fwrt=fread(jcn1, sizeof(int), insizeda[2], presolfile);
      if(fwrt== 0) printf("File read error.");
      fclose(presolfile);

    insize[j4*insizes]=nrow;
    insize[j4*insizes+1]=ncol;
    insize[j4*insizes+2]=nz;
    insize[j4*insizes+10]=rank;
    insize[j4*insizes+11]=j4;
    insize[j4*insizes+9]=laDi;
    insize[j4*insizes+15]=7+windx;//7+rank+jthrd+windx;//(j4+jthrd)%90+7;
    windx++;
    if(windx==eindx)windx=bindx;
      int *insized=(int *) calloc (6,sizeof(int));
      insized[0]=nrow;
      insized[1]=ncol;
      insized[2]=nz;
      insized[4]=laDi;
    ldsize=ceil((laDi/100.0)*nz);
      insized[5]=ldsize;
      if(cntl6==0&&SORD==0)cntl6in=0.3;
      else cntl6in=cntl6;
      spec51m_rank_(insized,&cntl6in,irn1,jcn1,vecbivi,irn,jcn,keep,w51,iw51);
      insized=realloc(insized,(5+nreg*insizes)*sizeof(int));
      for(i=0; i<nreg; i++) {
        for(j=0; j<insizes; j++) {
          insized[i*insizes+5+j]=insize[(i+j3*(nreg+1))*insizes+j];
        }
      }
      if(j3<10)strcpy(j1name,"000");
      if(j3<100&&j3>9)strcpy(j1name,"00");
      if(j3<1000&&j3>99)strcpy(j1name,"0");
      if(j3>=1000)j1name[0]='\0';
      sprintf(filename, "%d",j3);
      strcat(j1name,filename);
      strcpy(filename,"_rank");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "wb"))==NULL) {
        printf("Cannot open file.\n");
      }
      printf("Rank %d RRRRR! ncol %d nrow %d rank %d nz %d cntl6in %lf\n",rank,nrow,ncol,insized[3],nz,cntl6in);
      fwrt=fwrite(insized, sizeof(int), 5+nreg*insizes, presolfile);
      if(fwrt==0)printf("File write error.");
      fclose(presolfile);
      strcpy(filename,"_row");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "wb"))==NULL) {
        printf("Cannot open file.\n");
      }
      fwrt=fwrite(irn, sizeof(int), nrow, presolfile);
      if(fwrt== 0) printf("File write error.");
      fclose(presolfile);
      strcpy(filename,"_col");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "wb"))==NULL) {
        printf("Cannot open file.\n");
      }
      fwrt=fwrite(jcn, sizeof(int), ncol, presolfile);
      if(fwrt== 0) printf("File write error.");
      fclose(presolfile);
      free(insized);//1
    insize[j4*insizes+15]=0;
  }
    free(vecbivi);//1
    free(irn1);//1
    free(jcn1);//1
    free(irn);//1
    free(jcn);//1
    free(keep);
    free(w51);
    free(iw51);
  free(insizeda);
  }
  ierr = PetscFree(submatAij);//1
  printf("Rank %d End of preparation!\n",rank);
  MPI_Barrier(PETSC_COMM_WORLD);

  for(i=0; i<nmatin; i++) {
    ierr = PetscFree(submatCij[i]);//1
    ierr = PetscFree(submatBij[i]);//1
  }
  free(submatCij);//1
  free(submatBij);//1

  free(insize);//1
  return 0;
}


int HaNDBBDParSol(Mat A, Vec b, ha_cgetype *x1, uvadd VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int *ha_rows,int *ha_cols, uvadd ndblock,uvadd nreg,uvadd ntime,int *ha_ndblocks, uvadd *countvarintra1, uvadd *counteq, uvadd *counteqnoadd,uvdim laA,uvdim laDi,uvdim laD,PetscReal cntl3,PetscReal cntl6,PetscBool presol) {//,bool iter
  IS *rowindices=NULL,*colindices=NULL,*rowindicesbc=NULL,*colindicesbc1=NULL,*colindicesbc2=NULL,*Cindices=NULL,*Bindices=NULL,Cindicesc,Bindicesc,*BBindices=NULL,*BBiindices=NULL,*rowBBij=NULL,*colBBij=NULL;//,*colindicesbcpm,*colindicesbcpm1
  const PetscInt *nindices;
  PetscInt bfirst,bend,sblockin,nmatin,nmatinplus,nmatint,nmatinplust,nmatminust,nrowcolin,sumrowcolin,i,i1,j,j0,j1,j2,j3,j4,j5,j6,j7,l0,l1,l2,l3,l4,l5,rank1,proc1=0,nnzmax,j1nz,j1irnbs;
  Mat *submatAij=NULL,*submatC=NULL,*submatD=NULL,*submatBB=NULL,*submatBBij=NULL;//,*submatCij,*submatBij;,*submatB
  Vec vecxd;
  PetscInt *ai,*aj,*a1i,*a1j,*a2i,*a2j,*aic,*ajc;
  PetscReal cntl6in;
  uvadd ldsize;
  long int *biviindx1,*biviindx0,j7l;
  long int *obiviindx1= NULL;//(long int*)calloc(1,sizeof(long int));
  long int *obiviindx0= NULL;//(long int*)calloc(1,sizeof(long int));
  size_t freadresult,fwrt,frd;
  PetscInt nrow,ncol,nz,nrowc,ncolc,nrowb,ncolb,nzc,nzv;
  MPI_Status   status;
  clock_t timestr,timeend;
  int j2int,la1,colcut,rowcut,jthrd,nthrd=1,unequal=0;
  nthrd=omp_get_max_threads( );
  bool ifremove=false,iterstop;
  PetscScalar *vals,*val1s,*val2s,*valsc,vecval;
  ha_cgetype vval;
  PetscErrorCode ierr;
  PetscViewer viewer;
  FILE* fp1,*fp2,*fp3,*presolfile;
  int *begblock= (int *) calloc (mpisize,sizeof(int));
  int *ndblockinrank= (int *) calloc (mpisize,sizeof(int));
  char filename[1024],rankname[1024],j1name[1024],tempchar[1024];//,fn1[1024],fn2[1024],fn3[1024]
  if(rank<10)strcpy(rankname,"000");
  if(rank<100&&rank>9)strcpy(rankname,"00");
  if(rank<1000&&rank>99)strcpy(rankname,"0");
  if(rank>=1000)rankname[0]='\0';
  sprintf(filename, "%d",rank);
  strcat(rankname,filename);
  int nmatinBB=1;
  //reg matrix in one process only
  nmatint=(uvadd)ntime/mpisize;
  if(nmatint*mpisize<ntime)unequal=1;
  nmatminust=nmatint;
  if(nmatint*mpisize<ntime)nmatinplust=nmatint+1;
  else nmatinplust=nmatint;
  for(i=0; i<mpisize; i++)if(rank+1<=ntime-mpisize*nmatint)nmatint++;
  nmatin=(nreg+1)*nmatint;
  nmatinplus=(nreg+1)*nmatinplust;
  begblock[rank]=nmatin;
  printf("rank %d nmatin %d nmatint %d nmplus %d\n",rank,nmatin,nmatint,nmatinplus);
  for(i=0; i<mpisize; i++) {
    j=rank;
    j2=j;
    MPI_Bcast(&j2,1, MPI_INT,i, PETSC_COMM_WORLD);
    j=j2;
    MPI_Bcast(&begblock[j2],1, MPI_INT,i, PETSC_COMM_WORLD);
  }
  j2=begblock[0];
  begblock[0]=0;
  for(i=1; i<mpisize; i++) {
    j1=begblock[i];
    begblock[i]=begblock[i-1]+j2;
    j2=j1;
  }
  for(i=0; i<mpisize; i++)printf("rank %d beg block %d\n",rank,begblock[i]);
  j2=nmatin;
  MPI_Bcast(&j2,1, MPI_INT,0, PETSC_COMM_WORLD);
  ndblockinrank[0]=j2;
  for(i=1; i<mpisize; i++) {
    j2=nmatin;
    MPI_Bcast(&j2,1, MPI_INT,i, PETSC_COMM_WORLD);
    ndblockinrank[i]=ndblockinrank[i-1]+j2;
  }
  for(i=0; i<ndblock; i++) {
    j2=i-ndblockinrank[rank];
    if(j2<0&&j2>=-nmatin) {
      j1=i/(nreg+1);
      j=i-(nreg+1)*j1;
      if(rank==0)j3=0;
      else j3=ndblockinrank[rank-1];
    }
  }
  j1=0;
  for(j=0; j<mpisize; j++)if(j<rank)j1+=ndblockinrank[j];
  free(ndblockinrank);
  ierr = PetscMalloc(nmatin*sizeof(IS **),&rowindices);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatint*sizeof(IS **),&rowindicesbc);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatin*sizeof(IS **),&colindices);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatint*sizeof(IS **),&colindicesbc1);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatint*sizeof(IS **),&colindicesbc2);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatint*sizeof(IS **),&Cindices);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatint*sizeof(IS **),&Bindices);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatinBB*sizeof(IS **),&BBindices);
  CHKERRQ(ierr);
  Mat** submatBij = (Mat**)calloc(nmatin,sizeof(Mat*));
  Mat** submatCij = (Mat**)calloc(nmatin,sizeof(Mat*));
  Mat** submatB1 = (Mat**)calloc(nmatint,sizeof(Mat*));
  for (i=0; i<nmatint; i++) {
    submatB1[i] = (Mat*)calloc(1,sizeof(Mat));
  }
  Mat** submatB2 = (Mat**)calloc(nmatint,sizeof(Mat*));
  for (i=0; i<nmatint; i++) {
    submatB2[i] = (Mat*)calloc(1,sizeof(Mat));
  }

  ha_cgetype **yi2= (ha_cgetype**)calloc(nmatin,sizeof(ha_cgetype*));
  for (i=0; i<nmatin; i++) {
    yi2[i] = (ha_cgetype*)calloc(ha_ndblocks[i+begblock[rank]],sizeof(ha_cgetype));
  }
  sumrowcolin=0;
  for(j=0; j<ndblock; j++) {
    sumrowcolin+=ha_ndblocks[j];
  }
  printf("Vecs %ld sumrow %d!!!\n",VecSize,sumrowcolin);
  PetscInt *offblock= (PetscInt *) calloc (ndblock+1,sizeof(PetscInt));
  PetscInt *offblockrow= (PetscInt *) calloc (ndblock+1,sizeof(PetscInt));
  offblock[0]=0;
  for(j=0; j<ndblock; j++)offblock[j+1]=offblock[j]+countvarintra1[j+1]-countvarintra1[j]-ha_ndblocks[j];
  offblockrow[0]=0;
  for(j=0; j<ndblock; j++)offblockrow[j+1]=offblockrow[j]+counteqnoadd[j]-ha_ndblocks[j];

  char **fn01= (char**)calloc(nmatinplus,sizeof(char*));
  for (i=0; i<nmatinplus; i++) fn01[i] = (char*)calloc(1024,sizeof(char));
  char **fn02= (char**)calloc(nmatinplus,sizeof(char*));
  for (i=0; i<nmatinplus; i++) fn02[i] = (char*)calloc(1024,sizeof(char));
  char **fn03= (char**)calloc(nmatinplus,sizeof(char*));
  for (i=0; i<nmatinplus; i++) fn03[i] = (char*)calloc(1024,sizeof(char));

  timestr=clock();
  for(j1=0; j1<nmatinplus; j1++) {
    if(j1<nmatin) {
      if(j1<10)strcpy(j1name,"000");
      if(j1<100&&j1>9)strcpy(j1name,"00");
      if(j1<1000&&j1>99)strcpy(j1name,"0");
      if(j1>=1000)j1name[0]='\0';
      sprintf(filename, "%d",j1);
      strcat(j1name,filename);
      strcpy(filename,"_irnv");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      strcpy(fn01[j1],filename);
      strcpy(filename,"_keep");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      strcpy(fn02[j1],filename);
      strcpy(filename,"_vav");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      strcpy(fn03[j1],filename);
    }
  }
  printf("rank1w2 %d\n",rank);
  for(j=0; j<mpisize; j++) {
    for(j3=0; j3<nmatinplus; j3++) {
      if(j3>=nmatin)i=nmatin-1;
      else i=j3;
      bfirst=counteq[i+begblock[rank]];
      bend=ha_ndblocks[i+begblock[rank]];
      MPI_Bcast(&bfirst,1, MPI_INT,j, PETSC_COMM_WORLD);
      MPI_Bcast(&bend,1, MPI_INT,j, PETSC_COMM_WORLD);
      PetscInt *indices= (PetscInt *) calloc (bend,sizeof(PetscInt));
      ha_cgetype *yi0 = (ha_cgetype*)calloc(bend,sizeof(ha_cgetype));
      for(j1=0; j1<bend; j1++) { //ha_ndblocks[i+begblock[j]]
        j2=ha_rows[bfirst+j1];//+bfirst;
        if(j2>=Istart&&j2<Iend)indices[j1]=j2;
        else indices[j1]=-1;
      }
      VecGetValues(b,bend,indices,yi0);
      HaReduce(yi0,(forint)bend,mpisize,rank,j);
      if(rank==j&&i<nmatin) {
        memcpy (yi2[i],yi0,bend*sizeof(ha_cgetype));
      }
      free(indices);
      free(yi0);
    }
  }
  printf("rank1w2 %d\n",rank);
  ha_cgetype **yi1= (ha_cgetype**)calloc(nmatint,sizeof(ha_cgetype*));
  for (i=0; i<nmatint; i++) {
    j1=0;
    for (j=0; j<nreg+1; j++) {
      j1+=ha_ndblocks[i*(nreg+1)+j+begblock[rank]];
    }
    yi1[i] = (ha_cgetype*)calloc(j1,sizeof(ha_cgetype));
  }
  for (i=0; i<nmatint; i++) {
    j1=0;
    for (j=0; j<nreg+1; j++) {
      memcpy (&yi1[i][j1],yi2[i*(nreg+1)+j],ha_ndblocks[i*(nreg+1)+j+begblock[rank]]*sizeof(ha_cgetype));
      free(yi2[i*(nreg+1)+j]);
      j1+=ha_ndblocks[i*(nreg+1)+j+begblock[rank]];
    }
  }
  free(yi2);
  j1=0;
  for(i=0; i<nmatin; i++) {
    bfirst=counteq[i+begblock[rank]];
    PetscInt *indices= (PetscInt *) calloc (ha_ndblocks[i+begblock[rank]],sizeof(PetscInt));
    for(j=0; j<ha_ndblocks[i+begblock[rank]]; j++) {
      indices[j]=ha_rows[bfirst+j];//+bfirst;
    }

    ISCreateGeneral(PETSC_COMM_SELF,ha_ndblocks[i+begblock[rank]],indices,PETSC_COPY_VALUES,rowindices+i);
    bfirst=countvarintra1[i+begblock[rank]];
    for(j=0; j<ha_ndblocks[i+begblock[rank]]; j++) {
      indices[j]=ha_cols[bfirst+j];//+bfirst;
    }
    ISCreateGeneral(PETSC_COMM_SELF,ha_ndblocks[i+begblock[rank]],indices,PETSC_COPY_VALUES,colindices+i);
    free(indices);
  }
  for(i=0; i<nmatint; i++) {
    j2=0;
    for(j1=0; j1<nreg+1; j1++)j2+=ha_ndblocks[i*(nreg+1)+j1+begblock[rank]];
    PetscInt *indices= (PetscInt *) calloc (j2,sizeof(PetscInt));
    j2=0;
    for(j1=0; j1<nreg+1; j1++) {
      bfirst=counteq[i*(nreg+1)+j1+begblock[rank]];
      for(j=0; j<ha_ndblocks[i*(nreg+1)+j1+begblock[rank]]; j++) {
        indices[j2]=ha_rows[bfirst+j];//+bfirst;
        j2++;
      }
    }
    ISCreateGeneral(PETSC_COMM_SELF,j2,indices,PETSC_COPY_VALUES,rowindicesbc+i);
    j2=0;
    for(j1=0; j1<nreg; j1++) {
      bfirst=countvarintra1[i*(nreg+1)+j1+begblock[rank]];
      for(j=0; j<ha_ndblocks[i*(nreg+1)+j1+begblock[rank]]; j++) {
        indices[j2]=ha_cols[bfirst+j];//+bfirst;
        j2++;
      }
    }
    ISCreateGeneral(PETSC_COMM_SELF,j2,indices,PETSC_COPY_VALUES,colindicesbc1+i);
    j2=0;
    for(j1=nreg; j1<nreg+1; j1++) {
      bfirst=countvarintra1[i*(nreg+1)+j1+begblock[rank]];
      for(j=0; j<ha_ndblocks[i*(nreg+1)+j1+begblock[rank]]; j++) {
        indices[j2]=ha_cols[bfirst+j];//+bfirst;
        j2++;
      }
    }
    ISCreateGeneral(PETSC_COMM_SELF,j2,indices,PETSC_COPY_VALUES,colindicesbc2+i);
    free(indices);
  }
  PetscInt *indicesC= (PetscInt *) calloc (VecSize-sumrowcolin,sizeof(PetscInt));
  PetscInt *indicesB= (PetscInt *) calloc (VecSize-sumrowcolin,sizeof(PetscInt));
  for(i=0; i<ndblock; i++) {
    bend=countvarintra1[i+1];
    for(j=offblock[i]; j<offblock[i+1]; j++) {
      j1=j-offblock[i];
      indicesC[j]=ha_cols[bend-j1-1];//+bfirst;
    }
  }
  j1=offblock[ndblock];
  for(i=countvarintra1[ndblock]; i<VecSize; i++) {
    indicesC[j1]=i;
    j1++;
  }
  for(i=0; i<ndblock; i++) {
    bend=counteq[i]+counteqnoadd[i];
    for(j=offblockrow[i]; j<offblockrow[i+1]; j++) {
      j1=j-offblockrow[i];
      indicesB[j]=ha_rows[bend-j1-1];//+bfirst;
    }
  }
  j1=offblockrow[ndblock];
  for(i=counteq[ndblock-1]+counteqnoadd[ndblock-1]; i<VecSize; i++) {
    indicesB[j1]=i;
    j1++;
  }
  for(i=0; i<VecSize-sumrowcolin; i++)if(indicesB[i]>=Istart) {
      j=i;
      break;
    }
  for(j2=i; j2<VecSize-sumrowcolin; j2++)if(indicesB[j2]>=Iend) {
      j1=j2;
      break;
    }
  if(j2==VecSize-sumrowcolin)j1=VecSize-sumrowcolin;
  printf("rank %d j %d j1 %d istart %d iend %d\n",rank,j,j1,Istart,Iend);
  ha_cgetype *vecbiui= (ha_cgetype *) calloc (VecSize-sumrowcolin,sizeof(ha_cgetype));
  VecGetValues(b,j1-j,&indicesB[j],&vecbiui[j]);//implicite yd save mem
  ierr = VecDestroy(&b);
  CHKERRQ(ierr);
  ISCreateGeneral(PETSC_COMM_SELF,VecSize-sumrowcolin,indicesC,PETSC_COPY_VALUES,&Cindicesc);
  ISCreateGeneral(PETSC_COMM_SELF,VecSize-sumrowcolin,indicesB,PETSC_COPY_VALUES,&Bindicesc);
  for(i=0; i<nmatint; i++) {
    Cindices[i]=Cindicesc;
    Bindices[i]=Bindicesc;
  }
  free(indicesB);
  free(offblock);
  free(offblockrow);
  PetscInt *indexBB= (PetscInt *) calloc (VecSize,sizeof(PetscInt));
  for(i=0; i<VecSize; i++)indexBB[i]=i;
  ISCreateGeneral(PETSC_COMM_SELF,VecSize,indexBB,PETSC_COPY_VALUES,BBindices);
  printf("Begin preparation rank %d\n",rank);
  int insizes=17;
  int *insize=(int *) calloc (insizes*nmatin,sizeof(int));
  printf("rank %d A\n",rank);
  PetscInt *indicesbbij=(int *) calloc (1,sizeof(PetscInt));
  ierr = PetscMalloc(1*sizeof(IS **),&rowBBij);
  CHKERRQ(ierr);
  ierr = PetscMalloc(1*sizeof(IS **),&colBBij);
  CHKERRQ(ierr);

  for(i=0; i<nmatminust; i++) {
    j2=i*(nreg+1)+nreg;
    ierr = MatCreateSubMatrices(A,1,rowindices+j2,colindicesbc1+i,MAT_INITIAL_MATRIX,&submatBBij);
    nrow=submatBBij[0]->rmap->n;
    ncol=submatBBij[0]->cmap->n;
    if(nrow>ncol) {
      indicesbbij=realloc(indicesbbij,nrow*sizeof(PetscInt));
      for(j=0; j<nrow; j++)indicesbbij[j]=j;
    }
    else {
      indicesbbij=realloc(indicesbbij,ncol*sizeof(PetscInt));
      for(j=0; j<ncol; j++)indicesbbij[j]=j;
    }
    ISCreateGeneral(PETSC_COMM_SELF,nrow,indicesbbij,PETSC_COPY_VALUES,rowBBij);
    j3=0;
    for(j=0; j<nreg; j++) {
      j1=i*(nreg+1)+j;
      ierr = MatCreateSubMatrices(A,1,rowindices+j1,colindices+j2,MAT_INITIAL_MATRIX,&submatCij[j1]);
      ISCreateGeneral(PETSC_COMM_SELF,ha_ndblocks[j1+begblock[rank]],indicesbbij+j3,PETSC_COPY_VALUES,colBBij);
      ierr = MatCreateSubMatrices(submatBBij[0],1,rowBBij,colBBij,MAT_INITIAL_MATRIX,&submatBij[j1]);
      ierr = ISDestroy(&colBBij[0]);
      j3+=ha_ndblocks[j1+begblock[rank]];
    }
    ierr = MatDestroy(&submatBBij[0]);
    ierr = PetscFree(submatBBij);//1
    submatBBij=NULL;
    ierr = ISDestroy(&rowBBij[0]);
  }
  if(unequal==1){
  if(nmatint<nmatinplust) for(i=0; i<nreg; i++) {//to create again
      MatDestroy(&submatCij[(nmatint-1)*(nreg+1)+i][0]);
      ierr = PetscFree(submatCij[(nmatint-1)*(nreg+1)+i]);
      MatDestroy(&submatBij[(nmatint-1)*(nreg+1)+i][0]);
      ierr = PetscFree(submatBij[(nmatint-1)*(nreg+1)+i]);
    }
  j2=(nmatint-1)*(nreg+1)+nreg;
  ierr = MatCreateSubMatrices(A,1,rowindices+j2,colindicesbc1+(nmatint-1),MAT_INITIAL_MATRIX,&submatBBij);
  nrow=submatBBij[0]->rmap->n;
  ncol=submatBBij[0]->cmap->n;
  if(nrow>ncol) {
    indicesbbij=realloc(indicesbbij,nrow*sizeof(PetscInt));
    for(j=0; j<nrow; j++)indicesbbij[j]=j;
  }
  else {
    indicesbbij=realloc(indicesbbij,ncol*sizeof(PetscInt));
    for(j=0; j<ncol; j++)indicesbbij[j]=j;
  }
  ISCreateGeneral(PETSC_COMM_SELF,nrow,indicesbbij,PETSC_COPY_VALUES,rowBBij);
  j3=0;
  for(j=0; j<nreg; j++) {
    j1=(nmatint-1)*(nreg+1)+j;
    ierr = MatCreateSubMatrices(A,1,rowindices+j1,colindices+j2,MAT_INITIAL_MATRIX,&submatCij[j1]);
    ISCreateGeneral(PETSC_COMM_SELF,ha_ndblocks[j1+begblock[rank]],indicesbbij+j3,PETSC_COPY_VALUES,colBBij);
    ierr = MatCreateSubMatrices(submatBBij[0],1,rowBBij,colBBij,MAT_INITIAL_MATRIX,&submatBij[j1]);
    ierr = ISDestroy(&colBBij[0]);
    j3+=ha_ndblocks[j1+begblock[rank]];
  }
  ierr = MatDestroy(&submatBBij[0]);
  ierr = PetscFree(submatBBij);
  submatBBij=NULL;
  ierr = ISDestroy(&rowBBij[0]);
  }
  free(begblock);
  free(indicesbbij);
  PetscFree(rowBBij);
  PetscFree(colBBij);

  printf("rank %d A\n",rank);
  ierr = MatCreateSubMatrices(A,nmatint,rowindicesbc,Cindices,MAT_INITIAL_MATRIX,&submatC);
  printf("C\n");
  ierr = MatCreateSubMatrices(A,nmatinBB,Bindices,BBindices,MAT_INITIAL_MATRIX,&submatBB);
  ierr = MatDestroy(&A);
  CHKERRQ(ierr);
  PetscInt BBrow,BBcol;//,CCrow,CCcol;
  MatGetSize(submatBB[0],&BBrow,&BBcol);
  printf("rank %d bbrow %d bbcol %d\n",rank,BBrow,BBcol);//,CCrow,CCcol);
  PetscInt *indexBBi= (PetscInt *) calloc (BBrow,sizeof(PetscInt));
  for(i=0; i<BBrow; i++)indexBBi[i]=i;
  ierr = ISDestroy(&Bindicesc);
  ISCreateGeneral(PETSC_COMM_SELF,BBrow,indexBBi,PETSC_COPY_VALUES,&Bindicesc);
  for(i=0; i<nmatint; i++) {
    Bindices[i]=Bindicesc;
  }
  for(i=0; i<nmatint; i++) {
    ierr = MatCreateSubMatrix(submatBB[0],Bindices[i],colindicesbc1[i],MAT_INITIAL_MATRIX,&submatB1[i][0]);
    ierr = MatCreateSubMatrix(submatBB[0],Bindices[i],colindicesbc2[i],MAT_INITIAL_MATRIX,&submatB2[i][0]);
  }
  printf("B\n");
  if(rank==mpisize-1)proc1=1;
  printf("D\n");
  MatCreateSubMatrices(submatBB[0],proc1,Bindices,Cindices,MAT_INITIAL_MATRIX,&submatD);
  printf("End partitionning rank %d time %f\n",rank,((double)(clock()-timestr))/CLOCKS_PER_SEC);
  for (i=0; i<nmatin; i++) {
    ierr = ISDestroy(&rowindices[i]);
    CHKERRQ(ierr);
    ierr = ISDestroy(&colindices[i]);
    CHKERRQ(ierr);
  }
  ierr = ISDestroy(&Cindicesc);
  CHKERRQ(ierr);
  ierr = ISDestroy(&Bindicesc);
  CHKERRQ(ierr);
  PetscFree(colindices);
  ierr = ISDestroy(&BBindices[0]);
  PetscFree(rowindices);

  for (i=0; i<nmatint; i++) {
    ierr = ISDestroy(&rowindicesbc[i]);
    CHKERRQ(ierr);
  }
  PetscFree(rowindicesbc);

  PetscFree(Cindices);
  PetscFree(Bindices);
  PetscFree(BBindices);
  PetscFree(indexBB);
  PetscFree(indexBBi);
  ierr = MatDestroy(&submatBB[0]);
  PetscFree(submatBB);
  CHKERRQ(ierr);
  printf("OK1 rank %d sumrowcolin %d\n",rank,sumrowcolin);
  long int *bivinzrow=NULL;//(long int *) calloc (1,sizeof(long int));
  PetscInt *bivinzcol=NULL;//(PetscInt *) calloc (1,sizeof(PetscInt));
  PetscReal *ccolnorms=NULL;//(PetscReal *) calloc (1,sizeof(PetscReal));
  int bivirowsize=1,bivicolsize=1,bbrowij,ddrowi;
  long int vecbivisize=0,li,lj;//,halfvec;
  ha_cgetype *xi1point;
  long int *bivinzrow1=NULL;//(long int *) calloc (1,sizeof(long int));
  long int *bivinzcol1=NULL;//(long int *) calloc (1,sizeof(long int));
  uvadd xi1indx=0;
  long int nz0=0,nz1,nz3;//,nz2,halfj2;
  int fd1,fd2,fd3,frrsl1,frrsl2,frrsl3;
  omp_set_num_threads(smallthreads);
  #pragma omp parallel private(jthrd,nthrd,j4,j1name,filename,frd,fd1,nrow,ncol,i,j,presolfile,nz,nz1,cntl6in,fwrt,ldsize,j3) shared(insize)
  {
  int windx=0,bindx,eindx;
  long int insizeda0=0,insizeda1=0,insizeda2=0;
  jthrd=omp_get_thread_num();
  nthrd=omp_get_max_threads( );
  nthrd=90/nthrd;
  bindx=nthrd*jthrd;
  eindx=bindx+nthrd;
  windx=bindx;
  int *insizeda=(int *) calloc (5+nreg*insizes,sizeof(int));
  #pragma omp for  schedule (static)
  for(j3=0; j3<nmatint; j3++) {
    if(j3<10)strcpy(j1name,"000");
    if(j3<100&&j3>9)strcpy(j1name,"00");
    if(j3<1000&&j3>99)strcpy(j1name,"0");
    if(j3>=1000)j1name[0]='\0';
    sprintf(filename, "%d",j3);
    strcat(j1name,filename);
    strcpy(filename,"_rank");
    strcat(filename,rankname);
    strcat(filename,j1name);
    strcat(filename,".bin");
    #pragma omp critical
    {
    if((presolfile=fopen(filename, "r"))==NULL) {
      printf("Cannot open file. Please run presol by setting -presol 1!\n");
    }
    frd=fread(insizeda, sizeof(int), 5+nreg*insizes, presolfile);
    fclose(presolfile);
    if( frd== 0)printf("File read error. No presol!\n");
    }
    if(insizeda0<insizeda[0])insizeda0=insizeda[0];
    if(insizeda1<insizeda[1])insizeda1=insizeda[1];
    if(insizeda2<insizeda[2])insizeda2=insizeda[2];
  }
  ldsize=ceil((laDi/100.0)*insizeda2);
  ldsize+=10;
  insizeda2=ldsize;
  int *irn1=(int *) malloc(ldsize*sizeof(int));
  if(irn1==NULL)printf("Error!!!! No memory left for irn1!\n");
  int *jcn1=(int *) malloc(ldsize*sizeof(int));
  if(jcn1==NULL)printf("Error!!!! No memory left for jcn1!\n");
  ha_cgetype *vecbivi0= (ha_cgetype *) malloc (ldsize*sizeof(ha_cgetype));
  if(vecbivi0==NULL)printf("Error!!!! No memory left for vecbivi0!\n");
  int *irn1a=(int *) malloc (insizeda0*sizeof(int));
  if(irn1a==NULL)printf("Error!!!! No memory left for irn1a!\n");
  int *jcn1a=(int *) malloc (insizeda1*sizeof(int));
  if(jcn1a==NULL)printf("Error!!!! No memory left for jcn1a!\n");
  ha_cgetype *fw=(ha_cgetype *) malloc (insizeda0*sizeof(ha_cgetype));
  if(fw==NULL)printf("Error!!!! No memory left for fw!\n");
  int *fiw=(int *) malloc ((6*insizeda0+3*insizeda1)*sizeof(int));
  if(fiw==NULL)printf("Error!!!! No memory left for fiw!\n");
  int *fkeep=(int *) malloc ((insizeda0+5*insizeda1+4*insizeda1/1+7)*sizeof(int));//ICNTL(6)=1
  if(fkeep==NULL)printf("Error!!!! No memory left for fkeep!\n");
  #pragma omp for  schedule (static)
  for(j3=0; j3<nmatint; j3++) {
    j4=nreg+j3*(nreg+1);
    //read from pre
    #pragma omp critical
    {
    if(j3<10)strcpy(j1name,"000");
    if(j3<100&&j3>9)strcpy(j1name,"00");
    if(j3<1000&&j3>99)strcpy(j1name,"0");
    if(j3>=1000)j1name[0]='\0';
    sprintf(filename, "%d",j3);
    strcat(j1name,filename);
    strcpy(filename,"_rank");
    strcat(filename,rankname);
    strcat(filename,j1name);
    strcat(filename,".bin");
    if((presolfile=fopen(filename, "r"))==NULL) {
      printf("Cannot open file. Please run presol by setting -presol 1!\n");
    }
    frd=fread(insizeda, sizeof(int), 5+nreg*insizes, presolfile);
    fclose(presolfile);
    remove(filename);
    if( frd== 0)printf("File read error. No presol!\n");
    }
    nrow=insizeda[3];
    ncol=insizeda[3];
    #pragma omp critical
    {
    for(i=0; i<nreg; i++) {
      for(j=0; j<insizes; j++) {
        insize[(i+j3*(nreg+1))*insizes+j]=insizeda[i*insizes+5+j];
      }
    }
    strcpy(filename,"_row");
    strcat(filename,rankname);
    strcat(filename,j1name);
    strcat(filename,".bin");
    if((presolfile=fopen(filename, "r"))==NULL) {
      printf("Cannot open file.\n");
    }
    frd=fread(irn1a, sizeof(int), insizeda[0], presolfile);
    if(frd == 0)printf("File read error.");
    fclose(presolfile);
    remove(filename);
    strcpy(filename,"_col");
    strcat(filename,rankname);
    strcat(filename,j1name);
    strcat(filename,".bin");
    if((presolfile=fopen(filename, "r"))==NULL) {
      printf("Cannot open file.\n");
    }
    frd=fread(jcn1a, sizeof(int), insizeda[1], presolfile);
    if(frd== 0) printf("File read error.");
    fclose(presolfile);
    remove(filename);
    if(isLinux==0) {
      strcpy(filename,"_bivi");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "r"))==NULL) {
        printf("Cannot open file.\n");
      }
      frd=fread(vecbivi0, sizeof(ha_cgetype), insizeda[2], presolfile);
      if(frd== 0) printf("File read error %ld %d.\n",frd,insizeda[2]);
      fclose(presolfile);
      remove(filename);
      strcpy(filename,"_rbvi");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "r"))==NULL) {
        printf("Cannot open file.\n");
      }
      frd=fread(irn1, sizeof(int), insizeda[2], presolfile);
      if(frd== 0) printf("File read error.\n");
      fclose(presolfile);
      remove(filename);
      strcpy(filename,"_cbvi");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "r"))==NULL) {
        printf("Cannot open file.\n");
      }
      frd=fread(jcn1, sizeof(int), insizeda[2], presolfile);
      if(frd== 0) printf("File read error.");
      fclose(presolfile);
      remove(filename);
    }
    else {
      strcpy(filename,"_bivi");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      fd1 = open(filename, O_RDONLY);
      frd = read(fd1,vecbivi0, sizeof(ha_cgetype)*insizeda[2]);
      close(fd1);
      remove(filename);
      strcpy(filename,"_rbvi");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      fd1 = open(filename, O_RDONLY);
      frd = read(fd1,irn1, sizeof(int)*insizeda[2]);
      close(fd1);
      remove(filename);
      strcpy(filename,"_cbvi");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      fd1 = open(filename, O_RDONLY);
      frd = read(fd1,jcn1, sizeof(int)*insizeda[2]);
      close(fd1);
      remove(filename);
    }
    }
   
    if(insizeda[0]>insizeda[3]||insizeda[1]>insizeda[3]) {
      patio_mat_(insizeda,irn1,jcn1,vecbivi0,irn1a,jcn1a);
    }
    nz1=insizeda[2];
    //end read from pre
    nz=nz1;
    ldsize=ceil((laDi/100.0)*nz1);
    #pragma omp critical
    {
    insize[j4*insizes]=nrow;
    insize[j4*insizes+1]=ncol;
    insize[j4*insizes+2]=nz1;
    insize[j4*insizes+10]=rank;
    insize[j4*insizes+11]=j4;
    insize[j4*insizes+9]=laDi;
    insize[j4*insizes+16]=ldsize;
    insize[j4*insizes+15]=7+windx;//(j3+jthrd)%90+7;
    }
    windx++;
    if(windx==eindx)windx=bindx;
      prep48_alu1_(insize+j4*insizes,irn1,jcn1,vecbivi0,fw,fiw,fkeep);
    printf("rank %d ok0 nz %ld nzmax %ld ldsize %ld!!!!!!!!\n",rank,nz1,insizeda2,ldsize);
  }
  printf("rank %d ok1!!!!!!!!\n",rank);
  free(vecbivi0);
  printf("rank %d ok2!!!!!!!!\n",rank);
  free(irn1);
  printf("rank %d ok3!!!!!!!!\n",rank);
  free(jcn1);
  printf("rank %d ok4!!!!!!!!\n",rank);
  free(irn1a);
  printf("rank %d ok5!!!!!!!!\n",rank);
  free(jcn1a);
  printf("rank %d ok6!!!!!!!!\n",rank);
  free(fw);
  printf("rank %d ok7!!!!!!!!\n",rank);
  free(fiw);
  printf("rank %d ok8!!!!!!!!\n",rank);
  free(fkeep);
  printf("rank %d ok9!!!!!!!!\n",rank);
  free(insizeda);
  printf("rank %d ok10!!!!!!!!\n",rank);
  }
  omp_set_num_threads(mymaxnumthrd);
  if(presol)
    return 0;
  printf("Rank %d End of preparation!\n",rank);
  MPI_Barrier(PETSC_COMM_WORLD);

  ccolnorms=realloc(ccolnorms,BBrow*sizeof(PetscReal));
  bivinzrow=realloc(bivinzrow,BBrow*sizeof(long int));
  memset(bivinzrow,0,BBrow*sizeof(long int));
  bivinzcol=realloc(bivinzcol,BBrow*sizeof(PetscInt));
  memset(bivinzcol,0,BBrow*sizeof(PetscInt));
  printf("rank %d OKKK BBrow %d\n",rank,BBrow);
  for(j1=0; j1<nmatint; j1++) {
    MatGetColumnNorms(submatC[j1],NORM_1,ccolnorms);
    for(i=0; i<BBrow; i++)if(ccolnorms[i]>0)bivinzcol[i]++;
    Mat_SeqAIJ         *a1b=(Mat_SeqAIJ*)submatB1[j1][0]->data;//*aa=subA->data;
    ai= a1b->i;
    nz=a1b->nz;
    nrow=submatB1[j1][0]->rmap->n;
    for(i=0; i<nrow-1; i++)if(ai[i]!=ai[i+1])bivinzrow[i]++;
    if(ai[i]<nz)bivinzrow[i]++;

    Mat_SeqAIJ         *a2b=(Mat_SeqAIJ*)submatB2[j1][0]->data;//*aa=subA->data;
    ai= a2b->i;
    nz=a2b->nz;
    nrow=submatB2[j1][0]->rmap->n;
    for(j=0; j<nrow-1; j++)if(ai[j]!=ai[j+1])bivinzrow[j]++;
    if(ai[j]<nz)bivinzrow[j]++;
    nz=nz+a1b->nz;
  }
  free(ccolnorms);
  bivirowsize=0;
  bivicolsize=0;
  for(i=0; i<VecSize-sumrowcolin; i++) {
    if(bivinzrow[i]>0) {
      bivinzrow[i]=bivirowsize;
      bivirowsize++;
    }
    else bivinzrow[i]=-1;
    if(bivinzcol[i]>0) {
      bivinzcol[i]=bivicolsize;
      bivicolsize++;
    }
    else bivinzcol[i]=-1;
  }
  bivinzrow1=realloc (bivinzrow1,bivirowsize*sizeof(long int));//(int *) calloc (bivirowsize,sizeof(int));
  memset(bivinzrow1,0,bivirowsize*sizeof(long int));
  bivinzcol1=realloc (bivinzcol1,bivicolsize*sizeof(long int));//(int *) calloc (bivicolsize,sizeof(int));
  memset(bivinzcol1,0,bivicolsize*sizeof(long int));
  j=0;
  j2=0;
  long int lBBrow=BBrow;
  for(li=0; li<VecSize-sumrowcolin; li++) {
    if(bivinzrow[li]>=0) {
      bivinzrow1[j]=li*lBBrow;
      j++;
      bivinzrow[li]=bivinzrow[li]*bivicolsize;//faster MSOL
    }
    if(bivinzcol[li]>=0) {
      bivinzcol1[j2]=li;
      j2++;
    }
  }
  vecbivisize=((long int)bivirowsize)*((long int)bivicolsize);
  if(vecbivisize<0)printf("Overflow detected, the interface problem is too big!!!!!!\n");
  ha_cgetype *vecbivi= (ha_cgetype *) calloc (vecbivisize,sizeof(ha_cgetype));
  printf("vecbivisize %ld rank %d\n",vecbivisize,rank);
  xi1indx=0;
  int maxcolc1=0,maxcolc=0,maxrowc=0,maxrowc1=0;
  long int longsize;
  int *vecnrowc= (int*)calloc(nmatint+1,sizeof(int));
  for(j1=0; j1<nmatint; j1++) {
    MatGetSize(submatC[j1],&maxrowc1,&maxcolc1);
    vecnrowc[j1+1]=vecnrowc[j1]+maxrowc1;
    if(maxcolc<maxcolc1)maxcolc=maxcolc1;
    if(maxcolc<maxrowc1)maxcolc=maxrowc1;
    if(maxrowc<maxrowc1)maxrowc=maxrowc1;
  }
  int maxrowcij;
  ha_cgetype *xi1 = (ha_cgetype*)calloc(sumrowcolin,sizeof(ha_cgetype));
  #pragma omp parallel private(jthrd,timestr,aic,ajc,valsc,nrowc,ncolc,a1i,a1j,val1s,nz,a2i,a2j,val2s,nrowb,ncolb,i,j,j2,xi1point,xi1indx,maxrowcij,la1,fp1,fp2,fp3,freadresult,frrsl1,frrsl2,frrsl3,longsize,nzc) shared(submatC,submatB1,submatB2,xi1,submatCij,submatBij,insize,yi1,vecbivi,vecbiui)
  {
  int icntl[20],info[20];
  ha_cgetype cntl[10],rinfo[10],error1[3];
  ha_cgetype *w= (ha_cgetype*)calloc(4*maxcolc,sizeof(ha_cgetype));
  int *iw= (int*)calloc(maxcolc,sizeof(int));
  int **irnereg= (int**)calloc(nreg+1,sizeof(int*));
  int **keepreg= (int**)calloc(nreg+1,sizeof(int*));
  ha_cgetype **valereg = (ha_cgetype**)calloc(nreg+1,sizeof(ha_cgetype*));
  FILE** nfp1= (FILE**)calloc(nreg+1,sizeof(FILE*));
  FILE** nfp2= (FILE**)calloc(nreg+1,sizeof(FILE*));
  FILE** nfp3= (FILE**)calloc(nreg+1,sizeof(FILE*));
  int* nfd1= (int*)calloc(nreg+1,sizeof(int));
  int* nfd2= (int*)calloc(nreg+1,sizeof(int));
  int* nfd3= (int*)calloc(nreg+1,sizeof(int));
  Mat submatCT;
  jthrd=omp_get_thread_num();
  #pragma omp for
  for(j1=0; j1<nmatint; j1++) {
    xi1indx=vecnrowc[j1];
    timestr=clock();//time(&timestr);//=clock();
    MatTranspose(submatC[j1],MAT_INITIAL_MATRIX,&submatCT);
    Mat_SeqAIJ         *ac=(Mat_SeqAIJ*)submatCT->data;//*aa=subA->data;
    aic= ac->i;
    ajc= ac->j;
    valsc=ac->a;
    nzc=ac->nz;
    nrowc=submatCT->rmap->n;
    ncolc=submatCT->cmap->n;
    Mat_SeqAIJ         *a1b=(Mat_SeqAIJ*)submatB1[j1][0]->data;//*aa=subA->data;
    a1i= a1b->i;//irnbs+j1irnbs;
    a1j= a1b->j;//jcnbs+j1nz;
    val1s=a1b->a;//valbs+j1nz;
    nz=a1b->nz;//nzbs[j1];
    Mat_SeqAIJ         *a2b=(Mat_SeqAIJ*)submatB2[j1][0]->data;//*aa=subA->data;
    a2i= a2b->i;//irnbs+j1irnbs;
    a2j= a2b->j;//jcnbs+j1nz;
    val2s=a2b->a;//valbs+j1nz;
    nz=nz+a2b->nz;//nzbs[j1];
    nrowb=submatB1[j1][0]->rmap->n;
    ncolb=submatB1[j1][0]->cmap->n;
    PetscInt *ai12= (PetscInt*)calloc(nrowb+1,sizeof(PetscInt));
    PetscInt *aj12= (PetscInt*)calloc(nz,sizeof(PetscInt));
    PetscScalar *vals12= (PetscScalar*)calloc(nz,sizeof(PetscScalar));
    for(i=0; i<nrowb-1; i++) {
      ai12[i]=a1i[i]+a2i[i];
      memcpy(&aj12[ai12[i]],&a1j[a1i[i]],(a1i[i+1]-a1i[i])*sizeof(PetscInt));
      memcpy(&vals12[ai12[i]],&val1s[a1i[i]],(a1i[i+1]-a1i[i])*sizeof(PetscScalar));
      j2=ai12[i]+a1i[i+1]-a1i[i];
      memcpy(&vals12[j2],&val2s[a2i[i]],(a2i[i+1]-a2i[i])*sizeof(PetscScalar));
      for(j=a2i[i]; j<a2i[i+1]; j++)aj12[j2+j-a2i[i]]=a2j[j]+ncolb;
    }
    ai12[i]=a1i[i]+a2i[i];
    ai12[i+1]=nz;
    nz=a1b->nz;
    memcpy(&aj12[ai12[i]],&a1j[a1i[i]],(nz-a1i[i])*sizeof(PetscInt));
    memcpy(&vals12[ai12[i]],&val1s[a1i[i]],(nz-a1i[i])*sizeof(PetscScalar));
    j2=ai12[i]+nz-a1i[i];
    nz=a2b->nz;
    memcpy(&vals12[j2],&val2s[a2i[i]],(nz-a2i[i])*sizeof(PetscScalar));
    for(j=a2i[i]; j<nz; j++)aj12[j2+j-a2i[i]]=a2j[j]+ncolb;
    nz=nz+a1b->nz;
    ncolb=ncolb+submatB2[j1][0]->cmap->n;
    MatDestroy(&submatB1[j1][0]);//submatBT);
    free(submatB1[j1]);
    MatDestroy(&submatB2[j1][0]);//submatBT);
    free(submatB2[j1]);
    xi1point=xi1+xi1indx;
    printf("OK here j1 %d yi1 %ld xi1indx %ld a2j %d!\n",j1,sizeof(yi1[j1])/sizeof(PetscInt),xi1indx,ncolb);
    maxrowcij=0;
    #pragma omp critical
    {
    for(i=0; i<nreg+1; i++) {
      j2=j1*(nreg+1)+i;
      insize[j2*insizes+15]=j2%90+7;
      if(i!=nreg&&(submatCij[j2][0]->rmap->n)>maxrowcij)maxrowcij=submatCij[j2][0]->rmap->n;
      la1=ceil((insize[j2*insizes+9]/100.0)*insize[j2*insizes+2]);
      if(isLinux==0) {
        nfp1[i] = fopen(fn01[j2], "rb");
        if (nfp1[i]==NULL)printf("File opening error\n");
        irnereg[i] = realloc(irnereg[i],la1*sizeof(int));//(int*)calloc(la1,sizeof(int));
        freadresult=fread(irnereg[i],sizeof(int),la1,nfp1[i]);
        fclose(nfp1[i]);
        nfp2[i] = fopen(fn02[j2], "rb");
        if (nfp2[i]==NULL)printf("File opening error\n");
        keepreg[i] = realloc(keepreg[i],insize[j2*insizes+12]*sizeof(int));//(int*)calloc(insize[j1*insizes+12],sizeof(int));
        freadresult=fread(keepreg[i],sizeof(int),insize[j2*insizes+12],nfp2[i]);
        fclose(nfp2[i]);
        nfp3[i] = fopen(fn03[j2], "rb");
        if (nfp3[i]==NULL)printf("File opening error\n");
        valereg[i] = realloc(valereg[i],la1*sizeof(ha_cgetype));//(ha_cgetype*)calloc(la1,sizeof(ha_cgetype));
        freadresult=fread(valereg[i],sizeof(ha_cgetype),la1,nfp3[i]);
        fclose(nfp3[i]);
      }
      else {
        nfd1[i] = open(fn01[j2],O_RDONLY);
        nfd2[i] = open(fn02[j2],O_RDONLY);
        nfd3[i] = open(fn03[j2],O_RDONLY);
        la1=ceil((insize[j2*insizes+9]/100.0)*insize[j2*insizes+2]);
        irnereg[i] = realloc(irnereg[i],la1*sizeof(int));//(int*)calloc(la1,sizeof(int));
        keepreg[i] = realloc(keepreg[i],insize[j2*insizes+12]*sizeof(int));//(int*)calloc(insize[j1*insizes+12],sizeof(int));
        valereg[i] = realloc(valereg[i],la1*sizeof(ha_cgetype));//(ha_cgetype*)calloc(la1,sizeof(ha_cgetype));
        frrsl1=read(nfd1[i],irnereg[i],sizeof(int)*la1);
        frrsl2=read(nfd2[i],keepreg[i],insize[j2*insizes+12]*sizeof(int));
        frrsl3=read(nfd3[i],valereg[i],la1*sizeof(ha_cgetype));
        close(nfd1[i]);
        close(nfd2[i]);
        close(nfd3[i]);
      }
    }
    }
    printf("OK here11a nrowc %d ncolc %d maxrow %d!\n",nrowc,ncolc,maxrowcij);
    ha_cgetype *b02 = (ha_cgetype*)calloc(maxrowcij,sizeof(ha_cgetype));

    NDBBD_sol_nread1(rank,j1*(nreg+1),nreg,insize,insizes,submatCij,submatBij,yi1[j1],xi1point,irnereg,keepreg,valereg,cntl,rinfo,error1,icntl,info,w,iw,b02);
    ha_cgetype *bccol= (ha_cgetype*)calloc(ncolc,sizeof(ha_cgetype));
    longsize=ncolc*sizeof(ha_cgetype);
    for(i=0; i<nrowc-1; i++) {
      if(aic[i]<aic[i+1]) {
        memset(yi1[j1],0,longsize);
        for(j=aic[i]; j<aic[i+1]; j++) {
          yi1[j1][ajc[j]]=valsc[j];
        }
        NDBBD_sol_nread1(rank,j1*(nreg+1),nreg,insize,insizes,submatCij,submatBij,yi1[j1],bccol,irnereg,keepreg,valereg,cntl,rinfo,error1,icntl,info,w,iw,b02);
        spar_vbiviadd_(bccol,bivinzcol+i,bivinzrow,&vecbivisize,&nrowb,&ncolc,&nz,ai12,aj12,vals12,vecbivi);
      }
    }
    if(aic[i]<nzc) {
      memset(yi1[j1],0,longsize);
      for(j=aic[i]; j<nzc; j++) {
        yi1[j1][ajc[j]]=valsc[j];
      }
      NDBBD_sol_nread1(rank,j1*(nreg+1),nreg,insize,insizes,submatCij,submatBij,yi1[j1],bccol,irnereg,keepreg,valereg,cntl,rinfo,error1,icntl,info,w,iw,b02);
      spar_vbiviadd_(bccol,bivinzcol+i,bivinzrow,&vecbivisize,&nrowb,&ncolc,&nz,ai12,aj12,vals12,vecbivi);
    }
    free(bccol);
    free(b02);
    MatDestroy(&submatCT);
    free(yi1[j1]);
    //Multiply Bi by ui:
    spar_mulmin_(xi1point,&nrowb,&nz,ai12,aj12,vals12,vecbiui);
    free(ai12);
    free(aj12);
    free(vals12);
    timeend=clock();//(&timeend);
    printf("Submatrix %d rank %d calculation time %f\n",j1,rank,difftime(timeend,timestr)/CLOCKS_PER_SEC);
  }
  free(w);
  free(iw);
  for(i=0; i<nreg+1; i++) {
    free(irnereg[i]);//= (int*)calloc(1,sizeof(int));
    free(keepreg[i]);//= (int*)calloc(1,sizeof(int));
    free(valereg[i]);// = (ha_cgetype*)calloc(1,sizeof(ha_cgetype));
  }
  free(irnereg);//= (int**)calloc(nreg+1,sizeof(int*));
  free(keepreg);//= (int**)calloc(nreg+1,sizeof(int*));
  free(valereg);// = (ha_cgetype**)calloc(nreg+1,sizeof(ha_cgetype*));
  free(nfp1);
  free(nfp2);
  free(nfp3);
  free(nfd1);
  free(nfd2);
  free(nfd3);
  }
  printf("Completed calculation of partitioned matrices! Rank %d\n",rank);
  free(bivinzcol);
  free(bivinzrow);
  free(submatB1);
  free(submatB2);
  CHKERRQ(ierr);
  free(yi1);
  printf("Completed calculation of partitioned matrices! Rank %d\n",rank);
  timestr=clock();

  #pragma omp parallel private(li) reduction(+:lj)
  {
  lj=0;
  #pragma omp for
  for(li=0; li<vecbivisize; ++li){
    if(vecbivi[li]!=0)lj++;
  }
  }
  obiviindx1=realloc(obiviindx1,(lj+1)*sizeof(long int));
  obiviindx1[0]=-1;
  biviindx1=obiviindx1;
  biviindx1++;
  long int lj2=0;
  printf("vbivi %ld colsize %d\n",vecbivisize,bivicolsize);
  for(li=0; li<vecbivisize; ++li) {
    if(vecbivi[li]!=0) {
      vecbivi[lj2]=vecbivi[li];
      biviindx1[lj2]=bivinzcol1[li%bivicolsize]+bivinzrow1[li/bivicolsize];//[i-j1*bivicolsize]+bivinzrow1[j1];//*nrowb;
      ++lj2;
    }
  }
  free(bivinzrow1);
  free(bivinzcol1);
  printf("OK here11! j %ld j2 %ld bivi %ld\n",lj,lj2,biviindx1[lj2-1]);
  vecbivi=realloc(vecbivi,lj2*sizeof(ha_cgetype));
  printf("OK here11!\n");
  nz0=lj2;
  int group_size,ha_id,group_size1,group_size11,ha_id1,color;
  MPI_Comm_rank( HA_COMM, &ha_id);
  MPI_Comm_size(HA_COMM,&group_size);
  if(ha_id==group_size-1)color=1;
  else color=0;
  MPI_Comm_rank( HA1_COMM, &ha_id1);
  MPI_Comm_size(HA1_COMM,&group_size1);
  group_size11=group_size1;
  MPI_Bcast(&group_size11,1, MPI_INT,mpisize-1, PETSC_COMM_WORLD);
  long int sindx01,sindx02,sindx03;
  
  if(ha_id!=group_size-1) {
    MPI_Send(&lj2,1, MPI_LONG,group_size-1, 10, HA_COMM);
    sindx01=lj2/MAXSSIZE;
    sindx03=0;
    for(sindx02=0; sindx02<sindx01; sindx02++) {
    if(SORD==1)MPI_Send(vecbivi+sindx03,MAXSSIZE, MPI_DOUBLE,group_size-1, 11*sindx01+sindx02, HA_COMM);
    else MPI_Send(vecbivi+sindx03,MAXSSIZE, MPI_FLOAT,group_size-1, 11*sindx01+sindx02, HA_COMM);
    sindx03+=MAXSSIZE;
    }
    if(SORD==1)MPI_Send(vecbivi+sindx03,lj2-sindx01*MAXSSIZE, MPI_DOUBLE,group_size-1, 11, HA_COMM);
    else MPI_Send(vecbivi+sindx03,lj2-sindx01*MAXSSIZE, MPI_FLOAT,group_size-1, 11, HA_COMM);
    sindx03=0;
    for(sindx02=0; sindx02<sindx01; sindx02++) {
    MPI_Send(biviindx1+sindx03,MAXSSIZE, MPI_LONG,group_size-1, 12*sindx01+sindx02, HA_COMM);
    sindx03+=MAXSSIZE;
    }
    MPI_Send(biviindx1+sindx03,lj2-sindx01*MAXSSIZE, MPI_LONG,group_size-1, 12, HA_COMM);
    free(obiviindx1);
    obiviindx1=NULL;
    free(vecbivi);
    vecbivi=NULL;
  }
  if(ha_id==group_size-1) {
    for(j6=0; j6<group_size-1; j6++) {
      MPI_Recv(&lj2,1, MPI_LONG,  MPI_ANY_SOURCE,10, HA_COMM,&status);
      i=status.MPI_SOURCE;
      obiviindx0=realloc(obiviindx0,(lj2+1)*sizeof(long int));
      obiviindx0[0]=-1;
      biviindx0=obiviindx0;
      biviindx0++;
      ha_cgetype *vecbivi0= (ha_cgetype *) calloc (lj2,sizeof(ha_cgetype));
    sindx01=lj2/MAXSSIZE;
    sindx03=0;
    for(sindx02=0; sindx02<sindx01; sindx02++) {
      if(SORD==1)MPI_Recv(vecbivi0+sindx03,MAXSSIZE, MPI_DOUBLE,i,11*sindx01+sindx02, HA_COMM,&status);
      else MPI_Recv(vecbivi0+sindx03,MAXSSIZE, MPI_FLOAT,i,11*sindx01+sindx02, HA_COMM,&status);
      sindx03+=MAXSSIZE;
    }
      if(SORD==1)MPI_Recv(vecbivi0+sindx03,lj2-sindx01*MAXSSIZE, MPI_DOUBLE,i,11, HA_COMM,&status);
      else MPI_Recv(vecbivi0+sindx03,lj2-sindx01*MAXSSIZE, MPI_FLOAT,i,11, HA_COMM,&status);
    sindx03=0;
    for(sindx02=0; sindx02<sindx01; sindx02++) {
      MPI_Recv(biviindx0+sindx03,MAXSSIZE, MPI_LONG,i,12*sindx01+sindx02, HA_COMM,&status);
      sindx03+=MAXSSIZE;
    }
      MPI_Recv(biviindx0+sindx03,lj2-sindx01*MAXSSIZE, MPI_LONG,i,12, HA_COMM,&status);
      my_spar_compl_(biviindx1,&nz0,biviindx0,&lj2,&nz1);
      vecbivi=realloc(vecbivi,nz1*sizeof(ha_cgetype));
      obiviindx1=realloc(obiviindx1,(nz1+1)*sizeof(long int));
      obiviindx1[0]=-1;
      biviindx1=obiviindx1;
      biviindx1++;
      my_spar_add3l_(vecbivi,obiviindx1,&nz0,vecbivi0,obiviindx0,&lj2,&nz1);
      free(obiviindx0);
      obiviindx0=NULL;
      free(vecbivi0);
      vecbivi0=NULL;
      printf("!!!!!!!!!!!! rank11 %d nz0 %ld j2 %ld nz1 %ld\n",rank,nz0,lj2,nz1);
      nz0=nz1;
    }
  }

  if(group_size11>1&&color==1) {
    nz3=nz0;
    if(ha_id1!=group_size1-1) {
      lj2=nz3;
      MPI_Send(&lj2,1, MPI_LONG,group_size1-1, 10, HA1_COMM);
    sindx01=lj2/MAXSSIZE;
    sindx03=0;
    for(sindx02=0; sindx02<sindx01; sindx02++) {
      if(SORD==1)MPI_Send(vecbivi+sindx03,MAXSSIZE, MPI_DOUBLE,group_size1-1, 11*sindx01+sindx02, HA1_COMM);
      else MPI_Send(vecbivi+sindx03,MAXSSIZE, MPI_FLOAT,group_size1-1, 11*sindx01+sindx02, HA1_COMM);
      sindx03+=MAXSSIZE;
    }
      if(SORD==1)MPI_Send(vecbivi+sindx03,lj2-sindx01*MAXSSIZE, MPI_DOUBLE,group_size1-1, 11, HA1_COMM);
      else MPI_Send(vecbivi+sindx03,lj2-sindx01*MAXSSIZE, MPI_FLOAT,group_size1-1, 11, HA1_COMM);
    sindx03=0;
    for(sindx02=0; sindx02<sindx01; sindx02++) {
      MPI_Send(biviindx1+sindx03,MAXSSIZE, MPI_LONG,group_size1-1, 12*sindx01+sindx02, HA1_COMM);
      sindx03+=MAXSSIZE;
    }
      MPI_Send(biviindx1+sindx03,lj2-sindx01*MAXSSIZE, MPI_LONG,group_size1-1, 12, HA1_COMM);
      free(obiviindx1);
      obiviindx1=NULL;
      free(vecbivi);
      vecbivi=NULL;
    }
    if(ha_id1==group_size1-1) {
      for(j6=0; j6<group_size1-1; j6++) {
        MPI_Recv(&lj2,1, MPI_LONG,  MPI_ANY_SOURCE,10, HA1_COMM,&status);
        i=status.MPI_SOURCE;
        obiviindx0=realloc(obiviindx0,(lj2+1)*sizeof(long int));
        obiviindx0[0]=-1;
        biviindx0=obiviindx0;
        biviindx0++;
        ha_cgetype *vecbivi0= (ha_cgetype *) calloc (lj2,sizeof(ha_cgetype));
    sindx01=lj2/MAXSSIZE;
    sindx03=0;
    for(sindx02=0; sindx02<sindx01; sindx02++) {
        if(SORD==1)MPI_Recv(vecbivi0+sindx03,MAXSSIZE, MPI_DOUBLE,i,11*sindx01+sindx02, HA1_COMM,&status);
        else MPI_Recv(vecbivi0+sindx03,MAXSSIZE, MPI_FLOAT,i,11*sindx01+sindx02, HA1_COMM,&status);
        sindx03+=MAXSSIZE;
    }
        if(SORD==1)MPI_Recv(vecbivi0+sindx03,lj2-sindx01*MAXSSIZE, MPI_DOUBLE,i,11, HA1_COMM,&status);
        else MPI_Recv(vecbivi0+sindx03,lj2-sindx01*MAXSSIZE, MPI_FLOAT,i,11, HA1_COMM,&status);
    sindx03=0;
    for(sindx02=0; sindx02<sindx01; sindx02++) {
        MPI_Recv(biviindx0+sindx03,MAXSSIZE, MPI_LONG,i,12*sindx01+sindx02, HA1_COMM,&status);
        sindx03+=MAXSSIZE;
    }
        MPI_Recv(biviindx0+sindx03,lj2-sindx01*MAXSSIZE, MPI_LONG,i,12, HA1_COMM,&status);
        my_spar_compl_(biviindx1,&nz0,biviindx0,&lj2,&nz1);
        vecbivi=realloc(vecbivi,nz1*sizeof(ha_cgetype));
        obiviindx1=realloc(obiviindx1,(nz1+1)*sizeof(long int));
        obiviindx1[0]=-1;
        biviindx1=obiviindx1;
        biviindx1++;
        my_spar_add3l_(vecbivi,obiviindx1,&nz0,vecbivi0,obiviindx0,&lj2,&nz1);
        free(obiviindx0);
        obiviindx0=NULL;
        free(vecbivi0);
        vecbivi0=NULL;
        printf("!!!!!!!!!!!! rank22 %d nz0 %ld j2 %ld nz1 %ld\n",rank,nz0,lj2,nz1);
        nz0=nz1;
      }
    }
  }
  forint vecbiuisize=0;
  vecbiuisize=VecSize-sumrowcolin;
  printf("after bivi rank %d\n",rank);
  HaReduceNoComp(vecbiui,vecbiuisize,mpisize,rank,mpisize-1);
  printf("Completed MPI_Reduce Operation! Rank %d time %f\n",rank,((double)(clock()-timestr))/CLOCKS_PER_SEC);
  if(rank!=mpisize-1){
    free(vecbiui);
    vecbiui=NULL;
  }
  timestr=clock();
  ha_cgetype *xd;//= (ha_cgetype*)calloc(vecbiuisize,sizeof(ha_cgetype));
  ha_cgetype *x0;//,*vecbivi0;//= (ha_cgetype*)calloc(VecSize,sizeof(ha_cgetype));
  long int lnz;
  if(rank==mpisize-1) {
    Mat_SeqAIJ         *aa=(Mat_SeqAIJ*)submatD[0]->data;//*aa=subA->data;
    ai= aa->i;
    aj= aa->j;
    vals=aa->a;
    lnz=aa->nz;
    nrow=submatD[0]->rmap->n;
    ncol=submatD[0]->cmap->n;
    printf("D ncol %d nrow %d nz %ld\n",ncol,nrow,lnz);
    obiviindx0=realloc(obiviindx0,(lnz+1)*sizeof(long int));
    obiviindx0[0]=-1;
    biviindx0=obiviindx0;
    biviindx0++;
    lj2=0;
    lBBrow=ncol;
    for(li=0; li<nrow-1; li++) {
      j7l=li*lBBrow;
      for(j=ai[li]; j<ai[li+1]; j++) {
        biviindx0[lj2]=aj[j]+j7l;
        lj2++;
      }
    }
    j7l=li*lBBrow;
    for(j=ai[li]; j<lnz; j++) {
      biviindx0[lj2]=aj[j]+j7l;
      lj2++;
    }
    my_spar_compl_(biviindx1,&nz0,biviindx0,&lnz,&nz1);
    printf("nz0 %ld nz %ld nz1 %ld\n",nz0,lnz,nz1);
    ldsize=ceil((laD/100.0)*nz1);
    vecbivi=realloc(vecbivi,ldsize*sizeof(ha_cgetype));
    int *irn1=(int *) calloc (nz1,sizeof(int));
    int *jcn=(int *) calloc (ldsize,sizeof(int));
    printf("nz1 %ld obin %ld ond1 %ld\n",nz1,obiviindx1[0],obiviindx0[0]);
    my_spar_add4l_(vecbivi,obiviindx1,irn1,jcn,&nz0,vals,obiviindx0,&lnz,&nz1,&ncol);
    printf("nz1 %ld\n",nz1);
    MatDestroy(&submatD[0]);
    lnz=nz1;
    free(obiviindx1);
    obiviindx1=NULL;
    free(obiviindx0);
    obiviindx0=NULL;
    irn1=realloc(irn1,ldsize*sizeof(int));
    int *insizeD=(int *) calloc (6,sizeof(int));
    insizeD[0]=nrow;
    insizeD[1]=ncol;
    insizeD[2]=lnz;
    insizeD[3]=laD;
    insizeD[5]=ldsize;
    xd=(ha_cgetype *) calloc (vecbiuisize,sizeof(ha_cgetype));//realloc (xd,vecbiuisize*sizeof(ha_cgetype));
    spec48m_ssol2la_(insizeD,irn1,jcn,vecbivi,vecbiui,xd);
    free(vecbivi);
    vecbivi=NULL;
    free(vecbiui);
    vecbiui=NULL;
    free(irn1);
    free(jcn);
    free(insizeD);
    x0=(ha_cgetype *) calloc (VecSize,sizeof(ha_cgetype));//realloc (x0,*sizeof(ha_cgetype));
    for(i=0; i<vecbiuisize; i++) {
      x0[indicesC[i]]=xd[i];
    }
  }
  free(vecbiui);
  obiviindx1=realloc(obiviindx1,sizeof(long int));
  free(obiviindx1);
  obiviindx0=realloc(obiviindx0,sizeof(long int));
  free(obiviindx0);
  free(vecbivi);
  free(indicesC);
  if(rank!=mpisize-1) {
    xd=(ha_cgetype *) calloc (vecbiuisize,sizeof(ha_cgetype));//realloc (xd,vecbiuisize*sizeof(ha_cgetype));
    x0=(ha_cgetype *) calloc (VecSize,sizeof(ha_cgetype));//realloc (x0,*sizeof(ha_cgetype));
  }
  PetscFree(submatD);
  MPI_Barrier(PETSC_COMM_WORLD);
  if(SORD==1)MPI_Bcast(xd, vecbiuisize, MPI_DOUBLE,mpisize-1, PETSC_COMM_WORLD);
  else MPI_Bcast(xd, vecbiuisize, MPI_FLOAT,mpisize-1, PETSC_COMM_WORLD);
  printf("Interface problem solution Rank %d time %f\n",rank,((double)(clock()-timestr))/CLOCKS_PER_SEC);
  timestr=clock();
  xi1indx=0;

  omp_set_num_threads(smallthreads);
  #pragma omp parallel private(jthrd,j1,j2,j,i,timestr,ai,aj,vals,nz,nrow,xi1point,xi1indx,ierr,nindices) shared(submatC,submatCij,submatBij,insize,xi1,colindicesbc1,colindicesbc2,x0)
  {
  jthrd=omp_get_thread_num();
  ha_cgetype *biui0= (ha_cgetype *) calloc (maxrowc,sizeof(ha_cgetype));
  ha_cgetype *be0 = (ha_cgetype*)calloc(maxrowc,sizeof(ha_cgetype));
  #pragma omp for
  for(j1=0; j1<nmatint; j1++) {
    xi1indx=vecnrowc[j1];
    timestr=clock();
    j2=j1*(nreg+1);
    Mat_SeqAIJ         *ac=(Mat_SeqAIJ*)submatC[j1]->data;//*aa=subA->data;
    ai= ac->i;
    aj= ac->j;
    vals=ac->a;
    nz=ac->nz;
    nrow=submatC[j1]->rmap->n;
    spar_mulnoadd_(xd,&nrow,&nz,ai,aj,vals,be0);
    MatDestroy(&submatC[j1]);
    ifremove=true;
    NDBBD_sol(rank,j2,nreg,insize,insizes,submatCij,submatBij,be0,biui0,ifremove,fn01,fn02,fn03);
    for(i=0; i<nreg; i++) {
      MatDestroy(&submatCij[j2+i][0]);
      MatDestroy(&submatBij[j2+i][0]);
    }
    xi1point=xi1+xi1indx;
    ISGetIndices(colindicesbc1[j1],&nindices);
    ISGetSize(colindicesbc1[j1],&j);
    for(i=0; i<j; i++) {
      x0[nindices[i]]+=xi1point[i]-biui0[i];//vals[i]-biui0[i];//vecval;
    }
    ISRestoreIndices(colindicesbc1[j1],&nindices);
    ierr = ISDestroy(&colindicesbc1[j1]);
    ISGetIndices(colindicesbc2[j1],&nindices);
    for(j=i; j<nrow; j++) {
      x0[nindices[j-i]]+=xi1point[j]-biui0[j];//vals[i]-biui0[i];//vecval;
    }
    ISRestoreIndices(colindicesbc2[j1],&nindices);
    ierr = ISDestroy(&colindicesbc2[j1]);

    timeend=clock();//time(&timeend);
    printf("Submatrix %d rank %d calculation time %f\n",j1,rank,difftime(timeend,timestr)/CLOCKS_PER_SEC);
  }
  free(be0);
  free(biui0);
  }
  omp_set_num_threads(mymaxnumthrd);
  free(vecnrowc);
  for(i=0; i<nmatin; i++) {
    ierr = PetscFree(submatCij[i]);
    ierr = PetscFree(submatBij[i]);
  }
  ierr = PetscFree(colindicesbc1);
  ierr = PetscFree(colindicesbc2);
  free(submatCij);
  free(submatBij);
  for(j1=0; j1<nmatinplus; j1++) {
    free(fn01[j1]);
    free(fn02[j1]);
    free(fn03[j1]);
  }
  free(fn01);
  free(fn02);
  free(fn03);

  free(xi1);
  free(insize);
  ierr = PetscFree(submatC);
  printf("Solution calculation Rank %d time %f\n",rank,((double)(clock()-timestr))/CLOCKS_PER_SEC);
  MPI_Barrier(PETSC_COMM_WORLD);
  free(xd);
  timestr=clock();
  if(SORD==1)MPI_Allreduce(x0,x1,VecSize, MPI_DOUBLE, MPI_SUM,PETSC_COMM_WORLD);
  else MPI_Allreduce(x0,x1,VecSize, MPI_FLOAT, MPI_SUM,PETSC_COMM_WORLD);
  printf("Reduce solution Rank %d time %f\n",rank,((double)(clock()-timestr))/CLOCKS_PER_SEC);
  free(x0);
  printf("OK rui rank %d\n",rank);
  return 0;
}




bool NDBBD_sol(PetscInt rank, int begmat,int nreg,int * insize,int insizes, Mat **submatCij,Mat **submatBij,ha_cgetype *b,ha_cgetype *sol,bool ifremove,char** fn01,char** fn02, char** fn03) {
  FILE* fp1,*fp2,*fp3;
  PetscScalar *vals,*valsc,vecval;
  PetscInt *ai,*aj,*aic,*ajc,nrow,nz,maxrowcij;//,ncol
  PetscInt i,j,j1,j2,indx01,la1;
  ha_cgetype *b01,*b03,*sol1,*sol2;
  size_t freadresult;
  int *irne= NULL;//(int*)calloc(1,sizeof(int));
  int *keep= NULL;//(int*)calloc(1,sizeof(int));
  ha_cgetype *vale = NULL;//(ha_cgetype*)calloc(1,sizeof(ha_cgetype));
  maxrowcij=0;
  b01=b;
  sol1=sol;
  int fd1,fd2,fd3,frrsl1,frrsl2,frrsl3;
  for (j1=begmat; j1<nreg+begmat; j1++) {
    j2=j1*insizes;//-begmat;
    if((submatCij[j1][0]->rmap->n)>maxrowcij)maxrowcij=submatCij[j1][0]->rmap->n;
    #pragma omp critical(nsol)
    {
    if(isLinux==1) {
      la1=ceil((insize[j2+9]/100.0)*insize[j2+2]);
      irne = realloc(irne,la1*sizeof(int));//(int*)calloc(la1,sizeof(int));
      keep = realloc(keep,insize[j2+12]*sizeof(int));//(int*)calloc(insize[j1*insizes+12],sizeof(int));
      vale = realloc(vale,la1*sizeof(ha_cgetype));//(ha_cgetype*)calloc(la1,sizeof(ha_cgetype));
      fd1 = open(fn01[j1],O_RDONLY);
      frrsl1=read(fd1,irne,sizeof(int)*la1);
      close(fd1);
      fd2 = open(fn02[j1],O_RDONLY);
      frrsl2=read(fd2,keep,insize[j2+12]*sizeof(int));
      close(fd2);
      fd3 = open(fn03[j1],O_RDONLY);
      frrsl3=read(fd3,vale,la1*sizeof(ha_cgetype));
      close(fd3);
    }
    else {
      la1=ceil((insize[j2+9]/100.0)*insize[j2+2]);
      irne = realloc(irne,la1*sizeof(int));//(int*)calloc(la1,sizeof(int));
      keep = realloc(keep,insize[j2+12]*sizeof(int));//(int*)calloc(insize[j1*insizes+12],sizeof(int));
      vale = realloc(vale,la1*sizeof(ha_cgetype));//(ha_cgetype*)calloc(la1,sizeof(ha_cgetype));
      fp1 = fopen(fn01[j1], "rb");
      if (fp1==NULL)printf("File opening error\n");
      freadresult=fread(irne,sizeof(int),la1,fp1);
      fclose(fp1);
      fp2 = fopen(fn02[j1], "rb");
      if (fp2==NULL)printf("File opening error\n");
      freadresult=fread(keep,sizeof(int),insize[j2+12],fp2);
      fclose(fp2);
      fp3 = fopen(fn03[j1], "rb");
      if (fp3==NULL)printf("File opening error\n");
      freadresult=fread(vale,sizeof(ha_cgetype),la1,fp3);
      fclose(fp3);
    }
    }
    if(insize[j2+16]!=la1)insize[j2+16]=la1;
    spec48m_esol_(insize+j2,irne,vale,keep,b01,sol1);
    b01+=insize[j2];
    sol1+=insize[j2];
  }
  sol2=sol;
  for (j1=begmat; j1<nreg+begmat; j1++) {
    Mat_SeqAIJ         *ab=(Mat_SeqAIJ*)submatBij[j1][0]->data;
    ai= ab->i;
    aj= ab->j;
    vals=ab->a;
    nz=ab->nz;
    nrow=submatBij[j1][0]->rmap->n;
    spar_mulmin_(sol2,&nrow,&nz,ai,aj,vals,b01);
    sol2+=insize[j1*insizes];
  }
  j2=j1;//-begmat;
    #pragma omp critical(nsol)
    {
  if(isLinux==1) {
    la1=ceil((insize[j1*insizes+9]/100.0)*insize[j1*insizes+2]);
    irne = realloc(irne,la1*sizeof(int));//(int*)calloc(la1,sizeof(int));
    keep = realloc(keep,insize[j1*insizes+12]*sizeof(int));//(int*)calloc(insize[j1*insizes+12],sizeof(int));
    vale = realloc(vale,la1*sizeof(ha_cgetype));//(ha_cgetype*)calloc(la1,sizeof(ha_cgetype));
    fd1 = open(fn01[j2],O_RDONLY);
    frrsl1=read(fd1,irne,sizeof(int)*la1);
    close(fd1);
    fd2 = open(fn02[j2],O_RDONLY);
    frrsl2=read(fd2,keep,insize[j1*insizes+12]*sizeof(int));
    close(fd2);
    fd3 = open(fn03[j2],O_RDONLY);
    frrsl3=read(fd3,vale,la1*sizeof(ha_cgetype));
    close(fd3);
  }
  else {
    la1=ceil((insize[j1*insizes+9]/100.0)*insize[j1*insizes+2]);
    irne = realloc(irne,la1*sizeof(int));//(int*)calloc(la1,sizeof(int));
    keep = realloc(keep,insize[j1*insizes+12]*sizeof(int));//(int*)calloc(insize[j1*insizes+12],sizeof(int));
    vale = realloc(vale,la1*sizeof(ha_cgetype));//(ha_cgetype*)calloc(la1,sizeof(ha_cgetype));
    fp1 = fopen(fn01[j2], "rb");
    if (fp1==NULL)printf("File opening error\n");
    freadresult=fread(irne,sizeof(int),la1,fp1);
    fclose(fp1);
    fp2 = fopen(fn02[j2], "rb");
    if (fp2==NULL)printf("File opening error\n");
    freadresult=fread(keep,sizeof(int),insize[j1*insizes+12],fp2);
    fclose(fp2);
    fp3 = fopen(fn03[j2], "rb");
    if (fp3==NULL)printf("File opening error\n");
    freadresult=fread(vale,sizeof(ha_cgetype),la1,fp3);
    fclose(fp3);
  }
    }
  if(ifremove)remove(fn01[j2]);
  if(ifremove)remove(fn02[j2]);
  if(ifremove)remove(fn03[j2]);
  if(insize[j1*insizes+16]!=la1)insize[j1*insizes+16]=la1;
  spec48m_esol_(insize+j1*insizes,irne,vale,keep,b01,sol1);
  sol2=sol;
  b03=b;
  ha_cgetype *b02 = (ha_cgetype*)calloc(maxrowcij,sizeof(ha_cgetype));
  for (j1=begmat; j1<nreg+begmat; j1++) {
    j2=j1;//-begmat;
    #pragma omp critical(nsol)
    {
    if(isLinux==1) {
      la1=ceil((insize[j1*insizes+9]/100.0)*insize[j1*insizes+2]);
      irne = realloc(irne,la1*sizeof(int));//(int*)calloc(la1,sizeof(int));
      keep = realloc(keep,insize[j1*insizes+12]*sizeof(int));//(int*)calloc(insize[j1*insizes+12],sizeof(int));
      vale = realloc(vale,la1*sizeof(ha_cgetype));//(ha_cgetype*)calloc(la1,sizeof(ha_cgetype));
      fd1 = open(fn01[j2],O_RDONLY);
      frrsl1=read(fd1,irne,sizeof(int)*la1);
      close(fd1);
      fd2 = open(fn02[j2],O_RDONLY);
      frrsl2=read(fd2,keep,insize[j1*insizes+12]*sizeof(int));
      close(fd2);
      fd3 = open(fn03[j2],O_RDONLY);
      frrsl3=read(fd3,vale,la1*sizeof(ha_cgetype));
      close(fd3);
    }
    else {
      la1=ceil((insize[j1*insizes+9]/100.0)*insize[j1*insizes+2]);
      irne = realloc(irne,la1*sizeof(int));//(int*)calloc(la1,sizeof(int));
      keep = realloc(keep,insize[j1*insizes+12]*sizeof(int));//(int*)calloc(insize[j1*insizes+12],sizeof(int));
      vale = realloc(vale,la1*sizeof(ha_cgetype));//(ha_cgetype*)calloc(la1,sizeof(ha_cgetype));
      fp1 = fopen(fn01[j2], "rb");
      if (fp1==NULL)printf("File opening error\n");
      freadresult=fread(irne,sizeof(int),la1,fp1);
      fclose(fp1);
      fp2 = fopen(fn02[j2], "rb");
      if (fp2==NULL)printf("File opening error\n");
      freadresult=fread(keep,sizeof(int),insize[j1*insizes+12],fp2);
      fclose(fp2);
      fp3 = fopen(fn03[j2], "rb");
      if (fp3==NULL)printf("File opening error\n");
      freadresult=fread(vale,sizeof(ha_cgetype),la1,fp3);
      fclose(fp3);
    }
    }
    if(ifremove)remove(fn01[j2]);
    if(ifremove)remove(fn02[j2]);
    if(ifremove)remove(fn03[j2]);
    Mat_SeqAIJ         *ac=(Mat_SeqAIJ*)submatCij[j1][0]->data;//*aa=subA->data;
    ai= ac->i;
    aj= ac->j;
    vals=ac->a;
    nz=ac->nz;
    nrow=submatCij[j1][0]->rmap->n;
    spar_mulnoadd_(sol1,&nrow,&nz,ai,aj,vals,b02);
    if(insize[j1*insizes+16]!=la1)insize[j1*insizes+16]=la1;
    spec48m_esol_(insize+j1*insizes,irne,vale,keep,b02,b03);
    for(j=0; j<nrow; j++)sol2[j]-=b03[j];
    sol2+=nrow;
  }
  free(b02);
  free(irne);
  free(keep);
  free(vale);
  return true;
}


bool NDBBD_sol_nread1(PetscInt rank, int begmat,int nreg,int * insize,int insizes, Mat **submatCij,Mat **submatBij,ha_cgetype *b,ha_cgetype *sol,int** irnereg,int** keepreg,ha_cgetype** valereg,ha_cgetype *cntl,ha_cgetype *rinfo,ha_cgetype *error1,int *icntl,int *info,ha_cgetype *w,int *iw,ha_cgetype *b02) {
  PetscScalar *vals;//,*valsc;//,vecval;
  PetscInt *ai,*aj,nrow,nz;//,ncol,*aic,*ajc
  PetscInt i,j,j1,j2;//,indx01,la1;
  ha_cgetype *b01,*b03,*sol1,*sol2;//,*b02
  int * insize1;
  b01=b;
  sol1=sol;
  for (j1=begmat; j1<nreg+begmat; j1++) {
    j2=j1-begmat;//-begmat;
    j=j1*insizes;
    insize1=insize+j;
    insize1[16]=ceil((insize1[9]/100.0)*insize1[2]);
    spec48m_rpesol_(insize1,irnereg[j2],valereg[j2],keepreg[j2],b01,sol1,cntl,rinfo,error1,icntl,info,w,iw);//insize+j
    b01+=insize[j];
    sol1+=insize[j];
  }
  sol2=sol;
  for (j1=begmat; j1<nreg+begmat; j1++) {
    Mat_SeqAIJ         *ab=(Mat_SeqAIJ*)submatBij[j1][0]->data;
    ai= ab->i;
    aj= ab->j;
    vals=ab->a;
    nz=ab->nz;
    nrow=submatBij[j1][0]->rmap->n;
    spar_mulmin_(sol2,&nrow,&nz,ai,aj,vals,b01);
    sol2+=insize[j1*insizes];
  }
  j2=j1-begmat;//-begmat;
  insize1=insize+j1*insizes;
  insize1[16]=ceil((insize1[9]/100.0)*insize1[2]);
  spec48m_rpesol_(insize1,irnereg[j2],valereg[j2],keepreg[j2],b01,sol1,cntl,rinfo,error1,icntl,info,w,iw);//insize+j1*insizes
  sol2=sol;
  b03=b;
  for (j1=begmat; j1<nreg+begmat; j1++) {
    j2=j1-begmat;//-begmat;
    Mat_SeqAIJ         *ac=(Mat_SeqAIJ*)submatCij[j1][0]->data;//*aa=subA->data;
    ai= ac->i;
    aj= ac->j;
    vals=ac->a;
    nz=ac->nz;
    nrow=submatCij[j1][0]->rmap->n;
    spar_mulnoadd_(sol1,&nrow,&nz,ai,aj,vals,b02);
    insize1=insize+j1*insizes;
    insize1[16]=ceil((insize1[9]/100.0)*insize1[2]);
    spec48m_rpesol_(insize1,irnereg[j2],valereg[j2],keepreg[j2],b02,b03,cntl,rinfo,error1,icntl,info,w,iw);//insize+j1*insizes
    for(j=0; j<nrow; j++)sol2[j]-=b03[j];
    sol2+=nrow;
  }
  return true;
}

int HaReduce(ha_cgetype *vecbivi,forint vecbivisize,PetscInt mpisize,PetscInt rank,PetscInt targetrank) {
  if(mpisize==1)return 0;
  MPI_Status   status;
  int i,j,j1,j2,j3;
  j=0;
  if(rank!=targetrank){
  #pragma omp parallel private(i) reduction(+:j)
  {
  j=0;
  #pragma omp for
    for(i=0; i<vecbivisize; i++)if(vecbivi[i]!=0)j++;
  }
  }
  if(rank==targetrank) j=0;//(VecSize-sumrowcolin)*(VecSize-sumrowcolin);
  int *biviindx= (int *) calloc (j,sizeof(int));
  ha_cgetype *vecbivi0= (ha_cgetype *) calloc (1,sizeof(ha_cgetype));
  j2=0;
  if(rank!=targetrank)for(i=0; i<vecbivisize; i++) {
      if(vecbivi[i]!=0) {
        vecbivi[j2]=vecbivi[i];
        biviindx[j2]=i;
        j2++;
      }
    }
  MPI_Reduce(&j,&j2,1, MPI_INT, MPI_MAX,targetrank,PETSC_COMM_WORLD);
  if(rank==targetrank) {
    vecbivi0=realloc(vecbivi0,j2*sizeof(ha_cgetype));
    biviindx=realloc(biviindx,j2*sizeof(int));
  }

  for(j1=0; j1<targetrank; j1++) {
    j2=j;
    MPI_Bcast(&j2,1, MPI_INT,j1, PETSC_COMM_WORLD);
    j3=j1;
    MPI_Bcast(&j3,1, MPI_INT,j1, PETSC_COMM_WORLD);
    if(rank==j1) {
      if(SORD==1)MPI_Send(vecbivi,j, MPI_DOUBLE, targetrank, 10, PETSC_COMM_WORLD);
      else MPI_Send(vecbivi,j, MPI_FLOAT, targetrank, 10, PETSC_COMM_WORLD);
      MPI_Send(biviindx,j, MPI_INT, targetrank, 11, PETSC_COMM_WORLD);
    }
    if(rank==targetrank) {
      if(SORD==1)MPI_Recv(vecbivi0,j2, MPI_DOUBLE, j3, 10, PETSC_COMM_WORLD,&status);
      else MPI_Recv(vecbivi0,j2, MPI_FLOAT, j3, 10, PETSC_COMM_WORLD,&status);
      MPI_Recv(biviindx,j2, MPI_INT, j3, 11, PETSC_COMM_WORLD,&status);
        #pragma omp parallel private(i)
        {
        #pragma omp for
      for(i=0; i<j2; i++) {
        #pragma omp atomic
        vecbivi[biviindx[i]]+=vecbivi0[i];
      }
        }
    }
    MPI_Barrier(PETSC_COMM_WORLD);
  }
  for(j1=targetrank+1; j1<mpisize; j1++) {
    j2=j;
    MPI_Bcast(&j2,1, MPI_INT,j1, PETSC_COMM_WORLD);
    j3=j1;
    MPI_Bcast(&j3,1, MPI_INT,j1, PETSC_COMM_WORLD);
    if(rank==j1) {
      if(SORD==1)MPI_Send(vecbivi,j, MPI_DOUBLE, targetrank, 10, PETSC_COMM_WORLD);
      else MPI_Send(vecbivi,j, MPI_FLOAT, targetrank, 10, PETSC_COMM_WORLD);
      MPI_Send(biviindx,j, MPI_INT, targetrank, 11, PETSC_COMM_WORLD);
    }
    if(rank==targetrank) {
      if(SORD==1)MPI_Recv(vecbivi0,j2, MPI_DOUBLE, j3, 10, PETSC_COMM_WORLD,&status);
      else MPI_Recv(vecbivi0,j2, MPI_FLOAT, j3, 10, PETSC_COMM_WORLD,&status);
      MPI_Recv(biviindx,j2, MPI_INT, j3, 11, PETSC_COMM_WORLD,&status);
        #pragma omp parallel private(i)
        {
        #pragma omp for
      for(i=0; i<j2; i++) {
        #pragma omp atomic
        vecbivi[biviindx[i]]+=vecbivi0[i];
      }
        }
    }
    MPI_Barrier(PETSC_COMM_WORLD);
  }
  free(biviindx);
  free(vecbivi0);
  return 1;
}

int HaReduceNoComp(ha_cgetype *vecbivi,forint vecbivisize,PetscInt mpisize,PetscInt rank,PetscInt targetrank) {
  if(mpisize==1)return 0;
  MPI_Status   status;
  int i,j,j1,j2,j3;
  ha_cgetype *vecbivi0= (ha_cgetype *) calloc (1,sizeof(ha_cgetype));

  if(rank==targetrank) {
    vecbivi0=realloc(vecbivi0,vecbivisize*sizeof(ha_cgetype));
  }

  for(j1=0; j1<targetrank; j1++) {
    if(rank==j1) {
      if(SORD==1)MPI_Send(vecbivi,vecbivisize, MPI_DOUBLE, targetrank, 101, PETSC_COMM_WORLD);
      else MPI_Send(vecbivi,vecbivisize, MPI_FLOAT, targetrank, 101, PETSC_COMM_WORLD);
    }
    if(rank==targetrank) {
      if(SORD==1)MPI_Recv(vecbivi0,vecbivisize, MPI_DOUBLE, j1, 101, PETSC_COMM_WORLD,&status);
      else MPI_Recv(vecbivi0,vecbivisize, MPI_FLOAT, j1, 101, PETSC_COMM_WORLD,&status);
        #pragma omp parallel private(i)
        {
        #pragma omp for
      for(i=0; i<vecbivisize; i++) {
        #pragma omp atomic
        vecbivi[i]+=vecbivi0[i];
      }
        }
    }
    MPI_Barrier(PETSC_COMM_WORLD);
  }
  for(j1=targetrank+1; j1<mpisize; j1++) {
    if(rank==j1) {
      if(SORD==1)MPI_Send(vecbivi,vecbivisize, MPI_DOUBLE, targetrank, 102, PETSC_COMM_WORLD);
      else MPI_Send(vecbivi,vecbivisize, MPI_FLOAT, targetrank, 102, PETSC_COMM_WORLD);
    }
    if(rank==targetrank) {
      if(SORD==1)MPI_Recv(vecbivi0,vecbivisize, MPI_DOUBLE, j1, 102, PETSC_COMM_WORLD,&status);
      else MPI_Recv(vecbivi0,vecbivisize, MPI_FLOAT, j1, 102, PETSC_COMM_WORLD,&status);
        #pragma omp parallel private(i)
        {
        #pragma omp for
      for(i=0; i<vecbivisize; i++) {
        #pragma omp atomic
        vecbivi[i]+=vecbivi0[i];
      }
        }
    }
    MPI_Barrier(PETSC_COMM_WORLD);
  }
  free(vecbivi0);
  return 1;
}




bool spline(ha_cgetype* y,ha_cgetype* x,ha_cgetype sx0,ha_cgetype sxn,int size,ha_cgetype* w,int laA){//sizeof(x)=size+1
  uvadd lasize;
  ha_cgetype *matval= (ha_cgetype *) calloc (6*4*size,sizeof(ha_cgetype));
  ha_cgetype *bval= (ha_cgetype *) calloc (4*size,sizeof(ha_cgetype));
  int *irn= (int *) calloc (6*4*size,sizeof(int));
  int *jcn= (int *) calloc (6*4*size,sizeof(int));
  if(size<3){
    printf("Too few points!!!\n");
    return 0;
  }
  int i,j,mindx=0,bindx=0,i1,i2;
  int *insize=(int *) calloc (4,sizeof(int));
  insize[0]=4*size;
  insize[1]=4*size;
  
  for(i=1;i<size+1;i++){//ai+bixi+cixi2+dixi3=yi//indx 1:n
    matval[mindx]=1;
    irn[mindx]=i;
    jcn[mindx]=i;
    mindx++;
    matval[mindx]=x[i];
    irn[mindx]=i;
    jcn[mindx]=i+size;
    mindx++;
    matval[mindx]=x[i]*x[i];
    irn[mindx]=i;
    jcn[mindx]=i+2*size;
    mindx++;
    matval[mindx]=x[i]*x[i]*x[i];
    irn[mindx]=i;
    jcn[mindx]=i+3*size;
    mindx++;
    bval[bindx]=y[i];
    bindx++;
  }
  for(i=size;i<2*size;i++){//aip+bipxi+cipxi2+dipxi3=yi//indx 0:n-1
    i1=i+1;
    i2=i-size;
    matval[mindx]=1;
    irn[mindx]=i1;
    jcn[mindx]=i1-size;
    mindx++;
    matval[mindx]=x[i2];
    irn[mindx]=i1;
    jcn[mindx]=i1;
    mindx++;
    matval[mindx]=x[i2]*x[i2];
    irn[mindx]=i1;
    jcn[mindx]=i1+size;
    mindx++;
    matval[mindx]=x[i2]*x[i2]*x[i2];
    irn[mindx]=i1;
    jcn[mindx]=i1+2*size;
    mindx++;
    bval[bindx]=y[i2];
    bindx++;
  }
    matval[mindx]=1;
    irn[mindx]=2*size+1;
    jcn[mindx]=1+size;
    mindx++;
    matval[mindx]=2*x[0];
    irn[mindx]=2*size+1;
    jcn[mindx]=1+2*size;
    mindx++;
    matval[mindx]=3*x[0]*x[0];
    irn[mindx]=2*size+1;
    jcn[mindx]=1+3*size;
    mindx++;
    bval[bindx]=(y[1]-y[0])/(x[1]-x[0]);//sx0;//
    bindx++;
  for(i=2*size+2;i<3*size+1;i++){//indx 1:n-1
    matval[mindx]=1;
    irn[mindx]=i;
    jcn[mindx]=i-size-1;
    mindx++;
    matval[mindx]=-1;
    irn[mindx]=i;
    jcn[mindx]=i-size;
    mindx++;

    matval[mindx]=2*x[i-2*size-1];
    irn[mindx]=i;
    jcn[mindx]=i-1;
    mindx++;
    matval[mindx]=-2*x[i-2*size-1];
    irn[mindx]=i;
    jcn[mindx]=i;
    mindx++;

    matval[mindx]=3*x[i-2*size-1]*x[i-2*size-1];
    irn[mindx]=i;
    jcn[mindx]=i+size-1;
    mindx++;
    matval[mindx]=-3*x[i-2*size-1]*x[i-2*size-1];
    irn[mindx]=i;
    jcn[mindx]=i+size;
    mindx++;
    bval[bindx]=0;
    bindx++;
  }
  for(i=3*size+1;i<4*size;i++){//indx 1:n-1
    matval[mindx]=2;
    irn[mindx]=i;
    jcn[mindx]=i-size;
    mindx++;
    matval[mindx]=-2;
    irn[mindx]=i;
    jcn[mindx]=i-size+1;
    mindx++;

    matval[mindx]=6*x[i-3*size];
    irn[mindx]=i;
    jcn[mindx]=i;
    mindx++;
    matval[mindx]=-6*x[i-3*size];
    irn[mindx]=i;
    jcn[mindx]=i+1;
    mindx++;
    bval[bindx]=0;
    bindx++;
  }
    matval[mindx]=1;
    irn[mindx]=4*size;
    jcn[mindx]=2*size;
    mindx++;
    matval[mindx]=2*x[size];
    irn[mindx]=4*size;
    jcn[mindx]=3*size;
    mindx++;
    matval[mindx]=3*x[size]*x[size];
    irn[mindx]=4*size;
    jcn[mindx]=4*size;
    mindx++;
    bval[bindx]=(y[size]-y[size-1])/(x[size]-x[size-1]);//sxn;//
    printf("slop %lf app slop %lf y1 %lf y0 %lf\n",sxn,bval[bindx],y[size],y[size-1]);
    bindx++;
  insize[2]=mindx;
  lasize=ceil((laA*10/100.0)*mindx);
  insize[3]=laA*10;
  irn=realloc(irn,lasize*sizeof(int));
  jcn=realloc(jcn,lasize*sizeof(int));
  matval=realloc(matval,lasize*sizeof(ha_cgetype));
  spec48_ssol2la_(insize,irn,jcn,matval,bval,w);
    printf("bindx %d val %lf w1 %lf w2 %lf w3 %lf x %lf\n",bindx,w[size-1]+w[size-1+size]*x[size]+w[size-1+2*size]*x[size]*x[size]+w[size-1+3*size]*x[size]*x[size]*x[size],w[size-1],w[size-1+size],w[size-1+2*size],y[size]);
  free(insize);
  free(matval);
  free(bval);
  free(irn);
  free(jcn);
  return 1;
}

bool Johansen(PetscBool nohsl,PetscInt VecSize,Mat A,PetscInt dnz,PetscInt* dnnz,PetscInt onz,PetscInt* onnz,Mat B,PetscInt dnzB,PetscInt* dnnzB,PetscInt onzB,PetscInt* onnzB,Vec vecb,Vec vece,PetscInt rank,PetscInt rank_hsl,PetscInt mpisize,char* tabfile, char *commsyntax,ha_cgeset *ha_set,uvdim nset, ha_cgesetele *ha_setele, hcge_cof *ha_cof,uvadd ncof,hcge_cof *ha_var,uvadd nvar, ha_cgevar **ha_cofvar2,uvadd ncofvar,uvadd ncofele,uvadd nvarele,ha_cgeexovar **ha_cgeshock2,uvadd alltimeset,uvadd allregset,uvadd nintraeq,uvdim matsol,PetscInt Istart,PetscInt Iend,  uvadd nreg, uvadd ntime, uvadd *ha_eqadd, uvadd ndblock, uvadd *countvarintra1, uvadd *counteq, uvadd *counteqnoadd,uvdim laA,uvdim laDi,uvdim laD,PetscReal cntl3,PetscReal cntl6,PetscBool presol,uvdim nesteddbbd,int localsize,PetscInt *ndbbddrank1,forint* indata,uvdim mc66,forint *ptx,struct timeval begintime,ha_cgetype **xcf2){ //Johansen
  char tempfilenam[256],tempchar[256];
  PetscScalar value,*vals=NULL;
  PetscErrorCode ierr;
  PetscInt count,nz01,*ai=NULL,*aj=NULL;
  forint k=0,m=1;
  uvadd i,j,lasize;
  ha_cgetype *b1=NULL,*x0=NULL;
  bool IsIni;
  FILE* tempvar;
  ha_cgeexovar *ha_cgeshock;
  ha_cgeshock=*ha_cgeshock2;
  ha_cgetype *xcf;
  xcf=*xcf2;

  clock_t timestr,timeend,timemulti;
  struct timeval endtime;
  struct timespec gettime_now,gettime_beg,gettime_end;
  long int start_time=0;
  double rep_time;
  size_t freadresult;
  ha_cgevar *ha_cofvar;
  ha_cofvar=*ha_cofvar2;
  ha_cgevar *ha_cofvar1=NULL;
  
    if(nohsl) {
      MatCreate(PETSC_COMM_WORLD,&A);
    }
    else {
      MatCreate(PETSC_COMM_SELF,&A);
    }
    if(nesteddbbd==1)MatSetSizes(A,localsize,localsize,VecSize,VecSize);
    else MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,VecSize,VecSize);
    if(nohsl) {
      MatSetType(A,MATMPIAIJ);
    }
    else {
      MatSetType(A,MATSEQAIJ);
    }
    printf("OK1 nohsl %d!!!\n",nohsl);
    if(nohsl) {
      MatMPIAIJSetPreallocation(A,dnz,dnnz,onz,onnz);
    }
    else {
      MatSeqAIJSetPreallocation(A,dnz,dnnz);
    }

    ierr = MatSetOption(A,MAT_SYMMETRIC,PETSC_FALSE);
    CHKERRQ(ierr);


    if(nohsl) {
      MatCreate(PETSC_COMM_WORLD,&B);
    }
    else {
      MatCreate(PETSC_COMM_SELF,&B);
    }
    if(nesteddbbd==1)MatSetSizes(B,localsize,localsize,VecSize,VecSize);
    else MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,VecSize,VecSize);
    if(nohsl) {
      MatSetType(B,MATMPIAIJ);
    }
    else {
      MatSetType(B,MATSEQAIJ);
    }
    if(nohsl) {
      MatMPIAIJSetPreallocation(B,dnzB,dnnzB,onzB,onnzB);
    }
    else {
      MatSeqAIJSetPreallocation(B,dnzB,dnnzB);
    }

    gettimeofday(&endtime, NULL);
    if(rank==0)printf("Matrix preparation time %f\n",(endtime.tv_sec - begintime.tv_sec)+((double)(endtime.tv_usec - begintime.tv_usec))/ 1000000);
    
    if(rank==rank_hsl) {
      HaNewMatVal(tabfile,commsyntax,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele,ha_cgeshock,ndblock,alltimeset,allregset,ha_eqadd,counteq,nintraeq,A,B);
    }

    gettimeofday(&begintime, NULL);
    if(rank==0)printf("Matrix calculation time %f\n",(begintime.tv_sec - endtime.tv_sec)+((double)(begintime.tv_usec - endtime.tv_usec))/ 1000000);

    for (count=0; count<nvarele; count++) {
      if (ha_cgeshock[count].ShockId) {
        value = ha_cgeshock[count].ShockVal;
        dnz=ha_cgeshock[count].ExoIndx;
        VecSetValues(vece,1,&dnz,&value,INSERT_VALUES);
      }
    }
    MPI_Barrier(PETSC_COMM_WORLD);
    ierr = VecAssemblyBegin(vece);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(vece);
    CHKERRQ(ierr);
    if(rank==rank_hsl) {
      strcpy(tempfilenam,temdir);
      strcat(tempfilenam,"_tempshock");
      sprintf(tempchar, "%d",rank);
      strcat(tempfilenam,tempchar);
      strcat(tempfilenam,".bin");
      if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
        printf("Error opening file\n");
      }
      fwrite(ha_cgeshock, sizeof(ha_cgeexovar),nvarele, tempvar);
      fclose(tempvar);
      free(*ha_cgeshock2);//
      *ha_cgeshock2=NULL;//realloc (ha_cgeshock,1*sizeof(ha_cgeexovar));
      ha_cgeshock=*ha_cgeshock2;
    }
    if(rank==rank_hsl) {
      strcpy(tempfilenam,temdir);
      strcat(tempfilenam,"_tempvar");
      sprintf(tempchar, "%d",rank);
      strcat(tempfilenam,tempchar);
      strcat(tempfilenam,".bin");
      if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
        printf("Error opening file\n");
      }
      fwrite(ha_cofvar, sizeof(ha_cgevar),ncofele+nvarele, tempvar);
      fclose(tempvar);
      free(*ha_cofvar2);//
      *ha_cofvar2=NULL;//realloc (ha_cofvar,1*sizeof(ha_cgevar));
      ha_cofvar=*ha_cofvar2;
    }

    MPI_Barrier(PETSC_COMM_WORLD);
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    printf("OK11\n");
    ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    gettimeofday(&endtime, NULL);
    if(rank==0)printf("Matrix assembly time %f\n",(endtime.tv_sec - begintime.tv_sec)+((double)(endtime.tv_usec - begintime.tv_usec))/ 1000000);
    CHKERRQ(ierr);
    PetscViewer viewer;
    ierr = VecDuplicate(vece,&vecb);
    CHKERRQ(ierr);
    ierr = MatMult(B,vece,vecb);
    CHKERRQ(ierr);
    ierr = VecAssemblyBegin(vecb);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(vecb);
    CHKERRQ(ierr);
    ierr = MatDestroy(&B);
    CHKERRQ(ierr);
    ierr = VecDestroy(&vece);
    CHKERRQ(ierr);

    if(matsol>1) {
      gettimeofday(&begintime, NULL);
      clock_gettime(CLOCK_REALTIME, &gettime_beg);
      int *ha_rows= (int *) calloc (VecSize,sizeof(int));
      int *ha_cols= (int *) calloc (VecSize,sizeof(int));
      int *ha_ndblocks= (int *) calloc (ndblock,sizeof(int));
      if(matsol==2) {
        HaDBBDMatOder(A,VecSize,mpisize,rank,Istart,Iend,nvarele,ha_eqadd,ha_rows,ha_cols,ndblock,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,cntl6);
        x0=realloc (x0,VecSize*sizeof(ha_cgetype));
        HaDBBDParSol(A,vecb,x0,VecSize,mpisize,rank,Istart,Iend,ha_rows,ha_cols,ndblock,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,laD,cntl3);//,iter
      }
      if(matsol==3) {
        presol=1;
        if(presol){
        HaNDBBDMatOderPre(A,VecSize,mpisize,rank,Istart,Iend,nreg,ntime,nvarele,ha_eqadd,ha_rows,ha_cols,ndblock,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,laDi,cntl6,ndbbddrank1,presol);
        HaNDBBDParPre(A,vecb,x0,VecSize,mpisize,rank,Istart,Iend,ha_rows,ha_cols,ndblock,nreg,ntime,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,presol);//,iter
        }
        presol=0;
        HaNDBBDMatOder(A,VecSize,mpisize,rank,Istart,Iend,nreg,ntime,nvarele,ha_eqadd,ha_rows,ha_cols,ndblock,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,laDi,cntl6,ndbbddrank1,presol);
        x0=realloc (x0,VecSize*sizeof(ha_cgetype));
        HaNDBBDParSol(A,vecb,x0,VecSize,mpisize,rank,Istart,Iend,ha_rows,ha_cols,ndblock,nreg,ntime,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,presol);//,iter
      }
      time(&timeend);
      gettimeofday(&endtime, NULL);
      clock_gettime(CLOCK_REALTIME, &gettime_end);
      rep_time = ((double)(gettime_end.tv_nsec-gettime_beg.tv_nsec))/1000000000.0;
      if(rank==0)printf("One step calculation time %f\n",(endtime.tv_sec - begintime.tv_sec)+((double)(endtime.tv_usec - begintime.tv_usec))/ 1000000);
      if(rank==0)printf("One step calculation real time %lf\n",rep_time);
      free(ha_rows);
      free(ha_cols);
      free(ha_ndblocks);
      printf("rank %d\n",rank);
      MPI_Barrier(PETSC_COMM_WORLD);
    }
    else {

      /* This function should be called to be able to use PETSc routines
         from the FORTRAN subroutines needed by this program */
      gettimeofday(&begintime, NULL);
      MPI_Fint fcomm;
      fcomm = MPI_Comm_c2f(PETSC_COMM_WORLD);
      Mat_SeqAIJ         *aa=(Mat_SeqAIJ*)A->data;
      FILE *ofp;

      if(rank==rank_hsl) {
        ai= aa->i;
        aj= aa->j;
        vals=aa->a;
        nz01=aa->nz;
        count=0;
        for(i=0; i<nz01; i++) if(vals[i]!=0) {
            count++;
          }
        printf("count %d nz %d\n",count,nz01);
      }
      indata[0]=count;//.nz

      if(matsol==1) {
        int *irn=(int *) calloc (count,sizeof(int));
        int *irn1=(int *) calloc (nz01,sizeof(int));
        int *jcn=(int *) calloc (count,sizeof(int));
        ha_cgetype *values= (ha_cgetype *) calloc (count,sizeof(ha_cgetype));
        if(rank==rank_hsl) {
          for(i=0; i<VecSize-1; i++)for(j=ai[i]; j<ai[i+1]; j++) {
              irn1[j]=i+1;
            }
          for(j=ai[VecSize-1]; j<nz01; j++) {
            irn1[j]=VecSize;
          }
          j=0;
          for(i=0; i<nz01; i++) if(vals[i]!=0) {
              irn[j]=irn1[i];
              jcn[j]=aj[i]+1;
              values[j]=vals[i];
              j++;
            }
        }
        free(irn1);
        ierr = MatDestroy(&A);
        CHKERRQ(ierr);

        b1=realloc (b1,VecSize*sizeof(ha_cgetype));
        if(rank==rank_hsl) {
          VecGetArray(vecb,&vals);
          for(i=0; i<VecSize; i++) {
            b1[i]=vals[i];
          }
        }
        ierr = VecDestroy(&vecb);
        CHKERRQ(ierr);

        int *neleperrow= (int *) calloc (VecSize,sizeof(int));
        int *ai1= (int *) calloc (VecSize,sizeof(int));
        if(rank==rank_hsl) {
          j=1;
          for(i=1; i<count; i++) {
            if(irn[i]-irn[i-1]>0) {
              neleperrow[k]=j;
              ai1[k]=m;
              j=1;
              m=i+1;
              k++;
            }
            else {
              j++;
            }
          }
          neleperrow[k]=j;
          ai1[k]=ai1[k-1]+neleperrow[k-1];
        }
        x0=realloc (x0,VecSize*sizeof(ha_cgetype));
        if(mc66!=0)spec48_single_(ptx,irn,jcn,b1,values,x0,neleperrow,ai1,&fcomm);
        free(irn);
        if(mc66==0)spec48_nomc66_(ptx,jcn,b1,values,x0,neleperrow,&fcomm,counteq,countvarintra1);
        free(jcn);
        free(values);
        free(neleperrow);
        free(ai1);
        free(b1);//b1=realloc (b1,sizeof(ha_cgetype));
        b1=NULL;
      }
      else {
        x0=realloc (x0,VecSize*sizeof(ha_cgetype));
        lasize=ceil((laA/100.0)*count);
        int *irn=(int *) calloc (lasize,sizeof(int));
        int *irn1=(int *) calloc (nz01,sizeof(int));
        int *jcn=(int *) calloc (lasize,sizeof(int));
        ha_cgetype *values= (ha_cgetype *) calloc (lasize,sizeof(ha_cgetype));

        if(rank==rank_hsl) {
          for(i=0; i<VecSize-1; i++)for(j=ai[i]; j<ai[i+1]; j++) {
              irn1[j]=i+1;
            }
          for(j=ai[VecSize-1]; j<nz01; j++) {
            irn1[j]=VecSize;
          }
          j=0;
          for(i=0; i<nz01; i++) if(vals[i]!=0) {
              irn[j]=irn1[i];
              jcn[j]=aj[i]+1;
              values[j]=vals[i];
              j++;
            }
          VecGetArray(vecb,&vals);
        }
        free(irn1);
        PetscViewerDrawOpen(PETSC_COMM_WORLD,0,"",0,0,500,500,&viewer);
        ierr = MatView(A,viewer);
        PetscViewerDestroy(&viewer);
        ierr = MatDestroy(&A);
        CHKERRQ(ierr);

        int *insize=(int *) calloc (4,sizeof(int));
        insize[0]=VecSize;
        insize[1]=VecSize;
        insize[2]=count;
        insize[3]=laA;
        if(rank==rank_hsl)spec48_ssol2la_(insize,irn,jcn,values,vals,x0);
        free(insize);
        free(irn);
        free(jcn);
        free(values);
        ierr = VecDestroy(&vecb);
        CHKERRQ(ierr);
      }
      gettimeofday(&endtime, NULL);
      if(rank==0)printf("One step calculation time %f\n",(endtime.tv_sec - begintime.tv_sec)+((double)(endtime.tv_usec - begintime.tv_usec))/ 1000000);
    }
    if(rank==rank_hsl) {
      if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
        printf("Error opening file\n");
      }
      *ha_cofvar2=(ha_cgevar*)realloc (*ha_cofvar2,(ncofele+nvarele)*sizeof(ha_cgevar));
      freadresult=fread(*ha_cofvar2, sizeof(ha_cgevar),ncofele+nvarele, tempvar);
      fclose(tempvar);
      remove(tempfilenam);
      ha_cofvar=*ha_cofvar2;

      strcpy(tempfilenam,temdir);
      strcat(tempfilenam,"_tempshock");
      sprintf(tempchar, "%d",rank);
      strcat(tempfilenam,tempchar);
      strcat(tempfilenam,".bin");
      if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
        printf("Error opening file\n");
      }
      *ha_cgeshock2=(ha_cgeexovar*)realloc (*ha_cgeshock2,(nvarele)*sizeof(ha_cgeexovar));
      freadresult=fread(*ha_cgeshock2, sizeof(ha_cgeexovar),nvarele, tempvar);
      fclose(tempvar);
      remove(tempfilenam);
      ha_cgeshock=*ha_cgeshock2;      
    }
    *xcf2=(ha_cgetype*)realloc (*xcf2,nvarele*sizeof(ha_cgetype));
    xcf=*xcf2;
    printf("Hello world1!\n");
    if(rank==rank_hsl) {
      ha_cofvar1=ha_cofvar+ncofele;
      for(i=0; i<nvar; i++) {
        if(ha_var[i].change_real) {
          for(j=ha_var[i].begadd; j<ha_var[i].matsize+ha_var[i].begadd; j++) {
            if(ha_cgeshock[j].ShockId) {
              ha_cofvar1[j].var0=ha_cofvar1[j].varval;
              ha_cofvar1[j].varval+=ha_cgeshock[j].ShockVal;
              xcf[j]=ha_cgeshock[j].ShockVal;//varchange[j]
              ha_cofvar1[j].csolpupd=ha_cgeshock[j].ShockVal;
            }
            else {
              ha_cofvar1[j].var0=ha_cofvar1[j].varval;
              ha_cofvar1[j].varval+=x0[ha_cgeshock[j].ExoIndx];
              xcf[j]=x0[ha_cgeshock[j].ExoIndx];//varchange[j]
              ha_cofvar1[j].csolpupd=x0[ha_cgeshock[j].ExoIndx];
            }
          }
        }
        else {
          for(j=ha_var[i].begadd; j<ha_var[i].matsize+ha_var[i].begadd; j++) {
            if(ha_cgeshock[j].ShockId) {
              ha_cofvar1[j].var0=ha_cofvar1[j].varval;
              ha_cofvar1[j].varval+=ha_cgeshock[j].ShockVal*ha_cofvar1[j].var0/100;
              xcf[j]=ha_cgeshock[j].ShockVal;//varchange[j]
              ha_cofvar1[j].csolpupd=ha_cgeshock[j].ShockVal;
            }
            else {
              ha_cofvar1[j].var0=ha_cofvar1[j].varval;
              xcf[j]=x0[ha_cgeshock[j].ExoIndx];//varchange[j]
              ha_cofvar1[j].varval+=x0[ha_cgeshock[j].ExoIndx]/100*ha_cofvar1[j].varval;
              ha_cofvar1[j].csolpupd=x0[ha_cgeshock[j].ExoIndx];
            }
          }
        }
      }
      hnew_update(tabfile,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele);
      strcpy(commsyntax,"formula");
      IsIni=false;
      hnew_calcff(tabfile,commsyntax,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele,IsIni);

    }
    ha_cofvar1=NULL;
    free(x0);
    return 1;
  }

bool ModMidPoint(PetscBool nohsl,PetscInt VecSize,Mat* A1,PetscInt dnz,PetscInt* dnnz,PetscInt onz,PetscInt* onnz,Mat* B1,PetscInt dnzB,PetscInt* dnnzB,PetscInt onzB,PetscInt* onnzB,Vec* vecb1,Vec *vece1,PetscInt rank,PetscInt rank_hsl,PetscInt mpisize,char* tabfile, char *commsyntax,ha_cgeset *ha_set,uvdim nset, ha_cgesetele *ha_setele, hcge_cof *ha_cof,uvadd ncof,hcge_cof *ha_var,uvadd nvar, ha_cgevar **ha_cofvar2,uvadd ncofvar,uvadd ncofele,uvadd nvarele,ha_cgeexovar **ha_cgeshock2,uvadd alltimeset,uvadd allregset,uvadd nintraeq,uvdim matsol,PetscInt Istart,PetscInt Iend,  uvadd nreg, uvadd ntime, uvadd *ha_eqadd, uvadd ndblock, uvadd *countvarintra1, uvadd *counteq, uvadd *counteqnoadd,uvdim laA,uvdim laDi,uvdim laD,PetscReal cntl3,PetscReal cntl6,PetscBool presol,uvdim nesteddbbd,int localsize,PetscInt *ndbbddrank1,forint* indata,uvdim mc66,forint *ptx,struct timeval begintime,uvdim subints,MPI_Fint fcomm,ha_cgetype **xcf2,int Isbiupd){ //Modified midpoint Pearson 1991
  char tempfilenam[256],tempchar[256],solchar[255];
  PetscScalar value,*vals;
  PetscErrorCode ierr;
  PetscInt count,nz01,*ai,*aj;
  forint k=0,m=1;
  forint tindx1;//,tindx2;
  ha_cgetype temp1,temp2;
  uvadd i,j,lasize;
  uvdim subindx;
  ha_cgetype *b1=NULL;
  ha_cgetype *x1=NULL;
  ha_cgetype *xcf;
  xcf=*xcf2;
  bool IsIni;
  FILE* tempvar;
  PetscLogDouble time1,time0;
  clock_t timestr,timeend,timemulti;
  struct timeval endtime;
  struct timespec gettime_now,gettime_beg,gettime_end;
  long int start_time=0;
  double rep_time;
  size_t freadresult;
  ha_cgevar *ha_cofvar;
  ha_cofvar=*ha_cofvar2;
  ha_cgevar *ha_cofvar1;
  ha_cgeexovar *ha_cgeshock;
  ha_cgeshock=*ha_cgeshock2;
  Vec vece,vecb;
  Mat A,B;
  A=*A1;
  B=*B1;
  vece=*vece1,
  vecb=*vecb1;
  int stepcount;
  int nsteps=3;
  int sol;
  ha_cgetype vpercents=1.0,perprecis=0;
  FILE* solution;
  int maxsol=3;
              uvadd *counteqs= (uvadd *) calloc (ndblock+1,sizeof(uvadd));
              uvadd *counteqnoadds= (uvadd *) calloc (ndblock,sizeof(uvadd));
              uvadd *countvarintra1s= (uvadd *) calloc (ndblock+1,sizeof(uvadd));
              memcpy(counteqs,counteq,(ndblock+1)*sizeof(uvadd));
              memcpy(counteqnoadds,counteqnoadd,(ndblock)*sizeof(uvadd));
              memcpy(countvarintra1s,countvarintra1,(ndblock+1)*sizeof(uvadd));
    gettimeofday(&begintime, NULL);
    ha_cgetype *xc0= (ha_cgetype *) calloc (1,sizeof(ha_cgetype));
    ha_cgetype *xc12= (ha_cgetype *) calloc (1,sizeof(ha_cgetype));
    ha_cgetype *xc24= (ha_cgetype *) calloc (1,sizeof(ha_cgetype));
    int *xc124= (int *) calloc (1,sizeof(int));
    ha_cgetype *clag1= (ha_cgetype *) calloc (nvarele,sizeof(ha_cgetype));
    ha_cgetype *varchange= (ha_cgetype *) calloc (nvarele,sizeof(ha_cgetype));
    for(subindx=0; subindx<subints; subindx++) {
      for(sol=0; sol<maxsol; sol++) {
        if(sol==0)nsteps=step1;
        if(sol==1) nsteps=(int)step1*kindx1;
        if(sol==2) nsteps=(int)step1*kindx2;
        vpercents=(ha_cgetype)100/nsteps;
        for(stepcount=0; stepcount<nsteps; stepcount++) {
          printf("rank %d subint %d sol %d stepcount %d nsteps %d\n",rank,subindx,sol,stepcount,nsteps);
          MPI_Barrier(PETSC_COMM_WORLD);
          ierr = PetscGetCPUTime(&time0);
          CHKERRQ(ierr);
          if(stepcount==0) {
            MPI_Barrier(PETSC_COMM_WORLD);
            if(!(subindx==0&&sol==0&&stepcount==0)) {
              if(nohsl) {
                VecCreate(PETSC_COMM_WORLD,&vece);
              }
              else {
                VecCreate(PETSC_COMM_SELF,&vece);
              }
              if(nohsl) {
                VecSetType(vece,VECMPI);
              }
              else {
                VecSetType(vece,VECSEQ);
              }
              if(nesteddbbd==1)VecSetSizes(vece,localsize,VecSize);
              else VecSetSizes(vece,PETSC_DECIDE,VecSize);
              VecSetOption(vece, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
            }
            if(sol==0)for(i=0; i<ncofele; i++) {
                ha_cofvar[i].var0=ha_cofvar[i].varval;
              }
            else for(i=0; i<ncofele; i++) {
                ha_cofvar[i].varval=ha_cofvar[i].var0;
              }
            printf("rank %d OK!!!\n",rank);
            ha_cofvar1=ha_cofvar+ncofele;
            for(i=0; i<nvar; i++) {
              if(ha_var[i].change_real) {
                for(tindx1=ha_var[i].begadd; tindx1<ha_var[i].matsize+ha_var[i].begadd; tindx1++) {
                  if(ha_cgeshock[tindx1].ShockId) {
                    if(sol==0) {
                      ha_cofvar1[tindx1].var0=ha_cofvar1[tindx1].varval;
                    }
                    else {
                      ha_cofvar1[tindx1].varval=ha_cofvar1[tindx1].var0;
                    }
                    ha_cofvar1[tindx1].csolpupd=ha_cgeshock[tindx1].ShockVal/nsteps;
                    VecSetValue(vece,ha_cgeshock[tindx1].ExoIndx,ha_cofvar1[tindx1].csolpupd,INSERT_VALUES);
                  }
                  else {
                    if(sol==0) {
                      ha_cofvar1[tindx1].var0=ha_cofvar1[tindx1].varval;
                    }
                    else {
                      ha_cofvar1[tindx1].varval=ha_cofvar1[tindx1].var0;
                    }
                  }
                }
              }
              else {
                for(tindx1=ha_var[i].begadd; tindx1<ha_var[i].matsize+ha_var[i].begadd; tindx1++) {
                  if(ha_cgeshock[tindx1].ShockId) {
                    if(sol==0) {
                      ha_cofvar1[tindx1].var0=ha_cofvar1[tindx1].varval;
                    }
                    else {
                      ha_cofvar1[tindx1].varval=ha_cofvar1[tindx1].var0;
                    }
                    temp2=ha_cgeshock[tindx1].ShockVal;//subints;
                    ha_cofvar1[tindx1].csolpupd=(100+(subindx+1)*temp2)/(100+subindx*temp2)-1;//ha_cgeshock[ha_var[i].begadd+j].ShockVal/nsteps;//(exp(log(1+ha_cgeshock[ha_var[i].begadd+j].ShockVal/100)/nsteps)-1)*100;
                    ha_cofvar1[tindx1].csolpupd*=vpercents;//nsteps*100;
                    VecSetValue(vece,ha_cgeshock[tindx1].ExoIndx,ha_cofvar1[tindx1].csolpupd,INSERT_VALUES);
                  }
                  else {
                    if(sol==0) {
                      ha_cofvar1[tindx1].var0=ha_cofvar1[tindx1].varval;
                    }
                    else {
                      ha_cofvar1[tindx1].varval=ha_cofvar1[tindx1].var0;
                    }
                  }
                }
              }
            }
            MPI_Barrier(PETSC_COMM_WORLD);
            ierr = VecAssemblyBegin(vece);
            CHKERRQ(ierr);
            ierr = VecAssemblyEnd(vece);
            CHKERRQ(ierr);
          }
          printf("OK!!!\n");
          if(rank==rank_hsl) {

            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempclag1");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
              printf("Error opening file\n");
            }
            fwrite(clag1, sizeof(ha_cgetype),nvarele, tempvar);
            fclose(tempvar);
            free(clag1);
            clag1=NULL;

            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempvarchange");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
              printf("Error opening file\n");
            }
            fwrite(varchange, sizeof(ha_cgetype),nvarele, tempvar);
            fclose(tempvar);
            free(varchange);
            varchange=NULL;
          }
          MPI_Barrier(PETSC_COMM_WORLD);

          strcpy(commsyntax,"equation");

          if(nohsl) {
            MatCreate(PETSC_COMM_WORLD,&A);
          }
          else {
            MatCreate(PETSC_COMM_SELF,&A);
          }
          if(nesteddbbd==1)MatSetSizes(A,localsize,localsize,VecSize,VecSize);
          else MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,VecSize,VecSize);
          if(nohsl) {
            MatSetType(A,MATMPIAIJ);
          }
          else {
            MatSetType(A,MATSEQAIJ);
          }
          if(nohsl) {
            MatMPIAIJSetPreallocation(A,dnz,dnnz,onz,onnz);
          }
          else {
            MatSeqAIJSetPreallocation(A,dnz,dnnz);
          }

          if(nohsl) {
            MatCreate(PETSC_COMM_WORLD,&B);
          }
          else {
            MatCreate(PETSC_COMM_SELF,&B);
          }
          if(nesteddbbd==1)MatSetSizes(B,localsize,localsize,VecSize,VecSize);
          else MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,VecSize,VecSize);
          if(nohsl) {
            MatSetType(B,MATMPIAIJ);
          }
          else {
            MatSetType(B,MATSEQAIJ);
          }
          if(nohsl) {
            MatMPIAIJSetPreallocation(B,dnzB,dnnzB,onzB,onnzB);
          }
          else {
            MatSeqAIJSetPreallocation(B,dnzB,dnnzB);
          }
          printf("OKB!!!\n");

          if(rank==rank_hsl) {
            HaNewMatVal(tabfile,commsyntax,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele,ha_cgeshock,ndblock,alltimeset,allregset,ha_eqadd,counteq,nintraeq,A,B);
          }
          if(rank==rank_hsl) {
            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempvar");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
              printf("Error opening file\n");
            }
            fwrite(ha_cofvar, sizeof(ha_cgevar),ncofele+nvarele, tempvar);
            fclose(tempvar);
            free(*ha_cofvar2);
            *ha_cofvar2=NULL;
            ha_cofvar=*ha_cofvar2;

            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempshock");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
              printf("Error opening file\n");
            }
            fwrite(ha_cgeshock, sizeof(ha_cgeexovar),nvarele, tempvar);
            fclose(tempvar);
            free(*ha_cgeshock2);
            *ha_cgeshock2=NULL;
            ha_cgeshock=*ha_cgeshock2;
          }

          MPI_Barrier(PETSC_COMM_WORLD);
          ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
          CHKERRQ(ierr);
          ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
          CHKERRQ(ierr);
          ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);
          CHKERRQ(ierr);
          ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);
          CHKERRQ(ierr);
          ierr = VecDuplicate(vece,&vecb);
          CHKERRQ(ierr);
          if(rank==rank_hsl) {
            ierr = MatMult(B,vece,vecb);
            CHKERRQ(ierr);
          }
          ierr = VecDestroy(&vece);
          CHKERRQ(ierr);
          ierr = VecAssemblyBegin(vecb);
          CHKERRQ(ierr);
          ierr = VecAssemblyEnd(vecb);
          CHKERRQ(ierr);
          ierr = MatDestroy(&B);
          CHKERRQ(ierr);
          if(matsol>=2) {
            int *ha_rows= (int *) calloc (VecSize,sizeof(int));
            int *ha_cols= (int *) calloc (VecSize,sizeof(int));
            int *ha_ndblocks= (int *) calloc (ndblock,sizeof(int));
            time(&timestr);

            if(matsol==2) {
              HaDBBDMatOder(A,VecSize,mpisize,rank,Istart,Iend,nvarele,ha_eqadd,ha_rows,ha_cols,ndblock,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,cntl6);
              x1=realloc (x1,VecSize*sizeof(ha_cgetype));
              HaDBBDParSol(A,vecb,x1,VecSize,mpisize,rank,Istart,Iend,ha_rows,ha_cols,ndblock,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,laD,cntl3);//,iter
            }

            if(matsol==3) {
              presol=1;
              memcpy(counteq,counteqs,(ndblock+1)*sizeof(uvadd));
              memcpy(counteqnoadd,counteqnoadds,(ndblock)*sizeof(uvadd));
              memcpy(countvarintra1,countvarintra1s,(ndblock+1)*sizeof(uvadd));
              HaNDBBDMatOderPre(A,VecSize,mpisize,rank,Istart,Iend,nreg,ntime,nvarele,ha_eqadd,ha_rows,ha_cols,ndblock,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,laDi,cntl6,ndbbddrank1,presol);
              HaNDBBDParPre(A,vecb,x1,VecSize,mpisize,rank,Istart,Iend,ha_rows,ha_cols,ndblock,nreg,ntime,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,presol);//,iter
              presol=0;
              HaNDBBDMatOder(A,VecSize,mpisize,rank,Istart,Iend,nreg,ntime,nvarele,ha_eqadd,ha_rows,ha_cols,ndblock,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,laDi,cntl6,ndbbddrank1,presol);
              x1=realloc (x1,VecSize*sizeof(ha_cgetype));
              HaNDBBDParSol(A,vecb,x1,VecSize,mpisize,rank,Istart,Iend,ha_rows,ha_cols,ndblock,nreg,ntime,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,presol);//,iter
            }

            time(&timeend);
            MPI_Barrier(PETSC_COMM_WORLD);
            ierr = PetscGetCPUTime(&time1);
            ierr = PetscPrintf(PETSC_COMM_WORLD,"One step solution %f\n",time1-time0);
            if(rank==0)printf("One step calculation time %f\n",difftime(timeend,timestr));
            free(ha_rows);
            free(ha_cols);
            free(ha_ndblocks);
            printf("rank %d\n",rank);
            MPI_Barrier(PETSC_COMM_WORLD);
          }
          else {
            if(matsol==1) {
              if(rank==rank_hsl) {
                Mat_SeqAIJ *aa=(Mat_SeqAIJ*)A->data;
                ai= aa->i;
                aj= aa->j;
                vals=aa->a;
                nz01=aa->nz;
                count=0;
                for(i=0; i<nz01; i++) if(vals[i]!=0) {
                    count++;
                  }
              }
              int *irn=(int *) calloc (count,sizeof(int));
              int *irn1=(int *) calloc (nz01,sizeof(int));
              int *jcn=(int *) calloc (count,sizeof(int));
              ha_cgetype *values= (ha_cgetype *) calloc (count,sizeof(ha_cgetype));
              if(rank==rank_hsl) {
                for(i=0; i<VecSize-1; i++)for(j=ai[i]; j<ai[i+1]; j++) {
                    irn1[j]=i+1;
                  }
                for(j=ai[VecSize-1]; j<nz01; j++) {
                  irn1[j]=VecSize;
                }
                j=0;
                for(i=0; i<nz01; i++) if(vals[i]!=0) {
                    irn[j]=irn1[i];
                    jcn[j]=aj[i]+1;
                    values[j]=vals[i];
                    j++;
                  }
              }
              ierr = MatDestroy(&A);
              CHKERRQ(ierr);
              free(irn1);
              b1=realloc (b1,VecSize*sizeof(ha_cgetype));
              if(rank==rank_hsl) {
                VecGetArray(vecb,&vals);
                for(i=0; i<VecSize; i++) {
                  b1[i]=vals[i];
                }
              }
              ierr = VecDestroy(&vecb);
              CHKERRQ(ierr);
              int *neleperrow= (int *) calloc (VecSize,sizeof(int));
              int *ai1= (int *) calloc (VecSize,sizeof(int));
              if(rank==rank_hsl) {
                j=1;
                k=0,m=1;
                for(i=1; i<count; i++) {
                  if(irn[i]-irn[i-1]>0) {
                    neleperrow[k]=j;
                    ai1[k]=m;
                    j=1;
                    m=i+1;
                    k++;
                  }
                  else {
                    j++;
                  }
                }
                neleperrow[k]=j;
                ai1[k]=ai1[k-1]+neleperrow[k-1];
              }
              indata[1]=VecSize;//.m
              indata[0]=count;//.nz
              ptx = indata;
              x1=realloc (x1,VecSize*sizeof(ha_cgetype));
              ierr = PetscGetCPUTime(&time1);
              CHKERRQ(ierr);
              ierr = PetscPrintf(PETSC_COMM_WORLD,"Prepare time %f\n",time1-time0);
              CHKERRQ(ierr);
              ierr = PetscGetCPUTime(&time0);
              CHKERRQ(ierr);
              if(mc66!=0)spec48_single_(ptx,irn,jcn,b1,values,x1,neleperrow,ai1,&fcomm);
              free(irn);
              if(mc66==0)spec48_nomc66_(ptx,jcn,b1,values,x1,neleperrow,&fcomm,counteq,countvarintra1);
              ierr = PetscGetCPUTime(&time1);
              CHKERRQ(ierr);
              ierr = PetscPrintf(PETSC_COMM_WORLD,"LU time %f\n",time1-time0);
              CHKERRQ(ierr);
              ierr = PetscGetCPUTime(&time0);
              CHKERRQ(ierr);
              free(jcn);
              free(values);
              free(neleperrow);
              free(ai1);
              free(b1);
              b1=NULL;
            }
            else {
              if(rank==rank_hsl) {
                Mat_SeqAIJ *aa=(Mat_SeqAIJ*)A->data;
                ai= aa->i;
                aj= aa->j;
                vals=aa->a;
                nz01=aa->nz;
                count=0;
                for(i=0; i<nz01; i++) if(vals[i]!=0) {
                    count++;
                  }

              }
              lasize=ceil((laA/100.0)*count);
              int *irn=(int *) calloc (lasize,sizeof(int));
              int *irn1=(int *) calloc (nz01,sizeof(int));
              int *jcn=(int *) calloc (lasize,sizeof(int));
              ha_cgetype *values= (ha_cgetype *) calloc (lasize,sizeof(ha_cgetype));
              if(rank==rank_hsl) {
                for(i=0; i<VecSize-1; i++)for(j=ai[i]; j<ai[i+1]; j++) {
                    irn1[j]=i+1;
                  }
                for(j=ai[VecSize-1]; j<nz01; j++) {
                  irn1[j]=VecSize;
                }
                j=0;
                for(i=0; i<nz01; i++) if(vals[i]!=0) {
                    irn[j]=irn1[i];
                    jcn[j]=aj[i]+1;
                    values[j]=vals[i];
                    j++;
                  }
                VecGetArray(vecb,&vals);
              }
              ierr = MatDestroy(&A);
              CHKERRQ(ierr);
              free(irn1);
              x1=realloc (x1,VecSize*sizeof(ha_cgetype));
              ierr = PetscGetCPUTime(&time1);
              CHKERRQ(ierr);
              ierr = PetscPrintf(PETSC_COMM_WORLD,"Prepare time %f\n",time1-time0);
              CHKERRQ(ierr);
              ierr = PetscGetCPUTime(&time0);
              CHKERRQ(ierr);
              int *insize=(int *) calloc (4,sizeof(int));
              insize[0]=VecSize;
              insize[1]=VecSize;
              insize[2]=count;
              insize[3]=laA;
              if(rank==rank_hsl)spec48_ssol2la_(insize,irn,jcn,values,vals,x1);
              ierr = VecDestroy(&vecb);
              CHKERRQ(ierr);
              free(insize);
              ierr = PetscGetCPUTime(&time1);
              CHKERRQ(ierr);
              ierr = PetscPrintf(PETSC_COMM_WORLD,"LU time %f\n",time1-time0);
              CHKERRQ(ierr);
              ierr = PetscGetCPUTime(&time0);
              CHKERRQ(ierr);
              free(irn);
              free(jcn);
              free(values);
            }
          }
          if(rank==rank_hsl) {
            if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
              printf("Error opening file\n");
            }
            *ha_cgeshock2=(ha_cgeexovar*)realloc (*ha_cgeshock2,(nvarele)*sizeof(ha_cgeexovar));
            freadresult=fread(*ha_cgeshock2, sizeof(ha_cgeexovar),nvarele, tempvar);
            fclose(tempvar);
            remove(tempfilenam);
            ha_cgeshock=*ha_cgeshock2;

            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempvar");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
              printf("Error opening file\n");
            }
            *ha_cofvar2=(ha_cgevar*)realloc (*ha_cofvar2,(ncofele+nvarele)*sizeof(ha_cgevar));
            freadresult=fread(*ha_cofvar2, sizeof(ha_cgevar),ncofele+nvarele, tempvar);
            fclose(tempvar);
            remove(tempfilenam);
            ha_cofvar=*ha_cofvar2;

            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempclag1");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
              printf("Error opening file\n");
            }
            clag1=realloc (clag1,(nvarele)*sizeof(ha_cgetype));
            freadresult=fread(clag1, sizeof(ha_cgetype),nvarele, tempvar);
            fclose(tempvar);
            remove(tempfilenam);

            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempvarchange");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
              printf("Error opening file\n");
            }
            varchange=realloc (varchange,(nvarele)*sizeof(ha_cgetype));
            freadresult=fread(varchange, sizeof(ha_cgetype),nvarele, tempvar);
            fclose(tempvar);
            remove(tempfilenam);

          }
          printf("sol %d stepcount %d\n\n",sol,stepcount);
          MPI_Barrier(PETSC_COMM_WORLD);
          if(nohsl) {
            VecCreate(PETSC_COMM_WORLD,&vece);
          }
          else {
            VecCreate(PETSC_COMM_SELF,&vece);
          }
          if(nohsl) {
            VecSetType(vece,VECMPI);
          }
          else {
            VecSetType(vece,VECSEQ);
          }
          if(nesteddbbd==1)VecSetSizes(vece,localsize,VecSize);
          else VecSetSizes(vece,PETSC_DECIDE,VecSize);
          VecSetOption(vece, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
          ha_cofvar1=ha_cofvar+ncofele;
          if(stepcount==0) {
            for(i=0; i<nvar; i++) {
              if(ha_var[i].change_real) {
                for(tindx1=ha_var[i].begadd; tindx1<ha_var[i].matsize+ha_var[i].begadd; tindx1++) {
                  if(ha_cgeshock[tindx1].ShockId) {
                    ha_cofvar1[tindx1].varval+=ha_cofvar1[tindx1].csolpupd;
                    varchange[tindx1]=ha_cofvar1[tindx1].csolpupd;
                    VecSetValue(vece,ha_cgeshock[tindx1].ExoIndx,ha_cofvar1[tindx1].csolpupd,INSERT_VALUES);
                  }
                  else {
                    varchange[tindx1]=x1[ha_cgeshock[tindx1].ExoIndx];
                    ha_cofvar1[tindx1].varval+=x1[ha_cgeshock[tindx1].ExoIndx];
                    ha_cofvar1[tindx1].csolpupd=x1[ha_cgeshock[tindx1].ExoIndx];
                    clag1[tindx1]=0;
                  }
                }
              }
              else {
                for(tindx1=ha_var[i].begadd; tindx1<ha_var[i].matsize+ha_var[i].begadd; tindx1++) {
                  if(ha_cgeshock[tindx1].ShockId) {
                    varchange[tindx1]=ha_cofvar1[tindx1].csolpupd;
                    ha_cofvar1[tindx1].varval*=(1+ha_cofvar1[tindx1].csolpupd/100);
                    VecSetValue(vece,ha_cgeshock[tindx1].ExoIndx,ha_cofvar1[tindx1].csolpupd/(1+ha_cofvar1[tindx1].csolpupd/100),INSERT_VALUES);
                  }
                  else {
                    varchange[tindx1]=x1[ha_cgeshock[tindx1].ExoIndx];
                    ha_cofvar1[tindx1].csolpupd=x1[ha_cgeshock[tindx1].ExoIndx];
                    ha_cofvar1[tindx1].varval*=(1+ha_cofvar1[tindx1].csolpupd/100);
                    clag1[tindx1]=0;
                  }
                }
              }
            }
          }
          else {
            for(i=0; i<nvar; i++) {
              if(ha_var[i].change_real) {
                for(tindx1=ha_var[i].begadd; tindx1<ha_var[i].matsize+ha_var[i].begadd; tindx1++) {
                  if(ha_cgeshock[tindx1].ShockId) {
                    ha_cofvar1[tindx1].varval+=ha_cofvar1[tindx1].csolpupd;
                    varchange[tindx1]+=ha_cofvar1[tindx1].csolpupd;
                    VecSetValue(vece,ha_cgeshock[tindx1].ExoIndx,ha_cofvar1[tindx1].csolpupd,INSERT_VALUES);
                  }
                  else {
                    temp1=ha_cofvar1[tindx1].varval;//change;
                    varchange[tindx1]=clag1[tindx1]+2*x1[ha_cgeshock[tindx1].ExoIndx];//+=x1[ha_cgeshock[ha_var[i].begadd+j].ExoIndx];//
                    ha_cofvar1[tindx1].csolpupd=x1[ha_cgeshock[tindx1].ExoIndx];//ha_cofvar[ncofele+ha_var[i].begadd+j].varchange-temp1;
                    ha_cofvar1[tindx1].varval=clag1[tindx1]+2*x1[ha_cgeshock[tindx1].ExoIndx];//ha_cofvar[ncofele+ha_var[i].begadd+j].varchange-temp1;//ha_cofvar[ncofele+ha_var[i].begadd+j].csolpupd;
                    clag1[tindx1]=temp1;
                  }
                }
              }
              else {
                for(tindx1=ha_var[i].begadd; tindx1<ha_var[i].matsize+ha_var[i].begadd; tindx1++) {
                  if(ha_cgeshock[tindx1].ShockId) {
                    temp2=ha_cgeshock[tindx1].ShockVal;//subints;
                    temp1=(100+(subindx+1)*temp2)/(100+subindx*temp2)-1;
                    temp1*=vpercents;
                    ha_cofvar1[tindx1].csolpupd=temp1/(1+varchange[tindx1]/100);
                    varchange[tindx1]+=temp1;//*(1+ha_cofvar[ncofele+ha_var[i].begadd+j].varchange/100)
                    ha_cofvar1[tindx1].varval=(1+varchange[tindx1]/100)*ha_cofvar1[tindx1].var0;
                    VecSetValue(vece,ha_cgeshock[tindx1].ExoIndx,temp1/(1+varchange[tindx1]/100),INSERT_VALUES);
                  }
                  else {
                    temp1=varchange[tindx1];
                    varchange[tindx1]=clag1[tindx1]+2*x1[ha_cgeshock[tindx1].ExoIndx]*(100+temp1)/100;//+=x1[ha_cgeshock[ha_var[i].begadd+j].ExoIndx]*(1+temp1/100);//
                    ha_cofvar1[tindx1].csolpupd=x1[ha_cgeshock[tindx1].ExoIndx];//(ha_cofvar[ncofele+ha_var[i].begadd+j].varchange-temp1)/(1+temp1/100);
                    ha_cofvar1[tindx1].varval=(100+varchange[tindx1])/100*ha_cofvar1[tindx1].var0;
                    clag1[tindx1]=temp1;
                  }
                }
              }
            }
          }
          free(x1);
          x1=NULL;
          MPI_Barrier(PETSC_COMM_WORLD);
          ierr = VecAssemblyBegin(vece);
          CHKERRQ(ierr);
          MPI_Barrier(PETSC_COMM_WORLD);
          ierr = VecAssemblyEnd(vece);
          CHKERRQ(ierr);
          if(rank==rank_hsl) {
            if(stepcount==0) {
              hnew_update(tabfile,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele);
            }
            else {
              hnew_mupdate(tabfile,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele);
            }
            if(Isbiupd==1)hnew_biupd(rank,tabfile,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele,ha_cgeshock,nvarele,laA,subints,1,0,nsteps);
            if(Isbiupd==2)hnew_biupd(rank,tabfile,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele,ha_cgeshock,nvarele,laA,subints,1,2,nsteps);
            strcpy(commsyntax,"formula");
            IsIni=false;
            hnew_calcff(tabfile,commsyntax,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele,IsIni);
          }
          MPI_Barrier(PETSC_COMM_WORLD);
          ierr = PetscGetCPUTime(&time1);
          CHKERRQ(ierr);
          ierr = PetscPrintf(PETSC_COMM_WORLD,"Update time %f\n",time1-time0);
          CHKERRQ(ierr);
          ierr = PetscGetCPUTime(&time0);
          CHKERRQ(ierr);
          printf("OKKL!\n");
        }

        strcpy(commsyntax,"equation");
        if(rank==rank_hsl) {


          strcpy(tempfilenam,temdir);
          strcat(tempfilenam,"_tempclag1");
          sprintf(tempchar, "%d",rank);
          strcat(tempfilenam,tempchar);
          strcat(tempfilenam,".bin");
          if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
            printf("Error opening file\n");
          }
          fwrite(clag1, sizeof(ha_cgetype),nvarele, tempvar);
          fclose(tempvar);
          free(clag1);
          clag1=NULL;

          strcpy(tempfilenam,temdir);
          strcat(tempfilenam,"_tempvarchange");
          sprintf(tempchar, "%d",rank);
          strcat(tempfilenam,tempchar);
          strcat(tempfilenam,".bin");
          if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
            printf("Error opening file\n");
          }
          fwrite(varchange, sizeof(ha_cgetype),nvarele, tempvar);
          fclose(tempvar);
          free(varchange);
          varchange=NULL;

        }

        if(nohsl)MPI_Barrier(PETSC_COMM_WORLD);
        if(nohsl) {
          MatCreate(PETSC_COMM_WORLD,&A);
        }
        else {
          MatCreate(PETSC_COMM_SELF,&A);
        }
        if(nesteddbbd==1)MatSetSizes(A,localsize,localsize,VecSize,VecSize);
        else MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,VecSize,VecSize);
        if(nohsl) {
          MatSetType(A,MATMPIAIJ);
        }
        else {
          MatSetType(A,MATSEQAIJ);
        }
        if(nohsl) {
          MatMPIAIJSetPreallocation(A,dnz,dnnz,onz,onnz);
        }
        else {
          MatSeqAIJSetPreallocation(A,dnz,dnnz);
        }

        if(nohsl) {
          MatCreate(PETSC_COMM_WORLD,&B);
        }
        else {
          MatCreate(PETSC_COMM_SELF,&B);
        }
        if(nesteddbbd==1)MatSetSizes(B,localsize,localsize,VecSize,VecSize);
        else MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,VecSize,VecSize);
        if(nohsl) {
          MatSetType(B,MATMPIAIJ);
        }
        else {
          MatSetType(B,MATSEQAIJ);
        }
        if(nohsl) {
          MatMPIAIJSetPreallocation(B,dnzB,dnnzB,onzB,onnzB);
        }
        else {
          MatSeqAIJSetPreallocation(B,dnzB,dnnzB);
        }

        if(rank==rank_hsl) {
          HaNewMatVal(tabfile,commsyntax,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele,ha_cgeshock,ndblock,alltimeset,allregset,ha_eqadd,counteq,nintraeq,A,B);
        }

        if(rank==rank_hsl) {
          strcpy(tempfilenam,temdir);
          strcat(tempfilenam,"_tempvar");
          sprintf(tempchar, "%d",rank);
          strcat(tempfilenam,tempchar);
          strcat(tempfilenam,".bin");
          if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
            printf("Error opening file\n");
          }
          fwrite(ha_cofvar, sizeof(ha_cgevar),ncofele+nvarele, tempvar);
          fclose(tempvar);
          free(*ha_cofvar2);
          *ha_cofvar2=NULL;
          ha_cofvar=*ha_cofvar2;

          strcpy(tempfilenam,temdir);
          strcat(tempfilenam,"_tempshock");
          sprintf(tempchar, "%d",rank);
          strcat(tempfilenam,tempchar);
          strcat(tempfilenam,".bin");
          if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
            printf("Error opening file\n");
          }
          fwrite(ha_cgeshock, sizeof(ha_cgeexovar),nvarele, tempvar);
          fclose(tempvar);
          free(*ha_cgeshock2);
          *ha_cgeshock2=NULL;
          ha_cgeshock=*ha_cgeshock2;
        }

        MPI_Barrier(PETSC_COMM_WORLD);
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
        CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
        CHKERRQ(ierr);
        ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);
        CHKERRQ(ierr);
        ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);
        CHKERRQ(ierr);
        ierr = VecDuplicate(vece,&vecb);
        CHKERRQ(ierr);
        if(rank==rank_hsl) {
          ierr = MatMult(B,vece,vecb);
          CHKERRQ(ierr);
        }
        ierr = VecDestroy(&vece);
        CHKERRQ(ierr);
        ierr = VecAssemblyBegin(vecb);
        CHKERRQ(ierr);
        ierr = VecAssemblyEnd(vecb);
        CHKERRQ(ierr);
        ierr = MatDestroy(&B);
        CHKERRQ(ierr);
        if(matsol>=2) {
          int *ha_rows= (int *) calloc (VecSize,sizeof(int));
          int *ha_cols= (int *) calloc (VecSize,sizeof(int));
          int *ha_ndblocks= (int *) calloc (ndblock,sizeof(int));
          time(&timestr);

          if(matsol==2) {
            HaDBBDMatOder(A,VecSize,mpisize,rank,Istart,Iend,nvarele,ha_eqadd,ha_rows,ha_cols,ndblock,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,cntl6);
            x1=realloc (x1,VecSize*sizeof(ha_cgetype));
            HaDBBDParSol(A,vecb,x1,VecSize,mpisize,rank,Istart,Iend,ha_rows,ha_cols,ndblock,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,laD,cntl3);//,iter
          }

          if(matsol==3) {
            presol=1;
            memcpy(counteq,counteqs,(ndblock+1)*sizeof(uvadd));
            memcpy(counteqnoadd,counteqnoadds,(ndblock)*sizeof(uvadd));
            memcpy(countvarintra1,countvarintra1s,(ndblock+1)*sizeof(uvadd));
            HaNDBBDMatOderPre(A,VecSize,mpisize,rank,Istart,Iend,nreg,ntime,nvarele,ha_eqadd,ha_rows,ha_cols,ndblock,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,laDi,cntl6,ndbbddrank1,presol);
            HaNDBBDParPre(A,vecb,x1,VecSize,mpisize,rank,Istart,Iend,ha_rows,ha_cols,ndblock,nreg,ntime,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,presol);//,iter
            presol=0;
            HaNDBBDMatOder(A,VecSize,mpisize,rank,Istart,Iend,nreg,ntime,nvarele,ha_eqadd,ha_rows,ha_cols,ndblock,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,laDi,cntl6,ndbbddrank1,presol);
            x1=realloc (x1,VecSize*sizeof(ha_cgetype));
            HaNDBBDParSol(A,vecb,x1,VecSize,mpisize,rank,Istart,Iend,ha_rows,ha_cols,ndblock,nreg,ntime,ha_ndblocks,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,presol);//,iter
          }

          time(&timeend);
          MPI_Barrier(PETSC_COMM_WORLD);
          ierr = PetscGetCPUTime(&time1);
          ierr = PetscPrintf(PETSC_COMM_WORLD,"One step solution %f\n",time1-time0);
          if(rank==0)printf("One step calculation time %f\n",difftime(timeend,timestr));
          free(ha_rows);
          free(ha_cols);
          free(ha_ndblocks);
          printf("rank %d\n",rank);
          MPI_Barrier(PETSC_COMM_WORLD);
        }
        else {
          if(matsol==1) {
            if(rank==rank_hsl) {
              Mat_SeqAIJ *aa=(Mat_SeqAIJ*)A->data;
              ai= aa->i;
              aj= aa->j;
              vals=aa->a;
              nz01=aa->nz;
              count=0;
              for(i=0; i<nz01; i++) if(vals[i]!=0) {
                  count++;
                }
            }
            int *irn=(int *) calloc (count,sizeof(int));
            int *irn1=(int *) calloc (nz01,sizeof(int));
            int *jcn=(int *) calloc (count,sizeof(int));
            ha_cgetype *values= (ha_cgetype *) calloc (count,sizeof(ha_cgetype));
            if(rank==rank_hsl) {
              for(i=0; i<VecSize-1; i++)for(j=ai[i]; j<ai[i+1]; j++) {
                  irn1[j]=i+1;
                }
              for(j=ai[VecSize-1]; j<nz01; j++) {
                irn1[j]=VecSize;
              }
              j=0;
              for(i=0; i<nz01; i++) if(vals[i]!=0) {
                  irn[j]=irn1[i];
                  jcn[j]=aj[i]+1;
                  values[j]=vals[i];
                  j++;
                }
            }
            ierr = MatDestroy(&A);
            CHKERRQ(ierr);
            free(irn1);
            b1=realloc (b1,VecSize*sizeof(ha_cgetype));
            if(rank==rank_hsl) {
              VecGetArray(vecb,&vals);
              for(i=0; i<VecSize; i++) {
                b1[i]=vals[i];
              }
            }
            ierr = VecDestroy(&vecb);
            CHKERRQ(ierr);
            int *neleperrow= (int *) calloc (VecSize,sizeof(int));
            int *ai1= (int *) calloc (VecSize,sizeof(int));
            if(rank==rank_hsl) {
              j=1;
              k=0,m=1;
              for(i=1; i<count; i++) {
                if(irn[i]-irn[i-1]>0) {
                  neleperrow[k]=j;
                  ai1[k]=m;
                  j=1;
                  m=i+1;
                  k++;
                }
                else {
                  j++;
                }
              }
              neleperrow[k]=j;
              ai1[k]=ai1[k-1]+neleperrow[k-1];
            }
            indata[1]=VecSize;//.m
            indata[0]=count;//.nz
            ptx = indata;//&
            x1=realloc (x1,VecSize*sizeof(ha_cgetype));
            ierr = PetscGetCPUTime(&time1);
            CHKERRQ(ierr);
            ierr = PetscPrintf(PETSC_COMM_WORLD,"Prepare time %f\n",time1-time0);
            CHKERRQ(ierr);
            ierr = PetscGetCPUTime(&time0);
            CHKERRQ(ierr);
            if(mc66!=0)spec48_single_(ptx,irn,jcn,b1,values,x1,neleperrow,ai1,&fcomm);
            free(irn);
            if(mc66==0)spec48_nomc66_(ptx,jcn,b1,values,x1,neleperrow,&fcomm,counteq,countvarintra1);
            ierr = PetscGetCPUTime(&time1);
            CHKERRQ(ierr);
            ierr = PetscPrintf(PETSC_COMM_WORLD,"LU time %f\n",time1-time0);
            CHKERRQ(ierr);
            ierr = PetscGetCPUTime(&time0);
            CHKERRQ(ierr);
            free(jcn);
            free(values);
            free(neleperrow);
            free(ai1);
            free(b1);
            b1=NULL;
          }
          else {
            if(rank==rank_hsl) {
              Mat_SeqAIJ *aa=(Mat_SeqAIJ*)A->data;
              ai= aa->i;
              aj= aa->j;
              vals=aa->a;
              nz01=aa->nz;
              count=0;
              for(i=0; i<nz01; i++) if(vals[i]!=0) {
                  count++;
                }
            }
            lasize=ceil((laA/100.0)*count);
            int *irn=(int *) calloc (lasize,sizeof(int));
            int *irn1=(int *) calloc (nz01,sizeof(int));
            int *jcn=(int *) calloc (lasize,sizeof(int));
            ha_cgetype *values= (ha_cgetype *) calloc (lasize,sizeof(ha_cgetype));
            if(rank==rank_hsl) {
              for(i=0; i<VecSize-1; i++)for(j=ai[i]; j<ai[i+1]; j++) {
                  irn1[j]=i+1;
                }
              for(j=ai[VecSize-1]; j<nz01; j++) {
                irn1[j]=VecSize;
              }
              j=0;
              for(i=0; i<nz01; i++) if(vals[i]!=0) {
                  irn[j]=irn1[i];
                  jcn[j]=aj[i]+1;
                  values[j]=vals[i];
                  j++;
                }
              VecGetArray(vecb,&vals);
            }
            ierr = MatDestroy(&A);
            CHKERRQ(ierr);
            free(irn1);
            x1=realloc (x1,VecSize*sizeof(ha_cgetype));
            ierr = PetscGetCPUTime(&time1);
            CHKERRQ(ierr);
            ierr = PetscPrintf(PETSC_COMM_WORLD,"Prepare time %f\n",time1-time0);
            CHKERRQ(ierr);
            ierr = PetscGetCPUTime(&time0);
            CHKERRQ(ierr);
            int *insize=(int *) calloc (4,sizeof(int));
            insize[0]=VecSize;
            insize[1]=VecSize;
            insize[2]=count;
            insize[3]=laA;
            if(rank==rank_hsl)spec48_ssol2la_(insize,irn,jcn,values,vals,x1);
            free(insize);
            ierr = PetscGetCPUTime(&time1);
            CHKERRQ(ierr);
            ierr = PetscPrintf(PETSC_COMM_WORLD,"LU time %f\n",time1-time0);
            CHKERRQ(ierr);
            ierr = PetscGetCPUTime(&time0);
            CHKERRQ(ierr);
            ierr = VecDestroy(&vecb);
            CHKERRQ(ierr);
            free(irn);
            free(jcn);
            free(values);
          }
        }
        if(rank==rank_hsl) {
          if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
            printf("Error opening file\n");
          }
          *ha_cgeshock2=(ha_cgeexovar*)realloc (*ha_cgeshock2,(nvarele)*sizeof(ha_cgeexovar));
          freadresult=fread(*ha_cgeshock2, sizeof(ha_cgeexovar),nvarele, tempvar);
          fclose(tempvar);
          remove(tempfilenam);
          ha_cgeshock=*ha_cgeshock2;

          strcpy(tempfilenam,temdir);
          strcat(tempfilenam,"_tempvar");
          sprintf(tempchar, "%d",rank);
          strcat(tempfilenam,tempchar);
          strcat(tempfilenam,".bin");
          if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
            printf("Error opening file\n");
          }
          *ha_cofvar2=(ha_cgevar*)realloc (*ha_cofvar2,(ncofele+nvarele)*sizeof(ha_cgevar));
          freadresult=fread(*ha_cofvar2, sizeof(ha_cgevar),ncofele+nvarele, tempvar);
          fclose(tempvar);
          remove(tempfilenam);
          ha_cofvar=*ha_cofvar2;

          strcpy(tempfilenam,temdir);
          strcat(tempfilenam,"_tempclag1");
          sprintf(tempchar, "%d",rank);
          strcat(tempfilenam,tempchar);
          strcat(tempfilenam,".bin");
          if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
            printf("Error opening file\n");
          }
          clag1=realloc (clag1,(nvarele)*sizeof(ha_cgetype));
          freadresult=fread(clag1, sizeof(ha_cgetype),nvarele, tempvar);
          fclose(tempvar);
          remove(tempfilenam);

          strcpy(tempfilenam,temdir);
          strcat(tempfilenam,"_tempvarchange");
          sprintf(tempchar, "%d",rank);
          strcat(tempfilenam,tempchar);
          strcat(tempfilenam,".bin");
          if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
            printf("Error opening file\n");
          }
          varchange=realloc (varchange,(nvarele)*sizeof(ha_cgetype));
          freadresult=fread(varchange, sizeof(ha_cgetype),nvarele, tempvar);
          fclose(tempvar);
          remove(tempfilenam);

        }
        ha_cofvar1=ha_cofvar+ncofele;
        for(i=0; i<nvar; i++) {
          if(ha_var[i].change_real) {
            for(tindx1=ha_var[i].begadd; tindx1<ha_var[i].matsize+ha_var[i].begadd; tindx1++) {
              if(ha_cgeshock[tindx1].ShockId) {
                ha_cofvar1[tindx1].varval=0;
              }
              else {
                varchange[tindx1]=0.5*(varchange[tindx1]+clag1[tindx1]+x1[ha_cgeshock[tindx1].ExoIndx]);
                ha_cofvar1[tindx1].varval=0;//ha_cofvar[tindx2].var0+varchange[tindx1];//no distortion between steps
                clag1[tindx1]=0;
              }
            }
          }
          else {
            for(tindx1=ha_var[i].begadd; tindx1<ha_var[i].matsize+ha_var[i].begadd; tindx1++) {
              if(ha_cgeshock[tindx1].ShockId) {
                ha_cofvar1[tindx1].varval=0;
              }
              else {
                varchange[tindx1]=0.5*(varchange[tindx1]+clag1[tindx1]+x1[ha_cgeshock[tindx1].ExoIndx]*(1+varchange[tindx1]/100));
                ha_cofvar1[tindx1].varval=0;//ha_cofvar[tindx2].varval*varchange[tindx1]/100;
                clag1[tindx1]=0;
              }
            }
          }
        }
        if(rank==rank_hsl) {
          if(subindx!=0||sol!=0) {
            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempxcf");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
              printf("Error opening file\n");
            }
            *xcf2=(ha_cgetype*)realloc (*xcf2,(nvarele)*sizeof(ha_cgetype));
            xcf=*xcf2;
            freadresult=fread(xcf, sizeof(ha_cgetype),nvarele, tempvar);
            fclose(tempvar);
            remove(tempfilenam);

            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempxc12");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
              printf("Error opening file\n");
            }
            xc12=realloc (xc12,(nvarele)*sizeof(ha_cgetype));
            freadresult=fread(xc12, sizeof(ha_cgetype),nvarele, tempvar);
            fclose(tempvar);
            remove(tempfilenam);
            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempxc24");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
              printf("Error opening file\n");
            }
            xc24=realloc (xc24,(nvarele)*sizeof(ha_cgetype));
            freadresult=fread(xc24, sizeof(ha_cgetype),nvarele, tempvar);
            fclose(tempvar);
            remove(tempfilenam);
            
            xc0=realloc (xc0,(nvarele)*sizeof(ha_cgetype));
            if(subindx>0&&sol>0){
            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempxcO");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
              printf("Error opening file\n");
            }
            freadresult=fread(xc0, sizeof(ha_cgetype),nvarele, tempvar);
            fclose(tempvar);
            if(subindx==subints-1&&sol==maxsol-1)remove(tempfilenam);
            }

          }
          if(subindx==0&&sol==0) {
            xc0=realloc (xc0,nvarele*sizeof(ha_cgetype));
            *xcf2=(ha_cgetype*)realloc (*xcf2,nvarele*sizeof(ha_cgetype));
            xcf=*xcf2;
            for(i=0; i<nvarele; i++)xcf[i]=0;
          }
          if(sol==0)xc12=realloc (xc12,nvarele*sizeof(ha_cgetype));
          if(sol==0)xc24=realloc (xc24,nvarele*sizeof(ha_cgetype));
          if(subindx>0) {
            if(sol==0)for(i=0; i<nvarele; i++) xc0[i]=1+xcf[i]/100;//if(i==1287)printf("sol!!!!!!!!!!!!!!!!!! %d step %d xc %lf xc0 %lf k %d\n",sol,stepcount,1.0+xc[k]/100,xc0[i],i);}
            if(sol==0) {
              kval1=1.0/(kindx1*kindx1-1.0);
              kval2=1.0/(1-kindx1*kindx1)/(1.0-kindx2*kindx2);
              for(i=0; i<nvar; i++) {
                if(ha_var[i].change_real) {
                  for(k=ha_var[i].begadd; k<ha_var[i].matsize+ha_var[i].begadd; k++) {
                    xc12[k]=xcf[k]-varchange[k]*kval1;
                    xc24[k]=xcf[k];
                    xcf[k]+=varchange[k]*kval2;
                  }
                }
                else {
                  for(k=ha_var[i].begadd; k<ha_var[i].matsize+ha_var[i].begadd; k++) {
                    xc24[k]=xcf[k];
                    xc12[k]=xcf[k]-varchange[k]*xc0[k]*kval1;
                    xcf[k]+=varchange[k]*xc0[k]*kval2;//(100+xc0[k])*(100+varchange[k]/45)/100-100;//varchange[k]/45;
                  }
                }
              }
            }
            if(sol==1) {
              kval1=kindx1*kindx1/(kindx1*kindx1-1.0);
              kval2=kindx1*kindx1/(kindx2*kindx2-kindx1*kindx1);
              kval3=kindx1*kindx1*kindx1*kindx1/(kindx1*kindx1-kindx2*kindx2)/(1.0-kindx1*kindx1);
              for(i=0; i<nvar; i++) {
                if(ha_var[i].change_real) {
                  for(k=ha_var[i].begadd; k<ha_var[i].matsize+ha_var[i].begadd; k++) {
                    xc24[k]-=varchange[k]*kval2;
                    xc12[k]+=varchange[k]*kval1;
                    xcf[k]-=varchange[k]*kval3;
                  }
                }
                else {
                  for(k=ha_var[i].begadd; k<ha_var[i].matsize+ha_var[i].begadd; k++) {
                    xc24[k]-=varchange[k]*xc0[k]*kval2;
                    xc12[k]+=varchange[k]*xc0[k]*kval1;
                    xcf[k]-=varchange[k]*xc0[k]*kval3;//(100+xc0[k])*(100-20*varchange[k]/45)/100-100;
                  }
                }
              }
            }
            if(sol==2) {
              kval2=kindx2*kindx2/(kindx2*kindx2-kindx1*kindx1);
              kval3=kindx2*kindx2*kindx2*kindx2/(kindx1*kindx1-kindx2*kindx2)/(1.0-kindx2*kindx2);
              for(i=0; i<nvar; i++) {
                if(ha_var[i].change_real) {
                  for(k=ha_var[i].begadd; k<ha_var[i].matsize+ha_var[i].begadd; k++) {
                    xc24[k]+=varchange[k]*kval2;
                    xcf[k]+=varchange[k]*kval3;
                  }
                }
                else {
                  for(k=ha_var[i].begadd; k<ha_var[i].matsize+ha_var[i].begadd; k++) {
                    xc24[k]+=varchange[k]*xc0[k]*kval2;
                    xcf[k]+=varchange[k]*xc0[k]*kval3;//(100+xc0[k])*(100+64*varchange[k]/45)/100-100;
                  }
                }
              }
            }
          }
          else {
            if(sol==0){
              kval1=1.0/(kindx1*kindx1-1);
              kval2=1.0/(1-kindx1*kindx1)/(1-kindx2*kindx2);
              for(i=0; i<nvarele; i++) {
                xc12[i]=-varchange[i]*kval1;
                xcf[i]+=varchange[i]*kval2;
              }
            }
            if(sol==1) {
              kval1=kindx1*kindx1/(kindx1*kindx1-1.0);
              kval2=kindx1*kindx1/(kindx2*kindx2-kindx1*kindx1);
              kval3=kindx1*kindx1*kindx1*kindx1/(kindx1*kindx1-kindx2*kindx2)/(1.0-kindx1*kindx1);
              for(i=0; i<nvarele; i++) {
                xc24[i]=-varchange[i]*kval2;
                xc12[i]+=varchange[i]*kval1;
                xcf[i]-=varchange[i]*kval3;
              }
            }
            if(sol==2) {
              kval2=kindx2*kindx2/(kindx2*kindx2-kindx1*kindx1);
              kval3=kindx2*kindx2*kindx2*kindx2/(kindx1*kindx1-kindx2*kindx2)/(1.0-kindx2*kindx2);
              for(i=0; i<nvarele; i++) {
                xc24[i]+=varchange[i]*kval2;
                xcf[i]+=varchange[i]*kval3;
              }
            }
          }

          if(sol==maxsol-1){
          if(subindx==0){
          for(i=0; i<nvar; i++) {
            for(tindx1=ha_var[i].begadd; tindx1<ha_var[i].matsize+ha_var[i].begadd; tindx1++) {
               ha_cofvar1[tindx1].csolpupd=xcf[tindx1];
            }
          }
          }else{
          for(i=0; i<nvar; i++) {
            for(tindx1=ha_var[i].begadd; tindx1<ha_var[i].matsize+ha_var[i].begadd; tindx1++) {
               ha_cofvar1[tindx1].csolpupd=(100+xcf[tindx1])/xc0[tindx1]-100;
            }
          }
          }
          for(i=0; i<ncofele; i++) ha_cofvar[i].varval=ha_cofvar[i].var0;
          hnew_gupd(tabfile,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele);
          if(Isbiupd==1)hnew_biupd(rank,tabfile,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele,ha_cgeshock,nvarele,laA,subints,1,0,nsteps);
          if(Isbiupd==2)hnew_biupd(rank,tabfile,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele,ha_cgeshock,nvarele,laA,subints,1,2,nsteps);
          strcpy(commsyntax,"formula");
          IsIni=false;
          hnew_calcff(tabfile,commsyntax,ha_set,nset,ha_setele,ha_cof,ncof,ha_var,nvar,ha_cofvar,ncofele+nvarele,ncofele,IsIni);
          for(i=0; i<nvar; i++) {
            for(tindx1=ha_var[i].begadd; tindx1<ha_var[i].matsize+ha_var[i].begadd; tindx1++) {
               ha_cofvar1[tindx1].csolpupd=0;
            }
          }

          if(nohsl)MPI_Barrier(PETSC_COMM_WORLD);
          ierr = PetscGetCPUTime(&time1);
          CHKERRQ(ierr);
          ierr = PetscPrintf(PETSC_COMM_WORLD,"Last Update time %f\n",time1-time0);
          CHKERRQ(ierr);
          }

          
          if(!(subindx==subints-1&&sol==maxsol-1)) {
            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempxcf");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
              printf("Error opening file\n");
            }
            fwrite(xcf, sizeof(ha_cgetype),nvarele, tempvar);
            fclose(tempvar);
            free(*xcf2);
            *xcf2=NULL;
            xcf=*xcf2;
            
            if(subindx>0&&sol==0){            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempxcO");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
              printf("Error opening file\n");
            }
            fwrite(xc0, sizeof(ha_cgetype),nvarele, tempvar);
            fclose(tempvar);
            }
            free(xc0);
            xc0=NULL;
          }

    if(rank==rank_hsl&&sol==maxsol-1) {
            if(subindx==0){
              xc124=realloc (xc124,nvarele*sizeof(int));
              for(i=0; i<nvarele; i++)xc124[i]=6;
            }else{
            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempxc124");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
              printf("Error opening file\n");
            }
            xc124=realloc (xc124,nvarele*sizeof(int));
            freadresult=fread(xc124, sizeof(int),nvarele, tempvar);
            fclose(tempvar);
            remove(tempfilenam);
            }

      for(i=0; i<nvarele; i++){
        j=0;
        if(xc12[i]>0){
        while (xc12[i] >= 10){
          xc12[i] /= 10;
          j++;
        }
        }else{
        while (xc12[i] <= -10){
          xc12[i] /= 10;
          j++;
        }
        }
        xc24[i]/=pow(10,j);
        j=abs(floor((xc12[i]-xc24[i])*100000));
        if(j!=0){//}else {
          if(j<10){
            if(xc124[i]>5)xc124[i]=5;
          }else{
            if(j<100){
              if(xc124[i]>4)xc124[i]=4;
            }else {
              if(j<1000){
                if(xc124[i]>3)xc124[i]=3;
              }else {
                if(j<10000){
                  if(xc124[i]>2)xc124[i]=2;
                }else{
                  if(xc124[i]>1)xc124[i]=1;
                }
              }
            }
          }
        }
      }
      if(subindx!=subints-1){
            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempxc124");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
              printf("Error opening file\n");
            }
            fwrite(xc124, sizeof(int),nvarele, tempvar);
            fclose(tempvar);
            free(xc124);
            xc124=NULL;
      }
    }
            
          if(!(subindx==subints-1&&sol==maxsol-1)) {
            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempxc12");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
              printf("Error opening file\n");
            }
            fwrite(xc12, sizeof(ha_cgetype),nvarele, tempvar);
            fclose(tempvar);
            free(xc12);
            xc12=NULL;
            strcpy(tempfilenam,temdir);
            strcat(tempfilenam,"_tempxc24");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
              printf("Error opening file\n");
            }
            fwrite(xc24, sizeof(ha_cgetype),nvarele, tempvar);
            fclose(tempvar);
            free(xc24);
            xc24=NULL;
          }
          
          
          strcpy(solchar,temdir);
          if(subindx<10)strcat(solchar,"_tempsol0");
          else strcat(solchar,"_tempsol");
          sprintf(tempchar, "%d", subindx);
          strcat(solchar,tempchar);
          sprintf(tempchar, "%d", sol);
          strcat(solchar,tempchar);
          strcat(solchar,".bin");
          printf("solchar %s\n",solchar);
          if ( (solution = fopen(solchar, "wb")) == NULL ) {
            printf("Error opening file\n");
            return 1;
          }
          fwrite(xcf, sizeof(ha_cgetype),nvarele, solution);
          fclose(solution);
        }
        free(x1);
        x1=NULL;
      }
    }
    long int *precis= (long int *) calloc (6,sizeof(long int));
    
    if(rank==rank_hsl) {
      for(i=0; i<nvarele; i++){
   switch(xc124[i]) {
      case 6:
         precis[5]+=1;
         break;
      case 5:
         precis[4]+=1;
         break;
      case 4:
         precis[3]+=1;
         break;
      case 3:
         precis[2]+=1;
         break;
      case 2:
         precis[1]+=1;
         break;
      default :
         precis[0]+=1;
   }
    }
    if(rank==0)printf("Accurate at 6 digits        %ld\nAccurate at 5 digits        %ld\nAccurate at 4 digits        %ld\nAccurate at 3 digits        %ld\nAccurate at 2 digits        %ld\nAccurate at 1 digit or none %ld\n",precis[5],precis[4],precis[3],precis[2],precis[1],precis[0]);
    }
    free(precis);
    xc0=realloc (xc0,sizeof(ha_cgetype));
    free(xc0);
    xc0=NULL;
    free(xc12);
    xc12=NULL;
    free(xc24);
    xc24=NULL;
    free(xc124);
    xc124=NULL;
    free(clag1);
    free(varchange);
    gettimeofday(&endtime, NULL);
    if(rank==0)printf("Mmid method calculation time %f\n",(endtime.tv_sec - begintime.tv_sec)+((double)(endtime.tv_usec - begintime.tv_usec))/ 1000000);
              free(counteqs);
              free(counteqnoadds);
              free(countvarintra1s);
              ha_cofvar1=NULL;
              if(x1!=NULL)free(x1);
  }
