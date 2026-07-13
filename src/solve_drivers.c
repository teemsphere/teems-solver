#include <teems_solver.h>
#include <hsl_kernels.h>

bool solve_johansen(PetscBool nohsl,PetscInt VecSize,Mat A,PetscInt dnz,PetscInt* dnnz,PetscInt onz,PetscInt* onnz,Mat B,PetscInt dnzB,PetscInt* dnnzB,PetscInt onzB,PetscInt* onnzB,Vec vecb,Vec vece,PetscInt rank,PetscInt rank_hsl,PetscInt mpisize,char* tabfile, char *commsyntax,set_def *sets,dim_t nset, set_element *set_elems, array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar, elem_value **elem_vals2,offset_t ncofvar,offset_t ncofele,offset_t nvarele,closure_entry **closure_vals2,offset_t alltimeset,offset_t allregset,offset_t nintraeq,dim_t matsol,PetscInt Istart,PetscInt Iend,  offset_t nreg, offset_t ntime, offset_t *eq_addr, offset_t ndblock, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,dim_t laDi,dim_t laD,PetscReal cntl3,PetscReal cntl6,dim_t nesteddbbd,int localsize,PetscInt *ndbbddrank1,fortran_int* indata,dim_t mc66,fortran_int *ptx,struct timeval begintime,solve_real **xcf2){ //Johansen
  char tempfilenam[256],tempchar[256];
  PetscScalar value,*vals=NULL;
  PetscErrorCode ierr;
  PetscInt count,nz01,*ai=NULL,*aj=NULL;
  fortran_int k=0,m=1;
  offset_t i,j,lasize;
  solve_real *b1=NULL,*x0=NULL;
  PetscBool presol;/* NDBBD phase flag: presolve pass writes the interface files, solve pass consumes them */
  bool IsIni;
  FILE* tempvar;
  closure_entry *closure_vals;
  closure_vals=*closure_vals2;
  solve_real *xcf;
  xcf=*xcf2;

  clock_t timestr,timeend,timemulti;
  struct timeval endtime;
  struct timespec gettime_now,gettime_beg,gettime_end;
  long int start_time=0;
  double rep_time;
  size_t freadresult;
  elem_value *elem_vals;
  elem_vals=*elem_vals2;
  elem_value *elem_vals1=NULL;
  
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
    if(rank==0)logmsg(1,"Matrix preparation time %.2f s\n",(endtime.tv_sec - begintime.tv_sec)+((double)(endtime.tv_usec - begintime.tv_usec))/ 1000000);
    
    if(rank==rank_hsl) {
      jacobian_fill(tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,closure_vals,ndblock,alltimeset,allregset,eq_addr,counteq,nintraeq,A,B);
    }

    gettimeofday(&begintime, NULL);
    if(rank==0)logmsg(1,"Matrix calculation time %.2f s\n",(begintime.tv_sec - endtime.tv_sec)+((double)(begintime.tv_usec - endtime.tv_usec))/ 1000000);

    for (count=0; count<nvarele; count++) {
      if (closure_vals[count].is_exogenous) {
        value = closure_vals[count].shock_value;
        dnz=closure_vals[count].exo_index;
        VecSetValues(vece,1,&dnz,&value,INSERT_VALUES);
      }
    }
    MPI_Barrier(PETSC_COMM_WORLD);
    ierr = VecAssemblyBegin(vece);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(vece);
    CHKERRQ(ierr);
    if(rank==rank_hsl) {
      if(!inmemory){
      strcpy(tempfilenam,scratch_dir);
      strcat(tempfilenam,"_tempshock");
      sprintf(tempchar, "%d",rank);
      strcat(tempfilenam,tempchar);
      strcat(tempfilenam,".bin");
      if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
        printf("Error: cannot open %s for writing\n",tempfilenam);
      }
      fwrite(closure_vals, sizeof(closure_entry),nvarele, tempvar);
      fclose(tempvar);
      free(*closure_vals2);//
      *closure_vals2=NULL;//realloc (ha_cgeshock,1*sizeof(ha_cgeexovar));
      closure_vals=*closure_vals2;
      }
    }
    if(rank==rank_hsl) {
      if(!inmemory){
      strcpy(tempfilenam,scratch_dir);
      strcat(tempfilenam,"_tempvar");
      sprintf(tempchar, "%d",rank);
      strcat(tempfilenam,tempchar);
      strcat(tempfilenam,".bin");
      if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
        printf("Error: cannot open %s for writing\n",tempfilenam);
      }
      fwrite(elem_vals, sizeof(elem_value),ncofele+nvarele, tempvar);
      fclose(tempvar);
      free(*elem_vals2);//
      *elem_vals2=NULL;//realloc (ha_cofvar,1*sizeof(ha_cgevar));
      elem_vals=*elem_vals2;
      }
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

    gettimeofday(&endtime, NULL);
    if(rank==0)logmsg(1,"Matrix assembly time %.2f s\n",(endtime.tv_sec - begintime.tv_sec)+((double)(endtime.tv_usec - begintime.tv_usec))/ 1000000);
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

    if(matsol>=MM_DBBD) {
      gettimeofday(&begintime, NULL);
      clock_gettime(CLOCK_REALTIME, &gettime_beg);
      int *row_order= (int *) calloc (VecSize,sizeof(int));
      int *col_order= (int *) calloc (VecSize,sizeof(int));
      int *block_sizes= (int *) calloc (ndblock,sizeof(int));
      if(matsol==MM_DBBD) {
        dbbd_order(A,VecSize,mpisize,rank,Istart,Iend,nvarele,eq_addr,row_order,col_order,ndblock,block_sizes,countvarintra1,counteq,counteqnoadd,laA,cntl6);
        x0=realloc (x0,VecSize*sizeof(solve_real));
        dbbd_solve(A,vecb,x0,VecSize,mpisize,rank,Istart,Iend,row_order,col_order,ndblock,block_sizes,countvarintra1,counteq,counteqnoadd,laA,laD,cntl3);//,iter
      }
      if(matsol==MM_NDBBD) {
        presol=1;
        if(presol){
        ndbbd_order_presolve(A,VecSize,mpisize,rank,Istart,Iend,nreg,ntime,nvarele,eq_addr,row_order,col_order,ndblock,block_sizes,countvarintra1,counteq,counteqnoadd,laA,laDi,cntl6,ndbbddrank1,presol);
        ndbbd_presolve(A,vecb,x0,VecSize,mpisize,rank,Istart,Iend,row_order,col_order,ndblock,nreg,ntime,block_sizes,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,presol);//,iter
        }
        presol=0;
        ndbbd_order(A,VecSize,mpisize,rank,Istart,Iend,nreg,ntime,nvarele,eq_addr,row_order,col_order,ndblock,block_sizes,countvarintra1,counteq,counteqnoadd,laA,laDi,cntl6,ndbbddrank1,presol);
        x0=realloc (x0,VecSize*sizeof(solve_real));
        ndbbd_solve(A,vecb,x0,VecSize,mpisize,rank,Istart,Iend,row_order,col_order,ndblock,nreg,ntime,block_sizes,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,presol);//,iter
      }
      time(&timeend);
      gettimeofday(&endtime, NULL);
      clock_gettime(CLOCK_REALTIME, &gettime_end);
      rep_time = ((double)(gettime_end.tv_nsec-gettime_beg.tv_nsec))/1000000000.0;
      if(rank==0)logmsg(1,"Step time %.2f s\n",(endtime.tv_sec - begintime.tv_sec)+((double)(endtime.tv_usec - begintime.tv_usec))/ 1000000);
      if(rank==0)logmsg(1,"Step wall time %.2f s\n",rep_time);
      free(row_order);
      free(col_order);
      free(block_sizes);
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
        logmsg(2,"count %d nz %d\n",count,nz01);
      }
      indata[0]=count;//.nz

      if(matsol==MM_SBBD) {
        int *irn=(int *) calloc (count,sizeof(int));
        int *irn1=(int *) calloc (nz01,sizeof(int));
        int *jcn=(int *) calloc (count,sizeof(int));
        solve_real *values= (solve_real *) calloc (count,sizeof(solve_real));
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

        b1=realloc (b1,VecSize*sizeof(solve_real));
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
        x0=realloc (x0,VecSize*sizeof(solve_real));
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
        x0=realloc (x0,VecSize*sizeof(solve_real));
        lasize=ceil((laA/100.0)*count);
        int *irn=(int *) calloc (lasize,sizeof(int));
        int *irn1=(int *) calloc (nz01,sizeof(int));
        int *jcn=(int *) calloc (lasize,sizeof(int));
        solve_real *values= (solve_real *) calloc (lasize,sizeof(solve_real));

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
      if(rank==0)logmsg(1,"Step time %.2f s\n",(endtime.tv_sec - begintime.tv_sec)+((double)(endtime.tv_usec - begintime.tv_usec))/ 1000000);
    }
    if(rank==rank_hsl) {
      if(!inmemory){
      if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
        printf("Error: cannot open %s for reading\n",tempfilenam);
      }
      *elem_vals2=(elem_value*)realloc (*elem_vals2,(ncofele+nvarele)*sizeof(elem_value));
      freadresult=fread(*elem_vals2, sizeof(elem_value),ncofele+nvarele, tempvar);
      fclose(tempvar);
      remove(tempfilenam);
      }
      elem_vals=*elem_vals2;

      if(!inmemory){
      strcpy(tempfilenam,scratch_dir);
      strcat(tempfilenam,"_tempshock");
      sprintf(tempchar, "%d",rank);
      strcat(tempfilenam,tempchar);
      strcat(tempfilenam,".bin");
      if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
        printf("Error: cannot open %s for reading\n",tempfilenam);
      }
      *closure_vals2=(closure_entry*)realloc (*closure_vals2,(nvarele)*sizeof(closure_entry));
      freadresult=fread(*closure_vals2, sizeof(closure_entry),nvarele, tempvar);
      fclose(tempvar);
      remove(tempfilenam);
      }
      closure_vals=*closure_vals2;      
    }
    *xcf2=(solve_real*)realloc (*xcf2,nvarele*sizeof(solve_real));
    xcf=*xcf2;
    if(rank==rank_hsl) {
      elem_vals1=elem_vals+ncofele;
      for(i=0; i<nvar; i++) {
        if(vars[i].change_real) {
          for(j=vars[i].offset; j<vars[i].nelem+vars[i].offset; j++) {
            if(closure_vals[j].is_exogenous) {
              elem_vals1[j].initial=elem_vals1[j].value;
              elem_vals1[j].value+=closure_vals[j].shock_value;
              xcf[j]=closure_vals[j].shock_value;//varchange[j]
              elem_vals1[j].substep_base=closure_vals[j].shock_value;
            }
            else {
              elem_vals1[j].initial=elem_vals1[j].value;
              elem_vals1[j].value+=x0[closure_vals[j].exo_index];
              xcf[j]=x0[closure_vals[j].exo_index];//varchange[j]
              elem_vals1[j].substep_base=x0[closure_vals[j].exo_index];
            }
          }
        }
        else {
          for(j=vars[i].offset; j<vars[i].nelem+vars[i].offset; j++) {
            if(closure_vals[j].is_exogenous) {
              elem_vals1[j].initial=elem_vals1[j].value;
              elem_vals1[j].value+=closure_vals[j].shock_value*elem_vals1[j].initial/100;
              xcf[j]=closure_vals[j].shock_value;//varchange[j]
              elem_vals1[j].substep_base=closure_vals[j].shock_value;
            }
            else {
              elem_vals1[j].initial=elem_vals1[j].value;
              xcf[j]=x0[closure_vals[j].exo_index];//varchange[j]
              elem_vals1[j].value+=x0[closure_vals[j].exo_index]/100*elem_vals1[j].value;
              elem_vals1[j].substep_base=x0[closure_vals[j].exo_index];
            }
          }
        }
      }
      updates_apply(tabfile,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,0);
      strcpy(commsyntax,"formula");
      IsIni=false;
      formulas_execute(tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,IsIni);

    }
    elem_vals1=NULL;
    free(x0);
    return 1;
  }

bool solve_gragg(PetscBool nohsl,PetscInt VecSize,Mat* A1,PetscInt dnz,PetscInt* dnnz,PetscInt onz,PetscInt* onnz,Mat* B1,PetscInt dnzB,PetscInt* dnnzB,PetscInt onzB,PetscInt* onnzB,Vec* vecb1,Vec *vece1,PetscInt rank,PetscInt rank_hsl,PetscInt mpisize,char* tabfile, char *commsyntax,set_def *sets,dim_t nset, set_element *set_elems, array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar, elem_value **elem_vals2,offset_t ncofvar,offset_t ncofele,offset_t nvarele,closure_entry **closure_vals2,offset_t alltimeset,offset_t allregset,offset_t nintraeq,dim_t matsol,PetscInt Istart,PetscInt Iend,  offset_t nreg, offset_t ntime, offset_t *eq_addr, offset_t ndblock, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,dim_t laDi,dim_t laD,PetscReal cntl3,PetscReal cntl6,dim_t nesteddbbd,int localsize,PetscInt *ndbbddrank1,fortran_int* indata,dim_t mc66,fortran_int *ptx,struct timeval begintime,dim_t subints,MPI_Fint fcomm,solve_real **xcf2){ /* Gragg method (smoothed modified midpoint), Pearson 1991 eq. 6.1 / Alg. 7.1.2 */
  char tempfilenam[256],tempchar[256],solchar[255];
  PetscScalar value,*vals;
  PetscErrorCode ierr;
  PetscInt count,nz01,*ai,*aj;
  fortran_int k=0,m=1;
  fortran_int tindx1;//,tindx2;
  solve_real temp1,temp2;
  offset_t i,j,lasize;
  dim_t subindx;
  solve_real *b1=NULL;
  solve_real *x1=NULL;
  solve_real *xcf;
  xcf=*xcf2;
  PetscBool presol;/* NDBBD phase flag: presolve pass writes the interface files, solve pass consumes them */
  bool IsIni;
  FILE* tempvar;
  PetscLogDouble time1,time0;
  clock_t timestr,timeend,timemulti;
  struct timeval endtime;
  struct timespec gettime_now,gettime_beg,gettime_end;
  long int start_time=0;
  double rep_time;
  size_t freadresult;
  elem_value *elem_vals;
  elem_vals=*elem_vals2;
  elem_value *elem_vals1;
  closure_entry *closure_vals;
  closure_vals=*closure_vals2;
  Vec vece,vecb;
  Mat A,B;
  A=*A1;
  B=*B1;
  vece=*vece1,
  vecb=*vecb1;
  int stepcount;
  int nsteps=3;
  int sol;
  solve_real vpercents=1.0,perprecis=0;
  FILE* solution;
  int maxsol=3;
              offset_t *counteqs= (offset_t *) calloc (ndblock+1,sizeof(offset_t));
              offset_t *counteqnoadds= (offset_t *) calloc (ndblock,sizeof(offset_t));
              offset_t *countvarintra1s= (offset_t *) calloc (ndblock+1,sizeof(offset_t));
              memcpy(counteqs,counteq,(ndblock+1)*sizeof(offset_t));
              memcpy(counteqnoadds,counteqnoadd,(ndblock)*sizeof(offset_t));
              memcpy(countvarintra1s,countvarintra1,(ndblock+1)*sizeof(offset_t));
    gettimeofday(&begintime, NULL);
    solve_real *xc0= (solve_real *) calloc (1,sizeof(solve_real));
    solve_real *xc12= (solve_real *) calloc (1,sizeof(solve_real));
    solve_real *xc24= (solve_real *) calloc (1,sizeof(solve_real));
    int *xc124= (int *) calloc (1,sizeof(int));
    solve_real *clag1= (solve_real *) calloc (nvarele,sizeof(solve_real));
    solve_real *varchange= (solve_real *) calloc (nvarele,sizeof(solve_real));
    for(subindx=0; subindx<subints; subindx++) {
      for(sol=0; sol<maxsol; sol++) {
        if(sol==0)nsteps=steps1;
        if(sol==1) nsteps=(int)steps1*step_ratio2;
        if(sol==2) nsteps=(int)steps1*step_ratio3;
        vpercents=(solve_real)100/nsteps;
        for(stepcount=0; stepcount<nsteps; stepcount++) {
          logmsg(2,"rank %d subint %d sol %d stepcount %d nsteps %d\n",rank,subindx,sol,stepcount,nsteps);
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
                elem_vals[i].initial=elem_vals[i].value;
              }
            else for(i=0; i<ncofele; i++) {
                elem_vals[i].value=elem_vals[i].initial;
              }
            elem_vals1=elem_vals+ncofele;
            for(i=0; i<nvar; i++) {
              if(vars[i].change_real) {
                for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
                  if(closure_vals[tindx1].is_exogenous) {
                    if(sol==0) {
                      elem_vals1[tindx1].initial=elem_vals1[tindx1].value;
                    }
                    else {
                      elem_vals1[tindx1].value=elem_vals1[tindx1].initial;
                    }
                    elem_vals1[tindx1].substep_base=closure_vals[tindx1].shock_value/nsteps;
                    VecSetValue(vece,closure_vals[tindx1].exo_index,elem_vals1[tindx1].substep_base,INSERT_VALUES);
                  }
                  else {
                    if(sol==0) {
                      elem_vals1[tindx1].initial=elem_vals1[tindx1].value;
                    }
                    else {
                      elem_vals1[tindx1].value=elem_vals1[tindx1].initial;
                    }
                  }
                }
              }
              else {
                for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
                  if(closure_vals[tindx1].is_exogenous) {
                    if(sol==0) {
                      elem_vals1[tindx1].initial=elem_vals1[tindx1].value;
                    }
                    else {
                      elem_vals1[tindx1].value=elem_vals1[tindx1].initial;
                    }
                    temp2=closure_vals[tindx1].shock_value;//subints;
                    elem_vals1[tindx1].substep_base=(100+(subindx+1)*temp2)/(100+subindx*temp2)-1;//ha_cgeshock[ha_var[i].begadd+j].ShockVal/nsteps;//(exp(log(1+ha_cgeshock[ha_var[i].begadd+j].ShockVal/100)/nsteps)-1)*100;
                    elem_vals1[tindx1].substep_base*=vpercents;//nsteps*100;
                    VecSetValue(vece,closure_vals[tindx1].exo_index,elem_vals1[tindx1].substep_base,INSERT_VALUES);
                  }
                  else {
                    if(sol==0) {
                      elem_vals1[tindx1].initial=elem_vals1[tindx1].value;
                    }
                    else {
                      elem_vals1[tindx1].value=elem_vals1[tindx1].initial;
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
          if(rank==rank_hsl) {

            if(!inmemory){
            strcpy(tempfilenam,scratch_dir);
            strcat(tempfilenam,"_tempclag1");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
              printf("Error: cannot open %s for writing\n",tempfilenam);
            }
            fwrite(clag1, sizeof(solve_real),nvarele, tempvar);
            fclose(tempvar);
            free(clag1);
            clag1=NULL;
            }

            if(!inmemory){
            strcpy(tempfilenam,scratch_dir);
            strcat(tempfilenam,"_tempvarchange");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
              printf("Error: cannot open %s for writing\n",tempfilenam);
            }
            fwrite(varchange, sizeof(solve_real),nvarele, tempvar);
            fclose(tempvar);
            free(varchange);
            varchange=NULL;
            }
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
          if(rank==rank_hsl) {
            jacobian_fill(tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,closure_vals,ndblock,alltimeset,allregset,eq_addr,counteq,nintraeq,A,B);
          }
          if(rank==rank_hsl) {
            if(!inmemory){
            strcpy(tempfilenam,scratch_dir);
            strcat(tempfilenam,"_tempvar");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
              printf("Error: cannot open %s for writing\n",tempfilenam);
            }
            fwrite(elem_vals, sizeof(elem_value),ncofele+nvarele, tempvar);
            fclose(tempvar);
            free(*elem_vals2);
            *elem_vals2=NULL;
            elem_vals=*elem_vals2;
            }

            if(!inmemory){
            strcpy(tempfilenam,scratch_dir);
            strcat(tempfilenam,"_tempshock");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
              printf("Error: cannot open %s for writing\n",tempfilenam);
            }
            fwrite(closure_vals, sizeof(closure_entry),nvarele, tempvar);
            fclose(tempvar);
            free(*closure_vals2);
            *closure_vals2=NULL;
            closure_vals=*closure_vals2;
            }
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
          if(matsol>=MM_DBBD) {
            int *row_order= (int *) calloc (VecSize,sizeof(int));
            int *col_order= (int *) calloc (VecSize,sizeof(int));
            int *block_sizes= (int *) calloc (ndblock,sizeof(int));
            time(&timestr);

            if(matsol==MM_DBBD) {
              dbbd_order(A,VecSize,mpisize,rank,Istart,Iend,nvarele,eq_addr,row_order,col_order,ndblock,block_sizes,countvarintra1,counteq,counteqnoadd,laA,cntl6);
              x1=realloc (x1,VecSize*sizeof(solve_real));
              dbbd_solve(A,vecb,x1,VecSize,mpisize,rank,Istart,Iend,row_order,col_order,ndblock,block_sizes,countvarintra1,counteq,counteqnoadd,laA,laD,cntl3);//,iter
            }

            if(matsol==MM_NDBBD) {
              presol=1;
              memcpy(counteq,counteqs,(ndblock+1)*sizeof(offset_t));
              memcpy(counteqnoadd,counteqnoadds,(ndblock)*sizeof(offset_t));
              memcpy(countvarintra1,countvarintra1s,(ndblock+1)*sizeof(offset_t));
              ndbbd_order_presolve(A,VecSize,mpisize,rank,Istart,Iend,nreg,ntime,nvarele,eq_addr,row_order,col_order,ndblock,block_sizes,countvarintra1,counteq,counteqnoadd,laA,laDi,cntl6,ndbbddrank1,presol);
              ndbbd_presolve(A,vecb,x1,VecSize,mpisize,rank,Istart,Iend,row_order,col_order,ndblock,nreg,ntime,block_sizes,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,presol);//,iter
              presol=0;
              ndbbd_order(A,VecSize,mpisize,rank,Istart,Iend,nreg,ntime,nvarele,eq_addr,row_order,col_order,ndblock,block_sizes,countvarintra1,counteq,counteqnoadd,laA,laDi,cntl6,ndbbddrank1,presol);
              x1=realloc (x1,VecSize*sizeof(solve_real));
              ndbbd_solve(A,vecb,x1,VecSize,mpisize,rank,Istart,Iend,row_order,col_order,ndblock,nreg,ntime,block_sizes,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,presol);//,iter
            }

            time(&timeend);
            MPI_Barrier(PETSC_COMM_WORLD);
            ierr = PetscGetCPUTime(&time1);
            if(verbosity>=1){ierr = PetscPrintf(PETSC_COMM_WORLD,"One step solution %f\n",time1-time0);}
            if(rank==0)logmsg(1,"Step time %.2f s\n",difftime(timeend,timestr));
            free(row_order);
            free(col_order);
            free(block_sizes);
            MPI_Barrier(PETSC_COMM_WORLD);
          }
          else {
            if(matsol==MM_SBBD) {
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
              solve_real *values= (solve_real *) calloc (count,sizeof(solve_real));
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
              b1=realloc (b1,VecSize*sizeof(solve_real));
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
              x1=realloc (x1,VecSize*sizeof(solve_real));
              ierr = PetscGetCPUTime(&time1);
              CHKERRQ(ierr);
              if(verbosity>=1){ierr = PetscPrintf(PETSC_COMM_WORLD,"Prepare time %f\n",time1-time0);}
              CHKERRQ(ierr);
              ierr = PetscGetCPUTime(&time0);
              CHKERRQ(ierr);
              if(mc66!=0)spec48_single_(ptx,irn,jcn,b1,values,x1,neleperrow,ai1,&fcomm);
              free(irn);
              if(mc66==0)spec48_nomc66_(ptx,jcn,b1,values,x1,neleperrow,&fcomm,counteq,countvarintra1);
              ierr = PetscGetCPUTime(&time1);
              CHKERRQ(ierr);
              if(verbosity>=1){ierr = PetscPrintf(PETSC_COMM_WORLD,"LU time %f\n",time1-time0);}
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
              solve_real *values= (solve_real *) calloc (lasize,sizeof(solve_real));
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
              x1=realloc (x1,VecSize*sizeof(solve_real));
              ierr = PetscGetCPUTime(&time1);
              CHKERRQ(ierr);
              if(verbosity>=1){ierr = PetscPrintf(PETSC_COMM_WORLD,"Prepare time %f\n",time1-time0);}
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
              if(verbosity>=1){ierr = PetscPrintf(PETSC_COMM_WORLD,"LU time %f\n",time1-time0);}
              CHKERRQ(ierr);
              ierr = PetscGetCPUTime(&time0);
              CHKERRQ(ierr);
              free(irn);
              free(jcn);
              free(values);
            }
          }
          if(rank==rank_hsl) {
            if(!inmemory){
            if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
              printf("Error: cannot open %s for reading\n",tempfilenam);
            }
            *closure_vals2=(closure_entry*)realloc (*closure_vals2,(nvarele)*sizeof(closure_entry));
            freadresult=fread(*closure_vals2, sizeof(closure_entry),nvarele, tempvar);
            fclose(tempvar);
            remove(tempfilenam);
            }
            closure_vals=*closure_vals2;

            if(!inmemory){
            strcpy(tempfilenam,scratch_dir);
            strcat(tempfilenam,"_tempvar");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
              printf("Error: cannot open %s for reading\n",tempfilenam);
            }
            *elem_vals2=(elem_value*)realloc (*elem_vals2,(ncofele+nvarele)*sizeof(elem_value));
            freadresult=fread(*elem_vals2, sizeof(elem_value),ncofele+nvarele, tempvar);
            fclose(tempvar);
            remove(tempfilenam);
            }
            elem_vals=*elem_vals2;

            if(!inmemory){
            strcpy(tempfilenam,scratch_dir);
            strcat(tempfilenam,"_tempclag1");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
              printf("Error: cannot open %s for reading\n",tempfilenam);
            }
            clag1=realloc (clag1,(nvarele)*sizeof(solve_real));
            freadresult=fread(clag1, sizeof(solve_real),nvarele, tempvar);
            fclose(tempvar);
            remove(tempfilenam);
            }

            if(!inmemory){
            strcpy(tempfilenam,scratch_dir);
            strcat(tempfilenam,"_tempvarchange");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
              printf("Error: cannot open %s for reading\n",tempfilenam);
            }
            varchange=realloc (varchange,(nvarele)*sizeof(solve_real));
            freadresult=fread(varchange, sizeof(solve_real),nvarele, tempvar);
            fclose(tempvar);
            remove(tempfilenam);
            }

          }
          logmsg(2,"sol %d stepcount %d\n\n",sol,stepcount);
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
          elem_vals1=elem_vals+ncofele;
          if(stepcount==0) {
            for(i=0; i<nvar; i++) {
              if(vars[i].change_real) {
                for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
                  if(closure_vals[tindx1].is_exogenous) {
                    elem_vals1[tindx1].value+=elem_vals1[tindx1].substep_base;
                    varchange[tindx1]=elem_vals1[tindx1].substep_base;
                    VecSetValue(vece,closure_vals[tindx1].exo_index,elem_vals1[tindx1].substep_base,INSERT_VALUES);
                  }
                  else {
                    varchange[tindx1]=x1[closure_vals[tindx1].exo_index];
                    elem_vals1[tindx1].value+=x1[closure_vals[tindx1].exo_index];
                    elem_vals1[tindx1].substep_base=x1[closure_vals[tindx1].exo_index];
                    clag1[tindx1]=0;
                  }
                }
              }
              else {
                for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
                  if(closure_vals[tindx1].is_exogenous) {
                    varchange[tindx1]=elem_vals1[tindx1].substep_base;
                    elem_vals1[tindx1].value*=(1+elem_vals1[tindx1].substep_base/100);
                    VecSetValue(vece,closure_vals[tindx1].exo_index,elem_vals1[tindx1].substep_base/(1+elem_vals1[tindx1].substep_base/100),INSERT_VALUES);
                  }
                  else {
                    varchange[tindx1]=x1[closure_vals[tindx1].exo_index];
                    elem_vals1[tindx1].substep_base=x1[closure_vals[tindx1].exo_index];
                    elem_vals1[tindx1].value*=(1+elem_vals1[tindx1].substep_base/100);
                    clag1[tindx1]=0;
                  }
                }
              }
            }
          }
          else {
            for(i=0; i<nvar; i++) {
              if(vars[i].change_real) {
                for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
                  if(closure_vals[tindx1].is_exogenous) {
                    elem_vals1[tindx1].value+=elem_vals1[tindx1].substep_base;
                    varchange[tindx1]+=elem_vals1[tindx1].substep_base;
                    VecSetValue(vece,closure_vals[tindx1].exo_index,elem_vals1[tindx1].substep_base,INSERT_VALUES);
                  }
                  else {
                    temp1=elem_vals1[tindx1].value;//change;
                    varchange[tindx1]=clag1[tindx1]+2*x1[closure_vals[tindx1].exo_index];//+=x1[ha_cgeshock[ha_var[i].begadd+j].ExoIndx];//
                    elem_vals1[tindx1].substep_base=x1[closure_vals[tindx1].exo_index];//ha_cofvar[ncofele+ha_var[i].begadd+j].varchange-temp1;
                    elem_vals1[tindx1].value=clag1[tindx1]+2*x1[closure_vals[tindx1].exo_index];//ha_cofvar[ncofele+ha_var[i].begadd+j].varchange-temp1;//ha_cofvar[ncofele+ha_var[i].begadd+j].csolpupd;
                    clag1[tindx1]=temp1;
                  }
                }
              }
              else {
                for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
                  if(closure_vals[tindx1].is_exogenous) {
                    temp2=closure_vals[tindx1].shock_value;//subints;
                    temp1=(100+(subindx+1)*temp2)/(100+subindx*temp2)-1;
                    temp1*=vpercents;
                    elem_vals1[tindx1].substep_base=temp1/(1+varchange[tindx1]/100);
                    varchange[tindx1]+=temp1;//*(1+ha_cofvar[ncofele+ha_var[i].begadd+j].varchange/100)
                    elem_vals1[tindx1].value=(1+varchange[tindx1]/100)*elem_vals1[tindx1].initial;
                    VecSetValue(vece,closure_vals[tindx1].exo_index,temp1/(1+varchange[tindx1]/100),INSERT_VALUES);
                  }
                  else {
                    temp1=varchange[tindx1];
                    varchange[tindx1]=clag1[tindx1]+2*x1[closure_vals[tindx1].exo_index]*(100+temp1)/100;//+=x1[ha_cgeshock[ha_var[i].begadd+j].ExoIndx]*(1+temp1/100);//
                    elem_vals1[tindx1].substep_base=x1[closure_vals[tindx1].exo_index];//(ha_cofvar[ncofele+ha_var[i].begadd+j].varchange-temp1)/(1+temp1/100);
                    elem_vals1[tindx1].value=(100+varchange[tindx1])/100*elem_vals1[tindx1].initial;
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
              updates_apply(tabfile,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,0);
            }
            else {
              updates_apply(tabfile,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,1);
            }
            strcpy(commsyntax,"formula");
            IsIni=false;
            formulas_execute(tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,IsIni);
          }
          MPI_Barrier(PETSC_COMM_WORLD);
          ierr = PetscGetCPUTime(&time1);
          CHKERRQ(ierr);
          if(verbosity>=1){ierr = PetscPrintf(PETSC_COMM_WORLD,"Update time %f\n",time1-time0);}
          CHKERRQ(ierr);
          ierr = PetscGetCPUTime(&time0);
          CHKERRQ(ierr);
        }

        strcpy(commsyntax,"equation");
        if(rank==rank_hsl) {


          if(!inmemory){
          strcpy(tempfilenam,scratch_dir);
          strcat(tempfilenam,"_tempclag1");
          sprintf(tempchar, "%d",rank);
          strcat(tempfilenam,tempchar);
          strcat(tempfilenam,".bin");
          if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
            printf("Error: cannot open %s for writing\n",tempfilenam);
          }
          fwrite(clag1, sizeof(solve_real),nvarele, tempvar);
          fclose(tempvar);
          free(clag1);
          clag1=NULL;
          }

          if(!inmemory){
          strcpy(tempfilenam,scratch_dir);
          strcat(tempfilenam,"_tempvarchange");
          sprintf(tempchar, "%d",rank);
          strcat(tempfilenam,tempchar);
          strcat(tempfilenam,".bin");
          if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
            printf("Error: cannot open %s for writing\n",tempfilenam);
          }
          fwrite(varchange, sizeof(solve_real),nvarele, tempvar);
          fclose(tempvar);
          free(varchange);
          varchange=NULL;
          }

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
          jacobian_fill(tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,closure_vals,ndblock,alltimeset,allregset,eq_addr,counteq,nintraeq,A,B);
        }

        if(rank==rank_hsl) {
          if(!inmemory){
          strcpy(tempfilenam,scratch_dir);
          strcat(tempfilenam,"_tempvar");
          sprintf(tempchar, "%d",rank);
          strcat(tempfilenam,tempchar);
          strcat(tempfilenam,".bin");
          if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
            printf("Error: cannot open %s for writing\n",tempfilenam);
          }
          fwrite(elem_vals, sizeof(elem_value),ncofele+nvarele, tempvar);
          fclose(tempvar);
          free(*elem_vals2);
          *elem_vals2=NULL;
          elem_vals=*elem_vals2;
          }

          if(!inmemory){
          strcpy(tempfilenam,scratch_dir);
          strcat(tempfilenam,"_tempshock");
          sprintf(tempchar, "%d",rank);
          strcat(tempfilenam,tempchar);
          strcat(tempfilenam,".bin");
          if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
            printf("Error: cannot open %s for writing\n",tempfilenam);
          }
          fwrite(closure_vals, sizeof(closure_entry),nvarele, tempvar);
          fclose(tempvar);
          free(*closure_vals2);
          *closure_vals2=NULL;
          closure_vals=*closure_vals2;
          }
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
        if(matsol>=MM_DBBD) {
          int *row_order= (int *) calloc (VecSize,sizeof(int));
          int *col_order= (int *) calloc (VecSize,sizeof(int));
          int *block_sizes= (int *) calloc (ndblock,sizeof(int));
          time(&timestr);

          if(matsol==MM_DBBD) {
            dbbd_order(A,VecSize,mpisize,rank,Istart,Iend,nvarele,eq_addr,row_order,col_order,ndblock,block_sizes,countvarintra1,counteq,counteqnoadd,laA,cntl6);
            x1=realloc (x1,VecSize*sizeof(solve_real));
            dbbd_solve(A,vecb,x1,VecSize,mpisize,rank,Istart,Iend,row_order,col_order,ndblock,block_sizes,countvarintra1,counteq,counteqnoadd,laA,laD,cntl3);//,iter
          }

          if(matsol==MM_NDBBD) {
            presol=1;
            memcpy(counteq,counteqs,(ndblock+1)*sizeof(offset_t));
            memcpy(counteqnoadd,counteqnoadds,(ndblock)*sizeof(offset_t));
            memcpy(countvarintra1,countvarintra1s,(ndblock+1)*sizeof(offset_t));
            ndbbd_order_presolve(A,VecSize,mpisize,rank,Istart,Iend,nreg,ntime,nvarele,eq_addr,row_order,col_order,ndblock,block_sizes,countvarintra1,counteq,counteqnoadd,laA,laDi,cntl6,ndbbddrank1,presol);
            ndbbd_presolve(A,vecb,x1,VecSize,mpisize,rank,Istart,Iend,row_order,col_order,ndblock,nreg,ntime,block_sizes,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,presol);//,iter
            presol=0;
            ndbbd_order(A,VecSize,mpisize,rank,Istart,Iend,nreg,ntime,nvarele,eq_addr,row_order,col_order,ndblock,block_sizes,countvarintra1,counteq,counteqnoadd,laA,laDi,cntl6,ndbbddrank1,presol);
            x1=realloc (x1,VecSize*sizeof(solve_real));
            ndbbd_solve(A,vecb,x1,VecSize,mpisize,rank,Istart,Iend,row_order,col_order,ndblock,nreg,ntime,block_sizes,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,presol);//,iter
          }

          time(&timeend);
          MPI_Barrier(PETSC_COMM_WORLD);
          ierr = PetscGetCPUTime(&time1);
          if(verbosity>=1){ierr = PetscPrintf(PETSC_COMM_WORLD,"One step solution %f\n",time1-time0);}
          if(rank==0)logmsg(1,"Step time %.2f s\n",difftime(timeend,timestr));
          free(row_order);
          free(col_order);
          free(block_sizes);
          MPI_Barrier(PETSC_COMM_WORLD);
        }
        else {
          if(matsol==MM_SBBD) {
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
            solve_real *values= (solve_real *) calloc (count,sizeof(solve_real));
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
            b1=realloc (b1,VecSize*sizeof(solve_real));
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
            x1=realloc (x1,VecSize*sizeof(solve_real));
            ierr = PetscGetCPUTime(&time1);
            CHKERRQ(ierr);
            if(verbosity>=1){ierr = PetscPrintf(PETSC_COMM_WORLD,"Prepare time %f\n",time1-time0);}
            CHKERRQ(ierr);
            ierr = PetscGetCPUTime(&time0);
            CHKERRQ(ierr);
            if(mc66!=0)spec48_single_(ptx,irn,jcn,b1,values,x1,neleperrow,ai1,&fcomm);
            free(irn);
            if(mc66==0)spec48_nomc66_(ptx,jcn,b1,values,x1,neleperrow,&fcomm,counteq,countvarintra1);
            ierr = PetscGetCPUTime(&time1);
            CHKERRQ(ierr);
            if(verbosity>=1){ierr = PetscPrintf(PETSC_COMM_WORLD,"LU time %f\n",time1-time0);}
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
            solve_real *values= (solve_real *) calloc (lasize,sizeof(solve_real));
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
            x1=realloc (x1,VecSize*sizeof(solve_real));
            ierr = PetscGetCPUTime(&time1);
            CHKERRQ(ierr);
            if(verbosity>=1){ierr = PetscPrintf(PETSC_COMM_WORLD,"Prepare time %f\n",time1-time0);}
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
            if(verbosity>=1){ierr = PetscPrintf(PETSC_COMM_WORLD,"LU time %f\n",time1-time0);}
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
          if(!inmemory){
          if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
            printf("Error: cannot open %s for reading\n",tempfilenam);
          }
          *closure_vals2=(closure_entry*)realloc (*closure_vals2,(nvarele)*sizeof(closure_entry));
          freadresult=fread(*closure_vals2, sizeof(closure_entry),nvarele, tempvar);
          fclose(tempvar);
          remove(tempfilenam);
          }
          closure_vals=*closure_vals2;

          if(!inmemory){
          strcpy(tempfilenam,scratch_dir);
          strcat(tempfilenam,"_tempvar");
          sprintf(tempchar, "%d",rank);
          strcat(tempfilenam,tempchar);
          strcat(tempfilenam,".bin");
          if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
            printf("Error: cannot open %s for reading\n",tempfilenam);
          }
          *elem_vals2=(elem_value*)realloc (*elem_vals2,(ncofele+nvarele)*sizeof(elem_value));
          freadresult=fread(*elem_vals2, sizeof(elem_value),ncofele+nvarele, tempvar);
          fclose(tempvar);
          remove(tempfilenam);
          }
          elem_vals=*elem_vals2;

          if(!inmemory){
          strcpy(tempfilenam,scratch_dir);
          strcat(tempfilenam,"_tempclag1");
          sprintf(tempchar, "%d",rank);
          strcat(tempfilenam,tempchar);
          strcat(tempfilenam,".bin");
          if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
            printf("Error: cannot open %s for reading\n",tempfilenam);
          }
          clag1=realloc (clag1,(nvarele)*sizeof(solve_real));
          freadresult=fread(clag1, sizeof(solve_real),nvarele, tempvar);
          fclose(tempvar);
          remove(tempfilenam);
          }

          if(!inmemory){
          strcpy(tempfilenam,scratch_dir);
          strcat(tempfilenam,"_tempvarchange");
          sprintf(tempchar, "%d",rank);
          strcat(tempfilenam,tempchar);
          strcat(tempfilenam,".bin");
          if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
            printf("Error: cannot open %s for reading\n",tempfilenam);
          }
          varchange=realloc (varchange,(nvarele)*sizeof(solve_real));
          freadresult=fread(varchange, sizeof(solve_real),nvarele, tempvar);
          fclose(tempvar);
          remove(tempfilenam);
          }

        }
        elem_vals1=elem_vals+ncofele;
        for(i=0; i<nvar; i++) {
          if(vars[i].change_real) {
            for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
              if(closure_vals[tindx1].is_exogenous) {
                elem_vals1[tindx1].value=0;
              }
              else {
                varchange[tindx1]=0.5*(varchange[tindx1]+clag1[tindx1]+x1[closure_vals[tindx1].exo_index]);
                elem_vals1[tindx1].value=0;//ha_cofvar[tindx2].var0+varchange[tindx1];//no distortion between steps
                clag1[tindx1]=0;
              }
            }
          }
          else {
            for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
              if(closure_vals[tindx1].is_exogenous) {
                elem_vals1[tindx1].value=0;
              }
              else {
                varchange[tindx1]=0.5*(varchange[tindx1]+clag1[tindx1]+x1[closure_vals[tindx1].exo_index]*(1+varchange[tindx1]/100));
                elem_vals1[tindx1].value=0;//ha_cofvar[tindx2].varval*varchange[tindx1]/100;
                clag1[tindx1]=0;
              }
            }
          }
        }
        if(rank==rank_hsl) {
          if(subindx!=0||sol!=0) {
            if(!inmemory){
            strcpy(tempfilenam,scratch_dir);
            strcat(tempfilenam,"_tempxcf");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
              printf("Error: cannot open %s for reading\n",tempfilenam);
            }
            *xcf2=(solve_real*)realloc (*xcf2,(nvarele)*sizeof(solve_real));
            xcf=*xcf2;
            freadresult=fread(xcf, sizeof(solve_real),nvarele, tempvar);
            fclose(tempvar);
            remove(tempfilenam);
            }

            if(!inmemory){
            strcpy(tempfilenam,scratch_dir);
            strcat(tempfilenam,"_tempxc12");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
              printf("Error: cannot open %s for reading\n",tempfilenam);
            }
            xc12=realloc (xc12,(nvarele)*sizeof(solve_real));
            freadresult=fread(xc12, sizeof(solve_real),nvarele, tempvar);
            fclose(tempvar);
            remove(tempfilenam);
            }
            if(!inmemory){
            strcpy(tempfilenam,scratch_dir);
            strcat(tempfilenam,"_tempxc24");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
              printf("Error: cannot open %s for reading\n",tempfilenam);
            }
            xc24=realloc (xc24,(nvarele)*sizeof(solve_real));
            freadresult=fread(xc24, sizeof(solve_real),nvarele, tempvar);
            fclose(tempvar);
            remove(tempfilenam);
            }
            
            xc0=realloc (xc0,(nvarele)*sizeof(solve_real));
            if(subindx>0&&sol>0){
            if(!inmemory){
            strcpy(tempfilenam,scratch_dir);
            strcat(tempfilenam,"_tempxcO");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
              printf("Error: cannot open %s for reading\n",tempfilenam);
            }
            freadresult=fread(xc0, sizeof(solve_real),nvarele, tempvar);
            fclose(tempvar);
            if(subindx==subints-1&&sol==maxsol-1)remove(tempfilenam);
            }
            }

          }
          if(subindx==0&&sol==0) {
            xc0=realloc (xc0,nvarele*sizeof(solve_real));
            *xcf2=(solve_real*)realloc (*xcf2,nvarele*sizeof(solve_real));
            xcf=*xcf2;
            for(i=0; i<nvarele; i++)xcf[i]=0;
          }
          if(sol==0)xc12=realloc (xc12,nvarele*sizeof(solve_real));
          if(sol==0)xc24=realloc (xc24,nvarele*sizeof(solve_real));
          if(subindx>0) {
            if(sol==0)for(i=0; i<nvarele; i++) xc0[i]=1+xcf[i]/100;//if(i==1287)printf("sol!!!!!!!!!!!!!!!!!! %d step %d xc %lf xc0 %lf k %d\n",sol,stepcount,1.0+xc[k]/100,xc0[i],i);}
            if(sol==0) {
              extrap_w1=1.0/(step_ratio2*step_ratio2-1.0);
              extrap_w2=1.0/(1-step_ratio2*step_ratio2)/(1.0-step_ratio3*step_ratio3);
              for(i=0; i<nvar; i++) {
                if(vars[i].change_real) {
                  for(k=vars[i].offset; k<vars[i].nelem+vars[i].offset; k++) {
                    xc12[k]=xcf[k]-varchange[k]*extrap_w1;
                    xc24[k]=xcf[k];
                    xcf[k]+=varchange[k]*extrap_w2;
                  }
                }
                else {
                  for(k=vars[i].offset; k<vars[i].nelem+vars[i].offset; k++) {
                    xc24[k]=xcf[k];
                    xc12[k]=xcf[k]-varchange[k]*xc0[k]*extrap_w1;
                    xcf[k]+=varchange[k]*xc0[k]*extrap_w2;//(100+xc0[k])*(100+varchange[k]/45)/100-100;//varchange[k]/45;
                  }
                }
              }
            }
            if(sol==1) {
              extrap_w1=step_ratio2*step_ratio2/(step_ratio2*step_ratio2-1.0);
              extrap_w2=step_ratio2*step_ratio2/(step_ratio3*step_ratio3-step_ratio2*step_ratio2);
              extrap_w3=step_ratio2*step_ratio2*step_ratio2*step_ratio2/(step_ratio2*step_ratio2-step_ratio3*step_ratio3)/(1.0-step_ratio2*step_ratio2);
              for(i=0; i<nvar; i++) {
                if(vars[i].change_real) {
                  for(k=vars[i].offset; k<vars[i].nelem+vars[i].offset; k++) {
                    xc24[k]-=varchange[k]*extrap_w2;
                    xc12[k]+=varchange[k]*extrap_w1;
                    xcf[k]-=varchange[k]*extrap_w3;
                  }
                }
                else {
                  for(k=vars[i].offset; k<vars[i].nelem+vars[i].offset; k++) {
                    xc24[k]-=varchange[k]*xc0[k]*extrap_w2;
                    xc12[k]+=varchange[k]*xc0[k]*extrap_w1;
                    xcf[k]-=varchange[k]*xc0[k]*extrap_w3;//(100+xc0[k])*(100-20*varchange[k]/45)/100-100;
                  }
                }
              }
            }
            if(sol==2) {
              extrap_w2=step_ratio3*step_ratio3/(step_ratio3*step_ratio3-step_ratio2*step_ratio2);
              extrap_w3=step_ratio3*step_ratio3*step_ratio3*step_ratio3/(step_ratio2*step_ratio2-step_ratio3*step_ratio3)/(1.0-step_ratio3*step_ratio3);
              for(i=0; i<nvar; i++) {
                if(vars[i].change_real) {
                  for(k=vars[i].offset; k<vars[i].nelem+vars[i].offset; k++) {
                    xc24[k]+=varchange[k]*extrap_w2;
                    xcf[k]+=varchange[k]*extrap_w3;
                  }
                }
                else {
                  for(k=vars[i].offset; k<vars[i].nelem+vars[i].offset; k++) {
                    xc24[k]+=varchange[k]*xc0[k]*extrap_w2;
                    xcf[k]+=varchange[k]*xc0[k]*extrap_w3;//(100+xc0[k])*(100+64*varchange[k]/45)/100-100;
                  }
                }
              }
            }
          }
          else {
            if(sol==0){
              extrap_w1=1.0/(step_ratio2*step_ratio2-1);
              extrap_w2=1.0/(1-step_ratio2*step_ratio2)/(1-step_ratio3*step_ratio3);
              for(i=0; i<nvarele; i++) {
                xc12[i]=-varchange[i]*extrap_w1;
                xcf[i]+=varchange[i]*extrap_w2;
              }
            }
            if(sol==1) {
              extrap_w1=step_ratio2*step_ratio2/(step_ratio2*step_ratio2-1.0);
              extrap_w2=step_ratio2*step_ratio2/(step_ratio3*step_ratio3-step_ratio2*step_ratio2);
              extrap_w3=step_ratio2*step_ratio2*step_ratio2*step_ratio2/(step_ratio2*step_ratio2-step_ratio3*step_ratio3)/(1.0-step_ratio2*step_ratio2);
              for(i=0; i<nvarele; i++) {
                xc24[i]=-varchange[i]*extrap_w2;
                xc12[i]+=varchange[i]*extrap_w1;
                xcf[i]-=varchange[i]*extrap_w3;
              }
            }
            if(sol==2) {
              extrap_w2=step_ratio3*step_ratio3/(step_ratio3*step_ratio3-step_ratio2*step_ratio2);
              extrap_w3=step_ratio3*step_ratio3*step_ratio3*step_ratio3/(step_ratio2*step_ratio2-step_ratio3*step_ratio3)/(1.0-step_ratio3*step_ratio3);
              for(i=0; i<nvarele; i++) {
                xc24[i]+=varchange[i]*extrap_w2;
                xcf[i]+=varchange[i]*extrap_w3;
              }
            }
          }

          if(sol==maxsol-1){
          if(subindx==0){
          for(i=0; i<nvar; i++) {
            for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
               elem_vals1[tindx1].substep_base=xcf[tindx1];
            }
          }
          }else{
          for(i=0; i<nvar; i++) {
            for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
               elem_vals1[tindx1].substep_base=(100+xcf[tindx1])/xc0[tindx1]-100;
            }
          }
          }
          for(i=0; i<ncofele; i++) elem_vals[i].value=elem_vals[i].initial;
          updates_apply_product(tabfile,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele);
          strcpy(commsyntax,"formula");
          IsIni=false;
          formulas_execute(tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,IsIni);
          for(i=0; i<nvar; i++) {
            for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
               elem_vals1[tindx1].substep_base=0;
            }
          }

          if(nohsl)MPI_Barrier(PETSC_COMM_WORLD);
          ierr = PetscGetCPUTime(&time1);
          CHKERRQ(ierr);
          if(verbosity>=1){ierr = PetscPrintf(PETSC_COMM_WORLD,"Last Update time %f\n",time1-time0);}
          CHKERRQ(ierr);
          }

          
          if(!(subindx==subints-1&&sol==maxsol-1)) {
            if(!inmemory){
            strcpy(tempfilenam,scratch_dir);
            strcat(tempfilenam,"_tempxcf");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
              printf("Error: cannot open %s for writing\n",tempfilenam);
            }
            fwrite(xcf, sizeof(solve_real),nvarele, tempvar);
            fclose(tempvar);
            free(*xcf2);
            *xcf2=NULL;
            xcf=*xcf2;
            }
            
            if(!inmemory){
            if(subindx>0&&sol==0){            strcpy(tempfilenam,scratch_dir);
            strcat(tempfilenam,"_tempxcO");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
              printf("Error: cannot open %s for writing\n",tempfilenam);
            }
            fwrite(xc0, sizeof(solve_real),nvarele, tempvar);
            fclose(tempvar);
            }
            /* xc0 must stay resident in inmemory mode: it carries the
               cumulative multiplier into the sol>0 extrapolation passes
               of each subinterval (the disk path reloads it from
               _tempxcO). Freeing it unconditionally silently dropped the
               sol>0 subinterval contributions under -inmemory. */
            free(xc0);
            xc0=NULL;
            }
          }

    if(rank==rank_hsl&&sol==maxsol-1) {
            if(subindx==0){
              xc124=realloc (xc124,nvarele*sizeof(int));
              for(i=0; i<nvarele; i++)xc124[i]=6;
            }else{
            if(!inmemory){
            strcpy(tempfilenam,scratch_dir);
            strcat(tempfilenam,"_tempxc124");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ((tempvar = fopen(tempfilenam, "rb")) == NULL) {
              printf("Error: cannot open %s for reading\n",tempfilenam);
            }
            xc124=realloc (xc124,nvarele*sizeof(int));
            freadresult=fread(xc124, sizeof(int),nvarele, tempvar);
            fclose(tempvar);
            remove(tempfilenam);
            }
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
            if(!inmemory){
            strcpy(tempfilenam,scratch_dir);
            strcat(tempfilenam,"_tempxc124");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
              printf("Error: cannot open %s for writing\n",tempfilenam);
            }
            fwrite(xc124, sizeof(int),nvarele, tempvar);
            fclose(tempvar);
            free(xc124);
            xc124=NULL;
            }
      }
    }
            
          if(!(subindx==subints-1&&sol==maxsol-1)) {
            if(!inmemory){
            strcpy(tempfilenam,scratch_dir);
            strcat(tempfilenam,"_tempxc12");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
              printf("Error: cannot open %s for writing\n",tempfilenam);
            }
            fwrite(xc12, sizeof(solve_real),nvarele, tempvar);
            fclose(tempvar);
            free(xc12);
            xc12=NULL;
            }
            if(!inmemory){
            strcpy(tempfilenam,scratch_dir);
            strcat(tempfilenam,"_tempxc24");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ( (tempvar = fopen(tempfilenam, "wb")) == NULL ) {
              printf("Error: cannot open %s for writing\n",tempfilenam);
            }
            fwrite(xc24, sizeof(solve_real),nvarele, tempvar);
            fclose(tempvar);
            free(xc24);
            xc24=NULL;
            }
          }
          
          
          if(!inmemory){
          /* per-subinterval solution snapshot; nothing reads it back */
          strcpy(solchar,scratch_dir);
          if(subindx<10)strcat(solchar,"_tempsol0");
          else strcat(solchar,"_tempsol");
          sprintf(tempchar, "%d", subindx);
          strcat(solchar,tempchar);
          sprintf(tempchar, "%d", sol);
          strcat(solchar,tempchar);
          strcat(solchar,".bin");
          logmsg(2,"solchar %s\n",solchar);
          if ( (solution = fopen(solchar, "wb")) == NULL ) {
            printf("Error: cannot open %s for writing\n",solchar);
            return 1;
          }
          fwrite(xcf, sizeof(solve_real),nvarele, solution);
          fclose(solution);
          }
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
    xc0=realloc (xc0,sizeof(solve_real));
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
    if(rank==0)logmsg(1,"Gragg solve time %.2f s\n",(endtime.tv_sec - begintime.tv_sec)+((double)(endtime.tv_usec - begintime.tv_usec))/ 1000000);
              free(counteqs);
              free(counteqnoadds);
              free(countvarintra1s);
              elem_vals1=NULL;
              if(x1!=NULL)free(x1);
  }

