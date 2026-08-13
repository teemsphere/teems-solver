#include <teems_solver.h>
#include <hsl_kernels.h>

/* Runge-Kutta drivers (GEMPACK manual 26.5; Schiffmann 2022): explicit
   RK2/RK4 and embedded BoSha32/DoPri54 with component-by-component
   global error estimates and (embedded only) adaptive step size.

   One RK stage = one Jacobian refill + linear solve at the stage state,
   so every stage refactorizes exactly as the Gragg substeps do and all
   matrix methods (LU/SBBD/DBBD/NDBBD) work unchanged.  The driver
   integrates the cumulative solution X per variable element: the
   initial-based cumulative percent change for percent-change variables
   (RK.pdf eq. 7 linearizes stage combination in this scale: a stage
   solution V relative to its stage state contributes
   k = (1+X_stage/100)*V), or the cumulative ordinary change for
   change variables.  Exogenous elements follow the same
   uniform-in-levels path as the Gragg driver, which is linear in X, so
   their stage gradient is the constant h*S.

   No subintervals and no Richardson extrapolation: RK produces a single
   solution and (embedded flavors) the .acc error metrics replace the
   extrapolation accuracy census.  Value arrays stay resident for the
   whole solve — the Gragg !inmemory scratch spills are not mirrored
   here (the per-step snapshot the stages restore from must stay
   addressable between solves). */

extern void spec48_ssol2la_(int *INSIZE,int *IRN, int *JCN, solve_real *VA, solve_real *B, solve_real *X);
extern void spec48_single_(fortran_int *indata,int *irn, int *jcn,solve_real *b1, solve_real *values,solve_real *x1, int *neleperrow,int *ai1, MPI_Fint *fcomm);
extern void spec48_nomc66_(fortran_int *indata, int *jcn,solve_real *b1, solve_real *values,solve_real *x1, int *neleperrow, MPI_Fint *fcomm,fortran_int *rowptrin, fortran_int *colptrin);

/* Butcher tableaus (standard published coefficients; the GEMPACK paper
   names the schemes without reprinting them). a is row-major [s][s],
   lower triangle used. bhat==NULL marks a non-embedded scheme. order
   is the higher order p of the pair — the controller exponent (eq. 11). */
typedef struct {
  const char *name;
  int s;
  int order;
  const double *c;
  const double *a;
  const double *b;
  const double *bhat;
} rk_scheme;

static const double rk2_c[2]={0,0.5};
static const double rk2_a[4]={0,0, 0.5,0};
static const double rk2_b[2]={0,1.0};

static const double rk4_c[4]={0,0.5,0.5,1.0};
static const double rk4_a[16]={0,0,0,0, 0.5,0,0,0, 0,0.5,0,0, 0,0,1.0,0};
static const double rk4_b[4]={1.0/6,1.0/3,1.0/3,1.0/6};

static const double bosha_c[4]={0,1.0/2,3.0/4,1.0};
static const double bosha_a[16]={0,0,0,0,
                                 1.0/2,0,0,0,
                                 0,3.0/4,0,0,
                                 2.0/9,1.0/3,4.0/9,0};
static const double bosha_b[4]={2.0/9,1.0/3,4.0/9,0};
static const double bosha_bh[4]={7.0/24,1.0/4,1.0/3,1.0/8};

static const double dopri_c[7]={0,1.0/5,3.0/10,4.0/5,8.0/9,1.0,1.0};
static const double dopri_a[49]={
  0,0,0,0,0,0,0,
  1.0/5,0,0,0,0,0,0,
  3.0/40,9.0/40,0,0,0,0,0,
  44.0/45,-56.0/15,32.0/9,0,0,0,0,
  19372.0/6561,-25360.0/2187,64448.0/6561,-212.0/729,0,0,0,
  9017.0/3168,-355.0/33,46732.0/5247,49.0/176,-5103.0/18656,0,0,
  35.0/384,0,500.0/1113,125.0/192,-2187.0/6784,11.0/84,0};
static const double dopri_b[7]={35.0/384,0,500.0/1113,125.0/192,-2187.0/6784,11.0/84,0};
static const double dopri_bh[7]={5179.0/57600,0,7571.0/16695,393.0/640,-92097.0/339200,187.0/2100,1.0/40};

/* One linearized solve at the current elem_vals state: build A and B,
   multiply the stage shock vector into the RHS, dispatch to the matrix
   method, return the solution in x1[VecSize].  Consumes (destroys)
   vece.  Mirrors the solve block of the Gragg step loop, including the
   per-branch A/vecb ownership: the (N)DBBD solvers release them
   internally, the SBBD/LU paths destroy them here. */
static void rk_stage_solve(PetscBool nohsl,PetscInt VecSize,PetscInt BSize,
  PetscInt dnz,PetscInt* dnnz,PetscInt onz,PetscInt* onnz,
  PetscInt dnzB,PetscInt* dnnzB,PetscInt onzB,PetscInt* onnzB,
  PetscInt rank,PetscInt rank_hsl,PetscInt mpisize,
  char *tabfile,char *commsyntax,
  set_def *sets,dim_t nset,set_element *set_elems,
  array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar,
  elem_value *elem_vals,offset_t ncofele,offset_t nvarele,
  closure_entry *closure_vals,
  offset_t alltimeset,offset_t allregset,offset_t nintraeq,
  dim_t matsol,PetscInt Istart,PetscInt Iend,
  offset_t nreg,offset_t ntime,offset_t *eq_addr,offset_t ndblock,
  offset_t *countvarintra1,offset_t *counteq,offset_t *counteqnoadd,
  offset_t *countvarintra1s,offset_t *counteqs,offset_t *counteqnoadds,
  dim_t laA,dim_t laDi,dim_t laD,PetscReal cntl3,PetscReal cntl6,
  dim_t nesteddbbd,int localsize,PetscInt *ndbbddrank1,
  fortran_int *indata,dim_t mc66,fortran_int *ptx,MPI_Fint fcomm,
  Vec *vece,solve_real *x1) {
  PetscErrorCode ierr;
  PetscScalar *vals=NULL;
  PetscInt count=0,nz01=0,*ai=NULL,*aj=NULL;
  fortran_int k=0,m=1;
  /* -fastrefac: sequential LU keeps the MA48 pivot sequence across
     stages and refactorizes with JOB=2; option read once per process */
  static dim_t fastrefac=-1;
  if(fastrefac<0) {
    fastrefac=0;
    PetscOptionsGetInt(NULL,NULL,"-fastrefac",&fastrefac,NULL);
  }
  offset_t i,j;
  solve_real *b1=NULL;
  PetscBool presol;
  Mat A,B;
  Vec vecb;

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
    MatMPIAIJSetPreallocation(A,dnz,dnnz,onz,onnz);
  }
  else {
    MatSetType(A,MATSEQAIJ);
    MatSeqAIJSetPreallocation(A,dnz,dnnz);
  }

  if(nohsl) {
    MatCreate(PETSC_COMM_WORLD,&B);
  }
  else {
    MatCreate(PETSC_COMM_SELF,&B);
  }
  shock_mat_set_sizes(B,nesteddbbd,localsize,VecSize,BSize);
  if(nohsl) {
    MatSetType(B,MATMPIAIJ);
    MatMPIAIJSetPreallocation(B,dnzB,dnnzB,onzB,onnzB);
  }
  else {
    MatSetType(B,MATSEQAIJ);
    MatSeqAIJSetPreallocation(B,dnzB,dnnzB);
  }

  if(rank==rank_hsl) {
    jacobian_fill(tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,closure_vals,ndblock,alltimeset,allregset,eq_addr,counteq,nintraeq,A,B);
  }
  MPI_Barrier(PETSC_COMM_WORLD);
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
  ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
  ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
  /* vecb spans the equation rows (VecSize); vece may be wider (BSize) */
  if(nohsl) {
    VecCreate(PETSC_COMM_WORLD,&vecb);
    VecSetType(vecb,VECMPI);
  }
  else {
    VecCreate(PETSC_COMM_SELF,&vecb);
    VecSetType(vecb,VECSEQ);
  }
  if(nesteddbbd==1)VecSetSizes(vecb,localsize,VecSize);
  else VecSetSizes(vecb,PETSC_DECIDE,VecSize);
  /* the (N)DBBD back-solves probe vecb with -1 sentinels for rows other
     ranks own */
  VecSetOption(vecb, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
  if(rank==rank_hsl) {
    ierr = MatMult(B,*vece,vecb);CHKERRV(ierr);
  }
  ierr = VecDestroy(vece);CHKERRV(ierr);
  ierr = VecAssemblyBegin(vecb);CHKERRV(ierr);
  ierr = VecAssemblyEnd(vecb);CHKERRV(ierr);
  ierr = MatDestroy(&B);CHKERRV(ierr);

  if(matsol>=MM_DBBD) {
    int *row_order= (int *) calloc (VecSize,sizeof(int));
    int *col_order= (int *) calloc (VecSize,sizeof(int));
    int *block_sizes= (int *) calloc (ndblock,sizeof(int));
    if(matsol==MM_DBBD) {
      dbbd_order(A,VecSize,mpisize,rank,Istart,Iend,nvarele,eq_addr,row_order,col_order,ndblock,block_sizes,countvarintra1,counteq,counteqnoadd,laA,cntl6);
      dbbd_solve(A,vecb,x1,VecSize,mpisize,rank,Istart,Iend,row_order,col_order,ndblock,block_sizes,countvarintra1,counteq,counteqnoadd,laA,laD,cntl3);
    }
    if(matsol==MM_NDBBD) {
      presol=1;
      memcpy(counteq,counteqs,(ndblock+1)*sizeof(offset_t));
      memcpy(counteqnoadd,counteqnoadds,(ndblock)*sizeof(offset_t));
      memcpy(countvarintra1,countvarintra1s,(ndblock+1)*sizeof(offset_t));
      ndbbd_order_presolve(A,VecSize,mpisize,rank,Istart,Iend,nreg,ntime,nvarele,eq_addr,row_order,col_order,ndblock,block_sizes,countvarintra1,counteq,counteqnoadd,laA,laDi,cntl6,ndbbddrank1,presol);
      ndbbd_presolve(A,vecb,x1,VecSize,mpisize,rank,Istart,Iend,row_order,col_order,ndblock,nreg,ntime,block_sizes,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,presol);
      presol=0;
      ndbbd_order(A,VecSize,mpisize,rank,Istart,Iend,nreg,ntime,nvarele,eq_addr,row_order,col_order,ndblock,block_sizes,countvarintra1,counteq,counteqnoadd,laA,laDi,cntl6,ndbbddrank1,presol);
      ndbbd_solve(A,vecb,x1,VecSize,mpisize,rank,Istart,Iend,row_order,col_order,ndblock,nreg,ntime,block_sizes,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,cntl3,cntl6,presol);
    }
    free(row_order);
    free(col_order);
    free(block_sizes);
    MPI_Barrier(PETSC_COMM_WORLD);
  }
  else {
    if(matsol==MM_SBBD&&fastrefac) {
      /* -fastrefac: persistent MP48 instance, FACT_JOB=2 per stage */
      sbbd_fastrefac_solve(&A,&vecb,VecSize,rank,rank_hsl,indata,fcomm,counteq,countvarintra1,x1);
    }
    else if(matsol==MM_SBBD) {
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
      ierr = MatDestroy(&A);CHKERRV(ierr);
      free(irn1);
      b1=(solve_real *) calloc (VecSize,sizeof(solve_real));
      if(rank==rank_hsl) {
        VecGetArray(vecb,&vals);
        for(i=0; i<VecSize; i++) {
          b1[i]=vals[i];
        }
      }
      ierr = VecDestroy(&vecb);CHKERRV(ierr);
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
      indata[1]=VecSize;
      indata[0]=count;
      ptx = indata;
      if(mc66!=0)spec48_single_(ptx,irn,jcn,b1,values,x1,neleperrow,ai1,&fcomm);
      free(irn);
      if(mc66==0)spec48_nomc66_(ptx,jcn,b1,values,x1,neleperrow,&fcomm,counteq,countvarintra1);
      free(jcn);
      free(values);
      free(neleperrow);
      free(ai1);
      free(b1);
      b1=NULL;
    }
    else if(fastrefac) {
      /* -fastrefac: persistent pivot sequence, MA48B/BD JOB=2 per stage */
      if(rank==rank_hsl) {
        VecGetArray(vecb,&vals);
        lu_fastrefac_solve(A,VecSize,laA,vals,x1);
      }
      ierr = MatDestroy(&A);CHKERRV(ierr);
      ierr = VecDestroy(&vecb);CHKERRV(ierr);
    }
    else {
      /* A stays live through the factorize so a workspace-growth
         retry can re-stage the COO (MA48 clobbers it in place) */
      if(rank==rank_hsl) {
        VecGetArray(vecb,&vals);
        lu_grow_solve(A,VecSize,laA,vals,x1);
      }
      ierr = MatDestroy(&A);CHKERRV(ierr);
      ierr = VecDestroy(&vecb);CHKERRV(ierr);
    }
  }
}

/* Restore the step-base state and advance it by the X-space increment
   dx[]: coefficient values return to the snapshot, each variable
   element gets value/substep_base consistent with X = base_X + dx (the
   substep_base a following updates_apply_product reads is the percent (or
   ordinary) change relative to the restored base state).  stagex
   receives the resulting X. */
static void rk_state_set(array_def *vars,offset_t nvar,
  elem_value *elem_vals,offset_t ncofele,offset_t nvarele,
  store_real *base_vals,solve_real *base_X,solve_real *dx,
  solve_real *stagex) {
  offset_t i;
  fortran_int tindx1;
  elem_value *elem_vals1=elem_vals+ncofele;
  for(i=0; i<ncofele; i++)elem_vals[i].value=base_vals[i];
  for(i=0; i<nvar; i++) {
    if(vars[i].change_real) {
      for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
        stagex[tindx1]=base_X[tindx1]+dx[tindx1];
        elem_vals1[tindx1].value=base_vals[ncofele+tindx1]+(store_real)dx[tindx1];
        elem_vals1[tindx1].substep_base=dx[tindx1];
      }
    }
    else {
      for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
        stagex[tindx1]=base_X[tindx1]+dx[tindx1];
        elem_vals1[tindx1].value=(100+stagex[tindx1])/100*elem_vals1[tindx1].initial;
        elem_vals1[tindx1].substep_base=dx[tindx1]/(1+base_X[tindx1]/100);
      }
    }
  }
}

bool solve_rk(PetscBool nohsl,PetscInt VecSize,PetscInt dnz,PetscInt* dnnz,PetscInt onz,PetscInt* onnz,PetscInt dnzB,PetscInt* dnnzB,PetscInt onzB,PetscInt* onnzB,Vec *vece1,PetscInt rank,PetscInt rank_hsl,PetscInt mpisize,char* tabfile, char *commsyntax,set_def *sets,dim_t nset, set_element *set_elems, array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar, elem_value **elem_vals2,offset_t ncofele,offset_t nvarele,closure_entry **closure_vals2,offset_t alltimeset,offset_t allregset,offset_t nintraeq,dim_t matsol,PetscInt Istart,PetscInt Iend,offset_t nreg, offset_t ntime, offset_t *eq_addr, offset_t ndblock, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,dim_t laDi,dim_t laD,PetscReal cntl3,PetscReal cntl6,dim_t nesteddbbd,int localsize,PetscInt *ndbbddrank1,fortran_int* indata,dim_t mc66,fortran_int *ptx,struct timeval begintime,MPI_Fint fcomm,int solmethod,int adaptive,double epstol,double retryadj,int maxretries,solve_real **xcf2,solve_real **accmetric2) {
  PetscErrorCode ierr;
  offset_t i,j;
  fortran_int tindx1;
  int m,st;
  bool IsIni;
  struct timeval endtime;
  elem_value *elem_vals;
  elem_value *elem_vals1;
  closure_entry *closure_vals;
  elem_vals=*elem_vals2;
  closure_vals=*closure_vals2;
  Vec vece;
  vece=*vece1;
  solve_real *xcf;

  rk_scheme scheme;
  switch(solmethod) {
  case SM_RK2:
    scheme=(rk_scheme) {"RK2",2,2,rk2_c,rk2_a,rk2_b,NULL};
    break;
  case SM_RK4:
    scheme=(rk_scheme) {"RK4",4,4,rk4_c,rk4_a,rk4_b,NULL};
    break;
  case SM_BOSHA32:
    scheme=(rk_scheme) {"BoSha32",4,3,bosha_c,bosha_a,bosha_b,bosha_bh};
    break;
  default:
    scheme=(rk_scheme) {"DoPri54",7,5,dopri_c,dopri_a,dopri_b,dopri_bh};
    break;
  }
  int s=scheme.s;
  bool embedded=(scheme.bhat!=NULL);

  /* NDBBD's ordering mutates the block-count arrays; keep pristine
     copies to restore before each ordering pass (as the Gragg loop does) */
  offset_t *counteqs= (offset_t *) calloc (ndblock+1,sizeof(offset_t));
  offset_t *counteqnoadds= (offset_t *) calloc (ndblock,sizeof(offset_t));
  offset_t *countvarintra1s= (offset_t *) calloc (ndblock+1,sizeof(offset_t));
  memcpy(counteqs,counteq,(ndblock+1)*sizeof(offset_t));
  memcpy(counteqnoadds,counteqnoadd,(ndblock)*sizeof(offset_t));
  memcpy(countvarintra1s,countvarintra1,(ndblock+1)*sizeof(offset_t));

  gettimeofday(&begintime, NULL);

  /* exogenous columns run 0..nexo-1; under heavy condensation nexo can
     exceed VecSize, so vece and B's columns span BSize (see
     shock_vec_set_sizes / shock_mat_set_sizes for the layout rule) */
  PetscInt BSize;
  BSize=(PetscInt)(nvarele-VecSize-nbselems);      /* nexo */
  BSize=(BSize>VecSize)?BSize:VecSize;

  solve_real *x1= (solve_real *) calloc (VecSize,sizeof(solve_real));
  /* per-step state: X = cumulative solution per element (initial-based
     percent, or ordinary change), snapshots of every value at the step
     base, the s stage gradients in X-space, and the stage scratch */
  solve_real *base_X= (solve_real *) calloc (nvarele,sizeof(solve_real));
  solve_real *newx= (solve_real *) calloc (nvarele,sizeof(solve_real));
  solve_real *stagex= (solve_real *) calloc (nvarele,sizeof(solve_real));
  solve_real *dx= (solve_real *) calloc (nvarele,sizeof(solve_real));
  store_real *base_vals= (store_real *) calloc (ncofele+nvarele,sizeof(store_real));
  solve_real *kk[7];
  for(m=0; m<s; m++)kk[m]= (solve_real *) calloc (nvarele,sizeof(solve_real));
  solve_real *dcum=NULL,*errx=NULL;
  if(embedded) {
    dcum= (solve_real *) calloc (nvarele,sizeof(solve_real));
    errx= (solve_real *) calloc (nvarele,sizeof(solve_real));
  }
  /* backsolve recovery workspace (exo changes as fed to vece + the
     recovered per-stage changes of the backsolved elements) */
  solve_real *exo_z=NULL,*bsvals=NULL;
  if(nbselems>0) {
    exo_z= (solve_real *) calloc (nvarele,sizeof(solve_real));
    bsvals= (solve_real *) calloc (nbselems,sizeof(solve_real));
  }

  elem_vals1=elem_vals+ncofele;
  /* solution start: pin the pre-simulation levels and snapshot the base */
  for(i=0; i<ncofele; i++) {
    elem_vals[i].initial=elem_vals[i].value;
    base_vals[i]=elem_vals[i].value;
  }
  for(i=0; i<nvar; i++) {
    for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
      elem_vals1[tindx1].initial=elem_vals1[tindx1].value;
      base_vals[ncofele+tindx1]=elem_vals1[tindx1].value;
      base_X[tindx1]=0;
      stagex[tindx1]=0;
    }
  }

  double t=0,h=1.0/steps1,q;
  int stepno=0,checkretries=0;
  bool firstsolve=true,crossed_warned=false;
  double metric,stepdata[2];

  while(t<1.0-1e-12) {
    if(t+h>1.0)h=1.0-t;   /* land exactly on t=1 */
    for(m=0; m<s; m++) {
      if(m>0) {
        for(tindx1=0; tindx1<nvarele; tindx1++) {
          dx[tindx1]=0;
          for(st=0; st<m; st++)dx[tindx1]+=scheme.a[m*s+st]*kk[st][tindx1];
        }
        rk_state_set(vars,nvar,elem_vals,ncofele,nvarele,base_vals,base_X,dx,stagex);
        if(rank==rank_hsl) {
          /* product form, not updates_apply(0): the flag-0 update
             linearizes multi-factor product updates, an O(h^2)
             coefficient error per stage that caps every flavor at
             global first order; the product form is exact per factor
             when compounding from the step base */
          updates_apply_product(tabfile,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele);
          strcpy(commsyntax,"formula");
          IsIni=false;
          formulas_execute(tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,IsIni);
assertions_execute(tabfile,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,IsIni,teems_assertions_mode,0);
        }
      }
      /* stage shock vector: the full-substep exogenous change evaluated
         at the stage state.  The exogenous path is uniform in levels,
         which is linear in X, so the X-space gradient is the constant
         h*S per step and only the percent fed to the solver rescales
         with the stage state. */
      if(!firstsolve) {
        if(nohsl) {
          VecCreate(PETSC_COMM_WORLD,&vece);
          VecSetType(vece,VECMPI);
        }
        else {
          VecCreate(PETSC_COMM_SELF,&vece);
          VecSetType(vece,VECSEQ);
        }
        shock_vec_set_sizes(vece,nesteddbbd,localsize,VecSize,BSize);
        VecSetOption(vece, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
      }
      firstsolve=false;
      for(i=0; i<nvar; i++) {
        if(vars[i].change_real) {
          for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
            if(closure_vals[tindx1].is_exogenous) {
              solve_real ez=h*closure_vals[tindx1].shock_value;
              VecSetValue(vece,closure_vals[tindx1].exo_index,ez,INSERT_VALUES);
              if(exo_z!=NULL)exo_z[tindx1]=ez;
            }
          }
        }
        else {
          for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
            if(closure_vals[tindx1].is_exogenous) {
              solve_real ez=h*closure_vals[tindx1].shock_value/(1+stagex[tindx1]/100);
              VecSetValue(vece,closure_vals[tindx1].exo_index,ez,INSERT_VALUES);
              if(exo_z!=NULL)exo_z[tindx1]=ez;
            }
          }
        }
      }
      MPI_Barrier(PETSC_COMM_WORLD);
      ierr = VecAssemblyBegin(vece);CHKERRQ(ierr);
      ierr = VecAssemblyEnd(vece);CHKERRQ(ierr);

      rk_stage_solve(nohsl,VecSize,BSize,dnz,dnnz,onz,onnz,dnzB,dnnzB,onzB,onnzB,rank,rank_hsl,mpisize,tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele,nvarele,closure_vals,alltimeset,allregset,nintraeq,matsol,Istart,Iend,nreg,ntime,eq_addr,ndblock,countvarintra1,counteq,counteqnoadd,countvarintra1s,counteqs,counteqnoadds,laA,laDi,laD,cntl3,cntl6,nesteddbbd,localsize,ndbbddrank1,indata,mc66,ptx,fcomm,&vece,x1);

      /* recover the backsolved elements from their defining equations
         with this stage's solution (GEMPACK 14.1.3) */
      if(rank==rank_hsl&&nbselems>0)backsolve_recover(tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele,closure_vals,x1,exo_z,bsvals);

      /* stage gradient in X-space: a percent solved relative to the
         stage state scales by (1+X_stage/100) (RK.pdf eq. 7); the
         exogenous gradient is the constant h*S */
      for(i=0; i<nvar; i++) {
        if(vars[i].change_real) {
          for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
            if(closure_vals[tindx1].is_exogenous)kk[m][tindx1]=h*closure_vals[tindx1].shock_value;
            else if(closure_vals[tindx1].is_backsolved)kk[m][tindx1]=bsvals[closure_vals[tindx1].exo_index];
            else kk[m][tindx1]=x1[closure_vals[tindx1].exo_index];
          }
        }
        else {
          for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
            if(closure_vals[tindx1].is_exogenous)kk[m][tindx1]=h*closure_vals[tindx1].shock_value;
            else if(closure_vals[tindx1].is_backsolved)kk[m][tindx1]=bsvals[closure_vals[tindx1].exo_index]*(1+stagex[tindx1]/100);
            else kk[m][tindx1]=x1[closure_vals[tindx1].exo_index]*(1+stagex[tindx1]/100);
          }
        }
      }
    }

    /* combine: X_new = X_base + sum b_j k_j; embedded error
       eps = sum (b_j - bhat_j) k_j, already in X-space */
    bool crossed=false;
    for(i=0; i<nvar; i++) {
      for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
        solve_real acc=0,eacc=0;
        for(st=0; st<s; st++) {
          acc+=scheme.b[st]*kk[st][tindx1];
          if(embedded)eacc+=(scheme.b[st]-scheme.bhat[st])*kk[st][tindx1];
        }
        newx[tindx1]=base_X[tindx1]+acc;
        if(embedded)errx[tindx1]=eacc;
        if(!vars[i].change_real&&newx[tindx1]<=-100)crossed=true;
      }
    }
    metric=0;
    if(embedded) {
      for(tindx1=0; tindx1<nvarele; tindx1++) {
        double v=fabs(newx[tindx1]);
        double e=fabs(errx[tindx1])/((v>1)?v:1);
        if(e>metric)metric=e;
      }
    }
    /* the accept/reject branch must be collective: every rank follows
       rank 0's metric (Allreduce-equivalent; all ranks hold the full
       solution, HSL worker ranks hold stale copies) */
    stepdata[0]=metric;
    stepdata[1]=crossed?1:0;
    MPI_Bcast(stepdata,2,MPI_DOUBLE,0,PETSC_COMM_WORLD);
    metric=stepdata[0];
    crossed=(stepdata[1]!=0);

    if(adaptive&&embedded) {
      if(crossed&&adaptive==1) { /* adaptive=yes retries check failures;
                                    accuracy-only ignores them */
        checkretries++;
        if(checkretries>maxretries) {
          if(rank==0)printf("Error: a percent-change variable crossed -100%% and %d retries at reduced step size did not clear it; the simulation cannot continue (try more steps or a smaller shock)\n",maxretries);
          MPI_Abort(PETSC_COMM_WORLD,1);
        }
        h*=retryadj;
        if(rank==0)logmsg(1,"Step %d: percent variable crossed -100%%, retrying with step size %.4g\n",stepno+1,h);
        /* redo the step from its base: stage 0 reads only the
           coefficient values and stagex, so those are all that need
           restoring */
        for(i=0; i<ncofele; i++)elem_vals[i].value=base_vals[i];
        for(tindx1=0; tindx1<nvarele; tindx1++)stagex[tindx1]=base_X[tindx1];
        continue;
      }
      if(metric>epstol) {
        q=0.85*pow(epstol/metric,1.0/scheme.order);
        if(q<0.5)q=0.5;
        h*=q;
        if(h<1e-8) {
          if(rank==0)printf("Error: adaptive step size collapsed below 1e-8 without meeting epsTolerance %.3g (worst error metric %.3g); the problem may be too stiff for %s\n",epstol,metric,scheme.name);
          MPI_Abort(PETSC_COMM_WORLD,1);
        }
        if(rank==0)logmsg(1,"Step %d rejected (accuracy metric %.3g > %.3g), retrying with step size %.4g\n",stepno+1,metric,epstol,h);
        for(i=0; i<ncofele; i++)elem_vals[i].value=base_vals[i];
        for(tindx1=0; tindx1<nvarele; tindx1++)stagex[tindx1]=base_X[tindx1];
        continue;
      }
      checkretries=0;
    }
    else if(crossed&&!crossed_warned) {
      if(rank==0)printf("Warning: a percent-change variable crossed -100%% during the %s integration; results past this point are unreliable (use more steps or an adaptive method)\n",scheme.name);
      crossed_warned=true;
    }

    /* accept: propagate the cumulative error estimate (RK.pdf eq. 6/8;
       eps is already X-space so the (1+X/100) factor is folded in and
       the prior error scales with the step's multiplier) */
    if(embedded) {
      for(i=0; i<nvar; i++) {
        if(vars[i].change_real) {
          for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++)
            dcum[tindx1]+=fabs(errx[tindx1]);
        }
        else {
          for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
            solve_real vstep=(newx[tindx1]-base_X[tindx1])/(1+base_X[tindx1]/100);
            dcum[tindx1]=(1+vstep/100)*dcum[tindx1]+fabs(errx[tindx1]);
          }
        }
      }
    }
    /* advance the resident state to X_new and make it the next base */
    for(tindx1=0; tindx1<nvarele; tindx1++)dx[tindx1]=newx[tindx1]-base_X[tindx1];
    rk_state_set(vars,nvar,elem_vals,ncofele,nvarele,base_vals,base_X,dx,stagex);
    if(rank==rank_hsl) {
      updates_apply_product(tabfile,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele);
      strcpy(commsyntax,"formula");
      IsIni=false;
      formulas_execute(tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,IsIni);
assertions_execute(tabfile,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,IsIni,teems_assertions_mode,0);
    }
    for(i=0; i<ncofele; i++)base_vals[i]=elem_vals[i].value;
    for(i=0; i<nvar; i++) {
      for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
        base_vals[ncofele+tindx1]=elem_vals1[tindx1].value;
        base_X[tindx1]=newx[tindx1];
        stagex[tindx1]=newx[tindx1];
      }
    }
    t+=h;
    stepno++;
    if(rank==0)logmsg(2,"step %d done, t %.6f h %.4g metric %.3g\n",stepno,t,h,metric);
    if(adaptive&&embedded&&t<1.0-1e-12) {
      q=0.85*pow(epstol/((metric>1e-300)?metric:1e-300),1.0/scheme.order);
      if(q<0.5)q=0.5;
      if(q>2.0)q=2.0;
      h*=q;
    }
  }

  /* the cumulative solution is the solution: no extrapolation pass */
  *xcf2=(solve_real*)realloc (*xcf2,nvarele*sizeof(solve_real));
  xcf=*xcf2;
  for(tindx1=0; tindx1<nvarele; tindx1++)xcf[tindx1]=base_X[tindx1];

  if(embedded) {
    /* error metric per component, E = |Dcum/max(1,|X|)| (RK.pdf eq. 5),
       written to the .acc ride-along by main; face value in brackets
       of 0.02 on the worst metric */
    *accmetric2=(solve_real*)realloc (*accmetric2,nvarele*sizeof(solve_real));
    solve_real *accm=*accmetric2;
    double worst=0;
    for(tindx1=0; tindx1<nvarele; tindx1++) {
      double v=fabs(xcf[tindx1]);
      accm[tindx1]=fabs(dcum[tindx1])/((v>1)?v:1);
      if(accm[tindx1]>worst)worst=accm[tindx1];
    }
    int face=10-(int)floor(worst/0.02);
    if(face<1)face=1;
    if(face>10)face=10;
    if(rank==0)printf("Accuracy face value %d of 10 (worst accuracy metric %.3g)\n",face,worst);
  }

  /* final data pass: rebase the coefficients on the pre-simulation
     levels and apply the whole solution in one product update, exactly
     as the extrapolating drivers do */
  if(rank==rank_hsl) {
    for(i=0; i<nvar; i++) {
      for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
        elem_vals1[tindx1].substep_base=xcf[tindx1];
      }
    }
    for(i=0; i<ncofele; i++) elem_vals[i].value=elem_vals[i].initial;
    updates_apply_product(tabfile,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele);
    strcpy(commsyntax,"formula");
    IsIni=false;
    formulas_execute(tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,IsIni);
assertions_execute(tabfile,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,IsIni,teems_assertions_mode,0);
    for(i=0; i<nvar; i++) {
      for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
        elem_vals1[tindx1].substep_base=0;
      }
    }
  }

  free(x1);
  free(base_X);
  free(newx);
  free(stagex);
  free(dx);
  free(base_vals);
  for(m=0; m<s; m++)free(kk[m]);
  free(dcum);
  free(errx);
  free(exo_z);
  free(bsvals);
  free(counteqs);
  free(counteqnoadds);
  free(countvarintra1s);
  gettimeofday(&endtime, NULL);
  if(rank==0)logmsg(1,"%s solve time %.2f s (%d steps)\n",scheme.name,(endtime.tv_sec - begintime.tv_sec)+((double)(endtime.tv_usec - begintime.tv_usec))/ 1000000,stepno);
  return 1;
}

/* ------------------------------------------------------------------ */
/* C2 complementarity approximate run (design doc section 8; manual
   51.1.2/51.6/51.7.3): single-solution forward Euler over t in [0,1]
   with the E_$comp state machinery live. Per step: the per-component
   states and weight coefficients are evaluated from the step-base
   levels values (comp_states_set), the Newton-correction variable
   del_comp@ is shocked 1 IN FULL (NO_SPLIT semantics, 51.7.2 (d)),
   one Euler step is solved and applied, and the states are
   re-evaluated (comp_states_check); if any component changed state
   the step is redone from its base with a shorter length so the
   change lands just before the redone step's end (51.7.3), never
   shorter than redo_min_frac of the original step, and a redone step
   is accepted regardless. The accepted-step count can therefore
   exceed the request (51.6). The accurate run (51.7.1) is C3; this
   run's solution is the simulation result. */
bool solve_comp_approx(PetscBool nohsl,PetscInt VecSize,PetscInt dnz,PetscInt* dnnz,PetscInt onz,PetscInt* onnz,PetscInt dnzB,PetscInt* dnnzB,PetscInt onzB,PetscInt* onnzB,Vec *vece1,PetscInt rank,PetscInt rank_hsl,PetscInt mpisize,char* tabfile, char *commsyntax,set_def *sets,dim_t nset, set_element *set_elems, array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar, elem_value **elem_vals2,offset_t ncofele,offset_t nvarele,closure_entry **closure_vals2,offset_t alltimeset,offset_t allregset,offset_t nintraeq,dim_t matsol,PetscInt Istart,PetscInt Iend,offset_t nreg, offset_t ntime, offset_t *eq_addr, offset_t ndblock, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,dim_t laDi,dim_t laD,PetscReal cntl3,PetscReal cntl6,dim_t nesteddbbd,int localsize,PetscInt *ndbbddrank1,fortran_int* indata,dim_t mc66,fortran_int *ptx,struct timeval begintime,MPI_Fint fcomm,int napprox,int redo_steps,double redo_min_frac,solve_real **xcf2) {
  PetscErrorCode ierr;
  offset_t i;
  fortran_int tindx1;
  bool IsIni;
  struct timeval endtime;
  elem_value *elem_vals;
  elem_value *elem_vals1;
  closure_entry *closure_vals;
  elem_vals=*elem_vals2;
  closure_vals=*closure_vals2;
  Vec vece;
  vece=*vece1;
  solve_real *xcf;

  offset_t *counteqs= (offset_t *) calloc (ndblock+1,sizeof(offset_t));
  offset_t *counteqnoadds= (offset_t *) calloc (ndblock,sizeof(offset_t));
  offset_t *countvarintra1s= (offset_t *) calloc (ndblock+1,sizeof(offset_t));
  memcpy(counteqs,counteq,(ndblock+1)*sizeof(offset_t));
  memcpy(counteqnoadds,counteqnoadd,(ndblock)*sizeof(offset_t));
  memcpy(countvarintra1s,countvarintra1,(ndblock+1)*sizeof(offset_t));

  gettimeofday(&begintime, NULL);

  PetscInt BSize;
  BSize=(PetscInt)(nvarele-VecSize-nbselems);      /* nexo */
  BSize=(BSize>VecSize)?BSize:VecSize;

  solve_real *x1= (solve_real *) calloc (VecSize,sizeof(solve_real));
  solve_real *base_X= (solve_real *) calloc (nvarele,sizeof(solve_real));
  solve_real *newx= (solve_real *) calloc (nvarele,sizeof(solve_real));
  solve_real *stagex= (solve_real *) calloc (nvarele,sizeof(solve_real));
  solve_real *dx= (solve_real *) calloc (nvarele,sizeof(solve_real));
  store_real *base_vals= (store_real *) calloc (ncofele+nvarele,sizeof(store_real));
  solve_real *exo_z=NULL,*bsvals=NULL;
  if(nbselems>0) {
    exo_z= (solve_real *) calloc (nvarele,sizeof(solve_real));
    bsvals= (solve_real *) calloc (nbselems,sizeof(solve_real));
  }

  /* the Newton-correction variable's single element (51.7.2 (iii)).
     Only the model-owning rank validates: under HSL the worker ranks
     carry an empty model (nvar 0) and every deloff consumer below
     sits in an nvar-bounded loop, so deloff=-1 is inert there --
     the unconditional abort here deadlocked/killed SBBD comp runs */
  offset_t deloff=-1;
  for(i=0; i<nvar; i++)if(strcmp(vars[i].cofname,"del_comp@")==0) {
      deloff=vars[i].offset;
      break;
    }
  if(deloff<0&&rank==rank_hsl) {
    printf("Error: complementarity approximate run without a del_comp@ variable (internal)\n");
    MPI_Abort(PETSC_COMM_WORLD,1);
  }

  elem_vals1=elem_vals+ncofele;
  for(i=0; i<ncofele; i++) {
    elem_vals[i].initial=elem_vals[i].value;
    base_vals[i]=elem_vals[i].value;
  }
  for(i=0; i<nvar; i++) {
    for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
      elem_vals1[tindx1].initial=elem_vals1[tindx1].value;
      base_vals[ncofele+tindx1]=elem_vals1[tindx1].value;
      base_X[tindx1]=0;
      stagex[tindx1]=0;
    }
  }

  double t=0,hdef=1.0/napprox,h=hdef;
  int stepno=0,redoing=0;
  bool firstsolve=true;
  double stepdata[2];

  while(t<1.0-1e-12) {
    if(t+h>1.0)h=1.0-t;
    /* per-step states + E_$comp weights from the step-base values
       (a redo re-evaluates the same base: same states, shorter h) */
    if(rank==rank_hsl) {
      if(comp_states_set(sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals)<0)MPI_Abort(PETSC_COMM_WORLD,1);
    }
    if(!firstsolve) {
      if(nohsl) {
        VecCreate(PETSC_COMM_WORLD,&vece);
        VecSetType(vece,VECMPI);
      }
      else {
        VecCreate(PETSC_COMM_SELF,&vece);
        VecSetType(vece,VECSEQ);
      }
      shock_vec_set_sizes(vece,nesteddbbd,localsize,VecSize,BSize);
      VecSetOption(vece, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
    }
    firstsolve=false;
    for(i=0; i<nvar; i++) {
      if(vars[i].change_real) {
        for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
          if(closure_vals[tindx1].is_exogenous) {
            solve_real ez=(tindx1==deloff)?1.0:h*closure_vals[tindx1].shock_value;
            VecSetValue(vece,closure_vals[tindx1].exo_index,ez,INSERT_VALUES);
            if(exo_z!=NULL)exo_z[tindx1]=ez;
          }
        }
      }
      else {
        for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
          if(closure_vals[tindx1].is_exogenous) {
            solve_real ez=h*closure_vals[tindx1].shock_value/(1+stagex[tindx1]/100);
            VecSetValue(vece,closure_vals[tindx1].exo_index,ez,INSERT_VALUES);
            if(exo_z!=NULL)exo_z[tindx1]=ez;
          }
        }
      }
    }
    MPI_Barrier(PETSC_COMM_WORLD);
    ierr = VecAssemblyBegin(vece);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(vece);CHKERRQ(ierr);

    rk_stage_solve(nohsl,VecSize,BSize,dnz,dnnz,onz,onnz,dnzB,dnnzB,onzB,onnzB,rank,rank_hsl,mpisize,tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele,nvarele,closure_vals,alltimeset,allregset,nintraeq,matsol,Istart,Iend,nreg,ntime,eq_addr,ndblock,countvarintra1,counteq,counteqnoadd,countvarintra1s,counteqs,counteqnoadds,laA,laDi,laD,cntl3,cntl6,nesteddbbd,localsize,ndbbddrank1,indata,mc66,ptx,fcomm,&vece,x1);

    if(rank==rank_hsl&&nbselems>0)backsolve_recover(tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele,closure_vals,x1,exo_z,bsvals);

    /* the Euler increment in X-space */
    for(i=0; i<nvar; i++) {
      if(vars[i].change_real) {
        for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
          if(closure_vals[tindx1].is_exogenous)dx[tindx1]=(tindx1==deloff)?1.0:h*closure_vals[tindx1].shock_value;
          else if(closure_vals[tindx1].is_backsolved)dx[tindx1]=bsvals[closure_vals[tindx1].exo_index];
          else dx[tindx1]=x1[closure_vals[tindx1].exo_index];
        }
      }
      else {
        for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
          if(closure_vals[tindx1].is_exogenous)dx[tindx1]=h*closure_vals[tindx1].shock_value;
          else if(closure_vals[tindx1].is_backsolved)dx[tindx1]=bsvals[closure_vals[tindx1].exo_index]*(1+stagex[tindx1]/100);
          else dx[tindx1]=x1[closure_vals[tindx1].exo_index]*(1+stagex[tindx1]/100);
        }
      }
    }
    for(tindx1=0; tindx1<nvarele; tindx1++)newx[tindx1]=base_X[tindx1]+dx[tindx1];
    /* trial advance: values, updates, formulas -- then the state check */
    rk_state_set(vars,nvar,elem_vals,ncofele,nvarele,base_vals,base_X,dx,stagex);
    if(rank==rank_hsl) {
      updates_apply_product(tabfile,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele);
      strcpy(commsyntax,"formula");
      IsIni=false;
      formulas_execute(tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,IsIni);
      assertions_execute(tabfile,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,IsIni,teems_assertions_mode,0);
    }
    {
      offset_t nflip=0;
      double frac=1.0;
      if(rank==rank_hsl)comp_states_check(sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,&nflip,&frac);
      stepdata[0]=(double)nflip;
      stepdata[1]=frac;
      /* root 0, NOT rank_hsl: under nohsl rank_hsl==rank on every
         rank, so a rank_hsl root made every rank the broadcaster --
         the orphaned eager sends then mis-matched the raw
         PETSC_COMM_WORLD collectives inside the next step's
         (N)DBBD ordering (deadlock at approx step 2).  Rank 0 holds
         valid data in both modes: it is rank_hsl under HSL, and under
         nohsl every rank computes the state check itself. */
      MPI_Bcast(stepdata,2,MPI_DOUBLE,0,PETSC_COMM_WORLD);
      nflip=(offset_t)stepdata[0];
      frac=stepdata[1];
      if(nflip>0&&redo_steps&&!redoing) {
        double hnew=h*frac/0.995;
        if(hnew<redo_min_frac*h)hnew=redo_min_frac*h;
        if(hnew>h)hnew=h;
        /* restore the step base: rk_state_set with a zero increment
           returns every coefficient and variable slot to the base */
        for(tindx1=0; tindx1<nvarele; tindx1++)dx[tindx1]=0;
        rk_state_set(vars,nvar,elem_vals,ncofele,nvarele,base_vals,base_X,dx,stagex);
        if(rank==0)logmsg(1,"Step %d: %ld complementarity state change(s); redoing the step at %.4g of its length (51.7.3)\n",stepno+1,(long)nflip,hnew/h);
        h=hnew;
        redoing=1;
        continue;
      }
      redoing=0;
    }
    /* accept */
    for(i=0; i<ncofele; i++)base_vals[i]=elem_vals[i].value;
    for(i=0; i<nvar; i++) {
      for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
        base_vals[ncofele+tindx1]=elem_vals1[tindx1].value;
        base_X[tindx1]=newx[tindx1];
        stagex[tindx1]=newx[tindx1];
      }
    }
    t+=h;
    stepno++;
    if(rank==0)logmsg(2,"approx step %d done, t %.6f h %.4g\n",stepno,t,h);
    h=hdef;
  }

  *xcf2=(solve_real*)realloc (*xcf2,nvarele*sizeof(solve_real));
  xcf=*xcf2;
  for(tindx1=0; tindx1<nvarele; tindx1++)xcf[tindx1]=base_X[tindx1];

  /* final data pass: rebase the coefficients on the pre-simulation
     levels and apply the whole solution in one product update */
  if(rank==rank_hsl) {
    for(i=0; i<nvar; i++) {
      for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
        elem_vals1[tindx1].substep_base=xcf[tindx1];
      }
    }
    for(i=0; i<ncofele; i++) elem_vals[i].value=elem_vals[i].initial;
    updates_apply_product(tabfile,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele);
    strcpy(commsyntax,"formula");
    IsIni=false;
    formulas_execute(tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,IsIni);
    assertions_execute(tabfile,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,IsIni,teems_assertions_mode,0);
    for(i=0; i<nvar; i++) {
      for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
        elem_vals1[tindx1].substep_base=0;
      }
    }
    /* 51.7.5 post-simulation state check + 51.5.3-style change lines
       (the state runtime stays live: the C3 accurate run reads the
       final states from it; main frees it after the last pass) */
    comp_states_report(sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals);
  }

  free(x1);
  free(base_X);
  free(newx);
  free(stagex);
  free(dx);
  free(base_vals);
  free(exo_z);
  free(bsvals);
  free(counteqs);
  free(counteqnoadds);
  free(countvarintra1s);
  gettimeofday(&endtime, NULL);
  if(rank==0)logmsg(1,"Complementarity approximate run solve time %.2f s (%d Euler steps, %d requested)\n",(endtime.tv_sec - begintime.tv_sec)+((double)(endtime.tv_usec - begintime.tv_usec))/ 1000000,stepno,napprox);
  return 1;
}
