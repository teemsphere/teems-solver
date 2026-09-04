#include <teems_solver.h>
#include <hsl_kernels.h>

/* Runge-Kutta drivers (GEMPACK manual 26.5; Schiffmann 2022; Munthe-Kaas
   1999 for the chart): explicit Heun/RK2/RK4 and embedded BoSha32/DoPri54
   with component-by-component global error estimates and (embedded only)
   adaptive step size.

   One RK stage = one Jacobian refill + linear solve at the stage state,
   so every stage refactorizes exactly as the Gragg substeps do and all
   matrix methods (LU/SBBD/DBBD/NDBBD) work unchanged.

   CHART.  The driver integrates a per-element state Z:
     - log chart (default, -rkchart log): a percent-change element with a
       nonzero base level is carried as Z = log(level/level0), the
       exponential coordinate of the multiplicative group of positive
       levels.  Munthe-Kaas (1999, Thm 6/7) on an abelian group: the
       dexpinv correction is the identity, so every classical tableau
       keeps its order and level0*exp(Z) is positive for any tableau and
       any step -- a percent variable cannot reach -100 at any stage.  A
       stage solution V (percent, relative to the stage state) is the
       gradient V/100 in this chart; GEMPACK's accumulation
       X_i = X_{i-1} + (1+X_{i-1}/100) V_i is its first-order truncation.
     - percent chart (-rkchart percent): the original arithmetic -- the
       initial-based cumulative percent change X with stage gradients
       (1+X_stage/100)*V (RK.pdf eq. 7), kept for GEMPACK-orientation
       regression.
     Ordinary-change elements and percent elements with a zero base
     level use the additive (percent-chart) arithmetic in both charts.
   Exogenous elements follow the level-linear path in both charts:
   dZ/dt = S/(100*exp(Z)) in the log chart, the constant h*S in X.

   ECONOMY.  First-same-as-last (BoSha32, DoPri54): the last stage is
   evaluated at the accepted end state, which is the next step's base, so
   its gradient is reused as the next stage 0 (rescaled by h_new/h_old).
   A rejected step keeps its stage-0 gradient the same way.  Stage states
   are checked as they are realized (percent crossing -100 in the percent
   chart, coefficient range tests, assertions, the log-chart magnitude
   guard) and a failing step is abandoned before the remaining stages are
   solved; under -adaptive yes these are the manual 26.5.1 retry triggers
   1-4, with the step retried at h*retryadj.

   No subintervals and no Richardson extrapolation: RK produces a single
   solution and (embedded flavors) the per-element estimated error metric
   (.est ride-along) replaces the extrapolation accuracy census.  Value
   arrays stay resident for the whole solve. */

extern void spec48_ssol2la_(int *INSIZE,int *IRN, int *JCN, solve_real *VA, solve_real *B, solve_real *X);
extern void spec48_single_(fortran_int *indata,int *irn, int *jcn,solve_real *b1, solve_real *values,solve_real *x1, int *neleperrow,int *ai1, MPI_Fint *fcomm);
extern void spec48_nomc66_(fortran_int *indata, int *jcn,solve_real *b1, solve_real *values,solve_real *x1, int *neleperrow, MPI_Fint *fcomm,fortran_int *rowptrin, fortran_int *colptrin);

/* Butcher tableaus (standard published coefficients; the GEMPACK paper
   names the schemes without reprinting them). a is row-major [s][s],
   lower triangle used. bhat==NULL marks a non-embedded scheme. order
   is the higher order p of the pair -- the controller exponent (eq. 11).
   fsal marks a_{s,j} == b_j with c_s == 1 (last stage state == end
   state). */
typedef struct {
  const char *name;
  int s;
  int order;
  int fsal;
  const double *c;
  const double *a;
  const double *b;
  const double *bhat;
} rk_scheme;

static const double rk2_c[2]={0,0.5};
static const double rk2_a[4]={0,0, 0.5,0};
static const double rk2_b[2]={0,1.0};

/* Heun (explicit trapezoid): strong-stability-preserving with C = 1,
   which the midpoint tableau is not (Butcher 2014 sect. 8) */
static const double heun_c[2]={0,1.0};
static const double heun_a[4]={0,0, 1.0,0};
static const double heun_b[2]={0.5,0.5};

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
  offset_t nreg,offset_t ntime,PetscInt *eq_addr,offset_t ndblock,
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
  if(nohsl)MatSetType(A,MATMPIAIJ);
  else MatSetType(A,MATSEQAIJ);
  jac_mat_prealloc(A,"Jacobian",nohsl,rank==rank_hsl,Iend-Istart,dnz,dnnz,onz,onnz);

  if(nohsl) {
    MatCreate(PETSC_COMM_WORLD,&B);
  }
  else {
    MatCreate(PETSC_COMM_SELF,&B);
  }
  shock_mat_set_sizes(B,nesteddbbd,localsize,VecSize,BSize);
  if(nohsl)MatSetType(B,MATMPIAIJ);
  else MatSetType(B,MATSEQAIJ);
  jac_mat_prealloc(B,"exogenous block",nohsl,rank==rank_hsl,Iend-Istart,dnzB,dnnzB,onzB,onnzB);

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
      if(mc66==0) {
        sbbd_csr_solve(&A,&vecb,VecSize,rank,rank_hsl,indata,fcomm,counteq,countvarintra1,x1);
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

/* Restore the step-base state and advance it by the chart increment
   delta[]: coefficient values return to the snapshot, each variable
   element gets value/substep_base consistent with Z = base_Z + delta
   (the substep_base a following updates_apply_product reads is the
   percent (or ordinary) change relative to the restored base state).
   stagez receives the resulting Z.  mult[t] != 0 marks a log-chart
   element (Z = log(level/level0)); mult == NULL is the percent chart
   throughout (the complementarity driver's use).  delta may alias
   stagez: each element's increment is read before its state is
   written. */
static void rk_state_set(array_def *vars,offset_t nvar,
  elem_value *elem_vals,offset_t ncofele,offset_t nvarele,
  store_real *base_vals,solve_real *base_Z,solve_real *delta,
  solve_real *stagez,const char *mult) {
  offset_t i;
  fortran_int tindx1;
  elem_value *elem_vals1=elem_vals+ncofele;
  for(i=0; i<ncofele; i++)elem_vals[i].value=base_vals[i];
  for(i=0; i<nvar; i++) {
    if(vars[i].change_real) {
      for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
        solve_real d=delta[tindx1];
        stagez[tindx1]=base_Z[tindx1]+d;
        elem_vals1[tindx1].value=base_vals[ncofele+tindx1]+(store_real)d;
        elem_vals1[tindx1].substep_base=d;
      }
    }
    else {
      for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
        solve_real d=delta[tindx1];
        if(mult!=NULL&&mult[tindx1]) {
          stagez[tindx1]=base_Z[tindx1]+d;
          elem_vals1[tindx1].value=(store_real)(elem_vals1[tindx1].initial*exp(stagez[tindx1]));
          elem_vals1[tindx1].substep_base=100*(exp(d)-1);
        }
        else {
          stagez[tindx1]=base_Z[tindx1]+d;
          elem_vals1[tindx1].value=(100+stagez[tindx1])/100*elem_vals1[tindx1].initial;
          elem_vals1[tindx1].substep_base=d/(1+base_Z[tindx1]/100);
        }
      }
    }
  }
}

/* Stage-state checks after a realization (rank_hsl evaluates the
   updates/formulas/assertions; the chart checks run on the state
   arrays).  Returns the manual 26.5.1 trigger class: 0 clean,
   RK_FAIL_CROSSED (percent element <= -100 in the additive
   arithmetic), RK_FAIL_RANGE / RK_FAIL_ASSERT (counted by formula.c
   while teems_rk_stage_checks is set), RK_FAIL_GUARD (log-chart
   magnitude guard). */
#define RK_FAIL_CROSSED 1
#define RK_FAIL_RANGE   2
#define RK_FAIL_ASSERT  3
#define RK_FAIL_GUARD   4
#define RK_FAIL_SINGULAR 5

static int rk_state_check(array_def *vars,offset_t nvar,solve_real *stagez,const char *mult,double guard) {
  offset_t i;
  fortran_int tindx1;
  int code=0;
  for(i=0; i<nvar&&code==0; i++) {
    if(vars[i].change_real)continue;
    for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
      if(mult!=NULL&&mult[tindx1]) {
        if(fabs(stagez[tindx1])>guard) {
          code=RK_FAIL_GUARD;
          break;
        }
      }
      else if(stagez[tindx1]<=-100) {
        code=RK_FAIL_CROSSED;
        break;
      }
    }
  }
  if(code==0&&teems_check_viol_assert>0)code=RK_FAIL_ASSERT;
  if(code==0&&teems_check_viol_range>0)code=RK_FAIL_RANGE;
  return code;
}

/* name of the variable owning element index t, for log lines */
static const char *rk_elem_var(array_def *vars,offset_t nvar,fortran_int t) {
  offset_t i;
  for(i=0; i<nvar; i++)if(t>=vars[i].offset&&t<vars[i].offset+vars[i].nelem)return vars[i].cofname;
  return "?";
}

static const char *rk_fail_name(int code) {
  switch(code) {
  case RK_FAIL_CROSSED:
    return "a percent-change variable crossed -100%";
  case RK_FAIL_RANGE:
    return "a coefficient range test failed";
  case RK_FAIL_ASSERT:
    return "an assertion did not hold";
  case RK_FAIL_GUARD:
    return "a level moved beyond the -rkguard ratio";
  case RK_FAIL_SINGULAR:
    return "the Jacobian was singular at a stage state";
  }
  return "a check failed";
}

bool solve_rk(PetscBool nohsl,PetscInt VecSize,PetscInt dnz,PetscInt* dnnz,PetscInt onz,PetscInt* onnz,PetscInt dnzB,PetscInt* dnnzB,PetscInt onzB,PetscInt* onnzB,Vec *vece1,PetscInt rank,PetscInt rank_hsl,PetscInt mpisize,char* tabfile, char *commsyntax,set_def *sets,dim_t nset, set_element *set_elems, array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar, elem_value **elem_vals2,offset_t ncofele,offset_t nvarele,closure_entry **closure_vals2,offset_t alltimeset,offset_t allregset,offset_t nintraeq,dim_t matsol,PetscInt Istart,PetscInt Iend,offset_t nreg, offset_t ntime, PetscInt *eq_addr, offset_t ndblock, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,dim_t laDi,dim_t laD,PetscReal cntl3,PetscReal cntl6,dim_t nesteddbbd,int localsize,PetscInt *ndbbddrank1,fortran_int* indata,dim_t mc66,fortran_int *ptx,struct timeval begintime,MPI_Fint fcomm,int solmethod,int adaptive,double epstol,double retryadj,int maxretries,rk_options *rko,solve_real **xcf2,solve_real **accmetric2) {
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
    scheme=(rk_scheme) {"RK2",2,2,0,rk2_c,rk2_a,rk2_b,NULL};
    break;
  case SM_HEUN:
    scheme=(rk_scheme) {"Heun",2,2,0,heun_c,heun_a,heun_b,NULL};
    break;
  case SM_RK4:
    scheme=(rk_scheme) {"RK4",4,4,0,rk4_c,rk4_a,rk4_b,NULL};
    break;
  case SM_BOSHA32:
    scheme=(rk_scheme) {"BoSha32",4,3,1,bosha_c,bosha_a,bosha_b,bosha_bh};
    break;
  default:
    scheme=(rk_scheme) {"DoPri54",7,5,1,dopri_c,dopri_a,dopri_b,dopri_bh};
    break;
  }
  int s=scheme.s;
  bool embedded=(scheme.bhat!=NULL);
  bool chart_log=(rko->chart==RK_CHART_LOG);
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
  /* per-step state: Z = chart coordinate per element (log level ratio,
     initial-based percent, or ordinary change), the snapshot of every
     value at the step base, the s stage gradients in Z-space; stagez
     doubles as the stage/end-state scratch (the increment is formed in
     place) */
  solve_real *base_Z= (solve_real *) calloc (nvarele,sizeof(solve_real));
  solve_real *stagez= (solve_real *) calloc (nvarele,sizeof(solve_real));
  store_real *base_vals= (store_real *) calloc (ncofele+nvarele,sizeof(store_real));
  char *mult= (char *) calloc (nvarele,sizeof(char));
  /* elements that steer the accept test: percent-change variables
     (dimensionless metric); ordinary-change variables carry their own
     units, so at the start of the path max(1,|X|) reads their
     estimate in $ terms -- welfare-decomposition accumulators then
     hijack the step (measured on GTAP-RE); -rkscope all restores the
     GEMPACK rule */
  char *steer= (char *) calloc (nvarele,sizeof(char));
  offset_t nsteer=0;
  solve_real *kk[7];
  for(m=0; m<s; m++)kk[m]= (solve_real *) calloc (nvarele,sizeof(solve_real));
  solve_real *dcum=NULL,*errz=NULL;
  if(embedded) {
    dcum= (solve_real *) calloc (nvarele,sizeof(solve_real));
    errz= (solve_real *) calloc (nvarele,sizeof(solve_real));
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
  offset_t nmult=0;
  for(i=0; i<nvar; i++) {
    for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
      elem_vals1[tindx1].initial=elem_vals1[tindx1].value;
      base_vals[ncofele+tindx1]=elem_vals1[tindx1].value;
      base_Z[tindx1]=0;
      stagez[tindx1]=0;
      /* every ENDOGENOUS percent-change element rides the log chart:
         Z = log(1+X/100) needs no level value (a percent variable's
         value is the level only when a levels form declares one).
         Exogenous elements stay additive: their level-linear path is
         X = S*t exactly in the percent chart (constant gradient h*S),
         whereas integrating log(1+S*t/100) would give the shock itself
         a truncation error. */
      if(chart_log&&!vars[i].change_real&&!CL_EXO(tindx1)) {
        mult[tindx1]=1;
        nmult++;
      }
      if(!vars[i].change_real) {
        steer[tindx1]=1;
        nsteer++;
      }
    }
  }
  memset(&teems_rk_stats,0,sizeof(teems_rk_stats));
  teems_rk_stats.chart=rko->chart;
  teems_rk_stats.h_min=1;
  double t=0,h=1.0/steps1,q;
  if(rko->h0>0&&rko->h0<h)h=rko->h0;
  int stepno=0,checkretries=0;
  bool firstsolve=true,crossed_warned=false,have_k0=false,checks_suspended=false;
  double k0_h=0,metric,metric_prev=-1,stepdata[2];
  if(rank==0)logmsg(1,"%s: %s chart (%ld of %ld elements multiplicative), h0 %.4g, %s norm over %s (%ld elements), %s controller\n",scheme.name,chart_log?"log":"percent",(long)nmult,(long)nvarele,h,rko->norm==RK_NORM_RMS?"rms":"max",rko->scope==RK_SCOPE_PCT?"percent-change elements":"all elements",(long)(rko->scope==RK_SCOPE_PCT?nsteer:nvarele),rko->ctrl==RK_CTRL_PI?"PI":"standard");
  while(t<1.0-1e-12) {
    if(t+h>1.0)h=1.0-t;   /* land exactly on t=1 */
    int failcode=0;
    fortran_int worst_t=0;
    for(m=0; m<s; m++) {
      if(m>0) {
        for(tindx1=0; tindx1<nvarele; tindx1++) {
          solve_real d=0;
          for(st=0; st<m; st++)d+=scheme.a[m*s+st]*kk[st][tindx1];
          stagez[tindx1]=d;
        }
        rk_state_set(vars,nvar,elem_vals,ncofele,nvarele,base_vals,base_Z,stagez,stagez,mult);
        teems_check_viol_range=0;
        teems_check_viol_assert=0;
        teems_rk_stage_checks=(adaptive==1&&!checks_suspended)?1:0;
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
          failcode=rk_state_check(vars,nvar,stagez,mult,rko->guard);
        }
        teems_rk_stage_checks=0;
        MPI_Bcast(&failcode,1,MPI_INT,0,PETSC_COMM_WORLD);
        if(failcode&&adaptive==1)break;   /* abandon the remaining stages */
        if(failcode==RK_FAIL_CROSSED&&!crossed_warned) {
          if(rank==0)printf("Warning: a percent-change variable crossed -100%% during the %s integration; results past this point are unreliable (use more steps or an adaptive method)\n",scheme.name);
          crossed_warned=true;
        }
        failcode=0;
      }
      if(m==0&&have_k0) {
        /* stage 0 gradient known at this base (first-same-as-last after
           an accepted step, or the previous attempt of a retried step):
           rescale to the current step length instead of solving */
        if(h!=k0_h) {
          double r=h/k0_h;
          for(tindx1=0; tindx1<nvarele; tindx1++)kk[0][tindx1]*=r;
          k0_h=h;
        }
        teems_rk_stats.stage_solves_reused++;
        continue;
      }
      /* stage shock vector: the full-substep exogenous change evaluated
         at the stage state.  The exogenous path is uniform in levels:
         the percent fed to the solver rescales with the stage level. */
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
            if(CL_EXO(tindx1)) {
              solve_real ez=h*CL_SHOCK(tindx1);
              VecSetValue(vece,closure_vals[tindx1].exo_index,ez,INSERT_VALUES);
              if(exo_z!=NULL)exo_z[tindx1]=ez;
            }
          }
        }
        else {
          for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
            if(CL_EXO(tindx1)) {
              solve_real ez=mult[tindx1]?h*CL_SHOCK(tindx1)/exp(stagez[tindx1])
                            :h*CL_SHOCK(tindx1)/(1+stagez[tindx1]/100);
              VecSetValue(vece,closure_vals[tindx1].exo_index,ez,INSERT_VALUES);
              if(exo_z!=NULL)exo_z[tindx1]=ez;
            }
          }
        }
      }
      MPI_Barrier(PETSC_COMM_WORLD);
      ierr = VecAssemblyBegin(vece);CHKERRQ(ierr);
      ierr = VecAssemblyEnd(vece);CHKERRQ(ierr);
      /* a singular stage state may fail softly on the sequential LU
         path (the kernel returns instead of aborting); the bordered
         methods' kernels still abort */
      teems_stage_solve_failed=0;
      teems_rk_softfail=(adaptive==1&&matsol==MM_LU)?1:0;
      rk_stage_solve(nohsl,VecSize,BSize,dnz,dnnz,onz,onnz,dnzB,dnnzB,onzB,onnzB,rank,rank_hsl,mpisize,tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele,nvarele,closure_vals,alltimeset,allregset,nintraeq,matsol,Istart,Iend,nreg,ntime,eq_addr,ndblock,countvarintra1,counteq,counteqnoadd,countvarintra1s,counteqs,counteqnoadds,laA,laDi,laD,cntl3,cntl6,nesteddbbd,localsize,ndbbddrank1,indata,mc66,ptx,fcomm,&vece,x1);
      teems_rk_softfail=0;
      teems_rk_stats.stage_solves++;
      {
        int sf=teems_stage_solve_failed;
        MPI_Bcast(&sf,1,MPI_INT,0,PETSC_COMM_WORLD);
        if(sf) {
          failcode=RK_FAIL_SINGULAR;
          break;
        }
      }
      /* recover the backsolved elements from their defining equations
         with this stage's solution (GEMPACK 14.1.3) */
      if(rank==rank_hsl&&nbselems>0)backsolve_recover(tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele,closure_vals,x1,exo_z,bsvals);
      /* stage gradient in Z-space: log chart V/100; percent chart a
         percent solved relative to the stage state scales by
         (1+X_stage/100) (RK.pdf eq. 7); exogenous: h*S in X, and the
         level-linear path's log rate in the log chart */
      for(i=0; i<nvar; i++) {
        if(vars[i].change_real) {
          for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
            if(CL_EXO(tindx1))kk[m][tindx1]=h*CL_SHOCK(tindx1);
            else if(CL_BS(tindx1))kk[m][tindx1]=bsvals[closure_vals[tindx1].exo_index];
            else kk[m][tindx1]=x1[closure_vals[tindx1].exo_index];
          }
        }
        else {
          for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
            solve_real v;
            if(CL_EXO(tindx1))v=mult[tindx1]?h*CL_SHOCK(tindx1)/exp(stagez[tindx1]):h*CL_SHOCK(tindx1);
            else if(CL_BS(tindx1))v=mult[tindx1]?bsvals[closure_vals[tindx1].exo_index]:bsvals[closure_vals[tindx1].exo_index]*(1+stagez[tindx1]/100);
            else v=mult[tindx1]?x1[closure_vals[tindx1].exo_index]:x1[closure_vals[tindx1].exo_index]*(1+stagez[tindx1]/100);
            kk[m][tindx1]=mult[tindx1]?v/100:v;
          }
        }
      }
      if(m==0) {
        have_k0=true;
        k0_h=h;
        /* first step of an adaptive run: cap h from the initial
           gradient only against absurd first moves -- a level ratio
           beyond exp(5) in the log chart, or 90 points (the -100% wall)
           in the percent chart.  The error test, not this cap, sizes
           the step: a tighter cap (measured 2026-09-04 at half a level
           move) starts the run at h ~ 0.01-0.04 and the controller's
           2x growth clamp then needs dozens of steps to reach t = 1. */
        if(adaptive&&stepno==0&&checkretries==0&&rko->h0<=0&&t==0) {
          double hmax=h;
          for(i=0; i<nvar; i++) {
            if(vars[i].change_real)continue;
            for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
              double r=fabs(kk[0][tindx1])/h;   /* log rate, or percent per unit t */
              double cap=mult[tindx1]?5.0:90.0;
              if(r>0&&cap/r<hmax)hmax=cap/r;
            }
          }
          MPI_Bcast(&hmax,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
          if(hmax<h) {
            double hn=hmax;
            if(hn<1e-6)hn=1e-6;
            if(rank==0)logmsg(1,"Step 1: initial step size %.4g from the initial gradient (was %.4g)\n",hn,h);
        fflush(stdout);
            for(tindx1=0; tindx1<nvarele; tindx1++)kk[0][tindx1]*=hn/h;
            h=hn;
            k0_h=h;
          }
        }
      }
    }
    if(failcode==0) {
      /* combine: Z_new = Z_base + sum b_j k_j (formed in stagez);
         embedded error eps = sum (b_j - bhat_j) k_j, in Z-space */
      for(tindx1=0; tindx1<nvarele; tindx1++) {
        solve_real acc=0,eacc=0;
        for(st=0; st<s; st++) {
          acc+=scheme.b[st]*kk[st][tindx1];
          if(embedded)eacc+=(scheme.b[st]-scheme.bhat[st])*kk[st][tindx1];
        }
        stagez[tindx1]=base_Z[tindx1]+acc;
        if(embedded)errz[tindx1]=eacc;
      }
      bool crossed=false;
      for(i=0; i<nvar; i++) {
        if(vars[i].change_real)continue;
        for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++)
          if(!mult[tindx1]&&stagez[tindx1]<=-100)crossed=true;
      }
      metric=0;
      if(embedded) {
        /* per-element error metric of this step (RK.pdf eq. 5): the
           estimate in percent units over max(1,|X|); log chart:
           e_pct = 100*exp(Z)*|eps| (first order in eps) */
        double sumsq=0;
        for(tindx1=0; tindx1<nvarele; tindx1++) {
          double e,v;
          if(rko->scope==RK_SCOPE_PCT&&!steer[tindx1])continue;
          if(mult[tindx1]) {
            double ez=exp(stagez[tindx1]);
            v=fabs(100*(ez-1));
            e=100*ez*fabs(errz[tindx1]);
          }
          else {
            v=fabs(stagez[tindx1]);
            e=fabs(errz[tindx1]);
          }
          e/=((v>1)?v:1);
          if(e>metric) {
            metric=e;
            worst_t=tindx1;
          }
          sumsq+=e*e;
        }
        if(rko->norm==RK_NORM_RMS)metric=sqrt(sumsq/(double)(rko->scope==RK_SCOPE_PCT?nsteer:nvarele));
      }
      /* the accept/reject branch must be collective: every rank follows
         rank 0's metric (Allreduce-equivalent; all ranks hold the full
         solution, HSL worker ranks hold stale copies) */
      stepdata[0]=metric;
      stepdata[1]=crossed?1:0;
      MPI_Bcast(stepdata,2,MPI_DOUBLE,0,PETSC_COMM_WORLD);
      metric=stepdata[0];
      crossed=(stepdata[1]!=0);
      if(crossed&&adaptive==1)failcode=RK_FAIL_CROSSED;
      else if(adaptive&&embedded&&metric>epstol) {
        teems_rk_stats.rejects_accuracy++;
        q=0.85*pow(epstol/metric,1.0/scheme.order);
        if(q<0.5)q=0.5;
        h*=q;
        if(h<1e-8) {
          if(rank==0)printf("Error: adaptive step size collapsed below 1e-8 without meeting epsTolerance %.3g (worst error metric %.3g); the problem may be too stiff for %s\n",epstol,metric,scheme.name);
          MPI_Abort(PETSC_COMM_WORLD,1);
        }
        if(rank==0)logmsg(1,"Step %d rejected (accuracy metric %.3g > %.3g, worst element %s[%ld]), retrying with step size %.4g\n",stepno+1,metric,epstol,rk_elem_var(vars,nvar,worst_t),(long)worst_t,h);
        fflush(stdout);
        for(i=0; i<ncofele; i++)elem_vals[i].value=base_vals[i];
        for(tindx1=0; tindx1<nvarele; tindx1++)stagez[tindx1]=base_Z[tindx1];
        continue;
      }
      else if(crossed&&!crossed_warned) {
        if(rank==0)printf("Warning: a percent-change variable crossed -100%% during the %s integration; results past this point are unreliable (use more steps or an adaptive method)\n",scheme.name);
        crossed_warned=true;
      }
    }
    if(failcode==0) {
      /* realize the end state and run its checks before committing */
      for(tindx1=0; tindx1<nvarele; tindx1++)stagez[tindx1]-=base_Z[tindx1];
      rk_state_set(vars,nvar,elem_vals,ncofele,nvarele,base_vals,base_Z,stagez,stagez,mult);
      teems_check_viol_range=0;
      teems_check_viol_assert=0;
      teems_rk_stage_checks=(adaptive==1&&!checks_suspended)?1:0;
      if(rank==rank_hsl) {
        updates_apply_product(tabfile,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele);
        strcpy(commsyntax,"formula");
        IsIni=false;
        formulas_execute(tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,IsIni);
        assertions_execute(tabfile,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,IsIni,teems_assertions_mode,0);
        failcode=rk_state_check(vars,nvar,stagez,mult,rko->guard);
        if(failcode==RK_FAIL_CROSSED)failcode=0;   /* handled above */
      }
      teems_rk_stage_checks=0;
      MPI_Bcast(&failcode,1,MPI_INT,0,PETSC_COMM_WORLD);
      if(failcode&&adaptive!=1)failcode=0;   /* fixed-step: the checks warned or aborted already */
    }
    if(failcode) {
      /* manual 26.5.1: redo the step at a reduced size; the stage-0
         gradient survives (same base) */
      switch(failcode) {
      case RK_FAIL_CROSSED:
        teems_rk_stats.rejects_crossed++;
        break;
      case RK_FAIL_RANGE:
        teems_rk_stats.rejects_range++;
        break;
      case RK_FAIL_ASSERT:
        teems_rk_stats.rejects_assert++;
        break;
      case RK_FAIL_SINGULAR:
        teems_rk_stats.rejects_singular++;
        break;
      default:
        teems_rk_stats.rejects_guard++;
        break;
      }
      checkretries++;
      if(checkretries>maxretries) {
        bool fatal=(failcode==RK_FAIL_CROSSED||failcode==RK_FAIL_GUARD||failcode==RK_FAIL_SINGULAR||
                    (failcode==RK_FAIL_RANGE&&teems_range_test_updated==2)||
                    (failcode==RK_FAIL_ASSERT&&teems_assertions_mode==2));
        if(fatal) {
          if(rank==0)printf("Error: %s and %d retries at reduced step size did not clear it; the simulation cannot continue (try more steps or a smaller shock)\n",rk_fail_name(failcode),maxretries);
          MPI_Abort(PETSC_COMM_WORLD,1);
        }
        /* warn-level range/assertion violations: the retries are
           exhausted, so take the step with the warnings, as a
           fixed-step run would */
        if(rank==0)printf("Warning: %s and %d retries at reduced step size did not clear it; step %d proceeds with the warning\n",rk_fail_name(failcode),maxretries,stepno+1);
        checks_suspended=true;
      }
      else h*=retryadj;
      if(rank==0)logmsg(1,"Step %d: %s, retrying with step size %.4g\n",stepno+1,rk_fail_name(failcode),h);
        fflush(stdout);
      /* redo the step from its base: stage 0 reads only the
         coefficient values and stagez, so those are all that need
         restoring */
      for(i=0; i<ncofele; i++)elem_vals[i].value=base_vals[i];
      for(tindx1=0; tindx1<nvarele; tindx1++)stagez[tindx1]=base_Z[tindx1];
      continue;
    }
    checkretries=0;
    checks_suspended=false;
    /* accept: propagate the estimated error (RK.pdf eq. 6/8): additive
       in the log chart (relative errors compound), the (1+V/100)
       carry in the percent chart, plain sum for changes */
    if(embedded) {
      for(i=0; i<nvar; i++) {
        if(vars[i].change_real) {
          for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++)
            dcum[tindx1]+=fabs(errz[tindx1]);
        }
        else {
          for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
            if(mult[tindx1])dcum[tindx1]+=fabs(errz[tindx1]);
            else {
              solve_real vstep=(stagez[tindx1]-base_Z[tindx1])/(1+base_Z[tindx1]/100);
              dcum[tindx1]=(1+vstep/100)*dcum[tindx1]+fabs(errz[tindx1]);
            }
          }
        }
      }
    }
    /* commit the resident state as the next base */
    for(i=0; i<ncofele; i++)base_vals[i]=elem_vals[i].value;
    for(i=0; i<nvar; i++) {
      for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
        base_vals[ncofele+tindx1]=elem_vals1[tindx1].value;
        base_Z[tindx1]=stagez[tindx1];
      }
    }
    t+=h;
    stepno++;
    teems_rk_stats.steps++;
    if(h<teems_rk_stats.h_min)teems_rk_stats.h_min=h;
    if(h>teems_rk_stats.h_max)teems_rk_stats.h_max=h;
    if(metric>teems_rk_stats.worst_step_metric)teems_rk_stats.worst_step_metric=metric;
    if(rank==0)logmsg(1,"Step %d accepted: t %.5f, step size %.4g, accuracy metric %.3g (worst element %s[%ld])\n",stepno,t,h,metric,rk_elem_var(vars,nvar,worst_t),(long)worst_t);
    teems_rss_probe("step");
        fflush(stdout);
    /* first-same-as-last: the last stage was solved at the end state =
       the new base, so its gradient is the next stage 0 */
    if(scheme.fsal) {
      memcpy(kk[0],kk[s-1],nvarele*sizeof(solve_real));
      have_k0=true;
      k0_h=h;
    }
    else have_k0=false;
    if(adaptive&&embedded&&t<1.0-1e-12) {
      double mm=(metric>1e-300)?metric:1e-300;
      if(rko->ctrl==RK_CTRL_PI&&metric_prev>0) {
        /* PI control (Gustafsson): damp the growth with the previous
           step's estimate; wider clamp than the elementary rule */
        q=0.85*pow(epstol/mm,0.7/scheme.order)*pow(metric_prev/epstol,0.4/scheme.order);
        if(q<0.2)q=0.2;
        if(q>5.0)q=5.0;
      }
      else {
        q=0.85*pow(epstol/mm,1.0/scheme.order);
        if(q<0.5)q=0.5;
        if(q>2.0)q=2.0;
      }
      h*=q;
      metric_prev=mm;
    }
  }
  teems_rk_stage_checks=0;
  /* the cumulative solution is the solution: no extrapolation pass.
     Report in the solver's percent / change units. */
  *xcf2=(solve_real*)realloc (*xcf2,nvarele*sizeof(solve_real));
  xcf=*xcf2;
  for(tindx1=0; tindx1<nvarele; tindx1++)xcf[tindx1]=mult[tindx1]?100*(exp(base_Z[tindx1])-1):base_Z[tindx1];
  if(embedded) {
    /* estimated error metric per component, E = |Dcum/max(1,|X|)|
       (RK.pdf eq. 5; log chart: Dcum in percent units is
       100*exp(Z)*Dcum_log), written to the .est ride-along by main.
       This is the embedded pair's accumulated per-step estimate -- an
       indicator of the least-settled elements, NOT a bound (measured
       2026-09-03: it under-reads the true error on 30-46% of elements
       at scale). The face value is GEMPACK's 0.02-bracket summary of
       the worst estimate. */
    *accmetric2=(solve_real*)realloc (*accmetric2,nvarele*sizeof(solve_real));
    solve_real *accm=*accmetric2;
    double worst=0;
    for(tindx1=0; tindx1<nvarele; tindx1++) {
      double v=fabs(xcf[tindx1]);
      double d=mult[tindx1]?100*exp(base_Z[tindx1])*fabs(dcum[tindx1]):fabs(dcum[tindx1]);
      accm[tindx1]=d/((v>1)?v:1);
      if(accm[tindx1]>worst)worst=accm[tindx1];
    }
    int face=10-(int)floor(worst/0.02);
    if(face<1)face=1;
    if(face>10)face=10;
    teems_rk_stats.worst_est_metric=worst;
    if(rank==0)printf("Estimated accuracy face value %d of 10 (worst estimated accuracy metric %.3g; an indicator, not a bound)\n",face,worst);
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
  free(base_Z);
  free(stagez);
  free(base_vals);
  free(mult);
  free(steer);
  for(m=0; m<s; m++)free(kk[m]);
  free(dcum);
  free(errz);
  free(exo_z);
  free(bsvals);
  free(counteqs);
  free(counteqnoadds);
  free(countvarintra1s);
  gettimeofday(&endtime, NULL);
  if(rank==0)logmsg(1,"%s solve time %.2f s (%d steps, %ld stage solves, %ld reused; rejected %ld accuracy, %ld check)\n",scheme.name,(endtime.tv_sec - begintime.tv_sec)+((double)(endtime.tv_usec - begintime.tv_usec))/ 1000000,stepno,teems_rk_stats.stage_solves,teems_rk_stats.stage_solves_reused,teems_rk_stats.rejects_accuracy,teems_rk_stats.rejects_crossed+teems_rk_stats.rejects_range+teems_rk_stats.rejects_assert+teems_rk_stats.rejects_guard+teems_rk_stats.rejects_singular);
  teems_rss_probe("solve");
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
bool solve_comp_approx(PetscBool nohsl,PetscInt VecSize,PetscInt dnz,PetscInt* dnnz,PetscInt onz,PetscInt* onnz,PetscInt dnzB,PetscInt* dnnzB,PetscInt onzB,PetscInt* onnzB,Vec *vece1,PetscInt rank,PetscInt rank_hsl,PetscInt mpisize,char* tabfile, char *commsyntax,set_def *sets,dim_t nset, set_element *set_elems, array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar, elem_value **elem_vals2,offset_t ncofele,offset_t nvarele,closure_entry **closure_vals2,offset_t alltimeset,offset_t allregset,offset_t nintraeq,dim_t matsol,PetscInt Istart,PetscInt Iend,offset_t nreg, offset_t ntime, PetscInt *eq_addr, offset_t ndblock, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,dim_t laDi,dim_t laD,PetscReal cntl3,PetscReal cntl6,dim_t nesteddbbd,int localsize,PetscInt *ndbbddrank1,fortran_int* indata,dim_t mc66,fortran_int *ptx,struct timeval begintime,MPI_Fint fcomm,int napprox,int redo_steps,double redo_min_frac,solve_real **xcf2) {
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
          if(CL_EXO(tindx1)) {
            solve_real ez=(tindx1==deloff)?1.0:h*CL_SHOCK(tindx1);
            VecSetValue(vece,closure_vals[tindx1].exo_index,ez,INSERT_VALUES);
            if(exo_z!=NULL)exo_z[tindx1]=ez;
          }
        }
      }
      else {
        for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
          if(CL_EXO(tindx1)) {
            solve_real ez=h*CL_SHOCK(tindx1)/(1+stagex[tindx1]/100);
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
          if(CL_EXO(tindx1))dx[tindx1]=(tindx1==deloff)?1.0:h*CL_SHOCK(tindx1);
          else if(CL_BS(tindx1))dx[tindx1]=bsvals[closure_vals[tindx1].exo_index];
          else dx[tindx1]=x1[closure_vals[tindx1].exo_index];
        }
      }
      else {
        for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
          if(CL_EXO(tindx1))dx[tindx1]=h*CL_SHOCK(tindx1);
          else if(CL_BS(tindx1))dx[tindx1]=bsvals[closure_vals[tindx1].exo_index]*(1+stagex[tindx1]/100);
          else dx[tindx1]=x1[closure_vals[tindx1].exo_index]*(1+stagex[tindx1]/100);
        }
      }
    }
    for(tindx1=0; tindx1<nvarele; tindx1++)newx[tindx1]=base_X[tindx1]+dx[tindx1];
    /* trial advance: values, updates, formulas -- then the state check */
    rk_state_set(vars,nvar,elem_vals,ncofele,nvarele,base_vals,base_X,dx,stagex,NULL);
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
        rk_state_set(vars,nvar,elem_vals,ncofele,nvarele,base_vals,base_X,dx,stagex,NULL);
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
  teems_rss_probe("solve");
  return 1;
}
