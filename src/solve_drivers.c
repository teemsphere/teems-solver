#include <teems_solver.h>
#include <hsl_kernels.h>

/* --- MA48 workspace growth ---------------------------------------
   Every -3 (workspace too small) return in the solver funnels through
   ma48_grow_la.  The growth itself is MA48's suggested size with a
   doubling floor (which also repairs a suggestion that came back
   wrapped), but the essential part is the MA48_LA_MAX clamp: the HSL
   build is 32-bit, so once a request passes INT_MAX there is no larger
   workspace to ask for and every further doubling is a lie the caller
   cannot detect.  At the ceiling we stop and say so. */
offset_t ma48_grow_la(offset_t cur,offset_t suggested,offset_t nnz,const char *knob,long *used) {
  offset_t newla=suggested;
  if(newla<2*cur)newla=2*cur;
  if(newla>MA48_LA_MAX)newla=MA48_LA_MAX;
  if(newla<=cur) {
    errmsg("Error: factorizing %s needs a larger MA48 workspace than the %ld-element limit of the 32-bit HSL build (-%s is already at its %ld%% ceiling); this system is too large for one sequential factorization -- solve it with a bordered matrix_method (\"SBBD\" or \"DBBD\"), condense the model, or reduce its dimensions\n",
           probe_onfail_scope_label(),(long)MA48_LA_MAX,knob,
           (long)((100.0*MA48_LA_MAX)/nnz));
    fflush(stdout);
    MPI_Abort(PETSC_COMM_WORLD,1);
  }
  logmsg(1,"Note: MA48 workspace grown from %ld to %ld reals (equivalent -%s %ld)\n",
         (long)cur,(long)newla,knob,(long)ceil((100.0*newla)/nnz));
  {
    long eqpct=(long)ceil((100.0*newla)/nnz);
    #pragma omp critical(laused)
    if(eqpct>*used)*used=eqpct;
  }
  return newla;
}

/* Initial LA from a -la* percent, clamped to the 32-bit ceiling and
   formed in double: the product passes INT_MAX long before the
   percent looks unreasonable, and a wrapped initial LA is exactly as
   silent as a wrapped grown one. */
offset_t ma48_la_from_pct(dim_t pct,offset_t nnz) {
  double la=ceil((pct/100.0)*(double)nnz);
  if(la>(double)MA48_LA_MAX)return MA48_LA_MAX;
  if(la<0)return 0;
  return (offset_t)la;
}

/* A clamped LA still reaches tens of GB; these report the shortfall
   instead of segfaulting on an unchecked NULL. */
static void ma48_alloc_fail(offset_t n,size_t sz) {
  errmsg("Error: cannot allocate the %.1f GB MA48 workspace this factorization needs (%ld elements); free memory on this machine, use a bordered matrix_method (\"SBBD\" or \"DBBD\") to split the factorization, or condense the model\n",
         (n*(double)sz)/1073741824.0,(long)n);
  fflush(stdout);
  MPI_Abort(PETSC_COMM_WORLD,1);
}

void *ma48_alloc(offset_t n,size_t sz) {
  void *p=calloc(n,sz);
  if(p==NULL&&n>0)ma48_alloc_fail(n,sz);
  return p;
}

void *ma48_realloc(void *p,offset_t n,size_t sz) {
  void *q=realloc(p,n*sz);
  if(q==NULL&&n>0)ma48_alloc_fail(n,sz);
  return q;
}

/* --- Jacobian preallocation ceiling --------------------------------
   Mat*AIJSetPreallocation totals the per-row counts in a PetscInt.  For
   every HSL matrix_method A and B are MATSEQAIJ on PETSC_COMM_SELF and
   hold the whole system on rank_hsl, so on a 32-bit PetscInt build the
   total passes the ceiling near 2^31 nonzeros (~450M equations at 4.9
   per row); PETSc then returns an error and, unchecked, the matrix was
   dereferenced at the first insert (SEGV at address 0).  The total is
   formed here in 64-bit first so the message can name the count and
   the way out, and the PETSc return is checked so any other
   preallocation failure is reported instead of segfaulting. */
void jac_mat_prealloc(Mat M,const char *what,PetscBool mpi,int count,PetscInt nrows_local,PetscInt dnz,PetscInt *dnnz,PetscInt onz,PetscInt *onnz) {
  PetscErrorCode ierr;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  if(count&&dnnz!=NULL) {
    offset_t nnz=0;
    PetscInt i;
    for(i=0;i<nrows_local;i++)nnz+=(offset_t)dnnz[i]+((mpi&&onnz!=NULL)?(offset_t)onnz[i]:0);
    if(nnz>(offset_t)PETSC_INT_MAX) {
      if(mpi)
        errmsg("Error: assembling the %s for %s needs %ld nonzeros on rank %d, above the %ld-nonzero ceiling of the %d-bit PetscInt build; spread the system over more MPI ranks, condense the model, or reduce its dimensions\n",
               what,probe_onfail_scope_label(),(long)nnz,rank,(long)PETSC_INT_MAX,(int)(8*sizeof(PetscInt)));
      else
        errmsg("Error: assembling the %s for %s needs %ld nonzeros, above the %ld-nonzero ceiling of the %d-bit PetscInt build; every HSL matrix_method (\"LU\", \"SBBD\") keeps the whole system on one rank, so more ranks do not help -- condense the model, reduce its dimensions, or use the distributed matrix_method \"DBBD\"\n",
               what,probe_onfail_scope_label(),(long)nnz,(long)PETSC_INT_MAX,(int)(8*sizeof(PetscInt)));
      fflush(stdout);
      MPI_Abort(PETSC_COMM_WORLD,1);
    }
  }
  if(mpi)ierr=MatMPIAIJSetPreallocation(M,dnz,dnnz,onz,onnz);
  else ierr=MatSeqAIJSetPreallocation(M,dnz,dnnz);
  if(ierr) {
    errmsg("Error: PETSc could not preallocate the %s for %s (PETSc error %d)\n",what,probe_onfail_scope_label(),(int)ierr);
    fflush(stdout);
    MPI_Abort(PETSC_COMM_WORLD,1);
  }
}

/* Persistent-factor sequential LU (-fastrefac): the Jacobian's stored
   pattern is fixed across steps, so the MA48 pivot sequence is computed
   once and later steps only refactorize (MA48B/BD JOB=2) with fresh
   values.  Extraction keeps explicitly-stored zeros: an entry that is
   zero at the analyse state can become nonzero at a later step and must
   be present in the pattern.  rank_hsl only. */
static int *fr_irn=NULL,*fr_jcn=NULL;
static solve_real *fr_values=NULL;
static PetscInt fr_nz=-1;
static PetscInt *fr_pat=NULL;  /* the assembled column pattern the pivot sequence belongs to
                                  (the staged JCN is clobbered by MA48, so it cannot be
                                  compared): the assembly skips zero-valued entries, so the
                                  pattern is value-dependent and an equal NE does not prove
                                  it unchanged */
static offset_t fr_lasize=0;   /* current LA; grows on MA48 -3 returns and stays grown */
static int fr_ready=0;
/* ---- solve accuracy (manual 30.1.5, 30.6.1) and arithmetic errors
   (34.3) ----------------------------------------------------------
   After a solve A y = b the residual of every equation is compared with
   the sum of the absolute values of its terms (GEMPACK's residual
   ratio); the maximum across the simulation is logged and recorded in
   stats.json, and an equation whose ratio reaches 1e-4 (30.1.5: "somewhat
   concerned ... 0.0001 or larger") is reported as not satisfied very
   accurately. More than 100 such warnings make the run end with an
   error after all files are written, as GEMPACK does (30.6.1). The check
   needs A and b after the solve: the LU paths and -fastrefac SBBD keep
   them; one-shot SBBD, DBBD and NDBBD release A before or during the
   factorization, so their solves are counted as skipped unless A is
   kept (-residcheck 1, or the DBBD refinement step, which checks the
   refined solution). Every method's solution is scanned for NaN/Inf. */
double teems_resid_max=0.0;
long teems_resid_solves=0,teems_resid_warn=0,teems_resid_skipped=0;
#define TEEMS_RESID_WARN 1e-4
#define TEEMS_RESID_WARN_PRINT 10
/* A row whose terms all sit at rounding level (e.g. ps - po = 0 with
   both near 1e-14) has a residual as large as its terms and so a ratio
   of 1 whatever the solve's accuracy. The ratio's denominator is
   therefore at least TEEMS_RESID_FLOOR x the row's largest coefficient
   x the solution's largest magnitude: a residual below 1e-12 of that
   normwise scale (TEEMS_RESID_WARN x the floor) never warns. */
#define TEEMS_RESID_FLOOR 1e-8

static double resid_xmax(const solve_real *x, PetscInt n) {
  PetscInt i;
  double m=0.0;
  for(i=0; i<n; i++) if(fabs((double)x[i])>m) m=fabs((double)x[i]);
  return m;
}

void solve_x_check(const solve_real *x, PetscInt n, int doit) {
  PetscInt i;
  if(!doit||x==NULL) return;
  for(i=0; i<n; i++) if(teems_nonfinite((double)x[i])) {
      char lab[(MAXVARDIM+2)*NAMESIZE];
      probe_col_label(i,lab,sizeof(lab));
      errmsg("Error: the linear solve gave a value that is not finite (%s) for %s: the LHS matrix is singular or badly scaled at this step, or a coefficient overflowed (manual 34.1, 34.3)\n",teems_isnan_bits((double)x[i])?"NaN":"infinite",lab);
      MPI_Abort(PETSC_COMM_WORLD,1);
    }
}

static void solve_residual_check(Mat A, const PetscScalar *b, const solve_real *x, PetscInt n) {
  Mat_SeqAIJ *aa=(Mat_SeqAIJ*)A->data;
  PetscInt i,k;
  double worst=0.0,xmax;
  long nw=0;
  solve_x_check(x,n,1);
  xmax=resid_xmax(x,n);
  for(i=0; i<n; i++) {
    double sum=0.0,sabs=fabs((double)b[i]),r,q,amax=0.0,den;
    for(k=aa->i[i]; k<aa->i[i+1]; k++) {
      double t=(double)aa->a[k]*(double)x[aa->j[k]];
      sum+=t;
      sabs+=fabs(t);
      if(fabs((double)aa->a[k])>amax) amax=fabs((double)aa->a[k]);
    }
    r=fabs(sum-(double)b[i]);
    den=TEEMS_RESID_FLOOR*amax*xmax;
    if(sabs>den) den=sabs;
    q=(den>0.0)?r/den:0.0;
    if(q>worst) worst=q;
    if(q>=TEEMS_RESID_WARN) {
      nw++;
      if(nw<=TEEMS_RESID_WARN_PRINT) {
        char lab[(4*MAXVARDIM+2)*NAMESIZE];
        probe_row_label(i,lab,sizeof(lab));
        printf("Warning: equation %s is not satisfied very accurately (residual ratio %.3e: the sum of its terms is %.6g while the sum of their absolute values is %.6g); this may be because the LHS matrix is not really invertible (manual 30.6.1)\n",lab,q,sum-(double)b[i],sabs);
      }
    }
  }
  if(nw>TEEMS_RESID_WARN_PRINT) printf("Warning: ... %ld more equations of this solve are not satisfied very accurately (manual 30.6.1)\n",nw-TEEMS_RESID_WARN_PRINT);
  if(worst>teems_resid_max) teems_resid_max=worst;
  teems_resid_solves++;
  teems_resid_warn+=nw;
}

/* ---- factor once, solve many (Tier B3; API in teems_solver.h) ----
   What a method keeps when teems_fh_keep is set for its dispatch:
     LU one-shot     MA48 L/U (VA, IRN) + KEEP/controls in Fortran slot 1
     LU -fastrefac   nothing extra: the persistent factors are the handle
     SBBD one-shot   the MP48 instance (not torn down after its solve)
     SBBD -fastrefac nothing extra: the persistent instance
     DBBD            block factors, border blocks, interface factors
                     (block_solve.c, dfh)
   teems_resid_retain keeps the LHS matrix (a PETSc reference) and the
   RHS for the residual check of methods that release A before their
   solve completes (one-shot SBBD, DBBD, NDBBD). */
int teems_fh_selftest=0,teems_resid_all=0,teems_fh_keep=0,teems_resid_retain=0;
int teems_refine=1;
static int fh_refine_now=0;
static long fh_rf_solves=0;
static double fh_rf_before=0.0,fh_rf_after=0.0,fh_rf_secs=0.0;
static void fh_refine_cols(const solve_real *bloc,solve_real *x,int nc,PetscInt rank,PetscInt mpisize,double *before);
enum { FH_NONE=0, FH_LU, FH_LU_FR, FH_SBBD, FH_SBBD_FR, FH_DBBD };
static int fh_kind=FH_NONE;
static PetscInt fh_n=0;
static int fh_mpi=0,fh_resid_pending=0;
static int *fh_lu_irn=NULL;
static solve_real *fh_lu_va=NULL;
static fortran_int fh_sb_indata[5];
static Mat fh_A=NULL;
static solve_real *fh_b=NULL,*fh_bfull=NULL; /* RHS: this rank's rows (MPI) / whole on rank_hsl (seq); whole on rank 0 (self-test) */
static PetscInt fh_bstart=0,fh_bend=0;
static long fh_st_solves=0,fh_st_rhs=0,fh_st_bitdiff=0;
static double fh_st_maxabs=0.0,fh_st_maxrel=0.0;

void residual_note_skipped(PetscInt rank,PetscInt counting_rank) {
  if(teems_resid_retain)fh_resid_pending=1;
  else if(rank==counting_rank)teems_resid_skipped++;
}

void fh_request_check(dim_t matsol,dim_t mc66,PetscInt rank) {
  const char *who=NULL;
  if(!teems_fh_selftest)return;
  if(matsol==MM_NDBBD)who="matrix_method NDBBD";
  else if(matsol==MM_SBBD&&mc66!=0)who="SBBD with -withmc66 1";
  if(who==NULL)return;
  if(rank==0)errmsg("Error: this run solves extra right-hand sides with each step's factorization (-fhtest), which %s cannot keep yet; use matrix_method LU, SBBD or DBBD\n",who);
  MPI_Barrier(PETSC_COMM_WORLD);
  MPI_Abort(PETSC_COMM_WORLD,1);
}

/* Residual check (manual 30.1.5, 30.6.1) of a distributed solve: each
   rank checks its own rows of A against the whole solution x (every rank
   holds it after the bordered methods' final reduction); rank 0 prints
   the first warnings in row-ownership order and records the solve. */
static double solve_residual_check_mpi(Mat A,const solve_real *bloc,PetscInt bstart,PetscInt bend,const solve_real *x,PetscInt rank,PetscInt mpisize) {
  PetscInt rs,re,i,k,nc;
  const PetscInt *cols;
  const PetscScalar *va;
  double worst=0.0,gworst=0.0,xmax;
  long nw=0,gnw=0;
  double wbuf[4*TEEMS_RESID_WARN_PRINT];
  int nwp=0,*cnt=NULL,*disp=NULL,tot=0;
  double *all=NULL;
  MatGetOwnershipRange(A,&rs,&re);
  xmax=resid_xmax(x,A->cmap->N);
  for(i=rs; i<re; i++) {
    double bi=(i>=bstart&&i<bend)?(double)bloc[i-bstart]:0.0;
    double sum=0.0,sabs=fabs(bi),r,q,amax=0.0,den;
    MatGetRow(A,i,&nc,&cols,&va);
    for(k=0; k<nc; k++) {
      double t=(double)va[k]*(double)x[cols[k]];
      sum+=t;
      sabs+=fabs(t);
      if(fabs((double)va[k])>amax) amax=fabs((double)va[k]);
    }
    MatRestoreRow(A,i,&nc,&cols,&va);
    r=fabs(sum-bi);
    den=TEEMS_RESID_FLOOR*amax*xmax;
    if(sabs>den) den=sabs;
    q=(den>0.0)?r/den:0.0;
    if(q>worst)worst=q;
    if(q>=TEEMS_RESID_WARN) {
      nw++;
      if(nwp<TEEMS_RESID_WARN_PRINT) {
        wbuf[4*nwp]=(double)i;
        wbuf[4*nwp+1]=q;
        wbuf[4*nwp+2]=sum-bi;
        wbuf[4*nwp+3]=sabs;
        nwp++;
      }
    }
  }
  MPI_Reduce(&worst,&gworst,1,MPI_DOUBLE,MPI_MAX,0,PETSC_COMM_WORLD);
  MPI_Reduce(&nw,&gnw,1,MPI_LONG,MPI_SUM,0,PETSC_COMM_WORLD);
  if(rank==0) {
    cnt=(int *) calloc (mpisize,sizeof(int));
    disp=(int *) calloc (mpisize,sizeof(int));
  }
  nwp*=4;
  MPI_Gather(&nwp,1,MPI_INT,cnt,1,MPI_INT,0,PETSC_COMM_WORLD);
  if(rank==0) {
    for(i=0; i<mpisize; i++) {
      disp[i]=tot;
      tot+=cnt[i];
    }
    all=(double *) malloc ((tot>0?tot:1)*sizeof(double));
  }
  MPI_Gatherv(wbuf,nwp,MPI_DOUBLE,all,cnt,disp,MPI_DOUBLE,0,PETSC_COMM_WORLD);
  if(rank==0) {
    for(k=0; k<tot/4&&k<TEEMS_RESID_WARN_PRINT; k++) {
      char lab[(4*MAXVARDIM+2)*NAMESIZE];
      probe_row_label((PetscInt)all[4*k],lab,sizeof(lab));
      printf("Warning: equation %s is not satisfied very accurately (residual ratio %.3e: the sum of its terms is %.6g while the sum of their absolute values is %.6g); this may be because the LHS matrix is not really invertible (manual 30.6.1)\n",lab,all[4*k+1],all[4*k+2],all[4*k+3]);
    }
    if(gnw>TEEMS_RESID_WARN_PRINT)printf("Warning: ... %ld more equations of this solve are not satisfied very accurately (manual 30.6.1)\n",gnw-TEEMS_RESID_WARN_PRINT);
    if(gworst>teems_resid_max)teems_resid_max=gworst;
    teems_resid_solves++;
    teems_resid_warn+=gnw;
    free(cnt);
    free(disp);
    free(all);
  }
  return gworst;
}

/* ---- shock-group subtotals (Tier B3; GEMPACK manual 29; Harrison,
   Horridge and Pearson 1999, CoPS IP-73 sections 3.6, 4.2) ----
   At every solve the right-hand side B*dz splits by group: dz_r is the
   step's shock vector with every entry outside group r zeroed, and
   A p_r = B dz_r is solved with the step's kept factorization (the q_r
   sum to the step's RHS, so the p_r sum to its solution). Each group
   carries shadow copies of the driver's cumulative state and every
   update the driver makes to its own state is made to the group's with
   the same weights: ordinary changes add; a percent-change variable's
   step solution is weighted by the level at the start of the step
   relative to the pre-simulation level (1 + p/100, the MAIN solution's
   cumulative value; IP-73 eq. 7), the Gragg leapfrog, terminal
   smoothing, Richardson weights and subinterval compounding follow the
   driver's formulas. The columns therefore add up to the solution to
   rounding for any partition of the shocks, for Johansen, Euler and
   Gragg alike; how the total splits depends on the path (IP-73 3.4),
   which is the solver's straight line. Backsolved elements are
   recovered per group from p_r and the group's exogenous changes.
   Everything accumulates on rank 0. */
int teems_sub_active=0;
#define SUB_CHUNK 8
static solve_real *sub_rhs=NULL,*sub_x=NULL,*sub_bs=NULL,*sub_z=NULL;
static solve_real **sub_c=NULL,**sub_l=NULL,**sub_f=NULL;
static PetscInt sub_n=0;

/* the group right-hand sides B*dz_r, right after the step's own
   MatMult(B,vece,vecb) and before B is released; collective under the
   distributed methods, rank_hsl only otherwise */
void sub_step_rhs(Mat B,Vec vece,PetscInt VecSize,dim_t matsol,PetscInt rank,PetscInt rank_hsl) {
  int mpi=(matsol>=MM_DBBD),r;
  Vec tmp,out;
  PetscInt lo,hi;
  offset_t k;
  if(!teems_sub_active)return;
  if(!mpi&&rank!=rank_hsl)return;
  sub_n=VecSize;
  if(rank==0&&sub_rhs==NULL)sub_rhs=(solve_real *) malloc ((size_t)VecSize*teems_nsub*sizeof(solve_real));
  VecDuplicate(vece,&tmp);
  MatCreateVecs(B,NULL,&out);
  VecGetOwnershipRange(vece,&lo,&hi);
  for(r=0; r<teems_nsub; r++) {
    const PetscScalar *ev,*ov;
    PetscScalar *tv;
    VecSet(tmp,0.0);
    VecGetArrayRead(vece,&ev);
    VecGetArray(tmp,&tv);
    for(k=0; k<teems_subs[r].nmem; k++) {
      PetscInt e=(PetscInt)teems_subs[r].exo[k];
      if(e>=lo&&e<hi)tv[e-lo]=ev[e-lo];
    }
    VecRestoreArray(tmp,&tv);
    VecRestoreArrayRead(vece,&ev);
    MatMult(B,tmp,out);
    if(!mpi) {
      VecGetArrayRead(out,&ov);
      memcpy(sub_rhs+(size_t)r*VecSize,ov,VecSize*sizeof(solve_real));
      VecRestoreArrayRead(out,&ov);
    }
    else {
      VecScatter sct;
      Vec vz;
      VecScatterCreateToZero(out,&sct,&vz);
      VecScatterBegin(sct,out,vz,INSERT_VALUES,SCATTER_FORWARD);
      VecScatterEnd(sct,out,vz,INSERT_VALUES,SCATTER_FORWARD);
      if(rank==0) {
        VecGetArrayRead(vz,&ov);
        memcpy(sub_rhs+(size_t)r*VecSize,ov,VecSize*sizeof(solve_real));
        VecRestoreArrayRead(vz,&ov);
      }
      VecScatterDestroy(&sct);
      VecDestroy(&vz);
    }
  }
  VecDestroy(&tmp);
  VecDestroy(&out);
}

/* rank 0's full columns (n x nc, column-major) to this rank's rows of
   the kept A (nl x nc); collective */
static solve_real *fh_rows_scatter(const solve_real *full0,int nc,PetscInt rank,PetscInt mpisize) {
  const PetscInt *rng;
  PetscInt rs,re,n=fh_n;
  int *cnt=NULL,*disp=NULL,k,j;
  solve_real *loc;
  MatGetOwnershipRange(fh_A,&rs,&re);
  loc=(solve_real *) malloc ((re>rs?(size_t)(re-rs)*nc:1)*sizeof(solve_real));
  if(rank==0) {
    MatGetOwnershipRanges(fh_A,&rng);
    cnt=(int *) calloc (mpisize,sizeof(int));
    disp=(int *) calloc (mpisize,sizeof(int));
    for(j=0; j<mpisize; j++) {
      cnt[j]=(int)(rng[j+1]-rng[j]);
      disp[j]=(int)rng[j];
    }
  }
  for(k=0; k<nc; k++)MPI_Scatterv(rank==0?full0+(size_t)k*n:NULL,cnt,disp,SORD==1?MPI_DOUBLE:MPI_FLOAT,loc+(size_t)k*(re-rs),(int)(re-rs),SORD==1?MPI_DOUBLE:MPI_FLOAT,0,PETSC_COMM_WORLD);
  free(cnt);
  free(disp);
  return loc;
}

/* the group step solutions p_r, SUB_CHUNK columns at a time with the
   kept factorization (collective); each chunk's solutions overwrite its
   right-hand sides, so the step holds one n x G block on rank 0. Under
   -refine (DBBD) every column gets the same refinement step as the
   step's own solution: the step is linear in the right-hand side
   (x + A^-1 (b - A x) with x = A^-1 b), so the refined columns still sum
   to the refined solution, and a group of every shock reproduces it bit
   for bit (same kernels, same row sums) */
static void sub_step_solve(PetscInt rank,PetscInt n) {
  int c0,nc,mpisize;
  MPI_Comm_size(PETSC_COMM_WORLD,&mpisize);
  sub_x=sub_rhs;
  for(c0=0; c0<teems_nsub; c0+=SUB_CHUNK) {
    solve_real *xs;
    nc=teems_nsub-c0<SUB_CHUNK?teems_nsub-c0:SUB_CHUNK;
    xs=(solve_real *) malloc ((size_t)n*nc*sizeof(solve_real));
    if(teems_fh_solve(rank==0?sub_rhs+(size_t)c0*n:NULL,xs,nc)!=0) {
      if(rank==0)errmsg("Error: the subtotal solves found no kept factorization for this step\n");
      MPI_Abort(PETSC_COMM_WORLD,1);
    }
    if(fh_refine_now) {
      solve_real *bl=fh_rows_scatter(rank==0?sub_rhs+(size_t)c0*n:NULL,nc,rank,mpisize);
      fh_refine_cols(bl,xs,nc,rank,mpisize,NULL);
      free(bl);
    }
    if(rank==0)memcpy(sub_x+(size_t)c0*n,xs,(size_t)n*nc*sizeof(solve_real));
    free(xs);
  }
}

static void sub_alloc(offset_t nvarele,int multistep) {
  int r;
  if(sub_f==NULL) {
    sub_f=(solve_real **) calloc (teems_nsub,sizeof(solve_real*));
    for(r=0; r<teems_nsub; r++)sub_f[r]=(solve_real *) calloc (nvarele>0?nvarele:1,sizeof(solve_real));
  }
  if(multistep&&sub_c==NULL) {
    sub_c=(solve_real **) calloc (teems_nsub,sizeof(solve_real*));
    for(r=0; r<teems_nsub; r++)sub_c[r]=(solve_real *) calloc (nvarele>0?nvarele:1,sizeof(solve_real));
  }
  if((multistep==SUB_LEAP||multistep==SUB_SMOOTH)&&sub_l==NULL) {
    sub_l=(solve_real **) calloc (teems_nsub,sizeof(solve_real*));
    for(r=0; r<teems_nsub; r++)sub_l[r]=(solve_real *) calloc (nvarele>0?nvarele:1,sizeof(solve_real));
  }
  if(nbselems>0&&sub_bs==NULL) {
    sub_bs=(solve_real *) calloc (nbselems,sizeof(solve_real));
    sub_z=(solve_real *) calloc (nvarele>0?nvarele:1,sizeof(solve_real));
  }
}

/* group r's recovered backsolved changes: backsolve_recover on p_r and
   the step's exogenous changes restricted to group r */
static void sub_backsolve(int r,char *tabfile,char *commsyntax,set_def *sets,dim_t nset,set_element *set_elems,array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar,elem_value *elem_vals,offset_t ncofele,closure_entry *closure_vals,offset_t nvarele,const solve_real *exo_z) {
  offset_t k;
  if(nbselems<=0)return;
  memset(sub_z,0,nvarele*sizeof(solve_real));
  if(exo_z!=NULL)for(k=0; k<teems_subs[r].nmem; k++)sub_z[teems_subs[r].mem[k]]=exo_z[teems_subs[r].mem[k]];
  backsolve_recover(tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele,closure_vals,sub_x+(size_t)r*sub_n,sub_z,sub_bs);
}

/* mirror of one driver update on every group's state (rank 0). mode:
   SUB_JOHANSEN the one-step result, SUB_FIRST a pass's first step,
   SUB_EULER a forward step, SUB_LEAP a Gragg leapfrog step, SUB_SMOOTH
   Gragg's terminal smoothing. V is the driver's cumulative state before
   its own update (the level weight of percent-change variables). The
   exogenous members are set from the driver's state after its update
   (sub_exo_sync): a shock belongs wholly to its own group. */
void sub_update(int mode,char *tabfile,char *commsyntax,set_def *sets,dim_t nset,set_element *set_elems,array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar,elem_value *elem_vals,offset_t ncofele,closure_entry *closure_vals,offset_t nvarele,const solve_real *exo_z,const solve_real *V) {
  int r;
  offset_t i,j,k;
  sub_alloc(nvarele,mode);
  for(r=0; r<teems_nsub; r++) {
    solve_real *xr=sub_x+(size_t)r*sub_n,*c=(mode==SUB_JOHANSEN)?sub_f[r]:sub_c[r],*l=(sub_l==NULL)?NULL:sub_l[r];
    sub_backsolve(r,tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele,closure_vals,nvarele,exo_z);
    for(i=0; i<nvar; i++) {
      int chg=vars[i].change_real?1:0;
      for(j=vars[i].offset; j<vars[i].nelem+vars[i].offset; j++) {
        solve_real d,t;
        if(CL_EXO(j)) {
          if(mode==SUB_JOHANSEN||mode==SUB_FIRST)c[j]=0;
          continue;
        }
        d=CL_BS(j)?sub_bs[closure_vals[j].exo_index]:xr[closure_vals[j].exo_index];
        switch(mode) {
        case SUB_JOHANSEN:
        case SUB_FIRST:
          c[j]=d;
          if(l!=NULL)l[j]=0;
          break;
        case SUB_EULER:
          if(chg)c[j]+=d;
          else c[j]+=d*(100+V[j])/100;
          break;
        case SUB_LEAP:
          /* the driver's ordinary-change lag is the variable's stored
             value, i.e. the cumulative change at the storage precision
             (store_real); its percent-change lag is the double itself */
          if(chg) {
            t=(solve_real)(store_real)c[j];
            c[j]=l[j]+2*d;
          }
          else {
            t=c[j];
            c[j]=l[j]+2*d*(100+V[j])/100;
          }
          l[j]=t;
          break;
        case SUB_SMOOTH:
          if(chg)c[j]=0.5*(c[j]+l[j]+d);
          else c[j]=0.5*(c[j]+l[j]+d*(1+V[j]/100));
          l[j]=0;
          break;
        }
      }
    }
    if(mode==SUB_JOHANSEN)for(k=0; k<teems_subs[r].nmem; k++)c[teems_subs[r].mem[k]]=CL_SHOCK(teems_subs[r].mem[k]);
  }
}

void sub_exo_sync(const solve_real *varchange) {
  int r;
  offset_t k;
  for(r=0; r<teems_nsub; r++)for(k=0; k<teems_subs[r].nmem; k++)sub_c[r][teems_subs[r].mem[k]]=varchange[teems_subs[r].mem[k]];
}

/* a pass's Richardson contribution and the subinterval compounding, the
   driver's formulas on each group's pass state (xc0 = 1 + the main
   solution's cumulative percent change at the subinterval start) */
void sub_pass_end(int sol,dim_t subindx,array_def *vars,offset_t nvar,offset_t nvarele,const solve_real *xc0) {
  int r;
  offset_t i,k;
  for(r=0; r<teems_nsub; r++) {
    solve_real *c=sub_c[r],*f=sub_f[r];
    if(subindx>0) {
      for(i=0; i<nvar; i++) {
        for(k=vars[i].offset; k<vars[i].nelem+vars[i].offset; k++) {
          if(vars[i].change_real) {
            if(sol==0)f[k]+=c[k]*extrap_w2;
            if(sol==1)f[k]-=c[k]*extrap_w3;
            if(sol==2)f[k]+=c[k]*extrap_w3;
          }
          else {
            if(sol==0)f[k]+=c[k]*xc0[k]*extrap_w2;
            if(sol==1)f[k]-=c[k]*xc0[k]*extrap_w3;
            if(sol==2)f[k]+=c[k]*xc0[k]*extrap_w3;
          }
        }
      }
    }
    else {
      for(k=0; k<nvarele; k++) {
        if(sol==0)f[k]+=c[k]*extrap_w2;
        if(sol==1)f[k]-=c[k]*extrap_w3;
        if(sol==2)f[k]+=c[k]*extrap_w3;
      }
    }
  }
}

solve_real **sub_columns(void) {
  return sub_f;
}

void sub_free(void) {
  int r;
  for(r=0; r<teems_nsub; r++) {
    if(sub_f!=NULL)free(sub_f[r]);
    if(sub_c!=NULL)free(sub_c[r]);
    if(sub_l!=NULL)free(sub_l[r]);
  }
  free(sub_f);
  free(sub_c);
  free(sub_l);
  free(sub_rhs);
  free(sub_bs);
  free(sub_z);
  sub_f=sub_c=sub_l=NULL;
  sub_rhs=sub_x=sub_bs=sub_z=NULL;
}

void fh_step_begin(Mat A,Vec vecb,dim_t matsol,dim_t mc66,PetscInt VecSize,PetscInt rank,PetscInt rank_hsl) {
  fh_request_check(matsol,mc66,rank);
  /* -refine (DBBD): the refinement step needs the factors and A */
  fh_refine_now=(teems_refine&&matsol==MM_DBBD);
  teems_fh_keep=teems_fh_selftest||teems_sub_active||fh_refine_now;
  teems_resid_retain=(teems_resid_all&&matsol>=MM_SBBD)||fh_refine_now;
  fh_resid_pending=0;
  fh_kind=FH_NONE;
  fh_n=VecSize;
  fh_mpi=(matsol>=MM_DBBD);
  if(!teems_fh_keep&&!teems_resid_retain)return;
  /* the sequential methods know the system size on rank_hsl (0) only */
  if(teems_fh_keep)MPI_Bcast(&fh_n,1,MPIU_INT,0,PETSC_COMM_WORLD);
  if(teems_resid_retain) {
    PetscObjectReference((PetscObject)A);
    fh_A=A;
  }
  if(fh_mpi) {
    const PetscScalar *bv;
    VecGetOwnershipRange(vecb,&fh_bstart,&fh_bend);
    if(teems_resid_retain) {
      fh_b=(solve_real *) malloc ((fh_bend>fh_bstart?fh_bend-fh_bstart:1)*sizeof(solve_real));
      VecGetArrayRead(vecb,&bv);
      memcpy(fh_b,bv,(fh_bend-fh_bstart)*sizeof(solve_real));
      VecRestoreArrayRead(vecb,&bv);
    }
    if(teems_fh_selftest) {
      VecScatter sct;
      Vec vz;
      VecScatterCreateToZero(vecb,&sct,&vz);
      VecScatterBegin(sct,vecb,vz,INSERT_VALUES,SCATTER_FORWARD);
      VecScatterEnd(sct,vecb,vz,INSERT_VALUES,SCATTER_FORWARD);
      if(rank==0) {
        fh_bfull=(solve_real *) malloc (VecSize*sizeof(solve_real));
        VecGetArrayRead(vz,&bv);
        memcpy(fh_bfull,bv,VecSize*sizeof(solve_real));
        VecRestoreArrayRead(vz,&bv);
      }
      VecScatterDestroy(&sct);
      VecDestroy(&vz);
    }
  }
  else if(rank==rank_hsl) {
    const PetscScalar *bv;
    fh_b=(solve_real *) malloc (VecSize*sizeof(solve_real));
    VecGetArrayRead(vecb,&bv);
    memcpy(fh_b,bv,VecSize*sizeof(solve_real));
    VecRestoreArrayRead(vecb,&bv);
    fh_bfull=fh_b;
  }
}

int teems_fh_solve(const solve_real *rhs,solve_real *x,int nrhs) {
  int kind=fh_kind,rank=0,mpisize=1;
  PetscInt n=fh_n;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  MPI_Comm_size(PETSC_COMM_WORLD,&mpisize);
  if(kind==FH_NONE&&dbbd_fh_ready())kind=FH_DBBD;
  MPI_Bcast(&kind,1,MPI_INT,0,PETSC_COMM_WORLD);
  if(kind==FH_NONE||nrhs<1)return -1;
  switch(kind) {
  case FH_LU:
    if(rank==0) {
      int slot=1;
      spec48_keep_solve_(&slot,fh_lu_va,fh_lu_irn,&nrhs,(solve_real *)rhs,x);
    }
    break;
  case FH_LU_FR:
    if(rank==0) {
      int ins[4]={(int)n,(int)n,(int)fr_nz,(int)fr_lasize};
      spec48_persist_solve_(ins,fr_irn,fr_values,&nrhs,(solve_real *)rhs,x);
    }
    break;
  case FH_SBBD:
    spec48_nomc66_osolve_(fh_sb_indata,&nrhs,(solve_real *)rhs,x);
    break;
  case FH_SBBD_FR:
    spec48_nomc66_psolve_(fh_sb_indata,&nrhs,(solve_real *)rhs,x);
    break;
  case FH_DBBD:
    dbbd_fh_solve(rhs,x,nrhs,rank,mpisize);
    break;
  }
  if(kind!=FH_DBBD&&mpisize>1)MPI_Bcast(x,(int)(n*nrhs),SORD==1?MPI_DOUBLE:MPI_FLOAT,0,PETSC_COMM_WORLD);
  return 0;
}

void teems_fh_free(void) {
  int kind=fh_kind;
  if(kind==FH_NONE&&dbbd_fh_ready())kind=FH_DBBD;
  MPI_Bcast(&kind,1,MPI_INT,0,PETSC_COMM_WORLD);
  if(kind==FH_LU) {
    int slot=1;
    spec48_keep_free_(&slot);
    free(fh_lu_irn);
    free(fh_lu_va);
    fh_lu_irn=NULL;
    fh_lu_va=NULL;
  }
  if(kind==FH_SBBD)spec48_nomc66_ofree_();
  if(kind==FH_DBBD)dbbd_fh_free();
  fh_kind=FH_NONE;
}

/* ---- one step of iterative refinement (DBBD; -refine, default on) ----
   A DBBD solve is only as accurate as its diagonal blocks, which are
   square subsets of ill-conditioned rectangular region blocks: its
   componentwise backward error (the residual ratio) is 1e-9 to 1e-5 on
   real shocks, against about 1e-15 for LU, and a multi-step run carries
   those residuals into its updates (teems-dev fastrefac_dbbd_fix_report
   section 5). One step x += A^-1 (b - A x), the residual formed in
   double from the kept A and solved with the step's kept factors, brings
   it to LU's level. For nc columns (column-major; x valid on every rank,
   bloc = this rank's rows of b): each rank forms its rows of the
   residual, rank 0 gathers them, one teems_fh_solve takes all columns,
   and every rank adds the correction. *before (when given) receives this
   rank's worst residual ratio of the first column before the step, with
   the floor of the residual check. Collective. */
static void fh_refine_cols(const solve_real *bloc,solve_real *x,int nc,PetscInt rank,PetscInt mpisize,double *before) {
  PetscInt rs,re,n=fh_n,i,k,nc1;
  const PetscInt *rng,*cols;
  const PetscScalar *va;
  int *cnt=NULL,*disp=NULL,j,c;
  solve_real *rl,*rfull=NULL,*d;
  double worst=0.0;
  MatGetOwnershipRange(fh_A,&rs,&re);
  rl=(solve_real *) malloc ((re>rs?(size_t)(re-rs)*nc:1)*sizeof(solve_real));
  for(c=0; c<nc; c++) {
    const solve_real *xc=x+(size_t)c*n,*bc=bloc+(size_t)c*(re-rs);
    double xmax=(c==0&&before!=NULL)?resid_xmax(xc,n):0.0;
    for(i=rs; i<re; i++) {
      double sum=0.0,sabs=fabs((double)bc[i-rs]),amax=0.0,den;
      MatGetRow(fh_A,i,&nc1,&cols,&va);
      for(k=0; k<nc1; k++) {
        double t=(double)va[k]*(double)xc[cols[k]];
        sum+=t;
        sabs+=fabs(t);
        if(fabs((double)va[k])>amax)amax=fabs((double)va[k]);
      }
      MatRestoreRow(fh_A,i,&nc1,&cols,&va);
      rl[(size_t)c*(re-rs)+(i-rs)]=(solve_real)((double)bc[i-rs]-sum);
      if(c==0&&before!=NULL) {
        den=TEEMS_RESID_FLOOR*amax*xmax;
        if(sabs>den)den=sabs;
        if(den>0.0&&fabs((double)bc[i-rs]-sum)/den>worst)worst=fabs((double)bc[i-rs]-sum)/den;
      }
    }
  }
  if(before!=NULL)*before=worst;
  if(rank==0) {
    MatGetOwnershipRanges(fh_A,&rng);
    cnt=(int *) calloc (mpisize,sizeof(int));
    disp=(int *) calloc (mpisize,sizeof(int));
    for(j=0; j<mpisize; j++) {
      cnt[j]=(int)(rng[j+1]-rng[j]);
      disp[j]=(int)rng[j];
    }
    rfull=(solve_real *) malloc ((size_t)n*nc*sizeof(solve_real));
  }
  for(c=0; c<nc; c++)MPI_Gatherv(rl+(size_t)c*(re-rs),(int)(re-rs),SORD==1?MPI_DOUBLE:MPI_FLOAT,rank==0?rfull+(size_t)c*n:NULL,cnt,disp,SORD==1?MPI_DOUBLE:MPI_FLOAT,0,PETSC_COMM_WORLD);
  free(rl);
  free(cnt);
  free(disp);
  d=(solve_real *) malloc ((size_t)n*nc*sizeof(solve_real));
  if(teems_fh_solve(rfull,d,nc)!=0) {
    if(rank==0)errmsg("Error: the refinement step found no kept DBBD factorization for this solve\n");
    MPI_Abort(PETSC_COMM_WORLD,1);
  }
  for(i=0; i<(PetscInt)((size_t)n*nc); i++)x[i]+=d[i];
  free(d);
  free(rfull);
}

void fh_refine_totals(long *solves,double *before,double *after,double *secs) {
  *solves=fh_rf_solves;
  *before=fh_rf_before;
  *after=fh_rf_after;
  *secs=fh_rf_secs;
}

/* -fhtest: solve [b, 2b] again with the kept factorization and compare
   with the step's own solution x and 2x (the kit's check that every
   method can solve more right-hand sides with one factorization, at
   every solve of every driver) */
static void fh_selftest_run(const solve_real *x,PetscInt rank) {
  PetscInt VecSize=fh_n,i;
  solve_real *rhs=NULL,*xx=(solve_real *) malloc (2*VecSize*sizeof(solve_real));
  if(rank==0) {
    rhs=(solve_real *) calloc (2*VecSize,sizeof(solve_real));
    for(i=0; i<VecSize; i++) {
      rhs[i]=fh_bfull[i];
      rhs[VecSize+i]=2.0*fh_bfull[i];
    }
  }
  if(teems_fh_solve(rhs,xx,2)==0&&rank==0) {
    for(i=0; i<2*VecSize; i++) {
      double ref=(i<VecSize)?(double)x[i]:2.0*(double)x[i-VecSize];
      double d=fabs((double)xx[i]-ref),rl=(fabs(ref)>0.0)?d/fabs(ref):d;
      solve_real refr=(solve_real)ref;
      if(memcmp(&xx[i],&refr,sizeof(solve_real))!=0)fh_st_bitdiff++;
      if(d>fh_st_maxabs)fh_st_maxabs=d;
      if(rl>fh_st_maxrel)fh_st_maxrel=rl;
    }
    fh_st_solves++;
    fh_st_rhs+=2;
  }
  free(rhs);
  free(xx);
}

/* the consumers, in order: the self-test (on the step's own solution),
   the refinement step (DBBD), the residual check (of the refined
   solution, which gives the ratio after the step), the subtotal solves
   (refined the same way); then everything kept is released */
void fh_step_end(PetscInt VecSize,solve_real *x,PetscInt rank,PetscInt rank_hsl,PetscInt mpisize) {
  double before=0.0,gbefore=0.0,after=-1.0;
  struct timeval t0= {0,0},t1= {0,0};
  if(teems_fh_selftest)fh_selftest_run(x,rank);
  if(fh_refine_now) {
    gettimeofday(&t0,NULL);
    fh_refine_cols(fh_b,x,1,rank,mpisize,&before);
    MPI_Reduce(&before,&gbefore,1,MPI_DOUBLE,MPI_MAX,0,PETSC_COMM_WORLD);
    gettimeofday(&t1,NULL);
  }
  if(teems_resid_retain&&fh_resid_pending) {
    if(fh_mpi)after=solve_residual_check_mpi(fh_A,fh_b,fh_bstart,fh_bend,x,rank,mpisize);
    else if(rank==rank_hsl)solve_residual_check(fh_A,fh_b,x,VecSize);
  }
  if(fh_refine_now&&rank==0) {
    fh_rf_solves++;
    fh_rf_secs+=(t1.tv_sec-t0.tv_sec)+1e-6*(t1.tv_usec-t0.tv_usec);
    if(gbefore>fh_rf_before)fh_rf_before=gbefore;
    if(after>fh_rf_after)fh_rf_after=after;
    if(after>=0.0)logmsg(1,"Refinement step (DBBD): residual ratio %.3e before, %.3e after\n",gbefore,after);
    else logmsg(1,"Refinement step (DBBD): residual ratio %.3e before\n",gbefore);
  }
  if(teems_sub_active)sub_step_solve(rank,fh_n);
  if(teems_fh_keep)teems_fh_free();
  if(fh_A!=NULL)MatDestroy(&fh_A);
  if(fh_bfull!=fh_b)free(fh_bfull);
  free(fh_b);
  fh_b=NULL;
  fh_bfull=NULL;
  teems_fh_keep=0;
  teems_resid_retain=0;
  fh_resid_pending=0;
  fh_refine_now=0;
}

void fh_selftest_summary(PetscInt rank) {
  if(!teems_fh_selftest||rank!=0)return;
  printf("Factor handle self-test: %ld solves, %ld extra right-hand sides, %ld elements not bit-identical, max abs difference %.3e, max rel difference %.3e\n",fh_st_solves,fh_st_rhs,fh_st_bitdiff,fh_st_maxabs,fh_st_maxrel);
}


static void lu_fastrefac_extract(Mat A,PetscInt VecSize,dim_t laA) {
  Mat_SeqAIJ *aa=(Mat_SeqAIJ*)A->data;
  PetscInt i,j;
  offset_t floorla;
  free(fr_irn);
  free(fr_jcn);
  free(fr_values);
  free(fr_pat);
  fr_nz=aa->nz;
  fr_pat=(PetscInt *) ma48_alloc (fr_nz>0?fr_nz:1,sizeof(PetscInt));
  memcpy(fr_pat,aa->j,fr_nz*sizeof(PetscInt));
  floorla=ma48_la_from_pct(laA,fr_nz);
  /* a starved -laA (<100) must still stage all NE entries */
  if(floorla<fr_nz)floorla=fr_nz;
  if(fr_lasize<floorla)fr_lasize=floorla;
  fr_irn=(int *) ma48_alloc (fr_lasize,sizeof(int));
  fr_jcn=(int *) ma48_alloc (fr_lasize,sizeof(int));
  fr_values=(solve_real *) ma48_alloc (fr_lasize,sizeof(solve_real));
  for(i=0; i<VecSize; i++)for(j=aa->i[i]; j<aa->i[i+1]; j++) {
      fr_irn[j]=i+1;
      fr_jcn[j]=aa->j[j]+1;
      fr_values[j]=aa->a[j];
    }
  fr_ready=0;
}

void lu_fastrefac_solve(Mat A,PetscInt VecSize,dim_t laA,solve_real *rhs,solve_real *x) {
  int insize[6];
  Mat_SeqAIJ *aa=(Mat_SeqAIJ*)A->data;
  PetscInt i;
  int tries;
  if(!fr_ready||aa->nz!=fr_nz||memcmp(fr_pat,aa->j,fr_nz*sizeof(PetscInt))!=0) {
    if(fr_ready)printf("Note: the Jacobian pattern changed since the last analyse (an entry crossed zero); re-analysing the LU pivot sequence\n");
    lu_fastrefac_extract(A,VecSize,laA);
  }
  else {
    for(i=0; i<fr_nz; i++)fr_values[i]=aa->a[i];
  }
  probe_onfail_scope_set(A,VecSize,VecSize,"condensed system",-1,NULL,NULL,0,0,0,0);
  teems_condest_scope=teems_condest;
  for(tries=0; tries<6; tries++) {
    insize[0]=VecSize;
    insize[1]=VecSize;
    insize[2]=fr_nz;
    insize[3]=(int)fr_lasize;
    insize[4]=fr_ready;
    insize[5]=0;
    spec48_ssol2la_p_(insize,fr_irn,fr_jcn,fr_values,rhs,x);
    if(insize[4]==0) {
      fr_ready=1;
      teems_condest_scope=0;
      probe_onfail_scope_clear();
      solve_residual_check(A,rhs,x,VecSize);
      if(teems_fh_keep)fh_kind=FH_LU_FR;
      return;
    }
    if(insize[4]==-3) {
      /* MA48 workspace too small: grow to at least its suggested size
         (doubling floor guarantees progress) and redo the analyse */
      offset_t newla=ma48_grow_la(fr_lasize,insize[5],fr_nz,"laA",&teems_laA_used);
      fr_lasize=newla;
    }
    /* -3 or fast-factorize declined: fresh analyse on current values */
    lu_fastrefac_extract(A,VecSize,laA);
  }
  errmsg("Error: the MA48 workspace for %s did not converge after %d growth attempts; raise the initial workspace (laA/laD/laDi) or use a bordered matrix_method (\"SBBD\" or \"DBBD\")\n",probe_onfail_scope_label(),tries);
  MPI_Abort(PETSC_COMM_WORLD,1);
}

void lu_fastrefac_free(void) {
  free(fr_irn);
  free(fr_jcn);
  free(fr_values);
  free(fr_pat);
  fr_pat=NULL;
  fr_irn=NULL;
  fr_jcn=NULL;
  fr_values=NULL;
  fr_nz=-1;
  fr_lasize=0;
  fr_ready=0;
  spec48_persist_free_();
}

/* One-shot sequential LU with MA48 workspace growth: stages the
   nonzero COO from A -- which must stay live across the call, since
   MA48 clobbers the staged arrays in place -- and on a -3 workspace
   return grows LA to max(MA48's suggested size, 2x) and re-stages.
   rank_hsl only. */
static offset_t grow_hw=0;     /* largest LA a previous step needed */
static PetscInt grow_hw_n=-1;  /* the system size it applies to (a different system restarts
                                  the mark).  NOT the realized NE: the assembly skips
                                  zero-valued entries, so NE drifts by a few entries between
                                  steps as values cross zero, and keying on it restarted the
                                  mark at every Gragg step (I-long: 300->600 re-grown 17x) */

void lu_grow_solve(Mat A,PetscInt VecSize,dim_t laA,solve_real *rhs,solve_real *x) {
  Mat_SeqAIJ *aa=(Mat_SeqAIJ*)A->data;
  PetscInt i,j;
  offset_t count=0,k,lasize;
  int tries;
  int insize[8];
  for(i=0; i<aa->nz; i++) if(aa->a[i]!=0)count++;
  lasize=ma48_la_from_pct(laA,count);
  /* a starved -laA (<100) must still stage all NE entries; MA48 then
     returns -3 with its suggested size and the growth loop takes over */
  if(lasize<count)lasize=count;
  /* what a previous step had to grow to: the fill pattern is stable
     across steps, so starting from the high-water mark spares every
     later step the failed factorization that discovered it (oversizing
     is untouched pages, undersizing costs a whole factorization) */
  if(lasize<grow_hw&&VecSize==grow_hw_n)lasize=grow_hw;
  for(tries=0; tries<6; tries++) {
    int *irn=(int *) ma48_alloc (lasize,sizeof(int));
    int *jcn=(int *) ma48_alloc (lasize,sizeof(int));
    solve_real *values=(solve_real *) ma48_alloc (lasize,sizeof(solve_real));
    k=0;
    for(i=0; i<VecSize; i++)for(j=aa->i[i]; j<aa->i[i+1]; j++) if(aa->a[j]!=0) {
          irn[k]=i+1;
          jcn[k]=aa->j[j]+1;
          values[k]=aa->a[j];
          k++;
        }
    insize[0]=VecSize;
    insize[1]=VecSize;
    insize[2]=count;
    insize[3]=laA;
    insize[4]=0;
    insize[5]=(int)lasize;
    insize[6]=teems_fh_keep?1:0; /* Tier B3: keep the factorization for more solves */
    insize[7]=0;
    teems_condest_scope=teems_condest;
    spec48_ssol2la_(insize,irn,jcn,values,rhs,x);
    teems_condest_scope=0;
    free(jcn);
    if(teems_fh_keep&&insize[4]==0) {
      fh_lu_irn=irn;
      fh_lu_va=values;
      fh_kind=FH_LU;
    }
    else {
      free(irn);
      free(values);
    }
    if(insize[4]<0&&insize[4]!=-3) {
      /* soft failure (singular stage state under adaptive Runge-Kutta):
         hand the step back to the driver for a retry */
      teems_stage_solve_failed=1;
      return;
    }
    if(insize[4]!=-3) {
      /* record what this step actually needed so the next one starts
         there instead of rediscovering it with a failed factorization */
      if(VecSize!=grow_hw_n||lasize>grow_hw) {
        grow_hw=lasize;
        grow_hw_n=VecSize;
      }
      solve_residual_check(A,rhs,x,VecSize);
      return;
    }
    lasize=ma48_grow_la(lasize,insize[5],count,"laA",&teems_laA_used);
  }
  errmsg("Error: the MA48 workspace for %s did not converge after %d growth attempts; raise the initial workspace (laA/laD/laDi) or use a bordered matrix_method (\"SBBD\" or \"DBBD\")\n",probe_onfail_scope_label(),tries);
  MPI_Abort(PETSC_COMM_WORLD,1);
}

/* SBBD persistent MP48 instance (-fastrefac): border lists, per-block
   pivot sequences and factors live across steps; repeat steps refill
   VALUES/B and refactorize with FACT_JOB=2.  Same full-pattern
   extraction rule as the LU path.  Collective: every rank calls
   sbbd_fastrefac_solve with the redo decision broadcast from rank_hsl. */
static int sb_ready=0; /* a persistent instance exists (for the pattern-change note) */

/* SBBD one-shot solve (6.15(c)): MP48's host arrays are staged straight
   from PETSc's SeqAIJ CSR on the solving rank (stored zeros dropped, as
   before), the PETSc matrix and right-hand side are destroyed, and only
   then does MP48 factorize -- so rank_hsl holds one copy of the system
   during the factorization instead of the matrix, a 20-byte-per-entry
   COO staging and the host arrays together.  The on-failure probe scope
   points at the host arrays (registered from the Fortran staging).
   Collective: every rank enters both Fortran calls. */
void sbbd_csr_solve(Mat *A,Vec *vecb,PetscInt VecSize,PetscInt rank,PetscInt rank_hsl,
                    fortran_int *indata,MPI_Fint fcomm,
                    offset_t *counteq,offset_t *countvarintra1,solve_real *x) {
  Mat_SeqAIJ *aa=(Mat_SeqAIJ*)(*A)->data;
  PetscScalar *bv=NULL;
  PetscErrorCode ierr;
  indata[1]=VecSize;
  if(rank==rank_hsl)VecGetArray(*vecb,&bv);
  spec48_nomc66_stage_(indata,aa->i,aa->j,aa->a,bv,&fcomm,counteq,countvarintra1);
  if(rank==rank_hsl)VecRestoreArray(*vecb,&bv);
  ierr = VecDestroy(vecb);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = MatDestroy(A);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
  {
    int keep=teems_fh_keep?1:0;
    spec48_nomc66_run_(indata,x,&fcomm,&keep);
    if(keep) {
      memcpy(fh_sb_indata,indata,sizeof(fh_sb_indata));
      fh_kind=FH_SBBD;
    }
  }
  probe_onfail_scope_clear();
  solve_x_check(x,VecSize,rank==rank_hsl);
  residual_note_skipped(rank,rank_hsl);
}

/* -fastrefac SBBD (6.15(c) form): the persistent MP48 instance is staged
   and refilled straight from PETSc's CSR; no C-side copy of the pattern
   or values is kept.  The pattern test runs on the host against the
   instance's own EQVAR and its verdict is broadcast, so every rank calls
   the kernel with the same redo.  A and vecb stay live until after the
   call because a declined fast step is retried as a full rebuild from
   the same CSR. */
void sbbd_fastrefac_solve(Mat *A,Vec *vecb,PetscInt VecSize,PetscInt rank,PetscInt rank_hsl,
                          fortran_int *indata,MPI_Fint fcomm,
                          offset_t *counteq,offset_t *countvarintra1,solve_real *x) {
  int same=0,redo_io;
  PetscErrorCode ierr;
  PetscScalar *bv=NULL;
  Mat_SeqAIJ *aa=(Mat_SeqAIJ*)(*A)->data;
  indata[1]=VecSize;
  if(rank==rank_hsl) {
    indata[0]=aa->nz;
    spec48_nomc66_p_same_(indata,aa->i,aa->j,&same);
    if(!same&&sb_ready)printf("Note: the Jacobian pattern changed since the last analyse (an entry crossed zero); rebuilding the MP48 instance\n");
  }
  MPI_Bcast(&same,1,MPI_INT,rank_hsl,PETSC_COMM_WORLD);
  if(rank==rank_hsl)VecGetArray(*vecb,&bv);
  redo_io=same?1:0;
  spec48_nomc66_p_csr_(indata,aa->i,aa->j,aa->a,bv,x,&fcomm,counteq,countvarintra1,&redo_io);
  if(redo_io<0) {
    /* fast refactorize declined: rebuild the instance on current values */
    redo_io=0;
    spec48_nomc66_p_csr_(indata,aa->i,aa->j,aa->a,bv,x,&fcomm,counteq,countvarintra1,&redo_io);
    if(redo_io<0) {
      printf("MP48 instance rebuild failed with code %d\n",redo_io);
      /* -21 = structurally rank-deficient interface matrix (singular
         system; +2 warnings are mapped here by the kernel): name the
         defects before dying.  Other codes are setup/allocation/IO —
         a matching would be noise (the MC64/-3 lesson). */
      if(rank==rank_hsl&&redo_io==-21) {
        int zdiag=0;
        teems_onfail_diag_(&zdiag);
      }
      /* redo_io is replicated (the kernel sets it from pdata%ERROR on
         every rank), so all ranks are here: hold the siblings until
         the diagnosing rank has printed, or its output is truncated
         by their abort */
      MPI_Barrier(PETSC_COMM_WORLD);
      MPI_Abort(PETSC_COMM_WORLD,1);
    }
  }
  probe_onfail_scope_clear();
  if(rank==rank_hsl)solve_residual_check(*A,bv,x,VecSize);
  if(rank==rank_hsl)VecRestoreArray(*vecb,&bv);
  ierr = VecDestroy(vecb);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = MatDestroy(A);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
  sb_ready=1;
  if(teems_fh_keep) {
    memcpy(fh_sb_indata,indata,sizeof(fh_sb_indata));
    fh_kind=FH_SBBD_FR;
  }
}

void sbbd_fastrefac_free(void) {
  sb_ready=0;
  spec48_nomc66_pfree_();
}

/* Exogenous-side layout, shared by every site that (re)builds the shock
   vector or B (see the header comment).  A square B keeps the row split
   -- including NDBBD's time blocks, which the interface machinery needs
   -- while a wide B (nexo > VecSize under heavy condensation) takes
   PETSc's split of BSize on both, so MatMult(B,vece,vecb) stays
   layout-compatible. */
void shock_vec_set_sizes(Vec v,int nesteddbbd,PetscInt localsize,PetscInt VecSize,PetscInt BSize) {
  if(nesteddbbd==1&&BSize==VecSize)VecSetSizes(v,localsize,VecSize);
  else VecSetSizes(v,PETSC_DECIDE,BSize);
}

void shock_mat_set_sizes(Mat B,int nesteddbbd,PetscInt localsize,PetscInt VecSize,PetscInt BSize) {
  if(nesteddbbd==1)MatSetSizes(B,localsize,(BSize==VecSize)?localsize:PETSC_DECIDE,VecSize,BSize);
  else MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,VecSize,BSize);
}

bool solve_johansen(PetscBool nohsl,PetscInt VecSize,Mat A,PetscInt dnz,PetscInt* dnnz,PetscInt onz,PetscInt* onnz,Mat B,PetscInt dnzB,PetscInt* dnnzB,PetscInt onzB,PetscInt* onnzB,Vec vecb,Vec vece,PetscInt rank,PetscInt rank_hsl,PetscInt mpisize,char* tabfile, char *commsyntax,set_def *sets,dim_t nset, set_element *set_elems, array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar, elem_value **elem_vals2,offset_t ncofvar,offset_t ncofele,offset_t nvarele,closure_entry **closure_vals2,offset_t alltimeset,offset_t allregset,offset_t nintraeq,dim_t matsol,PetscInt Istart,PetscInt Iend,  offset_t nreg, offset_t ntime, PetscInt *eq_addr, offset_t ndblock, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,dim_t laDi,dim_t laD,PetscReal ma48_cntl4,dim_t nesteddbbd,int localsize,PetscInt *ndbbddrank1,fortran_int* indata,dim_t mc66,fortran_int *ptx,struct timeval begintime,solve_real **xcf2){ //Johansen
  char tempfilenam[256],tempchar[256];
  PetscScalar value,*vals=NULL;
  PetscErrorCode ierr;
  PetscInt count=0,nz01=0,*ai=NULL,*aj=NULL; /* stay 0 on ranks != rank_hsl */
  fortran_int k=0,m=1;
  offset_t i,j;
  solve_real *b1=NULL,*x0=NULL;
  PetscBool presol;/* NDBBD phase flag: presolve pass writes the interface files, solve pass consumes them */
  bool IsIni;
  FILE* tempvar;
  closure_entry *closure_vals;
  closure_vals=*closure_vals2;
  solve_real *xcf;
  xcf=*xcf2;

  clock_t timeend;
  struct timeval endtime;
  struct timespec gettime_beg,gettime_end;
  double rep_time;
  elem_value *elem_vals;
  elem_vals=*elem_vals2;
  elem_value *elem_vals1=NULL;
  /* backsolve recovery workspace: exogenous per-step changes as placed in
     vece, and the recovered changes of the backsolved elements */
  solve_real *exo_z=NULL,*bsvals=NULL;
  if(nbselems>0) {
    exo_z= (solve_real *) calloc (nvarele,sizeof(solve_real));
    bsvals= (solve_real *) calloc (nbselems,sizeof(solve_real));
  }
  /* exogenous columns run 0..nexo-1; under heavy condensation nexo can
     exceed VecSize, so vece and B's columns span BSize (see
     shock_vec_set_sizes / shock_mat_set_sizes for the layout rule) */
  PetscInt BSize;
  BSize=(PetscInt)(nvarele-VecSize-nbselems);      /* nexo */
  BSize=(BSize>VecSize)?BSize:VecSize;

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
    jac_mat_prealloc(A,"Jacobian",nohsl,rank==rank_hsl,Iend-Istart,dnz,dnnz,onz,onnz);

    ierr = MatSetOption(A,MAT_SYMMETRIC,PETSC_FALSE);
    CHKERRQ(ierr);


    if(nohsl) {
      MatCreate(PETSC_COMM_WORLD,&B);
    }
    else {
      MatCreate(PETSC_COMM_SELF,&B);
    }
    shock_mat_set_sizes(B,nesteddbbd,localsize,VecSize,BSize);
    if(nohsl) {
      MatSetType(B,MATMPIAIJ);
    }
    else {
      MatSetType(B,MATSEQAIJ);
    }
    jac_mat_prealloc(B,"exogenous block",nohsl,rank==rank_hsl,Iend-Istart,dnzB,dnnzB,onzB,onnzB);

    gettimeofday(&endtime, NULL);
    if(rank==0)logmsg(1,"Matrix preparation time %.2f s\n",(endtime.tv_sec - begintime.tv_sec)+((double)(endtime.tv_usec - begintime.tv_usec))/ 1000000);
    teems_rss_probe("matrix preparation");
    
    if(rank==rank_hsl) {
      jacobian_fill(tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,closure_vals,ndblock,alltimeset,allregset,eq_addr,counteq,nintraeq,A,B);
    }

    gettimeofday(&begintime, NULL);
    if(rank==0)logmsg(1,"Matrix calculation time %.2f s\n",(begintime.tv_sec - endtime.tv_sec)+((double)(begintime.tv_usec - endtime.tv_usec))/ 1000000);
    teems_rss_probe("matrix calculation");

    for (count=0; count<nvarele; count++) {
      if (CL_EXO(count)) {
        value = CL_SHOCK(count);
        dnz=closure_vals[count].exo_index;
        VecSetValues(vece,1,&dnz,&value,INSERT_VALUES);
        if(exo_z!=NULL)exo_z[count]=CL_SHOCK(count);
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
      if ( (tempvar = teems_fopen(tempfilenam, "wb")) == NULL ) {
        errmsg("Error: cannot open %s for writing\n",tempfilenam);
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
      if ( (tempvar = teems_fopen(tempfilenam, "wb")) == NULL ) {
        errmsg("Error: cannot open %s for writing\n",tempfilenam);
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
    teems_rss_probe("matrix assembly");
    CHKERRQ(ierr);
    PetscViewer viewer;
    /* vecb spans the equation rows (VecSize); vece may be wider (BSize) */
    {
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
      /* the (N)DBBD back-solves probe vecb with -1 sentinels for rows
         other ranks own; VecDuplicate(vece) used to inherit this option
         before vecb was created explicitly (the un-ignored -1 made
         VecGetValues error out and drop the rest of the batch) */
      VecSetOption(vecb, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
    }
    ierr = MatMult(B,vece,vecb);
    CHKERRQ(ierr);
    ierr = VecAssemblyBegin(vecb);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(vecb);
    CHKERRQ(ierr);
    sub_step_rhs(B,vece,VecSize,matsol,rank,rank_hsl);
    ierr = MatDestroy(&B);
    CHKERRQ(ierr);
    ierr = VecDestroy(&vece);
    CHKERRQ(ierr);
    fh_step_begin(A,vecb,matsol,mc66,VecSize,rank,rank_hsl);

    if(matsol>=MM_DBBD) {
      gettimeofday(&begintime, NULL);
      clock_gettime(CLOCK_REALTIME, &gettime_beg);
      int *row_order= (int *) calloc (VecSize,sizeof(int));
      int *col_order= (int *) calloc (VecSize,sizeof(int));
      int *block_sizes= (int *) calloc (ndblock,sizeof(int));
      if(matsol==MM_DBBD) {
        dbbd_order(A,VecSize,mpisize,rank,Istart,Iend,nvarele,eq_addr,row_order,col_order,ndblock,block_sizes,countvarintra1,counteq,counteqnoadd,laA,ma48_cntl4);
        x0=realloc (x0,VecSize*sizeof(solve_real));
        dbbd_solve(A,vecb,x0,VecSize,mpisize,rank,Istart,Iend,row_order,col_order,ndblock,block_sizes,countvarintra1,counteq,counteqnoadd,laA,laD);//,iter
      }
      if(matsol==MM_NDBBD) {
        presol=1;
        if(presol){
        teems_stage_mark("order-presolve");
        ndbbd_order_presolve(A,VecSize,mpisize,rank,Istart,Iend,nreg,ntime,nvarele,eq_addr,row_order,col_order,ndblock,block_sizes,countvarintra1,counteq,counteqnoadd,laA,laDi,ma48_cntl4,ndbbddrank1,presol);
        teems_stage_mark("presolve");
        ndbbd_presolve(A,vecb,x0,VecSize,mpisize,rank,Istart,Iend,row_order,col_order,ndblock,nreg,ntime,block_sizes,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,ma48_cntl4,presol);//,iter
        }
        presol=0;
        teems_stage_mark("order");
        ndbbd_order(A,VecSize,mpisize,rank,Istart,Iend,nreg,ntime,nvarele,eq_addr,row_order,col_order,ndblock,block_sizes,countvarintra1,counteq,counteqnoadd,laA,laDi,ma48_cntl4,ndbbddrank1,presol);
        x0=realloc (x0,VecSize*sizeof(solve_real));
        teems_stage_mark("solve");
        ndbbd_solve(A,vecb,x0,VecSize,mpisize,rank,Istart,Iend,row_order,col_order,ndblock,nreg,ntime,block_sizes,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,ma48_cntl4,presol);//,iter
        teems_stage_report("ndbbd");
      }
      time(&timeend);
      gettimeofday(&endtime, NULL);
      clock_gettime(CLOCK_REALTIME, &gettime_end);
      rep_time = ((double)(gettime_end.tv_nsec-gettime_beg.tv_nsec))/1000000000.0;
      if(rank==0)logmsg(1,"Step time %.2f s\n",(endtime.tv_sec - begintime.tv_sec)+((double)(endtime.tv_usec - begintime.tv_usec))/ 1000000);
      if(rank==0)logmsg(1,"Step wall time %.2f s\n",rep_time);
      teems_rss_probe("step");
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
        if(mc66==0) {
          x0=realloc (x0,VecSize*sizeof(solve_real));
          sbbd_csr_solve(&A,&vecb,VecSize,rank,rank_hsl,ptx,fcomm,counteq,countvarintra1,x0);
        }
        else {
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
        if(rank==rank_hsl)probe_onfail_scope_set_coo(irn,jcn,values,NULL,count,VecSize,VecSize,"SBBD system (MP48)",NULL,NULL);
        if(mc66!=0)spec48_single_(ptx,irn,jcn,b1,values,x0,neleperrow,ai1,&fcomm);
        if(mc66==0)spec48_nomc66_(ptx,jcn,b1,values,x0,neleperrow,&fcomm,counteq,countvarintra1);
              solve_x_check(x0,VecSize,rank==rank_hsl);
              residual_note_skipped(rank,rank_hsl);
        probe_onfail_scope_clear();
        free(irn);
        free(jcn);
        free(values);
        free(neleperrow);
        free(ai1);
        free(b1);//b1=realloc (b1,sizeof(ha_cgetype));
        b1=NULL;
        }
      }
      else {
        x0=realloc (x0,VecSize*sizeof(solve_real));
        if(rank==rank_hsl)VecGetArray(vecb,&vals);
        PetscViewerDrawOpen(PETSC_COMM_WORLD,0,"",0,0,500,500,&viewer);
        ierr = MatView(A,viewer);
        PetscViewerDestroy(&viewer);

        /* A stays live through the factorize for the on-failure
           diagnosis and the workspace-growth re-staging (the staged
           COO is MA48 workspace and is clobbered) */
        if(rank==rank_hsl) {
          probe_onfail_scope_set(A,VecSize,VecSize,"condensed system",-1,NULL,NULL,0,0,0,0);
          lu_grow_solve(A,VecSize,laA,vals,x0);
          probe_onfail_scope_clear();
        }
        ierr = MatDestroy(&A);
        CHKERRQ(ierr);
        ierr = VecDestroy(&vecb);
        CHKERRQ(ierr);
      }
      gettimeofday(&endtime, NULL);
      if(rank==0)logmsg(1,"Step time %.2f s\n",(endtime.tv_sec - begintime.tv_sec)+((double)(endtime.tv_usec - begintime.tv_usec))/ 1000000);
      teems_rss_probe("step");
    }
    fh_step_end(VecSize,x0,rank,rank_hsl,mpisize);
    if(rank==rank_hsl) {
      if(!inmemory){
      if ((tempvar = teems_fopen(tempfilenam, "rb")) == NULL) {
        errmsg("Error: cannot open %s for reading\n",tempfilenam);
      }
      *elem_vals2=(elem_value*)realloc (*elem_vals2,(ncofele+nvarele)*sizeof(elem_value));
      scratch_read(*elem_vals2, sizeof(elem_value),ncofele+nvarele,tempvar,tempfilenam);
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
      if ((tempvar = teems_fopen(tempfilenam, "rb")) == NULL) {
        errmsg("Error: cannot open %s for reading\n",tempfilenam);
      }
      *closure_vals2=(closure_entry*)realloc (*closure_vals2,(nvarele)*sizeof(closure_entry));
      scratch_read(*closure_vals2, sizeof(closure_entry),nvarele,tempvar,tempfilenam);
      fclose(tempvar);
      remove(tempfilenam);
      }
      closure_vals=*closure_vals2;      
    }
    *xcf2=(solve_real*)realloc (*xcf2,nvarele*sizeof(solve_real));
    xcf=*xcf2;
    if(rank==rank_hsl) {
      elem_vals1=elem_vals+ncofele;
      /* recover the backsolved elements from their defining equations
         before the updates read any variable's change (GEMPACK 14.1.3) */
      if(nbselems>0)backsolve_recover(tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele,closure_vals,x0,exo_z,bsvals);
      for(i=0; i<nvar; i++) {
        if(vars[i].change_real) {
          for(j=vars[i].offset; j<vars[i].nelem+vars[i].offset; j++) {
            if(CL_EXO(j)) {
              elem_vals1[j].initial=elem_vals1[j].value;
              elem_vals1[j].value+=CL_SHOCK(j);
              xcf[j]=CL_SHOCK(j);//varchange[j]
              elem_vals1[j].substep_base=CL_SHOCK(j);
            }
            else if(CL_BS(j)) {
              elem_vals1[j].initial=elem_vals1[j].value;
              elem_vals1[j].value+=bsvals[closure_vals[j].exo_index];
              xcf[j]=bsvals[closure_vals[j].exo_index];
              elem_vals1[j].substep_base=bsvals[closure_vals[j].exo_index];
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
            if(CL_EXO(j)) {
              elem_vals1[j].initial=elem_vals1[j].value;
              elem_vals1[j].value+=CL_SHOCK(j)*elem_vals1[j].initial/100;
              xcf[j]=CL_SHOCK(j);//varchange[j]
              elem_vals1[j].substep_base=CL_SHOCK(j);
            }
            else if(CL_BS(j)) {
              elem_vals1[j].initial=elem_vals1[j].value;
              xcf[j]=bsvals[closure_vals[j].exo_index];
              elem_vals1[j].value+=bsvals[closure_vals[j].exo_index]/100*elem_vals1[j].value;
              elem_vals1[j].substep_base=bsvals[closure_vals[j].exo_index];
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
      if(teems_sub_active&&rank==0)sub_update(SUB_JOHANSEN,tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele,closure_vals,nvarele,exo_z,NULL);
      updates_apply(tabfile,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,0);
      strcpy(commsyntax,"formula");
      IsIni=false;
      statements_execute(tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,IsIni,0);

    }
    elem_vals1=NULL;
    free(x0);
    free(exo_z);
    free(bsvals);
    return 1;
  }

/* -convrule 1: GEMPACK's treatment of poorly converging components
   (manual 26.2.5, 26.2.5.1; default off, as it changes reported
   results). Per element and subinterval the three pass results c1 c2
   c3 (step counts n1<n2<n3) are classified and the reported increment
   is
     - their average when they lie very close together or near zero
       (GEMPACK's MVC/OVC/MN0/ON0);
     - the extrapolation when the two two-pass extrapolations, from
       passes 1,2 and from passes 2,3, agree to at least one figure
       (manual 26.2.1's figures column; the confident codes CX/FCX);
     - the third result otherwise (MC?/MD?/MD!/OC?/OD?/OD!; the
       manual's GC10E1 example, figures 0, reports the 8-step result).
   The manual does not publish the code thresholds; "very close" is
   agreement to about six figures. The same rule applies to the
   (change)/(explicit) Update targets' extrapolated path values.
   conv_c holds the three passes, conv_w their extrapolation weights. */
static solve_real *conv_c[3]={NULL,NULL,NULL};
static double conv_w[3];
static offset_t conv_n[3];

int conv_pick(double c1,double c2,double c3,double q2,double q3,double E,double *R) {
  double sc=fmax(fabs(c1),fmax(fabs(c2),fabs(c3)));
  double e12=(q2*c2-c1)/(q2-1.0),e23=(q3*c3-q2*c2)/(q3-q2);
  if(sc<=1e-6||fabs(c3-c1)<=1e-6*sc) { *R=(c1+c2+c3)/3; return 1; }
  if(fabs(e12-e23)<0.1*fmax(fabs(e12),fabs(e23))) { *R=E; return 0; }
  *R=c3;
  return 2;
}

void conv_ratios(bool euler,double *q2,double *q3) {
  double p=euler?1.0:2.0;
  *q2=pow((double)llround(steps1*step_ratio2)/steps1,p);
  *q3=pow((double)llround(steps1*step_ratio3)/steps1,p);
}

static void conv_store(int sol, offset_t nvarele, const solve_real *varchange, double w) {
  conv_c[sol]=(solve_real *) realloc (conv_c[sol],(nvarele>0?nvarele:1)*sizeof(solve_real));
  memcpy(conv_c[sol],varchange,nvarele*sizeof(solve_real));
  conv_w[sol]=w;
}

static void conv_apply(array_def *vars, offset_t nvar, dim_t subindx, const solve_real *xc0, solve_real *xcf, bool euler) {
  offset_t i,k;
  double q2,q3;
  conv_ratios(euler,&q2,&q3);
  for(i=0; i<nvar; i++) for(k=vars[i].offset; k<vars[i].offset+vars[i].nelem; k++) {
      double c1=conv_c[0][k],c2=conv_c[1][k],c3=conv_c[2][k];
      double E=conv_w[0]*c1+conv_w[1]*c2+conv_w[2]*c3,R;
      int o=conv_pick(c1,c2,c3,q2,q3,E,&R);
      conv_n[o]++;
      if(o!=0) xcf[k]+=(R-E)*((!vars[i].change_real&&subindx>0)?xc0[k]:1);
    }
}

/* Extrapolation accuracy side-car (<stem>.xac, Tier B2; GEMPACK's XAC
   file, manual 26.2.3): per variable element the three pass solutions
   the Richardson extrapolation combines and the figures-of-accuracy
   code that the log only summarises. int64 header {version 1, nrow =
   nvarele, npass 3, nsubints}, then npass x nrow doubles (pass by pass,
   rows in .bin order), then nrow int32 codes: 6 = the passes agree to
   6 or more figures, 5..2 = that many, 1 = one figure or none (the
   minimum over subintervals). With subintervals a pass value is that
   pass over the last subinterval, compounded onto the extrapolated
   result of the ones before (the path the last subinterval starts
   from). Rank 0 only; a run that stops midway leaves a partial file
   that <stem>.outputs.json does not list. */
static FILE *xac_fp=NULL;
static char xac_path[TABREADLINE+8];
static solve_real *xac_base=NULL;

static void xac_pass(int sol,dim_t subindx,dim_t subints,array_def *vars,offset_t nvar,offset_t nvarele,const solve_real *varchange,const solve_real *xcf) {
  enum { CHUNK=1<<14 };
  solve_real buf[CHUNK];
  offset_t i,k,n,m;
  if(sol==0) {
    offset_t hdr[4]={1,nvarele,3,(offset_t)(teems_comp_nsub>1?teems_comp_nsub:subints)};
    if(xac_fp!=NULL)fclose(xac_fp);
    snprintf(xac_path,sizeof(xac_path),"%s.xac",teems_sol_stem);
    xac_fp=fopen(xac_path,"wb");
    if(xac_fp==NULL) {
      errmsg("Error: cannot open %s for writing: %s (the solver runs as uid %d)\n",xac_path,strerror(errno),(int)getuid());
      return;
    }
    fwrite(hdr,sizeof(offset_t),4,xac_fp);
    free(xac_base);
    xac_base=NULL;
    if(subindx>0) {
      xac_base=(solve_real *) malloc (nvarele*sizeof(solve_real));
      if(xac_base==NULL) {
        printf("Warning: no memory for the extrapolation accuracy file; %s not written\n",xac_path);
        fclose(xac_fp);
        xac_fp=NULL;
        remove(xac_path);
        return;
      }
      memcpy(xac_base,xcf,nvarele*sizeof(solve_real));
    }
    else if(teems_comp_xac_base!=NULL) {
      xac_base=(solve_real *) malloc (nvarele*sizeof(solve_real));
      if(xac_base!=NULL)memcpy(xac_base,teems_comp_xac_base,nvarele*sizeof(solve_real));
    }
  }
  if(xac_fp==NULL)return;
  for(i=0; i<nvar; i++) {
    if(vars[i].nelem<=0)continue;
    fseeko(xac_fp,(off_t)(4*sizeof(offset_t))+(off_t)sizeof(solve_real)*((off_t)sol*nvarele+vars[i].offset),SEEK_SET);
    for(k=0; k<vars[i].nelem; k+=n) {
      n=vars[i].nelem-k; if(n>CHUNK)n=CHUNK;
      for(m=0; m<n; m++) {
        offset_t e=vars[i].offset+k+m;
        if(xac_base==NULL)buf[m]=varchange[e];
        else if(vars[i].change_real)buf[m]=xac_base[e]+varchange[e];
        else buf[m]=xac_base[e]+varchange[e]*(1+xac_base[e]/100);
      }
      fwrite(buf,sizeof(solve_real),n,xac_fp);
    }
  }
}

static void xac_codes(offset_t nvarele,const int *codes) {
  if(xac_fp==NULL)return;
  fseeko(xac_fp,(off_t)(4*sizeof(offset_t))+(off_t)sizeof(solve_real)*3*nvarele,SEEK_SET);
  fwrite(codes,sizeof(int),nvarele,xac_fp);
  if(fclose(xac_fp)!=0)errmsg("Error: cannot write %s: %s\n",xac_path,strerror(errno));
  else outputs_note(xac_path,"extrapolation_accuracy",1);
  xac_fp=NULL;
  free(xac_base);
  xac_base=NULL;
}

bool solve_gragg(PetscBool nohsl,PetscInt VecSize,Mat* A1,PetscInt dnz,PetscInt* dnnz,PetscInt onz,PetscInt* onnz,Mat* B1,PetscInt dnzB,PetscInt* dnnzB,PetscInt onzB,PetscInt* onnzB,Vec* vecb1,Vec *vece1,PetscInt rank,PetscInt rank_hsl,PetscInt mpisize,char* tabfile, char *commsyntax,set_def *sets,dim_t nset, set_element *set_elems, array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar, elem_value **elem_vals2,offset_t ncofvar,offset_t ncofele,offset_t nvarele,closure_entry **closure_vals2,offset_t alltimeset,offset_t allregset,offset_t nintraeq,dim_t matsol,PetscInt Istart,PetscInt Iend,  offset_t nreg, offset_t ntime, PetscInt *eq_addr, offset_t ndblock, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,dim_t laDi,dim_t laD,PetscReal ma48_cntl4,dim_t nesteddbbd,int localsize,PetscInt *ndbbddrank1,fortran_int* indata,dim_t mc66,fortran_int *ptx,struct timeval begintime,dim_t subints,MPI_Fint fcomm,int solmethod,solve_real **xcf2){ /* multistep driver: Gragg (smoothed modified midpoint, Pearson 1991 eq. 6.1 / Alg. 7.1.2) or forward Euler, per solmethod */
  char tempfilenam[256],tempchar[256],solchar[255];
  PetscScalar *vals;
  PetscErrorCode ierr;
  PetscInt count=0,nz01=0,*ai=NULL,*aj=NULL; /* stay 0 on ranks != rank_hsl */
  fortran_int k=0,m=1;
  fortran_int tindx1;//,tindx2;
  solve_real temp1,temp2;
  offset_t i,j;
  dim_t subindx;
  solve_real *b1=NULL;
  solve_real *x1=NULL;
  solve_real *xcf;
  xcf=*xcf2;
  PetscBool presol;/* NDBBD phase flag: presolve pass writes the interface files, solve pass consumes them */
  bool IsIni;
  FILE* tempvar;
  PetscLogDouble time1,time0;
  clock_t timestr,timeend;
  struct timeval endtime;
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
  solve_real vpercents=1.0;
  FILE* solution;
  /* -single_run: one pass, no extrapolation triple */
  int maxsol=teems_single_run?1:(teems_two_run?2:3);
  /* Euler: forward step on every substep (no leapfrog, no terminal
     smoothing pass) and an h — not h^2 — truncation error series, so
     the Richardson weights below use the step ratios unsquared */
  bool euler=(solmethod==SM_EULER);
  /* midpoint: Gragg's leapfrog without the terminal smoothing pass --
     N passes for N steps, the result is the last leapfrog point (manual
     30.2); it never takes the exogenous variables past their end point */
  bool midpoint=(solmethod==SM_MIDPOINT);
  /* -fastrefac: sequential LU keeps the MA48 pivot sequence across
     steps and refactorizes with JOB=2 (analyse runs once per solve) */
  dim_t fastrefac=0;
  PetscOptionsGetInt(NULL,NULL,"-fastrefac",&fastrefac,NULL);
  dim_t convrule=0;
  PetscOptionsGetInt(NULL,NULL,"-convrule",&convrule,NULL);
  if(convrule&&(teems_single_run||teems_two_run||teems_sub_active)) {
    if(rank==0)errmsg("Error: -convrule 1 applies to a run extrapolating from three solutions without subtotals (manual 26.2.5)\n");
    MPI_Abort(PETSC_COMM_WORLD,1);
  }
  conv_n[0]=conv_n[1]=conv_n[2]=0;
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
    /* backsolve recovery workspace: exogenous per-step changes as placed
       in vece (captured at every vece fill site; used by the recovery of
       the following step), and the recovered changes of the backsolved
       elements.  Kept resident under !inmemory: recovery reads exo_z at
       the same point the update loops read x1. */
    solve_real *exo_z=NULL,*bsvals=NULL;
    if(nbselems>0) {
      exo_z= (solve_real *) calloc (nvarele,sizeof(solve_real));
      bsvals= (solve_real *) calloc (nbselems,sizeof(solve_real));
    }
    /* exogenous columns run 0..nexo-1; under heavy condensation nexo can
       exceed VecSize, so vece and B's columns span BSize (see
       shock_vec_set_sizes / shock_mat_set_sizes for the layout rule) */
    PetscInt BSize;
    BSize=(PetscInt)(nvarele-VecSize-nbselems);      /* nexo */
    BSize=(BSize>VecSize)?BSize:VecSize;
    for(subindx=0; subindx<subints; subindx++) {
      /* a new subinterval re-bases the data: the NDBBD cut is re-probed
         there (reuse is within a subinterval, across its passes/steps) */
      ndbbd_cut_cache_free();
      for(sol=0; sol<maxsol; sol++) {
        if(sol==0)nsteps=steps1;
        /* llround: truncation could drop a step when the ratio is not
           an integer multiple (Euler permits any increasing steps) */
        if(sol==1) nsteps=(int)llround(steps1*step_ratio2);
        if(sol==2) nsteps=(int)llround(steps1*step_ratio3);
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
              shock_vec_set_sizes(vece,nesteddbbd,localsize,VecSize,BSize);
              VecSetOption(vece, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
            }
            if(sol==0)for(i=0; i<ncofele; i++) {
                elem_vals[i].initial=elem_vals[i].value;
              }
            else for(i=0; i<ncofele; i++) {
                elem_vals[i].value=elem_vals[i].initial;
              }
            updates_path_restart();
            elem_vals1=elem_vals+ncofele;
            for(i=0; i<nvar; i++) {
              if(vars[i].change_real) {
                for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
                  if(CL_EXO(tindx1)) {
                    if(sol==0) {
                      elem_vals1[tindx1].initial=elem_vals1[tindx1].value;
                    }
                    else {
                      elem_vals1[tindx1].value=elem_vals1[tindx1].initial;
                    }
                    elem_vals1[tindx1].substep_base=CL_SHOCK(tindx1)/nsteps;
                    VecSetValue(vece,closure_vals[tindx1].exo_index,elem_vals1[tindx1].substep_base,INSERT_VALUES);
                    if(exo_z!=NULL)exo_z[tindx1]=elem_vals1[tindx1].substep_base;
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
                  if(CL_EXO(tindx1)) {
                    if(sol==0) {
                      elem_vals1[tindx1].initial=elem_vals1[tindx1].value;
                    }
                    else {
                      elem_vals1[tindx1].value=elem_vals1[tindx1].initial;
                    }
                    temp2=CL_SHOCK(tindx1);//subints;
                    elem_vals1[tindx1].substep_base=(100+(subindx+1)*temp2)/(100+subindx*temp2)-1;//ha_cgeshock[ha_var[i].begadd+j].ShockVal/nsteps;//(exp(log(1+ha_cgeshock[ha_var[i].begadd+j].ShockVal/100)/nsteps)-1)*100;
                    elem_vals1[tindx1].substep_base*=vpercents;//nsteps*100;
                    VecSetValue(vece,closure_vals[tindx1].exo_index,elem_vals1[tindx1].substep_base,INSERT_VALUES);
                    if(exo_z!=NULL)exo_z[tindx1]=elem_vals1[tindx1].substep_base;
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
            if ( (tempvar = teems_fopen(tempfilenam, "wb")) == NULL ) {
              errmsg("Error: cannot open %s for writing\n",tempfilenam);
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
            if ( (tempvar = teems_fopen(tempfilenam, "wb")) == NULL ) {
              errmsg("Error: cannot open %s for writing\n",tempfilenam);
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
          jac_mat_prealloc(A,"Jacobian",nohsl,rank==rank_hsl,Iend-Istart,dnz,dnnz,onz,onnz);

          if(nohsl) {
            MatCreate(PETSC_COMM_WORLD,&B);
          }
          else {
            MatCreate(PETSC_COMM_SELF,&B);
          }
          shock_mat_set_sizes(B,nesteddbbd,localsize,VecSize,BSize);
          if(nohsl) {
            MatSetType(B,MATMPIAIJ);
          }
          else {
            MatSetType(B,MATSEQAIJ);
          }
          jac_mat_prealloc(B,"exogenous block",nohsl,rank==rank_hsl,Iend-Istart,dnzB,dnnzB,onzB,onnzB);
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
            if ( (tempvar = teems_fopen(tempfilenam, "wb")) == NULL ) {
              errmsg("Error: cannot open %s for writing\n",tempfilenam);
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
            if ( (tempvar = teems_fopen(tempfilenam, "wb")) == NULL ) {
              errmsg("Error: cannot open %s for writing\n",tempfilenam);
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
          /* vecb spans the equation rows (VecSize); vece may be wider (BSize) */
          {
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
      /* the (N)DBBD back-solves probe vecb with -1 sentinels for rows
         other ranks own; VecDuplicate(vece) used to inherit this option
         before vecb was created explicitly (the un-ignored -1 made
         VecGetValues error out and drop the rest of the batch) */
      VecSetOption(vecb, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
    }
          if(rank==rank_hsl) {
            ierr = MatMult(B,vece,vecb);
            CHKERRQ(ierr);
          }
          sub_step_rhs(B,vece,VecSize,matsol,rank,rank_hsl);
          ierr = VecDestroy(&vece);
          CHKERRQ(ierr);
          ierr = VecAssemblyBegin(vecb);
          CHKERRQ(ierr);
          ierr = VecAssemblyEnd(vecb);
          CHKERRQ(ierr);
          ierr = MatDestroy(&B);
          CHKERRQ(ierr);
          fh_step_begin(A,vecb,matsol,mc66,VecSize,rank,rank_hsl);
          if(matsol>=MM_DBBD) {
            int *row_order= (int *) calloc (VecSize,sizeof(int));
            int *col_order= (int *) calloc (VecSize,sizeof(int));
            int *block_sizes= (int *) calloc (ndblock,sizeof(int));
            time(&timestr);

            if(matsol==MM_DBBD) {
              dbbd_order(A,VecSize,mpisize,rank,Istart,Iend,nvarele,eq_addr,row_order,col_order,ndblock,block_sizes,countvarintra1,counteq,counteqnoadd,laA,ma48_cntl4);
              x1=realloc (x1,VecSize*sizeof(solve_real));
              dbbd_solve(A,vecb,x1,VecSize,mpisize,rank,Istart,Iend,row_order,col_order,ndblock,block_sizes,countvarintra1,counteq,counteqnoadd,laA,laD);//,iter
            }

            if(matsol==MM_NDBBD) {
              presol=1;
              memcpy(counteq,counteqs,(ndblock+1)*sizeof(offset_t));
              memcpy(counteqnoadd,counteqnoadds,(ndblock)*sizeof(offset_t));
              memcpy(countvarintra1,countvarintra1s,(ndblock+1)*sizeof(offset_t));
              teems_stage_mark("order-presolve");
              ndbbd_order_presolve(A,VecSize,mpisize,rank,Istart,Iend,nreg,ntime,nvarele,eq_addr,row_order,col_order,ndblock,block_sizes,countvarintra1,counteq,counteqnoadd,laA,laDi,ma48_cntl4,ndbbddrank1,presol);
              teems_stage_mark("presolve");
              ndbbd_presolve(A,vecb,x1,VecSize,mpisize,rank,Istart,Iend,row_order,col_order,ndblock,nreg,ntime,block_sizes,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,ma48_cntl4,presol);//,iter
              presol=0;
              teems_stage_mark("order");
              ndbbd_order(A,VecSize,mpisize,rank,Istart,Iend,nreg,ntime,nvarele,eq_addr,row_order,col_order,ndblock,block_sizes,countvarintra1,counteq,counteqnoadd,laA,laDi,ma48_cntl4,ndbbddrank1,presol);
              x1=realloc (x1,VecSize*sizeof(solve_real));
              teems_stage_mark("solve");
              ndbbd_solve(A,vecb,x1,VecSize,mpisize,rank,Istart,Iend,row_order,col_order,ndblock,nreg,ntime,block_sizes,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,ma48_cntl4,presol);//,iter
              teems_stage_report("ndbbd");
            }

            time(&timeend);
            MPI_Barrier(PETSC_COMM_WORLD);
            ierr = PetscGetCPUTime(&time1);
            if(verbosity>=1){ierr = PetscPrintf(PETSC_COMM_WORLD,"One step solution %f\n",time1-time0);}
            if(rank==0)logmsg(1,"Step time %.2f s\n",difftime(timeend,timestr));
            teems_rss_probe("step");
            free(row_order);
            free(col_order);
            free(block_sizes);
            MPI_Barrier(PETSC_COMM_WORLD);
          }
          else {
            if(matsol==MM_SBBD&&fastrefac) {
              x1=realloc (x1,VecSize*sizeof(solve_real));
              sbbd_fastrefac_solve(&A,&vecb,VecSize,rank,rank_hsl,indata,fcomm,counteq,countvarintra1,x1);
            }
            else if(matsol==MM_SBBD) {
              if(mc66==0) {
                x1=realloc (x1,VecSize*sizeof(solve_real));
                sbbd_csr_solve(&A,&vecb,VecSize,rank,rank_hsl,ptx,fcomm,counteq,countvarintra1,x1);
                ierr = PetscGetCPUTime(&time1);
                CHKERRQ(ierr);
                if(verbosity>=1){ierr = PetscPrintf(PETSC_COMM_WORLD,"LU time %f\n",time1-time0);}
                CHKERRQ(ierr);
                ierr = PetscGetCPUTime(&time0);
                CHKERRQ(ierr);
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
              if(rank==rank_hsl)probe_onfail_scope_set_coo(irn,jcn,values,NULL,count,VecSize,VecSize,"SBBD system (MP48)",NULL,NULL);
              if(mc66!=0)spec48_single_(ptx,irn,jcn,b1,values,x1,neleperrow,ai1,&fcomm);
              if(mc66==0)spec48_nomc66_(ptx,jcn,b1,values,x1,neleperrow,&fcomm,counteq,countvarintra1);
              solve_x_check(x1,VecSize,rank==rank_hsl);
              residual_note_skipped(rank,rank_hsl);
              probe_onfail_scope_clear();
              free(irn);
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
            }
            else if(fastrefac) {
              x1=realloc (x1,VecSize*sizeof(solve_real));
              ierr = PetscGetCPUTime(&time1);
              CHKERRQ(ierr);
              if(verbosity>=1){ierr = PetscPrintf(PETSC_COMM_WORLD,"Prepare time %f\n",time1-time0);}
              CHKERRQ(ierr);
              ierr = PetscGetCPUTime(&time0);
              CHKERRQ(ierr);
              if(rank==rank_hsl) {
                VecGetArray(vecb,&vals);
                lu_fastrefac_solve(A,VecSize,laA,vals,x1);
              }
              ierr = MatDestroy(&A);
              CHKERRQ(ierr);
              ierr = VecDestroy(&vecb);
              CHKERRQ(ierr);
              ierr = PetscGetCPUTime(&time1);
              CHKERRQ(ierr);
              if(verbosity>=1){ierr = PetscPrintf(PETSC_COMM_WORLD,"LU time %f\n",time1-time0);}
              CHKERRQ(ierr);
              ierr = PetscGetCPUTime(&time0);
              CHKERRQ(ierr);
            }
            else {
              if(rank==rank_hsl)VecGetArray(vecb,&vals);
              x1=realloc (x1,VecSize*sizeof(solve_real));
              ierr = PetscGetCPUTime(&time1);
              CHKERRQ(ierr);
              if(verbosity>=1){ierr = PetscPrintf(PETSC_COMM_WORLD,"Prepare time %f\n",time1-time0);}
              CHKERRQ(ierr);
              ierr = PetscGetCPUTime(&time0);
              CHKERRQ(ierr);
              /* A stays live through the factorize for the on-failure
                 diagnosis and the workspace-growth re-staging (staged
                 COO is MA48 workspace) */
              if(rank==rank_hsl) {
                probe_onfail_scope_set(A,VecSize,VecSize,"condensed system",-1,NULL,NULL,0,0,0,0);
                lu_grow_solve(A,VecSize,laA,vals,x1);
                probe_onfail_scope_clear();
              }
              ierr = MatDestroy(&A);
              CHKERRQ(ierr);
              ierr = VecDestroy(&vecb);
              CHKERRQ(ierr);
              ierr = PetscGetCPUTime(&time1);
              CHKERRQ(ierr);
              if(verbosity>=1){ierr = PetscPrintf(PETSC_COMM_WORLD,"LU time %f\n",time1-time0);}
              CHKERRQ(ierr);
              ierr = PetscGetCPUTime(&time0);
              CHKERRQ(ierr);
            }
            teems_rss_probe("step");
          }
          fh_step_end(VecSize,x1,rank,rank_hsl,mpisize);
          if(rank==rank_hsl) {
            if(!inmemory){
            if ((tempvar = teems_fopen(tempfilenam, "rb")) == NULL) {
              errmsg("Error: cannot open %s for reading\n",tempfilenam);
            }
            *closure_vals2=(closure_entry*)realloc (*closure_vals2,(nvarele)*sizeof(closure_entry));
            scratch_read(*closure_vals2, sizeof(closure_entry),nvarele,tempvar,tempfilenam);
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
            if ((tempvar = teems_fopen(tempfilenam, "rb")) == NULL) {
              errmsg("Error: cannot open %s for reading\n",tempfilenam);
            }
            *elem_vals2=(elem_value*)realloc (*elem_vals2,(ncofele+nvarele)*sizeof(elem_value));
            scratch_read(*elem_vals2, sizeof(elem_value),ncofele+nvarele,tempvar,tempfilenam);
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
            if ((tempvar = teems_fopen(tempfilenam, "rb")) == NULL) {
              errmsg("Error: cannot open %s for reading\n",tempfilenam);
            }
            clag1=realloc (clag1,(nvarele)*sizeof(solve_real));
            scratch_read(clag1, sizeof(solve_real),nvarele,tempvar,tempfilenam);
            fclose(tempvar);
            remove(tempfilenam);
            }

            if(!inmemory){
            strcpy(tempfilenam,scratch_dir);
            strcat(tempfilenam,"_tempvarchange");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ((tempvar = teems_fopen(tempfilenam, "rb")) == NULL) {
              errmsg("Error: cannot open %s for reading\n",tempfilenam);
            }
            varchange=realloc (varchange,(nvarele)*sizeof(solve_real));
            scratch_read(varchange, sizeof(solve_real),nvarele,tempvar,tempfilenam);
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
          shock_vec_set_sizes(vece,nesteddbbd,localsize,VecSize,BSize);
          VecSetOption(vece, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
          elem_vals1=elem_vals+ncofele;
          /* recover the backsolved elements from their defining equations
             with this step's solution, before any update reads a
             variable's change (GEMPACK 14.1.3) */
          if(rank==rank_hsl&&nbselems>0)backsolve_recover(tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele,closure_vals,x1,exo_z,bsvals);
          if(teems_sub_active&&rank==0)sub_update(stepcount==0?SUB_FIRST:(euler?SUB_EULER:SUB_LEAP),tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele,closure_vals,nvarele,exo_z,varchange);
          if(stepcount==0) {
            for(i=0; i<nvar; i++) {
              if(vars[i].change_real) {
                for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
                  if(CL_EXO(tindx1)) {
                    elem_vals1[tindx1].value+=elem_vals1[tindx1].substep_base;
                    varchange[tindx1]=elem_vals1[tindx1].substep_base;
                    VecSetValue(vece,closure_vals[tindx1].exo_index,elem_vals1[tindx1].substep_base,INSERT_VALUES);
                    if(exo_z!=NULL)exo_z[tindx1]=elem_vals1[tindx1].substep_base;
                  }
                  else if(CL_BS(tindx1)) {
                    varchange[tindx1]=bsvals[closure_vals[tindx1].exo_index];
                    elem_vals1[tindx1].value+=bsvals[closure_vals[tindx1].exo_index];
                    elem_vals1[tindx1].substep_base=bsvals[closure_vals[tindx1].exo_index];
                    clag1[tindx1]=0;
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
                  if(CL_EXO(tindx1)) {
                    varchange[tindx1]=elem_vals1[tindx1].substep_base;
                    elem_vals1[tindx1].value*=(1+elem_vals1[tindx1].substep_base/100);
                    VecSetValue(vece,closure_vals[tindx1].exo_index,elem_vals1[tindx1].substep_base/(1+elem_vals1[tindx1].substep_base/100),INSERT_VALUES);
                    if(exo_z!=NULL)exo_z[tindx1]=elem_vals1[tindx1].substep_base/(1+elem_vals1[tindx1].substep_base/100);
                  }
                  else if(CL_BS(tindx1)) {
                    varchange[tindx1]=bsvals[closure_vals[tindx1].exo_index];
                    elem_vals1[tindx1].substep_base=bsvals[closure_vals[tindx1].exo_index];
                    elem_vals1[tindx1].value*=(1+elem_vals1[tindx1].substep_base/100);
                    clag1[tindx1]=0;
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
          else if(euler) {
            /* forward Euler substep: accumulate this step's change on the
               current state; clag1 stays unused (no leapfrog history) */
            for(i=0; i<nvar; i++) {
              if(vars[i].change_real) {
                for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
                  if(CL_EXO(tindx1)) {
                    elem_vals1[tindx1].value+=elem_vals1[tindx1].substep_base;
                    varchange[tindx1]+=elem_vals1[tindx1].substep_base;
                    VecSetValue(vece,closure_vals[tindx1].exo_index,elem_vals1[tindx1].substep_base,INSERT_VALUES);
                    if(exo_z!=NULL)exo_z[tindx1]=elem_vals1[tindx1].substep_base;
                  }
                  else if(CL_BS(tindx1)) {
                    varchange[tindx1]+=bsvals[closure_vals[tindx1].exo_index];
                    elem_vals1[tindx1].substep_base=bsvals[closure_vals[tindx1].exo_index];
                    elem_vals1[tindx1].value+=bsvals[closure_vals[tindx1].exo_index];
                  }
                  else {
                    varchange[tindx1]+=x1[closure_vals[tindx1].exo_index];
                    elem_vals1[tindx1].substep_base=x1[closure_vals[tindx1].exo_index];
                    elem_vals1[tindx1].value+=x1[closure_vals[tindx1].exo_index];
                  }
                }
              }
              else {
                for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
                  if(CL_EXO(tindx1)) {
                    temp2=CL_SHOCK(tindx1);//subints;
                    temp1=(100+(subindx+1)*temp2)/(100+subindx*temp2)-1;
                    temp1*=vpercents;
                    elem_vals1[tindx1].substep_base=temp1/(1+varchange[tindx1]/100);
                    varchange[tindx1]+=temp1;//*(1+ha_cofvar[ncofele+ha_var[i].begadd+j].varchange/100)
                    elem_vals1[tindx1].value=(1+varchange[tindx1]/100)*elem_vals1[tindx1].initial;
                    VecSetValue(vece,closure_vals[tindx1].exo_index,temp1/(1+varchange[tindx1]/100),INSERT_VALUES);
                    if(exo_z!=NULL)exo_z[tindx1]=temp1/(1+varchange[tindx1]/100);
                  }
                  else if(CL_BS(tindx1)) {
                    /* compound the per-step percent change onto the
                       cumulative one before rebasing on the initial */
                    varchange[tindx1]+=bsvals[closure_vals[tindx1].exo_index]*(100+varchange[tindx1])/100;
                    elem_vals1[tindx1].substep_base=bsvals[closure_vals[tindx1].exo_index];
                    elem_vals1[tindx1].value=(100+varchange[tindx1])/100*elem_vals1[tindx1].initial;
                  }
                  else {
                    varchange[tindx1]+=x1[closure_vals[tindx1].exo_index]*(100+varchange[tindx1])/100;
                    elem_vals1[tindx1].substep_base=x1[closure_vals[tindx1].exo_index];
                    elem_vals1[tindx1].value=(100+varchange[tindx1])/100*elem_vals1[tindx1].initial;
                  }
                }
              }
            }
          }
          else {
            for(i=0; i<nvar; i++) {
              if(vars[i].change_real) {
                for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
                  if(CL_EXO(tindx1)) {
                    elem_vals1[tindx1].value+=elem_vals1[tindx1].substep_base;
                    varchange[tindx1]+=elem_vals1[tindx1].substep_base;
                    VecSetValue(vece,closure_vals[tindx1].exo_index,elem_vals1[tindx1].substep_base,INSERT_VALUES);
                    if(exo_z!=NULL)exo_z[tindx1]=elem_vals1[tindx1].substep_base;
                  }
                  else if(CL_BS(tindx1)) {
                    temp1=elem_vals1[tindx1].value;
                    varchange[tindx1]=clag1[tindx1]+2*bsvals[closure_vals[tindx1].exo_index];
                    elem_vals1[tindx1].substep_base=bsvals[closure_vals[tindx1].exo_index];
                    elem_vals1[tindx1].value=clag1[tindx1]+2*bsvals[closure_vals[tindx1].exo_index];
                    clag1[tindx1]=temp1;
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
                  if(CL_EXO(tindx1)) {
                    temp2=CL_SHOCK(tindx1);//subints;
                    temp1=(100+(subindx+1)*temp2)/(100+subindx*temp2)-1;
                    temp1*=vpercents;
                    elem_vals1[tindx1].substep_base=temp1/(1+varchange[tindx1]/100);
                    varchange[tindx1]+=temp1;//*(1+ha_cofvar[ncofele+ha_var[i].begadd+j].varchange/100)
                    elem_vals1[tindx1].value=(1+varchange[tindx1]/100)*elem_vals1[tindx1].initial;
                    VecSetValue(vece,closure_vals[tindx1].exo_index,temp1/(1+varchange[tindx1]/100),INSERT_VALUES);
                    if(exo_z!=NULL)exo_z[tindx1]=temp1/(1+varchange[tindx1]/100);
                  }
                  else if(CL_BS(tindx1)) {
                    temp1=varchange[tindx1];
                    varchange[tindx1]=clag1[tindx1]+2*bsvals[closure_vals[tindx1].exo_index]*(100+temp1)/100;
                    elem_vals1[tindx1].substep_base=bsvals[closure_vals[tindx1].exo_index];
                    elem_vals1[tindx1].value=(100+varchange[tindx1])/100*elem_vals1[tindx1].initial;
                    clag1[tindx1]=temp1;
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
          if(teems_sub_active&&rank==0)sub_exo_sync(varchange);
          free(x1);
          x1=NULL;
          MPI_Barrier(PETSC_COMM_WORLD);
          ierr = VecAssemblyBegin(vece);
          CHKERRQ(ierr);
          MPI_Barrier(PETSC_COMM_WORLD);
          ierr = VecAssemblyEnd(vece);
          CHKERRQ(ierr);
          if(rank==rank_hsl) {
            if(stepcount==0||euler) {
              updates_apply(tabfile,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,0);
            }
            else {
              updates_apply(tabfile,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,1);
            }
            strcpy(commsyntax,"formula");
            IsIni=false;
            statements_execute(tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,IsIni,0);
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
        /* Gragg's terminal smoothing pass: one more Jacobian build and
           solve at the final state, then the half-sum correction. Euler
           and midpoint have no such pass — their accumulated varchange IS
           the solution, so skip straight to the state reset (nothing was
           spilled). */
        if(!euler&&!midpoint) {
        if(rank==rank_hsl) {


          if(!inmemory){
          strcpy(tempfilenam,scratch_dir);
          strcat(tempfilenam,"_tempclag1");
          sprintf(tempchar, "%d",rank);
          strcat(tempfilenam,tempchar);
          strcat(tempfilenam,".bin");
          if ( (tempvar = teems_fopen(tempfilenam, "wb")) == NULL ) {
            errmsg("Error: cannot open %s for writing\n",tempfilenam);
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
          if ( (tempvar = teems_fopen(tempfilenam, "wb")) == NULL ) {
            errmsg("Error: cannot open %s for writing\n",tempfilenam);
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
        jac_mat_prealloc(A,"Jacobian",nohsl,rank==rank_hsl,Iend-Istart,dnz,dnnz,onz,onnz);

        if(nohsl) {
          MatCreate(PETSC_COMM_WORLD,&B);
        }
        else {
          MatCreate(PETSC_COMM_SELF,&B);
        }
        shock_mat_set_sizes(B,nesteddbbd,localsize,VecSize,BSize);
        if(nohsl) {
          MatSetType(B,MATMPIAIJ);
        }
        else {
          MatSetType(B,MATSEQAIJ);
        }
        jac_mat_prealloc(B,"exogenous block",nohsl,rank==rank_hsl,Iend-Istart,dnzB,dnnzB,onzB,onnzB);

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
          if ( (tempvar = teems_fopen(tempfilenam, "wb")) == NULL ) {
            errmsg("Error: cannot open %s for writing\n",tempfilenam);
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
          if ( (tempvar = teems_fopen(tempfilenam, "wb")) == NULL ) {
            errmsg("Error: cannot open %s for writing\n",tempfilenam);
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
        /* vecb spans the equation rows (VecSize); vece may be wider (BSize) */
        {
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
      /* the (N)DBBD back-solves probe vecb with -1 sentinels for rows
         other ranks own; VecDuplicate(vece) used to inherit this option
         before vecb was created explicitly (the un-ignored -1 made
         VecGetValues error out and drop the rest of the batch) */
      VecSetOption(vecb, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
    }
        if(rank==rank_hsl) {
          ierr = MatMult(B,vece,vecb);
          CHKERRQ(ierr);
        }
        sub_step_rhs(B,vece,VecSize,matsol,rank,rank_hsl);
        ierr = VecDestroy(&vece);
        CHKERRQ(ierr);
        ierr = VecAssemblyBegin(vecb);
        CHKERRQ(ierr);
        ierr = VecAssemblyEnd(vecb);
        CHKERRQ(ierr);
        ierr = MatDestroy(&B);
        CHKERRQ(ierr);
        fh_step_begin(A,vecb,matsol,mc66,VecSize,rank,rank_hsl);
        if(matsol>=MM_DBBD) {
          int *row_order= (int *) calloc (VecSize,sizeof(int));
          int *col_order= (int *) calloc (VecSize,sizeof(int));
          int *block_sizes= (int *) calloc (ndblock,sizeof(int));
          time(&timestr);

          if(matsol==MM_DBBD) {
            dbbd_order(A,VecSize,mpisize,rank,Istart,Iend,nvarele,eq_addr,row_order,col_order,ndblock,block_sizes,countvarintra1,counteq,counteqnoadd,laA,ma48_cntl4);
            x1=realloc (x1,VecSize*sizeof(solve_real));
            dbbd_solve(A,vecb,x1,VecSize,mpisize,rank,Istart,Iend,row_order,col_order,ndblock,block_sizes,countvarintra1,counteq,counteqnoadd,laA,laD);//,iter
          }

          if(matsol==MM_NDBBD) {
            presol=1;
            memcpy(counteq,counteqs,(ndblock+1)*sizeof(offset_t));
            memcpy(counteqnoadd,counteqnoadds,(ndblock)*sizeof(offset_t));
            memcpy(countvarintra1,countvarintra1s,(ndblock+1)*sizeof(offset_t));
            teems_stage_mark("order-presolve");
            ndbbd_order_presolve(A,VecSize,mpisize,rank,Istart,Iend,nreg,ntime,nvarele,eq_addr,row_order,col_order,ndblock,block_sizes,countvarintra1,counteq,counteqnoadd,laA,laDi,ma48_cntl4,ndbbddrank1,presol);
            teems_stage_mark("presolve");
            ndbbd_presolve(A,vecb,x1,VecSize,mpisize,rank,Istart,Iend,row_order,col_order,ndblock,nreg,ntime,block_sizes,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,ma48_cntl4,presol);//,iter
            presol=0;
            teems_stage_mark("order");
            ndbbd_order(A,VecSize,mpisize,rank,Istart,Iend,nreg,ntime,nvarele,eq_addr,row_order,col_order,ndblock,block_sizes,countvarintra1,counteq,counteqnoadd,laA,laDi,ma48_cntl4,ndbbddrank1,presol);
            x1=realloc (x1,VecSize*sizeof(solve_real));
            teems_stage_mark("solve");
            ndbbd_solve(A,vecb,x1,VecSize,mpisize,rank,Istart,Iend,row_order,col_order,ndblock,nreg,ntime,block_sizes,countvarintra1,counteq,counteqnoadd,laA,laDi,laD,ma48_cntl4,presol);//,iter
            teems_stage_report("ndbbd");
          }

          time(&timeend);
          MPI_Barrier(PETSC_COMM_WORLD);
          ierr = PetscGetCPUTime(&time1);
          if(verbosity>=1){ierr = PetscPrintf(PETSC_COMM_WORLD,"One step solution %f\n",time1-time0);}
          if(rank==0)logmsg(1,"Step time %.2f s\n",difftime(timeend,timestr));
          teems_rss_probe("step");
          free(row_order);
          free(col_order);
          free(block_sizes);
          MPI_Barrier(PETSC_COMM_WORLD);
        }
        else {
          if(matsol==MM_SBBD&&fastrefac) {
            x1=realloc (x1,VecSize*sizeof(solve_real));
            sbbd_fastrefac_solve(&A,&vecb,VecSize,rank,rank_hsl,indata,fcomm,counteq,countvarintra1,x1);
          }
          else if(matsol==MM_SBBD) {
            if(mc66==0) {
              x1=realloc (x1,VecSize*sizeof(solve_real));
              sbbd_csr_solve(&A,&vecb,VecSize,rank,rank_hsl,ptx,fcomm,counteq,countvarintra1,x1);
              ierr = PetscGetCPUTime(&time1);
              CHKERRQ(ierr);
              if(verbosity>=1){ierr = PetscPrintf(PETSC_COMM_WORLD,"LU time %f\n",time1-time0);}
              CHKERRQ(ierr);
              ierr = PetscGetCPUTime(&time0);
              CHKERRQ(ierr);
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
            if(rank==rank_hsl)probe_onfail_scope_set_coo(irn,jcn,values,NULL,count,VecSize,VecSize,"SBBD system (MP48)",NULL,NULL);
            if(mc66!=0)spec48_single_(ptx,irn,jcn,b1,values,x1,neleperrow,ai1,&fcomm);
            if(mc66==0)spec48_nomc66_(ptx,jcn,b1,values,x1,neleperrow,&fcomm,counteq,countvarintra1);
              solve_x_check(x1,VecSize,rank==rank_hsl);
              residual_note_skipped(rank,rank_hsl);
            probe_onfail_scope_clear();
            free(irn);
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
          }
          else if(fastrefac) {
            x1=realloc (x1,VecSize*sizeof(solve_real));
            ierr = PetscGetCPUTime(&time1);
            CHKERRQ(ierr);
            if(verbosity>=1){ierr = PetscPrintf(PETSC_COMM_WORLD,"Prepare time %f\n",time1-time0);}
            CHKERRQ(ierr);
            ierr = PetscGetCPUTime(&time0);
            CHKERRQ(ierr);
            if(rank==rank_hsl) {
              VecGetArray(vecb,&vals);
              lu_fastrefac_solve(A,VecSize,laA,vals,x1);
            }
            ierr = MatDestroy(&A);
            CHKERRQ(ierr);
            ierr = VecDestroy(&vecb);
            CHKERRQ(ierr);
            ierr = PetscGetCPUTime(&time1);
            CHKERRQ(ierr);
            if(verbosity>=1){ierr = PetscPrintf(PETSC_COMM_WORLD,"LU time %f\n",time1-time0);}
            CHKERRQ(ierr);
            ierr = PetscGetCPUTime(&time0);
            CHKERRQ(ierr);
          }
          else {
            if(rank==rank_hsl)VecGetArray(vecb,&vals);
            x1=realloc (x1,VecSize*sizeof(solve_real));
            ierr = PetscGetCPUTime(&time1);
            CHKERRQ(ierr);
            if(verbosity>=1){ierr = PetscPrintf(PETSC_COMM_WORLD,"Prepare time %f\n",time1-time0);}
            CHKERRQ(ierr);
            ierr = PetscGetCPUTime(&time0);
            CHKERRQ(ierr);
            /* A stays live through the factorize for the on-failure
               diagnosis and the workspace-growth re-staging (staged
               COO is MA48 workspace) */
            if(rank==rank_hsl) {
              probe_onfail_scope_set(A,VecSize,VecSize,"condensed system",-1,NULL,NULL,0,0,0,0);
              lu_grow_solve(A,VecSize,laA,vals,x1);
              probe_onfail_scope_clear();
            }
            ierr = MatDestroy(&A);
            CHKERRQ(ierr);
            ierr = PetscGetCPUTime(&time1);
            CHKERRQ(ierr);
            if(verbosity>=1){ierr = PetscPrintf(PETSC_COMM_WORLD,"LU time %f\n",time1-time0);}
            CHKERRQ(ierr);
            ierr = PetscGetCPUTime(&time0);
            CHKERRQ(ierr);
            ierr = VecDestroy(&vecb);
            CHKERRQ(ierr);
          }
          teems_rss_probe("step");
        }
        fh_step_end(VecSize,x1,rank,rank_hsl,mpisize);
        if(rank==rank_hsl) {
          if(!inmemory){
          if ((tempvar = teems_fopen(tempfilenam, "rb")) == NULL) {
            errmsg("Error: cannot open %s for reading\n",tempfilenam);
          }
          *closure_vals2=(closure_entry*)realloc (*closure_vals2,(nvarele)*sizeof(closure_entry));
          scratch_read(*closure_vals2, sizeof(closure_entry),nvarele,tempvar,tempfilenam);
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
          if ((tempvar = teems_fopen(tempfilenam, "rb")) == NULL) {
            errmsg("Error: cannot open %s for reading\n",tempfilenam);
          }
          *elem_vals2=(elem_value*)realloc (*elem_vals2,(ncofele+nvarele)*sizeof(elem_value));
          scratch_read(*elem_vals2, sizeof(elem_value),ncofele+nvarele,tempvar,tempfilenam);
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
          if ((tempvar = teems_fopen(tempfilenam, "rb")) == NULL) {
            errmsg("Error: cannot open %s for reading\n",tempfilenam);
          }
          clag1=realloc (clag1,(nvarele)*sizeof(solve_real));
          scratch_read(clag1, sizeof(solve_real),nvarele,tempvar,tempfilenam);
          fclose(tempvar);
          remove(tempfilenam);
          }

          if(!inmemory){
          strcpy(tempfilenam,scratch_dir);
          strcat(tempfilenam,"_tempvarchange");
          sprintf(tempchar, "%d",rank);
          strcat(tempfilenam,tempchar);
          strcat(tempfilenam,".bin");
          if ((tempvar = teems_fopen(tempfilenam, "rb")) == NULL) {
            errmsg("Error: cannot open %s for reading\n",tempfilenam);
          }
          varchange=realloc (varchange,(nvarele)*sizeof(solve_real));
          scratch_read(varchange, sizeof(solve_real),nvarele,tempvar,tempfilenam);
          fclose(tempvar);
          remove(tempfilenam);
          }

        }
        elem_vals1=elem_vals+ncofele;
        /* recover the backsolved elements for the terminal smoothing
           solve (same pre-update evaluation point as the step loop) */
        if(rank==rank_hsl&&nbselems>0)backsolve_recover(tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele,closure_vals,x1,exo_z,bsvals);
        /* the subtotal columns smooth at the same pre-update point */
        if(teems_sub_active&&rank==0)sub_update(SUB_SMOOTH,tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele,closure_vals,nvarele,exo_z,varchange);
        /* the same terminal smoothing for the updated data of the
           (change)/(explicit) Update targets, (C[n-1] + C[n] + dC)/2 with
           dC from this pass's changes: their pass-end values feed the
           updated-data extrapolation only (the next pass restarts from
           the initial data), and for a change update linear in the
           variables it keeps the extrapolated value equal to the
           variables' own result */
        if(rank==rank_hsl&&updates_path_active()) {
          for(i=0; i<nvar; i++) {
            for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
              if(CL_EXO(tindx1)) {
                if(!vars[i].change_real) {
                  temp2=CL_SHOCK(tindx1);
                  temp1=((100+(subindx+1)*temp2)/(100+subindx*temp2)-1)*vpercents;
                  elem_vals1[tindx1].substep_base=temp1/(1+varchange[tindx1]/100);
                }
              }
              else if(CL_BS(tindx1)) elem_vals1[tindx1].substep_base=bsvals[closure_vals[tindx1].exo_index];
              else elem_vals1[tindx1].substep_base=x1[closure_vals[tindx1].exo_index];
            }
          }
          updates_apply(tabfile,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,2);
        }
        for(i=0; i<nvar; i++) {
          if(vars[i].change_real) {
            for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
              if(CL_EXO(tindx1)) {
                elem_vals1[tindx1].value=0;
              }
              else if(CL_BS(tindx1)) {
                varchange[tindx1]=0.5*(varchange[tindx1]+clag1[tindx1]+bsvals[closure_vals[tindx1].exo_index]);
                elem_vals1[tindx1].value=0;
                clag1[tindx1]=0;
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
              if(CL_EXO(tindx1)) {
                elem_vals1[tindx1].value=0;
              }
              else if(CL_BS(tindx1)) {
                varchange[tindx1]=0.5*(varchange[tindx1]+clag1[tindx1]+bsvals[closure_vals[tindx1].exo_index]*(1+varchange[tindx1]/100));
                elem_vals1[tindx1].value=0;
                clag1[tindx1]=0;
              }
              else {
                varchange[tindx1]=0.5*(varchange[tindx1]+clag1[tindx1]+x1[closure_vals[tindx1].exo_index]*(1+varchange[tindx1]/100));
                elem_vals1[tindx1].value=0;//ha_cofvar[tindx2].varval*varchange[tindx1]/100;
                clag1[tindx1]=0;
              }
            }
          }
        }
        }
        else {
          /* Euler/midpoint: the last substep created vece for a next fill that
             never comes; the smoothing pass destroyed it on the Gragg
             path */
          ierr = VecDestroy(&vece);
          CHKERRQ(ierr);
          /* reset the variable state for the next solution pass
             (mirrors the zeroing the smoothing loop performs) */
          elem_vals1=elem_vals+ncofele;
          for(i=0; i<nvar; i++) {
            for(tindx1=vars[i].offset; tindx1<vars[i].nelem+vars[i].offset; tindx1++) {
              elem_vals1[tindx1].value=0;
              clag1[tindx1]=0;
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
            if ((tempvar = teems_fopen(tempfilenam, "rb")) == NULL) {
              errmsg("Error: cannot open %s for reading\n",tempfilenam);
            }
            *xcf2=(solve_real*)realloc (*xcf2,(nvarele)*sizeof(solve_real));
            xcf=*xcf2;
            scratch_read(xcf, sizeof(solve_real),nvarele,tempvar,tempfilenam);
            fclose(tempvar);
            remove(tempfilenam);
            }

            if(!inmemory){
            strcpy(tempfilenam,scratch_dir);
            strcat(tempfilenam,"_tempxc12");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ((tempvar = teems_fopen(tempfilenam, "rb")) == NULL) {
              errmsg("Error: cannot open %s for reading\n",tempfilenam);
            }
            xc12=realloc (xc12,(nvarele)*sizeof(solve_real));
            scratch_read(xc12, sizeof(solve_real),nvarele,tempvar,tempfilenam);
            fclose(tempvar);
            remove(tempfilenam);
            }
            if(!inmemory){
            strcpy(tempfilenam,scratch_dir);
            strcat(tempfilenam,"_tempxc24");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ((tempvar = teems_fopen(tempfilenam, "rb")) == NULL) {
              errmsg("Error: cannot open %s for reading\n",tempfilenam);
            }
            xc24=realloc (xc24,(nvarele)*sizeof(solve_real));
            scratch_read(xc24, sizeof(solve_real),nvarele,tempvar,tempfilenam);
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
            if ((tempvar = teems_fopen(tempfilenam, "rb")) == NULL) {
              errmsg("Error: cannot open %s for reading\n",tempfilenam);
            }
            scratch_read(xc0, sizeof(solve_real),nvarele,tempvar,tempfilenam);
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
          {
            /* Richardson weights for the current solution pass. q is the
               step-ratio power matching the method's truncation error
               series: h^2 for Gragg (even-power expansion, Pearson 1991
               Thm 6.1), h for Euler. */
            double q2=euler?step_ratio2:step_ratio2*step_ratio2;
            double q3=euler?step_ratio3:step_ratio3*step_ratio3;
            if(sol==0) {
              extrap_w1=1.0/(q2-1.0);
              extrap_w2=1.0/(1-q2)/(1.0-q3);
              /* single pass: the pass itself is the result (weight 1) */
              if(teems_single_run) { extrap_w1=0.0; extrap_w2=1.0; }
              /* two solutions: (q c2 - c1)/(q - 1) */
              if(teems_two_run) extrap_w2=1.0/(1.0-q2);
            }
            if(sol==1) {
              extrap_w1=q2/(q2-1.0);
              extrap_w2=q2/(q3-q2);
              extrap_w3=q2*q2/(q2-q3)/(1.0-q2);
              if(teems_two_run) { extrap_w2=0.0; extrap_w3=-q2/(q2-1.0); }
            }
            if(sol==2) {
              extrap_w2=q3/(q3-q2);
              extrap_w3=q3*q3/(q2-q3)/(1.0-q3);
            }
          }
          /* updated data (manual 26.2): the (change)/(explicit) Update
             targets' values at the end of this pass enter the
             extrapolation with the weight the variables' final result
             gives this pass (sol 0: w2, sol 1: -w3, sol 2: w3) */
          if(rank==0&&maxsol==3&&subindx==subints-1)xac_pass(sol,subindx,subints,vars,nvar,nvarele,varchange,xcf);
          if(teems_sup&&rank==0&&subints==1&&(teems_sup==2||sol==maxsol-1)) {
            static const char *const udext[3]={".ud5",".ud6",".ud7"};
            coefficients_dump_phase(teems_sol_stem,udext[sol],5+sol,ncof,ncofele,elem_vals);
          }
          if(updates_path_active())updates_path_accumulate(coefs,ncof,elem_vals,(sol==0)?extrap_w2:((sol==1)?-extrap_w3:extrap_w3),sol==0);
          if(convrule&&updates_path_active()) {
            updates_path_conv_store(sol,coefs,ncof,elem_vals);
            if(sol==2) {
              double q2,q3;
              conv_ratios(euler,&q2,&q3);
              updates_path_conv_apply(q2,q3);
            }
          }
          if(subindx>0) {
            if(sol==0)for(i=0; i<nvarele; i++) xc0[i]=1+xcf[i]/100;//if(i==1287)printf("sol!!!!!!!!!!!!!!!!!! %d step %d xc %lf xc0 %lf k %d\n",sol,stepcount,1.0+xc[k]/100,xc0[i],i);}
            if(sol==0) {
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
              for(i=0; i<nvarele; i++) {
                xc12[i]=-varchange[i]*extrap_w1;
                xcf[i]+=varchange[i]*extrap_w2;
              }
            }
            if(sol==1) {
              for(i=0; i<nvarele; i++) {
                xc24[i]=-varchange[i]*extrap_w2;
                xc12[i]+=varchange[i]*extrap_w1;
                xcf[i]-=varchange[i]*extrap_w3;
              }
            }
            if(sol==2) {
              for(i=0; i<nvarele; i++) {
                xc24[i]+=varchange[i]*extrap_w2;
                xcf[i]+=varchange[i]*extrap_w3;
              }
            }
          }
          if(convrule&&rank==rank_hsl) {
            conv_store(sol,nvarele,varchange,sol==0?extrap_w2:(sol==1?-extrap_w3:extrap_w3));
            if(sol==2) conv_apply(vars,nvar,subindx,xc0,xcf,euler);
          }
          if(teems_sub_active&&rank==0)sub_pass_end(sol,subindx,vars,nvar,nvarele,xc0);

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
          /* product updates: exact one-shot form from the extrapolated
             totals; (change)/(explicit) targets: their extrapolated path
             values (SW5) -- not a one-shot recomputation, which is wrong
             for nonlinear change, explicit and counter updates */
          teems_upd_pathuse=1;
          updates_apply_product(tabfile,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele);
          teems_upd_pathuse=0;
          strcpy(commsyntax,"formula");
          IsIni=false;
          statements_execute(tabfile,commsyntax,sets,nset,set_elems,coefs,ncof,vars,nvar,elem_vals,ncofele+nvarele,ncofele,IsIni,0);
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
            if ( (tempvar = teems_fopen(tempfilenam, "wb")) == NULL ) {
              errmsg("Error: cannot open %s for writing\n",tempfilenam);
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
            if ( (tempvar = teems_fopen(tempfilenam, "wb")) == NULL ) {
              errmsg("Error: cannot open %s for writing\n",tempfilenam);
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
            if ((tempvar = teems_fopen(tempfilenam, "rb")) == NULL) {
              errmsg("Error: cannot open %s for reading\n",tempfilenam);
            }
            xc124=realloc (xc124,nvarele*sizeof(int));
            scratch_read(xc124, sizeof(int),nvarele,tempvar,tempfilenam);
            fclose(tempvar);
            remove(tempfilenam);
            }
            }

      /* the figures compare the two two-pass extrapolations, which a
         single- or two-solution run does not have (xc24 is unset) */
      if(maxsol==3) for(i=0; i<nvarele; i++){
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
        /* in double, clamped at the last bucket: two passes of different
           magnitude (one near zero) overflowed int, and the undefined
           conversion could land in the 5-digit bucket */
        {
          double dj=fabs(floor((xc12[i]-xc24[i])*100000));
          j=(dj<10000)?(int)dj:10000;
        }
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
      if(rank==0&&maxsol==3&&subindx==subints-1)xac_codes(nvarele,xc124);
      if(subindx!=subints-1){
            if(!inmemory){
            strcpy(tempfilenam,scratch_dir);
            strcat(tempfilenam,"_tempxc124");
            sprintf(tempchar, "%d",rank);
            strcat(tempfilenam,tempchar);
            strcat(tempfilenam,".bin");
            if ( (tempvar = teems_fopen(tempfilenam, "wb")) == NULL ) {
              errmsg("Error: cannot open %s for writing\n",tempfilenam);
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
            if ( (tempvar = teems_fopen(tempfilenam, "wb")) == NULL ) {
              errmsg("Error: cannot open %s for writing\n",tempfilenam);
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
            if ( (tempvar = teems_fopen(tempfilenam, "wb")) == NULL ) {
              errmsg("Error: cannot open %s for writing\n",tempfilenam);
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
            errmsg("Error: cannot open %s for writing: %s (the solver runs as uid %d)\n",solchar,strerror(errno),(int)getuid());
            return 1;
          }
          if(xcf!=NULL)fwrite(xcf, sizeof(solve_real),nvarele, solution);
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
    if(teems_single_run) {
      if(rank==0)printf("Accuracy estimates: not available for a single-pass run (-single_run 1; they need three multi-step solutions)\n");
    }
    else if(teems_two_run) {
      if(rank==0)printf("Accuracy estimates: not available from two solutions (-two_run 1; they need three multi-step solutions, manual 26.2.3)\n");
    }
    else if(rank==0)printf("Accurate at 6 digits        %ld\nAccurate at 5 digits        %ld\nAccurate at 4 digits        %ld\nAccurate at 3 digits        %ld\nAccurate at 2 digits        %ld\nAccurate at 1 digit or none %ld\n",precis[5],precis[4],precis[3],precis[2],precis[1],precis[0]);
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
    free(exo_z);
    free(bsvals);
    gettimeofday(&endtime, NULL);
    if(convrule&&rank==rank_hsl) {
      logmsg(1,"Convergence rule (-convrule 1, manual 26.2.5): %ld component result(s) extrapolated, %ld averaged (results very close or near zero), %ld taken from the %d-step solution (poor convergence)\n",(long)conv_n[0],(long)conv_n[1],(long)conv_n[2],(int)llround(steps1*step_ratio3));
      for(i=0; i<3; i++) { free(conv_c[i]); conv_c[i]=NULL; }
    }
    if(rank==0)logmsg(1,"%s solve time %.2f s\n",euler?"Euler":midpoint?"Midpoint":"Gragg",(endtime.tv_sec - begintime.tv_sec)+((double)(endtime.tv_usec - begintime.tv_usec))/ 1000000);
    teems_rss_probe("solve");
              free(counteqs);
              free(counteqnoadds);
              free(countvarintra1s);
              elem_vals1=NULL;
              if(x1!=NULL)free(x1);
    return true;
}

