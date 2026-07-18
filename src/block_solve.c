#include <teems_solver.h>
#include <hsl_kernels.h>

/* ==== NDBBD in-memory factor store (5.9) ====
   Per-block MA48 factor copies (VA/IRN/KEEP) held between the
   factorizations (regional blocks in ndbbd_presolve, interface blocks
   in ndbbd_solve's prep phase) and their consumers later in the same
   linear solve. Active only under -inmemory; disk mode keeps the
   legacy _vav/_irnv/_keep scratch files. Indexed by the rank-local
   block id (j4); puts from OpenMP loops touch distinct slots. */
static int ndbbd_fac_n=0;
static int **ndbbd_fac_irn=NULL,**ndbbd_fac_keep=NULL;
static solve_real **ndbbd_fac_va=NULL;
/* -fastrefac NDBBD extension: the regional (inner) blocks' factors
   persist across steps in this store — repeat steps refill VA(1:NE)
   and refactorize with MA48B/BD JOB=2 on the kept pivot sequence.
   ndbbd_fac_jcn holds the MA48 JCN mapping (needed for JOB=2) and
   ndbbd_fac_nz the staged pattern size; nz>=0 with a non-NULL jcn
   marks a persistent slot.  Per-TIME interface blocks keep the
   per-step flow (their pattern is assembled from value-dependent
   merges) and never set nz.  The flag forces the resident store
   regardless of -inmemory. */
static int **ndbbd_fac_jcn=NULL;
static PetscInt *ndbbd_fac_nz=NULL;
static dim_t ndbbd_fastrefac=-1;
static int nfr_flag(void) {
  if(ndbbd_fastrefac<0) {
    ndbbd_fastrefac=0;
    PetscOptionsGetInt(NULL,NULL,"-fastrefac",&ndbbd_fastrefac,NULL);
  }
  return (int)ndbbd_fastrefac;
}

/* ==== DBBD per-block persistent factors (-fastrefac) ====
   The diagonal-block COO is a raw copy of the stored CSR (no zero
   filter), so the pattern is structurally stable across steps: the
   MA48 analyse runs once per block per solve and later steps
   refactorize with MA48B/BD JOB=2 on the kept pivot sequence.
   IRN/JCN/VA/KEEP persist here (they double as the within-step
   factor handoff to the back-solve, replacing fac_* / scratch
   files).  Slots are indexed by the rank-local block id; the OMP
   factorization loop touches distinct slots. */
static int **dfr_irn=NULL,**dfr_jcn=NULL,**dfr_keep=NULL;
static solve_real **dfr_va=NULL;
static PetscInt *dfr_nz=NULL;
static int dfr_nblocks=0,dfr_ready=0;

/* -fastrefac DBBD extraction persistence: the five submatrix arrays
   and the per-block transposes survive across steps and are refilled
   with MAT_REUSE_MATRIX (the extraction structure is fixed while the
   orderings are; a row/col-ordering change triggers a full rebuild).
   The per-step IS/index construction stays as-is — REUSE only needs
   IS content, not object identity. */
static Mat *xfr_submatA=NULL,*xfr_submatC=NULL,*xfr_submatB=NULL,*xfr_submatBB=NULL,*xfr_submatD=NULL,*xfr_submatCT=NULL;
static int *xfr_order=NULL;
static PetscInt xfr_ordlen=0,xfr_nmatin=0;
static int xfr_proc1=0,xfr_ready=0;

void dbbd_fastextract_free(void) {
  PetscInt i;
  for(i=0; i<xfr_nmatin; i++) {
    if(xfr_submatA!=NULL)MatDestroy(&xfr_submatA[i]);
    if(xfr_submatC!=NULL)MatDestroy(&xfr_submatC[i]);
    if(xfr_submatB!=NULL)MatDestroy(&xfr_submatB[i]);
    if(xfr_submatCT!=NULL)MatDestroy(&xfr_submatCT[i]);
  }
  if(xfr_submatBB!=NULL)MatDestroy(&xfr_submatBB[0]);
  if(xfr_proc1&&xfr_submatD!=NULL)MatDestroy(&xfr_submatD[0]);
  if(xfr_submatA!=NULL)PetscFree(xfr_submatA);
  if(xfr_submatC!=NULL)PetscFree(xfr_submatC);
  if(xfr_submatB!=NULL)PetscFree(xfr_submatB);
  if(xfr_submatBB!=NULL)PetscFree(xfr_submatBB);
  if(xfr_submatD!=NULL)PetscFree(xfr_submatD);
  if(xfr_submatCT!=NULL)PetscFree(xfr_submatCT);
  xfr_submatA=NULL;
  xfr_submatC=NULL;
  xfr_submatB=NULL;
  xfr_submatBB=NULL;
  xfr_submatD=NULL;
  xfr_submatCT=NULL;
  free(xfr_order);
  xfr_order=NULL;
  xfr_ordlen=0;
  xfr_nmatin=0;
  xfr_proc1=0;
  xfr_ready=0;
}

void dbbd_fastrefac_free(void) {
  int i;
  for(i=0; i<dfr_nblocks; i++) {
    free(dfr_irn[i]);
    free(dfr_jcn[i]);
    free(dfr_keep[i]);
    free(dfr_va[i]);
  }
  free(dfr_irn);
  free(dfr_jcn);
  free(dfr_keep);
  free(dfr_va);
  free(dfr_nz);
  dfr_irn=NULL;
  dfr_jcn=NULL;
  dfr_keep=NULL;
  dfr_va=NULL;
  dfr_nz=NULL;
  dfr_nblocks=0;
  dfr_ready=0;
}

static void ndbbd_fac_init(int n) {
  int i;
  if(ndbbd_fac_n<n) {
    ndbbd_fac_irn=realloc(ndbbd_fac_irn,n*sizeof(int*));
    ndbbd_fac_keep=realloc(ndbbd_fac_keep,n*sizeof(int*));
    ndbbd_fac_va=realloc(ndbbd_fac_va,n*sizeof(solve_real*));
    ndbbd_fac_jcn=realloc(ndbbd_fac_jcn,n*sizeof(int*));
    ndbbd_fac_nz=realloc(ndbbd_fac_nz,n*sizeof(PetscInt));
    for(i=ndbbd_fac_n; i<n; i++) {
      ndbbd_fac_irn[i]=NULL;
      ndbbd_fac_keep[i]=NULL;
      ndbbd_fac_va[i]=NULL;
      ndbbd_fac_jcn[i]=NULL;
      ndbbd_fac_nz[i]=-1;
    }
    ndbbd_fac_n=n;
  }
}

static void ndbbd_fac_put(int idx,int *irn,int *keep,solve_real *va,long la,int t) {
  ndbbd_fac_irn[idx]=realloc(ndbbd_fac_irn[idx],la*sizeof(int));
  memcpy(ndbbd_fac_irn[idx],irn,la*sizeof(int));
  ndbbd_fac_keep[idx]=realloc(ndbbd_fac_keep[idx],t*sizeof(int));
  memcpy(ndbbd_fac_keep[idx],keep,t*sizeof(int));
  ndbbd_fac_va[idx]=realloc(ndbbd_fac_va[idx],la*sizeof(solve_real));
  memcpy(ndbbd_fac_va[idx],va,la*sizeof(solve_real));
}

static void ndbbd_fac_drop(int idx) {
  free(ndbbd_fac_irn[idx]);
  ndbbd_fac_irn[idx]=NULL;
  free(ndbbd_fac_keep[idx]);
  ndbbd_fac_keep[idx]=NULL;
  free(ndbbd_fac_va[idx]);
  ndbbd_fac_va[idx]=NULL;
  if(ndbbd_fac_jcn!=NULL) {
    free(ndbbd_fac_jcn[idx]);
    ndbbd_fac_jcn[idx]=NULL;
    ndbbd_fac_nz[idx]=-1;
  }
}

void ndbbd_fastrefac_free(void) {
  int i;
  for(i=0; i<ndbbd_fac_n; i++)ndbbd_fac_drop(i);
  free(ndbbd_fac_irn);
  free(ndbbd_fac_keep);
  free(ndbbd_fac_va);
  free(ndbbd_fac_jcn);
  free(ndbbd_fac_nz);
  ndbbd_fac_irn=NULL;
  ndbbd_fac_keep=NULL;
  ndbbd_fac_va=NULL;
  ndbbd_fac_jcn=NULL;
  ndbbd_fac_nz=NULL;
  ndbbd_fac_n=0;
}

/* hand one block's factors off: resident copy under -inmemory,
   otherwise the legacy scratch files the Fortran writers used to
   produce (same names, same bytes) */
static void ndbbd_fac_emit(int rank,int idx,int *irn,int *keep,solve_real *va,long la,int t) {
  if(inmemory||nfr_flag()) {
    ndbbd_fac_put(idx,irn,keep,va,la,t);
    return;
  }
  char fn[NAMESIZE+32];
  FILE *fp;
  sprintf(fn,"%s_vav%04d%04d.bin",scratch_dir,rank,idx);
  fp=fopen(fn,"wb");
  fwrite(va,sizeof(solve_real),la,fp);
  fclose(fp);
  sprintf(fn,"%s_irnv%04d%04d.bin",scratch_dir,rank,idx);
  fp=fopen(fn,"wb");
  fwrite(irn,sizeof(int),la,fp);
  fclose(fp);
  sprintf(fn,"%s_keep%04d%04d.bin",scratch_dir,rank,idx);
  fp=fopen(fn,"wb");
  fwrite(keep,sizeof(int),t,fp);
  fclose(fp);
}

int dbbd_solve(Mat A, Vec b, solve_real *x1, offset_t VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int *row_order,int *col_order, offset_t ndblock,int *block_sizes, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,dim_t laD,PetscReal cntl3) {//,bool iter
  IS *rowindices,*colindices,*Cindices,*Bindices,*BBindices,*BBiindices;
  const PetscInt *nindices;
  PetscInt bfirst,bend,sblockin,nmatin,nmatinplus,nrowcolin,sumrowcolin,i,i1,j,j0,j1,j2,j3,j4,j5,j6,l0,l1,l2,l3,l4,l5,rank1,proc1=0,nnzmax,j1nz,j1irnbs;
  Mat *submatA,*submatC,*submatB,*submatD,submatBT,submatCT,*submatBB;
  Vec vecxd;
  PetscInt *ai,*aj,*aic,*ajc;
  PetscInt nrow,ncol,nz,nrowc,ncolc,nrowb,ncolb,nzc,nzv;
  offset_t lasize,ldsize;
  MPI_Status   status;
  clock_t timestr,timeend;
  int j2int,la1;
  size_t freadresult;
  PetscScalar *vals,*valsc,vecval;
  solve_real vval;
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
  nmatin=(offset_t)ndblock/mpisize;
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
  /* -fastrefac: per-block persistent factors + persistent extraction;
     option read once */
  static dim_t dbbd_fastrefac=-1;
  if(dbbd_fastrefac<0) {
    dbbd_fastrefac=0;
    PetscOptionsGetInt(NULL,NULL,"-fastrefac",&dbbd_fastrefac,NULL);
  }
  if(dbbd_fastrefac&&dfr_nblocks!=nmatin) {
    dbbd_fastrefac_free();
    dfr_irn=(int**)calloc(nmatin,sizeof(int*));
    dfr_jcn=(int**)calloc(nmatin,sizeof(int*));
    dfr_keep=(int**)calloc(nmatin,sizeof(int*));
    dfr_va=(solve_real**)calloc(nmatin,sizeof(solve_real*));
    dfr_nz=(PetscInt*)calloc(nmatin,sizeof(PetscInt));
    for(i=0; i<nmatin; i++)dfr_nz[i]=-1;
    dfr_nblocks=nmatin;
  }
  /* the persisted extraction is valid only while the orderings are:
     any row/col-ordering change forces a fresh MAT_INITIAL build */
  int xfr_reuse=0;
  if(dbbd_fastrefac) {
    if(xfr_ready&&xfr_ordlen==2*(PetscInt)VecSize
       &&!memcmp(xfr_order,row_order,VecSize*sizeof(int))
       &&!memcmp(xfr_order+VecSize,col_order,VecSize*sizeof(int)))xfr_reuse=1;
    else {
      dbbd_fastextract_free();
      xfr_order=realloc(xfr_order,2*VecSize*sizeof(int));
      memcpy(xfr_order,row_order,VecSize*sizeof(int));
      memcpy(xfr_order+VecSize,col_order,VecSize*sizeof(int));
      xfr_ordlen=2*(PetscInt)VecSize;
    }
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
  /* MatCreateSubMatrices(MAT_INITIAL_MATRIX) allocates the submatrix
     arrays itself and overwrites the passed pointer — pre-allocating
     them here leaked one small array per call site per step */
  submatA=NULL;
  submatC=NULL;
  submatB=NULL;
  submatBB=NULL;
  PetscScalar **yi1= (PetscScalar**)calloc(nmatin,sizeof(PetscScalar*));
  for (i=0; i<nmatin; i++) yi1[i] = (PetscScalar*)calloc(block_sizes[i+begblock[rank]],sizeof(PetscScalar));
  sumrowcolin=0;
  for(j=0; j<ndblock; j++)sumrowcolin+=block_sizes[j];
  PetscInt *indicesC= (PetscInt *) calloc (VecSize-sumrowcolin,sizeof(PetscInt));
  PetscInt *indicesB= (PetscInt *) calloc (VecSize-sumrowcolin,sizeof(PetscInt));
  PetscInt *offblock= (PetscInt *) calloc (ndblock+1,sizeof(PetscInt));
  PetscInt *offblockrow= (PetscInt *) calloc (ndblock+1,sizeof(PetscInt));
  solve_real *vecbiui= (solve_real *) calloc (VecSize-sumrowcolin,sizeof(solve_real));
  offblock[0]=0;
  for(j=0; j<ndblock; j++)offblock[j+1]=offblock[j]+countvarintra1[j+1]-countvarintra1[j]-block_sizes[j];
  offblockrow[0]=0;
  for(j=0; j<ndblock; j++)offblockrow[j+1]=offblockrow[j]+counteqnoadd[j]-block_sizes[j];

  /* factor handoff between the block-factorization and back-solve
     phases: resident buffers under -inmemory, legacy _irnv/_keep/_vav
     scratch files otherwise (same bytes the Fortran writers produced) */
  int **fac_irn=NULL,**fac_keep=NULL;
  solve_real **fac_va=NULL;
  if(inmemory) {
    fac_irn=(int**)calloc(nmatinplus,sizeof(int*));
    fac_keep=(int**)calloc(nmatinplus,sizeof(int*));
    fac_va=(solve_real**)calloc(nmatinplus,sizeof(solve_real*));
  }
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
      strcpy(filename,scratch_dir);strcat(filename,"_irnv");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      strcpy(fn01[j1],filename);
      strcpy(filename,scratch_dir);strcat(filename,"_keep");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      strcpy(fn02[j1],filename);
      strcpy(filename,scratch_dir);strcat(filename,"_vav");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      strcpy(fn03[j1],filename);
    }
  }
  for(j=0; j<mpisize; j++) {
    for(j3=0; j3<nmatinplus; j3++) {
      if(j3==nmatin)i=j3-1;
      else i=j3;
      bfirst=counteq[i+begblock[rank]];
      bend=block_sizes[i+begblock[rank]];
      MPI_Bcast(&bfirst,1, MPI_INT,j, PETSC_COMM_WORLD);
      MPI_Bcast(&bend,1, MPI_INT,j, PETSC_COMM_WORLD);
      PetscInt *indices= (PetscInt *) calloc (bend,sizeof(PetscInt));
      solve_real *yi0 = (solve_real*)calloc(bend,sizeof(solve_real));
      for(j1=0; j1<bend; j1++) { //ha_ndblocks[i+begblock[j]]
        j2=row_order[bfirst+j1]+bfirst;
        if(j2>=Istart&&j2<Iend)indices[j1]=j2;
        else indices[j1]=-1;
      }
      VecGetValues(b,bend,indices,yi0);

      reduce_to_rank(yi0,(fortran_int)bend,mpisize,rank,j);
      if(rank==j&&i<nmatin) {
        memcpy (yi1[i],yi0,bend*sizeof(solve_real));
      }
      free(indices);
      free(yi0);
    }
  }
  j1=0;
  for(i=0; i<nmatin; i++) {
    bfirst=counteq[i+begblock[rank]];
    PetscInt *indices= (PetscInt *) calloc (block_sizes[i+begblock[rank]],sizeof(PetscInt));
    for(j=0; j<block_sizes[i+begblock[rank]]; j++) {
      indices[j]=row_order[bfirst+j]+bfirst;
    }

    bend=counteq[i+begblock[rank]]+counteqnoadd[i+begblock[rank]];
    for(j=offblockrow[i+begblock[rank]]; j<offblockrow[i+1+begblock[rank]]; j++) {
      j1=j-offblockrow[i+begblock[rank]];
      indicesB[j]=row_order[bend-j1-1]+bfirst;
    }
    ISCreateGeneral(PETSC_COMM_SELF,block_sizes[i+begblock[rank]],indices,PETSC_COPY_VALUES,rowindices+i);
    bfirst=countvarintra1[i+begblock[rank]];
    for(j=0; j<block_sizes[i+begblock[rank]]; j++) {
      indices[j]=col_order[bfirst+j]+bfirst;
    }
    bend=countvarintra1[i+1+begblock[rank]];
    for(j=offblock[i+begblock[rank]]; j<offblock[i+1+begblock[rank]]; j++) {
      j1=j-offblock[i+begblock[rank]];
      indicesC[j]=col_order[bend-j1-1]+bfirst;
    }

    ISCreateGeneral(PETSC_COMM_SELF,block_sizes[i+begblock[rank]],indices,PETSC_COPY_VALUES,colindices+i);
    free(indices);
  }
  j1=offblock[ndblock];
  for(i=countvarintra1[ndblock]; i<VecSize; i++) {
    indicesC[j1]=i;
    j1++;
  }
  j1=offblockrow[ndblock];
  for(i=counteq[ndblock-1]+counteqnoadd[ndblock-1]; i<VecSize; i++) {
    indicesB[j1]=i;
    j1++;
  }
  for(i=0; i<mpisize; i++) {
    j=offblockrow[begblock[rank]];
    j1=offblockrow[nmatin+begblock[rank]]-offblockrow[begblock[rank]];
    MPI_Bcast(&j,1, MPI_INT,i, PETSC_COMM_WORLD);
    MPI_Bcast(&j1,1, MPI_INT,i, PETSC_COMM_WORLD);
    MPI_Bcast(&indicesB[j],j1, MPI_INT,i, PETSC_COMM_WORLD);
  }
  /* this rank's [Istart,Iend) slice of the border rows; a rank owning
     no border rows must yield an empty range, not stale j/j1 */
  j=VecSize-sumrowcolin;
  for(i=0; i<VecSize-sumrowcolin; i++)if(indicesB[i]>=Istart) {
      j=i;
      break;
    }
  j1=VecSize-sumrowcolin;
  for(j2=j; j2<VecSize-sumrowcolin; j2++)if(indicesB[j2]>=Iend) {
      j1=j2;
      break;
    }
  logmsg(2,"rank %d j %d j1 %d istart %d iend %d\n",rank,j,j1,Istart,Iend);
  VecGetValues(b,j1-j,&indicesB[j],&vecbiui[j]);//implicite yd save mem
  ierr = VecDestroy(&b);
  CHKERRQ(ierr);
  for(i=0; i<mpisize; i++) {
    j=offblock[begblock[rank]];
    j1=offblock[nmatin+begblock[rank]]-offblock[begblock[rank]];
    MPI_Bcast(&j,1, MPI_INT,i, PETSC_COMM_WORLD);
    MPI_Bcast(&j1,1, MPI_INT,i, PETSC_COMM_WORLD);
    MPI_Bcast(&indicesC[j],j1, MPI_INT,i, PETSC_COMM_WORLD);
  }
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
  logmsg(2,"Begin partitioning rank %d\n",rank);
  if(dbbd_fastrefac&&xfr_reuse) {
    ierr = MatCreateSubMatrices(A,nmatin,rowindices,colindices,MAT_REUSE_MATRIX,&xfr_submatA);
    ierr = MatCreateSubMatrices(A,nmatin,rowindices,Cindices,MAT_REUSE_MATRIX,&xfr_submatC);
    ierr = MatCreateSubMatrices(A,nmatinBB,Bindices,BBindices,MAT_REUSE_MATRIX,&xfr_submatBB);
  }
  else {
    ierr = MatCreateSubMatrices(A,nmatin,rowindices,colindices,MAT_INITIAL_MATRIX,&submatA);
    ierr = MatCreateSubMatrices(A,nmatin,rowindices,Cindices,MAT_INITIAL_MATRIX,&submatC);
    ierr = MatCreateSubMatrices(A,nmatinBB,Bindices,BBindices,MAT_INITIAL_MATRIX,&submatBB);
    if(dbbd_fastrefac) {
      xfr_submatA=submatA;
      xfr_submatC=submatC;
      xfr_submatBB=submatBB;
    }
  }
  if(dbbd_fastrefac) {
    submatA=xfr_submatA;
    submatC=xfr_submatC;
    submatBB=xfr_submatBB;
  }
  ierr = MatDestroy(&A);
  CHKERRQ(ierr);
  PetscInt BBrow,BBcol,CCrow,CCcol;
  MatGetSize(submatBB[0],&BBrow,&BBcol);
  MatGetSize(submatC[0],&CCrow,&CCcol);
  logmsg(2,"bbrow %d bbcol %d ccrow %d cccol %d\n",BBrow,BBcol,BBcol,CCcol);
  PetscInt *indexBBi= (PetscInt *) calloc (BBrow,sizeof(PetscInt));
  for(i=0; i<BBrow; i++)indexBBi[i]=i;
  for(i=0; i<nmatin; i++) {
    /* re-created with new handles below; the first-generation IS was
       leaking here every step (only the second generation reached the
       destroy loop) */
    ierr = ISDestroy(&Bindices[i]);
    CHKERRQ(ierr);
    ISCreateGeneral(PETSC_COMM_SELF,BBrow,indexBBi,PETSC_COPY_VALUES,Bindices+i);
  }
  if(rank==mpisize-1)proc1=1;
  if(dbbd_fastrefac&&xfr_reuse) {
    ierr = MatCreateSubMatrices(submatBB[0],nmatin,Bindices,colindices,MAT_REUSE_MATRIX,&xfr_submatB);
    MatCreateSubMatrices(submatBB[0],proc1,Bindices,Cindices,MAT_REUSE_MATRIX,&xfr_submatD);
  }
  else {
    ierr = MatCreateSubMatrices(submatBB[0],nmatin,Bindices,colindices,MAT_INITIAL_MATRIX,&submatB);
    submatD=NULL;
    MatCreateSubMatrices(submatBB[0],proc1,Bindices,Cindices,MAT_INITIAL_MATRIX,&submatD);
    if(dbbd_fastrefac) {
      xfr_submatB=submatB;
      xfr_submatD=submatD;
      xfr_nmatin=nmatin;
      xfr_proc1=proc1;
    }
  }
  if(dbbd_fastrefac) {
    submatB=xfr_submatB;
    submatD=xfr_submatD;
  }
  logmsg(2,"End partitioning rank %d time %f\n",rank,((double)clock()-timestr)/CLOCKS_PER_SEC);
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
  if(!dbbd_fastrefac) {
  ierr = MatDestroy(&submatBB[0]);
  PetscFree(submatBB);
  }
  CHKERRQ(ierr);
  solve_real *xi1 = (solve_real*)calloc(sumrowcolin,sizeof(solve_real));
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
  if(vecbivisize<0)printf("Error: interface problem size overflow (negative vecbivisize); the border is too large\n");
  solve_real *vecbivi= (solve_real *) calloc (vecbivisize,sizeof(solve_real));
  logmsg(2,"vecbivisize %ld rank %d\n",vecbivisize,rank);

  solve_real *xi1point;
  offset_t xi1indx=0;
  int jthrd,nthrd=1;
  nthrd=omp_get_max_threads();
  offset_t *nthrds= (offset_t *) calloc (nthrd+1,sizeof(offset_t));
  offset_t *nthrds1= (offset_t *) calloc (nthrd+1,sizeof(offset_t));
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
      nthrds[jthrd+1]+=block_sizes[j1+begblock[rank]];
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
      int *irn=NULL,*jcn=NULL;
      solve_real *values=NULL;
      int dfr_redo=0;
      if(dbbd_fastrefac) {
        if(dfr_ready&&dfr_nz[j1]==nz) {
          /* pattern unchanged: refill values, keep the pivot sequence */
          memcpy (dfr_va[j1],vals,nz*sizeof(solve_real));
          dfr_redo=1;
        }
        else {
          free(dfr_irn[j1]);
          free(dfr_jcn[j1]);
          free(dfr_keep[j1]);
          free(dfr_va[j1]);
          dfr_nz[j1]=nz;
          dfr_irn[j1]=(int *) calloc (lasize,sizeof(int));
          dfr_jcn[j1]=(int *) calloc (lasize,sizeof(int));
          dfr_va[j1]=(solve_real *) calloc (lasize,sizeof(solve_real));
          dfr_keep[j1]=(int *) calloc (nrow+9*ncol+7,sizeof(int));
          memcpy (dfr_irn[j1]+nz,ai,(nrow+1)*sizeof(PetscInt));
          memcpy (dfr_jcn[j1],aj,nz*sizeof(PetscInt));
          memcpy (dfr_va[j1],vals,nz*sizeof(solve_real));
        }
      }
      else {
        irn=(int *) calloc (lasize,sizeof(int));
        jcn=(int *) calloc (lasize,sizeof(int));
        values= (solve_real *) calloc (lasize,sizeof(solve_real));
      }
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
      if(!dbbd_fastrefac) {
        memcpy (irn+nz,ai,(nrow+1)*sizeof(PetscInt));
        memcpy (jcn,aj,nz*sizeof(PetscInt));
        memcpy (values,vals,nz*sizeof(solve_real));
        MatDestroy(&submatA[j1]);
      }
      insize[j1*insizes+2]=nz;
      insize[j1*insizes+5]=nzc;
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
      int *keep=NULL;
      if(dbbd_fastrefac) {
        int redo_io=dfr_redo;
        spec48m_msol_p_(insize+j1*insizes,dfr_irn[j1],dfr_jcn[j1],dfr_va[j1],yi1[j1],xi1point,aic,ajc,valsc,ai,aj,vals,vecbivi,bivinzrow,bivinzcol,dfr_keep[j1],&redo_io);
        if(redo_io<0) {
          /* fast refactorize declined: re-stage the block from the
             still-live submatrix and redo the analyse */
          Mat_SeqAIJ *aa2=(Mat_SeqAIJ*)submatA[j1]->data;
          memcpy (dfr_irn[j1]+dfr_nz[j1],aa2->i,(nrow+1)*sizeof(PetscInt));
          memcpy (dfr_jcn[j1],aa2->j,dfr_nz[j1]*sizeof(PetscInt));
          memcpy (dfr_va[j1],aa2->a,dfr_nz[j1]*sizeof(solve_real));
          redo_io=0;
          spec48m_msol_p_(insize+j1*insizes,dfr_irn[j1],dfr_jcn[j1],dfr_va[j1],yi1[j1],xi1point,aic,ajc,valsc,ai,aj,vals,vecbivi,bivinzrow,bivinzcol,dfr_keep[j1],&redo_io);
        }
      }
      else {
      keep=(int *) calloc (nrow+9*ncol+7,sizeof(int));/* KEEP bound, ICNTL(6)=1; live length returned in insize[12] */
      spec48m_msol_(insize+j1*insizes,irn,jcn,values,yi1[j1],xi1point,aic,ajc,valsc,ai,aj,vals,vecbivi,bivinzrow,bivinzcol,keep);
      }
      MatDestroy(&submatCT);
      free(yi1[j1]);
      if(!dbbd_fastrefac) {
      free(jcn);
      if(inmemory) {
        fac_irn[j1]=irn;
        fac_keep[j1]=keep;
        fac_va[j1]=values;
      }
      else {
        FILE *ffp;
        ffp=fopen(fn03[j1],"wb");
        fwrite(values,sizeof(solve_real),lasize,ffp);
        fclose(ffp);
        ffp=fopen(fn01[j1],"wb");
        fwrite(irn,sizeof(int),lasize,ffp);
        fclose(ffp);
        ffp=fopen(fn02[j1],"wb");
        fwrite(keep,sizeof(int),insize[j1*insizes+12],ffp);
        fclose(ffp);
        free(irn);
        free(keep);
        free(values);
      }
      }
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
      if(!dbbd_fastrefac)MatDestroy(&submatB[j1]);//submatBT);
      xi1indx+=block_sizes[j1+begblock[rank]];
      time(&timeend);
      logmsg(2,"Submatrix %d rank %d thrd %d calculation time %f\n",j1,rank,jthrd,difftime(timeend,timestr));
    }
  }
  }
  if(dbbd_fastrefac) {
    dfr_ready=1;
    xfr_ready=1;
  }
  free(nthrds);
  free(nthrds1);
  free(bivinzcol);
  free(bivinzrow);
  if(!dbbd_fastrefac) {
  ierr = PetscFree(submatA);
  ierr = PetscFree(submatB);
  }
  CHKERRQ(ierr);
  free(yi1);
  logmsg(2,"Completed calculation of partitioned matrices, rank %d\n",rank);
  timestr=clock();

  #pragma omp parallel private(li) reduction(+:lj)
  {
  lj=0;
  #pragma omp for
  for(li=0; li<vecbivisize; ++li){
    if(vecbivi[li]!=0)lj++;
  }
  }
  /* first touch of these is realloc: they must start NULL (the NDBBD
     copy of this block already had the fix; an indeterminate pointer
     here is a latent glibc invalid-pointer abort) */
  long int *obiviindx1=NULL,*biviindx1=NULL,*obiviindx0=NULL,*biviindx0=NULL;
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
  vecbivi=realloc(vecbivi,lj2*sizeof(solve_real));
  nz0=lj2;
  logmsg(2,"lj %ld lj2 %ld rank %d nz0 %ld nz3 %ld\n",lj,lj2,rank,nz0,nz3);
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len,name_len_max,name_beg,class_size,color,group_size,node_rank,group_size1,group_size11,node_tail_rank;
  MPI_Comm_rank( node_comm, &node_rank);
  MPI_Comm_size(node_comm,&group_size);
  if(node_rank==group_size-1)color=1;
  else color=0;
  MPI_Comm_rank( node_tail_comm, &node_tail_rank);
  MPI_Comm_size(node_tail_comm,&group_size1);
  group_size11=group_size1;
  MPI_Bcast(&group_size11,1, MPI_INT,mpisize-1, PETSC_COMM_WORLD);
  int sindx01,sindx02,sindx03;

  if(node_rank!=group_size-1) {
    lj2=nz3;
    MPI_Send(&lj2,1, MPI_LONG,group_size-1, 10, node_comm);
    sindx01=lj2/MAXSSIZE;
    sindx03=0;
    for(sindx02=0; sindx02<sindx01; sindx02++) {
      if(SORD==1)MPI_Send(vecbivi+sindx03,MAXSSIZE, MPI_DOUBLE,group_size-1, 11*sindx01+sindx02, node_comm);
      else MPI_Send(vecbivi+sindx03,MAXSSIZE, MPI_FLOAT,group_size-1, 11*sindx01+sindx02, node_comm);
      sindx03+=MAXSSIZE;
    }
    if(SORD==1)MPI_Send(vecbivi+sindx03,lj2-sindx01*MAXSSIZE, MPI_DOUBLE,group_size-1, 11, node_comm);
    else MPI_Send(vecbivi+sindx03,lj2-sindx01*MAXSSIZE, MPI_FLOAT,group_size-1, 11, node_comm);
    free(vecbivi);
    vecbivi=NULL;
    sindx03=0;
    for(sindx02=0; sindx02<sindx01; sindx02++) {
      MPI_Send(biviindx1+sindx03,MAXSSIZE, MPI_LONG,group_size-1, 12*sindx01+sindx02, node_comm);
      sindx03+=MAXSSIZE;
    }
    MPI_Send(biviindx1+sindx03,lj2-sindx01*MAXSSIZE, MPI_LONG,group_size-1, 12, node_comm);
    free(obiviindx1);
    obiviindx1=NULL;
  }
  if(node_rank==group_size-1) {
    for(j6=0; j6<group_size-1; j6++) {
      MPI_Recv(&lj2,1, MPI_LONG,  MPI_ANY_SOURCE,10, node_comm,&status);
      i=status.MPI_SOURCE;
      obiviindx0=realloc(obiviindx0,(lj2+1)*sizeof(long int));
      obiviindx0[0]=-1;
      biviindx0=obiviindx0;
      biviindx0++;
      solve_real *vecbivi0= (solve_real *) calloc (lj2,sizeof(solve_real));
      sindx01=lj2/MAXSSIZE;
      sindx03=0;
      for(sindx02=0; sindx02<sindx01; sindx02++) {
        if(SORD==1)MPI_Recv(vecbivi0+sindx03,MAXSSIZE, MPI_DOUBLE,i,11*sindx01+sindx02, node_comm,&status);
        else MPI_Recv(vecbivi0+sindx03,MAXSSIZE, MPI_FLOAT,i,11*sindx01+sindx02, node_comm,&status);
        sindx03+=MAXSSIZE;
      }
      if(SORD==1)MPI_Recv(vecbivi0+sindx03,lj2-sindx01*MAXSSIZE, MPI_DOUBLE,i,11, node_comm,&status);
      else MPI_Recv(vecbivi0+sindx03,lj2-sindx01*MAXSSIZE, MPI_FLOAT,i,11, node_comm,&status);
      sindx03=0;
      for(sindx02=0; sindx02<sindx01; sindx02++) {
        MPI_Recv(biviindx0+sindx03,MAXSSIZE, MPI_LONG,i,12*sindx01+sindx02, node_comm,&status);
        sindx03+=MAXSSIZE;
      }
      MPI_Recv(biviindx0+sindx03,lj2-sindx01*MAXSSIZE, MPI_LONG,i,12, node_comm,&status);
      my_spar_compl_(biviindx1,&nz0,biviindx0,&lj2,&nz1);
      vecbivi=realloc(vecbivi,nz1*sizeof(solve_real));
      obiviindx1=realloc(obiviindx1,(nz1+1)*sizeof(long int));
      obiviindx1[0]=-1;
      biviindx1=obiviindx1;
      biviindx1++;
      my_spar_add3l_(vecbivi,obiviindx1,&nz0,vecbivi0,obiviindx0,&lj2,&nz1);
      free(obiviindx0);
      obiviindx0=NULL;
      free(vecbivi0);
      vecbivi0=NULL;
      logmsg(2,"rank %d nz0 %ld j2 %ld nz1 %ld\n",rank,nz0,lj2,nz1);
      nz0=nz1;
    }
  }

  if(group_size11>1&&color==1) {
    nz3=nz0;
    if(node_tail_rank!=group_size1-1) {
      lj2=nz3;
      MPI_Send(&lj2,1, MPI_LONG,group_size1-1, 10, node_tail_comm);
      sindx01=lj2/MAXSSIZE;
      sindx03=0;
      for(sindx02=0; sindx02<sindx01; sindx02++) {
        if(SORD==1) MPI_Send(vecbivi+sindx03,MAXSSIZE, MPI_DOUBLE,group_size1-1, 11*sindx01+sindx02, node_tail_comm);
        else MPI_Send(vecbivi+sindx03,MAXSSIZE, MPI_FLOAT,group_size1-1, 11*sindx01+sindx02, node_tail_comm);
        sindx03+=MAXSSIZE;
      }
      if(SORD==1) MPI_Send(vecbivi+sindx03,lj2-sindx01*MAXSSIZE, MPI_DOUBLE,group_size1-1, 11, node_tail_comm);
      else MPI_Send(vecbivi+sindx03,lj2-sindx01*MAXSSIZE, MPI_FLOAT,group_size1-1, 11, node_tail_comm);
      sindx03=0;
      for(sindx02=0; sindx02<sindx01; sindx02++) {
        MPI_Send(biviindx1+sindx03,MAXSSIZE, MPI_LONG,group_size1-1, 12*sindx01+sindx02, node_tail_comm);
        sindx03+=MAXSSIZE;
      }
      MPI_Send(biviindx1+sindx03,lj2-sindx01*MAXSSIZE, MPI_LONG,group_size1-1, 12, node_tail_comm);
      free(obiviindx1);
      obiviindx1=NULL;
      free(vecbivi);
      vecbivi=NULL;
    }
    if(node_tail_rank==group_size1-1) {
      for(j6=0; j6<group_size1-1; j6++) {
        MPI_Recv(&lj2,1, MPI_LONG,  MPI_ANY_SOURCE,10, node_tail_comm,&status);
        i=status.MPI_SOURCE;
        obiviindx0=realloc(obiviindx0,(lj2+1)*sizeof(long int));
        obiviindx0[0]=-1;
        biviindx0=obiviindx0;
        biviindx0++;
        solve_real *vecbivi0= (solve_real *) calloc (lj2,sizeof(solve_real));
        sindx01=lj2/MAXSSIZE;
        sindx03=0;
        for(sindx02=0; sindx02<sindx01; sindx02++) {
          if(SORD==1) MPI_Recv(vecbivi0+sindx03,MAXSSIZE, MPI_DOUBLE,i,11*sindx01+sindx02, node_tail_comm,&status);
          else MPI_Recv(vecbivi0+sindx03,MAXSSIZE, MPI_FLOAT,i,11*sindx01+sindx02, node_tail_comm,&status);
          sindx03+=MAXSSIZE;
        }
        if(SORD==1) MPI_Recv(vecbivi0+sindx03,lj2-sindx01*MAXSSIZE, MPI_DOUBLE,i,11, node_tail_comm,&status);
        else MPI_Recv(vecbivi0+sindx03,lj2-sindx01*MAXSSIZE, MPI_FLOAT,i,11, node_tail_comm,&status);
        sindx03=0;
        for(sindx02=0; sindx02<sindx01; sindx02++) {
          MPI_Recv(biviindx0+sindx03,MAXSSIZE, MPI_LONG,i,12*sindx01+sindx02, node_tail_comm,&status);
          sindx03+=MAXSSIZE;
        }
        MPI_Recv(biviindx0+sindx03,lj2-sindx01*MAXSSIZE, MPI_LONG,i,12, node_tail_comm,&status);
        my_spar_compl_(biviindx1,&nz0,biviindx0,&lj2,&nz1);
        vecbivi=realloc(vecbivi,nz1*sizeof(solve_real));
        obiviindx1=realloc(obiviindx1,(nz1+1)*sizeof(long int));
        obiviindx1[0]=-1;
        biviindx1=obiviindx1;
        biviindx1++;
        my_spar_add3l_(vecbivi,obiviindx1,&nz0,vecbivi0,obiviindx0,&lj2,&nz1);
        free(obiviindx0);
        obiviindx0=NULL;
        free(vecbivi0);
        logmsg(2,"rank %d nz0 %ld j2 %ld nz1 %ld\n",rank,nz0,lj2,nz1);
        nz0=nz1;
      }
    }
  }

  fortran_int vecbiuisize=0;
  vecbiuisize=VecSize-sumrowcolin;
  reduce_to_rank_nocompress(vecbiui,vecbiuisize,mpisize,rank,mpisize-1);
  logmsg(2,"Completed interface reduce, rank %d time %f\n",rank,((double)clock()-timestr)/CLOCKS_PER_SEC);
  if(rank!=mpisize-1){
    free(vecbiui);
    vecbiui=NULL;
  }
  timestr=clock();
  solve_real *xd;//= (ha_cgetype*)calloc(vecbiuisize,sizeof(ha_cgetype));
  solve_real *x0;//,*vecbivi0;//= (ha_cgetype*)calloc(VecSize,sizeof(ha_cgetype));
  long int j7,lnz;
  if(rank==mpisize-1) {
    Mat_SeqAIJ         *aa=(Mat_SeqAIJ*)submatD[0]->data;//*aa=subA->data;
    ai= aa->i;
    aj= aa->j;
    vals=aa->a;
    lnz=aa->nz;
    nrow=submatD[0]->rmap->n;
    ncol=submatD[0]->cmap->n;
    logmsg(2,"D ncol %d nrow %d nz %ld\n",ncol,nrow,lnz);
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
    logmsg(2,"nz0 %ld nz %ld nz1 %ld\n",nz0,lnz,nz1);



    logmsg(2,"nz1 %ld\n",nz1);
    ldsize=ceil((laD/100.0)*nz1);
    vecbivi=realloc(vecbivi,ldsize*sizeof(solve_real));
    int *irn1=(int *) calloc (nz1,sizeof(int));
    int *jcn=(int *) calloc (ldsize,sizeof(int));
    logmsg(2,"indx0 %ld indx1 %ld\n",obiviindx1[0],obiviindx0[0]);
    my_spar_add4l_(vecbivi,obiviindx1,irn1,jcn,&nz0,vals,obiviindx0,&lnz,&nz1,&ncol);//my_spar_add1l_



    if(!dbbd_fastrefac)MatDestroy(&submatD[0]);
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
    xd=(solve_real *) calloc (vecbiuisize,sizeof(solve_real));//realloc (xd,vecbiuisize*sizeof(ha_cgetype));
    spec48_ssol2la_(insizeD,irn1,jcn,vecbivi,vecbiui,xd);
    free(vecbivi);
    vecbivi=NULL;
    free(vecbiui);
    vecbiui=NULL;
    free(irn1);
    free(jcn);
    free(insizeD);
    x0=(solve_real *) calloc (VecSize,sizeof(solve_real));//realloc (x0,*sizeof(ha_cgetype));
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
    xd=(solve_real *) calloc (vecbiuisize,sizeof(solve_real));//realloc (xd,vecbiuisize*sizeof(ha_cgetype));
    x0=(solve_real *) calloc (VecSize,sizeof(solve_real));//realloc (x0,*sizeof(ha_cgetype));
  }
  if(!dbbd_fastrefac)PetscFree(submatD);
  MPI_Barrier(PETSC_COMM_WORLD);
  if(SORD==1)MPI_Bcast(xd, vecbiuisize, MPI_DOUBLE,mpisize-1, PETSC_COMM_WORLD);
  else MPI_Bcast(xd, vecbiuisize, MPI_FLOAT,mpisize-1, PETSC_COMM_WORLD);
  logmsg(2,"Interface problem solution rank %d time %f\n",rank,((double)clock()-timestr)/CLOCKS_PER_SEC);
  timestr=clock();
  xi1indx=0;
  for(j1=0; j1<nmatinplus; j1++) {
    if(j1<nmatin) {
      solve_real *biui0= (solve_real *) calloc (block_sizes[j1+begblock[rank]],sizeof(solve_real));
      int *irne,*keep;
      solve_real *vale;
      la1=ceil((insize[j1*insizes+9]/100.0)*insize[j1*insizes+2]);
      if(dbbd_fastrefac) {
        irne=dfr_irn[j1];
        keep=dfr_keep[j1];
        vale=dfr_va[j1];
      }
      else if(inmemory) {
        irne=fac_irn[j1];
        keep=fac_keep[j1];
        vale=fac_va[j1];
      }
      else {
        fp1 = fopen(fn01[j1], "rb");
        if (fp1==NULL)printf("Error: cannot open factor file %s\n",fn01[j1]);
        fp2 = fopen(fn02[j1], "rb");
        if (fp2==NULL)printf("Error: cannot open factor file %s\n",fn02[j1]);
        fp3 = fopen(fn03[j1], "rb");
        if (fp3==NULL)printf("Error: cannot open factor file %s\n",fn03[j1]);
        irne = (int*)calloc(la1,sizeof(int));
        keep = (int*)calloc(insize[j1*insizes+12],sizeof(int));
        vale = (solve_real*)calloc(la1,sizeof(solve_real));
        freadresult=fread(irne,sizeof(int),la1,fp1);
        freadresult=fread(keep,sizeof(int),insize[j1*insizes+12],fp2);
        freadresult=fread(vale,sizeof(solve_real),la1,fp3);
        fclose(fp1);
        fclose(fp2);
        fclose(fp3);
        remove(fn01[j1]);
        remove(fn02[j1]);
        remove(fn03[j1]);
      }
      Mat_SeqAIJ         *ac=(Mat_SeqAIJ*)submatC[j1]->data;//*aa=subA->data;
      ai= ac->i;
      aj= ac->j;
      vals=ac->a;
      nz=ac->nz;
      nrow=submatC[j1]->rmap->n;
      ncol=submatC[j1]->cmap->n;
      solve_real *be0 = (solve_real*)calloc(nrow,sizeof(solve_real));
      for(i=0; i<nrow-1; i++) {
        for(j=ai[i]; j<ai[i+1]; j++) {
          be0[i]+=vals[j]*xd[aj[j]];
        }
      }
      for(j=ai[i]; j<nz; j++) {
        be0[i]+=vals[j]*xd[aj[j]];
      }
      if(!dbbd_fastrefac)MatDestroy(&submatC[j1]);
      if(insize[j1*insizes+16]!=la1)insize[j1*insizes+16]=la1;
      spec48m_esol_(insize+j1*insizes,irne,vale,keep,be0,biui0);
      if(!dbbd_fastrefac) {
      free(irne);
      free(keep);
      free(vale);
      if(inmemory) {
        fac_irn[j1]=NULL;
        fac_keep[j1]=NULL;
        fac_va[j1]=NULL;
      }
      }
      free(be0);
      ISGetIndices(colindices[j1],&nindices);
      xi1point=xi1+xi1indx;
      for(i=0; i<block_sizes[j1+begblock[rank]]; i++) {
        x0[nindices[i]]+=xi1point[i]-biui0[i];//vals[i]-biui0[i];//vecval;
      }
      ISRestoreIndices(colindices[j1],&nindices);
      free(biui0);
      xi1indx+=block_sizes[j1+begblock[rank]];
    }
  }
  free(xi1);
  free(insize);
  if(!dbbd_fastrefac)ierr = PetscFree(submatC);
  free(begblock);
  logmsg(2,"Solution calculation rank %d time %f\n",rank,((double)clock()-timestr)/CLOCKS_PER_SEC);
  MPI_Barrier(PETSC_COMM_WORLD);
  free(xd);
  timestr=clock();
  if(SORD==1)MPI_Allreduce(x0,x1,VecSize, MPI_DOUBLE, MPI_SUM,PETSC_COMM_WORLD);
  else MPI_Allreduce(x0,x1,VecSize, MPI_FLOAT, MPI_SUM,PETSC_COMM_WORLD);
  logmsg(2,"Reduce solution rank %d time %f\n",rank,((double)clock()-timestr)/CLOCKS_PER_SEC);
  free(x0);
  for (i=0; i<nmatinplus; i++){
    free(fn01[i]);
    free(fn02[i]);
    free(fn03[i]);
  }
  free(fn01);
  free(fn02);
  free(fn03);
  if(inmemory) {
    free(fac_irn);
    free(fac_keep);
    free(fac_va);
  }
  for (i=0; i<nmatin; i++) {
    ierr = ISDestroy(&colindices[i]);
    CHKERRQ(ierr);
  }
  PetscFree(colindices);
  return 0;
}

int ndbbd_presolve(Mat A, Vec b, solve_real *x1, offset_t VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int *row_order,int *col_order, offset_t ndblock,offset_t nreg,offset_t ntime,int *block_sizes, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,dim_t laDi,dim_t laD,PetscReal cntl3,PetscReal cntl6,PetscBool presol) {//,bool iter
  IS *rowindices=NULL,*colindices=NULL,*colindicesbc1=NULL,*rowBBij=NULL,*colBBij=NULL;//,*colindicesbcpm,*colindicesbcpm1
  PetscInt bfirst,bend,sblockin,nmatin,nmatinplus,nmatint,nmatinplust,nmatminust,nrowcolin,i,i1,j,j0,j1,j2,j3,j4,j5,j6,j7,l0,l1,l2,l3,l4,l5,rank1,proc1=0,nnzmax,j1nz,j1irnbs;//,sumrowcolin
  Mat *submatAij=NULL,*submatBBij=NULL;//,*submatCij,*submatBij;,*submatB
  Vec vecxd;
  PetscInt *ai,*aj,*a1i,*a1j,*a2i,*a2j,*aic,*ajc;
  PetscReal cntl6in;
  offset_t lasize,ldsize;
  size_t freadresult,fwrt;
  PetscInt nrow,ncol,nz,nrowc,ncolc,nrowb,ncolb,nzc,nzv;
  MPI_Status   status;
  int j2int,la1,fd1,fd2,fd3,unequal=0;
  bool ifremove=false;
  PetscScalar *vals,*val1s,*val2s,*valsc,vecval;
  solve_real vval;
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
  nmatint=(offset_t)ntime/mpisize;
  if(nmatint*mpisize<ntime)unequal=1;
  nmatminust=nmatint;
  if(nmatint*mpisize<ntime)nmatinplust=nmatint+1;
  else nmatinplust=nmatint;
  for(i=0; i<mpisize; i++)if(rank+1<=ntime-mpisize*nmatint)nmatint++;
  nmatin=(nreg+1)*nmatint;
  nmatinplus=(nreg+1)*nmatinplust;
  if(inmemory||nfr_flag())ndbbd_fac_init(nmatin);
  begblock[rank]=nmatin;
  logmsg(2,"rank %d nmatin %d nmatint %d nmplus %d\n",rank,nmatin,nmatint,nmatinplus);
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
      strcpy(filename,scratch_dir);strcat(filename,"_irnv");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      strcpy(fn01[j1],filename);
      strcpy(filename,scratch_dir);strcat(filename,"_keep");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      strcpy(fn02[j1],filename);
      strcpy(filename,scratch_dir);strcat(filename,"_vav");
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
    PetscInt *indices= (PetscInt *) calloc (block_sizes[i+begblock[rank]],sizeof(PetscInt));
    for(j=0; j<block_sizes[i+begblock[rank]]; j++) {
      indices[j]=row_order[bfirst+j];//+bfirst;
    }

    ISCreateGeneral(PETSC_COMM_SELF,block_sizes[i+begblock[rank]],indices,PETSC_COPY_VALUES,rowindices+i);
    bfirst=countvarintra1[i+begblock[rank]];
    for(j=0; j<block_sizes[i+begblock[rank]]; j++) {
      indices[j]=col_order[bfirst+j];//+bfirst;
    }
    ISCreateGeneral(PETSC_COMM_SELF,block_sizes[i+begblock[rank]],indices,PETSC_COPY_VALUES,colindices+i);
    free(indices);//1
  }
  for(i=0; i<nmatint; i++) {
    j2=0;
    for(j1=0; j1<nreg+1; j1++)j2+=block_sizes[i*(nreg+1)+j1+begblock[rank]];
    PetscInt *indices= (PetscInt *) calloc (j2,sizeof(PetscInt));
    j2=0;
    for(j1=0; j1<nreg; j1++) {
      bfirst=countvarintra1[i*(nreg+1)+j1+begblock[rank]];
      for(j=0; j<block_sizes[i*(nreg+1)+j1+begblock[rank]]; j++) {
        indices[j2]=col_order[bfirst+j];//+bfirst;
        j2++;
      }
    }
    ISCreateGeneral(PETSC_COMM_SELF,j2,indices,PETSC_COPY_VALUES,colindicesbc1+i);
    free(indices);
  }
  logmsg(2,"Begin preparation rank %d\n",rank);
  ierr = MatCreateSubMatrices(A,nmatin,rowindices,colindices,MAT_INITIAL_MATRIX,&submatAij);
  int insizes=17;
  int *insize=(int *) calloc (insizes*nmatin,sizeof(int));
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
      ISCreateGeneral(PETSC_COMM_SELF,block_sizes[j1+begblock[rank]],indicesbbij+j3,PETSC_COPY_VALUES,colBBij);
      ierr = MatCreateSubMatrices(submatBBij[0],1,rowBBij,colBBij,MAT_INITIAL_MATRIX,&submatBij[j1]);
      ierr = ISDestroy(&colBBij[0]);
      j3+=block_sizes[j1+begblock[rank]];
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
    ISCreateGeneral(PETSC_COMM_SELF,block_sizes[j1+begblock[rank]],indicesbbij+j3,PETSC_COPY_VALUES,colBBij);
    ierr = MatCreateSubMatrices(submatBBij[0],1,rowBBij,colBBij,MAT_INITIAL_MATRIX,&submatBij[j1]);
    ierr = ISDestroy(&colBBij[0]);
    j3+=block_sizes[j1+begblock[rank]];
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
  solve_real *vecbivi=NULL;
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
    if(bivinzrow1==NULL)printf("Error: memory allocation failed for bivinzrow1\n");
    bivinzcol1=realloc (bivinzcol1,bivicolsizemax*sizeof(long int));
    if(bivinzcol1==NULL)printf("Error: memory allocation failed for bivinzcol1\n");
      vecbivi= (solve_real *) malloc (vecbivisizemax*sizeof(solve_real));
    if(vecbivi==NULL)printf("Error: memory allocation failed for vecbivi\n");
      irn=(int *) malloc (lasizemax*sizeof(int));
    if(irn==NULL)printf("Error: memory allocation failed for irn\n");
      jcn=(int *) malloc (lasizemax*sizeof(int));
    if(jcn==NULL)printf("Error: memory allocation failed for jcn\n");
      solve_real *values= (solve_real *) calloc (lasizemax,sizeof(solve_real));
    if(values==NULL)printf("Error: memory allocation failed for values\n");
      int* jcnb1=(int *) malloc (nebmax*sizeof(int));
    if(jcnb1==NULL)printf("Error: memory allocation failed for jcnb1\n");
      solve_real *sol48= (solve_real *) malloc (nrowamax*sizeof(solve_real));
    if(sol48==NULL)printf("Error: memory allocation failed for sol48\n");
      solve_real *b48= (solve_real *) malloc (nrowamax*sizeof(solve_real));
    if(b48==NULL)printf("Error: memory allocation failed for b48\n");
      int* keep=NULL,*iw51=NULL;
  keep=(int *) malloc ((nrowamax+5*ncolamax+4*ncolamax+7)*sizeof(int));//ICNTL(6)=1
    if(keep==NULL)printf("Error: memory allocation failed for keep\n");
  solve_real *w51=NULL;
  if(nrowamax>ncolamax) w51=(solve_real *) malloc ((5*nrowamax)*sizeof(solve_real));
  else w51=(solve_real *) malloc ((5*ncolamax)*sizeof(solve_real));
    if(w51==NULL)printf("Error: memory allocation failed for w51\n");
  iw51=(int *) malloc ((6*nrowamax+3*ncolamax)*sizeof(int));
    if(iw51==NULL)printf("Error: memory allocation failed for iw51\n");
  #pragma omp for schedule (static)
  for(j3=0; j3<nmatint; j3++) {
    bivirowsize=1;
    bivicolsize=1;
    bbrowij=submatBij[j3*(nreg+1)][0]->rmap->n;
    logmsg(2,"j3 %d bbrowij %djth %d windx %d\n",j3,bbrowij,jthrd,windx);
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
    logmsg(2,"rank %d vecbivisize %ld bivirowsize %d bivicolsize %d\n",rank,vecbivisize,bivirowsize,bivicolsize);
    if(vecbivisize<0)printf("Error: interface problem size overflow (negative vecbivisize); the border is too large\n");
    memset(vecbivi,0,vecbivisize*sizeof(solve_real));
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
      int nfr_redo=0;
      if(nfr_flag()) {
        if(ndbbd_fac_nz[j4]==nz&&ndbbd_fac_jcn[j4]!=NULL) {
          /* pattern unchanged: refill values, keep the pivot sequence */
          memcpy (ndbbd_fac_va[j4],vals,nz*sizeof(solve_real));
          nfr_redo=1;
        }
        else {
          ndbbd_fac_drop(j4);
          ndbbd_fac_irn[j4]=(int *) calloc (lasize,sizeof(int));
          ndbbd_fac_jcn[j4]=(int *) calloc (lasize,sizeof(int));
          ndbbd_fac_va[j4]=(solve_real *) calloc (lasize,sizeof(solve_real));
          ndbbd_fac_keep[j4]=(int *) calloc (nrow+9*ncol+7,sizeof(int));
          memcpy (ndbbd_fac_irn[j4]+nz,ai,(nrow+1)*sizeof(PetscInt));
          memcpy (ndbbd_fac_jcn[j4],aj,nz*sizeof(PetscInt));
          memcpy (ndbbd_fac_va[j4],vals,nz*sizeof(solve_real));
          ndbbd_fac_nz[j4]=nz;
        }
      }
      else {
        memcpy (irn+nz,ai,(nrow+1)*sizeof(PetscInt));
        memcpy (jcn,aj,nz*sizeof(PetscInt));
        memcpy (values,vals,nz*sizeof(solve_real));
        MatDestroy(&submatAij[j4]);//1
      }
      insize[j4*insizes+2]=nz;
      insize[j4*insizes+5]=nzc;
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
      if(nfr_flag()) {
        int redo_io=nfr_redo;
        prep48m_msol_p_(insize+j4*insizes,ndbbd_fac_irn[j4],ndbbd_fac_jcn[j4],ndbbd_fac_va[j4],aic,ajc,valsc,ai,aj,vals,vecbivi,bivinzrow,bivinzcol,jcnb1,sol48,b48,w51,iw51,ndbbd_fac_keep[j4],&redo_io);
        if(redo_io<0) {
          /* fast refactorize declined: re-stage from the still-live
             submatrix and redo the analyse */
          Mat_SeqAIJ *aa2=(Mat_SeqAIJ*)submatAij[j4]->data;
          memcpy (ndbbd_fac_irn[j4]+ndbbd_fac_nz[j4],aa2->i,(nrow+1)*sizeof(PetscInt));
          memcpy (ndbbd_fac_jcn[j4],aa2->j,ndbbd_fac_nz[j4]*sizeof(PetscInt));
          memcpy (ndbbd_fac_va[j4],aa2->a,ndbbd_fac_nz[j4]*sizeof(solve_real));
          redo_io=0;
          prep48m_msol_p_(insize+j4*insizes,ndbbd_fac_irn[j4],ndbbd_fac_jcn[j4],ndbbd_fac_va[j4],aic,ajc,valsc,ai,aj,vals,vecbivi,bivinzrow,bivinzcol,jcnb1,sol48,b48,w51,iw51,ndbbd_fac_keep[j4],&redo_io);
        }
        MatDestroy(&submatAij[j4]);//1
      }
      else {
      prep48m_msol_(insize+j4*insizes,irn,jcn,values,aic,ajc,valsc,ai,aj,vals,vecbivi,bivinzrow,bivinzcol,jcnb1,sol48,b48,w51,iw51,keep);
      ndbbd_fac_emit(rank,j4,irn,keep,values,insize[j4*insizes+16],insize[j4*insizes+12]);
      }
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
    if(irn1==NULL)printf("Error: memory allocation failed for irn1\n");
    jcn1=realloc(jcn1,irnmems*sizeof(int));
    if(jcn1==NULL)printf("Error: memory allocation failed for jcn1\n");

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
    logmsg(2,"D ncol %d nrow %d nz %d\n",ncol,nrow,nz);
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
    strcpy(filename,scratch_dir);strcat(filename,"_bivi");
    strcat(filename,rankname);
    strcat(filename,j1name);
    strcat(filename,".bin");
    if((presolfile=fopen(filename, "wb"))==NULL) {
      printf("Error: cannot open scratch file %s\n",filename);
    }
    fwrt=fwrite(vecbivi, sizeof(solve_real),nz1, presolfile);
    if(fwrt==0)printf("Error: short write on scratch file %s\n",filename);
    fclose(presolfile);
    strcpy(filename,scratch_dir);strcat(filename,"_rbvi");
    strcat(filename,rankname);
    strcat(filename,j1name);
    strcat(filename,".bin");
    if((presolfile=fopen(filename, "wb"))==NULL) {
      printf("Error: cannot open scratch file %s\n",filename);
    }
    fwrt=fwrite(irn1, sizeof(int), nz1, presolfile);
    if(fwrt== 0) printf("Error: short write on scratch file %s\n",filename);
    fclose(presolfile);
    strcpy(filename,scratch_dir);strcat(filename,"_cbvi");
    strcat(filename,rankname);
    strcat(filename,j1name);
    strcat(filename,".bin");
    if((presolfile=fopen(filename, "wb"))==NULL) {
      printf("Error: cannot open scratch file %s\n",filename);
    }
    fwrt=fwrite(jcn1, sizeof(int), nz1, presolfile);
    if(fwrt== 0) printf("Error: short write on scratch file %s\n",filename);
    fclose(presolfile);
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
      strcpy(filename,scratch_dir);strcat(filename,"_rank");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "wb"))==NULL) {
        printf("Error: cannot open scratch file %s\n",filename);
      }
      logmsg(2,"rank %d nrow %d ncol %d block rank %d nz %d cntl6 %lf\n",rank,nrow,ncol,insized[3],nz,cntl6);
      fwrt=fwrite(insized, sizeof(int), 5+nreg*insizes, presolfile);
      if(fwrt==0)printf("Error: short write on scratch file %s\n",filename);
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
  if(irn1==NULL)printf("Error: memory allocation failed for irn1\n");
  jcn1=realloc(jcn1,ldsize*sizeof(int));
  if(jcn1==NULL)printf("Error: memory allocation failed for jcn1\n");
  irn=realloc(irn,nrowmax*sizeof(int));
  if(irn==NULL)printf("Error: memory allocation failed for irn\n");
  jcn=realloc(jcn,ncolmax*sizeof(int));
  if(jcn==NULL)printf("Error: memory allocation failed for jcn\n");
  keep=(int *) malloc ((nrowmax+5*ncolmax+4*ncolmax+7)*sizeof(int));//ICNTL(6)=1
  if(keep==NULL)printf("Error: memory allocation failed for keep\n");
  if(nrowmax>ncolmax) w51=(solve_real *) malloc ((5*nrowmax)*sizeof(solve_real));
  else w51=(solve_real *) malloc ((5*ncolmax)*sizeof(solve_real));
  if(w51==NULL)printf("Error: memory allocation failed for w51\n");
  iw51=(int *) malloc ((6*nrowmax+3*ncolmax)*sizeof(int));
  if(iw51==NULL)printf("Error: memory allocation failed for iw51\n");
  vecbivi=realloc(vecbivi,ldsize*sizeof(solve_real));
  if(vecbivi==NULL)printf("Error: memory allocation failed for vecbivi\n");
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
    strcpy(filename,scratch_dir);strcat(filename,"_rank");
    strcat(filename,rankname);
    strcat(filename,j1name);
    strcat(filename,".bin");
    if((presolfile=fopen(filename, "r"))==NULL) {
      printf("Error: cannot open interface file %s; the presolve pass did not produce it (check scratch space)\n",filename);
    }
    fwrt=fread(insizeda, sizeof(int), 5+nreg*insizes, presolfile);
    fclose(presolfile);
    remove(filename);
    if( fwrt== 0)printf("Error: short read on interface file %s; the presolve pass did not produce it (check scratch space)\n",filename);
    nrow=insizeda[0];
    ncol=insizeda[1];
    nz1=insizeda[2];
    nz=nz1;
      strcpy(filename,scratch_dir);strcat(filename,"_bivi");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "r"))==NULL) {
        printf("Error: cannot open scratch file %s\n",filename);
      }
      fwrt=fread(vecbivi, sizeof(solve_real), insizeda[2], presolfile);
      if(fwrt== 0) printf("Error: short read on scratch file %s (%ld of %d records)\n",filename,fwrt,insizeda[2]);
      fclose(presolfile);
      strcpy(filename,scratch_dir);strcat(filename,"_rbvi");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "r"))==NULL) {
        printf("Error: cannot open scratch file %s\n",filename);
      }
      fwrt=fread(irn1, sizeof(int), insizeda[2], presolfile);
      if(fwrt== 0) printf("Error: short read on scratch file %s\n",filename);
      fclose(presolfile);
      strcpy(filename,scratch_dir);strcat(filename,"_cbvi");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "r"))==NULL) {
        printf("Error: cannot open scratch file %s\n",filename);
      }
      fwrt=fread(jcn1, sizeof(int), insizeda[2], presolfile);
      if(fwrt== 0) printf("Error: short read on scratch file %s\n",filename);
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
      strcpy(filename,scratch_dir);strcat(filename,"_rank");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "wb"))==NULL) {
        printf("Error: cannot open scratch file %s\n",filename);
      }
      logmsg(2,"rank %d nrow %d ncol %d block rank %d nz %d cntl6 %lf\n",rank,nrow,ncol,insized[3],nz,cntl6in);
      fwrt=fwrite(insized, sizeof(int), 5+nreg*insizes, presolfile);
      if(fwrt==0)printf("Error: short write on scratch file %s\n",filename);
      fclose(presolfile);
      strcpy(filename,scratch_dir);strcat(filename,"_row");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "wb"))==NULL) {
        printf("Error: cannot open scratch file %s\n",filename);
      }
      fwrt=fwrite(irn, sizeof(int), nrow, presolfile);
      if(fwrt== 0) printf("Error: short write on scratch file %s\n",filename);
      fclose(presolfile);
      strcpy(filename,scratch_dir);strcat(filename,"_col");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "wb"))==NULL) {
        printf("Error: cannot open scratch file %s\n",filename);
      }
      fwrt=fwrite(jcn, sizeof(int), ncol, presolfile);
      if(fwrt== 0) printf("Error: short write on scratch file %s\n",filename);
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
  logmsg(2,"End of preparation rank %d\n",rank);
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


int ndbbd_solve(Mat A, Vec b, solve_real *x1, offset_t VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int *row_order,int *col_order, offset_t ndblock,offset_t nreg,offset_t ntime,int *block_sizes, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,dim_t laDi,dim_t laD,PetscReal cntl3,PetscReal cntl6,PetscBool presol) {//,bool iter
  IS *rowindices=NULL,*colindices=NULL,*rowindicesbc=NULL,*colindicesbc1=NULL,*colindicesbc2=NULL,*Cindices=NULL,*Bindices=NULL,Cindicesc,Bindicesc,*BBindices=NULL,*BBiindices=NULL,*rowBBij=NULL,*colBBij=NULL;//,*colindicesbcpm,*colindicesbcpm1
  const PetscInt *nindices;
  PetscInt bfirst,bend,sblockin,nmatin,nmatinplus,nmatint,nmatinplust,nmatminust,nrowcolin,sumrowcolin,i,i1,j,j0,j1,j2,j3,j4,j5,j6,j7,l0,l1,l2,l3,l4,l5,rank1,proc1=0,nnzmax,j1nz,j1irnbs;
  Mat *submatAij=NULL,*submatC=NULL,*submatD=NULL,*submatBB=NULL,*submatBBij=NULL;//,*submatCij,*submatBij;,*submatB
  Vec vecxd;
  PetscInt *ai,*aj,*a1i,*a1j,*a2i,*a2j,*aic,*ajc;
  PetscReal cntl6in;
  offset_t ldsize;
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
  solve_real vval;
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
  nmatint=(offset_t)ntime/mpisize;
  if(nmatint*mpisize<ntime)unequal=1;
  nmatminust=nmatint;
  if(nmatint*mpisize<ntime)nmatinplust=nmatint+1;
  else nmatinplust=nmatint;
  for(i=0; i<mpisize; i++)if(rank+1<=ntime-mpisize*nmatint)nmatint++;
  nmatin=(nreg+1)*nmatint;
  nmatinplus=(nreg+1)*nmatinplust;
  if(inmemory||nfr_flag())ndbbd_fac_init(nmatin);
  begblock[rank]=nmatin;
  logmsg(2,"rank %d nmatin %d nmatint %d nmplus %d\n",rank,nmatin,nmatint,nmatinplus);
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
  for(i=0; i<mpisize; i++)logmsg(2,"rank %d beg block %d\n",rank,begblock[i]);
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

  solve_real **yi2= (solve_real**)calloc(nmatin,sizeof(solve_real*));
  for (i=0; i<nmatin; i++) {
    yi2[i] = (solve_real*)calloc(block_sizes[i+begblock[rank]],sizeof(solve_real));
  }
  sumrowcolin=0;
  for(j=0; j<ndblock; j++) {
    sumrowcolin+=block_sizes[j];
  }
  logmsg(2,"VecSize %ld sumrowcolin %d\n",VecSize,sumrowcolin);
  PetscInt *offblock= (PetscInt *) calloc (ndblock+1,sizeof(PetscInt));
  PetscInt *offblockrow= (PetscInt *) calloc (ndblock+1,sizeof(PetscInt));
  offblock[0]=0;
  for(j=0; j<ndblock; j++)offblock[j+1]=offblock[j]+countvarintra1[j+1]-countvarintra1[j]-block_sizes[j];
  offblockrow[0]=0;
  for(j=0; j<ndblock; j++)offblockrow[j+1]=offblockrow[j]+counteqnoadd[j]-block_sizes[j];

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
      strcpy(filename,scratch_dir);strcat(filename,"_irnv");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      strcpy(fn01[j1],filename);
      strcpy(filename,scratch_dir);strcat(filename,"_keep");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      strcpy(fn02[j1],filename);
      strcpy(filename,scratch_dir);strcat(filename,"_vav");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      strcpy(fn03[j1],filename);
    }
  }
  for(j=0; j<mpisize; j++) {
    for(j3=0; j3<nmatinplus; j3++) {
      if(j3>=nmatin)i=nmatin-1;
      else i=j3;
      bfirst=counteq[i+begblock[rank]];
      bend=block_sizes[i+begblock[rank]];
      MPI_Bcast(&bfirst,1, MPI_INT,j, PETSC_COMM_WORLD);
      MPI_Bcast(&bend,1, MPI_INT,j, PETSC_COMM_WORLD);
      PetscInt *indices= (PetscInt *) calloc (bend,sizeof(PetscInt));
      solve_real *yi0 = (solve_real*)calloc(bend,sizeof(solve_real));
      for(j1=0; j1<bend; j1++) { //ha_ndblocks[i+begblock[j]]
        j2=row_order[bfirst+j1];//+bfirst;
        if(j2>=Istart&&j2<Iend)indices[j1]=j2;
        else indices[j1]=-1;
      }
      VecGetValues(b,bend,indices,yi0);
      reduce_to_rank(yi0,(fortran_int)bend,mpisize,rank,j);
      if(rank==j&&i<nmatin) {
        memcpy (yi2[i],yi0,bend*sizeof(solve_real));
      }
      free(indices);
      free(yi0);
    }
  }
  solve_real **yi1= (solve_real**)calloc(nmatint,sizeof(solve_real*));
  for (i=0; i<nmatint; i++) {
    j1=0;
    for (j=0; j<nreg+1; j++) {
      j1+=block_sizes[i*(nreg+1)+j+begblock[rank]];
    }
    yi1[i] = (solve_real*)calloc(j1,sizeof(solve_real));
  }
  for (i=0; i<nmatint; i++) {
    j1=0;
    for (j=0; j<nreg+1; j++) {
      memcpy (&yi1[i][j1],yi2[i*(nreg+1)+j],block_sizes[i*(nreg+1)+j+begblock[rank]]*sizeof(solve_real));
      free(yi2[i*(nreg+1)+j]);
      j1+=block_sizes[i*(nreg+1)+j+begblock[rank]];
    }
  }
  free(yi2);
  j1=0;
  for(i=0; i<nmatin; i++) {
    bfirst=counteq[i+begblock[rank]];
    PetscInt *indices= (PetscInt *) calloc (block_sizes[i+begblock[rank]],sizeof(PetscInt));
    for(j=0; j<block_sizes[i+begblock[rank]]; j++) {
      indices[j]=row_order[bfirst+j];//+bfirst;
    }

    ISCreateGeneral(PETSC_COMM_SELF,block_sizes[i+begblock[rank]],indices,PETSC_COPY_VALUES,rowindices+i);
    bfirst=countvarintra1[i+begblock[rank]];
    for(j=0; j<block_sizes[i+begblock[rank]]; j++) {
      indices[j]=col_order[bfirst+j];//+bfirst;
    }
    ISCreateGeneral(PETSC_COMM_SELF,block_sizes[i+begblock[rank]],indices,PETSC_COPY_VALUES,colindices+i);
    free(indices);
  }
  for(i=0; i<nmatint; i++) {
    j2=0;
    for(j1=0; j1<nreg+1; j1++)j2+=block_sizes[i*(nreg+1)+j1+begblock[rank]];
    PetscInt *indices= (PetscInt *) calloc (j2,sizeof(PetscInt));
    j2=0;
    for(j1=0; j1<nreg+1; j1++) {
      bfirst=counteq[i*(nreg+1)+j1+begblock[rank]];
      for(j=0; j<block_sizes[i*(nreg+1)+j1+begblock[rank]]; j++) {
        indices[j2]=row_order[bfirst+j];//+bfirst;
        j2++;
      }
    }
    ISCreateGeneral(PETSC_COMM_SELF,j2,indices,PETSC_COPY_VALUES,rowindicesbc+i);
    j2=0;
    for(j1=0; j1<nreg; j1++) {
      bfirst=countvarintra1[i*(nreg+1)+j1+begblock[rank]];
      for(j=0; j<block_sizes[i*(nreg+1)+j1+begblock[rank]]; j++) {
        indices[j2]=col_order[bfirst+j];//+bfirst;
        j2++;
      }
    }
    ISCreateGeneral(PETSC_COMM_SELF,j2,indices,PETSC_COPY_VALUES,colindicesbc1+i);
    j2=0;
    for(j1=nreg; j1<nreg+1; j1++) {
      bfirst=countvarintra1[i*(nreg+1)+j1+begblock[rank]];
      for(j=0; j<block_sizes[i*(nreg+1)+j1+begblock[rank]]; j++) {
        indices[j2]=col_order[bfirst+j];//+bfirst;
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
      indicesC[j]=col_order[bend-j1-1];//+bfirst;
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
      indicesB[j]=row_order[bend-j1-1];//+bfirst;
    }
  }
  j1=offblockrow[ndblock];
  for(i=counteq[ndblock-1]+counteqnoadd[ndblock-1]; i<VecSize; i++) {
    indicesB[j1]=i;
    j1++;
  }
  /* this rank's [Istart,Iend) slice of the border rows; a rank owning
     no border rows must yield an empty range, not stale j/j1 */
  j=VecSize-sumrowcolin;
  for(i=0; i<VecSize-sumrowcolin; i++)if(indicesB[i]>=Istart) {
      j=i;
      break;
    }
  j1=VecSize-sumrowcolin;
  for(j2=j; j2<VecSize-sumrowcolin; j2++)if(indicesB[j2]>=Iend) {
      j1=j2;
      break;
    }
  logmsg(2,"rank %d j %d j1 %d istart %d iend %d\n",rank,j,j1,Istart,Iend);
  solve_real *vecbiui= (solve_real *) calloc (VecSize-sumrowcolin,sizeof(solve_real));
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
  logmsg(2,"Begin preparation rank %d\n",rank);
  int insizes=17;
  int *insize=(int *) calloc (insizes*nmatin,sizeof(int));
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
      ISCreateGeneral(PETSC_COMM_SELF,block_sizes[j1+begblock[rank]],indicesbbij+j3,PETSC_COPY_VALUES,colBBij);
      ierr = MatCreateSubMatrices(submatBBij[0],1,rowBBij,colBBij,MAT_INITIAL_MATRIX,&submatBij[j1]);
      ierr = ISDestroy(&colBBij[0]);
      j3+=block_sizes[j1+begblock[rank]];
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
    ISCreateGeneral(PETSC_COMM_SELF,block_sizes[j1+begblock[rank]],indicesbbij+j3,PETSC_COPY_VALUES,colBBij);
    ierr = MatCreateSubMatrices(submatBBij[0],1,rowBBij,colBBij,MAT_INITIAL_MATRIX,&submatBij[j1]);
    ierr = ISDestroy(&colBBij[0]);
    j3+=block_sizes[j1+begblock[rank]];
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
  ierr = MatCreateSubMatrices(A,nmatint,rowindicesbc,Cindices,MAT_INITIAL_MATRIX,&submatC);
  ierr = MatCreateSubMatrices(A,nmatinBB,Bindices,BBindices,MAT_INITIAL_MATRIX,&submatBB);
  ierr = MatDestroy(&A);
  CHKERRQ(ierr);
  PetscInt BBrow,BBcol;//,CCrow,CCcol;
  MatGetSize(submatBB[0],&BBrow,&BBcol);
  logmsg(2,"rank %d bbrow %d bbcol %d\n",rank,BBrow,BBcol);//,CCrow,CCcol);
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
  if(rank==mpisize-1)proc1=1;
  MatCreateSubMatrices(submatBB[0],proc1,Bindices,Cindices,MAT_INITIAL_MATRIX,&submatD);
  logmsg(2,"End partitioning rank %d time %f\n",rank,((double)(clock()-timestr))/CLOCKS_PER_SEC);
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
  long int *bivinzrow=NULL;//(long int *) calloc (1,sizeof(long int));
  PetscInt *bivinzcol=NULL;//(PetscInt *) calloc (1,sizeof(PetscInt));
  PetscReal *ccolnorms=NULL;//(PetscReal *) calloc (1,sizeof(PetscReal));
  int bivirowsize=1,bivicolsize=1,bbrowij,ddrowi;
  long int vecbivisize=0,li,lj;//,halfvec;
  solve_real *xi1point;
  long int *bivinzrow1=NULL;//(long int *) calloc (1,sizeof(long int));
  long int *bivinzcol1=NULL;//(long int *) calloc (1,sizeof(long int));
  offset_t xi1indx=0;
  long int nz0=0,nz1,nz3;//,nz2,halfj2;
  int fd1,fd2,fd3,frrsl1,frrsl2,frrsl3;
  omp_set_num_threads(section_threads);
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
    strcpy(filename,scratch_dir);strcat(filename,"_rank");
    strcat(filename,rankname);
    strcat(filename,j1name);
    strcat(filename,".bin");
    #pragma omp critical
    {
    if((presolfile=fopen(filename, "r"))==NULL) {
      printf("Error: cannot open interface file %s; the presolve pass did not produce it (check scratch space)\n",filename);
    }
    frd=fread(insizeda, sizeof(int), 5+nreg*insizes, presolfile);
    fclose(presolfile);
    if( frd== 0)printf("Error: short read on interface file %s; the presolve pass did not produce it (check scratch space)\n",filename);
    }
    if(insizeda0<insizeda[0])insizeda0=insizeda[0];
    if(insizeda1<insizeda[1])insizeda1=insizeda[1];
    if(insizeda2<insizeda[2])insizeda2=insizeda[2];
  }
  ldsize=ceil((laDi/100.0)*insizeda2);
  ldsize+=10;
  insizeda2=ldsize;
  int *irn1=(int *) malloc(ldsize*sizeof(int));
  if(irn1==NULL)printf("Error: memory allocation failed for irn1\n");
  int *jcn1=(int *) malloc(ldsize*sizeof(int));
  if(jcn1==NULL)printf("Error: memory allocation failed for jcn1\n");
  solve_real *vecbivi0= (solve_real *) malloc (ldsize*sizeof(solve_real));
  if(vecbivi0==NULL)printf("Error: memory allocation failed for vecbivi0\n");
  int *irn1a=(int *) malloc (insizeda0*sizeof(int));
  if(irn1a==NULL)printf("Error: memory allocation failed for irn1a\n");
  int *jcn1a=(int *) malloc (insizeda1*sizeof(int));
  if(jcn1a==NULL)printf("Error: memory allocation failed for jcn1a\n");
  solve_real *fw=(solve_real *) malloc (insizeda0*sizeof(solve_real));
  if(fw==NULL)printf("Error: memory allocation failed for fw\n");
  int *fiw=(int *) malloc ((6*insizeda0+3*insizeda1)*sizeof(int));
  if(fiw==NULL)printf("Error: memory allocation failed for fiw\n");
  int *fkeep=(int *) malloc ((insizeda0+5*insizeda1+4*insizeda1/1+7)*sizeof(int));//ICNTL(6)=1
  if(fkeep==NULL)printf("Error: memory allocation failed for fkeep\n");
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
    strcpy(filename,scratch_dir);strcat(filename,"_rank");
    strcat(filename,rankname);
    strcat(filename,j1name);
    strcat(filename,".bin");
    if((presolfile=fopen(filename, "r"))==NULL) {
      printf("Error: cannot open interface file %s; the presolve pass did not produce it (check scratch space)\n",filename);
    }
    frd=fread(insizeda, sizeof(int), 5+nreg*insizes, presolfile);
    fclose(presolfile);
    remove(filename);
    if( frd== 0)printf("Error: short read on interface file %s; the presolve pass did not produce it (check scratch space)\n",filename);
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
    strcpy(filename,scratch_dir);strcat(filename,"_row");
    strcat(filename,rankname);
    strcat(filename,j1name);
    strcat(filename,".bin");
    if((presolfile=fopen(filename, "r"))==NULL) {
      printf("Error: cannot open scratch file %s\n",filename);
    }
    frd=fread(irn1a, sizeof(int), insizeda[0], presolfile);
    if(frd == 0)printf("Error: short read on scratch file %s\n",filename);
    fclose(presolfile);
    remove(filename);
    strcpy(filename,scratch_dir);strcat(filename,"_col");
    strcat(filename,rankname);
    strcat(filename,j1name);
    strcat(filename,".bin");
    if((presolfile=fopen(filename, "r"))==NULL) {
      printf("Error: cannot open scratch file %s\n",filename);
    }
    frd=fread(jcn1a, sizeof(int), insizeda[1], presolfile);
    if(frd== 0) printf("Error: short read on scratch file %s\n",filename);
    fclose(presolfile);
    remove(filename);
    strcpy(filename,scratch_dir);strcat(filename,"_bivi");
    strcat(filename,rankname);
    strcat(filename,j1name);
    strcat(filename,".bin");
    if((presolfile=fopen(filename, "r"))==NULL) {
      printf("Error: cannot open scratch file %s\n",filename);
    }
    frd=fread(vecbivi0, sizeof(solve_real), insizeda[2], presolfile);
    if(frd== 0) printf("Error: short read on scratch file %s (%ld of %d records)\n",filename,frd,insizeda[2]);
    fclose(presolfile);
    remove(filename);
    strcpy(filename,scratch_dir);strcat(filename,"_rbvi");
    strcat(filename,rankname);
    strcat(filename,j1name);
    strcat(filename,".bin");
    if((presolfile=fopen(filename, "r"))==NULL) {
      printf("Error: cannot open scratch file %s\n",filename);
    }
    frd=fread(irn1, sizeof(int), insizeda[2], presolfile);
    if(frd== 0) printf("Error: short read on scratch file %s\n",filename);
    fclose(presolfile);
    remove(filename);
    strcpy(filename,scratch_dir);strcat(filename,"_cbvi");
    strcat(filename,rankname);
    strcat(filename,j1name);
    strcat(filename,".bin");
    if((presolfile=fopen(filename, "r"))==NULL) {
      printf("Error: cannot open scratch file %s\n",filename);
    }
    frd=fread(jcn1, sizeof(int), insizeda[2], presolfile);
    if(frd== 0) printf("Error: short read on scratch file %s\n",filename);
    fclose(presolfile);
    remove(filename);
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
    ndbbd_fac_emit(rank,j4,irn1,fkeep,vecbivi0,insize[j4*insizes+16],insize[j4*insizes+12]);
    logmsg(2,"rank %d nz %ld nzmax %ld ldsize %ld\n",rank,nz1,insizeda2,ldsize);
  }
  free(vecbivi0);
  free(irn1);
  free(jcn1);
  free(irn1a);
  free(jcn1a);
  free(fw);
  free(fiw);
  free(fkeep);
  free(insizeda);
  }
  omp_set_num_threads(max_threads);
  if(presol)
    return 0;
  logmsg(2,"End of preparation rank %d\n",rank);
  MPI_Barrier(PETSC_COMM_WORLD);

  ccolnorms=realloc(ccolnorms,BBrow*sizeof(PetscReal));
  bivinzrow=realloc(bivinzrow,BBrow*sizeof(long int));
  memset(bivinzrow,0,BBrow*sizeof(long int));
  bivinzcol=realloc(bivinzcol,BBrow*sizeof(PetscInt));
  memset(bivinzcol,0,BBrow*sizeof(PetscInt));
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
  if(vecbivisize<0)printf("Error: interface problem size overflow (negative vecbivisize); the border is too large\n");
  solve_real *vecbivi= (solve_real *) calloc (vecbivisize,sizeof(solve_real));
  logmsg(2,"vecbivisize %ld rank %d\n",vecbivisize,rank);
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
  solve_real *xi1 = (solve_real*)calloc(sumrowcolin,sizeof(solve_real));
  #pragma omp parallel private(jthrd,timestr,aic,ajc,valsc,nrowc,ncolc,a1i,a1j,val1s,nz,a2i,a2j,val2s,nrowb,ncolb,i,j,j2,xi1point,xi1indx,maxrowcij,la1,fp1,fp2,fp3,freadresult,frrsl1,frrsl2,frrsl3,longsize,nzc) shared(submatC,submatB1,submatB2,xi1,submatCij,submatBij,insize,yi1,vecbivi,vecbiui)
  {
  int icntl[20],info[20];
  solve_real cntl[10],rinfo[10],error1[3];
  solve_real *w= (solve_real*)calloc(4*maxcolc,sizeof(solve_real));
  int *iw= (int*)calloc(maxcolc,sizeof(int));
  int **irnereg= (int**)calloc(nreg+1,sizeof(int*));
  int **keepreg= (int**)calloc(nreg+1,sizeof(int*));
  solve_real **valereg = (solve_real**)calloc(nreg+1,sizeof(solve_real*));
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
    maxrowcij=0;
    #pragma omp critical
    {
    for(i=0; i<nreg+1; i++) {
      j2=j1*(nreg+1)+i;
      insize[j2*insizes+15]=j2%90+7;
      if(i!=nreg&&(submatCij[j2][0]->rmap->n)>maxrowcij)maxrowcij=submatCij[j2][0]->rmap->n;
      la1=ceil((insize[j2*insizes+9]/100.0)*insize[j2*insizes+2]);
      if(inmemory||nfr_flag()) {/* alias the resident factors; per-element frees are skipped below */
        irnereg[i]=ndbbd_fac_irn[j2];
        keepreg[i]=ndbbd_fac_keep[j2];
        valereg[i]=ndbbd_fac_va[j2];
      }
      else {
        nfp1[i] = fopen(fn01[j2], "rb");
        if (nfp1[i]==NULL)printf("Error: cannot open factor file %s\n",fn01[j2]);
        irnereg[i] = realloc(irnereg[i],la1*sizeof(int));//(int*)calloc(la1,sizeof(int));
        freadresult=fread(irnereg[i],sizeof(int),la1,nfp1[i]);
        fclose(nfp1[i]);
        nfp2[i] = fopen(fn02[j2], "rb");
        if (nfp2[i]==NULL)printf("Error: cannot open factor file %s\n",fn02[j2]);
        keepreg[i] = realloc(keepreg[i],insize[j2*insizes+12]*sizeof(int));//(int*)calloc(insize[j1*insizes+12],sizeof(int));
        freadresult=fread(keepreg[i],sizeof(int),insize[j2*insizes+12],nfp2[i]);
        fclose(nfp2[i]);
        nfp3[i] = fopen(fn03[j2], "rb");
        if (nfp3[i]==NULL)printf("Error: cannot open factor file %s\n",fn03[j2]);
        valereg[i] = realloc(valereg[i],la1*sizeof(solve_real));//(ha_cgetype*)calloc(la1,sizeof(ha_cgetype));
        freadresult=fread(valereg[i],sizeof(solve_real),la1,nfp3[i]);
        fclose(nfp3[i]);
      }
    }
    }
    solve_real *b02 = (solve_real*)calloc(maxrowcij,sizeof(solve_real));

    ndbbd_block_solve_mem(rank,j1*(nreg+1),nreg,insize,insizes,submatCij,submatBij,yi1[j1],xi1point,irnereg,keepreg,valereg,cntl,rinfo,error1,icntl,info,w,iw,b02);
    solve_real *bccol= (solve_real*)calloc(ncolc,sizeof(solve_real));
    longsize=ncolc*sizeof(solve_real);
    for(i=0; i<nrowc-1; i++) {
      if(aic[i]<aic[i+1]) {
        memset(yi1[j1],0,longsize);
        for(j=aic[i]; j<aic[i+1]; j++) {
          yi1[j1][ajc[j]]=valsc[j];
        }
        ndbbd_block_solve_mem(rank,j1*(nreg+1),nreg,insize,insizes,submatCij,submatBij,yi1[j1],bccol,irnereg,keepreg,valereg,cntl,rinfo,error1,icntl,info,w,iw,b02);
        spar_vbiviadd_(bccol,bivinzcol+i,bivinzrow,&vecbivisize,&nrowb,&ncolc,&nz,ai12,aj12,vals12,vecbivi);
      }
    }
    if(aic[i]<nzc) {
      memset(yi1[j1],0,longsize);
      for(j=aic[i]; j<nzc; j++) {
        yi1[j1][ajc[j]]=valsc[j];
      }
      ndbbd_block_solve_mem(rank,j1*(nreg+1),nreg,insize,insizes,submatCij,submatBij,yi1[j1],bccol,irnereg,keepreg,valereg,cntl,rinfo,error1,icntl,info,w,iw,b02);
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
    logmsg(2,"Submatrix %d rank %d calculation time %f\n",j1,rank,difftime(timeend,timestr)/CLOCKS_PER_SEC);
  }
  free(w);
  free(iw);
  for(i=0; i<nreg+1; i++) {
    if(inmemory||nfr_flag())break;/* elements alias the resident store */
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
  logmsg(2,"Completed calculation of partitioned matrices, rank %d\n",rank);
  free(bivinzcol);
  free(bivinzrow);
  free(submatB1);
  free(submatB2);
  CHKERRQ(ierr);
  free(yi1);
  logmsg(2,"Completed calculation of partitioned matrices, rank %d\n",rank);
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
  logmsg(2,"vbivi %ld colsize %d\n",vecbivisize,bivicolsize);
  for(li=0; li<vecbivisize; ++li) {
    if(vecbivi[li]!=0) {
      vecbivi[lj2]=vecbivi[li];
      biviindx1[lj2]=bivinzcol1[li%bivicolsize]+bivinzrow1[li/bivicolsize];//[i-j1*bivicolsize]+bivinzrow1[j1];//*nrowb;
      ++lj2;
    }
  }
  free(bivinzrow1);
  free(bivinzcol1);
  vecbivi=realloc(vecbivi,lj2*sizeof(solve_real));
  nz0=lj2;
  int group_size,node_rank,group_size1,group_size11,node_tail_rank,color;
  MPI_Comm_rank( node_comm, &node_rank);
  MPI_Comm_size(node_comm,&group_size);
  if(node_rank==group_size-1)color=1;
  else color=0;
  MPI_Comm_rank( node_tail_comm, &node_tail_rank);
  MPI_Comm_size(node_tail_comm,&group_size1);
  group_size11=group_size1;
  MPI_Bcast(&group_size11,1, MPI_INT,mpisize-1, PETSC_COMM_WORLD);
  long int sindx01,sindx02,sindx03;
  
  if(node_rank!=group_size-1) {
    MPI_Send(&lj2,1, MPI_LONG,group_size-1, 10, node_comm);
    sindx01=lj2/MAXSSIZE;
    sindx03=0;
    for(sindx02=0; sindx02<sindx01; sindx02++) {
    if(SORD==1)MPI_Send(vecbivi+sindx03,MAXSSIZE, MPI_DOUBLE,group_size-1, 11*sindx01+sindx02, node_comm);
    else MPI_Send(vecbivi+sindx03,MAXSSIZE, MPI_FLOAT,group_size-1, 11*sindx01+sindx02, node_comm);
    sindx03+=MAXSSIZE;
    }
    if(SORD==1)MPI_Send(vecbivi+sindx03,lj2-sindx01*MAXSSIZE, MPI_DOUBLE,group_size-1, 11, node_comm);
    else MPI_Send(vecbivi+sindx03,lj2-sindx01*MAXSSIZE, MPI_FLOAT,group_size-1, 11, node_comm);
    sindx03=0;
    for(sindx02=0; sindx02<sindx01; sindx02++) {
    MPI_Send(biviindx1+sindx03,MAXSSIZE, MPI_LONG,group_size-1, 12*sindx01+sindx02, node_comm);
    sindx03+=MAXSSIZE;
    }
    MPI_Send(biviindx1+sindx03,lj2-sindx01*MAXSSIZE, MPI_LONG,group_size-1, 12, node_comm);
    free(obiviindx1);
    obiviindx1=NULL;
    free(vecbivi);
    vecbivi=NULL;
  }
  if(node_rank==group_size-1) {
    for(j6=0; j6<group_size-1; j6++) {
      MPI_Recv(&lj2,1, MPI_LONG,  MPI_ANY_SOURCE,10, node_comm,&status);
      i=status.MPI_SOURCE;
      obiviindx0=realloc(obiviindx0,(lj2+1)*sizeof(long int));
      obiviindx0[0]=-1;
      biviindx0=obiviindx0;
      biviindx0++;
      solve_real *vecbivi0= (solve_real *) calloc (lj2,sizeof(solve_real));
    sindx01=lj2/MAXSSIZE;
    sindx03=0;
    for(sindx02=0; sindx02<sindx01; sindx02++) {
      if(SORD==1)MPI_Recv(vecbivi0+sindx03,MAXSSIZE, MPI_DOUBLE,i,11*sindx01+sindx02, node_comm,&status);
      else MPI_Recv(vecbivi0+sindx03,MAXSSIZE, MPI_FLOAT,i,11*sindx01+sindx02, node_comm,&status);
      sindx03+=MAXSSIZE;
    }
      if(SORD==1)MPI_Recv(vecbivi0+sindx03,lj2-sindx01*MAXSSIZE, MPI_DOUBLE,i,11, node_comm,&status);
      else MPI_Recv(vecbivi0+sindx03,lj2-sindx01*MAXSSIZE, MPI_FLOAT,i,11, node_comm,&status);
    sindx03=0;
    for(sindx02=0; sindx02<sindx01; sindx02++) {
      MPI_Recv(biviindx0+sindx03,MAXSSIZE, MPI_LONG,i,12*sindx01+sindx02, node_comm,&status);
      sindx03+=MAXSSIZE;
    }
      MPI_Recv(biviindx0+sindx03,lj2-sindx01*MAXSSIZE, MPI_LONG,i,12, node_comm,&status);
      my_spar_compl_(biviindx1,&nz0,biviindx0,&lj2,&nz1);
      vecbivi=realloc(vecbivi,nz1*sizeof(solve_real));
      obiviindx1=realloc(obiviindx1,(nz1+1)*sizeof(long int));
      obiviindx1[0]=-1;
      biviindx1=obiviindx1;
      biviindx1++;
      my_spar_add3l_(vecbivi,obiviindx1,&nz0,vecbivi0,obiviindx0,&lj2,&nz1);
      free(obiviindx0);
      obiviindx0=NULL;
      free(vecbivi0);
      vecbivi0=NULL;
      logmsg(2,"rank %d nz0 %ld j2 %ld nz1 %ld\n",rank,nz0,lj2,nz1);
      nz0=nz1;
    }
  }

  if(group_size11>1&&color==1) {
    nz3=nz0;
    if(node_tail_rank!=group_size1-1) {
      lj2=nz3;
      MPI_Send(&lj2,1, MPI_LONG,group_size1-1, 10, node_tail_comm);
    sindx01=lj2/MAXSSIZE;
    sindx03=0;
    for(sindx02=0; sindx02<sindx01; sindx02++) {
      if(SORD==1)MPI_Send(vecbivi+sindx03,MAXSSIZE, MPI_DOUBLE,group_size1-1, 11*sindx01+sindx02, node_tail_comm);
      else MPI_Send(vecbivi+sindx03,MAXSSIZE, MPI_FLOAT,group_size1-1, 11*sindx01+sindx02, node_tail_comm);
      sindx03+=MAXSSIZE;
    }
      if(SORD==1)MPI_Send(vecbivi+sindx03,lj2-sindx01*MAXSSIZE, MPI_DOUBLE,group_size1-1, 11, node_tail_comm);
      else MPI_Send(vecbivi+sindx03,lj2-sindx01*MAXSSIZE, MPI_FLOAT,group_size1-1, 11, node_tail_comm);
    sindx03=0;
    for(sindx02=0; sindx02<sindx01; sindx02++) {
      MPI_Send(biviindx1+sindx03,MAXSSIZE, MPI_LONG,group_size1-1, 12*sindx01+sindx02, node_tail_comm);
      sindx03+=MAXSSIZE;
    }
      MPI_Send(biviindx1+sindx03,lj2-sindx01*MAXSSIZE, MPI_LONG,group_size1-1, 12, node_tail_comm);
      free(obiviindx1);
      obiviindx1=NULL;
      free(vecbivi);
      vecbivi=NULL;
    }
    if(node_tail_rank==group_size1-1) {
      for(j6=0; j6<group_size1-1; j6++) {
        MPI_Recv(&lj2,1, MPI_LONG,  MPI_ANY_SOURCE,10, node_tail_comm,&status);
        i=status.MPI_SOURCE;
        obiviindx0=realloc(obiviindx0,(lj2+1)*sizeof(long int));
        obiviindx0[0]=-1;
        biviindx0=obiviindx0;
        biviindx0++;
        solve_real *vecbivi0= (solve_real *) calloc (lj2,sizeof(solve_real));
    sindx01=lj2/MAXSSIZE;
    sindx03=0;
    for(sindx02=0; sindx02<sindx01; sindx02++) {
        if(SORD==1)MPI_Recv(vecbivi0+sindx03,MAXSSIZE, MPI_DOUBLE,i,11*sindx01+sindx02, node_tail_comm,&status);
        else MPI_Recv(vecbivi0+sindx03,MAXSSIZE, MPI_FLOAT,i,11*sindx01+sindx02, node_tail_comm,&status);
        sindx03+=MAXSSIZE;
    }
        if(SORD==1)MPI_Recv(vecbivi0+sindx03,lj2-sindx01*MAXSSIZE, MPI_DOUBLE,i,11, node_tail_comm,&status);
        else MPI_Recv(vecbivi0+sindx03,lj2-sindx01*MAXSSIZE, MPI_FLOAT,i,11, node_tail_comm,&status);
    sindx03=0;
    for(sindx02=0; sindx02<sindx01; sindx02++) {
        MPI_Recv(biviindx0+sindx03,MAXSSIZE, MPI_LONG,i,12*sindx01+sindx02, node_tail_comm,&status);
        sindx03+=MAXSSIZE;
    }
        MPI_Recv(biviindx0+sindx03,lj2-sindx01*MAXSSIZE, MPI_LONG,i,12, node_tail_comm,&status);
        my_spar_compl_(biviindx1,&nz0,biviindx0,&lj2,&nz1);
        vecbivi=realloc(vecbivi,nz1*sizeof(solve_real));
        obiviindx1=realloc(obiviindx1,(nz1+1)*sizeof(long int));
        obiviindx1[0]=-1;
        biviindx1=obiviindx1;
        biviindx1++;
        my_spar_add3l_(vecbivi,obiviindx1,&nz0,vecbivi0,obiviindx0,&lj2,&nz1);
        free(obiviindx0);
        obiviindx0=NULL;
        free(vecbivi0);
        vecbivi0=NULL;
        logmsg(2,"rank %d nz0 %ld j2 %ld nz1 %ld\n",rank,nz0,lj2,nz1);
        nz0=nz1;
      }
    }
  }
  fortran_int vecbiuisize=0;
  vecbiuisize=VecSize-sumrowcolin;
  reduce_to_rank_nocompress(vecbiui,vecbiuisize,mpisize,rank,mpisize-1);
  logmsg(2,"Completed interface reduce, rank %d time %f\n",rank,((double)(clock()-timestr))/CLOCKS_PER_SEC);
  if(rank!=mpisize-1){
    free(vecbiui);
    vecbiui=NULL;
  }
  timestr=clock();
  solve_real *xd;//= (ha_cgetype*)calloc(vecbiuisize,sizeof(ha_cgetype));
  solve_real *x0;//,*vecbivi0;//= (ha_cgetype*)calloc(VecSize,sizeof(ha_cgetype));
  long int lnz;
  if(rank==mpisize-1) {
    Mat_SeqAIJ         *aa=(Mat_SeqAIJ*)submatD[0]->data;//*aa=subA->data;
    ai= aa->i;
    aj= aa->j;
    vals=aa->a;
    lnz=aa->nz;
    nrow=submatD[0]->rmap->n;
    ncol=submatD[0]->cmap->n;
    logmsg(2,"D ncol %d nrow %d nz %ld\n",ncol,nrow,lnz);
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
    logmsg(2,"nz0 %ld nz %ld nz1 %ld\n",nz0,lnz,nz1);
    ldsize=ceil((laD/100.0)*nz1);
    vecbivi=realloc(vecbivi,ldsize*sizeof(solve_real));
    int *irn1=(int *) calloc (nz1,sizeof(int));
    int *jcn=(int *) calloc (ldsize,sizeof(int));
    logmsg(2,"nz1 %ld obin %ld ond1 %ld\n",nz1,obiviindx1[0],obiviindx0[0]);
    my_spar_add4l_(vecbivi,obiviindx1,irn1,jcn,&nz0,vals,obiviindx0,&lnz,&nz1,&ncol);
    logmsg(2,"nz1 %ld\n",nz1);
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
    xd=(solve_real *) calloc (vecbiuisize,sizeof(solve_real));//realloc (xd,vecbiuisize*sizeof(ha_cgetype));
    spec48m_ssol2la_(insizeD,irn1,jcn,vecbivi,vecbiui,xd);
    free(vecbivi);
    vecbivi=NULL;
    free(vecbiui);
    vecbiui=NULL;
    free(irn1);
    free(jcn);
    free(insizeD);
    x0=(solve_real *) calloc (VecSize,sizeof(solve_real));//realloc (x0,*sizeof(ha_cgetype));
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
    xd=(solve_real *) calloc (vecbiuisize,sizeof(solve_real));//realloc (xd,vecbiuisize*sizeof(ha_cgetype));
    x0=(solve_real *) calloc (VecSize,sizeof(solve_real));//realloc (x0,*sizeof(ha_cgetype));
  }
  PetscFree(submatD);
  MPI_Barrier(PETSC_COMM_WORLD);
  if(SORD==1)MPI_Bcast(xd, vecbiuisize, MPI_DOUBLE,mpisize-1, PETSC_COMM_WORLD);
  else MPI_Bcast(xd, vecbiuisize, MPI_FLOAT,mpisize-1, PETSC_COMM_WORLD);
  logmsg(2,"Interface problem solution rank %d time %f\n",rank,((double)(clock()-timestr))/CLOCKS_PER_SEC);
  timestr=clock();
  xi1indx=0;

  omp_set_num_threads(section_threads);
  #pragma omp parallel private(jthrd,j1,j2,j,i,timestr,ai,aj,vals,nz,nrow,xi1point,xi1indx,ierr,nindices) shared(submatC,submatCij,submatBij,insize,xi1,colindicesbc1,colindicesbc2,x0)
  {
  jthrd=omp_get_thread_num();
  solve_real *biui0= (solve_real *) calloc (maxrowc,sizeof(solve_real));
  solve_real *be0 = (solve_real*)calloc(maxrowc,sizeof(solve_real));
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
    if(inmemory||nfr_flag()) {
      /* back-substitute from the resident store, then drop this
         window's factors (the ifremove equivalent) */
      int icntl2[20],info2[20],mrc=0;
      solve_real cntl2[10],rinfo2[10],error12[3];
      for(i=0; i<nreg; i++)if((submatCij[j2+i][0]->rmap->n)>mrc)mrc=submatCij[j2+i][0]->rmap->n;
      solve_real *b02=(solve_real*)calloc(mrc,sizeof(solve_real));
      solve_real *w2=(solve_real*)calloc(4*maxcolc,sizeof(solve_real));
      int *iw2=(int*)calloc(maxcolc,sizeof(int));
      ndbbd_block_solve_mem(rank,j2,nreg,insize,insizes,submatCij,submatBij,be0,biui0,ndbbd_fac_irn+j2,ndbbd_fac_keep+j2,ndbbd_fac_va+j2,cntl2,rinfo2,error12,icntl2,info2,w2,iw2,b02);
      free(b02);
      free(w2);
      free(iw2);
      for(i=0; i<nreg+1; i++) {
        if(nfr_flag()&&i!=nreg)continue;/* regional factors persist across steps */
        ndbbd_fac_drop(j2+i);
      }
    }
    else ndbbd_block_solve(rank,j2,nreg,insize,insizes,submatCij,submatBij,be0,biui0,ifremove,fn01,fn02,fn03);
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
    logmsg(2,"Submatrix %d rank %d calculation time %f\n",j1,rank,difftime(timeend,timestr)/CLOCKS_PER_SEC);
  }
  free(be0);
  free(biui0);
  }
  omp_set_num_threads(max_threads);
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
  logmsg(2,"Solution calculation rank %d time %f\n",rank,((double)(clock()-timestr))/CLOCKS_PER_SEC);
  MPI_Barrier(PETSC_COMM_WORLD);
  free(xd);
  timestr=clock();
  if(SORD==1)MPI_Allreduce(x0,x1,VecSize, MPI_DOUBLE, MPI_SUM,PETSC_COMM_WORLD);
  else MPI_Allreduce(x0,x1,VecSize, MPI_FLOAT, MPI_SUM,PETSC_COMM_WORLD);
  logmsg(2,"Reduce solution rank %d time %f\n",rank,((double)(clock()-timestr))/CLOCKS_PER_SEC);
  free(x0);
  return 0;
}




bool ndbbd_block_solve(PetscInt rank, int begmat,int nreg,int * insize,int insizes, Mat **submatCij,Mat **submatBij,solve_real *b,solve_real *sol,bool ifremove,char** fn01,char** fn02, char** fn03) {
  FILE* fp1,*fp2,*fp3;
  PetscScalar *vals,*valsc,vecval;
  PetscInt *ai,*aj,*aic,*ajc,nrow,nz,maxrowcij;//,ncol
  PetscInt i,j,j1,j2,indx01,la1;
  solve_real *b01,*b03,*sol1,*sol2;
  size_t freadresult;
  int *irne= NULL;//(int*)calloc(1,sizeof(int));
  int *keep= NULL;//(int*)calloc(1,sizeof(int));
  solve_real *vale = NULL;//(ha_cgetype*)calloc(1,sizeof(ha_cgetype));
  maxrowcij=0;
  b01=b;
  sol1=sol;
  int fd1,fd2,fd3,frrsl1,frrsl2,frrsl3;
  for (j1=begmat; j1<nreg+begmat; j1++) {
    j2=j1*insizes;//-begmat;
    if((submatCij[j1][0]->rmap->n)>maxrowcij)maxrowcij=submatCij[j1][0]->rmap->n;
    #pragma omp critical(nsol)
    {
    la1=ceil((insize[j2+9]/100.0)*insize[j2+2]);
    irne = realloc(irne,la1*sizeof(int));//(int*)calloc(la1,sizeof(int));
    keep = realloc(keep,insize[j2+12]*sizeof(int));//(int*)calloc(insize[j1*insizes+12],sizeof(int));
    vale = realloc(vale,la1*sizeof(solve_real));//(ha_cgetype*)calloc(la1,sizeof(ha_cgetype));
    fp1 = fopen(fn01[j1], "rb");
    if (fp1==NULL)printf("Error: cannot open factor file %s\n",fn01[j1]);
    freadresult=fread(irne,sizeof(int),la1,fp1);
    fclose(fp1);
    fp2 = fopen(fn02[j1], "rb");
    if (fp2==NULL)printf("Error: cannot open factor file %s\n",fn02[j1]);
    freadresult=fread(keep,sizeof(int),insize[j2+12],fp2);
    fclose(fp2);
    fp3 = fopen(fn03[j1], "rb");
    if (fp3==NULL)printf("Error: cannot open factor file %s\n",fn03[j1]);
    freadresult=fread(vale,sizeof(solve_real),la1,fp3);
    fclose(fp3);
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
  la1=ceil((insize[j1*insizes+9]/100.0)*insize[j1*insizes+2]);
  irne = realloc(irne,la1*sizeof(int));//(int*)calloc(la1,sizeof(int));
  keep = realloc(keep,insize[j1*insizes+12]*sizeof(int));//(int*)calloc(insize[j1*insizes+12],sizeof(int));
  vale = realloc(vale,la1*sizeof(solve_real));//(ha_cgetype*)calloc(la1,sizeof(ha_cgetype));
  fp1 = fopen(fn01[j2], "rb");
  if (fp1==NULL)printf("Error: cannot open factor file %s\n",fn01[j2]);
  freadresult=fread(irne,sizeof(int),la1,fp1);
  fclose(fp1);
  fp2 = fopen(fn02[j2], "rb");
  if (fp2==NULL)printf("Error: cannot open factor file %s\n",fn02[j2]);
  freadresult=fread(keep,sizeof(int),insize[j1*insizes+12],fp2);
  fclose(fp2);
  fp3 = fopen(fn03[j2], "rb");
  if (fp3==NULL)printf("Error: cannot open factor file %s\n",fn03[j2]);
  freadresult=fread(vale,sizeof(solve_real),la1,fp3);
  fclose(fp3);
    }
  if(ifremove)remove(fn01[j2]);
  if(ifremove)remove(fn02[j2]);
  if(ifremove)remove(fn03[j2]);
  if(insize[j1*insizes+16]!=la1)insize[j1*insizes+16]=la1;
  spec48m_esol_(insize+j1*insizes,irne,vale,keep,b01,sol1);
  sol2=sol;
  b03=b;
  solve_real *b02 = (solve_real*)calloc(maxrowcij,sizeof(solve_real));
  for (j1=begmat; j1<nreg+begmat; j1++) {
    j2=j1;//-begmat;
    #pragma omp critical(nsol)
    {
    la1=ceil((insize[j1*insizes+9]/100.0)*insize[j1*insizes+2]);
    irne = realloc(irne,la1*sizeof(int));//(int*)calloc(la1,sizeof(int));
    keep = realloc(keep,insize[j1*insizes+12]*sizeof(int));//(int*)calloc(insize[j1*insizes+12],sizeof(int));
    vale = realloc(vale,la1*sizeof(solve_real));//(ha_cgetype*)calloc(la1,sizeof(ha_cgetype));
    fp1 = fopen(fn01[j2], "rb");
    if (fp1==NULL)printf("Error: cannot open factor file %s\n",fn01[j2]);
    freadresult=fread(irne,sizeof(int),la1,fp1);
    fclose(fp1);
    fp2 = fopen(fn02[j2], "rb");
    if (fp2==NULL)printf("Error: cannot open factor file %s\n",fn02[j2]);
    freadresult=fread(keep,sizeof(int),insize[j1*insizes+12],fp2);
    fclose(fp2);
    fp3 = fopen(fn03[j2], "rb");
    if (fp3==NULL)printf("Error: cannot open factor file %s\n",fn03[j2]);
    freadresult=fread(vale,sizeof(solve_real),la1,fp3);
    fclose(fp3);
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


bool ndbbd_block_solve_mem(PetscInt rank, int begmat,int nreg,int * insize,int insizes, Mat **submatCij,Mat **submatBij,solve_real *b,solve_real *sol,int** irnereg,int** keepreg,solve_real** valereg,solve_real *cntl,solve_real *rinfo,solve_real *error1,int *icntl,int *info,solve_real *w,int *iw,solve_real *b02) {
  PetscScalar *vals;//,*valsc;//,vecval;
  PetscInt *ai,*aj,nrow,nz;//,ncol,*aic,*ajc
  PetscInt i,j,j1,j2;//,indx01,la1;
  solve_real *b01,*b03,*sol1,*sol2;//,*b02
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

int reduce_to_rank(solve_real *vecbivi,fortran_int vecbivisize,PetscInt mpisize,PetscInt rank,PetscInt targetrank) {
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
  solve_real *vecbivi0= (solve_real *) calloc (1,sizeof(solve_real));
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
    vecbivi0=realloc(vecbivi0,j2*sizeof(solve_real));
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

int reduce_to_rank_nocompress(solve_real *vecbivi,fortran_int vecbivisize,PetscInt mpisize,PetscInt rank,PetscInt targetrank) {
  if(mpisize==1)return 0;
  MPI_Status   status;
  int i,j,j1,j2,j3;
  solve_real *vecbivi0= (solve_real *) calloc (1,sizeof(solve_real));

  if(rank==targetrank) {
    vecbivi0=realloc(vecbivi0,vecbivisize*sizeof(solve_real));
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




