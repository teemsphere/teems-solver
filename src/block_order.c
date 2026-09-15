#include <teems_solver.h>
#include <hsl_kernels.h>

int dbbd_order(Mat A, offset_t VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend, offset_t nvarele, PetscInt *eq_addr,int *row_order,int *col_order, offset_t ndblock,int *block_sizes, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,solve_real cntl6) {
  IS *rowindices,*colindices;//,isrow,iscol;
  PetscInt bfirst,bend,nmatin,nrowcolin;
  Mat *submatA;
  PetscInt i,j,j1,j2,nrow,ncol,nz,nz1,*ai,*aj;
  PetscScalar *vals;
  PetscErrorCode ierr;
  offset_t lasize;
  int *row_order1= (int *) calloc (VecSize,sizeof(int));
  int *col_order1= (int *) calloc (VecSize,sizeof(int));
  int *block_sizes1= (int *) calloc (ndblock,sizeof(int));

  //*****New order test**********//
  nmatin=(offset_t)ndblock/mpisize;
  for(i=0; i<mpisize; i++)if(rank+1<=ndblock-mpisize*nmatin)nmatin++;
  int *begblock= (int *) calloc (mpisize,sizeof(int));
  ierr = PetscMalloc(nmatin*sizeof(IS **),&rowindices);
  CHKERRQ(ierr);
  ierr = PetscMalloc(nmatin*sizeof(IS **),&colindices);
  CHKERRQ(ierr);
  logmsg(2,"rank %d nmatin %d\n",rank,nmatin);
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
    logmsg(2,"rank %d nz %d\n",rank,nz);
    nz1=nz;
    if(nz1<nrow)nz1=nrow;
    if(nz1<ncol)nz1=ncol;
    lasize=ma48_la_from_pct(laA,nz1);
    /* a starved -laA (<100) must still stage all NE entries; MA48
       then returns -3 and the growth loop below takes over */
    if(lasize<nz1)lasize=nz1;
    int *irn=(int *) calloc (lasize,sizeof(int));
    int *jcn=(int *) calloc (lasize,sizeof(int));
    int *irn1=(int *)malloc(nrow*sizeof(int));
    int *jcn1=(int *)malloc(ncol*sizeof(int));
    solve_real *values= (solve_real *) calloc (lasize,sizeof(solve_real));
    int *keep=(int *) malloc ((nrow+5*ncol+4*ncol+7)*sizeof(int));//ICNTL(6)=1
  solve_real *w51;
  if(nrow>ncol) w51=(solve_real *) malloc ((5*nrow)*sizeof(solve_real));
  else w51=(solve_real *) malloc ((5*ncol)*sizeof(solve_real));
  int *iw51=(int *) malloc ((6*nrow+3*ncol)*sizeof(int));
    int *insize=(int *) calloc (8,sizeof(int));
    int tries;
    for(tries=0; tries<6; tries++) {
    insize[0]=nrow;
    insize[1]=ncol;
    insize[2]=nz;
    insize[4]=laA;
    insize[5]=lasize;
    insize[6]=0;
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
    if(insize[6]!=-3)break;
    {
      /* MA48 workspace too small for the rank probe: grow and
         re-stage (MA48 clobbers the staged triplet in place) */
      offset_t newla=ma48_grow_la(lasize,insize[7],nz1,"laA",&teems_laA_used);
      lasize=newla;
    }
    irn=ma48_realloc(irn,lasize,sizeof(int));
    jcn=ma48_realloc(jcn,lasize,sizeof(int));
    values=ma48_realloc(values,lasize,sizeof(solve_real));
    }
    if(insize[6]==-3) {
      printf("Error: the MA48 workspace for %s did not converge after %d growth attempts; raise the initial workspace (laA/laD/laDi) or use a bordered matrix_method (\"SBBD\" or \"DBBD\")\n",probe_onfail_scope_label(),tries);
      MPI_Abort(PETSC_COMM_WORLD,1);
    }
    for(i=0; i<nrow; i++) {
      row_order1[i+counteq[j1+begblock[rank]]]=irn1[i]-1;
    }
    for(i=0; i<ncol; i++) {
      col_order1[i+countvarintra1[j1+begblock[rank]]]=jcn1[i]-1;
    }
    logmsg(2,"rank %d j1 %d proc %d\n",insize[3],j1,rank);
    block_sizes1[j1+begblock[rank]]=insize[3];
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
  MPI_Allreduce(col_order1,col_order,VecSize, MPI_INT, MPI_SUM,PETSC_COMM_WORLD);
  free(col_order1);
  MPI_Allreduce(row_order1,row_order,VecSize, MPI_INT, MPI_SUM,PETSC_COMM_WORLD);
  free(row_order1);
  MPI_Allreduce(block_sizes1,block_sizes,ndblock, MPI_INT, MPI_SUM,PETSC_COMM_WORLD);
  free(block_sizes1);
  free(begblock);
  MPI_Barrier(PETSC_COMM_WORLD);
  return 1;
}



/* ---- NDBBD cut cache (see teems_solver.h) ---- */
static int ndcut_valid=0;
static offset_t ndcut_VecSize=-1,ndcut_ndblock=-1;
static int ndcut_ntime=-1;
static int *ndcut_row=NULL,*ndcut_col=NULL,*ndcut_bs=NULL;
static offset_t *ndcut_ceq=NULL,*ndcut_ceqna=NULL,*ndcut_cvi=NULL;
static PetscInt *ndcut_rank=NULL;
static int ndcut_niface=0;
static int *ndcut_if_rank=NULL,*ndcut_if_nrow=NULL,*ndcut_if_ncol=NULL,*ndcut_if_set=NULL;
static int **ndcut_if_irn=NULL,**ndcut_if_jcn=NULL;

void ndbbd_cut_cache_free(void) {
  int i;
  free(ndcut_row); free(ndcut_col); free(ndcut_bs);
  free(ndcut_ceq); free(ndcut_ceqna); free(ndcut_cvi); free(ndcut_rank);
  ndcut_row=ndcut_col=ndcut_bs=NULL;
  ndcut_ceq=ndcut_ceqna=ndcut_cvi=NULL;
  ndcut_rank=NULL;
  for(i=0; i<ndcut_niface; i++) {
    free(ndcut_if_irn[i]);
    free(ndcut_if_jcn[i]);
  }
  free(ndcut_if_irn); free(ndcut_if_jcn); free(ndcut_if_rank); free(ndcut_if_nrow); free(ndcut_if_ncol); free(ndcut_if_set);
  ndcut_if_irn=ndcut_if_jcn=NULL;
  ndcut_if_rank=ndcut_if_nrow=ndcut_if_ncol=ndcut_if_set=NULL;
  ndcut_niface=0;
  ndcut_valid=0;
  ndcut_VecSize=ndcut_ndblock=-1;
  ndcut_ntime=-1;
}

/* per-rank interface cut slots, one per local chain block; call before
   the presolve's interface loop (outside its parallel region) */
void ndbbd_cut_iface_init(int nmatint) {
  int i;
  if(ndcut_niface==nmatint)return;
  for(i=0; i<ndcut_niface; i++) {
    free(ndcut_if_irn[i]);
    free(ndcut_if_jcn[i]);
  }
  free(ndcut_if_irn); free(ndcut_if_jcn); free(ndcut_if_rank); free(ndcut_if_nrow); free(ndcut_if_ncol); free(ndcut_if_set);
  ndcut_niface=nmatint;
  ndcut_if_irn=(int **) calloc (nmatint,sizeof(int *));
  ndcut_if_jcn=(int **) calloc (nmatint,sizeof(int *));
  ndcut_if_rank=(int *) calloc (nmatint,sizeof(int));
  ndcut_if_nrow=(int *) calloc (nmatint,sizeof(int));
  ndcut_if_ncol=(int *) calloc (nmatint,sizeof(int));
  ndcut_if_set=(int *) calloc (nmatint,sizeof(int));
}

/* 1 = slot j3 holds a cut for a block of this shape (rank + permutations copied out) */
int ndbbd_cut_iface_get(int j3,int *rank_out,int *irn,int *jcn,int nrow,int ncol) {
  if(!teems_ndcutcache||teems_ndcutcache==2||j3<0||j3>=ndcut_niface||!ndcut_if_set[j3])return 0; /* 2 = regional cache only (bisect aid) */
  if(ndcut_if_nrow[j3]!=nrow||ndcut_if_ncol[j3]!=ncol)return 0;
  *rank_out=ndcut_if_rank[j3];
  memcpy(irn,ndcut_if_irn[j3],nrow*sizeof(int));
  memcpy(jcn,ndcut_if_jcn[j3],ncol*sizeof(int));
  return 1;
}

void ndbbd_cut_iface_put(int j3,int rank_val,const int *irn,const int *jcn,int nrow,int ncol) {
  if(!teems_ndcutcache||j3<0||j3>=ndcut_niface)return;
  free(ndcut_if_irn[j3]);
  free(ndcut_if_jcn[j3]);
  ndcut_if_irn[j3]=(int *) malloc ((nrow>0?nrow:1)*sizeof(int));
  ndcut_if_jcn[j3]=(int *) malloc ((ncol>0?ncol:1)*sizeof(int));
  memcpy(ndcut_if_irn[j3],irn,nrow*sizeof(int));
  memcpy(ndcut_if_jcn[j3],jcn,ncol*sizeof(int));
  ndcut_if_rank[j3]=rank_val;
  ndcut_if_nrow[j3]=nrow;
  ndcut_if_ncol[j3]=ncol;
  ndcut_if_set[j3]=1;
}

int ndbbd_order_presolve(Mat A, offset_t VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int nreg, int ntime, offset_t nvarele, PetscInt *eq_addr,int *row_order,int *col_order, offset_t ndblock,int *block_sizes, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,dim_t laDi,solve_real cntl6,PetscInt* ndbbdrank,PetscBool presol) {
  char filename[1024],rankname[1024];
  IS *rowindices=NULL,*colindices=NULL;//,isrow,iscol;
  PetscInt bfirst,bend,nmatin,nmatint;
  Mat *submatA=NULL;//,*submatD=NULL;
  PetscInt i,j,j1,j2,j3,j4,j5,j6,j7,j8,nrow,ncol,nz,nz1,*ai,*aj;
  offset_t lasize;
  solve_real cntl6in;
  PetscScalar *vals;
  PetscErrorCode ierr;
  MatInfo           matinfo;
  if(teems_ndcutcache&&teems_ndcutcache!=3&&ndcut_valid&&ndcut_VecSize==VecSize&&ndcut_ndblock==ndblock&&ndcut_ntime==ntime) { /* 3 = interface cache only (bisect aid) */
    /* the cut is structural: reuse step 1's instead of re-probing */
    memcpy(counteq,ndcut_ceq,(ndblock+1)*sizeof(offset_t));
    memcpy(counteqnoadd,ndcut_ceqna,ndblock*sizeof(offset_t));
    memcpy(countvarintra1,ndcut_cvi,(ndblock+1)*sizeof(offset_t));
    memcpy(row_order,ndcut_row,VecSize*sizeof(int));
    memcpy(col_order,ndcut_col,VecSize*sizeof(int));
    memcpy(block_sizes,ndcut_bs,ndblock*sizeof(int));
    memcpy(ndbbdrank,ndcut_rank,ntime*sizeof(PetscInt));
    return 1;
  }
  MatGetInfo(A,MAT_LOCAL,&matinfo);
  logmsg(2,"rank %d matinfo.nz_used %g\n",rank,matinfo.nz_used);
  int *block_sizes1= (int *) calloc (ndblock,sizeof(int));
  offset_t *counteq2= (offset_t *) calloc (ndblock+1,sizeof(offset_t));
  offset_t *counteqnoadd1= (offset_t *) calloc (ndblock,sizeof(offset_t));
  offset_t *counteqnoadd2= (offset_t *) calloc (ndblock,sizeof(offset_t));
  offset_t *countvarintra2= (offset_t *) calloc (ndblock+1,sizeof(offset_t));
  offset_t *countvarintra4= (offset_t *) calloc (ndblock+1,sizeof(offset_t));
  offset_t *countvarintra6= (offset_t *) calloc (ndblock+1,sizeof(offset_t));
  if(rank<10)strcpy(rankname,"000");
  if(rank<100&&rank>9)strcpy(rankname,"00");
  if(rank<1000&&rank>99)strcpy(rankname,"0");
  if(rank>=1000)rankname[0]='\0';
  sprintf(filename, "%d",rank);
  strcat(rankname,filename);
  //*****New order test**********//
  nmatint=(offset_t)ntime/mpisize;
  for(i=0; i<mpisize; i++)if(rank+1<=ntime-mpisize*nmatint)nmatint++;
  nmatin=(nreg+1)*nmatint;
  int *begblock= (int *) calloc (mpisize,sizeof(int));
  begblock[rank]=nmatin;
  logmsg(2,"rank %d nmatin %d nmatint %d ndblock %ld\n",rank,nmatin,nmatint,ndblock);
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
  int *row_order2= (int *) calloc (VecSize,sizeof(int));
  int *col_order2= (int *) calloc (VecSize,sizeof(int));
  if(cntl6==0&&SORD==0)cntl6in=1e-5;
  else cntl6in=cntl6;
  #pragma omp parallel private(j1,j2,j3,j4,j5,j6,j7,j8,ai,aj,vals,nz,nrow,ncol,nz1,i,j,bfirst,lasize) shared(cntl6in,row_order2,col_order2,block_sizes1,counteqnoadd1,countvarintra2,countvarintra1,submatA)
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
      lasize=ma48_la_from_pct(laA,nz1);
      /* starved -laA: staging still needs all NE entries per block */
      if(lasize<nz1)lasize=nz1;
      if(lasizemax<lasize)lasizemax=lasize;
      j8++;
    }
  }
  lasizemax+=10;
    long int alloc_la=lasizemax; /* current array length; grows on MA48 -3 */
    int *irn=(int *) calloc (lasizemax,sizeof(int));
    int *jcn=(int *) calloc (lasizemax,sizeof(int));
    int*irn1=NULL,*jcn1=NULL;
  irn1=realloc(irn1,nrowmax*sizeof(int));
  jcn1=realloc(jcn1,ncolmax*sizeof(int));
    solve_real *values= (solve_real *) calloc (lasizemax,sizeof(solve_real));
  int *keep=(int *) malloc ((nrowmax+5*ncolmax+4*ncolmax+7)*sizeof(int));//ICNTL(6)=1
  solve_real *w51;
  if(nrowmax>ncolmax) w51=(solve_real *) malloc ((5*nrowmax)*sizeof(solve_real));
  else w51=(solve_real *) malloc ((5*ncolmax)*sizeof(solve_real));
  int *iw51=(int *) malloc ((6*nrowmax+3*ncolmax)*sizeof(int));
    int *insize=(int *) calloc (8,sizeof(int));
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
      lasize=ma48_la_from_pct(laA,nz1);
      if(lasize<nz1)lasize=nz1;
      int tries;
      for(tries=0; tries<6; tries++) {
      insize[0]=nrow;
      insize[1]=ncol;
      insize[2]=nz;
      insize[4]=laA;
      insize[5]=lasize;
      insize[6]=0;
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
      if(insize[6]!=-3)break;
      {
        /* MA48 workspace too small for the rank probe: grow and
           re-stage (MA48 clobbers the staged triplet in place) */
        offset_t newla=ma48_grow_la(lasize,insize[7],nz1,"laA",&teems_laA_used);
        lasize=newla;
      }
      if(lasize>alloc_la) {
        alloc_la=lasize;
        irn=ma48_realloc(irn,alloc_la,sizeof(int));
        jcn=ma48_realloc(jcn,alloc_la,sizeof(int));
        values=ma48_realloc(values,alloc_la,sizeof(solve_real));
      }
      }
      if(insize[6]==-3) {
        printf("Error: the MA48 workspace for %s did not converge after %d growth attempts; raise the initial workspace (laA/laD/laDi) or use a bordered matrix_method (\"SBBD\" or \"DBBD\")\n",probe_onfail_scope_label(),tries);
        MPI_Abort(PETSC_COMM_WORLD,1);
      }
      bfirst=counteq[j3+begblock[rank]];
      for(i=0; i<nrow; i++) {
        row_order2[i+counteq[j3+begblock[rank]]]=irn1[i]-1+bfirst;
      }
      bfirst=countvarintra1[j3+begblock[rank]];
      for(i=0; i<ncol; i++) {
        col_order2[i+countvarintra1[j3+begblock[rank]]]=jcn1[i]-1+bfirst;
      }
      if(verbosity>=2) {
        unsigned long ph=1469598103934665603UL;
        for(i=0; i<nrow; i++)ph=(ph^(unsigned long)irn1[i])*1099511628211UL;
        for(i=0; i<ncol; i++)ph=(ph^(unsigned long)jcn1[i])*1099511628211UL;
        logmsg(2,"rank %d mat rank %d nrow %d ncol %d j3 %d proc %d permsum %lx\n",rank,insize[3],nrow,ncol,j3,rank,ph);
      }
      block_sizes1[j3+begblock[rank]]=insize[3];
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

  int *row_order1= (int *) calloc (VecSize,sizeof(int));
  int *col_order1= (int *) calloc (VecSize,sizeof(int));
  for(j1=0; j1<nmatint; j1++) {
    for(j2=0; j2<nreg; j2++) {
      j3=j1*(nreg+1)+j2;
      bfirst=counteq2[j3+begblock[rank]];
      for(i=bfirst; i<counteq2[j3+begblock[rank]]+counteqnoadd2[j3+begblock[rank]]; i++) {
        row_order1[i]=row_order2[i-bfirst+counteq[j3+begblock[rank]]];
      }
      bfirst=countvarintra4[j3+begblock[rank]];
      for(i=bfirst; i<countvarintra4[j3+1+begblock[rank]]; i++) {
        col_order1[i]=col_order2[i-bfirst+countvarintra1[j3+begblock[rank]]];
      }
    }
    j3=j1*(nreg+1)+j2;
    j4=counteq2[j3+begblock[rank]];
    j5=countvarintra4[j3+begblock[rank]];
    for(j2=0; j2<nreg; j2++) {
      j6=j1*(nreg+1)+j2;
      for(i=counteqnoadd[j6+begblock[rank]]+counteq[j6+begblock[rank]]-1; i>counteqnoadd2[j6+begblock[rank]]+counteq[j6+begblock[rank]]-1; i--) {
        row_order1[j4]=row_order2[i];
        j4++;
      }
      for(i=countvarintra1[j6+1+begblock[rank]]-1; i>countvarintra1[j6+begblock[rank]]+countvarintra6[j6+begblock[rank]]-1; i--) {
        col_order1[j5]=col_order2[i];
        j5++;
      }
    }
    for(i=j4; i<j4+counteqnoadd[j3+begblock[rank]]; i++) {
      row_order1[i]=i-j4+counteq[j3+begblock[rank]];
    }
    for(i=j5; i<j5+countvarintra1[j3+1+begblock[rank]]-countvarintra1[j3+begblock[rank]]; i++) {
      col_order1[i]=i-j5+countvarintra1[j3+begblock[rank]];
    }
  }
  free(col_order2);
  free(row_order2);
  for(j1=0; j1<nmatint; j1++) {
    j3=j1*(nreg+1)+nreg;
    nrow=counteqnoadd2[j3+begblock[rank]];
    ncol=countvarintra6[j3+begblock[rank]];
    j4=j1+begblock[rank]/(nreg+1);

    if(nrow<ncol)ndbbdrank[j4]=nrow;
    else ndbbdrank[j4]=ncol;

    if(ndbbdrank[j4]>0) {
      block_sizes1[j3+begblock[rank]]=ndbbdrank[j4];
    }

  }
  memcpy(counteq,counteq2,(ndblock+1)*sizeof(offset_t));
  memcpy(counteqnoadd,counteqnoadd2,(ndblock)*sizeof(offset_t));
  memcpy(countvarintra1,countvarintra4,(ndblock+1)*sizeof(offset_t));
  MPI_Barrier(PETSC_COMM_WORLD);
  MPI_Allreduce(col_order1,col_order,VecSize, MPI_INT, MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(row_order1,row_order,VecSize, MPI_INT, MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(block_sizes1,block_sizes,ndblock, MPI_INT, MPI_SUM,PETSC_COMM_WORLD);
  free(col_order1);
  free(row_order1);
  free(block_sizes1);
  free(begblock);
  free(counteq2);
  free(counteqnoadd1);
  free(counteqnoadd2);
  free(countvarintra2);
  free(countvarintra4);
  free(countvarintra6);
  if(teems_ndcutcache) {
    /* store the cut for the following steps */
    if(ndcut_VecSize!=VecSize||ndcut_ndblock!=ndblock||ndcut_ntime!=ntime) {
      free(ndcut_row); free(ndcut_col); free(ndcut_bs); free(ndcut_ceq); free(ndcut_ceqna); free(ndcut_cvi); free(ndcut_rank);
      ndcut_row=(int *) malloc (VecSize*sizeof(int));
      ndcut_col=(int *) malloc (VecSize*sizeof(int));
      ndcut_bs=(int *) malloc (ndblock*sizeof(int));
      ndcut_ceq=(offset_t *) malloc ((ndblock+1)*sizeof(offset_t));
      ndcut_ceqna=(offset_t *) malloc (ndblock*sizeof(offset_t));
      ndcut_cvi=(offset_t *) malloc ((ndblock+1)*sizeof(offset_t));
      ndcut_rank=(PetscInt *) malloc (ntime*sizeof(PetscInt));
      ndcut_VecSize=VecSize;
      ndcut_ndblock=ndblock;
      ndcut_ntime=ntime;
    }
    memcpy(ndcut_ceq,counteq,(ndblock+1)*sizeof(offset_t));
    memcpy(ndcut_ceqna,counteqnoadd,ndblock*sizeof(offset_t));
    memcpy(ndcut_cvi,countvarintra1,(ndblock+1)*sizeof(offset_t));
    memcpy(ndcut_row,row_order,VecSize*sizeof(int));
    memcpy(ndcut_col,col_order,VecSize*sizeof(int));
    memcpy(ndcut_bs,block_sizes,ndblock*sizeof(int));
    memcpy(ndcut_rank,ndbbdrank,ntime*sizeof(PetscInt));
    ndcut_valid=1;
  }
  return 1;
}

int ndbbd_order(Mat A, offset_t VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int nreg, int ntime, offset_t nvarele, PetscInt *eq_addr,int *row_order,int *col_order, offset_t ndblock,int *block_sizes, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,dim_t laDi,solve_real cntl6,PetscInt* ndbbdrank,PetscBool presol) {
  FILE *presolfile;
  char j1name[1024],filename[1024],rankname[1024];
  size_t frd;
  PetscInt bfirst,bend,nmatin,nmatint;
  PetscInt i,j,j1,j2,j3,j4;
  if(rank<10)strcpy(rankname,"000");
  if(rank<100&&rank>9)strcpy(rankname,"00");
  if(rank<1000&&rank>99)strcpy(rankname,"0");
  if(rank>=1000)rankname[0]='\0';
  sprintf(filename, "%d",rank);
  strcat(rankname,filename);
  //*****New order test**********//
  nmatint=(offset_t)ntime/mpisize;
  for(i=0; i<mpisize; i++)if(rank+1<=ntime-mpisize*nmatint)nmatint++;
  nmatin=(nreg+1)*nmatint;
  int *begblock= (int *) calloc (mpisize,sizeof(int));
  begblock[rank]=nmatin;
  logmsg(2,"rank %d nmatin %d nmatint %d ndblock %ld\n",rank,nmatin,nmatint,ndblock);
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

  int *row_order3= (int *) calloc (VecSize,sizeof(int));
  int *col_order3= (int *) calloc (VecSize,sizeof(int));
  for(i=0; i<nmatin; i++) {
    bfirst=counteq[i+begblock[rank]];
    bend=counteqnoadd[i+begblock[rank]]+counteq[i+begblock[rank]];
    for(j=bfirst; j<bend; j++)row_order3[j]=row_order[j];
    bfirst=countvarintra1[i+begblock[rank]];
    bend=countvarintra1[i+1+begblock[rank]];
    for(j=bfirst; j<bend; j++)col_order3[j]=col_order[j];
  }
  int *block_sizes2= (int *) calloc (ndblock,sizeof(int));
  for(j1=0; j1<nmatin; j1++)block_sizes2[j1+begblock[rank]]=block_sizes[j1+begblock[rank]];

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
      strcpy(filename,scratch_dir);strcat(filename,"_rank");
      strcat(filename,rankname);
      strcat(filename,j1name);
      strcat(filename,".bin");
      if((presolfile=fopen(filename, "r"))==NULL) {
        printf("Error: cannot open interface file %s; the presolve pass did not produce it (check scratch space)\n",filename);
      }
      frd=fread(insized, sizeof(int), 5, presolfile);
      fclose(presolfile);
      if( frd== 0)printf("Error: short read on interface file %s; falling back to the minimum rank assignment\n",filename);
      if(ndbbdrank[j4]>insized[3]) {
        int *irn1=(int *) calloc (insized[0],sizeof(int));
        int *jcn1=(int *) calloc (insized[1],sizeof(int));
        strcpy(filename,scratch_dir);strcat(filename,"_row");
        strcat(filename,rankname);
        strcat(filename,j1name);
        strcat(filename,".bin");
        if((presolfile=fopen(filename, "r"))==NULL) {
          printf("Error: cannot open interface file %s\n",filename);
        }
        frd=fread(irn1, sizeof(int), insized[0], presolfile);
        if(frd == 0)printf("Error: short read on interface file %s\n",filename);
        fclose(presolfile);
        strcpy(filename,scratch_dir);strcat(filename,"_col");
        strcat(filename,rankname);
        strcat(filename,j1name);
        strcat(filename,".bin");
        if((presolfile=fopen(filename, "r"))==NULL) {
          printf("Error: cannot open interface file %s\n",filename);
        }
        frd=fread(jcn1, sizeof(int), insized[0], presolfile);
        if( frd== 0) printf("Error: short read on interface file %s\n",filename);
        fclose(presolfile);
        int *indices= (int *) calloc (insized[0],sizeof(int));
        for(i=0; i<insized[0]; i++) {
          indices[i]=row_order3[counteq[j3+begblock[rank]]+irn1[i]-1];
        }
        memcpy(&row_order3[counteq[j3+begblock[rank]]],indices,insized[0]*sizeof(int));
        for(i=0; i<insized[0]; i++) {
          indices[i]=col_order3[countvarintra1[j3+begblock[rank]]+jcn1[i]-1];
        }
        memcpy(&col_order3[countvarintra1[j3+begblock[rank]]],indices,insized[0]*sizeof(int));
        ndbbdrank[j4]=insized[3];
        free(indices);
        free(irn1);
        free(jcn1);
      }
      free(insized);
    }

    if(ndbbdrank[j4]>0) {
      block_sizes2[j3+begblock[rank]]=ndbbdrank[j4];
    }

  }
  MPI_Barrier(PETSC_COMM_WORLD);
  MPI_Allreduce(col_order3,col_order,VecSize, MPI_INT, MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(row_order3,row_order,VecSize, MPI_INT, MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(block_sizes2,block_sizes,ndblock, MPI_INT, MPI_SUM,PETSC_COMM_WORLD);
  free(col_order3);
  free(row_order3);
  free(block_sizes2);
  free(begblock);
  return 1;
}


