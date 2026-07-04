#include <teems_solver.h>
#include <hsl_kernels.h>

int dbbd_order(Mat A, offset_t VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend, offset_t nvarele, offset_t *ha_eqadd,int *ha_rows,int *ha_cols, offset_t ndblock,int *ha_ndblocks, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,solve_real cntl6) {
  IS *rowindices,*colindices;//,isrow,iscol;
  PetscInt bfirst,bend,sblockin,nmatin,nmatinplus,nrowcolin,sumrowcolin;
  Mat *submatA;
  PetscInt i,j,j0,j1,j2,j3,j4,nrow,ncol,nz,nz1,*ai,*aj,la;
  PetscScalar *vals;
  PetscErrorCode ierr;
  PetscViewer viewer;
  offset_t lasize;
  int *ha_rows1= (int *) calloc (VecSize,sizeof(int));
  int *ha_cols1= (int *) calloc (VecSize,sizeof(int));
  int *ha_ndblocks1= (int *) calloc (ndblock,sizeof(int));

  //*****New order test**********//
  nmatin=(offset_t)ndblock/mpisize;
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
    solve_real *values= (solve_real *) calloc (lasize,sizeof(solve_real));
    int *keep=(int *) malloc ((nrow+5*ncol+4*ncol+7)*sizeof(int));//ICNTL(6)=1
  solve_real *w51;
  if(nrow>ncol) w51=(solve_real *) malloc ((5*nrow)*sizeof(solve_real));
  else w51=(solve_real *) malloc ((5*ncol)*sizeof(solve_real));
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


int ndbbd_order_presolve(Mat A, offset_t VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int nreg, int ntime, offset_t nvarele, offset_t *ha_eqadd,int *ha_rows,int *ha_cols, offset_t ndblock,int *ha_ndblocks, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,dim_t laDi,solve_real cntl6,PetscInt* ndbbdrank,PetscBool presol) {
  FILE *presolfile;
  char j1name[1024],filename[1024],rankname[1024];
  size_t frd;
  IS *rowindices=NULL,*colindices=NULL,*rowindicesD=NULL,*colindicesD=NULL;//,isrow,iscol;
  PetscInt bfirst,bend,sblockin,nmatin,nmatint,nmatinplus,nrowcolin,sumrowcolin,rfirst,cfirst;
  Mat *submatA=NULL;//,*submatD=NULL;
  PetscInt i,j,j0,j1,j2,j3,j4,j5,j6,j7,j8,nrow,ncol,nz,nz1,*ai,*aj,la;
  offset_t lasize;
  solve_real cntl6in;
  PetscScalar *vals;
  PetscErrorCode ierr;
  PetscViewer viewer;
  MatInfo           matinfo;
  MatGetInfo(A,MAT_LOCAL,&matinfo);
  printf("rank %d matinfo.nz_used %g\n",rank,matinfo.nz_used);
  int *ha_ndblocks1= (int *) calloc (ndblock,sizeof(int));
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
    solve_real *values= (solve_real *) calloc (lasizemax,sizeof(solve_real));
  int *keep=(int *) malloc ((nrowmax+5*ncolmax+4*ncolmax+7)*sizeof(int));//ICNTL(6)=1
  solve_real *w51;
  if(nrowmax>ncolmax) w51=(solve_real *) malloc ((5*nrowmax)*sizeof(solve_real));
  else w51=(solve_real *) malloc ((5*ncolmax)*sizeof(solve_real));
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
  memcpy(counteq,counteq2,(ndblock+1)*sizeof(offset_t));
  memcpy(counteqnoadd,counteqnoadd2,(ndblock)*sizeof(offset_t));
  memcpy(countvarintra1,countvarintra4,(ndblock+1)*sizeof(offset_t));
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

int ndbbd_order(Mat A, offset_t VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int nreg, int ntime, offset_t nvarele, offset_t *ha_eqadd,int *ha_rows,int *ha_cols, offset_t ndblock,int *ha_ndblocks, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,dim_t laDi,solve_real cntl6,PetscInt* ndbbdrank,PetscBool presol) {
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
  nmatint=(offset_t)ntime/mpisize;
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


