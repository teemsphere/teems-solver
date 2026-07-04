#ifndef HSL_KERNELS_H_INCLUDED
#define HSL_KERNELS_H_INCLUDED
/* Fortran kernels in hsl_kernels.f90 (wrappers around HSL MP48/MA48/MA51/MC66) */
#include <teems_solver.h>

extern void spec51m_rank_(int *INSIZE,solve_real *cntl6,int *IRN, int *JCN, solve_real *VA,int *IRNA, int *JCNA, int *KEEP,solve_real *w51, int *iw51);
extern void spec48_ssol2la_(int *INSIZE,int *IRN, int *JCN, solve_real *VA, solve_real *B, solve_real *X);
extern void spec48m_ssol2la_(int *INSIZE,int *IRN, int *JCN, solve_real *VA, solve_real *B, solve_real *X);
extern void spec48m_msol_(int *INSIZE,int *IRN, int *JCN, solve_real *VA, solve_real *B, solve_real *X, PetscInt *IRNC, PetscInt *JCNC, PetscScalar *VAC,int *IRNB,int *JCNB,PetscScalar *VALUESB,solve_real *VECBIVI,int *bivinzrow0,int *bivinzcol0);//, forint *IRNV, forint *JCNV, ha_cgetype *VAV
extern void spec48m_esol_(int *INSIZE,int *IRN, solve_real *VA,int *KEEP, solve_real *B, solve_real *X);
extern void spec48m_rpesol_(int *INSIZE,int *IRN, solve_real *VA,int *KEEP, solve_real *B, solve_real *X,solve_real *cntl,solve_real *rinfo,solve_real *error1,int *icntl,int *info,solve_real *w,int *iw);
extern void spec48_single_(fortran_int *indata,int *irn, int *jcn,solve_real *b1, solve_real *values,solve_real *x1, int *neleperrow,int *ai1, MPI_Fint *fcomm);
extern void spec48_nomc66_(fortran_int *indata, int *jcn,solve_real *b1, solve_real *values,solve_real *x1, int *neleperrow, MPI_Fint *fcomm,fortran_int *rowptrin, fortran_int *colptrin);
extern void my_spar_add3l_(solve_real *vecbivi, long int *biviindx,long int *nz1,solve_real *vecbivi0,long int *biviindx0,long int *nz0,long int *nz2);
extern void my_spar_add4l_(solve_real *vecbivi, long int *biviindx,int *irn, int *jcn,long int *nz1,solve_real *vecbivi0,long int *biviindx0,long int *nz0,long int *nz2,int *ncol);
extern void my_spar_compl_(long int *biviindx,long int *nz1,long int *biviindx0,long int *nz0,long int *nz2);
extern void prep48_alu1_(int *INSIZE,int *IRN,int* JCN,solve_real *VA,solve_real *W,int *IW,int *KEEP);
extern void prep48m_msol_(int *INSIZE,int *IRN, int *JCN, solve_real *VA, PetscInt *IRNC, PetscInt *JCNC, PetscScalar *VAC,int *IRNB,int *JCNB,PetscScalar *VALUESB,solve_real *VECBIVI,long int *bivinzrow0,int *bivinzcol0,int* jcnb1,solve_real *sol,solve_real *b,solve_real *w,int* iw,int* keep);
extern void spar_mulmin_(solve_real* sol,int* nrow,int* nz,int* irn,int* jcn,solve_real* va,solve_real* res);
extern void spar_mulnoadd_(solve_real* sol,int* nrow,int* nz,int* irn,int* jcn,solve_real* va,solve_real* res);
extern void spar_vbiviadd_(solve_real* sol,int* bvcol,long int* bvrow,long int* bvsize,int* nrow,int *ncol,int* nz,int* irn,int* jcn,solve_real* va,solve_real* res);
extern void patio_mat_(int* insizeda,int* IRN,int* JCN,solve_real* VBIVI,int* IRN1A,int* JCN1A);

#endif // HSL_KERNELS_H_INCLUDED
