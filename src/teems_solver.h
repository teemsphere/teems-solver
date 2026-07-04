#ifndef TEEMS_SOLVER_H_INCLUDED
#define TEEMS_SOLVER_H_INCLUDED

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>
#include <fcntl.h>
#include <petscksp.h>
#include <petscsys.h>
#include <petsctime.h>
#include <stdbool.h>
#include <aij.h> /* PETSc private SeqAIJ header; path supplied by makefile */
#include <petscdmda.h>
#include <omp.h>

#define NAMESIZE 256
#define DATREADLINE 150000
#define TABREADLINE 20000//2536
#define HEADERSIZE 5
#define TABLINESIZE 20000//81
#define NOPERTINSUM 5
#define MAXVARDIM 10 //maximum variable dimension
#define MAXSUPSET 12 //original
#define SORD 1
#define MAXSSIZE 187500000//1500000000/8
//To switch between single and double precision: SORD 0 single; SORD 1 double; change ha_cgetype and ha_floattype to float or double; change hsl_mp48ss.f90 accordingly. Note should change GUI too!!!
/* Definitions live in globals.c */
extern int isLinux;
extern int section_threads;
extern int max_threads;
extern double step_ratio2,step_ratio3,extrap_w1,extrap_w2,extrap_w3;
extern int steps1,steps2,steps3;
extern MPI_Comm node_comm,node_tail_comm;
extern char scratch_dir[NAMESIZE];

typedef double solve_real;
typedef int dim_t;
typedef long int offset_t;
typedef int exo_idx_t;//relax if nvarele greater than 2 billions
typedef long int fortran_int;
typedef float store_real;

//thu

//ha_cgeiof.c
typedef struct
{
  char logname[NAMESIZE];
  char filname[TABREADLINE];
} cmf_file_entry ;
int cmf_count_files(char *fname,char *comsyntax);

typedef struct
{
  dim_t dim1;
  char ch[NAMESIZE];
} datafile_labels ;
int datafile_read_labels(char *varname, char *filename,dim_t d1, datafile_labels *record);
int datafile_read_header_info(char *varname, char *filename,dim_t *vsize, char *longname,dim_t *d1);
int cmf_read(char *filename, int niodata, cmf_file_entry *iodata, char *tabfile, char *closure, char *shock);
int tab_preprocess(char *filename, char *newtabfile);
int tab_read_set_name(char *filename, char *varname, int indx, char *setname);

//ha_cgetab.c

typedef struct
{
  store_real value;
} sum_value ;
typedef struct
{
  store_real value;
  store_real initial;
  store_real substep_base;//cursol;//In stochastic case, with plus var, var0 holds x value.
} elem_value ;
char* str_rfind_any(char *line, char *finditems);
int str_rfind_ci(char *line, char *finditem);
int str_count_char(char *line, int finditem);
int str_count_ci(char *line, char *finditem);
char* str_rfind_toplevel(char *line, int finditem);
bool csv_read_ints(char *fileName,int* vec, int vecCol);

//ha_cgetab.c

typedef struct
{
  char header[HEADERSIZE];
  int fileid;
  char setname[NAMESIZE];
  char readele[TABREADLINE];
  offset_t offset;
  dim_t size;
  dim_t subsetid[MAXSUPSET];//Supersetid to be more precised
  bool intertemp;
  int intsup;
  bool regional;
  int regsup;
} set_def ;
typedef struct
{
  char setele[NAMESIZE];
  dim_t superset_pos[MAXSUPSET];
} set_element ;

typedef struct
{
  char cofname[NAMESIZE];
  offset_t offset;
  dim_t size;
  offset_t setid[MAXVARDIM];
  offset_t strides[MAXVARDIM];
  offset_t nelem;
  bool level_par;
  bool change_real;
  bool suplval;
  int gltype;//1 GE 2 GT 3 LE 4 LT
  store_real glval;
} array_def ;

typedef struct
{
  char sumname[NAMESIZE];
  offset_t summatsize;
  char sumindx[NAMESIZE];
  dim_t size;
  offset_t offset;
  offset_t sumsetid;
  char dimnames[MAXVARDIM][NAMESIZE];
  offset_t setid[MAXVARDIM];
  offset_t strides[MAXVARDIM];
} sum_def ;

typedef struct
{
  store_real value;
} elem_store ;


typedef struct
{
  exo_idx_t exo_index;
  bool is_exogenous;//1 exo
  store_real shock_value;
} closure_entry ;

typedef struct
{
  char index_name[NAMESIZE];
  offset_t setid;
  dim_t indx;
} quantifier ;


typedef struct
{
  char LinVarName[NAMESIZE];
  offset_t LinVarIndx;
  char dimnames[MAXVARDIM][NAMESIZE];
  char dimsetnames[MAXVARDIM][NAMESIZE];
  dim_t dimleadlag[MAXVARDIM];
  dim_t dimindx[MAXVARDIM];
} eq_var_ref ;

int sum_dedup_indices(char *formulain);
int tab_write_variables(char *filename, char *newtabfile,array_def *ha_var,offset_t nvar);

dim_t sets_count(char *fname);
int sets_read(char *fname, int niodata, cmf_file_entry *iodata, set_def *record,dim_t nset);
int sets_read_intertemporal(char *fname, int niodata, cmf_file_entry *iodata, set_def *record,dim_t nset);
dim_t set_union_named(set_element *ha_setele, set_def *ha_set,dim_t nset,dim_t i);
dim_t set_union_op(set_element *ha_setele, set_def *ha_set,dim_t nset,dim_t i);
dim_t set_difference(set_element *ha_setele, set_def *ha_set,dim_t nset,dim_t i);
dim_t subset_map_build(set_element *ha_setele, set_def *ha_set,dim_t nset,offset_t* contin);
char *closure_next_statement(char *commsyntax, FILE *filehandle, char *readline);
char *tab_next_statement(char *commsyntax, FILE *filehandle, char *readline,offset_t rlinesize);
char *tab_next_statement_resolved(char *commsyntax, FILE *filehandle, char *readline, elem_value *record, array_def *ha_cof,offset_t ncof,solve_real *zerodivide,offset_t rlinesize);
int str_find_ci(char *line, char *finditem);
char *str_replace_all(char *line, char *finditem, char *replitem);
int str_replace_char_all(char *line, int finditem, int replitem);
char *str_replace_all_bounded(char *line, char *finditem, char *replitem,dim_t nbuffer);
char *str_replace_first_bounded(char *line, char *finditem, char *replitem,dim_t nbuffer);
char *str_replace_first(char *line, char *finditem, char *replitem);
char *str_strip_comment(char *line, char *token);
int variables_read_defaults(char *fname, array_def *record, offset_t ncof);
offset_t variables_read(char *fname, char *commsyntax, array_def *record, offset_t ncof, set_def *ha_set,dim_t nset);
offset_t set_find_alltime(set_def *ha_set,dim_t nset);
offset_t tab_count_statements(char *fname, char *commsyntax);
offset_t closure_read(char *fname, char *commsyntax,closure_entry *ha_cgeshock, array_def *ha_var,offset_t nvar,set_def *ha_set,dim_t nset, set_element *ha_setele);
offset_t shocks_read(char *fname, char *commsyntax,closure_entry *ha_cgeshock,offset_t nvarele, array_def *ha_var,offset_t nvar,set_def *ha_set,dim_t nset, set_element *ha_setele,dim_t subints);
offset_t coefficients_read(char *fname, char *commsyntax, array_def *record, offset_t ncof, set_def *ha_set,offset_t nset);
int coef_resolve_sets(array_def *ha_cof,offset_t ncof, set_def *ha_set,dim_t nset, elem_store *ha_cofele);
offset_t data_read_files(char *fname, int niodata, cmf_file_entry *iodata, char *commsyntax,set_def *ha_set,dim_t nset, set_element *ha_setele,array_def *ha_cof,offset_t ncof, elem_store *ha_cofele,offset_t ncofele,array_def *ha_var,offset_t nvar, elem_store *ha_varele,offset_t nvarele);
int eq_replace_linvar(char *formulain,int linindx);
int eq_zero_linvar(char *formulain,int linindx);
int sum_count(char *formulain, char *commsyntax);
offset_t sum_parse(char *formulain, char *commsyntax, sum_def *sum_cof,quantifier *arSet,set_def *ha_set,dim_t nset,dim_t fdim,int j);
int sum_extract(char *formula);
int formula_normalize(char *fomulain);
offset_t subsets_read(char *fname, set_element *ha_setele, set_def *ha_set,dim_t nset);
char *str_replace_char(char *line, int finditem, int replitem);
//thu
//Moi

typedef struct
{
  dim_t Oper;//0 none; 1 *;2 /;3 +; 4 -;5 ^;71 if =;72 if >; 73 if <; 74 if <>; 75 if <=; 76 if >=;
  char TmpVarName[NAMESIZE];
  store_real TmpVarVal;
  dim_t Var1Type;//0 var cof;1 lin; 2 sum;3 lin;4 calvars; 5 number; 6 varchange; 41 id01, 42 abs, 43 loge
  offset_t Var1BegAdd;
  int Var1leadlag[MAXVARDIM];
  offset_t Var1SupSet[MAXVARDIM];
  offset_t Var1SSIndx[MAXVARDIM];
  offset_t Var1ADims[MAXVARDIM];
  store_real Var1Val;
  dim_t Var2Type;
  offset_t Var2BegAdd;
  int Var2leadlag[MAXVARDIM];
  offset_t Var2SupSet[MAXVARDIM];
  offset_t Var2SSIndx[MAXVARDIM];
  offset_t Var2ADims[MAXVARDIM];
  store_real Var2Val;

  dim_t Var3Type;
  offset_t Var3BegAdd;
  int Var3leadlag[MAXVARDIM];
  offset_t Var3SupSet[MAXVARDIM];
  offset_t Var3SSIndx[MAXVARDIM];
  offset_t Var3ADims[MAXVARDIM];
  store_real Var3Val;
} formula_op ;

solve_real formula_subst_scalar(char *var2, elem_value *record, array_def *ha_cof,offset_t ncof);
int formula_bind_operand(char *var2, set_def *ha_set,array_def *ha_cof,offset_t ncof, array_def *ha_var,offset_t nvar,offset_t ncofele,sum_def *sum_cof,int totalsum,formula_op *ha_calvar,int ha_calvarsize,quantifier *arSet,dim_t fdim,int varindex);
int leadlag_encode(char *line);
int parse_index_leadlag(char *p,int *leadlag);
int formula_compile(char *fomulain, set_def *ha_set,array_def *ha_cof, offset_t ncof, array_def *ha_var,offset_t nvar,offset_t ncofele,sum_def *sum_cof,dim_t totalsum,formula_op *ha_calvar,dim_t *ha_calvarsize,quantifier *arSet,dim_t fdim);
solve_real formula_eval(elem_value *record, set_def *ha_set,set_element *ha_setele,sum_value *ha_sumele,formula_op *ha_calvar,int ha_calvarsize,quantifier *arSet,dim_t fdim, solve_real zerodivide);
int formula_compile_pow(char *fomulain, set_def *ha_set,int npow,int ipar,array_def *ha_cof,offset_t ncof, array_def *ha_var,offset_t nvar,offset_t ncofele,sum_def *sum_cof,int totalsum,formula_op *ha_calvar,int *ha_calvarsize,quantifier *arSet,dim_t fdim);
int formula_compile_muldiv(char *fomulain, set_def *ha_set,int nmul,int ipar,array_def *ha_cof,offset_t ncof, array_def *ha_var,offset_t nvar,offset_t ncofele,sum_def *sum_cof,int totalsum,formula_op *ha_calvar,int *ha_calvarsize,quantifier *arSet,dim_t fdim);
int formula_compile_addsub(char *fomulain, set_def *ha_set,int nplu,int ipar,array_def *ha_cof,offset_t ncof, array_def *ha_var,offset_t nvar,offset_t ncofele,sum_def *sum_cof,int totalsum,formula_op *ha_calvar,int *ha_calvarsize,quantifier *arSet,dim_t fdim);
int formula_compile_if(char *fomulain, set_def *ha_set,int nif,int ipar,array_def *ha_cof,offset_t ncof, array_def *ha_var,offset_t nvar,offset_t ncofele,sum_def *sum_cof,int totalsum,formula_op *ha_calvar,int *ha_calvarsize,quantifier *arSet,dim_t fdim);
offset_t formulas_execute(char *fname, char *commsyntax,set_def *ha_set,dim_t nset, set_element *ha_setele, array_def *ha_cof,offset_t ncof,array_def *ha_var,offset_t nvar, elem_value *ha_cofvar,offset_t ncofvar,offset_t ncofele,bool IsIni);
int sum_eval(char *formulain, char *commsyntax,set_def *ha_set,dim_t nset, set_element *ha_setele,elem_value *ha_cofvar,offset_t ncofvar,offset_t ncofele, array_def *ha_cof,offset_t ncof, array_def *ha_var,offset_t nvar,sum_def *sum_cof,int totalsum,sum_value *ha_sumele,offset_t nsumele,formula_op *ha_calvar,quantifier *arSet1,dim_t fdim,int *sumindx,int j, solve_real zerodivide);
offset_t updates_apply(char *fname,set_def *ha_set,dim_t nset, set_element *ha_setele, array_def *ha_cof,offset_t ncof,array_def *ha_var,offset_t nvar, elem_value *ha_cofvar,offset_t ncofvar,offset_t ncofele,int midpoint);
offset_t updates_apply_product(char *fname,set_def *ha_set,dim_t nset, set_element *ha_setele, array_def *ha_cof,offset_t ncof,array_def *ha_var,offset_t nvar, elem_value *ha_cofvar,offset_t ncofvar,offset_t ncofele);
//lin

int eq_sum_parse(char *formulain, char *commsyntax, sum_def *sum_cof,quantifier *arSet,set_def *ha_set,dim_t nset,dim_t fdim,int j);
int eq_sum_eval(char *formulain, char *commsyntax,set_def *ha_set,dim_t nset, set_element *ha_setele,elem_value *ha_cofvar,offset_t ncofvar,offset_t ncofele, array_def *ha_cof,offset_t ncof,array_def *ha_var,offset_t nvar,sum_def *sum_cof,int totalsum,sum_value *ha_sumele,offset_t nsumele,formula_op *ha_calvar,quantifier *arSet1,dim_t fdim,int *sumindx,int j, solve_real zerodivide);
int eq_sum_replace(char *formulain, char *commsyntax,int LinIndx, eq_var_ref *LinVars,array_def *ha_var);
int jacobian_fill(char *fname, char *commsyntax,set_def *ha_set,offset_t nset, set_element *ha_setele, array_def *ha_cof,offset_t ncof,array_def *ha_var,offset_t nvar, elem_value *ha_cofvar,offset_t ncofvar,offset_t ncofele,closure_entry *ha_cgeshock,offset_t ndblock,offset_t alltimeset,offset_t allregset,offset_t *ha_eqadd,offset_t *counteq,offset_t nintraeq,Mat A,Mat B);
int eq_linvar_read(char *formulain,eq_var_ref *LinVars,int linindx,array_def *ha_var);
int jacobian_preallocate(char *fname, char *commsyntax,set_def *ha_set,dim_t nset,set_element *ha_setele,array_def *ha_cof,offset_t ncof,array_def *ha_var,offset_t nvar,elem_value *ha_cofvar,offset_t ncofvar,offset_t ncofele, offset_t nexo,closure_entry *ha_cgeshock,offset_t ndblock,offset_t alltimeset,offset_t allregset,bool *ha_eqint,offset_t *ha_eqadd,dim_t *ha_eqtime,dim_t *ha_eqreg,offset_t *counteq,offset_t nintraeq,bool *sbbd_overrid,PetscInt Istart,PetscInt Iend,PetscInt *dnz,PetscInt *dnnz,PetscInt *onz,PetscInt *onnz,PetscInt *dnzB,PetscInt *dnnzB,PetscInt *onzB,PetscInt *onnzB,int nesteddbbd);
int equation_order_read(char *fname, char *commsyntax,set_def *ha_set,dim_t nset,set_element *ha_setele,array_def *ha_cof,offset_t ncof,array_def *ha_var,offset_t nvar,elem_value *ha_cofvar,offset_t ncofvar,offset_t ncofele,closure_entry *ha_cgeshock,bool *var_inter,array_def *ha_eq,bool *ha_eqint,dim_t *eq_orderintra,dim_t *eq_orderreg,offset_t allregset,offset_t alltimeset,dim_t *orderintra,dim_t *orderreg);
int equation_order_read_nested(char *fname, char *commsyntax,set_def *ha_set,dim_t nset,set_element *ha_setele,array_def *ha_cof,offset_t ncof,array_def *ha_var,offset_t nvar,elem_value *ha_cofvar,offset_t ncofvar,offset_t ncofele,closure_entry *ha_cgeshock,bool *var_inter,array_def *ha_eq,bool *ha_eqint,dim_t *eq_orderintra,dim_t *eq_orderreg,offset_t allregset,offset_t alltimeset,dim_t *orderintra,dim_t *orderreg);
bool ndbbd_block_solve(PetscInt rank, int begmat,int nreg,int * insize,int insizes, Mat **submatCij,Mat **submatBij,solve_real *b,solve_real *sol,bool ifremove,char** fn01,char** fn02, char** fn03);
bool ndbbd_block_solve_mem(PetscInt rank, int begmat,int nreg,int * insize,int insizes, Mat **submatCij,Mat **submatBij,solve_real *b,solve_real *sol,int** irnereg,int** keepreg,solve_real** valereg,solve_real *cntl,solve_real *rinfo,solve_real *error1,int *icntl,int *info,solve_real *w,int *iw,solve_real *b02);

int ndbbd_order_presolve(Mat A, offset_t VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int nreg, int ntime, offset_t nvarele, offset_t *ha_eqadd,int *ha_rows,int *ha_cols, offset_t ndblock,int *ha_ndblocks, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,dim_t laDi,solve_real cntl6,PetscInt* ndbbdrank,PetscBool presol);
int dbbd_order(Mat A, offset_t VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend, offset_t nvarele, offset_t *ha_eqadd,int *ha_rows,int *ha_cols, offset_t ndblock,int *ha_ndblocks, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,solve_real cntl6);
int dbbd_solve(Mat A, Vec b, solve_real *x1, offset_t VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int *ha_rows,int *ha_cols, offset_t ndblock,int *ha_ndblocks, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,dim_t laD,PetscReal cntl3);//,bool iter
int ndbbd_order(Mat A, offset_t VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int nreg, int ntime, offset_t nvarele, offset_t *ha_eqadd,int *ha_rows,int *ha_cols, offset_t ndblock,int *ha_ndblocks, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,dim_t laDi,solve_real cntl6,PetscInt* ndbbdrank,PetscBool presol);
int ndbbd_presolve(Mat A, Vec b, solve_real *x1, offset_t VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int *ha_rows,int *ha_cols, offset_t ndblock,offset_t nreg,offset_t ntime,int *ha_ndblocks, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,dim_t laDi,dim_t laD,PetscReal cntl3,PetscReal cntl6,PetscBool presol);//,bool iter
int ndbbd_solve(Mat A, Vec b, solve_real *x1, offset_t VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int *ha_rows,int *ha_cols, offset_t ndblock,offset_t nreg,offset_t ntime,int *ha_ndblocks, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,dim_t laDi,dim_t laD,PetscReal cntl3,PetscReal cntl6,PetscBool presol);//,bool iter
int reduce_to_rank(solve_real *vecbivi,fortran_int vecbivisize,PetscInt mpisize,PetscInt rank,PetscInt targetrank);
int reduce_to_rank_nocompress(solve_real *vecbivi,fortran_int vecbivisize,PetscInt mpisize,PetscInt rank,PetscInt targetrank);
int outputs_write_csv(char *filename, char *newdatlogname, char *newdatfile,set_def *ha_set,dim_t nset, set_element *ha_setele,array_def *ha_cof,offset_t ncof,offset_t ncofele,array_def *ha_var,offset_t nvar,offset_t nvarele, elem_value *ha_cofvar);

offset_t subinterval_update(PetscInt rank,char *fname,set_def *ha_set,dim_t nset, set_element *ha_setele, array_def *ha_cof,offset_t ncof,array_def *ha_var,offset_t nvar, elem_value *ha_cofvar,offset_t ncofvar,offset_t ncofele,closure_entry *ha_cgeshock,offset_t nvarele,int laA,dim_t subints,bool IsIni,int IsSplint,int nsteps);
bool cubic_spline(solve_real* y,solve_real* x,solve_real sx0,solve_real sxn,int size,solve_real* w,int laA);

bool solve_johansen(PetscBool nohsl,PetscInt VecSize,Mat A,PetscInt dnz,PetscInt* dnnz,PetscInt onz,PetscInt* onnz,Mat B,PetscInt dnzB,PetscInt* dnnzB,PetscInt onzB,PetscInt* onnzB,Vec vecb,Vec vece,PetscInt rank,PetscInt rank_hsl,PetscInt mpisize,char* tabfile, char *commsyntax,set_def *ha_set,dim_t nset, set_element *ha_setele, array_def *ha_cof,offset_t ncof,array_def *ha_var,offset_t nvar, elem_value **ha_cofvar2,offset_t ncofvar,offset_t ncofele,offset_t nvarele,closure_entry **ha_cgeshock2,offset_t alltimeset,offset_t allregset,offset_t nintraeq,dim_t matsol,PetscInt Istart,PetscInt Iend,  offset_t nreg, offset_t ntime, offset_t *ha_eqadd, offset_t ndblock, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,dim_t laDi,dim_t laD,PetscReal cntl3,PetscReal cntl6,PetscBool presol,dim_t nesteddbbd,int localsize,PetscInt *ndbbddrank1,fortran_int* indata,dim_t mc66,fortran_int *ptx,struct timeval begintime,solve_real **xcf2);
bool solve_modified_midpoint(PetscBool nohsl,PetscInt VecSize,Mat* A,PetscInt dnz,PetscInt* dnnz,PetscInt onz,PetscInt* onnz,Mat* B,PetscInt dnzB,PetscInt* dnnzB,PetscInt onzB,PetscInt* onnzB,Vec* vecb,Vec *vece,PetscInt rank,PetscInt rank_hsl,PetscInt mpisize,char* tabfile, char *commsyntax,set_def *ha_set,dim_t nset, set_element *ha_setele, array_def *ha_cof,offset_t ncof,array_def *ha_var,offset_t nvar, elem_value **ha_cofvar2,offset_t ncofvar,offset_t ncofele,offset_t nvarele,closure_entry **ha_cgeshock2,offset_t alltimeset,offset_t allregset,offset_t nintraeq,dim_t matsol,PetscInt Istart,PetscInt Iend,  offset_t nreg, offset_t ntime, offset_t *ha_eqadd, offset_t ndblock, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,dim_t laDi,dim_t laD,PetscReal cntl3,PetscReal cntl6,PetscBool presol,dim_t nesteddbbd,int localsize,PetscInt *ndbbddrank1,fortran_int* indata,dim_t mc66,fortran_int *ptx,struct timeval begintime,dim_t subints,MPI_Fint fcomm,solve_real **xcf2,int Isbiupd);
#endif // TEEMS_SOLVER_H_INCLUDED

