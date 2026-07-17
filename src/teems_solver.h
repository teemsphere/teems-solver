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
#define SORD 1 /* 1 = double solve precision; to switch to single change
                  solve_real/store_real and FSORD in hsl_kernels.f90 */
#define MAXSSIZE 187500000//1500000000/8
/* Definitions live in globals.c */
extern int verbosity; /* -verbosity: 0 = errors/warnings + results only,
                         1 = phase progress (default), 2 = per-rank/per-block
                         debug detail. Errors, warnings, and the accuracy
                         summary teems-R parses are never gated. */
/* logmsg(level, ...): print when verbosity >= level. Not for errors or
   warnings (plain printf, "Error:"/"Warning:" prefix) and never inside
   per-element loops — hot-loop prints are removed, not gated. */
#define logmsg(lvl, ...) do{ if(verbosity>=(lvl)) printf(__VA_ARGS__); }while(0)
extern int inmemory; /* -inmemory: keep value arrays resident instead of spilling to scratch */
extern int section_threads;
extern int max_threads;
extern double step_ratio2,step_ratio3,extrap_w1,extrap_w2,extrap_w3;
extern int steps1,steps2,steps3;
extern MPI_Comm node_comm,node_tail_comm;
extern char scratch_dir[NAMESIZE];

typedef double solve_real;   /* linear-solve precision */
typedef int dim_t;           /* set sizes, dimension counts */
typedef long int offset_t;   /* element offsets into value arrays */
typedef int exo_idx_t;       /* exogenous-variable index; widen if nvarele > 2^31 */
typedef long int fortran_int;/* INTEGER(8) interop with hsl_kernels.f90 */
typedef float store_real;    /* coefficient storage precision (halves memory traffic) */

/* -matsol matrix method (Ha & Kompas 2016; Kompas & Ha 2019).
   Values match teems-R's matrix_method argument. */
enum matrix_method { MM_LU=0, MM_SBBD=1, MM_DBBD=2, MM_NDBBD=3 };
/* -solmed solution method (GEMPACK manual; Pearson 1991) */
enum solution_method { SM_GRAGG=1, SM_EULER=2, SM_JOHANSEN=10, SM_NOSOLVE=100 };
/* array_def.gltype: bound imposed on levels values */
enum bound_type { BT_NONE=0, BT_GE=1, BT_GT=2, BT_LE=3, BT_LT=4 };
/* formula_op.Oper: compiled formula operation */
enum op_code { OP_LOAD=0, OP_MUL=1, OP_DIV=2, OP_ADD=3, OP_SUB=4, OP_POW=5,
               OP_IF_EQ=71, OP_IF_GT=72, OP_IF_LT=73, OP_IF_NE=74,
               OP_IF_LE=75, OP_IF_GE=76 };
/* formula_op operand types (Var1Type/Var2Type/Var3Type) */
enum operand_type { OT_ARRAY=0, OT_LINVAR=1, OT_SUM=2, OT_LINVAR2=3,
                    OT_TEMP=4, OT_CONST=5, OT_CHANGE=6,
                    OT_TEMP_ID01=41, OT_TEMP_ABS=42, OT_TEMP_LOG=43 };


/* ================= cmf_io.c — command (CMF) file and data I/O ========== */

/* one "file <logname> <path>;" statement from the CMF */
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

/* ================= value records ======================================= */

/* one element of an evaluated SUM(...) */
typedef struct
{
  store_real value;
} sum_value ;
/* one element of a coefficient or variable array */
typedef struct
{
  store_real value;
  store_real initial;      /* pre-simulation (levels) value */
  store_real substep_base; /* base for midpoint update */
} elem_value ;
char* str_rfind_any(char *line, char *finditems);
int str_rfind_ci(char *line, char *finditem);
int str_count_char(char *line, int finditem);
int str_count_ci(char *line, char *finditem);
char* str_rfind_toplevel(char *line, int finditem);

/* ================= tab_parse.c — TAB-language model description ======== */

/* a SET statement */
typedef struct
{
  char header[HEADERSIZE];
  int fileid;
  char setname[NAMESIZE];
  char readele[TABREADLINE];
  offset_t offset;
  dim_t size;
  dim_t subsetid[MAXSUPSET]; /* ids of supersets this set maps into */
  bool intertemp;
  int intsup;
  bool regional;
  int regsup;
} set_def ;
/* one set element with its position in each superset */
typedef struct
{
  char setele[NAMESIZE];
  dim_t superset_pos[MAXSUPSET];
} set_element ;

/* a COEFFICIENT or VARIABLE declaration: name, dimensionality and layout */
typedef struct
{
  char cofname[NAMESIZE];
  offset_t offset;             /* start of this array in the value vector */
  dim_t size;                  /* number of dimensions */
  offset_t setid[MAXVARDIM];   /* set of each dimension */
  offset_t strides[MAXVARDIM]; /* row-major strides */
  offset_t nelem;
  bool level_par;              /* levels (not percentage-change) quantity */
  bool change_real;
  bool suplval;                /* values supplied by READ/FORMULA */
  int gltype;                  /* enum bound_type */
  store_real glval;            /* bound value */
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


/* closure and shock for one variable element.  exo_index is the matrix
   column for exogenous/endogenous elements; for backsolved elements it
   is the compact index into the per-step recovered-value array. */
typedef struct
{
  exo_idx_t exo_index;
  bool is_exogenous;
  bool is_backsolved;
  store_real shock_value;
} closure_entry ;

/* one "backsolve <var> using <eq> ;" statement (GEMPACK manual 10.16,
   14.1.3): the variable and its defining equation are eliminated from
   the condensed system; the solver recovers the variable's per-step
   values from the equation after each solve, before the data updates. */
typedef struct
{
  char eqname[NAMESIZE];   /* nominated defining equation */
  offset_t varindx;        /* index into vars[] */
  offset_t elem_base;      /* first slot in the recovered-value array */
} backsolve_def ;
extern backsolve_def *backsolves; /* nominated backsolve pairs, TAB order */
extern int nbacksolve;            /* number of backsolve statements */
extern offset_t nbselems;         /* total backsolved variable elements */
/* "equation" scan filter: SKIP excludes the nominated defining equations
   (every consumer of the scan then sees only the condensed system);
   ONLY inverts the filter for the recovery-program build. */
enum bs_scan_mode { BS_SCAN_SKIP=0, BS_SCAN_ONLY=1 };
extern int backsolve_scan_mode;
offset_t backsolve_read(char *fname, array_def *vars, offset_t nvar, closure_entry *closure_vals);
int backsolve_validate_refs(char *fname, array_def *vars);
int tab_equation_name(char *stmt, char *eqname);

/* one (all,index,SET) quantifier with its current position */
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
int tab_write_variables(char *filename, char *newtabfile,array_def *vars,offset_t nvar);

dim_t sets_count(char *fname);
int sets_read(char *fname, int niodata, cmf_file_entry *iodata, set_def *record,dim_t nset);
int sets_read_intertemporal(char *fname, int niodata, cmf_file_entry *iodata, set_def *record,dim_t nset);
dim_t set_union_named(set_element *set_elems, set_def *sets,dim_t nset,dim_t i);
dim_t set_expr_build(set_element *set_elems, set_def *sets,dim_t nset,dim_t i); /* "@<expr>" GEMPACK set expressions */
void set_equality_build(set_element *set_elems, set_def *sets,dim_t i); /* "=<idx>" SET <new> = <old>; */
dim_t set_union_op(set_element *set_elems, set_def *sets,dim_t nset,dim_t i);
dim_t set_difference(set_element *set_elems, set_def *sets,dim_t nset,dim_t i);
dim_t subset_map_build(set_element *set_elems, set_def *sets,dim_t nset,offset_t* contin);
char *closure_next_statement(char *commsyntax, FILE *filehandle, char *readline);
char *tab_next_statement(char *commsyntax, FILE *filehandle, char *readline,offset_t rlinesize);
char *tab_next_statement_resolved(char *commsyntax, FILE *filehandle, char *readline, elem_value *record, array_def *coefs,offset_t ncof,solve_real *zerodivide,offset_t rlinesize);
int str_find_ci(char *line, char *finditem);
char *str_replace_all(char *line, char *finditem, char *replitem);
int str_replace_char_all(char *line, int finditem, int replitem);
char *str_replace_all_bounded(char *line, char *finditem, char *replitem,dim_t nbuffer);
char *str_replace_first_bounded(char *line, char *finditem, char *replitem,dim_t nbuffer);
char *str_replace_first(char *line, char *finditem, char *replitem);
char *str_strip_comment(char *line, char *token);
int variables_read_defaults(char *fname, array_def *record, offset_t ncof);
offset_t variables_read(char *fname, char *commsyntax, array_def *record, offset_t ncof, set_def *sets,dim_t nset);
offset_t set_find_alltime(set_def *sets,dim_t nset);
offset_t tab_count_statements(char *fname, char *commsyntax);
offset_t closure_read(char *fname, char *commsyntax,closure_entry *closure_vals, array_def *vars,offset_t nvar,set_def *sets,dim_t nset, set_element *set_elems);
offset_t shocks_read(char *fname, char *commsyntax,closure_entry *closure_vals,offset_t nvarele, array_def *vars,offset_t nvar,set_def *sets,dim_t nset, set_element *set_elems,dim_t subints);
offset_t coefficients_read(char *fname, char *commsyntax, array_def *record, offset_t ncof, set_def *sets,offset_t nset);
int coef_resolve_sets(array_def *coefs,offset_t ncof, set_def *sets,dim_t nset, elem_store *coef_store);
offset_t data_read_files(char *fname, int niodata, cmf_file_entry *iodata, char *commsyntax,set_def *sets,dim_t nset, set_element *set_elems,array_def *coefs,offset_t ncof, elem_store *coef_store,offset_t ncofele,array_def *vars,offset_t nvar, elem_store *var_store,offset_t nvarele);
int eq_replace_linvar(char *formulain,int linindx);
int eq_zero_linvar(char *formulain,int linindx);
int sum_count(char *formulain, char *commsyntax);
offset_t sum_parse(char *formulain, char *commsyntax, sum_def *sum_cof,quantifier *arSet,set_def *sets,dim_t nset,dim_t fdim,int j);
int sum_extract(char *formula);
int formula_normalize(char *fomulain);
offset_t subsets_read(char *fname, set_element *set_elems, set_def *sets,dim_t nset);
char *str_replace_char(char *line, int finditem, int replitem);

/* ================= formula.c — FORMULA compile/evaluate, UPDATE ======== */

/* per-dimension addressing record of one operand: the eval loop reads
   all four fields together per dim, so they are interleaved (one record
   per dim) rather than held in four parallel arrays */
typedef struct
{
  offset_t ADims;            /* stride of this dim in the value array */
  int SupSet;                /* 1 = index via superset_pos, 0 = direct */
  int SSIndx;                /* column into superset_pos (< MAXSUPSET) */
  int leadlag;               /* intertemporal lead/lag shift */
} dim_addr ;

/* one operation of a compiled formula program (interpreted per element) */
typedef struct
{
  dim_t Oper;                /* enum op_code */
  store_real TmpVarVal;
  dim_t Var1Type;            /* enum operand_type */
  offset_t Var1BegAdd;
  dim_addr Var1Dims[MAXVARDIM];
  store_real Var1Val;
  dim_t Var2Type;
  offset_t Var2BegAdd;
  dim_addr Var2Dims[MAXVARDIM];
  store_real Var2Val;

  dim_t Var3Type;
  offset_t Var3BegAdd;
  dim_addr Var3Dims[MAXVARDIM];
  store_real Var3Val;
  /* compile-time only (generated temp name); kept last and short so the
     eval-hot fields above stay cache-dense */
  char TmpVarName[64];
} formula_op ;

solve_real formula_subst_scalar(char *var2, elem_value *record, array_def *coefs,offset_t ncof);
int formula_bind_operand(char *var2, set_def *sets,array_def *coefs,offset_t ncof, array_def *vars,offset_t nvar,offset_t ncofele,sum_def *sum_cof,int totalsum,formula_op *ops,int nops,quantifier *arSet,dim_t fdim,int varindex);
int leadlag_encode(char *line);
int parse_index_leadlag(char *p,int *leadlag);
int formula_compile(char *fomulain, set_def *sets,array_def *coefs, offset_t ncof, array_def *vars,offset_t nvar,offset_t ncofele,sum_def *sum_cof,dim_t totalsum,formula_op *ops,dim_t *nops,quantifier *arSet,dim_t fdim);
solve_real formula_eval(elem_value *record, set_def *sets,set_element *set_elems,sum_value *sum_vals,formula_op *ops,int nops,quantifier *arSet,dim_t fdim, solve_real zerodivide);
int formula_compile_pow(char *fomulain, set_def *sets,int npow,int ipar,array_def *coefs,offset_t ncof, array_def *vars,offset_t nvar,offset_t ncofele,sum_def *sum_cof,int totalsum,formula_op *ops,int *nops,quantifier *arSet,dim_t fdim);
int formula_compile_muldiv(char *fomulain, set_def *sets,int nmul,int ipar,array_def *coefs,offset_t ncof, array_def *vars,offset_t nvar,offset_t ncofele,sum_def *sum_cof,int totalsum,formula_op *ops,int *nops,quantifier *arSet,dim_t fdim);
int formula_compile_addsub(char *fomulain, set_def *sets,int nplu,int ipar,array_def *coefs,offset_t ncof, array_def *vars,offset_t nvar,offset_t ncofele,sum_def *sum_cof,int totalsum,formula_op *ops,int *nops,quantifier *arSet,dim_t fdim);
int formula_compile_if(char *fomulain, set_def *sets,int nif,int ipar,array_def *coefs,offset_t ncof, array_def *vars,offset_t nvar,offset_t ncofele,sum_def *sum_cof,int totalsum,formula_op *ops,int *nops,quantifier *arSet,dim_t fdim);
offset_t formulas_execute(char *fname, char *commsyntax,set_def *sets,dim_t nset, set_element *set_elems, array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar, elem_value *elem_vals,offset_t ncofvar,offset_t ncofele,bool IsIni);
int sum_eval(char *formulain, char *commsyntax,set_def *sets,dim_t nset, set_element *set_elems,elem_value *elem_vals,offset_t ncofvar,offset_t ncofele, array_def *coefs,offset_t ncof, array_def *vars,offset_t nvar,sum_def *sum_cof,int totalsum,sum_value *sum_vals,offset_t nsumele,formula_op *ops,quantifier *arSet1,dim_t fdim,int *sumindx,int j, solve_real zerodivide);
offset_t updates_apply(char *fname,set_def *sets,dim_t nset, set_element *set_elems, array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar, elem_value *elem_vals,offset_t ncofvar,offset_t ncofele,int midpoint);
offset_t updates_apply_product(char *fname,set_def *sets,dim_t nset, set_element *set_elems, array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar, elem_value *elem_vals,offset_t ncofvar,offset_t ncofele);

/* ============ jacobian.c — first-order derivative matrix assembly ======
   (Ha & Kompas 2016 §5; Kompas & Ha 2019) */

int eq_sum_parse(char *formulain, char *commsyntax, sum_def *sum_cof,quantifier *arSet,set_def *sets,dim_t nset,dim_t fdim,int j);
int eq_sum_replace(char *formulain, char *commsyntax,int LinIndx, eq_var_ref *LinVars,array_def *vars);
int jacobian_fill(char *fname, char *commsyntax,set_def *sets,offset_t nset, set_element *set_elems, array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar, elem_value *elem_vals,offset_t ncofvar,offset_t ncofele,closure_entry *closure_vals,offset_t ndblock,offset_t alltimeset,offset_t allregset,offset_t *eq_addr,offset_t *counteq,offset_t nintraeq,Mat A,Mat B);
void jacobian_cache_free(void); /* release the per-rank compiled-statement cache */
/* Recover the backsolved variables' per-step values from their retained
   defining equations (GEMPACK 14.1.3: after the condensed solve, before
   the data updates).  x = this step's solution vector; exo_z = per-element
   exogenous per-step changes captured at the vece fill sites; bsvals
   (nbselems) receives the recovered changes at each pair's elem_base.
   Returns -1 on a zero pivot or a malformed defining equation. */
int backsolve_recover(char *fname, char *commsyntax,set_def *sets,offset_t nset, set_element *set_elems, array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar, elem_value *elem_vals,offset_t ncofele,closure_entry *closure_vals,solve_real *x,solve_real *exo_z,solve_real *bsvals);
void backsolve_cache_free(void);
int eq_linvar_read(char *formulain,eq_var_ref *LinVars,int linindx,array_def *vars);
int jacobian_preallocate(char *fname, char *commsyntax,set_def *sets,dim_t nset,set_element *set_elems,array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar,elem_value *elem_vals,offset_t ncofvar,offset_t ncofele, offset_t nexo,closure_entry *closure_vals,offset_t ndblock,offset_t alltimeset,offset_t allregset,bool *eq_intertemp,offset_t *eq_addr,dim_t *eq_time,dim_t *eq_reg,offset_t *counteq,offset_t nintraeq,bool *sbbd_overrid,PetscInt Istart,PetscInt Iend,PetscInt *dnz,PetscInt *dnnz,PetscInt *onz,PetscInt *onnz,PetscInt *dnzB,PetscInt *dnnzB,PetscInt *onzB,PetscInt *onnzB,int nesteddbbd);
int equation_order_read(char *fname, char *commsyntax,set_def *sets,dim_t nset,set_element *set_elems,array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar,elem_value *elem_vals,offset_t ncofvar,offset_t ncofele,closure_entry *closure_vals,bool *var_inter,bool *ele_inter,array_def *eq_defs,bool *eq_intertemp,dim_t *eq_orderintra,dim_t *eq_orderreg,offset_t allregset,offset_t alltimeset,dim_t *orderintra,dim_t *orderreg);
int equation_order_read_nested(char *fname, char *commsyntax,set_def *sets,dim_t nset,set_element *set_elems,array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar,elem_value *elem_vals,offset_t ncofvar,offset_t ncofele,closure_entry *closure_vals,bool *var_inter,bool *ele_inter,array_def *eq_defs,bool *eq_intertemp,dim_t *eq_orderintra,dim_t *eq_orderreg,offset_t allregset,offset_t alltimeset,dim_t *orderintra,dim_t *orderreg);
/* ======= block_order.c / block_solve.c — (N)DBBD ordering and solve ====
   Doubly Bordered Block Diagonal decomposition per Ha & Kompas 2016 and
   Kompas & Ha 2019: reorder the Jacobian into diagonal blocks plus
   borders, LU-factor blocks in parallel (HSL MA48), form and solve the
   interface problem, then back-solve.  laA/laDi/laD control workspace
   sizing; cntl3/cntl6 are HSL control parameters. */

bool ndbbd_block_solve(PetscInt rank, int begmat,int nreg,int * insize,int insizes, Mat **submatCij,Mat **submatBij,solve_real *b,solve_real *sol,bool ifremove,char** fn01,char** fn02, char** fn03);
bool ndbbd_block_solve_mem(PetscInt rank, int begmat,int nreg,int * insize,int insizes, Mat **submatCij,Mat **submatBij,solve_real *b,solve_real *sol,int** irnereg,int** keepreg,solve_real** valereg,solve_real *cntl,solve_real *rinfo,solve_real *error1,int *icntl,int *info,solve_real *w,int *iw,solve_real *b02);

int ndbbd_order_presolve(Mat A, offset_t VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int nreg, int ntime, offset_t nvarele, offset_t *eq_addr,int *row_order,int *col_order, offset_t ndblock,int *block_sizes, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,dim_t laDi,solve_real cntl6,PetscInt* ndbbdrank,PetscBool presol);
int dbbd_order(Mat A, offset_t VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend, offset_t nvarele, offset_t *eq_addr,int *row_order,int *col_order, offset_t ndblock,int *block_sizes, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,solve_real cntl6);
int dbbd_solve(Mat A, Vec b, solve_real *x1, offset_t VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int *row_order,int *col_order, offset_t ndblock,int *block_sizes, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,dim_t laD,PetscReal cntl3);//,bool iter
int ndbbd_order(Mat A, offset_t VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int nreg, int ntime, offset_t nvarele, offset_t *eq_addr,int *row_order,int *col_order, offset_t ndblock,int *block_sizes, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,dim_t laDi,solve_real cntl6,PetscInt* ndbbdrank,PetscBool presol);
int ndbbd_presolve(Mat A, Vec b, solve_real *x1, offset_t VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int *row_order,int *col_order, offset_t ndblock,offset_t nreg,offset_t ntime,int *block_sizes, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,dim_t laDi,dim_t laD,PetscReal cntl3,PetscReal cntl6,PetscBool presol);//,bool iter
int ndbbd_solve(Mat A, Vec b, solve_real *x1, offset_t VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int *row_order,int *col_order, offset_t ndblock,offset_t nreg,offset_t ntime,int *block_sizes, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,dim_t laDi,dim_t laD,PetscReal cntl3,PetscReal cntl6,PetscBool presol);//,bool iter
int reduce_to_rank(solve_real *vecbivi,fortran_int vecbivisize,PetscInt mpisize,PetscInt rank,PetscInt targetrank);
int reduce_to_rank_nocompress(solve_real *vecbivi,fortran_int vecbivisize,PetscInt mpisize,PetscInt rank,PetscInt targetrank);
int outputs_write_csv(char *filename, char *newdatlogname, char *newdatfile,set_def *sets,dim_t nset, set_element *set_elems,array_def *coefs,offset_t ncof,offset_t ncofele,array_def *vars,offset_t nvar,offset_t nvarele, elem_value *elem_vals);

/* ============ solve_drivers.c — solution methods ======================= */

/* one-step (Johansen 1960) solution of the linearized system */
bool solve_johansen(PetscBool nohsl,PetscInt VecSize,Mat A,PetscInt dnz,PetscInt* dnnz,PetscInt onz,PetscInt* onnz,Mat B,PetscInt dnzB,PetscInt* dnnzB,PetscInt onzB,PetscInt* onnzB,Vec vecb,Vec vece,PetscInt rank,PetscInt rank_hsl,PetscInt mpisize,char* tabfile, char *commsyntax,set_def *sets,dim_t nset, set_element *set_elems, array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar, elem_value **elem_vals2,offset_t ncofvar,offset_t ncofele,offset_t nvarele,closure_entry **closure_vals2,offset_t alltimeset,offset_t allregset,offset_t nintraeq,dim_t matsol,PetscInt Istart,PetscInt Iend,  offset_t nreg, offset_t ntime, offset_t *eq_addr, offset_t ndblock, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,dim_t laDi,dim_t laD,PetscReal cntl3,PetscReal cntl6,dim_t nesteddbbd,int localsize,PetscInt *ndbbddrank1,fortran_int* indata,dim_t mc66,fortran_int *ptx,struct timeval begintime,solve_real **xcf2);
/* Multistep driver with Richardson extrapolation over
   steps1/steps2/steps3. solmethod selects the stepping scheme:
   SM_GRAGG (smoothed modified midpoint, Pearson 1991 eq. 6.1 /
   Alg. 7.1.2): Euler start, midpoint leapfrog, terminal smoothing
   pass, h^2 error series; SM_EULER: forward Euler on every substep,
   no smoothing pass, h error series (extrapolation weights differ
   accordingly). */
bool solve_gragg(PetscBool nohsl,PetscInt VecSize,Mat* A,PetscInt dnz,PetscInt* dnnz,PetscInt onz,PetscInt* onnz,Mat* B,PetscInt dnzB,PetscInt* dnnzB,PetscInt onzB,PetscInt* onnzB,Vec* vecb,Vec *vece,PetscInt rank,PetscInt rank_hsl,PetscInt mpisize,char* tabfile, char *commsyntax,set_def *sets,dim_t nset, set_element *set_elems, array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar, elem_value **elem_vals2,offset_t ncofvar,offset_t ncofele,offset_t nvarele,closure_entry **closure_vals2,offset_t alltimeset,offset_t allregset,offset_t nintraeq,dim_t matsol,PetscInt Istart,PetscInt Iend,  offset_t nreg, offset_t ntime, offset_t *eq_addr, offset_t ndblock, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,dim_t laDi,dim_t laD,PetscReal cntl3,PetscReal cntl6,dim_t nesteddbbd,int localsize,PetscInt *ndbbddrank1,fortran_int* indata,dim_t mc66,fortran_int *ptx,struct timeval begintime,dim_t subints,MPI_Fint fcomm,int solmethod,solve_real **xcf2);
#endif // TEEMS_SOLVER_H_INCLUDED

