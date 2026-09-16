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
#include "version.h"

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
/* errmsg(fmt, ...): every "Error:" line goes through here (printf
   drop-in: prints to stdout, flushes so the line survives a crash that
   follows, and counts). The count is the exit-status backstop: main()
   returns nonzero when any error was printed, so a site that reports
   and continues can never end in exit 0 (finding 3 of the 2026-09
   platform gate). Fatal sites still terminate on the spot
   (MPI_Abort / PetscFinalize + return 1). */
extern int teems_error_count;
int errmsg(const char *fmt, ...) __attribute__((format(printf,1,2)));
extern int inmemory; /* -inmemory: keep value arrays resident instead of spilling to scratch */
extern int section_threads;
extern int ndbbd_threads_used[4]; /* NDBBD team sizes chosen by the thread budget: presolve, interface-rank, interface-factor, schur; 0 = region not run */
double teems_mem_avail_bytes(void); /* bytes the node can still give this process: MemAvailable capped by the cgroup limit; -1 if unreadable */
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
/* Coefficient storage precision. float halves memory traffic and is the
   production default; -DTEEMS_STORE_F64 (the teems-solver-f64 binary)
   stores double for deep-ladder accuracy — f32 update rounding
   anti-converges with pass count (ROADMAP 6.3(f3)). In-memory only:
   file formats and MPI messages carry solve_real or structs of it. */
#ifdef TEEMS_STORE_F64
typedef double store_real;
#define TEEMS_STORE_PRECISION "double"
#else
typedef float store_real;
#define TEEMS_STORE_PRECISION "single"
#endif

/* -matsol matrix method (Ha & Kompas 2016; Kompas & Ha 2019).
   Values match teems-R's matrix_method argument. */
enum matrix_method { MM_LU=0, MM_SBBD=1, MM_DBBD=2, MM_NDBBD=3 };
/* -solmed solution method (GEMPACK manual; Pearson 1991; Schiffmann
   2022 / GEMPACK 26.5 for the Runge-Kutta flavors) */
enum solution_method { SM_GRAGG=1, SM_EULER=2, SM_RK2=3, SM_RK4=4, SM_BOSHA32=5, SM_DOPRI54=6, SM_HEUN=7, SM_JOHANSEN=10, SM_PROBE=100 };
/* array_def.gltype: bound imposed on levels values */
enum bound_type { BT_NONE=0, BT_GE=1, BT_GT=2, BT_LE=3, BT_LT=4 };
/* formula_op.Oper: compiled formula operation */
enum op_code { OP_LOAD=0, OP_MUL=1, OP_DIV=2, OP_ADD=3, OP_SUB=4, OP_POW=5,
               OP_MAXF=61, OP_MINF=62, OP_ID0VF=63, /* multi-arg intrinsics
               (manual 11.5/11.5.1): pairwise folds over compiled temps */
               OP_IF_EQ=71, OP_IF_GT=72, OP_IF_LT=73, OP_IF_NE=74,
               OP_IF_LE=75, OP_IF_GE=76 };
/* formula_op operand types (Var1Type/Var2Type/Var3Type) */
enum operand_type { OT_ARRAY=0, OT_LINVAR=1, OT_SUM=2, OT_LINVAR2=3,
                    OT_TEMP=4, OT_CONST=5, OT_CHANGE=6,
                    OT_TEMP_ID01=41, OT_TEMP_ABS=42, OT_TEMP_LOG=43,
                    OT_TEMP_EXP=44, OT_TEMP_SQRT=45, OT_TEMP_LOG10=46,
                    OT_TEMP_ROUND=47, OT_TEMP_TRUNC0=48, OT_TEMP_TRUNCB=49,
                    /* $POS(...) (manual 11.5.6): an OP_LOAD whose value is
                       a quantifier index's 1-based position -- Var1BegAdd
                       = quantifier slot, Var1Dims[0].MapId routes through
                       a mapping first, Var1Dims[0].SupSet/SSIndx then lift
                       into a named superset (Var1Dims[0].ADims = the set
                       the position is lifted from) */
                    OT_POS=50 };


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
/* conditional set builders `Set X = (all,i,SRC: cond);` (manual
   10.1.2): data-dependent conditions evaluated straight from the
   input files at transform time, the statement rewritten into an
   explicit element list + subset relation; -1 on error, no-op when
   the TAB has none */
int tab_setbuilder_transform(char *fname, cmf_file_entry *iodata, int niodata);
int tab_preprocess(char *filename, char *newtabfile);
int tab_read_set_name(char *filename, char *varname, int indx, char *setname);
void tab_decl_index_free(void);

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
/* set product SET = A x B (manual 10.1.6), parallel to sets[] (the
   set_def layout is the sol.set dump read by teems-R, so it stays):
   isprod marks a definition that is exactly one product, prod1/prod2
   are the factor set ids with the first varying fastest -- Mapping
   (project) projects onto a factor (10.13.2) */
extern bool *teems_set_isprod;
extern dim_t *teems_set_prod1;
extern dim_t *teems_set_prod2;
/* one set element with its position in each superset */
typedef struct
{
  char setele[NAMESIZE];
  dim_t superset_pos[MAXSUPSET];
} set_element ;

/* a MAPPING declaration (manual 11.9): a total many-to-one function
   from domain-set elements to codomain-set positions. Values arrive
   via Read (by_elements) (11.9.1a) and are immutable once used
   (11.9.9). */
typedef struct
{
  char mapname[NAMESIZE];
  dim_t fromset;   /* set id of the domain */
  dim_t toset;     /* set id of the codomain */
  bool onto;       /* (onto): every codomain element must be hit (11.9.3) */
  bool has_values;
  bool used;
  dim_t *values;   /* per domain element: position in the codomain */
  /* value sources beyond Read (by_elements): a Mapping (project)
     qualifier (10.13.2) fills the table from the domain's product
     structure at declaration; formula_assigned marks a mapping that
     some Formula assigns (10.13.1/11.9.1b) -- validated complete at
     first use instead of up front, with `assigned` tracking domain
     positions written so far */
  bool project;
  bool formula_assigned;
  unsigned char *assigned;
  dim_t nassigned;
} map_def ;

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
  int cond_mapid;              /* >0: mapping-equality condition on the summed index (11.4.11, M3) */
  char cond_rhs[NAMESIZE];     /* its RHS token: outer quantifier index or codomain element */
  /* coefficient-comparison condition (11.4.11; IF-survey gap 2):
     COEF(args) <op> <numeric const>, e.g. ENDOWFLAG(e,t) NE 0 */
  char cond_coef[NAMESIZE];    /* "" = none */
  char cond_cofargs[MAXVARDIM][NAMESIZE];
  dim_t cond_cofnargs;
  int cond_cofop;              /* 1..6 = eq/ne/gt/lt/ge/le */
  double cond_cofval;
} sum_def ;

/* a resolved coefficient-comparison sum condition: per-tuple value
   lookup bound to the evaluation frame + the summed index */
typedef struct
{
  offset_t cofid;              /* -1 = no coefficient condition */
  offset_t offset;             /* coefs[cofid].offset */
  dim_t nd;
  int bind[MAXVARDIM];         /* -2 the summed index; >=0 frame slot; -1 fixed */
  offset_t fix[MAXVARDIM];     /* fixed element position (bind -1) */
  offset_t strides[MAXVARDIM];
  int op;
  double cval;
} sum_cofcond ;

typedef struct
{
  store_real value;
} elem_store ;


/* closure and shock for one variable element.  exo_index is the matrix
   column for exogenous/endogenous elements; for backsolved elements it
   is the compact index into the per-step recovered-value array.  The
   closure flags and the shock value live in two side arrays indexed by
   the same element offset (6.16(b) struct diet: 4 + 1 + sizeof
   (store_real) bytes per element instead of a padded 12/16-byte
   struct); both are allocated with closure_vals, broadcast with it
   under nohsl, and stay resident when closure_vals spills to scratch.
   Access only through the CL_* macros. */
typedef struct
{
  exo_idx_t exo_index;
} closure_entry ;
#define CL_F_EXO ((unsigned char)1)
#define CL_F_BS  ((unsigned char)2)
extern unsigned char *teems_cl_flags; /* per element: CL_F_EXO | CL_F_BS */
extern store_real *teems_cl_shock;    /* per element: shock value (0 = none) */
#define CL_EXO(i) ((teems_cl_flags[(i)]&CL_F_EXO)!=0)
#define CL_BS(i)  ((teems_cl_flags[(i)]&CL_F_BS)!=0)
#define CL_SET_EXO(i,v) do{ if(v)teems_cl_flags[(i)]|=CL_F_EXO; else teems_cl_flags[(i)]&=(unsigned char)~CL_F_EXO; }while(0)
#define CL_SET_BS(i,v)  do{ if(v)teems_cl_flags[(i)]|=CL_F_BS;  else teems_cl_flags[(i)]&=(unsigned char)~CL_F_BS;  }while(0)
#define CL_SHOCK(i) (teems_cl_shock[(i)])

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
/* -assertions run switch (manual 25.3): 0 off, 1 warn, 2 fatal (default) */
extern int teems_assertions_mode;
/* declared-range checks (manual 25.4.4): 0 no, 1 warn (default), 2
   fatal; initial leg = formulas passes with IsIni, updated leg = the
   update executors and later formulas passes */
extern int teems_range_test_initial;
extern int teems_range_test_updated;
/* Runge-Kutta run controls (main parses -rkchart/-rknorm/-rkctrl/
   -rk_h0/-rkguard; solve_rk.c reads them) */
enum rk_chart_kind { RK_CHART_LOG=0, RK_CHART_PERCENT=1 };
enum rk_norm_kind { RK_NORM_MAX=0, RK_NORM_RMS=1 };
enum rk_ctrl_kind { RK_CTRL_STD=0, RK_CTRL_PI=1 };
enum rk_scope_kind { RK_SCOPE_PCT=0, RK_SCOPE_ALL=1 };
typedef struct {
  int chart;        /* rk_chart_kind */
  int norm;         /* rk_norm_kind: error-metric norm for the accept test */
  int ctrl;         /* rk_ctrl_kind: step-size controller */
  int scope;        /* rk_scope_kind: which elements steer the accept test */
  double h0;        /* initial step length (0 = 1/step1 capped by the initial gradient) */
  double guard;     /* log chart: |log(level/level0)| beyond which a stage is rejected */
} rk_options;
/* Runge-Kutta run record (solve_rk.c fills, main patches into stats.json) */
typedef struct {
  long stage_solves, stage_solves_reused, steps;
  long rejects_accuracy, rejects_crossed, rejects_range, rejects_assert, rejects_guard, rejects_singular;
  double h_min, h_max, worst_step_metric, worst_est_metric;
  int chart;
} teems_rk_stats_t;
extern teems_rk_stats_t teems_rk_stats;
/* set while a Runge-Kutta stage state is being realized under -adaptive
   yes: formula.c counts range-test and assertion violations into the
   two counters (and prints them as warnings) instead of aborting, so
   the driver can retry the step (manual 26.5.1 triggers 2 and 3) */
extern int teems_rk_stage_checks;
extern long teems_check_viol_range;
extern long teems_check_viol_assert;
/* nonzero while an LU stage solve may fail softly (a singular stage
   state under adaptive Runge-Kutta): the kernel returns and sets
   teems_stage_solve_failed instead of aborting */
extern int teems_rk_softfail;
extern int teems_stage_solve_failed;
/* (parameter)-qualified coefficients, parallel to coefs[] (F2) */
extern bool *teems_coef_is_param;
/* second range bound (one lower + one upper per declaration, manual
   10.19.1; audit A9), parallel to coefs[] -- slot 1 stays in the
   binary-locked array_def */
extern int *teems_coef_gltype2;
extern store_real *teems_coef_glval2;
/* satisfied Read (IfHeaderExists) targets, parallel to coefs[] (3.9) */
extern bool *teems_coef_ifhdr;
/* set mappings (manual 11.9), populated in main after the broadcast;
   consumed by the operand binder/eval and the statement guards */
extern map_def *teems_maps;
extern dim_t teems_nmap;
/* the set table, for name lookups inside the formula compiler ($POS
   set arguments) -- set once the elements are built */
extern set_def *teems_sets;
extern dim_t teems_nset;
extern set_element *teems_set_elems;
/* mark mappings assigned by Formula statements in fname (main TAB or
   the PostSim split) so mappings_validate defers their completeness
   check to first use */
void mapping_formula_scan(char *fname, map_def *maps, dim_t nmap);
int mapping_check_onto(map_def *maps, dim_t j, set_def *sets, set_element *set_elems);
void set_expr_mark_product(char *buf);
dim_t set_expr_bound(char **pp, set_def *record, dim_t nset, const char *owner, int *err);
/* GEMPACK dual-class zerodivide state (manual 10.11; plan A1): tracked
   positionally by the statement scanner, consulted by formula
   evaluation only under -gpzerodivide 1 (default 0 = the legacy single
   conflated default, bit-compatible). Initial GEMPACK state: 0/0 -> 0,
   nonzero/0 -> error. */
typedef struct {
  solve_real zbz_val;
  solve_real nbz_val;
  int zbz_on;
  int nbz_on;
} zdiv_state ;
extern zdiv_state teems_zdiv_scan;
extern int teems_gpzerodivide;
extern long teems_laA_used,teems_laDi_used,teems_laD_used; /* max grown -la* equivalent percent (la auto-sizing) */
/* MA48 workspace ceiling.  The HSL kernels are built with 32-bit
   integers -- INSIZE and MA48's own LA are integer(4) -- so a workspace
   request cannot exceed INT_MAX elements.  Past that the request wraps:
   MA48 is handed a negative or far-too-small LA while the caller has
   allocated the full 64-bit size, which is how an undersized run used
   to fail non-monotonically (a bigger retry could fail where a smaller
   one had succeeded) and then die in the allocator with no diagnosis. */
#define MA48_LA_MAX 2147483647L
/* Grow an MA48 workspace after a -3 return: applies MA48's suggestion
   with a doubling floor, clamps to MA48_LA_MAX, emits the growth note,
   records the -la* equivalent in *used, and aborts with a named
   diagnosis once the ceiling leaves nothing larger to ask for. */
offset_t ma48_grow_la(offset_t cur,offset_t suggested,offset_t nnz,const char *knob,long *used);
offset_t ma48_la_from_pct(dim_t pct,offset_t nnz); /* initial LA from a -la* percent, clamped to MA48_LA_MAX */
/* checked allocation for MA48 workspaces: a clamped LA still reaches
   tens of GB, and an unchecked NULL here is a segfault rather than a
   report */
void *ma48_alloc(offset_t n,size_t sz);
void *ma48_realloc(void *p,offset_t n,size_t sz);
/* Jacobian/exogenous-block preallocation with the PetscInt ceiling
   named: PETSc sums the row counts into a PetscInt, and the one-rank
   whole-system copy every HSL matrix_method keeps passes that near
   2^31 nonzeros (unchecked, the still-unallocated Mat segfaulted at
   the first insert).  count=1 on the rank that owns the row counts. */
void jac_mat_prealloc(Mat M,const char *what,PetscBool mpi,int count,PetscInt nrows_local,PetscInt dnz,PetscInt *dnnz,PetscInt onz,PetscInt *onnz);
/* -condest solve-quality diagnostics (MA60/MC71, sequential LU path);
   scope gate set by the LU wrappers, accumulators reduced into
   stats.json post-solve.  teems_condest_active_/report_ are the
   Fortran-visible entries used by SPEC48_SSOL2LA(_P). */
/* Phase resident-memory record (ROADMAP 6.16(a)): after each phase,
   step and solve line every rank reads its resident set
   (/proc/self/statm) and high-water mark (getrusage ru_maxrss); the
   values are reduced to rank 0 as max and sum over ranks, printed as a
   "memory:" line and kept per phase for stats.json ("rss_gb").
   Collective on PETSC_COMM_WORLD: call it on every rank.  A phase
   probed repeatedly (steps) keeps its maximum. */
/* Wall-clock stage record for the bordered drivers (NDBBD first):
   teems_stage_mark(name) closes the running stage and opens `name`
   (NULL closes only); teems_stage_report(prefix) reduces every stage's
   wall to its max over ranks and prints one line on rank 0, then
   resets.  Marks must sit at points every rank passes in the same
   order (the report is a collective); inside OpenMP regions only
   thread 0 records. */
/* NDBBD cut cache (-ndcutcache, default 1): the regional cut
   (ndbbd_order_presolve: rank-deficient block tails migrated to the
   time-interface blocks, counteq/countvarintra1/block_sizes, the
   row/col permutations, ndbbdrank) and the per-chain-block interface
   cut (ndbbd_presolve's rank probe: rank + row/col permutation) are
   fixed by the model's structure, yet were recomputed by a throwaway
   MA48 factorization of every block at every step (and every RK
   stage).  Step 1 computes them as before and stores them; later
   steps reuse them.  Keyed on VecSize/ndblock; freed on the
   complementarity re-entry (closure changes the partition). */
extern int teems_ndcutcache;
void ndbbd_cut_cache_free(void);
void ndbbd_cut_iface_init(int nmatint);
int ndbbd_cut_iface_get(int j3,int *rank_out,int *irn,int *jcn,int nrow,int ncol);
void ndbbd_cut_iface_put(int j3,int rank_val,const int *irn,const int *jcn,int nrow,int ncol);
#define TEEMS_STAGE_MAX 24
void teems_stage_mark(const char *name);
void teems_stage_report(const char *prefix);
#define TEEMS_RSS_MAX 16
typedef struct {
  const char *phase;
  double rss_max,rss_sum,hwm_max,hwm_sum;
  long count;
} teems_rss_entry;
extern teems_rss_entry teems_rss[TEEMS_RSS_MAX];
extern int teems_nrss;
void teems_rss_probe(const char *phase);
extern int teems_condest,teems_condest_scope;
extern double teems_condest_kw1max,teems_condest_kw2max,teems_condest_omegamax;
extern long teems_condest_solves,teems_condest_skips;
int teems_condest_active_(void);
void teems_condest_report_(int *status,double *omega1,double *omega2,double *erx,double *cond1,double *cond2,int *noiter);
void zdiv_scan_reset(void);
void zdiv_capture(void);
void zdiv_disable(void);
/* evaluate ASSERTION statements against the current values (manual
   10.14/25.3): rides each formulas_execute pass; (initial) assertions
   only when IsIni. postsim_pass: 0 = ordinary passes ((postsim)
   assertions skipped), 1 = the post-solve pass (only (postsim)
   assertions, initial/always ignored per manual 12.2.4). Returns the
   failure count (mode 1). */
offset_t assertions_execute(char *fname,set_def *sets,dim_t nset,set_element *set_elems,array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar,elem_value *elem_vals,offset_t ncofvar,offset_t ncofele,bool IsIni,int mode,int postsim_pass);
/* PostSim foundation F3: copy the composed solution (xcf, per
   variable element) into the variables' elem_vals slots so post-solve
   statements read simulation results as if variables were
   coefficients (manual ch.12) */
void postsim_expose_results(elem_value *elem_vals,offset_t ncofele,offset_t nvarele,solve_real *xcf);
int tab_has_postsim_assertions(char *fname);
/* Tier 0: split POSTSIM sections out of the preprocessed TAB (returns
   PostSim executable count; -1 on error); the -postsim run switch is
   parsed in main */
int tab_postsim_split(char *newtabfile, char *psfile);
/* C0: levels-statement transform -- expand Formula&Equation, pair
   levels variables with value coefficients + updates, linearize
   Equation (levels) by change differentiation (design doc
   mapping_complementarity_design.md section 5); no-op when the TAB
   has no levels statements; -1 on error */
int tab_levels_transform(char *fname);
/* C1: one Complementarity statement (manual 10.17/11.14; design doc
   section 7): parsed at transform time, set matching validated after
   set elements exist, closure integration after closure_read. Bound
   kinds: 0 none, 1 real constant, 2 levels variable, 3 parameter
   coefficient. Quantifier/argument set names are recorded textually
   by the transform (declarations are not yet read at that stage). */
typedef struct
{
  char name[NAMESIZE];              /* <= 10 chars (11.2.1) */
  char varname[NAMESIZE];           /* X: the complementarity variable */
  int lower_kind, upper_kind;
  double lower_const, upper_const;
  char lower_name[NAMESIZE], upper_name[NAMESIZE];
  dim_t nquant;
  char qidx[MAXVARDIM][NAMESIZE];   /* statement quantifier indices */
  char qset[MAXVARDIM][NAMESIZE];   /* statement quantifier sets */
  char xset[MAXVARDIM][NAMESIZE];   /* X's declared set per arg position */
  char lset[MAXVARDIM][NAMESIZE];   /* bound sets (kind 2/3 only) */
  char uset[MAXVARDIM][NAMESIZE];
  /* C2 state machinery: value-side names of the levels pairs (the
     pair coefficient name differs from the declared name for
     gen_lvN-renamed p_-leading pairs) and whether the linear pair is
     percent (column weight Xval/100) or change (weight 1). For a
     parameter bound the "value name" is the coefficient itself. */
  char xval[NAMESIZE], lval[NAMESIZE], uval[NAMESIZE];
  int xpct, lpct, upct;
} comp_def ;
extern comp_def *teems_comps;
extern dim_t teems_ncomp;
/* count of complementarity-variable components left ENDOGENOUS by the
   closure = components solved by the C2 approximate-run state
   machinery (exogenous components stay inert: their dummy comp@d is
   endogenous and absorbs the E_$comp row). Set by comp_closure_check
   on the closure-reading rank, broadcast in main. */
extern offset_t teems_comp_active;
int tab_complementarity_transform(char *fname);
int complementarities_validate(set_def *sets, dim_t nset, set_element *set_elems);
int comp_closure_check(closure_entry *closure_vals, array_def *vars, offset_t nvar, offset_t *nexo, set_def *sets, dim_t nset, set_element *set_elems);
/* C2 per-step state machinery (design doc section 8; manual 51.1.2/
   51.2/51.7.3/51.7.5). All three run on rank_hsl only (levels values
   are updated there); the driver broadcasts redo decisions.
   comp_states_set: evaluate the per-component state from the current
   levels values, write the E_$comp weight coefficients and store the
   step-start margins (lazy init on first call: resolve ids, pre-sim
   51.7.5 exactness check). comp_states_check: recompute states from
   the post-step values; returns the flip count and the smallest
   linear-interpolation crossing fraction. comp_states_report:
   post-sim 51.7.5 check + 51.5.3-style state-change log lines.
   comp_states_free: teardown. Return -1 on resolution failure. */
int comp_states_set(set_def *sets, dim_t nset, set_element *set_elems, array_def *coefs, offset_t ncof, array_def *vars, offset_t nvar, elem_value *elem_vals);
int comp_states_check(set_def *sets, dim_t nset, set_element *set_elems, array_def *coefs, offset_t ncof, array_def *vars, offset_t nvar, elem_value *elem_vals, offset_t *nflip, double *minfrac);
int comp_states_report(set_def *sets, dim_t nset, set_element *set_elems, array_def *coefs, offset_t ncof, array_def *vars, offset_t nvar, elem_value *elem_vals);
void comp_states_free(void);
/* C3 accurate run (manual 51.7.1/51.5.4; design doc section 8 tail).
   comp_accurate_prepare: capture the per-component target states
   from the CURRENT values (post-approximate normally; pre-sim under
   do_approx_run = no, with lazy init). comp_accurate_closure: the
   51.7.1 closure/shock modification, run on PRE-SIM-restored values
   (dummy endogenous + one component exogenized per active component;
   nexo unchanged). comp_verify_states: post-accurate 51.5.4/51.7.5
   checks; returns the violation count (caller maps warn/fatal). */
int comp_accurate_prepare(set_def *sets, dim_t nset, set_element *set_elems, array_def *coefs, offset_t ncof, array_def *vars, offset_t nvar, elem_value *elem_vals);
int comp_accurate_closure(closure_entry *closure_vals, array_def *vars, offset_t nvar, array_def *coefs, offset_t ncof, set_def *sets, dim_t nset, set_element *set_elems, elem_value *elem_vals);
offset_t comp_verify_states(set_def *sets, dim_t nset, set_element *set_elems, array_def *coefs, offset_t ncof, array_def *vars, offset_t nvar, elem_value *elem_vals);
/* the 51.6 run controls arrive as command-line flags (-comp_steps,
   -comp_redo, -comp_redo_min_frac, -comp_do_approx, -comp_do_acc,
   -comp_sberr_warn), passed by teems-R's ems_complementarity();
   there are no complementarity CMF statements */
/* LinVar token resolution incl. declared p_-/c_-leading names
   (design doc section 6); -1 when nothing matches */
offset_t linvar_resolve(char *vname, array_def *vars, offset_t nvar);
offset_t backsolve_read(char *fname, array_def *vars, offset_t nvar, closure_entry *closure_vals);
int backsolve_validate_refs(char *fname, array_def *vars);
int tab_equation_name(char *stmt, char *eqname);
int str_subst_all_bounded(char *line, const char *finditem, const char *replitem, size_t linesz);
void str_delete_char(char *s, char c); /* remove every occurrence of c in place, one pass */ /* forward-scanning replace-all within a buffer of linesz; -1 if the result would not fit */

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
  int dimmapid[MAXVARDIM];   /* >0: dim routes via teems_maps[id-1] (11.9.5, M2b) */
  int dimcondmap[MAXVARDIM]; /* >0: dim's enclosing sum has a mapping-equality condition (M3) */
  char dimcondrhs[MAXVARDIM][NAMESIZE]; /* its RHS token */
} eq_var_ref ;

/* per-equation-statement addressing metadata, captured by
   jacobian_preallocate (only when asked: eqmeta!=NULL, i.e. -solmed
   probe) so the structural diagnosis can invert a matrix row back to
   its named equation element: row -> eq_addr position -> statement +
   quantifier tuple.  Quantifier dims are in declared order, first
   slowest (the dcountdim1 layout).  var_ref/var_w record the
   statement's linear-variable references with their element-level
   incidence weights (entries written per reference, pre-merge) for
   the probe report's statement-level structure section. */
#define PROBE_MAXEQVARS 64
typedef struct
{
  char eqname[NAMESIZE];
  dim_t fdim;                     /* number of (all,) quantifiers */
  offset_t setid[4*MAXVARDIM];    /* quantifier set ids, declared order */
  offset_t base;                  /* first row element (matroworg) */
  offset_t nrows;                 /* quantifier-space size (nloops) */
  dim_t nvars_ref;                /* distinct variables referenced */
  offset_t var_ref[PROBE_MAXEQVARS]; /* vars[] indices */
  offset_t var_w[PROBE_MAXEQVARS];   /* incidence weights */
} eq_probe_meta ;

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
/* superset slot (0 = same set, >0 = superset_pos column, -1 = not a
   declared subset) and the named fatal for the -1 case */
dim_t set_supset_slot(set_def *sets, dim_t sub, dim_t sup);
void set_supset_fatal(const char *idx, const char *symname, const char *where, set_def *sets, dim_t sub, dim_t sup);
char *closure_next_statement(char *commsyntax, FILE *filehandle, char *readline);
char *tab_next_statement(char *commsyntax, FILE *filehandle, char *readline,offset_t rlinesize);
char *tab_next_statement_resolved(char *commsyntax, FILE *filehandle, char *readline, elem_value *record, array_def *coefs,offset_t ncof,solve_real *zerodivide,offset_t rlinesize);
int str_find_ci(char *line, char *finditem);
int str_find_token_ci(const char *base, const char *s, const char *pat);
int str_count_token_ci(const char *base, const char *s, const char *pat);
char *str_replace_all(char *line, char *finditem, char *replitem);
int str_replace_char_all(char *line, int finditem, int replitem);
char *str_replace_all_bounded(char *line, char *finditem, char *replitem,dim_t nbuffer);
char *str_replace_first_bounded(char *line, char *finditem, char *replitem,dim_t nbuffer);
char *str_replace_first(char *line, char *finditem, char *replitem);
char *str_strip_comment(char *line, char *token);
/* PostSim scope (Tier 0 residuals; manual 12.2.1-12.2.3) */
extern char (*teems_ps_coefnames)[NAMESIZE];
extern int teems_ps_ncoefs;
extern bool *teems_coef_is_ps;
extern int teems_ps_pass;
void postsim_mark_coefs(array_def *coefs, offset_t ncof);
/* 11.2.1 name uniqueness across coefficient/variable/set/mapping +
   reserved words (the 12.2.2 name-resolution spec pass) */
int names_validate(set_def *sets, dim_t nset, array_def *coefs, offset_t ncof, array_def *vars, offset_t nvar, map_def *maps, dim_t nmap);
/* MAPPING statements (manual 11.9): declarations, by_elements value
   reads, and the pre-use validation pass (range/onto/coverage) */
int mappings_read(char *fname, map_def *maps, dim_t nmap, set_def *sets, dim_t nset);
int mapping_values_read(char *fname, int niodata, cmf_file_entry *iodata, map_def *maps, dim_t nmap, set_def *sets, set_element *set_elems);
int mapping_use_guards(char *fname, map_def *maps, dim_t nmap);
void mapping_lower_calls(char *line);
void mapping_reject_in(char *line, const char *what);
void mapping_reject_lhs(char *line, const char *what);
char *mapping_token_split(char *p, int *mp);
char *sum_dim_identity(char *p);
char *sum_settok_extract(const char *sumtext);
char *sum_body_extract(const char *sumtext);
void sum_cond_parse(char *settok, const char *sumindx, int *cond_mapid, char *cond_rhs, sum_def *sc);
/* resolve/evaluate a coefficient-comparison sum condition against an
   evaluation frame (fatal on contract violations; no-op when the sum
   has no coefficient condition) */
void sum_cond_coef_resolve(sum_def *sc, quantifier *frame, dim_t nframe, set_def *sets, set_element *set_elems, array_def *coefs, offset_t ncof, sum_cofcond *out);
int sum_cofcond_test(const sum_cofcond *cc, elem_value *elem_vals, quantifier *frame, offset_t l1);
void sum_cond_domain_check(sum_def *sc, set_def *sets);
dim_t sum_cond_carry_rhs(sum_def *sc, quantifier *arSet, dim_t fdim, dim_t l3, char *interchar, set_def *sets);
void sum_cond_rhs_resolve(int cond_mapid, const char *cond_rhs, quantifier *frame, dim_t nframe, set_def *sets, set_element *set_elems, int *condpos, offset_t *condfix, dim_t *condss);
/* codomain position a mapping-equality condition compares against:
   the fixed element, or the RHS quantifier's current position lifted
   from its (sub)set into the codomain via superset_pos slot condss
   (0 = the quantifier ranges over the codomain itself) */
static inline offset_t sum_cond_target(int condpos, dim_t condss, offset_t condfix, const quantifier *frame, const set_def *sets, const set_element *set_elems) {
  if (condpos<0) return condfix;
  if (condss==0) return (offset_t)frame[condpos].indx;
  return (offset_t)set_elems[sets[frame[condpos].setid].offset+frame[condpos].indx].superset_pos[condss];
}

/* Read exactly n items from a solver scratch file. A short read means the
   file the earlier pass wrote is truncated or missing (scratch space
   exhausted, or the file removed under the run): abort with the file name
   rather than continue on whatever the buffer held. */
static inline void scratch_read(void *buf,size_t size,size_t n,FILE *fp,const char *fname) {
  if(fp==NULL||fread(buf,size,n,fp)!=n) {
    printf("Error: short read on scratch file %s; the pass that writes it did not complete (check free space under -tempdir)\n",fname);
    fflush(stdout);
    MPI_Abort(PETSC_COMM_WORLD,1);
  }
}
/* rewrite the word comparison operators (ge le gt lt ne eq, space-
   delimited, outside quotes) to their symbol forms ahead of
   whitespace stripping */
void tab_wordops_normalize(char *line);
int mappings_validate(map_def *maps, dim_t nmap, set_def *sets, set_element *set_elems);
offset_t postsim_reads_execute(char *psname, int niodata, cmf_file_entry *iodata, set_def *sets, dim_t nset, set_element *set_elems, array_def *coefs, offset_t ncof, offset_t ncofele, array_def *vars, offset_t nvar, offset_t nvarele, elem_value *elem_vals);
/* Default-statement helpers (manual 10.19; audit A6): positional
   semantics live in the readers; values are validated once up front */
int tab_default_value(char *line, char *out);
int tab_defaults_validate(char *fname);
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
  int MapId;                 /* >0: route via teems_maps[MapId-1] (11.9.4) */
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
int jacobian_fill(char *fname, char *commsyntax,set_def *sets,offset_t nset, set_element *set_elems, array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar, elem_value *elem_vals,offset_t ncofvar,offset_t ncofele,closure_entry *closure_vals,offset_t ndblock,offset_t alltimeset,offset_t allregset,PetscInt *eq_addr,offset_t *counteq,offset_t nintraeq,Mat A,Mat B);
void jacobian_cache_free(void); /* release the per-rank compiled-statement cache */
/* -fastrefac sequential-LU persistent refactorize (solve_drivers.c): the
   COO pattern and MA48 pivot sequence persist across steps/stages, later
   factorizations run MA48B/BD JOB=2 */
void lu_fastrefac_solve(Mat A,PetscInt VecSize,dim_t laA,solve_real *rhs,solve_real *x);
void lu_fastrefac_free(void); /* release persistent LU state after the solve dispatch */
void lu_grow_solve(Mat A,PetscInt VecSize,dim_t laA,solve_real *rhs,solve_real *x); /* one-shot LU w/ MA48 workspace grow-and-retry */
/* -fastrefac SBBD persistent MP48 instance (solve_drivers.c): border lists,
   per-block pivot sequences and factors persist across steps; repeat steps
   refactorize with FACT_JOB=2.  Collective — all ranks call both. */
void sbbd_fastrefac_solve(Mat *A,Vec *vecb,PetscInt VecSize,PetscInt rank,PetscInt rank_hsl,fortran_int *indata,MPI_Fint fcomm,offset_t *counteq,offset_t *countvarintra1,solve_real *x);
void sbbd_fastrefac_free(void);
/* SBBD one-shot solve staged straight from PETSc's CSR into MP48's host
   arrays (6.15(c)); destroys *A and *vecb before the factorization.
   Collective. */
void sbbd_csr_solve(Mat *A,Vec *vecb,PetscInt VecSize,PetscInt rank,PetscInt rank_hsl,fortran_int *indata,MPI_Fint fcomm,offset_t *counteq,offset_t *countvarintra1,solve_real *x);
/* -fastrefac DBBD per-block persistent factors (block_solve.c): the flag is
   read inside dbbd_solve, so all drivers inherit it */
void dbbd_fastrefac_free(void);
void dbbd_fastextract_free(void); /* persistent DBBD submatrix extraction (-fastrefac) */
/* -fastrefac NDBBD regional-block persistent factors (block_solve.c):
   flag read inside the ndbbd paths; frees the whole ndbbd_fac store */
void ndbbd_fastrefac_free(void);
/* Recover the backsolved variables' per-step values from their retained
   defining equations (GEMPACK 14.1.3: after the condensed solve, before
   the data updates).  x = this step's solution vector; exo_z = per-element
   exogenous per-step changes captured at the vece fill sites; bsvals
   (nbselems) receives the recovered changes at each pair's elem_base.
   Returns -1 on a zero pivot or a malformed defining equation. */
int backsolve_recover(char *fname, char *commsyntax,set_def *sets,offset_t nset, set_element *set_elems, array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar, elem_value *elem_vals,offset_t ncofele,closure_entry *closure_vals,solve_real *x,solve_real *exo_z,solve_real *bsvals);
void backsolve_cache_free(void);
int eq_linvar_read(char *formulain,eq_var_ref *LinVars,int linindx,array_def *vars);
/* Istart/Iend delimit this rank's rows of A (and, A being square with a
   matching layout, its diagonal-block columns); Cstart/Cend delimit the
   exogenous columns of B, which follow the shock vector's layout.  The
   two ranges coincide unless condensation left nexo > VecSize, where B
   is wider than it is tall and carries its own column split. */
int jacobian_preallocate(char *fname, char *commsyntax,set_def *sets,dim_t nset,set_element *set_elems,array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar,elem_value *elem_vals,offset_t ncofvar,offset_t ncofele, offset_t nexo,closure_entry *closure_vals,offset_t ndblock,offset_t alltimeset,offset_t allregset,bool *eq_intertemp,PetscInt *eq_addr,dim_t *eq_time,dim_t *eq_reg,offset_t *counteq,offset_t nintraeq,bool *sbbd_overrid,PetscInt VecSize,PetscInt Istart,PetscInt Iend,PetscInt Cstart,PetscInt Cend,PetscInt *dnz,PetscInt *dnnz,PetscInt *onz,PetscInt *onnz,PetscInt *dnzB,PetscInt *dnnzB,PetscInt *onzB,PetscInt *onnzB,int nesteddbbd,eq_probe_meta *eqmeta,offset_t *neqmeta);
/* Exogenous-side layout (the shock vector and B's columns).  Normally
   the exogenous columns mirror the equation rows and, under NDBBD,
   follow the same time-block split; when heavy condensation leaves more
   exogenous elements than unknowns the columns span BSize instead and
   take PETSc's own split, which the row blocks cannot express. */
void shock_vec_set_sizes(Vec v,int nesteddbbd,PetscInt localsize,PetscInt VecSize,PetscInt BSize);
void shock_mat_set_sizes(Mat B,int nesteddbbd,PetscInt localsize,PetscInt VecSize,PetscInt BSize);
/* -solmed probe: assemble the condensed Jacobian sequentially and run
   the HSL_MC79 maximum-matching / Dulmage-Mendelsohn structural
   diagnosis on it (full stored pattern + numerically realized
   pattern), naming defective variable and equation elements. */
int probe_structural(PetscInt VecSize,offset_t nvarele,offset_t ncofele,PetscInt dnz,PetscInt *dnnz,PetscInt dnzB,PetscInt *dnnzB,char *tabfile,char *commsyntax,set_def *sets,dim_t nset,set_element *set_elems,array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar,elem_value *elem_vals,closure_entry *closure_vals,offset_t ndblock,offset_t alltimeset,offset_t allregset,PetscInt *eq_addr,offset_t *counteq,offset_t nintraeq,eq_probe_meta *eqmeta,offset_t neqmeta,cmf_file_entry *iodata,int niodata,int noutdata,int nsoldata,int probefine,PetscInt mpisize,PetscInt rank);
/* on-failure diagnosis (part 3): the factorization kernels call
   teems_onfail_diag_() before their failure STOPs; solve paths
   register the naming context (once, from main) and the system about
   to be factorized (live Mat + local->condensed row/col maps). */
void probe_onfail_context(set_def *sets,set_element *set_elems,array_def *vars,offset_t nvar,closure_entry *closure_vals,offset_t nvarele,PetscInt *eq_addr,eq_probe_meta *eqmeta,offset_t neqmeta,PetscInt VecSize);
void probe_onfail_scope_set(Mat A,PetscInt m,PetscInt n,const char *label,int block_id,int *row_order,int *col_order,offset_t row_base,offset_t col_base,offset_t row_add,offset_t col_add);
void probe_onfail_scope_set_coo(const int *irn,const int *jcn,const solve_real *va,const int *rowlen,long nz,PetscInt m,PetscInt n,const char *label,int *row_map,int *col_map);
void probe_onfail_scope_set_csr(const int *rowptr,const int *jcn,const solve_real *va,long nz,PetscInt m,PetscInt n,const char *label); /* 1-based CSR staging (MP48 host arrays) */
void probe_onfail_scope_clear(void);
const char *probe_onfail_scope_label(void); /* label of the system registered for on-failure diagnosis, for errors raised outside probe.c */
void teems_onfail_diag_(int *info1);
void teems_onfail_abort_(void);
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

int ndbbd_order_presolve(Mat A, offset_t VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int nreg, int ntime, offset_t nvarele, PetscInt *eq_addr,int *row_order,int *col_order, offset_t ndblock,int *block_sizes, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,dim_t laDi,solve_real cntl6,PetscInt* ndbbdrank,PetscBool presol);
int dbbd_order(Mat A, offset_t VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend, offset_t nvarele, PetscInt *eq_addr,int *row_order,int *col_order, offset_t ndblock,int *block_sizes, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,solve_real cntl6);
int dbbd_solve(Mat A, Vec b, solve_real *x1, offset_t VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int *row_order,int *col_order, offset_t ndblock,int *block_sizes, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,dim_t laD,PetscReal cntl3);//,bool iter
int ndbbd_order(Mat A, offset_t VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int nreg, int ntime, offset_t nvarele, PetscInt *eq_addr,int *row_order,int *col_order, offset_t ndblock,int *block_sizes, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,dim_t laDi,solve_real cntl6,PetscInt* ndbbdrank,PetscBool presol);
int ndbbd_presolve(Mat A, Vec b, solve_real *x1, offset_t VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int *row_order,int *col_order, offset_t ndblock,offset_t nreg,offset_t ntime,int *block_sizes, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,dim_t laDi,dim_t laD,PetscReal cntl3,PetscReal cntl6,PetscBool presol);//,bool iter
int ndbbd_solve(Mat A, Vec b, solve_real *x1, offset_t VecSize, PetscInt mpisize, PetscInt rank, PetscInt Istart, PetscInt Iend,int *row_order,int *col_order, offset_t ndblock,offset_t nreg,offset_t ntime,int *block_sizes, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,dim_t laDi,dim_t laD,PetscReal cntl3,PetscReal cntl6,PetscBool presol);//,bool iter
int reduce_to_rank(solve_real *vecbivi,fortran_int vecbivisize,PetscInt mpisize,PetscInt rank,PetscInt targetrank);
int reduce_to_rank_nocompress(solve_real *vecbivi,fortran_int vecbivisize,PetscInt mpisize,PetscInt rank,PetscInt targetrank);
int outputs_write_csv(char *filename, char *newdatlogname, char *newdatfile,set_def *sets,dim_t nset, set_element *set_elems,array_def *coefs,offset_t ncof,offset_t ncofele,array_def *vars,offset_t nvar,offset_t nvarele, elem_value *elem_vals);

/* ============ solve_drivers.c — solution methods ======================= */

/* one-step (Johansen 1960) solution of the linearized system */
bool solve_johansen(PetscBool nohsl,PetscInt VecSize,Mat A,PetscInt dnz,PetscInt* dnnz,PetscInt onz,PetscInt* onnz,Mat B,PetscInt dnzB,PetscInt* dnnzB,PetscInt onzB,PetscInt* onnzB,Vec vecb,Vec vece,PetscInt rank,PetscInt rank_hsl,PetscInt mpisize,char* tabfile, char *commsyntax,set_def *sets,dim_t nset, set_element *set_elems, array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar, elem_value **elem_vals2,offset_t ncofvar,offset_t ncofele,offset_t nvarele,closure_entry **closure_vals2,offset_t alltimeset,offset_t allregset,offset_t nintraeq,dim_t matsol,PetscInt Istart,PetscInt Iend,  offset_t nreg, offset_t ntime, PetscInt *eq_addr, offset_t ndblock, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,dim_t laDi,dim_t laD,PetscReal cntl3,PetscReal cntl6,dim_t nesteddbbd,int localsize,PetscInt *ndbbddrank1,fortran_int* indata,dim_t mc66,fortran_int *ptx,struct timeval begintime,solve_real **xcf2);
/* Multistep driver with Richardson extrapolation over
   steps1/steps2/steps3. solmethod selects the stepping scheme:
   SM_GRAGG (smoothed modified midpoint, Pearson 1991 eq. 6.1 /
   Alg. 7.1.2): Euler start, midpoint leapfrog, terminal smoothing
   pass, h^2 error series; SM_EULER: forward Euler on every substep,
   no smoothing pass, h error series (extrapolation weights differ
   accordingly). */
bool solve_gragg(PetscBool nohsl,PetscInt VecSize,Mat* A,PetscInt dnz,PetscInt* dnnz,PetscInt onz,PetscInt* onnz,Mat* B,PetscInt dnzB,PetscInt* dnnzB,PetscInt onzB,PetscInt* onnzB,Vec* vecb,Vec *vece,PetscInt rank,PetscInt rank_hsl,PetscInt mpisize,char* tabfile, char *commsyntax,set_def *sets,dim_t nset, set_element *set_elems, array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar, elem_value **elem_vals2,offset_t ncofvar,offset_t ncofele,offset_t nvarele,closure_entry **closure_vals2,offset_t alltimeset,offset_t allregset,offset_t nintraeq,dim_t matsol,PetscInt Istart,PetscInt Iend,  offset_t nreg, offset_t ntime, PetscInt *eq_addr, offset_t ndblock, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,dim_t laDi,dim_t laD,PetscReal cntl3,PetscReal cntl6,dim_t nesteddbbd,int localsize,PetscInt *ndbbddrank1,fortran_int* indata,dim_t mc66,fortran_int *ptx,struct timeval begintime,dim_t subints,MPI_Fint fcomm,int solmethod,solve_real **xcf2);
/* solve_rk.c — Runge-Kutta drivers (GEMPACK 26.5; Schiffmann 2022):
   solmethod picks the flavor (SM_RK2/SM_RK4/SM_BOSHA32/SM_DOPRI54),
   steps1 the (initial) step count. adaptive: 0=no, 1=yes,
   2=accuracy-only; epstol/retryadj/maxretries tune the embedded
   controller. accmetric2 receives the per-element cumulative error
   metrics (embedded flavors only; main writes them to the .acc file). */
bool solve_rk(PetscBool nohsl,PetscInt VecSize,PetscInt dnz,PetscInt* dnnz,PetscInt onz,PetscInt* onnz,PetscInt dnzB,PetscInt* dnnzB,PetscInt onzB,PetscInt* onnzB,Vec *vece1,PetscInt rank,PetscInt rank_hsl,PetscInt mpisize,char* tabfile, char *commsyntax,set_def *sets,dim_t nset, set_element *set_elems, array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar, elem_value **elem_vals2,offset_t ncofele,offset_t nvarele,closure_entry **closure_vals2,offset_t alltimeset,offset_t allregset,offset_t nintraeq,dim_t matsol,PetscInt Istart,PetscInt Iend,offset_t nreg, offset_t ntime, PetscInt *eq_addr, offset_t ndblock, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,dim_t laDi,dim_t laD,PetscReal cntl3,PetscReal cntl6,dim_t nesteddbbd,int localsize,PetscInt *ndbbddrank1,fortran_int* indata,dim_t mc66,fortran_int *ptx,struct timeval begintime,MPI_Fint fcomm,int solmethod,int adaptive,double epstol,double retryadj,int maxretries,rk_options *rko,solve_real **xcf2,solve_real **accmetric2);
/* solve_rk.c -- C2 complementarity approximate run (design doc
   section 8; manual 51.1.2/51.7.3): single-solution forward Euler
   with per-step state evaluation, the del_comp@ Newton correction
   shocked 1 in full each step and step redo on state flips.
   napprox = requested Euler steps (51.6), redo_steps 0/1,
   redo_min_frac the minimum redone-step fraction (default 0.005). */
bool solve_comp_approx(PetscBool nohsl,PetscInt VecSize,PetscInt dnz,PetscInt* dnnz,PetscInt onz,PetscInt* onnz,PetscInt dnzB,PetscInt* dnnzB,PetscInt onzB,PetscInt* onnzB,Vec *vece1,PetscInt rank,PetscInt rank_hsl,PetscInt mpisize,char* tabfile, char *commsyntax,set_def *sets,dim_t nset, set_element *set_elems, array_def *coefs,offset_t ncof,array_def *vars,offset_t nvar, elem_value **elem_vals2,offset_t ncofele,offset_t nvarele,closure_entry **closure_vals2,offset_t alltimeset,offset_t allregset,offset_t nintraeq,dim_t matsol,PetscInt Istart,PetscInt Iend,offset_t nreg, offset_t ntime, PetscInt *eq_addr, offset_t ndblock, offset_t *countvarintra1, offset_t *counteq, offset_t *counteqnoadd,dim_t laA,dim_t laDi,dim_t laD,PetscReal cntl3,PetscReal cntl6,dim_t nesteddbbd,int localsize,PetscInt *ndbbddrank1,fortran_int* indata,dim_t mc66,fortran_int *ptx,struct timeval begintime,MPI_Fint fcomm,int napprox,int redo_steps,double redo_min_frac,solve_real **xcf2);
#endif // TEEMS_SOLVER_H_INCLUDED

