/* Single definitions of the program-wide globals declared extern in
 * teems_solver.h. */
#include <teems_solver.h>
#include <sys/resource.h>
#include <stdarg.h>

int verbosity = 1;
int teems_error_count = 0;
int errmsg(const char *fmt, ...) {
  va_list ap;
  int n;
  va_start(ap, fmt);
  n = vprintf(fmt, ap);
  va_end(ap);
  fflush(stdout);
  teems_error_count++;
  return n;
}
int inmemory;
int section_threads;
int max_threads;
double step_ratio2, step_ratio3, extrap_w1, extrap_w2, extrap_w3;
int steps1, steps2, steps3;
/* -single_run 1: one multi-step Euler/Gragg pass over -step1 steps, no
   Richardson extrapolation (GEMPACK "method = euler; steps = N;") */
int teems_single_run=0;
/* -two_run 1: extrapolate from two multi-step solutions (manual 26.1.2) */
int teems_two_run=0;
/* -sup (manual 26.8.1 SUP): 0 none, 1 last, 2 all multi-step solutions
   write their updated data (<stem>.ud5/.ud6/.ud7) */
int teems_sup=0;
/* -comp_do_acc 0 (51.5.6): the approximate run is the result */
int teems_comp_no_acc=0;
/* complementarity subintervals (51.7.4): their count, and the result
   compounded over the finished ones, onto which .xac pass values add */
int teems_comp_nsub=1;
/* loop control probe (tab_loop.c): assertions_execute counts the tuples
   of the "(loopctl)" statements in the window and the true ones */
int teems_loopctl_probe=0;
long teems_loopctl_n=0,teems_loopctl_true=0;
solve_real *teems_comp_xac_base=NULL;
/* -sui 1 (manual 26.8.2 SUI): .cof kind bit 2 marks the coefficients a
   Formula (Initial) sets, whose updated values are in the .cbin */
int teems_sui=0;
bool *teems_coef_is_fini = NULL;
/* RANDOM (manual 11.5.2): a value is a hash of the seed, the statement
   text, the occurrence in the statement and the element tuple, so every
   re-evaluation (steps, passes, ranks) draws the same number */
long teems_random_seed=1;
uint64_t teems_rand_stmt=0;
int teems_rand_count=0;
MPI_Comm node_comm, node_tail_comm;
char scratch_dir[NAMESIZE] = "/tmp/";
backsolve_def *backsolves = NULL;
int nbacksolve = 0;
offset_t nbselems = 0;
int backsolve_scan_mode = BS_SCAN_SKIP;

/* -assertions run switch (manual 25.3): 0 = off, 1 = warn,
   2 = fatal (default) -- consumed by assertions_execute */
int teems_assertions_mode = 2;

/* (parameter)-qualified coefficients (PostSim foundation F2): parallel
   to coefs[] -- array_def itself is binary-locked to sol.var */
bool *teems_coef_is_param = NULL;
/* INTEGER coefficients (manual 11.6.3): their Formulas default to INITIAL */
bool *teems_coef_is_int = NULL;
offset_t teems_n_int_coefs = 0;

/* second declared-range bound (audit A9): a declaration may carry one
   lower (GE/GT) and one upper (LE/LT) bound (manual 10.19.1); slot 1
   lives in array_def.gltype/glval (binary-locked to sol.var), the
   other-direction slot rides here, parallel to coefs[] */
int *teems_coef_gltype2 = NULL;
store_real *teems_coef_glval2 = NULL;

/* PostSim scope (Tier 0 residuals; manual 12.2.1-12.2.3): coefficient
   names declared inside POSTSIM sections, recorded by the split;
   is_ps parallels coefs[] (NULL when the TAB has no sections);
   ps_pass gates the PostSim-only formula rules */
char (*teems_ps_coefnames)[NAMESIZE] = NULL;
int teems_ps_ncoefs = 0;
bool *teems_coef_is_ps = NULL;
int teems_ps_pass = 0;
char (*teems_ps_wlogs)[NAMESIZE] = NULL;
int teems_ps_nwlogs = 0;
int teems_ps_ran = 0;

/* statement order (manual 10.1, 11.11.8, 12.2.1) */
long teems_stmt_start = -1;
long teems_ord_lo = -1, teems_ord_hi = -1;

/* set mappings (manual 11.9): filled in main once declarations and
   by_elements values are read and broadcast */
map_def *teems_maps = NULL;
dim_t teems_nmap = 0;
set_def *teems_sets = NULL;
dim_t teems_nset = 0;
set_element *teems_set_elems = NULL;
bool *teems_set_isprod = NULL;
dim_t *teems_set_prod1 = NULL;
dim_t *teems_set_prod2 = NULL;
/* intertemporal element stem of a p[a] - p[b] set ("" otherwise) and
   the number a of its first element (manual 16.2.1) */
char (*teems_set_itstem)[NAMESIZE] = NULL;
int *teems_set_itfirst = NULL;

/* complementarities (manual 10.17/11.14; design doc sections 7-8):
   records filled by tab_complementarity_transform on rank 0,
   broadcast in main alongside the mappings; consumed by
   complementarities_validate, comp_closure_check and the C2 state
   machinery (comp_states_*) */
comp_def *teems_comps = NULL;
dim_t teems_ncomp = 0;
/* endogenous complementarity-variable components (C2 active mode);
   set by comp_closure_check, broadcast after the closure section */
offset_t teems_comp_active = 0;

/* -range_test_initial/-range_test_updated run switches (manual
   25.4.4): 0 = off, 1 = warn (the GEMPACK default outside automatic
   accuracy), 2 = fatal */
int teems_range_test_initial = 1;
int teems_range_test_updated = 1;
teems_rk_stats_t teems_rk_stats;
int teems_rk_stage_checks = 0;
long teems_check_viol_range = 0;
long teems_check_viol_assert = 0;
int teems_rk_softfail = 0;
int teems_stage_solve_failed = 0;

/* dual-class zerodivide (manual 10.11): scanner-tracked state */
zdiv_state teems_zdiv_scan = { 0, 0, 1, 0 };

/* la* auto-sizing record: max grown -la* equivalent percent observed
   this run (0 = the configured size never grew); reduced across ranks
   and patched into stats.json after the solve */
long teems_laA_used = 0;
long teems_laDi_used = 0;
long teems_laD_used = 0;

/* closure side arrays (6.16(b)); see closure_entry in teems_solver.h */
unsigned char *teems_cl_flags = NULL;
int teems_nsub = 0;
sub_group *teems_subs = NULL;
store_real *teems_cl_shock = NULL;

int teems_ndcutcache = 1; /* -ndcutcache: reuse the NDBBD cuts across steps */

/* wall-clock stage record (bordered drivers); see teems_solver.h */
static const char *teems_stage_name[TEEMS_STAGE_MAX];
static double teems_stage_wall[TEEMS_STAGE_MAX];
static int teems_nstage=0;
static const char *teems_stage_cur=NULL;
static double teems_stage_t0=0.0;
void teems_stage_mark(const char *name) {
  int i;
  if(omp_in_parallel()&&omp_get_thread_num()!=0)return;
  double t=MPI_Wtime();
  if(teems_stage_cur!=NULL) {
    for(i=0; i<teems_nstage; i++)if(strcmp(teems_stage_name[i],teems_stage_cur)==0)break;
    if(i==teems_nstage&&teems_nstage<TEEMS_STAGE_MAX) {
      teems_stage_name[teems_nstage]=teems_stage_cur;
      teems_stage_wall[teems_nstage]=0.0;
      teems_nstage++;
    }
    if(i<TEEMS_STAGE_MAX)teems_stage_wall[i]+=t-teems_stage_t0;
  }
  teems_stage_cur=name;
  teems_stage_t0=t;
}
void teems_stage_report(const char *prefix) {
  double mx[TEEMS_STAGE_MAX];
  int i,rank=0,n=teems_nstage;
  teems_stage_mark(NULL);
  MPI_Allreduce(MPI_IN_PLACE,&n,1,MPI_INT,MPI_MIN,PETSC_COMM_WORLD);
  MPI_Allreduce(teems_stage_wall,mx,n>0?n:1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  if(rank==0&&n>0) {
    double tot=0.0;
    for(i=0; i<n; i++)tot+=mx[i];
    printf("%s stages (s, max over ranks; total %.2f):",prefix,tot);
    for(i=0; i<n; i++)printf("%s %s %.2f",i?",":"",teems_stage_name[i],mx[i]);
    printf("\n");
    fflush(stdout);
  }
  teems_nstage=0;
  teems_stage_cur=NULL;
}

/* phase resident-memory record (6.16(a)); see teems_solver.h */
teems_rss_entry teems_rss[TEEMS_RSS_MAX];
int teems_nrss = 0;
void teems_rss_probe(const char *phase) {
  double loc[2]={0.0,0.0},mx[2]={0.0,0.0},sm[2]={0.0,0.0};
  long pages=0,resident=0;
  int rank=0,size=1,i;
  FILE *fp=fopen("/proc/self/statm","r");
  if(fp!=NULL) {
    if(fscanf(fp,"%ld %ld",&pages,&resident)!=2)resident=0;
    fclose(fp);
  }
  loc[0]=(double)resident*(double)sysconf(_SC_PAGESIZE)/1073741824.0;
  struct rusage ru;
  if(getrusage(RUSAGE_SELF,&ru)==0)loc[1]=(double)ru.ru_maxrss*1024.0/1073741824.0;
  MPI_Reduce(loc,mx,2,MPI_DOUBLE,MPI_MAX,0,PETSC_COMM_WORLD);
  MPI_Reduce(loc,sm,2,MPI_DOUBLE,MPI_SUM,0,PETSC_COMM_WORLD);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  if(rank!=0)return;
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  logmsg(1,"memory: after %s, resident %.2f GB max per rank, %.2f GB over %d rank(s); high-water %.2f GB max, %.2f GB sum\n",phase,mx[0],sm[0],size,mx[1],sm[1]);
  for(i=0; i<teems_nrss; i++)if(strcmp(teems_rss[i].phase,phase)==0)break;
  if(i==teems_nrss) {
    if(teems_nrss==TEEMS_RSS_MAX)return;
    teems_nrss++;
    teems_rss[i].phase=phase;
    teems_rss[i].rss_max=teems_rss[i].rss_sum=teems_rss[i].hwm_max=teems_rss[i].hwm_sum=0.0;
    teems_rss[i].count=0;
  }
  teems_rss[i].count++;
  if(mx[0]>teems_rss[i].rss_max)teems_rss[i].rss_max=mx[0];
  if(sm[0]>teems_rss[i].rss_sum)teems_rss[i].rss_sum=sm[0];
  if(mx[1]>teems_rss[i].hwm_max)teems_rss[i].hwm_max=mx[1];
  if(sm[1]>teems_rss[i].hwm_sum)teems_rss[i].hwm_sum=sm[1];
}

/* -condest (MA60/MC71 solve-quality diagnostics, sequential LU path):
   the run flag, a scope gate the LU wrappers set around their kernel
   calls (SPEC48_SSOL2LA also serves the DBBD interface system, which
   must stay inert), and per-run accumulators reduced into stats.json.
   Diagnostic-only: the refined solution is never written back. */
int teems_condest = 0;
int teems_condest_scope = 0;
double teems_condest_kw1max = 0.0;
double teems_condest_kw2max = 0.0;
double teems_condest_omegamax = 0.0;
long teems_condest_solves = 0;
long teems_condest_skips = 0;

/* Fortran-visible: is condest active for the current kernel call? */
int teems_condest_active_(void) {
  return teems_condest_scope;
}

/* Fortran-visible reporter, called once per measured (or skipped)
   solve.  status: 0 = ok, 1 = zero right-hand side (nothing to
   measure), 2 = MA60 error return, 3 = refinement hit its iteration
   cap (estimates are the best obtained, still recorded). */
void teems_condest_report_(int *status,double *omega1,double *omega2,double *erx,double *cond1,double *cond2,int *noiter) {
  if(*status==1) {
    teems_condest_skips++;
    logmsg(1,"condest: zero right-hand side at this solve (null-shock step), nothing to measure\n");
    return;
  }
  if(*status==2) {
    logmsg(1,"condest: MA60 error return, no estimate for this solve\n");
    return;
  }
  teems_condest_solves++;
  {
    double om=(*omega1>*omega2)?*omega1:*omega2;
    if(om>teems_condest_omegamax)teems_condest_omegamax=om;
    if(*cond1>teems_condest_kw1max)teems_condest_kw1max=*cond1;
    if(*cond2>teems_condest_kw2max)teems_condest_kw2max=*cond2;
  }
  logmsg(1,"condest: backward error omega1 %.2e omega2 %.2e (%d refinement passes), forward error bound %.2e, kappa_w1 %.3e, kappa_w2 %.3e%s\n",
         *omega1,*omega2,*noiter,*erx,*cond1,*cond2,
         (*status==3)?" (iteration cap hit; best estimates kept)":"");
  if(*cond2>1e15||*cond1>1e15)
    logmsg(1,"condest: WARNING: the linear system is numerically near-singular at the current values (kappa_w2 %.1e): solutions are unreliable; the structural probe may pass (-solmed probe) -- look for near-zero data flows carried by the closure\n",*cond2);
}

void zdiv_scan_reset(void) {
  teems_zdiv_scan.zbz_val=0;
  teems_zdiv_scan.nbz_val=0;
  teems_zdiv_scan.zbz_on=1;
  teems_zdiv_scan.nbz_on=0;
}

/* Side-car output record (Tier B2). Rank 0 notes every file it writes;
   outputs_json_write lists them in <solfiles>.outputs.json as the very
   last write of a run that printed no Error line. A reader trusts a
   side-car only when this file lists it, so a file left by an earlier
   run in the same directory is never read against this run's solution,
   and the file's presence says the run completed. */
char teems_sol_stem[TABREADLINE] = "";
typedef struct {
  char *path;
  const char *kind;
  int version;
} outputs_entry;
static outputs_entry *outputs_list = NULL;
static long outputs_n = 0, outputs_cap = 0;

void outputs_note(const char *path,const char *kind,int format_version) {
  int rank=0;
  long i;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  if(rank!=0||path==NULL||path[0]=='\0')return;
  for(i=0; i<outputs_n; i++) {
    if(strcmp(outputs_list[i].path,path)==0) {
      outputs_list[i].kind=kind;
      outputs_list[i].version=format_version;
      return;
    }
  }
  if(outputs_n==outputs_cap) {
    long cap=outputs_cap?2*outputs_cap:64;
    outputs_entry *p=(outputs_entry *) realloc (outputs_list,cap*sizeof(outputs_entry));
    if(p==NULL)return;
    outputs_list=p;
    outputs_cap=cap;
  }
  outputs_list[outputs_n].path=strdup(path);
  if(outputs_list[outputs_n].path==NULL)return;
  outputs_list[outputs_n].kind=kind;
  outputs_list[outputs_n].version=format_version;
  outputs_n++;
}

static void json_str(FILE *fp,const char *s) {
  fputc('"',fp);
  for(; *s!='\0'; s++) {
    unsigned char c=(unsigned char)*s;
    if(c=='"'||c=='\\')fprintf(fp,"\\%c",c);
    else if(c<0x20)fprintf(fp,"\\u%04x",c);
    else fputc(c,fp);
  }
  fputc('"',fp);
}

int outputs_json_write(const char *stem,const char *run_id) {
  char path[TABREADLINE+32],tmp[TABREADLINE+40];
  const char *base;
  FILE *fp;
  long i;
  snprintf(path,sizeof(path),"%s.outputs.json",stem);
  snprintf(tmp,sizeof(tmp),"%s.outputs.json.tmp",stem);
  if((fp=fopen(tmp,"w"))==NULL) {
    errmsg("Error: cannot open %s for writing: %s (the solver runs as uid %d)\n",tmp,strerror(errno),(int)getuid());
    return 1;
  }
  fprintf(fp,"{\n  \"version\": 1,\n  \"solver_version\": \"%s\",\n  \"run_id\": ",TEEMS_SOLVER_VERSION);
  json_str(fp,run_id);
  fprintf(fp,",\n  \"complete\": true,\n  \"files\": [");
  for(i=0; i<outputs_n; i++) {
    base=strrchr(outputs_list[i].path,'/');
    base=base?base+1:outputs_list[i].path;
    fprintf(fp,"%s\n    {\"name\": ",i?",":"");
    json_str(fp,base);
    fprintf(fp,", \"path\": ");
    json_str(fp,outputs_list[i].path);
    fprintf(fp,", \"kind\": \"%s\", \"format_version\": %d}",outputs_list[i].kind,outputs_list[i].version);
  }
  fprintf(fp,"\n  ]\n}\n");
  if(fclose(fp)!=0||rename(tmp,path)!=0) {
    errmsg("Error: cannot write %s: %s\n",path,strerror(errno));
    remove(tmp);
    return 1;
  }
  return 0;
}
