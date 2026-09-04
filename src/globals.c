/* Single definitions of the program-wide globals declared extern in
 * teems_solver.h. */
#include <teems_solver.h>
#include <sys/resource.h>

int verbosity = 1;
int inmemory;
int section_threads;
int max_threads;
double step_ratio2, step_ratio3, extrap_w1, extrap_w2, extrap_w3;
int steps1, steps2, steps3;
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

/* Read (IfHeaderExists) satisfied per coefficient (manual 11.11.8;
   plan 3.9): under the engine's reads-then-formulas pass order, a
   satisfied conditional read supersedes formulas assigning the same
   coefficient -- reproducing GEMPACK's file-order outcome for the
   default-formula + conditional-read idiom */
bool *teems_coef_ifhdr = NULL;

/* set mappings (manual 11.9): filled in main once declarations and
   by_elements values are read and broadcast */
map_def *teems_maps = NULL;
dim_t teems_nmap = 0;

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

/* dual-class zerodivide (plan A1): scanner-tracked state + the
   GEMPACK-semantics switch (-gpzerodivide, default legacy) */
zdiv_state teems_zdiv_scan = { 0, 0, 1, 0 };
int teems_gpzerodivide = 0;

/* la* auto-sizing record: max grown -la* equivalent percent observed
   this run (0 = the configured size never grew); reduced across ranks
   and patched into stats.json after the solve */
long teems_laA_used = 0;
long teems_laDi_used = 0;
long teems_laD_used = 0;

/* closure side arrays (6.16(b)); see closure_entry in teems_solver.h */
unsigned char *teems_cl_flags = NULL;
store_real *teems_cl_shock = NULL;

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
