/* Single definitions of the program-wide globals declared extern in
 * teems_solver.h. */
#include <teems_solver.h>

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

/* CMF "Assertions = yes|no|warn" (manual 25.3): 0 = no, 1 = warn,
   2 = yes (default) -- consumed by assertions_execute */
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

/* CMF range-test modes (manual 25.4.4): 0 = no, 1 = warn (the
   GEMPACK default outside automatic accuracy), 2 = fatal */
int teems_range_test_initial = 1;
int teems_range_test_updated = 1;

/* dual-class zerodivide (plan A1): scanner-tracked state + the
   GEMPACK-semantics switch (-gpzerodivide, default legacy) */
zdiv_state teems_zdiv_scan = { 0, 0, 1, 0 };
int teems_gpzerodivide = 0;

void zdiv_scan_reset(void) {
  teems_zdiv_scan.zbz_val=0;
  teems_zdiv_scan.nbz_val=0;
  teems_zdiv_scan.zbz_on=1;
  teems_zdiv_scan.nbz_on=0;
}
