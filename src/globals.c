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
