/* Single definitions of the program-wide globals declared extern in
 * ha_cgeglobal.h. */
#include <ha_cgeglobal.h>

int isLinux;
int smallthreads;
int medthreads;
int mymaxnumthrd;
double kindx1, kindx2, kval1, kval2, kval3;
int step1, step2, step3;
MPI_Comm HA_COMM, HA1_COMM;
char temdir[NAMESIZE] = "/tmp/";
