/* The solver's one version constant. Surfaced by -version (exact line
   "teems-solver <version>", exit 0), the ungated startup banner (first
   line of every solver log), stats.json / probe.json ("solver_version"),
   and the image labels (the Dockerfiles' TEEMS_VERSION build arg; label
   = pre-pull hint, -version = ground truth). Bump here at tag time. */
#ifndef TEEMS_VERSION_H_INCLUDED
#define TEEMS_VERSION_H_INCLUDED
#define TEEMS_SOLVER_VERSION "1.1.0-dev.4"
#endif
