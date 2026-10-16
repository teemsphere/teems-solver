# teems-solver 1.1.0

## Breaking changes
* 1.1.0 is not compatible with 1.0.0 or with the teems R package 0.1.1, despite the minor-version step: the flags, the data layout and the binary's path all changed, and 1.1.0 is the first version that states its version (`-version`) and is checked by the R package before every solve (teems 1.0.0 requires solver ≥ 1.1.0). Users of 0.1.1 keep their current image until they upgrade both; see the README's "Updating an Existing Installation".
* The binary is renamed `hsl` → `teems-solver`
* Removed options: `-presol`, `-nesteddbbd` (derived from `-matsol`), `-regset`/`-enable_time` (the block structure is detected from the equations), `-ndbbd_bl_rank`, `-nestfile`, `-medthreads`, `-isLinux`, `-stoiter`
* `-solmed Mmid` is renamed `Gragg` (`Mmid` is a deprecated alias); the Stochastic/StoSim methods are removed (they never sampled); an unknown `-solmed` aborts
* Assertions, range tests, PostSim and complementarity controls are command-line flags, not CMF statements
* Shock and closure files are read strictly by the GEMPACK command-file rules (manual 20.7): unrecognized lines abort with their line number instead of being skipped
* Shock values are read in GEMPACK order (first subscript fastest) for every statement shape
* `-version` reports the solver version; the exit status follows a documented contract; teems-R 1.0.0 requires solver ≥ 1.1.0

## Bug fixes affecting results
* A uniform shock on a one-index variable restricted to an element or subset (`Shock qgdp("usa") = uniform x;`) shocked element 0 / the first n elements of the full set
* GTAPv7: the scalar `pxwwld` was dropped from `E_c1_cr`, `E_cnttotr` and `E_cntpinv` because a variable followed by a closing bracket was not recognized
* A shock with no argument list (`Shock pfactwld = uniform 1;`) shocked only the variable's first element
* DBBD solves left residuals of 1e-10 to 1e-5 that multi-step runs accumulated (on a 1.4M-equation GTAPv7 run, ~42,000 elements more than 1e-2 off LU); every DBBD solve now takes one refinement step with its kept factors
* DBBD/NDBBD could solve intertemporal models with sum-over-time or fixed-time references wrongly (border classification hole)
* Conditional formulas could read the wrong element when run on more than one thread (OpenMP race)
* Recursive formulas gave thread-count-dependent results; an offset on a Formula/Update left-hand side wrote element 0
* A later Update read earlier Updates' new values instead of start-of-step values (the GTAP models are unaffected)
* A repeated index (`T(d,d)`) read off-diagonal elements; a variable repeated inside a sum was counted |E| times
* `Set B = A;` on a set read from file skipped building the elements; `Formula & Equation` silently dropped the equation; set products produced empty sets
* A misspelled element or set in a closure marked the wrong elements exogenous; a scalar shock on an undeclared variable, tab-indented or BOM-prefixed shock lines, and malformed shock lines were dropped, and the model solved unshocked
* SBBD with no chain dimension exited 0 without a solution
* SBBD/DBBD wrote past the CGLOB array on every solve (results were correct by chance)
* Very large systems: MA48 workspace sizes wrapped at 32 bits and MP48 per-rank totals overflowed at 2^31

## Solution methods
* Euler (single run or Richardson extrapolation); Gragg with a single run (`-single_run`) or all-odd step counts
* Runge-Kutta: RK2, Heun, RK4, BoSha32, DoPri54, with adaptive step control, a log-level chart, and error estimates in `sol.est`
* `-solmed probe` (structural diagnosis) and `-solmed nosim` (formulas only); `-jacdump` exports the base-point Jacobian
* Subtotals using the kept factorization of every matrix method
* Complementarities (approximate/accurate runs) under all matrix methods
* Condensation: `Backsolve`/`Substitute` statements

## TABLO coverage
* Set expressions and equality, conditional set builders, set products, `$POS`, `Mapping` (in formulas, equations, sums, updates), levels variables and equations, `Complementarity`, `PostSim`, `Assertion`, `Read (IfHeaderExists)`, `Default`, ZERODIVIDE semantics, second range bound
* Intrinsics: EXP, SQRT, LOGE, LOG10, ROUND, ID0V, MAX, MIN, TRUNC0, TRUNCB, NORMAL, CUMNORMAL, LOGNORMAL, CUMLOGNORMAL, GPERF, GPERFC, RANDOM (`-random_seed`)
* Conditional sums and quantifier conditions comparing any two expressions; word comparison operators
* Case-insensitive element, header and name matching
* Malformed TAB, data, shock and closure input stops with a named error instead of a crash (sanitizer and fuzzing campaigns)

## Performance and resources
* Compiled equation programs are cached across Jacobian refills (wall −52% on small models)
* `-inmemory` (on by default except NDBBD) keeps value arrays resident
* `-fastrefac`: refactorize with the step-1 analysis for LU, SBBD and DBBD
* MA48/MP48 workspaces grow and retry instead of aborting; `-ma48u` sets the pivot threshold; `-condest` reports solve quality
* `-f64` image binary for double-precision coefficient storage
* NDBBD thread budget per region; OpenBLAS kernel pinned per architecture; images built per x86-64 ISA level
* Named aborts at the PetscInt and 32-bit MP48 ceilings

## Build
* HSL is built from the single libHSL tarball; `-Wall`-clean build
