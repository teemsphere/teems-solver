# TEEMS Solver — Development Roadmap

Working plan following the 2026 refactor (phases 0–4, logged below)
and phase 5's measured performance work. Completed items are never
deleted: phase-5 rows flip to **done** with findings attached, and
phases 0–4 are recorded in the completed-phases log at the end of this
file. `docs/solver-reference.md` describes the *resulting* system.

**Standing rules**
- Every change passes the golden suite (`.audit/verify.sh`, currently 9
  bit-identical gates: 5 method + 3 real-shock + 1 GEMPACK-orientation) and,
  where performance-relevant, interleaved A/B benchmarks on the
  calibrated rigs. Benchmark only on an idle machine; use paired runs.
- **Goal: full GEMPACK feature parity. Performance is the only veto.**
  Features must be zero-cost when unused. The GEMPACK manual is the
  parity spec; original TABs in `original_tabs/` (GTAP resource
  archives) and `teems_archive/tab_files/` (incl. `GDYNv3.6.tab`,
  `gtapep.tab`) are the test corpus — the teems-adapted TABs in
  `teems-R/inst/models/` are already reduced and cannot reveal gaps.

## Phase 5 — remaining performance/correctness work

| id | item | status / notes |
|---|---|---|
| 5.5 | Golden-suite widening (real asymmetric shocks: 2D probe, 4D, swap) | **done** — 8 goldens |
| 5.6 | Shock-value ordering aligned to GEMPACK standard (all statement forms); teems-R writer paired (cee6210); golden 9 `gmshk` pins orientation; teems-R custom-shock suite 18/18 | **done** |
| 5.7 | Netcut instrumentation: machine-readable netcut + block-size logging per run | **done** — `<solfiles>.stats.json` written by rank 0 at ordering time (pre-solve): size/method/netcut/border/per-block var+eq counts; excluded from golden manifests |
| 5.8 | `matrix_method = "auto"`: calibration sweeps (size × T/R × ranks), decision rule implemented in teems-R from deploy metadata | **done** — teems-R 3d561a2: deploy writes system_size/n_reg into metadata.rds; auto (now the default) = SBBD if intertemporal (single-rank SBBD also beats LU), else DBBD when n_tasks ≥ 2 and size ≥ 2M eq (or ≥ 1.5M with ≥ 100 regions), else LU. Static crossover data in solver-reference §9 (LU/DBBD ladders at 1.35M–8.3M eq × 33/163 blocks). Residual: NDBBD-vs-SBBD escalation (standing investigation); refine 1.5–2M × <100-region band if real models land there |
| 5.9 | In-memory block-factor handoff: redesign Fortran kernel interface so factors return via memory, not `_vav/_irnv/_keep` files | **done** (code + correctness) — kernels return factors via caller buffers (`KEEP` surfaced with the `M+9N+7` bound; `LA` was always known pre-call, so no grow-loop issue); C holds them resident under `-inmemory` or writes byte-identical legacy files otherwise. All 9 goldens pass in both modes for DBBD and NDBBD. NDBBD resident mode: zero factor traffic (was ~825MB/solve at 202k eq), +22MB RSS/rank, no scratch debris. **Open**: idle-machine A/B at bench-inter-L scale to decide flipping NDBBD's `-inmemory` default (its page-cache rationale is gone); residual `_bivi/_rbvi/_cbvi`+`_rank/_row/_col` interface spills (~65MB/solve, C-side both ends) as follow-up |
| 5.10 | Interpreter deep-compaction: interleaved per-dim operand records (`dim_addr`) in `formula_op`; follows the 1416→968B compaction (~1% win). `formula_eval` is 15.6% of wall at scale | **done** (code + correctness) — `dim_addr {ADims, SupSet, SSIndx, leadlag}` replaces the 4 parallel per-dim arrays per operand; the eval loop now reads one 24B record per dim instead of touching 4 arrays 96B apart. `formula_op` 968→848B; `SupSet`/`SSIndx` narrowed offset_t→int. All 9 goldens bit-identical; build warnings 110→102. **Open**: idle-machine A/B at scale to measure the `formula_eval` win. NB pre-existing latent bug preserved as-is (goldens are bit-identical gates): the Var2→Var3 dim-field copy self-assigns `SupSet` (`Var3Dims[i].SupSet=Var3Dims[i].SupSet`); ops arrays are calloc'd and slots written once, so operand 3's SupSet is always 0 there — superset indexing silently dropped if a superset-indexed operand lands in that slot (unexercised by the golden corpus). Fix as its own gated change |
| 5.11 | Small batch: verbosity-gated logging **including a full printf overhaul — professional wording (no `!!!!`/`OK012345`/`Hello world`/`RRRRR` debris), each message including errors as informative as possible (say which file/block/rank/phase and what to do about it) without impeding speed (hot-loop prints gated or removed; keep stdout parse contracts teems-R relies on, e.g. `Accurate`/netcut lines, or update both sides together)**; `-march` image variants; `ems_solve(inmemory=)` exposure; op-list caching (small-model benefit only); remove dead `isLinux` read path; NDBBD-without-`-regset` diagnostic (currently bare MPI abort) | pending (user emphasis 2026-07-06 on message quality) |

**Phase-5 findings that shape the above** (details in
`docs/solver-reference.md` §9): factorization ≈63% of wall at scale;
`formula_eval` ≈16%; per-step TAB re-parse only matters on small models.
`-inmemory` defaults: on for LU/SBBD/DBBD (+7% SBBD, neutral elsewhere,
no debris), off for NDBBD (its factor-file traffic needs the page cache
that spilling frees). Speed since refactor: SBBD ~5–8%; LU/DBBD neutral
+ cleanliness; NDBBD unchanged by design.

## Phase 6 — GEMPACK parity program (performance-vetoed)

| id | item | notes |
|---|---|---|
| 6.1 | **Parity inventory** — parse the original TAB corpus against the current front end; ranked checklist of unsupported statements/qualifiers with real-model frequencies. Refresh upstream docs (current GEMPACK, HSL, PETSc v3.25, MPICH) into `teems-docs/` | gates the rest of phase 6 |
| 6.2 | **Condensation** (substitute/backsolve/omit variables pre-factorization) — the largest speed lever in the project: shrinks the system that costs 63% of wall, and reduces netcut as a side effect. Design choice: teems-R-side TAB rewriting vs solver-side elimination. Must preserve the `.bin`/`.var` output contract (condensed variables recovered post-solve); interacts with closure/swap validation | highest-leverage item overall |
| 6.3 | **Stepping methods** (in cost-benefit order): (a) Euler multi-step — near-free (Mmid driver with midpoint flag = 0), immediate parity gap; (b) verify/add Gragg terminal smoothing pass (manual §30.2 distinguishes gragg vs midpoint; TEEMS does Euler-start + midpoint steps + 2-4-8 Richardson); (c) Runge-Kutta: RK2/RK4/BoSha32/DoPri54 via Butcher-tableau driver — stages reuse the existing refill+solve primitive, compatible with all matrix methods; adaptive accept/reject must be MPI-collective (Allreduce the embedded-error norm) and the controller tuned conservatively at scale (a rejected step wastes a full set of block factorizations); (d) Adams-Bashforth — ONE factorization per step at order k (vs RK4's 4): directly attacks the factorization wall; RK startup; stability caveats; (e) implicit stiff methods (Radau/BDF/ode15s) — parked: performance-negative for smooth CGE paths; revisit only if adaptive RK step-size collapse demonstrates stiffness | RK gate: fixed-accuracy wall-time vs Gragg 2-4-8 |
| 6.4 | Statement completion per 6.1 ranking: full IF semantics (opcode skeleton exists: `OP_IF_*`), POSTSIM, others | each zero-cost when unused |
| 6.5 | **Netcut enforcement** (border size = interface-problem size, governs bordered-method viability): E1 validate in `ems_model()` (error/warn on TIME×large-set variables in inter-period links); E2 auto-rewrite minimal intertemporal proxies during TAB assembly; E3 solver-side element-level border classification — E3 after 6.2, since condensation changes what reaches the border | E1 cheap, E2 sweet spot |
| 6.6 | Cross-repo hygiene: `hsl` → `teems-solver` binary rename (touch points documented in project notes); golden refresh | coordinate with teems-R + images |
| 6.7 | **HSL catalogue exploration** — assess additional HSL routines for speed/robustness across all solution × matrix methods; ranked shortlist with rationale in `teems-docs/hsl-catalogue-candidates.md` (top: MC64 pre-scaling, MA41/MA38 multifrontal alternative, MC33 bordered-block-triangular ordering vs netcut, HSL_MC79 Dulmage-Mendelsohn closure diagnosis). Retrieve specs/code for survivors; golden + A/B benchmark gates apply | user-requested 2026-07-05 |
| 6.8 | **Block-set auto-detection** (assessed 2026-07-06): drop the need for explicit `-regset`/`-enable_time`. Time side is already derivable — `(intertemporal)` TAB qualifier + `set_find_alltime`; the missing-flag case is even detected post-hoc (`sbbd_overrid` warning); NDBBD's timestep count self-derives from the alltime set (`ndbbdrank=ntime`). Regional side has no TAB qualifier and is name-convention only (teems-R hardcodes `REG`; unmatched names silently degrade to `allregset=-1`). All ordering code downstream is set-agnostic ("partition by set S"), so implement `-regset auto`: enumerate top-level non-intertemporal candidate sets, run the cheap pre-Jacobian ordering probe per candidate, score by netcut/border fraction, block count ≥ n_tasks, and block balance (the 5.7 stats fields), pick argmax or none (LU). Explicit flags stay as overrides (golden runs keep bit-identical behavior); log the choice in stats.json. Also fixes the silent-degradation footgun (warn when a named set is not found) | feeds `matrix_method` auto; pairs with 6.5/E3 |

## Completed phases 0–4 (2026 refactor log)

Reconstructed from git history; hashes are the anchor commits.

| phase | what was done | commits |
|---|---|---|
| 0 — verification harness | `.audit/verify.sh`: container rebuild (`teems-audit`) + golden solves checked bit-identical against manifests anchored to the pre-refactor binary; grew from 5 gates to 9 during phase 5. Benchmark rigs (`bench_run.sh` wall/RSS, `strace_run.sh` write-byte accounting, `bench-*` deployments). **Dev-machine only — not tracked in git**; described in solver-reference §12 | — |
| 1 — dead-code removal | 13 uncalled C functions + prototypes; unreachable solver-dispatch branches in `main()`; commented-out code stripped from C and Fortran; dead backup subroutine in `ha_mp48.f90`; dangling prototype; 15 uncalled Fortran kernels + externs; dead `-isLinux`/`-medthreads` options | `fdd2e76`, `08c1be8`, `4a05297`, `bab50d5`, `fa0b51b`, `b5669fe`, `7a403a5` |
| 2 — build hygiene | makefile repair (clean target, F77FLAGS typo, overridable `BUILD_DIR`/`OPT`); globals defined once, `-fcommon` dropped; PETSc private `aij.h` via include path; scratch dir configurable (`-tempdir`/`TMPDIR`); `mp48_mod.sh` sed edits replaced with committed patches; expedited build compiles current src (was stale snapshot); `BUILDING.md` added | `866fb85`, `6a958bc`, `1bfde3a`, `17f9b3a`, `74feb38`, `f775a4f`, `7a403a5` |
| 3 — de-duplication | inline Johansen and modified-midpoint blocks factored into `Johansen()`/`ModMidPoint()`; `hnew_mupdate` merged into `hnew_update` behind a midpoint flag | `15b9d5a`, `ffd7eb6`, `7a13606` |
| 4 — restructuring & renaming | `NAMING.md` conventions; types/fields/globals, functions, and source files renamed; `ha_newmfparse.c` split into `jacobian.c`/`block_order.c`/`block_solve.c`/`solve_drivers.c`; magic codes replaced with enums (`matrix_method`, `solution_method`, `bound_type`, `op_code`, `operand_type`); header documented by module section; residual `ha_*` parameters/locals/generated-name prefixes renamed | `f3ccb8e`, `fb441fd`, `92ec64c`, `b081034`, `e5a5158`, `ea9624f`, `e0f0e4c`, `4aff8be` |

Speed outcome across phases 0–4: neutral by design (golden-gated);
cleanliness and navigability were the deliverables. Performance work
began in phase 5.

## Standing investigations

- `-ndbbd_bl_rank`/`-nestfile` (expert overrides: hand-assign MPI ranks
  per NDBBD block from a CSV; default self-derives `ndbbdrank=ntime`):
  look into how they work and where they earn their keep — plausibly
  load-balancing heterogeneous block sizes (unbalanced regional blocks,
  the per-time remainder blocks visible in stats.json `block_nvar`) or
  pinning blocks on heterogeneous/cluster nodes. If a real use case
  exists, document it and consider folding into the auto-calibration
  (compute the assignment from measured block sizes instead of a
  hand-written CSV); if not, candidate for removal. NB the fallback
  `strcpy(filename,"./ndbbd_drank.csv")` in main.c is dead code — the
  CSV is only read when `-nestfile` is explicitly passed.

- NDBBD-vs-SBBD crossover at scale (medium ACTS, 8–16 ranks, cluster
  hardware) — dev-box sweep found no NDBBD win region.
- Stiffness monitoring once adaptive RK exists (see 6.3e).
- Ordering freedom in `equation_order_read_nested` (legacy note: "can
  try different varreg"); direct vs MC66 ordering netcut trade-off.
