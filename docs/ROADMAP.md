# TEEMS Solver — Development Roadmap

Working plan following the 2026 refactor (phases 0–4: verification
harness, dead-code removal, build hygiene, de-duplication, restructuring
and renaming — see `docs/solver-reference.md` §11–§13 and git history)
and phase 5's measured performance work.

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
| 5.9 | In-memory block-factor handoff: redesign Fortran kernel interface so factors return via memory, not `_vav/_irnv/_keep` files. Removes the page-cache competition that forces NDBBD's spill default. Constraint: factor sizes unknown pre-factorization (`laA` growth loop) → allocate-on-demand or grow-and-retry | pending (largest phase-5 effort) |
| 5.10 | Interpreter deep-compaction: interleaved per-dim operand records (`dim_addr`) in `formula_op`; follows the 1416→968B compaction (~1% win). `formula_eval` is 15.6% of wall at scale | pending (medium risk — hottest loop) |
| 5.11 | Small batch: verbosity-gated logging; `-march` image variants; `ems_solve(inmemory=)` exposure; op-list caching (small-model benefit only); remove dead `isLinux` read path; NDBBD-without-`-regset` diagnostic (currently bare MPI abort) | pending |

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

## Standing investigations

- NDBBD-vs-SBBD crossover at scale (medium ACTS, 8–16 ranks, cluster
  hardware) — dev-box sweep found no NDBBD win region.
- Stiffness monitoring once adaptive RK exists (see 6.3e).
- Ordering freedom in `equation_order_read_nested` (legacy note: "can
  try different varreg"); direct vs MC66 ordering netcut trade-off.
