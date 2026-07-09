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

## Phase 5 — performance/correctness work (**complete** 2026-07-09)

| id | item | status / notes |
|---|---|---|
| 5.5 | Golden-suite widening (real asymmetric shocks: 2D probe, 4D, swap) | **done** — 8 goldens |
| 5.6 | Shock-value ordering aligned to GEMPACK standard (all statement forms); teems-R writer paired (cee6210); golden 9 `gmshk` pins orientation; teems-R custom-shock suite 18/18 | **done** |
| 5.7 | Netcut instrumentation: machine-readable netcut + block-size logging per run | **done** — `<solfiles>.stats.json` written by rank 0 at ordering time (pre-solve): size/method/netcut/border/per-block var+eq counts; excluded from golden manifests |
| 5.8 | `matrix_method = "auto"`: calibration sweeps (size × T/R × ranks), decision rule implemented in teems-R from deploy metadata | **done** — teems-R 3d561a2: deploy writes system_size/n_reg into metadata.rds; auto (now the default) = SBBD if intertemporal (single-rank SBBD also beats LU), else DBBD when n_tasks ≥ 2 and size ≥ 2M eq (or ≥ 1.5M with ≥ 100 regions), else LU. Static crossover data in solver-reference §9 (LU/DBBD ladders at 1.35M–8.3M eq × 33/163 blocks). Residual: NDBBD-vs-SBBD escalation (standing investigation); refine 1.5–2M × <100-region band if real models land there |
| 5.9 | In-memory block-factor handoff: redesign Fortran kernel interface so factors return via memory, not `_vav/_irnv/_keep` files | **done** (code + correctness) — kernels return factors via caller buffers (`KEEP` surfaced with the `M+9N+7` bound; `LA` was always known pre-call, so no grow-loop issue); C holds them resident under `-inmemory` or writes byte-identical legacy files otherwise. All 9 goldens pass in both modes for DBBD and NDBBD. NDBBD resident mode: zero factor traffic (was ~825MB/solve at 202k eq), +22MB RSS/rank, no scratch debris. A/B closed 2026-07-09 (bench-inter-L NDBBD-2, 3 interleaved pairs): wall neutral (medians 385s im0 vs 381s im1; one 350s im1 run suggests mild upside), RSS 1.09→1.84GB/rank (+69%) for resident factors. **Default stays off** — new rationale is RAM cost, not page cache (that rationale died with 5.9); measured on the environment most favorable to im0 (fast local scratch), so slow/network-scratch clusters — NDBBD's habitat — are expected to favor `inmemory=TRUE`, now exposed via `ems_solve()`. Residual interface spills (~65MB/solve) folded into future work |
| 5.10 | Interpreter deep-compaction: interleaved per-dim operand records (`dim_addr`) in `formula_op`; follows the 1416→968B compaction (~1% win). `formula_eval` is 15.6% of wall at scale | **done** (code + correctness) — `dim_addr {ADims, SupSet, SSIndx, leadlag}` replaces the 4 parallel per-dim arrays per operand; the eval loop now reads one 24B record per dim instead of touching 4 arrays 96B apart. `formula_op` 968→848B; `SupSet`/`SSIndx` narrowed offset_t→int. All 9 goldens bit-identical; build warnings 110→102. A/B closed 2026-07-09 (bench-inter-L SBBD-2, eb7f278 vs 2be256f, 3 interleaved pairs): **speed-neutral** — median delta 0.4% inside ~10% session drift; arithmetic agrees (10% of a 16% `formula_eval` share ≈ 1.6% wall, unresolvable by wall-clock A/B). The change stands on its 120B/op footprint cut and layout clarity. NB pre-existing latent bug preserved as-is (goldens are bit-identical gates): the Var2→Var3 dim-field copy self-assigns `SupSet` (`Var3Dims[i].SupSet=Var3Dims[i].SupSet`); ops arrays are calloc'd and slots written once, so operand 3's SupSet is always 0 there — superset indexing silently dropped if a superset-indexed operand lands in that slot (unexercised by the golden corpus). Fix as its own gated change |
| 5.11 | Small batch: verbosity-gated logging **including a full printf overhaul — professional wording (no `!!!!`/`OK012345`/`Hello world`/`RRRRR` debris), each message including errors as informative as possible (say which file/block/rank/phase and what to do about it) without impeding speed (hot-loop prints gated or removed; keep stdout parse contracts teems-R relies on, e.g. `Accurate` lines, or update both sides together)**; `-march` image variants; `ems_solve(inmemory=)` exposure; op-list caching (small-model benefit only); remove dead `isLinux` read path; NDBBD-without-`-regset` diagnostic (currently bare MPI abort) | **done** — printf overhaul complete: `-verbosity {0,1,2}` (default 1; exported as `TEEMS_VERBOSITY` for the Fortran kernels), all 541 C printf sites triaged (debris deleted, per-rank/per-block detail → level 2, errors/warnings never gated and reworded with file/flag/remedy), MA48 duplicate-entry notes silenced below level 2 (`ICNTL(3)=1`), "no supplied values" warned once per array instead of per step (10,822 → 181 lines on golden-re; v0 log 12,358 → 187 lines), `Accurate` contract preserved at all levels; dead `isLinux` paths removed; NDBBD misconfiguration diagnostics added (missing/unmatched `-regset`, missing `-enable_time`; clean collective exit). `ems_solve(inmemory=, verbosity=)` exposed in teems-R (b98b31f; both default NULL = solver decides). `-march` image variants moved to 6.9 (2026-07-09). Op-list caching **measured and dropped** (2026-07-09): on the interactive-scale rig (GTAP-RE 10.5k eq, LU-1 Mmid, 1.16s wall) `formulas_execute` totals 0.04s across 15 per-step calls (3.4% of wall, parse+eval — caching removes only the parse slice); the real small-model cost is `jacobian_fill` (0.71s, 61% of wall), where equation parse interleaves with the mandatory per-step Jacobian refill — out of scope here and restructured wholesale by 6.2 anyway |

**Phase-5 findings** (details in `docs/solver-reference.md` §9):
factorization ≈63% of wall at scale; `formula_eval` ≈16% (compaction of
it measured speed-neutral — the wall is factorization); on small models
the per-step cost is `jacobian_fill`'s equation parse+refill (61% of a
1.16s interactive solve), not formulas (3.4%). `-inmemory` defaults: on
for LU/SBBD/DBBD (+7% SBBD, neutral elsewhere, no debris), off for
NDBBD on RAM-cost grounds (+69% RSS/rank for wall-neutral resident
factors on fast local scratch; slow/network scratch expected to favor
`inmemory=TRUE`, exposed in `ems_solve()`). Speed since refactor: SBBD
~5–8%; LU/DBBD neutral + cleanliness; NDBBD unchanged by design.

## Phase 6 — GEMPACK parity program (performance-vetoed)

**Work order (revised 2026-07-09):** (1) 6.1-scaled immediate language
gaps; (2) 6.5 netcut enforcement (E1 then E2); (3) 6.6 image fix batch
(rename + `-march` levels + golden refresh — a *needed fix*, not an
optimization: the published amd64 image is `-march=native`-built and
can SIGILL on older CPUs). Then 6.2 condensation and the rest.

| id | item | notes |
|---|---|---|
| 6.1 | **Immediate language gaps** (scaled back 2026-07-09 from the full parity inventory — that inventory is deferred, not dropped; the specs below are the documented limitations in `teems-manual/model_load.qmd`, whose workaround sections get deleted as each lands): (a) **set expressions** ✅ (2026-07-09, expanded to the full GEMPACK spec 10.1.1.1 per user guidance: UNION/INTERSECT/+/-/backslash over declared sets and quoted single elements, bracketed grouping, left-to-right; implied-SUBSET rules per the manual; '+' disjoint and '-' containment validated) — solver 0f59147 (recursive-descent set_expr_build; legacy pairwise encodings byte-identical, goldens 9/9; .audit/expr-test-kit exercises 3-term/brackets/INTERSECT/quoted forms), teems-R c3bc3bf (R/set_expr.R tokenizer+evaluator shared by parse and finalize; solve-equivalence test vs pairwise baseline passes), manual updated (teems-manual b880e49); (b) **set assignment** — `Set B # ... # = Set A;` without a duplicate Read statement; (c) **IF in Formula/Equation RHS** — `IF[c in MARG, VST(c,r)]`-style elemental conditionals (currently requires manually splitting the statement over subset/complement pairs). Each lands golden-gated + new test TABs exercising the form; update the manual and any teems-R-side validation in the same change | fixes real modeler friction now; full corpus inventory + upstream-docs refresh (GEMPACK, HSL, PETSc v3.25, MPICH into `teems-docs/`) resumes after |
| 6.2 | **Condensation** (substitute/backsolve/omit variables pre-factorization) — the largest speed lever in the project: shrinks the system that costs 63% of wall, and reduces netcut as a side effect. Design choice: teems-R-side TAB rewriting vs solver-side elimination. Must preserve the `.bin`/`.var` output contract (condensed variables recovered post-solve); interacts with closure/swap validation. NB also owns the small-model lever found 2026-07-09: `jacobian_fill` equation parse+refill is 61% of an interactive-scale solve — any condensation design that caches compiled equation programs (or eliminates re-parse) collects that win too | highest-leverage item overall |
| 6.3 | **Stepping methods** (in cost-benefit order): (a) Euler multi-step — near-free (Mmid driver with midpoint flag = 0), immediate parity gap; (b) verify/add Gragg terminal smoothing pass (manual §30.2 distinguishes gragg vs midpoint; TEEMS does Euler-start + midpoint steps + 2-4-8 Richardson); (c) Runge-Kutta: RK2/RK4/BoSha32/DoPri54 via Butcher-tableau driver — stages reuse the existing refill+solve primitive, compatible with all matrix methods; adaptive accept/reject must be MPI-collective (Allreduce the embedded-error norm) and the controller tuned conservatively at scale (a rejected step wastes a full set of block factorizations); (d) Adams-Bashforth — ONE factorization per step at order k (vs RK4's 4): directly attacks the factorization wall; RK startup; stability caveats; (e) implicit stiff methods (Radau/BDF/ode15s) — parked: performance-negative for smooth CGE paths; revisit only if adaptive RK step-size collapse demonstrates stiffness | RK gate: fixed-accuracy wall-time vs Gragg 2-4-8 |
| 6.4 | Statement completion per 6.1 ranking: full IF semantics (opcode skeleton exists: `OP_IF_*`), POSTSIM, others | each zero-cost when unused |
| 6.5 | **Netcut enforcement** (border size = interface-problem size, governs bordered-method viability): E1 validate in `ems_model()` (error/warn on TIME×large-set variables in inter-period links); E2 auto-rewrite minimal intertemporal proxies during TAB assembly; E3 solver-side element-level border classification — E3 after 6.2, since condensation changes what reaches the border | **priority raised 2026-07-09** — E1 cheap, E2 sweet spot |
| 6.6 | **Image fix batch** (absorbs 6.9, 2026-07-09): the current published approach is broken, not merely slow — amd64 base image built `-march=native` (SIGILL on CPUs older than the build host) while the solver's own objects carry no `-march` at all. One coordinated change: (a) `hsl` → `teems-solver` binary rename (touch points in project notes); (b) psABI-level image tags (x86-64-v2/v3, v4 only if it measures; arm64 stays portable `armv8-a`), the same level applied in both layers (base libs and solver makefile `OPT`); (c) teems-R selects the tag from host capability (`ld.so --help` lists supported levels; default v2) ; (d) golden refresh — required anyway since neither the rename nor the flag change preserves bit-identity, so re-anchor manifests once for the whole batch; per-level A/B on the bench rigs documents the v3/v4 speed effect | needed fix; coordinate teems-R + images; one retagging for everything |
| 6.7 | **HSL catalogue exploration** — assess additional HSL routines for speed/robustness across all solution × matrix methods; ranked shortlist with rationale in `teems-docs/hsl-catalogue-candidates.md` (top: MC64 pre-scaling, MA41/MA38 multifrontal alternative, MC33 bordered-block-triangular ordering vs netcut, HSL_MC79 Dulmage-Mendelsohn closure diagnosis). Retrieve specs/code for survivors; golden + A/B benchmark gates apply | user-requested 2026-07-05 |
| 6.8 | **Block-set auto-detection** (assessed 2026-07-06): drop the need for explicit `-regset`/`-enable_time`. Time side is already derivable — `(intertemporal)` TAB qualifier + `set_find_alltime`; the missing-flag case is even detected post-hoc (`sbbd_overrid` warning); NDBBD's timestep count self-derives from the alltime set (`ndbbdrank=ntime`). Regional side has no TAB qualifier and is name-convention only (teems-R hardcodes `REG`; unmatched names silently degrade to `allregset=-1`). All ordering code downstream is set-agnostic ("partition by set S"), so implement `-regset auto`: enumerate top-level non-intertemporal candidate sets, run the cheap pre-Jacobian ordering probe per candidate, score by netcut/border fraction, block count ≥ n_tasks, and block balance (the 5.7 stats fields), pick argmax or none (LU). Explicit flags stay as overrides (golden runs keep bit-identical behavior); log the choice in stats.json. Also fixes the silent-degradation footgun (warn when a named set is not found) | feeds `matrix_method` auto; pairs with 6.5/E3 |
| 6.9 | `-march` image variants | **folded into 6.6** (2026-07-09): reframed as a needed fix to the published image, not a performance item — see 6.6 for the batch |

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
