# TEEMS Solver — Development Roadmap

Working plan following the 2026 refactor (phases 0–4, logged below)
and phase 5's measured performance work. Completed items are never
deleted: phase-5 rows flip to **done** with findings attached, and
phases 0–4 are recorded in the completed-phases log at the end of this
file. `docs/solver-reference.md` describes the *resulting* system.

**Standing rules**
- Every change passes the golden suite (`.audit/verify.sh`, currently 14
  bit-identical gates: 5 method + 3 real-shock + 1 GEMPACK-orientation +
  4 subinterval gates spanning method × ranks × nsubints × inmemory,
  the LU one run in both -inmemory modes) and,
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
| 6.1 | **Immediate language gaps** (scaled back 2026-07-09 from the full parity inventory — that inventory is deferred, not dropped; the specs below are the documented limitations in `teems-manual/model_load.qmd`, whose workaround sections get deleted as each lands): (a) **set expressions** ✅ (2026-07-09, expanded to the full GEMPACK spec 10.1.1.1 per user guidance: UNION/INTERSECT/+/-/backslash over declared sets and quoted single elements, bracketed grouping, left-to-right; implied-SUBSET rules per the manual; '+' disjoint and '-' containment validated) — solver 0f59147 (recursive-descent set_expr_build; legacy pairwise encodings byte-identical, goldens 9/9; .audit/expr-test-kit exercises 3-term/brackets/INTERSECT/quoted forms), teems-R c3bc3bf (R/set_expr.R tokenizer+evaluator shared by parse and finalize; solve-equivalence test vs pairwise baseline passes), manual updated (teems-manual b880e49); (b) **set equality** ✅ (2026-07-09, GEMPACK spec 10.1.2.1: `Set B = A;` copies a previously defined non-intertemporal set and auto-generates the two SUBSET statements in both directions; the 13.3.1 `(Intertemporal)`/`(Non_Intertemporal)` conversion forms abort with a clear error) — solver: set_equality_build in tab_parse.c registers the mutual subsets (the legacy '='-copy branch did elements only; subset_map_build terminates on the cycle by adding self-entries), goldens 9/9, .audit/expr-test-kit TESTD probes formula addressing across the equality pair in both directions (sol.bin bit-identical); teems-R: bare `= A` no longer aborts, equality definitions keep their leading '=' so finalize's expression path (single-term) evaluates the alias, mutual subsets + cycle-safe visited-set transitive closure in tab_set_parse.R, e2e equation-quantifier solve-equivalence test passes; manual workaround section replaced with the supported spec; (c) **IF in Formula/Equation RHS** ✅ (2026-07-10, GEMPACK spec 11.4.5–11.4.7; rewrite-over-surgery: GEMPACK itself translates index-IN-set IFs to conditional sums, and the GTAP authors' hand adaptations are the exact target shapes) — teems-R `R/tab_if_rewrite.R` rewrites top-level additive `IF[cond, expr]` terms at parse time: Formulas → base statement with IF terms dropped + one accumulate formula per term over a narrowed domain (`idx in S` / `idx="elem"` → synthesized intersection set `IFS* = S & range`, whose all-INTERSECT implied subsets keep refs valid only over S in range per 11.4.7 rule 2; `coefref op const` → conditional quantifier); Equations → comparison terms become synthesized 0/1 indicator coefficients `IFC*` distributed over the value's top-level terms (the MAKESUNIT/E_qca shape), and one membership/element term splits the equation into complementary-domain equations `<name>A`/`<name>B` (indicator distribution would break rule 2's deemed range there); compound AND/OR/NOT, nested IF, coefref-vs-coefref comparisons, and >1 membership split per equation abort clearly. Solver: formula.c conditional-quantifier constant parse fixed — `atof` ran on the token still holding the operator, so ANY constant other than 0 silently read as 0; now parses past the operator (plus null-terminating the scalar-condvar strncpy), goldens 9/9, kit TESTE probes split VENDWMS at >1e6 / <=1e6 and recompose bit-exactly. Acceptance (test-ems_solve.R): original gtapv7.tab IF statements grafted over the shipped adaptations — VCB/VXW formulas solve bit-identically; E_qca/E_pca/E_pds equations bit-identical on common outputs (IFC1/IFC2 extra, consistent with MAKESUNIT appearing in base outputs); element probe `2 + IF[r="chn", 1]` value-verified under big3 aggregation; manual workaround section replaced with the supported spec. Remaining gaps → 6.4: IF nested inside expressions (needs solver-side evaluation), compound conditions, coefref-vs-coefref comparisons | fixes real modeler friction now; full corpus inventory + upstream-docs refresh (GEMPACK, HSL, PETSc v3.25, MPICH into `teems-docs/`) resumes after |
| 6.2 | **Condensation** (substitute/backsolve/omit variables pre-factorization) — the largest speed lever in the project: shrinks the system that costs 63% of wall, and reduces netcut as a side effect. Design choice: teems-R-side TAB rewriting vs solver-side elimination. Must preserve the `.bin`/`.var` output contract (condensed variables recovered post-solve); interacts with closure/swap validation. NB also owns the small-model lever found 2026-07-09: `jacobian_fill` equation parse+refill is 61% of an interactive-scale solve — any condensation design that caches compiled equation programs (or eliminates re-parse) collects that win too | highest-leverage item overall |
| 6.3 | **Stepping methods** (in cost-benefit order): (a) Euler multi-step — near-free (Mmid driver with midpoint flag = 0), immediate parity gap; (b) verify/add Gragg terminal smoothing pass (manual §30.2 distinguishes gragg vs midpoint; TEEMS does Euler-start + midpoint steps + 2-4-8 Richardson); (c) Runge-Kutta: RK2/RK4/BoSha32/DoPri54 via Butcher-tableau driver — stages reuse the existing refill+solve primitive, compatible with all matrix methods; adaptive accept/reject must be MPI-collective (Allreduce the embedded-error norm) and the controller tuned conservatively at scale (a rejected step wastes a full set of block factorizations); (d) Adams-Bashforth — ONE factorization per step at order k (vs RK4's 4): directly attacks the factorization wall; RK startup; stability caveats; (e) implicit stiff methods (Radau/BDF/ode15s) — parked: performance-negative for smooth CGE paths; revisit only if adaptive RK step-size collapse demonstrates stiffness | RK gate: fixed-accuracy wall-time vs Gragg 2-4-8 |
| 6.4 | Statement completion per 6.1 ranking: full IF semantics (opcode skeleton exists: `OP_IF_*`), POSTSIM, others | each zero-cost when unused |
| 6.5 | **Netcut enforcement** (border size = interface-problem size, governs bordered-method viability): E1 validate in `ems_model()` ✅ (2026-07-10, teems-R: `.check_netcut()` in `R/chk_netcut.R`, wired into `.process_tablo` after maths parse — mirrors the solver's border rule (jacobian.c `var_inter`: any lead/lag reference borders the variable's entire element block) by scanning Equation definitions for variable references carrying an index offset and warning (`model_wrn$netcut_inflation`) when such a variable is declared over ≥2 non-`(intertemporal)` sets — the minimal-proxy shape (one non-time dim, e.g. `kb(REG,TIME)`) passes silently and all shipped models are warning-free; set sizes are unknown at `ems_model()` time so dimension count is the structural signal; warning carries the proxy recipe, manual §netcut added, snapshot + no-warning tests in test-ems_model.R); E2 auto-rewrite minimal intertemporal proxies during TAB assembly ✅ (2026-07-10, teems-R: `.rewrite_tab_netcut()` in `R/tab_netcut_rewrite.R`, runs right after the IF rewrite in `.process_tablo` — lead/lag variable references that fix ≥1 dimension to a quoted element (the only statically-guaranteed reduction; corpus scan found *zero* violating models — INT_SCET/GTAP-RE/GTAP-INT all link periods via minimal REG×TIME variables already, so E2 guards new models) get a synthesized proxy variable `NCV*` over the remaining dims + linking equation `E_NCV*` inserted after the source declaration, the lead/lag moved onto the proxy, `orig_level` qualifiers dropped, other qualifiers inherited, proxies cached per (variable, element-slice), always-on with a cli inform listing rewrites; `.check_netcut()` skips synthesized proxies; full-set references remain E1-warned. Verified: rewritten slice-lead GTAP-RE graft solves bit-identically to a hand-proxied reference (SBBD-2, common outputs identical, proxy values == hand-proxy values; test-ems_solve.R) and netcut drops 207→36 on the probe (stats.json A/B with rewrite disabled via assignInNamespace). Solver untouched — quoted-element refs in equations were already supported, e.g. GTAP-INT `qo("capital",r,t)`); E3 solver-side element-level border classification — E3 after 6.2, since condensation changes what reaches the border | **priority raised 2026-07-09** — E1+E2 done, E3 waits on 6.2 |
| 6.6 | **Image fix batch** (absorbs 6.9, 2026-07-09): the current published approach is broken, not merely slow — amd64 base image built `-march=native` (SIGILL on CPUs older than the build host: the published binary carries 42,731 AVX/AVX-512 instructions from static `-march=native` PETSc) while the solver's own objects carry no `-march` at all. One coordinated change: (a) `hsl` → `teems-solver` binary rename ✅ (2026-07-10: makefile target, expedited/full Dockerfiles copy the new name and add a one-release `hsl` transition symlink in the runtime image, docs/README/verify.sh paths updated, teems-R invokes the new path (ee76e2a; cmd-string snapshots re-anchored); rename verified output-neutral — all 9 goldens bit-identical with the renamed binary; `teems:dev` built with rename + flag plumbing, GTAP-RE SBBD-2 solve OK; local `teems:latest` re-pointed at it, old image kept as `teems:pre-march-backup`. NB the pre-existing GTAP12 dynamic cross-method failure now manifests as numeric disagreement instead of a bound abort — new flags shift numerics on that borderline data; static sibling now passes; still owned by feature/GTAP12a); (b) psABI-level image tags ✅ **plumbing done** (2026-07-10): `MARCH` build arg in base/full Dockerfiles (amd64 default `x86-64-v2`, arm64 `armv8-a`, never `native`); resolved flags recorded in `${BUILD_DIR}/archflags` and consumed by the expedited layer (`make OPT="-Ofast $ARCH_FLAGS"`, HSL configures `-O2 $ARCH_FLAGS` — deliberately not `-Ofast` to keep MA48 pivoting semantics), so base/solver levels cannot mismatch; expedited gained `BASE_IMAGE` build arg for variant selection; MPICH now `--enable-fast=all,O3` + `MPICHLIB_*FLAGS` (arch flags kept out of the mpicc wrappers per install guide), PETSc `--with-memalign=64 --with-cxx=0 --with-fortran-bindings=0` (docs reviewed from teems-docs/mpich + teems-docs/petsc clones). Validated end-to-end: base v2 build clean, expedited picks up flags, binary has 0 ymm/zmm instructions (vs 42,731 in current image), GTAP-RE SBBD-2 deploy+solve OK via `docker_tag="v2-test"`. NB MA48's dgemm/dtrsm resolve to Debian OpenBLAS (runtime CPU dispatch) — BLAS3 already optimal per host; `OPENBLAS_NUM_THREADS=1` A/B queued for the batch. Maintainer guide at `docs/image-publishing.md` (variant matrix, buildx/publish/manifest commands, validation checklist, benchmark gates); (c) teems-R selects the tag from host capability (`ld.so --help` lists supported levels; default v2); (d) golden refresh — required anyway since neither the rename nor the flag change preserves bit-identity, so re-anchor manifests once for the whole batch; per-level A/B on the bench rigs documents the v3/v4 speed effect. REMAINING: rename (a), teems-R tag selection (c), golden refresh + A/Bs + actual publish (d) | needed fix; coordinate teems-R + images; one retagging for everything |
| 6.7 | **HSL catalogue exploration** — assess additional HSL routines for speed/robustness across all solution × matrix methods; ranked shortlist with rationale in `teems-docs/hsl-catalogue-candidates.md` (top: MC64 pre-scaling, MA41/MA38 multifrontal alternative, MC33 bordered-block-triangular ordering vs netcut, HSL_MC79 Dulmage-Mendelsohn closure diagnosis). Retrieve specs/code for survivors; golden + A/B benchmark gates apply | user-requested 2026-07-05 |
| 6.8 | **Structural partition detection** (assessed 2026-07-06 as "block-set auto-detection"; scope deepened 2026-07-12 per user directive: the solver must be *set-agnostic* — the structure of the equation system, not set names or invariant qualifiers, determines how each set is handled, and `-regset`/`-enable_time` are slated for removal) | **done** 2026-07-12 — both ordering dimensions now derive from the equations when the transitional flags are absent: (a) *chain dimension* — `chain_refs_scan` walks the equation statements (same preprocessing as `equation_order_read`) counting lead/lag offset references per declared dimension set; `chain_set_select` aggregates to top-level sets and picks the most-referenced (the `(intertemporal)` qualifier only licenses the syntax and is cross-checked: qualified-but-never-offset-referenced sets get no chain treatment, offsets on unqualified sets warn, multi-set offsets warn and majority wins). Confirmed name-agnostic against the uom corpus (`teems_tabs/uom_based/gtap_ebm.tab` chains on `allweek`, not `alltime`; its lat/long spatial grids couple via adjacency-coefficient sums, which the probe prices as border). (b) *diagonal-block partition* (DBBD/NDBBD) — every set with ≥2 elements outside the chain family indexing a variable dimension (directly or via subset, so subregional-style sets compete) is probed: `equation_order_read` + the extracted first counting pass (`block_var_count`, now shared by probe and live ordering) measure the real border/block sizes per candidate; viable = ≥ n_tasks nonempty blocks and border < half the system; smallest border wins, 2% near-ties broken by block balance (integer-compared — every rank reaches the same answer without communication). KEY SCORING LESSON: in nested (NDBBD) mode the per-chain-block interface column must count as border — it is the local border the interface solve pays for; scoring it intra made all candidates tie at netcut 27 on GTAP-RE and balance picked `acts` (half the system in the interface), which segfaulted the solve. With the fix `reg` wins 2.9% vs `demd` 15.8% / `comm` 34.9% / `acts` 46.6%. Acceptance: NDBBD, DBBD-intertemporal, DBBD-static, SBBD golden configs rerun with the flags REMOVED reproduce their manifests bit-identically (`.audit/accept_structural.sh`); `-regset auto` alias likewise; 14/14 goldens (explicit flags) bit-identical; LU-without-flags keeps no chain treatment (bit-compat) and gets a reworded structure-ignored warning (re-gated on `alltimeset<0`, was `-enable_time`-based). stats.json v2: `chain_source`/`partition_source` (`explicit`/`structural`/`none`), `chain_set`/`partition_set`, full `partition_auto` candidate table — the structural evidence for method auto. **`-solmed NoSol` verified as the structure-probe pathway**: 0.5 s wall on the GTAP-RE golden, runs full detection, writes the complete 13-candidate table to stats.json, no solve — the concrete mechanism by which `matrix_method` auto can later consult structure before choosing a method. Method-vs-structure aborts added: SBBD without a chain dimension previously handed HSL_MP48 zero blocks and finished with exit 0 and NO solution (silent-failure footgun, unreachable from teems-R which validates enable_time) — now a clean collective abort (verdict Bcast from rank 0, since HSL ranks≠0 never resolve alltimeset). Unmatched `-regset` names now warn + fall back to detection. NDBBD/DBBD abort with remedies when required dimensions resolve to none. NB the `acts`-partition segfault is a pre-existing NDBBD fragility for pathological partitions (never selectable now, but reproducible with explicit `-regset acts`) — noted under standing investigations | feeds `matrix_method` auto (solver-side method-vs-structure validation + `NoSol` probe run pathway); flag REMOVAL + teems-R dropping `-regset REG`/`-enable_time` awaits a golden re-anchor slot (fold into 6.6(d) or later); value-aware element-level classification (banded adjacency coefficients ≈ partitionable) pairs with 6.5/E3 |
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

- **RESOLVED 2026-07-11 — multi-subinterval bug was TWO defects** (fixed
  8778838 + the hsl_mp48d.patch revision; 14-gate golden suite green):
  (1) *xc0 unguarded free* (ours, phase 5.3 inmemory refactor,
  2026-07-05): `free(xc0)` in `solve_modified_midpoint` sat outside the
  `if(!inmemory)` guard while its reload stayed inside — under
  `-inmemory` the cumulative multiplier vanished for the sol>0
  extrapolation passes of every subinterval after the first, so the
  Richardson combination blended pass-inconsistent quantities (the
  exact 23/45 signature) and composed output violated the exogenous
  shock (pfactwld+1 → 0.5111). Symptoms varied by method/env: silent
  half-shock (LU/SBBD), a wild interface element (SBBD-2 on GTAP12a),
  bound aborts, ns=3 false singularity. Fix: free moved inside the
  guard; ns=1 outputs proven bit-identical.
  (2) *CGLOB overrun* (legacy, in the original project patch to
  HSL_MP48): the memory-reduction patch re-dimensioned
  `CGLOB(NEQ,NBLOCK)` → `CGLOB(MAXSBCOLS,NBLOCK)` with MAXSBCOLS
  estimated from block column pointers — but the border block's true
  column span (NCOL 5249 vs bound 1399 on the probe) exceeds any such
  estimate, so *every* SBBD/DBBD solve wrote ~15KB past the array
  (valgrind: invalid write, hsl_mp48d.f90:1293). Results stayed
  correct only because the overrun addresses were read back before
  reuse; ns≥2's allocation churn broke that and crashed
  (munmap_chunk) in some environments. Fix: hsl_mp48d.patch now
  measures the true per-block span in a pre-pass and sizes CGLOB from
  it; valgrind-clean; all goldens bit-identical (the round-trip had
  been output-neutral).
  Golden suite widened to 14 gates: + nsub-lu (Mmid ns=2 LU-1, run in
  BOTH -inmemory modes against one manifest — modes proven
  bit-identical), nsub-sbbd (ns=2 SBBD-2), nsub-ndbbd (ns=3 NDBBD-2),
  nsub-dbbd (static Mmid ns=2 DBBD-2); GTAP11c deployments with real
  numeraire shocks (.audit/deploy_golden_nsub.R). Coverage now spans
  solution methods × matrix methods × ranks × subintervals ×
  inmemory. teems-R test-ems_solve.R 40/40 — the "GTAP12 data class"
  cross-method failures were this bug, not the database.

## Standing investigations

- ~~`-ndbbd_bl_rank`/`-nestfile` per-block rank overrides~~ RESOLVED
  (2026-07-11): removed. Investigation corrected the semantics — "rank"
  is the NUMERICAL rank of each per-time interface block (the local
  border of the nested DBBD), not an MPI rank assignment. The options
  were a vestige of the research code's two-run workflow (`-presol`
  preparation run wrote per-block interface files recording each
  block's true rank from the rank-revealing MA51 factorization; a
  second run consumed them, with the CSV carrying the ranks between
  runs). Since the first TEEMS commit the presolve pass is hardcoded
  inline (`presol=1` at every NDBBD call site) and
  `ndbbd_order_presolve` unconditionally overwrites every array entry
  with min(nrow,ncol) before anything reads it — the CSV values could
  never influence a solve, on any code path; `-presol` was equally
  inert (parameter overwritten at all use sites; it was also read into
  an uninitialized bool when the flag was absent). The "fold into
  auto-calibration" idea is already the shipped behavior: the in-solve
  presolve measures each interface block's true rank and shrinks/
  permutes accordingly. Removed `-ndbbd_bl_rank`, `-nestfile`,
  `-presol`, `csv_read_ints()`, and the dead `./ndbbd_drank.csv`
  fallback; teems-R no longer passes `-presol 1` (compatible both
  directions — old binaries ignored the value anyway). Follow-up in the
  same sweep: `-nesteddbbd` removed too — it only ever restated
  `-matsol 3` (teems-R derived it from `matrix_method == "NDBBD"`; a
  mismatch in either direction was simply broken), so the solver now
  derives `nesteddbbd = (matsol==MM_NDBBD)`. NB solver-side this one
  is NOT ignorable by old binaries: an old solver given a new teems-R
  command runs NDBBD un-nested — runtime images must carry a solver
  at/after this change (teems:dev rebuilt).

- NDBBD-vs-SBBD crossover at scale (medium ACTS, 8–16 ranks, cluster
  hardware) — dev-box sweep found no NDBBD win region.
- NDBBD segfault under pathological partitions (found 2026-07-12 via
  6.8): `-matsol 3 -regset acts` on the GTAP-RE golden (interface
  column ≈ half the system) dies in SEGV during the solve. Structural
  detection now scores such partitions out, but the crash class
  predates 6.8 and is reachable with explicit flags; diagnose if NDBBD
  hardening becomes a priority.
- Stiffness monitoring once adaptive RK exists (see 6.3e).
- Ordering freedom in `equation_order_read_nested` (legacy note: "can
  try different varreg"); direct vs MC66 ordering netcut trade-off.
