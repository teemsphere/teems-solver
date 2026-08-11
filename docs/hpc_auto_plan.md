# HPC benchmark matrix & auto-setting calibration plan

Status: DRAFT 2026-08-11 (pre-instance). Target: single node, 32 CPU /
128 GB RAM. Consumes the gates other ROADMAP rows parked "on the HPC
matrix": bordered-method parallelism beyond 2 ranks, the RK-vs-Gragg
fixed-accuracy ranking, the `-fastrefac` default flip, TSan at
production thread counts. Adds the argument-by-argument auto-setting
program: the end state is that `ems_solve(cmf_path)` with no other
arguments does the right thing on both a laptop and this instance, and
every auto decision is reported (existing `matrix_method = "auto"`
info-message precedent) and recorded (stats.json `options`).

## 1. Ground rules

- **Interleaved A/B** (ab_*.sh precedent): alternate the compared
  configs run-by-run so thermal/cache drift cancels; 3 reps minimum,
  report median + spread. One benchmark at a time on the instance.
- **Truth references**: accuracy claims measured against f64 deep
  Gragg (the 6.3(f3) binary; `teems-solver-f64` Gragg 128-256-512 or
  deeper), never against another f32 run. Max-metric attribution rule
  stands: name the variables before judging (wev/ev/cnt* cancellation
  class rides level-scaled aggregates).
- **Instrumentation**: `/usr/bin/time -v` (wall, maxrss),
  sol.stats.json (netcut, border sizes, block counts, effective
  options), the solver accuracy summary, `.acc` where embedded
  methods run. Per-rank RSS via `mpiexec -n N /usr/bin/time -v` is
  aggregate-only; add a `smaps_rollup` sampler for the memory-per-rank
  curves.
- **Core accounting**: establish early whether the 32 CPUs are
  physical or SMT (lscpu); run the n_tasks ladder to 32 either way but
  interpret the 16→32 step against the topology.
- **Bit-reproducibility caveat**: `maxthreads > 1` changes reduction
  order — accuracy comparisons at >1 thread use tolerance bands, not
  bit equality. Re-run the TSan/Archer slice at the production thread
  counts chosen here before any threads-default change ships.

## 2. Model ladder

Existing rigs (checked into .audit, run out-of-the-box): bench-static
(~200k eq), bench-static-XL, bench-static-XXL, bench-inter,
bench-inter-L, plus -cond variants for the condensation interaction.
Add two instance-scale deployments built through teems-R from the full
GTAP database:

- **S-full**: static GTAPv7, full disaggregation (141 reg × 65
  sectors) — lands in the multi-million-equation band where the
  current auto DBBD thresholds (2e6 / 1.5e6 @ ≥100 reg) were
  extrapolated, and exercises the uncalibrated 1.5–2M × <100-region
  band via intermediate aggregations.
- **I-long**: GTAP-RE with the longest time horizon the data
  supports, sized to push toward the 128 GB ceiling at 8–16 ranks —
  the rig for memory-per-rank scaling and the inmemory/scratch
  decision at scale.

Riders on selected cells: one mapped model (bordered-map kit model)
and one complementarity model (comp kit) to confirm the 6.14 lifts
hold at rank counts > 2. Shock-magnitude ladder for the subinterval
and precision programs: pfactwld uniform {1, 10, 40, 100, 150}.

## 3. Core program: n_tasks × matrix_method scaling

The matrix everything else rides on:

- Methods {LU(n=1 only), SBBD, DBBD, NDBBD} × n_tasks {1, 2, 4, 8,
  16, 24, 32} × model ladder, Gragg 2-4-8, nsubints 1.
- Record wall, aggregate + per-rank RSS, netcut/border/block stats.
- Products:
  (a) parallel-efficiency curves and the knee per method × size —
  note ndblock caps useful ranks (35 region blocks static / 21 time
  blocks on current rigs; S-full raises the static cap to 141);
  (b) memory-per-rank model (HSL workspaces and border replication) →
  the `mem_budget` term of the n_tasks auto rule;
  (c) recalibrated matrix_method crossovers at real rank counts —
  the current auto was fit at ≤4 tasks;
  (d) oversubscription penalty (n_tasks 32 + n_threads 2, and 48
  ranks) so the auto can refuse gracefully.

## 4. Parked gates folded into the matrix

- **RK-vs-Gragg fixed-accuracy ranking** (6.3 gate): tune each of
  {Gragg, Euler, RK4, BoSha32, DoPri54 ±adaptive} to hit relative
  error targets {1e-4, 1e-6} vs f64 truth on {static-XL, inter-L,
  uniform-40 large shock} at n_tasks {1, 8, 32}; rank by wall. This
  is the calibration table the solution_method auto reads.
- **fastrefac flip** (6.13 condition): A/B `-fastrefac 1` at ranks
  {8, 16, 32} × {LU@1, SBBD, DBBD, NDBBD} × sizes; accuracy-neutral
  within f32 noise + no MA48 refactor-failure regressions → flip the
  default ON for LU/SBBD/DBBD (NDBBD separately: interaction with its
  disk-staged blocks). A flipped default is the purest auto: the flag
  disappears for users.
- **Precision at scale** (6.3(f3) follow-on): f32 vs f64 wall/RSS on
  S-full and I-long at {8, 32} ranks. Local result was RSS +7.4–8.6%,
  wall 0–7%; confirm the delta stays single-digit when coefficient
  arrays are tens of GB and memory bandwidth is contended at 32
  ranks. Also map the f32 accuracy floor on instance-scale models
  (deep-ladder ladders vs f64 truth) — this number parameterizes the
  precision auto threshold.
- **TSan at production threads**: re-run the Archer slice at the
  n_threads values the auto will actually pick.

## 5. Argument-by-argument auto assessment

Verdict key: AUTO = worth an automatic setting for an unsophisticated
user; POLICY = better served by a changed/confirmed default; EXPERT =
leave manual; N/A = workflow flag, no tuning content.

| Argument | Verdict | Decision signal | Calibration tests |
|---|---|---|---|
| solution_method | AUTO (via accuracy profile) | requested accuracy tier + shock magnitude + model class | §4 fixed-accuracy ranking; failure behavior on stiff/large shocks |
| matrix_method | AUTO (exists; recalibrate) | intertemporal flag, system size, n_reg — extend with netcut/border size from a prior stats.json or deploy metadata | §3 crossovers at 8–32 ranks; validate chosen-vs-best gap ≤10% across the ladder |
| n_tasks | AUTO (new; HPC-critical) | ndblock, system size, host cores, memory-per-rank estimate | §3 knees + per-rank RSS; rule `min(knee(size, method), ndblock, cores, mem_budget/rank_rss)`; oversubscription penalty |
| n_threads | POLICY (likely stay 1) | idle cores after n_tasks chosen | tasks×threads grid at fixed product 32 (32×1, 16×2, 8×4, 4×8) + `-smllthreads`; only promote if a grid cell wins reproducibly; TSan gate first |
| precision | AUTO | accuracy tier below f32 floor (~3e-6 rel), deep ladders, large shock × many subintervals | §4 precision-at-scale; floor map on S-full/I-long; post-solve advisory when the accuracy summary is floor-limited |
| steps | AUTO (with method) | accuracy tier; `.acc`/accuracy-summary feedback | step-count ladders on production rigs vs f64 truth per shock class; per-method steps table the profile reads |
| n_subintervals | AUTO (interim R heuristic; real fix = 6.3(f1) solver-side adaptive subintervals) | max shock magnitude from the shf, method | shock ladder × nsubints {1,2,4,8}: success/accuracy/wall; threshold table |
| adaptive | POLICY | method chosen by auto | adaptive-vs-fixed wall at equal achieved accuracy; retry-rate at 32 ranks (collective accept/reject validated at 2) — then default "yes" for embedded methods |
| eps_tolerance | AUTO (derived) | accuracy tier | epstol→achieved-error transfer curve at scale (toy kit shows `.acc` over-reads 2–30×; need the production mapping) |
| max_retries / retry_adjust | EXPERT (defaults confirmed) | — | stiff large-shock sweep to confirm defaults don't bind |
| laA / laD / laDi | AUTO (in-solver grow-and-retry; see note below the table) | MA48 reports the needed size itself: INFO(1) = −3 with the suggestion in INFO(3)/INFO(4), and the measured fill ratio is available even on success | port the fastrefac-LU grow loop (solve_drivers.c) to the fatal −3 sites; starved-la kit legs per method (goldens must stay bit-identical — LA is workspace, not pivot policy); required-ratio distribution across the ladder × methods for defaults |
| inmemory | POLICY (solver auto exists; widen at 128 GB) | solver memory estimate vs host RAM | inmemory 0/1 × method × size wall/RSS; estimator-vs-measured audit; tmpfs scratch A/B |
| verbosity | N/A | — | — |
| assertions / range_test_* / postsim | N/A (semantic switches; defaults already conservative) | — | — |
| complementarity | EXPERT (per-component logic already self-adjusting) | — | one comp rider cell at 8+ ranks to confirm scale behavior |
| pre_probe | AUTO (on-failure) | solve failed structurally-singular → auto-probe and re-report with named defects; optionally always-on below ~1e5 eq | probe cost curve vs size (claimed negligible <1e5, tens of seconds ~1e6 — verify at S-full scale, bordered + mapped included) |
| append_args: fastrefac | POLICY (flip default, §4) | — | §4 |
| append_args: nsbbdblocks / withmc66 | EXPERT (possible later auto from ndblock/ranks) | — | SBBD block-count sweep vs ranks, only if SBBD wins cells in §3 |
| append_args: cntl_3 / cntl_6 / tempdir / nowrites / gpzerodivide | EXPERT / N/A | — | — |
| suppress_outputs / terminal_run | N/A | — | — |

### la* sizing note (survey 2026-08-12)

All three parameters are the same quantity at different factorization
sites — workspace as a % of that matrix's nonzeros
(`lasize = ceil(la/100 * nz)`): laA at the main-system/diagonal-block
analyses (LU, SBBD, DBBD/NDBBD blocks), laD at the DBBD/NDBBD final
interface system, laDi at NDBBD's intermediate per-block interfaces.
The infrastructure for full auto-sizing already half-exists:

- The persistent LU kernel (fastrefac path) already does in-run
  reallocate-and-retry on INFO(1) = −3, growing to
  max(2×, MA48's suggested size) and logging the equivalent -laA
  (solve_drivers.c). MP48 (the SBBD block path) grows LA internally
  on FLAG = −3 as well (`LA = MAX(ICONTROL(6)*LA, INFO_MA60(3,...))`).
- Every remaining −3 site (one-shot LU, DBBD/NDBBD block and
  interface solves — 16 fatal MA48 sites in hsl_kernels.f90) aborts
  today with "increase -laA"; the −3/suggested-size return protocol
  (INSIZE(5)/INSIZE(6)) is already defined and just needs wiring to
  those callers. Growth at collective sites must broadcast the redo
  decision (sbbd_fastrefac precedent; rank_hsl-as-root gotchas
  apply), and NDBBD's disk-staged presolve blocks need a check that a
  regrow between stages doesn't invalidate staged sizes.
- On SUCCESS the analyse reports the measured fill ratio
  (INFO(3)/NE, currently printed at verbosity 2 only) — record
  la*_used per site in stats.json so repeat runs of a deployment
  warm-start at the measured ratio, and the HPC sweep can set
  cold-start defaults at the observed P95 instead of today's
  guessed 300/500/200.

With in-solver grow, la* stop being user decisions entirely:
undershoot costs a redone analyse (seconds), overshoot costs
transient RSS. The args stay as expert overrides. This supersedes
the earlier idea of an R-side catch-and-escalate loop, which would
pay a full pipeline re-run per retry.

## 6. Proposed user-facing surface

Not 20 knobs with auto values — two:

- `accuracy = c("standard", "exploratory", "high")` — resolves
  solution_method + steps + adaptive + eps_tolerance + precision +
  n_subintervals from the §4/§5 calibration tables. "standard" =
  today's Gragg 2-4-8 tier; "exploratory" = Johansen/coarse (model
  development loops); "high" = the fixed-accuracy winner at 1e-6 with
  precision="double" when the target sits below the f32 floor.
- `resources = "auto"` (or explicit n_tasks/n_threads as today) —
  resolves n_tasks, n_threads, inmemory, la* initial values from host
  inspection + deploy metadata.

Explicit arguments always override their resolved value (current
`matrix_method="auto"` semantics). Every resolved value is reported at
run time and recorded in stats.json options → model_diagnostics, so a
run is reproducible from its record regardless of how it was invoked.

## 7. Sequencing

1. Instance bring-up: topology audit, container runtime, rig sync,
   S-full/I-long deployment builds, smoke the harness (§1).
2. §3 core matrix (everything reads it).
3. §4 parked gates (reuse §3 cells where configs coincide).
4. §5 single-argument programs that need dedicated sweeps: la*
   distribution, subinterval shock ladder, epstol transfer, probe
   cost curve.
5. Fit + implement autos in teems-R (auto_method.R precedent:
   constants at top, calibration provenance in comments); flip
   defaults (fastrefac, adaptive) solver-side where §4 clears them.
6. Validation pass: `ems_solve(cmf)` + profiles on every ladder
   model; chosen-vs-best wall gap ≤10%; suite + goldens + kits green.

Exit criteria: an unsophisticated user on this instance gets ≥90% of
tuned performance and the correct accuracy tier with zero solver
arguments; every parked gate has a recorded verdict in its ROADMAP
row.

## 8. Open questions (user)

1. Single node only, or should the matrix reserve cells for eventual
   multi-node MPI (changes the n_tasks auto ceiling and NDBBD's
   presolve-file traffic story)?
2. Container runtime on the instance: docker as here, or
   singularity/apptainer (image conversion step)?
3. Is the full GTAP database licensed/available on the instance for
   S-full, or should S-full be built locally and shipped as a rig?
4. Default-flip appetite: ship cleared flips (fastrefac, adaptive,
   generous la* at high RAM) as they clear, or hold everything for
   one calibrated release with the profile surface?
