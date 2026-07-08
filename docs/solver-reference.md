# TEEMS Solver — Technical Reference

The TEEMS solver (`hsl` binary) is a parallel Computable General
Equilibrium (CGE) solver in the GEMPACK tradition: it reads a TAB-language
model description, base data, a closure, and shocks; linearizes the model;
solves the resulting first-order system with sparse direct methods; and
writes updated data and percentage-change solutions. It implements the
bordered block-diagonal ordering methods of Ha & Kompas (2016, SBBD) and
Kompas & Ha (2019, NDBBD) on top of PETSc, MPI, OpenMP, and the HSL
sparse linear algebra libraries (MA48, MA51, MC66, HSL_MP48).

Authored by Pham Van Ha and Tom Kompas; restructured and hardened in the
2026 refactor (see `src/NAMING.md` for the symbol lineage and
`BUILDING.md` for the build).

References used throughout:
- **[GM]** GEMPACK Manual, Horridge et al. 2024
- **[HK16]** Ha & Kompas, *Economic Modelling* 52 (2016) — SBBD
- **[KH19]** Kompas & Ha, *Economic Modelling* 80 (2019) — NDBBD
- **[D20]** Dixon (ed.), *Handbook of CGE Modeling*, ch. 20

---

## 1. Position in the TEEMS pipeline

teems-R prepares a self-contained run directory and invokes the solver in
a container:

```
docker run --rm --mount type=bind,src=<run dir>,dst=/opt/teems teems:<tag> \
  /opt/teems-solver/lib/mpi/bin/mpiexec -n <n_tasks> \
  /opt/teems-solver/solver/hsl -cmdfile /opt/teems/<model>.cmf <options>
```

The run directory contains:

| file | producer | role |
|---|---|---|
| `<model>.cmf` | `ems_deploy()` | command file: file bindings, tabfile/closure/shock paths |
| `<model>.tab` | `ems_model()` | TAB model description |
| `<model>.cls` | `ems_model()` | closure (exogenous variable list) |
| `*.shf` | `ems_shock()`/deploy | shock statements |
| `GTAPSETS/GTAPDATA/GTAPPARM.txt` | `ems_data()` | set elements and base-data arrays (text) |
| `out/` | solver | outputs (below) |

The CMF is a sequence of statements read by `cmf_read()`:
`file iodata "<logical>" "<path>";` binds model file names to paths;
`outdata` statements name per-set/per-coefficient CSV outputs; `soldata`
names the solution binary prefix; `tabfile`/`closure`/`shock` point at
the model inputs. All paths are absolute container paths when deployed
by teems-R.

Solver outputs, consumed by `ems_compose()`:

- `out/sets/*.csv`, `out/coefficients/*.csv` — post-simulation set
  listings and updated coefficient values.
- Solution binaries (prefix from `soldata`, default `solution`):

| file | contents |
|---|---|
| `.bin` | `nvarele × double` — percentage-change (or change) solution per variable element |
| `.var` | `nvar × array_def` — variable declarations (name, dims, set ids, strides) |
| `.set` | `nset × set_def` — set definitions |
| `.sel` | `nsetspace × set_element` — set elements with superset positions |
| `.mds` | 4 × long: `nsetspace, nvar, nvarele, nset` |
| `.stats.json` | per-run ordering statistics: system size, method, `netcut`, border sizes, per-block variable/equation counts (null/empty when no bordered ordering was built). Written before the solve, so failed runs still record their ordering; feeds `matrix_method` auto-calibration |

The structs are written raw; teems-R's `parse_solution.cpp` mirrors their
layout, so `array_def`/`set_def`/`set_element` are an ABI shared with
teems-R and must not be re-ordered or resized casually.

## 2. Model language

The solver reads a GEMPACK-style TAB subset (statement syntax per [GM]):

- `set` / `subset` — including unions (`+`), differences (`-`), and
  intertemporal set declarations (`(intertemporal)`), which mark the time
  dimension used by the bordered orderings.
- `coefficient` / `variable` — levels or percentage-change quantities
  (`(change)` for change variables; `(default=levels)` handled by
  `variables_read_defaults()`); optional bounds (`GE/GT/LE/LT`,
  `enum bound_type`).
- `read` — from the CMF-bound data files.
- `formula` — `(initial)` / `(always)` semantics as in [GM];
  `zerodivide` defaults honored by `tab_next_statement_resolved()`.
- `equation` — linearized equations; the left/right sides are compiled
  per equation block by `jacobian_fill()`.
- `update` — product-form and explicit updates (`updates_apply()`,
  `updates_apply_product()`); the modified-midpoint variant applies
  `substep_base + 2·Δ` (see §5).
- `swap` (closure) and shock statements — parsed by `closure_read()` /
  `shocks_read()` into per-element `closure_entry` records.

Shock statement values follow GEMPACK-standard ordering (first
subscript varies fastest; fixed dimensions collapsed to size 1).
Startup notes remind that intertemporal variables should be declared
with minimal dimensionality to keep the border (netcut) small.

Formulas are compiled once per statement into a `formula_op` program
(operator enum `op_code`, operand types `operand_type`) and evaluated
per element by the interpreter `formula_eval()` inside OpenMP loops.
Generated temporaries use the reserved prefixes `gen_sum`, `gen_par`,
`gen_pow`, `gen_mul`, `gen_add` in preprocessed TAB text.

## 3. Execution phases

`main()` orchestrates (all phases rank-parallel unless noted):

1. **Options and MPI topology.** PETSc init; per-node communicators
   (`node_comm`, `node_tail_comm`) derived from processor names.
2. **CMF read** (`cmf_read`), TAB preprocess (`tab_preprocess` writes the
   per-rank `_temp_tab_file<rank>.tab` working copy).
3. **Sets** (`sets_read`, `sets_read_intertemporal`, `subsets_read`,
   `subset_map_build`) — element lists plus superset position maps.
4. **Declarations** (`coefficients_read`, `variables_read`) — offsets and
   strides into the value vector.
5. **Data** (`data_read_files`) then **formulas** (`formulas_execute`) —
   fills `elem_vals[ncofele + nvarele]` (`store_real` = float storage).
6. **Closure and shocks** (`closure_read`, `shocks_read`) — the
   exogenous/endogenous partition must leave `VecSize − nexo` equations
   square; shocks divided across subintervals when `-nsubints > 1`.
7. **Equation ordering** (`equation_order_read[_nested]`) — assigns
   equations/variables to diagonal blocks by region/time per the matrix
   method; computes `netcut` (border size) [HK16 §5]. Rank 0 records the
   ordering to `<solfiles>.stats.json` (`ordering_stats_write`, §1).
8. **Jacobian preallocation and fill** (`jacobian_preallocate`,
   `jacobian_fill`) — PETSc AIJ matrices A (endogenous) and B
   (exogenous columns → RHS).
9. **Solve** (§5/§6), **updates**, iterate per step/subinterval.
10. **Outputs** — CSVs (`outputs_write_csv`) and solution binaries.

## 4. Source guide

Layout after the 2026 restructuring (`src/NAMING.md` maps every old name
to its new one, with the literature source of each term):

| file | lines | role |
|---|---|---|
| `main.c` | ~1,440 | orchestration: options, MPI topology, phases, solve dispatch, solution output |
| `teems_solver.h` | ~330 | umbrella header: constants, globals, typedefs, enums, structs, prototypes grouped per module |
| `globals.c` | 13 | single definitions of the program-wide globals |
| `str_util.c` | 109 | case-insensitive string search helpers |
| `cmf_io.c` | ~1,060 | CMF reading, TAB preprocessing, data-file headers, output CSVs |
| `tab_parse.c` | ~3,650 | TAB language: sets, declarations, data reads, closure, shocks |
| `formula.c` | ~3,280 | FORMULA compile/eval, UPDATE application, subinterval re-shocking |
| `jacobian.c` | ~2,510 | equation ordering, derivative-matrix preallocation and fill |
| `block_order.c` | 541 | (N)DBBD row/column ordering into bordered block form |
| `block_solve.c` | ~3,240 | (N)DBBD parallel factorization, interface problem, back-solve |
| `solve_drivers.c` | ~2,040 | Johansen and modified-midpoint drivers; spill/residency logic |
| `hsl_kernels.f90` | ~1,450 | Fortran wrappers around HSL MP48/MA48/MA51/MC66 plus sparse kernels |
| `hsl_kernels.h` | — | C extern declarations for the Fortran kernels |
| `makefile`, `mp48_mod.sh`, `patches/` | — | build; HSL source patching (see below) |

### main.c

A single `main()`: option parsing (§11); the node-level communicator
split (`node_comm`, `node_tail_comm`); CMF read and TAB preprocess;
phases 3–8 of §3 in order; the `-inmemory` default and residency
estimate; solve dispatch on `solution_method` (`solve_johansen`,
`solve_modified_midpoint`, and the `SM_STOCHASTIC`/`SM_STOSIM` repeat
loops); finally the output CSVs and the five solution binaries (§1).
Johansen LU solves run inline here through `spec48_single_` /
`spec48_nomc66_`.

### teems_solver.h / globals.c

Everything shared: buffer-size constants (`TABREADLINE` …); the globals
(`scratch_dir`, `inmemory`, thread counts, Gragg step state, node
communicators — defined once in `globals.c`); precision typedefs
(`solve_real` double / `store_real` float); index typedefs (`dim_t`,
`offset_t`); the enums of §5–§6 plus `bound_type`, `op_code`,
`operand_type`; and the model data structures — `set_def`/`set_element`,
`array_def` (offset + strides layout), `elem_value`
(value/initial/substep_base), `sum_def`/`sum_value`, `closure_entry`,
`quantifier`, `eq_var_ref`, `formula_op`. Prototypes are grouped under
per-module section banners.

### str_util.c

Stateless helpers: `str_rfind_any`, `str_rfind_ci`, `str_count_char`,
`str_count_ci`, `str_rfind_toplevel` (last occurrence of a character
with balanced parentheses after it — used to find a formula's top-level
`=`), `str_replace_char_all`. The replace family
(`str_replace_all[_bounded]`, `str_replace_first[_bounded]`,
`str_replace_char`, `str_find_ci`, `str_strip_comment`) lives at the
bottom of `tab_parse.c` for historical reasons.

### cmf_io.c

`cmf_count_files` / `cmf_read` parse the CMF (§1). `tab_preprocess`
writes the per-rank working TAB (`_temp_tab_file<rank>.tab`): statement
normalization, duplicate SUM-index renaming (`sum_dedup_indices`), and
variable-declaration append (`tab_write_variables`).
`datafile_read_header_info` / `datafile_read_labels` read text data-file
headers. `outputs_write_csv` writes the post-simulation
set/coefficient CSVs. `csv_read_ints` reads the optional `-nestfile`.

### tab_parse.c

The TAB front end, driven per statement type over the preprocessed file:

| area | functions |
|---|---|
| statement iteration | `tab_next_statement`, `tab_next_statement_resolved` (also inlines scalar coefficient values and applies `zerodivide`), `closure_next_statement`, `tab_count_statements` |
| sets | `sets_count`, `sets_read`, `sets_read_intertemporal`, `set_union_named` / `set_union_op` / `set_difference`, `subsets_read`, `subset_map_build` (fills `set_element.superset_pos`), `set_find_alltime` |
| declarations | `coefficients_read`, `variables_read`, `variables_read_defaults`, `coef_resolve_sets` |
| data | `data_read_files` (READ statements → element values) |
| closure/shocks | `closure_read` (`.cls`), `shocks_read` (`.shf`; divides shocks across subintervals) |
| formula surface | `formula_normalize`, `sum_parse` / `sum_count` / `sum_extract`, `eq_replace_linvar` / `eq_zero_linvar` (substitute or zero a linear variable inside an equation) |

### formula.c

The FORMULA/UPDATE engine. `formulas_execute` iterates formula
statements: parses quantifiers and conditions, resolves the target
array, compiles the RHS once via `formula_compile` (helpers:
`formula_compile_pow` / `_muldiv` / `_addsub` / `_if` precedence passes;
`formula_bind_operand` resolves each operand's array, strides, and
superset maps; `formula_subst_scalar` inlines scalar coefficients;
`leadlag_encode` / `parse_index_leadlag` handle `{t+1}` offsets), then
evaluates the resulting `formula_op` program per element with
`formula_eval` under OpenMP (per-thread program copies). `sum_eval`
materializes `SUM(...)` partials into `sum_value` arrays.
`updates_apply` (midpoint flag per §5) and `updates_apply_product`
implement UPDATE statements. `subinterval_update` re-applies shocks at
subinterval boundaries, optionally spline-interpolated (`cubic_spline`
from `solve_drivers.c`) — the `Isbiupd` path used by stochastic modes.

### jacobian.c

Linearized-system assembly. `equation_order_read` (flat) and
`equation_order_read_nested` (NDBBD) scan equations and variables,
assign them to diagonal blocks by region/time, and produce the
orderings plus `netcut` (border size) [HK16 §5]. `jacobian_preallocate`
counts nonzeros per row for exact PETSc preallocation (diagonal and
off-diagonal, for both A and B). `jacobian_fill` walks each equation
block, differentiates the linearized terms per element, and fills
A (endogenous columns) and B (exogenous columns, whose product with the
shock vector forms the RHS). `eq_sum_parse` / `eq_sum_eval` /
`eq_sum_replace` / `eq_linvar_read` handle SUM terms and
linear-variable references inside equations.

### block_order.c

`dbbd_order` (static/regional) and `ndbbd_order` / `ndbbd_order_presolve`
(intertemporal, nested) permute rows and columns into bordered block
form: per-block orders in `row_order`/`col_order`, block sizes in
`block_sizes`, block rank/singularity detection via `spec51m_rank_`
(MA51). With `-presol`, orderings persist in `_rank/_row/_col` scratch
files for reuse across runs.

### block_solve.c

The [KH19] steps 1–5 machinery:

- `dbbd_solve` — DBBD: per block, extract diagonal/border submatrices,
  factorize and form B·V (`spec48m_msol_`; factors handed off via the
  `_vav`/`_irnv`/`_keep` scratch files), merge interface contributions
  (`my_spar_add3l_` / `_add4l_` / `_compl_`, `reduce_to_rank*`), solve
  the interface problem (`spec48_ssol2la_` / `spec48m_ssol2la_`), then
  back-solve every block (`spec48m_esol_` / `spec48m_rpesol_`).
- `ndbbd_presolve` / `ndbbd_solve` — the nested variant: inner
  (regional) solves per time block with factors kept in caller arrays
  (`prep48_alu1_` / `prep48m_msol_`), then the outer interface.
- `ndbbd_block_solve` (file-staged) and `ndbbd_block_solve_mem`
  (array-staged) — back-solve of one nested block group.
- `reduce_to_rank` / `reduce_to_rank_nocompress` — chunked MPI
  reductions of dense interface vectors to a target rank.

### solve_drivers.c

`solve_johansen` and `solve_modified_midpoint` (§5): shock-vector
assembly, per-step Jacobian refill and solve (LU inline, or
`dbbd_solve` / `ndbbd_*` per `-matsol`), update application, Richardson
extrapolation and per-element precision accounting, the subinterval
loop. All class-1 spill/reload sites of §8 are here, guarded by
`inmemory`. `cubic_spline` supports `subinterval_update`.

### hsl_kernels.f90 / hsl_kernels.h

Fortran wrappers and kernels; factor-file paths take the
`TEEMS_SCRATCH` prefix (§8):

| subroutine | role | called from |
|---|---|---|
| `SPEC48_SINGLE` | SBBD solve via HSL_MP48 with MC66 ordering | main, solve_drivers |
| `SPEC48_NOMC66` | SBBD solve with the direct (precomputed) ordering [HK16] | main, solve_drivers |
| `SPEC51M_RANK` | MA51 rank/singularity detection per block | block_order, block_solve |
| `SPEC48M_MSOL` | factorize a block + multi-RHS solve (forms B·V); writes the factor files | block_solve |
| `SPEC48M_ESOL`, `SPEC48M_RPESOL` | back-solve from saved factors (plain / with refinement workspace) | block_solve |
| `SPEC48_SSOL2LA`, `SPEC48M_SSOL2LA` | single solve with workspace-doubling retry (interface problem) | drivers/main; block_solve |
| `PREP48_ALU1`, `PREP48M_MSOL` | factorize keeping factors in caller arrays (in-memory inner blocks) | block_solve |
| `MY_SPAR_ADD3L`, `MY_SPAR_ADD4L`, `MY_SPAR_COMPL` | merge/compress sorted sparse triplets during B·V accumulation (64-bit indices) | block_solve |
| `SPAR_MULMIN`, `SPAR_MULNOADD`, `SPAR_VBIVIADD` | sparse mat-vec kernels (y−=Ax; y=Ax; add B·V into a vector) | block_solve |
| `PATIO_MAT` | index packing during interface assembly | block_solve |

### Build and packaging files

`makefile` — PETSc-based; `BUILD_DIR` and `OPT ?= -Ofast` overridable;
serial `make` required for the HSL module dependencies. `mp48_mod.sh` +
`patches/*.patch` — TEEMS modifications applied to the staged HSL MP48
sources (64-bit duplicate-detection work arrays; `MA48→ZA48`-style
symbol renames so the statically patched copies cannot clash with the
dynamically linked libma48/libma51). `BUILDING.md` — build walkthrough.
`NAMING.md` — the old→new terminology map with literature citations.
Container builds live in `docker/`: `expedited_build` (reuses the
prebuilt `teems_base` image) and `full_build` (builds MPICH and PETSc
from source); both copy `./src` and stage the proprietary HSL packages
from tarball build-args.

## 5. Solution methods (`-solmed`, `enum solution_method`)

| value | method |
|---|---|
| `Johansen` | one-step solution of the linearized system [D20; Johansen 1960] |
| `Mmid` | Gragg's modified midpoint with Richardson extrapolation over `-step1/-step2/-step3` step counts (default 2-4-8), per subinterval [GM "Gragg"; Pearson 1991] |
| `Stochastic`, `StoSim` | Mmid-based stochastic variants (repeat solves, `-stoiter`) |
| `NoSol` | preparation only |

Mmid mechanics: for each step count `s`, the shock is applied in `s`
sub-steps; after each sub-step the Jacobian is refilled and solved, data
updated (midpoint rule: value advances from `substep_base` by twice the
computed change). The three solutions are Richardson-extrapolated with
weights `extrap_w1..3` derived from the step ratios; per-element error
codes are accumulated (`xc124`) and reported as precision counts.
teems-R maps `solution_method="mod_midpoint"` → `Mmid`; `"Johansen"`
forces `-nsubints 1`.

## 6. Matrix methods (`-matsol`, `enum matrix_method`)

| value | method | domain | notes |
|---|---|---|---|
| 0 `MM_LU` | serial MA48 LU | any | single rank does the factorization |
| 1 `MM_SBBD` | singly bordered block diagonal via HSL_MP48 [HK16] | **intertemporal** | blocks by time period; MC66 ordering optional (`-withmc66`), direct ordering by default; factors held in MP48 memory (sized by `-laA/-laDi`) |
| 2 `MM_DBBD` | doubly bordered block diagonal [KH19] | **static** (regional blocks) | C-side block factorization (MA48 kernels), factor handoff via scratch files |
| 3 `MM_NDBBD` | nested DBBD [KH19] | **intertemporal, T ≫ R** | regional blocks nested inside time blocks; requires `-regset`; `-nesteddbbd 1` |

Production guidance (measured, see §9): SBBD is the primary intertemporal
method; DBBD the primary static method; NDBBD targets many-timestep,
few-region models — its crossover point vs SBBD is still being mapped.

The bordered solve follows [KH19] steps 1–5: LU-factor diagonal blocks in
parallel; form interface contributions **B·V** (`vecbivi` — the paper's
B_i·V_i product); reduce to the interface problem (size = `netcut`);
solve it; back-solve all blocks. `reduce_to_rank()` performs the chunked
MPI reductions. MA51 (`spec51m_rank_`) detects rank/singularity per
block; singularity indicators appear in the log and are checked by
teems-R's `chk_solver_log`.

Workspace sizing: `laA`, `laDi`, `laD` scale factor/workspace sizes as
`ceil(la/100 × nnz)`; too small → factorization failure (increase and
re-run). teems-R defaults 300/500/200.

## 7. Parallelism

- **MPI** ranks via PETSc. With LU/SBBD (`nohsl == false`) only
  `rank_hsl = 0` holds the assembled system for HSL; with DBBD/NDBBD
  (`nohsl == true`) every rank owns blocks (`rank_hsl = rank`).
- Diagonal blocks are distributed round-robin over ranks; rank counts
  should divide the block count reasonably (blocks = time periods or
  regions × periods).
- **OpenMP** within rank: `-maxthreads` (global), `-smllthreads` (caps
  selected sections); formula evaluation and updates parallelize over
  elements with per-thread copies of the op-list.
- NDBBD's last (interface) block rank can be overridden via
  `-ndbbd_bl_rank` + `-nestfile <csv>`.

## 8. Scratch, memory, and the in-memory mode

All solver temp files go under `scratch_dir` — default `/tmp/`, override
with `-tempdir <dir>` or `TMPDIR`; exported to the Fortran kernels as
`TEEMS_SCRATCH`. User run directories receive only real outputs plus the
small per-rank preprocessed TAB copies.

Two spill classes exist:

1. **Driver array spills** (`solve_drivers.c`): before each factorization
   the value arrays (`elem_vals`, `closure_vals`) and Mmid step state
   (`clag1`, `varchange`, `xc*`) are written to `_temp*.bin` files and
   freed, then re-read after the solve.
2. **Block factor handoff** (DBBD/NDBBD): per-block LU factors
   (`VA/IRN/KEEP`). Since the 5.9 redesign the Fortran kernels
   (`SPEC48M_MSOL`, `PREP48M_MSOL`, `PREP48_ALU1`) return factors purely
   through caller-owned buffers (`KEEP` bound `M+9N+7`, live length in
   `INSIZE[12]`); the **C side** then either holds them resident
   (`-inmemory 1`: `fac_*` arrays in `dbbd_solve`, the `ndbbd_fac_*`
   store spanning the `ndbbd_presolve`→`ndbbd_solve` pair) or writes the
   legacy `_vav/_irnv/_keep` scratch files (`-inmemory 0`; same names
   and bytes). The C-written NDBBD interface products
   (`_bivi/_rbvi/_cbvi`, `_rank/_row/_col`) still spill in both modes —
   measured ~65MB per 202k-eq solve vs the ~825MB of factor traffic
   the resident mode eliminates; making them resident is a follow-up
   (their `_rank` rewrite-under-one-name state machine needs care).

`-inmemory 1` disables class 1 entirely (arrays stay resident), keeps
class-2 factors resident, and relocates scratch to tmpfs (`/dev/shm`)
unless `-tempdir`/`TMPDIR` is pinned. **Default: on for LU, SBBD, and
DBBD; off for NDBBD** — the NDBBD default predates the 5.9 factor
handoff (its rationale — factor-file traffic needing the page cache —
no longer applies in resident mode) and awaits an idle-machine A/B at
bench-inter-L scale before flipping. A startup check estimates
residency cost and falls back to spilling (with a warning) if it
exceeds half of `MemAvailable`. Explicit `-inmemory 0/1` always wins.

## 9. Measured characteristics

Benchmarks (GTAP12 database; 20-core/30GB node; 2 ranks unless noted):

| benchmark | size | config | wall | peak RSS/rank |
|---|---|---|---|---|
| static R32×medium | 1.35M eq | LU-1, Johansen | 12.4s | 0.41–0.48GB |
| static R32×medium | 1.35M eq | DBBD-2, Johansen | 20.4s | 0.61–0.66GB |
| intertemporal WB7×medium, T=21 | 4.36M eq | SBBD-2, Mmid | 107–115s | 1.1–1.6GB |
| intertemporal WB7×medium, T=21 | 4.36M eq | NDBBD-2, Mmid | 314–334s | 1.1–1.3GB |
| intertemporal WB7×medium, T=21 | 4.36M eq | DBBD-2, Mmid (off-domain) | ~420s | 2.1–2.6GB |

Scratch traffic before the in-memory work (202k-eq model, per solve):
SBBD 274MB, NDBBD 1.36GB — the motivation for §8.

T/R crossover sweep (WB7 and big3 × macro sectors, T ∈ {6,11,21,41},
min-of-2, production inmemory defaults): **SBBD won every shape** up to
T/R = 13.7; NDBBD-4 approaches SBBD-2 but never passes it, and LU trails
SBBD ~2× at larger T. At dev-box scale (≤800k equations, ≤4 ranks) no
NDBBD-over-SBBD region exists; the question escalates to larger models
and rank counts (see §13). SBBD also wins **single-rank** intertemporal
runs (202k eq: SBBD-1 6.3s vs LU-1 8.1s, min-of-3 interleaved).

Static LU/DBBD crossover (Johansen, min-of-reps, interleaved; 2026-07):

| rig | size | regions | LU-1 | DBBD-2 | DBBD-4 | DBBD-8 |
|---|---|---|---|---|---|---|
| R32×medium | 1.35M eq | 33 | **12.8s** | 20.4s | 14.8s | 19.1s |
| full×macro | 1.62M eq | 163 | 18.3s | 17.6s | **11.2s** | 14.4s |
| R32×full | 3.37M eq | 33 | 65.7s | — | **47.6s** | 59.1s |
| full×food | 8.32M eq | 163 | 450s | 357s | **284s** | 311s |

DBBD-4 is the sweet spot on this node (DBBD-8 always slower — border
serialization + contention). Crossover: DBBD overtakes LU by ~2M
equations at 33 regional blocks, and already by ~1.6M at 163 blocks.
These calibrate teems-R's `matrix_method = "auto"` rule (SBBD if
intertemporal; else DBBD when n_tasks ≥ 2 and the deployed system ≥ 2M
equations, or ≥ 1.5M with ≥ 100 regions; else LU — the 1.5–2M × <100
region band is uncalibrated and deliberately falls to LU).

Determinism: repeated same-binary runs are bit-identical for every
method (basis of the golden-run verification, below).

## 10. Numerical characteristics

- Solve precision `solve_real` = double; storage `store_real` = float
  (halves memory traffic for coefficient values; solution vectors and
  factorizations are double).
- Compiled `-Ofast` (includes `-ffast-math`: FP reassociation, no
  NaN/Inf guarantees). This is the tested production configuration;
  `make OPT=-O3` builds an IEEE-conformant binary for cross-checks.
- `zerodivide` defaults substitute configured values on 0/0 per [GM].
- Levels bounds (`bound_type`) are enforced during formula evaluation.

## 11. Command-line reference

| option | default | meaning |
|---|---|---|
| `-cmdfile <path>` | `./reg.cmf` | CMF file |
| `-matsol {0,1,2,3}` | 0 | matrix method (§6) |
| `-solmed <name>` | `Mmid` | solution method (§5) |
| `-step1/-step2/-step3` | 2/4/8 | Gragg step counts (all odd or all even) |
| `-nsubints n` | 1 | shock subintervals |
| `-regset <SET>` | — | regional set for ordering (required: LU output ordering, DBBD, NDBBD) |
| `-enable_time` | off | time-based SBBD ordering override |
| `-nesteddbbd {0,1}` | 0 | nested (NDBBD) ordering |
| `-presol {0,1}` | — | persist/reuse ordering & factor preparation |
| `-laA/-laDi/-laD n` | 2 (teems-R: 300/500/200) | workspace sizing, % of nnz |
| `-cntl_3 x` | — | HSL iterative-refinement threshold |
| `-cntl_6 x` | 0 | MA50 ordering control |
| `-withmc66 {0,1}` | 0 | MC66 ordering for SBBD |
| `-maxthreads n` | 1 | OpenMP threads |
| `-smllthreads n` | maxthreads | thread cap in selected sections |
| `-tempdir <dir>` | `/tmp/` (or `TMPDIR`) | scratch directory |
| `-inmemory {0,1}` | 1 except NDBBD | §8 |
| `-nsbbdblocks n` | 2 | SBBD block-count hint |
| `-ndbbd_bl_rank n` / `-nestfile <csv>` | — | NDBBD interface-block rank override |
| `-stoiter n` | 1 | stochastic repetitions |
| `-nowrites n` | 0 | suppress output writes |
| `-nox` | — | PETSc: no X output |

teems-R populates these from `ems_solve()` arguments
(`solution_method`, `matrix_method`, `n_tasks`, `steps`,
`n_subintervals`, `laA/laD/laDi`, `append_args` for anything else).

## 12. Verification and development infrastructure

- **Golden runs** (`.audit/verify.sh`, dev-machine): rebuilds from the
  working tree in the `teems-audit` container and checks nine solves
  bit-identically against manifests anchored to the pre-refactor binary:
  GTAPv7 static LU/Johansen; GTAP-RE intertemporal LU/Mmid; GTAP-RE
  SBBD, DBBD, NDBBD at 2 ranks; three real-shock runs (2D probe, 4D,
  swap); and a GEMPACK-orientation matrix shock (2D, both dims free).
  Any behavior change fails the gate.
- **Benchmark rigs**: `.audit/bench_run.sh` (wall/RSS via `time -v`),
  `.audit/strace_run.sh` (per-file write-byte accounting), deployments
  under `.audit/bench-*` produced by teems-R scripts.
- Build: see `BUILDING.md` (HSL staging, `src/patches/` applied by
  `mp48_mod.sh`, serial make requirement, `OPT` knob).

## 13. Known limitations and planned work

- **NDBBD without `-regset` aborts** with a bare MPI error (exit 59)
  instead of a diagnostic.
- **Error handling** is largely `printf` + `exit`/`return`; parse errors
  on one rank can abort non-collectively.
- Fixed-size line buffers (`TABREADLINE` = 20000) bound statement length;
  overflow warnings exist but are not hard guards everywhere.
- The binary is named `hsl` for historical reasons; renaming to
  `teems-solver` requires coordinated changes in teems-R and the images.
- Planned (measured motivation on file): formula op-list caching across
  steps (removes per-step TAB re-parsing); verbosity-gated logging
  (3,000+ unconditional printfs); NDBBD-vs-SBBD crossover mapping at
  scale. Done since this list was written (see ROADMAP 5.9/5.10):
  factor handoff through memory (`-inmemory` holds block factors
  resident); `formula_op` compaction (1416→968→848B, per-dim operand
  fields interleaved as `dim_addr` records).
