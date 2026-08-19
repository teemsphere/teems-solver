# TEEMS Solver — Technical Reference

The TEEMS solver (`teems-solver`; `teems-solver-f64` is the same
program with double-precision coefficient storage, §10) is a parallel
Computable General Equilibrium (CGE) solver in the GEMPACK tradition: it
reads a TAB-language model description, base data, a closure, and
shocks; linearizes the model; solves the resulting first-order system
with sparse direct methods; and writes updated data and
percentage-change solutions. It implements the bordered block-diagonal
ordering methods of Ha & Kompas (2016, SBBD) and Kompas & Ha (2019,
NDBBD) on top of PETSc, MPI, OpenMP, and the HSL sparse linear algebra
libraries (MA48, MA51, MC66, HSL_MP48; MC79 for the structural probe,
MA60/MC71 for the solve-quality diagnostics).

Authored by Pham Van Ha and Tom Kompas; restructured and hardened in the
2026 refactor (see `src/NAMING.md` for the symbol lineage and
`src/BUILDING.md` for the build). This document is the technical
reference for the solver; user-facing modelling rules (which TAB forms
the R front end accepts, rewrites or rejects) are in the teems manual
(`model_load`), and the solver's own build recipe is in the repository
`README.md`. Planning and design documents live in the separate
development repository, not here.

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
  /opt/teems-solver/solver/teems-solver -cmdfile /opt/teems/<model>.cmf <options>
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

- `out/sets/*.csv` — set listings; `out/coefficients/*.csv` /
  `out/postsim/*.csv` — updated coefficient values, one CSV per
  `Write` statement in the TAB (opt-in from teems-R since the
  coefficient dump below became the default transport; `%f`
  formatting, six fixed decimals).
- Solution binaries (prefix from `soldata`, default `solution`):

| file | contents |
|---|---|
| `.bin` | `nvarele × double` — percentage-change (or change) solution per variable element |
| `.var` | `nvar × array_def` — variable declarations (name, dims, set ids, strides) |
| `.set` | `nset × set_def` — set definitions |
| `.sel` | `nsetspace × set_element` — set elements with superset positions |
| `.mds` | 4 × long: `nsetspace, nvar, nvarele, nset` |
| `.cof` | coefficient dump header + declarations (`-cofdump`, default on): 4 × long `{version=1, ncof, ncofele, reserved}`, then `ncof × array_def` (same struct as `.var`), then `ncof × uint8 kind` (bit 0 = PostSim coefficient, bit 1 = `(parameter)`). Written after PostSim, so PostSim coefficients carry their computed values |
| `.cbin` | `ncofele × double` — updated (post-simulation) coefficient values, the coefficient slice of the value vector in `array_def.offset` order; the coefficient twin of `.bin`. Together with `.cof` this replaces the per-coefficient CSVs as teems-R's coefficient transport (2026-08) |
| `.stats.json` | per-run ordering statistics (v2): system size, method, `netcut`, border sizes, per-block variable/equation counts (null/empty when no bordered ordering was built), plus `chain_set`/`partition_set` with `chain_source`/`partition_source` (`explicit`/`structural`/`none`) and, when the partition probe ran, the full `partition_auto` candidate table (§6). Written before the solve, so failed runs still record their ordering; feeds `matrix_method` auto-calibration |

The structs are written raw; teems-R's `parse_solution.cpp` mirrors their
layout, so `array_def`/`set_def`/`set_element` are an ABI shared with
teems-R and must not be re-ordered or resized casually.

## 2. Model language

The solver reads a GEMPACK-style TAB subset (statement syntax per [GM]):

- `set` / `subset` — explicit lists, `read elements` from a data
  file, the full set-expression grammar of [GM] 10.1.1 (`+`, `-`, `\`,
  `UNION`, `INTERSECT`, quoted single elements, parentheses; `+`
  disjointness and `-` presence enforced), set equality, conditional
  set builders `Set X = (all,i,SRC: <cond>)` ([GM] 10.1.2 — evaluated
  from the input files ahead of set resolution by
  `tab_setbuilder_transform`; the condition coefficient must be
  file-Read), and intertemporal set declarations (`(intertemporal)`),
  which mark the time dimension used by the bordered orderings.
  Set products (`x`) are not supported.
- `mapping` ([GM] 11.9) — declared mappings between sets, values from
  `Read (by_elements)`, used in index expressions and in conditional
  sums (`sum{j,S: MAP(j) = i, ...}`); mapped equations solve under
  every matrix method.
- `coefficient` / `variable` — levels or percentage-change quantities
  (`(change)`, `(levels)` variables — levels equations and
  `Formula & Equation` pairs are linearized by `levels.c` at
  preprocess; `(default=…)` handled by `variables_read_defaults()`),
  `(parameter)`/`(non_parameter)`, `(integer)`; optional bounds
  (`GE/GT/LE/LT`, `enum bound_type`, two slots) with initial/updated
  range tests (`-range_test_*`).
- `read` — from the CMF-bound data files, including `(IfHeaderExists)`
  (absent header skipped) and `(by_elements)` mapping reads.
- `formula` — `(initial)` / `(always)` semantics as in [GM]; conditional
  quantifiers `(all,i,S: COEF(i) op c)` and general conditional sums;
  `zerodivide` defaults honored by `tab_next_statement_resolved()`
  (both GEMPACK zerodivide classes with `-gpzerodivide`); the [GM] 11.5
  intrinsics (ABS/MAX/MIN/SQRT/EXP/LOGE/LOG10/ID01/ID0V/ROUND/TRUNC0/
  TRUNCB, …; `$POS` is not supported). `IF[...]` conditional
  expressions are rewritten by the R front end into these forms.
- `assertion` (`-assertions` off/warn/fatal), `zerodivide` statements,
  `Default` statements ([GM] 10.19), `PostSim (Begin/End)` sections
  ([GM] ch. 12; `-postsim`, split by the preprocess and executed after
  the simulation), `write` (opt-in coefficient CSVs).
- `complementarity` ([GM] ch. 51) — approximate and accurate runs
  (`-comp_*` controls), under every matrix method.
- `equation` — linearized equations; the left/right sides are compiled
  per equation block by `jacobian_fill()`.
- `update` — product-form and explicit updates (`updates_apply()`,
  `updates_apply_product()`); the modified-midpoint variant applies
  `substep_base + 2·Δ` (see §5).
- `backsolve <var> using <eq>` — condensation ([GM] 10.16, 14.1.3): the
  named endogenous variable and its nominated defining equation are
  excluded from the solved system (`backsolve_read()` marks the elements;
  an "equation"-scan filter in `tab_parse.c` hides the equation from
  counting, ordering, preallocation, fill, and structural detection), and
  the variable's per-step changes are recovered from the retained
  equation after each solve, before the data updates
  (`backsolve_recover()`, riding the compiled statement programs of §9;
  zero pivot aborts the run).  The symbolic side of condensation —
  rewriting every other equation to eliminate the variable, `omit`,
  `substitute` — is teems-R's job (`ems_model(omit=, backsolve=)`);
  `omit`/`substitute` statements reaching the solver abort with a
  pointer back to model preparation.
- `swap` (closure) and shock statements — parsed by `closure_read()` /
  `shocks_read()` into per-element `closure_entry` records.

Shock statement values follow GEMPACK-standard ordering (first
subscript varies fastest; fixed dimensions collapsed to size 1).
Startup notes remind that intertemporal variables should be declared
with minimal dimensionality to keep the border (netcut) small;
element-level classification (§3 phase 7) borders only the elements a
lead/lag reference actually spans, so sliced references cost their
slice and full-set references cost the led periods rather than the
variable's whole element block.

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
   per-rank `_temp_tab_file<rank>.tab` working copy — statement
   normalisation, several statements per physical line split, block
   comments), then the preprocess passes over the working copy: the
   PostSim split, levels linearization (`levels.c`), mapping and
   set-builder transforms (`tab_setbuilder_transform` reads the
   condition coefficient straight from the data files, in their
   slice-per-blank-line text layout).
3. **Sets** (`sets_read`, `sets_read_intertemporal`, `subsets_read`,
   `subset_map_build`) — element lists plus superset position maps.
4. **Declarations** (`coefficients_read`, `variables_read`) — offsets and
   strides into the value vector.
5. **Data** (`data_read_files`) then **formulas** (`formulas_execute`) —
   fills `elem_vals[ncofele + nvarele]` (`store_real` = float storage).
6. **Closure and shocks** (`closure_read`, `shocks_read`) — the
   exogenous/endogenous partition must leave `VecSize − nexo` equations
   square; shocks divided across subintervals when `-nsubints > 1`.
   **Backsolve statements** (`backsolve_read`, `backsolve_validate_refs`)
   then mark condensed-out elements; `VecSize` shrinks by `nbselems` and
   the nominated defining equations vanish from every "equation" scan.
7. **Equation ordering** (`equation_order_read[_nested]`) — assigns
   equations/variables to diagonal blocks by region/time per the matrix
   method; computes `netcut` (border size) [HK16 §5]. References that
   are not block-safe on the chain dimension — a lead/lag offset, or a
   chain-dim index other than the equation's chain quantifier (sum
   indices, fixed elements) — border at element level
   (`border_mark_ref`): the mark spans exactly the referenced elements,
   each dim narrowed to its index set's positions (fixed dims arrive as
   tab_preprocess singleton subsets) shifted by that dim's offset and
   clipped. So `qo("capital",r,t+1)` borders only the capital slice at
   the led periods, and boundary-period elements of ordinary lead
   variables stay inside their diagonal blocks. A mark covering the
   whole variable collapses to the per-variable flag (`var_inter`),
   which the partition rules also read. Rank 0 records the
   ordering to `<solfiles>.stats.json` (`ordering_stats_write`, §1).
8. **Jacobian preallocation and fill** (`jacobian_preallocate`,
   `jacobian_fill`) — PETSc AIJ matrices A (endogenous) and B
   (exogenous columns → RHS).
9. **Solve** (§5/§6); **backsolve recovery** (`backsolve_recover`, per
   step, after the solve and before the updates — [GM] 14.1.3);
   **updates**; iterate per step/subinterval.
10. **Outputs** — CSVs (`outputs_write_csv`) and solution binaries.

## 4. Source guide

Layout after the 2026 restructuring (`src/NAMING.md` maps every old name
to its new one, with the literature source of each term):

| file | lines (2026-08) | role |
|---|---|---|
| `main.c` | ~2,640 | orchestration: options, MPI topology, phases, solve dispatch, solution output, `stats.json` |
| `teems_solver.h` | ~690 | umbrella header: constants, globals, typedefs, enums, structs, prototypes grouped per module |
| `globals.c` | ~130 | single definitions of the program-wide globals |
| `str_util.c` | 109 | case-insensitive string search helpers |
| `cmf_io.c` | ~2,180 | CMF reading, TAB preprocessing (incl. the PostSim split and set-builder transform), data-file headers, output CSVs, coefficient dump |
| `tab_parse.c` | ~5,760 | TAB language: sets and set expressions, mappings, declarations and qualifiers, data reads, closure, shocks, backsolve statements + equation-scan filter, name validation |
| `formula.c` | ~3,240 | FORMULA compile/eval, conditional quantifiers/sums, intrinsics, UPDATE application, subinterval re-shocking |
| `levels.c` | ~2,010 | levels variables / `Formula & Equation` linearization at preprocess |
| `jacobian.c` | ~3,060 | equation ordering, derivative-matrix preallocation and fill, backsolve recovery programs, complementarity state |
| `block_order.c` | 600 | (N)DBBD row/column ordering into bordered block form |
| `block_solve.c` | ~3,900 | (N)DBBD parallel factorization, interface problem, back-solve, persistent factors |
| `solve_drivers.c` | ~2,350 | Johansen, Gragg and Euler drivers; spill/residency logic; condest |
| `solve_rk.c` | ~960 | Runge–Kutta drivers (RK2/RK4, embedded BoSha32/DoPri54, adaptive control) |
| `probe.c` | ~820 | `-solmed probe`: MC79 structural diagnosis and `probe.json` |
| `hsl_kernels.f90` | ~2,530 | Fortran wrappers around HSL MP48/MA48/MA51/MC66/MC79/MA60 plus sparse kernels |
| `hsl_kernels.h` | — | C extern declarations for the Fortran kernels |
| `makefile`, `mp48_mod.sh`, `patches/` | — | build (`teems-solver` and `teems-solver-f64` targets); HSL source patching (see below) |

### main.c

A single `main()`: option parsing (§11); the node-level communicator
split (`node_comm`, `node_tail_comm`); CMF read and TAB preprocess;
phases 3–8 of §3 in order; the `-inmemory` default and residency
estimate; solve dispatch on `solution_method` (`solve_johansen`,
`solve_gragg`); finally the output CSVs and the five solution
binaries (§1).
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
set/coefficient CSVs.

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
implement UPDATE statements.

### jacobian.c

Linearized-system assembly. `equation_order_read` (flat) and
`equation_order_read_nested` (NDBBD) scan equations and variables,
assign them to diagonal blocks by region/time, and produce the
orderings plus `netcut` (border size) [HK16 §5]. `jacobian_preallocate`
counts nonzeros per row for exact PETSc preallocation (diagonal and
off-diagonal, for both A and B). A's diagonal block follows this rank's
row range; B's follows the exogenous column range the shock vector
carries, which parts company with the rows once condensation leaves
`nexo > VecSize` and B is wider than it is tall (`shock_vec_set_sizes` /
`shock_mat_set_sizes` set both layouts; NDBBD's time-block split applies
only while B stays square). The scan counts insertions rather than
distinct entries, so each count is capped at its block's row length
before PETSc sees it — a heavily substituted row can otherwise nominate
more entries than the block has columns, which PETSc rejects outright.
`jacobian_fill` walks each equation
block, differentiates the linearized terms per element, and fills
A (endogenous columns) and B (exogenous columns, whose product with the
shock vector forms the RHS). Since the 6.2 phase-0 refactor each
statement is compiled once per solve into a cached program
(`stmt_prog`: quantifier frames, per-occurrence coefficient ops from
`formula_compile`, SUM bodies as `sum_prog` records) and later fills
only re-evaluate it; the cache is per rank, keyed to the matrix
ownership range, and released by `jacobian_cache_free()`.
`eq_sum_parse` / `eq_sum_replace` / `eq_linvar_read` handle SUM terms
and linear-variable references during the build phase.

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

`solve_johansen` and `solve_gragg` (§5): shock-vector
assembly, per-step Jacobian refill and solve (LU inline, or
`dbbd_solve` / `ndbbd_*` per `-matsol`), update application, Richardson
extrapolation and per-element precision accounting, the subinterval
loop. All class-1 spill/reload sites of §8 are here, guarded by
`inmemory`.

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
| `Gragg` | Gragg's method (smoothed modified midpoint) with Richardson extrapolation over `-step1/-step2/-step3` step counts (default 2-4-8), per subinterval [GM "Gragg"; Pearson 1991 eq. 6.1/Alg. 7.1.2]; `Mmid` accepted as a deprecated alias (warns) |
| `Euler` | forward Euler multistep with Richardson extrapolation over the same three step counts [GM "Euler"] — shares the Gragg driver with the leapfrog and terminal smoothing disabled; the truncation error series is `h` (not `h²`), so the extrapolation weights use the step ratios unsquared and each extra solution gains one order (not two). Any strictly increasing step counts are allowed (no parity rule) |
| `RK2`, `RK4` | fixed-step explicit Runge–Kutta over `-step1` steps per subinterval (no extrapolation triple) |
| `BoSha32`, `DoPri54` | embedded Runge–Kutta pairs (Bogacki–Shampine 3(2), Dormand–Prince 5(4)); with `-adaptive yes` the embedded error estimate drives step control against `-epstol` (`accuracy-only` skips the check-failure retries; `-retryadj`/`-maxretries` set the retry policy) — the accuracy summary reports the per-step error metric |
| `probe` | preparation only — runs the full pre-solve pipeline (data, formulas, structural detection of both the chain dimension and the block partition irrespective of `-matsol`, ordering, `stats.json`) and skips the solve; a sub-second structure probe of a model (§6); `NoSol` accepted as a deprecated alias (warns). On a single rank the probe then assembles the condensed Jacobian and runs the HSL_MC79 maximum-matching / Dulmage-Mendelsohn structural diagnosis on (a) the full stored pattern — structural closure validity — and (b) the numerically realized pattern (entries nonzero at base data) — the zero-flow singularity class. Defects are reported as *named* variable and equation elements (under-determined variables = unmatched columns, over-constrained equations = unmatched rows, plus the coarse-DM entangled blocks), to the log and to `<solfiles>.probe.json` (version 2). The report also carries the statement-level equation-system structure — per statement its quantifier sets, row count and referenced variables with element-level incidence weights — and aggregates every defect/entangled set by variable and by equation statement, so the diagnosis stays readable at any scale. `-probefine 1` adds the fine decomposition's strongly-connected-component report: core-size histogram plus the composition (by equation and variable) of the largest simultaneous cores — the model's irreducible simultaneous structure vs its recursive remainder. Catches structural and zero-data singularity classes, not numerical near-singularity |

Gragg mechanics: for each step count `s`, the shock is applied in `s`
sub-steps — Euler first step, then midpoint leapfrog (value advances
from `substep_base` by twice the computed change), each preceded by a
Jacobian refill and solve — followed by one extra refill+solve at the
final state feeding the terminal smoothing pass
½(z_s + z_{s−1} + Δz), i.e. s+1 factorizations per pass (Pearson 1991
eq. 6.1). The three solutions are Richardson-extrapolated with weights
`extrap_w1..3` derived from the step ratios (valid for the even-power
error expansion Gragg's theorem guarantees for even `s`); per-element
error codes are accumulated (`xc124`) and reported as precision counts.

Euler mechanics: same driver and per-step refill+solve, but every
sub-step is a forward step from the current state (`updates_apply`
midpoint flag 0 throughout, no leapfrog history, no terminal
refill+solve), i.e. `s` factorizations per pass. Robustness properties
follow from using only local information: tolerates near-asymptote
levels, severe shocks, and exact −100 % percent-change shocks that the
midpoint family cannot start from (GM §"Gragg's method and the midpoint
method"). Step counts must be strictly increasing; parity is
unconstrained. Measured on GTAP-RE (40 % pfactwld, float64-coefficient
diagnostic build): extrapolated error contracts ~8× per step-count
doubling (the h³ rate), with Euler 8-16-32 matching Gragg 2-4-8.
teems-R maps `solution_method="Gragg"/"Euler"` → the same `-solmed`
values; `"Johansen"` forces `-nsubints 1`.

## 6. Matrix methods (`-matsol`, `enum matrix_method`)

| value | method | domain | notes |
|---|---|---|---|
| 0 `MM_LU` | serial MA48 LU | any | single rank does the factorization |
| 1 `MM_SBBD` | singly bordered block diagonal via HSL_MP48 [HK16] | **intertemporal** | blocks by time period; MC66 ordering optional (`-withmc66`), direct ordering by default; factors held in MP48 memory (sized by `-laA/-laDi`) |
| 2 `MM_DBBD` | doubly bordered block diagonal [KH19] | **static** (regional blocks) | C-side block factorization (MA48 kernels), factor handoff via scratch files |
| 3 `MM_NDBBD` | nested DBBD [KH19] | **intertemporal, T ≫ R** | partition blocks nested inside chain blocks; nested ordering auto-enabled |

Production guidance (measured, see §9): SBBD is the primary intertemporal
method; DBBD the primary static method; NDBBD targets many-timestep,
few-region models — its crossover point vs SBBD is still being mapped.

### Structural partition detection

The bordered orderings need two dimensions, and neither is a property
of a set's *name* — both are derived from the structure of the equation
system at run time:

- **Chain dimension** (conventionally time): the top-level set whose
  elements the equations couple through lead/lag index offsets
  (`x(...,t+1)`). Detected by scanning the equations for offset
  references; the `(intertemporal)` TAB qualifier only licenses the
  syntax and is cross-checked (a qualified set that no equation
  offset-references gets no chain treatment, offsets on an unqualified
  set draw a warning, offsets spanning several top-level sets draw a
  warning and the most-referenced set wins).
- **Diagonal-block partition** (conventionally regions; used by DBBD
  and NDBBD): every set with two or more elements outside the chain
  family that indexes at least one variable dimension — directly or
  through a subset, so finer subregional-style sets compete on equal
  terms — is probed: the pre-Jacobian ordering scan runs against it and
  the first counting pass measures the border and block sizes that
  partition would produce. A candidate is viable with ≥ `n_tasks`
  nonempty blocks and a border under half the system; the smallest
  border wins, near-ties (2 %) broken by block balance. In nested
  (NDBBD) mode the per-chain-block interface column counts as border —
  it is the local border the interface solve pays for. Sets coupled by
  cross-element sums (market clearing over regions, adjacency-matrix
  sums over spatial grids) price themselves out through the measured
  border rather than through any name- or qualifier-based rule.

The probe costs one equation-section scan per candidate, runs once
before the ordering, and is skipped entirely when the transitional
flags below resolve the dimensions explicitly. The full candidate table
(set, blocks, border, min/max block, viable) is logged at verbosity 1
and recorded in `stats.json` under `partition_auto`, with
`chain_source`/`partition_source` saying how each dimension was
resolved — the structural evidence teems-R's `matrix_method = "auto"`
consults for its `matrix_method = "auto"` decision. A `-solmed probe` run yields it without
solving, and independently of the `-matsol` it is launched with: the
probe always runs both detections, measures the nested (NDBBD)
partition geometry when a chain exists and the flat (DBBD) one
otherwise, and never aborts on a method-vs-structure mismatch (there is
no method to mismatch); the MC79 diagnosis it then runs is invariant to
the ordering, so the report is identical under any `-matsol`. Under a
real solve the detections are gated on the requested method (SBBD:
chain; DBBD/NDBBD: chain and partition). The scan sees
*declared* structure only; value-aware element-level classification
(e.g. recognising a banded adjacency coefficient as nearly
partitionable) is a planned refinement.

The transitional `-enable_time`/`-regset` overrides were removed once
structural detection reproduced their results bit-identically on the
golden suite; stale flags in old commands are silently ignored by the
options parser. Method-vs-structure
mismatches abort with a remedy: NDBBD when either dimension resolves
to none, DBBD when both do, and SBBD when no chain dimension exists
(previously HSL_MP48 received zero blocks, errored, and the run still
finished with exit 0 and no solution).

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

Persistent refactorize (`-fastrefac 1`; sequential LU, SBBD and DBBD): the
Jacobian's stored sparsity pattern is fixed across steps, so the
analyse work runs once per solve and every later step/stage
refactorizes on the kept pivot sequences. Sequential LU keeps the MA48
analyse (pivot-order search) and refactorizes with `MA48B/BD JOB=2`
(`ICNTL(11)=1` refactorizes any block whose entries turn unsuitable; a
declined fast factorize falls back to a fresh analyse). SBBD keeps the
whole MP48 instance alive across steps — border-column lists, per-block
pivot sequences and factors — and repeat steps refill `VALUES`/`B` on
the host and refactorize with `FACT_JOB=2` (numerical pivoting
retained; an MP48 error on the fast path rebuilds the instance). In
both cases the extraction keeps explicitly stored zeros — an entry that
is zero at the analyse state may become nonzero at a later step and
must stay in the pattern. Every MA48 factorization site is self-sizing: an MA48
`-3` (workspace too small) return grows the arrays to the size MA48
suggests (doubling floor, six-attempt cap) instead of aborting, the
grown size persists for the rest of the solve, and the sizes actually
used are recorded in `stats.json` (`la_used`) so a caller can warm-start
the next run. Measured (Gragg 2-4-8,
17 factorization cycles): 1.35M static LU −72% wall / +21% RSS; 4.36M
intertemporal SBBD-2 −36% wall / +28% RSS; 202k SBBD-2 −25%. DBBD keeps
per-diagonal-block persistent factors: each block's COO is a raw copy
of its stored CSR (pattern structurally stable), so the per-block
analyse runs once and later steps refactorize with `MA48B/BD JOB=2`;
the persistent arrays double as the within-step factor handoff to the
back-solve (no scratch files under the flag). NDBBD keeps its regional (inner) blocks' factors persistent the same
way, extending the `ndbbd_fac` store with the MA48 column mapping; the
flag forces the resident factor store even in NDBBD's default disk
mode. The interface problems are factorized fresh every step in all
bordered methods — their patterns are assembled from value-dependent
products and may legitimately change between steps. Off by default;
results shift only within factorization rounding (the analyse-state
pivot sequences applied at other step states), bounded by the usual
cross-method noise floors — on the DBBD and NDBBD acceptance rigs the
flag-on solutions came out exactly identical.

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
- NDBBD's per-time interface blocks are rank-revealed in-solve: the
  ordering presolve bounds each block's rank by min(nrow,ncol), the
  MA51 presolve factorization measures the true numerical rank, and
  `ndbbd_order` permutes any rank-deficient rows/columns out to the
  outer border before the solve proper.

## 8. Scratch, memory, and the in-memory mode

All solver temp files go under `scratch_dir` — default `/tmp/`, override
with `-tempdir <dir>` or `TMPDIR`; exported to the Fortran kernels as
`TEEMS_SCRATCH`. User run directories receive only real outputs plus the
small per-rank preprocessed TAB copies.

Two spill classes exist:

1. **Driver array spills** (`solve_drivers.c`): before each factorization
   the value arrays (`elem_vals`, `closure_vals`) and Gragg step state
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
| intertemporal WB7×medium, T=21 | 4.36M eq | SBBD-2, Gragg | 107–115s | 1.1–1.6GB |
| intertemporal WB7×medium, T=21 | 4.36M eq | NDBBD-2, Gragg | 314–334s | 1.1–1.3GB |
| intertemporal WB7×medium, T=21 | 4.36M eq | DBBD-2, Gragg (off-domain) | ~420s | 2.1–2.6GB |

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

Phase-5 closing A/Bs (2026-07-09, bench-inter-L, 3 interleaved pairs
each; raw walls in `.audit/ab_phase5_results.txt`):
- **`formula_op` compaction (5.10)** — speed-neutral (median delta 0.4%
  inside ~10% session drift; a 10% `formula_eval` win would be ~1.6% of
  wall, unresolvable by wall-clock A/B). Value is the 968→848B/op
  footprint and layout clarity.
- **NDBBD `-inmemory`** — wall neutral (medians 385s off vs 381s on);
  RSS 1.09→1.84GB/rank (+69%) for resident factors. Default stays off
  on RAM-cost grounds; on slow/network scratch the resident mode
  (`ems_solve(inmemory = TRUE)`) is expected to win — this box's fast
  local scratch + page cache is the best case for the file path.
- **Small-model cost structure** (GTAP-RE 10.5k eq, LU-1 Gragg, 1.16s
  wall): `jacobian_fill` = 0.71s (61%, equation parse + mandatory
  per-step refill); `formulas_execute` = 0.04s (3.4%). Formula op-list
  caching was dropped on these numbers; the equation path belongs to
  condensation.
- **Compiled-equation cache** (2026-07-16, condensation phase 0,
  `.audit/ab_eqcache.sh`, interleaved pairs, first run discarded):
  `jacobian_fill` split into a build phase (parse + `formula_compile`
  programs, run once per solve per rank) and an execute phase (SUM +
  element evaluation, run per step); later fills re-walk the statement
  file only to refresh `zerodivide` defaults. Wall medians old → new:
  GTAP-RE 10.5k LU-1 Gragg **1.60s → 0.76s (−52%)** — the small-model
  lever above, collected; 202k LU-1 Gragg 10.90s → 10.25s (−6%); 202k
  SBBD-2 Gragg 5.22s → 4.94s (−5%, new faster in every pair). RSS cost
  is the cache itself: +3–8MB/rank. Golden-bit-identical (14/14) and
  valgrind-clean; the remaining per-fill parse work is the statement
  walk, which stays because `zerodivide` may reference scalar
  coefficients whose values change between steps.

- **ISA level and BLAS threading** (2026-07-14, 6.6(d) A/Bs through
  the runtime images, `.audit/ab_march.sh`): `x86-64-v3` shows **no
  win** over `v2` — SBBD-2 on bench-inter-L medians 128.0s (v3) vs
  125.5s (v2), static LU-1 neutral — consistent with the BLAS3
  kernels already dispatching per-CPU at runtime and the rest being
  memory-bound; only the v2 base is published. The `x86-64-v2` flag
  scheme itself proved **bit-identical** to the historical
  no-`-march` build across the entire golden suite (v2 adds no
  FP-semantic instructions; FMA arrives at v3). Separately,
  `OPENBLAS_NUM_THREADS=1` beat the uncapped pthread pool by **~20%**
  on SBBD-2 at 4.36M eq (medians 103.5s vs 128.7s, decisive in every
  interleaved pair — one BLAS pool per rank oversubscribes the box);
  it is baked into the runtime images as an `ENV` and is
  golden-bit-identical.

Determinism: repeated same-binary runs are bit-identical for every
method (basis of the golden-run verification, below).

## 10. Numerical characteristics

- Solve precision `solve_real` = double; storage `store_real` = float
  in `teems-solver` (halves memory traffic for coefficient values;
  solution vectors and factorizations are double). `teems-solver-f64`
  is the same source built with double coefficient storage for
  accuracy-critical runs (teems-R `precision` argument); both binaries
  ship in every image.
- Compiled `-Ofast` (includes `-ffast-math`: FP reassociation, no
  NaN/Inf guarantees). This is the tested production configuration;
  `make OPT=-O3` builds an IEEE-conformant binary for cross-checks.
- `zerodivide` defaults substitute configured values on 0/0 per [GM].
- Levels bounds (`bound_type`) are enforced during formula evaluation.

Solve-quality diagnostics (`-condest 1`, sequential LU): each linear
solve is measured with HSL MA60 — iterative refinement plus the
Arioli–Demmel–Duff componentwise backward error ω₁/ω₂, a forward-error
bound, and the corresponding scaled condition numbers κω₁/κω₂ (Hager
1-norm estimation via MC71 internally). MA60 drives the estimate by
reverse communication, requesting `A⁻¹y` and `A⁻ᵀy` products that the
kept MA48 factors supply through `MA48CD` with `TRANS` set per `KASE`;
this transpose-solve requirement is why the flag is LU-only — the
bordered methods never factorize the composed system and have no
transpose path through their block+interface flows.

The measurement is *diagnostic-only*: refinement runs on a copy of the
solution and the caller's `X` is never updated, so a run's outputs are
bit-identical with the flag on or off (gated in both `-fastrefac`
states). Cost is a handful of back-solves per step — 0.18 s against a
2.7 s factorize at 1.35M — plus one pristine copy of the staged system
(the `-fastrefac` path already keeps one for its refill).

Reading the numbers: ω is the backward error (how far the solved
problem sits from the stated one); a backward-stable solve of a
near-singular system legitimately reports ω ≈ 0 *and* a near-zero
forward-error bound, so κω₂ — not `ERX` — is the conditioning signal.
Raw 1-norm condition numbers are useless here (CGE Jacobians mix
$-million and percent rows, so healthy systems already estimate ~1e18);
the scaled κω₂ has real dynamic range, and κω₂ > 1e15 raises a
near-singularity warning. That is precisely the class the MC79
structural probe cannot see: a system whose entries are all present and
nonzero, and structurally full-rank on both probe patterns, but
numerically degenerate at the base data (§5). Null-shock solves are
skipped with a note — with b = 0 there is nothing to measure. No
f32/f64 advisory is derived from κω₂: healthy models already exceed a
naive κ·ε test while validating to 1e-8 in f32, so a useful threshold
needs corpus calibration.

## 11. Command-line reference

| option | default | meaning |
|---|---|---|
| `-cmdfile <path>` | `./reg.cmf` | CMF file |
| `-matsol {0,1,2,3}` | 0 | matrix method (§6) |
| `-solmed <name>` | `Gragg` | solution method (§5) |
| `-step1/-step2/-step3` | 2/4/8 | Gragg step counts (all odd or all even) |
| `-nsubints n` | 1 | shock subintervals |
| `-laA/-laDi/-laD n` | 2 (teems-R: 300/500/200) | workspace sizing, % of nnz |
| `-fastrefac {0,1}` | 0 | all matrix methods: analyse once, fast refactorize per step (MA48 JOB=2 / MP48 FACT_JOB=2); LU auto-grows `laA` (§6) |
| `-cntl_3 x` | — | HSL iterative-refinement threshold |
| `-cntl_6 x` | 0 | MA50 ordering control |
| `-withmc66 {0,1}` | 0 | MC66 ordering for SBBD |
| `-maxthreads n` | 1 | OpenMP threads |
| `-smllthreads n` | maxthreads | thread cap in selected sections |
| `-tempdir <dir>` | `/tmp/` (or `TMPDIR`) | scratch directory |
| `-inmemory {0,1}` | 1 except NDBBD | §8 |
| `-nsbbdblocks n` | 2 | SBBD block-count hint |
| `-probefine {0,1}` | 0 | with `-solmed probe`: add the MC79 fine-DM strongly-connected-component report (§5) |
| `-condest {0,1}` | 0 | sequential LU only: per-solve quality diagnostics (MA60/MC71) — componentwise backward error ω₁/ω₂ with iterative refinement, forward-error bound, and the Arioli–Demmel–Duff scaled condition numbers κω₁/κω₂ — logged per linear solve and recorded as run maxima in `stats.json` (`condest` object). Diagnostic-only: solutions are bit-identical with the flag on or off (refinement runs on a copy). Null-shock (zero-rhs) solves are skipped with a note; κω₂ > 1e15 adds a numerically-near-singular warning — the class the structural probe cannot see. Ignored (with a warning) for the bordered methods, whose composed systems have no transpose-solve path |
| `-adaptive {no,yes,accuracy-only}` | no | embedded RK pairs only: error-controlled stepping against `-epstol` (default 0.1); `-retryadj` (0.5) / `-maxretries` (3) tune the check-failure retries |
| `-assertions {0,1,2}` | 1 | TAB `Assertion` statements: off / warn / fatal |
| `-range_test_initial`, `-range_test_updated {0,1,2}` | 2 / 1 | coefficient bound checks on initial and updated values: off / warn / fatal |
| `-gpzerodivide {0,1}` | 0 | GEMPACK's two-class zerodivide semantics (nonzero-by-zero vs zero-by-zero defaults) |
| `-postsim {0,1}` | 1 | execute the TAB's PostSim section after the simulation |
| `-comp_steps n`, `-comp_do_approx/-comp_do_acc {0,1}`, `-comp_redo {0,1}`, `-comp_redo_min_frac x`, `-comp_sberr_warn {0,1}` | see log | complementarity ([GM] ch. 51) approximate/accurate run controls |
| `-nowrites n` | 0 | suppress output writes |
| `-cofdump {0,1}` | 1 | write the `.cof`/`.cbin` coefficient dump (recorded in `stats.json` `options`) |
| `-verbosity {0,1,2}` | 1 | 0 = errors/warnings + accuracy summary only; 1 = phase progress and timings; 2 = per-rank/per-block debug detail (also exported as `TEEMS_VERBOSITY` for the Fortran kernels; MA48 duplicate-entry notes appear only at 2) |
| `-nox` | — | PETSc: no X output |

teems-R populates these from `ems_solve()`/`ems_RK()` arguments
(`solution_method`, `matrix_method`, `n_tasks`, `steps`,
`n_subintervals`, and the validated expert dots — `laA/laD/laDi`,
`fastrefac`, `condest`, `assertions`, `range_test_*`, RK controls,
complementarity controls). Options travel on the invocation, never in
the CMF (which is a file manifest only), and are echoed in
`stats.json`'s `options` object. Stale flags from older commands
(`-regset`, `-enable_time`, `-presol`, `-nesteddbbd`) are ignored.

## 12. Verification and development infrastructure

- **Golden runs** (`.audit/verify.sh`, dev-machine): rebuilds from the
  working tree in the `teems-audit` container (an `x86-64-v2` base
  since 2026-07-14; build mirrors the expedited image's
  `make OPT="-Ofast $ARCH_FLAGS"`) and checks 14 solves bit-identically
  against manifests anchored 2026-07-14 to the v2 flag scheme with
  flag-free commands (the original anchors traced back to the
  pre-refactor binary; the v2 re-anchor was proven to carry exactly
  the flag-removal delta — the v2 build itself changed no output bit):
  GTAPv7 static LU/Johansen; GTAP-RE intertemporal LU/Gragg; GTAP-RE
  SBBD, DBBD, NDBBD at 2 ranks; three real-shock runs (2D probe, 4D,
  swap); a GEMPACK-orientation matrix shock (2D, both dims free); and
  four subinterval gates spanning method × ranks × nsubints × inmemory.
  Any behavior change fails the gate.
- **Structural-detection acceptance** (`.audit/accept_structural.sh`):
  reruns the NDBBD, DBBD (intertemporal + static) and SBBD golden
  configurations flag-free and requires the structurally detected
  ordering to reproduce the golden manifests bit-identically (this
  gate is what licensed removing the transitional
  `-regset`/`-enable_time` overrides).
- **Benchmark rigs**: `.audit/bench_run.sh` (wall/RSS via `time -v`),
  `.audit/strace_run.sh` (per-file write-byte accounting), deployments
  under `.audit/bench-*` produced by teems-R scripts.
- **Feature kits** (`.audit/*-test-kit/run_*_tests.sh`): per-feature
  acceptance rigs on copies of the golden runs — set builders, set
  expressions, conditional sums, mapping, complementarity, levels,
  PostSim, qualifiers/defaults/bounds/names, intrinsics, zerodivide,
  condest, closure, methods (analytic and polynomial-exactness anchors
  for the stepping methods), probe. Each pins values with in-TAB
  `Assertion`s or against an independent computation of the same
  quantity.
- **Memory safety**: ASan/TSan builds of the golden suite
  (`.audit/verify-asan.sh`, `verify-tsan.sh`) and AFL++ fuzzing of the
  TAB/CMF/closure/shock readers (batches landed 2026-07); the
  `-fastrefac` and grow-and-retry paths were added under the same
  gates.
- Build: see `src/BUILDING.md` (HSL staging, `src/patches/` applied by
  `mp48_mod.sh`, serial make requirement, `OPT` knob) and the README's
  expedited image recipe. `.audit/` is developer tooling on the
  maintainer's machine, not part of the distributed tree.

## 13. Known limitations and planned work

- **Error handling** is largely `printf` + `exit`/`return`; parse errors
  on one rank can abort non-collectively. Messages carry
  `Error:`/`Warning:` prefixes with file/flag/remedy named, and teems-R
  maps them to typed conditions; the exit status is checked as well as
  the log.
- Fixed-size line buffers (`TABREADLINE` = 20000) bound statement length;
  overflow is diagnosed, not silently truncated.
- **Language forms rejected by design** (named errors): `$POS`, set
  products (`x`), formula-assigned mappings (a mapping needs a
  `Read (by_elements)`), formula-computed operands in conditional set
  builders, LOOP/DISPLAY/TRANSFER/BREAK/CYCLE, `(no_split)`,
  `linear_name=`, positional `Default` semantics through the R
  pipeline. Compound `AND/OR/NOT` conditions and `IF` nested inside
  expressions are not evaluated by the solver; the R front end
  rewrites the supported `IF` forms into conditional quantifiers,
  helper coefficients and domain splits before deployment.
- **Version handshake**: solver and teems-R are versioned together
  (v2.0.0 series; the `hsl` compatibility symlink is dropped at 2.0.0).
  A `-version` flag / banner / `stats.json` field surfacing one version
  constant, with a teems-R pre-flight, is planned for the 2.0.0
  release; until then the image label is the only version hint.
- **Calibration**: the `matrix_method`/`n_tasks` auto rules, the
  `-fastrefac` and MA48 pivot-threshold default flips, the NDBBD
  `-inmemory` default, RK-vs-Gragg rankings, and the structural-probe
  thresholds were fit at dev-box scale (≤ 4 ranks, ≤ 8M equations);
  a many-core calibration program on instance-scale rigs is the next
  step and will re-anchor those defaults.
