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

- names — coefficients, variables, sets, elements, indices, equations
  and logical files may contain `@` ([GM] 11.2.1, e.g. `X3@Y2`). The
  preprocess lowers a mapping call `MAP(i)` to the flat token `MAP~i`
  (`MAPMARK`, a character no GEMPACK name can contain), and the
  complementarity machinery's own names are recognised by their exact
  forms (`name_is_comp_derived`: `C@e/@d/@l/@u`, the `C@w*` weights,
  `e_C@…` equations for a declared complementarity `C`, and
  `del_comp@`), so a user `@` name is an ordinary name in closures,
  shocks, formulas and index positions. Text inside `# labels #` is
  opaque: a quoted element there is not lowered. Header names
  (`Read ... header "Nat"` against a data header `NAT`, set-element
  and `(by_elements)` mapping reads alike) and logical file names are
  case-insensitive, as in GEMPACK.
- `set` / `subset` — explicit lists, `read elements` from a data
  file, the full set-expression grammar of [GM] 10.1.1 (`+`, `-`, `\`,
  `UNION`, `INTERSECT`, quoted single elements, parentheses; `+`
  disjointness and `-` presence enforced), set equality, conditional
  set builders `Set X = (all,i,SRC: <cond>)` ([GM] 10.1.2 — evaluated
  from the input files ahead of set resolution by
  `tab_setbuilder_transform`; the source set may be declared by list,
  read, or derived — `A + B`, `A union B`, `A - B`, `A intersect B` —
  and may be a subset of the condition coefficient's dimension set;
  the condition coefficient must be
  file-Read or an indicator assigned constants only, its narrowing
  quantifier over a declared subset or the front end's synthesized
  `"ele" & RANGE` set, which the pre-pass resolves itself), and
  intertemporal set declarations (`(intertemporal)`),
  which mark the time dimension used by the bordered orderings, and
  set products `Set P = A x B` ([GM] 10.1.6; elements `a_b`, first
  factor fastest, the 11.7.11 compromise names when a pair would exceed
  12 characters, duplicates fatal). A set may be empty ([GM] 11.7.9) —
  an explicit `()`, a builder that selects nothing, a difference or
  intersection with nothing left: statements quantified over it have no
  tuples, sums over it are zero and variables over it have no elements.
- `mapping` ([GM] 11.9) — declared mappings between sets, values from
  `Read (by_elements)`, from `Mapping (project)` onto a factor of a
  set product ([GM] 10.13.2), or from Formulas ([GM] 10.13.1: a
  codomain position, typically `$POS(...)`, or a quoted codomain
  element under `(by_elements)`; a formula-assigned mapping must be
  complete before its first use and is frozen after it, 11.9.9); used
  in index expressions and in conditional sums
  (`sum{j,S: MAP(j) = i, ...}`, where `i` ranges over the codomain or
  a declared subset of it, or names a codomain element); the argument position may be declared over
  the codomain or a superset of it (`VST(COM2MAR(c),r)` with `MAR`
  a copy of `MARG_COMM`, `qo(i,S2R(s))` with `RGA ⊂ REG`), and the
  mapped index may range over the domain or a declared subset of it
  ([GM] 11.9.7) — the position routes through `superset_pos` on
  either side; unrelated sets stay a named fatal. Mappings may appear
  in formulas, equations, assertions and on the right-hand side of
  updates (`NCTAXLEV(r) = del_nctaxb(REGTOBLOC(r))`); a mapping in a
  left-hand-side argument list is a named fatal ([GM] 11.9.9); mapped
  equations solve under every matrix method, and a backsolved
  variable's defining equation may reference other variables through
  mappings (the recovery program routes the column like the Jacobian
  fill), while the backsolved variable's own occurrence must be plain
  (a mapped pivot is many-to-one: named fatal).
- `coefficient` / `variable` — levels or percentage-change quantities
  (`(change)`, `(levels)` variables — levels equations (with
  conditional sums — a coefficient condition becomes an `if()` factor
  on each linearized term, a mapping condition stays on the sum — and
  SQRT/EXP/LOGE/LOG10, differentiated by the chain rule), equations
  whose left side holds a conditioned sum (split at the top-level `=`),
  and `Formula & Equation` pairs are linearized by `levels.c` at
  preprocess; `(default=…)` handled by `variables_read_defaults()`),
  `(parameter)`/`(non_parameter)`, `(integer)`; optional bounds
  (`GE/GT/LE/LT`, `enum bound_type`, two slots) with initial/updated
  range tests (`-range_test_*`).
- `read` — from the CMF-bound data files, including `(IfHeaderExists)`
  (absent header skipped) and `(by_elements)` mapping reads.
- quoted element arguments, `V(c,"dom",i)` — lowered by the preprocess
  to a singleton subset and an index bound to it. The owning name is
  the identifier immediately before the enclosing `(`, whatever
  precedes it (`*[p(c,"dom")…]`, `{p(c,"dom")…}`, a name glued to its
  quantifier `(all,c,MAR)SALES(c,"dom")`); an equation argument that
  is neither a quantifier index, a sum index nor a lowered literal is a
  named fatal (it used to take its set from arbitrary statement text —
  a dropped or garbled term, or a SEGV). One synthesized set
  per (owning set, element) per statement, so a literal repeated in a
  statement shares one quantifier (E_regx0_A of ORANI-G carries ten
  `"dom"` arguments); an Equation or indexed Formula gets a leading
  `(all,iN,sub_iN)` quantifier, a scalar-target Formula has its
  right-hand side wrapped in `sum(iN,sub_iN, …)` instead. A qualifier
  group after the keyword (`equation (levels) NAME`, `formula
  (initial)`, `update (change)`, `formula & equation`) stays ahead of
  the synthesized quantifier, and an equation's name ahead of it too.
- repeated indices — `T(d,d)`, `TRADE(c,"dom",r,r)`, `v(e,e)` bind
  one frame index at several argument positions; each position adds
  its own stride (positions routed differently, e.g. a subset index on
  a coefficient declared over the subset and the full set, carry a
  second address), in formulas, equations (coefficients and variables,
  inside sums too) and updates. A repeated sum index is one loop, not
  one per position.
- file opens — every work, scratch and output file open is checked:
  a failure is a named fatal with the path, the access and the system
  reason (a leftover file owned by another uid in the run directory
  used to SEGV at startup); `stats.json` reads tolerate an absent file.
- names `C_*` — a coefficient may start with `c_` unless a variable of
  the tail name exists (then `c_NAME` is that variable's change
  column: named fatal); `c_` tokens that name a declared coefficient
  are not rewritten to `p_` in equations and updates.
- `formula` — `(initial)` / `(always)` semantics as in [GM]; a
  qualifier list `(initial, by_elements)` is split into its groups, and
  `write updated value to file F header "H"` is accepted with a warning
  (no separate updated-data file is written; the value is in the
  coefficient dump); conditional quantifiers
  `(all,i,S: <expr> op <expr>)` — the historical `COEF(i) op number`
  form directly, anything else (`$POS(c)=$POS(g)`, a coefficient or
  arithmetic on either side, several parenthesised groups) compiled as
  two expressions over the frame; a mapping on one side (`MAP(z) = r`,
  `MAP(z) = "ele"`, `A(z) = B(y)`) compares codomain positions, and a
  mapped argument inside a numeric comparison evaluates like any
  expression; a sum inside a quantifier condition is a named fatal;
  mapping Formulas also take an identity or copy right-hand side
  (`MAP(c) = c` with `c` over the codomain or a subset, `A(z) = B(z)`);
  a sum condition's coefficient may be indexed over subsets of its
  dimension sets (summed or quantifier index); an enclosing sum's
  condition no longer leaks into a nested sum's carried set; and
  general conditional sums,
  where the condition coefficient may also be a scalar (`(all,i,S: SC >
  1)`, `sum(i,S: SC > 1, …)`) and the word comparisons `EQ NE GT LT GE
  LE` are normalized to symbols while the condition is read; an
  undeclared condition name or a reference with the wrong argument
  count is a named fatal;
  `zerodivide` defaults honored by `tab_next_statement_resolved()`
  (both GEMPACK zerodivide classes with `-gpzerodivide`); the [GM] 11.5
  intrinsics (ABS/MAX/MIN/SQRT/EXP/LOGE/LOG10/ID01/ID0V/ROUND/TRUNC0/
  TRUNCB, …) and `$POS` in its five forms ([GM] 11.5.6: index,
  index in a superset, element literal in a set, mapped index, mapped
  index in a superset; compiled to a per-tuple position operand).
  Element-membership `IF[i in S, ...]` expressions are rewritten by
  the R front end into these forms; value-comparison `if(a op b, x)`
  compiles directly, with the comparison spelled either way
  (`< <= = <> > >=` or `lt le eq ne gt ge`).
- `assertion` (`-assertions` off/warn/fatal; conditional quantifiers
  `(all,i,S: C(i) > 0)` evaluate per tuple, each condition compiled as
  a residual through the formula engine — sums and `$POS` inside a
  condition still skip the assertion with a warning), `zerodivide` statements,
  `Default` statements ([GM] 10.19; `Equation (default=levels)` makes
  each following equation written without a qualifier a levels equation,
  applied by the levels transform until `Equation (default=linear)`),
  `PostSim (Begin/End)` sections
  ([GM] ch. 12; `-postsim`, split by the preprocess and executed after
  the simulation), `write` (opt-in coefficient CSVs).
- `complementarity` ([GM] ch. 51) — approximate and accurate runs
  (`-comp_*` controls), under every matrix method.
- `equation` — linearized equations; the left/right sides are compiled
  per equation block by `jacobian_fill()`.
- `update` — product-form, `(change)` and `(explicit)` updates
  (`updates_apply()`, `updates_apply_product()`); the modified-midpoint
  variant applies `substep_base + 2·Δ` (see §5). Every right-hand side
  reads the values at the start of the step and several Updates of one
  element apply in order, the later overriding ([GM] 11.12): a
  statement another statement could read, or which shares its target,
  holds its results until the pass ends (the rest — product updates of
  a coefficient nothing else touches — write in place). Final updated
  data after an extrapolated solve: product updates are the exact
  one-shot form from the extrapolated totals; a `(change)`/`(explicit)`
  target takes its per-pass path values (Gragg: terminal-smoothed like
  the variables; accumulated in double precision), Richardson-
  extrapolated with the variables' weights ([GM] 26.2). A step counter
  (`ITER = ITER + 1`) extrapolates meaninglessly, as in GEMPACK — run
  such models single-pass (`-single_run 1`). Quantifier conditions
  `(all,i,S: <expr> op <expr>)` on Updates are evaluated per element
  from start-of-step values; an element whose condition is false is
  left alone (no write, so it does not override an earlier Update).
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
   Every binding of a quantifier index to an argument position of a
   coefficient or variable resolves through `set_supset_slot`: the
   index's set must be the declared set or a declared/implied subset of
   it (manual 10.1.2), else `set_supset_fatal` names the index, the
   symbol and the missing `Subset` statement — in equations, formulas
   (operands, targets, IF conditions), partial Reads and set-qualified
   closure entries alike. The former fallback bound the index by its
   position within its own set and silently addressed the wrong
   elements of the declared set.
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
| `solve_rk.c` | ~1000 | Runge–Kutta drivers (RK2/Heun/RK4, embedded BoSha32/DoPri54, log-level chart, adaptive control with stage-level retries) |
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
communicators — defined once in `globals.c`); `version.h`
(`TEEMS_SOLVER_VERSION`, the one version constant); the two message
channels — `logmsg(level, …)` for progress gated by `-verbosity`, and
`errmsg(…)` for every `Error:` line (prints, flushes, counts; the
count is the exit-status backstop of §11); precision typedefs
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

Before an equation reaches the readers, every linear-variable
reference in it is prefixed `p_` (the internal linear-variable
marker, change and percent-change alike; levels variables are
excluded, their bare token is the value). The follower set that
recognises a token includes `)` and `}` since 2026-09-06: a scalar
variable closing a bracket group (GTAPv7's `[pxwcom(c) - pxwwld]`)
was never prefixed, so the equation reader bound it as a value and
dropped its column from `E_c1_cr`, `E_cnttotr` and `E_cntpinv`; the
linearity check in `jacobian.c` now names such a token as an unknown
name instead of letting it through.

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
superset maps, and aborts with a named error when a reference carries
a different number of indices than its declaration (a bare `ESUN` for
`ESUN(t)` used to hand a NULL index to the positional binder: SEGV in
an equation term, silent element-0 binding in a formula);
`formula_subst_scalar` inlines scalar coefficients;
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
and linear-variable references during the build phase. Before the
term splitter sees an equation, `eq_linearity_check` parses the
normalised `RHS-(LHS)` text by recursive descent and rejects, as a
named fatal citing manual 11.4.8, any use of a linear variable other
than coefficient-expression * variable: a variable inside a function
argument, as the base or exponent of `^`, multiplied by or dividing
another variable, in the conditional part of an IF or SUM, an unknown
function name, or an equation with no linear variable at all. The
splitter assumes that shape and silently mis-read anything else (a
power lost its exponent, a product bound one column and left MA48 a
structurally singular matrix). The nonlinear form belongs in
`Equation (Levels)` (see `levels.c`).

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

`makefile` — PETSc-based; `BUILD_DIR`, `OPT ?= -Ofast` and `WARN ?= -Wall` (our sources build warning-free; staged HSL sources are `-w`) overridable;
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
| `RK2`, `Heun`, `RK4` | fixed-step explicit Runge–Kutta over `-step1` steps (no extrapolation triple): midpoint, explicit trapezoid (strong-stability-preserving) and classic fourth order. Stages are combined in the log-level chart by default (`-rkchart log`: every endogenous percent-change element is carried as log(1+X/100), so no percent variable can reach −100 % at any stage and every tableau keeps its order — Munthe-Kaas on the multiplicative group; `-rkchart percent` is the GEMPACK-orientation arithmetic); exogenous elements stay on the exact level-linear path |
| `BoSha32`, `DoPri54` | embedded Runge–Kutta pairs (Bogacki–Shampine 3(2), Dormand–Prince 5(4)), first-same-as-last (an accepted step costs 3 / 6 new solves; a retried step reuses its stage-0 solve). With `-adaptive yes` the embedded estimate drives step control against `-epstol`, a step is redone at `-retryadj` when a stage state fails a coefficient range test, an assertion, the −100 % crossing (percent chart), the `-rkguard` level ratio or (LU) a singular factorization — the manual 26.5.1 triggers, abandoned at the failing stage; `accuracy-only` acts on the estimate alone. The accept test is steered by the percent-change elements (`-rkscope pct`; `all` restores GEMPACK's rule) in the `-rknorm max|rms` norm with the elementary or PI controller (`-rkctrl`). The accumulated estimate is written beside the solution as `sol.est` (one double per element, the paper's metric |Δ|/max(1,|X|)) — an indicator of the least-settled elements, measured NOT to be a bound; the run record (`runge_kutta` in stats.json) carries steps, stage solves, reuse and the rejection census |
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
The updated data follow the same scheme (§2 `update`): each pass's
`(change)`/`(explicit)` Update targets are smoothed like the variables
(`updates_apply` mode 2) and extrapolated with the same weights
(`updates_path_accumulate`); product updates are recomputed exactly
from the extrapolated totals at the end of each subinterval.

Single-pass runs (`-single_run 1`): the driver runs one pass (sol 0)
with weight 1 instead of the triple — Gragg keeps its terminal
smoothing (the GEMPACK single Gragg calculation), Euler is the plain
forward pass; subintervals chain passes as usual. An `ITER = ITER + 1`
counter then counts steps, and a RAS-style explicit update applied
per step gives exactly N iterations (pmfix kit `ras`).

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
- **NDBBD thread budget.** Four NDBBD regions own per-thread working
  sets that scale with the interface, not the rank count, so resident
  memory grows with threads × footprint (Q34, 234M equations, one rank:
  62 GB + 2.9 GB per thread measured at 4/8/16 threads). Each region is
  sized from its own structure just before its team starts and its team
  is capped so the ranks on the node fit 0.8 × the memory the node can
  still give: the presolve factor+schur team (dense interface product +
  Schur staging + MA48 workspace of the largest regional block), the
  interface-rank probe (the compressed interface staged at `-laDi`
  percent, 16 B per entry, + MA51 workspace), the interface
  factorization in `ndbbd_solve` (the same staging + a pristine copy for
  workspace-growth retries) and the outer Schur region (a chain block's
  factor set read back from scratch; aliased when resident). The cap
  only lowers `-maxthreads`/`-smllthreads`; block-to-thread assignment
  does not change results (goldens bit-identical at 1 and 4 threads).
  Each choice is logged once as `ndbbd threads: <region> rank r uses k of
  n (per thread X GB: <breakdown>; budget B GB = 0.8 x A GB available /
  m rank(s) on this node)` (again only when it moves; every step at
  `-verbosity 2`) and recorded in `stats.json` as `ndbbd_threads`
  (`asked`, `presolve`, `interface_rank`, `interface_factor`, `schur`).
  "Available" is `MemAvailable` capped by the container's cgroup limit
  minus its anonymous usage (`/sys/fs/cgroup/memory.max` or the v1
  `memory.limit_in_bytes`), so a Docker Desktop `--memory` cap is honoured
  even though `/proc/meminfo` reports the VM's total.
- **Reproducibility across thread counts.** LU, SBBD and NDBBD solutions
  are bit-identical run to run at any `-maxthreads` (measured 2026-09-14,
  golden GTAPv7 condensed rig, 2 ranks × 2 threads). DBBD is not: its
  border product and border right-hand side accumulate block
  contributions with OpenMP atomics (`SPEC48M_MSOL` VECBIVI, `dbbd_solve`
  vecbiui), so the floating-point summation order follows thread timing.
  A single Johansen step at 2 or 4 threads differs between two runs by at
  most ~1e-4 absolute on values of ~2e6 (7e-11 relative; 98 % of
  elements within 1e-9), numerical zeros flip sign, and Gragg 2-4-8 on a
  rough rig (2-digit accuracy on tens of thousands of elements) amplifies
  this to 1e-6..1e-3 relative on a sixth of the elements. Bit-for-bit
  reproducibility under DBBD needs `-maxthreads 1`; a fixed-order
  reduction (per-block buffers summed in block order) is the remedy if
  reproducibility at threads > 1 is ever required, at the cost of one
  border-sized buffer per block.
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
exceeds half of the available memory (`MemAvailable` capped by the
container's cgroup limit, see §7). Explicit `-inmemory 0/1` always wins.

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

- **ISA level re-checked** (2026-09-14, solver 1015811, `.audit/ab_march2.sh`,
  fresh v2/v3 base pair, real-shock rigs): `x86-64-v3` vs `v2` medians
  on the Linux box 32.9/33.0 s (static Johansen LU-1), 340.5/344.7 s
  (static Gragg 2-4-8 LU-1), 101.5/102.7 s (intertemporal Gragg SBBD-2,
  P-core pinned); on a Windows Docker Desktop laptop 76.0/76.6,
  632.0/637.3 and 383.3/414.0 s, v3 slower in every intertemporal
  pair. Verdict stands: no v3 win, v2 only. On hybrid CPUs pin
  multi-rank A/Bs to the performance cores (`AB_CPUSET`), otherwise
  rank placement adds ±20 % noise.

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
- **BLAS kernel pin.** Debian's `libopenblas` is built `DYNAMIC_ARCH`
  and selects its kernels from CPUID at run time, so an unpinned image
  produced different last digits per host: measured on the same image,
  a laptop dispatching `Haswell` against a workstation whose 13th-gen
  part the library does not recognise falling back to `Prescott`, with
  211 of 1,173 variables differing (forcing `Haswell` on both made all
  1,173 bit-identical). The runtime images therefore set
  `OPENBLAS_CORETYPE` to the family matching the ISA level they were
  compiled at — `NEHALEM` for the `x86-64-v2` base, `HASWELL` for
  `x86-64-v3`, `ARMV8` for `arm64` — recorded in
  `${BUILD_DIR}/blascore` by the base build, and the image build fails
  if the library does not actually select the pinned family (an
  unknown name is otherwise ignored silently). Every run records what
  was selected: `BLAS kernels:` in the log banner and `blas_core` in
  the `options` object of `<solfiles>.stats.json`. Speed is not a
  reason to pin or to unpin: the wall is dominated by the sparse
  factorization, and the same case measured 112.7/98.5 s under
  `Prescott` against 107.9/98.0 s under `Haswell`.
- **What the pin buys.** One-task runs of the same deployment are
  byte-identical on every host an image variant supports. At more than
  one task, cross-machine agreement is to rounding, not to the byte:
  the order in which ranks reduce is not fixed (2 of 562 variables at
  ~1e-14 in the measured pair). Different variants (v2, v3, arm64) are
  different kernel families and agree to the same rounding level, not
  bit for bit; and `-Ofast` reassociation means a build at another
  optimization level is a different numerical configuration again.

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
| `-version` (`--version`) | — | print exactly `teems-solver <version>` and exit 0, before MPI/PETSc start (no PETSc banner; the same constant as the log banner and `stats.json`). An image whose solver does not answer predates 1.1.0 |
| `-cmdfile <path>` | `./reg.cmf` | CMF file |
| `-matsol {0,1,2,3}` | 0 | matrix method (§6) |
| `-solmed <name>` | `Gragg` | solution method (§5) |
| `-step1/-step2/-step3` | 2/4/8 | Gragg step counts (all odd or all even) |
| `-single_run {0,1}` | 0 | Euler/Gragg: one pass of `-step1` steps, no Richardson extrapolation (GEMPACK `method = euler; steps = N;`); accuracy estimates are reported as unavailable, updated data are that pass's path values; ignored for Johansen, a named fatal with the Runge–Kutta methods; `stats.json` records `"single_run"` |
| `-nsubints n` | 1 | shock subintervals |
| `-laA/-laDi/-laD n` | 2 (teems-R: 300/500/200) | workspace sizing, % of nnz |
| `-fastrefac {0,1}` | 0 | all matrix methods: analyse once, fast refactorize per step (MA48 JOB=2 / MP48 FACT_JOB=2); LU auto-grows `laA` (§6) |
| `-ma48u x` | library defaults | MA48 / HSL_MP48 pivot threshold `CNTL(2)`, 0 < x ≤ 1, applied at every factorization site (sequential LU, DBBD/NDBBD blocks and rank probes, SBBD's MP48 instance). Absent = each library's own default (MA48 0.1, MP48 0.01), bit-identical to builds without the option; recorded in `stats.json` (`ma48u`, null when default). A calibration knob, not a tuning recommendation |
| `-cntl_3 x` | — | HSL iterative-refinement threshold |
| `-cntl_6 x` | 0 | MA50 ordering control |
| `-withmc66 {0,1}` | 0 | MC66 ordering for SBBD |
| `-maxthreads n` | 1 | OpenMP threads |
| `-smllthreads n` | maxthreads | thread cap in selected sections |
| `-tempdir <dir>` | `/tmp/` (or `TMPDIR`) | scratch directory |
| `-inmemory {0,1}` | 1 except NDBBD | §8 |
| `-ndcutcache {0,1,2,3}` | 1 | NDBBD only: reuse the step-1 cuts for the following steps and RK stages instead of re-probing every block with a throwaway MA48 factorization. The regional cut (rank-deficient block tails migrated into the time-interface blocks, permutations, block sizes) is structural and its reuse is bit-identical; the interface cut (each chain block's MA51 rank probe: rank plus row/column selection) is value-dependent, and reusing it moves outputs at the rounding level (≤1e-6 on the reference rigs, accuracy summaries identical). Dropped at every subinterval start and at the complementarity re-entry. 0 = re-probe every step (the previous behaviour; the remedy named in the MA48B failure message when a cached interface selection goes singular on a later step); 2/3 = regional-only / interface-only, bisect aids |
| `-nsbbdblocks n` | 2 | SBBD block-count hint |
| `-probefine {0,1}` | 0 | with `-solmed probe`: add the MC79 fine-DM strongly-connected-component report (§5) |
| `-probepattern {0,1}` | 0 | with `-solmed probe`: also write the assembled structural pattern with element names to `<solfiles>.probe.pattern` — one line per equation element, its variable elements space-separated (`=0` marks an entry that is zero at base data); the element-level view for tracing an entangled block by hand when the aggregated report is not enough |
| `-condest {0,1}` | 0 | sequential LU only: per-solve quality diagnostics (MA60/MC71) — componentwise backward error ω₁/ω₂ with iterative refinement, forward-error bound, and the Arioli–Demmel–Duff scaled condition numbers κω₁/κω₂ — logged per linear solve and recorded as run maxima in `stats.json` (`condest` object). Diagnostic-only: solutions are bit-identical with the flag on or off (refinement runs on a copy). Null-shock (zero-rhs) solves are skipped with a note; κω₂ > 1e15 adds a numerically-near-singular warning — the class the structural probe cannot see. Ignored (with a warning) for the bordered methods, whose composed systems have no transpose-solve path |
| `-adaptive {no,yes,accuracy-only}` | no | embedded RK pairs only: error-controlled stepping against `-epstol` (default 0.1; teems-R defaults to 0.01, which is where the log-chart driver matches the previous driver's accuracy at a third of the factorizations); `-retryadj` (0.5) / `-maxretries` (3) tune the check-failure retries (warn-level range/assertion violations are accepted with a warning once the retries are exhausted) |
| `-rkchart {log,percent}` | log | RK stage chart (§5) |
| `-rkscope {pct,all}` | pct | which elements steer the adaptive accept test (§5) |
| `-rknorm {max,rms}` | max | error-metric norm of the accept test |
| `-rkctrl {std,pi}` | std | step-size controller: elementary (0.85 safety, 0.5–2 clamp) or PI |
| `-rk_h0 h` | 0 | first step length; 0 = `1/step1` capped from the initial gradient (exp(5) level ratio in the log chart, 90 points in the percent chart) |
| `-rkguard r` | 1e30 | log chart: level ratio beyond which a stage state is rejected |
| `-assertions {0,1,2}` | 1 | TAB `Assertion` statements: off / warn / fatal |
| `-range_test_initial`, `-range_test_updated {0,1,2}` | 2 / 1 | coefficient bound checks on initial and updated values: off / warn / fatal |
| `-gpzerodivide {0,1}` | 0 | GEMPACK's two-class zerodivide semantics (nonzero-by-zero vs zero-by-zero defaults) |
| `-postsim {0,1}` | 1 | execute the TAB's PostSim section after the simulation |
| `-comp_steps n`, `-comp_do_approx/-comp_do_acc {0,1}`, `-comp_redo {0,1}`, `-comp_redo_min_frac x`, `-comp_sberr_warn {0,1}` | see log | complementarity ([GM] ch. 51) approximate/accurate run controls |
| `-nowrites n` | 0 | suppress output writes |
| `-cofdump {0,1}` | 1 | write the `.cof`/`.cbin` coefficient dump (recorded in `stats.json` `options`) |
| `-verbosity {0,1,2}` | 1 | 0 = errors/warnings + accuracy summary only; 1 = phase progress and timings; 2 = per-rank/per-block debug detail (also exported as `TEEMS_VERBOSITY` for the Fortran kernels; MA48 duplicate-entry notes appear only at 2) |
| `-nox` | — | PETSc: no X output |

**Exit status.** `0` means the run completed and no `Error:` line was
printed. `1` means an `Error:` line was printed: either a named abort
(option validation, TAB/closure/shock/data faults, factorization and
resource failures — `MPI_Abort(…, 1)` or `PetscFinalize` + `return 1`
on the spot) or a failure that was reported and the run continued to
the end (every `Error:` line goes through `errmsg()`, which counts,
and `main()` returns 1 when the count is nonzero; on multi-rank runs
`mpiexec` propagates any rank's nonzero status). Any other status is a
crash, a kill (137 = out of memory under a container limit) or an MPI
launch failure, none of which prints an `Error:` line. The first line
of every log is `teems-solver <version>`, unconditionally.

teems-R populates these from `ems_solve()`/`ems_RK()` arguments
(`solution_method`, `matrix_method`, `n_tasks`, `steps`,
`n_subintervals`, and the validated expert dots — `laA/laD/laDi`,
`fastrefac`, `condest`, `assertions`, `range_test_*`, RK controls,
complementarity controls). Options travel on the invocation, never in
the CMF (which is a file manifest only), and are echoed in
`stats.json`'s `options` object. Stale flags from older commands
(`-regset`, `-enable_time`, `-presol`, `-nesteddbbd`) are ignored.

## 12. Verifying a build

- **Determinism**: repeated runs of the same binary on the same inputs
  are bit-identical for every solution and matrix method and rank
  count, so a locally built image can be checked against a reference
  run by comparing output files byte for byte (the maintainers keep
  such golden manifests across methods × ranks × subintervals ×
  in-memory modes and gate every change on them).
- **Identity**: `teems-solver -version` (`docker run --rm teems:<tag>
  /opt/teems-solver/solver/teems-solver -version`) names the binary's
  version; the image label `org.opencontainers.image.version` is only
  the pre-pull hint set from the same constant at build time.
- **Smoke test**: the README's *Verification* section runs the shipped
  example through the image; the teems R package's test suite
  exercises the full pipeline against a local `teems:<tag>`
  (`ems_option_set(docker_tag=)`).
- **Structure before solving**: `-solmed probe` (or `ems_probe()`)
  reports the model's ordering statistics and any structural
  singularity by name without solving; `-condest 1` measures each
  linear solve's backward error and scaled condition number under
  sequential LU (§10).
- Build details: `src/BUILDING.md` (HSL staging, `src/patches/`
  applied by `mp48_mod.sh`, serial make requirement, `OPT` knob) and
  the README's expedited-image recipe.

## 13. Known limitations and planned work

- **DBBD with subintervals on large shocks can hit a structurally singular
  border Schur complement mid-path.** Seen on the 1.4M-equation GTAPv7
  R32 x medium rig with a uniform `aoall` +20 shock: `-matsol 2
  -nsubints 2` fails in the second subinterval with MA48 `INFO(1) = -4`
  on the 8,16x-row border system (rank one short; the on-failure probe
  names the under-determined variable and the over-constrained
  equation), while `-nsubints 1` and `4` solve. A flow reaches exactly
  zero along the subinterval path and the realized border pattern loses
  connectivity -- the same class as the value-dependent pattern the
  `-fastrefac` per-step test guards against. Reproduced twice on the
  calibration host, not on the development build with the identical
  inputs, so it is host-deterministic (2026-09). Remedies, as the
  diagnosis prints: change the subinterval count, or use a method that
  does not factorize that border (LU, or SBBD/NDBBD for intertemporal
  systems). Reproducer inputs are kept with the HPC archive.

- **Error handling** is largely `printf` + `exit`/`return`; parse errors
  on one rank can abort non-collectively. Messages carry
  `Error:`/`Warning:` prefixes with file/flag/remedy named, and teems-R
  maps them to typed conditions; the exit status is checked as well as
  the log.
- Fixed-size line buffers (`TABREADLINE` = 20000) bound statement length
  (the fatal names the statement's head; teems-R Backsolve/Substitute
  expansion is the usual cause);
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
- **Version handshake**: the solver image and the R package are
  versioned independently. The image moves 1.0.0 → 1.1.0: everything
  since 1.0.0 is additive for the package that drives it (the `hsl`
  compatibility symlink, deprecated aliases such as `-solmed Mmid`,
  and ignored legacy flags are kept precisely so that the released
  1.0.0-era package still runs bit-identically against it), while a
  newer R package may *require* a minimum image (the coefficient dump
  needs ≥ 1.1.0). The solver side of the handshake is in place: one
  constant (`src/version.h`) surfaced by `-version`, the ungated log
  banner, the `solver_version` field of `stats.json`/`probe.json`
  and the image label. The R package's pre-flight check against it
  (a minimum-version fence; no response = pre-1.1 image) is the
  remaining half. The symlink and aliases are removed only at a real
  2.0.0.
- **Calibration**: the `matrix_method`/`n_tasks` auto rules, the
  `-fastrefac` and MA48 pivot-threshold default flips, the NDBBD
  `-inmemory` default, RK-vs-Gragg rankings, and the structural-probe
  thresholds were fit at dev-box scale (≤ 4 ranks, ≤ 8M equations);
  a many-core calibration program on instance-scale rigs is the next
  step and will re-anchor those defaults.
