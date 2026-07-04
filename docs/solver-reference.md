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
  `substep_base + 2·Δ` (see §4).
- `swap` (closure) and shock statements — parsed by `closure_read()` /
  `shocks_read()` into per-element `closure_entry` records.

Known deviations from GEMPACK are printed at startup, most notably:
shocks to 2-dimensional variables use transposed orientation relative to
GEMPACK, and intertemporal variables should be declared with minimal
dimensionality to keep the border (netcut) small.

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
   method; computes `netcut` (border size) [HK16 §5].
8. **Jacobian preallocation and fill** (`jacobian_preallocate`,
   `jacobian_fill`) — PETSc AIJ matrices A (endogenous) and B
   (exogenous columns → RHS).
9. **Solve** (§4/§5), **updates**, iterate per step/subinterval.
10. **Outputs** — CSVs (`outputs_write_csv`) and solution binaries.

## 4. Solution methods (`-solmed`, `enum solution_method`)

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

## 5. Matrix methods (`-matsol`, `enum matrix_method`)

| value | method | domain | notes |
|---|---|---|---|
| 0 `MM_LU` | serial MA48 LU | any | single rank does the factorization |
| 1 `MM_SBBD` | singly bordered block diagonal via HSL_MP48 [HK16] | **intertemporal** | blocks by time period; MC66 ordering optional (`-withmc66`), direct ordering by default; factors held in MP48 memory (sized by `-laA/-laDi`) |
| 2 `MM_DBBD` | doubly bordered block diagonal [KH19] | **static** (regional blocks) | C-side block factorization (MA48 kernels), factor handoff via scratch files |
| 3 `MM_NDBBD` | nested DBBD [KH19] | **intertemporal, T ≫ R** | regional blocks nested inside time blocks; requires `-regset`; `-nesteddbbd 1` |

Production guidance (measured, see §8): SBBD is the primary intertemporal
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

## 6. Parallelism

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

## 7. Scratch, memory, and the in-memory mode

All solver temp files go under `scratch_dir` — default `/tmp/`, override
with `-tempdir <dir>` or `TMPDIR`; exported to the Fortran kernels as
`TEEMS_SCRATCH`. User run directories receive only real outputs plus the
small per-rank preprocessed TAB copies.

Two spill classes exist:

1. **Driver array spills** (`solve_drivers.c`): before each factorization
   the value arrays (`elem_vals`, `closure_vals`) and Mmid step state
   (`clag1`, `varchange`, `xc*`) are written to `_temp*.bin` files and
   freed, then re-read after the solve.
2. **Block factor handoff** (DBBD/NDBBD): the Fortran kernels write LU
   factors per block (`_vav/_irnv/_keep`), and the C side writes
   interface products (`_bivi/_rbvi/_cbvi`, `_rank/_row/_col`); re-read
   during back-solve, removed afterwards.

`-inmemory 1` disables class 1 entirely (arrays stay resident) and
relocates scratch to tmpfs (`/dev/shm`) unless `-tempdir`/`TMPDIR` is
pinned. **Default: on for LU, SBBD, and DBBD; off for NDBBD** — NDBBD's
factor traffic (measured 17× peak RSS) competes for the page cache that
array spilling frees, and forcing residency measured 6–11% slower at
4.4M equations on a 30GB node, while SBBD gained ~7% and DBBD was
neutral on both domains. A startup check estimates residency cost and
falls back to spilling (with a warning) if it exceeds half of
`MemAvailable`. Explicit `-inmemory 0/1` always wins.

## 8. Measured characteristics

Benchmarks (GTAP12 database; 20-core/30GB node; 2 ranks unless noted):

| benchmark | size | config | wall | peak RSS/rank |
|---|---|---|---|---|
| static R32×medium | 1.35M eq | LU-1, Johansen | 12.4s | 0.41–0.48GB |
| static R32×medium | 1.35M eq | DBBD-2, Johansen | 20.4s | 0.61–0.66GB |
| intertemporal WB7×medium, T=21 | 4.36M eq | SBBD-2, Mmid | 107–115s | 1.1–1.6GB |
| intertemporal WB7×medium, T=21 | 4.36M eq | NDBBD-2, Mmid | 314–334s | 1.1–1.3GB |
| intertemporal WB7×medium, T=21 | 4.36M eq | DBBD-2, Mmid (off-domain) | ~420s | 2.1–2.6GB |

Scratch traffic before the in-memory work (202k-eq model, per solve):
SBBD 274MB, NDBBD 1.36GB — the motivation for §7.

Determinism: repeated same-binary runs are bit-identical for every
method (basis of the golden-run verification, below).

## 9. Numerical characteristics

- Solve precision `solve_real` = double; storage `store_real` = float
  (halves memory traffic for coefficient values; solution vectors and
  factorizations are double).
- Compiled `-Ofast` (includes `-ffast-math`: FP reassociation, no
  NaN/Inf guarantees). This is the tested production configuration;
  `make OPT=-O3` builds an IEEE-conformant binary for cross-checks.
- `zerodivide` defaults substitute configured values on 0/0 per [GM].
- Levels bounds (`bound_type`) are enforced during formula evaluation.

## 10. Command-line reference

| option | default | meaning |
|---|---|---|
| `-cmdfile <path>` | `./reg.cmf` | CMF file |
| `-matsol {0,1,2,3}` | 0 | matrix method (§5) |
| `-solmed <name>` | `Mmid` | solution method (§4) |
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
| `-inmemory {0,1}` | 1 except NDBBD | §7 |
| `-nsbbdblocks n` | 2 | SBBD block-count hint |
| `-ndbbd_bl_rank n` / `-nestfile <csv>` | — | NDBBD interface-block rank override |
| `-stoiter n` | 1 | stochastic repetitions |
| `-nowrites n` | 0 | suppress output writes |
| `-nox` | — | PETSc: no X output |

teems-R populates these from `ems_solve()` arguments
(`solution_method`, `matrix_method`, `n_tasks`, `steps`,
`n_subintervals`, `laA/laD/laDi`, `append_args` for anything else).

## 11. Verification and development infrastructure

- **Golden runs** (`.audit/verify.sh`, dev-machine): rebuilds from the
  working tree in the `teems-audit` container and checks five solves
  bit-identically against manifests anchored to the pre-refactor binary:
  GTAPv7 static LU/Johansen; GTAP-RE intertemporal LU/Mmid; and GTAP-RE
  SBBD, DBBD, NDBBD at 2 ranks. Any behavior change fails the gate.
- **Benchmark rigs**: `.audit/bench_run.sh` (wall/RSS via `time -v`),
  `.audit/strace_run.sh` (per-file write-byte accounting), deployments
  under `.audit/bench-*` produced by teems-R scripts.
- Build: see `BUILDING.md` (HSL staging, `src/patches/` applied by
  `mp48_mod.sh`, serial make requirement, `OPT` knob).

## 12. Known limitations and planned work

- **NDBBD without `-regset` aborts** with a bare MPI error (exit 59)
  instead of a diagnostic.
- **Error handling** is largely `printf` + `exit`/`return`; parse errors
  on one rank can abort non-collectively.
- Fixed-size line buffers (`TABREADLINE` = 20000) bound statement length;
  overflow warnings exist but are not hard guards everywhere.
- The binary is named `hsl` for historical reasons; renaming to
  `teems-solver` requires coordinated changes in teems-R and the images.
- Planned (measured motivation on file): formula op-list caching across
  steps (removes per-step TAB re-parsing); `formula_op` compaction
  (~1.4KB/op today — cache-hostile in the evaluator); factor handoff
  through memory instead of files (would let DBBD/NDBBD run resident
  without page-cache competition); verbosity-gated logging (3,000+
  unconditional printfs); NDBBD-vs-SBBD crossover mapping at scale.
