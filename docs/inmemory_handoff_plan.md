# In-memory R↔solver handoffs (ROADMAP 6.13) — scoping report

Status: SCOPED 2026-08-16; §6 recommendation ACCEPTED by the user
2026-08-17 (dump all coefficients, `double`, value-only for now,
`ems_deploy(write_coefficients = FALSE)`), and the coefficient channel
LANDED the same day — see the ROADMAP 6.13 record for the shipped
shape and gates. §7 below records the outcome; §1–§5 are the scoping
evidence as measured before the change.

The 6.13 row asked three things before any build: (1) inventory the
solver read paths that assume a seekable `FILE*`; (2) measure the real
serialize/parse overhead on big deploys to size the win; (3) decide the
container transport. All three are answered here. Bench scripts lived
in the session scratchpad (bench613.R, ab_condlu.sh); the rigs are the
GTAP12 flexAgg data with the GTAPv7 internal model, plus two existing
`.audit/` deploys.

## 1. What the round trip actually costs

Machine: 20-CPU laptop, 30 GB, ext4 NVMe; teems:dev (solver 7b4703d),
teems-R working tree at 8ae6288 (post write-perf), Johansen 1-step LU
unless noted. Wall seconds, single warm runs (spread ≤10% where
repeated).

| rig | system size | R `ems_deploy` | deploy dir | solver input parse¹ | solver formulas¹ | docker+mpiexec launch | solve (step) | solver output write² | binaries | CSVs | R compose vars-only | R compose all |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| small (big3×macro_sector) | 3.5k eq | 0.25 | 0.4 MB (11 files) | <0.05 | | 0.5–0.7 | ~0.0 | ~0 | 1.5 MB | 0.09 MB / 264 | 0.23 | 0.65 |
| mid (R32×medium) | 1.41M eq | 0.39 | 4.8 MB | ~0.1 | ~0.1 | 0.5–0.7 | 16–31 | 0.2–0.5 | 15.6 MB | 25 MB / 264 | 0.24 | 5.1 |
| full-res (160×65, GTAP12) | 40.5M eq | 2.2–3.2 | 91.5 MB | 1.3 | 2.1 | 0.5 | (not run: hours) | — | — | — | — | — |
| `.audit/run-wideb-cond-lu` (condensed, 1.3M backsolved) | 56k eq | — | 4.5 MB | ~2.5 (incl. formulas) | | 0.4 | 7.0 | **4.6** | 14 MB | **277 MB / 792** | 0.59 | **33.1** |
| `.audit/run-wideb-inter-sbbd` (GTAP-RE, 21 steps) | — | — | 14.4 MB | — | — | — | — | — | 41 MB | **1.2 GB / 889** | 7.9³ | **243** (7 GB R object) |

¹ From a per-line timestamped log at full-res: CMF/TAB preprocess +
set parse 0–0.7 s, `data_read_files` over 91 MB of iodata text
0.7–2.0 s (~1.3 s), formulas 2.0–4.1 s. "Variable calculation time"
in the solver log is the sum. ² `-nowrites 0` vs `1`, interleaved 2×2
on the cond-lu rig: 22.9/22.9 vs 18.4/18.3 s wall — the CSV dump
(792 files, and one full re-scan of the preprocessed TAB per outdata
entry, `outputs_write_csv` cmf_io.c:988) costs the solver ~4.6 s.
³ vars-only 7.9 s on the intertemporal rig is R-side variable
composition (`CJ` over 265 variables × 21 periods), not parsing —
the five binaries parse in well under a second everywhere.

Two control measurements:

- **tmpfs vs disk makes no difference to the deploy write.** Full-res
  `ems_deploy` to `/dev/shm`: 2.87 / 2.21 s; to ext4: 2.63 / 2.24 s
  (interleaved). The write hop is CPU (formatting 91 MB of text with
  fwrite's formatter), not I/O; the page cache already absorbs the
  I/O. A memory-backed *file* transport therefore saves nothing on the
  R side — only skipping the text serialization altogether (binary
  arrays) would.
- **Container launch is a fixed ~0.5–0.7 s** (`docker run --rm true`
  0.31–0.42 s; `docker run … mpiexec -n 1 teems-solver -version`
  0.48–0.67 s). For the small rig that is 60–70% of the 0.8–0.95 s
  `ems_solve` wall; for anything the solver takes seconds on it is
  noise.

### Reading the table

- **R→solver (input hop)**: R serialize 0.25–3.2 s plus solver text
  parse ≤1.3 s, i.e. ≤4.5 s at 40M equations, ≤0.5 s at 1.4M. The
  matrix calc alone at full-res is 29 s and the LU factorization hours;
  at mid the single Johansen step is 16–31 s and a Gragg 2-4-6 run is
  ~12 steps. The input hop is **<5% of wall on any model where the
  solver time is noticeable**, and irrelevant on the small models
  where it isn't (0.3 s vs a 0.5 s container launch).
- **solver→R, binaries** (`.mds/.var/.set/.sel/.bin`): 15–41 MB,
  parsed by the cpp11 reader in ≤0.6 s. Not a cost. (`sol.set` is
  98% padding — 20 000-byte `readele` per set — but 1.3 MB is
  harmless.)
- **solver→R, coefficient CSVs**: **the** cost, by one to two orders
  of magnitude. 25 MB → 4.9 s of R; 277 MB → 4.6 s solver + 32.5 s R;
  1.2 GB → 235 s R and a 7 GB R object. R side is `readLines` +
  per-line `strsplit`/`as.numeric` in `.compose_coeff`/
  `.parse_coeff_block` (teems-R R/coeff_compose.R, coeff_block_parse.R).
  Also a **precision defect**: `outputs_write_csv` prints every value
  with `%f` (cmf_io.c:1094/1114/1116) — six fixed decimals, so
  coefficients below ~5e-7 come back as 0.000000 and small elasticities/
  shares lose relative precision; the binaries carry `store_real`
  intact.

## 2. Solver read-path inventory (item 1 of the scope)

Full inventory in the session record; the load-bearing facts:

- 181 `fopen` sites in `src/*.c` (132 read, 49 write); **zero**
  `open()/tmpfile/mkstemp/fmemopen/open_memstream/fscanf/getline`;
  Fortran does no file I/O. Of the 132 reads, 82 touch run inputs
  (CMF/TAB/iodata/closure/shock/_ps), 47 are binary scratch the solver
  wrote itself (`-inmemory`-gated spills in solve_drivers.c/
  block_solve.c, plus the NDBBD `_rank/_row/_col/_bivi…` presolve
  handoff which is *not* gated), 3 other.
- **Seekability is a non-issue.** Only 7 `fseek/rewind/ftell` sites
  exist; 5 are on solver-written outputs (`stats.json` read-modify-
  write patchers main.c:161/201, `_tempsol` spill) and 2 are the
  `rewind(fin)` in `tab_postsim_split` (cmf_io.c:1432/1460), which
  `fmemopen` supports anyway. Every other input reader is a forward
  `fgets` stream through `tab_next_statement*` (tab_parse.c:5215–5350,
  50 call sites) or one of 54 raw `fgets` sites (rewrite passes,
  `datafile_read_header_info/_labels`, the four `data_read_files`
  format branches).
- The real obstacles are (a) **repeat-open-by-name semantics** — the
  CMF is re-scanned 9× (`cmf_count_files` ×3 + `cmf_read` ×6), the
  preprocessed TAB is re-opened once per reader with no handle caching:
  ~35 setup scans plus `jacobian_fill`/`formulas_execute`/
  `updates_apply`/`assertions_execute`/`backsolve_recover` **per
  solve step** (a Gragg 2-4-8 run re-parses the TAB 50–80×), plus
  once per outdata entry in `outputs_write_csv` (268–889×); and (b) the
  **write-then-re-read rewrite chain** — `tab_preprocess` (2 passes) →
  `tab_defaults_validate` → `tab_postsim_split` → `tab_setbuilder_
  transform` → `tab_complementarity_transform` → `tab_levels_transform`
  → `tab_write_variables`, seven intermediate files and four
  `rename()` commits, all in the bind-mounted CWD (`_temp_tab_file000N
  .tab`/`_temp_tab_new_file000N.tab` are never removed). Under
  DBBD/NDBBD (`rank_hsl==rank`) **every rank** runs the CMF+TAB
  preprocess and writes its own per-rank scratch TAB. The C3
  comp-accurate re-entry (main.c:1679/2355) re-runs ordering +
  preallocation + the whole driver, so it re-opens the TAB again.
- Contract: everything is discovered from one PETSc option `-cmdfile`
  → `cmf_read` (cmf_io.c:193); the CMF is a pure file manifest
  (`iodata/outdata/soldata/tabfile/closure/shock`), controls are CLI.

Cost of an in-memory input transport, honestly sized: a name→buffer
registry handing out a fresh `fmemopen` per request (the 82 read
sites become one helper call each — mechanical), the rewrite chain
moved to `open_memstream` buffers (the actual work, all in cmf_io.c +
levels.c, and the C3 re-entry + per-rank duplication live there), and
a host↔container ingress for the buffers (§4). Weeks of solver work,
new failure surface in the most-fuzzed code in the tree (6.11 A(d)
lived in exactly these functions), for ≤4.5 s at 40M equations.

## 3. R-side transport inventory (what would have to change)

- `ems_deploy()` writes: `<model>.tab` (R-finalized, with auto
  `File (new)`/`Write` pairs for **every** set and coefficient — this
  is what makes 268–889 outdata lines), `.cls`, `<id>_HHMM.shf`,
  `.cmf`, one `<FILE>.txt` per TAB File statement (`.ragged_write`
  appends), `metadata.rds`, `model_diagnostics.txt`, empty `out/`
  tree. Returns only the CMF path (+ `tab_path` attr); everything is
  re-derived from disk later (`readLines` TAB for `(intertemporal)`,
  `readRDS(metadata.rds)`, `.retrieve_cmf`, compose re-parses the TAB
  twice for intertemporal).
- `ems_solve()`: one bind mount (deploy dir → `/opt/teems`), `docker
  run --rm … mpiexec -n N teems-solver -cmdfile … 2>&1 | tee
  out/solver_out_HHMM.txt`; **exit code discarded** (`system()` return
  ignored, and the pipeline's status would be `tee`'s) — completion
  and failure are detected solely by log grep in `.check_solver_log`.
- `ems_compose()`: five binaries via `src/parse_solution.cpp` (bare
  struct dumps of `array_def`/`set_def`/`set_element` at main.c:2424–
  2482, no versioning), `out/sets/*.csv` for canonical set casing +
  consistency, `out/coefficients/*.csv` + `out/postsim/*.csv` via
  pure-R line parsing, `metadata.rds` (`reference_year`), TAB (labels,
  set-name lists for column naming), `stats.json`.
- Post-solve, only the TAB, CMF, metadata.rds and `out/` are read
  back; `.cls`, `.shf` and the iodata `.txt` are solver-only
  (`.retrieve_tab_comp` computes an `input_files` attr nothing reads).
- Latent: `parse_solution.cpp` hardcodes `glval` as `float`; under
  the f64 binary it is `double` (same 456-byte `array_def` by
  padding luck), so `glval` is misread for `precision="double"` —
  harmless today because compose never consumes it.

## 4. Container transport (item 3 of the scope)

| option | verdict |
|---|---|
| Bind-mounted tmpfs / `/dev/shm` deploy dir | Already possible today via `ems_option_set(tempdir=)`; measured to change nothing (§1) because neither hop is I/O-bound. Not worth a default change. |
| stdin/socket/memfd protocol into the existing `docker run` | Only pays if the solver side is rebuilt per §2; still leaves the container launch (0.5–0.7 s) and the text-parse CPU. Reject. |
| Solver as a shared library loaded by R | Would remove the launch cost and the mpiexec boundary for n_tasks 1, but the container exists precisely to ship the PETSc/MPICH/HSL/gfortran toolchain (6.6 image contract, ISA-tagged images); an R-loadable .so means shipping that toolchain per platform, and HSL licensing terms constrain redistribution in a linked form. Reject for this horizon; revisit only if a per-run cost of ~0.5 s ever matters (many-run sensitivity loops), and then prefer a **persistent container + `docker exec`** which removes the same 0.35 s `docker run` cost with no build-contract change. |
| Keep files as transport, change **what** is transported | Recommended (§6). |

## 5. Ancillary findings (cheap fixes regardless of the decision)

1. `%f` coefficient CSVs (§1) — lossy; goes away with a binary channel,
   or `%.9g` as an interim (bit-changes every CSV → compose goldens).
2. Exit code discarded (`solve_implement.R:65–76`, `ems_probe.R:171–
   185`); the pipe hides the solver status. `set -o pipefail` in the
   `bash -c` string, or `${PIPESTATUS[0]}`, then check `system()`'s
   return in R.
3. Solver leaves `_temp_tab_file000N.tab`/`_temp_tab_new_file000N.tab`
   in the deploy dir; the `_ps` split file too.
4. `datafile_read_labels` (cmf_io.c:158) has no NULL check on `fopen`
   — a missing iodata file segfaults there rather than aborting cleanly.
5. Redundant TAB re-scans (§2) — a one-shot solver-internal cache
   (read the preprocessed TAB into memory once after `tab_write_
   variables`, hand `fmemopen` streams to every reader) would remove
   50–80 parses per run + 268–889 in `outputs_write_csv`, all behind
   the existing `FILE*` seam, no contract change. Worth measuring, but
   the per-scan cost is small (a 150–480 KB TAB) — likely ≤1 s even
   on the CSV-heavy runs; not a priority.

## 6. Recommendation (needs user decision)

The premise of 6.13 — that the filesystem round trip is the overhead —
does not survive measurement. The input hop is ≤5% of any run where
time matters and is CPU-bound formatting/parsing that a memory
transport would not remove; the solution binaries are cheap and
already ABI-mirrored. **The overhead is the coefficient CSV channel**,
which is 5–235 s of R time plus 0.3–4.6 s of solver time on the runs
measured, scales with model size × number of coefficients, is lossy,
and inflates R memory (7 GB compose object on the intertemporal rig).

Proposed reframing — **6.13 → "Binary + selective coefficient
channel"**, three slices, smallest first:

- **(a) Selective coefficient emission (R-only, days).** `ems_deploy()`
  (or `ems_solve()`) takes a coefficient selection; R writes
  `File (new)`/`Write` pairs and CMF `outdata` lines only for the
  selected ones (default: a curated short list, or none, with `"all"`
  available). Sets keep their CSVs (compose needs the casing).
  Removes most of the 792–889 files and the R parse time
  proportionally; no solver change; `-nowrites` stays as the blunt
  switch. This alone converts the 33 s / 243 s composes into
  seconds for typical use.
- **(b) Binary coefficient dump (solver + cpp11, ~1 week).** Mirror
  `.var/.bin`: `<stem>.cof` = `array_def` × ncof, `<stem>.cbin` =
  `store_real` × ncofele (coefficient slice of `elem_vals`, plus the
  postsim slots — same array), written next to the solution binaries
  under the same rank-0 guard; `parse_solution.cpp` gains a
  `parse_coefficients` twin (struct layout already mirrored) and
  `.compose_coeff` builds tibbles the way `.compose_var` does
  (`CJ` over set elements, `matsize` check). Exact `store_real`
  values, no `%f`; CSV emission becomes opt-in
  (`outputs_write_csv` kept for GEMPACK-style hand runs). Gate:
  compose goldens equal to CSV compose within `%f` rounding, then
  re-anchor. Fix the `glval` float/double mismatch in the same
  parser pass.
- **(c) Input hop — DEFER / reject as an implementation item.** Record
  the inventory (§2) so it is not re-done; revisit only if a use case
  appears where R↔solver latency (not solver time) dominates — e.g. a
  many-run loop on small models — and then the first lever is a
  persistent container (§4), not an in-memory input format.

Do (5.2) exit-code capture in the same teems-R pass as (a); (5.1)
falls out of (b); (5.3)/(5.4) are one-line solver hygiene.

Open questions for the user:
1. Accept the reframing (drop "no input files unless requested" as a
   goal; keep the deploy dir as the reproducibility artifact)?
2. Default coefficient selection for (a): none / curated per model /
   all-as-today?
3. Slice order — (a) first as a pure teems-R change, or (b) directly
   since it is the durable fix?

## 7. Outcome (2026-08-17)

Shipped as designed in §6 (b), with (a) subsumed: the solver dumps
every coefficient (`<stem>.cof` + `<stem>.cbin`, `-cofdump` default
on, written after PostSim) and `ems_compose(which=)` selects on read
through a cpp11 name-filtered seek — so nothing is selected at write
time and there is no per-file overhead to avoid. `ems_deploy()` no
longer emits coefficient `File`/`Write` pairs unless
`write_coefficients = TRUE`; sets keep their CSVs (compose uses them
for canonical casing). The CSV compose path stays as the fallback when
`sol.cof` is absent (older image; warns once).

Measured after the change (same rigs as §1): the 792-coefficient
condensed rig composes all coefficients in 5.4 s instead of 46 s (33 s
in §1 — the CSV path was slower on the re-run), the mid rig in 2.8 s
instead of 5.1 s, and the solver's post-solve CSV write (4.6 s on the
condensed rig) is gone by default. Values are exact `store_real`
instead of `%f`. The `initial` block, the exit-status ride-along
(shipped: `set -o pipefail` + `solve_err$solver_exit`), and the input
hop deferral are as decided.
