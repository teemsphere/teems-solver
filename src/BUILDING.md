# Building the TEEMS solver

The solver is built inside a container based on `matthewcantele/teems_base`
(Debian + MPICH + PETSc under `/opt/teems-solver/lib`). The proprietary HSL
packages (MA48, MA51, MC66, MP48) are supplied as tarballs at image build
time; see `docker/expedited_build/Dockerfile` for the canonical flow:

1. Build the four HSL libraries (`autoreconf && configure && make`).
   MP48 must be built with `make` **serially** — its Fortran module
   dependencies break under `make -j`.
2. Stage `hsl_mp48d.f90`, `ddeps.f`, and `hsl_mp01.mod` from the MP48 build
   into `src/`, plus `hsl_mc66d.f90` and `ddeps90.f90` from MC66.
3. Run `./mp48_mod.sh`, which applies `patches/*.patch` to the staged
   sources (64-bit work arrays; MA48/MA50/MC13/MC21/MC29/MC59/MC71 entry
   points renamed to `Z*` so the patched static copies cannot clash with
   the dynamically linked libma48/libma51). `patch --forward` fails loudly
   if the upstream HSL sources ever drift.
4. `make` in `src/` produces the `hsl` binary.

## Makefile knobs

- `BUILD_DIR` (default `/opt/teems-solver`) — root of the toolchain tree.
- `OPT` (default `-Ofast`) — optimization level. **Note:** `-Ofast` enables
  `-ffast-math` (FP reassociation, no NaN/Inf guarantees). This is the
  tested production configuration; use `make clean && make OPT=-O3` to
  produce an IEEE-conformant binary for numerical cross-checks.
- `make clean` removes objects and the binary (PETSc's clean also removes
  `*.mod`; the staged `hsl_mp01.mod` is restored automatically).

## Runtime notes

- Solver scratch files go to `/tmp/` by default; override with the
  `-tempdir <dir>` option or the `TMPDIR` environment variable.
- The `-smllthreads` option caps OpenMP threads in selected solver
  sections; `-maxthreads` sets the global thread count.

## -inmemory guidance

`-inmemory 1` keeps the value arrays resident (skipping the per-step
scratch spills) and, unless `-tempdir`/`TMPDIR` is set, places scratch on
tmpfs. Default is method-dependent: **on** for LU and SBBD (which write
nothing to scratch in this mode — MP48 holds SBBD factors in memory,
sized by laA/laDi), **off** for DBBD and NDBBD, whose block-factor file
traffic wants the page cache that array spilling frees. Measured on a
4.4M-equation intertemporal run (2 ranks, 30GB node): SBBD ~7% faster
with no scratch debris; DBBD neutral; NDBBD ~6-11% slower with residency
forced on. The sign for DBBD/NDBBD depends on factor-file traffic vs
free RAM, so on large-memory nodes forcing `-inmemory 1` may win there
too.
Override with an explicit `-inmemory 0/1`; a memory-availability check
falls back to spilling when the resident estimate exceeds half of
MemAvailable.
