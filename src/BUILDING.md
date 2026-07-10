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
4. `make` in `src/` produces the `teems-solver` binary.

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
scratch spills), keeps the DBBD/NDBBD block factors resident (since the
5.9 handoff redesign; `-inmemory 0` writes the legacy `_vav/_irnv/_keep`
scratch files instead), and, unless `-tempdir`/`TMPDIR` is set, places
scratch on tmpfs. Default is method-dependent: **on** for every method
except NDBBD. LU and SBBD write nothing to scratch in this mode (MP48
holds SBBD factors in memory, sized by laA/laDi) and gain ~7%; DBBD
measures neutral on both static (its production domain) and
intertemporal benchmarks and no longer touches scratch for factors.
NDBBD's **off** default predates the factor-handoff redesign (its
rationale — factor-file traffic competing for page cache — no longer
applies in resident mode, which eliminates that traffic entirely,
~825MB/solve at 202k equations, for +22MB RSS/rank); the default awaits
an idle-machine A/B at 4.4M-equation scale before flipping. NDBBD's
`_bivi/_rbvi/_cbvi` and `_rank/_row/_col` interface spills (~65MB/solve
at 202k eq) remain in both modes. The Johansen disk path also orphans
its spill files (~76MB/run at 1.35M equations) — residency avoids that.
Override with an explicit `-inmemory 0/1`; a memory-availability check
falls back to spilling when the resident estimate exceeds half of
MemAvailable.
