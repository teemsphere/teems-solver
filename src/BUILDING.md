# Building the TEEMS solver

The solver is built inside a container based on `matthewcantele/teems_base`
(Debian + MPICH + PETSc under `/opt/teems-solver/lib`). The proprietary HSL
packages (MA48, MA51, MA60, MC71, MC66, MC79, MP48) are supplied as a single
**libHSL** source snapshot (`libHSL.v<version>.tar.gz` from
https://licences.stfc.ac.uk) at image build time; see
`docker/expedited_build/Dockerfile` for the canonical flow:

1. `docker/stage_hsl.sh` stages the per-package sources out of the libHSL
   snapshot: the main decks are copied and the historical `ddeps` bundles
   are reproduced by concatenating the individual dependency decks (the
   recipes are byte-exact against the original per-package tarballs).
   Every deck used is checksum-pinned to the verified snapshot (2026.8.4);
   a snapshot that changes any used deck fails the build until it has been
   re-verified and its checksums recorded in `stage_hsl.sh`.
2. Build the two dynamically linked libraries — `libma48.so`, `libma51.so`
   (a direct `gfortran -shared` of deck + deps; same flags, soname, and
   `.libs/` layout the old autotools build produced) — and compile
   `common90.f90` (HSL_MP01) with `mpif90` for its module file.
3. Stage `hsl_mp48d.f90`, `ddeps.f`, and `hsl_mp01.mod` into `src/`, plus
   `hsl_mc66d.f90` and `ddeps90.f90` from MC66, plus `hsl_mc79i.f90`,
   `hsl_mc79i_ciface.f90` (F95 core + its C interface) and
   `include/hsl_mc79i.h` from MC79 (integer-only, dependency-free), plus
   `ma60d.f`, `mc71d.f`, and `fd15d.f` for `-condest`.
4. Run `./mp48_mod.sh`, which applies `patches/*.patch` to the staged
   sources (64-bit work arrays; MA48/MA50/MC13/MC21/MC29/MC59/MC71 entry
   points renamed to `Z*` so the patched static copies cannot clash with
   the dynamically linked libma48/libma51). `patch --forward` fails loudly
   if the upstream HSL sources ever drift.
5. `make` in `src/` produces the `teems-solver` binary.

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
