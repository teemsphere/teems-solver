# Building the TEEMS solver

The solver is built inside a container based on `matthewcantele/teems_base`
(Debian + MPICH + PETSc under `/opt/teems-solver/lib`). The proprietary HSL
packages (MA48, MA51, MA60, MC71, FD15, MC66, MC79, MP48) are supplied as a
single **libHSL** source snapshot (`libHSL.v<version>.tar.gz` from
https://licences.stfc.ac.uk) at image build time, passed to the build as
`--build-arg PATH_LIBHSL=<path>`. `docker/expedited_build/Dockerfile` is the
canonical flow; `docker/full_build` also builds MPICH and PETSc.

## ISA level

The base image records the ISA level it was compiled at in
`/opt/teems-solver/archflags`: `-march=x86-64-v2` for the published
`amd64` base, `-march=armv8-a` for `arm64`. The expedited build compiles the
HSL libraries and the solver at the same level (`make OPT="-Ofast
$(cat /opt/teems-solver/archflags)"`). The full build takes
`--build-arg MARCH=<level>` to compile everything for one machine; higher
levels showed no measurable gain, since the dense kernels dispatch at run
time through OpenBLAS. The runtime image pins `OPENBLAS_CORETYPE` to the
kernel family paired with this level (`--build-arg BLAS_CORETYPE`
overrides it for a base built at another level); see
`solver-reference.md` §10.

## Steps

1. `docker/stage_hsl.sh` stages the per-package sources out of the libHSL
   snapshot: the main decks are copied and the historical `ddeps` bundles
   are reproduced by concatenating the individual dependency decks (byte-exact
   against the original per-package tarballs). Every deck used is
   checksum-pinned to the verified snapshot (2026.8.4); a snapshot that
   changes any used deck fails the build until it has been re-verified and
   its checksums recorded in `stage_hsl.sh` (`ALLOW_HSL_DRIFT=1` overrides).
2. Build the two dynamically linked libraries, `libma48.so` and `libma51.so`
   (a direct `gfortran -O2 -fPIC -shared` of deck + deps, with the soname
   and `.libs/` layout of the old autotools build), and compile
   `common90.f90` (HSL_MP01) with `mpif90` for its module file.
3. Stage `hsl_mp48d.f90`, `ddeps.f` and `hsl_mp01.mod` into `src/` and run
   `src/mp48_mod.sh`, which applies `src/patches/hsl_mp48d.patch` and
   `src/patches/ddeps.patch` to those two sources: 64-bit
   duplicate-detection work arrays and the CGLOB allocation fix in MP48,
   and the MA48/MA50/MC13/MC21/MC29/MC59/MC71 entry points renamed to `Z*`
   so the statically patched copies cannot clash with the dynamically
   linked libma48/libma51. `patch --forward` fails loudly if the upstream
   HSL sources drift.
4. Stage the rest into `src/`: `hsl_mc66d.f90` and `ddeps90.f90` from MC66;
   `hsl_mc79i.f90`, `hsl_mc79i_ciface.f90` and `include/hsl_mc79i.h` from
   MC79 (the structural probe; integer-only, dependency-free); `ma60d.f`,
   `mc71d.f` and `fd15d.f` (FD15 is its own deck) for the `-condest`
   diagnostics.
5. `make teems-solver teems-solver-f64` in `src/`. Run `make` serially: the
   Fortran module dependencies of the HSL sources are not expressed for a
   parallel build.

The solver build in the `teems-audit` container used by the `.audit` kits
(`.audit/verify.sh`) follows the same steps against the libraries the
image already holds.

## Makefile (`src/makefile`)

- Targets: `teems-solver`, and `teems-solver-f64`, the same sources with
  double-precision coefficient storage (`precision = "f64"` in teems-R;
  the default binary stores coefficients in single precision and solves
  in double).
- `BUILD_DIR` (default `/opt/teems-solver`): root of the toolchain tree.
- `OPT` (default `-Ofast`): optimization level. `-Ofast` enables
  `-ffast-math` (floating-point reassociation, no NaN/Inf guarantees); the
  solver tests NaN/Inf on the bit patterns, so its arithmetic-error checks
  hold under it. This is the tested production configuration; `make clean
  && make OPT=-O3` gives an IEEE-conformant binary for numerical
  cross-checks.
- `WARN` (default `-Wall`): warning flags on top of `OPT` for our own
  sources; the tree builds warning-free and `.audit/verify.sh` holds the
  count at 0. The staged HSL sources are compiled `-w` as inputs.
- `make clean` removes objects and both binaries (PETSc's clean also
  removes `*.mod`; the staged `hsl_mp01.mod` is restored automatically).
