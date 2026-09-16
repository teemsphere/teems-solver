# TEEMS Solver Base Layer

[![License](https://img.shields.io/badge/License-AGPL_v3.0-blue.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-1.1.0-dev.1-green.svg)](https://github.com/matthewcantele/teems-solver/releases)

This Dockerfile builds the base layer for the TEEMS solver, containing all open-source dependencies required for the full solver build. This base layer approach significantly reduces build times for end users who only need to add their HSL libraries to complete the installation.

## Overview

The base layer includes:
- **MPICH** (v5.0.1) - Message Passing Interface implementation
- **PETSc** (v3.25.0) - Portable, Extensible Toolkit for Scientific Computation
- **Build tools** - gcc, gfortran, make, autotools

## Build Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `MARCH` | `x86-64-v2` (amd64), `armv8-a` (arm64) | psABI/ISA level everything in the image compiles at. Never `native` for a published image — that bakes in the build host's CPU |
| `BUILD_DIR` | `/opt/teems-solver` | Installation directory for solver components |
| `MPICH_VERSION` | `5.0.1` | MPICH version to install |
| `PETSC_ARCH` | `native` | PETSc architecture identifier |
| `PETSC_VERSION` | `v3.25.0` | PETSc version tag or branch |
| `TEEMS_TMP` | `/tmp/teems` | Temporary directory for builds (removed after build) |


## ISA level and BLAS kernels

After build, the image contains `archflags` and `blascore` — the ISA
level the layer compiled at, and the OpenBLAS kernel family that pairs
with it (`x86-64-v2` → `NEHALEM`, `x86-64-v3` → `HASWELL`, `armv8-a` →
`ARMV8`). Downstream builds compile at the recorded flags and pin
`OPENBLAS_CORETYPE` to the recorded family: Debian's `libopenblas` is
built `DYNAMIC_ARCH` and otherwise picks kernels from the host CPU at
run time, which made the same image give different last digits on
different machines. The runtime build fails if the two disagree, or if
the library does not actually select the pinned family.

## Directory Structure

After build, the image contains:

```
/opt/teems-solver/
├── archflags               # ISA level this layer compiled at
├── blascore                # OpenBLAS kernel family pairing with it
├── src/                    # TEEMS solver source code
├── lib/
│   ├── mpi/               # MPICH installation
│   └── petsc/             # PETSc build
└── LICENSES/
    └── petsc/             # PETSc license
```

## License

This project is licensed under the AGPL-3.0 License - see the [LICENSE](LICENSE) file for details.

Individual components retain their original licenses:
- MPICH: Custom open-source license
- PETSc: BSD 2-Clause License

## Code Authorship

- **Dockerfile Author**: Matthew Cantele
- **TEEMS Solver Source**: Tom Kompas, Ha Van Pham

## Contact

- **Maintainer**: [Matthew Cantele](mailto:matthew.cantele@protonmail.com)
- **Project Homepage**: [https://github.com/matthewcantele/teems-solver](https://github.com/matthewcantele/teems-solver)
- **Bug Reports**: [https://github.com/matthewcantele/teems-solver/issues](https://github.com/matthewcantele/teems-solver/issues)

## See Also

- [Main TEEMS Solver README](../../README.md)
- [Full Build Dockerfile](../full_build/Dockerfile)
- [Expedited Build Dockerfile](../expedited_build/Dockerfile)