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
| `BUILD_DIR` | `/opt/teems-solver` | Installation directory for solver components |
| `MPICH_VERSION` | `5.0.1` | MPICH version to install |
| `PETSC_ARCH` | `native` | PETSc architecture identifier |
| `PETSC_VERSION` | `v3.25.0` | PETSc version tag or branch |
| `TEEMS_TMP` | `/tmp/teems` | Temporary directory for builds (removed after build) |


## Directory Structure

After build, the image contains:

```
/opt/teems-solver/
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