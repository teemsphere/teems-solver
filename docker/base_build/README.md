# TEEMS Solver Base Layer

[![License](https://img.shields.io/badge/License-GPL-blue.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-0.92-green.svg)](https://github.com/matthewcantele/teems-solver/releases)

This Dockerfile builds the base layer for the TEEMS solver, containing all open-source dependencies required for the full solver build. This base layer approach significantly reduces build times for end users who only need to add their HSL libraries to complete the installation.

## Overview

The base layer includes:
- **OpenMPI** (v5.0.7) - Message Passing Interface implementation
- **MPICH** (v4.3.0) - Alternative MPI implementation
- **PETSc** (v3.24.4) - Portable, Extensible Toolkit for Scientific Computation
- **Build tools** - gcc, gfortran, make, autotools

## Build Arguments

The following build arguments can be customized:

| Argument | Default | Description |
|----------|---------|-------------|
| `BUILD_DIR` | `/opt/teems-solver` | Installation directory for solver components |
| `OPENMPI_VERSION` | `5.0.7` | OpenMPI version to install |
| `OPENMPI_BASE_VERSION` | `5.0` | OpenMPI major.minor version for download path |
| `MPICH_VERSION` | `4.3.0` | MPICH version to install |
| `PETSC_ARCH` | `native` | PETSc architecture identifier |
| `PETSC_VERSION` | `v3.24.4` | PETSc version tag or branch |
| `TEEMS_TMP` | `/tmp/teems` | Temporary directory for builds |


## Directory Structure

After build, the image contains:

```
/opt/teems-solver/
├── src/                    # TEEMS solver source code
├── lib/
│   ├── mpi/               # OpenMPI/MPICH installation
│   ├── petsc/             # PETSc installation
│   ├── ma48/              # Prepared for HSL MA48
│   ├── ma51/              # Prepared for HSL MA51
│   ├── mc66/              # Prepared for HSL MC66
│   └── mp48/              # Prepared for HSL MP48
└── LICENSES/
    ├── openmpi/           # OpenMPI license
    └── petsc/             # PETSc license
```

## License

This project is licensed under the GPLv3.0 License - see the [LICENSE](LICENSE) file for details.

Individual components retain their original licenses:
- OpenMPI: BSD 3-Clause License
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