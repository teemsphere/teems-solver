# TEEMS Solver

[![License](https://img.shields.io/badge/License-AGPL_v3.0-blue.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-1.0.0-green.svg)](https://github.com/matthewcantele/teems-solver/releases)
[![Docker](https://img.shields.io/badge/docker-required-2496ED.svg)](https://www.docker.com/)

High-performance optimization solver for the TEEMS (Trade and Environment Equilibrium Modeling System) R package. Built on PETSc with MPI parallelization and HSL sparse linear algebra libraries.

## Overview

The TEEMS solver is a specialized computational engine for large-scale computable general equilibrium (CGE) model runs. It combines:

- **PETSc** (v3.25.0) for portable, scalable scientific computation
- **MPICH** (v5.0.1) for parallel processing
- **HSL libraries** for efficient sparse matrix operations

The solver is distributed as a Docker image to ensure reproducibility and ease of deployment across different systems and platforms.

## Table of Contents

- [Prerequisites](#prerequisites)
- [Quick Start](#quick-start)
- [Installation](#installation)
  - [Docker Setup](#docker-setup)
  - [HSL Libraries](#hsl-libraries)
  - [Build Options](#build-options)
    - [Expedited Build (Recommended)](#expedited-build-recommended)
    - [Full Build](#full-build)
- [Verification](#verification)
- [Usage](#usage)
- [Troubleshooting](#troubleshooting)
- [References](#references)
- [License](#license)
- [Code Authorship](#code-authorship)
- [Contact](#contact)

## Prerequisites

### Required Software

- **Docker** (v20.10+) — [Installation guide](https://www.docker.com/get-started/)

### Required HSL Libraries

The HSL packages the solver uses are obtained as a single **libHSL** source
snapshot from the [STFC licensing portal](https://licences.stfc.ac.uk)
(available at no cost for academic use). Download the **source code tarball**
(`libHSL.v<version>.tar.gz`) — the prebuilt binary tarballs are not usable
by this build.

The build stages the following packages out of the snapshot:

| Package | Role |
|---------|------|
| MA48 | sparse LU |
| MA51 | block rank detection |
| MA60 | solve-quality diagnostics (`-condest`) |
| MC71 | condition estimation (with MA60) |
| HSL_MC66 | bordered ordering |
| HSL_MC79 | structural probe (`-solmed probe`) |
| HSL_MP48 | parallel bordered solve (SBBD) |

**Note**: the staged decks are checksum-pinned to the verified snapshot
(libHSL 2026.8.4); the build fails loudly if a newer snapshot changes any
package the solver uses, so that the change can be reviewed first.

### Linux-Specific Requirements

Linux users **must** configure Docker to run without `sudo`:
- Follow the [post-installation steps](https://docs.docker.com/engine/install/linux-postinstall/)

## Quick Start

```bash
# Clone repository
git clone --depth 1 https://github.com/matthewcantele/teems-solver.git
cd teems-solver
# Copy the libHSL source tarball into the empty hsl/ directory
```

```bash
# Expedited build (~5 minutes)
docker build --pull \
  -t teems:latest \
  --build-arg PATH_LIBHSL="hsl/libHSL.v2026.8.4.tar.gz" \
  -f ./docker/expedited_build/Dockerfile \
  .
```

## Installation

### Docker Setup

Install Docker following the official documentation:
- **Windows/Mac**: [Docker Desktop](https://www.docker.com/get-started/)
- **Linux**: [Docker Engine](https://docs.docker.com/engine/install/)

After installation, verify Docker is working:
```bash
docker --version
docker run hello-world
```

### HSL Libraries

1. Visit the [STFC licensing portal](https://licences.stfc.ac.uk) and request **libHSL**
2. Download the **source code tarball** (`libHSL.v<version>.tar.gz`)
3. Keep the downloaded `.tar.gz` file in its original format

### Build Options

Two build approaches are available:

| Build Type | Time | Description |
|------------|------|-------------|
| **Expedited** | ~5 min | Recommended for most users. Pulls a pre-built base image containing all open-source dependencies; only HSL compilation and the final solver link are performed locally. |
| **Full** | ~40 min | Compiles all dependencies from source. Use if you cannot pull from Dockerhub or require customization. |

#### Step 1: Clone Repository

```bash
git clone --depth 1 https://github.com/matthewcantele/teems-solver.git
cd teems-solver
```

#### Step 2: Prepare HSL Libraries

Copy the libHSL source tarball into the `hsl/` directory.

The `hsl/` directory should contain:
```
hsl/
└── libHSL.v2026.8.4.tar.gz
```

(Individual per-package HSL tarballs are no longer a supported input;
the build stops with a pointer to these instructions if given one.)

#### Step 3: Build Docker Image

##### Expedited Build (Recommended)

Pulls the pre-built base image (`matthewcantele/teems_base:latest`) from Dockerhub and adds HSL compilation and the final solver binary. The correct base image architecture (`amd64` or `arm64`) is selected automatically.

```bash
docker build --pull \
  -t teems:latest \
  --build-arg PATH_LIBHSL="hsl/libHSL.v2026.8.4.tar.gz" \
  -f ./docker/expedited_build/Dockerfile \
  .
```

**Build time**: ~5 minutes

To build against a specific published base (for example a pre-release tag), add `--build-arg BASE_IMAGE=matthewcantele/teems_base:<tag>`; the local image tag (`-t teems:<tag>`) is what `ems_option_set(docker_tag = "<tag>")` selects in the R package.

**ISA level.** The published `amd64` base is compiled at `x86-64-v2` (runs on every Intel/AMD CPU since ~2009) and the `arm64` base at `armv8-a` (all Apple M-series and Graviton); the expedited build compiles the HSL libraries and the solver at the same level, read from the base's `/opt/teems-solver/archflags`. Higher levels (`x86-64-v3` etc.) showed no measurable win — the dense kernels dispatch per CPU at run time through OpenBLAS — so only the portable levels are published. The full build accepts `--build-arg MARCH=<level>` if you want to compile everything for your own machine.

##### Full Build

Compiles all dependencies (MPICH, PETSc, HSL, solver) from source on a `debian:bookworm-slim` base:

```bash
docker build -t teems:latest \
  --build-arg PATH_LIBHSL="hsl/libHSL.v2026.8.4.tar.gz" \
  -f ./docker/full_build/Dockerfile \
  .
```

**Build time**: ~40 minutes (varies by system)

## Verification

After building, verify the image exists:

```bash
docker image ls | grep teems
```

Expected output:
```
REPOSITORY   TAG       IMAGE ID       CREATED          SIZE
teems        latest    abc123def456   2 minutes ago    ~175MB
```

## Usage

The TEEMS solver is designed to work seamlessly with the [teems R package](https://github.com/teemsphere/teems-R). Once the Docker image is built, the R package handles all solver invocation automatically via `ems_solve()`.

## Troubleshooting

### Docker Permission Denied (Linux)

**Error**: `permission denied while trying to connect to the Docker daemon socket`

**Solution**: Add your user to the docker group:
```bash
sudo usermod -aG docker $USER
newgrp docker
```

### Build Fails: HSL Libraries Not Found

**Error**: `COPY failed: file not found`

**Solution**: Verify HSL tarballs are in the `hsl/` directory with the exact filenames passed to `--build-arg`:
```bash
ls -l hsl/
```

### Getting Help

If you encounter issues not covered here:

1. Check [existing issues](https://github.com/matthewcantele/teems-solver/issues)
2. [Open a new issue](https://github.com/matthewcantele/teems-solver/issues/new) with system information (`uname -a`, `docker --version`, the solver version from `docker run --rm teems:<tag> /opt/teems-solver/solver/teems-solver -version`) and the complete error message

## References

The solver implements methods described in the following works.

**Tablo language and solution procedure**

- Horridge J.M., Jerie M., Mustakinov D. & Schiffmann F. (2018), *GEMPACK manual*, GEMPACK Software, Centre of Policy Studies, Victoria University, Melbourne, ISBN 978-1-921654-34-3. <https://ideas.repec.org/p/cop/wpaper/gpman.html>
- Pearson K.R. (1991), 'Solving Nonlinear Economic Models Accurately via a Linear Representation', Impact Preliminary Working Paper No. IP-55, Melbourne (July), pp.39. <http://www.copsmodels.com/elecpapr/ip-55.htm>

**Complementarities and subtotals**

- Harrison, W.J., Mark Horridge, K.R. Pearson and Glyn Wittwer (2002), 'A Practical Method for Explicitly Modeling Quotas and Other Complementarities', Computational Economics, June 2004, Vol. 23(4), pp. 325-341. [A preliminary version was Centre of Policy Studies and the Impact Project Preliminary Working Paper No. IP-78, Melbourne (April), pp.19.] <http://www.copsmodels.com/elecpapr/ip-78.htm>
- Harrison, W.J., J.M. Horridge and K.R. Pearson (2000), 'Decomposing Simulation Results with Respect to Exogenous Shocks', Computational Economics, vol.15, pp.227-249. [A preliminary version was Centre of Policy Studies and Impact Project Preliminary Working Paper No. IP-73, May 1999.] <http://www.copsmodels.com/elecpapr/ip-73.htm>

**Runge-Kutta methods**

- Schiffmann, F. (2022), 'Runge Kutta integrators for fast and accurate solutions in GEMPACK'.
- Munthe-Kaas, H. (1999), 'High order Runge-Kutta methods on manifolds', Applied Numerical Mathematics, vol. 29, pp. 115-127.
- Butcher, J.C. (2015), 'Runge-Kutta Methods for Ordinary Differential Equations', in M. Al-Baali, L. Grandinetti and A. Purnama (eds), *Numerical Analysis and Optimization*, Springer Proceedings in Mathematics & Statistics, vol. 134, pp. 37-58.
- Bogacki, P. and L.F. Shampine (1989), 'A 3(2) pair of Runge-Kutta formulas', Applied Mathematics Letters, vol. 2, pp. 321-325.
- Dormand, J.R. and P.J. Prince (1980), 'A family of embedded Runge-Kutta formulae', Journal of Computational and Applied Mathematics, vol. 6, pp. 19-26.
- Gustafsson, K. (1991), 'Control theoretic techniques for stepsize selection in explicit Runge-Kutta methods', ACM Transactions on Mathematical Software, vol. 17, pp. 533-554.

**Parallel matrix methods (SBBD, NDBBD)**

- Van Ha, P. and T. Kompas (2016), 'Solving intertemporal CGE models in parallel using a singly bordered block diagonal ordering technique', Economic Modelling, vol. 52, pp. 3-12. <https://doi.org/10.1016/j.econmod.2015.07.011>
- Kompas, T. and P. Van Ha (2019), 'The "curse of dimensionality" resolved: The effects of climate change and trade barriers in large dimensional modelling', Economic Modelling, vol. 80, pp. 103-110. <https://doi.org/10.1016/j.econmod.2018.08.011>

**Sparse linear algebra**

- Duff, I.S. and J.K. Reid (1993), 'MA48, a Fortran Code for Direct Solution of Sparse Unsymmetric Linear Systems of Equations', Rutherford Appleton Laboratory Report RAL-93-072, pp.62.

**Software**

- HSL (2026), A collection of Fortran codes for large scale scientific computation. <https://www.hsl.rl.ac.uk/>
- Balay, S. et al. (2026), *PETSc/TAO Users Manual*, Argonne National Laboratory, ANL-21/39 - Revision 3.25. <https://doi.org/10.2172/2998643>
- Balay, S. et al. (2026), *PETSc Web page*. <https://petsc.org/>

## License

This project is licensed under the **GNU Affero General Public License v3.0** — see the [LICENSE](LICENSE) file for details.

| Component | License |
|-----------|---------|
| TEEMS Solver | AGPL-3.0 |
| PETSc | BSD 2-Clause |
| MPICH | Custom Open Source |
| HSL Libraries | (user must obtain) |

## Code Authorship

| Component | Author(s) |
|-----------|-----------|
| Source code (`src/`) | Tom Kompas, Ha Van Pham, Matthew Cantele |
| Docker scripts | Matthew Cantele |

## Contact

- **Project Maintainer**: [Matthew Cantele](mailto:matthew.cantele@protonmail.com)
- **Project Homepage**: [https://github.com/matthewcantele/teems-solver](https://github.com/matthewcantele/teems-solver)
- **Bug Reports**: [https://github.com/matthewcantele/teems-solver/issues](https://github.com/matthewcantele/teems-solver/issues)
