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
- [License](#license)
- [Code Authorship](#code-authorship)
- [Contact](#contact)

## Prerequisites

### Required Software

- **Docker** (v20.10+) — [Installation guide](https://www.docker.com/get-started/)

### Required HSL Libraries

The following HSL libraries must be obtained directly from [HSL](https://www.hsl.rl.ac.uk/). These are available at no cost for academic use.

| Library | Version | Role |
|---------|---------|------|
| MA48 | 2.2.0 | sparse LU |
| MA51 | 1.0.0 | block rank detection |
| MA60 | 1.2.0 | solve-quality diagnostics (`-condest`) |
| MC71 | 1.0.0 | condition estimation (with MA60) |
| HSL_MC66 | 2.2.1 | bordered ordering |
| HSL_MC79 | 1.1.1 | structural probe (`-solmed probe`) |
| HSL_MP48 | 2.1.1 | parallel bordered solve (SBBD) |

**Note**: Backward compatibility with other HSL library versions is not guaranteed.

### Linux-Specific Requirements

Linux users **must** configure Docker to run without `sudo`:
- Follow the [post-installation steps](https://docs.docker.com/engine/install/linux-postinstall/)

## Quick Start

```bash
# Clone repository
git clone --depth 1 https://github.com/matthewcantele/teems-solver.git
cd teems-solver
# Copy HSL tarballs into the empty hsl/ directory
```

```bash
# Expedited build (~5 minutes)
docker build --pull \
  -t teems:latest \
  --build-arg PATH_HSL_MA48="hsl/ma48-2.2.0.tar.gz" \
  --build-arg PATH_HSL_MA51="hsl/ma51-1.0.0.tar.gz" \
  --build-arg PATH_HSL_MA60="hsl/ma60-1.2.0.tar.gz" \
  --build-arg PATH_HSL_MC66="hsl/hsl_mc66-2.2.1.tar.gz" \
  --build-arg PATH_HSL_MC71="hsl/mc71-1.0.0.tar.gz" \
  --build-arg PATH_HSL_MC79="hsl/hsl_mc79-1.1.1.tar.gz" \
  --build-arg PATH_HSL_MP48="hsl/hsl_mp48-2.1.1.tar.gz" \
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

1. Visit [HSL](https://www.hsl.rl.ac.uk/) and request the required libraries
2. Download the specific versions listed in [Prerequisites](#required-hsl-libraries)
3. Keep the downloaded `.tar.gz` files in their original format

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

Copy your HSL library tarballs into the `hsl/` directory.

The `hsl/` directory should contain:
```
hsl/
├── ma48-2.2.0.tar.gz
├── ma51-1.0.0.tar.gz
├── ma60-1.2.0.tar.gz
├── mc71-1.0.0.tar.gz
├── hsl_mc66-2.2.1.tar.gz
├── hsl_mc79-1.1.1.tar.gz
└── hsl_mp48-2.1.1.tar.gz
```

#### Step 3: Build Docker Image

##### Expedited Build (Recommended)

Pulls the pre-built base image (`matthewcantele/teems_base:latest`) from Dockerhub and adds HSL compilation and the final solver binary. The correct base image architecture (`amd64` or `arm64`) is selected automatically.

```bash
docker build --pull \
  -t teems:latest \
  --build-arg PATH_HSL_MA48="hsl/ma48-2.2.0.tar.gz" \
  --build-arg PATH_HSL_MA51="hsl/ma51-1.0.0.tar.gz" \
  --build-arg PATH_HSL_MA60="hsl/ma60-1.2.0.tar.gz" \
  --build-arg PATH_HSL_MC66="hsl/hsl_mc66-2.2.1.tar.gz" \
  --build-arg PATH_HSL_MC71="hsl/mc71-1.0.0.tar.gz" \
  --build-arg PATH_HSL_MC79="hsl/hsl_mc79-1.1.1.tar.gz" \
  --build-arg PATH_HSL_MP48="hsl/hsl_mp48-2.1.1.tar.gz" \
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
  --build-arg PATH_HSL_MA48="hsl/ma48-2.2.0.tar.gz" \
  --build-arg PATH_HSL_MA51="hsl/ma51-1.0.0.tar.gz" \
  --build-arg PATH_HSL_MA60="hsl/ma60-1.2.0.tar.gz" \
  --build-arg PATH_HSL_MC66="hsl/hsl_mc66-2.2.1.tar.gz" \
  --build-arg PATH_HSL_MC71="hsl/mc71-1.0.0.tar.gz" \
  --build-arg PATH_HSL_MC79="hsl/hsl_mc79-1.1.1.tar.gz" \
  --build-arg PATH_HSL_MP48="hsl/hsl_mp48-2.1.1.tar.gz" \
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
2. [Open a new issue](https://github.com/matthewcantele/teems-solver/issues/new) with system information (`uname -a`, `docker --version`) and the complete error message

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
