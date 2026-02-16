# TEEMS Solver

[![License](https://img.shields.io/badge/License-GPL_v3.0-blue.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-0.9-green.svg)](https://github.com/matthewcantele/teems-solver/releases)
[![Docker](https://img.shields.io/badge/docker-required-2496ED.svg)](https://www.docker.com/)

High-performance optimization solver for the TEEMS (Temporally Explicit Ecosystem Management Simulator) R package. Built on PETSc with MPI parallelization and HSL sparse linear algebra libraries.

## Overview

The TEEMS solver is a specialized computational engine designed for large-scale ecosystem management optimization problems. It combines:

- **PETSc** for portable, scalable scientific computation
- **MPI** (OpenMPI/MPICH) for parallel processing
- **HSL libraries** for efficient sparse matrix operations

The solver is distributed as Docker images to ensure reproducibility and ease of deployment across different systems.

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
- [Documentation](#documentation)
- [License](#license)
- [Code Authorship](#code-authorship)
- [Contact](#contact)

## Prerequisites

Before beginning installation, ensure you have:

### Required Software

- **Docker** (v20.10+) - [Installation guide](https://www.docker.com/get-started/)
- **Git** - For cloning the repository

### Required HSL Libraries

The following HSL libraries must be obtained directly from [HSL](https://www.hsl.rl.ac.uk/):

| Library | Version | Download Link |
|---------|---------|---------------|
| MA48 | 2.2.0 | [hsl.rl.ac.uk/catalogue/ma48.html](https://www.hsl.rl.ac.uk/catalogue/ma48.html) |
| MA51 | 1.0.0 | [hsl.rl.ac.uk/catalogue/ma51.html](https://www.hsl.rl.ac.uk/catalogue/ma51.html) |
| HSL_MC66 | 2.2.1 | [hsl.rl.ac.uk/catalogue/hsl_mc66.html](https://www.hsl.rl.ac.uk/catalogue/hsl_mc66.html) |
| HSL_MP48 | 2.1.1 | [hsl.rl.ac.uk/catalogue/hsl_mp48.html](https://www.hsl.rl.ac.uk/catalogue/hsl_mp48.html) |

**Note**: Backward compatibility with other HSL library versions is not guaranteed.

### Linux-Specific Requirements

Linux users **must** configure Docker to run without `sudo`:
- Follow the [post-installation steps](https://docs.docker.com/engine/install/linux-postinstall/)

## Quick Start

```bash
# Clone repository
git clone --depth 1 https://github.com/matthewcantele/teems-solver.git
cd teems-solver

# Copy HSL libraries to hsl/ directory
cp /path/to/your/hsl/*.tar.gz hsl/

# Build (expedited version - 5 minutes)
docker build -t teems:latest \
  --build-arg PATH_HSL_MA48="hsl/ma48-2.2.0.tar.gz" \
  --build-arg PATH_HSL_MA51="hsl/ma51-1.0.0.tar.gz" \
  --build-arg PATH_HSL_MC66="hsl/hsl_mc66-2.2.1.tar.gz" \
  --build-arg PATH_HSL_MP48="hsl/hsl_mp48-2.1.1.tar.gz" \
  -f ./docker/expedited_build/Dockerfile \
  .

# Verify installation
docker image ls | grep teems
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

| Build Type | Time | Use Case |
|------------|------|----------|
| **Expedited** | ~5 min | Recommended for most users |
| **Full** | ~40 min | Developers, customization needed |

The expedited build uses a pre-built base image containing all open-source dependencies, requiring only HSL library integration.

#### Step 1: Clone Repository

```bash
git clone --depth 1 https://github.com/matthewcantele/teems-solver.git
cd teems-solver
```

#### Step 2: Prepare HSL Libraries

Copy your HSL library files into the `hsl/` directory:

```bash
cp /path/to/downloads/*.tar.gz hsl/
```

The `hsl/` directory should contain:
```
hsl/
├── ma48-2.2.0.tar.gz
├── ma51-1.0.0.tar.gz
├── hsl_mc66-2.2.1.tar.gz
└── hsl_mp48-2.1.1.tar.gz
```

#### Step 3: Build Docker Image

##### Expedited Build (Recommended)

Uses pre-built base image with all open-source dependencies:

```bash
docker build -t teems:latest \
  --build-arg PATH_HSL_MA48="hsl/ma48-2.2.0.tar.gz" \
  --build-arg PATH_HSL_MA51="hsl/ma51-1.0.0.tar.gz" \
  --build-arg PATH_HSL_MC66="hsl/hsl_mc66-2.2.1.tar.gz" \
  --build-arg PATH_HSL_MP48="hsl/hsl_mp48-2.1.1.tar.gz" \
  -f ./docker/expedited_build/Dockerfile \
  .
```

**Build time**: ~5 minutes

##### Full Build

Compiles all dependencies from source:

```bash
docker build -t teems:latest \
  --build-arg PATH_HSL_MA48="hsl/ma48-2.2.0.tar.gz" \
  --build-arg PATH_HSL_MA51="hsl/ma51-1.0.0.tar.gz" \
  --build-arg PATH_HSL_MC66="hsl/hsl_mc66-2.2.1.tar.gz" \
  --build-arg PATH_HSL_MP48="hsl/hsl_mp48-2.1.1.tar.gz" \
  -f ./docker/full_build/Dockerfile \
  .
```

**Build time**: ~40 minutes (varies by system)

## Verification

After building, verify the image exists:

```bash
docker image ls
```

Expected output:
```
REPOSITORY   TAG       IMAGE ID       CREATED          SIZE
teems        latest    abc123def456   2 minutes ago    ~2GB
```

## Usage

The TEEMS solver is designed to work seamlessly with the TEEMS R package:

### With TEEMS R Package (Recommended)

```r
library(teems)

# The R package automatically interfaces with the Docker container
result <- teems_solve(model_data)
```

### In-Situ Solve

The TEEMS R package also supports in-container solving for enhanced performance:

## Troubleshooting

### Common Issues

#### Docker Permission Denied (Linux)

**Error**: `permission denied while trying to connect to the Docker daemon socket`

**Solution**: Add your user to the docker group:
```bash
sudo usermod -aG docker $USER
newgrp docker
```

#### Build Fails: HSL Libraries Not Found

**Error**: `COPY failed: file not found`

**Solution**: Verify HSL files are in the `hsl/` directory with exact filenames:
```bash
ls -l hsl/
```

### Getting Help

If you encounter issues not covered here:

1. Check [existing issues](https://github.com/matthewcantele/teems-solver/issues)
2. Review Docker logs: `docker logs <container_id>`
3. [Open a new issue](https://github.com/matthewcantele/teems-solver/issues/new) with:
   - System information (`uname -a`, `docker --version`)
   - Complete error message
   - Steps to reproduce

## Documentation

Additional documentation is available:

- [Base Layer Dockerfile](docker/base_layer/BASE_LAYER.md) - Details on the pre-built base image
- [Full Build Dockerfile](docker/full_build/README.md) - Complete build from source documentation
- [Expedited Build Dockerfile](docker/expedited_build/README.md) - Quick build documentation

## License

This project is licensed under the **GNU General Public License v3.0** - see the [LICENSE](LICENSE) file for details.

### Component Licenses

Individual components retain their original licenses:

| Component | License |
|-----------|---------|
| TEEMS Solver | GPL-3.0 |
| PETSc | BSD 2-Clause |
| OpenMPI | BSD 3-Clause |
| MPICH | Custom Open Source |
| HSL Libraries | (user must obtain) |

## Code Authorship

This work represents years of collaborative development:

| Component | Author(s) |
|-----------|-----------|
| C source code (`src/`) | Tom Kompas, Ha Van Pham |
| Docker scripts | Matthew Cantele |
| Build system | Matthew Cantele |

## Contact

- **Project Maintainer**: [Matthew Cantele](mailto:matthew.cantele@protonmail.com)
- **Project Homepage**: [https://github.com/matthewcantele/teems-solver](https://github.com/matthewcantele/teems-solver)
- **Bug Reports**: [https://github.com/matthewcantele/teems-solver/issues](https://github.com/matthewcantele/teems-solver/issues)

**Note**: This solver is actively maintained. For the latest updates and releases, watch this repository or check the [releases page](https://github.com/matthewcantele/teems-solver/releases).