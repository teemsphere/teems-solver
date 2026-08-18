FROM debian:bookworm-slim

LABEL org.opencontainers.image.name="TEEMS solver base layer Dockerfile"
LABEL org.opencontainers.image.version="2.0.0-dev.1"
LABEL org.opencontainers.image.description="Base layer Docker build for TEEMS solver (open-source components)"
LABEL org.opencontainers.image.author="Matthew Cantele <matthew.cantele@protonmail.com>"
LABEL org.opencontainers.image.maintainer="Matthew Cantele <matthew.cantele@protonmail.com>"
LABEL org.opencontainers.image.source="https://github.com/matthewcantele/teems-solver"
LABEL org.opencontainers.image.licenses="AGPL-3.0"

ARG TARGETARCH
# psABI/ISA level for amd64/arm64/arm (e.g. x86-64-v3); empty = portable
# default per arch. The resolved flags are recorded in
# ${BUILD_DIR}/archflags so downstream layers compile with the same level.
ARG MARCH=""
ARG BUILD_DIR=/opt/teems-solver
ARG MPICH_VERSION=5.0.1
ARG PETSC_ARCH=native
ARG PETSC_VERSION=v3.25.0
ARG TEEMS_TMP=/tmp/teems

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="${BUILD_DIR}/lib/mpi/bin:$PATH"

COPY ./src ${BUILD_DIR}/src

RUN echo 'APT::Install-Suggests "0";' >> /etc/apt/apt.conf.d/00-docker \
    && echo 'APT::Install-Recommends "0";' >> /etc/apt/apt.conf.d/00-docker \
    && apt-get update \
    && apt-get install -y \
        automake \
        autoconf \
        libtool \
        build-essential \
        gfortran \
        wget \
        ca-certificates \
        libopenblas-dev \
        git \
        python3 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && ARCH="${TARGETARCH:-}" \
    && if [ -z "$ARCH" ]; then \
         case "$(uname -m)" in \
           x86_64)        ARCH=amd64   ;; \
           aarch64|arm64) ARCH=arm64   ;; \
           armv7l)        ARCH=arm     ;; \
           ppc64le)       ARCH=ppc64le ;; \
           s390x)         ARCH=s390x   ;; \
           *)             ARCH=unknown ;; \
         esac; \
       fi \
    && case "$ARCH" in \
         amd64)   ARCH_FLAGS="-march=${MARCH:-x86-64-v2}"          ;; \
         arm64)   ARCH_FLAGS="-march=${MARCH:-armv8-a}"            ;; \
         arm)     ARCH_FLAGS="-march=${MARCH:-armv7-a} -mfpu=neon" ;; \
         ppc64le) ARCH_FLAGS="-mcpu=powerpc64le"                   ;; \
         *)       ARCH_FLAGS=""                                    ;; \
       esac \
    && if [ "$ARCH" = "unknown" ]; then OPT_FLAGS="-O2"; else OPT_FLAGS="-Ofast $ARCH_FLAGS"; fi \
    && printf "%s" "$ARCH_FLAGS" > ${BUILD_DIR}/archflags \
    && BLAS_LIB=$(find /usr/lib/$(gcc -print-multiarch) -name "libopenblas.so" | head -1) \
    # MPICH
    && mkdir -p ${TEEMS_TMP} \
    && cd ${TEEMS_TMP} \
    && mkdir mpich \
    && cd mpich \
    && wget https://www.mpich.org/static/downloads/${MPICH_VERSION}/mpich-${MPICH_VERSION}.tar.gz \
    && tar xfz mpich-${MPICH_VERSION}.tar.gz --strip-components=1 \
    # production MPICH: -O3 + NDEBUG internally (install guide "Compiler
    # Optimization Levels"); arch flags go in MPICHLIB_* so the mpicc/mpif90
    # wrappers stay clean — downstream compiles set their own levels
    && ./configure --prefix=${BUILD_DIR}/lib/mpi \
        --enable-fast=all,O3 \
        MPICHLIB_CFLAGS="$ARCH_FLAGS" \
        MPICHLIB_CXXFLAGS="$ARCH_FLAGS" \
        MPICHLIB_FFLAGS="$ARCH_FLAGS" \
        MPICHLIB_FCFLAGS="$ARCH_FLAGS" \
    && make \
    && make install \
    # PETSc
    && cd ${BUILD_DIR}/lib \
    && git clone -b ${PETSC_VERSION} https://gitlab.com/petsc/petsc.git petsc \
    && cd petsc \
    && mkdir -p ${BUILD_DIR}/LICENSES/petsc \
    && cp ${BUILD_DIR}/lib/petsc/LICENSE ${BUILD_DIR}/LICENSES/petsc \
    # memalign=64: SIMD-friendly array alignment (AVX2/AVX-512, SVE);
    # cxx/fortran bindings off: the solver calls PETSc from C only
    && python3 ./config/configure.py \
        PETSC_ARCH=${PETSC_ARCH} \
        PETSC_DIR=$PWD \
        --with-debugging=0 \
        --with-mpi-dir=${BUILD_DIR}/lib/mpi \
        COPTFLAGS="$OPT_FLAGS" \
        CXXOPTFLAGS="$OPT_FLAGS" \
        FOPTFLAGS="$OPT_FLAGS" \
        --with-memalign=64 \
        --with-cxx=0 \
        --with-fortran-bindings=0 \
        --with-blas-lib="$BLAS_LIB" \
        --with-lapack-lib="$BLAS_LIB" \
        --with-shared-libraries=0 \
        --with-openmp \
        --with-threadsafety \
        --with-log=0 \
        --with-x=0 \
        --force \
    && make PETSC_DIR=${PWD} PETSC_ARCH=${PETSC_ARCH} all \
    && rm -rf ${TEEMS_TMP}

CMD ["/bin/bash"]
