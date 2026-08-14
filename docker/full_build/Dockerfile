FROM debian:bookworm-slim AS builder

ARG TARGETARCH
# psABI/ISA level for amd64/arm64/arm (e.g. x86-64-v3); empty = portable
# default per arch. See docs/image-publishing.md.
ARG MARCH=""
ARG PATH_HSL_MA48
ARG PATH_HSL_MA51
ARG PATH_HSL_MA60
ARG PATH_HSL_MC66
ARG PATH_HSL_MC71
ARG PATH_HSL_MC79
ARG PATH_HSL_MP48
ARG BUILD_DIR=/opt/teems-solver
ARG MPICH_VERSION=5.0.1
ARG PETSC_ARCH=native
ARG PETSC_VERSION=v3.25.0
ARG TEEMS_TMP=/tmp/teems

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="${BUILD_DIR}/lib/mpi/bin:$PATH"

COPY ${PATH_HSL_MA48} ${BUILD_DIR}/hsl/ma48.tar.gz
COPY ${PATH_HSL_MA51} ${BUILD_DIR}/hsl/ma51.tar.gz
COPY ${PATH_HSL_MA60} ${BUILD_DIR}/hsl/ma60.tar.gz
COPY ${PATH_HSL_MC66} ${BUILD_DIR}/hsl/mc66.tar.gz
COPY ${PATH_HSL_MC71} ${BUILD_DIR}/hsl/mc71.tar.gz
COPY ${PATH_HSL_MC79} ${BUILD_DIR}/hsl/mc79.tar.gz
COPY ${PATH_HSL_MP48} ${BUILD_DIR}/hsl/mp48.tar.gz
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
    # HSL
    && mkdir -p ${BUILD_DIR}/LICENSES/hsl \
    && mkdir ${BUILD_DIR}/lib/ma48 \
    && tar xfz ${BUILD_DIR}/hsl/ma48.tar.gz --strip-components=1 -C ${BUILD_DIR}/lib/ma48 \
    && cp ${BUILD_DIR}/lib/ma48/LICENCE ${BUILD_DIR}/LICENSES/hsl \
    && mkdir ${BUILD_DIR}/lib/ma51 \
    && tar xfz ${BUILD_DIR}/hsl/ma51.tar.gz --strip-components=1 -C ${BUILD_DIR}/lib/ma51 \
    && mkdir ${BUILD_DIR}/lib/ma60 \
    && tar xfz ${BUILD_DIR}/hsl/ma60.tar.gz --strip-components=1 -C ${BUILD_DIR}/lib/ma60 \
    && mkdir ${BUILD_DIR}/lib/mc71 \
    && tar xfz ${BUILD_DIR}/hsl/mc71.tar.gz --strip-components=1 -C ${BUILD_DIR}/lib/mc71 \
    && mkdir ${BUILD_DIR}/lib/mc66 \
    && tar xfz ${BUILD_DIR}/hsl/mc66.tar.gz --strip-components=1 -C ${BUILD_DIR}/lib/mc66 \
    && mkdir ${BUILD_DIR}/lib/mc79 \
    && tar xfz ${BUILD_DIR}/hsl/mc79.tar.gz --strip-components=1 -C ${BUILD_DIR}/lib/mc79 \
    && mkdir ${BUILD_DIR}/lib/mp48 \
    && tar xfz ${BUILD_DIR}/hsl/mp48.tar.gz --strip-components=1 -C ${BUILD_DIR}/lib/mp48 \
    && cd ${BUILD_DIR}/lib/ma48 \
    && autoreconf -fi \
    && ./configure CFLAGS="-O2 $ARCH_FLAGS" FFLAGS="-O2 $ARCH_FLAGS" FCFLAGS="-O2 $ARCH_FLAGS" \
    && make \
    && cd ${BUILD_DIR}/lib/ma51 \
    && autoreconf -fi \
    && ./configure CFLAGS="-O2 $ARCH_FLAGS" FFLAGS="-O2 $ARCH_FLAGS" FCFLAGS="-O2 $ARCH_FLAGS" \
    && make \
    && cd ${BUILD_DIR}/lib/mp48 \
    && autoreconf -fi \
    && ./configure FFLAGS="-fallow-argument-mismatch -O2 $ARCH_FLAGS" FCFLAGS="-fallow-argument-mismatch -O2 $ARCH_FLAGS" \
    && make \
    && cp ${BUILD_DIR}/lib/mp48/src/hsl_mp48d.f90 ${BUILD_DIR}/src \
    && cp ${BUILD_DIR}/lib/mp48/src/ddeps.f ${BUILD_DIR}/src \
    && cp ${BUILD_DIR}/lib/mp48/src/hsl_mp01.mod ${BUILD_DIR}/src \
    && cd ${BUILD_DIR}/src \
    && chmod +x mp48_mod.sh \
    && ./mp48_mod.sh \
    && cp ${BUILD_DIR}/lib/mc66/src/hsl_mc66d.f90 ${BUILD_DIR}/src \
    && cp ${BUILD_DIR}/lib/mc66/src/ddeps90.f90 ${BUILD_DIR}/src \
    && cp ${BUILD_DIR}/lib/mc79/src/hsl_mc79i.f90 ${BUILD_DIR}/src \
    && cp ${BUILD_DIR}/lib/mc79/src/hsl_mc79i_ciface.f90 ${BUILD_DIR}/src \
    && cp ${BUILD_DIR}/lib/mc79/include/hsl_mc79i.h ${BUILD_DIR}/src \
    # -condest diagnostics (MA60/MC71): sources compile into the solver;
    # FD15AD is extracted alone -- ma60's ddeps.f bundles duplicate
    # MA48/MC59/MC71 and must never be compiled whole
    && cp ${BUILD_DIR}/lib/ma60/src/ma60d.f ${BUILD_DIR}/src \
    && cp ${BUILD_DIR}/lib/mc71/src/mc71d.f ${BUILD_DIR}/src \
    && awk '/DOUBLE PRECISION FUNCTION FD15AD/{f=1} f&&/^      SUBROUTINE/{exit} f' \
       ${BUILD_DIR}/lib/ma60/src/ddeps.f > ${BUILD_DIR}/src/fd15d.f \
    # final build — solver objects at the same ISA level as the libraries
    && cd ${BUILD_DIR}/src \
    && make OPT="-Ofast $ARCH_FLAGS" teems-solver teems-solver-f64 \
    && mv ${BUILD_DIR}/src ${BUILD_DIR}/solver \
    && rm -rf ${TEEMS_TMP} ${BUILD_DIR}/hsl \
    && mkdir -p ${BUILD_DIR}/mpi-runtime/lib ${BUILD_DIR}/mpi-runtime/bin \
    && find ${BUILD_DIR}/lib/mpi/lib -name "*.so*" \
           -exec cp -P {} ${BUILD_DIR}/mpi-runtime/lib/ \; \
    && for b in mpirun mpiexec mpiexec.hydra hydra_pmi_proxy hydra_nameserver hydra_persist; do \
           [ -f ${BUILD_DIR}/lib/mpi/bin/$b ] \
           && cp ${BUILD_DIR}/lib/mpi/bin/$b ${BUILD_DIR}/mpi-runtime/bin/; \
       done; true

# ── runtime image ──────────────────────────────────────────────────────────────
FROM debian:bookworm-slim

LABEL org.opencontainers.image.name="TEEMS solver Dockerfile (Debian bookworm-slim)"
LABEL org.opencontainers.image.version="1.0.0"
LABEL org.opencontainers.image.description="Debian bookworm-slim build for TEEMS solver"
LABEL org.opencontainers.image.author="Matthew Cantele <matthew.cantele@protonmail.com>"
LABEL org.opencontainers.image.maintainer="Matthew Cantele <matthew.cantele@protonmail.com>"
LABEL org.opencontainers.image.source="https://github.com/matthewcantele/teems-solver"
LABEL org.opencontainers.image.licenses="AGPL-3.0"

ARG BUILD_DIR=/opt/teems-solver
ARG USERNAME=teems
ARG USER_UID=1000
ARG USER_GID=$USER_UID

ARG DEBIAN_FRONTEND=noninteractive
ENV PATH="${BUILD_DIR}/lib/mpi/bin:$PATH"
# One BLAS thread per MPI rank: the pthread OpenBLAS otherwise opens a
# full thread pool per rank and oversubscribes (measured ~20% wall cost
# on SBBD-2 at 4.4M equations; golden outputs bit-identical either way)
ENV OPENBLAS_NUM_THREADS=1

# solver binary and LICENSES
COPY --from=builder ${BUILD_DIR}/solver/teems-solver     ${BUILD_DIR}/solver/teems-solver
COPY --from=builder ${BUILD_DIR}/solver/teems-solver-f64 ${BUILD_DIR}/solver/teems-solver-f64
COPY --from=builder ${BUILD_DIR}/LICENSES            ${BUILD_DIR}/LICENSES

# HSL shared libraries (MA48 and MA51 are linked dynamically)
COPY --from=builder \
    ${BUILD_DIR}/lib/ma48/src/.libs/libma48.so.0.0.0 \
    /usr/local/lib/libma48.so.0.0.0
COPY --from=builder \
    ${BUILD_DIR}/lib/ma51/src/.libs/libma51.so.0.0.0 \
    /usr/local/lib/libma51.so.0.0.0

# MPICH runtime — shared libs and launcher only
COPY --from=builder ${BUILD_DIR}/mpi-runtime/lib  ${BUILD_DIR}/lib/mpi/lib
COPY --from=builder ${BUILD_DIR}/mpi-runtime/bin  ${BUILD_DIR}/lib/mpi/bin

RUN echo 'APT::Install-Suggests "0";' >> /etc/apt/apt.conf.d/00-docker \
    && echo 'APT::Install-Recommends "0";' >> /etc/apt/apt.conf.d/00-docker \
    && apt-get update \
    && apt-get install -y \
        libopenblas0 \
        libgfortran5 \
        libgomp1 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && echo "${BUILD_DIR}/lib/mpi/lib" > /etc/ld.so.conf.d/teems-mpi.conf \
    && ldconfig \
    # transition symlink for scripts referencing the old binary name;
    # remove one release after the rename ships
    && ln -s teems-solver ${BUILD_DIR}/solver/hsl \
    && groupadd -g $USER_GID $USERNAME \
    && useradd -m -u $USER_UID -g $USER_GID -s /bin/bash $USERNAME

USER $USERNAME
WORKDIR /opt/teems
CMD ["/bin/bash"]
