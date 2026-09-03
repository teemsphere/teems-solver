# Select the base variant matching the target ISA level (see
# README, ISA levels); the solver and HSL libraries below are
# compiled at the level recorded in the base image's archflags file.
ARG BASE_IMAGE=matthewcantele/teems_base:latest
FROM ${BASE_IMAGE} AS builder

ARG PATH_LIBHSL
ARG BUILD_DIR=/opt/teems-solver

COPY ${PATH_LIBHSL} ${BUILD_DIR}/hsl/libHSL.tar.gz
COPY docker/stage_hsl.sh ${BUILD_DIR}/hsl/stage_hsl.sh
# Overwrite the src snapshot baked into teems_base with the current tree;
# without this the image silently ships stale solver code.
COPY ./src ${BUILD_DIR}/src

RUN ARCH_FLAGS="$(cat ${BUILD_DIR}/archflags 2>/dev/null || true)" \
    # HSL — staged from a single libHSL source snapshot; stage_hsl.sh
    # reconstructs the historical per-package layout from checksum-pinned
    # decks (see docker/stage_hsl.sh) and installs LICENSES/hsl
    && bash ${BUILD_DIR}/hsl/stage_hsl.sh ${BUILD_DIR}/hsl/libHSL.tar.gz ${BUILD_DIR} \
    # MA48/MA51 shared libraries (previously autotools/libtool; same
    # sources, same flags, same soname and .libs layout — double-precision
    # decks only: the solver references no single-precision entries)
    && mkdir -p ${BUILD_DIR}/lib/ma48/src/.libs ${BUILD_DIR}/lib/ma51/src/.libs \
    && cd ${BUILD_DIR}/lib/ma48/src \
    && gfortran -O2 $ARCH_FLAGS -fPIC -shared -Wl,-soname,libma48.so.0 \
           -o .libs/libma48.so.0.0.0 ma48d.f ddeps.f \
    && ln -sf libma48.so.0.0.0 .libs/libma48.so \
    && cd ${BUILD_DIR}/lib/ma51/src \
    && gfortran -O2 $ARCH_FLAGS -fPIC -shared -Wl,-soname,libma51.so.0 \
           -o .libs/libma51.so.0.0.0 ma51d.f ddeps.f \
    && ln -sf libma51.so.0.0.0 .libs/libma51.so \
    # hsl_mp01 module (previously a product of the MP48 autotools build)
    && cd ${BUILD_DIR}/lib/mp48/src \
    && ${BUILD_DIR}/lib/mpi/bin/mpif90 -fallow-argument-mismatch -O2 $ARCH_FLAGS -c common90.f90 \
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
    # -condest diagnostics (MA60/MC71/FD15): sources compile into the
    # solver; FD15 comes as its own deck from libHSL
    && cp ${BUILD_DIR}/lib/ma60/src/ma60d.f ${BUILD_DIR}/src \
    && cp ${BUILD_DIR}/lib/mc71/src/mc71d.f ${BUILD_DIR}/src \
    && cp ${BUILD_DIR}/lib/fd15/fd15d.f ${BUILD_DIR}/src \
    # final build — solver objects at the same ISA level as the base libs
    && cd ${BUILD_DIR}/src \
    && make OPT="-Ofast $ARCH_FLAGS" teems-solver teems-solver-f64 \
    && mv ${BUILD_DIR}/src ${BUILD_DIR}/solver \
    && rm -rf ${BUILD_DIR}/hsl \
    && mkdir -p ${BUILD_DIR}/mpi-runtime/lib ${BUILD_DIR}/mpi-runtime/bin \
    && find ${BUILD_DIR}/lib/mpi/lib -name "*.so*" \
           -exec cp -P {} ${BUILD_DIR}/mpi-runtime/lib/ \; \
    && for b in mpirun mpiexec mpiexec.hydra hydra_pmi_proxy hydra_nameserver hydra_persist; do \
           [ -f ${BUILD_DIR}/lib/mpi/bin/$b ] \
           && cp ${BUILD_DIR}/lib/mpi/bin/$b ${BUILD_DIR}/mpi-runtime/bin/; \
       done; true

# ── runtime image ──────────────────────────────────────────────────────────────
FROM debian:bookworm-slim

LABEL org.opencontainers.image.name="TEEMS solver expedited build Dockerfile"
LABEL org.opencontainers.image.version="1.1.0-dev.3"
LABEL org.opencontainers.image.description="Expedited Docker build for TEEMS solver"
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
