FROM ubuntu:24.04

LABEL org.opencontainers.image.name="TEEMS solver full Dockerfile"
LABEL org.opencontainers.image.version="0.93"
LABEL org.opencontainers.image.description="Full Docker build for TEEMS solver"
LABEL org.opencontainers.image.author="Matthew Cantele <matthew.cantele@protonmail.com>"
LABEL org.opencontainers.image.maintainer="Matthew Cantele <matthew.cantele@protonmail.com>"
LABEL org.opencontainers.image.source="https://github.com/matthewcantele/teems-solver"
LABEL org.opencontainers.image.licenses="AGPL-3.0"

ARG TARGETARCH
ARG PATH_HSL_MA48
ARG PATH_HSL_MA51
ARG PATH_HSL_MC66
ARG PATH_HSL_MP48
ARG BUILD_DIR=/opt/teems-solver
# ARG OPENMPI_VERSION=5.0.7
# ARG OPENMPI_BASE_VERSION=5.0
ARG MPICH_VERSION=4.3.0
ARG PETSC_ARCH=native
ARG PETSC_VERSION=v3.24.4
ARG USERNAME=teems
ARG USER_UID=1000
ARG USER_GID=$USER_UID
ARG TEEMS_TMP=/tmp/teems

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="${BUILD_DIR}/lib/mpi/bin:$PATH"

COPY ${PATH_HSL_MA48} ${TEEMS_TMP}/ma48.tar.gz
COPY ${PATH_HSL_MA51} ${TEEMS_TMP}/ma51.tar.gz
COPY ${PATH_HSL_MC66} ${TEEMS_TMP}/mc66.tar.gz
COPY ${PATH_HSL_MP48} ${TEEMS_TMP}/mp48.tar.gz
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
        libopenblas-dev \
        git \
        python3 \
        python3-pip \
        python3-setuptools \
        tzdata \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    # detect architecture and set optimization flags
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
         amd64)   OPT_FLAGS="-Ofast -march=native -mtune=native" ;; \
         arm64)   OPT_FLAGS="-Ofast -march=armv8-a"              ;; \
         arm)     OPT_FLAGS="-Ofast -march=armv7-a -mfpu=neon"   ;; \
         ppc64le) OPT_FLAGS="-Ofast -mcpu=powerpc64le"           ;; \
         s390x)   OPT_FLAGS="-Ofast"                             ;; \
         *)       OPT_FLAGS="-O2"                                ;; \
       esac \
    && BLAS_LIB=$(find /usr/lib/$(gcc -print-multiarch) -name "libopenblas.so" | head -1) \
    # MPICH
    && cd ${TEEMS_TMP} \
    && mkdir mpich \
    && cd mpich \
    && wget https://www.mpich.org/static/downloads/${MPICH_VERSION}/mpich-${MPICH_VERSION}.tar.gz \
    && tar xfz mpich-${MPICH_VERSION}.tar.gz --strip-components=1 \
    && ./configure --prefix=${BUILD_DIR}/lib/mpi \
    && make \
    && make install \
    # PETSC
    && cd ${BUILD_DIR}/lib \
    && git clone -b ${PETSC_VERSION} https://gitlab.com/petsc/petsc.git petsc \
    && cd petsc \
    && mkdir -p ${BUILD_DIR}/LICENSES/petsc \
    && cp ${BUILD_DIR}/lib/petsc/LICENSE ${BUILD_DIR}/LICENSES/petsc \
    && python3 ./config/configure.py \
        PETSC_ARCH=${PETSC_ARCH} \
        PETSC_DIR=$PWD \
        --with-debugging=0 \
        --with-mpi-dir=${BUILD_DIR}/lib/mpi \
        COPTFLAGS="$OPT_FLAGS" \
        CXXOPTFLAGS="$OPT_FLAGS" \
        FOPTFLAGS="$OPT_FLAGS" \
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
    && cd ${TEEMS_TMP} \
    && mkdir -p ${BUILD_DIR}/LICENSES/hsl \
    && mkdir ${BUILD_DIR}/lib/ma48 \
    && tar xfz ma48.tar.gz --strip-components=1 -C ${BUILD_DIR}/lib/ma48 \
    && cp ${BUILD_DIR}/lib/ma48/LICENCE ${BUILD_DIR}/LICENSES/hsl \
    && mkdir ${BUILD_DIR}/lib/ma51 \
    && tar xfz ${TEEMS_TMP}/ma51.tar.gz --strip-components=1 -C ${BUILD_DIR}/lib/ma51 \
    && mkdir ${BUILD_DIR}/lib/mc66 \
    && tar xfz ${TEEMS_TMP}/mc66.tar.gz --strip-components=1 -C ${BUILD_DIR}/lib/mc66 \
    && mkdir ${BUILD_DIR}/lib/mp48 \
    && tar xfz ${TEEMS_TMP}/mp48.tar.gz --strip-components=1 -C ${BUILD_DIR}/lib/mp48 \
    && cd ${BUILD_DIR}/lib/ma48 \
    && autoreconf -fi \
    && ./configure \
    && make \
    && ln -s ${BUILD_DIR}/lib/ma48/src/.libs/libma48.so.0.0.0 /usr/local/lib/libma48.so \
    && ln -s ${BUILD_DIR}/lib/ma48/src/.libs/libma48.so.0.0.0 /usr/local/lib/libma48.so.0 \
    && cd ${BUILD_DIR}/lib/ma51 \
    && autoreconf -fi \
    && ./configure \
    && make \
    && ln -s ${BUILD_DIR}/lib/ma51/src/.libs/libma51.so.0.0.0 /usr/local/lib/libma51.so \
    && ln -s ${BUILD_DIR}/lib/ma51/src/.libs/libma51.so.0.0.0 /usr/local/lib/libma51.so.0 \
    && cd ${BUILD_DIR}/lib/mp48 \
    && autoreconf -fi \
    && ./configure FFLAGS="-fallow-argument-mismatch" FCFLAGS="-fallow-argument-mismatch" \
    && make \
    && cp ${BUILD_DIR}/lib/mp48/src/hsl_mp48d.f90 ${BUILD_DIR}/src \
    && cp ${BUILD_DIR}/lib/mp48/src/ddeps.f ${BUILD_DIR}/src \
    && cp ${BUILD_DIR}/lib/mp48/src/hsl_mp01.mod ${BUILD_DIR}/src \
    && cd ${BUILD_DIR}/src \
    && chmod +x mp48_mod.sh \
    && ./mp48_mod.sh \
    && cp ${BUILD_DIR}/lib/mc66/src/hsl_mc66d.f90 ${BUILD_DIR}/src \
    && cp ${BUILD_DIR}/lib/mc66/src/ddeps90.f90 ${BUILD_DIR}/src \
    && ldconfig \
    # final build
    && cd ${BUILD_DIR}/src \
    && make \
    && mv ${BUILD_DIR}/src ${BUILD_DIR}/solver \
    && rm -rf ${TEEMS_TMP} \
    # rename the default ubuntu user (UID/GID 1000) that ships with Ubuntu 24.04
    && usermod -l $USERNAME ubuntu \
    && groupmod -n $USERNAME ubuntu \
    && usermod -d /home/$USERNAME -m $USERNAME \
    # cleanup
    && apt-get autoremove -y \
    && apt-get autoclean

USER $USERNAME
WORKDIR /opt/teems
CMD ["/bin/bash"]
