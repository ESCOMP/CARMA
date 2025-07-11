FROM fedora:latest

RUN dnf -y update \
    && dnf -y install \
        cmake \
        gcc \
        gcc-c++ \
        gcc-fortran \
        gdb \
        git \
        lapack-devel \
        make \
        netcdf-fortran-devel \
        python3 \
        python3-pip \
        valgrind \
    && dnf clean all

ENV LD_LIBRARY_PATH=/usr/local/lib64

# build CARMA
COPY . /carma
RUN mkdir build \
      && cd build \
      && cmake -D CMAKE_BUILD_TYPE=release \
               -D CARMA_INSTALL_INCLUDE_DIR=/usr/local/include \
               -D CARMA_INSTALL_MOD_DIR=/usr/local/include \
               /carma \
      && make install -j 8

WORKDIR /build