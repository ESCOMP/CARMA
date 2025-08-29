FROM fedora:latest

ARG BUILD_TYPE=Debug

RUN dnf -y update \
    && dnf -y install \
        cmake \
        gcc-c++ \
        gfortran \
        gdb \
        git \
        lapack-devel \
        lcov \
        make \
        netcdf-fortran-devel \
        python \
        valgrind \
        tree \
    && dnf clean all

# copy CARMA
COPY . /carma/

# build CARMA
RUN mkdir /build \
    && cd /build \
    && cmake \
       -D CMAKE_BUILD_TYPE=${BUILD_TYPE} \
       -D CARMA_ENABLE_MEMCHECK=ON \
       /carma \
    && make

WORKDIR /build