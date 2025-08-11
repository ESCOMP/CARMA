FROM fedora:latest

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
    && cmake -D \
       CARMA_ENABLE_MEMCHECK=ON \
       /carma \
    && make

WORKDIR /build