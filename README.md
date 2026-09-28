# Community Aerosol and Radiation Model for Atmospheres (CARMA)
[![License](https://img.shields.io/github/license/ESCOMP/CARMA.svg)](https://github.com/ESCOMP/CARMA/blob/main/LICENSE)
[![ubuntu](https://github.com/ESCOMP/CARMA/actions/workflows/Ubuntu.yml/badge.svg)](https://github.com/ESCOMP/CARMA/actions/workflows/Ubuntu.yml)

### Getting the code

```
git clone --recursive https://github.com/ESCOMP/CARMA.git
```

This populates source/base with source code from https://github.com/ESCOMP/CARMA_base.

___
## README

This project contains files to build version 4.0 of the Community Aerosol and
Radiation Model for Atmospheres (CARMA) that is based off of the 2.3 release,
but has been ported to Fortran 90 and repackaged so that it can be used as a
cloud and aerosol physics package embedded into GCMs. The project consists of
4 components that are each located in their own subdirectories:

  - bin    : the F90doc documentation generation tool
  - doc    : the change log and HTML based documentation
  - source : the CARMA microphysics layer
  - tests  : the test routines & benchmark results

Two script files have been provided to build and run the model or its
components. To build carma and the test cases, issue the following command
rom the root directory:
```
  make-carma.csh
```
This will build all the files in a subdirectory called build/carma.
Note: LAPACK external library is required for the CARMA base code
(This avoids the used of the copied LU matrix solve routines from copyrighted source)

To run a sample carma model, execute the following command from the root directory:
```
  run-carma.csh
```
This will copy the CARMA.exe executable to the directory run/carma and then
will execute it with all output going to the run/carma directory. The
dependency hierarchy is set up in the make files, so rebuilding should only
rebuild what is necessary. There are only two #defines that can be set in the
makefile to control the build. One specifies the precision (SINGLE) and the other
is used to build a version with extra debug information (DEBUG). The scripts can
also be used to build individual executables. For example, to build and run the
NUCTEST test routine and to have the run performed in a directory named nuctest
(this part is optional) execute:
```
  make-carma.csh NUCTEST.exe
  run-carma.csh NUCTEST.exe
```

For a legacy Makefile build with a separate debugger-friendly output directory,
use `make-carma-debug.csh`:

```
make-carma-debug.csh MIXEDSULFATE_2NC_TEST.exe
gdb ./build/carma-debug/MIXEDSULFATE_2NC_TEST.exe
```

This sets `CARMA_BUILD=carma-debug`, enables the `DEBUG` preprocessor symbol,
and retains the existing compiler debug flags for the selected Fortran compiler.
NOTE: bash and ksh users will need to use export rather than setenv.

Plotting routines have been created for each test in IDL. If you set CARMA_IDL
to be the path to your copy of IDL, then IDL should be launched automatical, and
the messages will indicate the program that you should run. On the Mac, you would
set a path like:
```
  setenv CARMA_IDL /Applications/itt/idl/idl80/bin/idl
```
After exectuing
```
 run-carma.csh NUCTEST.exe
```
you should see an IDL prompt that looks like this:

```
    ** Finished at Wed Jul 13 20:41:11 MDT 2011 **

  Running the IDL analysis routine read_nuctest.pro
    To run the test, in IDL you need to type the command: .r read_nuctest.pro
    To exit IDL, type the command: exit

  IDL Version 8.0, Mac OS X (darwin x86_64 m64). (c) 2010, ITT Visual Information Solutions
  Installation number: 95183-2518.
  Licensed for use by: ACD:acd-license:Linux.FL

  IDL>
```
Type:
```
  .r read_nuctest.pro
```
and you should get a plot generated showing the results of the test.


Documentation of the code is generated automatically by the make scripts using
the program f90doc, which is located in the bin directory. The generated
documentation is stored as html files in the directory doc/f90doc. To start
browsing the documentation, open the file doc/index.html.


The model supports OPEN/MP, which will allow the model to use multiple
processors if the model is called with multiple columns. You need to add
a compiler directive to enable OPEN/MP (-openmp for ifort) and you need
to specify the number of threads when you run. The  run script uses the
CARMA_THREADS flag to determine how many threads to allow during execution.
The dault is 1 thread. To run with 4 threads issue the command
```
  setenv CARMA_THREADS 4
```
before using the run script.

Several test cases have been developed to show how the CARMA model
can be used and to test that the model physics is performing correctly. The
source code for these tests are located in the tests subdirectory.

### Developing and debugging CMake tests

The CMake test targets are named `test_<name>` and are registered with CTest.
To add and debug a new test:

1. Add the test source file under `tests/`. The source must contain a Fortran
  program, for example `carma_mytest.F90`.
2. Choose a short CMake name. The name becomes both the CTest name and the
  executable suffix, so `mytest` produces `test_mytest`.
3. Register the test in `tests/CMakeLists.txt` with `create_standard_test`.
  Include `atmosphere_mod.F90` when the test follows the usual CARMA test
  pattern:

```
create_standard_test(NAME mytest SOURCES atmosphere_mod.F90 carma_mytest.F90)
```

4. Add every other source file that provides a module used by the test. For
  example, a NetCDF test commonly needs `test2nc_mod.F90`, `nc_types_mod.F90`,
  `ncio_mod.F90`, and `carmadiags_mod.F90`; a test using a shared local module
  must list that module source too. Put NetCDF-only registrations inside the
  existing `if(CARMA_ENABLE_NETCDF)` block.
5. Configure the debug build with NetCDF enabled. CMake 3.21 or newer is
  required. On the Cray environment, set the MPICH pkg-config path first so
  CMake can find the parallel NetCDF dependencies. If `cmake --version`
  reports an older version, use the supported CMake executable configured for
  this environment:

```sh
export PKG_CONFIG_PATH="/opt/cray/pe/mpich/8.1.27/ofi/gnu/9.1/lib/pkgconfig:${PKG_CONFIG_PATH}"
/mnt/lustre/a2fs-work4/work/y07/shared/utils/core/cmake/3.29.4/bin/cmake -S . -B build-debug -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=ftn "-DCMAKE_Fortran_FLAGS_DEBUG=-ffree-line-length-none -g -O0 -fbacktrace -fbounds-check" -DCARMA_ENABLE_NETCDF=ON
```

In the multiline command, each `\` must be the final character on its line;
do not put spaces after it. If copying multiline commands is inconvenient,
run the same command on one line instead.

6. Build only the new executable. Replace `mytest` with the `NAME` from step 3:

```sh
cmake --build build-debug --target test_mytest --parallel
```

7. Start debugging in VS Code by selecting `Debug CARMA test` in the Run and
  Debug view and entering `mytest` when prompted. The configuration builds
  `build-debug/test_mytest`, starts GDB, and uses `build-debug` as the working
  directory. Set breakpoints in the Fortran source before launching.

8. If the test writes files, check that its output appears in `build-debug`.
  You can also run the executable directly:

```sh
./build-debug/test_mytest
```

To run all of the tests interactively, use the command
```
  run-all.csh
```
To run all of the tests in the background and compare the results to
benchmark results, use the command
```
  run-regress.csh
```
The benchmarks are stored in tests/bench and were created on a Mac
using the Intel Fortran compiler. You may get slightly different
results on other platforms or with other compilers.


The entire project can be put into a tar file, using the command
```
  make-carma.csh tar
```
This will create a tar file called CARMA.tar in the current build directory,
which defaults to build/carma. NOTE: This tar file does not contain the
contents of the build, doc or run directories since they are
products of the make and run scripts.

CARMA has a ChangeLog in the doc directory. It contains the revision history
of CARMA. When changes are made, please prepend the ChangeLog with a
description of the change based upon filling out the ChangeLog_template
which is also stored in the doc directory.


Chuck Bardeen, Pete Colarco and Jamie Smith
Jul-2011

## CMake Build

To build CARMA with CMake locally and run the tests:
```
git clone https://github.com/ESCOMP/CARMA.git
cd CARMA
mkdir build
cd build
cmake ..
make
make test
```

Note that the CMake tests only run the test configurations, but do not verify results.
To include memory checking with Valgrind for the tests, add the `CARMA_ENABLE_MEMCHECK`
flag to the cmake step:
```
cmake -D CARMA_ENABLE_MEMCHECK=ON ..
```

To exclude NetCDF as a dependency:
```
cmake -D CARMA_ENABLE_NETCDF=OFF ..
```
Note that this will disable some of the tests.

## Docker Build

To build and run CARMA with Docker:
```
git clone https://github.com/ESCOMP/CARMA.git
cd CARMA
docker build -t carma . -f docker/Dockerfile
docker run carma bash
make test
```

Build args can be used to set CARMA build options
| Build Arg | Values |
|-----------|--------|
| BUILD_TYPE | Release / Debug |
| ENABLE_NETCDF | ON / OFF |
| ENABLE_MEMCHECK | ON / OFF |

To include the args modify the `docker build` step:
```
docker build -t carma . -f docker/Dockerfile --build-arg ENABLE_NETCDF=OFF --build-arg BUILD_TYPE=Debug
```
