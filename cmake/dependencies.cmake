find_package(PkgConfig REQUIRED)
include(FetchContent)

##################################################################################
# LAPACK

# Try to find BLAS first, then LAPACK
# Use shared libraries preference to avoid static library issues
set(BLA_STATIC OFF)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# ##############################################################################
# Memory check

if(CARMA_ENABLE_MEMCHECK)
  find_file(
    MEMCHECK_SUPPRESS_FILE
    DOC "Suppression file for memory checking"
    NAMES openmpi-valgrind.supp
    PATHS /usr/share/openmpi /usr/lib64/openmpi/share
          /usr/lib64/openmpi/share/openmpi /usr/share)
  if(MEMCHECK_SUPPRESS_FILE)
    set(MEMCHECK_SUPPRESS
        "--suppressions=${PROJECT_SOURCE_DIR}/test/valgrind.supp --suppressions=${MEMCHECK_SUPPRESS_FILE}"
    )
  else()
    set(MEMCHECK_SUPPRESS
        "--suppressions=${PROJECT_SOURCE_DIR}/test/valgrind.supp")
  endif()
endif()

##################################################################################
# NetCDF

find_package(PkgConfig REQUIRED)

pkg_check_modules(netcdff REQUIRED IMPORTED_TARGET netcdf-fortran)
pkg_check_modules(netcdfc REQUIRED IMPORTED_TARGET netcdf)

##################################################################################
# Google Test

if(PROJECT_IS_TOP_LEVEL)
  FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG be03d00f5f0cc3a997d1a368bee8a1fe93651f48)

  set(INSTALL_GTEST
      OFF
      CACHE BOOL "" FORCE)
  set(BUILD_GMOCK
      OFF
      CACHE BOOL "" FORCE)

  FetchContent_MakeAvailable(googletest)
endif()
