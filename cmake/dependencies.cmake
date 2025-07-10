find_package(PkgConfig REQUIRED)
include(FetchContent)

##################################################################################
# LAPACK

find_package(LAPACK REQUIRED)
find_package(LAPACKE REQUIRED)

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
