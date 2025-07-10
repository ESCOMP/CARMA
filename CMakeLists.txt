cmake_minimum_required(VERSION 3.21)

project(
  carma
  VERSION 4.0.0
  LANGUAGES Fortran
)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${PROJECT_SOURCE_DIR}/cmake")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

message ( STATUS "CMake build configuration for tuv-x(${CMAKE_BUILD_TYPE}) ${PROJECT_VERSION}" )

################################################################################
# Project wide setup options

option(CARMA_ENABLE_MEMCHECK "Enable memory checking in tests" OFF)

set(CARMA_MOD_DIR "${PROJECT_BINARY_DIR}/include" CACHE PATH "Directory to find Fortran module files during build")
set(CARMA_INSTALL_MOD_DIR "${INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}" CACHE PATH "Directory to install Fortran module files")
set(CARMA_INSTALL_INCLUDE_DIR "${INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}" CACHE PATH "Directory to install CARMA include files")

# Set up include and lib directories
set(CARMA_LIB_DIR "${PROJECT_BINARY_DIR}/lib")
include(GNUInstallDirs)
set(INSTALL_PREFIX "carma-${PROJECT_VERSION}")

################################################################################
# Dependencies

include(cmake/dependencies.cmake)

################################################################################
# CAMRA targets

add_subdirectory(src)

################################################################################
# Packaging

# Only include pacakaging if we are the top-level project
if(PROJECT_IS_TOP_LEVEL)
  add_subdirectory(packaging)
endif()
