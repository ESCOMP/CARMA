#! /bin/tcsh -f

# Build the legacy CARMA targets with debugger-friendly compiler flags.
# The output is isolated from the normal build/carma directory.
#
# Usage:
#   make-carma-debug.csh [build target]
#
# Examples:
#   make-carma-debug.csh MIXEDSULFATE_2NC_TEST.exe
#   make-carma-debug.csh all

if ($# > 1) then
  echo "Usage: make-carma-debug.csh [build target]"
  exit 2
endif

setenv CARMA_BUILD carma-debug
setenv CARMA_DEBUG 1

if ($# == 0) then
  ./make-carma.csh
else
  ./make-carma.csh "$1"
endif

set build_status=$status
exit $build_status