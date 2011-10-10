#! /bin/tcsh -f

# An entry point for build the CARMA code, which creates the build directory
# and populates it with the base makefile.
#
# NOTE: This script could be easily enhanced to manage multiple targets and
# multiple build directories, which may be useful as part of an automated test
# suite.
#
# Usage
#  make-carma.csh [build target]
#
#   build target - target label for the make
#
# Environment Variables
#   CARMA_BUILD [carma]
#     The subdirectory in which the build will be performed.

# Look for gmake first, but if not found then just use make.
setenv MAKETOOL "`which gmake`"
echo $MAKETOOL
if ("`echo $MAKETOOL | grep 'not found'`" != "") then
  setenv MAKETOOL "`which make`"
endif

echo "Using :  " $MAKETOOL
echo ""

# By default, build all of the targets in a directory named build.
set bldtgt=all

if ($# == 1) then
  set bldtgt="$1"
endif

if (! $?CARMA_BUILD ) then
  setenv CARMA_BUILD carma
endif

set blddir=build/$CARMA_BUILD
set docdir=doc/f90doc

# Create a directory for the build.
echo "Building the target $bldtgt in the directory $blddir"
mkdir -p $blddir

# Copy the makefile to the build directory.
cp Makefile $blddir/Makefile

# In MacOSX, the tar command will try to store off the resource fork as an
# extra file which is the same name as the original file with a ._ prefix.
# This flag stops that behavior
if ($bldtgt == tar) then
  setenv COPYFILE_DISABLE TRUE
endif

# Execute the make file in the build directory.
cd $blddir
$MAKETOOL $bldtgt

# Create the documentation.
if ($bldtgt != tar) then
  echo "Creating the documentation in the directory $docdir"

  # Create a directory for the build.
  cd ../..
  mkdir -p $docdir

  # Copy the makefile to the doc directory.
  cp Makefile $docdir/Makefile

  cd $docdir
  $MAKETOOL doc
endif

