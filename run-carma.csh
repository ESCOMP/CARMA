#! /bin/tcsh -f
# An entry point for running the CARMA code, which creates a run directory
# copies in the executable and the then runs the executable.
#
# NOTE: This script could be easily enhanced to manage mutliple executables and
# multiple run directories, which my be useful as part of an automated test
# suite.
#
# Usage
#  run-carma.csh [test]
#
#   test - the name of an executable test case
#

# Environment Variables
#   CARMA_BUILD [carma]
#     The subdirectory in which the build was performed.

# By default, run the single column model
if (! $?CARMA_BUILD ) then
  setenv CARMA_BUILD carma
endif

if (! $?CARMA_CASE ) then
  setenv CARMA_CASE $CARMA_BUILD
endif

if (! $?CARMA_THREADS ) then
  setenv CARMA_THREADS 1
endif

if (! $?CARMA_IDL ) then
  setenv CARMA_IDL idl
endif

set runtgt=CARMATEST.exe

if ($# == 1) then
  set runtgt="$1"
endif

set blddir=build/$CARMA_BUILD
set rundir=run/$CARMA_CASE
set testdir=tests
set idlfile="read_`echo $runtgt:r | tr '[A-Z]' '[a-z]'`.pro"
set outfile="carma_`echo $runtgt:r | tr '[A-Z]' '[a-z]'`.txt"

echo $idlfile
echo $outfile

# Create a directory for the build.
echo "Running $blddir/$runtgt in the directory $rundir using $CARMA_THREADS thread(s) ..."
mkdir -p $rundir

# Copy the executable to the run directory.
cp $blddir/$runtgt $rundir

if (-f $testdir/$idlfile) then
  
  # Don't overwrite the file in the run directory if it is newer than
  # the one in the test directory.
  #
  # NOTE: For the test on modification date to work, the copy must
  # preserve the modify date and time.
  if (-f $rundir/$idlfile) then
    if (-M "$rundir/$idlfile" > -M "$testdir/$idlfile") then
      setenv IDL_WARNING "  WARNING: $idlfile not copied, since $rundir/$idlfile is newer than $testdir/$idlfile"
    else
      cp -p $testdir/$idlfile $rundir
    endif
  else
    cp -p $testdir/$idlfile $rundir
  endif
endif

# Prepare for multiple threads, assuming Intel Compiler.
setenv OMP_NUM_THREADS $CARMA_THREADS
setenv KMP_STACKSIZE 128M

# Execute the make file in the build directory.
cd $rundir
echo "  ** Started at `date` **"
./$runtgt || echo '  *** Run Failed ***' && exit -1
echo "  ** Finished at `date` **"

if (-f $idlfile) then

  if (-f $outfile) then
    echo ""
    echo "Running the IDL analysis routine $idlfile"
  
    if ($?IDL_WARNING) then
      echo ""
      echo "$IDL_WARNING"
      echo ""
    endif
  
    echo "  To run the test, in IDL you need to type the command: .r $idlfile"
    echo "  To exit IDL, type the command: exit"
    
    # NOTE: If your invokation of IDL fails, check to see whether idl
    # is really on you path or if it is just an alias. Aliases don't work
    # properly in scripts, but this is how IDL is setup be default. You
    # can add the idl bin directory to your path so that this will work.
    echo ""
    $CARMA_IDL
  endif
endif


