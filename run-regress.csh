#! /bin/tcsh -f
# An entry point for running all of the carma regression tests. It creates
# a run directory, copies in the executables, runs all of the test executables,
# and then compares the results to the previous "benchmark" result.
#
# The benchmark results are run on a Mac using ifort
#
#
# Usage
#  run-carma.csh

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
set benchdir=tests/bench

# Create a directory for the build.
mkdir -p $rundir

# Copy the executable to the run directory.
cp $blddir/*TEST.exe $rundir

# Prepare for multiple threads, assuming Intel Compiler.
setenv OMP_NUM_THREADS $CARMA_THREADS
setenv KMP_STACKSIZE 128M

# Execute the tests.
cd $rundir

foreach runtgt (`ls -1 *TEST.exe`)
  echo ""
  echo ""
  echo "  ** Starting $runtgt at `date` **"
  ./$runtgt || echo '  *** Run Failed ***' && exit -1
  echo "  ** Finished at `date` **"
  echo ""
  
  set outfile="carma_`echo $runtgt:r | tr '[A-Z]' '[a-z]'`.txt"
  
  setenv FDIFF -sqw

  # The diff on AIX doesn't have the -q option ..."
  if (`uname` == AIX ) then
    setenv FDIFF -sw
  endif
  
  if (-f $outfile) diff $FDIFF $outfile ../../$benchdir/$outfile || exit(-1)
end

echo ""
echo ""
echo "All Tests Passed!"


