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

# Initialize string to capture user y/n response
set doall=""

# Initialize index to count number of pairs of differing files to 0
set ndiffs=0

# Clear the screen and advance two lines
clear
echo ""
echo ""

# Inquire if the user wants to to have the script stop on the first instance of
#	a difference found between a test file and its corresponding benchmark file
#
echo "     Do you want the script to stop on the first difference found between "
read -p "        a test file and its corresponding benchmark file? (y/n) " $doall

# Create a directory for the build.
mkdir -p $rundir

# Copy the executable to the run directory.
cp $blddir/*TEST.exe $rundir

# Prepare for multiple threads, assuming Intel Compiler.
setenv OMP_NUM_THREADS $CARMA_THREADS
setenv OMP_STACKSIZE 128M

# Execute the tests.
cd $rundir

foreach runtgt (`ls -1 *TEST.exe`)
  echo ""
  echo ""
  echo "  ** Starting $runtgt at `date` **"
  #Run the test program; if the test fails exit the script
  ./$runtgt || echo '  *** Run Failed ***' && exit -1
  echo "  ** Finished at `date` **"
  echo ""

  set outfile="carma_`echo $runtgt:r | tr '[A-Z]' '[a-z]'`.txt"

  setenv FDIFF -sqw

  # The diff on AIX doesn't have the -q option ..."
  if (`uname` == AIX ) then
    setenv FDIFF -sw
  endif

  # check for exisitence of output file and set up the diff-grep string
  if (-f $outfile) then
       set grepdiff=`diff $FDIFF $outfile ../../$benchdir/$outfile | grep -o differ`

       # if a difference is found between the test file and its benchmark file...
       if $grepdiff != "" then
           # ... then print out a message to that effect
           echo "     " $outfile "differs from benchmark file ../../"$benchdir/$outfile "     "
           echo ""
           #if the user wants the script to stop on the first difference, then do so...
           if ( $doall == "y" || $doall == "Y" ) then
			   echo "      Quitting."
			   echo ""
			   exit -1
           else
		       # otherwise, increment the number of differing file pairs found
		       set ndiffs = `expr $ndiffs + 1`
           endif
       endif
    endif
  endif
end

echo ""
echo "    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo ""
# How many differences were found between test and benchmark files?
echo "      There are " $ndiffs "test files that differ from their corresponding benchmark files."
echo ""
echo ""

