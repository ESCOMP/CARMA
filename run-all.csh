#! /bin/tcsh -f
# An entry point for running all of the carma tests. It creates a run directory
# copies in the executables and then runs all of the test executables.
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

# Create a directory for the build.
mkdir -p $rundir

# Copy the executable to the run directory.
cp $blddir/*TEST.exe $rundir


# Prepare for multiple threads, assuming Intel Compiler.
setenv OMP_NUM_THREADS $CARMA_THREADS
setenv KMP_STACKSIZE 128M

# Execute the tests.
cd $rundir

echo `ls -1 *TEST.exe`
foreach runtgt (`ls -1 *TEST.exe`)
  echo "  ** Started $runtgt at `date` **"
  ./$runtgt || echo '  *** Run Failed ***' && exit -1
  echo "  ** Finished at `date` **"
  echo ""
  
  set idlfile="read_`echo $runtgt:r | tr '[A-Z]' '[a-z]'`.pro"
  set outfile="carma_`echo $runtgt:r | tr '[A-Z]' '[a-z]'`.txt"
  
  if (! -f $idlfile) then
    cp -p ../../$testdir/$idlfile .
  endif

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
end


