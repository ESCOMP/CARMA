#! /bin/tcsh -f
# An entry point for viewing the results of the benchmark runs that
# are part of the CARMA distribution.
#
# Usage
#  view-bench.csh [view target]
#
#   build target - target label for the make
#

# Environment Variables
#   CARMA_BUILD [carma]
#     The subdirectory in which the build was performed.

if (! $?CARMA_IDL ) then
  setenv CARMA_IDL /Applications/itt/idl70/bin/idl
endif

set runtgt=CARMATEST.exe

if ($# == 1) then
  set runtgt="$1"
endif

set rundir=run/bench
set testdir=tests
set idlfile="read_`echo $runtgt:r | tr '[A-Z]' '[a-z]'`.pro"
set outfile="carma_`echo $runtgt:r | tr '[A-Z]' '[a-z]'`.txt"

echo $idlfile
echo $outfile

# Create a directory for the build.
echo "Viewing $testdir/bench/$outfile in the directory $rundir ..."
mkdir -p $rundir

# Copy the benchmark result to the run directory.
cp $testdir/bench/$outfile $rundir

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

# Execute the make file in the build directory.
cd $rundir

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


