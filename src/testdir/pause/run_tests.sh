#!/bin/sh

# Script to run PHAML pause tests.

THISDIR="pause"
echo "PHAML TEST:"
echo "PHAML TEST: $THISDIR: tests of the pause variables."

echo "PHAML TEST: Compile tests in $THISDIR"
make -s all
PASS=$?

if [ $PASS != 0 ]
then
   echo "PHAML TEST:"
   echo "PHAML TEST: FAILURE Unable to compile tests in $THISDIR."
   exit 1
fi

echo "PHAML TEST:"
echo "PHAML TEST: This program tests the pause variables in phaml_solve_pde."
echo "PHAML TEST: As the program runs you will be asked to 'press return to continue'"
echo "PHAML TEST: Whenever this prompt appears, press the Return or Enter key."
echo "PHAML TEST: This should occur 4 times."
echo "PHAML TEST:"
echo "PHAML TEST: Are you ready to run the test (y)? "
read ANS

ERR=0
RUN1="$RUNMPI -np 1 "
# TEMP stdin doesn't seem to be getting through LAM mpirun, but I know I can
# run a LAM master/slave program without mpirun.
if [ $PHAML_PARLIB = "lam" ]
then
   RUN1=
fi

if [ $PHAML_PARALLEL = "sequential" ]
then
   RUN1=
fi

for PROG in `ls test*.f90`;
do
   TESTN=`echo $PROG | sed -e s/\.f90//`
   echo "PHAML TEST: Run test $THISDIR/$TESTN"
   $RUN1 ./$TESTN.exe
   PASS=$?
   if [ $PASS != 0 ]
   then
      echo "PHAML TEST:"
      echo "PHAML TEST: FAILURE $THISDIR/$TESTN failed to run"
      ERR=1
   else
      echo "PHAML TEST: SUCCESS for $THISDIR/$TESTN"
   fi
done

exit $ERR
