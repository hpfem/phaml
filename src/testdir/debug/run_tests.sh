#!/bin/sh

# Script to run PHAML debug tests.

THISDIR="debug"
echo "PHAML TEST:"
echo "PHAML TEST: $THISDIR: tests of spawning with a debugger."

if [ $PHAML_PARALLEL != "messpass_spawn" ]
then
   echo "PHAML TEST: PHAML_PARALLEL is not messpass_spawn; spawning with debugger does not apply."
   echo "PHAML TEST: Skipping test."
   exit 0
fi

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
echo "PHAML TEST: The next programs test spawning the slaves under a debugger."
echo "PHAML TEST: As each program starts, two xterms will pop up, each running"
echo "PHAML TEST: a debugger.  In each window, perform whatever action the"
echo "PHAML TEST: debugger requires to start running the program; usually this"
echo "PHAML TEST: is entering 'run' at the debugger command prompt.  When the"
echo "PHAML TEST: program finishes you should be returned to the debugger"
echo "PHAML TEST: control in each window.  Terminate the debuggers, usually by"
echo "PHAML TEST: entering 'quit' or 'exit' at the command prompt, and the"
echo "PHAML TEST: xterm windows should go away."
echo "PHAML TEST:"
echo "PHAML TEST: Are you ready to run the tests (y)? "
read ANS

RUN1="$RUNMPI -np 1 "
# TEMP stdin doesn't seem to be getting through LAM mpirun, but I know I can
# run a LAM master/slave program without mpirun.
if [ $PHAML_PARLIB = "lam" ]
then
   RUN1=
fi

# run a test with gdb if it exists

echo "PHAML TEST: $THISDIR/test01 tests gdb"
which gdb > /dev/null >& /dev/null
PASS=$?
if [ $PASS != 0 ]
then
   echo "PHAML TEST: gdb does not exist; skipping $THISDIR/test01"
else
   echo "PHAML TEST: Do you want to run $THISDIR/test01 (y/n)? "
   read ANS
   if [ $ANS = "y" ]
   then
      echo "PHAML TEST: Run test $THISDIR/test01"
      $RUN1 ./test01.exe
      PASS=$?
      if [ $PASS != 0 ]
      then
         echo "PHAML TEST:"
         echo "PHAML TEST: FAILURE $THISDIR/test01 failed to run"
      else
         echo "PHAML TEST: SUCCESS for $THISDIR/test01"
      fi
   fi
fi

# run a test with the compiler supplied debugger if it exists

echo "PHAML TEST: $THISDIR/test02 tests a compiler supplied debugger"
case "$PHAML_F90" in
   absoft)
      DBG="Fx2" ;;
   intel)
      DBG="idb" ;;
   lahey)
      DBG="fdb" ;;
   nag)
      DBG="upsf95" ;;
   pathscale)
      DBG="pathdb" ;;
   pgi)
      DBG="pgdbg" ;;
   sgi|sun|xlf)
      DBG="dbx" ;;
   *)
      DBG="none" ;;
esac

if [ $DBG = "none" ]
then
   echo "PHAML TEST: Fortran compiler $PHAML_F90 does not supply a debugger."
   echo "PHAML TEST: Skipping $THISDIR/test02"
else
   which $DBG > /dev/null >& /dev/null
   PASS=$?
   if [ $PASS != 0 ]
   then
      echo "PHAML TEST: Fortran compiler is $PHAML_F90 but $DBG does not exist."
      echo "PHAML TEST: Skipping $THISDIR/test02"
   else
      echo "PHAML TEST: Fortran compiler is $PHAML_F90 and $DBG exists."
      echo "PHAML TEST: Do you want to run $THISDIR/test02 (y/n)? "
      read ANS
      if [ $ANS = "y" ]
      then
         echo "PHAML TEST: Run test $THISDIR/test02"
         echo $DBG > dbg2_command
         $RUN1 ./test02.exe
         PASS=$?
         if [ $PASS != 0 ]
         then
            echo "PHAML TEST:"
            echo "PHAML TEST: FAILURE $THISDIR/test02 failed to run"
         else
            echo "PHAML TEST: SUCCESS for $THISDIR/test02"
         fi
      fi
   fi
fi

# test the user's favorite debugger

echo "PHAML TEST:"
echo "PHAML TEST: You can also run a test with a debugger that you specify."
echo "PHAML TEST: Do you want to specify a debugger (y/n)? "
read ANS
echo "PHAML TEST:"

if [ $ANS = "y" ]
then
   echo "PHAML TEST: Enter the executable for the debugger: "
   read ANS
   which $ANS > /dev/null >& /dev/null
   PASS=$?
   if [ $PASS != 0 ]
   then
      echo "PHAML TEST: $ANS: Command not found.  Skipping $THISDIR/test03."
   else
      echo "PHAML TEST: Run test $THISDIR/test03"
      echo $ANS > dbg3_command
      $RUN1 ./test03.exe
      if [ $PASS != 0 ]
      then
         echo "PHAML TEST:"
         echo "PHAML TEST: FAILURE $THISDIR/test03 failed to run"
      else
         echo "PHAML TEST: SUCCESS for $THISDIR/test03"
      fi
   fi
fi

exit 0
