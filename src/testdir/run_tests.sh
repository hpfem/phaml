#!/bin/sh

# top level script to run PHAML test programs

# invoked by "make test [what=something]" in the top directory or testdir.
# The optional something can be clean, first, debug, graphics, pause,
# interactive, noninteractive, everything (default), or any of the test_*
# directories.  "test_" can be omitted from the directory name.

echo
echo "PHAML TEST: -------------------------- PHAML tests --------------------"

####################### get PHAML environment

PHAML_ARCH=$1
export PHAML_ARCH
shift
PHAML_OS=$1
export PHAML_OS
shift
PHAML_F90=$1
export PHAML_F90
shift
PHAML_C=$1
export PHAML_C
shift
PHAML_HASHSIZE=$1
export PHAML_HASHSIZE
shift
PHAML_PARALLEL=$1
export PHAML_PARALLEL
shift
PHAML_PARLIB=$1
export PHAML_PARLIB
shift
PHAML_GRAPHICS=$1
export PHAML_GRAPHICS
shift
PHAML_BLAS=$1
export PHAML_BLAS
shift
PHAML_LAPACK=$1
export PHAML_LAPACK
shift
PHAML_ARPACK=$1
export PHAML_ARPACK
shift
PHAML_BLOPEX=$1
export PHAML_BLOPEX
shift
PHAML_HYPRE=$1
export PHAML_HYPRE
shift
PHAML_MUMPS=$1
export PHAML_MUMPS
shift
PHAML_PETSC=$1
export PHAML_PETSC
shift
PHAML_SUPERLU=$1
export PHAML_SUPERLU
shift
PHAML_SYSTEM=$1
export PHAML_SYSTEM
shift
PHAML_ZOLTAN=$1
export PHAML_ZOLTAN
shift
PHAML_PARMETIS=$1
export PHAML_PARMETIS
shift
PHAML_JOSTLE=$1
export PHAML_JOSTLE
shift
PHAML_PATOH=$1
export PHAML_PATOH
shift
PHAML_PARKWAY=$1
export PHAML_PARKWAY
shift
PHAML_NEMESIS=$1
export PHAML_NEMESIS
shift
PHAML_DRUM=$1
export PHAML_DRUM
shift
RUNMPI=$1
export RUNMPI

shift
what=$1
if [ "x$what" = "x" ]
then
   what=everything
fi
export what

if [ $what = "clean" ]
then
   make clean
   exit 0
fi

# TEMP some situations not yet supported

if [ $PHAML_PARALLEL = "messpass_nospawn" ]
then
   echo "PHAML TEST: Tests have not yet been written for PARALLEL=messpass_nospawn."
   echo "PHAML TEST: Quitting."
   exit 1
fi

####################### setup

if [ "x$PBS_NODEFILE" != "x" ]
then

# if started from a PBS script, don't allow interactive tests

   if [ $what = "everything" -o $what = "interactive" ]
   then
      echo "PHAML TEST: Cannot run the interactive tests from a batch system."
      echo "PHAML TEST: Quitting."
      exit 1
   fi

# if not started from a PBS script, request that the user start a demon if
# one is needed

else

   if [ $what != "noninteractive" ]
   then
      echo "PHAML TEST:"
      case "$PHAML_PARLIB" in
         lam)
            echo "PHAML TEST: LAM requires that a LAM demon (lamd) be running on"
            echo "PHAML TEST: each processor.  It is usually started with lamboot."
            echo "PHAML TEST: Is the demon ready (y/n)? "
            read ANS ;;
         mpich2)
            echo "PHAML TEST: MPICH2 might require that the mpd demon be running on"
            echo "PHAML TEST: each processor.  It is usually started with mpdboot."
            echo "PHAML TEST: Is the demon ready if needed(y/n)? "
            read ANS ;;
         openmpi)
            ANS="y" ;;
         none)
            ANS="y" ;;
         *)
            echo "PHAML TEST: If you are using a parallel environment that requires"
            echo "PHAML TEST: starting a demon, please start the demon now."
            echo "PHAML TEST: Is the demon ready (y/n)? "
            read ANS ;;
      esac
      if [ $ANS != "y" ]
      then
         echo "PHAML TEST: terminating because the answer was not y"
         exit 1
      fi
   fi
fi

####################### first tests

# the first tests just make sure a program can be run

if [ $what = "everything" -o $what = "noninteractive" -o $what = "interactive" -o $what = "first" ]
then
   cd first
   ./run_tests.sh
   PASS=$?

   make alittleclean
   if [ $PASS != 0 ]
   then
      echo "PHAML_TEST:"
      echo "PHAML TEST: Tests to see if a program can be run failed."
      echo "PHAML TEST: Quitting."
      cd ..
      exit 1
   fi

   cd ..
fi

####################### pause

if [ $what = "everything" -o $what = "interactive" -o $what = "pause" ]
then
   cd pause
   ./run_tests.sh
   PASS=$?
   make alittleclean
   cd ..
fi

####################### debug

if [ $what = "everything" -o $what = "interactive" -o $what = "debug" ]
then
   cd debug
   ./run_tests.sh
   PASS=$?
   make alittleclean
   cd ..
fi

####################### graphics

if [ $what = "everything" -o $what = "interactive" -o $what = "graphics" ]
then
   cd graphics
   ./run_tests.sh
   PASS=$?
   make alittleclean
   cd ..
fi

####################### noninteractive tests

rm -f testresults
touch testresults
for DIR in `ls -1 | grep test_` ;
do
   DIRROOT=`echo $DIR | sed -e s/test_//`
   if [ $what = "everything" -o $what = "noninteractive" -o $what = $DIR -o $what = $DIRROOT ]
   then
      cd $DIR
      ./run_tests.sh | tee -a ../testresults
      PASS=$?
      make alittleclean
      cd ..
   fi
done

echo "PHAML TEST:"
echo "PHAML TEST: All tests have been completed.  You can see the results of"
echo "PHAML TEST: the tests in the file testdir/testresults."

exit 0
