#!/bin/sh

cd first
./mkmkfile.sh
cd ../pause
./mkmkfile.sh
cd ../debug
./mkmkfile.sh
cd ../graphics
./mkmkfile.sh
cd ..
for DIR in `ls -1 | grep test_` ;
do
   cd $DIR
   ./mkmkfile.sh
   cd ..
done

if [ -f Makefile ]
then
   mv -f Makefile Makefile.bak
fi

echo "nothing:" 1>>Makefile
echo "test:" 1>>Makefile
echo "	@ ./run_tests.sh $PHAML_ARCH $PHAML_OS $PHAML_F90 $PHAML_C $PHAML_HASHSIZE $PHAML_PARALLEL $PHAML_PARLIB $PHAML_GRAPHICS $PHAML_BLAS $PHAML_LAPACK $PHAML_ARPACK $PHAML_BLOPEX $PHAML_HYPRE $PHAML_MUMPS $PHAML_PETSC $PHAML_SUPERLU $PHAML_SYSTEM $PHAML_ZOLTAN $PHAML_PARMETIS $PHAML_JOSTLE $PHAML_PATOH $PHAML_PARKWAY $PHAML_NEMESIS $PHAML_DRUM $RUNMPI"' $(what)' 1>>Makefile
echo "clean:" 1>>Makefile
echo "	cd first; make clean; cd .." 1>>Makefile
echo "	cd pause; make clean; cd .." 1>>Makefile
echo "	cd debug; make clean; cd .." 1>>Makefile
echo "	cd graphics; make clean; cd .." 1>>Makefile
for DIR in `ls -1 | grep test_` ;
do
   echo "	cd $DIR; make clean; cd .." 1>>Makefile
done
echo "	rm -f testresults" 1>>Makefile
echo "reallyclean:" 1>>Makefile
echo "	cd first; make clean; cd .." 1>>Makefile
echo "	cd pause; make clean; cd .." 1>>Makefile
echo "	cd debug; make clean; cd .." 1>>Makefile
echo "	cd graphics; make clean; cd .." 1>>Makefile
for DIR in `ls -1 | grep test_` ;
do
   echo "	cd $DIR; make clean; cd .." 1>>Makefile
done
echo "	rm -f testresults" 1>>Makefile
echo "	cd first; rm -f Makefile Makefile.bak" 1>>Makefile
echo "	cd pause; rm -f Makefile Makefile.bak" 1>>Makefile
echo "	cd debug; rm -f Makefile Makefile.bak" 1>>Makefile
echo "	cd graphics; rm -f Makefile Makefile.bak" 1>>Makefile
for DIR in `ls -1 | grep test_` ;
do
   echo "	cd $DIR; rm -f Makefile Makefile.bak" 1>>Makefile
done
