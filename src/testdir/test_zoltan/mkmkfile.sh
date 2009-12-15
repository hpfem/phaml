#!/bin/bash

rm -f Makefile
f=Makefile

case "$PHAML_PARMETIS" in
   no)
      if [ -e test07.f90 ]
      then
         mv test07.f90 mtest07.f90
      fi ;;
   yes)
      if [ -e mtest07.f90 ]
      then
         mv mtest07.f90 test07.f90
      fi ;;
esac

case "$PHAML_DRUM" in
   no)
      if [ -e test08.f90 ]
      then
         mv test08.f90 mtest08.f90
      fi ;;
   yes)
      if [ -e mtest08.f90 ]
      then
         mv mtest08.f90 test08.f90
      fi ;;
esac

source ../mkinc/mkhead.inc

for PROG in `ls test*.f90` ;
do
   TESTN=`echo $PROG | sed -e s/\.f90//`
   echo "	$TESTN.exe \\" 1>>$f ;
done
echo "	phaml_slave" 1>>$f

source ../mkinc/mkslave.inc

for PROG in `ls test*.f90` ;
do
   TESTN=`echo $PROG | sed -e s/\.f90//`
   export TESTN
   source ../mkinc/mkmain.inc ;
done

source ../mkinc/mktail.inc
echo "	rm -f POWERFILE" 1>>$f
