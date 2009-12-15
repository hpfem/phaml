#!/bin/bash

rm -f Makefile
f=Makefile

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
