#!/bin/sh

# Script to run PHAML tests of printed output

THISDIR="test_print"
echo "PHAML TEST:"
echo "PHAML TEST: $THISDIR: tests of printed output"

echo "PHAML TEST: Compile tests in $THISDIR"
make -s all
PASS=$?

if [ $PASS != 0 ]
then
   echo "PHAML TEST:"
   echo "PHAML TEST: FAILURE Unable to compile tests in $THISDIR."
   exit 1
fi

ERR=0
RUN1="$RUNMPI -np 1 "
if [ $PHAML_PARALLEL = "sequential" ]
then
   RUN1=
fi

case "$PHAML_PARALLEL" in
   messpass_spawn) FORM="ms" ;;
   messpass_nospawn) FORM="spmd" ;;
   sequential) FORM="seq" ;;
esac

for PROG in `ls test*.f90`;
do
   TESTN=`echo $PROG | sed -e s/\.f90//`
   echo "PHAML TEST: Run test $THISDIR/$TESTN"
   $RUN1 ./$TESTN.exe > $TESTN.out
   PASS=$?
   if [ $PASS != 0 ]
   then
      echo "PHAML TEST:"
      echo "PHAML TEST: FAILURE $THISDIR/$TESTN failed to run"
      ERR=1
   else
      if [ $TESTN = "test06" -o $TESTN = "test07" ]
      then
         if [ $PHAML_PARALLEL != "sequential" ]
         then
            diff $TESTN.0.out $TESTN.0.$FORM.comp > $TESTN.0.diff
            PASS=$?
         fi
         diff $TESTN.1.out $TESTN.1.$FORM.comp >> $TESTN.1.diff
         if [ $PASS != 2 ]
         then
            PASS=$?
         fi
         if [ $PHAML_PARALLEL != "sequential" ]
         then
            diff $TESTN.2.out $TESTN.2.$FORM.comp >> $TESTN.2.diff
            if [ $PASS != 2 ]
            then
               PASS=$?
            fi
         fi
      else
         diff $TESTN.out $TESTN.$FORM.comp > $TESTN.diff
         PASS=$?
      fi
      if [ $PASS = 2 ]
      then
         echo "PHAML TEST:"
         echo "PHAML TEST: FAILURE $THISDIR/$TESTN diff failed"
         ERR=1
      else
         if [ $TESTN = "test06" -o $TESTN = "test07" ]
         then
            if [ $PHAML_PARALLEL != "sequential" ]
            then
               DIFFSIZE=`cat $TESTN.0.diff | wc -l`
            fi
            if [ $DIFFSIZE = 0 ]
            then
               DIFFSIZE=`cat $TESTN.1.diff | wc -l`
            fi
            if [ $DIFFSIZE = 0 ]
            then
               if [ $PHAML_PARALLEL != "sequential" ]
               then
                  DIFFSIZE=`cat $TESTN.2.diff | wc -l`
               fi
            fi
         else
            DIFFSIZE=`cat $TESTN.diff | wc -l`
         fi
         if [ $DIFFSIZE = 0 ]
         then
            echo "PHAML TEST: SUCCESS for $THISDIR/$TESTN"
         else
            echo "PHAML TEST:"
            echo "PHAML TEST: WARNING -- $THISDIR/$TESTN.out differs from the expected"
            echo "PHAML TEST:            results in $THISDIR/$TESTN.$FORM.comp."
            echo "PHAML TEST:            Examine $THISDIR/$TESTN.diff"
            echo "PHAML TEST:            to see if the difference is significant."
            echo "PHAML TEST:"
         fi
      fi
   fi
done

exit $ERR
