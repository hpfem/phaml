#!/bin/sh

# Script to run PHAML debug tests.

THISDIR="graphics"
echo "PHAML TEST:"
echo "PHAML TEST: $THISDIR: tests of graphics."

if [ $PHAML_GRAPHICS = "none" ]
then
   echo "PHAML TEST: PHAML_GRAPHICS is none."
   echo "PHAML TEST: Skipping tests."
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
echo "PHAML TEST: The next programs test the graphics.  In each test, one or"
echo "PHAML TEST: more PHAML graphics windows will appear.  The window can be"
echo "PHAML TEST: moved, iconified, resized, etc., in the usual manner."
echo "PHAML TEST: The left mouse button is used to rotate the image.  The middle"
echo "PHAML TEST: mouse button is used to zoom the image.  The arrow keys are"
echo "PHAML TEST: used to pan the image.  The right mouse button brings up"
echo "PHAML TEST: a menu of many things that can be done with the image.  You"
echo "PHAML TEST: should experiment with the items in the menu to verify"
echo "PHAML TEST: things work.  When done, press the Return or Enter key in"
echo "PHAML TEST: the window with the message 'press return to continue', which"
echo "PHAML TEST: will usually be this window.  This will terminate the graphics"
echo "PHAML TEST: window.  Do not terminate the graphics window by any other means."
echo "PHAML TEST: It may be a few seconds before an image is drawn in the window."
echo "PHAML TEST:"
echo "PHAML TEST: Are you ready to run the tests (y)? "
read ANS

RUN1="$RUNMPI -np 1 "
if [ $PHAML_PARALLEL = "sequential" ]
then
   RUN1=
fi
# TEMP stdin doesn't seem to be getting through LAM mpirun, but I know I can
# run a LAM master/slave program without mpirun.
if [ $PHAML_PARLIB = "lam" ]
then
   RUN1=
fi

echo "PHAML TEST: $THISDIR/test01 tests graphics from the master."
echo "PHAML TEST: Run test $THISDIR/test01"
if [ $PHAML_PARALLEL = "sequential" ]
then
   rm -f /tmp/phaml_message
   phaml_graphics &
fi
$RUN1 ./test01.exe
PASS=$?
if [ $PASS != 0 ]
then
   echo "PHAML TEST:"
   echo "PHAML TEST: FAILURE $THISDIR/test01 did not run."
else
   echo "PHAML TEST: SUCCESS for $THISDIR/test01."
fi

sleep 2
echo "PHAML TEST: $THISDIR/test02 tests graphics from the slaves."
echo "PHAML TEST: Run test $THISDIR/test02"
if [ $PHAML_PARALLEL = "sequential" ]
then
   rm -f /tmp/phaml_message
   phaml_graphics &
fi
$RUN1 ./test02.exe
PASS=$?
if [ $PASS != 0 ]
then
   echo "PHAML TEST:"
   echo "PHAML TEST: FAILURE $THISDIR/test02 did not run."
else
   echo "PHAML TEST: SUCCESS for $THISDIR/test02."
fi

sleep 2
echo "PHAML TEST: $THISDIR/test03 tests graphics from both the master and the slaves."
echo "PHAML TEST: Run test $THISDIR/test03"
if [ $PHAML_PARALLEL = "sequential" ]
then
   rm -f /tmp/phaml_message
   phaml_graphics &
fi
$RUN1 ./test03.exe
PASS=$?
if [ $PASS != 0 ]
then
   echo "PHAML TEST:"
   echo "PHAML TEST: FAILURE $THISDIR/test03 did not run."
else
   echo "PHAML TEST: SUCCESS for $THISDIR/test03."
fi

sleep 2
echo "PHAML TEST: $THISDIR/test04 tests graphics after each phase and pausing."
echo "PHAML TEST:"
echo "PHAML TEST: In this test you should press the Enter or Return key each"
echo "PHAML TEST: time the 'press return to continue' prompt appears.  You"
echo "PHAML TEST: may modify the graphics at any time."
echo "PHAML TEST:"
echo "PHAML TEST: Are you ready to run the test (y)? "
read ANS
echo "PHAML TEST: Run test $THISDIR/test04"
if [ $PHAML_PARALLEL = "sequential" ]
then
   rm -f /tmp/phaml_message
   phaml_graphics &
fi
$RUN1 ./test04.exe
PASS=$?
if [ $PASS != 0 ]
then
   echo "PHAML TEST:"
   echo "PHAML TEST: FAILURE $THISDIR/test04 did not run."
else
   echo "PHAML TEST: SUCCESS for $THISDIR/test04."
fi

sleep 2
echo "PHAML TEST: $THISDIR/test05 tests graphics with high order solution."
echo "PHAML TEST:"
echo "PHAML TEST: In this test set 'subelement resolution' in the graphics"
echo "PHAML TEST: menu to 2 for a smoother plot."
echo "PHAML TEST: To do this, move the mouse pointer over the graphics"
echo "PHAML TEST: window, press and hold down the right button, move over"
echo "PHAML TEST: 'subelement resolution' in the menu, move over '2' in the"
echo "PHAML TEST: submenu, and release the mouse button."
echo "PHAML TEST:"
echo "PHAML TEST: Are you ready to run the test (y)? "
read ANS
echo "PHAML TEST: Run test $THISDIR/test05"
if [ $PHAML_PARALLEL = "sequential" ]
then
   rm -f /tmp/phaml_message
   phaml_graphics &
fi
$RUN1 ./test05.exe
PASS=$?
if [ $PASS != 0 ]
then
   echo "PHAML TEST:"
   echo "PHAML TEST: FAILURE $THISDIR/test05 did not run."
else
   echo "PHAML TEST: SUCCESS for $THISDIR/test05."
fi

exit 0
