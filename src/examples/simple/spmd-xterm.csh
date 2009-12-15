#!/bin/csh -f

# WFM 9/2/05
# This script is modified from a script in some version of the LAM User's
# manual.  It allows launching some processes with an xterm, possibly
# running a debugger, and other processes without.
# Run as (for LAM)
# mpirun -w -np <number of procs> spmd-xterm.csh

# Set debugger to use, or leave blank for just an xterm

#set debugger=gdb
set debugger=

# Launch the master (rank 0) with an xterm and the other processes without.
# You could list others in here to also get an xterm, e.g.
#if ("$LAMRANK" == "0" || "$LAMRANK" == "1" || "$LAMRANK" == "2") then

if ("$LAMRANK" == "0") then
   xterm -e $debugger phaml
else
   phaml
endif

exit 0
