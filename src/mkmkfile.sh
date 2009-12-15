#!/bin/bash

# Make sure you are using bash instead of a Bourne shell.
# If you need to change the path to bash, you will also need to change it in
# mkmkfile.sh in directory src and all the directories under examples.

##############################################################################
#
# Instructions for modifying this file.
#
# Modification of this file for your system(s) consists of the following
# steps.  To find the places where these changes should be made,
# search for the string "Step"
#
# Step 1.  Set your default system configuration.
# Step 2.  Set paths, library names, etc.
# Step 3.  Reset paths, library names, etc. for specific system configurations
#
##############################################################################

##############################################################################

# Step 1.
# Set default system configuration.  These can be overridden by environment
# variables of the same name without "DEFAULT_".  Command line arguments
# of the same name without "DEFAULT_" override both these defaults and
# environment variables.  Command line arguments can also omit "PHAML_".
# The variable name and value should be separated by a space, e.g.
# "./mkmkfile.sh PARLIB mpi PHAML_F90 nag".
# Acceptable values are listed before each variable.  Hopefully they are
# self explanatory.  If not, look at the section on system specific settings.
# "./mkmkfile.sh help" will list all the variables, acceptable values and
# default values.
# You can also add your own values to the existing list if you think you
# know what you are doing.  Just use one of the existing values as a guide
# for what you need to add to this script.
# Naturally, not all combinations are valid.  There is very little checking
# of the validity of the combination.

# Architecture: origin rs6k sgi sun tflop x86
DEFAULT_PHAML_ARCH=x86

# OS: aix cougar irixn32 irix64 linux solaris
DEFAULT_PHAML_OS=linux

# Fortran 90 compiler: absoft g95 gfortran intel lahey nag pathscale pgi sgi sun xlf
DEFAULT_PHAML_F90=lahey

# C compiler: cc gcc
DEFAULT_PHAML_C=cc

# Hash size: 1 2
# (2 allows more levels of refinement, but uses more memory and time)
DEFAULT_PHAML_HASHSIZE=1

# Parallelism: messpass_spawn messpass_nospawn sequential
DEFAULT_PHAML_PARALLEL=messpass_spawn

# Parallel library: lam mpi mpich mpich2 myrinet openmpi pvm none
# (mpi is used for native MPI libraries, as opposed to lam or mpich)
DEFAULT_PHAML_PARLIB=lam

# Graphics: metro mesa none opengl
DEFAULT_PHAML_GRAPHICS=mesa

# BLAS: atlas compiler goto source standard vendor
DEFAULT_PHAML_BLAS=source

# LAPACK: atlas compiler source standard vendor
DEFAULT_PHAML_LAPACK=source

# ARPACK: no yes
DEFAULT_PHAML_ARPACK=no

# BLOPEX: no withhypre withpetsc yes
DEFAULT_PHAML_BLOPEX=no

# hypre: no yes
DEFAULT_PHAML_HYPRE=no

# MUMPS: no yes
DEFAULT_PHAML_MUMPS=no

# PETSc: no yes
DEFAULT_PHAML_PETSC=no

# SuperLU: no yes
DEFAULT_PHAML_SUPERLU=no

# Zoltan: no yes
DEFAULT_PHAML_ZOLTAN=no

# ParMETIS: no yes
DEFAULT_PHAML_PARMETIS=no

# JOSTLE: no yes
DEFAULT_PHAML_JOSTLE=no

# PaToH: no yes
DEFAULT_PHAML_PATOH=no

# ParKway: no yes
DEFAULT_PHAML_PARKWAY=no

# Nemesis: no yes
DEFAULT_PHAML_NEMESIS=no

# DRUM: no yes
DEFAULT_PHAML_DRUM=no

# Specific system: none dragon looneyjr looney octopus privet raritan sgis speedy suns tflop
DEFAULT_PHAML_SYSTEM=none

# end of Step 1; proceed to Step 2 below

##############################################################################

# Set the system configuration from environment variables when they are
# set, and the default when not.

if [ -z $PHAML_ARCH ]
then
   PHAML_ARCH=$DEFAULT_PHAML_ARCH
fi
if [ -z $PHAML_OS ]
then
   PHAML_OS=$DEFAULT_PHAML_OS
fi
if [ -z $PHAML_F90 ]
then
   PHAML_F90=$DEFAULT_PHAML_F90
fi
if [ -z $PHAML_C ]
then
   PHAML_C=$DEFAULT_PHAML_C
fi
if [ -z $PHAML_HASHSIZE ]
then
   PHAML_HASHSIZE=$DEFAULT_PHAML_HASHSIZE
fi
if [ -z $PHAML_PARALLEL ]
then
   PHAML_PARALLEL=$DEFAULT_PHAML_PARALLEL
fi
if [ -z $PHAML_PARLIB ]
then
   PHAML_PARLIB=$DEFAULT_PHAML_PARLIB
fi
if [ -z $PHAML_GRAPHICS ]
then
   PHAML_GRAPHICS=$DEFAULT_PHAML_GRAPHICS
fi
if [ -z $PHAML_BLAS ]
then
   PHAML_BLAS=$DEFAULT_PHAML_BLAS
fi
if [ -z $PHAML_LAPACK ]
then
   PHAML_LAPACK=$DEFAULT_PHAML_LAPACK
fi
if [ -z $PHAML_ARPACK ]
then
   PHAML_ARPACK=$DEFAULT_PHAML_ARPACK
fi
if [ -z $PHAML_BLOPEX ]
then
   PHAML_BLOPEX=$DEFAULT_PHAML_BLOPEX
fi
if [ -z $PHAML_HYPRE ]
then
   PHAML_HYPRE=$DEFAULT_PHAML_HYPRE
fi
if [ -z $PHAML_MUMPS ]
then
   PHAML_MUMPS=$DEFAULT_PHAML_MUMPS
fi
if [ -z $PHAML_PETSC ]
then
   PHAML_PETSC=$DEFAULT_PHAML_PETSC
fi
if [ -z $PHAML_SUPERLU ]
then
   PHAML_SUPERLU=$DEFAULT_PHAML_SUPERLU
fi
if [ -z $PHAML_ZOLTAN ]
then
   PHAML_ZOLTAN=$DEFAULT_PHAML_ZOLTAN
fi
if [ -z $PHAML_PARMETIS ]
then
   PHAML_PARMETIS=$DEFAULT_PHAML_PARMETIS
fi
if [ -z $PHAML_JOSTLE ]
then
   PHAML_JOSTLE=$DEFAULT_PHAML_JOSTLE
fi
if [ -z $PHAML_PATOH ]
then
   PHAML_PATOH=$DEFAULT_PHAML_PATOH
fi
if [ -z $PHAML_PARKWAY ]
then
   PHAML_PARKWAY=$DEFAULT_PHAML_PARKWAY
fi
if [ -z $PHAML_NEMESIS ]
then
   PHAML_NEMESIS=$DEFAULT_PHAML_NEMESIS
fi
if [ -z $PHAML_DRUM ]
then
   PHAML_DRUM=$DEFAULT_PHAML_DRUM
fi
if [ -z $PHAML_SYSTEM ]
then
   PHAML_SYSTEM=$DEFAULT_PHAML_SYSTEM
fi

##############################################################################

# Change system configuration based on command line arguments

while [ $# -gt 0 ]
do
   if [ "$1" = "help" ]
   then
      echo "mkmkfile.sh [ VAR val]* where VAR, choices for val, and (current setting) are:"
      echo "  ARCH origin rs6k sgi sun tflop x86 ($PHAML_ARCH)"
      echo "  OS aix cougar irixn32 irix64 linux solaris ($PHAML_OS)"
      echo "  F90 absoft g95 gfortran intel lahey nag pathscale pgi sgi sun xlf ($PHAML_F90)"
      echo "  C cc gcc ($PHAML_C)"
      echo "  HASHSIZE 1 2 ($PHAML_HASHSIZE)"
      echo "  PARALLEL messpass_spawn messpass_nospawn sequential ($PHAML_PARALLEL)"
      echo "  PARLIB lam mpi mpich mpich2 myrinet openmpi pvm none ($PHAML_PARLIB)"
      echo "  GRAPHICS metro mesa none opengl ($PHAML_GRAPHICS)"
      echo "  BLAS atlas compiler goto source standard vendor ($PHAML_BLAS)"
      echo "  LAPACK atlas compiler source standard vendor ($PHAML_LAPACK)"
      echo "  ARPACK no yes ($PHAML_ARPACK)"
      echo "  BLOPEX no withpetsc ($PHAML_BLOPEX)"
      echo "  HYPRE no yes ($PHAML_HYPRE)"
      echo "  MUMPS no yes ($PHAML_MUMPS)"
      echo "  PETSC no yes ($PHAML_PETSC)"
      echo "  SUPERLU no yes ($PHAML_SUPERLU)"
      echo "  ZOLTAN no yes ($PHAML_ZOLTAN)"
      echo "  PARMETIS no yes ($PHAML_PARMETIS)"
      echo "  JOSTLE no yes ($PHAML_JOSTLE)"
      echo "  PATOH no yes ($PHAML_PATOH)"
      echo "  PARKWAY no yes ($PHAML_PARKWAY)"
      echo "  NEMESIS no yes ($PHAML_NEMESIS)"
      echo "  DRUM no yes ($PHAML_DRUM)"
      echo "  SYSTEM none dragon octopus privet raritan looney looneyjr speedy sgis suns tflop ($PHAML_SYSTEM)"
      exit 0
   fi
   if [ $# -eq 1 ]
   then
      echo "mkmkfile.sh: value for configuration variable $1 is missing"
      exit 1
   else
      val=$2
   fi
   case "$1" in
      PHAML_ARCH|ARCH) PHAML_ARCH=$val ;;
      PHAML_OS|OS) PHAML_OS=$val ;;
      PHAML_F90|F90) PHAML_F90=$val ;;
      PHAML_C|C) PHAML_C=$val ;;
      PHAML_HASHSIZE|HASHSIZE) PHAML_HASHSIZE=$val ;;
      PHAML_PARALLEL|PARALLEL) PHAML_PARALLEL=$val ;;
      PHAML_PARLIB|PARLIB) PHAML_PARLIB=$val ;;
      PHAML_GRAPHICS|GRAPHICS) PHAML_GRAPHICS=$val ;;
      PHAML_BLAS|BLAS) PHAML_BLAS=$val ;;
      PHAML_LAPACK|LAPACK) PHAML_LAPACK=$val ;;
      PHAML_ARPACK|ARPACK) PHAML_ARPACK=$val ;;
      PHAML_BLOPEX|BLOPEX) PHAML_BLOPEX=$val ;;
      PHAML_HYPRE|HYPRE) PHAML_HYPRE=$val ;;
      PHAML_MUMPS|MUMPS) PHAML_MUMPS=$val ;;
      PHAML_PETSC|PETSC) PHAML_PETSC=$val ;;
      PHAML_SUPERLU|SUPERLU) PHAML_SUPERLU=$val ;;
      PHAML_ZOLTAN|ZOLTAN) PHAML_ZOLTAN=$val ;;
      PHAML_PARMETIS|PARMETIS) PHAML_PARMETIS=$val ;;
      PHAML_JOSTLE|JOSTLE) PHAML_JOSTLE=$val ;;
      PHAML_PATOH|PATOH) PHAML_PATOH=$val ;;
      PHAML_PARKWAY|PARKWAY) PHAML_PARKWAY=$val ;;
      PHAML_NEMESIS|NEMESIS) PHAML_NEMESIS=$val ;;
      PHAML_DRUM|DRUM) PHAML_DRUM=$val ;;
      PHAML_SYSTEM|SYSTEM) PHAML_SYSTEM=$val ;;
      *) echo "mkmkfile.sh: argument $1 is not a system configuration variable"
         exit 1 ;;
   esac
   shift
   shift
done

##############################################################################

# Check validity of system configuration

# Check that values are in the recognized lists

case "$PHAML_ARCH" in
   origin|rs6k|sgi|sun|tflop|x86) ;;
   *) echo "mkmkfile.sh: $PHAML_ARCH is not a valid value for PHAML_ARCH"
      echo "             use one of origin rs6k sgi sun tflop x86"
      exit 1 ;;
esac
case "$PHAML_OS" in
   aix|cougar|irixn32|irix64|linux|solaris) ;;
   *) echo "mkmkfile.sh: $PHAML_OS is not a valid value for PHAML_OS"
      echo "             use one of aix cougar irixn32 irix64 linux solaris"
      exit 1 ;;
esac
case "$PHAML_F90" in
   absoft|g95|gfortran|intel|lahey|nag|pathscale|pgi|sgi|sun|xlf) ;;
   *) echo "mkmkfile.sh: $PHAML_F90 is not a valid value for PHAML_F90"
      echo "             use one of absoft g95 gfortran intel lahey nag pathscale pgi sgi sun xlf"
      exit 1 ;;
esac
case "$PHAML_C" in
   cc|gcc) ;;
   *) echo "mkmkfile.sh: $PHAML_C is not a valid value for PHAML_C"
      echo "             use one of cc gcc"
      exit 1 ;;
esac
case "$PHAML_HASHSIZE" in
   1|2) ;;
   *) echo "mkmkfile.sh: $PHAML_HASHSIZE is not a valid value for PHAML_HASHSIZE"
      echo "             use one of 1 2"
      exit 1 ;;
esac
case "$PHAML_PARALLEL" in
   messpass_spawn|messpass_nospawn|sequential) ;;
   *) echo "mkmkfile.sh: $PHAML_PARALLEL is not a valid value for PHAML_PARALLEL"
      echo "             use one of messpass_spawn messpass_nospawn sequential"
      exit 1 ;;
esac
case "$PHAML_PARLIB" in
   lam|mpi|mpich|mpich2|myrinet|openmpi|pvm|none) ;;
   *) echo "mkmkfile.sh: $PHAML_PARLIB is not a valid value for PHAML_PARLIB"
      echo "             use one of lam mpi mpich mpich2 myrinet openmpi pvm none"
      exit 1 ;;
esac
case "$PHAML_GRAPHICS" in
   metro|mesa|none|opengl) ;;
   *) echo "mkmkfile.sh: $PHAML_GRAPHICS is not a valid value for PHAML_GRAPHICS"
      echo "             use one of metro mesa none opengl"
      exit 1 ;;
esac
case "$PHAML_BLAS" in
   atlas|compiler|goto|source|standard|vendor) ;;
   *) echo "mkmkfile.sh: $PHAML_BLAS is not a valid value for PHAML_BLAS"
      echo "             use one of atlas compiler goto source standard vendor"
      exit 1 ;;
esac
case "$PHAML_LAPACK" in
   atlas|compiler|source|standard|vendor) ;;
   *) echo "mkmkfile.sh: $PHAML_LAPACK is not a valid value for PHAML_LAPACK"
      echo "             use one of atlas compiler source standard vendor"
      exit 1 ;;
esac
case "$PHAML_ARPACK" in
   no|yes) ;;
   *) echo "mkmkfile.sh: $PHAML_ARPACK is not a valid value for PHAML_ARPACK"
      echo "             use one of no yes"
      exit 1 ;;
esac
case "$PHAML_BLOPEX" in
   no|withpetsc) ;;
   *) echo "mkmkfile.sh: $PHAML_BLOPEX is not a valid value for PHAML_BLOPEX"
      echo "             use one of no withpetsc"
      exit 1 ;;
esac
case "$PHAML_HYPRE" in
   no|yes) ;;
   *) echo "mkmkfile.sh: $PHAML_HYPRE is not a valid value for PHAML_HYPRE"
      echo "             use one of no yes"
      exit 1 ;;
esac
case "$PHAML_MUMPS" in
   no|yes) ;;
   *) echo "mkmkfile.sh: $PHAML_MUMPS is not a valid value for PHAML_MUMPS"
      echo "             use one of no yes"
      exit 1 ;;
esac
case "$PHAML_PETSC" in
   no|yes) ;;
   *) echo "mkmkfile.sh: $PHAML_PETSC is not a valid value for PHAML_PETSC"
      echo "             use one of no yes"
      exit 1 ;;
esac
case "$PHAML_SUPERLU" in
   no|yes) ;;
   *) echo "mkmkfile.sh: $PHAML_SUPERLU is not a valid value for PHAML_SUPERLU"
      echo "             use one of no yes"
      exit 1 ;;
esac
case "$PHAML_ZOLTAN" in
   no|yes) ;;
   *) echo "mkmkfile.sh: $PHAML_ZOLTAN is not a valid value for PHAML_ZOLTAN"
      echo "             use one of no yes"
      exit 1 ;;
esac
case "$PHAML_PARMETIS" in
   no|yes) ;;
   *) echo "mkmkfile.sh: $PHAML_PARMETIS is not a valid value for PHAML_PARMETIS"
      echo "             use one of no yes"
      exit 1 ;;
esac
case "$PHAML_JOSTLE" in
   no|yes) ;;
   *) echo "mkmkfile.sh: $PHAML_JOSTLE is not a valid value for PHAML_JOSTLE"
      echo "             use one of no yes"
      exit 1 ;;
esac
case "$PHAML_PATOH" in
   no|yes) ;;
   *) echo "mkmkfile.sh: $PHAML_PATOH is not a valid value for PHAML_PATOH"
      echo "             use one of no yes"
      exit 1 ;;
esac
case "$PHAML_PARKWAY" in
   no|yes) ;;
   *) echo "mkmkfile.sh: $PHAML_PARKWAY is not a valid value for PHAML_PARKWAY"
      echo "             use one of no yes"
      exit 1 ;;
esac
case "$PHAML_NEMESIS" in
   no|yes) ;;
   *) echo "mkmkfile.sh: $PHAML_NEMESIS is not a valid value for PHAML_NEMESIS"
      echo "             use one of no yes"
      exit 1 ;;
esac
case "$PHAML_DRUM" in
   no|yes) ;;
   *) echo "mkmkfile.sh: $PHAML_DRUM is not a valid value for PHAML_DRUM"
      echo "             use one of no yes"
      exit 1 ;;
esac
case "$PHAML_SYSTEM" in
   none|dragon|looneyjr|looney|octopus|privet|raritan|sgis|speedy|suns|tflop) ;;
   *) echo "mkmkfile.sh: $PHAML_SYSTEM is not a valid value for PHAML_SYSTEM"
      echo "             use one of none dragon looneyjr looney octopus privet raritan sgis speedy suns tflop"
      exit 1 ;;
esac

# check for easy-to-identify incorrect combinations

if [ $PHAML_PARALLEL == "sequential" -a $PHAML_PARLIB != "none" ]
then
   echo "mkmkfile.sh: PARALLEL==sequential requires PARLIB==none"
   exit 1
fi

if [ $PHAML_PARALLEL != "sequential" -a $PHAML_PARLIB == "none" ]
then
   echo "mkmkfile.sh: if PARALLEL is not sequential then PARLIB cannot be none"
   exit 1
fi

if [ $PHAML_PARLIB == "mpich" -a $PHAML_PARALLEL != "messpass_nospawn" ]
then
   echo "mkmkfile.sh: PARLIB==mpich requires PARALLEL==messpass_nospawn"
   exit 1
fi

if [ $PHAML_MUMPS != "no" -a $PHAML_PARLIB == "pvm" ]
then
   echo "mkmkfile.sh: MUMPS does not work with PVM"
   exit 1
fi

if [ $PHAML_PARLIB != "lam" -a $PHAML_PARLIB != "mpi" -a $PHAML_PARLIB != "mpich" -a $PHAML_PARLIB != "mpich2" -a $PHAML_PARLIB != "myrinet" -a $PHAML_PARLIB != "openmpi" ]
then
   if [ $PHAML_HYPRE != "no" ]
   then
      echo "mkmkfile.sh: hypre requires an MPI library"
      exit 1
   fi
   if [ $PHAML_PETSC != "no" ]
   then
      echo "mkmkfile.sh: PetSC requires an MPI library"
      exit 1
   fi
   if [ $PHAML_SUPERLU != "no" ]
   then
      echo "mkmkfile.sh: SuperLU requires an MPI library"
      exit 1
   fi
   if [ $PHAML_ZOLTAN != "no" ]
   then
      echo "mkmkfile.sh: Zoltan requires an MPI library"
      exit 1
   fi
fi

if [ $PHAML_BLOPEX == "yes" -a $PHAML_HYPRE == "yes" ]
then
   echo "mkmkfile.sh: BLOPEX==yes requires HYPRE==no"
   exit 1
fi

if [ $PHAML_BLOPEX == "withpetsc" -a $PHAML_HYPRE == "yes" ]
then
   echo "mkmkfile.sh: BLOPEX==withpetsc requires HYPRE==no"
   exit 1
fi

if [ $PHAML_BLOPEX == "withhypre" -a $PHAML_HYPRE == "no" ]
then
   echo "mkmkfile.sh: BLOPEX==withhypre requires HYPRE==yes"
   exit 1
fi

if [ $PHAML_BLOPEX == "withpetsc" -a $PHAML_PETSC == "no" ]
then
   echo "mkmkfile.sh: BLOPEX==withpetsc requires PETSC==yes"
   exit 1
fi

if [ $PHAML_ZOLTAN == "no" -a $PHAML_PARMETIS == "yes" ]
then
   echo "mkmkfile.sh: PARMETIS==yes requires ZOLTAN==yes"
   exit 1
fi

if [ $PHAML_ZOLTAN == "no" -a $PHAML_JOSTLE == "yes" ]
then
   echo "mkmkfile.sh: JOSTLE==yes requires ZOLTAN==yes"
   exit 1
fi

if [ $PHAML_ZOLTAN == "no" -a $PHAML_PATOH == "yes" ]
then
   echo "mkmkfile.sh: PATOH==yes requires ZOLTAN==yes"
   exit 1
fi

if [ $PHAML_ZOLTAN == "no" -a $PHAML_PARKWAY == "yes" ]
then
   echo "mkmkfile.sh: PARKWAY==yes requires ZOLTAN==yes"
   exit 1
fi

if [ $PHAML_ZOLTAN == "no" -a $PHAML_NEMESIS == "yes" ]
then
   echo "mkmkfile.sh: NEMESIS==yes requires ZOLTAN==yes"
   exit 1
fi

if [ $PHAML_ZOLTAN == "no" -a $PHAML_DRUM == "yes" ]
then
   echo "mkmkfile.sh: DRUM==yes requires ZOLTAN==yes"
   exit 1
fi

##############################################################################

# Misc variables

MAKELIB="ar rcv"
RANLIB=ranlib
PHAML_HOME=`pwd`

##############################################################################

# Step 2.  Set paths, library names, etc.
# In this section, set the variables that specify the specific names of
# libraries, paths to libraries, paths to module and include files,
# command names, etc.  Of course, you only need to set those you are
# actually going to use!  Most of these variables are used directly on the
# command line, so usually should include whatever compiler flags
# accompany them.  Follow the form used in the examples.  Some of these may
# appear to be unneccessarily complicated; this is because our machines have
# multiple compilers and libraries installed which would conflict if not
# handled carefully.  Your paths are likely to be much simpler, or even
# nonexistent in many cases (where the system already knows it).  Note the use
# of single quotes when a dollar sign is to be passed through.
# Values set here can be overwritten by step 3 if different system
# configurations require different values.
# Specific version numbers quoted here indicate the last version of a
# software package that PHAML was tested with.  This is not to imply
# that earlier or later versions will not work.

case "$PHAML_F90" in

# Absoft Fortran 95 Version 10.2.1

   absoft)
      F90=af95
# starting with Version 10.1 Absoft uses the lower case and attach underscore
# name mangling.  Earlier versions need -YEXT_NAMES=LCS -YEXT_SFX=_ -YCFRL=1
      FFLAGS="-O -w -YNO_CDEC"
      LINKER=$F90
      LINKFLAGS=
      MODFLAG="-p "
      MODSUFFIX=mod
      CPUSEC=cpusec_f95.f90
      CPUSEC_COMP='$(F90) $(FFLAGS)' ;;

# g95 4.0.1

   g95)
      F90=g95
      FFLAGS="-O -fno-second-underscore"
      LINKER=$F90
      LINKFLAGS=
      MODFLAG=-I
      MODSUFFIX=mod
      CPUSEC=cpusec_f95.f90
      CPUSEC_COMP='$(F90) $(FFLAGS)' ;;

# GNU Fortran 95 4.3.0 20080124

   gfortran)
      F90=gfortran
      FFLAGS="-O -fno-second-underscore"
      LINKER=$F90
      LINKFLAGS=
      MODFLAG=-I
      MODSUFFIX=mod
      CPUSEC=cpusec_f95.f90
      CPUSEC_COMP='$(F90) $(FFLAGS)' ;;

# Intel Fortran 11.1.046

   intel)
      F90=ifort
# as of 11.1 don't need to disable remark, but some of my systems are older
      FFLAGS="-O -w -diag-disable remark"
      LINKER=$F90
      LINKFLAGS=
      MODFLAG=-I
      MODSUFFIX=mod
      CPUSEC=cpusec_f95.f90
      CPUSEC_COMP='$(F90) $(FFLAGS)' ;;

# Lahey lf95 L6.20d

   lahey)
      F90=lf95
      FFLAGS=-O
      LINKER=$F90
      LINKFLAGS=--staticlink
      MODFLAG=-I
      MODSUFFIX=mod
      CPUSEC=cpusec_f95.f90
      CPUSEC_COMP='$(F90) $(FFLAGS)' ;;

# NAGWare Fortran 95 Release 5.2(662)

   nag)
      F90=nagfor
      FFLAGS='-O -w'
      if [ $PHAML_BLAS = "source" -o $PHAML_LAPACK = "source" ]
      then
         FFLAGS='-O -w -dusty -dcfuns'
      fi
      LINKER=$F90
      LINKFLAGS=
      MODFLAG=-I
      MODSUFFIX=mod
      CPUSEC=cpusec_f95.f90
      CPUSEC_COMP='$(F90) $(FFLAGS)' ;;

# PathScale EKOPath Fortran 95 Version 2.1, sold by Lahey

   pathscale)
      F90=pathf90
      FFLAGS=-O
      LINKER=$F90
      LINKFLAGS=
      MODFLAG=-I
      MODSUFFIX=mod
      CPUSEC=cpusec_f95.f90
      CPUSEC_COMP='$(F90) $(FFLAGS)' ;;

# PGI pgf90 6.0-2

   pgi)
      F90=pgf90
      FFLAGS=-O
      LINKER=$F90
      LINKFLAGS=
      MODFLAG=-I
      MODSUFFIX=mod
      CPUSEC=cpusec_f95.f90
      CPUSEC_COMP='$(F90) $(FFLAGS)' ;;

# SGI MIPSpro f90 7.4.4m

   sgi)
      F90=f90
      case "$PHAML_OS" in
         irixn32)
            FFLAGS="-O -n32"
            LINKFLAGS=-n32 ;;
         irix64)
            FFLAGS="-O -64"
            LINKFLAGS=-64 ;;
      esac
      LINKER=$F90
      MODFLAG=-I
      MODSUFFIX=mod
      CPUSEC=cpusec.etime.f
      CPUSEC_COMP='$(F90) $(FFLAGS)' ;;

# Sun Forte Developer 7 Fortran 95 7.0 2002/03/09

   sun)
      F90=f90
# This compiler crashes while compiling grid.f90 if -O is used
      FFLAGS=
      LINKER=$F90
      LINKFLAGS=
      MODFLAG=-M
      MODSUFFIX=mod
      CPUSEC=cpusec.etime.f
      CPUSEC_COMP='$(F90) $(FFLAGS)' ;;

# IBM XLF90 4.1.0.6

   xlf)
      F90=xlf90
      FFLAGS="-O -F:f90"
      LINKER=$F90
      LINKFLAGS=
      MODFLAG=-I
      MODSUFFIX=mod
      CPUSEC=cpusec.c1.c
      CPUSEC_COMP='$(CC) $(CFLAGS)' ;;

# invalid value for PHAML_F90

   *) echo "mkmkfile.sh: $PHAML_F90 is not a valid value for PHAML_F90"
      exit 1 ;;

esac

case "$PHAML_C" in

# native C compiler

   cc)
      CC=cc
      case "$PHAML_OS" in
         irix64)
            CFLAGS="-O -64" ;;
         *)
            CFLAGS=-O ;;
      esac ;;

# GCC 4.3.0 20080124

   gcc)
      CC=gcc
      CFLAGS=-O ;;

# invalid value for PHAML_C

   *) echo "mkmkfile.sh: $PHAML_C is not a valid value for PHAML_C"
      exit 1 ;;

esac

case "$PHAML_PARLIB" in

# LAM 7.1.4

   lam)
      RUNMPI="mpirun"
      MESSPASSINC=
      MESSPASSLIBS=
      MESSPASSPACK='mpipack1.o mpipack2.o'
      MPIF=$F90
      F90="mpif77"
      LINKER="mpif77"
      CC=mpicc ;;

# native MPI libraries; probably want to override in the system specific section

   mpi)
      RUNMPI="mpiexec"
      MESSPASSINC=""
      MESSPASSLIBS="-lmpi"
      MESSPASSPACK='mpipack1.o mpipack2.o' ;;

# MPICH 1.2.7

   mpich)
      RUNMPI="mpirun"
      MESSPASSINC='-I$(MPICH_HOME)'"/$PHAML_F90/include"
      MESSPASSLIBS='-L$(MPICH_HOME)'"/$PHAML_F90/lib -lmpich"
      MESSPASSPACK='mpipack1.o mpipack2.o' ;;

# MPICH2 1.1.1p1

   mpich2)
      RUNMPI="mpiexec"
      MESSPASSINC=
      MESSPASSLIBS=
      MESSPASSPACK='mpipack1.o mpipack2.o'
      F90="mpif90"
      LINKER="mpif90"
      CC=mpicc ;;

# Myrinet MPICH 1.2.4

   myrinet)
      RUNMPI="mpirun"
      MESSPASSINC='-I$(MPICH_HOME)/include'
      MESSPASSLIBS='-L$(MPICH_HOME)/lib -lpmpich -lmpich -L/home/pozo/MYRICOM/gm-1.2.3/binary/lib -lgm'
      MESSPASSPACK='mpipack1.o mpipack2.o' ;;

# Open MPI 1.1.2

   openmpi)

      RUNMPI="mpirun"
      MESSPASSINC=
      MESSPASSLIBS=
      MESSPASSPACK='mpipack1.o mpipack2.o'
      MPIF=$F90
      F90="mpif90"
      LINKER="mpif90"
      CC=mpicc ;;

# PVM 3.4.5

   pvm)
      RUNMPI=""
      MESSPASSINC='-I$(PVM_ROOT)/include'
      MESSPASSLIBS='-L$(PVM_ROOT)/lib/$(PVM_ARCH) -lfpvm3 -lpvm3'
      MESSPASSPACK='pvmpack1.o pvmpack2.o' ;;

# no message passing library

   none)
      RUNMPI="./"
      MESSPASSINC=""
      MESSPASSLIBS=""
      MESSPASSPACK="" ;;

# invalid value for PHAML_PARLIB

   *) echo "mkmkfile.sh: $PHAML_PARLIB is not a valid value for PHAML_PARLIB"
      exit 1 ;;

esac

# X11R6

XINCL=-I/usr/X11R6/include
XLIBS="-L/usr/X11R6/lib -lXaw -lXt -lXmu -lXi -lX11 -lXext -lm"

case "$PHAML_GRAPHICS" in

# MetroLink OpenGL 1.3, Glut 3.7.1, f90gl 1.2.0

   metro)
      OGLMODS="$MODFLAG"'$(OPENGL_HOME)/include/GL'
      OGLLIBS='-L$(OPENGL_HOME)/lib -L/usr/X11R6/lib -lf90glut -lf90GLU -lf90GL -lglut -lGLU -lGL' ;;

# Mesa 7.0.3, Glut MesaGLUT 7.0.3, f90gl 1.2.13

   mesa)
      OGLMODS="$MODFLAG"'$(MESA_HOME)/'"$PHAML_F90/include/GL"
      OGLLIBS='-L$(MESA_HOME)/lib -L$(MESA_HOME)/'"$PHAML_F90/lib -lf90glut -lf90GLU -lf90GL -lglut -lGLU -lGL" ;;

# No graphics library

   none)
      OGLMODS=
      OGLLIBS= ;;

# Native OpenGL, may need to redefine in the system specific section

   opengl)
      OGLMODS="$MODFLAG"'$(OPENGL_HOME)/include/GL'
      OGLLIBS='-L$(OPENGL_HOME)/lib -L/usr/X11R6/lib -lf90glut -lf90GLU -lf90GL -lglut -lGLU -lGL' ;;

# invalid value for PHAML_GRAPHICS

   *) echo "mkmkfile.sh: $PHAML_GRAPHICS is not a valid value for PHAML_GRAPHICS"
      exit 1 ;;

esac

# Zoltan 3.0

ZOLTANLIB='-L$(ZOLTAN_HOME)'"/$PHAML_F90/$PHAML_PARLIB/lib -lzoltan -lzoltan_comm -lzoltan_mem"
ZOLTANMOD="$MODFLAG"'$(ZOLTAN_HOME)'"/$PHAML_F90/$PHAML_PARLIB/include"
ZOLTANINC='-I$(ZOLTAN_HOME)/include -I$(ZOLTAN_HOME)/Utilities/Communication -I$(ZOLTAN_HOME)/Utilities/Memory -I$(ZOLTAN_HOME)/Utilities/DDirectory'

# ParMETIS 3.1.0

PARMETISLIB='-L$(PARMETIS_HOME)'"/$PHAML_PARLIB/lib -lparmetis -lmetis"

# Jostle untried

JOSTLELIB='-L$(JOSTLE_HOME) -ljostle'

# PaToH untried

PATOHLIB='-L$(PATOH_HOME) -lpatoh'

# ParKway untried

PARKWAYLIB='-L$(PARKWAY_HOME) -lparkway'

# Nemesis, Exodus and netcdf untried

NEMESISLIB='-L$(NEMESIS_HOME) -lnemIc -lexoIIv2c -lnetcdf'

# DRUM 2.0

DRUMLIB='-L$(DRUM_HOME)/lib/speedy -ldrum -lrt -lxml2'

case "$PHAML_ARPACK" in

# ARPACK no version number

# Since ARPACK does not support PVM, if PHAML_PARLIB is PVM, leave off the
# parallel library (and the program should specify nproc=1).

   yes)
      if [ $PHAML_PARALLEL = "sequential" -o $PHAML_PARLIB = "pvm" ]
      then
         ARPACKLIBS='-L$(ARPACK_HOME)/'"$PHAML_F90/$PHAML_PARLIB/lib -larpack"
      else
         ARPACKLIBS='-L$(ARPACK_HOME)/'"$PHAML_F90/$PHAML_PARLIB/lib -larpack -lparpack"
      fi ;;

   no)
      ARPACKLIBS= ;;

# invalid value for PHAML_ARPACK

   *) echo "mkmkfile.sh: $PHAML_ARPACK is not a valid value for PHAML_ARPACK"
      exit 1 ;;

esac

case "$PHAML_BLOPEX" in

# BLOPEX no version number

# stand alone BLOPEX

   yes)
      echo "mkmkfile.sh: BLOPEX is currently only supported via PETSc"
      exit 1
      BLOPEXINC=
      BLOPEXLIBS= ;;

# BLOPEX as an external package to PETSc.  The BLOPEX include files lobpcg.h,
# interpreter.h, multivector.h, and petsc-interface.h can be found
# in $PETSC_DIR/externalpackages/blopex_abstract/[krylov,multivector]
# and $PETSC_DIR/src/contrib/blopex/petsc-interface.  The BLOPEX library is
# also under blopex_abstract somewhere.  I relocate all these things.

   withpetsc)
      BLOPEXINC='-I$(PETSC_HOME)/blopex/include'
      BLOPEXLIBS='-L$(PETSC_HOME)/blopex/lib -lBLOPEX' ;;

# included with hypre

   withhypre)
      echo "mkmkfile.sh: BLOPEX is currently only supported via PETSc"
      exit 1
      BLOPEXINC=
      BLOPEXLIBS= ;;

   no)
      BLOPEXINC=
      BLOPEXLIBS= ;;

   *) echo "mkmkfile.sh: $PHAML_BLOPEX is not a valid value for PHAML_BLOPEX"
      exit 1 ;;

esac

case "$PHAML_BLAS" in

# ATLAS 3.6.0

   atlas)
      BLASLIBS='-L$(ATLAS_HOME) -lf77blas -lcblas -latlas' ;;

# some compilers come with prebuilt BLAS libraries

   compiler)
      case "$PHAML_F90" in
         pgi)
            BLASLIBS='-L$(PGI_HOME)/lib -lblas' ;;
         lahey)
            BLASLIBS='-L$(LAHEY_HOME)/lib -lblasmt' ;;
         intel)
            BLASLIBS='-L$(MKL_HOME) -lmkl_intel -lmkl_intel_thread -lmkl_core -liomp5' ;;
         *)
            BLASLIBS='-L/usr/lib -lblas' ;;
      esac ;;

   goto)
      BLASLIBS='-L$(GOTO_HOME) -lgoto -lpthread' ;;

# use source code supplied with PHAML instead of a precompiled library

   source)
      BLASLIBS= ;;

# use the standard BLAS library for this system

   standard)
      BLASLIBS='-L/usr/lib -lblas' ;;

# some hardware vendors supply BLAS libraries

   vendor)
      case "$PHAML_ARCH" in
         x86)
# MKL 11.1.046 supports intel and gfortran; otherwise use source
            if [ $PHAML_F90 = "intel" ]
            then
               BLASLIBS='-L$(MKL_HOME) -lmkl_intel -lmkl_intel_thread -lmkl_core -liomp5'
            else
               if [ $PHAML_F90 = "gfortran" ]
               then
                  BLASLIBS='-L$(MKL_HOME) -lmkl_gf -lmkl_gnu_thread -lmkl_core -liomp5'
               else
                  BLASLIBS=
               fi
            fi ;;

         sgi)
            BLASLIBS='-lscs' ;;
         sun)
            BLASLIBS='-xlic_lib=sunperf' ;;
         *)
            BLASLIBS='-L/usr/lib -lblas' ;;
      esac ;;

# invalid value for PHAML_BLAS

   *) echo "mkmkfile.sh: $PHAML_BLAS is not a valid value for PHAML_BLAS"
      exit 1 ;;

esac

case "$PHAML_LAPACK" in

# ATLAS 3.6.0

   atlas)
      LAPACKLIBS='-L$(ATLAS_HOME) -llapack' ;;

# some compilers come with prebuilt LAPACK libraries

   compiler)
      case "$PHAML_F90" in
         pgi)
            LAPACKLIBS='-L$(PGI_HOME)/lib -llapack' ;;
         lahey)
            LAPACKLIBS='-L$(LAHEY_HOME)/lib -llapackmt' ;;
         intel)
# Don't need them if BLAS already put them in
            if [ $PHAML_BLAS = "compiler" -o $PHAML_BLAS = "vendor" ]
            then
               LAPACKLIBS=
            else
               LAPACKLIBS='-L$(MKL_HOME) -lmkl_intel -lmkl_intel_thread -lmkl_core -liomp5'
            fi ;;
         *)
            LAPACKLIBS='-L/usr/lib -llapack' ;;
      esac ;;

# use source code supplied with PHAML instead of a precompiled library

   source)
      LAPACKLIBS= ;;

# use the standard LAPACK library for this system

   standard)
      LAPACKLIBS='-L/usr/lib -llapack' ;;

# some hardware vendors provide LAPACK libraries

   vendor)
      case "$PHAML_ARCH" in
         x86)
            if [ $PHAML_BLAS = "compiler" -o $PHAML_BLAS = "vendor" ]
            then
               LAPACKLIBS=
            else
               if [ $PHAML_F90 = "intel" ]
               then
                  LAPACKLIBS='-L$(MKL_HOME) -lmkl_intel -lmkl_intel_thread -lmkl_core -liomp5'
               else
                  if [ $PHAML_F90 = "gfortran" ]
                  then
                     LAPACKLIBS='-L$(MKL_HOME) -lmkl_gf -lmkl_gnu_thread -lmkl_core -liomp5'
                  else
                     LAPACKLIBS=
                  fi
               fi
            fi ;;
         sgi)
            LAPACKLIBS='-lscs' ;;
         sun)
            LAPACKLIBS='-xlic_lib=sunperf' ;;
         *)
            LAPACKLIBS='-L/usr/lib -llapack' ;;
      esac ;;

# invalid value for PHAML_LAPACK

   *) echo "mkmkfile.sh: $PHAML_LAPACK is not a valid value for PHAML_LAPACK"
      exit 1 ;;

esac

# PETSc 2.3.3 (Earlier versions may need changes in petsc_init.F90 and
#              petsc_interf.F90.  Search for "before" to find them.  Also,
#              select the right one for PETSCLIBS right here.)

case "$PHAML_PETSC" in

   yes)
# PETSc versions before 2.3.1
#      PETSCLIBS='-L$(PETSC_HOME)/'"$PHAML_F90/$PHAML_PARLIB"'/lib/$(PETSC_ARCH) -lpetscfortran -lpetscksp -lpetscts -lpetscsnes -lpetscdm -lpetscmat -lpetscvec -lpetsc'

# PETSc version 2.3.1 and later
      PETSCLIBS='-L$(PETSC_HOME)/'"$PHAML_F90/$PHAML_PARLIB"'/lib/$(PETSC_ARCH) -lpetscksp -lpetscts -lpetscsnes -lpetscdm -lpetscmat -lpetscvec -lpetsc'

      PETSCINC='-I$(PETSC_HOME)/'"$PHAML_F90/$PHAML_PARLIB"'/include -I$(PETSC_HOME)/'"$PHAML_F90/$PHAML_PARLIB"'/bmake/$(PETSC_ARCH) -I$(PETSC_HOME)/'"$PHAML_F90/$PHAML_PARLIB $MESSPASSINC" ;;

   no)
      PETSCLIBS=
      PETSCINC= ;;

   *) echo "mkmkfile.sh: $PHAML_PETSC is not a valid value for PHAML_PETSC"
      exit 1 ;;

esac

# hypre 2.0.0

case "$PHAML_HYPRE" in

   yes)

# use this for hypre version 2.0.0
      HYPRELIBS='-L$(HYPRE_HOME)/'"$PHAML_F90/$PHAML_PARLIB/lib -lHYPRE" ;;

# use this for hypre version 1.9.0b
#      HYPRELIBS='-L$(HYPRE_HOME)/'"$PHAML_F90/$PHAML_PARLIB/lib -lHYPRE_parcsr_ls -lHYPRE_DistributedMatrixPilutSolver -lHYPRE_ParaSails -lHYPRE_Euclid -lHYPRE_MatrixMatrix -lHYPRE_DistributedMatrix -lHYPRE_IJ_mv -lHYPRE_parcsr_mv -lHYPRE_seq_mv -lHYPRE_krylov -lHYPRE_utilities" ;;

# use this for hypre version 1.6.0
#      HYPRELIBS='-L$(HYPRE_HOME)/'"$PHAML_F90/$PHAML_PARLIB/lib -lHYPRE_parcsr_ls -lHYPRE_DistributedMatrixPilutSolver -lHYPRE_ParaSails -lHYPRE_Euclid -lHYPRE_MatrixMatrix -lHYPRE_DistributedMatrix -lHYPRE_IJ_mv -lHYPRE_parcsr_mv -lHYPRE_seq_mv -lkrylov -lHYPRE_utilities" ;;

   no)

      HYPRELIBS= ;;

   *) echo "mkmkfile.sh: $PHAML_HYPRE is not a valid value for PHAML_HYPRE"
      exit 1 ;;

esac

# BLACS 1.1 patch 3 (needed for MUMPS)

BLACSLIB='-L$(BLACS_HOME)/'"$PHAML_F90/$PHAML_PARLIB/lib -lblacs -lblacsF77init -lblacs"

# SCALAPACK 1.8 (needed for MUMPS)

SCALAPACKLIB='-L$(SCALAPACK_HOME)/'"$PHAML_F90/$PHAML_PARLIB/lib -lscalapack"

# MUMPS 4.9

case "$PHAML_MUMPS" in

   yes)

# sometime between versions 4.7.3 and 4.9 they added libmumps_common.  Pick
# the right form.
      MUMPSLIB='-L$(MUMPS_HOME)/'"$PHAML_F90/$PHAML_PARLIB/lib -ldmumps -lpord"
#      MUMPSLIB='-L$(MUMPS_HOME)/'"$PHAML_F90/$PHAML_PARLIB/lib -ldmumps -lmumps_common -lpord"
      MUMPSLIBS="$MUMPSLIB $SCALAPACKLIB $BLACSLIB"
      MUMPSINC='-I$(MUMPS_HOME)/'"$PHAML_F90/$PHAML_PARLIB/include" ;;

   no)
      MUMPSLIBS=
      MUMPSINC= ;;

   *) echo "mkmkfile.sh: $PHAML_MUMPS is not a valid value for PHAML_MUMPS"
      exit 1 ;;

esac

# SuperLU 2.0

case "$PHAML_SUPERLU" in

   yes)

      SUPERLULIBS='-L$(SUPERLU_HOME)/'"$PHAML_F90/$PHAML_PARLIB/lib -lsuperlu"
      SUPERLUINC='-I$(SUPERLU_HOME)/'"$PHAML_F90/$PHAML_PARLIB/include" ;;

   no)
      SUPERLULIBS=
      SUPERLUINC= ;;

   *) echo "mkmkfile.sh: $PHAML_SUPERLU is not a valid value for PHAML_SUPERLU"
      exit 1 ;;

esac

# If there are any other libraries that you know you need to link in,
# define them here (or in step 3).

OTHERLIBS=

# end of Step 2

##############################################################################

# Step 3.
# Here you can override the values set in step 2 for particular system
# configurations, for example if a particular library is located in
# different directories on different machines.  The variable
# PHAML_SYSTEM determines which section of changes are used.  You should
# make a new section for each of your systems where you need to override
# the values of step 2.  Follow the examples below.  The intention is that
# the value for PHAML_SYSTEM be a machine name, but in reality you can use
# whatever character string you want.  You probably want to add that string
# to the help line, too.  Search for speedy to see where you might want
# to put your system names.

case "$PHAML_SYSTEM" in

########
# raritan is a frontend to a Linux cluster

   raritan)

# raritan has a special environment setup

      case "$PHAML_PARLIB" in

         lam)

            MESSPASSINC='-I$(PREFERRED_INCDIR)'
            MESSPASSLIBS='-L$(PREFERRED_LIBDIR) -lmpi -lpthread'
            MESSPASSPACK='mpipack1.o mpipack2.o'
            LINKER='hf77'
# for lam/pgi, hf77 doesn't find the pgi libraries.  Don't use hf77 and
# explicitly link whatever libraries it needs for lam
            if [ $PHAML_F90 = "pgi" ]
            then
               MESSPASSLIBS='-L$(PREFERRED_LIBDIR) -llamf77mpi -lmpi -llam -lpthread -ldl'
               LINKER='pgf90'
            fi
            if [ $PHAML_MUMPS = "yes" ]
            then
               BLACSLIB='-L$(PREFERRED_LIBDIR) -lblacs -lblacsF77init -lblacs'
               SCALAPACKLIB='-L$(PREFERRED_LIBDIR) -lscalapack'
# the intel 64 lam blacs libraries have double underscore
# the intel 64 lam scalapack libraries have unresolved references
#            when I use my own blacs
               if [ $PREFERRED_ARCH = "64" ]
               then
                  BLACSLIB='-L$(BLACS_HOME)/intel/lam/lib -lblacs -lblacsF77init -lblacs'
                  SCALAPACKLIB='-L$(SCALAPACK_HOME)/intel/lam/lib -lscalapack'
               fi
               MUMPSLIBS="$MUMPSLIB $SCALAPACKLIB $BLACSLIB"
            fi ;;

         mpich)
            HOLD_F90="$F90"
            F90=mpif90
            LINKER=mpif90
            MESSPASSINC='-I$(PREFERRED_INCDIR)'
            MESSPASSLIBS=
            MESSPASSPACK='mpipack1.o mpipack2.o'
            if [ $PHAML_F90 = "lahey" ]
            then
               MESSPASSLIBS='-L$(PREFERRED_LIBDIR) -lfmpich -lmpichfsup'
            fi
# sse libraries (unresolved references in libphaml)
            if [ $PHAML_F90 = "pgi" ]
            then
               OTHERLIBS="$OTHERLIBS -lpgsse1 -lpgsse2"
            fi
            if [ $PHAML_MUMPS = "yes" ]
            then
               BLACSLIB='-L$(PREFERRED_LIBDIR) -lblacs -lblacsF77init -lblacs'
               SCALAPACKLIB='-L$(PREFERRED_LIBDIR) -lscalapack'
               MUMPSLIBS="$MUMPSLIB $SCALAPACKLIB $BLACSLIB"
            fi ;;

         openmpi)
            F90=mpif90
            LINKER=mpif90
            MESSPASSINC=
            MESSPASSLIBS= ;;

         mpich2)
            ;;

         *)
            if [ $PHAML_MUMPS = "yes" ]
            then
               BLACSLIB='-L$(PREFERRED_LIBDIR) -lblacs -lblacsF77init -lblacs'
               SCALAPACKLIB='-L$(PREFERRED_LIBDIR) -lscalapack'
               MUMPSLIBS="$MUMPSLIB $SCALAPACKLIB $BLACSLIB"
            fi ;;

      esac

# PETSc needs the new MESSPASSINC

      if [ $PHAML_PETSC = "yes" ]
      then
         PETSCINC='-I$(PETSC_HOME)/'"$PHAML_F90/$PHAML_PARLIB"'/include -I$(PETSC_HOME)/'"$PHAML_F90/$PHAML_PARLIB"'/bmake/$(PETSC_ARCH) -I$(PETSC_HOME)/'"$PHAML_F90/$PHAML_PARLIB $MESSPASSINC"
      fi

# Use the MPI library compiler commands for OpenMPI

      if [ $PHAML_PARLIB = "openmpi" ]
      then
         HOLD_F90="$F90"
         F90=mpif90
         LINKER=mpif90
         HOLD_CC="$CC"
         CC=mpicc
      fi

# if 64 bit chosen, use 64 bit X11, blas and lapack libraries
# There are other blas and lapack in the same directory as LAM which gets
# listed first, so use the libraries explicitly

      if [ $PREFERRED_ARCH = "64" ]
      then 
         XLIBS="-L/usr/X11R6/lib64 -lXaw -lXt -lXmu -lXi -lX11 -lXext -lm"
         if [ $PHAML_BLAS = "standard" ]
         then
            BLASLIBS='/usr/lib64/libblas.a'
         fi
         if [ $PHAML_LAPACK = "standard" ]
         then
            LAPACKLIBS='/usr/lib64/liblapack.a'
         fi
         if [ $PHAML_BLAS = "standard" -o $PHAML_LAPACK = "standard" ]
         then
            OTHERLIBS="$OTHERLIBS -L/usr/lib/gcc/x86_64-redhat-linux/3.4.3 -lg2c"
         fi
      fi

# Libraries that, at least for now, need g2c:
# atlas lapack, mkl

      if [ $PHAML_LAPACK = "atlas" -o $PHAML_BLAS = "vendor" -o $PHAML_LAPACK = "vendor" ]
      then
         OTHERLIBS="$OTHERLIBS -L/usr/lib/gcc/i386-redhat-linux/3.4.3 -lg2c"
      fi ;;

########
# speedy is a Linux laptop PC

   speedy)

# ARPACK needs etime, which can be obtained from libg2c

      if [ $PHAML_ARPACK = "yes" ]
      then
         OTHERLIBS="$OTHERLIBS -L/usr/lib/gcc/i386-redhat-linux/3.4.6 -lg2c"
      fi

# Changes if using the system LAM libraries instead of my own.
#
#      if [ $PHAML_PARLIB = "lam" ]
#      then
#         if [ $USE_MY_LAM == 0 ]
#         then
#            MESSPASSINC=
#            MESSPASSLIBS=
#            F90=mpif77
#            LINKER=mpif77
#            CC=mpicc
#         fi
#      fi

# using MUMPS 4.9 that needs mumps_common

      if [ $PHAML_MUMPS = "yes" ]
      then
         MUMPSLIB='-L$(MUMPS_HOME)/'"$PHAML_F90/$PHAML_PARLIB/lib -ldmumps -lmumps_common -lpord"
         MUMPSLIBS="$MUMPSLIB $SCALAPACKLIB $BLACSLIB"
      fi

# The Intel ifort compilation of PETSc on speedy requires _gfortran_getarg_i4
# and _gfortran_iargc from libgfortran, and dstebz_ and dstein_ from Lapack
# which are not in my lapack source.  -lgfortran automatically comes along
# with -llapack because /usr/bin/liblapack uses something from libgfortran.

      if [ $PHAML_F90 = "intel" -a $PHAML_PETSC = "yes" ]
      then
         if [ $PHAML_LAPACK = "source" ]
         then
            OTHERLIBS="$OTHERLIBS -llapack"
         else
            if [ $PHAML_LAPACK != "standard" ]
            then
               OTHERLIBS="$OTHERLIBS -L/local/apps/gfortran/irun/lib -lgfortran"
            fi
         fi
      fi ;;

########
# octopus is an 8 dual-core processor 64 bit Linux machine

   octopus)

# Mesa calls it lib64 on linux-x86-64

      OGLLIBS='-L$(MESA_HOME)/lib64 -L$(MESA_HOME)/'"$PHAML_F90/lib -lf90glut -lf90GLU -lf90GL -lglut -lGLU -lGL" ;;

########
# privet is a Linux PC with an Intel dual core processor

#   privet)

########
# looneyjr is a Linux laptop PC with several Fortran compilers

   looneyjr)

# The standard LAPACK library (and hence also ATLAS) needs some things from g2c.
# ARPACK needs etime, which can be obtained from libg2c
# Absoft F90 does not provide system, but it is in libg2c
# Also, if LAPACK is "compiler" and the compiler is not PGI or Lahey, the
# standard Lapack library is used.

      if [ $PHAML_LAPACK = "standard" -o $PHAML_LAPACK = "atlas" -o $PHAML_ARPACK = "yes" -o $PHAML_F90 = "absoft" ]
      then
         OTHERLIBS="$OTHERLIBS -L/usr/lib/gcc/i386-redhat-linux/3.4.6 -lg2c"
      else
         if [ $PHAML_LAPACK = "compiler" ]
         then
            if [ $PHAML_F90 != "lahey" -a $PHAML_F90 != "pgi" ]
            then
               OTHERLIBS="$OTHERLIBS -L/usr/lib/gcc/i386-redhat-linux/3.4.6 -lg2c"
            fi
         fi
      fi

# Use the system LAM libraries instead of my own.  Some compilers don't
# like some of the options that lamf77 puts on.  For those, use the
# list that showme gives with the offending options omitted.  This could
# change if LAM changes.

#      if [ $PHAML_PARLIB = "lam" ]
#      then
#         MESSPASSINC=
#         MESSPASSLIBS=
#         case "$PHAML_F90" in
#            absoft)
#               FFLAGS="$FFLAGS -I/usr/include/lam -I/usr/include/lam/32 -pthread" ;;
#            intel)
#               FFLAGS="$FFLAGS -I/usr/include/lam -I/usr/include/lam/32 -pthread" ;;
#            lahey)
#               FFLAGS="$FFLAGS -I/usr/include/lam -I/usr/include/lam/32 -m32" ;;
#            *)
#               FFLAGS="$FFLAGS `lamf77-32 -showme:compile`" ;;
#         esac
#         CFLAGS="$CFLAGS `lamcc-32 -showme:compile`"
#         LINKER=$F90
#         LINKFLAGS="$LINKFLAGS `lamf77-32 -showme:link`"
#         RUNMPI="lamrun"
#      fi

# use the system Open MPI instead of my own

#      if [ $PHAML_PARLIB = "openmpi" ]
#      then
#         MESSPASSINC=
#         MESSPASSLIBS=
#         case "$PHAML_F90" in
#            absoft)
#               FFLAGS="$FFLAGS -I/usr/include/openmpi -I/usr/include/openmpi/32 -pthread"
#               LINKFLAGS="$LINKFLAGS -pthread -L/usr/lib/openmpi -lmpi -lorte -lopal -ldl -Wl,--export-dynamic -lnsl -lutil -lm -ldl" ;;
#            intel)
#               FFLAGS="$FFLAGS -I/usr/include/openmpi -I/usr/include/openmpi/32 -pthread"
#               LINKFLAGS="$LINKFLAGS -pthread -L/usr/lib/openmpi -lmpi -lorte -lopal -ldl -Wl,--export-dynamic -lnsl -lutil -lm -ldl" ;;
#            lahey)
#               FFLAGS="$FFLAGS -I/usr/include/openmpi -I/usr/include/openmpi/32 -m32"
#               LINKFLAGS="$LINKFLAGS -m32 -L/usr/lib/openmpi -lmpi -lorte -lopal -ldl -Wl,--export-dynamic -lnsl -lutil -lm -ldl" ;;
#            *)
#               FFLAGS="$FFLAGS `mpif77 -showme:compile`"
#               LINKFLAGS="$LINKFLAGS `mpif77 -showme:link`" ;;
#         esac
#         CC=$HOLD_CC
#         CFLAGS="$CFLAGS `mpicc -showme:compile`"
#         F90=$HOLD_F90
#         LINKER=$F90
#      fi

# Absoft compiler doesn't link pthread by default

      if [ $PHAML_F90 = "absoft" ]
      then
         OTHERLIBS="$OTHERLIBS -lpthread"
      fi

# Absoft doesn't include it's own libU77.a?  That's crazy!

      if [ $PHAML_F90 = "absoft" -a $PHAML_PARLIB = "mpich" ]
      then
         OTHERLIBS="$OTHERLIBS -L/opt/absoft/lib -lU77"
      fi ;;

########
# looney is a 64 bit Linux PC

   looney)

# We need g2c to get subroutine system.  It is also needed by several other
# libraries.

      OTHERLIBS="$OTHERLIBS -L/usr/lib/gcc/x86_64-redhat-linux/3.4.6 -lg2c"

# If I'm not using gfortran but I use the system blas or lapack, I need to
# link in libgfortran

      if [ $PHAML_F90 != "gfortran" ]
      then
         if [ $PHAML_BLAS = "standard" -o $PHAML_LAPACK = "standard" ]
         then
            OTHERLIBS="$OTHERLIBS -lgfortran"
         fi
      fi

# using MUMPS 4.9 that needs mumps_common, unless it's the NAG compiler and
# either LAM or OPENMPI where it's an older version because of a SIGSEGV

      if [ $PHAML_MUMPS = "yes" ]
      then
         if [[ $PHAML_F90 != "nag" || ( $PHAML_PARLIB != "lam" && $PHAML_PARLIB != "openmpi" ) ]]
         then
            MUMPSLIB='-L$(MUMPS_HOME)/'"$PHAML_F90/$PHAML_PARLIB/lib -ldmumps -lmumps_common -lpord"
            MUMPSLIBS="$MUMPSLIB $SCALAPACKLIB $BLACSLIB"
         fi
      fi

# Mesa calls it lib64 on linux-x86-64

      OGLLIBS='-L$(MESA_HOME)/lib64 -L$(MESA_HOME)/'"$PHAML_F90/lib -lf90glut -lf90GLU -lf90GL -lglut -lGLU -lGL" ;;

########
# dragon is a Sun workstation

   dragon)

      MESSPASSINC="-I/Net/proj/zoltan/arch/solaris/include \
                   -I/export/home/mpich/mpich-1.2.4/include"
      MESSPASSLIBS="-L/export/home/mpich/mpich-1.2.4/lib \
                    -L/Net/local/gnu/lib/gcc-lib/sparc-sun-solaris2.7/3.0.3 \
                    -L/Net/local/lang/SUNWspro/lib \
                    -lmpich -lsocket -lnsl -laio -lgcc -lrt"
      PARMETISLIB="-L/home/u/kddevin/code/ParMETIS/ParMETIS2 -lparmetis -lmetis"
      JOSTLELIB="-L/Net/local/proj/zoltan/arch/solaris/lib -ljostle" ;;

########
# sgis is for our SGI workstations

   sgis)

# MESA has subdirectories mod32, mod64, lib32 and lib64 for 32 bit and 64 bit
      if [ $PHAML_GRAPHICS = "mesa" ]
      then
         if [ $PHAML_OS = "irixn32" ]
         then
            OGLMODS="$MODFLAG"'$(MESA_HOME)/mod32'
            OGLLIBS='-L$(MESA_HOME)/lib32'" -lf90glut -lMesaf90GLU -lMesaf90GL -lglut -lMesaGLU -lMesaGL"
         else
            OGLMODS="$MODFLAG"'$(MESA_HOME)/mod64'
            OGLLIBS='-L$(MESA_HOME)/lib64'" -lf90glut -lMesaf90GLU -lMesaf90GL -lglut -lMesaGLU -lMesaGL"
         fi
      fi

# Likewise for f90gl for SGI's OpenGL

      if [ $PHAML_GRAPHICS = "opengl" ]
      then
         if [ $PHAML_OS = "irixn32" ]
         then
            OGLMODS="$MODFLAG"'$(OPENGL_HOME)/mod32'
            OGLLIBS='-L$(OPENGL_HOME)/lib32'" -lf90glut -lf90GLU -lf90GL -lglut -lGLU -lGL"
         else
            OGLMODS="$MODFLAG"'$(OPENGL_HOME)/mod64'
            OGLLIBS='-L$(OPENGL_HOME)/lib64'" -lf90glut -lf90GLU -lf90GL -lglut -lGLU -lGL"
         fi
      fi ;;

########
# suns is for our Sun workstations running solaris, which may have more than
# one Fortran 90 compiler

   suns)

# The Sun F90 needs some extra libraries
      if [ $PHAML_F90 = "sun" ]
      then
         OTHERLIBS="$OTHERLIBS -lnsl -lsocket"
      fi ;;

########
# tflop is the ASCI Red machine at Sandia

   tflop)
      FFLAGS="-cougar -O"
      CFLAGS="-cougar -O"
      MESSPASSINC=
      MESSPASSLIBS='-lmpi'
      ZOLTANMOD="$MODFLAG"'$(HOME)/Zoltan/Obj_tflop'
      ZOLTANLIBS='-L$(HOME)/Zoltan/Obj_tflop -lzoltan -lzoltan_comm -lzoltan_mem -L/Net/proj/zoltan/arch/tflop/lib -lparmetis -lmetis' ;;

esac

# end of Step 3.  You should not need to change anything below here.

case "$PHAML_ZOLTAN" in

   no)
      ZOLTANMOD=
      ZOLTANLIBS= ;;

   yes)
      ZOLTANLIBS=$ZOLTANLIB
      if [ $PHAML_PARMETIS = "yes" ]
      then
         ZOLTANLIBS="$ZOLTANLIBS $PARMETISLIB"
      fi
      if [ $PHAML_JOSTLE = "yes" ]
      then
         ZOLTANLIBS="$ZOLTANLIBS $JOSTLELIB"
      fi
      if [ $PHAML_PATOH = "yes" ]
      then
         ZOLTANLIBS="$ZOLTANLIBS $PATOHLIB"
      fi
      if [ $PHAML_PARKWAY = "yes" ]
      then
         ZOLTANLIBS="$ZOLTANLIBS $PARKWAYLIB"
      fi
      if [ $PHAML_NEMESIS = "yes" ]
      then
         ZOLTANLIBS="$ZOLTANLIBS $NEMESISLIB"
      fi
      if [ $PHAML_DRUM = "yes" ]
      then
         ZOLTANLIBS="$DRUMLIB $ZOLTANLIBS"
      fi ;;

esac

if [ $PHAML_PARALLEL = "messpass_nospawn" ]
then
   case "$PHAML_GRAPHICS" in
      metro|mesa|opengl)
         PHAML_GETS_GRAPHICSLIBS='yes' ;;
   esac
fi

##############################################################################

# You should not have to change anything below here.

##############################################################################

# export the shell variables needed by mkmkfile.sh in subdirectories

export PHAML_ARCH
export PHAML_OS
export PHAML_F90
export PHAML_C
export PHAML_HASHSIZE
export PHAML_PARALLEL
export PHAML_PARLIB
export PHAML_GRAPHICS
export PHAML_BLAS
export PHAML_LAPACK
export PHAML_ARPACK
export PHAML_BLOPEX
export PHAML_HYPRE
export PHAML_MUMPS
export PHAML_PETSC
export PHAML_SUPERLU
export PHAML_SYSTEM
export PHAML_ZOLTAN
export PHAML_PARMETIS
export PHAML_JOSTLE
export PHAML_PATOH
export PHAML_PARKWAY
export PHAML_NEMESIS
export PHAML_DRUM
export F90
export FFLAGS
export CC
export CFLAGS
export LINKER
export LINKFLAGS
export MPIF
export MODFLAG
export MODSUFFIX
export MAKELIB
export RANLIB
export PHAML_HOME
export MESSPASSLIBS
export OGLLIBS
export XLIBS
export LAPACKLIBS
export BLASLIBS
export ZOLTANLIBS
export ARPACKLIBS
export BLOPEXINC
export BLOPEXLIBS
export PETSCLIBS
export HYPRELIBS
export MUMPSLIBS
export SUPERLULIBS
export OTHERLIBS
export PHAML_GETS_GRAPHICSLIBS
export ZOLTANMOD
export ZOLTANINC
export OGLMODS
export PETSCINC
export MUMPSINC
export SUPERLUINC
export MESSPASSINC
export RUNMPI
export CPUSEC
export CPUSEC_COMP

##############################################################################

# execute mkmkfile.sh for the src directory, testdir, and all the examples

echo "make makefile for src"
cd src
./mkmkfile.sh

echo "make makefile for testdir"
cd ../testdir
./mkmkfile.sh

cd ../examples
for dir in `ls`
do
   echo "make makefile for $dir"
   cd $dir
   ./mkmkfile.sh
   cd ..
done

cd ..

##############################################################################

# create the top level Makefile

if [ -f Makefile ]
then
   rm -f Makefile
fi
f=Makefile

echo "#---------------------------------------------------------------------!" 1>>$f
echo "#                                PHAML                                !" 1>>$f
echo "#                                                                     !" 1>>$f
echo "# The Parallel Hierarchical Adaptive MultiLevel code for solving      !" 1>>$f
echo "# linear elliptic partial differential equations of the form          !" 1>>$f
echo "# (PUx)x + (QUy)y + RU = F on 2D polygonal domains with mixed         !" 1>>$f
echo "# boundary conditions, and eigenvalue problems where F is lambda*U.   !" 1>>$f
echo "#                                                                     !" 1>>$f
echo "# PHAML is public domain software.  It was produced as part of work   !" 1>>$f
echo "# done by the U.S. Government, and is not subject to copyright in     !" 1>>$f
echo "# the United States.                                                  !" 1>>$f
echo "#                                                                     !" 1>>$f
echo "#     William F. Mitchell                                             !" 1>>$f
echo "#     Mathematical and Computational Sciences Division                !" 1>>$f
echo "#     National Institute of Standards and Technology                  !" 1>>$f
echo "#     william.mitchell@nist.gov                                       !" 1>>$f
echo "#     http://math.nist.gov/phaml                                      !" 1>>$f
echo "#                                                                     !" 1>>$f
echo "#---------------------------------------------------------------------!" 1>>$f
echo "" 1>>$f
echo "# Makefile created for system configuration:" 1>>$f
echo "#   Architecture:     " $PHAML_ARCH 1>>$f
echo "#   OS:               " $PHAML_OS 1>>$f
echo "#   F90 compiler:     " $PHAML_F90 1>>$f
echo "#   C compiler:       " $PHAML_C 1>>$f
echo "#   Hash size:        " $PHAML_HASHSIZE 1>>$f
echo "#   Parallel form:    " $PHAML_PARALLEL 1>>$f
echo "#   Parallel library: " $PHAML_PARLIB 1>>$f
echo "#   Graphics:         " $PHAML_GRAPHICS 1>>$f
echo "#   BLAS:             " $PHAML_BLAS 1>>$f
echo "#   LAPACK:           " $PHAML_LAPACK 1>>$f
echo "#   ARPACK:           " $PHAML_ARPACK 1>>$f
echo "#   BLOPEX:           " $PHAML_BLOPEX 1>>$f
echo "#   hypre:            " $PHAML_HYPRE 1>>$f
echo "#   MUMPS:            " $PHAML_MUMPS 1>>$f
echo "#   PETSc:            " $PHAML_PETSC 1>>$f
echo "#   SUPERLU:          " $PHAML_SUPERLU 1>>$f
echo "#   Zoltan:           " $PHAML_ZOLTAN 1>>$f
echo "#   ParMETIS:         " $PHAML_PARMETIS 1>>$f
echo "#   JOSTLE:           " $PHAML_JOSTLE 1>>$f
echo "#   PaToH:            " $PHAML_PATOH 1>>$f
echo "#   ParKway:          " $PHAML_PARKWAY 1>>$f
echo "#   Nemesis:          " $PHAML_NEMESIS 1>>$f
echo "#   DRUM:             " $PHAML_DRUM 1>>$f
echo "#   Specific system:  " $PHAML_SYSTEM 1>>$f
echo "" 1>>$f

echo "phaml:" 1>>$f
echo "	cd src; make lib" 1>>$f
echo "" 1>>$f

echo "clean:" 1>>$f
echo "	cd src; make clean" 1>>$f
echo "	cd testdir; make clean" 1>>$f
for dir in `ls examples`
do
echo "	cd examples/$dir; make clean" 1>>$f
done
echo "" 1>>$f

echo "reallyclean:" 1>>$f
echo "	cd src; make clean" 1>>$f
echo "	cd testdir; make reallyclean" 1>>$f
for dir in `ls examples`
do
echo "	cd examples/$dir; make clean" 1>>$f
done
echo "	rm -f lib/*" 1>>$f
echo "	rm -f modules/*" 1>>$f
echo "	cd src; rm -f Makefile Makefile.bak" 1>>$f
echo "	cd testdir; rm -f Makefile Makefile.bak" 1>>$f
for dir in `ls examples`
do
echo "	cd examples/$dir; rm -f Makefile Makefile.bak" 1>>$f
done
echo "" 1>>$f

echo "test:" 1>>$f
echo '	cd testdir; make test what=$(what)' 1>>$f
echo "" 1>>$f
