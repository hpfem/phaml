#!/bin/bash

# Make sure you are using bash instead of a Bourne shell.
# The following are some places where I had to change the shell line above.

# Sun at Sandia, cross compiling for tflop
#!/Net/local/gnu/bin/bash

# SGI at NIST
#!/usr/freeware/bin/bash

##############################################################################

# Create Makefile.
# Most shell variables are set in mkmkfile.sh in the top directory.

# back up old Makefile, and set f to be Makefile for brevity

if [ -f Makefile ]
then
   mv -f Makefile Makefile.bak
fi
f=Makefile

# write the header comments into Makefile

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
echo "#   HYPRE:            " $PHAML_HYPRE 1>>$f
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

# write makefile variables

echo "F90=$F90" 1>>$f
echo "FFLAGS=$FFLAGS" 1>>$f
echo "LINKER=$LINKER" 1>>$f
echo "LINKFLAGS=$LINKFLAGS" 1>>$f
if [ $PHAML_PARLIB = "lam" ]
then
echo "LAMHF77=$MPIF" 1>>$f
fi
if [ $PHAML_PARLIB = "openmpi" ]
then
echo "OMPI_F77=$MPIF" 1>>$f
echo "OMPI_FC=$MPIF" 1>>$f
echo "export OMPI_F77" 1>>$f
echo "export OMPI_FC" 1>>$f
fi
echo "" 1>>$f
echo "PHAML_HOME=$PHAML_HOME" 1>>$f
echo 'PHAML_MODDIR=$(PHAML_HOME)/modules' 1>>$f
echo 'PHAML_LIBDIR=$(PHAML_HOME)/lib' 1>>$f
echo 'PHAML_SRCDIR=$(PHAML_HOME)/src' 1>>$f
echo "" 1>>$f

# choose whether this is master/slave or spmd

if [ $PHAML_PARALLEL = "messpass_nospawn" ]
then 
MAIN=spmd
else
MAIN=master
fi

# write the main targets

if [ $PHAML_PARALLEL = "messpass_nospawn" ]
then 
echo "all: phaml" 1>>$f
elif [ $PHAML_PARALLEL = "sequential" ]
then 
   if [ $PHAML_GRAPHICS = "none" ]
   then
      echo "all: phaml" 1>>$f
   else 
      echo "all: phaml phaml_graphics" 1>>$f
   fi
elif [ $PHAML_GRAPHICS = "none" ]
then
echo "all: phaml phaml_slave" 1>>$f
else
echo "all: phaml phaml_slave phaml_graphics" 1>>$f
fi
echo "" 1>>$f

if [ ! $PHAML_PARALLEL = "messpass_nospawn" -a ! $PHAML_GRAPHICS = "none" ]
then
echo "phaml_graphics: "'\' 1>>$f
echo "	"'graphmain.o \' 1>>$f
echo "	pde.o " 1>>$f
echo "	"'$(LINKER) $(LINKFLAGS) -o phaml_graphics \' 1>>$f
echo "	"'graphmain.o \' 1>>$f
echo "	pde.o "'\' 1>>$f
echo "	"'-L$(PHAML_LIBDIR) -lphaml \' 1>>$f
if [ -n "$SUPERLULIBS" ]
then
echo "	$SUPERLULIBS "'\' 1>>$f
fi
if [ -n "$HYPRELIBS" ]
then
echo "	$HYPRELIBS "'\' 1>>$f
fi
if [ -n "$ARPACKLIBS" ]
then
echo "	$ARPACKLIBS "'\' 1>>$f
fi
if [ -n "$ZOLTANLIBS" ]
then
echo "	$ZOLTANLIBS "'\' 1>>$f
fi
if [ -n "$PETSCLIBS" ]
then
echo "	$PETSCLIBS "'\' 1>>$f
fi
if [ -n "$BLOPEXLIBS" ]
then
echo "	$BLOPEXLIBS "'\' 1>>$f
fi
if [ -n "$MUMPSLIBS" ]
then
echo "	$MUMPSLIBS "'\' 1>>$f
fi
echo "	$MESSPASSLIBS "'\' 1>>$f
echo "	$OGLLIBS "'\' 1>>$f
if [ -n "$LAPACKLIBS" ]
then
echo "	$LAPACKLIBS "'\' 1>>$f
fi
if [ -n "$BLASLIBS" ]
then
echo "	$BLASLIBS "'\' 1>>$f
fi
if [ -n "$OTHERLIBS" ]
then
echo "	$OTHERLIBS "'\' 1>>$f
fi
echo "	$XLIBS" 1>>$f
echo "" 1>>$f
fi

if [ ! $PHAML_PARALLEL = "sequential" -a ! $PHAML_PARALLEL = "messpass_nospawn" ]
then
echo "phaml_slave: "'\' 1>>$f
echo "	"'slave.o \' 1>>$f
echo "	pde.o" 1>>$f
echo "	"'$(LINKER) $(LINKFLAGS) -o phaml_slave \' 1>>$f
echo "	"'slave.o \' 1>>$f
echo "	pde.o "'\' 1>>$f
echo "	"'-L$(PHAML_LIBDIR) -lphaml \' 1>>$f
if [ -n "$HYPRELIBS" ]
then
echo "	$HYPRELIBS "'\' 1>>$f
fi
if [ -n "$ZOLTANLIBS" ]
then
echo "	$ZOLTANLIBS "'\' 1>>$f
fi
if [ -n "$ARPACKLIBS" ]
then
echo "	$ARPACKLIBS "'\' 1>>$f
fi
if [ -n "$PETSCLIBS" ]
then
echo "	$PETSCLIBS "'\' 1>>$f
fi
if [ -n "$BLOPEXLIBS" ]
then
echo "	$BLOPEXLIBS "'\' 1>>$f
fi
if [ -n "$MUMPSLIBS" ]
then
echo "	$MUMPSLIBS "'\' 1>>$f
fi
if [ -n "$SUPERLULIBS" ]
then
echo "	$SUPERLULIBS "'\' 1>>$f
fi
if [ -n "$LAPACKLIBS" -o -n "$BLASLIBS" -o -n "$OTHERLIBS" ]
then
   echo "	$MESSPASSLIBS "'\' 1>>$f
else
   echo "	$MESSPASSLIBS" 1>>$f
fi
if [ -n "$LAPACKLIBS" ]
then
   if [ -n "$BLASLIBS" -o -n "$OTHERLIBS" ]
   then
      echo "	$LAPACKLIBS "'\' 1>>$f
   else
      echo "	$LAPACKLIBS " 1>>$f
   fi
fi
if [ -n "$BLASLIBS" ]
then
   if [ -n "$OTHERLIBS" ]
   then
      echo "	$BLASLIBS "'\' 1>>$f
   else
      echo "	$BLASLIBS " 1>>$f
   fi
fi
if [ -n "$OTHERLIBS" ]
then
echo "	$OTHERLIBS" 1>>$f
fi
echo "" 1>>$f
fi

echo "phaml: "'\' 1>>$f
echo "	$MAIN.o pde.o" 1>>$f
echo "	"'$(LINKER) $(LINKFLAGS) -o phaml \' 1>>$f
echo "	$MAIN.o pde.o "'\' 1>>$f
if [ -n "$SUPERLULIBS" -o -n "$HYPRELIBS" -o -n "$ZOLTANLIBS" -o -n "$ARPACKLIBS" -o -n "$PETSCLIBS" -o -n "$BLOPEXLIBS" -o -n "$MUMPSLIBS" -o -n "$LAPACKLIBS" -o -n "$BLASLIBS" -o -n "$MESSPASSLIBS" -o -n "$PHAML_GETS_GRAPHICSLIBS" -o -n "$OTHERLIBS" ]
then
echo "	"'-L$(PHAML_LIBDIR) -lphaml \' 1>>$f
else
echo "	"'-L$(PHAML_LIBDIR) -lphaml' 1>>$f
fi
if [ -n "$SUPERLULIBS" ]
then
   if [ -n "$HYPRELIBS" -o -n "$ZOLTANLIBS" -o -n "$ARPACKLIBS" -o -n "$PETSCLIBS" -o -n "$BLOPEXLIBS" -o -n "$MUMPSLIBS" -o -n "$LAPACKLIBS" -o -n "$BLASLIBS" -o -n "$MESSPASSLIBS" -o -n "$PHAML_GETS_GRAPHICSLIBS" -o -n "$OTHERLIBS" ]
   then
echo "	$SUPERLULIBS "'\' 1>>$f
   else
echo "	$SUPERLULIBS " 1>>$f
   fi
fi
if [ -n "$HYPRELIBS" ]
then
   if [ -n "$ZOLTANLIBS" -o -n "$ARPACKLIBS" -o -n "$PETSCLIBS" -o -n "$BLOPEXLIBS" -o -n "$MUMPSLIBS" -o -n "$LAPACKLIBS" -o -n "$BLASLIBS" -o -n "$MESSPASSLIBS" -o -n "$PHAML_GETS_GRAPHICSLIBS" -o -n "$OTHERLIBS" ]
   then
echo "	$HYPRELIBS "'\' 1>>$f
   else
echo "	$HYPRELIBS " 1>>$f
   fi
fi
if [ -n "$ZOLTANLIBS" ]
then
   if [ -n "$ARPACKLIBS" -o -n "$PETSCLIBS" -o -n "$BLOPEXLIBS" -o -n "$MUMPSLIBS" -o -n "$LAPACKLIBS" -o -n "$BLASLIBS" -o -n "$MESSPASSLIBS" -o -n "$PHAML_GETS_GRAPHICSLIBS" -o -n "$OTHERLIBS" ]
   then
echo "	$ZOLTANLIBS "'\' 1>>$f
   else
echo "	$ZOLTANLIBS " 1>>$f
   fi
fi
if [ -n "$ARPACKLIBS" ]
then
   if [ -n "$PETSCLIBS" -o -n "$BLOPEXLIBS" -o -n "$MUMPSLIBS" -o -n "$LAPACKLIBS" -o -n "$BLASLIBS" -o -n "$MESSPASSLIBS" -o -n "$PHAML_GETS_GRAPHICSLIBS" -o -n "$OTHERLIBS" ]
   then
echo "	$ARPACKLIBS "'\' 1>>$f
   else
echo "	$ARPACKLIBS " 1>>$f
   fi
fi
if [ -n "$PETSCLIBS" ]
then
   if [ -n "$BLOPEXLIBS" -o -n "$MUMPSLIBS" -o -n "$LAPACKLIBS" -o -n "$BLASLIBS" -o -n "$MESSPASSLIBS" -o -n "$PHAML_GETS_GRAPHICSLIBS" -o -n "$OTHERLIBS" ]
   then
echo "	$PETSCLIBS "'\' 1>>$f
   else
echo "	$PETSCLIBS " 1>>$f
   fi
fi
if [ -n "$BLOPEXLIBS" ]
then
   if [ -n "$MUMPSLIBS" -o -n "$LAPACKLIBS" -o -n "$BLASLIBS" -o -n "$MESSPASSLIBS" -o -n "$PHAML_GETS_GRAPHICSLIBS" -o -n "$OTHERLIBS" ]
   then
echo "	$BLOPEXLIBS "'\' 1>>$f
   else
echo "	$BLOPEXLIBS " 1>>$f
   fi
fi
if [ -n "$MUMPSLIBS" ]
then
   if [ -n "$LAPACKLIBS" -o -n "$BLASLIBS" -o -n "$MESSPASSLIBS" -o -n "$PHAML_GETS_GRAPHICSLIBS" -o -n "$OTHERLIBS" ]
   then
echo "	$MUMPSLIBS "'\' 1>>$f
   else
echo "	$MUMPSLIBS " 1>>$f
   fi
fi
if [ -n "$MESSPASSLIBS" ]
then
   if [ -n "$LAPACKLIBS" -o -n "$BLASLIBS" -o -n "$PHAML_GETS_GRAPHICSLIBS" -o -n "$OTHERLIBS" ]
   then
echo "	$MESSPASSLIBS "'\' 1>>$f
   else
echo "	$MESSPASSLIBS " 1>>$f
   fi
fi
# OGLLIBS must come before anything that gives a directory containing glut or gl
if [ -n "$PHAML_GETS_GRAPHICSLIBS" ]
then
echo "  $OGLLIBS "'\' 1>>$f
fi
if [ -n "$LAPACKLIBS" ]
then
   if [ -n "$BLASLIBS" -o -n "$PHAML_GETS_GRAPHICSLIBS" -o -n "$OTHERLIBS" ]
   then
echo "	$LAPACKLIBS "'\' 1>>$f
   else
echo "	$LAPACKLIBS " 1>>$f
   fi
fi
if [ -n "$BLASLIBS" ]
then
   if [ -n "$PHAML_GETS_GRAPHICSLIBS" -o -n "$OTHERLIBS" ]
   then
echo "	$BLASLIBS "'\' 1>>$f
   else
echo "	$BLASLIBS " 1>>$f
   fi
fi
if [ -n "$OTHERLIBS" ]
then
   if [ -n "$PHAML_GETS_GRAPHICSLIBS" ]
   then
echo "	$OTHERLIBS "'\' 1>>$f
   else
echo "	$OTHERLIBS " 1>>$f
   fi
fi
if [ -n "$PHAML_GETS_GRAPHICSLIBS" ]
then
echo "	$XLIBS" 1>>$f
fi
echo "" 1>>$f

# write the rules for compiling files

if [ $PHAML_PARALLEL = "messpass_spawn" ]
then
echo 'slave.o: $(PHAML_SRCDIR)/slave.f90' 1>>$f
echo "	"'$(F90) $(FFLAGS) '"$MODFLAG"'$(PHAML_MODDIR)'" $ZOLTANMOD "'-o slave.o -c $(PHAML_SRCDIR)/slave.f90' 1>>$f
echo "" 1>>$f
fi

case "$PHAML_GRAPHICS" in
   metro|mesa|opengl)
      if [ ! $PHAML_PARALLEL = "messpass_nospawn" ]
      then
echo 'graphmain.o: $(PHAML_SRCDIR)/graphmain.f90' 1>>$f
echo "	"'$(F90) $(FFLAGS) '"$MODFLAG"'$(PHAML_MODDIR) -o graphmain.o -c $(PHAML_SRCDIR)/graphmain.f90' 1>>$f
echo "" 1>>$f
      fi ;;
esac

echo "$MAIN.o: $MAIN.f90" 1>>$f
echo "	"'$(F90) $(FFLAGS) '"$MODFLAG"'$(PHAML_MODDIR)'" $ZOLTANMOD -c $MAIN.f90" 1>>$f
echo "" 1>>$f

echo "pde.o: pde.f90" 1>>$f
echo "	"'$(F90) $(FFLAGS) '"$MODFLAG"'$(PHAML_MODDIR)'" $ZOLTANMOD -c pde.f90" 1>>$f
echo "" 1>>$f

echo "clean:" 1>>$f
echo "	rm -f *.o *.mod phaml phaml_slave phaml_graphics *.M *.stb" 1>>$f
