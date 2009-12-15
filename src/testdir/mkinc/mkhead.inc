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

# write the main targets

echo 'all: \' 1>>$f
