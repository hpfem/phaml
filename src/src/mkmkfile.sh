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
echo "#   hypre:            " $PHAML_HYPRE 1>>$f
echo "#   MUMPS:            " $PHAML_MUMPS 1>>$f
echo "#   PETSc:            " $PHAML_PETSC 1>>$f
echo "#   SuperLU:          " $PHAML_SUPERLU 1>>$f
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
echo "CC=$CC" 1>>$f
echo "CFLAGS=$CFLAGS" 1>>$f
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
echo "MAKELIB=$MAKELIB" 1>>$f
echo "RANLIB=$RANLIB" 1>>$f
echo "" 1>>$f

PHAML_MODDIR=$PHAML_HOME/modules
PHAML_LIBDIR=$PHAML_HOME/lib

# write the library target

echo "lib: "'\' 1>>$f
echo "	global.o cpusec.o phaml.o linsys.o linsys_util.o linsystype.o "'\' 1>>$f
echo "	eigen.o sort.o hbmg.o lapack_solve.o make_linsys.o linsys_io.o "'\' 1>>$f
echo "	quadrules.o basis.o errest.o evaluate.o sysdep.o "'\' 1>>$f
case "$PHAML_PARLIB" in
   lam|mpi|mpich|mpich2|myrinet|openmpi)
echo "	mpi_stringf.o mpi_stringc.o mpipack1.o mpipack2.o mpif.o "'\' 1>>$f
      if [ $PHAML_PARALLEL = "messpass_spawn" ]
      then
echo "	mpi2_stringf.o "'\' 1>>$f
      fi ;;
   pvm)
echo "	pvmpack1.o pvmpack2.o pvmf.o "'\' 1>>$f ;;
   none)
      : ;;
   *)
echo "*** mkmkfile.sh is missing case $PHAML_PARLIB for PHAML_PARLIB under phaml:" ;;
esac
case "$PHAML_ARPACK" in
   yes)
      case "$PHAML_PARLIB" in
         lam|mpi|mpich|mpich2|myrinet|openmpi)
            : ;;
         pvm|none)
echo "	parpack.o "'\' 1>>$f ;;
      esac ;;
   no)
echo "	arpack.o parpack.o "'\' 1>>$f ;;
   *)
echo "*** mkmkfile.sh is missing case $PHAML_ARPACK for PHAML_ARPACK under phaml:" ;;
esac
case "$PHAML_GRAPHICS" in
   metro|mesa|opengl)
echo "	graphics.o modview.o rendereps.o "'\' 1>>$f ;;
   none)
      if [ $PHAML_PARALLEL = "messpass_nospawn" ]
      then
echo "	graphics.o "'\' 1>>$f
      fi ;;
   *)
echo "*** mkmkfile.sh is missing case $PHAML_GRAPHICS for PHAML_GRAPHICS under phaml:" ;;
esac
if [ -z "$LAPACKLIBS" ]
then
echo "	lapack.o "'\' 1>>$f
fi
if [ -z "$BLASLIBS" ]
then
echo "	blas.o "'\' 1>>$f
fi
echo "	stopwatch.o hash.o hashext.o phamltype.o messpass.o grid.o "'\' 1>>$f
echo "	grid_io.o loadbal.o gridtype.o zoltan_interf.o petsc_interf.o "'\' 1>>$f
echo "	petsc_init.o petsctype.o calgo582.o dbinom.o mumps_interf.o mumps_struc.o "'\' 1>>$f
echo "	blopex_c.o krylov.o templates.o superlu_interf.o superlutype.o "'\' 1>>$f
if [ $PHAML_SUPERLU = "yes" ]
then
echo "	superlu_wrappers.o superlu_fwrap.o "'\' 1>>$f
fi
if [ $PHAML_ZOLTAN = "yes" ]
then
echo "	zoltanParams_read_file.o "'\' 1>>$f
fi
if [ $PHAML_HYPRE = "yes" ]
then
echo "	hypre_fix.o "'\' 1>>$f
fi
echo "	hypretype.o hypre_interf.o" 1>>$f
echo "	"'$(MAKELIB) libphaml.a \' 1>>$f
echo "	global.o cpusec.o phaml.o linsys.o linsystype.o linsys_util.o "'\' 1>>$f
echo "	eigen.o sort.o hbmg.o lapack_solve.o make_linsys.o linsys_io.o "'\' 1>>$f
echo "	quadrules.o basis.o errest.o evaluate.o sysdep.o "'\' 1>>$f
case "$PHAML_PARLIB" in
   lam|mpi|mpich|mpich2|myrinet|openmpi)
echo "	mpi_stringf.o mpi_stringc.o mpipack1.o mpipack2.o mpif.o "'\' 1>>$f
      if [ $PHAML_PARALLEL = "messpass_spawn" ]
      then
echo "	mpi2_stringf.o "'\' 1>>$f
      fi ;;
   pvm)
echo "	pvmpack1.o pvmpack2.o pvmf.o "'\' 1>>$f ;;
   none)
      : ;;
   *)
echo "*** mkmkfile.sh is missing case $PHAML_PARLIB for PHAML_PARLIB under phaml:" ;;
esac
case "$PHAML_ARPACK" in
   yes)
      case "$PHAML_PARLIB" in
         lam|mpi|mpich|mpich2|myrinet|openmpi)
            : ;;
         pvm|none)
echo "	parpack.o "'\' 1>>$f ;;
      esac ;;
   no)
echo "	arpack.o parpack.o "'\' 1>>$f ;;
   *)
echo "*** mkmkfile.sh is missing case $PHAML_ARPACK for PHAML_ARPACK under phaml:" ;;
esac
case "$PHAML_GRAPHICS" in
   metro|mesa|opengl)
echo "	graphics.o modview.o rendereps.o "'\' 1>>$f ;;
   none)
      if [ $PHAML_PARALLEL = "messpass_nospawn" ]
      then
echo "	graphics.o "'\' 1>>$f
      fi ;;
   *)
echo "*** mkmkfile.sh is missing case $PHAML_GRAPHICS for PHAML_GRAPHICS under phaml:" ;;
esac
if [ -z "$LAPACKLIBS" ]
then
echo "	lapack.o "'\' 1>>$f
fi
if [ -z "$BLASLIBS" ]
then
echo "	blas.o "'\' 1>>$f
fi
echo "	stopwatch.o hash.o hashext.o phamltype.o messpass.o grid.o "'\' 1>>$f
echo "	grid_io.o loadbal.o gridtype.o zoltan_interf.o petsc_interf.o "'\' 1>>$f
echo "	petsc_init.o petsctype.o calgo582.o dbinom.o mumps_interf.o mumps_struc.o "'\' 1>>$f
echo "	blopex_c.o krylov.o templates.o superlu_interf.o superlutype.o "'\' 1>>$f
if [ $PHAML_SUPERLU = "yes" ]
then
echo "	superlu_wrappers.o superlu_fwrap.o "'\' 1>>$f
fi
if [ $PHAML_ZOLTAN = "yes" ]
then
echo "	zoltanParams_read_file.o "'\' 1>>$f
fi
if [ $PHAML_HYPRE = "yes" ]
then
echo "	hypre_fix.o "'\' 1>>$f
fi
echo "	hypretype.o hypre_interf.o" 1>>$f
echo "	"'$(RANLIB) libphaml.a' 1>>$f
echo "	cp -f libphaml.a $PHAML_LIBDIR" 1>>$f
echo "	head -44 Makefile > $PHAML_LIBDIR/CONFIG" 1>>$f
echo "	cp -f *.$MODSUFFIX $PHAML_MODDIR" 1>>$f
echo "	head -44 Makefile > $PHAML_MODDIR/CONFIG" 1>>$f
echo "" 1>>$f

# write the rules for compiling files

echo "phaml.o: phaml.f90 "'\' 1>>$f
echo "        global.o linsys.o stopwatch.o messpass.o grid.o hash.o "'\' 1>>$f
echo "        phamltype.o gridtype.o grid_io.o loadbal.o linsys_io.o "'\' 1>>$f
echo "        zoltan_interf.o errest.o evaluate.o sysdep.o" 1>>$f
if [ $PHAML_PARLIB = "openmpi" ]
then
echo "	@echo" 1>>$f
echo "	@echo mpif90 is:" 1>>$f
echo "	@mpif90 --showme" 1>>$f
echo "	@echo" 1>>$f
fi
if [ $PHAML_PARLIB = "lam" ]
then
echo "	@echo" 1>>$f
echo "	@echo mpif77 is:" 1>>$f
echo "	@mpif77 --showme" 1>>$f
echo "	@echo" 1>>$f
fi
echo "	"'$(F90) $(FFLAGS)'" $ZOLTANMOD -c phaml.f90" 1>>$f
echo "" 1>>$f
echo "phamltype.o: phamltype.f90 "'\' 1>>$f
echo "        global.o messpass.o gridtype.o zoltan_interf.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" $ZOLTANMOD -c phamltype.f90" 1>>$f
echo "" 1>>$f
echo "grid.o: grid.f90 "'\' 1>>$f
echo "        global.o stopwatch.o messpass.o hash.o gridtype.o sort.o "'\' 1>>$f
echo "        make_linsys.o errest.o sysdep.o linsystype.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -c grid.f90" 1>>$f
echo "" 1>>$f
echo "gridtype.o: gridtype.f90 "'\' 1>>$f
echo "        global.o hash.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -c gridtype.f90" 1>>$f
echo "" 1>>$f
echo "grid_io.o: grid_io.f90 "'\' 1>>$f
echo "        global.o stopwatch.o messpass.o hash.o gridtype.o "'\' 1>>$f
echo "        zoltan_interf.o errest.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" $ZOLTANMOD -c grid_io.f90" 1>>$f
echo "" 1>>$f
echo "errest.o: errest.f90 "'\' 1>>$f
echo "        global.o messpass.o hash.o gridtype.o hash.o make_linsys.o "'\' 1>>$f
echo "        evaluate.o quadrules.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -c errest.f90" 1>>$f
echo "" 1>>$f
echo "evaluate.o: evaluate.f90 "'\' 1>>$f
echo "        global.o messpass.o gridtype.o basis.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -c evaluate.f90" 1>>$f
echo "" 1>>$f
echo "linsys.o: linsys.f90 "'\' 1>>$f
echo "        stopwatch.o global.o messpass.o linsystype.o linsys_util.o "'\' 1>>$f
echo "        hash.o gridtype.o eigen.o hbmg.o lapack_solve.o linsys_io.o "'\' 1>>$f
echo "        petsc_interf.o mumps_interf.o make_linsys.o krylov.o "'\' 1>>$f
echo "        hypre_interf.o superlu_interf.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -c linsys.f90" 1>>$f
echo "" 1>>$f
echo "linsystype.o: linsystype.f90 "'\' 1>>$f
echo "        global.o hash.o hashext.o messpass.o gridtype.o petsctype.o "'\' 1>>$f
echo "        mumps_struc.o hypretype.o superlutype.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -c linsystype.f90" 1>>$f
echo "" 1>>$f
echo "linsys_util.o: linsys_util.f90 "'\' 1>>$f
echo "        global.o messpass.o stopwatch.o linsystype.o hash.o hashext.o "'\' 1>>$f
echo "        sort.o gridtype.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -c linsys_util.f90" 1>>$f
echo "" 1>>$f
echo "krylov.o: krylov.f90 "'\' 1>>$f
echo "        global.o messpass.o stopwatch.o linsystype.o linsys_util.o "'\' 1>>$f
echo "        gridtype.o hbmg.o errest.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -c krylov.f90" 1>>$f
echo "" 1>>$f
echo "templates.o: templates.f" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -c templates.f" 1>>$f
echo "" 1>>$f
echo "eigen.o: eigen.f90 "'\' 1>>$f
echo "        global.o messpass.o linsystype.o linsys_util.o sort.o "'\' 1>>$f
echo "        gridtype.o hbmg.o lapack_solve.o petsc_interf.o krylov.o "'\' 1>>$f
echo "        mumps_interf.o mumps_struc.o hypre_interf.o superlu_interf.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -c eigen.f90" 1>>$f
echo "" 1>>$f
echo "sort.o: sort.f" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -c sort.f" 1>>$f
echo "" 1>>$f
echo "hbmg.o: hbmg.f90 "'\' 1>>$f
echo "        global.o messpass.o linsystype.o linsys_util.o "'\' 1>>$f
echo "        gridtype.o lapack_solve.o linsys_io.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -c hbmg.f90" 1>>$f
echo "" 1>>$f
echo "make_linsys.o: make_linsys.f90 "'\' 1>>$f
echo "        stopwatch.o global.o messpass.o linsystype.o evaluate.o "'\' 1>>$f
echo "        gridtype.o hash.o hashext.o linsys_util.o quadrules.o basis.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -c make_linsys.f90" 1>>$f
echo "" 1>>$f
echo "quadrules.o: quadrules.f90 "'\' 1>>$f
echo "        global.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -c quadrules.f90" 1>>$f
echo "" 1>>$f
echo "basis.o: basis.f90 "'\' 1>>$f
echo "        global.o quadrules.o messpass.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -c basis.f90" 1>>$f
echo "" 1>>$f
echo "linsys_io.o: linsys_io.f90 "'\' 1>>$f
echo "        stopwatch.o global.o messpass.o linsystype.o sysdep.o "'\' 1>>$f
echo "        gridtype.o make_linsys.o errest.o quadrules.o evaluate.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -c linsys_io.f90" 1>>$f
echo "" 1>>$f
echo "lapack_solve.o: lapack_solve.f90 "'\' 1>>$f
echo "        global.o gridtype.o messpass.o linsystype.o linsys_util.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -c lapack_solve.f90" 1>>$f
echo "" 1>>$f
echo "loadbal.o: loadbal.f90 "'\' 1>>$f
echo "        global.o stopwatch.o messpass.o hash.o "'\' 1>>$f
echo "        gridtype.o grid.o zoltan_interf.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" $ZOLTANMOD -c loadbal.f90" 1>>$f
echo "" 1>>$f
case "$PHAML_ZOLTAN" in
   yes)
echo "zoltan_interf.o: zoltan_interf.f90 "'\' 1>>$f
echo "        global.o messpass.o hash.o gridtype.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" $ZOLTANMOD -c zoltan_interf.f90" 1>>$f
echo "" 1>>$f
case "$PHAML_DRUM" in
   yes)
echo "zoltanParams_read_file.o: zoltanParams_read_file.c zoltanParams.h" 1>>$f
echo "	"'$(CC) $(CFLAGS)'" $ZOLTANINC $MESSPASSINC -DFORTRAN_INTERFACE -c zoltanParams_read_file.c" 1>>$f ;;
   no)
echo "zoltanParams_read_file.o: zoltanP_r_f_dum.f90" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -o zoltanParams_read_file.o -c zoltanP_r_f_dum.f90" 1>>$f;;
esac
echo "" 1>>$f ;;
   no)
echo "zoltan_interf.o: zoltan_interf.dum.f90 "'\' 1>>$f
echo "        global.o messpass.o hash.o gridtype.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -o zoltan_interf.o -c zoltan_interf.dum.f90" 1>>$f
echo "" 1>>$f ;;
esac

case "$PHAML_MUMPS" in
   yes)
echo "mumps_struc.o: mumps_struc.f" 1>>$f
      if [ $PHAML_F90 = "xlf" ]
      then
echo "	ln -s mumps_struc.f mumps_struc.f90" 1>>$f
echo "	"'$(F90) $(FFLAGS) -qfixed'" $MUMPSINC -c mumps_struc.f90" 1>>$f
echo "	rm -f mumps_struc.f90" 1>>$f
      else
echo "	"'$(F90) $(FFLAGS)'" $MUMPSINC -c mumps_struc.f" 1>>$f
      fi
echo "" 1>>$f
echo "mumps_interf.o: mumps_interf.f90 "'\' 1>>$f
echo "        global.o messpass.o linsystype.o linsys_util.o mumps_struc.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -c mumps_interf.f90" 1>>$f
echo "" 1>>$f ;;
   no)
echo "mumps_struc.o: mumps_struc_dum.f90" 1>>$f
echo "	"'$(F90) $(FFLAGS) -o mumps_struc.o -c mumps_struc_dum.f90' 1>>$f
echo "" 1>>$f
echo "mumps_interf.o: mumps_interf.dum.f90 "'\' 1>>$f
echo "        global.o messpass.o linsystype.o linsys_util.o mumps_struc.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -o mumps_interf.o -c mumps_interf.dum.f90" 1>>$f
echo "" 1>>$f ;;
esac

case "$PHAML_SUPERLU" in
   yes)
echo "superlu_interf.o: superlu_interf.f90 "'\' 1>>$f
echo "	global.o messpass.o linsystype.o superlutype.o "'\' 1>>$f
echo "	hash.o hashext.o linsys_util.o lapack_solve.o superlu_fwrap.o " 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -c superlu_interf.f90" 1>>$f
echo "" 1>>$f
echo "superlu_fwrap.o: superlu_fwrap.f90 "'\' 1>>$f
echo "	superlutype.o global.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -c superlu_fwrap.f90" 1>>$f
echo "" 1>>$f
echo "superlutype.o: superlutype.f90 "'\' 1>>$f
echo "	global.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -c superlutype.f90" 1>>$f
echo "" 1>>$f
echo "superlu_wrappers.o: superlu_wrappers.c" 1>>$f
echo "	"'$(CC) $(CFLAGS)'" $SUPERLUINC $MESSPASSINC -c superlu_wrappers.c" 1>>$f
echo "" 1>>$f ;;
   no)
echo "superlu_interf.o: superlu_interf.dum.f90 "'\' 1>>$f
echo "	global.o messpass.o linsystype.o superlutype.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -o superlu_interf.o -c superlu_interf.dum.f90" 1>>$f
echo "" 1>>$f
echo "superlutype.o: superlutype.f90 "'\' 1>>$f
echo "	global.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -c superlutype.f90" 1>>$f
echo "" 1>>$f ;;
esac

case "$PHAML_BLOPEX" in
   yes)
echo "blopex_c.o: blopex_c.c " 1>>$f
echo "	"'$(CC) $(CFLAGS)'" $BLOPEXINC -c blopex_c.c" 1>>$f
echo "" 1>>$f ;;
   withhypre)
echo "blopex_c.o: blopex_c.hypre.c " 1>>$f
echo "	"'$(CC) $(CFLAGS)'" $BLOPEXINC -o blopex_c.o -c blopex_c.hypre.c" 1>>$f
echo "" 1>>$f ;;
   withpetsc)
echo "blopex_c.o: blopex_c.petsc.c " 1>>$f
echo "	"'$(CC) $(CFLAGS)'" $BLOPEXINC $PETSCINC $MESSPASSINC -o blopex_c.o -c blopex_c.petsc.c" 1>>$f
echo "" 1>>$f ;;
   no)
echo "blopex_c.o: blopex_c.dum.f90 "'\' 1>>$f
echo "        messpass.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -o blopex_c.o -c blopex_c.dum.f90" 1>>$f
echo "" 1>>$f ;;
esac

case "$PHAML_PETSC" in
   yes)
echo "petsc_interf.o: petsc_interf.F90 "'\' 1>>$f
echo "        global.o hash.o hashext.o messpass.o petsctype.o gridtype.o hbmg.o lapack_solve.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" $PETSCINC -c petsc_interf.F90" 1>>$f
echo "" 1>>$f
echo "petsc_init.o: petsc_init.F90 "'\' 1>>$f
echo "        global.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" $PETSCINC -c petsc_init.F90" 1>>$f
echo "" 1>>$f
echo "petsctype.o: petsctype.F90 "'\' 1>>$f
echo "        global.o hash.o hashext.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" $PETSCINC -c petsctype.F90" 1>>$f
echo "" 1>>$f ;;
   no)
echo "petsc_interf.o: petsc_interf.dum.f90 "'\' 1>>$f
echo "        global.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -o petsc_interf.o -c petsc_interf.dum.f90" 1>>$f
echo "" 1>>$f
echo "petsc_init.o: petsc_init.dum.f90 "'\' 1>>$f
echo "        global.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -o petsc_init.o -c petsc_init.dum.f90" 1>>$f
echo "" 1>>$f
echo "petsctype.o: petsctype.dum.f90 "'\' 1>>$f
echo "        global.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -o petsctype.o -c petsctype.dum.f90" 1>>$f
echo "" 1>>$f ;;
esac

echo "hypretype.o: hypretype.f90" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -c hypretype.f90" 1>>$f
echo "" 1>>$f
if [ $PHAML_HYPRE = "yes" ]
then
   case "$PHAML_PARLIB" in
      lam|mpich2|openmpi)
echo "hypre_fix.o: hypre_fix.c" 1>>$f
echo "	"'$(CC) $(CFLAGS)'" $MESSPASSINC -DHYPRE_FIX=HFMPI2 -o hypre_fix.o -c hypre_fix.c" 1>>$f
echo "" 1>>$f ;;
      *)
echo "hypre_fix.o: hypre_fix.c" 1>>$f
echo "	"'$(CC) $(CFLAGS)'" $MESSPASSINC -DHYPRE_FIX=HFASIS -o hypre_fix.o -c hypre_fix.c" 1>>$f
echo "" 1>>$f ;;
   esac
echo "hypre_interf.o: hypre_interf.f90 "'\' 1>>$f
echo "        global.o messpass.o linsystype.o hash.o hashext.o "'\' 1>>$f
echo "        hypretype.o linsys_util.o lapack_solve.o gridtype.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -o hypre_interf.o -c hypre_interf.f90" 1>>$f
echo "" 1>>$f
else
echo "hypre_interf.o: hypre_interf.dum.f90 "'\' 1>>$f
echo "        global.o messpass.o linsystype.o hypretype.o gridtype.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -o hypre_interf.o -c hypre_interf.dum.f90" 1>>$f
echo "" 1>>$f
fi
echo "hash.o: hash$PHAML_HASHSIZE.f90 "'\' 1>>$f
echo "        global.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -o hash.o -c hash$PHAML_HASHSIZE.f90" 1>>$f
echo "" 1>>$f
echo "hashext.o: hashext$PHAML_HASHSIZE.f90 "'\' 1>>$f
echo "        global.o hash.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -o hashext.o -c hashext$PHAML_HASHSIZE.f90" 1>>$f
echo "" 1>>$f
echo "global.o: global.f90 "'\' 1>>$f
echo "        stopwatch.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -c global.f90" 1>>$f
echo "" 1>>$f
echo "sysdep.o: sysdep.f90" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -c sysdep.f90" 1>>$f
echo "" 1>>$f
echo "stopwatch.o: stopwatch.f90" 1>>$f
if [ $PHAML_PARLIB = "openmpi" ]
then
echo "	@echo" 1>>$f
echo "	@echo mpif90 is:" 1>>$f
echo "	@mpif90 --showme" 1>>$f
echo "	@echo" 1>>$f
fi
if [ $PHAML_PARLIB = "lam" ]
then
echo "	@echo" 1>>$f
echo "	@echo mpif77 is:" 1>>$f
echo "	@mpif77 --showme" 1>>$f
echo "	@echo" 1>>$f
fi
echo "	"'$(F90) $(FFLAGS)'" -c stopwatch.f90" 1>>$f
echo "" 1>>$f
echo "cpusec.o: $CPUSEC" 1>>$f
echo "	$CPUSEC_COMP -o cpusec.o -c $CPUSEC" 1>>$f
echo "" 1>>$f

case "$PHAML_ARPACK" in
   yes)
      case "$PHAML_PARLIB" in
         lam|mpi|mpich|mpich2|myrinet|openmpi)
            : ;;
         pvm|none)
echo "parpack.o: parpack.dum.f "'\' 1>>$f
echo "          global.o" 1>>$f
            if [ $PHAML_F90 = "xlf" ]
            then
echo "	ln -s parpack.dum.f parpack.dum.f90" 1>>$f
echo "	"'$(F90) $(FFLAGS) -qfixed'" -o parpack.o -c parpack.dum.f90" 1>>$f
echo "	rm -f parpack.dum.f90" 1>>$f
            else
echo "	"'$(F90) $(FFLAGS)'" -o parpack.o -c parpack.dum.f" 1>>$f
            fi
echo "" 1>>$f ;;
      esac ;;
   no)
echo "arpack.o: arpack.dum.f "'\' 1>>$f
echo "          global.o" 1>>$f
      if [ $PHAML_F90 = "xlf" ]
      then
echo "	ln -s arpack.dum.f arpack.dum.f90" 1>>$f
echo "	"'$(F90) $(FFLAGS) -qfixed'" -o arpack.o -c arpack.dum.f90" 1>>$f
echo "	rm -f arpack.dum.f90" 1>>$f
      else
echo "	"'$(F90) $(FFLAGS)'" -o arpack.o -c arpack.dum.f" 1>>$f
      fi
echo "" 1>>$f
echo "parpack.o: parpack.dum.f "'\' 1>>$f
echo "          global.o" 1>>$f
      if [ $PHAML_F90 = "xlf" ]
      then
echo "	ln -s parpack.dum.f parpack.dum.f90" 1>>$f
echo "	"'$(F90) $(FFLAGS) -qfixed'" -o parpack.o -c parpack.dum.f90" 1>>$f
echo "	rm -f parpack.dum.f90" 1>>$f
      else
echo "	"'$(F90) $(FFLAGS)'" -o parpack.o -c parpack.dum.f" 1>>$f
      fi
echo "" 1>>$f ;;
   *)
echo "*** mkmkfile.sh is missing case $PHAML_ARPACK for PHAML_ARPACK under arpack.o" ;;
esac

case "$PHAML_PARLIB" in
   lam|mpi|mpich|mpich2|myrinet|openmpi)
echo "mpi_stringf.o: mpi_stringf.f90" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -c mpi_stringf.f90" 1>>$f
echo "" 1>>$f
echo "mpi_stringc.o: mpi_stringc.c" 1>>$f
echo "	"'$(CC) $(CFLAGS)'" $MESSPASSINC -c mpi_stringc.c" 1>>$f
echo "" 1>>$f
echo "mpipack1.o: mpipack1.f90" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -c mpipack1.f90" 1>>$f
echo "" 1>>$f
echo "mpipack2.o: mpipack2.f90 "'\' 1>>$f
echo "        global.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -c mpipack2.f90" 1>>$f
echo "" 1>>$f
echo "mpif.o: mpif.f" 1>>$f
if [ $PHAML_F90 = "xlf" ]
then
echo "	ln -s mpif.f mpif.f90" 1>>$f
echo "	"'$(F90) $(FFLAGS) -qfixed'" $MESSPASSINC -c mpif.f90" 1>>$f
echo "	rm -f mpif.f90" 1>>$f
else
echo "	"'$(F90) $(FFLAGS)'" $MESSPASSINC -c mpif.f" 1>>$f
fi
echo "" 1>>$f
   case "$PHAML_PARALLEL" in
      messpass_nospawn)
echo "messpass.o: messpass.mpi1.f90 mpif.o "'\' 1>>$f
echo "        global.o hash.o hashext.o stopwatch.o petsc_init.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -o messpass.o -c messpass.mpi1.f90" 1>>$f
echo "" 1>>$f ;;
      messpass_spawn)
echo "messpass.o: messpass.mpi2.f90 mpif.o "'\' 1>>$f
echo "        global.o hash.o hashext.o stopwatch.o petsc_init.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -o messpass.o -c messpass.mpi2.f90" 1>>$f
echo "" 1>>$f
echo "mpi2_stringf.o: mpi2_stringf.f90" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -c mpi2_stringf.f90" 1>>$f
echo "" 1>>$f ;;
   esac ;;
   pvm)
echo "pvmpack1.o: pvmpack1.f90" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -c pvmpack1.f90" 1>>$f
echo "" 1>>$f
echo "pvmpack2.o: pvmpack2.f90 "'\' 1>>$f
echo "        global.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -c pvmpack2.f90" 1>>$f
echo "" 1>>$f
echo "pvmf.o: pvmf.f" 1>>$f
if [ $PHAML_F90 = "xlf" ]
then
echo "	ln -s pvmf.f pvmf.f90" 1>>$f
echo "	"'$(F90) $(FFLAGS) -qfixed'" $MESSPASSINC -c pvmf.f90" 1>>$f
echo "	rm -f pvmf.f90" 1>>$f
else
echo "	"'$(F90) $(FFLAGS)'" $MESSPASSINC -c pvmf.f" 1>>$f
fi
echo "" 1>>$f
echo "messpass.o: messpass.pvm.f90 pvmf.o "'\' 1>>$f
echo "        global.o hash.o hashext.o stopwatch.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -o messpass.o -c messpass.pvm.f90" 1>>$f
echo "" 1>>$f ;;
   none)
echo "messpass.o: messpass.dummy.f90 "'\' 1>>$f
echo "        global.o hash.o hashext.o stopwatch.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -o messpass.o -c messpass.dummy.f90" 1>>$f
echo "" 1>>$f ;;
   *)
echo "*** mkmkfile.sh is missing case $PHAML_PARLIB for PHAML_PARLIB under messpass.o" ;;
esac
echo "calgo582.o: calgo582.f" 1>>$f
if [ $PHAML_F90 = "xlf" ]
then
echo "	ln -s calgo582.f calgo582.f90" 1>>$f
echo "	"'$(F90) $(FFLAGS) -qfixed'" -c calgo582.f90" 1>>$f
echo "	rm -f calgo582.f90" 1>>$f
else
echo "	"'$(F90) $(FFLAGS)'" -c calgo582.f" 1>>$f
fi
echo "" 1>>$f
echo "dbinom.o: dbinom.f "'\' 1>>$f
echo "        messpass.o" 1>>$f
if [ $PHAML_F90 = "xlf" ]
then
echo "	ln -s dbinom.f dbinom.f90" 1>>$f
echo "	"'$(F90) $(FFLAGS) -qfixed'" -c dbinom.f90" 1>>$f
echo "	rm -f dbinom.f90" 1>>$f
else
echo "	"'$(F90) $(FFLAGS)'" -c dbinom.f" 1>>$f
fi
echo "" 1>>$f

if [ -z "$LAPACKLIBS" ]
then
echo "lapack.o: lapack.f" 1>>$f
   if [ $PHAML_F90 = "xlf" ]
   then
echo "	ln -s lapack.f lapack.f90" 1>>$f
echo "	"'$(F90) $(FFLAGS) -qfixed'" -c lapack.f90" 1>>$f
echo "	rm -f lapack.f90" 1>>$f
   else
echo "	"'$(F90) $(FFLAGS)'" -c lapack.f" 1>>$f
   fi
echo "" 1>>$f
fi
if [ -z "$BLASLIBS" ]
then
echo "blas.o: blas.f" 1>>$f
   if [ $PHAML_F90 = "xlf" ]
   then
echo "	ln -s blas.f blas.f90" 1>>$f
echo "	"'$(F90) $(FFLAGS) -qfixed'" -c blas.f90" 1>>$f
echo "	rm -f blas.f90" 1>>$f
   else
echo "	"'$(F90) $(FFLAGS)'" -c blas.f" 1>>$f
   fi
echo "" 1>>$f
fi

case "$PHAML_GRAPHICS" in
   metro|mesa|opengl)
echo "graphics.o: graphics.f90 "'\' 1>>$f
echo "        global.o messpass.o hash.o rendereps.o modview.o gridtype.o evaluate.o phamltype.o" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" $OGLMODS $ZOLTANMOD -c graphics.f90" 1>>$f
echo "" 1>>$f
echo "modview.o: modview.f90" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" $OGLMODS -c modview.f90" 1>>$f
echo "" 1>>$f
echo "rendereps.o: rendereps.f90" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" $OGLMODS -c rendereps.f90" 1>>$f
echo "" 1>>$f ;;
   none)
      if [ $PHAML_PARALLEL = "messpass_nospawn" ]
      then
echo "graphics.o: graphics.dum.f90" 1>>$f
echo "	"'$(F90) $(FFLAGS)'" -o graphics.o -c graphics.dum.f90" 1>>$f
echo "" 1>>$f
      fi ;;
esac

if [ $PHAML_PARLIB = "lam" -a $F90 = "mpif77" ]
then
echo "showme:" 1>>$f
echo "	"'$(F90) --showme' 1>>$f
fi

echo "clean:" 1>>$f
echo "	rm -f *.o *.mod libphaml.a *.M *.stb" 1>>$f
