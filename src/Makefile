#---------------------------------------------------------------------!
#                                PHAML                                !
#                                                                     !
# The Parallel Hierarchical Adaptive MultiLevel code for solving      !
# linear elliptic partial differential equations of the form          !
# (PUx)x + (QUy)y + RU = F on 2D polygonal domains with mixed         !
# boundary conditions, and eigenvalue problems where F is lambda*U.   !
#                                                                     !
# PHAML is public domain software.  It was produced as part of work   !
# done by the U.S. Government, and is not subject to copyright in     !
# the United States.                                                  !
#                                                                     !
#     William F. Mitchell                                             !
#     Mathematical and Computational Sciences Division                !
#     National Institute of Standards and Technology                  !
#     william.mitchell@nist.gov                                       !
#     http://math.nist.gov/phaml                                      !
#                                                                     !
#---------------------------------------------------------------------!

# Makefile created for system configuration:
#   Architecture:      x86
#   OS:                linux
#   F90 compiler:      gfortran
#   C compiler:        cc
#   Hash size:         1
#   Parallel form:     sequential
#   Parallel library:  none
#   Graphics:          none
#   BLAS:              source
#   LAPACK:            source
#   ARPACK:            no
#   BLOPEX:            no
#   hypre:             no
#   MUMPS:             no
#   PETSc:             no
#   SUPERLU:           no
#   Zoltan:            no
#   ParMETIS:          no
#   JOSTLE:            no
#   PaToH:             no
#   ParKway:           no
#   Nemesis:           no
#   DRUM:              no
#   Specific system:   none

phaml:
	cd src; make lib

clean:
	cd src; make clean
	cd testdir; make clean
	cd examples/all; make clean
	cd examples/curved; make clean
	cd examples/domains; make clean
	cd examples/eigenvalue; make clean
	cd examples/elliptic; make clean
	cd examples/L-domain; make clean
	cd examples/nonlinear; make clean
	cd examples/parabolic; make clean
	cd examples/periodic; make clean
	cd examples/rectangle; make clean
	cd examples/runtimepde; make clean
	cd examples/simple; make clean
	cd examples/system; make clean
	cd examples/system_ss; make clean

reallyclean:
	cd src; make clean
	cd testdir; make reallyclean
	cd examples/all; make clean
	cd examples/curved; make clean
	cd examples/domains; make clean
	cd examples/eigenvalue; make clean
	cd examples/elliptic; make clean
	cd examples/L-domain; make clean
	cd examples/nonlinear; make clean
	cd examples/parabolic; make clean
	cd examples/periodic; make clean
	cd examples/rectangle; make clean
	cd examples/runtimepde; make clean
	cd examples/simple; make clean
	cd examples/system; make clean
	cd examples/system_ss; make clean
	rm -f lib/*
	rm -f modules/*
	cd src; rm -f Makefile Makefile.bak
	cd testdir; rm -f Makefile Makefile.bak
	cd examples/all; rm -f Makefile Makefile.bak
	cd examples/curved; rm -f Makefile Makefile.bak
	cd examples/domains; rm -f Makefile Makefile.bak
	cd examples/eigenvalue; rm -f Makefile Makefile.bak
	cd examples/elliptic; rm -f Makefile Makefile.bak
	cd examples/L-domain; rm -f Makefile Makefile.bak
	cd examples/nonlinear; rm -f Makefile Makefile.bak
	cd examples/parabolic; rm -f Makefile Makefile.bak
	cd examples/periodic; rm -f Makefile Makefile.bak
	cd examples/rectangle; rm -f Makefile Makefile.bak
	cd examples/runtimepde; rm -f Makefile Makefile.bak
	cd examples/simple; rm -f Makefile Makefile.bak
	cd examples/system; rm -f Makefile Makefile.bak
	cd examples/system_ss; rm -f Makefile Makefile.bak

test:
	cd testdir; make test what=$(what)

