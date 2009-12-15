!---------------------------------------------------------------------!
!                                PHAML                                !
!                                                                     !
! The Parallel Hierarchical Adaptive MultiLevel code for solving      !
! linear elliptic partial differential equations of the form          !
! (PUx)x + (QUy)y + RU = F on 2D polygonal domains with mixed         !
! boundary conditions, and eigenvalue problems where F is lambda*U.   !
!                                                                     !
! PHAML is public domain software.  It was produced as part of work   !
! done by the U.S. Government, and is not subject to copyright in     !
! the United States.                                                  !
!                                                                     !
!     William F. Mitchell                                             !
!     Mathematical and Computational Sciences Division                !
!     National Institute of Standards and Technology                  !
!     william.mitchell@nist.gov                                       !
!     http://math.nist.gov/phaml                                      !
!                                                                     !
!---------------------------------------------------------------------!

module petsc_init_mod

!----------------------------------------------------
! This module contains dummy routines corresponding to those in
! petsc_init.F90 to satisfy references when PETSc is not used.
!----------------------------------------------------

implicit none
private
public petsc_init, petsc_finalize

contains

!          ----------
subroutine petsc_init(slaves_communicator)
!          ----------
!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: slaves_communicator
!----------------------------------------------------
! Local variables:

integer :: jerr
!----------------------------------------------------
! Begin executable code

jerr = slaves_communicator ! just to shut up picky compilers

end subroutine petsc_init

!          --------------
subroutine petsc_finalize
!          --------------

end subroutine petsc_finalize

end module petsc_init_mod
