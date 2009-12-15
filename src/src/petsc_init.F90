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
! This module contains routines that initialize and finalize PETSc.
! They cannot be in the petsc_interf module because petsc_interf
! uses message_passing, and message_passing calls these routines.
!----------------------------------------------------

!----------------------------------------------------
! Other modules used are:

use global
!----------------------------------------------------

implicit none
private
public petsc_init, petsc_finalize

!----------------------------------------------------
! The PETSc include files.  Note the use of preprocessor #include instead of
! the Fortran include statement, because the include files contain
! preprocessor directives.

#include "include/finclude/petsc.h"

!----------------------------------------------------
! The following parameters are defined:

!----------------------------------------------------

!----------------------------------------------------
! The following types are defined:

!----------------------------------------------------

!----------------------------------------------------
! The following variables are defined:

!----------------------------------------------------

contains

!          ----------
subroutine petsc_init(slaves_communicator)
!          ----------

!----------------------------------------------------
! This routine initializes PETSc
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: slaves_communicator
!----------------------------------------------------
! Local variables:

integer :: jerr
MPI_Comm petsc_world_communicator
!----------------------------------------------------
! Begin executable code

petsc_world_communicator = slaves_communicator

! For PETSc versions before 2.3.0
!call PetscSetCommWorld(slaves_communicator,jerr)

! For PETSc versons 2.3.0 and later
PETSC_COMM_WORLD = slaves_communicator

call PetscInitialize(PETSC_NULL_CHARACTER,jerr)

end subroutine petsc_init

!          --------------
subroutine petsc_finalize
!          --------------

!----------------------------------------------------
! This routine finalizes PETSc
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

!----------------------------------------------------
! Local variables:

integer :: jerr
!----------------------------------------------------
! Begin executable code

call PetscFinalize(jerr)

end subroutine petsc_finalize

end module petsc_init_mod
