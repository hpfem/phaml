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

module superlutype_mod

!----------------------------------------------------
! This module contains a data structure to hold the SuperLU form of a matrix.
!----------------------------------------------------

use global

implicit none
private
public superlu_matrix_type, superlu_ptr, dp

!----------------------------------------------------

integer, parameter :: dp = kind(0.0d0)

! kind of integer to hold a SuperLU pointer.  Use default integer.
! This might need to be changed on systems with large memory.
! If changed, be sure to change it in superlu_wrappers.c too.

integer, parameter :: superlu_ptr = kind(0)

! The SuperLU matrix data type

type superlu_matrix_type
   integer(superlu_ptr) :: grid            ! type gridinfo_t
   integer(superlu_ptr) :: options         ! type superlu_options_t
   integer(superlu_ptr) :: ScalePermstruct ! type ScalePermstruct_t
   integer(superlu_ptr) :: LUstruct        ! type LUstruct_t
   integer(superlu_ptr) :: SOLVEstruct     ! type SOLVEstruct_t
   integer(superlu_ptr) :: A               ! type SuperMatrix
   integer :: my_neq
   integer, pointer :: firsteq(:)
   integer, pointer :: global_eq(:)
   integer, pointer :: rowptr(:)
   integer, pointer :: colind(:)
   real(dp), pointer :: matval(:)
   real(dp), pointer :: b(:,:)
   real(my_real), pointer :: lapack_rhs(:)
end type superlu_matrix_type

end module superlutype_mod
