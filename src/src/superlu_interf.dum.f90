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

module superlu_interf

!----------------------------------------------------
! This module contains a dummhy version of the interface to the SuperLU package.
!----------------------------------------------------

!----------------------------------------------------
! Other modules used are:

use global
use message_passing
use superlutype_mod
use linsystype_mod

!----------------------------------------------------
implicit none
private
public create_superlu_linear_system, destroy_superlu_linear_system, &
       change_superlu_rhs, superlu_solve

contains

!          ----------------------------
subroutine create_superlu_linear_system(superlu_matrix, phaml_matrix, &
                                        procs, still_sequential)
!          ----------------------------
type(superlu_matrix_type), intent(inout) :: superlu_matrix
type(linsys_type), intent(inout) :: phaml_matrix
type(proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential

call warning("Dummy version of SuperLU interface routine called.")

end subroutine create_superlu_linear_system

!          -------------
subroutine superlu_solve(superlu_matrix,phaml_matrix,procs, &
                         still_sequential)
!          -------------
type(superlu_matrix_type), intent(inout) :: superlu_matrix
type(linsys_type), intent(inout) :: phaml_matrix
type(proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential

call warning("Dummy version of SuperLU interface routine called.")

end subroutine superlu_solve

!          -----------------------------
subroutine destroy_superlu_linear_system(superlu_matrix,procs)
!          -----------------------------
type(superlu_matrix_type), intent(inout) :: superlu_matrix
type(proc_info), intent(in) :: procs

call warning("Dummy version of SuperLU interface routine called.")

end subroutine destroy_superlu_linear_system

!          ------------------
subroutine change_superlu_rhs(superlu_matrix,phaml_matrix,rs,procs, &
                              still_sequential)
!          ------------------
type(superlu_matrix_type), intent(inout) :: superlu_matrix
type(linsys_type), intent(in) :: phaml_matrix
real(my_real), intent(in) :: rs(:)
type(proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential

call warning("Dummy version of SuperLU interface routine called.")

end subroutine change_superlu_rhs

end module superlu_interf
