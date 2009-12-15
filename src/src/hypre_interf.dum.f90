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

module hypre_interf

!----------------------------------------------------
! This module is a dummy version of the interface to hypre.
!----------------------------------------------------

!----------------------------------------------------
! Other modules used are:

use global
use message_passing
use hypretype_mod
use linsystype_mod
use gridtype_mod

!----------------------------------------------------
implicit none
private
public create_hypre_linear_system, destroy_hypre_linear_system, &
       change_hypre_rhs, zero_hypre_solution, hypre_solve, &
       hypre_lobpcg_solve_f
!----------------------------------------------------

contains

!          --------------------------
subroutine create_hypre_linear_system(hypre_matrix,phaml_matrix,procs, &
                                      still_sequential)
!          --------------------------
type(hypre_matrix_type), intent(inout) :: hypre_matrix
type(linsys_type), intent(inout) :: phaml_matrix
type(proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential
call warning("dummy hypre routine called")
end subroutine create_hypre_linear_system

!          ---------------------------
subroutine destroy_hypre_linear_system(hypre_matrix,procs)
!          ---------------------------
type(hypre_matrix_type), intent(inout) :: hypre_matrix
type(proc_info), intent(in) :: procs
call warning("dummy hypre routine called")
end subroutine destroy_hypre_linear_system

!          ----------------
subroutine change_hypre_rhs(hypre_matrix,phaml_matrix,rs,procs,still_sequential)
!          ----------------
type(hypre_matrix_type), intent(inout) :: hypre_matrix
type(linsys_type), intent(in) :: phaml_matrix
real(my_real), intent(in) :: rs(:)
type(proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential
call warning("dummy hypre routine called")
end subroutine change_hypre_rhs

!          -------------------
subroutine zero_hypre_solution(hypre_matrix,procs,still_sequential)
!          -------------------
type(hypre_matrix_type), intent(inout) :: hypre_matrix
type(proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential
call warning("dummy hypre routine called")
end subroutine zero_hypre_solution

!          -----------
subroutine hypre_solve(hypre_matrix,phaml_matrix,procs,solver_cntl, &
                       still_sequential)
!          -----------
type(hypre_matrix_type), intent(inout) :: hypre_matrix
type(linsys_type), intent(inout) :: phaml_matrix
type(solver_options), intent(in) :: solver_cntl
type(proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential
call warning("dummy hypre routine called")
end subroutine hypre_solve

!          --------------------
subroutine hypre_lobpcg_solve_f(phaml_matrix,full_matrix,hypre_matrix, &
                                solver_cntl,io_cntl,still_sequential,grid,procs)
!          --------------------

type(linsys_type), intent(inout), target :: phaml_matrix, full_matrix
type(hypre_matrix_type), intent(inout), target :: hypre_matrix
type(solver_options), intent(in), target :: solver_cntl
type(io_options), intent(in), target :: io_cntl
logical, intent(in) :: still_sequential
type(grid_type), intent(in), target :: grid
type(proc_info), intent(in), target :: procs

end subroutine hypre_lobpcg_solve_f

end module hypre_interf
