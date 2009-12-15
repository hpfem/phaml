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

module petsc_interf

!----------------------------------------------------
! This module contains dummy routines corresponding to those in
! petsc_interf.F90 to satisfy references when PETSc is not used.
!----------------------------------------------------

!----------------------------------------------------
! Other modules used are:

use global
use message_passing
use gridtype_mod
use petsc_type_mod
use linsystype_mod

!----------------------------------------------------

implicit none
private
public create_petsc_linear_system, create_petsc_linear_system_mf, &
       change_petsc_rhs, petsc_solve, destroy_petsc_linear_system, &
       petsc_lobpcg_solve_f

contains

!          --------------------------
subroutine create_petsc_linear_system(phaml_matrix,petsc_matrix, &
                                      still_sequential,procs)
!          --------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(in) :: phaml_matrix
type(petsc_matrix_type), intent(out), target :: petsc_matrix
logical, intent(in) :: still_sequential
type(proc_info), intent(in) :: procs
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

call warning("Dummy version of PETSc routine called.  Must compile with PETSc.")

! just to shut up picky compilers
petsc_matrix%nada=0

end subroutine create_petsc_linear_system

!          -----------------------------
subroutine create_petsc_linear_system_mf(phaml_matrix,petsc_matrix, &
                                         still_sequential)
!          -----------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(in) :: phaml_matrix
type(petsc_matrix_type), intent(out) :: petsc_matrix
logical, intent(in) :: still_sequential
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

call warning("Dummy version of PETSc routine called.  Must compile with PETSc.")

! just to shut up picky compilers
petsc_matrix%nada=0

end subroutine create_petsc_linear_system_mf

!          ----------------
subroutine change_petsc_rhs(rhs,petsc_matrix)
!          ----------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: rhs(:)
type(petsc_matrix_type), intent(inout) :: petsc_matrix
!----------------------------------------------------
! Local variables:

real(my_real) :: rvar
!----------------------------------------------------
! Begin executable code

call warning("Dummy version of PETSc routine called.  Must compile with PETSc.")

rvar = rhs(1)

end subroutine change_petsc_rhs

!          -----------
subroutine petsc_solve(phaml_matrix,petsc_matrix,solver_cntl,io_cntl, &
                       still_sequential,grid,procs)
!          -----------
!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(inout), target :: phaml_matrix
type(petsc_matrix_type), intent(inout), target :: petsc_matrix
type(solver_options), intent(in), target :: solver_cntl
type(io_options), intent(in), target :: io_cntl
logical, intent(in) :: still_sequential
type(grid_type), intent(in), target :: grid
type(proc_info), intent(in), target :: procs
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

call warning("Dummy version of PETSc routine called.  Must compile with PETSc.")

end subroutine petsc_solve

!          ---------------------------
subroutine destroy_petsc_linear_system(petsc_matrix)
!          ---------------------------

type(petsc_matrix_type), intent(inout) :: petsc_matrix

end subroutine destroy_petsc_linear_system

!          --------------------
subroutine petsc_lobpcg_solve_f(phaml_matrix,full_matrix,petsc_matrix, &
                                solver_cntl,io_cntl,still_sequential,grid,procs)
!          --------------------

type(linsys_type), intent(inout), target :: phaml_matrix, full_matrix
type(petsc_matrix_type), intent(inout), target :: petsc_matrix
type(solver_options), intent(in), target :: solver_cntl
type(io_options), intent(in), target :: io_cntl
logical, intent(in) :: still_sequential
type(grid_type), intent(in), target :: grid
type(proc_info), intent(in), target :: procs

end subroutine petsc_lobpcg_solve_f

end module petsc_interf
