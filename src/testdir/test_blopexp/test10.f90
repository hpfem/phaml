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

!----------------------------------------------------
! This file contains the main program supplied by the user.
! This version illustrates an eigenvalue problem, the Schroedinger
! Equation for a simple harmonic oscillator.
!----------------------------------------------------

!       ------------
program phaml_master
!       ------------

!----------------------------------------------------
! This is the main program.
!----------------------------------------------------

!----------------------------------------------------
! Modules used are:

use phaml
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Local variables:

integer :: nvert,num_eval,nproc
real(my_real) :: lambda0
type(phaml_solution_type) :: pde1

!----------------------------------------------------
! Begin executable code

nvert = 1000
lambda0 = 30.0_my_real
num_eval = 3
nproc=1
call phaml_create(pde1,nproc,eq_type=EIGENVALUE)
call phaml_solve_pde(pde1,                     &
                     max_vert=nvert,             &
                     sequential_vert=250,        &
                     print_grid_when=PHASES,     &
                     print_grid_who=MASTER,      &
                     print_error_when=PHASES,    &
                     print_error_who=MASTER,     &
                     print_errest_what=ENERGY_LINF_ERREST, &
                     print_eval_when=PHASES,     &
                     print_eval_who=MASTER,      &
                     print_warnings=.false.,     &
                     eigensolver=BLOPEX_SOLVER,  &
                     solver=PETSC_GMRES_SOLVER,  &
                     preconditioner=PETSC_SOR_PRECONDITION, &
                     lambda0_side=EIGEN_RIGHT,   &
                     transformation=SHIFT_INVERT,&
                     lambda0 = lambda0,          &
                     num_eval = num_eval)
call phaml_destroy(pde1)

end program phaml_master
