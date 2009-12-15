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
! This is a simple example that just solves one linear elliptic pde with
! periodic boundary conditions.
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
! Local variables

type(phaml_solution_type) :: soln
integer :: nproc
!----------------------------------------------------
! Begin executable code

nproc=1
!print *,"number of processors?"
!read *,nproc

call phaml_create(soln, nproc=nproc, triangle_files="period.1", &
                  draw_grid_who=MASTER  )

call phaml_solve_pde(soln,                   &
                     max_vert=1000,          &
                     mg_cycles=2,           &
                     refterm=DOUBLE_NVERT, &
                     pause_at_start=.true., &
                     pause_after_phases=.true., &
                     draw_grid_when=PHASES, &
                     print_grid_when=PHASES    , &
                     print_grid_who=MASTER,  &
                     print_linsys_when=PHASES, &
                     print_linsys_who=MASTER,  &
                     print_error_when=PHASES,&
                     print_error_what=ENERGY_LINF_ERR, &
                     print_errest_what=ENERGY_LINF_ERREST, &
                     print_error_who=MASTER)

call phaml_destroy(soln)

end program phaml_master
