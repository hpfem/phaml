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
! This is a simple example that just solves one linear elliptic pde.
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
character(len=32) :: filename
!----------------------------------------------------
! Begin executable code

! set mesh file name or read it from a prompt

filename = "superior.1"
print *,"Enter the root of the mesh file names (include .1 if they contain it)"
print *,"for example, superior.1"
read *, filename

call phaml_create(soln,nproc=1, &
                  triangle_files=filename, &
                  draw_grid_who=MASTER  )

call phaml_solve_pde(soln,                   &
                     max_vert=10000,          &
                     draw_grid_when = PHASES, &
                     pause_after_phases = .true., &
                     pause_at_start = .true., &
                     print_error_who = MASTER, &
                     print_error_when = PHASES, &
                     print_error_what=ENERGY_LINF_ERR, &
                     print_errest_what=ENERGY_LINF_ERREST, &
                     print_grid_when=PHASES    , &
                     print_grid_who=MASTER)

call phaml_destroy(soln)

end program phaml_master
