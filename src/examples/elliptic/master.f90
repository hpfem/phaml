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
use phaml_user_mod
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Local variables:

integer :: nvert,nproc,nloop
type(phaml_solution_type) :: pde1

!----------------------------------------------------
! Begin executable code

! Use either number of vertices or number of refine/solve loops as
! termination criterion, and either set it here or read it at run time

nvert = huge(0)
nvert = 5000
!print *,'number of vertices?'
!read *,nvert

nloop = huge(0)
!nloop = 13
!print *,'number of refine-solve loops?'
!read *,nloop

! Set the number of processors here, or read it at run time

nproc = 2
!print *,'number of processors?'
!read *,nproc

! If solving one of the examples that has a parameter, set it here or read it

prob_param = -20.0_my_real
!print *,"prob param"
!read *,prob_param

! create the pde object

call phaml_create(pde1,nproc, &
                  draw_grid_who = MASTER)

! Make sure the slaves have the problem parameter

call update_usermod(pde1)

! Solve the pde

call phaml_solve_pde(pde1, &
                     max_vert=nvert, &
                     max_refsolveloop=nloop, &
                     mg_cycles=5,  &
                     print_grid_when=PHASES    ,print_grid_who=MASTER  ,&
                     print_error_when=PHASES    ,print_error_who=MASTER, &
                     print_error_what=ENERGY_LINF_ERR, &
	             print_time_when=NEVER     ,print_time_who=MASTER    , &
	             draw_grid_when=PHASES    , &
                     pause_after_phases = .true.  ,&
                     pause_at_start = .false., &
                     pause_at_end=.false.)

call phaml_destroy(pde1)

end program phaml_master
