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
! Local variables:

integer :: nvert,nproc,nloop
type(phaml_solution_type) :: pde1

!----------------------------------------------------
! Begin executable code

! Use either number of vertices or number of refine/solve loops as
! termination criterion, and either set it here or read it at run time

nvert = huge(0)
!nvert = 20000
!print *,'number of vertices?'
!read *,nvert

!nloop = huge(0)
nloop = 7
!print *,'number of refine-solve loops?'
!read *,nloop

! Set the number of processors here, or read it at run time

nproc = 1
!print *,'number of processors?'
!read *,nproc

! create the pde object

call phaml_create(pde1,nproc,triangle_files="L.1", &
!                  spawn_form=DEBUG_SLAVE, debug_command="gdb", &
                  draw_grid_who = MASTER)

! Solve the pde

call phaml_solve_pde(pde1, &
                     max_vert=nvert, &
                     max_refsolveloop=nloop, &
                     print_grid_when=PHASES    ,print_grid_who=MASTER  ,&
                     print_error_when=PHASES   ,print_error_who=MASTER, &
                     print_error_what=ENERGY_LINF_ERR, &
                     print_errest_what=ENERGY_LINF_ERREST, &
	             draw_grid_when=PHASES     , &
                     pause_after_phases = .true.)

call phaml_destroy(pde1)

end program phaml_master
