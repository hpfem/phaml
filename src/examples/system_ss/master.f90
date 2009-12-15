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
!
! This version solves a system of two elliptic pdes.
!
! There are two strategies implemented in this example:
!
! 1) alternate between solving the two equations while refining the grids.
!    For this one, "loops for the first solution" determines how far the
!    grid is refined before starting the alternating process.  It should
!    be small, for example 4.  "loops after that" determines how many
!    times to alternate between the equations, while doubling the size
!    of the grid each time.  It should be fairly large, for example 10.

! 2) alternate between solving the two equations after refining the grids
!    to their final size.  For this one, "loops for the first solution"
!    determines the size of the grids.  It should be fairly large, for
!    example 12.  "loops after that" determines the number of times to
!    alternate between the equations on the final grids.  It should be
!    fairly small (for this example where convergence of this iteration
!    is fast), for example 5.
!
! Which strategy is used is determined by the value of refterm after the
! first solution.  If it is DOUBLE_NVERT, then the grid size is doubled in each
! solution (strategy 1).  If it is KEEP_NVERT, then the grid size stays
! the same in each subsequent solution.  Look for where refterm is changed at
! the end of the main loop to select the stragegy.
!
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

integer :: nproc,nloop,nloop2,i,init_form,refterm

!----------------------------------------------------
! Begin executable code

print *,'number of refine-solve loops for first solution?'
read *,nloop
print *,'number of refine-solve loops after that?'
read *,nloop2

nproc = 1
print *,'number of processors?'
read *,nproc

! create the phaml_solution variables

allocate(pde(2))
call phaml_create(pde(1),nproc, &
                  draw_grid_who = MASTER,id = 1)
call phaml_create(pde(2),nproc, &
                  draw_grid_who = MASTER,id = 2)

! connect the two solutions

call phaml_connect(1,2)

! set pde(2) to a zero function so it has something when pde(1) tries
! to evaluate it.  The zero is in subroutine icond.

call phaml_solve_pde(pde(2),                    &
                     max_vert=1,                &
                     refterm=DOUBLE_NVERT,      &
                     task=SET_INITIAL,          &
                     error_estimator=INITIAL_CONDITION, &
                     draw_grid_when=FINAL,      &
                     print_header_who=NO_ONE)

! loop to alternate between solving the two equations

refterm = DOUBLE_NVERT

do i=1,nloop2
   print *,"solving pde 1, iteration ",i
   call phaml_solve_pde(pde(1), &
                        max_refsolveloop=nloop, &
                        print_grid_when=FINAL,  &
                        print_grid_who=MASTER,  &
                        print_error_when=FINAL, &
                        print_error_who=MASTER, &
                        print_error_what=ENERGY_LINF_ERR, &
                        print_errest_what=ENERGY_LINF_ERREST, &
	                draw_grid_when=PHASES,  &
                        refterm=refterm,        &
                        mg_cycles=2,            &
                        print_header_who=NO_ONE,&
                        print_trailer_who=NO_ONE,&
                        pause_at_start = .true.)

   print *,"solving pde 2, iteration ",i
   call phaml_solve_pde(pde(2), &
                        max_refsolveloop=nloop, &
                        print_grid_when=FINAL,  &
                        print_grid_who=MASTER,  &
                        print_error_when=FINAL, &
                        print_error_who=MASTER, &
                        print_error_what=ENERGY_LINF_ERR, &
                        print_errest_what=ENERGY_LINF_ERREST, &
	                draw_grid_when=PHASES,  &
                        refterm=refterm,        &
                        mg_cycles=2,            &
                        print_header_who=NO_ONE,&
                        print_trailer_who=NO_ONE,&
                        pause_at_start = .true.)

   nloop = 1
! Reset refterm to KEEP_NVERT here to use strategy 2, or
! leave it alone for strategy 1.
!   refterm = KEEP_NVERT
end do

call phaml_destroy(pde(1),finalize_mpi=.false.)
call phaml_destroy(pde(2))
deallocate(pde)

end program phaml_master
