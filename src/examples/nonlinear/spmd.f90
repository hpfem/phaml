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
! This file contains a user main program for spawnless message passing mode.
!
! NOTE -- in the spawnless mode you can only call create_pde once on
!         each processor.  The number of slaves is determined by
!         the number of processes (usually specified on the command
!         line) and graphics options.  For example, for 2 slaves with
!         draw_grid_who=MASTER, there should be 4 processes (master,
!         one graphics and two slaves).
!----------------------------------------------------

!       ---------------
program phaml_spawnless
!       ---------------

!----------------------------------------------------
! Modules used are:

use phaml
use mpif_mod
implicit none

!----------------------------------------------------
! Local variables:

integer :: whodrawg
integer :: jerr
integer :: my_proc, total_nproc
integer :: nslave, subtract, divide

!----------------------------------------------------
! Begin executable code

! set the graphics options

whodrawg = MASTER

! initialize MPI, find out how many processors and what my rank is

call mpi_init(jerr)
call mpi_comm_size(MPI_COMM_WORLD,total_nproc,jerr)
call mpi_comm_rank(MPI_COMM_WORLD,my_proc,jerr)

! determine how many processors for slaves and graphics

subtract = 1
if (whodrawg == MASTER .or. whodrawg == EVERYONE) subtract = 2
divide = 1
if (whodrawg == SLAVES .or. whodrawg == EVERYONE) divide = 2
nslave = (total_nproc-subtract)/divide

! call the master, slave or graphics program depending on my rank

if (my_proc == 0) then
   call phaml_master(whodrawg,nslave)
elseif (my_proc <= nslave) then
   call phaml_slave
else
   call phaml_graphics
endif

stop
end program phaml_spawnless

!          ------------
subroutine phaml_master(whodrawg,nproc)
!          ------------

!----------------------------------------------------
! This is the main program for a master
!----------------------------------------------------

!----------------------------------------------------
! Modules used are:

use phaml
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

integer :: whodrawg, nproc
!----------------------------------------------------
! Local variables

integer :: nvert,i,j
real(my_real), dimension(101*101) :: x, y, u, uold
type(phaml_solution_type) :: soln

!----------------------------------------------------
! Begin executable code

! Use either number of vertices or number of refine/solve loops as
! termination criterion, and either set it here or read it at run time

nvert = 1000
!print *,'number of vertices?'
!read *,nvert

! create the pde object

call phaml_create(soln,nproc, &
                  draw_grid_who = whodrawg)

! set the initial condition

call phaml_solve_pde(soln,                      &
                     max_vert=nvert,            &
                     task=SET_INITIAL,          &
                     error_estimator=INITIAL_CONDITION, &
                     print_header_who=NO_ONE,   &
                     print_trailer_who=NO_ONE,  &
                     draw_grid_when=FINAL)

! set up a 101X101 grid on which to evaluate the solution to see how
! much the solution changes each iteration

x = (/ ((i/100.0_my_real, i=0,100), j=0,100) /)
y = (/ ((i/100.0_my_real, j=0,100), i=0,100) /)

! evaluate the initial guess

call phaml_evaluate(soln,x,y,u)

! solve the equation and evaluate the change in solution after each iteration

do i=1,6
   call phaml_copy_soln_to_old(soln)
   call phaml_solve_pde(soln,                    &
                        max_refsolveloop=1,      &
                        refterm=KEEP_NVERT,      &
                        mg_cycles=5,             &
                        draw_grid_when=FINAL,    &
                        print_header_who=NO_ONE, &
                        print_trailer_who=NO_ONE)
   uold = u
   call phaml_evaluate(soln,x,y,u)
   print *,"change in solution ",maxval(abs(u-uold))
end do

call phaml_destroy(soln)

end subroutine phaml_master
