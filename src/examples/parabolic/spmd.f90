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
integer :: my_processor, total_nproc
integer :: nslave, subtract, divide

!----------------------------------------------------
! Begin executable code

! set the graphics options

whodrawg = MASTER

! initialize MPI, find out how many processors and what my rank is

call mpi_init(jerr)
call mpi_comm_size(MPI_COMM_WORLD,total_nproc,jerr)
call mpi_comm_rank(MPI_COMM_WORLD,my_processor,jerr)

! determine how many processors for slaves and graphics

subtract = 1
if (whodrawg == MASTER .or. whodrawg == EVERYONE) subtract = 2
divide = 1
if (whodrawg == SLAVES .or. whodrawg == EVERYONE) divide = 2
nslave = (total_nproc-subtract)/divide

! call the master, slave or graphics program depending on my rank

if (my_processor == 0) then
   call phaml_master(whodrawg,nslave)
elseif (my_processor <= nslave) then
   call phaml_slave
else
   call phaml_graphics
endif

stop
end program phaml_spawnless

!          ------------
subroutine phaml_master(whodrawg,nslave)
!          ------------

!----------------------------------------------------
! This is the main program for a master
!----------------------------------------------------

!----------------------------------------------------
! Modules used are:

use global
use phaml
use phaml_user_mod
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

integer :: whodrawg, nslave
!----------------------------------------------------
! Local variables

integer :: nelem
real(my_real) :: finalt
type(phaml_solution_type) :: soln

!----------------------------------------------------
! Begin executable code

nelem = 256
!print *,"number of vertices for the grid (for example, 100)?"
!read *,nvert
finalt = 1.0_my_real
!print *,"final time (for example, 1.0)?"
!read *,finalt
deltat = .01_my_real
!print *,"time step (for example, .01)?"
!read *,deltat

! create the phaml_solution variables

t = 0.0_my_real
call phaml_create(soln,nslave,update_umod=.true., &
                  draw_grid_who = whodrawg)

! set the initial condition

call phaml_solve_pde(soln,                      &
                     max_elem=nelem,            &
                     task=SET_INITIAL,          &
                     error_estimator=INITIAL_CONDITION, &
                     print_header_who=NO_ONE,   &
                     print_trailer_who=NO_ONE,  &
                     draw_grid_when=FINAL)

! solve the equation until the final time is reached.

do
   t = t + deltat
   if (t > finalt) exit
   call update_usermod(soln)
   call phaml_copy_soln_to_old(soln)
   call phaml_solve_pde(soln,                    &
                        max_refsolveloop=1,      &
                        refterm=KEEP_NELEM,      &
                        mg_cycles=1,             &
                        draw_grid_when=FINAL ,   &
                        print_header_who=NO_ONE, &
                        print_trailer_who=NO_ONE)
   print *,"time = ",t
end do

call phaml_destroy(soln)

end subroutine phaml_master
