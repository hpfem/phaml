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

integer :: nvert,nloop
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

! create the pde object

call phaml_create(pde1,nproc, &
                  triangle_files="L.1", &
                  draw_grid_who = whodrawg)

! Solve the pde

call phaml_solve_pde(pde1, &
                     max_vert=nvert, &
                     max_refsolveloop=nloop, &
                     print_grid_when=PHASES    ,print_grid_who=MASTER  ,&
                     print_error_when=PHASES   ,print_error_who=MASTER, &
                     print_error_what=ENERGY_LINF_ERR, &
                     print_errest_what=ENERGY_LINF_ERREST, &
                     draw_grid_when=PHASES    , &
                     pause_after_phases = .false.)

call phaml_destroy(pde1)

end subroutine phaml_master
