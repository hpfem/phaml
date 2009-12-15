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
subroutine phaml_master(whodrawg,nslave)
!          ------------

!----------------------------------------------------
! This is the main program.
!----------------------------------------------------

!----------------------------------------------------
! Modules used are:

use phaml
use phaml_user_mod
!----------------------------------------------------

implicit none

! Dummy arguments

integer :: whodrawg, nslave
!----------------------------------------------------
! Local variables:

integer :: nvert,num_eval,print_err
real(my_real) :: lambda0
type(phaml_solution_type) :: pde1

!----------------------------------------------------
! Begin executable code

! Solve on a grid with 1000 vertices

nvert = 1000

! Find the eigenvalue closest to lambda0.  For the smallest eigenvalue,
! use -huge(0.0_my_real) for a more efficient solution.

lambda0 = -huge(0.0_my_real)
!print *,"lambda0?"
!read *,lambda0

! Set the number of eigenvalues to compute.  It will compute the
! num_eval eigenvalues closest to lambda0.

num_eval = 5
!print *,"number of eigenvalues?"
!read *,num_eval

! Input the parameter that determines the width of the well.
! The smallest eigenvalue is sqrt(prob_param*(2 sqrt(3))).
! 25 is a nice easy one (gives an eigenvalue of 9.6592583...)

prob_param = 25.0_my_real
!print *,"parameter that determines width of well (for example, 25.0)?"
!read *,prob_param

! Make sure prob_param is positive.

if (prob_param < 0.0_my_real) then
   print *,"the parameter must be positive; changing to ",abs(prob_param)
   prob_param = abs(prob_param)
endif

! The true solution is known for the smallest eigenvalue.  Print the error
! only if lambda0 is less than the smallest.  The error estimate, but not
! the actual error, will be printed for the other eigenvalues.

if (lambda0 < sqrt(2*prob_param)) then
   print_err = PHASES
else
   print_err = NEVER
endif

! Create the pde data structure

call phaml_create(pde1,nslave, &
                  eq_type = EIGENVALUE,       &
                  draw_grid_who = whodrawg)

! Send prob_param to the slaves

call update_usermod(pde1)

! Solve the problem

call phaml_solve_pde(pde1,                       &
                     max_vert=nvert,             &
                     sequential_vert=250,        &
                     print_grid_when=PHASES,     &
                     print_grid_who=MASTER,      &
                     print_error_when=print_err, &
                     print_error_who=MASTER,     &
                     print_error_what=ENERGY_LINF_ERR, &
                     print_errest_what=ENERGY_LINF_ERREST, &
                     print_eval_when=PHASES,     &
                     print_eval_who=MASTER,      &
                     draw_grid_when=PHASES,      &
                     mg_cycles = 10,             &
                     mg_tol=MG_NO_TOL,           &
                     inc_quad_order=3,           &
                     lambda0 = lambda0,          &
                     num_eval = num_eval)

! Destroy the pde data structure

call phaml_destroy(pde1)

end subroutine phaml_master
