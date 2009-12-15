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
! This version solves a parabolic equation.  See README.parabolic for
! a description of the implicit time stepping scheme.
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

use global
use phaml
use phaml_user_mod
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Local variables:

integer :: nproc,nvert
real(my_real) :: finalt
type(phaml_solution_type) :: soln

!----------------------------------------------------
! Begin executable code

nvert = 500
finalt = .05_my_real
deltat = .01_my_real

nproc = 4

! create the phaml_solution variables

t = 0.0_my_real
call phaml_create(soln,nproc,update_umod=.true.)

! set the initial condition

call phaml_solve_pde(soln,                      &
                     max_vert=nvert,            &
                     task=SET_INITIAL,          &
                     refterm=ONE_REF_HALF_ERRIND, &
                     error_estimator=INITIAL_CONDITION, &
                     print_header_who=NO_ONE,   &
                     print_trailer_who=NO_ONE)

! solve the equation until the final time is reached.

do
   t = t + deltat
   if (t > finalt) exit
   call update_usermod(soln)
   call phaml_copy_soln_to_old(soln)
   call phaml_solve_pde(soln,                    &
                        max_refsolveloop=1,      &
                        refterm=KEEP_NVERT,      &
                        mg_cycles=1,             &
                        print_grid_when=FINAL,   &
                        print_grid_who=MASTER,   &
                        print_error_when=FINAL,  &
                        print_error_who=MASTER,  &
                        print_errest_what=ENERGY_LINF_ERREST, &
                        print_header_who=NO_ONE, &
                        print_trailer_who=NO_ONE)
   write(6,"(A,SS,1P,E19.12E2)") "time = ",t
end do

call phaml_destroy(soln)

end program phaml_master
