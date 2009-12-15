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
! This version solves a nonlinear equation by Picard iteration.
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

integer :: nproc,nvert,i,j
real(my_real), dimension(101*101) :: x, y, u, uold
type(phaml_solution_type) :: soln

!----------------------------------------------------
! Begin executable code

print *,"number of vertices for the grid (for example, 1000)?"
read *,nvert
nproc=1
!print *,"number of processors?"
!read *,nproc

! create the phaml_solution variable

call phaml_create(soln,nproc,draw_grid_who = MASTER)

! set the initial condition

call phaml_solve_pde(soln,                      &
                     max_vert=nvert,            &
                     task=SET_INITIAL,          &
                     refterm=DOUBLE_NEQ,        &
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

! give time to adjust the graphics view

print *,"adjust graphics view and press return"
read *

! solve the equation and evaluate the change in solution after each iteration

do i=1,6
   call phaml_copy_soln_to_old(soln)
   call phaml_solve_pde(soln,                    &
                        max_refsolveloop=1,      &
                        refterm=KEEP_NVERT,      &
                        mg_cycles=10,            &
                        draw_grid_when=FINAL,    &
                        print_header_who=NO_ONE, &
                        print_trailer_who=NO_ONE)
   uold = u
   call phaml_evaluate(soln,x,y,u)
   print *,"change in solution ",maxval(abs(u-uold))
end do

print *,"press return to terminate program"
read *

call phaml_destroy(soln)

end program phaml_master
