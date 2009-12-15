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
! This example illustrates the use of all the procedures in the
! PHAML interface.
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
! Local variables

type(phaml_solution_type) :: soln
real(my_real), dimension(36) :: x, y, u
real(my_real) :: max_norm_err, soln_intgrl
integer :: i,j,nvert
!----------------------------------------------------
! Begin executable code

! create soln for solution using 2 processors

call phaml_create(soln,nproc=2)

! solve it, using 1000 vertices

call phaml_solve_pde(soln,                   &
                     max_vert=1000,          &
                     mg_cycles=5,            &
                     print_grid_when=PHASES, &
                     print_grid_who=MASTER,  &
                     print_error_when=PHASES,&
                     print_error_what=LINF_ERR, &
                     print_errest_what=LINF_ERREST, &
                     print_error_who=MASTER)

! evaluate the solution on a 6X6 grid

x = (/ ((i/5.0_my_real, i=0,5), j=0,5) /)
y = (/ ((i/5.0_my_real, j=0,5), i=0,5) /)
call phaml_evaluate(soln,x,y,u)
print *, "x, y, and solution for a 6X6 grid"
do i=1,36
   print *,x(i),y(i),u(i)
end do

! open output files to save the current state of soln

call phaml_popen(soln,unit=12,file="example_save.dat",form="FORMATTED")

! compress unrefined elements out of the data structures to reduce the
! size of the save file

call phaml_compress(soln)

! save soln

call phaml_store(soln,12)

! close the files

call phaml_pclose(soln,12)

! change the value of global variables in module phaml_user_mod and
! update them on the slaves with update_usermod.  (This doesn't actually
! affect this program, but is just here to show how to do it.)

iglobal1 = 10
iglobal2 = 20
realvar = 1.0_my_real
call update_usermod(soln)

! destroy soln

call phaml_destroy(soln,finalize_mpi=.false.)

! create soln again and read the saved data

call phaml_create(soln, nproc=2)
call phaml_popen(soln,13,"example_save.dat",form="FORMATTED")
call phaml_restore(soln,13)
call phaml_pclose(soln,13)

! get and print some information about soln

call phaml_query(soln, nvert=nvert,  linf_error=max_norm_err)
print *
print *,"grid has ",nvert," vertices"
print *,"max norm of error is ",max_norm_err

! refine the grid once, just to do something with the restored data

call phaml_solve_pde(soln,                        &
                     max_refsolveloop=1,          &
                     mg_cycles=5,                 &
                     print_grid_when=PHASES,      &
                     print_grid_who=MASTER,       &
                     print_error_when=PHASES,     &
                     print_error_who=MASTER,      &
                     print_error_what=ENERGY_LINF_ERR, &
                     print_errest_what=ENERGY_LINF_ERREST, &
                     print_header_who=NO_ONE,     &
                     print_trailer_who=NO_ONE     )

! compute the integral of the square of the solution, and scale the solution
! so that its L2 norm is 1.0

soln_intgrl = phaml_integrate(soln,1,p=2)
print *,"L2 norm of solution ",sqrt(soln_intgrl)
call phaml_scale(soln,1.0_my_real/sqrt(soln_intgrl))

! verify that the normalized solution has norm 1.0

soln_intgrl = phaml_integrate(soln,1,p=2)
print *,"L2 norm of normalized solution ",sqrt(soln_intgrl)

call phaml_destroy(soln)

! using subroutine phaml_connect does not make sense in this example.  See
! example system_ss for an example of using phaml_connect.

! using subroutines phaml_copy_soln_to_old and phaml_evaluate_old do not
! make sense in this example.  See example parabolic for an example of
! using them (phaml_copy_to_old is in master.f90 and phaml_evaluate_old
! is in pde.f90).

end program phaml_master
