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

nloop=1
nloop2=3

nproc = 2

! create the phaml_solution variables

allocate(pde(2))
call phaml_create(pde(1),nproc,id = 1)
call phaml_create(pde(2),nproc,id = 2)

! connect the two solutions

call phaml_connect(1,2)

! set pde(2) to a zero function so it has something when pde(1) tries
! to evaluate it.  The zero is in subroutine icond.

call phaml_solve_pde(pde(2),                    &
                     max_vert=1,                &
                     reftype=H_UNIFORM,         &
                     error_estimator=INITIAL_CONDITION, &
                     task=SET_INITIAL,          &
                     print_header_who=NO_ONE)

! loop to alternate between solving the two equations

refterm = DOUBLE_NVERT

do i=1,nloop2
   write(6,"(A,I11)") "solving pde 1, iteration ",i
   call phaml_solve_pde(pde(1), &
                        max_refsolveloop=nloop, &
                        sequential_vert=20,     &
                        print_grid_when=FINAL,  &
                        print_grid_who=MASTER,  &
                        print_error_when=FINAL, &
                        print_error_who=MASTER, &
                        print_error_what=ENERGY_LINF_ERR, &
                        print_errest_what=ENERGY_LINF_ERREST, &
                        refterm=refterm,        &
                        mg_cycles=2,            &
                        print_header_who=NO_ONE,&
                        print_trailer_who=NO_ONE)

   write(6,"(A,I11)") "solving pde 2, iteration ",i
   call phaml_solve_pde(pde(2), &
                        max_refsolveloop=nloop, &
                        sequential_vert=20,     &
                        print_grid_when=FINAL,  &
                        print_grid_who=MASTER,  &
                        print_error_when=FINAL, &
                        print_error_who=MASTER, &
                        print_error_what=ENERGY_LINF_ERR, &
                        print_errest_what=ENERGY_LINF_ERREST, &
                        refterm=refterm,        &
                        mg_cycles=2,            &
                        print_header_who=NO_ONE,&
                        print_trailer_who=NO_ONE)

end do

call phaml_destroy(pde(1),finalize_mpi=.false.)
call phaml_destroy(pde(2))
deallocate(pde)

end program phaml_master
