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

module mumps_interf

!----------------------------------------------------
! This module contains a dummy version of the interface to the mumps package.
!----------------------------------------------------

!----------------------------------------------------
! Other modules used are:

use global
use message_passing
use mumps_struc_mod
use linsystype_mod

!----------------------------------------------------
implicit none
private
public create_mumps_linear_system, destroy_mumps_linear_system, &
       change_mumps_rhs, mumps_solve

!----------------------------------------------------

contains

!          --------------------------
subroutine create_mumps_linear_system(mumps_matrix,phaml_matrix,procs, &
                                      solver_cntl,still_sequential)
!          --------------------------

type(mumps_matrix_type), intent(inout) :: mumps_matrix
type(linsys_type), intent(in) :: phaml_matrix
type(proc_info), intent(in) :: procs
type(solver_options), intent(in) :: solver_cntl
logical, intent(in) :: still_sequential

call fatal("Dummy version of create_mumps_linear_system called",procs=procs)

end subroutine create_mumps_linear_system

!          ---------------------------
subroutine destroy_mumps_linear_system(mumps_matrix,procs)
!          ---------------------------

type(mumps_matrix_type), intent(inout) :: mumps_matrix
type(proc_info), intent(in) :: procs

call fatal("Dummy version of destroy_mumps_linear_system called",procs=procs)

end subroutine destroy_mumps_linear_system

!          ----------------
subroutine change_mumps_rhs(mumps_matrix,phaml_matrix,rs,procs, &
                            still_sequential)
!          ----------------

type(mumps_matrix_type), intent(inout) :: mumps_matrix
type(linsys_type), intent(in) :: phaml_matrix
real(my_real), intent(in) :: rs(:)
type(proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential

call fatal("Dummy version of change_mumps_rhs called",procs=procs)

end subroutine change_mumps_rhs

!          -----------
subroutine mumps_solve(mumps_matrix,phaml_matrix,procs,still_sequential,&
                       noshadow)
!          -----------

type(mumps_matrix_type), intent(inout) :: mumps_matrix
type(linsys_type), intent(inout) :: phaml_matrix
type(proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential
logical, intent(in), optional :: noshadow

call fatal("Dummy version of mumps_solve called",procs=procs)

end subroutine mumps_solve

end module mumps_interf
