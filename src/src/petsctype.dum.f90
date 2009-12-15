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

module petsc_type_mod

!----------------------------------------------------
! This module contains a dummy version of the data structure to hold the
! PETSc form of a matrix.
!----------------------------------------------------

implicit none
private
public petsc_matrix_type, petsc_options, petsc_dummy

!----------------------------------------------------

type petsc_matrix_type
   integer :: nada
end type petsc_matrix_type

type petsc_options
   real(kind(0.0d0)) :: petsc_richardson_damping_factor
   real(kind(0.0d0)) :: petsc_chebychev_emin
   real(kind(0.0d0)) :: petsc_chebychev_emax
   real(kind(0.0d0)) :: petsc_rtol
   real(kind(0.0d0)) :: petsc_atol
   real(kind(0.0d0)) :: petsc_dtol
   real(kind(0.0d0)) :: petsc_ilu_dt
   real(kind(0.0d0)) :: petsc_ilu_dtcol
   real(kind(0.0d0)) :: petsc_sor_omega
   real(kind(0.0d0)) :: petsc_eisenstat_omega
   integer :: petsc_gmres_max_steps
   integer :: petsc_maxits
   integer :: petsc_ilu_levels
   integer :: petsc_icc_levels
   integer :: petsc_ilu_maxrowcount
   integer :: petsc_sor_its
   integer :: petsc_sor_lits
   integer :: petsc_asm_overlap
   logical :: petsc_eisenstat_nodiagscaling
end type petsc_options

type(petsc_options), parameter :: petsc_dummy = petsc_options(0.0d0, 0.0d0, &
   0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0, 0, 0, 0, 0, 0, &
   0, 0, .false.)

end module petsc_type_mod
