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

module hypretype_mod

!----------------------------------------------------
! This module contains a data structure to hold the hypre form of a matrix.
!----------------------------------------------------

use global

implicit none
private
public hypre_matrix_type, hypre_options, hypre_pointer

!----------------------------------------------------

integer, parameter :: dp = kind(0.0d0)

! kind of integer to hold a hypre pointer.  Documentation for hypre claims
! it should be integer*8 on most systems.

integer, parameter :: hypre_pointer = selected_int_kind(18)

! The hypre matrix data type

type hypre_matrix_type
   integer(hypre_pointer) :: ij_matrix      ! type HYPRE_IJMatrix
   integer(hypre_pointer) :: parcsr_matrix  ! type HYPRE_ParCSRMatrix
   integer(hypre_pointer) :: ij_rhs         ! type HYPRE_IJVector
   integer(hypre_pointer) :: par_rhs        ! type HYPRE_ParVector
   integer(hypre_pointer) :: ij_solution    ! type HYPRE_IJVector
   integer(hypre_pointer) :: par_solution   ! type HYPER_ParVector
   integer :: neq
   integer, pointer :: firsteq(:), global_eq(:)
   real(my_real), pointer :: lapack_rhs(:)
end type hypre_matrix_type

! Options to pass to hypre

type hypre_options
   integer :: BoomerAMG_MaxLevels
   integer :: BoomerAMG_MaxIter
   real(dp):: BoomerAMG_Tol
   real(dp):: BoomerAMG_StrongThreshold
   real(dp):: BoomerAMG_MaxRowSum
   integer :: BoomerAMG_CoarsenType
   integer :: BoomerAMG_MeasureType
   integer :: BoomerAMG_CycleType
   integer, pointer :: BoomerAMG_NumGridSweeps(:)
   integer, pointer :: BoomerAMG_GridRelaxType(:)
   integer, pointer :: BoomerAMG_GridRelaxPoints(:,:)
   real(dp), pointer:: BoomerAMG_RelaxWeight(:)
   integer :: BoomerAMG_DebugFlag
   real(dp):: ParaSails_thresh
   integer :: ParaSails_nlevels
   real(dp):: ParaSails_filter
   integer :: ParaSails_sym
   real(dp):: ParaSails_loadbal
   integer :: ParaSails_reuse
   integer :: ParaSails_logging
   real(dp):: PCG_Tol
   integer :: PCG_MaxIter
   integer :: PCG_TwoNorm
   integer :: PCG_RelChange
   integer :: PCG_Logging
   integer :: GMRES_KDim
   real(dp):: GMRES_Tol
   integer :: GMRES_MaxIter
   integer :: GMRES_Logging
end type hypre_options

end module hypretype_mod
