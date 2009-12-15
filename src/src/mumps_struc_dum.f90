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

! This is a dummy version of mumps_struc_mod for use when MUMPS is not used.
! It contains the definition of the dmumps_struc type with all the components
! that PHAML uses, and a dummy subroutine dmumps.  This is just to satisfy
! references in linsys.f90, and should never actually be used.

module mumps_struc_mod
use global

type dmumps_struc
   integer :: job, comm, sym, par, n, nz_loc
   integer :: icntl(40), info(40)
   integer, pointer :: irn_loc(:), jcn_loc(:)
   real(my_real), pointer :: a_loc(:), rhs(:)
end type dmumps_struc

type mumps_matrix_type
   integer :: dum
end type mumps_matrix_type

end module mumps_struc_mod
