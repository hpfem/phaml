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

! This is a fixed form file that wraps dmumps_struc.h into a module.  This
! provides access to all the entities in dmumps_struc.h and eliminates the
! problem of INCLUDEing a fixed form file into free form code.

      module mumps_struc_mod
      include "dmumps_struc.h"

! also define the mumps linear system data structure here

      type mumps_matrix_type
         type(dmumps_struc) :: mumps_var
         integer, pointer :: global_eq(:), firsteq(:)
      end type mumps_matrix_type

      end module mumps_struc_mod
