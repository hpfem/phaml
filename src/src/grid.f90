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

module grid_mod

!----------------------------------------------------
! This module contains subroutines to manipulate the grid.
!
! communication tags in this module are of the form 4xx
!----------------------------------------------------

!----------------------------------------------------
! Other modules used:

use global
!----------------------------------------------------

implicit none
private
!----------------------------------------------------
! Non-module procedures used are:

!----------------------------------------------------
! The following variables are defined:

real(my_real) :: n3pee(3,1), npq(3)

!----------------------------------------------------
! The following defined constants are defined:

!----------------------------------------------------

!contains

end module grid_mod
