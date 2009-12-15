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

module sysdep

!----------------------------------------------------
! This module contains routines that are not standard Fortran 90 and may
! need to be changed for different compilers.
!
! This version is known to work with the following compilers:
!  Lahey LF95 on Linux
!----------------------------------------------------

!----------------------------------------------------

implicit none
private
public my_system

!----------------------------------------------------

contains

!          ---------
subroutine my_system(cmd)
!          ---------

!----------------------------------------------------
! This routine executes a system command as if from the system command line
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

character(len=*), intent(in) :: cmd
!----------------------------------------------------
! Begin executable code

print *,"WARNING -- subroutine system was called, but does not exist in this compilation"

end subroutine my_system

end module sysdep
