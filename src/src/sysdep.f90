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
!----------------------------------------------------

!----------------------------------------------------

implicit none
private
public my_system, my_flush

!----------------------------------------------------

contains

!          ---------
subroutine my_system(cmd)
!          ---------

!----------------------------------------------------
! This routine executes a system command as if from the system command line.
! It is not supported by the standard, but I haven't had it fail on any
! compilers I have tried.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

character(len=*), intent(in) :: cmd
!----------------------------------------------------
! Begin executable code

call system(cmd)

end subroutine my_system

!          --------
subroutine my_flush(unit)
!          --------

!----------------------------------------------------
! This routine flushes the io buffer of the given unit.
!
! The flush statement is standard Fortran 2003.  According to
! ACM Fortran Forum August 2007, it is supported by Cray, gfortran,
! g95, IBM, and Intel, but not by NAG.  The article did not include
! Lahey, Absoft, PGI, Pathscale, or other compilers.
!
! Some compilers supply flush as a nonstandard intrinsic subroutine, so
! if the flush statement fails, you might try "call flush(unit)".
!
! If the compiler does not support flush, then just return without doing
! anything.
!
! Documentation for Lahey says it has flush as a subroutine.
!
! Documentation for Nag says it has a flush subroutine in a nonstandard
! intrinsic module, so might need
! use f90_unix_io, only:flush
!
! Documentation for Absoft makes no reference to flush.
!
! Documentation for PGI makes no reference to flush (except as an OpenMP
! directive), but the Pathscale documentation implies that PGI has it
! as a subroutine.
!
! Documentation for Pathscale says it has flush as a subroutine, but you might
! need -g77 or -pgi on the compiler command line.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: unit
!----------------------------------------------------
! Begin executable code

! Since more compilers support the nonstandard flush subroutine than the
! standard (2003) flush statement, I'll use the subroutine for now.

call flush(unit)

end subroutine my_flush

end module sysdep
