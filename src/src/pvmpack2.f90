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

! This contains four external routines that act as an interface to
! pvmfpack, pvmfunpack, pvmfpsend and pvmfprecv when the argument is type 
! real(my_real).
! This extra routine, in a separate file, is required to shut up picky
! fortran 90 compilers that won't allow sending different types of arguments
! to the same routine.

subroutine pvmfpack_my_real(type,xp,nitem,stride,info)
use global, only: my_real
integer :: type,nitem,stride,info
real(my_real) :: xp(*)
interface
!NAS$ ALIEN "F77 pvmfpack"
    subroutine pvmfpack(type,xp,nitem,stride,info)
    use global, only: my_real
    integer :: type,nitem,stride,info
    real(my_real) :: xp(*)
    end subroutine pvmfpack
end interface
call pvmfpack(type,xp,nitem,stride,info)
end subroutine pvmfpack_my_real
 
subroutine pvmfunpack_my_real(type,xp,nitem,stride,info)
use global, only: my_real
integer :: type,nitem,stride,info
real(my_real) :: xp(*)
interface
!NAS$ ALIEN "F77 pvmfunpack"
    subroutine pvmfunpack(type,xp,nitem,stride,info)
    use global, only: my_real
    integer :: type,nitem,stride,info
    real(my_real) :: xp(*)
    end subroutine pvmfunpack
end interface
call pvmfunpack(type,xp,nitem,stride,info)
end subroutine pvmfunpack_my_real

subroutine pvmfpsend_my_real(tid,msgtag,buf,len,datatype,info)
use global, only: my_real
integer :: tid,msgtag,len,datatype,info
real(my_real) :: buf(*)
interface
!NAS$ ALIEN "F77 pvmfpsend"
   subroutine pvmfpsend(tid,msgtag,buf,len,datatype,info)
   use global, only: my_real
   integer :: tid,msgtag,len,datatype,info
   real(my_real) :: buf(*)
   end subroutine pvmfpsend
end interface
call pvmfpsend(tid,msgtag,buf,len,datatype,info)
end subroutine pvmfpsend_my_real

subroutine pvmfprecv_my_real(tid,msgtag,buf,len,datatype,atid,atag,alen,info)
use global, only: my_real
integer :: tid,msgtag,len,datatype,atid,atag,alen,info
real(my_real) :: buf(*)
interface
!NAS$ ALIEN "F77 pvmfprecv"
   subroutine pvmfprecv(tid,msgtag,buf,len,datatype,atid,atag,alen,info)
   use global, only: my_real
   integer :: tid,msgtag,len,datatype,atid,atag,alen,info
   real(my_real) :: buf(*)
   end subroutine pvmfprecv
end interface
call pvmfprecv(tid,msgtag,buf,len,datatype,atid,atag,alen,info)
end subroutine pvmfprecv_my_real
