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

! MPI routines with string arguments; pass as integer array for portability

subroutine mpi_my_get_processor_name(name, l, ierr)
character(len=*) name
integer l, ierr

integer, allocatable :: int_name(:)
integer i, nchar

interface
!NAS$ ALIEN "C mpi_myc_get_processor_name"
   subroutine mpi_myc_get_processor_name(name,al,l,ierr)
   integer name(*),al,l,ierr
   end subroutine mpi_myc_get_processor_name
end interface

allocate(int_name(len(name)))
nchar = len(name)
call mpi_myc_get_processor_name(int_name, nchar, l, ierr)
name = " "
do i=1,nchar
   name(i:i) = char(int_name(i))
end do

deallocate(int_name)
end subroutine mpi_my_get_processor_name
