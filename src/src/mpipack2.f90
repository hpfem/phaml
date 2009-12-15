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

! This contains external routines that act as an interface to
! mpi routines that accept multiple types of arguments when the argument is type
! real(my_real).
! This extra routine, in a separate file, is required to shut up picky
! fortran 90 compilers that won't allow sending different types of arguments
! to the same routine.

subroutine mpi_pack_my_real(inbuf,incount,datatype,outbuf,outsize,position, &
                            comm,ierror)
use global, only: my_real
real(my_real) :: inbuf(*)
integer :: outbuf(*)
integer :: incount,datatype,outsize,position,comm,ierror
interface
!NAS$ ALIEN "C mpi_pack"
   subroutine mpi_pack(inbuf,incount,datatype,outbuf,outsize,position, &
                       comm,ierror)
   use global, only: my_real
   real(my_real) :: inbuf(*)
   integer :: outbuf(*)
   integer :: incount,datatype,outsize,position,comm,ierror
   end subroutine mpi_pack
end interface
call mpi_pack(inbuf,incount,datatype,outbuf,outsize,position,comm,ierror)
end subroutine mpi_pack_my_real

subroutine mpi_unpack_my_real(inbuf,insize,position,outbuf,outcount, &
                              datatype,comm,ierror)
use global, only: my_real
integer :: inbuf(*)
real(my_real) :: outbuf(*)
integer :: insize,position,outcount,datatype,comm,ierror
interface
!NAS$ ALIEN "C mpi_unpack"
   subroutine mpi_unpack(inbuf,insize,position,outbuf,outcount, &
                         datatype,comm,ierror)
   use global, only: my_real
   integer :: inbuf(*)
   real(my_real) :: outbuf(*)
   integer :: insize,position,outcount,datatype,comm,ierror
   end subroutine mpi_unpack
end interface
call mpi_unpack(inbuf,insize,position,outbuf,outcount,datatype,comm,ierror)
end subroutine mpi_unpack_my_real

subroutine mpi_reduce_my_real(sendbuf,recvbuf,count,datatype,op,root, &
                              comm,ierror)
use global, only: my_real
real(my_real) :: sendbuf(*)
real(my_real) :: recvbuf(*)
integer :: count,datatype,op,root,comm,ierror
interface
!NAS$ ALIEN "C mpi_reduce"
   subroutine mpi_reduce(sendbuf,recvbuf,count,datatype,op,root, &
                         comm,ierror)
   use global, only: my_real
   real(my_real) :: sendbuf(*)
   real(my_real) :: recvbuf(*)
   integer :: count,datatype,op,root,comm,ierror
   end subroutine mpi_reduce
end interface
call mpi_reduce(sendbuf,recvbuf,count,datatype,op,root, &
                comm,ierror)
end subroutine mpi_reduce_my_real

subroutine mpi_allreduce_my_real(sendbuf,recvbuf,count,datatype,op, &
                              comm,ierror)
use global, only: my_real
real(my_real) :: sendbuf(*)
real(my_real) :: recvbuf(*)
integer :: count,datatype,op,comm,ierror
interface
!NAS$ ALIEN "C mpi_allreduce"
   subroutine mpi_allreduce(sendbuf,recvbuf,count,datatype,op, &
                         comm,ierror)
   use global, only: my_real
   real(my_real) :: sendbuf(*)
   real(my_real) :: recvbuf(*)
   integer :: count,datatype,op,comm,ierror
   end subroutine mpi_allreduce
end interface
call mpi_allreduce(sendbuf,recvbuf,count,datatype,op, &
                   comm,ierror)
end subroutine mpi_allreduce_my_real

subroutine mpi_alltoallv_my_real(sendbuf, sendcounts, sdispls, sendtype, & 
                                recvbuf, recvcounts, rdispls, recvtype, comm, &
                                ierror)
use global, only: my_real
real(my_real) :: sendbuf(*), recvbuf(*)
integer :: sendcounts(*), sdispls(*)
integer :: recvcounts(*), rdispls(*)
integer :: sendtype, recvtype, comm, ierror
interface
!NAS$ ALIEN "C mpi_alltoallv"
   subroutine mpi_alltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf, &
                            recvcounts, rdispls, recvtype, comm, ierror)
   use global, only: my_real
   real(my_real) :: sendbuf(*), recvbuf(*)
   integer :: sendcounts(*), sdispls(*)
   integer :: recvcounts(*), rdispls(*)
   integer :: sendtype, recvtype, comm, ierror
   end subroutine mpi_alltoallv
end interface
call mpi_alltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf, &
                   recvcounts, rdispls, recvtype, comm, ierror)
end subroutine mpi_alltoallv_my_real

subroutine mpi_allgatherv_real(sendbuf, sendcount, sendtype, recvbuf, &
                               recvcounts, rdispls, recvtype, comm, ierror)
use global, only: my_real
real(my_real) :: sendbuf(*), recvbuf(*)
integer :: sendcount
integer :: recvcounts(*), rdispls(*)
integer :: sendtype, recvtype, comm, ierror
interface
!NAS$ ALIEN "C mpi_allgatherv"
   subroutine mpi_allgatherv(sendbuf, sendcount, sendtype, recvbuf, & 
                             recvcounts, rdispls, recvtype, comm, ierror)
   use global, only: my_real
   real(my_real) :: sendbuf(*), recvbuf(*)
   integer :: sendcount
   integer :: recvcounts(*), rdispls(*)
   integer :: sendtype, recvtype, comm, ierror
   end subroutine mpi_allgatherv
end interface
call mpi_allgatherv(sendbuf, sendcount, sendtype, recvbuf, &
                    recvcounts, rdispls, recvtype, comm, ierror)
end subroutine mpi_allgatherv_real
