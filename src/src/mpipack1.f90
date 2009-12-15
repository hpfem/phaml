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
! integer.
! This extra routine, in a separate file, is required to shut up picky
! fortran 90 compilers that won't allow sending different types of arguments
! to the same routine.

subroutine mpi_pack_int(inbuf,incount,datatype,outbuf,outsize,position, &
                        comm,ierror)
integer :: inbuf(*)
integer :: outbuf(*)
integer :: incount,datatype,outsize,position,comm,ierror
interface
!NAS$ ALIEN "C mpi_pack"
   subroutine mpi_pack(inbuf,incount,datatype,outbuf,outsize,position, &
                       comm,ierror)
   integer :: inbuf(*)
   integer :: outbuf(*)
   integer :: incount,datatype,outsize,position,comm,ierror
   end subroutine mpi_pack
end interface
call mpi_pack(inbuf,incount,datatype,outbuf,outsize,position,comm,ierror)
end subroutine mpi_pack_int

subroutine mpi_unpack_int(inbuf,insize,position,outbuf,outcount, &
                          datatype,comm,ierror)
integer :: inbuf(*)
integer :: outbuf(*)
integer :: insize,position,outcount,datatype,comm,ierror
interface
!NAS$ ALIEN "C mpi_unpack"
   subroutine mpi_unpack(inbuf,insize,position,outbuf,outcount, &
                         datatype,comm,ierror)
   integer :: inbuf(*)
   integer :: outbuf(*)
   integer :: insize,position,outcount,datatype,comm,ierror
   end subroutine mpi_unpack
end interface
call mpi_unpack(inbuf,insize,position,outbuf,outcount,datatype,comm,ierror)
end subroutine mpi_unpack_int

subroutine mpi_reduce_int(sendbuf,recvbuf,count,datatype,op,root, &
                             comm,ierror)
integer :: sendbuf(*)
integer :: recvbuf(*)
integer :: count,datatype,op,root,comm,ierror
interface
!NAS$ ALIEN "C mpi_reduce"
   subroutine mpi_reduce(sendbuf,recvbuf,count,datatype,op,root, &
                         comm,ierror)
   integer :: sendbuf(*)
   integer :: recvbuf(*)
   integer :: count,datatype,op,root,comm,ierror
   end subroutine mpi_reduce
end interface
call mpi_reduce(sendbuf,recvbuf,count,datatype,op,root, &
                comm,ierror)
end subroutine mpi_reduce_int

subroutine mpi_allreduce_int(sendbuf,recvbuf,count,datatype,op, &
                                comm,ierror)
integer :: sendbuf(*)
integer :: recvbuf(*)
integer :: count,datatype,op,comm,ierror
interface
!NAS$ ALIEN "C mpi_allreduce"
   subroutine mpi_allreduce(sendbuf,recvbuf,count,datatype,op, &
                         comm,ierror)
   integer :: sendbuf(*)
   integer :: recvbuf(*)
   integer :: count,datatype,op,comm,ierror
   end subroutine mpi_allreduce
end interface
call mpi_allreduce(sendbuf,recvbuf,count,datatype,op, &
                   comm,ierror)
end subroutine mpi_allreduce_int

subroutine mpi_alltoallv_int(sendbuf, sendcounts, sdispls, sendtype, &
                                recvbuf, recvcounts, rdispls, recvtype, comm, &
                                ierror)
integer :: sendbuf(*), sendcounts(*), sdispls(*)
integer :: recvbuf(*), recvcounts(*), rdispls(*)
integer :: sendtype, recvtype, comm, ierror
interface
!NAS$ ALIEN "C mpi_alltoallv"
   subroutine mpi_alltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf, &
                            recvcounts, rdispls, recvtype, comm, ierror)
   integer :: sendbuf(*), sendcounts(*), sdispls(*)
   integer :: recvbuf(*), recvcounts(*), rdispls(*)
   integer :: sendtype, recvtype, comm, ierror
   end subroutine mpi_alltoallv
end interface
call mpi_alltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf, &
                   recvcounts, rdispls, recvtype, comm, ierror)
end subroutine mpi_alltoallv_int

subroutine mpi_allgatherv_int(sendbuf, sendcount, sendtype, recvbuf, &
                              recvcounts, rdispls, recvtype, comm, ierror)
integer :: sendbuf(*), sendcount
integer :: recvbuf(*), recvcounts(*), rdispls(*)
integer :: sendtype, recvtype, comm, ierror
interface
!NAS$ ALIEN "C mpi_allgatherv"
   subroutine mpi_allgatherv(sendbuf, sendcount, sendtype, recvbuf, &
                             recvcounts, rdispls, recvtype, comm, ierror)
   integer :: sendbuf(*), sendcount
   integer :: recvbuf(*), recvcounts(*), rdispls(*)
   integer :: sendtype, recvtype, comm, ierror
   end subroutine mpi_allgatherv
end interface
call mpi_allgatherv(sendbuf, sendcount, sendtype, recvbuf, &
                    recvcounts, rdispls, recvtype, comm, ierror)
end subroutine mpi_allgatherv_int
