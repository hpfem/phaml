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

! MPI 2 routines with string arguments; pass as integer array for portability

! Unfortunately, we cannot use a C wrapper for safe character string passing
! with mpi_comm_spawn and mpi_info_set because of the handle arguments (info
! and communicators).  In Fortran the handle arguments are integers, but
! in C they are structures.  So the C side of these wrappers cannot call
! the C library functions with the Fortran handle arguments.  And they
! cannot call the Fortran library functions because the name mangling is
! processor dependent.  Therefore we must assume that the Fortran compiler
! used for compiling PHAML agrees with the compilers used for building
! the MPI library in how it passes character arguments, and hope for the best.
!
! But we still need an alternate form of mpi_comm_spawn because LAM (at 6.5.4
! at least) uses double complex for MPI_ARGV_NULL instead of character.  So
! use mpi_comm_spawn when passing MPI_ARGV_NULL and mpi_comm_spawn_char when
! passing a character argv.  Also note that MPI_ERRCODES_IGNORE (which PHAML
! always uses) is also double complex in LAM instead of an array of integers,
! so errcodes is incorrectly declared.

subroutine mpi_comm_spawn_char(command,argv,maxprocs,info,root,comm, &
                               intercomm,errcodes,ierror)
character(len=*) :: command, argv(*)
integer :: maxprocs, info, root, comm, intercomm, ierror
double complex errcodes

!interface
!!NAS$ ALIEN "C mpi_comm_spawn"
!   subroutine mpi_comm_spawn(command,argv,maxprocs,info,root,comm,intercomm, &
!                             errcodes,ierror)
!   character(len=*) :: command, argv(*)
!   integer :: maxprocs, info, root, comm, intercomm, errcodes(*), ierror
!   end subroutine mpi_comm_spawn
!end interface

call mpi_comm_spawn(command,argv,maxprocs,info,root,comm,intercomm,errcodes, &
                    ierror)
end subroutine mpi_comm_spawn_char
