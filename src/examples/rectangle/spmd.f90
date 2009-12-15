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

!----------------------------------------------------
! This file contains the main program supplied by the user.
! This is a simple example that just solves one linear elliptic pde.
!----------------------------------------------------

!       ---------------
program phaml_spawnless
!       ---------------

!----------------------------------------------------
! Modules used are:

use phaml
use mpif_mod
implicit none

!----------------------------------------------------
! Local variables:

integer :: whodrawg
integer :: jerr
integer :: my_proc, total_nproc
integer :: nslave, subtract, divide

!----------------------------------------------------
! Begin executable code

! set the graphics options

whodrawg = MASTER

! initialize MPI, find out how many processors and what my rank is

call mpi_init(jerr)
call mpi_comm_size(MPI_COMM_WORLD,total_nproc,jerr)
call mpi_comm_rank(MPI_COMM_WORLD,my_proc,jerr)

! determine how many processors for slaves and graphics

subtract = 1
if (whodrawg == MASTER .or. whodrawg == EVERYONE) subtract = 2
divide = 1
if (whodrawg == SLAVES .or. whodrawg == EVERYONE) divide = 2
nslave = (total_nproc-subtract)/divide

! call the master, slave or graphics program depending on my rank

if (my_proc == 0) then
   call phaml_master(whodrawg,nslave)
elseif (my_proc <= nslave) then
   call phaml_slave
else
   call phaml_graphics
endif

stop
end program phaml_spawnless

!          ------------
subroutine phaml_master(whodrawg,nslave)
!          ------------

!----------------------------------------------------
! This is the main program.
!----------------------------------------------------

!----------------------------------------------------
! Modules used are:

use phaml
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

integer :: whodrawg, nslave
!----------------------------------------------------
! Local variables

type(phaml_solution_type) :: soln
integer :: ngridx, ngridy
real :: xmin,xmax,ymin,ymax

interface
   subroutine make_data_files(nx,ny,xmin,xmax,ymin,ymax)
   integer, intent(in) :: nx, ny
   real, intent(in) :: xmin,xmax,ymin,ymax
   end subroutine make_data_files
end interface

!----------------------------------------------------
! Begin executable code

! create the triangle data files

xmin = 0.0
xmax = 0.5
ymin = 0.0
ymax = 1.0
ngridx = 4
ngridy = 8
call make_data_files(ngridx,ngridy,xmin,xmax,ymin,ymax)

call phaml_create(soln,nproc=nslave,draw_grid_who=whodrawg, &
                     triangle_files="rectangle.1")

call phaml_solve_pde(soln,                   &
                     max_vert=1000,         &
                     draw_grid_when=PHASES,  &
                     pause_after_phases=.false.,   &
                     print_grid_when=PHASES, &
                     print_grid_who=MASTER,  &
                     print_error_when=PHASES,&
                     print_error_what=LINF_ERR, &
                     print_error_who=MASTER)

call phaml_destroy(soln)

end subroutine phaml_master

!-----------------------------------------------------------------------

subroutine make_data_files(nx,ny,xmin,xmax,ymin,ymax)

! This subroutine creates a node file for an nx X ny grid on the rectangular
! domain (xmin,xmax) X (ymin,ymax), and runs triangle to create the
! required .node, .ele, .neigh and .edge files.
!
! Running triangle uses the common extension "call system", which is
! compiler dependent.  That call may need to be changed, or may not be
! supported at all, on some compilers.

implicit none

integer, intent(in) :: nx, ny
real, intent(in) :: xmin,xmax,ymin,ymax

integer :: i, j, count
real :: xvals(0:nx), yvals(0:ny)

! values for the grid lines

do i=0,nx
   xvals(i) = xmin + (xmax-xmin)*i/real(nx)
end do

do i=0,ny
   yvals(i) = ymin + (ymax-ymin)*i/real(ny)
end do

! write a .node file with the grid of vertices

open(unit=21,file="rectangle.node",status="replace")
write(21,*) (nx+1)*(ny+1), 2, 0, 0

count = 0
do i=0,nx
   do j=0,ny
      count = count + 1
      write(21,*) count, xvals(i), yvals(j)
   end do
end do

close(unit=21)

! run triangle
! NOTE: THIS IS COMPILER DEPENDENT.  You may need to change this statement.

call system("triangle -n -e rectangle.node")

end subroutine make_data_files
