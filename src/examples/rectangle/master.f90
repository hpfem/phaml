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

!       ------------
program phaml_master
!       ------------

!----------------------------------------------------
! This is the main program.
!----------------------------------------------------

!----------------------------------------------------
! Modules used are:

use phaml
!----------------------------------------------------

implicit none

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
xmax = 1.0
ymin = 0.0
ymax = 2.0
ngridx = 4
ngridy = 8
call make_data_files(ngridx,ngridy,xmin,xmax,ymin,ymax)

call phaml_create(soln,nproc=4,triangle_files="rectangle.1", &
draw_grid_who=MASTER)

call phaml_solve_pde(soln,                         &
                     max_vert=3200,                &
                     draw_grid_when=PHASES,        &
                     pause_at_start=.true.,        &
                     pause_after_phases=.true.,    &
                     print_grid_when=PHASES,       &
                     print_grid_who=MASTER,        &
                     print_error_when=PHASES,      &
                     print_error_what=ENERGY_LINF_ERR, &
                     print_errest_what=ENERGY_LINF_ERREST, &
                     print_error_who=MASTER)

call phaml_destroy(soln)

end program phaml_master

!-----------------------------------------------------------------------

subroutine make_data_files(nx,ny,xmin,xmax,ymin,ymax)

! This subroutine creates a node file for an nx X ny grid on the rectangular
! domain (xmin,xmax) X (ymin,ymax), and runs triangle to create the
! required .node, .ele, .neigh and .edge files.
!
! bmark is set to be 1, 3, 5 and 7 at the lower left, upper left, upper
! right and lower right vertices, respectively, and 2, 4, 6 and 8 at the
! left, top, right and bottom edges, respectively.
!
! Running triangle uses the common extension "call system", which is
! compiler dependent.  That call may need to be changed, or may not be
! supported at all, on some compilers.

implicit none

integer, intent(in) :: nx, ny
real, intent(in) :: xmin,xmax,ymin,ymax

integer :: i, j, count, bmark
real :: xvals(0:nx), yvals(0:ny)

! values for the grid lines

do i=0,nx
   xvals(i) = xmin + (xmax-xmin)*i/real(nx)
end do

do i=0,ny
   yvals(i) = ymin + (ymax-ymin)*i/real(ny)
end do

! write a .poly file with the grid of nodes, boundary edges, and bmark

open(unit=21,file="rectangle.poly",status="replace")

! nodes

write(21,*) (nx+1)*(ny+1), 2, 0, 1
count = 0
do i=0,nx
   do j=0,ny
      if (i==0) then
         if (j==0) then
            bmark = 1
         elseif (j==ny) then
            bmark = 3
         else
            bmark = 2
         endif
      elseif (i==nx) then
         if (j==0) then
            bmark = 7
         elseif (j==ny) then
            bmark = 5
         else
            bmark = 6
         endif
      else
         if (j==0) then
            bmark = 8
         elseif (j==ny) then
            bmark = 4
         else
            bmark = 0
         endif
      endif
      count = count + 1
      write(21,*) count, xvals(i), yvals(j), bmark
   end do
end do

! edges

count = 0
write(21,*) 2*nx+2*ny,1
do i=1,ny
   count = count + 1
   write(21,*) count,i,i+1,2
end do
do i=1,nx
   count = count + 1
   write(21,*) count,(i-1)*(ny+1)+1,i*(ny+1)+1,8
   count = count + 1
   write(21,*) count,i*(ny+1),(i+1)*(ny+1),4
end do
do i=1,ny
   count = count + 1
   write(21,*) count,nx*(ny+1)+i,nx*(ny+1)+i+1,6
end do

write(21,*) 0
close(unit=21)

! run triangle
! NOTE: THIS IS COMPILER DEPENDENT.  You may need to change this statement.

call system("triangle -neqj rectangle.poly")

end subroutine make_data_files
