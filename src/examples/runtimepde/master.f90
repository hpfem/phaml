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

!       ------------
program phaml_master
!       ------------

!----------------------------------------------------
! This is the main program.
!----------------------------------------------------

!----------------------------------------------------
! Modules used are:

use phaml
use pde_intf
use message_passing, only : PARALLEL
use global, only : SEQUENTIAL
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Local variables

!----------------------------------------------------
! Interface blocks for the external subroutines that define the PDEs

interface

subroutine pdecoef1(x,y,cxx,cxy,cyy,cx,cy,c,rs)
implicit none
double precision, intent(in) :: x,y
double precision, intent(out) :: cxx(:,:),cxy(:,:),cyy(:,:),cx(:,:),cy(:,:), &
                                 c(:,:),rs(:)
end subroutine pdecoef1

subroutine pdecoef2(x,y,cxx,cxy,cyy,cx,cy,c,rs)
implicit none
double precision, intent(in) :: x,y
double precision, intent(out) :: cxx(:,:),cxy(:,:),cyy(:,:),cx(:,:),cy(:,:), &
                                 c(:,:),rs(:)
end subroutine pdecoef2

subroutine dirich_bcond(x,y,bmark,itype,c,rs)
implicit none
double precision, intent(in) :: x,y
integer, intent(in) :: bmark
integer, intent(out) :: itype(:)
double precision, intent(out) :: c(:,:),rs(:)
end subroutine dirich_bcond

function true1(x,y,comp,eigen)
implicit none
double precision, intent(in) :: x,y
integer, intent(in) :: comp,eigen
double precision :: true1
end function true1

function true2(x,y,comp,eigen)
implicit none
double precision, intent(in) :: x,y
integer, intent(in) :: comp,eigen
double precision :: true2
end function true2

function truex1(x,y,comp,eigen)
implicit none
double precision, intent(in) :: x,y
integer, intent(in) :: comp,eigen
double precision :: truex1
end function truex1

function truex2(x,y,comp,eigen)
implicit none
double precision, intent(in) :: x,y
integer, intent(in) :: comp,eigen
double precision :: truex2
end function truex2

function truey1(x,y,comp,eigen)
implicit none
double precision, intent(in) :: x,y
integer, intent(in) :: comp,eigen
double precision :: truey1
end function truey1

end interface

!----------------------------------------------------
! Begin executable code

! Verify that the compilation of PHAML is sequential and double precision

if (PARALLEL /= SEQUENTIAL) then
   print *,'This example cannot be run in parallel.'
   print *,'Recompile PHAML with "PARALLEL sequential"'
   stop
endif

if (my_real /= kind(1.0d0)) then
   print *,'This example requires double precision.'
   print *,'Set my_real to kind(1.0d0) in global.f90'
   stop
endif

! Initialize the pde interface module and soln

call init_pde_intf
call phaml_create(1)

! Set the first pde

call set_pdecoefs(pdecoef1)
call set_bconds(dirich_bcond)
call set_trues(true1)
call set_truexs(truex1)
call set_trueys(truey1)

! Solve the first pde

print *,"solving the first pde"

call phaml_solve_pde(1,                      &
                     max_vert=1000,          &
                     pause_at_end=.true., &
                     print_grid_when=PHASES, &
                     print_grid_who=MASTER,  &
                     print_error_when=PHASES,&
                     print_error_what=ENERGY_LINF_ERR, &
                     print_errest_what=ENERGY_LINF_ERREST, &
                     print_error_who=MASTER)

! change the pde and true routines; bconds can stay the same because it
! uses true to set the Dirichlet boundary condition.  I also don't have to
! change truey because it just happens to be the same.

call set_pdecoefs(pdecoef2)
call set_trues(true2)
call set_truexs(truex2)

! solve the second pde, using the solution of the first pde as an initial guess

print *,"solving second pde"

call phaml_solve_pde(1,                      &
                     max_vert=40000,          &
                     print_grid_when=PHASES, &
                     print_grid_who=MASTER,  &
                     print_error_when=PHASES,&
                     print_error_what=ENERGY_LINF_ERR, &
                     print_errest_what=ENERGY_LINF_ERREST, &
                     print_error_who=MASTER)

call phaml_destroy(1)

end program phaml_master

! External routines to define the PDEs.  All other routines use the defaults.

subroutine pdecoef1(x,y,cxx,cxy,cyy,cx,cy,c,rs)
implicit none
double precision, intent(in) :: x,y
double precision, intent(out) :: cxx(:,:),cxy(:,:),cyy(:,:),cx(:,:),cy(:,:), &
                                 c(:,:),rs(:)
cxx = 1; cyy = 1
cxy = 0; cx = 0; cy = 0; c = 0
rs = -20*(x**3 + y**3)
end subroutine pdecoef1

subroutine pdecoef2(x,y,cxx,cxy,cyy,cx,cy,c,rs)
implicit none
double precision, intent(in) :: x,y
double precision, intent(out) :: cxx(:,:),cxy(:,:),cyy(:,:),cx(:,:),cy(:,:), &
                                 c(:,:),rs(:)
cxx = 1; cyy = 1
cxy = 0; cx = 0; cy = 0; c = 0
rs = -20*((1-x)**3 + y**3)
end subroutine pdecoef2

subroutine dirich_bcond(x,y,bmark,itype,c,rs)
use phaml
implicit none
double precision, intent(in) :: x,y
integer, intent(in) :: bmark
integer, intent(out) :: itype(:)
double precision, intent(out) :: c(:,:),rs(:)
interface
   function trues(x,y,comp,eigen)
   double precision, intent(in) :: x,y
   integer, intent(in) :: comp,eigen
   double precision :: trues
   end function trues
end interface
itype = DIRICHLET
c = 0
rs = trues(x,y,1,1)
end subroutine dirich_bcond

function true1(x,y,comp,eigen)
implicit none
double precision, intent(in) :: x,y
integer, intent(in) :: comp,eigen
double precision :: true1
true1 = x**5 + y**5
end function true1

function true2(x,y,comp,eigen)
implicit none
double precision, intent(in) :: x,y
integer, intent(in) :: comp,eigen
double precision :: true2
true2 = (1-x)**5 + y**5
end function true2

function truex1(x,y,comp,eigen)
implicit none
double precision, intent(in) :: x,y
integer, intent(in) :: comp,eigen
double precision :: truex1
truex1 = 5*x**4
end function truex1

function truex2(x,y,comp,eigen)
implicit none
double precision, intent(in) :: x,y
integer, intent(in) :: comp,eigen
double precision :: truex2
truex2 = -5*(1-x)**4
end function truex2

function truey1(x,y,comp,eigen)
implicit none
double precision, intent(in) :: x,y
integer, intent(in) :: comp,eigen
double precision :: truey1
truey1 = 5*y**4
end function truey1
