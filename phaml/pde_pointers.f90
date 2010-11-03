module pde_pointers

use phaml
implicit none

private
public set_pde_pointers, pdecoefs_ptr, bconds_ptr, iconds_ptr, trues_ptr, &
        truexs_ptr, trueys_ptr

interface

   subroutine pdecoefs_proc(x,y,cxx,cxy,cyy,cx,cy,c,rs)
   double precision, intent(in) :: x,y
   double precision, intent(out) :: cxx(:,:),cxy(:,:),cyy(:,:),cx(:,:), &
                                    cy(:,:),c(:,:),rs(:)
   end subroutine pdecoefs_proc

   subroutine bconds_proc(x,y,bmark,itype,c,rs)
   double precision, intent(in) :: x,y
   integer, intent(in) :: bmark
   integer, intent(out) :: itype(:)
   double precision, intent(out) :: c(:,:),rs(:)
   end subroutine bconds_proc

   function iconds_proc(x,y,comp,eigen)
   double precision, intent(in) :: x,y
   integer, intent(in) :: comp,eigen
   double precision :: iconds_proc
   end function iconds_proc

   function trues_proc(x,y,comp,eigen)
   double precision, intent(in) :: x,y
   integer, intent(in) :: comp,eigen
   double precision :: trues_proc
   end function trues_proc

   function truexs_proc(x,y,comp,eigen)
   double precision, intent(in) :: x,y
   integer, intent(in) :: comp,eigen
   double precision :: truexs_proc
   end function truexs_proc

   function trueys_proc(x,y,comp,eigen)
   double precision, intent(in) :: x,y
   integer, intent(in) :: comp,eigen
   double precision :: trueys_proc
   end function trueys_proc

end interface


procedure(pdecoefs_proc), pointer, save :: pdecoefs_ptr
procedure(bconds_proc), pointer, save :: bconds_ptr
procedure(iconds_proc), pointer, save :: iconds_ptr
procedure(trues_proc), pointer, save :: trues_ptr
procedure(truexs_proc), pointer, save :: truexs_ptr
procedure(trueys_proc), pointer, save :: trueys_ptr

contains

subroutine set_pde_pointers(pdecoefs, bconds, iconds, trues, truexs, trueys)
procedure(pdecoefs_proc) :: pdecoefs
procedure(bconds_proc) :: bconds
procedure(iconds_proc) :: iconds
procedure(trues_proc) :: trues
procedure(truexs_proc) :: truexs
procedure(trueys_proc) :: trueys
pdecoefs_ptr => pdecoefs
bconds_ptr => bconds
iconds_ptr => iconds
trues_ptr => trues
truexs_ptr => truexs
trueys_ptr => trueys
end subroutine

end module pde_pointers

subroutine pdecoefs(x,y,cxx,cxy,cyy,cx,cy,c,rs)
use phaml
use pde_pointers
real(my_real), intent(in) :: x,y
real(my_real), intent(out) :: cxx(:,:),cxy(:,:),cyy(:,:),cx(:,:),cy(:,:), &
                              c(:,:),rs(:)
call pdecoefs_ptr(x,y,cxx,cxy,cyy,cx,cy,c,rs)
end subroutine pdecoefs

subroutine bconds(x,y,bmark,itype,c,rs)
use phaml
use pde_pointers
real(my_real), intent(in) :: x,y
integer, intent(in) :: bmark
integer, intent(out) :: itype(:)
real(my_real), intent(out) :: c(:,:),rs(:)
call bconds_ptr(x,y,bmark,itype,c,rs)
end subroutine bconds

function iconds(x,y,comp,eigen)
use phaml
use pde_pointers
real(my_real), intent(in) :: x,y
integer, intent(in) :: comp,eigen
real(my_real) :: iconds
iconds = iconds_ptr(x,y,comp,eigen)
end function iconds

function trues(x,y,comp,eigen) ! real (my_real)
use phaml
use pde_pointers
real(my_real), intent(in) :: x,y
integer, intent(in) :: comp,eigen
real (my_real) :: trues
trues = trues_ptr(x,y,comp,eigen)
end function trues

function truexs(x,y,comp,eigen) ! real (my_real)
use phaml
use pde_pointers
real(my_real), intent(in) :: x,y
integer, intent(in) :: comp,eigen
real (my_real) :: truexs
truexs = truexs_ptr(x,y,comp,eigen)
end function truexs

function trueys(x,y,comp,eigen) ! real (my_real)
use phaml
use pde_pointers
real(my_real), intent(in) :: x,y
integer, intent(in) :: comp,eigen
real (my_real) :: trueys
trueys = trueys_ptr(x,y,comp,eigen)
end function trueys

!          --------------
subroutine boundary_point(ipiece,s,x,y)
!          --------------

!----------------------------------------------------
! This routine defines the boundary of the domain.  It returns the point
! (x,y) with parameter s on piece ipiece.
! If boundary_npiece <= 0 it is a dummy and the domain is given by triangle
! data files.
!----------------------------------------------------

!----------------------------------------------------
! Modules used are:

use phaml
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: ipiece
real(my_real), intent(in) :: s
real(my_real), intent(out) :: x,y

!----------------------------------------------------
! Begin executable code

! Dummy version

x = 0.0_my_real
y = 0.0_my_real

end subroutine boundary_point

!        ---------------
function boundary_npiece(hole)
!        ---------------

!----------------------------------------------------
! This routine gives the number of pieces in the boundary definition.
! If boundary_npiece <= 0 the domain is given by triangle data files.
!----------------------------------------------------

!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: hole
integer :: boundary_npiece
!----------------------------------------------------
! Begin executable code

boundary_npiece = 0

end function boundary_npiece

!          --------------
subroutine boundary_param(start,finish)
!          --------------

!----------------------------------------------------
! This routine gives the range of parameter values for each piece of the
! boundary.
! If boundary_npiece <= 0 it is a dummy and the domain is given by triangle
! data files.
!----------------------------------------------------

!----------------------------------------------------
! Modules used are:

use phaml
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(out) :: start(:), finish(:)
!----------------------------------------------------
! Begin executable code

! Dummy version

start = 0.0_my_real
finish = 0.0_my_real

end subroutine boundary_param

!          --------------
subroutine update_usermod(phaml_solution)
!          --------------

use phaml
type(phaml_solution_type), intent(in) :: phaml_solution

! Dummy version.  See examples/elliptic/usermod.f90 for a functional example.

end subroutine update_usermod

!        ---------------------
function phaml_integral_kernel(kernel,x,y)
!        ---------------------

use phaml
integer, intent(in) :: kernel
real(my_real), intent(in) :: x,y
real(my_real) :: phaml_integral_kernel

! Identity function

phaml_integral_kernel = 1.0

end function phaml_integral_kernel

!        ----------
function regularity(x,y)
!        ----------

use phaml
real(my_real), intent(in) :: x(3),y(3)
real(my_real) :: regularity

! Dummy version, assume infinitely differentiable everywhere.

regularity = huge(0.0_my_real)

end function regularity

