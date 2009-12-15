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
! This file contains the user supplied external subroutines that define
! the PDE(s) to be solved, and other required external subroutines.
!   pdecoefs bconds boundary_point boundary_npiece boundary_param iconds trues
!   truexs trueys update_usermod phaml_integral_kernel
!
! This version illustrates the subroutines for a single PDE, thus all the
! array arguments have dimension 1.  For a coupled system of PDEs, see
! examples/system.
!----------------------------------------------------

module phaml_user_mod
integer :: probno=0
end module phaml_user_mod

!          --------
subroutine pdecoefs(x,y,cxx,cxy,cyy,cx,cy,c,rs)
!          --------

!----------------------------------------------------
! This subroutine returns the coefficient and right hand side of the PDE
! at the point (x,y)
!
! The PDE is
!
!    -( cxx(x,y) * u  )  -( cyy(x,y) * u  ) + c(x,y) * u = rs(x,y)
!                   x  x                y  y
!
! For eigenvalue problems, the right hand side is lambda * u * rs(x,y)
!
! cxy, cx and cy are not currently used and should not be set.  They are
! for future expansion.
!
! NOTE: BE CAREFUL TO GET THE SIGNS RIGHT
! e.g. cxx=cyy=1 means rs=-(uxx+uyy)
!
!----------------------------------------------------

!----------------------------------------------------
! Modules used are:

use phaml
use phaml_user_mod
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x,y
real(my_real), intent(out) :: cxx(:,:),cxy(:,:),cyy(:,:),cx(:,:),cy(:,:), &
                              c(:,:),rs(:)
!----------------------------------------------------

!----------------------------------------------------
! Local variables

real(my_real) :: pi,a,t,tx,ty,txx,tyy,px,qy
!----------------------------------------------------

interface

   function trues(x,y,comp,eigen) ! real (my_real)
   use phaml
   real (my_real), intent(in) :: x,y
   integer, intent(in) :: comp,eigen
   real (my_real) :: trues
   end function trues
   function truexs(x,y,comp,eigen) ! real (my_real)
   use phaml
   real (my_real), intent(in) :: x,y
   integer, intent(in) :: comp,eigen
   real (my_real) :: truexs
   end function truexs
   function trueys(x,y,comp,eigen) ! real (my_real)
   use phaml
   real (my_real), intent(in) :: x,y
   integer, intent(in) :: comp,eigen
   real (my_real) :: trueys
   end function trueys

end interface

!----------------------------------------------------
! Begin executable code

cx=0; cy=0; cxy=0

select case(probno)
case(1)
   cxx(1,1) = 1.0_my_real
   cyy(1,1) = 1.0_my_real
   c(1,1) = 0.0_my_real
   pi = 4*atan(1.0_my_real)
   rs(1) = (4*pi**2*x**2 - 2)*sin(2*pi*y + 1)
case(2)
   cxx(1,1) = 1.0_my_real
   cyy(1,1) = 1.0_my_real
   c(1,1) = 0.0_my_real
   rs(1) = -4.0_my_real
case(3)
   cxx = 1 + x*x
   a = 4.0_my_real*y*y + 0.9_my_real
   cyy = 1 + a*a
   c = -(1.0_my_real+(8.0_my_real*y-x-4.0_my_real)**2)
   t = trues(x,y,1,1)
   tx = truexs(x,y,1,1)
   ty = trueys(x,y,1,1)
   txx = -sin(x)*sin(y)
   tyy = -sin(x)*sin(y)
   px = 2.0_my_real*x
   qy = 2.0_my_real*a*8.0_my_real*y
   rs = -cxx(1,1)*txx-px*tx -cyy(1,1)*tyy-qy*ty  +c(1,1)*t
case(4)
   cxx = 1.0_my_real
   cyy = 1.0_my_real
   cx = 1  + x*(1-x) + y*(1-y)
   cy = 1  + x*(1-x) + y*(1-y)
   c = 0.0_my_real
   rs = -90*(x**8+y**8) &
        + cx(1,1)*truexs(x,y,1,1) + cy(1,1)*trueys(x,y,1,1)
case(5)
   cxx = 1.0_my_real
   cyy = 1.0_my_real
   cxy = 1.0_my_real + x*(1-y)
   c = 0.0_my_real
   rs = -90*(x**8*y**10+x**10*y**8) &
           - cxy(1,1)*100*x**9*y**9 &
           - (1-y)*10*x**10*y**9

case default
   print *,"illegal problem number in pdecoefs"
   stop
end select

end subroutine pdecoefs

!          ------
subroutine bconds(x,y,bmark,itype,c,rs)
!          ------

!----------------------------------------------------
! This subroutine returns the boundary conditions at the point (x,y).
!
! Each boundary condition is either
!
!    u = rs(x,y) or  u  + c(x,y)*u = rs(x,y)
!                     n
!
! itype indicates which type of condition applies to each point, using
! symbolic constants from module phaml.  It is DIRICHLET for the first
! condition, NATURAL for the second condition with c==0, and MIXED for
! the second condition with c/=0.
!
!----------------------------------------------------

!----------------------------------------------------
! Modules used are:

use phaml
use phaml_user_mod
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x,y
integer, intent(in) :: bmark
integer, intent(out) :: itype(:)
real(my_real), intent(out) :: c(:,:),rs(:)
!----------------------------------------------------

!----------------------------------------------------
! Non-module procedures used are:

interface

   function trues(x,y,comp,eigen) ! real (my_real)
   use phaml
   real (my_real), intent(in) :: x,y
   integer, intent(in) :: comp,eigen
   real (my_real) :: trues
   end function trues
   function truexs(x,y,comp,eigen) ! real (my_real)
   use phaml
   real (my_real), intent(in) :: x,y
   integer, intent(in) :: comp,eigen
   real (my_real) :: truexs
   end function truexs
   function trueys(x,y,comp,eigen) ! real (my_real)
   use phaml
   real (my_real), intent(in) :: x,y
   integer, intent(in) :: comp,eigen
   real (my_real) :: trueys
   end function trueys

end interface

!----------------------------------------------------
! Local variables

!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

! Dirichlet boundary conditions

select case(probno)
case(1)
   select case(bmark)
   case (1)
      itype = DIRICHLET
      c = 0.0_my_real
      rs = trues(x,y,1,1)
   case (2)
      itype = NATURAL
      c = 0.0_my_real
      rs = 0.0_my_real
   case (-8, 8)
      itype = PERIODIC
      c = 0.0_my_real
      rs = 0.0_my_real
   case default
      print *,"illegal bmark passed to bconds"
   end select

case(2)
   select case(bmark)
   case(1,3,5,7,8)
      itype = DIRICHLET
      c = 0.0_my_real
      rs(1) = trues(x,y,1,1)
   case(2)
      itype = NATURAL
      c = 0.0_my_real
      rs(1) = -truexs(x,y,1,1)
   case(4)
      itype = MIXED
      c = 2.0_my_real
      rs(1) = trueys(x,y,1,1) + 2.0_my_real*trues(x,y,1,1)
   case(6)
      itype = NATURAL
      c = 0.0_my_real
      rs(1) = truexs(x,y,1,1)
   case default
      print *,"illegal bmark in bconds"
   end select

case (3,4,5)
   itype = DIRICHLET
   c = 0.0_my_real
   rs(1) = trues(x,y,1,1)

case default
   print *,"illegal problem number in bconds"
   stop
end select

end subroutine bconds

!        ------
function iconds(x,y,comp,eigen)
!        ------

!----------------------------------------------------
! This routine returns the initial condition for a time dependent problem.
! It can also be used for the initial guess for a nonlinear problem, or
! systems of equations solved by sucessive substitution.
! comp,eigen is which solution to use from a coupled system of PDEs or multiple
! eigenvectors from an eigenvalue problem, and is ignored in this example.
! For problems where there are no initial conditions, it is a dummy.
!----------------------------------------------------

!----------------------------------------------------
! Modules used are:

use phaml
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x,y
integer, intent(in) :: comp,eigen
real(my_real) :: iconds

!----------------------------------------------------
! Begin executable code

iconds = 0.0_my_real

end function iconds

!        -----
function trues(x,y,comp,eigen) ! real (my_real)
!        -----

!----------------------------------------------------
! This is the true solution of the differential equation, if known.
! comp,eigen is which solution to use from a coupled system of PDEs or multiple
! eigenvectors from an eigenvalue problem, and is ignored in this example.
!----------------------------------------------------

!----------------------------------------------------
! Modules used are:

use phaml
use phaml_user_mod
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x,y
integer, intent(in) :: comp,eigen
real (my_real) :: trues
!----------------------------------------------------

real(my_real) :: pi
!----------------------------------------------------
! Begin executable code

select case(probno)
case(1)
   pi = 4*atan(1.0_my_real)
   trues = x**2 * sin(2*pi*y + 1)
case(2)
   trues = x**2 + y**2
case(3)
   trues = sin(x)*sin(y)
case(4)
   trues = x**10 + y**10
case(5)
   trues = x**10 * y**10
case default
   print *,"illegal problem number in trues"
   stop
end select

end function trues

!        ------
function truexs(x,y,comp,eigen) ! real (my_real)
!        ------

!----------------------------------------------------
! This is the x derivative of the true solution of the differential
! equation, if known.
! comp,eigen is which solution to use from a coupled system of PDEs or multiple
! eigenvectors from an eigenvalue problem, and is ignored in this example.
!----------------------------------------------------

!----------------------------------------------------
! Modules used are:

use phaml
use phaml_user_mod
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x,y
integer, intent(in) :: comp,eigen
real (my_real) :: truexs
!----------------------------------------------------

real (my_real) :: pi
!----------------------------------------------------
! Begin executable code

select case(probno)
case(1)
   pi = 4*atan(1.0_my_real)
   truexs = 2*x * sin(2*pi*y + 1)
case(2)
   truexs = 2*x
case(3)
   truexs = cos(x)*sin(y)
case(4)
   truexs = 10*x**9
case(5)
   truexs = 10*x**9 * y**10
case default
   print *,"illegal problem number in truexs"
   stop
end select

end function truexs

!        ------
function trueys(x,y,comp,eigen) ! real (my_real)
!        ------

!----------------------------------------------------
! This is the y derivative of the true solution of the differential
! equation, if known.
! comp,eigen is which solution to use from a coupled system of PDEs or multiple
! eigenvectors from an eigenvalue problem, and is ignored in this example.
!----------------------------------------------------

!----------------------------------------------------
! Modules used are:

use phaml
use phaml_user_mod
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x,y
integer, intent(in) :: comp,eigen
real (my_real) :: trueys
!----------------------------------------------------

real(my_real) :: pi
!----------------------------------------------------
! Begin executable code

select case (probno)
case(1)
   pi = 4*atan(1.0_my_real)
   trueys = x**2 * 2*pi*cos(2*pi*y + 1)
case(2)
   trueys = 2*y
case(3)
   trueys = sin(x)*cos(y)
case(4)
   trueys = 10*y**9
case(5)
   trueys = 10*y**9 * x**10
case default
   print *,"illegal problem number in trueys"
   stop
end select

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
use phaml_user_mod
type(phaml_solution_type), intent(in) :: phaml_solution

integer :: iparam(1)
real(my_real) :: rparam(1)

iparam(1) = probno
rparam(1) = 0.0_my_real
call master_to_slaves(phaml_solution,iparam,rparam)
probno = iparam(1)

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
