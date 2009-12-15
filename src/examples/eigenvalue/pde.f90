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
! This version illustrates eigenvalue problems via the Schroedinger Equation
! for a simple harmonic oscillator.
!----------------------------------------------------

!          --------
subroutine pdecoefs(x,y,cxx,cxy,cyy,cx,cy,c,rs)
!          --------

!----------------------------------------------------
! This subroutine returns the coefficient and right hand side of the PDE
! at the point (x,y)
!
! The PDE is
!
!    -( cxx(x,y) * u  )  -( cyy(x,y) * u  ) + c(x,y) * u = rs(x,y) * lambda * u
!                   x  x                y  y
! with
!   cxx = cyy = 1.0
!   c = prob_param * (x^2+3y^2)/2
!   f = 1.0
!
! prob_param determines the radius of the well
!
! cxy, cx and cy are not currently used and should not be set.  They are
! for future expansion.
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

!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

cxx(1,1) = 1.0_my_real
cyy(1,1) = 1.0_my_real
c(1,1) = 0.5_my_real*prob_param*(x*x+3.0_my_real*y*y)
rs(1) = 1.0_my_real

cxy=0; cx=0; cy=0

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
! Local variables

!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

! Dirichlet boundary conditions

itype = DIRICHLET
c = 0.0_my_real
rs(1) = 0.0_my_real

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
! eigenvectors from an eigenvalue problem.
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

real (my_real) :: xmax,ymax,tmax,d,sq3
!----------------------------------------------------
! Begin executable code

! This solution is only valid for the smallest few eigenvalues.  We'll assume
! they are in order, but that might not be true if lambda0 /= -infinity.

sq3 = sqrt(3.0_my_real)
d = sqrt(prob_param/8)

select case (eigen)

case(1)

   trues = exp(-d*(x*x+sq3*y*y))

! Eigenvectors are unique up to a multiplicative constant, and are scaled to
! have M-norm 1.0; this is approximately the scaling factor for this solution

!   trues = 1.217006012276_my_real * trues

! TEMP disabling the 2nd and 3rd true solutions because 1) the scale factor
!      seems to be slightly off, 2) can't determine the sign of the scale factor

case(2)

   trues = huge(0.0_my_real)
!   xmax = sqrt(1/(2*d))
!   tmax = xmax*exp(-d*xmax**2)
!   trues = x*exp(-d*(x*x+sq3*y*y))
!   trues = 1.04009347278_my_real * trues / tmax

case(3)

   trues = huge(0.0_my_real)
!   ymax = sqrt(1/(2*sq3*d))
!   tmax = ymax*exp(-d*sq3*ymax**2)
!   trues = y*exp(-d*(x*x+sq3*y*y))
!   trues = 1.038366656409_my_real * trues / tmax

case default

   trues = huge(0.0_my_real)

end select

end function trues

!        ------
function truexs(x,y,comp,eigen) ! real (my_real)
!        ------

!----------------------------------------------------
! This is the x derivative of the true solution of the differential equation,
! if known.
! comp,eigen is which solution to use from a coupled system of PDEs or multiple
! eigenvectors from an eigenvalue problem.
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

real (my_real) :: xmax,ymax,tmax,d,sq3
!----------------------------------------------------
! Begin executable code

sq3 = sqrt(3.0_my_real)
d = sqrt(prob_param/8)

select case (eigen)

case(1)

   truexs = exp(-d*(x*x+sq3*y*y))
   truexs = -2*d*x*truexs
!   truexs = 1.217006012276_my_real * truexs

! TEMP disabling the 2nd and 3rd true solutions because 1) the scale factor
!      seems to be slightly off, 2) can't determine the sign of the scale factor

case(2)

   truexs = huge(0.0_my_real)
!   xmax = sqrt(1/(2*d))
!   tmax = xmax*exp(-d*xmax**2)
!   truexs = exp(-d*(x*x+sq3*y*y))
!   truexs = (1-2*d*x*x)*truexs
!   truexs = 1.04009347278_my_real * truexs / tmax

case(3)

   truexs = huge(0.0_my_real)
!   ymax = sqrt(1/(2*sq3*d))
!   tmax = ymax*exp(-d*sq3*ymax**2)
!   truexs = exp(-d*(x*x+sq3*y*y))
!   truexs = -2*d*x*y*truexs
!   truexs = 1.038366656409_my_real * truexs / tmax

case default

   truexs = huge(0.0_my_real)

end select

end function truexs

!        ------
function trueys(x,y,comp,eigen) ! real (my_real)
!        ------

!----------------------------------------------------
! This is the y derivative of the true solution of the differential equation,
! if known.
! comp,eigen is which solution to use from a coupled system of PDEs or multiple
! eigenvectors from an eigenvalue problem.
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

real (my_real) :: xmax,ymax,tmax,d,sq3
!----------------------------------------------------
! Begin executable code

sq3 = sqrt(3.0_my_real)
d = sqrt(prob_param/8)

select case (eigen)

case(1)

   trueys = exp(-d*(x*x+sq3*y*y))
   trueys = -2*sq3*d*y*trueys
!   trueys = 1.217006012276_my_real * trueys

! TEMP disabling the 2nd and 3rd true solutions because 1) the scale factor
!      seems to be slightly off, 2) can't determine the sign of the scale factor

case(2)

   trueys = huge(0.0_my_real)
!   xmax = sqrt(1/(2*d))
!   tmax = xmax*exp(-d*xmax**2)
!   trueys = exp(-d*(x*x+sq3*y*y))
!   trueys = -2*sq3*d*y*x*trueys
!   trueys = 1.04009347278_my_real * trueys / tmax

case(3)

   trueys = huge(0.0_my_real)
!   ymax = sqrt(1/(2*sq3*d))
!   tmax = ymax*exp(-d*sq3*ymax**2)
!   trueys = exp(-d*(x*x+sq3*y*y))
!   trueys = (1-2*d*sq3*y**2)*trueys
!   trueys = 1.038366656409_my_real * trueys / tmax

case default

   trueys = huge(0.0_my_real)

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

!----------------------------------------------------
! This routine updates the module variables on the slave processes by
! sending them from the master process
!----------------------------------------------------

use phaml
use phaml_user_mod

!----------------------------------------------------
! Dummy arguments


type(phaml_solution_type), intent(in) :: phaml_solution

!----------------------------------------------------
! Local variables:

! Declare these arrays big enough to hold the variables to be sent

integer :: iparam(1)
real(my_real) :: rparam(1)

!----------------------------------------------------
! Begin executable code 

! Copy the module variables into the arrays, putting integer variables
! into iparam and real variables into rparam.

   rparam(1) = prob_param

! Call the routine that performs the actual exchange.  Don't change this line.

   call master_to_slaves(phaml_solution,iparam,rparam)

! Copy the arrays into the module variables, using the same correspondence
! between module variable and array index as was used above.

   prob_param  = rparam(1)

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
