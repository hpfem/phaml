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
! This version illustrates several scalar linear elliptic PDEs.
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

real(my_real) :: e,a,t,tx,ty,txx,tyy,px,qy,rp117,rp5,r1,r2,r117,r500,r1000, &
                 r2000,r,dcxydx
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
! Begin executable code

cxy=0; cx=0; cy=0

select case (choice)

   case (1,2,3,4,5) ! Laplace problems
      cxx = 1.0_my_real
      cyy = 1.0_my_real
      c = 0.0_my_real
      rs = 0.0_my_real

   case (6,8) ! x**ipower+y**ipower
      cxx = 1.0_my_real
      cyy = 1.0_my_real
      c = 0.0_my_real
      rs = -ipower*(ipower-1)*(x**(ipower-2)+y**(ipower-2))

   case (7) ! peak
      cxx = 1.0_my_real
      cyy = 1.0_my_real
      c = 0.0_my_real
      rp117 = 0.117_my_real
      rp5 = 0.5_my_real
      r1 = 1.0_my_real 
      r2 = 2.0_my_real 
      r117 = 117.0_my_real
      r500 = 500.0_my_real
      r1000 = 1000.0_my_real
      r2000 = 2000.0_my_real
      e = exp(-r500*((x-rp5)**2+(y-rp117)**2))
      rs = -e*(-r2*y*(r1-y)+r2*(r1-x)*y*(r1-y)*(-r1000*x+r500)               &
           -r2*x*y*(r1-y)*(-r1000*x+r500)-r2000*x*(r1-x)*y*(r1-y)            &
           +x*(r1-x)*y*(r1-y)*(-r1000*x+r500)**2 -r2*x*(r1-x)                &
           +r2*x*(r1-x)*(r1-y)*(-r1000*y+r117)-r2*x*(r1-x)*y*(-r1000*y+r117) &
           +x*(r1-x)*y*(r1-y)*(-r1000*y+r117)**2)

   case (9) ! nonconstant coefficients
      cxx = 1 + x*x
      a = 4.0_my_real*y*y + 0.9_my_real
      cyy = 1 + a*a
      c = -(1.0_my_real+(8.0_my_real*y-x-4.0_my_real)**2)
      t = trues(x,y,1,1)
      tx = cos(x)*sin(y)
      ty = sin(x)*cos(y)
      txx = -sin(x)*sin(y)
      tyy = -sin(x)*sin(y)
      px = 2.0_my_real*x
      qy = 2.0_my_real*a*8.0_my_real*y
      rs = -cxx(1,1)*txx-px*tx -cyy(1,1)*tyy-qy*ty  +c(1,1)*t

   case (10) ! x(1-x)y(1-y)
      cxx = 1.0_my_real
      cyy = 1.0_my_real
      c = 0.0_my_real
      rs = 2.0_my_real*(x*(1.0_my_real-x)+y*(1.0_my_real-y))

   case (11) ! Helmholtz
      cxx = 1.0_my_real
      cyy = 1.0_my_real
      c = -prob_param**2
      rs = 0.0_my_real

   case (12) ! Poisson with same solution as case 11
      cxx = 1.0_my_real
      cyy = 1.0_my_real
      c = 0.0_my_real
      rs = prob_param**2 * trues(x,y,1,1)

   case (13) ! sin(prob_param*r)
      cxx = 1.0_my_real
      cyy = 1.0_my_real
      c = 0.0_my_real
      r = prob_param*sqrt(x**2+y**2)
      rs = prob_param**2*(sin(r) -cos(r)/r)

   case (14) ! x**ipower+y**ipower with Helmholtz
      cxx = 1.0_my_real
      cyy = 1.0_my_real
      c = prob_param
      rs = -ipower*(ipower-1)*(x**(ipower-2)+y**(ipower-2)) + &
            prob_param * (x**ipower + y**ipower)

   case (15) ! tanh wave front
      cxx = 1.0_my_real
      cyy = 1.0_my_real
      c = 0.0_my_real
      r = prob_param*(x**2 + y**2 - .5)
      rs = -4*prob_param/cosh(r)**2 + &
            8*prob_param**2*(x**2+y**2)*tanh(r)/cosh(r)**2

   case (16) ! mock eigenvalue problem
      cxx = 1.0_my_real
      cyy = 1.0_my_real
      c = 0.5_my_real*prob_param*(x*x+3.0_my_real*y*y)
      rs = sqrt(prob_param*(2.0_my_real+sqrt(3.0_my_real)))*trues(x,y,1,1)

   case (17) ! first order terms
      cxx = 1.0_my_real
      cyy = 1.0_my_real
      cx = 1  + x*(1-x) + y*(1-y)
      cy = 1  + x*(1-x) + y*(1-y)
      c = 0.0_my_real
      rs = -ipower*(ipower-1)*(x**(ipower-2)+y**(ipower-2)) &
           + cx(1,1)*truexs(x,y,1,1) + cy(1,1)*trueys(x,y,1,1)

   case (18) ! mixed derivative term
      cxx = 1.0_my_real
      cyy = 1.0_my_real
      cxy = 1+x*(1-y)
      dcxydx = (1-y)
      c = 0.0_my_real
      rs = -ipower*(ipower-1)*(x**(ipower-2)*y**ipower+x**ipower*y**(ipower-2))&
           - cxy(1,1)*ipower**2*x**(ipower-1)*y**(ipower-1) &
           - dcxydx*ipower*x**ipower*y**(ipower-1)

   case default
      print *,"ERROR -- pdecoefs called with unknown problem id"
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
! bmark is assumed to give the following values for the unit square:
!  1 - lower left corner
!  2 - left side
!  3 - upper left corner
!  4 - top
!  5 - upper right corner
!  6 - right side
!  7 - lower right corner
!  8 - bottom
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

end interface

!----------------------------------------------------
! Local variables

!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

select case (choice)

   case (1,2,3,6,7,9,10,11,12,13,14,15,16,17,18) ! DIRICHLET known solutions, use trues
      itype = DIRICHLET
      c = 0.0_my_real
      rs = trues(x,y,1,1)

   case (4) ! unknown solution
      itype = DIRICHLET
      c = 0.0_my_real
      select case (bmark)
         case (3,4,5) ! top
            rs = x*x*(1.0_my_real-x)*(1.0_my_real-x)
         case default
            rs = 0.0_my_real
      end select

   case (5) ! unknown solution
      itype = DIRICHLET
      c = 0.0_my_real
      select case (bmark)
         case (3,4,5) ! top
            if (x < prob_param) then
               rs = x/prob_param
            else
               rs = (1.0_my_real-x)/(1.0_my_real-prob_param)
            endif
         case default
            rs = 0.0_my_real
      end select

   case (8) ! natural and mixed conditions
      select case (bmark)
         case (2,3) ! left side
            itype = NATURAL
            c = 0.0_my_real
            if (ipower == 1) then
               rs = -1.0_my_real
            else
               rs = 0.0_my_real
            endif
         case (4) ! top
            itype = MIXED
            c = x
            rs = ipower*y**(ipower-1) + c(1,1)*trues(x,y,1,1)
         case (5,6) ! right side
            itype = NATURAL
            c = 0.0_my_real
            rs = ipower*x**(ipower-1)
         case (7,8,1) ! bottom
            itype = DIRICHLET
            c = 0.0_my_real
            rs = trues(x,y,1,1)
         case default
            print *,"ERROR -- bconds called with invalid bmark ",bmark
            stop
      end select

   case default
      print *,"ERROR -- bconds called with unknown problem id"
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
! Local variables

real (my_real) :: p,q,pi,r,theta1
!----------------------------------------------------
! Begin executable code

select case (choice)

   case (1) ! x+y
      trues = x+y

   case (2) ! something more complicated but still
            ! satisfies Laplace's equation
      p=x**4-6.0_my_real*x**2*y**2+y**4
      q=x**8-28.0_my_real*x**6*y**2+70.0_my_real*x**4*y**4-28.0_my_real*x**2*y**6+y**8
      trues = 1.1786_my_real-0.1801_my_real*p+0.006_my_real*q

   case (3) ! singularity at the origin, from the L domain problem
      pi=4.0_my_real*atan(1.0_my_real)
      r = sqrt(x*x+y*y)
      if (x.eq.0.0_my_real) then
         if(y.gt.0.0_my_real) then
            theta1=pi/2.0_my_real
         else
            theta1=3.0_my_real*pi/2.0_my_real
         endif
      else
         if (x.gt.0.0_my_real) then
            if(y.ge.00_my_real) then
              theta1=atan(y/x)
            else
              theta1=atan(y/x)+2.0_my_real*pi
            endif
         else 
            theta1=atan(y/x)+pi
         endif
      endif
      trues = r**(2.0_my_real/3.0_my_real)*sin(2.0_my_real*theta1/3.0_my_real)

   case (6,8,14,17) ! x**ipower+y**ipower
      trues = x**ipower + y**ipower

   case (7) ! peak
      trues = x*(1.0_my_real-x)*y*(1.0_my_real-y)* &
              exp(-500.0_my_real*((x-0.5_my_real)**2+(y-0.117_my_real)**2))

   case (4,5) ! unknown
      trues = huge(0.0_my_real)

   case (9) ! simpler than R&B #54
      trues = sin(x)*sin(y)

   case (10) ! x(1-x)y(1-y)
      trues = x*(1.0_my_real-x)*y*(1.0_my_real-y)

   case (11,12)
      trues = sin(prob_param*x) + sin(prob_param*y)

   case (13)
      trues = sin(prob_param*sqrt(x**2+y**2))

   case (15) ! tanh wave front
      trues = tanh(prob_param*(x**2+y**2-.5))

   case (16) ! mock eigenvalue problem, first eigenfunction
      trues = exp(-sqrt(prob_param/8)*(x*x+sqrt(3.0_my_real)*y*y))

   case (18) ! x**ipower*y**ipower
      trues = x**ipower * y**ipower

   case default
      print *,"ERROR -- trues called with unknown problem id"
      stop

end select

end function trues

!        ------
function truexs(x,y,comp,eigen) ! real (my_real)
!        ------

!----------------------------------------------------
! This is the x derivative of the true solution of the differential equation,
! if known.
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
! Local variables

real (my_real) :: p,q,pi,r,theta1
!----------------------------------------------------
! Begin executable code

select case (choice)

   case (6,8,14,17) ! x**ipower+y**ipower

      truexs = ipower*x**(ipower-1)

   case (9) ! simpler than R&B #54
      truexs = cos(x)*sin(y)

   case (10) ! x(1-x)y(1-y)
      truexs = (1.0_my_real-2.0_my_real*x)*y*(1.0_my_real-y)

   case (11,12) ! sin(a*x) + sin(a*y)
      truexs = prob_param*cos(prob_param*x)

   case (15) ! tanh wave front
      truexs = 2*prob_param*x/cosh(prob_param*(x**2+y**2-.5))**2

   case (16) ! mock eigenvalue problem, first eigenfunction
      truexs = exp(-sqrt(prob_param/8)*(x*x+sqrt(3.0_my_real)*y*y))
      truexs = -2*sqrt(prob_param/8)*x*truexs

   case (18) ! x**ipower*y**ipower
      truexs = ipower*x**(ipower-1) * y**ipower

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
! Local variables

real (my_real) :: p,q,pi,r,theta1
!----------------------------------------------------
! Begin executable code

select case (choice)

   case (6,8,14,17) ! x**ipower+y**ipower

      trueys = ipower*y**(ipower-1)

   case (9) ! simpler than R&B #54
      trueys = sin(x)*cos(y)

   case (10) ! x(1-x)y(1-y)
      trueys = (1.0_my_real-2.0_my_real*y)*x*(1.0_my_real-x)

   case (11,12) ! sin(a*x) + sin(a*y)
      trueys = prob_param*cos(prob_param*y)

   case (15) ! tanh wave front
      trueys = 2*prob_param*y/cosh(prob_param*(x**2+y**2-.5))**2

   case (16) ! mock eigenvalue problem, first eigenfunction
      trueys = exp(-sqrt(prob_param/8)*(x*x+sqrt(3.0_my_real)*y*y))
      trueys = -2*sqrt(3.0_my_real)*sqrt(prob_param/8)*y*trueys

   case (18) ! x**ipower*y**ipower
      trueys = ipower*y**(ipower-1) * x**ipower

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

   iparam(1) = ipower
   rparam(1) = prob_param

! Call the routine that performs the actual exchange.  Don't change this line.

   call master_to_slaves(phaml_solution,iparam,rparam)

! Copy the arrays into the module variables, using the same correspondence
! between module variable and array index as was used above.

   ipower = iparam(1)
   prob_param = rparam(1)

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
