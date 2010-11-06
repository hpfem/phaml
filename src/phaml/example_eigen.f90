module example_eigen

use phaml
implicit none

private
public setup_example_eigen

real(my_real), save :: prob_param=25.0_my_real

contains

subroutine setup_example_eigen()
use pde_pointers
call set_pde_pointers(pdecoefs, bconds, iconds, trues, truexs, trueys)
end subroutine

subroutine pdecoefs(x,y,cxx,cxy,cyy,cx,cy,c,rs)
real(my_real), intent(in) :: x,y
real(my_real), intent(out) :: cxx(:,:),cxy(:,:),cyy(:,:),cx(:,:),cy(:,:), &
                              c(:,:),rs(:)
cxx(1,1) = 1.0_my_real
cyy(1,1) = 1.0_my_real
c(1,1) = 0.5_my_real*prob_param*(x*x+3.0_my_real*y*y)
rs(1) = 1.0_my_real

cxy=0; cx=0; cy=0
end subroutine pdecoefs

subroutine bconds(x,y,bmark,itype,c,rs)
real(my_real), intent(in) :: x,y
integer, intent(in) :: bmark
integer, intent(out) :: itype(:)
real(my_real), intent(out) :: c(:,:),rs(:)
itype = DIRICHLET
c = 0.0_my_real
rs(1) = 0.0_my_real
end subroutine bconds

function iconds(x,y,comp,eigen)
real(my_real), intent(in) :: x,y
integer, intent(in) :: comp,eigen
real(my_real) :: iconds
iconds = 0.0_my_real
end function iconds

function trues(x,y,comp,eigen) ! real (my_real)
real(my_real), intent(in) :: x,y
integer, intent(in) :: comp,eigen
real (my_real) :: trues
real (my_real) :: xmax,ymax,tmax,d,sq3
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

function truexs(x,y,comp,eigen) ! real (my_real)
real(my_real), intent(in) :: x,y
integer, intent(in) :: comp,eigen
real (my_real) :: truexs

real (my_real) :: xmax,ymax,tmax,d,sq3
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

function trueys(x,y,comp,eigen) ! real (my_real)
real(my_real), intent(in) :: x,y
integer, intent(in) :: comp,eigen
real (my_real) :: trueys

real (my_real) :: xmax,ymax,tmax,d,sq3

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

end module
