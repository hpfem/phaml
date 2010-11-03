module example2

use phaml
implicit none

private
public setup_example2

real(my_real), save :: r0 = 0.7_my_real, &
                       xc = 0.0_my_real, &
                       yc = 0.0_my_real, &
                       alpha = 20.0_my_real

contains

subroutine setup_example2()
use pde_pointers
call set_pde_pointers(pdecoefs, bconds, iconds, trues, truexs, trueys)
end subroutine

subroutine pdecoefs(x,y,cxx,cxy,cyy,cx,cy,c,rs)
real(my_real), intent(in) :: x,y
real(my_real), intent(out) :: cxx(:,:),cxy(:,:),cyy(:,:),cx(:,:),cy(:,:), &
                              c(:,:),rs(:)

real(my_real) :: rx2, ry2, r, rr, denom, denom2

cxx(1,1) = 1.0_my_real
cyy(1,1) = 1.0_my_real
c(1,1) = 0.0_my_real

rx2 = (x-xc)**2
ry2 = (y-yc)**2
r = sqrt(rx2+ry2)
rr = r-r0
denom = ((1+(alpha*rr)**2))**2
denom2 = r*(1+(alpha*rr)**2)

rs(1) = 2*alpha**3*rr/denom - alpha/denom2

cxy=0; cx=0; cy=0

end subroutine pdecoefs

subroutine bconds(x,y,bmark,itype,c,rs)
real(my_real), intent(in) :: x,y
integer, intent(in) :: bmark
integer, intent(out) :: itype(:)
real(my_real), intent(out) :: c(:,:),rs(:)
itype = DIRICHLET
c = 0.0_my_real
rs = trues(x,y,1,1)
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
trues = atan(alpha*(sqrt((x-xc)**2 + (y-yc)**2)-r0))
end function trues

function truexs(x,y,comp,eigen) ! real (my_real)
real(my_real), intent(in) :: x,y
integer, intent(in) :: comp,eigen
real (my_real) :: truexs

real(my_real) :: r
r = sqrt((x-xc)**2+(y-yc)**2)
truexs = alpha*(x-xc)/(r*(1+(alpha*(r-r0))**2))
end function truexs

function trueys(x,y,comp,eigen) ! real (my_real)
real(my_real), intent(in) :: x,y
integer, intent(in) :: comp,eigen
real (my_real) :: trueys

real(my_real) :: r
r = sqrt((x-xc)**2+(y-yc)**2)
trueys = alpha*(y-yc)/(r*(1+(alpha*(r-r0))**2))
end function trueys


end module
