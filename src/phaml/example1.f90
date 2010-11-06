module example1

use phaml
implicit none

private
public setup_example1

contains

subroutine setup_example1()
use pde_pointers
call set_pde_pointers(pdecoefs, bconds, iconds, trues, truexs, trueys)
end subroutine

subroutine pdecoefs(x,y,cxx,cxy,cyy,cx,cy,c,rs)
real(my_real), intent(in) :: x,y
real(my_real), intent(out) :: cxx(:,:),cxy(:,:),cyy(:,:),cx(:,:),cy(:,:), &
                              c(:,:),rs(:)
cxx(1,1) = 1.0_my_real
cyy(1,1) = 1.0_my_real
c(1,1) = 0.0_my_real
rs(1) = -90*x**8 - 90*y**8

cxy=0; cx=0; cy=0
end subroutine pdecoefs

subroutine bconds(x,y,bmark,itype,c,rs)
real(my_real), intent(in) :: x,y
integer, intent(in) :: bmark
integer, intent(out) :: itype(:)
real(my_real), intent(out) :: c(:,:),rs(:)
itype = DIRICHLET
c = 0.0_my_real
rs(1) = trues(x,y,1,1)
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
trues = x**10 + y**10
end function trues

function truexs(x,y,comp,eigen) ! real (my_real)
real(my_real), intent(in) :: x,y
integer, intent(in) :: comp,eigen
real (my_real) :: truexs
truexs = 10*x
end function truexs

function trueys(x,y,comp,eigen) ! real (my_real)
real(my_real), intent(in) :: x,y
integer, intent(in) :: comp,eigen
real (my_real) :: trueys
trueys = 10*y
end function trueys


end module
