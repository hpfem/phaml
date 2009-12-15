program phaml_master
use phaml
use phaml_user_mod
implicit none
type(phaml_solution_type) :: soln
real(my_real) :: x(0:20),y(0:20),u(0:20),ux(0:20),uy(0:20),uxx(0:20),uyy(0:20)
integer :: i
call phaml_create(soln,nproc=4)
power=4
call update_usermod(soln)
call phaml_solve_pde(soln,max_vert=500,degree=3,mg_cycles=20)
x = (/ (i/20.0_my_real,i=0,20) /)
y = 0.4_my_real
call phaml_evaluate(soln,x,y,u,ux,uy,uxx,uyy)
do i=0,20
   write(6,"(SS,1P,2E19.12E2)") x(i),u(i)
end do
do i=0,20
   write(6,"(SS,1P,3E19.12E2)") x(i),ux(i),uy(i)
end do
do i=0,20
   write(6,"(SS,1P,3E19.12E2)") x(i),uxx(i),uyy(i)
end do
call phaml_destroy(soln)
end program phaml_master
