module python_gets

double precision, allocatable :: xvert(:), yvert(:)
integer, allocatable :: element_vertices(:,:), element_order(:)
integer :: nvert, nelem

contains

subroutine c_run(n, x, y, sol) bind(c)

use phaml
use iso_c_binding
use iso_c_utilities
implicit none

interface
   subroutine get_grid_params(soln)
   use global
   use gridtype_mod
   use phaml_type_mod
   !use python_defs
   type(phaml_solution_type), intent(in), target :: soln
   end subroutine get_grid_params
end interface
!----------------------------------------------------


integer(c_int) :: n
real(c_double) :: x(n), y(n), sol(n)
!f2py intent(in) x
!f2py intent(in) y
!f2py intent(out) sol
!----------------------------------------------------
! Local variables

type(phaml_solution_type) :: soln
integer :: i
! doesn't work:
!character(len=*) :: triangle_files = "domain"
!real(my_real), allocatable :: x(:), y(:), u(:)

!----------------------------------------------------
! Begin executable code

print*, "Start"

call phaml_create(soln,nproc=2,triangle_files="domain")

call phaml_solve_pde(soln,                   &
                     max_vert=100,          &
                     print_grid_when=PHASES, &
                     print_grid_who=MASTER,  &
                     print_error_when=PHASES,&
		     reftype=HP_ADAPTIVE, &
		     refterm=DOUBLE_NELEM, &
                     print_time_when=FINAL, print_time_who=MASTER, &
                     print_error_what=ENERGY_LINF_ERR, &
                     print_errest_what=ENERGY_LINF_ERREST, &
                     print_error_who=MASTER)

!allocate(x(3), y(3), u(3))
!x(1) = 0
!x(2) = 0.5
!x(3) = 1
!y(1) = 0
!y(2) = 0.5
!y(3) = 1

call phaml_evaluate(soln, x, y, sol)

print*, "x=", x
print*, "y=", y
print*, "u=", sol

call get_grid_params(soln)
print *,"xvert ",xvert(1:5)
print *,"yvert ",yvert(1:5)
do i=1,5
   print *,"elem verts ",element_vertices(i,:)
end do
print *,"elem deg ", element_order(1:5)
do i=nelem-5,size(element_vertices,dim=1)
   print *,"elem verts last ",element_vertices(i,:)
end do
print *,"elem deg last ", element_order(nelem-5:)

call phaml_destroy(soln)
end subroutine c_run

end module python_gets
