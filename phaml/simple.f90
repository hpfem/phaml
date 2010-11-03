module python_gets

use iso_c_utilities
use iso_c_binding
use phaml

implicit none

type(phaml_solution_type) :: this


contains

! c_phaml_* act like methods of a "Phaml" object, whose data is represented by
! the 'phaml_solution_type' object. The c_phaml_create() returns it, and then
! you pass it as the first parameter of all the other methods.

subroutine c_phaml_init(triangle_files, triangle_files_len, problem_number) &
    bind(c)
use example1
use example2
integer(c_int), intent(in) :: triangle_files_len
character(c_char), intent(in) :: triangle_files(triangle_files_len)
integer(c_int), intent(in) :: problem_number

select case(problem_number)
case (1); call setup_example1()
case (2); call setup_example2()
end select

call phaml_create(this, nproc=2, &
    triangle_files=char_array_to_string(triangle_files))
end subroutine

subroutine c_phaml_solve(term_energy_err) bind(c)
real(c_double), intent(in) :: term_energy_err
call phaml_solve_pde(this,                   &
                     term_energy_err=term_energy_err, &
                     print_grid_when=PHASES, &
                     print_grid_who=MASTER,  &
                     print_error_when=PHASES,&
                     reftype=HP_ADAPTIVE, &
                     refterm=DOUBLE_NELEM, &
                     print_time_when=FINAL, print_time_who=MASTER, &
                     print_error_what=ENERGY_LINF_ERR, &
                     print_errest_what=ENERGY_LINF_ERREST, &
                     print_error_who=MASTER)
end subroutine

subroutine c_phaml_get_mesh_info(n, nelem) bind(c)
integer(c_int), intent(out) :: n
integer(c_int), intent(out) :: nelem

n = size(this%grid%vertex)
nelem = this%grid%nelem_leaf
end subroutine

subroutine c_phaml_get_mesh(n, xvert, yvert, &
        nelem, element_vertices, orders) bind(c)
use global
use gridtype_mod
use phaml_type_mod
integer(c_int), intent(in) :: n
real(c_double), intent(out) :: xvert(n), yvert(n)
integer(c_int), intent(in) :: nelem
integer(c_int), intent(out) :: element_vertices(nelem, 3), orders(nelem)

type(phaml_solution_type), target :: soln
type(grid_type), pointer :: grid
integer :: ind, lev, elem

soln = this
grid => soln%grid

xvert = grid%vertex%coord%x
yvert = grid%vertex%coord%y

ind = 0
do lev = 1, grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (grid%element(elem)%isleaf) then
        ind = ind + 1
        element_vertices(ind,:) = grid%element(elem)%vertex
        orders(ind) = grid%element(elem)%degree
      endif
      elem = grid%element(elem)%next
   end do
end do

end subroutine


subroutine c_phaml_get_solution_values(n, x, y, values) bind(c)
integer(c_int) :: n
real(c_double) :: x(n), y(n), values(n)

integer :: i

call phaml_evaluate(this, x, y, values)
end subroutine


subroutine c_phaml_del() bind(c)
call phaml_destroy(this)
end subroutine

end module python_gets
