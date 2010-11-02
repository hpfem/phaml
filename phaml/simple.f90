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

subroutine c_phaml_init(triangle_files, triangle_files_len) bind(c)
integer(c_int) :: triangle_files_len
character(c_char) :: triangle_files(triangle_files_len)

call phaml_create(this, nproc=2, &
    triangle_files=char_array_to_string(triangle_files))
end subroutine

subroutine c_phaml_solve( n, x, y, sol) bind(c)
integer(c_int) :: n
real(c_double) :: x(n), y(n), sol(n)
call phaml_solve_pde(this,                   &
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
end subroutine


subroutine c_phaml_get_mesh() bind(c)
use global
use gridtype_mod
use phaml_type_mod
type(phaml_solution_type), target :: soln
type(grid_type), pointer :: grid
integer :: ind, lev, elem
double precision, allocatable :: xvert(:), yvert(:)
integer, allocatable :: element_vertices(:,:), element_order(:)
integer :: nvert, nelem


soln = this
grid => soln%grid

nvert = grid%nvert
nelem = grid%nelem_leaf

allocate(xvert(size(grid%vertex)),yvert(size(grid%vertex)), &
         element_vertices(nelem,3),element_order(nelem))

xvert = grid%vertex%coord%x
yvert = grid%vertex%coord%y

ind = 0
do lev = 1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (grid%element(elem)%isleaf) then
        ind = ind + 1
        element_vertices(ind,:) = grid%element(elem)%vertex
        element_order(ind) = grid%element(elem)%degree
        print *,"elem ind vertices ",elem,ind,element_vertices(ind,:)
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
