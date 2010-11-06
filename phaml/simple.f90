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
use example_eigen
integer(c_int), intent(in) :: triangle_files_len
character(c_char), intent(in) :: triangle_files(triangle_files_len)
integer(c_int), intent(in) :: problem_number

integer :: eq_type = ELLIPTIC

select case(problem_number)
case (1)
    call setup_example1()
case (2)
    call setup_example2()
case (3)
    call setup_example_eigen()
    eq_type = EIGENVALUE
case default; stop "Unknown problem_number"
end select

call phaml_create(this, nproc=1, eq_type=eq_type, &
    triangle_files=char_array_to_string(triangle_files))
end subroutine

subroutine c_phaml_solve(term_energy_err, max_eq, verbose, reftype, &
    hp_strategy, derefine, degree, error_estimator, lambda0, lambda_smallest, &
    num_eval) bind(c)
real(c_double), intent(in) :: term_energy_err
integer(c_int), intent(in) :: max_eq
integer(c_int), intent(in) :: verbose
integer(c_int), intent(in) :: reftype
integer(c_int), intent(in) :: hp_strategy
integer(c_int), intent(in) :: derefine
integer(c_int), intent(in) :: degree
integer(c_int), intent(in) :: error_estimator
real(c_double), intent(in) :: lambda0
integer(c_int), intent(in) :: lambda_smallest
integer(c_int), intent(in) :: num_eval

integer :: print_grid_when
integer :: print_error_when
integer :: print_header_who
integer :: print_trailer_who
integer :: print_time_when
real(my_real) :: lambda

if (verbose == 1) then
    print_grid_when = PHASES
    print_error_when = PHASES
    print_header_who = MASTER
    print_trailer_who = MASTER
    print_time_when = FINAL
else
    print_grid_when = NEVER
    print_error_when = NEVER
    print_header_who = NO_ONE
    print_trailer_who = NO_ONE
    print_time_when = NEVER
endif

if (lambda_smallest == 1) then
    lambda = -huge(0.0_my_real)
else
    lambda = lambda0
endif

call phaml_solve_pde(this,                       &
                     term_energy_err=term_energy_err, &
                     errtype=RELATIVE_ERROR, &
                     max_eq=max_eq, &
                     reftype=reftype, &
                     refterm=ONE_REF_HALF_ERRIND, &
                     error_estimator=error_estimator, &
                     derefine=(derefine==1), &
                     hp_strategy=hp_strategy, &
                     degree=degree, &
                     mg_cycles=100, &
                     max_lev=53, &
                     lambda0 = lambda,          &
                     num_eval = num_eval, &

                     print_grid_when=print_grid_when, &
                     print_grid_who=MASTER,  &
                     print_error_when=print_error_when,&
                     print_header_who=print_header_who, &
                     print_trailer_who=print_trailer_who, &
                     print_eval_when=PHASES, &
                     print_eval_who=MASTER, &
                     print_time_when=print_time_when, print_time_who=MASTER, &
                     print_error_what=ENERGY_LINF_L2_ERR, &
                     print_errest_what=ENERGY_LINF_L2_ERREST, &
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

call phaml_evaluate(this, x, y, values, eigen=1)
end subroutine


subroutine c_phaml_del() bind(c)
call phaml_destroy(this)
end subroutine

end module python_gets
