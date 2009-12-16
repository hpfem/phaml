subroutine foo (a)
integer a
print*, "Hello from Fortran!"
print*, "a=",a
end

subroutine run(n, x, y, sol)
use phaml
!----------------------------------------------------

implicit none

integer n
double precision x(n), y(n), sol(n)
!f2py intent(in) x
!f2py intent(in) y
!f2py intent(out) sol
!----------------------------------------------------
! Local variables

type(phaml_solution_type) :: soln
!real(my_real), allocatable :: x(:), y(:), u(:)

!----------------------------------------------------
! Begin executable code

print*, "Start"

call phaml_create(soln,nproc=2)

call phaml_solve_pde(soln,                   &
                     max_vert=1000,          &
                     print_grid_when=PHASES, &
                     print_grid_who=MASTER,  &
                     print_error_when=PHASES,&
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

call phaml_destroy(soln)
end
