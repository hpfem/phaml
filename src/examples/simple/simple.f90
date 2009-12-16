subroutine foo (a)
integer a
print*, "Hello from Fortran!"
print*, "a=",a
end

subroutine run()
use phaml
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Local variables

type(phaml_solution_type) :: soln
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

call phaml_destroy(soln)
end
