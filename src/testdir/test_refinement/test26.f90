program phaml_master
use phaml
implicit none
type(phaml_solution_type) :: soln
call phaml_create(soln,nproc=4)
call phaml_solve_pde(soln,print_grid_when=PHASES,print_grid_who=MASTER, &
                     max_eq=500, mg_cycles=2, &
                     inc_factor=4.0_my_real,refterm=DOUBLE_NVERT)
call phaml_destroy(soln)
end program phaml_master
