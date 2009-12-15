program phaml_master
use phaml
implicit none
type(phaml_solution_type) :: soln
call phaml_create(soln,nproc=4)
call phaml_solve_pde(soln,print_grid_when=PHASES,print_grid_who=MASTER, &
                     print_error_when=PHASES,print_error_what=ENERGY_L2_ERR, &
                     print_error_who=MASTER,max_vert=500, &
                     mg_tol=1.0e-8_my_real,mg_cycles=100)
call phaml_destroy(soln)
end program phaml_master
