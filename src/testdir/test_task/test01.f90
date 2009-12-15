program phaml_master
use phaml
implicit none
type(phaml_solution_type) :: soln
call phaml_create(soln,nproc=4)
call phaml_solve_pde(soln,print_grid_when=PHASES,print_grid_who=MASTER_ALL, &
                     print_error_when=PHASES,print_error_who=MASTER, &
                     print_error_what=ENERGY_LINF_L2_ERR,max_vert=100, &
                     sequential_vert=1000,task=SET_INITIAL, &
                     error_estimator=INITIAL_CONDITION)
call phaml_solve_pde(soln,print_grid_when=FINAL,print_grid_who=MASTER_ALL, &
                     sequential_vert=10,task=BALANCE_ONLY)
call phaml_solve_pde(soln,print_grid_when=PHASES,print_grid_who=MASTER_ALL, &
                     print_error_when=PHASES,print_error_who=MASTER, &
                     print_error_what=ENERGY_LINF_L2_ERR, &
                     task=REFINE_ONLY)
call phaml_solve_pde(soln,print_grid_when=PHASES,print_grid_who=MASTER, &
                     print_error_when=FINAL,print_error_who=MASTER, &
                     print_error_what=ENERGY_LINF_L2_ERR, &
                     task=SOLVE_ONLY,mg_cycles=5)
call phaml_destroy(soln)
end program phaml_master
