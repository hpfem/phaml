program phaml_master
use phaml
implicit none
type(phaml_solution_type) :: soln
call phaml_create(soln,nproc=4)
call phaml_solve_pde(soln,                   &
                     max_vert=500,           &
                     mg_cycles=2,            &
                     print_grid_when=PHASES, &
                     print_grid_who=MASTER_ALL  ,&
                     print_linsys_when=PHASES, &
                     print_linsys_who=MASTER_ALL, &
                     print_error_when=PHASES    , &
                     print_error_what=ENERGY_LINF_L2_ERR, &
                     print_errest_what=ENERGY_LINF_L2_ERREST, &
                     print_error_who=MASTER_ALL)
call phaml_destroy(soln)
end program phaml_master
