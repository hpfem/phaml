program phaml_master
use phaml
implicit none
type(phaml_solution_type) :: soln
call phaml_create(soln,nproc=4)
call phaml_solve_pde(soln,                   &
                     max_vert=500,           &
                     mg_cycles=2,            &
                     print_grid_when=FREQUENTLY, &
                     print_grid_who=MASTER  ,&
                     print_linsys_when=FREQUENTLY, &
                     print_linsys_who=MASTER, &
                     print_error_when=FREQUENTLY    , &
                     print_error_what=LINF_L2_ERR, &
                     print_errest_what=LINF_L2_ERREST, &
                     print_error_who=MASTER)
call phaml_destroy(soln)
end program phaml_master
