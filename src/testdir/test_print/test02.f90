program phaml_master
use phaml
implicit none
type(phaml_solution_type) :: soln
call phaml_create(soln,nproc=4)
call phaml_solve_pde(soln,                   &
                     max_vert=500,           &
                     mg_cycles=2,            &
                     print_grid_when=FINAL, &
                     print_grid_who=MASTER  ,&
                     print_linsys_when=NEVER, &
                     print_linsys_who=MASTER, &
                     print_header_who=NO_ONE, &
                     print_trailer_who=NO_ONE, &
                     print_error_when=FINAL    , &
                     print_error_what=ENERGY_LINF_ERR, &
                     print_errest_what=ENERGY_LINF_ERREST, &
                     print_error_who=MASTER)
call phaml_destroy(soln)
end program phaml_master
