program phaml_master
use phaml
implicit none
type(phaml_solution_type) :: soln
call phaml_create(soln,nproc=2,output_unit=11)
call phaml_popen(soln,11,"test07..out")
call phaml_solve_pde(soln,                   &
                     max_vert=500,           &
                     mg_cycles=2,            &
                     print_grid_when=PHASES, &
                     print_grid_who=EVERYONE  ,&
                     print_linsys_when=PHASES, &
                     print_linsys_who=EVERYONE, &
                     print_header_who=EVERYONE, &
                     print_trailer_who=EVERYONE, &
                     print_error_when=PHASES    , &
                     print_error_what=LINF_ERR, &
                     print_errest_what=LINF_ERREST, &
                     print_error_who=EVERYONE)
call phaml_pclose(soln,11)
call phaml_destroy(soln)
end program phaml_master
