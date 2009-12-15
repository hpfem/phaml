program phaml_master
use phaml
implicit none
type(phaml_solution_type) :: soln
call phaml_create(soln,nproc=2,output_unit=11)
call phaml_popen(soln,11,"test06..out")
call phaml_solve_pde(soln,                   &
                     max_vert=500,           &
                     mg_cycles=2,            &
                     print_grid_when=FINAL, &
                     print_grid_who=SLAVES  ,&
                     print_linsys_when=NEVER, &
                     print_linsys_who=SLAVES, &
                     print_header_who=SLAVES, &
                     print_trailer_who=SLAVES, &
                     print_error_when=FINAL    , &
                     print_error_what=L2_ERR, &
                     print_errest_what=L2_ERREST, &
                     print_error_who=SLAVES)
call phaml_pclose(soln,11)
call phaml_destroy(soln)
end program phaml_master
