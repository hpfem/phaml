program phaml_master
use phaml
implicit none
type(phaml_solution_type) :: soln
call phaml_create(soln,nproc=4,draw_grid_who=MASTER)
call phaml_solve_pde(soln,max_vert=500,print_grid_when=PHASES, &
                     print_grid_who=MASTER,draw_grid_when=PHASES, &
                     pause_after_draw=.true.,print_header_who=NO_ONE, &
                     print_error_when=PHASES,print_error_who=MASTER, &
                     print_error_what=LINF_ERR,mg_cycles=2, &
                     sequential_vert=20,inc_factor=8.0_my_real)
call phaml_destroy(soln)
end program phaml_master
