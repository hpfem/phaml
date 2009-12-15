program phaml_master
use phaml
implicit none
type(phaml_solution_type) :: soln
call phaml_create(soln,nproc=1,draw_grid_who=MASTER)
call phaml_solve_pde(soln,max_vert=1000,print_grid_when=FINAL, &
                     print_grid_who=MASTER,draw_grid_when=FINAL, &
                     pause_at_end=.true.,print_header_who=NO_ONE,mg_cycles=2)
call phaml_destroy(soln)
end program phaml_master
