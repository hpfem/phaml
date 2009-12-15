program phaml_master
use phaml
implicit none
type(phaml_solution_type) :: soln
call phaml_create(soln)
call phaml_solve_pde(soln,max_vert=100,print_grid_when=FINAL, &
                     print_grid_who=MASTER)
call phaml_destroy(soln)
end program phaml_master
