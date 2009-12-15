program phaml_master
use phaml
implicit none
type(phaml_solution_type) :: soln
call phaml_create(soln,nproc=4)
call phaml_solve_pde(soln,max_vert=1000,print_grid_when=FINAL, &
                     print_grid_who=MASTER)
call phaml_destroy(soln)
end program phaml_master
