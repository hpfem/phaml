program phaml_master
use phaml
implicit none
type(phaml_solution_type) :: soln
call phaml_create(soln,nproc=2,spawn_form=DEBUG_SLAVE)
call phaml_solve_pde(soln,max_refsolveloop=2,print_grid_when=PHASES, &
                     print_grid_who=SLAVES)
call phaml_destroy(soln)
end program phaml_master
