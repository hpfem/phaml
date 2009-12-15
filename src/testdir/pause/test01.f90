program phaml_master
use phaml
implicit none
type(phaml_solution_type) :: soln
call phaml_create(soln,nproc=4)
call phaml_solve_pde(soln,max_refsolveloop=2,print_grid_when=PHASES, &
                     print_grid_who=MASTER, pause_at_start=.true., &
                     pause_after_phases=.true.,pause_at_end=.true.)
call phaml_destroy(soln)
end program phaml_master
