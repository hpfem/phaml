program phaml_master
use phaml
implicit none
type(phaml_solution_type) :: soln
call phaml_create(soln,nproc=1)
call phaml_solve_pde(soln,print_grid_when=PHASES,print_grid_who=MASTER, &
                     max_eq=500,reftype=HP_ADAPTIVE, &
                     hp_strategy=HP_REFSOLN_ELEM)
call phaml_destroy(soln)
end program phaml_master
