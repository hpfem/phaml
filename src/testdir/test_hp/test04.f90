program phaml_master
use phaml
implicit none
type(phaml_solution_type) :: soln
call phaml_create(soln,nproc=4)
call phaml_solve_pde(soln,print_grid_when=PHASES,print_grid_who=MASTER, &
                     max_eq=4000,reftype=HP_ADAPTIVE,hp_strategy=HP_PRIOR2P_H1, &
                     degree=3,derefine=.false.)
call phaml_destroy(soln)
end program phaml_master
