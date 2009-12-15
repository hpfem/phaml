program phaml_master
use phaml
implicit none
type(phaml_solution_type) :: soln
call phaml_create(soln,nproc=4)
call phaml_solve_pde(soln,print_grid_when=FREQUENTLY,print_grid_who=MASTER_ALL, &
                     max_vert=1000, &
                     degree=4, &
                     prebalance=BALANCE_ELEMENTS, &
                     postbalance=BALANCE_EQUATIONS)
call phaml_destroy(soln)
end program phaml_master
