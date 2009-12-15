program phaml_master
use phaml
implicit none
type(phaml_solution_type) :: soln
call phaml_create(soln,nproc=4)
call phaml_solve_pde(soln,print_grid_when=PHASES,print_grid_who=MASTER_ALL, &
                     max_vert=1000, &
                     prebalance=BALANCE_ELEMENTS, &
                     postbalance=BALANCE_NONE)
call phaml_destroy(soln)
end program phaml_master
