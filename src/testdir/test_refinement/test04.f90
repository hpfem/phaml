program phaml_master
use phaml
implicit none
type(phaml_solution_type) :: soln
call phaml_create(soln,nproc=4)
call phaml_solve_pde(soln,print_grid_when=FINAL,print_grid_who=MASTER, &
                     max_refsolveloop=3)
call phaml_solve_pde(soln,print_grid_when=PHASES,print_grid_who=MASTER, &
                     max_eq=1000, &
                     reftype=P_ADAPTIVE,error_estimator=LOCAL_PROBLEM_P)
call phaml_destroy(soln)
end program phaml_master
