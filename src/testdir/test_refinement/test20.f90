program phaml_master
use phaml
implicit none
type(phaml_solution_type) :: soln
call phaml_create(soln,nproc=4)
call phaml_solve_pde(soln,print_grid_when=PHASES,print_grid_who=MASTER, &
                     max_eq=500,mg_cycles=2)
! TEMP have not decided how to handle refterm=KEEP_ERREST
!call phaml_solve_pde(soln,print_grid_when=PHASES,print_grid_who=MASTER, &
!                     max_refsolveloop=1, mg_cycles=2, &
!                     refterm=KEEP_ERREST)
call phaml_destroy(soln)
end program phaml_master
