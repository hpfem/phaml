program phaml_master
use phaml
use phaml_user_mod
implicit none
type(phaml_solution_type) :: soln
call phaml_create(soln,nproc=4)
power=4
call update_usermod(soln)
call phaml_solve_pde(soln,print_grid_when=FINAL,print_grid_who=MASTER, &
                     print_error_when=FINAL,print_error_who=MASTER, &
                     print_error_what=ENERGY_LINF_L2_ERR, &
                     max_vert=500,mg_cycles=2)
call phaml_destroy(soln)
end program phaml_master
