program phaml_master
use phaml
use phaml_user_mod
implicit none
type(phaml_solution_type) :: soln
probno=1
call phaml_create(soln,nproc=4,update_umod=.true.,triangle_files="period.1")
call phaml_solve_pde(soln,print_grid_when=PHASES,print_grid_who=MASTER, &
                     print_error_when=PHASES,print_error_who=MASTER, &
                     print_error_what=ENERGY_LINF_L2_ERR,max_vert=500, &
                     print_warnings=.false., &
                     refterm=DOUBLE_NVERT,mg_cycles=5)
call phaml_destroy(soln)
end program phaml_master
