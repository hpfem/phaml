program phaml_master
use phaml
use phaml_user_mod
implicit none
type(phaml_solution_type) :: soln
probno=4
call phaml_create(soln,nproc=1,update_umod=.true.)
call phaml_solve_pde(soln,print_grid_when=PHASES,print_grid_who=MASTER, &
                     print_error_when=PHASES,print_error_who=MASTER, &
                     print_error_what=ENERGY_LINF_L2_ERR,max_vert=500, &
                     ignore_quad_err=.true., &
                     solver=LAPACK_INDEFINITE_SOLVER)
call phaml_destroy(soln)
end program phaml_master
