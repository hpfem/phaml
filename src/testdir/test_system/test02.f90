program phaml_master
use phaml
implicit none
type(phaml_solution_type) :: soln
real(my_real) :: err1, err2
call phaml_create(soln,nproc=4,system_size=2)
call phaml_solve_pde(soln,                   &
                     max_vert=300,           &
                     degree=3,               &
                     solver=CG_SOLVER,       &
                     print_grid_when=PHASES, &
                     print_grid_who=MASTER  ,&
                     print_error_when=PHASES    , &
                     print_error_what=ENERGY_LINF_ERR, &
                     print_errest_what=ENERGY_LINF_ERREST, &
                     print_error_who=MASTER)
call phaml_query(soln,linf_error=err1,comp=1)
call phaml_query(soln,linf_error=err2,comp=2)
write(6,"(A,SS,1P,2E19.12E2)") "Linf errors from query ",err1,err2
call phaml_destroy(soln)
end program phaml_master
