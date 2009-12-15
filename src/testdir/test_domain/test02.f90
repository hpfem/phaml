program phaml_master
use phaml
implicit none
type(phaml_solution_type) :: soln
call phaml_create(soln,nproc=4,triangle_files="stadium",max_blen=.25_my_real, &
                  id=2)
call phaml_solve_pde(soln,print_grid_when=PHASES,print_grid_who=MASTER, &
                     print_error_when=FINAL,print_error_what=LINF_ERR, &
                     print_error_who=MASTER,max_elem=2000,mg_cycles=20)
call phaml_destroy(soln)
end program phaml_master
