program phaml_master
use phaml
implicit none
type(phaml_solution_type) :: soln
call phaml_create(soln,nproc=2)
call phaml_solve_pde(soln,max_vert=500,print_grid_when=FINAL, &
                     print_grid_who=MASTER)
call phaml_compress(soln)
call phaml_popen(soln,11,"savef.dat")
call phaml_store(soln,11)
call phaml_pclose(soln,11)
call phaml_destroy(soln,finalize_mpi=.false.)
call phaml_create(soln,nproc=2)
call phaml_popen(soln,12,"savef.dat")
call phaml_restore(soln,12)
call phaml_pclose(soln,12)
call phaml_solve_pde(soln,max_vert=500,print_grid_when=FINAL, &
                     print_grid_who=MASTER,inc_factor=1.1_my_real)
call phaml_compress(soln)
call phaml_popen(soln,13,"saveu.dat","UNFORMATTED")
call phaml_store(soln,13)
call phaml_pclose(soln,13)
call phaml_destroy(soln,finalize_mpi=.false.)
call phaml_create(soln,nproc=2)
call phaml_popen(soln,14,"saveu.dat","UNFORMATTED")
call phaml_restore(soln,14)
call phaml_pclose(soln,14)
call phaml_solve_pde(soln,max_vert=500,print_grid_when=FINAL, &
                     print_grid_who=MASTER,inc_factor=1.1_my_real)
call phaml_destroy(soln)

end program phaml_master
