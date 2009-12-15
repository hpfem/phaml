program phaml_master
use phaml
implicit none
type(phaml_solution_type) :: soln
call phaml_create(soln)
call phaml_solve_pde(soln,max_refsolveloop=1,degree=3,task=REFINE_ONLY)
open(unit=21,file="savestiff.dat")
open(unit=22,file="saverhs.dat")
call phaml_store_matrix(soln,stiffness_unit=21,rhs_unit=22)
close(21)
close(22)
call phaml_destroy(soln,finalize_mpi=.false.)
call phaml_create(soln,eq_type=EIGENVALUE)
call phaml_solve_pde(soln,max_refsolveloop=1,degree=3,task=REFINE_ONLY)
open(unit=23,file="savemass.dat")
call phaml_store_matrix(soln,mass_unit=23)
close(23)
call phaml_destroy(soln)

end program phaml_master
