program phaml_master
use phaml
implicit none
type(phaml_solution_type) :: soln
integer :: nvert,nelem,neq,nlev,min_degree,max_degree
integer, dimension(4) :: nvert_proc,nvert_own,nelem_proc,nelem_own, &
                         neq_proc,neq_own
real(my_real) :: linf_error,energy_error,l2_error, max_error_indicator, &
                 linf_error_estimate,energy_error_estimate, l2_error_estimate, &
                 linf_solution, l2_solution, energy_solution, linf_u, l2_u, &
                 energy_u, linf_true, l2_true, energy_true
call phaml_create(soln,nproc=4)
call phaml_solve_pde(soln,print_grid_when=FINAL,print_grid_who=MASTER, &
                     max_eq=500)
nvert_proc=0; nvert_own=0; nelem_proc=0; nelem_own=0; neq_proc=0; neq_own=0
call phaml_query(soln,nvert,nvert_proc,nvert_own,nelem,nelem_proc,nelem_own, &
                 neq,neq_proc,neq_own,nlev,min_degree,max_degree,linf_error, &
                 energy_error,l2_error,max_error_indicator, &
                 linf_error_estimate,energy_error_estimate,l2_error_estimate, &
                 linf_solution,l2_solution,energy_solution,linf_u,l2_u, &
                 energy_u,linf_true,l2_true,energy_true)
call phaml_destroy(soln)
write(6,100) "nvert ",nvert
write(6,101) "nvert_proc ",nvert_proc
write(6,101) "nvert_own ",nvert_own
write(6,100) "nelem ",nelem
write(6,101) "nelem_proc ",nelem_proc
write(6,101) "nelem_own ",nelem_own
write(6,100) "neq ",neq
write(6,101) "neq_proc ",neq_proc
write(6,101) "neq_own ",neq_own
write(6,100) "nlev ",nlev
write(6,100) "min_degree ",min_degree
write(6,100) "max_degree ",max_degree
write(6,102) "linf_error ",linf_error
write(6,102) "energy_error ",energy_error
write(6,102) "l2_error ",l2_error
write(6,102) "max_error_indicator ",max_error_indicator
write(6,102) "linf_error_estimate ",linf_error_estimate
write(6,102) "energy_error_estimate ",energy_error_estimate
write(6,102) "l2_error_estimate ",l2_error_estimate
write(6,102) "linf_solution ",linf_solution
write(6,102) "l2_solution ",l2_solution
write(6,102) "energy_solution ",energy_solution
write(6,102) "linf_u ",linf_u
write(6,102) "l2_u ",l2_u
write(6,102) "energy_u ",energy_u
write(6,102) "linf_true ",linf_true
write(6,102) "l2_true ",l2_true
write(6,102) "energy_true ",energy_true

100 format(A,I11)
101 format(A,4I11)
102 format(A,SS,1P,E19.12E2)
end program phaml_master
