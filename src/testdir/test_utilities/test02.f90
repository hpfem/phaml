program phaml_master
use phaml
implicit none
type(phaml_solution_type) :: soln
real(my_real) :: intgrl, invnorm
call phaml_create(soln,nproc=4)
call phaml_solve_pde(soln,print_grid_when=FINAL,print_grid_who=MASTER, &
                     max_vert=500)
intgrl = phaml_integrate(soln,1)
write(6,100) "integral of solution",intgrl
intgrl = phaml_integrate(soln,2)
write(6,100) "integral of x*solution",intgrl
intgrl = phaml_integrate(soln,1,p=2)
write(6,100) "integral of solution squared",intgrl
invnorm = sqrt(1/intgrl)
call phaml_scale(soln,invnorm)
intgrl = phaml_integrate(soln,1,p=2)
write(6,100) "integral of normalized solution squared",intgrl
call phaml_destroy(soln)
100 format(A,SS,1P,E19.12E2)
end program phaml_master
