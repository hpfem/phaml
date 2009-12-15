program phaml_master
use phaml
implicit none
type(phaml_solution_type) :: soln
character(len=128) :: dbg
open(10,file="dbg3_command")
read(10,*) dbg
close(10)
call phaml_create(soln,nproc=2,spawn_form=DEBUG_SLAVE,debug_command=dbg)
call phaml_solve_pde(soln,max_refsolveloop=2,print_grid_when=PHASES, &
                     print_grid_who=SLAVES)
call phaml_destroy(soln)
end program phaml_master
