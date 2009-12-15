   subroutine cpu_second(cpu,user,sys)
! f95 version
   implicit none
   real, intent(out) :: cpu, user, sys
   call cpu_time(cpu)
   user = -1.
   sys = -1.
   end subroutine cpu_second

