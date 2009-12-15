      subroutine cpu_second (cpu,user,sys)
      real cpu,user,sys
c
c  ------------------------------------------------
c  returns elapsed cpu time since start of job (sec)
c  ------------------------------------------------
c
c nil version; returns -1 for all clocks; can be used when there is no way
c to measure CPU time
c
      cpu = -1.
      user = -1.
      sys = -1.
c
      return
      end
