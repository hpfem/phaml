      subroutine cpu_second (cpu,user,sys)
      real cpu,user,sys
c
c  ------------------------------------------------
c  returns elapsed cpu time since start of job (sec)
C  ------------------------------------------------
C
C AIX version
C
      integer mclock
      cpu = mclock()/100.
      user = -1.
      sys = -1.
c
      return
      end
