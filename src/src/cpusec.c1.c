/* measure execution time in seconds */
/* Bill Mitchell 10/9/92             */
/*                                   */
/* Use unix routine"times"to measure */
/* user execution time between       */
/* calls to second in seconds.       */
/* Callable from FORTRAN and C.      */
/* Returns real                      */
/* in FORTRAN and float in C.        */

static int holdtime,flagtime=1;

/* FORTRAN callable version on systems that covert to upper case */
void CPU_SECOND(c,u,s)
float *c,*u,*s;
{
/* just call C version */
   void cpu_second();
   cpu_second(c,u,s);
   return;
}

/* FORTRAN callable version on systems that need added underscore */
void cpu_second_(c,u,s)
float *c,*u,*s;
{
/* just call C version */
   void cpu_second();
   cpu_second(c,u,s);
   return;
}

/* FORTRAN callable version on systems that need two underscores */
void cpu_second__(c,u,s)
float *c,*u,*s;
{
/* just call C version */
   void cpu_second();
   cpu_second(c,u,s);
   return;
}

/* C callable version and FORTRAN callable on systems that don't name mangle */
void cpu_second(c,u,s)
float *c,*u,*s;
{
#include <sys/types.h>
#include <sys/times.h>
#include <sys/param.h>
   clock_t times();
   clock_t t;
   struct tms t1;

/* call "times" */
   t=times(&t1);
   if (flagtime) {
      flagtime=0;
      holdtime=t;
   }
/* user time in 1/HZ seconds is in tms_utime */
/* HZ is in sys/param.h */
   *u = ((float)t1.tms_utime)/((float)HZ);
   *s = ((float)t1.tms_stime)/((float)HZ);
   *c = *u + *s;
   return;
}
