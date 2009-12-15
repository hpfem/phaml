/*--------------------------------------------------------------------!
!                                PHAML                                !
!                                                                     !
! The Parallel Hierarchical Adaptive MultiLevel code for solving      !
! linear elliptic partial differential equations of the form          !
! (PUx)x + (QUy)y + RU = F on 2D polygonal domains with mixed         !
! boundary conditions, and eigenvalue problems where F is lambda*U.   !
!                                                                     !
! PHAML is public domain software.  It was produced as part of work   !
! done by the U.S. Government, and is not subject to copyright in     !
! the United States.                                                  !
!                                                                     !
!     William F. Mitchell                                             !
!     Mathematical and Computational Sciences Division                !
!     National Institute of Standards and Technology                  !
!     william.mitchell@nist.gov                                       !
!     http://math.nist.gov/phaml                                      !
!                                                                     !
!--------------------------------------------------------------------*/

/* Different versions of hypre need different parts of the fixes in
   this file.  Look for "Section" to check whether that section should
   be used or commented out for your version.  I only know the answer for
   versions 1.6.0, 1.9.0b and 2.0.0.  For other versions you might have
   to experiment.  The current setup is for version 2.0.0
*/

#include <mpi.h>

/* Section on MPI communicator.  1.6.0, 1.9.0b and 2.0.0 need this section */

/* hypre casts the Fortran MPI communicator to a C communicator, which does
   not work with all MPI libraries.  MPI-2 provides a routine to do the
   conversion, but not all MPI libraries support MPI-2 yet.  This file
   provides an additional layer of wrappers for the relevant hypre subroutines
   to convert the communicator (HYPRE_FIX=HFMPI2) or leave it alone
   (HYPRE_FIX=HFASIS).  Other conversions can be added if needed.
*/

/* the second problem is that the fortran interface for
   HYPRE_ParCSRPCGSetLogging is missing in hypre 1.9.0b */

#define HFMPI2 1
#define HFASIS 2

/* I'm being lazy and just selecting the name mangling form here.  Define the
   string that selects the correct form for the Fortran compiler being used.
   I should set this in mkmkfile.sh. */

#define UNDERSCORE

#ifdef UNDERSCORE
#define my_hypre_ijmatrixcreate     my_hypre_ijmatrixcreate_
#define my_hypre_ijvectorcreate     my_hypre_ijvectorcreate_
#define my_hypre_parcsrpcgcreate    my_hypre_parcsrpcgcreate_
#define my_hypre_parcsrgmrescreate  my_hypre_parcsrgmrescreate_
#define my_hypre_parasailscreate    my_hypre_parasailscreate_
#define hypre_ijmatrixcreate        hypre_ijmatrixcreate_
#define hypre_ijvectorcreate        hypre_ijvectorcreate_
#define hypre_parcsrpcgcreate       hypre_parcsrpcgcreate_
#define hypre_parcsrgmrescreate     hypre_parcsrgmrescreate_
#define hypre_parasailscreate       hypre_parasailscreate_
#define my_hypre_boomeramgsetnumgridsweeps my_hypre_boomeramgsetnumgridsw_
#endif

#ifdef DOUBLE_UNDERSCORE
#define my_hypre_ijmatrixcreate     my_hypre_ijmatrixcreate__
#define my_hypre_ijvectorcreate     my_hypre_ijvectorcreate__
#define my_hypre_parcsrpcgcreate    my_hypre_parcsrpcgcreate__
#define my_hypre_parcsrgmrescreate  my_hypre_parcsrgmrescreate__
#define my_hypre_parasailscreate    my_hypre_parasailscreate__
#define hypre_ijmatrixcreate        hypre_ijmatrixcreate__
#define hypre_ijvectorcreate        hypre_ijvectorcreate__
#define hypre_parcsrpcgcreate       hypre_parcsrpcgcreate__
#define hypre_parcsrgmrescreate     hypre_parcsrgmrescreate__
#define hypre_parasailscreate       hypre_parasailscreate__
#define my_hypre_boomeramgsetnumgridsweeps my_hypre_boomeramgsetnumgridsw__
#endif

#ifdef CAPS
#define my_hypre_ijmatrixcreate     MY_HYPRE_IJMATRIXCREATE
#define my_hypre_ijvectorcreate     MY_HYPRE_IJVECTORCREATE
#define my_hypre_parcsrpcgcreate    MY_HYPRE_PARCSRPCGCREATE
#define my_hypre_parcsrgmrescreate  MY_HYPRE_PARCSRGMRESCREATE
#define my_hypre_parasailscreate    MY_HYPRE_PARASAILSCREATE
#define hypre_ijmatrixcreate        HYPRE_IJMATRIXCREATE
#define hypre_ijvectorcreate        HYPRE_IJVECTORCREATE
#define hypre_parcsrpcgcreate       HYPRE_PARCSRPCGCREATE
#define hypre_parcsrgmrescreate     HYPRE_PARCSRGMRESCREATE
#define hypre_parasailscreate       HYPRE_PARASAILSCREATE
#define my_hypre_boomeramgsetnumgridsweeps MY_HYPRE_BOOMERAMGSETNUMGRIDSW
#endif

MPI_Comm hypre_fix_comm_f2c(int *f_comm)
{
#if HYPRE_FIX==HFMPI2
   return MPI_Comm_f2c((MPI_Fint)(*f_comm));
#elif HYPRE_FIX==HFASIS
   return (MPI_Comm)(*f_comm);
#endif
}

void my_hypre_ijmatrixcreate(int *comm, int *ilower, int *iupper,
                             int *jlower, int *jupper, long int *matrix,
                             int *ierr)
{
   MPI_Comm new_comm;
   new_comm = hypre_fix_comm_f2c(comm);
   hypre_ijmatrixcreate(&new_comm, ilower, iupper, jlower, jupper, matrix,
                        ierr);
}

void my_hypre_ijvectorcreate(int *comm, int *jlower, int *jupper,
                             long int *vector, int *ierr)
{
   MPI_Comm new_comm;
   new_comm = hypre_fix_comm_f2c(comm);
   hypre_ijvectorcreate(&new_comm, jlower, jupper, vector, ierr);
}

void my_hypre_parcsrpcgcreate(int *comm, long int *solver, int *ierr)
{
   MPI_Comm new_comm;
   new_comm = hypre_fix_comm_f2c(comm);
   hypre_parcsrpcgcreate(&new_comm, solver, ierr);
}

void my_hypre_parcsrgmrescreate(int *comm, long int *solver, int *ierr)
{
   MPI_Comm new_comm;
   new_comm = hypre_fix_comm_f2c(comm);
   hypre_parcsrgmrescreate(&new_comm, solver, ierr);
}

void my_hypre_parasailscreate(int *comm, long int *solver, int *ierr)
{
   MPI_Comm new_comm;
   new_comm = hypre_fix_comm_f2c(comm);
   hypre_parasailscreate(&new_comm, solver, ierr);
}

/* Section on BoomerAMG_Set.  1.6.0, 1.9.0b and 2.0.0 need this section */

/* The array arguments in BoomerAMG_Set* are really wierd in hypre's F90* */

#define hypre_CTAlloc(type, count) \
( (type *)hypre_CAlloc((unsigned int)(count), (unsigned int)sizeof(type)) )

void my_hypre_boomeramgsetnumgridsweeps( long int *solver,
                                               int *num_grid_sweeps,
                                               int      *ierr             )
{
   int *sweeps;
   int i;
   sweeps = hypre_CTAlloc(int,4);
   for (i=0; i<4; i++) sweeps[i] = num_grid_sweeps[i];
   *ierr = (int) ( HYPRE_BoomerAMGSetNumGridSweeps(
/*                        (HYPRE_Solver) *solver, 
                        (int *)        num_grid_sweeps ) ); */
                                       *solver,
                                        sweeps) );
}

/* Section on PCGSetLogging.  1.9.0b needs this section */

/* Fortran interface for HYPRE_ParCSRPCGSetLogging is missing in 1.9.0b */

/*
void hypre_parcsrpcgsetlogging_(int *solver, int *logging, int *ierr)
{
   *ierr = (int) ( HYPRE_ParCSRPCGSetLogging( 
*/
                                          /*  (HYPRE_Solver) *solver,  */
/*
                                                             *solver,
                                              (int)          *logging ) );
}
*/
