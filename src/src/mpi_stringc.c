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

/* MPI routines with string arguments, to avoid system dependencies */

#include <stdlib.h>
#include <string.h>
#include <mpi.h>

void mpi_myc_get_processor_name_(name, al, l, ierr)
int *name, *al, *l, *ierr;
{
  char *ch_name;
  int i, string_length;

  ch_name = (char *)malloc(*al*sizeof(char));
  *ierr = MPI_Get_processor_name(ch_name, l);
  string_length = strlen(ch_name);
  if (*al>string_length) *al=string_length;
  for (i=0; i<*al; i++) name[i] = (unsigned int)ch_name[i];
  free(ch_name);
}

/* alternate names */

void mpi_myc_get_processor_name(name, al, l, ierr)
int *name, *al, *l, *ierr;
{ mpi_myc_get_processor_name_(name, al, l, ierr); }

void mpi_myc_get_processor_name__(name, al, l, ierr)
int *name, *al, *l, *ierr;
{ mpi_myc_get_processor_name_(name, al, l, ierr); }

void MPI_MYC_GET_PROCESSOR_NAME(name, al, l, ierr)
int *name, *al, *l, *ierr;
{ mpi_myc_get_processor_name_(name, al, l, ierr); }
