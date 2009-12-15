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

/*
Contains C routines for the interface with BLOPEX compiled with PETSc.
*/

/* Fortran name mangling choices */

#define LOWERCASE 1
#define UPPERCASE 2
#define UNDERSCORE 3
#define DOUBLEUNDERSCORE 4

/* Select the form of name mangling used by your Fortran compiler.  This
   can either be set here by changing the next statement, or on the C
   compiler command line with something like -DFMANGLE=UNDERSCORE
*/

#define FMANGLE UNDERSCORE

/* mangle the names of routines callable from Fortran */

#if FMANGLE==LOWERCASE
#define petsc_lobpcg_solve_c   petsc_lobpcg_solve_c
#elif FMANGLE==UPPERCASE
#define petsc_lobpcg_solve_c   PETSC_LOBPCG_SOLVE_C
#elif FMANGLE==UNDERSCORE
#define petsc_lobpcg_solve_c   petsc_lobpcg_solve_c_
#elif FMANGLE==DOUBLEUNDERSCORE
#define petsc_lobpcg_solve_c   petsc_lobpcg_solve_c__
#endif

/* include files for PETSc and BLOPEX */

#include "petscvec.h"
#include "petscksp.h"
#include "petscda.h"
#include "lobpcg.h"
#include "petsc-interface.h"
#include "interpreter.h"
#include "multivector.h"

/* auxillary data structure */

typedef struct
{
  mv_InterfaceInterpreter  ii;
} aux_data_struct;

/* pointers to the Fortran callback functions.  The functions are:
   opA multiplies by the matrix A in the eigenvalue problem Ax = lambda x
       or Ax = lambda Bx
   opB multiplies by the matrix B in the generalized eigenvalue problem
       Ax = lambda Bx
   opT is a preconditioner, for example an approximation to A^(-1)
   return_evec passes the eigenvectors from C to Fortran, one at a time
   initial_guess sets an initial guess for the eigenvectors, one at a time,
       in Fortran and passes them to C
*/

void (*hold_matmult_opA)(void *, void*, void*);
void (*hold_matmult_opB)(void *, void*, void*);
void (*hold_matmult_opT)(void *, void*, void*);
void (*hold_petsc_lobpcg_return_evec)(void *);
void (*hold_petsc_lobpcg_initial_guess)(void *);

/* C callback functions that call the Fortran callback functions with
   a single vector */

void OperatorASingleVector (void *data, void *x, void *y)
{
   hold_matmult_opA(data,&x,&y);
}

void OperatorBSingleVector (void *data, void *x, void *y)
{
   hold_matmult_opB(data,&x,&y);
}

void OperatorTSingleVector (void *data, void *x, void *y)
{
   hold_matmult_opT(data,&x,&y);
}

void petsc_lobpcg_initial_guess_SingleVector (void *data, void *x, void *y)
{
   hold_petsc_lobpcg_initial_guess(&x);
}

void petsc_lobpcg_return_evec_SingleVector (void *data, void *x, void *y)
{
   hold_petsc_lobpcg_return_evec(&x);
}

/* C callback functions that take a multivector and call the Interface
   Interpreter routine that calls the single-vector-callback with each vector */

void OperatorAMultiVector(void * data, void * x, void * y)
{
   ((aux_data_struct*)data)->ii.Eval(OperatorASingleVector, data, x, y);
}

void OperatorBMultiVector(void * data, void * x, void * y)
{
   ((aux_data_struct*)data)->ii.Eval(OperatorBSingleVector, data, x, y);
}

void OperatorTMultiVector(void * data, void * x, void * y)
{
   ((aux_data_struct*)data)->ii.Eval(OperatorTSingleVector, data, x, y);
}

void petsc_lobpcg_initial_guess_MultiVector(void * data, void * x, void * y)
{
   ((aux_data_struct*)data)->ii.Eval(petsc_lobpcg_initial_guess_SingleVector, data, x, y);
}

void petsc_lobpcg_return_evec_MultiVector(void * data, void * x, void * y)
{
   ((aux_data_struct*)data)->ii.Eval(petsc_lobpcg_return_evec_SingleVector, data, x, y);
}

/* the main routine called from Fortran to solve the eigenvalue problem */

void petsc_lobpcg_solve_c(
   Vec* u,                           /* prototype of a vector, not used   */
   int* my_own_eq,                   /* equations owned by this processor */
   int* num_eval,                    /* number of eigenvalues to compute  */
   int* maxit,                       /* maximum number of iterations      */
   double* atol,                     /* absolute error tolerance          */
   double* rtol,                     /* relative error tolerance          */
   double* eigenvalues,              /* computed eigenvalues              */
   void *matmult_opA,                /* Fortran routine for operator A    */
   void *matmult_opB,                /* Fortran routine for operator B    */
   void *matmult_opT,                /* Fortran routine for operator T    */
   void *petsc_lobpcg_return_evec,   /* Fortran routine gets eigenvectors */
   void *petsc_lobpcg_initial_guess, /* Fortran routine for initial guess */
   int* info)                        /* error code                        */
{

   PetscErrorCode             ierr;         /* for PETSc return code        */
   mv_MultiVectorPtr          eigenvectors; /* the eigenvectors             */
   double *                   eigs;         /* the eigenvalues              */
   double *                   eigs_hist;    /* history of eigenvalues       */
   double *                   resid;        /* the residuals                */
   double *                   resid_hist;   /* history of residuals         */
   int                        iterations;   /* number of iterations         */
   int                        n_eigs;       /* number of eigenvalues        */
   int                        i,j;
   PetscTruth                 outpt=PETSC_FALSE; /* print evals and resids  */
   lobpcg_Tolerance           lobpcg_tol;   /* residual tolerance           */
   mv_InterfaceInterpreter    ii;           /* Interface Interpreter        */
   lobpcg_BLASLAPACKFunctions blap_fn;      /* BLAS functions               */
   aux_data_struct            aux_data;     /* auxillary data               */

/* set the number of eigenvalues to compute */

   n_eigs = *num_eval;

/* set pointers to the Fortran callback functions */

   hold_matmult_opA = matmult_opA;
   hold_matmult_opB = matmult_opB;
   hold_matmult_opT = matmult_opT;
   hold_petsc_lobpcg_initial_guess = petsc_lobpcg_initial_guess;
   hold_petsc_lobpcg_return_evec = petsc_lobpcg_return_evec;

/* allocate memory for the eigenvalues, residuals and histories */

   ierr = PetscMalloc(sizeof(double)*n_eigs,&eigs);
   ierr = PetscMalloc(sizeof(double)*n_eigs*(*maxit+1),&eigs_hist);
   ierr = PetscMalloc(sizeof(double)*n_eigs,&resid);
   ierr = PetscMalloc(sizeof(double)*n_eigs*(*maxit+1),&resid_hist);

/* create the Interface Interpreter and put it in auxillary data */

   PETSCSetupInterpreter( &ii );
   aux_data.ii = ii;

/* set tolerances and BLAS routines */

   lobpcg_tol.absolute = *atol;
   lobpcg_tol.relative = *rtol;

   blap_fn.dpotrf = PETSC_dpotrf_interface;
   blap_fn.dsygv = PETSC_dsygv_interface;

/* create the multivector for eigenvectors */

   eigenvectors = mv_MultiVectorCreateFromSampleVector(&ii, n_eigs,*u);

/* set the initial guess.  The second instance of eigenvectors in this
   call isn't actually used, but something has to be passed in */

   petsc_lobpcg_initial_guess_MultiVector(&aux_data,
                                          mv_MultiVectorGetData(eigenvectors),
                                          mv_MultiVectorGetData(eigenvectors));

/* call the lobpcg solver from BLOPEX */

   ierr = lobpcg_solve( eigenvectors,
                        &aux_data,
                        OperatorAMultiVector,
                        &aux_data,
                        OperatorBMultiVector,
                        &aux_data,
                        OperatorTMultiVector,
                        NULL,
                        blap_fn,
                        lobpcg_tol,
                        *maxit,
                        0, /* verbosity, use 2 for debugging */
                        &iterations,
                        eigs,
                        eigs_hist,
                        n_eigs,
                        resid,
                        resid_hist,
                        n_eigs
   );

/* set the return error code to lobpcg's error code */

   *info = ierr;

/* copy the eigenvalues to the return variable */

   for (i=0;i<n_eigs;i++) eigenvalues[i] = eigs[i];

/* return the eigenvectors.  The second instance of eigenvectors isn't used
   here either */

   petsc_lobpcg_return_evec_MultiVector(&aux_data,
                                        mv_MultiVectorGetData(eigenvectors),
                                        mv_MultiVectorGetData(eigenvectors));

/* printed output, for debugging */

  if (outpt)
  {
      PetscPrintf(PETSC_COMM_WORLD,"Output from LOBPCG driver\n");
      PetscPrintf(PETSC_COMM_WORLD,"   iterations: %d\n",iterations);
      PetscPrintf(PETSC_COMM_WORLD,"   eigenvalues and residuals:\n");
      for (i=0;i<n_eigs;i++)
        {
                ierr = PetscPrintf(PETSC_COMM_WORLD,"%e %e\n",eigs[i],resid[i]);
        }

/*
      PetscPrintf(PETSC_COMM_WORLD,"Output from LOBPCG, eigenvalues history:\n");
      for (j=0; j<iterations+1; j++)
         for (i=0;i<n_eigs;i++)
         {
            ierr = PetscPrintf(PETSC_COMM_WORLD,"%e\n",*(eigs_hist+j*n_eigs+i));
         }
      PetscPrintf(PETSC_COMM_WORLD,"Output from LOBPCG, residual norms:\n");
      for (i=0;i<n_eigs;i++)
        {
                ierr = PetscPrintf(PETSC_COMM_WORLD,"%e\n",resid[i]);
        }

      PetscPrintf(PETSC_COMM_WORLD,"Output from LOBPCG, residual norms history:\n");
      for (j=0; j<iterations+1; j++)
         for (i=0;i<n_eigs;i++)
         {
            ierr = PetscPrintf(PETSC_COMM_WORLD,"%e\n",*(resid_hist+j*n_eigs+i));
         }
*/

   }

/* free work space */

   mv_MultiVectorDestroy(eigenvectors);
   ierr = PetscFree(eigs);
   ierr = PetscFree(eigs_hist);
   ierr = PetscFree(resid);
   ierr = PetscFree(resid_hist);
}
