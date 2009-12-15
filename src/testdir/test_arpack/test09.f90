!---------------------------------------------------------------------!
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
!---------------------------------------------------------------------!

!----------------------------------------------------
! This file contains the main program supplied by the user.
! This version illustrates an eigenvalue problem, the Schroedinger
! Equation for a simple harmonic oscillator.
!----------------------------------------------------

!       ------------
program phaml_master
!       ------------

!----------------------------------------------------
! This is the main program.
!----------------------------------------------------

!----------------------------------------------------
! Modules used are:

use phaml
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Local variables:

integer :: nvert,num_eval,nproc
real(my_real) :: lambda0
type(phaml_solution_type) :: pde1
real(my_real) :: evals(3), mlsresid, alsresid, el2resid(3), errest, &
                 eval_errest(3)
integer :: aiter, anconv, anumop, anumopb, anumreo, ainfo

!----------------------------------------------------
! Begin executable code

nvert = 1000
lambda0 = -huge(0.0_my_real)
num_eval = 3
nproc=4
call phaml_create(pde1,nproc,eq_type=EIGENVALUE)
call phaml_solve_pde(pde1,                     &
                     max_vert=nvert,             &
                     sequential_vert=250,        &
                     num_eval=num_eval,          &
                     lambda0=lambda0,            &
                     print_error_when=FINAL,     &
                     print_error_who=MASTER,     &
                     print_errest_what=LINF_ERREST, &
                     mg_tol=MG_NO_TOL,           &
                     mg_cycles=10)
call phaml_query(pde1, &
                 eigenvalues=evals, &
                 eigenvalue_error_estimate=eval_errest, &
                 error_estimator=LOCAL_PROBLEM_H, &
                 max_linsys_resid=mlsresid, &
                 ave_linsys_resid=alsresid, &
                 eigen_l2_resid=el2resid, &
                 arpack_iter=aiter, &
                 arpack_nconv=anconv, &
                 arpack_numop=anumop, &
                 arpack_numopb=anumopb, &
                 arpack_numreo=anumreo, &
                 arpack_info=ainfo)
call phaml_query(pde1,linf_error_estimate=errest,eigen=2)
write(6,102) "eigenvalues ",evals
write(6,102) "error estimates ",eval_errest
write(6,101) "max linear system residual ",mlsresid
write(6,101) "ave linear system residual ",alsresid
write(6,102) "eigensystem residuals ",el2resid
write(6,103) "arpack iterations ",aiter
write(6,103) "arpack number converged ",anconv
write(6,103) "arpack number of OP*x ",anumop
write(6,103) "arpack number of B*x ",anumopb
write(6,103) "arpack number of reorthogonalizatons ",anumreo
write(6,103) "arpack error flag ",ainfo
write(6,101) "error estimate for evec 2 ",errest
101 format(A,SS,1P,E19.12E2)
102 format(A,SS,1P,3E19.12E2)
103 format(A,I11)
call phaml_destroy(pde1)

end program phaml_master
