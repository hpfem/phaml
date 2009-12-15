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
!--------------------------------------------------------------------*/

! Contains a dummy version of the C routines for the interface with BLOPEX

subroutine petsc_lobpcg_solve_c(u, my_own_eq, num_eval, maxit, atol, rtol, &
                                eigenvalues, matmult_opA, matmult_opB, &
                                matmult_opT, petsc_lobpcg_return_evec, &
                                petsc_lobpcg_initial_guess, info)
use message_passing
double precision u(*), atol, rtol, eigenvalues(*)
integer my_own_eq, num_eval, maxit, info
external matmult_opA, matmult_opB, matmult_opT, petsc_lobpcg_return_evec, &
         petsc_lobpcg_initial_guess
call fatal("dummy version of petsc_lobpcg_solve_c called")
stop
end subroutine petsc_lobpcg_solve_c
