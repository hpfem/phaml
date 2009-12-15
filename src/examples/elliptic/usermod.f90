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
! This file contains the user supplied module phaml_user_mod
!----------------------------------------------------

module phaml_user_mod

!----------------------------------------------------
! This module contains user global data.
!
! In this version, that consists of a selection parameter to determine
! which problem to solve, and a couple parameters for those problems.
! This version just solves some simple scalar linear elliptic PDEs.
!----------------------------------------------------

!----------------------------------------------------
! Other modules used are:

use phaml
!----------------------------------------------------

implicit none

!----------------------------------------------------
! The following parameters are defined:

! Set choice to determine which problem to solve

integer, parameter :: choice = 6

! The choices are:
!  1) Laplace, solution x+y
!  2) more complicated Laplace solution
!  3) Laplace, singularity at the origin, from the popular L domain problem
!  4) Laplace with 0 boundary conditions on 3 sides and a
!     quartic on the top -- solution unknown
!  5) Laplace with 0.0 on 3 sides, piecewise linear on the top with the point
!     at which it is 1.0 specified by prob_param -- soln unknown
!  6) Poisson, solution x**ipower + y**ipower
!  7) Poisson, solution is a sharp peak
!  8) Poisson, x**ipower + y**ipower with Neuman and mixed boundary conditions
!  9) nonconstant coefficients
! 10) Poisson, solution x(1-x)y(1-y)
! 11) Helmholtz with parameter -a**2 and b.c. sin(a*x) + sin(a*y), a=prob_param.
!     Multigrid will fail if prob_param is sufficiently large.
! 12) Poisson, solution sin(a*x) + sin(a*y), a=prob_param
! 13) Poisson, solution sin(a*r), a=prob_param, r=sqrt(x**2+y**2)
! 14) Helmholtz with parameter a=prob_param, solution x**ipower + y**ipower
! 15) Poisson, solution tanh(prob_param*(x**2+y**2-.5)) circular wave front
! 16) Helmholtz with same coefficient and solution as the eigenvalue example
! 17) First order terms, solution x**ipower + y**ipower
! 18) Cross derivative term, solution x**ipower * y**ipower

integer, save :: ipower = 10 ! degree of polynomial for problem 6,8,14
real(my_real), save :: prob_param = 0.7_my_real
!----------------------------------------------------

end module phaml_user_mod
