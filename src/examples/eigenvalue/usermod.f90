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
! This version illustrates eigenvalue problems via the Schroedinger Equation
! for a simple harmonic oscillator where the cross section of the well is an
! ellipse with axes determined by prob_param.
! -uxx -uyy + prob_param*(x^2+3*y^2)/2 * u = lambda * u
! prob_param must be positive.
!
! For num_eval=1 and lambda0=-huge(0.), the eigenvalue is
! sqrt(prob_param*(2+sqrt(3)))
! and the solution (wave function) is
! exp(-sqrt(prob_param/8)*(x^2+sqrt(3)*y^2).
!
!----------------------------------------------------

!----------------------------------------------------
! Other modules used are:

use phaml
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Global variables

real(my_real), save :: prob_param
!----------------------------------------------------

end module phaml_user_mod
