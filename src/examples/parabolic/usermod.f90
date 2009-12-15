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
! This file contains module phaml_user_mod and subroutine update_usermod
!----------------------------------------------------

module phaml_user_mod

!----------------------------------------------------
! This module contains user global data.
!
! This version is for a parabolic equation.  It contains the variable that
! holds the time step.
!----------------------------------------------------

!----------------------------------------------------
! Other modules used are:

use global
!----------------------------------------------------

implicit none

!----------------------------------------------------

! The time step

real(my_real) :: t, deltat

end module phaml_user_mod

!          --------------
subroutine update_usermod(phaml_solution) 
!          --------------

!----------------------------------------------------
! This routine updates the module variables on the slave processes by
! sending them from the master process
!----------------------------------------------------

use phaml
use phaml_user_mod
   
!----------------------------------------------------
! Dummy arguments

type(phaml_solution_type), intent(in) :: phaml_solution

!----------------------------------------------------
! Local variables:

! Declare these arrays big enough to hold the variables to be sent

integer :: iparam(1)
real(my_real) :: rparam(2)

!----------------------------------------------------
! Begin executable code

! Copy the module variables into the arrays, putting integer variables
! into iparam and real variables into rparam.

   rparam(1) = deltat
   rparam(2) = t

! Call the routine that performs the actual exchange.  Don't change this line.

   call master_to_slaves(phaml_solution,iparam,rparam)

! Copy the arrays into the module variables, using the same correspondence
! between module variable and array index as was used above.
   
   deltat = rparam(1)
   t = rparam(2)
   
end subroutine update_usermod
