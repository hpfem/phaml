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
! This file contains a user defined module phaml_user_mod
!----------------------------------------------------

module phaml_user_mod

!----------------------------------------------------
! This module contains user global data.
!
! This data isn't actually used in this program, but is here just to
! illustrate how to change global data on the slaves from the master
! program main.f90 by using subroutine update_usermod.
!----------------------------------------------------

use phaml
implicit none

! Some global variables

integer, save :: iglobal1, iglobal2
real(my_real), save :: realvar

! power for solution

integer, save :: p = 5

end module phaml_user_mod

!          --------------
subroutine update_usermod(phaml_solution) 
!          --------------

!----------------------------------------------------
! This external subroutine updates the module variables on the slave processes
! by sending them from the master process.
!----------------------------------------------------

use phaml
use phaml_user_mod
   
!----------------------------------------------------
! Dummy arguments

type(phaml_solution_type), intent(in) :: phaml_solution

!----------------------------------------------------
! Local variables:

! Declare these arrays big enough to hold the variables to be sent

integer :: iparam(3)
real(my_real) :: rparam(1)

!----------------------------------------------------
! Begin executable code

! Copy the module variables into the arrays, putting integer variables
! into iparam and real variables into rparam.

   iparam(1) = iglobal1
   iparam(2) = iglobal2
   iparam(3) = p
   rparam(1) = realvar

! Call the routine that performs the actual exchange.

   call master_to_slaves(phaml_solution,iparam,rparam)

! Copy the arrays into the module variables, using the same correspondence
! between module variable and array index as was used above.
   
   iglobal1 = iparam(1)
   iglobal2 = iparam(2)
   p = iparam(3)
   realvar = rparam(1)
   
end subroutine update_usermod
