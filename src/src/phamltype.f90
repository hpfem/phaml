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

module phaml_type_mod

!----------------------------------------------------
! This module contains the definition of the primary type, phaml_solution_type.
! Although public, the internals should be considered to be private except
! in module phaml.
!----------------------------------------------------

use global
use message_passing
use gridtype_mod
use zoltan_interf

implicit none
private
public phaml_solution_type

!----------------------------------------------------
! Types defined are :

type phaml_solution_type
   type (grid_type) :: grid
   type (proc_info) :: procs
   integer :: outunit, errunit, pde_id, system_size, eq_type
   character(len=HOSTLEN) :: graphics_host
   type(Zoltan_Struct), pointer :: lb
   character(len=FN_LEN) :: zoltan_param_file
   logical :: i_draw_grid, master_draws_grid,still_sequential
end type phaml_solution_type

end module phaml_type_mod
