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

module zoltan_interf

!----------------------------------------------------
! This is a dummy version of zoltan_interf to satisfy references
!----------------------------------------------------

!----------------------------------------------------
! Other modules used:

use global
use message_passing
use hash_mod
use gridtype_mod

implicit none
private
public zoltan_create_lb, zoltan_destroy_lb, zoltan_init, zoltan_partition, &
       Zoltan_Struct
public zoltan_child_order

type Zoltan_Struct
   integer :: dummy
end type Zoltan_Struct

contains

!          ----------------
subroutine zoltan_create_lb(lb,procs)
!          ----------------

type(Zoltan_Struct), pointer :: lb
type(proc_info), intent(in), target :: procs
integer :: i

i = my_proc(procs) ! just to shut up the compiler warnings
!i = lb%dummy

return
end subroutine zoltan_create_lb

!          -----------------
subroutine zoltan_destroy_lb(lb,procs)
!          -----------------

type(Zoltan_Struct), pointer :: lb
type(proc_info), intent(in), target :: procs
integer :: i

i = my_proc(procs) ! just to shut up the compiler warnings
!i = lb%dummy

return
end subroutine zoltan_destroy_lb

!          -----------
subroutine zoltan_init(lb,method,zoltan_param_file,procs)
!          -----------

type(Zoltan_Struct), pointer :: lb
integer, intent(in) :: method
character(len=*), intent(in) :: zoltan_param_file
type(proc_info), intent(in) :: procs
integer :: i

! just to shut up compiler warnings
i = method
i = my_proc(procs)
i = ichar(zoltan_param_file(1:1))
!i = lb%dummy

return
end subroutine zoltan_init

!          ----------------
subroutine zoltan_partition(grid,procs,lb,partmeth, &
                            export_gid,export_part,first_call)
!          ----------------

type(grid_type), target, intent(in) :: grid
type(proc_info), target, intent(in) :: procs
type(Zoltan_Struct), pointer :: lb
integer, intent(in) :: partmeth
type(hash_key), pointer :: export_gid(:)
integer, pointer :: export_part(:)
logical, optional :: first_call
integer :: i

call warning("Zoltan called in zoltan_interf.dum")
! just to shut up the compiler warnings
if (first_call) then
   i = grid%nvert
   i = my_proc(procs)
   i = partmeth
endif

return
end subroutine zoltan_partition

subroutine zoltan_child_order(order,ierr,lb,grid,procs,still_sequential)
integer, intent(inout) :: order(:)
integer, intent(out) :: ierr
type(Zoltan_Struct), pointer :: lb
type (grid_type), intent(in) :: grid
type (proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential
ierr = order(1)
ierr = lb%dummy
ierr = 1 ! set to 1 to avoid changing order in draw_grid
end subroutine zoltan_child_order

end module zoltan_interf
