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

module message_passing

!----------------------------------------------------
! This module is a dummy version of the message passing library
! to satisfy references in a sequential compilation
!----------------------------------------------------

!----------------------------------------------------
! Other modules used:

use global
use hash_mod
use hash_eq_mod
use stopwatch

implicit none
private
public proc_info, init_comm, terminate_comm, increase_universe, &
       decrease_universe, phaml_send, phaml_recv, phaml_alltoall, &
       phaml_global_max, phaml_global_min, phaml_global_sum, PARALLEL, &
       phaml_barrier, pause_until_enter, warning, fatal, &
       my_proc, num_proc, log_nproc, hostname, graphics_proc, &
       this_processors_procs, pack_procs_size, pack_procs, unpack_procs, &
       cleanup_unpack_procs, slaves_communicator, max_real_buffer, &
       sequential_send, sequential_recv, &
       GRAPHICS_TERMINATE, GRAPHICS_INIT, GRAPHICS_GRID, &
       GRAPHICS_CLOSE, UNKNOWN_SENDER, &
       NORMAL_SPAWN, DEBUG_SLAVE, DEBUG_GRAPHICS, DEBUG_BOTH, &
       HOSTLEN

!----------------------------------------------------
! The following types are defined:

integer, parameter :: HOSTLEN = 64
type proc_info
   private
   integer :: nothing ! gotta have something in there, I guess
   logical :: do_graphics
end type proc_info

type message_type
   integer, pointer :: imess(:)
   real(my_real), pointer :: rmess(:)
   integer :: ni, nr
end type message_type

type message_list
   type(message_type), pointer :: message(:)
end type message_list

!----------------------------------------------------
! The following parameters are defined:

integer, parameter :: PARALLEL=SEQUENTIAL
integer, parameter :: UNKNOWN_SENDER = -2
integer, parameter :: GRAPHICS_TERMINATE = 1, &
                               GRAPHICS_INIT      = 2, &
                               GRAPHICS_GRID      = 3, &
                               GRAPHICS_CLOSE     = 5
integer, parameter :: NORMAL_SPAWN =   1, &
                               DEBUG_SLAVE  =   2, &
                               DEBUG_GRAPHICS = 3, &
                               DEBUG_BOTH     = 4
!----------------------------------------------------
! The following variables are defined:

integer, save :: max_real_buffer = huge(0)
type(proc_info), pointer :: this_processors_procs
!----------------------------------------------------
! The following generic interfaces are defined:
 
interface phaml_global_max
   module procedure phaml_global_max_real, phaml_global_max_int
end interface
 
interface phaml_global_min
   module procedure phaml_global_min_real, phaml_global_min_int
end interface
 
interface phaml_global_sum
   module procedure phaml_global_sum_real, phaml_global_sum_int, &
                    phaml_global_sum_realv,phaml_global_sum_intv
end interface

interface phaml_alltoall
   module procedure phaml_alltoall_int, phaml_alltoall_real, &
                    phaml_alltoall_intv, phaml_alltoall_realv
end interface
!----------------------------------------------------

contains

subroutine init_comm(procs,spawn_form,i_draw_grid,master_draws_grid, &
                     graphics_host,output_unit,error_unit,system_size,eq_type, &
                     pde_id,nproc,max_blen,triangle_files,update_umod, &
                     debug_command,display)
type (proc_info), target :: procs
integer, intent(inout) :: spawn_form
logical, intent(inout) :: i_draw_grid, master_draws_grid
character(len=*), intent(inout) :: graphics_host
integer, intent(inout) :: output_unit, error_unit, pde_id, system_size, eq_type
integer, intent(in) :: nproc
real(my_real), intent(inout) :: max_blen
character(len=*), intent(inout) :: triangle_files
logical, intent(inout) :: update_umod
character(len=*), intent(in), optional :: debug_command, display
logical :: l
integer :: i
procs%nothing = nproc
procs%do_graphics = i_draw_grid .or. master_draws_grid
this_processors_procs => procs
update_umod = .false.
! just to avoid compiler warnings that the dummy arguments aren't used
l = i_draw_grid
l = master_draws_grid
i = output_unit
i = error_unit
i = pde_id
return
end subroutine init_comm

!          ---------------
subroutine sequential_send(imess,ni,rmess,nr)
!          ---------------

!----------------------------------------------------
! This routine sends a message to graphics when running a sequential program
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: imess(:),ni,nr
real(my_real), intent(in) :: rmess(:)
!----------------------------------------------------
! Local variables:

logical :: exists, opened, firstcall=.true.
integer :: iounit
!----------------------------------------------------
! Begin executable code

! on the first call, if the message files exist, delete them

if (firstcall) then
   iounit = 11
   do
      inquire(unit=iounit,exist=exists,opened=opened)
      if (exists .and. .not. opened) exit
      iounit = iounit + 1
   end do
   inquire(file=trim(scratch)//"phaml_message",exist=exists)
   if (exists) then
      open(iounit,file=trim(scratch)//"phaml_message")
      close(iounit,status="DELETE")
   endif
   inquire(file=trim(scratch)//"phaml_lock",exist=exists)
   if (exists) then
      open(iounit,file=trim(scratch)//"phaml_lock")
      close(iounit,status="DELETE")
   endif
   firstcall = .false.
endif

! if the message file exists, wait for the graphics program to delete it

exists = .true.
do while (exists)
   inquire(file=trim(scratch)//"phaml_message",exist=exists)
end do 

! find an available io unit

iounit = 11
do
   inquire(unit=iounit,exist=exists,opened=opened)
   if (exists .and. .not. opened) exit
   iounit = iounit + 1
end do

! create a lock file to indicate data is being written

open(unit=iounit,file=trim(scratch)//"phaml_lock")
write(iounit,*) 0
close(iounit)

! write the message

open(unit=iounit,file=trim(scratch)//"phaml_message")
write(iounit,*) ni, nr
write(iounit,*) imess
write(iounit,*) rmess
close(iounit)

! delete the lock file

open(unit=iounit,file=trim(scratch)//"phaml_lock")
close(unit=iounit,status="DELETE")

end subroutine sequential_send

!          ---------------
subroutine sequential_recv(imess,ni,rmess,nr)
!          ---------------

!----------------------------------------------------
! This routine receives a message from the program when running sequentially
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, pointer :: imess(:)
real(my_real), pointer :: rmess(:)
integer, intent(out) :: ni, nr
!----------------------------------------------------
! Local variables:

logical :: exists, opened
integer :: allocstat, iounit
!----------------------------------------------------
! Begin executable code

! ni==0 indicates no message yet

ni = 0

! see if the message file exists

inquire(file=trim(scratch)//"phaml_message",exist=exists)
if (.not. exists) return

! see if the message is complete

inquire(file=trim(scratch)//"phaml_lock",exist=exists)
if (exists) return

! find an available io unit

iounit = 11
do
   inquire(unit=iounit,exist=exists,opened=opened)
   if (exists .and. .not. opened) exit
   iounit = iounit + 1
end do

! read the message file

open(unit=iounit,file=trim(scratch)//"phaml_message")
read(iounit,*) ni, nr
allocate(imess(ni),rmess(nr),stat=allocstat)
if (allocstat /= 0) then
   call fatal("allocation failed for receiving message",intlist=(/allocstat/))
   stop
endif
read(iounit,*) imess
read(iounit,*) rmess

! delete message file

close(unit=iounit,status="DELETE")

end subroutine sequential_recv

subroutine terminate_comm(procs,finalize_mpi)
type (proc_info), intent(inout) :: procs
logical, intent(in) :: finalize_mpi
if (procs%do_graphics) then
   call sequential_send((/GRAPHICS_TERMINATE/),1,(/0.0_my_real/),0)
endif
return
end subroutine terminate_comm

subroutine increase_universe(procs,recv_int)
type (proc_info), intent(inout) :: procs
integer, intent(in) :: recv_int(:)
end subroutine increase_universe

subroutine decrease_universe(procs)
type (proc_info), intent(inout) :: procs
end subroutine decrease_universe

subroutine phaml_send(procs,proc,int_message,ni,real_message,nr,tag)
type (proc_info), intent(in) :: procs
integer, intent(in) :: proc, tag,ni,nr
integer, intent(in) :: int_message(ni)
real(my_real), intent(in) :: real_message(nr)
integer i,ii
i=procs%nothing
i=proc; i=tag; ii=size(int_message); ii=size(real_message)
return
end subroutine phaml_send

subroutine phaml_recv(procs,proc,int_message,ni,real_message,nr,tag,noblock)
type (proc_info), intent(in) :: procs
integer, intent(in) :: tag
integer, intent(out) :: proc,ni,nr
integer, pointer :: int_message(:)
real(my_real), pointer :: real_message(:)
logical, optional, intent(in) :: noblock
! just to avoid compiler warnings that the dummy arguments aren't used
proc=procs%nothing
if (present(noblock)) then
   proc=tag; proc=0; ni=0; nr=0
endif
nullify(int_message,real_message)
return
end subroutine phaml_recv

subroutine phaml_alltoall_int(procs,send,nsend,recv,nrecv,tag)
type (proc_info), intent(in) :: procs
integer, intent(in) :: send(:)
integer, intent(in) :: nsend
integer, pointer :: recv(:)
integer, intent(out) :: nrecv(:)
integer, intent(in) :: tag
integer :: idum
if (size(nrecv)>0) nrecv=0
! just to avoid compiler warnings that the dummy arguments aren't used
idum=procs%nothing
if (size(send)>0) idum=nsend
if (size(nrecv)>0) idum=tag
nullify(recv)
end subroutine phaml_alltoall_int

subroutine phaml_alltoall_real(procs,send,nsend,recv,nrecv,tag)
type (proc_info), intent(in) :: procs
real(my_real), intent(in) :: send(:)
integer, intent(in) :: nsend
real(my_real), pointer :: recv(:)
integer, intent(out) :: nrecv(:)
integer, intent(in) :: tag
integer :: idum
if (size(nrecv)>0) nrecv=0
! just to avoid compiler warnings that the dummy arguments aren't used
idum=procs%nothing
if (size(send)>0) idum=nsend
if (size(nrecv)>0) idum=tag
nullify(recv)
end subroutine phaml_alltoall_real

subroutine phaml_alltoall_intv(procs,send,nsend,recv,nrecv,tag)
type (proc_info), intent(in) :: procs
integer, intent(in) :: send(:)
integer, intent(in) :: nsend(:)
integer, pointer :: recv(:)
integer, intent(out) :: nrecv(:)
integer, intent(in) :: tag
integer :: idum
if (size(nrecv)>0) nrecv=0
! just to avoid compiler warnings that the dummy arguments aren't used
if (size(send)>0 .and. size(nsend)>0) idum=procs%nothing
if (size(nrecv)>0) idum=tag
nullify(recv)
end subroutine phaml_alltoall_intv

subroutine phaml_alltoall_realv(procs,send,nsend,recv,nrecv,tag)
type (proc_info), intent(in) :: procs
real(my_real), intent(in) :: send(:)
integer, intent(in) :: nsend(:)
real(my_real), pointer :: recv(:)
integer, intent(out) :: nrecv(:)
integer, intent(in) :: tag
integer :: idum
if (size(nrecv)>0) nrecv=0
! just to avoid compiler warnings that the dummy arguments aren't used
if (size(send)>0 .and. size(nsend)>0) idum=procs%nothing
if (size(nrecv)>0) idum=tag
nullify(recv)
end subroutine phaml_alltoall_realv

function phaml_global_max_real(procs,var,tag,just1)
type (proc_info), intent(in) :: procs
real (my_real) phaml_global_max_real
real (my_real), intent(in) :: var
integer, intent(in) :: tag
logical, optional, intent(in) :: just1
integer :: itag
phaml_global_max_real = var
! just to avoid compiler warnings that the dummy arguments aren't used
itag = procs%nothing
if (present(just1)) then
   itag = tag
endif
return
end function phaml_global_max_real

function phaml_global_max_int(procs,var,tag,just1)
type (proc_info), intent(in) :: procs
integer phaml_global_max_int
integer, intent(in) :: var
integer, intent(in) :: tag
logical, optional, intent(in) :: just1
integer :: itag
phaml_global_max_int = var
! just to avoid compiler warnings that the dummy arguments aren't used
itag = procs%nothing
if (present(just1)) then
   itag = tag
endif
return
end function phaml_global_max_int

function phaml_global_min_real(procs,var,tag,just1)
type (proc_info), intent(in) :: procs
real (my_real) phaml_global_min_real
real (my_real), intent(in) :: var
integer, intent(in) :: tag
logical, optional, intent(in) :: just1
integer :: itag
phaml_global_min_real = var
! just to avoid compiler warnings that the dummy arguments aren't used
itag = procs%nothing
if (present(just1)) then
   itag = tag
endif
return
end function phaml_global_min_real

function phaml_global_min_int(procs,var,tag,just1)
type (proc_info), intent(in) :: procs
integer phaml_global_min_int
integer, intent(in) :: var
integer, intent(in) :: tag
logical, optional, intent(in) :: just1
integer :: itag
phaml_global_min_int = var
! just to avoid compiler warnings that the dummy arguments aren't used
itag = procs%nothing
if (present(just1)) then
   itag = tag
endif
return
end function phaml_global_min_int

function phaml_global_sum_real(procs,var,tag,just1)
type (proc_info), intent(in) :: procs
real (my_real) phaml_global_sum_real
real (my_real), intent(in) :: var
integer, intent(in) :: tag
logical, optional, intent(in) :: just1
integer :: itag
phaml_global_sum_real = var
! just to avoid compiler warnings that the dummy arguments aren't used
itag = procs%nothing
if (present(just1)) then
   itag = tag
endif
return
end function phaml_global_sum_real

function phaml_global_sum_int(procs,var,tag,just1)
type (proc_info), intent(in) :: procs
integer phaml_global_sum_int
integer, intent(in) :: var
integer, intent(in) :: tag
logical, optional, intent(in) :: just1
integer :: itag
phaml_global_sum_int = var
! just to avoid compiler warnings that the dummy arguments aren't used
itag = procs%nothing
if (present(just1)) then
   itag = tag
endif
return
end function phaml_global_sum_int

function phaml_global_sum_realv(procs,var,tag,just1)
type (proc_info), intent(in) :: procs
real(my_real), intent(in) :: var(:)
integer, intent(in) :: tag
logical, optional, intent(in) :: just1
real, dimension(size(var)) :: phaml_global_sum_realv
integer :: itag
phaml_global_sum_realv = var
! just to avoid compiler warnings that the dummy arguments aren't used
itag = procs%nothing
if (present(just1)) then
   itag = tag
endif
return
end function phaml_global_sum_realv

function phaml_global_sum_intv(procs,var,tag,just1)
type (proc_info), intent(in) :: procs
integer, intent(in) :: var(:)
integer, intent(in) :: tag
logical, optional, intent(in) :: just1
integer, dimension(size(var)) :: phaml_global_sum_intv
integer :: itag
phaml_global_sum_intv = var
! just to avoid compiler warnings that the dummy arguments aren't used
itag = procs%nothing
if (present(just1)) then
   itag = tag
endif
return
end function phaml_global_sum_intv

function my_proc(procs)
integer :: my_proc
type (proc_info), intent(in) :: procs
my_proc = procs%nothing
my_proc = 1
end function my_proc

function num_proc(procs)
integer :: num_proc
type (proc_info), intent(in) :: procs
num_proc = procs%nothing
num_proc = 1
end function num_proc

function log_nproc(procs)
integer :: log_nproc
type (proc_info), intent(in) :: procs
log_nproc = procs%nothing
log_nproc = 0
end function log_nproc
 
function hostname(procs)
character(len=HOSTLEN) hostname
type (proc_info), intent(in) :: procs
integer :: junk
junk = procs%nothing
hostname = " "
end function hostname

function graphics_proc(procs)
integer :: graphics_proc
type (proc_info), intent(in) :: procs
graphics_proc = procs%nothing
end function graphics_proc

function slaves_communicator(procs)
integer :: slaves_communicator ! note default integer
type (proc_info), intent(in) :: procs
slaves_communicator = 0
end function slaves_communicator

subroutine pack_procs_size(procs,ni,nr)
type(proc_info), intent(in) :: procs
integer, intent(out) :: ni,nr
ni = 0
nr = 0
end subroutine pack_procs_size

subroutine pack_procs(procs,send_int,first_int,send_real,first_real)
type(proc_info), intent(in) :: procs
integer, intent(inout) :: send_int(:)
real(my_real), intent(inout) :: send_real(:)
integer, intent(inout) :: first_int, first_real
end subroutine pack_procs

subroutine unpack_procs(procs,send_int,first_int,send_real,first_real)
type(proc_info), intent(out) :: procs
integer, intent(in) :: send_int(:)
real(my_real), intent(in) :: send_real(:)
integer, intent(inout) :: first_int, first_real
procs%nothing=0
end subroutine unpack_procs

subroutine cleanup_unpack_procs(procs)
type(proc_info), intent(inout) :: procs
end subroutine cleanup_unpack_procs

subroutine phaml_barrier(procs)
type(proc_info), intent(in) :: procs
end subroutine phaml_barrier

!          -----------------
subroutine pause_until_enter(procs,really_pause,just_one,dont_pausewatch)
!          -----------------
 
!----------------------------------------------------
! This routine pauses until the return (or enter) key is pressed.  The
! master and all slaves must call this routine to avoid deadlock.
! If really_pause is false, the pause does not occur.
! If just_one is present and true, the 'pause complete' message is sent
! from the master to only slave 1.
!----------------------------------------------------
 
!----------------------------------------------------
! Dummy arguments
 
type(proc_info), intent(in) :: procs
logical, intent(in) :: really_pause
logical, intent(in), optional :: just_one, dont_pausewatch
 
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

if (really_pause) then
   if (present(dont_pausewatch)) then
      if (.not. dont_pausewatch) call pause_watch(all_watches)
   else
      call pause_watch(all_watches)
   endif
   print *,"press return to continue"
   read *
   if (present(dont_pausewatch)) then
      if (.not. dont_pausewatch) call end_pause_watch(all_watches)
   else
      call end_pause_watch(all_watches)
   endif
endif

return
end subroutine pause_until_enter

!          -------
subroutine warning(msg,msg2,intlist,reallist)
!          -------

!----------------------------------------------------
! This routine handles warning messages
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

character(len=*), intent(in) :: msg
character(len=*), intent(in), optional :: msg2
integer, intent(in), optional, dimension(:) :: intlist
real(my_real), intent(in), optional, dimension(:) :: reallist
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

if (warn) then
   write(errunit,"(A)")
   write(errunit,"(A)") "------------------------------------------------------"
   write(errunit,"(3A)") "          PHAML Version ",version_number," WARNING"
   write(errunit,"(A)") msg
   if (present(msg2)) write(errunit,"(A)") msg2
   if (present(intlist)) write(errunit,"(7I11)") intlist
   if (present(reallist)) write(errunit,"(SS,1P,4E18.10E2)") reallist
   write(errunit,"(A)") "------------------------------------------------------"
   write(errunit,"(A)")
endif

end subroutine warning

!          -----
subroutine fatal(msg,msg2,intlist,reallist,procs)
!          -----

!----------------------------------------------------
! This routine handles fatal errors
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

character(len=*), intent(in) :: msg
character(len=*), intent(in), optional :: msg2
integer, intent(in), optional, dimension(:) :: intlist
real(my_real), intent(in), optional, dimension(:) :: reallist
type(proc_info), intent(in), optional :: procs
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

write(errunit,"(A)")
write(errunit,"(A)") "------------------------------------------------------"
write(errunit,"(3A)") "          PHAML Version ",version_number," ERROR"
write(errunit,"(A)") msg
if (present(msg2)) write(errunit,"(A)") msg2
if (present(intlist)) write(errunit,"(7I11)") intlist
if (present(reallist)) write(errunit,"(SS,1P,4E18.10E2)") reallist
write(errunit,"(A)") "------------------------------------------------------"
write(errunit,"(A)")

stop
end subroutine fatal

end module message_passing

!----------------------------------------------------
! Dummy PVM routines to satisfy external references to PVM in the hash module

subroutine pvmfpack(type,ints,nint,tag,info)
use message_passing
integer, intent(in) :: type,ints(:),nint,tag
integer, intent(out) :: info
! assignments just to shut up picky compilers
call warning("Dummy pvmfpack routine called")
info = type
info = ints(1)
info = nint
info = tag
info = 1
end subroutine pvmfpack

subroutine pvmfunpack(type,ints,nint,tag,info)
use message_passing
integer, intent(in) :: type,tag
integer, intent(out) :: ints(:),nint,info
call warning("Dummy pvmfpack routine called")
! assignments just to shut up picky compilers
ints = type
nint = tag
info = 1
end subroutine pvmfunpack
