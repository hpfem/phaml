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
! This module contains the phaml interface to a message passing library.
! This version is for PVM Release 3.4
!
! communication tags in this module are of the form 0xx
!----------------------------------------------------

!----------------------------------------------------
! Other modules used:

use global
use fpvm_mod
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
       GRAPHICS_CLOSE, &
       NORMAL_SPAWN, DEBUG_SLAVE, DEBUG_GRAPHICS, DEBUG_BOTH, &
       HOSTLEN

!----------------------------------------------------
! The following types are defined:

! NOTE changing this requires changing [un]pack_procs[_size]
integer, parameter :: HOSTLEN = 64
type proc_info
   private
   integer :: my_proc, nproc, log_nproc, graphics_proc
   character(len=HOSTLEN) :: hostname
   integer, pointer :: tids(:)
end type proc_info
! NOTE procs 1-n are the n slaves; proc 0 is the master; proc graphics_proc
! is my graphics (if it exists)

type message_type
   integer, pointer :: imess(:)
   real(my_real), pointer :: rmess(:)
   integer :: ni, nr
end type message_type

type message_list
   type(message_type), pointer :: message(:)
end type message_list

!----------------------------------------------------

!----------------------------------------------------
! The following parameters are defined:

integer, parameter :: PARALLEL           = PVM
integer, parameter :: UNKNOWN_SENDER     = -2
integer, parameter :: GRAPHICS_TERMINATE = 1, &
                      GRAPHICS_INIT      = 2, &
                      GRAPHICS_GRID      = 3, &
                      GRAPHICS_CLOSE     = 5
integer, parameter :: NORMAL_SPAWN       = 1, &
                      DEBUG_SLAVE        = 2, &
                      DEBUG_GRAPHICS     = 3, &
                      DEBUG_BOTH         = 4

logical, parameter :: print_comm = .false. ! for debugging and tuning
!----------------------------------------------------

!----------------------------------------------------
! The following variables are defined:

integer, save :: max_real_buffer = huge(0)
integer itype, rtype  ! the pvm type identifiers for integer and my_real
type(proc_info), pointer, save :: this_processors_procs

!----------------------------------------------------

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

interface

!NAS$ ALIEN "F77 pvmfmytid"
    subroutine pvmfmytid(tid)
    integer tid
    end subroutine pvmfmytid

!NAS$ ALIEN "F77 pvmfparent"
    subroutine pvmfparent(tid)
    integer tid
    end subroutine pvmfparent

!NAS$ ALIEN "F77 pvmfspawn"
    subroutine pvmfspawn(task,flag,where,ntask,tids,numt)
    character(len=*) :: task,where
    integer flag,ntask,tids(*),numt
    end subroutine pvmfspawn

!NAS$ ALIEN "F77 pvmfinitsend"
    subroutine pvmfinitsend(encoding,bufid)
    integer encoding,bufid
    end subroutine pvmfinitsend

!NAS$ ALIEN "F77 pvmfrecv"
    subroutine pvmfrecv(tid,msgtag,bufid)
    integer tid,msgtag,bufid
    end subroutine pvmfrecv

!NAS$ ALIEN "F77 pvmfnrecv"
    subroutine pvmfnrecv(tid,msgtag,bufid)
    integer tid,msgtag,bufid
    end subroutine pvmfnrecv
 
!NAS$ ALIEN "F77 pvmfexit"
   subroutine pvmfexit(info)
   integer info
   end subroutine pvmfexit

!NAS$ ALIEN "F77 pvmfsend"
   subroutine pvmfsend(tid,msgtag,info)
   integer tid,msgtag,info
   end subroutine pvmfsend

!NAS$ ALIEN "F77 pvmftidtohost"
   subroutine pvmftidtohost(mytid,mydtid)
   integer mytid,mydtid
   end subroutine pvmftidtohost

!NAS$ ALIEN "F77 pvmfconfig"
   subroutine pvmfconfig(nhost,narch,dtid,name,arch,speed,info)
   integer nhost,narch,dtid,speed,info
   character(len=*) :: name,arch
   end subroutine pvmfconfig

!NAS$ ALIEN "F77 pvmfbufinfo"
   subroutine pvmfbufinfo(bufid,bytes,msgtag,tid,info)
   integer bufid,bytes,msgtag,tid,info
   end subroutine pvmfbufinfo

end interface

!----------------------------------------------------

contains

!          ---------
subroutine init_comm(procs,spawn_form,i_draw_grid,master_draws_grid, &
                     graphics_host,output_unit,error_unit,system_size,eq_type, &
                     pde_id,nproc,max_blen,triangle_files,update_umod, &
                     debug_command,display)
!          ---------

!----------------------------------------------------
! This routine performs initialization for the communication.
! nproc should be provided for the master process which spawns others
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

type (proc_info), target :: procs
integer, intent(inout) :: spawn_form
logical, intent(inout) :: i_draw_grid, master_draws_grid
!(next two are intent(in) for master, intent(out) for slaves and graphics)
character(len=*), intent(inout) :: graphics_host
integer, intent(inout) :: output_unit, error_unit, pde_id, system_size, eq_type
integer, intent(in) :: nproc
real(my_real), intent(inout) :: max_blen
character(len=*), intent(inout) :: triangle_files
logical, intent(inout) :: update_umod
character(len=*), intent(in), optional :: debug_command, display
!----------------------------------------------------

!----------------------------------------------------
! Local variables:

integer :: from_proc,ni,nr,temp,proc,ngraph,i,j,info,my_type,allocstat
logical :: i_have_graphics
integer, pointer :: int_hostname(:)
integer, pointer :: irecv(:)
real(my_real), pointer :: rrecv(:)
real(my_real) :: rsend(1)
integer :: isend(14+HOSTLEN+nproc+1+len(triangle_files))
character(len=HOSTLEN) :: recv_name
character(len=HOSTLEN), allocatable :: from_hostname(:)
character(len=HOSTLEN) :: loc_graphics_host

integer :: mytid,parent_tid
integer nhost,narch,dtid,speed,mydtid
character(len=16) :: arch

!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

! set perminant pointer to procs

this_processors_procs => procs

! enroll in pvm

call pvmfmytid(mytid)

! identify the integer and real types

select case(bit_size(1)/8)
   case(2)
      itype = INTEGER2
   case(4)
      itype = INTEGER4
   case default
      call fatal("default integer kind does not match a pvm integer type")
      stop
end select

select case(kind(1.0_my_real))
   case(kind(1.0))
      rtype = REAL4
   case(kind(1.0d0))
      rtype = REAL8
   case default
      call fatal("kind my_real does not match a pvm real type")
      stop
end select

! different cases for master, slaves and graphics, determined by parent

call pvmfparent(parent_tid)

if (parent_tid < 0) then   ! I am the master process, start the slaves

   my_type = MASTER
   procs%my_proc = 0

! set the number of processors (slaves)

   procs%nproc = nproc
   procs%log_nproc = 0
   temp = procs%nproc
   do while (temp > 1)
      procs%log_nproc = procs%log_nproc + 1
      temp = temp/2
   end do
   allocate(procs%tids(0:procs%nproc+1),stat=allocstat)
   if (allocstat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in init_comm")
      return
   endif
   procs%tids = 0
   procs%tids(0) = mytid
   loc_graphics_host = trim(graphics_host)

! spawn the slaves

   select case(spawn_form)
   case (NORMAL_SPAWN,DEBUG_GRAPHICS)
      call pvmfspawn('phaml_slave',PVMDEFAULT,'*',procs%nproc, &
                     procs%tids(1:procs%nproc+1),info)
      if (info /= procs%nproc) call warning("pvmfspawn returned error code", &
                                            intlist=(/info/))
   case (DEBUG_SLAVE,DEBUG_BOTH)
      call pvmfspawn('phaml_slave',PVMTASKDEBUG,'*',procs%nproc, &
                     procs%tids(1:procs%nproc),info)
      if (info /= procs%nproc) call warning("pvmfspawn returned error code", &
                                            intlist=(/info/))
   case default
      call fatal("illegal value for spawn_form",intlist=(/spawn_form/))
      stop
   end select

! send slave code, nproc and graphics

   isend(1) = SLAVES
   isend(2) = procs%nproc
   if (i_draw_grid) then
      isend(3) = 1
   else
      isend(3) = 0
   endif
   if (master_draws_grid) then
      isend(4) = 1
   else
      isend(4) = 0
   endif
   isend(5) = 0   ! not used
   isend(6) = 0   ! not used
   isend(7) = pde_id
   isend(8) = spawn_form
   isend(9) = output_unit
   isend(10) = error_unit
   isend(11) = outunit
   isend(12) = system_size
   isend(13) = eq_type
   if (update_umod) then
      isend(14) = 1
   else
      isend(14) = 0
   endif
   isend(15:14+HOSTLEN) = (/ (ichar(loc_graphics_host(i:i)),i=1,HOSTLEN) /)
   isend(15+HOSTLEN:14+HOSTLEN+len(triangle_files)) = &
        (/ (ichar(triangle_files(i:i)),i=1,len(triangle_files)) /)
   isend(15+HOSTLEN+len(triangle_files):14+HOSTLEN+len(triangle_files)+procs%nproc+1) = &
         procs%tids(0:procs%nproc)
   rsend(1) = max_blen
   do proc=1,procs%nproc
      call phaml_send(procs,proc,isend, &
                      14+HOSTLEN+procs%nproc+1+len(triangle_files),rsend,1,10)
   end do

else ! I am either a slave or a graphics server

! receive the code that tells if I am slave or graphics, nproc and graphics

   procs%nproc = -2 ! prevents searching for from_proc before tids arrive
   call phaml_recv(procs,from_proc,irecv,ni,rrecv,nr,10)
   my_type = irecv(1)
   procs%nproc = irecv(2)
   i_draw_grid = irecv(3)==1
   master_draws_grid = irecv(4)==1
   pde_id = irecv(7)
   spawn_form = irecv(8)
   output_unit = irecv(9)
   error_unit = irecv(10)
   outunit = irecv(11)
   errunit = irecv(11)
   system_size = irecv(12)
   eq_type = irecv(13)
   update_umod = irecv(14)==1
   do i=1,min(HOSTLEN,len(graphics_host))
      graphics_host(i:i) = char(irecv(14+i))
   end do
   loc_graphics_host = trim(graphics_host)
   do i=1,len(triangle_files)
      triangle_files(i:i) = char(irecv(14+HOSTLEN+i))
   end do
   procs%log_nproc = 0
   temp = procs%nproc
   do while (temp > 1)
      procs%log_nproc = procs%log_nproc + 1
      temp = temp/2
   end do

   select case (my_type)

   case (SLAVES)
      allocate(procs%tids(0:procs%nproc+1),stat=allocstat)
      if (allocstat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in init_comm")
         return
      endif
      procs%tids(0:procs%nproc) = &
        irecv(15+HOSTLEN+len(triangle_files):14+HOSTLEN+len(triangle_files)+procs%nproc+1)
      procs%my_proc = -1
      do i=1,procs%nproc
         if (mytid == procs%tids(i)) then
            procs%my_proc = i
            exit
         endif
      end do
      if (procs%my_proc == -1) then
         call fatal("couldn't find my tid in the list of tids")
         stop
      endif

   case (GRAPHICS)
      allocate(procs%tids(0:0),stat=allocstat)
      if (allocstat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in init_comm")
         return
      endif
      procs%tids(0) = parent_tid
      procs%my_proc = -1
      procs%nproc = 0
      procs%graphics_proc = -1

   end select
   max_blen = rrecv(1)
   deallocate(irecv,rrecv,stat=allocstat)
endif

! if I'm the master or a slave, spawn graphics servers

if (my_type == MASTER .or. my_type == SLAVES) then

! determine number of servers and the index of my server

   procs%graphics_proc = -1
   i_have_graphics = .false.
   ngraph = 0
   if (master_draws_grid) then
      ngraph = ngraph + 1
      if (my_type == MASTER) i_have_graphics = .true.
   endif
   if (i_draw_grid) then
      ngraph = ngraph + procs%nproc
      if (my_type == SLAVES) i_have_graphics = .true.
   endif

   if (i_have_graphics) then
      if (trim(graphics_host) == "anywhere") then
         select case(spawn_form)
         case (NORMAL_SPAWN,DEBUG_SLAVE)
            call pvmfspawn('phaml_graphics',PVMDEFAULT,'*',1, &
                           procs%tids(procs%nproc+1:procs%nproc+1),info)
            if (info /= 1) call warning("pvmfspawn returned error code", &
                                        intlist=(/info/))
         case (DEBUG_GRAPHICS,DEBUG_BOTH)
            call pvmfspawn('phaml_graphics',PVMTASKDEBUG,'*',1, &
                           procs%tids(procs%nproc+1:procs%nproc+1),info)
            if (info /= 1) call warning("pvmfspawn returned error code", &
                                        intlist=(/info/))
         case default
            call fatal("illegal value for spawn_form",intlist=(/spawn_form/))
            stop
         end select
      else
         select case(spawn_form)
         case (NORMAL_SPAWN,DEBUG_SLAVE)
            call pvmfspawn('phaml_graphics',PVMTASKHOST,trim(graphics_host),&
                           1,procs%tids(procs%nproc+1:procs%nproc+1),info)
            if (info /= 1) call warning("pvmfspawn returned error code", &
                                        intlist=(/info/))
         case (DEBUG_GRAPHICS,DEBUG_BOTH)
            call pvmfspawn('phaml_graphics',PVMTASKDEBUG+PVMTASKHOST, &
                           trim(graphics_host),1, &
                           procs%tids(procs%nproc+1:procs%nproc+1),info)
            if (info /= 1) call warning("pvmfspawn returned error code", &
                                        intlist=(/info/))
         case default
            call fatal("illegal value for spawn_form",intlist=(/spawn_form/))
            stop
         end select
      endif
      procs%graphics_proc = procs%nproc + 1
   endif

! send graphics code, nproc and graphics

   if (procs%graphics_proc /= -1) then
      isend(1) = GRAPHICS
      isend(2) = procs%nproc
      if (my_type == SLAVES) then
         if (i_draw_grid) then
            isend(3) = 1
         else
            isend(3) = 0
         endif
      else
         if (master_draws_grid) then
            isend(3) = 1
         else
            isend(3) = 0
         endif
      endif
      if (master_draws_grid) then
         isend(4) = 1
      else
         isend(4) = 0
      endif
      isend(5) = 0    ! not used
      isend(6) = 0    ! not used
      isend(7) = pde_id
      isend(8) = spawn_form
      isend(9) = output_unit
      isend(10) = error_unit
      isend(11) = outunit
      isend(12) = system_size
      isend(13) = eq_type
      if (update_umod) then
         isend(14) = 1
      else
         isend(14) = 0
      endif
      isend(15:14+HOSTLEN) = (/ (ichar(loc_graphics_host(i:i)),i=1,HOSTLEN) /)
      isend(15+HOSTLEN:14+HOSTLEN+len(triangle_files)) = &
        (/ (ichar(triangle_files(i:i)),i=1,len(triangle_files)) /)
      rsend(1) = max_blen
      call phaml_send(procs,procs%graphics_proc,isend, &
                      14+HOSTLEN+len(triangle_files),rsend,1,10)
   endif

endif

! get the host name

call pvmftidtohost(mytid,mydtid)
procs%hostname = ''
do i=0,procs%nproc+1
   call pvmfconfig(nhost,narch,dtid,recv_name,arch,speed,info)
   if (mydtid == dtid) then
      procs%hostname=trim(recv_name)
      exit
   endif
end do
i = index(procs%hostname,".")
if (i /= 0) procs%hostname(i:) = " "

!write(outunit,"(A,I11,2A)") "I am processor number ",procs%my_proc," on host ",procs%hostname

! print host names from master

if (my_type == MASTER) then
   allocate(from_hostname(procs%nproc+1))
   do i=1,procs%nproc+1
      if (i==procs%nproc+1 .and. procs%tids(i)==0) exit
      call phaml_recv(procs,from_proc,int_hostname,ni,rrecv,nr,20)
      from_hostname(from_proc) = " "
      do j=1,ni
         from_hostname(from_proc)(j:j) = char(int_hostname(j))
      end do
      if (ni > 0) deallocate(int_hostname,stat=allocstat)
   end do
   do i=1,procs%nproc+1
!      write(outunit,"(A,I11,2A)") "Processor number ",i," is on host ",trim(from_hostname(i))
   end do
else
   allocate(int_hostname(len_trim(procs%hostname)),stat=allocstat)
   if (allocstat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in init_comm")
      return
   endif
   do i=1,len_trim(procs%hostname)
      int_hostname(i) = ichar(procs%hostname(i:i))
   end do
   call phaml_send(procs,MASTER,int_hostname,len_trim(procs%hostname), &
                   rsend,0,20)
   deallocate(int_hostname,stat=allocstat)
endif

return
end subroutine init_comm

!          --------------
subroutine terminate_comm(procs,finalize_mpi)
!          --------------

!----------------------------------------------------
! This routine terminates the communication package.
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

type (proc_info), intent(inout) :: procs
logical, intent(in) :: finalize_mpi
!----------------------------------------------------

!----------------------------------------------------
! Local variables:

integer :: i, allocstat
real (my_real) :: rmess(1)
integer :: imess(1)
integer :: ni, nr, from_proc, info
integer, pointer :: irecv(:)
real (my_real), pointer :: rrecv(:)
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

! have processor 0 wait for all others to complete

if (procs%my_proc == MASTER) then
   do i=1,procs%nproc
      call phaml_recv(procs,from_proc,irecv,ni,rrecv,nr,30)
   end do
   if (procs%graphics_proc /= -1) then
      imess(1)=GRAPHICS_TERMINATE
      call phaml_send(procs,procs%graphics_proc,imess,1,rmess,0,101)
   endif

! graphics engine; just quit

elseif (procs%my_proc == -1) then

   call pvmfexit(info)

! slave

else
   rmess(1)=float(procs%my_proc)
   call phaml_send(procs,MASTER,imess,0,rmess,1,30)
   if (procs%graphics_proc /= -1) then
      imess(1)=GRAPHICS_TERMINATE
      call phaml_send(procs,procs%graphics_proc,imess,1,rmess,0,101)
   endif
   call pvmfexit(info)
endif

! free up the space used for tids

deallocate(procs%tids,stat=allocstat)
nullify(this_processors_procs)

end subroutine terminate_comm

!          -----------------
subroutine increase_universe(procs,recv_int)
!          -----------------

!----------------------------------------------------
! Dummy version of a routine used by MPI 2
!----------------------------------------------------

type (proc_info), intent(inout) :: procs
integer, intent(in) :: recv_int(:)
end subroutine increase_universe

!          -----------------
subroutine decrease_universe(procs)
!          -----------------

!----------------------------------------------------
! Dummy version of a routine used by MPI 2
!----------------------------------------------------

type (proc_info), intent(inout) :: procs
end subroutine decrease_universe

!          ----------
subroutine phaml_send(procs,proc,int_message,ni,real_message,nr,tag)
!          ----------

!----------------------------------------------------
! This routine sends ni,nr,int_message,real_message to processor proc
! with a message tag tag.
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

type (proc_info), intent(in) :: procs
integer, intent(in) :: proc, tag,ni,nr
integer, intent(in) :: int_message(:)
real(my_real), intent(in) :: real_message(:)
!----------------------------------------------------

!----------------------------------------------------
! Local variables:

integer info, tag_default_kind, ni_default_kind, nr_default_kind
integer :: sizes(2)
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

if (print_comm) then
   write(outunit,"(3(A,I11))") "proc ",procs%my_proc," send to proc ",proc," tag ",tag
!   write(outunit,"(I11,A)")  ni," integers "
!   write(outunit,"(7I11))  int_message(1:ni)
!   write(outunit,"(I11,A)")  nr," reals "
!   write(outunit,"(4E18.10E2)")  real_message(1:nr)
endif

! initialize the buffer for sending

call pvmfinitsend(PVMDATADEFAULT,info)

! pack the data

sizes = (/ ni, nr /)
ni_default_kind = ni; nr_default_kind = nr

call pvmfpack_int(itype,sizes,2,1,info)
if (ni>0) call pvmfpack_int(itype,int_message,ni_default_kind,1,info)
if (nr>0) call pvmfpack_my_real(rtype,real_message,nr_default_kind,1,info)

! send it

tag_default_kind = tag
call pvmfsend(procs%tids(proc),tag_default_kind,info)

return
end subroutine phaml_send

!          ----------
subroutine phaml_recv(procs,proc,int_message,ni,real_message,nr,tag,noblock)
!          ----------

!----------------------------------------------------
! This routine receives a message with tag tag and unpacks ni,nr,
! int_message,real_message and proc
! If noblock is present and .true., then the receive is nonblocking.
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

type (proc_info), intent(in) :: procs
integer, intent(in) :: tag
integer, intent(out) :: proc,ni,nr
integer, pointer :: int_message(:)
real(my_real), pointer :: real_message(:)
logical, optional, intent(in) :: noblock
!----------------------------------------------------

!----------------------------------------------------
! Local variables:

integer info,bufid,bytes,from_tag,from_tid,i,tag_default_kind, &
        ni_default_kind,nr_default_kind,allocstat
integer :: sizes(2)
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

if (print_comm) then
   if (.not.present(noblock)) then
      write(outunit,"(2(A,I11))") "proc ",procs%my_proc," waiting to receive tag ",tag
   endif
endif

! receive the message

tag_default_kind = tag

! check for a message that some processor had a fatal error

call pvmfnrecv(-1,ERR_TAG,bufid)
if (bufid /= 0) then
   write(errunit,"(A)") "Some processor encountered a fatal error.  Halting."
endif

if (present(noblock)) then
   if (noblock) then
      call pvmfnrecv(-1,tag_default_kind,bufid)
   else
      call pvmfrecv(-1,tag_default_kind,bufid)
   endif
else
   call pvmfrecv(-1,tag_default_kind,bufid)
endif

! if bufid==0 it was an empty nonblocking receive

if (bufid==0) then
   proc=0; ni=0; nr=0
   return
endif

! identify the sender

call pvmfbufinfo(bufid,bytes,from_tag,from_tid,info)

proc=UNKNOWN_SENDER
do i=0,procs%nproc+1
   if (procs%tids(i) == from_tid) then
      proc = i
      exit
   endif
end do

! unpack the data, creating space along the way

call pvmfunpack_int(itype,sizes,2,1,info)
ni = sizes(1); nr=sizes(2)
ni_default_kind = ni; nr_default_kind = nr
if (ni>0) then
   allocate(int_message(ni),stat=allocstat)
   if (allocstat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in phaml_send")
      return
   endif
   call pvmfunpack_int(itype,int_message,ni_default_kind,1,info)
else
   nullify(int_message)
endif
if (nr>0) then
   allocate(real_message(nr),stat=allocstat)
   if (allocstat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in phaml_send")
      return
   endif
   call pvmfunpack_my_real(rtype,real_message,nr_default_kind,1,info)
else
   nullify(real_message)
endif

if (print_comm) then
   write(outunit,"(3(A,I11))") "proc ",procs%my_proc," received tag ",tag," from proc ",proc
!   write(outunit,"(I11,A)")  ni," integers "
!   write(outunit,"(7I11)")  int_message(1:ni)
!   write(outunit,"(I11,A)")  nr," reals "
!   write(outunit,"(4E18.10E2)")  real_message(1:nr)
endif

return
end subroutine phaml_recv

!          ------------------
subroutine phaml_alltoall_int(procs,send,nsend,recv,nrecv,tag)
!          ------------------

!----------------------------------------------------
! This routine performs all-to-all communication between the slaves.
! In this version, the integer array send, of size nsend, is sent to
! all slaves, including the sender.  The messages from all the processors are
! concatinated in recv, which is allocated in this routine and should be
! deallocated by the caller.  The messages are concatinated in order of
! processor number. nrecv, which should be passed in with size nproc, contains
! the size of each received message.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type (proc_info), intent(in) :: procs
integer, intent(in) :: send(:)
integer, intent(in) :: nsend
integer, pointer :: recv(:)
integer, intent(out) :: nrecv(:)
integer, intent(in) :: tag
!----------------------------------------------------
! Local variables:

integer :: proc, ni, nr, i
type arrays
   integer, pointer :: val(:)
end type arrays
type(arrays), allocatable :: temprecv(:)
real(my_real), pointer :: noreal(:)
!----------------------------------------------------
! Begin executable code

if (print_comm) then
   write(outunit,"(2(A,I11))") "proc ",procs%my_proc," in phaml_alltoall_int, tag ",tag
endif

! send the message to all processors except self

do proc=1,procs%nproc
   if (proc == procs%my_proc) cycle
   call phaml_send(procs,proc,send,nsend,(/0.0_my_real/),0,tag)
end do

! allocate place to keep received messages

allocate(temprecv(procs%nproc))

! receive messages.  I know I'm not getting one from myself, so use that
! pointer for the receive, and set the correct pointer after returning

do i=1,procs%nproc-1
   call phaml_recv(procs,proc,temprecv(procs%my_proc)%val,ni,noreal,nr,tag)
   nrecv(proc) = ni
   temprecv(proc)%val => temprecv(procs%my_proc)%val
end do
nrecv(procs%my_proc) = nsend

! allocate space for the return message

allocate(recv(sum(nrecv)))

! copy messages to return message

i=0
do proc=1,procs%nproc
   if (proc == procs%my_proc) then
      recv(i+1:i+nrecv(proc)) = send
   else
      recv(i+1:i+nrecv(proc)) = temprecv(proc)%val
   endif
   i = i + nrecv(proc)
end do

! free the individual messages

do proc=1,procs%nproc
   if (proc == procs%my_proc) cycle
   if (associated(temprecv(proc)%val)) deallocate(temprecv(proc)%val)
end do
deallocate(temprecv)

if (print_comm) then
   write(outunit,"(A,I11,A)") "proc ",procs%my_proc," done with phaml_alltoall_int"
endif

end subroutine phaml_alltoall_int

!          -------------------
subroutine phaml_alltoall_real(procs,send,nsend,recv,nrecv,tag)
!          -------------------

!----------------------------------------------------
! This routine performs all-to-all communication between the slaves.
! In this version, the real array send, of size nsend, is sent to
! all slaves, including the sender.  The messages from all the processors are
! concatinated in recv, which is allocated in this routine and should be
! deallocated by the caller.  The messages are concatinated in order of
! processor number. nrecv, which should be passed in with size nproc, contains
! the size of each received message.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type (proc_info), intent(in) :: procs
real(my_real), intent(in) :: send(:)
integer, intent(in) :: nsend
real(my_real), pointer :: recv(:)
integer, intent(out) :: nrecv(:)
integer, intent(in) :: tag
!----------------------------------------------------
! Local variables:

integer :: proc, ni, nr, i
type arrays
   real(my_real), pointer :: val(:)
end type arrays
type(arrays), allocatable :: temprecv(:)
integer, pointer :: noint(:)
!----------------------------------------------------
! Begin executable code

if (print_comm) then
   write(outunit,"(2(A,I11))") "proc ",procs%my_proc," in phaml_alltoall_real, tag ",tag
endif

! send the message to all processors except self

do proc=1,procs%nproc
   if (proc == procs%my_proc) cycle
   call phaml_send(procs,proc,(/0/),0,send,nsend,tag)
end do

! allocate place to keep received messages

allocate(temprecv(procs%nproc))

! receive messages.  I know I'm not getting one from myself, so use that
! pointer for the receive, and set the correct pointer after returning

do i=1,procs%nproc-1
   call phaml_recv(procs,proc,noint,ni,temprecv(procs%my_proc)%val,nr,tag)
   nrecv(proc) = nr
   temprecv(proc)%val => temprecv(procs%my_proc)%val
end do
nrecv(procs%my_proc) = nsend

! allocate space for the return message

allocate(recv(sum(nrecv)))

! copy messages to return message

i=0
do proc=1,procs%nproc
   if (proc == procs%my_proc) then
      recv(i+1:i+nrecv(proc)) = send
   else
      recv(i+1:i+nrecv(proc)) = temprecv(proc)%val
   endif
   i = i + nrecv(proc)
end do

! free the individual messages

do proc=1,procs%nproc
   if (proc == procs%my_proc) cycle
   if (associated(temprecv(proc)%val)) deallocate(temprecv(proc)%val)
end do
deallocate(temprecv)

if (print_comm) then
   write(outunit,"(A,I11,A)") "proc ",procs%my_proc," done with phaml_alltoall_real"
endif

end subroutine phaml_alltoall_real

!          -------------------
subroutine phaml_alltoall_intv(procs,send,nsend,recv,nrecv,tag)
!          -------------------

!----------------------------------------------------
! This routine performs all-to-all communication between the slaves.
! In this version, segments of the integer array send are sent to the
! other slaves.  nsend contains the size of each segment; the segments
! are contiguous and in processor-number order.  The messages from all
! processors are concatinated in recv, which is allocated in this routine
! and should be deallocated by the caller.  The messages are concatinated
! in order of processor number.  nrecv, which should be passed in with size
! nproc, contains the size of each received message.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type (proc_info), intent(in) :: procs
integer, intent(in) :: send(:)
integer, intent(in) :: nsend(:)
integer, pointer :: recv(:)
integer, intent(out) :: nrecv(:)
integer, intent(in) :: tag
!----------------------------------------------------
! Local variables:

integer :: proc, ni, nr, i
type arrays
   integer, pointer :: val(:)
end type arrays
type(arrays), allocatable :: temprecv(:)
real(my_real), pointer :: noreal(:)
!----------------------------------------------------
! Begin executable code

if (print_comm) then
   write(outunit,"(2(A,I11))") "proc ",procs%my_proc," in phaml_alltoall_intv, tag ",tag
endif

! allocate place to keep received messages

allocate(temprecv(procs%nproc))

! send the message to all processors except self

i=0
do proc=1,procs%nproc
   if (proc /= procs%my_proc) then
      call phaml_send(procs,proc,send(i+1:i+nsend(proc)),nsend(proc), &
                      (/0.0_my_real/),0,tag)
   endif
   i = i + nsend(proc)
end do

! receive messages.  I know I'm not getting one from myself, so use that
! pointer for the receive, and set the correct pointer after returning

do i=1,procs%nproc-1
   call phaml_recv(procs,proc,temprecv(procs%my_proc)%val,ni,noreal,nr,tag)
   nrecv(proc) = ni
   temprecv(proc)%val => temprecv(procs%my_proc)%val
end do
nrecv(procs%my_proc) = nsend(procs%my_proc)

! allocate space for the return message

allocate(recv(sum(nrecv)))

! copy messages to return message

i=0
do proc=1,procs%nproc
   if (proc == procs%my_proc) then
      recv(i+1:i+nrecv(proc)) = send(sum(nsend(1:procs%my_proc-1))+1:sum(nsend(1:procs%my_proc)))
   else
      recv(i+1:i+nrecv(proc)) = temprecv(proc)%val
   endif
   i = i + nrecv(proc)
end do

! free the individual messages

do proc=1,procs%nproc
   if (proc == procs%my_proc) cycle
   if (associated(temprecv(proc)%val)) deallocate(temprecv(proc)%val)
end do
deallocate(temprecv)

if (print_comm) then
   write(outunit,"(A,I11,A)") "proc ",procs%my_proc," done with phaml_alltoall_intv"
endif

end subroutine phaml_alltoall_intv

!          --------------------
subroutine phaml_alltoall_realv(procs,send,nsend,recv,nrecv,tag)
!          --------------------

!----------------------------------------------------
! This routine performs all-to-all communication between the slaves.
! In this version, segments of the real array send are sent to the
! other slaves.  nsend contains the size of each segment; the segments
! are contiguous and in processor-number order.  The messages from all
! processors are concatinated in recv, which is allocated in this routine
! and should be deallocated by the caller.  The messages are concatinated
! in order of processor number.  nrecv, which should be passed in with size
! nproc, contains the size of each received message.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type (proc_info), intent(in) :: procs
real(my_real), intent(in) :: send(:)
integer, intent(in) :: nsend(:)
real(my_real), pointer :: recv(:)
integer, intent(out) :: nrecv(:)
integer, intent(in) :: tag
!----------------------------------------------------
! Local variables:

integer :: proc, ni, nr, i
type arrays
   real(my_real), pointer :: val(:)
end type arrays
type(arrays), allocatable :: temprecv(:)
integer, pointer :: noint(:)
!----------------------------------------------------
! Begin executable code

if (print_comm) then
   write(outunit,"(2(A,I11))") "proc ",procs%my_proc," in phaml_alltoall_realv, tag ",tag
endif

! allocate place to keep received messages

allocate(temprecv(procs%nproc))

! send the message to all processors except self

i=0
do proc=1,procs%nproc
   if (proc /= procs%my_proc) then
      call phaml_send(procs,proc,(/0/),0,send(i+1:i+nsend(proc)),nsend(proc), &
                      tag)
   endif
   i = i + nsend(proc)
end do

! receive messages.  I know I'm not getting one from myself, so use that
! pointer for the receive, and set the correct pointer after returning

do i=1,procs%nproc-1
   call phaml_recv(procs,proc,noint,ni,temprecv(procs%my_proc)%val,nr,tag)
   nrecv(proc) = nr
   temprecv(proc)%val => temprecv(procs%my_proc)%val
end do
nrecv(procs%my_proc) = nsend(procs%my_proc)

! allocate space for the return message

allocate(recv(sum(nrecv)))

! copy messages to return message

i=0
do proc=1,procs%nproc
   if (proc == procs%my_proc) then
      recv(i+1:i+nrecv(proc)) = send(sum(nsend(1:procs%my_proc-1))+1:sum(nsend(1:procs%my_proc)))
   else
      recv(i+1:i+nrecv(proc)) = temprecv(proc)%val
   endif
   i = i + nrecv(proc)
end do

! free the individual messages

do proc=1,procs%nproc
   if (proc == procs%my_proc) cycle
   if (associated(temprecv(proc)%val)) deallocate(temprecv(proc)%val)
end do
deallocate(temprecv)

if (print_comm) then
   write(outunit,"(A,I11,A)") "proc ",procs%my_proc," done with phaml_alltoall_realv"
endif

end subroutine phaml_alltoall_realv

!        ---------------------
function phaml_global_max_real(procs,var,tag,just1)
!        ---------------------

!----------------------------------------------------
! This function performs a global reduction of a real variable.
! If just1 is present and .true., only processor 1 is given
! the result (other processors get the input returned); otherwise
! all processors get the result.
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

   real (my_real) phaml_global_max_real
   type (proc_info), intent(in) :: procs
   real (my_real), intent(in) :: var
   integer, intent(in) :: tag
   logical, optional, intent(in) :: just1
!----------------------------------------------------

!----------------------------------------------------
! Local variables:

integer i,ni,nr,p,allocstat
integer, pointer :: imess(:)
real (my_real), pointer :: rmess(:)
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

if (print_comm) then
   write(outunit,"(2(A,I11))") "proc ",procs%my_proc," in reduction operation, tag ",tag
endif

phaml_global_max_real = var

if (procs%my_proc == 1) then
! processor 1 receives messages from everyone and takes the max
   do i=2,procs%nproc
      call phaml_recv(procs,p,imess,ni,rmess,nr,tag)
      if (rmess(1) > phaml_global_max_real) phaml_global_max_real = rmess(1)
      deallocate(rmess,stat=allocstat)
   end do
! if .not.just1, processor 1 sends the result to the others
   if (.not.present(just1)) then
      allocate(imess(1),rmess(1),stat=allocstat)
      if (allocstat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in phaml_global_max_real")
         return
      endif
      rmess(1) = phaml_global_max_real
      do i=2,procs%nproc
         call phaml_send(procs,i,imess,0,rmess,1,tag)
      end do
      deallocate(imess,rmess,stat=allocstat)
   elseif (.not.just1) then
      allocate(imess(1),rmess(1),stat=allocstat)
      if (allocstat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in phaml_global_max_real")
         return
      endif
      rmess(1) = phaml_global_max_real
      do i=2,procs%nproc
         call phaml_send(procs,i,imess,0,rmess,1,tag)
      end do
      deallocate(imess,rmess,stat=allocstat)
   endif

else
! other processors send the value to processor 1
   allocate(imess(1),rmess(1),stat=allocstat)
   if (allocstat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in phaml_global_max_real")
      return
   endif
   rmess(1) = var
   call phaml_send(procs,1,imess,0,rmess,1,tag)
   deallocate(imess,rmess,stat=allocstat)
! if .not.just1, receive the result from processor 1
   if (.not.present(just1)) then
      call phaml_recv(procs,p,imess,ni,rmess,nr,tag)
      phaml_global_max_real = rmess(1)
      deallocate(rmess,stat=allocstat)
   elseif (.not.just1) then
      call phaml_recv(procs,p,imess,ni,rmess,nr,tag)
      phaml_global_max_real = rmess(1)
      deallocate(rmess,stat=allocstat)
   endif
endif

return
end function phaml_global_max_real

!        --------------------
function phaml_global_max_int(procs,var,tag,just1)
!        --------------------

!----------------------------------------------------
! This function performs a global reduction of an integer variable.
! If just1 is present and .true., only processor 1 is given
! the result (other processors get the input returned); otherwise
! all processors get the result.
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

   integer phaml_global_max_int
   type (proc_info), intent(in) :: procs
   integer, intent(in) :: var
   integer, intent(in) :: tag
   logical, optional, intent(in) :: just1
!----------------------------------------------------

!----------------------------------------------------
! Local variables:

integer i,ni,nr,p,allocstat
integer, pointer :: imess(:)
real (my_real), pointer :: rmess(:)
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

if (print_comm) then
   write(outunit,"(2(A,I11))") "proc ",procs%my_proc," in reduction operation, tag ",tag
endif

phaml_global_max_int = var

if (procs%my_proc == 1) then
! processor 1 receives messages from everyone and takes the max
   do i=2,procs%nproc
      call phaml_recv(procs,p,imess,ni,rmess,nr,tag)
      if (imess(1) > phaml_global_max_int) phaml_global_max_int = imess(1)
      deallocate(imess,stat=allocstat)
   end do
! if .not.just1, processor 1 sends the result to the others
   if (.not.present(just1)) then
      allocate(imess(1),rmess(1),stat=allocstat)
      if (allocstat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in phaml_global_max_int")
         return
      endif
      imess(1) = phaml_global_max_int
      do i=2,procs%nproc
         call phaml_send(procs,i,imess,1,rmess,0,tag)
      end do
      deallocate(imess,rmess,stat=allocstat)
   elseif (.not.just1) then
      allocate(imess(1),rmess(1),stat=allocstat)
      if (allocstat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in phaml_global_max_int")
         return
      endif
      imess(1) = phaml_global_max_int
      do i=2,procs%nproc
         call phaml_send(procs,i,imess,1,rmess,0,tag)
      end do
      deallocate(imess,rmess,stat=allocstat)
   endif

else
! other processors send the value to processor 1
   allocate(imess(1),rmess(1),stat=allocstat)
   if (allocstat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in phaml_global_max_int")
      return
   endif
   imess(1) = var
   call phaml_send(procs,1,imess,1,rmess,0,tag)
   deallocate(imess,rmess,stat=allocstat)
! if .not.just1, receive the result from processor 1
   if (.not.present(just1)) then
      call phaml_recv(procs,p,imess,ni,rmess,nr,tag)
      phaml_global_max_int = imess(1)
      deallocate(imess,stat=allocstat)
   elseif (.not.just1) then
      call phaml_recv(procs,p,imess,ni,rmess,nr,tag)
      phaml_global_max_int = imess(1)
      deallocate(imess,stat=allocstat)
   endif
endif

return
end function phaml_global_max_int

!        ---------------------
function phaml_global_min_real(procs,var,tag,just1)
!        ---------------------

!----------------------------------------------------
! This function performs a global reduction of a real variable.
! If just1 is present and .true., only processor 1 is given
! the result (other processors get the input returned); otherwise
! all processors get the result.
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

   real (my_real) phaml_global_min_real
   type (proc_info), intent(in) :: procs
   real (my_real), intent(in) :: var
   integer, intent(in) :: tag
   logical, optional, intent(in) :: just1
!----------------------------------------------------

!----------------------------------------------------
! Local variables:

integer i,ni,nr,p,allocstat
integer, pointer :: imess(:)
real (my_real), pointer :: rmess(:)
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

if (print_comm) then
   write(outunit,"(2(A,I11))") "proc ",procs%my_proc," in reduction operation, tag ",tag
endif

phaml_global_min_real = var

if (procs%my_proc == 1) then
! processor 1 receives messages from everyone and takes the min
   do i=2,procs%nproc
      call phaml_recv(procs,p,imess,ni,rmess,nr,tag)
      if (rmess(1) < phaml_global_min_real) phaml_global_min_real = rmess(1)
      deallocate(rmess,stat=allocstat)
   end do
! if .not.just1, processor 1 sends the result to the others
   if (.not.present(just1)) then
      allocate(imess(1),rmess(1),stat=allocstat)
      if (allocstat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in phaml_global_min_real")
         return
      endif
      rmess(1) = phaml_global_min_real
      do i=2,procs%nproc
         call phaml_send(procs,i,imess,0,rmess,1,tag)
      end do
      deallocate(imess,rmess,stat=allocstat)
   elseif (.not.just1) then
      allocate(imess(1),rmess(1),stat=allocstat)
      if (allocstat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in phaml_global_min_real")
         return
      endif
      rmess(1) = phaml_global_min_real
      do i=2,procs%nproc
         call phaml_send(procs,i,imess,0,rmess,1,tag)
      end do
      deallocate(imess,rmess,stat=allocstat)
   endif

else
! other processors send the value to processor 1
   allocate(imess(1),rmess(1),stat=allocstat)
   if (allocstat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in phaml_global_min_real")
      return
   endif
   rmess(1) = var
   call phaml_send(procs,1,imess,0,rmess,1,tag)
   deallocate(imess,rmess,stat=allocstat)
! if .not.just1, receive the result from processor 1
   if (.not.present(just1)) then
      call phaml_recv(procs,p,imess,ni,rmess,nr,tag)
      phaml_global_min_real = rmess(1)
      deallocate(rmess,stat=allocstat)
   elseif (.not.just1) then
      call phaml_recv(procs,p,imess,ni,rmess,nr,tag)
      phaml_global_min_real = rmess(1)
      deallocate(rmess,stat=allocstat)
   endif
endif

return
end function phaml_global_min_real

!        --------------------
function phaml_global_min_int(procs,var,tag,just1)
!        --------------------

!----------------------------------------------------
! This function performs a global reduction of an integer variable.
! If just1 is present and .true., only processor 1 is given
! the result (other processors get the input returned); otherwise
! all processors get the result.
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

   integer phaml_global_min_int
   type (proc_info), intent(in) :: procs
   integer, intent(in) :: var
   integer, intent(in) :: tag
   logical, optional, intent(in) :: just1
!----------------------------------------------------

!----------------------------------------------------
! Local variables:

integer i,ni,nr,p,allocstat
integer, pointer :: imess(:)
real (my_real), pointer :: rmess(:)
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

if (print_comm) then
   write(outunit,"(2(A,I11))") "proc ",procs%my_proc," in reduction operation, tag ",tag
endif

phaml_global_min_int = var

if (procs%my_proc == 1) then
! processor 1 receives messages from everyone and takes the min
   do i=2,procs%nproc
      call phaml_recv(procs,p,imess,ni,rmess,nr,tag)
      if (imess(1) < phaml_global_min_int) phaml_global_min_int = imess(1)
      deallocate(imess,stat=allocstat)
   end do
! if .not.just1, processor 1 sends the result to the others
   if (.not.present(just1)) then
      allocate(imess(1),rmess(1),stat=allocstat)
      if (allocstat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in phaml_global_min_int")
         return
      endif
      imess(1) = phaml_global_min_int
      do i=2,procs%nproc
         call phaml_send(procs,i,imess,1,rmess,0,tag)
      end do
      deallocate(imess,rmess,stat=allocstat)
   elseif (.not.just1) then
      allocate(imess(1),rmess(1),stat=allocstat)
      if (allocstat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in phaml_global_min_int")
         return
      endif
      imess(1) = phaml_global_min_int
      do i=2,procs%nproc
         call phaml_send(procs,i,imess,1,rmess,0,tag)
      end do
      deallocate(imess,rmess,stat=allocstat)
   endif

else
! other processors send the value to processor 1
   allocate(imess(1),rmess(1),stat=allocstat)
   if (allocstat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in phaml_global_min_int")
      return
   endif
   imess(1) = var
   call phaml_send(procs,1,imess,1,rmess,0,tag)
   deallocate(imess,rmess,stat=allocstat)
! if .not.just1, receive the result from processor 1
   if (.not.present(just1)) then
      call phaml_recv(procs,p,imess,ni,rmess,nr,tag)
      phaml_global_min_int = imess(1)
      deallocate(imess,stat=allocstat)
   elseif (.not.just1) then
      call phaml_recv(procs,p,imess,ni,rmess,nr,tag)
      phaml_global_min_int = imess(1)
      deallocate(imess,stat=allocstat)
   endif
endif

return
end function phaml_global_min_int

!        ---------------------
function phaml_global_sum_real(procs,var,tag,just1)
!        ---------------------

!----------------------------------------------------
! This function performs a global reduction of the variable.
! If just1 is present and .true., only processor 1 is given
! the result (other processors get the input returned); otherwise
! all processors get the result.
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

   real (my_real) phaml_global_sum_real
   type (proc_info), intent(in) :: procs
   real (my_real), intent(in) :: var
   integer, intent(in) :: tag
   logical, optional, intent(in) :: just1
!----------------------------------------------------

!----------------------------------------------------
! Local variables:

integer i,ni,nr,p,allocstat
integer, pointer :: imess(:)
real (my_real), pointer :: rmess(:)
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

if (print_comm) then
   write(outunit,"(2(A,I11))") "proc ",procs%my_proc," in reduction operation, tag ",tag
endif

phaml_global_sum_real = var

if (procs%my_proc == 1) then
! processor 1 receives messages from everyone and sums
   do i=2,procs%nproc
      call phaml_recv(procs,p,imess,ni,rmess,nr,tag)
      phaml_global_sum_real = phaml_global_sum_real + rmess(1)
      deallocate(rmess,stat=allocstat)
   end do
! if .not.just1, processor 1 sends the result to the others
   if (.not.present(just1)) then
      allocate(imess(1),rmess(1),stat=allocstat)
      if (allocstat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in phaml_global_sum_real")
         return
      endif
      rmess(1) = phaml_global_sum_real
      do i=2,procs%nproc
         call phaml_send(procs,i,imess,0,rmess,1,tag)
      end do
      deallocate(imess,rmess,stat=allocstat)
   elseif (.not.just1) then
      allocate(imess(1),rmess(1),stat=allocstat)
      if (allocstat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in phaml_global_sum_real")
         return
      endif
      rmess(1) = phaml_global_sum_real
      do i=2,procs%nproc
         call phaml_send(procs,i,imess,0,rmess,1,tag)
      end do
      deallocate(imess,rmess,stat=allocstat)
   endif

else
! other processors send the value to processor 1
   allocate(imess(1),rmess(1),stat=allocstat)
   if (allocstat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in phaml_global_sum_real")
      return
   endif
   rmess(1) = var
   call phaml_send(procs,1,imess,0,rmess,1,tag)
   deallocate(imess,rmess,stat=allocstat)
! if .not.just1, receive the result from processor 1
   if (.not.present(just1)) then
      call phaml_recv(procs,p,imess,ni,rmess,nr,tag)
      phaml_global_sum_real = rmess(1)
      deallocate(rmess,stat=allocstat)
   elseif (.not.just1) then
      call phaml_recv(procs,p,imess,ni,rmess,nr,tag)
      phaml_global_sum_real = rmess(1)
      deallocate(rmess,stat=allocstat)
   endif
endif

return
end function phaml_global_sum_real

!        --------------------
function phaml_global_sum_int(procs,var,tag,just1)
!        --------------------

!----------------------------------------------------
! This function performs a global reduction of the variable.
! If just1 is present and .true., only processor 1 is given
! the result (other processors get the input returned); otherwise
! all processors get the result.
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

   integer phaml_global_sum_int
   type (proc_info), intent(in) :: procs
   integer, intent(in) :: var
   integer, intent(in) :: tag
   logical, optional, intent(in) :: just1
!----------------------------------------------------

!----------------------------------------------------
! Local variables:

integer i,ni,nr,p,allocstat
integer, pointer :: imess(:)
real (my_real), pointer :: rmess(:)
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

if (print_comm) then
   write(outunit,"(2(A,I11))") "proc ",procs%my_proc," in reduction operation, tag ",tag
endif

phaml_global_sum_int = var

if (procs%my_proc == 1) then
! processor 1 receives messages from everyone and sums
   do i=2,procs%nproc
      call phaml_recv(procs,p,imess,ni,rmess,nr,tag)
      phaml_global_sum_int = phaml_global_sum_int + imess(1)
      deallocate(imess,stat=allocstat)
   end do
! if .not.just1, processor 1 sends the result to the others
   if (.not.present(just1)) then
      allocate(imess(1),rmess(1),stat=allocstat)
      if (allocstat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in phaml_global_sum_int")
         return
      endif
      imess(1) = phaml_global_sum_int
      do i=2,procs%nproc
         call phaml_send(procs,i,imess,1,rmess,0,tag)
      end do
      deallocate(imess,rmess,stat=allocstat)
   elseif (.not.just1) then
      allocate(imess(1),rmess(1),stat=allocstat)
      if (allocstat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in phaml_global_sum_int")
         return
      endif
      imess(1) = phaml_global_sum_int
      do i=2,procs%nproc
         call phaml_send(procs,i,imess,1,rmess,0,tag)
      end do
      deallocate(imess,rmess,stat=allocstat)
   endif

else
! other processors send the value to processor 1
   allocate(imess(1),rmess(1),stat=allocstat)
   if (allocstat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in phaml_global_sum_int")
      return
   endif
   imess(1) = var
   call phaml_send(procs,1,imess,1,rmess,0,tag)
   deallocate(imess,rmess,stat=allocstat)
! if .not.just1, receive the result from processor 1
   if (.not.present(just1)) then
      call phaml_recv(procs,p,imess,ni,rmess,nr,tag)
      phaml_global_sum_int = imess(1)
      deallocate(imess,stat=allocstat)
   elseif (.not.just1) then
      call phaml_recv(procs,p,imess,ni,rmess,nr,tag)
      phaml_global_sum_int = imess(1)
      deallocate(imess,stat=allocstat)
   endif
endif

return
end function phaml_global_sum_int

!        ----------------------
function phaml_global_sum_realv(procs,var,tag,just1)
!        ----------------------

!----------------------------------------------------
! This function performs a global reduction of the array of variables.
! If just1 is present and .true., only processor 1 is given
! the result (other processors get the input returned); otherwise
! all processors get the result.
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

   type (proc_info), intent(in) :: procs
   real(my_real), intent(in) :: var(:)
   integer, intent(in) :: tag
   logical, optional, intent(in) :: just1
   real(my_real), dimension(size(var)) :: phaml_global_sum_realv
!----------------------------------------------------

!----------------------------------------------------
! Local variables:

integer i,ni,nr,p,allocstat
integer, pointer :: imess(:)
real (my_real), pointer :: rmess(:)
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

if (print_comm) then
   write(outunit,"(2(A,I11))") "proc ",procs%my_proc," in reduction operation, tag ",tag
endif

phaml_global_sum_realv = var

if (procs%my_proc == 1) then
! processor 1 receives messages from everyone and sums
   do i=2,procs%nproc
      call phaml_recv(procs,p,imess,ni,rmess,nr,tag)
      phaml_global_sum_realv = phaml_global_sum_realv + rmess(1:nr)
      deallocate(rmess,stat=allocstat)
   end do
! if .not.just1, processor 1 sends the result to the others
   if (.not.present(just1)) then
      allocate(rmess(size(var)),imess(1),stat=allocstat)
      if (allocstat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in phaml_global_sum_realv")
         return
      endif
      rmess = phaml_global_sum_realv
      do i=2,procs%nproc
         call phaml_send(procs,i,imess,0,rmess,size(var),tag)
      end do
      deallocate(imess,rmess,stat=allocstat)
   elseif (.not.just1) then
      allocate(rmess(size(var)),imess(1),stat=allocstat)
      if (allocstat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in phaml_global_sum_realv")
         return
      endif
      rmess = phaml_global_sum_realv
      do i=2,procs%nproc
         call phaml_send(procs,i,imess,0,rmess,size(var),tag)
      end do
      deallocate(imess,rmess,stat=allocstat)
   endif

else
! other processors send the value to processor 1
   allocate(rmess(size(var)),imess(1),stat=allocstat)
   if (allocstat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in phaml_global_sum_realv")
      return
   endif
   rmess = var
   call phaml_send(procs,1,imess,0,rmess,size(var),tag)
   deallocate(imess,rmess,stat=allocstat)
! if .not.just1, receive the result from processor 1
   if (.not.present(just1)) then
      call phaml_recv(procs,p,imess,ni,rmess,nr,tag)
      phaml_global_sum_realv = rmess(1:nr)
      deallocate(rmess,stat=allocstat)
   elseif (.not.just1) then
      call phaml_recv(procs,p,imess,ni,rmess,nr,tag)
      phaml_global_sum_realv = rmess(1:nr)
      deallocate(rmess,stat=allocstat)
   endif
endif

return
end function phaml_global_sum_realv

!        ---------------------
function phaml_global_sum_intv(procs,var,tag,just1)
!        ---------------------

!----------------------------------------------------
! This function performs a global reduction of the array of variables.
! If just1 is present and .true., only processor 1 is given
! the result (other processors get the input returned); otherwise
! all processors get the result.
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

   type (proc_info), intent(in) :: procs
   integer, intent(in) :: var(:)
   integer, intent(in) :: tag
   logical, optional, intent(in) :: just1
   integer, dimension(size(var)) :: phaml_global_sum_intv
!----------------------------------------------------

!----------------------------------------------------
! Local variables:

integer i,ni,nr,p,allocstat
integer, pointer :: imess(:)
real (my_real), pointer :: rmess(:)
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

if (print_comm) then
   write(outunit,"(2(A,I11))") "proc ",procs%my_proc," in reduction operation, tag ",tag
endif

phaml_global_sum_intv = var

if (procs%my_proc == 1) then
! processor 1 receives messages from everyone and sums
   do i=2,procs%nproc
      call phaml_recv(procs,p,imess,ni,rmess,nr,tag)
      phaml_global_sum_intv = phaml_global_sum_intv + imess(1:ni)
      deallocate(imess,stat=allocstat)
   end do
! if .not.just1, processor 1 sends the result to the others
   if (.not.present(just1)) then
      allocate(imess(size(var)),rmess(1),stat=allocstat)
      if (allocstat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in phaml_global_sum_intv")
         return
      endif
      imess = phaml_global_sum_intv
      do i=2,procs%nproc
         call phaml_send(procs,i,imess,size(var),rmess,0,tag)
      end do
      deallocate(imess,rmess,stat=allocstat)
   elseif (.not.just1) then
      allocate(imess(size(var)),rmess(1),stat=allocstat)
      if (allocstat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in phaml_global_sum_intv")
         return
      endif
      imess = phaml_global_sum_intv
      do i=2,procs%nproc
         call phaml_send(procs,i,imess,size(var),rmess,0,tag)
      end do
      deallocate(imess,rmess,stat=allocstat)
   endif

else
! other processors send the value to processor 1
   allocate(imess(size(var)),rmess(1),stat=allocstat)
   if (allocstat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in phaml_global_sum_intv")
      return
   endif
   imess = var
   call phaml_send(procs,1,imess,size(var),rmess,0,tag)
   deallocate(imess,rmess,stat=allocstat)
! if .not.just1, receive the result from processor 1
   if (.not.present(just1)) then
      call phaml_recv(procs,p,imess,ni,rmess,nr,tag)
      phaml_global_sum_intv = imess(1:ni)
      deallocate(imess,stat=allocstat)
   elseif (.not.just1) then
      call phaml_recv(procs,p,imess,ni,rmess,nr,tag)
      phaml_global_sum_intv = imess(1:ni)
      deallocate(imess,stat=allocstat)
   endif
endif

return
end function phaml_global_sum_intv

!          ---------------
subroutine pack_procs_size(procs,ni,nr)
!          ---------------

!----------------------------------------------------
! This subroutine returns the number of integers and reals that must
! be allocated for a message containing proc_info
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(proc_info), intent(in) :: procs
integer, intent(out) :: ni,nr

!----------------------------------------------------
! Begin executable code

ni = 5 + HOSTLEN + size(procs%tids)
nr = 0

end subroutine pack_procs_size

!          ----------
subroutine pack_procs(procs,send_int,first_int,send_real,first_real)
!          ----------

!----------------------------------------------------
! This subroutine packs procs into messages send_int and send_real
! On input, first_* is the first location of send_* to be used for
! packing.  On output, they give the next available location.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(proc_info), intent(in) :: procs
integer, intent(inout) :: send_int(:)
real(my_real), intent(inout) :: send_real(:)
integer, intent(inout) :: first_int, first_real

!----------------------------------------------------
! Local variables:

integer :: i
real(my_real) :: r

!----------------------------------------------------
! Begin executable code

send_int(first_int) = procs%my_proc
send_int(first_int+1) = procs%nproc
send_int(first_int+2) = procs%log_nproc
send_int(first_int+3) = procs%graphics_proc
do i=1,HOSTLEN
   send_int(first_int+3+i) = ichar(procs%hostname(i:i))
end do
send_int(first_int+4+HOSTLEN) = size(procs%tids)
do i=0,size(procs%tids)-1
   send_int(first_int+5+HOSTLEN+i) = procs%tids(i)
end do
first_int = first_int + 5 + HOSTLEN + size(procs%tids)

! There are no real data
! Just to shut up verbose compiler warnings
r = first_real
i = size(send_real)

end subroutine pack_procs

!          ------------
subroutine unpack_procs(procs,send_int,first_int,send_real,first_real)
!          ------------

!----------------------------------------------------
! This subroutine unpacks procs from messages send_int and send_real
! On input, first_* is the first location of send_* to be used for
! unpacking.  On output, they give the next available location.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(proc_info), intent(out) :: procs
integer, intent(in) :: send_int(:)
real(my_real), intent(in) :: send_real(:)
integer, intent(inout) :: first_int, first_real

!----------------------------------------------------
! Local variables:

integer :: i,ntid,allocstat
real :: r

!----------------------------------------------------
! Begin executable code

procs%my_proc = send_int(first_int)
procs%nproc = send_int(first_int+1)
procs%log_nproc = send_int(first_int+2)
procs%graphics_proc = send_int(first_int+3)
do i=1,HOSTLEN
   procs%hostname(i:i) = char(send_int(first_int+3+i))
end do
ntid = send_int(first_int+4+HOSTLEN)
allocate(procs%tids(0:ntid-1),stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in unpack_procs")
   return
endif
do i=0,ntid-1
   procs%tids(i) = send_int(first_int+5+HOSTLEN+i)
end do
first_int = first_int + 5 + HOSTLEN + size(procs%tids)

! There are no real data
! Just to shut up verbose compiler warnings
r = first_real
i = size(send_real)

end subroutine unpack_procs

!          --------------------
subroutine cleanup_unpack_procs(procs)
!          --------------------

!----------------------------------------------------
! This routine deallocates the memory allocated by unpack_procs.  It
! should be called after you are done using the results of unpack_procs
! to avoid a memory leak.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(proc_info), intent(inout) :: procs
!----------------------------------------------------
! Local variables:

integer :: allocstat
!----------------------------------------------------
! Begin executable code

deallocate(procs%tids,stat=allocstat)

end subroutine cleanup_unpack_procs

!        -------
function my_proc(procs)
!        -------

!----------------------------------------------------
! This function returns my processor number
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

integer my_proc
type (proc_info), intent(in) :: procs
!----------------------------------------------------
! Begin executable code

my_proc = procs%my_proc
return
end function my_proc

!        --------
function num_proc(procs)
!        --------

!----------------------------------------------------
! This function returns the number of processors
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

integer num_proc
type (proc_info), intent(in) :: procs
!----------------------------------------------------
! Begin executable code

num_proc = procs%nproc
return
end function num_proc

!        ---------
function log_nproc(procs)
!        ---------

!----------------------------------------------------
! This function returns the base 2 logarithm of the number of processors
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

integer log_nproc
type (proc_info), intent(in) :: procs
!----------------------------------------------------
! Begin executable code

log_nproc = procs%log_nproc
return
end function log_nproc

!        --------
function hostname(procs)
!        --------

!----------------------------------------------------
! This function returns my hostname
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

character(len=HOSTLEN) hostname
type (proc_info), intent(in) :: procs
!----------------------------------------------------
! Begin executable code

hostname = trim(procs%hostname)
return
end function hostname

!        -------------
function graphics_proc(procs)
!        -------------

!----------------------------------------------------
! This function returns the processor number of my graphics server
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

integer :: graphics_proc
type (proc_info), intent(in) :: procs
!----------------------------------------------------
! Begin executable code

graphics_proc = procs%graphics_proc
return
end function graphics_proc

!        -------------------
function slaves_communicator(procs)
!        -------------------

!----------------------------------------------------
! This function returns the MPI communicator for the group of slaves
! In the PVM version, it returns 0 and is useless.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer :: slaves_communicator ! note default integer
type (proc_info), intent(in) :: procs

!----------------------------------------------------
! Local variables:

integer :: i

!----------------------------------------------------
! Begin executable code

slaves_communicator = 0
! just to shut up verbose compiler warnings
i = procs%my_proc

end function slaves_communicator

!          -------------
subroutine phaml_barrier(procs)
!          -------------

!----------------------------------------------------
! This routine waits until all processors (excluding master) get to it
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments
type(proc_info), intent(in) :: procs

!----------------------------------------------------
! Local variables:

integer :: i, p, no_int(1), ni, nr, allocstat
real(my_real) :: no_reals(1)
integer, pointer :: recv_int(:)
real(my_real), pointer :: recv_real(:)
!----------------------------------------------------
! Begin executable code

if (procs%my_proc == MASTER) return

! processor 1 waits for a message from all others and then sends them one

if (procs%my_proc == 1) then
   do i=2,procs%nproc
      call phaml_recv(procs,p,recv_int,ni,recv_real,nr,41)
      if (associated(recv_int)) deallocate(recv_int,stat=allocstat)
      if (associated(recv_real)) deallocate(recv_real,stat=allocstat)
   end do
   do i=2,procs%nproc
      call phaml_send(procs,i,no_int,0,no_reals,0,42)
   end do

! other processors send a message to processor 1 and wait for a reply

else
   call phaml_send(procs,1,no_int,0,no_reals,0,41)
   call phaml_recv(procs,p,recv_int,ni,recv_real,nr,42)
   if (associated(recv_int)) deallocate(recv_int,stat=allocstat)
   if (associated(recv_real)) deallocate(recv_real,stat=allocstat)

endif
return
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

integer :: proc,ni,nr,no_int(1),nsend
real(my_real) :: no_reals(1)
integer, pointer :: recv_int(:)
real(my_real), pointer :: recv_real(:)

!----------------------------------------------------
! Begin executable code

if (really_pause) then

   if (present(dont_pausewatch)) then
      if (.not. dont_pausewatch) call pause_watch(all_watches)
   else
      call pause_watch(all_watches)
   endif

! master waits for the return key and sends a message to the slaves

   if (procs%my_proc == MASTER) then
      write(outunit,"(A)") "press return to continue"
      read *
      if (present(just_one)) then
         if (just_one) then
            nsend = 1
         else
            nsend = procs%nproc
         endif
      else
         nsend = procs%nproc
      endif
      do proc=1,nsend
         call phaml_send(procs,proc,no_int,0,no_reals,0,50)
      end do

! slaves wait for message from master

   else
      call phaml_recv(procs,proc,recv_int,ni,recv_real,nr,50)
   endif

   if (present(dont_pausewatch)) then
      if (.not. dont_pausewatch) call end_pause_watch(all_watches)
   else
      call end_pause_watch(all_watches)
   endif

endif
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
   if (present(reallist)) write(errunit,"(4E18.10E2)") reallist
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

integer :: p
!----------------------------------------------------
! Begin executable code

write(errunit,"(A)")
write(errunit,"(A)") "------------------------------------------------------"
write(errunit,"(3A)") "          PHAML Version ",version_number," ERROR"
write(errunit,"(A)") msg
if (present(msg2)) write(errunit,"(A)") msg2
if (present(intlist)) write(errunit,"(7I11)") intlist
if (present(reallist)) write(errunit,"(4E18.10E2)") reallist
write(errunit,"(A)") "------------------------------------------------------"
write(errunit,"(A)")

if (present(procs)) then
   do p=0,num_proc(procs)
      call phaml_send(procs,p,(/ierr/),1,(/0.0_my_real/),0,ERR_TAG)
   end do
endif

stop
end subroutine fatal

!          ---------------
subroutine sequential_send(imess,ni,rmess,nr)
!          ---------------

!----------------------------------------------------
! dummy version of routine used for sequential programs
!----------------------------------------------------

integer, intent(in) :: imess(:),ni,nr
real(my_real), intent(in) :: rmess(:)

end subroutine sequential_send

!          ---------------
subroutine sequential_recv(imess,ni,rmess,nr)
!          ---------------

!----------------------------------------------------
! dummy version of routine used for sequential programs
!----------------------------------------------------

integer, pointer :: imess(:)
real(my_real), pointer :: rmess(:)
integer, intent(in) :: ni, nr

end subroutine sequential_recv

end module message_passing
