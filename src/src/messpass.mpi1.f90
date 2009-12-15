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
! This version is for spawnless MPI 1
!
! communication tags in this module are of the form 0xx
!----------------------------------------------------

!----------------------------------------------------
! Other modules used:

use global
use mpif_mod
use hash_mod
use hash_eq_mod
use stopwatch
use petsc_init_mod

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
   integer :: all_procs, slaves
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

integer, parameter :: PARALLEL           = MPI1
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

integer itype, rtype  ! the mpi type identifiers for integer and my_real
type(proc_info), pointer, save :: this_processors_procs
integer, save :: max_real_buffer ! largest number of reals safely sent
! size of messages in alltoall routines, in bytes.  Set to 0 to not limit size.
integer, parameter :: small_message = 16000

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
!NAS$ ALIEN

   subroutine mpi_comm_rank(communicator,myrank,ierror)
   integer communicator,myrank,ierror
   end subroutine mpi_comm_rank

   subroutine mpi_finalize(ierror)
   integer ierror
   end subroutine mpi_finalize

   subroutine mpi_get_processor_name(name,resultlen,ierror)
   character(len=*) name
   integer resultlen,ierror
   end subroutine mpi_get_processor_name

   subroutine mpi_get_count(status,datatype,count,ierror)
   integer status(*),datatype,count,ierror
   end subroutine mpi_get_count

   subroutine mpi_iprobe(source,tag,comm,flag,status,ierror)
   logical flag
   integer source,tag,comm,status(*),ierror
   end subroutine mpi_iprobe

   subroutine mpi_pack_size(incount,datatype,comm,size,ierror)
   integer incount,datatype,comm,size,ierror
   end subroutine mpi_pack_size

   subroutine mpi_probe(source,tag,comm,status,ierror)
   integer source,tag,comm,status(*),ierror
   end subroutine mpi_probe

   subroutine mpi_send(buf,count,datatype,dest,tag,comm,ierror)
   integer buf(*)
   integer count,datatype,dest,tag,comm,ierror
   end subroutine mpi_send

   subroutine mpi_recv(buf,count,datatype,source,tag,comm,status,ierror)
   integer buf(*)
   integer count,datatype,source,tag,comm,status(*),ierror
   end subroutine mpi_recv

   subroutine mpi_alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, &
                           recvtype, comm, ierror)
   integer :: sendbuf(*), sendcount, recvbuf(*), recvcount
   integer :: sendtype, recvtype, comm, ierror
   end subroutine mpi_alltoall

end interface

interface mpi_pack
   subroutine mpi_pack_int(inbuf,incount,datatype,outbuf,outsize,position, &
                           comm,ierror)
   use global
   integer :: inbuf(*)
   integer :: outbuf(*)
   integer :: incount,datatype,outsize,position,comm,ierror
   end subroutine mpi_pack_int

   subroutine mpi_pack_my_real(inbuf,incount,datatype,outbuf,outsize,position, &
                              comm,ierror)
   use global
   real(my_real) :: inbuf(*)
   integer :: outbuf(*)
   integer :: incount,datatype,outsize,position,comm,ierror
   end subroutine mpi_pack_my_real
end interface

interface mpi_unpack
   subroutine mpi_unpack_int(inbuf,insize,position,outbuf,outcount, &
                             datatype,comm,ierror)
   use global
   integer :: inbuf(*)
   integer :: outbuf(*)
   integer :: insize,position,outcount,datatype,comm,ierror
   end subroutine mpi_unpack_int

   subroutine mpi_unpack_my_real(inbuf,insize,position,outbuf,outcount, &
                                datatype,comm,ierror)
   use global
   integer :: inbuf(*)
   real(my_real) :: outbuf(*)
   integer :: insize,position,outcount,datatype,comm,ierror
   end subroutine mpi_unpack_my_real
end interface

interface mpi_reduce
   subroutine mpi_reduce_int(sendbuf,recvbuf,count,datatype,op,root, &
                             comm,ierror)
   use global
   integer :: sendbuf(*)
   integer :: recvbuf(*)
   integer :: count,datatype,op,root,comm,ierror
   end subroutine mpi_reduce_int

   subroutine mpi_reduce_my_real(sendbuf,recvbuf,count,datatype,op,root, &
                                 comm,ierror)
   use global
   real(my_real) :: sendbuf(*)
   real(my_real) :: recvbuf(*)
   integer :: count,datatype,op,root,comm,ierror
   end subroutine mpi_reduce_my_real
end interface

interface mpi_allreduce
   subroutine mpi_allreduce_int(sendbuf,recvbuf,count,datatype,op, &
                                comm,ierror)
   use global
   integer :: sendbuf(*)
   integer :: recvbuf(*)
   integer :: count,datatype,op,comm,ierror
   end subroutine mpi_allreduce_int

   subroutine mpi_allreduce_my_real(sendbuf,recvbuf,count,datatype,op, &
                                    comm,ierror)
   use global
   real(my_real) :: sendbuf(*)
   real(my_real) :: recvbuf(*)
   integer :: count,datatype,op,comm,ierror
   end subroutine mpi_allreduce_my_real
end interface

interface mpi_alltoallv
   subroutine mpi_alltoallv_int(sendbuf, sendcounts, sdispls, sendtype, &
                                recvbuf, recvcounts, rdispls, recvtype, &
                                   comm, ierror)
   use global
   integer :: sendbuf(*), sendcounts(*), sdispls(*)
   integer :: recvbuf(*), recvcounts(*), rdispls(*)
   integer :: sendtype, recvtype, comm, ierror
   end subroutine mpi_alltoallv_int

   subroutine mpi_alltoallv_my_real(sendbuf, sendcounts, sdispls, sendtype, &
                                   recvbuf, recvcounts, rdispls, recvtype, &
                                   comm, ierror)
   use global
   real(my_real) :: sendbuf(*), recvbuf(*)
   integer :: sendcounts(*), sdispls(*)
   integer :: recvcounts(*), rdispls(*)
   integer :: sendtype, recvtype, comm, ierror
   end subroutine mpi_alltoallv_my_real
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
integer :: isend(14+HOSTLEN+len(triangle_files))
character(len=MPI_MAX_PROCESSOR_NAME) :: recv_name
character(len=HOSTLEN), allocatable :: from_hostname(:)
character(len=HOSTLEN) :: loc_graphics_host

integer name_length
integer :: slave_group,everything_group

!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

! MPI must be initialized in the user code

! set perminant pointer to procs

this_processors_procs => procs

! identify the integer and real types
! Also set the maximum number of reals that can be safely sent without
! worrying about deadlock.  Assume the buffer size is 64KB, but if you know
! it is something else, change it here.  Also assume real is 4 bytes and
! double precision is 8 bytes.

select case(bit_size(1)/8)
   case(4)
      itype = MPI_INTEGER
   case default
      call fatal("default integer kind does not match an mpi integer type")
      stop
end select

select case(kind(1.0_my_real))
   case(kind(1.0))
      rtype = MPI_REAL
      max_real_buffer = 65536/4
   case(kind(1.0d0))
      rtype = MPI_DOUBLE_PRECISION
      max_real_buffer = 65536/8
   case default
      call fatal("kind my_real does not match an mpi real type")
      stop
end select

! different cases for master, slaves and graphics, determined by my rank

call mpi_comm_rank(MPI_COMM_WORLD,procs%my_proc,info)

if (procs%my_proc == MASTER) then ! I am the master process

   my_type = MASTER

! set the number of processors (slaves)

   procs%nproc = nproc
   procs%log_nproc = 0
   temp = procs%nproc
   do while (temp > 1)
      procs%log_nproc = procs%log_nproc + 1
      temp = temp/2
   end do
   loc_graphics_host = trim(graphics_host)

   procs%all_procs = MPI_COMM_WORLD

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
   rsend(1) = max_blen
   do proc=1,procs%nproc
      call phaml_send(procs,proc,isend,14+HOSTLEN+len(triangle_files),rsend, &
                      1,10)
   end do

else ! I am either a slave or a graphics server

   procs%all_procs = MPI_COMM_WORLD

! receive the code that tells if I am slave or graphics, nproc and graphics

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
   max_blen = rrecv(1)
   deallocate(irecv,rrecv,stat=allocstat)
   procs%log_nproc = 0
   temp = procs%nproc
   do while (temp > 1)
      procs%log_nproc = procs%log_nproc + 1
      temp = temp/2
   end do

endif

! if I'm the master or a slave, identify graphics server

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

   if (ngraph > 0) then
      if (i_have_graphics) then
         if (my_type == MASTER) then
            procs%graphics_proc = procs%nproc + 1
         elseif (my_type == SLAVES) then
            procs%graphics_proc = procs%my_proc + ngraph
         endif
      endif
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
      rsend(1) = max_blen
      call phaml_send(procs,procs%graphics_proc,isend, &
                      14+HOSTLEN+len(triangle_files),rsend,1,10)
   endif

endif

! create a communicator with just the slaves
! non-slaves will get MPI_COMM_NULL in procs%slaves

call mpi_comm_group(MPI_COMM_WORLD,everything_group,info)
call mpi_group_incl(everything_group,procs%nproc,(/(i,i=1,procs%nproc)/), &
                    slave_group,info)
call mpi_comm_create(MPI_COMM_WORLD,slave_group,procs%slaves,info)

! get the host name

call mpi_get_processor_name(recv_name,name_length,info)
if (name_length < MPI_MAX_PROCESSOR_NAME) &
   recv_name(name_length+1:MPI_MAX_PROCESSOR_NAME) = ' '
i = index(recv_name,".")
if (i /= 0) recv_name(i:) = " "
procs%hostname = recv_name

!write(outunit,"(A,I11,2A)") "I am processor number ",procs%my_proc," on host ",procs%hostname

! print host names from master

if (my_type == MASTER) then
   allocate(from_hostname(procs%nproc+ngraph))
   do i=1,procs%nproc+ngraph
      call phaml_recv(procs,from_proc,int_hostname,ni,rrecv,nr,20)
      from_hostname(from_proc) = " "
      do j=1,ni
         from_hostname(from_proc)(j:j) = char(int_hostname(j))
      end do
      if (ni > 0) deallocate(int_hostname,stat=allocstat)
   end do
   do i=1,procs%nproc+ngraph
      write(outunit,"(A,I11,2A)") "Processor number ",i," is on host ",trim(from_hostname(i))
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

! Initialize PETSc

if (my_type == SLAVES) then
   call petsc_init(procs%slaves)
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

integer info, i
real (my_real) rmess(1)
integer imess(1)
integer :: ni, nr, from_proc
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

   call mpi_finalize(info) ! get out of MPI

! graphics engine; just quit

elseif (procs%my_proc > procs%nproc) then

   call mpi_finalize(info) ! get out of MPI

! slave

else
   rmess(1)=float(procs%my_proc)
   call phaml_send(procs,MASTER,imess,0,rmess,1,30)
   if (procs%graphics_proc /= -1) then
      imess(1)=GRAPHICS_TERMINATE
      call phaml_send(procs,procs%graphics_proc,imess,1,rmess,0,101)
   endif
   call petsc_finalize ! terminate PETSc; only slaves call petsc_init
   call mpi_finalize(info) ! get out of MPI
endif

nullify(this_processors_procs)

return
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
integer :: part1, part2, part3, nbytes, pos, comm, dest
integer :: sizes(2)
integer, allocatable :: buffer(:)
integer :: allocstat
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

if (print_comm) then
   write(outunit,"(3(A,I11))") "proc ",procs%my_proc," send to proc ",proc," tag ",tag
!   write(outunit,"(I11,A)")  ni," integers "
!   write(outunit,"(7I11)")  int_message(1:ni)
!   write(outunit,"(I11,A)")  nr," reals "
!   write(outunit,"(4E18.10E2)")  real_message(1:nr)
endif

ni_default_kind = ni; nr_default_kind = nr

! set the communicator and destination; special case for graphics

comm = procs%all_procs
dest = proc

! initialize the buffer for sending

call mpi_pack_size( 2,itype,comm,part1,info)
call mpi_pack_size(ni,itype,comm,part2,info)
call mpi_pack_size(nr,rtype,comm,part3,info)
nbytes = part1+part2+part3

if (itype == MPI_INTEGER) then
   allocate(buffer((part1+part2+part3)/4),stat=allocstat)
   if (allocstat/=0) then
      ierr = ALLOC_FAILED
      call fatal("allocation of buffer failed in phaml_send", &
                  intlist=(/allocstat/))
      return
   endif
else
   call fatal("Didn't match itype in phaml_send")
   stop
endif

! pack the data

sizes = (/ ni, nr /)

pos = 0
call mpi_pack(sizes,2,itype,buffer,nbytes,pos,comm,info)
if (ni>0) call mpi_pack(int_message,ni_default_kind,itype,buffer,nbytes, &
                        pos,comm,info)
if (nr>0) call mpi_pack(real_message,nr_default_kind,rtype,buffer,nbytes, &
                        pos,comm,info)

! send it

tag_default_kind = tag
call mpi_send(buffer,pos,MPI_PACKED,dest,tag_default_kind,comm,info)

! free the buffer memory

deallocate(buffer,stat=allocstat)

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

integer info,tag_default_kind,ni_default_kind,nr_default_kind,allocstat, &
        status(MPI_STATUS_SIZE),nbyte,status2(MPI_STATUS_SIZE),pos
integer :: sizes(2)
integer, allocatable :: buffer(:)
logical :: message_waiting
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

if (print_comm) then
   if (.not.present(noblock)) then
      write(outunit,"(2(A,I11))") "proc ",procs%my_proc," waiting to receive tag ",tag
   endif
endif

! check for a message that some processor had a fatal error

call mpi_iprobe(MPI_ANY_SOURCE,ERR_TAG,procs%all_procs,message_waiting, &
                status,info)
if (message_waiting) then
   write(errunit,"(A)") "Some processor encountered a fatal error.  Halting."
endif

! probe the message to find the size for the buffer

tag_default_kind = tag

if (present(noblock)) then
   if (noblock) then
      call mpi_iprobe(MPI_ANY_SOURCE,tag_default_kind,procs%all_procs, &
                      message_waiting,status,info)
   else
      call mpi_probe(MPI_ANY_SOURCE,tag_default_kind,procs%all_procs, &
                     status,info)
      message_waiting = .true.
   endif
else
   call mpi_probe(MPI_ANY_SOURCE,tag_default_kind,procs%all_procs, &
                  status,info)
   message_waiting = .true.
endif

! if message_waiting is false, it was an empty nonblocking receive

if (.not. message_waiting) then
   proc=0; ni=0; nr=0
   return
endif

! determine size of buffer and allocate

call mpi_get_count(status,MPI_PACKED,nbyte,info)
if (itype == MPI_INTEGER) then
   allocate(buffer(nbyte/4),stat=allocstat)
   if (allocstat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in phaml_recv")
      return
   endif
else
   call fatal("Didn't identify itype in phaml_recv")
endif

! receive the message

call mpi_recv(buffer,nbyte,MPI_PACKED,status(MPI_SOURCE),tag_default_kind, &
              procs%all_procs,status2,info)

! identify the sender

proc = status(MPI_SOURCE)

! unpack the data, creating space along the way

pos = 0
call mpi_unpack(buffer,nbyte,pos,sizes,2,itype,procs%all_procs,info)
ni = sizes(1); nr=sizes(2)
ni_default_kind = ni; nr_default_kind = nr
if (ni>0) then
   allocate(int_message(ni),stat=allocstat)
   if (allocstat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in phaml_recv")
      return
   endif
   call mpi_unpack(buffer,nbyte,pos,int_message,ni_default_kind,itype, &
                   procs%all_procs,info)
else
   nullify(int_message)
endif
if (nr>0) then
   allocate(real_message(nr),stat=allocstat)
   if (allocstat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in phaml_recv")
      return
   endif
   call mpi_unpack(buffer,nbyte,pos,real_message,nr_default_kind,rtype, &
                   procs%all_procs,info)
else
   nullify(real_message)
endif

! free the buffer

deallocate(buffer,stat=allocstat)

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

integer :: i, jerr, astat, i1, i2, loop, isend, sm
integer :: displs(procs%nproc), nsenda(1), irecv(procs%nproc)
!----------------------------------------------------
! Begin executable code

if (print_comm) then
   write(outunit,"(2(A,I11))") "proc ",procs%my_proc," in phaml_alltoall_int, tag ",tag
endif

! use an all-gather to find out how much data each processor is sending

nsenda(1) = nsend
call mpi_allgather(nsenda,1,itype,nrecv,1,itype,procs%slaves,jerr)

! determine the starting point of each received message

displs(1) = 0
do i=2,procs%nproc
   displs(i) = displs(i-1) + nrecv(i-1)
end do

! allocate the receive buffer

allocate(recv(displs(procs%nproc)+nrecv(procs%nproc)),stat=astat)
if (astat/=0) then
   ierr = ALLOC_FAILED
   call fatal("allocation of buffer failed in phaml_alltoall",intlist=(/astat/))
   return
endif

! send/recv the data

sm = small_message/4 ! 4 bytes per integer
sm = sm/procs%nproc
if (small_message == 0 .or. maxval(nrecv) <= sm .or. &
    1+(maxval(nrecv)-1)/sm > 50) then

! all in one message if the messages are small enough or not limiting size
! or too many messages

   call mpi_allgatherv(send,nsend,itype,recv,nrecv,displs,itype, &
                       procs%slaves,jerr)

else ! send as several small messages

   do loop = 1,1+(maxval(nrecv)-1)/sm

! determine number sent and received from each processor this time

      isend = max(min(nsend-(loop-1)*sm,sm),0)
      do i=1,procs%nproc
         irecv(i) = max(min(nrecv(i)-(loop-1)*sm,sm),0)
      end do

! determine starting and ending points in this sent message

      i1 = min(nsend,1+(loop-1)*sm)
      i2 = i1+isend-1

! exchange messages

      call mpi_allgatherv(send(i1:i2),isend,itype,recv,irecv,displs,itype, &
                          procs%slaves,jerr)

! update the displacements for the next message

      displs = displs + irecv

   end do
endif

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

integer :: i, jerr, astat, i1, i2, loop, isend, sm
integer :: displs(procs%nproc), nsenda(1), irecv(procs%nproc)
!----------------------------------------------------
! Begin executable code

if (print_comm) then
   write(outunit,"(2(A,I11))") "proc ",procs%my_proc," in phaml_alltoall_real, tag ",tag
endif

! use an all-gather to find out how much data each processor is sending

nsenda(1) = nsend
call mpi_allgather(nsenda,1,itype,nrecv,1,itype,procs%slaves,jerr)

! determine the starting point of each received message

displs(1) = 0
do i=2,procs%nproc
   displs(i) = displs(i-1) + nrecv(i-1)
end do

! allocate the receive buffer

allocate(recv(displs(procs%nproc)+nrecv(procs%nproc)),stat=astat)
if (astat/=0) then
   ierr = ALLOC_FAILED
   call fatal("allocation of buffer failed in phaml_alltoall",intlist=(/astat/))
   return
endif

! send/recv the data

if (my_real == kind(0.0)) then
   sm = small_message/4
else
   sm = small_message/8
endif
sm = sm/procs%nproc
if (small_message == 0 .or. maxval(nrecv) <= sm .or. &
    1+(maxval(nrecv)-1)/sm > 50) then

! all in one message if the messages are small enough or not limiting size
! or too many messages

   call mpi_allgatherv(send,nsend,rtype,recv,nrecv,displs,rtype, &
                       procs%slaves,jerr)

else ! send as several small messages

   do loop = 1,1+(maxval(nrecv)-1)/sm

! determine number sent and received from each processor this time

      isend = max(min(nsend-(loop-1)*sm,sm),0)
      do i=1,procs%nproc
         irecv(i) = max(min(nrecv(i)-(loop-1)*sm,sm),0)
      end do

! determine starting and ending points in this sent message

      i1 = min(nsend,1+(loop-1)*sm)
      i2 = i1+isend-1

! exchange messages

      call mpi_allgatherv(send(i1:i2),isend,rtype,recv,irecv,displs,rtype, &
                          procs%slaves,jerr)


! update the displacements for the next message

      displs = displs + irecv

   end do
endif

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

integer :: i, jerr, astat, loop, sm, messsize(1)
integer :: sdispls(procs%nproc), rdispls(procs%nproc), isend(procs%nproc), &
           irecv(procs%nproc)
!----------------------------------------------------
! Begin executable code

if (print_comm) then
   write(outunit,"(2(A,I11))") "proc ",procs%my_proc," in phaml_alltoall_intv, tag ",tag
endif

! determine the starting point of each sent message

sdispls(1) = 0
do i=2,procs%nproc
   sdispls(i) = sdispls(i-1) + nsend(i-1)
end do

! use an all-to-all to find out how much data each processor is sending

call mpi_alltoall(nsend,1,itype,nrecv,1,itype,procs%slaves,jerr)

! determine the starting point of each received message

rdispls(1) = 0
do i=2,procs%nproc
   rdispls(i) = rdispls(i-1) + nrecv(i-1)
end do

! allocate the receive buffer

allocate(recv(rdispls(procs%nproc)+nrecv(procs%nproc)),stat=astat)
if (astat/=0) then
   ierr = ALLOC_FAILED
   call fatal("allocation of buffer failed in phaml_alltoall",intlist=(/astat/))
   return
endif

! determine the small message size and biggest actual message size

sm = small_message/4
call mpi_allreduce((/sum(nsend)/),messsize,1,itype,MPI_MAX,procs%slaves,jerr)

! send/recv the data

if (small_message == 0 .or. messsize(1) <= sm .or. &
    1+(messsize(1)-1)/sm > 50) then

! all in one message if the messages are small enough or not limiting size
! or too many messages

   call mpi_alltoallv(send,nsend,sdispls,itype,recv,nrecv,rdispls, &
                      itype,procs%slaves,jerr)

else ! send as several small messages

   do loop = 1,1+(messsize(1)-1)/sm

! determine number sent to and received from each processor this time

      do i=1,procs%nproc
         isend(i) = max(min(nsend(i)-(loop-1)*sm,sm),0)
         irecv(i) = max(min(nrecv(i)-(loop-1)*sm,sm),0)
      end do

! exchange messages

      call mpi_alltoallv(send,isend,sdispls,itype,recv,irecv,rdispls, &
                         itype,procs%slaves,jerr)

! update the displacements for the next message

      sdispls = sdispls + isend
      rdispls = rdispls + irecv

   end do
endif

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

integer :: i, jerr, astat, loop, sm, messsize(1)
integer :: sdispls(procs%nproc), rdispls(procs%nproc), isend(procs%nproc), &
           irecv(procs%nproc)
!----------------------------------------------------
! Begin executable code

if (print_comm) then
   write(outunit,"(2(A,I11))") "proc ",procs%my_proc," in phaml_alltoall_realv, tag ",tag
endif

! determine the starting point of each sent message

sdispls(1) = 0
do i=2,procs%nproc
   sdispls(i) = sdispls(i-1) + nsend(i-1)
end do

! use an all-to-all to find out how much data each processor is sending

call mpi_alltoall(nsend,1,itype,nrecv,1,itype,procs%slaves,jerr)

! determine the starting point of each received message

rdispls(1) = 0
do i=2,procs%nproc
   rdispls(i) = rdispls(i-1) + nrecv(i-1)
end do

! allocate the receive buffer

allocate(recv(rdispls(procs%nproc)+nrecv(procs%nproc)),stat=astat)
if (astat/=0) then
   ierr = ALLOC_FAILED
   call fatal("allocation of buffer failed in phaml_alltoall",intlist=(/astat/))
   return
endif

! determine the small message size and biggest actual message size

if (my_real == kind(0.0)) then
   sm = small_message/4
else
   sm = small_message/8
endif
call mpi_allreduce((/sum(nsend)/),messsize,1,itype,MPI_MAX,procs%slaves,jerr)

! send/recv the data

if (small_message == 0 .or. messsize(1) <= sm .or. &
    1+(messsize(1)-1)/sm > 50) then

! all in one message if the messages are small enough or not limiting size
! or too many messages

   call mpi_alltoallv(send,nsend,sdispls,rtype,recv,nrecv,rdispls, &
                      rtype,procs%slaves,jerr)

else ! send as several small messages

   do loop = 1,1+(messsize(1)-1)/sm

! determine number sent to and received from each processor this time

      do i=1,procs%nproc
         isend(i) = max(min(nsend(i)-(loop-1)*sm,sm),0)
         irecv(i) = max(min(nrecv(i)-(loop-1)*sm,sm),0)
      end do

! exchange messages

      call mpi_alltoallv(send,isend,sdispls,rtype,recv,irecv,rdispls, &
                         rtype,procs%slaves,jerr)

! update the displacements for the next message

      sdispls = sdispls + isend
      rdispls = rdispls + irecv

   end do
endif

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

real (my_real) :: result(1)
integer info
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

if (print_comm) then
   write(outunit,"(2(A,I11))") "proc ",procs%my_proc," in reduction operation, tag ",tag
endif

! for just1, use mpi_reduce and return

if (present(just1)) then; if (just1) then
   call mpi_reduce((/var/),result,1,rtype,MPI_MAX,0,procs%slaves,info)
   if (procs%my_proc == 1) then
      phaml_global_max_real = result(1)
   else
      phaml_global_max_real = var
   endif
   return
endif; endif

! if not just1, use mpi_allreduce

call mpi_allreduce((/var/),result,1,rtype,MPI_MAX,procs%slaves,info)
phaml_global_max_real = result(1)

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

integer :: result(1)
integer info
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

if (print_comm) then
   write(outunit,"(2(A,I11))") "proc ",procs%my_proc," in reduction operation, tag ",tag
endif

! for just1, use mpi_reduce and return

if (present(just1)) then; if (just1) then
   call mpi_reduce((/var/),result,1,itype,MPI_MAX,0,procs%slaves,info)
   if (procs%my_proc == 1) then
      phaml_global_max_int = result(1)
   else
      phaml_global_max_int = var
   endif
   return
endif; endif

! if not just1, use mpi_allreduce

call mpi_allreduce((/var/),result,1,itype,MPI_MAX,procs%slaves,info)
phaml_global_max_int = result(1)

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

real (my_real) :: result(1)
integer info
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

if (print_comm) then
   write(outunit,"(2(A,I11))") "proc ",procs%my_proc," in reduction operation, tag ",tag
endif

! for just1, use mpi_reduce and return

if (present(just1)) then; if (just1) then
   call mpi_reduce((/var/),result,1,rtype,MPI_MIN,0,procs%slaves,info)
   if (procs%my_proc == 1) then
      phaml_global_min_real = result(1)
   else
      phaml_global_min_real = var
   endif
   return
endif; endif

! if not just1, use mpi_allreduce

call mpi_allreduce((/var/),result,1,rtype,MPI_MIN,procs%slaves,info)
phaml_global_min_real = result(1)

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

integer :: result(1)
integer info
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

if (print_comm) then
   write(outunit,"(2(A,I11))") "proc ",procs%my_proc," in reduction operation, tag ",tag
endif

! for just1, use mpi_reduce and return

if (present(just1)) then; if (just1) then
   call mpi_reduce((/var/),result,1,itype,MPI_MIN,0,procs%slaves,info)
   if (procs%my_proc == 1) then
      phaml_global_min_int = result(1)
   else
      phaml_global_min_int = var
   endif
   return
endif; endif

! if not just1, use mpi_allreduce

call mpi_allreduce((/var/),result,1,itype,MPI_MIN,procs%slaves,info)
phaml_global_min_int = result(1)

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

real (my_real) :: result(1)
integer info
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

if (print_comm) then
   write(outunit,"(2(A,I11))") "proc ",procs%my_proc," in reduction operation, tag ",tag
endif

! for just1, use mpi_reduce and return

if (present(just1)) then; if (just1) then
   call mpi_reduce((/var/),result,1,rtype,MPI_SUM,0,procs%slaves,info)
   if (procs%my_proc == 1) then
      phaml_global_sum_real = result(1)
   else
      phaml_global_sum_real = var
   endif
   return
endif; endif

! if not just1, use mpi_allreduce

call mpi_allreduce((/var/),result,1,rtype,MPI_SUM,procs%slaves,info)
phaml_global_sum_real = result(1)

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

integer :: result(1)
integer info
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

if (print_comm) then
   write(outunit,"(2(A,I11))") "proc ",procs%my_proc," in reduction operation, tag ",tag
endif

! for just1, use mpi_reduce and return

if (present(just1)) then; if (just1) then
   call mpi_reduce((/var/),result,1,itype,MPI_SUM,0,procs%slaves,info)
   if (procs%my_proc == 1) then
      phaml_global_sum_int = result(1)
   else
      phaml_global_sum_int = var
   endif
   return
endif; endif

! if not just1, use mpi_allreduce

call mpi_allreduce((/var/),result,1,itype,MPI_SUM,procs%slaves,info)
phaml_global_sum_int = result(1)

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
   real, dimension(size(var)) :: phaml_global_sum_realv
!----------------------------------------------------

!----------------------------------------------------
! Local variables:

real(my_real) :: result(size(var))
integer info
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

if (print_comm) then
   write(outunit,"(2(A,I11))") "proc ",procs%my_proc," in reduction operation, tag ",tag
endif

! for just1, use mpi_reduce and return

if (present(just1)) then; if (just1) then
   call mpi_reduce(var,result,size(var),rtype,MPI_SUM,0,procs%slaves,info)
   if (procs%my_proc == 1) then
      phaml_global_sum_realv = result
   else
      phaml_global_sum_realv = var
   endif
   return
endif; endif

! if not just1, use mpi_allreduce

call mpi_allreduce(var,result,size(var),rtype,MPI_SUM,procs%slaves,info)
phaml_global_sum_realv = result

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

integer :: result(size(var))
integer info
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

if (print_comm) then
   write(outunit,"(2(A,I11))") "proc ",procs%my_proc," in reduction operation, tag ",tag
endif

! for just1, use mpi_reduce and return

if (present(just1)) then; if (just1) then
   call mpi_reduce(var,result,size(var),itype,MPI_SUM,0,procs%slaves,info)
   if (procs%my_proc == 1) then
      phaml_global_sum_intv = result
   else
      phaml_global_sum_intv = var
   endif
   return
endif; endif

! if not just1, use mpi_allreduce

call mpi_allreduce(var,result,size(var),itype,MPI_SUM,procs%slaves,info)
phaml_global_sum_intv = result

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

ni = procs%my_proc ! just to shut up verbose compilers
ni = 6 + HOSTLEN
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
send_int(first_int+4+HOSTLEN) = procs%all_procs
send_int(first_int+5+HOSTLEN) = procs%slaves

first_int = first_int + 6 + HOSTLEN

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

integer :: i
real(my_real) :: r

!----------------------------------------------------
! Begin executable code

procs%my_proc = send_int(first_int)
procs%nproc = send_int(first_int+1)
procs%log_nproc = send_int(first_int+2)
procs%graphics_proc = send_int(first_int+3)
do i=1,HOSTLEN
   procs%hostname(i:i) = char(send_int(first_int+3+i))
end do
procs%all_procs = send_int(first_int+4+HOSTLEN)
procs%slaves = send_int(first_int+5+HOSTLEN)

first_int = first_int + 6 + HOSTLEN

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

!----------------------------------------------------
! Begin executable code

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
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer :: slaves_communicator ! note default integer
type (proc_info), intent(in) :: procs
!----------------------------------------------------
! Begin executable code

slaves_communicator = procs%slaves

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
