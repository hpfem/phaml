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
! This module contains the interface to Zoltan, i.e., query functions
! and other useful stuff
!
! communication tags in this module are of the form 6xx
!----------------------------------------------------

!----------------------------------------------------
! Other modules used:

use global
use message_passing
use hash_mod
use gridtype_mod
use zoltan
use zoltan_types

implicit none
private
public zoltan_create_lb, zoltan_destroy_lb, zoltan_init, zoltan_partition, &
       Zoltan_Struct
public zoltan_child_order

! Pointers to the data used in the query functions.  Since all the Zoltan
! calls go through this module, I can just set these pointers before calling
! the load balancing routine instead of passing them through the data argument

type(grid_type), pointer :: global_grid
type(proc_info), pointer :: global_procs

! Element owner, needed by some Zoltan methods

integer, allocatable :: elem_owner(:)
logical, save :: need_owner = .false.

logical, save :: first_time

integer, save :: zoltan_method = -1

! interface blocks for zoltanParams

interface
   subroutine Zf90_zoltanparams_read_file(lb_addr,nbytes,ifile,nchar, &
                                          communicator,ierr)
   use zoltan
   integer(Zoltan_INT), intent(in) :: lb_addr(*), ifile(*)
   integer(Zoltan_INT), intent(in) :: nbytes, nchar, communicator
   integer(Zoltan_INT), intent(out) :: ierr
   end subroutine Zf90_zoltanparams_read_file
   subroutine zoltanparams_using_graph(retval)
   use zoltan
   integer(Zoltan_INT), intent(out) :: retval
   end subroutine zoltanparams_using_graph
end interface

contains

!        ---------
function z_num_obj(dummy, ierr)
!        ---------

!----------------------------------------------------
! This routine returns the number of leaf elements owned by this processor
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer(Zoltan_INT) :: z_num_obj
integer(Zoltan_INT), dimension(*), intent(in) :: dummy
integer(Zoltan_INT), intent(out) :: ierr

!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

ierr = ZOLTAN_OK

! on first call, pretend only processor 1 has the grid

if (first_time .and. my_proc(global_procs) /= 1) then
   z_num_obj = 0
   return
endif

z_num_obj = global_grid%nelem_leaf_own

return
end function z_num_obj

!          ----------
subroutine z_obj_list(dummy, nge, nle, gid, lid, wdim, wgts, ierr)
!          ----------

!----------------------------------------------------
! This routine returns the list of leaf elements owned by this processor.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer(Zoltan_INT), dimension(*), intent(in) :: dummy
integer(Zoltan_INT), intent(in) :: nge, nle
integer(Zoltan_INT), intent(out) :: gid(*)
integer(Zoltan_INT), intent(out) :: lid(*)
integer(Zoltan_INT), intent(in) :: wdim
real(Zoltan_FLOAT), intent(out) :: wgts(*)
integer(Zoltan_INT), intent(out) :: ierr

!----------------------------------------------------
! Local variables:

integer :: nleaf, i, loc_lid
type(hash_key) :: loc_gid

!----------------------------------------------------
! Begin executable code

! on first call, pretend only processor 1 has the grid

if (first_time .and. my_proc(global_procs) /= 1) then
   ierr = ZOLTAN_OK
   return
endif

! create the lists of IDs of owned leaf elements

nleaf = 0
call leaf_list(gid,lid,nleaf)

ierr = ZOLTAN_OK

! set the weights, if requested

if (wdim /= 0) then

! the list of gids just created has all the leaves owned by the processor's
! first partition, followed by all the leaves owned by its second partition,
! etc.  Go through the partitions, moving to the next partition when you
! encounter an unowned element.

   do i=1,nleaf
      loc_gid = hash_unpack_key(gid(1:KEY_SIZE),(i-1)*KEY_SIZE+1)
      loc_lid = hash_decode_key(loc_gid,global_grid%elem_hash)
      if (loc_lid == HASH_NOT_FOUND) then
         ierr = ZOLTAN_FATAL
         call fatal("Failed to find local ID in z_obj_list")
         return
      endif

      wgts((i-1)*wdim+1) = global_grid%element(loc_lid)%weight
      wgts((i-1)*wdim+2:i*wdim) = 0.0_Zoltan_FLOAT

   end do
endif

return
end subroutine z_obj_list

!          ---------
subroutine leaf_list(gid,lid,nleaf)
!          ---------

!----------------------------------------------------
! Calls leaf_list_recur with each child of the root
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(inout) :: gid(*)
integer, intent(inout) :: lid(*),nleaf

!----------------------------------------------------
! Local variables:

integer :: i

!----------------------------------------------------
! Begin executable code

i=global_grid%head_level_elem(1)
do while (i /= END_OF_LIST)
   call leaf_list_recur(i,gid,lid,nleaf)
   i = global_grid%element(i)%next
end do

end subroutine leaf_list

!                    ---------------
recursive subroutine leaf_list_recur(tree_node,gid,lid,nleaf)
!                    ---------------

!----------------------------------------------------
! This routine creates gid and lid lists of owned leaf elements
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: tree_node
integer, intent(inout) :: gid(*)
integer, intent(inout) :: lid(*),nleaf

!----------------------------------------------------
! Local variables:

integer :: child(MAX_CHILD), i, lgid(KEY_SIZE), allc(MAX_CHILD)

!----------------------------------------------------
! Begin executable code

allc = ALL_CHILDREN
child = get_child_lid(global_grid%element(tree_node)%gid,allc, &
                      global_grid%elem_hash)
if (child(1) == NO_CHILD) then
   if (global_grid%element(tree_node)%iown) then
      nleaf = nleaf + 1
      call hash_pack_key(global_grid%element(tree_node)%gid,lgid,1)
      gid((nleaf-1)*KEY_SIZE+1:nleaf*KEY_SIZE) = lgid
      lid(nleaf) = 0
   endif
else
   do i=1,MAX_CHILD
      if (child(i) == NO_CHILD) cycle
      call leaf_list_recur(child(i),gid,lid,nleaf)
   end do
endif

return
end subroutine leaf_list_recur

!        ----------
function z_num_geom(dummy, ierr)
!        ----------

!----------------------------------------------------
! This routine returns 2, the number of geometry coordinates
! RESTRICTION 2D
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer(Zoltan_INT) :: z_num_geom
integer(Zoltan_INT), dimension(*), intent(in) :: dummy
integer(Zoltan_INT), intent(out) :: ierr

!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

z_num_geom = 2
ierr = ZOLTAN_OK

return
end function z_num_geom

!          ------
subroutine z_geom(dummy,nge,nle,gid,lid,geom,ierr)
!          ------

!----------------------------------------------------
! This routine returns the midpoint of element gid
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer(Zoltan_INT), dimension(*), intent(in) :: dummy
integer(Zoltan_INT), intent(in) :: nge, nle
integer(Zoltan_INT), intent(in) :: gid(*)
integer(Zoltan_INT), intent(in) :: lid(*)
real(Zoltan_DOUBLE), intent(out) :: geom(*)
integer(Zoltan_INT), intent(out) :: ierr

!----------------------------------------------------
! Local variables:

integer :: i, loc_lid
type(hash_key) :: loc_gid
real(Zoltan_DOUBLE) :: sumit1, sumit2
!----------------------------------------------------
! Begin executable code

! find the element

loc_gid = hash_unpack_key(gid(1:KEY_SIZE),1)
loc_lid = hash_decode_key(loc_gid,global_grid%elem_hash)
if (loc_lid == HASH_NOT_FOUND) then
   ierr = ZOLTAN_FATAL
   call fatal("failed to find element in z_geom",procs=global_procs)
   return
endif

! RESTRICTION 2D
sumit1 = 0.0_Zoltan_DOUBLE
sumit2 = 0.0_Zoltan_DOUBLE
do i=1,VERTICES_PER_ELEMENT
   sumit1 = sumit1 + global_grid%vertex(global_grid%element(loc_lid)%vertex(i))%coord%x
   sumit2 = sumit2 + global_grid%vertex(global_grid%element(loc_lid)%vertex(i))%coord%y
end do
geom(1) = sumit1/VERTICES_PER_ELEMENT
geom(2) = sumit2/VERTICES_PER_ELEMENT
ierr = ZOLTAN_OK

return
end subroutine z_geom

!        ----------
function z_num_edge(dummy, nge, nle, gid, lid, ierr)
!        ----------

!----------------------------------------------------
! This routine returns the number of sides that are not boundary
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer(Zoltan_INT) :: z_num_edge
integer(Zoltan_INT), intent(in) :: nge, nle
integer(Zoltan_INT), dimension(*), intent(in) :: dummy
integer(Zoltan_INT), intent(in) :: gid(*)
integer(Zoltan_INT), intent(in) :: lid(*)
integer(Zoltan_INT), intent(out) :: ierr

!----------------------------------------------------
! Local variables:

integer :: i, count, parent, loc_lid
type(hash_key) :: loc_gid
!----------------------------------------------------
! Begin executable code

! find the element

loc_gid = hash_unpack_key(gid(1:KEY_SIZE),1)
loc_lid = hash_decode_key(loc_gid,global_grid%elem_hash)
if (loc_lid == HASH_NOT_FOUND) then
   ierr = ZOLTAN_FATAL
   call fatal("failed to find element in z_num_edge",procs=global_procs)
   return
endif

! look at the neighbors and count how many are not BOUNDARY

if (global_grid%element(loc_lid)%level == 1) then
   count = 0
   do i=1,3
      if (global_grid%initial_neighbor(i,loc_lid) /= BOUNDARY) count = count + 1
   end do
else
   count = 1
   if (.not. (global_grid%element(loc_lid)%mate == BOUNDARY)) count = count + 1
   parent=hash_decode_key(global_grid%element(loc_lid)%gid/2,global_grid%elem_hash)
   if (.not. (global_grid%element(parent)%mate == BOUNDARY)) count = count + 1
endif
z_num_edge = count
ierr = ZOLTAN_OK

return
end function z_num_edge

!          -------
subroutine z_edges(dummy, nge, nle, gid, lid, nbor, proc, wdim, wgt, ierr)
!          -------

!----------------------------------------------------
! This routine returns a list of sides that are not boundary
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer(Zoltan_INT), dimension(*), intent(in) :: dummy
integer(Zoltan_INT), intent(in) :: nge, nle
integer(Zoltan_INT), intent(in) :: gid(*)
integer(Zoltan_INT), intent(in) :: lid(*), wdim
integer(Zoltan_INT), intent(out) :: nbor(*), proc(*)
real(Zoltan_FLOAT), intent(out) :: wgt(*)
integer(Zoltan_INT), intent(out) :: ierr

!----------------------------------------------------
! Local variables:

integer :: i, count, neigh, loc_lid
integer :: loc_nbor(EDGES_PER_ELEMENT*KEY_SIZE), nbor_lid(EDGES_PER_ELEMENT)
type(hash_key) :: loc_gid
!----------------------------------------------------
! Begin executable code

if (wdim /= 0) then
   call fatal("edge weights for partition requested.  This requires", &
              "modifying subroutine z_edges in module zoltan_interf",procs=global_procs)
   stop
endif

! find the element

loc_gid = hash_unpack_key(gid(1:KEY_SIZE),1)
loc_lid = hash_decode_key(loc_gid,global_grid%elem_hash)
if (loc_lid == HASH_NOT_FOUND) then
   ierr = ZOLTAN_FATAL
   call fatal("failed to find element in z_edges",procs=global_procs)
   return
endif

! get the neighbors' lid

if (global_grid%element(loc_lid)%level == 1) then
   nbor_lid = global_grid%initial_neighbor(:,loc_lid)
else
   nbor_lid = get_neighbors(loc_lid,global_grid)
endif

! pack the gids

count = 0
do i=1,3
   neigh = nbor_lid(i)
   if (neigh /= BOUNDARY) then
      call hash_pack_key(global_grid%element(neigh)%gid,loc_nbor, &
                         count*KEY_SIZE+1)
      count = count + 1
      if (.not. allocated(elem_owner)) then
         call fatal("z_edges called without setting elem_owner",procs=global_procs)
         stop
      endif
      proc(count) = elem_owner(neigh) - 1 ! Zoltan starts at 0
   endif
end do
nbor(1:count*KEY_SIZE) = loc_nbor(1:count*KEY_SIZE)
ierr = ZOLTAN_OK

return
end subroutine z_edges

!        ------------
function z_num_coarse(dummy, ierr)
!        ------------

!----------------------------------------------------
! This routine returns the number of elements in the initial grid
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer(Zoltan_INT) :: z_num_coarse
integer(Zoltan_INT), dimension(*), intent(in) :: dummy
integer(Zoltan_INT), intent(out) :: ierr

!----------------------------------------------------
! Local variables:

integer :: i, elem
!----------------------------------------------------
! Begin executable code

if (.not. associated(global_grid%initial_neighbor)) then
   ierr = ZOLTAN_FATAL
   call fatal("z_num_coarse called before grid initialized",procs=global_procs)
else

! on first call, pretend only processor 1 has the grid

   if (first_time .and. my_proc(global_procs) /= 1) then
      z_num_coarse = 0
      ierr = ZOLTAN_OK

   else

      elem = global_grid%head_level_elem(1)
      i = 0
      do while (elem /= END_OF_LIST)
         i = i + 1
         elem = global_grid%element(elem)%next
      end do
      z_num_coarse = i
      ierr = ZOLTAN_OK
   endif

endif

return
end function z_num_coarse

!          -----------------
subroutine z_coarse_obj_list(dummy, nge, nle, gid, lid, mine, num_vert, verts, &
                             in_order, in_vert, out_vert, ierr)
!          -----------------

!----------------------------------------------------
! This routine returns the list of elements in the coarse grid.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer(Zoltan_INT), dimension(*), intent(in) :: dummy
integer(Zoltan_INT), intent(in) :: nge, nle
integer(Zoltan_INT), intent(out) :: gid(*)
integer(Zoltan_INT), intent(out) :: lid(*), mine(*), num_vert(*), &
                                verts(*), &
                                in_vert(*), out_vert(*)
integer(Zoltan_INT), intent(out) :: in_order, ierr

!----------------------------------------------------
! Local variables:

integer :: i, j, elem, last_elem, lgid(KEY_SIZE), &
           temp_verts(VERTICES_PER_ELEMENT*KEY_SIZE)

!----------------------------------------------------
! Begin executable code

in_order = 0

! on first call, pretend only processor 1 has the grid

if (first_time .and. my_proc(global_procs) /= 1) then
   ierr = ZOLTAN_OK
   return
endif

if (.not. associated(global_grid%initial_neighbor)) then
   ierr = ZOLTAN_FATAL
   call fatal("z_coarse_obj_list called before grid initialized",procs=global_procs)
else

   i = 0
   elem = global_grid%head_level_elem(1)
   do while (elem /= END_OF_LIST)
      i = i + 1
      call hash_pack_key(global_grid%element(elem)%gid,lgid,1)
      gid((i-1)*KEY_SIZE+1:i*KEY_SIZE) = lgid
      lid(i) = 0
      if (global_grid%element(elem)%iown) then
         mine(i) = 1
      else
         mine(i) = 0
      endif
      num_vert(i) = VERTICES_PER_ELEMENT
      do j=1,VERTICES_PER_ELEMENT
         call hash_pack_key(global_grid%vertex(global_grid%element(elem)%vertex(j))%gid, &
                            temp_verts,(j-1)*KEY_SIZE+1)
      end do
      verts((i-1)*VERTICES_PER_ELEMENT*KEY_SIZE+1:i*VERTICES_PER_ELEMENT*KEY_SIZE) = temp_verts
      call hash_pack_key(global_grid%vertex(global_grid%element(elem)%in)%gid,lgid,1)
      in_vert((i-1)*KEY_SIZE+1:i*KEY_SIZE) = lgid
      call hash_pack_key(global_grid%vertex(global_grid%element(elem)%out)%gid,lgid,1)
      out_vert((i-1)*KEY_SIZE+1:i*KEY_SIZE) = lgid
      last_elem = i
      elem = global_grid%element(elem)%next
   end do

   ierr = ZOLTAN_OK
endif

return
end subroutine z_coarse_obj_list

!        -----------
function z_num_child(dummy, nge, nle, gid, lid, ierr)
!        -----------

!----------------------------------------------------
! This routine returns the number of children of an element
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer(Zoltan_INT) :: z_num_child
integer(Zoltan_INT), intent(in) :: nge, nle
integer(Zoltan_INT), dimension(*), intent(in) :: dummy
integer(Zoltan_INT), intent(in) :: gid(*)
integer(Zoltan_INT), intent(in) :: lid(*)
integer(Zoltan_INT), intent(out) :: ierr

!----------------------------------------------------
! Local variables
integer :: loc_lid
type(hash_key) :: loc_gid
!----------------------------------------------------
! Begin executable code

! find the element

loc_gid = hash_unpack_key(gid(1:KEY_SIZE),1)
loc_lid = hash_decode_key(loc_gid,global_grid%elem_hash)
if (loc_lid == HASH_NOT_FOUND) then
   ierr = ZOLTAN_FATAL
   call fatal("failed to find element in z_num_child",procs=global_procs)
   return
endif

if (global_grid%element(loc_lid)%isleaf) then
   z_num_child = 0
else
! RESTRICTION bisection
   z_num_child = 2
endif
ierr = ZOLTAN_OK

return
end function z_num_child

!          ------------
subroutine z_child_list(dummy, nge, nle, pgid, plid, cgid, clid, mine, &
                        num_vert, verts, ref_type, in_vert, out_vert, ierr)
!          ------------

!----------------------------------------------------
! This routine returns the list of children of an element
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer(Zoltan_INT), dimension(*), intent(in) :: dummy
integer(Zoltan_INT), intent(in) :: nge, nle
integer(Zoltan_INT), intent(in) :: pgid(*)
integer(Zoltan_INT), intent(in) :: plid(*)
integer(Zoltan_INT), intent(out) :: cgid(*)
integer(Zoltan_INT), intent(out) :: clid(*), mine(*), num_vert(*), &
                                verts(*), in_vert(*), out_vert(*)
integer(Zoltan_INT), intent(out) :: ref_type, ierr

!----------------------------------------------------
! Local variables:

integer :: i, j, elem, lgid(KEY_SIZE), ind, loc_lid
integer :: child(MAX_CHILD), allc(MAX_CHILD)
type(hash_key) :: loc_gid

!----------------------------------------------------
! Begin executable code

! find the element

loc_gid = hash_unpack_key(pgid(1:KEY_SIZE),1)
loc_lid = hash_decode_key(loc_gid,global_grid%elem_hash)
if (loc_lid == HASH_NOT_FOUND) then
   ierr = ZOLTAN_FATAL
   call fatal("failed to find element in z_child_list",procs=global_procs)
   return
endif

allc = ALL_CHILDREN
child = get_child_lid(loc_gid,allc,global_grid%elem_hash)
if (child(1) == NO_CHILD) then
   ierr = ZOLTAN_FATAL
   call fatal("element in z_child_list has no children",procs=global_procs)
else
   ind = 1
! RESTRICTION bisection
   do i=1,2
      elem = child(i)
      call hash_pack_key(global_grid%element(elem)%gid,lgid,1)
      cgid((i-1)*KEY_SIZE+1:i*KEY_SIZE) = lgid
      clid(i) = 0
      if (global_grid%element(elem)%iown) then
         mine(i) = 1
      else
         mine(i) = 0
      endif
      num_vert(i) = VERTICES_PER_ELEMENT
      do j=1,VERTICES_PER_ELEMENT
         call hash_pack_key(global_grid%vertex(global_grid%element(elem)%vertex(j))%gid, &
                            verts(1:MAX_CHILD*KEY_SIZE),ind)
         ind = ind + KEY_SIZE
      end do
   end do
! RESTRICTION bisected triangles
   ref_type = ZOLTAN_TRI_BISECT
   ierr = ZOLTAN_OK
endif

return
end subroutine z_child_list

!          --------------
subroutine z_child_weight(data,nge,nle,global_id,local_id,wgt_dim,obj_wgt,ierr)
!          --------------

!----------------------------------------------------
! This routine returns the weight of the requested element
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer(Zoltan_INT), dimension(*), intent(in) :: data
integer(Zoltan_INT), intent(in) :: nge, nle
integer(Zoltan_INT), intent(in) :: global_id(*)
integer(Zoltan_INT), intent(in) :: local_id(*)
integer(Zoltan_INT), intent(in) :: wgt_dim
real(Zoltan_FLOAT), intent(out), dimension(*) :: obj_wgt
integer(Zoltan_INT), intent(out) :: ierr

!----------------------------------------------------
! Local variables

integer :: loc_lid
type(hash_key) :: loc_gid

!----------------------------------------------------
! Begin executable code

ierr = ZOLTAN_OK

! find the element

loc_gid = hash_unpack_key(global_id(1:KEY_SIZE),1)
loc_lid = hash_decode_key(loc_gid,global_grid%elem_hash)
if (loc_lid == HASH_NOT_FOUND) then
   ierr = ZOLTAN_FATAL
   call fatal("failed to find element in z_child_weight",procs=global_procs)
   return
endif

if (global_grid%element(loc_lid)%isleaf) then
   obj_wgt(1) = global_grid%element(loc_lid)%weight
else
   obj_wgt(1) = 0.0_lb_float
endif
obj_wgt(2:wgt_dim) = 0.0_lb_float

return
end subroutine z_child_weight

!        ----------
function z_partition(dummy, nge, nle, gid, lid, ierr)
!        ----------

!----------------------------------------------------
! This routine returns the current partition number of the given element
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer(Zoltan_INT) :: z_partition
integer(Zoltan_INT), intent(in) :: nge, nle
integer(Zoltan_INT), dimension(*), intent(in) :: dummy
integer(Zoltan_INT), intent(in) :: gid(*)
integer(Zoltan_INT), intent(in) :: lid(*)
integer(Zoltan_INT), intent(out) :: ierr

!----------------------------------------------------
! Local variables:

integer :: loc_lid
type(hash_key) :: loc_gid
!----------------------------------------------------
! Begin executable code

! For REFTREE I only need to get it right for elements I own.  If that
! doesn't hold for other methods, I might have to use get_elem_owner.

! on first call, pretend only processor 1 has the grid

if (first_time) then
   z_partition = 1

else
   loc_gid = hash_unpack_key(gid(1:KEY_SIZE),1)
   loc_lid = hash_decode_key(loc_gid,global_grid%elem_hash)

   if (loc_lid == HASH_NOT_FOUND) then
      z_partition = -1
   elseif (global_grid%element(loc_lid)%iown) then
      z_partition = my_proc(global_procs)
   else
      z_partition = -1
   endif

endif

! Zoltan is 0 based while PHAML is 1 based

z_partition = z_partition - 1

ierr = ZOLTAN_OK

end function z_partition

!          ----------------
subroutine zoltan_create_lb(lb,procs)
!          ----------------

!----------------------------------------------------
! This routine performs the setup of Zoltan for object lb
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(Zoltan_Struct), pointer :: lb
type(proc_info), intent(in), target :: procs

!----------------------------------------------------
! Local variables:

real(Zoltan_FLOAT) :: version
integer(Zoltan_INT) :: ierr, comm
character(len=1) :: gid_size
logical, save :: zoltan_initialized = .false.

!----------------------------------------------------
! Begin executable code

! no partitioning if 1 processor, and master does not participate

if (num_proc(procs)==1 .or. my_proc(procs)==MASTER) return

! don't initialize Zoltan more than once

if (.not. zoltan_initialized) then
   ierr = Zoltan_Initialize(version)
   if (ierr /= Zoltan_OK) then
      call fatal("Zoltan_Initialize returned error code",intlist=(/ierr/),procs=procs)
      stop
   endif
   zoltan_initialized = .true.
endif

! create the load balancing object

comm = slaves_communicator(procs)
lb => Zoltan_Create(comm)
if (.not. associated(lb)) then
   call fatal("Zoltan_Create failed",procs=procs)
   stop
endif

! set the Zoltan parameters that are independent of partitioning method

ierr = Zoltan_Set_Param(lb,"IMBALANCE_TOL","1.0")
ierr = Zoltan_Set_Param(lb,"OBJ_WEIGHT_DIM","1")
ierr = Zoltan_Set_Param(lb,"DEBUG_LEVEL","0")
if (KEY_SIZE > 9) then
! This is only because I am too lazy to deal with exceptions
! to the i1 format, and I doubt we would ever see keys so large
! that they need more than 9 integer words (OK, famous last words,
! but how big is 2**(9*32)?)
   call fatal("hash key size limited to 9 when using Zoltan",procs=procs)
   stop
endif
write(gid_size,"(i1)") KEY_SIZE
ierr = Zoltan_Set_Param(lb,"NUM_LID_ENTRIES","1")
ierr = Zoltan_Set_Param(lb,"NUM_GID_ENTRIES",gid_size)

ierr = Zoltan_Set_Num_Obj_Fn(lb,z_num_obj)
ierr = Zoltan_Set_Obj_List_Fn(lb,z_obj_list)
ierr = Zoltan_Set_Num_Geom_Fn(lb,z_num_geom)
ierr = Zoltan_Set_Geom_Fn(lb,z_geom)

return
end subroutine zoltan_create_lb

!          -----------------
subroutine zoltan_destroy_lb(lb,procs)
!          -----------------

!----------------------------------------------------
! This routine destroys Zoltan object lb
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(Zoltan_Struct), pointer :: lb
type(proc_info), intent(in), target :: procs

!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

! no partitioning if 1 processor, and master does not participate

if (num_proc(procs)==1 .or. my_proc(procs)==MASTER) return

call Zoltan_Destroy(lb)

return
end subroutine zoltan_destroy_lb

!          -----------
subroutine zoltan_init(lb,method,zoltan_param_file,procs)
!          -----------

!----------------------------------------------------
! This routine initializes Zoltan for the given method
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(Zoltan_Struct), pointer :: lb
integer, intent(in) :: method
character(len=*), intent(in) :: zoltan_param_file
type(proc_info), intent(in) :: procs

!----------------------------------------------------
! Local variables:

integer(Zoltan_INT) :: ierr, retval, comm

!----------------------------------------------------
! Begin executable code

if (num_proc(procs)==1 .or. my_proc(procs)==MASTER) return

select case(method)
case(ZOLTAN_RCB)
   ierr = Zoltan_Set_Param(lb,"LB_METHOD","RCB")
   ierr = Zoltan_Set_Param(lb,"CHECK_GEOM","0")
   ierr = Zoltan_Set_Param(lb,"RCB_OUTPUT_LEVEL","0")
case(ZOLTAN_OCT)
   ierr = Zoltan_Set_Param(lb,"LB_METHOD","OCTPART")
   ierr = Zoltan_Set_Param(lb,"OCT_DIM","2")
   ierr = Zoltan_Set_Param(lb,"OCT_OUTPUT_LEVEL","0")
case(ZOLTAN_METIS)
   ierr = Zoltan_Set_Num_Edges_Fn(lb,z_num_edge)
   ierr = Zoltan_Set_Edge_List_Fn(lb,z_edges)
   ierr = Zoltan_Set_Param(lb,"LB_METHOD","PARMETIS")
   ierr = Zoltan_Set_Param(lb,"CHECK_GRAPH","0")
   ierr = Zoltan_Set_Param(lb,"PARMETIS_OUTPUT_LEVEL","0")
   need_owner = .true.
case(ZOLTAN_REFTREE)
   ierr = Zoltan_Set_Num_Coarse_Obj_Fn(lb,z_num_coarse)
   ierr = Zoltan_Set_Coarse_Obj_List_Fn(lb,z_coarse_obj_list)
   ierr = Zoltan_Set_Num_Child_Fn(lb,z_num_child)
   ierr = Zoltan_Set_Child_List_Fn(lb,z_child_list)
   ierr = Zoltan_Set_Child_Weight_Fn(lb,z_child_weight)
   ierr = Zoltan_Set_Param(lb,"LB_METHOD","REFTREE")
case(ZOLTAN_RIB)
   ierr = Zoltan_Set_Param(lb,"LB_METHOD","RIB")
   ierr = Zoltan_Set_Param(lb,"CHECK_GEOM","0")
   ierr = Zoltan_Set_Param(lb,"RIB_OUTPUT_LEVEL","0")
case(ZOLTAN_HSFC)
   ierr = Zoltan_Set_Param(lb,"LB_METHOD","HSFC")
case(ZOLTAN_FILE)
   ierr = Zoltan_Set_Num_Coarse_Obj_Fn(lb,z_num_coarse)
   ierr = Zoltan_Set_Coarse_Obj_List_Fn(lb,z_coarse_obj_list)
   ierr = Zoltan_Set_Num_Child_Fn(lb,z_num_child)
   ierr = Zoltan_Set_Child_List_Fn(lb,z_child_list)
   ierr = Zoltan_Set_Child_Weight_Fn(lb,z_child_weight)
   ierr = Zoltan_Set_Num_Edges_Fn(lb,z_num_edge)
   ierr = Zoltan_Set_Edge_List_Fn(lb,z_edges)
   ierr = Zoltan_Set_Param(lb,"CHECK_GRAPH","0")
   comm = slaves_communicator(procs)
   ierr = zoltanparams_read_file(lb, trim(zoltan_param_file), comm)
   call zoltanparams_using_graph(retval)
   if (retval .eq. 1) then
      need_owner = .true.
   endif
end select

zoltan_method = method

end subroutine zoltan_init

!        ----------------------
function zoltanparams_read_file(lb,filename,communicator)
!        ----------------------

!----------------------------------------------------
! This routine call the Zoltan params library to read Zoltan parameters
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(Zoltan_Struct), pointer :: lb
character(len=*), intent(in) :: filename
integer(Zoltan_INT), intent(in) :: communicator
integer :: zoltanparams_read_file

!----------------------------------------------------
! Local variables:

integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: lb_addr
integer(Zoltan_INT), dimension(len(filename)) :: ifile
integer(Zoltan_INT) :: nbytes, nchar, i, ierr
!----------------------------------------------------
! Begin executable code

nbytes = Zoltan_PTR_LENGTH
call Zoltan_Get_Struct_Addr(lb,lb_addr)
nchar = len(filename)
do i=1,nchar
   ifile(i) = ichar(filename(i:i))
end do
call Zf90_zoltanparams_read_file(lb_addr,nbytes,ifile,nchar,communicator,ierr)
zoltanparams_read_file = ierr

end function zoltanparams_read_file

!          ----------------
subroutine zoltan_partition(grid,procs,lb,partmeth, &
                            export_gid,export_part,first_call)
!          ----------------

!----------------------------------------------------
! This routine uses Zoltan to partition the grid.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), target, intent(in) :: grid
type(proc_info), target, intent(in) :: procs
type(Zoltan_Struct), pointer :: lb
integer, intent(in) :: partmeth
type(hash_key), pointer :: export_gid(:)
integer, pointer :: export_part(:)
logical, optional :: first_call

!----------------------------------------------------
! Local variables:

integer(Zoltan_INT) :: ierr, num_imp, num_exp
integer(Zoltan_INT), pointer :: imp_gid(:),imp_lid(:),imp_proc(:),imp_part(:), &
                                exp_gid(:),exp_lid(:),exp_proc(:),exp_part(:)
integer :: i, astat, count, lid, my_part_num
logical :: changes
integer(Zoltan_INT) :: ngident, nlident
integer(Zoltan_INT), save :: last_npart = 0
type(hash_key), allocatable :: temp_gid(:)
integer, allocatable :: temp_part(:)
type(hash_key) :: loc_gid

!----------------------------------------------------
! Begin executable code

if (present(first_call)) then
   if (first_call) then
      first_time = .true.
   else
      first_time = .false.
   endif
else
   first_time = .false.
endif

! no partitioning if 1 processor, and master does not participate

if (num_proc(procs)==1 .or. my_proc(procs)==MASTER) return

! set pointers to data to pass to callback functions

global_grid => grid
global_procs => procs

! get the owner of each element, if needed

if (need_owner) then
   call get_elem_owner(grid,procs)
endif

! select the method for ParMETIS

if (partmeth == ZOLTAN_METIS .or. partmeth == ZOLTAN_FILE) then
   if (first_time) then
      ierr = Zoltan_Set_Param(lb,"PARMETIS_METHOD","PartKway")
   else
      ierr = Zoltan_Set_Param(lb,"PARMETIS_METHOD","RepartLDiffusion")
   endif
endif

nullify(imp_gid,imp_lid,imp_proc,imp_part,exp_gid,exp_lid,exp_proc,exp_part)

! call Zoltan

ierr = Zoltan_LB_Balance(lb,changes,ngident,nlident,num_imp,imp_gid,imp_lid, &
                         imp_proc,num_exp,exp_gid,exp_lid,exp_part)

if (need_owner) deallocate(elem_owner,stat=astat)

if (first_time) then

! the first time, only processor 1 computes the partition.  Send it to
! the others.

   call proc1_to_others(grid,procs,exp_gid,exp_part,num_exp)

endif

! copy results to phaml's arrays

if (num_exp > 0) then
   allocate(temp_gid(num_exp),temp_part(num_exp),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in zoltan_partition",procs=procs)
      return
   endif

   my_part_num = my_proc(procs)
   count = 0
   do i=1,num_exp
      loc_gid = hash_unpack_key(exp_gid,(i-1)*KEY_SIZE+1)
      lid = hash_decode_key(loc_gid,grid%elem_hash)
      if (lid /= HASH_NOT_FOUND) then
         if (grid%element(lid)%iown) then
            if (exp_part(i)+1 /= my_part_num) then
               count = count + 1
               temp_gid(count) = loc_gid
! Zoltan numbers partitions starting at 0, PHAML starts at 1
               temp_part(count) = exp_part(i) + 1
            endif
         endif
      endif
   end do
   if (count == 0) then
      nullify(export_gid, export_part)
   else
      allocate(export_gid(count),export_part(count),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in zoltan_partition",procs=procs)
         return
      endif
      export_gid = temp_gid(1:count)
      export_part = temp_part(1:count)
   endif
   deallocate(temp_gid,temp_part)
endif

! deallocate the arrays returned by Zoltan

if (first_time .and. my_proc(procs) /= 1) then
   deallocate(exp_gid,exp_part,stat=astat)
else
   ierr = Zoltan_LB_Free_Data(imp_gid,imp_lid,imp_proc,exp_gid,exp_lid,exp_part)
endif

! no longer flag it as the first call

first_time = .false.

end subroutine zoltan_partition

!          --------------
subroutine get_elem_owner(grid,procs)
!          --------------

!----------------------------------------------------
! This routine finds the owner of each element
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
type(proc_info), intent(in) :: procs
!----------------------------------------------------
! Local variables:

integer :: my_processor, nproc, count, lev, elem, part, i, lid, astat
! newcomm
integer :: nsend, scount, rcount, oscount
integer, allocatable :: nsendv(:), isend(:), nrecv(:)
integer, pointer :: irecv(:)
!----------------------------------------------------
! Begin executable code

my_processor = my_proc(procs)
nproc = num_proc(procs)

allocate(elem_owner(size(grid%element)),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in zoltan_partition",procs=procs)
   return
endif

allocate(nsendv(nproc), nrecv(nproc), stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in zoltan_partition",procs=procs)
   return
endif

! pass through the elements, counting the number that I don't own, and
! allocate space for a list of them, and space for the owner list

count = 0
do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (.not. grid%element(elem)%iown) count = count + 1
      elem = grid%element(elem)%next
   end do
end do

nsend = count*KEY_SIZE
allocate(isend(nsend),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in get_elem_owner",procs=procs)
   return
endif

elem_owner = 0

! pass through again, making a list of the ones I don't know and setting
! elem_owner for the ones I do

count = 1
do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (grid%element(elem)%iown) then
         elem_owner(elem) = my_processor
      else
         call hash_pack_key(grid%element(elem)%gid,isend,count)
         count = count + KEY_SIZE
      endif
      elem = grid%element(elem)%next
   end do
end do

! exchange lists with the other partitions

call phaml_alltoall(procs,isend,nsend,irecv,nrecv,610)

! scan through the list, counting the number that I own

count = 0
rcount = 1
do part=1,nproc
   do i=1,nrecv(part)/KEY_SIZE
      lid = hash_decode_key(hash_unpack_key(irecv,rcount),grid%elem_hash)
      if (lid /= HASH_NOT_FOUND) then
         if (grid%element(lid)%iown) count = count + 1
      endif
      rcount = rcount + KEY_SIZE
   end do
end do

! create the response with the list of gids I own

deallocate(isend,stat=astat)
allocate(isend(count*KEY_SIZE),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in get_elem_owner",procs=procs)
   return
endif
scount = 1
oscount = scount
rcount = 1
do part=1,nproc
   do i=1,nrecv(part)/KEY_SIZE
      lid = hash_decode_key(hash_unpack_key(irecv,rcount),grid%elem_hash)
      if (lid /= HASH_NOT_FOUND) then
         if (grid%element(lid)%iown) then
            isend(scount:scount+KEY_SIZE-1) = irecv(rcount:rcount+KEY_SIZE-1)
            scount = scount + KEY_SIZE
         endif
      endif
      rcount = rcount + KEY_SIZE
   end do
   nsendv(part) = scount - oscount
   oscount = scount
end do
if (associated(irecv)) deallocate(irecv,stat=astat)

! exchange responses

call phaml_alltoall(procs,isend,nsendv,irecv,nrecv,620)

! set elem_owner from the responses

count = 1
do part=1,nproc
   do i=1,nrecv(part)/KEY_SIZE
      lid = hash_decode_key(hash_unpack_key(irecv,count),grid%elem_hash)
      if (lid /= HASH_NOT_FOUND) then
         elem_owner(lid) = part
      endif
      count = count + KEY_SIZE
   end do
end do
if (associated(irecv)) deallocate(irecv,stat=astat)

deallocate(isend,nsendv,nrecv,stat=astat)

return
end subroutine get_elem_owner

!          ---------------
subroutine proc1_to_others(grid,procs,exp_gid,exp_proc,num_exp)
!          ---------------

!----------------------------------------------------
! This routine sends the partition from processor 1 to the other processors
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
type(proc_info), intent(in) :: procs
integer(Zoltan_INT), pointer :: exp_gid(:), exp_proc(:)
integer, intent(inout) :: num_exp

!----------------------------------------------------
! Local variables:

integer :: i, elem, next, ni, nr, my_processor, nproc, astat, loc_lid, ind
real(my_real) :: no_reals(1)
integer, allocatable :: export_to(:),send_int(:),num_kept(:)
integer, pointer :: recv_int(:)
real(my_real), pointer :: recv_real(:)
type(hash_key) :: loc_gid
!----------------------------------------------------
! Begin executable code

my_processor = my_proc(procs)
nproc = num_proc(procs)

! processor 1 builds the full list of exports and sends it to all other
! processors

if (my_processor == 1) then

! determine the new processor for each leaf.  This is needed because the
! export list doesn't tell which elements belong to processor 1

   allocate(export_to(size(grid%element)),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in proc1_to_others",procs=procs)
      return
   endif
   export_to = 0 ! Zoltan starts numbering at 0
   
   do i=1,num_exp
      loc_gid = hash_unpack_key(exp_gid,(i-1)*KEY_SIZE+1)
      loc_lid = hash_decode_key(loc_gid,grid%elem_hash)
      if (loc_lid == HASH_NOT_FOUND) then
         call warning("didn't find export element in proc1_to_others")
      else
         export_to(loc_lid) = exp_proc(i)
      endif
   end do

! copy the exports into a buffer for sending, by using a tree traversal

   allocate(send_int(grid%nelem_leaf*(KEY_SIZE+1)),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in proc1_to_others",procs=procs)
      return
   endif
   next = 1
   elem = grid%head_level_elem(1)
   do while (elem /= END_OF_LIST)
      call make_export_list(grid,export_to,elem,send_int,next)
      elem = grid%element(elem)%next
   end do

   deallocate(num_kept,export_to,stat=astat)

! send the buffer to the other processors

   do i=2,nproc
      call phaml_send(procs,i,send_int,next-1,no_reals,0,630)
   end do
   deallocate(send_int,stat=astat)

else

   if (associated(exp_gid)) deallocate(exp_gid)
   if (associated(exp_proc)) deallocate(exp_proc)

! processors other than processor 1 receive the export info from proc 1

   call phaml_recv(procs,i,recv_int,ni,recv_real,nr,630)

! determine how many exports there are and allocate the export list

   num_exp = ni/(KEY_SIZE+1)
   allocate(exp_gid(num_exp*KEY_SIZE),exp_proc(num_exp),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in proc1_to_others",procs=procs)
      return
   endif

! copy the exports from the message

   ind = 1
   do i=1,num_exp
      exp_gid((i-1)*KEY_SIZE+1:i*KEY_SIZE) = recv_int(ind:ind+KEY_SIZE-1)
      exp_proc(i) = recv_int(ind+KEY_SIZE)
      ind = ind + KEY_SIZE + 1
   end do

   deallocate(recv_int,stat=astat)

endif 

end subroutine proc1_to_others

!                    ----------------
recursive subroutine make_export_list(grid,export_to,elem,send_int,next)
!                    ----------------

!----------------------------------------------------
! This routine, called by processor 1, builds the export list to
! send to other processors
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: export_to(:)
integer, intent(in) :: elem
integer, intent(inout) :: send_int(:), next
!----------------------------------------------------
! Local variables:

integer :: children(MAX_CHILD), i, allc(MAX_CHILD)
!----------------------------------------------------
! Begin executable code

! if it is a leaf, add it to the list

allc = ALL_CHILDREN
children = get_child_lid(grid%element(elem)%gid,allc,grid%elem_hash)
if (children(1) == NO_CHILD) then
   call hash_pack_key(grid%element(elem)%gid,send_int,next)
   send_int(next+KEY_SIZE) = export_to(elem)
   next = next + KEY_SIZE + 1

! otherwise, traverse the children

else
   do i=1,MAX_CHILD
      if (children(i) /= NO_CHILD) then
         call make_export_list(grid,export_to,children(i),send_int,next)
      endif
   end do

endif

end subroutine make_export_list

!          ------------------
subroutine zoltan_child_order(order,ierr,lb,grid,procs,still_sequential)
!          ------------------
integer, intent(inout) :: order(:)
integer, intent(out) :: ierr
type(Zoltan_Struct), pointer :: lb
type(grid_type), target, intent(in) :: grid
type(proc_info), target, intent(in) :: procs
logical, intent(in) :: still_sequential

global_grid => grid
global_procs => procs
first_time = still_sequential

if (zoltan_method == ZOLTAN_REFTREE) then
   call Zoltan_Get_Child_Order(lb,order,ierr)
else
   ierr = Zoltan_FATAL
endif

end subroutine zoltan_child_order

end module zoltan_interf
