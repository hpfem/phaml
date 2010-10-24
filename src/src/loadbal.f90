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

module load_balance

!----------------------------------------------------
! This module contains routines for load balancing.
!
! communication tags in this module are of the form 5xx
!----------------------------------------------------

!----------------------------------------------------
! Other modules used are:

use global
use message_passing
use hash_mod
use gridtype_mod
use grid_util
use zoltan_interf
use stopwatch
use refine_elements
use error_estimators
use hp_strategies
use linsys_io

!----------------------------------------------------

implicit none
private
public partition, redistribute, set_weights, check_balance

!----------------------------------------------------
! The following parameters are defined:

! partition assignment for elements that are not in a single partition
integer, parameter :: NO_PARTITION = -1

!----------------------------------------------------

!----------------------------------------------------
! The following types are defined:

!----------------------------------------------------

!----------------------------------------------------
! The following variables are defined:

!----------------------------------------------------

contains

!          ---------
subroutine partition(grid,procs,lb,refcont,predictive_load_balance,form, &
                     balance_what,new_num_part,export_gid,export_proc, &
                     first_call)
!          ---------

!----------------------------------------------------
! This routine determines the partitioning of the grid.
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

type (grid_type), target, intent(inout) :: grid
type (proc_info), target, intent(in) :: procs
type (Zoltan_Struct), pointer :: lb
type (refine_options), intent(in) :: refcont
logical, intent(in) :: predictive_load_balance
integer, intent(in) :: form, balance_what, new_num_part
type(hash_key), pointer :: export_gid(:)
integer, pointer :: export_proc(:)
logical, optional :: first_call

!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

nullify(export_gid,export_proc)
if (my_proc(procs) == MASTER .or. PARALLEL == SEQUENTIAL) return

! start timing the partitioning process

call reset_watch((/ppartition,cppartition/))
call start_watch((/ppartition,tpartition/))

! set the weights

call set_weights(grid,predictive_load_balance,balance_what,refcont,procs)

! select method for partitioning.

select case(form)

   case(RTK)
      call reftree_kway(grid,procs,new_num_part,export_gid, &
                        export_proc,first_call)

   case(ZOLTAN_REFTREE, ZOLTAN_RCB, ZOLTAN_OCT, ZOLTAN_RIB, ZOLTAN_METIS, &
        ZOLTAN_HSFC, ZOLTAN_FILE)
      call zoltan_partition(grid,procs,lb,form,export_gid,export_proc, &
                            first_call)

   case default
      call fatal("unknown partition method specified",procs=procs)

end select

! stop timing the partitioning process

call stop_watch((/ppartition,tpartition/))

end subroutine partition

!          -----------
subroutine set_weights(grid,predictive,balance_what,refcont,procs)
!          -----------

!----------------------------------------------------
! This routine sets the weights for each element
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
logical, intent(in) :: predictive
integer, intent(in) :: balance_what
type (refine_options), intent(in) :: refcont
type(proc_info), intent(in) :: procs
!----------------------------------------------------
! Local variables:

integer :: assoc_verts(size(grid%element))
integer :: lev, vert, elem, elem_deg, side, edge, total
real(my_real) :: global_max_errind, reftol, normsoln, eta, eta_max, R, &
                 two_div_log2, loggamma
character(len=1) :: reftype(size(grid%element))
!----------------------------------------------------
! Begin executable code

! count the number of associated vertices for each element

assoc_verts = 0
do lev=1,grid%nlev
   vert = grid%head_level_vert(lev)
   do while (vert /= END_OF_LIST)
      assoc_verts(grid%vertex(vert)%assoc_elem) = &
         assoc_verts(grid%vertex(vert)%assoc_elem) + 1
      vert = grid%vertex(vert)%next
   end do
end do

! compute error indicators if needed

if (.not. grid%errind_up2date) then
   call all_error_indicators(grid,refcont%error_estimator)
endif

! compute maximum error indicator

global_max_errind = compute_global_max_errind(grid,procs)

! some constants for predictive load balancing

eta_max = global_max_errind
R = 100000.0_my_real
two_div_log2 = 2.0_my_real/log(2.0_my_real)
loggamma = log(0.125_my_real)

! for predictive load balancing, determine the refinement tolerance for the
! ONE_REF terminations

reftol = 0.0_my_real
if (predictive .and. refcont%refterm == ONE_REF) then
   call get_grid_info(grid,procs,.false.,3344, &
                      total_nelem_leaf=total,no_master=.true.)
   reftol = refcont%reftol/sqrt(real(total,my_real))
! TEMP only using first component of first eigenvalue
   if (grid%errtype == RELATIVE_ERROR .and. .not. (refcont%reftype == HP_ADAPTIVE .and. &
    (refcont%hp_strategy == HP_T3S .or. refcont%hp_strategy == HP_ALTERNATE))) then
      call norm_solution(grid,procs,.false.,1,1,energy=normsoln)
      if (normsoln /= 0.0_my_real) reftol = reftol*normsoln
   endif
endif
if (predictive .and. refcont%refterm == ONE_REF_HALF_ERRIND) then
   reftol = global_max_errind/refcont%inc_factor
endif

! for each element ...

do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (grid%element(elem)%isleaf .and. grid%element(elem)%iown) then

         elem_deg = grid%element(elem)%degree

         if (predictive) then

! for predictive load balancing, first make a guess as to whether or not this
! element will be refined and whether it will be refined by h or p

! For ONE_REF terminations, check the tolerance to determine if it will be
! refined.  If so, or if it is not a ONE_REF termination, determine h or p.
! TEMP this is a bit expensive because we are determining the h or p
!      refinement choice twice

            if ((refcont%refterm /= ONE_REF .and. &
                 refcont%refterm /= ONE_REF_HALF_ERRIND) .or. &
                 maxval(grid%element_errind(elem,:)) > reftol) then
               call mark_reftype_one(elem,grid,refcont,global_max_errind, &
                                     reftype,.false.)
            else
               reftype(elem) = "n"
            endif

            select case (reftype(elem))

            case ("n")

! if it will not be refined, just assign the weight to be the number of
! whatever we are balancing on associated with this element

               select case (balance_what)

               case (BALANCE_NONE)
                  grid%element(elem)%weight = 0.0_my_real

               case (BALANCE_ELEMENTS)
                  grid%element(elem)%weight = 1.0_my_real

               case (BALANCE_VERTICES)
                  grid%element(elem)%weight = assoc_verts(elem)

               case (BALANCE_EQUATIONS)
                  grid%element(elem)%weight = assoc_verts(elem) + &
                                       ((elem_deg-1)*(elem_deg-2))/2.0_my_real
                  do side=1,3
                     edge = grid%element(elem)%edge(side)
                     if (grid%edge(edge)%assoc_elem == elem) then
                        grid%element(elem)%weight = grid%element(elem)%weight +&
                           grid%edge(edge)%degree - 1
                     endif
                  end do

               case default
                  ierr = PHAML_INTERNAL_ERROR
                  call fatal("unknown value for balance_what in set_weights")
                  stop

               end select

            case ("h")

! if it will be refined by h, assign the number of entities there should be
! associated with the children

               select case (balance_what)

               case (BALANCE_NONE)
                  grid%element(elem)%weight = 0.0_my_real

               case (BALANCE_ELEMENTS)
                  grid%element(elem)%weight = 2.0_my_real

               case (BALANCE_VERTICES)
                  grid%element(elem)%weight = assoc_verts(elem) + 1

               case (BALANCE_EQUATIONS)
! the equations are: the associated vertices, the new vertex (guess that it
! will be associated, could be wrong), the face equations, the equations on
! the new edge, and the edge equations for associated edges (note the base
! will be refined and has twice as many, and guess they will be associated)
                  grid%element(elem)%weight = assoc_verts(elem)+1 + &
                                              (elem_deg-1)*(elem_deg-2) + &
                                              elem_deg-1
                  do side=1,3
                     edge = grid%element(elem)%edge(side)
                     if (grid%edge(edge)%assoc_elem == elem) then
                        if (side == 3) then
                           grid%element(elem)%weight = &
                                                   grid%element(elem)%weight + &
                                                   2*(grid%edge(edge)%degree-1)
                        else
                           grid%element(elem)%weight = &
                                                   grid%element(elem)%weight + &
                                                   grid%edge(edge)%degree - 1
                        endif
                     endif
                  end do

               case default
                  ierr = PHAML_INTERNAL_ERROR
                  call fatal("unknown value for balance_what in set_weights")
                  stop

               end select

            case ("p")

! if it will be refined by p, assign the number of entities there should be
! associated after p refinement

               select case (balance_what)

               case (BALANCE_NONE)
                  grid%element(elem)%weight = 0.0_my_real

               case (BALANCE_ELEMENTS)
                  grid%element(elem)%weight = 1.0_my_real

               case (BALANCE_VERTICES)
                  grid%element(elem)%weight = assoc_verts(elem)

               case (BALANCE_EQUATIONS)
! for associated edges, guess the edge degree is the same as this element
                  grid%element(elem)%weight = assoc_verts(elem) + &
                                              (elem_deg*(elem_deg-1))/2
                  do side=1,3
                     edge = grid%element(elem)%edge(side)
                     if (grid%edge(edge)%assoc_elem == elem) then
                        grid%element(elem)%weight = grid%element(elem)%weight +&
                                                    elem_deg
                     endif
                  end do

               case default
                  ierr = PHAML_INTERNAL_ERROR
                  call fatal("unknown value for balance_what in set_weights")
                  stop

               end select

            case default
               ierr = PHAML_INTERNAL_ERROR
               call fatal("unknown value for reftype in set_weights")
               stop

            end select

! For terminations that are not ONE_REF, multiply by a factor that will be
! larger in elements we expect to be refined more times.  For simplicity,
! we just multiply by an estimate of how many times we think the element
! would have to be refined to reach a certain error, rather than properly
! estimating the number of entities there would be after that many refinements.
! Given: the error estimate for the element, eta
!        the maximum error estimate over all elements, eta_max
!        some artificial factor by which to reduce the maximum error, R
!        a factor by which we believe p refinement reduces the error, gamma
!        the current element degree, p
! the number of refinements needed to reduce the error to eta_max/R is
! ceiling(A log(R*eta/eta_max)) where A is 2/(p*log(2)) for h refinement
! and -log(gamma) for p refinement.
! If eta is already less than eta_max/R, don't multiply by the factor.
! And don't worry about the ceiling.

            if (refcont%refterm /= ONE_REF .and. &
                 refcont%refterm /= ONE_REF_HALF_ERRIND) then
               eta = maxval(grid%element_errind(elem,:))
               if (eta > eta_max/R) then
                  if (reftype(elem) == "h") then
                     grid%element(elem)%weight = grid%element(elem)%weight * &
                      two_div_log2*log(eta*R/eta_max)/grid%element(elem)%degree
                  elseif (reftype(elem) == "p") then
                     grid%element(elem)%weight = grid%element(elem)%weight * &
                      (-loggamma)*log(eta*R/eta_max)
                  endif
               endif
            endif

         else ! not predictive

! if it is not predictive load balancing, just assign the weight to be
! the number of whatever we are balancing on associated with this element

            select case (balance_what)

            case (BALANCE_NONE)
               grid%element(elem)%weight = 0.0_my_real

            case (BALANCE_ELEMENTS)
               grid%element(elem)%weight = 1.0_my_real

            case (BALANCE_VERTICES)
               grid%element(elem)%weight = assoc_verts(elem)

            case (BALANCE_EQUATIONS)
               grid%element(elem)%weight = assoc_verts(elem) + &
                                       ((elem_deg-1)*(elem_deg-2))/2.0_my_real
               do side=1,3
                  edge = grid%element(elem)%edge(side)
                  if (grid%edge(edge)%assoc_elem == elem) then
                     grid%element(elem)%weight = grid%element(elem)%weight + &
                        grid%edge(edge)%degree - 1
                  endif
               end do

            case default
               ierr = PHAML_INTERNAL_ERROR
               call fatal("unknown value for balance_what in set_weights")
               stop

            end select

         endif ! predictive

      else ! not a leaf or I don't own it

         grid%element(elem)%weight = 0.0_my_real

      endif

      elem = grid%element(elem)%next
   end do
end do

end subroutine set_weights

!          ------------
subroutine reftree_kway(grid,procs,new_num_part,export_gid, &
                        export_part,first_call)
!          ------------

!----------------------------------------------------
! This routine partitions the grid using the k-way refinement tree method
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type (grid_type), intent(in) :: grid
type (proc_info), intent(in) :: procs
integer, intent(in) :: new_num_part
type(hash_key), pointer :: export_gid(:)
integer, pointer :: export_part(:)
logical, optional :: first_call
!----------------------------------------------------
! Local variables:

real(my_real) :: root_weight
real(my_real), allocatable :: weight(:)
integer :: elem, i, my_lid, part, num_exp, allocstat, count, nproc, &
           my_processor, issub, rssub, oissub, orssub, irsub, rrsub
real(my_real) :: part_size, current_size, cutoff
! newcomm
integer :: nisend
integer, allocatable :: isend(:),nisendv(:),nirecv(:),nrsendv(:),nrrecv(:)
integer, pointer :: irecv(:)
real(my_real), allocatable :: rsend(:)
real(my_real), pointer :: rrecv(:)
type(hash_key), allocatable :: leaf_list(:)
! 1 if all descendents are mine, 0 if no descendents are mine, -1 if some
! descendents are mine
integer, allocatable :: mine(:)
integer, allocatable :: new_partition(:)
logical :: skip_comm
!----------------------------------------------------
! Begin executable code

nproc = num_proc(procs)
my_processor = my_proc(procs)

allocate(nisendv(nproc),nirecv(nproc),nrsendv(nproc),nrrecv(nproc), &
         stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in reftree_kway",procs=procs)
   return
endif

allocate(mine(size(grid%element)),weight(size(grid%element)),stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in reftree_kway",procs=procs)
   return
endif

! sum the weights of the nodes that are assigned to this partition,
! and count the leaves that are not

count = 0
elem = grid%head_level_elem(1)
do while (elem /= END_OF_LIST)
   call sum_my_weights(grid,elem,weight,mine,count)
   elem = grid%element(elem)%next
end do

! on first call, all nodes are assigned to all partitions, so skip the
! communication part

if (present(first_call)) then
   if (first_call) then
      if (count == 0) then
         skip_comm = .true.
      else
         call warning("On first call to reftree_kway, found leaves not assigned to this processor.", &
                      "Partitions may be wrong.")
         skip_comm = .true.
      endif
   else
      skip_comm = .false.
   endif
else
   skip_comm = .false.
endif

if (.not. skip_comm) then

! make a list of the leaves that are not assigned to each partition
! TEMP may be more efficient to place them directly into isend, but it
!      is clearer to assign them into an array of gid's

   allocate(leaf_list(count),stat=allocstat)
   if (allocstat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in reftree_kway",procs=procs)
      return
   endif
   count = 0
   elem = grid%head_level_elem(1)
   do while (elem /= END_OF_LIST)
      call list_other_leaves(grid,elem,leaf_list,count)
      elem = grid%element(elem)%next
   end do

! get unknown leaf weights from other partitions

! copy leaf list to isend

   nisend = size(leaf_list)*KEY_SIZE
   allocate(isend(nisend),stat=allocstat)
   if (allocstat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in reftree_kway",procs=procs)
      return
   endif
   do i=1,size(leaf_list)
      call hash_pack_key(leaf_list(i),isend,1+(i-1)*KEY_SIZE)
   end do

! send my list of leaves for which I need partial sums, and get request list
! of other partitions

   call start_watch((/cppartition,ctpartition/))
   call phaml_alltoall(procs,isend,nisend,irecv,nirecv,510)
   call stop_watch((/cppartition,ctpartition/))

   deallocate(leaf_list,stat=allocstat)
   deallocate(isend,stat=allocstat)

! count the number of requests for which I have an answer, and allocate
! space for the next message

   count = 0
   irsub = 1
   do part=1,nproc
      do i=1,nirecv(part)/KEY_SIZE
         my_lid = hash_decode_key(hash_unpack_key(irecv,irsub),grid%elem_hash)
         if (my_lid /= HASH_NOT_FOUND) then
            if (weight(my_lid) /= 0.0_my_real) then
               count = count + 1
            endif
         endif
         irsub = irsub + KEY_SIZE
      end do
   end do

   allocate(isend(count*KEY_SIZE),rsend(count),stat=allocstat)
   if (allocstat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in reftree_kway",procs=procs)
      return
   endif

! build my reply to the requests

   irsub = 1
   issub = 1
   rssub = 1
   oissub = issub
   orssub = rssub
   do part=1,nproc
      do i=1,nirecv(part)/KEY_SIZE
         my_lid = hash_decode_key(hash_unpack_key(irecv,irsub),grid%elem_hash)
         if (my_lid /= HASH_NOT_FOUND) then
            if (weight(my_lid) /= 0.0_my_real) then
               isend(issub:issub+KEY_SIZE-1) = irecv(irsub:irsub+KEY_SIZE-1)
               issub = issub + KEY_SIZE
               rsend(rssub) = weight(my_lid)
               rssub = rssub + 1
            endif
         endif
         irsub = irsub + KEY_SIZE
      end do
      nisendv(part) = issub - oissub
      nrsendv(part) = rssub - orssub
      oissub = issub
      orssub = rssub
   end do

   if (associated(irecv)) deallocate(irecv,stat=allocstat)

! send the reply

   call start_watch((/cppartition,ctpartition/))
   call phaml_alltoall(procs,isend,nisendv,irecv,nirecv,520)
   call phaml_alltoall(procs,rsend,nrsendv,rrecv,nrrecv,521)
   call stop_watch((/cppartition,ctpartition/))
   deallocate(isend,rsend,stat=allocstat)

! add the replies into my weights

   irsub = 1
   rrsub = 1
   do part=1,nproc
      do i=1,nrrecv(part)
         my_lid = hash_decode_key(hash_unpack_key(irecv,irsub),grid%elem_hash)
         if (my_lid /= HASH_NOT_FOUND) then
            weight(my_lid) = weight(my_lid) + rrecv(rrsub)
         endif
         irsub = irsub + KEY_SIZE
         rrsub = rrsub + 1
      end do
   end do

! finish summing the weights

   elem = grid%head_level_elem(1)
   do while (elem /= END_OF_LIST)
      call sum_all_weights(grid,elem,weight,mine)
      elem = grid%element(elem)%next
   end do

   if (associated(irecv)) deallocate(irecv,stat=allocstat)
   if (associated(rrecv)) deallocate(rrecv,stat=allocstat)
   deallocate(nisendv,nirecv,nrsendv,nrrecv,stat=allocstat)

endif

! compute the root weight

root_weight = 0.0_my_real

elem = grid%head_level_elem(1)
do while (elem /= END_OF_LIST)
   root_weight = root_weight + weight(elem)
   elem = grid%element(elem)%next
end do

! determine partition sizes

part_size = root_weight/new_num_part

allocate(new_partition(size(grid%element)),stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in reftree_kway",procs=procs)
   return
endif

! traverse the tree to define partition and count the number of exports

num_exp = 0
part = 1
current_size = 0.0_my_real
cutoff = part_size

elem = grid%head_level_elem(1)
do while (elem /= END_OF_LIST)
   call define_partition(grid,elem,new_partition,weight, &
                         num_exp,part,current_size,cutoff,part_size, &
                         new_num_part,mine)
   elem = grid%element(elem)%next
end do

! if no exports, we're done

if (num_exp == 0) then
   nullify(export_gid, export_part)
else

! allocate space for the export lists

   allocate(export_gid(num_exp), export_part(num_exp), stat=allocstat)
   if (allocstat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in reftree_kway",procs=procs)
      return
   endif

! traverse the tree to make the export lists

   elem = grid%head_level_elem(1)
   count = 0
   do while (elem /= END_OF_LIST)
     call make_export_lists(grid,elem,mine,new_partition,export_gid, &
                            export_part,count)
     elem = grid%element(elem)%next
   end do

endif

deallocate(new_partition,weight,mine,stat=allocstat)

end subroutine reftree_kway

!                    --------------
recursive subroutine sum_my_weights(grid,subroot,weight,mine,count)
!                    --------------

!----------------------------------------------------
! This routine sums the weights of the elements owned by this partition, and
! counts the number of unowned leaves
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: subroot
real(my_real), intent(inout) :: weight(:)
integer, intent(inout) :: mine(:)
integer, intent(inout) :: count
!----------------------------------------------------
! Local variables:

integer :: children(MAX_CHILD), i, allc(MAX_CHILD)
logical :: none_assigned, all_assigned
!----------------------------------------------------
! Begin executable code

allc = ALL_CHILDREN
children = get_child_lid(grid%element(subroot)%gid,allc,grid%elem_hash)

if (children(1) == NO_CHILD) then

! if there are no children, then the sum is the weight of this node if it
! is assigned to this processor and 0.0 otherwise

   if (grid%element(subroot)%iown) then

      weight(subroot) = grid%element(subroot)%weight
      mine(subroot) = 1
   else
      weight(subroot) = 0.0_my_real
      mine(subroot) = 0
      count = count + 1
   endif

else

! if there are children, sum the weights of the children

! RESTRICTION no internal weights
   weight(subroot) = 0.0_my_real
   none_assigned = .true.
   all_assigned = .true.
   do i=1,MAX_CHILD
      if (children(i) /= NO_CHILD) then
         call sum_my_weights(grid,children(i),weight,mine,count)
         weight(subroot) = weight(subroot) + weight(children(i))
         if (mine(children(i))==1 .or. mine(children(i))==-1) then
            none_assigned = .false.
         endif
         if (mine(children(i))==0 .or. mine(children(i))==-1) then
            all_assigned = .false.
         endif
      endif
   end do
   if (none_assigned) then
      mine(subroot) = 0
   elseif (all_assigned) then
      mine(subroot) = 1
   else
      mine(subroot) = -1
   endif

endif

end subroutine sum_my_weights

!                    -----------------
recursive subroutine list_other_leaves(grid,subroot,leaf_list,count)
!                    -----------------

!----------------------------------------------------
! This routine creates a list of all the leaves that are not owned by this
! partition
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: subroot
type(hash_key), intent(inout) :: leaf_list(:)
integer, intent(inout) :: count
!----------------------------------------------------
! Local variables:

integer :: children(MAX_CHILD), i, allc(MAX_CHILD)
!----------------------------------------------------
! Begin executable code

allc = ALL_CHILDREN
children = get_child_lid(grid%element(subroot)%gid,allc,grid%elem_hash)

if (children(1) == NO_CHILD) then

   if (.not. grid%element(subroot)%iown) then
      count = count + 1
      leaf_list(count) = grid%element(subroot)%gid
   endif

else

   do i=1,MAX_CHILD
      if (children(i) /= NO_CHILD) then
         call list_other_leaves(grid,children(i),leaf_list,count)
      endif
   end do

endif

end subroutine list_other_leaves

!                    ---------------
recursive subroutine sum_all_weights(grid,subroot,weight,mine)
!                    ---------------

!----------------------------------------------------
! This routine sums the weights of all elements.  Elements for which mine is
! true already have the final sum; unowned leaves have received the sum from
! other partitions.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: subroot
real(my_real), intent(inout) :: weight(:)
integer, intent(inout) :: mine(:)
!----------------------------------------------------
! Local variables:

integer :: children(MAX_CHILD), i, allc(MAX_CHILD)
!----------------------------------------------------
! Begin executable code

allc = ALL_CHILDREN
children = get_child_lid(grid%element(subroot)%gid,allc,grid%elem_hash)

! if no children or I own all descendents, then done

if (children(1) /= NO_CHILD .and. mine(subroot) /= 1) then

! if there are children, sum the weights of the children

! RESTRICTION no internal weights
   weight(subroot) = 0.0_my_real
   do i=1,MAX_CHILD
      if (children(i) /= NO_CHILD) then
         call sum_all_weights(grid,children(i),weight,mine)
         weight(subroot) = weight(subroot) + weight(children(i))
      endif
   end do

endif

end subroutine sum_all_weights

!                    ----------------
recursive subroutine define_partition(grid,subroot,new_partition, &
                                      weight,num_exp,part,current_size,cutoff, &
                                      partition_size,new_num_part,mine)
!                    ----------------

!----------------------------------------------------
! This routine defines the partition and counts the number of exports
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type (grid_type), intent(in) :: grid
integer, intent(in) :: subroot
integer, intent(inout) :: new_partition(:)
real(my_real), intent(in) :: weight(:)
integer, intent(inout) :: num_exp, part
real(my_real), intent(inout) :: current_size, cutoff
real(my_real), intent(in) :: partition_size
integer, intent(in) :: new_num_part,mine(:)

!----------------------------------------------------
! Local variables:

real(my_real) :: newsize
integer :: children(MAX_CHILD), i, allc(MAX_CHILD)

!----------------------------------------------------
! Begin executable code

newsize = current_size + weight(subroot)
allc = ALL_CHILDREN
children = get_child_lid(grid%element(subroot)%gid,allc,grid%elem_hash)

if (newsize <= cutoff .or. part == new_num_part .or. &
    abs(newsize-cutoff) < 100*epsilon(0.0_my_real)) then

! this subtree fits in the current partition

   new_partition(subroot) = part
   current_size = newsize

! If this is this processor's partition, there are no exports below this
! node, so we don't have to traverse the subtree.
! If there are no leaves of this subtree assigned to this processor, there
! are no exports below this node.
! Otherwise, traverse the subtree setting partition and counting exports

   if (part /= grid%partition .and. mine(subroot) /= 0) then
      if (children(1) == NO_CHILD) then
         num_exp = num_exp + 1
      else
         do i=1,MAX_CHILD
           call mark_and_count(grid,children(i),part,mine,new_partition,num_exp)
         end do
      endif
   endif

else

! this subtree is too big for the current partition

   if (children(1) /= NO_CHILD) then

! if there are children, traverse them

      new_partition(subroot) = NO_PARTITION
      do i=1,MAX_CHILD
         call define_partition(grid, &
                               children(grid%element(subroot)%order(i)), &
                               new_partition,weight,num_exp,part, &
                               current_size,cutoff,partition_size, &
                               new_num_part,mine)
      end do

! if there are no children, move on to next partition

   else

      do while (newsize > cutoff)
         part = part + 1
         cutoff = part*partition_size
      end do
      new_partition(subroot) = part
      current_size = newsize
      if (part /= grid%partition .and. mine(subroot) /= 0) then
         num_exp = num_exp + 1
      endif
   endif
endif

end subroutine define_partition

!                    --------------
recursive subroutine mark_and_count(grid,subroot,part,mine,new_partition, &
                                    num_exp)
!                    --------------

!----------------------------------------------------
! This routine defines the new partition for elements in a subtree assigned
! to a different partition, and counts the exports
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: subroot,part,mine(:)
integer, intent(inout) :: new_partition(:),num_exp
!----------------------------------------------------
! Local variables:

integer :: children(MAX_CHILD), i, allc(MAX_CHILD)
!----------------------------------------------------
! Begin executable code

new_partition(subroot) = part
allc = ALL_CHILDREN
children = get_child_lid(grid%element(subroot)%gid,allc,grid%elem_hash)

if (children(1) == NO_CHILD) then
   if (mine(subroot) /= 0) num_exp = num_exp + 1
else
   do i=1,MAX_CHILD
      if (children(i) /= NO_CHILD) then
         call mark_and_count(grid,children(i),part,mine,new_partition,num_exp)
      endif
   end do

endif

end subroutine mark_and_count

!                    -----------------
recursive subroutine make_export_lists(grid,subroot,mine,new_partition, &
                                       export_gid,export_proc,count)
!                    -----------------

!----------------------------------------------------
! This routine creates the exports lists
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: subroot,mine(:),new_partition(:)
type(hash_key), intent(inout) :: export_gid(:)
integer, intent(inout) :: export_proc(:),count

!----------------------------------------------------
! Local variables:

integer :: children(MAX_CHILD), i, allc(MAX_CHILD)
!----------------------------------------------------
! Begin executable code

! if this subtree has no leaves assigned to this processor, or if the
! new partition of this subtree is this processor, then there can be no
! exports below it

if (mine(subroot) /= 0 .and. new_partition(subroot) /= grid%partition) then
   allc = ALL_CHILDREN
   children = get_child_lid(grid%element(subroot)%gid,allc,grid%elem_hash)
   if (children(1) == NO_CHILD) then

! if this is a leaf, put it on the export lists

      count = count + 1
      export_gid(count) = grid%element(subroot)%gid
      export_proc(count) = new_partition(subroot)

   else

! if it is not a leaf, traverse the children

      do i=1,MAX_CHILD
         call make_export_lists(grid,children(i),mine,new_partition, &
                                export_gid,export_proc,count)
      end do
   endif
endif

end subroutine make_export_lists

!          ------------
subroutine redistribute(grid,procs,refine_control,export_gid,export_part, &
                        first_call)
!          ------------

!----------------------------------------------------
! This routine redistributes the data structures over the partitions
! after 'partition' determines a new distribution.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type (grid_type), intent(inout) :: grid
type (proc_info), intent(in) :: procs
type(refine_options), intent(in) :: refine_control
type(hash_key), pointer :: export_gid(:)
integer, pointer :: export_part(:)
logical, intent(in), optional :: first_call
!----------------------------------------------------

!----------------------------------------------------
! Local variables:

integer, allocatable :: isend(:), nisend(:), nrsend(:), nirecv(:), nrrecv(:), &
                        iind(:), rind(:)
integer, pointer :: irecv(:)
real(my_real), allocatable :: rsend(:)
real(my_real), pointer :: rrecv(:)
type(hash_key) :: elemgid, gid
integer i, j, k, l, m, part, elemlid, errcode, allocstat, nproc, &
        iindex, rindex, d, to, elem, edge, vert, d1, d2, d3
logical :: skip_comm

!----------------------------------------------------
! Begin executable code

! no redistribution if not running in parallel

if (PARALLEL == SEQUENTIAL) return

grid%errind_up2date = .false.
grid_changed = .true.

! master does not participate

if (my_proc(procs) == MASTER) return

! start timing the redistribution process

call reset_watch((/pdistribute,cpdistribute/))
call start_watch((/pdistribute,tdistribute/))

nproc = num_proc(procs)

if (.not.associated(export_gid)) then
   allocate(export_gid(0),export_part(0),stat=allocstat)
   if (allocstat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in redistribute",procs=procs)
      return
   endif
endif

! remove myself as owner of the elements to be transfered

do i=1,size(export_gid)
   elemlid = hash_decode_key(export_gid(i),grid%elem_hash)
   if (elemlid == HASH_NOT_FOUND) then
      call warning("export list has a gid I don't have; skipping")
   else
      grid%element(elemlid)%iown = .false.
      call disown_parent(grid,export_gid(i))
      if (grid%element(elemlid)%degree >= 3) then
         grid%dof_own = grid%dof_own - grid%system_size * &
              ((grid%element(elemlid)%degree-2)*(grid%element(elemlid)%degree-1))/2
      endif
      do j=1,VERTICES_PER_ELEMENT
         if (grid%vertex(grid%element(elemlid)%vertex(j))%assoc_elem == elemlid) then
            grid%nvert_own = grid%nvert_own - 1
            grid%dof_own = grid%dof_own - grid%system_size
         endif
      end do
      do j=1,EDGES_PER_ELEMENT
         if (grid%edge(grid%element(elemlid)%edge(j))%assoc_elem == elemlid) then
            grid%nedge_own = grid%nedge_own - 1
            if (grid%edge(grid%element(elemlid)%edge(j))%degree >= 2) then
               grid%dof_own = grid%dof_own - grid%system_size * &
                          (grid%edge(grid%element(elemlid)%edge(j))%degree - 1)
            endif
         endif
      end do
   endif
end do
grid%nelem_leaf_own = grid%nelem_leaf_own - size(export_gid)

! on first call, all processors have all elements, so do not need to
! send the exported elements

if (present(first_call)) then
   if (first_call) then
      skip_comm = .true.
   else
      skip_comm = .false.
   endif
else
   skip_comm = .false.
endif

if (.not. skip_comm) then

   allocate(nisend(nproc),nrsend(nproc),nirecv(nproc),nrrecv(nproc), &
            iind(nproc),rind(nproc),stat=allocstat)
   if (allocstat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in redistribute",procs=procs)
      return
   endif

! pass through the export list to count how much will be sent to each
! processor.  see comments in the pass that copies the data

   nisend = 0
   nrsend = 0

   do i=1,size(export_gid)
      elemgid = export_gid(i)
      elemlid = hash_decode_key(elemgid,grid%elem_hash)
      to = export_part(i)
      nisend(to) = nisend(to) + KEY_SIZE
      nisend(to) = nisend(to) + 4
      nrsend(to) = nrsend(to) + VERTICES_PER_ELEMENT*grid%nsoln
      do j=1,EDGES_PER_ELEMENT
         d = grid%edge(grid%element(elemlid)%edge(j))%degree
         if (d > 1) nrsend(to) = nrsend(to) + grid%nsoln*(d-1)
      end do
      d = grid%element(elemlid)%degree
      if (d > 2) nrsend(to) = nrsend(to) + grid%nsoln*((d-1)*(d-2))/2
      if (grid%oldsoln_exists) then
         gid = elemgid
         do
            elem = hash_decode_key(gid,grid%elem_hash)
            if (grid%element(elem)%oldleaf) exit
            if (grid%element(elem)%level == 1) then
               ierr = PHAML_INTERNAL_ERROR
               call fatal("failed to find element with old solution")
               stop
            endif
            gid = gid/2
         end do
         nisend(to) = nisend(to) + KEY_SIZE
         nisend(to) = nisend(to) + 3
         if (associated(grid%element(elem)%oldsoln)) then
            nrsend(to) = nrsend(to) + size(grid%element(elem)%oldsoln)
         endif
         do j=1,EDGES_PER_ELEMENT
            edge = grid%element(elem)%edge(j)
            nisend(to) = nisend(to) + 3
            if (associated(grid%edge(edge)%oldsoln)) then
               nrsend(to) = nrsend(to) + size(grid%edge(edge)%oldsoln)
            endif
         end do
         do j=1,VERTICES_PER_ELEMENT
            vert = grid%element(elem)%vertex(j)
            nisend(to) = nisend(to) + 2
            if (associated(grid%vertex_oldsoln)) then
               nrsend(to) = nrsend(to) + size(grid%vertex_oldsoln,dim=2) * &
                                         size(grid%vertex_oldsoln,dim=3)
            endif
         end do
      endif
   end do ! next export gid

! allocate messages

   allocate(isend(sum(nisend)),rsend(sum(nrsend)),stat=allocstat)
   if (allocstat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in redistribute",procs=procs)
      return
   endif

! start the index counters at the beginning of each segment

   iind(1) = 1
   rind(1) = 1
   do i=2,nproc
      iind(i) = iind(i-1) + nisend(i-1)
      rind(i) = rind(i-1) + nrsend(i-1)
   end do

! for each element on the export list ...

   do i=1,size(export_gid)
      elemgid = export_gid(i)
      elemlid = hash_decode_key(elemgid,grid%elem_hash)
      to = export_part(i)

! pack the gid

      call hash_pack_key(elemgid,isend,iind(to))
      iind(to) = iind(to) + KEY_SIZE

! pack the edge and element degrees

      isend(iind(to):iind(to)+3)=(/grid%edge(grid%element(elemlid)%edge)%degree, &
                                   grid%element(elemlid)%degree /)
      iind(to) = iind(to) + 4

! pack the vertex solutions

      do j=1,VERTICES_PER_ELEMENT
         rsend(rind(to):rind(to)+grid%nsoln-1) = &
            reshape(grid%vertex_solution(grid%element(elemlid)%vertex(j),:,:), &
                    (/grid%nsoln/))
         rind(to) = rind(to) + grid%nsoln
      end do

! pack the edge solutions

      do j=1,EDGES_PER_ELEMENT
         d = grid%edge(grid%element(elemlid)%edge(j))%degree
         if (d > 1) then
            do k=1,grid%system_size
               do l=1,max(1,grid%num_eval)
                  rsend(rind(to):rind(to)+d-2) = &
                    grid%edge(grid%element(elemlid)%edge(j))%solution(1:d-1,k,l)
                  rind(to) = rind(to) + d-1
               end do
            end do
         endif
      end do

! pack the element solutions

      d = grid%element(elemlid)%degree
      if (d > 2) then
         do j=1,grid%system_size
            do k=1,max(1,grid%num_eval)
               rsend(rind(to):rind(to)+((d-1)*(d-2))/2-1) = &
                   grid%element(elemlid)%solution(1:((d-1)*(d-2))/2,j,k)
               rind(to) = rind(to) + ((d-1)*(d-2))/2
            end do
         end do
      endif

! if an old solution has been saved, for each export element find the element
! that has the old solution for that area (it can be the element or an ancestor)

      if (grid%oldsoln_exists) then

         gid = elemgid
         do
            elem = hash_decode_key(gid,grid%elem_hash)
            if (grid%element(elem)%oldleaf) exit
            if (grid%element(elem)%level == 1) then
               ierr = PHAML_INTERNAL_ERROR
               call fatal("failed to find element with old solution")
               stop
            endif
            gid = gid/2
         end do

! pack the gid of the element that has the old solution

         call hash_pack_key(gid,isend,iind(to))
         iind(to) = iind(to) + KEY_SIZE

! pack the element old solution and its dimensions

         if (associated(grid%element(elem)%oldsoln)) then
            isend(iind(to)  ) = size(grid%element(elem)%oldsoln,dim=1)
            isend(iind(to)+1) = size(grid%element(elem)%oldsoln,dim=2)
            isend(iind(to)+2) = size(grid%element(elem)%oldsoln,dim=3)
            iind(to) = iind(to) + 3
            do j=1,size(grid%element(elem)%oldsoln,dim=1)
               do k=1,size(grid%element(elem)%oldsoln,dim=2)
                  do l=1,size(grid%element(elem)%oldsoln,dim=3)
                     rsend(rind(to)) = grid%element(elem)%oldsoln(j,k,l)
                     rind(to) = rind(to) + 1
                   end do
               end do
            end do
         else
            isend(iind(to)  ) = 0
            isend(iind(to)+1) = 0
            isend(iind(to)+2) = 0
            iind(to) = iind(to) + 3
         endif

! pack the edge old solution and its dimensions

         do j=1,EDGES_PER_ELEMENT
            edge = grid%element(elem)%edge(j)
            if (associated(grid%edge(edge)%oldsoln)) then
               isend(iind(to)  ) = size(grid%edge(edge)%oldsoln,dim=1)
               isend(iind(to)+1) = size(grid%edge(edge)%oldsoln,dim=2)
               isend(iind(to)+2) = size(grid%edge(edge)%oldsoln,dim=3)
               iind(to) = iind(to) + 3
               do k=1,size(grid%edge(edge)%oldsoln,dim=1)
                  do l=1,size(grid%edge(edge)%oldsoln,dim=2)
                     do m=1,size(grid%edge(edge)%oldsoln,dim=3)
                        rsend(rind(to)) = grid%edge(edge)%oldsoln(k,l,m)
                        rind(to) = rind(to) + 1
                     end do
                  end do
               end do
            else
               isend(iind(to)  ) = 0
               isend(iind(to)+1) = 0
               isend(iind(to)+2) = 0
               iind(to) = iind(to) + 3
            endif
         end do

! pack the vertex old solution and its dimensions

         do j=1,VERTICES_PER_ELEMENT
            vert = grid%element(elem)%vertex(j)
            if (associated(grid%vertex_oldsoln)) then
               isend(iind(to)  ) = size(grid%vertex_oldsoln,dim=2)
               isend(iind(to)+1) = size(grid%vertex_oldsoln,dim=3)
               iind(to) = iind(to) + 2
               do k=1,size(grid%vertex_oldsoln,dim=2)
                  do l=1,size(grid%vertex_oldsoln,dim=3)
                     rsend(rind(to)) = grid%vertex_oldsoln(vert,k,l)
                     rind(to) = rind(to) + 1
                  end do
               end do
            else
               isend(iind(to)  ) = 0
               isend(iind(to)+1) = 0
               iind(to) = iind(to) + 2
            endif
         end do

      endif ! old solution exists

   end do ! next export gid

! exchange information with the other processors

   call start_watch((/cpdistribute,ctdistribute/))
   call phaml_alltoall(procs,isend,nisend,irecv,nirecv,530)
   call phaml_alltoall(procs,rsend,nrsend,rrecv,nrrecv,531)
   call stop_watch((/cpdistribute,ctdistribute/))

   deallocate(isend,rsend,stat=allocstat)

! derefine elements no longer needed

   call delete_unowned_elements(grid,refine_control)


! for each element received ...

   iindex = 1
   rindex = 1

   do part=1,nproc
      do while (iindex < sum(nirecv(1:part)))
         elemgid = hash_unpack_key(irecv,iindex)
         iindex = iindex + KEY_SIZE

! create the element, and any ancestors required by it

         call create_element(grid,elemgid,refine_control,errcode)
         if (errcode /= 0) then
            call fatal("grid full during redistribution",procs=procs)
            stop
         endif

! set myself as the owner

         elemlid = hash_decode_key(elemgid,grid%elem_hash)
         if (elemlid == HASH_NOT_FOUND) then
            ierr = PHAML_INTERNAL_ERROR
            call fatal("element just created doesn't have an elemlid in redistribute", &
                       procs=procs)
            stop
         endif
         if (.not. grid%element(elemlid)%iown) then
            grid%element(elemlid)%iown = .true.
            call own_parent(grid,elemgid)
            grid%nelem_leaf_own = grid%nelem_leaf_own + 1
            if (grid%element(elemlid)%degree >= 3) then
               grid%dof_own = grid%dof_own + grid%system_size * &
                  ((grid%element(elemlid)%degree-2)*(grid%element(elemlid)%degree-1))/2
            endif
            do j=1,VERTICES_PER_ELEMENT
               if (grid%vertex(grid%element(elemlid)%vertex(j))%assoc_elem == elemlid) then
                  grid%nvert_own = grid%nvert_own + 1
                  grid%dof_own = grid%dof_own + grid%system_size
               endif
            end do
            do j=1,EDGES_PER_ELEMENT
               if (grid%edge(grid%element(elemlid)%edge(j))%assoc_elem == elemlid) then
                  grid%nedge_own = grid%nedge_own - 1
                  if (grid%edge(grid%element(elemlid)%edge(j))%degree >= 2) then
                     grid%dof_own = grid%dof_own + grid%system_size * &
                           (grid%edge(grid%element(elemlid)%edge(j))%degree - 1)
                  endif
               endif
            end do
         endif

! set the edge and element degrees

         if (grid%element(elemlid)%isleaf) then
            do while (grid%element(elemlid)%degree < irecv(iindex+3))
               call p_refine_elem(grid,elemlid,refine_control)
            end do
         endif

! set the solution values of this element

! NOTE: for edges, it is possible that the current degree does not agree with
! the degree that was sent, because the p refinement of a neighboring element
! to bring up the degree may not have occured yet.  So use the minimum of
! edge degree and what was sent as the number to copy (the others will be
! copied when the neighbor is done), and increment the message pointer by
! the number that was sent.  This does not occur with vertices and elements.

! vertex solutions

         do j=1,VERTICES_PER_ELEMENT
            grid%vertex_solution(grid%element(elemlid)%vertex(j),:,:) = &
                   reshape(rrecv(rindex:rindex+grid%nsoln-1),&
                   (/grid%system_size,max(1,grid%num_eval)/))
            rindex = rindex + grid%nsoln
         end do

! edge solutions

         do j=1,EDGES_PER_ELEMENT
            d = min(grid%edge(grid%element(elemlid)%edge(j))%degree, &
                    irecv(iindex+j-1))
            do k=1,grid%system_size
               do l=1,max(1,grid%num_eval)
                  if (d > 1) then
                     grid%edge(grid%element(elemlid)%edge(j))%solution(1:d-1,k,l) = &
                            rrecv(rindex:rindex+d-2)
                  endif
                  rindex = rindex + irecv(iindex+j-1)-1
               end do
            end do
         end do

! element solutions

         d = min(grid%element(elemlid)%degree,irecv(iindex+3))
         if (d > 2) then
            do j=1,grid%system_size
               do k=1,max(1,grid%num_eval)
                  grid%element(elemlid)%solution(1:((d-1)*(d-2))/2,j,k) = &
                          rrecv(rindex:rindex+((d-1)*(d-2))/2-1)
                  rindex = rindex + ((irecv(iindex+3)-1)*(irecv(iindex+3)-2))/2
               end do
            end do
         endif

         iindex = iindex + 4

! if an old solution exists, set the old solution values sent

         if (grid%oldsoln_exists) then

            gid = hash_unpack_key(irecv,iindex)
            iindex = iindex + KEY_SIZE
            elem = hash_decode_key(gid,grid%elem_hash)
            if (elem == HASH_NOT_FOUND) then
               ierr = PHAML_INTERNAL_ERROR
               call fatal("element for old solution not found in redistribute",&
                          procs=procs)
               stop
            endif

! element old solution

            if (irecv(iindex) /= 0) then
               d1 = irecv(iindex  )
               d2 = irecv(iindex+1)
               d3 = irecv(iindex+2)
               if (associated(grid%element(elem)%oldsoln)) then
                  if (size(grid%element(elem)%oldsoln,dim=1) /= d1 .or. &
                      size(grid%element(elem)%oldsoln,dim=2) /= d2 .or. &
                      size(grid%element(elem)%oldsoln,dim=3) /= d3) then
                     deallocate(grid%element(elem)%oldsoln,stat=allocstat)
                  endif
               endif
               if (.not. associated(grid%element(elem)%oldsoln)) then
                  allocate(grid%element(elem)%oldsoln(d1,d2,d3),stat=allocstat)
                  if (allocstat /= 0) then
                     ierr = ALLOC_FAILED
                     call fatal("allocation failed in redistribute",procs=procs)
                     return
                  endif
               endif
               do j=1,d1
                  do k=1,d2
                     do l=1,d3
                        grid%element(elem)%oldsoln(j,k,l) = rrecv(rindex)
                        rindex = rindex + 1
                     end do
                  end do
               end do
            endif
            iindex = iindex + 3

! edge old solution

            do j=1,EDGES_PER_ELEMENT
               edge = grid%element(elem)%edge(j)
               if (irecv(iindex) /= 0) then
                  d1 = irecv(iindex  )
                  d2 = irecv(iindex+1)
                  d3 = irecv(iindex+2)
                  if (associated(grid%edge(edge)%oldsoln)) then
                     if (size(grid%edge(edge)%oldsoln,dim=1) /= d1 .or. &
                         size(grid%edge(edge)%oldsoln,dim=2) /= d2 .or. &
                         size(grid%edge(edge)%oldsoln,dim=3) /= d3) then
                        deallocate(grid%edge(edge)%oldsoln,stat=allocstat)
                     endif
                  endif
                  if (.not. associated(grid%edge(edge)%oldsoln)) then
                     allocate(grid%edge(edge)%oldsoln(d1,d2,d3),stat=allocstat)
                     if (allocstat /= 0) then
                        ierr = ALLOC_FAILED
                        call fatal("allocation failed in redistribute",procs=procs)
                        return
                     endif
                  endif
                  do k=1,d1
                     do l=1,d2
                        do m=1,d3
                           grid%edge(edge)%oldsoln(k,l,m) = rrecv(rindex)
                           rindex = rindex + 1
                        end do
                     end do
                  end do
               endif
               iindex = iindex + 3
            end do

! vertex old solution

            do j=1,VERTICES_PER_ELEMENT
               vert = grid%element(elem)%vertex(j)
               if (irecv(iindex) /= 0) then
                  d1 = irecv(iindex  )
                  d2 = irecv(iindex+1)
                  do k=1,d1
                     do l=1,d2
                        grid%vertex_oldsoln(vert,k,l) = rrecv(rindex)
                        rindex = rindex + 1
                     end do
                  end do
               endif
               iindex = iindex + 2
            end do

! set this element as an oldleaf

            grid%element(elem)%oldleaf = .true.

         endif ! old solution exists
      end do ! next imported element
   end do ! next processor

! clean up oldleaf so each point in the domain has only one oldleaf

   if (grid%oldsoln_exists) then
      call cleanup_oldleaf(grid)
   endif

! free memory

   if (associated(irecv)) deallocate(irecv,stat=allocstat)
   if (associated(rrecv)) deallocate(rrecv,stat=allocstat)
   deallocate(nirecv,nrrecv,nisend,nrsend,iind,rind,stat=allocstat)

   grid%errind_up2date = .false.

else ! skip_comm

! derefine elements no longer needed

   call delete_unowned_elements(grid,refine_control)

endif ! .not. skip_comm

! stop timing the redistribution process

call stop_watch((/pdistribute,tdistribute/))

end subroutine redistribute

!          -----------------------
subroutine delete_unowned_elements(grid,refcont)
!          -----------------------

!----------------------------------------------------
! This routine removes any elements that are not owned by this partition
! and are not needed for compatibility.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(refine_options), intent(in) :: refcont
!----------------------------------------------------
! Local variables:

integer :: elem
logical :: owned_child
!----------------------------------------------------
! Begin executable code

elem = grid%head_level_elem(1)
do while (elem /= END_OF_LIST)
   call delete_unowned_elements_recur(grid,elem,owned_child,refcont)
   elem = grid%element(elem)%next
end do

end subroutine delete_unowned_elements

!                    -----------------------------
recursive subroutine delete_unowned_elements_recur(grid,elem,owned_child, &
                                                   refcont)
!                    -----------------------------

!----------------------------------------------------
! This routine does the work of delete_unowned_elements
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: elem
logical, intent(out) :: owned_child
type(refine_options), intent(in) :: refcont

!----------------------------------------------------
! Local variables:

integer :: child(MAX_CHILD), i, mate, errcode, allc(MAX_CHILD)
logical :: this_owned_child
!----------------------------------------------------
! Begin executable code

allc = ALL_CHILDREN
child = get_child_lid(grid%element(elem)%gid,allc,grid%elem_hash)

! if this is a leaf, return with an indication of whether it is owned by
! this partition

if (child(1) == NO_CHILD) then
   owned_child = grid%element(elem)%iown
else

! otherwise, see if any of the descendents are owned

   owned_child = .false.
   do i=1,MAX_CHILD
      if (child(i) == NO_CHILD) cycle
      call delete_unowned_elements_recur(grid,child(i),this_owned_child,refcont)
      owned_child = owned_child .or. this_owned_child
   end do

! also check the descendents of the mate

   if (.not. owned_child) then
      if (.not. (grid%element(elem)%mate == BOUNDARY)) then
         mate = hash_decode_key(grid%element(elem)%mate,grid%elem_hash)
         if (mate /= HASH_NOT_FOUND) then
            child = get_child_lid(grid%element(mate)%gid,allc, &
                                        grid%elem_hash)
            do i=1,MAX_CHILD
              if (child(i) == NO_CHILD) cycle
              call delete_unowned_elements_recur(grid,child(i), &
                                                 this_owned_child,refcont)
              owned_child = owned_child .or. this_owned_child
            end do
         endif
      endif
   endif

! if not, derefine this element

   if (.not. owned_child) then
      call unbisect_triangle_pair(grid,elem,errcode,refcont)
   endif

endif

end subroutine delete_unowned_elements_recur

!                    -------------
recursive subroutine disown_parent(grid,child)
!                    -------------

!----------------------------------------------------
! This routine removes ownership of all ancestors of child
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(hash_key), intent(in) :: child
!----------------------------------------------------
! Local variables:

type(hash_key) :: parent
integer :: lid
!----------------------------------------------------
! Begin executable code

parent = child/MAX_CHILD
lid = hash_decode_key(parent,grid%elem_hash)
if (lid == HASH_NOT_FOUND) return

if (grid%element(lid)%iown) then
   grid%element(lid)%iown = .false.
   call disown_parent(grid,parent)
endif

end subroutine disown_parent

!                    ----------
recursive subroutine own_parent(grid,child)
!                    ----------

!----------------------------------------------------
! This routine adds ownership of ancestors of child for which
! all descendents are owned
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(hash_key), intent(in) :: child
!----------------------------------------------------
! Local variables:

type(hash_key) :: parent
integer :: children(MAX_CHILD), i, lid, allc(MAX_CHILD)
!----------------------------------------------------
! Begin executable code

parent = child/MAX_CHILD
lid = hash_decode_key(parent,grid%elem_hash)
if (lid == HASH_NOT_FOUND) return
allc = ALL_CHILDREN
children = get_child_lid(parent,allc,grid%elem_hash)

grid%element(lid)%iown = .true.
do i=1,MAX_CHILD
   if (children(i) /= NO_CHILD) then
      grid%element(lid)%iown = grid%element(lid)%iown .and. &
                               grid%element(children(i))%iown
   endif
end do

if (grid%element(lid)%iown) call own_parent(grid,parent)

end subroutine own_parent

!          ---------------
subroutine cleanup_oldleaf(grid)
!          ---------------

!----------------------------------------------------
! This routine cleans up oldleaf to make sure any point has only one
! element marked as a leaf for the old solution.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
!----------------------------------------------------
! Local variables:

integer :: elem
logical :: dum1, dum2
!----------------------------------------------------
! Begin executable code

! pass through all elements of the initial grid calling the worker recursive
! routine

elem = grid%head_level_elem(1)
do while (elem /= END_OF_LIST)
   call cleanup_oldleaf_recur(grid,elem,dum1,dum2)
   elem = grid%element(elem)%next
end do

end subroutine cleanup_oldleaf

!                    ---------------------
recursive subroutine cleanup_oldleaf_recur(grid,parent,any_oldleaf,all_oldleaf)
!                    ---------------------

!----------------------------------------------------
! This routine does the work of cleanup_oldleaf
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: parent
logical, intent(out) :: any_oldleaf, all_oldleaf
!----------------------------------------------------
! Local variables:

integer :: children(MAX_CHILD), i, allc(MAX_CHILD)
logical :: any_child, all_child
!----------------------------------------------------
! Begin executable code

! any_oldleaf for parent is true if any of its area is an oldleaf.
! all_oldleaf for parent is true if all of its area is covered by oldleafs

allc = ALL_CHILDREN
children = get_child_lid(grid%element(parent)%gid,allc,grid%elem_hash)

! if this is a leaf, then any_oldleaf and all_oldleaf are the same as oldleaf

if (children(1) == NO_CHILD) then
   any_oldleaf = grid%element(parent)%oldleaf
   all_oldleaf = grid%element(parent)%oldleaf

! otherwise, any_oldleaf is the or of the children and all_oldleaf is the and

else

   any_oldleaf = .false.
   all_oldleaf = .true.
   do i=1,MAX_CHILD
      call cleanup_oldleaf_recur(grid,children(i),any_child,all_child)
      any_oldleaf = any_oldleaf .or. any_child
      all_oldleaf = all_oldleaf .and. all_child
   end do

! if this non-leaf element is marked as an oldleaf, then if it is all_oldleaf
! unmark it as oldleaf, and if it is any_oldleaf but not all_oldleaf, then
! unmark any descendants that are oldleaf

   if (grid%element(parent)%oldleaf) then
      if (all_oldleaf) then
         grid%element(parent)%oldleaf = .false.
      elseif (any_oldleaf) then
         do i=1,MAX_CHILD
            call unset_oldleaf(grid,children(i))
         end do
      endif
      any_oldleaf = .true.
      all_oldleaf = .true.
   endif

endif

end subroutine cleanup_oldleaf_recur

!                    -------------
recursive subroutine unset_oldleaf(grid,parent)
!                    -------------

!----------------------------------------------------
! This routine unsets oldleaf in all descendents of parent and parent
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: parent
!----------------------------------------------------
! Local variables:

integer :: children(MAX_CHILD), i, allc(MAX_CHILD)
!----------------------------------------------------
! Begin executable code

grid%element(parent)%oldleaf = .false.

allc = ALL_CHILDREN
children = get_child_lid(grid%element(parent)%gid,allc,grid%elem_hash)

if (children(1) /= NO_CHILD) then
   do i=1,MAX_CHILD
      call unset_oldleaf(grid,children(i))
   end do
endif

end subroutine unset_oldleaf

!          -------------
subroutine check_balance(grid,procs,what,fileunit)
!          -------------

!----------------------------------------------------
! This routine prints information about how well "what" is balanced.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
type (proc_info), intent(in) :: procs
integer, intent(in) :: what
integer, intent(in), optional :: fileunit
!----------------------------------------------------
! Local variables:

integer :: nproc, i, p, ni, nr
integer, allocatable :: nentities(:), neq(:)
integer, pointer :: irecv(:)
real(my_real), pointer :: rrecv(:)
real(my_real) :: ave, aveeq
!----------------------------------------------------
! Begin executable code

if (what == BALANCE_NONE) return

nproc = num_proc(procs)

if (my_proc(procs) == MASTER) then

! receive the number of entities owned by each processor

   allocate(nentities(nproc),neq(nproc))
   do i=1,num_proc(procs)
      call phaml_recv(procs,p,irecv,ni,rrecv,nr,501)
      nentities(p) = irecv(1)
      neq(p) = irecv(2)
      deallocate(irecv)
      if (associated(rrecv)) deallocate(rrecv)
   end do

! print the number of entities owned by each

   write(outunit,"(A)")
   write(outunit,"(A)") "checking processor balance"
   write(outunit,"(A)")
   select case(what)
   case (BALANCE_ELEMENTS)
      write(outunit,"(A)") "number of elements owned by each processor"
   case (BALANCE_VERTICES)
      write(outunit,"(A)") "number of vertices owned by each processor"
   case (BALANCE_EQUATIONS)
      write(outunit,"(A)") "number of equations owned by each processor"
   case default
      call warning("bad value for what in check_balance")
   end select
! RESTRICTION no more than 1024 processors
   write(outunit,"(1024I11)") nentities

! print the biggest deviation from the mean

   ave = sum(nentities)/real(nproc,my_real)
   aveeq = sum(neq)/real(nproc,my_real)
   write(outunit,"(A,E18.10E2)") "Largest deviation from the mean (per cent)",&
                                 100*maxval(abs(nentities-ave))/ave
   if (present(fileunit)) write(fileunit,"(i10,2e18.10e2)") sum(neq), &
                                 100*maxval(abs(nentities-ave))/ave, &
                                 100*maxval(abs(neq-aveeq))/aveeq

   deallocate(nentities,neq)

! barrier so the slaves don't call this again before the master finishes

   call phaml_send(procs,1,(/0/),1,(/0.0_my_real/),0,502)

else ! SLAVE

! get the number of entities I own

   allocate(nentities(2))
   select case(what)
   case (BALANCE_ELEMENTS)
      nentities(1) = grid%nelem_leaf_own
      nentities(2) = count_dof(grid,just_owned=.true.)
   case (BALANCE_VERTICES)
      nentities(1) = grid%nvert_own
      nentities(2) = count_dof(grid,just_owned=.true.)
   case (BALANCE_EQUATIONS)
      nentities(1) = count_dof(grid,just_owned=.true.)
      nentities(2) = count_dof(grid,just_owned=.true.)
   case default
      call warning("bad value for what in check_balance")
      nentities = 0
   end select

! send it to the master

   call phaml_send(procs,MASTER,nentities,2,(/0.0_my_real/),0,501)

   deallocate(nentities)

! barrier so the slaves don't call this again before the master finishes

   if (my_proc(procs) == 1) then
      call phaml_recv(procs,p,irecv,ni,rrecv,nr,502)
      deallocate(irecv)
      if (associated(rrecv)) deallocate(rrecv)
   endif
   call phaml_barrier(procs)

endif ! MASTER or SLAVE

end subroutine check_balance

end module load_balance
