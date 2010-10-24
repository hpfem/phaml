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

module refine_adapt_mod

!----------------------------------------------------
! This module contains routines for the basic control of adaptive refinement.
!
! communication tags in this module are of the form 33xx
!----------------------------------------------------

!----------------------------------------------------
! Other modules used:

use global
use gridtype_mod
use grid_util
use linsystype_mod
use hp_strategies
use refine_elements
use message_passing
use hash_mod
use linsys_io
use load_balance
use error_estimators
!----------------------------------------------------

implicit none
private
public refine_adaptive

!----------------------------------------------------
! Non-module procedures used are:

!----------------------------------------------------
! The following defined constants are defined:

!----------------------------------------------------
! The following variables are defined:

!----------------------------------------------------

contains

!          ---------------
subroutine refine_adaptive(grid,procs,refine_control,solver_control,io_control,&
                           still_sequential,init_nvert,init_nelem,init_dof, &
                           loop,balance_what,predictive,no_time)
!          ---------------

!----------------------------------------------------
! This routine refines the grid in accordance with the parameters in
! refine_control.
! If derefine is true, the grid is first derefined to remove elements
! with very small error indicators.
! If refterm is DOUBLE_NVERT (NELEM, NEQ) then the number of vertices
! (elements, degrees of freedom) in the global grid is increased to
! init_nvert*inc_factor**loop. (init_nelem, init_dof)
! Each processor determines how many vertices (elements, dofs) by which to
! increase the grid based on the error indicators relative to those of other
! processors, unless load balancing is predictive in which case all processors
! get an equal share of the vertices (elements, dofs).
! If refterm is HALVE_ERREST then the maximum error indicator is reduced
! by the factor inc_factor.  However, if this leads to a very small change
! in the grid, additional refinement is performed.
! If refterm is KEEP_NVERT (NELEM, ERREST) then the number of vertices
! (elements) is kept constant, or the maximum error indicator is kept
! approximately constant.
! If refterm is ONE_REF then elements whose error indicator is larger than
! reftol are refined, and each element is refined at most once.
! If refterm is ONE_REF_HALF_ERRIND then elements whose error indicator is
! larger than maximum_error_indicator / inc_factor are refined once.
! If refterm ends with SMOOTH then after the criteria is met, refinement
! continues until the current error indicator bin is emptied.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type (proc_info), target, intent(in) :: procs
type(refine_options), intent(in) :: refine_control
type(solver_options), intent(in) :: solver_control
type(io_options), intent(in) :: io_control
integer, intent(in) :: init_nvert, init_nelem, init_dof, loop, balance_what
logical, intent(in) :: still_sequential, predictive
logical, intent(in), optional :: no_time

!----------------------------------------------------
! Local variables:

integer :: elem, errcode, astat, nproc, my_processor, target, mate, child(2)
real(my_real) :: global_max_errind
integer, pointer :: numhref(:), numpref(:)
! If len=1 is changed, also change it for temp_reftype in more_elements
character(len=1), pointer :: reftype(:)
type(errind_list) :: elist
logical :: return_to_elist, complete_elist, one_elist, target_met

!----------------------------------------------------
! Begin executable code

! MASTER doesn't participate

if (my_proc(procs) == MASTER) return

nproc = num_proc(procs)
my_processor = my_proc(procs)

! compute error indicators if needed

if (.not. grid%errind_up2date) then
   call all_error_indicators(grid,refine_control%error_estimator)
endif

! compute maximum error indicator

global_max_errind = compute_global_max_errind(grid,procs)

! set the target for terminating refinement

call set_target(target,grid,procs,refine_control,global_max_errind, &
                init_nvert,init_nelem,init_dof,loop,balance_what, &
                still_sequential,predictive)

! derefine the grid, if requested

if (refine_control%derefine) then
   call derefinement(grid,procs,refine_control,global_max_errind,target)
endif

! set the number of times to refine each element

allocate(numhref(size(grid%element)),numpref(size(grid%element)),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in refine",procs=procs)
   print *,"astat ",astat ! TEMP090814 for g95 problem
   return
endif

call set_numref(grid,procs,still_sequential,refine_control,numhref,numpref)

! create the lists that group elements into bins based on the error indicator

allocate(elist%next_errind(size(grid%element)), &
         elist%prev_errind(size(grid%element)),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in refine",procs=procs)
   return
endif

call make_elist(grid,procs,elist,refine_control,global_max_errind, &
                numhref,numpref,predictive,still_sequential)
elist%current_list = 1

! mark elements to be refined by h, by p or not, or leave undecided

allocate(reftype(size(grid%element)),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in refine",procs=procs)
   return
endif

call mark_reftype(grid,refine_control,global_max_errind,reftype,elist)

! set controls for the refine loop:
! return_to_elist, tells if children get put into an elist
! complete_elist, tells if you must finish the whole elist (SMOOTH)
! one_elist, tells if you quit after finishing the first elist

call set_loop_control(refine_control,return_to_elist,complete_elist,one_elist)

! main refinement loop

errcode = 0
do

! quit if the grid is full

   if (errcode /= 0) exit

! see if the target has been met

   target_met = check_target(grid,refine_control,elist,target)

! if the target is met and we don't complete lists, done

   if ((.not. complete_elist) .and. target_met) exit

! get the next element to refine

   elem = elist%head_errind(elist%current_list)

! if this is the end of the current list and either we do just the first list or
! the target has been met, done

   if (elem == END_OF_LIST .and. (one_elist .or. target_met)) exit

! if this is the end of the current list, find the beginning of the next
! nonempty list

   do while (elem==END_OF_LIST .and. elist%current_list<size(elist%head_errind))
      elist%current_list = elist%current_list + 1
      elem = elist%head_errind(elist%current_list)
   end do

! if all lists are empty, done

   if (elem == END_OF_LIST) exit

! if we return elements to the elists, complete the current elist, and have
! reached the last list, exit to avoid an infinite loop

   if (return_to_elist .and. complete_elist .and. &
       elist%current_list == size(elist%head_errind)) exit

! if we haven't determined the refinement type, do it now

   if (reftype(elem) == "u") then
      call mark_reftype_one(elem,grid,refine_control,global_max_errind, &
                            reftype,.false.)
   endif

! refine the element and determine the type of refinement of the new element(s)

   if (reftype(elem) == "h") then
      call bisect_triangle_pair(grid,elem,errcode,refine_control,elist,reftype,&
                                numhref,numpref,return_to_elist)
      reftype(elem) = "n"
      child = get_child_lid(grid%element(elem)%gid,(/1,2/),grid%elem_hash)
      call mark_reftype_one(child(1),grid,refine_control,1.0_my_real,reftype, &
                            .true.)
      call mark_reftype_one(child(2),grid,refine_control,1.0_my_real,reftype, &
                            .true.)
      if (grid%element(elem)%mate == BOUNDARY) then
         mate = BOUNDARY
      else
         mate = hash_decode_key(grid%element(elem)%mate,grid%elem_hash)
      endif
      if (mate /= BOUNDARY) then
         reftype(mate) = "n"
         child = get_child_lid(grid%element(mate)%gid,(/1,2/),grid%elem_hash)
         call mark_reftype_one(child(1),grid,refine_control,1.0_my_real, &
                               reftype,.true.)
         call mark_reftype_one(child(2),grid,refine_control,1.0_my_real, &
                               reftype,.true.)
      endif
   elseif (reftype(elem) == "p") then
      call p_refine_elem(grid,elem,refine_control,elist,numpref,return_to_elist)
      call mark_reftype_one(elem,grid,refine_control,1.0_my_real,reftype, &
                            .true.)
   else
      call remove_from_errind_list(elem,elist)
   endif

end do ! main refine loop

! free memory

deallocate(elist%next_errind,elist%prev_errind,reftype,numhref,numpref, &
           stat=astat)

end subroutine refine_adaptive

!          ----------
subroutine set_target(target,grid,procs,refine_control,global_max_errind, &
                      init_nvert,init_nelem,init_dof,loop,balance_what, &
                      still_sequential,predictive)
!          ----------

!----------------------------------------------------
! This routine determines the value to use as the target for terminating
! refinement
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(out) :: target
type(grid_type), intent(inout) :: grid
type (proc_info), target, intent(in) :: procs
type(refine_options), intent(in) :: refine_control
real(my_real), intent(in) :: global_max_errind
integer, intent(in) :: init_nvert, init_nelem, init_dof, loop, balance_what
logical, intent(in) :: still_sequential, predictive
!----------------------------------------------------
! Local variables:

real(my_real) :: my_fraction
integer :: my_num_big, lev, elem, total_num_big, total
!----------------------------------------------------
! Begin executable code

! cases where we don't need target

if (refine_control%refterm == ONE_REF .or. &
    refine_control%refterm == ONE_REF_HALF_ERRIND) then
   target = 0
   return
endif

if (refine_control%reftype == HP_ADAPTIVE .and. &
    refine_control%hp_strategy == HP_T3S) then
   if (refine_control%t3s_reftype == H_UNIFORM) then
      target = 0
      return
   endif
endif

! cases that don't need my_fraction

if (refine_control%refterm == HALVE_ERREST) then
   target = max(1,nint(log(refine_control%inc_factor)/log(binw)))
   return
endif

if (refine_control%refterm == KEEP_ERREST) then
! TEMP not doing KEEP_ERREST
   call warning("have not yet decided how to handle refterm==KEEP_ERREST")
   target = 0
endif

! set the target for the termination criterion

! if load balancing is not predictive, determine my fraction by comparing the
! number of elements with large error indicators on this processor with those
! on other processors.  For now, large means it will be in the first two error
! indicator lists

if (still_sequential) then
   my_fraction = 1.0_my_real
elseif (predictive) then
   my_fraction = my_weight_fraction(grid,procs,predictive, &
                                    balance_what,refine_control)
else
   my_num_big = 0
   do lev=1,grid%nlev
      elem = grid%head_level_elem(lev)
      do while (elem /= END_OF_LIST)
         if (grid%element(elem)%isleaf) then
            if (maxval(grid%element_errind(elem,:))/grid%element(elem)%work > &
                global_max_errind/(binw**2)) then
               my_num_big = my_num_big + 1
            endif
         endif
         elem = grid%element(elem)%next
      end do
   end do
   total_num_big = phaml_global_sum(procs,my_num_big,3332)
   my_fraction = my_num_big/real(total_num_big,my_real)
endif

! set the target value depending on what type of termination is used

select case (refine_control%refterm)

case (DOUBLE_NVERT, DOUBLE_NVERT_SMOOTH)

   total = init_nvert*refine_control%inc_factor**loop
   if (total > refine_control%max_vert) total = refine_control%max_vert
   target = ceiling(total*my_fraction)

case (DOUBLE_NELEM, DOUBLE_NELEM_SMOOTH)

   total = init_nelem*refine_control%inc_factor**loop
   if (total > refine_control%max_elem) total = refine_control%max_elem
   target = ceiling(total*my_fraction)

case (DOUBLE_NEQ, DOUBLE_NEQ_SMOOTH)

   total = init_dof*refine_control%inc_factor**loop
   if (total > refine_control%max_dof) total = refine_control%max_dof
   target = ceiling(total*my_fraction)

case (KEEP_NVERT, KEEP_NVERT_SMOOTH)

   call get_grid_info(grid,procs,still_sequential,3345,total_nvert=total,&
                      no_master=.true.)
   if (total > refine_control%max_vert) total = refine_control%max_vert
   target = total*my_fraction

case (KEEP_NELEM, KEEP_NELEM_SMOOTH)

   call get_grid_info(grid,procs,still_sequential,3345, &
                      total_nelem_leaf=total,no_master=.true.)
   if (total > refine_control%max_elem) total = refine_control%max_elem
   target = total*my_fraction

case (KEEP_NEQ, KEEP_NEQ_SMOOTH)

   call get_grid_info(grid,procs,still_sequential,3345, &
                      total_dof=total,no_master=.true.)
   if (total > refine_control%max_dof) total = refine_control%max_dof
   target = total*my_fraction

case default

   call fatal("illegal value for refterm",procs=procs)
   stop

end select

end subroutine set_target

!        ------------------
function my_weight_fraction(grid,procs,predictive,balance_what,refine_control)
!        ------------------

!----------------------------------------------------
! This routine computes the fraction of the total weight of leaf elements that
! belongs to this processor.  For homogeneous load balancing, it should be
! 1/nproc.  This is really only useful if load balancing is done for a
! heterogeneous computing system.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(proc_info), intent(in) :: procs
logical, intent(in) :: predictive
integer, intent(in) :: balance_what
type(refine_options), intent(in) :: refine_control
real(my_real) :: my_weight_fraction
!----------------------------------------------------
! Local variables:

real(my_real) :: my_total_weight, total_weight
integer :: lev, elem
!----------------------------------------------------
! Begin executable code

! set the current weights

call set_weights(grid,predictive,balance_what,refine_control,procs)

my_total_weight = 0.0_my_real
do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (grid%element(elem)%isleaf .and. grid%element(elem)%iown) then
         my_total_weight = my_total_weight + grid%element(elem)%weight
      endif
      elem = grid%element(elem)%next
   end do
end do

total_weight = phaml_global_sum(procs,my_total_weight,3321)

if (total_weight /= 0.0_my_real) then
   my_weight_fraction = my_total_weight/total_weight
else
   my_weight_fraction = 1.0_my_real
endif

end function my_weight_fraction

!          ------------
subroutine derefinement(grid,procs,refine_control,max_errind,target)
!          ------------

!----------------------------------------------------
! This routine performs h and p derefinement of the grid
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type (proc_info), intent(in) :: procs
type(refine_options), intent(in) :: refine_control
real(my_real), intent(in) :: max_errind
integer, intent(in) :: target
!----------------------------------------------------
! Local variables:

integer :: deref_list(2*size(grid%element)), degree_list(2*size(grid%element))
integer :: num_deref, lev, elem, parent, siblings(4), i, errcode, astat, &
           nproc, my_processor, p, next_elem, isub
! newcomm
integer :: nsend
integer, allocatable :: isend(:), nrecv(:)
integer, pointer :: irecv(:)
logical(small_logical) :: hcoarsen_checked(size(grid%element))
logical :: hcoarsen_possible, pcoarsen_possible, hcoarsen_occured, &
           pcoarsen_occured
real(my_real) :: hcoarsen_errind, pcoarsen_errind, small_errind
!----------------------------------------------------
! Begin executable code

! no derefinement for uniform refinement

if (refine_control%reftype == HP_ADAPTIVE .and. &
    refine_control%hp_strategy == HP_T3S) then
   if (refine_control%t3s_reftype == H_UNIFORM) then
      return
   endif
endif

! unrefine all unrefineable elements whose h or p coarsening error indicator
! is small enough.  If that's not enough, increase the size allowed for the
! coarsening indicator and repeat.

num_deref = 0
small_errind = max_errind/100

do

   hcoarsen_checked = .false.

! go through the elements from finest level to coarsest level to allow
! parents of unrefined elements to also be unrefined

   do lev=grid%nlev,1,-1
      elem = grid%head_level_elem(lev)
      do while (elem /= END_OF_LIST)

! relatives

         if (grid%element(elem)%level /= 1) then
            parent = hash_decode_key(grid%element(elem)%gid/2,grid%elem_hash)
            siblings(1:2) = get_child_lid(grid%element(parent)%gid,(/1,2/), &
                                          grid%elem_hash)
            if (grid%element(parent)%mate == BOUNDARY) then
               siblings(3:4) = NO_CHILD
            else
               siblings(3:4) = get_child_lid( &
                              grid%element(parent)%mate,(/1,2/),grid%elem_hash)
            endif
         endif

! Check conditions that allow h coarsening.
! If the element has already been checked for h coarsening (via a sibling)
! then don't need to check it again.
! Cannot h coarsen level 1 elements, elements that have children, elements
! whose siblings have children, or elements whose parent or parent's mate have
! children owned by different processors.  Only h coarsen elements that I own
! or I own a sibling.  Only do h coarsening with h or hp adaptive refinement.

         hcoarsen_possible = .not. hcoarsen_checked(elem)
         if (refine_control%reftype /= H_ADAPTIVE .and. &
             refine_control%reftype /= HP_ADAPTIVE) hcoarsen_possible = .false.
         if (grid%element(elem)%level == 1) hcoarsen_possible = .false.
         if (.not. grid%element(elem)%isleaf) hcoarsen_possible=.false.
         if (hcoarsen_possible) then
            do i=1,4
               if (i==3 .and. grid%element(parent)%mate == BOUNDARY) exit
               if (.not. grid%element(siblings(i))%isleaf) then
                  hcoarsen_possible = .false.
               endif
            end do
            if (grid%element(siblings(1))%iown .neqv. &
                grid%element(siblings(2))%iown) then
               hcoarsen_possible = .false.
            endif
            if (.not. grid%element(parent)%mate == BOUNDARY) then
               if (grid%element(siblings(3))%iown .neqv. &
                   grid%element(siblings(4))%iown) then
                  hcoarsen_possible = .false.
               endif
            endif
            if (grid%element(parent)%mate == BOUNDARY) then
               if (.not. grid%element(elem)%iown) hcoarsen_possible = .false.
            else
               if (.not. grid%element(elem)%iown .and. &
                .not. grid%element(siblings(3))%iown) hcoarsen_possible=.false.
            endif
         endif

! Check conditions that allow p coarsening.
! Cannot p coarsen an element that has p==1.  Only leaves need be p coarsened.
! Only p coarsen elements that I own.  Only do p coarsening with p or hp
! adaptive.

         pcoarsen_possible = grid%element(elem)%isleaf .and. &
                             grid%element(elem)%iown .and. &
                             grid%element(elem)%degree > 1
         if (refine_control%reftype /= P_ADAPTIVE .and. &
             refine_control%reftype /= HP_ADAPTIVE) pcoarsen_possible = .false.

! Compute coarsening error indicators

         if (hcoarsen_possible) then
            hcoarsen_errind = hcoarsen_indicator(grid,parent,siblings)
         else
            hcoarsen_errind = huge(0.0_my_real)
         endif
         if (pcoarsen_possible) then
            pcoarsen_errind = pcoarsen_indicator(grid,elem)
         else
            pcoarsen_errind = huge(0.0_my_real)
         endif

         hcoarsen_occured = .false.
         pcoarsen_occured = .false.

! If either coarsening indicator is small enough, perform the type of
! coarsening that has the least impact on the solution

         if (hcoarsen_errind < pcoarsen_errind .and. &
             hcoarsen_errind < small_errind) then

! find the next element that is not a sibling before doing h coarsening

            next_elem = grid%element(elem)%next
            do while (next_elem == siblings(1) .or. next_elem == siblings(2) &
                 .or. next_elem == siblings(3) .or. next_elem == siblings(4))
               next_elem = grid%element(next_elem)%next
            end do

! perform h coarsening

            call unbisect_triangle_pair(grid,parent,errcode,refine_control)
            if (errcode == 0) then
               num_deref = num_deref + 1
               deref_list(num_deref) = parent
               degree_list(num_deref) = -1
               hcoarsen_occured = .true.
            endif

         elseif (pcoarsen_errind < small_errind) then

! perform p coarsening

            call p_coarsen_elem(grid,elem,errcode,refine_control)
            if (errcode == 0) then
               if (num_deref == 0) then
                  num_deref = num_deref + 1
                  deref_list(num_deref) = elem
                  degree_list(num_deref) = -2
               elseif (deref_list(num_deref) /= elem .or. &
                       degree_list(num_deref) /= -2) then
                  num_deref = num_deref + 1
                  deref_list(num_deref) = elem
                  degree_list(num_deref) = -2
               endif
               pcoarsen_occured = .true.
            endif
         endif


! if p coarsening occured, do the same element again.
! otherwise, next element

         if (.not. pcoarsen_occured) then
            if (grid%element(elem)%level /= 1) then
               hcoarsen_checked(siblings(1)) = .true.
               hcoarsen_checked(siblings(2)) = .true.
               if (.not. grid%element(parent)%mate == BOUNDARY) then
                  hcoarsen_checked(siblings(3)) = .true.
                  hcoarsen_checked(siblings(4)) = .true.
               endif
            endif
            if (hcoarsen_occured) then
               elem = next_elem
            else
               elem = grid%element(elem)%next
            endif
         endif
      end do ! next element
   end do ! next level

! see if enough derefinement has occured

   select case (refine_control%refterm)
   case (DOUBLE_NVERT, KEEP_NVERT, DOUBLE_NVERT_SMOOTH, KEEP_NVERT_SMOOTH)
      if (grid%nvert_own <= target) exit
   case (DOUBLE_NELEM, KEEP_NELEM, DOUBLE_NELEM_SMOOTH, KEEP_NELEM_SMOOTH)
      if (grid%nelem_leaf_own <= target) exit
   case (DOUBLE_NEQ, DOUBLE_NEQ_SMOOTH, KEEP_NEQ, KEEP_NEQ_SMOOTH)
      if (grid%dof_own <= target) exit
   case default
      exit
   end select
   small_errind = 2*small_errind
   if (small_errind > max_errind) exit

end do

! Send the lists of derefined elements to other processors and send the
! current degree of elements that were p-derefined and -1 to indicate
! h-derefinement

nproc = num_proc(procs)
my_processor = my_proc(procs)

nsend = (KEY_SIZE+1)*num_deref
allocate(isend(nsend),nrecv(nproc),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in derefinement",procs=procs)
   stop
endif

do i=1,num_deref
   call hash_pack_key(grid%element(deref_list(i))%gid,isend, &
                  1+(i-1)*(KEY_SIZE+1))
   if (degree_list(i) == -1) then
      isend(i*(KEY_SIZE+1)) = -1
   else
      isend(i*(KEY_SIZE+1)) = grid%element(deref_list(i))%degree
   endif
enddo

call phaml_alltoall(procs,isend,nsend,irecv,nrecv,3350)

! Derefine elements that I have but were derefined by the owner.
! Also note which ones cannot be derefined.

num_deref = 0
isub = 1
do p=1,nproc
   if (p == my_processor) then
      isub = isub + nrecv(p)
      cycle
   endif
   do i=1,nrecv(p)/(KEY_SIZE+1)
      elem = hash_decode_key(hash_unpack_key(irecv,isub),grid%elem_hash)
      if (elem /= HASH_NOT_FOUND) then
         if (irecv(isub+KEY_SIZE) == -1) then
            call unbisect_triangle_pair(grid,elem,errcode,refine_control)
            if (errcode /= 0 .and. errcode /= -2) then
               num_deref = num_deref + 1
               deref_list(num_deref) = elem
               degree_list(num_deref) = -p
            endif
         else
            call p_coarsen_elem(grid,elem,errcode,refine_control)
            if (errcode /= 0) then
               num_deref = num_deref + 1
               deref_list(num_deref) = elem
               degree_list(num_deref) = 1
            endif
         endif
      endif
      isub = isub + KEY_SIZE+1
   end do
end do
if (associated(irecv)) deallocate(irecv)

! send the list of elements which could not be derefined.  For p refinement
! send the degree; for h refinement send -owner so only the owner rerefines

deallocate(isend)
nsend = (KEY_SIZE+1)*num_deref
allocate(isend(nsend),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in derefinement",procs=procs)
   stop
endif

do i=1,num_deref
   call hash_pack_key(grid%element(deref_list(i))%gid,isend, &
                      (1+(i-1)*(KEY_SIZE+1)))
   if (degree_list(i) < 0) then
      isend(i*(KEY_SIZE+1)) = degree_list(i)
   else
      isend(i*(KEY_SIZE+1)) = grid%element(deref_list(i))%degree
   endif
end do

call phaml_alltoall(procs,isend,nsend,irecv,nrecv,3351)

! rerefine elements that were returned as underefineable by some processor

isub = 1
do p=1,nproc
   do i=1,nrecv(p)/(KEY_SIZE+1)
      if (irecv(isub+KEY_SIZE) == -my_processor) then
         call create_element(grid,2*hash_unpack_key(irecv,isub), &
                             refine_control,errcode)
      elseif (irecv(isub+KEY_SIZE) > 0) then
         elem = hash_decode_key(hash_unpack_key(irecv,isub),grid%elem_hash)
         if (elem /= HASH_NOT_FOUND) then
            if (grid%element(elem)%isleaf) then
               do while (grid%element(elem)%degree < irecv(isub+KEY_SIZE))
                  call p_refine_elem(grid,elem,refine_control)
               end do
            endif
         endif
      endif
      isub = isub + KEY_SIZE+1
   end do
end do
if (associated(irecv)) deallocate(irecv)

deallocate(isend,nrecv)

end subroutine derefinement

!        ------------------
function pcoarsen_indicator(grid,elem)
!        ------------------

!----------------------------------------------------
! This routine computes the p coarsening error indicator for element elem
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: elem
real(my_real) :: pcoarsen_indicator
!----------------------------------------------------
! Local variables:

integer :: i, j, k, deg
real(my_real) :: temp
!----------------------------------------------------
! Begin executable code

! Compute the l2 norm of the p-hierarchic coefficients of the current
! degree of this element.  Don't worry about whether the edges are actually
! of the same degree.
! For systems or multiple eigenvectors, use the max of the individual solutions.

deg = grid%element(elem)%degree

if (deg <= 1) then
   pcoarsen_indicator = huge(0.0_my_real)
   return
endif

pcoarsen_indicator = 0

do j=1,grid%system_size
 do k=1,max(1,grid%num_eval)
   temp = 0
   do i = ((deg-2)*(deg-3))/2 + 1, ((deg-1)*(deg-2))/2
      temp = temp + grid%element(elem)%solution(i,j,k)**2
   end do
   do i=1,EDGES_PER_ELEMENT
      if (grid%edge(grid%element(elem)%edge(i))%degree >= deg) then
         temp = temp + grid%edge(grid%element(elem)%edge(i))%solution(deg-1,j,k)**2
      endif
   end do
   pcoarsen_indicator = max(pcoarsen_indicator,sqrt(temp))
 end do
end do

end function pcoarsen_indicator

!        ------------------
function hcoarsen_indicator(grid,parent,children)
!        ------------------

!----------------------------------------------------
! This routine computes the h coarsening error indicator for the elements in
! children, which should be siblings that are children of parent and its mate.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: parent,children(:)
real(my_real) :: hcoarsen_indicator
!----------------------------------------------------
! Local variables:

real(my_real), allocatable :: solution(:,:)
real(my_real) :: xvert(4),yvert(4),temp
integer :: i,j,k,degree,astat,isub
!----------------------------------------------------
! Begin executable code

! set degree to be the maximum degree, and allocate space for the local
! copy of the solution coefficients

degree = 0
do i=1,4
  if (i==3 .and. children(3)==NO_CHILD) exit
  do j=1,EDGES_PER_ELEMENT
     degree = max(degree,grid%edge(grid%element(children(i))%edge(j))%degree)
  end do
end do

allocate(solution((degree+1)**2+degree**2,grid%nsoln),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in hcoarsen_indicator")
   return
endif

! copy solution coefficients to the local variable.  see subroutine
! phier2nodal for the order of the solution components

solution = 0
solution(1,:) = reshape( &
                grid%vertex_solution(grid%element(children(1))%vertex(1),:,:), &
                (/grid%nsoln/))
solution(2,:) = reshape( &
                grid%vertex_solution(grid%element(children(2))%vertex(1),:,:), &
                (/grid%nsoln/))
solution(3,:) = reshape( &
                grid%vertex_solution(grid%element(children(1))%vertex(2),:,:), &
                (/grid%nsoln/))
if (children(3) /= NO_CHILD) then
   solution(4,:) = reshape( &
                grid%vertex_solution(grid%element(children(3))%vertex(2),:,:), &
                (/grid%nsoln/))
endif
isub = 5
do i=1,degree-1
   if (grid%edge(grid%element(children(1))%edge(1))%degree >= i+1) then
      solution(isub,:) = reshape(grid%edge(grid%element(children(1))%edge(1))%solution(i,:,:),(/grid%nsoln/))
   endif
   isub = isub+1
end do
solution(isub,:) = reshape( &
                grid%vertex_solution(grid%element(children(1))%vertex(3),:,:), &
                (/grid%nsoln/))
isub = isub+1
do i=1,degree-1
   if (grid%edge(grid%element(children(1))%edge(2))%degree >= i+1) then
      solution(isub,:) = reshape(grid%edge(grid%element(children(1))%edge(2))%solution(i,:,:),(/grid%nsoln/))
   endif
   isub = isub + 1
end do
do i=1,degree-1
   if (grid%edge(grid%element(children(1))%edge(3))%degree >= i+1) then
      solution(isub,:) = reshape(grid%edge(grid%element(children(1))%edge(3))%solution(i,:,:),(/grid%nsoln/))
   endif
   isub = isub + 1
end do
k = 1
do j=1,degree-2
   do i=1,j
      if (grid%element(children(1))%degree >= j+2) then
         solution(isub,:) = reshape(grid%element(children(1))%solution(k,:,:),(/grid%nsoln/))
      endif
      isub = isub + 1
      k = k+1
   end do
end do
do i=1,degree-1
   if (grid%edge(grid%element(children(2))%edge(2))%degree >= i+1) then
      solution(isub,:) = reshape(grid%edge(grid%element(children(2))%edge(2))%solution(i,:,:),(/grid%nsoln/))
   endif
   isub = isub + 1
end do
do i=1,degree-1
   if (grid%edge(grid%element(children(2))%edge(3))%degree >= i+1) then
      solution(isub,:) = reshape(grid%edge(grid%element(children(2))%edge(3))%solution(i,:,:),(/grid%nsoln/))
   endif
   isub = isub + 1
end do
k = 1
do j=1,degree-2
   do i=1,j
      if (grid%element(children(2))%degree >= j+2) then
         solution(isub,:) = reshape(grid%element(children(2))%solution(k,:,:),(/grid%nsoln/))
      endif
      isub = isub + 1
      k = k+1
   end do
end do
if (children(3) /= NO_CHILD) then
   do i=1,degree-1
      if (grid%edge(grid%element(children(3))%edge(1))%degree >= i+1) then
         solution(isub,:) = reshape(grid%edge(grid%element(children(3))%edge(1))%solution(i,:,:),(/grid%nsoln/))
      endif
      isub = isub + 1
   end do
   do i=1,degree-1
      if (grid%edge(grid%element(children(3))%edge(3))%degree >= i+1) then
         solution(isub,:) = reshape(grid%edge(grid%element(children(3))%edge(3))%solution(i,:,:),(/grid%nsoln/))
      endif
      isub = isub + 1
   end do
   k = 1
   do j=1,degree-2
      do i=1,j
         if (grid%element(children(3))%degree >= j+2) then
            solution(isub,:) = reshape(grid%element(children(3))%solution(k,:,:),(/grid%nsoln/))
         endif
         isub = isub + 1
         k = k+1
      end do
   end do
   do i=1,degree-1
      if (grid%edge(grid%element(children(4))%edge(3))%degree >= i+1) then
         solution(isub,:) = reshape(grid%edge(grid%element(children(4))%edge(3))%solution(i,:,:),(/grid%nsoln/))
      endif
      isub = isub + 1
   end do
   k = 1
   do j=1,degree-2
      do i=1,j
         if (grid%element(children(4))%degree >= j+2) then
            solution(isub,:) = reshape(grid%element(children(4))%solution(k,:,:),(/grid%nsoln/))
         endif
         isub = isub + 1
         k = k+1
      end do
   end do
endif

! set the outer vertices of the parents

xvert(1) = grid%vertex(grid%element(parent)%vertex(1))%coord%x
yvert(1) = grid%vertex(grid%element(parent)%vertex(1))%coord%y
xvert(2) = grid%vertex(grid%element(parent)%vertex(2))%coord%x
yvert(2) = grid%vertex(grid%element(parent)%vertex(2))%coord%y
xvert(3) = grid%vertex(grid%element(parent)%vertex(3))%coord%x
yvert(3) = grid%vertex(grid%element(parent)%vertex(3))%coord%y
if (children(3) == NO_CHILD) then
   xvert(4) = huge(0.0_my_real)
   yvert(4) = huge(0.0_my_real)
else
   xvert(4) = grid%vertex(grid%element(children(3))%vertex(2))%coord%x
   yvert(4) = grid%vertex(grid%element(children(3))%vertex(2))%coord%y
endif

! convert solution to nodal basis coefficients, and then to h-hierarchical
! basis coefficients

call phier2nodal(solution,xvert,yvert,degree)
call nodal2hhier(solution,xvert,yvert,degree)

! compute the l2 norm of the red node h-hierarchic coefficients; for systems
! and multiple eigenvectors use the max over all solutions
! See subroutine phier2nodal for the order of the nodes

hcoarsen_indicator = 0

do j=1,grid%nsoln
   temp = 0
! along edge between children(1) and children(2), including central node
   isub = 5
   do i=isub,isub+2*((degree-1)/2),2
      temp = temp + solution(i,j)**2
   end do
   if (degree > 1) then
! along edge between children(1) and children(3)
      isub = degree+5
      do i=isub,isub+2*((degree-2)/2),2
         temp = temp + solution(i,j)**2
      end do
! along edge between children(2) and children(4)
      isub = 3+3*degree + ((degree-1)*(degree-2))/2
      do i=isub,isub+2*((degree-2)/2),2
         temp = temp + solution(i,j)**2
      end do
! along edge between children(3) and children(4)
      if (children(3) /= NO_CHILD) then
         isub = 1 + 5*degree + (degree-1)*(degree-2)
         do i=isub,isub+2*((degree-2)/2),2
            temp = temp + solution(i,j)**2
         end do
      endif
      if (degree > 2) then
! interior of children(1)
         isub = 6 + 3*(degree-1)
         do k=1,(degree-1)/2
            do i=1,degree-2*k
               temp = temp + solution(isub,j)**2
               isub = isub + 1
            end do
            isub = isub + degree - 2*k - 1
         end do
! interior of children(2)
         isub = 1 + 5*degree + ((degree-1)*(degree-2))/2
         do k=1,(degree-1)/2
            do i=1,degree-2*k
               temp = temp + solution(isub,j)**2
               isub = isub + 1
            end do
            isub = isub + degree - 2*k - 1
         end do
         if (children(3) /= NO_CHILD) then
! interior of children(3)
            isub = -1 + 7*degree + (degree-1)*(degree-2)
            do k=1,(degree-1)/2
               do i=1,degree-2*k
                  temp = temp + solution(isub,j)**2
                  isub = isub + 1
               end do
               isub = isub + degree - 2*k - 1
            end do
! interior of children(4)
            isub = -2 + 8*degree + 3*((degree-1)*(degree-2))/2
            do k=1,(degree-1)/2
               do i=1,degree-2*k
                  temp = temp + solution(isub,j)**2
                  isub = isub + 1
               end do
               isub = isub + degree - 2*k - 1
            end do
         endif ! children(3) /= NO_CHILD
      endif ! degree > 2
   endif ! degree > 1

   hcoarsen_indicator = max(hcoarsen_indicator,sqrt(temp))
end do

deallocate(solution)

end function hcoarsen_indicator

!          ----------
subroutine make_elist(grid,procs,elist,refcont,global_max_errind, &
                      numhref,numpref,predictive_load_balance,still_sequential)
!          ----------

!----------------------------------------------------
! This routine groups the leaf elements into size(head_errind) lists.
! List 1 contains those with error indicator between maxerrind and
! maxerrind/binw, list 2 between maxerrind/binw and maxerrind/binw**2, etc.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type (proc_info), intent(in) :: procs
type(errind_list), intent(inout) :: elist
type(refine_options), intent(in) :: refcont
real(my_real), intent(in) :: global_max_errind
integer, intent(in) :: numhref(:), numpref(:)
logical, intent(in) :: predictive_load_balance, still_sequential

!----------------------------------------------------
! Local variables:

real(my_real), allocatable :: errind(:)
integer :: lev, elem, i, astat, total
real(my_real) :: max_errind, reftol, normsoln
logical :: is_uniform

!----------------------------------------------------
! Begin executable code

allocate(errind(size(grid%element)),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in make_elist",procs=procs)
   return
endif
errind = 0.0_my_real

! set error indicators and find maximum

do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (grid%element(elem)%iown) then
         if (grid%element(elem)%isleaf) then
            errind(elem) = &
                  maxval(grid%element_errind(elem,:))/grid%element(elem)%work
         endif
      endif
      elem = grid%element(elem)%next
   end do
end do

! use global max errind from before derefinement for computing bins

elist%max_errind = global_max_errind
max_errind = global_max_errind

if (refcont%refterm == ONE_REF) then
   call get_grid_info(grid,procs,still_sequential,3344, &
                      total_nelem_leaf=total,no_master=.true.)
   reftol = refcont%reftol/sqrt(real(total,my_real))
! TEMP only using first component of first eigenvalue
   if (grid%errtype == RELATIVE_ERROR .and. .not. (refcont%reftype == HP_ADAPTIVE .and. &
    (refcont%hp_strategy == HP_T3S .or. refcont%hp_strategy == HP_ALTERNATE))) then
      call norm_solution(grid,procs,still_sequential,1,1,energy=normsoln)
      if (normsoln /= 0.0_my_real) reftol = reftol*normsoln
   endif
endif
if (refcont%refterm == ONE_REF_HALF_ERRIND) then
   reftol = global_max_errind/refcont%inc_factor
endif

! create lists

elist%head_errind = END_OF_LIST
elist%tail_errind = END_OF_LIST
elist%next_errind = NOT_ON_LIST
elist%prev_errind = NOT_ON_LIST

! determine if this is uniform refinement

is_uniform = .false.
if (global_max_errind == 0.0_my_real) is_uniform = .true.
if (refcont%reftype == HP_ADAPTIVE .and. &
    refcont%hp_strategy == HP_T3S) then
   if (refcont%t3s_reftype == H_UNIFORM) is_uniform = .true.
endif

do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (grid%element(elem)%isleaf) then
! if numhref or numpref is positive, it goes in the first bin
         if (numhref(elem) > 0 .or. numpref(elem) > 0) then
            i = 1
! if both numhref and numpref are 0, it goes in the last bin
         elseif (numhref(elem) == 0 .and. numpref(elem) == 0) then
            i = size(elist%head_errind)
! uniform refinement puts everything in first bin
         elseif (is_uniform) then
            i = 1
         else
            select case (refcont%refterm)
! ONE_REF put large error indicators in first bin and everything else in second
            case (ONE_REF, ONE_REF_HALF_ERRIND)
               if (grid%element(elem)%iown .and. &
                   maxval(grid%element_errind(elem,:)) > reftol) then
                  i = 1
               else
                  i = 2
               endif
            case default
! otherwise find the right bin
               do i=1,size(elist%head_errind)-1
                  if (errind(elem) > max_errind/(binw**i)) exit
               end do
            end select
         endif
         if (elist%head_errind(i) == END_OF_LIST) then
            elist%head_errind(i) = elem
            elist%tail_errind(i) = elem
            elist%next_errind(elem) = END_OF_LIST
            elist%prev_errind(elem) = END_OF_LIST
         else
            elist%next_errind(elem) = elist%head_errind(i)
            elist%prev_errind(elist%head_errind(i)) = elem
            elist%head_errind(i) = elem
            elist%prev_errind(elem) = END_OF_LIST
         endif
      endif
      elem = grid%element(elem)%next
   end do
end do
deallocate(errind,stat=astat)

end subroutine make_elist

!          ----------------
subroutine set_loop_control(refine_control,return_to_elist,complete_elist, &
                            one_elist)
!          ----------------

!----------------------------------------------------
! This routine sets the variables that control the refine loop:
! return_to_elist, tells if children get put into an elist
! complete_elist, tells if you must finish the whole elist (SMOOTH)
! one_elist, tells if you quit after finishing the first elist
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(refine_options), intent(in) :: refine_control
logical, intent(out) :: return_to_elist, complete_elist, one_elist
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

select case (refine_control%reftype)
case (H_ADAPTIVE, P_ADAPTIVE)
   select case (refine_control%refterm)
   case (DOUBLE_NVERT, DOUBLE_NELEM, DOUBLE_NEQ, HALVE_ERREST, KEEP_NVERT, &
         KEEP_NELEM, KEEP_NEQ, KEEP_ERREST)
      return_to_elist = .true.
      complete_elist = .false.
      one_elist = .false.
   case (DOUBLE_NVERT_SMOOTH, DOUBLE_NELEM_SMOOTH, DOUBLE_NEQ_SMOOTH, &
         KEEP_NVERT_SMOOTH, KEEP_NELEM_SMOOTH, KEEP_NEQ_SMOOTH)
      return_to_elist = .true.
      complete_elist = .true.
      one_elist = .false.
   case (ONE_REF, ONE_REF_HALF_ERRIND)
      return_to_elist = .false.
      complete_elist = .true.
      one_elist = .true.
   end select
case (HP_ADAPTIVE)
   select case (refine_control%hp_strategy)
   case (HP_BIGGER_ERRIND, HP_APRIORI, HP_PRIOR2P_E, HP_PRIOR2P_H1, &
         HP_TYPEPARAM, HP_COEF_DECAY, HP_COEF_ROOT, HP_SMOOTH_PRED, &
         HP_NEXT3P)
      select case (refine_control%refterm)
      case (DOUBLE_NVERT, DOUBLE_NELEM, DOUBLE_NEQ, HALVE_ERREST, KEEP_NVERT, &
            KEEP_NELEM, KEEP_NEQ, KEEP_ERREST)
         return_to_elist = .true.
         complete_elist = .false.
         one_elist = .false.
      case (DOUBLE_NVERT_SMOOTH, DOUBLE_NELEM_SMOOTH, DOUBLE_NEQ_SMOOTH, &
            KEEP_NVERT_SMOOTH, KEEP_NELEM_SMOOTH, KEEP_NEQ_SMOOTH)
         return_to_elist = .true.
         complete_elist = .true.
         one_elist = .false.
      case (ONE_REF, ONE_REF_HALF_ERRIND)
         return_to_elist = .false.
         complete_elist = .true.
         one_elist = .true.
      end select
   case (HP_T3S)
      return_to_elist = .true.
      complete_elist = .true.
      one_elist = .true.
   case (HP_ALTERNATE)
      return_to_elist = .false.
      complete_elist = .true.
      one_elist = .true.
   case default
      ierr = PHAML_INTERNAL_ERROR
      call fatal("unexpected value for hp_strategy in set_loop_control")
      stop
   end select
case default
   ierr = PHAML_INTERNAL_ERROR
   call fatal("unexpected value for reftype in set_loop_control")
   stop
end select

end subroutine set_loop_control

!        ------------
function check_target(grid,refine_control,elist,target)
!        ------------

!----------------------------------------------------
! This routine checks to see if the target for terminating refinement is met
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
type(refine_options), intent(in) :: refine_control
type(errind_list), intent(in) :: elist
integer, intent(in) :: target
logical :: check_target
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

check_target = .false.

select case (refine_control%refterm)
case (DOUBLE_NVERT, KEEP_NVERT, DOUBLE_NVERT_SMOOTH, KEEP_NVERT_SMOOTH)
   if (grid%nvert_own >= target) check_target = .true.
case (DOUBLE_NELEM, KEEP_NELEM, DOUBLE_NELEM_SMOOTH, KEEP_NELEM_SMOOTH)
   if (grid%nelem_leaf_own >= target) check_target = .true.
case (DOUBLE_NEQ, KEEP_NEQ, DOUBLE_NEQ_SMOOTH, KEEP_NEQ_SMOOTH)
   if (grid%dof_own >= target) check_target = .true.
case (HALVE_ERREST, KEEP_ERREST)
   if (elist%current_list > target) check_target = .true.
case (ONE_REF, ONE_REF_HALF_ERRIND)
   check_target = .false.
case default
   call fatal("illegal value for refterm")
   stop
end select

end function check_target

end module refine_adapt_mod
