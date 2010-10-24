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

module refine_top

!----------------------------------------------------
! This module contains the top level refinement routines
!
! communication tags in this module are of the form 34xx
!----------------------------------------------------

!----------------------------------------------------
! Other modules used:

use global
use gridtype_mod
use grid_util
use linsystype_mod
use message_passing
use refine_uniform_mod
use refine_adapt_mod
use refine_elements
use hp_strategies
use hash_mod
!----------------------------------------------------

implicit none
private
public refine, reconcile

!----------------------------------------------------
! Non-module procedures used are:

!----------------------------------------------------
! The following defined constants are defined:

!----------------------------------------------------
! The following variables are defined:

!----------------------------------------------------

contains

!          ------
subroutine refine(grid,procs,refine_control,solver_control,io_control, &
                  still_sequential,init_nvert,init_nelem,init_dof,loop, &
                  balance_what,predictive,no_time)
!          ------

!----------------------------------------------------
! Top level refinement routine.  It mainly determines which refinement
! routine to call.
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

logical :: timeit
!----------------------------------------------------
! Begin executable code

! start timing the refinement process

if (present(no_time)) then
   timeit = .not. no_time
else
   timeit = .true.
endif
if (timeit) then
   call reset_watch(prefine)
   call start_watch((/prefine,trefine/))
endif

! The REFSOLN_EDGE and REFSOLN_ELEM hp-adaptive strategies have their own
! refine routine

if (refine_control%reftype == HP_ADAPTIVE .and. &
    (refine_control%hp_strategy == HP_REFSOLN_EDGE .or. &
     refine_control%hp_strategy == HP_REFSOLN_ELEM)) then
   call refine_refsoln(grid,procs,refine_control,solver_control, &
                       io_control,still_sequential,init_nvert,init_nelem, &
                       init_dof,loop,balance_what,predictive)

! NLP also has it's own routine

elseif (refine_control%reftype == HP_ADAPTIVE .and. &
        refine_control%hp_strategy == HP_NLP) then
   call refine_nlp(grid,procs,refine_control,still_sequential)

! uniform refinement

elseif (refine_control%reftype == P_UNIFORM) then
   call refine_uniform_p(grid,refine_control)
elseif (refine_control%reftype == H_UNIFORM) then
   call refine_uniform_h(grid,refine_control)

! all other adaptive refinements

else
   call refine_adaptive(grid,procs,refine_control,solver_control,io_control, &
                        still_sequential,init_nvert,init_nelem,init_dof,loop, &
                        balance_what,predictive,no_time)

endif

! error indicators should be recomputed and say we need to resend graphics data

grid%errind_up2date = .false.
grid_changed = .true.

! stop timing the refinement process

if (timeit) call stop_watch((/prefine,trefine/))

end subroutine refine

!          ---------
subroutine reconcile(grid,procs,refine_control,sequent)
!          ---------

!----------------------------------------------------
! This routine reconciles the refinements among the partitions
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(proc_info), intent(in) :: procs
type(refine_options), intent(in) :: refine_control
logical, intent(in) :: sequent

!----------------------------------------------------
! Local variables:

integer :: lev, elem, count, part, i, j, elemlid, errcode, vertlid, &
           astat, nproc
! newcomm
integer :: nsend
integer, allocatable :: isend(:), nrecv(:)
integer, pointer :: irecv(:)
type(hash_key) :: gid
logical, allocatable :: isboundary(:)
logical :: doit

!----------------------------------------------------
! Begin executable code

if (sequent .or. num_proc(procs) == 1) return

grid_changed = .true.

if (my_proc(procs) == MASTER) return

! start timing the reconciliation process

call reset_watch((/precon,cprecon/))
call start_watch((/precon,trecon/))

nproc = num_proc(procs)

allocate(nrecv(nproc),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in reconcile",procs=procs)
   return
endif

! make a list of the unowned elements I have refined along with the degree
! if it was p refined or -1 if it was h refined

count = 0
do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (grid%element(elem)%hrefined_unowned) count = count + 1
      if (grid%element(elem)%prefined_unowned) count = count + 1
      elem = grid%element(elem)%next
   end do
end do

nsend = (1+KEY_SIZE)*count
allocate(isend(nsend),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in reconcile",procs=procs)
   return
endif

count = 1
do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (grid%element(elem)%hrefined_unowned) then
         call hash_pack_key(grid%element(elem)%gid,isend,count)
         isend(count+KEY_SIZE) = -1
         grid%element(elem)%hrefined_unowned = .false.
         count = count + KEY_SIZE+1
      endif
      if (grid%element(elem)%prefined_unowned) then
         call hash_pack_key(grid%element(elem)%gid,isend,count)
         isend(count+KEY_SIZE) = grid%element(elem)%degree
         grid%element(elem)%prefined_unowned = .false.
         count = count + KEY_SIZE+1
      endif
      elem = grid%element(elem)%next
   end do
end do
 
! exchange the list with other partitions

call start_watch((/cprecon,ctrecon/))
call phaml_alltoall(procs,isend,nsend,irecv,nrecv,3470)
call stop_watch((/cprecon,ctrecon/))

errcode = 0
deallocate(isend,stat=astat)

! scan the list looking for elements that I own and have not refined,
! and refine those elements

count = 1
outer: do part=1,nproc
   do i=1,nrecv(part)/(KEY_SIZE+1)
      elemlid = hash_decode_key(hash_unpack_key(irecv,count),grid%elem_hash)
      if (elemlid == HASH_NOT_FOUND) then
         count = count + KEY_SIZE+1
         cycle
      endif
      if (.not. grid%element(elemlid)%iown) then
         count = count + KEY_SIZE+1
         cycle
      endif
      if (irecv(count+KEY_SIZE) == -1) then
         if (grid%element(elemlid)%isleaf) then
            call bisect_triangle_pair(grid,elemlid,errcode,refine_control)
            if (errcode /= 0) then
               call warning("refinement failed during reconciliation")
               exit outer
            endif
         endif
      else
         if (grid%element(elemlid)%isleaf) then
            do while (grid%element(elemlid)%degree < irecv(count+KEY_SIZE))
               call p_refine_elem(grid,elemlid,refine_control)
            end do
         endif
      endif
      count = count + KEY_SIZE+1
   end do
end do outer

if (associated(irecv)) deallocate(irecv,stat=astat)

! if doing adaptive p, p coarsen each unowned leaf as much as possible
! TEMP I may be p coarsening elements that need to be p refined for overlap,
!      which adds extra work and looses the solution

if (refine_control%reftype == P_ADAPTIVE .or. &
    refine_control%reftype == HP_ADAPTIVE) then
   do lev=1,grid%nlev
      elem = grid%head_level_elem(lev)
      do while (elem /= END_OF_LIST)
         if (grid%element(elem)%isleaf .and. &
             (.not. grid%element(elem)%iown)) then
            errcode = 0
            do while (errcode == 0)
               call p_coarsen_elem(grid,elem,errcode,refine_control)
            end do
         endif
         elem = grid%element(elem)%next
      end do
   end do
endif

! enforce overlap so that any element that touches my partition from the
! outside exists and has high enough degrees to match the one that owns it

! mark each vertex as being on the partition boundary or not by going through
! the elements and marking each vertex as being on the boundary if I own the
! element or vertex but not the other

allocate(isboundary(size(grid%vertex)),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in reconcile",procs=procs)
   return
endif
isboundary = .false.

do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      do i=1,VERTICES_PER_ELEMENT
         if (grid%element(elem)%iown .neqv. &
             grid%element(grid%vertex(grid%element(elem)%vertex(i))%assoc_elem)%iown) then
            isboundary(grid%element(elem)%vertex(i)) = .true.
         endif
      end do
      elem = grid%element(elem)%next
   end do
end do

! make a list of owned leaves that touch the partition boundary by going
! through the elements and finding those that have a marked vertex.  Send
! the element along with the vertices and the element and edge degrees

count = 0
do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (grid%element(elem)%iown .and. grid%element(elem)%isleaf) then
         do i=1,VERTICES_PER_ELEMENT
            if (isboundary(grid%element(elem)%vertex(i))) then
               count = count + 1
               exit
            endif
         end do
      endif
      elem = grid%element(elem)%next
   end do
end do

nsend=4*count*(KEY_SIZE+1)
allocate(isend(nsend),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in reconcile",procs=procs)
   return
endif

count = 1
do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (grid%element(elem)%iown .and. grid%element(elem)%isleaf) then
         do i=1,VERTICES_PER_ELEMENT
            if (isboundary(grid%element(elem)%vertex(i))) then
               call hash_pack_key(grid%element(elem)%gid,isend,count)
               call hash_pack_key(grid%vertex(grid%element(elem)%vertex(1))%gid, &
                                  isend,count+KEY_SIZE)
               call hash_pack_key(grid%vertex(grid%element(elem)%vertex(2))%gid, &
                                  isend,count+2*KEY_SIZE)
               call hash_pack_key(grid%vertex(grid%element(elem)%vertex(3))%gid, &
                                  isend,count+3*KEY_SIZE)
               isend(count+4*KEY_SIZE) = grid%element(elem)%degree
               isend(count+4*KEY_SIZE+1:count+4*KEY_SIZE+3) = &
                                  grid%edge(grid%element(elem)%edge)%degree
               count = count + 4*(KEY_SIZE+1)
               exit
            endif
         end do
      endif
      elem = grid%element(elem)%next
   end do
end do

! exchange lists with other processors

call start_watch((/cprecon,ctrecon/))
call phaml_alltoall(procs,isend,nsend,irecv,nrecv,3480)
call stop_watch((/cprecon,ctrecon/))

deallocate(isend,stat=astat)

! go through the received lists of elements.  If any vertex of the element is
! marked as being on the partition boundary, h refine to create the element if
! it doesn't exist, and p refine if the element or edge degrees are too small

count = 1
do part=1,nproc
   do i=1,nrecv(part)/(4*(KEY_SIZE+1))
      doit = .false.
      do j=1,VERTICES_PER_ELEMENT
         vertlid = hash_decode_key(hash_unpack_key(irecv,count+j*KEY_SIZE), &
                                   grid%vert_hash)
         if (vertlid == HASH_NOT_FOUND) cycle
         if (isboundary(vertlid)) then
            doit = .true.
            exit
         endif
      end do
      if (doit) then
         gid = hash_unpack_key(irecv,count)
         elemlid = hash_decode_key(gid,grid%elem_hash)
         if (elemlid == HASH_NOT_FOUND) then
            call create_element(grid,gid,refine_control,errcode)
            if (errcode /= 0) then
               call warning("refinement for overlap failed")
               count = count + 4*(KEY_SIZE+1)
               cycle
            endif
            elemlid = hash_decode_key(gid,grid%elem_hash)
         endif
         if (grid%element(elemlid)%isleaf) then
            do while (grid%element(elemlid)%degree < irecv(count+4*KEY_SIZE))
               call p_refine_elem(grid,elemlid,refine_control)
            end do
         endif
         do j=1,EDGES_PER_ELEMENT
            if (grid%edge(grid%element(elemlid)%edge(j))%degree < &
                irecv(count+4*KEY_SIZE+j)) then
               call adjust_edge_degree(grid,grid%element(elemlid)%edge(j), &
                                       irecv(count+4*KEY_SIZE+j))
            endif
         end do
      endif
      count = count + 4*(KEY_SIZE+1)
   end do
end do

if (associated(irecv)) deallocate(irecv,stat=astat)
deallocate(isboundary,stat=astat)

deallocate(nrecv,stat=astat)

grid%errind_up2date = .false.

! stop timing the reconciliation process

call stop_watch((/precon,trecon/))

end subroutine reconcile

!          ------------------
subroutine adjust_edge_degree(grid,edge,newdeg)
!          ------------------

!----------------------------------------------------
! This routine sets the edge degree to the given value.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: edge, newdeg
!----------------------------------------------------
! Local variables:

integer :: olddeg, edge_size, oldsize, astat, d1, d2, d3, i
real(my_real), pointer :: temp1(:,:,:)
!----------------------------------------------------
! Begin executable code

olddeg = grid%edge(edge)%degree

! make sure allocated memory is large enough.

edge_size = newdeg-1
if (associated(grid%edge(edge)%solution)) then
   oldsize = size(grid%edge(edge)%solution,dim=1)
else
   oldsize = 0
endif
if (oldsize < edge_size) then
   allocate(temp1(edge_size,grid%system_size,max(1,grid%num_eval)), &
            stat=astat)
   if (astat /= 0) then
      call fatal("allocation failed in adjust_edge_degree")
      stop
   endif
   temp1 = 0.0_my_real
   if (oldsize > 0) temp1(1:oldsize,:,:) = grid%edge(edge)%solution
   deallocate(grid%edge(edge)%solution, stat=astat)
   grid%edge(edge)%solution => temp1
   if (grid%have_true) then
      nullify(temp1)
      allocate(temp1(edge_size,grid%system_size,max(1,grid%num_eval)), &
               stat=astat)
      if (astat /= 0) then
         call fatal("allocation failed in adjust_edge_degree")
         stop
      endif
      temp1 = 0.0_my_real
      if (oldsize > 0) temp1(1:oldsize,:,:) = grid%edge(edge)%exact
      deallocate(grid%edge(edge)%exact, stat=astat)
      grid%edge(edge)%exact => temp1
   endif
endif
if (grid%oldsoln_exists) then
   if (associated(grid%edge(edge)%oldsoln)) then
      nullify(temp1)
      allocate(temp1(edge_size,grid%system_size,max(1,grid%num_eval)),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in adjust_edge_degree")
         stop
      endif
      temp1 = 0.0_my_real
      d1 = min(size(grid%edge(edge)%oldsoln,dim=1),size(temp1,dim=1))
      d2 = min(size(grid%edge(edge)%oldsoln,dim=2),size(temp1,dim=2))
      d3 = min(size(grid%edge(edge)%oldsoln,dim=3),size(temp1,dim=3))
      temp1(1:d1,1:d2,1:d3) = grid%edge(edge)%oldsoln(1:d1,1:d2,1:d3)
      deallocate(grid%edge(edge)%oldsoln)
      grid%edge(edge)%oldsoln => temp1
   endif
endif

! Set the edge degree.  Also set Dirichlet boundary conditions on Dirichlet
! edges that change degree and set solution to 0 for other new components

grid%edge(edge)%degree = newdeg
grid%dof = grid%dof + grid%system_size*(newdeg - olddeg)
! TEMP should grid%dof_own be adjusted here, too?
if (grid%have_true) then
   grid%edge(edge)%exact = 0.0_my_real
endif
do i=1,grid%system_size
   if (grid%edge_type(edge,i) == DIRICHLET) then
      call edge_exact(grid,edge,i,"d")
   else
      grid%edge(edge)%solution(newdeg-1,i,:) = 0.0_my_real
   endif
   if (grid%have_true) call edge_exact(grid,edge,i,"t")
end do

end subroutine adjust_edge_degree

end module refine_top
