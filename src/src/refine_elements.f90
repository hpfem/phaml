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

module refine_elements

!----------------------------------------------------
! This module contains routines with the details of element refinement.
!----------------------------------------------------

!----------------------------------------------------
! Other modules used:

use global
use gridtype_mod
use hash_mod
use grid_util
use error_estimators
use make_linsys
use message_passing, only: fatal, warning
!----------------------------------------------------

implicit none
private
public p_coarsen_elem, &            ! not thread safe; TEMP to be made so
       unbisect_triangle_pair,   &  ! not thread safe; TEMP to be made so
       p_refine_elem, &             ! not thread safe
       p_refine_element_interior, & ! conditionally thread safe
       enforce_edge_rule, &         ! conditionally thread safe
       before_h_refine, &           ! not thread safe
       bisect_triangle_pair, &      ! conditionally thread safe
       after_h_refine, &            ! not thread safe
       remove_from_errind_list, &   ! not thread safe
       create_element, &            ! not thread safe
       init_guess_ic, &             ! not thread safe
       init_guess_p                 ! not thread safe

!----------------------------------------------------
! Non-module procedures used are:

interface

   function trues(x,y,comp,eigen) ! real (my_real)
   use global
   real (my_real), intent(in) :: x,y
   integer, intent(in) :: comp, eigen
   real (my_real) :: trues
   end function trues

   subroutine bconds(x,y,bmark,itype,c,rs)
   use global
   real(my_real), intent(in) :: x,y
   integer, intent(in) :: bmark
   integer, intent(out) :: itype(:)
   real(my_real), intent(out) :: c(:,:),rs(:)
   end subroutine bconds

   function iconds(x,y,comp,eigen)
   use global
   real(my_real), intent(in) :: x,y
   integer, intent(in) :: comp,eigen
   real(my_real) :: iconds
   end function iconds

end interface

!----------------------------------------------------
! The following defined constants are defined:

!----------------------------------------------------
! The following variables are defined:

!----------------------------------------------------

contains

!          --------------
subroutine p_coarsen_elem(grid,elem,errcode,refine_control)
!          --------------

!----------------------------------------------------
! This routine reduces the degree of element elem by 1.
! errcode is 0 if successful
!            1 if an error occurred
!           -1 if cannot coarsen
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: elem
integer, intent(out) :: errcode
type(refine_options), intent(in) :: refine_control
!----------------------------------------------------
! Local variables:

integer :: j, deg, side, edge, neigh(EDGES_PER_ELEMENT), new_deg, &
           old_edge_deg(EDGES_PER_ELEMENT)
!----------------------------------------------------
! Begin executable code

errcode = 0

! if the degree is 1, cannot coarsen

if (grid%element(elem)%degree <= 1) then
   errcode = -1
   return
endif

! Reduce element degree and apply edge degree rule.
! Also set exact and Dirichlet boundary conditions on edges that change degree.

deg = grid%element(elem)%degree - 1
grid%element(elem)%degree = deg
if (deg >= 2) then
   grid%dof = grid%dof - grid%system_size*(deg - 1)
   if (grid%element(elem)%iown) grid%dof_own = grid%dof_own - &
                                               grid%system_size*(deg - 1)
endif

neigh = get_neighbors(elem,grid)
do side=1,EDGES_PER_ELEMENT
   edge = grid%element(elem)%edge(side)
   old_edge_deg(side) = grid%edge(edge)%degree
   if (neigh(side) == BOUNDARY) then
      new_deg = grid%element(elem)%degree
   elseif (refine_control%edge_rule == MINIMUM_RULE) then
      new_deg = min(grid%element(elem)%degree,grid%element(neigh(side))%degree)
   else ! (refine_control%edge_rule == MAXIMUM_RULE)
      new_deg = max(grid%element(elem)%degree,grid%element(neigh(side))%degree)
   endif
   if (new_deg /= old_edge_deg(side)) then
      grid%edge(edge)%degree = new_deg
      grid%dof = grid%dof + grid%system_size*(new_deg - old_edge_deg(side))
      if (grid%element(grid%edge(edge)%assoc_elem)%iown) then
         grid%dof_own = grid%dof_own + &
                        grid%system_size*(new_deg - old_edge_deg(side))
      endif
      if (grid%have_true) then
         grid%edge(edge)%exact = 0.0_my_real
      endif
      do j=1,grid%system_size
         if (grid%edge_type(edge,j) == DIRICHLET) then
            call edge_exact(grid,edge,j,"d")
         elseif (new_deg > old_edge_deg(side)) then
            grid%edge(edge)%solution(old_edge_deg(side):new_deg-1,j,:) = 0.0_my_real
         elseif (new_deg < old_edge_deg(side)) then
            grid%edge(edge)%solution(new_deg:old_edge_deg(side)-1,j,:) = 0.0_my_real
         endif
         if (grid%have_true) call edge_exact(grid,edge,j,"t")
      end do
   endif
end do ! next side

! fix exact for the face bases; other solution coefficients stay the same
! because the basis is p-hierarchical

if (grid%have_true) then
   do j=1,grid%system_size
      call elem_exact(grid,elem,j,"t")
   end do
endif

! compute error indicator

call error_indicator(grid,elem,refine_control%error_estimator, &
                     grid%element_errind(elem,:),grid%element(elem)%work)

end subroutine p_coarsen_elem

!          ----------------------
subroutine unbisect_triangle_pair(grid,parent,errcode,refcont)
!          ----------------------

!----------------------------------------------------
! This routine removes the bisection refinement of element parent and
! its mate.
! errcode is 0 if successful
!            1 if an error occurred
!           -1 if cannot unrefine because of grandchildren
!           -2 if nothing to derefine (no children)
!           -3 if cannot unrefine because children have different owners
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: parent
integer, intent(out) :: errcode
type(refine_options), intent(in) :: refcont

!----------------------------------------------------
! Local variables:

integer :: child(4), childedge(6), mate, i, j, vert1, vert2, par_deg, mate_deg,&
           old_par_deg, old_mate_deg, old_edge_deg(2*EDGES_PER_ELEMENT), &
           elem_size, oldsize, edge_size, astat, edge, edge_deg
real(my_real), pointer :: temp1(:,:,:)
!----------------------------------------------------
! Begin executable code

errcode = 0

! identify the children

child(1:2) = get_child_lid(grid%element(parent)%gid,(/1,2/),grid%elem_hash)

! if there are no children, nothing to derefine

if (child(1) == NO_CHILD) then
   errcode = -2
   return
endif

! if there are grandchildren, cannot derefine

if (.not.grid%element(child(1))%isleaf .or. &
    .not.grid%element(child(2))%isleaf) then
   errcode = -1
   return
endif

! prohibit derefining an element for which the children are owned by
! different partitions, because the other partition may not know that
! the element is no longer there for compatibility

if (grid%element(child(1))%iown .neqv. grid%element(child(2))%iown) then
   errcode = -3
   return
endif

! identify mate and mate's children, and check for grandchildren

if (grid%element(parent)%mate == BOUNDARY) then
   mate = BOUNDARY
else
   mate = hash_decode_key(grid%element(parent)%mate,grid%elem_hash)
   if (mate == HASH_NOT_FOUND) then
      errcode = 1
      call warning("mate does not exist for a triangle to be derefined")
      return
   endif
   child(3:4) = get_child_lid(grid%element(mate)%gid,(/1,2/),grid%elem_hash)
   if (.not.grid%element(child(3))%isleaf .or. &
       .not.grid%element(child(4))%isleaf) then
      errcode = -1
      return
   endif
   if (grid%element(child(3))%iown .neqv. grid%element(child(4))%iown) then
      errcode = -3
      return
   endif
endif

! identify the child edges

childedge(1) = grid%element(child(1))%edge(2)
childedge(2) = grid%element(child(2))%edge(2)
childedge(3) = grid%element(child(1))%edge(1)
if (mate /= BOUNDARY) then
   childedge(4) = grid%element(child(3))%edge(1)
   childedge(5) = grid%element(child(3))%edge(2)
   childedge(6) = grid%element(child(4))%edge(2)
endif

! determine if I own the unrefined elements

if (mate /= BOUNDARY) then
   if ((grid%element(child(3))%iown .and. grid%element(child(4))%iown) .neqv. grid%element(mate)%iown) then
      call warning("need to keep assignment of iown in unbisect_triangle_pair")
   endif
   grid%element(mate)%iown = grid%element(child(3))%iown .and. &
                             grid%element(child(4))%iown
endif

! reassign associated elements for the vertices

if (grid%vertex(grid%element(parent)%vertex(1))%assoc_elem == child(1)) then
   grid%vertex(grid%element(parent)%vertex(1))%assoc_elem = parent
endif
if (grid%vertex(grid%element(parent)%vertex(2))%assoc_elem == child(2)) then
   grid%vertex(grid%element(parent)%vertex(2))%assoc_elem = parent
endif
if (grid%vertex(grid%element(parent)%vertex(3))%assoc_elem == child(1)) then
   grid%vertex(grid%element(parent)%vertex(3))%assoc_elem = parent
endif
if (mate /= BOUNDARY) then
   if (grid%vertex(grid%element(mate)%vertex(1))%assoc_elem == child(3)) then
      grid%vertex(grid%element(mate)%vertex(1))%assoc_elem = mate
   endif
   if (grid%vertex(grid%element(mate)%vertex(2))%assoc_elem == child(4)) then
      grid%vertex(grid%element(mate)%vertex(2))%assoc_elem = mate
   endif
   if (grid%vertex(grid%element(mate)%vertex(3))%assoc_elem == child(3)) then
      grid%vertex(grid%element(mate)%vertex(3))%assoc_elem = mate
   endif
endif
do i=1,3
   if (any(grid%vertex_type(grid%element(parent)%vertex(i),:) == PERIODIC_SLAVE)) then
      grid%vertex(grid%element(parent)%vertex(i))%assoc_elem = &
         grid%vertex(grid%vertex(grid%element(parent)%vertex(i))%next)%assoc_elem
   endif
   if (any(grid%vertex_type(grid%element(parent)%vertex(i),:) == PERIODIC_MASTER)) then
      grid%vertex(grid%vertex(grid%element(parent)%vertex(i))%previous)%assoc_elem = &
         grid%vertex(grid%element(parent)%vertex(i))%assoc_elem
   endif
end do
if (mate /= BOUNDARY) then
 do i=1,3
   if (any(grid%vertex_type(grid%element(parent)%vertex(i),:) == PERIODIC_SLAVE)) then
      grid%vertex(grid%element(parent)%vertex(i))%assoc_elem = &
         grid%vertex(grid%vertex(grid%element(parent)%vertex(i))%next)%assoc_elem
   endif
   if (any(grid%vertex_type(grid%element(parent)%vertex(i),:) == PERIODIC_MASTER)) then
      grid%vertex(grid%vertex(grid%element(parent)%vertex(i))%previous)%assoc_elem = &
         grid%vertex(grid%element(parent)%vertex(i))%assoc_elem
   endif
 end do
endif

! reassign associated element for the edges of the parents

do i=1,EDGES_PER_ELEMENT
   if (grid%edge(grid%element(parent)%edge(i))%assoc_elem == child(1) .or. &
       grid%edge(grid%element(parent)%edge(i))%assoc_elem == child(2)) then
      grid%edge(grid%element(parent)%edge(i))%assoc_elem = parent
   endif
   if (mate /= BOUNDARY) then
      if (grid%edge(grid%element(mate)%edge(i))%assoc_elem == child(3) .or. &
          grid%edge(grid%element(mate)%edge(i))%assoc_elem == child(4)) then
         grid%edge(grid%element(mate)%edge(i))%assoc_elem = mate
      endif
   endif
end do

! determine the degree of the parent and mate

old_par_deg = grid%element(parent)%degree
par_deg = max(grid%element(child(1))%degree,grid%element(child(2))%degree)
if (mate == BOUNDARY) then
   old_mate_deg = old_par_deg
   mate_deg = par_deg
else
   old_mate_deg = grid%element(mate)%degree
   mate_deg = max(grid%element(child(3))%degree,grid%element(child(4))%degree)
endif

! set parent degree and make sure the parent has enough memory for solution

grid%element(parent)%degree = par_deg
elem_size = ((par_deg-1)*(par_deg-2))/2
if (associated(grid%element(parent)%solution)) then
   oldsize = size(grid%element(parent)%solution,dim=1)
else
   oldsize = 0
endif
if (oldsize < elem_size) then
   allocate(temp1(elem_size,grid%system_size,max(1,grid%num_eval)),stat=astat)
   if (astat /= 0) then
      call fatal("allocation failed in unbisect_triangle_pair")
      stop 
   endif
   temp1 = 0.0_my_real
   deallocate(grid%element(parent)%solution, stat=astat)
   grid%element(parent)%solution => temp1
   if (grid%have_true) then
      nullify(temp1)
      allocate(temp1(elem_size,grid%system_size,max(1,grid%num_eval)),stat=astat)
      if (astat /= 0) then
         call fatal("allocation failed in unbisect_triangle_pair")
         stop 
      endif
      temp1 = 0.0_my_real
      deallocate(grid%element(parent)%exact, stat=astat)
      grid%element(parent)%exact => temp1
   endif
endif

! set mate degree and make sure the mate has enough memory for solution

if (mate /= BOUNDARY) then
   grid%element(mate)%degree = mate_deg
   elem_size = ((mate_deg-1)*(mate_deg-2))/2
   if (associated(grid%element(mate)%solution)) then
      oldsize = size(grid%element(mate)%solution,dim=1)
   else
      oldsize = 0
   endif
   if (oldsize < elem_size) then
      allocate(temp1(elem_size,grid%system_size,max(1,grid%num_eval)), &
               stat=astat)
      if (astat /= 0) then
         call fatal("allocation failed in unbisect_triangle_pair")
         stop 
      endif
      temp1 = 0.0_my_real
      deallocate(grid%element(mate)%solution, stat=astat)
      grid%element(mate)%solution => temp1
      if (grid%have_true) then
         nullify(temp1)
         allocate(temp1(elem_size,grid%system_size,max(1,grid%num_eval)), &
                  stat=astat)
         if (astat /= 0) then
            call fatal("allocation failed in unbisect_triangle_pair")
            stop 
         endif
         temp1 = 0.0_my_real
         deallocate(grid%element(mate)%exact, stat=astat)
         grid%element(mate)%exact => temp1
      endif
   endif
endif

! for each edge, set the degree, and allocate and set solution and exact if
! necessary

do i=1,6
   if (i==4 .and. mate==BOUNDARY) exit
   if (i<=3) then
      edge = grid%element(parent)%edge(i)
   else
      edge = grid%element(mate)%edge(i-3)
   endif
   if (i<3) then
      edge_deg = par_deg
   elseif (i==3 .or. i==6) then
      select case (refcont%edge_rule)
      case (MINIMUM_RULE)
         edge_deg = min(par_deg,mate_deg)
      case (MAXIMUM_RULE)
         edge_deg = max(par_deg,mate_deg)
      end select
      grid%edge(edge)%degree = edge_deg
   else
      edge_deg = mate_deg
   endif
   if (i==3 .or. i==6) then
      old_edge_deg(i) = 1
   else
      old_edge_deg(i) = grid%edge(edge)%degree
   endif
   if (old_edge_deg(i) >= edge_deg) cycle
   edge_size = edge_deg-1
   if (associated(grid%edge(edge)%solution)) then
      oldsize = size(grid%edge(edge)%solution,dim=1)
   else
      oldsize = 0
   endif
   if (oldsize < edge_size) then
      allocate(temp1(edge_size,grid%system_size,max(1,grid%num_eval)), &
               stat=astat)
      if (astat /= 0) then
         call fatal("allocation failed in unbisect_triangle_pair")
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
            call fatal("allocation failed in unbisect_triangle_pair")
            stop
         endif
         temp1 = 0.0_my_real
         if (oldsize > 0) temp1(1:oldsize,:,:) = grid%edge(edge)%exact
         deallocate(grid%edge(edge)%exact, stat=astat)
         grid%edge(edge)%exact => temp1
      endif
   endif
   grid%edge(edge)%degree = edge_deg
   do j=1,grid%system_size
      if (grid%edge_type(edge,j) == DIRICHLET) then
         call edge_exact(grid,edge,j,"d")
      else
         grid%edge(edge)%solution(edge_deg-1,j,:) = 0
      endif
      if (grid%have_true) call edge_exact(grid,edge,j,"t")
   end do
end do

! set parent and mate element exact

if (par_deg > 2 .and. grid%have_true) then
   do i=1,grid%system_size
      call elem_exact(grid,parent,i,"t")
   end do
endif
if (mate /= BOUNDARY .and. mate_deg > 2 .and. grid%have_true) then
   do i=1,grid%system_size
      call elem_exact(grid,mate,i,"t")
   end do
endif

! don't say I have refined these elements

grid%element(parent)%hrefined_unowned = .false.
grid%element(parent)%isleaf = .true.
if (mate /= BOUNDARY) then
   grid%element(mate)%hrefined_unowned = .false.
   grid%element(mate)%isleaf = .true.
endif

! if the children are leaves in the old solution, the parents become old leaves

if (grid%element(child(1))%oldleaf) then
   grid%element(parent)%oldleaf = .true.
endif
if (mate /= BOUNDARY) then
   if (grid%element(child(3))%oldleaf) then
      grid%element(mate)%oldleaf = .true.
   endif
endif

! set initial solution

call hcoarsen_init_guess(grid,parent,mate,child)

! remove the children elements, edges and vertices from the level linked list,
! putting them on the free space linked list

vert1 = grid%element(child(1))%vertex(3)
if (any(grid%vertex_type(vert1,:) == PERIODIC_MASTER) .or. &
    any(grid%vertex_type(vert1,:) == PERIODIC_SLAVE)) then
   vert2 = grid%element(child(3))%vertex(3)
endif
do i=1,4
   if (i>2 .and. mate == BOUNDARY) exit
   if (associated(grid%element(child(i))%solution)) then
      deallocate(grid%element(child(i))%solution)
   endif
   if (associated(grid%element(child(i))%exact)) then
      deallocate(grid%element(child(i))%exact)
   endif
   if (grid%element(child(i))%previous == END_OF_LIST) then
      grid%head_level_elem(grid%element(child(i))%level) = grid%element(child(i))%next
   else
      grid%element(grid%element(child(i))%previous)%next = grid%element(child(i))%next
   endif
   if (grid%element(child(i))%next == END_OF_LIST) then
      ! nothing, unless I add end_level_elem to the grid data structure
   else
      grid%element(grid%element(child(i))%next)%previous = grid%element(child(i))%previous
   endif
   grid%element(child(i))%next = grid%next_free_elem
   grid%next_free_elem = child(i)
end do
do i=1,6
   if (i>3 .and. mate == BOUNDARY) exit
   if (i==5 .and. childedge(1)==childedge(5)) exit
   if (associated(grid%edge(childedge(i))%solution)) then
      deallocate(grid%edge(childedge(i))%solution)
   endif
   if (associated(grid%edge(childedge(i))%exact)) then
      deallocate(grid%edge(childedge(i))%exact)
   endif
   grid%edge(childedge(i))%next = grid%next_free_edge
   grid%next_free_edge = childedge(i)
end do
if (any(grid%vertex_type(vert1,:) == PERIODIC_MASTER) .or. &
    any(grid%vertex_type(vert1,:) == PERIODIC_SLAVE)) then
   if (grid%vertex(vert2)%previous == END_OF_LIST) then
      grid%head_level_vert(grid%element(child(1))%level) = grid%vertex(vert2)%next
   else
      grid%vertex(grid%vertex(vert2)%previous)%next = grid%vertex(vert2)%next
   endif
   if (grid%vertex(vert2)%next == END_OF_LIST) then
      ! nothing, unless I add end_level_vert to the grid data structure
   else
      grid%vertex(grid%vertex(vert2)%next)%previous = grid%vertex(vert2)%previous
   endif
   grid%vertex(vert2)%next = grid%next_free_vert
   grid%next_free_vert = vert2
endif
if (grid%vertex(vert1)%previous == END_OF_LIST) then
   grid%head_level_vert(grid%element(child(1))%level) = grid%vertex(vert1)%next
else
   grid%vertex(grid%vertex(vert1)%previous)%next = grid%vertex(vert1)%next
endif
if (grid%vertex(vert1)%next == END_OF_LIST) then
   ! nothing, unless I add end_level_vert to the grid data structure
else
   grid%vertex(grid%vertex(vert1)%next)%previous = grid%vertex(vert1)%previous
endif
grid%vertex(vert1)%next = grid%next_free_vert
grid%next_free_vert = vert1

! remove the children elements, edges and vertices from the hash tables

call hash_remove(grid%element(child(1))%gid,grid%elem_hash)
call hash_remove(grid%element(child(2))%gid,grid%elem_hash)
if (mate /= BOUNDARY) then
   call hash_remove(grid%element(child(3))%gid,grid%elem_hash)
   call hash_remove(grid%element(child(4))%gid,grid%elem_hash)
endif
call hash_remove(grid%edge(childedge(1))%gid,grid%edge_hash)
call hash_remove(grid%edge(childedge(2))%gid,grid%edge_hash)
call hash_remove(grid%edge(childedge(3))%gid,grid%edge_hash)
if (mate /= BOUNDARY) then
   call hash_remove(grid%edge(childedge(4))%gid,grid%edge_hash)
   if (childedge(1) /= childedge(5)) then
      call hash_remove(grid%edge(childedge(5))%gid,grid%edge_hash)
      call hash_remove(grid%edge(childedge(6))%gid,grid%edge_hash)
   endif
endif
call hash_remove(grid%vertex(vert1)%gid,grid%vert_hash)
if (any(grid%vertex_type(vert1,:) == PERIODIC_MASTER) .or. &
    any(grid%vertex_type(vert1,:) == PERIODIC_SLAVE)) then
   call hash_remove(grid%vertex(vert2)%gid,grid%vert_hash)
endif

! update scalars

if (mate == BOUNDARY) then
   grid%nelem = grid%nelem - 2
   grid%nelem_leaf = grid%nelem_leaf - 1
   grid%nedge = grid%nedge - 3
! TEMP also need to update nedge_own
else
   grid%nelem = grid%nelem - 4
   grid%nelem_leaf = grid%nelem_leaf - 2
   grid%nedge = grid%nedge - 4
! TEMP also need to update nedge_own
endif
if (any(grid%edge_type(grid%element(parent)%edge(3),:) == PERIODIC_MASTER) .or. &
    any(grid%edge_type(grid%element(parent)%edge(3),:) == PERIODIC_SLAVE)) then
   grid%nedge = grid%nedge - 2
! TEMP also need to update nedge_own
endif
grid%nvert = grid%nvert - 1
if (any(grid%vertex_type(vert1,:) == PERIODIC_MASTER) .or. &
    any(grid%vertex_type(vert1,:) == PERIODIC_SLAVE)) then
   grid%nvert = grid%nvert - 1
endif
do i=1,4
   if (i>2 .and. mate==BOUNDARY) exit
   if (grid%element(child(i))%iown) then
      grid%nelem_leaf_own = grid%nelem_leaf_own - 1
   endif
end do
if (grid%element(parent)%iown) then
   grid%nelem_leaf_own = grid%nelem_leaf_own + 1
endif
if (mate /= BOUNDARY) then
   if (grid%element(mate)%iown) then
      grid%nelem_leaf_own = grid%nelem_leaf_own + 1
   endif
endif
if (grid%element(grid%vertex(vert1)%assoc_elem)%iown) then
   grid%nvert_own = grid%nvert_own - 1
   if (any(grid%vertex_type(vert1,:) == PERIODIC_MASTER) .or. &
       any(grid%vertex_type(vert1,:) == PERIODIC_SLAVE)) then
      grid%nvert_own = grid%nvert_own - 1
   endif
endif
grid%dof = grid%dof - grid%system_size
if (grid%element(grid%vertex(vert1)%assoc_elem)%iown) then
   grid%dof_own = grid%dof_own - grid%system_size
endif
do i=1,4
   if (i>2 .and. mate==BOUNDARY) exit
   if (grid%element(child(i))%degree >= 3) then
      grid%dof = grid%dof - grid%system_size * &
        ((grid%element(child(i))%degree-1)*(grid%element(child(i))%degree-2))/2
      if (grid%element(child(i))%iown) then
       grid%dof_own = grid%dof_own - grid%system_size * &
        ((grid%element(child(i))%degree-1)*(grid%element(child(i))%degree-2))/2
      endif
   endif
end do
if (par_deg >= 3) then
   grid%dof = grid%dof + grid%system_size*((par_deg-1)*(par_deg-2))/2
   if (grid%element(parent)%iown) then
      grid%dof_own = grid%dof_own + grid%system_size*((par_deg-1)*(par_deg-2))/2
   endif
endif
if (mate /= BOUNDARY .and. mate_deg >= 3) then
   grid%dof = grid%dof + grid%system_size*((mate_deg-1)*(mate_deg-2))/2
   if (grid%element(mate)%iown) then
      grid%dof_own = grid%dof_own + grid%system_size*((mate_deg-1)*(mate_deg-2))/2
   endif
endif
do i=1,4
   if (i>3 .and. mate==BOUNDARY) exit
   grid%dof = grid%dof - grid%system_size*(grid%edge(childedge(i))%degree - 1)
   if (grid%element(grid%edge(childedge(i))%assoc_elem)%iown) then
      grid%dof_own = grid%dof_own - grid%system_size * &
                     (grid%edge(childedge(i))%degree - 1)
   endif
end do
do i=1,5
   if (i==4 .and. mate==BOUNDARY) exit
   if (i<=3) then
      edge = grid%element(parent)%edge(i)
   else
      edge = grid%element(mate)%edge(i-3)
   endif
   if (grid%edge(edge)%degree /= old_edge_deg(i)) then
      grid%dof = grid%dof - grid%system_size * &
                 (old_edge_deg(i) - grid%edge(edge)%degree)
      if (grid%element(grid%edge(edge)%assoc_elem)%iown) then
         grid%dof_own = grid%dof_own - grid%system_size * &
                        (old_edge_deg(i) - grid%edge(edge)%degree)
      endif
   endif
end do

! compute new error indicators

call error_indicator(grid,parent,refcont%error_estimator, &
                     grid%element_errind(parent,:), &
                     grid%element(parent)%work)

if (mate /= BOUNDARY) then
   call error_indicator(grid,mate,refcont%error_estimator, &
                        grid%element_errind(mate,:), &
                        grid%element(mate)%work)
endif

end subroutine unbisect_triangle_pair

!          -------------
subroutine p_refine_elem(grid,elem,refine_control,elist,numpref,return_to_elist)
!          -------------

!----------------------------------------------------
! This routine performs p refinement of element elem
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: elem
type(refine_options), intent(in) :: refine_control
type(errind_list), optional, intent(inout) :: elist
integer, optional, intent(inout) :: numpref(:)
logical, optional, intent(in) :: return_to_elist
!----------------------------------------------------
! Local variables:

integer :: i, j, errcode, old_edge_deg(EDGES_PER_ELEMENT)
logical :: add_to_list, need_new_errind
!----------------------------------------------------
! Begin executable code

! if the error indicator list is provided, remove the element from it

if (present(elist)) then
   call remove_from_errind_list(elem,elist)
endif

! refine the interior of the element

call p_refine_element_interior(grid,refine_control,elem,errcode)
if (errcode /= 0) return

! Enforce the edge rule on each edge.
! Also preserve the old edge degrees for p_init_guess.

do i=1,EDGES_PER_ELEMENT
   old_edge_deg(i) = grid%edge(grid%element(elem)%edge(i))%degree
   call enforce_edge_rule(grid,refine_control,grid%element(elem)%edge(i))
end do

! Set initial guess for new solution components.
! Don't set it for REFSOLN because it tries to evaluate a neighbor solution
! which hasn't been allocated yet.

if (refine_control%reftype /= HP_ADAPTIVE .or. &
    (refine_control%hp_strategy /= HP_REFSOLN_EDGE .and. &
     refine_control%hp_strategy /= HP_REFSOLN_ELEM)) then
   if (refine_control%error_estimator == INITIAL_CONDITION) then
      call init_guess_ic(grid,elem)
   else
      call init_guess_p(grid,elem,old_edge_deg)
   endif
endif

! compute new error indicators and
! add to error indicator list if present and control wants it

if (present(elist)) then

   if (present(numpref)) then
      if (numpref(elem) > 0) then
         numpref(elem) = numpref(elem)-1
         add_to_list = numpref(elem) > 0
         grid%element_errind(elem,:) = elist%max_errind*grid%element(elem)%work/binw + .01
         need_new_errind = .false.
      elseif (numpref(elem) < 0) then
         if (present(return_to_elist)) then
            add_to_list = return_to_elist
            need_new_errind = add_to_list
         else
            add_to_list = .false.
            need_new_errind = .false.
         endif
      else
         add_to_list = .false.
         need_new_errind = .false.
      endif
   elseif (present(return_to_elist)) then
      add_to_list = return_to_elist
      need_new_errind = add_to_list
   else
      add_to_list = .false.
      need_new_errind = .false.
   endif
   if (refine_control%hp_strategy == HP_SMOOTH_PRED) then
      need_new_errind=.true.
   endif
   if (need_new_errind) then
      call error_indicator(grid,elem,refine_control%error_estimator, &
                           grid%element_errind(elem,:),grid%element(elem)%work)
   endif
   if (add_to_list) then
      call add_to_errind_list(elem,elist,grid,refine_control)
   endif
endif

end subroutine p_refine_elem

!          -------------------------
subroutine p_refine_element_interior(grid,refine_control,elem,errcode, &
                                     delta_dof,delta_dof_own)
!          -------------------------

!----------------------------------------------------
! This routine performs p refinement of the interior of element elem,
! but not the edges.
!
! Thread safe if delta_dof and delta_dof_own are present.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(refine_options), intent(in) :: refine_control
integer, intent(in) :: elem
integer, intent(out) :: errcode
integer, intent(out), optional :: delta_dof, delta_dof_own
!----------------------------------------------------
! Local variables:

integer :: i, new_deg, new_size, old_size, d1, d2, d3
real(my_real), pointer :: temp(:,:,:)
!----------------------------------------------------
! Begin executable code

if (present(delta_dof) .neqv. present(delta_dof_own)) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("p_refine_element_interior: delta_* must all be present or absent")
   stop
endif

! initialize delta dof's in case we return early

if (present(delta_dof)) then
   delta_dof = 0
   delta_dof_own = 0
endif

! errcode = 0 success
!           1 did not refine

errcode = 1

! verify the element is a leaf

if (.not. grid%element(elem)%isleaf) return

! new degree for the element

new_deg = grid%element(elem)%degree + 1

! if the maximum degree would be exceeded, don't p-refine.

if (new_deg > refine_control%max_deg) then
   return
endif

errcode = 0

! make sure allocated memory is large enough.  If not, reallocate to
! degree+3 so we don't reallocate at every p refinement

if (associated(grid%element(elem)%solution)) then
   old_size = size(grid%element(elem)%solution,dim=1)
else
   old_size = 0
endif
new_size = ((new_deg-1)*(new_deg-2))/2

if (old_size < new_size) then
   new_size = ((new_deg+1)*new_deg)/2
   allocate(temp(new_size,grid%system_size,max(1,grid%num_eval)))
   temp = 0.0_my_real
   if (old_size > 0) then
      temp(1:old_size,:,:) = grid%element(elem)%solution
      deallocate(grid%element(elem)%solution)
   endif
   grid%element(elem)%solution => temp
   nullify(temp)
   if (grid%have_true) then
      allocate(temp(new_size,grid%system_size,max(1,grid%num_eval)))
      temp = 0.0_my_real
      if (old_size > 0) then
         temp(1:old_size,:,:) = grid%element(elem)%exact
         deallocate(grid%element(elem)%exact)
      endif
      grid%element(elem)%exact => temp
   endif
   nullify(temp)
   if (grid%oldsoln_exists) then
      if (associated(grid%element(elem)%oldsoln)) then
         allocate(temp(new_size,grid%system_size,max(1,grid%num_eval)))
         temp = 0.0_my_real
         d1 = min(size(grid%element(elem)%oldsoln,dim=1),size(temp,dim=1))
         d2 = min(size(grid%element(elem)%oldsoln,dim=2),size(temp,dim=2))
         d3 = min(size(grid%element(elem)%oldsoln,dim=3),size(temp,dim=3))
         temp(1:d1,1:d2,1:d3) = grid%element(elem)%oldsoln(1:d1,1:d2,1:d3)
         deallocate(grid%element(elem)%oldsoln)
         grid%element(elem)%oldsoln => temp
         nullify(temp)
      endif
   endif
endif

! Increment the element degree by one.

grid%element(elem)%degree = new_deg

if (present(delta_dof)) then
   if (new_deg >= 3) then
      delta_dof = grid%system_size*(new_deg-2)
      if (grid%element(elem)%iown) then
         delta_dof_own = grid%system_size*(new_deg-2)
      else
         delta_dof_own = 0
      endif
   else
      delta_dof = 0
      delta_dof_own = 0
   endif
else
   if (new_deg >= 3) then
      grid%dof = grid%dof + grid%system_size*(new_deg-2)
      if (grid%element(elem)%iown) grid%dof_own = grid%dof_own + &
                                   grid%system_size*(new_deg-2)
   endif
endif

! set the new solution coefficients to 0, and the exact solution

if (new_deg >= 3) then
   grid%element(elem)%solution(1+((new_deg-2)*(new_deg-3))/2: &
                               ((new_deg-1)*(new_deg-2))/2,:,:)= 0.0_my_real
   if (grid%have_true) then
      do i=1,grid%system_size
         call elem_exact(grid,elem,i,"t")
      end do
   endif
endif

! if I don't own it, set flag for reconciliation

if (.not. grid%element(elem)%iown) then
   grid%element(elem)%prefined_unowned = .true.
endif

! predicted error for SMOOTH_PRED hp strategy
! TEMP only using first eigenvector

grid%element(elem)%sp_eta_pred = &
   sqrt(refine_control%sp_gamma_p)*grid%element_errind(elem,1)

end subroutine p_refine_element_interior

!          -----------------
subroutine enforce_edge_rule(grid,refine_control,edge,delta_dof,delta_dof_own)
!          -----------------

!----------------------------------------------------
! This routine enforces either the minimum or maximum rule for edge
!
! Thread safe if:
!   delta_dof and delta_dof_own are present
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(refine_options), intent(in) :: refine_control
integer, intent(in) :: edge
integer, intent(out), optional :: delta_dof, delta_dof_own
!----------------------------------------------------
! Local variables:

integer :: elem(2), i, old_deg, new_deg, old_size, new_size, copy_size
real(my_real), pointer :: temp(:,:,:)
!----------------------------------------------------
! Begin executable code

if (present(delta_dof) .neqv. present(delta_dof_own)) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("enforce_edge_rule: delta_* must all be present or absent")
   stop
endif

! get the elements that share this edge

elem = get_edge_elements(grid,edge)

! determine what the degree of this element should be

if (elem(2) == BOUNDARY) then
   new_deg = grid%element(elem(1))%degree
elseif (refine_control%edge_rule == MINIMUM_RULE) then
   new_deg = min(grid%element(elem(1))%degree,grid%element(elem(2))%degree)
else ! (refine_control%edge_rule == MAXIMUM_RULE)
   new_deg = max(grid%element(elem(1))%degree,grid%element(elem(2))%degree)
endif

! assign the new degree

old_deg = grid%edge(edge)%degree
if (new_deg /= old_deg) then
   grid%edge(edge)%degree = new_deg
   if (present(delta_dof)) then
      delta_dof = grid%system_size*(new_deg - old_deg)
      if (grid%element(grid%edge(edge)%assoc_elem)%iown) then
        delta_dof_own = grid%system_size*(new_deg - old_deg)
      endif
   else
      grid%dof = grid%dof + grid%system_size*(new_deg - old_deg)
      if (grid%element(grid%edge(edge)%assoc_elem)%iown) then
        grid%dof_own = grid%dof_own + grid%system_size*(new_deg - old_deg)
      endif
   endif
else
   if (present(delta_dof)) then
      delta_dof = 0
      delta_dof_own = 0
   endif
endif

! See if the amount of memory for the solution is correct.
! To save on how often allocation/deallocation occurs,
! if more memory is needed, allocate for degree+2,
! don't deallocate until the new degree is 3 less than the current allocation

if (associated(grid%edge(edge)%solution)) then
   old_size = size(grid%edge(edge)%solution,dim=1)
else
   old_size = 0
endif

if (new_deg-1 < old_size-2) then
   new_size = new_deg-1
   copy_size = new_size
elseif (new_deg-1 > old_size) then
   new_size = new_deg+1
   copy_size = old_size
else
   new_size = old_size
   copy_size = 0
endif

! if necessary, reallocate solution memory

if (new_size /= old_size) then
   allocate(temp(new_size,grid%system_size,max(1,grid%num_eval)))
   temp = 0.0_my_real
   if (old_size > 0) then
      temp(1:copy_size,:,:) = grid%edge(edge)%solution(1:copy_size,:,:)
      deallocate(grid%edge(edge)%solution)
   endif
   grid%edge(edge)%solution => temp
   nullify(temp)
   if (grid%have_true) then
      allocate(temp(new_size,grid%system_size,max(1,grid%num_eval)))
      temp = 0.0_my_real
      if (old_size > 0) then
         temp(1:copy_size,:,:) = grid%edge(edge)%exact(1:copy_size,:,:)
         deallocate(grid%edge(edge)%exact)
      endif
      grid%edge(edge)%exact => temp
      nullify(temp)
   endif
   if (grid%oldsoln_exists) then
      allocate(temp(new_size,grid%system_size,max(1,grid%num_eval)))
      temp = 0.0_my_real
      if (old_size > 0) then
         temp(1:copy_size,:,:) = grid%edge(edge)%oldsoln(1:copy_size,:,:)
         deallocate(grid%edge(edge)%oldsoln)
      endif
      grid%edge(edge)%oldsoln => temp
      nullify(temp)
   endif
endif

! Set Dirichlet boundary conditions on Dirichlet edges, set solution to 0 for
! other new components, and set true solution.

if (new_deg /= old_deg) then
   if (grid%have_true) then
      grid%edge(edge)%exact = 0.0_my_real
   endif
   do i=1,grid%system_size
      if (grid%edge_type(edge,i) == DIRICHLET) then
         call edge_exact(grid,edge,i,"d")
      elseif (new_deg > old_deg) then
         grid%edge(edge)%solution(old_deg:new_size,i,:) = 0.0_my_real
      elseif (new_deg < old_deg) then
         grid%edge(edge)%solution(new_deg:new_size,i,:) = 0.0_my_real
      endif
      if (grid%have_true) call edge_exact(grid,edge,i,"t")
   end do
endif

end subroutine enforce_edge_rule

!          ---------------
subroutine before_h_refine(grid,element_list,nelem,elist,reftype,numpref, &
                           numhref,vert_lid,edge_lid,elem_lid)
!          ---------------

!----------------------------------------------------
! This routine performs parts of h refinement that must be done by a single
! OpenMP thread before the OpenMP-parallel refinement.
! All elements in element_list must be on the same level, and there must only
! be one from a compatibly divisible pair.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: element_list(:), nelem
type(errind_list), optional, intent(inout) :: elist
character(len=*), optional, pointer :: reftype(:)
integer, optional, pointer :: numhref(:), numpref(:)
integer, intent(out) :: vert_lid(:,:), edge_lid(:,:), elem_lid(:,:)
! second dim >= nelem; first dim 2              6              4
!----------------------------------------------------
! Local variables:

integer :: elem, mate, i, errcode, level, hold_vert_head, hold_elem_head, &
           hold_vert_prev, hold_elem_prev
!----------------------------------------------------
! Begin executable code

! if no elements in the list, nothing to do

if (nelem == 0) return

! assign -1 for excess lids, like when the mate is the boundary.
! It will also be returned as -1 if we run out of memory.

vert_lid = -1
edge_lid = -1
elem_lid = -1

level = grid%element(element_list(1))%level + 1

! see if maximum nlev needs to be increased

if (level > size(grid%head_level_elem)) then
   call extend_nlev(grid)
   if (ierr /= NO_ERROR) return
endif

! for each element on the list

do i=1,nelem

! identify the element and mate

   elem = element_list(i)
   if (grid%element(elem)%mate == BOUNDARY) then
      mate = BOUNDARY
   else
      mate = hash_decode_key(grid%element(elem)%mate,grid%elem_hash)
      if (mate == HASH_NOT_FOUND) then
         ierr = PHAML_INTERNAL_ERROR
         call fatal("before_h_refine: mate not found")
         stop
      endif
   endif

! keep track of the current head of the level linked lists in case we abort

   hold_elem_head = grid%head_level_elem(level)
   if (hold_elem_head /= END_OF_LIST) then
      hold_elem_prev = grid%element(hold_elem_head)%previous
   endif
   hold_vert_head = grid%head_level_vert(level)
   if (hold_vert_head /= END_OF_LIST) then
      hold_vert_prev = grid%element(hold_vert_head)%previous
   endif

! get the lids of the children and add them to the level linked lists

! element child 1

   if (grid%next_free_elem == END_OF_LIST) then
      call more_elements(grid,errcode,elist,reftype,numhref,numpref)
      if (errcode /= 0) then
         return
      endif
   endif
   elem_lid(1,i) = grid%next_free_elem
   grid%next_free_elem = grid%element(elem_lid(1,i))%next
   grid%element(elem_lid(1,i))%next = hold_elem_head
   if (hold_elem_head /= END_OF_LIST) then
      grid%element(hold_elem_head)%previous = elem_lid(1,i)
   endif

! element child 2

   if (grid%next_free_elem == END_OF_LIST) then
      call more_elements(grid,errcode,elist,reftype,numhref,numpref)
      if (errcode /= 0) then
         elem_lid = -1
         if (hold_elem_head /= END_OF_LIST) then
            grid%element(hold_elem_head)%previous = hold_elem_prev
         endif
         return
      endif
   endif
   elem_lid(2,i) = grid%next_free_elem
   grid%next_free_elem = grid%element(elem_lid(2,i))%next
   grid%element(elem_lid(2,i))%next = elem_lid(1,i)
   grid%element(elem_lid(1,i))%previous = elem_lid(2,i)
   grid%element(elem_lid(2,i))%previous = END_OF_LIST
   grid%head_level_elem(level) = elem_lid(2,i)

! element child 3

   if (mate /= BOUNDARY) then

      if (grid%next_free_elem == END_OF_LIST) then
         call more_elements(grid,errcode,elist,reftype,numhref,numpref)
         if (errcode /= 0) then
            elem_lid = -1
            if (hold_elem_head /= END_OF_LIST) then
               grid%element(hold_elem_head)%previous = hold_elem_prev
            endif
            grid%head_level_elem(level) = hold_elem_head
            return
         endif
      endif
      elem_lid(3,i) = grid%next_free_elem
      grid%next_free_elem = grid%element(elem_lid(3,i))%next
      grid%element(elem_lid(3,i))%next = elem_lid(2,i)
      grid%element(elem_lid(2,i))%previous = elem_lid(3,i)

! element child 4

      if (grid%next_free_elem == END_OF_LIST) then
         call more_elements(grid,errcode,elist,reftype,numhref,numpref)
         if (errcode /= 0) then
            elem_lid = -1
            if (hold_elem_head /= END_OF_LIST) then
               grid%element(hold_elem_head)%previous = hold_elem_prev
            endif
            grid%head_level_elem(level) = hold_elem_head
            return
         endif
      endif
      elem_lid(4,i) = grid%next_free_elem
      grid%next_free_elem = grid%element(elem_lid(4,i))%next
      grid%element(elem_lid(4,i))%next = elem_lid(3,i)
      grid%element(elem_lid(3,i))%previous = elem_lid(4,i)
      grid%element(elem_lid(4,i))%previous = END_OF_LIST
      grid%head_level_elem(level) = elem_lid(4,i)

   endif

! edge child 1

   if (grid%next_free_edge == END_OF_LIST) then
      call more_edges(grid,errcode)
      if (errcode /= 0) then
         elem_lid = -1
         if (hold_elem_head /= END_OF_LIST) then
            grid%element(hold_elem_head)%previous = hold_elem_prev
         endif
         grid%head_level_elem(level) = hold_elem_head
         return
      endif
   endif
   edge_lid(1,i) = grid%next_free_edge
   grid%next_free_edge = grid%edge(edge_lid(1,i))%next

! edge child 2

   if (grid%next_free_edge == END_OF_LIST) then
      call more_edges(grid,errcode)
      if (errcode /= 0) then
         elem_lid = -1
         edge_lid = -1
         if (hold_elem_head /= END_OF_LIST) then
            grid%element(hold_elem_head)%previous = hold_elem_prev
         endif
         grid%head_level_elem(level) = hold_elem_head
         return
      endif
   endif
   edge_lid(2,i) = grid%next_free_edge
   grid%next_free_edge = grid%edge(edge_lid(2,i))%next

! edge child 3

   if (grid%next_free_edge == END_OF_LIST) then
      call more_edges(grid,errcode)
      if (errcode /= 0) then
         elem_lid = -1
         edge_lid = -1
         if (hold_elem_head /= END_OF_LIST) then
            grid%element(hold_elem_head)%previous = hold_elem_prev
         endif
         grid%head_level_elem(level) = hold_elem_head
         return
      endif
   endif
   edge_lid(3,i) = grid%next_free_edge
   grid%next_free_edge = grid%edge(edge_lid(3,i))%next

! edge child 4

   if (mate /= BOUNDARY) then
      if (grid%next_free_edge == END_OF_LIST) then
         call more_edges(grid,errcode)
         if (errcode /= 0) then
            elem_lid = -1
            edge_lid = -1
            if (hold_elem_head /= END_OF_LIST) then
               grid%element(hold_elem_head)%previous = hold_elem_prev
            endif
            grid%head_level_elem(level) = hold_elem_head
            return
         endif
      endif
      edge_lid(4,i) = grid%next_free_edge
      grid%next_free_edge = grid%edge(edge_lid(4,i))%next
   endif

! edge child 5

   if (any(grid%edge_type(grid%element(elem)%edge(3),:)==PERIODIC_MASTER) .or. &
       any(grid%edge_type(grid%element(elem)%edge(3),:)==PERIODIC_SLAVE)) then

      if (grid%next_free_edge == END_OF_LIST) then
         call more_edges(grid,errcode)
         if (errcode /= 0) then
            elem_lid = -1
            edge_lid = -1
            if (hold_elem_head /= END_OF_LIST) then
               grid%element(hold_elem_head)%previous = hold_elem_prev
            endif
            grid%head_level_elem(level) = hold_elem_head
            return
         endif
      endif
      edge_lid(5,i) = grid%next_free_edge
      grid%next_free_edge = grid%edge(edge_lid(5,i))%next
      if (any(grid%edge_type(grid%element(elem)%edge(3),:)==PERIODIC_MASTER)) then
         grid%edge(edge_lid(5,i))%next = edge_lid(1,i)
      else
         grid%edge(edge_lid(1,i))%next = edge_lid(5,i)
      endif

! edge child 6

      if (grid%next_free_edge == END_OF_LIST) then
         call more_edges(grid,errcode)
         if (errcode /= 0) then
            elem_lid = -1
            edge_lid = -1
            if (hold_elem_head /= END_OF_LIST) then
               grid%element(hold_elem_head)%previous = hold_elem_prev
            endif
            grid%head_level_elem(level) = hold_elem_head
            return
         endif
      endif
      edge_lid(6,i) = grid%next_free_edge
      grid%next_free_edge = grid%edge(edge_lid(6,i))%next
      if (any(grid%edge_type(grid%element(elem)%edge(3),:)==PERIODIC_MASTER)) then
         grid%edge(edge_lid(6,i))%next = edge_lid(2,i)
      else
         grid%edge(edge_lid(2,i))%next = edge_lid(6,i)
      endif

   endif

! vertex child 1

   if (grid%next_free_vert == END_OF_LIST) then
      call more_verts(grid,errcode)
      if (errcode /= 0) then
         elem_lid = -1
         edge_lid = -1
         if (hold_elem_head /= END_OF_LIST) then
            grid%element(hold_elem_head)%previous = hold_elem_prev
         endif
         grid%head_level_elem(level) = hold_elem_head
         return
      endif
   endif
   vert_lid(1,i) = grid%next_free_vert
   grid%next_free_vert = grid%vertex(vert_lid(1,i))%next
   grid%vertex(vert_lid(1,i))%next = hold_vert_head
   grid%head_level_vert(level) = vert_lid(1,i)
   if (hold_vert_head /= END_OF_LIST) then
      grid%vertex(hold_vert_head)%previous = vert_lid(1,i)
   endif
   grid%vertex(vert_lid(1,i))%previous = END_OF_LIST

! vertex child 2

   if (any(grid%edge_type(grid%element(elem)%edge(3),:)==PERIODIC_MASTER) .or. &
       any(grid%edge_type(grid%element(elem)%edge(3),:)==PERIODIC_SLAVE)) then

      if (grid%next_free_vert == END_OF_LIST) then
         call more_verts(grid,errcode)
         if (errcode /= 0) then
            elem_lid = -1
            edge_lid = -1
            vert_lid = -1
            if (hold_elem_head /= END_OF_LIST) then
               grid%element(hold_elem_head)%previous = hold_elem_prev
               grid%head_level_elem(level) = hold_elem_head
            endif
            if (hold_vert_head /= END_OF_LIST) then
               grid%vertex(hold_vert_head)%previous = hold_vert_prev
            endif
            grid%head_level_vert(level) = hold_vert_head
            return
         endif
      endif
      vert_lid(2,i) = grid%next_free_vert
      grid%next_free_vert = grid%vertex(vert_lid(2,i))%next
      if (any(grid%edge_type(grid%element(elem)%edge(3),:)==PERIODIC_MASTER)) then
         grid%vertex(vert_lid(2,i))%next = vert_lid(1,i)
         grid%head_level_vert(level) = vert_lid(2,i)
         grid%vertex(vert_lid(1,i))%previous = vert_lid(2,i)
         grid%vertex(vert_lid(2,i))%previous = END_OF_LIST
      else
         grid%vertex(vert_lid(2,i))%next = grid%vertex(vert_lid(1,i))%next
         grid%vertex(vert_lid(1,i))%next = vert_lid(2,i)
         grid%vertex(vert_lid(2,i))%previous = vert_lid(1,i)
         if (grid%vertex(vert_lid(2,i))%next /= END_OF_LIST) then
            grid%vertex(hold_vert_head)%previous = vert_lid(2,i)
         endif
      endif

   endif

end do ! next element

end subroutine before_h_refine

!                    --------------------
recursive subroutine bisect_triangle_pair(grid,parent,errcode,refcont,elist, &
                                          reftype,numhref,numpref, &
                     return_to_elist,reduce_p,reduce_p_max,vert_child_lid, &
                     edge_child_lid,elem_child_lid,delta_dof,delta_dof_own, &
                     delta_nelem,delta_nelem_leaf,delta_nelem_leaf_own, &
                     delta_nedge,delta_nedge_own,delta_nvert,delta_nvert_own, &
                     max_nlev)
!                    --------------------

!----------------------------------------------------
! This routine refines triangle parent and its mate by bisection.
! errcode = 0 for success, 1 if the grid is full
!
! thread safe if:
!   all elements being refined in parallel have the same h-level and are
!     compatibly divisible
!   an element and its mate are not both in the list of elements being refined
!   elist is not present
!   vert_child_lid, edge_child_lid, elem_child_lid are present
!     (in this case, the tasks performed in before_h_refine and after_h_refine
!     are not performed here)
!   delta_* and max_nlev are present
!     (in this case, the grid scalars are not updated)
!   the solution at preexisting vertices does not change
!   if numhref is present, it is at most 1 for parent and mate
!   if numpref is present, it is 0 for parent and mate
!   if return_to_elist is present, it is false
!   TEMP reduce_p is not present
!   TEMP no periodic boundary conditions
!   some other conditions that I think hold; search for OPENMP
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: parent
integer, intent(inout) :: errcode
type(refine_options), intent(in) :: refcont
type(errind_list), optional, intent(inout) :: elist
character(len=*), optional, pointer :: reftype(:)
integer, optional, pointer :: numhref(:), numpref(:)
logical, optional, intent(in) :: return_to_elist
integer, optional, intent(in) :: reduce_p
integer, optional, intent(inout) :: reduce_p_max
integer, optional, intent(in) :: vert_child_lid(:), edge_child_lid(:), &
                                 elem_child_lid(:)
integer, optional, intent(out) :: delta_dof, delta_dof_own, delta_nelem, &
                                  delta_nelem_leaf, delta_nelem_leaf_own, &
                                  delta_nedge, delta_nedge_own, delta_nvert, &
                                  delta_nvert_own, max_nlev


!----------------------------------------------------
! Local variables:

integer :: child1, child2, child3, child4
integer :: edge1, edge2, edge3, edge4, edge5, edge6, vert1, vert2
integer :: mate, modval, i, j, astat, paredge, mateedge, masterparent, newdeg
type(hash_key) :: grandparent, grandparentmate, stepgrandparentmate, tempgid
integer :: bctype(grid%system_size)
real(my_real) :: bcrhs(grid%system_size), &
                 bccoef(grid%system_size,grid%system_size)
logical :: add_to_list, need_new_errind
!----------------------------------------------------
! Begin executable code

errcode = 0

! if the delta_* arguments are present, the change in the grid scalar is
! returned; otherwise the grid scalars are updated

if ((present(delta_dof) .neqv. present(delta_dof_own)) .or. &
    (present(delta_dof) .neqv. present(delta_nelem)) .or. &
    (present(delta_dof) .neqv. present(delta_nelem_leaf)) .or. &
    (present(delta_dof) .neqv. present(delta_nelem_leaf_own)) .or. &
    (present(delta_dof) .neqv. present(delta_nedge)) .or. &
    (present(delta_dof) .neqv. present(delta_nedge_own)) .or. &
    (present(delta_dof) .neqv. present(delta_nvert)) .or. &
    (present(delta_dof) .neqv. present(delta_nvert_own)) .or. &
    (present(delta_dof) .neqv. present(max_nlev))) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("bisect_triangle_pair: all or none of the deltas must be present")
   stop
endif

if (present(delta_dof)) then
   delta_dof = 0
   delta_dof_own = 0
   delta_nelem = 0
   delta_nelem_leaf = 0
   delta_nelem_leaf_own = 0
   delta_nedge = 0
   delta_nedge_own = 0
   delta_nvert = 0
   delta_nvert_own = 0
   max_nlev = 0
endif

! if the first elem_child_lid is -1, we ran out of memory

if (present(elem_child_lid)) then
   if (elem_child_lid(1) == -1) return
endif

! if the error indicator list is provided, remove the parent from it.
! OPENMP not thread safe; requires elist not present

if (present(elist)) then
   call remove_from_errind_list(parent,elist)
endif

! if the maximum number of levels would be exceeded or the hash key would
! overflow, don't h-refine.

if (grid%element(parent)%level >= refcont%max_lev .or. &
    hash_overflow(grid%edge(grid%element(parent)%edge(3))%gid,4,3) .or. &
    hash_overflow(grid%element(parent)%gid,2,1)) then
   return
endif

! see if refinement of this element requires extending the number of levels

if (.not. present(elem_child_lid)) then
   if (grid%element(parent)%level >= size(grid%head_level_elem)) then
      call extend_nlev(grid)
      if (ierr /= NO_ERROR) return
   endif
endif

! reduce_p and reduce_p_max should both be present or absent

if (present(reduce_p) .neqv. present(reduce_p_max)) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("reduce_p and reduce_p_max must both be present or absent in bisect_triangle_pair")
   stop
endif

! identify the mate, creating it if necessary.
! OPENMP not thread safe; requires compatibly divisible elements

if (present(reduce_p)) reduce_p_max = 0
if (grid%element(parent)%mate == BOUNDARY) then
   mate = BOUNDARY
else
   mate = hash_decode_key(grid%element(parent)%mate,grid%elem_hash)
   if (mate == HASH_NOT_FOUND) then
      mate = hash_decode_key(grid%element(parent)%mate/2,grid%elem_hash)
      call bisect_triangle_pair(grid,mate,errcode,refcont,elist,reftype, &
                                numhref,numpref,return_to_elist,reduce_p, &
                                reduce_p_max)
      if (errcode /= 0) return
      mate = hash_decode_key(grid%element(parent)%mate,grid%elem_hash)
      if (present(reftype)) reftype(mate) = "h"
   endif
   if (present(elist)) call remove_from_errind_list(mate,elist)
endif
if (present(reduce_p)) then
   if (reduce_p_max == 0) reduce_p_max = reduce_p
endif

! if the hash key would overflow, don't h-refine.

if (mate /= BOUNDARY) then
   if (hash_overflow(grid%element(mate)%gid,2,1)) then
      return
   endif
endif

! get the indices for the new elements, edges and vertex
! OPENMP thread safe iff elem_child_lid and friends are present

if ((present(elem_child_lid) .neqv. present(edge_child_lid)) .or. &
    (present(elem_child_lid) .neqv. present(vert_child_lid))) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("bisect_triangle_pair: all child_lid arrays must be present or absent")
   stop
endif

! indices were passed in

if (present(elem_child_lid)) then
   child1 = elem_child_lid(1)
   child2 = elem_child_lid(2)
   child3 = elem_child_lid(3)
   child4 = elem_child_lid(4)
   edge1 = edge_child_lid(1)
   edge2 = edge_child_lid(2)
   edge3 = edge_child_lid(3)
   edge4 = edge_child_lid(4)
   edge5 = edge_child_lid(5)
   edge6 = edge_child_lid(6)
   vert1 = vert_child_lid(1)
   vert2 = vert_child_lid(2)

! determine the indices here

else

   if (grid%next_free_elem == END_OF_LIST) then
      call more_elements(grid,errcode,elist,reftype,numhref,numpref)
      if (errcode /= 0) return
   endif
   child1 = grid%next_free_elem
   grid%next_free_elem = grid%element(child1)%next
   if (grid%next_free_elem == END_OF_LIST) then
      call more_elements(grid,errcode,elist,reftype,numhref,numpref)
      if (errcode /= 0) return
   endif
   child2 = grid%next_free_elem
   grid%next_free_elem = grid%element(child2)%next
   if (mate /= BOUNDARY) then
      if (grid%next_free_elem == END_OF_LIST) then
         call more_elements(grid,errcode,elist,reftype,numhref,numpref)
         if (errcode /= 0) return
      endif
      child3 = grid%next_free_elem
      grid%next_free_elem = grid%element(child3)%next
      if (grid%next_free_elem == END_OF_LIST) then
         call more_elements(grid,errcode,elist,reftype,numhref,numpref)
         if (errcode /= 0) return
      endif
      child4 = grid%next_free_elem
      grid%next_free_elem = grid%element(child4)%next
   endif
   if (grid%next_free_edge == END_OF_LIST) then
      call more_edges(grid,errcode)
      if (errcode /= 0) return
   endif
   edge1 = grid%next_free_edge
   grid%next_free_edge = grid%edge(edge1)%next
   if (grid%next_free_edge == END_OF_LIST) then
      call more_edges(grid,errcode)
      if (errcode /= 0) return
   endif
   edge2 = grid%next_free_edge
   grid%next_free_edge = grid%edge(edge2)%next
   if (grid%next_free_edge == END_OF_LIST) then
      call more_edges(grid,errcode)
      if (errcode /= 0) return
   endif
   edge3 = grid%next_free_edge
   grid%next_free_edge = grid%edge(edge3)%next
   if (mate /= BOUNDARY) then
      if (grid%next_free_edge == END_OF_LIST) then
         call more_edges(grid,errcode)
         if (errcode /= 0) return
      endif
      edge4 = grid%next_free_edge
      grid%next_free_edge = grid%edge(edge4)%next
   endif
   if (grid%next_free_vert == END_OF_LIST) then
      call more_verts(grid,errcode)
      if (errcode /= 0) return
   endif
   vert1 = grid%next_free_vert
   grid%next_free_vert = grid%vertex(vert1)%next

endif ! elem_child_lid present

if (grid%element(parent)%level /= 1) then
   grandparent = grid%element(parent)%gid/2
   grandparentmate = grid%element(hash_decode_key(grandparent,grid%elem_hash))%mate
endif

grid%element(parent)%isleaf = .false.

modval = mod(grid%element(parent)%gid,2)

! first child

grid%element(child1)%gid = 2*grid%element(parent)%gid
grid%element(child1)%edge = (/ edge3, edge1, &
                               grid%element(parent)%edge(2) /)
grid%element(child1)%vertex = (/ grid%element(parent)%vertex(1), &
                                 grid%element(parent)%vertex(3), &
                                 vert1 /)
if (grid%element(parent)%level == 1) then
   grid%element(child1)%mate = level2_mate(grid,parent,2)
elseif (grandparentmate == BOUNDARY) then
   grid%element(child1)%mate = BOUNDARY
else
   if (modval == 0) then
      grid%element(child1)%mate = 4*grandparentmate
   else
      grid%element(child1)%mate = 4*grandparentmate+2
   endif
endif
grid%element(child1)%level = grid%element(parent)%level+1
grid%element(child1)%iown = grid%element(parent)%iown
grid%element(child1)%hrefined_unowned = .false.
grid%element(child1)%prefined_unowned = .false.
if (.not. present(elem_child_lid)) then
   grid%element(child1)%next = grid%head_level_elem(grid%element(child1)%level)
   if (grid%element(child1)%next /= END_OF_LIST) then
      grid%element(grid%element(child1)%next)%previous = child1
   endif
endif
grid%element(child1)%isleaf = .true.
grid%element(child1)%oldleaf = .false.
grid%element(child1)%degree = grid%element(parent)%degree
! TEMP only using first eigenvector
grid%element(child1)%sp_eta_pred = refcont%sp_gamma_h*grid%element_errind(parent,1)/ &
   sqrt(2.0_my_real)**(grid%element(parent)%degree+2)
if (associated(grid%element(parent)%solution)) then
   allocate(grid%element(child1)%solution(size(grid%element(parent)%solution, &
            dim=1),grid%system_size,max(1,grid%num_eval)),stat=astat)
   grid%element(child1)%solution = 0.0_my_real
   if (grid%have_true) then
      allocate(grid%element(child1)%exact(size(grid%element(parent)%exact, &
               dim=1),grid%system_size,max(1,grid%num_eval)),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("memory allocation failed in bisect_triangle_pair")
         return
      endif
      grid%element(child1)%exact = 0.0_my_real
   endif
endif
if (grid%element(child1)%degree >= 3) then
   if (present(delta_dof)) then
      delta_dof = delta_dof + grid%system_size * &
           ((grid%element(child1)%degree-2)*(grid%element(child1)%degree-1))/2
   else
      grid%dof = grid%dof + grid%system_size * &
           ((grid%element(child1)%degree-2)*(grid%element(child1)%degree-1))/2
   endif
endif
if (grid%element(child1)%iown .and. grid%element(child1)%degree >= 3) then
   if (present(delta_dof)) then
      delta_dof_own = delta_dof_own + grid%system_size * &
           ((grid%element(child1)%degree-2)*(grid%element(child1)%degree-1))/2
   else
      grid%dof_own = grid%dof_own + grid%system_size * &
           ((grid%element(child1)%degree-2)*(grid%element(child1)%degree-1))/2
   endif
endif

if (.not. present(elem_child_lid)) then
   call hash_insert(grid%element(child1)%gid,child1,grid%elem_hash)
endif

! second child

grid%element(child2)%gid = grid%element(child1)%gid+1
grid%element(child2)%edge = (/ edge3, edge2, &
                               grid%element(parent)%edge(1) /)
grid%element(child2)%vertex = (/ grid%element(parent)%vertex(2), &
                                 grid%element(parent)%vertex(3), &
                                 vert1 /)
if (grid%element(parent)%level == 1) then
   grid%element(child2)%mate = level2_mate(grid,parent,1)
else
   if (modval == 0) then
      grid%element(child2)%mate = grid%element(child2)%gid+2
   else
      grid%element(child2)%mate = grid%element(child2)%gid-2
   endif
endif
grid%element(child2)%level = grid%element(parent)%level+1
grid%element(child2)%iown = grid%element(parent)%iown
if (.not. present(elem_child_lid)) then
   grid%element(child2)%next = child1
   grid%element(child1)%previous = child2
   grid%element(child2)%previous = END_OF_LIST
endif
grid%element(child2)%hrefined_unowned = .false.
grid%element(child2)%prefined_unowned = .false.
grid%element(child2)%isleaf = .true.
grid%element(child2)%oldleaf = .false.
grid%element(child2)%degree = grid%element(parent)%degree
grid%element(child2)%sp_eta_pred = grid%element(child1)%sp_eta_pred
if (associated(grid%element(parent)%solution)) then
   allocate(grid%element(child2)%solution(size(grid%element(parent)%solution, &
            dim=1),grid%system_size,max(1,grid%num_eval)),stat=astat)
   grid%element(child2)%solution = 0.0_my_real
   if (grid%have_true) then
      allocate(grid%element(child2)%exact(size(grid%element(parent)%exact, &
               dim=1),grid%system_size,max(1,grid%num_eval)),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("memory allocation failed in bisect_triangle_pair")
         return
      endif
      grid%element(child2)%exact = 0.0_my_real
   endif
endif

if (.not. present(elem_child_lid)) then
   call hash_insert(grid%element(child2)%gid,child2,grid%elem_hash)
endif

if (.not. grid%element(parent)%iown) then
   grid%element(parent)%hrefined_unowned = .true.
endif

if (mate == BOUNDARY) then

if (.not. present(elem_child_lid)) then
   grid%head_level_elem(grid%element(child2)%level) = child2
endif

else

if (grid%element(mate)%level /= 1) stepgrandparentmate = &
  grid%element(hash_decode_key(grid%element(parent)%mate/2,grid%elem_hash))%mate

modval = mod(grid%element(mate)%gid,2)

grid%element(mate)%isleaf = .false.

! third child

grid%element(child3)%gid = 2*grid%element(mate)%gid
grid%element(child3)%edge = (/ edge4, edge1, &
                               grid%element(mate)%edge(2) /)
grid%element(child3)%vertex = (/ grid%element(mate)%vertex(1), &
                                 grid%element(mate)%vertex(3), &
                                 vert1 /)
if (grid%element(mate)%level == 1) then
   grid%element(child3)%mate = level2_mate(grid,mate,2)
elseif (stepgrandparentmate == BOUNDARY) then
   grid%element(child3)%mate = BOUNDARY
else
   if (modval == 0) then
      grid%element(child3)%mate = 4*stepgrandparentmate
   else
      grid%element(child3)%mate = 4*stepgrandparentmate+2
   endif
endif
grid%element(child3)%level = grid%element(mate)%level+1
grid%element(child3)%iown = grid%element(mate)%iown
grid%element(child3)%hrefined_unowned = .false.
grid%element(child3)%prefined_unowned = .false.
if (.not. present(elem_child_lid)) then
   grid%element(child3)%next = child2
   grid%element(child2)%previous = child3
endif
grid%element(child3)%isleaf = .true.
grid%element(child3)%oldleaf = .false.
grid%element(child3)%degree = grid%element(mate)%degree
! TEMP only using first eigenvector
grid%element(child3)%sp_eta_pred = refcont%sp_gamma_h*grid%element_errind(mate,1)/ &
   sqrt(2.0_my_real)**(grid%element(mate)%degree+2)
if (associated(grid%element(mate)%solution)) then
   allocate(grid%element(child3)%solution(size(grid%element(mate)%solution, &
            dim=1),grid%system_size,max(1,grid%num_eval)),stat=astat)
   grid%element(child3)%solution = 0.0_my_real
   if (grid%have_true) then
      allocate(grid%element(child3)%exact(size(grid%element(mate)%exact, &
               dim=1),grid%system_size,max(1,grid%num_eval)),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("memory allocation failed in bisect_triangle_pair")
         return
      endif
      grid%element(child3)%exact = 0.0_my_real
   endif
endif
if (grid%element(child3)%degree >= 3) then
   if (present(delta_dof)) then
      delta_dof = delta_dof + grid%system_size * &
           ((grid%element(child3)%degree-2)*(grid%element(child3)%degree-1))/2
   else
      grid%dof = grid%dof + grid%system_size * &
           ((grid%element(child3)%degree-2)*(grid%element(child3)%degree-1))/2
   endif
endif
if (grid%element(child3)%iown .and. grid%element(child3)%degree >= 3) then
   if (present(delta_dof)) then
      delta_dof_own = delta_dof_own + grid%system_size * &
           ((grid%element(child3)%degree-2)*(grid%element(child3)%degree-1))/2
   else
      grid%dof_own = grid%dof_own + grid%system_size * &
           ((grid%element(child3)%degree-2)*(grid%element(child3)%degree-1))/2
   endif
endif


if (.not. present(elem_child_lid)) then
   call hash_insert(grid%element(child3)%gid,child3,grid%elem_hash)
endif

! fourth child

grid%element(child4)%gid = grid%element(child3)%gid+1
grid%element(child4)%edge = (/ edge4, edge2, &
                               grid%element(mate)%edge(1) /)
grid%element(child4)%vertex = (/ grid%element(mate)%vertex(2), &
                                 grid%element(mate)%vertex(3), &
                                 vert1 /)
if (grid%element(mate)%level == 1) then
   grid%element(child4)%mate = level2_mate(grid,mate,1)
else
   if (modval == 0) then
      grid%element(child4)%mate = grid%element(child4)%gid+2
   else
      grid%element(child4)%mate = grid%element(child4)%gid-2
   endif
endif
grid%element(child4)%level = grid%element(mate)%level+1
grid%element(child4)%iown = grid%element(mate)%iown
grid%element(child4)%hrefined_unowned = .false.
grid%element(child4)%prefined_unowned = .false.
if (.not. present(elem_child_lid)) then
   grid%element(child4)%next = child3
   grid%element(child3)%previous = child4
   grid%element(child4)%previous = END_OF_LIST
endif
grid%element(child4)%isleaf = .true.
grid%element(child4)%oldleaf = .false.
grid%element(child4)%degree = grid%element(mate)%degree
grid%element(child4)%sp_eta_pred = grid%element(child3)%sp_eta_pred
if (associated(grid%element(mate)%solution)) then
   allocate(grid%element(child4)%solution(size(grid%element(mate)%solution, &
            dim=1),grid%system_size,max(1,grid%num_eval)),stat=astat)
   grid%element(child4)%solution = 0.0_my_real
   if (grid%have_true) then
      allocate(grid%element(child4)%exact(size(grid%element(mate)%exact, &
               dim=1),grid%system_size,max(1,grid%num_eval)),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("memory allocation failed in bisect_triangle_pair")
         return
      endif
      grid%element(child4)%exact = 0.0_my_real
   endif
endif

if (.not. present(elem_child_lid)) then
   call hash_insert(grid%element(child4)%gid,child4,grid%elem_hash)
   grid%head_level_elem(grid%element(child4)%level) = child4
endif

if (.not. grid%element(mate)%iown) then
   grid%element(mate)%hrefined_unowned = .true.
endif

endif ! mate == BOUNDARY

! new vertex

if (mate == BOUNDARY) then
   grid%vertex(vert1)%gid = grid%element(child1)%gid
   grid%vertex(vert1)%assoc_elem = child1
elseif (grid%element(child1)%gid < grid%element(child3)%gid) then
   grid%vertex(vert1)%gid = grid%element(child1)%gid
   grid%vertex(vert1)%assoc_elem = child1
else
   grid%vertex(vert1)%gid = grid%element(child3)%gid
   grid%vertex(vert1)%assoc_elem = child3
endif
call point_on_edge(grid, &
                   grid%vertex(grid%element(parent)%vertex(1))%coord%x, &
                   grid%vertex(grid%element(parent)%vertex(1))%coord%y, &
                   grid%vertex(grid%element(parent)%vertex(1))%bmark, &
                   grid%vertex(grid%element(parent)%vertex(1))%bparam, &
                   grid%vertex(grid%element(parent)%vertex(2))%coord%x, &
                   grid%vertex(grid%element(parent)%vertex(2))%coord%y, &
                   grid%vertex(grid%element(parent)%vertex(2))%bmark, &
                   grid%vertex(grid%element(parent)%vertex(2))%bparam, &
                   grid%vertex(grid%element(parent)%vertex(3))%coord%x, &
                   grid%vertex(grid%element(parent)%vertex(3))%coord%y, &
                   grid%edge(grid%element(parent)%edge(3))%bmark, &
                   0.5_my_real, &
                   grid%vertex(vert1)%coord%x,grid%vertex(vert1)%coord%y, &
                   grid%vertex(vert1)%bparam)
grid%vertex(vert1)%bmark = grid%edge(grid%element(parent)%edge(3))%bmark
grid%vertex_type(vert1,:) = grid%edge_type(grid%element(parent)%edge(3),:)
if (any(grid%vertex_type(vert1,:) == DIRICHLET)) then
   call bconds(grid%vertex(vert1)%coord%x,grid%vertex(vert1)%coord%y, &
               grid%vertex(vert1)%bmark,bctype,bccoef,bcrhs)
endif
do i=1,grid%system_size
   if (grid%vertex_type(vert1,i) == DIRICHLET) then
      grid%vertex_solution(vert1,i,:) = bcrhs(i)
   else
      grid%vertex_solution(vert1,i,:) = 0.0_my_real
   endif
end do
if (.not. present(elem_child_lid)) then
   grid%vertex(vert1)%next = grid%head_level_vert(grid%element(child1)%level)
   grid%head_level_vert(grid%element(child1)%level) = vert1
   if (grid%vertex(vert1)%next /= END_OF_LIST) then
      grid%vertex(grid%vertex(vert1)%next)%previous = vert1
   endif
   grid%vertex(vert1)%previous = END_OF_LIST
endif
if (grid%have_true) then
   do j=1,max(1,grid%num_eval)
      do i=1,grid%system_size
         grid%vertex_exact(vert1,i,j)=trues(grid%vertex(vert1)%coord%x, &
                                            grid%vertex(vert1)%coord%y, &
                                            i,j)
      end do
   end do
endif
if (present(delta_dof)) then
   delta_dof = delta_dof + grid%system_size
   if (grid%element(grid%vertex(vert1)%assoc_elem)%iown) &
      delta_dof_own = delta_dof_own + grid%system_size
else
   grid%dof = grid%dof + grid%system_size
   if (grid%element(grid%vertex(vert1)%assoc_elem)%iown) &
      grid%dof_own = grid%dof_own + grid%system_size
endif

if (.not. present(elem_child_lid)) then
   call hash_insert(grid%vertex(vert1)%gid,vert1,grid%vert_hash)
endif

! second new vertex, if periodic boundary point(s)

masterparent = -1

if (any(grid%vertex_type(vert1,:) == PERIODIC_MASTER) .or. &
    any(grid%vertex_type(vert1,:) == PERIODIC_SLAVE)) then

   if (.not. present(elem_child_lid)) then
      if (grid%next_free_vert == END_OF_LIST) then
         call more_verts(grid,errcode)
         if (errcode /= 0) return
      endif
      vert2 = grid%next_free_vert
      grid%next_free_vert = grid%vertex(vert2)%next
   endif

   if (any(grid%vertex_type(vert1,:) == PERIODIC_MASTER)) then
      masterparent = parent
   else
      masterparent = mate
   endif

   grid%element(child3)%vertex(3) = vert2
   grid%element(child4)%vertex(3) = vert2

   if (grid%element(child2)%gid < grid%element(child4)%gid) then
      grid%vertex(vert2)%gid = grid%element(child2)%gid
   else
      grid%vertex(vert2)%gid = grid%element(child4)%gid
   endif
   if (masterparent == mate) then
      tempgid = grid%vertex(vert1)%gid
      grid%vertex(vert1)%gid = grid%vertex(vert2)%gid
      grid%vertex(vert2)%gid = tempgid
      if (.not. present(elem_child_lid)) then
         call hash_remove(grid%vertex(vert2)%gid,grid%vert_hash)
         call hash_insert(grid%vertex(vert1)%gid,vert1,grid%vert_hash)
      endif
      if (grid%vertex(vert1)%gid == grid%element(child2)%gid) then
         grid%vertex(vert1)%assoc_elem = child2
      elseif (grid%vertex(vert1)%gid == grid%element(child4)%gid) then
         grid%vertex(vert1)%assoc_elem = child4
      endif
   endif
   if (.not. present(elem_child_lid)) then
      call hash_insert(grid%vertex(vert2)%gid,vert2,grid%vert_hash)
   endif
   if (masterparent == parent) then
      grid%vertex(vert1)%assoc_elem = child1
   else
      grid%vertex(vert1)%assoc_elem = child3
   endif
   grid%vertex(vert2)%assoc_elem = grid%vertex(vert1)%assoc_elem
   call point_on_edge(grid, &
                      grid%vertex(grid%element(mate)%vertex(1))%coord%x, &
                      grid%vertex(grid%element(mate)%vertex(1))%coord%y, &
                      grid%vertex(grid%element(mate)%vertex(1))%bmark, &
                      grid%vertex(grid%element(mate)%vertex(1))%bparam, &
                      grid%vertex(grid%element(mate)%vertex(2))%coord%x, &
                      grid%vertex(grid%element(mate)%vertex(2))%coord%y, &
                      grid%vertex(grid%element(mate)%vertex(2))%bmark, &
                      grid%vertex(grid%element(mate)%vertex(2))%bparam, &
                      grid%vertex(grid%element(mate)%vertex(3))%coord%x, &
                      grid%vertex(grid%element(mate)%vertex(3))%coord%y, &
                      grid%edge(grid%element(mate)%edge(3))%bmark, &
                      0.5_my_real, &
                      grid%vertex(vert2)%coord%x,grid%vertex(vert2)%coord%y, &
                      grid%vertex(vert2)%bparam)
   grid%vertex(vert2)%bmark = grid%edge(grid%element(mate)%edge(3))%bmark
   grid%vertex_type(vert2,:) = grid%edge_type(grid%element(mate)%edge(3),:)
   if (any(grid%vertex_type(vert2,:) == DIRICHLET)) then
      call bconds(grid%vertex(vert2)%coord%x,grid%vertex(vert2)%coord%y, &
                  grid%vertex(vert2)%bmark,bctype,bccoef,bcrhs)
   endif
   do i=1,grid%system_size
      if (grid%vertex_type(vert2,i) == DIRICHLET) then
         grid%vertex_solution(vert2,i,:) = bcrhs(i)
      else
         grid%vertex_solution(vert2,i,:) = 0.0_my_real
      endif
   end do
   if (.not. present(elem_child_lid)) then
      if (any(grid%vertex_type(vert1,:) == PERIODIC_MASTER)) then
         grid%vertex(vert2)%next = vert1
         grid%head_level_vert(grid%element(child1)%level) = vert2
         grid%vertex(vert1)%previous = vert2
         grid%vertex(vert2)%previous = END_OF_LIST
      else
         grid%vertex(vert2)%next = grid%vertex(vert1)%next
         grid%vertex(vert1)%next = vert2
         grid%vertex(vert2)%previous = vert1
         if (grid%vertex(vert2)%next /= END_OF_LIST) then
            grid%vertex(grid%vertex(vert2)%next)%previous = vert2
         endif
      endif
   endif
   if (grid%have_true) then
      do j=1,max(1,grid%num_eval)
         do i=1,grid%system_size
            grid%vertex_exact(vert2,i,j)=trues(grid%vertex(vert2)%coord%x, &
                                               grid%vertex(vert2)%coord%y, &
                                               i,j)
         end do
      end do
   endif
   
endif

! new edges

paredge = grid%element(parent)%edge(3)
if (mate /= BOUNDARY) mateedge = grid%element(mate)%edge(3)

! edge 1

grid%edge(edge1)%gid = 4*grid%edge(paredge)%gid
grid%edge(edge1)%vertex = (/ grid%edge(paredge)%vertex(1),vert1 /)
grid%edge(edge1)%bmark = grid%edge(paredge)%bmark
grid%edge(edge1)%degree = grid%edge(paredge)%degree
grid%edge_type(edge1,:) = grid%edge_type(paredge,:)
if (mate == BOUNDARY) then
   grid%edge(edge1)%assoc_elem = child1
elseif (grid%element(child1)%gid < grid%element(child3)%gid) then
   grid%edge(edge1)%assoc_elem = child1
else
   grid%edge(edge1)%assoc_elem = child3
endif
if (grid%edge(edge1)%degree > 1) then
   allocate(grid%edge(edge1)%solution(grid%edge(edge1)%degree-1, &
            grid%system_size,max(1,grid%num_eval)),stat=astat)
   if (grid%have_true) then
      allocate(grid%edge(edge1)%exact(grid%edge(edge1)%degree-1, &
               grid%system_size,max(1,grid%num_eval)),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("memory allocation failed in bisect_triangle_pair")
         return
      endif
      grid%edge(edge1)%exact = 0.0_my_real
   else
      nullify(grid%edge(edge1)%exact)
   endif
   do i=1,grid%system_size
      if (grid%edge_type(edge1,i) == DIRICHLET) then
! TEMP OPENMP thread safe as long as the solution at existing vertices is
!             not changing (for example, by initial guess)
         call edge_exact(grid,edge1,i,"d")
      else
         grid%edge(edge1)%solution(:,i,:) = 0.0_my_real
      endif
      if (grid%have_true) call edge_exact(grid,edge1,i,"t")
   end do
else
   nullify(grid%edge(edge1)%solution)
   nullify(grid%edge(edge1)%exact)
endif

if (grid%edge(edge1)%degree >= 2) then
   if (present(delta_dof)) then
      delta_dof = delta_dof + grid%system_size*(grid%edge(edge1)%degree - 1)
      if (grid%element(grid%edge(edge1)%assoc_elem)%iown) &
         delta_dof_own = delta_dof_own + &
                         grid%system_size*(grid%edge(edge1)%degree - 1)
   else
      grid%dof = grid%dof + grid%system_size*(grid%edge(edge1)%degree - 1)
      if (grid%element(grid%edge(edge1)%assoc_elem)%iown) &
         grid%dof_own = grid%dof_own + &
                        grid%system_size*(grid%edge(edge1)%degree - 1)
   endif
endif

if (.not. present(elem_child_lid)) then
   call hash_insert(grid%edge(edge1)%gid,edge1,grid%edge_hash)
endif

! edge 2

grid%edge(edge2)%gid = 4*grid%edge(paredge)%gid + 1
grid%edge(edge2)%vertex = (/ grid%edge(paredge)%vertex(2),vert1 /)
grid%edge(edge2)%bmark = grid%edge(paredge)%bmark
grid%edge(edge2)%degree = grid%edge(paredge)%degree
grid%edge_type(edge2,:) = grid%edge_type(paredge,:)
if (mate == BOUNDARY) then
   grid%edge(edge2)%assoc_elem = child2
elseif (grid%element(child1)%gid < grid%element(child3)%gid) then
   grid%edge(edge2)%assoc_elem = child2
else
   grid%edge(edge2)%assoc_elem = child4
endif
if (grid%edge(edge2)%degree > 1) then
   allocate(grid%edge(edge2)%solution(grid%edge(edge2)%degree-1, &
            grid%system_size,max(1,grid%num_eval)),stat=astat)
   if (grid%have_true) then
      allocate(grid%edge(edge2)%exact(grid%edge(edge2)%degree-1, &
               grid%system_size,max(1,grid%num_eval)),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("memory allocation failed in bisect_triangle_pair")
         return
      endif
      grid%edge(edge2)%exact = 0.0_my_real
   else
      nullify(grid%edge(edge2)%exact)
   endif
   do i=1,grid%system_size
      if (grid%edge_type(edge2,i) == DIRICHLET) then
! TEMP OPENMP thread safe as long as the solution at existing vertices is
!             not changing (for example, by initial guess)
         call edge_exact(grid,edge2,i,"d")
      else
         grid%edge(edge2)%solution(:,i,:) = 0.0_my_real
      endif
      if (grid%have_true) call edge_exact(grid,edge2,i,"t")
   end do
else
   nullify(grid%edge(edge2)%solution)
   nullify(grid%edge(edge2)%exact)
endif
if (.not. present(elem_child_lid)) then
   call hash_insert(grid%edge(edge2)%gid,edge2,grid%edge_hash)
endif

! edge 3

if (mate /= BOUNDARY) then
   if (masterparent == parent) then
      if (grid%element(parent)%gid < grid%element(mate)%gid) then
         grid%edge(edge3)%gid = 4*grid%edge(paredge)%gid + 2
      else
         grid%edge(edge3)%gid = 4*grid%edge(paredge)%gid + 3
      endif
   else
      if (grid%element(parent)%gid < grid%element(mate)%gid) then
         grid%edge(edge3)%gid = 4*grid%edge(mateedge)%gid + 2
      else
         grid%edge(edge3)%gid = 4*grid%edge(mateedge)%gid + 3
      endif
   endif
else
   grid%edge(edge3)%gid = 4*grid%edge(paredge)%gid + 2
endif
grid%edge(edge3)%vertex = (/ grid%element(parent)%vertex(3),vert1 /)
grid%edge(edge3)%bmark = 0
grid%edge(edge3)%degree = grid%element(parent)%degree
grid%edge_type(edge3,:) = INTERIOR
grid%edge(edge3)%assoc_elem = child1
if (grid%edge(edge3)%degree > 1) then
   allocate(grid%edge(edge3)%solution(grid%edge(edge3)%degree-1, &
            grid%system_size,max(1,grid%num_eval)),stat=astat)
   grid%edge(edge3)%solution = 0.0_my_real
   if (grid%have_true) then
      allocate(grid%edge(edge3)%exact(grid%edge(edge3)%degree-1, &
               grid%system_size,max(1,grid%num_eval)),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("memory allocation failed in bisect_triangle_pair")
         return
      endif
      grid%edge(edge3)%exact = 0.0_my_real
      do i=1,grid%system_size
         call edge_exact(grid,edge3,i,"t")
      end do
   else
      nullify(grid%edge(edge3)%exact)
   endif
else
   nullify(grid%edge(edge3)%solution)
   nullify(grid%edge(edge3)%exact)
endif
if (grid%edge(edge3)%degree >= 2) then
   if (present(delta_dof)) then
      delta_dof = delta_dof + grid%system_size*(grid%edge(edge3)%degree - 1)
      if (grid%element(grid%edge(edge3)%assoc_elem)%iown) &
         delta_dof_own = delta_dof_own + &
                         grid%system_size*(grid%edge(edge3)%degree - 1)
   else
      grid%dof = grid%dof + grid%system_size*(grid%edge(edge3)%degree - 1)
      if (grid%element(grid%edge(edge3)%assoc_elem)%iown) &
         grid%dof_own = grid%dof_own + &
                        grid%system_size*(grid%edge(edge3)%degree - 1)
   endif
endif

if (.not. present(elem_child_lid)) then
   call hash_insert(grid%edge(edge3)%gid,edge3,grid%edge_hash)
endif

! edge 4

if (mate /= BOUNDARY) then
   if (masterparent == parent) then
      if (grid%element(parent)%gid < grid%element(mate)%gid) then
         grid%edge(edge4)%gid = 4*grid%edge(paredge)%gid + 3
      else
         grid%edge(edge4)%gid = 4*grid%edge(paredge)%gid + 2
      endif
   else
      if (grid%element(parent)%gid < grid%element(mate)%gid) then
         grid%edge(edge4)%gid = 4*grid%edge(mateedge)%gid + 3
      else
         grid%edge(edge4)%gid = 4*grid%edge(mateedge)%gid + 2
      endif
   endif
   if (any(grid%edge_type(edge1,:) == PERIODIC_MASTER) .or. &
       any(grid%edge_type(edge1,:) == PERIODIC_SLAVE)) then
      grid%edge(edge4)%vertex = (/ grid%element(mate)%vertex(3),vert2 /)
   else
      grid%edge(edge4)%vertex = (/ grid%element(mate)%vertex(3),vert1 /)
   endif
   grid%edge(edge4)%bmark = 0
   grid%edge(edge4)%degree = grid%element(mate)%degree
   grid%edge_type(edge4,:) = INTERIOR
   grid%edge(edge4)%assoc_elem = child3
   if (grid%edge(edge4)%degree > 1) then
      allocate(grid%edge(edge4)%solution(grid%edge(edge4)%degree-1, &
               grid%system_size,max(1,grid%num_eval)),stat=astat)
      grid%edge(edge4)%solution = 0.0_my_real
      if (grid%have_true) then
         allocate(grid%edge(edge4)%exact(grid%edge(edge4)%degree-1, &
                  grid%system_size,max(1,grid%num_eval)),stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("memory allocation failed in bisect_triangle_pair")
            return
         endif
         grid%edge(edge4)%exact = 0.0_my_real
         do i=1,grid%system_size
            call edge_exact(grid,edge4,i,"t")
         end do
      else
         nullify(grid%edge(edge4)%exact)
      endif
   else
      nullify(grid%edge(edge4)%solution)
      nullify(grid%edge(edge4)%exact)
   endif
   if (grid%edge(edge4)%degree >= 2) then
      if (present(delta_dof)) then
         delta_dof = delta_dof + grid%system_size*(grid%edge(edge4)%degree - 1)
         if (grid%element(grid%edge(edge4)%assoc_elem)%iown) &
            delta_dof_own = delta_dof_own + &
                            grid%system_size*(grid%edge(edge4)%degree - 1)
      else
         grid%dof = grid%dof + grid%system_size*(grid%edge(edge4)%degree - 1)
         if (grid%element(grid%edge(edge4)%assoc_elem)%iown) &
            grid%dof_own = grid%dof_own + &
                           grid%system_size*(grid%edge(edge4)%degree - 1)
      endif
   endif
   if (.not. present(elem_child_lid)) then
      call hash_insert(grid%edge(edge4)%gid,edge4,grid%edge_hash)
   endif
endif

! 5th and 6th edges if periodic boundary edge

if (any(grid%edge_type(edge1,:) == PERIODIC_MASTER) .or. &
    any(grid%edge_type(edge1,:) == PERIODIC_SLAVE)) then

   if (.not. present(elem_child_lid)) then
      if (grid%next_free_edge == END_OF_LIST) then
         call more_edges(grid,errcode)
         if (errcode /= 0) return
      endif
      edge5 = grid%next_free_edge
      grid%next_free_edge = grid%edge(edge5)%next
      if (grid%next_free_edge == END_OF_LIST) then
         call more_edges(grid,errcode)
         if (errcode /= 0) return
      endif
      edge6 = grid%next_free_edge
      grid%next_free_edge = grid%edge(edge6)%next
   endif

   grid%element(child3)%edge(2) = edge5
   grid%element(child4)%edge(2) = edge6

! edge 5

   grid%edge(edge5)%gid = 4*grid%edge(mateedge)%gid
   grid%edge(edge5)%vertex = (/ grid%edge(mateedge)%vertex(1),vert2 /)
   grid%edge(edge5)%bmark = grid%edge(mateedge)%bmark
   grid%edge(edge5)%degree = grid%edge(mateedge)%degree
   grid%edge_type(edge5,:) = grid%edge_type(mateedge,:)
   grid%edge(edge5)%assoc_elem = grid%edge(edge1)%assoc_elem
   if (.not. present(elem_child_lid)) then
      if (any(grid%edge_type(edge5,:) == PERIODIC_SLAVE)) then
         grid%edge(edge5)%next = edge1
      else
         grid%edge(edge1)%next = edge5
      endif
   endif
   if (grid%edge(edge5)%degree > 1) then
      allocate(grid%edge(edge5)%solution(grid%edge(edge5)%degree-1, &
               grid%system_size,max(1,grid%num_eval)),stat=astat)
      if (grid%have_true) then
         allocate(grid%edge(edge5)%exact(grid%edge(edge5)%degree-1, &
                  grid%system_size,max(1,grid%num_eval)),stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("memory allocation failed in bisect_triangle_pair")
            return
         endif
         grid%edge(edge5)%exact = 0.0_my_real
      else
         nullify(grid%edge(edge5)%exact)
      endif
      do i=1,grid%system_size
         if (grid%edge_type(edge5,i) == DIRICHLET) then
! TEMP OPENMP thread safe as long as the solution at existing vertices is
!             not changing (for example, by initial guess)
            call edge_exact(grid,edge5,i,"d")
         else
            grid%edge(edge5)%solution(:,i,:) = 0.0_my_real
         endif
         if (grid%have_true) call edge_exact(grid,edge5,i,"t")
      end do
   else
      nullify(grid%edge(edge5)%solution)
      nullify(grid%edge(edge5)%exact)
   endif
   if (.not. present(elem_child_lid)) then
      call hash_insert(grid%edge(edge5)%gid,edge5,grid%edge_hash)
   endif

! edge 6

   grid%edge(edge6)%gid = 4*grid%edge(mateedge)%gid + 1
   grid%edge(edge6)%vertex = (/ grid%edge(mateedge)%vertex(2),vert2 /)
   grid%edge(edge6)%bmark = grid%edge(mateedge)%bmark
   grid%edge(edge6)%degree = grid%edge(mateedge)%degree
   grid%edge_type(edge6,:) = grid%edge_type(mateedge,:)
   grid%edge(edge6)%assoc_elem = grid%edge(edge2)%assoc_elem
   if (.not. present(elem_child_lid)) then
      if (any(grid%edge_type(edge6,:) == PERIODIC_SLAVE)) then
         grid%edge(edge6)%next = edge2
      else
         grid%edge(edge2)%next = edge6
      endif
   endif
   if (grid%edge(edge6)%degree > 1) then
      allocate(grid%edge(edge6)%solution(grid%edge(edge6)%degree-1, &
               grid%system_size,max(1,grid%num_eval)),stat=astat)
      if (grid%have_true) then
         allocate(grid%edge(edge6)%exact(grid%edge(edge6)%degree-1, &
                  grid%system_size,max(1,grid%num_eval)),stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("memory allocation failed in bisect_triangle_pair")
            return
         endif
         grid%edge(edge6)%exact = 0.0_my_real
      else
         nullify(grid%edge(edge6)%exact)
      endif
      do i=1,grid%system_size
         if (grid%edge_type(edge6,i) == DIRICHLET) then
! TEMP OPENMP thread safe as long as the solution at existing vertices is
!             not changing (for example, by initial guess)
            call edge_exact(grid,edge6,i,"d")
         else
            grid%edge(edge6)%solution(:,i,:) = 0.0_my_real
         endif
         if (grid%have_true) call edge_exact(grid,edge6,i,"t")
      end do
   else
      nullify(grid%edge(edge6)%solution)
      nullify(grid%edge(edge6)%exact)
   endif
   if (.not. present(elem_child_lid)) then
      call hash_insert(grid%edge(edge6)%gid,edge6,grid%edge_hash)
   endif

endif

! if reduce_p is present, the children generated for compatibility should have
! their degree reduced by 1/sqrt(2).  This is for the REFSOLN_ELEM hp strategy.
! reduce_p is passed in as 0 and incremented with each recursive call.  Thus
! reduce_p==0 indicates the original call, so the reduction does not occur.
! reduce_p_max holds the largest value of reduce_p in this recursive chain.
! If reduce_p is not reduce_p_max, then the children of the mate are not
! reduced, because they were reduced when the mate's parent was refined.
! TEMP OPENMP I'm not sure how this works when I do one level at a time.  I may
!             need to make reduce_p array that gets set by the recursion in
!             creating the element list

if (present(reduce_p)) then
   if (reduce_p > 0) then
      newdeg = floor((grid%element(child1)%degree+1)/sqrt(2.0_my_real))
      do while (grid%element(child1)%degree > newdeg)
! TEMP OPENMP not thread safe until p_coarsen_elem is
         call p_coarsen_elem(grid,child1,errcode,refcont)
      end do
      do while (grid%element(child2)%degree > newdeg)
! TEMP OPENMP not thread safe until p_coarsen_elem is
         call p_coarsen_elem(grid,child2,errcode,refcont)
      end do
      if ((.not. mate == BOUNDARY) .and. reduce_p == reduce_p_max) then
         newdeg = floor((grid%element(child3)%degree+1)/sqrt(2.0_my_real))
         do while (grid%element(child3)%degree > newdeg)
! TEMP OPENMP not thread safe until p_coarsen_elem is
            call p_coarsen_elem(grid,child3,errcode,refcont)
         end do
         do while (grid%element(child4)%degree > newdeg)
! TEMP OPENMP not thread safe until p_coarsen_elem is
            call p_coarsen_elem(grid,child4,errcode,refcont)
         end do
      endif
   endif
endif

! set exact solution for element bases (now that edges are done)

if (grid%have_true) then
   do i=1,grid%system_size
      call elem_exact(grid,child1,i,"t")
      call elem_exact(grid,child2,i,"t")
      if (mate /= BOUNDARY) then
         call elem_exact(grid,child3,i,"t")
         call elem_exact(grid,child4,i,"t")
      endif
   end do
endif

! set the in and out vertices

call set_inout(grid,parent,child1,child2)
if (mate /= BOUNDARY) then
   call set_inout(grid,mate,child3,child4)
endif

! associated element for the vertices of the parent and mate
! TEMP OPENMP thread safe if assignment is atomic

if (grid%vertex(grid%element(parent)%vertex(1))%assoc_elem == parent) then
   grid%vertex(grid%element(parent)%vertex(1))%assoc_elem = child1
endif
if (grid%vertex(grid%element(parent)%vertex(2))%assoc_elem == parent) then
   grid%vertex(grid%element(parent)%vertex(2))%assoc_elem = child2
endif
if (grid%vertex(grid%element(parent)%vertex(3))%assoc_elem == parent) then
   grid%vertex(grid%element(parent)%vertex(3))%assoc_elem = child1
endif
if (mate /= BOUNDARY) then
   if (grid%vertex(grid%element(mate)%vertex(1))%assoc_elem == mate) then
      grid%vertex(grid%element(mate)%vertex(1))%assoc_elem = child3
   endif
   if (grid%vertex(grid%element(mate)%vertex(2))%assoc_elem == mate) then
      grid%vertex(grid%element(mate)%vertex(2))%assoc_elem = child4
   endif
   if (grid%vertex(grid%element(mate)%vertex(3))%assoc_elem == mate) then
      grid%vertex(grid%element(mate)%vertex(3))%assoc_elem = child3
   endif
endif

! associated elements of the edges of the parent and mate

if (grid%edge(grid%element(parent)%edge(1))%assoc_elem == parent) then
   grid%edge(grid%element(parent)%edge(1))%assoc_elem = child2
endif
if (grid%edge(grid%element(parent)%edge(2))%assoc_elem == parent) then
   grid%edge(grid%element(parent)%edge(2))%assoc_elem = child1
endif
if (grid%edge(grid%element(parent)%edge(3))%assoc_elem == parent) then
   grid%edge(grid%element(parent)%edge(3))%assoc_elem = child1
endif
if (mate /= BOUNDARY) then
   if (grid%edge(grid%element(mate)%edge(1))%assoc_elem == mate) then
      grid%edge(grid%element(mate)%edge(1))%assoc_elem = child4
   endif
   if (grid%edge(grid%element(mate)%edge(2))%assoc_elem == mate) then
      grid%edge(grid%element(mate)%edge(2))%assoc_elem = child3
   endif
   if (grid%edge(grid%element(mate)%edge(3))%assoc_elem == mate) then
      grid%edge(grid%element(mate)%edge(3))%assoc_elem = child3
   endif
endif

! check for periodic slave vertices, and set the assoc_elem to the master's
! TEMP OPENMP this might not be thread safe

do i=1,3
   if (any(grid%vertex_type(grid%element(parent)%vertex(i),:) == PERIODIC_SLAVE)) then
      grid%vertex(grid%element(parent)%vertex(i))%assoc_elem = &
        grid%vertex(grid%vertex(grid%element(parent)%vertex(i))%next)%assoc_elem
   endif
   if (mate /= BOUNDARY) then
    if (any(grid%vertex_type(grid%element(mate)%vertex(i),:) == PERIODIC_SLAVE)) then
      grid%vertex(grid%element(mate)%vertex(i))%assoc_elem = &
        grid%vertex(grid%vertex(grid%element(mate)%vertex(i))%next)%assoc_elem
    endif
   endif
end do

! grid scalars

if (present(delta_dof)) then

   if (mate == BOUNDARY) then
      delta_nelem = delta_nelem + 2
      delta_nelem_leaf = delta_nelem_leaf + 1
      if (grid%element(parent)%iown) &
         delta_nelem_leaf_own = delta_nelem_leaf_own + 1
      delta_nedge = delta_nedge + 3
! TEMP also need to change nedge_own
   else
      delta_nelem = delta_nelem + 4
      delta_nelem_leaf = delta_nelem_leaf + 2
      if (grid%element(parent)%iown) &
         delta_nelem_leaf_own = delta_nelem_leaf_own + 1
      if (grid%element(mate)%iown) &
         delta_nelem_leaf_own = delta_nelem_leaf_own + 1
      delta_nedge = delta_nedge + 4
! TEMP also need to change nedge_own
   endif
   delta_nvert = delta_nvert + 1
   if (grid%element(grid%vertex(vert1)%assoc_elem)%iown) then
      delta_nvert_own = delta_nvert_own + 1
   endif
   max_nlev = max(max_nlev,grid%element(child1)%level)
   if (any(grid%vertex_type(vert1,:) == PERIODIC_MASTER) .or. &
       any(grid%vertex_type(vert1,:) == PERIODIC_SLAVE)) then
      delta_nvert = delta_nvert + 1
      if (grid%element(grid%vertex(vert1)%assoc_elem)%iown) then
         delta_nvert_own = delta_nvert_own + 1
      endif
   endif
   if (any(grid%edge_type(grid%element(parent)%edge(3),:) == PERIODIC_MASTER) .or. &
       any(grid%edge_type(grid%element(parent)%edge(3),:) == PERIODIC_SLAVE)) then
      delta_nedge = delta_nedge + 2
! TEMP also need to change nedge_own
   endif

else ! deltas are not present

   if (mate == BOUNDARY) then
      grid%nelem = grid%nelem + 2
      grid%nelem_leaf = grid%nelem_leaf + 1
      if (grid%element(parent)%iown) &
         grid%nelem_leaf_own = grid%nelem_leaf_own + 1
      grid%nedge = grid%nedge + 3
! TEMP also need to change nedge_own
   else
      grid%nelem = grid%nelem + 4
      grid%nelem_leaf = grid%nelem_leaf + 2
      if (grid%element(parent)%iown) &
         grid%nelem_leaf_own = grid%nelem_leaf_own + 1
      if (grid%element(mate)%iown) &
         grid%nelem_leaf_own = grid%nelem_leaf_own + 1
      grid%nedge = grid%nedge + 4
! TEMP also need to change nedge_own
   endif
   grid%nvert = grid%nvert + 1
   if (grid%element(grid%vertex(vert1)%assoc_elem)%iown) then
      grid%nvert_own = grid%nvert_own + 1
   endif
   grid%nlev = max(grid%nlev,grid%element(child1)%level)
   if (any(grid%vertex_type(vert1,:) == PERIODIC_MASTER) .or. &
       any(grid%vertex_type(vert1,:) == PERIODIC_SLAVE)) then
      grid%nvert = grid%nvert + 1
      if (grid%element(grid%vertex(vert1)%assoc_elem)%iown) then
         grid%nvert_own = grid%nvert_own + 1
      endif
   endif
   if (any(grid%edge_type(grid%element(parent)%edge(3),:) == PERIODIC_MASTER) .or. &
       any(grid%edge_type(grid%element(parent)%edge(3),:) == PERIODIC_SLAVE)) then
      grid%nedge = grid%nedge + 2
! TEMP also need to change nedge_own
   endif

endif ! are deltas present

! set initial guess for solution

if (refcont%error_estimator == INITIAL_CONDITION) then
! TEMP OPENMP may require atomic assignment
   call init_guess_ic(grid,parent,(/child1,child2,child3,child4/))
else
! TEMP OPENMP thread safe if elemental_matrix is
   if (mate == BOUNDARY) then
      call init_guess_h(grid,parent,mate,(/child1,child2,BOUNDARY,BOUNDARY/))
   else
      call init_guess_h(grid,parent,mate,(/child1,child2,child3,child4/))
   endif
endif

! Set child error indicators and
! add children to the error indicator lists, if present and control wants it

if (present(elist)) then

   if (present(numhref)) then
      if (numhref(parent) > 0) then
         numhref(child1) = numhref(parent)-1
         numhref(child2) = numhref(parent)-1
         grid%element_errind(child1,:) = &
                      elist%max_errind*grid%element(child1)%work/binw + .01
         if (grid%element(child1)%mate == BOUNDARY) then
            grid%element(child1)%work = ((grid%element(child1)%degree)* &
                    (grid%element(child1)%degree+1))/2
         else
            grid%element(child1)%work = grid%element(child1)%degree**2
         endif
         grid%element_errind(child2,:) = &
                      elist%max_errind*grid%element(child2)%work/binw + .01
         if (grid%element(child2)%mate == BOUNDARY) then
            grid%element(child2)%work = ((grid%element(child2)%degree)* &
                    (grid%element(child2)%degree+1))/2
         else
            grid%element(child2)%work = grid%element(child2)%degree**2
         endif
         add_to_list = numhref(child1) > 0 .or. numpref(child1) > 0
         need_new_errind = .false.
      elseif (numhref(parent) < 0) then
         if (present(return_to_elist)) then
            add_to_list = return_to_elist
         else
            add_to_list = .false.
         endif
         need_new_errind = add_to_list
      else
         add_to_list = .false.
         need_new_errind = add_to_list
      endif
   elseif (present(return_to_elist)) then
      add_to_list = return_to_elist
      need_new_errind = add_to_list
   else
      add_to_list = .false.
      need_new_errind = add_to_list
   endif
   if (refcont%hp_strategy == HP_SMOOTH_PRED) then
      need_new_errind=.true.
   endif
   if (need_new_errind) then
! TEMP OPENMP haven't checked error_indicator for thread safety
!             EXPLICIT cannot be used during refinement because it needs the
!                      neighbors for the jump term, and they may or may not have
!                      been refined
!             Same should be true for EQUILIBRATED_RESIDUAL and LOCAL_P
!             HIERARCHICAL_COEF appears to be OK
!             I haven't checked the others.
      call error_indicator(grid,child1,refcont%error_estimator, &
                        grid%element_errind(child1,:),grid%element(child1)%work)
      call error_indicator(grid,child2,refcont%error_estimator, &
                        grid%element_errind(child2,:),grid%element(child2)%work)
   endif
! OPENMP not thread safe; requires add_to_list evaluate to false
   if (add_to_list) then
      call add_to_errind_list(child1,elist,grid,refcont)
      call add_to_errind_list(child2,elist,grid,refcont)
   endif

   if (.not. mate == BOUNDARY) then
      if (present(numhref)) then
         if (numhref(mate) > 0) then
            numhref(child3) = numhref(mate)-1
            numhref(child4) = numhref(mate)-1
            grid%element_errind(child3,:) = &
                        elist%max_errind*grid%element(child3)%work/binw + .01
            if (grid%element(child3)%mate == BOUNDARY) then
               grid%element(child3)%work = ((grid%element(child3)%degree)* &
                       (grid%element(child3)%degree+1))/2
            else
               grid%element(child3)%work = grid%element(child3)%degree**2
            endif
            grid%element_errind(child4,:) = &
                        elist%max_errind*grid%element(child4)%work/binw + .01
            if (grid%element(child4)%mate == BOUNDARY) then
               grid%element(child4)%work = ((grid%element(child4)%degree)* &
                       (grid%element(child4)%degree+1))/2
            else
               grid%element(child4)%work = grid%element(child4)%degree**2
            endif
            add_to_list = numhref(child3) > 0 .or. numpref(child3) > 0
            need_new_errind = .false.
         elseif (numhref(mate) < 0) then
            if (present(return_to_elist)) then
               add_to_list = return_to_elist
            else
               add_to_list = .false.
            endif
            need_new_errind = add_to_list
         else
            add_to_list = .false.
            need_new_errind = add_to_list
         endif
      elseif (present(return_to_elist)) then
         add_to_list = return_to_elist
         need_new_errind = add_to_list
      else
         add_to_list = .false.
         need_new_errind = add_to_list
      endif
      if (refcont%hp_strategy == HP_SMOOTH_PRED) then
         need_new_errind=.true.
      endif
      if (need_new_errind) then
! TEMP OPENMP haven't checked error_indicator for thread safety
         call error_indicator(grid,child3,refcont%error_estimator, &
                              grid%element_errind(child3,:), &
                              grid%element(child3)%work)
         call error_indicator(grid,child4,refcont%error_estimator, &
                              grid%element_errind(child4,:), &
                              grid%element(child4)%work)
      endif
! OPENMP not thread safe; requires add_to_list evaluate to false
      if (add_to_list) then
         call add_to_errind_list(child3,elist,grid,refcont)
         call add_to_errind_list(child4,elist,grid,refcont)
      endif
   endif

endif

! determine the type of refinement later

if (present(reftype)) then
   reftype(parent) = "u"
   reftype(child1) = "u"
   reftype(child2) = "u"
   if (.not. mate == BOUNDARY) then
      reftype(mate) = "u"
      reftype(child3) = "u"
      reftype(child4) = "u"
   endif
endif

end subroutine bisect_triangle_pair

!          --------------
subroutine after_h_refine(grid,element_list,nelem,vert_lid,edge_lid,elem_lid)
!          --------------

!----------------------------------------------------
! This routine performs parts of h refinement that must be done by a single
! OpenMP thread after the OpenMP-parallel refinement.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: element_list(:), nelem, vert_lid(:,:), edge_lid(:,:), &
                       elem_lid(:,:)
!----------------------------------------------------
! Local variables:

integer :: i, elem
!----------------------------------------------------
! Begin executable code

! for each element in the list ...

do i=1,nelem
   elem = element_list(i)

! add the children to the hash tables

   if (elem_lid(1,i) /= -1) &
      call hash_insert(grid%element(elem_lid(1,i))%gid,elem_lid(1,i), &
                       grid%elem_hash)
   if (elem_lid(2,i) /= -1) &
      call hash_insert(grid%element(elem_lid(2,i))%gid,elem_lid(2,i), &
                       grid%elem_hash)
   if (elem_lid(3,i) /= -1) &
      call hash_insert(grid%element(elem_lid(3,i))%gid,elem_lid(3,i), &
                       grid%elem_hash)
   if (elem_lid(4,i) /= -1) &
      call hash_insert(grid%element(elem_lid(4,i))%gid,elem_lid(4,i), &
                       grid%elem_hash)
   if (edge_lid(1,i) /= -1) &
      call hash_insert(grid%edge(edge_lid(1,i))%gid,edge_lid(1,i), &
                                 grid%edge_hash)
   if (edge_lid(2,i) /= -1) &
      call hash_insert(grid%edge(edge_lid(2,i))%gid,edge_lid(2,i), &
                                 grid%edge_hash)
   if (edge_lid(3,i) /= -1) &
      call hash_insert(grid%edge(edge_lid(3,i))%gid,edge_lid(3,i), &
                                 grid%edge_hash)
   if (edge_lid(4,i) /= -1) &
      call hash_insert(grid%edge(edge_lid(4,i))%gid,edge_lid(4,i), &
                                 grid%edge_hash)
   if (edge_lid(5,i) /= -1) &
      call hash_insert(grid%edge(edge_lid(5,i))%gid,edge_lid(5,i), &
                                 grid%edge_hash)
   if (edge_lid(6,i) /= -1) &
      call hash_insert(grid%edge(edge_lid(6,i))%gid,edge_lid(6,i), &
                                 grid%edge_hash)
   if (vert_lid(1,i) /= -1) &
      call hash_insert(grid%vertex(vert_lid(1,i))%gid,vert_lid(1,i), &
                       grid%vert_hash)
   if (vert_lid(2,i) /= -1) &
      call hash_insert(grid%vertex(vert_lid(2,i))%gid,vert_lid(2,i), &
                       grid%vert_hash)

end do ! elements

end subroutine after_h_refine

!          ---------
subroutine set_inout(grid,parent,child1,child2)
!          ---------

!----------------------------------------------------
! This routine sets the in and out vertices and order of child1 and child2
! RESTRICTION bisected triangles
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: parent, child1, child2
!----------------------------------------------------
! Local variables:

type(grid_type), intent(inout) :: grid
integer :: parent_in, parent_out, i
!----------------------------------------------------
! Begin executable code

! determine the elemental indices of the in and out vertices of the parent

parent_in = 0
parent_out = 0
do i=1,VERTICES_PER_ELEMENT
   if (grid%element(parent)%vertex(i) == grid%element(parent)%in) parent_in = i
   if (grid%element(parent)%vertex(i) == grid%element(parent)%out) parent_out = i
end do

! set the in and out of the children and order of the children based on
! which vertices are the in and out of the parent

select case(parent_in)
case (1)
   select case(parent_out)
   case(1)
      call fatal("element has the same vertex for in and out")
      stop
   case(2)
      grid%element(child1)%in  = grid%element(child1)%vertex(1)
      grid%element(child1)%out = grid%element(child1)%vertex(2)
      grid%element(child2)%in  = grid%element(child2)%vertex(2)
      grid%element(child2)%out = grid%element(child2)%vertex(1)
      grid%element(parent)%order = (/1,2/)
   case(3)
      grid%element(child1)%in  = grid%element(child1)%vertex(1)
      grid%element(child1)%out = grid%element(child1)%vertex(3)
      grid%element(child2)%in  = grid%element(child2)%vertex(3)
      grid%element(child2)%out = grid%element(child2)%vertex(2)
      grid%element(parent)%order = (/1,2/)
   case default
      call fatal("couldn't find the out vertex", &
                 intlist=(/parent,grid%element(parent)%out/))
      stop
   end select
case (2)
   select case(parent_out)
   case(1)
      grid%element(child2)%in  = grid%element(child2)%vertex(1)
      grid%element(child2)%out = grid%element(child2)%vertex(2)
      grid%element(child1)%in  = grid%element(child1)%vertex(2)
      grid%element(child1)%out = grid%element(child1)%vertex(1)
      grid%element(parent)%order = (/2,1/)
   case(2)
      call fatal("element has the same vertex for in and out")
      stop
   case(3)
      grid%element(child2)%in  = grid%element(child2)%vertex(1)
      grid%element(child2)%out = grid%element(child2)%vertex(3)
      grid%element(child1)%in  = grid%element(child1)%vertex(3)
      grid%element(child1)%out = grid%element(child1)%vertex(2)
      grid%element(parent)%order = (/2,1/)
   case default
      call fatal("couldn't find the out vertex", &
                 intlist=(/parent,grid%element(parent)%out/))
      stop
   end select
case (3)
   select case(parent_out)
   case(1)
      grid%element(child2)%in  = grid%element(child2)%vertex(2)
      grid%element(child2)%out = grid%element(child2)%vertex(3)
      grid%element(child1)%in  = grid%element(child1)%vertex(3)
      grid%element(child1)%out = grid%element(child1)%vertex(1)
      grid%element(parent)%order = (/2,1/)
   case(2)
      grid%element(child1)%in  = grid%element(child1)%vertex(2)
      grid%element(child1)%out = grid%element(child1)%vertex(3)
      grid%element(child2)%in  = grid%element(child2)%vertex(3)
      grid%element(child2)%out = grid%element(child2)%vertex(1)
      grid%element(parent)%order = (/1,2/)
   case(3)
      call fatal("element has the same vertex for in and out")
      stop
   case default
      call fatal("couldn't find the out vertex", &
                 intlist=(/parent,grid%element(parent)%out/))
      stop
   end select
case default
   call fatal("couldn't find the in vertex", &
              intlist=(/parent,grid%element(parent)%in/))
   stop
end select

! order of children shouldn't matter, but it's good to initialize

grid%element(child1)%order = (/1,2/)
grid%element(child2)%order = (/1,2/)

end subroutine set_inout

!                    --------------
recursive subroutine create_element(grid,elemgid,refine_control,err)
!                    --------------

!----------------------------------------------------
! This routine performs the refinements needed to create the element
! with global id elemgid
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(hash_key), intent(in) :: elemgid
type(refine_options), intent(in) :: refine_control
integer, intent(inout) :: err
!----------------------------------------------------
! Local variables:

integer :: parentlid
type(hash_key) :: parentgid
!----------------------------------------------------
! Begin executable code

err = 0

! return if the element already exists

if (hash_decode_key(elemgid,grid%elem_hash) /= HASH_NOT_FOUND) return

! identify the parent

parentgid = elemgid/MAX_CHILD
parentlid = hash_decode_key(parentgid,grid%elem_hash)

! if the parent doesn't exist, create it

if (parentlid == HASH_NOT_FOUND) then
   call create_element(grid,parentgid,refine_control,err)
endif

if (err /= 0) then
   call warning("refinement failed in create_element")
else

! refine the parent

   parentlid = hash_decode_key(parentgid,grid%elem_hash)
   call bisect_triangle_pair(grid,parentlid,err,refine_control)

endif

end subroutine create_element

!          -----------------------
subroutine remove_from_errind_list(elem,elist)
!          -----------------------

!----------------------------------------------------
! This routine removes elem from the error indicator lists
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: elem
type(errind_list) :: elist

!----------------------------------------------------
! Local variables:

integer :: i
!----------------------------------------------------
! Begin executable code

! make sure elem is on an error indicator list

if (elist%next_errind(elem) == NOT_ON_LIST) return

if (elist%prev_errind(elem) == END_OF_LIST) then
! elem is the head of a list
   do i=1,size(elist%head_errind)
      if (elist%head_errind(i) == elem) exit
   end do
   if (i > size(elist%head_errind)) then
      call warning("couldn't find appropriate an error indicator list in remove_from_errind_list")
      return
   else
      elist%head_errind(i) = elist%next_errind(elem)
      if (elist%next_errind(elem) == END_OF_LIST) then
! elem is both the head and the tail
         elist%tail_errind(i) = END_OF_LIST
      else
         elist%prev_errind(elist%next_errind(elem)) = END_OF_LIST
      endif
      elist%next_errind(elem) = NOT_ON_LIST
      elist%prev_errind(elem) = NOT_ON_LIST
   endif
elseif (elist%next_errind(elem) == END_OF_LIST) then
! elem is the tail of a list
   do i=1,size(elist%tail_errind)
      if (elist%tail_errind(i) == elem) exit
   end do
   if (i > size(elist%tail_errind)) then
      call warning("couldn't find appropriate an error indicator list in remove_from_errind_list")
      return
   else
      elist%tail_errind(i) = elist%prev_errind(elem)
      elist%next_errind(elist%prev_errind(elem)) = END_OF_LIST
      elist%next_errind(elem) = NOT_ON_LIST
      elist%prev_errind(elem) = NOT_ON_LIST
   endif
else
! elem is in the interior
   elist%next_errind(elist%prev_errind(elem)) = elist%next_errind(elem)
   elist%prev_errind(elist%next_errind(elem)) = elist%prev_errind(elem)
   elist%next_errind(elem) = NOT_ON_LIST
   elist%prev_errind(elem) = NOT_ON_LIST
endif

end subroutine remove_from_errind_list

!          ------------------
subroutine add_to_errind_list(elem,elist,grid,refcont)
!          ------------------

!----------------------------------------------------
! This routine adds element elem to the error indicator list
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: elem
type(errind_list), intent(inout) :: elist
type(grid_type), intent(in) :: grid
type(refine_options), intent(in) :: refcont
!----------------------------------------------------
! Local variables:

real(my_real) :: errind
integer :: i
!----------------------------------------------------
! Begin executable code

! don't add if I don't own it

if (.not. grid%element(elem)%iown) return

! make sure it is not already on a list

if (elist%next_errind(elem) /= NOT_ON_LIST) return

errind = maxval(grid%element_errind(elem,:))/grid%element(elem)%work
do i=elist%current_list,size(elist%head_errind)-1
   if (errind > elist%max_errind/(binw**i)) exit
end do
if (elist%head_errind(i) == END_OF_LIST) then
   elist%head_errind(i) = elem
   elist%tail_errind(i) = elem
   elist%next_errind(elem) = END_OF_LIST
   elist%prev_errind(elem) = END_OF_LIST
else
   elist%prev_errind(elem) = elist%tail_errind(i)
   elist%next_errind(elem) = END_OF_LIST
   elist%next_errind(elist%tail_errind(i)) = elem
   elist%tail_errind(i) = elem
endif

end subroutine add_to_errind_list

!          -------------
subroutine init_guess_ic(grid,elem,children)
!          -------------

!----------------------------------------------------
! This routine sets the solution in either element elem or the descendants
! of elem and its mate from the function in iconds.  The descendants are
! set if and only if children is present.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: elem
integer, optional, intent(in) :: children(4)
!----------------------------------------------------
! Local variables:

integer :: vert(5), edge(8), face(4)
integer :: nvert, nedge, nface, v, e, f, i, comp, eigen
!----------------------------------------------------
! Begin executable code

! set the entities to do

if (present(children)) then
   face(1:2) = children(1:2)
   edge(1:3) = grid%element(face(1))%edge
   edge(4:5) = grid%element(face(2))%edge(2:3)
   vert(1:3) = grid%element(elem)%vertex
   vert(4)   = grid%element(face(1))%vertex(3)
   if (grid%element(elem)%mate == BOUNDARY) then
      nvert = 4
      nedge = 5
      nface = 2
   else
      face(3:4) = children(3:4)
      edge(6:7) = grid%element(face(3))%edge(1:3:2)
      edge(8)   = grid%element(face(4))%edge(3)
      vert(5)   = grid%element(face(3))%vertex(2)
      nvert = 5
      nedge = 8
      nface = 4
   endif
else
   face(1)   = elem
   edge(1:3) = grid%element(elem)%edge
   vert(1:3) = grid%element(elem)%vertex
   nvert = 3
   nedge = 3
   nface = 1
endif

! set the vertex basis function coefficients

do i=1,nvert
   v = vert(i)
   do comp=1,grid%system_size
      if (grid%vertex_type(v,comp) /= DIRICHLET .and. &
          grid%vertex_type(v,comp) /= PERIODIC_MASTER_DIR .and. &
          grid%vertex_type(v,comp) /= PERIODIC_SLAVE_DIR) then
         do eigen=1,max(1,grid%num_eval)
            grid%vertex_solution(v,comp,eigen) = &
               iconds(grid%vertex(v)%coord%x, grid%vertex(v)%coord%y, &
                      comp,eigen)
         end do
      endif
   end do
end do

! set the edge basis function coefficients

do i=1,nedge
   e = edge(i)
   do comp=1,grid%system_size
      if (grid%edge_type(e,comp) /= DIRICHLET) then
         call edge_exact(grid,e,comp,"i")
      endif
   end do
end do

! set the face basis function coefficients

do i=1,nface
   f = face(i)
   do comp=1,grid%system_size
      call elem_exact(grid,f,comp,"i")
   end do
end do

end subroutine init_guess_ic

!          ------------
subroutine init_guess_p(grid,elem,old_edge_deg)
!          ------------

!----------------------------------------------------
! This routine sets the initial guess for the solution components associated
! with bases that have been added by the p refinement of element elem
! by solving a local Dirichlet residual problem for the red p-hierarchical bases
! over element elem using a domain consisting of elem and its three neighbors.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout), target :: grid
integer, intent(in) :: elem, old_edge_deg(:)
!----------------------------------------------------
! Local variables:

integer :: ss, i, j, k, degree(1+EDGES_PER_ELEMENT), bmark(EDGES_PER_ELEMENT), &
           edge_type(EDGES_PER_ELEMENT,grid%system_size), nred, n, astat, info,&
           neigh(3), v1, v2, nev, l
real(my_real) :: xvert(VERTICES_PER_ELEMENT), yvert(VERTICES_PER_ELEMENT), &
                 rmin
integer, allocatable :: ipiv(:)
real(my_real), allocatable :: mat(:,:), rs(:,:),elem_mat(:,:),elem_rs(:,:)
!----------------------------------------------------
! Begin executable code

ss = grid%system_size
nev = max(1,grid%num_eval)

! don't actually use rmin

rmin = 0.0_my_real

! the matrix is the elemental matrix for element elem plus a 1x1 system from
! each of the neighbors for the one red basis on each edge

xvert = grid%vertex(grid%element(elem)%vertex)%coord%x
yvert = grid%vertex(grid%element(elem)%vertex)%coord%y

degree(1+EDGES_PER_ELEMENT) = grid%element(elem)%degree
if (degree(1+EDGES_PER_ELEMENT) >= 3) then
   nred = degree(1+EDGES_PER_ELEMENT) - 2
else
   nred = 0
endif
do i=1,EDGES_PER_ELEMENT
   if (old_edge_deg(i) >= grid%edge(grid%element(elem)%edge(i))%degree) then
      degree(i) = 1 ! forces no red bases on this side
   elseif (old_edge_deg(i) == grid%edge(grid%element(elem)%edge(i))%degree-1) then
      degree(i) = grid%edge(grid%element(elem)%edge(i))%degree
      nred = nred + 1
   else
      ierr = PHAML_INTERNAL_ERROR
      call fatal("init_guess_p: edge degree changed by more than one")
      stop
   endif
   edge_type(i,:) = grid%edge_type(grid%element(elem)%edge(i),:)
   bmark(i) = grid%edge(grid%element(elem)%edge(i))%bmark
end do

! if nred is 0 then the element is degree 2 and all edges are greater than 2.
! There are no solution components to be initialized in this case.

if (nred == 0) return

n = nred*ss
allocate(mat(n,n),rs(n,nev),ipiv(n),elem_mat(ss,ss),elem_rs(ss,nev),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in error estimates")
   stop 
endif

if (nev > 1) then
   call elemental_matrix(grid, xvert, yvert, degree, edge_type, bmark, ss, &
                      0, .true., elem, rmin, "p", "r", mat, rs(:,1), &
                      loc_bconds_s=init_guess_dirich_bconds, extra_rhs=rs(:,2:))
else
   call elemental_matrix(grid, xvert, yvert, degree, edge_type, bmark, ss, &
                      0, .true., elem, rmin, "p", "r", mat, rs(:,1), &
                      loc_bconds_s=init_guess_dirich_bconds)
endif

! elemental 1x1 matrices over neighbor triangles by sides with red edge basis

neigh = get_neighbors(elem,grid)
bmark = huge(0)
edge_type(1,:) = INTERIOR
edge_type(2:3,:) = DIRICHLET
k = 1
do i=1,EDGES_PER_ELEMENT
   if (degree(i) == 1) cycle
   if (neigh(i) == BOUNDARY) then
      do j=1,ss
         mat(k+j-1,:) = 0
         mat(k+j-1,k+j-1) = 1
         rs(k+j-1,:) = 0.0_my_real
      end do
   else
      v1 = mod(i,3)+1
      v2 = mod(i+1,3)+1
      xvert(1) = grid%vertex(grid%element(elem)%vertex(v1))%coord%x
      yvert(1) = grid%vertex(grid%element(elem)%vertex(v1))%coord%y
      xvert(2) = grid%vertex(grid%element(elem)%vertex(v2))%coord%x
      yvert(2) = grid%vertex(grid%element(elem)%vertex(v2))%coord%y
      do j=1,3
         if (grid%element(neigh(i))%vertex(j) /= &
             grid%element(elem)%vertex(v1) .and. &
             grid%element(neigh(i))%vertex(j) /= &
             grid%element(elem)%vertex(v2)) then
            xvert(3) = grid%vertex(grid%element(neigh(i))%vertex(j))%coord%x
            yvert(3) = grid%vertex(grid%element(neigh(i))%vertex(j))%coord%y
            exit
         endif
      end do
      if (nev > 1) then
         call elemental_matrix(grid, xvert, yvert, (/1,1,degree(i),1/), &
                            edge_type, bmark, ss, 0, .true., &
                            neigh(i), rmin, "p", "r", elem_mat, elem_rs(:,1), &
                            loc_bconds_s=init_guess_dirich_bconds, &
                            extra_rhs=elem_rs(:,2:))
      else
         call elemental_matrix(grid, xvert, yvert, (/1,1,degree(i),1/), &
                            edge_type, bmark, ss, 0, .true., &
                            neigh(i), rmin, "p", "r", elem_mat, elem_rs(:,1), &
                            loc_bconds_s=init_guess_dirich_bconds)
      endif
      mat(k:k+ss-1,k:k+ss-1) = mat(k:k+ss-1,k:k+ss-1) + elem_mat
      rs(k:k+ss-1,:) = rs(k:k+ss-1,:) + elem_rs
   endif
   k = k+ss
end do
 
! solve the system

if (my_real == kind(0.0)) then
   call sgesv(n,nev,mat,n,ipiv,rs,n,info)
elseif (my_real == kind(0.0d0)) then
   call dgesv(n,nev,mat,n,ipiv,rs,n,info)
else
   ierr = USER_INPUT_ERROR
   call fatal("kind of real must be single or double to use LAPACK for error estimate")
   stop
endif

if (info /= 0) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("LAPACK returned error during init_guess_p",intlist=(/info/))
   stop
endif

! copy solution to grid%...%solution

do j=1,nev
   k = 1
   do i=1,EDGES_PER_ELEMENT
    if (degree(i) == 1) cycle
    do l=1,ss
      grid%edge(grid%element(elem)%edge(i))%solution(degree(i)-1,l,j) = &
         grid%edge(grid%element(elem)%edge(i))%solution(degree(i)-1,l,j) + rs(k+l-1,j)
    end do
    k = k+ss
   end do
   do i=1+((degree(4)-2)*(degree(4)-3))/2,((degree(4)-1)*(degree(4)-2))/2
    do l=1,ss
      grid%element(elem)%solution(i,l,j) = &
         grid%element(elem)%solution(i,l,j) + rs(k+l-1,j)
    end do
    k = k+ss
   end do
end do

! clean up memory

deallocate(mat,rs,ipiv,elem_mat,elem_rs,stat=astat)

end subroutine init_guess_p

!          ------------
subroutine init_guess_h(grid,elem,mate,child)
!          ------------

!----------------------------------------------------
! This routine sets the initial guess for the solution components associated
! with bases that have been added by the h refinement of element elem
! by solving a local Dirichlet problem using all the p-hierarchical bases
! over the children of elements elem and the mate of elem.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: elem,mate,child(4)
!----------------------------------------------------
! Local variables:

real(my_real) :: xvert(3), yvert(3), rmin
real(my_real), allocatable :: elem_mat(:,:),elem_rs(:,:),full_mat(:,:), &
                              full_rs(:,:)
integer :: degree(4), nbasis_elem, nbasis_full, elem_n, &
           full_n, astat, edge_type(3,grid%system_size), ss, i, j, k, info, &
           bmark(3), start_renum(17), nb, nev
integer, allocatable :: renum(:), ipiv(:)
!----------------------------------------------------
! Begin executable code

ss = grid%system_size
nev = max(1,grid%num_eval)

! number of p-hierarchical bases over the refinement of the parent and mate

if (mate == BOUNDARY) then
   nbasis_full = 4
else
   nbasis_full = 5
endif
do i=1,4
   if (i==3 .and. mate==BOUNDARY) exit
   do j=1,3
      if ((i==2.or.i==4) .and. j==1) cycle
      if (i>2 .and. j==2) cycle
      nbasis_full = nbasis_full + &
                    grid%edge(grid%element(child(i))%edge(j))%degree - 1
   end do
   nbasis_full = nbasis_full + &
         ((grid%element(child(i))%degree-1)*(grid%element(child(i))%degree-2))/2
end do

full_n = nbasis_full*ss

allocate(full_mat(full_n,full_n), full_rs(full_n,nev), ipiv(full_n), stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in error estimates")
   stop
endif

full_mat = 0.0_my_real
full_rs = 0.0_my_real

! beginning of renum for each segment of unknowns.  Segments are:
!
!                    1   elem vertex 1
!                  / | \
!                 /  |  \
!               10   6   12
!       elem   /  14 | 16  \   mate
!             /      |      \
!            3---8---5---9---4
!             \      |      /
!              \  15 | 17  /
!               11   7   13
!                 \  |  /
!                  \ | /
!                    2   elem vertex 2

start_renum(1) = 1
start_renum(2) = 2
start_renum(3) = 3
if (mate == BOUNDARY) then
   start_renum(4) = 3
else
   start_renum(4) = 4
endif
start_renum(5) = start_renum(4) + 1
start_renum(6) = start_renum(5) + 1
start_renum(7) = start_renum(6) + grid%edge(grid%element(child(1))%edge(2))%degree-1
start_renum(8) = start_renum(7) + grid%edge(grid%element(child(2))%edge(2))%degree-1
start_renum(9) = start_renum(8) + grid%edge(grid%element(child(1))%edge(1))%degree-1
if (mate == BOUNDARY) then
   start_renum(10) = start_renum(9)
else
   start_renum(10) = start_renum(9) + grid%edge(grid%element(child(3))%edge(1))%degree-1
endif
start_renum(11) = start_renum(10) + grid%edge(grid%element(child(1))%edge(3))%degree-1
start_renum(12) = start_renum(11) + grid%edge(grid%element(child(2))%edge(3))%degree-1
if (mate == BOUNDARY) then
   start_renum(13) = 0
else
   start_renum(13) = start_renum(12) + grid%edge(grid%element(child(3))%edge(3))%degree-1
endif
if (mate == BOUNDARY) then
   start_renum(14) = start_renum(12)
else
   start_renum(14) = start_renum(13) + grid%edge(grid%element(child(4))%edge(3))%degree-1
endif
start_renum(15) = start_renum(14) + ((grid%element(child(1))%degree-1)*(grid%element(child(1))%degree-2))/2
start_renum(16) = start_renum(15) + ((grid%element(child(2))%degree-1)*(grid%element(child(2))%degree-2))/2
if (mate == BOUNDARY) then
   start_renum(17) = 0
else
   start_renum(17) = start_renum(16) + ((grid%element(child(3))%degree-1)*(grid%element(child(3))%degree-2))/2
endif

! don't actually use rmin

rmin = 0.0_my_real

! for each of the children of elem and mate, compute the elemental matrix
! and right side, and assemble

! first child

degree = (/ grid%edge(grid%element(child(1))%edge(1))%degree, &
            grid%edge(grid%element(child(1))%edge(2))%degree, &
            grid%edge(grid%element(child(1))%edge(3))%degree, &
            grid%element(child(1))%degree /)
xvert = (/ grid%vertex(grid%element(child(1))%vertex(1))%coord%x, &
           grid%vertex(grid%element(child(1))%vertex(2))%coord%x, &
           grid%vertex(grid%element(child(1))%vertex(3))%coord%x /)
yvert = (/ grid%vertex(grid%element(child(1))%vertex(1))%coord%y, &
           grid%vertex(grid%element(child(1))%vertex(2))%coord%y, &
           grid%vertex(grid%element(child(1))%vertex(3))%coord%y /)

edge_type(1,:) = INTERIOR
bmark(1) = 0
if (mate == BOUNDARY) then
   edge_type(2,:) = grid%edge_type(grid%element(child(1))%edge(2),:)
   bmark(2) = grid%edge(grid%element(child(1))%edge(2))%bmark
else
   edge_type(2,:) = INTERIOR
   bmark(2) = 0
endif
edge_type(3,:) = DIRICHLET
bmark(3) = huge(0)

nbasis_elem = 3
do i=1,3
   nbasis_elem = nbasis_elem + degree(i)-1
end do
nbasis_elem = nbasis_elem + ((degree(4)-1)*(degree(4)-2))/2

elem_n = nbasis_elem*ss

allocate(elem_mat(elem_n,elem_n), elem_rs(elem_n,nev), renum(elem_n),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in error estimates")
   stop
endif

elem_mat = 0.0_my_real
elem_rs = 0.0_my_real

if (nev > 1) then
   call elemental_matrix(grid, xvert, yvert, degree, edge_type, bmark, ss, &
                         0, .false., elem, rmin, "p", "a", elem_mat, &
                         elem_rs(:,1), loc_bconds_s=init_guess_dirich_bconds, &
                         extra_rhs=elem_rs(:,2:))
else
   call elemental_matrix(grid, xvert, yvert, degree, edge_type, bmark, ss, &
                         0, .false., elem, rmin, "p", "a", elem_mat, &
                         elem_rs(:,1), loc_bconds_s=init_guess_dirich_bconds)
endif

renum(1:ss) = (/ (1+(start_renum(1)-1)*ss+j,j=0,ss-1) /)
renum(ss+1:2*ss) = (/ (1+(start_renum(3)-1)*ss+j,j=0,ss-1) /)
renum(2*ss+1:3*ss) = (/ (1+(start_renum(5)-1)*ss+j,j=0,ss-1) /)
i = 3*ss+1
nb = grid%edge(grid%element(child(1))%edge(1))%degree-1
renum(i:i+nb*ss-1) = (/ (1+(start_renum(8)-1)*ss+j,j=0,nb*ss-1) /)
i = i + nb*ss
nb = grid%edge(grid%element(child(1))%edge(2))%degree-1
renum(i:i+nb*ss-1) = (/ (1+(start_renum(6)-1)*ss+j,j=0,nb*ss-1) /)
i = i + nb*ss
nb = grid%edge(grid%element(child(1))%edge(3))%degree-1
renum(i:i+nb*ss-1) = (/ (1+(start_renum(10)-1)*ss+j,j=0,nb*ss-1) /)
i = i + nb*ss
nb = ((grid%element(child(1))%degree-1)*(grid%element(child(1))%degree-2))/2
renum(i:i+nb*ss-1) = (/ (1+(start_renum(14)-1)*ss+j,j=0,nb*ss-1) /)

do i=1,elem_n
   do j=1,elem_n
      full_mat(renum(i),renum(j)) = full_mat(renum(i),renum(j)) + elem_mat(i,j)
   end do
   full_rs(renum(i),:) = full_rs(renum(i),:) + elem_rs(i,:)
end do

! second child

degree = (/ grid%edge(grid%element(child(2))%edge(1))%degree, &
            grid%edge(grid%element(child(2))%edge(2))%degree, &
            grid%edge(grid%element(child(2))%edge(3))%degree, &
            grid%element(child(2))%degree /)
xvert = (/ grid%vertex(grid%element(child(2))%vertex(1))%coord%x, &
           grid%vertex(grid%element(child(2))%vertex(2))%coord%x, &
           grid%vertex(grid%element(child(2))%vertex(3))%coord%x /)
yvert = (/ grid%vertex(grid%element(child(2))%vertex(1))%coord%y, &
           grid%vertex(grid%element(child(2))%vertex(2))%coord%y, &
           grid%vertex(grid%element(child(2))%vertex(3))%coord%y /)

edge_type(1,:) = INTERIOR
bmark(1) = 0
if (mate == BOUNDARY) then
   edge_type(2,:) = grid%edge_type(grid%element(child(2))%edge(2),:)
   bmark(2) = grid%edge(grid%element(child(2))%edge(2))%bmark
else
   edge_type(2,:) = INTERIOR
   bmark(2) = 0
endif
edge_type(3,:) = DIRICHLET
bmark(3) = huge(0)

nbasis_elem = 3
do i=1,3
   nbasis_elem = nbasis_elem + degree(i)-1
end do
nbasis_elem = nbasis_elem + ((degree(4)-1)*(degree(4)-2))/2

elem_n = nbasis_elem*ss

if (elem_n /= size(renum)) then
   deallocate(elem_mat,elem_rs,renum,stat=astat)
   allocate(elem_mat(elem_n,elem_n), elem_rs(elem_n,nev), renum(elem_n), &
            stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in error estimates")
      stop
   endif
endif

elem_mat = 0.0_my_real
elem_rs = 0.0_my_real

if (nev > 1) then
   call elemental_matrix(grid, xvert, yvert, degree, edge_type, bmark, ss, &
                         0, .false., elem, rmin, "p", "a", elem_mat, &
                         elem_rs(:,1), loc_bconds_s=init_guess_dirich_bconds, &
                         extra_rhs=elem_rs(:,2:))
else
   call elemental_matrix(grid, xvert, yvert, degree, edge_type, bmark, ss, &
                         0, .false., elem, rmin, "p", "a", elem_mat, &
                         elem_rs(:,1), loc_bconds_s=init_guess_dirich_bconds)
endif

renum(1:ss) = (/ (1+(start_renum(2)-1)*ss+j,j=0,ss-1) /)
renum(ss+1:2*ss) = (/ (1+(start_renum(3)-1)*ss+j,j=0,ss-1) /)
renum(2*ss+1:3*ss) = (/ (1+(start_renum(5)-1)*ss+j,j=0,ss-1) /)
i = 3*ss+1
nb = grid%edge(grid%element(child(2))%edge(1))%degree-1
renum(i:i+nb*ss-1) = (/ (1+(start_renum(8)-1)*ss+j,j=0,nb*ss-1) /)
i = i + nb*ss
nb = grid%edge(grid%element(child(2))%edge(2))%degree-1
renum(i:i+nb*ss-1) = (/ (1+(start_renum(7)-1)*ss+j,j=0,nb*ss-1) /)
i = i + nb*ss
nb = grid%edge(grid%element(child(2))%edge(3))%degree-1
renum(i:i+nb*ss-1) = (/ (1+(start_renum(11)-1)*ss+j,j=0,nb*ss-1) /)
i = i + nb*ss
nb = ((grid%element(child(2))%degree-1)*(grid%element(child(2))%degree-2))/2
renum(i:i+nb*ss-1) = (/ (1+(start_renum(15)-1)*ss+j,j=0,nb*ss-1) /)

do i=1,elem_n
   do j=1,elem_n
      full_mat(renum(i),renum(j)) = full_mat(renum(i),renum(j)) + elem_mat(i,j)
   end do
   full_rs(renum(i),:) = full_rs(renum(i),:) + elem_rs(i,:)
end do

if (mate /= BOUNDARY) then

! third child

   degree = (/ grid%edge(grid%element(child(3))%edge(1))%degree, &
               grid%edge(grid%element(child(3))%edge(2))%degree, &
               grid%edge(grid%element(child(3))%edge(3))%degree, &
               grid%element(child(3))%degree /)
   xvert = (/ grid%vertex(grid%element(child(3))%vertex(1))%coord%x, &
              grid%vertex(grid%element(child(3))%vertex(2))%coord%x, &
              grid%vertex(grid%element(child(3))%vertex(3))%coord%x /)
   yvert = (/ grid%vertex(grid%element(child(3))%vertex(1))%coord%y, &
              grid%vertex(grid%element(child(3))%vertex(2))%coord%y, &
              grid%vertex(grid%element(child(3))%vertex(3))%coord%y /)

   edge_type(1,:) = INTERIOR
   bmark(1) = 0
   if (mate == BOUNDARY) then
      edge_type(2,:) = grid%edge_type(grid%element(child(3))%edge(2),:)
      bmark(2) = grid%edge(grid%element(child(3))%edge(2))%bmark
   else
      edge_type(2,:) = INTERIOR
      bmark(2) = 0
   endif
   edge_type(3,:) = DIRICHLET
   bmark(3) = huge(0)

   nbasis_elem = 3
   do i=1,3
      nbasis_elem = nbasis_elem + degree(i)-1
   end do
   nbasis_elem = nbasis_elem + ((degree(4)-1)*(degree(4)-2))/2

   elem_n = nbasis_elem*ss

   if (elem_n /= size(renum)) then
      deallocate(elem_mat,elem_rs,renum,stat=astat)
      allocate(elem_mat(elem_n,elem_n),elem_rs(elem_n,nev),renum(elem_n), &
               stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("memory allocation failed in error estimates")
         stop
      endif
   endif

   elem_mat = 0.0_my_real
   elem_rs = 0.0_my_real

   if (nev > 1) then
      call elemental_matrix(grid, xvert, yvert, degree, edge_type, bmark, ss, &
                            0, .false., mate, rmin, "p", "a", elem_mat, &
                          elem_rs(:,1), loc_bconds_s=init_guess_dirich_bconds, &
                            extra_rhs=elem_rs(:,2:))
   else
      call elemental_matrix(grid, xvert, yvert, degree, edge_type, bmark, ss, &
                            0, .false., mate, rmin, "p", "a", elem_mat, &
                            elem_rs(:,1), loc_bconds_s=init_guess_dirich_bconds)
   endif

   renum(1:ss) = (/ (1+(start_renum(1)-1)*ss+j,j=0,ss-1) /)
   renum(ss+1:2*ss) = (/ (1+(start_renum(4)-1)*ss+j,j=0,ss-1) /)
   renum(2*ss+1:3*ss) = (/ (1+(start_renum(5)-1)*ss+j,j=0,ss-1) /)
   i = 3*ss+1
   nb = grid%edge(grid%element(child(3))%edge(1))%degree-1
   renum(i:i+nb*ss-1) = (/ (1+(start_renum(9)-1)*ss+j,j=0,nb*ss-1) /)
   i = i + nb*ss
   nb = grid%edge(grid%element(child(3))%edge(2))%degree-1
   renum(i:i+nb*ss-1) = (/ (1+(start_renum(6)-1)*ss+j,j=0,nb*ss-1) /)
   i = i + nb*ss
   nb = grid%edge(grid%element(child(3))%edge(3))%degree-1
   renum(i:i+nb*ss-1) = (/ (1+(start_renum(12)-1)*ss+j,j=0,nb*ss-1) /)
   i = i + nb*ss
   nb = ((grid%element(child(3))%degree-1)*(grid%element(child(3))%degree-2))/2
   renum(i:i+nb*ss-1) = (/ (1+(start_renum(16)-1)*ss+j,j=0,nb*ss-1) /)

   do i=1,elem_n
      do j=1,elem_n
         full_mat(renum(i),renum(j)) = full_mat(renum(i),renum(j))+elem_mat(i,j)
      end do
      full_rs(renum(i),:) = full_rs(renum(i),:) + elem_rs(i,:)
   end do

! fourth child

   degree = (/ grid%edge(grid%element(child(4))%edge(1))%degree, &
               grid%edge(grid%element(child(4))%edge(2))%degree, &
               grid%edge(grid%element(child(4))%edge(3))%degree, &
               grid%element(child(4))%degree /)
   xvert = (/ grid%vertex(grid%element(child(4))%vertex(1))%coord%x, &
              grid%vertex(grid%element(child(4))%vertex(2))%coord%x, &
              grid%vertex(grid%element(child(4))%vertex(3))%coord%x /)
   yvert = (/ grid%vertex(grid%element(child(4))%vertex(1))%coord%y, &
              grid%vertex(grid%element(child(4))%vertex(2))%coord%y, &
              grid%vertex(grid%element(child(4))%vertex(3))%coord%y /)

   edge_type(1,:) = INTERIOR
   bmark(1) = 0
   if (mate == BOUNDARY) then
      edge_type(2,:) = grid%edge_type(grid%element(child(4))%edge(2),:)
      bmark(2) = grid%edge(grid%element(child(4))%edge(2))%bmark
   else
      edge_type(2,:) = INTERIOR
      bmark(2) = 0
   endif
   edge_type(3,:) = DIRICHLET
   bmark(3) = huge(0)

   nbasis_elem = 3
   do i=1,3
      nbasis_elem = nbasis_elem + degree(i)-1
   end do
   nbasis_elem = nbasis_elem + ((degree(4)-1)*(degree(4)-2))/2

   elem_n = nbasis_elem*ss

   if (elem_n /= size(renum)) then
      deallocate(elem_mat,elem_rs,renum,stat=astat)
      allocate(elem_mat(elem_n,elem_n),elem_rs(elem_n,nev),renum(elem_n), &
               stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("memory allocation failed in error estimates")
         stop
      endif
   endif

   elem_mat = 0.0_my_real
   elem_rs = 0.0_my_real

   if (nev > 1) then
      call elemental_matrix(grid, xvert, yvert, degree, edge_type, bmark, ss, &
                            0, .false., mate, rmin, "p", "a", elem_mat, &
                          elem_rs(:,1), loc_bconds_s=init_guess_dirich_bconds, &
                            extra_rhs=elem_rs(:,2:))
   else
      call elemental_matrix(grid, xvert, yvert, degree, edge_type, bmark, ss, &
                            0, .false., mate, rmin, "p", "a", elem_mat, &
                            elem_rs(:,1), loc_bconds_s=init_guess_dirich_bconds)
   endif

   renum(1:ss) = (/ (1+(start_renum(2)-1)*ss+j,j=0,ss-1) /)
   renum(ss+1:2*ss) = (/ (1+(start_renum(4)-1)*ss+j,j=0,ss-1) /)
   renum(2*ss+1:3*ss) = (/ (1+(start_renum(5)-1)*ss+j,j=0,ss-1) /)
   i = 3*ss+1
   nb = grid%edge(grid%element(child(4))%edge(1))%degree-1
   renum(i:i+nb*ss-1) = (/ (1+(start_renum(9)-1)*ss+j,j=0,nb*ss-1) /)
   i = i + nb*ss
   nb = grid%edge(grid%element(child(4))%edge(2))%degree-1
   renum(i:i+nb*ss-1) = (/ (1+(start_renum(7)-1)*ss+j,j=0,nb*ss-1) /)
   i = i + nb*ss
   nb = grid%edge(grid%element(child(4))%edge(3))%degree-1
   renum(i:i+nb*ss-1) = (/ (1+(start_renum(13)-1)*ss+j,j=0,nb*ss-1) /)
   i = i + nb*ss
   nb = ((grid%element(child(4))%degree-1)*(grid%element(child(4))%degree-2))/2
   renum(i:i+nb*ss-1) = (/ (1+(start_renum(17)-1)*ss+j,j=0,nb*ss-1) /)

   do i=1,elem_n
      do j=1,elem_n
         full_mat(renum(i),renum(j)) = full_mat(renum(i),renum(j))+elem_mat(i,j)
      end do
      full_rs(renum(i),:) = full_rs(renum(i),:) + elem_rs(i,:)
   end do

else ! mate is BOUNDARY

! replace rows corresponding to Dirichlet nodes with an identity row and
! the right side with the solution

! vertex 1

   do k=1,ss
      if (grid%vertex_type(grid%element(child(1))%vertex(1),k) == DIRICHLET) then
         full_rs((start_renum(1)-1)*ss+k,:) = &
            grid%vertex_solution(grid%element(child(1))%vertex(1),k,:)
         full_mat((start_renum(1)-1)*ss+k,:) = 0.0_my_real
         full_mat((start_renum(1)-1)*ss+k,(start_renum(1)-1)*ss+k) = 1.0_my_real
      endif
   end do

! vertex 2

   do k=1,ss
      if (grid%vertex_type(grid%element(child(2))%vertex(1),k) == DIRICHLET) then
         full_rs((start_renum(2)-1)*ss+k,:) = &
            grid%vertex_solution(grid%element(child(2))%vertex(1),k,:)
         full_mat((start_renum(2)-1)*ss+k,:) = 0.0_my_real
         full_mat((start_renum(2)-1)*ss+k,(start_renum(2)-1)*ss+k) = 1.0_my_real
      endif
   end do

! new vertex

   do k=1,ss
      if (grid%vertex_type(grid%element(child(1))%vertex(3),k) == DIRICHLET) then
         full_rs((start_renum(5)-1)*ss+k,:) = &
            grid%vertex_solution(grid%element(child(1))%vertex(3),k,:)
         full_mat((start_renum(5)-1)*ss+k,:) = 0.0_my_real
         full_mat((start_renum(5)-1)*ss+k,(start_renum(5)-1)*ss+k) = 1.0_my_real
      endif
   end do

! edge of first child

   do k=1,ss
      if (grid%edge_type(grid%element(child(1))%edge(2),k) == DIRICHLET) then
         do j=1,grid%edge(grid%element(child(1))%edge(2))%degree-1
            full_rs((start_renum(6)+j-2)*ss+k,:) = grid%edge(grid%element(child(1))%edge(2))%solution(j,k,:)
            full_mat((start_renum(6)+j-2)*ss+k,:) = 0.0_my_real
            full_mat((start_renum(6)+j-2)*ss+k,(start_renum(6)+j-2)*ss+k) = 1.0_my_real
         end do
      endif
   end do

! edge of second child

   do k=1,ss
      if (grid%edge_type(grid%element(child(2))%edge(2),k) == DIRICHLET) then
         do j=1,grid%edge(grid%element(child(2))%edge(2))%degree-1
            full_rs((start_renum(7)+j-2)*ss+k,:) = grid%edge(grid%element(child(2))%edge(2))%solution(j,k,:)
            full_mat((start_renum(7)+j-2)*ss+k,:) = 0.0_my_real
            full_mat((start_renum(7)+j-2)*ss+k,(start_renum(7)+j-2)*ss+k) = 1.0_my_real
         end do
      endif
   end do

endif

deallocate(elem_mat,elem_rs,renum,stat=astat)

! replace the equations on the boundary of parent U mate with Dirichlet

! vertex 1
do k=1,ss
   full_rs((start_renum(1)-1)*ss+k,:) = &
      grid%vertex_solution(grid%element(child(1))%vertex(1),k,:)
   full_mat((start_renum(1)-1)*ss+k,:) = 0.0_my_real
   full_mat((start_renum(1)-1)*ss+k,(start_renum(1)-1)*ss+k) = 1.0_my_real
end do
! vertex 2
do k=1,ss
   full_rs((start_renum(2)-1)*ss+k,:) = &
      grid%vertex_solution(grid%element(child(2))%vertex(1),k,:)
   full_mat((start_renum(2)-1)*ss+k,:) = 0.0_my_real
   full_mat((start_renum(2)-1)*ss+k,(start_renum(2)-1)*ss+k) = 1.0_my_real
end do
! vertex 3
do k=1,ss
   full_rs((start_renum(3)-1)*ss+k,:) = &
      grid%vertex_solution(grid%element(child(1))%vertex(2),k,:)
   full_mat((start_renum(3)-1)*ss+k,:) = 0.0_my_real
   full_mat((start_renum(3)-1)*ss+k,(start_renum(3)-1)*ss+k) = 1.0_my_real
end do
! vertex 4
if (mate /= BOUNDARY) then
   do k=1,ss
      full_rs((start_renum(4)-1)*ss+k,:) = &
         grid%vertex_solution(grid%element(child(3))%vertex(2),k,:)
      full_mat((start_renum(4)-1)*ss+k,:) = 0.0_my_real
      full_mat((start_renum(4)-1)*ss+k,(start_renum(4)-1)*ss+k) = 1.0_my_real
   end do
endif
! edge of child 1
do k=1,ss
   do j=1,grid%edge(grid%element(child(1))%edge(3))%degree-1
      full_rs((start_renum(10)+j-2)*ss+k,:) = grid%edge(grid%element(child(1))%edge(3))%solution(j,k,:)
      full_mat((start_renum(10)+j-2)*ss+k,:) = 0.0_my_real
      full_mat((start_renum(10)+j-2)*ss+k,(start_renum(10)+j-2)*ss+k) = 1.0_my_real
   end do
end do
! edge of child 2
do k=1,ss
   do j=1,grid%edge(grid%element(child(2))%edge(3))%degree-1
      full_rs((start_renum(11)+j-2)*ss+k,:) = grid%edge(grid%element(child(2))%edge(3))%solution(j,k,:)
      full_mat((start_renum(11)+j-2)*ss+k,:) = 0.0_my_real
      full_mat((start_renum(11)+j-2)*ss+k,(start_renum(11)+j-2)*ss+k) = 1.0_my_real
   end do
end do
if (mate /= BOUNDARY) then
! edge of child 3
   do k=1,ss
      do j=1,grid%edge(grid%element(child(3))%edge(3))%degree-1
         full_rs((start_renum(12)+j-2)*ss+k,:) = grid%edge(grid%element(child(3))%edge(3))%solution(j,k,:)
         full_mat((start_renum(12)+j-2)*ss+k,:) = 0.0_my_real
         full_mat((start_renum(12)+j-2)*ss+k,(start_renum(12)+j-2)*ss+k) = 1.0_my_real
      end do
   end do
! edge of child 4
   do k=1,ss
      do j=1,grid%edge(grid%element(child(4))%edge(3))%degree-1
         full_rs((start_renum(13)+j-2)*ss+k,:) = grid%edge(grid%element(child(4))%edge(3))%solution(j,k,:)
         full_mat((start_renum(13)+j-2)*ss+k,:) = 0.0_my_real
         full_mat((start_renum(13)+j-2)*ss+k,(start_renum(13)+j-2)*ss+k) = 1.0_my_real
      end do
   end do
endif

! solve the system

if (my_real == kind(0.0)) then
   call sgesv(full_n,nev,full_mat,full_n,ipiv,full_rs,full_n,info)
elseif (my_real == kind(0.0d0)) then
   call dgesv(full_n,nev,full_mat,full_n,ipiv,full_rs,full_n,info)
else
   ierr = USER_INPUT_ERROR
   call fatal("kind of real must be single or double to use LAPACK for error estimate")
   stop
endif

if (info /= 0) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("LAPACK returned error in init_guess_h",intlist=(/info/))
   stop
endif

! copy solution to grid solution

do k=1,nev
 do j=1,ss

   grid%vertex_solution(grid%element(child(1))%vertex(3),j,k) = &
      full_rs(j+ss*(start_renum(5)-1),k)
   if (mate /= BOUNDARY) then  ! just in case periodic b.c.
      grid%vertex_solution(grid%element(child(3))%vertex(3),j,k) = &
         full_rs(j+ss*(start_renum(5)-1),k)
   endif
   do i=1,grid%edge(grid%element(child(1))%edge(2))%degree-1
      grid%edge(grid%element(child(1))%edge(2))%solution(i,j,k) = &
         full_rs(j+ss*(start_renum(6)+i-2),k)
   end do
   do i=1,grid%edge(grid%element(child(2))%edge(2))%degree-1
      grid%edge(grid%element(child(2))%edge(2))%solution(i,j,k) = &
         full_rs(j+ss*(start_renum(7)+i-2),k)
   end do
   do i=1,grid%edge(grid%element(child(1))%edge(1))%degree-1
      grid%edge(grid%element(child(1))%edge(1))%solution(i,j,k) = &
         full_rs(j+ss*(start_renum(8)+i-2),k)
   end do
   if (mate /= BOUNDARY) then
      do i=1,grid%edge(grid%element(child(3))%edge(2))%degree-1
         grid%edge(grid%element(child(3))%edge(2))%solution(i,j,k) = &
            full_rs(j+ss*(start_renum(6)+i-2),k) ! just in case periodic b.c.
      end do
      do i=1,grid%edge(grid%element(child(4))%edge(2))%degree-1
         grid%edge(grid%element(child(4))%edge(2))%solution(i,j,k) = &
            full_rs(j+ss*(start_renum(7)+i-2),k) ! just in case periodic b.c.
      end do
      do i=1,grid%edge(grid%element(child(3))%edge(1))%degree-1
         grid%edge(grid%element(child(3))%edge(1))%solution(i,j,k)= &
            full_rs(j+ss*(start_renum(9)+i-2),k)
      end do
   endif
   do i=1,((grid%element(child(1))%degree-1)*(grid%element(child(1))%degree-2))/2
      grid%element(child(1))%solution(i,j,k) = &
         full_rs(j+ss*(start_renum(14)+i-2),k)
   end do
   do i=1,((grid%element(child(2))%degree-1)*(grid%element(child(2))%degree-2))/2
      grid%element(child(2))%solution(i,j,k) = &
         full_rs(j+ss*(start_renum(15)+i-2),k)
   end do
   if (mate /= BOUNDARY) then
      do i=1,((grid%element(child(3))%degree-1)*(grid%element(child(3))%degree-2))/2
         grid%element(child(3))%solution(i,j,k) = &
            full_rs(j+ss*(start_renum(16)+i-2),k)
      end do
      do i=1,((grid%element(child(4))%degree-1)*(grid%element(child(4))%degree-2))/2
         grid%element(child(4))%solution(i,j,k) = &
            full_rs(j+ss*(start_renum(17)+i-2),k)
      end do
   endif

 end do
end do

! clean up memory

deallocate(full_mat,full_rs,ipiv,stat=astat)

end subroutine init_guess_h

!          -------------------
subroutine hcoarsen_init_guess(grid,parent,mate,children)
!          -------------------

!----------------------------------------------------
! This routine sets the initial guess of the solution for the parents
! when a pair of triangles is unbisected
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: parent,mate,children(:)
!----------------------------------------------------
! Local variables:

real(my_real), allocatable :: solution(:,:), solution_c(:,:), &
                              oldsoln(:,:), oldsoln_c(:,:)
real(my_real) :: xvert(4),yvert(4)
integer :: i,j,k,degree,astat,isub,jsub
logical :: do_old, a1, a2, a3, a4
!----------------------------------------------------
! Begin executable code

! see if the old solution should also be set

if (grid%oldsoln_exists) then
   a1 = associated(grid%element(children(1))%oldsoln)
   a2 = associated(grid%element(children(2))%oldsoln)
   if (mate==BOUNDARY) then
      a3 = .false.
      a4 = .false.
   else
      a3 = associated(grid%element(children(3))%oldsoln)
      a4 = associated(grid%element(children(4))%oldsoln)
   endif
   do_old = a1 .or. a2 .or. a3 .or. a4
   if ((a1 .neqv. a2) .or. (a3 .neqv. a4)) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("In h coarsening, children are different w.r.t. existence of oldsoln.")
      stop
   endif
else
   do_old = .false.
endif

! set degree to be the maximum degree, and allocate space for the local
! copy of the solution coefficients

degree = 0
do i=1,4
  if (i==3 .and. mate==BOUNDARY) exit
  do j=1,EDGES_PER_ELEMENT
     degree = max(degree,grid%edge(grid%element(children(i))%edge(j))%degree)
  end do
end do

allocate(solution((degree+1)**2+degree**2,grid%nsoln), &
         solution_c((degree+1)**2,grid%nsoln),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in hcoarsen_initial_guess")
   return
endif

if (do_old) then
   allocate(oldsoln((degree+1)**2+degree**2,grid%nsoln), &
            oldsoln_c((degree+1)**2,grid%nsoln),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in hcoarsen_initial_guess")
      return
   endif
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
if (mate /= BOUNDARY) then
   solution(4,:) = reshape( &
                grid%vertex_solution(grid%element(children(3))%vertex(2),:,:), &
                (/grid%nsoln/))
endif
isub = 5
do i=1,degree-1
   if (grid%edge(grid%element(children(1))%edge(1))%degree >= i+1) then
      solution(isub,:) = reshape(grid%edge(grid%element(children(1))%edge(1))%solution(i,:,:), &
                                 (/grid%nsoln/))
   endif
   isub = isub+1
end do
solution(isub,:) = reshape( &
                grid%vertex_solution(grid%element(children(1))%vertex(3),:,:), &
                (/grid%nsoln/))
isub = isub+1
do i=1,degree-1
   if (grid%edge(grid%element(children(1))%edge(2))%degree >= i+1) then
      solution(isub,:) = reshape(grid%edge(grid%element(children(1))%edge(2))%solution(i,:,:), &
                                 (/grid%nsoln/))
   endif
   isub = isub + 1
end do
do i=1,degree-1
   if (grid%edge(grid%element(children(1))%edge(3))%degree >= i+1) then
      solution(isub,:) = reshape(grid%edge(grid%element(children(1))%edge(3))%solution(i,:,:), &
                                 (/grid%nsoln/))
   endif
   isub = isub + 1
end do
k = 1
do j=1,degree-2
   do i=1,j
      if (grid%element(children(1))%degree >= j+2) then
         solution(isub,:) = reshape(grid%element(children(1))%solution(k,:,:), &
                                    (/grid%nsoln/))
      endif
      isub = isub + 1
      k = k+1
   end do
end do
do i=1,degree-1
   if (grid%edge(grid%element(children(2))%edge(2))%degree >= i+1) then
      solution(isub,:) = reshape(grid%edge(grid%element(children(2))%edge(2))%solution(i,:,:), &
                                    (/grid%nsoln/))
   endif
   isub = isub + 1
end do
do i=1,degree-1
   if (grid%edge(grid%element(children(2))%edge(3))%degree >= i+1) then
      solution(isub,:) = reshape(grid%edge(grid%element(children(2))%edge(3))%solution(i,:,:), &
                                    (/grid%nsoln/))
   endif
   isub = isub + 1
end do
k = 1
do j=1,degree-2
   do i=1,j
      if (grid%element(children(2))%degree >= j+2) then
         solution(isub,:) = reshape(grid%element(children(2))%solution(k,:,:), &
                                    (/grid%nsoln/))
      endif
      isub = isub + 1
      k = k+1
   end do
end do
if (mate /= BOUNDARY) then
   do i=1,degree-1
      if (grid%edge(grid%element(children(3))%edge(1))%degree >= i+1) then
         solution(isub,:) = reshape(grid%edge(grid%element(children(3))%edge(1))%solution(i,:,:), &
                                    (/grid%nsoln/))
      endif
      isub = isub + 1
   end do
   do i=1,degree-1
      if (grid%edge(grid%element(children(3))%edge(3))%degree >= i+1) then
         solution(isub,:) = reshape(grid%edge(grid%element(children(3))%edge(3))%solution(i,:,:), &
                                    (/grid%nsoln/))
      endif
      isub = isub + 1
   end do
   k = 1
   do j=1,degree-2
      do i=1,j
         if (grid%element(children(3))%degree >= j+2) then
            solution(isub,:) = reshape(grid%element(children(3))%solution(k,:,:), &
                                 (/grid%nsoln/))
         endif
         isub = isub + 1
         k = k+1
      end do
   end do
   do i=1,degree-1
      if (grid%edge(grid%element(children(4))%edge(3))%degree >= i+1) then
         solution(isub,:) = reshape(grid%edge(grid%element(children(4))%edge(3))%solution(i,:,:), &
                                 (/grid%nsoln/))
      endif
      isub = isub + 1
   end do
   k = 1
   do j=1,degree-2
      do i=1,j
         if (grid%element(children(4))%degree >= j+2) then
            solution(isub,:) = reshape(grid%element(children(4))%solution(k,:,:), &
                                 (/grid%nsoln/))
         endif
         isub = isub + 1
         k = k+1
      end do
   end do
endif

if (do_old) then
   oldsoln = 0
   oldsoln(1,:) = reshape( &
                  grid%vertex_oldsoln(grid%element(children(1))%vertex(1),:,:),&
                  (/grid%nsoln/))
   oldsoln(2,:) = reshape( &
                  grid%vertex_oldsoln(grid%element(children(2))%vertex(1),:,:),&
                  (/grid%nsoln/))
   oldsoln(3,:) = reshape( &
                  grid%vertex_oldsoln(grid%element(children(1))%vertex(2),:,:),&
                  (/grid%nsoln/))
   if (mate /= BOUNDARY) then
      oldsoln(4,:) = reshape( &
                  grid%vertex_oldsoln(grid%element(children(3))%vertex(2),:,:),&
                  (/grid%nsoln/))
   endif
   isub = 5
   do i=1,degree-1
      if (grid%edge(grid%element(children(1))%edge(1))%degree >= i+1 .and. &
          associated(grid%edge(grid%element(children(1))%edge(1))%oldsoln)) then
         oldsoln(isub,:) = reshape( &
                   grid%edge(grid%element(children(1))%edge(1))%oldsoln(i,:,:),&
                   (/grid%nsoln/))
      endif
      isub = isub+1
   end do
   oldsoln(isub,:) = reshape( &
                  grid%vertex_oldsoln(grid%element(children(1))%vertex(3),:,:),&
                  (/grid%nsoln/))
   isub = isub+1
   do i=1,degree-1
      if (grid%edge(grid%element(children(1))%edge(2))%degree >= i+1 .and. &
          associated(grid%edge(grid%element(children(1))%edge(2))%oldsoln)) then
         oldsoln(isub,:) = reshape( &
                  grid%edge(grid%element(children(1))%edge(2))%oldsoln(i,:,:), &
                  (/grid%nsoln/))
      endif
      isub = isub + 1
   end do
   do i=1,degree-1
      if (grid%edge(grid%element(children(1))%edge(3))%degree >= i+1 .and. &
          associated(grid%edge(grid%element(children(1))%edge(3))%oldsoln)) then
         oldsoln(isub,:) = reshape( &
                   grid%edge(grid%element(children(1))%edge(3))%oldsoln(i,:,:),&
                   (/grid%nsoln/))
      endif
      isub = isub + 1
   end do
   k = 1
   do j=1,degree-2
      do i=1,j
         if (grid%element(children(1))%degree >= j+2 .and. &
             associated(grid%element(children(1))%oldsoln)) then
            oldsoln(isub,:) = reshape(grid%element(children(1))%oldsoln(k,:,:),&
                                      (/grid%nsoln/))
        
         endif
         isub = isub + 1
         k = k+1
      end do
   end do
   do i=1,degree-1
      if (grid%edge(grid%element(children(2))%edge(2))%degree >= i+1 .and. &
          associated(grid%edge(grid%element(children(2))%edge(2))%oldsoln)) then
         oldsoln(isub,:) = reshape( &
                  grid%edge(grid%element(children(2))%edge(2))%oldsoln(i,:,:), &
                  (/grid%nsoln/))
      endif
      isub = isub + 1
   end do
   do i=1,degree-1
      if (grid%edge(grid%element(children(2))%edge(3))%degree >= i+1 .and. &
          associated(grid%edge(grid%element(children(2))%edge(3))%oldsoln)) then
         oldsoln(isub,:) = reshape( &
                  grid%edge(grid%element(children(2))%edge(3))%oldsoln(i,:,:), &
                  (/grid%nsoln/))
      endif
      isub = isub + 1
   end do
   k = 1
   do j=1,degree-2
      do i=1,j
         if (grid%element(children(2))%degree >= j+2 .and. &
             associated(grid%element(children(2))%oldsoln)) then
            oldsoln(isub,:) = reshape(grid%element(children(2))%oldsoln(k,:,:),&
                                      (/grid%nsoln/))
         endif
         isub = isub + 1
         k = k+1
      end do
   end do
   if (mate /= BOUNDARY) then
      do i=1,degree-1
         if (grid%edge(grid%element(children(3))%edge(1))%degree >= i+1 .and. &
             associated(grid%edge(grid%element(children(3))%edge(1))%oldsoln)) then
            oldsoln(isub,:) = reshape( &
                   grid%edge(grid%element(children(3))%edge(1))%oldsoln(i,:,:),&
                   (/grid%nsoln/))
         endif
         isub = isub + 1
      end do
      do i=1,degree-1
         if (grid%edge(grid%element(children(3))%edge(3))%degree >= i+1 .and. &
             associated(grid%edge(grid%element(children(3))%edge(3))%oldsoln)) then
            oldsoln(isub,:) = reshape( &
                   grid%edge(grid%element(children(3))%edge(3))%oldsoln(i,:,:),&
                   (/grid%nsoln/))
         endif
         isub = isub + 1
      end do
      k = 1
      do j=1,degree-2
         do i=1,j
            if (grid%element(children(3))%degree >= j+2 .and. &
                associated(grid%element(children(3))%oldsoln)) then
               oldsoln(isub,:) = reshape( &
                                    grid%element(children(3))%oldsoln(k,:,:), &
                                    (/grid%nsoln/))
            endif
            isub = isub + 1
            k = k+1
         end do
      end do
      do i=1,degree-1
         if (grid%edge(grid%element(children(4))%edge(3))%degree >= i+1 .and. &
             associated(grid%edge(grid%element(children(4))%edge(3))%oldsoln)) then
            oldsoln(isub,:) = reshape( &
                   grid%edge(grid%element(children(4))%edge(3))%oldsoln(i,:,:),&
                   (/grid%nsoln/))
         endif
         isub = isub + 1
      end do
      k = 1
      do j=1,degree-2
         do i=1,j
            if (grid%element(children(4))%degree >= j+2 .and. &
                associated(grid%element(children(4))%oldsoln)) then
               oldsoln(isub,:) = reshape( &
                                 grid%element(children(4))%oldsoln(k,:,:), &
                                 (/grid%nsoln/))
            endif
            isub = isub + 1
            k = k+1
         end do
      end do
   endif
endif ! do_old

! set the outer vertices of the parents

xvert(1) = grid%vertex(grid%element(parent)%vertex(1))%coord%x
yvert(1) = grid%vertex(grid%element(parent)%vertex(1))%coord%y
xvert(2) = grid%vertex(grid%element(parent)%vertex(2))%coord%x
yvert(2) = grid%vertex(grid%element(parent)%vertex(2))%coord%y
xvert(3) = grid%vertex(grid%element(parent)%vertex(3))%coord%x
yvert(3) = grid%vertex(grid%element(parent)%vertex(3))%coord%y
if (mate == BOUNDARY) then
   xvert(4) = huge(0.0_my_real)
   yvert(4) = huge(0.0_my_real)
else
   xvert(4) = grid%vertex(grid%element(mate)%vertex(3))%coord%x
   yvert(4) = grid%vertex(grid%element(mate)%vertex(3))%coord%y
endif

! convert solution to nodal basis coefficients

call phier2nodal(solution,xvert,yvert,degree)
if (do_old) call phier2nodal(oldsoln,xvert,yvert,degree)

! extract solution for the parents, omiting those that are "red" nodes in
! the children.  See subroutines phier2nodal and nodal2phier for the order
! of the nodes

! vertices
solution_c(1:4,:) = solution(1:4,:)
if (mate == BOUNDARY) then
   isub = 4
else
   isub = 5
endif
! edge opposite vertex 1 of parent
jsub = 6 + 4*(degree-1) + ((degree-2)*(degree-1))/2
do i=1,degree-1
   solution_c(isub,:) = solution(jsub,:)
   isub = isub+1
   jsub = jsub+1
end do
! edge opposite vertex 2 of parent
jsub = 6 + 2*(degree-1)
do i=1,degree-1
   solution_c(isub,:) = solution(jsub,:)
   isub = isub+1
   jsub = jsub+1
end do
! common edge between parent and mate, half with vertex 1
jsub = degree+6
do i=1,(degree-1)/2
   solution_c(isub,:) = solution(jsub,:)
   isub = isub+1
   jsub = jsub+2
end do
! common edge between parent and mate, central vertex
if (2*(degree/2) == degree) then
   solution_c(isub,:) = solution(degree+4,:)
   isub = isub+1
endif
! common edge between parent and mate, half with vertex 2
jsub = 5 + 4*(degree-1) + ((degree-2)*(degree-1))/2
if (2*(degree/2) == degree) jsub = jsub-1
do i=1,(degree-1)/2
   solution_c(isub,:) = solution(jsub,:)
   jsub = jsub-2
   isub = isub+1
end do
! interior of parent
do j=1,degree-2
! from child 1
   jsub = 5 + 3*(degree-1) + degree-2 + j
   do i=1,(degree-j-1)/2
      solution_c(isub,:) = solution(jsub,:)
      jsub = jsub + 2*degree - 4*i - 3
      isub = isub+1
   end do
! from bisection edge
   if (2*((degree+j)/2) == degree+j) then
      jsub = 4+degree-j
      solution_c(isub,:) = solution(jsub,:)
      isub = isub+1
   endif
! from child 2
   jsub = 5 + 5*(degree-1) + (degree-2)*(degree-1) - (j*(j-1))/2
   if (2*((degree+j)/2) == degree+j) jsub = jsub - (j+1)
   do i=1,(degree-j-1)/2
      solution_c(isub,:) = solution(jsub,:)
      if (2*((degree+j)/2) == degree+j) then
         jsub = jsub - (2*j + 4*i + 1)
      else
         jsub = jsub - (2*j + 4*i - 1)
      endif
      isub = isub+1
   end do
end do
if (mate /= BOUNDARY) then
! edge opposite vertex 1 of mate
   jsub = 6 + 7*(degree-1) + 3*(((degree-2)*(degree-1))/2)
   do i=1,degree-1
      solution_c(isub,:) = solution(jsub,:)
      isub = isub+1
      jsub = jsub+1
   end do
! edge opposite vertex 2 of mate
   jsub = 6 + 6*(degree-1) + (degree-2)*(degree-1)
   do i=1,degree-1
      solution_c(isub,:) = solution(jsub,:)
      isub = isub+1
      jsub = jsub+1
   end do
! interior of mate
   do j=1,degree-2
! from child 3
      jsub = 5 + 7*(degree-1) + (degree-2)*(degree-1) + degree-2 + j
      do i=1,(degree-j-1)/2
         solution_c(isub,:) = solution(jsub,:)
         jsub = jsub + 2*degree - 4*i - 3
         isub = isub+1
      end do
! from bisection edge
      if (2*((degree+j)/2) == degree+j) then
         jsub = 5 + 6*(degree-1) + (degree-2)*(degree-1) + 1 - j
         solution_c(isub,:) = solution(jsub,:)
         isub = isub+1
      endif
! from child 4
      jsub = 5 + 8*(degree-1) + 2*(degree-2)*(degree-1) - (j*(j-1))/2
      if (2*((degree+j)/2) == degree+j) jsub = jsub - (j+1)
      do i=1,(degree-j-1)/2
         solution_c(isub,:) = solution(jsub,:)
         if (2*((degree+j)/2) == degree+j) then
            jsub = jsub - (2*j + 4*i + 1)
         else
            jsub = jsub - (2*j + 4*i - 1)
         endif
         isub = isub+1
      end do
   end do
endif

if (do_old) then
! vertices
   oldsoln_c(1:4,:) = oldsoln(1:4,:)
   if (mate == BOUNDARY) then
      isub = 4
   else
      isub = 5
   endif
! edge opposite vertex 1 of parent
   jsub = 6 + 4*(degree-1) + ((degree-2)*(degree-1))/2
   do i=1,degree-1
      oldsoln_c(isub,:) = oldsoln(jsub,:)
      isub = isub+1
      jsub = jsub+1
   end do
! edge opposite vertex 2 of parent
   jsub = 6 + 2*(degree-1)
   do i=1,degree-1
      oldsoln_c(isub,:) = oldsoln(jsub,:)
      isub = isub+1
      jsub = jsub+1
   end do
! common edge between parent and mate, half with vertex 1
   jsub = degree+6
   do i=1,(degree-1)/2
      oldsoln_c(isub,:) = oldsoln(jsub,:)
      isub = isub+1
      jsub = jsub+2
   end do
! common edge between parent and mate, central vertex
   if (2*(degree/2) == degree) then
      oldsoln_c(isub,:) = oldsoln(degree+4,:)
      isub = isub+1
   endif
! common edge between parent and mate, half with vertex 2
   jsub = 5 + 4*(degree-1) + ((degree-2)*(degree-1))/2
   if (2*(degree/2) == degree) jsub = jsub-1
   do i=1,(degree-1)/2
      oldsoln_c(isub,:) = oldsoln(jsub,:)
      jsub = jsub-2
      isub = isub+1
   end do
! interior of parent
   do j=1,degree-2
! from child 1
      jsub = 5 + 3*(degree-1) + degree-2 + j
      do i=1,(degree-j-1)/2
         oldsoln_c(isub,:) = oldsoln(jsub,:)
         jsub = jsub + 2*degree - 4*i - 3
         isub = isub+1
      end do
! from bisection edge
      if (2*((degree+j)/2) == degree+j) then
         jsub = 4+degree-j
         oldsoln_c(isub,:) = oldsoln(jsub,:)
         isub = isub+1
      endif
! from child 2
      jsub = 5 + 5*(degree-1) + (degree-2)*(degree-1) - (j*(j-1))/2
      if (2*((degree+j)/2) == degree+j) jsub = jsub - (j+1)
      do i=1,(degree-j-1)/2
         oldsoln_c(isub,:) = oldsoln(jsub,:)
         if (2*((degree+j)/2) == degree+j) then
            jsub = jsub - (2*j + 4*i + 1)
         else
            jsub = jsub - (2*j + 4*i - 1)
         endif
         isub = isub+1
      end do
   end do
   if (mate /= BOUNDARY) then
! edge opposite vertex 1 of mate
      jsub = 6 + 7*(degree-1) + 3*(((degree-2)*(degree-1))/2)
      do i=1,degree-1
         oldsoln_c(isub,:) = oldsoln(jsub,:)
         isub = isub+1
         jsub = jsub+1
      end do
! edge opposite vertex 2 of mate
      jsub = 6 + 6*(degree-1) + (degree-2)*(degree-1)
      do i=1,degree-1
         oldsoln_c(isub,:) = oldsoln(jsub,:)
         isub = isub+1
         jsub = jsub+1
      end do
! interior of mate
      do j=1,degree-2
! from child 3
         jsub = 5 + 7*(degree-1) + (degree-2)*(degree-1) + degree-2 + j
         do i=1,(degree-j-1)/2
            oldsoln_c(isub,:) = oldsoln(jsub,:)
            jsub = jsub + 2*degree - 4*i - 3
            isub = isub+1
         end do
! from bisection edge
         if (2*((degree+j)/2) == degree+j) then
            jsub = 5 + 6*(degree-1) + (degree-2)*(degree-1) + 1 - j
            oldsoln_c(isub,:) = oldsoln(jsub,:)
            isub = isub+1
         endif
! from child 4
         jsub = 5 + 8*(degree-1) + 2*(degree-2)*(degree-1) - (j*(j-1))/2
         if (2*((degree+j)/2) == degree+j) jsub = jsub - (j+1)
         do i=1,(degree-j-1)/2
            oldsoln_c(isub,:) = oldsoln(jsub,:)
            if (2*((degree+j)/2) == degree+j) then
               jsub = jsub - (2*j + 4*i + 1)
            else
               jsub = jsub - (2*j + 4*i - 1)
            endif
            isub = isub+1
         end do
      end do
   endif
endif ! do_old

! convert solution to p-hierarchical basis coefficients

call nodal2phier(solution_c,xvert,yvert,degree)
if (do_old) call nodal2phier(oldsoln_c,xvert,yvert,degree)

! copy solution to grid data structure

grid%vertex_solution(grid%element(parent)%vertex(1),:,:) = &
   reshape(solution_c(1,:),(/grid%system_size,max(1,grid%num_eval)/))
grid%vertex_solution(grid%element(parent)%vertex(2),:,:) = &
   reshape(solution_c(2,:),(/grid%system_size,max(1,grid%num_eval)/))
grid%vertex_solution(grid%element(parent)%vertex(3),:,:) = &
   reshape(solution_c(3,:),(/grid%system_size,max(1,grid%num_eval)/))
if (mate == boundary) then
   isub = 3
else
   grid%vertex_solution(grid%element(mate)%vertex(3),:,:) = &
      reshape(solution_c(4,:),(/grid%system_size,max(1,grid%num_eval)/))
   isub = 4
endif
do i=1,grid%edge(grid%element(parent)%edge(1))%degree-1
   grid%edge(grid%element(parent)%edge(1))%solution(i,:,:) = &
         reshape(solution_c(isub+i,:),(/grid%system_size,max(1,grid%num_eval)/))
end do
isub = isub + degree-1
do i=1,grid%edge(grid%element(parent)%edge(2))%degree-1
   grid%edge(grid%element(parent)%edge(2))%solution(i,:,:) = &
         reshape(solution_c(isub+i,:),(/grid%system_size,max(1,grid%num_eval)/))
end do
isub = isub + degree-1
do i=1,grid%edge(grid%element(parent)%edge(3))%degree-1
   grid%edge(grid%element(parent)%edge(3))%solution(i,:,:) = &
         reshape(solution_c(isub+i,:),(/grid%system_size,max(1,grid%num_eval)/))
end do
isub = isub + degree-1
do i=1,((grid%element(parent)%degree-2)*(grid%element(parent)%degree-1))/2
   grid%element(parent)%solution(i,:,:) = &
         reshape(solution_c(isub+i,:),(/grid%system_size,max(1,grid%num_eval)/))
end do
isub = isub + ((degree-2)*(degree-1))/2
if (mate /= BOUNDARY) then
   do i=1,grid%edge(grid%element(mate)%edge(1))%degree-1
      grid%edge(grid%element(mate)%edge(1))%solution(i,:,:) = &
         reshape(solution_c(isub+i,:),(/grid%system_size,max(1,grid%num_eval)/))
   end do
   isub = isub + degree-1
   do i=1,grid%edge(grid%element(mate)%edge(2))%degree-1
      grid%edge(grid%element(mate)%edge(2))%solution(i,:,:) = &
         reshape(solution_c(isub+i,:),(/grid%system_size,max(1,grid%num_eval)/))
   end do
   isub = isub + degree-1
   do i=1,((grid%element(mate)%degree-2)*(grid%element(mate)%degree-1))/2
      grid%element(mate)%solution(i,:,:) = &
         reshape(solution_c(isub+i,:),(/grid%system_size,max(1,grid%num_eval)/))
   end do
endif

if (do_old) then
   if (associated(grid%edge(grid%element(children(1))%edge(1))%oldsoln)) &
      deallocate(grid%edge(grid%element(children(1))%edge(1))%oldsoln)
   if (associated(grid%edge(grid%element(children(1))%edge(2))%oldsoln)) &
      deallocate(grid%edge(grid%element(children(1))%edge(2))%oldsoln)
   if (associated(grid%edge(grid%element(children(2))%edge(2))%oldsoln)) &
      deallocate(grid%edge(grid%element(children(2))%edge(2))%oldsoln)
   if (mate /= BOUNDARY) then
      if (associated(grid%edge(grid%element(children(3))%edge(1))%oldsoln)) &
         deallocate(grid%edge(grid%element(children(3))%edge(1))%oldsoln)
   endif
   if (associated(grid%element(children(1))%oldsoln)) &
      deallocate(grid%element(children(1))%oldsoln)
   if (associated(grid%element(children(2))%oldsoln)) &
      deallocate(grid%element(children(2))%oldsoln)
   if (mate /= BOUNDARY) then
      if (associated(grid%element(children(3))%oldsoln)) &
         deallocate(grid%element(children(3))%oldsoln)
      if (associated(grid%element(children(4))%oldsoln)) &
         deallocate(grid%element(children(4))%oldsoln)
   endif
   allocate(grid%edge(grid%element(parent)%edge(3))%oldsoln( &
      size(grid%edge(grid%element(parent)%edge(3))%solution,dim=1), &
      size(grid%edge(grid%element(parent)%edge(3))%solution,dim=2), &
      size(grid%edge(grid%element(parent)%edge(3))%solution,dim=3)))
   if (a1) then
      allocate(grid%element(parent)%oldsoln( &
         size(grid%element(parent)%solution,dim=1), &
         size(grid%element(parent)%solution,dim=2), &
         size(grid%element(parent)%solution,dim=3)))
   endif
   if (mate /= BOUNDARY .and. a3) then
      allocate(grid%element(mate)%oldsoln( &
         size(grid%element(mate)%solution,dim=1), &
         size(grid%element(mate)%solution,dim=2), &
         size(grid%element(mate)%solution,dim=3)))
   endif

   grid%vertex_oldsoln(grid%element(parent)%vertex(1),:,:) = &
              reshape(oldsoln_c(1,:),(/grid%system_size,max(1,grid%num_eval)/))
   grid%vertex_oldsoln(grid%element(parent)%vertex(2),:,:) = &
              reshape(oldsoln_c(2,:),(/grid%system_size,max(1,grid%num_eval)/))
   grid%vertex_oldsoln(grid%element(parent)%vertex(3),:,:) = &
              reshape(oldsoln_c(3,:),(/grid%system_size,max(1,grid%num_eval)/))
   if (mate == boundary) then
      isub = 3
   else
      grid%vertex_oldsoln(grid%element(mate)%vertex(3),:,:) = &
              reshape(oldsoln_c(4,:),(/grid%system_size,max(1,grid%num_eval)/))
      isub = 4
   endif
   if (a1) then
      do i=1,grid%edge(grid%element(parent)%edge(1))%degree-1
         grid%edge(grid%element(parent)%edge(1))%oldsoln(i,:,:) = &
          reshape(oldsoln_c(isub+i,:),(/grid%system_size,max(1,grid%num_eval)/))
      end do
   endif
   isub = isub + degree-1
   if (a1) then
      do i=1,grid%edge(grid%element(parent)%edge(2))%degree-1
         grid%edge(grid%element(parent)%edge(2))%oldsoln(i,:,:) = &
          reshape(oldsoln_c(isub+i,:),(/grid%system_size,max(1,grid%num_eval)/))
      end do
   endif
   isub = isub + degree-1
   do i=1,grid%edge(grid%element(parent)%edge(3))%degree-1
      grid%edge(grid%element(parent)%edge(3))%oldsoln(i,:,:) = &
          reshape(oldsoln_c(isub+i,:),(/grid%system_size,max(1,grid%num_eval)/))
   end do
   isub = isub + degree-1
   if (a1) then
      do i=1,((grid%element(parent)%degree-2)*(grid%element(parent)%degree-1))/2
         grid%element(parent)%oldsoln(i,:,:) = &
          reshape(oldsoln_c(isub+i,:),(/grid%system_size,max(1,grid%num_eval)/))
      end do
   endif
   isub = isub + ((degree-2)*(degree-1))/2
   if (mate /= BOUNDARY .and. a3) then
      do i=1,grid%edge(grid%element(mate)%edge(1))%degree-1
         grid%edge(grid%element(mate)%edge(1))%oldsoln(i,:,:) = &
          reshape(oldsoln_c(isub+i,:),(/grid%system_size,max(1,grid%num_eval)/))
      end do
      isub = isub + degree-1
      do i=1,grid%edge(grid%element(mate)%edge(2))%degree-1
         grid%edge(grid%element(mate)%edge(2))%oldsoln(i,:,:) = &
          reshape(oldsoln_c(isub+i,:),(/grid%system_size,max(1,grid%num_eval)/))
      end do
      isub = isub + degree-1
      do i=1,((grid%element(mate)%degree-2)*(grid%element(mate)%degree-1))/2
         grid%element(mate)%oldsoln(i,:,:) = &
          reshape(oldsoln_c(isub+i,:),(/grid%system_size,max(1,grid%num_eval)/))
      end do
   endif
endif ! do_old

deallocate(solution,solution_c,stat=astat)
if (do_old) deallocate(oldsoln,oldsoln_c,stat=astat)

end subroutine hcoarsen_init_guess

!          ------------------------
subroutine init_guess_dirich_bconds(x,y,bmark,itype,c,rs)
!          ------------------------

!----------------------------------------------------
! This routine returns boundary conditions for the initial guess local Dirichlet
! problem.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x,y
integer, intent(in) :: bmark
integer, intent(out) :: itype(:)
real(my_real), intent(out) :: c(:,:),rs(:)
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

! bmark = huge(0) indicates the special Dirichlet boundary of the local problem.
! Don't need to return the solution; it gets copied from grid%...%solution

if (bmark == huge(0)) then

   itype = DIRICHLET
   c = 0.0_my_real
   rs = 0.0_my_real

! otherwise evaluate the given boundary conditions

else

   call bconds(x,y,bmark,itype,c,rs)

endif

end subroutine init_guess_dirich_bconds

end module refine_elements
