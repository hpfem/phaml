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

module refine_uniform_mod

!----------------------------------------------------
! This module contains routines for uniform refinement of the grid.
!----------------------------------------------------

!----------------------------------------------------
! Other modules used:

use global
use gridtype_mod
use refine_elements
use grid_util
use hash_mod
use omp_lib
!----------------------------------------------------

implicit none
private
public refine_uniform_p, refine_uniform_h

!----------------------------------------------------
! Non-module procedures used are:

!----------------------------------------------------
! The following defined constants are defined:

!----------------------------------------------------
! The following variables are defined:

!----------------------------------------------------

contains

!          ----------------
subroutine refine_uniform_p(grid,refine_control)
!          ----------------

!----------------------------------------------------
! This routine performs a uniform p refinement of the grid
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(refine_options), intent(in) :: refine_control
!----------------------------------------------------
! Local variables:

integer :: element_list(size(grid%element)), edge_list(size(grid%edge)), &
           nelem, nedge, i, delta_dof, delta_dof_own, total_delta_dof, &
           total_delta_dof_own, errcode, old_edge_deg(EDGES_PER_ELEMENT)
!----------------------------------------------------
! Begin executable code

! make a list of all owned leaves

call list_elements(grid,element_list,nelem,level=0,own=.true.,leaf=.true.)

! refine all elements on the list

! refine the interior of the elements

total_delta_dof = 0
total_delta_dof_own = 0

!$omp parallel &
!$omp  default(shared) &
!$omp  private(i,errcode,delta_dof,delta_dof_own) &
!$omp  reduction(+ : total_delta_dof, total_delta_dof_own)

!$omp do
do i=1,nelem
   call p_refine_element_interior(grid,refine_control,element_list(i), &
                                  errcode,delta_dof,delta_dof_own)
   total_delta_dof = total_delta_dof + delta_dof
   total_delta_dof_own = total_delta_dof_own + delta_dof_own
end do
!$omp end do

! do things that must be done OpenMP-sequentially after the parallel loop
!$omp single

grid%dof = grid%dof + total_delta_dof
grid%dof_own = grid%dof_own + total_delta_dof_own

! make a list of all the edges that no longer satisfy the edge rule

call list_edges_without_rule(grid,refine_control,element_list,nelem,edge_list, &
                             nedge)

! enforce the edge rule on the listed edges

total_delta_dof = 0
total_delta_dof_own = 0

!$omp end single
!$omp do
do i=1,nedge
   call enforce_edge_rule(grid,refine_control,edge_list(i),delta_dof, &
                          delta_dof_own)
   total_delta_dof = total_delta_dof + delta_dof
   total_delta_dof_own = total_delta_dof_own + delta_dof_own
end do
!$omp end do

!$omp end parallel

! do things that must be done OpenMP-sequentially after the parallel loop

grid%dof = grid%dof + total_delta_dof
grid%dof_own = grid%dof_own + total_delta_dof_own

! For uniform refinement we don't need to set an initial guess, other than 0,
! for the new solution components, because we don't need a new error indicator
! to determine if the element should be refined again.
!
! TEMP NOTE: the adaptive refinement version of this will need to maintain elist
!       and conditionally compute error indicators, but uniform doesn't need it.

end subroutine refine_uniform_p

!          ----------------
subroutine refine_uniform_h(grid,refine_control)
!          ----------------

!----------------------------------------------------
! This routine performs a uniform h refinement of the grid
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(refine_options), intent(in) :: refine_control
!----------------------------------------------------
! Local variables:

character(len=1), pointer :: reftype(:)
integer :: lev, elem, errcode, mate, parent, i, nelem, &
           element_list(size(grid%element)), vert_lid(2,size(grid%element)), &
           edge_lid(6,size(grid%element)), elem_lid(4,size(grid%element)), &
           one_vert_lid(2),one_edge_lid(6),one_elem_lid(4), &
           delta_dof, total_delta_dof, delta_dof_own, total_delta_dof_own, &
           delta_nelem, total_delta_nelem, delta_nelem_leaf, &
           total_delta_nelem_leaf, delta_nelem_leaf_own, &
           total_delta_nelem_leaf_own, delta_nedge, total_delta_nedge, &
           delta_nedge_own, total_delta_nedge_own, delta_nvert, &
           total_delta_nvert, delta_nvert_own, total_delta_nvert_own, &
           max_nlev, max_max_nlev

!----------------------------------------------------
! Begin executable code

! mark elements that are owned leaves in the input grid, or are needed for
! compatibility, for refinement.  Only mark one in each compatibly divisible
! pair.

allocate(reftype(size(grid%element)))
reftype = "n"
do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (grid%element(elem)%isleaf .and. grid%element(elem)%iown) then
         if (grid%element(elem)%mate == BOUNDARY) then
            mate = BOUNDARY
         else
            mate = hash_decode_key(grid%element(elem)%mate,grid%elem_hash)
         endif
         if (mate == BOUNDARY) then
            reftype(elem) = "h"
         elseif (mate /= HASH_NOT_FOUND) then
            reftype(elem) = "h"
            reftype(mate) = "n"
         endif
         if (mate == HASH_NOT_FOUND) then
            reftype(elem) = "h"
            parent = elem
! traverse the compatibility chain
            do
               parent = hash_decode_key(grid%element(parent)%mate/2,grid%elem_hash)
               if (grid%element(parent)%mate == BOUNDARY) then
                  mate = BOUNDARY
               else
                  mate = hash_decode_key(grid%element(parent)%mate,grid%elem_hash)
               endif
               if (mate == BOUNDARY) then
                  reftype(parent) = "h"
                  exit
               elseif (mate /= HASH_NOT_FOUND) then
                  reftype(parent) = "h"
                  reftype(mate) = "n"
                  exit
               endif
               if (mate == HASH_NOT_FOUND) then
                  reftype(parent) = "h"
               endif
            end do
         endif
      endif
      elem = grid%element(elem)%next
   end do
end do

! for each level, starting with the coarsest level ...

do lev=1,grid%nlev

! create a list of elements to be refined

   nelem = 0
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (reftype(elem) == "h") then
         nelem = nelem + 1
         element_list(nelem) = elem
      endif
      elem = grid%element(elem)%next
   end do

! if the list is empty, more on to the next level

   if (nelem == 0) cycle

! get lids for the children, and other things that must be done
! OpenMP-sequentially before the parallel loop

   call before_h_refine(grid,element_list,nelem,reftype=reftype, &
                        vert_lid=vert_lid,edge_lid=edge_lid,elem_lid=elem_lid)

   total_delta_dof = 0
   total_delta_dof_own = 0
   total_delta_nelem = 0
   total_delta_nelem_leaf = 0
   total_delta_nelem_leaf_own = 0
   total_delta_nedge = 0
   total_delta_nedge_own = 0
   total_delta_nvert = 0
   total_delta_nvert_own = 0
   max_max_nlev = 0

!$omp parallel do &
!$omp  default(shared) &
!$omp  private(i,errcode,one_elem_lid,one_edge_lid,one_vert_lid, &
!$omp    delta_dof, delta_dof_own, delta_nelem, delta_nelem_leaf, &
!$omp    delta_nelem_leaf_own, delta_nedge, delta_nedge_own, &
!$omp    delta_nvert, delta_nvert_own, max_nlev) &
!$omp  reduction(+ : total_delta_dof, total_delta_dof_own, total_delta_nelem, &
!$omp    total_delta_nelem_leaf, total_delta_nelem_leaf_own, total_delta_nedge,&
!$omp    total_delta_nedge_own, total_delta_nvert, total_delta_nvert_own) &
!$omp  reduction(max : max_max_nlev)

   do i=1,nelem
! TEMP I should be able to pass these array sections, but ifort gets SIGSEGV
      one_elem_lid = elem_lid(:,i)
      one_edge_lid = edge_lid(:,i)
      one_vert_lid = vert_lid(:,i)
      call bisect_triangle_pair(grid,element_list(i),errcode,refine_control, &
                                reftype=reftype,vert_child_lid=one_vert_lid, &
                   edge_child_lid=one_edge_lid,elem_child_lid=one_elem_lid, &
                   delta_dof=delta_dof,delta_dof_own=delta_dof_own, &
                   delta_nelem=delta_nelem,delta_nelem_leaf=delta_nelem_leaf, &
                   delta_nelem_leaf_own=delta_nelem_leaf_own, &
                   delta_nedge=delta_nedge,delta_nedge_own=delta_nedge_own, &
                   delta_nvert=delta_nvert,delta_nvert_own=delta_nvert_own, &
                   max_nlev=max_nlev)

      total_delta_dof = total_delta_dof + delta_dof
      total_delta_dof_own = total_delta_dof_own + delta_dof_own
      total_delta_nelem = total_delta_nelem + delta_nelem
      total_delta_nelem_leaf = total_delta_nelem_leaf + delta_nelem_leaf
      total_delta_nelem_leaf_own = total_delta_nelem_leaf_own + delta_nelem_leaf_own
      total_delta_nedge = total_delta_nedge + delta_nedge
      total_delta_nedge_own = total_delta_nedge_own + delta_nedge_own
      total_delta_nvert = total_delta_nvert + delta_nvert
      total_delta_nvert_own = total_delta_nvert_own + delta_nvert_own
      max_max_nlev = max(max_max_nlev,max_nlev)

   end do
!$omp end parallel do

! do things that must be done OpenMP-sequentially after the parallel loop

   grid%dof = grid%dof + total_delta_dof
   grid%dof_own = grid%dof_own + total_delta_dof_own
   grid%nelem = grid%nelem + total_delta_nelem
   grid%nelem_leaf = grid%nelem_leaf + total_delta_nelem_leaf
   grid%nelem_leaf_own = grid%nelem_leaf_own + total_delta_nelem_leaf_own
   grid%nedge = grid%nedge + total_delta_nedge
   grid%nedge_own = grid%nedge_own + total_delta_nedge_own
   grid%nvert = grid%nvert + total_delta_nvert
   grid%nvert_own = grid%nvert_own + total_delta_nvert_own
   grid%nlev = max(grid%nlev,max_max_nlev)

   call after_h_refine(grid,element_list,nelem,vert_lid,edge_lid,elem_lid)

end do

end subroutine refine_uniform_h

end module refine_uniform_mod
