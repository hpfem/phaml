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

module nlp_vars
!----------------------------------------------------
! This module contains variables that need to be accessed in routines for
! the NLP solver in the NLP hp-adaptive strategy
!----------------------------------------------------
use global
use gridtype_mod
implicit none
type(grid_type), pointer :: hold_grid
type(refine_options), pointer :: hold_refcont
integer, allocatable :: nlp_elem_list(:), nlp_inverse_elem(:)
integer :: nlp_nelem
real(my_real), allocatable :: nlp_m(:)
real(my_real) :: nlp_tau
end module nlp_vars

module grid_mod

!----------------------------------------------------
! This module contains subroutines to manipulate the grid.
!
! communication tags in this module are of the form 4xx
!----------------------------------------------------

!----------------------------------------------------
! Other modules used:

use global
use message_passing
use hash_mod
use gridtype_mod
use stopwatch
use sort_mod
use make_linsys
use error_estimators
use evaluate
use basis_functions
use sysdep
use quadrature_rules
use linsys_io
use linsystype_mod
use linear_system
use nlp_vars
!----------------------------------------------------

implicit none
private
public allocate_grid, deallocate_grid, realloc_solution, refine, &
       reconcile, delete_unowned_elements, create_element, set_weights, &
       init_grid, reset_dirich_exact, p_refine_elem, count_memory, &
       ext_refsoln_errest_edge, ext_refsoln_errest_elem, &
       nlp_dof, nlp_ddof, nlp_errest, nlp_derrest

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

   subroutine boundary_point(ipiece,s,x,y)
   use global
   integer, intent(in) :: ipiece
   real(my_real), intent(in) :: s
   real(my_real), intent(out) :: x,y
   end subroutine boundary_point

   function boundary_npiece(hole)
   integer, intent(in) :: hole
   integer :: boundary_npiece
   end function boundary_npiece

   subroutine boundary_param(start,finish)
   use global
   real(my_real), intent(out) :: start(:), finish(:)
   end subroutine boundary_param

   function regularity(x,y)
   use global
   real(my_real), intent(in) :: x(3),y(3)
   real(my_real) :: regularity
   end function regularity
end interface

!----------------------------------------------------
! The following variables are defined:

real(my_real) :: n3pee(3,1), npq(3)

!----------------------------------------------------
! The following defined constants are defined:

! RESTRICTION no more than 16 triangles share a vertex in the triangle data
integer, parameter :: MAX_TD_VERT_NEIGH = 16
!----------------------------------------------------

contains

!---------------------------------------------------------------------
!  INITIALIZATION AND FINALIZATION ROUTINES
!---------------------------------------------------------------------

!          -------------
subroutine allocate_grid(grid,nvert,nev,type,degree)
!          -------------

!----------------------------------------------------
! This routine performs initial allocation of memory in the grid data structure
! and initializes it with an empty grid
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: nvert, nev, type, degree

!----------------------------------------------------
! Local variables:

integer :: i, astat1, astat2, astat3, dstat, npiece, nevdim

!----------------------------------------------------
! Begin executable code

nevdim = max(1,nev)

! grid data structure

allocate(grid%element(4*nvert),grid%edge(4*nvert),stat=astat1)
allocate(grid%vertex(nvert),stat=astat2)
allocate(grid%errest_energy(nevdim), &
         grid%errest_Linf(grid%system_size*nevdim), &
         grid%errest_L2(grid%system_size*nevdim), &
         grid%errest_eigenvalue(nevdim),stat=astat3)
if (astat1 /= 0 .or. astat2 /= 0 .or. astat3 /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in allocate_grid")
   deallocate(grid%element,grid%edge,grid%vertex,grid%errest_energy, &
              grid%errest_Linf, grid%errest_L2,grid%errest_eigenvalue, &
              stat=dstat)
   return
endif
if (type == EIGENVALUE) then
   allocate(grid%eigenvalue(nev),grid%eigenprob_l2_resid(nev), &
            grid%eigenprob_variance(nev),stat=astat1)
   if (astat1 /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in allocate_grid")
      deallocate(grid%element,grid%edge,grid%vertex,grid%errest_energy, &
                 grid%errest_Linf, grid%errest_L2,grid%eigenvalue, &
                 grid%eigenprob_l2_resid,grid%eigenprob_variance, &
                 grid%errest_eigenvalue, stat=dstat)
      return
   endif
   grid%eigenvalue = 0.0_my_real
   grid%eigenprob_l2_resid = 0.0_my_real
   grid%eigenprob_variance = 0.0_my_real
else
   nullify(grid%eigenvalue,grid%eigenprob_l2_resid,grid%eigenprob_variance)
endif
nullify(grid%initial_neighbor)
call hash_table_init(grid%elem_hash,size(grid%element))
if (ierr == ALLOC_FAILED) then
   deallocate(grid%element,grid%edge,grid%vertex,grid%errest_energy, &
              grid%errest_Linf,grid%errest_L2,grid%eigenvalue, &
              grid%eigenprob_l2_resid,grid%eigenprob_variance, &
              grid%errest_eigenvalue,stat=dstat)
   if (associated(grid%eigenvalue)) then
      deallocate(grid%eigenvalue,grid%eigenprob_l2_resid, &
                 grid%eigenprob_variance,stat=dstat)
   endif
   return
endif
call hash_table_init(grid%edge_hash,size(grid%edge))
if (ierr == ALLOC_FAILED) then
   call hash_table_destroy(grid%elem_hash)
   deallocate(grid%element,grid%edge,grid%vertex,grid%errest_energy, &
              grid%errest_Linf,grid%errest_L2,grid%eigenvalue, &
              grid%eigenprob_l2_resid,grid%eigenprob_variance, &
              grid%errest_eigenvalue,stat=dstat)
   if (associated(grid%eigenvalue)) then
      deallocate(grid%eigenvalue,grid%eigenprob_l2_resid, &
                 grid%eigenprob_variance,stat=dstat)
   endif
   return
endif
call hash_table_init(grid%vert_hash,size(grid%vertex))
if (ierr == ALLOC_FAILED) then
   call hash_table_destroy(grid%edge_hash)
   call hash_table_destroy(grid%elem_hash)
   deallocate(grid%element,grid%edge,grid%vertex,grid%errest_energy, &
              grid%errest_Linf,grid%errest_L2,grid%eigenvalue, &
              grid%eigenprob_l2_resid,grid%eigenprob_variance, &
              grid%errest_eigenvalue,stat=dstat)
   if (associated(grid%eigenvalue)) then
      deallocate(grid%eigenvalue,grid%eigenprob_l2_resid, &
                 grid%eigenprob_variance,stat=dstat)
   endif
   return
endif
grid%next_free_elem = 1
grid%next_free_edge = 1
grid%next_free_vert = 1
allocate(grid%head_level_elem(8),stat=astat1)
allocate(grid%head_level_vert(8),stat=astat2)
if (astat1 /= 0 .or. astat2 /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in allocate_grid")
   deallocate(grid%head_level_elem,grid%head_level_vert,stat=dstat)
   call hash_table_destroy(grid%elem_hash)
   call hash_table_destroy(grid%edge_hash)
   call hash_table_destroy(grid%vert_hash)
   deallocate(grid%element,grid%edge,grid%vertex,grid%errest_energy, &
              grid%errest_Linf,grid%errest_L2,grid%eigenvalue, &
              grid%eigenprob_l2_resid,grid%eigenprob_variance, &
              grid%errest_eigenvalue,stat=dstat)
   if (associated(grid%eigenvalue)) then
      deallocate(grid%eigenvalue,grid%eigenprob_l2_resid, &
                 grid%eigenprob_variance,stat=dstat)
   endif
   return
endif
grid%head_level_elem = END_OF_LIST
grid%head_level_vert = END_OF_LIST
grid%nelem = 0
grid%nelem_leaf = 0
grid%nelem_leaf_own = 0
grid%nedge = 0
grid%nedge_own = 0
grid%nvert = 0
grid%nvert_own = 0
grid%dof = 0
grid%dof_own = 0
grid%nlev = 0

! elements

grid%element%degree = 1
do i=1,size(grid%element)
   grid%element(i)%previous = i-1
   grid%element(i)%next = i+1
end do
grid%element(1)%previous = END_OF_LIST
grid%element(size(grid%element))%next = END_OF_LIST
do i=1,size(grid%element)
   nullify(grid%element(i)%solution,grid%element(i)%exact, &
           grid%element(i)%oldsoln)
end do
allocate(grid%element_errind(size(grid%element),max(1,nev)),stat=astat1)
if (astat1 /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in allocate_grid")
   deallocate(grid%element_errind,stat=dstat)
   deallocate(grid%head_level_elem,grid%head_level_vert,stat=dstat)
   call hash_table_destroy(grid%elem_hash)
   call hash_table_destroy(grid%edge_hash)
   call hash_table_destroy(grid%vert_hash)
   deallocate(grid%element,grid%edge,grid%vertex,grid%errest_energy, &
              grid%errest_Linf,grid%errest_L2,grid%eigenvalue, &
              grid%eigenprob_l2_resid,grid%eigenprob_variance, &
              grid%errest_eigenvalue,stat=dstat)
   if (associated(grid%eigenvalue)) then
      deallocate(grid%eigenvalue,grid%eigenprob_l2_resid, &
                 grid%eigenprob_variance,stat=dstat)
   endif
   return
endif
grid%element_errind = 0.0_my_real

! edges

grid%edge%degree = 1
do i=1,size(grid%edge)
   grid%edge(i)%next = i+1
end do
grid%edge(size(grid%edge))%next = END_OF_LIST
allocate(grid%edge_type(size(grid%edge),grid%system_size),stat=astat1)
if (astat1 /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in allocate_grid")
   deallocate(grid%edge_type,stat=dstat)
   deallocate(grid%element_errind,stat=dstat)
   deallocate(grid%head_level_elem,grid%head_level_vert,stat=dstat)
   call hash_table_destroy(grid%elem_hash)
   call hash_table_destroy(grid%edge_hash)
   call hash_table_destroy(grid%vert_hash)
   deallocate(grid%element,grid%edge,grid%vertex,grid%errest_energy, &
              grid%errest_Linf,grid%errest_L2,grid%eigenvalue, &
              grid%eigenprob_l2_resid,grid%eigenprob_variance, &
              grid%errest_eigenvalue,stat=dstat)
   if (associated(grid%eigenvalue)) then
      deallocate(grid%eigenvalue,grid%eigenprob_l2_resid, &
                 grid%eigenprob_variance,stat=dstat)
   endif
   return
endif
grid%edge_type = 0
do i=1,size(grid%edge)
   nullify(grid%edge(i)%solution,grid%edge(i)%exact,grid%edge(i)%oldsoln)
end do

! vertices

do i=1,size(grid%vertex)
   grid%vertex(i)%previous = i-1
   grid%vertex(i)%next = i+1
end do
grid%vertex(1)%previous = END_OF_LIST
grid%vertex(size(grid%vertex))%next = END_OF_LIST
allocate(grid%vertex_type(size(grid%vertex),grid%system_size),stat=astat1)
if (astat1 /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in allocate_grid")
   deallocate(grid%edge_type,stat=dstat)
   deallocate(grid%vertex_type,stat=dstat)
   deallocate(grid%element_errind,stat=dstat)
   deallocate(grid%head_level_elem,grid%head_level_vert,stat=dstat)
   call hash_table_destroy(grid%elem_hash)
   call hash_table_destroy(grid%edge_hash)
   call hash_table_destroy(grid%vert_hash)
   deallocate(grid%element,grid%edge,grid%vertex,grid%errest_energy, &
              grid%errest_Linf,grid%errest_L2,grid%eigenvalue, &
              grid%eigenprob_l2_resid,grid%eigenprob_variance, &
              grid%errest_eigenvalue,stat=dstat)
   if (associated(grid%eigenvalue)) then
      deallocate(grid%eigenvalue,grid%eigenprob_l2_resid, &
                 grid%eigenprob_variance,stat=dstat)
   endif
   return
endif
grid%vertex_type = 0
allocate(grid%vertex_solution(size(grid%vertex),grid%system_size,nevdim), &
         stat=astat1)
if (astat1 /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in allocate_grid")
   deallocate(grid%edge_type,grid%vertex_type,stat=dstat)
   deallocate(grid%element_errind,stat=dstat)
   deallocate(grid%head_level_elem,grid%head_level_vert,stat=dstat)
   call hash_table_destroy(grid%elem_hash)
   call hash_table_destroy(grid%edge_hash)
   call hash_table_destroy(grid%vert_hash)
   deallocate(grid%element,grid%edge,grid%vertex,grid%errest_energy, &
              grid%errest_Linf,grid%errest_L2,grid%eigenvalue, &
              grid%eigenprob_l2_resid,grid%eigenprob_variance, &
              grid%errest_eigenvalue,stat=dstat)
   if (associated(grid%eigenvalue)) then
      deallocate(grid%eigenvalue,grid%eigenprob_l2_resid, &
                 grid%eigenprob_variance,stat=dstat)
   endif
   return
endif
nullify(grid%vertex_exact,grid%vertex_oldsoln)
grid%vertex_solution = 0.0_my_real
grid%nsoln = grid%system_size*nevdim
grid%num_eval = nev
grid%have_true = .false.

! boundary parameter ranges, if using boundary subroutines

if (boundary_npiece(0) > 0) then
   npiece = boundary_npiece(0)
   i = 1
   do while (boundary_npiece(i) > 0)
      npiece = npiece + boundary_npiece(i)
      i = i+1
   end do
   allocate(grid%bp_start(npiece),grid%bp_finish(npiece),stat=astat1)
   if (astat1 /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in allocate_grid")
      stop
   endif
   call boundary_param(grid%bp_start,grid%bp_finish)
else
   nullify(grid%bp_start,grid%bp_finish)
endif

grid%errind_up2date = .false.
grid%oldsoln_exists = .false.

!call count_memory(grid)

end subroutine allocate_grid

!          ---------------
subroutine deallocate_grid(grid)
!          ---------------

!----------------------------------------------------
! This routine deallocates memory for a grid
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid

!----------------------------------------------------
! Local variables:

integer :: i, dstat

!----------------------------------------------------
! Begin executable code

if (associated(grid%vertex_type)) deallocate(grid%vertex_type,stat=dstat)
if (associated(grid%vertex_solution)) deallocate(grid%vertex_solution,stat=dstat)
if (associated(grid%vertex_exact)) deallocate(grid%vertex_exact,stat=dstat)
if (associated(grid%vertex_oldsoln)) deallocate(grid%vertex_oldsoln,stat=dstat)
if (associated(grid%edge_type)) deallocate(grid%edge_type,stat=dstat)
do i=1,size(grid%edge)
   if (associated(grid%edge(i)%solution)) deallocate(grid%edge(i)%solution,stat=dstat)
   if (associated(grid%edge(i)%exact)) deallocate(grid%edge(i)%exact,stat=dstat)
   if (associated(grid%edge(i)%oldsoln)) deallocate(grid%edge(i)%oldsoln,stat=dstat)
end do
if (associated(grid%element_errind)) deallocate(grid%element_errind,stat=dstat)
do i=1,size(grid%element)
   if (associated(grid%element(i)%solution)) deallocate(grid%element(i)%solution,stat=dstat)
   if (associated(grid%element(i)%exact)) deallocate(grid%element(i)%exact,stat=dstat)
   if (associated(grid%element(i)%oldsoln)) deallocate(grid%element(i)%oldsoln,stat=dstat)
end do
deallocate(grid%element, grid%edge, grid%vertex, grid%head_level_elem, &
           grid%head_level_vert, grid%errest_energy, grid%errest_Linf, &
           grid%errest_L2, grid%errest_eigenvalue,stat=dstat)
if (associated(grid%eigenvalue)) deallocate(grid%eigenvalue, &
                                            grid%eigenprob_variance, &
                                            grid%eigenprob_l2_resid,stat=dstat)
if (associated(grid%initial_neighbor)) &
   deallocate(grid%initial_neighbor, stat=dstat)
call hash_table_destroy(grid%elem_hash)
call hash_table_destroy(grid%edge_hash)
call hash_table_destroy(grid%vert_hash)
if (associated(grid%bp_start)) deallocate(grid%bp_start,grid%bp_finish,stat=dstat)

end subroutine deallocate_grid

!          ----------------
subroutine realloc_solution(grid,degree,system_size,neval)
!          ----------------

!----------------------------------------------------
! This routine reallocates the solution components of the grid for
! the given degree, system size and number of eigenvalues, keeping the first
! nsolut old solutions and filling in 0.0 if nsolut is bigger than
! the old number of solutions
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: degree, system_size, neval
!----------------------------------------------------
! Local variables:

integer :: i, ncopy, ncopyss, ncopyev, astat, nsolut, neq_edge, neq_elem, &
           neq_copy, deg, lev, elem, edge
real(my_real), pointer :: temp3(:,:,:)
!----------------------------------------------------
! Begin executable code

nsolut = system_size*neval
ncopy = min(nsolut,grid%nsoln)
ncopyss = min(size(grid%vertex_solution,2),system_size)
ncopyev = min(size(grid%vertex_solution,3),neval)

deallocate(grid%errest_energy, grid%errest_Linf, grid%errest_L2, &
           grid%errest_eigenvalue,stat=astat)
allocate(grid%errest_energy(neval), grid%errest_Linf(nsolut), &
         grid%errest_L2(nsolut), grid%errest_eigenvalue(neval),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in realloc_solution")
   return
endif

deallocate(grid%element_errind,stat=astat)
allocate(grid%element_errind(size(grid%element),neval),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in realloc_solution")
   return
endif

! vertices

if (size(grid%vertex_solution,2) /= system_size .or. &
    size(grid%vertex_solution,3) /= neval) then
   allocate(temp3(size(grid%vertex),system_size,neval),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in realloc_solution")
      return
   endif
   temp3 = 0.0_my_real
   temp3(:,1:ncopyss,1:ncopyev) = grid%vertex_solution(:,1:ncopyss,1:ncopyev)
   deallocate(grid%vertex_solution,stat=astat)
   grid%vertex_solution => temp3
   nullify(temp3)
endif

if (grid%have_true) then
 if (size(grid%vertex_exact,2) /= system_size .or. &
     size(grid%vertex_exact,3) /= neval) then
   allocate(temp3(size(grid%vertex),system_size,neval),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in realloc_solution")
      return
   endif
   temp3 = 0.0_my_real
   temp3(:,1:ncopyss,1:ncopyev) = grid%vertex_exact(:,1:ncopyss,1:ncopyev)
   deallocate(grid%vertex_exact,stat=astat)
   grid%vertex_exact => temp3
   nullify(temp3)
 endif
endif

do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)

! edges

      do i=1,EDGES_PER_ELEMENT
         edge = grid%element(elem)%edge(i)
         if (degree > 0) then
            deg = degree
            if (grid%edge(edge)%degree > 0) then
               grid%dof = grid%dof - (grid%edge(edge)%degree-1)
               if (grid%element(grid%edge(edge)%assoc_elem)%iown) then
                  grid%dof_own = grid%dof_own - (grid%edge(edge)%degree-1)
               endif
            endif
            grid%edge(edge)%degree = degree
            grid%dof = grid%dof + (grid%edge(edge)%degree-1)
            if (grid%element(grid%edge(edge)%assoc_elem)%iown) then
               grid%dof_own = grid%dof_own + (grid%edge(edge)%degree-1)
            endif
         else
            deg = grid%edge(edge)%degree
         endif
         neq_edge = deg-1
         if (deg <= 1) then
            if (associated(grid%edge(edge)%solution)) then
               deallocate(grid%edge(edge)%solution)
            endif
            if (associated(grid%edge(edge)%exact)) then
               deallocate(grid%edge(edge)%exact)
            endif
            nullify(grid%edge(edge)%solution,grid%edge(edge)%exact)
         else
            if (associated(grid%edge(edge)%solution)) then
               if (size(grid%edge(edge)%solution,dim=1) == neq_edge .and. &
                   size(grid%edge(edge)%solution,dim=2) == system_size .and. &
                   size(grid%edge(edge)%solution,dim=3) == neval) cycle
            endif
            if (associated(grid%edge(edge)%solution)) then
               neq_copy = min(neq_edge,size(grid%edge(edge)%solution,dim=1))
            else
               neq_copy = 0
            endif
            allocate(temp3(neq_edge,system_size,neval),stat=astat)
            if (astat /= 0) then
               ierr = ALLOC_FAILED
               call fatal("memory allocation failed in realloc_solution")
               return
            endif
            temp3 = 0.0_my_real
            if (neq_copy > 0) then
               temp3(1:neq_copy,1:ncopyss,1:ncopyev) = grid%edge(edge)%solution(1:neq_copy,1:ncopyss,1:ncopyev)
            endif
            if (associated(grid%edge(edge)%solution)) then
               deallocate(grid%edge(edge)%solution,stat=astat)
            endif
            grid%edge(edge)%solution => temp3
            nullify(temp3)
            if (grid%have_true) then
               if (associated(grid%edge(edge)%exact)) then
                  neq_copy = min(neq_edge,size(grid%edge(edge)%exact,dim=1))
               else
                  neq_copy = 0
               endif
               allocate(temp3(neq_edge,system_size,neval),stat=astat)
               if (astat /= 0) then
                  ierr = ALLOC_FAILED
                  call fatal("memory allocation failed in realloc_solution")
                  return
               endif
               temp3 = 0.0_my_real
               if (neq_copy > 0) then
                  temp3(1:neq_copy,1:ncopyss,1:ncopyev) = grid%edge(edge)%exact(1:neq_copy,1:ncopyss,1:ncopyev)
               endif
               if (associated(grid%edge(edge)%exact)) then
                  deallocate(grid%edge(edge)%exact,stat=astat)
               endif
               grid%edge(edge)%exact => temp3
               nullify(temp3)
            endif
         endif
      end do

! element

      if (degree > 0) then
         deg = degree
         if (grid%element(elem)%degree > 0) then
            grid%dof = grid%dof - &
                 ((grid%element(elem)%degree-2)*(grid%element(elem)%degree-1))/2
            if (grid%element(elem)%iown) then
               grid%dof_own = grid%dof_own - &
                 ((grid%element(elem)%degree-2)*(grid%element(elem)%degree-1))/2
            endif
         endif
         grid%element(elem)%degree = degree
         grid%dof = grid%dof + &
                 ((grid%element(elem)%degree-2)*(grid%element(elem)%degree-1))/2
         if (grid%element(elem)%iown) then
            grid%dof_own = grid%dof_own + &
                 ((grid%element(elem)%degree-2)*(grid%element(elem)%degree-1))/2
         endif
      else
         deg = grid%element(elem)%degree
      endif
      neq_elem = ((deg-2)*(deg-1))/2
      if (deg <= 2) then
         if (associated(grid%element(elem)%solution)) then
            deallocate(grid%element(elem)%solution)
         endif
         if (associated(grid%element(elem)%exact)) then
            deallocate(grid%element(elem)%exact)
         endif
         nullify(grid%element(elem)%solution,grid%element(elem)%exact)
      else
         if (associated(grid%element(elem)%solution)) then
            if (size(grid%element(elem)%solution,dim=1) == neq_elem .and. &
                size(grid%element(elem)%solution,dim=2) == system_size .and. &
                size(grid%element(elem)%solution,dim=3) == neval) then
               elem = grid%element(elem)%next
               cycle
            endif
         endif
         if (associated(grid%element(elem)%solution)) then
            neq_copy = min(neq_elem,size(grid%element(elem)%solution,dim=1))
         else
            neq_copy = 0
         endif
         allocate(temp3(neq_elem,system_size,neval),stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("memory allocation failed in realloc_solution")
            return
         endif
         temp3 = 0.0_my_real
         if (neq_copy > 0) then
            temp3(1:neq_copy,1:ncopyss,1:ncopyev) = grid%element(elem)%solution(1:neq_copy,1:ncopyss,1:ncopyev)
         endif
         if (associated(grid%element(elem)%solution)) then
            deallocate(grid%element(elem)%solution,stat=astat)
         endif
         grid%element(elem)%solution => temp3
         nullify(temp3)
         if (grid%have_true) then
            if (associated(grid%element(elem)%exact)) then
               neq_copy = min(neq_elem,size(grid%element(elem)%exact,dim=1))
            else
               neq_copy = 0
            endif
            allocate(temp3(neq_elem,system_size,neval),stat=astat)
            if (astat /= 0) then
               ierr = ALLOC_FAILED
               call fatal("memory allocation failed in realloc_solution")
               return
            endif
            temp3 = 0.0_my_real
            if (neq_copy > 0) then
               temp3(1:neq_copy,1:ncopyss,1:ncopyev) = grid%element(elem)%exact(1:neq_copy,1:ncopyss,1:ncopyev)
            endif
            if (associated(grid%element(elem)%exact)) then
               deallocate(grid%element(elem)%exact,stat=astat)
            endif
            grid%element(elem)%exact => temp3
            nullify(temp3)
         endif
      endif
      elem = grid%element(elem)%next
   end do
end do

grid%nsoln = nsolut
grid%errind_up2date = .false.

end subroutine realloc_solution

!          ------------------
subroutine reset_dirich_exact(grid)
!          ------------------

!----------------------------------------------------
! This routine resets Dirichlet boundary conditions and the exact solution
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
!----------------------------------------------------
! Local variables:

integer :: lev, vert, edge, elem, i, j
integer :: bctype(grid%system_size)
real(my_real) :: bcrhs(grid%system_size), &
                 bccoef(grid%system_size,grid%system_size)
logical(small_logical) :: edge_done(size(grid%edge))
!----------------------------------------------------
! Begin executable code

! for each vertex, set Dirichlet b.c. and exact

do lev=1,grid%nlev
   vert = grid%head_level_vert(lev)
   do while (vert /= END_OF_LIST)
      if (any(grid%vertex_type(vert,:) == DIRICHLET) .or. &
          any(grid%vertex_type(vert,:) == PERIODIC_SLAVE_DIR) .or. &
          any(grid%vertex_type(vert,:) == PERIODIC_MASTER_DIR)) then
         call bconds(grid%vertex(vert)%coord%x,grid%vertex(vert)%coord%y, &
                     grid%vertex(vert)%bmark,bctype,bccoef,bcrhs)
         do i=1,grid%system_size
           if (grid%vertex_type(vert,i) == DIRICHLET .or. &
               grid%vertex_type(vert,i) == PERIODIC_SLAVE_DIR .or. &
               grid%vertex_type(vert,i) == PERIODIC_MASTER_DIR) then
            grid%vertex_solution(vert,i,:) = bcrhs(i)
           endif
         end do
      endif
      if (grid%have_true) then
         do j=1,max(1,grid%num_eval)
            do i=1,grid%system_size
               grid%vertex_exact(vert,i,j) = trues(grid%vertex(vert)%coord%x, &
                                                   grid%vertex(vert)%coord%y, &
                                                   i,j)
            end do
         end do
      endif
      vert = grid%vertex(vert)%next
   end do
end do

! for each leaf element
!   for each edge of element that has not already been done
!      set Dirichlet b.c. and exact
!   set exact

edge_done = .false.
do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (grid%element(elem)%isleaf) then
         do j=1,EDGES_PER_ELEMENT
            edge = grid%element(elem)%edge(j)
            if (.not. edge_done(edge)) then
               do i=1,grid%system_size
                  if (grid%edge_type(edge,i) == DIRICHLET) then
                     call edge_exact(grid,edge,i,"d")
                  endif
                  if (grid%have_true) call edge_exact(grid,edge,i,"t")
               end do
            endif
         end do
         if (grid%have_true) then
            do i=1,grid%system_size
               call elem_exact(grid,elem,i,"t")
            end do
         endif
      endif
      elem = grid%element(elem)%next
   end do
end do

end subroutine reset_dirich_exact

!          ---------
subroutine copy_grid(old_grid,new_grid)
!          ---------

!----------------------------------------------------
! This routine makes an exact copy of old_grid in new_grid
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: old_grid
type(grid_type), intent(out) :: new_grid
!----------------------------------------------------
! Local variables:

integer :: i
!----------------------------------------------------
! Begin executable code

! elements

if (associated(old_grid%element)) then
   allocate(new_grid%element(size(old_grid%element)))
   do i=1,size(old_grid%element)
      new_grid%element(i)%gid = old_grid%element(i)%gid
      new_grid%element(i)%mate = old_grid%element(i)%mate
      new_grid%element(i)%weight = old_grid%element(i)%weight
      if (associated(old_grid%element(i)%solution)) then
         allocate(new_grid%element(i)%solution( &
                  size(old_grid%element(i)%solution,dim=1), &
                  size(old_grid%element(i)%solution,dim=2), &
                  size(old_grid%element(i)%solution,dim=3)))
         new_grid%element(i)%solution = old_grid%element(i)%solution
      else
         nullify(new_grid%element(i)%solution)
      endif
      if (associated(old_grid%element(i)%exact)) then
         allocate(new_grid%element(i)%exact( &
                  size(old_grid%element(i)%exact,dim=1), &
                  size(old_grid%element(i)%exact,dim=2), &
                  size(old_grid%element(i)%exact,dim=3)))
         new_grid%element(i)%exact = old_grid%element(i)%exact
      else
         nullify(new_grid%element(i)%exact)
      endif
      if (associated(old_grid%element(i)%oldsoln)) then
         allocate(new_grid%element(i)%oldsoln( &
                  size(old_grid%element(i)%oldsoln,dim=1), &
                  size(old_grid%element(i)%oldsoln,dim=2), &
                  size(old_grid%element(i)%oldsoln,dim=3)))
         new_grid%element(i)%oldsoln = old_grid%element(i)%oldsoln
      else
         nullify(new_grid%element(i)%oldsoln)
      endif
      new_grid%element(i)%work = old_grid%element(i)%work
      new_grid%element(i)%sp_eta_pred = old_grid%element(i)%sp_eta_pred
      new_grid%element(i)%vertex = old_grid%element(i)%vertex
      new_grid%element(i)%edge = old_grid%element(i)%edge
      new_grid%element(i)%degree = old_grid%element(i)%degree
      new_grid%element(i)%level = old_grid%element(i)%level
      new_grid%element(i)%in = old_grid%element(i)%in
      new_grid%element(i)%out = old_grid%element(i)%out
      new_grid%element(i)%order = old_grid%element(i)%order
      new_grid%element(i)%next = old_grid%element(i)%next
      new_grid%element(i)%previous = old_grid%element(i)%previous
      new_grid%element(i)%isleaf = old_grid%element(i)%isleaf
      new_grid%element(i)%oldleaf = old_grid%element(i)%oldleaf
      new_grid%element(i)%iown = old_grid%element(i)%iown
      new_grid%element(i)%hrefined_unowned = old_grid%element(i)%hrefined_unowned
      new_grid%element(i)%prefined_unowned = old_grid%element(i)%prefined_unowned
   end do
else
   nullify(new_grid%element)
endif

! edges

if (associated(old_grid%edge)) then
   allocate(new_grid%edge(size(old_grid%edge)))
   do i=1,size(old_grid%edge)
      new_grid%edge(i)%gid = old_grid%edge(i)%gid
      new_grid%edge(i)%vertex = old_grid%edge(i)%vertex
      new_grid%edge(i)%bmark = old_grid%edge(i)%bmark
      new_grid%edge(i)%degree = old_grid%edge(i)%degree
      new_grid%edge(i)%assoc_elem = old_grid%edge(i)%assoc_elem
      if (associated(old_grid%edge(i)%solution)) then
         allocate(new_grid%edge(i)%solution( &
                  size(old_grid%edge(i)%solution,dim=1), &
                  size(old_grid%edge(i)%solution,dim=2), &
                  size(old_grid%edge(i)%solution,dim=3)))
         new_grid%edge(i)%solution = old_grid%edge(i)%solution
      else
         nullify(new_grid%edge(i)%solution)
      endif
      if (associated(old_grid%edge(i)%exact)) then
         allocate(new_grid%edge(i)%exact( &
                  size(old_grid%edge(i)%exact,dim=1), &
                  size(old_grid%edge(i)%exact,dim=2), &
                  size(old_grid%edge(i)%exact,dim=3)))
         new_grid%edge(i)%exact = old_grid%edge(i)%exact
      else
         nullify(new_grid%edge(i)%exact)
      endif
      if (associated(old_grid%edge(i)%oldsoln)) then
         allocate(new_grid%edge(i)%oldsoln( &
                  size(old_grid%edge(i)%oldsoln,dim=1), &
                  size(old_grid%edge(i)%oldsoln,dim=2), &
                  size(old_grid%edge(i)%oldsoln,dim=3)))
         new_grid%edge(i)%oldsoln = old_grid%edge(i)%oldsoln
      else
         nullify(new_grid%edge(i)%oldsoln)
      endif
      new_grid%edge(i)%next = old_grid%edge(i)%next
   end do
else
   nullify(new_grid%edge)
endif

! vertices

if (associated(old_grid%vertex)) then
   allocate(new_grid%vertex(size(old_grid%vertex)))
   do i=1,size(old_grid%vertex)
      new_grid%vertex(i)%gid = old_grid%vertex(i)%gid
      new_grid%vertex(i)%coord = old_grid%vertex(i)%coord
      new_grid%vertex(i)%bparam = old_grid%vertex(i)%bparam
      new_grid%vertex(i)%bmark = old_grid%vertex(i)%bmark
      new_grid%vertex(i)%assoc_elem = old_grid%vertex(i)%assoc_elem
      new_grid%vertex(i)%next = old_grid%vertex(i)%next
      new_grid%vertex(i)%previous = old_grid%vertex(i)%previous
   end do
else
   nullify(new_grid%vertex)
endif

! remaining fields

call hash_table_copy(old_grid%elem_hash,new_grid%elem_hash)
call hash_table_copy(old_grid%edge_hash,new_grid%edge_hash)
call hash_table_copy(old_grid%vert_hash,new_grid%vert_hash)
new_grid%boundbox_min = old_grid%boundbox_min
new_grid%boundbox_max = old_grid%boundbox_max
if (associated(old_grid%vertex_solution)) then
   allocate(new_grid%vertex_solution( &
            size(old_grid%vertex_solution,dim=1), &
            size(old_grid%vertex_solution,dim=2), &
            size(old_grid%vertex_solution,dim=3)))
   new_grid%vertex_solution = old_grid%vertex_solution
else
   nullify(new_grid%vertex_solution)
endif
if (associated(old_grid%vertex_exact)) then
   allocate(new_grid%vertex_exact( &
            size(old_grid%vertex_exact,dim=1), &
            size(old_grid%vertex_exact,dim=2), &
            size(old_grid%vertex_exact,dim=3)))
   new_grid%vertex_exact = old_grid%vertex_exact
else
   nullify(new_grid%vertex_exact)
endif
if (associated(old_grid%vertex_oldsoln)) then
   allocate(new_grid%vertex_oldsoln( &
            size(old_grid%vertex_oldsoln,dim=1), &
            size(old_grid%vertex_oldsoln,dim=2), &
            size(old_grid%vertex_oldsoln,dim=3)))
   new_grid%vertex_oldsoln = old_grid%vertex_oldsoln
else
   nullify(new_grid%vertex_oldsoln)
endif
if (associated(old_grid%element_errind)) then
   allocate(new_grid%element_errind( &
            size(old_grid%element_errind,dim=1), &
            size(old_grid%element_errind,dim=2)))
   new_grid%element_errind = old_grid%element_errind
else
   nullify(new_grid%element_errind)
endif
if (associated(old_grid%eigenvalue)) then
   allocate(new_grid%eigenvalue(size(old_grid%eigenvalue)))
   new_grid%eigenvalue = old_grid%eigenvalue
else
   nullify(new_grid%eigenvalue)
endif
new_grid%eigen_linsys_max_l2_resid = old_grid%eigen_linsys_max_l2_resid
new_grid%eigen_linsys_ave_l2_resid = old_grid%eigen_linsys_ave_l2_resid
if (associated(old_grid%eigenprob_l2_resid)) then
   allocate(new_grid%eigenprob_l2_resid(size(old_grid%eigenprob_l2_resid)))
   new_grid%eigenprob_l2_resid = old_grid%eigenprob_l2_resid
else
   nullify(new_grid%eigenprob_l2_resid)
endif
if (associated(old_grid%eigenprob_variance)) then
   allocate(new_grid%eigenprob_variance(size(old_grid%eigenprob_variance)))
   new_grid%eigenprob_variance = old_grid%eigenprob_variance
else
   nullify(new_grid%eigenprob_variance)
endif
if (associated(old_grid%errest_energy)) then
   allocate(new_grid%errest_energy(size(old_grid%errest_energy)))
   new_grid%errest_energy = old_grid%errest_energy
else
   nullify(new_grid%errest_energy)
endif
if (associated(old_grid%errest_Linf)) then
   allocate(new_grid%errest_Linf(size(old_grid%errest_Linf)))
   new_grid%errest_Linf = old_grid%errest_Linf
else
   nullify(new_grid%errest_Linf)
endif
if (associated(old_grid%errest_L2)) then
   allocate(new_grid%errest_L2(size(old_grid%errest_L2)))
   new_grid%errest_L2 = old_grid%errest_L2
else
   nullify(new_grid%errest_L2)
endif
if (associated(old_grid%errest_eigenvalue)) then
   allocate(new_grid%errest_eigenvalue(size(old_grid%errest_eigenvalue)))
   new_grid%errest_eigenvalue = old_grid%errest_eigenvalue
else
   nullify(new_grid%errest_eigenvalue)
endif
new_grid%refsoln_errest = old_grid%refsoln_errest
new_grid%max_blen = old_grid%max_blen
if (associated(old_grid%bp_start)) then
   allocate(new_grid%bp_start(size(old_grid%bp_start)))
   new_grid%bp_start = old_grid%bp_start
else
   nullify(new_grid%bp_start)
endif
if (associated(old_grid%bp_finish)) then
   allocate(new_grid%bp_finish(size(old_grid%bp_finish)))
   new_grid%bp_finish = old_grid%bp_finish
else
   nullify(new_grid%bp_finish)
endif
if (associated(old_grid%edge_type)) then
   allocate(new_grid%edge_type( &
            size(old_grid%edge_type,dim=1), &
            size(old_grid%edge_type,dim=2)))
   new_grid%edge_type = old_grid%edge_type
else
   nullify(new_grid%edge_type)
endif
if (associated(old_grid%vertex_type)) then
   allocate(new_grid%vertex_type( &
            size(old_grid%vertex_type,dim=1), &
            size(old_grid%vertex_type,dim=2)))
   new_grid%vertex_type = old_grid%vertex_type
else
   nullify(new_grid%vertex_type)
endif
if (associated(old_grid%initial_neighbor)) then
   allocate(new_grid%initial_neighbor( &
            size(old_grid%initial_neighbor,dim=1), &
            size(old_grid%initial_neighbor,dim=2)))
   new_grid%initial_neighbor = old_grid%initial_neighbor
else
   nullify(new_grid%initial_neighbor)
endif
if (associated(old_grid%head_level_elem)) then
   allocate(new_grid%head_level_elem(size(old_grid%head_level_elem)))
   new_grid%head_level_elem = old_grid%head_level_elem
else
   nullify(new_grid%head_level_elem)
endif
if (associated(old_grid%head_level_vert)) then
   allocate(new_grid%head_level_vert(size(old_grid%head_level_vert)))
   new_grid%head_level_vert = old_grid%head_level_vert
else
   nullify(new_grid%head_level_vert)
endif
new_grid%next_free_elem = old_grid%next_free_elem
new_grid%next_free_edge = old_grid%next_free_edge
new_grid%next_free_vert = old_grid%next_free_vert
new_grid%partition = old_grid%partition
new_grid%system_size = old_grid%system_size
new_grid%num_eval = old_grid%num_eval
new_grid%nsoln = old_grid%nsoln
new_grid%nelem = old_grid%nelem
new_grid%nelem_leaf = old_grid%nelem_leaf
new_grid%nelem_leaf_own = old_grid%nelem_leaf_own
new_grid%nedge = old_grid%nedge
new_grid%nedge_own = old_grid%nedge_own
new_grid%nvert = old_grid%nvert
new_grid%nvert_own = old_grid%nvert_own
new_grid%nlev = old_grid%nlev
new_grid%dof = old_grid%dof
new_grid%dof_own = old_grid%dof_own
new_grid%arpack_iter = old_grid%arpack_iter
new_grid%arpack_nconv = old_grid%arpack_nconv
new_grid%arpack_numop = old_grid%arpack_numop
new_grid%arpack_numopb = old_grid%arpack_numopb
new_grid%arpack_numreo = old_grid%arpack_numreo
new_grid%arpack_info = old_grid%arpack_info
new_grid%errtype = old_grid%errtype
new_grid%errind_up2date = old_grid%errind_up2date
new_grid%oldsoln_exists = old_grid%oldsoln_exists
new_grid%have_true = old_grid%have_true
new_grid%triangle_files = old_grid%triangle_files

end subroutine copy_grid

!---------------------------------------------------------------------
!  ROUTINES FOR CREATING THE INITIAL GRID
!---------------------------------------------------------------------

!          ---------
subroutine init_grid(grid,procs,degree,partition,set_iown)
!          ---------

!----------------------------------------------------
! This routine initializes the grid starting with the triangulation given by
! .node, .ele, .edge and .neigh files in the format generated by triangle.  The
! root name of the files is given in grid%triangle_files.  If the domain is
! given by boundary subroutines, it first creates the triangle data files.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(proc_info), intent(in) :: procs
integer, intent(in) :: degree,partition
logical, intent(in) :: set_iown

!----------------------------------------------------
! Local variables:

type(triangle_data) :: td
integer :: i, j, k, neigh, stat, nedge, edge, v1s, v2s, v1m, v2m, neighbmark, &
           igid
integer, allocatable :: mate(:)
integer :: elemedge_bmark(EDGES_PER_ELEMENT,size(grid%element))
integer :: bctype(grid%system_size)
real(my_real) :: bccoef(grid%system_size,grid%system_size), &
                 bcrhs(grid%system_size)
integer :: ni, nr
integer, pointer :: imess(:)
real(my_real), pointer :: rmess(:)

integer, parameter :: NO_MATE = BOUNDARY-1
!----------------------------------------------------
! Begin executable code

! create the triangle data files if the boundary is given by subroutines

if (boundary_npiece(0) > 0) then
   if (.false.) then
! TEMP this would cause multiple creations of the same file, in parallel
      call create_triangle_files(grid)
   else
! TEMP this requires all processes share the same file system
      if (my_proc(procs) == MASTER .or. PARALLEL==SEQUENTIAL) then
         call create_triangle_files(grid)
         do i=1,num_proc(procs)
            call phaml_send(procs,i,(/0/),0,(/0.0_my_real/),0,401)
         end do
      else
         call phaml_recv(procs,j,imess,ni,rmess,nr,401)
      endif
   endif
endif

! get the starting grid from the data files

call read_triangle_data(grid,td)

! included master for creating triangle data files and checking on opening
! triangle data files, but don't generate a grid for the master

if (my_proc(procs) == MASTER) then
   return
endif

! space for PERIODIC_MASTER vertex for PERIODIC_SLAVEs

allocate(td%vert_master(td%nvert),stat=stat)
if (stat /= 0) then
   call fatal("allocation failed for mates in init_grid")
   stop
endif

td%vert_master = NO_MATE

! match up sides with periodic boundary conditions

call match_periodic(grid,td,NO_MATE)

! pair up most triangles

allocate(mate(td%ntri),stat=stat)
if (stat /= 0) then
   call fatal("allocation failed for mates in init_grid")
   stop
endif

mate = NO_MATE
call pair_triangles(td,mate,NO_MATE)

! Refine the starting grid by bisecting paired triangles and trisecting others.
! This defines these components of grid: vertex%coord, element%gid, &
! element%mate, element%vertex, nelem, nvert, initial_neighbor
! It also sets vertex%assoc_elem to be the master vertex for vertices that
! are PERIODIC_SLAVE (negative bmark) or endpoints of a PERIODIC SLAVE side,
! and NO_MATE for all other vertices. This is just a convenient component to
! use for temporary storage.

call refine_start(grid,td,mate,elemedge_bmark,NO_MATE)

! Set the vertex linked list.  Need to do vertices on sides that have
! periodic boundary conditions first so that
! next(PERIODIC_SLAVE) = PERIODIC_MASTER.  And need to have nonpeaks before
! peaks so that ancestors of level 2 equations come before parents in the
! numbering of equations (see make_linear_system).  Since the linked list is
! built from the tail to the head, the order of building it is
!   peak, periodic
!   peak, nonperiodic
!   nonpeak, periodic
!   nonpeak, nonperiodic
! assoc_elem was used as a place to store the master associated with a slave.
! Also set vertex type to a form of PERIODIC or INTERIOR as an initialization.

grid%head_level_vert(1) = END_OF_LIST
do i=td%nvert+1,grid%nvert
   if (grid%vertex(i)%assoc_elem /= NO_MATE .and. &
       grid%vertex(i)%assoc_elem /= NO_MATE-1) then
      grid%vertex(grid%vertex(i)%assoc_elem)%next = grid%head_level_vert(1)
      grid%vertex(i)%next = grid%vertex(i)%assoc_elem
      grid%head_level_vert(1) = i
      grid%vertex_type(i,:) = PERIODIC_SLAVE
      grid%vertex_type(grid%vertex(i)%assoc_elem,:) = PERIODIC_MASTER
      grid%vertex(grid%vertex(i)%assoc_elem)%assoc_elem = NO_MATE-1
   endif
end do

do i=td%nvert+1,grid%nvert
   if (grid%vertex(i)%assoc_elem == NO_MATE) then
      grid%vertex(i)%next = grid%head_level_vert(1)
      grid%head_level_vert(1) = i
      grid%vertex_type(i,:) = INTERIOR
   endif
end do

do i=1,td%nvert
   if (grid%vertex(i)%assoc_elem /= NO_MATE .and. &
       grid%vertex(i)%assoc_elem /= NO_MATE-1) then
      grid%vertex(grid%vertex(i)%assoc_elem)%next = grid%head_level_vert(1)
      grid%vertex(i)%next = grid%vertex(i)%assoc_elem
      grid%head_level_vert(1) = i
      grid%vertex_type(i,:) = PERIODIC_SLAVE
      grid%vertex_type(grid%vertex(i)%assoc_elem,:) = PERIODIC_MASTER
      grid%vertex(grid%vertex(i)%assoc_elem)%assoc_elem = NO_MATE-1
   endif
end do

do i=1,td%nvert
   if (grid%vertex(i)%assoc_elem == NO_MATE) then
      grid%vertex(i)%next = grid%head_level_vert(1)
      grid%head_level_vert(1) = i
      grid%vertex_type(i,:) = INTERIOR
   endif
end do

i = grid%head_level_vert(1)
grid%vertex(i)%previous = END_OF_LIST
do while (i /= END_OF_LIST)
   j = grid%vertex(i)%next
   if (j /= END_OF_LIST) grid%vertex(j)%previous = i
   i = j
end do

! set remaining components of grid
! also set vertex 1 to be the smaller of vertex 1 and vertex 2 so that
! pairs agree

! define the elements

grid%head_level_elem(1) = 1
do i=1,grid%nelem
   if (grid%element(i)%vertex(1) > grid%element(i)%vertex(2)) then
      j = grid%element(i)%vertex(1)
      grid%element(i)%vertex(1) = grid%element(i)%vertex(2)
      grid%element(i)%vertex(2) = j
      j = grid%initial_neighbor(1,i)
      grid%initial_neighbor(1,i) = grid%initial_neighbor(2,i)
      grid%initial_neighbor(2,i) = j
      j = elemedge_bmark(1,i)
      elemedge_bmark(1,i) = elemedge_bmark(2,i)
      elemedge_bmark(2,i) = j
   endif
   grid%element(i)%degree = degree
   if (degree > 2) then
      allocate(grid%element(i)%solution(((degree-1)*(degree-2))/2,grid%system_size,max(1,grid%num_eval)),&
               stat=stat)
      if (stat /= 0) then
         call fatal("allocation failed for mates in init_grid")
         stop
      endif
      grid%element(i)%solution = 0.0_my_real
      if (grid%have_true) then
         allocate(grid%element(i)%exact(((degree-1)*(degree-2))/2,grid%system_size,max(1,grid%num_eval)),&
                  stat=stat)
         if (stat /= 0) then
            call fatal("allocation failed for mates in init_grid")
            stop
         endif
         grid%element(i)%exact = 0.0_my_real
      endif
   endif
   grid%element(i)%in = grid%element(i)%vertex(1)
   grid%element(i)%out = grid%element(i)%vertex(2)
   grid%element(i)%order = (/1,2/)
   grid%element(i)%level = 1
   grid%element(i)%next = i+1
   grid%element(i)%iown = set_iown
   grid%element(i)%hrefined_unowned = .false.
   grid%element(i)%prefined_unowned = .false.
   grid%element(i)%isleaf = .true.
   grid%element(i)%oldleaf = .false.
   grid%element(i)%sp_eta_pred = 0.0_my_real
   call hash_insert(grid%element(i)%gid,i,grid%elem_hash)
   do j=1,3
      grid%vertex(grid%element(i)%vertex(j))%assoc_elem = i
      if (grid%initial_neighbor(j,i) == BOUNDARY) then
         where (grid%vertex_type(grid%element(i)%vertex(1+mod(j,3)),:) == INTERIOR) &
            grid%vertex_type(grid%element(i)%vertex(1+mod(j,3)),:) = DIRICHLET
         where (grid%vertex_type(grid%element(i)%vertex(1+mod(j+1,3)),:) == INTERIOR) &
            grid%vertex_type(grid%element(i)%vertex(1+mod(j+1,3)),:) = DIRICHLET
      endif
   end do
end do
grid%element(grid%nelem)%next = END_OF_LIST

grid%next_free_elem = grid%nelem+1
grid%nelem_leaf = grid%nelem
if (set_iown) then
   grid%nelem_leaf_own = grid%nelem
else
   grid%nelem_leaf_own = 0
endif

! define the vertices

grid%boundbox_min = point(huge(0.0_my_real),huge(0.0_my_real))
grid%boundbox_max = point(-huge(0.0_my_real),-huge(0.0_my_real))
do i=1,grid%nvert
   grid%vertex(i)%gid = i
   if (grid%vertex_type(i,1) == INTERIOR) then
      grid%vertex_solution(i,:,:) = 0.0_my_real
   else
      call bconds(grid%vertex(i)%coord%x,grid%vertex(i)%coord%y, &
                  grid%vertex(i)%bmark,bctype,bccoef,bcrhs)
      do j=1,grid%system_size
         if (grid%vertex_type(i,j) == PERIODIC_SLAVE .or. &
             grid%vertex_type(i,j) == PERIODIC_MASTER) then
            select case(bctype(j))
            case (DIRICHLET)
               if (grid%vertex_type(i,j) == PERIODIC_SLAVE) then
                  grid%vertex_type(i,j) = PERIODIC_SLAVE_DIR
               else
                  grid%vertex_type(i,j) = PERIODIC_MASTER_DIR
               endif
            case (NATURAL)
               if (grid%vertex_type(i,j) == PERIODIC_SLAVE) then
                  grid%vertex_type(i,j) = PERIODIC_SLAVE_NAT
               else
                  grid%vertex_type(i,j) = PERIODIC_MASTER_NAT
               endif
            case (MIXED)
               if (grid%vertex_type(i,j) == PERIODIC_SLAVE) then
                  grid%vertex_type(i,j) = PERIODIC_SLAVE_MIX
               else
                  grid%vertex_type(i,j) = PERIODIC_MASTER_MIX
               endif
            case (PERIODIC)
               if (grid%vertex(i)%bmark < 0) then
                  grid%vertex_type(i,j) = PERIODIC_SLAVE
               else
                  grid%vertex_type(i,j) = PERIODIC_MASTER
               endif
            end select
         else
            grid%vertex_type(i,j) = bctype(j)
         endif
         if (grid%vertex_type(i,j) == DIRICHLET .or. &
             grid%vertex_type(i,j) == PERIODIC_SLAVE_DIR .or. &
             grid%vertex_type(i,j) == PERIODIC_MASTER_DIR) then
            grid%vertex_solution(i,j,:) = bcrhs(j)
         else
            grid%vertex_solution(i,j,:) = 0.0_my_real
         endif
      end do
   endif
   call hash_insert(grid%vertex(i)%gid,i,grid%vert_hash)
   if (grid%vertex(i)%coord%x < grid%boundbox_min%x) then
      grid%boundbox_min%x = grid%vertex(i)%coord%x
   endif
   if (grid%vertex(i)%coord%x > grid%boundbox_max%x) then
      grid%boundbox_max%x = grid%vertex(i)%coord%x
   endif
   if (grid%vertex(i)%coord%y < grid%boundbox_min%y) then
      grid%boundbox_min%y = grid%vertex(i)%coord%y
   endif
   if (grid%vertex(i)%coord%y > grid%boundbox_max%y) then
      grid%boundbox_max%y = grid%vertex(i)%coord%y
   endif
end do
grid%next_free_vert = grid%nvert+1
if (set_iown) then
   grid%nvert_own = grid%nvert
else
   grid%nvert_own = 0
endif

! check for true solution known

if (trues((grid%boundbox_max%x+grid%boundbox_min%x)/2, &
          (grid%boundbox_max%y+grid%boundbox_min%y)/2, 1, 1) &
    == huge(0.0_my_real)) then
   grid%have_true = .false.
else
   grid%have_true = .true.
   allocate(grid%vertex_exact(size(grid%vertex),grid%system_size,max(1,grid%num_eval)), &
            stat=stat)
   if (stat /= 0) then
      call fatal("allocation failed for vertex_exact in init_grid")
      stop
   endif
endif

! define the edges

nedge = 0
do i=1,grid%nelem
   do j=1,EDGES_PER_ELEMENT
      neigh = grid%initial_neighbor(j,i)
      neighbmark = 1
      if (neigh /= BOUNDARY) then
         do k=1,EDGES_PER_ELEMENT
            if (grid%initial_neighbor(k,neigh) == i) then
               neighbmark = elemedge_bmark(k,neigh)
            endif
         end do
      endif
      if (neigh == BOUNDARY .or. neigh > i .or. elemedge_bmark(j,i) < 0 .or. &
          neighbmark < 0) then
         nedge = nedge + 1
         grid%element(i)%edge(j) = nedge
         grid%edge(nedge)%vertex(1) = &
            grid%element(i)%vertex(1+mod(j,EDGES_PER_ELEMENT))
         grid%edge(nedge)%vertex(2) = &
            grid%element(i)%vertex(1+mod(j+1,EDGES_PER_ELEMENT))
         if (grid%vertex(grid%edge(nedge)%vertex(2))%gid < &
             grid%vertex(grid%edge(nedge)%vertex(1))%gid) then
            k = grid%edge(nedge)%vertex(1)
            grid%edge(nedge)%vertex(1) = grid%edge(nedge)%vertex(2)
            grid%edge(nedge)%vertex(2) = k
         endif
         grid%edge(nedge)%bmark = elemedge_bmark(j,i)
         if (neigh == BOUNDARY .or. elemedge_bmark(j,i) < 0 .or. &
             neighbmark < 0) then
! point is not on curved boundary, but I'm only after bctype
            call bconds((grid%vertex(grid%edge(nedge)%vertex(1))%coord%x + &
                         grid%vertex(grid%edge(nedge)%vertex(2))%coord%x)/2, &
                        (grid%vertex(grid%edge(nedge)%vertex(1))%coord%y + &
                         grid%vertex(grid%edge(nedge)%vertex(2))%coord%y)/2, &
                        grid%edge(nedge)%bmark,bctype,bccoef,bcrhs)
            grid%edge_type(nedge,:) = bctype
            do k=1,grid%system_size
               if (grid%edge_type(nedge,k) == PERIODIC) then
                  if (grid%edge(nedge)%bmark < 0) then
                     grid%edge_type(nedge,k) = PERIODIC_SLAVE
                  else
                     grid%edge_type(nedge,k) = PERIODIC_MASTER
                  endif
               endif
            end do
         else
            grid%edge_type(nedge,:) = INTERIOR
         endif
         grid%edge(nedge)%assoc_elem = i
         grid%edge(nedge)%degree = degree
         if (degree > 1) then
            allocate(grid%edge(nedge)%solution(degree-1,grid%system_size,max(1,grid%num_eval)),&
                     stat=stat)
            if (stat /= 0) then
               call fatal("allocation failed for mates in init_grid")
               stop
            endif
            if (grid%have_true) then
               allocate(grid%edge(nedge)%exact(degree-1,grid%system_size,max(1,grid%num_eval)),&
                        stat=stat)
               if (stat /= 0) then
                  call fatal("allocation failed for mates in init_grid")
                  stop
               endif
               grid%edge(nedge)%exact = 0.0_my_real
            endif
            do k=1,grid%system_size
               if (grid%edge_type(nedge,k) == DIRICHLET) then
                  call edge_exact(grid,nedge,k,"d")
               else
                  grid%edge(nedge)%solution(:,k,:) = 0.0_my_real
               endif
            end do
         endif
      else
         do k=1,EDGES_PER_ELEMENT
            if (grid%initial_neighbor(k,neigh) == i) then
               grid%element(i)%edge(j) = grid%element(neigh)%edge(k)
               exit
            endif
         end do
      endif
   end do
end do

! Initial edge gids start at ceiling(nedge/3)
igid = nedge/3
if (3*igid /= nedge) igid= igid+1
do i=1,nedge
   grid%edge(i)%gid = igid + i-1
   call hash_insert(grid%edge(i)%gid,i,grid%edge_hash)
end do
grid%next_free_edge = nedge+1

grid%nedge = nedge
if (set_iown) then
   grid%nedge_own = grid%nedge
else
   grid%nedge_own = 0
endif
grid%dof = grid%nvert
if (degree >= 2) grid%dof = grid%dof + nedge*(degree-1)
if (degree >= 3) grid%dof = grid%dof + grid%nelem*(((degree-2)*(degree-1))/2)
grid%dof_own = grid%dof

! find the master edge for PERIODIC_SLAVE edges

do i=1,grid%nelem
   do j=1,EDGES_PER_ELEMENT
      edge = grid%element(i)%edge(j)
      if (any(grid%edge_type(edge,:) == PERIODIC_SLAVE)) then
         neigh = grid%initial_neighbor(j,i)
         v1s = grid%edge(edge)%vertex(1)
         v2s = grid%edge(edge)%vertex(2)
         do k=1,EDGES_PER_ELEMENT
            v1m = grid%edge(grid%element(neigh)%edge(k))%vertex(1)
            v2m = grid%edge(grid%element(neigh)%edge(k))%vertex(2)
            if ((v1m == grid%vertex(v1s)%next .and. v2m == grid%vertex(v2s)%next) .or. &
                (v2m == grid%vertex(v1s)%next .and. v1m == grid%vertex(v2s)%next)) then
               grid%edge(edge)%next = grid%element(neigh)%edge(k)
               exit
            endif
         end do
      endif
   end do
end do

grid%partition = partition
grid%nlev = 1
grid%errind_up2date = .false.

! smooth the triangle shapes

! TEMP071217 optional no smoothing; for the battery example

if (.not. dont_smooth) then
   call smooth_grid(grid,procs)
endif

! compute exact after smoothing, because vertices moved

if (grid%have_true) then
   do i=1,grid%nvert
      do k=1,max(1,grid%num_eval)
         do j=1,grid%system_size
            grid%vertex_exact(i,j,k) = trues(grid%vertex(i)%coord%x, &
                                             grid%vertex(i)%coord%y, &
                                             j,k)
         end do
      end do
   end do
endif

if (grid%have_true) then
   do i=1,grid%nelem
      do k=1,grid%system_size
         do j=1,EDGES_PER_ELEMENT
            if (grid%initial_neighbor(j,i) > i) cycle
            call edge_exact(grid,grid%element(i)%edge(j),k,"t")
         end do
         call elem_exact(grid,i,k,"t")
      end do
   end do
endif

! set the path through the initial grid for RTK.  If using Zoltan, Zoltan will
! reset it.

call init_path(grid)

deallocate(td%vert_coord,td%tri_edge,td%tri_vert,td%tri_neigh,td%edge_tri, &
           td%edge_vert,td%vert_tri,td%vert_edge,td%vert_bmark, &
           td%vert_bparam,td%edge_bmark,td%vert_master,stat=stat)

end subroutine init_grid

!          ---------------------
subroutine create_triangle_files(grid)
!          ---------------------

!----------------------------------------------------
! This routine creates the triangle data files for the domain given by the
! domain boundary subroutines.
! grid%max_blen is the maximum length of the line between adjacent boundary
! vertices.
! grid%triangle_files is the root of the created triangle data files.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
!----------------------------------------------------
! Local variables:

logical :: exists, opened
logical, allocatable :: singleton(:)
character(len=20) :: carea
integer :: total_npiece, k, astat, iounit, nhole, nvert, hole, vert, &
           edge, piece, first_vert
integer, allocatable :: npiece(:), end_hole(:), last_vert(:)
real(my_real) :: area, x1, y1, x2, y2
real(my_real), allocatable :: mid(:)
!----------------------------------------------------
! Begin executable code

! count the number of holes

nhole = 1
do while (boundary_npiece(nhole) > 0)
   nhole = nhole + 1
end do
nhole = nhole - 1
allocate(npiece(0:nhole),end_hole(-1:nhole),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in create_triangle_files")
   stop
endif

! get number of boundary pieces, the piece that ends each hole, and
! parameter limits

total_npiece = 0
end_hole(-1) = 0
do hole=0,nhole
   npiece(hole) = boundary_npiece(hole)
   total_npiece = total_npiece + npiece(hole)
   end_hole(hole) = end_hole(hole-1) + npiece(hole)
end do
allocate(mid(total_npiece),last_vert(total_npiece),singleton(0:total_npiece), &
         stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in create_triangle_files")
   stop
endif
mid = (grid%bp_start+grid%bp_finish)/2

! find an available unit and open file for writing

iounit = 11
do
   inquire(unit=iounit,exist=exists,opened=opened)
   if (exists .and. .not. opened) exit
   iounit = iounit + 1
end do
open(unit=iounit,file=trim(grid%triangle_files)//".poly",status="replace")

! count the number of vertices and edges that will be written, and identify
! boundary pieces that are just a single vertex

nvert = 0
singleton(0) = .false.
do piece=1,total_npiece
   singleton(piece) = grid%bp_start(piece)==grid%bp_finish(piece)
   call create_triangle_files_recur(piece,grid%bp_start(piece),mid(piece), &
                                    grid%max_blen,nvert,.not.singleton(piece-1))
   if (.not. singleton(piece)) then
      call create_triangle_files_recur(piece,mid(piece),grid%bp_finish(piece), &
                                       grid%max_blen,nvert,.true.)
   endif
end do

! write the vertex count line

write(iounit,"(I12,3I2)") nvert,2,0,1

! make sure the end of each piece is the beginning of the next, and that
! no part has a singleton at the end

do hole=0,nhole
   do piece=end_hole(hole-1)+1,end_hole(hole)
      if (piece == end_hole(hole)) then
         k = end_hole(hole-1) + 1
      else
         k = piece+1
      endif
      call boundary_point(piece,grid%bp_finish(piece),x1,y1)
      call boundary_point(k,grid%bp_start(k),x2,y2)
      if (sqrt((x1-x2)**2 + (y1-y2)**2) > 10000*epsilon(1.0_my_real)) then
         ierr = USER_INPUT_ERROR
         call fatal("End of boundary piece does not match beginning of next piece.", &
                    reallist=(/x1,y1,x2,y2/),intlist=(/piece,k/))
         stop
      endif
   end do
   if (singleton(end_hole(hole))) then
      ierr = USER_INPUT_ERROR
      call fatal("Final boundary piece of outer boundary or any hole cannot be a single point.", &
                 "Make that point be the first piece, instead.")
      stop
   endif
end do

! write the vertices
! when a piece is a singleton, don't do the second half and don't do the
! first point of the next piece

nvert = 0
do piece=1,total_npiece
   call create_triangle_files_recur(piece,grid%bp_start(piece),mid(piece), &
                                    grid%max_blen,nvert, &
                                    .not.singleton(piece-1),iounit)
   if (.not. singleton(piece)) then
      call create_triangle_files_recur(piece,mid(piece),grid%bp_finish(piece), &
                                       grid%max_blen,nvert,.true.,iounit)
   endif
   last_vert(piece) = nvert
end do

! write the edge count line

write(iounit,"(I12,I2)") nvert,1

! write the edges

hole = 0
first_vert = 1
piece = 1
do while(singleton(piece))
   piece = piece+1
end do
do vert=1,nvert
   edge = vert
   if (vert == last_vert(piece)) then
      if (piece == end_hole(hole)) then
         write(iounit,"(4I12)") edge,vert,first_vert,piece
         hole = hole+1
         first_vert = vert+1
      else
         write(iounit,"(4I12)") edge,vert,vert+1,piece
      endif
      piece = piece+1
      if (piece > total_npiece) exit
      do while(singleton(piece))
         piece = piece+1
      end do
   else
      write(iounit,"(4I12)") edge,vert,vert+1,piece
   endif
end do

! write the hole count line

write(iounit,"(I12)") nhole

! for each hole compute a point inside the hole as the average of the
! piece endpoints and midpoints.  this can fail if the hole is too concave.

do hole=1,nhole
   x2 = 0.0_my_real
   y2 = 0.0_my_real
   k = 0
   do piece=end_hole(hole-1)+1,end_hole(hole)
      call boundary_point(piece,grid%bp_start(piece),x1,y1)
      x2 = x2 + x1
      y2 = y2 + y1
      call boundary_point(piece,(grid%bp_start(piece)+grid%bp_finish(piece))/2,&
                          x1,y1)
      x2 = x2 + x1
      y2 = y2 + y1
      k = k+2
   end do
   x2 = x2/k
   y2 = y2/k
   write(iounit,"(I12,SS,1P,2E21.13E2)") hole,x2,y2
end do

close(iounit)

deallocate(mid,last_vert,singleton,npiece,end_hole,stat=astat)

! run Triangle to create the other triangle data files

if (grid%max_blen == huge(0.0_my_real)) then
   call my_system("triangle -pneqjQIY "//trim(grid%triangle_files)//".poly")
else
   area = grid%max_blen**2*sqrt(3.0_my_real)/4
   write(carea,"(f20.10)") area
   k = index(carea," ",back=.true.)
   call my_system("triangle -pneq28 -jQI -a"//carea(k+1:)//" "//trim(grid%triangle_files)//".poly")
endif

end subroutine create_triangle_files

!          ------
recursive subroutine create_triangle_files_recur(ipiece,s1,s2,max_blen,count, &
                                                inc_start,iounit)
!          ------

!----------------------------------------------------
! This routine increments count by the number of vertices and edges needed
! for boundary piece i between parameters s1 and s2 with maximum segment
! length max_blen, and optionally write the vertices to unit iounit if present.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: ipiece
real(my_real), intent(in) :: s1, s2, max_blen
integer, intent(inout) :: count
logical, intent(in) :: inc_start
integer, intent(in), optional :: iounit
!----------------------------------------------------
! Local variables:

real(my_real) :: x1, x2, y1, y2, sm
!----------------------------------------------------
! Begin executable code

! get the coordinates of the endpoints

call boundary_point(ipiece,s1,x1,y1)
call boundary_point(ipiece,s2,x2,y2)

if (sqrt((x2-x1)**2 + (y2-y1)**2) > max_blen) then

! if this segment is too long, bisect it

   sm = (s1+s2)/2
   call create_triangle_files_recur(ipiece,s1,sm,max_blen,count,inc_start,iounit)
   call create_triangle_files_recur(ipiece,sm,s2,max_blen,count,.true.,iounit)

! if it is short enough, increment the count and write the starting endpoint

else

   if (inc_start) then
      count = count+1
      if (present(iounit)) write(iounit,"(I12,SS,1P,2E21.13E2,I12)") count,x1,y1,ipiece
   endif

endif

end subroutine create_triangle_files_recur

!          ------------------
subroutine read_triangle_data(grid,td)
!          ------------------

!----------------------------------------------------
! This routine reads data from .node, .ele, .edge and .neigh files in the format
! of Jonathan Richard Shewchuk's mesh generation program "triangle".
!
! NOTE: I assume there are no comment lines before the end of the data.
!       Triangle 1.5 seems to obey this.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
type(triangle_data), intent(out) :: td
!----------------------------------------------------
! Local variables:

integer :: i, j, k, stat, td_dim, td_natt, td_nbm, bmark, vert, td_npt, tri, &
           v1, v2, v3, td_ntri2, td_nneigh, edge, end1, end2, iounit, nset
logical, allocatable :: used(:)
real(my_real) :: x, y
logical :: exists, opened
!----------------------------------------------------
! Begin executable code

! find an available i/o unit number

iounit = 11
do
   inquire(unit=iounit,exist=exists,opened=opened)
   if (exists .and. .not. opened) exit
   iounit = iounit + 1
end do

! read the node (vertex) data from the triangle data file

open(unit=iounit,file=trim(grid%triangle_files)//".node",status="old", &
     action="read",iostat=stat)
if (stat /= 0) then
   call fatal("open failed for file "//trim(grid%triangle_files)//".node", &
              "iostat is ",intlist=(/stat/))
   stop
endif

read(iounit,*) td%nvert, td_dim, td_natt, td_nbm
if (td_natt /= 0) then
   call fatal("number of attributes in .node file must be 0")
   stop
endif
allocate(used(td%nvert),stat=stat)
if (stat /= 0) then
   call fatal("memory allocation for vertices from .node file failed", &
              intlist=(/stat,td%nvert/))
   stop
endif

allocate(td%vert_tri(MAX_TD_VERT_NEIGH,td%nvert), &
         td%vert_edge(MAX_TD_VERT_NEIGH,td%nvert), &
         td%vert_bmark(td%nvert), td%vert_bparam(td%nvert), &
         td%vert_coord(td%nvert), &
         stat=stat)
if (stat /= 0) then
   call fatal("memory allocation for vertices from .node file failed", &
              intlist=(/stat,td%nvert/))
   stop
endif

if (td_nbm == 0) then
   call fatal("boundary markers are required in data from Triangle")
   stop
endif

do i=1,td%nvert
   read(iounit,*) vert,x,y,bmark
   if (vert < 1) then
      call fatal("vertices in .node file must be numbered starting at 1")
      stop
   endif
   if (vert > td%nvert) then
      call fatal("vertex number in .node file is larger than stated number of vertices", &
                 intlist=(/vert,td%nvert/))
      stop
   endif
   td%vert_coord(vert)%x = x
   td%vert_coord(vert)%y = y
   td%vert_bmark(vert) = bmark
   td%vert_bparam(vert) = find_point_bparam(bmark,x,y,grid)
end do

close(unit=iounit)

! read the element data from the .ele file

open(unit=iounit,file=trim(grid%triangle_files)//".ele",status="old", &
     action="read",iostat=stat)
if (stat /= 0) then
   call fatal("open failed for file "//trim(grid%triangle_files)//".ele", &
              "iostat is ",intlist=(/stat/))
   stop
endif

read(iounit,*) td%ntri, td_npt, td_natt
allocate(td%tri_edge(EDGES_PER_ELEMENT,td%ntri), &
         td%tri_vert(VERTICES_PER_ELEMENT,td%ntri), &
         stat=stat)
if (stat /= 0) then
   call fatal("memory allocation for triangles from .ele file failed", &
              intlist=(/stat,td%ntri/))
   stop
endif

used = .false.
do i=1,td%ntri
   read(iounit,*) tri,v1,v2,v3
   if (tri < 1) then
      call fatal("triangles in .ele file must be numbered starting at 1")
      stop
   endif
   if (tri > td%ntri) then
      call fatal("triangle number in .ele file is larger than stated number of triangles", &
                 intlist=(/tri,td%ntri/))
      stop
   endif
   td%tri_vert(1,tri) = v1
   td%tri_vert(2,tri) = v2
   td%tri_vert(3,tri) = v3
   used(v1) = .true.
   used(v2) = .true.
   used(v3) = .true.
end do

close(unit=iounit)

if (.not. all(used)) then
   ierr = USER_INPUT_ERROR
   call fatal("There are unused nodes in the .node file.", &
              "Use the -j flag when running triangle.")
   stop
endif

deallocate(used)

! read the neighbor data from the .neigh file

open(unit=iounit,file=trim(grid%triangle_files)//".neigh",status="old", &
     action="read",iostat=stat)
if (stat /= 0) then
   call fatal("open failed for file "//trim(grid%triangle_files)//".neigh", &
              "iostat is ",intlist=(/stat/))
   stop
endif

read(iounit,*) td_ntri2, td_nneigh
if (td_ntri2 /= td%ntri) then
   call fatal("number of triangles in .neigh file is not the same as number in .ele file", &
              intlist=(/td_ntri2,td%ntri/))
   stop
endif
allocate(td%tri_neigh(3,td%ntri),stat=stat)
if (stat /= 0) then
   call fatal("memory allocation for neighbors from .neigh file failed", &
              intlist=(/stat,td%ntri/))
   stop
endif

do i=1,td%ntri
   read(iounit,*) tri,v1,v2,v3
   if (tri < 1) then
      call fatal("triangles in .neigh file must be numbered starting at 1")
      stop
   endif
   if (tri > td%ntri) then
      call fatal("triangle number in .neigh file is larger than stated number of triangles", &
                 intlist=(/tri,td%ntri/))
      stop
   endif
   td%tri_neigh(1,tri) = v1
   td%tri_neigh(2,tri) = v2
   td%tri_neigh(3,tri) = v3
end do

close(unit=iounit)

! read the edge data from the triangle edge file

open(unit=iounit,file=trim(grid%triangle_files)//".edge",status="old", &
     action="read",iostat=stat)
if (stat /= 0) then
   call fatal("open failed for file "//trim(grid%triangle_files)//".edge", &
              "iostat is ",intlist=(/stat/))
   stop
endif

read(iounit,*) td%nedge, td_nbm
allocate(td%edge_tri(2,td%nedge), td%edge_vert(2,td%nedge), &
         td%edge_bmark(td%nedge), stat=stat)
if (stat /= 0) then
   call fatal("memory allocation for vertices from .edge file failed", &
              intlist=(/stat,td%nedge/))
   stop
endif

if (td_nbm == 0) then
   call fatal("boundary markers are required in data from Triangle")
   stop
endif

do i=1,td%nedge
   read(iounit,*) edge,end1,end2,bmark
   if (edge < 1) then
      call fatal("edges in .edge file must be numbered starting at 1")
      stop
   endif
   if (edge > td%nedge) then
      call fatal("edge number in .edge file is larger than stated number of edges", &
                 intlist=(/edge,td%nedge/))
      stop
   endif
   td%edge_vert(1,edge) = end1
   td%edge_vert(2,edge) = end2
   td%edge_bmark(edge) = bmark
end do

close(unit=iounit)

! NEW
! derive other components of triangle data

! set the triangle list for each vertex

td%vert_tri = -1
do i=1,td%ntri
   do j=1,3
      do k=1,MAX_TD_VERT_NEIGH
         if (td%vert_tri(k,td%tri_vert(j,i)) == -1) exit
      end do
      if (k == MAX_TD_VERT_NEIGH+1) then
         call fatal("too many neighbors of a vertex in triangle data")
         stop
      endif
      td%vert_tri(k,td%tri_vert(j,i)) = i
   end do
end do

! set the edge list for each vertex

td%vert_edge = -1
do i=1,td%nedge
   do j=1,2
      do k=1,MAX_TD_VERT_NEIGH
         if (td%vert_edge(k,td%edge_vert(j,i)) == -1) exit
      end do
      if (k == MAX_TD_VERT_NEIGH+1) then
         call fatal("too many neighbors of a vertex in triangle data")
         stop
      endif
      td%vert_edge(k,td%edge_vert(j,i)) = i
   end do
end do

! set the edge list of each triangle, and triangle list of each edge

td%tri_edge = -1
td%edge_tri = -1

! for each triangle
do i=1,td%ntri
   nset = 0
! for each vertex of the triangle
   do j=1,3
! search the edges of the vertex for any that contain another vertex of the
! triangle
! for each edge of this vertex
      do k=1,MAX_TD_VERT_NEIGH
         if (td%vert_edge(k,td%tri_vert(j,i)) == -1) exit
! if the first vertex of the edge is this vertex, see if the second vertex of
! the edge is a vertex of the triangle
         if (td%edge_vert(1,td%vert_edge(k,td%tri_vert(j,i))) == td%tri_vert(j,i)) then
            if (td%edge_vert(2,td%vert_edge(k,td%tri_vert(j,i))) == &
                td%tri_vert(1,i) .or. &
                td%edge_vert(2,td%vert_edge(k,td%tri_vert(j,i))) == &
                td%tri_vert(2,i) .or. &
                td%edge_vert(2,td%vert_edge(k,td%tri_vert(j,i))) == &
                td%tri_vert(3,i)) then
! if so make it an edge of this triangle, and make this triangle a triangle
! of that edge, unless it has already been set
               if (td%edge_tri(1,td%vert_edge(k,td%tri_vert(j,i))) /= i .and. &
                   td%edge_tri(2,td%vert_edge(k,td%tri_vert(j,i))) /= i) then
                  nset = nset + 1
                  td%tri_edge(nset,i) = td%vert_edge(k,td%tri_vert(j,i))
                  if (td%edge_tri(1,td%vert_edge(k,td%tri_vert(j,i))) == -1) then
                     td%edge_tri(1,td%vert_edge(k,td%tri_vert(j,i))) = i
                  elseif (td%edge_tri(2,td%vert_edge(k,td%tri_vert(j,i))) == -1) then
                     td%edge_tri(2,td%vert_edge(k,td%tri_vert(j,i))) = i
                  else
                     call fatal("too many triangles neighboring an edge in read_trianlge_data")
                     stop
                  endif
               endif
            endif
! if the second vertex of the edge is this vertex, see if the first vertex of
! the edge is a vertex of the triangle
         elseif (td%edge_vert(2,td%vert_edge(k,td%tri_vert(j,i))) == td%tri_vert(j,i)) then
            if (td%edge_vert(1,td%vert_edge(k,td%tri_vert(j,i))) == &
                td%tri_vert(1,i) .or. &
                td%edge_vert(1,td%vert_edge(k,td%tri_vert(j,i))) == &
                td%tri_vert(2,i) .or. &
                td%edge_vert(1,td%vert_edge(k,td%tri_vert(j,i))) == &
                td%tri_vert(3,i)) then
! if so make it an edge of this triangle, and make this triangle a triangle
! of that edge, unless it has already been set
               if (td%edge_tri(1,td%vert_edge(k,td%tri_vert(j,i))) /= i .and. &
                   td%edge_tri(2,td%vert_edge(k,td%tri_vert(j,i))) /= i) then
                  nset = nset + 1
                  td%tri_edge(nset,i) = td%vert_edge(k,td%tri_vert(j,i))
                  if (td%edge_tri(1,td%vert_edge(k,td%tri_vert(j,i))) == -1) then
                     td%edge_tri(1,td%vert_edge(k,td%tri_vert(j,i))) = i
                  elseif (td%edge_tri(2,td%vert_edge(k,td%tri_vert(j,i))) == -1) then
                     td%edge_tri(2,td%vert_edge(k,td%tri_vert(j,i))) = i
                  else
                     call fatal("too many triangles neighboring an edge in read_trianlge_data")
                     stop
                  endif
               endif
            endif
         endif
      end do
   end do
end do

! verify that all triangles have 3 edges and all edges have 2 triangles or
! are boundary

do i=1,td%ntri
   if (td%tri_edge(3,i) == -1) then
      call fatal("didn't assign 3 edges to all triangles in read_triangle_data")
      stop
   endif
end do
! TEMP must verify that bmark/=0 iff vertex is on boundary.  Might need to
! change the documentation
do i=1,td%nedge
   if (td%edge_tri(1,i) == -1 .or. &
       (td%edge_bmark(i) == 0 .and. td%edge_tri(2,i) == -1)) then
      call fatal("didn't assign 2 triangles or 1 triangle and boundary mark to all edges in read_triangle_data")
      stop
   endif
end do

end subroutine read_triangle_data

!        -----------------
function find_point_bparam(piece,x,y,grid)
!        -----------------

!----------------------------------------------------
! This routine finds the parameter for the point (x,y) on boundary piece piece.
!
! The point must lie on this piece of the boundary or the routine will get
! stuck.  It uses a secant root finder to find the root of
! (x(s)-x)**2 - (y(s)-y)**2
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: piece
real(my_real), intent(in) :: x,y
real(my_real) :: find_point_bparam
type(grid_type), intent(in) :: grid
!----------------------------------------------------
! Local variables:

integer :: nstart, i
real(my_real) :: s_start, s_finish, s, s0, s1, s2, f0, f1, f2, xtry, ytry, &
                 xymag
real(my_real), parameter :: small_enough = 100*epsilon(0.0_my_real)
!----------------------------------------------------
! Begin executable code

! no parameter if the domain is defined from triangle data files, or the
! point is interior

if (boundary_npiece(0) <= 0 .or. piece == 0) then
   find_point_bparam = 0.0_my_real
   return
endif

! check the starting and finishing points

s_start = grid%bp_start(piece)
s_finish = grid%bp_finish(piece)
call boundary_point(piece,s_start,xtry,ytry)
xymag = max(1.0_my_real,sqrt(x*x+y*y))
if (sqrt((x-xtry)**2+(y-ytry)**2)/xymag < small_enough) then
   find_point_bparam = s_start
   return
endif

call boundary_point(piece,s_finish,xtry,ytry)
if (sqrt((x-xtry)**2+(y-ytry)**2) < small_enough) then
   find_point_bparam = s_finish
   return
endif

! begin by going through nstart equally spaced parameters to find the one
! that minimizes the target function

nstart = 100
f0 = huge(0.0_my_real)
do i=1,nstart-1
   s = s_start + i*((s_finish-s_start)/nstart)
   call boundary_point(piece,s,xtry,ytry)
   if ((x-xtry)**2+(y-ytry)**2 < f0) then
      f0 = (x-xtry)**2+(y-ytry)**2
      s0 = s
   endif
end do

! if found one that's good enough, return it

if (f0 < small_enough) then
   find_point_bparam = s0
   return
endif

! get a second starting point

s1 = s0 + (s_finish-s_start)/(10*nstart)
call boundary_point(piece,s1,xtry,ytry)
f1 = (x-xtry)**2+(y-ytry)**2

! perform secant method iteration until we are close enough.  The method
! fails if this goes into an infinite loop or a secant iteration finds a
! parameter that is out of the range of this piece's parameters.

do
   s2 = s1 - f1*(s1-s0)/(f1-f0)
   if (s2 < s_start .or. s2 > s_finish) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("secant method failed in find_point_bparam")
      stop
   endif
   call boundary_point(piece,s2,xtry,ytry)
   f2 = (x-xtry)**2+(y-ytry)**2
   if (f2 < small_enough) then
      find_point_bparam = s2
      return
   endif
   s0 = s1; s1 = s2; f0 = f1; f1 = f2
end do

end function find_point_bparam

!          -------------
subroutine point_on_edge(grid,x1,y1,bmark1,bparam1,x2,y2,bmark2,bparam2,x3,y3, &
                         emark,f,x,y,bparam)
!          -------------

!----------------------------------------------------
! This routine determines the (x,y) coordinate and bparam for a point on the
! (possibly curved) edge [(x1,y1),(x2,y2)] of the triangle with that edge and
! vertex (x3,y3).  The bmarks and bparams of the endpoints, and edge bparam
! are also input.  The point is 0 < f < 1 of the way from point 1 to 2.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
real(my_real), intent(in) :: x1,y1,bparam1,x2,y2,bparam2,x3,y3
integer, intent(in) :: bmark1,bmark2,emark
real(my_real), intent(in) :: f
real(my_real), intent(out) :: x,y,bparam
!----------------------------------------------------
! Local variables:

integer :: piece
real(my_real) :: p1, p2, pm, temp, a1, a2, a3, xm, ym
!----------------------------------------------------
! Begin executable code

! if the boundary subroutines are not used or it is an interior edge,
! then return the midpoint and bparam doesn't matter

if (boundary_npiece(0) <= 0 .or. emark == 0) then
   x = (1-f)*x1 + f*x2
   y = (1-f)*y1 + f*y2
   bparam = 0
   return
endif

! determine parameters for the endpoints on piece emark

! first point is on same piece as edge
if (bmark1 == emark) then
   p1 = bparam1
! first point is a singleton.  A singleton cannot be the last piece of a
! hole, so we can look for the next larger index that is not a singleton.
! This will be emark iff bmark1 comes before emark.
elseif (bparam1 == grid%bp_start(bmark1) .and. &
        bparam1 == grid%bp_finish(bmark1)) then
   piece = bmark1+1
   do while (grid%bp_start(piece) == grid%bp_finish(piece))
      piece = piece + 1
   end do
   if (piece == emark) then
      p1 = grid%bp_start(emark)
   else
      p1 = grid%bp_finish(emark)
   endif
! first point is the beginning of it's piece, so also end of edge piece
elseif (bparam1 == grid%bp_start(bmark1)) then
   p1 = grid%bp_finish(emark)
! first point is the end of it's piece, so also the begining of edge piece
elseif (bparam1 == grid%bp_finish(bmark1)) then
   p1 = grid%bp_start(emark)
! failed to find what it is
else
   ierr = PHAML_INTERNAL_ERROR
   call fatal("Failed to find edge endpoint in edge piece.")
endif

if (bmark2 == emark) then
   p2 = bparam2
elseif (bparam2 == grid%bp_start(bmark2) .and. &
        bparam2 == grid%bp_finish(bmark2)) then
   piece = bmark2+1
   do while (grid%bp_start(piece) == grid%bp_finish(piece))
      piece = piece + 1
   end do
   if (piece == emark) then
      p2 = grid%bp_start(emark)
   else
      p2 = grid%bp_finish(emark)
   endif
elseif (bparam2 == grid%bp_start(bmark2)) then
   p2 = grid%bp_finish(emark)
elseif (bparam2 == grid%bp_finish(bmark2)) then
   p2 = grid%bp_start(emark)
else
   ierr = PHAML_INTERNAL_ERROR
   call fatal("Failed to find edge endpoint in edge piece.")
endif

! find the equation of the line that goes through the desired point of the edge
! and the opposite vertex, in the form a1*X + a2*Y + a3 = 0

xm = (1-f)*x1 + f*x2
ym = (1-f)*y1 + f*y2
a1 = y3-ym
a2 = xm-x3
a3 = x3*ym - y3*xm

! if vertex 1 in the line equation is positive, swap p1 and p2 so that
! p1 goes with the negative one

if (a1*x1+a2*y1+a3 > 0.0_my_real) then
   temp = p1
   p1 = p2
   p2 = temp
endif

! use a bisection root finder to find the point where the line intersects
! the boundary

do
   pm = (p1+p2)/2
   call boundary_point(emark,pm,xm,ym)
   if (abs(a1*xm+a2*ym+a3) < 100*epsilon(0.0_my_real)) then
      x = xm
      y = ym
      bparam = pm
      return
   elseif (a1*xm+a2*ym+a3 < 0.0_my_real) then
      p1 = pm
   else
      p2 = pm
   endif
end do

end subroutine point_on_edge

!          --------------
subroutine match_periodic(grid,td,NO_MATE)
!          --------------

!----------------------------------------------------
! This routine matches pairs of sides that are periodic bounday so that
! the vertices align.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
type(triangle_data), intent(inout) :: td
!----------------------------------------------------
! Local variables:

integer :: mark, astat, astat2, jerr, itemp, edge, i, m, edge1, edge2, vert1, &
           vert2, vert3, endvert(2,2), nseg(2), end1, end2, tri1, tri2, syssize
real(my_real) :: acclen1, acclen2, seglen1, seglen2, newacclen1, newacclen2, &
                 slope1, slope2, frac1, frac2, lenside(2)
integer :: first_edge(2), last_edge(2)
integer, pointer :: next_edge(:),prev_edge(:),assoc_tri(:)
logical, pointer :: swap(:)
logical :: rev(2), straight1, straight2
integer, intent(in) :: NO_MATE
!----------------------------------------------------

syssize = grid%system_size

! allocate space for linked lists

allocate(next_edge(td%nedge),prev_edge(td%nedge),assoc_tri(td%nedge), &
         stat=astat)

! swap indicates if the endpoints of an edge need to be visited in opposite
! order while traversing the boundary side

allocate(swap(td%nedge),stat=astat2)

if (astat /= 0 .or. astat2 /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in match_periodic")
   stop
endif

! one of the pair of periodic edges is required to have a negative bmark,
! while the matching edge is its absolute value.  Repeat for each
! negative bmark.

mark = -huge(0)
outer: do
   mark = minval(td%edge_bmark, mask=(td%edge_bmark > mark) )
   if (mark >= 0) exit

! for the boundary side with mark mark and its periodic companion with
! mark abs(mark) (aka -mark) create a linked list of consecutive boundary
! edges, and corresponding triangles.

   do m = mark, -mark, -2*mark
      if (m==mark) then
         i=1
      else
         i=2
      endif

      call list_mark(m,syssize,td,first_edge(i),last_edge(i),next_edge, &
                     prev_edge,assoc_tri,jerr) 

! if jerr=-1, then the negative mark was not used to designate a periodic side
! so nothing to do

      if (jerr == -1) cycle outer
   end do

! find the vertices at the ends of the boundary segments

   do i=1,2
      if (next_edge(first_edge(i)) == END_OF_LIST) then
         endvert(1,i) = td%edge_vert(1,first_edge(i))
      elseif (td%edge_vert(1,first_edge(i))== &
              td%edge_vert(1,next_edge(first_edge(i))) .or. &
              td%edge_vert(1,first_edge(i))== &
              td%edge_vert(2,next_edge(first_edge(i)))) then
         endvert(1,i) = td%edge_vert(2,first_edge(i))
      else
         endvert(1,i) = td%edge_vert(1,first_edge(i))
      endif
      if (prev_edge(last_edge(i)) == END_OF_LIST) then
         endvert(2,i) = td%edge_vert(2,last_edge(i))
      elseif (td%edge_vert(1,last_edge(i))== &
              td%edge_vert(1,prev_edge(last_edge(i))) .or. &
              td%edge_vert(1,last_edge(i))== &
              td%edge_vert(2,prev_edge(last_edge(i)))) then
         endvert(2,i) = td%edge_vert(2,last_edge(i))
      else
         endvert(2,i) = td%edge_vert(1,last_edge(i))
      endif
   end do

! attempt to find the direction to traverse the boundary segments so that the
! periodic conditions match
! TEMP this is not fool proof

   rev = .false.

! if the first and last x coordinates differ for both of them, then go from
! the small x to the large x

   if (td%vert_coord(endvert(1,1))%x /= td%vert_coord(endvert(2,1))%x .and. &
       td%vert_coord(endvert(1,2))%x /= td%vert_coord(endvert(2,2))%x) then
      if (td%vert_coord(endvert(1,1))%x > td%vert_coord(endvert(2,1))%x) then
         rev(1) = .true.
      endif
      if (td%vert_coord(endvert(1,2))%x > td%vert_coord(endvert(2,2))%x) then
         rev(2) = .true.
      endif

! if that failed, same idea with the y coordinate

   elseif (td%vert_coord(endvert(1,1))%y /= td%vert_coord(endvert(2,1))%y .and. &
           td%vert_coord(endvert(1,2))%y /= td%vert_coord(endvert(2,2))%y) then
      if (td%vert_coord(endvert(1,1))%y > td%vert_coord(endvert(2,1))%y) then
         rev(1) = .true.
      endif
      if (td%vert_coord(endvert(1,2))%y > td%vert_coord(endvert(2,2))%y) then
         rev(2) = .true.
      endif

! if that failed, put the end closer to the origin first

   else
      if (td%vert_coord(endvert(1,1))%x**2+td%vert_coord(endvert(1,1))%y**2 > &
          td%vert_coord(endvert(2,1))%x**2+td%vert_coord(endvert(2,1))%y**2) then
         rev(1) = .true.
      endif
      if (td%vert_coord(endvert(1,2))%x**2+td%vert_coord(endvert(1,2))%y**2 > &
          td%vert_coord(endvert(2,2))%x**2+td%vert_coord(endvert(2,2))%y**2) then
         rev(2) = .true.
      endif

   endif

! if either list wants to be reversed, reverse it

   do i=1,2
      if (rev(i)) then
         itemp = endvert(1,i)
         endvert(1,i) = endvert(2,i)
         endvert(2,i) = itemp
         edge = first_edge(i)
         itemp = first_edge(i)
         first_edge(i) = last_edge(i)
         last_edge(i) = itemp
         do
            if (edge == END_OF_LIST) exit
            itemp = next_edge(edge)
            next_edge(edge) = prev_edge(edge)
            prev_edge(edge) = itemp
            edge = itemp
         end do
      endif
   end do

! traverse each list to determine whether the first or second vertex of each
! edge matches the preceeding edge

   do i=1,2

! check starting point

      if (td%edge_vert(1,first_edge(i)) == endvert(1,i)) then
         swap(first_edge(i)) = .false.
      elseif (td%edge_vert(2,first_edge(i)) == endvert(1,i)) then
         swap(first_edge(i)) = .true.
      else
         ierr = PHAML_INTERNAL_ERROR
         call fatal("first vertex does not match endvert")
         stop
      endif

! traverse the list

      edge = first_edge(i)
      do
         if (edge == last_edge(i)) exit
         if (swap(edge)) then
            if (td%edge_vert(1,next_edge(edge)) == td%edge_vert(1,edge)) then
               swap(next_edge(edge)) = .false.
            elseif (td%edge_vert(2,next_edge(edge))==td%edge_vert(1,edge)) then
               swap(next_edge(edge)) = .true.
            else
               ierr = PHAML_INTERNAL_ERROR
               call fatal("first vertex does not match either vertex in next edge")
               stop
            endif
         else
            if (td%edge_vert(1,next_edge(edge)) == td%edge_vert(2,edge)) then
               swap(next_edge(edge)) = .false.
            elseif (td%edge_vert(2,next_edge(edge))==td%edge_vert(2,edge)) then
               swap(next_edge(edge)) = .true.
            else
               ierr = PHAML_INTERNAL_ERROR
               call fatal("second vertex does not match either vertex in next edge")
               stop
            endif
         endif

         edge = next_edge(edge)
      end do

! verify the last point

      if ((swap(edge).and.td%edge_vert(1,edge)/=endvert(2,i)) .or. &
          (.not.swap(edge).and.td%edge_vert(2,edge)/=endvert(2,i))) then
         ierr = PHAML_INTERNAL_ERROR
         call fatal("last point did not match second endvert")
         stop
      endif

   end do ! next boundary side

! count the number of edge segments in each list and compute the total
! length of each boundary side

   do i=1,2
      nseg(i) = 0
      lenside(i) = 0.0_my_real
      edge = first_edge(i)
      do
         if (edge == END_OF_LIST) exit
         end1 = td%edge_vert(1,edge)
         end2 = td%edge_vert(2,edge)
         lenside(i) = lenside(i) + &
                      (td%vert_coord(end1)%x - td%vert_coord(end2)%x)**2 + &
                      (td%vert_coord(end1)%y - td%vert_coord(end2)%y)**2
         nseg(i) = nseg(i) + 1
         edge = next_edge(edge)
      end do
   end do

! traverse the two boundary sides lining up the vertices according to
! fraction of total length.  If the number of segments are different,
! insert a new vertex when there is a large difference in the accumulated
! length

   edge1 = first_edge(1)
   edge2 = first_edge(2)
   acclen1 = 0.0_my_real
   acclen2 = 0.0_my_real

   do

! check for end of lists

      if (edge1 == END_OF_LIST) then
         if (edge2 == END_OF_LIST) then
            exit
         else
            ierr = PHAML_INTERNAL_ERROR
            call fatal("did not reach end of lists at same time when matching periodic boundaries")
            stop
         endif
      elseif (edge2 == END_OF_LIST) then
         ierr = PHAML_INTERNAL_ERROR
         call fatal("did not reach end of lists at same time when matching periodic boundaries")
         stop
      endif

! determine segment scaled lengths and new accumulations

      end1 = td%edge_vert(1,edge1)
      end2 = td%edge_vert(2,edge1)
      seglen1 = ((td%vert_coord(end1)%x - td%vert_coord(end2)%x)**2 + &
                 (td%vert_coord(end1)%y - td%vert_coord(end2)%y)**2)/lenside(1)
      newacclen1 = acclen1 + seglen1
      end1 = td%edge_vert(1,edge2)
      end2 = td%edge_vert(2,edge2)
      seglen2 = ((td%vert_coord(end1)%x - td%vert_coord(end2)%x)**2 + &
                 (td%vert_coord(end1)%y - td%vert_coord(end2)%y)**2)/lenside(2)
      newacclen2 = acclen2 + seglen2

! if no difference in the new accumulations, take it

      if (abs(newacclen1-newacclen2) < 10*epsilon(0.0_my_real)) then
         edge1 = next_edge(edge1)
         edge2 = next_edge(edge2)
         acclen1 = newacclen1
         acclen2 = newacclen2
         cycle
      endif

! if both lists are on the last section, take it.  If only one is on its
! last section, insert the required number of vertices in the other one

      if (next_edge(edge1) == END_OF_LIST) then
         if (next_edge(edge2) == END_OF_LIST) then
            edge1 = next_edge(edge1)
            edge2 = next_edge(edge2)
            acclen1 = newacclen1
            acclen2 = newacclen2
            cycle
         else
            call break_edge(grid,1,edge1,assoc_tri(edge1),nseg(2)-nseg(1),m, &
                            td,first_edge(1),last_edge(1),next_edge,prev_edge, &
                            assoc_tri,swap,NO_MATE)
            nseg(1) = nseg(2)
            cycle
         endif
      elseif (next_edge(edge2) == END_OF_LIST) then
         call break_edge(grid,2,edge2,assoc_tri(edge2),nseg(1)-nseg(2),m,td, &
                         first_edge(2),last_edge(2),next_edge,prev_edge, &
                         assoc_tri,swap,NO_MATE)
         nseg(2) = nseg(1)
         cycle
      endif

! note whether the current sections and the sections that follow form
! a straight line

      end1 = td%edge_vert(1,edge1)
      end2 = td%edge_vert(2,edge1)
      if (td%vert_coord(end2)%x == td%vert_coord(end1)%x) then
         slope1 = huge(0.0_my_real)/2
      else
         slope1 = (td%vert_coord(end2)%y-td%vert_coord(end1)%y) / &
                  (td%vert_coord(end2)%x-td%vert_coord(end1)%x)
      endif
      end1 = end2
      end2 = td%edge_vert(1,next_edge(edge1))
      if (td%vert_coord(end2)%x == td%vert_coord(end1)%x) then
         slope2 = huge(0.0_my_real)/2
      else
         slope2 = (td%vert_coord(end2)%y-td%vert_coord(end1)%y) / &
                  (td%vert_coord(end2)%x-td%vert_coord(end1)%x)
      endif
      straight1 = abs(slope1-slope2)<100*epsilon(0.0_my_real)

      end1 = td%edge_vert(1,edge2)
      end2 = td%edge_vert(2,edge2)
      if (td%vert_coord(end2)%x == td%vert_coord(end1)%x) then
         slope1 = huge(0.0_my_real)/2
      else
         slope1 = (td%vert_coord(end2)%y-td%vert_coord(end1)%y) / &
                  (td%vert_coord(end2)%x-td%vert_coord(end1)%x)
      endif
      end1 = end2
      end2 = td%edge_vert(1,next_edge(edge2))
      if (td%vert_coord(end2)%x == td%vert_coord(end1)%x) then
         slope2 = huge(0.0_my_real)/2
      else
         slope2 = (td%vert_coord(end2)%y-td%vert_coord(end1)%y) / &
                  (td%vert_coord(end2)%x-td%vert_coord(end1)%x)
      endif
      straight2 = abs(slope1-slope2)<100*epsilon(0.0_my_real)

! if the new accumulated lengths are well separated and the number of
! segments is not equal, then insert a new vertex in the longer one if
! the longer one has fewer segments

      if (abs(newacclen1-newacclen2) > max(seglen1,seglen2)/2 .and. &
          nseg(1) /= nseg(2)) then

         if (nseg(1) > nseg(2) .and. newacclen2 > newacclen1) then
            call break_edge(grid,2,edge2,assoc_tri(edge2),1,m,td, &
                            first_edge(2),last_edge(2),next_edge,prev_edge, &
                            assoc_tri,swap,NO_MATE)
            nseg(2) = nseg(2) + 1
            cycle
         elseif (nseg(2) > nseg(1) .and. newacclen1 > newacclen2) then
            call break_edge(grid,1,edge1,assoc_tri(edge1),1,m,td, &
                            first_edge(1),last_edge(1),next_edge,prev_edge, &
                            assoc_tri,swap,NO_MATE)
            nseg(1) = nseg(1) + 1
            cycle
         endif
      endif

! otherwise, move the vertices at the end of the current segments to agree,
! depending on which are on straight lines

! TEMP is there a possibility of degenerating triangles here?

      if (straight1) then
         if (straight2) then
            newacclen1 = (newacclen1+newacclen2)/2
            newacclen2 = newacclen1
            frac1 = (newacclen1-acclen1)/seglen1
            frac2 = (newacclen2-acclen2)/seglen2
            seglen1 = newacclen1 - acclen1
            seglen2 = newacclen2 - acclen2
         else
            newacclen1 = newacclen2
            frac1 = (newacclen1-acclen1)/seglen1
            seglen1 = newacclen1 - acclen1
            frac2 = 1
         endif
      elseif (straight2) then
         newacclen2 = newacclen1
         frac2 = (newacclen2-acclen2)/seglen2
         seglen2 = newacclen2 - acclen2
         frac1 = 1
      endif

      tri1 = assoc_tri(edge1)
      do i=1,3
         if (td%tri_vert(i,tri1) /= td%edge_vert(1,edge1) .and. &
             td%tri_vert(i,tri1) /= td%edge_vert(2,edge1)) then
            vert3 = td%tri_vert(i,tri1)
         endif
      end do
      call point_on_edge(grid, &
                         td%vert_coord(td%edge_vert(1,edge1))%x, &
                         td%vert_coord(td%edge_vert(1,edge1))%y, &
                         td%vert_bmark(td%edge_vert(1,edge1)), &
                         td%vert_bparam(td%edge_vert(1,edge1)), &
                         td%vert_coord(td%edge_vert(2,edge1))%x, &
                         td%vert_coord(td%edge_vert(2,edge1))%y, &
                         td%vert_bmark(td%edge_vert(2,edge1)), &
                         td%vert_bparam(td%edge_vert(2,edge1)), &
                         td%vert_coord(vert3)%x, td%vert_coord(vert3)%y, &
                         td%edge_bmark(edge1), frac1, &
                         td%vert_coord(td%edge_vert(2,edge1))%x, &
                         td%vert_coord(td%edge_vert(2,edge1))%y, &
                         td%vert_bparam(td%edge_vert(2,edge1)))
      tri1 = assoc_tri(edge2)
      do i=1,3
         if (td%tri_vert(i,tri1) /= td%edge_vert(1,edge2) .and. &
             td%tri_vert(i,tri1) /= td%edge_vert(2,edge2)) vert3 = td%tri_vert(i,tri1)
      end do
      call point_on_edge(grid, &
                         td%vert_coord(td%edge_vert(1,edge2))%x, &
                         td%vert_coord(td%edge_vert(1,edge2))%y, &
                         td%vert_bmark(td%edge_vert(1,edge2)), &
                         td%vert_bparam(td%edge_vert(1,edge2)), &
                         td%vert_coord(td%edge_vert(2,edge2))%x, &
                         td%vert_coord(td%edge_vert(2,edge2))%y, &
                         td%vert_bmark(td%edge_vert(2,edge2)), &
                         td%vert_bparam(td%edge_vert(2,edge2)), &
                         td%vert_coord(vert3)%x, td%vert_coord(vert3)%y, &
                         td%edge_bmark(edge2), frac2, &
                         td%vert_coord(td%edge_vert(2,edge2))%x, &
                         td%vert_coord(td%edge_vert(2,edge2))%y, &
                         td%vert_bparam(td%edge_vert(2,edge2)))

      edge1 = next_edge(edge1)
      edge2 = next_edge(edge2)
      acclen1 = newacclen1
      acclen2 = newacclen2

   end do ! next segment of both sides

! traverse the edge lists to make the corresponding elements neighbors

   edge1 = first_edge(1)
   edge2 = first_edge(2)

   if (swap(edge1)) then
      vert1 = td%edge_vert(2,edge1)
   else
      vert1 = td%edge_vert(1,edge1)
   endif
   if (swap(edge2)) then
      vert2 = td%edge_vert(2,edge2)
   else
      vert2 = td%edge_vert(1,edge2)
   endif
   if (td%edge_bmark(edge1) < 0) td%vert_master(vert1) = vert2
   if (td%edge_bmark(edge2) < 0) td%vert_master(vert2) = vert1

   do
      if (edge1 == END_OF_LIST) then
         if (edge2 == END_OF_LIST) then
            exit
         else
            ierr = PHAML_INTERNAL_ERROR
            call fatal("did not reach end of lists at same time when matching periodic boundaries")
            stop
         endif
      elseif (edge2 == END_OF_LIST) then
         ierr = PHAML_INTERNAL_ERROR
         call fatal("did not reach end of lists at same time when matching periodic boundaries")
         stop
      endif

      tri1 = assoc_tri(edge1)
      tri2 = assoc_tri(edge2)
      do i=1,3
         if (td%tri_vert(i,tri1) /= td%edge_vert(1,edge1) .and. &
             td%tri_vert(i,tri1) /= td%edge_vert(2,edge1)) then
            td%tri_neigh(i,tri1) = tri2
            exit
         endif
      end do
      do i=1,3
         if (td%tri_vert(i,tri2) /= td%edge_vert(1,edge2) .and. &
             td%tri_vert(i,tri2) /= td%edge_vert(2,edge2)) then
            td%tri_neigh(i,tri2) = tri1
            exit
         endif
      end do

      if (swap(edge1)) then
         vert1 = td%edge_vert(1,edge1)
      else
         vert1 = td%edge_vert(2,edge1)
      endif
      if (swap(edge2)) then
         vert2 = td%edge_vert(1,edge2)
      else
         vert2 = td%edge_vert(2,edge2)
      endif
      if (td%edge_bmark(edge1) < 0) td%vert_master(vert1) = vert2
      if (td%edge_bmark(edge2) < 0) td%vert_master(vert2) = vert1

      edge1 = next_edge(edge1)
      edge2 = next_edge(edge2)
   end do

end do outer ! next negative mark

deallocate(swap,next_edge,prev_edge,assoc_tri,stat=astat)

end subroutine match_periodic

!          ---------
subroutine list_mark(m,syssize,td,first_edge,last_edge,next_edge,prev_edge, &
                     assoc_tri,jerr)
!          ---------

!----------------------------------------------------
! This routine creates a linked list of edges that have bmark m, and the
! triangles containing those edges.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: m,syssize
type(triangle_data), intent(in) :: td
integer, intent(out) :: first_edge, last_edge
integer, intent(out) :: next_edge(:),prev_edge(:),assoc_tri(:)
integer, intent(out) :: jerr
!----------------------------------------------------
! Local variables:

integer :: tri, edge, vert, itype(syssize), tri0, edge0, dir, &
           vert1, vert2, i, j, k, newedge, newtri
real(my_real) :: c(syssize,syssize), rs(syssize)
!----------------------------------------------------
! Begin executable code

jerr = 0

! find some edge with bmark m

first_edge = END_OF_LIST
do edge=1,td%nedge
   if (td%edge_bmark(edge) == m) then
      first_edge = edge
      last_edge = edge
      next_edge(edge) = END_OF_LIST
      prev_edge(edge) = END_OF_LIST
      exit
   endif
end do

if (first_edge == END_OF_LIST) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("failed to find an edge with bmark ",intlist=(/m/))
   stop
endif

! find a triangle containing this edge

vert1 = td%edge_vert(1,edge)
vert2 = td%edge_vert(2,edge)
assoc_tri(edge) = -1
outer: do i=1,td%ntri
   do j=1,3
      if (td%tri_vert(j,i) == vert1) then
         do k=1,3
            if (td%tri_vert(k,i) == vert2) then
               assoc_tri(edge) = i
               exit outer
            endif
         end do
      end if
   end do
end do outer

if (assoc_tri(edge) == -1) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("failed to find a triangle with the edge")
   stop
endif

! check the boundary condition type for this edge.  If it is not PERIODIC and
! bmark is negative, then they simply used a negative bmark.  But if it is
! not PERIODIC and bmark is positive, then they did not designate the matching
! of periodic sides with +-bmark

! point may not be on curved boundary, but I'm only after itype
call bconds((td%vert_coord(td%edge_vert(1,edge))%x + &
             td%vert_coord(td%edge_vert(2,edge))%x)/2, &
            (td%vert_coord(td%edge_vert(1,edge))%y + &
             td%vert_coord(td%edge_vert(2,edge))%y)/2, &
            m,itype,c,rs)
if (all(itype /= PERIODIC)) then
   if (m < 0) then
      jerr = -1
      return
   else
      ierr = USER_INPUT_ERROR
      call fatal("boundary side with given bmark does not have PERIODIC", &
                "boundary conditions but the negative of it does",intlist=(/m/))
      stop
   endif
endif

! traverse the boundary in both directions from the found edge until
! the bmark changes, building a linked list of edges and triangles

tri0 = assoc_tri(edge)
edge0 = edge

do dir=1,2
   tri = tri0
   edge = edge0
   vert = td%edge_vert(dir,edge0)

   do

! find the other boundary edge that contains vert.

      do newedge=1,td%nedge
         if (td%edge_bmark(newedge) == 0) cycle ! interior edge
         if (newedge == edge) cycle ! same edge
         if (td%edge_vert(1,newedge) == vert .or. &
             td%edge_vert(2,newedge) == vert) exit
      end do

      if (newedge > td%nedge) then
         ierr = PHAML_INTERNAL_ERROR
         call fatal("couldn't find second edge with given vertex in list_mark",&
                    intlist=(/vert/))
         stop
      endif

! if bmark is m, put this on the list and continue.  if not, we have passed
! this end of the boundary segment with bmark m

      if (td%edge_bmark(newedge) == m) then
         if (dir == 1) then
            prev_edge(edge) = newedge
            next_edge(newedge) = edge
            prev_edge(newedge) = END_OF_LIST
            first_edge = newedge
         else
            next_edge(edge) = newedge
            prev_edge(newedge) = edge
            next_edge(newedge) = END_OF_LIST
            last_edge = newedge
         endif
         vert1 = td%edge_vert(1,newedge)
         vert2 = td%edge_vert(2,newedge)
         assoc_tri(newedge) = -1
outer3:  do i=1,td%ntri
            do j=1,3
               if (td%tri_vert(j,i) == vert1) then
                  do k=1,3
                     if (td%tri_vert(k,i) == vert2) then
                        newtri = i
                        exit outer3
                     endif
                  end do
               end if
            end do
         end do outer3
         assoc_tri(newedge) = newtri
         tri = newtri
         edge = newedge
         if (td%edge_vert(1,edge) == vert) then
            vert = td%edge_vert(2,edge)
         else
            vert = td%edge_vert(1,edge)
         endif
      else
         exit
      endif
   end do
end do

end subroutine list_mark

!          ----------
subroutine break_edge(grid,bside,edge,tri,nbreak,m,td,first_edge,last_edge, &
                      next_edge,prev_edge,assoc_tri,swap,NO_MATE)
!          ----------

!----------------------------------------------------
! This routine breaks an edge into nbreak equal pieces by refining triangle
! tri into nbreak+1 triangles by dividing from the vertex opposite edge edge

! TEMP some of the components of td (tri_edge, edge_tri, vert_tri, vert_edge)
! do not get redefined in here.  I'm not sure if that matters.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: bside, edge, tri, nbreak, m
type(triangle_data), intent(inout) :: td
integer, intent(inout) :: first_edge, last_edge
integer, pointer :: next_edge(:), prev_edge(:), assoc_tri(:)
logical, pointer :: swap(:)
integer, intent(in) :: NO_MATE
!----------------------------------------------------
! Local variables:

integer :: i, new_nvert, new_ntri, new_nedge, oppvert, astat, vert, last, j, &
           lastneigh, lastloc
real(my_real) :: f
type(point), pointer :: new_vert_coord(:)
integer, pointer :: new_tri_edge(:,:), new_tri_vert(:,:), new_tri_neigh(:,:), new_edge_tri(:,:), new_edge_vert(:,:), &
                    new_vert_tri(:,:), new_vert_edge(:,:), new_vert_bmark(:), new_edge_bmark(:), new_next(:), &
                    new_prev(:), new_assoc(:), new_vert_master(:)
real(my_real), pointer :: new_vert_bparam(:)
logical, pointer :: new_swap(:)
real(my_real) :: xtemp, ytemp
!----------------------------------------------------
! Begin executable code

! reallocate space for new entities

new_nvert = td%nvert + nbreak
new_ntri = td%ntri + nbreak
new_nedge = td%nedge + 2*nbreak

allocate(new_vert_coord(new_nvert), new_tri_edge(3,new_ntri), &
         new_tri_vert(3,new_ntri), new_tri_neigh(3,new_ntri), &
         new_edge_tri(2,new_nedge), new_edge_vert(2,new_nedge), &
         new_vert_bmark(new_nvert), new_vert_tri(MAX_TD_VERT_NEIGH,new_nvert), &
         new_vert_edge(MAX_TD_VERT_NEIGH,new_nvert), &
         new_vert_bparam(new_nvert), new_edge_bmark(new_nedge), &
         new_next(new_nedge), new_prev(new_nedge), new_assoc(new_nedge), &
         new_swap(new_nedge), new_vert_master(new_nvert),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in break_edge")
   stop
endif

! copy existing data to new location

new_vert_coord(1:td%nvert) = td%vert_coord
new_tri_edge(:,1:td%ntri) = td%tri_edge
new_tri_vert(:,1:td%ntri) = td%tri_vert
new_tri_neigh(:,1:td%ntri) = td%tri_neigh
new_edge_tri(:,1:td%nedge) = td%edge_tri
new_edge_vert(:,1:td%nedge) = td%edge_vert
new_vert_tri(:,1:td%nvert) = td%vert_tri
new_vert_edge(:,1:td%nvert) = td%vert_edge
new_vert_bmark(1:td%nvert) = td%vert_bmark
new_vert_bparam(1:td%nvert) = td%vert_bparam
new_edge_bmark(1:td%nedge) = td%edge_bmark
new_next(1:td%nedge) = next_edge
new_prev(1:td%nedge) = prev_edge
new_assoc(1:td%nedge) = assoc_tri
new_vert_master(1:td%nvert) = td%vert_master
new_vert_master(td%nvert+1:new_nvert) = NO_MATE
new_swap(1:td%nedge) = swap

! free old memory

deallocate(td%vert_coord,td%tri_edge,td%tri_vert,td%tri_neigh,td%edge_tri, &
           td%edge_vert,td%vert_tri,td%vert_edge,td%vert_bmark, &
           td%vert_bparam,td%edge_bmark,next_edge,prev_edge,assoc_tri, &
           td%vert_master,swap,stat=astat)
if (astat/=0) then
   call warning("deallocation failed in break_edge")
endif

! set standard names to new memory

td%vert_coord => new_vert_coord
td%tri_edge => new_tri_edge
td%tri_vert => new_tri_vert
td%tri_neigh => new_tri_neigh
td%edge_tri => new_edge_tri
td%edge_vert => new_edge_vert
td%vert_tri => new_vert_tri
td%vert_edge => new_vert_edge
td%vert_bmark => new_vert_bmark
td%vert_bparam => new_vert_bparam
td%edge_bmark => new_edge_bmark
next_edge => new_next
prev_edge => new_prev
assoc_tri => new_assoc
td%vert_master => new_vert_master
swap => new_swap

! find the opposite vertex in the associated triangle

do vert=1,3
   if (td%tri_vert(vert,tri) /= td%edge_vert(1,edge) .and. &
       td%tri_vert(vert,tri) /= td%edge_vert(2,edge)) exit
end do
oppvert = td%tri_vert(vert,tri)

! define the new vertices

do i=1,nbreak
   f = i/(nbreak+1.0_my_real)
   if (swap(edge)) f = 1-f
   call point_on_edge(grid, &
                      td%vert_coord(td%edge_vert(1,edge))%x, &
                      td%vert_coord(td%edge_vert(1,edge))%y, &
                      td%vert_bmark(td%edge_vert(1,edge)), &
                      td%vert_bparam(td%edge_vert(1,edge)), &
                      td%vert_coord(td%edge_vert(2,edge))%x, &
                      td%vert_coord(td%edge_vert(2,edge))%y, &
                      td%vert_bmark(td%edge_vert(2,edge)), &
                      td%vert_bparam(td%edge_vert(2,edge)), &
                      td%vert_coord(oppvert)%x, td%vert_coord(oppvert)%y, m,f, &
                      xtemp, ytemp, &
                      td%vert_bparam(td%nvert+i))
   td%vert_coord(td%nvert+i)%x = xtemp
   td%vert_coord(td%nvert+i)%y = ytemp
   td%vert_bmark(td%nvert+i) = m
end do

! define the new edges that replace the edge being broken.  The first one
! overwrites the existing one.

if (swap(edge)) then
   last = td%edge_vert(1,edge)
   td%edge_vert(1,edge) = td%nvert + 1
else
   last = td%edge_vert(2,edge)
   td%edge_vert(2,edge) = td%nvert + 1
endif

do i=1,nbreak-1
   td%edge_vert(1,td%nedge+i) = td%nvert + i
   td%edge_vert(2,td%nedge+i) = td%nvert + i + 1
   swap(td%nedge+i) = .false.
   td%edge_bmark(td%nedge+i) = m
end do

td%edge_vert(1,td%nedge+nbreak) = td%nvert + nbreak
td%edge_vert(2,td%nedge+nbreak) = last
swap(td%nedge+nbreak) = .false.
td%edge_bmark(td%nedge+nbreak) = m

! define the edges that slice through the triangle being broken

do i=1,nbreak
   td%edge_vert(1,td%nedge+nbreak+i) = oppvert
   td%edge_vert(2,td%nedge+nbreak+i) = td%nvert + i
   swap(td%nedge+nbreak+i) = .false.
   td%edge_bmark(td%nedge+nbreak+i) = m
end do

! define the new triangles; the first one overwrites the existing one

do i=1,3
   if (td%tri_vert(i,tri) == last) td%tri_vert(i,tri) = td%nvert + 1
end do
do j=1,nbreak
   td%tri_vert(1,td%ntri+j) = oppvert
   td%tri_vert(2,td%ntri+j) = td%nvert + j
   td%tri_vert(3,td%ntri+j) = td%nvert + j + 1
end do
td%tri_vert(3,td%ntri+nbreak) = last

! set the new neighbors

lastneigh = -1
outer1: do i=1,3
   lastneigh = td%tri_neigh(i,tri)
   if (lastneigh <= 0) cycle
   do j=1,3
      lastloc = j
      if (td%tri_vert(j,lastneigh) == last) exit outer1
   enddo
end do outer1

if (i > 3) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("didn't find last neighbor")
   stop
endif

td%tri_neigh(lastloc,tri) = td%ntri + 1
do i=1,nbreak
   td%tri_neigh(1,td%ntri+i) = -1
   td%tri_neigh(2,td%ntri+i) = td%ntri + i + 1
   td%tri_neigh(3,td%ntri+i) = td%ntri + i - 1
end do
td%tri_neigh(3,td%ntri+1) = tri
td%tri_neigh(2,td%ntri+nbreak) = lastneigh

! set the linked lists and associated triangle

last = next_edge(edge)
do i=1,nbreak
   prev_edge(td%nedge+i) = td%nedge + i - 1
   next_edge(td%nedge+i) = td%nedge + i + 1
   assoc_tri(td%nedge+i) = tri
end do
next_edge(edge) = td%nedge + 1
prev_edge(td%nedge+1) = edge
next_edge(td%nedge+nbreak) = last
if (last /= END_OF_LIST) then
   prev_edge(last) = td%nedge+nbreak
else
   last_edge = td%nedge+nbreak
endif

! set new counts

td%nvert = new_nvert
td%ntri = new_ntri
td%nedge = new_nedge

end subroutine break_edge

!          --------------
subroutine pair_triangles(td,mate,NO_MATE)
!          --------------

!----------------------------------------------------
! This routine pairs up most of the triangles in the grid given by td.
!
! In this version, make one pass through the triangles pairing each
! unpaired triangle with the unpaired neighbor opposite the largest angle.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(triangle_data), intent(in) :: td
integer, intent(inout) :: mate(:)
integer, intent(in) :: NO_MATE
!----------------------------------------------------
! Local variables:

real(my_real) :: cos_angle(3,size(mate))
integer :: i, j, maxang
real(my_real) :: x1, x2, x3, y1, y2, y3, dx1, dx2, dy1, dy2, denom, mincos
!----------------------------------------------------
! Begin executable code

! compute the cosines of the angles of all the triangles

do i=1,td%ntri
   do j=1,3
      x1 = td%vert_coord(td%tri_vert(j,i))%x
      y1 = td%vert_coord(td%tri_vert(j,i))%y
      x2 = td%vert_coord(td%tri_vert(1+mod(j+1,3),i))%x
      y2 = td%vert_coord(td%tri_vert(1+mod(j+1,3),i))%y
      x3 = td%vert_coord(td%tri_vert(1+mod(j,3),i))%x
      y3 = td%vert_coord(td%tri_vert(1+mod(j,3),i))%y
      dx1 = x2-x1
      dy1 = y2-y1
      dx2 = x3-x1
      dy2 = y3-y1
      denom = sqrt((dx1*dx1+dy1*dy1)*(dx2*dx2+dy2*dy2))
      cos_angle(j,i) = (dx1*dx2+dy1*dy2)/denom
   end do
end do

! for each triangle not already assigned a mate

do i=1,td%ntri
   if (mate(i) /= NO_MATE) cycle


! find the largest angle (smallest cosine) opposite an unassigned or boundary
! neighbor

   maxang = -1
   mincos = 2.0_my_real
   do j=1,3
      if (td%tri_neigh(j,i) == -1) then
         if (cos_angle(j,i) < mincos .and. &
             abs(cos_angle(j,i)-mincos) > 100*epsilon(0.0_my_real)) then
            mincos = cos_angle(j,i)
            maxang = j
         endif
      elseif (mate(td%tri_neigh(j,i)) == NO_MATE) then
         if (cos_angle(j,i) < mincos .and. &
             abs(cos_angle(j,i)-mincos) > 100*epsilon(0.0_my_real)) then
            mincos = cos_angle(j,i)
            maxang = j
         endif
      endif
   end do

! set mate to indicate which of the three angles the mate is opposite

   if (maxang /= -1) then
      mate(i) = maxang

! if the mate is not the boundary, then find the angle in mate opposite
! triangle i and set its mate

      if (td%tri_neigh(maxang,i) /= -1) then
         do j=1,3
            if (td%tri_neigh(j,td%tri_neigh(maxang,i)) == i) then
               mate(td%tri_neigh(maxang,i)) = j
            endif
         end do
      endif

   endif

end do

end subroutine pair_triangles

!          ------------
subroutine refine_start(grid,td,mate,elemedge_bmark,NO_MATE)
!          ------------

!----------------------------------------------------
! This routine refines the starting grid in mate and td
! by bisecting those that have mates and trisecting those that don't.  It
! creates the following components of grid:
! vertex%coord, element%gid, element%mate, element%vertex, nelem, nvert,
! initial_neighbor
! It also sets vertex%assoc_elem to be the master vertex for vertices that
! are PERIODIC_SLAVE (negative bmark) or the endpoint of a PERIODIC_SLAVE
! side, and NO_MATE for all other vertices. This is just a convenient component
! to use for temporary storage.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(triangle_data), intent(in) :: td
integer, intent(inout) :: mate(:)
integer, intent(in) :: NO_MATE
integer, intent(out) :: elemedge_bmark(:,:)
!----------------------------------------------------
! Local variables:

integer :: i, j, k, next_vert, next_elem, &
           next_gid, tri, neigh, my_neigh, peak, old1, old2, peakm, old1m, &
           old2m, my_mate, stat, edge, bound_periodic_nvert, refedge, refedgem
integer :: refined_to(3,size(mate)), vert_edges(32,td%nvert)
!----------------------------------------------------
! Begin executable code

! determine the number of vertices and elements in the initial grid

bound_periodic_nvert = 0
grid%nvert = td%nvert
grid%nelem = 0
do i=1,td%ntri
   if (mate(i) == NO_MATE) then
      grid%nvert = grid%nvert + 1
      grid%nelem = grid%nelem + 3
   elseif (td%tri_neigh(mate(i),i) == -1) then ! boundary mate
      grid%nvert = grid%nvert + 1
      grid%nelem = grid%nelem + 2
   elseif (td%tri_neigh(mate(i),i) > i) then ! count pair when encountering first
      grid%nvert = grid%nvert + 1
      grid%nelem = grid%nelem + 4
   endif
   do j=1,EDGES_PER_ELEMENT
      if (td%vert_bmark(td%tri_vert(j,i)) < 0) then
         bound_periodic_nvert = bound_periodic_nvert + 1
      endif
   end do
end do

! make sure initial grid was allocated big enough and allocate initial_neighbor

if (size(grid%vertex) < grid%nvert + bound_periodic_nvert) then
   deallocate(grid%vertex)
   allocate(grid%vertex(2*grid%nvert),stat=stat)
   if (stat /= 0) then
      call fatal("allocation failed for grid%vertex in refine_start", &
                 intlist=(/grid%nvert/))
      stop
   endif
endif
if (size(grid%element) < grid%nelem) then
   deallocate(grid%element)
   allocate(grid%element(2*grid%nelem),stat=stat)
   if (stat /= 0) then
      call fatal("allocation failed for grid%element in refine_start", &
                 intlist=(/grid%nelem/))
      stop
   endif
endif
allocate(grid%initial_neighbor(EDGES_PER_ELEMENT,grid%nelem),stat=stat)
if (stat /= 0) then
   call fatal("allocation failed for grid%initial_neighbor in refine_start", &
              intlist=(/grid%nelem/))
   stop
endif

next_gid = grid%nelem

! copy vertex coordinates, boundary markers, and masters for PERIODIC_SLAVEs
! to grid data structure.  assoc_elem is used as temporary storage for masters.

do i=1,td%nvert
   grid%vertex(i)%coord = td%vert_coord(i)
   grid%vertex(i)%bmark = td%vert_bmark(i)
   grid%vertex(i)%bparam = td%vert_bparam(i)
   grid%vertex(i)%assoc_elem = td%vert_master(i)
end do
next_vert = td%nvert+1

! invert the mapping between edges and vertices

vert_edges = -1
do i=1,td%nedge
   do k=1,2
      do j=1,size(vert_edges,dim=1)
         if (vert_edges(j,td%edge_vert(k,i)) == -1) exit
      end do
      if (j > size(vert_edges,dim=1)) then
         call fatal("need bigger dimension for vert_edges")
         stop
      endif
      vert_edges(j,td%edge_vert(k,i)) = i
   end do
end do

! pass through elements, refining them and creating remainder of grid data

refined_to = -1
next_elem = 1
do tri=1,td%ntri

! see if it is already refined

   if (refined_to(1,tri) /= -1) cycle

! no mate, trisect it

   if (mate(tri) == NO_MATE) then

      grid%vertex(next_vert)%coord%x = sum(td%vert_coord(td%tri_vert(:,tri))%x)/3
      grid%vertex(next_vert)%coord%y = sum(td%vert_coord(td%tri_vert(:,tri))%y)/3
      grid%vertex(next_vert)%bmark = 0
      grid%vertex(next_vert)%assoc_elem = NO_MATE

      grid%element(next_elem  )%vertex(1) = td%tri_vert(2,tri)
      grid%element(next_elem  )%vertex(2) = td%tri_vert(3,tri)
      grid%element(next_elem  )%vertex(3) = next_vert

      grid%element(next_elem+1)%vertex(1) = td%tri_vert(3,tri)
      grid%element(next_elem+1)%vertex(2) = td%tri_vert(1,tri)
      grid%element(next_elem+1)%vertex(3) = next_vert

      grid%element(next_elem+2)%vertex(1) = td%tri_vert(1,tri)
      grid%element(next_elem+2)%vertex(2) = td%tri_vert(2,tri)
      grid%element(next_elem+2)%vertex(3) = next_vert

      grid%element(next_elem  )%gid = next_gid
      grid%element(next_elem+1)%gid = next_gid+1
      grid%element(next_elem+2)%gid = next_gid+2

      refined_to(1,tri) = next_elem
      refined_to(2,tri) = next_elem+1
      refined_to(3,tri) = next_elem+2

      do j=0,2
         elemedge_bmark(:,next_elem+j) = 0
         do i=1,size(vert_edges,dim=1)
            edge = vert_edges(i,grid%element(next_elem+j)%vertex(1))
            if (edge == -1) exit
            if (td%edge_vert(1,edge) == grid%element(next_elem+j)%vertex(2) .or. &
                td%edge_vert(2,edge) == grid%element(next_elem+j)%vertex(2)) exit
         end do
         if (edge == -1 .or. i > size(vert_edges,dim=1)) then
            call fatal("failed to find edge")
            stop
         endif
         elemedge_bmark(3,next_elem+j) = td%edge_bmark(edge)
      end do

      grid%initial_neighbor(1,next_elem  ) = next_elem+1
      grid%initial_neighbor(2,next_elem  ) = next_elem+2
      grid%initial_neighbor(3,next_elem  ) = -1

      grid%initial_neighbor(1,next_elem+1) = next_elem+2
      grid%initial_neighbor(2,next_elem+1) = next_elem
      grid%initial_neighbor(3,next_elem+1) = -1

      grid%initial_neighbor(1,next_elem+2) = next_elem
      grid%initial_neighbor(2,next_elem+2) = next_elem+1
      grid%initial_neighbor(3,next_elem+2) = -1

      do j=1,3
         neigh = td%tri_neigh(j,tri)

         if (neigh == -1) then
            grid%initial_neighbor(3,next_elem+j-1) = BOUNDARY
            grid%element(next_elem+j-1)%mate = BOUNDARY
            cycle
         endif

         if (refined_to(1,neigh) == -1) cycle ! wait until neighbor is refined

         my_neigh = find_neigh(refined_to(:,neigh), &
                               grid%element(next_elem+j-1)%vertex, &
                               grid%element,grid%vertex)
         grid%initial_neighbor(3,next_elem+j-1) = my_neigh
         grid%element(next_elem+j-1)%mate = grid%element(my_neigh)%gid
         grid%initial_neighbor(3,my_neigh) = next_elem+j-1
         grid%element(my_neigh)%mate = grid%element(next_elem+j-1)%gid
      end do

      next_vert = next_vert + 1
      next_elem = next_elem + 3
      next_gid = next_gid + 3

! mate is boundary, bisect tri

   elseif (td%tri_neigh(mate(tri),tri) == -1) then

      peak = mate(tri)
      old1 = 1
      old2 = 2
      if (peak == 1) old1 = 3
      if (peak == 2) old2 = 3

      do i=1,size(vert_edges,dim=1)
         edge = vert_edges(i,td%tri_vert(old1,tri))
         if (edge == -1) exit
         if (td%edge_vert(1,edge) == td%tri_vert(old2,tri) .or. &
             td%edge_vert(2,edge) == td%tri_vert(old2,tri)) exit
      end do
      if (edge == -1 .or. i > size(vert_edges,dim=1)) then
         call fatal("failed to find edge")
         stop
      endif

      call point_on_edge(grid, &
                         td%vert_coord(td%tri_vert(old1,tri))%x, &
                         td%vert_coord(td%tri_vert(old1,tri))%y, &
                         td%vert_bmark(td%tri_vert(old1,tri)), &
                         td%vert_bparam(td%tri_vert(old1,tri)), &
                         td%vert_coord(td%tri_vert(old2,tri))%x, &
                         td%vert_coord(td%tri_vert(old2,tri))%y, &
                         td%vert_bmark(td%tri_vert(old2,tri)), &
                         td%vert_bparam(td%tri_vert(old2,tri)), &
                         td%vert_coord(td%tri_vert(peak,tri))%x, &
                         td%vert_coord(td%tri_vert(peak,tri))%y, &
                         td%edge_bmark(edge), &
                         0.5_my_real, &
                         grid%vertex(next_vert)%coord%x, &
                         grid%vertex(next_vert)%coord%y, &
                         grid%vertex(next_vert)%bparam)
      grid%vertex(next_vert)%bmark = td%edge_bmark(edge)
      grid%vertex(next_vert)%assoc_elem = NO_MATE

      grid%element(next_elem  )%vertex(1) = td%tri_vert(peak,tri)
      grid%element(next_elem  )%vertex(2) = td%tri_vert(old1,tri)
      grid%element(next_elem  )%vertex(3) = next_vert

      grid%element(next_elem+1)%vertex(1) = td%tri_vert(peak,tri)
      grid%element(next_elem+1)%vertex(2) = td%tri_vert(old2,tri)
      grid%element(next_elem+1)%vertex(3) = next_vert

      grid%element(next_elem  )%gid = next_gid
      grid%element(next_elem+1)%gid = next_gid+1

      refined_to(1,tri) = next_elem
      refined_to(2,tri) = next_elem+1
      refined_to(3,tri) = -1

      do j=0,1
         elemedge_bmark(:,next_elem+j) = 0
         elemedge_bmark(1,next_elem+j) = grid%vertex(next_vert)%bmark
         do i=1,size(vert_edges,dim=1)
            edge = vert_edges(i,grid%element(next_elem+j)%vertex(1))
            if (edge == -1) exit
            if (td%edge_vert(1,edge) == grid%element(next_elem+j)%vertex(2) .or. &
                td%edge_vert(2,edge) == grid%element(next_elem+j)%vertex(2)) exit
         end do
         if (edge == -1 .or. i > size(vert_edges,dim=1)) then
            call fatal("failed to find edge")
            stop
         endif
         elemedge_bmark(3,next_elem+j) = td%edge_bmark(edge)
      end do

      grid%initial_neighbor(1,next_elem  ) = BOUNDARY
      grid%initial_neighbor(2,next_elem  ) = next_elem+1
      grid%initial_neighbor(3,next_elem  ) = -1

      grid%initial_neighbor(1,next_elem+1) = BOUNDARY
      grid%initial_neighbor(2,next_elem+1) = next_elem
      grid%initial_neighbor(3,next_elem+1) = -1

      do j=1,2
         if (j == 1) then
            neigh = td%tri_neigh(old2,tri)
         else
            neigh = td%tri_neigh(old1,tri)
         endif

         if (neigh == -1) then
            grid%initial_neighbor(3,next_elem+j-1) = BOUNDARY
            grid%element(next_elem+j-1)%mate = BOUNDARY
            cycle
         endif

         if (refined_to(1,neigh) == -1) cycle ! wait until neighbor is refined

         my_neigh = find_neigh(refined_to(:,neigh), &
                               grid%element(next_elem+j-1)%vertex, &
                               grid%element,grid%vertex)
         grid%initial_neighbor(3,next_elem+j-1) = my_neigh
         grid%element(next_elem+j-1)%mate = grid%element(my_neigh)%gid
         grid%initial_neighbor(3,my_neigh) = next_elem+j-1
         grid%element(my_neigh)%mate = grid%element(next_elem+j-1)%gid
      end do

      next_vert = next_vert + 1
      next_elem = next_elem + 2
      next_gid = next_gid + 2

! mate is another triangle, bisect pair of triangles

   else

      my_mate = td%tri_neigh(mate(tri),tri)
      if (td%tri_neigh(mate(my_mate),my_mate) /= tri) then
         call fatal("mates are not symmetric in refine_start")
         stop
      endif

      peak = mate(tri)
      old1 = 1
      old2 = 2
      if (peak == 1) old1 = 3
      if (peak == 2) old2 = 3

      do i=1,size(vert_edges,dim=1)
         refedge = vert_edges(i,td%tri_vert(old1,tri))
         if (refedge == -1) exit
         if (td%edge_vert(1,refedge) == td%tri_vert(old2,tri) .or. &
             td%edge_vert(2,refedge) == td%tri_vert(old2,tri)) exit
      end do
      if (refedge == -1 .or. i > size(vert_edges,dim=1)) then
         call fatal("failed to find edge")
         stop
      endif

      peakm = mate(my_mate)
      old1m = 1 ! temporary for finding edge; will set after edge is determined
      old2m = 2
      if (peakm == 1) old1m = 3
      if (peakm == 2) old2m = 3

      do i=1,size(vert_edges,dim=1)
         refedgem = vert_edges(i,td%tri_vert(old1m,my_mate))
         if (refedgem == -1) exit
         if (td%edge_vert(1,refedgem) == td%tri_vert(old2m,my_mate) .or. &
             td%edge_vert(2,refedgem) == td%tri_vert(old2m,my_mate)) exit
      end do
      if (refedgem == -1 .or. i > size(vert_edges,dim=1)) then
         call fatal("failed to find edge")
         stop
      endif

      old1m = -1
      old2m = -1
      do j=1,3
         if (td%edge_bmark(refedge) < 0 .and. &
             td%vert_master(td%tri_vert(old1,tri)) /= NO_MATE) then
            if (td%tri_vert(j,my_mate) == td%vert_master(td%tri_vert(old1,tri))) old1m = j
         elseif (td%edge_bmark(refedgem) < 0 .and. &
                 td%vert_master(td%tri_vert(j,my_mate)) /= NO_MATE) then
            if (td%vert_master(td%tri_vert(j,my_mate)) == td%tri_vert(old1,tri)) old1m = j
         else
            if (td%tri_vert(j,my_mate) == td%tri_vert(old1,tri)) old1m = j
         endif
         if (td%edge_bmark(refedge) < 0 .and. &
             td%vert_master(td%tri_vert(old2,tri)) /= NO_MATE) then
            if (td%tri_vert(j,my_mate) == td%vert_master(td%tri_vert(old2,tri))) old2m = j
         elseif (td%edge_bmark(refedgem) < 0 .and. &
                 td%vert_master(td%tri_vert(j,my_mate)) /= NO_MATE) then
            if (td%vert_master(td%tri_vert(j,my_mate)) == td%tri_vert(old2,tri)) old2m = j
         else
            if (td%tri_vert(j,my_mate) == td%tri_vert(old2,tri)) old2m = j
         endif
      end do
      if (old1m == -1 .or. old2m == -1) then
         call fatal("didn't find matching vertices in mate in refine_start")
         stop
      endif

      call point_on_edge(grid, &
                         td%vert_coord(td%tri_vert(old1,tri))%x, &
                         td%vert_coord(td%tri_vert(old1,tri))%y, &
                         td%vert_bmark(td%tri_vert(old1,tri)), &
                         td%vert_bparam(td%tri_vert(old1,tri)), &
                         td%vert_coord(td%tri_vert(old2,tri))%x, &
                         td%vert_coord(td%tri_vert(old2,tri))%y, &
                         td%vert_bmark(td%tri_vert(old2,tri)), &
                         td%vert_bparam(td%tri_vert(old2,tri)), &
                         td%vert_coord(td%tri_vert(peak,tri))%x, &
                         td%vert_coord(td%tri_vert(peak,tri))%y, &
                         td%edge_bmark(refedge), &
                         0.5_my_real, &
                         grid%vertex(next_vert)%coord%x, &
                         grid%vertex(next_vert)%coord%y, &
                         grid%vertex(next_vert)%bparam)
      grid%vertex(next_vert)%bmark = td%edge_bmark(refedge)
      grid%vertex(next_vert)%assoc_elem = NO_MATE

      grid%element(next_elem  )%vertex(1) = td%tri_vert(peak,tri)
      grid%element(next_elem  )%vertex(2) = td%tri_vert(old1,tri)
      grid%element(next_elem  )%vertex(3) = next_vert

      grid%element(next_elem+1)%vertex(1) = td%tri_vert(peak,tri)
      grid%element(next_elem+1)%vertex(2) = td%tri_vert(old2,tri)
      grid%element(next_elem+1)%vertex(3) = next_vert

      if (td%edge_bmark(refedge) < 0 .or. td%edge_bmark(refedgem) < 0) then
         next_vert = next_vert + 1
         grid%nvert = grid%nvert + 1
         call point_on_edge(grid, &
                            td%vert_coord(td%tri_vert(old1m,my_mate))%x, &
                            td%vert_coord(td%tri_vert(old1m,my_mate))%y, &
                            td%vert_bmark(td%tri_vert(old1m,my_mate)), &
                            td%vert_bparam(td%tri_vert(old1m,my_mate)), &
                            td%vert_coord(td%tri_vert(old2m,my_mate))%x, &
                            td%vert_coord(td%tri_vert(old2m,my_mate))%y, &
                            td%vert_bmark(td%tri_vert(old2m,my_mate)), &
                            td%vert_bparam(td%tri_vert(old2m,my_mate)), &
                            td%vert_coord(td%tri_vert(peakm,my_mate))%x, &
                            td%vert_coord(td%tri_vert(peakm,my_mate))%y, &
                            td%edge_bmark(refedgem), &
                            0.5_my_real, &
                            grid%vertex(next_vert)%coord%x, &
                            grid%vertex(next_vert)%coord%y, &
                            grid%vertex(next_vert)%bparam)
         grid%vertex(next_vert)%bmark = td%edge_bmark(refedgem)
         grid%vertex(next_vert)%assoc_elem = NO_MATE
         if (td%edge_bmark(refedge) < 0) then
            grid%vertex(next_vert-1)%assoc_elem = next_vert
         elseif (td%edge_bmark(refedgem) < 0) then
            grid%vertex(next_vert)%assoc_elem = next_vert-1
         endif
      endif

      grid%element(next_elem+2)%vertex(1) = td%tri_vert(peakm,my_mate)
      grid%element(next_elem+2)%vertex(2) = td%tri_vert(old1m,my_mate)
      grid%element(next_elem+2)%vertex(3) = next_vert

      grid%element(next_elem+3)%vertex(1) = td%tri_vert(peakm,my_mate)
      grid%element(next_elem+3)%vertex(2) = td%tri_vert(old2m,my_mate)
      grid%element(next_elem+3)%vertex(3) = next_vert

      grid%element(next_elem  )%gid = next_gid
      grid%element(next_elem+1)%gid = next_gid+1
      grid%element(next_elem+2)%gid = next_gid+2
      grid%element(next_elem+3)%gid = next_gid+3

      refined_to(1,tri) = next_elem
      refined_to(2,tri) = next_elem+1
      refined_to(3,tri) = -1
      refined_to(1,my_mate) = next_elem+2
      refined_to(2,my_mate) = next_elem+3
      refined_to(3,my_mate) = -1

      do j=0,3
         elemedge_bmark(2,next_elem+j) = 0
         do i=1,size(vert_edges,dim=1)
            edge = vert_edges(i,grid%element(next_elem+j)%vertex(1))
            if (edge == -1) exit
            if (td%edge_vert(1,edge) == grid%element(next_elem+j)%vertex(2) .or. &
                td%edge_vert(2,edge) == grid%element(next_elem+j)%vertex(2)) exit
         end do
         if (edge == -1 .or. i > size(vert_edges,dim=1)) then
            call fatal("failed to find edge")
            stop
         endif
         elemedge_bmark(3,next_elem+j) = td%edge_bmark(edge)
      end do
      elemedge_bmark(1,next_elem  ) = td%edge_bmark(refedge)
      elemedge_bmark(1,next_elem+1) = td%edge_bmark(refedge)
      elemedge_bmark(1,next_elem+2) = td%edge_bmark(refedgem)
      elemedge_bmark(1,next_elem+3) = td%edge_bmark(refedgem)

      grid%initial_neighbor(1,next_elem  ) = next_elem+2
      grid%initial_neighbor(2,next_elem  ) = next_elem+1
      grid%initial_neighbor(3,next_elem  ) = -1

      grid%initial_neighbor(1,next_elem+1) = next_elem+3
      grid%initial_neighbor(2,next_elem+1) = next_elem
      grid%initial_neighbor(3,next_elem+1) = -1

      grid%initial_neighbor(1,next_elem+2) = next_elem
      grid%initial_neighbor(2,next_elem+2) = next_elem+3
      grid%initial_neighbor(3,next_elem+2) = -1

      grid%initial_neighbor(1,next_elem+3) = next_elem+1
      grid%initial_neighbor(2,next_elem+3) = next_elem+2
      grid%initial_neighbor(3,next_elem+3) = -1

      do j=1,4
         select case(j)
         case(1)
            neigh = td%tri_neigh(old2,tri)
         case(2)
            neigh = td%tri_neigh(old1,tri)
         case(3)
            neigh = td%tri_neigh(old2m,my_mate)
         case(4)
            neigh = td%tri_neigh(old1m,my_mate)
         end select

         if (neigh == -1) then
            grid%initial_neighbor(3,next_elem+j-1) = BOUNDARY
            grid%element(next_elem+j-1)%mate = BOUNDARY
            cycle
         endif

         if (refined_to(1,neigh) == -1) cycle ! wait until neighbor is refined

         my_neigh = find_neigh(refined_to(:,neigh), &
                               grid%element(next_elem+j-1)%vertex, &
                               grid%element,grid%vertex)
         grid%initial_neighbor(3,next_elem+j-1) = my_neigh
         grid%element(next_elem+j-1)%mate = grid%element(my_neigh)%gid
         grid%initial_neighbor(3,my_neigh) = next_elem+j-1
         grid%element(my_neigh)%mate = grid%element(next_elem+j-1)%gid
      end do

      next_vert = next_vert + 1
      next_elem = next_elem + 4
      next_gid = next_gid + 4

   endif ! no mate, boundary or paired

end do

end subroutine refine_start

!        ----------
function find_neigh(search_in,search_for,elements,vertices)
!        ----------

!----------------------------------------------------
! This routine looks through the elements indexed by search_in to find one with
! the first two vertices in search_for
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: search_in(3),search_for(:)
type(element_t), intent(in) :: elements(:)
type(vertex_t), intent(in) :: vertices(:)
integer :: find_neigh
!----------------------------------------------------
! Local variables:

integer :: j
!----------------------------------------------------
! Begin executable code

do j=1,3
   if (search_in(j) == -1) cycle
   if ((search_for(1) == elements(search_in(j))%vertex(1) .or. &
        search_for(1) == elements(search_in(j))%vertex(2) .or. &
        search_for(1) == elements(search_in(j))%vertex(3)) .and. &
       (search_for(2) == elements(search_in(j))%vertex(1) .or. &
        search_for(2) == elements(search_in(j))%vertex(2) .or. &
        search_for(2) == elements(search_in(j))%vertex(3))) then
      find_neigh = search_in(j)
      return
   endif
   if ((search_for(1) == vertices(elements(search_in(j))%vertex(1))%assoc_elem .or. &
        search_for(1) == vertices(elements(search_in(j))%vertex(2))%assoc_elem .or. &
        search_for(1) == vertices(elements(search_in(j))%vertex(3))%assoc_elem) .and. &
       (search_for(2) == vertices(elements(search_in(j))%vertex(1))%assoc_elem .or. &
        search_for(2) == vertices(elements(search_in(j))%vertex(2))%assoc_elem .or. &
        search_for(2) == vertices(elements(search_in(j))%vertex(3))%assoc_elem)) then
      find_neigh = search_in(j)
      return
   endif
   if ((vertices(search_for(1))%assoc_elem == elements(search_in(j))%vertex(1) .or. &
        vertices(search_for(1))%assoc_elem == elements(search_in(j))%vertex(2) .or. &
        vertices(search_for(1))%assoc_elem == elements(search_in(j))%vertex(3)) .and. &
       (vertices(search_for(2))%assoc_elem == elements(search_in(j))%vertex(1) .or. &
        vertices(search_for(2))%assoc_elem == elements(search_in(j))%vertex(2) .or. &
        vertices(search_for(2))%assoc_elem == elements(search_in(j))%vertex(3))) then
      find_neigh = search_in(j)
      return
   endif
end do

call fatal("didn't find neighbor that shares two vertices in find_neigh")
stop

end function find_neigh

!          -----------
subroutine smooth_grid(grid,procs)
!          -----------

!----------------------------------------------------
! This routine moves non-boundary vertices to smooth the triangle shapes
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(proc_info), intent(in) :: procs
!----------------------------------------------------
! Local variables:

integer, allocatable :: neigh_vert(:,:)
integer :: astat, elem, i, j, k, vert, vert2
real(my_real) :: x,y

!----------------------------------------------------
! Begin executable code

if (my_proc(procs) == MASTER) return

! For each vertex, construct a list of surrounding vertices
! Surely there won't be 32 neighbors

allocate(neigh_vert(32,grid%nvert),stat=astat)
if (astat /= 0) then
   call fatal("allocation failed for neigh_vert in smooth_grid")
   stop
endif
neigh_vert = -1

elem = grid%head_level_elem(1)
do while (elem /= END_OF_LIST)
   do vert = 1,3
      do vert2 = 1,3
         if (vert == vert2) cycle
         do i=1,32
            if (neigh_vert(i,grid%element(elem)%vertex(vert)) == -1) exit
            if (neigh_vert(i,grid%element(elem)%vertex(vert)) == &
                grid%element(elem)%vertex(vert2)) exit
         end do
         if (i==32) then
            call fatal("more than 32 neighbors in initial grid")
            stop
         endif
         if (neigh_vert(i,grid%element(elem)%vertex(vert)) == -1) then
            neigh_vert(i,grid%element(elem)%vertex(vert)) = &
               grid%element(elem)%vertex(vert2)
         endif
      end do
   end do
   elem = grid%element(elem)%next
end do

! Simple Laplacian smoothing.  Each interior vertex gets moved to the
! geometric center of its neighbors.

do k=1,10 ! do ten iterations
   do i=1,grid%nvert
      if (grid%vertex_type(i,1) /= INTERIOR) cycle
      x = 0.0_my_real
      y = 0.0_my_real
      do j=1,32
         if (neigh_vert(j,i) == -1) exit
         x = x + grid%vertex(neigh_vert(j,i))%coord%x
         y = y + grid%vertex(neigh_vert(j,i))%coord%y
      end do
      grid%vertex(i)%coord = point(x/(j-1),y/(j-1))
   end do
end do

end subroutine smooth_grid

!          ---------
subroutine init_path(grid)
!          ---------

!----------------------------------------------------
! This routine creates a path through the elements of the initial grid,
! along with in- and out- vertices, using a Sierpinski space filling curve.
! It is possible for the path to have discontinuities.
!
! This routine sets grid%head_level_elem(1) and, for each element of the
! initial grid, element%next, element%previous, element%in and element%out.
! It assumes the grid%nelem elements are in the first grid%nelem entries
! of grid%element, 
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
!----------------------------------------------------
! Local variables:

real(my_real) :: xcent, ycent, sfccoord(grid%nelem)
real(my_real) :: xmin, xmax, ymin, ymax
integer :: iperm(grid%nelem)
integer :: i, elem, prev, prevprev, jerr
!----------------------------------------------------
! Begin executable code

! For each element, determine the SFC mapping from its center to 1D

xmin = grid%boundbox_min%x
xmax = grid%boundbox_max%x
ymin = grid%boundbox_min%y
ymax = grid%boundbox_max%y

do elem=1,grid%nelem
   xcent = 0.0_my_real
   ycent = 0.0_my_real
   do i=1,VERTICES_PER_ELEMENT
      xcent = xcent + grid%vertex(grid%element(elem)%vertex(i))%coord%x
      ycent = ycent + grid%vertex(grid%element(elem)%vertex(i))%coord%y
   end do
   xcent = xcent/3
   ycent = ycent/3
   xcent = (xcent-xmin)/(xmax-xmin)
   ycent = (ycent-ymin)/(ymax-ymin)
   sfccoord(elem) = invsierpinski2d(xcent,ycent)
end do

! Sort the SFC coordinates to get the order of the elements

call sort(sfccoord,grid%nelem,iperm,1,jerr)

! Pass through the elements in order setting previous/next and looking
! for in/out vertices

elem = iperm(1)
grid%head_level_elem(1) = elem
grid%element(elem)%previous = END_OF_LIST
grid%element(elem)%in = grid%element(elem)%vertex(1)
prev = elem

do i=2,grid%nelem
   elem = iperm(i)
   grid%element(prev)%next = elem
   grid%element(elem)%previous = prev
   if (i == 2) then
      prevprev = -1
   else
      prevprev = iperm(i-2)
   endif
   call find_inout(elem,prev,prevprev,grid)
   prev = elem
end do

grid%element(elem)%next = END_OF_LIST
if (grid%element(elem)%in == grid%element(elem)%vertex(1)) then
   grid%element(elem)%out = grid%element(elem)%vertex(2)
else
   grid%element(elem)%out = grid%element(elem)%vertex(1)
endif

end subroutine init_path

!          ----------
subroutine find_inout(elem,prev,prevprev,grid)
!          ----------

!----------------------------------------------------
! This routine finds and sets an in-vertex for eleme and out-vertex for prev.
! Usually they will be the same vertex, but if necessary they will be different
! which will cause the Hamiltonian path to be disconnected.  If necessary and
! possible, this routine will change the in-vertex of prev and out-vertex
! of prevprev to make in and out be the same.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: elem, prev, prevprev
type(grid_type), intent(inout) :: grid
!----------------------------------------------------
! Local variables:

integer :: i, j
logical :: looking, shared
!----------------------------------------------------
! Begin executable code

! First look for a shared vertex of elem and prev that is not the
! in-vertex of prev

looking = .true.
do i=1,VERTICES_PER_ELEMENT
   do j=1,VERTICES_PER_ELEMENT
      if (grid%element(elem)%vertex(i) == grid%element(prev)%vertex(j) .and. &
          grid%element(prev)%vertex(j) /= grid%element(prev)%in) then
         grid%element(prev)%out = grid%element(prev)%vertex(j)
         grid%element(elem)%in  = grid%element(elem)%vertex(i)
         looking = .false.
         exit
      endif
   end do
   if (.not. looking) exit
end do

! If that failed then the in-vertex of prev must be shared, or else there
! are no shared vertices between elem and prev

if (looking) then
   shared = .false.
   do i=1,VERTICES_PER_ELEMENT
      if (grid%element(elem)%vertex(i) == grid%element(prev)%in) then
         shared = .true.
         exit
      endif
   end do

! If the in-vertex of prev is not shared, then prev and elem are not adjacent,
! so pick any vertices for the in and out

   if (.not. shared) then
      grid%element(elem)%in = grid%element(elem)%vertex(1)
      if (grid%element(prev)%in == grid%element(prev)%vertex(1)) then
         grid%element(prev)%out = grid%element(prev)%vertex(2)
      else
         grid%element(prev)%out = grid%element(prev)%vertex(1)
      endif
      looking = .false.
   endif
endif

! If the in-vertex of prev is shared, then it must be the only shared vertex
! between elem and prev (otherwise the first approach would have worked).
! Try to change the in-vertex of prev and use it as the out-vertex of prev
! and in-vertex of elem.

! First, if prev is the beginning of the path then set in-vertex to be
! any other vertex

if (looking) then
   if (prevprev == -1) then
      grid%element(elem)%in = grid%element(prev)%in
      grid%element(prev)%out = grid%element(prev)%in
      if (grid%element(prev)%in == grid%element(prev)%vertex(1)) then
         grid%element(prev)%in = grid%element(prev)%vertex(2)
      else
         grid%element(prev)%in = grid%element(prev)%vertex(1)
      endif
      looking = .false.
   endif
endif

! Second, if prev and prevprev are not adjacent or contain a broken link,
! i.e., prevprev%out is not prev%in, then change the in-vertex of prev to
! any other vertex

if (looking) then
   if (grid%element(prev)%in /= grid%element(prevprev)%out) then
      grid%element(elem)%in = grid%element(prev)%in
      grid%element(prev)%out = grid%element(prev)%in
      if (grid%element(prev)%in == grid%element(prev)%vertex(1)) then
         grid%element(prev)%in = grid%element(prev)%vertex(2)
      else
         grid%element(prev)%in = grid%element(prev)%vertex(1)
      endif
      looking = .false.
   endif
endif

! Third and final, look for a shared vertex between prev and prevprev that
! is not the in-vertex of either prev or prevprev and use that as the new
! out-vertex of prevprev and in-vertex of prev.

if (looking) then
   do i=1,VERTICES_PER_ELEMENT
      if (grid%element(prev)%in /= grid%element(prev)%vertex(i)) then
         do j=1,VERTICES_PER_ELEMENT
            if (grid%element(prevprev)%vertex(j)==grid%element(prev)%vertex(i) &
         .and. grid%element(prevprev)%vertex(j)/=grid%element(prevprev)%in) then
               grid%element(elem)%in = grid%element(prev)%in
               grid%element(prev)%out = grid%element(prev)%in
               grid%element(prev)%in = grid%element(prev)%vertex(i)
               grid%element(prevprev)%out = grid%element(prev)%vertex(i)
               looking = .false.
               exit
            endif
         end do
      endif
      if (.not. looking) exit
   end do
endif

! If that failed, give up and put in a broken link by using the in-vertex
! of prev (the only shared vertex between elem and prev) as the in-vertex
! of elem and any other vertex of prev as the out-vertex of prev.

if (looking) then
   grid%element(elem)%in = grid%element(prev)%in
   if (grid%element(prev)%in == grid%element(prev)%vertex(1)) then
      grid%element(prev)%out = grid%element(prev)%vertex(2)
   else
      grid%element(prev)%out = grid%element(prev)%vertex(1)
   endif
endif

end subroutine find_inout

!        ---------------
function invsierpinski2d(x,y)
!        ---------------

!----------------------------------------------------
! This routine returns the Sierpinski key in [0,1] for the coordinates
! (x,y) in [0,1]X[0,1]
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x,y
real(my_real) :: invsierpinski2d
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

! sanity check for input arguments

if (x < 0 .or. x > 1 .or. y < 0 .or. y > 1) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("Input for Sierpinski out of range.")
   stop
endif

! Begin recursion that computes key

invsierpinski2d = sier2d(x, y, 0, 0, 0.5_my_real, 0.5_my_real, 0.0_my_real, &
                         0.0_my_real)

end function invsierpinski2d

!                  ------
recursive function sier2d(x, y, state, level, addf, addp, peakx, peaky) result(res)
!                  ------

!----------------------------------------------------
! This routine recursively computes the Sierpinski mapping from the unit square
! to the unit interval.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x, y, addf, addp, peakx, peaky
integer, intent(in) :: state, level
real(my_real) :: res
!----------------------------------------------------
! Local variables:

integer, parameter :: MAXLEV = 24
!----------------------------------------------------
! Begin executable code

if (level >= MAXLEV) then
   res = 0.0_my_real
   return
endif

select case(state)

case(0)
   if (y > x) then
      res = sier2d(x,y,1,level+1,addf/2,addp,0.0_my_real,1.0_my_real)
   else
      res = sier2d(x,y,2,level+1,addf/2,addp,1.0_my_real,0.0_my_real) + addf
   endif

case(1)
   if ( y-peaky < -(x-peakx) ) then
      res = sier2d(x,y,5,level+1,addf/2,addp,peakx+addp,peaky-addp)
   else
      res = sier2d(x,y,6,level+1,addf/2,addp,peakx+addp,peaky-addp) + addf
   endif

case(2)
   if ( y-peaky >= -(x-peakx) ) then
      res = sier2d(x,y,7,level+1,addf/2,addp,peakx-addp,peaky+addp)
   else
      res = sier2d(x,y,8,level+1,addf/2,addp,peakx-addp,peaky+addp) + addf
   endif

case(3)
   if ( y-peaky <= x-peakx ) then
      res = sier2d(x,y,8,level+1,addf/2,addp,peakx+addp,peaky+addp)
   else
      res = sier2d(x,y,5,level+1,addf/2,addp,peakx+addp,peaky+addp) + addf
   endif

case(4)
   if ( y-peaky > x-peakx ) then
      res = sier2d(x,y,6,level+1,addf/2,addp,peakx-addp,peaky-addp)
   else
      res = sier2d(x,y,7,level+1,addf/2,addp,peakx-addp,peaky-addp) + addf
   endif

case(5)
   if ( y < peaky ) then
      res = sier2d(x,y,1,level+1,addf/2,addp/2,peakx-addp,peaky)
   else
      res = sier2d(x,y,3,level+1,addf/2,addp/2,peakx-addp,peaky) + addf
   endif

case(6)
   if ( x < peakx ) then
      res = sier2d(x,y,4,level+1,addf/2,addp/2,peakx,peaky+addp)
   else
      res = sier2d(x,y,1,level+1,addf/2,addp/2,peakx,peaky+addp) + addf
   endif

case(7)
   if ( y >= peaky ) then
      res = sier2d(x,y,2,level+1,addf/2,addp/2,peakx+addp,peaky)
   else
      res = sier2d(x,y,4,level+1,addf/2,addp/2,peakx+addp,peaky) + addf
   endif

case(8)
   if ( x >= peakx ) then
      res = sier2d(x,y,3,level+1,addf/2,addp/2,peakx,peaky-addp)
   else
      res = sier2d(x,y,2,level+1,addf/2,addp/2,peakx,peaky-addp) + addf
   endif

end select

end function sier2d

!---------------------------------------------------------------------
!  GRID REFINEMENT ROUTINES
!---------------------------------------------------------------------

!                    ------
recursive subroutine refine(grid,procs,refine_control,solver_control, &
                            io_control,still_sequential,init_nvert,init_nelem, &
                            init_dof,loop,balance_what,predictive,no_time)
!                    ------

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

integer :: elem, errcode, astat, nproc, my_processor, target
real(my_real) :: global_max_errind
integer, pointer :: numhref(:), numpref(:)
! If len=1 is changed, also change it for temp_reftype in more_elements
character(len=1), pointer :: reftype(:)
type(errind_list) :: elist
logical :: return_to_elist, complete_elist, one_elist, target_met, timeit

!----------------------------------------------------
! Begin executable code

grid_changed = .true. ! says we need to resend graphics data

! The REFSOLN_EDGE and REFSOLN_ELEM hp-adaptive strategies have their own
! refine routine

if (refine_control%reftype == HP_ADAPTIVE .and. &
    (refine_control%hp_strategy == HP_REFSOLN_EDGE .or. &
     refine_control%hp_strategy == HP_REFSOLN_ELEM)) then
   call refine_refsoln(grid,procs,refine_control,solver_control, &
                       io_control,still_sequential,init_nvert,init_nelem, &
                       init_dof,loop,balance_what,predictive)
   return
endif

! NLP also has it's own routine

if (refine_control%reftype == HP_ADAPTIVE .and. &
    refine_control%hp_strategy == HP_NLP) then
   call refine_nlp(grid,procs,refine_control,still_sequential)
   return
endif

! MASTER doesn't participate

if (my_proc(procs) == MASTER) return

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
   if (timeit) call stop_watch((/prefine,trefine/))
   return
endif

call set_numref(grid,procs,still_sequential,refine_control,numhref,numpref)

! create the lists that group elements into bins based on the error indicator

allocate(elist%next_errind(size(grid%element)), &
         elist%prev_errind(size(grid%element)),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in refine",procs=procs)
   if (timeit) call stop_watch((/prefine,trefine/))
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
   if (timeit) call stop_watch((/prefine,trefine/))
   return
endif

call mark_reftype(grid,refine_control,global_max_errind,reftype,numhref, &
                  numpref,elist)

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
                            reftype,numhref,numpref,.false.)
   endif

! refine the element

   if (reftype(elem) == "h") then
      call bisect_triangle_pair(grid,elem,errcode,refine_control,elist,reftype,&
                                numhref,numpref,return_to_elist)
   elseif (reftype(elem) == "p") then
      call p_refine_elem(grid,elem,refine_control,elist,reftype,numhref, &
                         numpref,return_to_elist)
   else
      call remove_from_errind_list(elem,elist)
   endif

end do ! main refine loop

! free memory

deallocate(elist%next_errind,elist%prev_errind,reftype,numhref,numpref, &
           stat=astat)

! error indicators should be recomputed TEMP or should they?

grid%errind_up2date = .false.

! stop timing the refinement process

if (timeit) call stop_watch((/prefine,trefine/))

end subroutine refine

!        -------------------------
function compute_global_max_errind(grid,procs)
!        -------------------------

!----------------------------------------------------
! This routine computes the global maximum error indicator
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
type(proc_info), intent(in) :: procs
real(my_real) :: compute_global_max_errind
!----------------------------------------------------
! Local variables:

integer :: lev, elem
!----------------------------------------------------
! Begin executable code

compute_global_max_errind = 0.0_my_real
do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (grid%element(elem)%iown .and. grid%element(elem)%isleaf) then
         compute_global_max_errind = max(compute_global_max_errind, &
                    maxval(grid%element_errind(elem,:))/grid%element(elem)%work)
      endif
      elem = grid%element(elem)%next
   end do
end do

compute_global_max_errind = phaml_global_max(procs,compute_global_max_errind,440)

end function compute_global_max_errind

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

if (refine_control%reftype == H_UNIFORM .or. &
    refine_control%reftype == P_UNIFORM .or. &
    refine_control%refterm == ONE_REF .or. &
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
   total_num_big = phaml_global_sum(procs,my_num_big,432)
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

   call get_grid_info(grid,procs,still_sequential,445,total_nvert=total,&
                      no_master=.true.)
   if (total > refine_control%max_vert) total = refine_control%max_vert
   target = total*my_fraction

case (KEEP_NELEM, KEEP_NELEM_SMOOTH)

   call get_grid_info(grid,procs,still_sequential,445, &
                      total_nelem_leaf=total,no_master=.true.)
   if (total > refine_control%max_elem) total = refine_control%max_elem
   target = total*my_fraction

case (KEEP_NEQ, KEEP_NEQ_SMOOTH)

   call get_grid_info(grid,procs,still_sequential,445, &
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

call set_weights(grid,predictive,balance_what,refine_control)

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

total_weight = phaml_global_sum(procs,my_total_weight,421)

if (total_weight /= 0.0_my_real) then
   my_weight_fraction = my_total_weight/total_weight
else
   my_weight_fraction = 1.0_my_real
endif

end function my_weight_fraction

!        --------------------
function prior2p_e_regularity(grid,elem)
!        --------------------

!----------------------------------------------------
! This routine computes the regularity constant for the PRIOR2P
! hp-adaptive strategy using the energy norm.
!
! TEMP only using the first eigensolution
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: elem
real(my_real) :: prior2p_e_regularity
!----------------------------------------------------
! Local variables:

real(my_real) :: xvert(VERTICES_PER_ELEMENT), yvert(VERTICES_PER_ELEMENT)
real(my_real), allocatable :: lmat(:,:), lrhs(:), solut(:)
real(my_real) :: rmin, errest2, errest1, au
integer :: degree(EDGES_PER_ELEMENT+1), bmark(EDGES_PER_ELEMENT), &
           edge_type(EDGES_PER_ELEMENT,grid%system_size)
integer, allocatable :: rowdeg(:)
integer :: i, j, nleq, ss, isub, allocstat
!----------------------------------------------------
! Begin executable code

! need p at least 3, so do p refinement if p < 3

if (grid%element(elem)%degree < 3) then
   prior2p_e_regularity = huge(0.0_my_real)
   return
endif

ss = grid%system_size

! count the number of equations associated with this element

nleq = VERTICES_PER_ELEMENT
do i=1,EDGES_PER_ELEMENT
   nleq = nleq + grid%edge(grid%element(elem)%edge(i))%degree-1
end do
nleq = nleq + ((grid%element(elem)%degree-1)*(grid%element(elem)%degree-2))/2
nleq = nleq*ss

allocate(rowdeg(nleq),lmat(nleq,nleq),lrhs(nleq),solut(nleq),stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in prior2p_e_regularity")
   return
endif

! compute the element matrix to be used in computing the energy norm

do i=1,VERTICES_PER_ELEMENT
   xvert(i) = grid%vertex(grid%element(elem)%vertex(i))%coord%x
   yvert(i) = grid%vertex(grid%element(elem)%vertex(i))%coord%y
end do

do i=1,EDGES_PER_ELEMENT
   degree(i) = grid%edge(grid%element(elem)%edge(i))%degree
   edge_type(i,:) = grid%edge_type(grid%element(elem)%edge(i),:)
   bmark(i) = grid%edge(grid%element(elem)%edge(i))%bmark
end do
degree(EDGES_PER_ELEMENT+1) = grid%element(elem)%degree

rmin = 0

call elemental_matrix(grid, xvert, yvert, degree, edge_type, bmark, &
                      grid%system_size, 0, .false., elem, rmin, &
                      "p","a",lmat, lrhs, loc_bconds_s=bconds)

! set the degree and solution associated with each row

isub = 1
do i=1,VERTICES_PER_ELEMENT
   rowdeg(isub:isub+ss-1) = 1
   solut(isub:isub+ss-1) = grid%vertex_solution(grid%element(elem)%vertex(i),:,1)
   isub = isub + ss
end do
do i=1,EDGES_PER_ELEMENT
   do j = 2,degree(i)
      rowdeg(isub:isub+ss-1) = j
      solut(isub:isub+ss-1) = grid%edge(grid%element(elem)%edge(i))%solution(j-1,:,1)
      isub = isub + ss
   end do
end do
do i=3,degree(4)
   do j=1,i-2
      rowdeg(isub:isub+ss-1) = i
      solut(isub:isub+ss-1) = grid%element(elem)%solution(j+((i-2)*(i-3))/2,:,1)
      isub = isub + ss
   end do
end do

! compute the error estimates for subapproximations of degree p-2 and p-1,
! i.e. the energy norm of the parts of the approximation due to bases of
! degree p-1 and p

errest2 = 0
errest1 = 0
do i=1,nleq
   if (rowdeg(i) == degree(4)-1) then
      au = 0
      do j=1,nleq
         if (rowdeg(j) == degree(4)-1) then
            au = au + lmat(i,j)*solut(j)
         endif
      enddo
      errest2 = errest2 + au*solut(i)
   endif
   if (rowdeg(i) == degree(4)) then
      au = 0
      do j=1,nleq
         if (rowdeg(j) == degree(4)) then
            au = au + lmat(i,j)*solut(j)
         endif
      enddo
      errest1 = errest1 + au*solut(i)
   endif
end do
errest2 = sqrt(errest2)
errest1 = sqrt(errest1)

! the regularity comes from ||e_p-1||/||e_p-2|| = ((p-1)/(p-2))^-(m-1)

if (errest1 == 0.0_my_real .or. errest2 == 0.0_my_real) then
   prior2p_e_regularity = degree(4)
else
   prior2p_e_regularity = 1-log(errest1/errest2)/log(real(degree(4)-1)/(degree(4)-2))
endif

deallocate(rowdeg,lmat,lrhs,solut,stat=allocstat)

end function prior2p_e_regularity

!        ---------------------
function prior2p_h1_regularity(grid,elem)
!        ---------------------

!----------------------------------------------------
! This routine computes the regularity constant for the PRIOR2P
! hp-adaptive strategy using the H^1 norm.
!
! TEMP only using the first eigensolution and component
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: elem
real(my_real) :: prior2p_h1_regularity
!----------------------------------------------------
! Local variables:

real(my_real) :: xvert(VERTICES_PER_ELEMENT), yvert(VERTICES_PER_ELEMENT)
real(my_real), allocatable :: vert_soln(:,:,:), edg1_soln(:,:,:), &
                              edg2_soln(:,:,:), edg3_soln(:,:,:), &
                              elem_soln(:,:,:), ux(:,:,:), uy(:,:,:)
real(my_real), pointer :: qweights(:), xquad(:), yquad(:)
real(my_real) :: errest1, errest2
integer :: i, j, deg, astat, nqpoints, jerr
!----------------------------------------------------
! Begin executable code

! need p at least 3, so do p refinement if p < 3

if (grid%element(elem)%degree < 3) then
   prior2p_h1_regularity = huge(0.0_my_real)
   return
endif

! degree of this element, p

deg = grid%element(elem)%degree

! save the solutions associated with this element

allocate(vert_soln(VERTICES_PER_ELEMENT, &
                   size(grid%vertex_solution,dim=2), &
                   size(grid%vertex_solution,dim=3)))
if (associated(grid%edge(grid%element(elem)%edge(1))%solution)) &
allocate(edg1_soln(size(grid%edge(grid%element(elem)%edge(1))%solution,dim=1),&
                   size(grid%edge(grid%element(elem)%edge(1))%solution,dim=2),&
                   size(grid%edge(grid%element(elem)%edge(1))%solution,dim=3)))
if (associated(grid%edge(grid%element(elem)%edge(2))%solution)) &
allocate(edg2_soln(size(grid%edge(grid%element(elem)%edge(2))%solution,dim=1),&
                   size(grid%edge(grid%element(elem)%edge(2))%solution,dim=2),&
                   size(grid%edge(grid%element(elem)%edge(2))%solution,dim=3)))
if (associated(grid%edge(grid%element(elem)%edge(3))%solution)) &
allocate(edg3_soln(size(grid%edge(grid%element(elem)%edge(3))%solution,dim=1),&
                   size(grid%edge(grid%element(elem)%edge(3))%solution,dim=2),&
                   size(grid%edge(grid%element(elem)%edge(3))%solution,dim=3)))
if (associated(grid%element(elem)%solution)) &
allocate(elem_soln(size(grid%element(elem)%solution,dim=1), &
                   size(grid%element(elem)%solution,dim=2), &
                   size(grid%element(elem)%solution,dim=3)))

do i=1,VERTICES_PER_ELEMENT
   vert_soln(i,:,:) = grid%vertex_solution(grid%element(elem)%vertex(i),:,:)
end do
if (associated(grid%edge(grid%element(elem)%edge(1))%solution)) then
   edg1_soln = grid%edge(grid%element(elem)%edge(1))%solution
endif
if (associated(grid%edge(grid%element(elem)%edge(2))%solution)) then
   edg2_soln = grid%edge(grid%element(elem)%edge(2))%solution
endif
if (associated(grid%edge(grid%element(elem)%edge(3))%solution)) then
   edg3_soln = grid%edge(grid%element(elem)%edge(3))%solution
endif
if (associated(grid%element(elem)%solution)) then
   elem_soln = grid%element(elem)%solution
endif

! get a quadrature rule that is exact for polynomial of degree (p-1)^2
! TEMP can we get away with a less accurate integral?

do i=1,VERTICES_PER_ELEMENT
   xvert(i) = grid%vertex(grid%element(elem)%vertex(i))%coord%x
   yvert(i) = grid%vertex(grid%element(elem)%vertex(i))%coord%y
end do

call quadrature_rule_tri(min(MAX_QUAD_ORDER_TRI,(deg-1)**2),xvert,yvert, &
                         nqpoints,qweights,xquad,yquad,jerr,stay_in=.true.)

if (jerr /= 0) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("error returned by quadrature_rule in prior2p_h1_regularity")
   stop
endif
   
! estimate the error in the degree p-2 solution by using only the degree
! p-1 and p parts of the solution, i.e., set all the solution coefficients
! below degree p-1 to 0
! TEMP it has to be much more efficient to just work with the degree p-1
! and degree p parts of the solution rather than setting coefficients to
! 0 and using the evaluate routine.  This approach is just to get something
! programmed quickly to see if the strategy is viable.

do i=1,VERTICES_PER_ELEMENT
   grid%vertex_solution(grid%element(elem)%vertex(i),:,:) = 0
end do
do i=2,deg-2
   do j=1,EDGES_PER_ELEMENT
      if (grid%edge(grid%element(elem)%edge(j))%degree >= i) then
         grid%edge(grid%element(elem)%edge(j))%solution(i-1,:,:) = 0
      endif
   end do
end do
if (deg > 4) then
   do i=1,((deg-4)*(deg-3))/2
      grid%element(elem)%solution(i,:,:) = 0
   end do
endif

! evaluate the derivatives of the solution at the quadrature points

allocate(ux(1,1,nqpoints),uy(1,1,nqpoints),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in prior2p_h1_regularity")
   return
endif

call evaluate_soln_local(grid,xquad,yquad,elem,(/1/),(/1/),ux=ux,uy=uy)

! compute the H^1 norm of the error estimate for the degree p-2 solution
! TEMP does the H^1 norm also include the derivative jumps on the edges?

errest2 = 0
do i=1,nqpoints
   errest2 = errest2 + qweights(i)*(ux(1,1,i)**2+uy(1,1,i)**2)
end do
errest2 = sqrt(errest2)

! now compute the estimate of the error in the degree p-1 approximations
! using only the degree p part fo the solution

do i=1,EDGES_PER_ELEMENT
   if (grid%edge(grid%element(elem)%edge(i))%degree >= deg-1) then
      grid%edge(grid%element(elem)%edge(i))%solution(deg-2,:,:) = 0
   endif
end do

if (deg > 3) then
   do i=1+((deg-4)*(deg-3))/2, ((deg-3)*(deg-2))/2
      grid%element(elem)%solution(i,:,:) = 0
   end do
endif

call evaluate_soln_local(grid,xquad,yquad,elem,(/1/),(/1/),ux=ux,uy=uy)

errest1 = 0
do i=1,nqpoints
   errest1 = errest1 + qweights(i)*(ux(1,1,i)**2+uy(1,1,i)**2)
end do
errest1 = sqrt(errest1)

! the regularity comes from ||e_p-1||/||e_p-2|| = ((p-1)/(p-2))^-(m-1)

if (errest1 == 0.0_my_real .or. errest2 == 0.0_my_real) then
   prior2p_h1_regularity = deg
else
   prior2p_h1_regularity = 1-log(errest1/errest2)/log(real(deg-1)/(deg-2))
endif

! restore the solution

do i=1,VERTICES_PER_ELEMENT
   grid%vertex_solution(grid%element(elem)%vertex(i),:,:) = vert_soln(i,:,:)
end do
if (associated(grid%edge(grid%element(elem)%edge(1))%solution)) &
   grid%edge(grid%element(elem)%edge(1))%solution = edg1_soln
if (associated(grid%edge(grid%element(elem)%edge(2))%solution)) &
   grid%edge(grid%element(elem)%edge(2))%solution = edg2_soln
if (associated(grid%edge(grid%element(elem)%edge(3))%solution)) &
   grid%edge(grid%element(elem)%edge(3))%solution = edg3_soln
if (associated(grid%element(elem)%solution)) &
   grid%element(elem)%solution = elem_soln

! free memory

deallocate(vert_soln, qweights, xquad, yquad, ux, uy, stat=astat)
if (allocated(edg1_soln)) deallocate(edg1_soln)
if (allocated(edg2_soln)) deallocate(edg2_soln)
if (allocated(edg3_soln)) deallocate(edg3_soln)
if (allocated(elem_soln)) deallocate(elem_soln)

end function prior2p_h1_regularity

!        -----------------
function next3p_regularity(grid,elem)
!        -----------------

!----------------------------------------------------
! This routine computes the regularity for the NEXT3P strategy.
! TEMP It will be much more efficient to compute the equilibrated residual
!      error estimates for the whole grid at once, when possible.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: elem
real(my_real) :: next3p_regularity
!----------------------------------------------------
! Local variables:

integer :: deg, i
double precision :: a, b, t, alpha
! only need C2 and normpsi2 for debugging
!double precision :: normpsi2, t1, t2, C2
!----------------------------------------------------
! Begin executable code

! compute the three error estimates

call equilibrated_residual_ei_one(grid,elem,multi_p=n3pee)

deg = grid%element(elem)%degree
npq(1) = deg+1
npq(2) = deg+2
npq(3) = deg+3

alpha = -1
a = .001d0
b = 0
do i=1,10
   b = b+1
   t = next3pfunc(b)
   if (next3pfunc(a)*next3pfunc(b) < 0.0d0) exit
end do
if (next3pfunc(a)*next3pfunc(b) > 0.0d0) then
!   print *,elem," couldn't bracket solution"
   alpha = 0.0d0
endif

if (alpha /= 0.0d0) then
   t = 1.0d-13
   alpha = dzero(next3pfunc,a,b,t)
! only need C2 and normpsi2 for debugging
!   t1 = n3pee(3,1)**2 - n3pee(1,1)**2
!   t2 = npq(3)**(-2*alpha) - npq(1)**(-2*alpha)
!   C2 = -t1/t2
!   normpsi2 = C2*npq(1)**(-2*alpha) + n3pee(1,1)**2
   next3p_regularity = 1+alpha
else
   next3p_regularity = huge(0.0_my_real)
endif

end function next3p_regularity

!        ----------
function next3pfunc(x)
!        ----------

!----------------------------------------------------
! The function who's root gives alpha for next3p
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x
real(my_real) :: next3pfunc
!----------------------------------------------------
! Local variables:

real(my_real) :: t1, t2, t3, t4
!----------------------------------------------------
! Begin executable code

t1 = n3pee(3,1)**2 - n3pee(1,1)**2
t2 = npq(3)**(-2*x) - npq(1)**(-2*x)
t3 = n3pee(2,1)**2 - n3pee(1,1)**2
t4 = npq(2)**(-2*x) - npq(1)**(-2*x)

next3pfunc = t1/t2 - t3/t4

end function next3pfunc

!              -----
      FUNCTION DZERO(F,A,B,T)
!              -----
      REAL(MY_REAL) DZERO

!----------------------------------------------------
! A slight modification (mostly for free format) of DZERO from netlib
!----------------------------------------------------
!
!  FINDS THE REAL ROOT OF THE FUNCTION F LYING BETWEEN A AND B
!  TO WITHIN A TOLERANCE OF
!
!         6*D1MACH(3) * DABS(DZERO) + 2 * T
!
!  F(A) AND F(B) MUST HAVE OPPOSITE SIGNS
!
!  THIS IS BRENTS ALGORITHM
!
!  A, STORED IN SA, IS THE PREVIOUS BEST APPROXIMATION (I.E. THE OLD B)
!  B, STORED IN SB, IS THE CURRENT BEST APPROXIMATION
!  C IS THE MOST RECENTLY COMPUTED POINT SATISFYING F(B)*F(C) .LT. 0
!  D CONTAINS THE CORRECTION TO THE APPROXIMATION
!  E CONTAINS THE PREVIOUS VALUE OF D
!  M CONTAINS THE BISECTION QUANTITY (C-B)/2
!
      REAL(MY_REAL) F,A,B,T,TT,SA,SB,C,D,E,FA,FB,FC,TOL,M,P,Q,R,S
      EXTERNAL F
      DOUBLE PRECISION D1MACH

      TT = T
      IF (T .LE. 0.0D0) TT = 10.D0*D1MACH(1)

      SA = A
      SB = B
      FA = F(SA)
      FB = F(SB)
      IF (FA .NE. 0.0D0) GO TO 5
      DZERO = SA
      RETURN
  5   IF (FB .EQ. 0.0D0) GO TO 140
      IF (DSIGN(FA,FB) .EQ. FA) print *, ' DZERO - F(A) AND F(B) ARE NOT OF OPPOSITE SIGN'

 10   C  = SA
      FC = FA
      E  = SB-SA
      D  = E

!  INTERCHANGE B AND C IF DABS F(C) .LT. DABS F(B)

 20   IF (DABS(FC).GE.DABS(FB)) GO TO 30
      SA = SB
      SB = C
      C  = SA
      FA = FB
      FB = FC
      FC = FA

 30   TOL = 2.0D0*D1MACH(4)*DABS(SB)+TT
      M = 0.5D0*(C-SB)

!  SUCCESS INDICATED BY M REDUCES TO UNDER TOLERANCE OR
!  BY F(B) = 0

      IF ((DABS(M).LE.TOL).OR.(FB.EQ.0.0D0)) GO TO 140

!  A BISECTION IS FORCED IF E, THE NEXT-TO-LAST CORRECTION
!  WAS LESS THAN THE TOLERANCE OR IF THE PREVIOUS B GAVE
!  A SMALLER F(B).  OTHERWISE GO TO 40.

      IF ((DABS(E).GE.TOL).AND.(DABS(FA).GE.DABS(FB))) GO TO 40
      E = M
      D = E
      GO TO 100
 40   S = FB/FA

!  QUADRATIC INTERPOLATION CAN ONLY BE DONE IF A (IN SA)
!  AND C ARE DIFFERENT POINTS.
!  OTHERWISE DO THE FOLLOWING LINEAR INTERPOLATION

      IF (SA.NE.C) GO TO 50
      P = 2.0D0*M*S
      Q = 1.0D0-S
      GO TO 60

!  INVERSE QUADRATIC INTERPOLATION

 50   Q = FA/FC
      R = FB/FC
      P = S*(2.0D0*M*Q*(Q-R)-(SB-SA)*(R-1.0D0))
      Q = (Q-1.0D0)*(R-1.0D0)*(S-1.0D0)
 60   IF (P.LE.0.0D0) GO TO 70
      Q = -Q
      GO TO 80
 70   P = -P

!  UPDATE THE QUANTITIES USING THE NEWLY COMPUTED
!  INTERPOLATE UNLESS IT WOULD EITHER FORCE THE
!  NEW POINT TOO FAR TO ONE SIDE OF THE INTERVAL
!  OR WOULD REPRESENT A CORRECTION GREATER THAN
!  HALF THE PREVIOUS CORRECTION.

!  IN THESE LAST TWO CASES - DO THE BISECTION
!  BELOW (FROM STATEMENT 90 TO 100)

 80   S = E
      E = D
      IF ((2.0D0*P.GE.3.0D0*M*Q-DABS(TOL*Q)).OR.(P.GE.DABS(0.5D0*S*Q))) GO TO 90
      D = P/Q
      GO TO 100
 90   E = M
      D = E

!  SET A TO THE PREVIOUS B

 100  SA = SB
      FA = FB

!  IF THE CORRECTION TO BE MADE IS SMALLER THAN
!  THE TOLERANCE, JUST TAKE A  DELTA STEP  (DELTA=TOLERANCE)
!         B = B + DELTA * SIGN(M)

      IF (DABS(D).LE.TOL) GO TO 110
      SB = SB+D
      GO TO 130

 110  IF (M.LE.0.0D0) GO TO 120
      SB = SB+TOL
      GO TO 130

 120  SB = SB-TOL
 130  FB = F(SB)

!  IF F(B) AND F(C) HAVE THE SAME SIGN ONLY
!  LINEAR INTERPOLATION (NOT INVERSE QUADRATIC)
!  CAN BE DONE

      IF ((FB.GT.0.0D0).AND.(FC.GT.0.0D0)) GO TO 10
      IF ((FB.LE.0.0D0).AND.(FC.LE.0.0D0)) GO TO 10
      GO TO 20

!***SUCCESS***
 140  DZERO = SB
      RETURN
      END FUNCTION DZERO

!        -----------
function typeparam_R(grid,elem)
!        -----------

!----------------------------------------------------
! This routine computes the ratio e(t,p)/e(t,p-1) for the type parameter
! comparison for element elem
! TEMP only using first eigenvector
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: elem
real(my_real) :: typeparam_R
!----------------------------------------------------
! Local variables:

real(my_real) :: xvert(VERTICES_PER_ELEMENT), yvert(VERTICES_PER_ELEMENT)
real(my_real), allocatable :: lmat(:,:), lrhs(:), solut(:)
real(my_real) :: rmin, errest2(1), errest1, au
integer :: degree(EDGES_PER_ELEMENT+1), bmark(EDGES_PER_ELEMENT), &
           edge_type(EDGES_PER_ELEMENT,grid%system_size)
integer, allocatable :: rowdeg(:)
integer :: i, j, nleq, ss, isub, allocstat
!----------------------------------------------------
! Begin executable code

! if the degree is 1, then R is 0 to force p refinement

if (grid%element(elem)%degree == 1) then
   typeparam_R = 0.0_my_real

else

! compute the error estimate for degree p-1.  It is the norm of the p'th
! degree part of the solution

   ss = grid%system_size

! count the number of equations associated with this element

   nleq = VERTICES_PER_ELEMENT
   do i=1,EDGES_PER_ELEMENT
      nleq = nleq + grid%edge(grid%element(elem)%edge(i))%degree-1
   end do
   nleq = nleq + ((grid%element(elem)%degree-1)*(grid%element(elem)%degree-2))/2
   nleq = nleq*ss

   allocate(rowdeg(nleq),lmat(nleq,nleq),lrhs(nleq),solut(nleq),stat=allocstat)
   if (allocstat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in typeparam_R")
      return
   endif

! compute the element matrix to be used in computing the energy norm

   do i=1,VERTICES_PER_ELEMENT
      xvert(i) = grid%vertex(grid%element(elem)%vertex(i))%coord%x
      yvert(i) = grid%vertex(grid%element(elem)%vertex(i))%coord%y
   end do

   do i=1,EDGES_PER_ELEMENT
      degree(i) = grid%edge(grid%element(elem)%edge(i))%degree
      edge_type(i,:) = grid%edge_type(grid%element(elem)%edge(i),:)
      bmark(i) = grid%edge(grid%element(elem)%edge(i))%bmark
   end do
   degree(EDGES_PER_ELEMENT+1) = grid%element(elem)%degree

   rmin = 0

! TEMP should this be just "r"?

   call elemental_matrix(grid, xvert, yvert, degree, edge_type, bmark, &
                         grid%system_size, 0, .false., elem, rmin, &
                         "p","a",lmat, lrhs, loc_bconds_s=bconds)

! set the degree and solution associated with each row

   isub = 1
   do i=1,VERTICES_PER_ELEMENT
      rowdeg(isub:isub+ss-1) = 1
      solut(isub:isub+ss-1) = grid%vertex_solution(grid%element(elem)%vertex(i),:,1)
      isub = isub + ss
   end do
   do i=1,EDGES_PER_ELEMENT
      do j = 2,degree(i)
         rowdeg(isub:isub+ss-1) = j
         solut(isub:isub+ss-1) = grid%edge(grid%element(elem)%edge(i))%solution(j-1,:,1)
         isub = isub + ss
      end do
   end do
   do i=3,degree(4)
      do j=1,i-2
         rowdeg(isub:isub+ss-1) = i
         solut(isub:isub+ss-1) = grid%element(elem)%solution(j+((i-2)*(i-3))/2,:,1)
         isub = isub + ss
      end do
   end do

! compute the error estimates for subapproximations of degree p-1,
! i.e. the energy norm of the part of the approximation due to bases of
! degree p

   errest1 = 0
   do i=1,nleq
      if (rowdeg(i) == degree(4)) then
         au = 0
         do j=1,nleq
            if (rowdeg(j) == degree(4)) then
               au = au + lmat(i,j)*solut(j)
            endif
         enddo
         errest1 = errest1 + au*solut(i)
      endif
   end do
   errest1 = sqrt(errest1)

! if the p-1 error estimate is 0, then R is 0

   if (errest1 == 0.0_my_real) then
      typeparam_R = 0.0_my_real

   else

! otherwise compute the p error estimate as a local Neuman error estimate

      call error_indicator(grid,elem,LOCAL_PROBLEM_P,energy=errest2)
      typeparam_R = errest2(1)/errest1

   endif
endif ! p==1

end function typeparam_R

!        --------------------
function coef_root_regularity(grid,elem)
!        --------------------

!----------------------------------------------------
! This routine estimates the regularity by the root test on the basis coefs
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: elem
real(my_real) :: coef_root_regularity
!----------------------------------------------------
! Local variables:

! RESTRICTION max_deg <= 50
real(my_real) :: lcoefs(50)

real(my_real) :: a_i
integer :: i, deg
!----------------------------------------------------
! Begin executable code

! if the degree is 1, log(deg) in the denominator doesn't work, so
! return 10 to force p refinement

if (grid%element(elem)%degree == 1) then
   coef_root_regularity = 10.0_my_real

else

! get the basis coefficients

   call basis_coefs(grid,elem,lcoefs)

! find the last nonzero basis coefficients

   a_i   = -1.0_my_real
   do i = grid%element(elem)%degree,1,-1
      if (lcoefs(i) /= 0.0_my_real) then
         a_i = lcoefs(i)
         exit
      endif
   end do

! if we didn't find a nonzero coefficients, force p refinement

   if (a_i < 0.0_my_real) then
      coef_root_regularity = huge(0.0_my_real)

   else

! compute the regularity estimate

      deg = grid%element(elem)%degree
      coef_root_regularity = log((2*deg+1)/(2*a_i**2))/(2*log(real(deg,my_real))) - 0.5_my_real

   endif
endif

end function coef_root_regularity

!        ---------------
function coef_decay_rate(grid,elem)
!        ---------------

!----------------------------------------------------
! This routine estimates the decay rate of the basis coefficients.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: elem
real(my_real) :: coef_decay_rate
!----------------------------------------------------
! Local variables:

! RESTRICTION max_deg <= 50
real(my_real) :: lcoefs(50)

integer, parameter :: nfit = 4
real(my_real) :: sx, sx2, sy, sxy
integer :: i, s, e
!----------------------------------------------------
! Begin executable code

! if the degree is 1, there's not enough data to fit the exponential so
! return 10 to force p refinement

if (grid%element(elem)%degree == 1) then
   coef_decay_rate = 10.0_my_real

else

! get the basis coefficients

   call basis_coefs(grid,elem,lcoefs)

! fit the last nfit points to c*exp(-sigma*p)

   e = grid%element(elem)%degree
   s = max(1,e-nfit+1)

   sx  = 0.0_my_real
   sx2 = 0.0_my_real
   sy  = 0.0_my_real
   sxy = 0.0_my_real
   do i=s,e
      sx  = sx + i
      sx2 = sx2 + i*i
      sy  = sy + log(lcoefs(i))
      sxy = sxy + i*log(lcoefs(i))
   end do

   coef_decay_rate = - ((e-s+1)*sxy - sx*sy) / ((e-s+1)*sx2 - sx*sx)

endif

end function coef_decay_rate

!          -----------
subroutine basis_coefs(grid,elem,lcoefs)
!          -----------

!----------------------------------------------------
! This routine puts a measure of the size of the basis coefficients into lcoefs.
! lcoefs(i) is the sum of the absolute value of the coefficients of bases of
! degree i in elem.
! TEMP only using first component of first eigenvector
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: elem
real(my_real), intent(out) :: lcoefs(:)
!----------------------------------------------------
! Local variables:

integer :: i,j,deg
!----------------------------------------------------
! Begin executable code

lcoefs=0.0_my_real

! vertex bases

do i=1,3
   lcoefs(1)=lcoefs(1) + &
             abs(grid%vertex_solution(grid%element(elem)%vertex(i),1,1))
end do

! edge bases

do i=1,grid%element(elem)%degree-1
   do j=1,3
      if (grid%edge(grid%element(elem)%edge(j))%degree >= i+1) then
         lcoefs(i+1)=lcoefs(i+1) + &
                     abs(grid%edge(grid%element(elem)%edge(j))%solution(i,1,1))
      endif
   end do
end do

! face bases

do i=1,((grid%element(elem)%degree-2)*(grid%element(elem)%degree-1))/2
   deg = ceiling((3+sqrt(8.0*i+1))/2)
   lcoefs(deg)=lcoefs(deg) + abs(grid%element(elem)%solution(i,1,1))
end do

end subroutine basis_coefs

!          ----------
subroutine refine_nlp(grid,procs,refine_control,still_sequential)
!          ----------

!----------------------------------------------------
! This routine performs one refinement phase under the NLP hp strategy
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(proc_info), intent(in) :: procs
type(refine_options), intent(in) :: refine_control
logical, intent(in) :: still_sequential
!----------------------------------------------------
! Local variables:

integer :: i, lev, elem, maxdeg
integer, allocatable :: new_level(:), new_degree(:)
real(my_real), allocatable :: new_params(:)
logical :: changed
!----------------------------------------------------
! Begin executable code

! TEMP not parallel yet

if (num_proc(procs) > 1) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("NLP not yet implemented in parallel")
   stop
endif

! algencan requires double precision

if (my_real /= kind(0.0d0)) then
   ierr = USER_INPUT_ERROR
   call fatal("NLP requires my_real to be double precision")
   stop
endif

! MASTER doesn't participate

if (my_proc(procs) == MASTER) return

! start timing the refinement process

call reset_watch(prefine)
call start_watch((/prefine,trefine/))

! initializations for the NLP hp-adaptive strategy

call nlp_initialize(grid,refine_control)

! make sure tau is such that we can get a feasible solution by making sure
! it is at least 1.5 times the error estimate of a maximally refined grid

allocate(new_params(2*nlp_nelem))
do i=1,nlp_nelem
   new_params(i) = grid%element(nlp_elem_list(i))%level + &
                   refine_control%nlp_max_h_inc
   new_params(i+nlp_nelem) = grid%element(nlp_elem_list(i))%degree + &
                             refine_control%nlp_max_p_inc
end do

nlp_tau = max(nlp_tau,1.5_my_real*nlp_errest(new_params))

! believe that the algorithm has been doing the right thing so far, and
! set the initial guess such that elements that have been more h and/or p
! refined are more refined that way in the initial guess

maxdeg = 0
do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (grid%element(elem)%iown .and. grid%element(elem)%isleaf) then
         maxdeg=max(maxdeg,grid%element(elem)%degree)
      endif
      elem = grid%element(elem)%next
   end do
end do
do i=1,nlp_nelem
   new_params(i) = grid%element(nlp_elem_list(i))%level + max(0, &
                   refine_control%nlp_max_h_inc - &
                   (grid%nlev - grid%element(nlp_elem_list(i))%level))
   new_params(i+nlp_nelem) = grid%element(nlp_elem_list(i))%degree + max(0, &
                             refine_control%nlp_max_p_inc - &
                             (maxdeg - grid%element(nlp_elem_list(i))%degree))
end do

! call the optimization routine

call algencanma(new_params,nlp_tau/100)

! integerize the new parameters

allocate(new_level(nlp_nelem),new_degree(nlp_nelem))
do i=1,nlp_nelem
   new_level(i) = new_params(i) + 0.5_my_real
   new_degree(i) = new_params(i+nlp_nelem) + 0.5_my_real
end do

! if no change, uniform refinement to avoid stalling

changed = .false.
do i=1,nlp_nelem
   if (new_level(i) /= grid%element(nlp_elem_list(i))%level) changed = .true.
   if (new_degree(i) /= grid%element(nlp_elem_list(i))%degree) changed = .true.
   if (changed) exit
end do
if (.not. changed) then
   new_level = new_level + 1
   new_degree = new_degree + 1
endif

! perform the refinements

call nlp_refine_grid(grid,new_level,new_degree,refine_control)

! free memory

deallocate(new_params,new_level,new_degree)

! finish the NLP hp-adaptive strategy

call nlp_finalize()

! error indicators should be recomputed

grid%errind_up2date = .false.

! stop timing the refinement process

call stop_watch((/prefine,trefine/))

end subroutine refine_nlp

!          --------------
subroutine nlp_initialize(grid,refine_control)
!          --------------

!----------------------------------------------------
! This routine performs initializations for the NLP hp-adaptive strategy
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout), target :: grid
type(refine_options), intent(in), target :: refine_control
!----------------------------------------------------
! Local variables:

integer :: i
!----------------------------------------------------
! Begin executable code

! create a list of the leaf elements.  This is the order in which the
! candidate level and degree parameters will be listed.

call create_nlp_elem_list(grid)

! we will need to use grid and refine_control in the objective function and
! constraints, but can't pass it through the NLP software

hold_grid => grid
hold_refcont => refine_control

! make sure the error indicators are up to date

if (.not. grid%errind_up2date) then
   call all_error_indicators(grid,refine_control%error_estimator)
endif

! set the error tolerance to 1/4 the current error estimate
! TEMP only using first eigenvector

nlp_tau = 0
do i=1,nlp_nelem
   nlp_tau = nlp_tau + grid%element_errind(nlp_elem_list(i),1)**2
end do
nlp_tau = nlp_tau/16

! approximate the smoothness the same as in the NEXT3P strategy, but don't
! allow m<1.5 because it can cause stalling because of a large tau, or
! bigger than 10.  If the degree is 1 then next3p return infinity, so use
! 1.5 in that case.

allocate(nlp_m(nlp_nelem))
do i=1,nlp_nelem
   nlp_m(i) = max(1.5_my_real,next3p_regularity(grid,nlp_elem_list(i)))
   if (grid%element(nlp_elem_list(i))%degree == 1) nlp_m(i) = 1.5_my_real
   if (nlp_m(i) > 10.0_my_real) nlp_m(i) = 10
end do

end subroutine nlp_initialize

!          ------------
subroutine nlp_finalize()
!          ------------

!----------------------------------------------------
! This routine finishes the NLP hp-adaptive strategy
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

! deallocate the nlp_elem_list

call destroy_nlp_elem_list

! free memory from m

deallocate(nlp_m)

end subroutine nlp_finalize

!          --------------------
subroutine create_nlp_elem_list(grid)
!          --------------------

!----------------------------------------------------
! This routine creates the list of leaf elements in nlp_elem_list
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
!----------------------------------------------------
! Local variables:

integer :: lev, elem
!----------------------------------------------------
! Begin executable code

allocate(nlp_elem_list(size(grid%element)),nlp_inverse_elem(size(grid%element)))

nlp_nelem = 0
do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (grid%element(elem)%isleaf) then
         nlp_nelem = nlp_nelem+1
         nlp_elem_list(nlp_nelem) = elem
         nlp_inverse_elem(elem) = nlp_nelem
      endif
      elem = grid%element(elem)%next
   end do
end do

end subroutine create_nlp_elem_list

!          ---------------------
subroutine destroy_nlp_elem_list()
!          ---------------------

!----------------------------------------------------
! This routine frees the memory of nlp_elem_list
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

deallocate(nlp_elem_list,nlp_inverse_elem)

end subroutine destroy_nlp_elem_list

!        -------
function nlp_dof(candidate_params)
!        -------

!----------------------------------------------------
! This routine computes the (approximate) number of freedom in the grid
! with the candidate refinement for NLP
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: candidate_params(:)
real(my_real) :: nlp_dof
!----------------------------------------------------
! Local variables:

integer :: i, elem, l
real(my_real) :: lhat, phat
!----------------------------------------------------
! Begin executable code

nlp_dof = 0.0_my_real
do i=1,nlp_nelem
   elem = nlp_elem_list(i)
   l = hold_grid%element(elem)%level
   lhat = candidate_params(i)
   phat = candidate_params(i+nlp_nelem)
   nlp_dof = nlp_dof + 2**(lhat-l) * phat**2
end do

nlp_dof = nlp_dof/2

end function nlp_dof

!          --------
subroutine nlp_ddof(candidate_params,ddof)
!          --------

!----------------------------------------------------
! This routine computes the derivative of the dof function with respect to
! each of the parameters, and returns them in ddof.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: candidate_params(:)
real(my_real), intent(out) :: ddof(:)
!----------------------------------------------------
! Local variables:

integer :: i, elem, l
real(my_real) :: lhat, phat
!----------------------------------------------------
! Begin executable code

do i=1,nlp_nelem
   elem = nlp_elem_list(i)
   l = hold_grid%element(elem)%level
   lhat = candidate_params(i)
   phat = candidate_params(i+nlp_nelem)
   ddof(i) = log(2.0_my_real) * 2**(lhat-l) * phat**2 / 2
   ddof(i+nlp_nelem) = 2**(lhat-l) * phat
end do

end subroutine nlp_ddof

!        ----------
function nlp_errest(candidate_params)
!        ----------

!----------------------------------------------------
! This routine computes the square of the error estimate for NLP
! TEMP only using first eigenvalue
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: candidate_params(:)
real(my_real) :: nlp_errest
!----------------------------------------------------
! Local variables:

integer :: i, elem, l
real(my_real) :: m, p, eta, lhat, phat
!----------------------------------------------------
! Begin executable code

nlp_errest = 0.0_my_real
do i=1,nlp_nelem
   elem = nlp_elem_list(i)

! convenience names

   l = hold_grid%element(elem)%level
   lhat = candidate_params(i)
   p = hold_grid%element(elem)%degree
   phat = candidate_params(i+nlp_nelem)
   eta = hold_grid%element_errind(elem,1)
   m = nlp_m(i)

! contribution of this element to the error estimate

   nlp_errest = nlp_errest + &
                2**(min(p,m-1)*(l-lhat)) * (p/phat)**(2*(m-1)) * eta**2

end do

nlp_errest=nlp_errest/nlp_tau

end function nlp_errest

!          -----------
subroutine nlp_derrest(candidate_params,derrest)
!          -----------

!----------------------------------------------------
! This routine computes the derivative of the errest function with respect to
! each of the parameters, and returns them in derrest.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: candidate_params(:)
real(my_real), intent(out) :: derrest(:)
!----------------------------------------------------
! Local variables:

integer :: i, elem, l
real(my_real) :: m, p, eta, lhat, phat
!----------------------------------------------------
! Begin executable code

do i=1,nlp_nelem
   elem = nlp_elem_list(i)

! convenience names

   l = hold_grid%element(elem)%level
   lhat = candidate_params(i)
   p = hold_grid%element(elem)%degree
   phat = candidate_params(i+nlp_nelem)
   eta = hold_grid%element_errind(elem,1)
   m = nlp_m(i)

! derivatives w.r.t. the l and p of this element

   derrest(i) = -log(2.0_my_real)*min(p,m-1)*2**(min(p,m-1)*(l-lhat)) * &
                    (p/phat)**(2*(m-1)) * eta**2
   derrest(i+nlp_nelem) = -2*(m-1) * &
                2**(min(p,m-1)*(l-lhat)) * (p/phat)**(2*(m-1)) * eta**2 / phat

end do

derrest = derrest/nlp_tau

end subroutine nlp_derrest

!          ---------------
subroutine nlp_refine_grid(grid,new_level,new_degree,refine_control)
!          ---------------

!----------------------------------------------------
! This routine performs the refinements determined to be optimal by NLP
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: new_level(:), new_degree(:)
type(refine_options), intent(in) :: refine_control
!----------------------------------------------------
! Local variables:

integer :: i, j, elem, errcode, nrefine(nlp_nelem), parent, elem_check
type(hash_key) :: elem_gid, parent_gid, list_gid(nlp_nelem)
!----------------------------------------------------
! Begin executable code


! first perform all the p coarsening/refinement.  The correct p will then
! propagate to the parents/children during h coarsening/refinement.
! Also, save a list of the gids for use in h refinement.

do i=1,nlp_nelem
   elem = nlp_elem_list(i)
   list_gid(i) = grid%element(elem)%gid
   do while (grid%element(elem)%degree > new_degree(i))
      call p_coarsen_elem(grid,elem,errcode,refine_control)
      if (errcode /= 0) then
         ierr = PHAML_INTERNAL_ERROR
         call fatal("could not p coarsen in nlp_refine_grid")
         stop
      endif
   end do
   do while (grid%element(elem)%degree < new_degree(i))
      call p_refine_elem(grid,elem,refine_control)
   end do
end do

! determine how many times elements should be h refined or coarsened before
! we lose information to coarsening

do i=1,nlp_nelem
   nrefine(i) = new_level(i) - grid%element(nlp_elem_list(i))%level
end do

! next perform all the h coarsenings.  Note that if an element is not found
! in the element hash table, then it must have already been coarsened as the
! mate of some other element that was coarsened.

do i=1,nlp_nelem
   elem = nlp_elem_list(i)
   elem_gid = grid%element(elem)%gid
   do j=1,-nrefine(i)
      parent_gid = elem_gid/2
      elem_check = hash_decode_key(elem_gid,grid%elem_hash)
      if (elem_check /= HASH_NOT_FOUND) then
         parent = hash_decode_key(parent_gid,grid%elem_hash)
         if (parent == HASH_NOT_FOUND) then
            ierr = PHAML_INTERNAL_ERROR
            call fatal("didn't find parent in hash table in nlp_refine_grid")
            stop
         endif
         call unbisect_triangle_pair(grid,parent,errcode,refine_control)
      endif
      elem_gid = parent_gid
   end do
end do

! finally, perform h refinements.  Since elements on nlp_elem_list may no
! longer exist, get their local id from the saved global id and create
! them if necessary

do i=1,nlp_nelem
   if (nrefine(i) > 0) then
      elem = hash_decode_key(list_gid(i),grid%elem_hash)
      if (elem == HASH_NOT_FOUND) then
         call create_element(grid,list_gid(i),refine_control,errcode)
         elem = hash_decode_key(list_gid(i),grid%elem_hash)
      endif

! then use a recursive routine to perform nrefine refinements making sure the
! whole binary tree of descendants is created

      call nlp_refine_grid_recur(grid,elem,refine_control,nrefine(i))
   endif
end do

end subroutine nlp_refine_grid

!                    ---------------------
recursive subroutine nlp_refine_grid_recur(grid,elem,refine_control,nrefine)
!                    ---------------------

!----------------------------------------------------
! This routine h-refines element elem nrefine times
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: elem
type(refine_options), intent(in) :: refine_control
integer, intent(in) :: nrefine
!----------------------------------------------------
! Local variables:

integer :: child(MAX_CHILD), errcode
!----------------------------------------------------
! Begin executable code

if (nrefine > 0) then

! refine this element

   if (grid%element(elem)%isleaf) then
      call bisect_triangle_pair(grid,elem,errcode,refine_control)
   endif

! refine the children

   if (nrefine > 1) then
      child = get_child_lid(grid%element(elem)%gid,ALL_CHILDREN,grid%elem_hash)
      if (child(1) /= NO_CHILD) then
         call nlp_refine_grid_recur(grid,child(1),refine_control,nrefine-1)
         call nlp_refine_grid_recur(grid,child(2),refine_control,nrefine-1)
      endif
   endif
endif

end subroutine nlp_refine_grid_recur

! interface to algencan


!     ******************************************************************
!     ******************************************************************

      subroutine algencanma(new_params,newepsfeas)

      implicit none

      real(my_real), intent(out) :: new_params(:)
      real(my_real), intent(in) :: newepsfeas

!     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )

!     LOCAL SCALARS
      logical checkder
      integer inform,iprint,m,n,ncomp
      double precision cnorm,epsfeas,epsopt,f,nlpsupn,snorm

!     LOCAL ARRAYS
      logical coded(10),equatn(mmax),linear(mmax)
      double precision l(nmax),lambda(mmax),u(nmax),x(nmax)

!     EXTERNAL SUBROUTINES
!      external algencan,param

!     SET UP PROBLEM DATA

      call param(epsfeas,epsopt,iprint,ncomp)
      epsfeas = newepsfeas

      call inip(n,x,l,u,m,lambda,equatn,linear,coded,checkder,new_params)

      call algencan(epsfeas,epsopt,iprint,ncomp,n,x,l,u,m,lambda,equatn, &
      linear,coded,checkder,f,cnorm,snorm,nlpsupn,inform)

      call endp(n,x,l,u,m,lambda,equatn,linear)

      new_params(1:2*nlp_nelem) = x(1:2*nlp_nelem)

      end subroutine

!     ******************************************************************
!     ******************************************************************

      subroutine param(epsfeas,epsopt,iprint,ncomp)

!     SCALAR ARGUMENTS
      integer iprint,ncomp
      double precision epsfeas,epsopt

      epsfeas  = 1.0d-12
      epsopt   = 1.0d-2

      iprint   = 10
      ncomp    = 6

      end subroutine

!          -----------
subroutine set_weights(grid,predictive,balance_what,refcont)
!          -----------

!----------------------------------------------------
! This routine sets the weights for each element based on what is to be
! balanced and whether or not the balancing will be done before refinement
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
logical, intent(in) :: predictive
integer, intent(in) :: balance_what
type (refine_options), intent(in) :: refcont
!----------------------------------------------------
! Local variables:

integer :: i, lev, vert, elem, astat, edge
real(my_real), allocatable :: assoc_verts(:)
real(my_real) :: energy_errind(grid%system_size), work, err, err2
!----------------------------------------------------
! Begin executable code

! count number of associated vertices for each element if balancing vertices,
! or number of associated equations if balancing equations

if (balance_what==BALANCE_VERTICES .or. balance_what==BALANCE_EQUATIONS) then

   allocate(assoc_verts(size(grid%element)),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in set_weights")
      return
   endif

   assoc_verts = 0

   if (balance_what == BALANCE_EQUATIONS) then

      do lev=1,grid%nlev
         elem = grid%head_level_elem(lev)
         do while (elem /= END_OF_LIST)
            if (grid%element(elem)%isleaf) then
! equations associated with the edges
               do i=1,EDGES_PER_ELEMENT
                  edge = grid%element(elem)%edge(i)
                  vert = grid%edge(edge)%vertex(2)
                  if (grid%edge_type(edge,1) == INTERIOR) then
                     assoc_verts(grid%vertex(vert)%assoc_elem) = &
                           assoc_verts(grid%vertex(vert)%assoc_elem) + &
                           (grid%edge(edge)%degree-1)/2.0_my_real
                  else
                     assoc_verts(grid%vertex(vert)%assoc_elem) = &
                           assoc_verts(grid%vertex(vert)%assoc_elem) + &
                           grid%edge(edge)%degree - 1
                  endif
               end do
! equations associated with the elements
               if (grid%element(elem)%degree > 2) then
                  assoc_verts(elem) = assoc_verts(elem) + &
                 ((grid%element(elem)%degree-1)*(grid%element(elem)%degree-2))/2
               endif
            endif
            elem = grid%element(elem)%next
         end do
      end do

   endif

! the vertices, or equations associated with them

   do lev=1,grid%nlev
      vert = grid%head_level_vert(lev)
      do while (vert /= END_OF_LIST)
         assoc_verts(grid%vertex(vert)%assoc_elem) = &
            assoc_verts(grid%vertex(vert)%assoc_elem) + 1
         vert = grid%vertex(vert)%next
      end do
   end do

endif

! if predictive make sure the error indicators are up to date

if (predictive) then
   if (refcont%reftype==HP_ADAPTIVE .and. &
       (refcont%hp_strategy==HP_REFSOLN_EDGE .or. &
        refcont%hp_strategy==HP_REFSOLN_ELEM)) then
      call all_error_indicators(grid,EXPLICIT_ERRIND)
   elseif (.not. grid%errind_up2date) then
      call all_error_indicators(grid,refcont%error_estimator)
   endif
endif

do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)

! set weight based on what is to be balanced

      select case (balance_what)

      case (BALANCE_NONE)

         grid%element(elem)%weight = 0.0_my_real

      case (BALANCE_ELEMENTS)

         if (grid%element(elem)%isleaf .and. grid%element(elem)%iown) then
            grid%element(elem)%weight = 1.0_my_real
         else
            grid%element(elem)%weight = 0.0_my_real
         endif

      case (BALANCE_VERTICES, BALANCE_EQUATIONS)

         if (grid%element(elem)%isleaf .and. grid%element(elem)%iown) then
            grid%element(elem)%weight = assoc_verts(elem)
         else
            grid%element(elem)%weight = 0.0_my_real
         endif

      end select

! if predictive, scale by the error indicator which is only defined for leaves

      if (predictive .and. grid%element(elem)%isleaf) then
         select case (refcont%reftype)
         case (H_UNIFORM, H_ADAPTIVE, P_UNIFORM, P_ADAPTIVE)
            grid%element(elem)%weight = grid%element(elem)%weight * &
              maxval(grid%element_errind(elem,:))/grid%element(elem)%work
         case (HP_ADAPTIVE)
            select case (refcont%hp_strategy)
            case (HP_BIGGER_ERRIND)
               if (grid%element(elem)%level < refcont%max_lev) then
                  call error_indicator(grid,elem,LOCAL_PROBLEM_H, &
                                       energy=energy_errind,work=work)
                  err = maxval(energy_errind/work)
               else
                  err = 0.0_my_real
               endif
               if (grid%element(elem)%degree < refcont%max_deg) then
                  call error_indicator(grid,elem,LOCAL_PROBLEM_P, &
                                       energy=energy_errind,work=work)
                  err2 = maxval(energy_errind/work)
               else
                  err2 = 0.0_my_real
               endif
               grid%element(elem)%weight = grid%element(elem)%weight * &
                                           max(err,err2)
            case (HP_APRIORI, HP_PRIOR2P_E, HP_PRIOR2P_H1, HP_T3S, HP_NLP, &
                  HP_TYPEPARAM, HP_ALTERNATE, HP_COEF_DECAY, HP_COEF_ROOT, &
                  HP_SMOOTH_PRED, HP_NEXT3P, HP_REFSOLN_EDGE, HP_REFSOLN_ELEM)
               grid%element(elem)%weight = grid%element(elem)%weight * &
                 maxval(grid%element_errind(elem,:))/grid%element(elem)%work
            case default
               ierr = PHAML_INTERNAL_ERROR
               call fatal("bad case for hp_strategy in setting partitioning weights")
               stop
            end select
         end select
      endif

      elem = grid%element(elem)%next
   end do
end do

if (balance_what == BALANCE_VERTICES) then
   deallocate(assoc_verts)
endif

end subroutine set_weights

!          -------------
subroutine p_refine_elem(grid,elem,refine_control,elist,reftype,numhref, &
                         numpref,return_to_elist)
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
character(len=*), optional, intent(inout) :: reftype(:)
integer, optional, intent(inout) :: numhref(:),numpref(:)
logical, optional, intent(in) :: return_to_elist
!----------------------------------------------------
! Local variables:

integer :: i, j, deg, side, edge, edge_size, elem_size, oldsize, astat, d1,&
           d2, d3, old_edge_deg(EDGES_PER_ELEMENT), new_deg, &
           neigh(EDGES_PER_ELEMENT)
real(my_real), pointer :: temp1(:,:,:)
logical :: add_to_list, need_new_errind
!----------------------------------------------------
! Begin executable code

! if the error indicator list is provided, remove the element from it

if (present(elist)) then
   call remove_from_errind_list(elem,elist)
endif

! verify the element is a leaf

if (.not. grid%element(elem)%isleaf) return

! new degree for the element

deg = grid%element(elem)%degree + 1

! if the maximum degree would be exceeded, don't p-refine.

if (deg > refine_control%max_deg) then
   return
endif

! make sure allocated memory is large enough.  If not, reallocate to
! degree+3 so we don't reallocate at every p refinement

elem_size = ((deg-1)*(deg-2))/2
if (associated(grid%element(elem)%solution)) then
   oldsize = size(grid%element(elem)%solution,dim=1)
else
   oldsize = 0
endif
if (oldsize < elem_size) then
   elem_size = ((deg+1)*deg)/2
   allocate(temp1(elem_size,grid%system_size,max(1,grid%num_eval)),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in p_refine_elem")
      stop
   endif
   temp1 = 0.0_my_real
   if (oldsize > 0) temp1(1:oldsize,:,:) = grid%element(elem)%solution
   deallocate(grid%element(elem)%solution, stat=astat)
   grid%element(elem)%solution => temp1
   if (grid%have_true) then
      nullify(temp1)
      allocate(temp1(elem_size,grid%system_size,max(1,grid%num_eval)),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in p_refine_elem")
         stop
      endif
      temp1 = 0.0_my_real
      if (oldsize > 0) temp1(1:oldsize,:,:) = grid%element(elem)%exact
      deallocate(grid%element(elem)%exact, stat=astat)
      grid%element(elem)%exact => temp1
   endif
   if (grid%oldsoln_exists) then
      if (associated(grid%element(elem)%oldsoln)) then
         nullify(temp1)
         allocate(temp1(elem_size,grid%system_size,max(1,grid%num_eval)),stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("allocation failed in p_refine_elem")
            stop
         endif
         temp1 = 0.0_my_real
         d1 = min(size(grid%element(elem)%oldsoln,dim=1),size(temp1,dim=1))
         d2 = min(size(grid%element(elem)%oldsoln,dim=2),size(temp1,dim=2))
         d3 = min(size(grid%element(elem)%oldsoln,dim=3),size(temp1,dim=3))
         temp1(1:d1,1:d2,1:d3) = grid%element(elem)%oldsoln(1:d1,1:d2,1:d3)
         deallocate(grid%element(elem)%oldsoln)
         grid%element(elem)%oldsoln => temp1
      endif
   endif
endif

do side=1,EDGES_PER_ELEMENT
   edge_size = deg-1
   edge = grid%element(elem)%edge(side)
   if (associated(grid%edge(edge)%solution)) then
      oldsize = size(grid%edge(edge)%solution,dim=1)
   else
      oldsize = 0
   endif
   if (oldsize < edge_size) then
      edge_size = deg+1
      allocate(temp1(edge_size,grid%system_size,max(1,grid%num_eval)), &
               stat=astat)
      if (astat /= 0) then
         call fatal("allocation failed in p_refine_elem")
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
            call fatal("allocation failed in p_refine_elem")
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
            call fatal("allocation failed in p_refine_elem")
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
end do

! Increment the element degree by one and apply the appropriate rule for
! the edge degrees.  Also set Dirichlet boundary conditions on Dirichlet
! edges that change degree and set solution to 0 for other new components.

grid%element(elem)%degree = deg
if (deg >= 3) then
   grid%dof = grid%dof + deg-2
   if (grid%element(elem)%iown) grid%dof_own = grid%dof_own + deg-2
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
      grid%dof = grid%dof + new_deg - old_edge_deg(side)
      if (grid%element(grid%edge(edge)%assoc_elem)%iown) then
         grid%dof_own = grid%dof_own + new_deg - old_edge_deg(side)
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

if (deg > 2) then
   grid%element(elem)%solution(1+((deg-2)*(deg-3))/2:((deg-1)*(deg-2))/2,:,:) = 0.0_my_real
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

! Set initial guess for new solution components.
! Don't set it for REFSOLN because it tries to evaluate a neighbor solution
! which hasn't been allocated yet.

if (refine_control%reftype /= HP_ADAPTIVE .or. &
    (refine_control%hp_strategy /= HP_REFSOLN_EDGE .and. &
     refine_control%hp_strategy /= HP_REFSOLN_ELEM)) then
   if (refine_control%error_estimator == INITIAL_CONDITION) then
      call init_guess_ic(grid,elem,.false.)
   else
      call init_guess_p(grid,elem,old_edge_deg)
   endif
endif

! predicted error for SMOOTH_PRED hp strategy
! TEMP only using first eigenvector

grid%element(elem)%sp_eta_pred = sqrt(refine_control%sp_gamma_p)*grid%element_errind(elem,1)

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

if (present(reftype)) then
   if (.not. present(numhref) .or. .not. present(numpref)) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("must have numhref and numpref in p_refine_elem if reftype is present")
      stop
   endif
   call mark_reftype_one(elem,grid,refine_control,1.0_my_real,reftype,numhref, &
                         numpref,.true.)
endif

end subroutine p_refine_elem

!                    --------------------
recursive subroutine bisect_triangle_pair(grid,parent,errcode,refcont,elist, &
                                          reftype,numhref,numpref, &
                                          return_to_elist,reduce_p,reduce_p_max)
!                    --------------------

!----------------------------------------------------
! This routine refines triangle parent and its mate by bisection.
! errcode = 0 for success, 1 if the grid is full
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

! if the error indicator list is provided, remove the parent from it

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

if (grid%element(parent)%level >= size(grid%head_level_elem)) then
   call extend_nlev(grid)
   if (ierr /= NO_ERROR) return
endif

! reduce_p and reduce_p_max should both be present or absent

if (present(reduce_p) .neqv. present(reduce_p_max)) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("reduce_p and reduce_p_max must both be present or absent in bisect_triangle_pair")
   stop
endif

! identify the mate, creating it if necessary

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

! if the mate was scheduled for p refinement, perform that p refinement first

if (mate /= BOUNDARY .and. present(reftype)) then
   if (reftype(mate)=="p") then
      call p_refine_elem(grid,mate,refcont,elist,reftype,numhref, &
                         numpref,return_to_elist)
   endif
endif

! get the indices for the new elements, edges and vertex

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
grid%element(child1)%next = grid%head_level_elem(grid%element(child1)%level)
if (grid%element(child1)%next /= END_OF_LIST) then
   grid%element(grid%element(child1)%next)%previous = child1
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
   grid%dof = grid%dof + &
           ((grid%element(child1)%degree-2)*(grid%element(child1)%degree-1))/2
endif
if (grid%element(child1)%iown .and. grid%element(child1)%degree >= 3) then
   grid%dof_own = grid%dof_own + &
           ((grid%element(child1)%degree-2)*(grid%element(child1)%degree-1))/2
endif

call hash_insert(grid%element(child1)%gid,child1,grid%elem_hash)

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
grid%element(child2)%next = child1
grid%element(child2)%hrefined_unowned = .false.
grid%element(child2)%prefined_unowned = .false.
grid%element(child1)%previous = child2
grid%element(child2)%previous = END_OF_LIST
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

call hash_insert(grid%element(child2)%gid,child2,grid%elem_hash)

if (.not. grid%element(parent)%iown) then
   grid%element(parent)%hrefined_unowned = .true.
endif

if (mate == BOUNDARY) then

grid%head_level_elem(grid%element(child2)%level) = child2

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
grid%element(child3)%next = child2
grid%element(child2)%previous = child3
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
   grid%dof = grid%dof + &
           ((grid%element(child3)%degree-2)*(grid%element(child3)%degree-1))/2
endif
if (grid%element(child3)%iown .and. grid%element(child3)%degree >= 3) then
   grid%dof_own = grid%dof_own + &
           ((grid%element(child3)%degree-2)*(grid%element(child3)%degree-1))/2
endif

call hash_insert(grid%element(child3)%gid,child3,grid%elem_hash)

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
grid%element(child4)%next = child3
grid%element(child3)%previous = child4
grid%element(child4)%previous = END_OF_LIST
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

call hash_insert(grid%element(child4)%gid,child4,grid%elem_hash)

grid%head_level_elem(grid%element(child4)%level) = child4

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
grid%vertex(vert1)%next = grid%head_level_vert(grid%element(child1)%level)
grid%head_level_vert(grid%element(child1)%level) = vert1
if (grid%vertex(vert1)%next /= END_OF_LIST) then
   grid%vertex(grid%vertex(vert1)%next)%previous = vert1
endif
grid%vertex(vert1)%previous = END_OF_LIST
if (grid%have_true) then
   do j=1,max(1,grid%num_eval)
      do i=1,grid%system_size
         grid%vertex_exact(vert1,i,j)=trues(grid%vertex(vert1)%coord%x, &
                                            grid%vertex(vert1)%coord%y, &
                                            i,j)
      end do
   end do
endif
grid%dof = grid%dof + 1
if (grid%element(grid%vertex(vert1)%assoc_elem)%iown) grid%dof_own = grid%dof_own + 1

call hash_insert(grid%vertex(vert1)%gid,vert1,grid%vert_hash)

! second new vertex, if periodic boundary point(s)

masterparent = -1

if (any(grid%vertex_type(vert1,:) == PERIODIC_MASTER) .or. &
    any(grid%vertex_type(vert1,:) == PERIODIC_SLAVE)) then

   if (grid%next_free_vert == END_OF_LIST) then
      call more_verts(grid,errcode)
      if (errcode /= 0) return
   endif
   vert2 = grid%next_free_vert
   grid%next_free_vert = grid%vertex(vert2)%next

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
      call hash_remove(grid%vertex(vert1)%gid,grid%vert_hash)
      tempgid = grid%vertex(vert1)%gid
      grid%vertex(vert1)%gid = grid%vertex(vert2)%gid
      grid%vertex(vert2)%gid = tempgid
      call hash_insert(grid%vertex(vert1)%gid,vert1,grid%vert_hash)
      if (grid%vertex(vert1)%gid == grid%element(child2)%gid) then
         grid%vertex(vert1)%assoc_elem = child2
      elseif (grid%vertex(vert1)%gid == grid%element(child4)%gid) then
         grid%vertex(vert1)%assoc_elem = child4
      endif
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
   if (grid%have_true) then
      do j=1,max(1,grid%num_eval)
         do i=1,grid%system_size
            grid%vertex_exact(vert2,i,j)=trues(grid%vertex(vert2)%coord%x, &
                                               grid%vertex(vert2)%coord%y, &
                                               i,j)
         end do
      end do
   endif

   call hash_insert(grid%vertex(vert2)%gid,vert2,grid%vert_hash)
   
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
   grid%dof = grid%dof + grid%edge(edge1)%degree - 1
   if (grid%element(grid%edge(edge1)%assoc_elem)%iown) then
      grid%dof_own = grid%dof_own + grid%edge(edge1)%degree - 1
   endif
endif

call hash_insert(grid%edge(edge1)%gid,edge1,grid%edge_hash)

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
call hash_insert(grid%edge(edge2)%gid,edge2,grid%edge_hash)

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
   grid%dof = grid%dof + grid%edge(edge3)%degree - 1
   if (grid%element(grid%edge(edge3)%assoc_elem)%iown) then
      grid%dof_own = grid%dof_own + grid%edge(edge3)%degree - 1
   endif
endif

call hash_insert(grid%edge(edge3)%gid,edge3,grid%edge_hash)

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
      grid%dof = grid%dof + grid%edge(edge4)%degree - 1
      if (grid%element(grid%edge(edge4)%assoc_elem)%iown) then
         grid%dof_own = grid%dof_own + grid%edge(edge4)%degree - 1
      endif
   endif
   call hash_insert(grid%edge(edge4)%gid,edge4,grid%edge_hash)
endif

! 5th and 6th edges if periodic boundary edge

if (any(grid%edge_type(edge1,:) == PERIODIC_MASTER) .or. &
    any(grid%edge_type(edge1,:) == PERIODIC_SLAVE)) then

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

   grid%element(child3)%edge(2) = edge5
   grid%element(child4)%edge(2) = edge6

! edge 5

   grid%edge(edge5)%gid = 4*grid%edge(mateedge)%gid
   grid%edge(edge5)%vertex = (/ grid%edge(mateedge)%vertex(1),vert2 /)
   grid%edge(edge5)%bmark = grid%edge(mateedge)%bmark
   grid%edge(edge5)%degree = grid%edge(mateedge)%degree
   grid%edge_type(edge5,:) = grid%edge_type(mateedge,:)
   grid%edge(edge5)%assoc_elem = grid%edge(edge1)%assoc_elem
   if (any(grid%edge_type(edge5,:) == PERIODIC_SLAVE)) then
      grid%edge(edge5)%next = edge1
   else
      grid%edge(edge1)%next = edge5
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
   call hash_insert(grid%edge(edge5)%gid,edge5,grid%edge_hash)

! edge 6

   grid%edge(edge6)%gid = 4*grid%edge(mateedge)%gid + 1
   grid%edge(edge6)%vertex = (/ grid%edge(mateedge)%vertex(2),vert2 /)
   grid%edge(edge6)%bmark = grid%edge(mateedge)%bmark
   grid%edge(edge6)%degree = grid%edge(mateedge)%degree
   grid%edge_type(edge6,:) = grid%edge_type(mateedge,:)
   grid%edge(edge6)%assoc_elem = grid%edge(edge2)%assoc_elem
   if (any(grid%edge_type(edge6,:) == PERIODIC_SLAVE)) then
      grid%edge(edge6)%next = edge2
   else
      grid%edge(edge2)%next = edge6
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
   call hash_insert(grid%edge(edge6)%gid,edge6,grid%edge_hash)

endif

! if reduce_p is present, the children generated for compatibility should have
! their degree reduced by 1/sqrt(2).  This is for the REFSOLN_ELEM hp strategy.
! reduce_p is passed in as 0 and incremented with each recursive call.  Thus
! reduce_p==0 indicates the original call, so the reduction does not occur.
! reduce_p_max holds the largest value of reduce_p in this recursive chain.
! If reduce_p is not reduce_p_max, then the children of the mate are not
! reduced, because they were reduced when the mate's parent was refined.

if (present(reduce_p)) then
   if (reduce_p > 0) then
      newdeg = floor((grid%element(child1)%degree+1)/sqrt(2.0_my_real))
      do while (grid%element(child1)%degree > newdeg)
         call p_coarsen_elem(grid,child1,errcode,refcont)
      end do
      do while (grid%element(child2)%degree > newdeg)
         call p_coarsen_elem(grid,child2,errcode,refcont)
      end do
      if ((.not. mate == BOUNDARY) .and. reduce_p == reduce_p_max) then
         newdeg = floor((grid%element(child3)%degree+1)/sqrt(2.0_my_real))
         do while (grid%element(child3)%degree > newdeg)
            call p_coarsen_elem(grid,child3,errcode,refcont)
         end do
         do while (grid%element(child4)%degree > newdeg)
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
if (grid%element(child1)%level > grid%nlev) then
   grid%nlev = grid%element(child1)%level
endif
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

! set initial guess for solution

if (refcont%error_estimator == INITIAL_CONDITION) then
   call init_guess_ic(grid,parent,.true.)
else
   call init_guess_h(grid,parent)
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
      call error_indicator(grid,child1,refcont%error_estimator, &
                        grid%element_errind(child1,:),grid%element(child1)%work)
      call error_indicator(grid,child2,refcont%error_estimator, &
                        grid%element_errind(child2,:),grid%element(child2)%work)
   endif
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
         call error_indicator(grid,child3,refcont%error_estimator, &
                              grid%element_errind(child3,:), &
                              grid%element(child3)%work)
         call error_indicator(grid,child4,refcont%error_estimator, &
                              grid%element_errind(child4,:), &
                              grid%element(child4)%work)
      endif
      if (add_to_list) then
         call add_to_errind_list(child3,elist,grid,refcont)
         call add_to_errind_list(child4,elist,grid,refcont)
      endif
   endif

endif

if (present(reftype)) then
   if (.not. present(numhref) .or. .not. present(numpref)) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("must have numhref and numpref in bisect_triangle_pair if reftype is present")
      stop
   endif
   reftype(parent) = "n"
   call mark_reftype_one(child1,grid,refcont,1.0_my_real,reftype,numhref, &
                         numpref,.true.)
   call mark_reftype_one(child2,grid,refcont,1.0_my_real,reftype,numhref, &
                         numpref,.true.)
   if (mate /= BOUNDARY) then
      reftype(mate) = "n"
      call mark_reftype_one(child3,grid,refcont,1.0_my_real,reftype,numhref, &
                         numpref,.true.)
      call mark_reftype_one(child4,grid,refcont,1.0_my_real,reftype,numhref, &
                         numpref,.true.)
   endif
endif

end subroutine bisect_triangle_pair

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

!        -----------
function level2_mate(grid,parent,k)
!        -----------

!----------------------------------------------------
! This routine determines the mate of a triangle on level 2, which
! should be a child of the k'th neighbor of the level 1 parent
! RESTRICTION triangles
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: parent,k
type(hash_key) :: level2_mate
!----------------------------------------------------
! Local variables:

integer :: neigh, i
!----------------------------------------------------
! Begin executable code

neigh = grid%initial_neighbor(k,parent)

! easy case -- boundary

if (neigh == BOUNDARY) then
   level2_mate = BOUNDARY
   return
endif

! determine which neighbor parent is of the neighbor

do i=1,3
   if (grid%initial_neighbor(i,neigh) == parent) exit
end do
if (i > 3) then
   call fatal("initial_neighbor is not reflexive")
   stop
endif
if (i == 3) then
   call fatal("initial_neighbor mates are not reflexive")
   stop
endif

! compute the mate

level2_mate = 2*grid%element(neigh)%gid + (2-i)

end function level2_mate

!          ------------
subroutine init_guess_h(grid,elem)
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
integer, intent(in) :: elem
!----------------------------------------------------
! Local variables:

real(my_real) :: xvert(3), yvert(3), rmin
real(my_real), allocatable :: elem_mat(:,:),elem_rs(:,:),full_mat(:,:), &
                              full_rs(:,:)
integer :: mate, child(4), degree(4), nbasis_elem, nbasis_full, elem_n, &
           full_n, astat, edge_type(3,grid%system_size), ss, i, j, k, info, &
           bmark(3), start_renum(17), nb, nev
integer, allocatable :: renum(:), ipiv(:)
!----------------------------------------------------
! Begin executable code

ss = grid%system_size
nev = max(1,grid%num_eval)

! lid of the mate

if (grid%element(elem)%mate == BOUNDARY) then
   mate = BOUNDARY
else
   mate = hash_decode_key(grid%element(elem)%mate,grid%elem_hash)
endif

! lids of the children

child(1) = hash_decode_key(2*grid%element(elem)%gid,grid%elem_hash)
child(2) = hash_decode_key(2*grid%element(elem)%gid+1,grid%elem_hash)
if (mate == BOUNDARY) then
   child(3) = BOUNDARY
   child(4) = BOUNDARY
else
   child(3) = hash_decode_key(2*grid%element(mate)%gid,grid%elem_hash)
   child(4) = hash_decode_key(2*grid%element(mate)%gid+1,grid%elem_hash)
endif
if (any(child == HASH_NOT_FOUND)) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("didn't find child in init_guess_h")
   stop
endif

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
                            0, .false., elem, rmin, "p", "a", elem_mat, &
                          elem_rs(:,1), loc_bconds_s=init_guess_dirich_bconds, &
                            extra_rhs=elem_rs(:,2:))
   else
      call elemental_matrix(grid, xvert, yvert, degree, edge_type, bmark, ss, &
                            0, .false., elem, rmin, "p", "a", elem_mat, &
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
                            0, .false., elem, rmin, "p", "a", elem_mat, &
                          elem_rs(:,1), loc_bconds_s=init_guess_dirich_bconds, &
                            extra_rhs=elem_rs(:,2:))
   else
      call elemental_matrix(grid, xvert, yvert, degree, edge_type, bmark, ss, &
                            0, .false., elem, rmin, "p", "a", elem_mat, &
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

   grid%vertex_solution(grid%element(child(1))%vertex(1),j,k) = &
      full_rs(j+ss*(start_renum(1)-1),k)
   grid%vertex_solution(grid%element(child(2))%vertex(1),j,k) = &
      full_rs(j+ss*(start_renum(2)-1),k)
   grid%vertex_solution(grid%element(child(1))%vertex(2),j,k) = &
      full_rs(j+ss*(start_renum(3)-1),k)
   if (mate /= BOUNDARY) then
      grid%vertex_solution(grid%element(child(3))%vertex(2),j,k) = &
         full_rs(j+ss*(start_renum(4)-1),k)
   endif
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
   do i=1,grid%edge(grid%element(child(1))%edge(3))%degree-1
      grid%edge(grid%element(child(1))%edge(3))%solution(i,j,k) = &
         full_rs(j+ss*(start_renum(10)+i-2),k)
   end do
   do i=1,grid%edge(grid%element(child(2))%edge(3))%degree-1
      grid%edge(grid%element(child(2))%edge(3))%solution(i,j,k) = &
         full_rs(j+ss*(start_renum(11)+i-2),k)
   end do
   if (mate /= BOUNDARY) then
      do i=1,grid%edge(grid%element(child(3))%edge(3))%degree-1
         grid%edge(grid%element(child(3))%edge(3))%solution(i,j,k)= &
            full_rs(j+ss*(start_renum(12)+i-2),k)
      end do
      do i=1,grid%edge(grid%element(child(4))%edge(3))%degree-1
         grid%edge(grid%element(child(4))%edge(3))%solution(i,j,k)= &
            full_rs(j+ss*(start_renum(13)+i-2),k)
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
nred = degree(1+EDGES_PER_ELEMENT) - 2
do i=1,EDGES_PER_ELEMENT
   if (old_edge_deg(i) >= degree(1+EDGES_PER_ELEMENT)) then
      degree(i) = 1 ! forces no red bases on this side
   else
      degree(i) = degree(1+EDGES_PER_ELEMENT)
      nred = nred + 1
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

!          -------------
subroutine init_guess_ic(grid,elem,isparent)
!          -------------

!----------------------------------------------------
! This routine sets the solution in either element elem or the descendants
! of elem and its mate from the function in iconds.  isparent determines
! if elem is the element or the parent.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: elem
logical, intent(in) :: isparent
!----------------------------------------------------
! Local variables:

integer :: vert(5), edge(8), face(4)
integer :: nvert, nedge, nface, v, e, f, i, comp, eigen
!----------------------------------------------------
! Begin executable code

! set the entities to do

if (isparent) then
   face(1:2) = get_child_lid(grid%element(elem)%gid,(/1,2/),grid%elem_hash)
   edge(1:3) = grid%element(face(1))%edge
   edge(4:5) = grid%element(face(2))%edge(2:3)
   vert(1:3) = grid%element(elem)%vertex
   vert(4)   = grid%element(face(1))%vertex(3)
   if (grid%element(elem)%mate == BOUNDARY) then
      nvert = 4
      nedge = 5
      nface = 2
   else
      face(3:4) = get_child_lid(grid%element(elem)%mate,(/1,2/), &
                  grid%elem_hash)
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

!          --------------
subroutine compute_normal(grid,x,y,elem,edge,normal)
!          --------------

!----------------------------------------------------
! This routine calculates the outward unit normal for the edge of
! element elem containing point (x,y).  The index of the edge in
! elem is also returned.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
real(my_real), intent(in) :: x,y
integer, intent(in) :: elem
integer, intent(out) :: edge
real(my_real), intent(out) :: normal(2)
!----------------------------------------------------
! Local variables:

integer :: iarr(1)
real(my_real) :: x1, x2, x3, y1, y2, y3, zeta(3)
logical :: clockwise
!----------------------------------------------------
! Begin executable code

! determine which edge of the triangle this point is on by computing the
! barycentric coordinates and looking for one that is 0

zeta = barycentric(x,y,elem,grid,no_det=.true.)

iarr = minloc(abs(zeta))
edge = iarr(1)

x1 = grid%vertex(grid%element(elem)%vertex(1))%coord%x
x2 = grid%vertex(grid%element(elem)%vertex(2))%coord%x
x3 = grid%vertex(grid%element(elem)%vertex(3))%coord%x
y1 = grid%vertex(grid%element(elem)%vertex(1))%coord%y
y2 = grid%vertex(grid%element(elem)%vertex(2))%coord%y
y3 = grid%vertex(grid%element(elem)%vertex(3))%coord%y

! determine if the vertices are clockwise by checking the sign of the cross
! product

clockwise = (x2-x1)*(y3-y2)-(y2-y1)*(x3-x2) < 0.0_my_real

! change the use of x1 etc to be the two endpoints of the edge

if (edge == 1) then
   x1 = x2
   y1 = y2
   x2 = x3
   y2 = y3
elseif (edge == 2) then
   x2 = x1
   y2 = y1
   x1 = x3
   y1 = y3
endif

! compute the outward unit normal vector for this edge

if (clockwise) then
   normal(1) = -(y2-y1)/sqrt((x2-x1)**2+(y2-y1)**2)
   normal(2) =  (x2-x1)/sqrt((x2-x1)**2+(y2-y1)**2)
else
   normal(1) =  (y2-y1)/sqrt((x2-x1)**2+(y2-y1)**2)
   normal(2) = -(x2-x1)/sqrt((x2-x1)**2+(y2-y1)**2)
endif

end subroutine compute_normal

!        -----------
function barycentric(x,y,elem,grid,no_det) result(zeta)
!        -----------

!----------------------------------------------------
! This function returns the barycentric coordinates of (x,y) in element elem.
! If no_det is present and .true., the coordinates are not divided by the
! determinant (area of the triangle).  This is useful if you are only looking
! at the signs of the coordinates, where it can be wrong if the point is
! very close to an edge and the area is very small so roundoff can be very bad.
! RESTRICTION triangles
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x,y
integer, intent(in) :: elem
type(grid_type), intent(in) :: grid
real(my_real), dimension(VERTICES_PER_ELEMENT) :: zeta
logical, optional, intent(in) :: no_det

!----------------------------------------------------
! Local variables:

real(quad_real) :: x1,x2,x3,y1,y2,y3,det
real(quad_real) :: xy1,xy2,xy3,yx1,yx2,yx3,x1y2,x2y3,x3y1,x1y3,x2y1,x3y2
real(quad_real) :: sump,summ
logical :: local_no_det

!----------------------------------------------------
! Begin executable code

! check for request for no determinant

if (present(no_det)) then
   local_no_det = no_det
else
   local_no_det = .false.
endif

! local variables for the vertices, to make the code easier to read

x1 = real(grid%vertex(grid%element(elem)%vertex(1))%coord%x,quad_real)
x2 = real(grid%vertex(grid%element(elem)%vertex(2))%coord%x,quad_real)
x3 = real(grid%vertex(grid%element(elem)%vertex(3))%coord%x,quad_real)
y1 = real(grid%vertex(grid%element(elem)%vertex(1))%coord%y,quad_real)
y2 = real(grid%vertex(grid%element(elem)%vertex(2))%coord%y,quad_real)
y3 = real(grid%vertex(grid%element(elem)%vertex(3))%coord%y,quad_real)

! compute the barycentric coordinates of the point

! reduce roundoff by summing all the positive parts and negative parts
! separately, and then adding the two partial sums

! products needed for the sums

xy1 = real(x,quad_real)*y1
xy2 = real(x,quad_real)*y2
xy3 = real(x,quad_real)*y3
yx1 = real(y,quad_real)*x1
yx2 = real(y,quad_real)*x2
yx3 = real(y,quad_real)*x3
x1y2 = x1*y2
x2y3 = x2*y3
x3y1 = x3*y1
x1y3 = x1*y3
x2y1 = x2*y1
x3y2 = x3*y2

! compute the determinant

! det = x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2)

sump=0.0_quad_real; summ=0.0_quad_real
if (x1y2 > 0.0_quad_real) then; sump=sump+x1y2; else; summ=summ+x1y2; endif
if (x1y3 < 0.0_quad_real) then; sump=sump-x1y3; else; summ=summ-x1y3; endif
if (x2y3 > 0.0_quad_real) then; sump=sump+x2y3; else; summ=summ+x2y3; endif
if (x2y1 < 0.0_quad_real) then; sump=sump-x2y1; else; summ=summ-x2y1; endif
if (x3y1 > 0.0_quad_real) then; sump=sump+x3y1; else; summ=summ+x3y1; endif
if (x3y2 < 0.0_quad_real) then; sump=sump-x3y2; else; summ=summ-x3y2; endif
det = sump + summ

! if the request is to not divide by the determinant, only use its sign

if (local_no_det) then
   det = sign(1.0_quad_real,det)
endif

! compute the coordinates

! zeta(1) = (x*(y2-y3) + y*(x3-x2) + (x2*y3-x3*y2))/det

sump=0.0_quad_real; summ=0.0_quad_real
if (xy2 > 0.0_quad_real) then; sump=sump+xy2; else; summ=summ+xy2; endif
if (xy3 < 0.0_quad_real) then; sump=sump-xy3; else; summ=summ-xy3; endif
if (yx3 > 0.0_quad_real) then; sump=sump+yx3; else; summ=summ+yx3; endif
if (yx2 < 0.0_quad_real) then; sump=sump-yx2; else; summ=summ-yx2; endif
if (x2y3 > 0.0_quad_real) then; sump=sump+x2y3; else; summ=summ+x2y3; endif
if (x3y2 < 0.0_quad_real) then; sump=sump-x3y2; else; summ=summ-x3y2; endif
zeta(1) = sump/det + summ/det

! zeta(2) = (x*(y3-y1) + y*(x1-x3) + (x3*y1-x1*y3))/det

sump=0.0_quad_real; summ=0.0_quad_real
if (xy3 > 0.0_quad_real) then; sump=sump+xy3; else; summ=summ+xy3; endif
if (xy1 < 0.0_quad_real) then; sump=sump-xy1; else; summ=summ-xy1; endif
if (yx1 > 0.0_quad_real) then; sump=sump+yx1; else; summ=summ+yx1; endif
if (yx3 < 0.0_quad_real) then; sump=sump-yx3; else; summ=summ-yx3; endif
if (x3y1 > 0.0_quad_real) then; sump=sump+x3y1; else; summ=summ+x3y1; endif
if (x1y3 < 0.0_quad_real) then; sump=sump-x1y3; else; summ=summ-x1y3; endif
zeta(2) = sump/det + summ/det

! zeta(3) = (x*(y1-y2) + y*(x2-x1) + (x1*y2-x2*y1))/det

sump=0.0_quad_real; summ=0.0_quad_real
if (xy1 > 0.0_quad_real) then; sump=sump+xy1; else; summ=summ+xy1; endif
if (xy2 < 0.0_quad_real) then; sump=sump-xy2; else; summ=summ-xy2; endif
if (yx2 > 0.0_quad_real) then; sump=sump+yx2; else; summ=summ+yx2; endif
if (yx1 < 0.0_quad_real) then; sump=sump-yx1; else; summ=summ-yx1; endif
if (x1y2 > 0.0_quad_real) then; sump=sump+x1y2; else; summ=summ+x1y2; endif
if (x2y1 < 0.0_quad_real) then; sump=sump-x2y1; else; summ=summ-x2y1; endif
zeta(3) = sump/det + summ/det

if (VERTICES_PER_ELEMENT /= 3) then
   zeta(4:VERTICES_PER_ELEMENT) = 0.0_my_real
   call warning("function barycentric needs to be changed for nontriangles")
endif

end function barycentric

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

if (refine_control%reftype == H_UNIFORM .or. &
    refine_control%reftype == P_UNIFORM) then
   return
endif

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

call phaml_alltoall(procs,isend,nsend,irecv,nrecv,450)

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

call phaml_alltoall(procs,isend,nsend,irecv,nrecv,451)

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
   grid%dof = grid%dof - deg + 1
   if (grid%element(elem)%iown) grid%dof_own = grid%dof_own - deg + 1
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
      grid%dof = grid%dof + new_deg - old_edge_deg(side)
      if (grid%element(grid%edge(edge)%assoc_elem)%iown) then
         grid%dof_own = grid%dof_own + new_deg - old_edge_deg(side)
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
grid%dof = grid%dof - 1
if (grid%element(grid%vertex(vert1)%assoc_elem)%iown) then
   grid%dof_own = grid%dof_own - 1
endif
do i=1,4
   if (i>2 .and. mate==BOUNDARY) exit
   if (grid%element(child(i))%degree >= 3) then
      grid%dof = grid%dof - &
        ((grid%element(child(i))%degree-1)*(grid%element(child(i))%degree-2))/2
      if (grid%element(child(i))%iown) then
       grid%dof_own = grid%dof_own - &
        ((grid%element(child(i))%degree-1)*(grid%element(child(i))%degree-2))/2
      endif
   endif
end do
if (par_deg >= 3) then
   grid%dof = grid%dof + ((par_deg-1)*(par_deg-2))/2
   if (grid%element(parent)%iown) then
      grid%dof_own = grid%dof_own + ((par_deg-1)*(par_deg-2))/2
   endif
endif
if (mate /= BOUNDARY .and. mate_deg >= 3) then
   grid%dof = grid%dof + ((mate_deg-1)*(mate_deg-2))/2
   if (grid%element(mate)%iown) then
      grid%dof_own = grid%dof_own + ((mate_deg-1)*(mate_deg-2))/2
   endif
endif
do i=1,4
   if (i>3 .and. mate==BOUNDARY) exit
   grid%dof = grid%dof - grid%edge(childedge(i))%degree + 1
   if (grid%element(grid%edge(childedge(i))%assoc_elem)%iown) then
      grid%dof_own = grid%dof_own - grid%edge(childedge(i))%degree + 1
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
      grid%dof = grid%dof - old_edge_deg(i) + grid%edge(edge)%degree
      if (grid%element(grid%edge(edge)%assoc_elem)%iown) then
         grid%dof_own = grid%dof_own - old_edge_deg(i) + grid%edge(edge)%degree
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

!          -----------
subroutine phier2nodal(solution,xvert,yvert,degree)
!          -----------

!----------------------------------------------------
! This routine converts the values in solution from coefficients of a
! p-hierarchical basis to coefficients of a nodal basis, over four triangles
! that are siblings in an h refinement with the outer vertices given in
! (xvert,yvert).  The order of the coefficients (using cubics as an example) is

!                  vertex 1
!                      1
!                    / | \
!                   10 8  20
!                  / 12|22 \
!                 11   9    21
!                /     |      \
!    vertex 3   3-5-6--7-19-18-4  vertex 4
!                \     14     /
!                 16   |    24
!                  \ 17|25 /
!                   15 13 23
!                    \ | /
!                      2
!                  vertex 2

! For the p-hierarchical basis, those on the same edge or face functions of
! the same element are in order of degree, as in the basis function routines.
! All elements and edges must have the same degree.  The third and fourth
! triangles can be omitted by setting the 4th vertex to be huge(0.0_my_real).

! If the third and fourth triangles are omitted, the indices remain the same.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(inout) :: solution(:,:)
real(my_real), intent(in) :: xvert(:),yvert(:)
integer, intent(in) :: degree
!----------------------------------------------------
! Local variables:

real(my_real) :: xnode(size(solution,dim=1)),ynode(size(solution,dim=1)), &
                 basis(((degree+1)*(degree+2))/2,size(solution,dim=1)), &
                 new_soln(size(solution,dim=1),size(solution,dim=2))
real(my_real) :: rdeg,xmid,ymid,fact,fact1,fact2,fact3
integer :: i,j,isub,nsub,nbasis,last
integer :: jsub(((degree+1)*(degree+2))/2), ksub(((degree+1)*(degree+2))/2)
!----------------------------------------------------
! Begin executable code


rdeg = real(degree,my_real)
nbasis = ((degree+1)*(degree+2))/2

! midpoint of the edge between vertices 1 and 2

xmid = (xvert(1)+xvert(2))/2
ymid = (yvert(1)+yvert(2))/2

! the coefficients of the vertex bases (1-4 and the midpoint 4+degree) are
! the same in both bases

new_soln = 0
new_soln(1:4,:) = solution(1:4,:)
new_soln(4+degree,:) = solution(4+degree,:)

! define the coordinates of the as yet unaddressed nodes in the first triangle,
! and set jsub to the index corresponding to those nodes,
! and set ksub to the index of the basis functions

ksub(1) = 1
ksub(2) = 3
ksub(3) = 4+degree
isub = 1
do i=1,degree-1
   fact = i/rdeg
   xnode(isub) = fact*xmid + (1-fact)*xvert(3)
   ynode(isub) = fact*ymid + (1-fact)*yvert(3)
   jsub(isub) = 4+i
   ksub(isub+3) = jsub(isub)
   isub = isub+1
end do
do i=1,degree-1
   fact = i/rdeg
   xnode(isub) = fact*xmid + (1-fact)*xvert(1)
   ynode(isub) = fact*ymid + (1-fact)*yvert(1)
   jsub(isub) = 4+degree+i
   ksub(isub+3) = jsub(isub)
   isub = isub+1
end do
do i=1,degree-1
   fact = i/rdeg
   xnode(isub) = fact*xvert(3) + (1-fact)*xvert(1)
   ynode(isub) = fact*yvert(3) + (1-fact)*yvert(1)
   jsub(isub) = 3+2*degree+i
   ksub(isub+3) = jsub(isub)
   isub = isub+1
end do
do i=1,degree-2
   fact3 = i/rdeg
   do j=1,degree-1-i
      fact2 = (1-fact3)*j/(rdeg-1)
      fact1 = 1-fact3-fact2
      xnode(isub) = fact1*xvert(1) + fact2*xvert(3) + fact3*xmid
      ynode(isub) = fact1*yvert(1) + fact2*yvert(3) + fact3*ymid
      jsub(isub) = 5+isub
      ksub(isub+3) = jsub(isub)
      isub = isub+1
   end do
end do

nsub = 3*(degree-1) + ((degree-2)*(degree-1))/2
if (nsub == 0) then
   last = 0
else
   last = jsub(nsub)
endif

! evaluate the p-hierarchical basis functions at the nodes

if (nsub /= 0) then
   call p_hier_basis_func(xnode(1:nsub),ynode(1:nsub), &
                          (/xvert(1),xvert(3),xmid/), &
                          (/yvert(1),yvert(3),ymid/), &
                          (/degree,degree,degree,degree/),"a",basis(:,1:nsub))
endif

! evaluate the solution at the nodes

do i=1,nsub
   new_soln(jsub(i),:) = 0
   do j=1,nbasis
      new_soln(jsub(i),:) = new_soln(jsub(i),:) + solution(ksub(j),:)*basis(j,i)
   end do
end do

! same process for the second triangle

ksub(1) = 2
ksub(2) = 3
ksub(3) = 4+degree
do i=1,degree-1
   ksub(i+3) = i+4
end do
isub = 1
do i=1,degree-1
   fact = i/rdeg
   xnode(isub) = fact*xmid + (1-fact)*xvert(2)
   ynode(isub) = fact*ymid + (1-fact)*yvert(2)
   jsub(isub) = last + isub
   ksub(isub+degree+2) = jsub(isub)
   isub = isub+1
end do
do i=1,degree-1
   fact = i/rdeg
   xnode(isub) = fact*xvert(3) + (1-fact)*xvert(2)
   ynode(isub) = fact*yvert(3) + (1-fact)*yvert(2)
   jsub(isub) = last + isub
   ksub(isub+degree+2) = jsub(isub)
   isub = isub+1
end do
do i=1,degree-2
   fact3 = i/rdeg
   do j=1,degree-1-i
      fact2 = (1-fact3)*j/(rdeg-1)
      fact1 = 1-fact3-fact2
      xnode(isub) = fact1*xvert(2) + fact2*xvert(3) + fact3*xmid
      ynode(isub) = fact1*yvert(2) + fact2*yvert(3) + fact3*ymid
      jsub(isub) = last + isub
      ksub(isub+degree+2) = jsub(isub)
      isub = isub+1
   end do
end do
nsub = 2*(degree-1) + ((degree-2)*(degree-1))/2
if (nsub == 0) then
   last = 0
else
   last = jsub(nsub)
endif
if (nsub /= 0) then
   call p_hier_basis_func(xnode(1:nsub),ynode(1:nsub), &
                          (/xvert(2),xvert(3),xmid/), &
                          (/yvert(2),yvert(3),ymid/), &
                          (/degree,degree,degree,degree/),"a",basis(:,1:nsub))
endif
do i=1,nsub
   new_soln(jsub(i),:) = 0
   do j=1,nbasis
      new_soln(jsub(i),:) = new_soln(jsub(i),:) + solution(ksub(j),:)*basis(j,i)
   end do
end do

! same process for the third triangle if it exists

if (xvert(4) /= huge(0.0_my_real)) then
   ksub(1) = 1
   ksub(2) = 4
   ksub(3) = 4+degree
   isub = 1
   do i=1,degree-1
      fact = i/rdeg
      xnode(isub) = fact*xmid + (1-fact)*xvert(4)
      ynode(isub) = fact*ymid + (1-fact)*yvert(4)
      jsub(isub) = last + isub
      ksub(isub+3) = jsub(isub)
      isub = isub+1
   end do
   do i=1,degree-1
      ksub(isub+2+i) = 4+degree+i
   end do
   do i=1,degree-1
      fact = i/rdeg
      xnode(isub) = fact*xvert(4) + (1-fact)*xvert(1)
      ynode(isub) = fact*yvert(4) + (1-fact)*yvert(1)
      jsub(isub) = last + isub
      ksub(isub+degree+2) = jsub(isub)
      isub = isub+1
   end do
   do i=1,degree-2
      fact3 = i/rdeg
      do j=1,degree-1-i
         fact2 = (1-fact3)*j/(rdeg-1)
         fact1 = 1-fact3-fact2
         xnode(isub) = fact1*xvert(1) + fact2*xvert(4) + fact3*xmid
         ynode(isub) = fact1*yvert(1) + fact2*yvert(4) + fact3*ymid
         jsub(isub) = last + isub
         ksub(isub+degree+2) = jsub(isub)
         isub = isub+1
      end do
   end do
   nsub = 2*(degree-1) + ((degree-2)*(degree-1))/2
   if (nsub == 0) then
      last = 0
   else
      last = jsub(nsub)
   endif
   if (nsub /= 0) then
      call p_hier_basis_func(xnode(1:nsub),ynode(1:nsub), &
                             (/xvert(1),xvert(4),xmid/), &
                             (/yvert(1),yvert(4),ymid/), &
                            (/degree,degree,degree,degree/),"a",basis(:,1:nsub))
   endif
   do i=1,nsub
      new_soln(jsub(i),:) = 0
      do j=1,nbasis
         new_soln(jsub(i),:) = new_soln(jsub(i),:) + solution(ksub(j),:)*basis(j,i)
      end do
   end do

! same process for the fourth triangle if it exists

   ksub(1) = 2
   ksub(2) = 4
   ksub(3) = 4+degree
   do i=1,degree-1
      ksub(i+3) = (degree+1)**2 + 1 + i
   end do
   do i=1,degree-1
      ksub(degree+2+i) = 2 + ((degree+1)*(degree+2))/2 + i
   end do
   isub = 1
   do i=1,degree-1
      fact = i/rdeg
      xnode(isub) = fact*xvert(4) + (1-fact)*xvert(2)
      ynode(isub) = fact*yvert(4) + (1-fact)*yvert(2)
      jsub(isub) = last + isub
      ksub(isub+2*degree+1) = jsub(isub)
      isub = isub+1
   end do
   do i=1,degree-2
      fact3 = i/rdeg
      do j=1,degree-1-i
         fact2 = (1-fact3)*j/(rdeg-1)
         fact1 = 1-fact3-fact2
         xnode(isub) = fact1*xvert(2) + fact2*xvert(4) + fact3*xmid
         ynode(isub) = fact1*yvert(2) + fact2*yvert(4) + fact3*ymid
         jsub(isub) = last + isub
         ksub(isub+2*degree+1) = jsub(isub)
         isub = isub+1
      end do
   end do
   nsub = (degree-1) + ((degree-2)*(degree-1))/2
   if (nsub /= 0) then
      call p_hier_basis_func(xnode(1:nsub),ynode(1:nsub), &
                             (/xvert(2),xvert(4),xmid/), &
                             (/yvert(2),yvert(4),ymid/), &
                            (/degree,degree,degree,degree/),"a",basis(:,1:nsub))
   endif
   do i=1,nsub
      new_soln(jsub(i),:) = 0
      do j=1,nbasis
         new_soln(jsub(i),:) = new_soln(jsub(i),:) + solution(ksub(j),:)*basis(j,i)
      end do
   end do
endif

! copy new coefficients to return variable

solution = new_soln

end subroutine phier2nodal

!          -----------
subroutine nodal2phier(solution,xvert,yvert,degree)
!          -----------

!----------------------------------------------------
! This routine converts the values in solution from coefficients of a
! nodal basis to coefficients of a p-hierarchical basis, over two triangles
! that are mates with the outer vertices given in (xvert,yvert).  The order
! of the coefficients (using cubics as an example) is

!              vertex 1
!                  1
!                / | \
!               7  |  14
!              /   |   \
!             8    9    15
!            /     |     \
! vertex 3  3  11  |  16  4  vertex 4
!            \     |     /
!             6    10   13
!              \   |   /
!               5  |  12
!                \ | /
!                  2
!              vertex 2

! For the p-hierarchical basis, those on the same edge or face functions of
! the same element are in order of degree, as in the basis function routines.
! All elements and edges must have the same degree.  The second triangle
! can be omitted by setting the 4th vertex to be huge(0.0_my_real).

! If the second triangle is omitted, then all the indices after 3 are
! decremented by 1.

!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(inout) :: solution(:,:)
real(my_real), intent(in) :: xvert(:),yvert(:)
integer, intent(in) :: degree
!----------------------------------------------------
! Local variables:

real(my_real) :: xnode((degree+1)**2),ynode((degree+1)**2), &
                 basis((degree+1)**2,(degree+1)**2), &
                 basistemp(((degree+1)*(degree+2))/2,((degree+1)*(degree+2))/2)
real(my_real) :: rdeg,fact,fact1,fact2,fact3
integer :: index1(((degree+1)*(degree+2))/2),index2(((degree+1)*(degree+2))/2),&
           ipiv((degree+1)**2)
integer :: i,j,isub,jsub,nnode,info
!----------------------------------------------------
! Begin executable code

rdeg = real(degree,my_real)

! define the nodes

! nodes that are vertices (will fix 4 later if only one triangle)

xnode(1:4) = xvert(1:4)
ynode(1:4) = yvert(1:4)
index1(1:3) = (/1,2,3/)
index2(1:3) = (/1,2,4/)

! nodes in first triangle

isub = 5
jsub = 4
do i=1,degree-1
   fact = i/rdeg
   xnode(isub) = fact*xvert(3) + (1-fact)*xvert(2)
   ynode(isub) = fact*yvert(3) + (1-fact)*yvert(2)
   index1(jsub) = isub
   jsub = jsub+1
   isub = isub+1
end do
do i=1,degree-1
   fact = i/rdeg
   xnode(isub) = fact*xvert(3) + (1-fact)*xvert(1)
   ynode(isub) = fact*yvert(3) + (1-fact)*yvert(1)
   index1(jsub) = isub
   jsub = jsub+1
   isub = isub+1
end do
do i=1,degree-1
   fact = i/rdeg
   xnode(isub) = fact*xvert(2) + (1-fact)*xvert(1)
   ynode(isub) = fact*yvert(2) + (1-fact)*yvert(1)
   index1(jsub) = isub
   jsub = jsub+1
   isub = isub+1
end do
do i=1,degree-2
   fact3 = i/rdeg
   do j=1,degree-1-i
      fact2 = (1-fact3)*j/(rdeg-1)
      fact1 = 1-fact3-fact2
      xnode(isub) = fact1*xvert(1) + fact2*xvert(2) + fact3*xvert(3)
      ynode(isub) = fact1*yvert(1) + fact2*yvert(2) + fact3*yvert(3)
      index1(jsub) = isub
      jsub = jsub+1
      isub = isub+1
   end do
end do

! if only one triangle, done, but remove vertex 4

if (xvert(4) == huge(0.0_my_real)) then
   nnode = isub-2
   xnode(4:nnode) = xnode(5:nnode+1)
   ynode(4:nnode) = ynode(5:nnode+1)
   index1(4:nnode) = index1(4:nnode)-1

! otherwise, nodes in second triangle

else

   jsub = 4
   do i=1,degree-1
      fact = i/rdeg
      xnode(isub) = fact*xvert(4) + (1-fact)*xvert(2)
      ynode(isub) = fact*yvert(4) + (1-fact)*yvert(2)
      index2(jsub) = isub
      jsub = jsub+1
      isub = isub+1
   end do
   do i=1,degree-1
      fact = i/rdeg
      xnode(isub) = fact*xvert(4) + (1-fact)*xvert(1)
      ynode(isub) = fact*yvert(4) + (1-fact)*yvert(1)
      index2(jsub) = isub
      jsub = jsub+1
      isub = isub+1
   end do
   do i=1,degree-1
      index2(jsub) = 2 + 2*degree + i
      jsub = jsub+1
   end do
   do i=1,degree-2
      fact3 = i/rdeg
      do j=1,degree-1-i
         fact2 = (1-fact3)*j/(rdeg-1)
         fact1 = 1-fact3-fact2
         xnode(isub) = fact1*xvert(1) + fact2*xvert(2) + fact3*xvert(4)
         ynode(isub) = fact1*yvert(1) + fact2*yvert(2) + fact3*yvert(4)
         index2(jsub) = isub
         jsub = jsub+1
         isub = isub+1
      end do
   end do
   nnode = isub-1

endif

! evaluate the p-hierarchical basis functions at the nodes

basis = 0
call p_hier_basis_func(xnode(index1),ynode(index1),xvert(1:3),yvert(1:3), &
                       (/degree,degree,degree,degree/),"a",basistemp)
basis(index1,index1) = basistemp
if (xvert(4) /= huge(0.0_my_real)) then
   call p_hier_basis_func(xnode(index2),ynode(index2), &
                          (/xvert(1),xvert(2),xvert(4)/), &
                          (/yvert(1),yvert(2),yvert(4)/), &
                          (/degree,degree,degree,degree/),"a",basistemp)
   basis(index2,index2) = basistemp
endif

! if n is the vector of nodal basis coefficients, p is the vector of
! p-hierarchical basis coefficients, and A is the matrix with A_ij having
! the value of the jth p-hierarchical basis at the ith node, then
! n = Ap.  Therefore, to get p from n, compute p = A^-1 n.  Note however that
! the matrix basis is A^T.

basis = transpose(basis)

if (my_real == kind(0.0)) then
   call sgesv(nnode,size(solution,dim=2),basis,size(basis,dim=1),ipiv, &
              solution,size(solution,dim=1),info)
else
   call dgesv(nnode,size(solution,dim=2),basis,size(basis,dim=1),ipiv, &
              solution,size(solution,dim=1),info)
endif

if (info /= 0) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("lapack failed in nodal2phier",intlist=(/info/))
   stop
endif

end subroutine nodal2phier

!          -----------
subroutine nodal2hhier(solution,xvert,yvert,degree)
!          -----------

!----------------------------------------------------
! This routine converts the values in solution from coefficients of a
! nodal basis to coefficients of a h-hierarchical basis, over four triangles
! that are siblings in an h refinement with the outer vertices given in
! (xvert,yvert).  The order of the coefficients is the same as in phier2nodal.
! All elements and edges must have the same degree.  The third and fourth
! triangles can be omitted by setting the 4th vertex to be huge(0.0_my_real).
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(inout) :: solution(:,:)
real(my_real), intent(in) :: xvert(:),yvert(:)
integer, intent(in) :: degree
!----------------------------------------------------
! Local variables:

integer :: nred1, nblack1, nbasis1, nred2, nblack2, nbasis2, nred4, nblack4, &
           nbasis4, astat, i, j, k, rsub, bsub, ind
logical :: red
real(my_real) :: xvert5, yvert5, rdeg, w1, w2, w3
real(my_real) :: s(size(solution,dim=1),size(solution,dim=1)), &
                 prod(size(solution,dim=1),size(solution,dim=2))
real(my_real), allocatable :: xnode(:), ynode(:), basis(:,:)
integer, allocatable :: ind_black(:), ind_red(:)
!----------------------------------------------------
! Begin executable code

! some useful constants

rdeg = real(degree,my_real)
xvert5 = (xvert(1)+xvert(2))/2
yvert5 = (yvert(1)+yvert(2))/2

! number of basis functions / nodes in 1 triangle and all 2 or 4 triangles

nbasis1 = ((degree+1)*(degree+2))/2
if (2*(degree/2) == degree) then
   nred1 = (degree/2)*(degree/2+1)
else
   nred1 = ((degree+1)/2)**2
endif
nblack1 = nbasis1 - nred1
nred2 = (degree*(degree+1))/2
nblack2 = nbasis1
nbasis2 = nred2 + nblack2
nred4 = degree**2
nblack4 = (degree+1)**2
nbasis4 = nred4 + nblack4
if (nbasis4 /= size(solution,dim=1)) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("size mismatch in nodal2hhier")
   stop
endif

! allocate memory based on number of bases/nodes

allocate(xnode(nred2),ynode(nred2),basis(nblack2,nred2),ind_black(nblack2), &
         ind_red(nred2),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in nodal2hhier")
   return
endif

! The (i,j)th entry of the matrix that converts nodal coefficients to
! hierarchical coefficients is the negative of the jth h-hierarchical basis
! function at the ith node if j is black and i is red, and identity on the
! diagonal.  The black h-hierarchical bases are the nodal bases of the parent
! triangles.  Go through the 1 or 2 parents evaluating the nodal bases at
! the red nodes and inserting the value in the right place.

s = 0
do i=1,nbasis4
   s(i,i) = 1
end do

! parent

! define the red nodes, red node indices, and black basis indices; go along
! lines parallel to (vert1,vert2)

rsub = 1
bsub = 1
! vertex 1
ind_black(bsub) = 1
bsub = bsub + 1
! nodes along line between child 1 and child 3
red = .true.
ind = 5 + degree-1 + 1
do i=1,degree-1
   if (red) then
      w3 = i/rdeg
      xnode(rsub) = (1-w3)*xvert(1) + w3*xvert5
      ynode(rsub) = (1-w3)*yvert(1) + w3*yvert5
      ind_red(rsub) = ind
      rsub = rsub + 1
   else
      ind_black(bsub) = ind
      bsub = bsub + 1
   endif
   ind = ind + 1
   red = .not. red
end do
! central node
if (red) then
   xnode(rsub) = xvert5
   ynode(rsub) = yvert5
   ind_red(rsub) = degree + 4
   rsub = rsub + 1
else
   ind_black(bsub) = degree + 4
   bsub = bsub + 1
endif
red = .not. red
! nodes along edge between child 2 and child 4, in reverse order
ind = 5 + 4*(degree-1) + ((degree-2)*(degree-1))/2
do i=1,degree-1
   if (red) then
      w1 = i/rdeg
      xnode(rsub) = w1*xvert(2) + (1-w1)*xvert5
      ynode(rsub) = w1*yvert(2) + (1-w1)*yvert5
      ind_red(rsub) = ind
      rsub = rsub + 1
   else
      ind_black(bsub) = ind
      bsub = bsub + 1
   endif
   ind = ind - 1
   red = .not. red
end do
! vertex 2
ind_black(bsub) = 2
bsub = bsub + 1
! for each line parallel to (vert1, vert2)
do k=1,degree-1
   w2 = k/rdeg
! node on edge between vert1 and vert3
   ind_black(bsub) = 5 + 2*(degree-1) + k
   bsub = bsub + 1
! interior nodes in child 1
   red = .true.
   ind = 5 + 3*(degree-1) + k
   do j=1,degree-1-k
      w3 = (1-w2)*j/real(degree-k,my_real)
      w1 = 1-w2-w3
      if (red) then
         xnode(rsub) = w1*xvert(1) + w2*xvert(3) + w3*xvert5
         ynode(rsub) = w1*yvert(1) + w2*yvert(3) + w3*yvert5
         ind_red(rsub) = ind
         rsub = rsub + 1
      else
         ind_black(bsub) = ind
         bsub = bsub + 1
      endif
      ind = ind + degree - j - 1
      red = .not. red
   end do
! line between child 1 and child 2
   if (red) then
      xnode(rsub) = w2*xvert(3) + (1-w2)*xvert5
      ynode(rsub) = w2*yvert(3) + (1-w2)*yvert5
      ind_red(rsub) = 5 + degree-1 - k
      rsub = rsub + 1
   else
      ind_black(bsub) = 5 + degree-1 - k
      bsub = bsub + 1
   endif
   red = .not. red
! interior nodes in child 2
   ind = 5 + 5*(degree-1) + (degree-2)*(degree-1) - ((k-1)*k)/2
   do j=1,degree-1-k
      w1 = (1-w2)*j/real(degree-k,my_real)
      w3 = 1-w2-w1
      if (red) then
         xnode(rsub) = w1*xvert(2) + w2*xvert(3) + w3*xvert5
         ynode(rsub) = w1*yvert(2) + w2*yvert(3) + w3*yvert5
         ind_red(rsub) = ind
         rsub = rsub + 1
      else
         ind_black(bsub) = ind
         bsub = bsub + 1
      endif
      ind = ind - (k+j)
      red = .not. red
   end do
! node on line between vert2 and vert3
   ind_black(bsub) = 5 + 4*(degree-1) + ((degree-2)*(degree-1))/2 + k
   bsub = bsub + 1
end do ! lines parallel to (vert1, vert2)
! vertex 3
ind_black(bsub) = 3
bsub = bsub + 1

! evaluate the nodal basis functions of the parent at the red nodes

call nodal_basis_func(xnode,ynode,(/xvert(1),xvert(2),xvert(3)/), &
                      (/yvert(1),yvert(2),yvert(3)/),degree,"a",basis)

! copy the values into the matrix

do j=1,nblack2
   do i=1,nred2
      s(ind_red(i),ind_black(j)) = -basis(j,i)
   end do
end do

if (xvert(4) /= huge(0.0_my_real)) then

! mate

! define the red nodes, red node indices, and black basis indices; go along
! lines parallel to (vert1,vert2)

   rsub = 1
   bsub = 1
! vertex 1
   ind_black(bsub) = 1
   bsub = bsub + 1
! nodes along line between child 1 and child 3
   red = .true.
   ind = 5 + degree-1 + 1
   do i=1,degree-1
      if (red) then
         w3 = i/rdeg
         xnode(rsub) = (1-w3)*xvert(1) + w3*xvert5
         ynode(rsub) = (1-w3)*yvert(1) + w3*yvert5
         ind_red(rsub) = ind
         rsub = rsub + 1
      else
         ind_black(bsub) = ind
         bsub = bsub + 1
      endif
      ind = ind + 1
      red = .not. red
   end do
! central node
   if (red) then
      xnode(rsub) = xvert5
      ynode(rsub) = yvert5
      ind_red(rsub) = degree + 4
      rsub = rsub + 1
   else
      ind_black(bsub) = degree + 4
      bsub = bsub + 1
   endif
   red = .not. red
! nodes along edge between child 2 and child 4, in reverse order
   ind = 5 + 4*(degree-1) + ((degree-2)*(degree-1))/2
   do i=1,degree-1
      if (red) then
         w1 = i/rdeg
         xnode(rsub) = w1*xvert(2) + (1-w1)*xvert5
         ynode(rsub) = w1*yvert(2) + (1-w1)*yvert5
         ind_red(rsub) = ind
         rsub = rsub + 1
      else
         ind_black(bsub) = ind
         bsub = bsub + 1
      endif
      ind = ind - 1
      red = .not. red
   end do
! vertex 2
   ind_black(bsub) = 2
   bsub = bsub + 1
! for each line parallel to (vert1, vert2)
   do k=1,degree-1
      w2 = k/rdeg
! node on edge between vert1 and vert4
      ind_black(bsub) = 5 + 6*(degree-1) + (degree-2)*(degree-1) + k
      bsub = bsub + 1
! interior nodes in child 3
      red = .true.
      ind = 5 + 7*(degree-1) + (degree-2)*(degree-1) + k
      do j=1,degree-1-k
         w3 = (1-w2)*j/real(degree-k,my_real)
         w1 = 1-w2-w3
         if (red) then
            xnode(rsub) = w1*xvert(1) + w2*xvert(4) + w3*xvert5
            ynode(rsub) = w1*yvert(1) + w2*yvert(4) + w3*yvert5
            ind_red(rsub) = ind
            rsub = rsub + 1
         else
            ind_black(bsub) = ind
            bsub = bsub + 1
         endif
         ind = ind + degree - j - 1
         red = .not. red
      end do
! line between child 3 and child 4
      if (red) then
         xnode(rsub) = w2*xvert(4) + (1-w2)*xvert5
         ynode(rsub) = w2*yvert(4) + (1-w2)*yvert5
         ind_red(rsub) = 6 + 6*(degree-1) + (degree-2)*(degree-1) - k
         rsub = rsub + 1
      else
         ind_black(bsub) = 6 + 6*(degree-1) + (degree-2)*(degree-1) - k
         bsub = bsub + 1
      endif
      red = .not. red
! interior nodes in child 4
      ind = 5 + 8*(degree-1) + 2*(degree-2)*(degree-1) - ((k-1)*k)/2
      do j=1,degree-1-k
         w1 = (1-w2)*j/real(degree-k,my_real)
         w3 = 1-w2-w1
         if (red) then
            xnode(rsub) = w1*xvert(2) + w2*xvert(4) + w3*xvert5
            ynode(rsub) = w1*yvert(2) + w2*yvert(4) + w3*yvert5
            ind_red(rsub) = ind
            rsub = rsub + 1
         else
            ind_black(bsub) = ind
            bsub = bsub + 1
         endif
         ind = ind - (k+j)
         red = .not. red
      end do
! node on line between vert2 and vert4
      ind_black(bsub) = 5 + 7*(degree-1) + 3*((degree-2)*(degree-1))/2 + k
      bsub = bsub + 1
   end do ! lines parallel to (vert1, vert2)
! vertex 4
   ind_black(bsub) = 4
   bsub = bsub + 1

! evaluate the nodal basis functions of the mate at the red nodes

   call nodal_basis_func(xnode,ynode,(/xvert(1),xvert(2),xvert(4)/), &
                         (/yvert(1),yvert(2),yvert(4)/),degree,"a",basis)

! copy the values into the matrix

   do j=1,nblack2
      do i=1,nred2
         s(ind_red(i),ind_black(j)) = -basis(j,i)
      end do
   end do

endif ! xvert(4) /= huge

! convert values in solution to h-hierarchical

if (my_real == kind(1.0)) then
   call sgemm("N","N",nbasis4,size(solution,dim=2),nbasis4,1.0, &
              s,nbasis4,solution,nbasis4,0.0,prod,nbasis4)
elseif (my_real == kind(1.0d0)) then
   call dgemm("N","N",nbasis4,size(solution,dim=2),nbasis4,1.0d0, &
              s,nbasis4,solution,nbasis4,0.0d0,prod,nbasis4)
else
   ierr = PHAML_INTERNAL_ERROR
   call fatal("my_real is neither single nor double precision. Can't call GEMM")
   stop
endif

! copy result back to solution

solution = prod

deallocate(xnode,ynode,basis,ind_black,ind_red)

end subroutine nodal2hhier

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
   call get_grid_info(grid,procs,still_sequential,444, &
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
if (refcont%reftype == H_UNIFORM .or. &
    refcont%reftype == P_UNIFORM .or. &
    global_max_errind == 0.0_my_real) is_uniform = .true.
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

!          ------------
subroutine mark_reftype(grid,refine_control,global_max_errind,reftype,numhref,&
                        numpref,elist)
!          ------------

!----------------------------------------------------
! This routine marks each element for h, p or no refinement
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(refine_options), intent(in) :: refine_control
real(my_real), intent(in) :: global_max_errind
character(len=*), intent(inout) :: reftype(:)
integer, intent(in) :: numhref(:), numpref(:)
type(errind_list), intent(in) :: elist

!----------------------------------------------------
! Local variables:

integer :: lev, elem
!----------------------------------------------------
! Begin executable code

! first pass, mark all elements, but delay the decision for expensive strategies

do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      call mark_reftype_one(elem,grid,refine_control,global_max_errind,reftype,&
                            numhref,numpref,.true.)
      elem = grid%element(elem)%next
   end do ! next element
end do ! next level

! second pass, mark all delayed elements in bin 1 without delay

elem = elist%head_errind(1)
do while (elem /= END_OF_LIST)
   if (reftype(elem) == "u") then
      call mark_reftype_one(elem,grid,refine_control,global_max_errind,reftype,&
                            numhref,numpref,.false.)
   endif
   elem = elist%next_errind(elem)
end do

end subroutine mark_reftype

!          ----------------
subroutine mark_reftype_one(elem,grid,refine_control,global_max_errind,reftype,&
                            numhref,numpref,delay)
!          ----------------

!----------------------------------------------------
! This routine marks one element for h, p or no refinement.
! If delay is true, don't make the decision now if it is expensive.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: elem
type(grid_type), intent(inout) :: grid
type(refine_options), intent(in) :: refine_control
real(my_real), intent(in) :: global_max_errind
character(len=*), intent(inout) :: reftype(:)
integer, intent(in) :: numhref(:),numpref(:)
logical, intent(in) :: delay

!----------------------------------------------------
! Local variables:

integer :: add
real(my_real) :: energy_errind(grid%system_size), err, err2, reg
!----------------------------------------------------
! Begin executable code

! initalize to no refinement

reftype(elem) = "n"

! only refine leaves that I own

if (grid%element(elem)%isleaf .and. grid%element(elem)%iown) then

! if the global maximum error indicator is 0, refine everything by h

   if (global_max_errind == 0.0_my_real) then
      reftype(elem) = "h"
   else

      select case (refine_control%reftype)

! h refinement refines all elements by h

      case (H_UNIFORM, H_ADAPTIVE)
         reftype(elem) = "h"

! p refinement refines all elements by p

      case (P_UNIFORM, P_ADAPTIVE)
         reftype(elem) = "p"

! hp strategies

      case (HP_ADAPTIVE)

! for all strategies, use p if h is maxed out or visa versa

         if (grid%element(elem)%level >= refine_control%max_lev .and. &
             grid%element(elem)%degree >= refine_control%max_deg) then
            reftype(elem) = "n"
         elseif (grid%element(elem)%level >= refine_control%max_lev) then
            reftype(elem) = "p"
         elseif (grid%element(elem)%degree >= refine_control%max_deg) then
            reftype(elem) = "h"
         else

            select case(refine_control%hp_strategy)

! bigger errind computes the local problem h and p error indicators and
! refines by the one that is larger.  Experimentally I found that multiplying
! the p error indicator by 2 is better than just taking the largest and
! better than weighting by the increase in degrees of freedom on most of
! the test problems.  (Weighting by increase in dof was better on elasticity
! and peak.)

            case (HP_BIGGER_ERRIND)

               if (delay) then
                  reftype(elem) = "u"
               else
                  call error_indicator(grid,elem,LOCAL_PROBLEM_H, &
                                       energy=energy_errind)
                  err = maxval(energy_errind)
! use this to weight by increase in dof
!                  add = ((grid%element(elem)%degree-1)* &
!                         (grid%element(elem)%degree-2))/2 + &
!                        2*grid%element(elem)%degree-1
!                  if (.not. grid%element(elem)%mate == BOUNDARY) then
!                     mate = hash_decode_key(grid%element(elem)%mate, &
!                                            grid%elem_hash)
!                     add = add + ((grid%element(mate)%degree-1)* &
!                                  (grid%element(mate)%degree-2))/2 + &
!                                 grid%element(mate)%degree-1
!                  endif
!                  err = err/add
! end of weighting by increase in dof
                  call error_indicator(grid,elem,LOCAL_PROBLEM_P, &
                                       energy=energy_errind)
                  err2 = maxval(energy_errind)
! use this to weight by increase in dof
!                  add = 3+max(0,grid%element(elem)%degree-2)
!                  err2 = err2/add
! end of weighting by increase in dof
! use this to bias toward p-refinement
                  err2 = 2*err2
! end of biasing toward p-refinement
                  if (err > err2) then
                     reftype(elem) = "h"
                  else
                     reftype(elem) = "p"
                  endif
               endif

! regularity-based strategies compute a (strategy dependent) estimate of
! the regularity and use p refinement where the estimate is sufficiently
! large

            case (HP_APRIORI, HP_PRIOR2P_E, HP_PRIOR2P_H1, HP_COEF_ROOT, &
                  HP_NEXT3P)

               if (delay) then
                  reftype(elem) = "u"
               else

                  select case(refine_control%hp_strategy)
                  case (HP_APRIORI)
                     reg = regularity( &
                             grid%vertex(grid%element(elem)%vertex)%coord%x, &
                             grid%vertex(grid%element(elem)%vertex)%coord%y)
                     add = 1
                  case (HP_PRIOR2P_E)
                     reg = prior2p_e_regularity(grid,elem)
                     add = 1
                  case (HP_PRIOR2P_H1)
                     reg = prior2p_h1_regularity(grid,elem)
                     add = 1
                  case (HP_COEF_ROOT)
                     reg = coef_root_regularity(grid,elem)
                     add = 1
                  case (HP_NEXT3P)
                     reg = next3p_regularity(grid,elem)
                     add = 1
                  end select

                  if (grid%element(elem)%degree + add <= reg) then
                     reftype(elem) = "p"
                  else
                     reftype(elem) = "h"
                  endif

               endif

! Type parameter uses h if e(t,p)/e(t,p-1) > gamma and p otherwise

            case (HP_TYPEPARAM)

               if (delay) then
                  reftype(elem) = "u"
               else
                  if (typeparam_R(grid,elem) > refine_control%tp_gamma) then
                     reftype(elem) = "h"
                  else
                     reftype(elem) = "p"
                  endif
               endif

! Texas 3 Step and Alternate do everything the same; which one was determined
! back in phaml_solve_pde

            case (HP_T3S, HP_ALTERNATE)
               if (refine_control%t3s_reftype == H_UNIFORM) then
                  reftype(elem) = "h"
               elseif (refine_control%t3s_reftype == H_ADAPTIVE) then
                  reftype(elem) = "h"
               elseif (refine_control%t3s_reftype == P_ADAPTIVE) then
                  reftype(elem) = "p"
               else
                  ierr = PHAML_INTERNAL_ERROR
                  call fatal("bad value for t3s_reftype")
                  stop
               endif

! Coefficient decay uses p if the decay rate is bigger than 1

            case (HP_COEF_DECAY)
               if (delay) then
                  reftype(elem) = "u"
               else
                  if (coef_decay_rate(grid,elem) > 1.0_my_real) then
                     reftype(elem) = "p"
                  else
                     reftype(elem) = "h"
                  endif
               endif

! Smooth prediction uses h if the error estimate is larger than the
! predicted error estimate
! TEMP only using first eigenvector

            case (HP_SMOOTH_PRED)
               if (grid%element_errind(elem,1) > grid%element(elem)%sp_eta_pred) then
                  reftype(elem) = "h"
               else
                  reftype(elem) = "p"
               endif

            case default
               ierr = PHAML_INTERNAL_ERROR
               call fatal("bad case of hp_strategy under adaptive refinement")
               stop
            end select

         endif ! reached max h or max p

      case default
         ierr = PHAML_INTERNAL_ERROR
         call fatal("bad case of reftype under adaptive refinement")
         stop
      end select
   endif ! global_max_errind == 0
endif ! owned leaf

end subroutine mark_reftype_one

!          ----------
subroutine set_numref(grid,procs,still_sequential,refine_control,numhref, &
                      numpref)
!          ----------

!----------------------------------------------------
! This routine sets the number of times each element should be refined
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential
type(refine_options), intent(in) :: refine_control
integer, intent(inout) :: numhref(:), numpref(:)
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

! -1 indicates the chosen scheme does not use this.

if (refine_control%reftype /= HP_ADAPTIVE) then
   numhref = -1
   numpref = -1

else

   select case (refine_control%hp_strategy)

   case default
      numhref = -1
      numpref = -1

   case (HP_T3S)

      call t3s_set_numref(grid,procs,refine_control,numhref,numpref)

   end select

endif

end subroutine set_numref

!          --------------
subroutine t3s_set_numref(grid,procs,refine_control,numhref,numpref)
!          --------------

!----------------------------------------------------
! This routine sets the number of times each element should be refined for
! Texas 3 Step
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(proc_info), intent(in) :: procs
type(refine_options), intent(in) :: refine_control
integer, intent(inout) :: numhref(:), numpref(:)
!----------------------------------------------------
! Local variables:

real(my_real) :: theta_I, log2, r, p_0_k, theta_I_k, theta_T, nu_k
real(my_real), dimension(size(grid%element)) :: mu_k, theta_0_k
integer :: N_I, N_I_own, lev, elem, i
!----------------------------------------------------
! Begin executable code

! The uniform refinement phase of Texas 3 Step doesn't use it

if (refine_control%t3s_reftype == H_UNIFORM) then
   numhref = -1
   numpref = -1

! Part I: determine the number of h refinements to do
! ------

elseif (refine_control%t3s_reftype == H_ADAPTIVE) then

! set numhref to 0, just to be safe about nonleaves and to not refine
! unowned elements.  numpref is 0

   numhref = 0
   numpref = 0

! compute error indicators if needed

   if (.not. grid%errind_up2date) then
      call all_error_indicators(grid,refine_control%error_estimator)
   endif

! set constants that are independent of element

   theta_I = refine_control%t3s_h_target

! for first guess at numhref, assume N_I is twice N_0

   N_I_own = 2*grid%nelem_leaf_own
   N_I = phaml_global_sum(procs,N_I_own,410)
   N_I_own = 0

! for each owned leaf element ...

   do lev=1,grid%nlev
      elem = grid%head_level_elem(lev)
      do while (elem /= END_OF_LIST)
         if (grid%element(elem)%isleaf .and. &
             grid%element(elem)%iown) then

! set the constants for this element

! get the regularity, which should be min(s+2,1+alpha) where rhs is in H^s(T)
! and alpha = pi/(angle of reentrant corner)

            r = regularity(grid%vertex(grid%element(elem)%vertex)%coord%x, &
                           grid%vertex(grid%element(elem)%vertex)%coord%y)

            mu_k(elem) = min(real(grid%element(elem)%degree,my_real),r-1)
! TEMP only using first eigenvector
            theta_0_k(elem) = grid%element_errind(elem,1)

            numhref(elem) = max(1.0_my_real, &
                        (theta_0_k(elem)**2*N_I/theta_I**2)**(1/(1+mu_k(elem))))
            N_I_own = N_I_own + numhref(elem)

         endif ! is owned leaf
         elem = grid%element(elem)%next
      end do ! next element
   end do ! next level

! iterate to balance N_I and sum(numhref)

   do i=1,5
      N_I = phaml_global_sum(procs,N_I_own,410+i)
      N_I_own = 0
      do lev=1,grid%nlev
         elem = grid%head_level_elem(lev)
         do while (elem /= END_OF_LIST)
            if (grid%element(elem)%isleaf .and. &
                grid%element(elem)%iown) then
               numhref(elem) = max(1.0_my_real, &
                        (theta_0_k(elem)**2*N_I/theta_I**2)**(1/(1+mu_k(elem))))
               N_I_own = N_I_own + numhref(elem)
            endif ! is owned leaf
            elem = grid%element(elem)%next
         end do ! next element
      end do ! next level
   end do

! number of refinements is the base-2 log of the number of desired children.
! round and cap it at t3s_maxref

   log2 = log(2.0_my_real)
   do lev=1,grid%nlev
      elem = grid%head_level_elem(lev)
      do while (elem /= END_OF_LIST)
         if (grid%element(elem)%isleaf .and. &
             grid%element(elem)%iown) then
            if (numhref(elem) < 2) then
               numhref(elem) = 0
            else
               numhref(elem) = 0.5 + log(real(numhref(elem),my_real))/log2
            endif
            numhref(elem) = min(numhref(elem),refine_control%t3s_maxref)
         endif ! is owned leaf
         elem = grid%element(elem)%next
      end do ! next element
   end do ! next level

! Part II: determine the number of p refinement to do
! -------

elseif (refine_control%t3s_reftype == P_ADAPTIVE) then

! set numpref to 0, just to be safe about nonleaves and to not refine
! unowned elements.  numhref is 0

   numhref = 0
   numpref = 0

! compute error indicators if needed

   if (.not. grid%errind_up2date) then
      call all_error_indicators(grid,refine_control%error_estimator)
   endif

! set constants that are independent of element

   theta_T = refine_control%t3s_p_target
   N_I_own = grid%nelem_leaf_own
   N_I = phaml_global_sum(procs,N_I_own,460)

! for each owned leaf element ...

   do lev=1,grid%nlev
      elem = grid%head_level_elem(lev)
      do while (elem /= END_OF_LIST)
         if (grid%element(elem)%isleaf .and. &
             grid%element(elem)%iown) then

! set the constants for this element

! get the regularity, which should be min(s+2,1+alpha) where rhs is in H^s(T)
! and alpha = pi/(angle of reentrant corner)

            r = regularity(grid%vertex(grid%element(elem)%vertex)%coord%x, &
                           grid%vertex(grid%element(elem)%vertex)%coord%y)

! TEMP for smooth areas, say r>10, use nu=1, because I really don't understand

            if (r > 10.0_my_real) then
               nu_k = 1
            else
               nu_k = r-1
            endif
! TEMP only using first eigenvector
            theta_I_k = grid%element_errind(elem,1)
            p_0_k = grid%element(elem)%degree

            numpref(elem) = p_0_k*(N_I*theta_I_k**2/theta_T**2)**(1/(2*nu_k))

         endif ! is owned leaf
         elem = grid%element(elem)%next
      end do ! next element
   end do ! next level

! number of refinements is the difference between the old and new p_k

   do lev=1,grid%nlev
      elem = grid%head_level_elem(lev)
      do while (elem /= END_OF_LIST)
         if (grid%element(elem)%isleaf .and. &
             grid%element(elem)%iown) then
            numpref(elem) = max(numpref(elem)-grid%element(elem)%degree,0)
            numpref(elem) = min(numpref(elem),refine_control%t3s_maxdeginc)
         endif
         elem = grid%element(elem)%next
      end do
   end do

! bad value for t3s_reftype

else
   ierr = PHAML_INTERNAL_ERROR
   call fatal("bad value for t3s_reftype in t3s_set_numref")
   stop
endif

end subroutine t3s_set_numref

!        ----------------
function element_diameter(grid,elem)
!        ----------------

!----------------------------------------------------
! This routine computes the diameter of element elem
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: elem
real(my_real) :: element_diameter
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

! take the diameter to be the longest side length

element_diameter = &
   sqrt((grid%vertex(grid%element(elem)%vertex(2))%coord%x - &
         grid%vertex(grid%element(elem)%vertex(1))%coord%x)**2 + &
        (grid%vertex(grid%element(elem)%vertex(2))%coord%y - &
         grid%vertex(grid%element(elem)%vertex(1))%coord%y)**2)

element_diameter = max(element_diameter, &
   sqrt((grid%vertex(grid%element(elem)%vertex(3))%coord%x - &
         grid%vertex(grid%element(elem)%vertex(2))%coord%x)**2 + &
        (grid%vertex(grid%element(elem)%vertex(3))%coord%y - &
         grid%vertex(grid%element(elem)%vertex(2))%coord%y)**2))

element_diameter = max(element_diameter, &
   sqrt((grid%vertex(grid%element(elem)%vertex(1))%coord%x - &
         grid%vertex(grid%element(elem)%vertex(3))%coord%x)**2 + &
        (grid%vertex(grid%element(elem)%vertex(1))%coord%y - &
         grid%vertex(grid%element(elem)%vertex(3))%coord%y)**2))

end function element_diameter

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
case (H_UNIFORM, P_UNIFORM)
   return_to_elist = .false.
   complete_elist = .true.
   one_elist = .true.
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
   end select
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

!          -----------
subroutine extend_nlev(grid)
!          -----------

!----------------------------------------------------
! This routine increases the number of refinement levels supported in
! the grid data structure
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid

!----------------------------------------------------
! Local variables:

integer, pointer :: old_array(:)
integer :: oldlev, newlev, astat, dstat

!----------------------------------------------------
! Begin executable code

oldlev = size(grid%head_level_elem)
newlev = 2*oldlev
old_array => grid%head_level_elem
allocate(grid%head_level_elem(newlev),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in extend_nlev")
   grid%head_level_elem => old_array
   return
endif
grid%head_level_elem(1:oldlev) = old_array
grid%head_level_elem(oldlev+1:newlev) = END_OF_LIST
deallocate(old_array,stat=dstat)
old_array => grid%head_level_vert
allocate(grid%head_level_vert(newlev),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in extend_nlev")
   grid%head_level_vert => old_array
   return
endif
grid%head_level_vert(1:oldlev) = old_array
grid%head_level_vert(oldlev+1:newlev) = END_OF_LIST
deallocate(old_array,stat=dstat)

end subroutine extend_nlev

!          -------------
subroutine more_elements(grid,errcode,elist,reftype,numhref,numpref)
!          -------------

!----------------------------------------------------
! This routine increases the space for the number of elements
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(out) :: errcode
type(errind_list), optional, intent(inout) :: elist
character(len=*), optional, pointer :: reftype(:)
integer, optional, pointer :: numhref(:), numpref(:)
!----------------------------------------------------
! Local variables:

type(element_t), pointer :: temp_elem(:)
integer, pointer :: temp_eprev(:), temp_enext(:)
character(len=1), pointer :: temp_reftype(:)
integer, pointer :: temp_numref(:)
real(my_real), pointer :: temp_errind(:,:)
integer :: allocstat, oldsize, newsize, i
!----------------------------------------------------
! Begin executable code

oldsize = size(grid%element)
newsize = int(1.5*oldsize)

errcode = 0
temp_elem => grid%element
nullify(grid%element)
allocate(grid%element(newsize),stat=allocstat)
if (allocstat /= 0) then
   errcode = 1
   call fatal("increased allocation failed")
   grid%element => temp_elem
   return
endif

temp_errind => grid%element_errind
nullify(grid%element_errind)
allocate(grid%element_errind(newsize,max(1,grid%num_eval)),stat=allocstat)
if (allocstat /= 0) then
   errcode = 1
   call fatal("increased allocation failed")
   grid%element => temp_elem
   grid%element_errind => temp_errind
   return
endif

if (present(elist)) then
   temp_eprev => elist%prev_errind
   temp_enext => elist%next_errind
   nullify(elist%prev_errind, elist%next_errind)
   allocate(elist%prev_errind(newsize),elist%next_errind(newsize),stat=allocstat)
   if (allocstat /= 0) then
      errcode = 1
      call fatal("increased allocation failed")
      elist%prev_errind => temp_eprev
      elist%next_errind => temp_enext
      deallocate(grid%element,stat=allocstat)
      grid%element => temp_elem
      grid%element_errind => temp_errind
      return
   endif
   elist%prev_errind(1:oldsize) = temp_eprev
   elist%next_errind(1:oldsize) = temp_enext
   elist%prev_errind(oldsize+1:newsize) = NOT_ON_LIST
   elist%next_errind(oldsize+1:newsize) = NOT_ON_LIST
   deallocate(temp_eprev,temp_enext,stat=allocstat)
endif

if (present(reftype)) then
   temp_reftype => reftype
   nullify(reftype)
   allocate(reftype(newsize),stat=allocstat)
   if (allocstat /= 0) then
      errcode = 1
      call fatal("increased allocation failed")
      stop
   endif
   reftype(1:oldsize) = temp_reftype
   reftype(oldsize+1:newsize) = "n"
   deallocate(temp_reftype,stat=allocstat)
endif

if (present(numhref)) then
   temp_numref => numhref
   nullify(numhref)
   allocate(numhref(newsize),stat=allocstat)
   if (allocstat /= 0) then
      errcode = 1
      call fatal("increased allocation failed")
      stop
   endif
   numhref(1:oldsize) = temp_numref
   numhref(oldsize+1:newsize) = -1
   deallocate(temp_numref,stat=allocstat)
endif

if (present(numpref)) then
   temp_numref => numpref
   nullify(numpref)
   allocate(numpref(newsize),stat=allocstat)
   if (allocstat /= 0) then
      errcode = 1
      call fatal("increased allocation failed")
      stop
   endif
   numpref(1:oldsize) = temp_numref
   numpref(oldsize+1:newsize) = -1
   deallocate(temp_numref,stat=allocstat)
endif

grid%element(1:oldsize) = temp_elem
grid%element(oldsize+1:newsize)%degree = 1
do i=oldsize+1,newsize
   nullify(grid%element(i)%solution,grid%element(i)%exact, &
           grid%element(i)%oldsoln)
end do
grid%element_errind(1:oldsize,:) = temp_errind
grid%element_errind(oldsize+1:newsize,:) = 0.0_my_real

deallocate(temp_elem,temp_errind,stat=allocstat)
grid%next_free_elem = oldsize+1
grid%element(oldsize+1:newsize)%previous = (/ (i,i=oldsize,newsize-1) /)
grid%element(oldsize+1:newsize-1)%next = (/ (i,i=oldsize+2,newsize) /)
grid%element(newsize)%next = END_OF_LIST

!call count_memory(grid)

end subroutine more_elements

!          ----------
subroutine more_edges(grid,errcode)
!          ----------

!----------------------------------------------------
! This routine increases the space for the number of edges
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(out) :: errcode
!----------------------------------------------------
! Local variables:

type(edge_t), pointer :: temp_edge(:)
integer, pointer :: temp_edge_type(:,:)
integer :: allocstat, oldsize, newsize, i
!----------------------------------------------------
! Begin executable code

oldsize = size(grid%edge)
newsize = int(1.5*oldsize)

errcode = 0
temp_edge => grid%edge
temp_edge_type => grid%edge_type
nullify(grid%edge,grid%edge_type)
allocate(grid%edge(newsize),grid%edge_type(newsize,grid%system_size), &
         stat=allocstat)
if (allocstat /= 0) then
   errcode = 1
   call fatal("increased allocation failed")
   grid%edge => temp_edge
   grid%edge_type => temp_edge_type
   return
endif
grid%edge(1:oldsize) = temp_edge
grid%edge(oldsize+1:newsize)%degree = 1
do i=oldsize+1,newsize
   nullify(grid%edge(i)%solution,grid%edge(i)%exact,grid%edge(i)%oldsoln)
end do
grid%edge_type(1:oldsize,:) = temp_edge_type
grid%edge_type(oldsize+1:newsize,:) = 0

deallocate(temp_edge,temp_edge_type,stat=allocstat)
grid%next_free_edge = oldsize+1
grid%edge(oldsize+1:newsize-1)%next = (/ (i,i=oldsize+2,newsize) /)
grid%edge(newsize)%next = END_OF_LIST

!call count_memory(grid)

end subroutine more_edges

!          ----------
subroutine more_verts(grid,errcode)
!          ----------

!----------------------------------------------------
! This routine increases the space for the number of vertices
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(out) :: errcode
!----------------------------------------------------
! Local variables:

type(vertex_t), pointer :: temp_vertex(:)
integer, pointer :: temp_vertex_type(:,:)
real(my_real), pointer :: temp_vertex_solution(:,:,:), &
                          temp_vertex_exact(:,:,:), temp_vertex_oldsoln(:,:,:)
integer :: astat, astat2, astat3, oldsize, newsize, i
!----------------------------------------------------
! Begin executable code

oldsize = size(grid%vertex)
newsize = int(1.5*oldsize)

errcode = 0
temp_vertex => grid%vertex
temp_vertex_type => grid%vertex_type
temp_vertex_solution => grid%vertex_solution
temp_vertex_exact => grid%vertex_exact
temp_vertex_oldsoln => grid%vertex_oldsoln
nullify(grid%vertex,grid%vertex_type)
allocate(grid%vertex(newsize),grid%vertex_type(newsize,grid%system_size), &
         grid%vertex_solution(newsize,grid%system_size,max(1,grid%num_eval)), &
         stat=astat)
if (grid%have_true) then
   allocate(grid%vertex_exact(newsize,grid%system_size,max(1,grid%num_eval)), &
            stat=astat2)
else
   astat2 = 0
endif
if (associated(grid%vertex_oldsoln)) then
   allocate(grid%vertex_oldsoln(newsize,grid%system_size,max(1,grid%num_eval)), &
            stat=astat3)
else
   astat3 = 0
endif
if (astat /= 0 .or. astat2 /= 0 .or. astat3 /= 0) then
   errcode = 1
   call fatal("increased allocation failed")
   grid%vertex => temp_vertex
   grid%vertex_type => temp_vertex_type
   grid%vertex_solution => temp_vertex_solution
   grid%vertex_exact => temp_vertex_exact
   grid%vertex_oldsoln => temp_vertex_oldsoln
   return
endif
grid%vertex(1:oldsize) = temp_vertex
grid%vertex_type(1:oldsize,:) = temp_vertex_type
grid%vertex_type(oldsize+1:newsize,:) = 0
grid%vertex_solution(1:oldsize,:,:) = temp_vertex_solution
grid%vertex_solution(oldsize+1:newsize,:,:) = 0
if (grid%have_true) then
   grid%vertex_exact(1:oldsize,:,:) = temp_vertex_exact
   grid%vertex_exact(oldsize+1:newsize,:,:) = 0
endif
if (associated(grid%vertex_oldsoln)) then
   grid%vertex_oldsoln(1:oldsize,:,:) = temp_vertex_oldsoln
   grid%vertex_oldsoln(oldsize+1:newsize,:,:) = 0
endif

deallocate(temp_vertex,temp_vertex_type,temp_vertex_solution,stat=astat)
if (associated(temp_vertex_exact)) deallocate(temp_vertex_exact,stat=astat)
if (associated(temp_vertex_oldsoln)) deallocate(temp_vertex_oldsoln,stat=astat)
grid%next_free_vert = oldsize+1
grid%vertex(oldsize+1:newsize)%previous = (/ (i,i=oldsize,newsize-1) /)
grid%vertex(oldsize+1:newsize-1)%next = (/ (i,i=oldsize+2,newsize) /)
grid%vertex(newsize)%next = END_OF_LIST

!call count_memory(grid)

end subroutine more_verts

!          ------------
subroutine count_memory(grid)
!          ------------

!----------------------------------------------------
! This routine adds up the amount of memory used in the grid data structure.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
!----------------------------------------------------
! Local variables:

integer :: isize, rsize, lsize, lssize, csize, p0size, p1size, p2size, &
           p3size, p4size, hsize
integer :: mem, i
!----------------------------------------------------
! Begin executable code

isize = bit_size(1)/8

! TEMP need a good portable way to determine the number of bytes in real.
! I think the following should always work, at least on the computers of 2006.
! The standard intrinsic function digits returns the number of binary digits in
! the real model for the type of real given.  The number of binary digits cannot
! possibly be bigger than the number of bits, so we look for the first bit
! size that can hold that number of binary digits.  I assume the number of
! bytes is a power of 2, it's at least 4 and no more than 64.

rsize = digits(1.0_my_real)
if (rsize <= 32) then
   rsize = 4
elseif (rsize <= 64) then
   rsize = 8
elseif (rsize <= 128) then
   rsize = 16
elseif (rsize <= 256) then
   rsize = 32
else
   rsize = 64
endif

! TEMP assume a default logical is 4 bytes and a small logical is 1 byte

lsize = 4
lssize = 1

! TEMP assume a default character is 1 byte

csize = 1

! The size of pointer descriptors varies among compilers.  There is also
! lost memory for small pointer arrays because of alignment requirements;
! this can be very large (lf95 can loose up to 192 bytes for every allocated
! pointer).  I do not include the memory lost to alignment in the memory
! count in this routine.  The following byte counts for pointer descriptors
! and alignment requirements were determined experimentally for some Fortran
! compilers on Linux.
!          Lahey   Intel   g95   Absoft   Nag
! Rank 0     8       4      4      24      4
! Rank 1    28      36     28      36     20
! Rank 2    44      48     40      48     32
! Rank 3    60      60     52      60     44
! Rank 4    76      72     64      72     56
! alignment 192     16     32      96      8

p0size = 8
p1size = 28
p2size = 44
p3size = 60
p4size = 76

! hash key size

hsize = KEY_SIZE*isize

mem = 0

! entries in grid_type

mem = mem + 3*p1size  & ! pointers for element, edge and vertex
          + 3*0       & ! hash tables ! TEMP look for size
          + 2*2*rsize & ! bounding box
          + 3*p3size  & ! solutions
          + 2*p2size  & ! error indicators
          + p1size    & ! pointer for eigenvalue
          + 2*rsize   & ! eigen linsys resids
          + 8*p1size  & ! pointers for eigenprob and errests
          + rsize     & ! max_blen
          + 2*p1size  & ! pointers for bp
          + 2*p2size  & ! pointers for type
          + p2size    & ! pointer for initial neighbor
          + 2*p1size  & ! pointer for linked list heads
          + 23*isize  & ! whole mess of integers
          + 2*lsize   & ! 2 logicals
          + FN_LEN*csize

! allocated pointers in grid_type, except element, edge and vertex

if (associated(grid%vertex_solution)) mem = mem + rsize*size(grid%vertex_solution)
if (associated(grid%vertex_exact)) mem = mem + rsize*size(grid%vertex_exact)
if (associated(grid%vertex_oldsoln)) mem = mem + rsize*size(grid%vertex_oldsoln)
if (associated(grid%element_errind)) mem = mem + rsize*size(grid%element_errind)
if (associated(grid%eigenvalue)) mem = mem + rsize*size(grid%eigenvalue)
if (associated(grid%eigenprob_l2_resid)) mem = mem + rsize*size(grid%eigenprob_l2_resid)
if (associated(grid%eigenprob_variance)) mem = mem + rsize*size(grid%eigenprob_variance)
if (associated(grid%errest_energy)) mem = mem + rsize*size(grid%errest_energy)
if (associated(grid%errest_Linf)) mem = mem + rsize*size(grid%errest_Linf)
if (associated(grid%errest_L2)) mem = mem + rsize*size(grid%errest_L2)
if (associated(grid%errest_eigenvalue)) &
   mem = mem + rsize*size(grid%errest_eigenvalue)
if (associated(grid%bp_start)) mem = mem + rsize*size(grid%bp_start)
if (associated(grid%bp_finish)) mem = mem + rsize*size(grid%bp_finish)
if (associated(grid%edge_type)) mem = mem + isize*size(grid%edge_type)
if (associated(grid%vertex_type)) mem = mem + isize*size(grid%vertex_type)
if (associated(grid%initial_neighbor)) mem = mem + isize*size(grid%initial_neighbor)
if (associated(grid%head_level_elem)) mem = mem + isize*size(grid%head_level_elem)
if (associated(grid%head_level_vert)) mem = mem + isize*size(grid%head_level_vert)

! elements

if (associated(grid%element)) then
   mem = mem + size(grid%element)*( &
                  hsize*2                    & ! gid and mate
                + rsize                      & ! weight
                + p2size*3                   & ! pointers for solutions
                + p1size*2                   & ! pointers for error indicators
                + rsize*2                    & ! works
                + isize*VERTICES_PER_ELEMENT & ! vertex
                + isize*EDGES_PER_ELEMENT    & ! edge
                + isize*(4+MAX_CHILD+2)      & ! bunch of integers
                + lssize*5                   & ! small logicals
               )

   do i=1,size(grid%element)
      if (associated(grid%element(i)%solution)) mem = mem + rsize*size(grid%element(i)%solution)
      if (associated(grid%element(i)%exact)) mem = mem + rsize*size(grid%element(i)%exact)
      if (associated(grid%element(i)%oldsoln)) mem = mem + rsize*size(grid%element(i)%oldsoln)
   end do

endif

! edges

if (associated(grid%edge)) then
   mem = mem + size(grid%edge)*( &
                  hsize &    ! gid
                + isize*5 &  ! integers
                + p1size &   ! pointer for type
                + p2size*3 & ! pointers for solutions
                + isize    & ! next
               )

   do i=1,size(grid%edge)
      if (associated(grid%edge(i)%solution)) mem = mem + rsize*size(grid%edge(i)%solution)
      if (associated(grid%edge(i)%exact)) mem = mem + rsize*size(grid%edge(i)%exact)
      if (associated(grid%edge(i)%oldsoln)) mem = mem + rsize*size(grid%edge(i)%oldsoln)
   end do

endif

! vertices

if (associated(grid%vertex)) then
   mem = mem + size(grid%vertex)*( &
                  hsize &    ! gid
                + rsize*2 &  ! coord
                + rsize &    ! bparam
                + p1size &   ! pointer for type
                + isize*4 &  ! integers
               )
endif

write(outunit,"(A,I12,A)") "Memory for grid data structure ",mem," bytes"

end subroutine count_memory

!---------------------------------------------------------------------
!  ROUTINES COMMON TO REFSOLN_EDGE AND REFSOLN_ELEM
!---------------------------------------------------------------------

!          --------------
subroutine refine_refsoln(grid,procs,refine_control,solver_control, &
                          io_control,still_sequential,init_nvert, &
                          init_nelem,init_dof,loop,balance_what,predictive)
!          --------------

!----------------------------------------------------
! This routine is the top level refine routine for the REFSOLN_* hp strategies.
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

!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

select case (refine_control%hp_strategy)
case (HP_REFSOLN_EDGE)
   call refine_refsoln_edge(grid,procs,refine_control,solver_control, &
                          io_control,still_sequential,init_nvert, &
                          init_nelem,init_dof,loop,balance_what,predictive)
case (HP_REFSOLN_ELEM)
   call refine_refsoln_elem(grid,procs,refine_control,solver_control, &
                          io_control,still_sequential,init_nvert, &
                          init_nelem,init_dof,loop,balance_what,predictive)
case default
   ierr = PHAML_INTERNAL_ERROR
   call fatal("bad case of hp_strategy in refine_refsoln")
   stop
end select

end subroutine refine_refsoln

!---------------------------------------------------------------------
!  ROUTINES FOR THE REFSOLN_EDGE HP-ADAPTIVE STRATEGY
!---------------------------------------------------------------------

!          -------------------
subroutine refine_refsoln_edge(grid,procs,refine_control,solver_control, &
                               io_control,still_sequential,init_nvert, &
                               init_nelem,init_dof,loop,balance_what,predictive)
!          -------------------

!----------------------------------------------------
! This routine is the top level refine routine for the REFSOLN_EDGE hp strategy
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

!----------------------------------------------------
! Local variables:

type(grid_type) :: old_grid, ref_soln
type(refine_options) :: loc_refcont
type(io_options) :: loc_iocont
!----------------------------------------------------
! Begin executable code

! TEMP only using first component

   if (grid%system_size /= 1) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("REFSOLN_EDGE not yet implemented for multicomponent systems")
      stop
   endif

! TEMP only using first eigenvector

   if (grid%num_eval > 1) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("REFSOLN_EDGE not yet implemented for multiple eigenvalues")
      stop
   endif

! TEMP not parallel

   if (num_proc(procs) > 1) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("REFSOLN_EDGE not yet implemented in parallel")
      stop
   endif

! start timing the refinement process

   call reset_watch(prefine)
   call start_watch((/prefine,trefine/))

! Apply the Dirichlet kludge

   call dirichlet_kludge(grid,refine_control)

! Make a copy of the grid as it currently exists.  Note the copy is in old_grid.

   call copy_grid(grid,old_grid)

! Create the reference solution by performing uniform h and p refinements
! and solving on the fine grid.

   loc_refcont = refine_control
   loc_refcont%error_estimator = HIERARCHICAL_COEFFICIENT
   loc_refcont%reftype = H_UNIFORM
   loc_iocont = io_control
   loc_iocont%print_linsys_when = NEVER
   loc_iocont%print_error_when = NEVER
   call refine(grid,procs,loc_refcont,solver_control,loc_iocont, &
               still_sequential,init_nvert,init_nelem,init_dof,loop, &
               balance_what,predictive,no_time=.true.)
   loc_refcont%reftype = P_UNIFORM
   call refine(grid,procs,loc_refcont,solver_control,loc_iocont, &
               still_sequential,init_nvert,init_nelem,init_dof,loop, &
               balance_what,predictive,no_time=.true.)
   call solve(grid,procs,loc_iocont,solver_control,still_sequential,.false., &
              no_time=.true.)
   call copy_grid(grid,ref_soln)

! MASTER doesn't participate further; it was only here to monitor solve

if (my_proc(procs) == MASTER) then
   call deallocate_grid(old_grid)
   call deallocate_grid(ref_soln)
   grid%errind_up2date = .false.
   call stop_watch((/prefine,trefine/))
   return
endif

! Compute the error estimate.

   grid%refsoln_errest = compute_refsoln_errest_edge(ref_soln,old_grid)

! If the error estimate is large enough, determine the optimal refinements
! and unrefine grid to meet them.

   if (grid%refsoln_errest > refine_control%term_energy_err) then
      call det_and_perf_opt_ref_edge(grid,old_grid,ref_soln,refine_control)
   else
      call deallocate_grid(grid)
      call copy_grid(old_grid,grid)
   endif

! Destroy the old grid and reference solution

   call deallocate_grid(old_grid)
   call deallocate_grid(ref_soln)

! Error indicators have not been (and won't really be) set

   grid%errind_up2date = .false.

! stop timing the refinement process

   call stop_watch((/prefine,trefine/))

end subroutine refine_refsoln_edge

!          ----------------
subroutine dirichlet_kludge(grid,refine_control)
!          ----------------

!----------------------------------------------------
! If a first degree element has its base on a Dirichlet boundary with a linear
! function as the boundary condition, and is in a region where only p refinement
! will occur, then it will never be refined.  This is because the edge will
! not be selected for refinement because the solution is exact on it, the
! other edges won't be selected because they are not bases and won't become
! bases unless the neighbor is h refined, and will keep degree 1 because of
! the minimum rule, and the element won't have its degree increased because
! increasing from 1 to 2 doesn't add any face bases and the edges are unchanged
! so the element rate is 0.
!
! This routine removes this situation by performing h refinement of any edge
! that is on a Dirichlet boundary with the boundary condition at the midpoint
! of the edge equal to the average of the solution at the endpoints before
! the refsoln_edge algorithm even starts.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(refine_options), intent(in) :: refine_control
!----------------------------------------------------
! Local variables:

integer :: lev, elem, edge, vert1, vert2, itype(grid%system_size), errcode
real(my_real) :: xmid, ymid, umid, rs(grid%system_size), &
                 c(grid%system_size,grid%system_size)
!----------------------------------------------------
! Begin executable code

! for each leaf element

do lev=1,grid%nlev
   elem=grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (grid%element(elem)%isleaf) then

! if the base is a Dirichlet edge

         edge = grid%element(elem)%edge(3)
         if (grid%edge_type(edge,1) == DIRICHLET) then

! get the boundary condition at the midpoint

            vert1 = grid%edge(edge)%vertex(1)
            vert2 = grid%edge(edge)%vertex(2)
            xmid = (grid%vertex(vert1)%coord%x + grid%vertex(vert2)%coord%x)/2
            ymid = (grid%vertex(vert1)%coord%y + grid%vertex(vert2)%coord%y)/2
            call bconds(xmid,ymid,grid%edge(edge)%bmark,itype,c,rs)

! if it is the average of the solution at the endpoints

            umid = (grid%vertex_solution(vert1,1,1) + &
                    grid%vertex_solution(vert2,1,1))/2
            if (abs(rs(1)-umid) < 100*epsilon(0.0_my_real)) then

! h-refine the element

               call bisect_triangle_pair(grid,elem,errcode,refine_control)

! next element

            endif
         endif
      endif
      elem = grid%element(elem)%next
   end do
end do

end subroutine dirichlet_kludge

!        ---------------------------
function compute_refsoln_errest_edge(ref_soln,old_grid)
!        ---------------------------

!----------------------------------------------------
! This routine computes the error estimate
!    |u_{h,p} - u_{h/sqrt(2),p+1}| / |u_{h/sqrt(2),p+1}|
! where |.| is the H^1 seminorm
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: ref_soln, old_grid
real(my_real) :: compute_refsoln_errest_edge
!----------------------------------------------------
! Local variables:

integer :: lev, elem, nqpoints, jerr, parent, old_elem, grandparent
real(my_real) :: norm_fine
real(my_real), pointer :: qweights(:), xquad(:), yquad(:)
real(my_real), allocatable :: ux_fine(:,:,:), uy_fine(:,:,:), ux_old(:,:,:), &
                              uy_old(:,:,:)
type(hash_key) :: parent_gid, grandparent_gid
!----------------------------------------------------
! Begin executable code

! initalize the norm of the difference and solution

   compute_refsoln_errest_edge = 0
   norm_fine = 0

!  for each leaf element of the fine grid

   do lev=1,ref_soln%nlev
      elem = ref_soln%head_level_elem(lev)
      do while (elem /= END_OF_LIST)
         if (ref_soln%element(elem)%iown .and. ref_soln%element(elem)%isleaf) then

! determine a quadrature rule for p+1 (which is p in the refined grid)

            call quadrature_rule_tri(ref_soln%element(elem)%degree, &
                       ref_soln%vertex(ref_soln%element(elem)%vertex)%coord%x, &
                       ref_soln%vertex(ref_soln%element(elem)%vertex)%coord%y, &
                       nqpoints,qweights,xquad,yquad,jerr,.true.)

! identify the element in old_grid that contains this element in the fine grid,
! as the grandparent if it is a leaf and parent otherwise

            parent_gid = ref_soln%element(elem)%gid/2
            parent = hash_decode_key(parent_gid,ref_soln%elem_hash)
            if (ref_soln%element(parent)%level == 1) then
               old_elem = parent
            else
               grandparent_gid = ref_soln%element(parent)%gid/2
               grandparent = hash_decode_key(grandparent_gid,ref_soln%elem_hash)
               if (old_grid%element(grandparent)%isleaf) then
                  old_elem = grandparent
               else
                  old_elem = parent
               endif
            endif

! evaluate both solutions

            allocate(ux_fine(1,1,nqpoints),uy_fine(1,1,nqpoints), &
                     ux_old(1,1,nqpoints),uy_old(1,1,nqpoints))
            call evaluate_soln_local(ref_soln,xquad,yquad,elem,(/1/),(/1/), &
                                     ux=ux_fine,uy=uy_fine)
            call evaluate_soln_local(old_grid,xquad,yquad,old_elem,(/1/),(/1/),&
                                     ux=ux_old,uy=uy_old)

! add contributions to integrals

            compute_refsoln_errest_edge = compute_refsoln_errest_edge + &
                           sum(qweights*((ux_old(1,1,:)-ux_fine(1,1,:))**2 + &
                                         (uy_old(1,1,:)-uy_fine(1,1,:))**2))
            norm_fine = norm_fine + &
                        sum(qweights*(ux_fine(1,1,:)**2+uy_fine(1,1,:)**2))

            deallocate(ux_fine,uy_fine,ux_old,uy_old)

! next element

         endif
         elem = ref_soln%element(elem)%next
      end do
   end do

! square roots and normalization
! Use abs because negative quadrature weights can cause negative squared norm

   compute_refsoln_errest_edge = sqrt(abs(compute_refsoln_errest_edge))
   if (norm_fine /= 0.0_my_real) then
      compute_refsoln_errest_edge = compute_refsoln_errest_edge/sqrt(abs(norm_fine))
   endif

end function compute_refsoln_errest_edge

!          -----------------------
subroutine ext_refsoln_errest_edge(grid,procs,refine_control,solver_control, &
                                   io_control,still_sequential)
!          -----------------------

!----------------------------------------------------
! This routine can be called from outside the module to compute the
! error estimate used by REFSOLN_EDGE
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type (proc_info), target, intent(in) :: procs
type(refine_options), intent(in) :: refine_control
type(solver_options), intent(in) :: solver_control
type(io_options), intent(in) :: io_control
logical, intent(in) :: still_sequential
!----------------------------------------------------
! Local variables:

type(grid_type) :: ref_soln
type(refine_options) :: loc_refcont
type(io_options) :: loc_iocont
!----------------------------------------------------
! Begin executable code

! TEMP only using first component

if (grid%system_size /= 1) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("REFSOLN_EDGE not yet implemented for multicomponent systems")
   stop
endif

! TEMP only using first eigenvector

if (grid%num_eval > 1) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("REFSOLN_EDGE not yet implemented for multiple eigenvalues")
   stop
endif

! TEMP not parallel

if (num_proc(procs) > 1) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("REFSOLN_EDGE not yet implemented in parallel")
   stop
endif

! Make a copy of the grid as it currently exists.

call copy_grid(grid,ref_soln)

! Create the reference solution by performing uniform h and p refinements
! and solving on the fine grid.

loc_refcont = refine_control
loc_refcont%error_estimator = HIERARCHICAL_COEFFICIENT
loc_refcont%reftype = H_UNIFORM
loc_iocont = io_control
loc_iocont%print_linsys_when = NEVER
loc_iocont%print_error_when = NEVER
call refine(ref_soln,procs,loc_refcont,solver_control,loc_iocont, &
            still_sequential,0,0,0,0,0,.false.,no_time=.true.)
loc_refcont%reftype = P_UNIFORM
call refine(ref_soln,procs,loc_refcont,solver_control,loc_iocont, &
            still_sequential,0,0,0,0,0,.false.,no_time=.true.)
call solve(ref_soln,procs,loc_iocont,solver_control,still_sequential,.false., &
           no_time=.true.)

! MASTER doesn't participate further; it was only here to monitor solve

if (my_proc(procs) == MASTER) then
   call deallocate_grid(ref_soln)
   grid%errind_up2date = .true.
   return
endif

! Compute the error estimate.

grid%refsoln_errest = compute_refsoln_errest_edge(ref_soln,grid)
call deallocate_grid(ref_soln)
grid%errind_up2date = .true.

end subroutine ext_refsoln_errest_edge

!          -------------------------
subroutine det_and_perf_opt_ref_edge(grid,old_grid,ref_soln,refine_control)
!          -------------------------

!----------------------------------------------------
! This is the top level routine for determining the optimal refinements to
! perform, and performing those refinements.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(grid_type), intent(in) :: old_grid, ref_soln
type(refine_options), intent(in) :: refine_control
!----------------------------------------------------
! Local variables:

character(len=1) :: edge_reftype(size(old_grid%edge))
real(my_real) :: guaranteed_edge_rate(size(old_grid%edge)), &
                 diff_old(size(old_grid%edge)), edge_rate_max
integer :: best_h(size(old_grid%edge))
integer :: lev, elem, edge
!----------------------------------------------------
! Begin executable code

! initially set every edge refinement type to no and guaranteed edge rate
! to a number smaller than any that will be computed

   edge_reftype = "n"
   guaranteed_edge_rate = -huge(0.0_my_real)

! for each base edge in old_grid

   do lev=1,old_grid%nlev
      elem = old_grid%head_level_elem(lev)
      do while (elem /= END_OF_LIST)
         if (old_grid%element(elem)%isleaf) then
            edge = old_grid%element(elem)%edge(3)

! decide if this edge should be refined by p or h (if it is refined),
! identify the best h refinement, and compute the guaranteed edge rate.
! If edge_reftype is not "n", this edge has already been done.

            if (edge_reftype(edge) == "n") then
               call choose_edge_p_or_h(edge,elem,ref_soln,old_grid, &
                                       refine_control, &
                                       edge_reftype(edge),best_h(edge), &
                                       diff_old(edge))
               call compute_guaranteed_edge_rate(ref_soln,old_grid,edge,elem, &
                                                 best_h(edge),diff_old(edge), &
                                                 guaranteed_edge_rate(edge))
            endif

! next edge

         endif
         elem = old_grid%element(elem)%next
      end do
   end do

! determine the maximum edge rate

   edge_rate_max = maxval(guaranteed_edge_rate)

! determine which edges should be refined

   call determine_edges_to_refine(edge_reftype,guaranteed_edge_rate, &
                                  edge_rate_max)

! in the fine grid, un-h-refine elements for which the base is not marked
! for h refinement

   call unrefine_edge_h_refinements(grid,old_grid,edge_reftype,refine_control)

! set the optimal degrees of the edges

   call det_and_set_opt_deg_edges(grid,old_grid,ref_soln,refine_control, &
                                  edge_reftype,best_h,edge_rate_max,diff_old)

! set the optimal degrees of the elements

   call det_and_set_opt_deg_elem(grid,old_grid,ref_soln,refine_control, &
                                 edge_rate_max)

! enforce the minimum rule for edge degree

   call impose_minimum_rule(grid,refine_control)

end subroutine det_and_perf_opt_ref_edge

!          ------------------
subroutine choose_edge_p_or_h(edge,elem,ref_soln,old_grid,refine_control, &
                              edge_reftype,best_h,diff_old)
!          ------------------

!----------------------------------------------------
! This routine determines whether an edge should be refined by p or h,
! and identifies the best h refinement.
! It also returns the norm of the difference between the old and reference
! solutions so it doesn't have to be computed again for guaranteed rate.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: edge, elem
type(grid_type), intent(in) :: ref_soln, old_grid
type(refine_options), intent(in) :: refine_control
character(len=1), intent(out) :: edge_reftype
integer, intent(out) :: best_h
real(my_real), intent(out) :: diff_old
!----------------------------------------------------
! Local variables:

integer :: p(2), i
real(my_real) :: rate_p, rate_h, rate_i
real(my_real) :: diff_hp(3), diff_p(3), diff_h(3)
!----------------------------------------------------
! Begin executable code

! compute |u_fine - w_{hp}|^2 where |.| is the weighted H^1 seminorm
! and w_{hp} is the projection based interpolant on edge of the old grid

   p(1) = old_grid%edge(edge)%degree
   p(2) = 0
   call weighted_H1_seminorm_diff_squared(ref_soln,old_grid,edge,elem,p,diff_hp)

! Save the norm to return it.
! Note the norm of the difference is in the third component (the first two
! are used for child contributions)

   diff_old = diff_hp(3)

! compute |u_fine - w_p|^2 where w_p is the projection based interpolant
! on the p refined edge

   p(1) = p(1) + 1
   call weighted_H1_seminorm_diff_squared(ref_soln,old_grid,edge,elem,p,diff_p)

! compute the error decrease rate for the p refined edge

   rate_p = diff_hp(3) - diff_p(3)

! for each of the p possible h-refinements, compute the error decrease rate
! and identify the best h refinement

   rate_h = -huge(0.0_my_real)
   do i=1,old_grid%edge(edge)%degree

! compute rate

      p(1) = i
      p(2) = old_grid%edge(edge)%degree + 1 - i
      call weighted_H1_seminorm_diff_squared(ref_soln,old_grid,edge,elem,p,diff_h)
      rate_i = diff_hp(3) - diff_h(3)

! check for largest rate

      if (rate_i > rate_h) then
         best_h = i
         rate_h = rate_i
      endif

! next possible h-refinement

   end do

! determine the edge refinement type based on largest error decrease rate

   if (old_grid%element(elem)%level >= refine_control%max_lev .and. &
       old_grid%edge(edge)%degree >= refine_control%max_deg) then
      edge_reftype = "n"
   elseif (old_grid%element(elem)%level >= refine_control%max_lev) then
      edge_reftype = "p"
   elseif (old_grid%edge(edge)%degree >= refine_control%max_deg) then
      edge_reftype = "h"
   elseif (rate_p > rate_h) then
      edge_reftype = "p"
   else
      edge_reftype = "h"
   endif

end subroutine choose_edge_p_or_h

!          ---------------------------------
subroutine weighted_H1_seminorm_diff_squared(ref_soln,old_grid,edge,elem,p,diff)
!          ---------------------------------

!----------------------------------------------------
! This routine computes |u_fine - w|^2 where |.| is the weighted H1 seminorm
! on edge edge in old_grid and w is the projection based interpolation of the
! reference solution in grid onto edge refined as specified in p.  If
! p(2)==0, then the projection is onto a single edge of degree p(1).  Otherwise
! it is onto two edges, formed by bisecting edge, with degrees p(1) and p(2).
! The edge is the base of elem in old_grid.
!
! On return, diff(3) contains the norm.  If p(2) /= 0, diff(1) and diff(2)
! contain the child contributions.
!
!
! The interpolant w = u_1 + u_2 is defined by solving the variational problem
! on each child edge e_i, i=1..2, with the understanding that e_1 is edge if
! p(2)==0:
!   u_2 in V(e_i)
!   <u_2,v>_e_i = <(u_fine-u_1),v>_e_i  forall v in V(e_i)
! where V(e_i) is the space of bubble function over e_i up to degree p(i),
! and the inner product is the weighted H^1 seminorm (xi is the Greek letter)
! <u,v>_e = integral_0^1 du/dxi dv/dxi dxi
! with x(xi) = x_1 + (x_2-x_1)xi, y(xi) = y_1 + (y_2-y_1)xi
! where (x_1,y_1) and (x_2,y_2) are the endpoints of e_i.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: ref_soln, old_grid
integer, intent(in) :: edge, elem, p(2)
real(my_real), intent(out) :: diff(3)
!----------------------------------------------------
! Local variables:

integer :: endpt(3), n, nmax, qorder, nqpoints, jerr, loop, qp, i, j, &
           halfqp, child1, child2
integer, allocatable :: ipiv(:)
real(my_real) :: du1dxi, xline(2), yline(2), xtri(3), ytri(3), deltax, deltay, &
                 integrand
real(my_real), allocatable :: a(:,:), b(:,:), dphidxi(:,:), dudx(:,:,:), &
                              dudy(:,:,:)
real(my_real), pointer :: qweights(:), xquad(:), yquad(:), xquadu(:), &
                          yquadu(:), qweightshalf(:), xquadhalf(:), yquadhalf(:)
!----------------------------------------------------
! Begin executable code

! initialize result

   diff = 0.0_my_real

! determine the elements in the fine grid that contain the children of edge

! first guess is the children of elem

   child1 = hash_decode_key(2*ref_soln%element(elem)%gid,ref_soln%elem_hash)
   child2 = hash_decode_key(2*ref_soln%element(elem)%gid+1,ref_soln%elem_hash)
   if (child1 == HASH_NOT_FOUND .or. child2 == HASH_NOT_FOUND) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("didn't find child element in weighted_H1_seminorm_diff_squared")
      stop
   endif

! Endpoints of the edges.  1 is first vertex of edge, 2 is midpoint of edge,
! 3 is second vertex of edge

   endpt(1) = old_grid%edge(edge)%vertex(1)
   endpt(2) = ref_soln%element(child1)%vertex(3)
   endpt(3) = old_grid%edge(edge)%vertex(2)

! if either of the children was further refined (for compatibility), then
! the element to use is its first child

   if (.not. ref_soln%element(child1)%isleaf) then
      child1 = hash_decode_key(2*ref_soln%element(child1)%gid, &
                                ref_soln%elem_hash)
   endif
   if (.not. ref_soln%element(child2)%isleaf) then
      child2 = hash_decode_key(2*ref_soln%element(child2)%gid, &
                                ref_soln%elem_hash)
   endif
   if (child1 == HASH_NOT_FOUND .or. child2 == HASH_NOT_FOUND) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("didn't find child element in weighted_H1_seminorm_diff_squared")
      stop
   endif

! (xtri,ytri) is a triangle containing the unit interval as the first side,
! used for evaluating basis functions as a function of (Greek letter) xi

   xtri(1) = 0.0_my_real
   ytri(1) = 0.0_my_real
   xtri(2) = 1.0_my_real
   ytri(2) = 0.0_my_real
   xtri(3) = 0.5_my_real
   ytri(3) = 1.0_my_real

! Determine the maximum size of the linear systems to solve and allocate
! space for them.

   nmax = max(p(1)-1,p(2)-1)
   allocate(a(nmax,nmax),b(nmax,1),ipiv(nmax))

! compute the integrals over each subinterval, or the whole interval

   do loop=1,2
      if (loop==2 .and. p(2)==0) exit

! size of the linear system for this (sub)interval

      n = p(loop)-1

! Note that in the following, p(loop) = n+1 allows us to have p(2)=p(1) when
! p(2)==0
! quadratures for the matrix involve polynomials of degree 2*(p(loop)-1), and
! quadratures for the right side involve polynomials of degree
! (1+degree(u_fine)-1)+(p(loop)-1).
! quadratures for the final integral involve polynomials of degree
! 2*max(degree(u_fine),p(loop)-1)
! Pick a quadrature rule that will be exact for all of them

      qorder = max(2*n,(old_grid%edge(edge)%degree+1)+n)
      qorder = max(qorder,2*(old_grid%edge(edge)%degree+1))

! the quadrature rule is exact for polynomials up to order 2*qorder-1, so
! get the order from the degree (currently in qorder) by (qorder+1)/2, except
! add 2 because of truncation in integer division

      qorder = min((qorder+2)/2,MAX_QUAD_ORDER_LINE)

! u_1 is the linear vertex interpolant of u_fine given by
!   u_1(xi) = u_fine(x_1,y_1) + [u_fine(x_2,y_2)-u_fine(x_1,y_1)]xi
! where points 1 and 2 are the endpoints of the subinterval if p(2)/=0
! and the whole interval if p(2)==0.
! But we only need du_1/dxi which is the expression in square brackets


      if (p(2)==0) then
         du1dxi = ref_soln%vertex_solution(endpt(3),1,1) - &
                  ref_soln%vertex_solution(endpt(1),1,1)
      else
         if (loop==1) then
            du1dxi = ref_soln%vertex_solution(endpt(2),1,1) - &
                     ref_soln%vertex_solution(endpt(1),1,1)
         else
            du1dxi = ref_soln%vertex_solution(endpt(3),1,1) - &
                     ref_soln%vertex_solution(endpt(2),1,1)
         endif
      endif

! matrix entries are given by
! a_ij = integral dphi_i/dxi dphi_j/dxi
! right hand side is given by
! b_i = integral [(x_2-x_1)du_fine/dx + (y_2-y_1)du_fine/dy - du1dxi] dphi_i/dxi
! because du_fine/dxi = du_fine/dx dx/dxi + du_fine/dy dy/dxi
! and dx/dxi = x_2-x_1, dy/dxi = y_2-y_1

! get the quadrature rules.

! For the whole interval (p(2)==0), need to do the quadrature of each half to
! get the exact integral, so make a double length quadrature rule by
! concatinating the rules for each interval.  For the integral w.r.t. xi, get
! the rule for the unit interval and double it up into each half of the interval

      if (p(2)==0) then
         xline(1) = ref_soln%vertex(endpt(1))%coord%x
         xline(2) = ref_soln%vertex(endpt(2))%coord%x
         yline(1) = ref_soln%vertex(endpt(1))%coord%y
         yline(2) = ref_soln%vertex(endpt(2))%coord%y
         call quadrature_rule_line(qorder,xline,yline,nqpoints,qweightshalf, &
                                   xquadhalf,yquadhalf,jerr)
         allocate(qweights(2*nqpoints),xquadu(2*nqpoints),yquadu(2*nqpoints), &
                  xquad(2*nqpoints),yquad(2*nqpoints))
         xquadu(1:nqpoints) = xquadhalf
         yquadu(1:nqpoints) = yquadhalf
         deallocate(qweightshalf,xquadhalf,yquadhalf)
         xline(1) = ref_soln%vertex(endpt(2))%coord%x
         xline(2) = ref_soln%vertex(endpt(3))%coord%x
         yline(1) = ref_soln%vertex(endpt(2))%coord%y
         yline(2) = ref_soln%vertex(endpt(3))%coord%y
         call quadrature_rule_line(qorder,xline,yline,nqpoints,qweightshalf, &
                                   xquadhalf,yquadhalf,jerr)
         xquadu(nqpoints+1:2*nqpoints) = xquadhalf
         yquadu(nqpoints+1:2*nqpoints) = yquadhalf
         deallocate(qweightshalf,xquadhalf,yquadhalf)
         call quadrature_rule_line(qorder,xtri(1:2),ytri(1:2),nqpoints, &
                                   qweightshalf,xquadhalf,yquadhalf,jerr)
         xquad(1:nqpoints) = xquadhalf/2
         yquad(1:nqpoints) = 0
         xquad(nqpoints+1:2*nqpoints) = xquadhalf/2 + 0.5_my_real
         yquad(nqpoints+1:2*nqpoints) = 0
         qweights(1:nqpoints) = qweightshalf
         qweights(nqpoints+1:2*nqpoints) = qweightshalf
! since the quadrature weights include the interval length, and now each half
! of the quadrature points form an integral over half the unit interval, the
! weights should be divided by 2
! TEMP but I get strong p refinement at the singularity with this, whereas I
!      get a pretty good strategy without it.  maybe I'm wrong about halving it
!      On the other hand, by dividing by 2 the integral is very close to the
!      integral by using a single quadrature rule over the whole interval
!         qweights = qweights/2
         deallocate(qweightshalf,xquadhalf,yquadhalf)
         nqpoints = 2*nqpoints
         deltax = ref_soln%vertex(endpt(3))%coord%x - &
                  ref_soln%vertex(endpt(1))%coord%x
         deltay = ref_soln%vertex(endpt(3))%coord%y - &
                  ref_soln%vertex(endpt(1))%coord%y

! for the half intervals, get the quadrature rule for this interval

      else
         if (loop==1) then
            xline(1) = ref_soln%vertex(endpt(1))%coord%x
            xline(2) = ref_soln%vertex(endpt(2))%coord%x
            yline(1) = ref_soln%vertex(endpt(1))%coord%y
            yline(2) = ref_soln%vertex(endpt(2))%coord%y
            call quadrature_rule_line(qorder,xline,yline, &
                                      nqpoints,qweights,xquadu,yquadu,jerr)
            deallocate(qweights)
            call quadrature_rule_line(qorder,xtri(1:2),ytri(1:2),nqpoints, &
                                      qweights,xquad,yquad,jerr)
         else
            xline(1) = ref_soln%vertex(endpt(2))%coord%x
            xline(2) = ref_soln%vertex(endpt(3))%coord%x
            yline(1) = ref_soln%vertex(endpt(2))%coord%y
            yline(2) = ref_soln%vertex(endpt(3))%coord%y
            call quadrature_rule_line(qorder,xline,yline, &
                                      nqpoints,qweights,xquadu,yquadu,jerr)
            deallocate(qweights)
            call quadrature_rule_line(qorder,xtri(1:2),ytri(1:2),nqpoints, &
                                      qweights,xquad,yquad,jerr)
         endif
         deltax = xline(2) - xline(1)
         deltay = yline(2) - yline(1)
      endif

! evaluate the derivatives of the bubble basis functions at the quadrature points

      if (n > 0) then
         allocate(dphidxi(p(loop)+2,nqpoints))
         call p_hier_basis_func(xquad,yquad,xtri,ytri,(/0,0,p(loop),0/),"a", &
                                basisx=dphidxi)
         dphidxi(1:n,:) = dphidxi(4:p(loop)+2,:)
      endif

! evaluate the derivatives of u_fine at the quadrature points

      allocate(dudx(1,1,nqpoints),dudy(1,1,nqpoints))

! If there is a single interval (p(2)==0), then the quadrature points are
! split among the two children of edge.  Those for which the quadrature
! points on the unit interval are greater than 1/2 are in the second child.
! If there are two intervals, the quadrature points are all in the first or
! second child.

      if (p(2)==0) then
         halfqp = nqpoints/2
         call evaluate_soln_local(ref_soln,xquadu(1:halfqp),yquadu(1:halfqp), &
                                  child1,(/1/),(/1/), &
                                  ux=dudx(:,:,1:halfqp),uy=dudy(:,:,1:halfqp))
         call evaluate_soln_local(ref_soln,xquadu(halfqp+1:),yquadu(halfqp+1:),&
                                 child2,(/1/),(/1/), &
                                 ux=dudx(:,:,halfqp+1:),uy=dudy(:,:,halfqp+1:))
      else
         if (loop==1) then
            call evaluate_soln_local(ref_soln,xquadu,yquadu,child1,(/1/),(/1/),&
                                     ux=dudx,uy=dudy)
         else
            call evaluate_soln_local(ref_soln,xquadu,yquadu,child2,(/1/),(/1/),&
                                     ux=dudx,uy=dudy)
         endif
      endif

! compute the integrals

      a = 0
      b = 0
      do qp=1,nqpoints
         do i=1,n
            do j=1,n
               a(i,j) = a(i,j) + qweights(qp)*dphidxi(i,qp)*dphidxi(j,qp)
            end do
            b(i,1) = b(i,1) + qweights(qp)* &
                  (deltax*dudx(1,1,qp)+deltay*dudy(1,1,qp)-du1dxi)*dphidxi(i,qp)
         end do
      end do

! solve the linear system

      if (n > 0) then
         if (my_real == kind(1.0)) then
            call sgesv(n,1,a,nmax,ipiv,b,nmax,jerr)
         elseif (my_real == kind(1.0d0)) then
            call dgesv(n,1,a,nmax,ipiv,b,nmax,jerr)
         else
            ierr = PHAML_INTERNAL_ERROR
            call fatal("in weighted_H1_seminorm_diff_squared, LAPACK requires single or double precision")
            stop
         endif

         if (jerr /= 0) then
            ierr = PHAML_INTERNAL_ERROR
            call fatal("in weighted_H1_seminorm_diff_squared, dgesv failed")
            stop
         endif
      endif

! compute the integral of (du_fine/dxi - du_1/dxi - du_2/dxi)**2

      do qp=1,nqpoints
         integrand = deltax*dudx(1,1,qp)+deltay*dudy(1,1,qp)-du1dxi
         do i=1,n
            integrand = integrand - b(i,1)*dphidxi(i,qp)
         end do
         diff(loop) = diff(loop) + qweights(qp)*integrand**2
      end do

! deallocate memory allocated in this loop

      deallocate(qweights,xquad,yquad,xquadu,yquadu,dudx,dudy)
      if (n > 0) deallocate(dphidxi)

! next subinterval

   end do

! free workspace

   deallocate(a,b,ipiv)

! compute total integral

   diff(3) = diff(1) + diff(2)

end subroutine weighted_H1_seminorm_diff_squared

!          ------------------------------
subroutine compute_guaranteed_edge_rate(ref_soln,old_grid,edge,elem,best_h, &
                                        diff_old,guaranteed_edge_rate, &
                                        edge_rate_max,new_degree)
!          ------------------------------

!----------------------------------------------------
! This routine computes the guaranteed edge rate for edge.
! If edge_rate_max and new_degree are present, it, instead, computes the new
! degrees for an h-refined edge by quitting at 1/3 edge_rate_max.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: ref_soln, old_grid
integer, intent(in) :: edge, elem, best_h
real(my_real), intent(in) :: diff_old
real(my_real), intent(out) :: guaranteed_edge_rate
real(my_real), intent(in), optional :: edge_rate_max
integer, intent(out), optional :: new_degree(2)
!----------------------------------------------------
! Local variables:

! TEMP081203 try best_degree and best_rate
integer :: p, degree(2), best_degree(2)
real(my_real) :: diff_new(3), edge_rate, best_rate
!----------------------------------------------------
! Begin executable code

! edge_rate_max and new_degree must both be present or absent

   if (present(edge_rate_max) .neqv. present(new_degree)) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("in compute_guaranteed_edge_rate, edge_rate_max and new_degree must both be present or both absent")
      stop
   endif

! initialize result

   guaranteed_edge_rate = 0.0_my_real

! begin with the best h refinement

   p = old_grid%edge(edge)%degree
   degree(1) = best_h
   degree(2) = p+1 - best_h
   best_degree = degree ! TEMP081203
   best_rate = -huge(0.0_my_real) ! TEMP081203

! set new_degree to the best h refinement to return that if the edge rate
! condition is never met
! TEMP I should go back to the book and justify that to myself

   if (present(new_degree)) new_degree = degree

! do until both child edges have degree p+1

   do

! compute |u_fine - w_new|^2 where |.| is the weighted H^1 seminorm and
! w_new is the projection based interpolant on the current h refinement of edge

      call weighted_H1_seminorm_diff_squared(ref_soln,old_grid,edge,elem, &
                                             degree,diff_new)

!                       |u_fine-w_{hp}|^2 - |u_fine-w_new|^2
! compute edge_rate =  --------------------------------------
!                              (p_1+p_2-1)-(p-1)

      edge_rate = (diff_old-diff_new(3)) / (degree(1)+degree(2)-p)
      guaranteed_edge_rate = max(edge_rate, guaranteed_edge_rate)

! TEMP081203
! if this is the best rate we have seen, keep it

      if (edge_rate > best_rate) then
         best_rate = edge_rate
         best_degree = degree
      endif
! end TEMP081203

! if edge_rate_max is present, we are computing the new degrees

      if (present(edge_rate_max)) then
         if (edge_rate < edge_rate_max/3) then
            new_degree = degree
            new_degree = best_degree ! TEMP081203
            exit
         endif
      endif

! p refine any child with p_i < p+1 and rate > 70% of maximum rate

      if (degree(1) >= p+1 .and. degree(2) >= p+1) exit
      if (degree(1) >= p+1) then
         degree(2) = degree(2) + 1
      elseif (degree(2) >= p+1) then
         degree(1) = degree(1) + 1
      else
         if (degree(1) < p+1 .and. &
             diff_new(1) >= 0.7_my_real*max(diff_new(1),diff_new(2))) then
            degree(1) = degree(1) + 1
         endif
         if (degree(2) < p+1 .and. &
             diff_new(2) >= 0.7_my_real*max(diff_new(1),diff_new(2))) then
            degree(2) = degree(2) + 1
         endif
      endif

! next noncompetative h-refinement

   end do

end subroutine compute_guaranteed_edge_rate

!          -------------------------
subroutine determine_edges_to_refine(reftype,guaranteed_edge_rate,edge_rate_max)
!          -------------------------

!----------------------------------------------------
! This routine determines which edges should be refined
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

character(len=1), intent(inout) :: reftype(:)
real(my_real), intent(in) :: guaranteed_edge_rate(:)
real(my_real), intent(in) :: edge_rate_max
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

! mark edges for which the guaranteed edge rate is too small as not to refine

   where (guaranteed_edge_rate < edge_rate_max/3) reftype = "n"

end subroutine determine_edges_to_refine

!          ---------------------------
subroutine unrefine_edge_h_refinements(grid,old_grid,edge_reftype, &
                                       refine_control)
!          ---------------------------

!----------------------------------------------------
! This routine un-h-refines elements in the fine grid for which the base
! of the parent in the coarse grid is not marked for h-refinement, and
! remarks edges that need to be h-refined for compatibility
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(grid_type), intent(in) :: old_grid
character(len=1), intent(inout) :: edge_reftype(:)
type(refine_options), intent(in) :: refine_control
!----------------------------------------------------
! Local variables:

integer :: lev, elem, errcode
!----------------------------------------------------
! Begin executable code

! for each leaf element in old_grid (in reverse order of levels)

   do lev=old_grid%nlev,1,-1
      elem = old_grid%head_level_elem(lev)
      do while (elem /= END_OF_LIST)
         if (old_grid%element(elem)%isleaf) then

! if the base of this element is not marked for h-refinement, then unrefine
! it in the fine grid

            if (edge_reftype(old_grid%element(elem)%edge(3)) /= "h") then
               call unbisect_triangle_pair(grid,elem,errcode,refine_control)

! if unrefinement returned error code -1, there are grandchildren of elem so
! the base will actually have to be h-refined

               if (errcode == -1) then
                  edge_reftype(old_grid%element(elem)%edge(3)) = "h"
               endif
            endif

! next element

         endif
         elem = old_grid%element(elem)%next
      end do
   end do

end subroutine unrefine_edge_h_refinements

!          --------------------------------------
subroutine det_and_set_opt_deg_edges(grid,old_grid,ref_soln,refine_control, &
                                     edge_reftype,best_h,edge_rate_max,diff_old)
!          --------------------------------------

!----------------------------------------------------
! This routine determines the optimal edge degrees and sets them.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(grid_type), intent(in) :: old_grid, ref_soln
type(refine_options), intent(in) :: refine_control
character(len=1), intent(in) :: edge_reftype(:)
integer, intent(in) :: best_h(:)
real(my_real), intent(in) ::  edge_rate_max
real(my_real), intent(in) :: diff_old(:)
!----------------------------------------------------
! Local variables:

logical :: visited(size(old_grid%edge))
integer :: lev, elem, side, edge, p, new_degree(2), children(2), mate
real(my_real) :: guaranteed_edge_rate
!----------------------------------------------------
! Begin executable code

! for each edge in old_grid

   visited = .false.
   do lev=1,old_grid%nlev
      elem = old_grid%head_level_elem(lev)
      do while (elem /= END_OF_LIST)
         if (old_grid%element(elem)%isleaf) then
            do side=1,3
               edge = old_grid%element(elem)%edge(side)
               if (.not. visited(edge)) then
                  visited(edge) = .true.
                  p = old_grid%edge(edge)%degree

! if the edge is not to be refined, it keeps the same degree

                  if (edge_reftype(edge) == "n") then
                     call set_edge_degree(grid,refine_control,edge,p)

! if the edge is to be p-refined, increase p by one

                  elseif (edge_reftype(edge) == "p") then
                     call set_edge_degree(grid,refine_control,edge,p+1)

! if the edge is h-refined, determine the degrees of the child edges in the
! same way the guaranteed edge rate was computed, except quitting when the
! edge rate drops below 1/3 the maximum edge rate, and set the degree for the
! other two edges formed by bisecting two triangles to 1, to be set by the
! minimum rule later
! this only applies to the base of the current element

                  elseif (side==3) then ! (edge_reftype(edge) == "h" and base)
                     call compute_guaranteed_edge_rate(ref_soln,old_grid,edge, &
                                                       elem,best_h(edge), &
                                                       diff_old(edge), &
                                                       guaranteed_edge_rate, &
                                                       edge_rate_max,new_degree)
                     children = get_child_lid(grid%element(elem)%gid, &
                                              ALL_CHILDREN,grid%elem_hash)
                     call set_edge_degree(grid,refine_control, &
                                grid%element(children(1))%edge(2),new_degree(1))
                     call set_edge_degree(grid,refine_control, &
                                grid%element(children(2))%edge(2),new_degree(2))
                     call set_edge_degree(grid,refine_control, &
                                grid%element(children(1))%edge(1),1)
                     if (.not.(grid%element(elem)%mate == BOUNDARY)) then
                        mate = hash_decode_key(grid%element(elem)%mate, &
                                               grid%elem_hash)
                        children = get_child_lid(grid%element(mate)%gid, &
                                                 ALL_CHILDREN,grid%elem_hash)
                        call set_edge_degree(grid,refine_control, &
                                            grid%element(children(1))%edge(1),1)
                     endif

! remaining case of h refinement but not the base is handled by the element
! for which it is a base, so it has not been visited yet

                  else
                     visited(edge) = .false.

                  endif

! next edge
               endif
            end do
         endif
         elem = old_grid%element(elem)%next
      end do
   end do

end subroutine det_and_set_opt_deg_edges

!          ------------------------
subroutine det_and_set_opt_deg_elem(grid,old_grid,ref_soln,refine_control, &
                                    edge_rate_max)
!          ------------------------

!----------------------------------------------------
! This is a high level routine to determine and set the degree of the elements
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(grid_type), intent(in) :: old_grid, ref_soln
type(refine_options), intent(in) :: refine_control
real(my_real), intent(in) :: edge_rate_max
!----------------------------------------------------
! Local variables:

real(my_real) :: guaranteed_element_rate(size(old_grid%element))
real(my_real) :: element_rate_max, rate_max
integer :: lev, elem, ndescendants
integer :: descendants(4), new_degree(4)
!----------------------------------------------------
! Begin executable code

! initialize the guaranteed element rate at a number smaller than will be computed

   guaranteed_element_rate = -huge(0.0_my_real)

! for each leaf element in old_grid

   do lev=1,old_grid%nlev
      elem = old_grid%head_level_elem(lev)
      do while (elem /= END_OF_LIST)
         if (old_grid%element(elem)%isleaf) then

! determine the set of leaf descendants of elem in the h-refined grid.  This
! could be the element itself, or the two children, or one child and two
! grandchildren, or four grandchildren

            call get_leaf_descendants(grid,elem,descendants,ndescendants)

! initialize the degrees of the descendants according to the minimum rule

            call initialize_element_degrees(grid,refine_control,descendants, &
                                            ndescendants)

! compute the guaranteed element rate for elem

            call compute_guaranteed_element_rate(grid,old_grid,ref_soln,elem, &
                                                 descendants,ndescendants, &
                                                 guaranteed_element_rate(elem))

! next element

         endif
         elem = old_grid%element(elem)%next
      end do
   end do

! determine the maximum rate

   element_rate_max = maxval(guaranteed_element_rate)
   rate_max = max(edge_rate_max,element_rate_max)

! for each leaf element in old_grid

   do lev=1,old_grid%nlev
      elem = old_grid%head_level_elem(lev)
      do while (elem /= END_OF_LIST)
         if (old_grid%element(elem)%isleaf) then

! determine and set the degree for the descendants of this element

            call get_leaf_descendants(grid,elem,descendants,ndescendants)
            call compute_guaranteed_element_rate(grid,old_grid,ref_soln,elem, &
                                                descendants,ndescendants, &
                                                guaranteed_element_rate(elem), &
                                                rate_max,new_degree)
            call set_element_degree(grid,refine_control,descendants, &
                                    ndescendants,new_degree)

! next element

         endif
         elem = old_grid%element(elem)%next
      end do
   end do

end subroutine det_and_set_opt_deg_elem

!          --------------------
subroutine get_leaf_descendants(grid,elem,descendants,ndescendants)
!          --------------------

!----------------------------------------------------
! This routine returns a list of the leaf descendants of elem where elem is
! a leaf in old_grid.  Consequently, the descendants go at most two
! generations.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: elem
integer, intent(out) :: descendants(4), ndescendants
!----------------------------------------------------
! Local variables:

integer :: children(2)
!----------------------------------------------------
! Begin executable code

! if the element is a leaf, it is the only descendant

   if (grid%element(elem)%isleaf) then
      descendants(1) = elem
      ndescendants = 1
   else

! get the children of elem

      children = get_child_lid(grid%element(elem)%gid,ALL_CHILDREN, &
                               grid%elem_hash)

! if the first child is a leaf, it is the first descendant

      if (grid%element(children(1))%isleaf) then
         descendants(1) = children(1)
         ndescendants = 1

! otherwise its children are the first two descendants

      else
         descendants(1:2) = get_child_lid(grid%element(children(1))%gid, &
                                          ALL_CHILDREN,grid%elem_hash)
         ndescendants = 2
      endif

! if the second child is a leaf, it is the next descendant

      if (grid%element(children(2))%isleaf) then
         descendants(ndescendants+1) = children(2)
         ndescendants = ndescendants+1

! otherwise its children are the remaining descendants

      else
         descendants(ndescendants+1:ndescendants+2) = &
            get_child_lid(grid%element(children(2))%gid,ALL_CHILDREN, &
                          grid%elem_hash)
         ndescendants = ndescendants + 2
      endif
   endif

end subroutine get_leaf_descendants

!          --------------------------
subroutine initialize_element_degrees(grid,refine_control,descendants, &
                                      ndescendants)
!          --------------------------

!----------------------------------------------------
! This routine sets the degree of the elements in descendants to
! satisfy the minimum rule.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(refine_options), intent(in) :: refine_control
integer, intent(in) :: descendants(:), ndescendants
!----------------------------------------------------
! Local variables:

integer :: i, elem, olddeg, newdeg, oldsize, newsize, d1, d2, d3, astat
real(my_real), pointer :: temp1(:,:,:)
!----------------------------------------------------
! Begin executable code

! for each decendant

   do i=1,ndescendants
      elem = descendants(i)

! degree(descendant) = maximum of descendant edge degrees

      olddeg = grid%element(elem)%degree
      newdeg = maxval(grid%edge(grid%element(elem)%edge)%degree)
      grid%element(elem)%degree = newdeg

! adjust degrees of freedom count

      if (olddeg >= 3) grid%dof = grid%dof - ((olddeg-2)*(olddeg-1))/2
      if (newdeg >= 3) grid%dof = grid%dof + ((newdeg-2)*(newdeg-1))/2
      if (grid%element(elem)%iown) then
         if (olddeg >= 3) grid%dof_own=grid%dof_own - ((olddeg-2)*(olddeg-1))/2
         if (newdeg >= 3) grid%dof_own=grid%dof_own + ((newdeg-2)*(newdeg-1))/2
      endif

! make sure allocated memory is large enough

      newsize = ((newdeg-1)*(newdeg-2))/2
      if (associated(grid%element(elem)%solution)) then
         oldsize = size(grid%element(elem)%solution,dim=1)
      else
         oldsize = 0
      endif
      if (oldsize < newsize) then
         allocate(temp1(newsize,grid%system_size,max(1,grid%num_eval)))
         temp1 = 0.0_my_real
         if (oldsize > 0) temp1(1:oldsize,:,:) = grid%element(elem)%solution
         deallocate(grid%element(elem)%solution, stat=astat)
         grid%element(elem)%solution => temp1
         if (grid%have_true) then
            nullify(temp1)
            allocate(temp1(newsize,grid%system_size,max(1,grid%num_eval)))
            temp1 = 0.0_my_real
            if (oldsize > 0) temp1(1:oldsize,:,:) = grid%element(elem)%exact
            deallocate(grid%element(elem)%exact, stat=astat)
            grid%element(elem)%exact => temp1
         endif
         if (grid%oldsoln_exists) then
            if (associated(grid%element(elem)%oldsoln)) then
               nullify(temp1)
               allocate(temp1(newsize,grid%system_size,max(1,grid%num_eval)))
               temp1 = 0.0_my_real
               d1=min(size(grid%element(elem)%oldsoln,dim=1),size(temp1,dim=1))
               d2=min(size(grid%element(elem)%oldsoln,dim=2),size(temp1,dim=2))
               d3=min(size(grid%element(elem)%oldsoln,dim=3),size(temp1,dim=3))
               temp1(1:d1,1:d2,1:d3)=grid%element(elem)%oldsoln(1:d1,1:d2,1:d3)
               deallocate(grid%element(elem)%oldsoln)
               grid%element(elem)%oldsoln => temp1
            endif
         endif
      endif

   end do

end subroutine initialize_element_degrees

!          -------------------------------
subroutine compute_guaranteed_element_rate(grid,old_grid,ref_soln,elem, &
                                           descendants,ndescendants, &
                                           guaranteed_element_rate, &
                                           rate_max,new_degree)
!          -------------------------------

!----------------------------------------------------
! This routine computes the guaranteed element rate for elem.
! If rate_max and new_degree are present, it, instead, computes the new
! degrees for the descendants by quitting at 1/3 rate_max.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid, old_grid, ref_soln
integer, intent(in) :: elem, descendants(:), ndescendants
real(my_real), intent(out) :: guaranteed_element_rate
real(my_real), intent(in), optional :: rate_max
integer, intent(out), optional :: new_degree(:)
!----------------------------------------------------
! Local variables:

! TEMP081203 try best_degree and best_rate
integer :: p, degree(ndescendants), i, N_new, N_old, best_degree(ndescendants)
real(my_real) :: diff_hp(2), diff_new(ndescendants+1), element_rate, best_rate
!----------------------------------------------------
! Begin executable code

! rate_max and new_degree must both be present or absent

   if (present(rate_max) .neqv. present(new_degree)) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("in compute_guaranteed_element_rate, rate_max and new_degree must both be present or both absent")
      stop
   endif

! initialize result

   guaranteed_element_rate = 0.0_my_real

! begin with degrees set in initialize_element_degrees

   p = old_grid%element(elem)%degree
   do i = 1,ndescendants
      degree(i) = grid%element(descendants(i))%degree
   end do
   best_degree = degree ! TEMP081203
   best_rate = -huge(0.0_my_real) ! TEMP081203

! put that in new_degree also, in case the rate condition is never met
! TEMP justify?

   if (present(new_degree)) new_degree(1:ndescendants) = degree

! compute |u_fine - w_{hp}|^2 where |.| is the H^1 seminorm over elem and
! w_{hp} is the projection based interpolant on elem on the old grid

   call H1_elem_seminorm_diff_squared(ref_soln,grid,(/elem/),1,(/p/),diff_hp)
   N_old = ((p-1)*(p-2))/2

! do until all descendants have degree >= p+1

   do

! compute |u_fine - w_new|^2 where |.| is the H^1 seminorm over elem and
! w_new is the projection based interpolant on the descendants of elem with
! the current degrees

      call H1_elem_seminorm_diff_squared(ref_soln,grid,descendants, &
                                         ndescendants,degree,diff_new)

!                          |u_fine-w_{hp}|^2 - |u_fine-w_new|^2
! compute element_rate =  --------------------------------------
!                                    N_new - N_hp

! N_new is the number of interior degrees of freedom, but also make sure it
! is greater than N_old

      N_new = 0
      do i=1,ndescendants
         N_new = N_new + ((degree(i)-1)*(degree(i)-2))/2
      end do
      N_new = max(N_old+1,N_new)
      element_rate = (diff_hp(2)-diff_new(ndescendants+1)) / (N_new - N_old)
      guaranteed_element_rate = max(element_rate, guaranteed_element_rate)

! TEMP081203
! if this is the best rate we have seen, keep it

      if (element_rate > best_rate) then
         best_rate = element_rate
         best_degree = degree
      endif
! end TEMP081203
! if rate_max is present, we are computing the new degrees

      if (present(rate_max)) then
         if (element_rate >= 0.0_my_real .and. element_rate < rate_max/3) then
            new_degree(1:ndescendants) = degree
            new_degree(1:ndescendants) = best_degree ! TEMP081203
            exit
         endif
      endif

! p refine any descendant with p_i < p+1 and rate > 70% of maximum rate
! but don't include those that are already at p+1

      if (all(degree >= p+1)) exit
      where (degree >= p+1) diff_new(1:ndescendants) = 0.0_my_real
      do i=1,ndescendants
         if (degree(i) < p+1 .and. &
             diff_new(i) >= 0.7_my_real*maxval(diff_new(1:ndescendants))) then
            degree(i) = degree(i) + 1
         endif
      end do

! next assignment

   end do

end subroutine compute_guaranteed_element_rate

!          -----------------------------
subroutine H1_elem_seminorm_diff_squared(ref_soln,grid,descendants, &
                                         ndescendants,p,diff)
!          -----------------------------

!----------------------------------------------------
! This routine computes |u_fine - w|^2 where |.| is the H1 seminorm on the
! elements listed in descendants and w is the projection based interpolation
! of the reference solution onto the elements with degrees specified in p.
!
! On return, diff(1:ndescendants) contain the contribution of each element,
! and diff(ndescendants+1) contains the sum of them.
!
!
! On each descendant, the interpolant w = u_1 + u_2 + u_3,
! with u_2 = u_2,1 + u_2,2 + u_2,3, is defined by solving four variational
! problems.
! u_1 is the linear interpolant of the reference solution at the three vertices.
! For each edge e=1,2,3 we solve for the edge interpolant exactly as in
! subroutine weighted_H1_seminorm_diff_squared.
! Finally, the interior bubble is given by solving
!   u_3 in V(K)
!   <u_3,v>_K = <(u_fine-u_1-u_2,v>_K  forall v in V(K)
! where V(K) is the space of bubble functions over element K up to degree p(K)
! and the inner product is the H^1 seminorm over element K.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: ref_soln, grid
integer, intent(in) :: ndescendants, descendants(ndescendants), p(ndescendants)
real(my_real), intent(out) :: diff(ndescendants+1)
!----------------------------------------------------
! Local variables:

integer :: descend, nbasis, side, degree(4), endpt(2), elem, n_sofar, edge, &
           n, qorder, nqpoints, qp, i, j, jerr, have_it, container
integer, allocatable :: ipiv(:)
real(my_real) :: xtri(3), ytri(3), xline(2), yline(2), deltax, deltay, du1dxi, &
                 xvert(3), yvert(3), integrandx, integrandy
real(my_real), allocatable :: interp(:), a(:,:), b(:,:), dphidxi(:,:), &
                              dudx(:,:,:), dudy(:,:,:), phi(:,:), dphidx(:,:), &
                              dphidy(:,:), du1u2dx(:), du1u2dy(:)
real(my_real), pointer :: qweights(:), xquad(:), yquad(:), xquadu(:), yquadu(:)
!----------------------------------------------------
! Begin executable code

! (xtri,ytri) is a triangle containing the unit interval on the first edge, used
! for evaluating edge basis functions as a function of (Greek letter) xi

   xtri(1) = 0.0_my_real
   ytri(1) = 0.0_my_real
   xtri(2) = 1.0_my_real
   ytri(2) = 0.0_my_real
   xtri(3) = 0.5_my_real
   ytri(3) = 1.0_my_real

! For each descendant

   do descend = 1,ndescendants
      elem = descendants(descend)

! identify the degree of the edges and interior of this element

      degree(1:3) = grid%edge(grid%element(elem)%edge)%degree
      degree(4) = p(descend)

! interp will contain the coefficients of the interpolant, in the same order
! that the basis functions are returned

      nbasis = 3
      do side=1,3
         nbasis = nbasis + degree(side)-1
      end do
      nbasis = nbasis + ((degree(4)-1)*(degree(4)-2))/2
      allocate(interp(nbasis))

! the linear interpolant is just the solution at the vertices

      interp(1:3) = ref_soln%vertex_solution(ref_soln%element(elem)%vertex,1,1)
      n_sofar = 3

! for each edge

      do side=1,3
         edge = ref_soln%element(elem)%edge(side)
         endpt = ref_soln%edge(edge)%vertex

! allocate linear system for the coefficients of these edge basis functions

         n = degree(side)-1
         allocate(a(n,n),b(n,1),ipiv(n))

! Endpoints of the interval over which to integrate.

         xline = ref_soln%vertex(endpt)%coord%x
         yline = ref_soln%vertex(endpt)%coord%y
         deltax = xline(2) - xline(1)
         deltay = yline(2) - yline(1)

! u_1 is the linear vertex interpolant of u_fine given by
!   u_1(xi) = u_fine(x_1,y_1) + [u_fine(x_2,y_2)-u_fine(x_1,y_1)]xi
! but we only need du_1/dxi which is the expression in square brackets

         du1dxi = ref_soln%vertex_solution(endpt(2),1,1) - &
                  ref_soln%vertex_solution(endpt(1),1,1)

! matrix entries are given by
! a_ij = integral dphi_i/dxi dphi_j/dxi
! right hand side is given by
! b_i = integral [(x_2-x_1)du_fine/dx + (y_2-y_1)du_fine/dy - du1dxi] dphi_i/dxi
! because du_fine/dxi = du_fine/dx dx/dxi + du_fine/dy dy/dxi
! and dx/dxi = x_2-x_1, dy/dxi = y_2-y_1

! quadratures for the matrix involve polynomials of degree 2*(degree(side)-1),
! and quadratures for the right side involve polynomials of degree
! (1+degree(u_fine)-1)+(degree(side)-1).
! Pick a quadrature rule that will be exact for both of them

         qorder = max(2*(degree(side)-1), &
                      (ref_soln%edge(edge)%degree+1)+(degree(side)-1))
         qorder = min((qorder+2)/2,MAX_QUAD_ORDER_LINE)
! TEMP if edge is not a leaf in ref_soln, the integral will not be exact, but
!      will be real close with the maximum quadrature order
         qorder = MAX_QUAD_ORDER_LINE ! TEMP081125
         call quadrature_rule_line(qorder,xline,yline, &
                                   nqpoints,qweights,xquadu,yquadu,jerr)
         deallocate(qweights)
         call quadrature_rule_line(qorder,xtri(1:2),ytri(1:2),nqpoints, &
                                   qweights,xquad,yquad,jerr)

! evaluate the derivatives of the bubble basis functions at quadrature points

         if (n > 0) then
            allocate(dphidxi(degree(side)+2,nqpoints))
            call p_hier_basis_func(xquad,yquad,xtri,ytri, &
                                   (/0,0,degree(side),0/),"a",basisx=dphidxi)
            dphidxi(1:n,:) = dphidxi(4:degree(side)+2,:)
         endif

! evaluate the derivatives of u_fine at the quadrature points
! If elem is not a leaf in ref_soln, then the correct containing element must
! be found to pass to evaluate_soln_local.  Since we don't know the clustering
! of quadrature points in the descendants of elem, we have to do one quadrature
! point at a time.
! TEMP not so true for edges, I was thinking of quadrature over the elements.
!     This should still work but is less efficient than doing the two half edges

         allocate(dudx(1,1,nqpoints),dudy(1,1,nqpoints))
         if (ref_soln%element(elem)%isleaf) then
            call evaluate_soln_local(ref_soln,xquadu,yquadu,elem,(/1/),(/1/), &
                                     ux=dudx,uy=dudy)
         else
            do qp=1,nqpoints
               call find_containing_leaf(xquadu(qp),yquadu(qp),elem,have_it, &
                                         container,ref_soln)
! TEMP when I parallelize this, note that have_it=0 and elem=0 if I don't
!      own the element that contains the point.  for now, just check that its
!      not zero, which it shouldn't be if nproc==1
               if (have_it==0) then
                  ierr = PHAML_INTERNAL_ERROR
                  call fatal("couldn't find containing element in H1_elem_seminorm_diff_squared")
                  stop
               endif
               call evaluate_soln_local(ref_soln,xquadu(qp:qp),yquadu(qp:qp), &
                                        container,(/1/),(/1/), &
                                        ux=dudx(:,:,qp:qp),uy=dudy(:,:,qp:qp))
            enddo
         endif

! compute the integrals

         a = 0
         b = 0
         do qp=1,nqpoints
            do i=1,n
               do j=1,n
                  a(i,j) = a(i,j) + qweights(qp)*dphidxi(i,qp)*dphidxi(j,qp)
               end do
               b(i,1) = b(i,1) + qweights(qp)* &
                  (deltax*dudx(1,1,qp)+deltay*dudy(1,1,qp)-du1dxi)*dphidxi(i,qp)
            end do
         end do

! solve the linear system

         if (n > 0) then
            if (my_real == kind(1.0)) then
               call sgesv(n,1,a,n,ipiv,b,n,jerr)
            elseif (my_real == kind(1.0d0)) then
               call dgesv(n,1,a,n,ipiv,b,n,jerr)
            else
               ierr = PHAML_INTERNAL_ERROR
               call fatal("in H1_elem_seminorm_diff_squared, LAPACK requires single or double precision")
               stop
            endif

            if (jerr /= 0) then
               ierr = PHAML_INTERNAL_ERROR
               call fatal("in H1_elem_seminorm_diff_squared, dgesv failed")
               stop
            endif
         endif

! copy the results to interp

         interp(n_sofar+1:n_sofar+n) = b(:,1)
         n_sofar = n_sofar + n

! deallocate memory allocated in this loop

         deallocate(a,b,ipiv,qweights,xquad,yquad,xquadu,yquadu,dudx,dudy)
         if (n > 0) deallocate(dphidxi)

! next edge

      end do

! vertices of the element

      xvert = ref_soln%vertex(ref_soln%element(elem)%vertex)%coord%x
      yvert = ref_soln%vertex(ref_soln%element(elem)%vertex)%coord%y

! polynomial degrees involved in integrals are:
! matrix: 2*(degree(4)-1)
! rhs:    degree(4)+max(degree(u_fine),degree(1:3))
! norm:   2*max(degree(u_fine),degree(1:4))
! Pick a quadrature rule that is exact for the largest of these, which happens
! to be the norm, assuming the edge degrees are not larger than the interior

      qorder = 2*max(ref_soln%element(elem)%degree,maxval(degree))
      qorder = min(qorder,MAX_QUAD_ORDER_TRI)
      call quadrature_rule_tri(qorder,xvert,yvert,nqpoints,qweights, &
                               xquad,yquad,jerr,.true.)

! evaluate the derivatives of all basis functions at the quadrature points

      allocate(phi(nbasis,nqpoints),dphidx(nbasis,nqpoints), &
               dphidy(nbasis,nqpoints))
      call p_hier_basis_func(xquad,yquad,xvert,yvert,degree,"a", &
                             basis=phi,basisx=dphidx,basisy=dphidy)

! evaluate the derivatives of u_1+u_2 at the quadrature points

      allocate(du1u2dx(nqpoints),du1u2dy(nqpoints))
      du1u2dx = 0.0_my_real
      du1u2dy = 0.0_my_real
      do qp=1,nqpoints
         do i=1,n_sofar
            du1u2dx(qp) = du1u2dx(qp) + interp(i)*dphidx(i,qp)
            du1u2dy(qp) = du1u2dy(qp) + interp(i)*dphidy(i,qp)
         end do
      end do

! evaluate the derivatives of u_fine at the quadrature points
! If elem is not a leaf in ref_soln, then the correct containing element must
! be found to pass to evaluate_soln_local.  Since we don't know the clustering
! of quadrature points in the descendants of elem, we have to do one quadrature
! point at a time.

         allocate(dudx(1,1,nqpoints),dudy(1,1,nqpoints))
         if (ref_soln%element(elem)%isleaf) then
            call evaluate_soln_local(ref_soln,xquad,yquad,elem,(/1/),(/1/), &
                                     ux=dudx,uy=dudy)
         else
            do qp=1,nqpoints
               call find_containing_leaf(xquad(qp),yquad(qp),elem,have_it, &
                                         container,ref_soln)
! TEMP when I parallelize this, note that have_it=0 and elem=0 if I don't
!      own the element that contains the point.  for now, just check that its
!      not zero, which it shouldn't be if nproc==1
               if (have_it==0) then
                  ierr = PHAML_INTERNAL_ERROR
                  call fatal("couldn't find containing element in H1_elem_seminorm_diff_squared")
                  stop
               endif
               call evaluate_soln_local(ref_soln,xquad(qp:qp),yquad(qp:qp), &
                                        container,(/1/),(/1/), &
                                        ux=dudx(:,:,qp:qp),uy=dudy(:,:,qp:qp))
            enddo
         endif


! evaluate the integrals
! a_ij = integral dphi_i/dx*dphi_j/dx + dphi_i/dy*dphi_j/dy
! b_i = integral (d(u_fine-u_1-u_2)/dx + d(u_fine-u_1-u_2)/dy) * phi_i

      n = nbasis - n_sofar
      allocate(a(n,n),b(n,1),ipiv(n))
      a = 0
      b = 0
      do qp=1,nqpoints
         do i=1,n
            do j=1,n
               a(i,j) = a(i,j) + qweights(qp)* &
                                (dphidx(i+n_sofar,qp)*dphidx(j+n_sofar,qp) + &
                                 dphidy(i+n_sofar,qp)*dphidy(j+n_sofar,qp))
            end do
! TEMP081209 maybe it's supposed to be phi', not phi
!            b(i,1) = b(i,1) + qweights(qp)*(dudx(1,1,qp)-du1u2dx(qp) + &
!                                     dudy(1,1,qp)-du1u2dy(qp))*phi(i+n_sofar,qp)
            b(i,1) = b(i,1) + qweights(qp)*( &
                   (dudx(1,1,qp)-du1u2dx(qp))*dphidx(i+n_sofar,qp) + &
                   (dudy(1,1,qp)-du1u2dy(qp))*dphidy(i+n_sofar,qp))

         end do
      end do

! solve the linear system

      if (n > 0) then
         if (my_real == kind(1.0)) then
            call sgesv(n,1,a,n,ipiv,b,n,jerr)
         elseif (my_real == kind(1.0d0)) then
            call dgesv(n,1,a,n,ipiv,b,n,jerr)
         else
            ierr = PHAML_INTERNAL_ERROR
            call fatal("in H1_elem_seminorm_diff_squared, LAPACK requires single or double precision")
            stop
         endif

         if (jerr /= 0) then
            ierr = PHAML_INTERNAL_ERROR
            call fatal("in H1_elem_seminorm_diff_squared, dgesv failed")
            stop
         endif
      endif

! copy the results to interp

      interp(n_sofar+1:nbasis) = b(:,1)

! compute the integral of d(u_fine-w)dx**2 + d(u_fine-w)/dy**2

      diff(descend) = 0.0_my_real
      do qp = 1,nqpoints
         integrandx = dudx(1,1,qp)
         integrandy = dudy(1,1,qp)
         do i=1,nbasis
            integrandx = integrandx - interp(i)*dphidx(i,qp)
            integrandy = integrandy - interp(i)*dphidy(i,qp)
         end do
         diff(descend) = diff(descend) + &
                         qweights(qp)*(integrandx**2+integrandy**2)
      end do

! deallocate memory

      deallocate(a,b,ipiv,dudx,dudy,du1u2dx,du1u2dy,phi,dphidx,dphidy, &
                 qweights,xquad,yquad,interp)

! next descendant

   end do

! compute the total seminorm over all descendants

   diff(ndescendants+1) = sum(diff(1:ndescendants))

end subroutine H1_elem_seminorm_diff_squared

!          -------------------
subroutine impose_minimum_rule(grid,refine_control)
!          -------------------

!----------------------------------------------------
! This routine sets all edge degrees in accordance with the minimum rule
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(refine_options), intent(in) :: refine_control
!----------------------------------------------------
! Local variables:

logical :: visited(size(grid%edge))
integer :: neighbors(3), lev, elem, side, edge, newdeg
!----------------------------------------------------
! Begin executable code

! for each leaf edge

   visited = .false.
   do lev=1,grid%nlev
      elem = grid%head_level_elem(lev)
      do while (elem /= END_OF_LIST)
         if (grid%element(elem)%isleaf) then
            neighbors = get_neighbors(elem,grid)
            do side = 1,3
               edge = grid%element(elem)%edge(side)
               if (.not. visited(edge)) then
                  visited(edge) = .true.

! the new degree is the minimum of the degree of elem and the neighbor

                  newdeg = grid%element(elem)%degree
                  if (neighbors(side) /= BOUNDARY) then
                     newdeg = min(newdeg,grid%element(neighbors(side))%degree)
                  endif

! if the new degree is not the current degree, change it

                  if (newdeg /= grid%edge(edge)%degree) then
                     call set_edge_degree(grid,refine_control,edge,newdeg)
                  endif

! next edge

               endif
            end do
         endif
         elem = grid%element(elem)%next
      end do
   end do

end subroutine impose_minimum_rule

!          ---------------
subroutine set_edge_degree(grid,refine_control,edge,new_deg)
!          ---------------

!----------------------------------------------------
! This routine sets the degree of edge to new_deg
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(refine_options), intent(in) :: refine_control
integer, intent(in) :: edge, new_deg
!----------------------------------------------------
! Local variables:

integer :: olddeg, newdeg, edge_size, oldsize, d1, d2, d3, j, astat
real(my_real), pointer :: temp1(:,:,:)
!----------------------------------------------------
! Begin executable code

! Remember the old degree and set the new degree.  If they're the same
! then there is nothing to do.

   olddeg = grid%edge(edge)%degree
   newdeg = min(new_deg,refine_control%max_deg)
   if (olddeg == newdeg) return
   grid%edge(edge)%degree = newdeg

! make sure there is enough memory for solution etc.

   edge_size = newdeg-1
   if (associated(grid%edge(edge)%solution)) then
      oldsize = size(grid%edge(edge)%solution,dim=1)
   else
      oldsize = 0
   endif
   if (oldsize < edge_size) then
      edge_size = newdeg+1
      allocate(temp1(edge_size,grid%system_size,max(1,grid%num_eval)))
      temp1 = 0.0_my_real
      if (oldsize > 0) temp1(1:oldsize,:,:) = grid%edge(edge)%solution
      deallocate(grid%edge(edge)%solution,stat=astat)
      grid%edge(edge)%solution => temp1
      if (grid%have_true) then
         nullify(temp1)
         allocate(temp1(edge_size,grid%system_size,max(1,grid%num_eval)))
         temp1 = 0.0_my_real
         if (oldsize > 0) temp1(1:oldsize,:,:) = grid%edge(edge)%exact
         deallocate(grid%edge(edge)%exact,stat=astat)
         grid%edge(edge)%exact => temp1
      endif
   endif
   if (grid%oldsoln_exists) then
      if (associated(grid%edge(edge)%oldsoln)) then
         nullify(temp1)
         allocate(temp1(edge_size,grid%system_size,max(1,grid%num_eval)))
         temp1 = 0.0_my_real
         d1 = min(size(grid%edge(edge)%oldsoln,dim=1),size(temp1,dim=1))
         d2 = min(size(grid%edge(edge)%oldsoln,dim=2),size(temp1,dim=2))
         d3 = min(size(grid%edge(edge)%oldsoln,dim=3),size(temp1,dim=3))
         temp1(1:d1,1:d2,1:d3) = grid%edge(edge)%oldsoln(1:d1,1:d2,1:d3)
         deallocate(grid%edge(edge)%oldsoln,stat=astat)
         grid%edge(edge)%oldsoln => temp1
      endif
   endif

! adjust the degree of freedom count for grid

   grid%dof = grid%dof + newdeg - olddeg
   if (grid%element(grid%edge(edge)%assoc_elem)%iown) then
      grid%dof_own = grid%dof_own + newdeg - olddeg
   endif

! if it is a Dirichlet boundary edge, set the solution.  Otherwise, zero
! the new parts of the solution or no longer used parts of the solution.
! Also set the exact solution if true was given.

   if (grid%have_true .and. associated(grid%edge(edge)%exact)) then
      grid%edge(edge)%exact = 0.0_my_real
   endif
   do j=1,grid%system_size
      if (grid%edge_type(edge,j) == DIRICHLET) then
         call edge_exact(grid,edge,j,"d")
      elseif (newdeg > olddeg .and. associated(grid%edge(edge)%solution)) then
         grid%edge(edge)%solution(olddeg:newdeg-1,j,:) = 0.0_my_real
      elseif (newdeg < olddeg .and. associated(grid%edge(edge)%solution)) then
         grid%edge(edge)%solution(newdeg:olddeg-1,j,:) = 0.0_my_real
      endif
      if (grid%have_true) call edge_exact(grid,edge,j,"t")
   end do

end subroutine set_edge_degree

!          ------------------
subroutine set_element_degree(grid,refine_control,descendants,ndescendants, &
                              new_degree)
!          ------------------

!----------------------------------------------------
! This routine sets the degree of element descendants(i) to new_degree(i)
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(refine_options), intent(in) :: refine_control
integer, intent(in) :: descendants(:), ndescendants, new_degree(:)
!----------------------------------------------------
! Local variables:

integer :: i, jerr
!----------------------------------------------------
! Begin executable code

! just use the p_refine and p_coarsen routines to adjust the degree

   do i=1,ndescendants

      do while (grid%element(descendants(i))%degree > new_degree(i))
         call p_coarsen_elem(grid,descendants(i),jerr,refine_control)
      end do

      do while (grid%element(descendants(i))%degree < new_degree(i))
         call p_refine_elem(grid,descendants(i),refine_control)
      end do

   end do

end subroutine set_element_degree

!---------------------------------------------------------------------
!  ROUTINES FOR THE REFSOLN_ELEM HP-ADAPTIVE STRATEGY
!---------------------------------------------------------------------

!          -------------------
subroutine refine_refsoln_elem(grid,procs,refine_control,solver_control, &
                               io_control,still_sequential,init_nvert, &
                               init_nelem,init_dof,loop,balance_what,predictive)
!          -------------------

!----------------------------------------------------
! This routine is the top level refine routine for the REFSOLN_ELEM hp strategy
! TEMP only using first eigenvector
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

!----------------------------------------------------
! Local variables:

type(grid_type) :: ref_soln
type(refine_options) :: loc_refcont
type(io_options) :: loc_iocont
!----------------------------------------------------
! Begin executable code

! TEMP not parallel

   if (num_proc(procs) > 1) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("REFSOLN_ELEM not yet implemented in parallel")
      stop
   endif

! start timing the refinement process

   call reset_watch(prefine)
   call start_watch((/prefine,trefine/))

! Create the reference solution by copying the current grid, performing
! uniform h and p refinements, and solving on the fine grid.

   call copy_grid(grid,ref_soln)
   loc_refcont = refine_control
   loc_refcont%error_estimator = HIERARCHICAL_COEFFICIENT
   loc_refcont%reftype = H_UNIFORM
   loc_iocont = io_control
   loc_iocont%print_linsys_when = NEVER
   loc_iocont%print_error_when = NEVER
   call refine(ref_soln,procs,loc_refcont,solver_control,loc_iocont, &
               still_sequential,init_nvert,init_nelem,init_dof,loop, &
               balance_what,predictive,no_time=.true.)
   loc_refcont%reftype = P_UNIFORM
   call refine(ref_soln,procs,loc_refcont,solver_control,loc_iocont, &
               still_sequential,init_nvert,init_nelem,init_dof,loop, &
               balance_what,predictive,no_time=.true.)
   call solve(ref_soln,procs,loc_iocont,solver_control,still_sequential, &
              .false.,no_time=.true.)

! MASTER doesn't participate further; it was only here to monitor solve

if (my_proc(procs) == MASTER) then
   call deallocate_grid(ref_soln)
   grid%errind_up2date = .false.
   call stop_watch((/prefine,trefine/))
   return
endif

! Compute the error estimate.

   grid%refsoln_errest = compute_refsoln_errest_elem(ref_soln,grid)

! If the error estimate is large enough, determine and perform refinements

   if (grid%refsoln_errest > refine_control%term_energy_err) then
      call det_and_perf_opt_ref_elem(grid,ref_soln,refine_control)
   endif

! Destroy the reference solution

   call deallocate_grid(ref_soln)

! Error indicators have not been (and won't really be) set

   grid%errind_up2date = .false.

! stop timing the refinement process

   call stop_watch((/prefine,trefine/))

end subroutine refine_refsoln_elem

!        ---------------------------
function compute_refsoln_errest_elem(ref_soln,old_grid)
!        ---------------------------

!----------------------------------------------------
! This routine computes the global error estimate
!    ||u_{h,p} - u_{h/sqrt(2),p+1}||_H^1
! and sets the elemental error estimates in old_grid
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: ref_soln
type(grid_type), intent(inout) :: old_grid
real(my_real) :: compute_refsoln_errest_elem
!----------------------------------------------------
! Local variables:

integer :: lev, elem, nqpoints, jerr, parent, old_elem, grandparent, ss, i, j
real(my_real) :: contribution
real(my_real), pointer :: qweights(:), xquad(:), yquad(:)
real(my_real), allocatable :: ux_fine(:,:,:), uy_fine(:,:,:), ux_old(:,:,:), &
                              uy_old(:,:,:), u_fine(:,:,:), u_old(:,:,:)
type(hash_key) :: parent_gid, grandparent_gid
!----------------------------------------------------
! Begin executable code

   ss = old_grid%system_size

! initalize the norm of the difference and solution, and elemental error estimates

   compute_refsoln_errest_elem = 0
   old_grid%element_errind = 0.0_my_real

!  for each leaf element of the fine grid

   do lev=1,ref_soln%nlev
      elem = ref_soln%head_level_elem(lev)
      do while (elem /= END_OF_LIST)
         if (ref_soln%element(elem)%iown .and. ref_soln%element(elem)%isleaf) then

! determine a quadrature rule for p+1 (which is p in the refined grid)

            call quadrature_rule_tri(ref_soln%element(elem)%degree, &
                       ref_soln%vertex(ref_soln%element(elem)%vertex)%coord%x, &
                       ref_soln%vertex(ref_soln%element(elem)%vertex)%coord%y, &
                       nqpoints,qweights,xquad,yquad,jerr,.true.)

! identify the element in old_grid that contains this element in the fine grid,
! as the grandparent if it is a leaf and parent otherwise

            parent_gid = ref_soln%element(elem)%gid/2
            parent = hash_decode_key(parent_gid,ref_soln%elem_hash)
            if (ref_soln%element(parent)%level == 1) then
               old_elem = parent
            else
               grandparent_gid = ref_soln%element(parent)%gid/2
               grandparent = hash_decode_key(grandparent_gid,ref_soln%elem_hash)
               if (old_grid%element(grandparent)%isleaf) then
                  old_elem = grandparent
               else
                  old_elem = parent
               endif
            endif

! evaluate both solutions

            allocate(ux_fine(ss,1,nqpoints),uy_fine(ss,1,nqpoints), &
                     ux_old(ss,1,nqpoints),uy_old(ss,1,nqpoints), &
                     u_fine(ss,1,nqpoints),u_old(ss,1,nqpoints))
            call evaluate_soln_local(ref_soln,xquad,yquad,elem,(/(i,i=1,ss)/), &
                                     (/1/),u=u_fine,ux=ux_fine,uy=uy_fine)
            call evaluate_soln_local(old_grid,xquad,yquad,old_elem, &
                                     (/(i,i=1,ss)/),(/1/),u=u_old,ux=ux_old, &
                                     uy=uy_old)

! compute the contribution of this element

            contribution = 0.0_my_real
            do i = 1,nqpoints
               do j=1,ss
                  contribution = contribution + &
                              qweights(i)*((ux_old(j,1,i)-ux_fine(j,1,i))**2 + &
                                           (uy_old(j,1,i)-uy_fine(j,1,i))**2 + &
                                           ( u_old(j,1,i)- u_fine(j,1,i))**2)
               end do
            end do

! add contributions to integrals

            compute_refsoln_errest_elem = compute_refsoln_errest_elem + contribution
            old_grid%element_errind(old_elem,1) = &
                  old_grid%element_errind(old_elem,1) + contribution

            deallocate(u_fine,u_old,ux_fine,uy_fine,ux_old,uy_old)

! next element

         endif
         elem = ref_soln%element(elem)%next
      end do
   end do

! square roots
! Because quadrature weights can be negative, it is possible to get
! negative squared error estimates.  Use abs to avoid NaN.

   compute_refsoln_errest_elem = sqrt(abs(compute_refsoln_errest_elem))
   old_grid%element_errind = sqrt(abs(old_grid%element_errind))

end function compute_refsoln_errest_elem

!          ------------------------------------
subroutine ext_refsoln_errest_elem(grid,procs,refine_control,solver_control, &
                                   io_control,still_sequential)
!          ------------------------------------

!----------------------------------------------------
! This routine can be called from outside the module to compute the
! error estimate used by REFSOLN_ELEM
! TEMP only using first eigenvector
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type (proc_info), target, intent(in) :: procs
type(refine_options), intent(in) :: refine_control
type(solver_options), intent(in) :: solver_control
type(io_options), intent(in) :: io_control
logical, intent(in) :: still_sequential
!----------------------------------------------------
! Local variables:

type(grid_type) :: ref_soln
type(refine_options) :: loc_refcont
type(io_options) :: loc_iocont
!----------------------------------------------------
! Begin executable code

! TEMP not parallel

if (num_proc(procs) > 1) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("REFSOLN not yet implemented in parallel")
   stop
endif

! Make a copy of the grid as it currently exists.

call copy_grid(grid,ref_soln)

! Create the reference solution by performing uniform h and p refinements
! and solving on the fine grid.

loc_refcont = refine_control
loc_refcont%error_estimator = HIERARCHICAL_COEFFICIENT
loc_refcont%reftype = H_UNIFORM
loc_iocont = io_control
loc_iocont%print_linsys_when = NEVER
loc_iocont%print_error_when = NEVER
call refine(ref_soln,procs,loc_refcont,solver_control,loc_iocont, &
            still_sequential,0,0,0,0,0,.false.,no_time=.true.)
loc_refcont%reftype = P_UNIFORM
call refine(ref_soln,procs,loc_refcont,solver_control,loc_iocont, &
            still_sequential,0,0,0,0,0,.false.,no_time=.true.)
call solve(ref_soln,procs,loc_iocont,solver_control,still_sequential,.false., &
           no_time=.true.)

! MASTER doesn't participate further; it was only here to monitor solve

if (my_proc(procs) == MASTER) then
   call deallocate_grid(ref_soln)
   grid%errind_up2date = .true.
   return
endif

! Compute the error estimate.

grid%refsoln_errest = compute_refsoln_errest_elem(ref_soln,grid)
call deallocate_grid(ref_soln)
grid%errind_up2date = .true.

end subroutine ext_refsoln_errest_elem

!          -------------------------
subroutine det_and_perf_opt_ref_elem(grid,ref_soln,refine_control)
!          -------------------------

!----------------------------------------------------
! This is the top level routine for determining the optimal refinements to
! perform, and performing those refinements, with the element based version
! of refsoln.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(grid_type), intent(in) :: ref_soln
type(refine_options), intent(in) :: refine_control
!----------------------------------------------------
! Local variables:

integer :: opt_ref(size(grid%element),2)
!----------------------------------------------------
! Begin executable code

! determine optimal refinements

call det_opt_ref_elem(grid,ref_soln,opt_ref,refine_control)

! perform the refinements

call perform_refinements_elem(grid,opt_ref,refine_control)

end subroutine det_and_perf_opt_ref_elem

!          ----------------
subroutine det_opt_ref_elem(grid,ref_soln,opt_ref,refine_control)
!          ----------------

!----------------------------------------------------
! This routine determines the optimal refinement of each element, from the
! choices of p-refine once or twice and h-refine with child degrees of all
! combinations of p0, p0+1 and p0+2 where p0 = floor((p+1)/sqrt(2)).
! opt_ref(e,1:2) is (-1,-1) if e should not be refined
!                   ( p,-1) if e should be p-refined to degree p
!                   (p1,p2) if e should be h-refined with degrees p1 and p2
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid, ref_soln
integer, intent(out) :: opt_ref(:,:)
type(refine_options), intent(in) :: refine_control
!----------------------------------------------------
! Local variables:

integer :: lev, elem
real(my_real) :: max_errest
!----------------------------------------------------
! Begin executable code

! begin by marking everything as not to be refined

opt_ref = -1

! determine the maximum error estimate

max_errest = 0.0_my_real
do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (grid%element(elem)%isleaf) then
         max_errest = max(max_errest,grid%element_errind(elem,1))
      endif
      elem = grid%element(elem)%next
   end do
end do

! for each leaf element ...

do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (grid%element(elem)%isleaf) then

! if the error estimate is large enough, determine the optimal refinement

         if (grid%element_errind(elem,1) > max_errest/refine_control%inc_factor) then
            opt_ref(elem,1:2) = det_one_opt_ref_elem(elem,grid,ref_soln, &
                                                     refine_control)
         endif
      endif

! next element

      elem = grid%element(elem)%next
   end do
end do

end subroutine det_opt_ref_elem

!        --------------------
function det_one_opt_ref_elem(elem,grid,ref_soln,refcont)
!        --------------------

!----------------------------------------------------
! This routine determines the optimal refinement for element elem
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: elem
type(grid_type), intent(in) :: grid,ref_soln
integer, dimension(2) :: det_one_opt_ref_elem
type(refine_options), intent(in) :: refcont
!----------------------------------------------------
! Local variables:

integer :: p0, p_u, p_ref, child1, child2, nqpoints_child1, nqpoints_child2, &
           jerr, nbasis_elem, nbasis_child, i, j, k, l, nbasis(11), nsubset, &
           nbasis_children, adim, candidate, nbasis_u, P, Q, R, nkept, winner, &
           ss, i2
integer, allocatable :: subset(:), ipiv(:), renum(:), invrenum(:)
real(my_real) :: xvert_elem(3), yvert_elem(3), xvert_child1(3), &
                 yvert_child1(3), xvert_child2(3), yvert_child2(3), &
                 projerr(11), proj, projx, projy, projerr_u, logprojerr_u, &
                 logprojerr(11), ave, sd, biggest, reduction
real(my_real), pointer :: qweights_child1(:),xquad_child1(:),yquad_child1(:), &
                          qweights_child2(:),xquad_child2(:),yquad_child2(:)
real(my_real), allocatable ::  u_ref_child1(:,:,:), ux_ref_child1(:,:,:), &
                              uy_ref_child1(:,:,:),  u_ref_child2(:,:,:), &
                              ux_ref_child2(:,:,:), uy_ref_child2(:,:,:), &
                               basis_elem_child1(:,:), basisx_elem_child1(:,:),&
                              basisy_elem_child1(:,:),  basis_elem_child2(:,:),&
                              basisx_elem_child2(:,:), basisy_elem_child2(:,:),&
                               basis_child1(:,:), basisx_child1(:,:), &
                              basisy_child1(:,:),  basis_child2(:,:), &
                              basisx_child2(:,:), basisy_child2(:,:), &
                              ip_a(:,:),ip_b(:,:),a(:,:),b(:,:)
!----------------------------------------------------
! Begin executable code

ss = grid%system_size

! base degree to use for h-refinement

p0 = floor((grid%element(elem)%degree+1)/sqrt(2.0_my_real))

! degree of reference solution and solution

p_u = grid%element(elem)%degree
p_ref = p_u+1

! children of elem in the reference solution

child1 = get_child_lid(ref_soln%element(elem)%gid,1,ref_soln%elem_hash)
child2 = get_child_lid(ref_soln%element(elem)%gid,2,ref_soln%elem_hash)

! vertices of elem and it's children

xvert_elem = grid%vertex(grid%element(elem)%vertex)%coord%x
yvert_elem = grid%vertex(grid%element(elem)%vertex)%coord%y
xvert_child1 = (/xvert_elem(1),xvert_elem(3),(xvert_elem(1)+xvert_elem(2))/2/)
yvert_child1 = (/yvert_elem(1),yvert_elem(3),(yvert_elem(1)+yvert_elem(2))/2/)
xvert_child2 = (/xvert_elem(2),xvert_elem(3),(xvert_elem(1)+xvert_elem(2))/2/)
yvert_child2 = (/yvert_elem(2),yvert_elem(3),(yvert_elem(1)+yvert_elem(2))/2/)

! quadratures are with polynomials up to degree (p_u+2)**2
! get quadrature rules for both of the children; quadratures on elem will be
! done over the two children

call quadrature_rule_tri(min((p_u+2)**2,MAX_QUAD_ORDER_TRI), &
                         xvert_child1,yvert_child1,nqpoints_child1, &
                         qweights_child1,xquad_child1,yquad_child1,jerr, &
                         stay_in=.true.)
call quadrature_rule_tri(min((p_u+2)**2,MAX_QUAD_ORDER_TRI), &
                         xvert_child2,yvert_child2,nqpoints_child2, &
                         qweights_child2,xquad_child2,yquad_child2,jerr, &
                         stay_in=.true.)

! evaluate the reference solution and first derivatives at the quadrature points

allocate( u_ref_child1(ss,1,nqpoints_child1), &
         ux_ref_child1(ss,1,nqpoints_child1), &
         uy_ref_child1(ss,1,nqpoints_child1), &
          u_ref_child2(ss,1,nqpoints_child2), &
         ux_ref_child2(ss,1,nqpoints_child2), &
         uy_ref_child2(ss,1,nqpoints_child2))
call evaluate_soln_local(ref_soln,xquad_child1,yquad_child1,child1, &
                         (/(i,i=1,ss)/),(/1/),u_ref_child1,ux_ref_child1, &
                         uy_ref_child1)
call evaluate_soln_local(ref_soln,xquad_child2,yquad_child2,child2, &
                         (/(i,i=1,ss)/),(/1/),u_ref_child2,ux_ref_child2, &
                         uy_ref_child2)

! evaluate the basis functions of elem up to degree p+2

nbasis_elem = ((p_u+3)*(p_u+4))/2
allocate( basis_elem_child1(nbasis_elem,nqpoints_child1), &
         basisx_elem_child1(nbasis_elem,nqpoints_child1), &
         basisy_elem_child1(nbasis_elem,nqpoints_child1), &
          basis_elem_child2(nbasis_elem,nqpoints_child2), &
         basisx_elem_child2(nbasis_elem,nqpoints_child2), &
         basisy_elem_child2(nbasis_elem,nqpoints_child2))
call p_hier_basis_func(xquad_child1,yquad_child1,xvert_elem,yvert_elem, &
                       (/p_u+2,p_u+2,p_u+2,p_u+2/),"a",basis_elem_child1, &
                       basisx_elem_child1,basisy_elem_child1)
call p_hier_basis_func(xquad_child2,yquad_child2,xvert_elem,yvert_elem, &
                       (/p_u+2,p_u+2,p_u+2,p_u+2/),"a",basis_elem_child2, &
                       basisx_elem_child2,basisy_elem_child2)

! evaluate the basis functions of the children up to degree p0+2

nbasis_child = ((p0+3)*(p0+4))/2 ! each child
nbasis_children = (p0+3)**2      ! both children, some bases are in common
allocate( basis_child1(nbasis_child,nqpoints_child1), &
         basisx_child1(nbasis_child,nqpoints_child1), &
         basisy_child1(nbasis_child,nqpoints_child1), &
          basis_child2(nbasis_child,nqpoints_child2), &
         basisx_child2(nbasis_child,nqpoints_child2), &
         basisy_child2(nbasis_child,nqpoints_child2))
call p_hier_basis_func(xquad_child1,yquad_child1,xvert_child1,yvert_child1, &
                       (/p0+2,p0+2,p0+2,p0+2/),"a",basis_child1, &
                       basisx_child1,basisy_child1)
call p_hier_basis_func(xquad_child2,yquad_child2,xvert_child2,yvert_child2, &
                       (/p0+2,p0+2,p0+2,p0+2/),"a",basis_child2, &
                       basisx_child2,basisy_child2)

! compute the inner products needed for H^1 projections onto elem

adim = max(nbasis_elem,nbasis_children)
allocate(ip_a(adim,adim),ip_b(adim,ss))
ip_a = 0.0_my_real
ip_b = 0.0_my_real
do i=1,nbasis_elem
   do j=1,nbasis_elem
      do k=1,nqpoints_child1
         ip_a(i,j) = ip_a(i,j) + qweights_child1(k) * &
                            (basisx_elem_child1(i,k)*basisx_elem_child1(j,k) + &
                             basisy_elem_child1(i,k)*basisy_elem_child1(j,k) + &
                              basis_elem_child1(i,k)* basis_elem_child1(j,k))
      end do
      do k=1,nqpoints_child2
         ip_a(i,j) = ip_a(i,j) + qweights_child2(k) * &
                            (basisx_elem_child2(i,k)*basisx_elem_child2(j,k) + &
                             basisy_elem_child2(i,k)*basisy_elem_child2(j,k) + &
                              basis_elem_child2(i,k)* basis_elem_child2(j,k))
      end do
   end do
   do k=1,nqpoints_child1
      do l=1,ss
         ip_b(i,l) = ip_b(i,l) + qweights_child1(k) * &
                               (ux_ref_child1(l,1,k)*basisx_elem_child1(i,k) + &
                                uy_ref_child1(l,1,k)*basisy_elem_child1(i,k) + &
                                 u_ref_child1(l,1,k)* basis_elem_child1(i,k))
      end do
   end do
   do k=1,nqpoints_child2
      do l=1,ss
         ip_b(i,l) = ip_b(i,l) + qweights_child2(k) * &
                               (ux_ref_child2(l,1,k)*basisx_elem_child2(i,k) + &
                                uy_ref_child2(l,1,k)*basisy_elem_child2(i,k) + &
                                 u_ref_child2(l,1,k)* basis_elem_child2(i,k))
      end do
   end do
end do

! determine the subset of bases over elem used for the original element (elem,p)

allocate(subset(adim))
subset(1:p_u+2) = (/(l,l=1,p_u+2)/)
subset(p_u+3:2*p_u+1) = (/(l,l=p_u+5,2*p_u+3)/)
subset(2*p_u+2:3*p_u) = (/(l,l=2*p_u+6,3*p_u+4)/)
subset(3*p_u+1:((p_u+1)*(p_u+2))/2) = (/(l,l=3*p_u+7,((p_u+3)*(p_u+4))/2-2*p_u+1)/)

! set the size of the space (elem,p+1)

nbasis_u = ((p_u+1)*(p_u+2))/2

! copy the inner products needed for the projection on (elem,p)

allocate(a(adim,adim),b(adim,ss),ipiv(adim))
do i=1,nbasis_u
   do j=1,nbasis_u
      a(i,j) = ip_a(subset(i),subset(j))
   end do
   do j=1,ss
      b(i,j) = ip_b(subset(i),j)
   end do
end do

! solve for the projection on (elem,p)

if (my_real == kind(0.0)) then
   call sgesv(nbasis_u,ss,a,adim,ipiv,b,adim,jerr)
elseif (my_real == kind(0.0d0)) then
   call dgesv(nbasis_u,ss,a,adim,ipiv,b,adim,jerr)
else
   ierr = USER_INPUT_ERROR
   call fatal("kind of real must be single or double to use LAPACK for refsoln_elem")
   stop
endif
if (jerr /= 0) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("LAPACK returned error during H^1 projection for refsoln_elem",intlist=(/0,jerr/))
   stop
endif

! compute the projection error on (elem,p)

projerr_u = 0.0_my_real
do k=1,nqpoints_child1
   do j=1,ss
      proj  = 0.0_my_real
      projx = 0.0_my_real
      projy = 0.0_my_real
      do i=1,nbasis_u
         proj  = proj  + b(i,j)* basis_elem_child1(subset(i),k)
         projx = projx + b(i,j)*basisx_elem_child1(subset(i),k)
         projy = projy + b(i,j)*basisy_elem_child1(subset(i),k)
      end do
      projerr_u = projerr_u + qweights_child1(k) * ( &
                              (ux_ref_child1(j,1,k)-projx)**2 + &
                              (uy_ref_child1(j,1,k)-projy)**2 + &
                              ( u_ref_child1(j,1,k)-proj )**2)
   end do
end do
do k=1,nqpoints_child2
   do j=1,ss
      proj  = 0.0_my_real
      projx = 0.0_my_real
      projy = 0.0_my_real
      do i=1,nbasis_u
         proj  = proj  + b(i,j)* basis_elem_child2(subset(i),k)
         projx = projx + b(i,j)*basisx_elem_child2(subset(i),k)
         projy = projy + b(i,j)*basisy_elem_child2(subset(i),k)
      end do
      projerr_u = projerr_u + qweights_child2(k) * ( &
                              (ux_ref_child2(j,1,k)-projx)**2 + &
                              (uy_ref_child2(j,1,k)-projy)**2 + &
                              ( u_ref_child2(j,1,k)-proj )**2)
   end do
end do
projerr_u = sqrt(projerr_u)

! determine the subset of bases over elem used for (elem,p+1)

subset(1:p_u+3) = (/(l,l=1,p_u+3)/)
subset(p_u+4:2*p_u+3) = (/(l,l=p_u+5,2*p_u+4)/)
subset(2*p_u+4:3*p_u+3) = (/(l,l=2*p_u+6,3*p_u+5)/)
subset(3*p_u+4:((p_u+2)*(p_u+3))/2) = (/(l,l=3*p_u+7,((p_u+3)*(p_u+4))/2-p_u)/)

! set the size of the space (elem,p+1)

candidate = 1
nbasis(candidate) = ((p_u+2)*(p_u+3))/2
nsubset = nbasis(candidate)

! copy the inner products needed for the projection on (elem,p+1)

do i=1,nsubset
   do j=1,nsubset
      a(i,j) = ip_a(subset(i),subset(j))
   end do
   do j=1,ss
      b(i,j) = ip_b(subset(i),j)
   end do
end do

! solve for the projection on (elem,p+1)

if (my_real == kind(0.0)) then
   call sgesv(nsubset,ss,a,adim,ipiv,b,adim,jerr)
else ! my_real == kind(0.0d0)
   call dgesv(nsubset,ss,a,adim,ipiv,b,adim,jerr)
endif
if (jerr /= 0) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("LAPACK returned error during H^1 projection for refsoln_elem",intlist=(/candidate,jerr/))
   stop
endif

! compute the projection error on (elem,p+1)

projerr = 0.0_my_real
do k=1,nqpoints_child1
   do j=1,ss
      proj  = 0.0_my_real
      projx = 0.0_my_real
      projy = 0.0_my_real
      do i=1,nsubset
         proj  = proj  + b(i,j)* basis_elem_child1(subset(i),k)
         projx = projx + b(i,j)*basisx_elem_child1(subset(i),k)
         projy = projy + b(i,j)*basisy_elem_child1(subset(i),k)
      end do
      projerr(candidate) = projerr(candidate) + qweights_child1(k) * ( &
                                             (ux_ref_child1(j,1,k)-projx)**2 + &
                                             (uy_ref_child1(j,1,k)-projy)**2 + &
                                             ( u_ref_child1(j,1,k)-proj )**2)
   end do
end do
do k=1,nqpoints_child2
   do j=1,ss
      proj  = 0.0_my_real
      projx = 0.0_my_real
      projy = 0.0_my_real
      do i=1,nsubset
         proj  = proj  + b(i,j)* basis_elem_child2(subset(i),k)
         projx = projx + b(i,j)*basisx_elem_child2(subset(i),k)
         projy = projy + b(i,j)*basisy_elem_child2(subset(i),k)
      end do
      projerr(candidate) = projerr(candidate) + qweights_child2(k) * ( &
                                             (ux_ref_child2(j,1,k)-projx)**2 + &
                                             (uy_ref_child2(j,1,k)-projy)**2 + &
                                             ( u_ref_child2(j,1,k)-proj )**2)
   end do
end do
projerr(candidate) = sqrt(projerr(candidate))

! set the size of the space (elem,p+2)

candidate = candidate + 1
nbasis(candidate) = ((p_u+3)*(p_u+4))/2
nsubset = nbasis(candidate)

! determine the subset of bases over elem used for (elem,p+2)

subset(1:nsubset) = (/(l,l=1,nsubset)/)

! copy the inner products needed for the projection on (elem,p+2)

a = ip_a
b = ip_b

! solve for the projection on (elem,p+2)

if (my_real == kind(0.0)) then
   call sgesv(nsubset,ss,a,adim,ipiv,b,adim,jerr)
else ! my_real == kind(0.0d0)
   call dgesv(nsubset,ss,a,adim,ipiv,b,adim,jerr)
endif
if (jerr /= 0) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("LAPACK returned error during H^1 projection for refsoln_elem",intlist=(/candidate,jerr/))
   stop
endif

! compute the projection error on (elem,p+2)

do k=1,nqpoints_child1
   do j=1,ss
      proj  = 0.0_my_real
      projx = 0.0_my_real
      projy = 0.0_my_real
      do i=1,nsubset
         proj  = proj  + b(i,j)* basis_elem_child1(subset(i),k)
         projx = projx + b(i,j)*basisx_elem_child1(subset(i),k)
         projy = projy + b(i,j)*basisy_elem_child1(subset(i),k)
      end do
      projerr(candidate) = projerr(candidate) + qweights_child1(k) * ( &
                                             (ux_ref_child1(j,1,k)-projx)**2 + &
                                             (uy_ref_child1(j,1,k)-projy)**2 + &
                                             ( u_ref_child1(j,1,k)-proj )**2)
   end do
end do
do k=1,nqpoints_child2
   do j=1,ss
      proj  = 0.0_my_real
      projx = 0.0_my_real
      projy = 0.0_my_real
      do i=1,nsubset
         proj  = proj  + b(i,j)* basis_elem_child2(subset(i),k)
         projx = projx + b(i,j)*basisx_elem_child2(subset(i),k)
         projy = projy + b(i,j)*basisy_elem_child2(subset(i),k)
      end do
      projerr(candidate) = projerr(candidate) + qweights_child2(k) * ( &
                                             (ux_ref_child2(j,1,k)-projx)**2 + &
                                             (uy_ref_child2(j,1,k)-projy)**2 + &
                                             ( u_ref_child2(j,1,k)-proj )**2)
   end do
end do
projerr(candidate) = sqrt(projerr(candidate))

deallocate(basis_elem_child1,basisx_elem_child1,basisy_elem_child1, &
           basis_elem_child2,basisx_elem_child2,basisy_elem_child2)

! the basis function indices for bases of degree p0+2 are (where child1 is on
! the left, p is p0 and P=(p+3)(p+4)/2)
!
!                                   2
!                                  / \
!                                /  |  \
!                              /    |    \
!                            /      |<---------- 4..p+4
!                          /        |        \
!            2p+6..3p+6  /          |          \  P+p+3+..P+2p+3
!                      /  3p+7..P   | P+2p+4..   \
!                    /              | (p+3)^2      \
!                  /                |                \
!                 -------------------------------------
!                1   p+5..2p+5       3  P+2..P+p+2    P+1
!
! The indices for the subset of bases where the first child has degree p0+i
! and the second child has degree p0+j, i,j=0..2, with k=min(i,j), p is p0,
! Q = (p+i-2)(p+i-1)/2 + 3p+k+2i, R = Q+2p+2j + (p+j-2)(p+j-1)/2 - 1
!
!                                   2
!                                  / \
!                                /  |  \
!                              /    |    \
!                            /      |<---------- 4..p+k+2
!                          /        |        \
!     2p+k+i+2..3p+k+2i  /          |          \  Q+p+j+1..Q+2p+2j-1
!                      /  3p+k+2i+1 | Q+2p+2j    \
!                    /    ..Q       | ..R          \
!                  /                |                \
!                 -------------------------------------
!                1  p+k+3..2p+k+i+1  3  Q+2..Q+p+j    Q+1
!
! determine a renumbering of the degree p0+2 bases in child 2 from the local
! numbering returned by p_hier_basis_func to the global numbering of the
! above diagram, and its inverse from the global numbering to the local
! numbering using -1 for global indices that don't have a local basis

allocate(renum(nbasis_child),invrenum(nbasis_children))
P = ((p0+3)*(p0+4))/2
renum(1) = P+1
renum(2:p0+4) = (/(l,l=2,p0+4)/)
renum(p0+5:P) = (/(l,l=P+2,(p0+3)**2)/)
invrenum = -1
do i=1,nbasis_child
   invrenum(renum(i)) = i
end do

! compute the inner products of the bases in the children
! note that bases 2, 3, ..., p0+4 of the two children are the same basis
! functions (the bases associated with the common vertices and common edge)

! compute the inner products needed for H^1 projections onto the children

ip_a = 0.0_my_real
ip_b = 0.0_my_real
do i=1,nbasis_child
   do j=1,nbasis_child
      do k=1,nqpoints_child1
         ip_a(i,j) = ip_a(i,j) + qweights_child1(k) * &
                                 (basisx_child1(i,k)*basisx_child1(j,k) + &
                                  basisy_child1(i,k)*basisy_child1(j,k) + &
                                   basis_child1(i,k)* basis_child1(j,k))
      end do
   end do
   do k=1,nqpoints_child1
      do j=1,ss
         ip_b(i,j) = ip_b(i,j) + qweights_child1(k) * &
                          (ux_ref_child1(j,1,k)*basisx_child1(i,k) + &
                           uy_ref_child1(j,1,k)*basisy_child1(i,k) + &
                            u_ref_child1(j,1,k)* basis_child1(i,k))
      end do
   end do
end do
do i=1,nbasis_child
   do j=1,nbasis_child
      do k=1,nqpoints_child2
         ip_a(renum(i),renum(j)) = ip_a(renum(i),renum(j)) + &
                                   qweights_child2(k) * &
                                   (basisx_child2(i,k)*basisx_child2(j,k) + &
                                    basisy_child2(i,k)*basisy_child2(j,k) + &
                                     basis_child2(i,k)* basis_child2(j,k))
      end do
   end do
   do k=1,nqpoints_child2
      do j=1,ss
         ip_b(renum(i),j) = ip_b(renum(i),j) + qweights_child2(k) * &
                             (ux_ref_child2(j,1,k)*basisx_child2(i,k) + &
                              uy_ref_child2(j,1,k)*basisy_child2(i,k) + &
                               u_ref_child2(j,1,k)* basis_child2(i,k))
      end do
   end do
end do

! for each combination of degrees on the children ..

do i=0,2
   do j=0,2
      candidate = candidate + 1

! reject any h-refinements where both children have degree > p_u

      if (p0+i > p_u .and. p0+j > p_u) then
         projerr(candidate) = -1.0_my_real
         cycle
      endif

      k = min(i,j)
      Q = ((p0+i-2)*(p0+i-1))/2 + 3*p0+k+2*i
      R = Q+2*p0+2*j + ((p0+j-2)*(p0+j-1))/2 - 1

! determine the subset of bases used for (p0+i,p0+j)

      subset(1:p0+k+2) = (/(l,l=1,p0+k+2)/)
      subset(p0+k+3:2*p0+k+i+1) = (/(l,l=p0+5,2*p0+3+i)/)
      subset(2*p0+k+i+2:3*p0+k+2*i) = (/(l,l=2*p0+6,3*p0+4+i)/)
      subset(3*p0+k+2*i+1:Q) = (/(l,l=3*p0+7,P)/)
      subset(Q+1) = P+1
      subset(Q+2:Q+p0+j) = (/(l,l=P+2,P+p0+j)/)
      subset(Q+p0+j+1:Q+2*p0+2*j-1) = (/(l,l=P+p0+3,P+2*p0+j+1)/)
      subset(Q+2*p0+2*j:R) = (/(l,l=P+2*p0+4,(p0+3)**2+j-2)/)

! set the size of the space for (p0+i,p0+j)

      nbasis(candidate) = R
      nsubset = nbasis(candidate)

! copy the inner products needed for the projection on (p0+i,p0+j)

      do k=1,nsubset
         do l=1,nsubset
            a(k,l) = ip_a(subset(k),subset(l))
         end do
         do l=1,ss
            b(k,l) = ip_b(subset(k),l)
         end do
      end do

! solve for the projection on (p0+i,p0+j)

      if (my_real == kind(0.0)) then
         call sgesv(nsubset,ss,a,adim,ipiv,b,adim,jerr)
      else ! my_real == kind(0.0d0)
         call dgesv(nsubset,ss,a,adim,ipiv,b,adim,jerr)
      endif
      if (jerr /= 0) then
         ierr = PHAML_INTERNAL_ERROR
         call warning("LAPACK returned error during H^1 projection for refsoln_elem",intlist=(/candidate,jerr/))
         projerr(candidate) = -1.0_my_real
      endif

! compute the projection error on (p0+i,p0+j)

      do k=1,nqpoints_child1
         do i2=1,ss
            proj  = 0.0_my_real
            projx = 0.0_my_real
            projy = 0.0_my_real
            do l=1,nsubset
               if (subset(l) <= P) then
                  proj  = proj  + b(l,i2)* basis_child1(subset(l),k)
                  projx = projx + b(l,i2)*basisx_child1(subset(l),k)
                  projy = projy + b(l,i2)*basisy_child1(subset(l),k)
               endif
            end do
            projerr(candidate) = projerr(candidate) + qweights_child1(k) * ( &
                                            (ux_ref_child1(i2,1,k)-projx)**2 + &
                                            (uy_ref_child1(i2,1,k)-projy)**2 + &
                                            ( u_ref_child1(i2,1,k)-proj )**2)
         end do
      end do
      do k=1,nqpoints_child2
         do i2=1,ss
            proj  = 0.0_my_real
            projx = 0.0_my_real
            projy = 0.0_my_real
            do l=1,nsubset
               if (invrenum(subset(l)) /= -1) then
                  proj  = proj  + b(l,i2)* basis_child2(invrenum(subset(l)),k)
                  projx = projx + b(l,i2)*basisx_child2(invrenum(subset(l)),k)
                  projy = projy + b(l,i2)*basisy_child2(invrenum(subset(l)),k)
               endif
            end do
            projerr(candidate) = projerr(candidate) + qweights_child2(k) * ( &
                                            (ux_ref_child2(i2,1,k)-projx)**2 + &
                                            (uy_ref_child2(i2,1,k)-projy)**2 + &
                                            ( u_ref_child2(i2,1,k)-proj )**2)
         end do
      end do
      projerr(candidate) = sqrt(projerr(candidate))

! next candidate

   end do
end do

deallocate(u_ref_child1,ux_ref_child1, uy_ref_child1,u_ref_child2, &
           ux_ref_child2,uy_ref_child2,basis_child1,basisx_child1, &
           basisy_child1,basis_child2,basisx_child2,basisy_child2, &
           ip_a,ip_b,subset,a,b,ipiv,renum,invrenum)

! discard candidates that do not reduce the projection error; flag them
! by setting the projection error to -1
! Also discard any that are close to 0 to avoid problems with log.

where (projerr > projerr_u) projerr = -1.0_my_real
where (abs(projerr) < 1.0e-12) projerr = -1.0_my_real
nkept = 11 - count(projerr==-1.0_my_real)

! if none were kept, no refinement

if (nkept <= 0) then
   det_one_opt_ref_elem = (/-1,-1/)
else

! compute the logarithms of the projection errors

   logprojerr_u = log10(projerr_u)
   where (projerr /= -1.0_my_real) logprojerr = log10(projerr)

! compute the average and standard deviation of the logarithms of projerr

   ave = 0.0_my_real
   do candidate=1,11
      if (projerr(candidate) /= -1.0_my_real) then
         ave = ave + logprojerr(candidate)
      endif
   end do
   ave = ave / nkept

   sd = 0.0_my_real
   do candidate=1,11
      if (projerr(candidate) /= -1.0_my_real) then
         sd = sd + (logprojerr(candidate)-ave)**2
      endif
   end do
   sd = sqrt(sd/nkept)
   
! select the winning candidate

   winner = 0
   biggest = -huge(0.0_my_real)
   do candidate=1,11
      if (projerr(candidate) == -1.0_my_real .or. &
          logprojerr(candidate) >= ave+sd .or. &
          nbasis(candidate) <= nbasis_u) cycle
      reduction = (logprojerr_u-logprojerr(candidate))/(nbasis(candidate)-nbasis_u)
      if (candidate==1 .or. candidate==2) reduction=refcont%refsoln_pbias*reduction
      if (reduction > biggest) then
         winner = candidate
         biggest = reduction
      endif
   end do

! return the degree(s) of the winner

   select case (winner)
   case (0)
      det_one_opt_ref_elem = (/-1,-1/)
   case (1)
      det_one_opt_ref_elem = (/p_u+1,-1/)
   case (2)
      det_one_opt_ref_elem = (/p_u+2,-1/)
   case default
      j = winner/3 - 1
      i = winner - 3*(j+1)
      det_one_opt_ref_elem = (/p0+i,p0+j/)
   end select

endif

end function det_one_opt_ref_elem

!          ------------------------
subroutine perform_refinements_elem(grid,opt_ref,refine_control)
!          ------------------------

!----------------------------------------------------
! This routine performs the refinements indicated in opt_ref.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(inout) :: opt_ref(:,:)
type(refine_options), intent(in) :: refine_control
!----------------------------------------------------
! Local variables:

integer :: lev, elem, err, children(MAX_CHILD), reduce_p_max
!----------------------------------------------------
! Begin executable code

! do all h refinements first, because some h refinements of mates might
! remove p refinements

! for each element

do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)

! no refinement if first opt_ref is -1

      if (opt_ref(elem,1) /= -1) then

! if second opt_ref is not -1, h-refinement is indicated, but make sure it
! hasn't already been h-refined as a mate

         if (opt_ref(elem,2) /= -1) then
            if (grid%element(elem)%isleaf) then
               call bisect_triangle_pair(grid,elem,err,refine_control, &
                                         reduce_p=0,reduce_p_max=reduce_p_max)
            endif

! use p-refinement to set the degrees of the children

            children = get_child_lid(grid%element(elem)%gid,ALL_CHILDREN, &
                                     grid%elem_hash)
            do while (grid%element(children(1))%degree > opt_ref(elem,1))
               call p_coarsen_elem(grid,children(1),err,refine_control)
            end do
            do while (grid%element(children(2))%degree > opt_ref(elem,2))
               call p_coarsen_elem(grid,children(2),err,refine_control)
            end do
            do while (grid%element(children(1))%degree < opt_ref(elem,1) .and. &
                      grid%element(children(1))%degree < refine_control%max_deg)
               call p_refine_elem(grid,children(1),refine_control)
            end do
            do while (grid%element(children(2))%degree < opt_ref(elem,2) .and. &
                      grid%element(children(2))%degree < refine_control%max_deg)
               call p_refine_elem(grid,children(2),refine_control)
            end do

         endif ! opt_ref(2) /= -1
      endif ! opt_ref(1) /= -1

! next element

      elem = grid%element(elem)%next
   end do
end do

! now do p refinements

! for each element

do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)

! no refinement if first opt_ref is -1
! Note that all of opt_ref was initialized at -1, so children of h-refined
! elements have the first opt_ref at -1

      if (opt_ref(elem,1) /= -1) then

! if the second opt_ref is -1, p-refinement is indicated, but if it was
! h-refined for compatibility, don't do it

         if (opt_ref(elem,2) == -1) then

            if (grid%element(elem)%isleaf) then
               do while (grid%element(elem)%degree < opt_ref(elem,1) .and. &
                         grid%element(elem)%degree < refine_control%max_deg)
                  call p_refine_elem(grid,elem,refine_control)
               end do
            endif

         endif ! opt_ref(2) == -1
      endif ! opt_ref(1) /= -1

! next element

      elem = grid%element(elem)%next
   end do
end do
end subroutine perform_refinements_elem

!---------------------------------------------------------------------
! INTER-PARTITION COMMUNICATION ROUTINES
!---------------------------------------------------------------------

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
call phaml_alltoall(procs,isend,nsend,irecv,nrecv,470)
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
call phaml_alltoall(procs,isend,nsend,irecv,nrecv,480)
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
grid%dof = grid%dof + newdeg - olddeg
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

end module grid_mod

!---------------------------------------------------------------------
! EXTERNAL ROUTINES THAT INTERFACE TO ALGENCAN FOR HP_NLP
!---------------------------------------------------------------------

!     =================================================================
!     File: toyprob.f
!     =================================================================

!     =================================================================
!     Module: Subroutines that define the problem
!     =================================================================

!     Last update of any of the component of this module: 
 
!     March 4, 2008.

!     Users are encouraged to download periodically updated versions of 
!     this code at the TANGO home page:
 
!     www.ime.usp.br/~egbirgin/tango/ 

!     ******************************************************************
!     ******************************************************************

      subroutine inip(n,x,l,u,m,lambda,equatn,linear,coded,checkder,init_guess)

      use global
      use nlp_vars
      implicit none

!     SCALAR ARGUMENTS
      integer m,n
      logical checkder

!     ARRAY ARGUMENTS
      logical coded(10),equatn(*),linear(*)
      double precision l(*),lambda(*),u(*),x(*)
      real(my_real) :: init_guess(*)

!     This subroutine must set some problem data. For achieving this 
!     objective YOU MUST MODIFY it according to your problem. See below 
!     where your modifications must be inserted.
!     
!     Parameters of the subroutine:
!
!     On Entry:
!
!     This subroutine has no input parameters.
!
!     On Return
!
!     n        integer,
!              number of variables,
!
!     x        double precision x(n),
!              initial point,
!
!     l        double precision l(n),
!              lower bounds on x,
!
!     u        double precision u(n),
!              upper bounds on x,
!
!     m        integer,
!              number of constraints (excluding the bounds),
!
!     lambda   double precision lambda(m),
!              initial estimation of the Lagrange multipliers,
!
!     equatn   logical equatn(m)
!              for each constraint j, set equatn(j) = .true. if it is an 
!              equality constraint of the form c_j(x) = 0, and set 
!              equatn(j) = .false. if it is an inequality constraint of 
!              the form c_j(x) <= 0,
!
!     linear   logical linear(m)
!              for each constraint j, set linear(j) = .true. if it is a 
!              linear constraint, and set linear(j) = .false. if it is a
!              nonlinear constraint.

!     LOCAL SCALARS
      integer i

!     ******************************************************************
!     FROM HERE ON YOU MUST MODIFY THE SUBROUTINE TO SET YOUR PROBLEM
!     DATA:
!     ******************************************************************

!     Number of variables

      n = 2*nlp_nelem

!     Initial point

      do i = 1,n
          x(i) = init_guess(i)
      end do

!     Lower and upper bounds

      do i = 1,nlp_nelem
          l(i) = max(1, &
                     hold_grid%element(nlp_elem_list(i))%level - &
                     hold_refcont%nlp_max_h_dec)
          u(i) = min(hold_refcont%max_lev, &
                     hold_grid%element(nlp_elem_list(i))%level + &
                     hold_refcont%nlp_max_h_inc)
          l(i+nlp_nelem) = max(2, &
                    hold_grid%element(nlp_elem_list(i))%degree - &
                     hold_refcont%nlp_max_p_dec)
          u(i+nlp_nelem) = min(hold_refcont%max_deg, &
                    hold_grid%element(nlp_elem_list(i))%degree + &
                     hold_refcont%nlp_max_p_inc)
      end do

!     Number of constraints (equalities plus inequalities)

      m = 3*nlp_nelem + 1

!     Lagrange multipliers approximation. Most users prefer to use the 
!     null initial Lagrange multipliers estimates. However, if the 
!     problem that you are solving is "slightly different" from a 
!     previously solved problem of which you know the correct Lagrange 
!     multipliers, we encourage you to set these multipliers as initial 
!     estimates. Of course, in this case you are also encouraged to use 
!     the solution of the previous problem as initial estimate of the 
!     solution. Similarly, most users prefer to use rho = 10 as initial 
!     penalty parameters. But in the case mentioned above (good 
!     estimates of solution and Lagrange multipliers) larger values of 
!     the penalty parameters (say, rho = 1000) may be more useful. More 
!     warm-start procedures are being elaborated.

      do i = 1,m
          lambda(i) =  0.0d0
      end do

!     For each constraint i, set equatn(i) = .true. if it is an equality
!     constraint of the form c_i(x) = 0, and set equatn(i) = .false. if 
!     it is an inequality constraint of the form c_i(x) <= 0.

      do i = 1,m
          equatn(i) = .false.
      end do

!     For each constraint i, set linear(i) = .true. if it is a linear
!     constraint, otherwise set linear(i) = .false.

      do i = 1,m-1
         linear(i) = .true.
      end do
      linear(m) =  .false.

!     Indicate which subroutines did you code.

      coded( 1) = .true.  ! evalf    
      coded( 2) = .true.  ! evalg   
      coded( 3) = .true.  ! evalh   
      coded( 4) = .true.  ! evalc   
      coded( 5) = .true.  ! evaljac
      coded( 6) = .true.  ! evalhc  
      coded( 7) = .false. ! evalfc
      coded( 8) = .false. ! evalgjac
      coded( 9) = .false. ! evalhl
      coded(10) = .false. ! evalhlp

!     Set checkder = TRUE if you code some derivatives and you would
!     like them to be tested by finite differences. It is highly 
!     recommended.

      checkder = .false.

      end subroutine

!     ******************************************************************
!     ******************************************************************

      subroutine evalf(n,x,f,flag)

      use grid_mod
      implicit none

!     SCALAR ARGUMENTS
      integer flag,n
      double precision f

!     ARRAY ARGUMENTS
      double precision x(n)

!     This subroutine must compute the objective function. For achieving 
!     this objective YOU MUST MODIFY it according to your problem. See 
!     below where your modifications must be inserted.
!     
!     Parameters of the subroutine:
!
!     On Entry:
!
!     n        integer,
!              number of variables,
!
!     x        double precision x(n),
!              current point,
!
!     On Return
!
!     f        double precision,
!              objective function value at x,
!
!     flag     integer,
!              You must set it to any number different of 0 (zero) if 
!              some error ocurred during the evaluation of the objective 
!              function. (For example, trying to compute the square root 
!              of a negative number, dividing by zero or a very small 
!              number, etc.) If everything was o.k. you must set it 
!              equal to zero.

!     Objective function

!     ******************************************************************
!     FROM HERE ON YOU MUST MODIFY THE SUBROUTINE TO SET YOUR OBJECTIVE
!     FUNCTION:
!     ******************************************************************

      flag = 0

      f = nlp_dof(x)

!     ******************************************************************
!     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALF.
!     ******************************************************************

      end subroutine

!     ******************************************************************
!     ******************************************************************
 
      subroutine evalg(n,x,g,flag)

      use grid_mod
      implicit none

!     SCALAR ARGUMENTS
      integer flag,n

!     ARRAY ARGUMENTS
      double precision g(n),x(n)

!     This subroutine must compute the gradient vector of the objective 
!     function. For achieving these objective YOU MUST MODIFY it in the 
!     way specified below. However, if you decide to use numerical 
!     derivatives (we dont encourage this option at all!) you dont need
!     to modify evalg.
!
!     Parameters of the subroutine:
!
!     On Entry:
!
!     n        integer,
!              number of variables,
!
!     x        double precision x(n),
!              current point,
!
!     On Return
!
!     g        double precision g(n),
!              gradient vector of the objective function evaluated at x,
!
!     flag     integer,
!              You must set it to any number different of 0 (zero) if 
!              some error ocurred during the evaluation of any component 
!              of the gradient vector. (For example, trying to compute 
!              the square root of a negative number, dividing by zero or 
!              a very small number, etc.) If everything was o.k. you 
!              must set it equal to zero.

!     Gradient vector

!     ******************************************************************
!     FROM HERE ON YOU MUST MODIFY THE SUBROUTINE TO SET THE GRADIENT
!     VECTOR OF YOUR OBJECTIVE FUNCTION: 
!     ******************************************************************

      flag = 0

      call nlp_ddof(x,g)

!     ******************************************************************
!     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALG. 
!     ******************************************************************
 
      end subroutine

!     ******************************************************************
!     ******************************************************************
 
      subroutine evalh(n,x,hlin,hcol,hval,hnnz,flag)

      use nlp_vars
      implicit none

!     SCALAR ARGUMENTS
      integer flag,n,hnnz

!     ARRAY ARGUMENTS
      integer hcol(*),hlin(*)
      double precision hval(*),x(n)

!     This subroutine might compute the Hessian matrix of the objective 
!     function. For achieving this objective YOU MAY MODIFY it according 
!     to your problem. To modify this subroutine IS NOT MANDATORY. See 
!     below where your modifications must be inserted.
!     
!     Parameters of the subroutine:
!
!     On Entry:
!
!     n        integer,
!              number of variables,
!
!     x        double precision x(n),
!              current point,
!
!     On Return
!
!     hnnz     integer,
!              number of perhaps-non-null elements of the computed 
!              Hessian,
!
!     hlin     integer hlin(hnnz),
!              see below,
!
!     hcol     integer hcol(hnnz),
!              see below,
!
!     hval     double precision hval(hnnz),
!              the non-null value of the (hlin(k),hcol(k)) position 
!              of the Hessian matrix of the objective function must 
!              be saved at hval(k). Just the lower triangular part of
!              Hessian matrix must be computed,
!
!     flag     integer,
!              You must set it to any number different of 0 (zero) if 
!              some error ocurred during the evaluation of the Hessian
!              matrix of the objective funtion. (For example, trying 
!              to compute the square root of a negative number, 
!              dividing by zero or a very small number, etc.) If 
!              everything was o.k. you must set it equal to zero.

!     ******************************************************************
!     FROM HERE ON YOU MAY (OPTIONALY) MODIFY THE SUBROUTINE TO SET THE 
!     HESSIAN MATRIX OF YOUR OBJECTIVE FUNCTION: 
!     ******************************************************************

      integer i, elem, l
      double precision lhat, phat

      flag = 0

      hnnz = 3*nlp_nelem

      do i=1,nlp_nelem
         elem = nlp_elem_list(i)
         l = hold_grid%element(elem)%level
         lhat = x(i)
         phat = x(i+nlp_nelem)
         hlin(i) = i
         hcol(i) = i
         hval(i) = log(2.0d0)**2 * 2**(lhat-l) * phat**2 / 2
         hlin(i+nlp_nelem) = i+nlp_nelem
         hcol(i+nlp_nelem) = i
         hval(i+nlp_nelem) =  log(2.0d0) * 2**(lhat-l) * phat
         hlin(i+2*nlp_nelem) = i+nlp_nelem
         hcol(i+2*nlp_nelem) = i+nlp_nelem
         hval(i+2*nlp_nelem) = 2**(lhat-l)
      end do

!     ******************************************************************
!     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALH. 
!     ******************************************************************
 
      end subroutine

!     ******************************************************************
!     ******************************************************************
 
      subroutine evalc(n,x,ind,c,flag)

      use gridtype_mod
      use nlp_vars
      use grid_mod
      implicit none

!     SCALAR ARGUMENTS
      integer ind,flag,n
      double precision c

!     ARRAY ARGUMENTS
      double precision x(n)

!     This subroutine must compute the ind-th constraint of your 
!     problem. For achieving this objective YOU MUST MOFIFY it 
!     according to your problem. See below the places where your 
!     modifications must be inserted.
!
!     Parameters of the subroutine:
!
!     On Entry:
!
!     n        integer,
!              number of variables,
!
!     x        double precision x(n),
!              current point,
!
!     ind      integer,
!              index of the constraint to be computed,
!
!     On Return
!
!     c        double precision,
!              ind-th constraint evaluated at x,
!
!     flag     integer
!              You must set it to any number different of 0 (zero) if 
!              some error ocurred during the evaluation of the 
!              constraint. (For example, trying to compute the square 
!              root of a negative number, dividing by zero or a very 
!              small number, etc.) If everything was o.k. you must set 
!              it equal to zero.
 
!     LOCAL VARIABLES
      integer elem, neigh_index, all_neighbors(3), neigh

!     Constraints

!     ******************************************************************
!     FROM HERE ON YOU MUST MODIFY THE SUBROUTINE TO SET YOUR CONSTRAINTS 
!     ******************************************************************

      flag = 0

! the nonlinear error estimate constraint

      if (ind == 1) then
         c = nlp_errest(x) - 1.0d0

! a linear constraint that neighbors differ in level by at most 1
! If the neighbor is the boundary, then put in a (redundant) constraint
! for max_lev, because it would be too hard to number the constraints otherwise

      else
         elem = nlp_elem_list((ind+1)/3)
         neigh_index = 1+mod(ind-2,3)
         all_neighbors = get_neighbors(elem,hold_grid)
         neigh = all_neighbors(neigh_index)
         if (neigh == BOUNDARY) then
            c = x(nlp_inverse_elem(elem)) - hold_refcont%max_lev
         else
            c = x(nlp_inverse_elem(elem)) - x(nlp_inverse_elem(neigh)) - 1
         endif

      endif

!     ******************************************************************
!     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALC. 
!     ******************************************************************
 
      end subroutine

!     ******************************************************************
!     ******************************************************************
 
      subroutine evaljac(n,x,ind,jcvar,jcval,jcnnz,flag)

      use nlp_vars
      use grid_mod
      implicit none

!     SCALAR ARGUMENTS
      integer flag,ind,jcnnz,n

!     ARRAY ARGUMENTS
      integer jcvar(n)
      double precision x(n),jcval(n)

!     This subroutine must compute the gradient of the ind-th constraint.
!     For achieving these objective YOU MUST MODIFY it in the way 
!     specified below.
!
!     Parameters of the subroutine:
!
!     On Entry:
!
!     n        integer,
!              number of variables,
!
!     x        double precision x(n),
!              current point,
!
!     ind      integer,
!              index of the constraint whose gradient will be computed,
!
!     On Return
!
!     jcnnz    integer,
!              number of perhaps-non-null elements of the computed 
!              gradient,
!
!     jcvar    integer jcvar(jcnnz),
!              see below,
!
!     jcval    double precision jcval(jcnnz),
!              the non-null value of the partial derivative of the 
!              ind-th constraint with respect to the jcvar(k)-th 
!              variable must be saved at jcval(k).
!
!     flag     integer
!              You must set it to any number different of 0 (zero) if 
!              some error ocurred during the evaluation of the 
!              constraint. (For example, trying to compute the square 
!              root of a negative number, dividing by zero or a very 
!              small number, etc.) If everything was o.k. you must set 
!              it equal to zero.

!     Sparse gradient vector of the ind-th constraint

!     ******************************************************************
!     FROM HERE ON YOU MUST MODIFY THE SUBROUTINE TO SET THE GRADIENTS  
!     OF YOUR CONSTRAINTS: 
!     ******************************************************************

      integer :: i, elem, neigh_index, all_neighbors(3), neigh

      flag = 0

! the nonlinear error estimate constraint

      if (ind == 1) then
         jcnnz = n
         jcvar = (/(i,i=1,n)/)
         call nlp_derrest(x,jcval)

! a linear constraint that neighbors differ in level by at most 1

      else
         elem = nlp_elem_list((ind+1)/3)
         neigh_index = 1+mod(ind-2,3)
         all_neighbors = get_neighbors(elem,hold_grid)
         neigh = all_neighbors(neigh_index)
         if (neigh == BOUNDARY) then
            jcnnz = 1
            jcvar(1) = nlp_inverse_elem(elem)
            jcval(1) = 1.0d0
         else
            jcnnz = 2
            jcvar(1) = nlp_inverse_elem(elem)
            jcval(1) = 1.0d0
            jcvar(2) = nlp_inverse_elem(neigh)
            jcval(2) = -1.0d0
         endif

      endif

!     ******************************************************************
!     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALJAC. 
!     ******************************************************************
 
      end subroutine

!     ******************************************************************
!     ******************************************************************
 
      subroutine evalhc(n,x,ind,hclin,hccol,hcval,hcnnz,flag)

      use nlp_vars
      implicit none

!     SCALAR ARGUMENTS
      integer flag,hcnnz,ind,n

!     ARRAY ARGUMENTS
      integer hccol(*),hclin(*)
      double precision hcval(*),x(n)

!     This subroutine might compute the Hessian matrix of the ind-th
!     constraint. For achieving this objective YOU MAY MODIFY it 
!     according to your problem. To modify this subroutine IS NOT 
!     MANDATORY. See below where your modifications must be inserted.
!     
!     Parameters of the subroutine:
!
!     On Entry:
!
!     n        integer,
!              number of variables,
!
!     x        double precision x(n),
!              current point,
!
!     ind      integer,
!              index of the constraint whose Hessian will be computed,
!
!     On Return
!
!     hcnnz    integer,
!              number of perhaps-non-null elements of the computed 
!              Hessian,
!
!     hclin    integer hclin(hcnnz),
!              see below,
!
!     hccol    integer hccol(hcnnz),
!              see below,
!
!     hcval    double precision hcval(hcnnz),
!              the non-null value of the (hclin(k),hccol(k)) position 
!              of the Hessian matrix of the ind-th constraint must 
!              be saved at hcval(k). Just the lower triangular part of
!              Hessian matrix must be computed,
!
!     flag     integer,
!              You must set it to any number different of 0 (zero) if 
!              some error ocurred during the evaluation of the Hessian
!              matrix of the ind-th constraint. (For example, trying 
!              to compute the square root of a negative number, 
!              dividing by zero or a very small number, etc.) If 
!              everything was o.k. you must set it equal to zero.

!     ******************************************************************
!     FROM HERE ON YOU MAY (OPTIONALY) MODIFY THE SUBROUTINE TO SET THE 
!     HESSIANS OF YOUR CONSTRAINTS: 
!     ******************************************************************

      integer i, elem, l
      double precision m, p, eta, lhat, phat

      flag = 0

      if (ind == 1) then

         hcnnz = 3*nlp_nelem

         do i=1,nlp_nelem
            elem = nlp_elem_list(i)
            l = hold_grid%element(elem)%level
            lhat = x(i)
            p = hold_grid%element(elem)%degree
            phat = x(i+nlp_nelem)
            eta = hold_grid%element_errind(elem,1)
            m = nlp_m(i)
            hclin(i) = i
            hccol(i) = i
            hcval(i) = (log(2.0d0)*min(p,m-1))**2 * &
                       2**(min(p,m-1)*(l-lhat)) * (p/phat)**(2*(m-1)) * eta**2
            hclin(i+nlp_nelem) = i+nlp_nelem
            hccol(i+nlp_nelem) = i
            hcval(i+nlp_nelem) = log(2.0d0)*min(p,m-1)*2*(m-1) * &
                 2**(min(p,m-1)*(l-lhat)) * (p/phat)**(2*(m-1)) * eta**2 / phat
            hclin(i+2*nlp_nelem) = i+nlp_nelem
            hccol(i+2*nlp_nelem) = i+nlp_nelem
            hcval(i+2*nlp_nelem) = (2-2*m)*(2-2*m-1) * &
              2**(min(p,m-1)*(l-lhat)) * (p/phat)**(2*(m-1)) * eta**2 / phat**2
         end do
         hcval(1:3*nlp_nelem) = hcval(1:3*nlp_nelem)/nlp_tau

      else

         hcnnz = 0

      endif

!     ******************************************************************
!     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALHC. 
!     ******************************************************************
 
      end subroutine

!     ******************************************************************
!     ******************************************************************

      subroutine evalfc(n,x,f,m,c,flag)

      implicit none

!     SCALAR ARGUMENTS
      integer flag,m,n
      double precision f

!     ARRAY ARGUMENTS
      double precision c(m),x(n)

      flag = - 1

      end subroutine

!     ******************************************************************
!     ******************************************************************
 
      subroutine evalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,flag)

      implicit none

!     SCALAR ARGUMENTS
      integer flag,jcnnz,m,n

!     ARRAY ARGUMENTS
      integer jcfun(*),jcvar(*)
      double precision g(n),jcval(*),x(n)

      flag = - 1

      end subroutine

!     ******************************************************************
!     ******************************************************************

      subroutine evalhl(n,x,m,lambda,scalef,scalec,hllin,hlcol,hlval, &
      hlnnz,flag)

      implicit none

!     SCALAR ARGUMENTS
      integer flag,hlnnz,m,n
      double precision scalef

!     ARRAY ARGUMENTS
      integer hlcol(*),hllin(*)
      double precision hlval(*),lambda(m),scalec(m),x(n)

      flag = - 1

      end subroutine

!     ******************************************************************
!     ******************************************************************

      subroutine evalhlp(n,x,m,lambda,sf,sc,p,hp,goth,flag)

      implicit none

!     SCALAR ARGUMENTS
      logical goth
      integer flag,m,n
      double precision sf

!     ARRAY ARGUMENTS
      double precision hp(n),lambda(m),p(n),sc(m),x(n)

!     This subroutine might compute the product of the Hessian of the
!     Lagrangian times vector p (just the Hessian of the objective 
!     function in the unconstrained or bound-constrained case). 
!     
!     Parameters of the subroutine:
!
!     On Entry:
!
!     n        integer,
!              number of variables,
!
!     x        double precision x(n),
!              current point,
!
!     m        integer,
!              number of constraints,
!
!     lambda   double precision lambda(m),
!              vector of Lagrange multipliers,
!
!     p        double precision p(n),
!              vector of the matrix-vector product,
!
!     goth     logical,
!              can be used to indicate if the Hessian matrices were
!              computed at the current point. It is set to .false.
!              by the optimization method every time the current
!              point is modified. Sugestion: if its value is .false. 
!              then compute the Hessians, save them in a common 
!              structure and set goth to .true.. Otherwise, just use 
!              the Hessians saved in the common block structure,
!
!     On Return
!
!     hp       double precision hp(n),
!              Hessian-vector product,
!
!     goth     logical,
!              see above,
!              
!     flag     integer,
!              You must set it to any number different of 0 (zero) if 
!              some error ocurred during the evaluation of the 
!              Hessian-vector product. (For example, trying to compute 
!              the square root of a negative number, dividing by zero 
!              or a very small number, etc.) If everything was o.k. you 
!              must set it equal to zero.

      flag = - 1

      end subroutine

!     ******************************************************************
!     ******************************************************************

      subroutine endp(n,x,l,u,m,lambda,equatn,linear)
      implicit none
      integer m,n
      logical equatn(m),linear(m)
      double precision l(n),lambda(m),u(n),x(n)
      end subroutine endp

      subroutine algencan(epsfeas,epsopt,iprint,ncomp,n,x,l,u,m,lambda, &
                     equatn,linear,coded,checkder,f,cnorm,snorm,nlpsupn,inform)
      use message_passing
      use global
      logical checkder
      integer inform,iprint,m,n,ncomp
      double precision cnorm,epsfeas,epsopt,f,nlpsupn,snorm
      logical coded(10),equatn(m),linear(m)
      double precision l(n),lambda(m),u(n),x(n)

! dummy routine to indicate the real routine must be obtained

      ierr = USER_INPUT_ERROR
      call fatal("HP_NLP requires a third party optimization routine that is not included with the PHAML distribution.", &
                 "Contact phaml@nist.gov for instructions on obtaining and configuring it.")
      stop
      end subroutine algencan
