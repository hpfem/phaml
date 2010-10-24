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

module grid_util

!----------------------------------------------------
! This module contains grid utility routines.
! communication tags in this module are of the form 35xx
!----------------------------------------------------

!----------------------------------------------------
! Other modules used:

use global
use hash_mod
use gridtype_mod
use quadrature_rules
use basis_functions
use message_passing
!----------------------------------------------------

implicit none
private
public allocate_grid, &        ! not thread safe
       deallocate_grid, &      ! not thread safe
       realloc_solution, &     ! not thread safe; TEMP OPENMP can be
       reset_dirich_exact, &   ! not thread safe; TEMP OPENMP can be
       copy_grid, &            ! not thread safe
       edge_exact, &           ! thread safe if vertex not changing
       elem_exact, &           ! thread safe if vertex and edge not changing
       point_on_edge, &        ! thread safe
       level2_mate, &          ! thread safe
       get_edge_elements, &    ! thread safe
       extend_nlev, &          ! not thread safe
       more_verts, &           ! not thread safe
       more_edges, &           ! not thread safe
       more_elements, &        ! not thread safe
       phier2nodal, &          ! thread safe
       nodal2phier, &          ! thread safe
       nodal2hhier, &          ! thread safe
       element_diameter, &     ! thread safe
       list_elements, &        ! not thread safe
       list_edges_without_rule, & ! not thread safe
       count_dof, &            ! not thread safe
       compute_global_max_errind ! not thread safe

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

end interface

!----------------------------------------------------
! The following defined constants are defined:

!----------------------------------------------------
! The following variables are defined:

!----------------------------------------------------

contains

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

!call count_grid_memory(grid)

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
               grid%dof = grid%dof - (grid%edge(edge)%degree-1)*system_size
               if (grid%element(grid%edge(edge)%assoc_elem)%iown) then
                  grid%dof_own = grid%dof_own - &
                                 (grid%edge(edge)%degree-1)*system_size
               endif
            endif
            grid%edge(edge)%degree = degree
            grid%dof = grid%dof + (grid%edge(edge)%degree-1)*system_size
            if (grid%element(grid%edge(edge)%assoc_elem)%iown) then
               grid%dof_own = grid%dof_own + &
                              (grid%edge(edge)%degree-1)*system_size
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
            grid%dof = grid%dof - system_size * &
                 ((grid%element(elem)%degree-2)*(grid%element(elem)%degree-1))/2
            if (grid%element(elem)%iown) then
               grid%dof_own = grid%dof_own - system_size * &
                 ((grid%element(elem)%degree-2)*(grid%element(elem)%degree-1))/2
            endif
         endif
         grid%element(elem)%degree = degree
         grid%dof = grid%dof + system_size * &
                 ((grid%element(elem)%degree-2)*(grid%element(elem)%degree-1))/2
         if (grid%element(elem)%iown) then
            grid%dof_own = grid%dof_own + system_size * &
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

!          ----------
subroutine edge_exact(grid,edge,sysrank,what)
!          ----------

!----------------------------------------------------
! This routine sets "exact" solutions on an edge.  what indicates what is
! being set, and can be "d" for Dirichlet boundary conditions, "i" for
! initial conditions, or "t" for the true solution.
! It assumes the coefficients for the linear elements at the endpoints of
! edge are already set, and computes the coefficients of the high order
! elements along that edge as a least squares fit to the function.
! For multiple eigenvalue problems, all eigenfunctions are set.  It is
! assumed that they all satisfy the same boundary conditions (bconds is
! defined for component=1..system_size) but have different initial conditions
! and true solutions (iconds and trues are defined for component=1..system_size
! and eigen=1..num_eval)
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: edge, sysrank
character(len=1), intent(in) :: what
!----------------------------------------------------
! Local variables:

integer :: i, j, nqp, nev, jerr, deg(4), astat, n, ss
real(my_real) :: xv(3),yv(3)
real(my_real), pointer :: weight(:), xq(:), yq(:)
real(my_real), allocatable :: basis(:,:), a(:), b(:,:), true(:,:)
integer, allocatable :: ipiv(:)
real(my_real) :: c(grid%system_size,grid%system_size),rs(grid%system_size)
integer :: itype(grid%system_size)

!----------------------------------------------------
! Non-module procedures used are:

interface

   function trues(x,y,comp,eigen) ! real (my_real)
   use global
   real (my_real), intent(in) :: x,y
   integer, intent(in) :: comp,eigen
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
   integer, intent(in) :: comp, eigen
   real(my_real) :: iconds
   end function iconds

end interface

!----------------------------------------------------
! Begin executable code

! nothing to do if only linear bases on this edge

if (grid%edge(edge)%degree < 2) return

! Dirichlet boundary conditions override initial conditions

if (what == "i" .and. grid%edge_type(edge,sysrank) == DIRICHLET) return

! useful constants

nev = max(1,grid%num_eval)
ss = grid%system_size

! local copy of the vertex coordinates, plus a fake third vertex to make
! a triangle

xv(1) = grid%vertex(grid%edge(edge)%vertex(1))%coord%x
xv(2) = grid%vertex(grid%edge(edge)%vertex(2))%coord%x
xv(3) = (xv(1)+xv(2))/2
yv(1) = grid%vertex(grid%edge(edge)%vertex(1))%coord%y
yv(2) = grid%vertex(grid%edge(edge)%vertex(2))%coord%y
yv(3) = (yv(1)+yv(2))/2
if (abs(xv(2)-xv(1)) > abs(yv(2)-yv(1))) then
   yv(3) = yv(3) + 1
else
   xv(3) = xv(3) + 1
endif

! get the quadrature rule for the edge

call quadrature_rule_line(min(MAX_QUAD_ORDER_LINE,grid%edge(edge)%degree+1), &
                          xv(1:2), yv(1:2), nqp, weight, xq, yq, jerr)
if (jerr /= 0) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("Error getting quadrature rule in edge_exact",intlist=(/jerr/))
   stop
endif

! evaluate basis functions at the quadrature points

allocate(basis(3+grid%edge(edge)%degree-1,nqp),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in edge_exact")
   stop
endif

deg = 1
deg(3) = grid%edge(edge)%degree

call p_hier_basis_func(xq,yq,xv,yv,deg,"a",basis)

! evaluate the function at the quadrature points

allocate(true(nqp,nev),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in edge_exact")
   stop
endif

select case (what)
case ("d")
   do i=1,nqp
      call bconds(xq(i),yq(i),grid%edge(edge)%bmark,itype,c,rs)
      true(i,:) = rs(sysrank)
   end do
case ("i")
   do j=1,nev
      do i=1,nqp
         true(i,j) = iconds(xq(i),yq(i),sysrank,j)
      end do
   end do
case ("t")
   do j=1,nev
      do i=1,nqp
         true(i,j) = trues(xq(i),yq(i),sysrank,j)
         if (true(i,j) == huge(0.0_my_real)) true(i,j) = 0.0_my_real
      end do
   end do
end select

! take off the linear bases part of the function

do j=1,nev
 if (what == "t") then
  if (associated(grid%vertex_exact)) then
   true(:,j) = true(:,j) - &
    grid%vertex_exact(grid%edge(edge)%vertex(1),sysrank,j)*basis(1,:)-&
    grid%vertex_exact(grid%edge(edge)%vertex(2),sysrank,j)*basis(2,:)
  endif
 else
  true(:,j) = true(:,j) - &
   grid%vertex_solution(grid%edge(edge)%vertex(1),sysrank,j)*basis(1,:) - &
   grid%vertex_solution(grid%edge(edge)%vertex(2),sysrank,j)*basis(2,:)
 endif
end do

! set up least squares linear system, storing the matrix in LAPACK symmetric
! packed form

n = grid%edge(edge)%degree - 1

allocate(a((n*(n+1))/2),b(n,nev),ipiv(n),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in edge_exact")
   stop
endif

do i=1,n
   do j=i,n
      a(i+(j*(j-1))/2) = sum(weight*basis(3+i,:)*basis(3+j,:))
   end do
   do j=1,nev
      b(i,j) = sum(weight*basis(3+i,:)*true(:,j))
   end do
end do

! solve the least squares system

if (my_real == kind(1.0e0)) then
   call sspsv("U",n,nev,a,ipiv,b,n,jerr)
elseif (my_real == kind(1.0d0)) then
   call dspsv("U",n,nev,a,ipiv,b,n,jerr)
else
   ierr = PHAML_INTERNAL_ERROR
      call fatal("my_real is neither single nor double precision. Can't call LAPACK")
   stop 
endif

! copy the solution to the edge data structure

do j=1,nev
   if (what == "t") then
      if (grid%have_true) then
         do i=1,n
            grid%edge(edge)%exact(i,sysrank,j) = b(i,j)
         end do
      endif
   else
      do i=1,n
         grid%edge(edge)%solution(i,sysrank,j) = b(i,j)
      end do
   endif
end do

deallocate(weight,xq,yq,basis,a,b,true,ipiv)

end subroutine edge_exact

!          ----------
subroutine elem_exact(grid,elem,sysrank,what)
!          ----------

!----------------------------------------------------
! This routine sets "exact" solutions for face basis functions.
! what indicates what is being set, and can be "i" for initial conditions
! or "t" for the true solution.
! It assumes the coefficients for the linear elements at the vertices of
! elem and edge basis functions are already set, and computes the
! coefficients of the high order face bases in elem as a least squares
! fit to the function.
! For multiple eigenvalue problems, all eigenfunctions are set.  It is
! assumed that they have different initial conditions
! and true solutions (iconds and trues are defined for component=1..system_size
! and eigen=1..num_eval)
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: elem, sysrank
character(len=1), intent(in) :: what
!----------------------------------------------------
! Local variables:

integer :: i, j, k, nqp, nev, ss, jerr, deg(4), astat, n, degree, nbasis, isub
real(my_real) :: xv(3),yv(3)
real(my_real), pointer :: weight(:), xq(:), yq(:)
real(my_real), allocatable :: basis(:,:), wbasis(:,:), a(:), b(:,:), true(:,:),&
                              afull(:,:)
integer, allocatable :: ipiv(:)

!----------------------------------------------------
! Non-module procedures used are:

interface

   function trues(x,y,comp,eigen) ! real (my_real)
   use global
   real (my_real), intent(in) :: x,y
   integer, intent(in) :: comp,eigen
   real (my_real) :: trues
   end function trues

   function iconds(x,y,comp,eigen)
   use global
   real(my_real), intent(in) :: x,y
   integer, intent(in) :: comp, eigen
   real(my_real) :: iconds
   end function iconds

end interface

!----------------------------------------------------
! Begin executable code

! must be at least cubics to have face bases

if (grid%element(elem)%degree < 3) return

! useful constants

nev = max(1,grid%num_eval)
ss = grid%system_size

! local copy of the vertex coordinates

xv = grid%vertex(grid%element(elem)%vertex)%coord%x
yv = grid%vertex(grid%element(elem)%vertex)%coord%y

! get the quadrature rule for the triangle

degree = grid%element(elem)%degree
do i=1,EDGES_PER_ELEMENT
   degree = max(degree,grid%edge(grid%element(elem)%edge(i))%degree)
end do
if (degree > 1) then
   degree = min(MAX_QUAD_ORDER_TRI,2*degree-2)
endif
call quadrature_rule_tri(degree,xv,yv,nqp,weight,xq,yq,jerr,stay_in=.true.)
if (jerr /= 0) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("Error getting quadrature rule in elem_exact",intlist=(/jerr/))
   stop
endif

! evaluate basis functions at the quadrature points

nbasis = 3
do i=1,EDGES_PER_ELEMENT 
   deg(i) = grid%edge(grid%element(elem)%edge(i))%degree
   nbasis = nbasis + max(0,deg(i)-1)
end do
deg(4) = grid%element(elem)%degree
nbasis = nbasis + ((deg(4)-1)*(deg(4)-2))/2
allocate(basis(nbasis,nqp),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in elem_exact")
   stop
endif

call p_hier_basis_func(xq,yq,xv,yv,deg,"a",basis)

! evaluate the function at the quadrature points

allocate(true(nqp,nev),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in elem_exact")
   stop
endif

do j=1,nev
   do i=1,nqp
      if (what == "t") then
         true(i,j) = trues(xq(i),yq(i),sysrank,j)
         if (true(i,j) == huge(0.0_my_real)) true(i,j) = 0.0_my_real
      else
         true(i,j) = iconds(xq(i),yq(i),sysrank,j)
      endif
   end do
end do

! take off the vertex- and edge-basis parts of the function

do i=1,nev
 if (what == "t") then
  if (associated(grid%vertex_exact)) then
   true(:,i) = true(:,i) - &
        grid%vertex_exact(grid%element(elem)%vertex(1),sysrank,i)*basis(1,:) - &
        grid%vertex_exact(grid%element(elem)%vertex(2),sysrank,i)*basis(2,:) - &
        grid%vertex_exact(grid%element(elem)%vertex(3),sysrank,i)*basis(3,:)
   isub = 3
   do j=1,EDGES_PER_ELEMENT
      do k=1,grid%edge(grid%element(elem)%edge(j))%degree-1
         isub = isub + 1
         if (grid%have_true) then
            true(:,i) = true(:,i) - &
               grid%edge(grid%element(elem)%edge(j))%exact(k,sysrank,i)*basis(isub,:)
         endif
      end do
   end do
  endif
 else
   true(:,i) = true(:,i) - &
      grid%vertex_solution(grid%element(elem)%vertex(1),sysrank,i)*basis(1,:) -&
      grid%vertex_solution(grid%element(elem)%vertex(2),sysrank,i)*basis(2,:) -&
      grid%vertex_solution(grid%element(elem)%vertex(3),sysrank,i)*basis(3,:)
   isub = 3
   do j=1,EDGES_PER_ELEMENT
      do k=1,grid%edge(grid%element(elem)%edge(j))%degree-1
        isub = isub + 1
        true(:,i) = true(:,i) - &
         grid%edge(grid%element(elem)%edge(j))%solution(k,sysrank,i)*basis(isub,:)
      end do
   end do
 endif
end do

! set up least squares linear system, storing the matrix in LAPACK symmetric
! packed form

n = ((grid%element(elem)%degree-1)*(grid%element(elem)%degree-2))/2

allocate(a((n*(n+1))/2),b(n,nev),ipiv(n),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in elem_exact")
   stop
endif

isub = 3
do i=1,EDGES_PER_ELEMENT 
   isub = isub + max(0,grid%edge(grid%element(elem)%edge(i))%degree-1)
end do

allocate(afull(n,n),wbasis(nqp,n),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in elem_exact")
   stop
endif
afull = 0.0_my_real
basis(1:n,:) = basis(isub+1:nbasis,:)
do i=1,nqp
   wbasis(i,1:n) = weight(i)*basis(1:n,i)
end do

if (my_real == kind(1.0)) then
   call sgemm("N","N",n,n,nqp,1.0_my_real,basis,nbasis,wbasis,nqp, &
              0.0_my_real,afull,n)
elseif (my_real == kind(1.0d0)) then
   call dgemm("N","N",n,n,nqp,1.0_my_real,basis,nbasis,wbasis,nqp, &
              0.0_my_real,afull,n)
else
   ierr = PHAML_INTERNAL_ERROR
   call fatal("my_real is neither single nor double precision. Can't call GEMM")
   stop 
endif

do i=1,n
   do j=i,n
      a(i+(j*(j-1))/2) = afull(i,j)
   end do
end do

if (my_real == kind(1.0)) then
   call sgemm("T","N",n,nev,nqp,1.0_my_real,wbasis,nqp,true,nqp,0.0_my_real,b,n)
elseif (my_real == kind(1.0d0)) then
   call dgemm("T","N",n,nev,nqp,1.0_my_real,wbasis,nqp,true,nqp,0.0_my_real,b,n)
endif

! solve the least squares system

if (my_real == kind(1.0e0)) then
   call sspsv("U",n,nev,a,ipiv,b,n,jerr)
elseif (my_real == kind(1.0d0)) then
   call dspsv("U",n,nev,a,ipiv,b,n,jerr)
else
   ierr = PHAML_INTERNAL_ERROR
      call fatal("my_real is neither single nor double precision. Can't call LAPACK")
   stop
endif

! copy the solution to the element data structure

do j=1,nev
   if (what == "t") then
      if (grid%have_true) then
         do i=1,n
            grid%element(elem)%exact(i,sysrank,j) = b(i,j)
         end do
      endif
   else
      do i=1,n
         grid%element(elem)%solution(i,sysrank,j) = b(i,j)
      end do
   endif
end do

deallocate(weight,xq,yq,basis,wbasis,a,b,true,ipiv,afull)

end subroutine elem_exact

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

!        -----------------
function get_edge_elements(grid,edge)
!        -----------------

!----------------------------------------------------
! This routine returns the two elements that share edge.  If the edge is
! on the boundary, the second one will be BOUNDARY.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: edge
integer, dimension(2) :: get_edge_elements
!----------------------------------------------------
! Local variables:

integer :: elem1, side, i, neigh(EDGES_PER_ELEMENT)
!----------------------------------------------------
! Begin executable code

! one of them is the associated element

elem1 = grid%edge(edge)%assoc_elem

! determine which edge of the first element it is

side = 0
do i=1,EDGES_PER_ELEMENT
   if (grid%element(elem1)%edge(i) == edge) then
      side = i
      exit
   endif
end do

if (side == 0) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("get_edge_elements: failed to find edge on associated element")
   stop
endif

! get the neighbor elements

neigh = get_neighbors(elem1,grid)

! the i'th neighbor shares the i'th side

get_edge_elements(1) = elem1
get_edge_elements(2) = neigh(side)

end function get_edge_elements

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

!call count_grid_memory(grid)

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

!call count_grid_memory(grid)

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

!call count_grid_memory(grid)

end subroutine more_verts

!          -----------------
subroutine count_grid_memory(grid)
!          -----------------

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

end subroutine count_grid_memory

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
   call sgemm("N","N",nbasis4,size(solution,dim=2),nbasis4,1.0_my_real, &
              s,nbasis4,solution,nbasis4,0.0,prod,nbasis4)
elseif (my_real == kind(1.0d0)) then
   call dgemm("N","N",nbasis4,size(solution,dim=2),nbasis4,1.0_my_real, &
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

!          -------------
subroutine list_elements(grid,list,nelem,level,own,leaf)
!          -------------

!----------------------------------------------------
! This routine creates a list of elements in the beginning of list.
! nelem is the number of returned elements.  level is the h-refinement
! level to use; if level=0, all levels are used.  If own is true, only
! owned elements are listed.  If leaf is true, only leaves are listed.
! list must be large enough to hold the list.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(out) :: list(:), nelem
integer, intent(in) :: level
logical, intent(in) :: own, leaf
!----------------------------------------------------
! Local variables:

integer :: lolim, hilim, lev, elem
!----------------------------------------------------
! Begin executable code

if (level == 0) then
   lolim = 1
   hilim = grid%nlev
else
   lolim = level
   hilim = level
endif

nelem = 0
do lev=lolim,hilim
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (((.not. own) .or. grid%element(elem)%iown) .and. &
          ((.not. leaf) .or. grid%element(elem)%isleaf)) then
         nelem = nelem + 1
         list(nelem) = elem
      endif
      elem = grid%element(elem)%next
   end do
end do

end subroutine list_elements

!          -----------------------
subroutine list_edges_without_rule(grid,refine_control,element_list,nelem, &
                                   edge_list,nedge)
!          -----------------------

!----------------------------------------------------
! This routine makes a list of the edges of elements in element_list(1:nelem)
! that do not satisfy the edge rule.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
type(refine_options), intent(in) :: refine_control
integer, intent(in) :: element_list(:), nelem
integer, intent(out) :: edge_list(:), nedge
!----------------------------------------------------
! Local variables:

integer :: i, j, edge_elements(2), elem, edge
logical :: on_list(size(grid%edge))
!----------------------------------------------------
! Begin executable code

on_list = .false.

nedge = 0
do i=1,nelem
   elem = element_list(i)
   do j=1,3
      edge = grid%element(elem)%edge(j)
      if (on_list(edge)) cycle
      edge_elements = get_edge_elements(grid,edge)
      if (edge_elements(2) == BOUNDARY) then
         if (grid%edge(edge)%degree /= &
             grid%element(edge_elements(1))%degree) then
            nedge = nedge + 1
            edge_list(nedge) = edge
            on_list(edge) = .true.
         endif
      elseif (refine_control%edge_rule == MINIMUM_RULE) then
         if (grid%edge(edge)%degree /= min( &
                grid%element(edge_elements(1))%degree, &
                grid%element(edge_elements(2))%degree)) then
            nedge = nedge + 1
            edge_list(nedge) = edge
            on_list(edge) = .true.
         endif
      else ! (refine_control%edge_rule == MAXIMUM_RULE)
         if (grid%edge(edge)%degree /= max( &
                grid%element(edge_elements(1))%degree, &
                grid%element(edge_elements(2))%degree)) then
            nedge = nedge + 1
            edge_list(nedge) = edge
            on_list(edge) = .true.
         endif
      endif
   end do
end do

end subroutine list_edges_without_rule

!        ---------
function count_dof(grid,procs,just_owned)
!        ---------

!----------------------------------------------------
! This routine counts the degrees of freedom in grid.  If procs is present,
! it returns the dof of the global grid.  If just_owned is present and true, it
! returns the owned dof of this processor's grid.  Otherwise it returns all
! the dof of this processor's grid.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
type(proc_info), intent(in), optional :: procs
logical, intent(in), optional :: just_owned
integer :: count_dof
!----------------------------------------------------
! Local variables:

logical :: visited(size(grid%edge)), only_owned
integer :: lev, vert, elem, deg, side, edge
!----------------------------------------------------
! Begin executable code

! can't have both procs and just_owned

if (present(procs) .and. present(just_owned)) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("Can't give both procs and just_owned in count_dof.")
   stop
endif

! determine if the count is of all dof or only the owned dof

if (present(procs)) then
   only_owned = .true.
elseif (present(just_owned)) then
   only_owned = just_owned
else
   only_owned = .false.
endif

! count the vertices; each has system_size dofs

count_dof = 0
do lev=1,grid%nlev
   vert = grid%head_level_vert(lev)
   do while (vert /= END_OF_LIST)
      if (grid%element(grid%vertex(vert)%assoc_elem)%iown .or. &
          .not. only_owned) then
         count_dof = count_dof + count(grid%vertex_type(vert,:)/=PERIODIC_SLAVE)
      endif
      vert = grid%vertex(vert)%next
   end do
end do

! go through the leaf elements counting the degree-dependent dofs and counting
! the dofs on the edges that have not yet been visited

visited = .false.
do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (grid%element(elem)%isleaf .and. &
          (grid%element(elem)%iown .or. .not. only_owned)) then
         deg = grid%element(elem)%degree
         if (deg >= 3) then
            count_dof = count_dof + grid%system_size*((deg-1)*(deg-2))/2
         endif
         do side=1,EDGES_PER_ELEMENT
            edge = grid%element(elem)%edge(side)
            if (.not. visited(edge) .and. &
                (grid%element(grid%edge(edge)%assoc_elem)%iown .or. &
                 .not. only_owned)) then
               count_dof = count_dof + (grid%edge(edge)%degree-1)* &
                           count(grid%edge_type(edge,:)/=PERIODIC_SLAVE)
               visited(edge) = .true.
            endif
         end do
      endif
      elem = grid%element(elem)%next
   end do
end do

! sum over all processors

if (present(procs)) then
   count_dof = phaml_global_sum(procs,count_dof,3501)
endif

end function count_dof

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

compute_global_max_errind = phaml_global_max(procs,compute_global_max_errind, &
                                             3540)

end function compute_global_max_errind

end module grid_util
