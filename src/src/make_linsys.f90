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

module make_linsys

!----------------------------------------------------
! This module contains routines to create and destroy the linear system
!
! communication tags in this module are of the form 11xx
!----------------------------------------------------

!----------------------------------------------------
! Other modules used are:

use global
use message_passing
use stopwatch
use hash_mod
use hash_eq_mod
use linsystype_mod
use gridtype_mod
use linsys_util
use quadrature_rules
use basis_functions
use evaluate

!----------------------------------------------------

implicit none
private
public create_linear_system, destroy_linear_system, edge_exact, elem_exact, &
       elemental_matrix

!----------------------------------------------------
! Non-module procedures used are:

interface

   subroutine bconds(x,y,bmark,itype,c,rs)
   use global
   real(my_real), intent(in) :: x,y
   integer, intent(in) :: bmark
   integer, intent(out) :: itype(:)
   real(my_real), intent(out) :: c(:,:),rs(:)
   end subroutine bconds

   subroutine pdecoefs(x,y,cxx,cxy,cyy,cx,cy,c,rs)
   use global
   real(my_real), intent(in) :: x,y
   real(my_real), intent(out) :: cxx(:,:),cxy(:,:),cyy(:,:),cx(:,:),cy(:,:), &
                                 c(:,:),rs(:)
   end subroutine pdecoefs

end interface

!----------------------------------------------------

contains

!          --------------------
subroutine create_linear_system(linear_system,grid,procs,solver_cntl, &
                                io_cntl,still_sequential,notime)
!          --------------------

!----------------------------------------------------
! This routine creates the linear system.
!
! RESTRICTION 1) if a vertex in the initial grid is a peak, it is a peak of all
!             its triangles in the initial grid
!             2) in the linked list passage through the vertices of the initial
!             grid, non-peaks come before peaks
!             These are required for h-hierarchical basis changes, but
!             not necessary otherwise.  init_grid from Triangle data files
!             satisfies this.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(out) :: linear_system
type(grid_type), intent(in) :: grid
type(proc_info), intent(in) :: procs
type(solver_options), intent(in) :: solver_cntl
type(io_options), intent(in) :: io_cntl
logical, intent(in) :: still_sequential
logical, intent(in), optional :: notime

!----------------------------------------------------
! Local variables:

logical :: timeit
integer :: nelem_leaf, maxeq_per_elem, allocstat
integer, allocatable :: eqlist(:), leaf_element(:)
real(my_real) :: lambda0
real(my_real), allocatable :: tempvec(:)
integer :: i, j

!----------------------------------------------------
! Begin executable code

! MASTER doesn't have a linear system, but needs some things set

if (my_proc(procs) == MASTER) then
   linear_system%neq      = 0
   linear_system%neq_vert = 0
   linear_system%neq_edge = 0
   linear_system%neq_face = 0
   nullify(linear_system%begin_row)
   nullify(linear_system%nn_comm_end_of_send)
   return
endif

! Start timing the assembly process.

if (present(notime)) then
   timeit = .not. notime
else
   timeit = .true.
endif

if (timeit) then
   call reset_watch((/passemble,cpassemble/))
   call start_watch((/passemble,tassemble/))
endif

! Initialize minimum value of r(x,y).

linear_system%rmin = huge(0.0_my_real)

! Create a list of all the leaf elements.

allocate(leaf_element(grid%nelem_leaf),stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in create_linear_system",procs=procs)
   return
endif
call make_leaf_list(grid,leaf_element,nelem_leaf)

! Set the easy scalars.

linear_system%coarse_band_exists = .false.
linear_system%lapack_symm_band_exists = .false.
linear_system%lapack_gen_band_exists = .false.
linear_system%petsc_matrix_exists = .false.
linear_system%mumps_matrix_exists = .false.
linear_system%superlu_matrix_exists = .false.
linear_system%hypre_matrix_exists = .false.

linear_system%system_size = solver_cntl%system_size
!linear_system%nlev = grid%nlev
linear_system%nlev = phaml_global_max(procs,grid%nlev,1101)
linear_system%maxdeg = 1

! Determine number of vertex, edge and face equations and maximum number of
! equations related to any element.

call determine_neq(grid,linear_system,leaf_element,nelem_leaf, &
                   maxeq_per_elem)
linear_system%maxdeg = phaml_global_max(procs,linear_system%maxdeg,1100)

! Create the hash table.

call hash_table_init(linear_system%eq_hash,linear_system%neq)

! Allocate space for linear_system.

call allocate_linsys(linear_system,solver_cntl,procs)

! Create the matrices for converting linear basis functions between
! the nodal basis and the h-hierarchical basis.

call create_basis_change(linear_system)

! Determine the order of the equations.  In the process, set begin_level,
! equation_type, iown, and equation gids, copy solution from the grid, and
! insert equations into the hash table.

call equation_order(linear_system,grid,leaf_element,nelem_leaf)

! Create the column_index array.  This allocates column_index, and sets
! column_index, begin_row, end_row, end_row_linear, and end_row_edge

allocate(eqlist(maxeq_per_elem),stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in create_linear_system",procs=procs)
   return
endif
call create_column_index(linear_system,grid,procs,eqlist,leaf_element, &
                         nelem_leaf)

! Allocate the matrix value arrays, and initialize them and the right hand
! side arrays to 0.0

call allocate_matrix(linear_system,solver_cntl,io_cntl,procs)

! Compute the matrix and right hand side.

call compute_matrix_rhs(linear_system,grid,procs,solver_cntl,maxeq_per_elem, &
                        leaf_element,nelem_leaf,eqlist,timeit,still_sequential)

deallocate(eqlist, leaf_element, stat=allocstat)
if (allocstat /= 0) then
   call warning("deallocation failed in create_linear_system")
endif

! Determine the linear system communication map, i.e., which equation ids
! need to be communicated under various situations

call make_communication_map(linear_system,grid,procs,still_sequential,timeit)

! TEMP080114 charged particles
! Undocumented feature for computing the electric field of charged particles
if (charged_particles) then
   call charged_particles_rhs(linear_system,grid)
endif

! Undocumented feature for using Crank-Nicholson to solve the time-dependent
! Schroedinger equation.  Subtract A*u_old from the right hand side.

if (crank_nicholson) then
   allocate(tempvec(linear_system%neq))
   linear_system%matrix_val => linear_system%stiffness
   call matrix_times_vector(linear_system%solution(1:),tempvec,linear_system, &
                            procs,still_sequential,10001,10002,10003,10004, &
                            10005,10006)
   linear_system%rhs = linear_system%rhs - tempvec
   where (linear_system%equation_type == DIRICHLET) linear_system%rhs = 0.0_my_real
   deallocate(tempvec)
endif

! For BLOPEX, set Dirichlet columns and rows to the identity in both stiffness
! and mass.  Note that this restricts Dirichlet boundary conditions to be
! homogeneous, but I believe that is required for well-posedness of the
! problem anyway.

if (solver_cntl%eq_type == EIGENVALUE .and. &
    solver_cntl%eigensolver == BLOPEX_SOLVER) then
   do i=1,linear_system%neq
      do j=linear_system%begin_row(i),linear_system%end_row(i)
         if (linear_system%column_index(j) == NO_ENTRY) cycle
         if (linear_system%equation_type(i) == DIRICHLET .or. &
             linear_system%equation_type(linear_system%column_index(j)) == DIRICHLET) then
            if (linear_system%column_index(j) == i) then
               linear_system%stiffness(j) = 1.0_my_real
               linear_system%mass(j) = 1.0_my_real
            else
               linear_system%stiffness(j) = 0.0_my_real
               linear_system%mass(j) = 0.0_my_real
            endif
         endif
      end do
   end do
endif

! If we are using BLOPEX and computing eigenvalues to the right of lambda0,
! then we solve -Ax = (-lambda)Mx transformed to
! M (-A-lambda0M)^(-1) M x = 1/(-lambda-lambda0) M x
! so negate the stiffness matrix and lambda0.

if (solver_cntl%eq_type == EIGENVALUE .and. &
    solver_cntl%eigensolver == BLOPEX_SOLVER .and. &
    solver_cntl%lambda0_side == EIGEN_RIGHT .and. &
    solver_cntl%lambda0 /= -huge(0.0_my_real)) then
   linear_system%stiffness = - linear_system%stiffness
endif

! For eigenvalue problems, set the shifted matrix A-lambda0*M.
! For BLOPEX, if lambda0 = -inf, don't shift.

if (solver_cntl%eq_type == EIGENVALUE) then
   lambda0 = solver_cntl%lambda0
   if (lambda0 == -huge(0.0_my_real)) then
      if (solver_cntl%eigensolver == BLOPEX_SOLVER) then
         linear_system%rmin = 0.0_my_real
         lambda0 = 0.0_my_real
      else
         linear_system%rmin = phaml_global_min(procs,linear_system%rmin,1120)
         lambda0 = linear_system%rmin
      endif
   else
      linear_system%rmin = 0.0_my_real
   endif
   if (solver_cntl%eigensolver == BLOPEX_SOLVER .and. &
       solver_cntl%lambda0_side == EIGEN_RIGHT) then
      linear_system%shifted = linear_system%stiffness + lambda0*linear_system%mass
   else
      linear_system%shifted = linear_system%stiffness - lambda0*linear_system%mass
   endif
else
   linear_system%rmin = 0.0_my_real
endif

! If the residual is to be computed while solving a boundary value problem
! (print_error_when is FREQUENTLY or TOO_MUCH) keep the stiffness by copying
! it to condensed.

if ((io_cntl%print_error_when == FREQUENTLY .or. &
     io_cntl%print_error_when == TOO_MUCH) .and. &
     solver_cntl%eq_type /= EIGENVALUE) then
   linear_system%condensed = linear_system%stiffness
   linear_system%rhs_cond = linear_system%rhs_nocond
endif

! Perform static condensation to eliminate the face bases

call static_condensation(linear_system,grid,procs)

! stop timing the assembly process

if (timeit) then
   call stop_watch((/passemble,tassemble/))
endif

end subroutine create_linear_system

!          --------------
subroutine make_leaf_list(grid,leaf_element,nelem)
!          --------------

!----------------------------------------------------
! This routine makes a list of the leaf elements
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(inout) :: leaf_element(:), nelem

!----------------------------------------------------
! Local variables:

integer :: elem

!----------------------------------------------------
! Begin executable code

nelem = 0
elem = grid%head_level_elem(1)
do while (elem /= END_OF_LIST)
   call make_leaf_list1(grid,leaf_element,nelem,elem)
   elem = grid%element(elem)%next
end do
if (nelem /= grid%nelem_leaf) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("number of leaf elements found is not equal to number of leaf elements in grid%nelem_leaf in make_leaf_list", &
              intlist=(/nelem,grid%nelem_leaf/))
   return
endif

end subroutine make_leaf_list

!                    ---------------
recursive subroutine make_leaf_list1(grid,leaf_element,nelem,elem)
!                    ---------------

!----------------------------------------------------
! This routine makes a list of the leaf elements below element elem
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(inout) :: leaf_element(:), nelem
integer, intent(in) :: elem

!----------------------------------------------------
! Local variables:

integer :: i
integer :: child(MAX_CHILD), allc(MAX_CHILD)

!----------------------------------------------------
! Begin executable code

allc = ALL_CHILDREN
child = get_child_lid(grid%element(elem)%gid,allc,grid%elem_hash)
if (child(1) /= NO_CHILD) then

! for non-leaves, move on to the children

   do i=1,MAX_CHILD
      if (child(i) /= NO_CHILD) then
         call make_leaf_list1(grid,leaf_element,nelem,child(i))
      endif
   end do

else

! add leaves to the list

   nelem = nelem + 1
   leaf_element(nelem) = elem

endif

end subroutine make_leaf_list1

!          -------------
subroutine determine_neq(grid,linear_system,leaf_element,nelem_leaf, &
                         maxeq_per_elem)
!          -------------

!----------------------------------------------------
! This routine determines the number of vertex, edge and face equations
! and maximum number of equations related to any element
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(inout) :: linear_system
type(grid_type), intent(in) :: grid
integer, intent(in) :: leaf_element(:), nelem_leaf
integer, intent(out) :: maxeq_per_elem
!----------------------------------------------------
! Local variables:

logical :: visited(size(grid%edge))
integer :: i, j, k, elem, edge, vert, degree, eq_this_elem, lev, syssize, &
           neqvert, neqedge, neqface
!----------------------------------------------------
! Begin executable code

syssize = linear_system%system_size
visited = .false.
neqvert = grid%nvert
neqedge = 0
neqface = 0
maxeq_per_elem = 0
do j=1,nelem_leaf
   elem = leaf_element(j)
   eq_this_elem = VERTICES_PER_ELEMENT
   do i=1,EDGES_PER_ELEMENT
      edge = grid%element(elem)%edge(i)
      degree = grid%edge(edge)%degree
      if (degree > 1) then
         eq_this_elem = eq_this_elem + degree-1
         if (visited(edge)) cycle
         neqedge = neqedge + degree - 1
         visited(edge) = .true.
      endif
   end do
   degree = grid%element(elem)%degree
   linear_system%maxdeg = max(degree,linear_system%maxdeg)
   if (degree > 2) then
      eq_this_elem = eq_this_elem + ((degree-1)*(degree-2))/2
      neqface = neqface + ((degree-1)*(degree-2))/2
   endif
   maxeq_per_elem = max(maxeq_per_elem,eq_this_elem)
end do
neqvert = neqvert*syssize
neqedge = neqedge*syssize
neqface = neqface*syssize

! remove slave periodic entities from equation count

do lev=1,grid%nlev
   vert = grid%head_level_vert(lev)
   do while (vert /= END_OF_LIST)
      do j=1,syssize
         if (grid%vertex_type(vert,j) == PERIODIC_SLAVE .or. &
             grid%vertex_type(vert,j) == PERIODIC_SLAVE_DIR .or. &
             grid%vertex_type(vert,j) == PERIODIC_SLAVE_NAT .or. &
             grid%vertex_type(vert,j) == PERIODIC_SLAVE_MIX) neqvert=neqvert-1
      end do
      vert = grid%vertex(vert)%next
   end do
end do
do k=1,nelem_leaf
   elem = leaf_element(k)
   do i=1,EDGES_PER_ELEMENT
      edge = grid%element(elem)%edge(i)
      if (grid%edge(edge)%degree <= 1) cycle
      do j=1,syssize
         if (grid%edge_type(edge,j) == PERIODIC_SLAVE .or. &
             grid%edge_type(edge,j) == PERIODIC_SLAVE_DIR .or. &
             grid%edge_type(edge,j) == PERIODIC_SLAVE_NAT .or. &
             grid%edge_type(edge,j) == PERIODIC_SLAVE_MIX) then
            neqedge = neqedge - (grid%edge(edge)%degree-1)
         endif
      end do
   end do
end do

linear_system%neq = neqvert + neqedge + neqface
linear_system%neq_vert = neqvert
linear_system%neq_edge = neqedge
linear_system%neq_face = neqface

end subroutine determine_neq

!          ---------------
subroutine allocate_linsys(linear_system,solver_cntl,procs)
!          ---------------

!----------------------------------------------------
! This routine allocates memory in linear_system
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(inout) :: linear_system
type(solver_options), intent(in) :: solver_cntl
type(proc_info), intent(in) :: procs
!----------------------------------------------------
! Local variables:

integer :: neq, syssize, allocstat
!----------------------------------------------------
! Begin executable code

neq = linear_system%neq
syssize = linear_system%system_size

! allocate space based on neq, nlev, etc.
! Any changes to this list should be matched in the deallocate statement
! in destroy_linear_system

! TEMP it might be possible to save space by not storing the face equations
!      (just the factorized blocks) if not solving an eigenvalue problem

! TEMP some of these can probably omit face equations, and maybe even edge
!      equations, to reduce memory requirements

allocate(linear_system%begin_row(neq+1), &
         linear_system%end_row_linear(neq+1), &
         linear_system%end_row_edge(neq+1), &
         linear_system%end_row_face(neq+1), &
         linear_system%begin_level(linear_system%nlev+3), &
         linear_system%rhs_nocond(neq), &
         linear_system%r_mine(neq), &
         linear_system%r_others(neq), &
         linear_system%need_r_others(neq), &
         linear_system%iown(neq), &
         linear_system%solution(0:neq), &
         linear_system%equation_type(neq), &
         linear_system%gid(neq), &
         linear_system%hold_dirich(neq), &
         linear_system%s_int(5*syssize,5*syssize), &
         linear_system%s_int_inv(5*syssize,5*syssize), &
         linear_system%s_bnd(4*syssize,4*syssize), &
         linear_system%s_bnd_inv(4*syssize,4*syssize), &
         stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in create_linear_system",procs=procs)
   return
endif
linear_system%end_row => linear_system%end_row_face
linear_system%rhs => linear_system%rhs_nocond
if (solver_cntl%num_eval > 1) then
   allocate(linear_system%evecs(neq,solver_cntl%num_eval-1),stat=allocstat)
   if (allocstat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in create_linear_system",procs=procs)
      return
   endif
else
   nullify(linear_system%evecs)
endif
nullify(linear_system%matrix_val,linear_system%column_index, &
        linear_system%stiffness,linear_system%mass, &
        linear_system%shifted)

end subroutine allocate_linsys

!          -------------------
subroutine create_basis_change(linear_system)
!          -------------------

!----------------------------------------------------
! This routine creates the matrices for converting linear basis functions
! between ! the nodal basis and the h-hierarchical basis
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(inout) :: linear_system
!----------------------------------------------------
! Local variables:

integer :: i, syssize
!----------------------------------------------------
! Begin executable code

syssize = linear_system%system_size

linear_system%s_int = 0.0_my_real
linear_system%s_int_inv = 0.0_my_real
linear_system%s_bnd = 0.0_my_real
linear_system%s_bnd_inv = 0.0_my_real

do i=1,5*syssize
   linear_system%s_int(i,i) = 1.0_my_real
   linear_system%s_int_inv(i,i) = 1.0_my_real
end do
do i=1,4*syssize
   linear_system%s_bnd(i,i) = 1.0_my_real
   linear_system%s_bnd_inv(i,i) = 1.0_my_real
end do

do i=1,syssize
   linear_system%s_int(4*syssize+i,i) = 0.5_my_real
   linear_system%s_int(4*syssize+i,i+syssize) = 0.5_my_real
   linear_system%s_int_inv(4*syssize+i,i) = -0.5_my_real
   linear_system%s_int_inv(4*syssize+i,i+syssize) = -0.5_my_real
   linear_system%s_bnd(3*syssize+i,i) = 0.5_my_real
   linear_system%s_bnd(3*syssize+i,i+syssize) = 0.5_my_real
   linear_system%s_bnd_inv(3*syssize+i,i) = -0.5_my_real
   linear_system%s_bnd_inv(3*syssize+i,i+syssize) = -0.5_my_real
end do

end subroutine create_basis_change

!          --------------
subroutine equation_order(linear_system,grid,leaf_element,nelem_leaf)
!          --------------

!----------------------------------------------------
! This routine orders the equations with vertex bases (linear bases) first,
! ordered by level, and then edge bases, and then face bases.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(inout) :: linear_system
type(grid_type), intent(in) :: grid
integer, intent(in) :: leaf_element(:), nelem_leaf
!----------------------------------------------------
! Local variables:

logical :: visited(size(grid%edge))
integer :: i, j, k, l, eq, lev, vert, edge, elem, degree, syssize
!----------------------------------------------------
! Begin executable code

! First pass through the vertices to determine the linear-basis equation order,
! global IDs, initial solution, equation type and ownership, and the start of
! equations for each level

syssize = linear_system%system_size
linear_system%solution(0) = 0.0_my_real
eq = 1
do lev=1,grid%nlev
   linear_system%begin_level(lev) = eq
   vert = grid%head_level_vert(lev)
   do while (vert /= END_OF_LIST)
      do j=1,syssize
         if (grid%vertex_type(vert,j) == PERIODIC_SLAVE .or. &
             grid%vertex_type(vert,j) == PERIODIC_SLAVE_DIR .or. &
             grid%vertex_type(vert,j) == PERIODIC_SLAVE_NAT .or. &
             grid%vertex_type(vert,j) == PERIODIC_SLAVE_MIX) cycle
         call grid_to_eq(grid,linear_system,VERTEX_ID,1,j, &
                        grid%vertex(vert)%gid,eqn_gid=linear_system%gid(eq))
         call hash_insert(linear_system%gid(eq),eq,linear_system%eq_hash)
         linear_system%solution(eq) = grid%vertex_solution(vert,j,1)
         select case (grid%vertex_type(vert,j))
         case (PERIODIC_MASTER_DIR)
            linear_system%equation_type(eq) = DIRICHLET
         case (PERIODIC_MASTER_NAT)
            linear_system%equation_type(eq) = NATURAL
         case (PERIODIC_MASTER_MIX)
            linear_system%equation_type(eq) = MIXED
         case default
            linear_system%equation_type(eq) = grid%vertex_type(vert,j)
         end select
         linear_system%iown(eq) =  grid%element(grid%vertex(vert)%assoc_elem)%iown
         eq = eq + 1
      end do
      vert = grid%vertex(vert)%next
   end do
end do

do lev=grid%nlev+1,linear_system%nlev
   linear_system%begin_level(lev) = eq
end do

if (eq-1 /= linear_system%neq_vert) then
   call warning("number of vertex equations created is not equal to neqvert")
endif
linear_system%begin_level(linear_system%nlev+1) = eq

! Next pass through the elements setting the same components of the
! edge basis equations.  Note the edge bases all fall between
! begin_level(nlev+1) and begin_level(nlev+2)-1.

visited = .false.
do l=1,nelem_leaf
   elem = leaf_element(l)
   do i=1,EDGES_PER_ELEMENT
      edge = grid%element(elem)%edge(i)
      if (visited(edge)) cycle
      visited(edge) = .true.
      degree = grid%edge(edge)%degree
      if (degree < 2) cycle
      do k=1,degree-1
         do j=1,syssize
            if (grid%edge_type(edge,j) == PERIODIC_SLAVE) cycle
            call grid_to_eq(grid,linear_system,EDGE_ID,k,j, &
                          grid%edge(edge)%gid,eqn_gid=linear_system%gid(eq))
            call hash_insert(linear_system%gid(eq),eq,linear_system%eq_hash)
            linear_system%solution(eq) = grid%edge(edge)%solution(k,j,1)
            linear_system%equation_type(eq) = grid%edge_type(edge,j)
            linear_system%iown(eq) = &
               grid%element(grid%edge(edge)%assoc_elem)%iown
            eq = eq + 1
         end do
      end do
   end do
end do

if (eq-1 /= linear_system%neq_vert + linear_system%neq_edge) then
   call warning("number of edge equations is not equal to neqedge")
endif
linear_system%begin_level(linear_system%nlev+2) = eq

! Then pass through the elements setting the same components of the
! face basis equations.  Note the face bases all fall between
! begin_level(nlev+2) and begin_level(nlev+3)-1.

do l=1,nelem_leaf
   elem = leaf_element(l)
   degree = grid%element(elem)%degree
   if (degree > 2) then
      do k=1,((degree-1)*(degree-2))/2
         do j=1,syssize
            call grid_to_eq(grid,linear_system,ELEMENT_ID,k,j, &
                 grid%element(elem)%gid,eqn_gid=linear_system%gid(eq+j-1))
            call hash_insert(linear_system%gid(eq+j-1),eq+j-1, &
                             linear_system%eq_hash)
         end do
         linear_system%solution(eq:eq+syssize-1) = &
                                     grid%element(elem)%solution(k,1:syssize,1)
         linear_system%equation_type(eq:eq+syssize-1) = INTERIOR
         linear_system%iown(eq:eq+syssize-1) = grid%element(elem)%iown
         eq = eq + syssize
      end do
   endif
end do

if (eq-1 /= linear_system%neq_vert + linear_system%neq_edge + &
            linear_system%neq_face) then
   call warning("number of face equations is not equal to neqface")
endif
linear_system%begin_level(linear_system%nlev+3) = eq

end subroutine equation_order

!          -------------------
subroutine create_column_index(linear_system,grid,procs,eqlist, &
                               leaf_element,nelem_leaf)
!          -------------------

!----------------------------------------------------
! This routine creates column_index and related arrays.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(inout) :: linear_system
type(grid_type), intent(in) :: grid
type(proc_info), intent(in) :: procs
integer, intent(out) :: eqlist(:)
integer, intent(in) :: leaf_element(:), nelem_leaf
!----------------------------------------------------
! Local variables:

logical :: visited(size(grid%edge))
integer :: i, j, k, l, m, lev, vert, edge, elem, eq, eq1, neq, nentry, degree, &
           indx, syssize, objtype, brank, srank, allocstat, max_entries_in_row
integer, allocatable :: adjacencies(:,:)
!----------------------------------------------------
! Begin executable code

neq = linear_system%neq
syssize = linear_system%system_size

! Determine the maximum number of entries in any row of the matrix

max_entries_in_row = find_max_entries_in_row(grid)

! Pass through the elements to determine which bases have
! overlapping support.  For each equation build a temporary list
! of equations from this pass, and then create the column_index array

! create and initialize the lists

allocate(adjacencies(max_entries_in_row,neq),stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation for adjacencies failed",procs=procs)
   return
endif
! starts with itself, for the diagonal
adjacencies = 0
adjacencies(1,:) = (/ (eq,eq=1,neq) /)

! Pass through the elements building lists and counting entries.
! The list only contains the first equation of each block if system_size > 1.

nentry = neq/syssize
do i=1,nelem_leaf
   eq = 1
   do j=1,VERTICES_PER_ELEMENT
      call grid_to_eq(grid,linear_system,VERTEX_ID,1,1, &
          grid%vertex(grid%element(leaf_element(i))%vertex(j))%gid,eqlist(eq))
      eq = eq + 1
   end do
   do j=1,EDGES_PER_ELEMENT
      do brank = 1,grid%edge(grid%element(leaf_element(i))%edge(j))%degree-1
         call grid_to_eq(grid,linear_system,EDGE_ID,brank,1, &
               grid%edge(grid%element(leaf_element(i))%edge(j))%gid,eqlist(eq))
         eq = eq + 1
      end do
   end do
   degree = grid%element(leaf_element(i))%degree
   do j=1,((degree-1)*(degree-2))/2
      call grid_to_eq(grid,linear_system,ELEMENT_ID,j,1, &
                      grid%element(leaf_element(i))%gid,eqlist(eq))
      eq = eq + 1
   end do
   call add_adjacencies(eqlist,eq-1,adjacencies,nentry)
end do

! Allocate column_index.
! Any changes to this list should be matched in the deallocate statement
! in destroy_linear_system

allocate(linear_system%column_index(nentry*syssize*syssize), stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation for matrix failed",procs=procs)
   return
endif

! Build column_index, begin_row, end_row_linear, end_row_edge and end_row.

! first the equations from vertices, level by level

indx = 1
do lev=1,grid%nlev
   vert = grid%head_level_vert(lev)
   do while (vert /= END_OF_LIST)
      call grid_to_eq(grid,linear_system,VERTEX_ID,1,1, &
                      grid%vertex(vert)%gid,eq1)
      eq = eq1 - 1
      do i=1,syssize
         if (grid%vertex_type(vert,i) == PERIODIC_SLAVE .or. &
             grid%vertex_type(vert,i) == PERIODIC_SLAVE_DIR .or. &
             grid%vertex_type(vert,i) == PERIODIC_SLAVE_NAT .or. &
             grid%vertex_type(vert,i) == PERIODIC_SLAVE_MIX) cycle
         eq = eq + 1
         linear_system%begin_row(eq) = indx
         do m=1,max_entries_in_row
            if (adjacencies(m,eq1) == 0) exit
            linear_system%column_index(indx:indx+syssize-1) = &
                                       (/(adjacencies(m,eq1)+j,j=0,syssize-1)/)
! for the diagonal block (which comes first in the row storage), swap so that
! the diagonal entry is the first entry of the row
            if (indx < linear_system%begin_row(eq) + syssize) then
              linear_system%column_index(indx) = adjacencies(m,eq1)+i-1
              linear_system%column_index(indx+i-1) = adjacencies(m,eq1)
            endif
! if this block is from a linear basis, set end_row_linear so it is the last
! linear basis column in the end
            call eq_to_grid(linear_system,linear_system%gid(adjacencies(m,eq1)), &
                            objtype,brank,srank)
            if (objtype == VERTEX_ID) then
               linear_system%end_row_linear(eq) = indx + syssize - 1
            endif
            if (objtype == VERTEX_ID .or. objtype == EDGE_ID) then
               linear_system%end_row_edge(eq) = indx + syssize - 1
            endif
            indx = indx + syssize
         end do
         linear_system%end_row(eq) = indx - 1
      end do
      vert = grid%vertex(vert)%next
   end do
end do

! then the same for the edge bases

visited = .false.
do l=1,nelem_leaf
   elem = leaf_element(l)
   do j=1,EDGES_PER_ELEMENT
      edge = grid%element(elem)%edge(j)
      if (visited(edge)) cycle
      visited(edge) = .true.
      degree = grid%edge(edge)%degree
      if (degree < 2) cycle
      do k=1,degree-1
         call grid_to_eq(grid,linear_system,EDGE_ID,k,1,grid%edge(edge)%gid,eq1)
         eq = eq1 - 1
         do i=1,syssize
            if (grid%edge_type(edge,i) == PERIODIC_SLAVE) cycle
            eq = eq + 1
            linear_system%begin_row(eq) = indx
            do m=1,max_entries_in_row
               if (adjacencies(m,eq1) == 0) exit
               linear_system%column_index(indx:indx+syssize-1) = &
                  (/(adjacencies(m,eq1)+j,j=0,syssize-1)/)
               if (indx < linear_system%begin_row(eq) + syssize) then
                  linear_system%column_index(indx) = adjacencies(m,eq1)+i-1
                  linear_system%column_index(indx+i-1) = adjacencies(m,eq1)
               endif
               call eq_to_grid(linear_system, &
                              linear_system%gid(adjacencies(m,eq1)), &
                              objtype,brank,srank)
               if (objtype == VERTEX_ID) then
                  linear_system%end_row_linear(eq) = indx + syssize - 1
               endif
               if (objtype == VERTEX_ID .or. objtype == EDGE_ID) then
                  linear_system%end_row_edge(eq) = indx + syssize - 1
               endif
               indx = indx + syssize
            end do ! next column (block)
            linear_system%end_row(eq) = indx - 1
         end do ! next PDE
      end do ! next basis
   end do ! next edge
end do ! next element

! finally the same for the face bases

do l=1,nelem_leaf
   elem = leaf_element(l)
   degree = grid%element(elem)%degree
   if (degree < 3) cycle
   do k=1,((degree-1)*(degree-2))/2
      call grid_to_eq(grid,linear_system,ELEMENT_ID,k,1,grid%element(elem)%gid,&
                      eq1)
      eq = eq1 - 1
      do i=1,syssize
         eq = eq + 1
         linear_system%begin_row(eq) = indx
         do m=1,max_entries_in_row
            if (adjacencies(m,eq1) == 0) exit
            linear_system%column_index(indx:indx+syssize-1) = &
                  (/(adjacencies(m,eq1)+j,j=0,syssize-1)/)
            if (indx < linear_system%begin_row(eq) + syssize) then
               linear_system%column_index(indx) = adjacencies(m,eq1)+i-1
               linear_system%column_index(indx+i-1) = adjacencies(m,eq1)
            endif
            call eq_to_grid(linear_system, &
                           linear_system%gid(adjacencies(m,eq1)), &
                           objtype,brank,srank)
            if (objtype == VERTEX_ID) then
               linear_system%end_row_linear(eq) = indx + syssize - 1
            endif
            if (objtype == VERTEX_ID .or. objtype == EDGE_ID) then
               linear_system%end_row_edge(eq) = indx + syssize - 1
            endif
            indx = indx + syssize
         end do ! next column (block)
         linear_system%end_row(eq) = indx - 1
      end do ! next PDE
   end do ! next basis
end do ! next element

linear_system%begin_row(neq+1) = indx
linear_system%end_row(neq+1) = indx

if (indx /= nentry*syssize*syssize+1) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("miscounted number of entries in the matrix.", &
              intlist=(/indx,nentry*syssize*syssize+1/),procs=procs)
   return
endif

! free up the memory used by adjacencies

deallocate(adjacencies,stat=allocstat)
if (allocstat /= 0) then
   call warning("deallocation of adjacencies failed in create_linear_system")
endif

end subroutine create_column_index

!        -----------------------
function find_max_entries_in_row(grid)
!        -----------------------

!----------------------------------------------------
! This routine determines the maximum number of entries in any row when
! system_size is 1 (otherwise multiply by system_size).
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer :: find_max_entries_in_row
!----------------------------------------------------
! Local variables:

integer :: num_adjacent(size(grid%vertex))
integer :: lev, elem, i, j, deg
integer :: v(VERTICES_PER_ELEMENT), neigh(EDGES_PER_ELEMENT)
!----------------------------------------------------
! Begin executable code

! The maximum will occur at a row corresponding to a vertex, so add up the
! number of overlapping bases for each vertex basis.  Pass through the elements
! and for each vertex of the element add 1) twice the number of face bases,
! 2) the number of edge bases on adjacent edges, or twice the number if the
! edge is on the boundary, 3) twice the number of edge bases on the opposite
! edge, and 4) 1 for each adjacent vertex, or 2 if the vertex is on the
! boundary.  Then divide by 2.  And add 1 for the diagonal.
! The business with 2 takes care of duplications over adjacent elements.

num_adjacent = 0

! for each element
do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (grid%element(elem)%isleaf) then
         v = grid%element(elem)%vertex
         neigh = get_neighbors(elem,grid)
! for each vertex of this element
         do i=1,VERTICES_PER_ELEMENT
            deg = grid%element(elem)%degree
! face bases
            num_adjacent(v(i)) = num_adjacent(v(i)) + (deg-1)*(deg-2)
            do j=1,EDGES_PER_ELEMENT
               deg = grid%edge(grid%element(elem)%edge(j))%degree
               if (j==i) then
! opposite edge
                  num_adjacent(v(i)) = num_adjacent(v(i)) + 2*(deg-1)
               else
! adjacent edges
                  if (neigh(j)==BOUNDARY) then
                     num_adjacent(v(i)) = num_adjacent(v(i)) + 2*(deg-1)
                  else
                     num_adjacent(v(i)) = num_adjacent(v(i)) + deg-1
                  endif
! adjacent vertices
                  if (neigh(6-(i+j))==BOUNDARY) then
                     num_adjacent(v(i)) = num_adjacent(v(i)) + 2
                  else
                     num_adjacent(v(i)) = num_adjacent(v(i)) + 1
                  endif
               endif
            end do
         end do
      endif
      elem = grid%element(elem)%next
   end do
end do

find_max_entries_in_row = maxval(num_adjacent)/2 + 1

end function find_max_entries_in_row

!          ---------------
subroutine add_adjacencies(list,n,adjacencies,nentry)
!          ---------------

!----------------------------------------------------
! This routine takes a list of equations and makes sure all pairs are on the
! lists of adjacencies
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: list(:), n
integer, intent(inout) :: adjacencies(:,:)
integer, intent(inout) :: nentry

!----------------------------------------------------
! Local variables:

integer :: i, j, list_eq, adj_eq, entry, temp
!----------------------------------------------------
! Begin executable code

! For each ordered pair of list entries

do i=1,n
   list_eq = list(i)
   do j=1,n
      adj_eq = list(j)
      if (list_eq == adj_eq) cycle ! skip diagonal

! look for adjacency with larger node number in list

      entry = 2
      do
! if 0, reached end; add it to the end
         if (adjacencies(entry,list_eq) == 0) then
            adjacencies(entry,list_eq) = adj_eq
            nentry = nentry + 1
            exit
         endif
! if equal, already on list
         if (adjacencies(entry,list_eq) == adj_eq) exit
! if greater, insert here and move remaining entries up one
         if (adjacencies(entry,list_eq) > adj_eq) then
            do
               temp = adjacencies(entry,list_eq)
               if (temp == 0) exit
               adjacencies(entry,list_eq) = adj_eq
               adj_eq = temp
               entry = entry + 1
            end do
            adjacencies(entry,list_eq) = adj_eq
            nentry = nentry + 1
            exit
         endif
! if less, move on to next
         entry = entry + 1
      end do

   end do
end do

end subroutine add_adjacencies

!          ------
subroutine make_communication_map(linear_system,grid,procs,still_sequential, &
                                  timeit)
!          ------

!----------------------------------------------------
! This routine creates the communication map for linear_system, i.e, list
! of equation ids that need to be communicated in various situations
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(inout) :: linear_system
type(grid_type), intent(in) :: grid
type(proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential, timeit
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

call make_nn_comm_map(linear_system,grid,procs,still_sequential,timeit)
call make_fudop_comm_map(linear_system,grid,procs,still_sequential,timeit)

end subroutine make_communication_map

!          ----------------
subroutine make_nn_comm_map(linear_system,grid,procs,still_sequential,timeit)
!          ----------------

!----------------------------------------------------
! This routine creates the nearest neighbor communication map for
! linear_system, i.e, list of equation ids that are adjacent to owned equation
! IDS that need to be communicated
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(inout) :: linear_system
type(grid_type), intent(in) :: grid
type(proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential, timeit
!----------------------------------------------------
! Local variables:

integer :: nproc, eq, i, counter, ind, proc, lid, object_type, basis_rank, &
           system_rank, deg, nint, astat
integer, allocatable :: send_int(:), nsend(:), nrecv(:)
integer, pointer :: recv_int(:)
logical :: owned_neigh
type(hash_key_eq) :: gid
!----------------------------------------------------
! Begin executable code

nproc = num_proc(procs)

nullify(linear_system%nn_comm_end_of_send)
allocate(linear_system%nn_comm_remote_neigh(linear_system%neq))
linear_system%nn_comm_remote_neigh = .false.

! cases where this is not necessary

if (my_proc(procs) == MASTER .or. still_sequential .or. nproc==1) return

! make lists of nearest neighbor data

! make a list of all equations that I don't own but are adjacent to one I
! do own, i.e., unowned rows with at least one owned column.

! count the equations

counter = 0
do eq=1,linear_system%neq
   if (.not. linear_system%iown(eq)) then
      owned_neigh = .false.
      do i=linear_system%begin_row(eq)+1,linear_system%end_row(eq)
         if (linear_system%column_index(i) == NO_ENTRY) cycle
         if (linear_system%iown(linear_system%column_index(i))) then
            owned_neigh = .true.
            exit
         endif
      end do
      if (owned_neigh) counter = counter + KEY_SIZE+1
   endif
end do

allocate(send_int(counter),nsend(nproc),nrecv(nproc),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in make_nn_comm_map",procs=procs)
   return 
endif
   
! Make a list of the equations I don't own but are a neighbor of one I do own.
! Also flag the equations that have a neighbor I don't own.

counter = 0 
do eq=1,linear_system%neq
   if (.not. linear_system%iown(eq)) then
      owned_neigh = .false.
      do i=linear_system%begin_row(eq)+1,linear_system%end_row(eq)
         if (linear_system%column_index(i) == NO_ENTRY) cycle
         if (linear_system%iown(linear_system%column_index(i))) then
            owned_neigh = .true.
            linear_system%nn_comm_remote_neigh(linear_system%column_index(i)) = .true.
         endif
      end do
      if (owned_neigh) then
         call hash_pack_key(linear_system%gid(eq),send_int,counter+1)
         counter = counter + KEY_SIZE+1
      endif
   endif
end do

! perform all_to_all exchange of the unowned equation IDs

if (timeit) call start_watch((/cpassemble,ctassemble/))
call phaml_alltoall(procs,send_int,counter,recv_int,nrecv,1140)
if (timeit) call stop_watch((/cpassemble,ctassemble/))
deallocate(send_int)

! for each processor, create a list of global IDs that I own to return to
! the processor, and a list of local IDs to use when data is exchanged.  The
! lists of IDs for vertices and edges are ordered by degree so we can send
! messages containing data for bases up to a given degree, and then followed
! by IDs for faces so they can be included when sending for the full
! uncondensed matrix.

! count the number I have of each degree for each processor

allocate(linear_system%nn_comm_end_of_send(nproc,linear_system%maxdeg+1))
linear_system%nn_comm_end_of_send = 0
ind = 1
do proc=1,nproc
   do eq=1,nrecv(proc)/(KEY_SIZE+1)
      gid = hash_unpack_key(recv_int,ind,extended=.true.)
      lid = hash_decode_key(gid,linear_system%eq_hash)
      if (lid /= HASH_NOT_FOUND) then
         if (linear_system%iown(lid)) then
            call eq_to_grid(linear_system,gid,object_type,basis_rank, &
                            system_rank)
            if (object_type == VERTEX_ID) then
               deg = 1
            elseif (object_type == EDGE_ID) then
               deg = basis_rank+1
            else ! ELEMENT_ID
               deg = linear_system%maxdeg+1
            endif
            linear_system%nn_comm_end_of_send(proc,deg) = &
                 linear_system%nn_comm_end_of_send(proc,deg) + 1
         endif
      endif
      ind = ind + KEY_SIZE+1
   end do
end do

! set the end of lists of each degree

do proc=1,nproc
   do deg=2,linear_system%maxdeg+1
      linear_system%nn_comm_end_of_send(proc,deg) = &
                 linear_system%nn_comm_end_of_send(proc,deg) + &
                 linear_system%nn_comm_end_of_send(proc,deg-1)
   end do
end do

! allocate the lists

allocate(linear_system%nn_comm_send_lid(nproc))
nint = nproc
do proc=1,nproc
   allocate(linear_system%nn_comm_send_lid(proc)%lid(linear_system%nn_comm_end_of_send(proc,linear_system%maxdeg+1)))
   nint = nint + linear_system%nn_comm_end_of_send(proc,linear_system%maxdeg+1)*(KEY_SIZE+1)
end do
allocate(send_int(nint))

! create the lists to keep

! reset the end of list to the beginning and increment as we go

do proc=1,nproc
   do deg = linear_system%maxdeg+1,2,-1
      linear_system%nn_comm_end_of_send(proc,deg) = linear_system%nn_comm_end_of_send(proc,deg-1)
   end do
   linear_system%nn_comm_end_of_send(proc,1) = 0
end do

ind = 1
do proc=1,nproc
   do eq=1,nrecv(proc)/(KEY_SIZE+1)
      gid = hash_unpack_key(recv_int,ind,extended=.true.)
      lid = hash_decode_key(gid,linear_system%eq_hash)
      if (lid /= HASH_NOT_FOUND) then
         if (linear_system%iown(lid)) then
            call eq_to_grid(linear_system,gid,object_type,basis_rank, &
                            system_rank)
            if (object_type == VERTEX_ID) then
               deg = 1
            elseif (object_type == EDGE_ID) then
               deg = basis_rank+1
            else ! ELEMENT_ID
               deg = linear_system%maxdeg+1
            endif
            linear_system%nn_comm_end_of_send(proc,deg) = &
                       linear_system%nn_comm_end_of_send(proc,deg) + 1
            linear_system%nn_comm_send_lid(proc)%lid(linear_system%nn_comm_end_of_send(proc,deg)) = lid
         endif
      endif
      ind = ind + KEY_SIZE+1
   end do
end do

! create the list to send

ind = 1
do proc=1,nproc
   send_int(ind) = linear_system%nn_comm_end_of_send(proc,linear_system%maxdeg+1)
   ind = ind + 1
   do eq = 1,linear_system%nn_comm_end_of_send(proc,linear_system%maxdeg+1)
      gid = linear_system%gid(linear_system%nn_comm_send_lid(proc)%lid(eq))
      call hash_pack_key(gid,send_int,ind)
      ind = ind + KEY_SIZE+1
   end do
   nsend(proc) = 1 + linear_system%nn_comm_end_of_send(proc,linear_system%maxdeg+1)*(KEY_SIZE+1)
end do

! exchange the ID lists

deallocate(recv_int)
if (timeit) call start_watch((/cpassemble,ctassemble/))
call phaml_alltoall(procs,send_int,nsend,recv_int,nrecv,1141)
if (timeit) call stop_watch((/cpassemble,ctassemble/))

! for each processor, create a list of local IDs of data I will be receiving
! from that processor

allocate(linear_system%nn_comm_recv_lid(nproc))
ind = 1
do proc=1,nproc
   allocate(linear_system%nn_comm_recv_lid(proc)%lid(recv_int(ind)))
   ind = ind + 1
   do eq=1,nrecv(proc)/(KEY_SIZE+1)
      gid = hash_unpack_key(recv_int,ind,extended=.true.)
      ind = ind + KEY_SIZE+1
      lid = hash_decode_key(gid,linear_system%eq_hash)
      linear_system%nn_comm_recv_lid(proc)%lid(eq) = lid
   end do
end do

! free memory

deallocate(send_int,recv_int,nsend,nrecv)

! The following components of linear system were allocated in this routine
! and will be deallocated in destroy_linear_system.  If others are added in
! here, be sure to destroy them, too.

!   nn_comm_end_of_send
!   nn_comm_send_lid
!   nn_comm_send_lid()%lid
!   nn_comm_recv_lid
!   nn_comm_recv_lid()%lid

end subroutine make_nn_comm_map

!          -------------------
subroutine make_fudop_comm_map(linear_system,grid,procs,still_sequential,timeit)
!          -------------------

!----------------------------------------------------
! This routine creates the fudop communication map for linear_system, i.e, list
! of equation ids that are unowned and need to be communicated, and equation
! ids for which r_others needs to be communicated
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(inout) :: linear_system
type(grid_type), intent(in) :: grid
type(proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential, timeit
!----------------------------------------------------
! Local variables:

type(hash_key_eq) :: gid
integer :: nproc, my_processor, i, j, k, p, q, neigh, counter, eq, astat, ind, &
           proc, lid, nint, step, nunowned, nrothers
integer, allocatable :: send_int(:), nrecv(:), power_map(:,:), old_power(:,:), &
                        nsend(:), all_proxy(:,:), location(:,:), &
                        mess_size(:,:,:), nchunk(:,:)
integer, pointer :: nn_map(:), recv_int(:), all_sizes(:)
logical(small_logical), allocatable :: need_r_others(:)
logical :: do_resid
!----------------------------------------------------
! Begin executable code

nproc = num_proc(procs)
my_processor = my_proc(procs)
nullify(linear_system%fudop_comm_proxy)

! cases where this is not necessary

if (my_proc(procs) == MASTER .or. still_sequential .or. nproc==1) return

! get the full nearest neighbor map of all processors

allocate(send_int(nproc),nrecv(nproc))
do p=1,nproc
   send_int(p) = min(1,linear_system%nn_comm_end_of_send(p,linear_system%maxdeg+1))
end do
if (timeit) call start_watch((/cpassemble,ctassemble/))
call phaml_alltoall(procs,send_int,nproc,nn_map,nrecv,1150)
if (timeit) call stop_watch((/cpassemble,ctassemble/))
deallocate(send_int)

! determine the proxy for each processor pair, i.e. a nearest neighbor to which
! a message is sent to be delivered to the other processor

allocate(all_proxy(nproc,nproc))
all_proxy = -1

! first, nearest neighbors are their own proxy

do p=1,nproc
   all_proxy(p,p) = p
   do q=1,nproc
      if (nn_map(nproc*(p-1)+q) /= 0) then
         all_proxy(p,q) = q
      endif
   end do
end do

if (any(all_proxy==-1)) then

! By raising the map to the k'th power, we get a map of processors that are
! within a distance of k via nearest neighbor connections.  Start with the
! nearest neighbor map.

   allocate(power_map(nproc,nproc),old_power(nproc,nproc))
   power_map = reshape(nn_map,(/nproc,nproc/))

! iterate until all proxies have been determined

   do

! for each processor pair that does not yet have a proxy

      do p=1,nproc
         do q=1,nproc
            if (all_proxy(p,q) /= -1) cycle

! check each of p's neighbors to see if any of them knows where to send it

            do neigh=1,nproc
               if (nn_map(nproc*(p-1)+neigh) == 0) cycle
               if (power_map(neigh,q) /= 0) then
                  all_proxy(p,q) = neigh
                  exit
               endif
            end do
         end do
      end do

! see if all have been assigned

      if (all(all_proxy/=-1)) exit

! next power of the nearest neighbor map

      old_power = power_map
      do i=1,nproc
         do j=1,nproc
            power_map(i,j) = 0
            do k=1,nproc
               power_map(i,j) = power_map(i,j) + nn_map(nproc*(i-1)+k)*old_power(k,j)
            end do
         end do
      end do

   end do ! until all proxies have been determined
   deallocate(power_map,old_power)
endif

deallocate(nn_map)

! keep my proxies

allocate(linear_system%fudop_comm_proxy(nproc))
linear_system%fudop_comm_proxy = all_proxy(my_processor,:)

! determine which equations need r_others

allocate(need_r_others(linear_system%neq))
call make_need_r_others(linear_system)
! TEMP090212 eventually remove need_r_others from the linsys data structure
! and return the result directly in need_r_others
need_r_others = linear_system%need_r_others

! make lists of fudop data

! make a list of all equations that I don't own followed by a list of all
! equations for which I need r_others

! count the equations

counter = 0
do eq=1,linear_system%neq
   if (.not. linear_system%iown(eq)) then
      counter = counter + KEY_SIZE+1
   endif
   if (need_r_others(eq)) then
      counter = counter + KEY_SIZE+1
   endif
end do
counter = counter + KEY_SIZE+1

allocate(send_int(counter),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in make_fudop_comm_map",procs=procs)
   return
endif

! make the list

counter = 0
do eq=1,linear_system%neq
   if (.not. linear_system%iown(eq)) then
      call hash_pack_key(linear_system%gid(eq),send_int,counter+1)
      counter = counter + KEY_SIZE+1
   endif
end do
call hash_pack_key(NULL_KEY_EQ,send_int,counter+1)
counter = counter + KEY_SIZE+1
do eq=1,linear_system%neq
   if (need_r_others(eq)) then
      call hash_pack_key(linear_system%gid(eq),send_int,counter+1)
      counter = counter + KEY_SIZE+1
   endif
end do

! perform all_to_all exchange of the equation IDs

if (timeit) call start_watch((/cpassemble,ctassemble/))
call phaml_alltoall(procs,send_int,counter,recv_int,nrecv,1151)
if (timeit) call stop_watch((/cpassemble,ctassemble/))
deallocate(send_int)

! for each processor, create a list of global IDs that I own to return to
! the processor, and a list of local IDs to use when data is exchanged, for
! the unowned equations and r_others

! count the number I have for each processor

allocate(linear_system%fudop_comm_nsend(nproc), &
         linear_system%resid_comm_nsend(nproc))
linear_system%fudop_comm_nsend = 0
linear_system%resid_comm_nsend = 0
ind = 1
do proc=1,nproc
   do_resid = .false.
   do eq=1,nrecv(proc)/(KEY_SIZE+1)
      gid = hash_unpack_key(recv_int,ind,extended=.true.)
      if (gid == NULL_KEY_EQ) then
         do_resid = .true.
         ind = ind + KEY_SIZE+1
         cycle
      endif
      lid = hash_decode_key(gid,linear_system%eq_hash)
      if (lid /= HASH_NOT_FOUND) then
         if (proc /= my_processor) then
            if (do_resid) then
               linear_system%resid_comm_nsend(proc) = &
                    linear_system%resid_comm_nsend(proc) + 1
            else
               if (linear_system%iown(lid)) then
                  linear_system%fudop_comm_nsend(proc) = &
                       linear_system%fudop_comm_nsend(proc) + 1
               endif
            endif
         endif
      endif
      ind = ind + KEY_SIZE+1
   end do
end do

! allocate the lists

allocate(linear_system%fudop_comm_send_lid(nproc), &
         linear_system%resid_comm_send_lid(nproc))
nint = 2*nproc
do proc=1,nproc
   allocate(linear_system%fudop_comm_send_lid(proc)%lid(linear_system%fudop_comm_nsend(proc)), &
            linear_system%resid_comm_send_lid(proc)%lid(linear_system%resid_comm_nsend(proc)))
   nint = nint + (linear_system%fudop_comm_nsend(proc) + &
                  linear_system%resid_comm_nsend(proc))*(KEY_SIZE+1)
end do
allocate(send_int(nint),nsend(nproc))

! create the lists to keep

ind = 1
do proc=1,nproc
   counter = 0
   do_resid = .false.
   do eq=1,nrecv(proc)/(KEY_SIZE+1)
      gid = hash_unpack_key(recv_int,ind,extended=.true.)
      if (gid == NULL_KEY_EQ) then
         do_resid = .true.
         ind = ind + KEY_SIZE+1
         counter = 0
         cycle
      endif
      lid = hash_decode_key(gid,linear_system%eq_hash)
      if (lid /= HASH_NOT_FOUND) then
         if (proc /= my_processor) then
            if (do_resid) then
               counter = counter + 1
               linear_system%resid_comm_send_lid(proc)%lid(counter) = lid
            else
               if (linear_system%iown(lid)) then
                  counter = counter + 1
                  linear_system%fudop_comm_send_lid(proc)%lid(counter) = lid
               endif
            endif
         endif
      endif
      ind = ind + KEY_SIZE+1
   end do
end do

! create the list to send

ind = 1
do proc=1,nproc
   send_int(ind) = linear_system%fudop_comm_nsend(proc)
   send_int(ind+1) = linear_system%resid_comm_nsend(proc)
   ind = ind + 2
   do eq = 1,linear_system%fudop_comm_nsend(proc)
      gid = linear_system%gid(linear_system%fudop_comm_send_lid(proc)%lid(eq))
      call hash_pack_key(gid,send_int,ind)
      ind = ind + KEY_SIZE+1
   end do
   do eq = 1,linear_system%resid_comm_nsend(proc)
      gid = linear_system%gid(linear_system%resid_comm_send_lid(proc)%lid(eq))
      call hash_pack_key(gid,send_int,ind)
      ind = ind + KEY_SIZE+1
   end do
   nsend(proc) = 2 + (linear_system%fudop_comm_nsend(proc) + &
                      linear_system%resid_comm_nsend(proc))*(KEY_SIZE+1)
end do

! exchange the ID lists

deallocate(recv_int)
if (timeit) call start_watch((/cpassemble,ctassemble/))
call phaml_alltoall(procs,send_int,nsend,recv_int,nrecv,1152)
if (timeit) call stop_watch((/cpassemble,ctassemble/))

! for each processor, create a list of local IDs of data I will be receiving
! from that processor, for both unowned equations and r_others

allocate(linear_system%fudop_comm_recv_lid(nproc), &
         linear_system%resid_comm_recv_lid(nproc))
ind = 1
do proc=1,nproc
   nunowned = recv_int(ind)
   nrothers = recv_int(ind+1)
   allocate(linear_system%fudop_comm_recv_lid(proc)%lid(nunowned), &
            linear_system%resid_comm_recv_lid(proc)%lid(nrothers))
   ind = ind + 2
   do eq=1,nunowned
      gid = hash_unpack_key(recv_int,ind,extended=.true.)
      ind = ind + KEY_SIZE+1
      lid = hash_decode_key(gid,linear_system%eq_hash)
      linear_system%fudop_comm_recv_lid(proc)%lid(eq) = lid
   end do
   do eq=1,nrothers
      gid = hash_unpack_key(recv_int,ind,extended=.true.)
      ind = ind + KEY_SIZE+1
      lid = hash_decode_key(gid,linear_system%eq_hash)
      linear_system%resid_comm_recv_lid(proc)%lid(eq) = lid
   end do
end do

deallocate(send_int,recv_int,nsend)

! for each step of the fudop exchange, determine how big a message I send to
! each processor, and how many processors I will receive a message from

! get everyone's list of message lengths

call phaml_alltoall(procs,linear_system%fudop_comm_nsend,nproc,all_sizes, &
                    nrecv,1153)

! location(i,j) tells what processor contains the message from processor i
! to processor j during the current step of the exchange.  Initially the
! messages are on the senders.  0 indicates there is no message or the
! message has reached its destination.

allocate(location(nproc,nproc))
do i=1,nproc
   do j=1,nproc
      if (all_sizes(nproc*(i-1)+j) == 0) then
         location(i,j) = 0
      else
         location(i,j) = i
      endif
   end do
end do

allocate(mess_size(nproc,nproc,nproc),nchunk(nproc,nproc))
mess_size = 0
nchunk = 0

! for each step of the exchange ...

step = 1
do

! determine the message sizes and move the messages

   do i=1,nproc
      do j=1,nproc
         if (location(i,j) /= 0) then
            mess_size(location(i,j),all_proxy(location(i,j),j),step) = &
               mess_size(location(i,j),all_proxy(location(i,j),j),step) + &
               all_sizes(nproc*(i-1)+j)
            if (location(i,j)==my_processor) then
               nchunk(all_proxy(location(i,j),j),step) = &
                  nchunk(all_proxy(location(i,j),j),step) + 1
            endif
            location(i,j) = all_proxy(location(i,j),j)
            if (location(i,j) == j) location(i,j) = 0
         endif
      end do
   end do

! see if we need another step

   if (all(location == 0)) exit

   step = step + 1
end do

deallocate(all_sizes)

! keep my message sizes, how many pieces make up each message from this
! processor, and how many messages I receive in each step

allocate(linear_system%fudop_comm_mess_size(nproc,step), &
         linear_system%fudop_comm_nchunk(nproc,step), &
         linear_system%fudop_comm_nproc_recv(step))

linear_system%fudop_comm_mess_size = mess_size(my_processor,:,1:step)
linear_system%fudop_comm_nchunk = nchunk(:,1:step)
do i=1,step
   linear_system%fudop_comm_nproc_recv(i) = count(mess_size(:,my_processor,i)/=0)
end do

! repeat using message lengths that include r_others

! get everyone's list of message lengths

call phaml_alltoall(procs, &
              linear_system%fudop_comm_nsend + linear_system%resid_comm_nsend, &
                    nproc,all_sizes,nrecv,1154)

! location(i,j) tells what processor contains the message from processor i
! to processor j during the current step of the exchange.  Initially the
! messages are on the senders.  0 indicates there is no message or the
! message has reached its destination.

do i=1,nproc
   do j=1,nproc
      if (all_sizes(nproc*(i-1)+j) == 0) then
         location(i,j) = 0
      else
         location(i,j) = i
      endif
   end do
end do

mess_size = 0
nchunk = 0

! for each step of the exchange ...

step = 1
do

! determine the message sizes and move the messages

   do i=1,nproc
      do j=1,nproc
         if (location(i,j) /= 0) then
            mess_size(location(i,j),all_proxy(location(i,j),j),step) = &
               mess_size(location(i,j),all_proxy(location(i,j),j),step) + &
               all_sizes(nproc*(i-1)+j)
            if (location(i,j)==my_processor) then
               nchunk(all_proxy(location(i,j),j),step) = &
                  nchunk(all_proxy(location(i,j),j),step) + 1
            endif
            location(i,j) = all_proxy(location(i,j),j)
            if (location(i,j) == j) location(i,j) = 0
         endif
      end do
   end do

! see if we need another step

   if (all(location == 0)) exit

   step = step + 1
end do

deallocate(all_sizes,location,all_proxy,nrecv)

! keep my message sizes, how many pieces make up each message from this
! processor, and how many messages I receive in each step

allocate(linear_system%resid_comm_mess_size(nproc,step), &
         linear_system%resid_comm_nchunk(nproc,step), &
         linear_system%resid_comm_nproc_recv(step))

linear_system%resid_comm_mess_size = mess_size(my_processor,:,1:step)
linear_system%resid_comm_nchunk = nchunk(:,1:step)
do i=1,step
   linear_system%resid_comm_nproc_recv(i) = count(mess_size(:,my_processor,i)/=0)
end do

deallocate(mess_size, nchunk)

! The following components of linear system were allocated in this routine
! and will be deallocated in destroy_linear_system.  If others are added in
! here, be sure to destroy them, too.

! fudop_comm_proxy
! fudop_comm_nsend
! fudop_comm_send_lid
! fudop_comm_send_lid()%lid
! fudop_comm_recv_lid
! fudop_comm_recv_lid()%lid
! fudop_comm_mess_size
! fudop_comm_nproc_recv
! fudop_comm_nchunk
! resid_comm_nsend
! resid_comm_send_lid
! resid_comm_send_lid()%lid
! resid_comm_recv_lid
! resid_comm_recv_lid()%lid
! resid_comm_mess_size
! resid_comm_nproc_recv
! resid_comm_nchunk

end subroutine make_fudop_comm_map

!          ---------------
subroutine allocate_matrix(linear_system,solver_cntl,io_cntl,procs)
!          ---------------

!----------------------------------------------------
! This routine allocates the matrix value arrays and initializes them and
! the right hand side arrays to 0.0
! Any changes to this list should be matched in the deallocate statement
! in destroy_linear_system
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(inout) :: linear_system
type(solver_options), intent(in) :: solver_cntl
type(io_options), intent(in) :: io_cntl
type(proc_info), intent(in) :: procs
!----------------------------------------------------
! Local variables:

integer :: siz, allocstat, neq
!----------------------------------------------------
! Begin executable code

siz = size(linear_system%column_index)
neq = linear_system%neq

! Allocate stiffness matrix.

allocate(linear_system%stiffness(siz), stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation for matrix failed",procs=procs)
   return
endif
linear_system%stiffness  = 0.0_my_real

! Allocate mass matrix and shifted matrix if this is an eigenvalue problem.

if (solver_cntl%eq_type == EIGENVALUE) then
   allocate(linear_system%mass(siz),linear_system%shifted(siz), stat=allocstat)
   if (allocstat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation for matrix failed",procs=procs)
      return
   endif
   linear_system%mass = 0.0_my_real
   linear_system%shifted = 0.0_my_real
else
   nullify(linear_system%mass,linear_system%shifted)
endif

! If this is an eigenvalue problem, the matrix to be statically condensed
! is the shifted matrix; otherwise it is the stiffness matrix.  But if
! the residual is to be computed while solving a boundary value problem
! (print_error_when is FREQUENTLY or TOO_MUCH) we need to keep the stiffness
! matrix and use different space for condensed.

if (solver_cntl%eq_type == EIGENVALUE) then
   linear_system%condensed => linear_system%shifted
   linear_system%rhs_cond => linear_system%rhs_nocond
else
   if (io_cntl%print_error_when == FREQUENTLY .or. &
       io_cntl%print_error_when == TOO_MUCH) then
      allocate(linear_system%condensed(siz),linear_system%rhs_cond(neq), &
               stat=allocstat)
      if (allocstat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("memory allocation for matrix failed",procs=procs)
         return
      endif
   else
      linear_system%condensed => linear_system%stiffness
      linear_system%rhs_cond => linear_system%rhs_nocond
   endif
endif

! Initialize rhs related arrays to 0

linear_system%rhs        = 0.0_my_real
linear_system%r_mine     = 0.0_my_real
linear_system%r_others   = 0.0_my_real
linear_system%need_r_others = .false.

end subroutine allocate_matrix

!          ------------------
subroutine compute_matrix_rhs (linear_system,grid,procs,solver_cntl, &
                               maxeq_per_elem,leaf_element,nelem_leaf,eqlist, &
                               timeit,still_sequential)
!          ------------------

!----------------------------------------------------
! This routine computes the matrix and rhs values
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(inout) :: linear_system
type(grid_type), intent(in) :: grid
type(proc_info), intent(in) :: procs
type(solver_options), intent(in) :: solver_cntl
integer, intent(in) :: maxeq_per_elem,leaf_element(:),nelem_leaf
integer, intent(inout) :: eqlist(:)
logical, intent(in) :: timeit, still_sequential
!----------------------------------------------------
! Local variables:

integer :: i, siz, elem, nleq, allocstat
real(my_real), allocatable :: local_stiffness(:,:), local_mass(:,:), &
                              local_rhs(:)
!----------------------------------------------------
! Begin executable code

! space for elemental matrices

siz = maxeq_per_elem*linear_system%system_size
allocate(local_stiffness(siz,siz),local_mass(siz,siz),local_rhs(siz), &
         stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation for elemental matrix failed",procs=procs)
   return
endif

! pass through the elements

do i=1,nelem_leaf
   elem = leaf_element(i)

! compute the local matrix and right hand side

   call compute_local_matrix(grid,linear_system,solver_cntl,elem, &
                             local_stiffness,local_rhs,local_mass,eqlist,nleq)

! unless we are ignoring quadrature errors, remove the linear basis
! contributions if I don't own the element.  They will be gotten from the
! owner of the element in fix_quad_err.

   if (.not. solver_cntl%ignore_quad_err .and. &
      (.not. grid%element(elem)%iown)) then
      local_stiffness(1:3*linear_system%system_size,1:3*linear_system%system_size) = 0.0_my_real
      local_mass(1:3*linear_system%system_size,1:3*linear_system%system_size) = 0.0_my_real
      local_rhs(1:3*linear_system%system_size) = 0.0_my_real
   endif

! assemble local values into global matrix etc.

   call assemble(linear_system,eqlist,nleq,local_stiffness,local_rhs, &
                 local_mass)

end do

deallocate(local_stiffness, local_mass, local_rhs, stat=allocstat)
if (allocstat /= 0) then
   call warning("deallocation failed in create_linear_system")
endif

! fix the quadrature errors

if (.not. solver_cntl%ignore_quad_err) then
   call fix_quad_error(grid,procs,linear_system,still_sequential,timeit)
endif

end subroutine compute_matrix_rhs

!          --------------------
subroutine compute_local_matrix(grid,linear_system,solver_cntl,elem,lmat, &
                                lrhs,lmassmat,eqlist,nleq)
!          --------------------

!----------------------------------------------------
! This routine computes the local stiffness and mass matrices and right side
! over element elem.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
type(linsys_type), intent(inout) :: linear_system
type(solver_options), intent(in) :: solver_cntl
integer, intent(in) :: elem
real(my_real), intent(out) :: lmat(:,:),lrhs(:),lmassmat(:,:)
integer, intent(out) :: eqlist(:), nleq

!----------------------------------------------------
! Local variables:

real(my_real) :: xvert(VERTICES_PER_ELEMENT), yvert(VERTICES_PER_ELEMENT)
integer :: i,j,k,deg
integer :: degree(EDGES_PER_ELEMENT+1), bmark(EDGES_PER_ELEMENT), &
           edge_type(EDGES_PER_ELEMENT,linear_system%system_size)

!----------------------------------------------------
! Begin executable code

! determine local IDs of the equations

nleq = 1
do j=1,VERTICES_PER_ELEMENT
   call grid_to_eq(grid,linear_system,VERTEX_ID,1,1, &
                   grid%vertex(grid%element(elem)%vertex(j))%gid,eqlist(nleq))
   nleq = nleq + 1
end do
do j=1,EDGES_PER_ELEMENT
   do k = 1,grid%edge(grid%element(elem)%edge(j))%degree-1
      call grid_to_eq(grid,linear_system,EDGE_ID,k,1, &
                      grid%edge(grid%element(elem)%edge(j))%gid,eqlist(nleq))
      nleq = nleq + 1
   end do
end do
deg = grid%element(elem)%degree
do j=1,((deg-1)*(deg-2))/2
   call grid_to_eq(grid,linear_system,ELEMENT_ID,j,1, &
                   grid%element(elem)%gid,eqlist(nleq))
   nleq = nleq + 1
end do
nleq = nleq - 1

! get vertex coordinates

do i=1,VERTICES_PER_ELEMENT
   xvert(i) = grid%vertex(grid%element(elem)%vertex(i))%coord%x
   yvert(i) = grid%vertex(grid%element(elem)%vertex(i))%coord%y
end do

! set degree, edge_type and bmark from the element

do i=1,EDGES_PER_ELEMENT
   degree(i) = grid%edge(grid%element(elem)%edge(i))%degree
   edge_type(i,:) = grid%edge_type(grid%element(elem)%edge(i),:)
   bmark(i) = grid%edge(grid%element(elem)%edge(i))%bmark
end do
degree(EDGES_PER_ELEMENT+1) = grid%element(elem)%degree

if (associated(linear_system%mass)) then
   call elemental_matrix(grid, xvert, yvert, degree, edge_type, bmark, &
                         linear_system%system_size, solver_cntl%inc_quad_order,&
                         .false., elem, linear_system%rmin, &
                         "p","a",lmat, lrhs, loc_bconds_s=bconds, &
                         lmassmat=lmassmat)
else
   call elemental_matrix(grid, xvert, yvert, degree, edge_type, bmark, &
                         linear_system%system_size, solver_cntl%inc_quad_order,&
                         .false., elem, linear_system%rmin, &
                         "p","a",lmat, lrhs, loc_bconds_s=bconds)
endif

end subroutine compute_local_matrix

!          --------------------
subroutine elemental_matrix(grid, xvert, yvert, degree, edge_type, bmark, ss, &
                            inc_quad_order, resid_prob, elem, &
                            rmin, which_basis, which_set, lmat, lrhs, &
                            loc_bconds_a, loc_bconds_s, lmassmat, extra_rhs)
!          --------------------

!----------------------------------------------------
! This routine computes the local stiffness matrix and right hand side, and
! optionally mass matrix over the element with the given vertices and
! degrees, using the pde coefficient functions and given boundary conditions.
! extra_rhs is used for multiple right hand sides, currently supports only
! residual problems for multiple eigenvalues. Dimensions of extra_rhs are
! (system_size*nbasis,num_eval-1).  The first eigenvalue goes in lrhs and
! eigenvalues 2 to num_eval are in order in extra_rhs.
! which_basis is "n" for nodal, "h" for h-hierarchical, "p" for p-hierarchical
! which_set is "r" for red bases, "b" for black bases and "a" for all
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
real(my_real), intent(in) :: xvert(:), yvert(:)
integer, intent(in) :: degree(4), edge_type(:,:), bmark(:)
integer, intent(in) :: ss ! system_size
integer, intent(in) :: inc_quad_order, elem
logical, intent(in) :: resid_prob
real(my_real), intent(inout) :: rmin
character(len=1), intent(in) :: which_basis, which_set
real(my_real), intent(out) :: lmat(:,:),lrhs(:)
interface
   subroutine loc_bconds_a(x,y,bmark,itype,c,rs,extra_rs)
   use global
   real(my_real), intent(in) :: x(:),y(:)
   integer, intent(in) :: bmark
   integer, intent(out) :: itype(:)
   real(my_real), intent(out) :: c(:,:,:),rs(:,:)
   real(my_real), intent(out), optional :: extra_rs(:,:,:)
   end subroutine loc_bconds_a
   subroutine loc_bconds_s(x,y,bmark,itype,c,rs)
   use global
   real(my_real), intent(in) :: x,y
   integer, intent(in) :: bmark
   integer, intent(out) :: itype(:)
   real(my_real), intent(out) :: c(:,:),rs(:)
   end subroutine loc_bconds_s
end interface
optional :: loc_bconds_a, loc_bconds_s
real(my_real), intent(out), optional :: lmassmat(:,:), extra_rhs(:,:)


!----------------------------------------------------
! Local variables:

integer :: bctype(ss)
real(my_real) :: cxx(ss,ss), cxy(ss,ss), cyy(ss,ss), cx(ss,ss), cy(ss,ss), &
                 c(ss,ss), rs(ss)
real(my_real), allocatable :: basis(:,:),basisx(:,:),basisy(:,:), &
                              u(:,:,:),ux(:,:,:),uy(:,:,:), &
                              bcrs(:,:),bcc(:,:,:),extra_bcrs(:,:,:)
real(my_real) :: xend(2),yend(2)
integer :: i,j,k,l,qp,side,nquad_pts,jerr,astat,nbasis,qdegree,edim
logical :: doit
logical, save :: warned=.false.
real(my_real), pointer :: xquad(:), yquad(:), quad_weight(:)
!----------------------------------------------------
! Begin executable code

! nodal and h-hierarchical bases require the same degree throughout

if (which_basis == "n" .or. which_basis == "h") then
   if (degree(1) /= degree(4) .or. degree(2) /= degree(4) .or. &
       degree(3) /= degree(4)) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("degree must be the same for edges and face for nodal or h-hierarchical basis", &
                 intlist = degree)
      stop
   endif
endif

! determine the degree of the quadrature rule that is exact for
! polynomials of the same degree as the approximating space.  degree has
! polynomial degree for the edge opposite vertices 1, 2, and 3, and
! the interior of the triangle, in that order.

qdegree = 2*maxval(degree)
if (qdegree > MAX_QUAD_ORDER_TRI) then
   if (.not. warned) then
      call warning("Element degree requires quadrature rule larger than available.", &
                   "Answers may not be as accurate as expected.")
      warned = .true.
   endif
   qdegree = MAX_QUAD_ORDER_TRI
endif
qdegree = qdegree + inc_quad_order
if (qdegree > MAX_QUAD_ORDER_TRI) qdegree = MAX_QUAD_ORDER_TRI
if (qdegree < 1) qdegree = 1

! compute element quadrature points and weights

call quadrature_rule_tri(qdegree,xvert,yvert,nquad_pts,quad_weight,xquad, &
                         yquad,jerr,stay_in=.true.)
if (jerr /= 0) then
   select case(jerr)
   case (1)
      ierr = PHAML_INTERNAL_ERROR
      call fatal("No rule for quadrature of requested order.",intlist=(/qdegree/))
      stop
   case(2)
      ierr = ALLOC_FAILED
      call fatal("allocation failed in quadrature_rules")
      stop
   case default
      ierr = PHAML_INTERNAL_ERROR
      call fatal("Quadrature rule returned error ",intlist=(/jerr/))
      stop
   end select
endif

! number of basis functions over this element

select case (which_set)

case("a")
   nbasis = 3
   do i=1,EDGES_PER_ELEMENT 
      nbasis = nbasis + max(0,degree(i)-1)
   end do
   nbasis = nbasis + max(0,((degree(4)-1)*(degree(4)-2))/2)

case("b")
   if (which_basis == "p") then
      nbasis = 3
      do i=1,EDGES_PER_ELEMENT
         nbasis = nbasis + max(0,degree(i)-2)
      end do
      if (degree(4) > 3) then
         nbasis = nbasis + ((degree(4)-3)*(degree(4)-2))/2
      endif
   else
      if (2*(degree(1)/2) == degree(1)) then
         nbasis = ((degree(1)+2)*(degree(1)+2))/4
      else
         nbasis = ((degree(1)+1)*(degree(1)+3))/4
      endif
   endif

case("r")
   if (which_basis == "p") then
      nbasis = 0
      do i=1,EDGES_PER_ELEMENT
         if (degree(i) >= 2) nbasis = nbasis + 1
      end do
      if (degree(4) >= 3) nbasis = nbasis + degree(4)-2
   else
      if (2*(degree(1)/2) == degree(1)) then
         nbasis = (degree(1)*(degree(1)+2))/4
      else
         nbasis = ((degree(1)+1)*(degree(1)+1))/4
      endif
   endif

end select

! compute the values of the basis functions at the quadrature points

allocate(basis(nbasis,nquad_pts),basisx(nbasis,nquad_pts), &
         basisy(nbasis,nquad_pts), stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in elemental_matrix")
   return
endif

select case (which_basis)
case("n")
   call nodal_basis_func(xquad,yquad,xvert,yvert,degree(1),which_set, &
                         basis,basisx,basisy)
case("h")
   call h_hier_basis_func(xquad,yquad,xvert,yvert,degree(1),which_set, &
                          basis,basisx,basisy)
case("p")
   call p_hier_basis_func(xquad,yquad,xvert,yvert,degree,which_set, &
                          basis,basisx,basisy)
end select

! initialize sums to 0.0

lmat = 0.0_my_real
lrhs = 0.0_my_real
if (present(lmassmat)) lmassmat = 0.0_my_real
if (present(extra_rhs)) extra_rhs = 0.0_my_real

if (present(extra_rhs)) then
   edim = 1+size(extra_rhs,2)
else
   edim = 1
endif

! for residual problems, evaluate the solution and derivatives, and for
! eigenvalue problems, evaluate the solution

if (resid_prob) then
   allocate(u(ss,edim,nquad_pts),ux(ss,edim,nquad_pts),uy(ss,edim,nquad_pts), &
            stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in elemental_matrix")
      return
   endif
   call evaluate_soln_local(grid,xquad,yquad,elem,(/(i,i=1,ss)/), &
                            (/(j,j=1,edim)/),u,ux,uy)

elseif (grid%num_eval > 0) then
   allocate(u(ss,edim,nquad_pts),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in elemental_matrix")
      return
   endif
   call evaluate_soln_local(grid,xquad,yquad,elem,(/(i,i=1,ss)/), &
                            (/(j,j=1,edim)/),u)
endif

! for each quadrature point

do qp = 1,nquad_pts

! evaluate pde coefficients

   call pdecoefs(xquad(qp),yquad(qp),cxx,cxy,cyy,cx,cy,c,rs)

   if (ss == 1) then

! special BLAS based code for scalar PDEs

      if (my_real == kind(1.0)) then
         call sgemm("N","N",nbasis,nbasis,1,quad_weight(qp)*cxx(1,1), &
                    basisx(1,qp),nbasis,basisx(1,qp),1,1.0_my_real,lmat,size(lmat,1))
         call sgemm("N","N",nbasis,nbasis,1,quad_weight(qp)*cyy(1,1), &
                    basisy(1,qp),nbasis,basisy(1,qp),1,1.0_my_real,lmat,size(lmat,1))
         call sgemm("N","N",nbasis,nbasis,1,quad_weight(qp)*cxy(1,1), &
                    basisx(1,qp),nbasis,basisy(1,qp),1,1.0_my_real,lmat,size(lmat,1))
         call sgemm("N","N",nbasis,nbasis,1,quad_weight(qp)*cx(1,1), &
                    basis(1,qp),nbasis,basisx(1,qp),1,1.0_my_real,lmat,size(lmat,1))
         call sgemm("N","N",nbasis,nbasis,1,quad_weight(qp)*cy(1,1), &
                    basis(1,qp),nbasis,basisy(1,qp),1,1.0_my_real,lmat,size(lmat,1))
         call sgemm("N","N",nbasis,nbasis,1,quad_weight(qp)*c(1,1), &
                    basis(1,qp),nbasis,basis(1,qp),1,1.0_my_real,lmat,size(lmat,1))
         if (present(lmassmat)) then
            call sgemm("N","N",nbasis,nbasis,1,quad_weight(qp)*rs(1), &
                       basis(1,qp),nbasis,basis(1,qp),1,1.0_my_real,lmassmat, &
                       size(lmassmat,1))
         endif
      elseif (my_real == kind(1.0d0)) then
         call dgemm("N","N",nbasis,nbasis,1,quad_weight(qp)*cxx(1,1), &
                    basisx(1,qp),nbasis,basisx(1,qp),1,1.0d0,lmat,size(lmat,1))
         call dgemm("N","N",nbasis,nbasis,1,quad_weight(qp)*cyy(1,1), &
                    basisy(1,qp),nbasis,basisy(1,qp),1,1.0d0,lmat,size(lmat,1))
         call dgemm("N","N",nbasis,nbasis,1,quad_weight(qp)*cxy(1,1), &
                    basisx(1,qp),nbasis,basisy(1,qp),1,1.0d0,lmat,size(lmat,1))
         call dgemm("N","N",nbasis,nbasis,1,quad_weight(qp)*cx(1,1), &
                    basis(1,qp),nbasis,basisx(1,qp),1,1.0d0,lmat,size(lmat,1))
         call dgemm("N","N",nbasis,nbasis,1,quad_weight(qp)*cy(1,1), &
                    basis(1,qp),nbasis,basisy(1,qp),1,1.0d0,lmat,size(lmat,1))
         call dgemm("N","N",nbasis,nbasis,1,quad_weight(qp)*c(1,1), &
                    basis(1,qp),nbasis,basis(1,qp),1,1.0d0,lmat,size(lmat,1))
         if (present(lmassmat)) then
            call dgemm("N","N",nbasis,nbasis,1,quad_weight(qp)*rs(1), &
                       basis(1,qp),nbasis,basis(1,qp),1,1.0d0,lmassmat, &
                       size(lmassmat,1))
         endif
      else
         ierr = PHAML_INTERNAL_ERROR
         call fatal("my_real is neither single nor double precision. Can't call GEMM")
         stop 
      endif
      do j=1,nbasis
         if (grid%num_eval <= 0) then
            lrhs(j) = lrhs(j) + quad_weight(qp) * rs(1)*basis(j,qp)
         else
            lrhs(j) = lrhs(j) + quad_weight(qp) * &
                                  rs(1)*grid%eigenvalue(1)*u(1,1,qp)*basis(j,qp)
         endif
         if (present(extra_rhs)) then
            do i=1,size(extra_rhs,2)
               extra_rhs(j,i) = extra_rhs(j,i) + quad_weight(qp) * &
                               rs(1)*grid%eigenvalue(i+1)*u(1,i+1,qp)*basis(j,qp)
            end do
         endif
         if (resid_prob) then
            lrhs(j) = lrhs(j) + quad_weight(qp) * &
                           (-cxx(1,1)*ux(1,1,qp)*basisx(j,qp) &
                            -cyy(1,1)*uy(1,1,qp)*basisy(j,qp) &
                            -cxy(1,1)*uy(1,1,qp)*basisx(j,qp) &
                            - cx(1,1)*ux(1,1,qp)*basis(j,qp) &
                            - cy(1,1)*uy(1,1,qp)*basis(j,qp) &
                              -c(1,1)*u(1,1,qp)*basis(j,qp))
            if (present(extra_rhs)) then
               do i=1,size(extra_rhs,2)
                  extra_rhs(j,i) = extra_rhs(j,i) + quad_weight(qp) * &
                                       (-cxx(1,1)*ux(1,i+1,qp)*basisx(j,qp) &
                                        -cyy(1,1)*uy(1,i+1,qp)*basisy(j,qp) &
                                        -cxy(1,1)*uy(1,i+1,qp)*basisx(j,qp) &
                                        - cx(1,1)*ux(1,i+1,qp)*basis(j,qp) &
                                        - cy(1,1)*uy(1,i+1,qp)*basis(j,qp) &
                                          -c(1,1)*u(1,i+1,qp)*basis(j,qp))
               end do
            endif
         endif
      end do

   else ! ss /= 1

! for each basis pair

      do j=1,nbasis
       do i=1,nbasis

! contribution of this quadrature point to this integrals

         do l=1,ss
            do k=1,ss
               lmat(ss*(i-1)+k,ss*(j-1)+l) = &
                             lmat(ss*(i-1)+k,ss*(j-1)+l) + quad_weight(qp) * &
                                     (cxx(k,l)*basisx(i,qp)*basisx(j,qp) + &
                                      cyy(k,l)*basisy(i,qp)*basisy(j,qp) + &
                                      cxy(k,l)*basisx(i,qp)*basisy(j,qp)+ &
                                       cx(k,l)*basis (i,qp)*basisx(j,qp) + &
                                       cy(k,l)*basis (i,qp)*basisy(j,qp) + &
                                        c(k,l)*basis (i,qp)*basis (j,qp))
               if (present(lmassmat) .and. k==l) then
                  lmassmat(ss*(i-1)+k,ss*(j-1)+l) = &
                                  lmassmat(ss*(i-1)+k,ss*(j-1)+l) + &
                                  quad_weight(qp)*basis(i,qp)*basis(j,qp)*rs(k)
               endif
            end do
         end do
       end do

! right hand side integral

       do k=1,ss
         if (grid%num_eval <= 0) then
             lrhs(ss*(j-1)+k) = lrhs(ss*(j-1)+k) + &
                                   quad_weight(qp) * rs(k)*basis(j,qp)
         else
             lrhs(ss*(j-1)+k) = lrhs(ss*(j-1)+k) + &
               rs(k)* quad_weight(qp) * grid%eigenvalue(1)*u(k,1,qp)*basis(j,qp)
         endif
         if (present(extra_rhs)) then
            do i=1,size(extra_rhs,2)
               extra_rhs(ss*(j-1)+k,i) = extra_rhs(ss*(j-1)+k,i) + &
                               quad_weight(qp) * &
                         rs(k)*grid%eigenvalue(i+1)*u(k,i+1,qp)*basis(j,qp)
            end do
         endif

! if setting up a residual problem, residual part of right hand side

          if (resid_prob) then
             do l=1,ss
                lrhs(ss*(j-1)+k) = lrhs(ss*(j-1)+k) + quad_weight(qp) * &
                                  (-cxx(k,l)*ux(l,1,qp)*basisx(j,qp) &
                                   -cyy(k,l)*uy(l,1,qp)*basisy(j,qp) &
                                   -cxy(k,l)*uy(l,1,qp)*basisx(j,qp) &
                                   - cx(k,l)*ux(l,1,qp)*basis(j,qp) &
                                   - cy(k,l)*uy(l,1,qp)*basis(j,qp) &
                                     -c(k,l)*u(l,1,qp)*basis(j,qp))
             end do
             if (present(extra_rhs)) then
                do i=1,size(extra_rhs,2)
                   do l=1,ss
                      extra_rhs(ss*(j-1)+k,i) = extra_rhs(ss*(j-1)+k,i) + &
                              quad_weight(qp) * &
                              (-cxx(k,l)*ux(l,i+1,qp)*basisx(j,qp) &
                               -cyy(k,l)*uy(l,i+1,qp)*basisy(j,qp) &
                               -cxy(k,l)*uy(l,i+1,qp)*basisx(j,qp) &
                               - cx(k,l)*ux(l,i+1,qp)*basis(j,qp) &
                               - cy(k,l)*uy(l,i+1,qp)*basis(j,qp) &
                                 -c(k,l)*u(l,i+1,qp)*basis(j,qp))
                   end do
                end do
             endif
          endif
       end do
     end do

   endif ! ss == 1

! Determine the minimum value of coefu, to be used as lambda0 when
! solving an eigenvalue problem for the smallest eigenvalues.
! If the equation has been multiplied through by some value (e.g. when
! using cylindrical coordinates it is multiplied by x) that value is
! supposed to be in coefrhs.  Divide coefu by that to get the correct
! lower bound on eigenvalues.

   do i=1,ss
      do j=1,ss
         if (rs(j) /= 0.0_my_real) then
            rmin = min(rmin,c(i,j)/rs(j))
         endif
      end do
   end do

end do

if (resid_prob) then
   deallocate(u,ux,uy)
elseif (grid%num_eval > 0) then
   deallocate(u)
endif
deallocate(quad_weight,xquad,yquad,basis,basisx,basisy)

! compute boundary integrals

! experiments showed some cases needed one more degree of quadrature

qdegree = min(MAX_QUAD_ORDER_LINE,qdegree+1)

do side=1,EDGES_PER_ELEMENT

! do sides that have at least one component that is NATURAL or MIXED

   doit = .false.
   do j=1,ss
     if (edge_type(side,j) == NATURAL .or. edge_type(side,j) == MIXED) then
        doit = .true.
        exit
     endif
   end do

   if (.not. doit) cycle

! compute side quadrature points and weights

   select case (side)
   case (1)
      xend(1) = xvert(2)
      xend(2) = xvert(3)
      yend(1) = yvert(2)
      yend(2) = yvert(3)
   case (2)
      xend(1) = xvert(3)
      xend(2) = xvert(1)
      yend(1) = yvert(3)
      yend(2) = yvert(1)
   case (3)
      xend(1) = xvert(1)
      xend(2) = xvert(2)
      yend(1) = yvert(1)
      yend(2) = yvert(2)
   end select
   call quadrature_rule_line(qdegree,xend,yend,nquad_pts,quad_weight,xquad, &
                             yquad,jerr)
   if (jerr /= 0) then
      select case(jerr)
      case (1)
         ierr = PHAML_INTERNAL_ERROR
         call fatal("No rule for quadrature of requested order.", &
                    intlist=(/qdegree/))
         stop
      case(2)
         ierr = ALLOC_FAILED
         call fatal("allocation failed in quadrature_rules")
         stop
      case default
         ierr = PHAML_INTERNAL_ERROR
         call fatal("Quadrature rule returned error ",intlist=(/jerr/))
         stop
      end select
   endif

! set the values of the basis functions at the quadrature points

   if (present(extra_rhs)) then
      allocate(basis(nbasis,nquad_pts),bcrs(ss,nquad_pts),bcc(ss,ss,nquad_pts),&
               extra_bcrs(ss,size(extra_rhs,dim=2),nquad_pts),stat=astat)
   else
      allocate(basis(nbasis,nquad_pts),bcrs(ss,nquad_pts),bcc(ss,ss,nquad_pts),&
               stat=astat)
   endif
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in elemental_matrix")
      return
   endif

   select case(which_basis)
   case("n")
      call nodal_basis_func(xquad,yquad,xvert,yvert,degree(1),which_set,basis)
   case("h")
      call h_hier_basis_func(xquad,yquad,xvert,yvert,degree(1),which_set,basis)
   case("p")
      call p_hier_basis_func(xquad,yquad,xvert,yvert,degree,which_set,basis)
   end select

! evaluate boundary conditions at all quadrature points

   if (present(loc_bconds_a)) then
      if (present(extra_rhs)) then
         call loc_bconds_a(xquad,yquad,bmark(side),bctype,bcc,bcrs,extra_bcrs)
      else
         call loc_bconds_a(xquad,yquad,bmark(side),bctype,bcc,bcrs)
      endif
   elseif (present(loc_bconds_s)) then
      do qp=1,nquad_pts
         call loc_bconds_s(xquad(qp),yquad(qp),bmark(side),bctype,bcc(:,:,qp), &
                           bcrs(:,qp))
         if (present(extra_rhs)) then
            do i=1,size(extra_bcrs,dim=2)
               extra_bcrs(:,i,qp) = bcrs(:,qp)
            end do
         endif
      end do
   else
      ierr = PHAML_INTERNAL_ERROR
      call fatal("must give one of loc_bconds forms in elemental_matrix")
      stop
   endif

! for each quadrature point

   do qp=1,nquad_pts

! for each pair of basis functions

      do j=1,nbasis
         do i=1,nbasis
            do l=1,ss
               do k=1,ss
                  if (bctype(k) /= MIXED) cycle
                  lmat(ss*(i-1)+k,ss*(j-1)+l) = &
                     lmat(ss*(i-1)+k,ss*(j-1)+l) + &
                     quad_weight(qp)*bcc(k,l,qp)*basis(i,qp)*basis(j,qp)
               end do
            end do
         end do
         do k=1,ss
            if (bctype(k) /= NATURAL .and. bctype(k) /= MIXED) cycle
            lrhs(ss*(j-1)+k) = lrhs(ss*(j-1)+k) + &
                                    quad_weight(qp)*bcrs(k,qp)*basis(j,qp)
            if (present(extra_rhs)) then
               do i=1,size(extra_rhs,2)
                  extra_rhs(ss*(j-1)+k,i) = extra_rhs(ss*(j-1)+k,i) + &
                                  quad_weight(qp)*extra_bcrs(k,i,qp)*basis(j,qp)
               end do
            endif
         end do
      end do
   end do

   deallocate(quad_weight,xquad,yquad,basis,bcrs,bcc)
   if (present(extra_rhs)) then
      deallocate(extra_bcrs)
   endif

end do

end subroutine elemental_matrix

!          --------
subroutine assemble(linear_system,eq_num,nleq,local_stiffness,local_rhs, &
                    local_mass)
!          --------

!----------------------------------------------------
! This routine adds the local matrix from one element to the global matrix.
! The equations for that element are given in eq_num.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(inout) :: linear_system
integer, intent(in) :: eq_num(:), nleq
real(my_real), intent(in) :: local_stiffness(:,:),local_rhs(:),local_mass(:,:)

!----------------------------------------------------
! Local variables:

integer :: eq, col, i, leq, lcol, syssize, eq1, eq2, s1, s2

!----------------------------------------------------
! Begin executable code

syssize = linear_system%system_size

! for each equation ...

do eq1=1,nleq
   do s1=1,syssize
      leq = syssize*(eq1-1)+s1
      eq = eq_num(eq1)+s1-1

! add contribution to right hand side

      linear_system%rhs(eq) = linear_system%rhs(eq) + local_rhs(leq)

! for each column

      do eq2=1,nleq
         do s2=1,syssize
            lcol = syssize*(eq2-1)+s2
            col = eq_num(eq2)+s2-1

! find the column in the global matrix

            do i=linear_system%begin_row(eq),linear_system%end_row(eq)
               if (linear_system%column_index(i)==col) exit
            end do
            if (linear_system%column_index(i) /= col) then
               ierr = PHAML_INTERNAL_ERROR
               call fatal("failed to find column in assemble")
               return
            endif

! add contribution to the matrix

            linear_system%stiffness(i) = linear_system%stiffness(i) + &
                                         local_stiffness(leq,lcol)
            if (associated(linear_system%mass)) then
               linear_system%mass(i)=linear_system%mass(i)+local_mass(leq,lcol)
            endif

         end do ! next PDE
      end do ! next column
   end do ! next PDE
end do ! next equation

end subroutine assemble

!          --------------
subroutine fix_quad_error(grid,procs,linear_system,still_sequential,timeit)
!          --------------

!----------------------------------------------------
! This routine gets the contributions to the matrix and right hand side
! from elements owned by other processors.  This is used instead of
! computing them myself so that the quadrature errors are on the order of
! the fine grid elements instead of the element I have.
! This only works for linear elements
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
type(proc_info), intent(in) :: procs
type(linsys_type), intent(inout) :: linear_system
logical, intent(in) :: still_sequential, timeit

!----------------------------------------------------
! Local variables:

integer :: i, lev, eq, indx, allocstat, my_processor, nonzero, &
           nproc,p, j, k, col, elem, nisend, isub, rsub, sisub, srsub
logical :: rhs_next, mass_next
logical(small_logical) :: need_fix(linear_system%neq_vert)
type(hash_key_eq) :: gid
integer, allocatable :: send_int(:), nirecv(:), nrrecv(:),nisendv(:), nrsendv(:)
integer, pointer :: recv_int(:)
real(my_real), allocatable :: send_real(:)
real(my_real), pointer :: recv_real(:)

!----------------------------------------------------
! Begin executable code

nproc = num_proc(procs)
my_processor = my_proc(procs)

! master doesn't do this and nothing to do if sequential

if (my_processor == MASTER .or.  still_sequential .or. &
    nproc == 1) return

! some memory for messages

allocate(nirecv(nproc),nrrecv(nproc),nisendv(nproc),nrsendv(nproc), &
         stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in fix_quad_error",procs=procs)
   return
endif

! determine which equations need fixing

need_fix = .false.
elem=grid%head_level_elem(1)
do while (elem /= END_OF_LIST)
   call mark_need_fix(grid,linear_system,need_fix,elem)
   elem = grid%element(elem)%next
end do

! convert to the hierarchical basis

linear_system%matrix_val => linear_system%stiffness
do lev=linear_system%nlev+2,2,-1
   call basis_change(lev,TO_HIER,linear_system,do_mass=.true.)
end do

! count the number of equations that need fixing

nisend = 0
do eq=1,linear_system%neq_vert
   if (need_fix(eq)) then
      nisend = nisend + KEY_SIZE+1
   endif
end do

allocate(send_int(nisend),stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in fix_quad_error",procs=procs)
   return
endif

! make a list of the equations that need fixing

nisend = 0
do eq=1,linear_system%neq_vert
   if (need_fix(eq)) then
      call hash_pack_key(linear_system%gid(eq),send_int,nisend+1)
      nisend = nisend + KEY_SIZE+1
   endif
end do

! send the list

if (timeit) call start_watch((/cpassemble,ctassemble/))
call phaml_alltoall(procs,send_int,nisend,recv_int,nirecv,1102)
if (timeit) call stop_watch((/cpassemble,ctassemble/))

deallocate(send_int,stat=allocstat)

! count the number of entries to be sent to each processor

nisendv = 0
nrsendv = 0
isub = 1
do p=1,nproc
   if (p == my_processor) then
      isub = isub + nirecv(p)
      cycle
   endif
   do i=1,nirecv(p)/(KEY_SIZE+1)
      gid = hash_unpack_key(recv_int,isub,.true.)
      isub = isub + KEY_SIZE+1
      eq = hash_decode_key(gid,linear_system%eq_hash)
      if (eq /= HASH_NOT_FOUND) then

! how many nonzeros in the matrix row

         nonzero = 0
         do j=linear_system%begin_row(eq),linear_system%end_row_linear(eq)
            if (linear_system%column_index(j) == NO_ENTRY) cycle
            if (associated(linear_system%mass)) then
               if (linear_system%stiffness(j) /= 0.0_my_real .or. &
                   linear_system%mass(j) /= 0.0_my_real) nonzero = nonzero + 2
            else
               if (linear_system%stiffness(j) /= 0.0_my_real) then
                  nonzero = nonzero + 1
               endif
            endif
         end do

! count the rhs if there are any nonzeros in the matrix or the rhs is nonzero,
! plus also one for a NULL_KEY marker to end the row

         if (nonzero > 0 .or. linear_system%rhs(eq) /= 0.0_my_real) then
            nisendv(p) = nisendv(p) + (nonzero+2)*(KEY_SIZE+1)
            nrsendv(p) = nrsendv(p) + nonzero+1
         else
            nisendv(p) = nisendv(p) + nonzero*(KEY_SIZE+1)
            nrsendv(p) = nrsendv(p) + nonzero
         endif

      endif
   end do
end do

! allocate memory for the messages

allocate(send_int(sum(nisendv)),send_real(sum(nrsendv)),stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in fix_quad_error",procs=procs)
   return
endif

! pack the data
! the data for each row consists of
!   first the equation gid in the int and rhs in the real
!   then the stiffness matrix and mass matrix values in the real and column
!   gid in the int whereever either the stiffness matrix or mass matrix
!   is nonzero, finally a NULL hash key to indicate end-of-row.

isub = 1
sisub = 1
srsub = 1
do p=1,nproc
   if (p == my_processor) then
      isub = isub + nirecv(p)
      cycle
   endif
   do i=1,nirecv(p)/(KEY_SIZE+1)
      gid = hash_unpack_key(recv_int,isub,.true.)
      isub = isub + KEY_SIZE+1
      eq = hash_decode_key(gid,linear_system%eq_hash)
      if (eq /= HASH_NOT_FOUND) then

! how many nonzeros in the matrix row

         nonzero = 0
         do j=linear_system%begin_row(eq),linear_system%end_row_linear(eq)
            if (linear_system%column_index(j) == NO_ENTRY) cycle
            if (associated(linear_system%mass)) then
               if (linear_system%stiffness(j) /= 0.0_my_real .or. &
                   linear_system%mass(j) /= 0.0_my_real) then
                  nonzero = 2
                  exit
               endif
            else
               if (linear_system%stiffness(j) /= 0.0_my_real) then
                  nonzero = 1
                  exit
               endif
            endif
         end do

! pack the data

         if (nonzero > 0 .or. linear_system%rhs(eq) /= 0.0_my_real) then
            send_int(sisub:sisub+KEY_SIZE) = recv_int(isub-(KEY_SIZE+1):isub-1)
            sisub = sisub + KEY_SIZE+1
            send_real(srsub) = linear_system%rhs(eq)
            srsub = srsub + 1
            do j=linear_system%begin_row(eq),linear_system%end_row_linear(eq)
               if (linear_system%column_index(j) == NO_ENTRY) cycle
               if (associated(linear_system%mass)) then
                  if (linear_system%stiffness(j) /= 0.0_my_real .or. &
                      linear_system%mass(j) /= 0.0_my_real) then
                     call hash_pack_key(linear_system%gid(linear_system%column_index(j)), &
                                        send_int,sisub)
                     sisub = sisub + KEY_SIZE+1
                     send_real(srsub) = linear_system%stiffness(j)
                     srsub = srsub + 1
! TEMP I really only need the column index once for mass and stiffness
                     call hash_pack_key(linear_system%gid(linear_system%column_index(j)), &
                                        send_int,sisub)
                     sisub = sisub + KEY_SIZE+1
                     send_real(srsub) = linear_system%mass(j)
                     srsub = srsub + 1
                  endif
               else
                  if (linear_system%stiffness(j) /= 0.0_my_real) then
                     call hash_pack_key(linear_system%gid(linear_system%column_index(j)), &
                                        send_int,sisub)
                     sisub = sisub + KEY_SIZE+1
                     send_real(srsub) = linear_system%stiffness(j)
                     srsub = srsub + 1
                  endif
               endif
            end do
            call hash_pack_key(NULL_KEY_EQ,send_int,sisub)
            sisub = sisub + KEY_SIZE+1
         endif
      endif
   end do
end do

if (associated(recv_int)) deallocate(recv_int,stat=allocstat)

! send the lists to the other processors

if (timeit) call start_watch((/cpassemble,ctassemble/))
call phaml_alltoall(procs,send_int,nisendv,recv_int,nirecv,1103)
call phaml_alltoall(procs,send_real,nrsendv,recv_real,nrrecv,1104)
if (timeit) call stop_watch((/cpassemble,ctassemble/))

deallocate(send_int,send_real,stat=allocstat)

! receive the lists and add contributions to my linear system

isub = 1
rsub = 1
do p=1,nproc
   if (p == my_processor) cycle
   rhs_next = .true.
   mass_next = .false.
   do j=1,nirecv(p)/(KEY_SIZE+1)
      gid = hash_unpack_key(recv_int,isub,.true.)
      isub = isub + KEY_SIZE+1
      if (gid == NULL_KEY_EQ) then
         rhs_next = .true.
      elseif (rhs_next) then
         rhs_next = .false.
         eq = hash_decode_key(gid,linear_system%eq_hash)
         if (eq == HASH_NOT_FOUND) then
            call warning("received data for equation I don't have in fix_quad_error")
            call hash_print_key(gid,errunit)
            cycle
         endif
         linear_system%rhs(eq) = linear_system%rhs(eq) + recv_real(rsub)
         rsub = rsub + 1
      elseif (mass_next) then
! indx was computed in the else clause of the previous time through the loop
         if (indx /= 0) then
            linear_system%mass(indx) = linear_system%mass(indx) +recv_real(rsub)
            rsub = rsub + 1
         endif
         mass_next = .false.
      else
         if (associated(linear_system%mass)) mass_next = .true.
         if (eq == HASH_NOT_FOUND) cycle
         col = hash_decode_key(gid,linear_system%eq_hash)
         if (col == HASH_NOT_FOUND) then
            call warning("received data for column I don't have in fix_quad_error")
            call hash_print_key(gid,errunit)
            cycle
         endif
! find the location of this column, if I have it
         indx = 0
         do k=linear_system%begin_row(eq),linear_system%end_row_linear(eq)
            if (linear_system%column_index(k) == col) then
               indx = k
               exit
            endif
         end do
         if (indx /= 0) then
            linear_system%stiffness(indx) = linear_system%stiffness(indx) + &
                                            recv_real(rsub)
            rsub = rsub + 1
         else
            call warning("couldn't find column in fix_quad_error")
         endif
      endif
   end do
end do

! convert back to nodal basis

do lev=2,linear_system%nlev+2
   call basis_change(lev,TO_NODAL,linear_system,do_mass=.true.)
end do

! free message memory

if (associated(recv_int)) deallocate(recv_int,stat=allocstat)
if (associated(recv_real)) deallocate(recv_real,stat=allocstat)
deallocate(nisendv,nrsendv,nirecv,nrrecv,stat=allocstat)

end subroutine fix_quad_error

!                    -------------
recursive subroutine mark_need_fix(grid,linear_system,need_fix,elem)
!                    -------------

!----------------------------------------------------
! This routine determines which equations need quadrature contributions
! from other processors
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
type(linsys_type), intent(in) :: linear_system
logical(small_logical), intent(inout) :: need_fix(:)
integer, intent(in) :: elem
!----------------------------------------------------
! Local variables:

integer :: i,j,eq,child(MAX_CHILD),allc(MAX_CHILD)
!----------------------------------------------------
! Begin executable code

! recursively traverse the refinement tree.  When I own an element, quit
! traversing that branch.  Otherwise, mark the equations associated with
! the three vertices as needing to be fixed.

if (grid%element(elem)%iown) return

do i=1,VERTICES_PER_ELEMENT
   do j=1,linear_system%system_size
      call grid_to_eq(grid,linear_system,VERTEX_ID,1,j, &
                      grid%vertex(grid%element(elem)%vertex(i))%gid,eq)
      need_fix(eq) = .true.
   end do
end do

allc = ALL_CHILDREN
child = get_child_lid(grid%element(elem)%gid,allc,grid%elem_hash)

if (child(1) == NO_CHILD) return

do i=1,MAX_CHILD
   call mark_need_fix(grid,linear_system,need_fix,child(i))
end do

end subroutine mark_need_fix

!          -------------------
subroutine static_condensation(linear_system,grid,procs)
!          -------------------

!----------------------------------------------------
! This routine performs static condensation to eliminate the face bases.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(inout) :: linear_system
type(grid_type), intent(in) :: grid
type(proc_info), intent(in) :: procs
!----------------------------------------------------
! Local variables:

integer :: i, j, k, row, col, eq, objtype, brank, srank, lid, degree, loc_neq, &
           syssize, rowlen, allocstat, info
real(my_real) :: temp
real(my_real), allocatable :: row_block(:,:), row_block_t(:,:)
real(my_real), pointer :: block(:,:)
integer, pointer :: ipiv(:)
!----------------------------------------------------
! Begin executable code

syssize = linear_system%system_size

! factor the blocks of edge bases from one edge and face bases from one element

nullify(linear_system%edge_block,linear_system%elem_block)

! allocate the structures that contain the factored blocks

allocate(linear_system%edge_block(size(grid%edge)), &
         linear_system%elem_block(size(grid%element)),stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation for element block failed",procs=procs)
   return
endif
do i=1,size(linear_system%edge_block)
   nullify(linear_system%edge_block(i)%matrix, &
           linear_system%edge_block(i)%ipiv)
   linear_system%edge_block(i)%neq = 0
end do
do i=1,size(linear_system%elem_block)
   nullify(linear_system%elem_block(i)%matrix, &
           linear_system%elem_block(i)%ipiv)
   linear_system%elem_block(i)%neq = 0
end do
      
! go through the high order equations, identifying each block corresponding to
! an edge or element

eq = linear_system%begin_level(linear_system%nlev+1)
do
   if (eq >= linear_system%begin_level(linear_system%nlev+3)) exit

! find the edge or element corresponding to the block this row starts

   call eq_to_grid(linear_system,linear_system%gid(eq),objtype,brank,srank, &
                   lid,grid)

! make sure this is the beginning of a block

   if (brank /= 1 .or. srank /= 1) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("not at beginning of block in factorization of high order equations", &
                 intlist=(/eq,brank,srank,objtype,i/),procs=procs)
      stop
   endif

! get the edge/element degree, determine number of equations, and allocate mem

   select case(objtype)

   case (EDGE_ID)
      degree = grid%edge(lid)%degree
      loc_neq = (degree-1)*syssize
         
      allocate(linear_system%edge_block(lid)%matrix(loc_neq,loc_neq), &
               linear_system%edge_block(lid)%ipiv(loc_neq),stat=allocstat)
      if (allocstat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("memory allocation for edge block failed",procs=procs)
         return
      endif

      linear_system%edge_block(lid)%neq = loc_neq
      block => linear_system%edge_block(lid)%matrix
      ipiv => linear_system%edge_block(lid)%ipiv

   case (ELEMENT_ID)
      degree = grid%element(lid)%degree
      loc_neq = syssize*((degree-1)*(degree-2))/2
         
      allocate(linear_system%elem_block(lid)%matrix(loc_neq,loc_neq), &
               linear_system%elem_block(lid)%ipiv(loc_neq),stat=allocstat)
      if (allocstat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("memory allocation for element block failed",procs=procs)
         return
      endif
      linear_system%elem_block(lid)%neq = loc_neq
      block => linear_system%elem_block(lid)%matrix
      ipiv => linear_system%elem_block(lid)%ipiv

   case (VERTEX_ID)
      ierr = PHAML_INTERNAL_ERROR
      call fatal("got a vertex ID when examining high order equation", &
                  procs=procs)
      stop

   case default
      ierr = PHAML_INTERNAL_ERROR
      call fatal("object type invalid  when examining high order equation", &
                 intlist=(/objtype/),procs=procs)
      stop

   end select

! copy the entries for this edge or face into block, and the rows and rhs
! into row_block and row_block_t (transpose)

   rowlen = linear_system%end_row(eq)-linear_system%begin_row(eq)+1
   allocate(row_block(loc_neq,rowlen+1),row_block_t(rowlen+1,loc_neq), &
            stat=allocstat)
   if (allocstat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation in make_linear_system failed",procs=procs)
      return
   endif
   row_block = 0.0_my_real
   row_block_t = 0.0_my_real

   block = 0.0_my_real
   do i=1,loc_neq
      do j=linear_system%begin_row(eq+i-1),linear_system%end_row(eq+i-1)
         col = linear_system%column_index(j)
         if (col >= eq .and. col < eq+loc_neq) then
            block(i,col-eq+1) = linear_system%condensed(j)
         endif
         row_block(i,j-linear_system%begin_row(eq+i-1)+1) = &
            linear_system%condensed(j)
         row_block_t(j-linear_system%begin_row(eq+i-1)+1,i) = &
             get_matval(linear_system,linear_system%condensed,col,eq+i-1)
      end do
      row_block(i,rowlen+1) = linear_system%rhs(eq+i-1)
      row_block_t(rowlen+1,i) = linear_system%rhs(eq+i-1)
   end do

! factor the block to PLU

   if (my_real == kind(0.0)) then
      call sgetrf(loc_neq,loc_neq,block,loc_neq,ipiv,info)
   elseif (my_real == kind(0.0d0)) then
      call dgetrf(loc_neq,loc_neq,block,loc_neq,ipiv,info)
   else
      ierr = PHAML_INTERNAL_ERROR
      call fatal("lapack requires my_real is either default single or double precision")
      return
   endif
   if (info /= 0) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("lapack sgetrf failed for factoring high order bases", &
                 intlist=(/info/),procs=procs)
      stop
   endif

! form the Schur complement A1 - A21 A2^-1 A12 for the condensed equations
! (and also the rhs) where A1 is vertex and edge bases and A2 is face bases.

   if (objtype == ELEMENT_ID) then

! multiply the rows and rhs by A2^-1

      if (my_real == kind(0.0)) then
         call sgetrs("N",loc_neq,rowlen+1,block,loc_neq,ipiv,row_block, &
                     loc_neq,info)
      else
         call dgetrs("N",loc_neq,rowlen+1,block,loc_neq,ipiv,row_block, &
                     loc_neq,info)
      endif
      if (info /= 0) then
         ierr = PHAML_INTERNAL_ERROR
         call fatal("lapack sgetrs solution failed for static condensation",&
                    intlist=(/info/),procs=procs)
         stop
      endif

! form the Schur complement

      do i=1,rowlen
         row = linear_system%column_index(linear_system%begin_row(eq)+i-1)
         if (row >= eq .and. row < eq+loc_neq) cycle
         do j=1,rowlen
            col = linear_system%column_index(linear_system%begin_row(eq)+j-1)
            if (col >= eq .and. col < eq+loc_neq) cycle
! TEMP this should be done with level 3 blas
            temp = 0.0_my_real
            do k=1,loc_neq
               temp = temp + row_block_t(i,k)*row_block(k,j)
            end do
            do k=linear_system%begin_row(row),linear_system%end_row(row)
               if (linear_system%column_index(k) == col) exit
            end do
            if (k > linear_system%end_row(row)) then
               ierr = PHAML_INTERNAL_ERROR
               call fatal("couldn't find col in row during static condensation", &
                          intlist=(/col,row/),procs=procs)
               stop
            endif
            linear_system%condensed(k) = linear_system%condensed(k) - temp
         end do
         temp = 0.0_my_real
         do k=1,loc_neq
            temp = temp + row_block_t(i,k)*row_block(k,rowlen+1)
         end do
         linear_system%rhs_cond(row) = linear_system%rhs_cond(row) - temp
      end do
   endif ! face basis

   deallocate(row_block,row_block_t,stat=allocstat)
   eq = eq + loc_neq
end do ! blocks of high order equations

! after static condensation, the matrix rows should not contain the face bases,
! the number of equations should not include face bases, and the matrix values
! and rhs should be the Schur complement

linear_system%end_row => linear_system%end_row_edge
linear_system%neq = linear_system%neq_vert + linear_system%neq_edge
linear_system%matrix_val => linear_system%condensed
linear_system%rhs => linear_system%rhs_cond

end subroutine static_condensation

!        ----------
function get_matval(linear_system,matval,row,col)
!        ----------

!----------------------------------------------------
! This routine get the (row,col) matrix value in matval
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(in) :: linear_system
real(my_real), intent(in) :: matval(:)
integer, intent(in) :: row,col
real(my_real) :: get_matval
!----------------------------------------------------
! Local variables:

integer :: i

!----------------------------------------------------
! Begin executable code

! TEMP this is used to build the transpose block in static condensation.
!      If that's all it is used for, this can be made more efficient when
!      the matrix is symmetric, and maybe other times, too

! seach row for column col to get the value; return 0 if col is not found

get_matval = 0.0_my_real
do i=linear_system%begin_row(row), linear_system%end_row(row)
   if (linear_system%column_index(i) == col) then
      get_matval = matval(i)
      exit
   endif
end do

end function get_matval

!          ---------------------
subroutine destroy_linear_system(linear_system)
!          ---------------------

!----------------------------------------------------
! This routine frees the memory used by the linear system
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type) :: linear_system

!----------------------------------------------------
! local variables

integer :: i, allocstat

!----------------------------------------------------
! Begin executable code

if (associated(linear_system%begin_row)) then
   if (associated(linear_system%condensed,linear_system%stiffness)) then
      nullify(linear_system%condensed)
   endif
   if (associated(linear_system%condensed,linear_system%shifted)) then
      nullify(linear_system%condensed)
   endif
   if (associated(linear_system%rhs_cond,linear_system%rhs_nocond)) then
      nullify(linear_system%rhs_cond)
   endif
   deallocate(linear_system%stiffness, &
              linear_system%column_index, &
              linear_system%begin_row, &
              linear_system%end_row_face, &
              linear_system%end_row_linear, &
              linear_system%end_row_edge, &
              linear_system%begin_level, &
              linear_system%rhs_nocond, &
              linear_system%r_mine, &
              linear_system%r_others, &
              linear_system%need_r_others, &
              linear_system%iown, &
              linear_system%solution, &
              linear_system%equation_type, &
              linear_system%gid, &
              linear_system%hold_dirich, &
              linear_system%s_int, &
              linear_system%s_bnd, &
              linear_system%s_int_inv, &
              linear_system%s_bnd_inv, &
              stat=allocstat)
   if (allocstat /= 0) then
      call warning("deallocation failed in destroy_linear_system")
   endif
   nullify(linear_system%end_row)
   if (associated(linear_system%evecs)) then
      deallocate(linear_system%evecs,stat=allocstat)
      if (allocstat /= 0) then
         call warning("deallocation failed in destroy_linear_system")
      endif
   endif
   nullify(linear_system%matrix_val)
   if (associated(linear_system%mass)) then
      deallocate(linear_system%mass, linear_system%shifted, stat=allocstat)
      if (allocstat /= 0) then
         call warning("deallocation failed in destroy_linear_system (mass)")
      endif
   endif
   if (associated(linear_system%condensed)) then
      deallocate(linear_system%condensed,linear_system%rhs_cond,stat=allocstat)
      if (allocstat /= 0) then
         call warning("deallocation failed in destroy_linear_system (condensed)")
      endif
   endif
   call hash_table_destroy(linear_system%eq_hash)
   if (associated(linear_system%edge_block)) then
      do i=1,size(linear_system%edge_block)
         if (associated(linear_system%edge_block(i)%matrix)) then
            deallocate(linear_system%edge_block(i)%matrix, &
                       linear_system%edge_block(i)%ipiv,stat=allocstat)
            if (allocstat /= 0) then
               call warning("deallocation failed in destroy_linear_system (edge_block component)")
            endif
         endif
      end do
      deallocate(linear_system%edge_block,stat=allocstat)
      if (allocstat /= 0) then
         call warning("deallocation failed in destroy_linear_system (edge_block)")
      endif
   endif
   if (associated(linear_system%elem_block)) then
      do i=1,size(linear_system%elem_block)
         if (associated(linear_system%elem_block(i)%matrix)) then
            deallocate(linear_system%elem_block(i)%matrix, &
                       linear_system%elem_block(i)%ipiv,stat=allocstat)
            if (allocstat /= 0) then
               call warning("deallocation failed in destroy_linear_system (elem_block component)")
            endif
         endif
      end do
      deallocate(linear_system%elem_block,stat=allocstat)
      if (allocstat /= 0) then
         call warning("deallocation failed in destroy_linear_system (elem_block)")
      endif
   endif
   if (associated(linear_system%nn_comm_remote_neigh)) then
      deallocate(linear_system%nn_comm_remote_neigh)
   endif
   if (associated(linear_system%nn_comm_end_of_send)) then
      do i=1,size(linear_system%nn_comm_send_lid)
         deallocate(linear_system%nn_comm_send_lid(i)%lid, &
                    linear_system%nn_comm_recv_lid(i)%lid)
      end do
      deallocate(linear_system%nn_comm_end_of_send, &
                 linear_system%nn_comm_send_lid, &
                 linear_system%nn_comm_recv_lid)
   endif
   if (associated(linear_system%fudop_comm_proxy)) then
      deallocate(linear_system%fudop_comm_proxy)
      do i=1,size(linear_system%fudop_comm_send_lid)
         deallocate(linear_system%fudop_comm_send_lid(i)%lid, &
                    linear_system%fudop_comm_recv_lid(i)%lid, &
                    linear_system%resid_comm_send_lid(i)%lid, &
                    linear_system%resid_comm_recv_lid(i)%lid)
      end do
      deallocate(linear_system%fudop_comm_nsend, &
                 linear_system%fudop_comm_send_lid, &
                 linear_system%fudop_comm_recv_lid, &
                 linear_system%fudop_comm_mess_size, &
                 linear_system%fudop_comm_nchunk, &
                 linear_system%fudop_comm_nproc_recv, &
                 linear_system%resid_comm_nsend, &
                 linear_system%resid_comm_send_lid, &
                 linear_system%resid_comm_recv_lid, &
                 linear_system%resid_comm_mess_size, &
                 linear_system%resid_comm_nchunk, &
                 linear_system%resid_comm_nproc_recv)
   endif
endif

end subroutine destroy_linear_system

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

! TEMP080114 charged particles
! Undocumented feature for computing the electric field of charged particles
!          ---------------------
subroutine charged_particles_rhs(linsys,grid)
!          ---------------------

!----------------------------------------------------
! This routine changes the right hand side of the linear system for an
! equation that computes the electric field due to charged particles.
! It is assumed that the particle charge density is a weighted sum of
! delta functions which are centered at interior vertices of the initial
! grid.  The weights are the charges.  The right hand side of the PDE
! (rs in subroutine pdecoefs) must return the charge if (x,y) is the
! location of a particle, and 0.0 otherwise.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(inout) :: linsys
type(grid_type), intent(in) :: grid
!----------------------------------------------------
! Local variables:

! TEMP assume the equation is scalar, i.e. system_size == 1
real(my_real) :: cxx(1,1),cxy(1,1),cyy(1,1),cx(1,1),cy(1,1),c(1,1),rs(1),x,y
integer :: eq, object_type, basis_rank, system_rank, grid_lid
!----------------------------------------------------
! Begin executable code

! TEMP assume the equation is scalar, i.e. system_size == 1

if (linsys%system_size /= 1) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("systems of equations for charged particles is not implemented")
   stop
endif

! initialize the right hand side to be 0, except at Dirichlet boundary
! points where it keeps the boundary value previously assigned.

where (linsys%equation_type /= DIRICHLET) rs = 0.0_my_real

! go through the vertices of the initial grid (i.e. equations on level 1)
! skipping the Dirichlet points

do eq=linsys%begin_level(1),linsys%begin_level(2)-1
   if (linsys%equation_type(eq) == DIRICHLET) cycle

! get the grid point associated with this equation, and its coordinates

   call eq_to_grid(linsys,linsys%gid(eq),object_type,basis_rank,system_rank, &
                   grid_lid,grid)
   x = grid%vertex(grid_lid)%coord%x
   y = grid%vertex(grid_lid)%coord%y

! evaluate the PDE right hand side at this point and assign it to the
! linear system right hand side

   call pdecoefs(x,y,cxx,cxy,cyy,cx,cy,c,rs)
   linsys%rhs(eq) = rs(1)

end do

end subroutine charged_particles_rhs

end module make_linsys
