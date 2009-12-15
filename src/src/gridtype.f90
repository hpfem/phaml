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

module gridtype_mod

!----------------------------------------------------
! This module contains data structures for the grid.
!----------------------------------------------------

!----------------------------------------------------
! Other modules used:

use global
use hash_mod
use message_passing
!----------------------------------------------------

implicit none
private
public VERTICES_PER_ELEMENT, EDGES_PER_ELEMENT, MAX_CHILD, ALL_CHILDREN, &
       END_OF_LIST, NOT_ON_LIST, NO_CHILD, BOUNDARY, point, element_t, &
       edge_t, vertex_t, grid_type, refine_options, errind_list, &
       triangle_data, get_child_lid, get_child_gid, get_neighbors, binw, &
       get_grid_info

!----------------------------------------------------
! The following parameters are defined:

! values that depend on the kind of elements, polynomial degree, etc.
! RESTRICTION This version for bisected triangles.

integer, parameter :: VERTICES_PER_ELEMENT = 3, &
                      EDGES_PER_ELEMENT    = 3, &
                      MAX_CHILD            = 2

! argument for get_child when all children are wanted

integer, private :: i
integer, parameter :: ALL_CHILDREN(MAX_CHILD) = (/ (i,i=1,MAX_CHILD) /)

! flags

integer, parameter :: END_OF_LIST = -10, & ! end of linked list
                      NOT_ON_LIST = -11, & ! link value when not on any list
                      NO_CHILD    = -12    ! child does not exist

! neighbor or mate of an element when that side is on the boundary

integer, parameter :: BOUNDARY = -6

!----------------------------------------------------
! The following types are defined:

type point
   real(my_real) :: x,y ! RESTRICTION 2D, add z for 3D
end type point

type element_t
   type(hash_key) :: gid
   type(hash_key) :: mate ! RESTRICTION 2D, multiple mates in 3D
   real(my_real) :: weight
! subscripts are basis rank, system rank, eigenvalue
   real(my_real), pointer :: solution(:,:,:), exact(:,:,:), oldsoln(:,:,:)
   real(my_real) :: work
   real(my_real) :: sp_eta_pred
   integer :: vertex(VERTICES_PER_ELEMENT)
   integer :: edge(EDGES_PER_ELEMENT)
   integer :: degree
   integer :: level
   integer :: in, out
   integer :: order(MAX_CHILD)
   integer :: next, previous ! links for available memory list when not in use
                             ! and elements of one level when in use
   logical(small_logical) :: isleaf, oldleaf
   logical(small_logical) :: iown ! true if I own all leaves below this element
   logical(small_logical) :: hrefined_unowned, prefined_unowned
end type element_t

! edge owner is the same as it's vertex 2
! there is no linked list of edges in use.  next is used for linked list of
! free memory and pointer from PERIODIC_SLAVE to PERIODIC_MASTER

type edge_t
   type(hash_key) :: gid
   integer :: vertex(2)
   integer :: bmark, degree, assoc_elem
! subscripts are basis rank, system rank, eigenvalue
   real(my_real), pointer :: solution(:,:,:), exact(:,:,:), oldsoln(:,:,:)
   integer :: next
end type edge_t

type vertex_t
   type(hash_key) :: gid
   type(point) :: coord
   real(my_real) :: bparam
   integer :: bmark
   integer :: assoc_elem
   integer :: next, previous
end type vertex_t

! if not solving an eigenproblem, num_eval is 0 and nsoln is system_size.

type grid_type
   type(element_t), pointer :: element(:)
   type(edge_t), pointer :: edge(:)
   type(vertex_t), pointer :: vertex(:)
   type(hash_table) :: elem_hash, edge_hash, vert_hash
   type(point) :: boundbox_min, boundbox_max
! subscripts are lid, system rank, eigenvalue
   real(my_real), pointer :: vertex_solution(:,:,:), vertex_exact(:,:,:), &
                             vertex_oldsoln(:,:,:)
! subscripts are element, eigenvalue
   real(my_real), pointer :: element_errind(:,:)
   real(my_real), pointer :: eigenvalue(:)
   real(my_real) :: eigen_linsys_max_l2_resid, eigen_linsys_ave_l2_resid
   real(my_real), pointer :: eigenprob_l2_resid(:), eigenprob_variance(:)
   real(my_real), pointer :: errest_energy(:), errest_Linf(:), errest_L2(:), &
                             errest_eigenvalue(:)
   real(my_real) :: refsoln_errest
   real(my_real) :: max_blen
   real(my_real), pointer :: bp_start(:), bp_finish(:)
   integer, pointer :: edge_type(:,:), vertex_type(:,:)
   integer, pointer :: initial_neighbor(:,:) ! (EDGES_PER_ELEMENT,nelem_init)
   integer, pointer :: head_level_elem(:), head_level_vert(:)
   integer :: next_free_elem, next_free_edge, next_free_vert
   integer :: partition
   integer :: system_size, num_eval, nsoln
   integer :: nelem, nelem_leaf, nelem_leaf_own, nedge, nedge_own, &
              nvert, nvert_own, nlev, dof, dof_own
   integer :: arpack_iter, arpack_nconv, arpack_numop, arpack_numopb, &
              arpack_numreo, arpack_info
   integer :: errtype ! really belongs in io_options but that doesn't get
                      ! passed where it is needed
   logical :: errind_up2date, oldsoln_exists, have_true
   character(len=FN_LEN) :: triangle_files
end type grid_type

! data structure for data from triangle files, read and derived

type triangle_data
   type(point), pointer :: vert_coord(:)
   real(my_real), pointer :: vert_bparam(:)
   integer, pointer :: tri_edge(:,:), tri_vert(:,:), tri_neigh(:,:), &
                       edge_tri(:,:), edge_vert(:,:), edge_bmark(:), &
                       vert_tri(:,:), vert_edge(:,:), vert_bmark(:), &
                       vert_master(:)
   integer :: ntri, nedge, nvert
end type triangle_data

type refine_options
   real(my_real) :: inc_factor, reftol, term_energy_err, term_Linf_err, &
                    term_L2_err, t3s_gamma, t3s_eta, t3s_h_target, &
                    t3s_p_target, tp_gamma, sp_gamma_h, sp_gamma_p, &
                    refsoln_pbias
   integer :: error_estimator, reftype, refterm, edge_rule, hp_strategy, &
              max_vert, max_elem, max_dof, max_lev, max_deg, t3s_nunif, &
              t3s_maxref, t3s_maxdeginc, t3s_reftype, nlp_max_h_dec, &
              nlp_max_h_inc, nlp_max_p_dec, nlp_max_p_inc
   logical :: derefine
end type refine_options

! error indicator lists, to determine which elements get refined next during
! adaptive refinement.  List k contains leaf elements in this partition for
! which m/binw**(k-1) > e > m/binw**k where m is max_errind, binw is the
! bin width (2 for cutoffs at 1, 1/2, 1/4, 1/8, ...)  and e is the error
! indicator for the element, with the last set going down to 0.0
! Currently using 16 lists of width 4th root of 2.

type errind_list
   integer :: current_list
   integer :: head_errind(16), tail_errind(16)
   integer, pointer :: next_errind(:), prev_errind(:)
   real(my_real) :: max_errind
end type

!real(my_real), parameter :: binw = sqrt(sqrt(2.0_my_real))
real(my_real), parameter :: binw = 1.1892071150027_my_real
!----------------------------------------------------
! Generic procedures

private get_child_lid_scalar, get_child_lid_array, &
        get_child_gid_scalar, get_child_gid_array

interface get_child_lid
   module procedure get_child_lid_scalar, get_child_lid_array
end interface

interface get_child_gid
   module procedure get_child_gid_scalar, get_child_gid_array
end interface

contains

!        --------------------
function get_child_gid_scalar(gid,child)
!        --------------------

!----------------------------------------------------
! This function returns the global ID of child number child of
! element number gid
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_key), intent(in) :: gid
integer, intent(in) :: child
type(hash_key) :: get_child_gid_scalar

!----------------------------------------------------
! Local variables: 

!----------------------------------------------------
! Begin executable code

get_child_gid_scalar = MAX_CHILD*gid+(child-1)

end function get_child_gid_scalar

!        -------------------
function get_child_gid_array(gid,child)
!        -------------------

!----------------------------------------------------
! This function returns the global ID of the children of element
! number gid listed in child
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_key), intent(in) :: gid
integer, intent(in) :: child(:)
type(hash_key) :: get_child_gid_array(size(child))

!----------------------------------------------------
! Local variables: 

integer :: i

!----------------------------------------------------
! Begin executable code

do i=1,size(child)
   get_child_gid_array(i) = MAX_CHILD*gid+(child(i)-1)
end do

end function get_child_gid_array

!        --------------------
function get_child_lid_scalar(gid,child,table)
!        --------------------

!----------------------------------------------------
! This function returns the local ID of child number child of
! element number gid
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_key), intent(in) :: gid
integer, intent(in) :: child
type(hash_table), intent(in) :: table
integer :: get_child_lid_scalar

!----------------------------------------------------
! Local variables: 

!----------------------------------------------------
! Begin executable code

if (hash_overflow(gid,MAX_CHILD,child-1)) then
   get_child_lid_scalar = NO_CHILD
else
   get_child_lid_scalar = hash_decode_key(MAX_CHILD*gid+(child-1),table)
   if (get_child_lid_scalar == HASH_NOT_FOUND) then
      get_child_lid_scalar = NO_CHILD
   endif
endif

end function get_child_lid_scalar

!        -------------------
function get_child_lid_array(gid,child,table)
!        -------------------

!----------------------------------------------------
! This function returns the local ID of the children of element
! number gid listed in child
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_key), intent(in) :: gid
integer, intent(in) :: child(:)
type(hash_table), intent(in) :: table
integer :: get_child_lid_array(size(child))

!----------------------------------------------------
! Local variables: 

integer :: i

!----------------------------------------------------
! Begin executable code

do i=1,size(child)
   if (hash_overflow(gid,MAX_CHILD,child(i)-1)) then
      get_child_lid_array(i) = NO_CHILD
   else
      get_child_lid_array(i)=hash_decode_key(MAX_CHILD*gid+(child(i)-1), &
                                             table)
      if (get_child_lid_array(i) == HASH_NOT_FOUND) then
         get_child_lid_array(i) = NO_CHILD
      endif
   endif
end do

end function get_child_lid_array

!        -------------
function get_neighbors(lid,grid)
!        -------------

!----------------------------------------------------
! This function returns the local IDs of the neighbors of element lid.
! Order is (sibling,child of parent's mate,mate) i.e. opposite vert 1,2,3.
! On level 1, the order is the same as in initial_neighbor.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: lid
type(grid_type), intent(in) :: grid
integer :: get_neighbors(EDGES_PER_ELEMENT)

!----------------------------------------------------
! Local variables: 

integer :: parent_lid, mate_lid, child_lid
type(hash_key) :: gid, parent, child1, parent_mate_gid, gid_neighbor
type(hash_key) :: children(MAX_CHILD)
integer :: allc(MAX_CHILD)
integer :: i, j, count

!----------------------------------------------------
! Begin executable code

if (grid%element(lid)%level == 1) then
! neighbors of a level 1 element
   get_neighbors = grid%initial_neighbor(:,lid)
   do i=1,EDGES_PER_ELEMENT
! if the initial neighbor was refined, get the child that shares two vertices
! RESTRICTION bisected triangles
      if (get_neighbors(i) /= BOUNDARY) then
         gid_neighbor = grid%element(get_neighbors(i))%gid
         if (hash_overflow(gid_neighbor,2,1)) then
            child_lid = HASH_NOT_FOUND
         else
            allc = ALL_CHILDREN
            children = get_child_gid(gid_neighbor,allc)
            child_lid = hash_decode_key(children(1),grid%elem_hash)
         endif
         if (child_lid /= HASH_NOT_FOUND) then
           count = 0
           do j=1,3
             if (grid%element(child_lid)%vertex(1)==grid%element(lid)%vertex(j).or.&
                 grid%element(child_lid)%vertex(2)==grid%element(lid)%vertex(j).or.&
                 grid%element(child_lid)%vertex(3)==grid%element(lid)%vertex(j))then
                count = count + 1
              endif
            end do
            if (count == 2) then
               get_neighbors(i) = child_lid
            else
               get_neighbors(i) = hash_decode_key(children(2),grid%elem_hash)
            endif
         endif
      endif
   end do
else
! RESTRICTION bisected triangles
   gid = grid%element(lid)%gid
   parent = gid/2
   parent_lid = hash_decode_key(parent,grid%elem_hash)
   child1 = 2*parent
! sibling
   if (child1 == gid) then
      gid_neighbor = child1+1
   else
      gid_neighbor = child1
   endif
! if the sibling was refined, get the child whose vertex 1 is my vertex 2
   if (hash_overflow(gid_neighbor,2,1)) then
      child_lid = HASH_NOT_FOUND
   else
      allc = ALL_CHILDREN
      children = get_child_gid(gid_neighbor,allc)
      child_lid = hash_decode_key(children(1),grid%elem_hash)
   endif
   if (child_lid /= HASH_NOT_FOUND) then
      if (grid%element(child_lid)%vertex(1) == grid%element(lid)%vertex(2)) then
         get_neighbors(1) = child_lid
         gid_neighbor = BOUNDARY ! so I don't hash decode it later
      else
         gid_neighbor = children(2)
      endif
   endif
   if (.not. gid_neighbor == BOUNDARY) then
      get_neighbors(1) = hash_decode_key(gid_neighbor,grid%elem_hash)
   endif
! child of parent's mate
   parent_mate_gid = grid%element(parent_lid)%mate
   if (parent_mate_gid == BOUNDARY) then
      get_neighbors(2) = BOUNDARY
      gid_neighbor = BOUNDARY
   else
      if (child1 == gid) then
         gid_neighbor = 2*parent_mate_gid
      else
         gid_neighbor = 2*parent_mate_gid+1
      endif
! if the neighbor was refined, get the child whose vertex 1 is my vertex 1
      if (hash_overflow(gid_neighbor,2,1)) then
         child_lid = HASH_NOT_FOUND
      else
         allc = ALL_CHILDREN
         children = get_child_gid(gid_neighbor,allc)
         child_lid = hash_decode_key(children(1),grid%elem_hash)
      endif
      if (child_lid /= HASH_NOT_FOUND) then
         if (grid%element(child_lid)%vertex(1)==grid%element(lid)%vertex(1)) then
            get_neighbors(2) = child_lid
            gid_neighbor = BOUNDARY
         else
            gid_neighbor = children(2)
         endif
      endif
   endif
   if (.not. gid_neighbor == BOUNDARY) then
      get_neighbors(2) = hash_decode_key(gid_neighbor,grid%elem_hash)
   endif
! mate
   if (grid%element(lid)%mate == BOUNDARY) then
      get_neighbors(3) = BOUNDARY
   else
      mate_lid = hash_decode_key(grid%element(lid)%mate,grid%elem_hash)
      if (mate_lid == HASH_NOT_FOUND) then
         get_neighbors(3)=hash_decode_key(grid%element(lid)%mate/2,grid%elem_hash)
      else
         get_neighbors(3) = mate_lid
      endif
   endif
endif

end function get_neighbors

!          -------------
subroutine get_grid_info(grid,procs,still_sequential,tag,nelem,nlev,nvert, &
                         nvert_own,nelem_own,nelem_leaf,nelem_leaf_own,     &
                         dof,dof_own,total_nvert,total_nelem_leaf,total_dof, &
                         max_nlev,mindeg,maxdeg,no_master)
!          -------------

!----------------------------------------------------
! This routine returns information about the grid
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

type (grid_type), intent(in) :: grid
type (proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential
integer, intent(in) :: tag
integer, intent(out), optional :: nelem,nlev,nvert,nvert_own, &
                                  nelem_own,nelem_leaf,nelem_leaf_own, &
                                  dof,dof_own,total_nvert,total_nelem_leaf, &
                                  total_dof,max_nlev,mindeg,maxdeg
logical, intent(in), optional :: no_master
!----------------------------------------------------

integer :: proc,ni,nr,astat,lev,elem,min_mindeg,max_maxdeg
integer, pointer :: recv_int(:)
real(my_real), pointer :: recv_real(:)
real(my_real) :: no_reals(1)
logical :: send_to_master

!----------------------------------------------------
! Begin executable code

if (present(no_master)) then
   if (no_master) then
      send_to_master = .false.
   else
      send_to_master = .true.
   endif
else
   send_to_master = .true.
endif

if (my_proc(procs)==MASTER) then

   if (present(nelem)) nelem = 0
   if (present(nlev)) nlev = 0
   if (present(nvert)) nvert = 0
   if (present(nvert_own)) nvert_own = 0
   if (present(nelem_own)) nelem_own = 0
   if (present(nelem_leaf)) nelem_leaf = 0
   if (present(nelem_leaf_own)) nelem_leaf_own = 0
   if (present(dof)) dof = 0
   if (present(dof_own)) dof_own = 0
   if (present(total_nvert)) then
      call phaml_recv(procs,proc,recv_int,ni,recv_real,nr,10*tag+4)
      total_nvert = recv_int(1)
      deallocate(recv_int,stat=astat)
   end if
   if (present(total_nelem_leaf)) then
      call phaml_recv(procs,proc,recv_int,ni,recv_real,nr,10*tag+6)
      total_nelem_leaf = recv_int(1)
      deallocate(recv_int,stat=astat)
   end if
   if (present(total_dof)) then
      call phaml_recv(procs,proc,recv_int,ni,recv_real,nr,10*tag+9)
      total_dof = recv_int(1)
      deallocate(recv_int,stat=astat)
   end if
   if (present(max_nlev)) then
      call phaml_recv(procs,proc,recv_int,ni,recv_real,nr,10*tag+8)
      max_nlev = recv_int(1)
      deallocate(recv_int,stat=astat)
   end if
   if (present(mindeg)) then
      call phaml_recv(procs,proc,recv_int,ni,recv_real,nr,10*tag+11)
      mindeg = recv_int(1)
      deallocate(recv_int,stat=astat)
   end if
   if (present(maxdeg)) then
      call phaml_recv(procs,proc,recv_int,ni,recv_real,nr,10*tag+13)
      maxdeg = recv_int(1)
      deallocate(recv_int,stat=astat)
   end if

else ! not processor 0

   if (present(nelem)) nelem = grid%nelem
   if (present(nlev)) nlev = grid%nlev
   if (present(nvert)) nvert = grid%nvert
   if (present(nvert_own)) nvert_own = grid%nvert_own
   if (present(nelem_own)) nelem_own = 0 ! was grid%nelem_own
   if (present(nelem_leaf)) nelem_leaf = grid%nelem_leaf
   if (present(nelem_leaf_own)) nelem_leaf_own = grid%nelem_leaf_own
   if (present(dof)) dof = grid%dof
   if (present(dof_own)) dof_own = grid%dof_own

   if (present(total_nvert)) then
      if (still_sequential) then
         total_nvert = grid%nvert_own
      else
         total_nvert = phaml_global_sum(procs,grid%nvert_own,10*tag+1)
      endif
      if (my_proc(procs)==1 .and. send_to_master) then
         call phaml_send(procs,MASTER,(/total_nvert/),1,no_reals,0,10*tag+4)
      endif
   end if
   if (present(total_nelem_leaf)) then
      if (still_sequential) then
         total_nelem_leaf = grid%nelem_leaf_own
      else
         total_nelem_leaf = phaml_global_sum(procs,grid%nelem_leaf_own,10*tag+3)
      endif
      if (my_proc(procs)==1 .and. send_to_master) then
         call phaml_send(procs,MASTER,(/total_nelem_leaf/),1,no_reals,0,10*tag+6)
      endif
   end if
   if (present(total_dof)) then
      if (still_sequential) then
         total_dof = grid%dof_own
      else
         total_dof = phaml_global_sum(procs,grid%dof_own,10*tag+2)
      endif
      if (my_proc(procs)==1 .and. send_to_master) then
         call phaml_send(procs,MASTER,(/total_dof/),1,no_reals,0,10*tag+9)
      endif
   end if
   if (present(max_nlev)) then
      if (still_sequential) then
         max_nlev = grid%nlev
      else
         max_nlev = phaml_global_max(procs,grid%nlev,10*tag+7)
      endif
      if (my_proc(procs)==1 .and. send_to_master) then
         call phaml_send(procs,MASTER,(/max_nlev/),1,no_reals,0,10*tag+8)
      endif
   endif
   if (present(mindeg) .or. present(maxdeg)) then
      if (present(mindeg)) mindeg = huge(0)
      if (present(maxdeg)) maxdeg = 0
      do lev=1,grid%nlev
         elem = grid%head_level_elem(lev)
         do while (elem /= END_OF_LIST)
            if (grid%element(elem)%iown .and. grid%element(elem)%isleaf) then
               if (present(mindeg)) mindeg=min(mindeg,grid%element(elem)%degree)
               if (present(maxdeg)) maxdeg=max(maxdeg,grid%element(elem)%degree)
            endif
            elem = grid%element(elem)%next
         end do
      end do
      if (send_to_master) then
         if (present(mindeg)) then
            min_mindeg = phaml_global_min(procs,mindeg,10*tag+10)
            if (my_proc(procs)==1) then
               call phaml_send(procs,MASTER,(/min_mindeg/),1,no_reals,0,10*tag+11)
            endif
         endif
         if (present(maxdeg)) then
            max_maxdeg = phaml_global_max(procs,maxdeg,10*tag+12)
            if (my_proc(procs)==1) then
               call phaml_send(procs,MASTER,(/max_maxdeg/),1,no_reals,0,10*tag+13)
            endif
         endif
      endif
   endif

endif

end subroutine get_grid_info

end module gridtype_mod
