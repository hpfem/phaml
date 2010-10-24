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

module linsys_util

!----------------------------------------------------
! This module contains utilities used by linear_system and related modules
!
! communication tags in this module are of the form 12xx
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
use sort_mod

!----------------------------------------------------

implicit none
private
public basis_change, exchange_neigh_vect, exchange_fudop_vect, &
       send_neigh_vect, recv_neigh_vect, make_need_r_others, &
       exchange_fudop_soln_residual, sum_fudop_vect, all_to_mine, mine_to_all, &
       grid_to_eq, eq_to_grid, matrix_times_vector, linsys_residual, &
       test_precon, TO_NODAL, TO_HIER, VERTEX_ID, EDGE_ID, ELEMENT_ID

!----------------------------------------------------
! Parameters

! gives the direction of a basis change
integer, parameter :: TO_NODAL = 1, TO_HIER = 2

! the shift (2^bits) for defining the rank field of an equation gid
integer, parameter :: SHIFT = 4

! specifies the object type for equation gids

integer, parameter :: VERTEX_ID  = 1, &
                      EDGE_ID    = 2, &
                      ELEMENT_ID = 3

contains

!          ----------
subroutine grid_to_eq(grid,linear_system,object_type,basis_rank, &
                      system_rank,grid_gid,eqn_lid,eqn_gid)
!          ----------

!----------------------------------------------------
! This routine takes a vertex, edge or element gid, a rank within that
! object, and a rank within a system of PDEs and returns an equation lid
! and/or gid.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
type(linsys_type), intent(in) :: linear_system
integer, intent(in) :: object_type, basis_rank, system_rank
type(hash_key), intent(in) :: grid_gid
integer, intent(out), optional :: eqn_lid
type(hash_key_eq), intent(out), optional :: eqn_gid
!----------------------------------------------------
! Local variables:

type(hash_key_eq) :: loc_gid
type(hash_key) :: gid
integer :: array_gid(KEY_SIZE+1), lid
!----------------------------------------------------
! Begin executable code

gid = grid_gid
if (object_type == VERTEX_ID .or. object_type == EDGE_ID) then
   select case (object_type)
   case (VERTEX_ID)
      lid = hash_decode_key(grid_gid,grid%vert_hash)
      if (lid /= HASH_NOT_FOUND) then
         if (grid%vertex_type(lid,system_rank) == PERIODIC_SLAVE .or. &
             grid%vertex_type(lid,system_rank) == PERIODIC_SLAVE_DIR .or. &
             grid%vertex_type(lid,system_rank) == PERIODIC_SLAVE_NAT .or. &
             grid%vertex_type(lid,system_rank) == PERIODIC_SLAVE_MIX) then
            gid = grid%vertex(grid%vertex(lid)%next)%gid
         endif
      endif
   case (EDGE_ID)
      lid = hash_decode_key(grid_gid,grid%edge_hash)
      if (lid /= HASH_NOT_FOUND) then
         if (grid%edge_type(lid,system_rank) == PERIODIC_SLAVE) then
            gid = grid%edge(grid%edge(lid)%next)%gid
         endif
      endif
   end select
endif

call hash_pack_key(gid,array_gid,1)
array_gid(KEY_SIZE+1) = object_type + &
                   SHIFT*(basis_rank*(linear_system%system_size+1)+system_rank)
loc_gid = hash_unpack_key(array_gid,1,.true.)
if (present(eqn_gid)) eqn_gid = loc_gid
if (present(eqn_lid)) eqn_lid = hash_decode_key(loc_gid,linear_system%eq_hash)

end subroutine grid_to_eq

!          ----------
subroutine eq_to_grid(linear_system,eqn_gid,object_type,basis_rank, &
                      system_rank,grid_lid,grid,grid_gid)
!          ----------

!----------------------------------------------------
! This routine takes an equation gid and returns the vertex, edge or
! element lid and/or gid, rank on that object of the corresponding basis
! function, and rank in a system of PDEs.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(in) :: linear_system
type(hash_key_eq), intent(in) :: eqn_gid
integer, intent(out) :: object_type, basis_rank, system_rank
integer, intent(out), optional :: grid_lid
type(grid_type), intent(in), optional :: grid
type(hash_key), intent(out), optional :: grid_gid
!----------------------------------------------------
! Local variables:

type(hash_key) :: loc_gid
integer :: array_gid(KEY_SIZE+1), ranks
!----------------------------------------------------
! Begin executable code

if (present(grid_lid) .neqv. present(grid)) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("one of grid_lid or grid is present without the other in eq_to_gid")
   stop
endif
call hash_pack_key(eqn_gid,array_gid,1)
ranks = array_gid(KEY_SIZE+1)/SHIFT
object_type = array_gid(KEY_SIZE+1) - ranks*SHIFT
basis_rank = ranks/(linear_system%system_size+1)
system_rank = ranks - basis_rank*(linear_system%system_size+1)
loc_gid = hash_unpack_key(array_gid,1)
if (present(grid_gid)) grid_gid = loc_gid
if (present(grid_lid)) then
   select case(object_type)
   case (VERTEX_ID)
      grid_lid = hash_decode_key(loc_gid,grid%vert_hash)
   case (EDGE_ID)
      grid_lid = hash_decode_key(loc_gid,grid%edge_hash)
   case (ELEMENT_ID)
      grid_lid = hash_decode_key(loc_gid,grid%elem_hash)
   case default
      ierr = PHAML_INTERNAL_ERROR
      call fatal("illegal object type in eq_to_grid",intlist=(/object_type/))
      grid_lid = HASH_NOT_FOUND
   end select
endif

end subroutine eq_to_grid

!          ------------
subroutine basis_change(lev,direction,phaml_matrix,skip_matrix,do_mass)
!          ------------

!----------------------------------------------------
! This routine converts the linear subsystem of level 1 through lev between
! nodal and the two level hierarchical basis.  direction is either TO_NODAL
! or TO_HIER to determine which way to convert.  It converts matrix_val,
! the column indices, and the vectors rhs, r_mine, r_others and solution.
! If skip_matrix is present and .true., then only the vectors and column
! indices are converted, not the matrix values.
! If do_mass is present and .true. and the mass matrix exists, it is also
! changed.
! RESTRICTION linear bases
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: lev, direction
type(linsys_type), intent(inout), target :: phaml_matrix
logical, intent(in), optional :: skip_matrix, do_mass
!----------------------------------------------------
! Local variables:

integer :: ss, ssm1, i, j, v, a1, a2, p1, p2, vrow, a1row, a2row, p1row, &
           p2row, va1, va2, vp1, vp2, a1a2, a1p1, a1p2, a1v, a2a1, a2p1, a2p2, &
           a2v, p1a1, p1a2, p1v, p2a1, p2a2, p2v, vrowend, a1rowend, a2rowend, &
           p1rowend, p2rowend, vrowlen, a1rowlen, a2rowlen, p1rowlen, &
           p2rowlen, vsep, a1sep, a2sep, p1sep, p2sep, noent1, noent2
logical :: do_matrix, boundary_vert, loc_do_mass
!----------------------------------------------------
! Begin executable code

! for p-hierarchical high order bases, just change the end of row between
! the whole row and end of edge or linear bases

if (lev > phaml_matrix%nlev) then
   select case (direction)
   case (TO_HIER)
      if (lev == phaml_matrix%nlev+2) then
         phaml_matrix%end_row => phaml_matrix%end_row_edge
      elseif (lev == phaml_matrix%nlev+1) then
         phaml_matrix%end_row => phaml_matrix%end_row_linear
      endif
   case (TO_NODAL)
      if (lev == phaml_matrix%nlev+2) then
         phaml_matrix%end_row => phaml_matrix%end_row_face
      elseif (lev == phaml_matrix%nlev+1) then
         phaml_matrix%end_row => phaml_matrix%end_row_edge
      endif
   end select
   return
endif

! for systems of equations, only work with the first component.  To get at
! all the components, use blocks of size system size

! some useful constants

ss = phaml_matrix%system_size
ssm1 = ss-1
do_matrix = .true.
if (present(skip_matrix)) then
   if (skip_matrix) then
      do_matrix = .false.
   endif
endif
if (present(do_mass)) then
   if (do_mass .and. associated(phaml_matrix%mass)) then
      loc_do_mass = .true.
   else
      loc_do_mass = .false.
   endif
else
   loc_do_mass = .false.
endif

! v is the beginning of a block of rows associated with a red vertex

do v = phaml_matrix%begin_level(lev), phaml_matrix%begin_level(lev+1)-1, ss

! determine if the vertex is a boundary vertex

   select case (phaml_matrix%equation_type(v))
   case (INTERIOR, PERIODIC, PERIODIC_MASTER, PERIODIC_SLAVE)
      boundary_vert = .false.
   case (DIRICHLET, NATURAL, MIXED)
      boundary_vert = .true.
   case default
      call fatal("unrecognized equation_type in basis_change")
   end select

! start of row for red vertex

   vrow = phaml_matrix%begin_row(v)

! identify local ids of parents and ancestors

   a1 = phaml_matrix%column_index(vrow+  ss)
   a2 = phaml_matrix%column_index(vrow+2*ss)
   p1 = phaml_matrix%column_index(vrow+3*ss)
   if (boundary_vert) then
      p2 = 0
   else
      p2 = phaml_matrix%column_index(vrow+4*ss)
   endif

! start of row for each neighbor

   a1row = phaml_matrix%begin_row(a1)
   a2row = phaml_matrix%begin_row(a2)
   p1row = phaml_matrix%begin_row(p1)
   if (.not. boundary_vert) then
      p2row = phaml_matrix%begin_row(p2)
   endif

! ends of rows and length of rows

   vrowend  = phaml_matrix%end_row(v)
   a1rowend = phaml_matrix%end_row(a1)
   a2rowend = phaml_matrix%end_row(a2)
   p1rowend = phaml_matrix%end_row(p1)
   if (.not. boundary_vert) then
      p2rowend = phaml_matrix%end_row(p2)
   endif

   vrowlen  = vrowend  - vrow  + 1
   a1rowlen = a1rowend - a1row + 1
   a2rowlen = a2rowend - a2row + 1
   p1rowlen = p1rowend - p1row + 1
   if (.not. boundary_vert) then
      p2rowlen = p2rowend - p2row + 1
   endif

! convert vectors

   if (direction==TO_HIER) then
      phaml_matrix%rhs(a1:a1+ssm1) = phaml_matrix%rhs(a1:a1+ssm1) + &
                                   phaml_matrix%rhs(v:v+ssm1)/2
      phaml_matrix%rhs(a2:a2+ssm1) = phaml_matrix%rhs(a2:a2+ssm1) + &
                                   phaml_matrix%rhs(v:v+ssm1)/2
      phaml_matrix%r_mine(a1:a1+ssm1) = phaml_matrix%r_mine(a1:a1+ssm1) + &
                                      phaml_matrix%r_mine(v:v+ssm1)/2
      phaml_matrix%r_mine(a2:a2+ssm1) = phaml_matrix%r_mine(a2:a2+ssm1) + &
                                      phaml_matrix%r_mine(v:v+ssm1)/2
      phaml_matrix%r_others(a1:a1+ssm1) = phaml_matrix%r_others(a1:a1+ssm1) + &
                                        phaml_matrix%r_others(v:v+ssm1)/2
      phaml_matrix%r_others(a2:a2+ssm1) = phaml_matrix%r_others(a2:a2+ssm1) + &
                                        phaml_matrix%r_others(v:v+ssm1)/2
      phaml_matrix%solution(v:v+ssm1) = phaml_matrix%solution(v:v+ssm1) -&
            (phaml_matrix%solution(a1:a1+ssm1)+phaml_matrix%solution(a2:a2+ssm1))/2
   else
      phaml_matrix%rhs(a1:a1+ssm1) = phaml_matrix%rhs(a1:a1+ssm1) - &
                                   phaml_matrix%rhs(v:v+ssm1)/2
      phaml_matrix%rhs(a2:a2+ssm1) = phaml_matrix%rhs(a2:a2+ssm1) - &
                                   phaml_matrix%rhs(v:v+ssm1)/2
      phaml_matrix%r_mine(a1:a1+ssm1) = phaml_matrix%r_mine(a1:a1+ssm1) - &
                                      phaml_matrix%r_mine(v:v+ssm1)/2
      phaml_matrix%r_mine(a2:a2+ssm1) = phaml_matrix%r_mine(a2:a2+ssm1) - &
                                      phaml_matrix%r_mine(v:v+ssm1)/2
      phaml_matrix%r_others(a1:a1+ssm1) = phaml_matrix%r_others(a1:a1+ssm1) - &
                                        phaml_matrix%r_others(v:v+ssm1)/2
      phaml_matrix%r_others(a2:a2+ssm1) = phaml_matrix%r_others(a2:a2+ssm1) - &
                                        phaml_matrix%r_others(v:v+ssm1)/2
      phaml_matrix%solution(v:v+ssm1) = phaml_matrix%solution(v:v+ssm1) +&
            (phaml_matrix%solution(a1:a1+ssm1)+phaml_matrix%solution(a2:a2+ssm1))/2
   endif

! convert matrix

! In general, xy is the location for vertex y in the row of vertex x.
! If direction is TO_HIER, a1a2 and a2a1 give the location to start placing
! the new entries for ancestor-ancestor inner products.
! If direction is TO_NODAL, then a1v, a2v, p1v and p2v give the location
! to start placing the new entries for the red vertex inner products.

! For those entries that are new after the basis change, identify where
! they should go in the rows to keep indices increasing across the rows,
! and the location of existing or new empty space, and shift the row accordingly

   if (direction == TO_HIER) then

! location to place a2 in row of a1
      a1a2 = 0
      do i=ss+1,a1rowlen,ss
         if (phaml_matrix%column_index(a1row+i-1) > a2) then
            a1a2 = a1row+i-1
            exit
         endif
      end do
      if (a1a2 == 0) then
         ierr = PHAML_INTERNAL_ERROR
         call fatal("basis_change:failed to find place for ancestor 2 in ancestor 1")
         stop
      endif

! location of v in a1, now considered empty space
      a1v = 0
      do i=2*ss+1,a1rowlen,ss
         if (phaml_matrix%column_index(a1row+i-1) == v) then
            a1v = a1row+i-1
            exit
         endif
      end do
      if (a1v == 0) then
         ierr = PHAML_INTERNAL_ERROR
         call fatal("basis_change:failed to find red vertex in ancestor 1")
         stop
      endif

! location to place a1 in row of a2
      a2a1 = 0
      do i=ss+1,a2rowlen,ss
         if (phaml_matrix%column_index(a2row+i-1) > a1) then
            a2a1 = a2row+i-1
            exit
         endif
      end do
      if (a2a1 == 0) then
         ierr = PHAML_INTERNAL_ERROR
         call fatal("basis_change:failed to find place for ancestor 2 in ancestor 1")
         stop
      endif

! location of v in a2, now considered empty space
      a2v = 0
      do i=2*ss+1,a2rowlen,ss
         if (phaml_matrix%column_index(a2row+i-1) == v) then
            a2v = a2row+i-1
            exit
         endif
      end do
      if (a2v == 0) then
         ierr = PHAML_INTERNAL_ERROR
         call fatal("basis_change:failed to find red vertex in ancestor 2")
         stop
      endif

! perform the shift
      do i=1,ss
         vsep  = phaml_matrix%begin_row(v +i-1) - phaml_matrix%begin_row(v )
         a1sep = phaml_matrix%begin_row(a1+i-1) - phaml_matrix%begin_row(a1)
         a2sep = phaml_matrix%begin_row(a2+i-1) - phaml_matrix%begin_row(a2)
         p1sep = phaml_matrix%begin_row(p1+i-1) - phaml_matrix%begin_row(p1)
         if (.not. boundary_vert) then
          p2sep = phaml_matrix%begin_row(p2+i-1) - phaml_matrix%begin_row(p2)
         endif
         phaml_matrix%column_index(a1a2+a1sep+ss:a1v+a1sep+ssm1) = &
            phaml_matrix%column_index(a1a2+a1sep:a1v+a1sep-1)
         phaml_matrix%column_index(a1a2+a1sep:a1a2+a1sep+ssm1) = &
            (/ ((a2+j),j=0,ssm1) /)
         phaml_matrix%column_index(a2a1+a2sep+ss:a2v+a2sep+ssm1) = &
         phaml_matrix%column_index(a2a1+a2sep:a2v+a2sep-1)
         phaml_matrix%column_index(a2a1+a2sep:a2a1+a2sep+ssm1) = &
            (/ ((a1+j),j=0,ssm1) /)
         if (do_matrix) then
            phaml_matrix%matrix_val(a1a2+a1sep+ss:a1v+a1sep+ssm1) = &
               phaml_matrix%matrix_val(a1a2+a1sep:a1v+a1sep-1)
            phaml_matrix%matrix_val(a2a1+a2sep+ss:a2v+a2sep+ssm1) = &
               phaml_matrix%matrix_val(a2a1+a2sep:a2v+a2sep-1)
            if (loc_do_mass) then
               phaml_matrix%mass(a1a2+a1sep+ss:a1v+a1sep+ssm1) = &
                  phaml_matrix%mass(a1a2+a1sep:a1v+a1sep-1)
               phaml_matrix%mass(a2a1+a2sep+ss:a2v+a2sep+ssm1) = &
                  phaml_matrix%mass(a2a1+a2sep:a2v+a2sep-1)
            endif
         endif
      end do

   else ! direction == TO_NODAL

! location to place v in row of a1
      a1v = 0
      do i=2*ss+1,a1rowlen,ss
         if (phaml_matrix%column_index(a1row+i-1) == NO_ENTRY .and. &
             a1v == 0) then
            a1v = a1row+i-1
         endif
         if (phaml_matrix%column_index(a1row+i-1) > v) then
            a1v = a1row+i-1-ss
            exit
         endif
      end do
      if (a1v == 0) then
         a1v = a1rowend - ss + 1
      endif

! location of a2 in a1, now considered empty space
      a1a2 = 0
      do i=ss+1,a1rowlen,ss
         if (phaml_matrix%column_index(a1row+i-1) == a2) then
            a1a2 = a1row+i-1
            exit
         endif
      end do
      if (a1a2 == 0) then
         ierr = PHAML_INTERNAL_ERROR
         call fatal("basis_change:failed to find ancestor 2 in ancestor 1")
         stop
      endif

! location to place v in row of a2
      a2v = 0
      do i=2*ss+1,a2rowlen,ss
         if (phaml_matrix%column_index(a2row+i-1) == NO_ENTRY .and. &
             a2v == 0) then
            a2v = a2row+i-1
         endif
         if (phaml_matrix%column_index(a2row+i-1) > v) then
            a2v = a2row+i-1-ss
            exit
         endif
      end do
      if (a2v == 0) then
         a2v = a2rowend - ss + 1
      endif

! location of a2 in a1, now considered empty space
      a2a1 = 0
      do i=ss+1,a2rowlen,ss
         if (phaml_matrix%column_index(a2row+i-1) == a1) then
            a2a1 = a2row+i-1
            exit
         endif
      end do
      if (a2a1 == 0) then
         ierr = PHAML_INTERNAL_ERROR
         call fatal("basis_change:failed to find ancestor 1 in ancestor 2")
         stop
      endif

! location to place v in row of p1
      p1v = 0
      do i=3*ss+1,p1rowlen,ss
         if (phaml_matrix%column_index(p1row+i-1) == NO_ENTRY .and. &
             p1v == 0) then
            p1v = p1row+i-1
         endif
         if (phaml_matrix%column_index(p1row+i-1) > v) then
            p1v = p1row+i-1-ss
            exit
         endif
      end do
      if (phaml_matrix%column_index(p1v) == NO_ENTRY) then
         noent1 = 0
      else
         noent1 = p1v+1
         do while (phaml_matrix%column_index(p1v) /= NO_ENTRY)
            noent1 = noent1 + 1
         end do
      endif

! location to place v in row of p2
      if (.not. boundary_vert) then
         p2v = 0
         do i=3*ss+1,p2rowlen,ss
            if (phaml_matrix%column_index(p2row+i-1)==NO_ENTRY .and. &
                p2v==0) then
               p2v = p2row+i-1
            endif
            if (phaml_matrix%column_index(p2row+i-1) > v) then
               p2v = p2row+i-1-ss
               exit
            endif
         end do
         if (phaml_matrix%column_index(p2v) == NO_ENTRY) then
            noent2 = 0
         else
            noent2 = p2v+1
            do while (phaml_matrix%column_index(p2v) /= NO_ENTRY)
               noent2 = noent2 + 1
            end do
         endif
      endif

! perform the shifts
      do i=1,ss
         vsep  = phaml_matrix%begin_row(v +i-1) - phaml_matrix%begin_row(v )
         a1sep = phaml_matrix%begin_row(a1+i-1) - phaml_matrix%begin_row(a1)
         a2sep = phaml_matrix%begin_row(a2+i-1) - phaml_matrix%begin_row(a2)
         p1sep = phaml_matrix%begin_row(p1+i-1) - phaml_matrix%begin_row(p1)
         if (.not. boundary_vert) then
            p2sep = phaml_matrix%begin_row(p2+i-1) - phaml_matrix%begin_row(p2)
         endif
         if (do_matrix) then
            phaml_matrix%matrix_val(a1a2+a1sep:a1v+a1sep-1) = &
               phaml_matrix%matrix_val(a1a2+a1sep+ss:a1v+a1sep+ssm1)
            if (loc_do_mass) then
               phaml_matrix%mass(a1a2+a1sep:a1v+a1sep-1) = &
                  phaml_matrix%mass(a1a2+a1sep+ss:a1v+a1sep+ssm1)
            endif
         endif
         phaml_matrix%column_index(a1a2+a1sep:a1v+a1sep-1) = &
            phaml_matrix%column_index(a1a2+a1sep+ss:a1v+a1sep+ssm1)
         if (do_matrix) then
            phaml_matrix%matrix_val(a2a1+a2sep:a2v+a2sep-1) = &
               phaml_matrix%matrix_val(a2a1+a2sep+ss:a2v+a2sep+ssm1)
            if (loc_do_mass) then
               phaml_matrix%mass(a2a1+a2sep:a2v+a2sep-1) = &
                  phaml_matrix%mass(a2a1+a2sep+ss:a2v+a2sep+ssm1)
            endif
         endif
         phaml_matrix%column_index(a2a1+a2sep:a2v+a2sep-1) = &
            phaml_matrix%column_index(a2a1+a2sep+ss:a2v+a2sep+ssm1)
         if (phaml_matrix%column_index(p1v) /= NO_ENTRY) then
            if (do_matrix) then
               phaml_matrix%matrix_val(p1v+p1sep+ss:noent1+p1sep+ssm1) = &
                  phaml_matrix%matrix_val(p1v+p1sep:noent1+p1sep-1)
               if (loc_do_mass) then
                  phaml_matrix%mass(p1v+p1sep+ss:noent1+p1sep+ssm1) = &
                     phaml_matrix%mass(p1v+p1sep:noent1+p1sep-1)
               endif
            endif
            phaml_matrix%column_index(p1v+p1sep+ss:noent1+p1sep+ssm1) = &
               phaml_matrix%column_index(p1v+p1sep:noent1+p1sep-1)
         endif
         if (.not. boundary_vert) then
            if (phaml_matrix%column_index(p2v) /= NO_ENTRY) then
               if (do_matrix) then
                  phaml_matrix%matrix_val(p2v+p2sep+ss:noent2+p2sep+ssm1) = &
                     phaml_matrix%matrix_val(p2v+p2sep:noent2+p2sep-1)
                  if (loc_do_mass) then
                     phaml_matrix%mass(p2v+p2sep+ss:noent2+p2sep+ssm1) = &
                        phaml_matrix%mass(p2v+p2sep:noent2+p2sep-1)
                  endif
               endif
               phaml_matrix%column_index(p2v+p2sep+ss:noent2+p2sep+ssm1) = &
                  phaml_matrix%column_index(p2v+p2sep:noent2+p2sep-1)
            endif
         endif
      end do

   endif ! direction

! end of shifting to make space for new entries

! Identify the start of each vertex in each of the other vertices rows when
! they already exist

! location of other vertices in row of v
   va1 = vrow +   ss
   va2 = vrow + 2*ss
   vp1 = vrow + 3*ss
   vp2 = vrow + 4*ss

   if (do_matrix) then

! location of p1 in row of a1
      a1p1 = 0
      do i=ss+1,a1rowlen,ss
         if (phaml_matrix%column_index(a1row+i-1) == p1) then
            a1p1 = a1row+i-1
            exit
         endif
      end do
      if (a1p1 == 0) then
         ierr = PHAML_INTERNAL_ERROR
         call fatal("basis_change:failed to find parent 1 in ancestor 1")
         stop
      endif

! location of p2 in row of a1
      if (.not. boundary_vert) then
         a1p2 = 0
         do i=2*ss+1,a1rowlen,ss
            if (phaml_matrix%column_index(a1row+i-1) == p2) then
               a1p2 = a1row+i-1
               exit
            endif
         end do
         if (a1p2 == 0) then
            ierr = PHAML_INTERNAL_ERROR
            call fatal("basis_change:failed to find parent 2 in ancestor 1")
            stop
         endif
      endif

! location of p1 in row of a2
      a2p1 = 0
      do i=ss+1,a2rowlen,ss
         if (phaml_matrix%column_index(a2row+i-1) == p1) then
            a2p1 = a2row+i-1
            exit
         endif
      end do
      if (a2p1 == 0) then
         ierr = PHAML_INTERNAL_ERROR
         call fatal("basis_change:failed to find parent 1 in ancestor 2")
         stop
      endif

! location of p2 in row of a2
      if (.not. boundary_vert) then
         a2p2 = 0
         do i=2*ss+1,a2rowlen,ss
            if (phaml_matrix%column_index(a2row+i-1) == p2) then
               a2p2 = a2row+i-1
               exit
            endif
         end do
         if (a2p2 == 0) then
            ierr = PHAML_INTERNAL_ERROR
            call fatal("basis_change:failed to find parent 2 in ancestor 2")
            stop
         endif
      endif

! location of a1 in row of p1
      p1a1 = 0
      do i=ss+1,p1rowlen,ss
         if (phaml_matrix%column_index(p1row+i-1) == a1) then
            p1a1 = p1row+i-1
            exit
         endif
      end do
      if (p1a1 == 0) then
         ierr = PHAML_INTERNAL_ERROR
         call fatal("basis_change:failed to find ancestor 1 in parent 1")
         stop
      endif

! location of a2 in row of p1
      p1a2 = 0
      do i=ss+1,p1rowlen,ss
         if (phaml_matrix%column_index(p1row+i-1) == a2) then
            p1a2 = p1row+i-1
            exit
         endif
      end do
      if (p1a2 == 0) then
         ierr = PHAML_INTERNAL_ERROR
         call fatal("basis_change:failed to find ancestor 2 in parent 1")
         stop
      endif

   endif ! do_matrix

! location of v in p1
   if (direction==TO_HIER) then
      p1v = 0
      do i=3*ss+1,p1rowlen,ss
         if (phaml_matrix%column_index(p1row+i-1) == v) then
            p1v = p1row+i-1
            exit
         endif
      end do
      if (p1v == 0) then
         ierr = PHAML_INTERNAL_ERROR
         call fatal("basis_change:failed to find red vertex in parent 1")
         stop
      endif
   endif

   if (.not. boundary_vert) then

      if (do_matrix) then

! location of a1 in row of p2
         p2a1 = 0
         do i=ss+1,p2rowlen,ss
            if (phaml_matrix%column_index(p2row+i-1) == a1) then
               p2a1 = p2row+i-1
               exit
            endif
         end do
         if (p2a1 == 0) then
            ierr = PHAML_INTERNAL_ERROR
            call fatal("basis_change:failed to find ancestor 1 in parent 2")
            stop
         endif

! location of a2 in row of p2
         p2a2 = 0
         do i=ss+1,p2rowlen,ss
            if (phaml_matrix%column_index(p2row+i-1) == a2) then
               p2a2 = p2row+i-1
               exit
            endif
         end do
         if (p2a2 == 0) then
            ierr = PHAML_INTERNAL_ERROR
            call fatal("basis_change:failed to find ancestor 2 in parent 2")
            stop
         endif

      endif ! do_matrix

! location of v in row of p2
      if (direction==TO_HIER) then
         p2v = 0
         do i=3*ss+1,p2rowlen,ss
            if (phaml_matrix%column_index(p2row+i-1) == v) then
               p2v = p2row+i-1
               exit
            endif
         end do
      endif
      if (p2v == 0) then
         ierr = PHAML_INTERNAL_ERROR
         call fatal("basis_change:failed to find red vertex in parent 2")
         stop
      endif
   endif

! end of finding each vertex in the rows of the other vertices

! transform matrix values

   if (direction == TO_HIER) then

      do i=1,ss
         vsep  = phaml_matrix%begin_row(v +i-1) - phaml_matrix%begin_row(v )
         a1sep = phaml_matrix%begin_row(a1+i-1) - phaml_matrix%begin_row(a1)
         a2sep = phaml_matrix%begin_row(a2+i-1) - phaml_matrix%begin_row(a2)
         p1sep = phaml_matrix%begin_row(p1+i-1) - phaml_matrix%begin_row(p1)
         if (.not. boundary_vert) then
          p2sep = phaml_matrix%begin_row(p2+i-1) - phaml_matrix%begin_row(p2)
         endif
! <a1,a1>
         if (do_matrix) then
            phaml_matrix%matrix_val(a1row+a1sep:a1row+a1sep+ssm1) = &
               phaml_matrix%matrix_val(a1row+a1sep:a1row+a1sep+ssm1) + &
               phaml_matrix%matrix_val(vrow+vsep:vrow+vsep+ssm1)/4
            do j=1,ss
               if (j==1) then
                  phaml_matrix%matrix_val(a1row+a1sep+j-1) = &
                     phaml_matrix%matrix_val(a1row+a1sep+j-1) + &
                     phaml_matrix%matrix_val(va1+vsep+i-1)
               elseif (j==i) then
                  phaml_matrix%matrix_val(a1row+a1sep+j-1) = &
                     phaml_matrix%matrix_val(a1row+a1sep+j-1) + &
                     phaml_matrix%matrix_val(va1+vsep)
               else
                  phaml_matrix%matrix_val(a1row+a1sep+j-1) = &
                     phaml_matrix%matrix_val(a1row+a1sep+j-1) + &
                     phaml_matrix%matrix_val(va1+vsep+j-1)
               endif
            end do
! <a2,a2>
            phaml_matrix%matrix_val(a2row+a2sep:a2row+a2sep+ssm1) = &
               phaml_matrix%matrix_val(a2row+a2sep:a2row+a2sep+ssm1) + &
               phaml_matrix%matrix_val(vrow+vsep:vrow+vsep+ssm1)/4
            do j=1,ss
               if (j==1) then
                  phaml_matrix%matrix_val(a2row+a2sep+j-1) = &
                     phaml_matrix%matrix_val(a2row+a2sep+j-1) + &
                     phaml_matrix%matrix_val(va2+vsep+i-1)
               elseif (j==i) then
                  phaml_matrix%matrix_val(a2row+a2sep+j-1) = &
                     phaml_matrix%matrix_val(a2row+a2sep+j-1) + &
                     phaml_matrix%matrix_val(va2+vsep)
               else
                  phaml_matrix%matrix_val(a2row+a2sep+j-1) = &
                     phaml_matrix%matrix_val(a2row+a2sep+j-1) + &
                     phaml_matrix%matrix_val(va2+vsep+j-1)
               endif
            end do
! <a1,a2>
            phaml_matrix%matrix_val(a1a2+a1sep:a1a2+a1sep+ssm1) = &
               phaml_matrix%matrix_val(va1+vsep:va1+vsep+ssm1)/2 + &
               phaml_matrix%matrix_val(va2+vsep:va2+vsep+ssm1)/2
            do j=1,ss
               if (j==1) then
                  phaml_matrix%matrix_val(a1a2+a1sep+j-1) = &
                     phaml_matrix%matrix_val(a1a2+a1sep+j-1) + &
                     phaml_matrix%matrix_val(vrow+vsep+i-1)/4
               elseif (j==i) then
                  phaml_matrix%matrix_val(a1a2+a1sep+j-1) = &
                     phaml_matrix%matrix_val(a1a2+a1sep+j-1) + &
                     phaml_matrix%matrix_val(vrow+vsep)/4
               else
                  phaml_matrix%matrix_val(a1a2+a1sep+j-1) = &
                     phaml_matrix%matrix_val(a1a2+a1sep+j-1) + &
                     phaml_matrix%matrix_val(vrow+vsep+j-1)/4
               endif
            end do
! <a2,a1>
            phaml_matrix%matrix_val(a2a1+a2sep:a2a1+a2sep+ssm1) = &
               phaml_matrix%matrix_val(a1a2+a1sep:a1a2+a1sep+ssm1)
! <a1,p1>
            phaml_matrix%matrix_val(a1p1+a1sep:a1p1+a1sep+ssm1) = &
               phaml_matrix%matrix_val(a1p1+a1sep:a1p1+a1sep+ssm1) + &
               phaml_matrix%matrix_val(vp1+vsep:vp1+vsep+ssm1)/2
! <p1,a1>
            phaml_matrix%matrix_val(p1a1+p1sep:p1a1+p1sep+ssm1) = &
               phaml_matrix%matrix_val(a1p1+a1sep:a1p1+a1sep+ssm1)
! <a2,p1>
            phaml_matrix%matrix_val(a2p1+a2sep:a2p1+a2sep+ssm1) = &
               phaml_matrix%matrix_val(a2p1+a2sep:a2p1+a2sep+ssm1) + &
               phaml_matrix%matrix_val(vp1+vsep:vp1+vsep+ssm1)/2
! <p1,a2>
            phaml_matrix%matrix_val(p1a2+p1sep:p1a2+p1sep+ssm1) = &
               phaml_matrix%matrix_val(a2p1+a2sep:a2p1+a2sep+ssm1)
            if (.not. boundary_vert) then
! <a1,p2>
               phaml_matrix%matrix_val(a1p2+a1sep:a1p2+a1sep+ssm1) = &
                  phaml_matrix%matrix_val(a1p2+a1sep:a1p2+a1sep+ssm1) + &
                  phaml_matrix%matrix_val(vp2+vsep:vp2+vsep+ssm1)/2
! <p2,a1>
               phaml_matrix%matrix_val(p2a1+p2sep:p2a1+p2sep+ssm1) = &
                  phaml_matrix%matrix_val(a1p2+a1sep:a1p2+a1sep+ssm1)
! <a2,p2>
               phaml_matrix%matrix_val(a2p2+a2sep:a2p2+a2sep+ssm1) = &
                  phaml_matrix%matrix_val(a2p2+a2sep:a2p2+a2sep+ssm1) + &
                  phaml_matrix%matrix_val(vp2+vsep:vp2+vsep+ssm1)/2
! <p2,a2>
               phaml_matrix%matrix_val(p2a2+p2sep:p2a2+p2sep+ssm1) = &
                  phaml_matrix%matrix_val(a2p2+a2sep:a2p2+a2sep+ssm1)
            endif
! <v,a1>
            do j=1,ss
               if (j==1) then
                  phaml_matrix%matrix_val(va1+vsep+j-1) = &
                     phaml_matrix%matrix_val(va1+vsep+j-1) + &
                     phaml_matrix%matrix_val(vrow+vsep+i-1)/2
               elseif (j==i) then
                  phaml_matrix%matrix_val(va1+vsep+j-1) = &
                     phaml_matrix%matrix_val(va1+vsep+j-1) + &
                     phaml_matrix%matrix_val(vrow+vsep)/2
               else
                  phaml_matrix%matrix_val(va1+vsep+j-1) = &
                     phaml_matrix%matrix_val(va1+vsep+j-1) + &
                     phaml_matrix%matrix_val(vrow+vsep+j-1)/2
               endif
            end do
! <v,a2>
            do j=1,ss
               if (j==1) then
                  phaml_matrix%matrix_val(va2+vsep+j-1) = &
                     phaml_matrix%matrix_val(va2+vsep+j-1) + &
                     phaml_matrix%matrix_val(vrow+vsep+i-1)/2
               elseif (j==i) then
                  phaml_matrix%matrix_val(va2+vsep+j-1) = &
                     phaml_matrix%matrix_val(va2+vsep+j-1) + &
                     phaml_matrix%matrix_val(vrow+vsep)/2
               else
                  phaml_matrix%matrix_val(va2+vsep+j-1) = &
                     phaml_matrix%matrix_val(va2+vsep+j-1) + &
                     phaml_matrix%matrix_val(vrow+vsep+j-1)/2
               endif
            end do
         endif ! do_matrix
! remove red vertex from parents' column indices
         phaml_matrix%column_index(p1v+p1sep:p1v+p1sep+ssm1) = NO_ENTRY
         if (.not. boundary_vert) then
            phaml_matrix%column_index(p2v+p2sep:p2v+p2sep+ssm1) = NO_ENTRY
         endif
      end do

   else ! direction == TO_NODAL

      do i=1,ss
         vsep  = phaml_matrix%begin_row(v +i-1) - phaml_matrix%begin_row(v )
         a1sep = phaml_matrix%begin_row(a1+i-1) - phaml_matrix%begin_row(a1)
         a2sep = phaml_matrix%begin_row(a2+i-1) - phaml_matrix%begin_row(a2)
         p1sep = phaml_matrix%begin_row(p1+i-1) - phaml_matrix%begin_row(p1)
         if (.not. boundary_vert) then
          p2sep = phaml_matrix%begin_row(p2+i-1) - phaml_matrix%begin_row(p2)
         endif
! <a1,a1>
         if (do_matrix) then
            phaml_matrix%matrix_val(a1row+a1sep:a1row+a1sep+ssm1) = &
               phaml_matrix%matrix_val(a1row+a1sep:a1row+a1sep+ssm1) + &
               phaml_matrix%matrix_val(vrow+vsep:vrow+vsep+ssm1)/4
            do j=1,ss
               if (j==1) then
                  phaml_matrix%matrix_val(a1row+a1sep+j-1) = &
                     phaml_matrix%matrix_val(a1row+a1sep+j-1) - &
                     phaml_matrix%matrix_val(va1+vsep+i-1)
               elseif (j==i) then
                  phaml_matrix%matrix_val(a1row+a1sep+j-1) = &
                     phaml_matrix%matrix_val(a1row+a1sep+j-1) - &
                     phaml_matrix%matrix_val(va1+vsep)
               else
                  phaml_matrix%matrix_val(a1row+a1sep+j-1) = &
                     phaml_matrix%matrix_val(a1row+a1sep+j-1) - &
                     phaml_matrix%matrix_val(va1+vsep+j-1)
               endif
            end do
! <a2,a2>
            phaml_matrix%matrix_val(a2row+a2sep:a2row+a2sep+ssm1) = &
               phaml_matrix%matrix_val(a2row+a2sep:a2row+a2sep+ssm1) + &
               phaml_matrix%matrix_val(vrow+vsep:vrow+vsep+ssm1)/4
            do j=1,ss
               if (j==1) then
                  phaml_matrix%matrix_val(a2row+a2sep+j-1) = &
                     phaml_matrix%matrix_val(a2row+a2sep+j-1) - &
                     phaml_matrix%matrix_val(va2+vsep+i-1)
               elseif (j==i) then
                  phaml_matrix%matrix_val(a2row+a2sep+j-1) = &
                     phaml_matrix%matrix_val(a2row+a2sep+j-1) - &
                     phaml_matrix%matrix_val(va2+vsep)
               else
                  phaml_matrix%matrix_val(a2row+a2sep+j-1) = &
                     phaml_matrix%matrix_val(a2row+a2sep+j-1) - &
                     phaml_matrix%matrix_val(va2+vsep+j-1)
               endif
            end do
! <a1,p1>
            phaml_matrix%matrix_val(a1p1+a1sep:a1p1+a1sep+ssm1) = &
               phaml_matrix%matrix_val(a1p1+a1sep:a1p1+a1sep+ssm1) - &
               phaml_matrix%matrix_val(vp1+vsep:vp1+vsep+ssm1)/2
! <p1,a1>
            phaml_matrix%matrix_val(p1a1+p1sep:p1a1+p1sep+ssm1) = &
               phaml_matrix%matrix_val(a1p1+a1sep:a1p1+a1sep+ssm1)
! <a2,p1>
            phaml_matrix%matrix_val(a2p1+a2sep:a2p1+a2sep+ssm1) = &
               phaml_matrix%matrix_val(a2p1+a2sep:a2p1+a2sep+ssm1) - &
               phaml_matrix%matrix_val(vp1+vsep:vp1+vsep+ssm1)/2
! <p1,a2>
            phaml_matrix%matrix_val(p1a2+p1sep:p1a2+p1sep+ssm1) = &
               phaml_matrix%matrix_val(a2p1+a2sep:a2p1+a2sep+ssm1)
            if (.not. boundary_vert) then
! <a1,p2>
               phaml_matrix%matrix_val(a1p2+a1sep:a1p2+a1sep+ssm1) = &
                  phaml_matrix%matrix_val(a1p2+a1sep:a1p2+a1sep+ssm1) - &
                  phaml_matrix%matrix_val(vp2+vsep:vp2+vsep+ssm1)/2
! <p2,a1>
               phaml_matrix%matrix_val(p2a1+p2sep:p2a1+p2sep+ssm1) = &
                  phaml_matrix%matrix_val(a1p2+a1sep:a1p2+a1sep+ssm1)
! <a2,p2>
               phaml_matrix%matrix_val(a2p2+a2sep:a2p2+a2sep+ssm1) = &
                  phaml_matrix%matrix_val(a2p2+a2sep:a2p2+a2sep+ssm1) - &
                  phaml_matrix%matrix_val(vp2+vsep:vp2+vsep+ssm1)/2
! <p2,a2>
               phaml_matrix%matrix_val(p2a2+p2sep:p2a2+p2sep+ssm1) = &
                  phaml_matrix%matrix_val(a2p2+a2sep:a2p2+a2sep+ssm1)
            endif
! <v,a1>
            do j=1,ss
               if (j==1) then
                  phaml_matrix%matrix_val(va1+vsep+j-1) = &
                     phaml_matrix%matrix_val(va1+vsep+j-1) - &
                     phaml_matrix%matrix_val(vrow+vsep+i-1)/2
               elseif (j==i) then
                  phaml_matrix%matrix_val(va1+vsep+j-1) = &
                     phaml_matrix%matrix_val(va1+vsep+j-1) - &
                     phaml_matrix%matrix_val(vrow+vsep)/2
               else
                  phaml_matrix%matrix_val(va1+vsep+j-1) = &
                     phaml_matrix%matrix_val(va1+vsep+j-1) - &
                     phaml_matrix%matrix_val(vrow+vsep+j-1)/2
               endif
            end do
! <v,a2>
            do j=1,ss
               if (j==1) then
                  phaml_matrix%matrix_val(va2+vsep+j-1) = &
                     phaml_matrix%matrix_val(va2+vsep+j-1) - &
                     phaml_matrix%matrix_val(vrow+vsep+i-1)/2
               elseif (j==i) then
                  phaml_matrix%matrix_val(va2+vsep+j-1) = &
                     phaml_matrix%matrix_val(va2+vsep+j-1) - &
                     phaml_matrix%matrix_val(vrow+vsep)/2
               else
                  phaml_matrix%matrix_val(va2+vsep+j-1) = &
                     phaml_matrix%matrix_val(va2+vsep+j-1) - &
                     phaml_matrix%matrix_val(vrow+vsep+j-1)/2
               endif
            end do
! <a1,v>
            phaml_matrix%matrix_val(a1v+a1sep:a1v+a1sep+ssm1) = &
               phaml_matrix%matrix_val(va1+vsep:va1+vsep+ssm1)
! <a2,v>
            phaml_matrix%matrix_val(a2v+a2sep:a2v+a2sep+ssm1) = &
               phaml_matrix%matrix_val(va2+vsep:va2+vsep+ssm1)
! <p1,v>
            phaml_matrix%matrix_val(p1v+p1sep:p1v+p1sep+ssm1) = &
               phaml_matrix%matrix_val(vp1+vsep:vp1+vsep+ssm1)
! <p2,v>
            if (.not. boundary_vert) then
               phaml_matrix%matrix_val(p2v+p2sep:p2v+p2sep+ssm1) = &
                  phaml_matrix%matrix_val(vp2+vsep:vp2+vsep+ssm1)
               phaml_matrix%column_index(p2v+p2sep:p2v+p2sep+ssm1) = &
                  phaml_matrix%column_index(a1v+a1sep:a1v+a1sep+ssm1)
            endif
         endif ! do_matrix
! column index changes
! <a1,v>
         phaml_matrix%column_index(a1v+a1sep:a1v+a1sep+ssm1) = &
            (/ ((v+j),j=0,ssm1) /)
! <a2,v>
         phaml_matrix%column_index(a2v+a2sep:a2v+a2sep+ssm1) = &
            phaml_matrix%column_index(a1v+a1sep:a1v+a1sep+ssm1)
! <p1,v>
         phaml_matrix%column_index(p1v+p1sep:p1v+p1sep+ssm1) = &
            phaml_matrix%column_index(a1v+a1sep:a1v+a1sep+ssm1)
! <p2,v>
         if (.not. boundary_vert) then
            phaml_matrix%column_index(p2v+p2sep:p2v+p2sep+ssm1) = &
               phaml_matrix%column_index(a1v+a1sep:a1v+a1sep+ssm1)
         endif
      end do

   endif ! direction

! end of transforming matrix values

! transform mass matrix values

   if (do_matrix .and. loc_do_mass) then
      if (direction == TO_HIER) then

         do i=1,ss
            vsep  = phaml_matrix%begin_row(v +i-1) - phaml_matrix%begin_row(v )
            a1sep = phaml_matrix%begin_row(a1+i-1) - phaml_matrix%begin_row(a1)
            a2sep = phaml_matrix%begin_row(a2+i-1) - phaml_matrix%begin_row(a2)
            p1sep = phaml_matrix%begin_row(p1+i-1) - phaml_matrix%begin_row(p1)
            if (.not. boundary_vert) then
             p2sep = phaml_matrix%begin_row(p2+i-1) - phaml_matrix%begin_row(p2)
            endif
! <a1,a1>
            phaml_matrix%mass(a1row+a1sep:a1row+a1sep+ssm1) = &
               phaml_matrix%mass(a1row+a1sep:a1row+a1sep+ssm1) + &
               phaml_matrix%mass(vrow+vsep:vrow+vsep+ssm1)/4
            do j=1,ss
               if (j==1) then
                  phaml_matrix%mass(a1row+a1sep+j-1) = &
                     phaml_matrix%mass(a1row+a1sep+j-1) + &
                     phaml_matrix%mass(va1+vsep+i-1)
               elseif (j==i) then
                  phaml_matrix%mass(a1row+a1sep+j-1) = &
                     phaml_matrix%mass(a1row+a1sep+j-1) + &
                     phaml_matrix%mass(va1+vsep)
               else
                  phaml_matrix%mass(a1row+a1sep+j-1) = &
                     phaml_matrix%mass(a1row+a1sep+j-1) + &
                     phaml_matrix%mass(va1+vsep+j-1)
               endif
            end do
! <a2,a2>
            phaml_matrix%mass(a2row+a2sep:a2row+a2sep+ssm1) = &
               phaml_matrix%mass(a2row+a2sep:a2row+a2sep+ssm1) + &
               phaml_matrix%mass(vrow+vsep:vrow+vsep+ssm1)/4
            do j=1,ss
               if (j==1) then
                  phaml_matrix%mass(a2row+a2sep+j-1) = &
                     phaml_matrix%mass(a2row+a2sep+j-1) + &
                     phaml_matrix%mass(va2+vsep+i-1)
               elseif (j==i) then
                  phaml_matrix%mass(a2row+a2sep+j-1) = &
                     phaml_matrix%mass(a2row+a2sep+j-1) + &
                     phaml_matrix%mass(va2+vsep)
               else
                  phaml_matrix%mass(a2row+a2sep+j-1) = &
                     phaml_matrix%mass(a2row+a2sep+j-1) + &
                     phaml_matrix%mass(va2+vsep+j-1)
               endif
            end do
! <a1,a2>
            phaml_matrix%mass(a1a2+a1sep:a1a2+a1sep+ssm1) = &
               phaml_matrix%mass(va1+vsep:va1+vsep+ssm1)/2 + &
               phaml_matrix%mass(va2+vsep:va2+vsep+ssm1)/2
            do j=1,ss
               if (j==1) then
                  phaml_matrix%mass(a1a2+a1sep+j-1) = &
                     phaml_matrix%mass(a1a2+a1sep+j-1) + &
                     phaml_matrix%mass(vrow+vsep+i-1)/4
               elseif (j==i) then
                  phaml_matrix%mass(a1a2+a1sep+j-1) = &
                     phaml_matrix%mass(a1a2+a1sep+j-1) + &
                     phaml_matrix%mass(vrow+vsep)/4
               else
                  phaml_matrix%mass(a1a2+a1sep+j-1) = &
                     phaml_matrix%mass(a1a2+a1sep+j-1) + &
                     phaml_matrix%mass(vrow+vsep+j-1)/4
               endif
            end do
! <a2,a1>
            phaml_matrix%mass(a2a1+a2sep:a2a1+a2sep+ssm1) = &
               phaml_matrix%mass(a1a2+a1sep:a1a2+a1sep+ssm1)
! <a1,p1>
            phaml_matrix%mass(a1p1+a1sep:a1p1+a1sep+ssm1) = &
               phaml_matrix%mass(a1p1+a1sep:a1p1+a1sep+ssm1) + &
               phaml_matrix%mass(vp1+vsep:vp1+vsep+ssm1)/2
! <p1,a1>
            phaml_matrix%mass(p1a1+p1sep:p1a1+p1sep+ssm1) = &
               phaml_matrix%mass(a1p1+a1sep:a1p1+a1sep+ssm1)
! <a2,p1>
            phaml_matrix%mass(a2p1+a2sep:a2p1+a2sep+ssm1) = &
               phaml_matrix%mass(a2p1+a2sep:a2p1+a2sep+ssm1) + &
               phaml_matrix%mass(vp1+vsep:vp1+vsep+ssm1)/2
! <p1,a2>
            phaml_matrix%mass(p1a2+p1sep:p1a2+p1sep+ssm1) = &
               phaml_matrix%mass(a2p1+a2sep:a2p1+a2sep+ssm1)
            if (.not. boundary_vert) then
! <a1,p2>
               phaml_matrix%mass(a1p2+a1sep:a1p2+a1sep+ssm1) = &
                  phaml_matrix%mass(a1p2+a1sep:a1p2+a1sep+ssm1) + &
                  phaml_matrix%mass(vp2+vsep:vp2+vsep+ssm1)/2
! <p2,a1>
               phaml_matrix%mass(p2a1+p2sep:p2a1+p2sep+ssm1) = &
                  phaml_matrix%mass(a1p2+a1sep:a1p2+a1sep+ssm1)
! <a2,p2>
               phaml_matrix%mass(a2p2+a2sep:a2p2+a2sep+ssm1) = &
                  phaml_matrix%mass(a2p2+a2sep:a2p2+a2sep+ssm1) + &
                  phaml_matrix%mass(vp2+vsep:vp2+vsep+ssm1)/2
! <p2,a2>
               phaml_matrix%mass(p2a2+p2sep:p2a2+p2sep+ssm1) = &
                  phaml_matrix%mass(a2p2+a2sep:a2p2+a2sep+ssm1)
            endif
! <v,a1>
            do j=1,ss
               if (j==1) then
                  phaml_matrix%mass(va1+vsep+j-1) = &
                     phaml_matrix%mass(va1+vsep+j-1) + &
                     phaml_matrix%mass(vrow+vsep+i-1)/2
               elseif (j==i) then
                  phaml_matrix%mass(va1+vsep+j-1) = &
                     phaml_matrix%mass(va1+vsep+j-1) + &
                     phaml_matrix%mass(vrow+vsep)/2
               else
                  phaml_matrix%mass(va1+vsep+j-1) = &
                     phaml_matrix%mass(va1+vsep+j-1) + &
                     phaml_matrix%mass(vrow+vsep+j-1)/2
               endif
            end do
! <v,a2>
            do j=1,ss
               if (j==1) then
                  phaml_matrix%mass(va2+vsep+j-1) = &
                     phaml_matrix%mass(va2+vsep+j-1) + &
                     phaml_matrix%mass(vrow+vsep+i-1)/2
               elseif (j==i) then
                  phaml_matrix%mass(va2+vsep+j-1) = &
                     phaml_matrix%mass(va2+vsep+j-1) + &
                     phaml_matrix%mass(vrow+vsep)/2
               else
                  phaml_matrix%mass(va2+vsep+j-1) = &
                     phaml_matrix%mass(va2+vsep+j-1) + &
                     phaml_matrix%mass(vrow+vsep+j-1)/2
               endif
            end do
         end do

      else ! direction == TO_NODAL

         do i=1,ss
            vsep  = phaml_matrix%begin_row(v +i-1) - phaml_matrix%begin_row(v )
            a1sep = phaml_matrix%begin_row(a1+i-1) - phaml_matrix%begin_row(a1)
            a2sep = phaml_matrix%begin_row(a2+i-1) - phaml_matrix%begin_row(a2)
            p1sep = phaml_matrix%begin_row(p1+i-1) - phaml_matrix%begin_row(p1)
            if (.not. boundary_vert) then
             p2sep = phaml_matrix%begin_row(p2+i-1) - phaml_matrix%begin_row(p2)
            endif
! <a1,a1>
            phaml_matrix%mass(a1row+a1sep:a1row+a1sep+ssm1) = &
               phaml_matrix%mass(a1row+a1sep:a1row+a1sep+ssm1) + &
               phaml_matrix%mass(vrow+vsep:vrow+vsep+ssm1)/4
            do j=1,ss
               if (j==1) then
                  phaml_matrix%mass(a1row+a1sep+j-1) = &
                     phaml_matrix%mass(a1row+a1sep+j-1) - &
                     phaml_matrix%mass(va1+vsep+i-1)
               elseif (j==i) then
                  phaml_matrix%mass(a1row+a1sep+j-1) = &
                     phaml_matrix%mass(a1row+a1sep+j-1) - &
                     phaml_matrix%mass(va1+vsep)
               else
                  phaml_matrix%mass(a1row+a1sep+j-1) = &
                     phaml_matrix%mass(a1row+a1sep+j-1) - &
                     phaml_matrix%mass(va1+vsep+j-1)
               endif
            end do
! <a2,a2>
            phaml_matrix%mass(a2row+a2sep:a2row+a2sep+ssm1) = &
               phaml_matrix%mass(a2row+a2sep:a2row+a2sep+ssm1) + &
               phaml_matrix%mass(vrow+vsep:vrow+vsep+ssm1)/4
            do j=1,ss
               if (j==1) then
                  phaml_matrix%mass(a2row+a2sep+j-1) = &
                     phaml_matrix%mass(a2row+a2sep+j-1) - &
                     phaml_matrix%mass(va2+vsep+i-1)
               elseif (j==i) then
                  phaml_matrix%mass(a2row+a2sep+j-1) = &
                     phaml_matrix%mass(a2row+a2sep+j-1) - &
                     phaml_matrix%mass(va2+vsep)
               else
                  phaml_matrix%mass(a2row+a2sep+j-1) = &
                     phaml_matrix%mass(a2row+a2sep+j-1) - &
                     phaml_matrix%mass(va2+vsep+j-1)
               endif
            end do
! <a1,p1>
            phaml_matrix%mass(a1p1+a1sep:a1p1+a1sep+ssm1) = &
               phaml_matrix%mass(a1p1+a1sep:a1p1+a1sep+ssm1) - &
               phaml_matrix%mass(vp1+vsep:vp1+vsep+ssm1)/2
! <p1,a1>
            phaml_matrix%mass(p1a1+p1sep:p1a1+p1sep+ssm1) = &
               phaml_matrix%mass(a1p1+a1sep:a1p1+a1sep+ssm1)
! <a2,p1>
            phaml_matrix%mass(a2p1+a2sep:a2p1+a2sep+ssm1) = &
               phaml_matrix%mass(a2p1+a2sep:a2p1+a2sep+ssm1) - &
               phaml_matrix%mass(vp1+vsep:vp1+vsep+ssm1)/2
! <p1,a2>
            phaml_matrix%mass(p1a2+p1sep:p1a2+p1sep+ssm1) = &
               phaml_matrix%mass(a2p1+a2sep:a2p1+a2sep+ssm1)
            if (.not. boundary_vert) then
! <a1,p2>
               phaml_matrix%mass(a1p2+a1sep:a1p2+a1sep+ssm1) = &
                  phaml_matrix%mass(a1p2+a1sep:a1p2+a1sep+ssm1) - &
                  phaml_matrix%mass(vp2+vsep:vp2+vsep+ssm1)/2
! <p2,a1>
               phaml_matrix%mass(p2a1+p2sep:p2a1+p2sep+ssm1) = &
                  phaml_matrix%mass(a1p2+a1sep:a1p2+a1sep+ssm1)
! <a2,p2>
               phaml_matrix%mass(a2p2+a2sep:a2p2+a2sep+ssm1) = &
                  phaml_matrix%mass(a2p2+a2sep:a2p2+a2sep+ssm1) - &
                  phaml_matrix%mass(vp2+vsep:vp2+vsep+ssm1)/2
! <p2,a2>
               phaml_matrix%mass(p2a2+p2sep:p2a2+p2sep+ssm1) = &
                  phaml_matrix%mass(a2p2+a2sep:a2p2+a2sep+ssm1)
            endif
! <v,a1>
            do j=1,ss
               if (j==1) then
                  phaml_matrix%mass(va1+vsep+j-1) = &
                     phaml_matrix%mass(va1+vsep+j-1) - &
                     phaml_matrix%mass(vrow+vsep+i-1)/2
               elseif (j==i) then
                  phaml_matrix%mass(va1+vsep+j-1) = &
                     phaml_matrix%mass(va1+vsep+j-1) - &
                     phaml_matrix%mass(vrow+vsep)/2
               else
                  phaml_matrix%mass(va1+vsep+j-1) = &
                     phaml_matrix%mass(va1+vsep+j-1) - &
                     phaml_matrix%mass(vrow+vsep+j-1)/2
               endif
            end do
! <v,a2>
            do j=1,ss
               if (j==1) then
                  phaml_matrix%mass(va2+vsep+j-1) = &
                     phaml_matrix%mass(va2+vsep+j-1) - &
                     phaml_matrix%mass(vrow+vsep+i-1)/2
               elseif (j==i) then
                  phaml_matrix%mass(va2+vsep+j-1) = &
                     phaml_matrix%mass(va2+vsep+j-1) - &
                     phaml_matrix%mass(vrow+vsep)/2
               else
                  phaml_matrix%mass(va2+vsep+j-1) = &
                     phaml_matrix%mass(va2+vsep+j-1) - &
                     phaml_matrix%mass(vrow+vsep+j-1)/2
               endif
            end do
! <a1,v>
            phaml_matrix%mass(a1v+a1sep:a1v+a1sep+ssm1) = &
               phaml_matrix%mass(va1+vsep:va1+vsep+ssm1)
! <a2,v>
            phaml_matrix%mass(a2v+a2sep:a2v+a2sep+ssm1) = &
               phaml_matrix%mass(va2+vsep:va2+vsep+ssm1)
! <p1,v>
            phaml_matrix%mass(p1v+p1sep:p1v+p1sep+ssm1) = &
               phaml_matrix%mass(vp1+vsep:vp1+vsep+ssm1)
! <p2,v>
            if (.not. boundary_vert) then
               phaml_matrix%mass(p2v+p2sep:p2v+p2sep+ssm1) = &
                  phaml_matrix%mass(vp2+vsep:vp2+vsep+ssm1)
            endif
         end do

      endif ! direction
   endif ! do_mass

! end of transforming mass matrix values

end do ! next vertex

end subroutine basis_change

!          -------------------
subroutine exchange_neigh_vect(vector,procs,linear_system,tag1,tag2,tag3, &
                               max_lev,notime)
!          -------------------

!----------------------------------------------------
! This routine exchanges the given vector between processors
! at equations that are immediate neighbors, i.e., nonzero coefficient in matrix
! If max_lev is present, it only exchanges points with level up to max_lev.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(inout) :: vector(:)
type(proc_info), intent(in) :: procs
type(linsys_type), intent(in) :: linear_system
integer, intent(in) :: tag1, tag2, tag3
integer, optional :: max_lev
logical, optional :: notime

!----------------------------------------------------
! Local variables:

type(hash_key_eq) :: gid
real(my_real), allocatable :: send_real(:)
real(my_real), pointer :: recv_real(:)
integer, allocatable :: send_int(:), nsend(:), nrecv(:)
integer, pointer :: recv_int(:)
integer :: i, counter, eq, allocstat, ind, p, lid, rind, sindi, sindr, nproc, &
           hilim
logical :: timeit, owned_neigh
!----------------------------------------------------
! Begin executable code

! if there is only one processor, nothing to exchange

nproc = num_proc(procs)
if (nproc == 1) return

! time communication part of the solve

if (present(notime)) then
   timeit = .not. notime
else
   timeit = .true. 
endif
if (timeit) then
   call start_watch((/cpsolve,ctsolve/))
endif

! determine how many equations I don't own but are a neighbor of one I do own

if (present(max_lev)) then
   hilim = linear_system%begin_level(max_lev+1)-1
else
   hilim = linear_system%neq
endif

counter = 0
do eq=1,hilim
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

! allocate space for the first messages

allocate(send_int(counter),nsend(nproc),nrecv(nproc),stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in exchange_neigh_vect",procs=procs)
   return
endif

! make a list of the equations I don't own but are a neighbor of one I do own

counter = 0
do eq=1,hilim
   if (.not. linear_system%iown(eq)) then
      owned_neigh = .false.
      do i=linear_system%begin_row(eq)+1,linear_system%end_row(eq)
         if (linear_system%column_index(i) == NO_ENTRY) cycle
         if (linear_system%iown(linear_system%column_index(i))) then
            owned_neigh = .true.
            exit
         endif
      end do
      if (owned_neigh) then
         call hash_pack_key(linear_system%gid(eq),send_int,counter+1)
         counter = counter + KEY_SIZE+1
      endif
   endif 
end do

! exchange lists with other processors

call phaml_alltoall(procs,send_int,counter,recv_int,nrecv,tag1)
deallocate(send_int,stat=allocstat)

! count the number I own in each of the other processors' requests

nsend = 0
ind = 1
do p=1,nproc
   do eq=1,nrecv(p)/(KEY_SIZE+1)
      gid = hash_unpack_key(recv_int,ind,extended=.true.)
      lid = hash_decode_key(gid,linear_system%eq_hash)
      if (lid /= HASH_NOT_FOUND) then
         if (linear_system%iown(lid)) then
            nsend(p) = nsend(p) + 1
         endif
      endif
      ind = ind + KEY_SIZE+1
   end do
end do

! allocate space for next messages

allocate(send_int(sum(nsend)*(KEY_SIZE+1)),send_real(sum(nsend)), &
         stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in exchange_neigh_vect",procs=procs)
   return
endif

! create the reply with values of vector for equations I own

rind = 1
sindi = 1
sindr = 1
do p=1,nproc
   do eq=1,nrecv(p)/(KEY_SIZE+1)
      gid = hash_unpack_key(recv_int,rind,extended=.true.)
      lid = hash_decode_key(gid,linear_system%eq_hash)
      if (lid /= HASH_NOT_FOUND) then
         if (linear_system%iown(lid)) then
            send_int(sindi:sindi+KEY_SIZE) = recv_int(rind:rind+KEY_SIZE)
            send_real(sindr) = vector(lid)
            sindi = sindi + KEY_SIZE+1
            sindr = sindr + 1
         endif
      endif
      rind = rind + KEY_SIZE+1
   end do
end do
if (associated(recv_int)) deallocate(recv_int,stat=allocstat)

! send the reply

call phaml_alltoall(procs,send_int,nsend*(KEY_SIZE+1),recv_int,nrecv,tag2)
call phaml_alltoall(procs,send_real,nsend,recv_real,nrecv,tag3)

deallocate(send_int,send_real,nsend,stat=allocstat)

! copy the replies into vector

do eq=1,sum(nrecv)
   gid = hash_unpack_key(recv_int,1+(eq-1)*(KEY_SIZE+1),extended=.true.)
   lid = hash_decode_key(gid,linear_system%eq_hash)
   if (lid == HASH_NOT_FOUND) then
      call warning("received reply for an equation I don't have in exchange_neigh_vect")
   else
      vector(lid) = recv_real(eq)
   endif
end do

deallocate(recv_int,recv_real,nrecv,stat=allocstat)

if (timeit) then
   call stop_watch((/cpsolve,ctsolve/))
endif

end subroutine exchange_neigh_vect

!          -------------------
subroutine old_exchange_fudop_vect(vector,procs,linear_system,tag1,tag2,tag3, &
                               only_lev,max_lev,notime)
!          -------------------

!----------------------------------------------------
! This routine exchanges the given vector between processors
! at all FuDoP equations.  If only_lev is present, it only exchanges
! points at level only_lev.  If max_lev is present, it only exchanges
! points with level up to max_lev.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(inout) :: vector(:)
type(proc_info), intent(in) :: procs
type(linsys_type), intent(in) :: linear_system
integer, intent(in) :: tag1, tag2, tag3
integer, intent(in), optional :: only_lev, max_lev
logical, optional, intent(in) :: notime

!----------------------------------------------------
! Local variables:

type(hash_key_eq) :: gid
real(my_real), allocatable :: send_real(:)
real(my_real), pointer :: recv_real(:)
integer, allocatable :: send_int(:), nsend(:), nrecv(:)
integer, pointer :: recv_int(:)
integer :: counter, lolim, hilim, eq, allocstat, ind, p, lid, rind, sindi, &
           sindr, nproc
logical :: timeit
!----------------------------------------------------
! Begin executable code

! if there is only one processor, nothing to exchange

nproc = num_proc(procs)
if (nproc == 1) return

! time communication part of the solve

if (present(notime)) then
   timeit = .not. notime
else
   timeit = .true. 
endif
if (timeit) then
   call start_watch((/cpsolve,ctsolve/))
endif

! determine how many equations I have but don't own

counter = 0
if (present(only_lev)) then
   lolim = linear_system%begin_level(only_lev)
   hilim = linear_system%begin_level(only_lev+1)-1
elseif (present(max_lev)) then
   lolim = 1
   hilim = linear_system%begin_level(max_lev+1)-1
else
   lolim = 1
   hilim = linear_system%neq
endif
do eq=lolim,hilim
   if (.not. linear_system%iown(eq)) then
      counter = counter + KEY_SIZE+1
   endif
end do

! allocate space for the first messages

allocate(send_int(counter),nsend(nproc),nrecv(nproc),stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in exchange_fudop_vect",procs=procs)
   return
endif

! make a list of the equations I don't own

counter = 0
do eq=lolim,hilim
   if (.not. linear_system%iown(eq)) then
      call hash_pack_key(linear_system%gid(eq),send_int,counter+1)
      counter = counter + KEY_SIZE+1
   endif 
end do

! exchange lists of unowned equations with other processors

call phaml_alltoall(procs,send_int,counter,recv_int,nrecv,tag1)
deallocate(send_int,stat=allocstat)

! count the number I own in each of the other processors' requests

nsend = 0
ind = 1
do p=1,nproc
   do eq=1,nrecv(p)/(KEY_SIZE+1)
      gid = hash_unpack_key(recv_int,ind,extended=.true.)
      lid = hash_decode_key(gid,linear_system%eq_hash)
      if (lid /= HASH_NOT_FOUND) then
         if (linear_system%iown(lid)) then
            nsend(p) = nsend(p) + 1
         endif
      endif
      ind = ind + KEY_SIZE+1
   end do
end do

! allocate space for next messages

allocate(send_int(sum(nsend)*(KEY_SIZE+1)),send_real(sum(nsend)), &
         stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in exchange_fudop_vect",procs=procs)
   return
endif

! create the reply with values of vector for equations I own

rind = 1
sindi = 1
sindr = 1
do p=1,nproc
   do eq=1,nrecv(p)/(KEY_SIZE+1)
      gid = hash_unpack_key(recv_int,rind,extended=.true.)
      lid = hash_decode_key(gid,linear_system%eq_hash)
      if (lid /= HASH_NOT_FOUND) then
         if (linear_system%iown(lid)) then
            send_int(sindi:sindi+KEY_SIZE) = recv_int(rind:rind+KEY_SIZE)
            send_real(sindr) = vector(lid)
            sindi = sindi + KEY_SIZE+1
            sindr = sindr + 1
         endif
      endif
      rind = rind + KEY_SIZE+1
   end do
end do
if (associated(recv_int)) deallocate(recv_int,stat=allocstat)

! send the reply

call phaml_alltoall(procs,send_int,nsend*(KEY_SIZE+1),recv_int,nrecv,tag2)
call phaml_alltoall(procs,send_real,nsend,recv_real,nrecv,tag3)

deallocate(send_int,send_real,nsend,stat=allocstat)

! copy the replies into vector

do eq=1,sum(nrecv)
   gid = hash_unpack_key(recv_int,1+(eq-1)*(KEY_SIZE+1),extended=.true.)
   lid = hash_decode_key(gid,linear_system%eq_hash)
   if (lid == HASH_NOT_FOUND) then
      call warning("received reply for an equation I don't have in exchange_fudop_vect")
   else
      vector(lid) = recv_real(eq)
   endif
end do

deallocate(recv_int,recv_real,nrecv,stat=allocstat)

if (timeit) then
   call stop_watch((/cpsolve,ctsolve/))
endif

end subroutine old_exchange_fudop_vect

!          -------------------
subroutine exchange_fudop_vect(vector,procs,linear_system,tag1,tag2,tag3, &
                               notime)
!          -------------------

!----------------------------------------------------
! This routine exchanges the given vector between processors
! at all FuDoP equations.  If only_lev is present, it only exchanges
! points at level only_lev.  If max_lev is present, it only exchanges
! points with level up to max_lev.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(inout) :: vector(:)
type(proc_info), intent(in) :: procs
type(linsys_type), intent(in) :: linear_system
! TEMP090210 when all is said and done, I only need 2 tag
integer, intent(in) :: tag1, tag2, tag3
logical, optional, intent(in) :: notime
!----------------------------------------------------
! Local variables:

type int_arrays
   integer, pointer :: val(:)
end type int_arrays
type real_arrays
   real(my_real), pointer :: val(:)
end type real_arrays
integer :: nproc, my_processor, p, q, ni, nr, i, j, step, loop, recv_rind, &
           to, from, mess_size, proxy, nstep, ndelay
logical :: timeit
integer, allocatable :: send_iind(:), send_rind(:), send_iind_next(:), &
                        send_rind_next(:)
integer, pointer :: recv_int(:)
real(my_real), pointer :: recv_real(:)
type(int_arrays), allocatable :: send_int(:), send_int_next(:)
type(real_arrays), allocatable :: send_real(:), send_real_next(:)
!----------------------------------------------------
! Begin executable code

if (.not. new_comm) then ! TEMP090209
   call old_exchange_fudop_vect(vector,procs,linear_system,tag1,tag2,tag3,notime=notime) ! TEMP090209
   return ! TEMP090209
endif ! TEMP090209

nproc = num_proc(procs)
my_processor = my_proc(procs)

! conditions under which there is nothing to do

if (nproc==1 .or. my_processor==MASTER) return

! time communication part of the solve

if (present(notime)) then
   timeit = .not. notime
else
   timeit = .true. 
endif
if (timeit) then
   call start_watch((/cpsolve,ctsolve/))
endif

allocate(send_int(nproc),send_real(nproc),send_int_next(nproc), &
         send_real_next(nproc),send_iind(nproc), send_rind(nproc), &
         send_iind_next(nproc),send_rind_next(nproc))

! create and send the first messages

do p=1,nproc
   if (linear_system%fudop_comm_mess_size(p,1) == 0) cycle
   allocate(send_int(p)%val(3*linear_system%fudop_comm_nchunk(p,1)), &
            send_real(p)%val(linear_system%fudop_comm_mess_size(p,1)))
   send_iind(p) = 0
   send_rind(p) = 0
   do q=1,nproc
      if (linear_system%fudop_comm_nsend(q) /= 0) then
         if (linear_system%fudop_comm_proxy(q) == p) then
            send_int(p)%val(send_iind(p)+1) = q
            send_int(p)%val(send_iind(p)+2) = my_processor
            send_int(p)%val(send_iind(p)+3) = linear_system%fudop_comm_nsend(q)
            send_iind(p) = send_iind(p) + 3
            do i=1,linear_system%fudop_comm_nsend(q)
               send_real(p)%val(send_rind(p)+i) = vector(linear_system%fudop_comm_send_lid(q)%lid(i))
            end do
            send_rind(p) = send_rind(p) + linear_system%fudop_comm_nsend(q)
         endif
      endif
   end do
end do

! multiple steps of sending messages to nearest neighbors to get to distant
! locations

nstep=size(linear_system%fudop_comm_nproc_recv)
do step = 1,nstep

! send the next messages, except if the message is large and the receiver
! has a smaller ID than mine then delay the message to avoid deadlock

   do p=1,nproc
      if (linear_system%fudop_comm_mess_size(p,step) /= 0) then
         if (p > my_processor .or. &
             linear_system%fudop_comm_mess_size(p,step) <= max_real_buffer) then
            call phaml_send(procs,p,send_int(p)%val,send_iind(p), &
                            send_real(p)%val,send_rind(p),tag1+13*step)
            deallocate(send_int(p)%val,send_real(p)%val)
         else
            call phaml_send(procs,p,(/0,0,-1/),3,(/0.0_my_real/),0, &
                            tag1+13*step)
         endif
      endif
   end do

! allocate space for the next message

   if (step /= nstep) then
      do p=1,nproc
         if (linear_system%fudop_comm_mess_size(p,step+1) /= 0) then
            allocate(send_int_next(p)%val(3*linear_system%fudop_comm_nchunk(p,step+1)),&
                    send_real_next(p)%val(linear_system%fudop_comm_mess_size(p,step+1)))
         endif
      end do
   endif

   send_iind_next = 0
   send_rind_next = 0
   ndelay = 0

! for the number of messages to be received ...

   do loop=1,linear_system%fudop_comm_nproc_recv(step)

! receive a message

      call phaml_recv(procs,p,recv_int,ni,recv_real,nr,tag1+13*step)

! if recv_int(3), which normally contains the message size of the first chunk,
! is negative, this message has been delayed

      if (recv_int(3) < 0) then
         ndelay = ndelay + 1
         cycle
      endif

! chunks for which I am the destination get processed; the others are passed on

      recv_rind = 0
      do i=1,ni,3
         to = recv_int(i)
         from = recv_int(i+1)
         mess_size = recv_int(i+2)
         if (to == my_processor) then
            do j=1,mess_size
               vector(linear_system%fudop_comm_recv_lid(from)%lid(j)) = &
                  recv_real(recv_rind+j)
            end do
         else
            proxy = linear_system%fudop_comm_proxy(to)
            send_int_next(proxy)%val(send_iind_next(proxy)+1) = to
            send_int_next(proxy)%val(send_iind_next(proxy)+2) = from
            send_int_next(proxy)%val(send_iind_next(proxy)+3) = mess_size
            send_iind_next(proxy) = send_iind_next(proxy) + 3
            send_real_next(proxy)%val(send_rind_next(proxy)+1:send_rind_next(proxy)+mess_size) = &
               recv_real(recv_rind+1:recv_rind+mess_size)
            send_rind_next(proxy) = send_rind_next(proxy) + mess_size
         end if
         recv_rind = recv_rind + mess_size
      end do
      if (associated(recv_int)) deallocate(recv_int)
      if (associated(recv_real)) deallocate(recv_real)
   end do

! send the messages that were delayed

   do p=1,my_processor-1
      if (linear_system%fudop_comm_mess_size(p,step) > max_real_buffer) then
         call phaml_send(procs,p,send_int(p)%val,send_iind(p), &
                         send_real(p)%val,send_rind(p),tag2+13*step)
         deallocate(send_int(p)%val,send_real(p)%val)
      endif
   end do

! receive and unpack or forward the delayed messages

   do loop=1,ndelay
      call phaml_recv(procs,p,recv_int,ni,recv_real,nr,tag2+13*step)
      recv_rind = 0
      do i=1,ni,3
         to = recv_int(i)
         from = recv_int(i+1)
         mess_size = recv_int(i+2)
         if (to == my_processor) then
            do j=1,mess_size
               vector(linear_system%fudop_comm_recv_lid(from)%lid(j)) = &
                  recv_real(recv_rind+j)
            end do
         else
            proxy = linear_system%fudop_comm_proxy(to)
            send_int_next(proxy)%val(send_iind_next(proxy)+1) = to
            send_int_next(proxy)%val(send_iind_next(proxy)+2) = from
            send_int_next(proxy)%val(send_iind_next(proxy)+3) = mess_size
            send_iind_next(proxy) = send_iind_next(proxy) + 3
            send_real_next(proxy)%val(send_rind_next(proxy)+1:send_rind_next(proxy)+mess_size) = &
               recv_real(recv_rind+1:recv_rind+mess_size)
            send_rind_next(proxy) = send_rind_next(proxy) + mess_size
         end if
         recv_rind = recv_rind + mess_size
      end do
      if (associated(recv_int)) deallocate(recv_int)
      if (associated(recv_real)) deallocate(recv_real)
   end do

! move the next messages to the current messages for the next loop

   if (step /= nstep) then
      do p=1,nproc
         if (linear_system%fudop_comm_mess_size(p,step+1) /= 0) then
            send_int(p)%val => send_int_next(p)%val
            send_real(p)%val => send_real_next(p)%val
            nullify(send_int_next(p)%val,send_real_next(p)%val)
            send_iind(p) = send_iind_next(p)
            send_rind(p) = send_rind_next(p)
         endif
      end do
   end if

end do ! step

deallocate(send_iind,send_rind,send_int,send_real,send_iind_next, &
           send_rind_next,send_int_next,send_real_next)

if (timeit) then
   call stop_watch((/cpsolve,ctsolve/))
endif

end subroutine exchange_fudop_vect

!          ----------------------------
subroutine o_exchange_fudop_soln_residual(procs,linear_system,tag1,tag2,tag3, &
                                        no_soln,max_lev)
!          ----------------------------

!----------------------------------------------------
! This routine exchanges the solution and residual between processors
! at all FuDoP equations
! If no_soln is present and true, only the residual is passed.
! If max_lev is present, only entities of level max_lev or less are passed.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(proc_info), intent(in) :: procs
type(linsys_type), intent(inout) :: linear_system
integer, intent(in) :: tag1, tag2, tag3
logical, intent(in), optional :: no_soln
integer, intent(in), optional :: max_lev

!----------------------------------------------------
! Local variables:

type(hash_key_eq) :: gid
real(my_real), allocatable :: send_real(:)
real(my_real), pointer :: recv_real(:)
integer, allocatable :: send_int(:), nsend(:), nrecv(:)
integer, pointer :: recv_int(:)
integer :: i, counter, eq, allocstat, ind, p, lid, rind, sindi, sindr, nproc, &
           my_processor, hilim
logical :: do_resid, do_soln
!----------------------------------------------------
! Begin executable code

! if there is only one processor, nothing to exchange

nproc = num_proc(procs)
if (nproc == 1) return

! time communication part of the solve

call start_watch((/cpsolve,ctsolve/))

my_processor = my_proc(procs)
linear_system%r_others = 0.0_my_real

if (present(no_soln)) then
   do_soln = .not. no_soln
else
   do_soln = .true.
endif

! determine the size of the first message

if (present(max_lev)) then
   hilim = linear_system%begin_level(max_lev+1)-1
else
   hilim = linear_system%neq
endif

counter = 0
do eq=1,hilim
   if (linear_system%need_r_others(eq)) then
      counter = counter + KEY_SIZE+1
   endif
   if (do_soln .and. .not. linear_system%iown(eq)) then
      counter = counter + KEY_SIZE+1
   endif
end do
counter = counter + KEY_SIZE+1

! allocate space for the first messages

allocate(send_int(counter),nsend(nproc),nrecv(nproc),stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in exchange_fudop_soln_residual",procs=procs)
   return
endif

! make a list of all the equations where I need r_others, and a list
! of all equations I don't own for which I need the solution

counter = 0
do eq=1,hilim
   if (linear_system%need_r_others(eq)) then
      call hash_pack_key(linear_system%gid(eq),send_int,counter+1)
      counter = counter + KEY_SIZE+1
   endif
end do

call hash_pack_key(NULL_KEY_EQ,send_int,counter+1)
counter = counter + KEY_SIZE+1

if (do_soln) then
   do eq=1,hilim
      if (.not. linear_system%iown(eq)) then
         call hash_pack_key(linear_system%gid(eq),send_int,counter+1)
         counter = counter + KEY_SIZE+1
      endif 
   end do
endif

! exchange lists with other processors

call phaml_alltoall(procs,send_int,counter,recv_int,nrecv,tag1)
deallocate(send_int,stat=allocstat)

! from other partition's requests, count the residuals I have and solutions
! I own

nsend = 0
ind = 1
do p=1,nproc
   do_resid = .true.
   do eq=1,nrecv(p)/(KEY_SIZE+1)
      gid = hash_unpack_key(recv_int,ind,extended=.true.)
      if (gid == NULL_KEY_EQ) then
         do_resid = .false.
         nsend(p) = nsend(p) + 1
         ind = ind + KEY_SIZE+1
         cycle
      endif
      lid = hash_decode_key(gid,linear_system%eq_hash)
      if (lid /= HASH_NOT_FOUND) then
         if (do_resid .or. linear_system%iown(lid)) then
            nsend(p) = nsend(p) + 1
         endif
      endif
      ind = ind + KEY_SIZE+1
   end do
end do

! allocate space for next messages

allocate(send_int(sum(nsend)*(KEY_SIZE+1)),send_real(sum(nsend)), &
         stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in exchange_fudop_soln_residual",procs=procs)
   return
endif

! create the reply with the residuals and solutions

rind = 1
sindi = 1
sindr = 1
do p=1,nproc
   do_resid = .true.
   do eq=1,nrecv(p)/(KEY_SIZE+1)
      gid = hash_unpack_key(recv_int,rind,extended=.true.)
      if (gid == NULL_KEY_EQ) then
         do_resid = .false.
         call hash_pack_key(NULL_KEY_EQ,send_int,sindi)
         sindi = sindi + KEY_SIZE+1
         send_real(sindr) = 0.0_my_real
         sindr = sindr + 1
         rind = rind + KEY_SIZE+1
         cycle
      endif
      lid = hash_decode_key(gid,linear_system%eq_hash)
      if (lid /= HASH_NOT_FOUND) then
         if (do_resid .or. linear_system%iown(lid)) then
            send_int(sindi:sindi+KEY_SIZE) = recv_int(rind:rind+KEY_SIZE)
            if (do_resid) then
               send_real(sindr) = linear_system%r_mine(lid)
            else
               send_real(sindr) = linear_system%solution(lid)
            endif
            sindi = sindi + KEY_SIZE+1
            sindr = sindr + 1
         endif
      endif
      rind = rind + KEY_SIZE+1
   end do
end do
if (associated(recv_int)) deallocate(recv_int,stat=allocstat)

! send the reply

call phaml_alltoall(procs,send_int,nsend*(KEY_SIZE+1),recv_int,nrecv,tag2)
call phaml_alltoall(procs,send_real,nsend,recv_real,nrecv,tag3)

deallocate(send_int,send_real,nsend,stat=allocstat)

! add the replies into r_others or assign to solution

eq = 0
do p=1,nproc
   do_resid = .true.
   do i=1,nrecv(p)
      eq = eq + 1
      gid = hash_unpack_key(recv_int,1+(eq-1)*(KEY_SIZE+1),extended=.true.)
      if (gid == NULL_KEY_EQ) then
         do_resid = .false.
         cycle
      endif
      lid = hash_decode_key(gid,linear_system%eq_hash)
      if (lid == HASH_NOT_FOUND) then
         call warning("received reply for an equation I don't have in exchange_fudop_soln_residual")
      else
         if (p /= my_processor) then
            if (do_resid) then
               linear_system%r_others(lid) = linear_system%r_others(lid) + recv_real(eq)
            else
               linear_system%solution(lid) = recv_real(eq)
            endif
         endif
      endif
   end do
end do

deallocate(recv_int,recv_real,nrecv,stat=allocstat)

call stop_watch((/cpsolve,ctsolve/))

end subroutine o_exchange_fudop_soln_residual

!          ----------------------------
subroutine exchange_fudop_soln_residual(procs,linear_system,tag1,tag2,tag3, &
                                        no_soln,max_lev,notime)
!          ----------------------------

!----------------------------------------------------
! This routine exchanges the given vector between processors
! at all FuDoP equations.  If only_lev is present, it only exchanges
! points at level only_lev.  If max_lev is present, it only exchanges
! points with level up to max_lev.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(proc_info), intent(in) :: procs
type(linsys_type), intent(inout) :: linear_system
! TEMP090210 when all is said and done, I only need 2 tags
integer, intent(in) :: tag1, tag2, tag3
logical, intent(in), optional :: no_soln
integer, intent(in), optional :: max_lev
logical, intent(in), optional :: notime
!----------------------------------------------------
! Local variables:

type int_arrays
   integer, pointer :: val(:)
end type int_arrays
type real_arrays
   real(my_real), pointer :: val(:)
end type real_arrays
integer :: nproc, my_processor, p, q, ni, nr, i, j, step, loop, recv_rind, &
           to, from, mess_size, proxy, ind, lasteq, nstep, ndelay
logical :: timeit, do_soln
integer, allocatable :: send_iind(:), send_rind(:), send_iind_next(:), &
                        send_rind_next(:)
integer, pointer :: recv_int(:)
real(my_real), pointer :: recv_real(:)
type(int_arrays), allocatable :: send_int(:), send_int_next(:)
type(real_arrays), allocatable :: send_real(:), send_real_next(:)
!----------------------------------------------------
! Begin executable code

if (.not. new_comm) then ! TEMP090209
   call o_exchange_fudop_soln_residual(procs,linear_system,tag1,tag2,tag3, &
                                        no_soln,max_lev) ! TEMP090209
   return ! TEMP090209
endif ! TEMP090209

nproc = num_proc(procs)
my_processor = my_proc(procs)

! conditions under which there is nothing to do

if (nproc==1 .or. my_processor==MASTER) return

! time communication part of the solve

if (present(notime)) then
   timeit = .not. notime
else
   timeit = .true. 
endif
if (timeit) then
   call start_watch((/cpsolve,ctsolve/))
endif

if (present(no_soln)) then
   do_soln = .not. no_soln
else
   do_soln = .true.
endif

if (present(max_lev)) then
   lasteq = linear_system%begin_level(max_lev+1) - 1
else
   lasteq = linear_system%neq
endif

linear_system%r_others = 0.0_my_real
allocate(send_int(nproc),send_real(nproc),send_int_next(nproc), &
         send_real_next(nproc),send_iind(nproc), send_rind(nproc), &
         send_iind_next(nproc),send_rind_next(nproc))


! create the first messages

do p=1,nproc
   if (linear_system%resid_comm_mess_size(p,1) == 0) cycle
   allocate(send_int(p)%val(3*linear_system%resid_comm_nchunk(p,1)), &
            send_real(p)%val(linear_system%resid_comm_mess_size(p,1)))
   send_iind(p) = 0
   send_rind(p) = 0
   do q=1,nproc
      if (linear_system%resid_comm_nsend(q) /= 0) then
         if (linear_system%fudop_comm_proxy(q) == p) then
            send_int(p)%val(send_iind(p)+1) = q
            send_int(p)%val(send_iind(p)+2) = my_processor
            send_int(p)%val(send_iind(p)+3) = linear_system%fudop_comm_nsend(q) + &
                                              linear_system%resid_comm_nsend(q)
            send_iind(p) = send_iind(p) + 3
            do i=1,linear_system%fudop_comm_nsend(q)
               send_real(p)%val(send_rind(p)+i) = linear_system%solution(linear_system%fudop_comm_send_lid(q)%lid(i))
            end do
            send_rind(p) = send_rind(p) + linear_system%fudop_comm_nsend(q)
            do i=1,linear_system%resid_comm_nsend(q)
               send_real(p)%val(send_rind(p)+i) = linear_system%r_mine(linear_system%resid_comm_send_lid(q)%lid(i))
            end do
            send_rind(p) = send_rind(p) + linear_system%resid_comm_nsend(q)
         endif
      endif
   end do
end do

! multiple steps of sending messages to nearest neighbors to get to distant
! locations

nstep=size(linear_system%resid_comm_nproc_recv)
do step = 1,nstep

! send the next messages, except if the message is large and the receiver
! has a smaller ID than mine then delay the message to avoid deadlock

   do p=1,nproc
      if (linear_system%resid_comm_mess_size(p,step) /= 0) then
         if (p > my_processor .or. &
             linear_system%resid_comm_mess_size(p,step) <= max_real_buffer) then
            call phaml_send(procs,p,send_int(p)%val,send_iind(p), &
                            send_real(p)%val,send_rind(p),tag1+13*step)
            deallocate(send_int(p)%val,send_real(p)%val)
         else
            call phaml_send(procs,p,(/0,0,-1/),3,(/0.0_my_real/),0, &
                            tag1+13*step)
         endif
      endif
   end do

! allocate space for the next messages

   if (step /= nstep) then
      do p=1,nproc
         if (linear_system%resid_comm_mess_size(p,step+1) /= 0) then
            allocate(send_int_next(p)%val(3*linear_system%resid_comm_nchunk(p,step+1)),&
                     send_real_next(p)%val(linear_system%resid_comm_mess_size(p,step+1)))
         endif
      end do
   endif

   send_iind_next = 0
   send_rind_next = 0
   ndelay = 0

! for the number of messages to be received ...

   do loop=1,linear_system%resid_comm_nproc_recv(step)

! receive a message

      call phaml_recv(procs,p,recv_int,ni,recv_real,nr,tag1+13*step)

! if recv_int(3), which normally contains the message size of the first chunk,
! is negative, this message has been delayed

      if (recv_int(3) < 0) then
         ndelay = ndelay + 1
         cycle
      endif

! chunks for which I am the destination get processed; the others are passed on

      recv_rind = 0
      do i=1,ni,3
         to = recv_int(i)
         from = recv_int(i+1)
         mess_size = recv_int(i+2)
         if (to == my_processor) then
            do j=1,mess_size
               if (j <= size(linear_system%fudop_comm_recv_lid(from)%lid)) then
                  if (do_soln) then
                     ind = linear_system%fudop_comm_recv_lid(from)%lid(j)
                     if (ind <= lasteq) then
                        linear_system%solution(ind) = recv_real(recv_rind+j)
                     endif
                  endif
               else
                  ind = j - size(linear_system%fudop_comm_recv_lid(from)%lid)
                  ind = linear_system%resid_comm_recv_lid(from)%lid(ind)
                  if (ind <= lasteq) then
                     linear_system%r_others(ind) = linear_system%r_others(ind)+&
                                                   recv_real(recv_rind+j)
                  endif
               endif
            end do
         else
            proxy = linear_system%fudop_comm_proxy(to)
            send_int_next(proxy)%val(send_iind_next(proxy)+1) = to
            send_int_next(proxy)%val(send_iind_next(proxy)+2) = from
            send_int_next(proxy)%val(send_iind_next(proxy)+3) = mess_size
            send_iind_next(proxy) = send_iind_next(proxy) + 3
            send_real_next(proxy)%val(send_rind_next(proxy)+1:send_rind_next(proxy)+mess_size) = &
               recv_real(recv_rind+1:recv_rind+mess_size)
            send_rind_next(proxy) = send_rind_next(proxy) + mess_size
         end if
         recv_rind = recv_rind + mess_size
      end do
      if (associated(recv_int)) deallocate(recv_int)
      if (associated(recv_real)) deallocate(recv_real)
   end do

! send the messages that were delayed

   do p=1,my_processor-1
      if (linear_system%resid_comm_mess_size(p,step) > max_real_buffer) then
         call phaml_send(procs,p,send_int(p)%val,send_iind(p), &
                         send_real(p)%val,send_rind(p),tag2+13*step)
         deallocate(send_int(p)%val,send_real(p)%val)
      endif
   end do

! receive and unpack or forward the delayed messages

   do loop=1,ndelay
      call phaml_recv(procs,p,recv_int,ni,recv_real,nr,tag2+13*step)
      recv_rind = 0
      do i=1,ni,3
         to = recv_int(i)
         from = recv_int(i+1)
         mess_size = recv_int(i+2)
         if (to == my_processor) then
            do j=1,mess_size
               if (j <= size(linear_system%fudop_comm_recv_lid(from)%lid)) then
                  if (do_soln) then
                     ind = linear_system%fudop_comm_recv_lid(from)%lid(j)
                     if (ind <= lasteq) then
                        linear_system%solution(ind) = recv_real(recv_rind+j)
                     endif
                  endif
               else
                  ind = j - size(linear_system%fudop_comm_recv_lid(from)%lid)
                  ind = linear_system%resid_comm_recv_lid(from)%lid(ind)
                  if (ind <= lasteq) then
                     linear_system%r_others(ind) = linear_system%r_others(ind)+&
                                                   recv_real(recv_rind+j)
                  endif
               endif
            end do
         else
            proxy = linear_system%fudop_comm_proxy(to)
            send_int_next(proxy)%val(send_iind_next(proxy)+1) = to
            send_int_next(proxy)%val(send_iind_next(proxy)+2) = from
            send_int_next(proxy)%val(send_iind_next(proxy)+3) = mess_size
            send_iind_next(proxy) = send_iind_next(proxy) + 3
            send_real_next(proxy)%val(send_rind_next(proxy)+1:send_rind_next(proxy)+mess_size) = &
               recv_real(recv_rind+1:recv_rind+mess_size)
            send_rind_next(proxy) = send_rind_next(proxy) + mess_size
         end if
         recv_rind = recv_rind + mess_size
      end do
      if (associated(recv_int)) deallocate(recv_int)
      if (associated(recv_real)) deallocate(recv_real)
   end do

! move the next messages to the current messages for the next loop

   if (step /= nstep) then
      do p=1,nproc
         if (linear_system%resid_comm_mess_size(p,step+1) /= 0) then
            send_int(p)%val => send_int_next(p)%val
            send_real(p)%val => send_real_next(p)%val
            nullify(send_int_next(p)%val,send_real_next(p)%val)
            send_iind(p) = send_iind_next(p)
            send_rind(p) = send_rind_next(p)
         endif
      end do
   end if

end do ! step

deallocate(send_iind,send_rind,send_int,send_real,send_iind_next, &
           send_rind_next,send_int_next,send_real_next)

if (timeit) then
   call stop_watch((/cpsolve,ctsolve/))
endif

end subroutine exchange_fudop_soln_residual

!          --------------
! TEMP090204
subroutine old_sum_fudop_vect(vector,procs,linear_system,tag1,tag2,tag3, &
                          only_lev)
!          --------------

!----------------------------------------------------
! This routine sums the given vector over all processors
! at all FuDoP equations.  If only_lev is present, it only sums
! points at level only_lev
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(inout) :: vector(:)
type(proc_info), intent(in) :: procs
type(linsys_type), intent(in) :: linear_system
integer, intent(in) :: tag1, tag2, tag3
integer, intent(in), optional :: only_lev 

!----------------------------------------------------
! Local variables:

type(hash_key_eq) :: gid
real(my_real), allocatable :: send_real(:)
real(my_real), pointer :: recv_real(:)
integer, allocatable :: send_int(:), nsend(:), nrecv(:)
integer, pointer :: recv_int(:)
integer :: counter, lolim, hilim, eq, allocstat, ind, p, lid, rind, sindi, &
           sindr, nproc
!----------------------------------------------------
! Begin executable code

! if there is only one processor, nothing to exchange

nproc = num_proc(procs)
if (nproc == 1) return

! time communication part of the solve

call start_watch((/cpsolve,ctsolve/))

! determine how many equations I have

if (present(only_lev)) then
   lolim = linear_system%begin_level(only_lev)
   hilim = linear_system%begin_level(only_lev+1)-1
else
   lolim = 1
   hilim = linear_system%neq
endif
counter = (hilim - lolim + 1)*(KEY_SIZE+1)

! allocate space for the first messages

allocate(send_int(counter),nsend(nproc),nrecv(nproc),stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in sum_fudop_vect",procs=procs)
   return
endif

! make a list of the equations I have

counter = 0
do eq=lolim,hilim
   call hash_pack_key(linear_system%gid(eq),send_int,counter+1)
   counter = counter + KEY_SIZE+1
end do

! exchange lists of unowned equations with other processors

call phaml_alltoall(procs,send_int,counter,recv_int,nrecv,tag1)
deallocate(send_int,stat=allocstat)

! count the number I have in each of the other processors' requests

nsend = 0
ind = 1
do p=1,nproc
   do eq=1,nrecv(p)/(KEY_SIZE+1)
      gid = hash_unpack_key(recv_int,ind,extended=.true.)
      lid = hash_decode_key(gid,linear_system%eq_hash)
      if (lid /= HASH_NOT_FOUND) then
         nsend(p) = nsend(p) + 1
      endif
      ind = ind + KEY_SIZE+1
   end do
end do

! allocate space for next messages

allocate(send_int(sum(nsend)*(KEY_SIZE+1)),send_real(sum(nsend)), &
         stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in sum_fudop_vect",procs=procs)
   return
endif

! create the reply with values of vector for equations I have

rind = 1
sindi = 1
sindr = 1
do p=1,nproc
   do eq=1,nrecv(p)/(KEY_SIZE+1)
      gid = hash_unpack_key(recv_int,rind,extended=.true.)
      lid = hash_decode_key(gid,linear_system%eq_hash)
      if (lid /= HASH_NOT_FOUND) then
         send_int(sindi:sindi+KEY_SIZE) = recv_int(rind:rind+KEY_SIZE)
         send_real(sindr) = vector(lid)
         sindi = sindi + KEY_SIZE+1
         sindr = sindr + 1
      endif
      rind = rind + KEY_SIZE+1
   end do
end do
if (associated(recv_int)) deallocate(recv_int,stat=allocstat)

! send the reply

call phaml_alltoall(procs,send_int,nsend*(KEY_SIZE+1),recv_int,nrecv,tag2)
call phaml_alltoall(procs,send_real,nsend,recv_real,nrecv,tag3)

deallocate(send_int,send_real,nsend,stat=allocstat)

! sum the replies in vector

vector = 0.0_my_real
do eq=1,sum(nrecv)
   gid = hash_unpack_key(recv_int,1+(eq-1)*(KEY_SIZE+1),extended=.true.)
   lid = hash_decode_key(gid,linear_system%eq_hash)
   if (lid == HASH_NOT_FOUND) then
      call warning("received reply for an equation I don't have in sum_fudop_vect")
   else
      vector(lid) = vector(lid) + recv_real(eq)
   endif
end do

deallocate(recv_int,recv_real,nrecv,stat=allocstat)

call stop_watch((/cpsolve,ctsolve/))

! TEMP090204
end subroutine old_sum_fudop_vect

! TEMP090204 new version
!          --------------
subroutine sum_fudop_vect(vector,procs,linear_system,tag1,tag2,tag3, &
                          only_lev)
!          --------------

!----------------------------------------------------
! This routine sums the given vector over all processors
! at all FuDoP equations.  If only_lev is present, it only sums
! points at level only_lev
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(inout) :: vector(:)
type(proc_info), intent(in) :: procs
type(linsys_type), intent(in) :: linear_system
! TEMP090210 when all is said and done, I only need 2 tags
integer, intent(in) :: tag1, tag2, tag3
integer, intent(in), optional :: only_lev

!----------------------------------------------------
! Local variables:

type(hash_key_eq) :: gid
real(my_real), allocatable :: send_real(:)
real(my_real), pointer :: recv_real(:)
integer, allocatable :: send_int(:), nrecv(:)
integer, pointer :: recv_int(:)
integer :: counter, counter2, lolim, hilim, eq, allocstat, ind, lid, nproc
!----------------------------------------------------
! Begin executable code

! TEMP090204
if (.not. new_comm) then
   call old_sum_fudop_vect(vector,procs,linear_system,tag1,tag2,tag3,only_lev)
   return
endif
! end TEMP090204

! if there is only one processor, nothing to exchange

nproc = num_proc(procs)
if (nproc == 1) return

! time communication part of the solve

call start_watch((/cpsolve,ctsolve/))

! determine how many equations I have

if (present(only_lev)) then
   lolim = linear_system%begin_level(only_lev)
   hilim = linear_system%begin_level(only_lev+1)-1
else
   lolim = 1
   hilim = linear_system%neq
endif
counter = (hilim - lolim + 1)*(KEY_SIZE+1)
counter2 = hilim - lolim + 1

! allocate space for the messages

allocate(send_int(counter),send_real(counter2),nrecv(nproc), &
         stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in sum_fudop_vect",procs=procs)
   return
endif

! make a list of the equations I have and the values of vector for them

counter = 0
counter2 = 0
do eq=lolim,hilim
   call hash_pack_key(linear_system%gid(eq),send_int,counter+1)
   counter = counter + KEY_SIZE+1
   counter2 = counter2 + 1
   send_real(counter2) = vector(eq)
end do

! exchange lists of unowned equations with other processors

call phaml_alltoall(procs,send_int,counter,recv_int,nrecv,tag1)
call phaml_alltoall(procs,send_real,counter2,recv_real,nrecv,tag2)
deallocate(send_int,send_real,stat=allocstat)

! sum the values received from each processor for equations that I have

vector(lolim:hilim) = 0
ind = 1
do eq=1,sum(nrecv)
   gid = hash_unpack_key(recv_int,ind,extended=.true.)
   lid = hash_decode_key(gid,linear_system%eq_hash)
   if (lid /= HASH_NOT_FOUND) then
      vector(lid) = vector(lid) + recv_real(eq)
   endif
   ind = ind + KEY_SIZE+1
end do

deallocate(recv_int,recv_real,nrecv,stat=allocstat)

call stop_watch((/cpsolve,ctsolve/))

end subroutine sum_fudop_vect

!          ---------------
subroutine send_neigh_vect(vector,procs,linear_system,tag,max_deg,notime)
!          ---------------

!----------------------------------------------------
! This routine sends nearest neighbor entries of vector to nearest neighbor
! processors.  The vector is assumed to be associated with rows, or IDs, of
! the matrix, and only those associated with bases of degree less than or
! equal to max_deg are sent.  If max_deg = linear_system%maxdeg+1 then
! equations associated with face bases are included; there is no capability of
! including only including face bases up to a given degree.
! The messages should be received by recv_neigh_vect.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(inout) :: vector(:)
type(proc_info), intent(in) :: procs
type(linsys_type), intent(in) :: linear_system
integer, intent(in) :: tag, max_deg
logical, intent(in), optional :: notime
!----------------------------------------------------
! Local variables:

integer :: nproc, p, i, nr
real(my_real), allocatable :: send_real(:)
logical :: timeit
!----------------------------------------------------
! Begin executable code

nproc = num_proc(procs)

! conditions under which there is nothing to do

if (nproc==1 .or. my_proc(procs)==MASTER) return

! time communication part of the solve

if (present(notime)) then
   timeit = .not. notime
else
   timeit = .true. 
endif
if (timeit) then
   call start_watch((/cpsolve,ctsolve/))
endif

! for each processor that is a nearest neighbor

do p=1,nproc
   if (linear_system%nn_comm_end_of_send(p,linear_system%maxdeg+1) == 0) cycle

! create the message, send it, and destroy it

   nr = linear_system%nn_comm_end_of_send(p,max_deg)
   allocate(send_real(nr))
   do i=1,nr
      send_real(i) = vector(linear_system%nn_comm_send_lid(p)%lid(i))
   end do
   call phaml_send(procs,p,(/0/),0,send_real,nr,tag)
   deallocate(send_real)
end do

if (timeit) then
   call stop_watch((/cpsolve,ctsolve/))
endif
 
end subroutine send_neigh_vect

!          ---------------
subroutine recv_neigh_vect(vector,procs,linear_system,tag,notime)
!          ---------------

!----------------------------------------------------
! This routine receives nearest neighbor entries of vector from nearest neighbor
! processors.  The vector is assumed to be associated with rows, or IDs, of
! the matrix.  The messages should have been sent by send_neigh_vect.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(inout) :: vector(:)
type(proc_info), intent(in) :: procs
type(linsys_type), intent(in) :: linear_system
integer, intent(in) :: tag
logical, intent(in), optional :: notime
!----------------------------------------------------
! Local variables:

integer :: nproc, p, i, nr, proc, ni
integer, pointer :: recv_int(:)
real(my_real), pointer :: recv_real(:)
logical :: timeit
!----------------------------------------------------
! Begin executable code

nproc = num_proc(procs)

! conditions under which there is nothing to do

if (nproc==1 .or. my_proc(procs)==MASTER) return

! time communication part of the solve

if (present(notime)) then
   timeit = .not. notime
else
   timeit = .true. 
endif
if (timeit) then
   call start_watch((/cpsolve,ctsolve/))
endif

! for each processor that is a nearest neighbor

do p=1,nproc
   if (size(linear_system%nn_comm_recv_lid(p)%lid) == 0) cycle

! receive a message (not necessarily from p), copy it and destroy it

   call phaml_recv(procs,proc,recv_int,ni,recv_real,nr,tag)
   do i=1,nr
      vector(linear_system%nn_comm_recv_lid(proc)%lid(i)) = recv_real(i)
   end do
   if (associated(recv_real)) deallocate(recv_real)
end do
 
if (timeit) then
   call stop_watch((/cpsolve,ctsolve/))
endif
 
end subroutine recv_neigh_vect

!          ------------------
subroutine make_need_r_others(linear_system)
!          ------------------

!----------------------------------------------------
! This routine identifies which equations need to receive r_others, i.e.,
! rows of the matrix that have a column not owned by this processor
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(inout) :: linear_system
!----------------------------------------------------
! Local variables:

integer :: lev, i, j, eq, nblack, neigh
logical :: has_face, has_edge
integer, allocatable :: black_list(:)
logical(small_logical), allocatable :: on_black_list(:)
!----------------------------------------------------
! Begin executable code

! if the matrix includes face bases, change to ending with edge bases

has_face = associated(linear_system%end_row,linear_system%end_row_face)
if (has_face) call basis_change(linear_system%nlev+2,TO_HIER,linear_system, &
                                skip_matrix=.true.)

! if the matrix includes edge bases, pass through all vertex and edge rows
! to identify any rows that need r_others in relax_ho, and change to ending
! with linear bases

has_edge = associated(linear_system%end_row,linear_system%end_row_edge)
if (has_edge .and. linear_system%maxdeg > 1) then
   do eq = linear_system%begin_level(1),linear_system%begin_level(linear_system%nlev+2)-1
      do i=linear_system%begin_row(eq)+1,linear_system%end_row(eq)
         if (linear_system%column_index(i) /= NO_ENTRY) then
            if (.not. linear_system%iown(linear_system%column_index(i))) then
               linear_system%need_r_others(eq) = .true.
            endif
         endif
      end do
   end do
endif
if (has_edge) call basis_change(linear_system%nlev+1,TO_HIER,linear_system, &
                                skip_matrix=.true.)

! simulate the first half of a V cycle to identify the equations that need
! r_others at each level in relax.  See hbmg.f90/relax for comments.

allocate(on_black_list(linear_system%neq),black_list(linear_system%neq))
do lev=linear_system%nlev,2,-1
   nblack = 0
   on_black_list = .false.
   do eq = linear_system%begin_level(lev), linear_system%begin_level(lev+1)-1
      do i=linear_system%begin_row(eq)+1,linear_system%end_row(eq)
         neigh = linear_system%column_index(i)
         if (neigh == NO_ENTRY) cycle
         if (.not. linear_system%iown(neigh)) then
            linear_system%need_r_others(eq) = .true.
         endif
         if (on_black_list(neigh)) cycle
         if (neigh >= linear_system%begin_level(lev)) cycle
         on_black_list(neigh) = .true.
         nblack = nblack + 1
         black_list(nblack) = neigh
      end do
   end do
   do i = 1,nblack
      eq = black_list(i)
      do j=linear_system%begin_row(eq)+1,linear_system%end_row(eq)
         if (linear_system%column_index(j) /= NO_ENTRY) then
            if (.not. linear_system%iown(linear_system%column_index(j))) then
               linear_system%need_r_others(eq) = .true.
            endif
         endif
      end do
   end do
   call basis_change(lev,TO_HIER,linear_system,skip_matrix=.true.)
end do
deallocate(on_black_list,black_list)

! need r_others for all level 1 vertices

do eq = linear_system%begin_level(1), linear_system%begin_level(2)-1
   linear_system%need_r_others(eq) = .true.
end do

! convert back to the nodal basis

do lev=2,linear_system%nlev
   call basis_change(lev,TO_NODAL,linear_system,skip_matrix=.true.)
end do
if (has_edge) call basis_change(linear_system%nlev+1,TO_NODAL,linear_system, &
                                skip_matrix=.true.)
if (has_face) call basis_change(linear_system%nlev+2,TO_NODAL,linear_system, &
                                skip_matrix=.true.)

end subroutine make_need_r_others

!          -----------
subroutine all_to_mine(matrix,invec,outvec)
!          -----------

!----------------------------------------------------
! This routine extracts the owned entries of invec, putting them in outvec
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(in) :: matrix
real(my_real), intent(in) :: invec(:)
real(my_real), intent(out) :: outvec(:)
!----------------------------------------------------
! Local variables:

integer :: i,j
!----------------------------------------------------
! Begin executable code

j = 0
do i=1,matrix%neq
   if (matrix%iown(i)) then
      j = j+1
      outvec(j) = invec(i)
   endif
end do
outvec(j+1:) = 0.0_my_real

end subroutine all_to_mine

!          -----------
subroutine mine_to_all(matrix,invec,outvec,fill)
!          -----------

!----------------------------------------------------
! This routine spreads the entries of invec into the owned entries of
! outvec. If fill is present, the rest of outvec is filled with fill,
! otherwise it is left alone.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(in) :: matrix
real(my_real), intent(in) :: invec(:)
real(my_real), intent(inout) :: outvec(:)
real(my_real), intent(in), optional :: fill
!----------------------------------------------------
! Local variables:

integer :: i,j
!----------------------------------------------------
! Begin executable code

j = 0
do i=1,matrix%neq
   if (matrix%iown(i)) then
      j = j+1
      outvec(i) = invec(j)
   else
      if (present(fill)) outvec(i) = fill
   endif
end do

end subroutine mine_to_all

!          -------------------
subroutine matrix_times_vector(x,y,mat,procs,still_sequential,tag1,tag2, &
                               tag3,tag4,tag5,tag6,nodirch,nocomm1,nocomm2, &
                               comm2nn,natural,notime)
!          -------------------

!----------------------------------------------------
! This routine multiplies the matrix mat times the vector x and returns the
! result in the vector y.  The result is the _global_ matrix times x, not
! this processor's matrix.  Neighboring elements of x are communicated before
! the product, and the result for unowned equations are communicated after.
! If nodirch is present and true, entries in the matrix corresponding to
! Dirichlet points are not included.
! If nocomm1 is present and true, x is not communicated from the neighboring
! points.  Use this if you know the immediate neighbor shadows are current.
! If nocomm2 is present and true, y is not communicated from the fudop shadows.
! Use this if you only need the product at owned points.
! If comm2nn is present and true, y is only communicated at neighbor shadows.
! If natural is present and true, Dirichlet boundary rows are treated like
! natural boundary conditions, i.e. the actual row of the matrix is used for
! the product instead of the identity row of Dirichlet boundary conditions.
! If end_row indicates rows have been shortened to exclude face or edge bases,
! the face and/or edge bases rows are omitted in the product.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(inout), target :: x(:)
real(my_real), intent(out), target :: y(:)
type(linsys_type), intent(in) :: mat
type(proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential
! TEMP090210 when all is said and done, I only need 2 tags
integer, intent(in) :: tag1, tag2, tag3, tag4, tag5, tag6
logical, optional, intent(in) :: nodirch, nocomm1, nocomm2, comm2nn, natural, &
                                 notime

!----------------------------------------------------
! Local variables:

integer :: i, j, lasteq, loop, lastcomm
logical :: loc_nodirch, comm1, comm2, nncomm2, loc_natural

!----------------------------------------------------
! Begin executable code

if (present(nodirch)) then
   loc_nodirch = nodirch
else
   loc_nodirch = .false.
endif

if (present(nocomm1)) then
   comm1 = .not. nocomm1
else
   comm1 = .true.
endif

if (present(nocomm2)) then
   comm2 = .not. nocomm2
else
   comm2 = .true.
endif

if (present(comm2nn)) then
   nncomm2 = comm2nn
else
   nncomm2 = .false.
endif

if (present(natural)) then
   loc_natural = natural
else
   loc_natural = .false.
endif

if (associated(mat%end_row,mat%end_row_linear)) then
   lasteq = mat%begin_level(mat%nlev+1)-1
   lastcomm = 1
elseif (associated(mat%end_row,mat%end_row_edge)) then
   lasteq = mat%begin_level(mat%nlev+2)-1
   lastcomm = mat%maxdeg
else
   lasteq = mat%neq
   lastcomm = mat%maxdeg+1
endif

! On first loop, send neighbor data and do rows that don't need neighbor data.
! On second loop, receive neighbor data and do remaining rows.

y = 0.0_my_real
do loop=1,2
   if (.not. new_comm .and. loop==2) cycle ! TEMP090128

! exchange neighbor values

   if (.not. still_sequential .and. comm1) then
    if (new_comm) then ! TEMP090128
      if (loop==1) then
         call send_neigh_vect(x,procs,mat,tag1,lastcomm,notime=notime)
      else
         call recv_neigh_vect(x,procs,mat,tag1,notime=notime)
      endif
    else ! TEMP090128
      call exchange_neigh_vect(x,procs,mat,tag1,tag2,tag3,notime=notime) ! TEMP090128
    endif ! TEMP090128
   endif

! compute product for rows I own

   do i=1,lasteq
      if (.not. mat%iown(i)) cycle
      if (new_comm) then ! TEMP090127
      if ((loop==1 .and. mat%nn_comm_remote_neigh(i)) .or. &
          (loop==2 .and. .not. mat%nn_comm_remote_neigh(i))) then
         cycle
      endif
      endif ! TEMP090127
      if (mat%equation_type(i) == DIRICHLET .and. .not.loc_natural) then
         y(i) = x(i)
      else
         y(i) = 0.0_my_real
         do j=mat%begin_row(i),mat%end_row(i)
            if (mat%column_index(j) == NO_ENTRY) cycle
            if (loc_nodirch.and.mat%equation_type(mat%column_index(j))==DIRICHLET) cycle
            y(i) = y(i) + mat%matrix_val(j)*x(mat%column_index(j))
         end do
      endif
   end do
end do

! exchange result values for shadow copies

if (.not. new_comm) comm2 = .true. ! TEMP090203
if (.not. still_sequential .and. comm2) then
   if (nncomm2) then
     if (new_comm) then ! TEMP090202
      call send_neigh_vect(y,procs,mat,tag4,lastcomm,notime=notime)
      call recv_neigh_vect(y,procs,mat,tag4,notime=notime)
     else ! TEMP090202
      call exchange_fudop_vect(y,procs,mat,tag4,tag5,tag6,notime=notime) ! TEMP090202
     endif ! TEMP090202
   else
      call exchange_fudop_vect(y,procs,mat,tag4,tag5,tag6,notime=notime)
   endif
endif

return
end subroutine matrix_times_vector

!          ---------------
subroutine linsys_residual(linear_system,procs,still_sequential,iter,tag, &
                           printit,start_anew,resid,relresid,no_master)
!          ---------------

!----------------------------------------------------
! This routine computes and prints and/or returns the l2 norm of the residual
! of the linear system solution, and reduction factor.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(inout) :: linear_system
type(proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential
integer, intent(in) :: iter, tag
logical, intent(in) :: printit,start_anew
real(my_real), optional :: resid, relresid
logical, intent(in), optional :: no_master
!----------------------------------------------------
! Local variables:

real(my_real), save :: last_time = 1.0_my_real
integer :: i, proc, ni, nr
real(my_real) :: l2_norm, rhs_norm, temp(linear_system%neq), loc_resid
integer, pointer :: irecv(:)
real(my_real), pointer :: rrecv(:)
logical :: yes_master
!----------------------------------------------------
! Begin executable code

if (present(no_master)) then
   yes_master = .not. no_master
else
   yes_master = .true.
endif

! master

if (my_proc(procs) == MASTER) then

! sum up the partial sums received from the slaves

   l2_norm = 0.0_my_real
   rhs_norm = 0.0_my_real
   do i=1,num_proc(procs)
      call phaml_recv(procs,proc,irecv,ni,rrecv,nr,tag)
      l2_norm = l2_norm + rrecv(1)
      rhs_norm = rhs_norm + rrecv(2)
      deallocate(rrecv)
      if (still_sequential) exit
   end do
   l2_norm = sqrt(l2_norm)
   rhs_norm = sqrt(rhs_norm)
   if (rhs_norm == 0.0_my_real) rhs_norm = 1.0_my_real
   if (last_time == 0.0_my_real) last_time = 1.0_my_real

! print the residual

   if (printit) then
      if (start_anew) then
         write(outunit,"(A)")
         write(outunit,"(A)") "L2 norm of linear system residual"
         write(outunit,"(A)")
         write(outunit,"(1x,A)") "     iteration    absolute residual    relative residual      reduction"
         write(outunit,"(A)")
         write(outunit,"(SS,1P,1X,I11,2E18.10E2)") 0,l2_norm,l2_norm/rhs_norm
      else
         write(outunit,"(SS,1P,1X,I11,3E18.10E2)") iter,l2_norm,l2_norm/rhs_norm,l2_norm/last_time
      endif
      last_time = l2_norm
   endif

! return the residual

   if (present(resid)) resid = l2_norm
   if (present(relresid)) relresid = l2_norm/rhs_norm

! slave

else

! compute this processor's part of the residual

   call matrix_times_vector(linear_system%solution(1:),temp,linear_system, &
                            procs,still_sequential,10,20,30,0,0,0, &
                            nocomm2=.true.)
   temp = linear_system%rhs(1:linear_system%neq) + &
          linear_system%r_mine(1:linear_system%neq) + &
          linear_system%r_others(1:linear_system%neq) - temp

! compute square of l2 norm of this processor's part of residual and rhs

   l2_norm = 0.0_my_real
   rhs_norm = 0.0_my_real
   do i=1,linear_system%neq
      if (linear_system%equation_type(i) /= DIRICHLET .and. &
          linear_system%iown(i)) then
         l2_norm = l2_norm + temp(i)**2
         rhs_norm = rhs_norm + linear_system%rhs(i)**2
      endif
   end do

! send the partial norms to the master

   if ((.not. still_sequential .or. my_proc(procs) == 1) .and. yes_master) then
      call phaml_send(procs,MASTER,(/0/),0,(/l2_norm,rhs_norm/),2,tag)
   endif

! return the residual

   if (present(resid) .or. present(relresid)) then
      if (still_sequential) then
         loc_resid = sqrt(l2_norm)
      else
         loc_resid = sqrt(phaml_global_sum(procs,l2_norm,tag+5))
      endif
      if (present(resid)) then
         resid = loc_resid
      endif
      if (present(relresid)) then
         if (still_sequential) then
            rhs_norm = sqrt(rhs_norm)
         else
            rhs_norm = sqrt(phaml_global_sum(procs,rhs_norm,tag+6))
         endif
         if (rhs_norm == 0.0_my_real) rhs_norm = 1.0_my_real
         relresid = loc_resid/rhs_norm
      endif
   endif

endif

end subroutine linsys_residual

!          -----------
subroutine test_precon(invec,outvec,matrix,grid,procs,solver_cntl, &
                       still_sequential)
!          -----------

!----------------------------------------------------
! This routine provides a place to implement and test a new preconditioner
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: invec(:)
real(my_real), intent(out) :: outvec(:)
type(linsys_type), intent(inout) :: matrix
type(grid_type), intent(inout) :: grid
type(proc_info), intent(in) :: procs
type(solver_options), intent(in) :: solver_cntl
logical, intent(in) :: still_sequential
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

call warning("no test preconditioner currently implemented; using identity")

outvec = invec

end subroutine test_precon

end module linsys_util
