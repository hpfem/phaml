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

module grid_io

!----------------------------------------------------
! This module contains printed and graphical i/o for grids
!
! communication tags in this module are of the form 7xx
!----------------------------------------------------

use global
use message_passing
use hash_mod
use gridtype_mod
use stopwatch
use zoltan_interf
use error_estimators
!----------------------------------------------------

implicit none
private
public draw_grid, print_element, print_vertex, print_edge, &
       print_all_elements, print_all_vertices, print_all_edges, &
       print_entire_grid, print_grid_info

contains

!          ------------------
subroutine print_all_elements(grid)
!          ------------------

!----------------------------------------------------
! This subroutine prints the data structure for all the elements
!----------------------------------------------------
 
!----------------------------------------------------
! Dummy arguments
 
type (grid_type), intent(in) :: grid
!----------------------------------------------------
! Local variables

integer :: lev, elem

!----------------------------------------------------
! Begin executable code

do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      call print_element(grid,elem)
      elem = grid%element(elem)%next
   end do
end do

end subroutine print_all_elements

!          ---------------
subroutine print_all_edges(grid)
!          ---------------

!----------------------------------------------------
! This subroutine prints the data structure for all the edges
!----------------------------------------------------
 
!----------------------------------------------------
! Dummy arguments
 
type (grid_type), intent(in) :: grid
!----------------------------------------------------
! Local variables

integer :: lev, elem, edge, i
logical :: printed(size(grid%edge))

!----------------------------------------------------
! Begin executable code

! go through elements and print any edges that haven't already been printed

printed = .false.
do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      do i=1,EDGES_PER_ELEMENT
         edge = grid%element(elem)%edge(i)
         if (.not. printed(edge)) then
            call print_edge(grid,edge)
            printed(edge) = .true.
         endif
      end do
      elem = grid%element(elem)%next
   end do
end do

end subroutine print_all_edges

!          ------------------
subroutine print_all_vertices(grid)
!          ------------------

!----------------------------------------------------
! This subroutine prints the data structure for all the vertices
!----------------------------------------------------
 
!----------------------------------------------------
! Dummy arguments
 
type (grid_type), intent(in) :: grid
!----------------------------------------------------
! Local variables

integer lev, vert
!----------------------------------------------------
! Begin executable code

do lev=1,grid%nlev
   vert = grid%head_level_vert(lev)
   do while (vert /= END_OF_LIST)
      call print_vertex(grid,vert)
      vert = grid%vertex(vert)%next
   end do
end do

end subroutine print_all_vertices

!          -----------------
subroutine print_entire_grid(grid)
!          -----------------

!----------------------------------------------------
! This subroutine prints the data structures for the whole grid
!----------------------------------------------------
 
!----------------------------------------------------
! Dummy arguments
 
type (grid_type), intent(in) :: grid
!----------------------------------------------------
! Begin executable code

call print_all_vertices(grid)
call print_all_edges(grid)
call print_all_elements(grid)

end subroutine print_entire_grid

!          -------------
subroutine print_element(grid,elem)
!          -------------

!----------------------------------------------------
! This subroutine prints the data structure for element elem
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type (grid_type), intent(in) :: grid
integer, intent(in) :: elem

!----------------------------------------------------
! Begin executable code

write(outunit,"(A)") '--------------------'
write(outunit,"(A,I11)") 'Element number ',elem
write(outunit,"(A)") '--------------------'
call hash_print_key(grid%element(elem)%gid,outunit,"Global ID: ")
write(outunit,"(A,3I11)") 'Vertices:  ',grid%element(elem)%vertex
write(outunit,"(A,3I11)") 'Edges:     ',grid%element(elem)%edge
call hash_print_key(grid%element(elem)%mate,outunit,"Mate GID: ")
write(outunit,"(A,I11)") 'Level:     ',grid%element(elem)%level
write(outunit,"(A,I11)") 'Degree:    ',grid%element(elem)%degree
if (associated(grid%element(elem)%solution)) then
   write(outunit,"(SS,1P,A,100E18.10E2)") 'Solution:  ',grid%element(elem)%solution
else
   write(outunit,"(A)") 'Solution:  disassociated'
endif
if (associated(grid%element(elem)%exact)) then
   write(outunit,"(SS,1P,A,100E18.10E2)") 'Exact:     ',grid%element(elem)%exact
else
   write(outunit,"(A)") 'Exact:     disassociated'
endif
if (associated(grid%element(elem)%oldsoln)) then
   write(outunit,"(SS,1P,A,100E18.10E2)") 'Oldsoln:   ',grid%element(elem)%oldsoln
else
   write(outunit,"(A)") 'Oldsoln:   disassociated'
endif
write(outunit,"(A,L1)") 'Owner:     ',grid%element(elem)%iown
write(outunit,"(A,2I11)") 'Next/Prev: ',grid%element(elem)%next,grid%element(elem)%previous
write(outunit,"(A,2I11)") 'In/Out:    ',grid%element(elem)%in,grid%element(elem)%out
if (grid%element(elem)%level == 1) then
   write(outunit,"(A,3I11)") 'Neighbors: ',grid%initial_neighbor(:,elem)
endif

end subroutine print_element

!          ----------
subroutine print_edge(grid,edge)
!          ----------

!----------------------------------------------------
! This subroutine prints the data structure for edge edge
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type (grid_type), intent(in) :: grid
integer, intent(in) :: edge

!----------------------------------------------------
! Begin executable code

write(outunit,"(A)") '------------------'
write(outunit,"(A,I11)") 'Edge number ',edge
write(outunit,"(A)") '------------------'
call hash_print_key(grid%edge(edge)%gid,outunit,"Global ID: ")
write(outunit,"(A,2I11)") 'Vertices:           ',grid%edge(edge)%vertex
if (associated(grid%edge(edge)%solution)) then
   write(outunit,"(SS,1P,A,100E18.10E2)") 'Solution:           ',grid%edge(edge)%solution
else
   write(outunit,"(A)") 'Solution:           disassociated'
endif
if (associated(grid%edge(edge)%exact)) then
   write(outunit,"(SS,1P,A,100E18.10E2)") 'Exact:              ',grid%edge(edge)%exact
else
   write(outunit,"(A)") 'Exact:              disassociated'
endif
if (associated(grid%edge(edge)%oldsoln)) then
   write(outunit,"(SS,1P,A,100E18.10E2)") 'Oldsoln:            ',grid%edge(edge)%oldsoln
else
   write(outunit,"(A)") 'Oldsoln:            disassociated'
endif
write(outunit,"(A,I11)") 'Bmark:              ',grid%edge(edge)%bmark
write(outunit,"(A,100I11)") 'Type:               ',grid%edge_type(edge,:)
write(outunit,"(A,I11)") 'Degree:             ',grid%edge(edge)%degree
write(outunit,"(A,I11)") 'Associated element: ',grid%edge(edge)%assoc_elem
write(outunit,"(A,I11)") 'Next:               ',grid%edge(edge)%next

return
end subroutine print_edge

!          ------------
subroutine print_vertex(grid,vert)
!          ------------

!----------------------------------------------------
! This subroutine prints the data structure for vertex vertex
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type (grid_type), intent(in) :: grid
integer, intent(in) :: vert

!----------------------------------------------------
! Begin executable code

write(outunit,"(A)") '--------------------'
write(outunit,"(A,I11)") 'Vertex number ',vert
write(outunit,"(A)") '--------------------'
call hash_print_key(grid%vertex(vert)%gid,outunit,"Global ID: ")
write(outunit,"(SS,1P,A,2E18.10E2)") 'Coordinates:        ',grid%vertex(vert)%coord
write(outunit,"(SS,1P,A,100E18.10E2)") 'Solution:           ',grid%vertex_solution(vert,:,:)
if (associated(grid%vertex_exact)) then
   write(outunit,"(SS,1P,A,100E18.10E2)") 'Exact:              ',grid%vertex_exact(vert,:,:)
else
   write(outunit,"(A)") 'Exact:              disassociated'
endif
if (associated(grid%vertex_oldsoln)) then
   write(outunit,"(SS,1P,A,100E18.10E2)") 'Oldsoln:            ',grid%vertex_oldsoln(vert,:,:)
else
   write(outunit,"(A)") 'Oldsoln:            disassociated'
endif
write(outunit,"(A,I11)") 'Bmark:              ',grid%vertex(vert)%bmark
write(outunit,"(A,100I11)") 'Type:               ',grid%vertex_type(vert,:)
write(outunit,"(A,I11)") 'Associated element: ',grid%vertex(vert)%assoc_elem
write(outunit,"(A,2I11)") 'Next/Previous:      ',grid%vertex(vert)%next,grid%vertex(vert)%previous

return
end subroutine print_vertex

!          ---------------
subroutine print_grid_info(grid,procs,io_control,still_sequential,this_time, &
                           tag)
!          ---------------

!----------------------------------------------------
! This routine prints information about the grid
!----------------------------------------------------
 
implicit none
 
!----------------------------------------------------
! Dummy arguments
 
type (grid_type), intent(in) :: grid
type (proc_info), intent(in) :: procs
type (io_options), intent(in) :: io_control
logical, intent(in) :: still_sequential
integer, intent(in) :: this_time(:),tag
!----------------------------------------------------
 
!----------------------------------------------------
! Local variables:
 
integer :: i,proc,np,who,when,astat
integer :: nelem,nlev,nvert,nvert_own,nelem_own,nelem_leaf,nelem_leaf_own, &
           mindeg,maxdeg,dof,dof_own
integer, allocatable :: ivert(:),ielem(:),ilev(:),ielem_leaf(:), &
                        ivert_own(:),ielem_own(:),ielem_leaf_own(:), &
                        imindeg(:),imaxdeg(:),idof(:),idof_own(:)
integer :: send_int(11),ni,nr
real (my_real) :: no_reals(1)
integer, pointer :: recv_int(:)
real (my_real), pointer :: recv_real(:)
 
!----------------------------------------------------
 
!----------------------------------------------------
! Begin executable code
 
! If this is not the right time to print, return

who = io_control%print_grid_who
when = io_control%print_grid_when

if (.not. any(this_time == when) .or. who == NO_ONE) return

! If I'm the master and only slaves print the grid, return

if (my_proc(procs) == MASTER .and. who == SLAVES) return

! stop the clocks

call pause_watch(all_watches)

! If the master will print, wait for the master to get here so the slaves
! don't run away from it

if (who == MASTER .or. who == MASTER_ALL .or. who == EVERYONE) then
   if (my_proc(procs) == MASTER) then
      do proc=1,num_proc(procs)
         call phaml_send(procs,proc,(/1/),1,no_reals,0,711)
      end do
   else
      call phaml_recv(procs,proc,recv_int,ni,recv_real,nr,711)
      if (associated(recv_int)) deallocate(recv_int)
   endif
endif

! If I'm not the master, collect my grid info

if (my_proc(procs) /= MASTER) then
   call get_grid_info(grid,procs,still_sequential,710,nelem,nlev,nvert, &
                      nvert_own,nelem_own,nelem_leaf,nelem_leaf_own, &
                      dof,dof_own,mindeg=mindeg,maxdeg=maxdeg,no_master=.true.)
endif

! If the master will print, get it the info

if (who == MASTER .or. who == MASTER_ALL .or. who == EVERYONE) then
   if (my_proc(procs) /= MASTER) then
      send_int(1) = nvert
      send_int(2) = nelem
      send_int(3) = nlev
      send_int(4) = nelem_leaf
      send_int(5) = nvert_own
      send_int(6) = nelem_own
      send_int(7) = nelem_leaf_own
      send_int(8) = mindeg
      send_int(9) = maxdeg
      send_int(10) = dof
      send_int(11) = dof_own
      ni = 11
      nr = 0
      call phaml_send(procs,MASTER,send_int,ni,no_reals,nr,tag)
   else
      np = num_proc(procs)
      allocate(ivert(np),ielem(np),ilev(np),ielem_leaf(np),ivert_own(np), &
               ielem_own(np),ielem_leaf_own(np),imindeg(np),imaxdeg(np), &
               idof(np),idof_own(np),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in print_grid_info",procs=procs)
         return
      endif
      do i=1,np
         call phaml_recv(procs,proc,recv_int,ni,recv_real,nr,tag)
         ivert(proc) = recv_int(1)
         ielem(proc) = recv_int(2)
         ilev(proc) = recv_int(3)
         ielem_leaf(proc) = recv_int(4)
         ivert_own(proc) = recv_int(5)
         ielem_own(proc) = recv_int(6)
         ielem_leaf_own(proc) = recv_int(7)
         imindeg(proc) = recv_int(8)
         imaxdeg(proc) = recv_int(9)
         idof(proc) = recv_int(10)
         idof_own(proc) = recv_int(11)
         deallocate(recv_int,stat=astat)
      end do
      if (still_sequential) then
         nvert = ivert_own(1)
         nelem = ielem_own(1)
         nlev = ilev(1)
         nelem_leaf = ielem_leaf_own(1)
         mindeg = imindeg(1)
         maxdeg = imaxdeg(1)
         dof = idof(1)
      else
         nvert = sum(ivert_own)
         nelem = sum(ielem_own)
         nlev = maxval(ilev)
         nelem_leaf = sum(ielem_leaf_own)
         mindeg = minval(imindeg)
         maxdeg = maxval(imaxdeg)
         dof = sum(idof_own)
      endif
   endif
endif

! print the info, if requested

if (my_proc(procs) /= MASTER) then
 if (who == SLAVES .or. who == EVERYONE) then
   write(outunit,"(A)")
   write(outunit,"(A)") 'My Grid:'
   write(outunit,"(A,I11)") '   number of vertices       = ',nvert
   write(outunit,"(A,I11)") '   number of vertices I own = ',nvert_own
   write(outunit,"(A,I11)") '   number of elements       = ',nelem
   write(outunit,"(A,I11)") '   number of leaf elements  = ',nelem_leaf
   write(outunit,"(A,I11)") '   leaf elements I own      = ',nelem_leaf_own
   write(outunit,"(A,I11)") '   number of levels         = ',nlev
   write(outunit,"(A,I11)") '   degrees of freedom       = ',dof
   write(outunit,"(A,I11)") '   degrees of freedom I own = ',dof_own
   write(outunit,"(A,2I11)") '   min and max degree       = ',mindeg,maxdeg
 endif
else
 if (who == MASTER_ALL) then
   write(outunit,"(A)")
   write(outunit,"(A)") 'Individual Grids:'
   write(outunit,"(A,128I11)") '   number of vertices       = ',ivert
   write(outunit,"(A,128I11)") '   number of vertices I own = ',ivert_own
   write(outunit,"(A,128I11)") '   number of elements       = ',ielem
   write(outunit,"(A,128I11)") '   number of leaf elements  = ',ielem_leaf
   write(outunit,"(A,128I11)") '   leaf elements I own      = ',ielem_leaf_own
   write(outunit,"(A,128I11)") '   number of levels         = ',ilev
   write(outunit,"(A,128I11)") '   degrees of freedom       = ',idof
   write(outunit,"(A,128I11)") '   degrees of freedom I own = ',idof_own
   write(outunit,"(A,128I11)") '   min degree               = ',imindeg
   write(outunit,"(A,128I11)") '   max degree               = ',imaxdeg
 endif
 if (who == MASTER .or. who == EVERYONE .or. who == MASTER_ALL) then
   write(outunit,"(A)")
   write(outunit,"(A)") 'Total Grid:'
   write(outunit,"(A,I11)") '   number of vertices      = ',nvert
   write(outunit,"(A,I11)") '   number of leaf elements = ',nelem_leaf
   write(outunit,"(A,I11)") '   number of levels        = ',nlev
   write(outunit,"(A,I11)") '   degrees of freedom      = ',dof
   write(outunit,"(A,2I11)") '   min and max degree      = ',mindeg,maxdeg
 endif
endif

if (allocated(ivert)) then
   deallocate(ivert,ielem,ilev,ielem_leaf,ivert_own,ielem_own,ielem_leaf_own, &
              imindeg,imaxdeg,idof,idof_own,stat=astat)
endif

call end_pause_watch(all_watches)

return
end subroutine print_grid_info

!          ---------
subroutine draw_grid(grid,procs,io_control,ref_control,i_draw,master_draws, &
                     still_sequential,this_time,partition_method,lb)
!          ---------

!----------------------------------------------------
! This subroutine sends a message to the graphics process to redraw the grid.
! If the grid data has changed, it also sends the new data.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type (grid_type), intent(inout) :: grid
type (proc_info), intent(in) :: procs
type (io_options), intent(in) :: io_control
type (refine_options), intent(in) :: ref_control
logical, intent(in) :: i_draw, master_draws, still_sequential
integer, intent(in) :: this_time(:)
integer, intent(in) :: partition_method
type(Zoltan_Struct), pointer :: lb
!----------------------------------------------------

!----------------------------------------------------
! Local variables:

integer :: ni, nr, astat
integer, pointer :: imess(:)
real (my_real), pointer :: rmess(:)
integer :: imess1(1)
real (my_real) :: noreals(1)
integer, allocatable :: hold_next(:), hold_prev(:)
integer :: hold_head
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

! see if I should do anything

if (.not. (i_draw .or. master_draws)) return
if (.not. any(this_time == io_control%draw_grid_when)) return

call pause_watch(all_watches)

! slave

if (my_proc(procs) /= MASTER) then

! if the grid has changed since last sent to the graphics process ...

   if (grid_changed) then

! if using ZOLTAN_REFTREE,  preserve the level-by-level element linked list
! and get the child order from Zoltan

      if (partition_method == ZOLTAN_REFTREE) then

         allocate(hold_next(size(grid%element)), &
                  hold_prev(size(grid%element)),stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("allocation failed in draw_grid",procs=procs)
            return
         endif
         hold_next = grid%element%next
         hold_prev = grid%element%previous
         hold_head = grid%head_level_elem(1)

         call child_order_from_zoltan(grid,procs,lb,partition_method, &
                                      still_sequential)
      endif

! if slaves are plotting, send the grid to the graphics process

      if (i_draw) then
         allocate(imess(3+scalar_imess_size(grid) + &
                        elem_imess_size(grid)*grid%nelem + &
                        edge_imess_size(grid)*grid%nedge + &
                        vert_imess_size(grid)*grid%nvert), &
                  rmess(scalar_rmess_size(grid) + &
                        elem_rmess_size(grid)*grid%nelem + &
                        edge_rmess_size(grid)*grid%nedge + &
                        vert_rmess_size(grid)*grid%nvert), &
                  stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("allocation failed in draw_grid",procs=procs)
            return
         endif
         call pack_grid(grid,procs,ref_control,.false.,still_sequential, &
                        imess,rmess,ni,nr)
         if (PARALLEL == SEQUENTIAL) then
            call sequential_send(imess(1:ni),ni,rmess(1:nr),nr)
         else
            call phaml_send(procs,graphics_proc(procs),imess,ni,rmess,nr,101)
         endif
         deallocate(imess,rmess,stat=astat)
      endif

! if master is plotting, send the grid to the master

      if (master_draws) then
         allocate(imess(3+scalar_imess_size(grid) + &
                        elem_imess_size(grid)*grid%nelem + &
                        edge_imess_size(grid)*grid%nedge + &
                        vert_imess_size(grid)*grid%nvert), &
                  rmess(scalar_rmess_size(grid) + &
                        elem_rmess_size(grid)*grid%nelem + &
                        edge_rmess_size(grid)*grid%nedge + &
                        vert_rmess_size(grid)*grid%nvert), &
                  stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("allocation failed in draw_grid",procs=procs)
            return
         endif
         call pack_grid(grid,procs,ref_control,.true.,still_sequential, &
                        imess,rmess,ni,nr)
         if (.not.still_sequential .or. my_proc(procs) == 1) then
            call phaml_send(procs,MASTER,imess,ni,rmess,nr,720)
         endif
         deallocate(imess,rmess,stat=astat)
      endif

! if using ZOLTAN_REFTREE, restore the level-by-level element linked list

      if (partition_method == ZOLTAN_REFTREE) then
         grid%element%next = hold_next
         grid%element%previous = hold_prev
         grid%head_level_elem(1) = hold_head
         deallocate(hold_next, hold_prev, stat=astat)
      endif

! if the grid has not changed, just send a message to redraw

   else

      if (i_draw) then
         imess1(1) = GRAPHICS_GRID
         if (PARALLEL == SEQUENTIAL) then
            call sequential_send(imess1,1,noreals,0)
         else
            call phaml_send(procs,graphics_proc(procs),imess1,1,noreals,0,101)
         endif
      endif

   endif

! master

else

   if (master_draws) then

! if the grid has changed since last sent to the graphics process ...

      if (grid_changed) then

! receive the individual grids from all the slaves and build the composite grid

         call merge_grids(grid,procs,imess,ni,rmess,nr,still_sequential)

! send the composite grid to the master's graphics process

         if (PARALLEL == SEQUENTIAL) then
            call sequential_send(imess(1:ni),ni,rmess(1:nr),nr)
         else
            call phaml_send(procs,graphics_proc(procs),imess,ni,rmess,nr,101)
         endif
         deallocate(imess,rmess,stat=astat)

! if the grid has not changed, just send a message to redraw

      else

         imess1(1) = GRAPHICS_GRID
         if (PARALLEL == SEQUENTIAL) then
            call sequential_send(imess1,1,noreals,0)
         else
            call phaml_send(procs,graphics_proc(procs),imess1,1,noreals,0,101)
         endif

      endif
   endif

endif ! slave or master

grid_changed = .false. ! says the graphics data has been sent since last change

! pause, if requested

call pause_until_enter(procs,io_control%pause_after_draw,dont_pausewatch=.true.)

call end_pause_watch(all_watches)

end subroutine draw_grid

!          ---------
subroutine pack_grid(grid,procs,ref_control,for_master,still_sequential, &
                     imess,rmess,ni,nr)
!          ---------

!----------------------------------------------------
! This subroutine packs the message containing the grid and flags for
! the graphics processor
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type (proc_info), intent(in) :: procs
type (refine_options), intent(in) :: ref_control
logical, intent(in) :: for_master, still_sequential
integer, intent(out) :: imess(:)
real (my_real), intent(out) :: rmess(:)
integer, intent(out) :: ni,nr

!----------------------------------------------------
! Local variables:

logical(small_logical), allocatable :: sendit(:), sendite(:), seen_edge(:)
character(len=HOSTLEN) :: host
integer :: i,k,ind,lev,elem,edge,rind,vert,elem_per_proc,edge_per_proc, &
           vert_per_proc,nproc,my_processor,melem,medge,mvert,astat,ssize, &
           neigh_lid(EDGES_PER_ELEMENT)
type(hash_key) :: boundary_gid

!----------------------------------------------------
! Begin executable code

if (.not. grid%errind_up2date) then
   call all_error_indicators(grid,ref_control%error_estimator)
endif

boundary_gid = BOUNDARY
my_processor = my_proc(procs)
nproc = num_proc(procs)

if (still_sequential .and. for_master .and. my_processor/=1) return

! if packing for the master, determine which elements and edges to send to the
! master as all that I own or have a descendent that I own

if (for_master) then
   allocate(sendit(size(grid%element)),sendite(size(grid%edge)),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in pack grid",procs=procs)
      return
   endif
   call set_sendit(grid,sendit,sendite,my_processor)
endif

! pack the command number

imess(1) = GRAPHICS_GRID

! determine the dimension for elements, edges and vertices needed on
! the graphics server

melem = 0
medge = 0
mvert = 0
do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      melem = max(melem,elem)
      do i=1,EDGES_PER_ELEMENT
         medge = max(medge,grid%element(elem)%edge(i))
      end do
      elem = grid%element(elem)%next
   end do
   vert = grid%head_level_vert(lev)
   do while (vert /= END_OF_LIST)
      mvert = max(mvert,vert)
      vert = grid%vertex(vert)%next
   end do
end do

if (for_master) then
   if (still_sequential) then
      elem_per_proc=10**int(log10(real(melem,my_real))+1)
      edge_per_proc=10**int(log10(real(medge,my_real))+1)
      vert_per_proc=10**int(log10(real(mvert,my_real))+1)
   else
      elem_per_proc = &
              10**int(log10(phaml_global_max(procs,real(melem,my_real),730))+1)
      edge_per_proc = &
              10**int(log10(phaml_global_max(procs,real(medge,my_real),735))+1)
      vert_per_proc = &
              10**int(log10(phaml_global_max(procs,real(mvert,my_real),740))+1)
   endif
   imess(2) = (nproc+1)*elem_per_proc
   imess(3) = (nproc+1)*edge_per_proc
   imess(4) = (nproc+1)*vert_per_proc
else
   imess(2) = melem
   imess(3) = medge
   imess(4) = mvert
endif
imess(5) = size(grid%vertex_solution,2)
imess(6) = size(grid%vertex_solution,3)

! pack the processor info

imess(7) = nproc
imess(8) = my_processor
ind = 8
host = hostname(procs)
do i=1,HOSTLEN
   imess(ind+ i) = ichar(host(i:i))
end do
ind = ind+HOSTLEN

! pack grid_type variables

imess(ind+1) = grid%nsoln/grid%system_size ! number of solutions or eigenvectors
imess(ind+2) = grid%system_size
ind = ind+2

! pack number of children of the root

ind = scalar_imess_size(grid)
rind = 0

! pack the elements

do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (for_master) then
         if (.not. sendit(elem)) then
            elem = grid%element(elem)%next
            cycle
         endif
         imess(ind+1) = elem + my_processor*elem_per_proc
      else
         imess(ind+1) = elem
      endif
      ind = ind + 1
      if (associated(grid%element(elem)%solution)) then
         ssize = size(grid%element(elem)%solution)
      else
         ssize = 0
      endif
      imess(ind+1) = ssize
      ind = ind + 1
      call hash_pack_key(grid%element(elem)%gid,imess,ind+1)
      ind = ind + KEY_SIZE
      rmess(rind+1) = maxval(grid%element_errind(elem,:))
      rmess(rind+2) = grid%element(elem)%work
      rind = rind + 2
      if (ssize /= 0) then
         rmess(rind+1:rind+ssize) = reshape(grid%element(elem)%solution,(/ssize/))
         rind = rind + ssize
         if (grid%have_true) then
            rmess(rind+1:rind+ssize) = reshape(grid%element(elem)%exact,(/ssize/))
         else
            rmess(rind+1:rind+ssize) = 0.0_my_real
         endif
         where (rmess(rind+1:rind+ssize) == huge(0.0_my_real))
            rmess(rind+1:rind+ssize) = 0.0_my_real
         endwhere
         rind = rind + ssize
      endif
      do i=1,VERTICES_PER_ELEMENT
         call hash_pack_key(grid%vertex(grid%element(elem)%vertex(i))%gid, &
                            imess,ind+1+(i-1)*KEY_SIZE)
      end do
      ind = ind + VERTICES_PER_ELEMENT*KEY_SIZE
      do i=1,EDGES_PER_ELEMENT
         call hash_pack_key(grid%edge(grid%element(elem)%edge(i))%gid,imess, &
                           ind+1+(i-1)*KEY_SIZE)
      end do
      ind = ind + EDGES_PER_ELEMENT*KEY_SIZE
      imess(ind+1) = grid%element(elem)%degree
      ind = ind + 1
      imess(ind+1) = grid%element(elem)%level
      ind = ind + 1
      call hash_pack_key(grid%vertex(grid%element(elem)%in)%gid,imess, &
                         ind+1)
      ind = ind + KEY_SIZE
      call hash_pack_key(grid%vertex(grid%element(elem)%out)%gid,imess, &
                         ind+1)
      ind = ind + KEY_SIZE
      imess(ind+1:ind+MAX_CHILD) = grid%element(elem)%order
      ind = ind + MAX_CHILD
      if (grid%element(elem)%isleaf) then
         imess(ind+1) = 1
      else
         imess(ind+1) = 0
      endif
      ind = ind + 1
      if (grid%element(elem)%iown) then
         imess(ind+1) = my_proc(procs)
      else
         imess(ind+1) = 0
      endif
      ind = ind + 1
      neigh_lid = get_neighbors(elem,grid)
      do i=1,EDGES_PER_ELEMENT
         if (neigh_lid(i) == BOUNDARY) then
            call hash_pack_key(boundary_gid,imess,ind+1+(i-1)*KEY_SIZE)
         else
            call hash_pack_key(grid%element(neigh_lid(i))%gid,imess, &
                               ind+1+(i-1)*KEY_SIZE)
         endif
      end do
      ind = ind + EDGES_PER_ELEMENT*KEY_SIZE
      elem = grid%element(elem)%next
   end do
end do

imess(ind+1) = END_OF_ELEMENTS
ind = ind+1

! pack the edges

allocate(seen_edge(medge),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in pack grid",procs=procs)
   return
endif
seen_edge = .false.

do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      do k=1,EDGES_PER_ELEMENT
         edge = grid%element(elem)%edge(k)
         if (seen_edge(edge)) cycle
         seen_edge(edge) = .true.
         if (associated(grid%edge(edge)%solution)) then
            ssize = size(grid%edge(edge)%solution)
         else
            ssize = 0
         endif
         if (for_master) then
            if (.not. sendite(edge)) cycle
            imess(ind+1) = edge + edge_per_proc*my_processor
            imess(ind+2) = ssize
            imess(ind+3) = grid%edge(edge)%assoc_elem +my_processor*elem_per_proc
         else
            imess(ind+1) = edge
            imess(ind+2) = ssize
            imess(ind+3) = grid%edge(edge)%assoc_elem
         endif
         ind = ind + 3
         call hash_pack_key(grid%edge(edge)%gid,imess,ind+1)
         ind = ind + KEY_SIZE
         imess(ind+1) = grid%edge(edge)%degree
         ind = ind + 1
         if (ssize /= 0) then
            rmess(rind+1:rind+ssize) = reshape(grid%edge(edge)%solution,(/ssize/))
            rind = rind + ssize
            if (grid%have_true) then
               rmess(rind+1:rind+ssize) = reshape(grid%edge(edge)%exact,(/ssize/))
            else
               rmess(rind+1:rind+ssize) = 0.0_my_real
            endif
            where (rmess(rind+1:rind+ssize) == huge(0.0_my_real))
               rmess(rind+1:rind+ssize) = 0.0_my_real
            endwhere
            rind = rind + ssize
         endif
      end do
      elem = grid%element(elem)%next
   end do
end do

deallocate(seen_edge,stat=astat)
imess(ind+1) = END_OF_EDGES
ind = ind+1

! pack the vertices

do lev=1,grid%nlev
   vert = grid%head_level_vert(lev)
   do while (vert /= END_OF_LIST)
      if (for_master) then
         if (.not. grid%element(grid%vertex(vert)%assoc_elem)%iown) then
            vert = grid%vertex(vert)%next
            cycle
         endif
         imess(ind+1) = vert + vert_per_proc*my_processor
         imess(ind+2) = grid%vertex(vert)%assoc_elem +my_processor*elem_per_proc
      else
         imess(ind+1) = vert
         imess(ind+2) = grid%vertex(vert)%assoc_elem
      endif
      ind = ind + 2
      call hash_pack_key(grid%vertex(vert)%gid,imess,ind+1)
      ind = ind + KEY_SIZE
      rmess(rind+1) = grid%vertex(vert)%coord%x
      rmess(rind+2) = grid%vertex(vert)%coord%y
      rind = rind + 2
      rmess(rind+1:rind+grid%nsoln) = reshape(grid%vertex_solution(vert,:,:), &
                                              (/grid%nsoln/))
      rind = rind + grid%nsoln
      if (associated(grid%vertex_exact)) then
         where (reshape(grid%vertex_exact(vert,:,:),(/grid%nsoln/)) == huge(0.0_my_real))
            rmess(rind+1:rind+grid%nsoln) = 0.0_my_real
         elsewhere
            rmess(rind+1:rind+grid%nsoln) = reshape(grid%vertex_exact(vert,:,:), &
                                                    (/grid%nsoln/))
         endwhere
      else
         rmess(rind+1:rind+grid%nsoln) = 0.0_my_real
      endif
      rind = rind + grid%nsoln
      vert = grid%vertex(vert)%next
   end do
end do
imess(ind+1) = END_OF_VERTICES

ni = ind+1
nr = rind

if (ni > size(imess) .or. nr > size(rmess)) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("Incorrect allocation for message to graphics process.", &
              intlist=(/ni,size(imess),nr,size(rmess)/))
   stop
endif

if (for_master) deallocate(sendit,sendite,stat=astat)

end subroutine pack_grid

!          -----------------
! function scalar_imess_size and friends
!          -----------------
! These functions give the size of the message containing the grid.  They
! give the amount of integer and real space used for scalar information and
! for each element, edge and vertex.

integer function scalar_imess_size(grid)
type(grid_type), intent(in) :: grid
scalar_imess_size = 10 + HOSTLEN
end function scalar_imess_size

integer function elem_imess_size(grid)
type(grid_type), intent(in) :: grid
elem_imess_size = 6 + MAX_CHILD + &
                 (3 + 2*EDGES_PER_ELEMENT + VERTICES_PER_ELEMENT)*KEY_SIZE
end function elem_imess_size

integer function edge_imess_size(grid)
type(grid_type), intent(in) :: grid
edge_imess_size = 4 + KEY_SIZE
end function edge_imess_size

integer function vert_imess_size(grid)
type(grid_type), intent(in) :: grid
vert_imess_size = 2 + KEY_SIZE
end function vert_imess_size

integer function scalar_rmess_size(grid)
type(grid_type), intent(in) :: grid
logical(small_logical) :: seen_edge(size(grid%edge))
integer :: lev, elem, edge, i
seen_edge = .false.
scalar_rmess_size = 0
do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (associated(grid%element(elem)%solution)) then
         scalar_rmess_size = scalar_rmess_size + &
                             2*size(grid%element(elem)%solution)
      endif
      do i=1,EDGES_PER_ELEMENT
         edge = grid%element(elem)%edge(i)
         if (seen_edge(edge)) cycle
         seen_edge(edge) = .true.
         if (associated(grid%edge(edge)%solution)) then
            scalar_rmess_size = scalar_rmess_size + &
                                2*size(grid%edge(edge)%solution)
         endif
      end do
      elem = grid%element(elem)%next
   end do
end do
end function scalar_rmess_size

integer function elem_rmess_size(grid)
type(grid_type), intent(in) :: grid
elem_rmess_size = 2
end function elem_rmess_size

integer function edge_rmess_size(grid)
type(grid_type), intent(in) :: grid
edge_rmess_size = 0
end function edge_rmess_size

integer function vert_rmess_size(grid)
type(grid_type), intent(in) :: grid
vert_rmess_size = 2 + 2*grid%nsoln
end function vert_rmess_size

!          -----------
subroutine merge_grids(grid,procs,imess,ni,rmess,nr,still_sequential)
!          -----------

!----------------------------------------------------
! This routine is called by the master to merge the grids from each slave
! into a single grid for the master's graphics processor
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
type(proc_info), intent(in) :: procs
integer, pointer :: imess(:)
real (my_real), pointer :: rmess(:)
integer, intent(out) :: ni,nr
logical, intent(in) :: still_sequential
!----------------------------------------------------

!----------------------------------------------------
! Local variables:

character(len=HOSTLEN) :: host
integer :: i, j, iind, rind, piind, prind, np, isize, rsize, astat, ssize
integer, allocatable :: ord(:)
type recv_type
   integer, pointer :: imess(:)
   real(my_real), pointer :: rmess(:)
   integer :: ni,nr,proc,iind,rind
end type recv_type
type (recv_type), allocatable :: recv(:)
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

! receive all the slave messages

if (still_sequential) then
   np = 1
else
   np = num_proc(procs)
endif

allocate(recv(np),ord(np),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in merge grids",procs=procs)
   return
endif
do i=1,np
   call phaml_recv(procs,recv(i)%proc,recv(i)%imess,recv(i)%ni,recv(i)%rmess, &
                   recv(i)%nr,720)
   ord(recv(i)%proc) = i
end do

allocate(imess(sum(recv%ni)),rmess(sum(recv%nr)),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in merge grids",procs=procs)
   return
endif

! pack the scalar information from processor 1; the others should be the same

imess(1:scalar_imess_size(grid)) = recv(1)%imess(1:scalar_imess_size(grid))

! change the processor to this processor

imess(8) = my_proc(procs)
iind = 8
host = hostname(procs)
do i=1,HOSTLEN
   imess(iind+ i) = ichar(host(i:i))
end do
iind = scalar_imess_size(grid)
rind = 0

! pack elements from each processor

isize = elem_imess_size(grid)
rsize = elem_rmess_size(grid)

do j=1,np
   i = ord(j)
   piind = scalar_imess_size(grid)
   prind = 0
   do while (recv(i)%imess(piind+1) /= END_OF_ELEMENTS)
      imess(iind+1:iind+isize) = recv(i)%imess(piind+1:piind+isize)
      ssize = imess(iind+2)
      iind = iind   + isize
      piind = piind + isize
      rmess(rind+1:rind+rsize+2*ssize) = recv(i)%rmess(prind+1:prind+rsize+2*ssize)
      rind = rind   + rsize + 2*ssize
      prind = prind + rsize + 2*ssize
   end do
   recv(i)%iind = piind + 1 ! bookmark for when I copy the edges
   recv(i)%rind = prind
end do
imess(iind+1) = END_OF_ELEMENTS
iind = iind + 1

! pack edges from each processor

isize = edge_imess_size(grid)
rsize = edge_rmess_size(grid)

do j=1,np
   i = ord(j)
   piind = recv(i)%iind
   prind = recv(i)%rind
   do while (recv(i)%imess(piind+1) /= END_OF_EDGES)
      imess(iind+1:iind+isize) = recv(i)%imess(piind+1:piind+isize)
      ssize = imess(iind+2)
      iind = iind   + isize
      piind = piind + isize
      rmess(rind+1:rind+rsize+2*ssize) = recv(i)%rmess(prind+1:prind+rsize+2*ssize)
      rind = rind   + rsize + 2*ssize
      prind = prind + rsize + 2*ssize
   end do
   recv(i)%iind = piind + 1 ! bookmark for when I copy the vertices
   recv(i)%rind = prind
end do
imess(iind+1) = END_OF_EDGES
iind = iind + 1

! pack vertices from each processor

do j=1,np
   i = ord(j)
   piind = recv(i)%iind
   prind = recv(i)%rind
   imess(iind+1:iind+recv(i)%ni-1-piind) = recv(i)%imess(piind+1:recv(i)%ni-1)
   rmess(rind+1:rind+recv(i)%nr-prind) = recv(i)%rmess(prind+1:recv(i)%nr)
   iind = iind+recv(i)%ni-1-piind
   rind = rind+recv(i)%nr-prind
end do
imess(iind+1) = END_OF_VERTICES
ni = iind + 1
nr = rind

do i=1,np
   if (recv(i)%ni /= 0) deallocate(recv(i)%imess,stat=astat)
   if (recv(i)%nr /= 0) deallocate(recv(i)%rmess,stat=astat)
end do
deallocate(recv,ord,stat=astat)

return
end subroutine merge_grids

!          ----------
subroutine set_sendit(grid,sendit,sendite,my_processor)
!          ----------

!----------------------------------------------------
! This routine sets flags to indicate whether or not to send an element
! to the master for the composite grid graphics.  Send it if it or any
! of its descendents are owned by my_processor
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
logical(small_logical), intent(inout) :: sendit(:), sendite(:)
integer, intent(in) :: my_processor
!----------------------------------------------------
! Local variables:

integer :: elem
!----------------------------------------------------
! Begin executable code

elem = grid%head_level_elem(1)
do while (elem /= END_OF_LIST)
   call set_sendit_recur(grid,sendit,sendite,my_processor,elem)
   elem = grid%element(elem)%next
end do
end subroutine set_sendit

!                    ----------------
recursive subroutine set_sendit_recur(grid,sendit,sendite,my_processor,subroot)
!                    ----------------

!----------------------------------------------------
! This routine recursively traverses the tree to do the work of set_sendit
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
logical(small_logical), intent(inout) :: sendit(:), sendite(:)
integer, intent(in) :: my_processor, subroot
!----------------------------------------------------
! Local variables:

integer :: child(MAX_CHILD), i, allc(MAX_CHILD), edge, cedge
!----------------------------------------------------
! Begin executable code

allc = ALL_CHILDREN
child = get_child_lid(grid%element(subroot)%gid,allc,grid%elem_hash)

if (child(1) == NO_CHILD) then

! for a leaf, return whether or not I own it

   sendit(subroot) = grid%element(subroot)%iown
   do i=1,EDGES_PER_ELEMENT
      edge = grid%element(subroot)%edge(i)
      sendite(edge) = grid%element(grid%edge(edge)%assoc_elem)%iown
   end do

else

! otherwise, see if any of the children are owned
! RESTRICTION bisected triangles; only the third edge needs children checked

   sendit(subroot) = .false.
   edge = grid%element(subroot)%edge(3)
   sendite(edge) = .false.
   do i=1,MAX_CHILD
      if (child(i) /= NO_CHILD) then
         call set_sendit_recur(grid,sendit,sendite,my_processor,child(i))
         sendit(subroot) = sendit(subroot) .or. sendit(child(i))
         cedge = grid%element(child(i))%edge(2)
         if (sendite(cedge)) sendite(edge) = .true.
      endif
   end do

endif

end subroutine set_sendit_recur

!          -----------------------
subroutine child_order_from_zoltan(grid,procs,lb,partition_method, &
                                   still_sequential)
!          -----------------------

!----------------------------------------------------
! This routine gets the child order (for the reftree sfc) from Zoltan.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type (proc_info), intent(in) :: procs
type(Zoltan_Struct), pointer :: lb
integer, intent(in) :: partition_method
logical, intent(in) :: still_sequential
!----------------------------------------------------
! Local variables:

integer :: i, nelem_init, maxelem_init, child1, child2, parent, &
           astat, jerr, allc(MAX_CHILD), children(MAX_CHILD)
integer, allocatable :: order(:)
type(hash_key) :: tempkey
!----------------------------------------------------
! Begin executable code

! must be using Zoltan REFTREE as the partition method

if (partition_method /= ZOLTAN_REFTREE) return

! count number of elements in the initial grid

nelem_init = 0
maxelem_init = 0
i = grid%head_level_elem(1)
do while (i /= END_OF_LIST)
   nelem_init = nelem_init + 1
   maxelem_init = max(maxelem_init,i)
   i = grid%element(i)%next
end do

! space for order returned by zoltan and retaining element linked list

allocate(order(KEY_SIZE*(1 + 3*nelem_init + 7*(grid%nelem-grid%nelem_leaf))), &
         stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in child_order_from_zoltan",procs=procs)
   return
endif

! get the order from zoltan

order = -1
call zoltan_child_order(order,jerr,lb,grid,procs,still_sequential)
if (jerr /= 0) then
   call warning("zoltan_child_order returned error.", &
                "Retaining PHAML's child order.",(/jerr/))
   deallocate(order,stat=astat)
   return
endif

! the first 1+nelem_init are the root and it's children (initial grid). Skip
! the root and reorder the initial grid by resetting previous/next

i = KEY_SIZE+1
tempkey = hash_unpack_key(order,i)
child1 = hash_decode_key(tempkey, grid%elem_hash)
if (child1 == HASH_NOT_FOUND) then
   call fatal("zoltan child order returned unknown GID",procs=procs)
   stop
endif
grid%head_level_elem(1) = child1
grid%element(child1)%previous = END_OF_LIST
tempkey = hash_unpack_key(order,i+KEY_SIZE)
grid%element(child1)%in = hash_decode_key(tempkey, grid%vert_hash)
if (grid%element(child1)%in == HASH_NOT_FOUND) then
   call fatal("zoltan child order returned unknown GID for in",procs=procs)
   stop
endif
tempkey = hash_unpack_key(order,i+2*KEY_SIZE)
grid%element(child1)%out = hash_decode_key(tempkey, grid%vert_hash)
if (grid%element(child1)%out == HASH_NOT_FOUND) then
   call fatal("zoltan child order returned unknown GID for out",procs=procs)
   stop
endif
child2 = child1
do i=4*KEY_SIZE+1,(3*nelem_init+1)*KEY_SIZE,3*KEY_SIZE
   tempkey = hash_unpack_key(order,i)
   child1 = hash_decode_key(tempkey, grid%elem_hash)
   if (child1 == HASH_NOT_FOUND) then
      call fatal("zoltan child order returned unknown GID",procs=procs)
      stop
   endif
   grid%element(child2)%next = child1
   grid%element(child1)%previous = child2
   tempkey = hash_unpack_key(order,i+KEY_SIZE)
   grid%element(child1)%in = hash_decode_key(tempkey, grid%vert_hash)
   if (grid%element(child1)%in == HASH_NOT_FOUND) then
      call fatal("zoltan child order returned unknown GID for in",procs=procs)
      stop
   endif
   tempkey = hash_unpack_key(order,i+2*KEY_SIZE)
   grid%element(child1)%out = hash_decode_key(tempkey, grid%vert_hash)
   if (grid%element(child1)%out == HASH_NOT_FOUND) then
      call fatal("zoltan child order returned unknown GID for out",procs=procs)
      stop
   endif
   child2 = child1
end do
grid%element(child2)%next = END_OF_LIST
i = (3*nelem_init+1)*KEY_SIZE + 1

! don't get refinements if still sequential and not processor 1
if (still_sequential .and. my_proc(procs) /= 1) i = size(order)+1

! go through the rest of the list returned by Zoltan and reset order
! to correspond to the order of the children

do while (i < size(order))
! RESTRICTION bisection
   tempkey = hash_unpack_key(order,i)
   if (tempkey == -1) then
      i = i+7*KEY_SIZE
      cycle
   endif
   parent = hash_decode_key(tempkey, grid%elem_hash)
   if (parent == HASH_NOT_FOUND) then
      call fatal("zoltan child order returned unknown GID for parent",procs=procs)
      stop
   endif
   tempkey = hash_unpack_key(order,i+KEY_SIZE)
   child1 = hash_decode_key(tempkey, grid%elem_hash)
   if (child1 == HASH_NOT_FOUND) then
      call fatal("zoltan child order returned unknown GID for child1",procs=procs)
      stop
   endif
   tempkey = hash_unpack_key(order,i+2*KEY_SIZE)
   grid%element(child1)%in = hash_decode_key(tempkey, grid%vert_hash)
   if (grid%element(child1)%in == HASH_NOT_FOUND) then
      call fatal("zoltan child order returned unknown GID for in",procs=procs)
      stop
   endif
   tempkey = hash_unpack_key(order,i+3*KEY_SIZE)
   grid%element(child1)%out = hash_decode_key(tempkey, grid%vert_hash)
   if (grid%element(child1)%out == HASH_NOT_FOUND) then
      call fatal("zoltan child order returned unknown GID for out",procs=procs)
      stop
   endif
   tempkey = hash_unpack_key(order,i+4*KEY_SIZE)
   child2 = hash_decode_key(tempkey, grid%elem_hash)
   if (child2 == HASH_NOT_FOUND) then
      call fatal("zoltan child order returned unknown GID for child2",procs=procs)
      stop
   endif
   tempkey = hash_unpack_key(order,i+5*KEY_SIZE)
   grid%element(child2)%in = hash_decode_key(tempkey, grid%vert_hash)
   if (grid%element(child2)%in == HASH_NOT_FOUND) then
      call fatal("zoltan child order returned unknown GID for in",procs=procs)
      stop
   endif
   tempkey = hash_unpack_key(order,i+6*KEY_SIZE)
   grid%element(child2)%out = hash_decode_key(tempkey, grid%vert_hash)
   if (grid%element(child2)%out == HASH_NOT_FOUND) then
      call fatal("zoltan child order returned unknown GID for out",procs=procs)
      stop
   endif
   allc = ALL_CHILDREN
   children = get_child_lid(grid%element(parent)%gid,allc,grid%elem_hash)
   if (child1 == children(1) .and. child2 == children(2)) then
      grid%element(parent)%order = (/1,2/)
   elseif (child2 == children(1) .and. child1 == children(2)) then
      grid%element(parent)%order = (/2,1/)
   else
      call warning("in draw_grid zoltan did not return the right children", &
                   "parent, children, zoltan", &
                   intlist=(/parent,children,child1,child2/))
   endif
   i = i+7*KEY_SIZE
end do

deallocate(order,stat=astat)

end subroutine child_order_from_zoltan

end module grid_io
