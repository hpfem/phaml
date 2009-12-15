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

!----------------------------------------------------
! This file contains the program for graphics processes.  It consists of
! a module with the routines that start/terminate graphics windows and
! routines that draw the graphics, followed by a short main program, which
! is actually a subroutine.
!
! RESTRICTION 2D
!----------------------------------------------------

module graphics_mod

!----------------------------------------------------
! Other modules used:

use global
use message_passing
use opengl_gl
use opengl_glu
use opengl_glut
use hash_mod
use view_modifier
use gridtype_mod
use evaluate
use phaml_type_mod
!----------------------------------------------------

implicit none

!----------------------------------------------------
! The following types are defined:

!----------------------------------------------------
! The following parameters are defined:

! values for what function to plot on the grid

integer(glcint), parameter :: DRAW_NO_FUNCTION    = 20, &
                              DRAW_SOLUTION       = 21, &
                              DRAW_TRUE           = 22, &
                              DRAW_ERROR          = 23, &
                              DRAW_LEVELS         = 24, &
                              DRAW_ERRIND         = 25

! values for preprocessing a function to draw

integer(glcint), parameter :: PREPROC_NONE = 0, &
                              PREPROC_ABS  = 1, &
                              PREPROC_LOG  = 2, &
                              PREPROC_SQ   = 3, &
                              PREPROC_NEG  = 4

! values for what color to use

integer(glcint), parameter :: COLOR_TRANSPARENT = 0, &
                              COLOR_WHITE       = 1, &
                              COLOR_BLACK       = 2, &
                              COLOR_SOLUTION    = 3, &
                              COLOR_TRUE        = 4, &
                              COLOR_ERROR       = 5, &
                              COLOR_OWNER       = 6, &
                              COLOR_VERT_OWNER  = 7, &
                              COLOR_SIZE        = 8, &
                              COLOR_PART_BOUND  = 9, &
                              COLOR_DEGREE      = 10, &
                              COLOR_ERRIND      = 11

! values for color schemes

integer(glcint), parameter :: SCHEME_RAINBOW        = 0, &
                              SCHEME_GRAY           = 1, &
                              SCHEME_STRIPE         = 2, &
                              SCHEME_DOUBLE_RAINBOW = 3, &
                              SCHEME_STEP_SEQ       = 4

! values for contour location

integer(glcint), parameter :: CONT_XYPLANE = 0, &
                              CONT_SURFACE = 1

! values for what label to place on elements and vertices

integer(glcint), parameter :: LABEL_NOLABEL   = 0, &
                              LABEL_LID       = 1, &
                              LABEL_GID       = 2, &
                              LABEL_PARTITION = 3
! other parameters

integer, parameter :: REFTREE_ROOT = 0, &
                      NO_OWNER = 0
real (my_real), parameter :: myzero  = 0.0_my_real
real (glfloat), parameter :: glfzero  = 0.0_glfloat, &
                             glfhalf  = 0.5_glfloat, &
                             glfone   = 1.0_glfloat
real(gldouble), parameter :: gldzero = 0.0_gldouble, &
                             gldhalf = 0.5_gldouble, &
                             gldone  = 1.0_gldouble
real(gldouble), parameter :: MAX_DEGREE = 23

!----------------------------------------------------
! The following variables are defined:

! the grid data structure

type (grid_type) :: grid

! convenience pointers into the grid data structure

type (element_t), pointer :: element(:)
type (edge_t), pointer ::    edge(:)
type (vertex_t), pointer ::  vertex(:)
real(my_real), pointer :: vertex_solution(:,:,:), vertex_exact(:,:,:)

integer, allocatable ::         elem_owner(:)
type(hash_key), allocatable ::  neighbors(:,:)

! graphics choices

integer :: label_elements = LABEL_NOLABEL, &
           label_verts    = LABEL_NOLABEL, &
           label_edges    = LABEL_NOLABEL, &
           color_elements = COLOR_SOLUTION, &
           color_lines    = COLOR_BLACK, &
           draw_func      = DRAW_SOLUTION, &
           draw_cont      = DRAW_NO_FUNCTION, &
           preproc_func   = PREPROC_NONE, &
           color_scheme   = SCHEME_RAINBOW, &
           step_scheme_steps = 4, &
           step_scheme_hues  = 4, &
           subelement_resolution = 0, &
           contour_location = CONT_XYPLANE
logical :: vert_associated_element=.false., edge_associated_element=.false.

! cropping region

real(my_real) :: xcrop1 = -huge(myzero), &
                 xcrop2 =  huge(myzero), &
                 ycrop1 = -huge(myzero), &
                 ycrop2 =  huge(myzero)

! contour plot parameters

logical :: contour_values_given  = .false.
integer :: num_contour = 21
real(my_real), allocatable :: actual_contours(:)

! parameters for exploding the grid into partitions

type(point), allocatable :: part_cent(:), explode_shift(:)
real(gldouble) :: old_explode_factor = 1.0_gldouble

! which lights are used

logical :: leftlighton=.false., rightlighton=.true., toplighton=.false., &
           bottomlighton=.false., movelighton=.false.

! useful min and max values

real(my_real) :: xmin, xmax, ymin, ymax, zmin, zmax, &
                 minsolut, maxsolut, maxabssolut, mintrue,  maxtrue, &
                 maxabstrue, minerror, maxerror, minerrind, maxerrind, &
                 minabserr, maxabserr, minsize, maxsize, &
                 maxdomain
real(my_real), allocatable :: all_minsolut(:,:),all_maxsolut(:,:), &
                              all_mintrue (:,:),all_maxtrue (:,:), &
                              all_minerror(:,:),all_maxerror(:,:), &
                              all_maxabserr(:,:)

! components are scaled individually or all the same

logical :: indiv_compnt_scale = .true.

! other parameters for control of drawing the grid

logical :: draw_axes_flag  = .false., &
           draw_sfc_flag   = .false., &
           draw_sfcio_flag = .false., &
           grid_newview    = .true.
integer :: eigen          = 1, &
           compnt         = 1
real(glfloat) :: offset   = 10.0_glfloat

! GLUT identifiers

integer(glcint) :: grid_win, grid_menu, eigen_menu=-1, component_menu=-1, &
                   scheme_menu, eigen_submenu, eigen_10s_menus(99)
integer(gluint) :: grid_list=1

! misc variables

integer :: my_processor, nproc, neigen, ncompnt, menu_neigen, menu_ncompnt
character (len=HOSTLEN) :: host
type (proc_info), pointer :: procs
real(gldouble), allocatable :: owner_color(:,:)
real(gldouble) :: lookat_x, lookat_y, lookat_z
logical :: window_initialized = .false.
integer(glcint) :: intzero=0

! TEMP080109 for exact solution of battery problem.  I should only need this
!            for debugging, i.e. have pde(1) draw the error or pde(2) the true
!            soln, but it is always needed for the battery problem because of
!            the determination of max and min true when unpacking the grid.
!            battery/master.f90 tells graphics to set this to .true.
logical :: TEMP_battery = .false.

!----------------------------------------------------
! Non-module procedures used are:

interface

   function trues(x,y,comp,eigen) ! real (my_real)
   use global
   real (my_real), intent(in) :: x,y
   integer, intent(in) :: comp,eigen
   real (my_real) :: trues
   end function trues

end interface

contains

!          ---------------
subroutine process_message
!          ---------------

!----------------------------------------------------
! This routine checks for a message and processes it
!----------------------------------------------------

!----------------------------------------------------
! Local variables:

integer :: proc, ni, nr, allocstat
integer, pointer :: imess(:)
real (my_real), pointer :: rmess(:)
type(phaml_solution_type) :: phaml_solution
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

! receive a message if there is one

if (PARALLEL == SEQUENTIAL) then
   call sequential_recv(imess,ni,rmess,nr)
else
   call phaml_recv(procs,proc,imess,ni,rmess,nr,101,noblock=.true.)
endif

! first entry in imess gives the job

if (ni > 0) then
   select case (imess(1))

! initialize graphics

   case (GRAPHICS_INIT)
      allocate(grid%head_level_elem(1))
      call hash_table_init(grid%elem_hash,imess(2))
      call hash_table_init(grid%edge_hash,imess(3))
      call hash_table_init(grid%vert_hash,imess(4))
      xmin = rmess(1)
      xmax = rmess(2)
      ymin = rmess(3)
      ymax = rmess(4)
      maxdomain = max(abs(xmax-xmin),abs(ymax-ymin))
      deallocate(imess,rmess,stat=allocstat)
      if (allocstat /= 0) then
         call warning("deallocation failed",intlist=(/allocstat/))
      endif
      call init_windows

! draw grid

   case (GRAPHICS_GRID)
      if (ni > 1) then ! new grid data sent
         call reallocate_grid(imess(2),imess(3),imess(4),imess(5),imess(6))
         call unpack_grid(imess,rmess)
         call set_owner_color
      endif
      call draw_grid
      deallocate(imess,stat=allocstat)
      if (allocstat /= 0) then
         call warning("deallocation failed",intlist=(/allocstat/))
      endif
      if (nr > 0) then
         deallocate(rmess,stat=allocstat)
         if (allocstat /= 0) then
            call warning("deallocation failed",intlist=(/allocstat/))
         endif
      endif

! terminate

   case (GRAPHICS_TERMINATE)

      call reallocate_grid(0,0,0,0,0)

      if (allocated(part_cent)) deallocate(part_cent,stat=allocstat)
      if (allocated(explode_shift)) deallocate(explode_shift,stat=allocstat)
      call glutdestroywindow(grid_win)
      call terminate_comm(procs,.true.)
      deallocate(procs,stat=allocstat)
      deallocate(imess,stat=allocstat)
      stop

   case (9) ! update usermod

      phaml_solution%procs = procs
      call update_usermod(phaml_solution)
      if (ni > 0) deallocate(imess)

   case (10) ! increase process universe

      call increase_universe(procs,imess)
      deallocate(imess,stat=allocstat)

   case (11) ! decrease process universe

      call decrease_universe(procs)
      deallocate(imess,stat=allocstat)

   case (12) ! save as postscript
      call write_postscript(intzero)
      deallocate(imess,stat=allocstat)

   case(13) ! TEMP080109 for exact solution of battery problem.
      TEMP_battery = .true.
      deallocate(imess,stat=allocstat)

   case(21) ! TEMP080814 scale for time dependent Schroedinger
      tds_scale = imess(2)
      deallocate(imess,stat=allocstat)

   end select
endif

return
end subroutine process_message
 
!          ---------------
subroutine reallocate_grid(nelem,nedge,nvert,soln2,soln3)
!          ---------------

!----------------------------------------------------
! This routine reallocates the grid arrays
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: nelem, nedge, nvert, soln2, soln3
!----------------------------------------------------
! Local variables:

integer :: i, allocstat
!----------------------------------------------------
! Begin executable code

! Nullify convenience pointers

nullify(element,edge,vertex,vertex_solution,vertex_exact)

! Deallocate memory that is already allocated

if (associated(grid%element_errind)) &
   deallocate(grid%element_errind,stat=allocstat)
if (associated(grid%element)) then
   do i=1,size(grid%element)
      if (associated(grid%element(i)%solution)) &
         deallocate(grid%element(i)%solution,stat=allocstat)
      if (associated(grid%element(i)%exact)) &
         deallocate(grid%element(i)%exact,stat=allocstat)
! TEMP071219 for exact solution for the battery problem
      if (TEMP_battery) then
         if (associated(grid%element(i)%oldsoln)) &
            deallocate(grid%element(i)%oldsoln,stat=allocstat)
         nullify(grid%element(i)%oldsoln)
      endif
! end TEMP071219
   end do
   deallocate(grid%element,stat=allocstat)
endif
if (associated(grid%edge)) then
   do i=1,size(grid%edge)
      if (associated(grid%edge(i)%solution)) &
         deallocate(grid%edge(i)%solution,stat=allocstat)
      if (associated(grid%edge(i)%exact)) &
         deallocate(grid%edge(i)%exact,stat=allocstat)
! TEMP071219 for exact solution for the battery problem
      if (TEMP_battery) then
         if (associated(grid%edge(i)%oldsoln)) &
            deallocate(grid%edge(i)%oldsoln,stat=allocstat)
         nullify(grid%edge(i)%oldsoln)
      endif
! end TEMP071219
   end do
   deallocate(grid%edge,stat=allocstat)
endif
if (associated(grid%vertex)) deallocate(grid%vertex,stat=allocstat)
if (associated(grid%vertex_solution)) deallocate(grid%vertex_solution)
if (associated(grid%vertex_exact)) deallocate(grid%vertex_exact)
! TEMP071219 for exact solution for the battery problem
if (TEMP_battery) then
   if (associated(grid%vertex_oldsoln)) deallocate(grid%vertex_oldsoln)
   nullify(grid%vertex_oldsoln)
endif
! end TEMP071219
if (allocated(elem_owner)) deallocate(elem_owner,stat=allocstat)
if (allocated(neighbors)) deallocate(neighbors,stat=allocstat)

! allocate to new sizes, if the new size is not 0

if (nelem > 0) then
   allocate(grid%element(nelem),elem_owner(nelem), &
            neighbors(EDGES_PER_ELEMENT,nelem), &
            stat=allocstat)
   if (allocstat /= 0) then
      call fatal("allocation failed for drawing grid",intlist=(/allocstat/))
      stop
   endif
   allocate(grid%element_errind(nelem,1),stat=allocstat)
   if (allocstat /= 0) then
      call fatal("allocation failed for drawing grid",intlist=(/allocstat/))
      stop
   endif
endif

if (nedge > 0) then
   allocate(grid%edge(nedge),stat=allocstat)
   if (allocstat /= 0) then
      call fatal("allocation failed for drawing grid",intlist=(/allocstat/))
      stop
   endif
endif

if (nvert > 0) then
   allocate(grid%vertex(nvert),grid%vertex_solution(nvert,soln2,soln3), &
            grid%vertex_exact(nvert,soln2,soln3),stat=allocstat)
   if (allocstat /= 0) then
      call fatal("allocation failed for drawing grid",intlist=(/allocstat/))
      stop
   endif
endif

! nullify components to be allocated later

do i=1,nelem
   nullify(grid%element(i)%solution,grid%element(i)%exact)
! TEMP071219 for exact solution for the battery problem
   if (TEMP_battery) nullify(grid%element(i)%oldsoln)
! end TEMP071219
end do
do i=1,nedge
   nullify(grid%edge(i)%solution,grid%edge(i)%exact)
! TEMP071219 for exact solution for the battery problem
   if (TEMP_battery) nullify(grid%edge(i)%oldsoln)
! end TEMP071219
end do

! replace hash tables

call hash_table_destroy(grid%elem_hash)
call hash_table_destroy(grid%edge_hash)
call hash_table_destroy(grid%vert_hash)
if (nelem > 0) call hash_table_init(grid%elem_hash,nelem)
if (nedge > 0) call hash_table_init(grid%edge_hash,nedge)
if (nvert > 0) call hash_table_init(grid%vert_hash,nvert)

! set convenience pointers into grid data structrue

if (nelem > 0) element => grid%element
if (nedge > 0) edge => grid%edge
if (nvert > 0) then
   vertex => grid%vertex
   vertex_solution => grid%vertex_solution
   vertex_exact => grid%vertex_exact
endif

end subroutine reallocate_grid

!          -----------
subroutine unpack_grid(imess,rmess)
!          -----------

!----------------------------------------------------
! This subroutine unpacks the message containing the grid and flags 
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: imess(:)
real (my_real), intent(in) :: rmess(:)
!----------------------------------------------------
! Local variables:

integer :: ind,rind,elem,vert,i,allocstat,other_lid,part,oldvert,newvert, &
           num_elem,iedge,ssize
real (my_real) :: perim,xcent,ycent,valcent(1)
logical :: duplicate
type(hash_key), allocatable :: vert_gid(:,:), edge_gid(:,:), invert(:), &
                               outvert(:)
real(my_real), allocatable :: sum_perim(:)
integer, pointer :: head_temp(:)

!----------------------------------------------------
! Begin executable code

! initialize min and max values

xmin      =  huge(myzero)
xmax      = -huge(myzero)
ymin      =  huge(myzero)
ymax      = -huge(myzero)
zmin      =  myzero
zmax      =  myzero
minsolut  =  huge(myzero)
maxsolut  = -huge(myzero)
mintrue   =  huge(myzero)
maxtrue   = -huge(myzero)
minerror  =  huge(myzero)
maxerror  = -huge(myzero)
minerrind =  huge(myzero)
maxerrind = -huge(myzero)
minsize   =  huge(myzero)
maxsize   = -huge(myzero)

allocate(vert_gid(VERTICES_PER_ELEMENT,size(element)), &
         edge_gid(EDGES_PER_ELEMENT,size(element)), &
         invert(size(element)), outvert(size(element)), stat=allocstat)
if (allocstat /= 0) then
   call fatal("allocation failed for unpacking grid",intlist=(/allocstat/))
   stop
endif

! unpack the processor info

nproc = imess(7)
my_processor  = imess(8)
ind = 8
host = " "
do i=1,HOSTLEN
   host(i:i) = char(imess(ind+i))
end do
ind = ind + HOSTLEN

! unpack number of eigenfunctions and components (system_size)

neigen = imess(ind+1)
if (eigen > neigen) eigen = neigen
ncompnt = imess(ind+2)
grid%system_size = ncompnt
ind = ind+2

! determine the eigenfunction and component to use

if (compnt > ncompnt) compnt = ncompnt
if (allocated(all_minsolut)) deallocate(all_minsolut,stat=allocstat)
if (allocated(all_maxsolut)) deallocate(all_maxsolut,stat=allocstat)
if (allocated(all_mintrue)) deallocate(all_mintrue,stat=allocstat)
if (allocated(all_maxtrue)) deallocate(all_maxtrue,stat=allocstat)
if (allocated(all_minerror)) deallocate(all_minerror,stat=allocstat)
if (allocated(all_maxerror)) deallocate(all_maxerror,stat=allocstat)
if (allocated(all_maxabserr)) deallocate(all_maxabserr,stat=allocstat)
allocate(all_minsolut(ncompnt,neigen), all_maxsolut(ncompnt,neigen), &
         all_mintrue(ncompnt,neigen), all_maxtrue(ncompnt,neigen), &
         all_minerror(ncompnt,neigen), all_maxerror(ncompnt,neigen), &
         all_maxabserr(ncompnt,neigen), stat=allocstat)
if (allocstat /= 0) then
   call fatal("allocation failed for min/max solut",intlist=(/allocstat/))
   stop
endif
all_minsolut  =  huge(myzero)
all_mintrue   =  huge(myzero)
all_minerror  =  huge(myzero)
all_maxsolut  = -huge(myzero)
all_maxtrue   = -huge(myzero)
all_maxerror  = -huge(myzero)
all_maxabserr = -huge(myzero)

! unpack the elements

element%level = -999
num_elem = 0
grid%head_level_elem = END_OF_LIST
grid%nlev = 0
rind = 0
elem = imess(ind+1)
do while (elem /= END_OF_ELEMENTS)
   ind = ind + 1
   ssize = imess(ind+1)
   ind = ind + 1
   element(elem)%gid = hash_unpack_key(imess,ind+1)
   ind = ind + KEY_SIZE
   call hash_insert(element(elem)%gid,elem,grid%elem_hash)
   grid%element_errind(elem,:) = rmess(rind+1)/rmess(rind+2)
   rind = rind + 2
   allocate(element(elem)%solution(ssize/(neigen*ncompnt),ncompnt,neigen), &
            element(elem)%exact(ssize/(neigen*ncompnt),ncompnt,neigen), &
            stat=allocstat)
   if (allocstat /= 0) then
      call fatal("allocation failed for solution",intlist=(/allocstat/))
      stop
   endif
   element(elem)%solution=reshape(rmess(rind+1:rind+ssize), &
                                  (/ssize/(neigen*ncompnt),ncompnt,neigen/))
   rind = rind + ssize
   element(elem)%exact=reshape(rmess(rind+1:rind+ssize), &
                               (/ssize/(neigen*ncompnt),ncompnt,neigen/))
   rind = rind + ssize
   do i=1,VERTICES_PER_ELEMENT
      vert_gid(i,elem) = hash_unpack_key(imess,ind+1+(i-1)*KEY_SIZE)
   end do
   ind = ind + VERTICES_PER_ELEMENT*KEY_SIZE
   do i=1,EDGES_PER_ELEMENT
      edge_gid(i,elem) = hash_unpack_key(imess,ind+1+(i-1)*KEY_SIZE)
   end do
   ind = ind + EDGES_PER_ELEMENT*KEY_SIZE
   element(elem)%degree = imess(ind+1)
   ind = ind + 1
   element(elem)%level = imess(ind+1)
   ind = ind + 1
   invert(elem) = hash_unpack_key(imess,ind+1)
   ind = ind + KEY_SIZE
   outvert(elem) = hash_unpack_key(imess,ind+1)
   ind = ind + KEY_SIZE
   element(elem)%order = imess(ind+1:ind+MAX_CHILD)
   ind = ind + MAX_CHILD
   element(elem)%isleaf = (imess(ind+1) == 1)
   ind = ind + 1
   elem_owner(elem) = imess(ind+1)
   ind = ind + 1
   do i=1,EDGES_PER_ELEMENT
      neighbors(i,elem) = hash_unpack_key(imess,ind+1+(i-1)*KEY_SIZE)
   end do
   ind = ind + EDGES_PER_ELEMENT*KEY_SIZE
   grid%nlev = max(grid%nlev,element(elem)%level)
   num_elem = max(num_elem,elem)
   if (element(elem)%level > size(grid%head_level_elem)) then
      head_temp => grid%head_level_elem
      allocate(grid%head_level_elem(element(elem)%level))
      grid%head_level_elem = END_OF_LIST
      grid%head_level_elem(1:size(head_temp)) = head_temp
      deallocate(head_temp)
   endif
! avoid duplicates from different processors
   duplicate = .false.
   i = grid%head_level_elem(element(elem)%level)
   do while (i /= END_OF_LIST)
      if (element(i)%gid == element(elem)%gid) then
         duplicate = .true.
         exit
      endif
      i = element(i)%next
   end do
   if (.not. duplicate) then
      element(elem)%next = grid%head_level_elem(element(elem)%level)
      grid%head_level_elem(element(elem)%level) = elem
   endif
   elem = imess(ind+1)
end do
ind = ind + 1

! unpack the edges

iedge = imess(ind+1)
do while (iedge /= END_OF_EDGES)
   ind = ind + 1
   ssize = imess(ind+1)
   ind = ind + 1
   edge(iedge)%assoc_elem = imess(ind+1)
   ind = ind + 1
   edge(iedge)%gid = hash_unpack_key(imess,ind+1)
   other_lid = hash_decode_key(edge(iedge)%gid,grid%edge_hash)
! insert in hash table if this is the first instance or I am the owner
   if (other_lid == HASH_NOT_FOUND) then
      call hash_insert(edge(iedge)%gid,iedge,grid%edge_hash)
   elseif (elem_owner(edge(iedge)%assoc_elem) == my_processor) then
      call hash_remove(edge(iedge)%gid,grid%edge_hash)
      call hash_insert(edge(iedge)%gid,iedge,grid%edge_hash)
   endif
   ind = ind + KEY_SIZE
   edge(iedge)%degree = imess(ind+1)
   ind = ind + 1
   allocate(edge(iedge)%solution(ssize/(neigen*ncompnt),ncompnt,neigen), &
           edge(iedge)%exact(ssize/(neigen*ncompnt),ncompnt,neigen),stat=allocstat)
   if (allocstat /= 0) then
      call fatal("allocation failed for solution",intlist=(/allocstat/))
      stop
   endif
   edge(iedge)%solution=reshape(rmess(rind+1:rind+ssize), &
                                (/ssize/(neigen*ncompnt),ncompnt,neigen/))
   rind = rind + ssize
   edge(iedge)%exact=reshape(rmess(rind+1:rind+ssize), &
                             (/ssize/(neigen*ncompnt),ncompnt,neigen/))
   rind = rind + ssize
   iedge = imess(ind+1)
end do
ind = ind+1

! unpack the vertices

vert = imess(ind+1)
do while (vert /= END_OF_VERTICES)
   ind = ind + 1
   vertex(vert)%assoc_elem = imess(ind+1)
   ind = ind + 1
   vertex(vert)%gid = hash_unpack_key(imess,ind+1)
   other_lid = hash_decode_key(vertex(vert)%gid,grid%vert_hash)
! insert in hash table if this is the first instance or I am the owner
   if (other_lid == HASH_NOT_FOUND) then
      call hash_insert(vertex(vert)%gid,vert,grid%vert_hash)
   elseif (elem_owner(vertex(vert)%assoc_elem) == my_processor) then
      call hash_remove(vertex(vert)%gid,grid%vert_hash)
      call hash_insert(vertex(vert)%gid,vert,grid%vert_hash)
   endif
   ind = ind + KEY_SIZE
   vertex(vert)%coord%x = rmess(rind+1)
   vertex(vert)%coord%y = rmess(rind+2)
   rind = rind + 2
   vertex_solution(vert,:,:) = reshape(rmess(rind+1:rind+neigen*ncompnt), &
                                      (/ncompnt,neigen/))
   rind = rind + neigen*ncompnt
   vertex_exact(vert,:,:) = reshape(rmess(rind+1:rind+neigen*ncompnt), &
                                    (/ncompnt,neigen/))
   rind = rind + neigen*ncompnt
   xmin     = min(xmin,vertex(vert)%coord%x)
   xmax     = max(xmax,vertex(vert)%coord%x)
   ymin     = min(ymin,vertex(vert)%coord%y)
   ymax     = max(ymax,vertex(vert)%coord%y)
   all_minsolut = min(all_minsolut,vertex_solution(vert,:,:))
   all_maxsolut = max(all_maxsolut,vertex_solution(vert,:,:))
   all_mintrue = min(all_mintrue,vertex_exact(vert,:,:))
   all_maxtrue = max(all_maxtrue,vertex_exact(vert,:,:))
   all_minerror = min(all_minerror,vertex_solution(vert,:,:)-vertex_exact(vert,:,:))
   all_maxerror = max(all_maxerror,vertex_solution(vert,:,:)-vertex_exact(vert,:,:))
   vert  = imess(ind+1)
end do

maxdomain = max(abs(xmax-xmin),abs(ymax-ymin))

! TEMP071026 keep a constant z scale, for movie of time dependent Schroedinger
!tds_scale = -1 ! TEMP081007
if (tds_scale == -1) then
   all_minsolut = -3.0_my_real
   all_maxsolut =  3.0_my_real
elseif (tds_scale /= 0) then
   all_maxsolut = max(abs(all_minsolut(1,tds_scale)), &
                      abs(all_maxsolut(1,tds_scale)))
   all_minsolut = -all_maxsolut
endif

xcrop1 = max(xcrop1,xmin-.01_my_real)
xcrop2 = min(xcrop2,xmax+.01_my_real)
ycrop1 = max(ycrop1,ymin-.01_my_real)
ycrop2 = min(ycrop2,ymax+.01_my_real)

! set the edges of the elements from the edge gids,
! set the vertices of the elements from the vertex gids,
! decode the in and out vertices
! find min and max error indicator from owned leaf elements
! find min and max error at element midpoints
! find the min and max element edge size

do elem=1,num_elem
   if (element(elem)%level == -999) cycle
   do i=1,EDGES_PER_ELEMENT
      iedge = hash_decode_key(edge_gid(i,elem),grid%edge_hash)
      if (iedge == HASH_NOT_FOUND) then
         call warning("edge hash code not found; using edge 1; element", &
                      intlist=(/elem/))
         call hash_print_key(edge_gid(i,elem),errunit)
         element(elem)%edge(i) = 1
      else
         element(elem)%edge(i) = iedge
      endif
   end do
   do i=1,VERTICES_PER_ELEMENT
      vert = hash_decode_key(vert_gid(i,elem),grid%vert_hash)
      if (vert == HASH_NOT_FOUND) then
         call warning("vertex hash code not found; using vertex 1; element", &
                      intlist=(/elem/))
         call hash_print_key(vert_gid(i,elem),errunit)
         element(elem)%vertex(i) = 1
      else
         element(elem)%vertex(i) = vert
      endif
   end do
   if (element(elem)%isleaf) then
      if (elem_owner(elem) == my_processor .or.  my_processor==MASTER) then
         minerrind = min(minerrind,grid%element_errind(elem,1))
         maxerrind = max(maxerrind,grid%element_errind(elem,1))
      endif
   endif
   element(elem)%in      = hash_decode_key(invert(elem),grid%vert_hash)
   element(elem)%out     = hash_decode_key(outvert(elem),grid%vert_hash)
   if (element(elem)%isleaf) then
      minsize = min(minsize,sqrt( &
         (vertex(element(elem)%vertex(1))%coord%x-vertex(element(elem)%vertex(2))%coord%x)**2 + &
         (vertex(element(elem)%vertex(1))%coord%y-vertex(element(elem)%vertex(2))%coord%y)**2))
      maxsize = max(maxsize,sqrt( &
         (vertex(element(elem)%vertex(1))%coord%x-vertex(element(elem)%vertex(2))%coord%x)**2 + &
         (vertex(element(elem)%vertex(1))%coord%y-vertex(element(elem)%vertex(2))%coord%y)**2))
      minsize = min(minsize,sqrt( &
         (vertex(element(elem)%vertex(2))%coord%x-vertex(element(elem)%vertex(3))%coord%x)**2 + &
         (vertex(element(elem)%vertex(2))%coord%y-vertex(element(elem)%vertex(3))%coord%y)**2))
      maxsize = max(maxsize,sqrt( &
         (vertex(element(elem)%vertex(2))%coord%x-vertex(element(elem)%vertex(3))%coord%x)**2 + &
         (vertex(element(elem)%vertex(2))%coord%y-vertex(element(elem)%vertex(3))%coord%y)**2))
      minsize = min(minsize,sqrt( &
         (vertex(element(elem)%vertex(3))%coord%x-vertex(element(elem)%vertex(1))%coord%x)**2 + &
         (vertex(element(elem)%vertex(3))%coord%y-vertex(element(elem)%vertex(1))%coord%y)**2))
      maxsize = max(maxsize,sqrt( &
         (vertex(element(elem)%vertex(3))%coord%x-vertex(element(elem)%vertex(1))%coord%x)**2 + &
         (vertex(element(elem)%vertex(3))%coord%y-vertex(element(elem)%vertex(1))%coord%y)**2))
   endif
end do
minsize = log(minsize)
maxsize = log(maxsize)

! TEMP071219 for exact solution for the battery problem
if (TEMP_battery) call exact_to_oldsoln

! Set max and min solution, true and error.
! Evaluate errors at element midpoints to avoid 0 error for initial condition.

do elem=1,num_elem
   if (element(elem)%level == -999) cycle
   if (.not. element(elem)%isleaf) cycle
   if (elem_owner(elem) /= my_processor .and.  my_processor/=MASTER) cycle
   xcent = sum(vertex(element(elem)%vertex)%coord%x)/3
   ycent = sum(vertex(element(elem)%vertex)%coord%y)/3
   call graphics_evaluate_soln((/xcent/),(/ycent/),elem,valcent)
! TEMP only makes the adjustment for the first eigen,compnt
   valcent(1) = valcent(1) - graphics_trues(xcent,ycent)
   if (valcent(1) < all_minerror(1,1)) all_minerror(1,1) = valcent(1)
   if (valcent(1) > all_maxerror(1,1)) all_maxerror(1,1) = valcent(1)
end do

all_maxabserr = max(abs(all_minerror),abs(all_maxerror))
call set_max_min

deallocate(vert_gid,edge_gid,invert,outvert,stat=allocstat)
if (allocstat /= 0) then
   call warning("deallocation failed for vert_gid",intlist=(/allocstat/))
endif

! determine a center of each partition as a weighted average of the
! centers of the elements, using the square of the perimeter of the
! elements as the weights

if (allocated(part_cent)) deallocate(part_cent,stat=allocstat)
if (allocated(explode_shift)) deallocate(explode_shift,stat=allocstat)
allocate(part_cent(nproc),explode_shift(NO_OWNER:nproc), &
         stat=allocstat)
if (allocstat /= 0) then
   call fatal("allocation failed for explode data structures", &
              intlist=(/allocstat/))
   stop
endif
part_cent = point(myzero,myzero)
explode_shift = point(myzero,myzero)
allocate(sum_perim(nproc),stat=allocstat)
if (allocstat /= 0) then
   call fatal("allocation failed for sum_perim",intlist=(/allocstat/))
   stop
endif
sum_perim = myzero
do elem=1,num_elem
   if (element(elem)%level == -999) cycle
   if (.not. element(elem)%isleaf) cycle
   part = elem_owner(elem)
   if (part == NO_OWNER) cycle
   oldvert = element(elem)%vertex(VERTICES_PER_ELEMENT)
   perim = myzero; xcent = myzero; ycent = myzero
   do i=1,VERTICES_PER_ELEMENT
      newvert = element(elem)%vertex(i)
      perim = perim + sqrt( &
               (vertex(newvert)%coord%x-vertex(oldvert)%coord%x)**2 + &
               (vertex(newvert)%coord%y-vertex(oldvert)%coord%y)**2)
      xcent = xcent + vertex(newvert)%coord%x
      ycent = ycent + vertex(newvert)%coord%y
      oldvert = newvert
   end do
   part_cent(part)%x = part_cent(part)%x+perim*perim*xcent/VERTICES_PER_ELEMENT
   part_cent(part)%y = part_cent(part)%y+perim*perim*ycent/VERTICES_PER_ELEMENT
   sum_perim(part) = sum_perim(part) + perim*perim
end do
where(sum_perim /= myzero) part_cent%x = part_cent%x / sum_perim
where(sum_perim /= myzero) part_cent%y = part_cent%y / sum_perim
deallocate(sum_perim,stat=allocstat)
if (allocstat /= 0) then
   call warning("deallocation failed for sum_perim",intlist=(/allocstat/))
endif

! change the menu if the number of eigenfunctions has changed

if (menu_neigen /= neigen) then
   call update_eigen_menu
endif
if (menu_ncompnt /= ncompnt+2) then
   call update_compnt_menu
endif

grid_newview = .true.

return
end subroutine unpack_grid

!          -----------
subroutine set_max_min()
!          -----------

!----------------------------------------------------
! This routine set max and min values for the current eigen and compnt
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

if (allocated(all_minsolut)) then

   select case(compnt)

   case (-1)
      minsolut = 0.0_my_real
      maxsolut = sum(abs(all_maxsolut(:,eigen)))
      maxabssolut = maxsolut
      if (maxabssolut == 0.0_my_real) maxabssolut = 1.0_my_real
      mintrue = 0.0_my_real
      maxtrue = sum(abs(all_maxtrue(:,eigen)))
      maxabstrue = maxtrue
      if (maxabstrue == 0.0_my_real) maxabstrue = 1.0_my_real
      minerror = 0.0_my_real
      maxerror = sum(abs(all_maxerror(:,eigen)))
      maxabserr = sum(abs(all_maxabserr(:,eigen)))
      if (maxabserr == 0.0_my_real) maxabserr = 1.0_my_real

   case (-2)
      minsolut = 0.0_my_real
      maxsolut = sum(all_maxsolut(:,eigen)**2)
      maxabssolut = maxsolut
      if (maxabssolut == 0.0_my_real) maxabssolut = 1.0_my_real
      mintrue = 0.0_my_real
      maxtrue = sum(all_maxtrue(:,eigen)**2)
      maxabstrue = maxtrue
      if (maxabstrue == 0.0_my_real) maxabstrue = 1.0_my_real
      minerror = 0.0_my_real
      maxerror = sum(all_maxerror(:,eigen)**2)
      maxabserr = sum(all_maxabserr(:,eigen)**2)
      if (maxabserr == 0.0_my_real) maxabserr = 1.0_my_real

   case default
      if (indiv_compnt_scale) then
         minsolut = all_minsolut(compnt,eigen)
         maxsolut = all_maxsolut(compnt,eigen)
         maxabssolut = max(abs(minsolut),abs(maxsolut))
         if (maxabssolut == 0.0_my_real) maxabssolut = 1.0_my_real
         mintrue = all_mintrue(compnt,eigen)
         maxtrue = all_maxtrue(compnt,eigen)
         maxabstrue = max(abs(mintrue),abs(maxtrue))
         if (maxabstrue == 0.0_my_real) maxabstrue = 1.0_my_real
         minerror = all_minerror(compnt,eigen)
         maxerror = all_maxerror(compnt,eigen)
         maxabserr = all_maxabserr(compnt,eigen)
         if (maxabserr == 0.0_my_real) maxabserr = 1.0_my_real
      else
         minsolut = minval(all_minsolut(:,eigen))
         maxsolut = maxval(all_maxsolut(:,eigen))
         maxabssolut = max(abs(minsolut),abs(maxsolut))
         if (maxabssolut == 0.0_my_real) maxabssolut = 1.0_my_real
         mintrue = minval(all_mintrue(:,eigen))
         maxtrue = maxval(all_maxtrue(:,eigen))
         maxabstrue = max(abs(mintrue),abs(maxtrue))
         if (maxabstrue == 0.0_my_real) maxabstrue = 1.0_my_real
         minerror = minval(all_minerror(:,eigen))
         maxerror = maxval(all_maxerror(:,eigen))
         maxabserr = maxval(all_maxabserr(:,eigen))
         if (maxabserr == 0.0_my_real) maxabserr = 1.0_my_real
      endif

   end select
endif

end subroutine set_max_min

!          ----------------
subroutine update_eigen_menu
!          ----------------

!----------------------------------------------------
! This routine changes the number of entries in the eigenfunction selection
! menu when the number of eigenfunctions changes
!----------------------------------------------------

!----------------------------------------------------
! Local variables:

integer(glcint) :: hold
!----------------------------------------------------
! Begin executable code

if (eigen_menu == -1) return ! no eigenfunction menu

hold = glutgetmenu()
call glutsetmenu(eigen_menu)

if (menu_neigen > neigen) then
   call decrease_eigen_menu(menu_neigen,neigen)
elseif (menu_neigen < neigen) then
   call increase_eigen_menu(menu_neigen,neigen)
endif

menu_neigen = neigen
call glutsetmenu(hold)
end subroutine update_eigen_menu

!          ------------------
subroutine increase_eigen_menu(oldn,newn)
!          ------------------

!----------------------------------------------------
! This routine increases the number of entries in the eigenfunction selection menu
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments:

integer, intent(in) :: oldn, newn

!----------------------------------------------------
! Local variables:

integer(glcint) :: i, j, oldone, newone
character(len=17) :: str
!----------------------------------------------------
! Begin executable code

! add new entries in primary menu

do i=oldn+1,min(newn,9)
   write(str,"(a14,i3)") "eigenfunction ",i
   call glutaddmenuentry(trim(str),i)
end do

! create submenu for selection of 10's if needed and not already there

if (oldn < 10 .and. newn >= 10) then
   eigen_submenu = glutcreatemenu(select_eigen)
   call glutsetmenu(eigen_menu)
   call glutaddsubmenu("more",eigen_submenu)
endif

! create menu for each of the 10's that are needed

do i=1+oldn/10,newn/10
   eigen_10s_menus(i) = glutcreatemenu(select_eigen)
end do
do i=1+oldn/10,newn/10
   call glutsetmenu(eigen_submenu)
   str = "  0's"
   write(str(1:2),"(i2)") i
   call glutaddsubmenu(trim(str),eigen_10s_menus(i))
end do

! create the menu entries in the 10's menus

oldone = oldn - 10*(oldn/10)
newone = newn - 10*(newn/10)

! the old last 10's menu

if (oldn >= 10) then
   i = oldn/10
   call glutsetmenu(eigen_10s_menus(i))
   do j = oldone+1, min(9,newn-10*i)
      write(str,"(a14,i3)") "eigenfunction ",10*i+j
      call glutaddmenuentry(trim(str),10*i+j)
   end do
endif

! menus before the new last menu

do i=1+oldn/10,newn/10-1
   call glutsetmenu(eigen_10s_menus(i))
   do j=0,9
      write(str,"(a14,i3)") "eigenfunction ",10*i+j
      call glutaddmenuentry(trim(str),10*i+j)
   end do
end do

! new last menu

if (oldn/10 /= newn/10) then
   i = newn/10
   call glutsetmenu(eigen_10s_menus(i))
   do j=0,newone
      write(str,"(a14,i3)") "eigenfunction ",10*i+j
      call glutaddmenuentry(trim(str),10*i+j)
   end do
endif

end subroutine increase_eigen_menu

!          ------------------
subroutine decrease_eigen_menu(oldn,newn)
!          ------------------

!----------------------------------------------------
! This routine decreases the number of entries in the eigenfunction selection menu
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments:

integer, intent(in) :: oldn, newn

!----------------------------------------------------
! Local variables:

integer(glcint) :: i, j, oldone, newone
!----------------------------------------------------
! Begin executable code

oldone = oldn - 10*(oldn/10)
newone = newn - 10*(newn/10)

! remove last menu if not needed

if (oldn >= 10 .and. oldn/10 /= newn/10) then
   i = oldn/10
   call glutsetmenu(eigen_10s_menus(i))
   do j=oldone,0,-1
      call glutremovemenuitem(j+1)
   end do
   call glutsetmenu(eigen_submenu)
   call glutremovemenuitem(i)
   call glutsetmenu(eigen_submenu)
endif

! remove menus between the old and new last menu

do i=oldn/10-1,newn/10+1,-1
   call glutsetmenu(eigen_10s_menus(i))
   do j=9,0,-1
      call glutremovemenuitem(j+1)
   end do
   call glutsetmenu(eigen_submenu)
   call glutremovemenuitem(i)
   call glutdestroymenu(eigen_10s_menus(i))
end do

! remove menu items in the new last menu

if (newn >= 10) then
   i = newn/10
   call glutsetmenu(eigen_10s_menus(i))
   do j = min(9,oldn-10*i), newone+1, -1
      call glutremovemenuitem(j+1)
   end do
endif

! remove submenu for selection of 10's if there and no longer needed

if (oldn >= 10 .and. newn < 10) then
   call glutsetmenu(eigen_menu)
   call glutremovemenuitem(10)
   call glutdestroymenu(eigen_submenu)
endif

! remove entries in primary menu

do i=min(9,oldn),newn+1,-1
   call glutremovemenuitem(i)
end do

end subroutine decrease_eigen_menu

!          ----------------
subroutine update_compnt_menu
!          ----------------

!----------------------------------------------------
! This routine changes the number of entries in the component selection
! menu when the number of components changes
!----------------------------------------------------

!----------------------------------------------------
! Local variables:

integer(glcint) :: i, hold
character(len=15) :: str
!----------------------------------------------------
! Begin executable code

if (component_menu == -1) return ! no component menu

hold = glutgetmenu()
call glutsetmenu(component_menu)

str = ""
if (menu_ncompnt > ncompnt+2) then
   do i=menu_ncompnt,ncompnt+3,-1
      call glutremovemenuitem(i-2)
   end do

elseif (menu_ncompnt < ncompnt+2) then
   do i=menu_ncompnt+1,ncompnt+2
      write(str,"(a10,i3)") "component ",i-2
      call glutaddmenuentry(trim(str),i-2)
   end do
endif

menu_ncompnt = ncompnt+2
call glutsetmenu(hold)
return
end subroutine update_compnt_menu

!          ---------------
subroutine set_owner_color
!          ---------------

!----------------------------------------------------
! This routine defines the color to use for coloring by owner
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

!----------------------------------------------------
! Local variables:

integer :: proc, allocstat
real(gldouble) :: max_val
!----------------------------------------------------
! Begin executable code

if (allocated(owner_color)) then
   deallocate(owner_color,stat=allocstat)
   if (allocstat /= 0) then
      call warning("deallocation failed",intlist=(/allocstat/))
   endif
endif
allocate(owner_color(4,nproc),stat=allocstat)
if (allocstat /= 0) then
   call fatal("allocation failed for owner colors", &
              intlist=(/allocstat/))
   stop
endif

if (color_scheme == SCHEME_STEP_SEQ) then
   max_val = step_scheme_steps*step_scheme_hues
else
   max_val = nproc
endif
do proc=1,nproc
   call get_rainbow(real(proc,gldouble),1.0_gldouble,max_val,owner_color(:,proc))
end do

return
end subroutine set_owner_color

!          -----------------
subroutine preproc_and_scale(selector,val,ppval,ppssval,ppmin,ppssmin,ppmax, &
                             ppssmax)
!          -----------------

!----------------------------------------------------
! This routine converts val to the preprocessed value ppval and the
! preprocessed, scaled and shifted value ppssval
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: selector
real(my_real), intent(in) :: val
real(my_real), intent(out) :: ppval, ppssval, ppmin, ppssmin, ppmax, ppssmax
!----------------------------------------------------
! Local variables:

real(my_real) :: minval, maxval
!----------------------------------------------------
! Begin executable code

select case (preproc_func)

case (PREPROC_NONE)
   ppval = val
   select case (selector)
   case (DRAW_NO_FUNCTION)
      ppmin = 0.0_my_real
      ppmax = 0.0_my_real
      ppssval = 0.5_my_real
      ppssmin = 0.5_my_real
      ppssmax = 0.5_my_real
   case (DRAW_SOLUTION, COLOR_SOLUTION)
      ppmin = minsolut
      ppmax = maxsolut
      if ((maxsolut-minsolut) == 0.0_my_real) then
         ppssval = ppval
         ppssmin = ppmin
         ppssmax = ppmax
      else
         ppssval = (ppval-minsolut)/(maxsolut-minsolut)
         ppssmin = 0.0_my_real
         ppssmax = 1.0_my_real
      endif
   case (DRAW_TRUE, COLOR_TRUE)
      ppmin = mintrue
      ppmax = maxtrue
      if ((maxtrue-mintrue) == 0.0_my_real) then
         ppssval = ppval
         ppssmin = ppmin
         ppssmax = ppmax
      else
         ppssval = (ppval-mintrue)/(maxtrue-mintrue)
         ppssmin = 0.0_my_real
         ppssmax = 1.0_my_real
      endif
   case (DRAW_ERROR, COLOR_ERROR)
      ppmin = minerror
      ppmax = maxerror
      if ((maxerror-minerror) == 0.0_my_real) then
         ppssval = ppval
         ppssmin = ppmin
         ppssmax = ppmax
      else
         ppssval = (ppval-minerror)/(maxerror-minerror)
         ppssmin = 0.0_my_real
         ppssmax = 1.0_my_real
      endif
   case (DRAW_LEVELS)
      ppmin = 0.0_my_real
      ppmax = grid%nlev
      ppssval = val
      ppssmin = ppmin
      ppssmax = ppmax
   case (DRAW_ERRIND, COLOR_ERRIND)
      ppmin = minerrind
      ppmax = maxerrind
      if ((maxerrind-minerrind) == 0.0_my_real) then
         ppssval = ppval
         ppssmin = ppmin
         ppssmax = ppmax
      else
         ppssval = (ppval-minerrind)/(maxerrind-minerrind)
         ppssmin = 0.0_my_real
         ppssmax = 1.0_my_real
      endif
   end select

case (PREPROC_ABS)
   ppval = abs(val)
   select case (selector)
   case (DRAW_NO_FUNCTION)
      ppmin = 0.0_my_real
      ppmax = 0.0_my_real
      ppssval = 0.5_my_real
      ppssmin = 0.5_my_real
      ppssmax = 0.5_my_real
   case (DRAW_SOLUTION, COLOR_SOLUTION)
      if ((maxsolut-minsolut) == 0.0_my_real) then
         ppmin = abs(minsolut)
         ppmax = ppmin
         ppssval = ppval
         ppssmin = ppmin
         ppssmax = ppmin
      else
         if (minsolut >= 0.0_my_real) then
            ppmin = minsolut
            ppmax = maxsolut
            ppssval = (ppval-minsolut)/(maxsolut-minsolut)
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         elseif (maxsolut <= 0.0_my_real) then
            ppmin = -maxsolut
            ppmax = -minsolut
            ppssval = (ppval+maxsolut)/(maxsolut-minsolut)
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         else
            ppmin = 0.0_my_real
            ppmax = max(abs(minsolut),abs(maxsolut))
            ppssval = ppval/max(abs(minsolut),abs(maxsolut))
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         endif
      endif
   case (DRAW_TRUE, COLOR_TRUE)
      if ((maxtrue-mintrue) == 0.0_my_real) then
         ppmin = abs(mintrue)
         ppmax = ppmin
         ppssval = ppval
         ppssmin = ppmin
         ppssmax = ppmin
      else
         if (mintrue >= 0.0_my_real) then
            ppmin = mintrue
            ppmax = maxtrue
            ppssval = (ppval-mintrue)/(maxtrue-mintrue)
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         elseif (maxsolut <= 0.0_my_real) then
            ppmin = -maxtrue
            ppmax = -mintrue
            ppssval = (ppval+maxtrue)/(maxtrue-mintrue)
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         else
            ppmin = 0.0_my_real
            ppmax = max(abs(mintrue),abs(maxtrue))
            ppssval = ppval/max(abs(mintrue),abs(maxtrue))
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         endif
      endif
   case (DRAW_ERROR, COLOR_ERROR)
      if ((maxerror-minerror) == 0.0_my_real) then
         ppmin = abs(minerror)
         ppmax = ppmin
         ppssval = ppval
         ppssmin = ppmin
         ppssmax = ppmin
      else
         if (minerror >= 0.0_my_real) then
            ppmin = minerror
            ppmax = maxerror
            ppssval = (ppval-minerror)/(maxerror-minerror)
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         elseif (maxerror <= 0.0_my_real) then
            ppmin = -maxerror
            ppmax = -minerror
            ppssval = (ppval+maxerror)/(maxerror-minerror)
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         else
            ppmin = 0.0_my_real
            ppmax = max(abs(minerror),abs(maxerror))
            ppssval = ppval/max(abs(minerror),abs(maxerror))
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         endif
      endif
   case (DRAW_LEVELS)
      ppmin = 0.0_my_real
      ppmax = grid%nlev
      ppssval = val
      ppssmin = ppmin
      ppssmax = ppmax
   case (DRAW_ERRIND, COLOR_ERRIND)
      if ((maxerrind-minerrind) == 0.0_my_real) then
         ppmin = abs(minerrind)
         ppmax = ppmin
         ppssval = ppval
         ppssmin = ppmin
         ppssmax = ppmin
      else
         if (minerrind >= 0.0_my_real) then
            ppmin = minerrind
            ppmax = maxerrind
            ppssval = (ppval-minerrind)/(maxerrind-minerrind)
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         elseif (maxerrind <= 0.0_my_real) then
            ppmin = -maxerrind
            ppmax = -minerrind
            ppssval = (ppval+maxerrind)/(maxerrind-minerrind)
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         else
            ppmin = 0.0_my_real
            ppmax = max(abs(minerrind),abs(maxerrind))
            ppssval = ppval/max(abs(minerrind),abs(maxerrind))
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         endif
      endif
   end select

case (PREPROC_NEG)
   ppval = -val
   select case (selector)
   case (DRAW_NO_FUNCTION)
      ppmin = 0.0_my_real
      ppmax = 0.0_my_real
      ppssval = 0.5_my_real
      ppssmin = 0.5_my_real
      ppssmax = 0.5_my_real
   case (DRAW_SOLUTION, COLOR_SOLUTION)
      if ((maxsolut-minsolut) == 0.0_my_real) then
         ppmin = -minsolut
         ppmax = ppmin
         ppssval = ppval
         ppssmin = ppmin
         ppssmax = ppmin
      else
         ppmin = -maxsolut
         ppmax = -minsolut
         ppssval = (ppval+maxsolut)/(maxsolut-minsolut)
         ppssmin = 0.0_my_real
         ppssmax = 1.0_my_real
      endif
   case (DRAW_TRUE, COLOR_TRUE)
      if ((maxtrue-mintrue) == 0.0_my_real) then
         ppmin = -mintrue
         ppmax = ppmin
         ppssval = ppval
         ppssmin = ppmin
         ppssmax = ppmin
      else
         ppmin = -maxtrue
         ppmax = -mintrue
         ppssval = (ppval+maxtrue)/(maxtrue-mintrue)
         ppssmin = 0.0_my_real
         ppssmax = 1.0_my_real
      endif
   case (DRAW_ERROR, COLOR_ERROR)
      if ((maxerror-minerror) == 0.0_my_real) then
         ppmin = -minerror
         ppmax = ppmin
         ppssval = ppval
         ppssmin = ppmin
         ppssmax = ppmin
      else
         ppmin = -maxerror
         ppmax = -minerror
         ppssval = (ppval+maxerror)/(maxerror-minerror)
         ppssmin = 0.0_my_real
         ppssmax = 1.0_my_real
      endif
   case (DRAW_LEVELS)
      ppmin = 0.0_my_real
      ppmax = grid%nlev
      ppssval = val
      ppssmin = ppmin
      ppssmax = ppmax
   case (DRAW_ERRIND, COLOR_ERRIND)
      if ((maxerrind-minerrind) == 0.0_my_real) then
         ppmin = -minerrind
         ppmax = ppmin
         ppssval = ppval
         ppssmin = ppmin
         ppssmax = ppmin
      else
         ppmin = -maxerrind
         ppmax = -minerrind
         ppssval = (ppval+maxerrind)/(maxerrind-minerrind)
         ppssmin = 0.0_my_real
         ppssmax = 1.0_my_real
      endif
   end select

case (PREPROC_SQ)
   ppval = val**2
   select case (selector)
   case (DRAW_NO_FUNCTION)
      ppmin = 0.0_my_real
      ppmax = 0.0_my_real
      ppssval = 0.5_my_real
      ppssmin = 0.5_my_real
      ppssmax = 0.5_my_real
   case (DRAW_SOLUTION, COLOR_SOLUTION)
      if ((maxsolut-minsolut) == 0.0_my_real) then
         ppmin = minsolut**2
         ppmax = ppmin
         ppssval = ppval
         ppssmin = ppmin
         ppssmax = ppmin
      else
         if (minsolut > 0.0_my_real) then
            ppmin = minsolut**2
            ppmax = maxsolut**2
            ppssval = (ppval-minsolut**2)/(maxsolut**2-minsolut**2)
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         elseif (maxsolut < 0.0_my_real) then
            ppmin = maxsolut**2
            ppmax = minsolut**2
            ppssval = (ppval-maxsolut**2)/(maxsolut**2-minsolut**2)
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         else
            ppmin = 0.0_my_real
            ppmax = max(maxsolut**2,minsolut**2)
            ppssval = ppval/max(maxsolut**2,minsolut**2)
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         endif
      endif
   case (DRAW_TRUE, COLOR_TRUE)
      if ((maxtrue-mintrue) == 0.0_my_real) then
         ppmin = mintrue**2
         ppmax = ppmin
         ppssval = ppval
         ppssmin = ppmin
         ppssmax = ppmin
      else
         if (mintrue > 0.0_my_real) then
            ppmin = mintrue**2
            ppmax = maxtrue**2
            ppssval = (ppval-mintrue**2)/(maxtrue**2-mintrue**2)
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         elseif (maxtrue < 0.0_my_real) then
            ppmin = maxtrue**2
            ppmax = mintrue**2
            ppssval = (ppval-maxtrue**2)/(maxtrue**2-mintrue**2)
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         else
            ppmin = 0.0_my_real
            ppmax = max(maxtrue**2,mintrue**2)
            ppssval = ppval/max(maxtrue**2,mintrue**2)
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         endif
      endif
   case (DRAW_ERROR, COLOR_ERROR)
      if ((maxerror-minerror) == 0.0_my_real) then
         ppmin = minerror**2
         ppmax = ppmin
         ppssval = ppval
         ppssmin = ppmin
         ppssmax = ppmin
      else
         if (minerror > 0.0_my_real) then
            ppmin = minerror**2
            ppmax = maxerror**2
            ppssval = (ppval-minerror**2)/(maxerror**2-minerror**2)
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         elseif (maxerror < 0.0_my_real) then
            ppmin = maxerror**2
            ppmax = minerror**2
            ppssval = (ppval-maxerror**2)/(maxerror**2-minerror**2)
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         else
            ppmin = 0.0_my_real
            ppmax = max(maxerror**2,minerror**2)
            ppssval = ppval/max(maxerror**2,minerror**2)
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         endif
      endif
   case (DRAW_LEVELS)
      ppmin = 0.0_my_real
      ppmax = grid%nlev
      ppssval = val
      ppssmin = ppmin
      ppssmax = ppmax
   case (DRAW_ERRIND, COLOR_ERRIND)
      if ((maxerrind-minerrind) == 0.0_my_real) then
         ppmin = minerrind**2
         ppmax = ppmin
         ppssval = ppval
         ppssmin = ppmin
         ppssmax = ppmin
      else
         if (minerrind > 0.0_my_real) then
            ppmin = minerrind**2
            ppmax = maxerrind**2
            ppssval = (ppval-minerrind**2)/(maxerrind**2-minerrind**2)
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         elseif (maxerrind < 0.0_my_real) then
            ppmin = maxerrind**2
            ppmax = minerrind**2
            ppssval = (ppval-maxerrind**2)/(maxerrind**2-minerrind**2)
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         else
            ppmin = 0.0_my_real
            ppmax = max(maxerrind**2,minerrind**2)
            ppssval = ppval/max(maxerrind**2,minerrind**2)
            ppssmin = 0.0_my_real
            ppssmax = 1.0_my_real
         endif
      endif
   end select

case (PREPROC_LOG)
   if (abs(val) < 1.0e-14_my_real) then
      ppval = -14.0_my_real
   else
      ppval = log10(abs(val))
   endif
   select case (selector)
   case (DRAW_NO_FUNCTION)
      ppval = 0.0_my_real
      minval = 1.0_my_real
      maxval = 1.0_my_real
   case (DRAW_SOLUTION, COLOR_SOLUTION)
      minval = minsolut
      maxval = maxsolut
   case (DRAW_TRUE, COLOR_TRUE)
      minval = mintrue
      maxval = maxtrue
   case (DRAW_ERROR, COLOR_ERROR)
      minval = minerror
      maxval = maxerror
   case (DRAW_ERRIND, COLOR_ERRIND)
      minval = minerrind
      maxval = maxerrind
   end select
   if (minval < 0.0_my_real) then
      if (maxval > 0.0_my_real) then
         ppmin = -14
         ppmax = max(log10(abs(minval)),log10(abs(maxval)))
      else
         if (maxval > -1.0e-14_my_real) then
            ppmin = -14
            if (minval > -1.0e-14_my_real) then
               ppmax = -14
            else
               ppmax = log10(abs(minval))
            endif
         else
            ppmin = log10(abs(maxval))
            ppmax = log10(abs(minval))
         endif
      endif
   elseif (minval < 1.0e-14_my_real) then
      ppmin = -14
      if (maxval < 1.0e-14_my_real) then
         ppmax = -14
      else
         ppmax = log10(abs(maxval))
      endif
   else
      ppmin = log10(abs(minval))
      ppmax = log10(abs(maxval))
   endif
   if (ppmax == ppmin) then
      ppssval = ppval
      ppssmin = ppmin
      ppssmax = ppmax
   else
      ppssval = (ppval-ppmin)/(ppmax-ppmin)
      ppssmin = 0.0_my_real
      ppssmax = 1.0_my_real
   endif

end select

ppssval = 2*maxdomain*(ppssval-0.5_my_real)
ppssmin = 2*maxdomain*(ppssmin-0.5_my_real)
ppssmax = 2*maxdomain*(ppssmax-0.5_my_real)

end subroutine preproc_and_scale

!          ---------
subroutine draw_grid
!          ---------

!----------------------------------------------------
! This subroutine draws the grid.  Optionally, the elements and/or vertices
! can be labeled, the associated element can be marked, and the elements,
! vertices, and/or element edges can be colored.
!----------------------------------------------------

!----------------------------------------------------
! Local variables:

real(glfloat) :: fblack(4) = (/glfzero,glfzero,glfzero,glfone/)
logical :: nolight
character(len=14+HOSTLEN) :: title
real(my_real) :: contour_value(num_contour), contour_ps(num_contour), &
                 minz, maxz, pval, psval, psmin, psmax
integer :: i, allocstat
logical(small_logical), allocatable :: labeled(:)

!----------------------------------------------------
! Begin executable code

! reset the graphics window

call reset_view
if (grid_newview) call recalcview

call glutsetwindow(grid_win)
write(title(12:13),'(i2)') my_processor
title(1:11)='PHAML grid '
title(14:14)=' '
title(15:14+HOSTLEN)=host
call glutSetWindowTitle(title)
call gldeletelists(grid_list, 1_glsizei)
call glnewlist(grid_list, gl_compile_and_execute)

! set lighting

if (draw_func == DRAW_NO_FUNCTION .or. draw_func == DRAW_LEVELS .or. &
    draw_func == DRAW_ERRIND .or. color_elements == COLOR_TRANSPARENT .or. &
    color_elements == COLOR_WHITE .or. color_scheme == SCHEME_STRIPE) then
   call gldisable(gl_lighting)
   nolight = .true.
else
   call glenable(gl_lighting)
   nolight = .false.
endif

! amount to shift each partition for exploding

do i=1,nproc
   explode_shift(i)%x = (part_cent(i)%x-lookat_x)*(explode_factor-gldone)
   explode_shift(i)%y = (part_cent(i)%y-lookat_y)*(explode_factor-gldone)
end do

! recursively draw the elements solid

if (color_elements /= COLOR_TRANSPARENT) then

   call glenable(GL_POLYGON_OFFSET_FILL)
   call glpolygonmode(gl_front_and_back, gl_fill)
   call glbegin(gl_triangles)

   call draw_subgrid_solid(REFTREE_ROOT,nolight)

   call glend
   call gldisable(GL_POLYGON_OFFSET_FILL)

endif

! recursively draw the element sides

if (color_lines /= COLOR_TRANSPARENT) then

   call glpolygonmode(gl_front_and_back, gl_line)
   call glbegin(gl_lines)

   call draw_subgrid_lines(REFTREE_ROOT,nolight)

   call glend

endif

! recursively label elements

if (label_elements /= LABEL_NOLABEL .or. label_verts /= LABEL_NOLABEL .or. &
    label_edges /= LABEL_NOLABEL .or. vert_associated_element .or. &
    edge_associated_element) then

   allocate(labeled(size(vertex)),stat=allocstat)
   if (allocstat /= 0) then
      call fatal("allocation failed for labeled",intlist=(/allocstat/))
      stop
   endif

   if (nolight) then
      call glcolor3d(gldzero, gldzero, gldzero)
   else
      call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse,fblack)
   endif

   labeled = .false.

   call draw_subgrid_label(REFTREE_ROOT,labeled,nolight)

   deallocate(labeled,stat=allocstat)
   if (allocstat /= 0) then
      call warning("deallocation failed for labeled",intlist=(/allocstat/))
   endif

endif

! recursively draw contour plot

if (draw_cont /= DRAW_NO_FUNCTION) then

   call glpolygonmode(gl_front_and_back, gl_line)
   call glbegin(gl_lines)

   if (nolight) then
      call glcolor3d(gldzero, gldzero, gldzero)
   else
      call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse,fblack)
   endif

! set the contour values

   if (contour_values_given) then
      contour_value = actual_contours
   else
      call preproc_and_scale(draw_cont,1.0_my_real,pval,psval,minz,psmin,maxz, &
                             psmax)
      do i=1,num_contour
         contour_value(i) = minz+(maxz-minz)*(i-1)/real(num_contour-1,my_real)
         contour_ps(i) = psmin+(psmax-psmin)*(i-1)/real(num_contour-1,my_real)
      end do
   endif

   call draw_contour(REFTREE_ROOT,contour_value,contour_ps)

   call glend

endif

! recursively draw space filling curve

if (draw_sfc_flag) then

   call glpolygonmode(gl_front_and_back, gl_line)
   call glbegin(gl_line_strip)

   if (nolight) then
      call glcolor3d(gldzero, gldzero, gldzero)
   else
      call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse,fblack)
   endif

   call draw_sfc(REFTREE_ROOT)

   call glend
endif

! label in and out vertices

if (draw_sfcio_flag) then

   if (nolight) then
      call glcolor3d(gldzero, gldzero, gldzero)
   else
      call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse,fblack)
   endif

   call draw_sfcio(REFTREE_ROOT,1.0_gldouble)

endif

! draw axes

if (draw_axes_flag) then

   if (nolight) then
      call glcolor3d(gldzero, gldzero, gldzero)
   else
      call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse,fblack)
   endif

   call draw_axes

endif

call glendlist
call glutpostredisplay

return
end subroutine draw_grid

!                    ------------------
recursive subroutine draw_subgrid_lines(elem,nolight)
!                    ------------------

!----------------------------------------------------
! This routine draws the element edges for the subtree below element elem
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: elem
logical, intent(in) :: nolight
!----------------------------------------------------

!----------------------------------------------------
! Local variables:

logical :: part_bound(EDGES_PER_ELEMENT)
integer :: i
integer :: allc(MAX_CHILD), children(MAX_CHILD)
real(my_real) :: x(VERTICES_PER_ELEMENT), y(VERTICES_PER_ELEMENT)
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

! if elem is the root, then do all the initial elements

if (elem == REFTREE_ROOT) then
   i = grid%head_level_elem(1)
   do while (i /= END_OF_LIST)
      call draw_subgrid_lines(i,nolight)
      i = element(i)%next
   end do
   return
endif

! draw this element only if it is a leaf or we are drawing all levels

if (element(elem)%isleaf .or. draw_func == DRAW_LEVELS) then

! draw this element only if one of the vertices is inside the croping range

 if (inside_crop(element(elem)%vertex)) then

! draw the element edges

   if (color_lines /= COLOR_TRANSPARENT) then

! set the vertices into x,y arrays

      x = vertex(element(elem)%vertex(1:VERTICES_PER_ELEMENT))%coord%x
      y = vertex(element(elem)%vertex(1:VERTICES_PER_ELEMENT))%coord%y

! determine if each edge is on a partition boundary

      do i=1,EDGES_PER_ELEMENT
         if (neighbors(i,elem) == BOUNDARY) then
            part_bound(i) = .true.
         elseif (elem_owner(hash_decode_key(neighbors(i,elem),grid%elem_hash)) &
                 /= elem_owner(elem)) then
            part_bound(i) = .true.
         else
            part_bound(i) = .false.
         endif
      end do

! draw the element

      call draw_element_lines(grid,x,y,elem, &
  elem_owner(vertex(element(elem)%vertex(1:VERTICES_PER_ELEMENT))%assoc_elem), &
  elem_owner(edge(element(elem)%edge(1:EDGES_PER_ELEMENT))%assoc_elem), &
                              part_bound,(/.true.,.true.,.true./),nolight,0)

   endif
 endif
endif

! draw the elements in the subtree below the children

allc = ALL_CHILDREN
children = get_child_lid(element(elem)%gid,allc,grid%elem_hash)
do i=1,MAX_CHILD
   if (children(i) /= NO_CHILD) then
      call draw_subgrid_lines(children(i),nolight)
   endif
end do

return
end subroutine draw_subgrid_lines

!                    ------------------
recursive subroutine draw_element_lines(grid,xmr,ymr,elem,vown,eown,part_bound,&
                                        sides,nolight,recur_lev)
!                    ------------------

!----------------------------------------------------
! This routine draws part of the edges of element elem.  If recur_lev
! is smaller than subelement_resolution, it replaces the triangle with
! four triangles and draws the edges of those that are part of elem's boundary.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
real(my_real), intent(in) :: xmr(:),ymr(:)
integer, intent(in) :: elem,vown(:),eown(:)
logical, intent(in) :: part_bound(:),sides(:),nolight
integer, intent(in) :: recur_lev
!----------------------------------------------------
! Local variables:

logical :: csides(EDGES_PER_ELEMENT),cpart_bound(EDGES_PER_ELEMENT)
integer :: i, j, cvown(VERTICES_PER_ELEMENT), ceown(EDGES_PER_ELEMENT)
real(my_real) :: xmid(EDGES_PER_ELEMENT),ymid(EDGES_PER_ELEMENT), &
                 valmr(VERTICES_PER_ELEMENT),val1(1),ppval,ppssval,ppmin, &
                 ppssmin,ppmax,ppssmax
real(gldouble) :: x(VERTICES_PER_ELEMENT),y(VERTICES_PER_ELEMENT), &
                 val(VERTICES_PER_ELEMENT),normal(VERTICES_PER_ELEMENT),val2, &
                 elem_size,max_val
real(gldouble) :: color(3,4), &
                 grey(4)  = (/gldhalf,gldhalf,gldhalf,gldone/), &
                 white(4) = (/gldone, gldone, gldone, gldone/), &
                 black(4) = (/gldzero,gldzero,gldzero,gldone/)
real(glfloat) :: fcolor(3,4)
!----------------------------------------------------
! Begin executable code

if (recur_lev < subelement_resolution) then

! draw the children that have a piece of elem's boundary for better resolution

! define midpoints of the sides

   xmid(1) = (xmr(2)+xmr(3))/2
   ymid(1) = (ymr(2)+ymr(3))/2
   xmid(2) = (xmr(3)+xmr(1))/2
   ymid(2) = (ymr(3)+ymr(1))/2
   xmid(3) = (xmr(1)+xmr(2))/2
   ymid(3) = (ymr(1)+ymr(2))/2

! for each child, determine which sides are on the boundary and the owner of
! the vertices and edges, and draw if any side is on elem's boundary

   csides(1) = .false.
   csides(2) = sides(2)
   csides(3) = sides(3)
   cpart_bound(1) = .false.
   cpart_bound(2) = part_bound(2)
   cpart_bound(3) = part_bound(3)
   cvown = vown(1)
   ceown = eown
   if (any(csides)) then
      call draw_element_lines(grid,(/xmr(1),xmid(3),xmid(2)/), &
                              (/ymr(1),ymid(3),ymid(2)/),elem,cvown,ceown, &
                              cpart_bound,csides,nolight,recur_lev+1)
   endif

   csides(1) = .false.
   csides(2) = sides(3)
   csides(3) = sides(1)
   cpart_bound(1) = .false.
   cpart_bound(2) = part_bound(3)
   cpart_bound(3) = part_bound(1)
   cvown = vown(2)
   ceown(1) = 0
   ceown(2) = eown(3)
   ceown(3) = eown(1)
   if (any(csides)) then
      call draw_element_lines(grid,(/xmr(2),xmid(1),xmid(3)/), &
                              (/ymr(2),ymid(1),ymid(3)/),elem,cvown,ceown, &
                              cpart_bound,csides,nolight,recur_lev+1)
   endif

   csides(1) = .false.
   csides(2) = sides(1)
   csides(3) = sides(2)
   cpart_bound(1) = .false.
   cpart_bound(2) = part_bound(1)
   cpart_bound(3) = part_bound(2)
   cvown = vown(3)
   ceown(1) = 0
   ceown(2) = eown(1)
   ceown(3) = eown(2)
   if (any(csides)) then
      call draw_element_lines(grid,(/xmr(3),xmid(2),xmid(1)/), &
                              (/ymr(3),ymid(2),ymid(1)/),elem,cvown,ceown, &
                              cpart_bound,csides,nolight,recur_lev+1)
   endif

else ! draw the triangle

! (xmr,ymr) are used for evaluating solution; (x,y) are used as vertices

   x = xmr + explode_shift(elem_owner(elem))%x
   y = ymr + explode_shift(elem_owner(elem))%y

! if coloring lines by size, determine the length of the longest side

   if (color_lines == COLOR_SIZE) then
      elem_size = maxval( &
         (/ sqrt((x(1)-x(2))**2 + (y(1)-y(2))**2), &
            sqrt((x(2)-x(3))**2 + (y(2)-y(3))**2), &
            sqrt((x(3)-x(1))**2 + (y(3)-y(1))**2) /) )
      elem_size = elem_size*2**subelement_resolution
      elem_size = log(elem_size)
   endif

! determine the z value for the function being drawn

   select case (draw_func)
      case (DRAW_SOLUTION)
         call graphics_evaluate_soln(xmr,ymr,elem,valmr)
      case (DRAW_TRUE)
         do i=1,VERTICES_PER_ELEMENT
            valmr(i) = graphics_trues(xmr(i),ymr(i))
         end do
      case (DRAW_ERROR)
         call graphics_evaluate_soln(xmr,ymr,elem,valmr)
         do i=1,VERTICES_PER_ELEMENT
            valmr(i) = (valmr(i) - graphics_trues(xmr(i),ymr(i)))
         end do
      case (DRAW_ERRIND)
         if (elem_owner(elem)==my_processor .or.  my_processor==MASTER) then
            valmr = grid%element_errind(elem,1)
         else
            valmr = gldzero
         endif
      case (DRAW_LEVELS)
         valmr = element(elem)%level
      case default
         valmr = gldzero
   end select
   do i=1,VERTICES_PER_ELEMENT
      call preproc_and_scale(draw_func,valmr(i),ppval,ppssval,ppmin,ppssmin, &
                             ppmax,ppssmax)
      val(i) = ppssval
   end do
   normal = normcrossprod(x,y,val)

! find the color of each vertex (in some cases, color of edge)

   do i=1,VERTICES_PER_ELEMENT

      select case (color_lines)
      case (COLOR_VERT_OWNER)
         if (vown(i) /= NO_OWNER) then
            color(i,:) = owner_color(:,vown(i))
         else
            color(i,:) = grey
         endif

      case (COLOR_OWNER)
         if (eown(i) /= NO_OWNER) then
            color(i,:) = owner_color(:,eown(i))
         else
            color(i,:) = grey
         endif

      case (COLOR_BLACK, COLOR_PART_BOUND)
         color(i,:) = black

      case (COLOR_SOLUTION)
         call graphics_evaluate_soln(xmr(i:i),ymr(i:i),elem,val1)
         call preproc_and_scale(color_lines,val1(1),ppval,ppssval,ppmin, &
                                ppssmin,ppmax,ppssmax)
         val2 = ppval
         call get_rainbow(val2,real(ppmin,gldouble),real(ppmax,gldouble), &
                          color(i,:))

      case (COLOR_TRUE)
         val1(1) = graphics_trues(xmr(i),ymr(i))
         call preproc_and_scale(color_lines,val1(1),ppval,ppssval,ppmin, &
                                ppssmin,ppmax,ppssmax)
         val2 = ppval
         call get_rainbow(val2,real(ppmin,gldouble),real(ppmax,gldouble), &
                          color(i,:))

      case (COLOR_ERROR)
         call graphics_evaluate_soln(xmr(i:i),ymr(i:i),elem,val1)
         val1(1) = val1(1) - graphics_trues(xmr(i),ymr(i))
         call preproc_and_scale(color_lines,val1(1),ppval,ppssval,ppmin, &
                                ppssmin,ppmax,ppssmax)
         val2 = ppval
         call get_rainbow(val2,real(ppmin,gldouble),real(ppmax,gldouble), &
                          color(i,:))

      case (COLOR_SIZE)
         call get_rainbow(elem_size,real(minsize,gldouble), &
                          real(maxsize,gldouble),color(i,:))

      case (COLOR_DEGREE)
         if (color_scheme == SCHEME_STEP_SEQ) then
            max_val = step_scheme_steps*step_scheme_hues
         else
            max_val = MAX_DEGREE
         endif
         call get_rainbow(real(edge(element(elem)%edge(i))%degree,gldouble), &
                          1.0_gldouble,max_val,color(i,:))

      end select

   end do

! draw each edge of the triangle

   do i=1,EDGES_PER_ELEMENT

      if (.not. sides(i)) cycle

! add the vertex to the display list

      select case (color_lines)

! when coloring by edge owner or edge degree, draw edge i with color i

      case (COLOR_OWNER, COLOR_DEGREE)

         if (nolight) then
            call glcolor3dv(color(i,:))
         else
            fcolor = color
            call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse,fcolor(i,:))
         endif

         do j=1,3
            if (j==i) cycle
            call glnormal3dv(normal)
            call glvertex3d(x(j),y(j),val(j))
         end do

! for drawing partition boundary, draw opposite line if opposite element is
! boundary or owned by a different processor

      case (COLOR_PART_BOUND)

         if (part_bound(i)) then

            if (nolight) then
               call glcolor3dv(color(i,:))
            else
               fcolor = color
               call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse,fcolor(i,:))
            endif

            do j=1,3
               if (j==i) cycle
               call glnormal3dv(normal)
               call glvertex3d(x(j),y(j),val(j))
            end do

         endif

! for most cases, draw the edge with colors from the vertices

      case default

         do j=1,3
            if (j==i) cycle
            if (nolight) then
               call glcolor3dv(color(j,:))
            else
               fcolor = color
               call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse,fcolor(j,:))
            endif
            call glnormal3dv(normal)
            call glvertex3d(x(j),y(j),val(j))
         end do

      end select

   end do

endif

end subroutine draw_element_lines

!                    ------------------
recursive subroutine draw_subgrid_label(elem,labeled,nolight)
!                    ------------------

!----------------------------------------------------
! This routine labels elements, vertices, etc.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: elem
logical(small_logical), intent(inout) :: labeled(:)
logical, intent(in) :: nolight
!----------------------------------------------------

!----------------------------------------------------
! Local variables:

integer :: j,i,int_key(KEY_SIZE)
integer :: allc(MAX_CHILD), children(MAX_CHILD)
real(gldouble) :: x(VERTICES_PER_ELEMENT), y(VERTICES_PER_ELEMENT)
real(gldouble) :: xx,yy,zz,h
real(my_real) :: ppval,ppssval,ppmin,ppssmin,ppmax,ppssmax
real, parameter :: delta=.010
real(gldouble) :: elem_lc(4) = (/gldzero,gldzero,gldhalf,gldone/), &
                 edge_lc(4) = (/gldzero,gldhalf,gldzero,gldone/), &
                 vert_lc(4) = (/gldhalf,gldzero,gldzero,gldone/), &
                 black(4)   = (/gldzero,gldzero,gldzero,gldone/)
real(glfloat) :: felem_lc(4) = (/glfzero,glfzero,glfhalf,glfone/), &
                 fedge_lc(4) = (/glfzero,glfhalf,glfzero,glfone/), &
                 fvert_lc(4) = (/glfhalf,glfzero,glfzero,glfone/), &
                 fblack(4)   = (/glfzero,glfzero,glfzero,glfone/)
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

! if elem is the root, then do all the initial elements

if (elem == REFTREE_ROOT) then
   i = grid%head_level_elem(1)
   do while (i /= END_OF_LIST)
      call draw_subgrid_label(i,labeled,nolight)
      i = element(i)%next
   end do
   return
endif

! label this element only if it is a leaf or we are drawing all levels

if (element(elem)%isleaf .or. draw_func == DRAW_LEVELS) then

! label this element only if one of the vertices is inside the croping range

 if (inside_crop(element(elem)%vertex)) then

! number the element and vertices

   if (label_elements /= LABEL_NOLABEL .or. label_verts /= LABEL_NOLABEL .or. &
       label_edges /= LABEL_NOLABEL .or. vert_associated_element .or. &
       edge_associated_element) then

! determine the location of the label

      x = vertex(element(elem)%vertex(1:VERTICES_PER_ELEMENT))%coord%x + &
          explode_shift(elem_owner(elem))%x
      y = vertex(element(elem)%vertex(1:VERTICES_PER_ELEMENT))%coord%y + &
          explode_shift(elem_owner(elem))%y
      h = sqrt((x(1)-x(2))**2 + (y(1)-y(2))**2) + &
          sqrt((x(2)-x(3))**2 + (y(2)-y(3))**2) + &
          sqrt((x(3)-x(1))**2 + (y(3)-y(1))**2)
      xx = sum(x)/size(x)
      yy = sum(y)/size(y)
      if (draw_func == DRAW_LEVELS) then
         zz = element(elem)%level
      else
         zz = gldzero
      endif
      call preproc_and_scale(draw_func,zz,ppval,ppssval,ppmin,ppssmin, &
                             ppmax,ppssmax)
      zz = ppssval

! element label

      select case (label_elements)
      case (LABEL_NOLABEL)
      case (LABEL_LID)
         if (nolight) then
            call glcolor3d(elem_lc(1),elem_lc(2),elem_lc(3))
         else
            call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse,felem_lc)
         endif
         call number(elem,xx,yy,zz,h)
      case (LABEL_GID)
         if (nolight) then
            call glcolor3d(elem_lc(1),elem_lc(2),elem_lc(3))
         else
            call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse,felem_lc)
         endif
         call hash_pack_key(element(elem)%gid,int_key,1)
         call number(int_key(1),xx,yy,zz,h)
      end select

! vertex label

      select case (label_verts)
      case (LABEL_NOLABEL)
      case (LABEL_LID)
         if (nolight) then
            call glcolor3d(vert_lc(1),vert_lc(2),vert_lc(3))
         else
            call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse,fvert_lc)
         endif
         do j=1,VERTICES_PER_ELEMENT
           if (.not. labeled(element(elem)%vertex(j))) then
            xx = x(j)
            yy = y(j)
            call number(element(elem)%vertex(j),xx,yy,zz,h)
            labeled(element(elem)%vertex(j)) = .true.
           endif
         end do
      case (LABEL_GID)
         if (nolight) then
            call glcolor3d(vert_lc(1),vert_lc(2),vert_lc(3))
         else
            call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse,fvert_lc)
         endif
         do j=1,VERTICES_PER_ELEMENT
           if (.not. labeled(element(elem)%vertex(j))) then
            xx = x(j)
            yy = y(j)
            call hash_pack_key(vertex(element(elem)%vertex(j))%gid,int_key,1)
            call number(int_key(1),xx,yy,zz,h)
            labeled(element(elem)%vertex(j)) = .true.
           endif
         end do
      end select

! draw a line segment from the vertices into the associated element

      if (vert_associated_element) then
         if (nolight) then
            call glcolor3d(gldzero,gldzero,gldzero)
         else
            call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse,fblack)
         endif
         xx = sum(x)/size(x)
         yy = sum(y)/size(y)
         do j=1,VERTICES_PER_ELEMENT
            if (vertex(element(elem)%vertex(j))%assoc_elem == elem) then
               call gllinewidth(2.0_glfloat)
               call glbegin(gl_lines)
               call glVertex3d(.25_gldouble*xx+.75_gldouble*x(j), &
                               .25_gldouble*yy+.75_gldouble*y(j),zz)
               call glVertex3d(x(j),y(j),zz)
               call glend
               call gllinewidth(glfone)
            endif
         end do
      endif

! draw a line segment from the edges into the associated element

      if (edge_associated_element) then
         if (nolight) then
            call glcolor3d(gldzero,gldzero,gldzero)
         else
            call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse,fblack)
         endif
         do j=1,EDGES_PER_ELEMENT
            if (edge(element(elem)%edge(j))%assoc_elem == elem) then
               xx = (x(mod(j,3)+1) + x(mod(j+1,3)+1))/2
               yy = (y(mod(j,3)+1) + y(mod(j+1,3)+1))/2
               call gllinewidth(2.0_glfloat)
               call glbegin(gl_lines)
               call glVertex3d(.90_gldouble*xx+.10_gldouble*x(j), &
                               .90_gldouble*yy+.10_gldouble*y(j),zz)
               call glVertex3d(xx,yy,zz)
               call glend
               call gllinewidth(glfone)
            endif
         end do
      endif

! edge label

      select case (label_edges)
      case (LABEL_NOLABEL)
      case (LABEL_LID)
         if (nolight) then
            call glcolor3d(edge_lc(1),edge_lc(2),edge_lc(3))
         else
            call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse,fedge_lc)
         endif
         do j=1,EDGES_PER_ELEMENT
            xx = (x(mod(j,3)+1) + x(mod(j+1,3)+1))/2
            yy = (y(mod(j,3)+1) + y(mod(j+1,3)+1))/2
            h = 2.*sqrt( (x(mod(j,3)+1) - x(mod(j+1,3)+1))**2 + &
                         (y(mod(j,3)+1) - y(mod(j+1,3)+1))**2)
            call number(element(elem)%edge(j),xx,yy,zz,h)
         end do
      case (LABEL_GID)
         if (nolight) then
            call glcolor3d(edge_lc(1),edge_lc(2),edge_lc(3))
         else
            call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse,fedge_lc)
         endif
         do j=1,EDGES_PER_ELEMENT
            xx = (x(mod(j,3)+1) + x(mod(j+1,3)+1))/2
            yy = (y(mod(j,3)+1) + y(mod(j+1,3)+1))/2
            h = 2.*sqrt( (x(mod(j,3)+1) - x(mod(j+1,3)+1))**2 + &
                         (y(mod(j,3)+1) - y(mod(j+1,3)+1))**2)
            call hash_pack_key(edge(element(elem)%edge(j))%gid,int_key,1)
            call number(int_key(1),xx,yy,zz,h)
         end do
      end select

   endif

 endif
endif

! draw the elements in the subtree below the children

allc = ALL_CHILDREN
children = get_child_lid(element(elem)%gid,allc,grid%elem_hash)
do i=1,MAX_CHILD
   if (children(i) /= NO_CHILD) then
      call draw_subgrid_label(children(i),labeled,nolight)
   endif
end do

return
end subroutine draw_subgrid_label

!                    ------------------
recursive subroutine draw_subgrid_solid(elem,nolight)
!                    ------------------

!----------------------------------------------------
! This routine draws the solid elements for the subtree below element elem
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: elem
logical, intent(in) :: nolight
!----------------------------------------------------

!----------------------------------------------------
! Local variables:

integer :: i
integer :: allc(MAX_CHILD), children(MAX_CHILD)
real(my_real) :: x(VERTICES_PER_ELEMENT), y(VERTICES_PER_ELEMENT)
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

! if elem is the root, then do all the initial elements

if (elem == REFTREE_ROOT) then
   i = grid%head_level_elem(1)
   do while (i /= END_OF_LIST)
      call draw_subgrid_solid(i,nolight)
      i = element(i)%next
   end do
   return
endif

! only draw leaf elements, unless the function is DRAW_LEVELS

if (element(elem)%isleaf .or. draw_func == DRAW_LEVELS) then

! draw this element only if one of the vertices is inside the croping range

 if (inside_crop(element(elem)%vertex)) then

! set the vertices into x,y arrays

   x = vertex(element(elem)%vertex(1:VERTICES_PER_ELEMENT))%coord%x
   y = vertex(element(elem)%vertex(1:VERTICES_PER_ELEMENT))%coord%y

! draw the element

   call draw_element_solid(grid,x,y,elem,nolight,0)

 endif
endif

! draw the elements in the subtree below the children

allc = ALL_CHILDREN
children = get_child_lid(element(elem)%gid,allc,grid%elem_hash)
do i=1,MAX_CHILD
   if (children(i) /= NO_CHILD) then
      call draw_subgrid_solid(children(i),nolight)
   endif
end do

return
end subroutine draw_subgrid_solid

!                    ------------------
recursive subroutine draw_element_solid(grid,xmr,ymr,elem,nolight,recur_lev)
!                    ------------------

!----------------------------------------------------
! This routine draws a triangle with vertices (xmr,ymr) solid.  If recur_lev
! is smaller than subelement_resolution, it replaces the triangle with
! four triangles by connecting the edge midpoints.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
real(my_real), intent(in) :: xmr(:),ymr(:)
integer, intent(in) :: elem, recur_lev
logical, intent(in) :: nolight
!----------------------------------------------------
! Local variables:

integer :: i
real(my_real) :: xmid(EDGES_PER_ELEMENT),ymid(EDGES_PER_ELEMENT), &
                 valmr(VERTICES_PER_ELEMENT),val1(1),ppval,ppssval,ppmin, &
                 ppssmin,ppmax,ppssmax
real(gldouble) :: x(VERTICES_PER_ELEMENT),y(VERTICES_PER_ELEMENT), &
                 val(VERTICES_PER_ELEMENT),normal(VERTICES_PER_ELEMENT),val2, &
                 elem_size,max_val
real(gldouble) :: color(4), &
                 grey(4)  = (/gldhalf,gldhalf,gldhalf,gldone/), &
                 white(4) = (/gldone, gldone, gldone, gldone/)
real(glfloat) :: fcolor(4), &
                 fgrey(4)  = (/glfhalf,glfhalf,glfhalf,glfone/), &
                 fwhite(4) = (/glfone, glfone, glfone, glfone/)
!----------------------------------------------------
! Begin executable code

if (recur_lev < subelement_resolution) then

! draw four children of this triangle for better resolution

! define midpoints of the sides

   xmid(1) = (xmr(2)+xmr(3))/2
   ymid(1) = (ymr(2)+ymr(3))/2
   xmid(2) = (xmr(3)+xmr(1))/2
   ymid(2) = (ymr(3)+ymr(1))/2
   xmid(3) = (xmr(1)+xmr(2))/2
   ymid(3) = (ymr(1)+ymr(2))/2

! draw the children

   call draw_element_solid(grid, &
                           (/xmr(1),xmid(3),xmid(2)/), &
                           (/ymr(1),ymid(3),ymid(2)/), &
                           elem,nolight,recur_lev+1)
   call draw_element_solid(grid, &
                           (/xmr(2),xmid(1),xmid(3)/), &
                           (/ymr(2),ymid(1),ymid(3)/), &
                           elem,nolight,recur_lev+1)
   call draw_element_solid(grid, &
                           (/xmr(3),xmid(2),xmid(1)/), &
                           (/ymr(3),ymid(2),ymid(1)/), &
                           elem,nolight,recur_lev+1)
   call draw_element_solid(grid, &
                           (/xmid(1),xmid(2),xmid(3)/), &
                           (/ymid(1),ymid(2),ymid(3)/), &
                           elem,nolight,recur_lev+1)

else ! draw the triangle

! (xmr,ymr) are used for evaluating solution; (x,y) are used as vertices

   x = xmr + explode_shift(elem_owner(elem))%x
   y = ymr + explode_shift(elem_owner(elem))%y

! if coloring elements by size, determine the length of the longest side
! of the parent element

   if (color_elements == COLOR_SIZE) then
      elem_size = maxval( &
         (/ sqrt((x(1)-x(2))**2 + (y(1)-y(2))**2), &
            sqrt((x(2)-x(3))**2 + (y(2)-y(3))**2), &
            sqrt((x(3)-x(1))**2 + (y(3)-y(1))**2) /) )
      elem_size = elem_size*2**subelement_resolution
      elem_size = log(elem_size)
   endif

! determine the z value for the function being drawn

   select case (draw_func)
      case (DRAW_SOLUTION)
         call graphics_evaluate_soln(xmr,ymr,elem,valmr)
      case (DRAW_TRUE)
         do i=1,VERTICES_PER_ELEMENT
            valmr(i) = graphics_trues(xmr(i),ymr(i))
         end do
      case (DRAW_ERROR)
         call graphics_evaluate_soln(xmr,ymr,elem,valmr)
         do i=1,VERTICES_PER_ELEMENT
            valmr(i) = valmr(i) - graphics_trues(xmr(i),ymr(i))
         end do
      case (DRAW_ERRIND)
         if (elem_owner(elem) == my_processor .or.  my_processor == MASTER) then
            valmr = grid%element_errind(elem,1)
         else
            valmr = gldzero
         endif
      case (DRAW_LEVELS)
         valmr = element(elem)%level
      case default
         valmr = gldzero
   end select
   do i=1,VERTICES_PER_ELEMENT
      call preproc_and_scale(draw_func,valmr(i),ppval,ppssval,ppmin,ppssmin, &
                             ppmax,ppssmax)
      val(i) = ppssval
   end do
   normal = normcrossprod(xmr,ymr,val)
   call glnormal3dv(normal)

! find the color to draw the element

   do i=1,VERTICES_PER_ELEMENT

      select case (color_elements)
      case (COLOR_OWNER)
         if (elem_owner(elem) /= NO_OWNER) then
            if(nolight) then
               call glcolor4dv(owner_color(:,elem_owner(elem)))
            else
               fcolor = owner_color(:,elem_owner(elem))
               call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse, &
                                 fcolor)
            endif
         else
            if (nolight) then
               call glcolor3d(gldhalf, gldhalf, gldhalf)
            else
               call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse, &
                                    fgrey)
            endif
         endif

      case (COLOR_WHITE)
         if (nolight) then
            call glcolor3d(gldone,gldone,gldone)
         else
            call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse,fwhite)
         endif

      case (COLOR_SOLUTION)
         call graphics_evaluate_soln(xmr(i:i),ymr(i:i),elem,val1)
         call preproc_and_scale(color_elements,val1(1),ppval,ppssval,ppmin, &
                                ppssmin,ppmax,ppssmax)
         val2 = ppval
         call get_rainbow(val2,real(ppmin,gldouble),real(ppmax,gldouble),color)
         if (nolight) then
            call glcolor3dv(color)
         else
            fcolor = color
            call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse, &
                              fcolor)
         endif

      case (COLOR_TRUE)
         val1(1) = graphics_trues(xmr(i),ymr(i))
         call preproc_and_scale(color_elements,val1(1),ppval,ppssval,ppmin, &
                                ppssmin,ppmax,ppssmax)
         val2 = ppval
         call get_rainbow(val2,real(ppmin,gldouble),real(ppmax,gldouble),color)
         if (nolight) then
            call glcolor3dv(color)
         else
            fcolor = color
            call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse, &
                              fcolor)
         endif

      case (COLOR_SIZE)
         call get_rainbow(elem_size,real(minsize,gldouble), &
                          real(maxsize,gldouble),color)
         if (nolight) then
            call glcolor3dv(color)
         else
            fcolor = color
            call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse, &
                              fcolor)
         endif

      case (COLOR_DEGREE)
         if (color_scheme == SCHEME_STEP_SEQ) then
            max_val = step_scheme_steps*step_scheme_hues
         else
            max_val = MAX_DEGREE
         endif
         call get_rainbow(real(element(elem)%degree,gldouble),1.0_gldouble, &
                          max_val,color)
         if (nolight) then
            call glcolor3dv(color)
         else
            fcolor = color
            call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse, &
                              fcolor)
         endif

      case (COLOR_ERROR)
         call graphics_evaluate_soln(xmr(i:i),ymr(i:i),elem,val1)
         val1(1) = val1(1) - graphics_trues(xmr(i),ymr(i))
         call preproc_and_scale(color_elements,val1(1),ppval,ppssval,ppmin, &
                                ppssmin,ppmax,ppssmax)
         val2 = ppval
         call get_rainbow(val2,real(ppmin,gldouble),real(ppmax,gldouble),color)
         if (nolight) then
            call glcolor3dv(color)
         else
            fcolor = color
            call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse, &
                              fcolor)
         endif

      case (COLOR_ERRIND)
         if (elem_owner(elem) /= NO_OWNER) then
            val1(1) = grid%element_errind(elem,1)
            call preproc_and_scale(color_elements,val1(1),ppval,ppssval,ppmin, &
                                   ppssmin,ppmax,ppssmax)
            val2 = ppval
            call get_rainbow(val2,real(ppmin,gldouble),real(ppmax,gldouble),color)
            if (nolight) then
               call glcolor3dv(color)
            else
               fcolor = color
               call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse, &
                                 fcolor)
            endif
         else
            if (nolight) then
               call glcolor3d(gldhalf, gldhalf, gldhalf)
            else
               call glmaterialfv(gl_front_and_back,gl_ambient_and_diffuse,fgrey)
            endif
         endif

      end select

! add this vertex to the drawing of the triangle

      call glvertex3d(x(i),y(i),val(i))

   end do ! next vertex

endif ! check for recursion

end subroutine draw_element_solid

!                    ------------
recursive subroutine draw_contour(elem,contour_value,contour_ps)
!                    ------------

!----------------------------------------------------
! This routine draws the contour plot
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: elem
real(my_real) :: contour_value(:),contour_ps(:)
!----------------------------------------------------

!----------------------------------------------------
! Local variables:

real(my_real) :: x(3), y(3)
integer :: i
integer :: allc(MAX_CHILD), children(MAX_CHILD)
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

! if elem is the root, then do all the initial elements

if (elem == REFTREE_ROOT) then
   i = grid%head_level_elem(1)
   do while (i /= END_OF_LIST)
      call draw_contour(i,contour_value,contour_ps)
      i = element(i)%next
   end do
   return
endif

! draw contour lines in this element only if it is a leaf

if (element(elem)%isleaf) then

! draw contour lines in this element only if one of the vertices is inside the croping range

   if (inside_crop(element(elem)%vertex)) then

! set the vertices

      x = vertex(element(elem)%vertex(1:VERTICES_PER_ELEMENT))%coord%x
      y = vertex(element(elem)%vertex(1:VERTICES_PER_ELEMENT))%coord%y

! draw the element

      call draw_element_contour(grid,x,y,contour_value,contour_ps,elem,0)

   endif ! not cropped
endif ! is a leaf

! draw the contours in the subtree below the children

allc = ALL_CHILDREN
children = get_child_lid(element(elem)%gid,allc,grid%elem_hash)
do i=1,MAX_CHILD
   if (children(i) /= NO_CHILD) then
      call draw_contour(children(i),contour_value,contour_ps)
   endif
end do

return
end subroutine draw_contour

!                    --------------------
recursive subroutine draw_element_contour(grid,xmr,ymr,contour_value, &
                                          contour_ps,elem,recur_lev)
!                    --------------------

!----------------------------------------------------
! This routine draws the contours in a triangle with vertices (xmr,ymr).
! If recur_lev is smaller than subelement_resolution, it replaces the
! triangle with four triangles by connecting the edge midpoints.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
real(my_real), intent(in) :: xmr(:),ymr(:), contour_value(:), contour_ps(:)
integer, intent(in) :: elem, recur_lev
!----------------------------------------------------
! Local variables:

integer :: i
real(my_real) :: xmid(EDGES_PER_ELEMENT),ymid(EDGES_PER_ELEMENT), &
                 zmr(VERTICES_PER_ELEMENT)
real(gldouble) :: x(VERTICES_PER_ELEMENT),y(VERTICES_PER_ELEMENT), &
                 z(VERTICES_PER_ELEMENT)
real(gldouble) :: xt, yt, zt, frac, xcross1, xcross2, ycross1, ycross2, zcross
real(my_real) :: val,ppval,ppssval,ppmin,ppssmin,ppmax,ppssmax
real(my_real) :: delta=.05_my_real
!----------------------------------------------------
! Begin executable code

if (recur_lev < subelement_resolution) then

! draw four children of this triangle for better resolution

! define midpoints of the sides

   xmid(1) = (xmr(2)+xmr(3))/2
   ymid(1) = (ymr(2)+ymr(3))/2
   xmid(2) = (xmr(3)+xmr(1))/2
   ymid(2) = (ymr(3)+ymr(1))/2
   xmid(3) = (xmr(1)+xmr(2))/2
   ymid(3) = (ymr(1)+ymr(2))/2

! draw the children

   call draw_element_contour(grid, &
                           (/xmr(1),xmid(3),xmid(2)/), &
                           (/ymr(1),ymid(3),ymid(2)/), &
                           contour_value,contour_ps,elem,recur_lev+1)
   call draw_element_contour(grid, &
                           (/xmr(2),xmid(1),xmid(3)/), &
                           (/ymr(2),ymid(1),ymid(3)/), &
                           contour_value,contour_ps,elem,recur_lev+1)
   call draw_element_contour(grid, &
                           (/xmr(3),xmid(2),xmid(1)/), &
                           (/ymr(3),ymid(2),ymid(1)/), &
                           contour_value,contour_ps,elem,recur_lev+1)
   call draw_element_contour(grid, &
                           (/xmid(1),xmid(2),xmid(3)/), &
                           (/ymid(1),ymid(2),ymid(3)/), &
                           contour_value,contour_ps,elem,recur_lev+1)

else ! draw the triangle

! (xmr,ymr) are used for evaluating solution; (x,y) are used as vertices

   x = xmr + explode_shift(elem_owner(elem))%x
   y = ymr + explode_shift(elem_owner(elem))%y

   select case (draw_cont)
   case (DRAW_SOLUTION)
      call graphics_evaluate_soln(xmr,ymr,elem,zmr)
      z = zmr
   case (DRAW_TRUE)
      do i=1,VERTICES_PER_ELEMENT
         z(i) = graphics_trues(xmr(i),ymr(i))
      end do
   case (DRAW_ERROR)
      call graphics_evaluate_soln(xmr,ymr,elem,zmr)
      do i=1,VERTICES_PER_ELEMENT
         z(i) = (zmr(i) - graphics_trues(xmr(i),ymr(i)))
      end do
   case default
      z(i) = gldzero
   end select
   do i=1,VERTICES_PER_ELEMENT
      val = z(i)
      call preproc_and_scale(draw_cont,val,ppval,ppssval,ppmin,ppssmin,ppmax, &
                             ppssmax)
      z(i) = ppval
   end do

! order the vertices by z value

   xt = x(1); yt = y(1); zt = z(1)
   if (z(2) < z(1)) then
      xt = x(1); yt = y(1); zt = z(1)
      if (z(3) < z(1)) then
         if (z(2) < z(3)) then
            x(1) = x(2); y(1) = y(2); z(1) = z(2)
            x(2) = x(3); y(2) = y(3); z(2) = z(3)
         else
            x(1) = x(3); y(1) = y(3); z(1) = z(3)
         endif
         x(3) = xt; y(3) = yt; z(3) = zt
      else
         x(1) = x(2); y(1) = y(2); z(1) = z(2)
         x(2) = xt; y(2) = yt; z(2) = zt
      endif
   elseif (z(3) < z(1)) then
      x(1) = x(3); y(1) = y(3); z(1) = z(3)
      x(3) = x(2); y(3) = y(2); z(3) = z(2)
      x(2) = xt; y(2) = yt; z(2) = zt
   elseif (z(3) < z(2)) then
      xt = x(2); yt = y(2); zt = z(2)
      x(2) = x(3); y(2) = y(3); z(2) = z(3)
      x(3) = xt; y(3) = yt; z(3) = zt
   endif

! if z(1)==z(3), the function is flat in the element and has no contours

   if (z(1) /= z(3)) then

! for each contour value

      do i = 1,num_contour

! see if it passes through this element

         if (contour_value(i) < z(1)) cycle
         if (contour_value(i) > z(3)) exit

! see where it crosses the 1-3 edge

         frac = (contour_value(i)-z(1))/(z(3)-z(1))
         xcross1 = (1.0_gldouble - frac)*x(1) + frac*x(3)
         ycross1 = (1.0_gldouble - frac)*y(1) + frac*y(3)

! see where it crosses one of the other edges

         if (contour_value(i) == z(2)) then
            xcross2 = x(2)
            ycross2 = y(2)
         elseif (contour_value(i) < z(2)) then
            frac = (contour_value(i)-z(1))/(z(2)-z(1))
            xcross2 = (1.0_gldouble - frac)*x(1) + frac*x(2)
            ycross2 = (1.0_gldouble - frac)*y(1) + frac*y(2)
         else
            frac = (contour_value(i)-z(2))/(z(3)-z(2))
            xcross2 = (1.0_gldouble - frac)*x(2) + frac*x(3)
            ycross2 = (1.0_gldouble - frac)*y(2) + frac*y(3)
         endif

! add the line segment to the display list

         select case(contour_location)
         case (CONT_XYPLANE)
            val = 0
            call preproc_and_scale(draw_func,val,ppval,ppssval,ppmin,ppssmin, &
                                   ppmax,ppssmax)
            if (draw_func == DRAW_NO_FUNCTION) then
               zcross = ppssval
            else
               zcross = ppssmin-delta*(ppssmax-ppssmin)
            endif
         case (CONT_SURFACE)
             zcross = contour_ps(i)
         end select
         call glvertex3d(xcross1,ycross1,zcross)
         call glvertex3d(xcross2,ycross2,zcross)

      end do ! next contour
   endif ! not flat element

endif ! recursive call

end subroutine draw_element_contour

!                    --------
recursive subroutine draw_sfc(elem)
!                    --------

!----------------------------------------------------
! This routine draws the space filling curve used for partitioning
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: elem
!----------------------------------------------------

!----------------------------------------------------
! Local variables:

real(gldouble) :: x, y, z
integer :: i
integer :: allc(MAX_CHILD), children(MAX_CHILD)
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

! if elem is the root, then do all the initial elements

if (elem == REFTREE_ROOT) then
   i = grid%head_level_elem(1)
   do while (i /= END_OF_LIST)
      call draw_sfc(i)
      i = element(i)%next
   end do
   return
endif

! add the center of this element to the SFC only if it is a leaf

if (element(elem)%isleaf) then

! set the midpoint

   x = sum(vertex(element(elem)%vertex(1:VERTICES_PER_ELEMENT))%coord%x)/VERTICES_PER_ELEMENT + &
       explode_shift(elem_owner(elem))%x
   y = sum(vertex(element(elem)%vertex(1:VERTICES_PER_ELEMENT))%coord%y)/VERTICES_PER_ELEMENT + &
       explode_shift(elem_owner(elem))%y
   z = gldzero

! add the midpoint to the display list

   call glvertex3d(x,y,z)

endif

! draw the sfc in the subtree below the children

allc = ALL_CHILDREN
children = get_child_lid(element(elem)%gid,allc,grid%elem_hash)
do i=1,MAX_CHILD
   if (children(element(elem)%order(i)) /= NO_CHILD) then
      call draw_sfc(children(element(elem)%order(i)))
   endif
end do

return
end subroutine draw_sfc

!                    ----------
recursive subroutine draw_sfcio(elem,h)
!                    ----------

!----------------------------------------------------
! This routine labels the in and out vertices for the space filling curve
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: elem
real(gldouble), intent(in) :: h
!----------------------------------------------------

!----------------------------------------------------
! Local variables:

real(gldouble) :: x, y, z, hh
real(gldouble), parameter :: shiftx = .1_gldouble, shifty = .01_gldouble
integer :: i
integer :: allc(MAX_CHILD), children(MAX_CHILD)
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

! if elem is the root, then do all the initial elements

if (elem == REFTREE_ROOT) then
   i = grid%head_level_elem(1)
   do while (i /= END_OF_LIST)
      hh = 8*elem_diam(i)
      call draw_sfcio(i,hh)
      i = element(i)%next
   end do
   return
endif

! add the labels to this element only if it is a leaf

if (element(elem)%isleaf) then

! set the midpoint

   x = sum(vertex(element(elem)%vertex(1:VERTICES_PER_ELEMENT))%coord%x)/VERTICES_PER_ELEMENT + &
       explode_shift(elem_owner(elem))%x
   y = sum(vertex(element(elem)%vertex(1:VERTICES_PER_ELEMENT))%coord%y)/VERTICES_PER_ELEMENT + &
       explode_shift(elem_owner(elem))%y
   z = gldzero

! label the in vertex

   call number(1,(x+vertex(element(elem)%in)%coord%x)/2 - h*shiftx, &
                 (y+vertex(element(elem)%in)%coord%y)/2 - h*shifty, &
                 z,h)

! label the out vertex

   call number(0,(x+vertex(element(elem)%out)%coord%x)/2 - h*shiftx, &
                 (y+vertex(element(elem)%out)%coord%y)/2 - h*shifty, &
                 z,h)

endif

! label the elements in the subtree below the children

allc = ALL_CHILDREN
children = get_child_lid(element(elem)%gid,allc,grid%elem_hash)
do i=1,MAX_CHILD
   if (children(element(elem)%order(i)) /= NO_CHILD) then
     call draw_sfcio(children(element(elem)%order(i)),h/sqrt(2.0_gldouble))
   endif
end do

return
end subroutine draw_sfcio

!          ----------------------
subroutine graphics_evaluate_soln(x,y,elem,z)
!          ----------------------

!----------------------------------------------------
! This routine calls evaluate solution with the current eigen and
! compnt taking care when compnt is the L1 or L2 sum
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x(:), y(:)
integer, intent(in) :: elem
real(my_real), intent(out) :: z(:)
!----------------------------------------------------
! Local variables:

integer :: i
real(my_real) :: zloc(ncompnt,1,size(z)), z1(1,1,size(z))
!----------------------------------------------------
! Begin executable code

select case (compnt)

case(-1)
   call evaluate_soln_local(grid,x,y,elem,(/(i,i=1,ncompnt)/),(/eigen/),zloc)
   z = 0.0_my_real
   do i=1,ncompnt
      z = z + abs(zloc(i,1,:))
   end do

case(-2)
   call evaluate_soln_local(grid,x,y,elem,(/(i,i=1,ncompnt)/),(/eigen/),zloc)
   z = 0.0_my_real
   do i=1,ncompnt
      z = z + zloc(i,1,:)**2
   end do

case default
   call evaluate_soln_local(grid,x,y,elem,(/compnt/),(/eigen/),z1)
   z = z1(1,1,:)

end select

end subroutine graphics_evaluate_soln

!        --------------
function graphics_trues(x,y)
!        --------------

!----------------------------------------------------
! This routine calls trues for the current eigen and compnt taking
! care when compnt is the L1 or L2 sum
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x,y
real(my_real) :: graphics_trues
!----------------------------------------------------
! Local variables:

integer :: i
real(my_real) :: zloc(ncompnt)
!----------------------------------------------------
! Begin executable code

select case (compnt)

case(-1)
   do i=1,ncompnt
      zloc(i) = trues(x,y,i,eigen)
   end do
   graphics_trues = sum(abs(zloc))

case(-2)
   do i=1,ncompnt
      zloc(i) = trues(x,y,i,eigen)
   end do
   graphics_trues = sum(zloc**2)

case default
! TEMP071219 for exact solution for the battery problem
   if (TEMP_battery) then
      call evaluate_oldsoln_local(x,y,zloc(1),comp=compnt,eigen=eigen)
      graphics_trues = zloc(1)
   else
      graphics_trues = trues(x,y,compnt,eigen)
   endif

end select

end function graphics_trues

!        ---------
function elem_diam(elem)
!        ---------

!----------------------------------------------------
! This routine computes the diameter of an element, defined to be
! the minimum distance from the center to one of the vertices.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: elem
real(gldouble) :: elem_diam
!----------------------------------------------------
! Local variables:

real(gldouble) :: x, y, mindist
integer :: i
!----------------------------------------------------
! Begin executable code

! set the midpoint

x = sum(vertex(element(elem)%vertex(1:VERTICES_PER_ELEMENT))%coord%x)/VERTICES_PER_ELEMENT
y = sum(vertex(element(elem)%vertex(1:VERTICES_PER_ELEMENT))%coord%y)/VERTICES_PER_ELEMENT

! compute the distance to each vertex, keeping the minimum

mindist = huge(0.0_gldouble)

do i=1,VERTICES_PER_ELEMENT
   mindist = min(mindist,sqrt((vertex(element(elem)%vertex(i))%coord%x-x)**2 + &
                              (vertex(element(elem)%vertex(i))%coord%y-y)**2))
end do

elem_diam = mindist

end function elem_diam

!          ---------
subroutine get_rainbow(val,minval,maxval,c)
!          ---------

!----------------------------------------------------
! This routine sets the color for rainbow plots
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments:

real(gldouble), intent(in) :: val,maxval,minval
real(gldouble), intent(out) :: c(4)

!----------------------------------------------------
! Local variables:

real(gldouble) :: f,h,s,v,sbase
integer :: hbase

!----------------------------------------------------
! Begin executable code

if (maxval > minval) then
   f = (val-minval)/(maxval-minval)
else ! probably maxval==minval
   f = gldhalf
endif
if (f < 0.0) f = 0
if (f > 1.0) f = 1

if (color_scheme == SCHEME_STRIPE) then
   f = real(floor(f*(num_contour-1)),gldouble)/(num_contour-1)
endif

if (color_scheme == SCHEME_STEP_SEQ .and. f==1.0_gldouble) f=0.9999_gldouble

select case(color_scheme)

case (SCHEME_RAINBOW,SCHEME_STRIPE)
   if (f < .25) then
      c(1) = 0.0_gldouble
      c(2) = 4.0_gldouble * f
      c(3) = 1.0_gldouble
      c(4) = 1.0_gldouble
   elseif (f < .5) then
      c(1) = 0.0_gldouble
      c(2) = 1.0_gldouble
      c(3) = 2.0_gldouble - 4.0_gldouble*f
      c(4) = 1.0_gldouble
   elseif (f < .75) then
      c(1) = 4.0_gldouble * f - 2.0_gldouble
      c(2) = 1.0_gldouble
      c(3) = 0.0_gldouble
      c(4) = 1.0_gldouble
   else
      c(1) = 1.0_gldouble
      c(2) = 4.0_gldouble - 4.0_gldouble*f
      c(3) = 0.0_gldouble
      c(4) = 1.0_gldouble
   endif

case (SCHEME_DOUBLE_RAINBOW)
   f = 2*f
   if (f < .25) then
      c(1) = 0.0_gldouble
      c(2) = 4.0_gldouble * f
      c(3) = 1.0_gldouble
      c(4) = 1.0_gldouble
   elseif (f < .5) then
      c(1) = 0.0_gldouble
      c(2) = 1.0_gldouble
      c(3) = 2.0_gldouble - 4.0_gldouble*f
      c(4) = 1.0_gldouble
   elseif (f < .75) then
      c(1) = 4.0_gldouble * f - 2.0_gldouble
      c(2) = 1.0_gldouble
      c(3) = 0.0_gldouble
      c(4) = 1.0_gldouble
   elseif (f < 1.0) then
      c(1) = 1.0_gldouble
      c(2) = 4.0_gldouble - 4.0_gldouble*f
      c(3) = 0.0_gldouble
      c(4) = 1.0_gldouble
   elseif (f < 1.25) then
      c(1) = 0.5_gldouble
      c(2) = 0.5_gldouble + 2.0_gldouble * (f-1)
      c(3) = 1.0_gldouble
      c(4) = 1.0_gldouble
   elseif (f < 1.5) then
      c(1) = 0.5_gldouble
      c(2) = 1.0_gldouble
      c(3) = 0.5_gldouble + (2.0_gldouble - 4.0_gldouble*(f-1))/2
      c(4) = 1.0_gldouble
   elseif (f < 1.75) then
      c(1) = .5_gldouble+(4.0_gldouble * (f-1) - 2.0_gldouble)/2
      c(2) = 1.0_gldouble
      c(3) = 0.5_gldouble
      c(4) = 1.0_gldouble
   else
      c(1) = 1.0_gldouble
      c(2) = .2_gldouble+(4.0_gldouble - 4.0_gldouble*(f-1))/2
      c(3) = 0.5_gldouble
      c(4) = 1.0_gldouble
   endif

case (SCHEME_GRAY)
   c = f
   c(4) = 1.0_gldouble

case (SCHEME_STEP_SEQ)
   c(4) = 1.0_gldouble
   hbase = floor(f*step_scheme_hues)
   sbase = f*step_scheme_hues - hbase
   h = (360.0_gldouble*hbase)/step_scheme_hues
   s = 0.2_gldouble + 0.8_gldouble*(1-sbase)
   v = 0.8_gldouble + 0.2_gldouble*sbase
   call hsv2rgb(h,s,v,c(1),c(2),c(3))

! stepped sequential scheme from
! http://geography.uoregon.edu/datagraphics/color_scales.htm
!   c(4) = 1.0_gldouble
!   select case (floor(25*f))
!   case (0)
!      c(1) = 0.600_gldouble; c(2) = 0.060_gldouble; c(3) = 0.060_gldouble
!   case (1)
!      c(1) = 0.700_gldouble; c(2) = 0.175_gldouble; c(3) = 0.175_gldouble
!   case (2)
!      c(1) = 0.800_gldouble; c(2) = 0.320_gldouble; c(3) = 0.320_gldouble
!   case (3)
!      c(1) = 0.900_gldouble; c(2) = 0.495_gldouble; c(3) = 0.495_gldouble
!   case (4)
!      c(1) = 1.000_gldouble; c(2) = 0.700_gldouble; c(3) = 0.700_gldouble
!   case (5)
!      c(1) = 0.600_gldouble; c(2) = 0.330_gldouble; c(3) = 0.060_gldouble
!   case (6)
!      c(1) = 0.700_gldouble; c(2) = 0.438_gldouble; c(3) = 0.175_gldouble
!   case (7)
!      c(1) = 0.800_gldouble; c(2) = 0.560_gldouble; c(3) = 0.320_gldouble
!   case (8)
!      c(1) = 0.900_gldouble; c(2) = 0.697_gldouble; c(3) = 0.495_gldouble
!   case (9)
!      c(1) = 1.000_gldouble; c(2) = 0.850_gldouble; c(3) = 0.700_gldouble
!   case (10)
!      c(1) = 0.420_gldouble; c(2) = 0.600_gldouble; c(3) = 0.060_gldouble
!   case (11)
!      c(1) = 0.525_gldouble; c(2) = 0.700_gldouble; c(3) = 0.175_gldouble
!   case (12)
!      c(1) = 0.640_gldouble; c(2) = 0.800_gldouble; c(3) = 0.320_gldouble
!   case (13)
!      c(1) = 0.765_gldouble; c(2) = 0.900_gldouble; c(3) = 0.495_gldouble
!   case (14)
!      c(1) = 0.900_gldouble; c(2) = 1.000_gldouble; c(3) = 0.700_gldouble
!   case (15)
!      c(1) = 0.060_gldouble; c(2) = 0.420_gldouble; c(3) = 0.600_gldouble
!   case (16)
!      c(1) = 0.175_gldouble; c(2) = 0.525_gldouble; c(3) = 0.700_gldouble
!   case (17)
!      c(1) = 0.320_gldouble; c(2) = 0.640_gldouble; c(3) = 0.800_gldouble
!   case (18)
!      c(1) = 0.495_gldouble; c(2) = 0.765_gldouble; c(3) = 0.900_gldouble
!   case (19)
!      c(1) = 0.700_gldouble; c(2) = 0.900_gldouble; c(3) = 1.000_gldouble
!   case (20)
!      c(1) = 0.150_gldouble; c(2) = 0.060_gldouble; c(3) = 0.600_gldouble
!   case (21)
!      c(1) = 0.262_gldouble; c(2) = 0.175_gldouble; c(3) = 0.700_gldouble
!   case (22)
!      c(1) = 0.400_gldouble; c(2) = 0.320_gldouble; c(3) = 0.800_gldouble
!   case (23)
!      c(1) = 0.562_gldouble; c(2) = 0.495_gldouble; c(3) = 0.900_gldouble
!   case (24)
!      c(1) = 0.750_gldouble; c(2) = 0.700_gldouble; c(3) = 1.000_gldouble
!   case default
!      c(1) = 0.000_gldouble; c(2) = 0.000_gldouble; c(3) = 0.000_gldouble
!   end select

end select

return
end subroutine get_rainbow

!          -------
subroutine hsv2rgb(h,s,v,r,g,b)
!          -------

!----------------------------------------------------
! This routine converts (h,s,v) color representation to (r,g,b).
! The formulas come from the Wikipedia page on HSL and HSV
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(gldouble), intent(in) :: h,s,v
real(gldouble), intent(out) :: r,g,b
!----------------------------------------------------
! Local variables:

real(gldouble) :: f, p, q, t
!----------------------------------------------------
! Begin executable code

f = h/60 - floor(h/60)
p = v*(1-s)
q = v*(1-f*s)
t = v*(1-(1-f)*s)
select case (mod(floor(h/60),6))
case (0)
   r=v; g=t; b=p
case (1)
   r=q; g=v; b=p
case (2)
   r=p; g=v; b=t
case (3)
   r=p; g=q; b=v
case (4)
   r=t; g=p; b=v
case (5)
   r=v; g=p; b=q
case default
   print *,"bad case in hsv2rgb"
   r=0; g=0; b=0
end select

end subroutine hsv2rgb

!        -------------
function normcrossprod(x,y,z)
!        -------------

!----------------------------------------------------
! The function returns the normalized cross product of the vectors
! (p1,p2) and (p1,p3) with the coordinates of the p's in x,y,z.  If the
! points are counterclockwise, the cross product is negated to point in
! the direction of a clockwise orientation of the points.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(gldouble), dimension(:) :: x,y,z
real(gldouble), dimension(3) :: normcrossprod

!----------------------------------------------------
! Local variables

real(gldouble) :: t1(3),t2(3),norm

!----------------------------------------------------
! Begin executable code

t1(1) = x(2) - x(1)
t1(2) = y(2) - y(1)
t1(3) = z(2) - z(1)
t2(1) = x(3) - x(1)
t2(2) = y(3) - y(1)
t2(3) = z(3) - z(1)

normcrossprod(1) = t1(2)*t2(3) - t1(3)*t2(2)
normcrossprod(2) = t1(3)*t2(1) - t1(1)*t2(3)
normcrossprod(3) = t1(1)*t2(2) - t1(2)*t2(1)

if (normcrossprod(3) > 0.0_gldouble) normcrossprod = -normcrossprod

norm = sqrt(dot_product(normcrossprod,normcrossprod))
if (norm /= gldzero) normcrossprod = normcrossprod/norm

return
end function normcrossprod

!                -----------
logical function inside_crop(vert)
!                -----------

!----------------------------------------------------
! This routine determines if at least one of the vertices is inside
! the croping range
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: vert(:)

!----------------------------------------------------
! Local variables:
integer :: i
real(my_real) :: x,y

!----------------------------------------------------
! Begin executable code

inside_crop = .false.
do i=1,VERTICES_PER_ELEMENT
   x = vertex(vert(i))%coord%x
   y = vertex(vert(i))%coord%y
   if (x>xcrop1 .and. x<xcrop2 .and. y>ycrop1 .and. y<ycrop2) &
      inside_crop = .true.
end do

return
end function inside_crop

!          ---------
subroutine draw_axes
!          ---------

!----------------------------------------------------
! This routine draws x, y and z axes
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

!----------------------------------------------------
! Local variables:

real(gldouble), parameter :: delta = .05_gldouble, half=.5_gldouble
real(gldouble) :: xstart,xend,ystart,yend,zstart,zend
real(gldouble) :: xmind,xmaxd,ymind,ymaxd,zmind,zmaxd,zminlabel,zmaxlabel
real(gldouble) :: sdelta,temp1,temp2
!----------------------------------------------------
! Begin executable code

xmind = xmin
xmaxd = xmax
ymind = ymin
ymaxd = ymax

sdelta = delta*min(xmax-xmin,ymax-ymin)

! set the endpoints for the axes

xstart = xmin - delta*(xmax-xmin)
xend   = xmax + delta*(xmax-xmin)
ystart = ymin - delta*(ymax-ymin)
yend   = ymax + delta*(ymax-ymin)
call preproc_and_scale(draw_func,1.0_my_real,temp1,temp2,zminlabel,zmind, &
                       zmaxlabel,zmaxd)
zstart = zmind-delta*(zmaxd-zmind)
zend   = zmaxd+delta*(zmaxd-zmind)

! draw the axes

call glbegin(gl_lines)
call glvertex3d(xstart,ystart,zstart) ! x axis
call glvertex3d(xend,ystart,zstart)
call glvertex3d(xstart,ystart,zstart) ! y axis
call glvertex3d(xstart,yend,zstart)
call glvertex3d(xstart,ystart,zstart) ! z axis
call glvertex3d(xstart,ystart,zend)

! add ticks

call glvertex3d(xmind,ystart-half*sdelta,zstart) ! x ticks
call glvertex3d(xmind,ystart+half*sdelta,zstart)
call glvertex3d(half*(xmind+xmaxd),ystart-half*sdelta,zstart)
call glvertex3d(half*(xmind+xmaxd),ystart+half*sdelta,zstart)
call glvertex3d(xmaxd,ystart-half*sdelta,zstart)
call glvertex3d(xmaxd,ystart+half*sdelta,zstart)

call glvertex3d(xstart-half*sdelta,ymind,zstart) ! y ticks
call glvertex3d(xstart+half*sdelta,ymind,zstart)
call glvertex3d(xstart-half*sdelta,half*(ymind+ymaxd),zstart)
call glvertex3d(xstart+half*sdelta,half*(ymind+ymaxd),zstart)
call glvertex3d(xstart-half*sdelta,ymaxd,zstart)
call glvertex3d(xstart+half*sdelta,ymaxd,zstart)

if (zmax-zmin > 10._my_real*tiny(myzero)) then
   call glvertex3d(xstart,ystart-half*sdelta,zmind) ! z ticks
   call glvertex3d(xstart,ystart+half*sdelta,zmind)
   call glvertex3d(xstart,ystart-half*sdelta,half*(zmind+zmaxd))
   call glvertex3d(xstart,ystart+half*sdelta,half*(zmind+zmaxd))
   call glvertex3d(xstart,ystart-half*sdelta,zmaxd)
   call glvertex3d(xstart,ystart+half*sdelta,zmaxd)
endif

call glend

! label the axes

call real_number(xmin,xmind,ystart-sdelta,zstart,10.*sdelta)
call real_number((xmin+xmax)/2.,half*(xmind+xmaxd),ystart-sdelta, &
                 zstart,10.*sdelta)
call real_number(xmax,xmaxd,ystart-sdelta,zstart,10.*sdelta)

call real_number(ymin,xstart-sdelta,ymind,zstart,10.*sdelta)
call real_number((ymin+ymax)/2.,xstart-sdelta,half*(ymind+ymaxd), &
                 zstart,10.*sdelta)
call real_number(ymax,xstart-sdelta,ymaxd,zstart,10.*sdelta)

if (zmax-zmin > 10._my_real*tiny(0._my_real)) then
   call real_number(zminlabel,xstart,ystart-sdelta,zmind,10.*sdelta)
   call real_number((zminlabel+zmaxlabel)/2.,xstart,ystart-sdelta, &
                    half*(zmind+zmaxd),10.*sdelta)
   call real_number(zmaxlabel,xstart,ystart-sdelta,zmaxd,10.*sdelta)
endif

return
end subroutine draw_axes

!          ------
subroutine number(n,x,y,z,h)
!          ------

!----------------------------------------------------
! This routine draws the number n at the point (x,y,z) using h as a
! scaling factor for the height of the number
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: n
real(gldouble), intent(in) :: x,y,z,h

!----------------------------------------------------
! Local variables

character(len=9) :: chn
integer :: i

!----------------------------------------------------
! Begin executable code

write(chn,'(i9)') n
if (chn /= " ") then
   do while (chn(1:1) == " ")
      chn = chn(2:9)
   end do
endif
call glmatrixmode(gl_modelview)
call glpushmatrix
call gltranslated(x,y,z)
call glscaled(h/2500.0_gldouble,h/2500.0_gldouble,1._gldouble)
do i=1,9
   call glutstrokecharacter(glut_stroke_roman,ichar(chn(i:i)))
end do
call glpopmatrix

return
end subroutine number

!          -----------
subroutine real_number(n,x,y,z,h)
!          -----------

!----------------------------------------------------
! This routine draws the real number n at the point (x,y,z) using h as a
! scaling factor for the height of the number
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: n
real(gldouble), intent(in) :: x,y,z,h

!----------------------------------------------------
! Local variables

character(len=8) chn
integer :: i

!----------------------------------------------------
! Begin executable code

write(chn,'(g8.2)') n
call glmatrixmode(gl_modelview)
call glpushmatrix
call gltranslated(x-h/8.,y,z)
call glscaled(h/2500.0_gldouble,h/2500.0_gldouble,1._gldouble)
do i=1,8
   call glutstrokecharacter(glut_stroke_roman,ichar(chn(i:i)))
end do
call glpopmatrix

return
end subroutine real_number

!          ------------
subroutine grid_display
!          ------------
if (.not. window_initialized) return
call glutsetwindow(grid_win)
call reset_view
if (grid_newview) call recalcview
if (explode_factor /= old_explode_factor) then
   call draw_grid
   old_explode_factor = explode_factor
endif
call glclear(ior(gl_color_buffer_bit,gl_depth_buffer_bit))
call glcalllist(grid_list)
call glutswapbuffers
end subroutine grid_display

!          ----------
subroutine recalcview
!          ----------
real, parameter :: pi = 3.1415926

select case(draw_func)
   case(DRAW_NO_FUNCTION)
      zmin = myzero
      zmax = myzero
   case(DRAW_SOLUTION)
      zmin = minsolut
      zmax = maxsolut
   case(DRAW_TRUE)
      zmin = mintrue
      zmax = maxtrue
   case(DRAW_ERROR)
      zmin = minerror
      zmax = maxerror
   case(DRAW_ERRIND)
      zmin = minerrind
      zmax = maxerrind
   case(DRAW_LEVELS)
      zmin = myzero
      zmax = grid%nlev
   case default
      zmin = myzero
      zmax = myzero
end select
grid_newview = .false.
end subroutine recalcview

!          --------------
subroutine set_draw_lines(selection)
!          --------------
integer(glcint), intent(in out) :: selection
color_lines = selection
call draw_grid
return
end subroutine set_draw_lines

!          -------------
subroutine set_draw_elem(selection)
!          -------------
integer(glcint), intent(in out) :: selection
color_elements = selection
call draw_grid
return
end subroutine set_draw_elem

!          --------
subroutine set_func(selection)
!          --------
integer(glcint), intent(in out) :: selection
draw_func = selection
grid_newview = .true.
call draw_grid
return
end subroutine set_func

!          ---------------
subroutine set_subelem_res(selection)
!          ---------------
integer(glcint), intent(in out) :: selection
select case (selection)
case (-2)
   if (subelement_resolution > 0) then
      subelement_resolution = subelement_resolution - 1
   endif
case (-1)
   subelement_resolution = subelement_resolution + 1
case default
   subelement_resolution = selection
end select
grid_newview = .true.
call draw_grid
return
end subroutine set_subelem_res

!          ----------------
subroutine set_color_scheme(selection)
!          ----------------
integer(glcint), intent(in out) :: selection
color_scheme = selection
call set_owner_color
call draw_grid
return
end subroutine set_color_scheme

!          ---------------
subroutine set_step_scheme(selection)
!          ---------------
integer(glcint), intent(in out) :: selection
integer(glcint) :: hold
character(len=32) :: str
select case(selection)
case(1)
   step_scheme_steps=5
case(2)
   step_scheme_steps=10
case(3)
   step_scheme_steps=16
case(4)
   step_scheme_steps=max(1,step_scheme_steps-1)
case(5)
   step_scheme_steps=step_scheme_steps+1
case(6)
   step_scheme_hues=5
case(7)
   step_scheme_hues=10
case(8)
   step_scheme_hues=16
case(9)
   step_scheme_hues=max(1,step_scheme_hues-1)
case(10)
   step_scheme_hues=step_scheme_hues+1
end select
hold = glutgetmenu()
call glutsetmenu(scheme_menu)
write(str,"(A18,2I4)") "stepped sequential",step_scheme_steps,step_scheme_hues
call glutchangetomenuentry(5,trim(str),SCHEME_STEP_SEQ)
call glutsetmenu(hold)
call set_owner_color
call draw_grid
end subroutine set_step_scheme

!          -----------
subroutine set_contour(selection)
!          -----------
integer(glcint), intent(in out) :: selection
draw_cont = selection
grid_newview = .true.
call draw_grid
end subroutine set_contour

!          -----------------
subroutine set_contour_param(selection)
!          -----------------
integer(glcint), intent(in out) :: selection
integer :: allocstat

select case(selection)
   case(1)
      print *,"enter number of uniform contour lines"
      read *,num_contour
      contour_values_given = .false.
   case(2)
      print *,"enter number of contours:"
      read *, num_contour
      if (allocated(actual_contours)) deallocate(actual_contours,stat=allocstat)
      allocate(actual_contours(num_contour),stat=allocstat)
      if (allocstat /= 0) then
         call fatal("allocation failed for actual_contours", &
                    intlist=(/allocstat/))
         stop
      endif
      print *,"enter ",num_contour," contour values:"
      read *,actual_contours
      contour_values_given = .true.
   case(3)
      num_contour = num_contour + 1
      contour_values_given = .false.
   case(4)
      num_contour = num_contour - 1
      if (num_contour < 2) num_contour = 2
      contour_values_given = .false.
   case(5)
      num_contour = num_contour + 10
      contour_values_given = .false.
   case(6)
      num_contour = num_contour - 10
      if (num_contour < 2) num_contour = 2
      contour_values_given = .false.
   case(7)
      num_contour = num_contour*2 - 1
      contour_values_given = .false.
   case(8)
      num_contour = (num_contour+1)/2
      if (num_contour < 2) num_contour = 2
      contour_values_given = .false.
end select

grid_newview = .true.
call draw_grid
end subroutine set_contour_param

!          --------------------
subroutine set_contour_location(selection)
!          --------------------
integer(glcint), intent(inout) :: selection
contour_location = selection
grid_newview = .true.
call draw_grid
end subroutine set_contour_location

!          ----------------
subroutine set_preproc_func(selection)
!          ----------------
integer(glcint), intent(in out) :: selection
preproc_func = selection
grid_newview = .true.
call draw_grid
return
end subroutine set_preproc_func

!          -------
subroutine set_sfc(selection)
!          -------
integer(glcint), intent(in out) :: selection

select case(selection)
   case(0)
      draw_sfc_flag = .not. draw_sfc_flag
   case(1)
      draw_sfcio_flag = .not. draw_sfcio_flag
end select

grid_newview = .true.
call draw_grid
end subroutine set_sfc

!          ----------
subroutine set_offset(selection)
!          ----------
integer(glcint), intent(in out) :: selection

offset = offset + selection
if (offset < 0.0_glfloat) offset = 0.0_glfloat
call glpolygonoffset(offset,0.0_glfloat)

grid_newview = .true.
call draw_grid
end subroutine set_offset

!          -------------
subroutine toggle_lights(selection)
!          -------------
integer(glcint), intent(in out) :: selection

select case(selection)
   case(0)
      if (movelighton) then
         call gldisable(gl_light0)
         movelighton = .false.
      else
         call glenable(gl_light0)
         movelighton = .true.
      endif
   case(1)
      if (rightlighton) then
         call gldisable(gl_light1)
         rightlighton = .false.
      else
         call glenable(gl_light1)
         rightlighton = .true.
      endif
   case(2)
      if (leftlighton) then
         call gldisable(gl_light2)
         leftlighton = .false.
      else
         call glenable(gl_light2)
         leftlighton = .true.
      endif
   case(3)
      if (toplighton) then
         call gldisable(gl_light3)
         toplighton = .false.
      else
         call glenable(gl_light3)
         toplighton = .true.
      endif
   case(4)
      if (bottomlighton) then
         call gldisable(gl_light4)
         bottomlighton = .false.
      else
         call glenable(gl_light4)
         bottomlighton = .true.
      endif
end select

return
end subroutine toggle_lights

!          --------------
subroutine set_elem_label(selection)
!          --------------
integer(glcint), intent(in out) :: selection
label_elements = selection
call draw_grid
return
end subroutine set_elem_label

!          -------------
subroutine set_vert_label(selection)
!          -------------
integer(glcint), intent(in out) :: selection
label_verts = selection
call draw_grid
return
end subroutine set_vert_label

!          -------------
subroutine set_edge_label(selection)
!          -------------
integer(glcint), intent(in out) :: selection
label_edges = selection
call draw_grid
return
end subroutine set_edge_label

!          ----------------
subroutine set_compnt_scale(selection)
!          ----------------
integer(glcint), intent(in out) :: selection
select case(selection)
case(0)
   indiv_compnt_scale = .true.
case(1)
   indiv_compnt_scale = .false.
end select
call set_max_min
call draw_grid
return
end subroutine set_compnt_scale

!          -------------
subroutine set_assoc_elem(selection)
!          -------------
integer(glcint), intent(in out) :: selection
select case(selection)
case(1) ! toggle vertex associated element
   vert_associated_element = .not. vert_associated_element
case(2) ! toggle edge associated element
   edge_associated_element = .not. edge_associated_element
end select
call draw_grid
return
end subroutine set_assoc_elem

!          ------------
subroutine select_eigen(selection)
!          ------------
integer(glcint), intent(in out) :: selection
eigen = selection
call set_max_min
call draw_grid
return
end subroutine select_eigen

!          -------------
subroutine select_compnt(selection)
!          -------------
integer(glcint), intent(in out) :: selection
compnt = selection
call set_max_min
call draw_grid
return
end subroutine select_compnt

!          ----------------
subroutine write_postscript(selection)
!          ----------------
use rendereps
integer(glcint), intent(in out) :: selection
integer :: bsize
character(len=11) :: filename

! just a colored triangle plus triangle edges seems to require about
! 69*number of leaf elements for the buffer size.  100*total elements
! should be plenty.  If it is not, then the output postscript file
! contains just a white square.  Note that with subelement resolution, the
! number of elements drawn is 4**subelement_resolution times number of elements.

filename = "render .eps"
write(filename(7:7),"(i1)") my_processor
bsize = 100*size(element)*4**subelement_resolution
call glutsetcursor(glut_cursor_wait)
select case(selection)
case(0)
   call outputeps(bsize,.true.,filename)
case(1)
   call outputeps(bsize,.false.,filename)
end select
call glutsetcursor(glut_cursor_inherit)

return
end subroutine write_postscript

!          ------------
subroutine menu_handler(selection)
!          ------------
integer(glcint), intent(in out) :: selection

select case(selection)
case(12)
   draw_axes_flag = .not. draw_axes_flag
   grid_newview = .true.
   call draw_grid
case(13)
   print *,"enter crop region as xmin, xmax, ymin, ymax"
   read *,xcrop1, xcrop2, ycrop1, ycrop2
   grid_newview = .true.
   call draw_grid

end select

return
end subroutine menu_handler

!          ------------
subroutine init_windows
!          ------------
real(gldouble) :: lookfrom_x, lookfrom_y, lookfrom_z
integer(glcint) :: lines, inter, func, lights, elem_label, vert_label, &
                   edge_label, mark_av, contour, contour_ncont, contour_vals, &
                   contour_place, postscript, sfc, off, preproc, subelem, &
                   step_scheme_params, component_scale
integer(glcint) :: menuid
character(len=32) :: str

! set the grid window view and callback functions

   call glutsetwindow(grid_win)
   call glpolygonoffset(offset,0.0_glfloat)
   lookat_x = (xmax+xmin)/2
   lookat_y = (ymax+ymin)/2
   lookat_z = 0
   lookfrom_x = lookat_x + 3.0_gldouble*max(xmax-xmin,ymax-ymin)
   lookfrom_y = lookat_y - 9.0_gldouble*max(xmax-xmin,ymax-ymin)
   lookfrom_z = lookat_z + 4.0_gldouble*max(xmax-xmin,ymax-ymin)
   menuid = view_modifier_init(lookfrom_x, lookfrom_y, lookfrom_z, &
                               lookat_x,   lookat_y,   lookat_z)

! create menu for grid window

   lines = glutcreatemenu(set_draw_lines)
   call glutaddmenuentry("no lines",COLOR_TRANSPARENT)
   call glutaddmenuentry("black",COLOR_BLACK)
   call glutaddmenuentry("edge owner",COLOR_OWNER)
   call glutaddmenuentry("vertex owner",COLOR_VERT_OWNER)
   call glutaddmenuentry("computed solution",COLOR_SOLUTION)
   call glutaddmenuentry("true solution",COLOR_TRUE)
   call glutaddmenuentry("error",COLOR_ERROR)
   call glutaddmenuentry("size",COLOR_SIZE)
   call glutaddmenuentry("degree",COLOR_DEGREE)
   call glutaddmenuentry("partition boundary black",COLOR_PART_BOUND)
   inter = glutcreatemenu(set_draw_elem)
   call glutaddmenuentry("transparent",COLOR_TRANSPARENT)
   call glutaddmenuentry("white",COLOR_WHITE)
   call glutaddmenuentry("owner",COLOR_OWNER)
   call glutaddmenuentry("computed solution",COLOR_SOLUTION)
   call glutaddmenuentry("true solution",COLOR_TRUE)
   call glutaddmenuentry("error",COLOR_ERROR)
   call glutaddmenuentry("size",COLOR_SIZE)
   call glutaddmenuentry("degree",COLOR_DEGREE)
   call glutaddmenuentry("error indicator",COLOR_ERRIND)
   func = glutcreatemenu(set_func)
   call glutaddmenuentry("no function",DRAW_NO_FUNCTION)
   call glutaddmenuentry("computed solution",DRAW_SOLUTION)
   call glutaddmenuentry("true solution",DRAW_TRUE)
   call glutaddmenuentry("error",DRAW_ERROR)
   call glutaddmenuentry("levels",DRAW_LEVELS)
   call glutaddmenuentry("error indicator",DRAW_ERRIND)
   preproc = glutcreatemenu(set_preproc_func)
   call glutaddmenuentry("none",PREPROC_NONE)
   call glutaddmenuentry("-f",PREPROC_NEG)
   call glutaddmenuentry("abs(f)",PREPROC_ABS)
   call glutaddmenuentry("f**2",PREPROC_SQ)
   call glutaddmenuentry("log(abs(f))",PREPROC_LOG)
   lights = glutcreatemenu(toggle_lights)
   call glutaddmenuentry("movable light",0_glcint)
   call glutaddmenuentry("right light",1_glcint)
   call glutaddmenuentry("left light",2_glcint)
   call glutaddmenuentry("top light",3_glcint)
   call glutaddmenuentry("bottom light",4_glcint)
   elem_label = glutcreatemenu(set_elem_label)
   call glutaddmenuentry("none",LABEL_NOLABEL)
   call glutaddmenuentry("local id",LABEL_LID)
   call glutaddmenuentry("global id",LABEL_GID)
   vert_label = glutcreatemenu(set_vert_label)
   call glutaddmenuentry("none",LABEL_NOLABEL)
   call glutaddmenuentry("local id",LABEL_LID)
   call glutaddmenuentry("global id",LABEL_GID)
   edge_label = glutcreatemenu(set_edge_label)
   call glutaddmenuentry("none",LABEL_NOLABEL)
   call glutaddmenuentry("local id",LABEL_LID)
   call glutaddmenuentry("global id",LABEL_GID)
   subelem = glutcreatemenu(set_subelem_res)
   call glutaddmenuentry("0",0_glcint)
   call glutaddmenuentry("1",1_glcint)
   call glutaddmenuentry("2",2_glcint)
   call glutaddmenuentry("3",3_glcint)
   call glutaddmenuentry("increase",-1_glcint)
   call glutaddmenuentry("decrease",-2_glcint)
   step_scheme_params = glutcreatemenu(set_step_scheme)
   call glutaddmenuentry("steps=5",1_glcint)
   call glutaddmenuentry("steps=10",2_glcint)
   call glutaddmenuentry("steps=16",3_glcint)
   call glutaddmenuentry("decrement steps",4_glcint)
   call glutaddmenuentry("increment steps",5_glcint)
   call glutaddmenuentry("hues=5",6_glcint)
   call glutaddmenuentry("hues=10",7_glcint)
   call glutaddmenuentry("hues=16",8_glcint)
   call glutaddmenuentry("decrement hues",9_glcint)
   call glutaddmenuentry("increment hues",10_glcint)
   scheme_menu = glutcreatemenu(set_color_scheme)
   call glutaddmenuentry("rainbow",SCHEME_RAINBOW)
   call glutaddmenuentry("double rainbow",SCHEME_DOUBLE_RAINBOW)
   call glutaddmenuentry("gray scale",SCHEME_GRAY)
   call glutaddmenuentry("striped",SCHEME_STRIPE)
! if this is changed so that stepped sequential is not the 5th entry, need to also
! make a change in set_step_scheme
   write(str,"(A18,2I4)") "stepped sequential",step_scheme_steps,step_scheme_hues
   call glutaddmenuentry(trim(str),SCHEME_STEP_SEQ)
   call glutaddsubmenu("change stepped sequential",step_scheme_params)
   mark_av = glutcreatemenu(set_assoc_elem)
   call glutaddmenuentry("toggle vertex associated element",1_glcint)
   call glutaddmenuentry("toggle edge associated element",2_glcint)
   contour_ncont = glutcreatemenu(set_contour_param)
   call glutaddmenuentry("increment by 1",3_glcint)
   call glutaddmenuentry("decrement by 1",4_glcint)
   call glutaddmenuentry("increment by 10",5_glcint)
   call glutaddmenuentry("decrement by 10",6_glcint)
   call glutaddmenuentry("double",7_glcint)
   call glutaddmenuentry("cut in half",8_glcint)
   call glutaddmenuentry("enter number in debug window",1_glcint)
   contour_vals = glutcreatemenu(set_contour_param)
   call glutaddmenuentry("enter values in debug window",2_glcint)
   contour_place = glutcreatemenu(set_contour_location)
   call glutaddmenuentry("x-y plane",CONT_XYPLANE)
   call glutaddmenuentry("surface",CONT_SURFACE)
   contour = glutcreatemenu(set_contour)
   call glutaddmenuentry("no contour plot",DRAW_NO_FUNCTION)
   call glutaddmenuentry("computed solution",DRAW_SOLUTION)
   call glutaddmenuentry("true solution",DRAW_TRUE)
   call glutaddmenuentry("error",DRAW_ERROR)
   call glutaddsubmenu("set number of uniform lines",contour_ncont)
   call glutaddsubmenu("set nonuniform lines",contour_vals)
   call glutaddsubmenu("contour location",contour_place)
   postscript = glutcreatemenu(write_postscript)
   call glutaddmenuentry("Write out Encapsulated PS (sorted)",0_glcint)
   call glutaddmenuentry("Write out Encapsulated PS (unsorted)",1_glcint)
   sfc = glutcreatemenu(set_sfc)
   call glutaddmenuentry("toggle space filling curve",0_glcint)
   call glutaddmenuentry("toggle in/out vertex label",1_glcint)
   off = glutcreatemenu(set_offset)
   call glutaddmenuentry("decrease by 10",-10_glcint)
   call glutaddmenuentry("decrease by 1",-1_glcint)
   call glutaddmenuentry("increase by 1",1_glcint)
   call glutaddmenuentry("increase by 10",10_glcint)
   eigen_menu = glutcreatemenu(select_eigen)
   menu_neigen = 1
   write(str,"(a14,i3)") "eigenfunction ",1
   call glutaddmenuentry(trim(str),1_glcint)
   component_menu = glutcreatemenu(select_compnt)
   menu_ncompnt = 3
   write(str,"(a6)") "L1 sum"
   call glutaddmenuentry(trim(str),-1_glcint)
   write(str,"(a6)") "L2 sum"
   call glutaddmenuentry(trim(str),-2_glcint)
   write(str,"(a10,i3)") "component ",1
   call glutaddmenuentry(trim(str),1_glcint)
   component_scale = glutcreatemenu(set_compnt_scale)
   call glutaddmenuentry("individual",0_glcint)
   call glutaddmenuentry("all the same",1_glcint)
   grid_menu = glutcreatemenu(menu_handler)
   call glutaddsubmenu("view modifier",menuid)
   call glutaddsubmenu("element edge color",lines)
   call glutaddsubmenu("element interior color",inter)
   call glutaddsubmenu("function",func)
   call glutaddsubmenu("contour plots",contour)
   call glutaddsubmenu("preprocess function",preproc)
   call glutaddsubmenu("subelement resolution",subelem)
   call glutaddsubmenu("color scheme",scheme_menu)
   call glutaddsubmenu("toggle lights",lights)
   call glutaddsubmenu("element label",elem_label)
   call glutaddsubmenu("edge label", edge_label)
   call glutaddsubmenu("vertex label", vert_label)
   call glutaddsubmenu("associated element", mark_av)
   call glutaddsubmenu("eigenfunction to use",eigen_menu)
   call glutaddsubmenu("component to use",component_menu)
   call glutaddsubmenu("component scale",component_scale)
   call glutaddsubmenu("space filling curve",sfc)
   call glutaddsubmenu("grid offset",off)
   call glutaddmenuentry("crop (debug window)",13_glcint)
   call glutaddmenuentry("toggle axes",12_glcint)
   call glutaddsubmenu("write postscript",postscript)
   call glutattachmenu(glut_right_button)

! set up lighting conditions for grid

   call glclearcolor(1.0_glclampf, 1.0_glclampf, 1.0_glclampf, 1.0_glclampf)
   call gllightfv(gl_light0, gl_diffuse, (/glfone,glfone,glfone,glfone/))
   call gllightfv(gl_light1, gl_diffuse, (/glfone,glfone,glfone,glfone/))
   call gllightfv(gl_light2, gl_diffuse, (/glfone,glfone,glfone,glfone/))
   call gllightfv(gl_light3, gl_diffuse, (/glfone,glfone,glfone,glfone/))
   call gllightfv(gl_light4, gl_diffuse, (/glfone,glfone,glfone,glfone/))
!   call gllightfv(gl_light1, gl_position, (/-2.5, -.5, -2.0, 0.0/))
   call gllightfv(gl_light1, gl_position, (/-7.1268,1.68,-6.81, 0.0/))
   call gllightfv(gl_light2, gl_position, (/1.5,-.5,-2.,0./))
   call gllightfv(gl_light3, gl_position, (/-.5,-.5,-2.,0./))
   call gllightfv(gl_light4, gl_position, (/-.5,-.5, 2.,0./))
   call glenable(gl_lighting)
   if (movelighton) call glenable(gl_light0)
   if (rightlighton) call glenable(gl_light1)
   if (leftlighton) call glenable(gl_light2)
   if (toplighton) call glenable(gl_light3)
   if (bottomlighton) call glenable(gl_light4)
   call gllightmodelfv(gl_light_model_ambient, (/glfhalf,glfhalf,glfhalf,glfone/))
   call gldepthfunc(gl_lequal)
   call glenable(gl_depth_test)

   window_initialized = .true.

end subroutine init_windows

! some compiler doesn't let timer pass itself, so use 2 timers
! that pass each other

!                    -----
recursive subroutine timer(selection)
!                    -----
integer(glcint), intent(in out) :: selection
integer :: i

i = selection ! just to shut up compilers that warn about unused arguments

call process_message
call gluttimerfunc(10_glcuint,timer2,0_glcint)

return
end subroutine timer

!                    -----
recursive subroutine timer2(selection)
!                    -----
integer(glcint), intent(in out) :: selection
integer :: i

i = selection ! just to shut up compilers that warn about unused arguments

call process_message
call gluttimerfunc(10_glcuint,timer,0_glcint)

return
end subroutine timer2

! TEMP071219 for exact solution for the battery problem, copy exact to oldsoln

!          ----------------
subroutine exact_to_oldsoln
!          ----------------

!----------------------------------------------------
! Local variables:

integer :: lev, elem, i, edge
!----------------------------------------------------
! Begin executable code

! initialize data structures for oldsoln; this puts solution in there

allocate(grid%initial_neighbor(3,size(grid%element)))
grid%initial_neighbor = BOUNDARY
call set_grid_for_old_soln(grid)
call copy_old(grid)

! copy exact into oldsoln

do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (grid%element(elem)%oldleaf) then
         if (associated(grid%element(elem)%exact)) then
            grid%element(elem)%oldsoln = grid%element(elem)%exact
         endif
         do i=1,EDGES_PER_ELEMENT
            edge = grid%element(elem)%edge(i)
            if (associated(grid%edge(edge)%exact)) then
               grid%edge(edge)%oldsoln = grid%edge(edge)%exact
            endif
         end do ! next edge
      endif
      elem = grid%element(elem)%next
   end do ! next element
end do ! next level
grid%vertex_oldsoln = grid%vertex_exact

end subroutine exact_to_oldsoln

! end TEMP071219

end module graphics_mod

!-------------------------------------------
! The main program
!-------------------------------------------

!          --------------
subroutine phaml_graphics
!          --------------
use graphics_mod
use opengl_gl
use opengl_glut
implicit none

integer :: spawn_form=0, allocstat, system_size, eq_type
character(len=32) :: dummy_char
logical :: junk1,junk2,update_umod
type(phaml_solution_type) :: phaml_solution
integer :: proc,ni,nr
integer, pointer :: recv_int(:)
real(my_real), pointer :: recv_real(:)
!----------------------------------------------------
! Begin executable code

! nullify grid data structures pointers so they can be tested for associated

nullify (grid%element,grid%edge,grid%vertex)

! initialize communication

junk1 = .false.
junk2 = .false.
allocate(procs,stat=allocstat)
if (allocstat /= 0) then
   call fatal("allocation failed for procs",intlist=(/allocstat/))
   stop
endif
call init_comm(procs,spawn_form,junk1,junk2,dummy_char,outunit, &
               errunit,system_size,eq_type,my_pde_id,nproc,grid%max_blen, &
               grid%triangle_files,update_umod)
if (update_umod) then
   if (PARALLEL /= SEQUENTIAL) then
      call phaml_recv(procs,proc,recv_int,ni,recv_real,nr,101)
      if (recv_int(1) /= 9) then
         call fatal("In graphics start up, received code that is not update_usermod.")
         stop
      endif
      phaml_solution%procs = procs
      if (associated(recv_int)) deallocate(recv_int)
      if (associated(recv_real)) deallocate(recv_real)
      call update_usermod(phaml_solution)
   endif
endif

! use unit 6 for all output from the graphics process

outunit = 6
errunit = 6

! initialize GLUT

call glutinit
call glutinitdisplaymode(ior(ior(glut_double,glut_rgb),glut_depth))
call glutinitwindowsize(500_glcint,500_glcint)
call gluttimerfunc(10_glcuint,timer,0_glcint)

! create the window

grid_win = glutcreatewindow("PHAML grid")
call glutdisplayfunc(grid_display)

! go into the GLUT infinite loop

call glutmainloop

end subroutine phaml_graphics
