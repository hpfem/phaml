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

module evaluate

!----------------------------------------------------
! This module contains routines for evaluation of the solution.
!
! communication tags in this module are of the form 9xx
!----------------------------------------------------

!----------------------------------------------------
! Other modules used are:

use global
use message_passing
use gridtype_mod
use basis_functions
!----------------------------------------------------

implicit none
private
public evaluate_soln, evaluate_soln_slave, evaluate_soln_local, &
       evaluate_oldsoln_local, copy_old, set_grid_for_old_soln, &
       find_containing_leaf

!----------------------------------------------------
! The following variables are defined:

type(grid_type), pointer, save :: grid_for_old_soln
integer, save :: roundoff_fudge = 1000

!----------------------------------------------------

contains

!          -------------
subroutine evaluate_soln(procs,x,y,u,ux,uy,uxx,uyy,comp,eigen)
!          -------------

!----------------------------------------------------
! This routine evaluates the solution at the points in the arrays (x,y)
! and returns them in the array soln.  The return value is 0.0 for
! points that are outside the domain.  This would only be called by
! the master or a slave requesting an evaluation by another slave.
! If comp and/or eigen is present, it tells which component and/or which
! eigenfunction to return.  The default is 1 for each.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(proc_info), intent(in) :: procs
real(my_real), intent(in) :: x(:),y(:)
real(my_real), optional, intent(out) :: u(:),ux(:),uy(:),uxx(:),uyy(:)
integer, optional, intent(in) :: comp,eigen

!----------------------------------------------------
! Local variables:

integer, allocatable :: send_int(:)
real(my_real), allocatable :: send_real(:)
integer, pointer :: recv_int(:)
real(my_real), pointer :: recv_real(:)
integer :: i, j, ni, nr, proc, first_int, first_real, allocstat

!----------------------------------------------------
! Begin executable code

! Invoke the slaves to evaluate the solution at the points in
! elements they own, and merge the results into a single array

! send the message to evaluate the solution and the points at which to evaluate

call pack_procs_size(this_processors_procs,ni,nr)
allocate(send_int(8+ni),send_real(nr+2*size(x)),stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in evaluate_soln",procs=procs)
   return
endif
send_int(1) = 3 ! code for "evaluate"
if (present(comp)) then
   send_int(2) = comp
else
   send_int(2) = 1
endif
if (present(eigen)) then
   send_int(3) = eigen
else
   send_int(3) = 1
endif
send_int(4:8) = 0
if (present(u  )) send_int(4) = 1
if (present(ux )) send_int(5) = 1
if (present(uy )) send_int(6) = 1
if (present(uxx)) send_int(7) = 1
if (present(uyy)) send_int(8) = 1
ni = ni+8
first_int = 9
first_real = 1
call pack_procs(this_processors_procs,send_int,first_int,send_real,first_real)
send_real(nr+1:nr+size(x)) = x
send_real(nr+size(x)+1:nr+2*size(x)) = y
nr = nr+2*size(x)
do proc=1,num_proc(procs)
   call phaml_send(procs,proc,send_int,ni,send_real,nr,101)
end do
deallocate(send_real,stat=allocstat)

! receive the solution and flags indicating which points were evaluated
! from each slave

if (present(u  )) u   = 0.0_my_real
if (present(ux )) ux  = 0.0_my_real
if (present(uy )) uy  = 0.0_my_real
if (present(uxx)) uxx = 0.0_my_real
if (present(uyy)) uyy = 0.0_my_real

do i=1,num_proc(procs)
   call phaml_recv(procs,proc,recv_int,ni,recv_real,nr,910)
   j = 0
   if (present(u  )) then
      where (recv_int == 1) u   = recv_real(j+1:j+size(x))
      j = j+size(x)
   endif
   if (present(ux )) then
      where (recv_int == 1) ux  = recv_real(j+1:j+size(x))
      j = j+size(x)
   endif
   if (present(uy )) then
      where (recv_int == 1) uy  = recv_real(j+1:j+size(x))
      j = j+size(x)
   endif
   if (present(uxx)) then
      where (recv_int == 1) uxx = recv_real(j+1:j+size(x))
      j = j+size(x)
   endif
   if (present(uyy)) then
      where (recv_int == 1) uyy = recv_real(j+1:j+size(x))
      j = j+size(x)
   endif
   deallocate(recv_int,recv_real,stat=allocstat)
end do

end subroutine evaluate_soln

!          -------------------
subroutine evaluate_soln_slave(grid,invoker_procs,x,y,comp,eigen, &
                               u,ux,uy,uxx,uyy,which)
!          -------------------

!----------------------------------------------------
! This is the slaves version of evaluate_soln.  If the u's are present, the
! solution is returned in them.  Otherwise, it sends the solutions indicated
! by which to the invoker.  Use of u's or which should be exclusive.
! comp and eigen tell which component and eigenfunction to return.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
type(proc_info), intent(in) :: invoker_procs
real(my_real), intent(in) :: x(:),y(:)
integer, intent(in) :: comp, eigen
real(my_real), intent(out), optional :: u(:),ux(:),uy(:),uxx(:),uyy(:)
integer, optional, intent(in) :: which(:)

!----------------------------------------------------
! Local variables:

real(my_real) :: loc_u(size(x)), loc_ux(size(x)), loc_uy(size(x)), &
                 loc_uxx(size(x)),loc_uyy(size(x)),xc(3), yc(3)
integer :: have_it(size(x))
integer :: i,j,k,elem,nbasis,astat,isub,nr
real(my_real), allocatable :: basis(:),basisx(:),basisy(:),basisxx(:),basisyy(:)
real(my_real), allocatable :: send_real(:)
integer :: loc_which(5)

!----------------------------------------------------
! Begin executable code

! determine which solutions to evaluate

if (present(which)) then
   loc_which = which
else
   loc_which = 0
   if (present(u  )) loc_which(1) = 1
   if (present(ux )) loc_which(2) = 1
   if (present(uy )) loc_which(3) = 1
   if (present(uxx)) loc_which(4) = 1
   if (present(uyy)) loc_which(5) = 1
endif

have_it = 0

! look for points in my partition

! For each point ...

do i=1,size(x)

! find the element it is in.  have_it==1 if I own the element.
! elem is 0 if I do not own it.

   call find_containing_element(x(i),y(i),have_it(i),elem,grid)

   if (have_it(i)==1) then

! evaluate the basis functions at this point

      nbasis = VERTICES_PER_ELEMENT
      do j=1,EDGES_PER_ELEMENT
         if (grid%edge(grid%element(elem)%edge(j))%degree > 1) then
            nbasis = nbasis + grid%edge(grid%element(elem)%edge(j))%degree - 1
         endif
      end do
      if (grid%element(elem)%degree > 2) then
         nbasis = nbasis + &
               ((grid%element(elem)%degree-2)*(grid%element(elem)%degree-1))/2
      endif
      allocate(basis(nbasis),basisx(nbasis),basisy(nbasis),basisxx(nbasis), &
               basisyy(nbasis),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("memory allocation failed in evaluate_soln_slave")
         stop
      endif
      xc = grid%vertex(grid%element(elem)%vertex)%coord%x
      yc = grid%vertex(grid%element(elem)%vertex)%coord%y
      if (loc_which(4)==1 .or. loc_which(5)==1) then
         call p_hier_basis_func(x(i),y(i),xc,yc, &
                                (/grid%edge(grid%element(elem)%edge)%degree, &
                                  grid%element(elem)%degree /), "a", &
                                  basis,basisx,basisy,basisxx,basisyy)
      elseif (loc_which(2)==1 .or. loc_which(3)==1) then
         call p_hier_basis_func(x(i),y(i),xc,yc, &
                                (/grid%edge(grid%element(elem)%edge)%degree, &
                                  grid%element(elem)%degree /), "a", &
                                  basis,basisx,basisy)
      else
         call p_hier_basis_func(x(i),y(i),xc,yc, &
                                (/grid%edge(grid%element(elem)%edge)%degree, &
                                  grid%element(elem)%degree /), "a", &
                                  basis)
      endif

! compute the solution at this point

      loc_u  (i) = 0.0_my_real
      loc_ux (i) = 0.0_my_real
      loc_uy (i) = 0.0_my_real
      loc_uxx(i) = 0.0_my_real
      loc_uyy(i) = 0.0_my_real

      do j=1,VERTICES_PER_ELEMENT
         if (loc_which(1)==1) loc_u  (i) = loc_u  (i) + &
              basis  (j)*grid%vertex_solution(grid%element(elem)%vertex(j), &
              comp,eigen)
         if (loc_which(2)==1) loc_ux (i) = loc_ux (i) + &
              basisx (j)*grid%vertex_solution(grid%element(elem)%vertex(j), &
              comp,eigen)
         if (loc_which(3)==1) loc_uy (i) = loc_uy (i) + &
              basisy (j)*grid%vertex_solution(grid%element(elem)%vertex(j), &
              comp,eigen)
         if (loc_which(4)==1) loc_uxx(i) = loc_uxx(i) + &
              basisxx(j)*grid%vertex_solution(grid%element(elem)%vertex(j), &
              comp,eigen)
         if (loc_which(5)==1) loc_uyy(i) = loc_uyy(i) + &
              basisyy(j)*grid%vertex_solution(grid%element(elem)%vertex(j), &
              comp,eigen)
      end do
      isub = VERTICES_PER_ELEMENT
      do j=1,EDGES_PER_ELEMENT
         if (grid%edge(grid%element(elem)%edge(j))%degree < 2) cycle
         do k=1,grid%edge(grid%element(elem)%edge(j))%degree-1
            isub = isub + 1
            if (loc_which(1)==1) loc_u  (i) = loc_u  (i) + &
               basis  (isub)*grid%edge(grid%element(elem)%edge(j))%solution(k, &
               comp,eigen)
            if (loc_which(2)==1) loc_ux (i) = loc_ux (i) + &
               basisx (isub)*grid%edge(grid%element(elem)%edge(j))%solution(k, &
               comp,eigen)
            if (loc_which(3)==1) loc_uy (i) = loc_uy (i) + &
               basisy (isub)*grid%edge(grid%element(elem)%edge(j))%solution(k, &
               comp,eigen)
            if (loc_which(4)==1) loc_uxx(i) = loc_uxx(i) + &
               basisxx(isub)*grid%edge(grid%element(elem)%edge(j))%solution(k, &
               comp,eigen)
            if (loc_which(5)==1) loc_uyy(i) = loc_uyy(i) + &
               basisyy(isub)*grid%edge(grid%element(elem)%edge(j))%solution(k, &
               comp,eigen)
         end do
      end do
      if (grid%element(elem)%degree > 2) then
         do k=1,((grid%element(elem)%degree-2)*(grid%element(elem)%degree-1))/2
            isub = isub + 1
            if (loc_which(1)==1) loc_u  (i) = loc_u  (i) + &
               basis  (isub)*grid%element(elem)%solution(k,comp,eigen)
            if (loc_which(2)==1) loc_ux (i) = loc_ux (i) + &
               basisx (isub)*grid%element(elem)%solution(k,comp,eigen)
            if (loc_which(3)==1) loc_uy (i) = loc_uy (i) + &
               basisy (isub)*grid%element(elem)%solution(k,comp,eigen)
            if (loc_which(4)==1) loc_uxx(i) = loc_uxx(i) + &
               basisxx(isub)*grid%element(elem)%solution(k,comp,eigen)
            if (loc_which(5)==1) loc_uyy(i) = loc_uyy(i) + &
               basisyy(isub)*grid%element(elem)%solution(k,comp,eigen)
         end do
      endif

   else

      loc_u  (i) = 0.0_my_real
      loc_ux (i) = 0.0_my_real
      loc_uy (i) = 0.0_my_real
      loc_uxx(i) = 0.0_my_real
      loc_uyy(i) = 0.0_my_real

   endif
   deallocate(basis,basisx,basisy,basisxx,basisyy,stat=astat)
end do

! Copy the solution to the u's or send the solution to the invoker.
! Note that my_proc(invoker_procs) is not me, it is the invoker.

if (present(which)) then
   nr = size(x)*count(loc_which==1)
   allocate(send_real(nr),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in evaluate_soln_slave")
      stop
   endif
   j = 0
   if (loc_which(1) == 1) then
      send_real(j+1:j+size(x)) = loc_u
      j = j + size(x)
   endif
   if (loc_which(2) == 1) then
      send_real(j+1:j+size(x)) = loc_ux
      j = j + size(x)
   endif
   if (loc_which(3) == 1) then
      send_real(j+1:j+size(x)) = loc_uy
      j = j + size(x)
   endif
   if (loc_which(4) == 1) then
      send_real(j+1:j+size(x)) = loc_uxx
      j = j + size(x)
   endif
   if (loc_which(5) == 1) then
      send_real(j+1:j+size(x)) = loc_uyy
      j = j + size(x)
   endif
   call phaml_send(invoker_procs,my_proc(invoker_procs),have_it,size(x), &
                   send_real,nr,910)
   deallocate(send_real)
else
   if (present(u  )) u   = loc_u
   if (present(ux )) ux  = loc_ux
   if (present(uy )) uy  = loc_uy
   if (present(uxx)) uxx = loc_uxx
   if (present(uyy)) uyy = loc_uyy
endif

end subroutine evaluate_soln_slave

!          -----------------------
subroutine find_containing_element(x,y,have_it,elem,grid)
!          -----------------------

!----------------------------------------------------
! This routine looks for an element containing the point (x,y).  If it
! determines that (x,y) is not in an element owned by this processor, it
! quits looking and returns (have_it=0,elem=0).  Otherwise it returns
! have_it=1 and the element index in elem.  If the point falls on an
! element boundary, it is indeterminant as to which of the containing
! elements it returns.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x,y
integer, intent(out) :: have_it,elem
type(grid_type), intent(in) :: grid

!----------------------------------------------------
! Local variables:

integer :: e,root,i,neigh(EDGES_PER_ELEMENT)
real(my_real) :: bc(VERTICES_PER_ELEMENT)
!----------------------------------------------------
! Begin executable code

! find an element in the initial grid that contains (x,y)

! first look for it by moving from a triangle to a neighbor in the direction
! of the point, which is characterized by a negative barycentric coordinate

root = -10
e = grid%head_level_elem(1)
do
! if all barycentric coordinates are positive, it's in there
   bc = barycentric(x,y,e,grid,no_det=.true.)
   if (all(bc >= -roundoff_fudge*epsilon(0.0_my_real))) then
      root = e
      exit
   endif
   neigh = grid%initial_neighbor(:,e)
   do i=1,VERTICES_PER_ELEMENT
      if (bc(i) < 0) then
         e = neigh(i)
         if (e /= BOUNDARY) exit
      endif
   end do
   if (e == BOUNDARY) exit
end do

! if that failed, go through all the initial elements

if (root == -10) then
   e = grid%head_level_elem(1)
   do while (e /= END_OF_LIST)
      if (all(barycentric(x,y,e,grid,no_det=.true.) >= &
              -roundoff_fudge*epsilon(0.0_my_real))) then
         root = e
         exit
      endif
      e = grid%element(e)%next
   end do
endif

! go down the refinement tree to find the leaf element that contains (x,y)

if (root /= -10) then
   call find_containing_leaf(x,y,root,have_it,elem,grid)
else
   have_it = 0
   elem = 0
endif

end subroutine find_containing_element

!                    --------------------
recursive subroutine find_containing_leaf(x,y,root,have_it,elem,grid)
!                    --------------------

!----------------------------------------------------
! This routine recursively goes down the refinement tree, starting at root,
! to find a leaf element containing (x,y) as long as (x,y) may be in an
! element I own.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x,y
integer, intent(in) :: root
integer, intent(out) :: have_it,elem
type(grid_type), intent(in) :: grid

!----------------------------------------------------
! Local variables:

integer :: children(MAX_CHILD), i, allc(MAX_CHILD)
!----------------------------------------------------
! Begin executable code

! if root is a leaf, we've found it; but verify ownership

allc = ALL_CHILDREN
children = get_child_lid(grid%element(root)%gid,allc,grid%elem_hash)
if (children(1) == NO_CHILD) then
   if (grid%element(root)%iown) then
      have_it = 1
      elem = root
   else
      have_it = 0
      elem = 0
   endif
   return
endif

! otherwise, look at the barycentric coordinates of (x,y) in each child,
! until one is found where they are all positive

have_it = 0
elem = 0
do i=1,MAX_CHILD
   if (all(barycentric(x,y,children(i),grid,no_det=.true.) >= &
           -roundoff_fudge*epsilon(0.0_my_real))) then
      call find_containing_leaf(x,y,children(i),have_it,elem,grid)
      exit
   endif
end do

end subroutine find_containing_leaf

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

!          -------------------
subroutine evaluate_soln_local(grid,x,y,elem,comp,eigen,u,ux,uy,uxx,uyy,uxy)
!          -------------------

!----------------------------------------------------
! This routine evaluates the solution and/or derivatives of the solution
! on this processor.  The points (x,y) must all be in element elem.
! comp and eigen tell which components and eigenfunctions to return.
! The three dimensions of u correspond to comp, eigen and point.
! If ux or uy is present, then both must be present.  If uxx or uyy is
! present, then all 5 must be present.  This is only to reduce
! the number of forms of the call to basis_function, so this restriction can
! easily be removed.  Likewise, if uxy is present then u, ux and uy must be
! present and uxx and uyy either both present or not present.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
real(my_real), intent(in) :: x(:),y(:)
integer, intent(in) :: elem
integer, intent(in) :: comp(:),eigen(:)
real(my_real), optional, intent(out) :: u(:,:,:),ux(:,:,:),uy(:,:,:), &
                                        uxx(:,:,:),uyy(:,:,:),uxy(:,:,:)
! dimensions are (comp,eigen,point)

!----------------------------------------------------
! Local variables:

real(my_real) :: xc(3), yc(3)
real(my_real), allocatable :: basis(:,:),basisx(:,:),basisy(:,:),basisxx(:,:), &
                              basisyy(:,:),basisxy(:,:),solnvect(:,:)
integer :: i,j,k,l,p,nbasis,astat,isub
logical :: useblas

!----------------------------------------------------
! Begin executable code

! evaluate the basis functions at these points

nbasis = VERTICES_PER_ELEMENT
do j=1,EDGES_PER_ELEMENT
   if (grid%edge(grid%element(elem)%edge(j))%degree > 1) then
      nbasis = nbasis + grid%edge(grid%element(elem)%edge(j))%degree - 1
   endif
end do
if (grid%element(elem)%degree > 2) then
   nbasis = nbasis + &
         ((grid%element(elem)%degree-2)*(grid%element(elem)%degree-1))/2
endif
allocate(basis(nbasis,size(x)),basisx(nbasis,size(x)),basisy(nbasis,size(x)), &
         basisxx(nbasis,size(x)),basisyy(nbasis,size(x)), &
         basisxy(nbasis,size(x)),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in evaluate_soln_local")
   stop
endif
if (present(ux) .neqv. present(uy)) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("need both or neither of ux and uy in evaluate_soln_local")
   stop
endif
if (present(uxx) .or. present(uyy)) then
   if (.not. present(u) .or. .not. present(ux) .or. .not. present(uy) .or. &
       .not. present(uxx) .or. .not. present(uyy)) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("if uxx or uyy is present in evaluate_soln_local then all args except uxy must be present")
      stop
   endif
endif
if (present(uxy)) then
   if (.not. present(u) .or. .not. present(ux) .or. .not. present(uy) .or. &
       (present(uxx) .neqv. present(uyy))) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("evaluate_soln_local: uxy present requires u, ux and uy present and uxx and uyy both or neither present")
      stop
   endif
endif
xc = grid%vertex(grid%element(elem)%vertex)%coord%x
yc = grid%vertex(grid%element(elem)%vertex)%coord%y
if (present(u) .and. .not. present(ux)) then
   call p_hier_basis_func(x,y,xc,yc, &
                          (/grid%edge(grid%element(elem)%edge)%degree, &
                            grid%element(elem)%degree /),"a", &
                          basis)
elseif (.not. present(u) .and. present(ux)) then
   call p_hier_basis_func(x,y,xc,yc, &
                          (/grid%edge(grid%element(elem)%edge)%degree, &
                            grid%element(elem)%degree /),"a", &
                          basisx=basisx,basisy=basisy)
elseif (present(uxy) .and. present(uxx)) then
   call p_hier_basis_func(x,y,xc,yc, &
                          (/grid%edge(grid%element(elem)%edge)%degree, &
                            grid%element(elem)%degree /),"a", &
                          basis,basisx,basisy,basisxx,basisyy,basisxy)
elseif (present(uxy) .and. .not. present(uxx)) then
   call p_hier_basis_func(x,y,xc,yc, &
                          (/grid%edge(grid%element(elem)%edge)%degree, &
                            grid%element(elem)%degree /),"a", &
                          basis,basisx,basisy,basisxy=basisxy)
elseif (present(uxx)) then
   call p_hier_basis_func(x,y,xc,yc, &
                          (/grid%edge(grid%element(elem)%edge)%degree, &
                            grid%element(elem)%degree /),"a", &
                          basis,basisx,basisy,basisxx,basisyy)
else
   call p_hier_basis_func(x,y,xc,yc, &
                          (/grid%edge(grid%element(elem)%edge)%degree, &
                           grid%element(elem)%degree /),"a", &
                          basis,basisx,basisy)
endif

! compute the solution at this point

if (present(u)) u = 0.0_my_real
if (present(ux)) ux = 0.0_my_real
if (present(uy)) uy = 0.0_my_real
if (present(uxx)) uxx = 0.0_my_real
if (present(uyy)) uyy = 0.0_my_real
if (present(uxy)) uxy = 0.0_my_real

! special BLAS code for size(comp)==size(eigen)==1.  Copy solution values to
! an array and use GEMM.  TEMP Probably can do something similar with /=1.
! Also need the first 2 dimensions of the u's to be 1 since we're
! (illegally!) passing a rank 3 array to a rank 2 dummy argument.

useblas = size(comp)==1 .and. size(eigen)==1
if (present(u)) useblas = useblas .and. size(u,dim=1) == 1 .and. &
                          size(u,dim=2) == 1
if (present(ux)) useblas = useblas .and. size(ux,dim=1) == 1 .and. &
                          size(ux,dim=2) == 1
if (present(uy)) useblas = useblas .and. size(uy,dim=1) == 1 .and. &
                          size(uy,dim=2) == 1
if (present(uxx)) useblas = useblas .and. size(uxx,dim=1) == 1 .and. &
                          size(uxx,dim=2) == 1
if (present(uyy)) useblas = useblas .and. size(uyy,dim=1) == 1 .and. &
                          size(uyy,dim=2) == 1
if (present(uxy)) useblas = useblas .and. size(uxy,dim=1) == 1 .and. &
                          size(uxy,dim=2) == 1

if (useblas) then
   allocate(solnvect(1,nbasis),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in evaluate_soln_local")
      stop
   endif

   do j=1,VERTICES_PER_ELEMENT
      solnvect(1,j) = grid%vertex_solution(grid%element(elem)%vertex(j), &
                                           comp(1),eigen(1))
   end do
   isub = VERTICES_PER_ELEMENT
   do j=1,EDGES_PER_ELEMENT
      if (grid%edge(grid%element(elem)%edge(j))%degree < 2) cycle
      do k=1,grid%edge(grid%element(elem)%edge(j))%degree-1
         isub = isub + 1
         solnvect(1,isub) = grid%edge(grid%element(elem)%edge(j))%solution(k, &
                            comp(1),eigen(1))
      end do
   end do
   if (grid%element(elem)%degree > 2) then
      do k=1,((grid%element(elem)%degree-2)*(grid%element(elem)%degree-1))/2
         isub = isub + 1
         solnvect(1,isub) = grid%element(elem)%solution(k,comp(1),eigen(1))
      end do
   endif

   if (my_real == kind(1.0)) then
      if (present(u)) call sgemm("N","N",1,size(x),nbasis,1.0,solnvect,1, &
                                 basis,nbasis,1.0,u,1)
      if (present(ux)) call sgemm("N","N",1,size(x),nbasis,1.0,solnvect,1, &
                                  basisx,nbasis,1.0,ux,1)
      if (present(uy)) call sgemm("N","N",1,size(x),nbasis,1.0,solnvect,1, &
                                  basisy,nbasis,1.0,uy,1)
      if (present(uxx)) call sgemm("N","N",1,size(x),nbasis,1.0,solnvect,1, &
                                  basisxx,nbasis,1.0,uxx,1)
      if (present(uyy)) call sgemm("N","N",1,size(x),nbasis,1.0,solnvect,1, &
                                  basisyy,nbasis,1.0,uyy,1)
      if (present(uxy)) call sgemm("N","N",1,size(x),nbasis,1.0,solnvect,1, &
                                  basisxy,nbasis,1.0,uxy,1)
   elseif (my_real == kind(1.0d0)) then
      if (present(u)) call dgemm("N","N",1,size(x),nbasis,1.0d0,solnvect,1, &
                                 basis,nbasis,1.0d0,u,1)
      if (present(ux)) call dgemm("N","N",1,size(x),nbasis,1.0d0,solnvect,1, &
                                  basisx,nbasis,1.0d0,ux,1)
      if (present(uy)) call dgemm("N","N",1,size(x),nbasis,1.0d0,solnvect,1, &
                                  basisy,nbasis,1.0d0,uy,1)
      if (present(uxx)) call dgemm("N","N",1,size(x),nbasis,1.0d0,solnvect,1, &
                                  basisxx,nbasis,1.0d0,uxx,1)
      if (present(uyy)) call dgemm("N","N",1,size(x),nbasis,1.0d0,solnvect,1, &
                                  basisyy,nbasis,1.0d0,uyy,1)
      if (present(uxy)) call dgemm("N","N",1,size(x),nbasis,1.0d0,solnvect,1, &
                                  basisxy,nbasis,1.0d0,uxy,1)
   else
      ierr = PHAML_INTERNAL_ERROR
      call fatal("my_real is neither single nor double precision. Can't call GEMM")
      stop
   endif

   deallocate(solnvect,stat=astat)

else ! general case with size(comp)/=1 .or. size(eigen)/=1

   do j=1,VERTICES_PER_ELEMENT
      do p=1,size(x)
       do i=1,size(comp)
        do l=1,size(eigen)
         if (present(u)) u(i,l,p) = u(i,l,p) + &
            basis(j,p)*grid%vertex_solution(grid%element(elem)%vertex(j), &
                       comp(i),eigen(l))
         if (present(ux)) ux(i,l,p) = ux(i,l,p) + &
            basisx(j,p)*grid%vertex_solution(grid%element(elem)%vertex(j), &
                        comp(i),eigen(l))
         if (present(uy)) uy(i,l,p) = uy(i,l,p) + &
            basisy(j,p)*grid%vertex_solution(grid%element(elem)%vertex(j), &
                        comp(i),eigen(l))
         if (present(uxx)) uxx(i,l,p) = uxx(i,l,p) + &
            basisxx(j,p)*grid%vertex_solution(grid%element(elem)%vertex(j), &
                        comp(i),eigen(l))
         if (present(uyy)) uyy(i,l,p) = uyy(i,l,p) + &
            basisyy(j,p)*grid%vertex_solution(grid%element(elem)%vertex(j), &
                        comp(i),eigen(l))
         if (present(uxy)) uxy(i,l,p) = uxy(i,l,p) + &
            basisxy(j,p)*grid%vertex_solution(grid%element(elem)%vertex(j), &
                        comp(i),eigen(l))
        end do
       end do
      end do
   end do
   isub = VERTICES_PER_ELEMENT
   do j=1,EDGES_PER_ELEMENT
      if (grid%edge(grid%element(elem)%edge(j))%degree < 2) cycle
      do k=1,grid%edge(grid%element(elem)%edge(j))%degree-1
         isub = isub + 1
         do p=1,size(x)
          do i=1,size(comp)
           do l=1,size(eigen)
            if (present(u)) u(i,l,p) = u(i,l,p) + &
               basis(isub,p)*grid%edge(grid%element(elem)%edge(j))%solution(k, &
                             comp(i),eigen(l))
            if (present(ux)) ux(i,l,p) = ux(i,l,p) + &
               basisx(isub,p)*grid%edge(grid%element(elem)%edge(j))%solution(k,&
                              comp(i),eigen(l))
            if (present(uy)) uy(i,l,p) = uy(i,l,p) + &
               basisy(isub,p)*grid%edge(grid%element(elem)%edge(j))%solution(k,&
                              comp(i),eigen(l))
            if (present(uxx)) uxx(i,l,p) = uxx(i,l,p) + &
              basisxx(isub,p)*grid%edge(grid%element(elem)%edge(j))%solution(k,&
                              comp(i),eigen(l))
            if (present(uyy)) uyy(i,l,p) = uyy(i,l,p) + &
              basisyy(isub,p)*grid%edge(grid%element(elem)%edge(j))%solution(k,&
                              comp(i),eigen(l))
            if (present(uxy)) uxy(i,l,p) = uxy(i,l,p) + &
              basisxy(isub,p)*grid%edge(grid%element(elem)%edge(j))%solution(k,&
                              comp(i),eigen(l))
           end do
          end do
         end do
      end do
   end do
   if (grid%element(elem)%degree > 2) then
      do k=1,((grid%element(elem)%degree-2)*(grid%element(elem)%degree-1))/2
         isub = isub + 1
         do p=1,size(x)
          do i=1,size(comp)
           do l=1,size(eigen)
            if (present(u)) u(i,l,p) = u(i,l,p) + &
               basis(isub,p)*grid%element(elem)%solution(k,comp(i),eigen(l))
            if (present(ux)) ux(i,l,p) = ux(i,l,p) + &
               basisx(isub,p)*grid%element(elem)%solution(k,comp(i),eigen(l))
            if (present(uy)) uy(i,l,p) = uy(i,l,p) + &
               basisy(isub,p)*grid%element(elem)%solution(k,comp(i),eigen(l))
            if (present(uxx)) uxx(i,l,p) = uxx(i,l,p) + &
               basisxx(isub,p)*grid%element(elem)%solution(k,comp(i),eigen(l))
            if (present(uyy)) uyy(i,l,p) = uyy(i,l,p) + &
               basisyy(isub,p)*grid%element(elem)%solution(k,comp(i),eigen(l))
            if (present(uxy)) uxy(i,l,p) = uxy(i,l,p) + &
               basisxy(isub,p)*grid%element(elem)%solution(k,comp(i),eigen(l))
           end do
          end do
         end do
      end do
   endif

endif ! size(comp)==size(eigen)==1

deallocate(basis,basisx,basisy,basisxx,basisyy,basisxy,stat=astat)

end subroutine evaluate_soln_local

!          ----------------------
subroutine evaluate_oldsoln_local(x,y,u,ux,uy,uxx,uyy,comp,eigen)
!          ----------------------

!----------------------------------------------------
! This routine evaluates the old solution and/or derivatives of the old solution
! on this processor.
! If present, comp and/or eigen tell which component and eigenfunction to return
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x,y
real(my_real), optional, intent(out) :: u,ux,uy,uxx,uyy
integer, optional, intent(in) :: comp, eigen

!----------------------------------------------------
! Local variables:

integer :: j,k,nbasis,astat,isub,elem,loc_comp,loc_eigen
real(my_real), allocatable :: basis(:),basisx(:),basisy(:),basisxx(:),basisyy(:)
real(my_real) :: xc(3), yc(3), solnval
type(grid_type), pointer :: grid

!----------------------------------------------------
! Begin executable code

! make sure the old solution exists

grid => grid_for_old_soln
if (.not. grid%oldsoln_exists) then
   call warning("Old solution has not been saved, and hence cannot be evaluated", &
                "Returning 0.0 as old solution.")
   if (present(u)) u=0.0_my_real
   if (present(ux)) ux=0.0_my_real
   if (present(uy)) uy=0.0_my_real
   if (present(uxx)) uxx=0.0_my_real
   if (present(uyy)) uyy=0.0_my_real
   return
endif

! determine the component and eigenfunction

if (present(comp)) then
   if (comp <= 0) then
      call fatal("phaml_evaluate_oldsoln: comp must be a positive integer.  It is ", &
                  intlist=(/comp/))
      stop
   endif
   if (comp > grid%system_size) then
      call fatal("phaml_evaluate_oldsoln: comp must be no larger than system size.  They are ", &
                 intlist=(/comp,grid%system_size/))
      stop
   endif
   loc_comp = comp
else
   loc_comp = 1
endif
if (present(eigen)) then
   if (eigen <= 0) then
      call fatal("phaml_evaluate_oldsoln: eigen must be a positive integer.  It is ", &
                 intlist=(/eigen/))
      stop
   endif
   if (eigen > max(1,grid%num_eval)) then
      call fatal("phaml_evaluate_oldsoln: eigen must be no larger than the number of eigenvalues computed.  They are ", &
                 intlist=(/eigen,grid%num_eval/))
      stop
   endif
   loc_eigen = eigen
else
   loc_eigen = 1 
endif

! find the element containing (x,y)

call find_old_containing_element(x,y,elem,grid)
if (elem == 0) then
   if (present(u)) u=0.0_my_real
   if (present(ux)) ux=0.0_my_real
   if (present(uy)) uy=0.0_my_real
   if (present(uxx)) uxx=0.0_my_real
   if (present(uyy)) uyy=0.0_my_real
   return
endif

! evaluate the basis functions at the point

nbasis = VERTICES_PER_ELEMENT
do j=1,EDGES_PER_ELEMENT
   if (grid%edge(grid%element(elem)%edge(j))%degree > 1) then
      nbasis = nbasis + grid%edge(grid%element(elem)%edge(j))%degree - 1
   endif
end do
if (grid%element(elem)%degree > 2) then
   nbasis = nbasis + &
         ((grid%element(elem)%degree-2)*(grid%element(elem)%degree-1))/2
endif
allocate(basis(nbasis),basisx(nbasis),basisy(nbasis),basisxx(nbasis), &
         basisyy(nbasis),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in evaluate_oldsoln_local")
   stop
endif
xc = grid%vertex(grid%element(elem)%vertex)%coord%x
yc = grid%vertex(grid%element(elem)%vertex)%coord%y
if (present(uxx) .or. present(uyy)) then
   call p_hier_basis_func(x,y,xc,yc, &
                          (/grid%edge(grid%element(elem)%edge)%degree, &
                           grid%element(elem)%degree /),"a", &
                          basis,basisx,basisy,basisxx,basisyy)
elseif (present(ux) .or. present(uy)) then
   call p_hier_basis_func(x,y,xc,yc, &
                          (/grid%edge(grid%element(elem)%edge)%degree, &
                            grid%element(elem)%degree /),"a", &
                          basis,basisx,basisy)
else
   call p_hier_basis_func(x,y,xc,yc, &
                          (/grid%edge(grid%element(elem)%edge)%degree, &
                           grid%element(elem)%degree /),"a", &
                          basis)
endif

! compute the solution at this point

if (present(u)) u = 0.0_my_real
if (present(ux)) ux = 0.0_my_real
if (present(uy)) uy = 0.0_my_real
if (present(uxx)) uxx = 0.0_my_real
if (present(uyy)) uyy = 0.0_my_real

do j=1,VERTICES_PER_ELEMENT
   solnval = grid%vertex_oldsoln(grid%element(elem)%vertex(j),loc_comp,loc_eigen)
   if (present(u)) u = u + solnval*basis(j)
   if (present(ux)) ux = ux + solnval*basisx(j)
   if (present(uy)) uy = uy + solnval*basisy(j)
   if (present(uxx)) uxx = uxx + solnval*basisxx(j)
   if (present(uyy)) uyy = uyy + solnval*basisyy(j)
end do
isub = VERTICES_PER_ELEMENT
do j=1,EDGES_PER_ELEMENT
   if (grid%edge(grid%element(elem)%edge(j))%degree < 2) cycle
   do k=1,grid%edge(grid%element(elem)%edge(j))%degree-1
      isub = isub + 1
      if (.not. associated(grid%edge(grid%element(elem)%edge(j))%oldsoln)) cycle
      solnval = grid%edge(grid%element(elem)%edge(j))%oldsoln(k,loc_comp,loc_eigen)
      if (present(u)) u = u + solnval*basis(isub)
      if (present(ux)) ux = ux + solnval*basisx(isub)
      if (present(uy)) uy = uy + solnval*basisy(isub)
      if (present(uxx)) uxx = uxx + solnval*basisxx(isub)
      if (present(uyy)) uyy = uyy + solnval*basisyy(isub)
   end do
end do
if (grid%element(elem)%degree > 2) then
   do k=1,((grid%element(elem)%degree-2)*(grid%element(elem)%degree-1))/2
      isub = isub + 1
      if (.not. associated(grid%element(elem)%oldsoln)) cycle
      solnval = grid%element(elem)%oldsoln(k,loc_comp,loc_eigen)
      if (present(u)) u = u + solnval*basis(isub)
      if (present(ux)) ux = ux + solnval*basisx(isub)
      if (present(uy)) uy = uy + solnval*basisy(isub)
      if (present(uxx)) uxx = uxx + solnval*basisxx(isub)
      if (present(uyy)) uyy = uyy + solnval*basisyy(isub)
   end do
endif

deallocate(basis,basisx,basisy,basisxx,basisyy,stat=astat)

end subroutine evaluate_oldsoln_local

!          ---------------------------
subroutine find_old_containing_element(x,y,elem,grid)
!          ---------------------------

!----------------------------------------------------
! This routine looks for an "effective old leaf" element containing the point
! (x,y).  An element is an effective old leaf if oldsoln is allocated.
! If the point falls on an element boundary, it is indeterminant as to which
! of the containing elements it returns.
! If the point does not fall in any element, 0 is returned.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x,y
integer, intent(out) :: elem
type(grid_type), intent(in) :: grid

!----------------------------------------------------
! Local variables:

integer :: e,root,i,neigh(EDGES_PER_ELEMENT)
real(my_real) :: bc(VERTICES_PER_ELEMENT)
!----------------------------------------------------
! Begin executable code

! find an element in the initial grid that contains (x,y)

! first look for it by moving from a triangle to a neighbor in the direction
! of the point, which is characterized by a negative barycentric coordinate

root = -10
e = grid%head_level_elem(1)
do
! if all barycentric coordinates are positive, it's in there
   bc = barycentric(x,y,e,grid,no_det=.true.)
   if (all(bc >= -roundoff_fudge*epsilon(0.0_my_real))) then
      root = e
      exit
   endif
   neigh = grid%initial_neighbor(:,e)
   do i=1,VERTICES_PER_ELEMENT
      if (bc(i) < 0) then
         e = neigh(i)
         if (e /= BOUNDARY) exit
      endif
   end do
   if (e == BOUNDARY) exit
end do

! if that failed, go through all the initial elements

if (root == -10) then
   e = grid%head_level_elem(1)
   do while (e /= END_OF_LIST)
      if (all(barycentric(x,y,e,grid,no_det=.true.) >= &
              -roundoff_fudge*epsilon(0.0_my_real))) then
         root = e
         exit
      endif
      e = grid%element(e)%next
   end do
endif

! go down the refinement tree to find the effective leaf element that
! contains (x,y)

if (root /= -10) then
   call find_old_containing_leaf(x,y,root,elem,grid)
else
   elem = 0
endif

end subroutine find_old_containing_element

!                    ------------------------
recursive subroutine find_old_containing_leaf(x,y,root,elem,grid)
!                    ------------------------

!----------------------------------------------------
! This routine recursively goes down the refinement tree, starting at root,
! to find a element containing (x,y) that is marked as oldleaf
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x,y
integer, intent(in) :: root
integer, intent(out) :: elem
type(grid_type), intent(in) :: grid

!----------------------------------------------------
! Local variables:

integer :: children(MAX_CHILD), i, allc(MAX_CHILD)
!----------------------------------------------------
! Begin executable code

! if root is oldleaf, we've found it

if (grid%element(root)%oldleaf) then
   elem = root
   return
endif

! otherwise, check the children

allc = ALL_CHILDREN
children = get_child_lid(grid%element(root)%gid,allc,grid%elem_hash)
if (children(1) == NO_CHILD) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("evaluate_oldsoln could not find a containing element marked as oldleaf.")
   stop
endif

! look at the barycentric coordinates of (x,y) in each child,
! until one is found where they are all positive

elem = 0
do i=1,MAX_CHILD
   if (all(barycentric(x,y,children(i),grid,no_det=.true.) >= &
           -roundoff_fudge*epsilon(0.0_my_real))) then
      call find_old_containing_leaf(x,y,children(i),elem,grid)
      exit
   endif
end do

end subroutine find_old_containing_leaf

!          --------
subroutine copy_old(grid)
!          --------

!----------------------------------------------------
! This routine copies solution into oldsoln
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
!----------------------------------------------------
! Local variables:

integer :: lev, elem, astat, i, edge
!----------------------------------------------------
! Begin executable code

do lev=1,grid%nlev

! for each element

   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      if (grid%element(elem)%isleaf) then
         if (associated(grid%element(elem)%solution)) then

! if the element is a leaf and has a solution, then allocate oldsoln to the
! right size and copy solution to it

            if (associated(grid%element(elem)%oldsoln)) then
               if (size(grid%element(elem)%oldsoln,dim=1) /= &
                   size(grid%element(elem)%solution,dim=1) .or. &
                   size(grid%element(elem)%oldsoln,dim=2) /= &
                   size(grid%element(elem)%solution,dim=2) .or. &
                   size(grid%element(elem)%oldsoln,dim=3) /= &
                   size(grid%element(elem)%solution,dim=3)) then
                  deallocate(grid%element(elem)%oldsoln,stat=astat)
               endif
            endif
            if (.not. associated(grid%element(elem)%oldsoln)) then
               allocate(grid%element(elem)%oldsoln( &
                          size(grid%element(elem)%solution,dim=1), &
                          size(grid%element(elem)%solution,dim=2), &
                          size(grid%element(elem)%solution,dim=3)),stat=astat)
               if (astat /= 0) then
                  ierr = ALLOC_FAILED
                  call fatal("memory allocation failed in copy_old")
                  stop
               endif
            endif
            grid%element(elem)%oldsoln = grid%element(elem)%solution
            grid%element(elem)%oldleaf = .true.

         else

! if element is a leaf and does not have a solution, make sure it does not
! have oldsoln either

            if (associated(grid%element(elem)%oldsoln)) then
               deallocate(grid%element(elem)%oldsoln,stat=astat)
            endif
            grid%element(elem)%oldleaf = .true.

         endif
      else

! if element is not a leaf, make sure it does not have an oldsoln or say
! that it is an oldleaf

         if (associated(grid%element(elem)%oldsoln)) then
            deallocate(grid%element(elem)%oldsoln,stat=astat)
         endif
         grid%element(elem)%oldleaf = .false.

      endif

! for each edge

      do i=1,EDGES_PER_ELEMENT
         edge = grid%element(elem)%edge(i)

         if (grid%element(elem)%isleaf) then
            if (associated(grid%edge(edge)%solution)) then

! if the element is a leaf and the edge solution is allocated, make sure the
! edge oldsoln is the same size and copy solution to it

               if (associated(grid%edge(edge)%oldsoln)) then
                  if (size(grid%edge(edge)%oldsoln,dim=1) /= &
                      size(grid%edge(edge)%solution,dim=1) .or. &
                      size(grid%edge(edge)%oldsoln,dim=2) /= &
                      size(grid%edge(edge)%solution,dim=2) .or. &
                      size(grid%edge(edge)%oldsoln,dim=3) /= &
                      size(grid%edge(edge)%solution,dim=3)) then
                     deallocate(grid%edge(edge)%oldsoln,stat=astat)
                  endif
               endif
               if (.not. associated(grid%edge(edge)%oldsoln)) then
                  allocate(grid%edge(edge)%oldsoln( &
                             size(grid%edge(edge)%solution,dim=1), &
                             size(grid%edge(edge)%solution,dim=2), &
                             size(grid%edge(edge)%solution,dim=3)),stat=astat)
                  if (astat /= 0) then
                     ierr = ALLOC_FAILED
                     call fatal("memory allocation failed in copy_old")
                     stop
                  endif
               endif
               grid%edge(edge)%oldsoln = grid%edge(edge)%solution

            else

! if element is a leaf and edge does not have a solution, make sure it does not
! have oldsoln either

               if (associated(grid%edge(edge)%oldsoln)) then
                  deallocate(grid%edge(edge)%oldsoln,stat=astat)
               endif
            endif

! if element is not a leaf, don't change the edge because the neighbor might
! be a leaf and want the edge set, and it won't cause any damage to leave it

         endif
      end do ! next edge

      elem = grid%element(elem)%next
   end do ! next element

end do ! next level

! make sure the oldsoln is the same size as solution and copy

if (associated(grid%vertex_oldsoln)) then
   if (size(grid%vertex_oldsoln,dim=1) /= &
       size(grid%vertex_solution,dim=1) .or. &
       size(grid%vertex_oldsoln,dim=2) /= &
       size(grid%vertex_solution,dim=2) .or. &
       size(grid%vertex_oldsoln,dim=3) /= &
       size(grid%vertex_solution,dim=3)) then
      deallocate(grid%vertex_oldsoln,stat=astat)
   endif
endif
if (.not. associated(grid%vertex_oldsoln)) then
   allocate(grid%vertex_oldsoln(size(grid%vertex_solution,1), &
            size(grid%vertex_solution,2),size(grid%vertex_solution,3)), &
            stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in copy_old")
      stop
   endif
endif
grid%vertex_oldsoln = grid%vertex_solution

grid%oldsoln_exists = .true.

end subroutine copy_old

!          ---------------------
subroutine set_grid_for_old_soln(grid)
!          ---------------------

!----------------------------------------------------
! This routine sets grid_for_old_soln so that the grid is accessible even
! though it cannot be passed through the user routines pdecoef, etc.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), target :: grid

!----------------------------------------------------
! Begin executable code

grid_for_old_soln => grid

end subroutine set_grid_for_old_soln

end module evaluate
