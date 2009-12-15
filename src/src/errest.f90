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

module error_estimators

!----------------------------------------------------
! This module contains subroutines to compute error indicators and estimates.
! Many of these indicators are defined in
! Mark Ainsworth and J. Tinsley Oden, A Posteriori Error Estimation in
! Finite Element Analysis, John Wiley & Sons, New York, 2000.

! communication tags in this module are of the form 8xx
!----------------------------------------------------

!----------------------------------------------------
! Other modules used:

use global
use message_passing
use gridtype_mod
use hash_mod
use make_linsys
use evaluate
use basis_functions
use quadrature_rules
!----------------------------------------------------

implicit none
private
public all_error_indicators, error_indicator, error_estimate, set_edge_mass, &
       equilibrated_residual_ei_one

!----------------------------------------------------
! Non-module procedures used are:

interface

   subroutine pdecoefs(x,y,cxx,cxy,cyy,cx,cy,c,rs)
   use global
   real(my_real), intent(in) :: x,y
   real(my_real), intent(out) :: cxx(:,:),cxy(:,:),cyy(:,:),cx(:,:),cy(:,:), &
                                 c(:,:),rs(:)
   end subroutine pdecoefs

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

   function trues(x,y,comp,eigen)
   use global
   real (my_real), intent(in) :: x,y
   integer, intent(in) :: comp,eigen
   real (my_real) :: trues
   end function trues

   function truexs(x,y,comp,eigen)
   use global
   real (my_real), intent(in) :: x,y
   integer, intent(in) :: comp,eigen
   real (my_real) :: truexs
   end function truexs

   function trueys(x,y,comp,eigen)
   use global
   real (my_real), intent(in) :: x,y
   integer, intent(in) :: comp,eigen
   real (my_real) :: trueys
   end function trueys
end interface

!----------------------------------------------------
! The following variables are defined:

type(grid_type), pointer :: hold_grid
integer :: hold_elem

! for equilibrated residual patch type
integer, parameter :: INTERIOR_VERTEX     = 1, &
                      NEUMANN_NEUMANN     = 2, &
                      DIRICHLET_NEUMANN   = 3, &
                      DIRICHLET_DIRICHLET = 4
real(my_real), allocatable :: mu_K(:,:,:,:,:)
type edge_mass_type
   real(my_real), pointer :: M(:,:)
end type edge_mass_type
type(edge_mass_type), allocatable :: edge_mass(:)
!----------------------------------------------------

contains

!          ---------------
subroutine error_indicator(grid,elem,method,energy,work,energy_est,Linf,L2)
!          ---------------

!----------------------------------------------------
! This subroutine computes the error indicators and estimates for element elem.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: elem, method
real(my_real), intent(out), optional :: energy(:), work, energy_est(:), &
                                        Linf(:,:), L2(:,:)

!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

select case(method)

case (EXPLICIT_ERRIND)
   call explicit_ei(grid,elem,energy,work,energy_est,Linf,L2)

case (HIERARCHICAL_COEFFICIENT)
   call hier_coef_ei(grid,elem,energy,work,energy_est,Linf,L2)

case (TRUE_DIFF)
   call true_diff_ei(grid,elem,energy,work,energy_est,Linf,L2)

case (LOCAL_PROBLEM_H)
   call loc_prob_h_ei(grid,elem,method,energy,work,energy_est,Linf,L2)

case (LOCAL_PROBLEM_P)
   call loc_prob_p_ei(grid,elem,method,energy,work,energy_est,Linf,L2)

case (EQUILIBRATED_RESIDUAL)
   call equilibrated_residual_ei_one(grid,elem,energy,work,energy_est,Linf,L2)

case (INITIAL_CONDITION)
   call init_cond_ei(grid,elem,energy,work,energy_est,Linf,L2)

case (REFSOLN_ERREST)
   if (present(energy)) energy = 0.0_my_real
   if (present(work)) work = 1.0_my_real
   if (present(energy_est)) energy_est = 0.0_my_real
   if (present(Linf)) Linf = 0.0_my_real
   if (present(L2)) L2 = 0.0_my_real

case default
   call fatal("illegal value for error_indicator choice")
   stop

end select

end subroutine error_indicator

!          ------------
subroutine hier_coef_ei(grid,elem,energy,work,energy_est,Linf,L2)
!          ------------

!----------------------------------------------------
! This routine computes the hierarchical coefficient error indicator.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: elem
real(my_real), intent(out), optional :: energy(:), energy_est(:), Linf(:,:), &
                                        L2(:,:), work
!----------------------------------------------------
! Local variables:

integer :: i, d, n1, n2, n3, par, eigen, comp
real(my_real) :: error_indic(grid%system_size,max(1,grid%num_eval))
real(my_real) :: loc_energy, area, xvert(3), yvert(3)
!----------------------------------------------------
! Begin executable code

! for linear elements, determine the coefficient of the h-hierarchical basis
! that was added when this element was created

if (grid%element(elem)%degree == 1) then

   if (grid%element(elem)%level == 1) then ! cludge on level 1
      n1 = grid%element(elem)%vertex(1)
      n2 = grid%element(elem)%vertex(2)
   else
      par = hash_decode_key(grid%element(elem)%gid/2,grid%elem_hash)
      n1 = grid%element(par)%vertex(1)
      n2 = grid%element(par)%vertex(2)
   endif
   n3 = grid%element(elem)%vertex(3)
   error_indic = abs(grid%vertex_solution(n3,:,:) - &
          (grid%vertex_solution(n1,:,:) + grid%vertex_solution(n2,:,:))/2)

! for high order elements, use the coefficients of the highest order
! p-hierarchical bases

else

   error_indic = 0.0_my_real
   n1 = 0
   d = grid%element(elem)%degree
   if (d > 2) then
      do i=((d-2)*(d-3))/2+1,((d-1)*(d-2))/2
         error_indic = error_indic + grid%element(elem)%solution(i,:,:)**2
         n1 = n1 + 1
      end do
   endif

   do i=1,3
      d = grid%edge(grid%element(elem)%edge(i))%degree
      if (d >= 2) then
         error_indic = error_indic + grid%edge(grid%element(elem)%edge(i))%solution(d-1,:,:)**2
         n1 = n1 + 1
      endif
   end do

   error_indic = sqrt(error_indic)
   if (n1 > 0) error_indic = error_indic/n1

endif

do eigen=1,max(1,grid%num_eval)
   if (present(energy) .or. present(energy_est)) loc_energy = 0.0_my_real
   do comp=1,grid%system_size

! For energy norm, sum over the components of a system

      if (present(energy) .or. present(energy_est)) then
         loc_energy = loc_energy + error_indic(comp,eigen)**2
      endif

! L infinity norm is just the max of one value

      if (present(Linf)) Linf(comp,eigen) = error_indic(comp,eigen)

! L2 norm is the sqrt of the integral of the error squared over the triangle(s).
! Using a midpoint rule it is sqrt(area of triangles) * error(vert3)/3
! TEMP that's for linear elements.  Should do something different for
! high order elements.

      if (present(L2)) then
         xvert = grid%vertex(grid%element(elem)%vertex)%coord%x
         yvert = grid%vertex(grid%element(elem)%vertex)%coord%y
         area = abs(xvert(1)*(yvert(2)-yvert(3)) + &
                    xvert(2)*(yvert(3)-yvert(1)) + &
                    xvert(3)*(yvert(1)-yvert(2))) / 2
         L2(comp,eigen) = sqrt(area)*error_indic(comp,eigen)/3
      endif
   end do

   if (present(energy)) energy(eigen) = sqrt(loc_energy)
   if (present(energy_est)) energy_est(eigen) = sqrt(loc_energy)
end do

if (present(work)) then
   if (grid%element(elem)%mate == BOUNDARY) then
      work = ((grid%element(elem)%degree)*(grid%element(elem)%degree+1))/2
   else
      work = grid%element(elem)%degree**2
   endif
endif

end subroutine hier_coef_ei

!          ------------
subroutine true_diff_ei(grid,elem,energy,work,energy_est,Linf,L2)
!          ------------

!----------------------------------------------------
! This routine uses the difference between the computed solution and true
! solution as an error indicator.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: elem
real(my_real), intent(out), optional :: energy(:), energy_est(:), Linf(:,:), &
                                        L2(:,:), work
!----------------------------------------------------
! Local variables:

real (my_real), pointer :: qw(:),xq(:),yq(:)
real (my_real), allocatable :: err(:,:,:),errx(:,:,:),erry(:,:,:)
real(my_real) :: cxx(grid%system_size,grid%system_size), &
                 cxy(grid%system_size,grid%system_size), &
                 cyy(grid%system_size,grid%system_size), &
                 cx(grid%system_size,grid%system_size),  &
                 cy(grid%system_size,grid%system_size),  &
                 c(grid%system_size,grid%system_size),   &
                 rs(grid%system_size)
real(my_real) :: xc(3), yc(3), loc_energy(max(1,grid%num_eval))
integer :: nqp, jerr, ss, neigen, astat, comp, comp2, eigen, qp
!----------------------------------------------------
! Begin executable code

ss = grid%system_size
neigen = max(1,grid%num_eval)

! determine the quadrature rule

xc = grid%vertex(grid%element(elem)%vertex)%coord%x
yc = grid%vertex(grid%element(elem)%vertex)%coord%y
call quadrature_rule_tri(6,xc,yc,nqp,qw,xq,yq,jerr,stay_in=.true.)

! evaluate the solution at the quadrature points

allocate(err(ss,neigen,nqp),errx(ss,neigen,nqp),erry(ss,neigen,nqp),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in true_diff_ei")
   stop
endif
call evaluate_soln_local(grid,xq,yq,elem,(/(comp,comp=1,ss)/), &
                         (/(eigen,eigen=1,neigen)/),err,errx,erry)

! evaluate the error at the quadrature points

do eigen=1,neigen
   do comp=1,ss
      do qp=1,nqp
         err (comp,eigen,qp) = err (comp,eigen,qp) - &
                               trues (xq(qp),yq(qp),comp,eigen)
         errx(comp,eigen,qp) = errx(comp,eigen,qp) - &
                               truexs(xq(qp),yq(qp),comp,eigen)
         erry(comp,eigen,qp) = erry(comp,eigen,qp) - &
                               trueys(xq(qp),yq(qp),comp,eigen)
      end do
   end do
end do

! compute the energy norm of the error

if (present(energy) .or. present(energy_est)) then
   do qp=1,nqp
      call pdecoefs(xq(qp),yq(qp),cxx,cxy,cyy,cx,cy,c,rs)
      do eigen=1,neigen
         loc_energy(eigen) = 0.0_my_real
         do comp=1,ss
            do comp2=1,ss
               loc_energy(eigen) = loc_energy(eigen) + qw(qp)* &
                  (cxx(comp,comp2)*errx(comp,eigen,qp)*errx(comp2,eigen,qp) + &
                   cyy(comp,comp2)*erry(comp,eigen,qp)*erry(comp2,eigen,qp) + &
                   cxy(comp,comp2)*erry(comp,eigen,qp)*errx(comp2,eigen,qp) + &
                    cx(comp,comp2)*errx(comp,eigen,qp)* err(comp2,eigen,qp) + &
                    cy(comp,comp2)*erry(comp,eigen,qp)* err(comp2,eigen,qp) + &
                     c(comp,comp2)* err(comp,eigen,qp)* err(comp2,eigen,qp))
            end do
         end do
      end do
   end do
   if (present(energy)) energy = sqrt(abs(loc_energy))
   if (present(energy_est)) energy_est = sqrt(abs(loc_energy))
endif

! L-infinity norm

if (present(Linf)) then
   Linf = 0.0_my_real
   do eigen=1,neigen
      do comp=1,ss
         do qp=1,nqp
            Linf(comp,eigen) = max(Linf(comp,eigen),abs(err(comp,eigen,qp)))
         end do
      end do
   end do
endif

! L2 norm

if (present(L2)) then
   L2 = 0.0_my_real
   do eigen=1,neigen
      do comp=1,ss
         do qp=1,nqp
            L2(comp,eigen) = L2(comp,eigen) + qw(qp)*err(comp,eigen,qp)**2
         end do
      end do
   end do
   L2 = sqrt(L2)
endif

if (present(work)) then
   if (grid%element(elem)%mate == BOUNDARY) then
      work = ((grid%element(elem)%degree)*(grid%element(elem)%degree+1))/2
   else
      work = grid%element(elem)%degree**2
   endif
endif

deallocate(qw,xq,yq,err,errx,erry)

end subroutine true_diff_ei

!          -----------
subroutine explicit_ei(grid,elem,energy,work,energy_est,Linf,L2)
!          -----------

!----------------------------------------------------
! This routine computes the explicit error indicator given in chapter 2
! of Ainsworth and Oden.  The energy, L2 and Linf estimates are given on
! pages 22, 34 and 38 respectively.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: elem
real(my_real), intent(out), optional :: energy(:), work, energy_est(:), &
                                        Linf(:,:), L2(:,:)
!----------------------------------------------------
! Local variables:

real(my_real) :: xvert(VERTICES_PER_ELEMENT), yvert(VERTICES_PER_ELEMENT)
real(my_real), pointer :: qweights(:), xquad(:), yquad(:)
real(my_real), allocatable :: resid(:,:,:)
real(my_real) :: int_resid_L2sq_comps(grid%system_size,max(1,grid%num_eval)), &
                 int_resid_Linf_comps(grid%system_size,max(1,grid%num_eval)), &
                 bnd_resid_L2sq_comps(EDGES_PER_ELEMENT, &
                                      grid%system_size,max(1,grid%num_eval)), &
                 bnd_resid_Linf_comps(EDGES_PER_ELEMENT, &
                                      grid%system_size,max(1,grid%num_eval)), &
                 loc_energy(max(1,grid%num_eval))
real(my_real) :: elem_diam, edge_len(EDGES_PER_ELEMENT)
integer :: qorder, jerr, nqpoints, astat, i, e, edge, eigen, comp, soln, pow(4)
integer :: all_comp(grid%system_size), all_eigen(max(1,grid%num_eval))
!----------------------------------------------------
! Begin executable code

all_comp  = (/ (i,i=1,grid%system_size) /)
all_eigen = (/ (i,i=1,max(1,grid%num_eval)) /)

! use 4th order quadrature rules

qorder = 4

! copy the vertex coordinates to local arrays

xvert = grid%vertex(grid%element(elem)%vertex)%coord%x
yvert = grid%vertex(grid%element(elem)%vertex)%coord%y

! get a quadrature rule for the triangle

call quadrature_rule_tri(qorder,xvert,yvert,nqpoints,qweights,xquad,yquad,jerr, &
                         stay_in=.true.)

! compute the interior residual at the quadrature points

allocate(resid(grid%system_size,max(1,grid%num_eval),nqpoints),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in explicit_ei")
   stop
endif

call interior_residual(grid,xquad,yquad,elem,all_comp,all_eigen,resid)

! compute the L infinity and the square of the L2 norm of each comp,eigen
! interior residual

int_resid_L2sq_comps(:,:) = 0.0_my_real
int_resid_Linf_comps(:,:) = 0.0_my_real
do i=1,nqpoints
   int_resid_L2sq_comps = int_resid_L2sq_comps + qweights(i)*resid(:,:,i)**2
   int_resid_Linf_comps = max(int_resid_Linf_comps,abs(resid(:,:,i)))
end do

deallocate(qweights, xquad, yquad, resid, stat=astat)

! for each edge

do e=1,EDGES_PER_ELEMENT
   edge = grid%element(elem)%edge(e)

! copy the vertex coordinates and degrees to local arrays

   xvert(1:2) = grid%vertex(grid%edge(edge)%vertex)%coord%x
   yvert(1:2) = grid%vertex(grid%edge(edge)%vertex)%coord%y
   pow(e) = grid%edge(edge)%degree

! get a quadrature rule for the edge

   call quadrature_rule_line(qorder,xvert(1:2),yvert(1:2),nqpoints,qweights, &
                             xquad,yquad,jerr)

! compute the boundary residual at the quadrature points

   if (e==1) then
      allocate(resid(grid%system_size,max(1,grid%num_eval),nqpoints),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("memory allocation failed in explicit_ei")
         stop
      endif
   endif

   call boundary_residual(grid,xquad,yquad,elem,edge,all_comp,all_eigen,resid)

! compute the L infinity and the square of the L2 norm of each comp,eigen
! boundary residual

   bnd_resid_L2sq_comps(e,:,:) = 0.0_my_real
   bnd_resid_Linf_comps(e,:,:) = 0.0_my_real
   do i=1,nqpoints
      bnd_resid_L2sq_comps(e,:,:) = bnd_resid_L2sq_comps(e,:,:) + &
                                    qweights(i)*resid(:,:,i)**2
      bnd_resid_Linf_comps(e,:,:) = max(bnd_resid_Linf_comps(e,:,:), &
                                        abs(resid(:,:,i)))
   end do

   deallocate(qweights, xquad, yquad, stat=astat)

! compute the length of this edge

   edge_len(e) = sqrt((xvert(2)-xvert(1))**2 + (yvert(2)-yvert(1))**2)

end do

deallocate(resid,stat=astat)

! set the element diameter as the longest edge length

elem_diam = maxval(edge_len)

! compute error indicators and estimates

soln = 0
pow(4) = grid%element(elem)%degree
do eigen=1,max(1,grid%num_eval)
   if (present(energy) .or. present(energy_est)) then
      loc_energy(eigen) = elem_diam**2*sum(int_resid_L2sq_comps(:,eigen)) + &
         sum(elem_diam*sum(bnd_resid_L2sq_comps(:,:,eigen),dim=2))
   endif
   do comp=1,grid%system_size
      soln = soln + 1
! NOTE: this will not give exactly the same estimate as on page 38 of
! Ainsworth & Oden because the max of the individual parts should be taken
! outside this routine and then summed
      if (present(Linf)) then
         Linf(comp,eigen) = elem_diam**2*int_resid_Linf_comps(comp,eigen) + &
               maxval(edge_len*bnd_resid_Linf_comps(:,comp,eigen))
      endif
      if (present(L2)) then
         L2(comp,eigen) = elem_diam**4*int_resid_L2sq_comps(comp,eigen) + &
                 sum(elem_diam**3*bnd_resid_L2sq_comps(:,comp,eigen))
      endif
   end do
end do

if (present(energy)) energy = sqrt(loc_energy)
if (present(energy_est)) energy_est = sqrt(loc_energy)
if (present(L2)) L2 = sqrt(L2)

! the arbitrary constant in the error bound; these seem reasonable for
! many of the example and test programs

if (present(energy)) energy = energy/20
if (present(energy_est)) energy_est = energy_est/20
if (present(Linf)) Linf = Linf/100
if (present(L2)) L2 = L2/100

if (present(work)) then
   if (grid%element(elem)%mate == BOUNDARY) then
      work = ((grid%element(elem)%degree)*(grid%element(elem)%degree+1))/2
   else
      work = grid%element(elem)%degree**2
   endif
endif

end subroutine explicit_ei

!          -------------
subroutine loc_prob_h_ei(grid,elem,method,energy,work,energy_est,Linf,L2)
!          -------------

!----------------------------------------------------
! Set up and solve a local residual problem using h refinement of element elem
! and it's (possibly non-existent) mate, using only the red bases for the
! approximate solution, and define the h error indicator and estimates from
! the result.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in), target :: grid
integer, intent(in) :: elem, method
real(my_real), intent(out), optional :: energy(:),work,energy_est(:), &
                                        Linf(:,:),L2(:,:)
!----------------------------------------------------
! Local variables:

real(my_real) :: xvert(3), yvert(3), rmin, alpha, x, y, area, &
                 solut(grid%system_size,1,1), &
                 c(grid%system_size,grid%system_size), rs(grid%system_size)
real(my_real), allocatable :: elem_mat(:,:),elem_rs(:,:),full_mat(:,:), &
                              copy_mat(:,:),full_rs(:,:),temp(:,:)
integer :: mate, deg, degree(4), nred_elem, nred_full, elem_n, full_n, astat, &
           edge_type(3,grid%system_size), ss, i, i2, j, k, k1, k2, k3, k4, &
           info, bmark(3), itype(grid%system_size), nev
integer, allocatable :: renum(:), ipiv(:)
logical :: compatibly_divisible
!----------------------------------------------------
! Begin executable code

ss = grid%system_size
nev = max(1,grid%num_eval)
hold_grid => grid
hold_elem = elem

! lid of the mate, and flag to indicate if elem is compatibly divisible

if (grid%element(elem)%mate == BOUNDARY) then
   mate = BOUNDARY
   compatibly_divisible = .true.
else
   mate = hash_decode_key(grid%element(elem)%mate,grid%elem_hash)
   if (mate == HASH_NOT_FOUND) then
      mate = hash_decode_key(grid%element(elem)%mate/2,grid%elem_hash)
      compatibly_divisible = .false.
   else
      compatibly_divisible = .true.
   endif
endif

! degree of the elements

deg = grid%element(elem)%degree
if (mate /= BOUNDARY) deg = min(deg,grid%element(mate)%degree)
degree = (/ deg, deg, deg, deg /)

! number of red h-hierarchical bases over the refinement of triangles
! of this degree; both individual elements and all children

if (2*(deg/2) == deg) then
   nred_elem = (deg*(deg+2))/4
else
   nred_elem = ((deg+1)*(deg+1))/4
endif
if (mate == BOUNDARY) then
   nred_full = (deg*(deg+1))/2
else
   nred_full = deg**2
endif

elem_n = nred_elem*ss
full_n = nred_full*ss
allocate(elem_mat(elem_n,elem_n), elem_rs(elem_n,nev), full_mat(full_n,full_n),&
         full_rs(full_n,nev), copy_mat(full_n,full_n), renum(elem_n), &
         ipiv(full_n), temp(full_n,nev), stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in error estimates")
   stop
endif

full_mat = 0.0_my_real
full_rs = 0.0_my_real
elem_mat = 0.0_my_real
elem_rs = 0.0_my_real

! don't actually use rmin

rmin = 0.0_my_real

! for each of the children of elem and mate, compute the elemental matrix
! and right side, and assemble

! first child

xvert = (/ grid%vertex(grid%element(elem)%vertex(1))%coord%x, &
           grid%vertex(grid%element(elem)%vertex(3))%coord%x, &
          (grid%vertex(grid%element(elem)%vertex(1))%coord%x + &
           grid%vertex(grid%element(elem)%vertex(2))%coord%x)/2 /)
yvert = (/ grid%vertex(grid%element(elem)%vertex(1))%coord%y, &
           grid%vertex(grid%element(elem)%vertex(3))%coord%y, &
          (grid%vertex(grid%element(elem)%vertex(1))%coord%y + &
           grid%vertex(grid%element(elem)%vertex(2))%coord%y)/2 /)

edge_type(1,:) = INTERIOR
bmark(1) = 0
if (mate == BOUNDARY) then
   edge_type(2,:) = grid%edge_type(grid%element(elem)%edge(3),:)
   bmark(2) = grid%edge(grid%element(elem)%edge(3))%bmark
else
   edge_type(2,:) = INTERIOR
   bmark(2) = 0
endif
edge_type(3,:) = DIRICHLET
bmark(3) = huge(0)

if (nev > 1) then
   call elemental_matrix(grid, xvert, yvert, degree, edge_type, bmark, ss, &
                         0, .true., elem, rmin, "h", "r", elem_mat, &
                         elem_rs(:,1), loc_bconds_a=resid_dirich_bconds, &
                         extra_rhs=elem_rs(:,2:))
else
   call elemental_matrix(grid, xvert, yvert, degree, edge_type, bmark, ss, &
                         0, .true., elem, rmin, "h", "r", elem_mat, &
                         elem_rs(:,1), loc_bconds_a=resid_dirich_bconds)
endif

full_mat(1:elem_n,1:elem_n) = elem_mat
full_rs(1:elem_n,:) = elem_rs

! second child

xvert = (/ grid%vertex(grid%element(elem)%vertex(2))%coord%x, &
           grid%vertex(grid%element(elem)%vertex(3))%coord%x, &
          (grid%vertex(grid%element(elem)%vertex(1))%coord%x + &
           grid%vertex(grid%element(elem)%vertex(2))%coord%x)/2 /)
yvert = (/ grid%vertex(grid%element(elem)%vertex(2))%coord%y, &
           grid%vertex(grid%element(elem)%vertex(3))%coord%y, &
          (grid%vertex(grid%element(elem)%vertex(1))%coord%y + &
           grid%vertex(grid%element(elem)%vertex(2))%coord%y)/2 /)

edge_type(1,:) = INTERIOR
bmark(1) = 0
if (mate == BOUNDARY) then
   edge_type(2,:) = grid%edge_type(grid%element(elem)%edge(3),:)
   bmark(2) = grid%edge(grid%element(elem)%edge(3))%bmark
else
   edge_type(2,:) = INTERIOR
   bmark(2) = 0
endif
edge_type(3,:) = DIRICHLET
bmark(3) = huge(0)

if (nev > 1) then
   call elemental_matrix(grid, xvert, yvert, degree, edge_type, bmark, ss, &
                         0, .true., elem, rmin, "h", "r", elem_mat, &
                         elem_rs(:,1), loc_bconds_a=resid_dirich_bconds, &
                         extra_rhs=elem_rs(:,2:))
else
   call elemental_matrix(grid, xvert, yvert, degree, edge_type, bmark, ss, &
                         0, .true., elem, rmin, "h", "r", elem_mat, &
                         elem_rs(:,1), loc_bconds_a=resid_dirich_bconds)
endif

k1 = 1
k2 = elem_n/ss + 1
k3 = deg
do i=1,(deg+1)/2
   do j=1,deg-2*i+1
      renum(k1:k1+ss-1) = (/ (ss*(k2-1)+i2, i2=1,ss) /)
      k1 = k1+ss
      k2 = k2+1
   end do
   renum(k1:k1+ss-1) = (/ (ss*(k3-1)+i2, i2=1,ss) /)
   k1 = k1+ss
   k3 = k3+deg-2*i
end do

do i=1,elem_n
   do j=1,elem_n
      full_mat(renum(i),renum(j)) = full_mat(renum(i),renum(j)) + elem_mat(i,j)
   end do
   full_rs(renum(i),:) = full_rs(renum(i),:) + elem_rs(i,:)
end do

if (mate /= BOUNDARY) then

! third child

   if (compatibly_divisible) then

      xvert = (/ grid%vertex(grid%element(mate)%vertex(1))%coord%x, &
                 grid%vertex(grid%element(mate)%vertex(3))%coord%x, &
                (grid%vertex(grid%element(mate)%vertex(1))%coord%x + &
                 grid%vertex(grid%element(mate)%vertex(2))%coord%x)/2 /)
      yvert = (/ grid%vertex(grid%element(mate)%vertex(1))%coord%y, &
                 grid%vertex(grid%element(mate)%vertex(3))%coord%y, &
                (grid%vertex(grid%element(mate)%vertex(1))%coord%y + &
                 grid%vertex(grid%element(mate)%vertex(2))%coord%y)/2 /)

   else

      xvert = (/ grid%vertex(grid%element(elem)%vertex(1))%coord%x, &
                (grid%vertex(grid%element(mate)%vertex(1))%coord%x + &
                 grid%vertex(grid%element(mate)%vertex(2))%coord%x)/2, &
                (grid%vertex(grid%element(elem)%vertex(1))%coord%x + &
                 grid%vertex(grid%element(elem)%vertex(2))%coord%x)/2 /)
      yvert = (/ grid%vertex(grid%element(elem)%vertex(1))%coord%y, &
                (grid%vertex(grid%element(mate)%vertex(1))%coord%y + &
                 grid%vertex(grid%element(mate)%vertex(2))%coord%y)/2, &
                (grid%vertex(grid%element(elem)%vertex(1))%coord%y + &
                 grid%vertex(grid%element(elem)%vertex(2))%coord%y)/2 /)

   endif

   edge_type(1,:) = INTERIOR
   bmark(1) = 0
   edge_type(2,:) = INTERIOR
   bmark(2) = 0
   edge_type(3,:) = DIRICHLET
   bmark(3) = huge(0)

   if (nev > 1) then
      call elemental_matrix(grid, xvert, yvert, degree, edge_type, bmark, ss, &
                            0, .true., mate, rmin, "h", "r", elem_mat, &
                            elem_rs(:,1), loc_bconds_a=resid_dirich_bconds, &
                            extra_rhs=elem_rs(:,2:))
   else
      call elemental_matrix(grid, xvert, yvert, degree, edge_type, bmark, ss, &
                            0, .true., mate, rmin, "h", "r", elem_mat, &
                            elem_rs(:,1), loc_bconds_a=resid_dirich_bconds)
   endif

   k1 = 1
   k2 = (deg*(deg+1))/2 + 1
   k3 = 1
   do i=1,(deg+1)/2
      renum(k1:k1+ss-1) = (/ (ss*(k3-1)+i2, i2=1,ss) /)
      k1 = k1+ss
      k3 = k3+deg-2*i+2
      do j=1,deg-2*i+1
         renum(k1:k1+ss-1) = (/ (ss*(k2-1)+i2, i2=1,ss) /)
         k1 = k1+ss
         k2 = k2+1
      end do
   end do

   do i=1,elem_n
      do j=1,elem_n
        full_mat(renum(i),renum(j)) = full_mat(renum(i),renum(j))+elem_mat(i,j)
      end do
      full_rs(renum(i),:) = full_rs(renum(i),:) + elem_rs(i,:)
   end do

! fourth child

   if (compatibly_divisible) then

      xvert = (/ grid%vertex(grid%element(mate)%vertex(2))%coord%x, &
                 grid%vertex(grid%element(mate)%vertex(3))%coord%x, &
                (grid%vertex(grid%element(mate)%vertex(1))%coord%x + &
                 grid%vertex(grid%element(mate)%vertex(2))%coord%x)/2 /)
      yvert = (/ grid%vertex(grid%element(mate)%vertex(2))%coord%y, &
                 grid%vertex(grid%element(mate)%vertex(3))%coord%y, &
                (grid%vertex(grid%element(mate)%vertex(1))%coord%y + &
                 grid%vertex(grid%element(mate)%vertex(2))%coord%y)/2 /)

   else

      xvert = (/ grid%vertex(grid%element(elem)%vertex(2))%coord%x, &
                (grid%vertex(grid%element(mate)%vertex(1))%coord%x + &
                 grid%vertex(grid%element(mate)%vertex(2))%coord%x)/2, &
                (grid%vertex(grid%element(elem)%vertex(1))%coord%x + &
                 grid%vertex(grid%element(elem)%vertex(2))%coord%x)/2 /)
      yvert = (/ grid%vertex(grid%element(elem)%vertex(2))%coord%y, &
                (grid%vertex(grid%element(mate)%vertex(1))%coord%y + &
                 grid%vertex(grid%element(mate)%vertex(2))%coord%y)/2, &
                (grid%vertex(grid%element(elem)%vertex(1))%coord%y + &
                 grid%vertex(grid%element(elem)%vertex(2))%coord%y)/2 /)

   endif

   edge_type(1,:) = INTERIOR
   bmark(1) = 0
   edge_type(2,:) = INTERIOR
   bmark(2) = 0
   edge_type(3,:) = DIRICHLET
   bmark(3) = huge(0)

   if (nev > 1) then
      call elemental_matrix(grid, xvert, yvert, degree, edge_type, bmark, ss, &
                            0, .true., mate, rmin, "h", "r", elem_mat, &
                            elem_rs(:,1), loc_bconds_a=resid_dirich_bconds, &
                            extra_rhs=elem_rs(:,2:))
   else
      call elemental_matrix(grid, xvert, yvert, degree, edge_type, bmark, ss, &
                            0, .true., mate, rmin, "h", "r", elem_mat, &
                            elem_rs(:,1), loc_bconds_a=resid_dirich_bconds)
   endif

   k1 = 1
   if (2*(deg/2) == deg) then
      k2 = (deg*(3*deg+2))/4 + 1
   else
      k2 = ((deg+1)*(3*deg-1))/4 + 1
   endif
   k3 = elem_n/ss + 1
   k4 = (deg*(deg+1))/2 + deg-1
   do i=1,deg/2
      renum(k1:k1+ss-1) = (/ (ss*(k3-1)+i2, i2=1,ss) /)
      k1 = k1+ss
      k3 = k3+deg-2*i+1
      do j=1,deg-2*i
         renum(k1:k1+ss-1) = (/ (ss*(k2-1)+i2, i2=1,ss) /)
         k1 = k1+ss
         k2 = k2+1
      end do
      renum(k1:k1+ss-1) = (/ (ss*(k4-1)+i2, i2=1,ss) /)
      k1 = k1+ss
      k4 = k4+deg-2*i-1
   end do
   if (2*(deg/2) /= deg) then
      renum(k1:k1+ss-1) = (/ (elem_n-ss+i2, i2=1,ss) /)
   endif

   do i=1,elem_n
      do j=1,elem_n
        full_mat(renum(i),renum(j)) = full_mat(renum(i),renum(j))+elem_mat(i,j)
      end do
      full_rs(renum(i),:) = full_rs(renum(i),:) + elem_rs(i,:)
   end do

   copy_mat = full_mat

else ! mate is BOUNDARY

! replace rows corresponding to Dirichlet nodes with an identity row and
! the right side with the error

   copy_mat = full_mat

! points in first triangle

   k1 = 1
   do i=1,(deg+1)/2
      alpha = real(2*i-1,my_real)/real(2*deg,my_real)
      x = (1-alpha)*grid%vertex(grid%element(elem)%vertex(1))%coord%x + &
             alpha *grid%vertex(grid%element(elem)%vertex(2))%coord%x
      y = (1-alpha)*grid%vertex(grid%element(elem)%vertex(1))%coord%y + &
             alpha *grid%vertex(grid%element(elem)%vertex(2))%coord%y
      if (any(grid%edge_type(grid%element(elem)%edge(3),:) == DIRICHLET)) then
         call evaluate_soln_local(grid,(/x/),(/y/),elem,(/(j,j=1,ss)/),(/1/), &
                                  solut)
         call bconds(x,y,grid%edge(grid%element(elem)%edge(3))%bmark,itype,c,rs)
      endif
      do j=1,ss
         if (grid%edge_type(grid%element(elem)%edge(3),j) == DIRICHLET) then
            full_rs(k1+j-1,:) = rs(j) - solut(j,1,1)
            full_mat(k1+j-1,:) = 0.0_my_real
            full_mat(k1+j-1,k1+j-1) = 1.0_my_real
         endif
      end do
      k1 = k1 + ss*(deg-2*(i-1))
   end do

! points in second triangle

   k1 = elem_n + 1
   do i=1,deg/2
      alpha = real(2*i-1,my_real)/real(2*deg,my_real)
      x = (1-alpha)*grid%vertex(grid%element(elem)%vertex(2))%coord%x + &
             alpha *grid%vertex(grid%element(elem)%vertex(1))%coord%x
      y = (1-alpha)*grid%vertex(grid%element(elem)%vertex(2))%coord%y + &
             alpha *grid%vertex(grid%element(elem)%vertex(1))%coord%y
      if (any(grid%edge_type(grid%element(elem)%edge(3),:) == DIRICHLET)) then
         call evaluate_soln_local(grid,(/x/),(/y/),elem,(/(j,j=1,ss)/),(/1/), &
                                  solut)
         call bconds(x,y,grid%edge(grid%element(elem)%edge(3))%bmark,itype,c,rs)
      endif
      do j=1,ss
         if (grid%edge_type(grid%element(elem)%edge(3),j) == DIRICHLET) then
            full_rs(k1+j-1,:) = rs(j) - solut(j,1,1)
            full_mat(k1+j-1,:) = 0.0_my_real
            full_mat(k1+j-1,k1+j-1) = 1.0_my_real
         endif
      end do
      k1 = k1 + ss*(deg+1-2*i)
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
   call fatal("LAPACK returned error during h-hierarchical error estimate computation",intlist=(/info/))
   stop
endif

! set error indicators and estimates

!                                 T
! energy error indicator is sqrt(e Ae)

if (my_real == kind(1.0)) then
   call sgemm("N","N",full_n,nev,full_n,1.0_my_real,copy_mat, &
             full_n,full_rs,full_n,0.0_my_real,temp,full_n)
elseif (my_real == kind(1.0d0)) then
   call dgemm("N","N",full_n,nev,full_n,1.0_my_real,copy_mat, &
             full_n,full_rs,full_n,0.0_my_real,temp,full_n)
else
   ierr = PHAML_INTERNAL_ERROR
   call fatal("my_real is neither single nor double precision. Can't call GEMM")
   stop
endif

do k=1,nev

   if (present(energy)) &
      energy(k) = sqrt(abs(dot_product(temp(:,k),full_rs(:,k))))
   if (present(energy_est)) &
      energy_est(k) = sqrt(abs(dot_product(temp(:,k),full_rs(:,k))))

end do

! use the maximum for L infinity norm

if (present(Linf)) then
   do k=1,nev
      do j=1,ss
         Linf(j,k) = maxval(abs(full_rs(j:full_n:ss,k)))
      end do
   end do
endif

! for L2 estimate, use the maximum as function value in L2 integral

if (present(L2)) then
   if (grid%element(elem)%iown) then
      xvert = grid%vertex(grid%element(elem)%vertex)%coord%x
      yvert = grid%vertex(grid%element(elem)%vertex)%coord%y
      area = abs(xvert(1)*(yvert(2)-yvert(3)) + &
                 xvert(2)*(yvert(3)-yvert(1)) + &
                 xvert(3)*(yvert(1)-yvert(2))) / 2
   else
      area = 0.0_my_real
   endif
   do k=1,nev
      do j=1,ss
         L2(j,k) = sqrt(area)*maxval(abs(full_rs(j:full_n:ss,k)))
      end do
   end do
endif

! clean up memory

deallocate(elem_mat,elem_rs,full_mat,full_rs,copy_mat,renum,ipiv,temp,stat=astat)

if (present(work)) then
   if (grid%element(elem)%mate == BOUNDARY) then
      work = ((grid%element(elem)%degree)*(grid%element(elem)%degree+1))/2
   else
      work = grid%element(elem)%degree**2
   endif
endif

end subroutine loc_prob_h_ei

!          -------------
subroutine loc_prob_p_ei(grid,elem,method,energy,work,energy_est,Linf,L2)
!          -------------

!----------------------------------------------------
! Set up and solve a local residual problem using p refinement of element elem,
! using only the red bases for the approximate solution, and define the p
! error indicator and estimates from the result.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in), target :: grid
integer, intent(in) :: elem, method
real(my_real), intent(out), optional :: energy(:),work,energy_est(:), &
                                        Linf(:,:),L2(:,:)
!----------------------------------------------------
! Local variables:

integer :: ss, i, j, k, degree(1+EDGES_PER_ELEMENT), bmark(EDGES_PER_ELEMENT), &
           edge_type(EDGES_PER_ELEMENT,grid%system_size), nred, n, astat, &
           nev, info
real(my_real) :: xvert(VERTICES_PER_ELEMENT), yvert(VERTICES_PER_ELEMENT), &
                 rmin, area
integer, allocatable :: ipiv(:)
real(my_real), allocatable :: mat(:,:), copy_mat(:,:), rs(:,:), temp(:), &
                              temp2(:,:)
!----------------------------------------------------
! Begin executable code

ss = grid%system_size
nev = max(1,grid%num_eval)

hold_grid => grid
hold_elem = elem

! don't actually use rmin

rmin = 0.0_my_real

! the matrix is the elemental matrix for element elem, with modified natural
! boundary conditions and the degree increased by 1

xvert = grid%vertex(grid%element(elem)%vertex)%coord%x
yvert = grid%vertex(grid%element(elem)%vertex)%coord%y

degree(1+EDGES_PER_ELEMENT) = grid%element(elem)%degree + 1
nred = degree(1+EDGES_PER_ELEMENT) - 2
do i=1,EDGES_PER_ELEMENT
   if (grid%edge(grid%element(elem)%edge(i))%degree >= degree(1+EDGES_PER_ELEMENT)) then
      degree(i) = 1 ! forces no red bases on this side
   else
      degree(i) = degree(1+EDGES_PER_ELEMENT)
      nred = nred + 1
   endif
   edge_type(i,:) = grid%edge_type(grid%element(elem)%edge(i),:)
   bmark(i) = grid%edge(grid%element(elem)%edge(i))%bmark
end do
where (edge_type(:,1) == INTERIOR) bmark = huge(0)
where (edge_type == INTERIOR) edge_type = NATURAL

! special case -- if nred is 0 then the element is degree 1 and all edges
! are greater than 1.  In this case use degree 2 for the edges so that an
! error estimate can be computed

if (nred == 0) then
   degree = 2
   nred = 3
endif

! another special case -- if the element is degree 1 and the only edges that
! are not greater than 1 are on the Dirichlet boundary, use degree 2

if (grid%element(elem)%degree == 1 .and. &
    (edge_type(1,1) == DIRICHLET .or. degree(1) == 1) .and. &
    (edge_type(2,1) == DIRICHLET .or. degree(2) == 1) .and. &
    (edge_type(3,1) == DIRICHLET .or. degree(3) == 1)) then
   degree = 2
   nred = 3
endif

n = nred*ss
allocate(mat(n,n),copy_mat(n,n),rs(n,nev),ipiv(n),temp(n),temp2(n,nev), &
         stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in error estimates")
   stop 
endif

if (nev > 1) then
   call elemental_matrix(grid, xvert, yvert, degree, edge_type, bmark, ss, &
                         0, .true., elem, rmin, "p", "r", mat, &
                         rs(:,1), loc_bconds_a=resid_natural_bconds, &
                         extra_rhs=rs(:,2:))
else
   call elemental_matrix(grid, xvert, yvert, degree, edge_type, bmark, ss, &
                         0, .true., elem, rmin, "p", "r", mat, &
                         rs(:,1), loc_bconds_a=resid_natural_bconds)
endif

! adjust for Dirichlet boundary conditions on the real boundary

k = 1
do i=1,EDGES_PER_ELEMENT
   if (degree(i) >=2) then
      do j=1,ss
         if (edge_type(i,j) == DIRICHLET) then
            mat(k,:) = 0.0_my_real
            mat(k,k) = 1.0_my_real
            rs(k,:) = 0.0_my_real
         endif
         k = k+1
      end do
   endif
end do

! solve the system

copy_mat = mat
if (my_real == kind(0.0)) then
   call sgesv(n,nev,copy_mat,n,ipiv,rs,n,info)
elseif (my_real == kind(0.0d0)) then
   call dgesv(n,nev,copy_mat,n,ipiv,rs,n,info)
else
   ierr = USER_INPUT_ERROR
   call fatal("kind of real must be single or double to use LAPACK for error estimate")
   stop
endif

if (info /= 0) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("LAPACK returned error during p-hierarchical error estimate computation",intlist=(/info/))
   stop
endif

! set error indicators and estimates

!                                 T
! energy error indicator is sqrt(e Ae)

if (my_real == kind(1.0)) then
   call sgemm("N","N",n,nev,n,1.0_my_real,mat,n,rs,n,0.0_my_real,temp2,n)
elseif (my_real == kind(1.0d0)) then
   call dgemm("N","N",n,nev,n,1.0_my_real,mat,n,rs,n,0.0_my_real,temp2,n)
else
   ierr = PHAML_INTERNAL_ERROR
   call fatal("my_real is neither single nor double precision. Can't call GEMM")
   stop
endif

do k=1,nev

   if (present(energy)) &
      energy(k) = sqrt(abs(dot_product(temp2(:,k),rs(:,k))))
   if (present(energy_est)) &
      energy_est(k) = sqrt(abs(dot_product(temp2(:,k),rs(:,k))))

end do

! use midpoint value for L infinity and L2 norms

if (present(Linf) .or. present(L2)) then
   call p_hier_basis_func((xvert(1)+xvert(2)+xvert(3))/3, &
                          (yvert(1)+yvert(2)+yvert(3))/3, &
                          xvert,yvert,degree,"r",temp)
endif

if (present(Linf)) then
   do k=1,nev
      do j=1,ss
         Linf(j,k) = abs(sum(rs(j:n:ss,k)*temp(1:n/ss)))
      end do
   end do
endif

if (present(L2)) then
   area = abs(xvert(1)*(yvert(2)-yvert(3)) + &
              xvert(2)*(yvert(3)-yvert(1)) + &
              xvert(3)*(yvert(1)-yvert(2))) / 2
   do k=1,nev
      do j=1,ss
         L2(j,k) = sqrt(area)*abs(sum(rs(j:n:ss,k)*temp(1:n/ss)))
      end do
   end do
endif

! clean up memory

deallocate(mat,copy_mat,rs,ipiv,temp,temp2,stat=astat)

if (present(work)) then
   if (grid%element(elem)%mate == BOUNDARY) then
      work = ((grid%element(elem)%degree)*(grid%element(elem)%degree+1))/2
   else
      work = grid%element(elem)%degree**2
   endif
endif

end subroutine loc_prob_p_ei

!          -------------------
subroutine resid_dirich_bconds(x,y,bmark,itype,c,rs,extra_rs)
!          -------------------

!----------------------------------------------------
! This routine returns boundary conditions for the residual local Dirichlet
! problem.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x(:),y(:)
integer, intent(in) :: bmark 
integer, intent(out) :: itype(:)
real(my_real), intent(out) :: c(:,:,:),rs(:,:)
real(my_real), intent(out), optional :: extra_rs(:,:,:)
!----------------------------------------------------
! Local variables:

type(grid_type), pointer :: grid
integer :: elem, ss, edge, i, j, k
real(my_real) ::  u(hold_grid%system_size,max(1,hold_grid%num_eval),size(x)), &
                 ux(hold_grid%system_size,max(1,hold_grid%num_eval),size(x)), &
                 uy(hold_grid%system_size,max(1,hold_grid%num_eval),size(x)), normal(2), &
                 cxx(hold_grid%system_size,hold_grid%system_size,size(x)), &
                 cxy(hold_grid%system_size,hold_grid%system_size,size(x)), &
                 cyy(hold_grid%system_size,hold_grid%system_size,size(x)), &
                 cx (hold_grid%system_size,hold_grid%system_size,size(x)), &
                 cy (hold_grid%system_size,hold_grid%system_size,size(x)), &
                 cc (hold_grid%system_size,hold_grid%system_size,size(x)), &
                 crs(hold_grid%system_size,size(x))
!----------------------------------------------------
! Begin executable code

! bmark = huge(0) indicates the special Dirichlet boundary of the local problem

if (bmark == huge(0)) then

   itype = DIRICHLET
   c = 0.0_my_real
   rs = 0.0_my_real
   if (present(extra_rs)) extra_rs = 0.0_my_real

else

! otherwise evaluate the given boundary conditions and change Neuman and
! mixed right sides to the residual.  Dirichlet conditions don't matter because
! they get reset in residual_h.

   grid => hold_grid
   elem = hold_elem

   ss = grid%system_size

   do i=1,size(x)
      call bconds(x(i),y(i),bmark,itype,c(:,:,i),rs(:,i))
   end do
   if (present(extra_rs)) then
      do i=1,size(extra_rs,dim=2)
         extra_rs(:,i,:) = rs
      end do
   endif

   if (any(itype==NATURAL) .or. any(itype==MIXED)) then
      do i=1,size(x)
         call pdecoefs(x(i),y(i),cxx(:,:,i),cxy(:,:,i),cyy(:,:,i),cx(:,:,i), &
                       cy(:,:,i),cc(:,:,i),crs(:,i))
      end do
      call evaluate_soln_local(grid,x,y,elem,(/(i,i=1,ss)/),(/(i,i=1,max(1,grid%num_eval))/), &
                               u=u,ux=ux,uy=uy)
      call old_compute_normal(grid,x(1),y(1),elem,edge,normal)
      do i=1,ss
         if (itype(i)==NATURAL) then
            rs(i,:) = rs(i,:) - (normal(1)*(ux(i,1,:)*cxx(i,i,:)+uy(i,1,:)*cxy(i,i,:)) &
                               + normal(2)*uy(i,1,:)*cyy(i,i,:))
            if (present(extra_rs)) then
               do k=1,size(extra_rs,dim=2)
                  extra_rs(i,k,:) = extra_rs(i,k,:) - &
                     (normal(1)*(ux(i,k+1,:)*cxx(i,i,:) + uy(i,k+1,:)*cxy(i,i,:)) &
                    + normal(2)*uy(i,k+1,:)*cyy(i,i,:))
               end do
            endif
         elseif (itype(i)==MIXED) then
            rs(i,:) = rs(i,:) - (normal(1)*(ux(i,1,:)*cxx(i,i,:)+uy(i,1,:)*cxy(i,i,:)) &
                               + normal(2)*uy(i,1,:)*cyy(i,i,:))
            do j=1,ss
               rs(i,:) = rs(i,:) - cc(i,j,:)*u(j,1,:)
            end do
            if (present(extra_rs)) then
               do k=1,size(extra_rs,dim=2)
                  extra_rs(i,k,:) = extra_rs(i,k,:) - &
                     (normal(1)*(ux(i,k+1,:)*cxx(i,i,:)+uy(i,k+1,:)*cxy(i,i,:)) &
                    + normal(2)*uy(i,k+1,:)*cyy(i,i,:))
                  do j=1,ss
                     extra_rs(i,k,:) = extra_rs(i,k,:) - cc(i,j,:)*u(j,k+1,:)
                  end do
               end do
            endif
         endif
      end do
   endif

endif

end subroutine resid_dirich_bconds

!          --------------------
subroutine resid_natural_bconds(x,y,bmark,itype,c,rs,extra_rs)
!          --------------------

!----------------------------------------------------
! This routine returns boundary conditions for the residual local Neuman
! problem.  All points (x,y) must be on the same edge.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x(:),y(:)
integer, intent(in) :: bmark 
integer, intent(out) :: itype(:)
real(my_real), intent(out) :: c(:,:,:),rs(:,:)
real(my_real), intent(out), optional :: extra_rs(:,:,:)
!----------------------------------------------------
! Local variables:

type(grid_type), pointer :: grid
integer :: elem, ss, i, j, k, edge, neigh, neighbors(3)
real(my_real) :: u(hold_grid%system_size,max(1,hold_grid%num_eval),size(x)), &
                 ux(hold_grid%system_size,max(1,hold_grid%num_eval),size(x)), &
                 uy(hold_grid%system_size,max(1,hold_grid%num_eval),size(x)), normal(2), &
                 cxx(hold_grid%system_size,hold_grid%system_size,size(x)), &
                 cxy(hold_grid%system_size,hold_grid%system_size,size(x)), &
                 cyy(hold_grid%system_size,hold_grid%system_size,size(x)), &
                 cx (hold_grid%system_size,hold_grid%system_size,size(x)), &
                 cy (hold_grid%system_size,hold_grid%system_size,size(x))
real(my_real) :: xx(size(x)),yy(size(y)) ! TEMP081103 battery cludge
!----------------------------------------------------
! Begin executable code

grid => hold_grid
elem = hold_elem

ss = grid%system_size

! determine which edge contains (x,y) and compute the outward unit normal

call old_compute_normal(grid,x(1),y(1),elem,edge,normal)

! bmark = huge(0) indicates the special Neuman boundary of the local problem

if (bmark == huge(0)) then

   itype = NATURAL
   c = 0.0_my_real
   rs = 0.0_my_real
   if (present(extra_rs)) extra_rs = 0.0_my_real

! TEMP081103 battery cludge; move the points slightly inside the triangle
!            in case this is an edge with discontinuous coefficients

   if (battery_cludge) then
      xx = x - normal(1)*1.0e-13_my_real
      yy = y - normal(2)*1.0e-13_my_real
   else
      xx = x
      yy = y
   endif

! evaluate the coefficients of the pde to use in the natural b.c.

   do i=1,size(x)
      call pdecoefs(xx(i),yy(i),cxx(:,:,i),cxy(:,:,i),cyy(:,:,i),cx(:,:,i), &
                    cy(:,:,i),c(:,:,i),rs(:,i))
   end do

! find the neighbor that shares the (x,y) edge; it has the same index as edge

   neighbors = get_neighbors(elem,grid)
   neigh = neighbors(edge)

! compute the derivatives of the solution in this element

   call evaluate_soln_local(grid,x,y,elem,(/(i,i=1,ss)/),(/(i,i=1,max(1,grid%num_eval))/), &
                            ux=ux,uy=uy)

! add normal derivative of solution to right side

   do i=1,hold_grid%system_size
      rs(i,:) = normal(1)*(ux(i,1,:)*cxx(i,i,:)+uy(i,1,:)*cxy(i,i,:)) &
              + normal(2)*uy(i,1,:)*cyy(i,i,:)
      if (present(extra_rs)) then
         do k=1,size(extra_rs,dim=2)
            extra_rs(i,k,:) = normal(1)*(ux(i,k+1,:)*cxx(i,i,:)+uy(i,k+1,:)*cxy(i,i,:)) &
                            + normal(2)*uy(i,k+1,:)*cyy(i,i,:)
         end do
      endif
   end do

! TEMP081103 battery cludge; move the points slightly outside the triangle
!            and reevaluate the coefficients, in case this is an edge with
!            discontinuous coefficients

   if (battery_cludge) then
      xx = x + normal(1)*1.0e-13_my_real
      yy = y + normal(2)*1.0e-13_my_real
      do i=1,size(x)
         call pdecoefs(xx(i),yy(i),cxx(:,:,i),cxy(:,:,i),cyy(:,:,i),cx(:,:,i), &
                       cy(:,:,i),c(:,:,i),rs(:,i))
      end do
   endif

! compute the derivatives of the solution in the neighbor

   call evaluate_soln_local(grid,x,y,neigh,(/(i,i=1,ss)/),(/(i,i=1,max(1,grid%num_eval))/), &
                            ux=ux,uy=uy)

! set rs to be -1/2 jump
! TEMP I don't know why this should be -1/2 the jump, but adding the extra rs
!      is the only thing that works

   do i=1,hold_grid%system_size
      rs(i,:) = rs(i,:) - &
                (normal(1)*(ux(i,1,:)*cxx(i,i,:)+uy(i,1,:)*cxy(i,i,:)) &
               + normal(2)*uy(i,1,:)*cyy(i,i,:) - rs(i,:))/2
      if (present(extra_rs)) then
         do k=1,size(extra_rs,dim=2)
            extra_rs(i,k,:) = extra_rs(i,k,:) &
                            - (normal(1)*(ux(i,k+1,:)*cxx(i,i,:) + &
                                          uy(i,k+1,:)*cxy(i,i,:)) + &
                               normal(2)*uy(i,k+1,:)*cyy(i,i,:) - extra_rs(i,k,:))/2
         end do
      endif
   end do

else

! If this is not a special Neuman boundary,  evaluate the given boundary
! conditions and subtract the normal derivative of the solution on natural
! boundaries to get the residual

   do i=1,size(x)
      call bconds(x(i),y(i),bmark,itype,c(:,:,i),rs(:,i))
   end do
   if (present(extra_rs)) then
      do i=1,size(extra_rs,dim=2)
         extra_rs(:,i,:) = rs
      end do
   endif

   if (any(itype==MIXED)) then
      call evaluate_soln_local(grid,x,y,elem,(/(j,j=1,ss)/),(/(i,i=1,max(1,grid%num_eval))/),u)
   endif

   do i=1,ss
      if (itype(i) == NATURAL .or. itype(i) == MIXED) then
         call evaluate_soln_local(grid,x,y,elem,(/i/),(/(i,i=1,max(1,grid%num_eval))/), &
                                  ux=ux,uy=uy)
         rs(i,:) = rs(i,:) - normal(1)*(ux(1,1,:)*cxx(i,i,:)+uy(1,1,:)*cxy(i,i,:)) &
                           - normal(2)*uy(1,1,:)*cyy(i,i,:)
         if (present(extra_rs)) then
            do k=1,size(extra_rs,dim=2)
               extra_rs(i,k,:) = extra_rs(i,k,:) - &
                              normal(1)*(ux(1,k+1,:)*cxx(i,i,:) + &
                                         uy(1,k+1,:)*cxy(i,i,:)) &
                            - normal(2)*uy(1,k+1,:)*cyy(i,i,:)
            end do
         endif
      endif
      if (itype(i) == MIXED) then
         do j=1,ss
            rs(i,:) = rs(i,:) - c(i,j,:)*u(j,1,:)
         end do
         if (present(extra_rs)) then
            do k=1,size(extra_rs,dim=2)
               do j=1,ss
                  extra_rs(i,k,:) = extra_rs(i,k,:) - c(i,j,:)*u(j,k+1,:)
               end do
            end do
         endif
      endif
   end do

endif

end subroutine resid_natural_bconds

!          ------------------
subroutine old_compute_normal(grid,x,y,elem,edge,normal)
!          ------------------

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

end subroutine old_compute_normal

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
subroutine init_cond_ei(grid,elem,energy,work,energy_est,Linf,L2)
!          ------------

!----------------------------------------------------
! This routine computes an error indicator as the difference between the
! computed solution and the function in iconds.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: elem
real(my_real), intent(out), optional :: energy(:), energy_est(:), Linf(:,:), &
                                        L2(:,:), work
!----------------------------------------------------
! Local variables:

real (my_real), pointer :: qw(:),xq(:),yq(:)
real (my_real), allocatable :: err(:,:,:)
real(my_real) :: xc(3), yc(3), loc_energy(max(1,grid%num_eval))
integer :: nqp, jerr, ss, neigen, astat, comp, eigen, qp
!----------------------------------------------------
! Begin executable code

ss = grid%system_size
neigen = max(1,grid%num_eval)

! determine the quadrature rule

xc = grid%vertex(grid%element(elem)%vertex)%coord%x
yc = grid%vertex(grid%element(elem)%vertex)%coord%y
call quadrature_rule_tri(6,xc,yc,nqp,qw,xq,yq,jerr,stay_in=.true.)

! evaluate the solution at the quadrature points

allocate(err(ss,neigen,nqp),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in init_cond_ei")
   stop
endif
call evaluate_soln_local(grid,xq,yq,elem,(/(comp,comp=1,ss)/), &
                         (/(eigen,eigen=1,neigen)/),err)

! evaluate the error at the quadrature points

do eigen=1,neigen
   do comp=1,ss
      do qp=1,nqp
         err (comp,eigen,qp) = err (comp,eigen,qp) - &
                               iconds(xq(qp),yq(qp),comp,eigen)
      end do
   end do
end do

! we don't have derivatives to do an energy norm, so just use the
! discrete L2 norm

if (present(energy) .or. present(energy_est)) then
   do qp=1,nqp
      do eigen=1,neigen
         loc_energy(eigen) = 0.0_my_real
         do comp=1,ss
            loc_energy(eigen) = loc_energy(eigen) + err(comp,eigen,qp)**2
         end do
      end do
   end do
   if (present(energy)) energy = sqrt(loc_energy)
   if (present(energy_est)) energy_est = sqrt(loc_energy)
endif

! L-infinity norm

if (present(Linf)) then
   Linf = 0.0_my_real
   do eigen=1,neigen
      do comp=1,ss
         do qp=1,nqp
            Linf(comp,eigen) = max(Linf(comp,eigen),abs(err(comp,eigen,qp)))
         end do
      end do
   end do
endif

! L2 norm

if (present(L2)) then
   L2 = 0.0_my_real
   do eigen=1,neigen
      do comp=1,ss
         do qp=1,nqp
            L2(comp,eigen) = L2(comp,eigen) + qw(qp)*err(comp,eigen,qp)**2
         end do
      end do
   end do
   L2 = sqrt(L2)
endif

if (present(work)) then
   if (grid%element(elem)%mate == BOUNDARY) then
      work = ((grid%element(elem)%degree)*(grid%element(elem)%degree+1))/2
   else
      work = grid%element(elem)%degree**2
   endif
endif

deallocate(qw,xq,yq,err)

end subroutine init_cond_ei

!          --------------------
subroutine all_error_indicators(grid,method)
!          --------------------

!----------------------------------------------------
! This routine computes the error indicators for all elements, and the
! global error estimates
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
integer, intent(in) :: method
!----------------------------------------------------
! Local variables:

integer :: lev, elem, mate, eval, comp, nev, ss, nelem, astat
real(my_real) :: energy(max(1,grid%num_eval)), &
                 Linf(grid%system_size,max(1,grid%num_eval)), &
                 L2(grid%system_size,max(1,grid%num_eval)), &
                 diam,xvert(3),yvert(3),area,domain_area
real(my_real), allocatable :: all_energy(:,:),all_work(:),all_energy_est(:,:), &
                              all_Linf(:,:,:),all_L2(:,:,:)
logical :: compatibly_divisible
!----------------------------------------------------
! Begin executable code

ss = grid%system_size
nev = max(1,grid%num_eval)
nelem = size(grid%element)

! For the equilibrated residual error estimate, compute all of the error
! indicators at once

if (method == EQUILIBRATED_RESIDUAL) then
   allocate(all_energy(nev,nelem),all_work(nelem),all_energy_est(nev,nelem), &
            all_Linf(ss,nev,nelem),all_L2(ss,nev,nelem),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in all_error_indicators")
      stop 
   endif
   call equilibrated_residual_ei_all(grid,all_energy,all_work,all_energy_est, &
                                     all_Linf,all_L2)
endif

domain_area = (grid%boundbox_max%x-grid%boundbox_min%x)*(grid%boundbox_max%y-grid%boundbox_min%y)

! initialize error estimates to 0

grid%errest_energy = 0.0_my_real
grid%errest_Linf = 0.0_my_real
grid%errest_L2 = 0.0_my_real
grid%errest_eigenvalue = 0.0_my_real

! for each leaf element ...

do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      grid%element_errind(elem,:) = 0.0_my_real
      grid%element(elem)%work = 1.0_my_real
      if (grid%element(elem)%isleaf) then

! identify the mate and whether or not elem is compatibly divisible

         compatibly_divisible = .true.
         if (grid%element(elem)%mate == BOUNDARY) then
            mate = BOUNDARY
         else
            mate = hash_decode_key(grid%element(elem)%mate,grid%elem_hash)
            if (mate == HASH_NOT_FOUND) then
               mate = hash_decode_key(grid%element(elem)%mate/2,grid%elem_hash)
               compatibly_divisible = .false.
            endif
         endif

! compute error indicators

         if (method == EQUILIBRATED_RESIDUAL) then
            grid%element_errind(elem,:) = all_energy(:,elem)
            grid%element(elem)%work = all_work(elem)
            energy = all_energy_est(:,elem)
            Linf = all_Linf(:,:,elem)
            L2 = all_L2(:,:,elem)
         else
            call error_indicator(grid,elem,method, &
                                 grid%element_errind(elem,:), &
                                 grid%element(elem)%work,energy,Linf,L2)
         endif

! if I own the element, add to error estimate.

         if (grid%element(elem)%iown) then
            grid%errest_energy = grid%errest_energy + energy**2
            grid%errest_Linf = &
               max(grid%errest_Linf,reshape(Linf,(/grid%nsoln/)))
            grid%errest_L2 = grid%errest_L2 + reshape(L2**2,(/grid%nsoln/))
            if (grid%num_eval > 0) then
               xvert = grid%vertex(grid%element(elem)%vertex)%coord%x
               yvert = grid%vertex(grid%element(elem)%vertex)%coord%y
               area = abs(xvert(1)*(yvert(2)-yvert(3)) + &
                          xvert(2)*(yvert(3)-yvert(1)) + &
                          xvert(3)*(yvert(1)-yvert(2))) / 2
               diam = sqrt(area/domain_area)
               do eval=1,grid%num_eval
                  do comp=1,grid%system_size
                     grid%errest_eigenvalue(eval)=grid%errest_eigenvalue(eval) &
                      + ((diam**(grid%element(elem)%degree-1))*L2(comp,eval))**2
                  end do
               end do
            endif
         endif

      endif

      elem = grid%element(elem)%next
   end do
end do

grid%errest_energy = sqrt(grid%errest_energy)
grid%errest_L2 = sqrt(grid%errest_L2)
grid%errest_eigenvalue = sqrt(grid%errest_eigenvalue)
! remove the arbitrary constant that was applied to the L2 error indicators
if (method == EXPLICIT_ERRIND) then
   grid%errest_eigenvalue = grid%errest_eigenvalue*10
endif
if (method == REFSOLN_ERREST) then
   grid%errest_energy = grid%refsoln_errest
endif


if (method == EQUILIBRATED_RESIDUAL) then
   deallocate(all_energy,all_work,all_energy_est,all_Linf,all_L2)
endif

grid%errind_up2date = .true.

end subroutine all_error_indicators

!          --------------
subroutine error_estimate(grid,procs,method,soln,errest_energy, &
                          errest_Linf,errest_L2,errest_eigenvalue)
!          --------------

!----------------------------------------------------
! This routine returns requested estimates of the error over this partition.
! If soln is present then it indicates which solution to use. Default is
! maximum over all solutions.  For errest_eigenvalue, soln is which
! eigenvalue to estimate.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(proc_info), intent(in) :: procs
integer, intent(in), optional :: method
integer, intent(in), optional :: soln
real(my_real), intent(out), optional :: errest_energy, errest_Linf, errest_L2, &
                                        errest_eigenvalue
!----------------------------------------------------
! Local variables:

integer :: lo, hi1, hi2
!----------------------------------------------------
! Begin executable code

if (my_proc(procs) == MASTER) then
   if (present(errest_energy)) errest_energy = 0.0_my_real
   if (present(errest_Linf)) errest_Linf = 0.0_my_real
   if (present(errest_L2)) errest_L2 = 0.0_my_real
   if (present(errest_eigenvalue)) errest_eigenvalue = 0.0_my_real
   return
endif

! for REFSOLN_ERREST, the error estimate is stored in the grid data
! structure.  The places where this routine is called without method present
! (hbmg and krylov solvers to get error estimate for termination criteria)
! are OK to use a different method.
! The norm is actually the H1 seminorm, but we'll return it in energy anyway.

if (present(method)) then
   if (method == REFSOLN_ERREST) then
      if (present(errest_energy)) errest_energy = grid%refsoln_errest
      return
   endif
endif

! compute error indicators if they are not up to date

if (.not. grid%errind_up2date) then
   if (present(method)) then
      call all_error_indicators(grid,method)
   else
      call all_error_indicators(grid,EXPLICIT_ERRIND)
   endif
endif

if (present(soln)) then
   lo = soln
   hi1 = soln
   hi2 = soln
else
   lo = 1
   hi1 = max(1,grid%num_eval)
   hi2 = grid%nsoln
endif

if (present(errest_energy)) errest_energy=maxval(grid%errest_energy(lo:hi1))
if (present(errest_Linf)) errest_Linf = maxval(grid%errest_Linf(lo:hi2))
if (present(errest_L2)) errest_L2 = maxval(grid%errest_L2(lo:hi2))
if (present(errest_eigenvalue)) errest_eigenvalue = &
                                   maxval(grid%errest_eigenvalue(lo:hi1))

end subroutine error_estimate

!          -----------------
subroutine interior_residual(grid,x,y,elem,comp,eigen,resid)
!          -----------------

!----------------------------------------------------
! This routine computes the interior residual
!    r = f + dp/dx du/dx + p d2u/dx2 + dq/dy du/dy + q d2u/dy2 - ru
!        - cx du/dx - cy du/dy - cxy d2u/dxdy
! at the given points (x,y) which lie in the interior of element elem where
! u is the approximate solution.  comp and eigen are any subset of the
! components and eigenfunctions.  Subscripts of resid are comp, eigen, point.
! dp/dx and dq/dy are estimated with h = diam(elem)/10.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
real(my_real), intent(in) :: x(:), y(:)
integer, intent(in) :: elem, comp(:), eigen(:)
real(my_real), intent(out) :: resid(:,:,:)
!----------------------------------------------------
! Local variables:

real(my_real) ::   u(size(comp),size(eigen),size(x)), &
                  ux(size(comp),size(eigen),size(x)), &
                  uy(size(comp),size(eigen),size(x)), &
                 uxy(size(comp),size(eigen),size(x)), &
                 uxx(size(comp),size(eigen),size(x)), &
                 uyy(size(comp),size(eigen),size(x))
real(my_real) :: cxx(grid%system_size,grid%system_size), &
                 dxx(grid%system_size,grid%system_size), &
                 cxy(grid%system_size,grid%system_size), &
                 dxy(grid%system_size,grid%system_size), &
                 cyy(grid%system_size,grid%system_size), &
                 dyy(grid%system_size,grid%system_size), &
                  cx(grid%system_size,grid%system_size), &
                  cy(grid%system_size,grid%system_size), &
                   c(grid%system_size,grid%system_size), &
                  rs(grid%system_size), &
                  rse(grid%system_size,max(1,grid%num_eval))
real(my_real) :: diam, hx, hy, xvert(3), yvert(3)
integer :: i, j, p
!----------------------------------------------------
! Begin executable code

! evaluate the solution and derivatives at the points (x,y)

call evaluate_soln_local(grid,x,y,elem,comp,eigen,u,ux,uy,uxx,uyy,uxy)

! determine the diameter of the element as the longest edge length

xvert = grid%vertex(grid%element(elem)%vertex)%coord%x
yvert = grid%vertex(grid%element(elem)%vertex)%coord%y
diam =          sqrt((xvert(2)-xvert(1))**2 + (yvert(2)-yvert(1))**2)
diam = max(diam,sqrt((xvert(3)-xvert(2))**2 + (yvert(3)-yvert(2))**2))
diam = max(diam,sqrt((xvert(1)-xvert(3))**2 + (yvert(1)-yvert(3))**2))

! for each point

do p=1,size(x)

! evaluate cxx and cxy at the point displaced in x

   if (abs(x(p)-grid%boundbox_min%x) > abs(x(p)-grid%boundbox_max%x)) then
      hx = -diam/10.0_my_real
   else
      hx =  diam/10.0_my_real
   endif
   call pdecoefs(x(p)+hx,y(p),dxx,dxy,cyy,cx,cy,c,rs)

! evaluate cyy at the point displaced in y

   if (abs(y(p)-grid%boundbox_min%y) > abs(y(p)-grid%boundbox_max%y)) then
      hy = -diam/10.0_my_real
   else
      hy =  diam/10.0_my_real
   endif
   call pdecoefs(x(p),y(p)+hy,cxx,cxy,dyy,cx,cy,c,rs)

! evaluate the pde data at (x,y)

   call pdecoefs(x(p),y(p),cxx,cxy,cyy,cx,cy,c,rs)

! replace dxx, dxy and dyy with the approximate derivative of p, cxy and q

   dxx = (dxx-cxx)/hx
   dxy = (dxy-cxy)/hx
   dyy = (dyy-cyy)/hy

! change rhs for eigenvalue problems

   if(grid%num_eval > 0) then
      do i=1,grid%num_eval
         rse(:,i) = rs * grid%eigenvalue(i) * u(:,i,p)
      end do
   else
      rse(:,1) = rs
   endif

! evaluate the residual

   do i=1,size(comp)
      do j=1,size(eigen)
         resid(i,j,p) = rse(i,j) + &
              sum(dxx(i,:)*ux(:,j,p) + cxx(i,:)*uxx(:,j,p) + &
                  dyy(i,:)*uy(:,j,p) + cyy(i,:)*uyy(:,j,p) + &
                  dxy(i,:)*uy(:,j,p) + cxy(i,:)*uxy(:,j,p) - &
                  cx(i,:)*ux(:,j,p) - cy(i,:)*uy(:,j,p) - &
                  c(i,:)*u(:,j,p) )
      end do
   end do

end do

end subroutine interior_residual

!          -----------------
subroutine boundary_residual(grid,x,y,elem,edge,comp,eigen,resid)
!          -----------------

!----------------------------------------------------
! This routine computes the boundary residual
!       { g-du/dn-cu   if edge is a Natural or Mixed boundary edge
!   R = { -[du/dn]     if edge is an interior edge or Periodic boundary edge
!       { 0            if edge is a Dirichlet boundary edge
! at the given points (x,y) which lie on the given edge of element elem, where
! u is the approximate solution, n is the unit outward (from element elem)
! normal, and [] is the jump discontinuity across the edge.
! du/dn is actually more complicated, involving the pde and bc coefficients
! comp and eigen are any subset of the components and eigenfunctions.
! Subscripts of resid are comp, eigen, point.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
real(my_real), intent(in) :: x(:),y(:)
integer, intent(in) :: elem, edge, comp(:), eigen(:)
real(my_real), intent(out) :: resid(:,:,:)
!----------------------------------------------------
! Local variables:

real(my_real), allocatable :: u(:,:,:),ux(:,:,:),uy(:,:,:),cxx(:,:,:), &
                              cxy(:,:,:),cyy(:,:,:),cx(:,:,:),cy(:,:,:), &
                              c(:,:,:),rs(:,:),bcc(:,:,:),bcrs(:,:)
real(my_real) :: normal(2)
integer :: itype(grid%system_size)
integer :: i, j, p, iedge, neigh, neighbors(3)
logical :: first_loop, first_loop2
logical, save :: first_call = .true. ! TEMP for periodic problem
!----------------------------------------------------
! Begin executable code

allocate(u(size(comp),size(eigen),size(x)), &
         ux(size(comp),size(eigen),size(x)), &
         uy(size(comp),size(eigen),size(x)), &
         cxx(grid%system_size,grid%system_size,size(x)), &
         cxy(grid%system_size,grid%system_size,size(x)), &
         cyy(grid%system_size,grid%system_size,size(x)), &
         cx (grid%system_size,grid%system_size,size(x)), &
         cy (grid%system_size,grid%system_size,size(x)), &
         c  (grid%system_size,grid%system_size,size(x)), &
         rs (grid%system_size,size(x)), &
         bcc(grid%system_size,grid%system_size,size(x)), &
         bcrs(grid%system_size,size(x)))

first_loop = .true.
first_loop2 = .true.

do i=1,size(comp)

   if (grid%edge_type(edge,i) == DIRICHLET) then

      resid(i,:,:) = 0.0_my_real

   elseif (grid%edge_type(edge,i) == NATURAL .or. &
           grid%edge_type(edge,i) == MIXED .or. &
           grid%edge_type(edge,i) == INTERIOR .or. &
           grid%edge_type(edge,i) == PERIODIC_MASTER .or. &
           grid%edge_type(edge,i) == PERIODIC_SLAVE) then

      if (first_loop) then
         first_loop = .false.

! determine which edge of elem is edge

         if (grid%element(elem)%edge(1) == edge) then
            iedge = 1
         elseif (grid%element(elem)%edge(2) == edge) then
            iedge = 2
         elseif (grid%element(elem)%edge(3) == edge) then
            iedge = 3
         else
            ierr = PHAML_INTERNAL_ERROR
            call fatal("in boundary_residual, edge does not belong to elem")
            stop
         endif

! compute the outward unit normal

         call compute_normal(grid,elem,iedge,normal)

! evaluate the coefficients of the pde to use in the natural b.c.

         do p=1,size(x)
            call pdecoefs(x(p),y(p),cxx(:,:,p),cxy(:,:,p),cyy(:,:,p),cx(:,:,p),&
                          cy(:,:,p),c(:,:,p),rs(:,p))
         end do

! evaluate the boundary conditions

         if (grid%edge_type(edge,i) == NATURAL .or. &
              grid%edge_type(edge,i) == MIXED) then
            do p=1,size(x)
               call bconds(x(p),y(p),grid%edge(edge)%bmark,itype,bcc(:,:,p), &
                           bcrs(:,p))
            end do
         endif

! compute the solution and derivatives in this element

         call evaluate_soln_local(grid,x,y,elem,comp,eigen,u=u,ux=ux,uy=uy)

      endif ! first_loop

      if (grid%edge_type(edge,i) == NATURAL .or. &
           grid%edge_type(edge,i) == MIXED) then

! compute the residual of the boundary conditions

         do j=1,size(eigen)
            do p=1,size(x)
               resid(i,j,p) = bcrs(i,p) - &
                              sum(normal(1)*(cxx(i,:,p)*ux(:,j,p) + &
                                             cxy(i,:,p)*uy(:,j,p)) + &
                                  normal(2)*cyy(i,:,p)*uy(:,j,p) + &
                                  bcc(i,:,p)*u(:,j,p) )
            end do
         end do

      else

! on an interior or periodic edge, compute the first term of the jump

         do j=1,size(eigen)
            do p=1,size(x)
               resid(i,j,p) = sum(-normal(1)*(cxx(i,:,p)*ux(:,j,p) &
                                             +cxy(i,:,p)*uy(:,j,p)) &
                                  +normal(2)*cyy(i,:,p)*uy(:,j,p) )
            end do
         end do

! determine the neighbor triangle

         if (first_loop2) then
            first_loop2 = .false.

            if (grid%edge_type(edge,i) == INTERIOR) then
               neighbors = get_neighbors(elem,grid)
               neigh = neighbors(iedge)

            elseif (grid%edge_type(edge,i) == PERIODIC_SLAVE .or. &
                    grid%edge_type(edge,i) == PERIODIC_MASTER) then
! TEMP for periodic problem
               if (first_call) then
                  call warning("Have not determined how to find the neighbor across a periodic boundary in", &
                               "boundary_residual.  Omitting jump term.")
                  first_call = .false.
               endif
               resid(i,:,:) = 0.0_my_real
               cycle
      
            endif

! compute the solution derivatives in the neighbor

            call evaluate_soln_local(grid,x,y,neigh,comp,eigen,ux=ux,uy=uy)

         endif ! first_loop2

! subtract the normal derivative to get the jump

         do j=1,size(eigen)
            do p=1,size(x)
               resid(i,j,p) = resid(i,j,p) + &
                       sum(normal(1)*(cxx(i,:,p)*ux(:,j,p) &
                                     +cxy(i,:,p)*uy(:,j,p)) &
                         - normal(2)*cyy(i,:,p)*uy(:,j,p) )
            end do
         end do

      endif

   else
      ierr = PHAML_INTERNAL_ERROR
      call fatal("unrecognized edge type in boundary_residual", &
                 intlist=(/grid%edge_type(edge,i)/))
      stop

   endif

end do ! i=component

deallocate(u,ux,uy,cxx,cxy,cyy,cx,cy,c,rs,bcc,bcrs)

end subroutine boundary_residual

!          --------------
subroutine compute_normal(grid,elem,iedge,normal)
!          --------------

!----------------------------------------------------
! This routine calculates the outward unit normal for the iedge'th edge in
! element elem.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: elem, iedge
real(my_real), intent(out) :: normal(2)
!----------------------------------------------------
! Local variables:

real(my_real) :: x1, x2, x3, y1, y2, y3
logical :: clockwise
!----------------------------------------------------
! Begin executable code

! determine if the vertices are clockwise by checking the sign of the cross
! product of two edges

x1 = grid%vertex(grid%element(elem)%vertex(1))%coord%x
x2 = grid%vertex(grid%element(elem)%vertex(2))%coord%x
x3 = grid%vertex(grid%element(elem)%vertex(3))%coord%x
y1 = grid%vertex(grid%element(elem)%vertex(1))%coord%y
y2 = grid%vertex(grid%element(elem)%vertex(2))%coord%y
y3 = grid%vertex(grid%element(elem)%vertex(3))%coord%y

clockwise = (x2-x1)*(y3-y2)-(y2-y1)*(x3-x2) < 0.0_my_real

! change the use of x1 etc to be the two endpoints of the edge

if (iedge == 1) then
   x1 = x2
   y1 = y2
   x2 = x3
   y2 = y3
elseif (iedge == 2) then
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

!          ----------------------------
subroutine equilibrated_residual_ei_all(grid,energy,work,energy_est,Linf,L2, &
                                        multi_p)
!          ----------------------------

!----------------------------------------------------
! This routine computes the equilibrated residual error estimate for all
! elements.  See Ainsworth & Oden, Chapter 6.
! If multi_p is present it computes multiple energy error estimates using
! spaces of degree p+1, p+2, ... p+size(multi_p,dim=1)
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
real(my_real), intent(out), optional :: energy(:,:), work(:), energy_est(:,:), &
                                        Linf(:,:,:), L2(:,:,:), multi_p(:,:,:)
!----------------------------------------------------
! Local variables:

integer :: lev, vert, K, ss, nev, astat
!----------------------------------------------------
! Begin executable code

ss = grid%system_size
nev = max(1,grid%num_eval)

! allocate space for the flux moments

allocate(mu_K(size(grid%element),3,2,ss,nev),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in equilibrated_residual_ei_all")
   stop
endif

! For each vertex ...

do lev = 1,grid%nlev
   vert = grid%head_level_vert(lev)
   do while (vert /= END_OF_LIST)
      call eq_resid_worker1(grid,vert)
      vert = grid%vertex(vert)%next
   end do
end do

! for each element K ...

do lev=1,grid%nlev
   K = grid%head_level_elem(lev)
   do while (K /= END_OF_LIST)
      if (.not. grid%element(K)%isleaf) then
         K = grid%element(K)%next
         cycle
      endif
      call eq_resid_worker2_all(grid,K,energy,work,energy_est,Linf,L2,multi_p)
      K = grid%element(K)%next
   end do
end do

deallocate(mu_K)

end subroutine equilibrated_residual_ei_all

!          ----------------------------
subroutine equilibrated_residual_ei_one(grid,K,energy,work,energy_est,Linf,L2, &
                                        multi_p)
!          ----------------------------

!----------------------------------------------------
! This routine computes the equilibrated residual error estimate for element K.
! See Ainsworth & Oden, Chapter 6.
! If multi_p is present it computes multiple energy error estimates using
! spaces of degree p+1, p+2, ... p+size(multi_p,dim=1)
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: K
real(my_real), intent(out), optional :: energy(:), work, energy_est(:), &
                                        Linf(:,:), L2(:,:), multi_p(:,:)
!----------------------------------------------------
! Local variables:

integer :: i, vert, ss, nev, astat
!----------------------------------------------------
! Begin executable code

ss = grid%system_size
nev = max(1,grid%num_eval)

! allocate space for the flux moments

allocate(mu_K(size(grid%element),3,2,ss,nev),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in equilibrated_residual_ei_one")
   stop
endif

do i=1,3
   vert = grid%element(K)%vertex(i)
   call eq_resid_worker1(grid,vert)
end do

call eq_resid_worker2_one(grid,K,energy,work,energy_est,Linf,L2,multi_p)

deallocate(mu_K)

end subroutine equilibrated_residual_ei_one

!          ----------------
subroutine eq_resid_worker1(grid,vert)
!          ----------------

!----------------------------------------------------
! This routine performs the work of the equilibrated residual error
! indicator inside the "for each vertex" loop.  It ultimately ends with
! mu_K defined around the given vertex.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: vert
!----------------------------------------------------
! Local variables:

integer :: i1, i2, K, num_patch_elements, astat, vert_local_index, &
           edge_local_index, vert_local_edge_index, gamma, ss, nev, comp, &
           patch_type(grid%system_size)
integer, pointer :: patch_element(:), patch_edge(:)
real(my_real), allocatable :: sigma(:,:,:), mu_K_tilde(:,:,:,:)
!----------------------------------------------------
! Begin executable code

ss = grid%system_size
nev = max(1,grid%num_eval)

! Identify the patch of elements around the vertex

call build_patch(grid,vert,patch_element,patch_edge,num_patch_elements, &
                 patch_type)

! allocate memory for sigma, the solution of the topological matrix, which
! initially holds the right hand side tilde{Delta}_K(theta_n), and mu_K_tilde

allocate(sigma(num_patch_elements,nev,ss), &
         mu_K_tilde(num_patch_elements,2,ss,nev), stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in eq_resid_worker1")
   stop
endif

! For each patch element, K ...

do i1 = 1,num_patch_elements
   K = patch_element(i1)

! determine the local index of this vertex in this element

   if (grid%element(K)%vertex(1) == vert) then
      vert_local_index = 1
   elseif (grid%element(K)%vertex(2) == vert) then
      vert_local_index = 2
   elseif (grid%element(K)%vertex(3) == vert) then
      vert_local_index = 3
   else
      ierr = PHAML_INTERNAL_ERROR
      call fatal("didn't find vertex in patch element in eq_resid_worker1")
   endif

! set the right hand side for the topology matrix system
! tilde{Delta}_K(theta_n) = B_K(u_X,theta_n) - (f,theta_n)_K -
! integral_partialK < frac{partial u_X}{partial n_K} > theta_n ds
! where
!        { 1/2 n_K dot {(del u_X)_K + (del u_X)_K'} on partialK cap partialK'
! <..> = { n_K dot (del u_x)_K                      on partialK cap Gamma_D
!        { g                                        on partialK cap Gamma_N

   sigma(i1,:,:) = Delta_K_tilde(grid,K,vert_local_index)

! next patch element

end do

! solve for the sigma vector

call solve_topmat(sigma,num_patch_elements,patch_type)

! For each patch edge, gamma ...

do i1 = 1,num_patch_elements+1
   gamma = patch_edge(i1)

! For elements on each side of gamma

   do i2 = 1,2

! no element before the first edge or after the last edge

      if ((i1==1 .and. i2==1) .or. &
          (i1==num_patch_elements+1 .and. i2==2)) then
         cycle
      endif

      K = patch_element(i1+i2-2)

! Compute the approximate flux moments
! tilde{mu}_{K,n}^gamma = integral_gamma theta_n n_K dot del u_X|_K ds

      mu_K_tilde(i1+i2-2,3-i2,:,:) = compute_mu_K_tilde(grid,K,gamma,vert)

! next element; next patch edge

   end do
end do

! For each patch element, K ...

do i1 = 1,num_patch_elements
   K = patch_element(i1)

! For each element edge, gamma ...

   do i2 = 1,2
      gamma = patch_edge(i1+i2-1)

! determine the local index of the edge in the element

      if (grid%element(K)%edge(1) == gamma) then
         edge_local_index = 1
      elseif (grid%element(K)%edge(2) == gamma) then
         edge_local_index = 2
      elseif (grid%element(K)%edge(3) == gamma) then
         edge_local_index = 3
      else
         ierr = PHAML_INTERNAL_ERROR
         call fatal("didn't find edge in element in eq_resid_worker1")
         stop
      endif

! determine the local index of the vertex in the edge

      if (grid%edge(gamma)%vertex(1) == vert) then
         vert_local_edge_index = 1
      elseif (grid%edge(gamma)%vertex(2) == vert) then
         vert_local_edge_index = 2
      else
         ierr = PHAML_INTERNAL_ERROR
         call fatal("didn't find vertex in edge in eq_resid_worker1")
         stop
      endif

! set the flux moments
!                 { 1/2 (sigma_{K,n} - sigma_{K',n} + tilde{mu}_{K,n}^gamma -
!                 { tilde{mu}_{K',n}^gamma), gamma = partialK cap partialK'
!                 {
! mu_{K,n}^gamma ={ integral_gamma g theta_n ds, gamma = partialK cap Gamma_N
!                 {
!                 { sigma_{K,n} + tilde{mu}_{K,n}^gamma, gamma = partialK cap Gamma_D

! different components of a system of equations may have different types

      do comp=1,ss

! first element, first edge

         if (i1==1 .and. i2==1) then

!    Neumann
            if (patch_type(comp)==NEUMANN_NEUMANN) then

               mu_K(K,edge_local_index,vert_local_edge_index,comp,:) = &
                  integral_g_theta_n(grid,vert,gamma,K,comp)
!    Dirichlet
            elseif (patch_type(comp)==DIRICHLET_NEUMANN .or. &
                    patch_type(comp)==DIRICHLET_DIRICHLET) then

               mu_K(K,edge_local_index,vert_local_edge_index,comp,:) = &
                  sigma(i1,:,comp) + mu_K_tilde(i1,i2,comp,:)
!    interior
            else

               mu_K(K,edge_local_index,vert_local_edge_index,comp,:) = &
                  (sigma(i1,:,comp)-sigma(num_patch_elements,:,comp) + &
                   mu_K_tilde(i1,i2,comp,:) - &
                   mu_K_tilde(num_patch_elements,2,comp,:))/2
            endif

! last element, second edge

         elseif (i1==num_patch_elements .and. i2==2) then

!    Neumann
            if (patch_type(comp)==NEUMANN_NEUMANN .or. &
                patch_type(comp)==DIRICHLET_NEUMANN) then

               mu_K(K,edge_local_index,vert_local_edge_index,comp,:) = &
                  integral_g_theta_n(grid,vert,gamma,K,comp)

!    Dirichlet
            elseif (patch_type(comp)==DIRICHLET_DIRICHLET) then

               mu_K(K,edge_local_index,vert_local_edge_index,comp,:) = &
                  sigma(i1,:,comp) + mu_K_tilde(i1,i2,comp,:)

!    interior
            else

               mu_K(K,edge_local_index,vert_local_edge_index,comp,:) = &
                  (sigma(i1,:,comp)-sigma(1,:,comp) + &
                   mu_K_tilde(i1,i2,comp,:)-mu_K_tilde(1,1,comp,:))/2
            endif

! edge between element and previous element

         elseif (i2==1) then

            mu_K(K,edge_local_index,vert_local_edge_index,comp,:) = &
               (sigma(i1,:,comp)-sigma(i1-1,:,comp) + &
                mu_K_tilde(i1,i2,comp,:)-mu_K_tilde(i1-1,2,comp,:))/2

! edge between element and next element

         else

            mu_K(K,edge_local_index,vert_local_edge_index,comp,:) = &
               (sigma(i1,:,comp)-sigma(i1+1,:,comp) + &
                mu_K_tilde(i1,i2,comp,:)-mu_K_tilde(i1+1,1,comp,:))/2

         endif

      end do ! next component
   end do ! next edge
end do ! next patch element

! deallocate the patch, sigma and mu_K_tilde for this vertex

deallocate(patch_element, patch_edge, sigma, mu_K_tilde)

end subroutine eq_resid_worker1

!          --------------------
subroutine eq_resid_worker2_all(grid,K,energy,work,energy_est,Linf,L2,multi_p)
!          --------------------

!----------------------------------------------------
! This routine performs the work of the equilibrated residual error
! indicator inside the "for each element" loop.  It ultimately ends with
! defining the error indicators and estimates for element K.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in), target :: grid
integer, intent(in) :: K
real(my_real), intent(out), optional :: energy(:,:), work(:), energy_est(:,:), &
                                        Linf(:,:,:), L2(:,:,:), multi_p(:,:,:)
!----------------------------------------------------
! Local variables:

integer :: i, j, deg, degree(4), edge_type(3,grid%system_size), d, &
           bmark(3), n, info, nev, ss, astat, rank, lwork, inc_deg, qorder
real(my_real), allocatable :: matrix(:,:), rhs(:,:), factors(:,:), Ax(:,:), &
                              temp(:), S(:), rwork(:), err(:,:)
real(my_real) :: xvert(3), yvert(3), rmin, area
!----------------------------------------------------
! Begin executable code

ss = grid%system_size
nev = max(1,grid%num_eval)

! set up the local Neumann problem
! B_K(phi_K,v) = (f,v)_K - B_K(u_X,v) + integral_{partialK} g_K v ds all v in V_K

xvert(1) = grid%vertex(grid%element(K)%vertex(1))%coord%x
xvert(2) = grid%vertex(grid%element(K)%vertex(2))%coord%x
xvert(3) = grid%vertex(grid%element(K)%vertex(3))%coord%x
yvert(1) = grid%vertex(grid%element(K)%vertex(1))%coord%y
yvert(2) = grid%vertex(grid%element(K)%vertex(2))%coord%y
yvert(3) = grid%vertex(grid%element(K)%vertex(3))%coord%y

degree(1) = grid%edge(grid%element(K)%edge(1))%degree
degree(2) = grid%edge(grid%element(K)%edge(2))%degree
degree(3) = grid%edge(grid%element(K)%edge(3))%degree
degree(4) = grid%element(K)%degree
if (present(multi_p)) then
   inc_deg = size(multi_p,dim=1)
else
   inc_deg = 1
endif
deg = maxval(degree) + inc_deg
degree = deg

edge_type = NATURAL
bmark = 0
rmin = 0.0_my_real

n = ss*(3*deg + ((deg-2)*(deg-1))/2)
allocate(matrix(n,n),factors(n,n),rhs(n,nev),err(n,nev),Ax(n,nev),temp(n), &
         stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in eq_resid_worker2")
   stop
endif

hold_elem = K
hold_grid => grid

qorder = 0
if (nev > 1) then
   call elemental_matrix(grid, xvert, yvert, degree, edge_type, bmark, &
                         ss, qorder, .true., K, rmin, "p", "a", matrix, &
                         rhs(:,1), equil_resid_bconds,extra_rhs=rhs(:,2:))
else
   call elemental_matrix(grid, xvert, yvert, degree, edge_type, bmark, &
                         ss, qorder, .true., K, rmin, "p", "a", matrix, &
                         rhs(:,1), equil_resid_bconds)
endif

! if multi_p is present, compute error estimates for degree+size(multi_p),
! degree+size(multi_p)-1,...,degree+1 by solving the linear system,
! computing sqrt(e^T A e), removing the rows and columns corresponding to
! degree+size(multi_p) from A, solving the linear system again,
! computing sqrt(e^T A e) again, etc.
! At the end of this, A, rhs and err will be for the degree+1 system, which
! is correct for the single error estimates.

do i=1,inc_deg

   factors = matrix
   err = rhs

! use the LAPACK solver for the minimum norm solution of the least squares
! problem

   lwork = 2*(3*n+max(2*n,nev))
   allocate(S(n),rwork(lwork),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in eq_resid_worker2")
      stop
   endif
   if (my_real == kind(1.0)) then
      call sgelss(n,n,nev,factors,n,err,n,S,1.0e-5_my_real,rank,rwork, &
                  lwork,info)
   elseif (my_real == kind(1.0d0)) then
      call dgelss(n,n,nev,factors,n,err,n,S,1.0e-12_my_real,rank,rwork, &
                  lwork,info)
   else
      ierr = PHAML_INTERNAL_ERROR
      call fatal("my_real is neither single nor double precision. Can't call LAPACK routines in eq_resid_worker2")
      stop
   endif
   if (info /= 0) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("dgelss failed in eq_resid_worker2")
      stop
   endif
   deallocate(S,rwork)

   if (present(multi_p)) then

! sqrt(e^T A e)

      if (my_real == kind(1.0)) then
         call sgemm("N","N",n,nev,n,1.0_my_real,matrix,n,err,n,0.0_my_real,Ax,n)
      else
         call dgemm("N","N",n,nev,n,1.0_my_real,matrix,n,err,n,0.0_my_real,Ax,n)
      endif
      do j=1,nev
            multi_p(inc_deg-i+1,j,K) = sqrt(abs(dot_product(Ax(:,j),err(:,j))))
      end do

! decrease the system by one degree

      if (i /= inc_deg) then
         d = deg-i+1

! edge 1 equations of degree d

         do j=(d+1)*ss+1,(d+1)*ss+ss
            matrix(:,j) = 0.0_my_real
            matrix(j,:) = 0.0_my_real
            matrix(j,j) = 1.0_my_real
            rhs(j,:)    = 0.0_my_real
         end do

! edge 2 equations of degree d

         do j=(d+deg)*ss+1,(d+deg)*ss+ss
            matrix(:,j) = 0.0_my_real
            matrix(j,:) = 0.0_my_real
            matrix(j,j) = 1.0_my_real
            rhs(j,:)    = 0.0_my_real
         end do

! edge 3 equations of degree d

         do j=(d+2*deg-1)*ss+1,(d+2*deg-1)*ss+ss
            matrix(:,j) = 0.0_my_real
            matrix(j,:) = 0.0_my_real
            matrix(j,j) = 1.0_my_real
            rhs(j,:)    = 0.0_my_real
         end do

! face equations of degree d

         do j=(3*deg+(d-2)*(d-3)/2)*ss+1,(3*deg+(d-1)*(d-2)/2)*ss
            matrix(:,j) = 0.0_my_real
            matrix(j,:) = 0.0_my_real
            matrix(j,j) = 1.0_my_real
            rhs(j,:)    = 0.0_my_real
         end do

      endif
   endif
end do

! set error indicators and estimates

! energy error indicator is sqrt(e^T A e)

if (my_real == kind(1.0)) then
   call sgemm("N","N",n,nev,n,1.0_my_real,matrix,n,err,n,0.0_my_real,Ax,n)
else
   call dgemm("N","N",n,nev,n,1.0_my_real,matrix,n,err,n,0.0_my_real,Ax,n)
endif

do i=1,nev

   if (present(energy)) &
      energy(i,K) = sqrt(abs(dot_product(Ax(:,i),err(:,i))))
   if (present(energy_est)) &
      energy_est(i,K) = sqrt(abs(dot_product(Ax(:,i),err(:,i))))

end do

! use midpoint value for L infinity and L2 norms

if (present(Linf) .or. present(L2)) then
   call p_hier_basis_func((xvert(1)+xvert(2)+xvert(3))/3, &
                          (yvert(1)+yvert(2)+yvert(3))/3, &
                          xvert,yvert,degree,"a",temp)
endif

if (present(Linf)) then
   do i=1,nev
      do j=1,ss
         Linf(j,i,K) = abs(sum(err(j:n:ss,i)*temp(1:n/ss)))
      end do
   end do
endif

if (present(L2)) then
   area = abs(xvert(1)*(yvert(2)-yvert(3)) + &
              xvert(2)*(yvert(3)-yvert(1)) + &
              xvert(3)*(yvert(1)-yvert(2))) / 2
   do i=1,nev
      do j=1,ss
         L2(j,i,K) = sqrt(area)*abs(sum(err(j:n:ss,i)*temp(1:n/ss)))
      end do
   end do
endif

if (present(work)) then
   if (grid%element(K)%mate == BOUNDARY) then
      work(K) = ((grid%element(K)%degree)*(grid%element(K)%degree+1))/2
   else
      work(K) = grid%element(K)%degree**2
   endif
endif

! free memory

deallocate(matrix,factors,rhs,err,Ax,temp)

end subroutine eq_resid_worker2_all

!          --------------------
subroutine eq_resid_worker2_one(grid,K,energy,work,energy_est,Linf,L2,multi_p)
!          --------------------

!----------------------------------------------------
! This routine performs the work of the equilibrated residual error
! indicator inside the "for each element" loop.  It ultimately ends with
! defining the error indicators and estimates for element K.
! This version is for one element at a time; the only difference is the rank
! of the output arrays.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in), target :: grid
integer, intent(in) :: K
real(my_real), intent(out), optional :: energy(:), work, energy_est(:), &
                                        Linf(:,:), L2(:,:), multi_p(:,:)
!----------------------------------------------------
! Local variables:

integer :: i, j, deg, degree(4), edge_type(3,grid%system_size), d, &
           bmark(3), n, info, nev, ss, astat, rank, lwork, inc_deg, qorder
real(my_real), allocatable :: matrix(:,:), rhs(:,:), factors(:,:), Ax(:,:), &
                              temp(:), S(:), rwork(:), err(:,:)
real(my_real) :: xvert(3), yvert(3), rmin, area
!----------------------------------------------------
! Begin executable code

ss = grid%system_size
nev = max(1,grid%num_eval)

! set up the local Neumann problem
! B_K(phi_K,v) = (f,v)_K - B_K(u_X,v) + integral_{partialK} g_K v ds all v in V_K

xvert(1) = grid%vertex(grid%element(K)%vertex(1))%coord%x
xvert(2) = grid%vertex(grid%element(K)%vertex(2))%coord%x
xvert(3) = grid%vertex(grid%element(K)%vertex(3))%coord%x
yvert(1) = grid%vertex(grid%element(K)%vertex(1))%coord%y
yvert(2) = grid%vertex(grid%element(K)%vertex(2))%coord%y
yvert(3) = grid%vertex(grid%element(K)%vertex(3))%coord%y

degree(1) = grid%edge(grid%element(K)%edge(1))%degree
degree(2) = grid%edge(grid%element(K)%edge(2))%degree
degree(3) = grid%edge(grid%element(K)%edge(3))%degree
degree(4) = grid%element(K)%degree
if (present(multi_p)) then
   inc_deg = size(multi_p,dim=1)
else
   inc_deg = 1
endif
deg = maxval(degree) + inc_deg
degree = deg

edge_type = NATURAL
bmark = 0
rmin = 0.0_my_real

n = ss*(3*deg + ((deg-2)*(deg-1))/2)
allocate(matrix(n,n),factors(n,n),rhs(n,nev),err(n,nev),Ax(n,nev),temp(n), &
         stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in eq_resid_worker2")
   stop
endif

hold_elem = K
hold_grid => grid

qorder = 0
if (nev > 1) then
   call elemental_matrix(grid, xvert, yvert, degree, edge_type, bmark, &
                         ss, qorder, .true., K, rmin, "p", "a", matrix, &
                         rhs(:,1), equil_resid_bconds,extra_rhs=rhs(:,2:))
else
   call elemental_matrix(grid, xvert, yvert, degree, edge_type, bmark, &
                         ss, qorder, .true., K, rmin, "p", "a", matrix, &
                         rhs(:,1), equil_resid_bconds)
endif

! if multi_p is present, compute error estimates for degree+size(multi_p),
! degree+size(multi_p)-1,...,degree+1 by solving the linear system,
! computing sqrt(e^T A e), removing the rows and columns corresponding to
! degree+size(multi_p) from A, solving the linear system again,
! computing sqrt(e^T A e) again, etc.
! At the end of this, A, rhs and err will be for the degree+1 system, which
! is correct for the single error estimates.

do i=1,inc_deg

   factors = matrix
   err = rhs

! use the LAPACK solver for the minimum norm solution of the least squares
! problem

   lwork = 2*(3*n+max(2*n,nev))
   allocate(S(n),rwork(lwork),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in eq_resid_worker2")
      stop
   endif
   if (my_real == kind(1.0)) then
      call sgelss(n,n,nev,factors,n,err,n,S,1.0e-5_my_real,rank,rwork, &
                  lwork,info)
   elseif (my_real == kind(1.0d0)) then
      call dgelss(n,n,nev,factors,n,err,n,S,1.0e-12_my_real,rank,rwork, &
                  lwork,info)
   else
      ierr = PHAML_INTERNAL_ERROR
      call fatal("my_real is neither single nor double precision. Can't call LAPACK routines in eq_resid_worker2")
      stop
   endif
   if (info /= 0) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("dgelss failed in eq_resid_worker2")
      stop
   endif
   deallocate(S,rwork)

   if (present(multi_p)) then

! sqrt(e^T A e)

      if (my_real == kind(1.0)) then
         call sgemm("N","N",n,nev,n,1.0_my_real,matrix,n,err,n,0.0_my_real,Ax,n)
      else
         call dgemm("N","N",n,nev,n,1.0_my_real,matrix,n,err,n,0.0_my_real,Ax,n)
      endif
      do j=1,nev
            multi_p(inc_deg-i+1,j) = sqrt(abs(dot_product(Ax(:,j),err(:,j))))
      end do

! decrease the system by one degree

      if (i /= inc_deg) then
         d = deg-i+1

! edge 1 equations of degree d

         do j=(d+1)*ss+1,(d+1)*ss+ss
            matrix(:,j) = 0.0_my_real
            matrix(j,:) = 0.0_my_real
            matrix(j,j) = 1.0_my_real
            rhs(j,:)    = 0.0_my_real
         end do

! edge 2 equations of degree d

         do j=(d+deg)*ss+1,(d+deg)*ss+ss
            matrix(:,j) = 0.0_my_real
            matrix(j,:) = 0.0_my_real
            matrix(j,j) = 1.0_my_real
            rhs(j,:)    = 0.0_my_real
         end do

! edge 3 equations of degree d

         do j=(d+2*deg-1)*ss+1,(d+2*deg-1)*ss+ss
            matrix(:,j) = 0.0_my_real
            matrix(j,:) = 0.0_my_real
            matrix(j,j) = 1.0_my_real
            rhs(j,:)    = 0.0_my_real
         end do

! face equations of degree d

         do j=(3*deg+(d-2)*(d-3)/2)*ss+1,(3*deg+(d-1)*(d-2)/2)*ss
            matrix(:,j) = 0.0_my_real
            matrix(j,:) = 0.0_my_real
            matrix(j,j) = 1.0_my_real
            rhs(j,:)    = 0.0_my_real
         end do

      endif
   endif
end do

! set error indicators and estimates

! energy error indicator is sqrt(e^T A e)

if (my_real == kind(1.0)) then
   call sgemm("N","N",n,nev,n,1.0_my_real,matrix,n,err,n,0.0_my_real,Ax,n)
else
   call dgemm("N","N",n,nev,n,1.0_my_real,matrix,n,err,n,0.0_my_real,Ax,n)
endif

do i=1,nev

   if (present(energy)) &
      energy(i) = sqrt(abs(dot_product(Ax(:,i),err(:,i))))
   if (present(energy_est)) &
      energy_est(i) = sqrt(abs(dot_product(Ax(:,i),err(:,i))))

end do

! use midpoint value for L infinity and L2 norms

if (present(Linf) .or. present(L2)) then
   call p_hier_basis_func((xvert(1)+xvert(2)+xvert(3))/3, &
                          (yvert(1)+yvert(2)+yvert(3))/3, &
                          xvert,yvert,degree,"a",temp)
endif

if (present(Linf)) then
   do i=1,nev
      do j=1,ss
         Linf(j,i) = abs(sum(err(j:n:ss,i)*temp(1:n/ss)))
      end do
   end do
endif

if (present(L2)) then
   area = abs(xvert(1)*(yvert(2)-yvert(3)) + &
              xvert(2)*(yvert(3)-yvert(1)) + &
              xvert(3)*(yvert(1)-yvert(2))) / 2
   do i=1,nev
      do j=1,ss
         L2(j,i) = sqrt(area)*abs(sum(err(j:n:ss,i)*temp(1:n/ss)))
      end do
   end do
endif

if (present(work)) then
   if (grid%element(K)%mate == BOUNDARY) then
      work = ((grid%element(K)%degree)*(grid%element(K)%degree+1))/2
   else
      work = grid%element(K)%degree**2
   endif
endif

! free memory

deallocate(matrix,factors,rhs,err,Ax,temp)

end subroutine eq_resid_worker2_one

!          ------------
subroutine solve_topmat(rhs,n,ptype)
!          ------------

!----------------------------------------------------
! This routine solves the topology matrix systems for the equilibrated
! residual error estimator.  rhs is the right hand side on input and the
! solution on output.  n is the size of the sytem.  ptype is the patch_type.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(inout) :: rhs(:,:,:)
integer, intent(in) :: n, ptype(:)
!----------------------------------------------------
! Local variables:

type topmat_type
   real(my_real), pointer :: A(:,:)
end type topmat_type

! RESTRICTION no more than 20 elements sharing a vertex
integer, parameter :: maxelem = 20
type(topmat_type), save :: topmat(maxelem,4)
logical, save :: first_call = .true.

integer :: i, j, astat, info, lwork, nrhs, ss, comp, type, rank
real(my_real), allocatable :: T(:,:), rwork(:), S(:)
!----------------------------------------------------
! Begin executable code

ss = size(rhs,dim=3)
nrhs = size(rhs,dim=2)

! on the first call, nullify all the topological matrices so we can test
! them for allocation on subsequent calls

if (first_call) then
   do j=1,4
      do i=1,maxelem
         nullify(topmat(i,j)%A)
      end do
   end do
   first_call = .false.
endif

! check input

if (n > maxelem) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("too many elements surrounding a vertex in solve_topmat")
   stop
endif
if (size(rhs,dim=1) < n) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("size of rhs not big enough for n in solve_topmat")
   stop
endif

! do each component separately, as they may have different types

do comp=1,ss
   type = ptype(comp)

! on first call for this size and type, allocate and set the matrix

   if (.not. associated(topmat(n,type)%A)) then

! allocations

      allocate(topmat(n,type)%A(n,n),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in solve_topmat")
         stop
      endif

! set the topology matrix

      topmat(n,type)%A = 0
      do i=1,n
         if (i /= 1) topmat(n,type)%A(i,i-1) = -0.5_my_real
         topmat(n,type)%A(i,i) = 1
         if (i /= n) topmat(n,type)%A(i,i+1) = -0.5_my_real
      end do
      select case (type)
      case (INTERIOR_VERTEX)
         topmat(n,type)%A(1,n) = -0.5_my_real
         topmat(n,type)%A(n,1) = -0.5_my_real
      case (NEUMANN_NEUMANN)
         topmat(n,type)%A(1,1) = 0.5_my_real
         topmat(n,type)%A(n,n) = 0.5_my_real
      case (DIRICHLET_NEUMANN)
         topmat(n,type)%A(n,n) = 0.5_my_real
         topmat(n,type)%A(1,1) = 1.5_my_real
      case (DIRICHLET_DIRICHLET)
         topmat(n,type)%A(1,1) = 1.5_my_real
         topmat(n,type)%A(n,n) = 1.5_my_real
      end select

! factor A if it has a Dirichlet type

      if (type==DIRICHLET_NEUMANN .or. type==DIRICHLET_DIRICHLET) then
         if (my_real == kind(1.0)) then
            call spotrf("U",n,topmat(n,type)%A,n,info)
         elseif (my_real == kind(1.0d0)) then
            call dpotrf("U",n,topmat(n,type)%A,n,info)
         else
            ierr = PHAML_INTERNAL_ERROR
            call fatal("real is not single or double precision.  cannot call LAPACK in solve_topmat")
            stop
         endif
      else
         info = 0
      endif

      if (info /= 0) then
         ierr = PHAML_INTERNAL_ERROR
         call fatal("dpotrf failed in solve_topmat")
         stop
      endif

   endif ! first call for this type and size

   select case (type)

   case(INTERIOR_VERTEX,NEUMANN_NEUMANN)

! non-Dirichlet types are singular.  Solve for the linear least squares
! solution.  dgelss does not reuse the factors

      lwork = 2*(3*n+max(2*n,nrhs))
      allocate(T(n,n),S(n),rwork(lwork),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in solve_topmat")
         stop
      endif

      T = topmat(n,type)%A

      if (my_real == kind(1.0)) then
         call sgelss(n,n,nrhs,T,n,rhs(:,:,comp),n,S,1.0e-5_my_real,rank,rwork, &
                     lwork,info)
      elseif (my_real == kind(1.0d0)) then
         call dgelss(n,n,nrhs,T,n,rhs(:,:,comp),n,S,1.0e-12_my_real,rank,rwork,&
                     lwork,info)
      else
         ierr = PHAML_INTERNAL_ERROR
         call fatal("my_real is neither single nor double precision. Can't call LAPACK routines in solve_topmat")
         stop
      endif
      if (info /= 0) then
         ierr = PHAML_INTERNAL_ERROR
         call fatal("dgelss failed in solve_topmat")
         stop
      endif
      deallocate(T,S,rwork)

   case(DIRICHLET_NEUMANN,DIRICHLET_DIRICHLET)

! Dirichlet types are nonsingular.  Solve using factored matrix.

      if (my_real == kind(1.0)) then
         call spotrs("U",n,nrhs,topmat(n,type)%A,n,rhs(:,:,comp),n,info)
      else
         call dpotrs("U",n,nrhs,topmat(n,type)%A,n,rhs(:,:,comp),n,info)
      endif

      if (info /= 0) then
         ierr = PHAML_INTERNAL_ERROR
         call fatal("dpotrs failed in solve_topmat")
         stop
      endif

   case default
      ierr = PHAML_INTERNAL_ERROR
      call fatal("bad value for type in solve_topmat")
      stop

   end select

end do ! comp

end subroutine solve_topmat

!          -----------
subroutine build_patch(grid,vert,patch_element,patch_edge,num_patch_elements, &
                       patch_type)
!          -----------

!----------------------------------------------------
! This routine builds a patch consisting of the elements around vertex vert.
! patch_element is allocated in this routine and should be deallocated when
! no longer needed.  Sequential elements in the array are adjacent elements
! in the grid.  If vert is on the boundary, the first and last elements
! are on the boundary; if one is a Dirichlet boundary and the other is a
! Neumann boundary, the Dirichlet is the first element.  If vert is interior,
! the first element is the associated element and the last element is adjacent
! to the first element.
! patch_edge(i) is the edge between patch_element(i-1) and patch_element(i).
! If the vertex is interior, patch_edge(1) and patch_edge(num_patch_elements+1)
! are both the edge between patch_element(num_patch_elements) and
! patch_element(1).  If the vertex is on the boundary, patch_edge(1) and
! patch_edge(num_patch_elements+1) are the boundary edges.
! num_patch_elements is returned as the number of elements in the patch.
! patch type is the type of patch
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: vert
integer, pointer :: patch_element(:), patch_edge(:)
integer, intent(out) :: num_patch_elements, patch_type(:)
!----------------------------------------------------
! Local variables:

integer, parameter :: max_patch_elements = 50
integer :: elements(max_patch_elements), edges(max_patch_elements+1),neigh(3)
integer :: i, boundary_loc, bedge1, bedge2, comp
logical :: found_second_boundary, first_loop
!----------------------------------------------------
! Begin executable code

! initializations

boundary_loc = -1 ! gives the position in elements of the first boundary element

! start with the associated element

num_patch_elements = 1
elements(num_patch_elements) = grid%vertex(vert)%assoc_elem

! get the neighboring elements

neigh = get_neighbors(elements(num_patch_elements),grid)

! see if this element is on the boundary

do i = 1,3
   if (neigh(i) == BOUNDARY) then
      if (grid%edge(grid%element(elements(num_patch_elements))%edge(i))%vertex(1) == vert .or. &
          grid%edge(grid%element(elements(num_patch_elements))%edge(i))%vertex(2) == vert) then
         boundary_loc = num_patch_elements
         bedge1 = grid%element(elements(num_patch_elements))%edge(i)
         exit
      endif
   endif
end do

! find a neighbor that contains vert

do i=1,3
   if (neigh(i) == BOUNDARY) cycle
   if (grid%element(neigh(i))%vertex(1) == vert .or. &
       grid%element(neigh(i))%vertex(2) == vert .or. &
       grid%element(neigh(i))%vertex(3) == vert) then
      num_patch_elements = num_patch_elements + 1
      elements(num_patch_elements) = neigh(i)
      edges(num_patch_elements) = grid%element(elements(num_patch_elements-1))%edge(i)
      exit
   endif
end do

! if we didn't find one, then the vertex must be between two boundary edges
! on the same element

if (num_patch_elements == 1) then
   allocate(patch_element(1),patch_edge(2))
   patch_element(1) = elements(1)
   patch_edge(1) = bedge1
   do i = 1,3
      if (neigh(i) == BOUNDARY .and. &
          grid%element(elements(1))%edge(i) /= bedge1) then
         if (grid%edge(grid%element(elements(1))%edge(i))%vertex(1) == vert .or. &
             grid%edge(grid%element(elements(1))%edge(i))%vertex(2) == vert) then
            patch_edge(2) = grid%element(elements(1))%edge(i)
            exit
         endif
      endif
   end do
   do comp=1,grid%system_size
      if (grid%edge_type(patch_edge(1),comp) == DIRICHLET .and. &
          grid%edge_type(patch_edge(2),comp) == DIRICHLET) then
         patch_type(comp) = DIRICHLET_DIRICHLET
      elseif (grid%edge_type(patch_edge(1),comp) == DIRICHLET) then
         patch_type(comp) = DIRICHLET_NEUMANN
      elseif (grid%edge_type(patch_edge(2),comp) == DIRICHLET) then
         patch_edge(1:2) = patch_edge(2:1:-1)
         patch_type(comp) = DIRICHLET_NEUMANN
      else
         patch_type(comp) = NEUMANN_NEUMANN
      endif
   end do
   return
endif

! repeat until we hit a boundary or return to the starting element

do

   if (boundary_loc /= -1) exit
   if (elements(num_patch_elements) == elements(1)) then
      edges(1) = edges(num_patch_elements)
      num_patch_elements = num_patch_elements - 1
      exit
   endif

! get the neighboring elements

   neigh = get_neighbors(elements(num_patch_elements),grid)

! find one that contains vert and isn't the previous element

   do i = 1,3
      if (neigh(i) == elements(num_patch_elements-1)) cycle
      if (neigh(i) == BOUNDARY) then
         if (grid%edge(grid%element(elements(num_patch_elements))%edge(i))%vertex(1) == vert .or. &
             grid%edge(grid%element(elements(num_patch_elements))%edge(i))%vertex(2) == vert) then
            boundary_loc = num_patch_elements
            bedge1 = grid%element(elements(num_patch_elements))%edge(i)
            exit
         endif
         cycle
      endif
      if (grid%element(neigh(i))%vertex(1) == vert .or. &
          grid%element(neigh(i))%vertex(2) == vert .or. &
          grid%element(neigh(i))%vertex(3) == vert) then
         num_patch_elements = num_patch_elements + 1
         if (num_patch_elements > max_patch_elements) then
            ierr = PHAML_INTERNAL_ERROR
            call fatal("too many elements in patch in build_patch")
            stop
         endif
         elements(num_patch_elements) = neigh(i)
         edges(num_patch_elements) = grid%element(elements(num_patch_elements-1))%edge(i)
         exit
      endif
   end do

end do

! if a boundary was found, go to the other side of the first element and
! repeat until the other boundary is found

if (boundary_loc /= -1) then

   found_second_boundary = .false.

! if the first element is not the boundary ...

   if (boundary_loc /= 1) then

! get the neighboring elements of the first element

      neigh = get_neighbors(elements(1),grid)

! find one that contains vert and isn't the second element

      do i = 1,3
         if (neigh(i) == elements(2)) cycle
         if (neigh(i) == BOUNDARY) then
            if (grid%edge(grid%element(elements(1))%edge(i))%vertex(1) == vert .or. &
                grid%edge(grid%element(elements(1))%edge(i))%vertex(2) == vert) then
               found_second_boundary = .true.
               bedge2 = grid%element(elements(1))%edge(i)
               exit
            endif
            cycle
         endif
         if (grid%element(neigh(i))%vertex(1) == vert .or. &
             grid%element(neigh(i))%vertex(2) == vert .or. &
             grid%element(neigh(i))%vertex(3) == vert) then
            num_patch_elements = num_patch_elements + 1
            if (num_patch_elements > max_patch_elements) then
               ierr = PHAML_INTERNAL_ERROR
               call fatal("too many elements in patch in build_patch")
               stop
            endif
            elements(num_patch_elements) = neigh(i)
            edges(num_patch_elements) = grid%element(elements(1))%edge(i)
            exit
         endif
      end do

   endif

! repeat until the second boundary is found

   first_loop = .true.
   do

      if (found_second_boundary) exit

! get the neighboring elements

      neigh = get_neighbors(elements(num_patch_elements),grid)

! find one that contains vert and isn't the previous element

      do i = 1,3
         if (first_loop) then
            if (neigh(i) == elements(1)) cycle
         else
            if (neigh(i) == elements(num_patch_elements-1)) cycle
         endif
         if (neigh(i) == BOUNDARY) then
            if (grid%edge(grid%element(elements(num_patch_elements))%edge(i))%vertex(1) == vert .or. &
                grid%edge(grid%element(elements(num_patch_elements))%edge(i))%vertex(2) == vert) then
               found_second_boundary = .true.
               bedge2 = grid%element(elements(num_patch_elements))%edge(i)
               exit
            endif
            cycle
         endif
         if (grid%element(neigh(i))%vertex(1) == vert .or. &
             grid%element(neigh(i))%vertex(2) == vert .or. &
             grid%element(neigh(i))%vertex(3) == vert) then
            num_patch_elements = num_patch_elements + 1
            if (num_patch_elements > max_patch_elements) then
               ierr = PHAML_INTERNAL_ERROR
               call fatal("too many elements in patch in build_patch")
               stop
            endif
            elements(num_patch_elements) = neigh(i)
            edges(num_patch_elements) = grid%element(elements(num_patch_elements-1))%edge(i)
            exit
         endif
      end do
      first_loop = .false.

   end do

endif

! copy the elements into patch_element

allocate(patch_element(num_patch_elements),patch_edge(num_patch_elements+1))

! if there was a boundary, go from the boundary back to the beginning, and
! then after the boundary until the end

if (boundary_loc /= -1) then
   do i=boundary_loc,1,-1
      patch_element(boundary_loc-i+1) = elements(i)
      if (i /= 1) patch_edge(boundary_loc-i+2) = edges(i)
   end do
   do i=boundary_loc+1,num_patch_elements
      patch_element(i) = elements(i)
      patch_edge(i) = edges(i)
   end do
   patch_edge(1) = bedge1
   patch_edge(num_patch_elements+1) = bedge2

! otherwise just copy it

else
   do i=1,num_patch_elements
      patch_element(i) = elements(i)
      patch_edge(i) = edges(i)
   end do
   patch_edge(num_patch_elements+1) = edges(num_patch_elements+1)
endif

! determine the patch type

do comp=1,grid%system_size
   if (boundary_loc == -1) then
      patch_type(comp) = INTERIOR_VERTEX
   elseif (grid%edge_type(bedge1,comp) == DIRICHLET .and. &
           grid%edge_type(bedge2,comp) == DIRICHLET) then
      patch_type(comp) = DIRICHLET_DIRICHLET
   elseif (grid%edge_type(bedge1,comp) == DIRICHLET) then
      patch_type(comp) = DIRICHLET_NEUMANN
   elseif (grid%edge_type(bedge2,comp) == DIRICHLET) then
      patch_element(1:num_patch_elements) = patch_element(num_patch_elements:1:-1)
      patch_edge(1:num_patch_elements+1) = patch_edge(num_patch_elements+1:1:-1)
      patch_type(comp) = DIRICHLET_NEUMANN
   else
      patch_type(comp) = NEUMANN_NEUMANN
   endif
end do

end subroutine build_patch

!        -------------
function Delta_K_tilde(grid,K,vert_local_index)
!        -------------

!----------------------------------------------------
! This routine computes the right hand side for the topology matrix system
! tilde{Delta}_K(theta_n) = B_K(u_X,theta_n) - (f,theta_n)_K -
! integral_partialK < frac{partial u_X}{partial n_K} > theta_n ds
! where
!        { 1/2 n_K dot {(del u_X)_K + (del u_X)_K'} on partialK cap partialK'
! <..> = { n_K dot (del u_x)_K                      on partialK cap Gamma_D
!        { g                                        on partialK cap Gamma_N

!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: K, vert_local_index
real(my_real) :: Delta_K_tilde(max(1,grid%num_eval),grid%system_size)
!----------------------------------------------------
! Local variables:

integer :: i, j, k2, side, edge, qorder, nqpoints, jerr, astat, ss, nev
real(my_real) :: xvert(5), yvert(5)
real(my_real), pointer :: xquad(:),yquad(:), qweights(:)
real(my_real), allocatable :: cxx(:,:,:), cxy(:,:,:), cyy(:,:,:), cx(:,:,:), &
                              cy(:,:,:), c(:,:,:), rs(:,:)
real(my_real), allocatable :: basis(:,:), basisx(:,:), basisy(:,:), &
                              u(:,:,:), ux(:,:,:), uy(:,:,:), &
                              dudn(:,:,:)
!----------------------------------------------------
! Begin executable code

ss = grid%system_size
nev = max(1,grid%num_eval)

! get the coordinates of the vertices of element K, including an extra
! copy of the first and second to make it easy to identify the edge endpoints

do i=1,3
   xvert(i) = grid%vertex(grid%element(K)%vertex(i))%coord%x
   yvert(i) = grid%vertex(grid%element(K)%vertex(i))%coord%y
end do
xvert(4) = xvert(1)
yvert(4) = yvert(1)
xvert(5) = xvert(2)
yvert(5) = yvert(2)

! theta_n is the piecewise linear basis at vertex vert, phi, u_X is the current
! approximate solution, U.  So B_K(u_X,theta_n) - (f,theta_n)_K is
! integral_K p*U_x*phi_x + q*U_y*phi*y + (r*U-f)*phi dx dy

! get a quadrature rule

qorder = min(1+grid%element(K)%degree,MAX_QUAD_ORDER_TRI)
call quadrature_rule_tri(qorder,xvert(1:3),yvert(1:3),nqpoints,qweights, &
                         xquad,yquad,jerr,stay_in=.true.)

! evaluate the pde coefficients at the quadrature points

allocate(cxx(ss,ss,nqpoints),cxy(ss,ss,nqpoints),cyy(ss,ss,nqpoints), &
         cx(ss,ss,nqpoints),cy(ss,ss,nqpoints),c(ss,ss,nqpoints), &
         rs(ss,nqpoints),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in Delta_K_tilde")
   stop
endif

do i = 1,nqpoints
   call pdecoefs(xquad(i),yquad(i),cxx(:,:,i),cxy(:,:,i),cyy(:,:,i), &
                 cx(:,:,i),cy(:,:,i),c(:,:,i),rs(:,i))
end do

! evaluate the basis functions at the quadrature points

allocate(basis(3,nqpoints),basisx(3,nqpoints),basisy(3,nqpoints),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in Delta_K_tilde")
   stop
endif

call p_hier_basis_func(xquad,yquad,xvert(1:3),yvert(1:3),(/1,1,1,1/),"a", &
                       basis,basisx,basisy)

! evaluate the solution at the quadrature points

allocate(u(ss,nev,nqpoints),ux(ss,nev,nqpoints),uy(ss,nev,nqpoints), &
         stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in Delta_K_tilde")
   stop
endif

call evaluate_soln_local(grid,xquad,yquad,K,(/(i,i=1,ss)/),(/(i,i=1,nev)/), &
                         u,ux,uy)

! evaluate integral_K p*U_x*phi_x + q*U_y*phi_y + (r*U-f)*phi dx dy

do j=1,ss
   do i=1,nev
      if (grid%num_eval <= 0) then
         Delta_K_tilde(i,j) = -sum(qweights*rs(j,:)*basis(vert_local_index,:))
         do k2=1,ss
            Delta_K_tilde(i,j) = Delta_K_tilde(i,j) + sum(qweights * &
                          (cxx(j,k2,:)*ux(k2,i,:)*basisx(vert_local_index,:) + &
                           cyy(j,k2,:)*uy(k2,i,:)*basisy(vert_local_index,:) + &
                           cxy(j,k2,:)*uy(k2,i,:)*basisx(vert_local_index,:) + &
                            cx(j,k2,:)*ux(k2,i,:)*basis (vert_local_index,:) + &
                            cy(j,k2,:)*uy(k2,i,:)*basis (vert_local_index,:) + &
                             c(j,k2,:)* u(k2,i,:)*basis (vert_local_index,:)))
         end do
      else
         Delta_K_tilde(i,j) = -sum(qweights* &
                 rs(j,:)*grid%eigenvalue(i)*u(j,i,:)*basis(vert_local_index,:))
         do k2=1,ss
            Delta_K_tilde(i,j) = Delta_K_tilde(i,j) + sum(qweights * &
                          (cxx(j,k2,:)*ux(k2,i,:)*basisx(vert_local_index,:) + &
                           cyy(j,k2,:)*uy(k2,i,:)*basisy(vert_local_index,:) + &
                           cxy(j,k2,:)*uy(k2,i,:)*basisx(vert_local_index,:) + &
                            cx(j,k2,:)*ux(k2,i,:)*basis (vert_local_index,:) + &
                            cy(j,k2,:)*uy(k2,i,:)*basis (vert_local_index,:) + &
                             c(j,k2,:)* u(k2,i,:)*basis (vert_local_index,:)))
         end do
      endif
   end do
end do

! free memory for triangle quadrature rule, bases, and solution

deallocate(xquad,yquad,qweights,cxx,cxy,cyy,cx,cy,c,rs,basis,basisx,basisy,u, &
           ux,uy)

! compute the integral around the boundary.  Go through the three boundary
! sides, skipping the one that is opposite vert since theta_n is zero there.

do side=1,3
   if (side == vert_local_index) cycle
   edge = grid%element(K)%edge(side)

! get a quadrature rule

   qorder = min(2*grid%edge(edge)%degree,MAX_QUAD_ORDER_LINE)
   call quadrature_rule_line(qorder,xvert(side+1:side+2),yvert(side+1:side+2), &
                             nqpoints,qweights,xquad,yquad,jerr)

! evaluate the <...> function at the quadrature points

   allocate(dudn(ss,nev,nqpoints),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in Delta_K_tilde")
      stop
   endif

   call bracket_dudn(grid,xquad,yquad,K,side,dudn)

! evaluate the basis function at the quadrature points

   allocate(basis(3,nqpoints),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in Delta_K_tilde")
      stop
   endif

   call p_hier_basis_func(xquad,yquad,xvert(1:3),yvert(1:3),(/1,1,1,1/),"a", &
                          basis)

! subtract the integral along this side

   do j=1,ss
      do i=1,nev
         Delta_K_tilde(i,j) = Delta_K_tilde(i,j) - &
                            sum(qweights*dudn(j,i,:)*basis(vert_local_index,:))
      end do
   end do

! free memory for quadrature points, du/dn and basis

   deallocate(xquad,yquad,qweights,dudn,basis)

! next side

end do

end function Delta_K_tilde

!          ------------
subroutine bracket_dudn(grid,x,y,K,side,dudn)
!          ------------

!----------------------------------------------------
! This routine computes
!           { 1/2 n_K dot {(del u_X)_K + (del u_X)_K'} on partialK cap partialK'
! <du/dn> = { n_K dot (del u_x)_K                      on partialK cap Gamma_D
!           { g                                        on partialK cap Gamma_N
! at (x,y) in dudn.  x, y and dudn must all have the same size.
! all (x,y) must lie on the same side of element K.
! side is the local index of side in element K (1, 2 or 3)
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
real(my_real), intent(in) :: x(:), y(:)
integer, intent(in) :: side, K
real(my_real), intent(out) :: dudn(:,:,:)
!----------------------------------------------------
! Local variables:

real(my_real) :: u_K(grid%system_size,max(1,grid%num_eval),size(x)), &
                 ux_K(grid%system_size,max(1,grid%num_eval),size(x)), &
                 ux_Kprime(grid%system_size,max(1,grid%num_eval),size(x)), &
                 uy_K(grid%system_size,max(1,grid%num_eval),size(x)), &
                 uy_Kprime(grid%system_size,max(1,grid%num_eval),size(x)), &
                 normal(2)
real(my_real) :: c(grid%system_size,grid%system_size),rs(grid%system_size)
integer :: itype(grid%system_size)
integer :: i, j, neigh(3), edge, Kprime, ss, nev, comp
!----------------------------------------------------
! Begin executable code

ss = grid%system_size
nev = max(1,grid%num_eval)

edge = grid%element(K)%edge(side)

! evaluate the normal to this side, if needed

if (any(grid%edge_type(edge,:) == INTERIOR) .or. &
    any(grid%edge_type(edge,:) == DIRICHLET)) then
   call compute_normal(grid,K,side,normal)
endif

! evaluate the derivatives of the solution in K at the points, if needed

if (any(grid%edge_type(edge,:) == INTERIOR) .or. &
    any(grid%edge_type(edge,:) == DIRICHLET).or. &
    any(grid%edge_type(edge,:) == MIXED)) then
   call evaluate_soln_local(grid,x,y,K,(/(i,i=1,ss)/),(/(i,i=1,nev)/), &
                            u=u_K,ux=ux_K,uy=uy_K)
endif

! determine the element, Kprime, that shares this side, and which side it
! is in the local index of Kprime, and evaluate the derivatives of the
! solution in Kprime at the points, if needed

if (any(grid%edge_type(edge,:) == INTERIOR)) then
   neigh = get_neighbors(K,grid)
   Kprime = -1
   do i=1,3
      if (neigh(i) == BOUNDARY) cycle
      do j=1,3
         if (grid%element(neigh(i))%edge(j) == edge) then
            Kprime = neigh(i)
            exit
         endif
      end do
   end do

   if (Kprime == -1) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("didn't find neighbor sharing side in bracket_dudn")
      stop
   endif

   call evaluate_soln_local(grid,x,y,Kprime,(/(i,i=1,ss)/),(/(i,i=1,nev)/), &
                            ux=ux_Kprime,uy=uy_Kprime)

endif

! do each component

do comp=1,ss

! Three cases

   select case(grid%edge_type(edge,comp))

   case (INTERIOR)

! compute the function

      do i=1,nev
         dudn(comp,i,:) = (normal(1)*(ux_K(comp,i,:)+ux_Kprime(comp,i,:)) + &
                           normal(2)*(uy_K(comp,i,:)+uy_Kprime(comp,i,:)))/2
      end do

   case (DIRICHLET)

! compute the function

      do i=1,nev
         dudn(comp,i,:) = normal(1)*ux_K(comp,i,:) + normal(2)*uy_K(comp,i,:)
      end do

   case (NATURAL,MIXED)

! evaluate the right hand side of the boundary conditions at this point.
! for mixed, also add c*u which should be part of B_K(u_X,theta_n)

      do i = 1,size(x)
         call bconds(x(i),y(i),grid%edge(edge)%bmark,itype,c,rs)
         dudn(comp,:,i) = rs(comp)
         if (grid%edge_type(edge,comp) == MIXED) then
            dudn(comp,:,i) = dudn(comp,:,i) + c(comp,comp)*u_K(comp,:,i)
         endif
      end do

   case default
      ierr = PHAML_INTERNAL_ERROR
      call fatal("unknown edge type in bracket_dudn")

   end select

end do ! comp

end subroutine bracket_dudn

!        ------------------
function compute_mu_K_tilde(grid,K,gamma,vert)
!        ------------------

!----------------------------------------------------
! This routine computes the approximate flux moments
! tilde{mu}_{K,n}^gamma = integral_gamma theta_n n_K dot del u_X|_K ds
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: K, gamma, vert
real(my_real) :: compute_mu_K_tilde(grid%system_size,max(1,grid%num_eval))
!----------------------------------------------------
! Local variables:

real(my_real) :: xverte(2), yverte(2), xvertt(3), yvertt(3), normal(2)
integer :: i, j, vert_local_index, qorder, nqpoints, jerr, edge_local_index, &
           ss, nev, astat
real(my_real), pointer :: qweights(:), xquad(:), yquad(:)
real(my_real), allocatable :: basis(:,:), ux(:,:,:), uy(:,:,:), dudn(:,:,:)
!----------------------------------------------------
! Begin executable code

ss = grid%system_size
nev = max(1,grid%num_eval)

! get the vertex coordinates

xverte(1) = grid%vertex(grid%edge(gamma)%vertex(1))%coord%x
xverte(2) = grid%vertex(grid%edge(gamma)%vertex(2))%coord%x
yverte(1) = grid%vertex(grid%edge(gamma)%vertex(1))%coord%y
yverte(2) = grid%vertex(grid%edge(gamma)%vertex(2))%coord%y
xvertt(1) = grid%vertex(grid%element(K)%vertex(1))%coord%x
xvertt(2) = grid%vertex(grid%element(K)%vertex(2))%coord%x
xvertt(3) = grid%vertex(grid%element(K)%vertex(3))%coord%x
yvertt(1) = grid%vertex(grid%element(K)%vertex(1))%coord%y
yvertt(2) = grid%vertex(grid%element(K)%vertex(2))%coord%y
yvertt(3) = grid%vertex(grid%element(K)%vertex(3))%coord%y

! identify the local index of the vertex in the element

if (grid%element(K)%vertex(1) == vert) then
   vert_local_index = 1
elseif (grid%element(K)%vertex(2) == vert) then
   vert_local_index = 2
elseif (grid%element(K)%vertex(3) == vert) then
   vert_local_index = 3
else
   ierr = PHAML_INTERNAL_ERROR
   call fatal("didn't find vertex in element in compute_mu_K_tilde")
endif

! determine the local index of the edge in the element

if (grid%element(K)%edge(1) == gamma) then
   edge_local_index = 1
elseif (grid%element(K)%edge(2) == gamma) then
   edge_local_index = 2
elseif (grid%element(K)%edge(3) == gamma) then
   edge_local_index = 3
else
   ierr = PHAML_INTERNAL_ERROR
   call fatal("didn't find edge in element in compute_mu_K_tilde")
   stop
endif

! get a quadrature rule

qorder = min(grid%edge(gamma)%degree+2,MAX_QUAD_ORDER_LINE)
call quadrature_rule_line(qorder,xverte,yverte,nqpoints,qweights, &
                          xquad,yquad,jerr)

! allocate memory based on number of quadrature points

allocate(basis(3,nqpoints),ux(ss,nev,nqpoints),uy(ss,nev,nqpoints), &
         dudn(ss,nev,nqpoints),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in compute_mu_K_tilde")
   stop
endif

! evaluate the basis function at the quadrature points

call p_hier_basis_func(xquad,yquad,xvertt,yvertt,(/1,1,1,1/),"a",basis)

! evaluate the normal to this side

call compute_normal(grid,K,edge_local_index,normal)

! evaluate the derivatives of the solution in K at the quadrature points

call evaluate_soln_local(grid,xquad,yquad,K,(/(i,i=1,ss)/),(/(i,i=1,nev)/), &
                         ux=ux,uy=uy)

! compute the normal derivative

do j=1,ss
   do i=1,nev
      dudn(j,i,:) = normal(1)*ux(j,i,:) + normal(2)*uy(j,i,:)
   end do
end do

! compute the integral

do j=1,ss
   do i=1,nev
      compute_mu_K_tilde(j,i) = &
         sum(qweights*basis(vert_local_index,:)*dudn(j,i,:))
   end do
end do

! free memory

deallocate(qweights,xquad,yquad,basis,ux,uy,dudn)

end function compute_mu_K_tilde

!        ------------------
function integral_g_theta_n(grid,vert,edge,elem,comp)
!        ------------------

!----------------------------------------------------
! This routine computes the integral of the right hand side of component
! comp's natural boundary condition times the linear basis function centered
! at vert, over the given edge.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: vert, edge, elem, comp
real(my_real) :: integral_g_theta_n
!----------------------------------------------------
! Local variables:

real(my_real) :: xverte(2),yverte(2),xvertt(3),yvertt(3)
real(my_real), pointer :: qweights(:), xquad(:), yquad(:)
real(my_real), allocatable :: basis(:,:), g(:)
integer :: i, qorder, jerr, nqpoints, vert_local_index, astat
real(my_real) :: c(grid%system_size,grid%system_size),rs(grid%system_size)
integer :: itype(grid%system_size)
!----------------------------------------------------
! Begin executable code

! get the vertex coordinates

xverte(1) = grid%vertex(grid%edge(edge)%vertex(1))%coord%x
xverte(2) = grid%vertex(grid%edge(edge)%vertex(2))%coord%x
yverte(1) = grid%vertex(grid%edge(edge)%vertex(1))%coord%y
yverte(2) = grid%vertex(grid%edge(edge)%vertex(2))%coord%y
xvertt(1) = grid%vertex(grid%element(elem)%vertex(1))%coord%x
xvertt(2) = grid%vertex(grid%element(elem)%vertex(2))%coord%x
xvertt(3) = grid%vertex(grid%element(elem)%vertex(3))%coord%x
yvertt(1) = grid%vertex(grid%element(elem)%vertex(1))%coord%y
yvertt(2) = grid%vertex(grid%element(elem)%vertex(2))%coord%y
yvertt(3) = grid%vertex(grid%element(elem)%vertex(3))%coord%y

! identify the local index of the vertex in the element

if (grid%element(elem)%vertex(1) == vert) then
   vert_local_index = 1
elseif (grid%element(elem)%vertex(2) == vert) then
   vert_local_index = 2
elseif (grid%element(elem)%vertex(3) == vert) then
   vert_local_index = 3
else
   ierr = PHAML_INTERNAL_ERROR
   call fatal("didn't find vertex in patch element in integral_g_theta_n")
endif

! get a quadrature rule

qorder = min(grid%edge(edge)%degree+2,MAX_QUAD_ORDER_LINE)
call quadrature_rule_line(qorder,xverte,yverte,nqpoints,qweights, &
                          xquad,yquad,jerr)

! allocate memory based on number of quadrature points

allocate(basis(3,nqpoints),g(nqpoints),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in integral_g_theta_n")
   stop
endif

! evaluate the basis function at the quadrature points

call p_hier_basis_func(xquad,yquad,xvertt,yvertt,(/1,1,1,1/),"a",basis)

! evaluate the boundary condition at the quadrature points

do i = 1,nqpoints
   call bconds(xquad(i),yquad(i),grid%edge(edge)%bmark,itype,c,rs)
   g(i) = rs(comp)
end do

! compute integral

integral_g_theta_n = sum(qweights*g*basis(vert_local_index,:))

! free memory

deallocate(qweights,xquad,yquad,basis)

end function integral_g_theta_n

!          ------------------
subroutine equil_resid_bconds(x,y,bmark,itype,c,rs,extra_rs)
!          ------------------

!----------------------------------------------------
! This routine computes the boundary conditions for the equilibrated
! residual error estimator
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x(:),y(:)
integer, intent(in) :: bmark
integer, intent(out) :: itype(:)
real(my_real), intent(out) :: c(:,:,:),rs(:,:)
real(my_real), intent(out), optional :: extra_rs(:,:,:)
!----------------------------------------------------
! Local variables:

integer :: i, j, k, k2, iarr(1), side, edge, l, rr, degree(4), deg, qorder, &
           jerr, nqpoints, ss, nev, edegree(4), astat
real(my_real), pointer :: qweights(:),xquad(:),yquad(:)
real(my_real) :: zeta(3), xvert(5), yvert(5), side_length
real(my_real), allocatable :: cxx(:,:,:), cxy(:,:,:), cyy(:,:,:), cx(:,:,:), &
                              cy(:,:,:), r(:,:,:), f(:,:), basis(:,:), &
                              basisx(:,:), basisy(:,:), u(:,:,:), ux(:,:,:), &
                              uy(:,:,:), alpha(:,:,:)
!----------------------------------------------------
! Begin executable code

ss = hold_grid%system_size
if (present(extra_rs)) then
   nev = 1 + size(extra_rs,dim=2)
else
   nev = 1
endif

itype = NATURAL
c = 0.0_my_real

! determine which side of the triangle these points are on by computing the
! barycentric coordinates and looking for one that is closest to 0

zeta = barycentric(x(1),y(1),hold_elem,hold_grid,no_det=.true.)
iarr = minloc(abs(zeta))
side = iarr(1)
edge = hold_grid%element(hold_elem)%edge(side)

! determine the vertex index in the element for the first (l for left) and
! second (rr) vertices of the side

if (hold_grid%element(hold_elem)%vertex(1) == hold_grid%edge(edge)%vertex(1)) then
   l = 1
elseif (hold_grid%element(hold_elem)%vertex(2) == hold_grid%edge(edge)%vertex(1)) then
   l = 2
else
   l = 3
endif
if (hold_grid%element(hold_elem)%vertex(1) == hold_grid%edge(edge)%vertex(2)) then
   rr = 1
elseif (hold_grid%element(hold_elem)%vertex(2) == hold_grid%edge(edge)%vertex(2)) then
   rr = 2
else
   rr = 3
endif

! vertex coordinates

xvert(1) = hold_grid%vertex(hold_grid%element(hold_elem)%vertex(1))%coord%x
xvert(2) = hold_grid%vertex(hold_grid%element(hold_elem)%vertex(2))%coord%x
xvert(3) = hold_grid%vertex(hold_grid%element(hold_elem)%vertex(3))%coord%x
yvert(1) = hold_grid%vertex(hold_grid%element(hold_elem)%vertex(1))%coord%y
yvert(2) = hold_grid%vertex(hold_grid%element(hold_elem)%vertex(2))%coord%y
yvert(3) = hold_grid%vertex(hold_grid%element(hold_elem)%vertex(3))%coord%y
xvert(4) = xvert(1)
yvert(4) = yvert(1)
xvert(5) = xvert(2)
yvert(5) = yvert(2)

! length of the side the points are on

side_length = sqrt((xvert(l)-xvert(rr))**2 + (yvert(l)-yvert(rr))**2)

! degrees

degree(1) = hold_grid%edge(hold_grid%element(hold_elem)%edge(1))%degree
degree(2) = hold_grid%edge(hold_grid%element(hold_elem)%edge(2))%degree
degree(3) = hold_grid%edge(hold_grid%element(hold_elem)%edge(3))%degree
degree(4) = hold_grid%element(hold_elem)%degree
deg = maxval(degree)

! allocate the solution of the edge mass system

allocate(alpha(deg+1,ss,nev),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in equil_resid_bconds")
   stop
endif

! the first two entries of the right side are from the mu_K from linear bases

alpha(1,:,:) = mu_K(hold_elem,side,1,:,:)
alpha(2,:,:) = mu_K(hold_elem,side,2,:,:)

! other entries are Delta_K(theta_n) for the high order bases on this side,
! i.e. B_K(u_X,theta_n) - (f,theta_n)_K

! get a quadrature rule

qorder = min(2*deg,MAX_QUAD_ORDER_TRI)
call quadrature_rule_tri(qorder,xvert(1:3),yvert(1:3),nqpoints,qweights,xquad,yquad,jerr,stay_in=.true.)

! evaluate the pde coefficients at the quadrature points

allocate(cxx(ss,ss,nqpoints),cxy(ss,ss,nqpoints),cyy(ss,ss,nqpoints), &
         cx(ss,ss,nqpoints),cy(ss,ss,nqpoints),r(ss,ss,nqpoints), &
         f(ss,nqpoints),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in equil_resid_bconds")
   stop
endif

do i = 1,nqpoints
   call pdecoefs(xquad(i),yquad(i),cxx(:,:,i),cxy(:,:,i),cyy(:,:,i), &
                 cx(:,:,i),cy(:,:,i),r(:,:,i),f(:,i))
end do

! evaluate the edge basis functions at the quadrature points

allocate(basis(2+deg,nqpoints),basisx(2+deg,nqpoints),basisy(2+deg,nqpoints), &
         stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in equil_resid_bconds")
   stop
endif

edegree = 0
edegree(side) = deg
call p_hier_basis_func(xquad,yquad,xvert(1:3),yvert(1:3),edegree,"a", &
                       basis,basisx,basisy)

! evaluate the solution at the quadrature points

allocate(u(ss,nev,nqpoints),ux(ss,nev,nqpoints),uy(ss,nev,nqpoints), &
         stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in equil_resid_bconds")
   stop
endif

call evaluate_soln_local(hold_grid,xquad,yquad,hold_elem,(/(i,i=1,ss)/), &
                         (/(i,i=1,nev)/),u,ux,uy)

! evaluate integral_K p*U_x*phi_x + q*U_y*phi_y + (r*U-f)*phi dx dy, a.k.a.
! Delta_K(theta_n) a.k.a. mu_K a.k.a. the right side of the edge mass system

do i=3,deg+1
   do k=1,ss
      do j=1,nev
         if (hold_grid%num_eval <= 0) then
            alpha(i,k,j) = -sum(qweights*f(k,:)*basis(i+1,:))
            do k2=1,ss
               alpha(i,k,j) = alpha(i,k,j) + sum(qweights * &
                                (cxx(k,k2,:)*ux(k2,j,:)*basisx(i+1,:) + &
                                 cyy(k,k2,:)*uy(k2,j,:)*basisy(i+1,:) + &
                                 cxy(k,k2,:)*uy(k2,j,:)*basisx(i+1,:) + &
                                  cx(k,k2,:)*ux(k2,j,:)*basis (i+1,:) + &
                                  cy(k,k2,:)*uy(k2,j,:)*basis (i+1,:) + &
                                   r(k,k2,:)* u(k2,j,:)*basis (i+1,:)))
            end do
         else
            alpha(i,k,j) = -sum(qweights* &
                           f(k,:)*hold_grid%eigenvalue(j)*u(k,j,:)*basis(i+1,:))
            do k2=1,ss
               alpha(i,k,j) = alpha(i,k,j) + sum(qweights * &
                             (cxx(k,k2,:)*ux(k2,j,:)*basisx(i+1,:) + &
                              cyy(k,k2,:)*uy(k2,j,:)*basisy(i+1,:) + &
                              cxy(k,k2,:)*uy(k2,j,:)*basisx(i+1,:) + &
                               cx(k,k2,:)*ux(k2,j,:)*basis (i+1,:) + &
                               cy(k,k2,:)*uy(k2,j,:)*basis (i+1,:) + &
                                r(k,k2,:)* u(k2,j,:)*basis (i+1,:)))
            end do
         endif 
      end do
   end do
end do

! free memory used to get to here

deallocate(qweights,xquad,yquad,cxx,cxy,cyy,cx,cy,r,f,basis,basisx,basisy,u, &
           ux,uy)

! for mixed boundary conditions, also need integral_{R\cupK} bcc*U*phi

do k=1,ss
   if (hold_grid%edge_type(edge,k) == MIXED) then

! get a quadrature rule

      qorder = min(2*hold_grid%edge(edge)%degree,MAX_QUAD_ORDER_LINE)
      call quadrature_rule_line(qorder,xvert(side+1:side+2),yvert(side+1:side+2), &
                                nqpoints,qweights,xquad,yquad,jerr)

! evaluate the boundary conditions at the quadrature points

      allocate(r(ss,ss,nqpoints),f(ss,nqpoints),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in equil_resid_bconds")
         stop
      endif

      do i = 1,nqpoints
         call bconds(xquad(i),yquad(i),hold_grid%edge(edge)%bmark,itype,r(:,:,i),f(:,i))
      end do

! evaluate the basis functions at the quadrature points

      allocate(basis(2+deg,nqpoints),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in Delta_K_tilde")
         stop
      endif

      call p_hier_basis_func(xquad,yquad,xvert(1:3),yvert(1:3),edegree,"a",basis)

! evaluate the solution at the quadrature points

      allocate(u(1,nev,nqpoints),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in equil_resid_bconds")
         stop
      endif

      call evaluate_soln_local(hold_grid,xquad,yquad,hold_elem,(/k/),(/(i,i=1,nev)/),u)

! evaluate integral

      do i=3,deg+1
         do j=1,nev
            alpha(i,k,j) = alpha(i,k,j) + sum(qweights*r(k,k,:)*u(k,j,:)*basis(i+1,:))
         end do
      end do

! free memory

      deallocate(qweights,xquad,yquad,r,f,basis,u)

! next component

   endif
end do

! solve for alpha

if (my_real == kind(0.0)) then
   call spotrs("U",deg+1,ss*nev,edge_mass(deg)%M,deg+1,alpha,deg+1,jerr)
else
   call dpotrs("U",deg+1,ss*nev,edge_mass(deg)%M,deg+1,alpha,deg+1,jerr)
endif

! the mass matrix is scaled by the side length

alpha = alpha/side_length

! evaluate the basis functions at the given points

allocate(basis(deg+2,size(x)))
call p_hier_basis_func(x,y,xvert,yvert,edegree,"a",basis)

! remove the basis of the opposite vertex

basis(side:deg+1,:) = basis(side+1:deg+2,:)

! recover the flux g_K|_gamma = alpha dot basis

do k=1,ss
   do i=1,size(x)
      rs(k,i) = sum(alpha(:,k,1)*basis(:deg+1,i))
   end do
   if (present(extra_rs)) then
      do j=2,nev
         do i=1,size(x)
            extra_rs(k,j-1,i) = sum(alpha(:,k,j)*basis(:,i))
         end do
      end do
   endif
end do

deallocate(basis,alpha)

end subroutine equil_resid_bconds

!          -------------
subroutine set_edge_mass(deg)
!          -------------

!----------------------------------------------------
! This routine computes and factorizes the edge mass matrices used by the
! equilibrated residual error estimator, up to degree deg
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: deg
!----------------------------------------------------
! Local variables:

integer :: i, j, qorder, jerr, astat, nqpoints
real(my_real) :: xvert(3),yvert(3)
real(my_real), pointer :: qweights(:), xquad(:), yquad(:), basis(:,:)
!----------------------------------------------------
! Begin executable code

! see if it has already been set

if (allocated(edge_mass)) then

! if it has already been set to at least degree deg, don't have to do anything

   if (size(edge_mass) >= deg) return

! if it has already been set but not to degree deg, throw it away

   do i=1,size(edge_mass)
      deallocate(edge_mass(i)%M)
   end do
   deallocate(edge_mass)

endif

! allocate the edge mass matrices

allocate(edge_mass(deg),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in set_edge_mass")
   stop
endif

do i=1,deg
   allocate(edge_mass(i)%M(i+1,i+1),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in set_edge_mass")
      stop
   endif
end do

! compute the integrals on the unit line; scale by side length when used

xvert(1) = 0.0_my_real
xvert(2) = 1.0_my_real
xvert(3) = 0.5_my_real
yvert(1) = 0.0_my_real
yvert(2) = 0.0_my_real
yvert(3) = 0.5_my_real

! get a quadrature rule

qorder = min(deg+1,MAX_QUAD_ORDER_LINE)
call quadrature_rule_line(qorder,xvert,yvert,nqpoints,qweights,xquad,yquad,jerr)

! evaluate the basis functions at the quadrature points

allocate(basis(deg+2,nqpoints),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in set_edge_mass")
   stop
endif
call p_hier_basis_func(xquad,yquad,xvert,yvert,(/0,0,deg,0/),"a",basis)

! remove the extraneous linear basis of the third vertex

basis(3:deg+1,:) = basis(4:deg+2,:)

! compute the mass matrix for the highest degree

do i=1,deg+1
   do j=1,i-1
      edge_mass(deg)%M(i,j) = edge_mass(deg)%M(j,i)
   end do
   do j=i,deg+1
      edge_mass(deg)%M(i,j) = sum(qweights*basis(i,:)*basis(j,:))
   end do
end do

! for each degree, extract the submatrix up to that degree and factor it

do i=1,deg
   edge_mass(i)%M(1:i+1,1:i+1) = edge_mass(deg)%M(1:i+1,1:i+1)
   if (my_real == kind(0.0)) then
      call spotrf("U",i+1,edge_mass(i)%M,i+1,jerr)
   elseif (my_real == kind(0.0d0)) then
      call dpotrf("U",i+1,edge_mass(i)%M,i+1,jerr)
   else
      ierr = PHAML_INTERNAL_ERROR
      call fatal("real is not single or double precision.  cannot call LAPACK in set_edge_mass")
      stop
   endif
   if (jerr /= 0) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("dpotrf failed in set_edge_mass")
      stop
   endif
end do

! free memory

deallocate(qweights,xquad,yquad,basis)

end subroutine set_edge_mass

end module error_estimators
