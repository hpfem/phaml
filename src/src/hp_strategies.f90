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

module hp_strategies

!----------------------------------------------------
! This module contains routines that define the hp-adaptive strategies.
!
! communication tags in this module are of the form 32xx
!----------------------------------------------------

!----------------------------------------------------
! Other modules used:

use global
use gridtype_mod
use grid_util
use refine_uniform_mod
use refine_elements
use message_passing
use hash_mod
use quadrature_rules
use basis_functions
use evaluate
use error_estimators
use make_linsys
use linsystype_mod
use linear_system
use nlp_vars
!----------------------------------------------------

implicit none
private
public mark_reftype, mark_reftype_one, set_numref, refine_refsoln, &
       ext_refsoln_errest_edge, ext_refsoln_errest_elem, &
       refine_nlp, nlp_dof, nlp_ddof, nlp_errest, nlp_derrest

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

   function regularity(x,y)
   use global
   real(my_real), intent(in) :: x(3),y(3)
   real(my_real) :: regularity
   end function regularity

end interface

!----------------------------------------------------
! The following defined constants are defined:

!----------------------------------------------------
! The following variables are defined:

real(my_real) :: n3pee(3,1), npq(3)

!----------------------------------------------------

contains

!===========================================================================
! Routines common to many strategies
!===========================================================================

!          ------------
subroutine mark_reftype(grid,refine_control,global_max_errind,reftype,elist)
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
                            .true.)
      elem = grid%element(elem)%next
   end do ! next element
end do ! next level

! second pass, mark all delayed elements in bin 1 without delay

elem = elist%head_errind(1)
do while (elem /= END_OF_LIST)
   if (reftype(elem) == "u") then
      call mark_reftype_one(elem,grid,refine_control,global_max_errind,reftype,&
                            .false.)
   endif
   elem = elist%next_errind(elem)
end do

end subroutine mark_reftype

!          ----------------
subroutine mark_reftype_one(elem,grid,refine_control,global_max_errind,reftype,&
                            delay)
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

! for all strategies, use p if h is maxed out or visa versa, unless we're
! stopping when we reach the maximum

         if (grid%element(elem)%level >= refine_control%max_lev .and. &
             grid%element(elem)%degree >= refine_control%max_deg) then
            reftype(elem) = "n"
         elseif (grid%element(elem)%level >= refine_control%max_lev) then
            if (refine_control%stop_on_maxlev) then
               reftype(elem) = "n"
            else
               reftype(elem) = "p"
            endif
         elseif (grid%element(elem)%degree >= refine_control%max_deg) then
            if (refine_control%stop_on_maxdeg) then
               reftype(elem) = "n"
            else
               reftype(elem) = "h"
            endif
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

!===========================================================================
! T3S, TEXAS 3 STEP
!===========================================================================

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
   N_I = phaml_global_sum(procs,N_I_own,3210)
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
      N_I = phaml_global_sum(procs,N_I_own,3210+i)
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
   N_I = phaml_global_sum(procs,N_I_own,3260)

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

!===========================================================================
! PRIOR2P
!===========================================================================

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

!===========================================================================
! NEXT3P
!===========================================================================

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

!===========================================================================
! TYPEPARAM
!===========================================================================

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

!===========================================================================
! COEF_ROOT
!===========================================================================

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

!===========================================================================
! COEF_DECAY
!===========================================================================

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

!===========================================================================
! common to COEF_ROOT and COEF_DECAY
!===========================================================================

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

!===========================================================================
! Routines common to REFSOLN_EDGE and REFSOLN_ELEM
!===========================================================================

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

!===========================================================================
! REFSOLN_EDGE
!===========================================================================

!          -------------------
subroutine refine_refsoln_edge(grid,procs,refine_control,solver_control, &
                               io_control,still_sequential,init_nvert, &
                               init_nelem,init_dof,loop,balance_what,predictive)
!          -------------------

!----------------------------------------------------
! This routine is the top level refine routine for the REFSOLN_EDGE hp strategy
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

type(grid_type) :: old_grid, ref_soln
type(refine_options) :: loc_refcont
type(io_options) :: loc_iocont
!----------------------------------------------------
! Begin executable code

! TEMP not parallel

   if (num_proc(procs) > 1) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("REFSOLN_EDGE not yet implemented in parallel")
      stop
   endif

! Make a copy of the grid as it currently exists.  Note the copy is in old_grid.

   call copy_grid(grid,old_grid)

! Create the reference solution by performing uniform h and p refinements
! and solving on the fine grid.

   loc_refcont = refine_control
   loc_refcont%error_estimator = HIERARCHICAL_COEFFICIENT
   loc_refcont%reftype = P_UNIFORM
   loc_iocont = io_control
   loc_iocont%print_linsys_when = NEVER
   loc_iocont%print_error_when = NEVER
   call refine_uniform_p(grid,loc_refcont)
   loc_refcont%reftype = H_UNIFORM
   call refine_uniform_h(grid,loc_refcont)
   call solve(grid,procs,loc_iocont,solver_control,still_sequential,.false., &
              no_time=.true.)
   call copy_grid(grid,ref_soln)

! MASTER doesn't participate further; it was only here to monitor solve

if (my_proc(procs) == MASTER) then
   call deallocate_grid(old_grid)
   call deallocate_grid(ref_soln)
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

end subroutine refine_refsoln_edge

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

integer :: lev, elem, nqpoints, jerr, parent, old_elem, grandparent, ss, i, j
real(my_real) :: norm_fine, contribution
real(my_real), pointer :: qweights(:), xquad(:), yquad(:)
real(my_real), allocatable :: ux_fine(:,:,:), uy_fine(:,:,:), ux_old(:,:,:), &
                              uy_old(:,:,:)
type(hash_key) :: parent_gid, grandparent_gid
!----------------------------------------------------
! Begin executable code

   ss = old_grid%system_size

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

            allocate(ux_fine(ss,1,nqpoints),uy_fine(ss,1,nqpoints), &
                     ux_old(ss,1,nqpoints),uy_old(ss,1,nqpoints))
            call evaluate_soln_local(ref_soln,xquad,yquad,elem,(/(i,i=1,ss)/), &
                                     (/1/),ux=ux_fine,uy=uy_fine)
            call evaluate_soln_local(old_grid,xquad,yquad,old_elem, &
                                     (/(i,i=1,ss)/),(/1/),ux=ux_old,uy=uy_old)

! compute the contribution of this element

            contribution = 0.0_my_real
            do i = 1,nqpoints
               do j=1,ss
                  contribution = contribution + &
                              qweights(i)*((ux_old(j,1,i)-ux_fine(j,1,i))**2 + &
                                           (uy_old(j,1,i)-uy_fine(j,1,i))**2)
               end do
            end do

! add contributions to integrals

            compute_refsoln_errest_edge = compute_refsoln_errest_edge + contribution

! same for the norm of the fine solution

            contribution = 0.0_my_real
            do i = 1,nqpoints
               do j=1,ss
                  contribution = contribution + &
                             qweights(i)*(ux_fine(j,1,i)**2 + uy_fine(j,1,i)**2)
               end do
            end do
            norm_fine = norm_fine + contribution

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
loc_refcont%reftype = P_UNIFORM
loc_iocont = io_control
loc_iocont%print_linsys_when = NEVER
loc_iocont%print_error_when = NEVER
call refine_uniform_p(ref_soln,loc_refcont)
loc_refcont%reftype = H_UNIFORM
call refine_uniform_h(ref_soln,loc_refcont)
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
! perform, and performing those refinements, with the edge based version
! of refsoln.
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
                                  edge_rate_max,refine_control)

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
   call weighted_H1_seminorm_diff_sq(ref_soln,old_grid,edge,elem,p,diff_hp)

! Save the norm to return it.
! Note the norm of the difference is in the third component (the first two
! are used for child contributions)

   diff_old = diff_hp(3)

! compute |u_fine - w_p|^2 where w_p is the projection based interpolant
! on the p refined edge

   p(1) = p(1) + 1
   call weighted_H1_seminorm_diff_sq(ref_soln,old_grid,edge,elem,p,diff_p)

! compute the error decrease rate for the p refined edge

   rate_p = diff_hp(3) - diff_p(3)

! for each of the p possible h-refinements, compute the error decrease rate
! and identify the best h refinement

   rate_h = -huge(0.0_my_real)
   do i=1,old_grid%edge(edge)%degree

! compute rate

      p(1) = i
      p(2) = old_grid%edge(edge)%degree + 1 - i
      call weighted_H1_seminorm_diff_sq(ref_soln,old_grid,edge,elem,p,diff_h)
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

!          ----------------------------
subroutine weighted_H1_seminorm_diff_sq(ref_soln,old_grid,edge,elem,p,diff)
!          ----------------------------

!----------------------------------------------------
! This routine computes |u_fine - w|^2 where |.| is the weighted H1 seminorm
! on edge edge in old_grid and w is the projection based interpolation of the
! reference solution in grid onto edge refined as specified in p.  If
! p(2)==0, then the projection is onto a single edge of degree p(1).  Otherwise
! it is onto two edges, formed by bisecting edge, with degrees p(1) and p(2).
! edge must not be a leaf in ref_soln.  If p(2) /= 0, the edge must be the
! base of elem.  If p(2) == 0, the edge can be any side of elem.
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
! <u,v>_e_i = integral_e_i du/dt dv/dt d(x,y)
! where d/dt is the derivative along the tangent direction of e_i.
! If (x_1,y_1) and (x_2,y_2) are the starting and ending endpoints of e_i,
! and r = sqrt((x2-x1)**2 + (y2-y1)**2), then for any F(x,y)
! dF/dt = (x2-x1)/r dF/dx + (y2-y1)/r dF/dy
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: ref_soln, old_grid
integer, intent(in) :: edge, elem, p(2)
real(my_real), intent(out) :: diff(3)
!----------------------------------------------------
! Local variables:

integer :: endpt(3), n, nmax, qorder, nqpoints, jerr, loop, qp, i, j, &
           halfqp, child1, child2, child3, child4, degree(4), ss, comp
integer, allocatable :: ipiv(:)
real(my_real) :: du1dt, xline(2), yline(2), xtri(3), ytri(3), deltax, deltay, &
                 edge_length, integrand
real(my_real), allocatable :: a(:,:), b(:,:), dphidx(:,:), dphidy(:,:), &
                              dphidt(:,:),dudx(:,:,:), dudy(:,:,:), dudt(:,:,:)
real(my_real), pointer :: qweights(:), xquad(:), yquad(:), yquadu(:), &
                          qweightshalf(:), xquadhalf(:), yquadhalf(:)
!----------------------------------------------------
! Begin executable code

   ss = old_grid%system_size

! initialize result

   diff = 0.0_my_real

! determine the elements in the fine grid that contain the children of edge

! get the children of elem in the reference grid.  These are used for
! evaluating basis functions along the children edges.

   child1 = hash_decode_key(2*ref_soln%element(elem)%gid,ref_soln%elem_hash)
   child2 = hash_decode_key(2*ref_soln%element(elem)%gid+1,ref_soln%elem_hash)
   if (child1 == HASH_NOT_FOUND .or. child2 == HASH_NOT_FOUND) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("didn't find child element in weighted_H1_seminorm_diff_sq")
      stop
   endif

! Endpoints of the edges.  1 is first vertex of edge, 2 is midpoint of edge,
! 3 is second vertex of edge

   endpt(1) = ref_soln%edge(edge)%vertex(1)
   endpt(2) = ref_soln%element(child1)%vertex(3)
   endpt(3) = ref_soln%edge(edge)%vertex(2)

! Determine the elements to use for evaluating the solution in the reference
! grid.  If either of the children was further refined (for compatibility), then
! the element to use is its first child.

   if (ref_soln%element(child1)%isleaf) then
      child3 = child1
   else
      child3 = hash_decode_key(2*ref_soln%element(child1)%gid, &
                                ref_soln%elem_hash)
   endif
   if (ref_soln%element(child2)%isleaf) then
      child4 = child2
   else
      child4 = hash_decode_key(2*ref_soln%element(child2)%gid, &
                                ref_soln%elem_hash)
   endif
   if (child3 == HASH_NOT_FOUND .or. child4 == HASH_NOT_FOUND) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("didn't find child element in weighted_H1_seminorm_diff_sq")
      stop
   endif

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

      qorder = max(2*n,(ref_soln%edge(edge)%degree)+n)
      qorder = max(qorder,2*(ref_soln%edge(edge)%degree))

! the quadrature rule is exact for polynomials up to order 2*qorder-1, so
! get the order from the degree (currently in qorder) by (qorder+1)/2, except
! add 2 because of truncation in integer division

      qorder = min((qorder+2)/2,MAX_QUAD_ORDER_LINE)

! for multicomponent solutions, each component is computed independent of the
! other components; the squares of the norms are summed to get the whole norm

      do comp=1,ss

! u_1 is the linear vertex interpolant of u_fine.  The tangential
! derivative, i.e. the slope of the linear interpolant, is given by the
! change in solution divided by the length of the line, which is the
! subinterval if p(2)/=0 and the whole interval if p(2)==0.


         if (p(2)==0) then
            deltax = ref_soln%vertex(endpt(3))%coord%x - &
                     ref_soln%vertex(endpt(1))%coord%x
            deltay = ref_soln%vertex(endpt(3))%coord%y - &
                     ref_soln%vertex(endpt(1))%coord%y
            edge_length = sqrt(deltax**2+deltay**2)
            du1dt = (ref_soln%vertex_solution(endpt(3),comp,1) - &
                     ref_soln%vertex_solution(endpt(1),comp,1)) / edge_length
         else
            if (loop==1) then
               deltax = ref_soln%vertex(endpt(2))%coord%x - &
                        ref_soln%vertex(endpt(1))%coord%x
               deltay = ref_soln%vertex(endpt(2))%coord%y - &
                        ref_soln%vertex(endpt(1))%coord%y
               edge_length = sqrt(deltax**2+deltay**2)
               du1dt = (ref_soln%vertex_solution(endpt(2),comp,1) - &
                        ref_soln%vertex_solution(endpt(1),comp,1)) / edge_length
            else
            deltax = ref_soln%vertex(endpt(3))%coord%x - &
                        ref_soln%vertex(endpt(2))%coord%x
               deltay = ref_soln%vertex(endpt(3))%coord%y - &
                        ref_soln%vertex(endpt(2))%coord%y
               edge_length = sqrt(deltax**2+deltay**2)
               du1dt = (ref_soln%vertex_solution(endpt(3),comp,1) - &
                        ref_soln%vertex_solution(endpt(2),comp,1)) / edge_length
            endif
         endif

! matrix entries are given by
! a_ij = integral dphi_i/dt dphi_j/dt
!  where dphi_k/dt = deltax/edge_length dphi_k/dx + deltay/edge_length dphi_k/dy
! and the right hand side is given by
! b_i = integral [du_fine/dt - du1dt] dphi_i/dt

! get the quadrature rules.

! For the whole interval (p(2)==0), need to do the quadrature of each half to
! get the exact integral, so make a double length quadrature rule by
! concatinating the rules for each interval.

         if (p(2)==0) then
            xline(1) = ref_soln%vertex(endpt(1))%coord%x
            xline(2) = (ref_soln%vertex(endpt(1))%coord%x + &
                        ref_soln%vertex(endpt(3))%coord%x) / 2
            yline(1) = ref_soln%vertex(endpt(1))%coord%y
            yline(2) = (ref_soln%vertex(endpt(1))%coord%y + &
                        ref_soln%vertex(endpt(3))%coord%y) / 2
            call quadrature_rule_line(qorder,xline,yline,nqpoints,qweightshalf,&
                                      xquadhalf,yquadhalf,jerr)
            allocate(qweights(2*nqpoints),xquad(2*nqpoints),yquad(2*nqpoints))
            xquad(1:nqpoints) = xquadhalf
            yquad(1:nqpoints) = yquadhalf
            qweights(1:nqpoints) = qweightshalf
            deallocate(qweightshalf,xquadhalf,yquadhalf)
            xline(1) = (ref_soln%vertex(endpt(1))%coord%x + &
                        ref_soln%vertex(endpt(3))%coord%x) / 2
            xline(2) = ref_soln%vertex(endpt(3))%coord%x
            yline(1) = (ref_soln%vertex(endpt(1))%coord%y + &
                        ref_soln%vertex(endpt(3))%coord%y) / 2
            yline(2) = ref_soln%vertex(endpt(3))%coord%y
            call quadrature_rule_line(qorder,xline,yline,nqpoints,qweightshalf,&
                                      xquadhalf,yquadhalf,jerr)
            xquad(nqpoints+1:2*nqpoints) = xquadhalf
            yquad(nqpoints+1:2*nqpoints) = yquadhalf
            qweights(nqpoints+1:2*nqpoints) = qweightshalf
            deallocate(qweightshalf,xquadhalf,yquadhalf)
            nqpoints = 2*nqpoints

! for the half intervals, get the quadrature rule for this interval

         else
            if (loop==1) then
               xline(1) = ref_soln%vertex(endpt(1))%coord%x
               xline(2) = ref_soln%vertex(endpt(2))%coord%x
               yline(1) = ref_soln%vertex(endpt(1))%coord%y
               yline(2) = ref_soln%vertex(endpt(2))%coord%y
               call quadrature_rule_line(qorder,xline,yline, &
                                         nqpoints,qweights,xquad,yquad,jerr)
            else
               xline(1) = ref_soln%vertex(endpt(2))%coord%x
               xline(2) = ref_soln%vertex(endpt(3))%coord%x
               yline(1) = ref_soln%vertex(endpt(2))%coord%y
               yline(2) = ref_soln%vertex(endpt(3))%coord%y
               call quadrature_rule_line(qorder,xline,yline, &
                                         nqpoints,qweights,xquad,yquad,jerr)
            endif
         endif

! evaluate the derivatives of the bubble basis functions at the quadrature points

         if (p(2) == 0) then
            xtri(1) = ref_soln%vertex(ref_soln%element(elem)%vertex(1))%coord%x
            xtri(2) = ref_soln%vertex(ref_soln%element(elem)%vertex(2))%coord%x
            xtri(3) = ref_soln%vertex(ref_soln%element(elem)%vertex(3))%coord%x
            ytri(1) = ref_soln%vertex(ref_soln%element(elem)%vertex(1))%coord%y
            ytri(2) = ref_soln%vertex(ref_soln%element(elem)%vertex(2))%coord%y
            ytri(3) = ref_soln%vertex(ref_soln%element(elem)%vertex(3))%coord%y
            if (ref_soln%element(elem)%edge(1) == edge) then
               degree = (/p(loop),0,0,0/)
            elseif (ref_soln%element(elem)%edge(2) == edge) then
               degree = (/0,p(loop),0,0/)
            elseif (ref_soln%element(elem)%edge(3) == edge) then
               degree = (/0,0,p(loop),0/)
            else
               ierr = PHAML_INTERNAL_ERROR
               call fatal("weighted_H1_seminorm_diff_sq: didn't find edge in element")
            endif
         else
            if (loop==1) then
               xtri(1)= ref_soln%vertex(ref_soln%element(child1)%vertex(1))%coord%x
               xtri(2)= ref_soln%vertex(ref_soln%element(child1)%vertex(2))%coord%x
               xtri(3)= ref_soln%vertex(ref_soln%element(child1)%vertex(3))%coord%x
               ytri(1)= ref_soln%vertex(ref_soln%element(child1)%vertex(1))%coord%y
               ytri(2)= ref_soln%vertex(ref_soln%element(child1)%vertex(2))%coord%y
               ytri(3)= ref_soln%vertex(ref_soln%element(child1)%vertex(3))%coord%y
            else
               xtri(1)= ref_soln%vertex(ref_soln%element(child2)%vertex(1))%coord%x
               xtri(2)= ref_soln%vertex(ref_soln%element(child2)%vertex(2))%coord%x
               xtri(3)= ref_soln%vertex(ref_soln%element(child2)%vertex(3))%coord%x
               ytri(1)= ref_soln%vertex(ref_soln%element(child2)%vertex(1))%coord%y
               ytri(2)= ref_soln%vertex(ref_soln%element(child2)%vertex(2))%coord%y
               ytri(3)= ref_soln%vertex(ref_soln%element(child2)%vertex(3))%coord%y
            endif
            degree = (/0,p(loop),0,0/)
         endif

         if (n > 0) then
            allocate(dphidx(p(loop)+2,nqpoints),dphidy(p(loop)+2,nqpoints), &
                     dphidt(p(loop)+2,nqpoints))
            call p_hier_basis_func(xquad,yquad,xtri,ytri,degree,"a", &
                                   basisx=dphidx,basisy=dphidy)
            dphidx(1:n,:) = dphidx(4:p(loop)+2,:)
            dphidy(1:n,:) = dphidy(4:p(loop)+2,:)
            dphidt(1:n,:) = dphidx(1:n,:)*deltax/edge_length + &
                            dphidy(1:n,:)*deltay/edge_length
         endif

! evaluate the derivatives of u_fine at the quadrature points

         allocate(dudx(1,1,nqpoints),dudy(1,1,nqpoints),dudt(1,1,nqpoints))

! If there is a single interval (p(2)==0) and edge is the base, then the
! quadrature points are split among the two children of edge.  Those for which
! the quadrature points on the unit interval are greater than 1/2 are in the
! second child.  If there is a single interval that is not the base, the
! points are all in the child that contains edge.  If there are two intervals,
! the quadrature points are all in the first or second child depending on loop.

         if (p(2)==0) then
            if (ref_soln%element(elem)%edge(1) == edge) then
               call evaluate_soln_local(ref_soln,xquad,yquad,child2,(/comp/), &
                                        (/1/),ux=dudx,uy=dudy)
            elseif (ref_soln%element(elem)%edge(2) == edge) then
               call evaluate_soln_local(ref_soln,xquad,yquad,child1,(/comp/), &
                                        (/1/),ux=dudx,uy=dudy)
            else
              halfqp = nqpoints/2
              call evaluate_soln_local(ref_soln,xquad(1:halfqp), &
                                       yquad(1:halfqp),child3,(/comp/),(/1/), &
                                    ux=dudx(:,:,1:halfqp),uy=dudy(:,:,1:halfqp))
              call evaluate_soln_local(ref_soln,xquad(halfqp+1:), &
                                       yquad(halfqp+1:),child4,(/comp/),(/1/), &
                                  ux=dudx(:,:,halfqp+1:),uy=dudy(:,:,halfqp+1:))
            endif
         else
            if (loop==1) then
               call evaluate_soln_local(ref_soln,xquad,yquad,child3,(/comp/), &
                                        (/1/),ux=dudx,uy=dudy)
            else
               call evaluate_soln_local(ref_soln,xquad,yquad,child4,(/comp/), &
                                        (/1/),ux=dudx,uy=dudy)
            endif
         endif
         dudt = dudx*deltax/edge_length + dudy*deltay/edge_length

! compute the integrals

         a = 0
         b = 0
         do qp=1,nqpoints
            do i=1,n
               do j=1,n
                  a(i,j) = a(i,j) + qweights(qp)*dphidt(i,qp)*dphidt(j,qp)
               end do
               b(i,1) = b(i,1) + qweights(qp)*(dudt(1,1,qp)-du1dt)*dphidt(i,qp)
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
               call fatal("in weighted_H1_seminorm_diff_sq, LAPACK requires single or double precision")
               stop
            endif

            if (jerr /= 0) then
               ierr = PHAML_INTERNAL_ERROR
               call fatal("in weighted_H1_seminorm_diff_sq, dgesv failed")
               stop
            endif
         endif

! compute the integral of (du_fine/dt - du_1/dt - du_2/dt)**2

         do qp=1,nqpoints
            integrand = dudt(1,1,qp)-du1dt
            do i=1,n
               integrand = integrand - b(i,1)*dphidt(i,qp)
            end do
            diff(loop) = diff(loop) + qweights(qp)*integrand**2
         end do

! deallocate memory allocated in this loop

         deallocate(qweights,xquad,yquad,dudx,dudy,dudt)
         if (n > 0) deallocate(dphidx,dphidy,dphidt)

! next component

      end do

! next subinterval

   end do

! free workspace

   deallocate(a,b,ipiv)

! compute total integral

   diff(3) = diff(1) + diff(2)

end subroutine weighted_H1_seminorm_diff_sq

!          ----------------------------
subroutine compute_guaranteed_edge_rate(ref_soln,old_grid,edge,elem,best_h, &
                                        diff_old,guaranteed_edge_rate, &
                                        edge_rate_max,new_degree)
!          ----------------------------

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

integer :: p, degree(2), best_degree(2), side, otheredge, mate
real(my_real) :: diff_new(3), edge_rate, best_rate, diff_other(3)
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

! First determine the rate for p refinement of each of the three
! edges of elem.  The other two edges are considered because we get no
! refinement if edge is along an isocontour of the true solution and the
! approximate solution is sufficiently accurate (probably even more generally
! the solution is a polynomial of degree less than or equal to the degree of
! the edge).  This can probably also be justified in the same sense as using
! noncompetitive assignment of p's to the h refinement.

   do side=1,3
      otheredge = old_grid%element(elem)%edge(side)
      degree(1) = old_grid%edge(otheredge)%degree
      degree(2) = 0
      call weighted_H1_seminorm_diff_sq(ref_soln,old_grid,otheredge,elem, &
                                        degree,diff_other)
      degree(1) = old_grid%edge(otheredge)%degree + 1
      call weighted_H1_seminorm_diff_sq(ref_soln,old_grid,otheredge,elem, &
                                        degree,diff_new)
      edge_rate = diff_other(3)-diff_new(3)
      guaranteed_edge_rate = max(edge_rate, guaranteed_edge_rate)
   end do

   if (.not.(ref_soln%element(elem)%mate == BOUNDARY)) then
      mate = hash_decode_key(ref_soln%element(elem)%mate,ref_soln%elem_hash)
      do side=1,2
         otheredge = ref_soln%element(mate)%edge(side)
         degree(1) = ref_soln%edge(otheredge)%degree - 1
         degree(2) = 0
         call weighted_H1_seminorm_diff_sq(ref_soln,old_grid,otheredge,mate, &
                                           degree,diff_other)
         degree(1) = ref_soln%edge(otheredge)%degree
         call weighted_H1_seminorm_diff_sq(ref_soln,old_grid,otheredge,mate, &
                                           degree,diff_new)
         edge_rate = diff_other(3)-diff_new(3)
         guaranteed_edge_rate = max(edge_rate, guaranteed_edge_rate)
      end do
   endif

! begin with the best h refinement

   p = old_grid%edge(edge)%degree
   degree(1) = best_h
   degree(2) = p+1 - best_h
   best_degree = degree
   best_rate = -huge(0.0_my_real)

! set new_degree to the best h refinement to return that if the edge rate
! condition is never met

   if (present(new_degree)) new_degree = degree

! do until both child edges have degree p+1

   do

! compute |u_fine - w_new|^2 where |.| is the weighted H^1 seminorm and
! w_new is the projection based interpolant on the current h refinement of edge

      call weighted_H1_seminorm_diff_sq(ref_soln,old_grid,edge,elem, &
                                        degree,diff_new)

!                       |u_fine-w_{hp}|^2 - |u_fine-w_new|^2
! compute edge_rate =  --------------------------------------
!                              (p_1+p_2-1)-(p-1)

      edge_rate = (diff_old-diff_new(3)) / (degree(1)+degree(2)-p)
      guaranteed_edge_rate = max(edge_rate, guaranteed_edge_rate)

! if this is the best rate we have seen, keep it

      if (edge_rate > best_rate) then
         best_rate = edge_rate
         best_degree = degree
      endif

! if edge_rate_max is present, we are computing the new degrees

      if (present(edge_rate_max)) then
         if (edge_rate < edge_rate_max/3) then
            new_degree = degree
            new_degree = best_degree
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

! scale guaranteed edge rate by the edge length

   guaranteed_edge_rate = guaranteed_edge_rate * &
            sqrt((ref_soln%vertex(ref_soln%edge(edge)%vertex(1))%coord%x - &
                  ref_soln%vertex(ref_soln%edge(edge)%vertex(2))%coord%x)**2 + &
                 (ref_soln%vertex(ref_soln%edge(edge)%vertex(1))%coord%y - &
                  ref_soln%vertex(ref_soln%edge(edge)%vertex(2))%coord%y)**2)

end subroutine compute_guaranteed_edge_rate

!          -------------------------
subroutine determine_edges_to_refine(reftype,guaranteed_edge_rate, &
                                     edge_rate_max,refine_control)
!          -------------------------

!----------------------------------------------------
! This routine determines which edges should be refined
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

character(len=1), intent(inout) :: reftype(:)
real(my_real), intent(in) :: guaranteed_edge_rate(:)
real(my_real), intent(in) :: edge_rate_max
type(refine_options), intent(in) :: refine_control
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

! mark edges for which the guaranteed edge rate is too small as not to refine

   where (guaranteed_edge_rate < edge_rate_max/refine_control%inc_factor) &
      reftype = "n"

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

      if (olddeg >= 3) grid%dof = grid%dof - &
                                  grid%system_size*((olddeg-2)*(olddeg-1))/2
      if (newdeg >= 3) grid%dof = grid%dof + &
                                  grid%system_size*((newdeg-2)*(newdeg-1))/2
      if (grid%element(elem)%iown) then
         if (olddeg >= 3) grid%dof_own=grid%dof_own - &
                                  grid%system_size*((olddeg-2)*(olddeg-1))/2
         if (newdeg >= 3) grid%dof_own=grid%dof_own + &
                                  grid%system_size*((newdeg-2)*(newdeg-1))/2
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

integer :: p, degree(ndescendants), i, N_new, N_old, best_degree(ndescendants)
real(my_real) :: diff_hp(2), diff_new(ndescendants+1), element_rate, best_rate
real(my_real) :: fudge=1.0d-14
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
   best_degree = degree
   best_rate = -huge(0.0_my_real)

! put that in new_degree also, in case the rate condition is never met

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
      if (abs(element_rate) < fudge) element_rate=0
      guaranteed_element_rate = max(element_rate, guaranteed_element_rate)

! if this is the best rate we have seen, keep it

      if (element_rate > best_rate+fudge) then
         best_rate = element_rate
         best_degree = degree
      endif

! if rate_max is present, we are computing the new degrees

      if (present(rate_max)) then
         if (element_rate >= 0.0_my_real .and. element_rate < rate_max/3) then
            new_degree(1:ndescendants) = degree
            new_degree(1:ndescendants) = best_degree
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
! subroutine weighted_H1_seminorm_diff_sq.
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
           n, qorder, sub_nqpoints, nqpoints, qp, i, j, jerr, have_it, &
           container, num_leaf_descend, ss, comp
integer, allocatable :: ipiv(:)
integer, pointer :: leaf_descend(:)
real(my_real) :: xtri(3), ytri(3), xline(2), yline(2), deltax, deltay, du1dxi, &
                 xvert(3), yvert(3), integrandx, integrandy, xvert_elem(3), &
                 yvert_elem(3)
real(my_real), allocatable :: interp(:), a(:,:), b(:,:), dphidxi(:,:), &
                              dudx(:,:,:), dudy(:,:,:), phi(:,:), dphidx(:,:), &
                              dphidy(:,:), du1u2dx(:), du1u2dy(:)
real(my_real), pointer :: sub_qweights(:), sub_xquad(:), sub_yquad(:), &
                          xquadu(:), yquadu(:)
real(my_real), allocatable :: qweights(:), xquad(:), yquad(:)
!----------------------------------------------------
! Begin executable code

   ss = grid%system_size

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

! vertices of the element

      xvert_elem=ref_soln%vertex(ref_soln%element(elem)%vertex)%coord%x
      yvert_elem=ref_soln%vertex(ref_soln%element(elem)%vertex)%coord%y

! interp will contain the coefficients of the interpolant, in the same order
! that the basis functions are returned

      nbasis = 3
      do side=1,3
         nbasis = nbasis + degree(side)-1
      end do
      nbasis = nbasis + ((degree(4)-1)*(degree(4)-2))/2
      allocate(interp(nbasis))

! make a list of the leaf descendants in ref_soln of this descendant.
! Integrals must be done over the leaves because they involve solution
! derivatives so we cannot have quadrature points on the leaf element edges

      nullify(leaf_descend)
      call list_leaf_descendants(ref_soln,elem,leaf_descend,num_leaf_descend)

! for multicomponent solutions, each component is computed independent of the
! other components; the squares of the norms are summed to get the whole norm

      do comp=1,ss

! the linear interpolant is just the solution at the vertices

         interp(1:3) = ref_soln%vertex_solution(ref_soln%element(elem)%vertex, &
                                                comp,1)
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

            du1dxi = ref_soln%vertex_solution(endpt(2),comp,1) - &
                     ref_soln%vertex_solution(endpt(1),comp,1)

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

! odd order quadrature rules include the midpoint, which might be on the
! boundary between two elements where the derivative is not defined

            if (2*(qorder/2) /= qorder) qorder = qorder + 1

! if edge is not a leaf in ref_soln, the integral will not be exact, but
! will be real close with the maximum quadrature order
            qorder = MAX_QUAD_ORDER_LINE

            call quadrature_rule_line(qorder,xline,yline, &
                                      nqpoints,sub_qweights,xquadu,yquadu,jerr)
            deallocate(sub_qweights)
            call quadrature_rule_line(qorder,xtri(1:2),ytri(1:2),nqpoints, &
                                      sub_qweights,sub_xquad,sub_yquad,jerr)
            allocate(qweights(nqpoints),xquad(nqpoints),yquad(nqpoints))
            qweights = sub_qweights
            xquad = sub_xquad
            yquad = sub_yquad

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

            allocate(dudx(1,1,nqpoints),dudy(1,1,nqpoints))
            if (ref_soln%element(elem)%isleaf) then
               call evaluate_soln_local(ref_soln,xquadu,yquadu,elem,(/comp/), &
                                        (/1/),ux=dudx,uy=dudy)
            else
               do qp=1,nqpoints
                  call find_containing_leaf(xquadu(qp),yquadu(qp),elem,have_it,&
                                            container,ref_soln)
! TEMP when I parallelize this, note that have_it=0 and elem=0 if I don't
!      own the element that contains the point.  for now, just check that its
!      not zero, which it shouldn't be if nproc==1
                  if (have_it==0) then
                     ierr = PHAML_INTERNAL_ERROR
                     call fatal("couldn't find containing element in H1_elem_seminorm_diff_squared")
                     stop
                  endif
                  call evaluate_soln_local(ref_soln,xquadu(qp:qp), &
                                           yquadu(qp:qp),container,(/comp/), &
                                           (/1/),ux=dudx(:,:,qp:qp), &
                                           uy=dudy(:,:,qp:qp))
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
                  b(i,1) = b(i,1) + qweights(qp)*dphidxi(i,qp)* &
                           (deltax*dudx(1,1,qp)+deltay*dudy(1,1,qp)-du1dxi)
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

            deallocate(a,b,ipiv,sub_qweights,sub_xquad,sub_yquad,qweights, &
                       xquad,yquad,xquadu,yquadu,dudx,dudy)
            if (n > 0) deallocate(dphidxi)

! next edge

         end do

! element interior

! polynomial degrees involved in integrals are:
! matrix: 2*(degree(4)-1)
! rhs:    degree(4)+max(degree(u_fine),degree(1:3))
! norm:   2*max(degree(u_fine),degree(1:4))
! Pick a quadrature rule that is exact for the largest of these, which happens
! to be the norm, assuming the edge degrees are not larger than the interior

         qorder = 2*max(ref_soln%element(elem)%degree,maxval(degree))
         qorder = min(qorder,MAX_QUAD_ORDER_TRI)

! determine a set of quadrature points using quadrature rules in each of
! the leaf descendants of elem

         do i=1,num_leaf_descend

! vertices of the leaf descendant, or subelement

            xvert=ref_soln%vertex(ref_soln%element(leaf_descend(i))%vertex)%coord%x
            yvert=ref_soln%vertex(ref_soln%element(leaf_descend(i))%vertex)%coord%y

! quadrature points in this subelement

            call quadrature_rule_tri(qorder,xvert,yvert,sub_nqpoints, &
                                   sub_qweights,sub_xquad,sub_yquad,jerr,.true.)

! at the first subelement, allocate space for composite quadrature rule and
! basis functions

            if (i==1) then
               nqpoints = sub_nqpoints * num_leaf_descend
               allocate(qweights(nqpoints),xquad(nqpoints),yquad(nqpoints), &
                        phi(nbasis,nqpoints),dphidx(nbasis,nqpoints), &
                        dphidy(nbasis,nqpoints))
            endif

! copy subelement quadrature rule into composite quadrature rule

            qweights(1+(i-1)*sub_nqpoints:i*sub_nqpoints) = sub_qweights
            xquad(1+(i-1)*sub_nqpoints:i*sub_nqpoints) = sub_xquad
            yquad(1+(i-1)*sub_nqpoints:i*sub_nqpoints) = sub_yquad

! evaluate the derivatives of all basis functions at the quadrature points

            call p_hier_basis_func(sub_xquad,sub_yquad,xvert_elem,yvert_elem, &
                         degree,"a", &
                         basis=phi(:,1+(i-1)*sub_nqpoints:i*sub_nqpoints), &
                         basisx=dphidx(:,1+(i-1)*sub_nqpoints:i*sub_nqpoints), &
                         basisy=dphidy(:,1+(i-1)*sub_nqpoints:i*sub_nqpoints))

         end do

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

         allocate(dudx(1,1,nqpoints),dudy(1,1,nqpoints))
         do i=1,num_leaf_descend
            call evaluate_soln_local(ref_soln, &
                             xquad(1+(i-1)*sub_nqpoints:i*sub_nqpoints), &
                             yquad(1+(i-1)*sub_nqpoints:i*sub_nqpoints), &
                             leaf_descend(i),(/comp/),(/1/), &
                             ux=dudx(:,:,1+(i-1)*sub_nqpoints:i*sub_nqpoints), &
                             uy=dudy(:,:,1+(i-1)*sub_nqpoints:i*sub_nqpoints))

         end do


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
                    sub_qweights,sub_xquad,sub_yquad,qweights,xquad,yquad)

! next component

      end do
      deallocate(interp,leaf_descend)

! next descendant

   end do

! compute the total seminorm over all descendants

   diff(ndescendants+1) = sum(diff(1:ndescendants))

end subroutine H1_elem_seminorm_diff_squared

!                    ---------------------
recursive subroutine list_leaf_descendants(grid,elem,leaf_descend, &
                                           num_leaf_descend)
!                    ---------------------

!----------------------------------------------------
! This routine returns a list of all the leaf descendants of elem in grid.
! On input, the pointer status of leaf_descend must not be undefined.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: elem
integer, pointer :: leaf_descend(:)
integer, intent(out) :: num_leaf_descend
!----------------------------------------------------
! Local variables:

integer :: child(MAX_CHILD),i
integer, pointer :: temp(:)
!----------------------------------------------------
! Begin executable code

child = get_child_lid(grid%element(elem)%gid,ALL_CHILDREN,grid%elem_hash)

! if elem is a leaf, put it on the list

if (child(1) == NO_CHILD) then
   if (associated(leaf_descend)) then
      num_leaf_descend = size(leaf_descend) + 1
      allocate(temp(num_leaf_descend))
      temp(1:num_leaf_descend-1) = leaf_descend(1:num_leaf_descend-1)
      deallocate(leaf_descend)
      leaf_descend => temp
      leaf_descend(num_leaf_descend) = elem
   else
      allocate(leaf_descend(1))
      leaf_descend(1) = elem
      num_leaf_descend = 1
   endif

! otherwise, make the leaf descendant lists of the children

else
   do i=1,MAX_CHILD
      call list_leaf_descendants(grid,child(i),leaf_descend,num_leaf_descend)
   end do
endif

end subroutine list_leaf_descendants

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

   grid%dof = grid%dof + grid%system_size*(newdeg - olddeg)
   if (grid%element(grid%edge(edge)%assoc_elem)%iown) then
      grid%dof_own = grid%dof_own + grid%system_size*(newdeg - olddeg)
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

!===========================================================================
! REFSOLN_ELEM
!===========================================================================

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

! Create the reference solution by copying the current grid, performing
! uniform h and p refinements, and solving on the fine grid.

   call copy_grid(grid,ref_soln)
   loc_refcont = refine_control
   loc_refcont%error_estimator = HIERARCHICAL_COEFFICIENT
   loc_refcont%reftype = H_UNIFORM
   loc_iocont = io_control
   loc_iocont%print_linsys_when = NEVER
   loc_iocont%print_error_when = NEVER
   call refine_uniform_h(ref_soln,loc_refcont)
   loc_refcont%reftype = P_UNIFORM
   call refine_uniform_p(ref_soln,loc_refcont)
   call solve(ref_soln,procs,loc_iocont,solver_control,still_sequential, &
              .false.,no_time=.true.)

! MASTER doesn't participate further; it was only here to monitor solve

if (my_proc(procs) == MASTER) then
   call deallocate_grid(ref_soln)
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
call refine_uniform_h(ref_soln,loc_refcont)
loc_refcont%reftype = P_UNIFORM
call refine_uniform_p(ref_soln,loc_refcont)
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
! The projection onto a linear element is singular if h is small enough that
! the u term in the A inner products is below machine epsilon relative to the
! gradient term.  In that case, use the value of the linear part of the
! reference solution as the projection.
   if (p_u==1) then
      b(1,1)=ref_soln%vertex_solution(grid%element(elem)%vertex(1),1,1)
      b(2,1)=ref_soln%vertex_solution(grid%element(elem)%vertex(2),1,1)
      b(3,1)=ref_soln%vertex_solution(grid%element(elem)%vertex(3),1,1)
      jerr = 0
   else
      ierr = PHAML_INTERNAL_ERROR
      call fatal("LAPACK returned error during H^1 projection for refsoln_elem",intlist=(/0,jerr/))
      stop
   endif
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
      subset(3*p0+k+2*i+1:Q) = (/(l,l=3*p0+7,Q-k-2*i+6)/)
      subset(Q+1) = P+1
      subset(Q+2:Q+p0+j) = (/(l,l=P+2,P+p0+j)/)
      subset(Q+p0+j+1:Q+2*p0+2*j-1) = (/(l,l=P+p0+3,P+2*p0+j+1)/)
      subset(Q+2*p0+2*j:R) = (/(l,l=P+2*p0+4,P+R-Q-2*j+4)/)

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

! discard candidates for which the projection is singular

      if (jerr /= 0) then
         projerr(candidate) = -1.0_my_real
      else
         projerr(candidate) = 0.0_my_real
      endif

! compute the projection error on (p0+i,p0+j)

      if (projerr(candidate) == 0.0_my_real) then
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
               projerr(candidate) = projerr(candidate) + qweights_child1(k)*( &
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
                     proj  = proj + b(l,i2)* basis_child2(invrenum(subset(l)),k)
                     projx = projx+ b(l,i2)*basisx_child2(invrenum(subset(l)),k)
                     projy = projy+ b(l,i2)*basisy_child2(invrenum(subset(l)),k)
                  endif
               end do
               projerr(candidate) = projerr(candidate) + qweights_child2(k)*( &
                                            (ux_ref_child2(i2,1,k)-projx)**2 + &
                                            (uy_ref_child2(i2,1,k)-projy)**2 + &
                                            ( u_ref_child2(i2,1,k)-proj )**2)
            end do
         end do
         projerr(candidate) = sqrt(projerr(candidate))
      endif

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
! TEMP091215 raise N to 1/3; or not; or don't do log of errors
         reduction = (logprojerr_u-logprojerr(candidate))/(nbasis(candidate)-nbasis_u)
!         reduction = (projerr_u-projerr(candidate))/(nbasis(candidate)-nbasis_u)
!         reduction = (logprojerr_u-logprojerr(candidate))/(nbasis(candidate)**.33333d0-nbasis_u**.33333d0)
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

!===========================================================================
! NLP
!===========================================================================

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

! initializations for the NLP hp-adaptive strategy

call nlp_initialize(grid,refine_control)

! make sure tau is such that we can get a feasible solution by making sure
! it is at least 1.5 times the error estimate of a maximally refined grid.
! Note that nlp_errest scales its result by nlp_tau, so undo that.

allocate(new_params(2*nlp_nelem))
do i=1,nlp_nelem
   new_params(i) = grid%element(nlp_elem_list(i))%level + &
                   refine_control%nlp_max_h_inc
   new_params(i+nlp_nelem) = grid%element(nlp_elem_list(i))%degree + &
                             refine_control%nlp_max_p_inc
end do

nlp_tau = max(nlp_tau,nlp_tau*1.5_my_real*nlp_errest(new_params))

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

! approximate the smoothness the same as in the PRIOR2P strategy, but don't
! allow m<1.5 because it can cause stalling because of a large tau, or bigger
! than 10.

allocate(nlp_m(nlp_nelem))
do i=1,nlp_nelem
   nlp_m(i) = min(10.0_my_real, &
                  max(1.5_my_real,prior2p_h1_regularity(grid,nlp_elem_list(i))))
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

end module hp_strategies

!===========================================================================
! External routines that interface to algencan for NLP
!===========================================================================

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

      use hp_strategies
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

      use hp_strategies
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
      use hp_strategies
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
      use hp_strategies
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
