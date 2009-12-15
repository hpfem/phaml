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

module basis_functions

!----------------------------------------------------
! This module contains routines to compute the basis functions
!----------------------------------------------------


!----------------------------------------------------
! Other modules used are:

use global
use message_passing, only: fatal
!----------------------------------------------------

implicit none
private
public p_hier_basis_func, h_hier_basis_func, nodal_basis_func

!----------------------------------------------------
! The following generic interfaces are defined:
   
interface p_hier_basis_func
   module procedure p_hier_basis_func_a, p_hier_basis_func_s
end interface

interface nodal_basis_func
   module procedure nodal_basis_func_a, nodal_basis_func_s
end interface

interface h_hier_basis_func
   module procedure h_hier_basis_func_a, h_hier_basis_func_s
end interface

!----------------------------------------------------
! Defined constants:

! choice between the Szabo-Babuska and Carnevali basis
integer, parameter :: SZABO = 1, CARNEVALI = 2
integer, parameter :: use_basis = SZABO

!----------------------------------------------------

contains

!-------------------------------------------------------
! ROUTINES FOR P-HIERARCHICAL BASIS
!-------------------------------------------------------

!          -------------------
subroutine p_hier_basis_func_a(x,y,xvert,yvert,degree,set,basis,basisx,basisy, &
                               basisxx,basisyy,basisxy)
!          -------------------

!----------------------------------------------------
! This routine computes p-hierarchical basis functions and derivatives at (x,y)
! for the triangle given by (xvert,yvert).  basis, basisx, basisy, basisxx,
! basisyy and basisxy are all optional; only those present are computed.
! degree determines the polynomial degree of the basis functions; degree(1:3)
! is the degree of the edge bases opposite the corresponding vertex, degree(4)
! is the degree of the face bases.
!
! set determines which subset of the bases is returned: "a" requests all
! basis functions, "b" requests only the black bases, and "r" the red bases.
!
!  The order in which the bases are returned:
!  - the vertex bases at vertex 1, then 2, then 3
!  - the edge bases at the side opposite vertex 1, then 2, then 3
!    - for each side, the quadratic basis, then the cubic, ...
!  - the face bases, in order of degree
!
! The basis functions are as defined in
! Adjerid, S., Aiffa, M. and Flaherty, J.E., Hierarchical Finite Element
! Bases for Triangular and Tetrahedral Elements, Comp. Meth. Appl. Mech. Eng.,
! 190, (1999), pp 2925-2941.
! except that for the Carnevali edge bases in equation 2.3d, i is replaced
! by i-1 everywhere, which gives the right polynomial degree and I think
! agrees with the expression in the thesis of Whiting (another RPI student).
!
! There are two choices for the basis functions, selected by a defined
! constant above.  Either the Szabo-Babuska bases defined in
! Szabo, B. and Babuska, I., Introduction to Finite Element Analysis,
! John Wiley and Sons, New York, 1989.
! [Note: I was unable to find this book as cited in Adjerid et al.  Instead
!  I found Finite Element Analysis, John Wiley and Sons, 1991]
! or the Carnevali bases defined in
! Carnevali, P., Morris, R. B., Tsuji, Y. and Taylor, G., New Basis Functions
! and Computational Procedures for the p-Version Finite Element Analysis,
! Int. J. Num. Meth. Engng., 36, (1993), pp. 3759-3779.
!
! Choosing the Szabo-Babuska basis and performing block eliminations of the
! rows/columns corresponding to related face bases in the stiffness matrix
! should be equivalent to using the Adjerid basis for Poisson equations.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x(:), y(:), xvert(3), yvert(3)
integer, intent(in) :: degree(4)
character(len=1), intent(in) :: set
real(my_real), intent(out), optional :: basis(:,:), basisx(:,:), basisy(:,:), &
                                        basisxx(:,:), basisyy(:,:), basisxy(:,:)

!----------------------------------------------------
! Local variables:

real(my_real) :: zeta(size(x),3), dzdx(3), dzdy(3), zetav1v2(size(x)), &
                 dzdxzeta12(size(x)), dzdyzeta12(size(x)), const, &
                 zetaprod(size(x)), forbasisx(size(x)), forbasisy(size(x)), &
                 forbasisxx(size(x)), forbasisyy(size(x)), forbasisxy(size(x))
real(my_real), allocatable :: legendre_d0(:,:), legendre_d1(:,:), &
                              legendre_d2(:,:), legendre_d3(:,:), &
                              legendre2_d0(:,:), legendre2_d1(:,:), &
                              legendre2_d2(:,:), c(:), lbasis(:,:)
real(my_real) :: Fr1r2(size(x)),dFdL1(size(x)),dFdL2(size(x)),dFdL3(size(x)), &
                 d2FdL(size(x),3,3)
integer :: i, j, nbasis, nedge, npoint, side, v1, v2, isub, astat1, astat2, &
           astat3, r, r1, r2, ik, k, p, nface, nedger, nedgeb, nfacer, nfaceb, &
           nvert

!----------------------------------------------------
! External routines:

interface
   double precision function dbinom(n,m)
   integer :: n,m
   end function dbinom
end interface

!----------------------------------------------------
! Begin executable code

if (present(basis)) basis=0.0_my_real
if (present(basisx)) basisx=0.0_my_real
if (present(basisy)) basisy=0.0_my_real
if (present(basisxx)) basisxx=0.0_my_real
if (present(basisyy)) basisyy=0.0_my_real
if (present(basisxy)) basisxy=0.0_my_real

! check that the size of the return arrays is big enough

npoint = size(x)
nedgeb = 0
nedger = 0
do i=1,3
   if (degree(i) > 1) then
      nedgeb = nedgeb + degree(i)-2
      nedger = nedger + 1
   endif
end do
if (degree(4) > 2) then
   nfaceb = ((degree(4)-2)*(degree(4)-3))/2
   nfacer = degree(4)-2
else
   nfaceb = 0
   nfacer = 0
endif

! always consider the vertex bases to be black, even with degree=1
select case(set)
case ("a")
   nvert = 3
   nedge = nedgeb + nedger
   nface = nfaceb + nfacer
case ("b")
   nvert = 3
   nedge = nedgeb
   nface = nfaceb
case ("r")
   nvert = 0
   nedge = nedger
   nface = nfacer
end select

nbasis = nvert + nedge + nface

if (present(basis)) then
   if (size(basis,1) < nbasis .or. size(basis,2) < npoint) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("Size of basis array is too small in p_hier_basis_func.", &
                 intlist = (/size(basis,1),size(basis,2),nbasis,npoint/))
      stop
   endif
endif

if (present(basisx)) then
   if (size(basisx,1) < nbasis .or. size(basisx,2) < npoint) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("Size of basisx array is too small in p_hier_basis_func.", &
                 intlist = (/size(basisx,1),size(basisx,2),nbasis,npoint/))
      stop
   endif
endif

if (present(basisy)) then
   if (size(basisy,1) < nbasis .or. size(basisy,2) < npoint) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("Size of basisy array is too small in p_hier_basis_func.", &
                 intlist = (/size(basisy,1),size(basisy,2),nbasis,npoint/))
      stop
   endif
endif

if (present(basisxx)) then
   if (size(basisxx,1) < nbasis .or. size(basisxx,2) < npoint) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("Size of basisxx array is too small in p_hier_basis_func.", &
                 intlist = (/size(basisxx,1),size(basisxx,2),nbasis,npoint/))
      stop
   endif
endif

if (present(basisyy)) then
   if (size(basisyy,1) < nbasis .or. size(basisyy,2) < npoint) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("Size of basisyy array is too small in p_hier_basis_func.", &
                 intlist = (/size(basisyy,1),size(basisyy,2),nbasis,npoint/))
      stop
   endif
endif

if (present(basisxy)) then
   if (size(basisxy,1) < nbasis .or. size(basisxy,2) < npoint) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("Size of basisxy array is too small in p_hier_basis_func.", &
                 intlist = (/size(basisxy,1),size(basisxy,2),nbasis,npoint/))
      stop
   endif
endif

if (present(basis)) then
   allocate(lbasis(size(basis,1),size(basis,2)),stat=astat1)
elseif (present(basisx)) then
   allocate(lbasis(size(basisx,1),size(basisx,2)),stat=astat1)
elseif (present(basisy)) then
   allocate(lbasis(size(basisy,1),size(basisy,2)),stat=astat1)
elseif (present(basisxx)) then
   allocate(lbasis(size(basisxx,1),size(basisxx,2)),stat=astat1)
elseif (present(basisyy)) then
   allocate(lbasis(size(basisyy,1),size(basisyy,2)),stat=astat1)
elseif (present(basisxy)) then
   allocate(lbasis(size(basisxy,1),size(basisxy,2)),stat=astat1)
else
   ierr = PHAML_INTERNAL_ERROR
   call fatal("no return argument requested in p_hier_basis_func_a")
   stop
endif
if (astat1 /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in p_hier_basis_func_a")
   stop
endif

! compute the barycentric coordinates of each point

call barycentric(x,y,xvert,yvert,zeta,dzdx,dzdy)

!-------------
! Vertex bases

! These are the usual linear bases

if (set == "a" .or. set == "b") then

   lbasis(1:3,:) = transpose(zeta)

   if (present(basisx)) then
      basisx(1,:) = dzdx(1)
      basisx(2,:) = dzdx(2)
      basisx(3,:) = dzdx(3)
   endif

   if (present(basisy)) then
      basisy(1,:) = dzdy(1)
      basisy(2,:) = dzdy(2)
      basisy(3,:) = dzdy(3)
   endif

   if (present(basisxx)) basisxx = 0.0_my_real
   if (present(basisyy)) basisyy = 0.0_my_real
   if (present(basisxy)) basisxy = 0.0_my_real

endif

!-----------
! Edge bases

if (maxval(degree) > 1) then

 select case (use_basis)

! the Szabo-Babuska basis

 case (SZABO)

   allocate(legendre_d1(size(x),maxval(degree)-1),c(maxval(degree)-1), &
            stat=astat1)
   if (present(basisx) .or. present(basisy) .or. present(basisxx) .or. &
       present(basisyy) .or. present(basisxy)) then
      allocate(legendre_d2(size(x),maxval(degree)-1),stat=astat2)
   else
      astat2 = 0
   endif
   if (present(basisxx) .or. present(basisyy) .or. present(basisxy)) then
      allocate(legendre_d3(size(x),maxval(degree)-1),stat=astat3)
   else
      astat3 = 0
   endif
   if (astat1 /= 0 .or. astat2 /= 0 .or. astat3 /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in p_hier_basis_func_a")
      stop
   endif

   do i=1,maxval(degree)-1
      c(i) = -8*sqrt(4*i+2.0_my_real)/(i*(i+1))
   end do

! for each edge ...

   isub = nvert
   do side=1,3

      if (degree(side) > 1) then

! identify the endpoints of the edge

         if (side == 2) then
            v2 = mod(side  ,3)+1
            v1 = mod(side+1,3)+1
         else
            v1 = mod(side  ,3)+1
            v2 = mod(side+1,3)+1
         endif

! evaluate the derivatives of the Legendre polynomials

         if (present(basisxx) .or. present(basisyy) .or. present(basisxy)) then
            call dlegendre(zeta(:,v2)-zeta(:,v1),degree(side)-1, &
                           deriv1=legendre_d1,deriv2=legendre_d2, &
                           deriv3=legendre_d3)
         elseif (present(basisx) .or. present(basisy)) then
            call dlegendre(zeta(:,v2)-zeta(:,v1),degree(side)-1, &
                           deriv1=legendre_d1,deriv2=legendre_d2)
         else
            call dlegendre(zeta(:,v2)-zeta(:,v1),degree(side)-1, &
                           deriv1=legendre_d1)
         endif

         zetav1v2 = zeta(:,v1)*zeta(:,v2)
         dzdxzeta12 = dzdx(v1)*zeta(:,v2)+dzdx(v2)*zeta(:,v1)
         dzdyzeta12 = dzdy(v1)*zeta(:,v2)+dzdy(v2)*zeta(:,v1)

! compute the basis functions and derivatives as requested

         if (set == "a" .or. set == "b") then
            if (present(basis)) then
               do i=1,degree(side)-2
                  lbasis(isub+i,:) = c(i)*zetav1v2*legendre_d1(:,i)
               end do
            endif

            if (present(basisx)) then
               do i=1,degree(side)-2
                  basisx(isub+i,:) = c(i) * ( dzdxzeta12*legendre_d1(:,i) + &
                  zetav1v2*legendre_d2(:,i)*(dzdx(v2)-dzdx(v1)) )
               end do
            endif

            if (present(basisy)) then
               do i=1,degree(side)-2
                  basisy(isub+i,:) = c(i) * ( dzdyzeta12*legendre_d1(:,i) + &
                  zetav1v2*legendre_d2(:,i)*(dzdy(v2)-dzdy(v1)) )
               end do
            endif

            if (present(basisxx)) then
               do i=1,degree(side)-2
                  basisxx(isub+i,:) = c(i) * ( &
                          2*dzdx(v1)*dzdx(v2)*legendre_d1(:,i) + &
                          2*dzdxzeta12*legendre_d2(:,i)*(dzdx(v2)-dzdx(v1)) + &
                          zetav1v2*legendre_d3(:,i)*(dzdx(v2)-dzdx(v1))**2 )
               end do
            endif

            if (present(basisyy)) then
               do i=1,degree(side)-2
                  basisyy(isub+i,:) = c(i) * ( &
                          2*dzdy(v1)*dzdy(v2)*legendre_d1(:,i) + &
                          2*dzdyzeta12*legendre_d2(:,i)*(dzdy(v2)-dzdy(v1)) + &
                          zetav1v2*legendre_d3(:,i)*(dzdy(v2)-dzdy(v1))**2 )
               end do
            endif

            if (present(basisxy)) then
               do i=1,degree(side)-2
                  basisxy(isub+i,:) = c(i) * ( &
                          (dzdx(v1)*dzdy(v2)+dzdy(v1)*dzdx(v2))*legendre_d1(:,i) + &
                          (dzdxzeta12*(dzdy(v2)-dzdy(v1)) + dzdyzeta12*(dzdx(v2)-dzdx(v1)))*legendre_d2(:,i) + &
                          zetav1v2*legendre_d3(:,i)*(dzdx(v2)-dzdx(v1))*(dzdy(v2)-dzdy(v1)))
               end do
            endif

            isub = isub + degree(side)-2
         endif

         if (set == "a" .or. set == "r") then
            i=degree(side)-1
            if (present(basis)) then
               lbasis(isub+1,:) = zetav1v2*legendre_d1(:,i)*c(i)
            endif
            if (present(basisx)) then
               basisx(isub+1,:) = c(i) * ( dzdxzeta12*legendre_d1(:,i) + &
                  zetav1v2*legendre_d2(:,i)*(dzdx(v2)-dzdx(v1)) )
            endif
            if (present(basisy)) then
               basisy(isub+1,:) = c(i) * ( dzdyzeta12*legendre_d1(:,i) + &
                  zetav1v2*legendre_d2(:,i)*(dzdy(v2)-dzdy(v1)) )
            endif
            if (present(basisxx)) then
               basisxx(isub+1,:) = c(i) * ( &
                       2*dzdx(v1)*dzdx(v2)*legendre_d1(:,i) + &
                       2*dzdxzeta12*legendre_d2(:,i)*(dzdx(v2)-dzdx(v1)) + &
                       zetav1v2*legendre_d3(:,i)*(dzdx(v2)-dzdx(v1))**2 )
            endif
            if (present(basisyy)) then
               basisyy(isub+1,:) = c(i) * ( &
                       2*dzdy(v1)*dzdy(v2)*legendre_d1(:,i) + &
                       2*dzdyzeta12*legendre_d2(:,i)*(dzdy(v2)-dzdy(v1)) + &
                       zetav1v2*legendre_d3(:,i)*(dzdy(v2)-dzdy(v1))**2 )
            endif
            if (present(basisxy)) then
               basisxy(isub+1,:) = c(i) * ( &
                       (dzdx(v1)*dzdy(v2)+dzdy(v1)*dzdx(v2))*legendre_d1(:,i) + &
                       (dzdxzeta12*(dzdy(v2)-dzdy(v1)) + dzdyzeta12*(dzdx(v2)-dzdx(v1)))*legendre_d2(:,i) + &
                       zetav1v2*legendre_d3(:,i)*(dzdx(v2)-dzdx(v1))*(dzdy(v2)-dzdy(v1)))
            endif

            isub = isub + 1
         endif
      endif
   end do

   if (present(basisxx) .or. present(basisyy) .or. present(basisxy)) then
      deallocate(legendre_d3)
   endif
   if (present(basisx) .or. present(basisy) .or. present(basisxx) .or. &
       present(basisyy) .or. present(basisxy)) then
      deallocate(legendre_d2)
   endif
   deallocate(legendre_d1,c)

! the Carnevali basis

 case (CARNEVALI)

   if (present(basisxx) .or. present(basisyy) .or. present(basisxy)) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("Second derivatives for Carnevali basis not yet written.")
      stop
   endif

! for each edge ...

   isub = nvert
   do side=1,3

      if (degree(side) > 1) then

! identify the endpoints of the edge

         if (side == 2) then
            v2 = mod(side  ,3)+1
            v1 = mod(side+1,3)+1
         else
            v1 = mod(side  ,3)+1
            v2 = mod(side+1,3)+1
         endif

! compute the basis functions and derivatives as requested

         if (set == "a" .or. set == "b") then
            do i=1,degree(side)-2
               lbasis(isub+i,:) = 0.0_my_real
               if (present(basisx)) then
                  basisx(isub+i,:) = 0.0_my_real
               endif
               if (present(basisy)) then
                  basisy(isub+i,:) = 0.0_my_real
               endif
               do k=0,i-1
                  const = (-1)**k * dbinom(i-1,k) * dbinom(i,k) / (k+1)

                     lbasis(isub+i,:) = lbasis(isub+i,:) + &
                                    const * zeta(:,v1)**k * zeta(:,v2)**(i-1-k)

                  if (present(basisx)) then
                     if (k /= 0) then
                        basisx(isub+i,:) = basisx(isub+i,:) + const * &
                             k*dzdx(v1)*zeta(:,v1)**(k-1)*zeta(:,v2)**(i-1-k)
                     endif
                     if (k /= i-1) then
                        basisx(isub+i,:) = basisx(isub+i,:) + const * &
                             (i-1-k)*dzdx(v2)*zeta(:,v2)**(i-k-2)*zeta(:,v1)**k
                     endif
                  endif

                  if (present(basisy)) then
                     if (k /= 0) then
                        basisy(isub+i,:) = basisy(isub+i,:) + const * &
                             k*dzdy(v1)*zeta(:,v1)**(k-1)*zeta(:,v2)**(i-1-k)
                     endif
                     if (k /= i-1) then
                        basisy(isub+i,:) = basisy(isub+i,:) + const * &
                             (i-1-k)*dzdy(v2)*zeta(:,v2)**(i-k-2)*zeta(:,v1)**k
                     endif
                  endif
               end do
               if (present(basisx)) then
                  basisx(isub+i,:) = dzdx(v1)*zeta(:,v2)*lbasis(isub+i,:) + &
                                     dzdx(v2)*zeta(:,v1)*lbasis(isub+i,:) + &
                                     zeta(:,v1)*zeta(:,v2)*basisx(isub+i,:)
               endif
               if (present(basisy)) then 
                  basisy(isub+i,:) = dzdy(v1)*zeta(:,v2)*lbasis(isub+i,:) + &
                                     dzdy(v2)*zeta(:,v1)*lbasis(isub+i,:) + &
                                     zeta(:,v1)*zeta(:,v2)*basisy(isub+i,:)
               endif
               lbasis(isub+i,:) = lbasis(isub+i,:)*zeta(:,v1)*zeta(:,v2)
            end do
            isub = isub + degree(side)-2
         endif

         if (set == "a" .or. set == "r") then
            i=degree(side)-1
            lbasis(isub+1,:) = 0.0_my_real
            if (present(basisx)) then
               basisx(isub+1,:) = 0.0_my_real
            endif
            if (present(basisy)) then
               basisy(isub+1,:) = 0.0_my_real
            endif
            do k=0,i-1
               const = (-1)**k * dbinom(i-1,k) * dbinom(i,k) / (k+1)
               lbasis(isub+1,:) = lbasis(isub+1,:) + &
                                    const * zeta(:,v1)**k * zeta(:,v2)**(i-1-k)
               if (present(basisx)) then
                  if (k /= 0) then
                     basisx(isub+1,:) = basisx(isub+1,:) + const * &
                          k*dzdx(v1)*zeta(:,v1)**(k-1)*zeta(:,v2)**(i-1-k)
                  endif
                  if (k /= i-1) then
                     basisx(isub+1,:) = basisx(isub+1,:) + const * &
                          (i-1-k)*dzdx(v2)*zeta(:,v2)**(i-k-2)*zeta(:,v1)**k
                  endif
               endif
               if (present(basisy)) then
                  if (k /= 0) then
                     basisy(isub+1,:) = basisy(isub+1,:) + const * &
                          k*dzdy(v1)*zeta(:,v1)**(k-1)*zeta(:,v2)**(i-1-k)
                  endif
                  if (k /= i-1) then
                     basisy(isub+1,:) = basisy(isub+1,:) + const * &
                          (i-1-k)*dzdy(v2)*zeta(:,v2)**(i-k-2)*zeta(:,v1)**k
                  endif
               endif
            end do
            if (present(basisx)) then
               basisx(isub+1,:) = dzdx(v1)*zeta(:,v2)*lbasis(isub+1,:) + &
                                  dzdx(v2)*zeta(:,v1)*lbasis(isub+1,:) + &
                                  zeta(:,v1)*zeta(:,v2)*basisx(isub+1,:)
            endif
            if (present(basisy)) then
               basisy(isub+1,:) = dzdy(v1)*zeta(:,v2)*lbasis(isub+1,:) + &
                                  dzdy(v2)*zeta(:,v1)*lbasis(isub+1,:) + &
                                  zeta(:,v1)*zeta(:,v2)*basisy(isub+1,:)
            endif
            lbasis(isub+1,:) = lbasis(isub+1,:)*zeta(:,v1)*zeta(:,v2)
            isub = isub + 1
         endif
      endif
   end do
 end select
endif

!-----------
! Face bases

! no face functions for linear and quadratic bases

if (degree(4) > 2) then

! evaluate the derivatives of the Legendre polynomials

if (use_basis == SZABO) then
   allocate(legendre_d0(size(x),degree(4)-3),legendre2_d0(size(x),degree(4)-3),&
            stat=astat1)
   if (present(basisx) .or. present(basisy) .or. present(basisxx) .or. &
       present(basisyy) .or. present(basisxy)) then
      allocate(legendre_d1(size(x),degree(4)-3), &
               legendre2_d1(size(x),degree(4)-3),stat=astat2)
   else
      astat2 = 0
   endif
   if (present(basisxx) .or. present(basisyy) .or. present(basisxy)) then
      allocate(legendre_d2(size(x),degree(4)-3), &
               legendre2_d2(size(x),degree(4)-3),stat=astat3)
   else
      astat3 = 0
   endif
   if (astat1 /= 0 .or. astat2 /= 0 .or. astat3 /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in p_hier_basis_func_a")
      stop
   endif
   if (present(basisxx) .or. present(basisyy) .or. present(basisxy)) then
      call dlegendre(zeta(:,2)-zeta(:,1),degree(4)-3,deriv0=legendre_d0, &
                     deriv1=legendre_d1,deriv2=legendre_d2)
      call dlegendre(2*zeta(:,3)-1,degree(4)-3,deriv0=legendre2_d0, &
                     deriv1=legendre2_d1,deriv2=legendre2_d2)
   elseif (present(basisx) .or. present(basisy)) then
      call dlegendre(zeta(:,2)-zeta(:,1),degree(4)-3,deriv0=legendre_d0, &
                     deriv1=legendre_d1)
      call dlegendre(2*zeta(:,3)-1,degree(4)-3,deriv0=legendre2_d0, &
                     deriv1=legendre2_d1)
   else
      call dlegendre(zeta(:,2)-zeta(:,1),degree(4)-3,deriv0=legendre_d0)
      call dlegendre(2*zeta(:,3)-1,degree(4)-3,deriv0=legendre2_d0)
   endif
endif

! common products

zetaprod = zeta(:,1)*zeta(:,2)*zeta(:,3)
forbasisx = dzdx(1)*zeta(:,2)*zeta(:,3) + dzdx(2)*zeta(:,1)*zeta(:,3) + &
            dzdx(3)*zeta(:,1)*zeta(:,2)
forbasisy = dzdy(1)*zeta(:,2)*zeta(:,3) + dzdy(2)*zeta(:,1)*zeta(:,3) + &
            dzdy(3)*zeta(:,1)*zeta(:,2)
forbasisxx = 2*(dzdx(1)*dzdx(2)*zeta(:,3) + &
                dzdx(1)*zeta(:,2)*dzdx(3) + &
                zeta(:,1)*dzdx(2)*dzdx(3))
forbasisyy = 2*(dzdy(1)*dzdy(2)*zeta(:,3) + &
                dzdy(1)*zeta(:,2)*dzdy(3) + &
                zeta(:,1)*dzdy(2)*dzdy(3))
forbasisxy = dzdx(1)*dzdy(2)*zeta(:,3) + dzdx(1)*zeta(:,2)*dzdy(3) + &
             dzdy(1)*dzdx(2)*zeta(:,3) + zeta(:,1)*dzdx(2)*dzdy(3) + &
             dzdy(1)*zeta(:,2)*dzdx(3) + zeta(:,1)*dzdy(2)*dzdx(3)

! evaluate the face bases and derivatives

do k=3,degree(4)
   ik = ((k-3)*(k-2))/2
   if (set == "r") ik = 0
   r = nedge + nvert
   do r1 = 0,k-3
      r2 = k-3-r1
      r = r+1

      select case (use_basis)

      case (SZABO)
         if (r1==0 .and. r2==0) then
            Fr1r2 = 1
         elseif (r1==0) then
            Fr1r2 = legendre2_d0(:,r2)
         elseif (r2==0) then
            Fr1r2 = legendre_d0(:,r1)
         else
            Fr1r2 = legendre_d0(:,r1)*legendre2_d0(:,r2)
         endif
         if (present(basisx) .or. present(basisy) .or. present(basisxx) .or. &
             present(basisyy) .or. present(basisxy)) then
            if (r1==0 .and. r2==0) then
               dFdL1 = 0.0_my_real
               dFdL2 = 0.0_my_real
               dFdL3 = 0.0_my_real
            elseif (r1==0) then
               dFdL1 = 0.0_my_real
               dFdL2 = 0.0_my_real
               dFdL3 = 2*legendre2_d1(:,r2)
            elseif (r2==0) then
               dFdL1 = -legendre_d1(:,r1)
               dFdL2 =  legendre_d1(:,r1)
               dFdL3 = 0.0_my_real
            else
               dFdL1 = -legendre_d1(:,r1)*legendre2_d0(:,r2)
               dFdL2 = -dFdL1
               dFdL3 = 2*legendre_d0(:,r1)*legendre2_d1(:,r2)
            endif
         endif
         if (present(basisxx) .or. present(basisyy) .or. present(basisxy)) then
            if (r1==0 .and. r2==0) then
               d2FdL = 0.0_my_real
            elseif (r1==0) then
               d2FdL = 0.0_my_real
               d2FdL(:,3,3) = 4*legendre2_d2(:,r2)
            elseif (r2==0) then
               d2FdL = 0.0_my_real
               d2FdL(:,1,1) = legendre_d2(:,r1)
               d2FdL(:,2,2) = d2FdL(:,1,1)
               d2FdL(:,1,2) = -d2FdL(:,1,1)
               d2FdL(:,2,1) = d2FdL(:,1,2)
            else
               d2FdL(:,1,1) = legendre_d2(:,r1)*legendre2_d0(:,r2)
               d2FdL(:,1,2) = -d2FdL(:,1,1)
               d2FdL(:,2,1) = d2FdL(:,1,2)
               d2FdL(:,2,2) = d2FdL(:,1,1)
               d2FdL(:,1,3) = -2*legendre_d1(:,r1)*legendre2_d1(:,r2)
               d2FdL(:,3,1) = d2FdL(:,1,3)
               d2FdL(:,2,3) = -d2FdL(:,1,3)
               d2FdL(:,3,2) = d2FdL(:,2,3)
               d2FdL(:,3,3) = 4*legendre_d0(:,r1)*legendre2_d2(:,r2)
            endif
         endif

      case (CARNEVALI)
         Fr1r2 = 0.0_my_real
         dFdL1 = 0.0_my_real
         dFdL2 = 0.0_my_real
         dFdL3 = 0.0_my_real
         do i=0,r2
            do j=0,r1
               const = 1.0_my_real
               do p=1,i+j
                  const = const * (p*(r1+r2+2) - (p*(p-1))/2)
               end do
               const = (-0.5_my_real)**(i+j) * dbinom(r1,j) * dbinom(r1+1,j) * &
                       dbinom(r2,i) * dbinom(r2+1,i) * factorial(i) * &
                       factorial(j) * factorial(i+j) / const

               Fr1r2 = Fr1r2 + const * zeta(:,1)**(r1-j) * zeta(:,2)**(r2-i)
               if (present(basisx) .or. present(basisy)) then
                  if (j /= r1) then
                     dFdL1 = dFdL1 + const * &
                                  (r1-j)*zeta(:,1)**(r1-j-1) * zeta(:,2)**(r2-i)
                  endif
                  if (i /= r2) then
                     dFdL2 = dFdL2 + const * &
                                  zeta(:,1)**(r1-j) * (r2-i)*zeta(:,2)**(r2-i-1)
                  endif
               endif
            end do
         end do
      end select

      if (set == "a" .or. (set == "b" .and. k < degree(4)) .or. &
          (set == "r" .and. k == degree(4))) then
         if (present(basis)) then
            lbasis(ik+r,:) = zetaprod*Fr1r2
         endif
         if (present(basisx)) then
            basisx(ik+r,:) = Fr1r2*forbasisx + &
                            zetaprod*(dFdL1*dzdx(1)+dFdL2*dzdx(2)+dFdL3*dzdx(3))
         endif
         if (present(basisy)) then
            basisy(ik+r,:) = Fr1r2*forbasisy + &
                            zetaprod*(dFdL1*dzdy(1)+dFdL2*dzdy(2)+dFdL3*dzdy(3))
         endif
         if (present(basisxx)) then
            basisxx(ik+r,:) = forbasisxx*Fr1r2 + &
                        2*forbasisx*(dFdL1*dzdx(1)+dFdL2*dzdx(2)+dFdL3*dzdx(3))
            do i=1,size(x)
               basisxx(ik+r,i) = basisxx(ik+r,i) + &
                         zetaprod(i)*dot_product(dzdx,matmul(d2FdL(i,:,:),dzdx))
            end do
         endif
         if (present(basisyy)) then
            basisyy(ik+r,:) = forbasisyy*Fr1r2 + &
                        2*forbasisy*(dFdL1*dzdy(1)+dFdL2*dzdy(2)+dFdL3*dzdy(3))
            do i=1,size(x)
               basisyy(ik+r,i) = basisyy(ik+r,i) + &
                         zetaprod(i)*dot_product(dzdy,matmul(d2FdL(i,:,:),dzdy))
            end do
         endif
         if (present(basisxy)) then
            basisxy(ik+r,:) = forbasisxy*Fr1r2 + &
                       forbasisx*(dFdL1*dzdy(1)+dFdL2*dzdy(2)+dFdL3*dzdy(3)) + &
                       forbasisy*(dFdL1*dzdx(1)+dFdL2*dzdx(2)+dFdL3*dzdx(3))
            do i=1,size(x)
               basisxy(ik+r,i) = basisxy(ik+r,i) + &
                         zetaprod(i)*dot_product(dzdx,matmul(d2FdL(i,:,:),dzdy))
            end do
         endif
      endif
   end do
end do

if (use_basis == SZABO) then
   deallocate(legendre_d0,legendre2_d0)
   if (present(basisx) .or. present(basisy) .or. present(basisxx) .or. &
       present(basisyy) .or. present(basisxy)) then
      deallocate(legendre_d1,legendre2_d1)
   endif
   if (present(basisxx) .or. present(basisyy) .or. present(basisxy)) then
      deallocate(legendre_d2,legendre2_d2)
   endif
endif

endif ! degree(4) > 2

if (present(basis)) basis = lbasis
deallocate(lbasis)

end subroutine p_hier_basis_func_a

!          -------------------
subroutine p_hier_basis_func_s(x,y,xvert,yvert,degree,set,basis,basisx,basisy,&
                               basisxx,basisyy,basisxy)
!          -------------------

!----------------------------------------------------
! This routine is a version of p_hier_basis_func for a single point.  It turns
! it into an array of size 1 and calls the version for an array of points.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x, y, xvert(3), yvert(3)
integer, intent(in) :: degree(4)
character(len=1), intent(in) :: set
real(my_real), intent(out), optional :: basis(:), basisx(:), basisy(:), &
                                        basisxx(:), basisyy(:), basisxy(:)
!----------------------------------------------------
! Local variables:

real(my_real) :: x_a(1), y_a(1)
real(my_real), allocatable :: basis_a(:,:), basisx_a(:,:), basisy_a(:,:), &
                              basisxx_a(:,:), basisyy_a(:,:), basisxy_a(:,:)
integer :: code, astat1, astat2, astat3
!----------------------------------------------------
! Begin executable code

x_a = x
y_a = y

code = 0
astat1 = 0
astat2 = 0
astat3 = 0

if (present(basis)) then
   allocate(basis_a(size(basis),1),stat=astat1)
   code = code + 1
endif

if (present(basisx) .and. astat1 == 0) then
   allocate(basisx_a(size(basisx),1),stat=astat2)
   code = code + 2
endif

if (present(basisy) .and. astat1 == 0 .and. astat2 == 0) then
   allocate(basisy_a(size(basisy),1),stat=astat3)
   code = code + 4
endif

if (astat1 /= 0 .or. astat2 /= 0 .or. astat3 /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in p_hier_basis_func_s")
   stop
endif

if (present(basisxx)) then
   allocate(basisxx_a(size(basisxx),1),stat=astat1)
   code = code + 8
endif

if (present(basisyy) .and. astat1 == 0) then
   allocate(basisyy_a(size(basisyy),1),stat=astat2)
   code = code + 16
endif

if (present(basisxy) .and. astat1 == 0 .and. astat2 == 0) then
   allocate(basisxy_a(size(basisxy),1),stat=astat3)
   code = code + 32
endif

if (astat1 /= 0 .or. astat2 /= 0 .or. astat3 /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in p_hier_basis_func_s")
   stop
endif

select case (code)
case(0)
case(1)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basis=basis_a)
case(2)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basisx=basisx_a)
case(3)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basis=basis_a, &
                            basisx=basisx_a)
case(4)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basisy=basisy_a)
case(5)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basis=basis_a, &
                            basisy=basisy_a)
case(6)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basisx=basisx_a, &
                            basisy=basisy_a)
case(7)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basis=basis_a, &
                            basisx=basisx_a,basisy=basisy_a)
case(8)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basisxx=basisxx_a)
case(9)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basis=basis_a, &
                            basisxx=basisxx_a)
case(10)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basisx=basisx_a, &
                            basisxx=basisxx_a)
case(11)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basis=basis_a, &
                            basisx=basisx_a,basisxx=basisxx_a)
case(12)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basisy=basisy_a, &
                            basisxx=basisxx_a)
case(13)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basis=basis_a, &
                            basisy=basisy_a,basisxx=basisxx_a)
case(14)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basisx=basisx_a, &
                            basisy=basisy_a,basisxx=basisxx_a)
case(15)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basis=basis_a, &
                            basisx=basisx_a,basisy=basisy_a,basisxx=basisxx_a)
case(16)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basisyy=basisyy_a)
case(17)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basis=basis_a, &
                            basisyy=basisyy_a)
case(18)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basisx=basisx_a, &
                            basisyy=basisyy_a)
case(19)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basis=basis_a, &
                            basisx=basisx_a,basisyy=basisyy_a)
case(20)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basisy=basisy_a, &
                            basisyy=basisyy_a)
case(21)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basis=basis_a, &
                            basisy=basisy_a,basisyy=basisyy_a)
case(22)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basisx=basisx_a, &
                            basisy=basisy_a,basisyy=basisyy_a)
case(23)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basis=basis_a, &
                            basisx=basisx_a,basisy=basisy_a,basisyy=basisyy_a)
case(24)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basisxx=basisxx_a, &
                            basisyy=basisyy_a)
case(25)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basis=basis_a, &
                            basisxx=basisxx_a,basisyy=basisyy_a)
case(26)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basisx=basisx_a, &
                            basisxx=basisxx_a,basisyy=basisyy_a)
case(27)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basis=basis_a, &
                            basisx=basisx_a,basisxx=basisxx_a,basisyy=basisyy_a)
case(28)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basisy=basisy_a, &
                            basisxx=basisxx_a,basisyy=basisyy_a)
case(29)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basis=basis_a, &
                            basisy=basisy_a,basisxx=basisxx_a,basisyy=basisyy_a)
case(30)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basisx=basisx_a, &
                            basisy=basisy_a,basisxx=basisxx_a,basisyy=basisyy_a)
case(31)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basis=basis_a, &
                            basisx=basisx_a,basisy=basisy_a,basisxx=basisxx_a, &
                            basisyy=basisyy_a)
case(32)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basisxy=basisxy_a)
case(33)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basis=basis_a, &
                            basisxy=basisxy_a)
case(34)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basisx=basisx_a, &
                            basisxy=basisxy_a)
case(35)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basis=basis_a, &
                            basisx=basisx_a,basisxy=basisxy_a)
case(36)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basisy=basisy_a, &
                            basisxy=basisxy_a)
case(37)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basis=basis_a, &
                            basisy=basisy_a,basisxy=basisxy_a)
case(38)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basisx=basisx_a, &
                            basisy=basisy_a,basisxy=basisxy_a)
case(39)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basis=basis_a, &
                            basisx=basisx_a,basisy=basisy_a,basisxy=basisxy_a)
case(40)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basisxx=basisxx_a, &
                            basisxy=basisxy_a)
case(41)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basis=basis_a, &
                            basisxx=basisxx_a,basisxy=basisxy_a)
case(42)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basisx=basisx_a, &
                            basisxx=basisxx_a,basisxy=basisxy_a)
case(43)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basis=basis_a, &
                            basisx=basisx_a,basisxx=basisxx_a,basisxy=basisxy_a)
case(44)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basisy=basisy_a, &
                            basisxx=basisxx_a,basisxy=basisxy_a)
case(45)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basis=basis_a, &
                            basisy=basisy_a,basisxx=basisxx_a,basisxy=basisxy_a)
case(46)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basisx=basisx_a, &
                            basisy=basisy_a,basisxx=basisxx_a,basisxy=basisxy_a)
case(47)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basis=basis_a, &
                            basisx=basisx_a,basisy=basisy_a,basisxx=basisxx_a, &
                            basisxy=basisxy_a)
case(48)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basisyy=basisyy_a, &
                            basisxy=basisxy_a)
case(49)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basis=basis_a, &
                            basisyy=basisyy_a,basisxy=basisxy_a)
case(50)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basisx=basisx_a, &
                            basisyy=basisyy_a,basisxy=basisxy_a)
case(51)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basis=basis_a, &
                            basisx=basisx_a,basisyy=basisyy_a,basisxy=basisxy_a)
case(52)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basisy=basisy_a, &
                            basisyy=basisyy_a,basisxy=basisxy_a)
case(53)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basis=basis_a, &
                            basisy=basisy_a,basisyy=basisyy_a,basisxy=basisxy_a)
case(54)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basisx=basisx_a, &
                            basisy=basisy_a,basisyy=basisyy_a,basisxy=basisxy_a)
case(55)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basis=basis_a, &
                            basisx=basisx_a,basisy=basisy_a,basisyy=basisyy_a, &
                            basisxy=basisxy_a)
case(56)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basisxx=basisxx_a, &
                            basisyy=basisyy_a,basisxy=basisxy_a)
case(57)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basis=basis_a, &
                            basisxx=basisxx_a,basisyy=basisyy_a, &
                            basisxy=basisxy_a)
case(58)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basisx=basisx_a, &
                            basisxx=basisxx_a,basisyy=basisyy_a, &
                            basisxy=basisxy_a)
case(59)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basis=basis_a, &
                            basisx=basisx_a,basisxx=basisxx_a, &
                            basisyy=basisyy_a,basisxy=basisxy_a)
case(60)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basisy=basisy_a, &
                            basisxx=basisxx_a,basisyy=basisyy_a, &
                            basisxy=basisxy_a)
case(61)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basis=basis_a, &
                            basisy=basisy_a,basisxx=basisxx_a, &
                            basisyy=basisyy_a,basisxy=basisxy_a)
case(62)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basisx=basisx_a, &
                            basisy=basisy_a,basisxx=basisxx_a, &
                            basisyy=basisyy_a,basisxy=basisxy_a)
case(63)
   call p_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basis=basis_a, &
                            basisx=basisx_a,basisy=basisy_a,basisxx=basisxx_a, &
                            basisyy=basisyy_a,basisxy=basisxy_a)
end select

if (present(basis) ) basis  = basis_a(:,1)
if (present(basisx)) basisx = basisx_a(:,1)
if (present(basisy)) basisy = basisy_a(:,1)
if (present(basisxx)) basisxx = basisxx_a(:,1)
if (present(basisyy)) basisyy = basisyy_a(:,1)
if (present(basisxy)) basisxy = basisxy_a(:,1)

end subroutine p_hier_basis_func_s

!          -----------
subroutine barycentric(x,y,xvert,yvert,zeta,dzdx,dzdy)
!          -----------

!----------------------------------------------------
! This subroutine returns the barycentric coordinates of each (x,y) and
! their derivatives.
! RESTRICTION triangles
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x(:),y(:),xvert(:),yvert(:)
real(my_real), intent(out) :: zeta(:,:),dzdx(:),dzdy(:)

!----------------------------------------------------
! Local variables:

real(quad_real) :: x1,x2,x3,y1,y2,y3,det,ssump,ssumm
real(quad_real), dimension(size(x)) :: xy1,xy2,xy3,yx1,yx2,yx3,sump,summ
real(quad_real) :: x1y2,x2y3,x3y1,x1y3,x2y1,x3y2

!----------------------------------------------------
! Begin executable code

! put the vertices in local variables to shorten expressions

x1 = real(xvert(1),quad_real)
x2 = real(xvert(2),quad_real)
x3 = real(xvert(3),quad_real)
y1 = real(yvert(1),quad_real)
y2 = real(yvert(2),quad_real)
y3 = real(yvert(3),quad_real)

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

ssump=0.0_quad_real; ssumm=0.0_quad_real
if (x1y2 > 0.0_quad_real) then; ssump=ssump+x1y2; else; ssumm=ssumm+x1y2; endif
if (x1y3 < 0.0_quad_real) then; ssump=ssump-x1y3; else; ssumm=ssumm-x1y3; endif
if (x2y3 > 0.0_quad_real) then; ssump=ssump+x2y3; else; ssumm=ssumm+x2y3; endif
if (x2y1 < 0.0_quad_real) then; ssump=ssump-x2y1; else; ssumm=ssumm-x2y1; endif
if (x3y1 > 0.0_quad_real) then; ssump=ssump+x3y1; else; ssumm=ssumm+x3y1; endif
if (x3y2 < 0.0_quad_real) then; ssump=ssump-x3y2; else; ssumm=ssumm-x3y2; endif
det = ssump + ssumm

! compute the coordinates

! zeta(1) = (x*(y2-y3) + y*(x3-x2) + (x2*y3-x3*y2))/det

sump=0.0_quad_real; summ=0.0_quad_real
where (xy2 > 0.0_quad_real)
   sump=sump+xy2
elsewhere
   summ=summ+xy2
endwhere
where (xy3 < 0.0_quad_real)
   sump=sump-xy3
elsewhere
   summ=summ-xy3
endwhere
where (yx3 > 0.0_quad_real)
   sump=sump+yx3
elsewhere
   summ=summ+yx3
endwhere
where (yx2 < 0.0_quad_real)
   sump=sump-yx2
elsewhere
   summ=summ-yx2
endwhere
if (x2y3 > 0.0_quad_real) then; sump=sump+x2y3; else; summ=summ+x2y3; endif
if (x3y2 < 0.0_quad_real) then; sump=sump-x3y2; else; summ=summ-x3y2; endif
zeta(:,1) = sump/det + summ/det

! zeta(2) = (x*(y3-y1) + y*(x1-x3) + (x3*y1-x1*y3))/det

sump=0.0_quad_real; summ=0.0_quad_real
where (xy3 > 0.0_quad_real)
   sump=sump+xy3
elsewhere
   summ=summ+xy3
endwhere
where (xy1 < 0.0_quad_real)
   sump=sump-xy1
elsewhere
   summ=summ-xy1
endwhere
where (yx1 > 0.0_quad_real)
   sump=sump+yx1
elsewhere
   summ=summ+yx1
endwhere
where (yx3 < 0.0_quad_real)
   sump=sump-yx3
elsewhere
   summ=summ-yx3
endwhere
if (x3y1 > 0.0_quad_real) then; sump=sump+x3y1; else; summ=summ+x3y1; endif
if (x1y3 < 0.0_quad_real) then; sump=sump-x1y3; else; summ=summ-x1y3; endif
zeta(:,2) = sump/det + summ/det

! zeta(3) = (x*(y1-y2) + y*(x2-x1) + (x1*y2-x2*y1))/det

sump=0.0_quad_real; summ=0.0_quad_real
where (xy1 > 0.0_quad_real)
   sump=sump+xy1
elsewhere
   summ=summ+xy1
endwhere
where (xy2 < 0.0_quad_real)
   sump=sump-xy2
elsewhere
   summ=summ-xy2
endwhere
where (yx2 > 0.0_quad_real)
   sump=sump+yx2
elsewhere
   summ=summ+yx2
endwhere
where (yx1 < 0.0_quad_real)
   sump=sump-yx1
elsewhere
   summ=summ-yx1
endwhere
if (x1y2 > 0.0_quad_real) then; sump=sump+x1y2; else; summ=summ+x1y2; endif
if (x2y1 < 0.0_quad_real) then; sump=sump-x2y1; else; summ=summ-x2y1; endif
zeta(:,3) = sump/det + summ/det

! derivatives

dzdx(1) = (y2-y3)/det
dzdx(2) = (y3-y1)/det
dzdx(3) = (y1-y2)/det
dzdy(1) = (x3-x2)/det
dzdy(2) = (x1-x3)/det
dzdy(3) = (x2-x1)/det

end subroutine barycentric

!          ---------
subroutine dlegendre(x,degree,deriv0,deriv1,deriv2,deriv3)
!          ---------

!----------------------------------------------------
! This routine optionally returns the Legendre polynomials and first and
! second and third derivatives from degree 1 to degree.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x(:)
integer, intent(in) :: degree
real(my_real), intent(out), optional :: deriv0(size(x),degree), &
                            deriv1(size(x),degree), deriv2(size(x),degree), &
                            deriv3(size(x),degree)
!----------------------------------------------------
! Local variables:

real(my_real) :: legendre(size(x),0:degree)
real(my_real) :: c1,c2
integer :: d
!----------------------------------------------------
! Begin executable code

if (degree == 0) return

legendre(:,0) = 1
legendre(:,1) = x
if (present(deriv1)) deriv1(:,1) = 1
if (present(deriv2)) deriv2(:,1) = 0
if (present(deriv3)) deriv3(:,1) = 0
if (degree > 1) then
   legendre(:,2) = (3*x*x-1)/2
   if (present(deriv1)) deriv1(:,2) = 3*x
   if (present(deriv2)) deriv2(:,2) = 3
   if (present(deriv3)) deriv3(:,2) = 0
endif

do d=3,degree
   c1 = (2*(d-1.0_my_real)+1)/d
   c2 = (d-1.0_my_real)/d
   legendre(:,d) = c1*x*legendre(:,d-1)
   if (my_real == kind(0.0)) then
      call saxpy(size(x),-c2,legendre(1,d-2),1,legendre(1,d),1)
      if (present(deriv1)) then
         deriv1(:,d) = c1*x*deriv1(:,d-1)
         call saxpy(size(x),c1,legendre(1,d-1),1,deriv1(1,d),1)
         call saxpy(size(x),-c2,deriv1(1,d-2),1,deriv1(1,d),1)
      endif
      if (present(deriv2)) then
         deriv2(:,d) = c1*x*deriv2(:,d-1)
         call saxpy(size(x),2*c1,deriv1(1,d-1),1,deriv2(1,d),1)
         call saxpy(size(x),-c2,deriv2(1,d-2),1,deriv2(1,d),1)
      endif
      if (present(deriv3)) then
         deriv3(:,d) = c1*x*deriv3(:,d-1)
         call saxpy(size(x),3*c1,deriv2(1,d-1),1,deriv3(1,d),1)
         call saxpy(size(x),-c2,deriv3(1,d-2),1,deriv3(1,d),1)
      endif
   elseif (my_real == kind(0.0d0)) then
      call daxpy(size(x),-c2,legendre(1,d-2),1,legendre(1,d),1)
      if (present(deriv1)) then
         deriv1(:,d) = c1*x*deriv1(:,d-1)
         call daxpy(size(x),c1,legendre(1,d-1),1,deriv1(1,d),1)
         call daxpy(size(x),-c2,deriv1(1,d-2),1,deriv1(1,d),1)
      endif
      if (present(deriv2)) then
         deriv2(:,d) = c1*x*deriv2(:,d-1)
         call daxpy(size(x),2*c1,deriv1(1,d-1),1,deriv2(1,d),1)
         call daxpy(size(x),-c2,deriv2(1,d-2),1,deriv2(1,d),1)
      endif
      if (present(deriv3)) then
         deriv3(:,d) = c1*x*deriv3(:,d-1)
         call daxpy(size(x),3*c1,deriv2(1,d-1),1,deriv3(1,d),1)
         call daxpy(size(x),-c2,deriv3(1,d-2),1,deriv3(1,d),1)
      endif
   else
      ierr = PHAML_INTERNAL_ERROR
      call fatal("need single or double precision for BLAS")
      stop
   endif
end do

if (present(deriv0)) deriv0 = legendre(:,1:)

end subroutine dlegendre

!        ---------
function factorial(n)
!        ---------

!----------------------------------------------------
! This routine computes n!
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: n
integer :: factorial
!----------------------------------------------------
! Local variables:

integer :: i
!----------------------------------------------------
! Begin executable code

if (n < 0) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("attempt to compute factorial of negative number",intlist=(/n/))
   stop
endif

factorial = 1
do i=2,n
   factorial = factorial*i
end do

end function factorial

!-------------------------------------------------------
! ROUTINES FOR NODAL BASIS
!-------------------------------------------------------

!          ------------------
subroutine nodal_basis_func_a(x,y,xvert,yvert,degree,set,basis,basisx,basisy)
!          ------------------

!----------------------------------------------------
! This routine computes nodal basis functions and derivatives at (x,y)
! for the triangle given by (xvert,yvert).  basis, basisx and basisy are all
! optional; only those present are computed.  degree determines the polynomial
! degree of the basis functions;
!
! set determines which subset of the bases is returned: "a" requests all
! basis functions, "b" requests only the black bases, and "r" the red bases.
!
!  The order in which the bases are returned:
!  - those along the edge from vertex 1 to vertex 2
!  - subsequent lines of nodes going towards vertex 3
!  For example, for cubics
! vertex 1  1
!           |\
!           2 5
!           |  \
!           3 6 8
!           |    \
! vertex 2  4-7-9-10  vertex 3
!
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x(:), y(:), xvert(3), yvert(3)
integer, intent(in) :: degree
character(len=1), intent(in) :: set
real(my_real), intent(out), optional :: basis(:,:), basisx(:,:), basisy(:,:)
!----------------------------------------------------
! Local variables:

real(my_real) :: zeta(size(x),3), dzdx(3), dzdy(3), b(size(x)), bx(size(x)), &
                 by(size(x)), factor(size(x)), nf
integer :: i, j, bas, nbasis, npoint, node_bary(3)
!----------------------------------------------------
! Begin executable code

! check that the size of the return arrays is big enough

npoint = size(x)
select case(set)
case("a")
   nbasis = ((degree+1)*(degree+2))/2
case("b")
   if (2*(degree/2) == degree) then
      nbasis = ((degree+2)*(degree+2))/4
   else
      nbasis = ((degree+1)*(degree+3))/4
   endif
case("r")
   if (2*(degree/2) == degree) then
      nbasis = (degree*(degree+2))/4
   else
      nbasis = ((degree+1)*(degree+1))/4
   endif
end select

if (present(basis)) then
   if (size(basis,1) < nbasis .or. size(basis,2) < npoint) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("Size of basis array is too small in nodal_basis_func.", &
                 intlist = (/size(basis,1),size(basis,2),nbasis,npoint/))
      stop
   endif
endif

if (present(basisx)) then
   if (size(basisx,1) < nbasis .or. size(basisx,2) < npoint) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("Size of basisx array is too small in nodal_basis_func.", &
                 intlist = (/size(basisx,1),size(basisx,2),nbasis,npoint/))
      stop
   endif
endif

if (present(basisy)) then
   if (size(basisy,1) < nbasis .or. size(basisy,2) < npoint) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("Size of basisy array is too small in nodal_basis_func.", &
                 intlist = (/size(basisy,1),size(basisy,2),nbasis,npoint/))
      stop
   endif
endif

! compute the barycentric coordinates of (x,y)

call barycentric(x,y,xvert,yvert,zeta,dzdx,dzdy)

! for each basis function, characterized by the node barycentric coordinates
! scaled by degree

if (set == "r") then
   node_bary(1) = degree
   node_bary(2) = -1
   node_bary(3) = 1
else
   node_bary(1) = degree + 1
   node_bary(2) = -1
   node_bary(3) = 0
endif

do bas=1,nbasis

! barycentric coordinates of node of next basis

   node_bary(1) = node_bary(1) - 1
   node_bary(2) = node_bary(2) + 1
   if (node_bary(1) < 0) then
      if (set == "a") then
         node_bary(3) = node_bary(3) + 1
      else
         node_bary(3) = node_bary(3) + 2
      endif
      node_bary(2) = 0
      node_bary(1) = degree - node_bary(3)
   endif

! compute basis function and derivatives

   b = 1.0_my_real
   if (present(basisx)) bx = 0.0_my_real
   if (present(basisy)) by = 0.0_my_real
   nf = 1.0_my_real

   do i=1,3
      do j=0,degree-1
         if (j < node_bary(i)) then
            nf = nf*real(node_bary(i)-j,my_real)/degree
            factor = zeta(:,i) - real(j,my_real)/degree
            if (present(basisx)) bx = bx*factor + b*dzdx(i)
            if (present(basisy)) by = by*factor + b*dzdy(i)
            b = b*factor
         endif
      end do
   end do

   if (present(basis)) basis(bas,:) = b/nf
   if (present(basisx)) basisx(bas,:) = bx/nf
   if (present(basisy)) basisy(bas,:) = by/nf

end do

end subroutine nodal_basis_func_a

!          ------------------
subroutine nodal_basis_func_s(x,y,xvert,yvert,degree,set,basis,basisx,basisy)
!          ------------------

!----------------------------------------------------
! This routine is a version of nodal_basis_func for a single point.  It turns
! it into an array of size 1 and calls the version for an array of points.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x, y, xvert(3), yvert(3)
integer, intent(in) :: degree
character(len=1), intent(in) :: set
real(my_real), intent(out), optional :: basis(:), basisx(:), basisy(:)
!----------------------------------------------------
! Local variables:

real(my_real) :: x_a(1), y_a(1)
real(my_real), allocatable :: basis_a(:,:), basisx_a(:,:), basisy_a(:,:)
integer :: code, astat1, astat2, astat3
!----------------------------------------------------
! Begin executable code

x_a = x
y_a = y

code = 0
astat1 = 0
astat2 = 0
astat3 = 0

if (present(basis)) then
   allocate(basis_a(size(basis),1),stat=astat1)
   code = code + 1
endif

if (present(basisx) .and. astat1 == 0) then
   allocate(basisx_a(size(basisx),1),stat=astat2)
   code = code + 2
endif

if (present(basisy) .and. astat1 == 0 .and. astat2 == 0) then
   allocate(basisy_a(size(basisy),1),stat=astat3)
   code = code + 4
endif

if (astat1 /= 0 .or. astat2 /= 0 .or. astat3 /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in nodal_basis_func_s")
   stop
endif

select case (code)
case(0)
case(1)
   call nodal_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basis_a)
case(2)
   call nodal_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basisx=basisx_a)
case(3)
   call nodal_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basis_a,basisx_a)
case(4)
   call nodal_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basisy=basisy_a)
case(5)
   call nodal_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basis_a,basisy=basisy_a)
case(6)
   call nodal_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basisx_a,basisy=basisy_a)
case(7)
   call nodal_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basis_a,basisx_a,basisy_a)
end select

if (present(basis) ) basis  = basis_a(:,1)
if (present(basisx)) basisx = basisx_a(:,1)
if (present(basisy)) basisy = basisy_a(:,1)

end subroutine nodal_basis_func_s

!-------------------------------------------------------
! ROUTINES FOR H-HIERARCHICAL BASIS
!-------------------------------------------------------

!          -------------------
subroutine h_hier_basis_func_a(x,y,xvert,yvert,degree,set,basis,basisx,basisy)
!          -------------------

!----------------------------------------------------
! This routine computes h-hierarchical basis functions and derivatives at (x,y)
! for the triangle given by (xvert,yvert).  basis, basisx and basisy are all
! optional; only those present are computed.  degree determines the polynomial
! degree of the basis functions;
!
! set determines which subset of the bases is returned: "a" requests all
! basis functions, "b" requests only the black bases, and "r" the red bases.
!
!  The order in which the bases are returned:
!  - those along the edge from vertex 1 to vertex 2
!  - subsequent lines of black nodes (every other line) going towards vertex 3
!  - the red nodes going along lines parallel to the line from vertex 1 to 2
!  For example, for cubics
! vertex 1  1
!           |\
!           2 7
!           |  \
!           3 8 5
!           |    \
! vertex 2  4-9-6-10  vertex 3
!
! vertex 3 is assumed to be the newest node from a bisection, hence the
! parent triangle would have vertex 1 and 2, and a vertex that is twice
! the distance from vertex 1 as vertex 3
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x(:), y(:), xvert(3), yvert(3)
integer, intent(in) :: degree
character(len=1), intent(in) :: set
real(my_real), intent(out), optional :: basis(:,:), basisx(:,:), basisy(:,:)
!----------------------------------------------------
! Local variables:

real(my_real) :: xvert_par(3), yvert_par(3), &
                 basis_loc(((degree+1)*(degree+2))/2,size(x)), &
                 basisx_loc(((degree+1)*(degree+2))/2,size(x)), &
                 basisy_loc(((degree+1)*(degree+2))/2,size(x))
integer :: i, j, npoint, nbasis, nblack, nred, code, n1, n2
!----------------------------------------------------
! Begin executable code

! number of red and black nodes

if (2*(degree/2) == degree) then
   nblack = ((degree+2)*(degree+2))/4
   nred = (degree*(degree+2))/4
else
   nblack = ((degree+1)*(degree+3))/4
   nred = ((degree+1)*(degree+1))/4
endif

! check that the size of the return arrays is big enough

npoint = size(x)
select case(set)
case("a")
   nbasis = ((degree+1)*(degree+2))/2
case("b")
   nbasis = nblack
case("r")
   nbasis = nred
end select

if (present(basis)) then
   if (size(basis,1) < nbasis .or. size(basis,2) < npoint) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("Size of basis array is too small in h_hier_basis_func.", &
                 intlist = (/size(basis,1),size(basis,2),nbasis,npoint/))
      stop
   endif
endif

if (present(basisx)) then
   if (size(basisx,1) < nbasis .or. size(basisx,2) < npoint) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("Size of basisx array is too small in h_hier_basis_func.", &
                 intlist = (/size(basisx,1),size(basisx,2),nbasis,npoint/))
      stop
   endif
endif

if (present(basisy)) then
   if (size(basisy,1) < nbasis .or. size(basisy,2) < npoint) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("Size of basisy array is too small in h_hier_basis_func.", &
                 intlist = (/size(basisy,1),size(basisy,2),nbasis,npoint/))
      stop
   endif
endif

! set the code that determines which basis arguments are passed to nodal_basis

code = 0
if (present(basis)) code = code + 1
if (present(basisx)) code = code + 2
if (present(basisy)) code = code + 4

if (set == "a" .or. set == "b") then

! to get the bases associated with black nodes, get the nodal bases from
! the parent triangle

   xvert_par = xvert
   yvert_par = yvert
   xvert_par(3) = 2*xvert(3) - xvert(1)
   yvert_par(3) = 2*yvert(3) - yvert(1)

   select case(code)
   case(0)
   case (1)
      call nodal_basis_func(x,y,xvert_par,yvert_par,degree,"a",basis_loc)
   case (2)
      call nodal_basis_func(x,y,xvert_par,yvert_par,degree,"a", &
                            basisx=basisx_loc)
   case (3)
      call nodal_basis_func(x,y,xvert_par,yvert_par,degree,"a",basis_loc, &
                            basisx_loc)
   case (4)
      call nodal_basis_func(x,y,xvert_par,yvert_par,degree,"a", &
                            basisy=basisy_loc)
   case (5)
      call nodal_basis_func(x,y,xvert_par,yvert_par,degree,"a",basis_loc, &
                            basisy=basisy_loc)
   case (6)
      call nodal_basis_func(x,y,xvert_par,yvert_par,degree,"a", &
                            basisx=basisx_loc,basisy=basisy_loc)
   case (7)
      call nodal_basis_func(x,y,xvert_par,yvert_par,degree,"a",basis_loc, &
                            basisx_loc,basisy_loc)
   end select

! extract those that are in this triangle

   n1 = 0
   n2 = 0

   do i=1,(degree+2)/2
      do j=1,degree+3-2*i
         n1 = n1+1
         n2 = n2+1
         if (present(basis)) basis(n1,1:size(x)) = basis_loc(n2,:)
         if (present(basisx)) basisx(n1,1:size(x)) = basisx_loc(n2,:)
         if (present(basisy)) basisy(n1,1:size(x)) = basisy_loc(n2,:)
      end do
      n2 = n2 + i-1
   end do

endif

if (set == "a" .or. set == "r") then

! to get the bases associated with red nodes, get the nodal bases from
! for the red nodes

   select case(code)
   case(0)
   case (1)
      call nodal_basis_func(x,y,xvert,yvert,degree,"r",basis_loc)
   case (2)
      call nodal_basis_func(x,y,xvert,yvert,degree,"r",basisx=basisx_loc)
   case (3)
      call nodal_basis_func(x,y,xvert,yvert,degree,"r",basis_loc,basisx_loc)
   case (4)
      call nodal_basis_func(x,y,xvert,yvert,degree,"r",basisy=basisy_loc)
   case (5)
      call nodal_basis_func(x,y,xvert,yvert,degree,"r",basis_loc, &
                            basisy=basisy_loc)
   case (6)
      call nodal_basis_func(x,y,xvert,yvert,degree,"r", &
                            basisx=basisx_loc,basisy=basisy_loc)
   case (7)
      call nodal_basis_func(x,y,xvert,yvert,degree,"r",basis_loc, &
                            basisx_loc,basisy_loc)
   end select

! copy them to the result arrays

   if (set == "a") then
      n1 = nblack
   else
      n1 = 0
   endif

   do i=1,nred
      if (present(basis)) basis(i+n1,1:size(x)) = basis_loc(i,:)
      if (present(basisx)) basisx(i+n1,1:size(x)) = basisx_loc(i,:)
      if (present(basisy)) basisy(i+n1,1:size(x)) = basisy_loc(i,:)
   end do

end if

end subroutine h_hier_basis_func_a

!          -------------------
subroutine h_hier_basis_func_s(x,y,xvert,yvert,degree,set,basis,basisx,basisy)
!          -------------------

!----------------------------------------------------
! This routine is a version of h_hier_basis_func for a single point.  It turns
! it into an array of size 1 and calls the version for an array of points.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x, y, xvert(3), yvert(3)
integer, intent(in) :: degree
character(len=1), intent(in) :: set
real(my_real), intent(out), optional :: basis(:), basisx(:), basisy(:)
!----------------------------------------------------
! Local variables:

real(my_real) :: x_a(1), y_a(1)
real(my_real), allocatable :: basis_a(:,:), basisx_a(:,:), basisy_a(:,:)
integer :: code, astat1, astat2, astat3
!----------------------------------------------------
! Begin executable code

x_a = x
y_a = y

code = 0
astat1 = 0
astat2 = 0
astat3 = 0

if (present(basis)) then
   allocate(basis_a(size(basis),1),stat=astat1)
   code = code + 1
endif

if (present(basisx) .and. astat1 == 0) then
   allocate(basisx_a(size(basisx),1),stat=astat2)
   code = code + 2
endif

if (present(basisy) .and. astat1 == 0 .and. astat2 == 0) then
   allocate(basisy_a(size(basisy),1),stat=astat3)
   code = code + 4
endif

if (astat1 /= 0 .or. astat2 /= 0 .or. astat3 /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in nodal_basis_func_s")
   stop
endif

select case (code)
case(0)
case(1)
   call h_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basis_a)
case(2)
   call h_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basisx=basisx_a)
case(3)
   call h_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basis_a,basisx_a)
case(4)
   call h_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basisy=basisy_a)
case(5)
   call h_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basis_a,basisy=basisy_a)
case(6)
   call h_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basisx_a,basisy=basisy_a)
case(7)
   call h_hier_basis_func_a(x_a,y_a,xvert,yvert,degree,set,basis_a,basisx_a,basisy_a)
end select

if (present(basis) ) basis  = basis_a(:,1)
if (present(basisx)) basisx = basisx_a(:,1)
if (present(basisy)) basisy = basisy_a(:,1)

end subroutine h_hier_basis_func_s

end module basis_functions
