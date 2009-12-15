module quadrature_rules

!----------------------------------------------------
! This module contains quadrature rules for integration over
!  - a line in the plane, up to order 24
!  - a triangle in the plane, up to order 45
!----------------------------------------------------

!----------------------------------------------------
! Quadrature rules for the line in the plane are taken from
! Abramowitz, M. and Stegun, I.A., "Handbook of Mathematical Functions",
!   Dover Publications, New York, ninth printing (1972).
! Some of the rules are not given (order=11,13,14,15,17,18,19,21,22,23).  In
! these cases, the first available larger order is returned.
! Integration of x**p over the unit interval shows the k order rule to be
! exact for polynomials up to degree 2k-1 and the error to be O(h^(2k))
!----------------------------------------------------

!----------------------------------------------------
! Quadrature rules for the triangle in the plane up to order 20 are taken from
! Dunavant, D.A., "High Degree Efficient Symmetrical Gaussian Quadrature
!   Rules for the Triangle", Int. J. Num. Meth. Eng., 21 (1985) pp 1129-1148.
! Integration of x**p + y**p over the unit square shows the k order rule to
! be exact for polynomials up to degree k (k+1 if k is even) and the error
! to be O(h^(k+1)) (k+2 if k is even).
!
! Quadrature rules for the triangle in the plane of order greater than 20 are
! computed by a product rule from the rules for lines along with the Duffy
! transformation from the square to the triangle.  The order of the line
! quadrature rule is chosen to give at least the same order as the trend of
! the Dunavant rules.  Odd degree rules are exact for polynomials up to
! degree k+1 and the error is O(h^(k+2)).  Even degree rules are exact for
! polynomials up to degree k+2 and the error is O(h^(k+3)).  Some of the line
! quadrature rules are not given.  In these cases, the first available
! larger order is returned.

!----------------------------------------------------
! The public procedures are:

! These routines take the order of the quadrature rule and the vertices
! that define the domain of integration as input.  They return the
! number of quadrature points, the weights, the quadrature points, and
! an error status.  The size of the domain of integration is included in
! the weights, so the integral is given by
!
! sum i=1,nqpoints [qweights(i) * integrand(xquad(i),...)]
!
! The weights and quadrature points are allocated in these routines and
! should be deallocated by the calling routine.
!
! The error status, ierr is:
!   0 - sucess
!   1 - no rule for quadrature of the requested order
!   2 - memory allocation failed
!
!---
!  subroutine quadrature_rule_line(order,xvert,yvert,nqpoints,qweights, &
!                                  xquad,yquad,ierr)
!  integer, intent(in) :: order
!  real(RKIND), intent(in) :: xvert(2),yvert(2)
!  integer, intent(out) :: nqpoints
!  real(RKIND), pointer :: qweights(:)
!  real(RKIND), pointer :: xquad(:),yquad(:)
!  integer, intent(out) :: ierr
!  end subroutine quadrature_rule_line
!
!---
!  subroutine quadrature_rule_tri(order,xvert,yvert,nqpoints,qweights, &
!                                 xquad,yquad,ierr,stay_in)
!  integer, intent(in) :: order
!  real(RKIND), intent(in) :: xvert(3),yvert(3)
!  integer, intent(out) :: nqpoints
!  real(RKIND), pointer :: qweights(:)
!  real(RKIND), pointer :: xquad(:),yquad(:)
!  integer, intent(out) :: ierr
!  logical, intent(in), optional :: stay_in
!  end subroutine quadrature_rule_tri
!
!----------------------------------------------------
! Other modules used are:

use global, only: my_real
!----------------------------------------------------

implicit none
private
public RKIND, MAX_QUAD_ORDER_LINE, MAX_QUAD_ORDER_TRI, &
       quadrature_rule_line, quadrature_rule_tri

!----------------------------------------------------
! The following parameters are defined:

! kind for real dummy arguments

integer, parameter :: RKIND = my_real

! maximum order of quadrature rule currently supported

integer, parameter :: MAX_QUAD_ORDER_LINE = 24, &
                      MAX_QUAD_ORDER_TRI  = 45, &
                      MAX_DUNAVANT_ORDER  = 20
!----------------------------------------------------

!----------------------------------------------------
! The following types are defined:

type double_arrays
   real(kind(0.0d0)), pointer :: val(:)
end type double_arrays

!----------------------------------------------------
! The following variables are defined:

! number of Gauss quadrature points

integer, save :: ng_line(MAX_QUAD_ORDER_LINE), &
                 ng_tri(MAX_DUNAVANT_ORDER)

! quadrature weights

type(double_arrays), save :: weight_line(MAX_QUAD_ORDER_LINE), &
                             weight_tri(MAX_DUNAVANT_ORDER)

! quadrature points, in barycentric coordinates

type(double_arrays), save ::  alpha_line(MAX_QUAD_ORDER_LINE), &
                              alpha_tri(MAX_DUNAVANT_ORDER),   &
                              beta_tri(MAX_DUNAVANT_ORDER),    &
                              gamma_tri(MAX_DUNAVANT_ORDER)

contains

!          --------------------
subroutine quadrature_rule_line(order,xvert,yvert,nqpoints,qweights, &
                                xquad,yquad,ierr)
!          --------------------

!----------------------------------------------------
! This routine returns quadrature points and weights for a quadrature rule
! of order order over the line in 2D defined by (xvert,yvert). The integral is
! given by the sum of the weights times the function at the quadrature points
! (the weights include the domain size).
!
! ierr is the return status:
!   0 - sucess
!   1 - no rule for quadrature of the requested order
!   2 - memory allocation failed
!
! This routine allocates qweights, xquad and yquad; they should be deallocated
! by the caller.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: order
real(RKIND), intent(in) :: xvert(2),yvert(2)
integer, intent(out) :: nqpoints
real(RKIND), pointer :: qweights(:)
real(RKIND), pointer :: xquad(:),yquad(:)
integer, intent(inout) :: ierr
!----------------------------------------------------
! Local variables:

integer :: astat, ord, i
real(RKIND) :: area, alpha
logical, save :: first_call = .true.
!----------------------------------------------------
! Begin executable code

ierr = 0

if (order < 1 .or. order > MAX_QUAD_ORDER_LINE) then
   ierr = 1
   return
endif

! on first call, set up the quadrature rule data base

if (first_call) then
   call setup_line(ierr)
   first_call = .false.
endif

! some rules are missing; use the first available rule that is larger

ord = order
if (order == 11) ord = 12
if (order > 12 .and. order < 16) ord = 16
if (order > 16 .and. order < 20) ord = 20
if (order > 20 .and. order < 24) ord = 24

! length of line

area = sqrt(abs((xvert(2)-xvert(1))**2 + (yvert(2)-yvert(1))**2))

! copy/compute quadrature rule

nqpoints = ng_line(ord)
allocate(qweights(nqpoints),xquad(nqpoints),yquad(nqpoints),stat=astat)
if (astat /= 0) then
   ierr = 2
   return
endif
do i=1,ord/2
   qweights(i) = area * weight_line(ord)%val((ord+1)/2-i+1)/2
   alpha = (1-alpha_line(ord)%val((ord+1)/2-i+1))/2
   xquad(i) = alpha*xvert(2) + (1-alpha)*xvert(1)
   yquad(i) = alpha*yvert(2) + (1-alpha)*yvert(1)
end do
do i=ord/2+1,ord
   qweights(i) = area * weight_line(ord)%val(i-ord/2)/2
   alpha = (1+alpha_line(ord)%val(i-ord/2))/2
   xquad(i) = alpha*xvert(2) + (1-alpha)*xvert(1)
   yquad(i) = alpha*yvert(2) + (1-alpha)*yvert(1)
end do

end subroutine quadrature_rule_line

!          -------------------
subroutine quadrature_rule_tri(order,xvert,yvert,nqpoints,qweights, &
                               xquad,yquad,ierr,stay_in)
!          -------------------

!----------------------------------------------------
! This routine returns quadrature points and weights for a quadrature rule
! of order order over the triangle defined by (xvert,yvert). The integral is
! given by the sum of the weights times the function at the quadrature points
! (the weights include the domain size).
!
! Some of the quadrature rules contain quadrature points that are outside
! the triangle.  If stay_in is present and true, a quadrature rule that
! uses points only in the triangle is returned.  This may be a rule of higher
! order than requested.
!
! For some orders a rule is not given.  In these cases the first higher order
! rule is returned.
!
! ierr is the return status:
!   0 - sucess
!   1 - no rule for quadrature of the requested order
!   2 - memory allocation failed
!
! This routine allocates qweights, xquad and yquad; they should be deallocated
! by the caller.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: order
real(RKIND), intent(in) :: xvert(3),yvert(3)
integer, intent(out) :: nqpoints
real(RKIND), pointer :: qweights(:)
real(RKIND), pointer :: xquad(:),yquad(:)
integer, intent(inout) :: ierr
logical, intent(in), optional :: stay_in
!----------------------------------------------------
! Local variables:

logical :: loc_stay_in
!----------------------------------------------------
! Begin executable code

ierr = 0

! check for valid order

if (order < 1 .or. order > MAX_QUAD_ORDER_TRI) then
   ierr = 1
   return
endif

! Select between Dunvant rules and product rules.  Special case for
! order=20 with stay_in true.

if (present(stay_in)) then
   loc_stay_in = stay_in
else
   loc_stay_in = .false.
endif
if (order == 20 .and. loc_stay_in) then
   call quadrature_rule_tri_product(order,xvert,yvert,nqpoints,qweights, &
                                    xquad,yquad,ierr,stay_in)
elseif (order > MAX_DUNAVANT_ORDER) then
   call quadrature_rule_tri_product(order,xvert,yvert,nqpoints,qweights, &
                                    xquad,yquad,ierr,stay_in)
else
   call quadrature_rule_tri_dunavant(order,xvert,yvert,nqpoints,qweights, &
                                     xquad,yquad,ierr,stay_in)
endif

end subroutine quadrature_rule_tri

!          ----------------------------
subroutine quadrature_rule_tri_dunavant(order,xvert,yvert,nqpoints,qweights, &
                                        xquad,yquad,ierr,stay_in)
!          ----------------------------

!----------------------------------------------------
! This routine returns quadrature points and weights as defined by
! the Dunavant rules.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: order
real(RKIND), intent(in) :: xvert(3),yvert(3)
integer, intent(out) :: nqpoints
real(RKIND), pointer :: qweights(:)
real(RKIND), pointer :: xquad(:),yquad(:)
integer, intent(inout) :: ierr
logical, intent(in), optional :: stay_in
!----------------------------------------------------
! Local variables:

integer :: astat, loc_order
real(RKIND) :: area
logical, save :: first_call = .true.
!----------------------------------------------------
! Begin executable code

! determine actual order to use

if (present(stay_in)) then
   if (stay_in) then
      select case(order)
      case (11)
         loc_order = 12
      case (15, 16)
         loc_order = 17
      case (18, 20)
         loc_order = 19
      case default
         loc_order = order
      end select
   else
      loc_order = order
   endif
else
   loc_order = order
endif

! on first call, set up the quadrature rule data base

if (first_call) then
   call setup_tri(ierr)
   first_call = .false.
endif

! area of the triangle

area = abs(xvert(1)*(yvert(2)-yvert(3)) + &
           xvert(2)*(yvert(3)-yvert(1)) + &
           xvert(3)*(yvert(1)-yvert(2))) / 2

! copy quadrature rule

nqpoints = ng_tri(loc_order)
allocate(qweights(nqpoints),xquad(nqpoints),yquad(nqpoints),stat=astat)
if (astat /= 0) then
   ierr = 2
   return
endif
qweights = weight_tri(loc_order)%val * area
xquad = xvert(1)*alpha_tri(loc_order)%val + &
        xvert(2)*beta_tri(loc_order)%val  + &
        xvert(3)*gamma_tri(loc_order)%val
yquad = yvert(1)*alpha_tri(loc_order)%val + &
        yvert(2)*beta_tri(loc_order)%val  + &
        yvert(3)*gamma_tri(loc_order)%val

end subroutine quadrature_rule_tri_dunavant

!          ---------------------------
subroutine quadrature_rule_tri_product(order,xvert,yvert,nqpoints,qweights, &
                                       xquad,yquad,ierr,stay_in)
!          ---------------------------

!----------------------------------------------------
! This routine returns quadrature points and weights for the triangle by
! using a product rule and a mapping from the square to the triangle.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: order
real(RKIND), intent(in) :: xvert(3),yvert(3)
integer, intent(out) :: nqpoints
real(RKIND), pointer :: qweights(:)
real(RKIND), pointer :: xquad(:),yquad(:)
integer, intent(inout) :: ierr
logical, intent(in), optional :: stay_in
!----------------------------------------------------
! Local variables:

integer :: nqpoints_1D, astat, i, j, p
real(RKIND), pointer :: qweights_1D(:), xquad_1D(:), yquad_1D(:)
real(RKIND) :: zeta1, zeta2, zeta3, area
!----------------------------------------------------
! Begin executable code

! get the quadrature rule for the unit line with an order that gives
! the right convergence rate

call quadrature_rule_line(2+order/2,(/0.0_RKIND,1.0_RKIND/), &
                          (/0.0_RKIND,0.0_RKIND/),nqpoints_1D,qweights_1D, &
                          xquad_1D,yquad_1D,ierr)
if (ierr /= 0) return

! allocate memory for the rules

nqpoints = nqpoints_1D**2
allocate(qweights(nqpoints),xquad(nqpoints),yquad(nqpoints),stat=astat)
if (astat /= 0) then
   ierr = 2
   return
endif

! The Duffy transformation says
!
! int_0_1 int_0_x f(x,y) dx dy = int_0_1 int_0_1 x f(x,xt) dt dx
!
! so from the product rule on the square is mapped to the reference triangle
! (0,0), (1,1), (1,0) by mapping the square's quadrature point (x,y) to the
! triangle's point (x,xy) and multiplying the weight by x.  The reference
! triangle is then mapped to the requested triangle via the barycentric
! coordinates noting that the barycentric coordinates for the reference
! triangle are given by 1-x, y and x-y.  Also multiply the weights by the
! area of the triangle.

! area of the triangle

area = abs(xvert(1)*(yvert(2)-yvert(3)) + &
           xvert(2)*(yvert(3)-yvert(1)) + &
           xvert(3)*(yvert(1)-yvert(2)))

p = 0
do i=1,nqpoints_1D
   do j=1,nqpoints_1D
      p = p+1
      xquad(p) = xquad_1D(i)
      yquad(p) = xquad_1d(j)*xquad_1D(i)
      qweights(p) = qweights_1D(i)*qweights_1D(j)*xquad_1D(i)*area
      zeta1 = 1-xquad(p)
      zeta2 = yquad(p)
      zeta3 = xquad(p) - yquad(p)
      xquad(p) = zeta1*xvert(1) + zeta2*xvert(2) + zeta3*xvert(3)
      yquad(p) = zeta1*yvert(1) + zeta2*yvert(2) + zeta3*yvert(3)
   end do
end do

! free up memory from the 1D rules

deallocate(qweights_1D,xquad_1D,yquad_1D,stat=astat)

end subroutine quadrature_rule_tri_product

!          ----------
subroutine setup_line(ierr)
!          ----------

!----------------------------------------------------
! This routine is called the first time a quadrature rule over a line
! is requested, and defines all the weights and quadrature points.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(inout) :: ierr
!----------------------------------------------------
! Local variables:

integer :: i, astat, dim
!----------------------------------------------------
! Begin executable code

! define the number of quadrature points

ng_line = (/(i,i=1,MAX_QUAD_ORDER_LINE)/)

! allocate space for the quadrature rules

do i=1,MAX_QUAD_ORDER_LINE
   dim = i/2
   if (2*(i/2)/=i) dim = dim+1
   allocate(weight_line(i)%val(dim), alpha_line(i)%val(dim),stat=astat)
   if (astat /= 0) then
      ierr = 2
      return
   endif
end do

! For each rule, define the weights and points

weight_line(1)%val = (/ &
 2.000000000000000D+00 /)
alpha_line(1)%val = (/ &
 0.000000000000000D+00 /)

weight_line(2)%val = (/ &
 1.000000000000000D+00 /)
alpha_line(2)%val = (/ &
 0.577350269189626D+00 /)

weight_line(3)%val = (/ &
 0.888888888888889D+00, 0.555555555555556D+00 /)
alpha_line(3)%val = (/ &
 0.000000000000000D+00, 0.774596669241483D+00 /)

weight_line(4)%val = (/ &
 0.652145154862546D+00, 0.347854845137454D+00 /)
alpha_line(4)%val = (/ &
 0.339981043584856D+00, 0.861136311594053D+00 /)

weight_line(5)%val = (/ &
 0.568888888888889D+00, 0.478628670499366D+00, 0.236926885056189D+00 /)
alpha_line(5)%val = (/ &
 0.000000000000000D+00, 0.538469310105683D+00, 0.906179845938664D+00 /)

weight_line(6)%val = (/ &
 0.467913934572691D+00, 0.360761573048139D+00, 0.171324492379170D+00 /)
alpha_line(6)%val = (/ &
 0.238619186083197D+00, 0.661209386466265D+00, 0.932469514203152D+00 /)

weight_line(7)%val = (/ &
 0.417959183673469D+00, 0.381830050505119D+00, 0.279705391489277D+00, &
 0.129484966168870D+00 /)
alpha_line(7)%val = (/ &
 0.000000000000000D+00, 0.405845151377397D+00, 0.741531185599394D+00, &
 0.949107912342759D+00 /)

weight_line(8)%val = (/ &
 0.362683783378362D+00, 0.313706645877887D+00, 0.222381034453374D+00, &
 0.101228536290376D+00 /)
alpha_line(8)%val = (/ &
 0.183434642495650D+00, 0.525532409916329D+00, 0.796666477413627D+00, &
 0.960289856497536D+00 /)

weight_line(9)%val = (/ &
 0.330239355001260D+00, 0.312347077040003D+00, 0.260610696402935D+00, &
 0.180648160694857D+00, 0.081274388361574D+00 /)
alpha_line(9)%val = (/ &
 0.000000000000000D+00, 0.324253423403809D+00, 0.613371432700590D+00, &
 0.836031107326636D+00, 0.968160239507626D+00 /)

weight_line(10)%val = (/ &
 0.295524224714753D+00, 0.269266719309996D+00, 0.219086362515982D+00, &
 0.149451349150581D+00, 0.066671344308688D+00 /)
alpha_line(10)%val = (/ &
 0.148874338981631D+00, 0.433395394129247D+00, 0.679409568299024D+00, &
 0.865063366688985D+00, 0.973906528517172D+00 /)

weight_line(12)%val = (/ &
 0.249147045813403D+00, 0.233492536538355D+00, 0.203167426723066D+00, &
 0.160078328543346D+00, 0.106939325995318D+00, 0.047175336386512D+00 /)
alpha_line(12)%val = (/ &
 0.125233408511469D+00, 0.367831498998180D+00, 0.587317954286617D+00, &
 0.769902674194305D+00, 0.904117256370475D+00, 0.981560634246719D+00 /)

weight_line(16)%val = (/ &
 0.189450610455068496285D+00, 0.182603415044923588867D+00, &
 0.169156519395002538189D+00, 0.149595988816576732081D+00, &
 0.124628971255533872052D+00, 0.095158511682492784810D+00, &
 0.062253523938647892863D+00, 0.027152459411754094852D+00 /)
alpha_line(16)%val = (/ &
 0.095012509837637440185D+00, 0.281603550779258913230D+00, &
 0.458016777657227386342D+00, 0.617876244402643748447D+00, &
 0.755404408355003033895D+00, 0.865631202387831743880D+00, &
 0.944575023073232576078D+00, 0.989400934991649932596D+00 /)

weight_line(20)%val = (/ &
 0.152753387130725850698D+00, 0.149172986472603746788D+00, &
 0.142096109318382051329D+00, 0.131688638449176626898D+00, &
 0.118194531961518417312D+00, 0.101930119817240435037D+00, &
 0.083276741576704748725D+00, 0.062672048334109063570D+00, &
 0.040601429800386941331D+00, 0.017614007139152118312D+00 /)
alpha_line(20)%val = (/ &
 0.076526521133497333755D+00, 0.227785851141645078080D+00, &
 0.373706088715419560673D+00, 0.510867001950827098004D+00, &
 0.636053680726515025453D+00, 0.746331906460150792614D+00, &
 0.839116971822218823395D+00, 0.912234428251325905868D+00, &
 0.963971927277913791268D+00, 0.993128599185094924786D+00 /)

weight_line(24)%val = (/ &
 0.127938195346752156974D+00, 0.125837456346828296121D+00, &
 0.121670472927803391204D+00, 0.115505668053725601353D+00, &
 0.107444270115965634783D+00, 0.097618652104113888270D+00, &
 0.086190161531953275917D+00, 0.073346481411080305734D+00, &
 0.059298584915436780746D+00, 0.044277438817419806169D+00, &
 0.028531388628933663181D+00, 0.012341229799987199547D+00 /)
alpha_line(24)%val = (/ &
 0.064056892862605626085D+00, 0.191118867473616309159D+00, &
 0.315042679696163374387D+00, 0.433793507626045138487D+00, &
 0.545421471388839535658D+00, 0.648093651936975569252D+00, &
 0.740124191578554364244D+00, 0.820001985973902921954D+00, &
 0.886415527004401034213D+00, 0.938274552002732758524D+00, &
 0.974728555971309498198D+00, 0.995187219997021360180D+00 /)
end subroutine setup_line

!          ---------
subroutine setup_tri(ierr)
!          ---------

!----------------------------------------------------
! This routine is called the first time a quadrature rule over a triangle
! is requested, and defines all the weights and quadrature points.
! These quadrature rules come from
! Dunavant, D.A., "High Degree Efficient Symmetrical Gaussian Quadrature
!   Rules for the Triangle", Int. J. Num. Meth. Eng., 21 (1985) pp 1129-1148.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(inout) :: ierr
!----------------------------------------------------
! Local variables:

integer :: i, astat
!----------------------------------------------------
! Begin executable code

! define the number of quadrature points

ng_tri = (/1,3,4,6,7,12,13,16,19,25,27,33,37,42,48,52,61,70,73,79/)

! allocate space for the quadrature rules

do i=1,MAX_DUNAVANT_ORDER
   allocate(weight_tri(i)%val(ng_tri(i)), alpha_tri(i)%val(ng_tri(i)), &
            beta_tri(i)%val(ng_tri(i)), gamma_tri(i)%val(ng_tri(i)), stat=astat)
   if (astat /= 0) then
      ierr = 2
      return
   endif
end do

! For each rule, define the weights and points

weight_tri(1)%val = (/ &
 1.000000000000000D+00 /)
alpha_tri(1)%val = (/ &
 0.333333333333333D+00 /)
beta_tri(1)%val = (/ &
 0.333333333333333D+00 /)
gamma_tri(1)%val = (/ &
 0.333333333333333D+00 /)

weight_tri(2)%val = (/ &
 0.333333333333333D+00, 0.333333333333333D+00, 0.333333333333333D+00 /)
alpha_tri(2)%val = (/ &
 0.666666666666667D+00, 0.166666666666667D+00, 0.166666666666667D+00 /)
beta_tri(2)%val = (/ &
 0.166666666666667D+00, 0.666666666666667D+00, 0.166666666666667D+00 /)
gamma_tri(2)%val = (/ &
 0.166666666666667D+00, 0.166666666666667D+00, 0.666666666666667D+00 /)

weight_tri(3)%val = (/ &
-0.562500000000000D+00, 0.520833333333333D+00, 0.520833333333333D+00, &
 0.520833333333333D+00 /)
alpha_tri(3)%val = (/ &
 0.333333333333333D+00, 0.600000000000000D+00, 0.200000000000000D+00, &
 0.200000000000000D+00 /)
beta_tri(3)%val = (/ &
 0.333333333333333D+00, 0.200000000000000D+00, 0.600000000000000D+00, &
 0.200000000000000D+00 /)
gamma_tri(3)%val = (/ &
 0.333333333333333D+00, 0.200000000000000D+00, 0.200000000000000D+00, &
 0.600000000000000D+00 /)

weight_tri(4)%val = (/ &
 0.223381589678011D+00, 0.223381589678011D+00, 0.223381589678011D+00, &
 0.109951743655322D+00, 0.109951743655322D+00, 0.109951743655322D+00 /)
alpha_tri(4)%val = (/ &
 0.108103018168070D+00, 0.445948490915965D+00, 0.445948490915965D+00, &
 0.816847572980459D+00, 0.091576213509771D+00, 0.091576213509771D+00 /)
beta_tri(4)%val = (/ &
 0.445948490915965D+00, 0.108103018168070D+00, 0.445948490915965D+00, &
 0.091576213509771D+00, 0.816847572980459D+00, 0.091576213509771D+00 /)
gamma_tri(4)%val = (/ &
 0.445948490915965D+00, 0.445948490915965D+00, 0.108103018168070D+00, &
 0.091576213509771D+00, 0.091576213509771D+00, 0.816847572980459D+00 /)

weight_tri(5)%val = (/ &
 0.225000000000000D+00, 0.132394152788506D+00, 0.132394152788506D+00, &
 0.132394152788506D+00, 0.125939180544827D+00, 0.125939180544827D+00, &
 0.125939180544827D+00 /)
alpha_tri(5)%val = (/ &
 0.333333333333333D+00, 0.059715871789770D+00, 0.470142064105115D+00, &
 0.470142064105115D+00, 0.797426985353087D+00, 0.101286507323456D+00, &
 0.101286507323456D+00 /)
beta_tri(5)%val = (/ &
 0.333333333333333D+00, 0.470142064105115D+00, 0.059715871789770D+00, &
 0.470142064105115D+00, 0.101286507323456D+00, 0.797426985353087D+00, &
 0.101286507323456D+00 /)
gamma_tri(5)%val = (/ &
 0.333333333333333D+00, 0.470142064105115D+00, 0.470142064105115D+00, &
 0.059715871789770D+00, 0.101286507323456D+00, 0.101286507323456D+00, &
 0.797426985353087D+00 /)

weight_tri(6)%val = (/ &
 0.116786275726379D+00, 0.116786275726379D+00, 0.116786275726379D+00, &
 0.050844906370207D+00, 0.050844906370207D+00, 0.050844906370207D+00, &
 0.082851075618374D+00, 0.082851075618374D+00, 0.082851075618374D+00, &
 0.082851075618374D+00, 0.082851075618374D+00, 0.082851075618374D+00 /)
alpha_tri(6)%val = (/ &
 0.501426509658179D+00, 0.249286745170910D+00, 0.249286745170910D+00, &
 0.873821971016996D+00, 0.063089014491502D+00, 0.063089014491502D+00, &
 0.053145049844817D+00, 0.053145049844817D+00, 0.310352451033784D+00, &
 0.310352451033784D+00, 0.636502499121399D+00, 0.636502499121399D+00 /)
beta_tri(6)%val = (/ &
 0.249286745170910D+00, 0.501426509658179D+00, 0.249286745170910D+00, &
 0.063089014491502D+00, 0.873821971016996D+00, 0.063089014491502D+00, &
 0.310352451033784D+00, 0.636502499121399D+00, 0.053145049844817D+00, &
 0.636502499121399D+00, 0.310352451033784D+00, 0.053145049844817D+00 /)
gamma_tri(6)%val = (/ &
 0.249286745170910D+00, 0.249286745170910D+00, 0.501426509658179D+00, &
 0.063089014491502D+00, 0.063089014491502D+00, 0.873821971016996D+00, &
 0.636502499121399D+00, 0.310352451033784D+00, 0.636502499121399D+00, &
 0.053145049844817D+00, 0.053145049844817D+00, 0.310352451033784D+00 /)

weight_tri(7)%val = (/ &
-0.149570044467682D+00, 0.175615257433208D+00, 0.175615257433208D+00, &
 0.175615257433208D+00, 0.053347235608838D+00, 0.053347235608838D+00, &
 0.053347235608838D+00, 0.077113760890257D+00, 0.077113760890257D+00, &
 0.077113760890257D+00, 0.077113760890257D+00, 0.077113760890257D+00, &
 0.077113760890257D+00 /)
alpha_tri(7)%val = (/ &
 0.333333333333333D+00, 0.479308067841920D+00, 0.260345966079040D+00, &
 0.260345966079040D+00, 0.869739794195568D+00, 0.065130102902216D+00, &
 0.065130102902216D+00, 0.048690315425316D+00, 0.048690315425316D+00, &
 0.312865496004874D+00, 0.312865496004874D+00, 0.638444188569810D+00, &
 0.638444188569810D+00 /)
beta_tri(7)%val = (/ &
 0.333333333333333D+00, 0.260345966079040D+00, 0.479308067841920D+00, &
 0.260345966079040D+00, 0.065130102902216D+00, 0.869739794195568D+00, &
 0.065130102902216D+00, 0.312865496004874D+00, 0.638444188569810D+00, &
 0.048690315425316D+00, 0.638444188569810D+00, 0.312865496004874D+00, &
 0.048690315425316D+00 /)
gamma_tri(7)%val = (/ &
 0.333333333333333D+00, 0.260345966079040D+00, 0.260345966079040D+00, &
 0.479308067841920D+00, 0.065130102902216D+00, 0.065130102902216D+00, &
 0.869739794195568D+00, 0.638444188569810D+00, 0.312865496004874D+00, &
 0.638444188569810D+00, 0.048690315425316D+00, 0.048690315425316D+00, &
 0.312865496004874D+00 /)

weight_tri(8)%val = (/ &
 0.144315607677787D+00, 0.095091634267285D+00, 0.095091634267285D+00, &
 0.095091634267285D+00, 0.103217370534718D+00, 0.103217370534718D+00, &
 0.103217370534718D+00, 0.032458497623198D+00, 0.032458497623198D+00, &
 0.032458497623198D+00, 0.027230314174435D+00, 0.027230314174435D+00, &
 0.027230314174435D+00, 0.027230314174435D+00, 0.027230314174435D+00, &
 0.027230314174435D+00 /)
alpha_tri(8)%val = (/ &
 0.333333333333333D+00, 0.081414823414554D+00, 0.459292588292723D+00, &
 0.459292588292723D+00, 0.658861384496480D+00, 0.170569307751760D+00, &
 0.170569307751760D+00, 0.898905543365938D+00, 0.050547228317031D+00, &
 0.050547228317031D+00, 0.008394777409958D+00, 0.008394777409958D+00, &
 0.263112829634638D+00, 0.263112829634638D+00, 0.728492392955404D+00, &
 0.728492392955404D+00 /)
beta_tri(8)%val = (/ &
 0.333333333333333D+00, 0.459292588292723D+00, 0.081414823414554D+00, &
 0.459292588292723D+00, 0.170569307751760D+00, 0.658861384496480D+00, &
 0.170569307751760D+00, 0.050547228317031D+00, 0.898905543365938D+00, &
 0.050547228317031D+00, 0.263112829634638D+00, 0.728492392955404D+00, &
 0.008394777409958D+00, 0.728492392955404D+00, 0.263112829634638D+00, &
 0.008394777409958D+00 /)
gamma_tri(8)%val = (/ &
 0.333333333333333D+00, 0.459292588292723D+00, 0.459292588292723D+00, &
 0.081414823414554D+00, 0.170569307751760D+00, 0.170569307751760D+00, &
 0.658861384496480D+00, 0.050547228317031D+00, 0.050547228317031D+00, &
 0.898905543365938D+00, 0.728492392955404D+00, 0.263112829634638D+00, &
 0.728492392955404D+00, 0.008394777409958D+00, 0.008394777409958D+00, &
 0.263112829634638D+00 /)

weight_tri(9)%val = (/ &
 0.097135796282799D+00, 0.031334700227139D+00, 0.031334700227139D+00, &
 0.031334700227139D+00, 0.077827541004774D+00, 0.077827541004774D+00, &
 0.077827541004774D+00, 0.079647738927210D+00, 0.079647738927210D+00, &
 0.079647738927210D+00, 0.025577675658698D+00, 0.025577675658698D+00, &
 0.025577675658698D+00, 0.043283539377289D+00, 0.043283539377289D+00, &
 0.043283539377289D+00, 0.043283539377289D+00, 0.043283539377289D+00, &
 0.043283539377289D+00 /)
alpha_tri(9)%val = (/ &
 0.333333333333333D+00, 0.020634961602525D+00, 0.489682519198738D+00, &
 0.489682519198738D+00, 0.125820817014127D+00, 0.437089591492937D+00, &
 0.437089591492937D+00, 0.623592928761935D+00, 0.188203535619033D+00, &
 0.188203535619033D+00, 0.910540973211095D+00, 0.044729513394453D+00, &
 0.044729513394453D+00, 0.036838412054736D+00, 0.036838412054736D+00, &
 0.221962989160766D+00, 0.221962989160766D+00, 0.741198598784498D+00, &
 0.741198598784498D+00 /)
beta_tri(9)%val = (/ &
 0.333333333333333D+00, 0.489682519198738D+00, 0.020634961602525D+00, &
 0.489682519198738D+00, 0.437089591492937D+00, 0.125820817014127D+00, &
 0.437089591492937D+00, 0.188203535619033D+00, 0.623592928761935D+00, &
 0.188203535619033D+00, 0.044729513394453D+00, 0.910540973211095D+00, &
 0.044729513394453D+00, 0.221962989160766D+00, 0.741198598784498D+00, &
 0.036838412054736D+00, 0.741198598784498D+00, 0.221962989160766D+00, &
 0.036838412054736D+00 /)
gamma_tri(9)%val = (/ &
 0.333333333333333D+00, 0.489682519198738D+00, 0.489682519198738D+00, &
 0.020634961602525D+00, 0.437089591492937D+00, 0.437089591492937D+00, &
 0.125820817014127D+00, 0.188203535619033D+00, 0.188203535619033D+00, &
 0.623592928761935D+00, 0.044729513394453D+00, 0.044729513394453D+00, &
 0.910540973211095D+00, 0.741198598784498D+00, 0.221962989160766D+00, &
 0.741198598784498D+00, 0.036838412054736D+00, 0.036838412054736D+00, &
 0.221962989160766D+00 /)

weight_tri(10)%val = (/ &
 0.908179903827540D-01, 0.367259577564670D-01, 0.367259577564670D-01, &
 0.367259577564670D-01, 0.453210594355280D-01, 0.453210594355280D-01, &
 0.453210594355280D-01, 0.727579168454200D-01, 0.727579168454200D-01, &
 0.727579168454200D-01, 0.727579168454200D-01, 0.727579168454200D-01, &
 0.727579168454200D-01, 0.283272425310570D-01, 0.283272425310570D-01, &
 0.283272425310570D-01, 0.283272425310570D-01, 0.283272425310570D-01, &
 0.283272425310570D-01, 0.942166696373300D-02, 0.942166696373300D-02, &
 0.942166696373300D-02, 0.942166696373300D-02, 0.942166696373300D-02, &
 0.942166696373300D-02 /)
alpha_tri(10)%val = (/ &
 0.333333333333333D+00, 0.288447332326850D-01, 0.485577633383657D+00, &
 0.485577633383657D+00, 0.781036849029926D+00, 0.109481575485037D+00, &
 0.109481575485037D+00, 0.141707219414880D+00, 0.141707219414880D+00, &
 0.307939838764121D+00, 0.307939838764121D+00, 0.550352941820999D+00, &
 0.550352941820999D+00, 0.250035347626860D-01, 0.250035347626860D-01, &
 0.246672560639903D+00, 0.246672560639903D+00, 0.728323904597411D+00, &
 0.728323904597411D+00, 0.954081540029900D-02, 0.954081540029900D-02, &
 0.668032510122000D-01, 0.668032510122000D-01, 0.923655933587500D+00, &
 0.923655933587500D+00 /)
beta_tri(10)%val = (/ &
 0.333333333333333D+00, 0.485577633383657D+00, 0.288447332326850D-01, &
 0.485577633383657D+00, 0.109481575485037D+00, 0.781036849029926D+00, &
 0.109481575485037D+00, 0.307939838764121D+00, 0.550352941820999D+00, &
 0.141707219414880D+00, 0.550352941820999D+00, 0.307939838764121D+00, &
 0.141707219414880D+00, 0.246672560639903D+00, 0.728323904597411D+00, &
 0.250035347626860D-01, 0.728323904597411D+00, 0.246672560639903D+00, &
 0.250035347626860D-01, 0.668032510122000D-01, 0.923655933587500D+00, &
 0.954081540029900D-02, 0.923655933587500D+00, 0.668032510122000D-01, &
 0.954081540029900D-02 /)
gamma_tri(10)%val = (/ &
 0.333333333333333D+00, 0.485577633383657D+00, 0.485577633383657D+00, &
 0.288447332326850D-01, 0.109481575485037D+00, 0.109481575485037D+00, &
 0.781036849029926D+00, 0.550352941820999D+00, 0.307939838764121D+00, &
 0.550352941820999D+00, 0.141707219414880D+00, 0.141707219414880D+00, &
 0.307939838764121D+00, 0.728323904597411D+00, 0.246672560639903D+00, &
 0.728323904597411D+00, 0.250035347626860D-01, 0.250035347626860D-01, &
 0.246672560639903D+00, 0.923655933587500D+00, 0.668032510122000D-01, &
 0.923655933587500D+00, 0.954081540029900D-02, 0.954081540029900D-02, &
 0.668032510122000D-01 /)

weight_tri(11)%val = (/ &
 0.927006328961000D-03, 0.927006328961000D-03, 0.927006328961000D-03, &
 0.771495349148130D-01, 0.771495349148130D-01, 0.771495349148130D-01, &
 0.593229773807740D-01, 0.593229773807740D-01, 0.593229773807740D-01, &
 0.361845405034180D-01, 0.361845405034180D-01, 0.361845405034180D-01, &
 0.136597310026780D-01, 0.136597310026780D-01, 0.136597310026780D-01, &
 0.523371119622040D-01, 0.523371119622040D-01, 0.523371119622040D-01, &
 0.523371119622040D-01, 0.523371119622040D-01, 0.523371119622040D-01, &
 0.207076596391410D-01, 0.207076596391410D-01, 0.207076596391410D-01, &
 0.207076596391410D-01, 0.207076596391410D-01, 0.207076596391410D-01 /)
alpha_tri(11)%val = (/ &
 -.692220965415170D-01, 0.534611048270758D+00, 0.534611048270758D+00, &
 0.202061394068290D+00, 0.398969302965855D+00, 0.398969302965855D+00, &
 0.593380199137435D+00, 0.203309900431282D+00, 0.203309900431282D+00, &
 0.761298175434837D+00, 0.119350912282581D+00, 0.119350912282581D+00, &
 0.935270103777448D+00, 0.323649481112760D-01, 0.323649481112760D-01, &
 0.501781383104950D-01, 0.501781383104950D-01, 0.356620648261293D+00, &
 0.356620648261293D+00, 0.593201213428213D+00, 0.593201213428213D+00, &
 0.210220165361660D-01, 0.210220165361660D-01, 0.171488980304042D+00, &
 0.171488980304042D+00, 0.807489003159792D+00, 0.807489003159792D+00 /)
beta_tri(11)%val = (/ &
 0.534611048270758D+00, -.692220965415170D-01, 0.534611048270758D+00, &
 0.398969302965855D+00, 0.202061394068290D+00, 0.398969302965855D+00, &
 0.203309900431282D+00, 0.593380199137435D+00, 0.203309900431282D+00, &
 0.119350912282581D+00, 0.761298175434837D+00, 0.119350912282581D+00, &
 0.323649481112760D-01, 0.935270103777448D+00, 0.323649481112760D-01, &
 0.356620648261293D+00, 0.593201213428213D+00, 0.501781383104950D-01, &
 0.593201213428213D+00, 0.356620648261293D+00, 0.501781383104950D-01, &
 0.171488980304042D+00, 0.807489003159792D+00, 0.210220165361660D-01, &
 0.807489003159792D+00, 0.171488980304042D+00, 0.210220165361660D-01 /)
gamma_tri(11)%val = (/ &
 0.534611048270758D+00, 0.534611048270758D+00, -.692220965415170D-01, &
 0.398969302965855D+00, 0.398969302965855D+00, 0.202061394068290D+00, &
 0.203309900431282D+00, 0.203309900431282D+00, 0.593380199137435D+00, &
 0.119350912282581D+00, 0.119350912282581D+00, 0.761298175434837D+00, &
 0.323649481112760D-01, 0.323649481112760D-01, 0.935270103777448D+00, &
 0.593201213428213D+00, 0.356620648261293D+00, 0.593201213428213D+00, &
 0.501781383104950D-01, 0.501781383104950D-01, 0.356620648261293D+00, &
 0.807489003159792D+00, 0.171488980304042D+00, 0.807489003159792D+00, &
 0.210220165361660D-01, 0.210220165361660D-01, 0.171488980304042D+00 /)

weight_tri(12)%val = (/ &
 0.257310664404550D-01, 0.257310664404550D-01, 0.257310664404550D-01, &
 0.436925445380380D-01, 0.436925445380380D-01, 0.436925445380380D-01, &
 0.628582242178850D-01, 0.628582242178850D-01, 0.628582242178850D-01, &
 0.347961129307090D-01, 0.347961129307090D-01, 0.347961129307090D-01, &
 0.616626105155900D-02, 0.616626105155900D-02, 0.616626105155900D-02, &
 0.403715577663810D-01, 0.403715577663810D-01, 0.403715577663810D-01, &
 0.403715577663810D-01, 0.403715577663810D-01, 0.403715577663810D-01, &
 0.223567732023030D-01, 0.223567732023030D-01, 0.223567732023030D-01, &
 0.223567732023030D-01, 0.223567732023030D-01, 0.223567732023030D-01, &
 0.173162311086590D-01, 0.173162311086590D-01, 0.173162311086590D-01, &
 0.173162311086590D-01, 0.173162311086590D-01, 0.173162311086590D-01 /)
alpha_tri(12)%val = (/ &
 0.235652204523900D-01, 0.488217389773805D+00, 0.488217389773805D+00, &
 0.120551215411079D+00, 0.439724392294460D+00, 0.439724392294460D+00, &
 0.457579229975768D+00, 0.271210385012116D+00, 0.271210385012116D+00, &
 0.744847708916828D+00, 0.127576145541586D+00, 0.127576145541586D+00, &
 0.957365299093579D+00, 0.213173504532100D-01, 0.213173504532100D-01, &
 0.115343494534698D+00, 0.115343494534698D+00, 0.275713269685514D+00, &
 0.275713269685514D+00, 0.608943235779788D+00, 0.608943235779788D+00, &
 0.228383322222570D-01, 0.228383322222570D-01, 0.281325580989940D+00, &
 0.281325580989940D+00, 0.695836086787803D+00, 0.695836086787803D+00, &
 0.257340505483300D-01, 0.257340505483300D-01, 0.116251915907597D+00, &
 0.116251915907597D+00, 0.858014033544073D+00, 0.858014033544073D+00 /)
beta_tri(12)%val = (/ &
 0.488217389773805D+00, 0.235652204523900D-01, 0.488217389773805D+00, &
 0.439724392294460D+00, 0.120551215411079D+00, 0.439724392294460D+00, &
 0.271210385012116D+00, 0.457579229975768D+00, 0.271210385012116D+00, &
 0.127576145541586D+00, 0.744847708916828D+00, 0.127576145541586D+00, &
 0.213173504532100D-01, 0.957365299093579D+00, 0.213173504532100D-01, &
 0.275713269685514D+00, 0.608943235779788D+00, 0.115343494534698D+00, &
 0.608943235779788D+00, 0.275713269685514D+00, 0.115343494534698D+00, &
 0.281325580989940D+00, 0.695836086787803D+00, 0.228383322222570D-01, &
 0.695836086787803D+00, 0.281325580989940D+00, 0.228383322222570D-01, &
 0.116251915907597D+00, 0.858014033544073D+00, 0.257340505483300D-01, &
 0.858014033544073D+00, 0.116251915907597D+00, 0.257340505483300D-01 /)
gamma_tri(12)%val = (/ &
 0.488217389773805D+00, 0.488217389773805D+00, 0.235652204523900D-01, &
 0.439724392294460D+00, 0.439724392294460D+00, 0.120551215411079D+00, &
 0.271210385012116D+00, 0.271210385012116D+00, 0.457579229975768D+00, &
 0.127576145541586D+00, 0.127576145541586D+00, 0.744847708916828D+00, &
 0.213173504532100D-01, 0.213173504532100D-01, 0.957365299093579D+00, &
 0.608943235779788D+00, 0.275713269685514D+00, 0.608943235779788D+00, &
 0.115343494534698D+00, 0.115343494534698D+00, 0.275713269685514D+00, &
 0.695836086787803D+00, 0.281325580989940D+00, 0.695836086787803D+00, &
 0.228383322222570D-01, 0.228383322222570D-01, 0.281325580989940D+00, &
 0.858014033544073D+00, 0.116251915907597D+00, 0.858014033544073D+00, &
 0.257340505483300D-01, 0.257340505483300D-01, 0.116251915907597D+00 /)

weight_tri(13)%val = (/ &
 0.525209234008020D-01, 0.112801452093300D-01, 0.112801452093300D-01, &
 0.112801452093300D-01, 0.314235183624540D-01, 0.314235183624540D-01, &
 0.314235183624540D-01, 0.470725025041940D-01, 0.470725025041940D-01, &
 0.470725025041940D-01, 0.473635865363550D-01, 0.473635865363550D-01, &
 0.473635865363550D-01, 0.311675290457940D-01, 0.311675290457940D-01, &
 0.311675290457940D-01, 0.797577146507400D-02, 0.797577146507400D-02, &
 0.797577146507400D-02, 0.368484027287320D-01, 0.368484027287320D-01, &
 0.368484027287320D-01, 0.368484027287320D-01, 0.368484027287320D-01, &
 0.368484027287320D-01, 0.174014633038220D-01, 0.174014633038220D-01, &
 0.174014633038220D-01, 0.174014633038220D-01, 0.174014633038220D-01, &
 0.174014633038220D-01, 0.155217868390450D-01, 0.155217868390450D-01, &
 0.155217868390450D-01, 0.155217868390450D-01, 0.155217868390450D-01, &
 0.155217868390450D-01 /)
alpha_tri(13)%val = (/ &
 0.333333333333333D+00, 0.990363012059100D-02, 0.495048184939705D+00, &
 0.495048184939705D+00, 0.625667297808520D-01, 0.468716635109574D+00, &
 0.468716635109574D+00, 0.170957326397447D+00, 0.414521336801277D+00, &
 0.414521336801277D+00, 0.541200855914337D+00, 0.229399572042831D+00, &
 0.229399572042831D+00, 0.771151009607340D+00, 0.114424495196330D+00, &
 0.114424495196330D+00, 0.950377217273082D+00, 0.248113913634590D-01, &
 0.248113913634590D-01, 0.948538283795790D-01, 0.948538283795790D-01, &
 0.268794997058761D+00, 0.268794997058761D+00, 0.636351174561660D+00, &
 0.636351174561660D+00, 0.181007732788070D-01, 0.181007732788070D-01, &
 0.291730066734288D+00, 0.291730066734288D+00, 0.690169159986905D+00, &
 0.690169159986905D+00, 0.222330766740900D-01, 0.222330766740900D-01, &
 0.126357385491669D+00, 0.126357385491669D+00, 0.851409537834241D+00, &
 0.851409537834241D+00 /)
beta_tri(13)%val = (/ &
 0.333333333333333D+00, 0.495048184939705D+00, 0.990363012059100D-02, &
 0.495048184939705D+00, 0.468716635109574D+00, 0.625667297808520D-01, &
 0.468716635109574D+00, 0.414521336801277D+00, 0.170957326397447D+00, &
 0.414521336801277D+00, 0.229399572042831D+00, 0.541200855914337D+00, &
 0.229399572042831D+00, 0.114424495196330D+00, 0.771151009607340D+00, &
 0.114424495196330D+00, 0.248113913634590D-01, 0.950377217273082D+00, &
 0.248113913634590D-01, 0.268794997058761D+00, 0.636351174561660D+00, &
 0.948538283795790D-01, 0.636351174561660D+00, 0.268794997058761D+00, &
 0.948538283795790D-01, 0.291730066734288D+00, 0.690169159986905D+00, &
 0.181007732788070D-01, 0.690169159986905D+00, 0.291730066734288D+00, &
 0.181007732788070D-01, 0.126357385491669D+00, 0.851409537834241D+00, &
 0.222330766740900D-01, 0.851409537834241D+00, 0.126357385491669D+00, &
 0.222330766740900D-01 /)
gamma_tri(13)%val = (/ &
 0.333333333333333D+00, 0.495048184939705D+00, 0.495048184939705D+00, &
 0.990363012059100D-02, 0.468716635109574D+00, 0.468716635109574D+00, &
 0.625667297808520D-01, 0.414521336801277D+00, 0.414521336801277D+00, &
 0.170957326397447D+00, 0.229399572042831D+00, 0.229399572042831D+00, &
 0.541200855914337D+00, 0.114424495196330D+00, 0.114424495196330D+00, &
 0.771151009607340D+00, 0.248113913634590D-01, 0.248113913634590D-01, &
 0.950377217273082D+00, 0.636351174561660D+00, 0.268794997058761D+00, &
 0.636351174561660D+00, 0.948538283795790D-01, 0.948538283795790D-01, &
 0.268794997058761D+00, 0.690169159986905D+00, 0.291730066734288D+00, &
 0.690169159986905D+00, 0.181007732788070D-01, 0.181007732788070D-01, &
 0.291730066734288D+00, 0.851409537834241D+00, 0.126357385491669D+00, &
 0.851409537834241D+00, 0.222330766740900D-01, 0.222330766740900D-01, &
 0.126357385491669D+00 /)

weight_tri(14)%val = (/ &
 0.218835813694290D-01, 0.218835813694290D-01, 0.218835813694290D-01, &
 0.327883535441250D-01, 0.327883535441250D-01, 0.327883535441250D-01, &
 0.517741045072920D-01, 0.517741045072920D-01, 0.517741045072920D-01, &
 0.421625887369930D-01, 0.421625887369930D-01, 0.421625887369930D-01, &
 0.144336996697770D-01, 0.144336996697770D-01, 0.144336996697770D-01, &
 0.492340360240000D-02, 0.492340360240000D-02, 0.492340360240000D-02, &
 0.246657532125640D-01, 0.246657532125640D-01, 0.246657532125640D-01, &
 0.246657532125640D-01, 0.246657532125640D-01, 0.246657532125640D-01, &
 0.385715107870610D-01, 0.385715107870610D-01, 0.385715107870610D-01, &
 0.385715107870610D-01, 0.385715107870610D-01, 0.385715107870610D-01, &
 0.144363081135340D-01, 0.144363081135340D-01, 0.144363081135340D-01, &
 0.144363081135340D-01, 0.144363081135340D-01, 0.144363081135340D-01, &
 0.501022883850100D-02, 0.501022883850100D-02, 0.501022883850100D-02, &
 0.501022883850100D-02, 0.501022883850100D-02, 0.501022883850100D-02 /)
alpha_tri(14)%val = (/ &
 0.220721792756430D-01, 0.488963910362179D+00, 0.488963910362179D+00, &
 0.164710561319092D+00, 0.417644719340454D+00, 0.417644719340454D+00, &
 0.453044943382323D+00, 0.273477528308839D+00, 0.273477528308839D+00, &
 0.645588935174913D+00, 0.177205532412543D+00, 0.177205532412543D+00, &
 0.876400233818255D+00, 0.617998830908730D-01, 0.617998830908730D-01, &
 0.961218077502598D+00, 0.193909612487010D-01, 0.193909612487010D-01, &
 0.571247574036480D-01, 0.571247574036480D-01, 0.172266687821356D+00, &
 0.172266687821356D+00, 0.770608554774996D+00, 0.770608554774996D+00, &
 0.929162493569720D-01, 0.929162493569720D-01, 0.336861459796345D+00, &
 0.336861459796345D+00, 0.570222290846683D+00, 0.570222290846683D+00, &
 0.146469500556540D-01, 0.146469500556540D-01, 0.298372882136258D+00, &
 0.298372882136258D+00, 0.686980167808088D+00, 0.686980167808088D+00, &
 0.126833093287200D-02, 0.126833093287200D-02, 0.118974497696957D+00, &
 0.118974497696957D+00, 0.879757171370171D+00, 0.879757171370171D+00 /)
beta_tri(14)%val = (/ &
 0.488963910362179D+00, 0.220721792756430D-01, 0.488963910362179D+00, &
 0.417644719340454D+00, 0.164710561319092D+00, 0.417644719340454D+00, &
 0.273477528308839D+00, 0.453044943382323D+00, 0.273477528308839D+00, &
 0.177205532412543D+00, 0.645588935174913D+00, 0.177205532412543D+00, &
 0.617998830908730D-01, 0.876400233818255D+00, 0.617998830908730D-01, &
 0.193909612487010D-01, 0.961218077502598D+00, 0.193909612487010D-01, &
 0.172266687821356D+00, 0.770608554774996D+00, 0.571247574036480D-01, &
 0.770608554774996D+00, 0.172266687821356D+00, 0.571247574036480D-01, &
 0.336861459796345D+00, 0.570222290846683D+00, 0.929162493569720D-01, &
 0.570222290846683D+00, 0.336861459796345D+00, 0.929162493569720D-01, &
 0.298372882136258D+00, 0.686980167808088D+00, 0.146469500556540D-01, &
 0.686980167808088D+00, 0.298372882136258D+00, 0.146469500556540D-01, &
 0.118974497696957D+00, 0.879757171370171D+00, 0.126833093287200D-02, &
 0.879757171370171D+00, 0.118974497696957D+00, 0.126833093287200D-02 /)
gamma_tri(14)%val = (/ &
 0.488963910362179D+00, 0.488963910362179D+00, 0.220721792756430D-01, &
 0.417644719340454D+00, 0.417644719340454D+00, 0.164710561319092D+00, &
 0.273477528308839D+00, 0.273477528308839D+00, 0.453044943382323D+00, &
 0.177205532412543D+00, 0.177205532412543D+00, 0.645588935174913D+00, &
 0.617998830908730D-01, 0.617998830908730D-01, 0.876400233818255D+00, &
 0.193909612487010D-01, 0.193909612487010D-01, 0.961218077502598D+00, &
 0.770608554774996D+00, 0.172266687821356D+00, 0.770608554774996D+00, &
 0.571247574036480D-01, 0.571247574036480D-01, 0.172266687821356D+00, &
 0.570222290846683D+00, 0.336861459796345D+00, 0.570222290846683D+00, &
 0.929162493569720D-01, 0.929162493569720D-01, 0.336861459796345D+00, &
 0.686980167808088D+00, 0.298372882136258D+00, 0.686980167808088D+00, &
 0.146469500556540D-01, 0.146469500556540D-01, 0.298372882136258D+00, &
 0.879757171370171D+00, 0.118974497696957D+00, 0.879757171370171D+00, &
 0.126833093287200D-02, 0.126833093287200D-02, 0.118974497696957D+00 /)

weight_tri(15)%val = (/ &
 0.191687564284900D-02, 0.191687564284900D-02, 0.191687564284900D-02, &
 0.442490272711450D-01, 0.442490272711450D-01, 0.442490272711450D-01, &
 0.511865487188520D-01, 0.511865487188520D-01, 0.511865487188520D-01, &
 0.236877358706880D-01, 0.236877358706880D-01, 0.236877358706880D-01, &
 0.132897756900210D-01, 0.132897756900210D-01, 0.132897756900210D-01, &
 0.474891660819200D-02, 0.474891660819200D-02, 0.474891660819200D-02, &
 0.385500725995930D-01, 0.385500725995930D-01, 0.385500725995930D-01, &
 0.385500725995930D-01, 0.385500725995930D-01, 0.385500725995930D-01, &
 0.272158143206240D-01, 0.272158143206240D-01, 0.272158143206240D-01, &
 0.272158143206240D-01, 0.272158143206240D-01, 0.272158143206240D-01, &
 0.218207736679700D-02, 0.218207736679700D-02, 0.218207736679700D-02, &
 0.218207736679700D-02, 0.218207736679700D-02, 0.218207736679700D-02, &
 0.215053198477310D-01, 0.215053198477310D-01, 0.215053198477310D-01, &
 0.215053198477310D-01, 0.215053198477310D-01, 0.215053198477310D-01, &
 0.767394263104900D-02, 0.767394263104900D-02, 0.767394263104900D-02, &
 0.767394263104900D-02, 0.767394263104900D-02, 0.767394263104900D-02 /)
alpha_tri(15)%val = (/ &
 -.139458337164860D-01, 0.506972916858243D+00, 0.506972916858243D+00, &
 0.137187291433955D+00, 0.431406354283023D+00, 0.431406354283023D+00, &
 0.444612710305711D+00, 0.277693644847144D+00, 0.277693644847144D+00, &
 0.747070217917492D+00, 0.126464891041254D+00, 0.126464891041254D+00, &
 0.858383228050628D+00, 0.708083859746860D-01, 0.708083859746860D-01, &
 0.962069659517853D+00, 0.189651702410730D-01, 0.189651702410730D-01, &
 0.133734161966621D+00, 0.133734161966621D+00, 0.261311371140087D+00, &
 0.261311371140087D+00, 0.604954466893291D+00, 0.604954466893291D+00, &
 0.363666773969170D-01, 0.363666773969170D-01, 0.388046767090269D+00, &
 0.388046767090269D+00, 0.575586555512814D+00, 0.575586555512814D+00, &
 -.101748831265710D-01, -.101748831265710D-01, 0.285712220049916D+00, &
 0.285712220049916D+00, 0.724462663076655D+00, 0.724462663076655D+00, &
 0.368438698758780D-01, 0.368438698758780D-01, 0.215599664072284D+00, &
 0.215599664072284D+00, 0.747556466051838D+00, 0.747556466051838D+00, &
 0.124598093311990D-01, 0.124598093311990D-01, 0.103575616576386D+00, &
 0.103575616576386D+00, 0.883964574092416D+00, 0.883964574092416D+00 /)
beta_tri(15)%val = (/ &
 0.506972916858243D+00, -.139458337164860D-01, 0.506972916858243D+00, &
 0.431406354283023D+00, 0.137187291433955D+00, 0.431406354283023D+00, &
 0.277693644847144D+00, 0.444612710305711D+00, 0.277693644847144D+00, &
 0.126464891041254D+00, 0.747070217917492D+00, 0.126464891041254D+00, &
 0.708083859746860D-01, 0.858383228050628D+00, 0.708083859746860D-01, &
 0.189651702410730D-01, 0.962069659517853D+00, 0.189651702410730D-01, &
 0.261311371140087D+00, 0.604954466893291D+00, 0.133734161966621D+00, &
 0.604954466893291D+00, 0.261311371140087D+00, 0.133734161966621D+00, &
 0.388046767090269D+00, 0.575586555512814D+00, 0.363666773969170D-01, &
 0.575586555512814D+00, 0.388046767090269D+00, 0.363666773969170D-01, &
 0.285712220049916D+00, 0.724462663076655D+00, -.101748831265710D-01, &
 0.724462663076655D+00, 0.285712220049916D+00, -.101748831265710D-01, &
 0.215599664072284D+00, 0.747556466051838D+00, 0.368438698758780D-01, &
 0.747556466051838D+00, 0.215599664072284D+00, 0.368438698758780D-01, &
 0.103575616576386D+00, 0.883964574092416D+00, 0.124598093311990D-01, &
 0.883964574092416D+00, 0.103575616576386D+00, 0.124598093311990D-01 /)
gamma_tri(15)%val = (/ &
 0.506972916858243D+00, 0.506972916858243D+00, -.139458337164860D-01, &
 0.431406354283023D+00, 0.431406354283023D+00, 0.137187291433955D+00, &
 0.277693644847144D+00, 0.277693644847144D+00, 0.444612710305711D+00, &
 0.126464891041254D+00, 0.126464891041254D+00, 0.747070217917492D+00, &
 0.708083859746860D-01, 0.708083859746860D-01, 0.858383228050628D+00, &
 0.189651702410730D-01, 0.189651702410730D-01, 0.962069659517853D+00, &
 0.604954466893291D+00, 0.261311371140087D+00, 0.604954466893291D+00, &
 0.133734161966621D+00, 0.133734161966621D+00, 0.261311371140087D+00, &
 0.575586555512814D+00, 0.388046767090269D+00, 0.575586555512814D+00, &
 0.363666773969170D-01, 0.363666773969170D-01, 0.388046767090269D+00, &
 0.724462663076655D+00, 0.285712220049916D+00, 0.724462663076655D+00, &
 -.101748831265710D-01, -.101748831265710D-01, 0.285712220049916D+00, &
 0.747556466051838D+00, 0.215599664072284D+00, 0.747556466051838D+00, &
 0.368438698758780D-01, 0.368438698758780D-01, 0.215599664072284D+00, &
 0.883964574092416D+00, 0.103575616576386D+00, 0.883964574092416D+00, &
 0.124598093311990D-01, 0.124598093311990D-01, 0.103575616576386D+00 /)

weight_tri(16)%val = (/ &
 0.468756974276420D-01, 0.640587857858500D-02, 0.640587857858500D-02, &
 0.640587857858500D-02, 0.417102967393870D-01, 0.417102967393870D-01, &
 0.417102967393870D-01, 0.268914842500640D-01, 0.268914842500640D-01, &
 0.268914842500640D-01, 0.421325227616500D-01, 0.421325227616500D-01, &
 0.421325227616500D-01, 0.300002668427730D-01, 0.300002668427730D-01, &
 0.300002668427730D-01, 0.142000989250240D-01, 0.142000989250240D-01, &
 0.142000989250240D-01, 0.358246235127300D-02, 0.358246235127300D-02, &
 0.358246235127300D-02, 0.327731474606270D-01, 0.327731474606270D-01, &
 0.327731474606270D-01, 0.327731474606270D-01, 0.327731474606270D-01, &
 0.327731474606270D-01, 0.152983062484410D-01, 0.152983062484410D-01, &
 0.152983062484410D-01, 0.152983062484410D-01, 0.152983062484410D-01, &
 0.152983062484410D-01, 0.238624419283900D-02, 0.238624419283900D-02, &
 0.238624419283900D-02, 0.238624419283900D-02, 0.238624419283900D-02, &
 0.238624419283900D-02, 0.190847927558990D-01, 0.190847927558990D-01, &
 0.190847927558990D-01, 0.190847927558990D-01, 0.190847927558990D-01, &
 0.190847927558990D-01, 0.685005454654200D-02, 0.685005454654200D-02, &
 0.685005454654200D-02, 0.685005454654200D-02, 0.685005454654200D-02, &
 0.685005454654200D-02 /)
alpha_tri(16)%val = (/ &
 0.333333333333333D+00, 0.523891610312300D-02, 0.497380541948438D+00, &
 0.497380541948438D+00, 0.173061122901295D+00, 0.413469438549352D+00, &
 0.413469438549352D+00, 0.590828018660170D-01, 0.470458599066991D+00, &
 0.470458599066991D+00, 0.518892500060958D+00, 0.240553749969521D+00, &
 0.240553749969521D+00, 0.704068411554854D+00, 0.147965794222573D+00, &
 0.147965794222573D+00, 0.849069624685052D+00, 0.754651876574740D-01, &
 0.754651876574740D-01, 0.966807194753950D+00, 0.165964026230250D-01, &
 0.165964026230250D-01, 0.103575692245252D+00, 0.103575692245252D+00, &
 0.296555596579887D+00, 0.296555596579887D+00, 0.599868711174861D+00, &
 0.599868711174861D+00, 0.200834116554160D-01, 0.200834116554160D-01, &
 0.337723063403079D+00, 0.337723063403079D+00, 0.642193524941505D+00, &
 0.642193524941505D+00, -.434100261413900D-02, -.434100261413900D-02, &
 0.204748281642812D+00, 0.204748281642812D+00, 0.799592720971327D+00, &
 0.799592720971327D+00, 0.419417864680100D-01, 0.419417864680100D-01, &
 0.189358492130623D+00, 0.189358492130623D+00, 0.768699721401368D+00, &
 0.768699721401368D+00, 0.143173202306810D-01, 0.143173202306810D-01, &
 0.852836156826570D-01, 0.852836156826570D-01, 0.900399064086661D+00, &
 0.900399064086661D+00 /)
beta_tri(16)%val = (/ &
 0.333333333333333D+00, 0.497380541948438D+00, 0.523891610312300D-02, &
 0.497380541948438D+00, 0.413469438549352D+00, 0.173061122901295D+00, &
 0.413469438549352D+00, 0.470458599066991D+00, 0.590828018660170D-01, &
 0.470458599066991D+00, 0.240553749969521D+00, 0.518892500060958D+00, &
 0.240553749969521D+00, 0.147965794222573D+00, 0.704068411554854D+00, &
 0.147965794222573D+00, 0.754651876574740D-01, 0.849069624685052D+00, &
 0.754651876574740D-01, 0.165964026230250D-01, 0.966807194753950D+00, &
 0.165964026230250D-01, 0.296555596579887D+00, 0.599868711174861D+00, &
 0.103575692245252D+00, 0.599868711174861D+00, 0.296555596579887D+00, &
 0.103575692245252D+00, 0.337723063403079D+00, 0.642193524941505D+00, &
 0.200834116554160D-01, 0.642193524941505D+00, 0.337723063403079D+00, &
 0.200834116554160D-01, 0.204748281642812D+00, 0.799592720971327D+00, &
 -.434100261413900D-02, 0.799592720971327D+00, 0.204748281642812D+00, &
 -.434100261413900D-02, 0.189358492130623D+00, 0.768699721401368D+00, &
 0.419417864680100D-01, 0.768699721401368D+00, 0.189358492130623D+00, &
 0.419417864680100D-01, 0.852836156826570D-01, 0.900399064086661D+00, &
 0.143173202306810D-01, 0.900399064086661D+00, 0.852836156826570D-01, &
 0.143173202306810D-01 /)
gamma_tri(16)%val = (/ &
 0.333333333333333D+00, 0.497380541948438D+00, 0.497380541948438D+00, &
 0.523891610312300D-02, 0.413469438549352D+00, 0.413469438549352D+00, &
 0.173061122901295D+00, 0.470458599066991D+00, 0.470458599066991D+00, &
 0.590828018660170D-01, 0.240553749969521D+00, 0.240553749969521D+00, &
 0.518892500060958D+00, 0.147965794222573D+00, 0.147965794222573D+00, &
 0.704068411554854D+00, 0.754651876574740D-01, 0.754651876574740D-01, &
 0.849069624685052D+00, 0.165964026230250D-01, 0.165964026230250D-01, &
 0.966807194753950D+00, 0.599868711174861D+00, 0.296555596579887D+00, &
 0.599868711174861D+00, 0.103575692245252D+00, 0.103575692245252D+00, &
 0.296555596579887D+00, 0.642193524941505D+00, 0.337723063403079D+00, &
 0.642193524941505D+00, 0.200834116554160D-01, 0.200834116554160D-01, &
 0.337723063403079D+00, 0.799592720971327D+00, 0.204748281642812D+00, &
 0.799592720971327D+00, -.434100261413900D-02, -.434100261413900D-02, &
 0.204748281642812D+00, 0.768699721401368D+00, 0.189358492130623D+00, &
 0.768699721401368D+00, 0.419417864680100D-01, 0.419417864680100D-01, &
 0.189358492130623D+00, 0.900399064086661D+00, 0.852836156826570D-01, &
 0.900399064086661D+00, 0.143173202306810D-01, 0.143173202306810D-01, &
 0.852836156826570D-01 /)

weight_tri(17)%val = (/ &
 0.334371992908030D-01, 0.509341544050700D-02, 0.509341544050700D-02, &
 0.509341544050700D-02, 0.146708645276380D-01, 0.146708645276380D-01, &
 0.146708645276380D-01, 0.243508783536720D-01, 0.243508783536720D-01, &
 0.243508783536720D-01, 0.311075508689690D-01, 0.311075508689690D-01, &
 0.311075508689690D-01, 0.312571112186200D-01, 0.312571112186200D-01, &
 0.312571112186200D-01, 0.248156543396650D-01, 0.248156543396650D-01, &
 0.248156543396650D-01, 0.140560730705570D-01, 0.140560730705570D-01, &
 0.140560730705570D-01, 0.319467617377900D-02, 0.319467617377900D-02, &
 0.319467617377900D-02, 0.811965531899300D-02, 0.811965531899300D-02, &
 0.811965531899300D-02, 0.811965531899300D-02, 0.811965531899300D-02, &
 0.811965531899300D-02, 0.268057422831630D-01, 0.268057422831630D-01, &
 0.268057422831630D-01, 0.268057422831630D-01, 0.268057422831630D-01, &
 0.268057422831630D-01, 0.184599932108220D-01, 0.184599932108220D-01, &
 0.184599932108220D-01, 0.184599932108220D-01, 0.184599932108220D-01, &
 0.184599932108220D-01, 0.847686853432800D-02, 0.847686853432800D-02, &
 0.847686853432800D-02, 0.847686853432800D-02, 0.847686853432800D-02, &
 0.847686853432800D-02, 0.182927967700250D-01, 0.182927967700250D-01, &
 0.182927967700250D-01, 0.182927967700250D-01, 0.182927967700250D-01, &
 0.182927967700250D-01, 0.666563200416500D-02, 0.666563200416500D-02, &
 0.666563200416500D-02, 0.666563200416500D-02, 0.666563200416500D-02, &
 0.666563200416500D-02 /)
alpha_tri(17)%val = (/ &
 0.333333333333333D+00, 0.565891888645200D-02, 0.497170540556774D+00, &
 0.497170540556774D+00, 0.356473547507510D-01, 0.482176322624625D+00, &
 0.482176322624625D+00, 0.995200619584370D-01, 0.450239969020782D+00, &
 0.450239969020782D+00, 0.199467521245206D+00, 0.400266239377397D+00, &
 0.400266239377397D+00, 0.495717464058095D+00, 0.252141267970953D+00, &
 0.252141267970953D+00, 0.675905990683077D+00, 0.162047004658461D+00, &
 0.162047004658461D+00, 0.848248235478508D+00, 0.758758822607460D-01, &
 0.758758822607460D-01, 0.968690546064356D+00, 0.156547269678220D-01, &
 0.156547269678220D-01, 0.101869288269190D-01, 0.101869288269190D-01, &
 0.334319867363658D+00, 0.334319867363658D+00, 0.655493203809423D+00, &
 0.655493203809423D+00, 0.135440871671036D+00, 0.135440871671036D+00, &
 0.292221537796944D+00, 0.292221537796944D+00, 0.572337590532020D+00, &
 0.572337590532020D+00, 0.544239242905830D-01, 0.544239242905830D-01, &
 0.319574885423190D+00, 0.319574885423190D+00, 0.626001190286228D+00, &
 0.626001190286228D+00, 0.128685608336370D-01, 0.128685608336370D-01, &
 0.190704224192292D+00, 0.190704224192292D+00, 0.796427214974071D+00, &
 0.796427214974071D+00, 0.671657824135240D-01, 0.671657824135240D-01, &
 0.180483211648746D+00, 0.180483211648746D+00, 0.752351005937729D+00, &
 0.752351005937729D+00, 0.146631822248280D-01, 0.146631822248280D-01, &
 0.807113136795640D-01, 0.807113136795640D-01, 0.904625504095608D+00, &
 0.904625504095608D+00 /)
beta_tri(17)%val = (/ &
 0.333333333333333D+00, 0.497170540556774D+00, 0.565891888645200D-02, &
 0.497170540556774D+00, 0.482176322624625D+00, 0.356473547507510D-01, &
 0.482176322624625D+00, 0.450239969020782D+00, 0.995200619584370D-01, &
 0.450239969020782D+00, 0.400266239377397D+00, 0.199467521245206D+00, &
 0.400266239377397D+00, 0.252141267970953D+00, 0.495717464058095D+00, &
 0.252141267970953D+00, 0.162047004658461D+00, 0.675905990683077D+00, &
 0.162047004658461D+00, 0.758758822607460D-01, 0.848248235478508D+00, &
 0.758758822607460D-01, 0.156547269678220D-01, 0.968690546064356D+00, &
 0.156547269678220D-01, 0.334319867363658D+00, 0.655493203809423D+00, &
 0.101869288269190D-01, 0.655493203809423D+00, 0.334319867363658D+00, &
 0.101869288269190D-01, 0.292221537796944D+00, 0.572337590532020D+00, &
 0.135440871671036D+00, 0.572337590532020D+00, 0.292221537796944D+00, &
 0.135440871671036D+00, 0.319574885423190D+00, 0.626001190286228D+00, &
 0.544239242905830D-01, 0.626001190286228D+00, 0.319574885423190D+00, &
 0.544239242905830D-01, 0.190704224192292D+00, 0.796427214974071D+00, &
 0.128685608336370D-01, 0.796427214974071D+00, 0.190704224192292D+00, &
 0.128685608336370D-01, 0.180483211648746D+00, 0.752351005937729D+00, &
 0.671657824135240D-01, 0.752351005937729D+00, 0.180483211648746D+00, &
 0.671657824135240D-01, 0.807113136795640D-01, 0.904625504095608D+00, &
 0.146631822248280D-01, 0.904625504095608D+00, 0.807113136795640D-01, &
 0.146631822248280D-01 /)
gamma_tri(17)%val = (/ &
 0.333333333333333D+00, 0.497170540556774D+00, 0.497170540556774D+00, &
 0.565891888645200D-02, 0.482176322624625D+00, 0.482176322624625D+00, &
 0.356473547507510D-01, 0.450239969020782D+00, 0.450239969020782D+00, &
 0.995200619584370D-01, 0.400266239377397D+00, 0.400266239377397D+00, &
 0.199467521245206D+00, 0.252141267970953D+00, 0.252141267970953D+00, &
 0.495717464058095D+00, 0.162047004658461D+00, 0.162047004658461D+00, &
 0.675905990683077D+00, 0.758758822607460D-01, 0.758758822607460D-01, &
 0.848248235478508D+00, 0.156547269678220D-01, 0.156547269678220D-01, &
 0.968690546064356D+00, 0.655493203809423D+00, 0.334319867363658D+00, &
 0.655493203809423D+00, 0.101869288269190D-01, 0.101869288269190D-01, &
 0.334319867363658D+00, 0.572337590532020D+00, 0.292221537796944D+00, &
 0.572337590532020D+00, 0.135440871671036D+00, 0.135440871671036D+00, &
 0.292221537796944D+00, 0.626001190286228D+00, 0.319574885423190D+00, &
 0.626001190286228D+00, 0.544239242905830D-01, 0.544239242905830D-01, &
 0.319574885423190D+00, 0.796427214974071D+00, 0.190704224192292D+00, &
 0.796427214974071D+00, 0.128685608336370D-01, 0.128685608336370D-01, &
 0.190704224192292D+00, 0.752351005937729D+00, 0.180483211648746D+00, &
 0.752351005937729D+00, 0.671657824135240D-01, 0.671657824135240D-01, &
 0.180483211648746D+00, 0.904625504095608D+00, 0.807113136795640D-01, &
 0.904625504095608D+00, 0.146631822248280D-01, 0.146631822248280D-01, &
 0.807113136795640D-01 /)

weight_tri(18)%val = (/ &
 0.308099399376470D-01, 0.907243667940400D-02, 0.907243667940400D-02, &
 0.907243667940400D-02, 0.187613169395940D-01, 0.187613169395940D-01, &
 0.187613169395940D-01, 0.194410979854770D-01, 0.194410979854770D-01, &
 0.194410979854770D-01, 0.277539486108100D-01, 0.277539486108100D-01, &
 0.277539486108100D-01, 0.322562253514570D-01, 0.322562253514570D-01, &
 0.322562253514570D-01, 0.250740326169220D-01, 0.250740326169220D-01, &
 0.250740326169220D-01, 0.152719279718320D-01, 0.152719279718320D-01, &
 0.152719279718320D-01, 0.679392202296300D-02, 0.679392202296300D-02, &
 0.679392202296300D-02, -.222309872992000D-02, -.222309872992000D-02, &
 -.222309872992000D-02, 0.633191407640600D-02, 0.633191407640600D-02, &
 0.633191407640600D-02, 0.633191407640600D-02, 0.633191407640600D-02, &
 0.633191407640600D-02, 0.272575380491380D-01, 0.272575380491380D-01, &
 0.272575380491380D-01, 0.272575380491380D-01, 0.272575380491380D-01, &
 0.272575380491380D-01, 0.176767856494650D-01, 0.176767856494650D-01, &
 0.176767856494650D-01, 0.176767856494650D-01, 0.176767856494650D-01, &
 0.176767856494650D-01, 0.183794846380700D-01, 0.183794846380700D-01, &
 0.183794846380700D-01, 0.183794846380700D-01, 0.183794846380700D-01, &
 0.183794846380700D-01, 0.810473280819200D-02, 0.810473280819200D-02, &
 0.810473280819200D-02, 0.810473280819200D-02, 0.810473280819200D-02, &
 0.810473280819200D-02, 0.763412907072500D-02, 0.763412907072500D-02, &
 0.763412907072500D-02, 0.763412907072500D-02, 0.763412907072500D-02, &
 0.763412907072500D-02, 0.461876607940000D-04, 0.461876607940000D-04, &
 0.461876607940000D-04, 0.461876607940000D-04, 0.461876607940000D-04, &
 0.461876607940000D-04 /)
alpha_tri(18)%val = (/ &
 0.333333333333333D+00, 0.133103827381570D-01, 0.493344808630921D+00, &
 0.493344808630921D+00, 0.615788115160860D-01, 0.469210594241957D+00, &
 0.469210594241957D+00, 0.127437208225989D+00, 0.436281395887006D+00, &
 0.436281395887006D+00, 0.210307658653168D+00, 0.394846170673416D+00, &
 0.394846170673416D+00, 0.500410862393686D+00, 0.249794568803157D+00, &
 0.249794568803157D+00, 0.677135612512315D+00, 0.161432193743843D+00, &
 0.161432193743843D+00, 0.846803545029257D+00, 0.765982274853710D-01, &
 0.765982274853710D-01, 0.951495121293100D+00, 0.242524393534500D-01, &
 0.242524393534500D-01, 0.913707265566071D+00, 0.431463672169650D-01, &
 0.431463672169650D-01, 0.843053620242000D-02, 0.843053620242000D-02, &
 0.358911494940944D+00, 0.358911494940944D+00, 0.632657968856636D+00, &
 0.632657968856636D+00, 0.131186551737188D+00, 0.131186551737188D+00, &
 0.294402476751957D+00, 0.294402476751957D+00, 0.574410971510855D+00, &
 0.574410971510855D+00, 0.502031515656750D-01, 0.502031515656750D-01, &
 0.325017801641814D+00, 0.325017801641814D+00, 0.624779046792512D+00, &
 0.624779046792512D+00, 0.663292638109160D-01, 0.663292638109160D-01, &
 0.184737559666046D+00, 0.184737559666046D+00, 0.748933176523037D+00, &
 0.748933176523037D+00, 0.119961945662360D-01, 0.119961945662360D-01, &
 0.218796800013321D+00, 0.218796800013321D+00, 0.769207005420443D+00, &
 0.769207005420443D+00, 0.148581005901250D-01, 0.148581005901250D-01, &
 0.101179597136408D+00, 0.101179597136408D+00, 0.883962302273467D+00, &
 0.883962302273467D+00, -.352220152879490D-01, -.352220152879490D-01, &
 0.208747552825860D-01, 0.208747552825860D-01, 0.101434726000536D+01, &
 0.101434726000536D+01 /)
beta_tri(18)%val = (/ &
 0.333333333333333D+00, 0.493344808630921D+00, 0.133103827381570D-01, &
 0.493344808630921D+00, 0.469210594241957D+00, 0.615788115160860D-01, &
 0.469210594241957D+00, 0.436281395887006D+00, 0.127437208225989D+00, &
 0.436281395887006D+00, 0.394846170673416D+00, 0.210307658653168D+00, &
 0.394846170673416D+00, 0.249794568803157D+00, 0.500410862393686D+00, &
 0.249794568803157D+00, 0.161432193743843D+00, 0.677135612512315D+00, &
 0.161432193743843D+00, 0.765982274853710D-01, 0.846803545029257D+00, &
 0.765982274853710D-01, 0.242524393534500D-01, 0.951495121293100D+00, &
 0.242524393534500D-01, 0.431463672169650D-01, 0.913707265566071D+00, &
 0.431463672169650D-01, 0.358911494940944D+00, 0.632657968856636D+00, &
 0.843053620242000D-02, 0.632657968856636D+00, 0.358911494940944D+00, &
 0.843053620242000D-02, 0.294402476751957D+00, 0.574410971510855D+00, &
 0.131186551737188D+00, 0.574410971510855D+00, 0.294402476751957D+00, &
 0.131186551737188D+00, 0.325017801641814D+00, 0.624779046792512D+00, &
 0.502031515656750D-01, 0.624779046792512D+00, 0.325017801641814D+00, &
 0.502031515656750D-01, 0.184737559666046D+00, 0.748933176523037D+00, &
 0.663292638109160D-01, 0.748933176523037D+00, 0.184737559666046D+00, &
 0.663292638109160D-01, 0.218796800013321D+00, 0.769207005420443D+00, &
 0.119961945662360D-01, 0.769207005420443D+00, 0.218796800013321D+00, &
 0.119961945662360D-01, 0.101179597136408D+00, 0.883962302273467D+00, &
 0.148581005901250D-01, 0.883962302273467D+00, 0.101179597136408D+00, &
 0.148581005901250D-01, 0.208747552825860D-01, 0.101434726000536D+01, &
 -.352220152879490D-01, 0.101434726000536D+01, 0.208747552825860D-01, &
 -.352220152879490D-01 /)
gamma_tri(18)%val = (/ &
 0.333333333333333D+00, 0.493344808630921D+00, 0.493344808630921D+00, &
 0.133103827381570D-01, 0.469210594241957D+00, 0.469210594241957D+00, &
 0.615788115160860D-01, 0.436281395887006D+00, 0.436281395887006D+00, &
 0.127437208225989D+00, 0.394846170673416D+00, 0.394846170673416D+00, &
 0.210307658653168D+00, 0.249794568803157D+00, 0.249794568803157D+00, &
 0.500410862393686D+00, 0.161432193743843D+00, 0.161432193743843D+00, &
 0.677135612512315D+00, 0.765982274853710D-01, 0.765982274853710D-01, &
 0.846803545029257D+00, 0.242524393534500D-01, 0.242524393534500D-01, &
 0.951495121293100D+00, 0.431463672169650D-01, 0.431463672169650D-01, &
 0.913707265566071D+00, 0.632657968856636D+00, 0.358911494940944D+00, &
 0.632657968856636D+00, 0.843053620242000D-02, 0.843053620242000D-02, &
 0.358911494940944D+00, 0.574410971510855D+00, 0.294402476751957D+00, &
 0.574410971510855D+00, 0.131186551737188D+00, 0.131186551737188D+00, &
 0.294402476751957D+00, 0.624779046792512D+00, 0.325017801641814D+00, &
 0.624779046792512D+00, 0.502031515656750D-01, 0.502031515656750D-01, &
 0.325017801641814D+00, 0.748933176523037D+00, 0.184737559666046D+00, &
 0.748933176523037D+00, 0.663292638109160D-01, 0.663292638109160D-01, &
 0.184737559666046D+00, 0.769207005420443D+00, 0.218796800013321D+00, &
 0.769207005420443D+00, 0.119961945662360D-01, 0.119961945662360D-01, &
 0.218796800013321D+00, 0.883962302273467D+00, 0.101179597136408D+00, &
 0.883962302273467D+00, 0.148581005901250D-01, 0.148581005901250D-01, &
 0.101179597136408D+00, 0.101434726000536D+01, 0.208747552825860D-01, &
 0.101434726000536D+01, -.352220152879490D-01, -.352220152879490D-01, &
 0.208747552825860D-01 /)

weight_tri(19)%val = (/ &
 0.329063313889190D-01, 0.103307318912720D-01, 0.103307318912720D-01, &
 0.103307318912720D-01, 0.223872472630160D-01, 0.223872472630160D-01, &
 0.223872472630160D-01, 0.302661258694680D-01, 0.302661258694680D-01, &
 0.302661258694680D-01, 0.304909678021980D-01, 0.304909678021980D-01, &
 0.304909678021980D-01, 0.241592127416410D-01, 0.241592127416410D-01, &
 0.241592127416410D-01, 0.160508035868010D-01, 0.160508035868010D-01, &
 0.160508035868010D-01, 0.808458026178400D-02, 0.808458026178400D-02, &
 0.808458026178400D-02, 0.207936202748500D-02, 0.207936202748500D-02, &
 0.207936202748500D-02, 0.388487690498100D-02, 0.388487690498100D-02, &
 0.388487690498100D-02, 0.388487690498100D-02, 0.388487690498100D-02, &
 0.388487690498100D-02, 0.255741606120220D-01, 0.255741606120220D-01, &
 0.255741606120220D-01, 0.255741606120220D-01, 0.255741606120220D-01, &
 0.255741606120220D-01, 0.888090357333800D-02, 0.888090357333800D-02, &
 0.888090357333800D-02, 0.888090357333800D-02, 0.888090357333800D-02, &
 0.888090357333800D-02, 0.161245467617310D-01, 0.161245467617310D-01, &
 0.161245467617310D-01, 0.161245467617310D-01, 0.161245467617310D-01, &
 0.161245467617310D-01, 0.249194181749100D-02, 0.249194181749100D-02, &
 0.249194181749100D-02, 0.249194181749100D-02, 0.249194181749100D-02, &
 0.249194181749100D-02, 0.182428401189510D-01, 0.182428401189510D-01, &
 0.182428401189510D-01, 0.182428401189510D-01, 0.182428401189510D-01, &
 0.182428401189510D-01, 0.102585637361990D-01, 0.102585637361990D-01, &
 0.102585637361990D-01, 0.102585637361990D-01, 0.102585637361990D-01, &
 0.102585637361990D-01, 0.379992885530200D-02, 0.379992885530200D-02, &
 0.379992885530200D-02, 0.379992885530200D-02, 0.379992885530200D-02, &
 0.379992885530200D-02 /)
alpha_tri(19)%val = (/ &
 0.333333333333333D+00, 0.207800258539870D-01, 0.489609987073006D+00, &
 0.489609987073006D+00, 0.909262146042150D-01, 0.454536892697893D+00, &
 0.454536892697893D+00, 0.197166638701138D+00, 0.401416680649431D+00, &
 0.401416680649431D+00, 0.488896691193805D+00, 0.255551654403098D+00, &
 0.255551654403098D+00, 0.645844115695741D+00, 0.177077942152130D+00, &
 0.177077942152130D+00, 0.779877893544096D+00, 0.110061053227952D+00, &
 0.110061053227952D+00, 0.888942751496321D+00, 0.555286242518400D-01, &
 0.555286242518400D-01, 0.974756272445543D+00, 0.126218637772290D-01, &
 0.126218637772290D-01, 0.361141784841200D-02, 0.361141784841200D-02, &
 0.395754787356943D+00, 0.395754787356943D+00, 0.600633794794645D+00, &
 0.600633794794645D+00, 0.134466754530780D+00, 0.134466754530780D+00, &
 0.307929983880436D+00, 0.307929983880436D+00, 0.557603261588784D+00, &
 0.557603261588784D+00, 0.144460257761150D-01, 0.144460257761150D-01, &
 0.264566948406520D+00, 0.264566948406520D+00, 0.720987025817365D+00, &
 0.720987025817365D+00, 0.469335788381780D-01, 0.469335788381780D-01, &
 0.358539352205951D+00, 0.358539352205951D+00, 0.594527068955871D+00, &
 0.594527068955871D+00, 0.286112035056700D-02, 0.286112035056700D-02, &
 0.157807405968595D+00, 0.157807405968595D+00, 0.839331473680839D+00, &
 0.839331473680839D+00, 0.223861424097916D+00, 0.223861424097916D+00, &
 0.750505969759110D-01, 0.750505969759110D-01, 0.701087978926173D+00, &
 0.701087978926173D+00, 0.346470748167600D-01, 0.346470748167600D-01, &
 0.142421601113383D+00, 0.142421601113383D+00, 0.822931324069857D+00, &
 0.822931324069857D+00, 0.101611192962780D-01, 0.101611192962780D-01, &
 0.654946280829380D-01, 0.654946280829380D-01, 0.924344252620784D+00, &
 0.924344252620784D+00 /)
beta_tri(19)%val = (/ &
 0.333333333333333D+00, 0.489609987073006D+00, 0.207800258539870D-01, &
 0.489609987073006D+00, 0.454536892697893D+00, 0.909262146042150D-01, &
 0.454536892697893D+00, 0.401416680649431D+00, 0.197166638701138D+00, &
 0.401416680649431D+00, 0.255551654403098D+00, 0.488896691193805D+00, &
 0.255551654403098D+00, 0.177077942152130D+00, 0.645844115695741D+00, &
 0.177077942152130D+00, 0.110061053227952D+00, 0.779877893544096D+00, &
 0.110061053227952D+00, 0.555286242518400D-01, 0.888942751496321D+00, &
 0.555286242518400D-01, 0.126218637772290D-01, 0.974756272445543D+00, &
 0.126218637772290D-01, 0.395754787356943D+00, 0.600633794794645D+00, &
 0.361141784841200D-02, 0.600633794794645D+00, 0.395754787356943D+00, &
 0.361141784841200D-02, 0.307929983880436D+00, 0.557603261588784D+00, &
 0.134466754530780D+00, 0.557603261588784D+00, 0.307929983880436D+00, &
 0.134466754530780D+00, 0.264566948406520D+00, 0.720987025817365D+00, &
 0.144460257761150D-01, 0.720987025817365D+00, 0.264566948406520D+00, &
 0.144460257761150D-01, 0.358539352205951D+00, 0.594527068955871D+00, &
 0.469335788381780D-01, 0.594527068955871D+00, 0.358539352205951D+00, &
 0.469335788381780D-01, 0.157807405968595D+00, 0.839331473680839D+00, &
 0.286112035056700D-02, 0.839331473680839D+00, 0.157807405968595D+00, &
 0.286112035056700D-02, 0.750505969759110D-01, 0.701087978926173D+00, &
 0.223861424097916D+00, 0.701087978926173D+00, 0.750505969759110D-01, &
 0.223861424097916D+00, 0.142421601113383D+00, 0.822931324069857D+00, &
 0.346470748167600D-01, 0.822931324069857D+00, 0.142421601113383D+00, &
 0.346470748167600D-01, 0.654946280829380D-01, 0.924344252620784D+00, &
 0.101611192962780D-01, 0.924344252620784D+00, 0.654946280829380D-01, &
 0.101611192962780D-01 /)
gamma_tri(19)%val = (/ &
 0.333333333333333D+00, 0.489609987073006D+00, 0.489609987073006D+00, &
 0.207800258539870D-01, 0.454536892697893D+00, 0.454536892697893D+00, &
 0.909262146042150D-01, 0.401416680649431D+00, 0.401416680649431D+00, &
 0.197166638701138D+00, 0.255551654403098D+00, 0.255551654403098D+00, &
 0.488896691193805D+00, 0.177077942152130D+00, 0.177077942152130D+00, &
 0.645844115695741D+00, 0.110061053227952D+00, 0.110061053227952D+00, &
 0.779877893544096D+00, 0.555286242518400D-01, 0.555286242518400D-01, &
 0.888942751496321D+00, 0.126218637772290D-01, 0.126218637772290D-01, &
 0.974756272445543D+00, 0.600633794794645D+00, 0.395754787356943D+00, &
 0.600633794794645D+00, 0.361141784841200D-02, 0.361141784841200D-02, &
 0.395754787356943D+00, 0.557603261588784D+00, 0.307929983880436D+00, &
 0.557603261588784D+00, 0.134466754530780D+00, 0.134466754530780D+00, &
 0.307929983880436D+00, 0.720987025817365D+00, 0.264566948406520D+00, &
 0.720987025817365D+00, 0.144460257761150D-01, 0.144460257761150D-01, &
 0.264566948406520D+00, 0.594527068955871D+00, 0.358539352205951D+00, &
 0.594527068955871D+00, 0.469335788381780D-01, 0.469335788381780D-01, &
 0.358539352205951D+00, 0.839331473680839D+00, 0.157807405968595D+00, &
 0.839331473680839D+00, 0.286112035056700D-02, 0.286112035056700D-02, &
 0.157807405968595D+00, 0.701087978926173D+00, 0.750505969759110D-01, &
 0.701087978926173D+00, 0.223861424097916D+00, 0.223861424097916D+00, &
 0.750505969759110D-01, 0.822931324069857D+00, 0.142421601113383D+00, &
 0.822931324069857D+00, 0.346470748167600D-01, 0.346470748167600D-01, &
 0.142421601113383D+00, 0.924344252620784D+00, 0.654946280829380D-01, &
 0.924344252620784D+00, 0.101611192962780D-01, 0.101611192962780D-01, &
 0.654946280829380D-01 /)

weight_tri(20)%val = (/ &
 0.330570555416240D-01, 0.867019185663000D-03, 0.867019185663000D-03, &
 0.867019185663000D-03, 0.116600527164480D-01, 0.116600527164480D-01, &
 0.116600527164480D-01, 0.228769363564210D-01, 0.228769363564210D-01, &
 0.228769363564210D-01, 0.304489826739380D-01, 0.304489826739380D-01, &
 0.304489826739380D-01, 0.306248917253550D-01, 0.306248917253550D-01, &
 0.306248917253550D-01, 0.243680576768000D-01, 0.243680576768000D-01, &
 0.243680576768000D-01, 0.159974320320240D-01, 0.159974320320240D-01, &
 0.159974320320240D-01, 0.769830181560200D-02, 0.769830181560200D-02, &
 0.769830181560200D-02, -.632060497488000D-03, -.632060497488000D-03, &
 -.632060497488000D-03, 0.175113430119300D-02, 0.175113430119300D-02, &
 0.175113430119300D-02, 0.164658391895760D-01, 0.164658391895760D-01, &
 0.164658391895760D-01, 0.164658391895760D-01, 0.164658391895760D-01, &
 0.164658391895760D-01, 0.483903354048500D-02, 0.483903354048500D-02, &
 0.483903354048500D-02, 0.483903354048500D-02, 0.483903354048500D-02, &
 0.483903354048500D-02, 0.258049065346500D-01, 0.258049065346500D-01, &
 0.258049065346500D-01, 0.258049065346500D-01, 0.258049065346500D-01, &
 0.258049065346500D-01, 0.847109105444100D-02, 0.847109105444100D-02, &
 0.847109105444100D-02, 0.847109105444100D-02, 0.847109105444100D-02, &
 0.847109105444100D-02, 0.183549141062800D-01, 0.183549141062800D-01, &
 0.183549141062800D-01, 0.183549141062800D-01, 0.183549141062800D-01, &
 0.183549141062800D-01, 0.704404677908000D-03, 0.704404677908000D-03, &
 0.704404677908000D-03, 0.704404677908000D-03, 0.704404677908000D-03, &
 0.704404677908000D-03, 0.101126849274620D-01, 0.101126849274620D-01, &
 0.101126849274620D-01, 0.101126849274620D-01, 0.101126849274620D-01, &
 0.101126849274620D-01, 0.357390938595000D-02, 0.357390938595000D-02, &
 0.357390938595000D-02, 0.357390938595000D-02, 0.357390938595000D-02, &
 0.357390938595000D-02 /)
alpha_tri(20)%val = (/ &
 0.333333333333333D+00, -.190092870440000D-02, 0.500950464352200D+00, &
 0.500950464352200D+00, 0.235740841305430D-01, 0.488212957934729D+00, &
 0.488212957934729D+00, 0.897266360994350D-01, 0.455136681950283D+00, &
 0.455136681950283D+00, 0.196007481363421D+00, 0.401996259318289D+00, &
 0.401996259318289D+00, 0.488214180481157D+00, 0.255892909759421D+00, &
 0.255892909759421D+00, 0.647023488009788D+00, 0.176488255995106D+00, &
 0.176488255995106D+00, 0.791658289326483D+00, 0.104170855336758D+00, &
 0.104170855336758D+00, 0.893862072318140D+00, 0.530689638409300D-01, &
 0.530689638409300D-01, 0.916762569607942D+00, 0.416187151960290D-01, &
 0.416187151960290D-01, 0.976836157186356D+00, 0.115819214068220D-01, &
 0.115819214068220D-01, 0.487415836648390D-01, 0.487415836648390D-01, &
 0.344855770229001D+00, 0.344855770229001D+00, 0.606402646106160D+00, &
 0.606402646106160D+00, 0.631411594860500D-02, 0.631411594860500D-02, &
 0.377843269594854D+00, 0.377843269594854D+00, 0.615842614456541D+00, &
 0.615842614456541D+00, 0.134316520547348D+00, 0.134316520547348D+00, &
 0.306635479062357D+00, 0.306635479062357D+00, 0.559048000390295D+00, &
 0.559048000390295D+00, 0.139738939623920D-01, 0.139738939623920D-01, &
 0.249419362774742D+00, 0.249419362774742D+00, 0.736606743262866D+00, &
 0.736606743262866D+00, 0.755491329097640D-01, 0.755491329097640D-01, &
 0.212775724802802D+00, 0.212775724802802D+00, 0.711675142287434D+00, &
 0.711675142287434D+00, -.836815320822700D-02, -.836815320822700D-02, &
 0.146965436053239D+00, 0.146965436053239D+00, 0.861402717154987D+00, &
 0.861402717154987D+00, 0.266860632587140D-01, 0.266860632587140D-01, &
 0.137726978828923D+00, 0.137726978828923D+00, 0.835586957912363D+00, &
 0.835586957912363D+00, 0.105477192941410D-01, 0.105477192941410D-01, &
 0.596961091490070D-01, 0.596961091490070D-01, 0.929756171556853D+00, &
 0.929756171556853D+00 /)
beta_tri(20)%val = (/ &
 0.333333333333333D+00, 0.500950464352200D+00, -.190092870440000D-02, &
 0.500950464352200D+00, 0.488212957934729D+00, 0.235740841305430D-01, &
 0.488212957934729D+00, 0.455136681950283D+00, 0.897266360994350D-01, &
 0.455136681950283D+00, 0.401996259318289D+00, 0.196007481363421D+00, &
 0.401996259318289D+00, 0.255892909759421D+00, 0.488214180481157D+00, &
 0.255892909759421D+00, 0.176488255995106D+00, 0.647023488009788D+00, &
 0.176488255995106D+00, 0.104170855336758D+00, 0.791658289326483D+00, &
 0.104170855336758D+00, 0.530689638409300D-01, 0.893862072318140D+00, &
 0.530689638409300D-01, 0.416187151960290D-01, 0.916762569607942D+00, &
 0.416187151960290D-01, 0.115819214068220D-01, 0.976836157186356D+00, &
 0.115819214068220D-01, 0.344855770229001D+00, 0.606402646106160D+00, &
 0.487415836648390D-01, 0.606402646106160D+00, 0.344855770229001D+00, &
 0.487415836648390D-01, 0.377843269594854D+00, 0.615842614456541D+00, &
 0.631411594860500D-02, 0.615842614456541D+00, 0.377843269594854D+00, &
 0.631411594860500D-02, 0.306635479062357D+00, 0.559048000390295D+00, &
 0.134316520547348D+00, 0.559048000390295D+00, 0.306635479062357D+00, &
 0.134316520547348D+00, 0.249419362774742D+00, 0.736606743262866D+00, &
 0.139738939623920D-01, 0.736606743262866D+00, 0.249419362774742D+00, &
 0.139738939623920D-01, 0.212775724802802D+00, 0.711675142287434D+00, &
 0.755491329097640D-01, 0.711675142287434D+00, 0.212775724802802D+00, &
 0.755491329097640D-01, 0.146965436053239D+00, 0.861402717154987D+00, &
 -.836815320822700D-02, 0.861402717154987D+00, 0.146965436053239D+00, &
 -.836815320822700D-02, 0.137726978828923D+00, 0.835586957912363D+00, &
 0.266860632587140D-01, 0.835586957912363D+00, 0.137726978828923D+00, &
 0.266860632587140D-01, 0.596961091490070D-01, 0.929756171556853D+00, &
 0.105477192941410D-01, 0.929756171556853D+00, 0.596961091490070D-01, &
 0.105477192941410D-01 /)
gamma_tri(20)%val = (/ &
 0.333333333333333D+00, 0.500950464352200D+00, 0.500950464352200D+00, &
 -.190092870440000D-02, 0.488212957934729D+00, 0.488212957934729D+00, &
 0.235740841305430D-01, 0.455136681950283D+00, 0.455136681950283D+00, &
 0.897266360994350D-01, 0.401996259318289D+00, 0.401996259318289D+00, &
 0.196007481363421D+00, 0.255892909759421D+00, 0.255892909759421D+00, &
 0.488214180481157D+00, 0.176488255995106D+00, 0.176488255995106D+00, &
 0.647023488009788D+00, 0.104170855336758D+00, 0.104170855336758D+00, &
 0.791658289326483D+00, 0.530689638409300D-01, 0.530689638409300D-01, &
 0.893862072318140D+00, 0.416187151960290D-01, 0.416187151960290D-01, &
 0.916762569607942D+00, 0.115819214068220D-01, 0.115819214068220D-01, &
 0.976836157186356D+00, 0.606402646106160D+00, 0.344855770229001D+00, &
 0.606402646106160D+00, 0.487415836648390D-01, 0.487415836648390D-01, &
 0.344855770229001D+00, 0.615842614456541D+00, 0.377843269594854D+00, &
 0.615842614456541D+00, 0.631411594860500D-02, 0.631411594860500D-02, &
 0.377843269594854D+00, 0.559048000390295D+00, 0.306635479062357D+00, &
 0.559048000390295D+00, 0.134316520547348D+00, 0.134316520547348D+00, &
 0.306635479062357D+00, 0.736606743262866D+00, 0.249419362774742D+00, &
 0.736606743262866D+00, 0.139738939623920D-01, 0.139738939623920D-01, &
 0.249419362774742D+00, 0.711675142287434D+00, 0.212775724802802D+00, &
 0.711675142287434D+00, 0.755491329097640D-01, 0.755491329097640D-01, &
 0.212775724802802D+00, 0.861402717154987D+00, 0.146965436053239D+00, &
 0.861402717154987D+00, -.836815320822700D-02, -.836815320822700D-02, &
 0.146965436053239D+00, 0.835586957912363D+00, 0.137726978828923D+00, &
 0.835586957912363D+00, 0.266860632587140D-01, 0.266860632587140D-01, &
 0.137726978828923D+00, 0.929756171556853D+00, 0.596961091490070D-01, &
 0.929756171556853D+00, 0.105477192941410D-01, 0.105477192941410D-01, &
 0.596961091490070D-01 /)

end subroutine setup_tri

end module quadrature_rules
