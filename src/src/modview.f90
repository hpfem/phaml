module view_modifier

! This module provides facilities to modify the view in an OpenGL window.
! The mouse buttons and keyboard arrow keys can be used to zoom, pan,
! rotate and change the scale.  A menu or submenu can be used to select which
! buttons perform which function and to reset the view to the initial settings.
! Handles up to MAXWIN windows simultaneously; reset this symbolic constant
! below to increase the number of windows.

! William F. Mitchell
! william.mitchell@nist.gov
! Mathematical and Computational Sciences Division
! National Institute of Standards and Technology
! April, 1998
! October, 2000 added multiwindow capability

! To use this module:
!
! 1) put a USE view_modifier statement in any program unit that calls a
!    procedure in this module
!
! 2) set the initial operation assignments and scale below the
!    "Initial configuration" comment below
!
! 3) call view_modifier_init after glutCreateWindow
!    It takes six real(gldouble) arguments to give the point you are looking
!    from and the point you are looking at.
!    This is a function that returns integer(kind=glcint) menuid.  The menuid
!    is the ID returned by glutCreateMenu.  You can either use the view_modifier
!    menu as your menu by calling glutAttachMenu immediately after
!    view_modifier_init, as in
!       menuid = view_modifier_init(lookfrom_x, lookfrom_y, lookfrom_z, &
!                                   lookat_x,   lookat_y,   lookat_z)
!       call glutAttachMenu(GLUT_RIGHT_BUTTON)
!    or by using the menuid to attach a submenu to your own menu, as in
!       call glutAddSubMenu("View Modifier",menuid)
!
! 4) in any callback functions that update the display, put
!       call reset_view
!    as the first executable statement
!
! Note that view_modifier_init sets the callback functions for glutMouseFunc,
! glutMotionFunc and glutSpecialFunc, so don't call these yourself
!
! The menu allows you to select what operation is attached to the left and
! middle mouse buttons and arrow keys, reset to the initial view, and quit.
! The right mouse button should be used for the menu.
!
! The public real(GLfloat) variable explode_factor contains a value
! in [1.0,infinity) which can be changed by moving the mouse forward or
! back while holding down a button selected for explode.  It is up to
! the application to figure out what to do with it.  If there are
! multiple windows, the same explode_factor applies to all.

use opengl_gl
use opengl_glu
use opengl_glut
implicit none
private
public :: view_modifier_init, reset_view, explode_factor
private :: ZOOM, PAN, ROTATE, SCALEX, SCALEY, SCALEZ, EXPLODE, LIGHT, RESET, &
           ABOVE, ABOVE_ORIGIN, QUIT, PI, &
           left_button_func, middle_button_func, arrow_key_func, &
           init_lookat, init_lookfrom, abszmax, init_lightpos, &
           init_xscale_factor, init_yscale_factor, init_zscale_factor, &
           angle, shift, xscale_factor, yscale_factor, zscale_factor, &
           lightpos, moving_left, moving_middle, begin_left, begin_middle, &
           cart2sphere, sphere2cart, cart3D_plus_cart3D, cart3D_minus_cart3D, &
           reset_to_init, mouse, motion, arrows, &
           menu_handler, set_left_button, set_middle_button, set_arrow_keys

! Maximum number of windows that view_modifier can handle simultaneously.
integer, parameter :: MAXWIN = 4

integer(kind=glcint), parameter :: ZOOM = 1, PAN = 2, ROTATE = 3, SCALEX = 4, &
                      SCALEY = 5, SCALEZ = 6, EXPLODE = 7, LIGHT = 8
integer(kind=glcint), parameter :: RESET = 10, ABOVE = 11, &
                      ABOVE_ORIGIN = 12, QUIT = 13
real(kind=gldouble), parameter :: PI = 3.141592653589793_gldouble

type, private :: cart2D ! 2D cartesian coordinates
   real(kind=gldouble) :: x, y
end type cart2D

type, private :: cart3D ! 3D cartesian coordinates
   real(kind=gldouble) :: x, y, z
end type cart3D

type, private :: sphere3D ! 3D spherical coordinates
   real(kind=gldouble) :: theta, phi, rho
end type sphere3D

type(cart2D), save, dimension(MAXWIN) :: angle
type(cart3D), save, dimension(MAXWIN) :: shift
real(kind=gldouble), save, dimension(MAXWIN) :: xscale_factor, yscale_factor, &
                                                zscale_factor, &
                                                init_xscale_factor, &
                                                init_yscale_factor, &
                                                init_zscale_factor
real(kind=glfloat), save :: explode_factor = 1.0_glfloat
logical, save :: moving_left, moving_middle
type(cart2D), save :: begin_left, begin_middle
integer(kind=glcint), save, dimension(MAXWIN) :: winids = -1
type(cart3D), save, dimension(MAXWIN) :: init_lookfrom, init_lookat
type(sphere3D), save, dimension(MAXWIN) :: init_lightpos, lightpos

interface operator(+)
   module procedure cart3D_plus_cart3D
end interface
interface operator(-)
   module procedure cart3D_minus_cart3D
end interface

! ------- Initial configuration -------

! Set the initial operation performed by each button and the arrow keys.
! The operations are ZOOM, PAN, ROTATE, SCALEX, SCALEY, SCALEZ, EXPLODE and
!   LIGHT

integer, save, dimension(MAXWIN) ::   left_button_func = ROTATE, &
                                    middle_button_func = ZOOM, &
                                        arrow_key_func = PAN

! Set an estimate of the maximum absolute value of z divided by the diameter
! of the (x,y) domain, for initializing the z scale factor

real(kind=gldouble), parameter :: &
   abszmax = 1.0_gldouble

! -------- end of Initial configuration ------

contains

!          -------------
subroutine reset_to_init(win)
!          -------------

integer, intent(in) :: win

! This resets the view to the initial configuration

type(sphere3D) :: slookfrom
type(cart3D) :: cpos
real(kind=glfloat) :: apos(4)

slookfrom = cart2sphere(init_lookfrom(win)-init_lookat(win))
angle(win)%x = -180.0_gldouble*slookfrom%theta/PI - 90.0_gldouble
angle(win)%y = -180.0_gldouble*slookfrom%phi/PI
shift(win)%x = 0.0_gldouble
shift(win)%y = 0.0_gldouble
shift(win)%z = -slookfrom%rho
xscale_factor(win) = init_xscale_factor(win)
yscale_factor(win) = init_yscale_factor(win)
zscale_factor(win) = init_zscale_factor(win)
lightpos(win) = init_lightpos(win)
cpos = sphere2cart(lightpos(win))
apos(1) = cpos%x
apos(2) = cpos%y
apos(3) = cpos%z
apos(4) = 0
call reset_view
call gllightfv(gl_light0, gl_position, apos)

call glutPostRedisplay

return
end subroutine reset_to_init

!          ---------------
subroutine view_from_above(win)
!          ---------------

integer, intent(in) :: win

! This sets the view to be from straight above

type(sphere3D) :: slookfrom
type(cart3D) :: cpos
real(kind=glfloat) :: apos(4)

slookfrom = cart2sphere(cart3D(0.0,0.0,1.0))
angle(win)%x = -180.0_gldouble*slookfrom%theta/PI
angle(win)%y = -180.0_gldouble*slookfrom%phi/PI
cpos = sphere2cart(lightpos(win))
apos(1) = cpos%x
apos(2) = cpos%y
apos(3) = cpos%z
apos(4) = 0
call reset_view
call gllightfv(gl_light0, gl_position, apos)

call glutPostRedisplay

return
end subroutine view_from_above

!          ----------------------
subroutine view_from_above_origin(win)
!          ----------------------

integer, intent(in) :: win

! This sets the view to be from straight above the origin

type(sphere3D) :: slookfrom
type(cart3D) :: cpos
real(kind=glfloat) :: apos(4)

slookfrom = cart2sphere(cart3D(0.0,0.0,1.0))
angle(win)%x = -180.0_gldouble*slookfrom%theta/PI
angle(win)%y = -180.0_gldouble*slookfrom%phi/PI
shift(win)%x = init_lookat(win)%x
shift(win)%y = init_lookat(win)%y
cpos = sphere2cart(lightpos(win))
apos(1) = cpos%x
apos(2) = cpos%y
apos(3) = cpos%z
apos(4) = 0
call reset_view
call gllightfv(gl_light0, gl_position, apos)

call glutPostRedisplay

return
end subroutine view_from_above_origin

!          ----------
subroutine reset_view
!          ----------

! This routine resets the view to the current orientation and scale

integer :: win

win = find_win()

call glMatrixMode(GL_MODELVIEW)
call glPopMatrix
call glPushMatrix
call glTranslated(shift(win)%x, shift(win)%y, shift(win)%z)
call glRotated(angle(win)%x, 0.0_gldouble, 0.0_gldouble, 1.0_gldouble)
call glRotated(angle(win)%y, cos(PI*angle(win)%x/180.0_gldouble), &
               -sin(PI*angle(win)%x/180.0_gldouble), 0.0_gldouble)
call glTranslated(-init_lookat(win)%x, -init_lookat(win)%y, -init_lookat(win)%z)
call glScaled(xscale_factor(win),yscale_factor(win),zscale_factor(win))

return
end subroutine reset_view

!          -----
subroutine mouse(button, state, x, y)
!          -----
integer(kind=glcint), intent(in out) :: button, state, x, y

! This gets called when a mouse button changes
 
  if (button == GLUT_LEFT_BUTTON .and. state == GLUT_DOWN) then
    moving_left = .true.
    begin_left = cart2D(x,y)
  endif
  if (button == GLUT_LEFT_BUTTON .and. state == GLUT_UP) then
    moving_left = .false.
  endif
  if (button == GLUT_MIDDLE_BUTTON .and. state == GLUT_DOWN) then
    moving_middle = .true.
    begin_middle = cart2D(x,y)
  endif
  if (button == GLUT_MIDDLE_BUTTON .and. state == GLUT_UP) then
    moving_middle = .false.
  endif
end subroutine mouse

!          ------
subroutine motion(x, y)
!          ------
integer(kind=glcint), intent(in out) :: x, y

! This gets called when the mouse moves

integer :: button_function, win
type(cart2D) :: begin
type(cart3D) :: cpos
real(kind=glfloat) :: apos(4)
real(kind=gldouble) :: factor

win = find_win()

! Determine and apply the button function

if (moving_left) then
   button_function = left_button_func(win)
   begin = begin_left
else if(moving_middle) then
   button_function = middle_button_func(win)
   begin = begin_middle
end if

select case(button_function)
case (ZOOM)
   if (y < begin%y) then
      factor = 1.0_gldouble/(1.0_gldouble + .002_gldouble*(begin%y-y))
   else if (y > begin%y) then
      factor = 1.0_gldouble + .002_gldouble*(y-begin%y)
   else
      factor = 1.0_gldouble
   end if
   shift(win)%z = factor*shift(win)%z
case (PAN)
   shift(win)%x = shift(win)%x - .001*shift(win)%z*(x - begin%x)
   shift(win)%y = shift(win)%y + .001*shift(win)%z*(y - begin%y)
case (ROTATE)
   angle(win)%x = angle(win)%x + (x - begin%x)
   angle(win)%y = angle(win)%y + (y - begin%y)
   cpos = sphere2cart(lightpos(win))
   apos(1) = cpos%x
   apos(2) = cpos%y
   apos(3) = cpos%z
   apos(4) = 0
   call gllightfv(gl_light0, gl_position, apos)
case (SCALEX)
   if (y < begin%y) then
      factor = 1.0_gldouble + .002_gldouble*(begin%y-y)
   else if (y > begin%y) then
      factor = 1.0_gldouble/(1.0_gldouble + .002_gldouble*(y-begin%y))
   else
      factor = 1.0_gldouble
   end if
   xscale_factor(win) = xscale_factor(win) * factor
case (SCALEY)
   if (y < begin%y) then
      factor = 1.0_gldouble + .002_gldouble*(begin%y-y)
   else if (y > begin%y) then
      factor = 1.0_gldouble/(1.0_gldouble + .002_gldouble*(y-begin%y))
   else
      factor = 1.0_gldouble
   end if
   yscale_factor(win) = yscale_factor(win) * factor
case (SCALEZ)
   if (y < begin%y) then
      factor = 1.0_gldouble + .002_gldouble*(begin%y-y)
   else if (y > begin%y) then
      factor = 1.0_gldouble/(1.0_gldouble + .002_gldouble*(y-begin%y))
   else
      factor = 1.0_gldouble
   end if
   zscale_factor(win) = zscale_factor(win) * factor
case (EXPLODE)
   if (y < begin%y) then
      factor = 1.0_gldouble + .002_gldouble*(begin%y-y)
   else if (y > begin%y) then
      factor = 1.0_gldouble/(1.0_gldouble + .002_gldouble*(y-begin%y))
   else
      factor = 1.0_gldouble
   end if
   explode_factor = explode_factor * factor
   if (explode_factor < 1.0_glfloat) explode_factor = 1.0_glfloat
case (LIGHT)
   lightpos(win)%theta = lightpos(win)%theta + .01_gldouble*(x - begin%x)
   lightpos(win)%phi = lightpos(win)%phi + .01_gldouble*(y - begin%y)
   cpos = sphere2cart(lightpos(win))
   apos(1) = cpos%x
   apos(2) = cpos%y
   apos(3) = cpos%z
   apos(4) = 0
   call gllightfv(gl_light0, gl_position, apos)
end select

! update private variables and redisplay

if (moving_left) then
   begin_left = cart2D(x,y)
else if(moving_middle) then
   begin_middle = cart2D(x,y)
endif

if (moving_left .or. moving_middle) then
   call glutPostRedisplay
endif

return
end subroutine motion

!          ------
subroutine arrows(key, x, y)
!          ------
integer(glcint), intent(in out) :: key, x, y

! This routine handles the arrow key operations

integer :: win
real(kind=gldouble) :: factor
type(cart3D) :: cpos
real(kind=glfloat) :: apos(4)

win = find_win()

select case(arrow_key_func(win))
case(ZOOM)
   select case(key)
   case(GLUT_KEY_DOWN)
      factor = 1.0_gldouble + .02_gldouble
   case(GLUT_KEY_UP)
      factor = 1.0_gldouble/(1.0_gldouble + .02_gldouble)
   case default
      factor = 1.0_gldouble
   end select
   shift(win)%z = factor*shift(win)%z
case(PAN)
   select case(key)
   case(GLUT_KEY_LEFT)
      shift(win)%x = shift(win)%x + .002*shift(win)%z
   case(GLUT_KEY_RIGHT)
      shift(win)%x = shift(win)%x - .002*shift(win)%z
   case(GLUT_KEY_DOWN)
      shift(win)%y = shift(win)%y + .002*shift(win)%z
   case(GLUT_KEY_UP)
      shift(win)%y = shift(win)%y - .002*shift(win)%z
   end select
case(ROTATE)
   select case(key)
   case(GLUT_KEY_LEFT)
      angle(win)%x = angle(win)%x - 1.0_gldouble
   case(GLUT_KEY_RIGHT)
      angle(win)%x = angle(win)%x + 1.0_gldouble
   case(GLUT_KEY_DOWN)
      angle(win)%y = angle(win)%y - 1.0_gldouble
   case(GLUT_KEY_UP)
      angle(win)%y = angle(win)%y + 1.0_gldouble
   end select
   cpos = sphere2cart(lightpos(win))
   apos(1) = cpos%x
   apos(2) = cpos%y
   apos(3) = cpos%z
   apos(4) = 0
   call gllightfv(gl_light0, gl_position, apos)
case(SCALEX)
   select case(key)
   case(GLUT_KEY_DOWN)
      factor = 1.0_gldouble/(1.0_gldouble + .02_gldouble)
   case(GLUT_KEY_UP)
      factor = 1.0_gldouble + .02_gldouble
   case default
      factor = 1.0_gldouble
   end select
   xscale_factor(win) = xscale_factor(win) * factor
case(SCALEY)
   select case(key)
   case(GLUT_KEY_DOWN)
      factor = 1.0_gldouble/(1.0_gldouble + .02_gldouble)
   case(GLUT_KEY_UP)
      factor = 1.0_gldouble + .02_gldouble
   case default
      factor = 1.0_gldouble
   end select
   yscale_factor(win) = yscale_factor(win) * factor
case(SCALEZ)
   select case(key)
   case(GLUT_KEY_DOWN)
      factor = 1.0_gldouble/(1.0_gldouble + .02_gldouble)
   case(GLUT_KEY_UP)
      factor = 1.0_gldouble + .02_gldouble
   case default
      factor = 1.0_gldouble
   end select
   zscale_factor(win) = zscale_factor(win) * factor
case(EXPLODE)
   select case(key)
   case(GLUT_KEY_DOWN)
      factor = 1.0_gldouble/(1.0_gldouble + .02_gldouble)
   case(GLUT_KEY_UP)
      factor = 1.0_gldouble + .02_gldouble
   case default
      factor = 1.0_gldouble
   end select
   explode_factor = explode_factor * factor
   if (explode_factor < 1.0_glfloat) explode_factor = 1.0_glfloat
case (LIGHT)
   select case(key)
   case(GLUT_KEY_LEFT)
      lightpos(win)%theta = lightpos(win)%theta - .01_gldouble
   case(GLUT_KEY_RIGHT)
      lightpos(win)%theta = lightpos(win)%theta + .01_gldouble
   case(GLUT_KEY_DOWN)
      lightpos(win)%phi = lightpos(win)%phi - .01_gldouble
   case(GLUT_KEY_UP)
      lightpos(win)%phi = lightpos(win)%phi + .01_gldouble
   end select
   cpos = sphere2cart(lightpos(win))
   apos(1) = cpos%x
   apos(2) = cpos%y
   apos(3) = cpos%z
   apos(4) = 0
   call gllightfv(gl_light0, gl_position, apos)

end select
   
call glutPostRedisplay

return
end subroutine arrows

!          ------------
subroutine menu_handler(value)
!          ------------
integer(kind=glcint), intent(in out) :: value

! This routine handles the first level entries in the menu

integer :: win

win = find_win()

select case(value)

case(RESET)
   call reset_to_init(win)
case(ABOVE)
   call view_from_above(win)
case(ABOVE_ORIGIN)
   call view_from_above_origin(win)
case(QUIT)
   stop

end select

return
end subroutine menu_handler

!          ---------------
subroutine set_left_button(value)
!          ---------------
integer(kind=glcint), intent(in out) :: value

! This routine sets the function of the left button as given by menu selection

integer :: win

win = find_win()
left_button_func(win) = value

return
end subroutine set_left_button

!          -----------------
subroutine set_middle_button(value)
!          -----------------
integer(kind=glcint), intent(in out) :: value

! This routine sets the function of the middle button as given by menu selection

integer :: win

win = find_win()
middle_button_func(win) = value

return
end subroutine set_middle_button

!          --------------
subroutine set_arrow_keys(value)
!          --------------
integer(kind=glcint), intent(in out) :: value

! This routine sets the function of the arrow keys as given by menu selection

integer :: win

win = find_win()
arrow_key_func(win) = value

return
end subroutine set_arrow_keys

!        ------------------
function view_modifier_init(lookfrom_x, lookfrom_y, lookfrom_z, &
                            lookat_x,   lookat_y,   lookat_z) result(menuid)
!        ------------------
real(kind=gldouble), intent(in) :: lookfrom_x, lookfrom_y, lookfrom_z, &
                                   lookat_x,   lookat_y,   lookat_z
integer(kind=glcint) :: menuid

! This initializes the view modifier variables and sets initial view.
! It should be called immediately after glutCreateWindow

integer(kind=glcint) :: button_left, button_middle, arrow_keys
integer :: win

! find an available window number

do win = 1,MAXWIN
   if (winids(win) == -1) exit
end do

if (win > MAXWIN) then
   print *,"ERROR in view_modifier: maximum number of windows exceeded"
   menuid = -1
   return
endif

winids(win) = glutGetWindow()

! set the initial view parameters

init_lookfrom(win) = cart3D(lookfrom_x, lookfrom_y, lookfrom_z)
init_lookat(win) = cart3D(lookat_x, lookat_y, lookat_z)
init_xscale_factor(win) = 1.0_gldouble
init_yscale_factor(win) = 1.0_gldouble
init_zscale_factor(win) = 0.5_gldouble/abszmax
init_lightpos(win) = sphere3D(2.91_gldouble,2.32_gldouble,10.0_gldouble)

! set the callback functions

call glutMouseFunc(mouse)
call glutMotionFunc(motion)
call glutSpecialFunc(arrows)

! create the menu

button_left = glutCreateMenu(set_left_button)
call glutAddMenuEntry("rotate",ROTATE)
call glutAddMenuEntry("zoom",ZOOM)
call glutAddMenuEntry("pan",PAN)
call glutAddMenuEntry("scale x",SCALEX)
call glutAddMenuEntry("scale y",SCALEY)
call glutAddMenuEntry("scale z",SCALEZ)
call glutAddMenuEntry("explode",EXPLODE)
call glutAddMenuEntry("move light",LIGHT)
button_middle = glutCreateMenu(set_middle_button)
call glutAddMenuEntry("rotate",ROTATE)
call glutAddMenuEntry("zoom",ZOOM)
call glutAddMenuEntry("pan",PAN)
call glutAddMenuEntry("scale x",SCALEX)
call glutAddMenuEntry("scale y",SCALEY)
call glutAddMenuEntry("scale z",SCALEZ)
call glutAddMenuEntry("explode",EXPLODE)
call glutAddMenuEntry("move light",LIGHT)
arrow_keys = glutCreateMenu(set_arrow_keys)
call glutAddMenuEntry("rotate",ROTATE)
call glutAddMenuEntry("zoom",ZOOM)
call glutAddMenuEntry("pan",PAN)
call glutAddMenuEntry("scale x",SCALEX)
call glutAddMenuEntry("scale y",SCALEY)
call glutAddMenuEntry("scale z",SCALEZ)
call glutAddMenuEntry("explode",EXPLODE)
call glutAddMenuEntry("move light",LIGHT)
menuid = glutCreateMenu(menu_handler)
call glutAddSubMenu("left mouse button",button_left)
call glutAddSubMenu("middle mouse button",button_middle)
call glutAddSubMenu("arrow keys",arrow_keys)
call glutAddMenuEntry("reset to initial view",RESET)
call glutAddMenuEntry("view from above",ABOVE)
call glutAddMenuEntry("view from above origin",ABOVE_ORIGIN)
call glutAddMenuEntry("quit",QUIT)

! set the perspective
! smaller values of zNear (third argument) allow zooming in further, and larger
! values of zFar (fourth argument) allowing zooming out further.  However,
! a large ratio affects the depth buffer precision.  Roughly log_2(zFar/zNear)
! bits are lost.

call glMatrixMode(GL_PROJECTION)
!call gluPerspective(10.0_gldouble, 1.0_gldouble, 0.1_gldouble, 200.0_gldouble)
!call gluPerspective(10.0_gldouble, 1.0_gldouble, 0.0001_gldouble, 1000.0_gldouble)
call gluPerspective(10.0_gldouble, 1.0_gldouble, 0.05_gldouble, 500.0_gldouble)

! set the initial view

call glPushMatrix
call reset_to_init(win)

return
end function view_modifier_init

!        --------
function find_win() result(win)
!        --------
integer :: win

! Determines the window number of the current window.  Returns -1 if
! if couldn't find it.

integer(kind=glcint) :: winid

winid = glutGetWindow()

do win = 1,MAXWIN
   if (winids(win) == winid) exit
end do

if (win > MAXWIN) then
   print *,"ERROR in view_modifier: current window is not registered"
   win = -1
endif

return
end function find_win

!        -----------
function sphere2cart(spoint) result(cpoint)
!        -----------
type(sphere3D), intent(in) :: spoint
type(cart3D) :: cpoint

! This converts a 3D point from spherical to cartesean coordinates

real(kind=gldouble) :: t,p,r

t=spoint%theta
p=spoint%phi
r=spoint%rho

cpoint%x = r*cos(t)*sin(p)
cpoint%y = r*sin(t)*sin(p)
cpoint%z = r*cos(p)

return
end function sphere2cart

!        -----------
function cart2sphere(cpoint) result(spoint)
!        -----------
type(cart3D), intent(in) :: cpoint
type(sphere3D) :: spoint

! This converts a 3D point from cartesean to spherical coordinates

real(kind=gldouble) :: x,y,z

x=cpoint%x
y=cpoint%y
z=cpoint%z

spoint%rho = sqrt(x*x+y*y+z*z)
if (x==0.0_gldouble .and. y==0.0_gldouble) then
   spoint%theta = 0.0_gldouble
else
   spoint%theta = atan2(y,x)
end if
if (spoint%rho == 0.0_gldouble) then
   spoint%phi = 0.0_gldouble
else
   spoint%phi = acos(z/spoint%rho)
endif

return
end function cart2sphere

!        ------------------
function cart3D_plus_cart3D(cart1,cart2) result(cart3)
!        ------------------
type(cart3D), intent(in) :: cart1, cart2
type(cart3D) :: cart3

! Compute the sum of two 3D cartesean points

cart3%x = cart1%x + cart2%x
cart3%y = cart1%y + cart2%y
cart3%z = cart1%z + cart2%z

return
end function cart3D_plus_cart3D

!        -------------------
function cart3D_minus_cart3D(cart1,cart2) result(cart3)
!        -------------------
type(cart3D), intent(in) :: cart1, cart2
type(cart3D) :: cart3

! Compute the difference of two 3D cartesean points

cart3%x = cart1%x - cart2%x
cart3%y = cart1%y - cart2%y
cart3%z = cart1%z - cart2%z

return
end function cart3D_minus_cart3D

end module view_modifier
