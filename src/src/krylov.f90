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

module krylov

!----------------------------------------------------
! This module contains Krylov space iterative solvers.
!
! communication tags in this module are of the form 30xx
!----------------------------------------------------

!----------------------------------------------------
! Other modules used are:

use global
use message_passing
use stopwatch
use gridtype_mod
use linsystype_mod
use linsys_util
use hbmg
use error_estimators
!----------------------------------------------------

implicit none
private
public phaml_cg, phaml_gmres

!----------------------------------------------------
! Symbolic constants are:

integer, parameter :: DO_MATVEC = 1, &
                      DO_PRECON = 2, &
                      QUIT = 3
!----------------------------------------------------
! Module variables are:

type(linsys_type), save, pointer :: loc_linsys
type(proc_info), save, pointer :: loc_procs
type(grid_type), save, pointer :: loc_grid
type(io_options), save, pointer :: loc_io_cntl
type(solver_options), save, pointer :: loc_solver_cntl
logical, save, pointer :: loc_still_sequential

type arrays
   real(my_real), pointer :: val(:)
   integer :: n
end type arrays

type(arrays), save, allocatable :: each_rhs(:)

integer, save :: myneq, neq
real(my_real), save, allocatable :: tempvec(:)

!----------------------------------------------------
! Interfaces are

interface

   function sdot(n,x,incx,y,incy)
   use global
   integer :: n, incx, incy
   real(my_real) :: x(n), y(n)
   real(my_real) :: sdot
   end function sdot

   function ddot(n,x,incx,y,incy)
   use global
   integer :: n, incx, incy
   real(my_real) :: x(n), y(n)
   real(my_real) :: ddot
   end function ddot

end interface

contains

!          -----------
subroutine phaml_gmres(linsys,procs,grid,io_cntl,solver_cntl,still_sequential)
!          -----------

!----------------------------------------------------
! This routine routine solves the linear system using GMRES.
! Processor 1 calls the Templates GMRES subroutine.  The other
! processors call a routine that waits for processor 1 to request help with
! matvecs, etc.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(inout), target :: linsys
type(proc_info), intent(in), target :: procs
type(grid_type), intent(inout), target :: grid
type(io_options), intent(in), target :: io_cntl
type(solver_options), intent(in), target :: solver_cntl
logical, intent(in), target :: still_sequential
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

! set module pointers to variables that are needed in the matvec and
! preconditioner routines

loc_linsys => linsys
loc_procs => procs
loc_grid => grid
loc_io_cntl => io_cntl
loc_solver_cntl => solver_cntl
loc_still_sequential => still_sequential

! processor 1 calls a routine that gathers information, calls the Templates
! CG subroutine and scatters the result

if (my_proc(loc_procs) == 1) then

   call phaml_gmres1

! master calls multigrid when the others do

elseif (my_proc(loc_procs) == MASTER) then

   call master_waits

! other processors call a routine that helps with matvecs, etc.

else ! not processor 1

   call help_solve

endif

end subroutine phaml_gmres

!          ------------
subroutine phaml_gmres1
!          ------------

!----------------------------------------------------
! This routine is called by processor 1 to gather information from the other
! processors, call the Templates GMRES routine, and scatter the result.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

!----------------------------------------------------
! Local variables:

real(my_real), allocatable :: rhs(:), solution(:), holdsoln(:)
real(kind(0.0d0)), allocatable :: work(:), work2(:)
real(my_real), pointer :: rrecv(:)
real(my_real) :: resid
integer, pointer :: irecv(:)
integer :: nproc, astat, icount, p, i, j, ni, nr, iter, info, restart, &
           objtype, brank, srank, lid
!----------------------------------------------------
! Begin executable code

! Only supporting double precision at this time

if (my_real /= kind(1.0d0)) then
   ierr = UNCLASSIFIED_ERROR
   call fatal("GMRES_SOLVER currently only supports double precision", &
              "Change definition of my_real in global.f90")
   stop
endif

! Determine the total number of equations

myneq = count(loc_linsys%iown(1:loc_linsys%neq) .and. &
              loc_linsys%equation_type(1:loc_linsys%neq) /= DIRICHLET)
if (loc_still_sequential) then
   neq = myneq
else
   call start_watch((/cpsolve,ctsolve/))
   neq = phaml_global_sum(loc_procs,myneq,3001,.true.)
   call stop_watch((/cpsolve,ctsolve/))
endif

! some constants and workspace

nproc = num_proc(loc_procs)
allocate(tempvec(loc_linsys%neq),holdsoln(loc_linsys%neq),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in phaml_gmres1",procs=loc_procs)
   stop
endif
iter = loc_solver_cntl%krylov_iter
restart = max(1,min(neq-2,loc_solver_cntl%krylov_restart))
resid = loc_solver_cntl%krylov_tol
if (resid == KRYLOV_ERREST_TOL) then
! TEMP for systems of equations, should be a combination not the max
   call error_estimate(loc_grid,loc_procs,errest_L2=resid)
   if (.not. loc_still_sequential) then
      call start_watch((/cpsolve,ctsolve/))
      resid = sqrt(phaml_global_sum(loc_procs,resid**2,3005,.true.))
      call stop_watch((/cpsolve,ctsolve/))
   endif
   resid = max(100*epsilon(1.0_my_real),resid/100)
endif

! save the current solution to restore unowned parts later

holdsoln = loc_linsys%solution(1:loc_linsys%neq)

! Start the right hand side with the entries from this processor, subtracting
! off Dirichlet contributions

allocate(rhs(neq),solution(neq),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in phaml_gmres1",procs=loc_procs)
   stop
endif

icount = 0
do i=1,loc_linsys%neq
   if (loc_linsys%iown(i) .and. loc_linsys%equation_type(i) /= DIRICHLET) then
      icount = icount + 1
      rhs(icount) = loc_linsys%rhs(i)
      do j=loc_linsys%begin_row(i),loc_linsys%end_row(i)
         if (loc_linsys%column_index(j) == NO_ENTRY) cycle
         if (loc_linsys%equation_type(loc_linsys%column_index(j)) == DIRICHLET) then
            rhs(icount) = rhs(icount) - &
               loc_linsys%matrix_val(j)*loc_linsys%solution(loc_linsys%column_index(j))
         endif
      end do
   endif
end do

! Gather the rest of the right hand side from the other processors

if (.not. loc_still_sequential) then

   allocate(each_rhs(nproc),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in phaml_gmres1",procs=loc_procs)
      stop
   endif

   call start_watch((/cpsolve,ctsolve/))
   do i=2,nproc
      call phaml_recv(loc_procs,p,irecv,ni,rrecv,nr,3002)
      each_rhs(p)%n = nr
      each_rhs(p)%val => rrecv
      if (associated(irecv)) deallocate(irecv,stat=astat)
   end do
   call stop_watch((/cpsolve,ctsolve/))

   do i=2,nproc
      if (associated(each_rhs(i)%val)) then
         rhs(icount+1:icount+each_rhs(i)%n) = each_rhs(i)%val
         icount = icount + each_rhs(i)%n
         deallocate(each_rhs(i)%val,stat=astat)
      endif
   end do

   if (icount /= neq) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("count and neq are not the same in phaml_gmres1", &
                 intlist=(/icount,neq/),procs=loc_procs)
      stop
   endif

endif ! not still_sequential

! Set initial guess to be existing solution.  Start with entries on this proc.

solution = 0.0_my_real
icount = 0
do i=1,loc_linsys%neq
   if (loc_linsys%iown(i) .and. loc_linsys%equation_type(i) /= DIRICHLET) then
      icount = icount + 1
      solution(icount) = loc_linsys%solution(i)
   endif
end do

! Gather the rest of the solution from the other processors

if (.not. loc_still_sequential) then

   call start_watch((/cpsolve,ctsolve/))
   do i=2,nproc
      call phaml_recv(loc_procs,p,irecv,ni,rrecv,nr,3003)
      each_rhs(p)%val => rrecv
      if (associated(irecv)) deallocate(irecv,stat=astat)
   end do
   call stop_watch((/cpsolve,ctsolve/))

   do i=2,nproc
      if (associated(each_rhs(i)%val)) then
         solution(icount+1:icount+each_rhs(i)%n) = each_rhs(i)%val
         icount = icount + each_rhs(i)%n
         deallocate(each_rhs(i)%val,stat=astat)
      endif
   end do

endif ! not still_sequential

! workspace

allocate(work((restart+6)*neq),work2(2*restart**2+2),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in phaml_gmres1",procs=loc_procs)
   stop
endif
work = 0
work2 = 0

! Call the Templates GMRES routine

call gmres(neq,rhs,solution,restart,work,neq,work2,restart,iter,resid, &
        matvec,precond,info, &
        loc_io_cntl%print_error_when==FREQUENTLY.or.loc_io_cntl%print_error_when==TOO_MUCH)

if (info > 0) then
   call warning("GMRES did not converge to tolerance",reallist=(/resid/),&
                intlist=(/info/))
elseif (info < 0) then
   ierr = UNCLASSIFIED_ERROR
   call fatal("GMRES returned error code",intlist=(/info/),procs=loc_procs)
   stop
endif

! Scatter the solution to the other processors

call start_watch((/cpsolve,ctsolve/))
if (loc_still_sequential) then
   do p=2,nproc
      call phaml_send(loc_procs,p,(/QUIT/),1,solution,neq,3000)
   end do
else
   icount = myneq
   do p=2,nproc
      call phaml_send(loc_procs,p,(/QUIT/),1, &
                     solution(icount+1:icount+each_rhs(p)%n),each_rhs(p)%n,3000)
      icount = icount + each_rhs(p)%n
   end do
endif
call phaml_send(loc_procs,MASTER,(/QUIT/),1,(/0.0_my_real/),0,3000)
call stop_watch((/cpsolve,ctsolve/))

! Copy this processor's part of the solution, getting Dirichlet points from grid

icount = 0
do i=1,loc_linsys%neq
   if (loc_linsys%iown(i) .and. loc_linsys%equation_type(i) /= DIRICHLET) then
      icount = icount + 1
      loc_linsys%solution(i) = solution(icount)
   elseif (loc_linsys%equation_type(i)==DIRICHLET) then
      call eq_to_grid(loc_linsys,loc_linsys%gid(i),objtype,brank, &
                      srank,lid,loc_grid)
      select case(objtype)
      case (VERTEX_ID)
         loc_linsys%solution(i) = loc_grid%vertex_solution(lid,srank,1)
      case (EDGE_ID)
         loc_linsys%solution(i) = loc_grid%edge(lid)%solution(brank,srank,1)
      end select
   endif
end do

! restore previous solution at unowned equations.  This is to get the edge
! solutions that the owner doesn't have because it's a different level

do i=1,loc_linsys%neq
   if (.not. loc_linsys%iown(i)) loc_linsys%solution(i) = holdsoln(i)
end do

! exchange solution at unowned equations

if (.not. loc_still_sequential) then
   call exchange_fudop_vect(loc_linsys%solution(1:),loc_procs,loc_linsys,3010,&
                            3011,3012)
endif

deallocate(tempvec,holdsoln,rhs,solution,work,work2,stat=astat)
if (.not. loc_still_sequential) deallocate(each_rhs,stat=astat)

end subroutine phaml_gmres1

!          ----------
subroutine help_solve
!          ----------

!----------------------------------------------------
! This routine is called by all processors except processor 1 to help
! with matvecs, etc.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

!----------------------------------------------------
! Local variables:

real(my_real), allocatable :: rhs(:), holdsoln(:)
real(my_real), pointer :: rrecv(:)
real(my_real) :: resid
integer, pointer :: irecv(:)
integer :: neq, dum, astat, i, j, icount, ni, nr, p, lev, objtype, brank, &
           srank, lid
type(io_options) :: io_cntl_noprint
!----------------------------------------------------
! Begin executable code

io_cntl_noprint = loc_io_cntl
io_cntl_noprint%print_error_when = NEVER

allocate(tempvec(loc_linsys%neq),holdsoln(loc_linsys%neq),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in help_solve",procs=loc_procs)
   stop
endif

! save the current solution to restore unowned parts later

holdsoln = loc_linsys%solution(1:loc_linsys%neq)

if (.not. loc_still_sequential) then

! Determine the total number of equations

   neq = count(loc_linsys%iown(1:loc_linsys%neq) .and. &
               loc_linsys%equation_type(1:loc_linsys%neq) /= DIRICHLET)
   call start_watch((/cpsolve,ctsolve/))
   dum = phaml_global_sum(loc_procs,neq,3001,.true.)
   call stop_watch((/cpsolve,ctsolve/))

! help with error estimate

   if (loc_solver_cntl%krylov_tol == KRYLOV_ERREST_TOL) then
! TEMP for systems of equations, should be a combination not the max
      call error_estimate(loc_grid,loc_procs,errest_L2=resid)
      call start_watch((/cpsolve,ctsolve/))
      resid = sqrt(phaml_global_sum(loc_procs,resid**2,3005,.true.))
      call stop_watch((/cpsolve,ctsolve/))
   endif

! Send this processor's part of the right hand side, with Dirichlet
! contributions removed, to processor 1

   allocate(rhs(neq),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in help_solve",procs=loc_procs)
      stop
   endif

   icount = 0
   do i=1,loc_linsys%neq
      if (loc_linsys%iown(i) .and. loc_linsys%equation_type(i) /= DIRICHLET) then
         icount = icount + 1
         rhs(icount) = loc_linsys%rhs(i)
         do j=loc_linsys%begin_row(i),loc_linsys%end_row(i)
            if (loc_linsys%column_index(j) == NO_ENTRY) cycle
            if (loc_linsys%equation_type(loc_linsys%column_index(j)) == DIRICHLET) then
               rhs(icount) = rhs(icount) - &
                  loc_linsys%matrix_val(j)*loc_linsys%solution(loc_linsys%column_index(j))
            endif
         end do
      endif
   end do

   if (icount /= neq) then
      ierr = PHAML_INTERNAL_ERROR
      call fatal("count and neq are not the same in help_solve", &
                 intlist=(/icount,neq/),procs=loc_procs)
      stop
   endif

   call start_watch((/cpsolve,ctsolve/))
   call phaml_send(loc_procs,1,(/0/),0,rhs,neq,3002)
   call stop_watch((/cpsolve,ctsolve/))

! Send this processor's part of the solution for the initial guess; use rhs
! for space

   icount = 0
   do i=1,loc_linsys%neq
      if (loc_linsys%iown(i) .and. loc_linsys%equation_type(i) /= DIRICHLET) then
         icount = icount + 1
         rhs(icount) = loc_linsys%solution(i)
      endif
   end do

   call start_watch((/cpsolve,ctsolve/))
   call phaml_send(loc_procs,1,(/0/),0,rhs,neq,3003)
   call stop_watch((/cpsolve,ctsolve/))

endif ! not still_sequential

! go into a loop to process requests from processor 1.  It it exited when
! the quit request arrives

do

! do not count this receive as communication time
   call phaml_recv(loc_procs,p,irecv,ni,rrecv,nr,3000)

   select case(irecv(1))

   case (DO_MATVEC)

! copy the vector to a temporary

      icount = 0
      tempvec = 0.0_my_real
      do i=1,loc_linsys%neq
         if (loc_linsys%iown(i) .and. loc_linsys%equation_type(i) /= DIRICHLET) then
            icount = icount + 1
            tempvec(i) = rrecv(icount)
         endif
      end do

      deallocate(irecv,stat=astat)
      if (associated(rrecv)) deallocate(rrecv,stat=astat)

! perform the matvec with the result in loc_linsys%solution

      call matrix_times_vector(tempvec,loc_linsys%solution(1:loc_linsys%neq), &
                               loc_linsys,loc_procs,loc_still_sequential, &
                               3021,3022,3023,3024,3025,3026,nodirch=.true., &
                               nocomm2=.true.)

! return this processor's part of the product to processor 1; use rhs for space

      icount = 0
      do i=1,loc_linsys%neq
         if (loc_linsys%iown(i) .and. loc_linsys%equation_type(i) /= DIRICHLET) then
            icount = icount + 1
            rhs(icount) = loc_linsys%solution(i)
         endif
      end do

      call start_watch((/cpsolve,ctsolve/))
      call phaml_send(loc_procs,1,(/0/),0,rhs,neq,3030)
      call stop_watch((/cpsolve,ctsolve/))

   case (DO_PRECON)

      select case(loc_solver_cntl%preconditioner)

      case (NO_PRECONDITION)

         ierr = PHAML_INTERNAL_ERROR
         call fatal("Processor other than 1 received NO_PRECONDITION in help_solve")
         stop

      case (MG_PRECONDITION)

! preserve the right hand side

         tempvec = loc_linsys%rhs(1:loc_linsys%neq)

! copy the vector to the right hand side

         icount = 0
         loc_linsys%rhs(1:loc_linsys%neq) = 0.0_my_real
         do i=1,loc_linsys%neq
            if (loc_linsys%iown(i) .and. loc_linsys%equation_type(i) /= DIRICHLET) then
               icount = icount + 1
               loc_linsys%rhs(i) = rrecv(icount)
            endif
         end do

! get unowned components of right hand side

         if (.not. loc_still_sequential) then
            do lev=loc_linsys%nlev+1,2,-1
               call basis_change(lev,TO_HIER,loc_linsys)
            end do
            call sum_fudop_vect(loc_linsys%rhs(1:loc_linsys%neq),loc_procs, &
                                loc_linsys,3041,3042,3043)
            do lev=2,loc_linsys%nlev+1
               call basis_change(lev,TO_NODAL,loc_linsys)
            end do
         endif

! apply multigrid

         loc_linsys%solution(1:loc_linsys%neq) = 0.0_my_real
         call multigrid(loc_grid,loc_procs,loc_linsys,io_cntl_noprint, &
                        loc_solver_cntl,loc_still_sequential)

! return this processor's part of the preconditioned vector to processor 1;
! use rhs for space

         if (.not. loc_still_sequential) then
            icount = 0
            do i=1,loc_linsys%neq
               if (loc_linsys%iown(i) .and. loc_linsys%equation_type(i) /= DIRICHLET) then
                  icount = icount + 1
                  rhs(icount) = loc_linsys%solution(i)
               endif
            end do

            call start_watch((/cpsolve,ctsolve/))
            call phaml_send(loc_procs,1,(/0/),0,rhs,neq,3040)
            call stop_watch((/cpsolve,ctsolve/))
         endif

! restore the original right hand side

         loc_linsys%rhs(1:loc_linsys%neq) = tempvec

      case default

         ierr = USER_INPUT_ERROR
         call fatal("Illegal choice of preconditioner for CG_SOLVER or GMRES_SOLVER")
         stop

      end select ! preconditioner

      deallocate(irecv,stat=astat)
      if (associated(rrecv)) deallocate(rrecv,stat=astat)

   case (QUIT)

! copy solution to linsys, getting Dirichlet points from grid

      icount = 0
      do i=1,loc_linsys%neq
         if ((loc_still_sequential .or. loc_linsys%iown(i)) .and. &
              loc_linsys%equation_type(i) /= DIRICHLET) then
            icount = icount + 1
            loc_linsys%solution(i) = rrecv(icount)
         elseif (loc_linsys%equation_type(i)==DIRICHLET) then
            call eq_to_grid(loc_linsys,loc_linsys%gid(i),objtype,brank, &
                            srank,lid,loc_grid)
            select case(objtype)
            case (VERTEX_ID)
               loc_linsys%solution(i) = loc_grid%vertex_solution(lid,srank,1)
            case (EDGE_ID)
               loc_linsys%solution(i) = loc_grid%edge(lid)%solution(brank,srank,1)
            end select
         endif
      end do

      deallocate(irecv,stat=astat)
      if (associated(rrecv)) deallocate(rrecv,stat=astat)

! restore previous solution at unowned equations.  This is to get the edge
! solutions that the owner doesn't have because it's a different level

      do i=1,loc_linsys%neq
         if (.not. loc_linsys%iown(i)) loc_linsys%solution(i) = holdsoln(i)
      end do

! exchange solution at unowned equations

      if (.not. loc_still_sequential) then
         call exchange_fudop_vect(loc_linsys%solution(1:),loc_procs, &
                                  loc_linsys,3010,3011,3012)
      endif

      exit

   end select

end do

deallocate(tempvec,holdsoln,stat=astat)
if (allocated(rhs)) deallocate(rhs,stat=astat)

end subroutine help_solve

!          ------------
subroutine master_waits
!          ------------

!----------------------------------------------------
! This routine is called by the master process.  The master needs to be
! involved in calls to multgrid, but other than that it is just waiting.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

!----------------------------------------------------
! Local variables:

integer :: p,ni,nr,astat
integer, pointer :: irecv(:)
real(my_real), pointer :: rrecv(:)
type(io_options) :: io_cntl_noprint
!----------------------------------------------------
! Begin executable code

io_cntl_noprint = loc_io_cntl
io_cntl_noprint%print_error_when = NEVER

do
   call phaml_recv(loc_procs,p,irecv,ni,rrecv,nr,3000)
   select case(irecv(1))
   case (DO_PRECON)
      select case(loc_solver_cntl%preconditioner)
      case (MG_PRECONDITION)
         call multigrid(loc_grid,loc_procs,loc_linsys,io_cntl_noprint, &
                        loc_solver_cntl,loc_still_sequential)
      end select
      deallocate(irecv,stat=astat)
      if (associated(rrecv)) deallocate(rrecv,stat=astat)
   case (QUIT)
      deallocate(irecv,stat=astat)
      if (associated(rrecv)) deallocate(rrecv,stat=astat)
      exit
   end select
end do

end subroutine master_waits

!          ------
subroutine matvec(alpha,x,beta,y)
!          ------

!----------------------------------------------------
! This routine is called on processor 1 by GMRES to perform the matrix-vector
! multiply y := alpha*A*x + beta*y where A is the stiffness matrix
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(kind(0.0d0)) :: alpha,x(*),beta,y(*)
!----------------------------------------------------
! Local variables:

real(my_real), pointer :: rrecv(:)
integer, pointer :: irecv(:)
integer :: nproc, p, icount, i, ni, nr, astat
!----------------------------------------------------
! Begin executable code

nproc = num_proc(loc_procs)

! inform the other processors that a matvec is needed and send the vector

if (.not. loc_still_sequential) then
   call start_watch((/cpsolve,ctsolve/))
   icount = myneq
   do p=2,nproc
      call phaml_send(loc_procs,p,(/DO_MATVEC/),1, &
                      x(icount+1:icount+each_rhs(p)%n),each_rhs(p)%n,3000)
      icount = icount + each_rhs(p)%n
   end do
   call stop_watch((/cpsolve,ctsolve/))
endif

! copy the vector to a temporary

icount = 0
tempvec = 0.0_my_real
do i=1,loc_linsys%neq
   if (loc_linsys%iown(i) .and. loc_linsys%equation_type(i) /= DIRICHLET) then
      icount = icount + 1
      tempvec(i) = x(icount)
   endif
end do

! perform the matvec with the result in loc_linsys%solution

call matrix_times_vector(tempvec,loc_linsys%solution(1:loc_linsys%neq), &
                         loc_linsys,loc_procs,loc_still_sequential,3021,3022, &
                         3023,3024,3025,3026,nodirch=.true.,nocomm2=.true.)

! get this processor's part of the result

icount = 0
do i=1,loc_linsys%neq
   if (loc_linsys%iown(i) .and. loc_linsys%equation_type(i) /= DIRICHLET) then
      icount = icount + 1
      y(icount) = alpha*loc_linsys%solution(i) + beta*y(icount)
   endif
end do

! Gather the rest of the result from the other processors

if (.not. loc_still_sequential) then
   call start_watch((/cpsolve,ctsolve/))
   do i=2,nproc
      call phaml_recv(loc_procs,p,irecv,ni,rrecv,nr,3030)
      each_rhs(p)%val => rrecv
      if (associated(irecv)) deallocate(irecv,stat=astat)
   end do
   call stop_watch((/cpsolve,ctsolve/))

   do i=2,nproc
      if (associated(each_rhs(i)%val)) then
         y(icount+1:icount+each_rhs(i)%n) = alpha*each_rhs(i)%val + &
                                           beta*y(icount+1:icount+each_rhs(i)%n)
         icount = icount + each_rhs(i)%n
         deallocate(each_rhs(i)%val,stat=astat)
      endif
   end do
endif

end subroutine matvec

!          -------
subroutine precond(x,b)
!          -------

!----------------------------------------------------
! This routine performs preconditioning for GMRES
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(kind(0.0d0)) :: x(*), b(*)
!----------------------------------------------------
! Local variables:

integer :: i, p, ni, nr, astat, icount, nproc, lev
integer, pointer :: irecv(:)
real(my_real), pointer :: rrecv(:)
type(io_options) :: io_cntl_noprint
!----------------------------------------------------
! Begin executable code

select case(loc_solver_cntl%preconditioner)

! no preconditioning

case (NO_PRECONDITION)

   do i=1,neq
      x(i) = b(i)
   end do

! HBMG multigrid preconditioning

case (MG_PRECONDITION)

   nproc = num_proc(loc_procs)

   io_cntl_noprint = loc_io_cntl
   io_cntl_noprint%print_error_when = NEVER

! request preconditioning of the other processors and send unowned parts of b

   call start_watch((/cpsolve,ctsolve/))
   if (loc_still_sequential) then
      do p=2,nproc
         call phaml_send(loc_procs,p,(/DO_PRECON/),1,b(1:neq),neq,3000)
      end do
      call phaml_send(loc_procs,MASTER,(/DO_PRECON/),1,(/0.0_my_real/),0,3000)
   else
      icount = myneq
      do p=2,nproc
         call phaml_send(loc_procs,p,(/DO_PRECON/),1, &
                         b(icount+1:icount+each_rhs(p)%n),each_rhs(p)%n,3000)
         icount = icount + each_rhs(p)%n
      end do
      call phaml_send(loc_procs,MASTER,(/DO_PRECON/),1,(/0.0_my_real/),0,3000)
   endif
   call stop_watch((/cpsolve,ctsolve/))

! preserve right hand side of linear system

   tempvec = loc_linsys%rhs(1:loc_linsys%neq)

! copy my part of b to linsys%rhs

   icount = 0
   loc_linsys%rhs(1:loc_linsys%neq) = 0.0_my_real
   do i=1,loc_linsys%neq
      if (loc_linsys%iown(i) .and. loc_linsys%equation_type(i) /= DIRICHLET) then
         icount = icount + 1
         loc_linsys%rhs(i) = b(icount)
      endif
   end do

! get unowned right hand side entries

   if (.not. loc_still_sequential) then
      do lev=loc_linsys%nlev+1,2,-1
         call basis_change(lev,TO_HIER,loc_linsys)
      end do
      call sum_fudop_vect(loc_linsys%rhs(1:loc_linsys%neq),loc_procs, &
                          loc_linsys,3041,3042,3043)
      do lev=2,loc_linsys%nlev+1
         call basis_change(lev,TO_NODAL,loc_linsys)
      end do
   endif

! apply multigrid

   loc_linsys%solution(1:loc_linsys%neq) = 0.0_my_real
   call multigrid(loc_grid,loc_procs,loc_linsys,io_cntl_noprint, &
                  loc_solver_cntl,loc_still_sequential)

! get this processor's part of the result

   icount = 0
   do i=1,loc_linsys%neq
      if (loc_linsys%iown(i) .and. loc_linsys%equation_type(i) /= DIRICHLET) then
         icount = icount + 1
         x(icount) = loc_linsys%solution(i)
      endif
   end do

! Gather the rest of the result from the other processors

   if (.not. loc_still_sequential) then
      call start_watch((/cpsolve,ctsolve/))
      do i=2,nproc
         call phaml_recv(loc_procs,p,irecv,ni,rrecv,nr,3040)
         each_rhs(p)%val => rrecv
         if (associated(irecv)) deallocate(irecv,stat=astat)
      end do
      call stop_watch((/cpsolve,ctsolve/))

      do i=2,nproc
         if (associated(each_rhs(i)%val)) then
            x(icount+1:icount+each_rhs(i)%n) = each_rhs(i)%val
            icount = icount + each_rhs(i)%n
            deallocate(each_rhs(i)%val,stat=astat)
         endif
      end do
   endif

! restore original rhs

   loc_linsys%rhs(1:loc_linsys%neq) = tempvec

case default

   ierr = USER_INPUT_ERROR
   call fatal("Illegal choice of preconditioner for GMRES_SOLVER")
   stop

end select

end subroutine precond

!          --------
subroutine phaml_cg(linsys,procs,grid,io_cntl,solver_cntl,still_sequential)
!          --------

!----------------------------------------------------
! This routine routine solves the linear system using conjugate gradients.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(inout) :: linsys
type(proc_info), intent(in) :: procs
type(grid_type), intent(inout) :: grid
type(io_options), intent(in) :: io_cntl
type(solver_options), intent(in) :: solver_cntl
logical, intent(in) :: still_sequential
!----------------------------------------------------
! Local variables:

integer, pointer :: recv_int(:)
real(my_real), pointer :: recv_real(:)
real(my_real), allocatable :: x(:), r(:), p(:), q(:), &
                 myr(:), myp(:), myq(:), &
                 holdsoln(:), holdrhs(:)
real(my_real) :: rho, oldrho, alpha, beta, normr, pdotq, tol, normrhs, &
                 last_time
integer :: loop, i, j, neq, ni, nr, lev
logical :: dp
type(io_options) :: io_cntl_noprint
!----------------------------------------------------
! Begin executable code

! master just prints error, if requested, until told it's done

if (my_proc(procs) == MASTER) then
   if (io_cntl%print_error_when == FREQUENTLY .or. &
       io_cntl%print_error_when == TOO_MUCH) then
      loop=0
      do
         loop = loop+1
         call phaml_recv(procs,i,recv_int,ni,recv_real,nr,3040+loop)
         if (recv_int(1) == 0) exit
         if (loop == 1) then
            write(outunit,"(A)")
            write(outunit,"(A)") "L2 norm of linear system residual"
            write(outunit,"(A)")
            write(outunit,"(1x,A)") "     iteration    absolute residual    relative residual      reduction"
            write(outunit,"(A)")
            write(outunit,"(SS,1P,1X,I11,2E18.10E2)") 0,recv_real(1),recv_real(1)/recv_real(2)
         else
            write(outunit,"(SS,1P,1X,I11,3E18.10E2)") loop-1,recv_real(1),recv_real(1)/recv_real(2),recv_real(1)/last_time
         endif
         last_time = recv_real(1)
         if (associated(recv_int)) deallocate(recv_int)
         if (associated(recv_real)) deallocate(recv_real)
      end do
      if (associated(recv_int)) deallocate(recv_int)
      if (associated(recv_real)) deallocate(recv_real)
   endif
   return
endif

allocate(x(linsys%neq), r(linsys%neq), p(linsys%neq), q(linsys%neq), &
         myr(linsys%neq), myp(linsys%neq), myq(linsys%neq), &
         holdsoln(linsys%neq), holdrhs(linsys%neq))

! handy constants

neq = linsys%neq
dp = my_real == kind(1.0d0)
if (solver_cntl%krylov_tol == KRYLOV_ERREST_TOL) then
   call error_estimate(grid,procs,errest_L2=tol)
   if (.not. still_sequential) tol = sqrt(phaml_global_sum(procs,tol**2,3001))
   tol = tol/100
else
   tol = solver_cntl%krylov_tol
endif
io_cntl_noprint = io_cntl
io_cntl_noprint%print_error_when = NEVER

! hold solution to restore Dirichlet boundary conditions later, and rhs

holdsoln = linsys%solution(1:neq)
holdrhs  = linsys%rhs(1:neq)

! compute the norm of the rhs of the equations CG is solving

where (linsys%equation_type(1:neq) /= DIRICHLET) linsys%solution(1:neq) = 0
call matrix_times_vector(linsys%solution(1:neq),r,linsys,procs, &
                         still_sequential,3080,3081,3082,3083,3084,3085, &
                         nocomm2=.true.,nodirch=.true.)
r = linsys%rhs(1:neq) - r
do i=1,neq
   if (linsys%equation_type(i) /= DIRICHLET) then
      do j=linsys%begin_row(i),linsys%end_row(i)
         if (linsys%column_index(j) == NO_ENTRY) cycle
         if (linsys%equation_type(linsys%column_index(j)) == DIRICHLET) then
            r(i) = r(i) - linsys%matrix_val(j)*holdsoln(linsys%column_index(j))
         endif
      end do
   endif
end do
where (linsys%equation_type(1:neq) == DIRICHLET) r = 0.0_my_real
where (.not. linsys%iown(1:neq)) r = 0.0_my_real
if (dp) then
   normrhs = ddot(neq,r,1,r,1)
else
   normrhs = sdot(neq,r,1,r,1)
endif
if (.not. still_sequential) then
   call start_watch((/cpsolve,ctsolve/))
   normrhs = phaml_global_sum(procs,normrhs,3006)
   call stop_watch((/cpsolve,ctsolve/))
endif
normrhs = sqrt(normrhs)

! compute initial residual r <-- b - Ax

linsys%solution(1:neq) = holdsoln
linsys%rhs(1:neq) = holdrhs

x = linsys%solution(1:neq)
call matrix_times_vector(x,r,linsys,procs,still_sequential,3010,3011,3012, &
                         3013,3014,3015,nocomm2=.true.,nodirch=.true.)
r = linsys%rhs(1:neq) - r

! subtract out Dirichlet boundary conditions, and set them to 0 so they
! aren't used in the preconditioner

do i=1,neq
   if (linsys%equation_type(i) /= DIRICHLET) then
      do j=linsys%begin_row(i),linsys%end_row(i)
         if (linsys%column_index(j) == NO_ENTRY) cycle
         if (linsys%equation_type(linsys%column_index(j)) == DIRICHLET) then
            r(i) = r(i) - linsys%matrix_val(j)*holdsoln(linsys%column_index(j))
         endif
      end do
   else
      linsys%solution(i) = 0.0_my_real
   endif
end do
where (linsys%equation_type(1:neq) == DIRICHLET) r = 0.0_my_real
myr = r
where (.not. linsys%iown(1:neq)) myr = 0.0_my_real

! repeat until converged

rho = 0.0_my_real
loop = 0

do
   loop = loop+1

! test for convergence

   if (dp) then
      normr = ddot(neq,myr,1,myr,1)
   else
      normr = sdot(neq,myr,1,myr,1)
   endif
   if (.not. still_sequential) then
      call start_watch((/cpsolve,ctsolve/))
      normr = phaml_global_sum(procs,normr,3002)
      call stop_watch((/cpsolve,ctsolve/))
   endif
   normr = sqrt(normr)
   if (io_cntl%print_error_when == FREQUENTLY .or. &
       io_cntl%print_error_when == TOO_MUCH) then
      if (my_proc(procs) == 1) then
         call phaml_send(procs,MASTER,(/1/),1,(/normr,normrhs/),2,3040+loop)
      endif
   endif
   if (normr/normrhs < tol) exit
   if (loop > solver_cntl%krylov_iter) then
      call warning("CG did not converge to tolerance",reallist=(/tol,normr/))
      exit
   endif

! apply preconditioner; solve Mz=r  where z stays in linsys%solution

   select case(solver_cntl%preconditioner)

   case (MG_PRECONDITION)
      linsys%rhs(1:neq) = myr
      if (.not. still_sequential) then
         do lev=linsys%nlev+1,2,-1
            call basis_change(lev,TO_HIER,linsys)
         end do
         call sum_fudop_vect(linsys%rhs(1:neq),procs,linsys,3071,3072,3073)
         do lev=2,linsys%nlev+1
            call basis_change(lev,TO_NODAL,linsys)
         end do
      endif
      linsys%solution(1:neq) = 0.0_my_real
      call multigrid(grid,procs,linsys,io_cntl_noprint,solver_cntl, &
                     still_sequential,no_master=.true.)
      if (.not. still_sequential) then
        if (.not. new_comm) then ! TEMP090204
         call exchange_fudop_vect(linsys%solution(1:neq),procs,linsys,3091,3092,3093) ! TEMP090204
        endif ! TEMP090204
      endif

   case (NO_PRECONDITION)
      linsys%solution(1:neq) = r

   end select

! rho <-- (r,z)

   oldrho = rho
   if (dp) then
      rho = ddot(neq,myr,1,linsys%solution(1:neq),1)
   else
      rho = sdot(neq,myr,1,linsys%solution(1:neq),1)
   endif
   if (.not. still_sequential) then
      call start_watch((/cpsolve,ctsolve/))
      rho = phaml_global_sum(procs,rho,3003)
      call stop_watch((/cpsolve,ctsolve/))
   endif

! p <-- z + beta p where beta=0 first time and rho/oldrho after that

   if (loop == 1) then
      p = linsys%solution(1:neq)
   else
      beta = rho/oldrho
      if (dp) then
         call daxpy(neq,beta,p,1,linsys%solution(1:neq),1)
      else
         call saxpy(neq,beta,p,1,linsys%solution(1:neq),1)
      endif
      p = linsys%solution(1:neq)
   endif
   where (linsys%equation_type(1:neq) == DIRICHLET) p = 0.0_my_real
   myp = p
   where (.not. linsys%iown(1:neq)) myp = 0.0_my_real

! q <-- Ap

   call matrix_times_vector(p,q,linsys,procs,still_sequential,3020,3021,3022, &
                            0,0,0,nodirch=.true.,nocomm2=.true.)
   where (linsys%equation_type(1:neq) == DIRICHLET) q = 0.0_my_real
   myq = q
   where (.not. linsys%iown(1:neq)) myq = 0.0_my_real

! alpha <-- rho / (p,q)

   if (dp) then
      pdotq = ddot(neq,myp,1,myq,1)
   else
      pdotq = sdot(neq,myp,1,myq,1)
   endif
   if (.not. still_sequential) then
      call start_watch((/cpsolve,ctsolve/))
      pdotq = phaml_global_sum(procs,pdotq,3004)
      call stop_watch((/cpsolve,ctsolve/))
   endif
   alpha = rho/pdotq

! x <-- x + alpha p

   if (dp) then
      call daxpy(neq,alpha,p,1,x,1)
   else
      call saxpy(neq,alpha,p,1,x,1)
   endif

! r <-- r - alpha q

   if (dp) then
      call daxpy(neq,-alpha,q,1,r,1)
   else
      call saxpy(neq,-alpha,q,1,r,1)
   endif
   where (linsys%equation_type(1:neq) == DIRICHLET) r = 0.0_my_real
   myr = r
   where (.not. linsys%iown(1:neq)) myr = 0.0_my_real

end do

! tell the master we're done

if (io_cntl%print_error_when == FREQUENTLY .or. &
    io_cntl%print_error_when == TOO_MUCH) then
   if (my_proc(procs) == 1) then
      call phaml_send(procs,MASTER,(/0/),1,(/0.0_my_real/),0,3040+loop+1)
   endif
endif

! copy solution

linsys%solution(1:neq) = x

! restore Dirichlet entries of solution and previous solution at unowned
! equations to get the edge solutions that the owner doesn't have because
! it's a different level

where (linsys%equation_type(1:neq) == DIRICHLET .or. &
       (.not. linsys%iown(1:neq))) linsys%solution(1:neq) = holdsoln

! exchange solution at unowned equations

if (.not. still_sequential) then
   call exchange_fudop_vect(linsys%solution(1:),procs,linsys,3031,3032,3033)
endif

! restore right hand side

linsys%rhs(1:neq) = holdrhs

deallocate(x, r, p, q, myr, myp, myq, holdsoln, holdrhs)
end subroutine phaml_cg

end module krylov
