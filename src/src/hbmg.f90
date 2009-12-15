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

module hbmg

!----------------------------------------------------
! This module contains routines for the hierarchical basis multigrid method
!
! communication tags in this module are of the form 14xx and 16xx
! TEMP090127 new communication approach
!----------------------------------------------------

!----------------------------------------------------
! Other modules used are:

use global
use message_passing
use gridtype_mod
use linsystype_mod
use linsys_util
use lapack_solve
use linsys_io
use error_estimators

!----------------------------------------------------

implicit none
private
public multigrid, mg_precon

!----------------------------------------------------
! The following variables are defined:

!----------------------------------------------------

contains

!          ---------
subroutine multigrid(grid,procs,linear_system,io_cntl,solver_cntl, &
                     still_sequential,maxlev,no_master)
!          ---------

!----------------------------------------------------
! This routine solves the linear system via multigrid.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(proc_info), intent(in) :: procs
type(linsys_type), intent(inout) :: linear_system
type(io_options), intent(in) :: io_cntl
type(solver_options), intent(in) :: solver_cntl
logical, intent(in) :: still_sequential
integer, intent(in), optional :: maxlev
logical, intent(in), optional :: no_master

!----------------------------------------------------
! Local variables:

integer :: ncyc, cycle, lev, nu1, nu2, nu1ho, nu2ho, nlev, proc, ni, nr
real(my_real) :: mg_tol, resid
integer, pointer :: irecv(:)
real(my_real), pointer :: rrecv(:)
logical :: fudop_comm, conventional_comm, yes_master

!----------------------------------------------------
! Begin executable code

if (present(no_master)) then
   yes_master = .not. no_master
else
   yes_master = .true.
endif

fudop_comm = solver_cntl%mg_comm == MGCOMM_FUDOP .and. &
             num_proc(procs) > 1 .and. &
             .not. still_sequential
conventional_comm = solver_cntl%mg_comm == MGCOMM_CONVENTIONAL .and. &
             num_proc(procs) > 1 .and. &
             .not. still_sequential
ncyc = solver_cntl%ncycle
nu1 = solver_cntl%prerelax
nu2 = solver_cntl%postrelax
nu1ho = solver_cntl%prerelax_ho
nu2ho = solver_cntl%postrelax_ho
if (present(maxlev)) then
   nlev = maxlev
else
   if (linear_system%maxdeg == 1) then
      nlev = linear_system%nlev
   else
      nlev = linear_system%nlev+1
   endif
endif

! determine tolerance of residual for termination

mg_tol = solver_cntl%mg_tol
if (mg_tol == MG_ERREST_TOL) then
   if (my_proc(procs) == MASTER) then
      call phaml_recv(procs,proc,irecv,ni,rrecv,nr,1400)
      mg_tol = rrecv(1)
      deallocate(rrecv)
   else
! TEMP for systems of equations, should be a combination not the max
      call error_estimate(grid,procs,errest_L2=mg_tol)
      mg_tol = sqrt(phaml_global_sum(procs,mg_tol**2,1401))
      mg_tol = max(100*epsilon(1.0_my_real),mg_tol/100)
      if (my_proc(procs) == 1 .and. yes_master) then
         call phaml_send(procs,MASTER,(/0/),0,(/mg_tol/),1,1400)
      endif
   endif
endif

resid = 0.0_my_real

! master process just prints the error (if requested) and returns

if (my_proc(procs) == MASTER) then
   if (io_cntl%print_error_when == FREQUENTLY .or. &
       io_cntl%print_error_when == TOO_MUCH) then
      call linsys_residual(linear_system,procs,still_sequential,0,1410,.true., &
                           .true.)
   endif
   do cycle=1,ncyc
      if (io_cntl%print_error_when == FREQUENTLY .or. &
          io_cntl%print_error_when == TOO_MUCH) then
         call linsys_residual(linear_system,procs,still_sequential,cycle, &
                              1410+cycle,.true.,.false.,relresid=resid)
      else
         if (mg_tol /= MG_NO_TOL) then
            call linsys_residual(linear_system,procs,still_sequential,cycle, &
                                 1410+cycle,.false.,.false.,relresid=resid)
         endif
      endif
      if (resid < mg_tol) exit
   end do
   return
endif

! make sure the solution at unowned points is current

if (fudop_comm .or. conventional_comm) then
   call exchange_fudop_vect(linear_system%solution(1:),procs, &
                            linear_system,1402,1403,1404)
endif

! print the error before solution (if requested)

if (io_cntl%print_error_when == FREQUENTLY .or. &
    io_cntl%print_error_when == TOO_MUCH) then
   call linsys_residual(linear_system,procs,still_sequential,0,1410,.true., &
                        .true.,no_master=no_master)
endif

! initialize r_other

if (fudop_comm .or. conventional_comm) then
   call init_r_other(linear_system%solution(1:))
else
   linear_system%r_mine = 0.0_my_real
   linear_system%r_others = 0.0_my_real
endif

! repeat ncycle times or until L2 norm of residual < mg_tol

do cycle=1,ncyc

! for each level from finest to coarsest, relaxation and restriction

   do lev=nlev,2,-1
      if (lev == linear_system%nlev+1) then
         call relax_ho(lev,nu1ho,linear_system,-1,procs, &
                       fudop_comm.or.conventional_comm)
      else
         call relax(lev,nu1,linear_system,conventional_comm)
      endif
      call basis_change(lev,TO_HIER,linear_system)
      call rhs_minus_Axh(linear_system,lev)
      if (conventional_comm) then
         call exchange_fudop_soln_residual(procs,linear_system, &
                                           1420+cycle+lev,1430+cycle+lev, &
                                           1440+cycle+lev,max_lev=lev)
      endif
   end do

! fix solution values and residuals via information from other processors

   if (fudop_comm) then
      call exchange_fudop_soln_residual(procs,linear_system, &
                                        1420+cycle,1430+cycle,1440+cycle)
   endif

! solve on coarsest grid

   call coarse_solve(linear_system)

! for each level from coarsest to finest, prolongation and relaxation

   do lev=2,nlev
      call rhs_plus_Axh(linear_system,lev)
      if (conventional_comm) then
         call exchange_fudop_soln_residual(procs,linear_system, &
                                           1450+cycle+lev,1460+cycle+lev, &
                                           1470+cycle+lev,max_lev=lev)
      endif
      call basis_change(lev,TO_NODAL,linear_system)
      if (lev == linear_system%nlev+1) then
         call relax_ho(lev,nu2ho,linear_system,1,procs, &
                       fudop_comm.or.conventional_comm)
      else
         call relax(lev,nu2,linear_system,conventional_comm)
      endif
   end do

! fix solution values via information from other processors

   if (fudop_comm .or. conventional_comm) then
      call exchange_fudop_vect(linear_system%solution(1:),procs, &
                               linear_system,1450+cycle,1460+cycle, &
                               1470+cycle)
   endif

! print error after this cycle (if requested) and test for small residual

   if (io_cntl%print_error_when == FREQUENTLY .or. &
       io_cntl%print_error_when == TOO_MUCH) then
      call linsys_residual(linear_system,procs,still_sequential,cycle, &
                           1410+cycle,.true.,.false.,relresid=resid,no_master=no_master)
   else
      if (mg_tol /= MG_NO_TOL) then
         call linsys_residual(linear_system,procs,still_sequential,cycle, &
                              1410+cycle,.false.,.false.,relresid=resid,no_master=no_master)
      endif
   endif
   if (resid < mg_tol) exit

end do ! next cycle

return

contains

!          ------------
subroutine init_r_other(solution)
!          ------------

!----------------------------------------------------
! This routine initializes r_other before the first V-cycle
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(inout) :: solution(:)
!----------------------------------------------------
! Local variables:

real(my_real) :: hold_soln(size(solution))
!----------------------------------------------------
! Begin executable code

! clear out residuals

linear_system%r_mine = 0.0_my_real
linear_system%r_others = 0.0_my_real

! hold the solution so that this doesn't change it

hold_soln = solution

! do the first half of a V cycle to generate the residual

do lev=nlev,2,-1
   if (lev == linear_system%nlev+1) then
      call relax_ho(lev,nu1ho,linear_system,-1,procs, &
                    fudop_comm.or.conventional_comm)
   else
      call relax(lev,nu1,linear_system,conventional_comm)
   endif
   call basis_change(lev,TO_HIER,linear_system)
   call rhs_minus_Axh(linear_system,lev)
   if (conventional_comm) then
      call exchange_fudop_soln_residual(procs,linear_system,1600+lev,1610+lev, &
                                        1620+lev,max_lev=lev)
   endif
end do

! exchange the residual

if (fudop_comm) then
   call exchange_fudop_soln_residual(procs,linear_system,1600,1610,1620)
endif

! second half of V cycle to restore to nodal form, but don't need relaxation

do lev=2,nlev
   call rhs_plus_Axh(linear_system,lev)
   if (conventional_comm) then
      call exchange_fudop_soln_residual(procs,linear_system, &
                                        1630+lev,1640+lev,1650+lev, &
                                        no_soln=.true.,max_lev=lev)
   endif
   call basis_change(lev,TO_NODAL,linear_system)
end do

! restore the solution

solution = hold_soln

return
end subroutine init_r_other

end subroutine multigrid

!          -----
subroutine relax(lev,nu,matrix,conventional_comm)
!          -----

!----------------------------------------------------
! This routine perform nu/2 iterations of red-(local)black relaxation, where
! the red points are the level lev points.
! Half an iteration means red relaxation.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: lev, nu
type(linsys_type), intent(inout) :: matrix
logical, intent(in) :: conventional_comm

!----------------------------------------------------
! Local variables:

integer :: iter, eq, nblack, i, neigh, allocstat
integer :: j
logical :: red_only
integer, allocatable :: black_list(:)
logical(small_logical), allocatable :: on_black_list(:)
!----------------------------------------------------
! Begin executable code
 
if (nu > 1) then
   allocate(on_black_list(matrix%neq),black_list(matrix%neq),stat=allocstat)
   if (allocstat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed for black_list in relax")
      return
   endif
   nblack = 0
   on_black_list = .false.
endif

! repeat nu/2 times, but make sure the red get relaxed when nu is odd

do iter = 1, (nu+1)/2

! will this be red relaxation only?

!                nu is odd           last iteration
   red_only = (2*(nu/2) /= nu .and. iter == (nu+1)/2)

! for each non-Dirichlet equation of level lev

   do eq = matrix%begin_level(lev), matrix%begin_level(lev+1)-1

      if (matrix%equation_type(eq) == DIRICHLET) cycle

! red relaxation

      matrix%solution(eq) = matrix%rhs(eq) + matrix%r_mine(eq) + &
                            matrix%r_others(eq)
      do i=matrix%begin_row(eq)+1,matrix%end_row(eq)
         matrix%solution(eq) = matrix%solution(eq) - &
                   matrix%matrix_val(i)*matrix%solution(matrix%column_index(i))
      end do
      matrix%solution(eq) = matrix%solution(eq) / &
                            matrix%matrix_val(matrix%begin_row(eq))

! first time, add neighbors to black lists, if there will be black relaxation

      if ((.not.red_only) .and. iter==1) then
         do i=matrix%begin_row(eq)+1,matrix%end_row(eq)
            neigh = matrix%column_index(i)
            if (neigh == NO_ENTRY) cycle
            if (on_black_list(neigh)) cycle
            if (neigh >= matrix%begin_level(lev)) cycle ! must be a lower level (black)
            on_black_list(neigh) = .true.
            nblack = nblack + 1
            black_list(nblack) = neigh
         end do
      endif

   end do ! next equation

! black relaxation, if desired

   if (red_only) exit

   do i = 1,nblack
      eq = black_list(i)
      if (matrix%equation_type(eq) == DIRICHLET) cycle
      if (.not. red_only) then
         matrix%solution(eq) = matrix%rhs(eq) + matrix%r_mine(eq) + &
                               matrix%r_others(eq)
      endif
      do j=matrix%begin_row(eq)+1,matrix%end_row(eq)
         if (.not. red_only) then
            matrix%solution(eq) = matrix%solution(eq) - &
                matrix%matrix_val(j)*matrix%solution(matrix%column_index(j))
         endif
      end do
      if (.not. red_only) then
         matrix%solution(eq) = matrix%solution(eq) / &
                               matrix%matrix_val(matrix%begin_row(eq))
      endif
   end do

end do ! next iteration

! free memory

if (nu > 1) then
   if (allocated(black_list)) then
      deallocate(on_black_list,black_list,stat=allocstat)
      if (allocstat /= 0) then
         call warning("deallocation failed in relax")
      endif
   endif
endif

return
end subroutine relax

!          --------
subroutine relax_ho(lev,nu,matrix,dir,procs,communicate)
!          --------

!----------------------------------------------------
! This routine performs one direction of the p-multigrid cycle.  If dir==-1
! it performs nu Gauss-Seidel iterations with the entire (condensed) matrix,
! then with the matrix up to maxdeg-1, then maxdeg-2, ... 2.  If dir==1 it
! starts at 2 and works up to maxdeg.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: lev, nu
type(linsys_type), intent(inout) :: matrix
integer, intent(in) :: dir
type(proc_info), intent(in) :: procs
logical, intent(in) :: communicate

!----------------------------------------------------
! Local variables:

integer :: iter, eq, i, j, gsiter, loop
integer :: object_type,basis_rank,system_rank,eqdeg
integer :: countdeg(matrix%maxdeg), startdeg(matrix%maxdeg+1), &
           order(matrix%begin_level(lev+1)-1)
logical :: isneigh(matrix%begin_level(lev+1)-1)
!----------------------------------------------------
! Begin executable code
 
! TEMP Determine an order for the equations such that all equations of
!      the same degree are together. This should be done once and for all when
!      the matrix is created, maybe even actually placing them in that order.

countdeg = 0
do eq=1,matrix%begin_level(lev+1)-1
   call eq_to_grid(matrix,matrix%gid(eq),object_type,basis_rank, &
                   system_rank)
   if (object_type == VERTEX_ID) then
      eqdeg = 1
   else
      eqdeg = basis_rank+1
   endif
   countdeg(eqdeg) = countdeg(eqdeg)+1
end do
startdeg(1) = 1
do i=2,matrix%maxdeg+1
   startdeg(i) = startdeg(i-1)+countdeg(i-1)
end do
countdeg = startdeg(1:matrix%maxdeg)
do eq=1,matrix%begin_level(lev+1)-1
   call eq_to_grid(matrix,matrix%gid(eq),object_type,basis_rank, &
                   system_rank)
   if (object_type == VERTEX_ID) then
      eqdeg = 1
   else
      eqdeg = basis_rank+1
   endif
   order(countdeg(eqdeg)) = eq
   countdeg(eqdeg) = countdeg(eqdeg)+1
end do

! end TEMP

if (dir == -1) then

   do iter = matrix%maxdeg,2,-1
    isneigh = .false. ! TEMP090401
    do gsiter=1,nu
     do loop=1,2
      if (.not. new_comm .and. loop==2) cycle ! TEMP090127

      if (communicate) then
        if (new_comm) then ! TEMP090127
         if (loop==1) then
            call send_neigh_vect(matrix%solution(1:),procs,matrix,1660+iter, &
                                 min(iter+1,matrix%maxdeg))
         else
            call recv_neigh_vect(matrix%solution(1:),procs,matrix,1660+iter)
         endif
        else ! TEMP090127
         call exchange_neigh_vect(matrix%solution(1:),procs, & ! TEMP090127
                                  matrix,1660+iter,1670+iter,1680+iter) ! TEMP090127
        endif ! TEMP090127
      endif
      do i = matrix%begin_level(lev+1)-1,1,-1
         eq = order(i)
         if (new_comm) then ! TEMP090127
         if ((loop==1 .and. matrix%nn_comm_remote_neigh(eq)) .or. &
             (loop==2 .and. .not. matrix%nn_comm_remote_neigh(eq))) then
            cycle
         endif
         endif ! TEMP090127
         if (matrix%equation_type(eq) == DIRICHLET) cycle
         call eq_to_grid(matrix,matrix%gid(eq),object_type,basis_rank, &
                         system_rank)
         if (object_type == VERTEX_ID) then
            eqdeg = 1
         else
            eqdeg = basis_rank+1
         endif
! TEMP not a very efficient way to do this
         if (eqdeg > iter) cycle
         if (eqdeg < iter .and. .not. isneigh(eq)) cycle
         if (.not. communicate .or. matrix%iown(eq)) then
            matrix%solution(eq) = matrix%rhs(eq) + matrix%r_mine(eq) + &
                                  matrix%r_others(eq)
         endif
         do j=matrix%begin_row(eq)+1,matrix%end_row(eq)
            if (matrix%column_index(j) == NO_ENTRY) cycle
            if (eqdeg==iter) isneigh(matrix%column_index(j)) = .true.
            if (.not. communicate .or. matrix%iown(eq)) then
               matrix%solution(eq) = matrix%solution(eq) - &
                   matrix%matrix_val(j)*matrix%solution(matrix%column_index(j))
            endif
         end do
         if (.not. communicate .or. matrix%iown(eq)) then
            matrix%solution(eq) = matrix%solution(eq) / &
                                  matrix%matrix_val(matrix%begin_row(eq))
         endif
      end do

     end do ! all neighbors local or not
    end do ! next Gauss-Seidel iteration
   end do ! next iteration
   if (communicate) then
      call exchange_fudop_vect(matrix%solution(1:),procs, &
                               matrix,1661,1671,1681)
   endif

elseif (dir == 1) then

   do iter = 2,matrix%maxdeg
    isneigh = .false.
    do gsiter=1,nu
     do loop=1,2
      if (.not. new_comm .and. loop==2) cycle ! TEMP090127

      if (communicate) then
       if (new_comm) then ! TEMP090127
         if (loop==1) then
            call send_neigh_vect(matrix%solution(1:),procs,matrix,1660+iter, &
                                 iter-1)
         else
            call recv_neigh_vect(matrix%solution(1:),procs,matrix,1660+iter)
          endif
       else ! TEMP090127
         call exchange_neigh_vect(matrix%solution(1:),procs, & ! TEMP090127
                                  matrix,1480+iter,1490+iter,1690+iter) ! TEMP090127
       endif ! TEMP090127
      endif
      do i = matrix%begin_level(lev+1)-1,1,-1
         eq = order(i)
         if (new_comm) then ! TEMP090127
         if ((loop==1 .and. matrix%nn_comm_remote_neigh(eq)) .or. &
             (loop==2 .and. .not. matrix%nn_comm_remote_neigh(eq))) then
            cycle
         endif
         endif ! TEMP090127
         if (matrix%equation_type(eq) == DIRICHLET) cycle
         if (communicate .and. .not. matrix%iown(eq)) cycle
         call eq_to_grid(matrix,matrix%gid(eq),object_type,basis_rank, &
                         system_rank)
         if (object_type == VERTEX_ID) then
            eqdeg = 1
         else
            eqdeg = basis_rank+1
         endif
         if (eqdeg > iter) cycle
         if (eqdeg < iter .and. .not. isneigh(eq)) cycle
         matrix%solution(eq) = matrix%rhs(eq) + matrix%r_mine(eq) + &
                               matrix%r_others(eq)
         do j=matrix%begin_row(eq)+1,matrix%end_row(eq)
            if (matrix%column_index(j) == NO_ENTRY) cycle
            if (eqdeg==iter) isneigh(matrix%column_index(j)) = .true.
            matrix%solution(eq) = matrix%solution(eq) - &
                    matrix%matrix_val(j)*matrix%solution(matrix%column_index(j))
         end do
         matrix%solution(eq) = matrix%solution(eq) / &
                               matrix%matrix_val(matrix%begin_row(eq))
      end do

     end do ! all neighbors local or not
    end do ! next Gauss-Seidel iteration
   end do ! next iteration

else

   call fatal("bad value for dir in relax_ho",intlist=(/dir/))
   stop

endif

end subroutine relax_ho

!          ------------
subroutine rhs_plus_Axh(matrix,lev)
!          ------------

!----------------------------------------------------
! This routine computes the addition of the coarse-fine block of the
! hierarchical matrix times the fine part of the hierarchical solution
! vector and removes it from r_mine and r_others
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(inout) :: matrix
integer, intent(in) :: lev

!----------------------------------------------------
! Local variables:

integer :: eq, col

!----------------------------------------------------
! Begin executable code

do eq = matrix%begin_level(lev),matrix%begin_level(lev+1)-1
   do col = matrix%begin_row(eq)+1,matrix%end_row(eq)
      if (matrix%column_index(col) == NO_ENTRY) cycle
      if (matrix%column_index(col) >= matrix%begin_level(lev)) cycle
      if (matrix%iown(eq)) then
         matrix%r_mine(matrix%column_index(col)) = &
                                     matrix%r_mine(matrix%column_index(col)) + &
                                     matrix%matrix_val(col)*matrix%solution(eq)
      else
         matrix%r_others(matrix%column_index(col)) = &
                                   matrix%r_others(matrix%column_index(col)) + &
                                   matrix%matrix_val(col)*matrix%solution(eq)
      endif
   end do
end do

return
end subroutine rhs_plus_Axh

!          -------------
subroutine rhs_minus_Axh(matrix,lev)
!          -------------

!----------------------------------------------------
! This routine computes the subtraction of the coarse-fine block of the
! hierarchical matrix times the fine part of the hierarchical solution
! vector and stores it in r_mine and r_others
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(inout) :: matrix
integer, intent(in) :: lev

!----------------------------------------------------
! Local variables:

integer :: eq, col

!----------------------------------------------------
! Begin executable code

do eq = matrix%begin_level(lev),matrix%begin_level(lev+1)-1
   do col = matrix%begin_row(eq)+1,matrix%end_row(eq)
      if (matrix%column_index(col) == NO_ENTRY) cycle
      if (matrix%column_index(col) >= matrix%begin_level(lev)) cycle
      if (matrix%iown(eq)) then
         matrix%r_mine(matrix%column_index(col)) = &
                                     matrix%r_mine(matrix%column_index(col)) - &
                                     matrix%matrix_val(col)*matrix%solution(eq)
      else
         matrix%r_others(matrix%column_index(col)) = &
                                   matrix%r_others(matrix%column_index(col)) - &
                                   matrix%matrix_val(col)*matrix%solution(eq)
      endif
   end do
end do

return
end subroutine rhs_minus_Axh

!          ------------
subroutine coarse_solve(linear_system)
!          ------------

!----------------------------------------------------
! This routine solves the coarse grid problem using the LAPACK routines for
! symmetric positive definite band matricies.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(inout) :: linear_system
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

! All partitions have the same coarse grid.  Let each processor do the
! coarse grid for itself.

! If the coarse grid matrix has not been put into band form and
! factorized, do so.

if (.not. linear_system%coarse_band_exists) then
   call make_lapack_symm_band(1,linear_system,linear_system%coarse_matrix)
   if (ierr /= NO_ERROR) return
   linear_system%coarse_band_exists = .true.
endif

! Solve the coarse grid system

call lapack_spd(1,linear_system,linear_system%coarse_matrix)

end subroutine coarse_solve

!          ---------
subroutine mg_precon(invec,outvec,matrix,grid,procs,io_cntl, &
                     solver_cntl,still_sequential)
!          ---------

!----------------------------------------------------
! This routine applies a multigrid V cycle as a preconditioner.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: invec(:)
real(my_real), intent(out) :: outvec(:)
type(linsys_type), intent(inout) :: matrix
type(grid_type), intent(inout) :: grid
type(proc_info), intent(in) :: procs
type(io_options), intent(in) :: io_cntl
type(solver_options), intent(in) :: solver_cntl
logical, intent(in) :: still_sequential
!----------------------------------------------------
! Local variables:

real(my_real) :: holdrhs(size(matrix%rhs)), holdsoln(0:size(matrix%rhs))
!----------------------------------------------------
! Begin executable code

! Keep rhs and solution

holdrhs = matrix%rhs
holdsoln = matrix%solution

! Copy the invec to rhs; the size of invec should be the same as rhs

matrix%rhs = invec

! Set the initial guess to 0.0

matrix%solution(1:) = 0.0_my_real

! Set Dirichlet points

where (matrix%equation_type == DIRICHLET) matrix%solution(1:) = matrix%rhs

! Apply a multigrid cycle

call multigrid(grid,procs,matrix,io_cntl,solver_cntl,still_sequential)

! Copy solution (which now contains the preconditioner times invec) to outvec

outvec = matrix%solution(1:)

! Restore rhs and solution

matrix%rhs = holdrhs
matrix%solution = holdsoln

end subroutine mg_precon

end module hbmg
