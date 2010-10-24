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

module linsys_io

!----------------------------------------------------
! This module contains routines for I/O related to the linear system,
! including computing and printing the error.
!
! communication tags in this module are of the form 13xx
!----------------------------------------------------

!----------------------------------------------------
! Other modules used are:

use global
use message_passing
use stopwatch
use hash_eq_mod
use gridtype_mod
use linsystype_mod
use linsys_util
use make_linsys
use error_estimators
use quadrature_rules
use evaluate
use sysdep

!----------------------------------------------------

implicit none
private
public print_linear_system, print_linsys_info, print_error_info, norm_error, &
       norm_solution, norm_true, store_matrix

!----------------------------------------------------
! Non-module procedures used are:

interface

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

   subroutine pdecoefs(x,y,cxx,cxy,cyy,cx,cy,c,rs)
   use global
   real(my_real), intent(in) :: x,y
   real(my_real), intent(out) :: cxx(:,:),cxy(:,:),cyy(:,:),cx(:,:),cy(:,:), &
                                 c(:,:),rs(:)
   end subroutine pdecoefs

end interface

!----------------------------------------------------
! Undocumented feature.  Write the number of degrees of freedom and
! various norms of the error and error estimates and computation time
! to a file whenever the error is printed and the number of dof has changed.
! You should open and close unit convfileunit in the master program to use this.

! The variable save_convergence is defined in global so it can be changed
! in the master main program.

integer, save :: last_dof=0
logical, save :: first_conv=.true.
!----------------------------------------------------

contains

!----------------------------------------------------------------
! CODE FOR I/O
!----------------------------------------------------------------

!          -------------------
subroutine print_linear_system(ls,procs,grid)
!          -------------------

!----------------------------------------------------
! This routine prints the linear system
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(in) :: ls
type (proc_info), intent(in) :: procs
type(grid_type), intent(in), optional :: grid
!----------------------------------------------------
! Local variables:

integer :: i,objtype,brank,srank,glid

!----------------------------------------------------
! Begin executable code

if (my_proc(procs) == MASTER) return

write(outunit,"(A)")  "--------------"
write(outunit,"(A)")  "Linear System:"
write(outunit,"(A)")  "--------------"
write(outunit,"(A)")

write(outunit,"(A,I11)")  "  number of levels ",ls%nlev
write(outunit,"(A)")  "  beginning of each level"
write(outunit,"(7I11)")  ls%begin_level
write(outunit,"(A)")  "  beginning of each row"
write(outunit,"(7I11)")  ls%begin_row
write(outunit,"(A)")  "  end of each row"
write(outunit,"(7I11)")  ls%end_row

if (present(grid)) then
   write(outunit,"(A)")  "  equation correspondence: eqn, grid point, type, basis rank, system rank"
   write(outunit,"(A,I2,A,I2,A,I2)")  "    types are ",VERTEX_ID," vertex; ",EDGE_ID," edge; ", &
                     ELEMENT_ID," element;"
   do i=1,ls%neq_vert + ls%neq_edge + ls%neq_face
      call eq_to_grid(ls,ls%gid(i),objtype,brank,srank,glid,grid)
      write(outunit,"(A,5I11)") "    ",i,glid,objtype,brank,srank
   end do
endif

write(outunit,"(A)")  "  matrix values and column indices"
do i=1,ls%neq_vert + ls%neq_edge + ls%neq_face
   write(outunit,"(SS,1P,4E18.10E2)")  ls%matrix_val(ls%begin_row(i):ls%end_row(i))
   write(outunit,"(4I21)")  ls%column_index(ls%begin_row(i):ls%end_row(i))
end do

if (associated(ls%mass)) then
   write(outunit,"(A)")  "  mass matrix values"
   do i=1,ls%neq_vert + ls%neq_edge + ls%neq_face
      write(outunit,"(SS,1P,4E18.10E2)")  ls%mass(ls%begin_row(i):ls%end_row(i))
      write(outunit,"(4I21)")  ls%column_index(ls%begin_row(i):ls%end_row(i))
   end do
endif

write(outunit,"(A)")  "  right hand side"
write(outunit,"(SS,1P,4E18.10E2)")  ls%rhs

write(outunit,"(A)")  "  solution(0)"
write(outunit,"(SS,1P,E18.10E2)")  ls%solution(0)
write(outunit,"(A)")  "  solution vector"
write(outunit,"(SS,1P,4E18.10E2)")  ls%solution(1:)

write(outunit,"(A)")  "  equation type, key is"
write(outunit,"(A,I2)")  "  INTERIOR ",INTERIOR
write(outunit,"(A,I2)")  "  DIRICHLET ",DIRICHLET
write(outunit,"(A,I2)")  "  NATURAL ",NATURAL
write(outunit,"(A,I2)")  "  MIXED ",MIXED
write(outunit,"(A,I2)")  "  PERIODIC ",PERIODIC_MASTER
write(outunit,"(7I11)")  ls%equation_type

return
end subroutine print_linear_system

!          ----------------
subroutine print_error_info(grid,procs,io_control,still_sequential,this_time, &
                            tag,title,reduction,errest)
!          ----------------

!----------------------------------------------------
! This routine prints norms of the error
!----------------------------------------------------

implicit none

!----------------------------------------------------
! Dummy arguments

type (grid_type), intent(inout) :: grid
type (proc_info), intent(in) :: procs
type(io_options), intent(in) :: io_control
logical, intent(in) :: still_sequential
integer, intent(in) :: this_time(:)
integer, intent(in) :: tag
character (len=*), optional, intent(in) :: title
logical, optional, intent(in) :: reduction
integer, intent(in), optional :: errest
!----------------------------------------------------

!----------------------------------------------------
! Local variables:

integer :: no_ints(1),ni,nr,proc,i,np,allocstat,who,when,comp,nevec,evec, &
           nsoln,soln
integer :: ints(4)
real (my_real), allocatable :: send_real(:)
integer, pointer :: recv_int(:)
real (my_real), pointer :: recv_real(:)
real (my_real), save, allocatable :: old_energy_error(:), old_Linf_error(:), &
                                     old_L2_error(:), old_energy_errest(:), &
                                     old_Linf_errest(:), old_L2_errest(:)
real (my_real), allocatable :: new_energy_error(:), new_Linf_error(:), &
                               new_L2_error(:), new_energy_errest(:), &
                               new_Linf_errest(:), new_L2_errest(:), &
                               normte(:), normtinf(:), normt2(:), &
                               normse(:), normsinf(:), norms2(:)
real (my_real), allocatable :: iold_energy_error(:,:), iold_Linf_error(:,:), &
                               iold_L2_error(:,:), iold_energy_errest(:,:), &
                               iold_Linf_errest(:,:), iold_L2_errest(:,:), &
                               inew_energy_error(:,:), inew_Linf_error(:,:), &
                               inew_L2_error(:,:), inew_energy_errest(:,:), &
                               inew_Linf_errest(:,:), inew_L2_errest(:,:), &
                               inormte(:,:), inormtinf(:,:), inormt2(:,:), &
                               inormse(:,:), inormsinf(:,:), inorms2(:,:)
logical :: true_provided, printit, printall, print_energy_error, &
           print_Linf_error, print_L2_error, print_energy_errest, &
           print_Linf_errest, print_L2_errest
character(len=22) :: fmt1, fmt2
real, pointer :: times(:,:)
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

! If this is not the right time to print, return

who = io_control%print_error_who
when = io_control%print_error_when

if (.not. any(this_time == when) .or. who == NO_ONE) return

! If I'm the master and only slaves print the error, return

if (my_proc(procs) == MASTER .and. who == SLAVES) return

! stop the clocks

call pause_watch(all_watches)

nevec = max(1,grid%num_eval)
nsoln = grid%nsoln
if (nevec*grid%system_size /= nsoln) then
   call warning("nsoln is not neval*syssize; printing errors for first solution only", &
                intlist=(/nsoln,nevec,grid%system_size/))
   nevec = 1
endif
np = num_proc(procs)
fmt1 = "(SS,1P,A,E18.10E2)"
fmt2 = "(SS,1P,A, 128E18.10E2)"

! make space for holding errors from one call to another

if (allocated(old_energy_error)) then
   if (size(old_energy_error) /= nevec) then
      deallocate(old_energy_error,old_Linf_error,old_L2_error, &
                 old_energy_errest,old_Linf_errest,old_L2_errest, &
                 stat=allocstat)
   endif
endif
if (.not. allocated(old_energy_error)) then
   allocate(old_energy_error(nevec),old_Linf_error(nsoln),old_L2_error(nsoln), &
            old_energy_errest(nevec),old_Linf_errest(nsoln), &
            old_L2_errest(nsoln),stat=allocstat)
   if (allocstat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in print_error_info",procs=procs)
      return
   endif
   old_energy_error = 0.0_my_real
   old_Linf_error = 0.0_my_real
   old_L2_error = 0.0_my_real
   old_energy_errest = 0.0_my_real
   old_Linf_errest = 0.0_my_real
   old_L2_errest = 0.0_my_real
endif

allocate(new_energy_error(nevec),new_Linf_error(nsoln),new_L2_error(nsoln), &
      new_energy_errest(nevec),new_Linf_errest(nsoln),new_L2_errest(nsoln), &
      normte(nevec), normtinf(nsoln), normt2(nsoln), &
      normse(nevec), normsinf(nsoln), norms2(nsoln), &
      stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in print_error_info",procs=procs)
   return
endif

new_energy_error = 1.0_my_real
new_Linf_error = 1.0_my_real
new_L2_error = 1.0_my_real
new_energy_errest = 0.0_my_real
new_Linf_errest = 0.0_my_real
new_L2_errest = 0.0_my_real
normte = 1.0_my_real
normtinf = 1.0_my_real
normt2 = 1.0_my_real
normse = 1.0_my_real
normsinf = 1.0_my_real
norms2 = 1.0_my_real

! compute all the needed errors and error estimates

if (my_proc(procs) /= MASTER) then
   soln = 0
   do evec=1,nevec
      do comp=1,grid%system_size
         soln = soln + 1
         true_provided=(trues(0.5_my_real,0.5_my_real,comp,evec)/=huge(0.0_my_real))
         if (true_provided) then
            select case (io_control%print_error_what)
            case (ENERGY_ERR)
               if (comp==1) then
                  call norm_error(grid,procs,still_sequential,comp,evec, &
                                  my_energy_norm=new_energy_error(evec))
                  if (grid%errtype == RELATIVE_ERROR) then
                     call norm_true(grid,procs,still_sequential,comp,evec, &
                                    energy=normte(evec))
                  endif
               endif
            case (ENERGY_LINF_ERR, ENERGY_LINF_L2_ERR)
               if (comp==1) then
                  call norm_error(grid,procs,still_sequential,comp,evec, &
                                  my_energy_norm=new_energy_error(evec), &
                                  my_Linf_norm=new_Linf_error(soln), &
                                  my_L2_norm=new_L2_error(soln))
                  if (grid%errtype == RELATIVE_ERROR) then
                     call norm_true(grid,procs,still_sequential,comp,evec, &
                                    linf=normtinf(soln),l2=normt2(soln), &
                                    energy=normte(evec))
                  endif
               else
                  call norm_error(grid,procs,still_sequential,comp,evec, &
                                  my_Linf_norm=new_Linf_error(soln), &
                                  my_L2_norm=new_L2_error(soln))
                  if (grid%errtype == RELATIVE_ERROR) then
                     call norm_true(grid,procs,still_sequential,comp,evec, &
                                    linf=normtinf(soln),l2=normt2(soln))
                  endif
               endif
            case (ENERGY_L2_ERR)
               if (comp==1) then
                  call norm_error(grid,procs,still_sequential,comp,evec, &
                                  my_energy_norm=new_energy_error(evec), &
                                  my_L2_norm=new_L2_error(soln))
                  if (grid%errtype == RELATIVE_ERROR) then
                     call norm_true(grid,procs,still_sequential,comp,evec, &
                                    l2=normt2(soln),energy=normte(evec))
                  endif
               else
                  call norm_error(grid,procs,still_sequential,comp,evec, &
                                  my_L2_norm=new_L2_error(soln))
                  if (grid%errtype == RELATIVE_ERROR) then
                     call norm_true(grid,procs,still_sequential,comp,evec, &
                                    l2=normt2(soln))
                  endif
               endif
            case (LINF_ERR, LINF_L2_ERR)
               call norm_error(grid,procs,still_sequential,comp,evec, &
                               my_Linf_norm=new_Linf_error(soln), &
                               my_L2_norm=new_L2_error(soln))
                  if (grid%errtype == RELATIVE_ERROR) then
                     call norm_true(grid,procs,still_sequential,comp,evec, &
                                    linf=normtinf(soln),l2=normt2(soln))
                  endif
            case (L2_ERR)
               call norm_error(grid,procs,still_sequential,comp,evec, &
                               my_L2_norm=new_L2_error(soln))
                  if (grid%errtype == RELATIVE_ERROR) then
                     call norm_true(grid,procs,still_sequential,comp,evec, &
                                    l2=normt2(soln))
                  endif
            end select
         endif
         if (present(errest)) then
            call error_estimate(grid,procs,errest,soln, &
                                errest_Linf=new_Linf_errest(soln), &
                                errest_L2=new_L2_errest(soln))
            if (grid%errtype == RELATIVE_ERROR) then
               call norm_solution(grid,procs,still_sequential,comp,evec, &
                                  linf=normsinf(soln),l2=norms2(soln))
            endif
            if (comp==1) then
               call error_estimate(grid,procs,errest,evec, &
                                   errest_energy=new_energy_errest(evec))
               if (grid%errtype == RELATIVE_ERROR) then
                  call norm_solution(grid,procs,still_sequential,comp,evec, &
                                     energy=normse(evec))
               endif
            endif
         endif
      end do
   end do
endif

! If the master will print, get it the info

if (who /= SLAVES) then
   if (my_proc(procs) /= MASTER) then
      allocate(send_real(12*nsoln+6*nevec),stat=allocstat)
      if (allocstat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in print_error_info",procs=procs)
         return
      endif
      send_real(                 1:   nsoln        ) = new_Linf_error
      send_real(   nsoln+        1: 2*nsoln        ) = old_Linf_error
      send_real( 2*nsoln+        1: 3*nsoln        ) = new_Linf_errest
      send_real( 3*nsoln+        1: 4*nsoln        ) = old_Linf_errest
      send_real( 4*nsoln+        1: 5*nsoln        ) = new_L2_error
      send_real( 5*nsoln+        1: 6*nsoln        ) = old_L2_error
      send_real( 6*nsoln+        1: 7*nsoln        ) = new_L2_errest
      send_real( 7*nsoln+        1: 8*nsoln        ) = old_L2_errest
      send_real( 8*nsoln+        1: 9*nsoln        ) = normtinf
      send_real( 9*nsoln+        1:10*nsoln        ) = normt2
      send_real(10*nsoln+        1:11*nsoln        ) = normsinf
      send_real(11*nsoln+        1:12*nsoln        ) = norms2
      send_real(12*nsoln+        1:12*nsoln+  nevec) = new_energy_error
      send_real(12*nsoln+  nevec+1:12*nsoln+2*nevec) = old_energy_error
      send_real(12*nsoln+2*nevec+1:12*nsoln+3*nevec) = new_energy_errest
      send_real(12*nsoln+3*nevec+1:12*nsoln+4*nevec) = old_energy_errest
      send_real(12*nsoln+4*nevec+1:12*nsoln+5*nevec) = normte
      send_real(12*nsoln+5*nevec+1:12*nsoln+6*nevec) = normse
      ni = 0
      nr = 12*nsoln+6*nevec
      call phaml_send(procs,MASTER,no_ints,ni,send_real,nr,tag)
      call phaml_recv(procs,proc,recv_int,ni,recv_real,nr,tag+2)
      if (associated(recv_int)) then
         save_convergence = recv_int(1)==1
         deallocate(recv_int,stat=allocstat)
      endif
      if (save_convergence) then
         call get_grid_info(grid,procs,still_sequential,2929,total_dof=ints(1),&
                            mindeg=ints(2),maxdeg=ints(3),nlev=ints(4))
         call read_watch(times,(/trefine,tassemble,tsolve,ttotal/),(/"wall"/))
         send_real(1:4) = times(1:4,1)
         deallocate(times)
         call phaml_send(procs,MASTER,ints,4,send_real,4,tag+1)
      endif
      deallocate(send_real,stat=allocstat)
      if (allocstat /= 0) then
         call warning("deallocation of messge failed in print_error_info")
      endif
   else ! MASTER
      allocate(iold_energy_error(np,nevec),iold_Linf_error(np,nsoln), &
               iold_L2_error(np,nsoln),iold_energy_errest(np,nevec), &
               iold_Linf_errest(np,nsoln),iold_L2_errest(np,nsoln), &
               inew_energy_error(np,nevec),inew_Linf_error(np,nsoln), &
               inew_L2_error(np,nsoln),inew_energy_errest(np,nevec), &
               inew_Linf_errest(np,nsoln),inew_L2_errest(np,nsoln), &
               inormte(np,nevec),inormtinf(np,nsoln),inormt2(np,nsoln), &
               inormse(np,nevec),inormsinf(np,nsoln),inorms2(np,nsoln), &
               stat=allocstat)
      if (allocstat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in print_error_info",procs=procs)
         return
      endif
      do i=1,np
         call phaml_recv(procs,proc,recv_int,ni,recv_real,nr,tag)
         if (save_convergence) then
            call phaml_send(procs,proc,(/1/),1,(/0.0_my_real/),0,tag+2)
         else
            call phaml_send(procs,proc,(/0/),1,(/0.0_my_real/),0,tag+2)
         endif
         inew_Linf_error(proc,:) = &
                           recv_real(                 1:   nsoln        )
         iold_Linf_error(proc,:) = &
                           recv_real(   nsoln+        1: 2*nsoln        )
         inew_Linf_errest(proc,:) = &
                           recv_real( 2*nsoln+        1: 3*nsoln        )
         iold_Linf_errest(proc,:) = &
                           recv_real( 3*nsoln+        1: 4*nsoln        )
         inew_L2_error(proc,:) = &
                           recv_real( 4*nsoln+        1: 5*nsoln        )
         iold_L2_error(proc,:) = &
                           recv_real( 5*nsoln+        1: 6*nsoln        )
         inew_L2_errest(proc,:) = &
                           recv_real( 6*nsoln+        1: 7*nsoln        )
         iold_L2_errest(proc,:) = &
                           recv_real( 7*nsoln+        1: 8*nsoln        )
         inormtinf(proc,:) = &
                           recv_real( 8*nsoln+        1: 9*nsoln        )
         inormt2(proc,:) = &
                           recv_real( 9*nsoln+        1:10*nsoln        )
         inormsinf(proc,:) = &
                           recv_real(10*nsoln+        1:11*nsoln        )
         inorms2(proc,:) = &
                           recv_real(11*nsoln+        1:12*nsoln        )
         inew_energy_error(proc,:) = &
                           recv_real(12*nsoln+        1:12*nsoln+  nevec)
         iold_energy_error(proc,:) = &
                           recv_real(12*nsoln+  nevec+1:12*nsoln+2*nevec)
         inew_energy_errest(proc,:) = &
                           recv_real(12*nsoln+2*nevec+1:12*nsoln+3*nevec)
         iold_energy_errest(proc,:) = &
                           recv_real(12*nsoln+3*nevec+1:12*nsoln+4*nevec)
         inormte(proc,:) = &
                           recv_real(12*nsoln+4*nevec+1:12*nsoln+5*nevec)
         inormse(proc,:) = &
                           recv_real(12*nsoln+5*nevec+1:12*nsoln+6*nevec)
         deallocate(recv_real,stat=allocstat)
         if (allocstat /= 0) then
            call warning("deallocation of messge failed in print_error_info")
         endif
      end do
      if (save_convergence) then
         ints(2) = huge(0)
         ints(3) = 0
         ints(4) = 0
         allocate(times(4,1))
         times = 0
         do i=1,np
            call phaml_recv(procs,proc,recv_int,ni,recv_real,nr,tag+1)
            ints(1) = recv_int(1)
            ints(2) = min(ints(2),recv_int(2))
            ints(3) = max(ints(3),recv_int(3))
            ints(4) = max(ints(4),recv_int(4))
            times(1,1) = max(real(times(1,1),my_real),recv_real(1))
            times(2,1) = max(real(times(2,1),my_real),recv_real(2))
            times(3,1) = max(real(times(3,1),my_real),recv_real(3))
            times(4,1) = max(real(times(4,1),my_real),recv_real(4))
            deallocate(recv_int,recv_real)
         end do
      endif
      new_Linf_error = maxval(inew_Linf_error,dim=1)
      old_Linf_error = maxval(iold_Linf_error,dim=1)
      new_Linf_errest = maxval(inew_Linf_errest,dim=1)
      old_Linf_errest = maxval(iold_Linf_errest,dim=1)
      normtinf = maxval(inormtinf,dim=1)
      normsinf = maxval(inormsinf,dim=1)
      if (still_sequential) then
         new_L2_error = inew_L2_error(1,:)
         old_L2_error = iold_L2_error(1,:)
         new_L2_errest = inew_L2_errest(1,:)
         old_L2_errest = iold_L2_errest(1,:)
         normt2 = inormt2(1,:)
         norms2 = inorms2(1,:)
         new_energy_error = inew_energy_error(1,:)
         old_energy_error = iold_energy_error(1,:)
         new_energy_errest = inew_energy_errest(1,:)
         old_energy_errest = iold_energy_errest(1,:)
         normte = inormte(1,:)
         normse = inormte(1,:)
      else
         new_L2_error = sqrt(sum(inew_L2_error**2,dim=1))
         old_L2_error = sqrt(sum(iold_L2_error**2,dim=1))
         new_L2_errest = sqrt(sum(inew_L2_errest**2,dim=1))
         old_L2_errest = sqrt(sum(iold_L2_errest**2,dim=1))
         normt2 = sqrt(sum(inormt2**2,dim=1))
         norms2 = sqrt(sum(inorms2**2,dim=1))
         new_energy_error = sqrt(sum(inew_energy_error**2,dim=1))
         old_energy_error = sqrt(sum(iold_energy_error**2,dim=1))
         new_energy_errest = sqrt(sum(inew_energy_errest**2,dim=1))
         old_energy_errest = sqrt(sum(iold_energy_errest**2,dim=1))
         normte = sqrt(sum(inormte**2,dim=1))
         normse = sqrt(sum(inormse**2,dim=1))
      endif
   endif
endif

! change to relative error if so requested.  Note that slaves divide their
! part of the error by their part of the true solution, which is not right

if (grid%errtype == RELATIVE_ERROR) then
   where (normte /= 0.0_my_real) new_energy_error = new_energy_error/normte
   where (normtinf /= 0.0_my_real) new_Linf_error = new_Linf_error/normtinf
   where (normt2 /= 0.0_my_real) new_L2_error = new_L2_error/normt2
   where (normse /= 0.0_my_real) new_energy_errest = new_energy_errest/normse
   where (normsinf /= 0.0_my_real) new_Linf_errest = new_Linf_errest/normsinf
   where (norms2 /= 0.0_my_real) new_L2_errest = new_L2_errest/norms2
   if (my_proc(procs) == MASTER) then
      where (inormte /= 0.0_my_real) inew_energy_error = inew_energy_error/inormte
      where (inormtinf /= 0.0_my_real) inew_Linf_error = inew_Linf_error/inormtinf
      where (inormt2 /= 0.0_my_real) inew_L2_error = inew_L2_error/inormt2
      where (inormse /= 0.0_my_real) inew_energy_errest = inew_energy_errest/inormse
      where (inormsinf /= 0.0_my_real) inew_Linf_errest = inew_Linf_errest/inormsinf
      where (inorms2 /= 0.0_my_real) inew_L2_errest = inew_L2_errest/inorms2
   endif
endif

if (PARALLEL==SEQUENTIAL .or. my_proc(procs) == MASTER) then
   if (PARALLEL==SEQUENTIAL .and. save_convergence) then
      call get_grid_info(grid,procs,still_sequential,2929,total_dof=ints(1),&
                         mindeg=ints(2),maxdeg=ints(3),nlev=ints(4))
      call read_watch(times,(/trefine,tassemble,tsolve,ttotal/),(/"wall"/))
   endif
   if (save_convergence .and. ints(1) /= last_dof) then
      if (first_conv) then
         write(convfileunit,"(A)") "#"
         write(convfileunit,"(A)") "# Columns are: ndof, energy err, energy errest, Linf err, Linf errest, L2 err,"
         write(convfileunit,"(A)") "#              L2 errest, nlev, mindeg, maxdeg, refine time, discretize time,"
         write(convfileunit,"(A)") "#              solve time, total time"
         write(convfileunit,"(A)") "#"
         first_conv = .false.
      endif
      write(convfileunit,"(I11,SS,1P,6E23.15E2,3I11,4E15.7E2)") ints(1), &
                  new_energy_error(1), new_energy_errest(1), &
                  new_Linf_error(1), new_Linf_errest(1), new_L2_error(1), &
                  new_L2_errest(1), ints(4), ints(2), ints(3), times
      deallocate(times)
      call my_flush(convfileunit)
      last_dof = ints(1)
   endif
endif

! determine whether or not this processor prints

printit = ((my_proc(procs)==MASTER .and. (who==EVERYONE .or. who==MASTER)) &
         .or.                                                              &
           (my_proc(procs)/=MASTER .and. (who==EVERYONE .or. who==SLAVES)))

printall = (my_proc(procs)==MASTER .and. who==MASTER_ALL)

print_energy_error = io_control%print_error_what == ENERGY_ERR .or. &
                     io_control%print_error_what == ENERGY_LINF_ERR .or. &
                     io_control%print_error_what == ENERGY_LINF_L2_ERR .or. &
                     io_control%print_error_what == ENERGY_L2_ERR
print_Linf_error = io_control%print_error_what == ENERGY_LINF_ERR .or. &
                   io_control%print_error_what == ENERGY_LINF_L2_ERR .or. &
                   io_control%print_error_what == LINF_L2_ERR .or. &
                   io_control%print_error_what == LINF_ERR
print_L2_error = io_control%print_error_what == ENERGY_L2_ERR .or. &
                 io_control%print_error_what == ENERGY_LINF_L2_ERR .or. &
                 io_control%print_error_what == LINF_L2_ERR .or. &
                 io_control%print_error_what == L2_ERR
print_energy_errest = io_control%print_errest_what == ENERGY_ERREST .or. &
                 io_control%print_errest_what == ENERGY_LINF_ERREST .or. &
                 io_control%print_errest_what == ENERGY_LINF_L2_ERREST .or. &
                 io_control%print_errest_what == ENERGY_L2_ERREST
print_Linf_errest = io_control%print_errest_what == LINF_ERREST .or. &
                 io_control%print_errest_what == ENERGY_LINF_ERREST .or. &
                 io_control%print_errest_what == ENERGY_LINF_L2_ERREST .or. &
                 io_control%print_errest_what == LINF_L2_ERREST
print_L2_errest = io_control%print_errest_what == L2_ERREST .or. &
                 io_control%print_errest_what == ENERGY_L2_ERREST .or. &
                 io_control%print_errest_what == ENERGY_LINF_L2_ERREST .or. &
                 io_control%print_errest_what == LINF_L2_ERREST

! print the errors and/or error estimates

if (printit .or. printall) then
   write(outunit,"(A)")
   if (present(title)) then
      write(outunit,"(A)") title
   else
      write(outunit,"(A)") "Norms of error:"
   endif

   soln = 0
   do evec=1,nevec
      if (associated(grid%eigenvalue)) then
         write(outunit,"(A,I11)") " Eigenvector ",evec
      endif
      true_provided=(trues(0.5_my_real,0.5_my_real,1,evec)/=huge(0.0_my_real))
      if (true_provided) then
        if (print_energy_error) then
         if (printit) then
            write(outunit,fmt1) "   energy norm of error     ",new_energy_error(evec)
         else
            write(outunit,fmt2) "   energy norm of error     ",inew_energy_error(:,evec)
         endif
         if (present(reduction)) then
            if (printit) then
               if (reduction .and. old_energy_error(evec)/=0.0_my_real) then
                  write(outunit,fmt1) "      reduced by            ",new_energy_error(evec)/old_energy_error(evec)
               endif
            else
               if (reduction .and. all(iold_energy_error(:,evec)/=0.0_my_real)) then
                  write(outunit,fmt2) "      reduced by            ",inew_energy_error(:,evec)/iold_energy_error(:,evec)
               endif
            endif
         endif
        endif
      endif
      if (print_energy_errest) then
         if (printit) then
            write(outunit,fmt1) "   energy error estimate    ",new_energy_errest(evec)
         else
            write(outunit,fmt2) "   energy error estimate    ",inew_energy_errest(:,evec)
         endif
         if (present(reduction)) then
            if (printit) then
               if (reduction .and. old_energy_errest(evec)/=0.0_my_real) then
                  write(outunit,fmt1) "      reduced by            ",new_energy_errest(evec)/old_energy_errest(evec)
               endif
            else
               if (reduction .and. all(iold_energy_errest(:,evec)/=0.0_my_real)) then
                  write(outunit,fmt2) "      reduced by            ",inew_energy_errest(:,evec)/iold_energy_errest(:,evec)
               endif
            endif
         endif
         if (true_provided) then
            if (printit) then
               if (new_energy_error(evec) /= 0.0_my_real) then
                  write(outunit,fmt1) "      effectivity           ",new_energy_errest(evec)/new_energy_error(evec)
               endif
            else
               if (all(inew_energy_error(:,evec) /= 0.0_my_real)) then
                  write(outunit,fmt2) "      effectivity           ",inew_energy_errest(:,evec)/inew_energy_error(:,evec)
               endif
            endif
         endif
      endif
      do comp=1,grid%system_size
         soln = soln + 1
         if (grid%system_size > 1) then
            write(outunit,"(A,I11)") "  Component ",comp
         endif
         true_provided=(trues(0.5_my_real,0.5_my_real,comp,evec)/=huge(0.0_my_real))
         if (true_provided) then
           if (print_Linf_error) then
            if (printit) then
               write(outunit,fmt1) "   Linf norm of error       ",new_Linf_error(soln)
            else
               write(outunit,fmt2) "   Linf norm of error       ",inew_Linf_error(:,soln)
            endif
            if (present(reduction)) then
               if (printit) then
                  if (reduction .and. old_Linf_error(soln)/=0.0_my_real) then
                     write(outunit,fmt1) "      reduced by            ",new_Linf_error(soln)/old_Linf_error(soln)
                  endif
               else
                  if (reduction .and. all(iold_Linf_error(:,soln)/=0.0_my_real)) then
                     write(outunit,fmt2) "      reduced by            ",inew_Linf_error(:,soln)/iold_Linf_error(:,soln)
                  endif
               endif
            endif
           endif
         endif
         if (print_Linf_errest) then
            if (printit) then
               write(outunit,fmt1) "   Linf error estimate      ",new_Linf_errest(soln)
            else
               write(outunit,fmt2) "   Linf error estimate      ",inew_Linf_errest(:,soln)
            endif
            if (present(reduction)) then
               if (printit) then
                  if (reduction .and. old_Linf_errest(soln)/=0.0_my_real) then
                     write(outunit,fmt1) "      reduced by            ",new_Linf_errest(soln)/old_Linf_errest(soln)
                  endif
               else
                  if (reduction .and. all(iold_Linf_errest(:,soln)/=0.0_my_real)) then
                     write(outunit,fmt2) "      reduced by            ",inew_Linf_errest(:,soln)/iold_Linf_errest(:,soln)
                  endif
               endif
            endif
            if (true_provided) then
               if (printit) then
                  if (new_Linf_error(soln) /= 0.0_my_real) then
                     write(outunit,fmt1) "      effectivity           ",new_Linf_errest(soln)/new_Linf_error(soln)
                  endif
               else
                  if (all(inew_Linf_error(:,soln) /= 0.0_my_real)) then
                     write(outunit,fmt2) "      effectivity          ",inew_Linf_errest(:,soln)/inew_Linf_error(:,soln)
                  endif
               endif
            endif
         endif
         if (true_provided) then
           if (print_L2_error) then
            if (printit) then
               write(outunit,fmt1) "   L2 norm of error         ",new_L2_error(soln)
            else
               write(outunit,fmt2) "   L2 norm of error         ",inew_L2_error(:,soln)
            endif
            if (present(reduction)) then
               if (printit) then
                  if (reduction .and. old_L2_error(soln)/=0.0_my_real) then
                     write(outunit,fmt1) "      reduced by            ",new_L2_error(soln)/old_L2_error(soln)
                  endif
               else
                  if (reduction .and. all(iold_L2_error(:,soln)/=0.0_my_real)) then
                     write(outunit,fmt2) "      reduced by            ",inew_L2_error(:,soln)/iold_L2_error(:,soln)
                  endif
               endif
            endif
           endif
         endif
         if (print_L2_errest) then
            if (printit) then
               write(outunit,fmt1) "   L2 error estimate        ",new_L2_errest(soln)
            else
               write(outunit,fmt2) "   L2 error estimate        ",inew_L2_errest(:,soln)
            endif
            if (present(reduction)) then
               if (printit) then
                  if (reduction .and. old_L2_errest(soln)/=0.0_my_real) then
                     write(outunit,fmt1) "      reduced by            ",new_L2_errest(soln)/old_L2_errest(soln)
                  endif
               else
                  if (reduction .and. all(iold_L2_errest(:,soln)/=0.0_my_real)) then
                     write(outunit,fmt2) "      reduced by            ",inew_L2_errest(:,soln)/iold_L2_errest(:,soln)
                  endif
               endif
            endif
            if (true_provided) then
               if (printit) then
                  if (new_L2_error(soln) /= 0.0_my_real) then
                     write(outunit,fmt1) "      effectivity           ",new_L2_errest(soln)/new_L2_error(soln)
                  endif
               else
                  if (all(inew_L2_error(:,soln) /= 0.0_my_real)) then
                     write(outunit,fmt2) "      effectivity          ",inew_L2_errest(:,soln)/inew_L2_error(:,soln)
                  endif
               endif
            endif
         endif
      end do
   end do
endif

old_energy_error = new_energy_error
old_energy_errest = new_energy_errest
old_Linf_error = new_Linf_error
old_Linf_errest = new_Linf_errest
old_L2_error = new_L2_error
old_L2_errest = new_L2_errest

deallocate(new_energy_error,new_energy_errest, new_Linf_error,new_Linf_errest, &
           new_L2_error,new_L2_errest,normte,normtinf,normt2,normse,normsinf, &
           norms2,stat=allocstat)
if (allocstat /= 0) then
   call warning("deallocation failed in print_error_info")
endif
if (allocated(inew_energy_error)) then
   deallocate(inew_energy_error,inew_energy_errest, inew_Linf_error, &
              inew_Linf_errest, inew_L2_error,inew_L2_errest, &
              iold_energy_error,iold_energy_errest, iold_Linf_error, &
              iold_Linf_errest, iold_L2_error,iold_L2_errest,inormte,inormtinf,&
              inormt2,inormse,inormsinf,inorms2,stat=allocstat)
   if (allocstat /= 0) then
      call warning("deallocation failed in print_error_info")
   endif
endif

call end_pause_watch(all_watches)

end subroutine print_error_info

!          -----------------
subroutine print_linsys_info(linsys,procs,io_control,still_sequential, &
                             this_time,tag)
!          -----------------

!----------------------------------------------------
! This routine prints information about the linear system
!----------------------------------------------------
 
implicit none
 
!----------------------------------------------------
! Dummy arguments
 
type (linsys_type), intent(in) :: linsys
type (proc_info), intent(in) :: procs
type (io_options), intent(in) :: io_control
logical, intent(in) :: still_sequential
integer, intent(in) :: this_time(:),tag
!----------------------------------------------------
 
!----------------------------------------------------
! Local variables:
 
integer :: i,proc,np,who,when,astat,lev,eq,col
integer :: neq,neq_own,nonzero,nonzero_own,coarse_neq,coarse_bandwidth, &
           lapack_neq,lapack_bandwidth,neq_cond,neq_cond_own,nonzero_cond, &
           nonzero_cond_own
integer, allocatable :: ineq(:),ineq_own(:),inonzero(:),inonzero_own(:), &
                        ineq_cond(:),ineq_cond_own(:),inonzero_cond(:), &
                        inonzero_cond_own(:)
integer :: send_int(12),ni,nr
real (my_real) :: no_reals(1)
integer, pointer :: recv_int(:)
real (my_real), pointer :: recv_real(:)
 
!----------------------------------------------------
 
!----------------------------------------------------
! Begin executable code
 
! If this is not the right time to print, return

who = io_control%print_linsys_who
when = io_control%print_linsys_when

if (.not. any(this_time == when) .or. who == NO_ONE) return

! If I'm the master and only slaves print, return

if (my_proc(procs) == MASTER .and. who == SLAVES) return

! stop the clocks

call pause_watch(all_watches)

! If I'm not the master, collect my linear system info

if (my_proc(procs) /= MASTER) then
   neq = linsys%neq_vert + linsys%neq_edge + linsys%neq_face
   neq_own = 0
   nonzero = 0
   nonzero_own = 0
   neq_cond = linsys%neq_vert + linsys%neq_edge
   neq_cond_own = 0
   nonzero_cond = 0
   nonzero_cond_own = 0
   do lev=1,linsys%nlev+2
      do eq=linsys%begin_level(lev),linsys%begin_level(lev+1)-1
         if (linsys%iown(eq)) then
            neq_own = neq_own + 1
            if (lev < linsys%nlev+2) neq_cond_own = neq_cond_own + 1
         endif
         do col=linsys%begin_row(eq),linsys%end_row_face(eq)
            if (linsys%column_index(col) == NO_ENTRY) cycle
            nonzero = nonzero + 1
            if (linsys%iown(eq)) nonzero_own = nonzero_own + 1
            if (lev < linsys%nlev+2 .and. col <= linsys%end_row_edge(eq)) then
               nonzero_cond = nonzero_cond + 1
               if (linsys%iown(eq)) nonzero_cond_own = nonzero_cond_own + 1
            endif
         end do
      end do
   end do
   if (linsys%coarse_band_exists) then
      coarse_neq = linsys%coarse_matrix%neq
      coarse_bandwidth = linsys%coarse_matrix%halfbandwidth
   else
      coarse_neq = -1
      coarse_bandwidth = -1
   endif
   if (linsys%lapack_gen_band_exists .or. &
       linsys%lapack_symm_band_exists) then
      lapack_neq = linsys%lapack_mat%neq
      lapack_bandwidth = linsys%lapack_mat%halfbandwidth
   else
      lapack_neq = -1
      lapack_bandwidth = -1
   endif
endif

! If the master will print, get it the info

if (who == MASTER .or. who == MASTER_ALL .or. who == EVERYONE) then
   if (my_proc(procs) /= MASTER) then
      send_int(1) = neq
      send_int(2) = neq_own
      send_int(3) = nonzero
      send_int(4) = nonzero_own
      send_int(5) = coarse_neq
      send_int(6) = coarse_bandwidth
      send_int(7) = lapack_neq
      send_int(8) = lapack_bandwidth
      send_int(9) = neq_cond
      send_int(10) = neq_cond_own
      send_int(11) = nonzero_cond
      send_int(12) = nonzero_cond_own
      ni = 12
      nr = 0
      call phaml_send(procs,MASTER,send_int,ni,no_reals,nr,tag)
   else
      np = num_proc(procs)
      allocate(ineq(np),ineq_own(np),inonzero(np),inonzero_own(np), &
               ineq_cond(np),ineq_cond_own(np),inonzero_cond(np), &
               inonzero_cond_own(np),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in print_linsys_info",procs=procs)
         return
      endif
      do i=1,np
         call phaml_recv(procs,proc,recv_int,ni,recv_real,nr,tag)
         ineq(proc) = recv_int(1)
         ineq_own(proc) = recv_int(2)
         inonzero(proc) = recv_int(3)
         inonzero_own(proc) = recv_int(4)
         coarse_neq = recv_int(5)
         coarse_bandwidth = recv_int(6)
         lapack_neq = recv_int(7)
         lapack_bandwidth = recv_int(8)
         ineq_cond(proc) = recv_int(9)
         ineq_cond_own(proc) = recv_int(10)
         inonzero_cond(proc) = recv_int(11)
         inonzero_cond_own(proc) = recv_int(12)
         deallocate(recv_int,stat=astat)
      end do
      if (still_sequential) then
         neq = ineq(1)
         nonzero = inonzero(1)
         neq_cond = ineq_cond(1)
         nonzero_cond = inonzero_cond(1)
      else
         neq = sum(ineq_own)
         nonzero = sum(inonzero_own)
         neq_cond = sum(ineq_cond_own)
         nonzero_cond = sum(inonzero_cond_own)
      endif
   endif
endif

! print the info, if requested

if (my_proc(procs) /= MASTER) then
 if (who == SLAVES .or. who == EVERYONE) then
   write(outunit,"(A)")
   write(outunit,"(A)") 'My Linear System:'
   write(outunit,"(A,I11)") '   number of equations = ',neq
   write(outunit,"(A,I11)") '   number of equations I own = ',neq_own
   write(outunit,"(A,I11)") '   number of nonzeroes = ',nonzero
   write(outunit,"(A,I11)") '   condensed matrix number of equations = ',neq_cond
   write(outunit,"(A,I11)") '   condensed matrix number of equations I own = ',neq_cond_own
   write(outunit,"(A,I11)") '   condensed matrix number of nonzeroes = ',nonzero_cond
   if (coarse_neq /= -1) then
      write(outunit,"(A,I11)") '   coarse matrix neq = ',coarse_neq
      write(outunit,"(A,I11)") '   coarse (half) bandwidth = ',coarse_bandwidth
   endif
   if (lapack_neq /= -1) then
      write(outunit,"(A,I11)") '   lapack matrix neq = ',lapack_neq
      write(outunit,"(A,I11)") '   lapack (half) bandwidth = ',lapack_bandwidth
   endif
 endif
else
 if (who == MASTER_ALL) then
   write(outunit,"(A)")
   write(outunit,"(A)") 'Individual Linear Systems:'
   write(outunit,"(A,128I11)") '   number of equations = ',ineq
   write(outunit,"(A,128I11)") '   number of equations I own = ',ineq_own
   write(outunit,"(A,128I11)") '   number of nonzeroes = ',inonzero
   write(outunit,"(A,128I11)") '   condensed matrix number of equations = ',ineq_cond
   write(outunit,"(A,128I11)") '   condensed matrix number of equations I own = ',ineq_cond_own
   write(outunit,"(A,128I11)") '   condensed matrix number of nonzeroes = ',inonzero_cond
   if (coarse_neq /= -1) then
      write(outunit,"(A,128I11)") '   coarse matrix neq = ',coarse_neq
      write(outunit,"(A,128I11)") '   coarse (half) bandwidth = ',coarse_bandwidth
   endif
   if (lapack_neq /= -1) then
      write(outunit,"(A,128I11)") '   lapack matrix neq = ',lapack_neq
      write(outunit,"(A,128I11)") '   lapack (half) bandwidth = ',lapack_bandwidth
   endif
 endif
 if (who == MASTER .or. who == EVERYONE .or. who == MASTER_ALL) then
   write(outunit,"(A)")
   write(outunit,"(A)") 'Total Linear System:'
   write(outunit,"(A,I11)") '   number of equations = ',neq
   write(outunit,"(A,I11)") '   number of nonzeroes = ',nonzero
   write(outunit,"(A,I11)") '   condensed matrix number of equations = ',neq_cond
   write(outunit,"(A,I11)") '   condensed matrix number of nonzeroes = ',nonzero_cond
   if (coarse_neq /= -1) then
      write(outunit,"(A,I11)") '   coarse matrix neq = ',coarse_neq
      write(outunit,"(A,I11)") '   coarse (half) bandwidth = ',coarse_bandwidth
   endif
   if (lapack_neq /= -1) then
      write(outunit,"(A,I11)") '   lapack matrix neq = ',lapack_neq
      write(outunit,"(A,I11)") '   lapack (half) bandwidth = ',lapack_bandwidth
   endif
 endif
endif

if (allocated(ineq)) then
   deallocate(ineq,ineq_own,inonzero,inonzero_own,stat=astat)
endif

call end_pause_watch(all_watches)

end subroutine print_linsys_info

!----------------------------------------------------------------
! CODE FOR COMPUTING ERRORS
!----------------------------------------------------------------

!          ----------
subroutine norm_error(grid,procs,still_sequential,comp,eigen, &
                      my_energy_norm,my_Linf_norm,my_L2_norm,energy_norm, &
                      Linf_norm,L2_norm)
!          ----------

!----------------------------------------------------
! This routine measures the error (defined as solution-exact)
! in the requested norms for the requested component of the solution.
! If the energy norm is requested, the component should be the beginning
! of an eigensolution; if not then the eigensolution containing that component
! is used.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments
 
type (grid_type), intent(in) :: grid
type (proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential
integer, intent(in) :: comp,eigen
real (my_real), optional, intent(out) :: my_energy_norm, my_Linf_norm, &
                                     my_L2_norm, energy_norm, Linf_norm, L2_norm
!----------------------------------------------------
 
!----------------------------------------------------
! Local variables:
 
real (my_real), pointer :: qw(:),xq(:),yq(:)
real (my_real), allocatable :: u(:,:,:),ux(:,:,:),uy(:,:,:)
real(my_real) :: cxx(grid%system_size,grid%system_size), &
                 cxy(grid%system_size,grid%system_size), &
                 cyy(grid%system_size,grid%system_size), &
                 cx(grid%system_size,grid%system_size),  &
                 cy(grid%system_size,grid%system_size),  &
                 c(grid%system_size,grid%system_size),   &
                 rs(grid%system_size)
integer :: lev,vert,elem,astat,i,j,k,nqp,ss
real (my_real) :: loc_Linf_norm, loc_L2_norm, loc_energy_norm_squared, &
                  t, xc(3), yc(3)
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

if (my_proc(procs)==MASTER) return

ss = grid%system_size

! energy norm

if (present(my_energy_norm) .or. present(energy_norm)) then

   loc_energy_norm_squared = 0.0_my_real
   do lev=1,grid%nlev
      elem = grid%head_level_elem(lev)
      do while (elem /= END_OF_LIST)
         if (.not.grid%element(elem)%isleaf .or. &
             .not.grid%element(elem)%iown) then
            elem = grid%element(elem)%next
            cycle
         endif
         xc = grid%vertex(grid%element(elem)%vertex)%coord%x
         yc = grid%vertex(grid%element(elem)%vertex)%coord%y
! quadrature order determined experimentally with u=x**10+y**10 on uniform grid
         call quadrature_rule_tri(6,xc,yc,nqp,qw,xq,yq,i,stay_in=.true.)
         allocate(u(ss,1,nqp),ux(ss,1,nqp),uy(ss,1,nqp), &
                  stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("allocation failed in norm_error",procs=procs)
            return
         endif
         call evaluate_soln_local(grid,xq,yq,elem,(/(j,j=1,ss)/),(/eigen/), &
                                  u,ux,uy)
         do i=1,nqp
            call pdecoefs(xq(i),yq(i),cxx,cxy,cyy,cx,cy,c,rs)
            do j=1,ss
               u(j,1,i) = trues(xq(i),yq(i),j,eigen) - u(j,1,i)
               ux(j,1,i) = truexs(xq(i),yq(i),j,eigen) - ux(j,1,i)
               uy(j,1,i) = trueys(xq(i),yq(i),j,eigen) - uy(j,1,i)
            end do
            do j=1,ss
               do k=1,ss
                  loc_energy_norm_squared = loc_energy_norm_squared + &
                     qw(i)*(cxx(j,k)*ux(j,1,i)*ux(k,1,i) + &
                            cyy(j,k)*uy(j,1,i)*uy(k,1,i) + &
                            cxy(j,k)*ux(j,1,i)*uy(k,1,i) + &
                             cx(j,k)*ux(j,1,i)*u(k,1,i) + &
                             cy(j,k)*uy(j,1,i)*u(k,1,i) + &
                              c(j,k)*u(j,1,i)*u(k,1,i))
               end do
            end do
         end do
         deallocate(qw,xq,yq,u,ux,uy)
         elem = grid%element(elem)%next
      end do
   end do

   loc_energy_norm_squared = abs(loc_energy_norm_squared)
   if (present(my_energy_norm)) my_energy_norm = sqrt(loc_energy_norm_squared)
   if (present(energy_norm)) then
      if (still_sequential) then
         energy_norm = sqrt(loc_energy_norm_squared)
      else
         energy_norm = phaml_global_sum(procs,loc_energy_norm_squared,1302)
         energy_norm = sqrt(energy_norm)
      endif
   endif

endif

! vertex contribution to Linf norm

loc_Linf_norm = 0.0_my_real

if (present(my_Linf_norm) .or. present(Linf_norm)) then

! go through the vertices of the grid

   do lev=1,grid%nlev
      vert = grid%head_level_vert(lev)
      do while (vert /= END_OF_LIST)
         if (grid%element(grid%vertex(vert)%assoc_elem)%iown) then
           loc_Linf_norm = max(loc_Linf_norm, &
                               abs(grid%vertex_solution(vert,comp,eigen) - &
                                   grid%vertex_exact(vert,comp,eigen)))
         endif
         vert = grid%vertex(vert)%next
      end do
   end do

endif

! L2 norm and quadrature point contributions to Linf norm

if (present(my_Linf_norm) .or. present(Linf_norm) .or. &
    present(my_L2_norm) .or. present(L2_norm)) then

   loc_L2_norm = 0.0_my_real

! go through the elements of the grid

   do lev=1,grid%nlev
      elem = grid%head_level_elem(lev)
      do while (elem /= END_OF_LIST)
         if (grid%element(elem)%isleaf .and. grid%element(elem)%iown) then

! compute the quadrature points, evaluate the solution, and compute the
! integral and max over the quadrature points

            xc = grid%vertex(grid%element(elem)%vertex)%coord%x
            yc = grid%vertex(grid%element(elem)%vertex)%coord%y
            call quadrature_rule_tri(6,xc,yc,nqp,qw,xq,yq,i,stay_in=.true.)
            allocate(u(1,1,nqp),stat=astat)
            if (astat /= 0) then
               ierr = ALLOC_FAILED
               call fatal("allocation failed in norm_error",procs=procs)
               return
            endif
            call evaluate_soln_local(grid,xq,yq,elem,(/comp/),(/eigen/),u)
            do i=1,nqp
               t = trues(xq(i),yq(i),comp,eigen)
               loc_Linf_norm = max(loc_Linf_norm,abs(u(1,1,i)-t))
               loc_L2_norm = loc_L2_norm + qw(i)*(u(1,1,i)-t)**2
            end do
            deallocate(qw,xq,yq,u)

         endif
         elem = grid%element(elem)%next
      end do
   end do

   if (present(my_Linf_norm)) my_Linf_norm = loc_Linf_norm
   if (present(Linf_norm)) then
      if (still_sequential) then
         Linf_norm = loc_Linf_norm
      else
         Linf_norm = phaml_global_max(procs,loc_Linf_norm,1302)
      endif
   endif

   if (present(my_L2_norm)) my_L2_norm = sqrt(abs(loc_L2_norm))
   if (present(L2_norm)) then
      if (still_sequential) then
         L2_norm = sqrt(abs(loc_L2_norm))
      else
         L2_norm = sqrt(abs(phaml_global_sum(procs,loc_L2_norm,1303)))
      endif
   endif

endif

end subroutine norm_error

!          -------------
subroutine norm_solution(grid,procs,still_sequential,comp,eigen, &
                         discrete_energy,linf,l2,energy)
!          -------------

!----------------------------------------------------
! This routine computes the requested norms of the requested component of the
! solution over the region owned by this processor.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments
 
type (grid_type), intent(in) :: grid
type (proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential
integer, intent(in) :: comp,eigen
real (my_real), optional, intent(out) :: discrete_energy, linf, l2, energy
!----------------------------------------------------
 
!----------------------------------------------------
! Local variables:
 
real (my_real), pointer :: qw(:),xq(:),yq(:)
real (my_real), allocatable :: u(:,:,:),ux(:,:,:),uy(:,:,:)
real(my_real) :: cxx(grid%system_size,grid%system_size), &
                 cxy(grid%system_size,grid%system_size), &
                 cyy(grid%system_size,grid%system_size), &
                 cx(grid%system_size,grid%system_size),  &
                 cy(grid%system_size,grid%system_size),  &
                 c(grid%system_size,grid%system_size),   &
                 rs(grid%system_size)
integer :: lev,vert,elem,astat,i,j,k,nqp,ss
integer :: objtype,brank,srank,objlid
real (my_real) :: xc(3), yc(3)
type(linsys_type) :: linear_system
type(solver_options) :: solver_cntl
type(io_options) :: io_cntl
real(my_real), pointer :: x(:)
real(my_real), allocatable :: y(:)
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

if (my_proc(procs)==MASTER) return

ss = grid%system_size

if (present(discrete_energy)) then

! discrete energy norm

   solver_cntl%system_size = grid%system_size
   solver_cntl%num_eval = grid%num_eval
   solver_cntl%eq_type = ELLIPTIC
   solver_cntl%lambda0 = -huge(0.0_my_real)
   solver_cntl%ignore_quad_err = .true.
   solver_cntl%inc_quad_order = 0
   io_cntl = io_options(NEVER,NO_ONE,NEVER,NO_ONE,TOO_MUCH,NO_ONE,NEVER,NEVER, &
                        NEVER,NO_ONE,PHASES,.false.)
   call create_linear_system(linear_system,grid,procs,solver_cntl,io_cntl, &
                             still_sequential,notime=.true.)
   linear_system%end_row => linear_system%end_row_face
   linear_system%neq = linear_system%neq_vert + linear_system%neq_edge + &
                       linear_system%neq_face
   linear_system%matrix_val => linear_system%stiffness
   linear_system%rhs => linear_system%rhs_nocond
   if (eigen > 1) then
      do i=1,linear_system%neq
         call eq_to_grid(linear_system,linear_system%gid(i),objtype,brank, &
                         srank,objlid,grid)
         select case (objtype)
         case (VERTEX_ID)
            linear_system%solution(i) = grid%vertex_solution(objlid,srank,eigen)
         case (EDGE_ID)
            linear_system%solution(i) = grid%edge(objlid)%solution(brank,srank,eigen)
         case (ELEMENT_ID)
            linear_system%solution(i) = grid%element(objlid)%solution(brank,srank,eigen)
         end select
      end do
   endif
   x => linear_system%solution(1:)
   allocate(y(linear_system%neq))
   call matrix_times_vector(x,y,linear_system,procs,still_sequential, &
                            1310,1311,1312,1313,1314,1315,natural=.true., &
                            notime=.true.,nocomm2=.true.)
   discrete_energy = 0.0_my_real
   do i=1,linear_system%neq
      if (linear_system%iown(i)) then
         discrete_energy = discrete_energy + x(i)*y(i)
      endif
   end do
   discrete_energy = sqrt(abs(discrete_energy))

endif

! energy norm

if (present(energy)) then

   energy = 0.0_my_real
   do lev=1,grid%nlev
      elem = grid%head_level_elem(lev)
      do while (elem /= END_OF_LIST)
         if (.not.grid%element(elem)%isleaf .or. &
             .not.grid%element(elem)%iown) then
            elem = grid%element(elem)%next
            cycle
         endif
         xc = grid%vertex(grid%element(elem)%vertex)%coord%x
         yc = grid%vertex(grid%element(elem)%vertex)%coord%y
! quadrature order determined experimentally with u=x**10+y**10 on uniform grid
         call quadrature_rule_tri(6,xc,yc,nqp,qw,xq,yq,i,stay_in=.true.)
         allocate(u(ss,1,nqp),ux(ss,1,nqp),uy(ss,1,nqp), &
                  stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("allocation failed in norm_solution",procs=procs)
            return
         endif
         call evaluate_soln_local(grid,xq,yq,elem,(/(j,j=1,ss)/),(/eigen/), &
                                  u,ux,uy)
         do i=1,nqp
            call pdecoefs(xq(i),yq(i),cxx,cxy,cyy,cx,cy,c,rs)
            do j=1,ss
               do k=1,ss
                  energy = energy + qw(i)*(cxx(j,k)*ux(j,1,i)*ux(k,1,i) + &
                                           cyy(j,k)*uy(j,1,i)*uy(k,1,i) + &
                                           cxy(j,k)*ux(j,1,i)*uy(k,1,i) + &
                                            cx(j,k)*ux(j,1,i)*u(k,1,i) + &
                                            cy(j,k)*uy(j,1,i)*u(k,1,i) + &
                                             c(j,k)*u(j,1,i)*u(k,1,i))
               end do
            end do
         end do
         deallocate(qw,xq,yq,u,ux,uy)
         elem = grid%element(elem)%next
      end do
   end do

   energy = sqrt(abs(energy))

endif

! Linf norm

if (present(Linf)) then

   Linf = 0.0_my_real

! go through the vertices of the grid

   do lev=1,grid%nlev
      vert = grid%head_level_vert(lev)
      do while (vert /= END_OF_LIST)
         if (grid%element(grid%vertex(vert)%assoc_elem)%iown) then
           Linf = max(Linf,abs(grid%vertex_solution(vert,comp,eigen)))
         endif
         vert = grid%vertex(vert)%next
      end do
   end do

endif

! L2 norm and midpoints for L infinity norm

if (present(Linf) .or. present(L2)) then

   if (present(L2)) L2 = 0.0_my_real

! go through the elements of the grid

   do lev=1,grid%nlev
      elem = grid%head_level_elem(lev)
      do while (elem /= END_OF_LIST)
         if (grid%element(elem)%isleaf .and. grid%element(elem)%iown) then

! compute the quadrature points, evaluate the solution, and compute the
! integral and max over the quadrature points

            xc = grid%vertex(grid%element(elem)%vertex)%coord%x
            yc = grid%vertex(grid%element(elem)%vertex)%coord%y
            call quadrature_rule_tri(6,xc,yc,nqp,qw,xq,yq,i,stay_in=.true.)
            allocate(u(1,1,nqp),stat=astat)
            if (astat /= 0) then
               ierr = ALLOC_FAILED
               call fatal("allocation failed in norm_solution",procs=procs)
               return
            endif
            call evaluate_soln_local(grid,xq,yq,elem,(/comp/),(/eigen/),u)
            do i=1,nqp
               if (present(Linf)) Linf = max(Linf,abs(u(1,1,i)))
               if (present(L2)) L2 = L2 + qw(i)*u(1,1,i)**2
            end do
            deallocate(qw,xq,yq,u)

         endif
         elem = grid%element(elem)%next
      end do
   end do

   if (present(L2)) L2 = sqrt(abs(L2))

endif

end subroutine norm_solution

!          ---------
subroutine norm_true(grid,procs,still_sequential,comp,eigen,linf,l2,energy)
!          ---------

!----------------------------------------------------
! This routine computes the requested norms of the requested component of the
! true solution over the region owned by this processor.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments
 
type (grid_type), intent(in) :: grid
type (proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential
integer, intent(in) :: comp,eigen
real (my_real), optional, intent(out) :: linf, l2, energy
!----------------------------------------------------
 
!----------------------------------------------------
! Local variables:
 
real (my_real), pointer :: qw(:),xq(:),yq(:)
real (my_real), allocatable :: u(:,:,:),ux(:,:,:),uy(:,:,:)
real(my_real) :: cxx(grid%system_size,grid%system_size), &
                 cxy(grid%system_size,grid%system_size), &
                 cyy(grid%system_size,grid%system_size), &
                 cx(grid%system_size,grid%system_size),  &
                 cy(grid%system_size,grid%system_size),  &
                 c(grid%system_size,grid%system_size),   &
                 rs(grid%system_size)
integer :: lev,vert,elem,astat,i,j,k,nqp,ss
real (my_real) :: xc(3), yc(3)
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

if (my_proc(procs)==MASTER) return

ss = grid%system_size

! energy norm

if (present(energy)) then

   energy = 0.0_my_real
   do lev=1,grid%nlev
      elem = grid%head_level_elem(lev)
      do while (elem /= END_OF_LIST)
         if (.not.grid%element(elem)%isleaf .or. &
             .not.grid%element(elem)%iown) then
            elem = grid%element(elem)%next
            cycle
         endif
         xc = grid%vertex(grid%element(elem)%vertex)%coord%x
         yc = grid%vertex(grid%element(elem)%vertex)%coord%y
! quadrature order determined experimentally with u=x**10+y**10 on uniform grid
         call quadrature_rule_tri(6,xc,yc,nqp,qw,xq,yq,i,stay_in=.true.)
         allocate(u(ss,1,nqp),ux(ss,1,nqp),uy(ss,1,nqp),stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("allocation failed in norm_true",procs=procs)
            return
         endif
         do i=1,nqp
            call pdecoefs(xq(i),yq(i),cxx,cxy,cyy,cx,cy,c,rs)
            do j=1,ss
               u(j,1,i) = trues(xq(i),yq(i),j,eigen)
               ux(j,1,i) = truexs(xq(i),yq(i),j,eigen)
               uy(j,1,i) = trueys(xq(i),yq(i),j,eigen)
            end do
            do j=1,ss
               do k=1,ss
                  energy = energy + qw(i)*(cxx(j,k)*ux(j,1,i)*ux(k,1,i) + &
                                           cyy(j,k)*uy(j,1,i)*uy(k,1,i) + &
                                           cxy(j,k)*ux(j,1,i)*uy(k,1,i) + &
                                            cx(j,k)*ux(j,1,i)*u(k,1,i) + &
                                            cy(j,k)*uy(j,1,i)*u(k,1,i) + &
                                             c(j,k)*u(j,1,i)*u(k,1,i))
               end do
            end do
         end do
         deallocate(qw,xq,yq,u,ux,uy)
         elem = grid%element(elem)%next
      end do
   end do

   energy = sqrt(abs(energy))

endif

! Linf norm

if (present(Linf)) then

   Linf = 0.0_my_real

! go through the vertices of the grid

   do lev=1,grid%nlev
      vert = grid%head_level_vert(lev)
      do while (vert /= END_OF_LIST)
         if (grid%element(grid%vertex(vert)%assoc_elem)%iown) then
           Linf = max(Linf, &
                    trues(grid%vertex(vert)%coord%x,grid%vertex(vert)%coord%y, &
                          comp,eigen))
         endif
         vert = grid%vertex(vert)%next
      end do
   end do

endif

! L2 norm and midpoints for L infinity norm

if (present(Linf) .or. present(L2)) then

   if (present(L2)) L2 = 0.0_my_real

! go through the elements of the grid

   do lev=1,grid%nlev
      elem = grid%head_level_elem(lev)
      do while (elem /= END_OF_LIST)
         if (grid%element(elem)%isleaf .and. grid%element(elem)%iown) then

! compute the quadrature points, evaluate the solution, and compute the
! integral and max over the quadrature points

            xc = grid%vertex(grid%element(elem)%vertex)%coord%x
            yc = grid%vertex(grid%element(elem)%vertex)%coord%y
            call quadrature_rule_tri(6,xc,yc,nqp,qw,xq,yq,i,stay_in=.true.)
            allocate(u(1,1,nqp),stat=astat)
            if (astat /= 0) then
               ierr = ALLOC_FAILED
               call fatal("allocation failed in norm_true",procs=procs)
               return
            endif
            do i=1,nqp
               u(1,1,i) = trues(xq(i),yq(i),comp,eigen)
               if (present(Linf)) Linf = max(Linf,abs(u(1,1,i)))
               if (present(L2)) L2 = L2 + qw(i)*u(1,1,i)**2
            end do
            deallocate(qw,xq,yq,u)

         endif
         elem = grid%element(elem)%next
      end do
   end do

   if (present(L2)) L2 = sqrt(abs(L2))

endif

end subroutine norm_true

!          ------------
subroutine store_matrix(grid,procs,still_sequential,system_size,eq_type, &
                        inc_quad_order,sunit,runit,munit)
!          ------------

!----------------------------------------------------
! This routine writes the stiffness matrix, mass matrix and/or right hand
! side to files in Matrix Market format.  The rhs is an Nx1 matrix.
! The master collects the matrices from the slaves and writes to the
! given unit numbers, which must be opened when phaml_store_matrix is
! called by the user.  The presence of the unit arguments determine which
! matrices are written.  sunit is for the stiffness matrix, runit for the
! right hand side, and munit for the mass matrix.  The order of the rows
! is to merge the matrices of the slaves by level, and within each level
! order groups by the processor number, and within the group the order in
! which the processor has them.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
type(proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential
integer, intent(in) :: system_size, eq_type, inc_quad_order
integer, intent(in), optional :: sunit, runit, munit
!----------------------------------------------------
! Local variables:

integer, allocatable :: isend(:), nsend(:), nrecv(:)
real(my_real), allocatable :: rsend(:)
integer, pointer :: irecv(:)
real(my_real), pointer :: rrecv(:)
type(linsys_type) :: linsys
type(solver_options) :: solver_cntl
type(io_options) :: io_cntl
type(hash_key_eq) :: gid
integer, allocatable :: lid2global(:), nown(:,:), start_row(:,:), &
                        hold_begin(:)
integer :: nproc, my_processor, whichmat, iounit, p, proc, nlev, astat, &
           ni, nr, lev, neq, nnonzero, eq, col, geq, counter, i, lid
!----------------------------------------------------
! Begin executable code

! if compilation is for sequential, the lone processor writes the file

if (PARALLEL == SEQUENTIAL) then

! set up solver control and i/o control as needed to create linear system

   solver_cntl%system_size = system_size
   solver_cntl%eq_type = eq_type
   solver_cntl%lambda0 = 0.0_my_real
   solver_cntl%num_eval = 0
   solver_cntl%ignore_quad_err = .true.
   solver_cntl%inc_quad_order = inc_quad_order

   io_cntl%print_error_when = FREQUENTLY ! prevents condensed from overwriting stiffness

! create the linear system

   call create_linear_system(linsys,grid,procs,solver_cntl, &
                             io_cntl,still_sequential,notime=.true.)

! remove static condensation

   linsys%neq = linsys%neq_vert + linsys%neq_edge + &
                linsys%neq_face
   linsys%end_row => linsys%end_row_face
   linsys%matrix_val => linsys%stiffness
   linsys%rhs => linsys%rhs_nocond

! count the nonzeroes

   nnonzero = 0
   do eq=1,linsys%neq
      if (linsys%equation_type(eq) == DIRICHLET) then
         nnonzero = nnonzero + 1
      else
         do col=linsys%begin_row(eq),linsys%end_row(eq)
            if (linsys%column_index(col) == NO_ENTRY) cycle
            nnonzero = nnonzero + 1
         end do
      endif
   end do

! for each matrix to be stored ...

   do whichmat=1,3

! set the output unit

      select case(whichmat)
      case (1)
         if (.not. present(sunit)) cycle
         iounit = sunit
      case (2)
         if (.not. present(runit)) cycle
         iounit = runit
      case (3)
         if (.not. present(munit)) cycle
         iounit = munit
      end select

! write the header information

      write(iounit,"(A)") "%%MatrixMarket matrix coordinate real general"
      if (whichmat == 2) then
         write(iounit,"(3I11)") linsys%neq, 1, linsys%neq ! rhs
      else
         write(iounit,"(3I11)") linsys%neq, linsys%neq, nnonzero ! stiffness and mass
      endif

! for each equation ....

      do eq=1,linsys%neq

         if (whichmat == 2) then

! rhs writes the rhs or the solution if it is a Dirichlet eq

            if (linsys%equation_type(eq) == DIRICHLET) then
               write(iounit,"(2I11,SS,1P,E18.10E2)") eq,1,linsys%solution(eq)
            else
               write(iounit,"(2I11,SS,1P,E18.10E2)") eq,1,linsys%rhs(eq)
            endif

         else

! matrices go through the row, or just 1 on the diagonal for Dirichlet equations

            if (linsys%equation_type(eq) == DIRICHLET) then
               write(iounit,"(2I11,SS,1P,E18.10E2)") eq,eq,1.0_my_real
            else
               do col=linsys%begin_row(eq),linsys%end_row(eq)
                  if (linsys%column_index(col) == NO_ENTRY) cycle
                  if (whichmat == 1) then
                     write(iounit,"(2I11,SS,1P,E18.10E2)") eq, &
                                                   linsys%column_index(col), &
                                                   linsys%stiffness(col)
                  else
                     write(iounit,"(2I11,SS,1P,E18.10E2)") eq, &
                                                   linsys%column_index(col), &
                                                   linsys%mass(col)
                  endif
               end do
            endif

         endif
      end do
   end do

! destroy the linear system

   call destroy_linear_system(linsys)

else ! not sequential compilation (dropping one level of indentation)

! useful values

nproc = num_proc(procs)
my_processor = my_proc(procs)

! if still sequential, processor 1 has the entire grid.
! CAUTION: if create_linear_system needs communication among the slaves
!          (which it doesn't at the time of this writing), should return
!          after creating (and destroying) the linear system

if (still_sequential) then
   if (my_processor > 1) return
   nproc = 1
endif

! code for MASTER

if (my_processor == MASTER) then

! receive the number of levels

   call phaml_recv(procs,proc,irecv,ni,rrecv,nr,1321)
   nlev = irecv(1)
   deallocate(irecv,stat=astat)

! receive from each processor the number of owned equations on each level

   allocate(nown(nproc,nlev+2),start_row(nproc,nlev+2),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in store_matrix",procs=procs)
      return
   endif
   nown = 0
   do p=1,nproc
      call phaml_recv(procs,proc,irecv,ni,rrecv,nr,1322)
      nown(proc,1:ni) = irecv
      deallocate(irecv,stat=astat)
   end do

! determine the global matrix numbering of the beginning of each level on
! each processor

   start_row(1,1) = 1
   do proc=2,nproc
      start_row(proc,1) = start_row(proc-1,1) + nown(proc-1,1)
   end do
   do lev=2,nlev+2
      start_row(1,lev) = start_row(nproc,lev-1) + nown(nproc,lev-1)
      do proc=2,nproc
         start_row(proc,lev) = start_row(proc-1,lev) + nown(proc-1,lev)
      end do
   end do

   neq = start_row(nproc,nlev+2) + nown(nproc,nlev+2) - 1

! receive the number of nonzeroes

   nnonzero = 0
   do p=1,nproc
      call phaml_recv(procs,proc,irecv,ni,rrecv,nr,1323)
      nnonzero = nnonzero + irecv(1)
      deallocate(irecv,stat=astat)
   end do

! send each processor it's global matrix numbering at the start of each level

   do proc=1,nproc
      call phaml_send(procs,proc,start_row(proc,:),nlev+2,(/0.0_my_real/),0, &
                      1324)
   end do

! for each matrix to be stored ...

   do whichmat=1,3

! set the output unit

      select case(whichmat)
      case (1)
         if (.not. present(sunit)) cycle
         iounit = sunit
      case (2)
         if (.not. present(runit)) cycle
         iounit = runit
      case (3)
         if (.not. present(munit)) cycle
         iounit = munit
      end select

! write the header information

      write(iounit,"(A)") "%%MatrixMarket matrix coordinate real general"
      if (whichmat == 2) then
         write(iounit,"(3I11)") neq, 1, neq ! rhs
      else
         write(iounit,"(3I11)") neq, neq, nnonzero ! stiffness and mass
      endif

! for each level ...

      do lev=1,nlev+2

! request the equations of this level from the slaves

         do p=1,nproc
            call phaml_send(procs,p,(/whichmat,lev/),2,(/0.0_my_real/),0,1327)
         end do

! for each processor ...

         do p=1,nproc

! receive the equations

            call phaml_recv(procs,proc,irecv,ni,rrecv,nr,1328)

! write them to the output file

            do i=1,nr
               if (whichmat == 2) then ! rhs
                  write(iounit,"(2I11,SS,1P,E18.10E2)") irecv(i),1,rrecv(i)
               else ! stiffness and mass
                  write(iounit,"(2I11,SS,1P,E18.10E2)") irecv(2*(i-1)+1), &
                                                      irecv(2*(i-1)+2),rrecv(i)
               endif
            end do

            if (associated(irecv)) deallocate(irecv,stat=astat)
            if (associated(rrecv)) deallocate(rrecv,stat=astat)

         end do ! next processor
      end do ! next level
   end do ! next matrix

   deallocate(nown,start_row,stat=astat)

! tell the slaves to quit
   
   do p=1,nproc
      call phaml_send(procs,p,(/4/),1,(/0.0_my_real/),0,1327)
   end do

! code for SLAVES

else ! slave

! set up solver control and i/o control as needed to create linear system

   solver_cntl%system_size = system_size
   solver_cntl%eq_type = eq_type
   solver_cntl%lambda0 = 0.0_my_real
   solver_cntl%num_eval = 0
   solver_cntl%ignore_quad_err = .true.
   solver_cntl%inc_quad_order = inc_quad_order

   io_cntl%print_error_when = FREQUENTLY ! prevents condensed from overwriting stiffness

! create the linear system

   call create_linear_system(linsys,grid,procs,solver_cntl, &
                             io_cntl,still_sequential,notime=.true.)

! remove static condensation

   linsys%neq = linsys%neq_vert + linsys%neq_edge + &
                linsys%neq_face
   linsys%end_row => linsys%end_row_face
   linsys%matrix_val => linsys%stiffness
   linsys%rhs => linsys%rhs_nocond

! send the number of levels to the master

   if (still_sequential) then
      nlev = linsys%nlev
   else
      nlev = phaml_global_max(procs,linsys%nlev,1320)
   endif
   if (my_processor == 1) call phaml_send(procs,MASTER,(/nlev/),1, &
                                          (/0.0_my_real/),0,1321)

! change the beginning of the levels so that every processor uses the
! same value for nlev=1 and nlev+2

   if (linsys%nlev /= nlev) then
      allocate(hold_begin(linsys%nlev+3))
      hold_begin = linsys%begin_level
      deallocate(linsys%begin_level)
      allocate(linsys%begin_level(nlev+3))
      linsys%begin_level(1:linsys%nlev+3) = hold_begin
      deallocate(hold_begin)
      linsys%begin_level(nlev+3) = linsys%begin_level(linsys%nlev+3)
      linsys%begin_level(nlev+2) = linsys%begin_level(linsys%nlev+2)
      linsys%begin_level(nlev+1) = linsys%begin_level(linsys%nlev+1)
      do lev=linsys%nlev+1,nlev
         linsys%begin_level(lev) = linsys%begin_level(nlev+1)
      end do
      linsys%nlev = nlev
   endif

! count the number of owned equations on each level and number of nonzeros

   allocate(nown(linsys%nlev+2,1),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in store_matrix",procs=procs)
      return
   endif

   nown = 0
   nnonzero = 0
   lev = 1
   do eq=1,linsys%neq
      do while (eq >= linsys%begin_level(lev+1))
         lev = lev+1
      end do
      if (linsys%iown(eq)) then
         nown(lev,1) = nown(lev,1) + 1
         if (linsys%equation_type(eq) == DIRICHLET) then
            nnonzero = nnonzero + 1
         else
            do col=linsys%begin_row(eq),linsys%end_row(eq)
               if (linsys%column_index(col) /= NO_ENTRY) nnonzero = nnonzero + 1
            end do
         endif
      endif
   end do

! send the number of owned equations on each level to the master

   call phaml_send(procs,MASTER,nown(:,1),linsys%nlev+2,(/0.0_my_real/),0,1322)

   deallocate(nown,stat=astat)

! send the number of nonzeroes in the rows owned by this processor to the master

   call phaml_send(procs,MASTER,(/nnonzero/),1,(/0.0_my_real/),0,1323)

! receive the global matrix numbering of the beginning of each level

   call phaml_recv(procs,proc,irecv,ni,rrecv,nr,1324)

! set the mapping from my local id to the global matrix numbering

   allocate(lid2global(linsys%neq),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in store_matrix",procs=procs)
      return
   endif

   lid2global = 0
   lev = 0
   do eq=1,linsys%neq
      do while (eq >= linsys%begin_level(lev+1))
         lev = lev+1
         geq = irecv(lev)
      end do
      if (linsys%iown(eq)) then
         lid2global(eq) = geq
         geq = geq + 1
      endif
   end do

   deallocate(irecv,stat=astat)

! request global matrix numbering for unowned columns from other processors

   if (.not. still_sequential .and. nproc > 1 .and. parallel /= SEQUENTIAL) then

      counter = 0
      do eq=1,linsys%neq
         if (linsys%iown(eq)) then
            do col=linsys%begin_row(eq),linsys%end_row(eq)
               if (linsys%column_index(col) /= NO_ENTRY) then
                  if (.not. linsys%iown(linsys%column_index(col))) then
                     counter = counter + KEY_SIZE+1
                  endif
               endif
            end do
         end if
      end do

      allocate(isend(counter),nrecv(nproc),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in store_matrix",procs=procs)
         return
      endif

      counter = 1
      do eq=1,linsys%neq
         if (linsys%iown(eq)) then
            do col=linsys%begin_row(eq),linsys%end_row(eq)
               if (linsys%column_index(col) /= NO_ENTRY) then
                  if (.not. linsys%iown(linsys%column_index(col))) then
                     call hash_pack_key(linsys%gid(linsys%column_index(col)), &
                                        isend,counter)
                     counter = counter + KEY_SIZE+1
                  endif
               endif
            end do
         end if
      end do

      call phaml_alltoall(procs,isend,counter-1,irecv,nrecv,1325)

      deallocate(isend,stat=astat)

! construct response for those that I own

      counter = 1
      ni = 0
      do p=1,nproc
         if (p == my_processor) then
            counter = counter + nrecv(p)
            cycle
         endif
         do i=1,nrecv(p)/(KEY_SIZE+1)
            gid = hash_unpack_key(irecv,counter,extended=.true.)
            lid = hash_decode_key(gid,linsys%eq_hash)
            if (lid /= HASH_NOT_FOUND) then
               if (linsys%iown(lid)) then
                  ni = ni + KEY_SIZE+2
               endif
            endif
            counter = counter + KEY_SIZE+1
         end do
      end do

      allocate(isend(ni),nsend(nproc),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in store_matrix",procs=procs)
         return
      endif

      counter = 1
      ni = 1
      do p=1,nproc
         if (p == my_processor) then
            counter = counter + nrecv(p)
            nsend(p) = 0
            cycle
         endif
         do i=1,nrecv(p)/(KEY_SIZE+1)
            gid = hash_unpack_key(irecv,counter,extended=.true.)
            lid = hash_decode_key(gid,linsys%eq_hash)
            if (lid /= HASH_NOT_FOUND) then
               if (linsys%iown(lid)) then
                  call hash_pack_key(gid,isend,ni)
                  ni = ni + KEY_SIZE+1
                  isend(ni) = lid2global(lid)
                  ni = ni + 1
               endif
            endif
            counter = counter + KEY_SIZE+1
         end do
         nsend(p) = ni - sum(nsend(1:p-1)) - 1
      end do

      if (associated(irecv)) deallocate(irecv,stat=astat)

      call phaml_alltoall(procs,isend,nsend,irecv,nrecv,1326)

      deallocate(isend,nsend,stat=astat)

! set the mapping from my local id to the global matrix numbering received

      counter = 1
      do p=1,nproc
         do i=1,nrecv(p)/(KEY_SIZE+2)
            gid = hash_unpack_key(irecv,counter,extended=.true.)
            lid = hash_decode_key(gid,linsys%eq_hash)
            counter = counter + KEY_SIZE+1
            if (lid /= HASH_NOT_FOUND) then
               lid2global(lid) = irecv(counter)
            endif
            counter = counter + 1
         end do
      end do

      if (associated(irecv)) deallocate(irecv,stat=astat)
      deallocate(nrecv,stat=astat)

   endif ! not still sequential

! loop to process requests from the master

   do

! receive request.  First integer is which matrix to send, or quit.
! Second integer is level to send.

      call phaml_recv(procs,proc,irecv,ni,rrecv,nr,1327)

      whichmat = irecv(1)
      select case (whichmat)

! send stiffness or mass matrix

      case(1,3)

         lev = irecv(2)
         deallocate(irecv,stat=astat)

         counter = 0
         do eq=linsys%begin_level(lev),linsys%begin_level(lev+1)-1
            if (.not. linsys%iown(eq)) cycle
            geq = lid2global(eq)
            if (linsys%equation_type(eq) == DIRICHLET) then
               counter = counter + 1
            else
               do col=linsys%begin_row(eq),linsys%end_row(eq)
                  if (linsys%column_index(col) /= NO_ENTRY) then
                     counter = counter + 1
                  endif
               end do
            endif
         end do

         allocate(isend(2*counter),rsend(counter),stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("allocation failed in store_matrix",procs=procs)
            return
         endif

         counter = 0
         do eq=linsys%begin_level(lev),linsys%begin_level(lev+1)-1
            if (.not. linsys%iown(eq)) cycle
            geq = lid2global(eq)
            if (linsys%equation_type(eq) == DIRICHLET) then
               isend(2*counter+1) = geq
               isend(2*counter+2) = geq
               rsend(counter+1) = 1.0_my_real
               counter = counter + 1
            else
               do col=linsys%begin_row(eq),linsys%end_row(eq)
                  if (linsys%column_index(col) /= NO_ENTRY) then
                     isend(2*counter+1) = geq
                     isend(2*counter+2) = lid2global(linsys%column_index(col))
                     if (whichmat == 1) then
                        rsend(counter+1) = linsys%stiffness(col)
                     else
                        rsend(counter+1) = linsys%mass(col)
                     endif
                     counter = counter + 1
                  endif
               end do
            endif
         end do

         call phaml_send(procs,MASTER,isend,2*counter,rsend,counter,1328)

         deallocate(isend,rsend,stat=astat)

! send right hand side

      case(2)

         lev = irecv(2)
         deallocate(irecv,stat=astat)

         counter = 0
         do eq=linsys%begin_level(lev),linsys%begin_level(lev+1)-1
            if (linsys%iown(eq)) counter = counter + 1
         end do

         allocate(isend(counter),rsend(counter),stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("allocation failed in store_matrix",procs=procs)
            return
         endif

         counter = 0
         do eq=linsys%begin_level(lev),linsys%begin_level(lev+1)-1
            if (.not. linsys%iown(eq)) cycle
            counter = counter + 1
            isend(counter) = lid2global(eq)
            if (linsys%equation_type(eq) == DIRICHLET) then
               rsend(counter) = linsys%solution(eq)
            else
               rsend(counter) = linsys%rhs(eq)
            endif
         end do

         call phaml_send(procs,MASTER,isend,counter,rsend,counter,1328)

         deallocate(isend,rsend,stat=astat)

! quit

      case(4)

         deallocate(irecv,stat=astat)
         exit

      end select
   end do

! destroy the linear system

   call destroy_linear_system(linsys)

   deallocate(lid2global,stat=astat)

endif ! master or slave

endif ! sequential compilation

end subroutine store_matrix

end module linsys_io
