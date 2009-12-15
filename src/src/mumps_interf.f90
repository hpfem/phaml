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

module mumps_interf

!----------------------------------------------------
! This module contains the interface to the mumps package.
!
! communication tags in this module are of the form 15xx
!----------------------------------------------------

!----------------------------------------------------
! Other modules used are:

use global
use message_passing
use mumps_struc_mod
use hash_mod
use hash_eq_mod
use linsystype_mod
use linsys_util

!----------------------------------------------------
implicit none
private
public create_mumps_linear_system, destroy_mumps_linear_system, &
       change_mumps_rhs, mumps_solve

!----------------------------------------------------
! The following types are defined:

!----------------------------------------------------
! The following parameters are defined:

!----------------------------------------------------
! The following variables are defined:

integer, save :: fillin_size = 25
!----------------------------------------------------

contains

!          --------------------------
subroutine create_mumps_linear_system(mumps_matrix,phaml_matrix,procs, &
                                      solver_cntl,still_sequential)
!          --------------------------

!----------------------------------------------------
! This routine initializes the mumps data structure, creates and factors
! the matrix, and creates the rhs.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(mumps_matrix_type), intent(inout) :: mumps_matrix
type(linsys_type), intent(in) :: phaml_matrix
type(proc_info), intent(in) :: procs
type(solver_options), intent(in) :: solver_cntl
logical, intent(in) :: still_sequential

!----------------------------------------------------
! Local variables:

integer :: allocstat
!----------------------------------------------------
! Begin executable code

! master does not participate

if (my_proc(procs) == MASTER) return

! allocate space to keep track of equation numbering

allocate(mumps_matrix%firsteq(num_proc(procs)+1), &
         mumps_matrix%global_eq(phaml_matrix%neq),stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in create_mumps_linear_system",procs=procs)
   return
endif

! start global_eq at -1 so any reference to an unset one will generate
! an error

mumps_matrix%global_eq = -1

! initialize mumps_var

call initialize_mumps(mumps_matrix%mumps_var,procs,solver_cntl)

! copy the matrix to MUMPS format

call make_mumps_matrix(mumps_matrix,phaml_matrix,procs,still_sequential)

! copy the right hand side to MUMPS format

call make_mumps_rhs(mumps_matrix,phaml_matrix,phaml_matrix%rhs,procs, &
                    still_sequential)

! factor the matrix

! if no equations, nothing to factor

if (mumps_matrix%mumps_var%n > 0) then

! with 1 or 2 equations, PORD fails.  Force MUMPS to use Approximate
! Minimum Degree in that case.

   if (mumps_matrix%mumps_var%n < 3) then
      mumps_matrix%mumps_var%icntl(7) = 0
   endif

! analysis and factor in one call

   mumps_matrix%mumps_var%job = 4
   do
      call dmumps(mumps_matrix%mumps_var)
      if (mumps_matrix%mumps_var%info(1) /= -9) exit
! if info is -9, increase space for fill-in and try again
      fillin_size = 1.5*mumps_matrix%mumps_var%icntl(14)
      mumps_matrix%mumps_var%icntl(14) = fillin_size
      call warning("increasing space for MUMPS fill-in",intlist=(/fillin_size/))
   end do
   if (mumps_matrix%mumps_var%info(1) < 0) then
      call fatal("MUMPS returned error.  Code is", &
                 intlist=(/mumps_matrix%mumps_var%info(1), &
                           mumps_matrix%mumps_var%info(2)/),procs=procs)
      stop
   elseif (mumps_matrix%mumps_var%info(1) > 0) then
      call warning("MUMPS returned warning. Code is", &
                    intlist=(/mumps_matrix%mumps_var%info(1), &
                              mumps_matrix%mumps_var%info(2)/))
   endif

endif

end subroutine create_mumps_linear_system

!          ---------------------------
subroutine destroy_mumps_linear_system(mumps_matrix,procs)
!          ---------------------------

!----------------------------------------------------
! This routine deallocates memory in mumps_var and calls mumps destroy
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(mumps_matrix_type), intent(inout) :: mumps_matrix
type(proc_info), intent(in) :: procs
!----------------------------------------------------
! Local variables:

integer :: allocstat
!----------------------------------------------------
! Begin executable code

! master does not participate

if (my_proc(procs) == MASTER) return

if (associated(mumps_matrix%firsteq)) &
   deallocate(mumps_matrix%firsteq,stat=allocstat)
if (associated(mumps_matrix%global_eq)) &
   deallocate(mumps_matrix%global_eq,stat=allocstat)
if (associated(mumps_matrix%mumps_var%irn_loc)) &
   deallocate(mumps_matrix%mumps_var%irn_loc,stat=allocstat)
if (associated(mumps_matrix%mumps_var%jcn_loc)) &
   deallocate(mumps_matrix%mumps_var%jcn_loc,stat=allocstat)
if (associated(mumps_matrix%mumps_var%a_loc)) &
   deallocate(mumps_matrix%mumps_var%a_loc,stat=allocstat)
if (associated(mumps_matrix%mumps_var%rhs)) &
   deallocate(mumps_matrix%mumps_var%rhs,stat=allocstat)

mumps_matrix%mumps_var%job = -2
call dmumps(mumps_matrix%mumps_var)
if (mumps_matrix%mumps_var%info(1) < 0) then
   call fatal("MUMPS returned error.  Code is", &
              intlist=(/mumps_matrix%mumps_var%info(1), &
                        mumps_matrix%mumps_var%info(2)/),procs=procs)
   stop
elseif (mumps_matrix%mumps_var%info(1) > 0) then
   call warning("MUMPS returned warning. Code is", &
                 intlist=(/mumps_matrix%mumps_var%info(1), &
                           mumps_matrix%mumps_var%info(2)/))
endif

end subroutine destroy_mumps_linear_system

!          ----------------
subroutine change_mumps_rhs(mumps_matrix,phaml_matrix,rs,procs,still_sequential)
!          ----------------

!----------------------------------------------------
! This routine changes the rhs in the mumps linear system
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(mumps_matrix_type), intent(inout) :: mumps_matrix
type(linsys_type), intent(in) :: phaml_matrix
real(my_real), intent(in) :: rs(:)
type(proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential

!----------------------------------------------------
! Begin executable code

! master does not participate

if (my_proc(procs) == MASTER) return

! same thing as making the rhs in the first place

call make_mumps_rhs(mumps_matrix,phaml_matrix,rs,procs,still_sequential)

end subroutine change_mumps_rhs

!          -----------
subroutine mumps_solve(mumps_matrix,phaml_matrix,procs,still_sequential, &
                       noshadow)
!          -----------

!----------------------------------------------------
! This routine solves the linear system with an already factored MUMPS matrix.
! If noshadow is present and true, the solution at the shadow vertices is not
! communicated.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(mumps_matrix_type), intent(inout) :: mumps_matrix
type(linsys_type), intent(inout) :: phaml_matrix
type(proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential
logical, intent(in), optional :: noshadow

!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

! master does not participate

if (my_proc(procs) == MASTER) return

! if no equations, nothing to solve

if (mumps_matrix%mumps_var%n > 0) then

! call mumps to solve

   mumps_matrix%mumps_var%job = 3
   call dmumps(mumps_matrix%mumps_var)
   if (mumps_matrix%mumps_var%info(1) < 0) then
      call fatal("MUMPS returned error.  Code is", &
                 intlist=(/mumps_matrix%mumps_var%info(1), &
                           mumps_matrix%mumps_var%info(2)/),procs=procs)
      stop
   elseif (mumps_matrix%mumps_var%info(1) > 0) then
      call warning("MUMPS returned warning. Code is", &
                    intlist=(/mumps_matrix%mumps_var%info(1), &
                              mumps_matrix%mumps_var%info(2)/))
   endif

! distribute the solution over the processors

   call mumps_distrib_soln(mumps_matrix,phaml_matrix,procs, &
                           still_sequential,noshadow)

endif

end subroutine mumps_solve

!          ----------------
subroutine initialize_mumps(mumps_var,procs,solver_cntl)
!          ----------------

!----------------------------------------------------
! This routine initializes the mumps data structure.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(dmumps_struc), intent(inout) :: mumps_var
type(proc_info), intent(in) :: procs
type(solver_options), intent(in) :: solver_cntl
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

! tell MUMPS to use the slaves communicator

mumps_var%comm = slaves_communicator(procs)

! select SPD or general symmetric solver

select case (solver_cntl%solver)
case (MUMPS_NONSYM_SOLVER)
   mumps_var%sym = 0
case (MUMPS_SPD_SOLVER)
   mumps_var%sym = 1
case (MUMPS_GEN_SOLVER)
   mumps_var%sym = 2
case default
   ierr = PHAML_INTERNAL_ERROR
   call fatal("mumps_solve called, but solver is not MUMPS",procs=procs)
   return
end select

! slave 1 particpates in the solution process

mumps_var%par = 1

! have MUMPS initialize mumps_var

mumps_var%job = -1
call dmumps(mumps_var)
if (mumps_var%info(1) < 0) then
   call fatal("MUMPS returned error.  Code is",intlist=(/mumps_var%info(1), &
               mumps_var%info(2)/),procs=procs)
   stop
elseif (mumps_var%info(1) > 0) then
   call warning("MUMPS returned warning. Code is",intlist=(/mumps_var%info(1), &
                mumps_var%info(2)/))
endif

! the rest of my nondefault selections

mumps_var%icntl(1) = errunit ! unit for MUMPS error messages
mumps_var%icntl(3) = 0 ! suppress printing diagnostics
mumps_var%icntl(4) = 2 ! print only errors and warnings
mumps_var%icntl(5) = 0 ! the matrix is in assembled form
mumps_var%icntl(10) = 5 ! perform up to 5 steps of iterative refinement
mumps_var%icntl(14) = fillin_size ! extra space for fill-in
mumps_var%icntl(18) = 3 ! the matrix is already distributed
mumps_var%cntl(2) = 1.0e-13_my_real ! requirement to stop iterative refinement
mumps_var%cntl(3) = 0  ! default value can erroneously report singular matrix

end subroutine initialize_mumps

!          -----------------
subroutine make_mumps_matrix(mumps_matrix,phaml_matrix,procs, &
                             still_sequential)
!          -----------------

!----------------------------------------------------
! This routine copies the matrix to MUMPS format in mumps_var.
! The components of mumps_var that are set are:
!    n, nz_loc, irn_loc, jcn_loc, a_loc
! Also returns the global equation number assigned to each equation.  Dirichlet
! points get equation number DIRICHLET.
! Also returns the first equation, in global numbering, on each processor.
! This should be dimensioned at nproc+1, as an extra entry gets N+1 so
! the section for a processor is given by firsteq(p):firsteq(p+1)-1.
! When still_sequential, firsteq is set to 1 except in p+1 where it is N+1
!  
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(mumps_matrix_type), intent(inout), target :: mumps_matrix
type(linsys_type), intent(in) :: phaml_matrix
type(proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential
!----------------------------------------------------
! Local variables:

integer :: i, nproc, my_processor, my_neq, my_first_eq, counter, count_other, &
           not_own, allocstat, p, lid, j, col, isub, scounter, &
           rcounter, oscounter, KEY_SIZE_EQ
integer, allocatable :: neq_all(:)
! newcomm
integer, allocatable :: isend(:),nsendv(:),nrecv(:)
integer, pointer :: irecv(:)
type(dmumps_struc), pointer :: mumps_var
!----------------------------------------------------
! Begin executable code

nproc = num_proc(procs)
my_processor = my_proc(procs)
mumps_var => mumps_matrix%mumps_var
KEY_SIZE_EQ = KEY_SIZE+1

! if still sequential, processors other than 1 set up an empty matrix

if (still_sequential) then
   if (my_processor /= 1) then
      mumps_var%nz_loc = 0
      allocate(mumps_var%irn_loc(1),mumps_var%jcn_loc(1),mumps_var%a_loc(1))
      counter = 0
      do i=1,phaml_matrix%neq
         if (phaml_matrix%equation_type(i) == DIRICHLET) then
            mumps_matrix%global_eq(i) = DIRICHLET
         else
            counter = counter+1
            mumps_matrix%global_eq(i) = counter
         endif
      end do
      mumps_var%n = counter
      mumps_matrix%firsteq = 1
      mumps_matrix%firsteq(my_processor+1) = counter+1
      return
   endif
endif

! allocate space for received messages

allocate(nsendv(nproc), nrecv(nproc), stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in make_mumps_matrix",procs=procs)
   return
endif

! count the number of equations owned by this processor; do not include
! Dirichlet points, which are pre-eliminated.  Also count the number of
! equations owned by other processors.

my_neq = 0
not_own = 0
do i=1,phaml_matrix%neq
   if (phaml_matrix%iown(i)) then
      if (phaml_matrix%equation_type(i) /= DIRICHLET) my_neq = my_neq + 1
   else
      not_own = not_own + 1
   endif
end do

! determine how many equations are owned by each processor

allocate(neq_all(nproc))
neq_all = 0
neq_all(my_processor) = my_neq
if (.not. still_sequential) then
   neq_all = phaml_global_sum(procs,neq_all,1501)
endif

! set the total number of equations

mumps_var%n = sum(neq_all)

! set the global number of the first equation on each processor to be the
! sum of the number of equations on lower processor, plus 1.  Include
! an extra one so firsteq(nproc+1)-firsteq(nproc) works.
! Special case for still sequential, when only processor 1 gets here.

if (still_sequential) then
   mumps_matrix%firsteq = 1
   mumps_matrix%firsteq(my_processor+1) = mumps_var%n+1
else
   mumps_matrix%firsteq(1) = 1
   do i=2,nproc+1
      mumps_matrix%firsteq(i) = mumps_matrix%firsteq(i-1) + neq_all(i-1)
   end do
endif

deallocate(neq_all,stat=allocstat)
my_first_eq = mumps_matrix%firsteq(my_processor)

! determine the global number for each equation on this processor.
! First set the owned equations in the order they are found, using
! DIRICHLET (which is negative so it can't be a real equation number) for
! the Dirichlet points, and make a list of the gids of unowned equations.

allocate(isend(not_own*KEY_SIZE_EQ),stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in make_mumps_matrix",procs=procs)
   return
endif
counter = 0
count_other = 1
do i=1,phaml_matrix%neq
   if (phaml_matrix%iown(i)) then
      if (phaml_matrix%equation_type(i) == DIRICHLET) then
         mumps_matrix%global_eq(i) = DIRICHLET
      else
         mumps_matrix%global_eq(i) = my_first_eq + counter
         counter = counter + 1
      endif
   else
      call hash_pack_key(phaml_matrix%gid(i),isend,count_other)
      count_other = count_other + KEY_SIZE_EQ
   endif
end do

! send the list of unowned gids to the other processors

if (.not. still_sequential) then
   call phaml_alltoall(procs,isend,count_other-1,irecv,nrecv,1502)
endif

deallocate(isend, stat=allocstat) 

! if still sequential, there are no unowned equations

if (.not. still_sequential) then

! from the other processors' requests, send the global_eq of the ones I own

   counter = 0
   isub = 1
   do p=1,nproc
      if (p == my_processor) then
         isub = isub + nrecv(p)
         cycle
      endif
      do i=1,nrecv(p)/KEY_SIZE_EQ
         lid = hash_decode_key(hash_unpack_key(irecv,isub,.true.), &
                               phaml_matrix%eq_hash)
         if (lid /= HASH_NOT_FOUND) then
            if (phaml_matrix%iown(lid)) then
               counter = counter + 1
            endif
         endif
         isub = isub + KEY_SIZE_EQ
      end do
   end do

   allocate(isend(counter*(KEY_SIZE_EQ+1)),stat=allocstat)
   if (allocstat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in make_mumps_matrix",procs=procs)
      return
   endif

   scounter = 1
   rcounter = 1
   oscounter = scounter
   do p=1,nproc
      if (p == my_processor) then
         rcounter = rcounter + nrecv(p)
         nsendv(p) = 0
         cycle
      endif
      do i=1,nrecv(p)/KEY_SIZE_EQ
         lid = hash_decode_key(hash_unpack_key(irecv,rcounter,.true.), &
                               phaml_matrix%eq_hash)
         if (lid /= HASH_NOT_FOUND) then
            if (phaml_matrix%iown(lid)) then
               isend(scounter:scounter+KEY_SIZE_EQ-1) = &
                        irecv(rcounter:rcounter+KEY_SIZE_EQ-1)
               isend(scounter+KEY_SIZE_EQ) = mumps_matrix%global_eq(lid)
               scounter = scounter + KEY_SIZE_EQ+1
            endif
         endif
         rcounter = rcounter + KEY_SIZE_EQ
      end do
      nsendv(p) = scounter - oscounter
      oscounter = scounter
   end do

   if (associated(irecv)) deallocate(irecv,stat=allocstat)

! send the replies

   call phaml_alltoall(procs,isend,nsendv,irecv,nrecv,1503)
   deallocate(isend,stat=allocstat)

! copy the replies into global_eq

   counter=1
   do p=1,nproc
      if (p == my_processor) cycle
      do i=1,nrecv(p)/(KEY_SIZE_EQ+1)
         lid = hash_decode_key(hash_unpack_key(irecv,counter,.true.), &
                               phaml_matrix%eq_hash)
         if (lid == HASH_NOT_FOUND) then
            call warning("received reply for an equation I don't have in make_mumps_matrix")
         else
            mumps_matrix%global_eq(lid) = irecv(counter+KEY_SIZE_EQ)
         endif
         counter = counter + KEY_SIZE_EQ+1
      end do
   end do

   if (associated(irecv)) deallocate(irecv,stat=allocstat)
   deallocate(nsendv,nrecv,stat=allocstat)

endif ! .not. still_sequential

! count the number of nonzeros in the rows owned by this processor, omitting
! Dirichlet points and only counting the upper triangle if its symmetric

counter = 0
do i=1,phaml_matrix%neq
   if (.not. phaml_matrix%iown(i)) cycle
   if (phaml_matrix%equation_type(i) == DIRICHLET) cycle
   do j=phaml_matrix%begin_row(i),phaml_matrix%end_row(i)
      if (phaml_matrix%column_index(j) == NO_ENTRY) cycle
      col = mumps_matrix%global_eq(phaml_matrix%column_index(j))
      if (col == DIRICHLET) cycle
      if (col < mumps_matrix%global_eq(i) .and. mumps_var%sym /= 0) cycle
      counter = counter + 1
   end do
end do

mumps_var%nz_loc = counter

! copy the rows owned by this processor to MUMPS distributed matrix structure,
! only including the upper triangle

allocate(mumps_var%irn_loc(counter),mumps_var%jcn_loc(counter), &
         mumps_var%a_loc(counter),stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in make_mumps_matrix",procs=procs)
   return
endif

counter = 0
do i=1,phaml_matrix%neq
   if (.not. phaml_matrix%iown(i)) cycle
   if (phaml_matrix%equation_type(i) == DIRICHLET) cycle
   do j=phaml_matrix%begin_row(i),phaml_matrix%end_row(i)
      if (phaml_matrix%column_index(j) == NO_ENTRY) cycle
      col = mumps_matrix%global_eq(phaml_matrix%column_index(j))
      if (col == DIRICHLET) cycle
      if (col < mumps_matrix%global_eq(i) .and. mumps_var%sym /= 0) cycle
      counter = counter + 1
      mumps_var%irn_loc(counter) = mumps_matrix%global_eq(i)
      mumps_var%jcn_loc(counter) = col
      mumps_var%a_loc(counter) = phaml_matrix%matrix_val(j)
   end do
end do

end subroutine make_mumps_matrix

!          --------------
subroutine make_mumps_rhs(mumps_matrix,phaml_matrix,rs,procs,still_sequential)
!          --------------

!----------------------------------------------------
! This routine copies rs into the mumps data structure.  The entire
! right hand side goes on processor 1.  make_mumps_matrix should be
! called before this.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(mumps_matrix_type), intent(inout) :: mumps_matrix
type(linsys_type), intent(in) :: phaml_matrix
real(my_real), intent(in) :: rs(:)
type(proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential
!----------------------------------------------------
! Local variables:

integer :: send_int(phaml_matrix%neq)
real(my_real) :: send_real(phaml_matrix%neq)
integer, pointer :: recv_int(:)
real(my_real), pointer :: recv_real(:)
integer :: i, j, counter, p, ni, nr, allocstat
!----------------------------------------------------
! Begin executable code

if (my_proc(procs) /= 1) then ! processors other than 1

   nullify(mumps_matrix%mumps_var%rhs)

! if still sequential, processor 1 has the whole rhs

   if (still_sequential) return

! send the rhs entries this processor owns, and the corresponding global eq,
! to processor 1, after removing Dirichlet contributions

   counter = 0
   do i=1,phaml_matrix%neq
      if (.not. phaml_matrix%iown(i)) cycle
      if (phaml_matrix%equation_type(i) == DIRICHLET) cycle
      counter = counter + 1
      send_real(counter) = rs(i)
      send_int(counter) = mumps_matrix%global_eq(i)
      do j=phaml_matrix%begin_row(i),phaml_matrix%end_row(i)
         if (phaml_matrix%column_index(j) == NO_ENTRY) cycle
         if (phaml_matrix%equation_type(phaml_matrix%column_index(j)) == DIRICHLET) then
            send_real(counter) = send_real(counter) - &
               phaml_matrix%matrix_val(j)*phaml_matrix%solution(phaml_matrix%column_index(j))
         endif
      end do
   end do

   call phaml_send(procs,1,send_int,counter,send_real,counter,1504)

else ! processor 1

   allocate(mumps_matrix%mumps_var%rhs(mumps_matrix%mumps_var%n),stat=allocstat)
   if (allocstat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in make_mumps_rhs",procs=procs)
      return
   endif

! copy the rs that this processor owns, and remove Dirichlet contributions

   do i=1,phaml_matrix%neq
      if (.not. phaml_matrix%iown(i)) cycle
      if (phaml_matrix%equation_type(i) == DIRICHLET) cycle
      mumps_matrix%mumps_var%rhs(mumps_matrix%global_eq(i)) = rs(i)
      do j=phaml_matrix%begin_row(i),phaml_matrix%end_row(i)
         if (phaml_matrix%column_index(j) == NO_ENTRY) cycle
         if (phaml_matrix%equation_type(phaml_matrix%column_index(j)) == DIRICHLET) then
            mumps_matrix%mumps_var%rhs(mumps_matrix%global_eq(i)) =  &
               mumps_matrix%mumps_var%rhs(mumps_matrix%global_eq(i)) - &
               phaml_matrix%matrix_val(j)*phaml_matrix%solution(phaml_matrix%column_index(j))
         endif
      end do
   end do

! receive other rhs from other processors

   if (.not. still_sequential) then
      do i=2,num_proc(procs)
         call phaml_recv(procs,p,recv_int,ni,recv_real,nr,1504)
         do j=1,ni
            mumps_matrix%mumps_var%rhs(recv_int(j)) = recv_real(j)
         end do
         if (associated(recv_int)) deallocate(recv_int,recv_real,stat=allocstat)
      end do
   endif

endif ! processor 1

end subroutine make_mumps_rhs

!          ------------------
subroutine mumps_distrib_soln(mumps_matrix,phaml_matrix,procs, &
                              still_sequential,noshadow)
!          ------------------

!----------------------------------------------------
! This routine distributes the solution, contained in mumps_var%rhs on
! processor 1, to solution on each processor
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(mumps_matrix_type), intent(in) :: mumps_matrix
type(linsys_type), intent(inout) :: phaml_matrix
type(proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential
logical, intent(in), optional :: noshadow
!----------------------------------------------------
! Local variables:

integer :: p, i, allocstat, ni, nr, first, lasteq
integer, pointer :: recv_int(:)
real(my_real), pointer :: recv_real(:)
logical :: shadow
!----------------------------------------------------
! Begin executable code

if (present(noshadow)) then
   shadow = .not. noshadow
else
   shadow = .true.
endif

if (my_proc(procs) == 1) then

! processor 1 sends the section of solution that belongs to each processor
! to the owning processor

   do p=2,num_proc(procs)
      if (still_sequential) then
         first = 1
         lasteq = mumps_matrix%mumps_var%n
      else
         first = mumps_matrix%firsteq(p)
         lasteq = mumps_matrix%firsteq(p+1)-1
      endif
      call phaml_send(procs,p,(/0/),0,mumps_matrix%mumps_var%rhs(first:lasteq),&
                      lasteq-first+1,1505)
   end do

! and copies it's own section to recv_real

   allocate(recv_real(mumps_matrix%firsteq(2)-1),stat=allocstat)
   if (allocstat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in mumps_distrib_soln",procs=procs)
      return
   endif
   recv_real = mumps_matrix%mumps_var%rhs(1:mumps_matrix%firsteq(2)-1)

else

! other processors receive the solution from processor 1

   call phaml_recv(procs,p,recv_int,ni,recv_real,nr,1505)

endif

! copy into solution.  global_eq(i) gives the global equation number for
! this processor's local equation i (or DIRICHLET).  Subtracting this
! processor's first equation number and adding 1 gives the index into
! this processor's section of the solution, which is in recv_real

p = my_proc(procs)
do i=1,phaml_matrix%neq
   if (mumps_matrix%global_eq(i) == DIRICHLET) cycle
   if (.not. phaml_matrix%iown(i)) cycle
   phaml_matrix%solution(i) = recv_real(mumps_matrix%global_eq(i)-mumps_matrix%firsteq(p)+1)
end do

if (associated(recv_real)) deallocate(recv_real,stat=allocstat)

! exchange solution to get unowned points

if (.not.still_sequential .and. shadow) then
   call exchange_fudop_vect(phaml_matrix%solution(1:),procs,phaml_matrix, &
                            1510,1511,1512)
endif

end subroutine mumps_distrib_soln

end module mumps_interf
