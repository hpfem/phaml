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

module hypre_interf

!----------------------------------------------------
! This module contains the interface to the hypre package.
!
! communication tags in this module are of the form 18xx
!----------------------------------------------------

!----------------------------------------------------
! Other modules used are:

use global
use message_passing
use hypretype_mod
use hash_mod
use hash_eq_mod
use linsystype_mod
use linsys_util
use lapack_solve

!----------------------------------------------------
implicit none
private
public create_hypre_linear_system, destroy_hypre_linear_system, &
       change_hypre_rhs, zero_hypre_solution, hypre_solve

!----------------------------------------------------
! Parameters from hypre; make sure these agree with hypre's C include files

integer, parameter :: HYPRE_PARCSR = 5555

! These constants come from the Fortran wrappers for the SetPrecond
! routines for PCG and GMRES

integer, parameter :: NO_PRECOND_ID        = 0, &
                      DS_PRECOND_ID        = 1, &
                      AMG_PRECOND_ID       = 2, &
                      PARASAILS_PRECOND_ID = 4
!----------------------------------------------------

contains

!          --------------------------
subroutine create_hypre_linear_system(hypre_matrix,phaml_matrix,procs, &
                                      still_sequential)
!          --------------------------

!----------------------------------------------------
! This routine creates the linear system in hypre format
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hypre_matrix_type), intent(inout) :: hypre_matrix
type(linsys_type), intent(inout) :: phaml_matrix
type(proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential
!----------------------------------------------------
! Local variables:

integer :: my_processor, nproc
integer :: allocstat
!----------------------------------------------------
! Begin executable code

! useful information about the processors

my_processor = my_proc(procs)
nproc = num_proc(procs)

! master does not participate

if (my_processor == MASTER) return

! copy the matrix to hypre format

call make_hypre_matrix(hypre_matrix,phaml_matrix,procs,still_sequential)

! copy the right hand side to hypre format

call make_hypre_rhs(hypre_matrix,phaml_matrix,phaml_matrix%rhs,procs, &
                    still_sequential)

! copy the solution into the hypre solution vector as initial guess

call make_hypre_soln(hypre_matrix,phaml_matrix,procs,still_sequential)

end subroutine create_hypre_linear_system

!          ---------------------------
subroutine destroy_hypre_linear_system(hypre_matrix,procs)
!          ---------------------------

!----------------------------------------------------
! This routine deallocates memory associated with hypre_matrix
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments
   
type(hypre_matrix_type), intent(inout) :: hypre_matrix
type(proc_info), intent(in) :: procs
!----------------------------------------------------
! Local variables:

integer :: allocstat, herr
!----------------------------------------------------
! Begin executable code

! master does not participate

if (my_proc(procs) == MASTER) return

if (associated(hypre_matrix%firsteq)) &
   deallocate(hypre_matrix%firsteq,stat=allocstat)
if (associated(hypre_matrix%global_eq)) &
   deallocate(hypre_matrix%global_eq,stat=allocstat)
if (associated(hypre_matrix%lapack_rhs)) &
   deallocate(hypre_matrix%lapack_rhs,stat=allocstat)
if (hypre_matrix%ij_matrix /= -1) then
   call HYPRE_IJMatrixDestroy(hypre_matrix%ij_matrix,herr)
   call HYPRE_IJVectorDestroy(hypre_matrix%ij_rhs,herr)
   call HYPRE_IJVectorDestroy(hypre_matrix%ij_solution,herr)
endif

end subroutine destroy_hypre_linear_system

!          -----------------
subroutine make_hypre_matrix(hypre_matrix,phaml_matrix,procs, &
                             still_sequential)
!          -----------------

!----------------------------------------------------
! This routine copies the matrix to hypre format.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hypre_matrix_type), intent(inout) :: hypre_matrix
type(linsys_type), intent(inout) :: phaml_matrix
type(proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential
!----------------------------------------------------
! Local variables:

integer :: my_processor, nproc, comm
integer :: i, my_neq, not_own, allocstat, my_first_eq, counter, count_other, &
           p, lid, my_last_eq, herr, row, col, j, isub, scounter, &
           rcounter, oscounter, KEY_SIZE_EQ
integer, allocatable :: neq_all(:), owned_col(:), unowned_col(:), cols(:)
real(kind(0.0d0)), allocatable :: values(:)
! newcomm
integer :: nsend
integer, allocatable :: isend(:),nsendv(:),nrecv(:)
integer, pointer :: irecv(:)

!----------------------------------------------------
! Begin executable code

! useful information about the processors

comm = slaves_communicator(procs)
my_processor = my_proc(procs)
nproc = num_proc(procs)
KEY_SIZE_EQ = KEY_SIZE+1

! if still sequential, use LAPACK solver on each processor instead

if (still_sequential) then
   if (phaml_matrix%lapack_symm_band_exists) then
      call destroy_lapack_band(phaml_matrix%lapack_mat)
      phaml_matrix%lapack_symm_band_exists = .false.
   endif
   if (.not. phaml_matrix%lapack_gen_band_exists) then
      if (phaml_matrix%neq == phaml_matrix%neq_vert+phaml_matrix%neq_edge) then
         call make_lapack_gen_band(phaml_matrix%nlev+1,phaml_matrix, &
                                   phaml_matrix%lapack_mat)
      else
         call make_lapack_gen_band(phaml_matrix%nlev+2,phaml_matrix, &
                                   phaml_matrix%lapack_mat)
      endif
      phaml_matrix%lapack_gen_band_exists = .true.
   endif
   nullify(hypre_matrix%firsteq)
   nullify(hypre_matrix%global_eq)
   hypre_matrix%ij_matrix = -1
   return
endif

! allocate space for received messages
   
allocate(nsendv(nproc), nrecv(nproc), stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in make_hypre_matrix",procs=procs)
   return
endif

! allocate space to keep track of equation numbering

allocate(hypre_matrix%firsteq(nproc+1), &
         hypre_matrix%global_eq(phaml_matrix%neq),stat=allocstat)
nullify(hypre_matrix%lapack_rhs)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in make_hypre_matrix",procs=procs)
   return
endif

! edge and face equations that I don't own might not exist on the processor
! that owns the region containing them, because the edge or face has been
! refined.  Start global_eq at a negative number so that it doesn't fall
! in this processor's range if it doesn't get set by some processor

hypre_matrix%global_eq = -10

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

allocate(neq_all(nproc),stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in make_hypre_matrix",procs=procs)
   return
endif
neq_all = 0
neq_all(my_processor) = my_neq
neq_all = phaml_global_sum(procs,neq_all,1801)

! set the total number of equations

hypre_matrix%neq = sum(neq_all)

! set the global number of the first equation on each processor to be the
! sum of the number of equations on lower processor, plus 1.  Include
! an extra one so firsteq(nproc+1)-firsteq(nproc) works.

hypre_matrix%firsteq(1) = 1
do i=2,nproc+1
   hypre_matrix%firsteq(i) = hypre_matrix%firsteq(i-1) + neq_all(i-1)
end do

deallocate(neq_all,stat=allocstat)
my_first_eq = hypre_matrix%firsteq(my_processor)
my_last_eq  = hypre_matrix%firsteq(my_processor+1)-1

! determine the global number for each equation on this processor.
! First set the owned equations in the order they are found, using
! DIRICHLET (which is negative so it can't be a real equation number) for
! the Dirichlet points, and make a list of the gids of unowned equations.

allocate(isend(not_own*KEY_SIZE_EQ),stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in make_hypre_matrix",procs=procs)
   return
endif
counter = 0 
count_other = 1
do i=1,phaml_matrix%neq
   if (phaml_matrix%iown(i)) then
      if (phaml_matrix%equation_type(i) == DIRICHLET) then
         hypre_matrix%global_eq(i) = DIRICHLET
      else
         hypre_matrix%global_eq(i) = my_first_eq + counter
         counter = counter + 1
      endif
   else
      call hash_pack_key(phaml_matrix%gid(i),isend,count_other)
      count_other = count_other + KEY_SIZE_EQ
   endif
end do

! send the list of unowned gids to the other processors

call phaml_alltoall(procs,isend,count_other-1,irecv,nrecv,1802)

deallocate(isend, stat=allocstat)

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
   call fatal("allocation failed in make_hypre_matrix",procs=procs)
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
            isend(scounter+KEY_SIZE_EQ) = hypre_matrix%global_eq(lid)
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

call phaml_alltoall(procs,isend,nsendv,irecv,nrecv,1803)
deallocate(isend,stat=allocstat)

! copy the replies into global_eq

counter=1
do p=1,nproc
   if (p == my_processor) cycle
   do i=1,nrecv(p)/(KEY_SIZE_EQ+1)
      lid = hash_decode_key(hash_unpack_key(irecv,counter,.true.), &
                            phaml_matrix%eq_hash)
      if (lid == HASH_NOT_FOUND) then
         call warning("received reply for an equation I don't have in make_hypre_matrix")
      else
         hypre_matrix%global_eq(lid) = irecv(counter+KEY_SIZE_EQ)
      endif
      counter = counter + KEY_SIZE_EQ+1
   end do
end do

if (associated(irecv)) deallocate(irecv,stat=allocstat)
deallocate(nsendv,nrecv,stat=allocstat)

! count the number of owned and unowned entries in the owned rows

allocate(owned_col(my_neq),unowned_col(my_neq),stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in make_hypre_matrix",procs=procs)
   return
endif

owned_col = 0; unowned_col = 0

do i=1,phaml_matrix%neq
   row = hypre_matrix%global_eq(i)
   if (row < my_first_eq .or. row > my_last_eq) cycle ! not own or Dirichlet
   row = row - my_first_eq + 1
   do j=phaml_matrix%begin_row(i),phaml_matrix%end_row(i)
      if (phaml_matrix%column_index(j) == NO_ENTRY) cycle
      col = hypre_matrix%global_eq(phaml_matrix%column_index(j))
      if (col == DIRICHLET) cycle
      if (col < my_first_eq .or. col > my_last_eq) then
         unowned_col(row) = unowned_col(row) + 1
      else
         owned_col(row) = owned_col(row) + 1
      endif
   end do
end do

! create the hypre matrix

call my_HYPRE_IJMatrixCreate(comm, my_first_eq, my_last_eq, my_first_eq, &
                             my_last_eq, hypre_matrix%ij_matrix, herr)
call HYPRE_IJMatrixSetObjectType(hypre_matrix%ij_matrix, HYPRE_PARCSR, herr)
call HYPRE_IJMatrixSetRowSizes(hypre_matrix%ij_matrix, owned_col+unowned_col, &
                               herr)
call HYPRE_IJMatrixSetDiagOffdSizes(hypre_matrix%ij_matrix, owned_col, &
                                    unowned_col,herr)
call HYPRE_IJMatrixInitialize(hypre_matrix%ij_matrix, herr)

! copy the matrix information into the data structures that get passed
! to hypre to copy them into the hypre matrix

counter = 0
do i=1,phaml_matrix%neq
   row = hypre_matrix%global_eq(i)
   if (row < my_first_eq .or. row > my_last_eq) cycle ! not own or Dirichlet
   do j=phaml_matrix%begin_row(i),phaml_matrix%end_row(i)
      if (phaml_matrix%column_index(j) == NO_ENTRY) cycle
      col = hypre_matrix%global_eq(phaml_matrix%column_index(j))
      if (col == DIRICHLET) cycle
      counter = counter + 1
   end do
end do

allocate(cols(counter),values(counter),stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in make_hypre_matrix",procs=procs)
   return
endif

counter = 0
do i=1,phaml_matrix%neq
   row = hypre_matrix%global_eq(i)
   if (row < my_first_eq .or. row > my_last_eq) cycle ! not own or Dirichlet
   do j=phaml_matrix%begin_row(i),phaml_matrix%end_row(i)
      if (phaml_matrix%column_index(j) == NO_ENTRY) cycle
      col = hypre_matrix%global_eq(phaml_matrix%column_index(j))
      if (col == DIRICHLET) cycle
      counter = counter + 1
      cols(counter) = col
      values(counter) = phaml_matrix%matrix_val(j)
   end do
end do

! place the values in the hypre_matrix

call HYPRE_IJMatrixSetValues(hypre_matrix%ij_matrix, my_neq, &
                             owned_col+unowned_col, &
                             (/ (i,i=my_first_eq,my_last_eq) /), &
                             cols, values, herr)

deallocate(cols,values,stat=allocstat)
deallocate(owned_col,unowned_col,stat=allocstat)

! finalize the matrix

call HYPRE_IJMatrixAssemble(hypre_matrix%ij_matrix, herr)
call HYPRE_IJMatrixGetObject(hypre_matrix%ij_matrix, &
                             hypre_matrix%parcsr_matrix, herr)

end subroutine make_hypre_matrix

!          --------------
subroutine make_hypre_rhs(hypre_matrix,phaml_matrix,rs,procs,still_sequential)
!          --------------

!----------------------------------------------------
! This routine copies rs into the hypre data structure. 
! make_hypre_matrix must be called before this
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hypre_matrix_type), intent(inout) :: hypre_matrix
type(linsys_type), intent(in) :: phaml_matrix
real(my_real), intent(in) :: rs(:)
type(proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential
!----------------------------------------------------
! Local variables:

integer :: my_processor, comm
integer :: herr, i, j, my_first_eq, my_last_eq, allocstat, counter, eq
real(kind(1.0d0)), allocatable :: my_rhs(:)
!----------------------------------------------------
! Begin executable code

! if still_sequential, put the rhs in a holding place

if (still_sequential) then
   allocate(hypre_matrix%lapack_rhs(size(rs)),stat=allocstat)
   if (allocstat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in make_hypre_rhs",procs=procs)
      return
   endif
   hypre_matrix%lapack_rhs = rs
   return
endif

! useful information about the processors

comm = slaves_communicator(procs)
my_processor = my_proc(procs)

my_first_eq = hypre_matrix%firsteq(my_processor)
my_last_eq  = hypre_matrix%firsteq(my_processor+1)-1

! Create the rhs object

call my_HYPRE_IJVectorCreate(comm, my_first_eq, my_last_eq, hypre_matrix%ij_rhs, &
                          herr)
call HYPRE_IJVectorSetObjectType(hypre_matrix%ij_rhs, HYPRE_PARCSR, herr)
call HYPRE_IJVectorInitialize(hypre_matrix%ij_rhs, herr)

! Set the rhs of equations that this processor owns, including the
! removal of Dirichlet boundary conditions

allocate(my_rhs(my_last_eq-my_first_eq+1),stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in make_hypre_rhs",procs=procs)
   return
endif

counter = 0
do i=1,phaml_matrix%neq
   eq = hypre_matrix%global_eq(i)
   if (eq < my_first_eq .or. eq > my_last_eq) cycle ! not own or Dirichlet
   counter = counter + 1
   my_rhs(counter) = rs(i)
   do j=phaml_matrix%begin_row(i),phaml_matrix%end_row(i)
      if (phaml_matrix%column_index(j) == NO_ENTRY) cycle
      if (phaml_matrix%equation_type(phaml_matrix%column_index(j)) == DIRICHLET) then
         my_rhs(counter) = my_rhs(counter) - &
            phaml_matrix%matrix_val(j)*phaml_matrix%solution(phaml_matrix%column_index(j))
      endif
   end do
end do

! copy the rhs to the hypre matrix

call HYPRE_IJVectorSetValues(hypre_matrix%ij_rhs, size(my_rhs), &
                             (/ (i,i=my_first_eq,my_last_eq) /), my_rhs, herr)

deallocate(my_rhs,stat=allocstat)

! finalize the rhs

call HYPRE_IJVectorAssemble(hypre_matrix%ij_rhs, herr)
call HYPRE_IJVectorGetObject(hypre_matrix%ij_rhs, hypre_matrix%par_rhs, herr)

end subroutine make_hypre_rhs

!          ---------------
subroutine make_hypre_soln(hypre_matrix,phaml_matrix,procs,still_sequential)
!          ---------------

!----------------------------------------------------
! This routine copies phaml_matrix%solution into the hypre data structure. 
! make_hypre_matrix must be called before this
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hypre_matrix_type), intent(inout) :: hypre_matrix
type(linsys_type), intent(in) :: phaml_matrix
type(proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential
!----------------------------------------------------
! Local variables:

integer :: my_processor, comm
integer :: herr, i, j, my_first_eq, my_last_eq, allocstat, counter, eq
real(kind(1.0d0)), allocatable :: my_soln(:)
!----------------------------------------------------
! Begin executable code

! nothing to do if still sequential

if (still_sequential) return

! useful information about the processors

comm = slaves_communicator(procs)
my_processor = my_proc(procs)

my_first_eq = hypre_matrix%firsteq(my_processor)
my_last_eq  = hypre_matrix%firsteq(my_processor+1)-1

! Create the solution object

call my_HYPRE_IJVectorCreate(comm, my_first_eq, my_last_eq, hypre_matrix%ij_solution, &
                          herr)
call HYPRE_IJVectorSetObjectType(hypre_matrix%ij_solution, HYPRE_PARCSR, herr)
call HYPRE_IJVectorInitialize(hypre_matrix%ij_solution, herr)

! Set the solution of equations that this processor owns

allocate(my_soln(my_last_eq-my_first_eq+1),stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in make_hypre_soln",procs=procs)
   return
endif

counter = 0
do i=1,phaml_matrix%neq
   eq = hypre_matrix%global_eq(i)
   if (eq < my_first_eq .or. eq > my_last_eq) cycle ! not own or Dirichlet
   counter = counter + 1
   my_soln(counter) = phaml_matrix%solution(i)
end do

! copy the solution to the hypre matrix

call HYPRE_IJVectorSetValues(hypre_matrix%ij_solution, size(my_soln), &
                             (/ (i,i=my_first_eq,my_last_eq) /), my_soln, herr)

deallocate(my_soln,stat=allocstat)

! finalize the solution

call HYPRE_IJVectorAssemble(hypre_matrix%ij_solution, herr)
call HYPRE_IJVectorGetObject(hypre_matrix%ij_solution, hypre_matrix%par_solution, herr)

end subroutine make_hypre_soln

!          ----------------
subroutine change_hypre_rhs(hypre_matrix,phaml_matrix,rs,procs,still_sequential)
!          ----------------

!----------------------------------------------------
! This routine copies rs into the existing hypre data structure. 
! make_hypre_rhs must be called before this
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hypre_matrix_type), intent(inout) :: hypre_matrix
type(linsys_type), intent(in) :: phaml_matrix
real(my_real), intent(in) :: rs(:)
type(proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential
!----------------------------------------------------
! Local variables:

integer :: my_processor
integer :: herr, i, j, my_first_eq, my_last_eq, allocstat, counter, eq
real(kind(1.0d0)), allocatable :: my_rhs(:)
!----------------------------------------------------
! Begin executable code

! useful information about the processors

my_processor = my_proc(procs)
if (my_processor == MASTER) return

! if still sequential, copy rhs to a holding place allocated in make_rhs

if (still_sequential) then
   hypre_matrix%lapack_rhs = rs
   return
endif

my_first_eq = hypre_matrix%firsteq(my_processor)
my_last_eq  = hypre_matrix%firsteq(my_processor+1)-1

! Reinitialize the rhs object

call HYPRE_IJVectorInitialize(hypre_matrix%ij_rhs, herr)

! Set the rhs of equations that this processor owns, including the
! removal of Dirichlet boundary conditions

allocate(my_rhs(my_last_eq-my_first_eq+1),stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in make_hypre_rhs",procs=procs)
   return
endif

counter = 0
do i=1,phaml_matrix%neq
   eq = hypre_matrix%global_eq(i)
   if (eq < my_first_eq .or. eq > my_last_eq) cycle ! not own or Dirichlet
   counter = counter + 1
   my_rhs(counter) = rs(i)
   do j=phaml_matrix%begin_row(i),phaml_matrix%end_row(i)
      if (phaml_matrix%column_index(j) == NO_ENTRY) cycle
      if (phaml_matrix%equation_type(phaml_matrix%column_index(j)) == DIRICHLET) then
         my_rhs(counter) = my_rhs(counter) - &
            phaml_matrix%matrix_val(j)*phaml_matrix%solution(phaml_matrix%column_index(j))
      endif
   end do
end do

! copy the rhs to the hypre matrix

call HYPRE_IJVectorSetValues(hypre_matrix%ij_rhs, size(my_rhs), &
                             (/ (i,i=my_first_eq,my_last_eq) /), my_rhs, herr)

deallocate(my_rhs,stat=allocstat)

! finalize the rhs

call HYPRE_IJVectorAssemble(hypre_matrix%ij_rhs, herr)
call HYPRE_IJVectorGetObject(hypre_matrix%ij_rhs, hypre_matrix%par_rhs, herr)

end subroutine change_hypre_rhs

!          -------------------
subroutine zero_hypre_solution(hypre_matrix,procs,still_sequential)
!          -------------------

!----------------------------------------------------
! This routine zeroes out solution in the hypre data structure. 
! make_hypre_soln must be called before this
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hypre_matrix_type), intent(inout) :: hypre_matrix
type(proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential
!----------------------------------------------------
! Local variables:

integer :: my_processor
integer :: herr, i, my_first_eq, my_last_eq, allocstat
!----------------------------------------------------
! Begin executable code

! nothing to do if still sequential

if (still_sequential) return

! useful information about the processors

my_processor = my_proc(procs)
if (my_processor == MASTER) return

my_first_eq = hypre_matrix%firsteq(my_processor)
my_last_eq  = hypre_matrix%firsteq(my_processor+1)-1

! Reinitialize the solution object

call HYPRE_IJVectorInitialize(hypre_matrix%ij_solution, herr)

! set zeroes in the solution of the hypre matrix

call HYPRE_IJVectorSetValues(hypre_matrix%ij_solution, &
                             my_last_eq-my_first_eq+1, &
                             (/ (i,i=my_first_eq,my_last_eq) /), &
                             (/ (0.0d0,i=my_first_eq,my_last_eq) /), herr)

! finalize the solution

call HYPRE_IJVectorAssemble(hypre_matrix%ij_solution, herr)
call HYPRE_IJVectorGetObject(hypre_matrix%ij_solution, hypre_matrix%par_solution, herr)

end subroutine zero_hypre_solution

!          -----------
subroutine hypre_solve(hypre_matrix,phaml_matrix,procs,solver_cntl, &
                       still_sequential)
!          -----------

!----------------------------------------------------
! This routine solves the linear system using hypre
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hypre_matrix_type), intent(inout) :: hypre_matrix
type(linsys_type), intent(inout) :: phaml_matrix
type(solver_options), intent(in) :: solver_cntl
type(proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential

!----------------------------------------------------
! Local variables:

integer(hypre_pointer) :: solver, precond
integer :: comm, herr, my_processor, precond_id, allocstat
real(my_real), allocatable :: hold_rhs(:)
!----------------------------------------------------
! Begin executable code

my_processor = my_proc(procs)
if (my_processor == MASTER) return
comm = slaves_communicator(procs)

! if still sequential, use lapack solver

if (still_sequential) then
   allocate(hold_rhs(size(phaml_matrix%rhs)),stat=allocstat)
   if (allocstat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in hypre_solve",procs=procs)
      return
   endif
   hold_rhs = phaml_matrix%rhs
   phaml_matrix%rhs = hypre_matrix%lapack_rhs
   if (phaml_matrix%neq == phaml_matrix%neq_vert+phaml_matrix%neq_edge) then
     call lapack_indef(phaml_matrix%nlev+1,phaml_matrix,phaml_matrix%lapack_mat)
   else
     call lapack_indef(phaml_matrix%nlev+2,phaml_matrix,phaml_matrix%lapack_mat)
   endif
   phaml_matrix%rhs = hold_rhs
   deallocate(hold_rhs,stat=allocstat)
   return
endif

! select the solver to use

select case(solver_cntl%solver)

! BoomerAMG

case(HYPRE_BOOMERAMG_SOLVER)

! create solver

   call HYPRE_BoomerAMGCreate(solver, herr)

! set user provided options

   if (solver_cntl%hypre_cntl%BoomerAMG_MaxLevels /= huge(0)) then
      call HYPRE_BoomerAMGSetMaxLevels(solver, &
         solver_cntl%hypre_cntl%BoomerAMG_MaxLevels,herr)
   endif

   if (solver_cntl%hypre_cntl%BoomerAMG_MaxIter /= huge(0)) then
      call HYPRE_BoomerAMGSetMaxIter(solver, &
         solver_cntl%hypre_cntl%BoomerAMG_MaxIter,herr)
   endif

   if (solver_cntl%hypre_cntl%BoomerAMG_Tol /= huge(0.0d0)) then
      call HYPRE_BoomerAMGSetTol(solver, &
         solver_cntl%hypre_cntl%BoomerAMG_Tol,herr)
   endif

   if (solver_cntl%hypre_cntl%BoomerAMG_StrongThreshold /= huge(0.0d0)) then
      call HYPRE_BoomerAMGSetStrongThrshld(solver, &
         solver_cntl%hypre_cntl%BoomerAMG_StrongThreshold,herr)
   endif

   if (solver_cntl%hypre_cntl%BoomerAMG_MaxRowSum /= huge(0.0d0)) then
      call HYPRE_BoomerAMGSetMaxRowSum(solver, &
         solver_cntl%hypre_cntl%BoomerAMG_MaxRowSum,herr)
   endif

   if (solver_cntl%hypre_cntl%BoomerAMG_CoarsenType /= huge(0)) then
      call HYPRE_BoomerAMGSetCoarsenType(solver, &
         solver_cntl%hypre_cntl%BoomerAMG_CoarsenType,herr)
   endif

   if (solver_cntl%hypre_cntl%BoomerAMG_MeasureType /= huge(0)) then
      call HYPRE_BoomerAMGSetMeasureType(solver, &
         solver_cntl%hypre_cntl%BoomerAMG_MeasureType,herr)
   endif

   if (solver_cntl%hypre_cntl%BoomerAMG_CycleType /= huge(0)) then
      call HYPRE_BoomerAMGSetCycleType(solver, &
         solver_cntl%hypre_cntl%BoomerAMG_CycleType,herr)
   endif

   if (associated(solver_cntl%hypre_cntl%BoomerAMG_NumGridSweeps)) then
      call my_HYPRE_BoomerAMGSetNumGridSw(solver, &
         solver_cntl%hypre_cntl%BoomerAMG_NumGridSweeps,herr)
   endif

   if (associated(solver_cntl%hypre_cntl%BoomerAMG_GridRelaxType)) then
      call HYPRE_BoomerAMGSetGridRelaxType(solver, &
         solver_cntl%hypre_cntl%BoomerAMG_GridRelaxType,herr)
   endif

   if (associated(solver_cntl%hypre_cntl%BoomerAMG_GridRelaxPoints)) then
      call HYPRE_BoomerAMGSetGridRelaxPnts(solver, &
         solver_cntl%hypre_cntl%BoomerAMG_GridRelaxPoints,herr)
   endif

   if (associated(solver_cntl%hypre_cntl%BoomerAMG_RelaxWeight)) then
      call HYPRE_BoomerAMGSetRelaxWeight(solver, &
         solver_cntl%hypre_cntl%BoomerAMG_RelaxWeight,herr)
   endif

   if (solver_cntl%hypre_cntl%BoomerAMG_DebugFlag /= huge(0)) then
      call HYPRE_BoomerAMGSetDebugFlag(solver, &
         solver_cntl%hypre_cntl%BoomerAMG_DebugFlag,herr)
   endif

! solve

   call HYPRE_BoomerAMGSetup(solver, hypre_matrix%parcsr_matrix, &
                             hypre_matrix%par_rhs, hypre_matrix%par_solution, &
                             herr)
   call HYPRE_BoomerAMGSolve(solver, hypre_matrix%parcsr_matrix, &
                             hypre_matrix%par_rhs, hypre_matrix%par_solution, &
                             herr)

! get the solution

   call get_solution

! destroy solver

   call HYPRE_BoomerAMGDestroy(solver, herr)

! PCG

case(HYPRE_PCG_SOLVER)

! create solver

   call my_HYPRE_ParCSRPCGCreate(comm, solver, herr)

! set user provided options

   if (solver_cntl%hypre_cntl%PCG_Tol /= huge(0.0d0)) then
      call HYPRE_ParCSRPCGSetTol(solver, &
         solver_cntl%hypre_cntl%PCG_Tol,herr)
   endif

   if (solver_cntl%hypre_cntl%PCG_MaxIter /= huge(0)) then
      call HYPRE_ParCSRPCGSetMaxIter(solver, &
         solver_cntl%hypre_cntl%PCG_MaxIter,herr)
   endif

   if (solver_cntl%hypre_cntl%PCG_TwoNorm /= huge(0)) then
      call HYPRE_ParCSRPCGSetTwoNorm(solver, &
         solver_cntl%hypre_cntl%PCG_TwoNorm,herr)
   endif

   if (solver_cntl%hypre_cntl%PCG_RelChange /= huge(0)) then
      call HYPRE_ParCSRPCGSetRelChange(solver, &
         solver_cntl%hypre_cntl%PCG_RelChange,herr)
   endif

   if (solver_cntl%hypre_cntl%PCG_Logging /= huge(0)) then
      call HYPRE_ParCSRPCGSetLogging(solver, &
         solver_cntl%hypre_cntl%PCG_Logging,herr)
   endif

! set preconditioner

   call set_precond
   call HYPRE_ParCSRPCGSetPrecond(solver, precond_id, precond, herr)

! solve

   call HYPRE_ParCSRPCGSetup(solver, hypre_matrix%parcsr_matrix, &
                             hypre_matrix%par_rhs, hypre_matrix%par_solution, &
                             herr)
   call HYPRE_ParCSRPCGSolve(solver, hypre_matrix%parcsr_matrix, &
                             hypre_matrix%par_rhs, hypre_matrix%par_solution, &
                             herr)

! get the solution

   call get_solution

! destroy preconditioner and solver

   select case(solver_cntl%preconditioner)
   case(HYPRE_BOOMERAMG_PRECONDITION)
      call HYPRE_BoomerAMGDestroy(precond,herr)
   case(HYPRE_PARASAILS_PRECONDITION)
      call HYPRE_ParaSailsDestroy(precond,herr)
   end select

   call HYPRE_ParCSRPCGDestroy(solver, herr)

! GMRES

case(HYPRE_GMRES_SOLVER)

! create solver

   call my_HYPRE_ParCSRGMRESCreate(comm, solver, herr)

! set user provided options

   if (solver_cntl%hypre_cntl%GMRES_KDim /= huge(0)) then
      call HYPRE_ParCSRGMRESSetKDim(solver, &
         solver_cntl%hypre_cntl%GMRES_KDim,herr)
   endif

   if (solver_cntl%hypre_cntl%GMRES_Tol /= huge(0.0d0)) then
      call HYPRE_ParCSRGMRESSetTol(solver, &
         solver_cntl%hypre_cntl%GMRES_Tol,herr)
   endif

   if (solver_cntl%hypre_cntl%GMRES_MaxIter /= huge(0)) then
      call HYPRE_ParCSRGMRESSetMaxIter(solver, &
         solver_cntl%hypre_cntl%GMRES_MaxIter,herr)
   endif

   if (solver_cntl%hypre_cntl%GMRES_Logging /= huge(0)) then
      call HYPRE_ParCSRGMRESSetLogging(solver, &
         solver_cntl%hypre_cntl%GMRES_Logging,herr)
   endif

! set preconditioner

   call set_precond
   call HYPRE_ParCSRGMRESSetPrecond(solver, precond_id, precond, herr)

! solve

   call HYPRE_ParCSRGMRESSetup(solver, hypre_matrix%parcsr_matrix, &
                             hypre_matrix%par_rhs, hypre_matrix%par_solution, &
                             herr)
   call HYPRE_ParCSRGMRESSolve(solver, hypre_matrix%parcsr_matrix, &
                             hypre_matrix%par_rhs, hypre_matrix%par_solution, &
                             herr)

! get the solution

   call get_solution

! destroy preconditioner and solver

   select case(solver_cntl%preconditioner)
   case(HYPRE_BOOMERAMG_PRECONDITION)
      call HYPRE_BoomerAMGDestroy(precond,herr)
   case(HYPRE_PARASAILS_PRECONDITION)
      call HYPRE_ParaSailsDestroy(precond,herr)
   end select

   call HYPRE_ParCSRGMRESDestroy(solver, herr)

case default

   ierr = USER_INPUT_ERROR
   call fatal("Illegal choice of solver in hypre_solve",procs=procs)
   return

end select

return

contains

! internal subroutine to set preconditioner

subroutine set_precond

real(kind(0.0d0)) :: thresh
integer :: nlevels

select case(solver_cntl%preconditioner)

case(NO_PRECONDITION)

   precond_id = NO_PRECOND_ID

case(HYPRE_DS_PRECONDITION)

   precond_id = DS_PRECOND_ID

case(HYPRE_BOOMERAMG_PRECONDITION)

   precond_id = AMG_PRECOND_ID

   call HYPRE_BoomerAMGCreate(precond, herr)

   if (solver_cntl%hypre_cntl%BoomerAMG_MaxLevels /= huge(0)) then
      call HYPRE_BoomerAMGSetMaxLevels(precond, &
         solver_cntl%hypre_cntl%BoomerAMG_MaxLevels,herr)
   endif

   if (solver_cntl%hypre_cntl%BoomerAMG_MaxIter /= huge(0)) then
      call HYPRE_BoomerAMGSetMaxIter(precond, &
         solver_cntl%hypre_cntl%BoomerAMG_MaxIter,herr)
   endif

   if (solver_cntl%hypre_cntl%BoomerAMG_Tol /= huge(0.0d0)) then
      call HYPRE_BoomerAMGSetTol(precond, &
         solver_cntl%hypre_cntl%BoomerAMG_Tol,herr)
   endif

   if (solver_cntl%hypre_cntl%BoomerAMG_StrongThreshold /= huge(0.0d0)) then
      call HYPRE_BoomerAMGSetStrongThrshld(precond, &
         solver_cntl%hypre_cntl%BoomerAMG_StrongThreshold,herr)
   endif

   if (solver_cntl%hypre_cntl%BoomerAMG_MaxRowSum /= huge(0.0d0)) then
      call HYPRE_BoomerAMGSetMaxRowSum(precond, &
         solver_cntl%hypre_cntl%BoomerAMG_MaxRowSum,herr)
   endif

   if (solver_cntl%hypre_cntl%BoomerAMG_CoarsenType /= huge(0)) then
      call HYPRE_BoomerAMGSetCoarsenType(precond, &
         solver_cntl%hypre_cntl%BoomerAMG_CoarsenType,herr)
   endif

   if (solver_cntl%hypre_cntl%BoomerAMG_MeasureType /= huge(0)) then
      call HYPRE_BoomerAMGSetMeasureType(precond, &
         solver_cntl%hypre_cntl%BoomerAMG_MeasureType,herr)
   endif

   if (solver_cntl%hypre_cntl%BoomerAMG_CycleType /= huge(0)) then
      call HYPRE_BoomerAMGSetCycleType(precond, &
         solver_cntl%hypre_cntl%BoomerAMG_CycleType,herr)
   endif

   if (associated(solver_cntl%hypre_cntl%BoomerAMG_NumGridSweeps)) then
      call my_HYPRE_BoomerAMGSetNumGridSw(precond, &
         solver_cntl%hypre_cntl%BoomerAMG_NumGridSweeps,herr)
   endif

   if (associated(solver_cntl%hypre_cntl%BoomerAMG_GridRelaxType)) then
      call HYPRE_BoomerAMGSetGridRelaxType(precond, &
         solver_cntl%hypre_cntl%BoomerAMG_GridRelaxType,herr)
   endif

   if (associated(solver_cntl%hypre_cntl%BoomerAMG_GridRelaxPoints)) then
      call HYPRE_BoomerAMGSetGridRelaxPnts(precond, &
         solver_cntl%hypre_cntl%BoomerAMG_GridRelaxPoints,herr)
   endif

   if (associated(solver_cntl%hypre_cntl%BoomerAMG_RelaxWeight)) then
      call HYPRE_BoomerAMGSetRelaxWeight(precond, &
         solver_cntl%hypre_cntl%BoomerAMG_RelaxWeight,herr)
   endif

   if (solver_cntl%hypre_cntl%BoomerAMG_DebugFlag /= huge(0)) then
      call HYPRE_BoomerAMGSetDebugFlag(precond, &
         solver_cntl%hypre_cntl%BoomerAMG_DebugFlag,herr)
   endif

case(HYPRE_PARASAILS_PRECONDITION)

   precond_id = PARASAILS_PRECOND_ID

   call my_HYPRE_ParaSailsCreate(comm,precond,herr)

   if (solver_cntl%hypre_cntl%ParaSails_thresh /= huge(0.0d0) .or. &
       solver_cntl%hypre_cntl%ParaSails_nlevels /= huge(0)) then
      if (solver_cntl%hypre_cntl%ParaSails_thresh == huge(0.0d0)) then
         thresh = 0.1d0
      else
         thresh = solver_cntl%hypre_cntl%ParaSails_thresh
      endif
      if (solver_cntl%hypre_cntl%ParaSails_nlevels == huge(0)) then
         nlevels = 1
      else
         nlevels = solver_cntl%hypre_cntl%ParaSails_nlevels
      endif
      call HYPRE_ParaSailsSetParams(precond,thresh,nlevels,herr)
   endif

   if (solver_cntl%hypre_cntl%ParaSails_filter /= huge(0.0d0)) then
      call HYPRE_ParaSailsSetFilter(precond, &
         solver_cntl%hypre_cntl%ParaSails_filter,herr)
   endif

   if (solver_cntl%hypre_cntl%ParaSails_sym /= huge(0)) then
      call HYPRE_ParaSailsSetSym(precond, &
         solver_cntl%hypre_cntl%ParaSails_sym,herr)
   endif

! TEMP SetLoadBal and SetReuse are missing from the hypre Fortran interface

   if (solver_cntl%hypre_cntl%ParaSails_loadbal /= huge(0.0d0)) then
      call warning("SetLoadbal is missing from the hypre Fortran interface.", &
                   "Not changing its value.")
!      call HYPRE_ParaSailsSetLoadbal(precond, &
!         solver_cntl%hypre_cntl%ParaSails_loadbal,herr)
   endif

   if (solver_cntl%hypre_cntl%ParaSails_reuse /= huge(0)) then
      call warning("SetReuse is missing from the hypre Fortran interface.", &
                   "Not changing its value.")
!      call HYPRE_ParaSailsSetReuse(precond, &
!         solver_cntl%hypre_cntl%ParaSails_reuse,herr)
   endif

   if (solver_cntl%hypre_cntl%ParaSails_logging /= huge(0)) then
      call HYPRE_ParaSailsSetLogging(precond, &
         solver_cntl%hypre_cntl%ParaSails_logging,herr)
   endif

case default

   ierr = USER_INPUT_ERROR
   call fatal("Illegal choice for preconditioner for hypre.",procs=procs)
   return

end select

end subroutine set_precond

! internal routine to copy the hypre solution to phaml_matrix

subroutine get_solution

real(kind(0.0d0)), allocatable :: my_soln(:)
integer :: i, counter, my_first_eq, my_last_eq

my_first_eq = hypre_matrix%firsteq(my_processor)
my_last_eq  = hypre_matrix%firsteq(my_processor+1)-1

allocate(my_soln(my_last_eq-my_first_eq+1),stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in hypre_solve",procs=procs)
   return
endif

call HYPRE_IJVectorGetValues(hypre_matrix%ij_solution, size(my_soln), &
                             (/ (i,i=my_first_eq,my_last_eq) /), my_soln, herr)

counter = 0 
do i=1,phaml_matrix%neq
   if (phaml_matrix%iown(i)) then
      if (phaml_matrix%equation_type(i) /= DIRICHLET) then
         counter = counter + 1
         phaml_matrix%solution(i) = my_soln(counter)
      endif
   endif
end do

call exchange_fudop_vect(phaml_matrix%solution(1:),procs,phaml_matrix, &
                         1810,1811,1812)

end subroutine get_solution

end subroutine hypre_solve

end module hypre_interf
