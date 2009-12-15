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

module superlu_interf

!----------------------------------------------------
! This module contains the interface to the SuperLU package.
!
! communication tags in this module are of the form 19xx
!----------------------------------------------------

!----------------------------------------------------
! Other modules used are:

use global, AVOID_NATURAL=>NATURAL
use message_passing
use linsystype_mod
use superlutype_mod
use superlu_mod
use hash_mod
use hash_eq_mod
use linsys_util
use lapack_solve

!----------------------------------------------------
implicit none
private
public create_superlu_linear_system, destroy_superlu_linear_system, &
       change_superlu_rhs, superlu_solve

!----------------------------------------------------
! The following parameters are defined:

! These values come from superlu_defs.h.  If the values in there change with
! the version of SuperLU, then they need to be changed here, too.

integer, parameter :: NO                      = 0, & ! yes_no_t
                      YES                     = 1, &
                      DOFACT                  = 0, & ! fact_t
                      SamePattern             = 1, &
                      SamePattern_SameRowPerm = 2, &
                      FACTORED                = 3, &
                      NOROWPERM               = 0, & ! rowperm_t
                      LargeDiag               = 1, &
                      MY_PERMR                = 2, &
                      NATURAL                 = 0, & ! colperm_t
                      MMD_ATA                 = 1, &
                      MMD_AT_PLUS_A           = 2, &
                      COLAMD                  = 3, &
                      MY_PERMC                = 4, &
                      NOTRANS                 = 0, & ! trans_t
                      TRANS                   = 1, &
                      CONJ                    = 2, &
                      NOEQUIL                 = 0, & ! DiagScale_t  Need?
                      ROW                     = 1, &
                      COL                     = 2, &
                      BOTH                    = 3, &
                      NOREFINE                = 0, & ! IterRefine_t
                      DOUBLE                  = 1, &
                      EXTRA                   = 2, &
                      LUSUP                   = 0, & ! MemType  Need?
                      UCOL                    = 1, &
                      LSUB                    = 2, &
                      USUB                    = 3, &
                      SYSTEM                  = 0, & ! LU_space_t  Need?
                      USER                    = 1
integer, parameter :: SLU_NC                  = 0, & ! Stype_t
                      SLU_NR                  = 1, &
                      SLU_SC                  = 2, &
                      SLU_SR                  = 3, &
                      SLU_NCP                 = 4, &
                      SLU_DN                  = 5, &
                      SLU_NR_loc              = 6, &
                      SLU_S                   = 0, & ! Dtype_t
                      SLU_D                   = 1, &
                      SLU_C                   = 2, &
                      SLU_Z                   = 3, &
                      SLU_GE                  = 0, & ! Mtype_t
                      SLU_TRLU                = 1, &
                      SLU_TRUU                = 2, &
                      SLU_TRL                 = 3, &
                      SLU_TRU                 = 4, &
                      SLU_SYL                 = 5, &
                      SLU_SYU                 = 6, &
                      SLU_HEL                 = 7, &
                      SLU_HEU                 = 8


!----------------------------------------------------

contains

!          ----------------------------
subroutine create_superlu_linear_system(superlu_matrix, phaml_matrix, &
                                        procs, still_sequential)
!          ----------------------------

!----------------------------------------------------
! This routine initializes the SuperLU data structures, creates
! the matrix, and creates the rhs.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(superlu_matrix_type), intent(inout) :: superlu_matrix
type(linsys_type), intent(inout) :: phaml_matrix
type(proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential

!----------------------------------------------------
! Local variables:

integer :: nproc, nprow, npcol, comm, m, n
!----------------------------------------------------
! Begin executable code

! master does not participate

if (my_proc(procs) == MASTER) return

! If still sequential, just make the matrix and rhs

if (still_sequential) then
   call make_superlu_matrix(superlu_matrix,phaml_matrix,procs, &
                            still_sequential)
   call make_superlu_rhs(superlu_matrix,phaml_matrix,phaml_matrix%rhs,procs, &
                         still_sequential)
   return
endif

! create C structures in superlu_matrix_type

call f_create_gridinfo(superlu_matrix%grid)
call f_create_options(superlu_matrix%options)
call f_create_ScalePermstruct(superlu_matrix%ScalePermstruct)
call f_create_LUstruct(superlu_matrix%LUstruct)
call f_create_SOLVEstruct(superlu_matrix%SOLVEstruct)
call f_create_SuperMatrix(superlu_matrix%A)

! determine number of process rows and columns that use all slave processors
! and are as close to square as possible

nproc = num_proc(procs)
nprow = ceiling(sqrt(real(nproc)))
do
   npcol = nproc/nprow
   if (nprow*npcol == nproc) exit
   nprow = nprow - 1
end do

! initialize the SuperLU process grid

comm = slaves_communicator(procs)
call f_superlu_gridinit(comm, nprow, npcol, superlu_matrix%grid)

! make the matrix

call make_superlu_matrix(superlu_matrix,phaml_matrix,procs, &
                         still_sequential)

! make the rhs

call make_superlu_rhs(superlu_matrix,phaml_matrix,phaml_matrix%rhs,procs, &
                      still_sequential)

! set the default input options

call f_set_default_options_dist(superlu_matrix%options)

! don't want it to print statistics

call set_superlu_options(superlu_matrix%options,PrintStat=NO)

! initialize ScalePermstruct and LUstruct

call get_SuperMatrix(superlu_matrix%A,nrow=m,ncol=n)

call f_ScalePermstructInit(m, n, superlu_matrix%ScalePermstruct)
call f_LUstructInit(m, n, superlu_matrix%LUstruct)

end subroutine create_superlu_linear_system

!          -------------
subroutine superlu_solve(superlu_matrix,phaml_matrix,procs, &
                         still_sequential)
!          -------------

!----------------------------------------------------
! This routine solves the linear system using SuperLU
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(superlu_matrix_type), intent(inout) :: superlu_matrix
type(linsys_type), intent(inout) :: phaml_matrix
type(proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential
!----------------------------------------------------
! Local variables:

integer(superlu_ptr) :: stat ! SuperLUStat_t
real(dp) :: berr(1)
integer :: info, p, i, nr, ni, allocstat
integer, pointer :: recv_int(:)
real(my_real), pointer :: recv_real(:)
real(my_real), allocatable :: hold_rhs(:)
!----------------------------------------------------
! Begin executable code

! if still sequential, use lapack solver

if (still_sequential) then
   allocate(hold_rhs(size(phaml_matrix%rhs)),stat=allocstat)
   if (allocstat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in superlu_solve",procs=procs)
      return
   endif
   hold_rhs = phaml_matrix%rhs
   phaml_matrix%rhs = superlu_matrix%lapack_rhs
   if (phaml_matrix%neq == phaml_matrix%neq_vert+phaml_matrix%neq_edge) then
     call lapack_indef(phaml_matrix%nlev+1,phaml_matrix,phaml_matrix%lapack_mat)
   else
     call lapack_indef(phaml_matrix%nlev+2,phaml_matrix,phaml_matrix%lapack_mat)
   endif
   phaml_matrix%rhs = hold_rhs
   deallocate(hold_rhs,stat=allocstat)
   return
endif

! initialize the statistics variables

call f_create_SuperLUStat(stat)
call f_PStatInit(stat)

! call the linear equation solver

call f_pdgssvx(superlu_matrix%options, superlu_matrix%A, &
               superlu_matrix%ScalePermstruct, superlu_matrix%b, &
               superlu_matrix%my_neq, 1, superlu_matrix%grid, &
               superlu_matrix%LUstruct, superlu_matrix%SOLVEstruct, berr, &
               stat, info)

if (info /= 0) then
   call warning("SuperLU solution returned nonzero info",intlist=(/info/))
endif

! indicate that the matrix is now factored

call set_superlu_options(superlu_matrix%options,Fact=FACTORED)

! copy the solution into phaml_matrix solution.  Start by copying the
! non-Dirichlet values for points this processor owns.

p = my_proc(procs)
if (p == 1 .or. .not. still_sequential) then
   do i=1,phaml_matrix%neq
      if (superlu_matrix%global_eq(i) == DIRICHLET) cycle
      if (.not. phaml_matrix%iown(i)) cycle
      phaml_matrix%solution(i)=superlu_matrix%b(superlu_matrix%global_eq(i) - &
                                                superlu_matrix%firsteq(p)+1,1)
   end do
endif

! if still sequential, send solution from proc 1 to the others.
! if not, exchange solutions

if (still_sequential) then

   if (my_proc(procs) == 1) then
      do p=2,num_proc(procs)
         call phaml_send(procs,p,(/0/),0,phaml_matrix%solution, &
                         size(phaml_matrix%solution),1915)
      end do
   else
      call phaml_recv(procs,p,recv_int,ni,recv_real,nr,1915)
      phaml_matrix%solution = recv_real(1:size(phaml_matrix%solution))
      deallocate(recv_real)
   endif
else

   call exchange_fudop_vect(phaml_matrix%solution(1:),procs,phaml_matrix, &
                            1910,1911,1912)

endif

! free memory

call f_PStatFree(stat)
call f_destroy_SuperLUStat(stat)

end subroutine superlu_solve

!          -----------------------------
subroutine destroy_superlu_linear_system(superlu_matrix,procs)
!          -----------------------------

!----------------------------------------------------
! This routine deallocates all memory associated with superlu_matrix.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(superlu_matrix_type), intent(inout) :: superlu_matrix
type(proc_info), intent(in) :: procs

!----------------------------------------------------
! Local variables:

integer :: n, init, allocstat
!----------------------------------------------------
! Begin executable code

if (my_proc(procs) == MASTER) return

! If lapack_rhs is allocated, then this would have been created under
! still_sequential, and deallocating lapack_rhs is all that should be done

if (associated(superlu_matrix%lapack_rhs)) then
   deallocate(superlu_matrix%lapack_rhs, stat=allocstat)
   return
endif

! deallocate SuperLU allocated storage

! Do not Destroy_CompRowLoc_Matrix; it deallocates rowptr, colind and matval
!call f_Destroy_CompRowLoc_Matrix_dis(superlu_matrix%A)
call f_ScalePermstructFree(superlu_matrix%ScalePermstruct)
call get_SuperMatrix(superlu_matrix%A,ncol=n)
! TEMP Destroy_LU was crashing in free
!call f_Destroy_LU(n, superlu_matrix%grid, superlu_matrix%LUstruct)
call f_LUstructFree(superlu_matrix%LUstruct)
call get_superlu_options(superlu_matrix%options, SolveInitialized=init)
if (init == YES) then
   call f_dSolveFinalize(superlu_matrix%options, superlu_matrix%SOLVEstruct)
endif

! deallocate Fortran memory in superlu_matrix

if (associated(superlu_matrix%firsteq)) &
   deallocate(superlu_matrix%firsteq, stat=allocstat)
if (associated(superlu_matrix%global_eq)) &
   deallocate(superlu_matrix%global_eq, stat=allocstat)
if (associated(superlu_matrix%rowptr)) &
   deallocate(superlu_matrix%rowptr, stat=allocstat)
if (associated(superlu_matrix%colind)) &
   deallocate(superlu_matrix%colind, stat=allocstat)
if (associated(superlu_matrix%matval)) &
   deallocate(superlu_matrix%matval, stat=allocstat)
if (associated(superlu_matrix%b)) &
   deallocate(superlu_matrix%b, stat=allocstat)

! release the SuperLU process grid

call f_superlu_gridexit(superlu_matrix%grid)

! destroy C structures in superlu_matrix_type

call f_destroy_gridinfo(superlu_matrix%grid)
call f_destroy_options(superlu_matrix%options)
call f_destroy_ScalePermstruct(superlu_matrix%ScalePermstruct)
call f_destroy_LUstruct(superlu_matrix%LUstruct)
call f_destroy_SOLVEstruct(superlu_matrix%SOLVEstruct)
call f_destroy_SuperMatrix(superlu_matrix%A)

end subroutine destroy_superlu_linear_system

!          -------------------
subroutine make_superlu_matrix(superlu_matrix,phaml_matrix,procs, &
                               still_sequential)
!          -------------------

!----------------------------------------------------
! This routine copies the matrix to SuperLU format in superlu_matrix.
! Also returns the global equation number assigned to each equation.  Dirichlet
! points get equation number DIRICHLET.
! Also returns the first equation, in global numbering, on each processor.
! This should be dimensioned at nproc+1, as an extra entry gets N+1 so
! the section for a processor is given by firsteq(p):firsteq(p+1)-1.
! When still_sequential, the "distributed" matrix resides entirely on proc 1.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(superlu_matrix_type), intent(inout) :: superlu_matrix
type(linsys_type), intent(inout) :: phaml_matrix
type(proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential
!----------------------------------------------------
! Local variables:

integer :: i, nproc, my_processor, my_neq, my_first_eq, counter, count_other, &
           not_own, allocstat, p, lid, j, col, total_neq, nnz_loc, row, isub, &
           scounter, rcounter, oscounter, KEY_SIZE_EQ
integer, allocatable :: neq_all(:)
! newcomm
integer, allocatable :: isend(:),nsendv(:),nrecv(:)
integer, pointer :: irecv(:)

! Interface block for SuperLU's Create routine with the arrays whose
! addresses are being saved declared as targets, to prevent a
! compiler from passing them by value/return instead of by reference.
! Also, use a local array for the allocation rather than allocating
! directly in the structure.  Note, however, this prevents having more
! than one SuperLU matrix at a time.

real(dp), pointer, save :: matval(:)
integer, pointer, save :: colind(:), rowptr(:)

interface
subroutine f_dCreate_CompRowLoc_Matrix_dis(A, m, n, nnz_loc, m_loc, fst_row, &
                                    nzval, colind, rowptr, stype, dtype, mtype)
use superlutype_mod
integer(superlu_ptr) :: A
integer :: m, n, nz_loc, m_loc, fst_row
real(dp), target :: nzval(nnz_loc)
integer, target :: colind(nnz_loc), rowptr(m_loc+1)
integer :: stype, dtype, mtype
end subroutine f_dCreate_CompRowLoc_Matrix_dis
end interface

!----------------------------------------------------
! Begin executable code

nproc = num_proc(procs)
my_processor = my_proc(procs)
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
         call make_lapack_gen_band(phaml_matrix%nlev+1,phaml_matrix, &
                                   phaml_matrix%lapack_mat)
      endif
      phaml_matrix%lapack_gen_band_exists = .true.
   endif
   nullify(superlu_matrix%firsteq)
   nullify(superlu_matrix%global_eq)
   nullify(superlu_matrix%rowptr)
   nullify(superlu_matrix%colind)
   nullify(superlu_matrix%matval)
   nullify(superlu_matrix%b)
   return
endif

! allocate space for received messages

allocate(nsendv(nproc), nrecv(nproc), stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in make_superlu_matrix",procs=procs)
   return
endif

! allocate space to keep track of equation numbering

allocate(superlu_matrix%firsteq(nproc+1), &
         superlu_matrix%global_eq(phaml_matrix%neq),stat=allocstat)
nullify(superlu_matrix%lapack_rhs)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in make_superlu_matrix",procs=procs)
   return
endif

! if still sequential, processors other than 1 set up an empty matrix

if (still_sequential) then
   if (my_processor /= 1) then
      counter = 0
      do i=1,phaml_matrix%neq
         if (phaml_matrix%equation_type(i) == DIRICHLET) then
            superlu_matrix%global_eq(i) = DIRICHLET
         else
            counter = counter+1
            superlu_matrix%global_eq(i) = counter
         endif
      end do
      superlu_matrix%firsteq = counter+1
      superlu_matrix%firsteq(1) = 1
      superlu_matrix%my_neq = 0
      allocate(superlu_matrix%rowptr(1),superlu_matrix%colind(1), &
               superlu_matrix%matval(1),stat=allocstat)
      if (allocstat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in make_superlu_matrix",procs=procs)
         return
      endif
      call f_dCreate_CompRowLoc_Matrix_dis(superlu_matrix%A, counter, &
               counter, 0, 0, counter+1, superlu_matrix%matval, &
               superlu_matrix%colind, superlu_matrix%rowptr, SLU_NR_loc, &
               SLU_D, SLU_GE);

      return
   endif
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
superlu_matrix%my_neq = my_neq

! determine how many equations are owned by each processor

allocate(neq_all(nproc))
neq_all = 0
neq_all(my_processor) = my_neq
if (.not. still_sequential) then
   neq_all = phaml_global_sum(procs,neq_all,1901)
endif

! set the total number of equations

total_neq = sum(neq_all)

! set the global number of the first equation on each processor to be the
! sum of the number of equations on lower processor, plus 1.  Include
! an extra one so firsteq(nproc+1)-firsteq(nproc) works.
! Special case for still sequential, when only processor 1 gets here.

if (still_sequential) then
   superlu_matrix%firsteq = 1
   superlu_matrix%firsteq(my_processor+1) = total_neq+1
else
   superlu_matrix%firsteq(1) = 1
   do i=2,nproc+1
      superlu_matrix%firsteq(i) = superlu_matrix%firsteq(i-1) + neq_all(i-1)
   end do
endif

deallocate(neq_all,stat=allocstat)
my_first_eq = superlu_matrix%firsteq(my_processor)

! determine the global number for each equation on this processor.
! First set the owned equations in the order they are found, using
! DIRICHLET (which is negative so it can't be a real equation number) for
! the Dirichlet points, and make a list of the gids of unowned equations.

allocate(isend(not_own*KEY_SIZE_EQ),stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in make_superlu_matrix",procs=procs)
   return
endif
counter = 0
count_other = 1
do i=1,phaml_matrix%neq
   if (phaml_matrix%iown(i)) then
      if (phaml_matrix%equation_type(i) == DIRICHLET) then
         superlu_matrix%global_eq(i) = DIRICHLET
      else
         superlu_matrix%global_eq(i) = my_first_eq + counter
         counter = counter + 1
      endif
   else
      call hash_pack_key(phaml_matrix%gid(i),isend,count_other)
      count_other = count_other + KEY_SIZE_EQ
   endif
end do

! send the list of unowned gids to the other processors

if (.not. still_sequential) then
   call phaml_alltoall(procs,isend,count_other-1,irecv,nrecv,1902)
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
      call fatal("allocation failed in make_superlu_matrix",procs=procs)
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
               isend(scounter+KEY_SIZE_EQ) = superlu_matrix%global_eq(lid)
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

   call phaml_alltoall(procs,isend,nsendv,irecv,nrecv,1903)
   deallocate(isend,stat=allocstat)

! copy the replies into global_eq

   counter=1
   do p=1,nproc
      if (p == my_processor) cycle
      do i=1,nrecv(p)/(KEY_SIZE_EQ+1)
         lid = hash_decode_key(hash_unpack_key(irecv,counter,.true.), &
                               phaml_matrix%eq_hash)
         if (lid == HASH_NOT_FOUND) then
            call warning("received reply for an equation I don't have in make_superlu_matrix")
         else
            superlu_matrix%global_eq(lid) = irecv(counter+KEY_SIZE_EQ)
         endif
         counter = counter + KEY_SIZE_EQ+1
      end do
   end do

   if (associated(irecv)) deallocate(irecv,stat=allocstat)
   deallocate(nsendv,nrecv,stat=allocstat)

endif ! .not. still_sequential

! count the number of nonzeros in the rows owned by this processor, omitting
! Dirichlet points

counter = 0
do i=1,phaml_matrix%neq
   if (.not. phaml_matrix%iown(i)) cycle
   if (phaml_matrix%equation_type(i) == DIRICHLET) cycle
   do j=phaml_matrix%begin_row(i),phaml_matrix%end_row(i)
      if (phaml_matrix%column_index(j) == NO_ENTRY) cycle
      col = superlu_matrix%global_eq(phaml_matrix%column_index(j))
      if (col == DIRICHLET) cycle
      counter = counter + 1
   end do
end do

nnz_loc = counter

! copy the rows owned by this processor to superlu_matrix matrix structure,

allocate(rowptr(my_neq+1),colind(counter),matval(counter),stat=allocstat)

if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in make_superlu_matrix",procs=procs)
   return
endif

superlu_matrix%rowptr => rowptr
superlu_matrix%colind => colind
superlu_matrix%matval => matval

counter = 0
row = 0
do i=1,phaml_matrix%neq
   if (.not. phaml_matrix%iown(i)) cycle
   if (phaml_matrix%equation_type(i) == DIRICHLET) cycle
   row = row+1
   superlu_matrix%rowptr(row) = counter+1
   do j=phaml_matrix%begin_row(i),phaml_matrix%end_row(i)
      if (phaml_matrix%column_index(j) == NO_ENTRY) cycle
      col = superlu_matrix%global_eq(phaml_matrix%column_index(j))
      if (col == DIRICHLET) cycle
      counter = counter + 1
      superlu_matrix%colind(counter) = col
      superlu_matrix%matval(counter) = phaml_matrix%matrix_val(j)
   end do
end do
row = row+1
superlu_matrix%rowptr(row) = counter+1

! SuperLU has 0 based arrays, so all indices need to be decremented

superlu_matrix%rowptr = superlu_matrix%rowptr - 1
superlu_matrix%colind = superlu_matrix%colind - 1
my_first_eq = my_first_eq - 1

! create the SuperLU compressed row distributed matrix

call f_dCreate_CompRowLoc_Matrix_dis(superlu_matrix%A, total_neq, total_neq, &
               nnz_loc, my_neq, my_first_eq, matval, &
               colind, rowptr, SLU_NR_loc, &
               SLU_D, SLU_GE);

end subroutine make_superlu_matrix

!          ----------------
subroutine make_superlu_rhs(superlu_matrix,phaml_matrix,rs,procs,still_sequential)
!          ----------------

!----------------------------------------------------
! This routine creates the rhs vector and places rs in it
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(superlu_matrix_type), intent(inout) :: superlu_matrix
type(linsys_type), intent(in) :: phaml_matrix
real(my_real), intent(in) :: rs(:)
type(proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential
!----------------------------------------------------
! Local variables:

integer :: allocstat
!----------------------------------------------------
! Begin executable code

if (my_proc(procs) == MASTER) return

! if still_sequential, put the rhs in a holding place

if (still_sequential) then
   allocate(superlu_matrix%lapack_rhs(size(rs)),stat=allocstat)
   if (allocstat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in make_superlu_rhs",procs=procs)
      return
   endif
   superlu_matrix%lapack_rhs = rs
   return
endif

! allocate memory for b

allocate(superlu_matrix%b(superlu_matrix%my_neq,1),stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in make_superlu_rhs",procs=procs)
   return
endif

! setting the values is the same as changing the values

call change_superlu_rhs(superlu_matrix,phaml_matrix,rs,procs,still_sequential)

end subroutine make_superlu_rhs

!          ------------------
subroutine change_superlu_rhs(superlu_matrix,phaml_matrix,rs,procs, &
                              still_sequential)
!          ------------------

!----------------------------------------------------
! This routine copies rs to superlu_matrix%b.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(superlu_matrix_type), intent(inout) :: superlu_matrix
type(linsys_type), intent(in) :: phaml_matrix
real(my_real), intent(in) :: rs(:)
type(proc_info), intent(in) :: procs
logical, intent(in) :: still_sequential
!----------------------------------------------------
! Local variables:

integer :: row, i, j
!----------------------------------------------------
! Begin executable code

if (my_proc(procs) == MASTER) return

! if still sequential, copy rhs to a holding place allocated in make_rhs

if (still_sequential) then
   superlu_matrix%lapack_rhs = rs
   return
endif

if (my_proc(procs) == 1 .or. .not. still_sequential) then
   row = 0
   do i=1,phaml_matrix%neq
      if (.not. phaml_matrix%iown(i)) cycle
      if (phaml_matrix%equation_type(i) == DIRICHLET) cycle
      row = row+1
      superlu_matrix%b(row,1) = rs(i)
      do j=phaml_matrix%begin_row(i),phaml_matrix%end_row(i)
         if (phaml_matrix%column_index(j) == NO_ENTRY) cycle
         if (phaml_matrix%equation_type(phaml_matrix%column_index(j)) == DIRICHLET) then
            superlu_matrix%b(row,1) = superlu_matrix%b(row,1) - &
               phaml_matrix%matrix_val(j)*phaml_matrix%solution(phaml_matrix%column_index(j))
         endif
      end do
   end do
endif

end subroutine change_superlu_rhs

end module superlu_interf
