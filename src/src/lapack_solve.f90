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

module lapack_solve

!----------------------------------------------------
! This module contains routines for solving linear systems with LAPACK
!
! communication tags in this module are of the form 17xx
!----------------------------------------------------

!----------------------------------------------------
! Other modules used are:

use global
use message_passing
use linsystype_mod
use linsys_util

!----------------------------------------------------

implicit none
private
public make_lapack_symm_band, make_lapack_gen_band, destroy_lapack_band, &
       lapack_spd, lapack_indef, lapack_precon

contains

!          ---------------------
subroutine make_lapack_symm_band(use_nlev,phaml_matrix,lapack_matrix)
!          ---------------------

!----------------------------------------------------
! This routine makes and factorizes a LAPACK symmetric band matrix from the
! first use_nlev refinement levels of the PHAML matrix, which should be in
! nodal form for level use_nlev.  Dirichlet equations are omitted; make
! sure they are accounted for when setting up the right hand side.
! halfbandwidth does not include the diagonal.  lapack_matrix%rhs is
! allocated but not set.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: use_nlev
type(linsys_type), intent(in) :: phaml_matrix
type(lapack_band_matrix), intent(out) :: lapack_matrix

!----------------------------------------------------
! Local variables:

integer :: i, j, k, num_nonzero, workd, astat, counter, space, jerr, profil, &
           orig_bandwd, row, col, info
integer :: nodir_renum(phaml_matrix%neq), nodir_inv_renum(phaml_matrix%neq), &
           degree(phaml_matrix%neq), rstart(phaml_matrix%neq), &
           minband_renum(phaml_matrix%neq)
integer, allocatable :: connec(:), work(:)
!----------------------------------------------------
! Begin executable code

! count the number of equations without Dirichlet points, and the number
! of nonzeroes in the non-Dirichlet rows of mat, and create a renumbering
! that omits Dirichlet points.

lapack_matrix%neq = 0
num_nonzero = 0
nodir_renum = -1
do i=1,use_nlev
   do j=phaml_matrix%begin_level(i),phaml_matrix%begin_level(i+1)-1
      if (phaml_matrix%equation_type(j) == DIRICHLET) cycle
      lapack_matrix%neq = lapack_matrix%neq + 1
      nodir_renum(j) = lapack_matrix%neq
      nodir_inv_renum(lapack_matrix%neq) = j
      do k=phaml_matrix%begin_row(j),phaml_matrix%end_row(j)
         if (phaml_matrix%column_index(k) == NO_ENTRY) cycle
         if (phaml_matrix%equation_type(phaml_matrix%column_index(k)) == DIRICHLET) cycle
         if (phaml_matrix%column_index(k) >= phaml_matrix%begin_level(use_nlev+1)) cycle
         num_nonzero = num_nonzero + 1
      end do
   end do
end do

! nullify all array components of lapack_matrix

nullify(lapack_matrix%matrix,lapack_matrix%rhs,lapack_matrix%renum, &
        lapack_matrix%inv_renum,lapack_matrix%ipiv)

! if no equations, never use lapack_matrix so don't continue

if (lapack_matrix%neq == 0) return

! create column pointers for the matrix with Dirichlet points removed

allocate(connec(num_nonzero-lapack_matrix%neq),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in make_lapack_symm_band")
   return
endif

orig_bandwd = 0
counter = 0
do i=1,lapack_matrix%neq
   j = nodir_inv_renum(i)
   rstart(i) = counter + 1
   do k=phaml_matrix%begin_row(j)+1,phaml_matrix%end_row(j)
      if (phaml_matrix%column_index(k) == NO_ENTRY) cycle
      if (phaml_matrix%equation_type(phaml_matrix%column_index(k)) == DIRICHLET) cycle
      if (phaml_matrix%column_index(k) >= phaml_matrix%begin_level(use_nlev+1)) cycle
      counter = counter + 1
      connec(counter) = nodir_renum(phaml_matrix%column_index(k))
      orig_bandwd = max(orig_bandwd,abs(connec(counter)-i))
   end do
   degree(i) = counter - rstart(i) + 1
end do

! Use CALGO 582 (Gibbs, Poole and Stockmeyer) to find a bandwidth reduction
! ordering

workd = 6*lapack_matrix%neq+3
allocate(work(workd),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in make_lapack_symm_band")
   return
endif
minband_renum = (/ (i,i=1,phaml_matrix%neq) /)
call gpskca(lapack_matrix%neq,degree,rstart,connec,.false.,workd,minband_renum,&
            work,lapack_matrix%halfbandwidth,profil,jerr,space)
if (jerr /= 0) then
   work(1:lapack_matrix%neq) = (/ (i,i=1,lapack_matrix%neq) /) ! contains minband_inv_renum
   minband_renum(1:lapack_matrix%neq) = (/ (i,i=1,lapack_matrix%neq) /)
   lapack_matrix%halfbandwidth = orig_bandwd
   call warning("Bandwidth reduction reordering routine failed.", &
                "Rediculously large bandwidth may cause allocation failure or long run time.")
endif

! done with column pointers for reduced array

deallocate(connec)

! compose nodir and minband renumberings

! renumbering and inverse renumbering of the equations.  renum(i) is the
! equation number in the symmetric matrix corresponding to original equation i;
! dimension renum(neq(use_nlev)); renum(i) = -1 if i is a Dirichlet point.
! inv_renum(i) is the original equation number corresponding to the symmetric
! matrix equation i; dimension inv_renum(symm_neq).

allocate(lapack_matrix%renum(phaml_matrix%neq), &
         lapack_matrix%inv_renum(lapack_matrix%neq),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in make_lapack_symm_band")
   return
endif

do i=1,phaml_matrix%neq
   if (nodir_renum(i) == -1) then
      lapack_matrix%renum(i) = -1
   else
      lapack_matrix%renum(i) = minband_renum(nodir_renum(i))
   endif
end do

do i=1,lapack_matrix%neq
   lapack_matrix%inv_renum(i) = nodir_inv_renum(work(i))
end do

! done with work space for gpskca

deallocate(work)

! Copy matrix values to symmetric band form, upper triangle version

allocate(lapack_matrix%matrix(lapack_matrix%halfbandwidth+1,lapack_matrix%neq),&
         stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in make_lapack_symm_band")
   return
endif
lapack_matrix%matrix = 0.0_my_real

do i=1,use_nlev
   do j=phaml_matrix%begin_level(i),phaml_matrix%begin_level(i+1)-1
      if (phaml_matrix%equation_type(j) == DIRICHLET) cycle
      row = lapack_matrix%renum(j)
      if (row == -1) cycle
      do k=phaml_matrix%begin_row(j),phaml_matrix%end_row(j)
         if (phaml_matrix%column_index(k) == NO_ENTRY) cycle
         if (phaml_matrix%equation_type(phaml_matrix%column_index(k)) == DIRICHLET) cycle
         if (phaml_matrix%column_index(k) >= phaml_matrix%begin_level(use_nlev+1)) cycle
         col = lapack_matrix%renum(phaml_matrix%column_index(k))
         if (col == -1) cycle
         if (col < row) cycle
         lapack_matrix%matrix(lapack_matrix%halfbandwidth+1+row-col,col) = &
            phaml_matrix%matrix_val(k)
      end do
   end do
end do

! Factor the matrix

if (my_real == kind(0.0e0)) then
   call spbtrf("U",lapack_matrix%neq,lapack_matrix%halfbandwidth, &
               lapack_matrix%matrix,lapack_matrix%halfbandwidth+1,info)
elseif (my_real == kind(0.0d0)) then
   call dpbtrf("U",lapack_matrix%neq,lapack_matrix%halfbandwidth, &
               lapack_matrix%matrix,lapack_matrix%halfbandwidth+1,info)
else
   ierr = PHAML_INTERNAL_ERROR
   call fatal("my_real is neither default single nor double precision")
   return
endif
if (info /= 0) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("LAPACK SPBTRF failed.",intlist=(/info/))
   return
endif

! allocate the rhs

allocate(lapack_matrix%rhs(lapack_matrix%neq,1),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in make_lapack_symm_band")
   return
endif

end subroutine make_lapack_symm_band

!          --------------------
subroutine make_lapack_gen_band(use_nlev,phaml_matrix,lapack_matrix)
!          --------------------

!----------------------------------------------------
! This routine makes and factorizes a LAPACK general band matrix from the
! first use_nlev refinement levels of the PHAML matrix, which should be in
! nodal form for level use_nlev.  Dirichlet equations are omitted; make
! sure they are accounted for when setting up the right hand side.
! halfbandwidth does not include the diagonal.  lapack_matrix%rhs is
! allocated but not set.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: use_nlev
type(linsys_type), intent(in) :: phaml_matrix
type(lapack_band_matrix), intent(out) :: lapack_matrix

!----------------------------------------------------
! Local variables:

integer :: i, j, k, num_nonzero, workd, astat, counter, space, jerr, profil, &
           orig_bandwd, row, col, info
integer :: nodir_renum(phaml_matrix%neq), nodir_inv_renum(phaml_matrix%neq), &
           degree(phaml_matrix%neq), rstart(phaml_matrix%neq), &
           minband_renum(phaml_matrix%neq)
integer, allocatable :: connec(:), work(:)
!----------------------------------------------------
! Begin executable code

! count the number of equations without Dirichlet points, and the number
! of nonzeroes in the non-Dirichlet rows of mat, and create a renumbering
! that omits Dirichlet points.

lapack_matrix%neq = 0
num_nonzero = 0
nodir_renum = -1
do i=1,use_nlev
   do j=phaml_matrix%begin_level(i),phaml_matrix%begin_level(i+1)-1
      if (phaml_matrix%equation_type(j) == DIRICHLET) cycle
      lapack_matrix%neq = lapack_matrix%neq + 1
      nodir_renum(j) = lapack_matrix%neq
      nodir_inv_renum(lapack_matrix%neq) = j
      do k=phaml_matrix%begin_row(j),phaml_matrix%end_row(j)
         if (phaml_matrix%column_index(k) == NO_ENTRY) cycle
         if (phaml_matrix%equation_type(phaml_matrix%column_index(k)) == DIRICHLET) cycle
         if (phaml_matrix%column_index(k) >= phaml_matrix%begin_level(use_nlev+1)) cycle
         num_nonzero = num_nonzero + 1
      end do
   end do
end do

! nullify all array components of lapack_matrix

nullify(lapack_matrix%matrix,lapack_matrix%rhs,lapack_matrix%renum, &
        lapack_matrix%inv_renum,lapack_matrix%ipiv)

! if no equations, never use lapack_matrix so don't continue

if (lapack_matrix%neq == 0) return

! create column pointers for the matrix with Dirichlet points removed

allocate(connec(num_nonzero-lapack_matrix%neq),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in make_lapack_gen_band")
   return
endif

orig_bandwd = 0
counter = 0
do i=1,lapack_matrix%neq
   j = nodir_inv_renum(i)
   rstart(i) = counter + 1
   do k=phaml_matrix%begin_row(j)+1,phaml_matrix%end_row(j)
      if (phaml_matrix%column_index(k) == NO_ENTRY) cycle
      if (phaml_matrix%equation_type(phaml_matrix%column_index(k)) == DIRICHLET) cycle
      if (phaml_matrix%column_index(k) >= phaml_matrix%begin_level(use_nlev+1)) cycle
      counter = counter + 1
      connec(counter) = nodir_renum(phaml_matrix%column_index(k))
      orig_bandwd = max(orig_bandwd,abs(connec(counter)-i))
   end do
   degree(i) = counter - rstart(i) + 1
end do

! Use CALGO 582 (Gibbs, Poole and Stockmeyer) to find a bandwidth reduction
! ordering

workd = 6*lapack_matrix%neq+3
allocate(work(workd),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in make_lapack_gen_band")
   return
endif
work = 0 ! BUG workaround for bug in gpskca; seems to have assumed initial 0's
minband_renum = (/ (i,i=1,phaml_matrix%neq) /)
call gpskca(lapack_matrix%neq,degree,rstart,connec,.false.,workd,minband_renum,&
            work,lapack_matrix%halfbandwidth,profil,jerr,space)
if (jerr /= 0) then
   work(1:lapack_matrix%neq) = (/ (i,i=1,lapack_matrix%neq) /) ! contains minband_inv_renum
   minband_renum(1:lapack_matrix%neq) = (/ (i,i=1,lapack_matrix%neq) /)
   lapack_matrix%halfbandwidth = orig_bandwd
   call warning("Bandwidth reduction reordering routine failed.", &
                "Rediculously large bandwidth may cause allocation failure or long run time.")
endif

! done with column pointers for reduced array

deallocate(connec)

! compose nodir and minband renumberings

! renumbering and inverse renumbering of the equations.  renum(i) is the
! equation number in the band matrix corresponding to original equation i;
! dimension renum(neq(use_nlev)); renum(i) = -1 if i is a Dirichlet point.
! inv_renum(i) is the original equation number corresponding to the band
! matrix equation i; dimension inv_renum(gen_neq).

allocate(lapack_matrix%renum(phaml_matrix%neq),lapack_matrix%inv_renum(lapack_matrix%neq),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in make_lapack_gen_band")
   return
endif

do i=1,phaml_matrix%neq
   if (nodir_renum(i) == -1) then
      lapack_matrix%renum(i) = -1
   else
      lapack_matrix%renum(i) = minband_renum(nodir_renum(i))
   endif
end do

do i=1,lapack_matrix%neq
   lapack_matrix%inv_renum(i) = nodir_inv_renum(work(i))
end do

! done with work space for gpskca

deallocate(work)

! Copy matrix values to general band form with extra space for factorization

allocate(lapack_matrix%matrix(3*lapack_matrix%halfbandwidth+1,lapack_matrix%neq), &
         lapack_matrix%ipiv(lapack_matrix%neq),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in make_lapack_gen_band")
   return
endif
lapack_matrix%matrix = 0.0_my_real

do i=1,use_nlev
   do j=phaml_matrix%begin_level(i),phaml_matrix%begin_level(i+1)-1
      if (phaml_matrix%equation_type(j) == DIRICHLET) cycle
      row = lapack_matrix%renum(j)
      if (row == -1) cycle
      do k=phaml_matrix%begin_row(j),phaml_matrix%end_row(j)
         if (phaml_matrix%column_index(k) == NO_ENTRY) cycle
         if (phaml_matrix%equation_type(phaml_matrix%column_index(k)) == DIRICHLET) cycle
         if (phaml_matrix%column_index(k) >= phaml_matrix%begin_level(use_nlev+1)) cycle
         col = lapack_matrix%renum(phaml_matrix%column_index(k))
         if (col == -1) cycle
         lapack_matrix%matrix(2*lapack_matrix%halfbandwidth+1+row-col,col) = &
            phaml_matrix%matrix_val(k)
      end do
   end do
end do

! Factor the matrix

if (my_real == kind(0.0e0)) then
   call sgbtrf(lapack_matrix%neq,lapack_matrix%neq,lapack_matrix%halfbandwidth,&
               lapack_matrix%halfbandwidth,lapack_matrix%matrix, &
               3*lapack_matrix%halfbandwidth+1,lapack_matrix%ipiv,info)
elseif (my_real == kind(0.0d0)) then
   call dgbtrf(lapack_matrix%neq,lapack_matrix%neq,lapack_matrix%halfbandwidth,&
               lapack_matrix%halfbandwidth,lapack_matrix%matrix, &
               3*lapack_matrix%halfbandwidth+1,lapack_matrix%ipiv,info)
else
   ierr = PHAML_INTERNAL_ERROR
   call fatal("my_real is neither default single nor double precision")
   return
endif
if (info /= 0) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("LAPACK SGBTRF failed.",intlist=(/info/))
   return
endif

! allocate the rhs

allocate(lapack_matrix%rhs(lapack_matrix%neq,1),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in make_lapack_gen_band")
   return
endif

end subroutine make_lapack_gen_band

!          -------------------
subroutine destroy_lapack_band(lapack_matrix)
!          -------------------

!----------------------------------------------------
! This routine gets rid of the band storage of a matrix
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(lapack_band_matrix), intent(inout) :: lapack_matrix
!----------------------------------------------------
! Local variables:

integer :: astat
!----------------------------------------------------
! Begin executable code

if (associated(lapack_matrix%matrix)) deallocate(lapack_matrix%matrix, stat=astat)
if (associated(lapack_matrix%rhs)) deallocate(lapack_matrix%rhs, stat=astat)
if (associated(lapack_matrix%renum)) deallocate(lapack_matrix%renum, stat=astat)
if (associated(lapack_matrix%inv_renum)) deallocate(lapack_matrix%inv_renum, stat=astat)
if (associated(lapack_matrix%ipiv)) deallocate(lapack_matrix%ipiv, stat=astat)

end subroutine destroy_lapack_band

!          ----------
subroutine lapack_spd(use_nlev,phaml_matrix,lapack_matrix)
!          ----------

!----------------------------------------------------
! This routine solves the linear system using the LAPACK routines for
! symmetric positive definite band matricies.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: use_nlev
type(linsys_type), intent(inout) :: phaml_matrix
type(lapack_band_matrix), intent(inout) :: lapack_matrix
!----------------------------------------------------
! Local variables:

integer :: info,i,j,k,lev
!----------------------------------------------------
! Begin executable code

! if the number of equations is 0, there is nothing to solve

if (lapack_matrix%neq == 0) return

! Create the right hand side, by moving to the new numbering and eliminating
! Dirichlet boundary conditions

do lev=1,use_nlev
   do i=phaml_matrix%begin_level(lev),phaml_matrix%begin_level(lev+1)-1
      if (phaml_matrix%equation_type(i) == DIRICHLET) cycle
      j = lapack_matrix%renum(i)
      if (j == -1) cycle
      lapack_matrix%rhs(j,1) = phaml_matrix%rhs(i) + phaml_matrix%r_mine(i) + &
                               phaml_matrix%r_others(i)
      do k=phaml_matrix%begin_row(i),phaml_matrix%end_row(i)
         if (phaml_matrix%column_index(k) == NO_ENTRY) cycle
         if (phaml_matrix%equation_type(phaml_matrix%column_index(k)) == DIRICHLET) then
            lapack_matrix%rhs(j,1)=lapack_matrix%rhs(j,1) - &
               phaml_matrix%matrix_val(k)*phaml_matrix%solution(phaml_matrix%column_index(k))
         endif
      end do
   end do
end do

! Solve the system

if (my_real == kind(0.0e0)) then
   call spbtrs("U",lapack_matrix%neq,lapack_matrix%halfbandwidth,1, &
               lapack_matrix%matrix,lapack_matrix%halfbandwidth+1, &
               lapack_matrix%rhs,lapack_matrix%neq,info)
elseif (my_real == kind(0.0d0)) then
   call dpbtrs("U",lapack_matrix%neq,lapack_matrix%halfbandwidth,1, &
               lapack_matrix%matrix,lapack_matrix%halfbandwidth+1, &
               lapack_matrix%rhs,lapack_matrix%neq,info)
else
   ierr = PHAML_INTERNAL_ERROR
   call fatal("my_real is neither default single nor double precision")
   return
endif
if (info /= 0) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("LAPACK SPBTRS failed.",intlist=(/info/))
   return
endif

! Copy the solution into the solution vector

do i=1,lapack_matrix%neq
   phaml_matrix%solution(lapack_matrix%inv_renum(i)) = lapack_matrix%rhs(i,1)
end do

end subroutine lapack_spd

!          ------------
subroutine lapack_indef(use_nlev,phaml_matrix,lapack_matrix)
!          ------------

!----------------------------------------------------
! This routine solves the linear system using the LAPACK routines for
! general band matricies.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: use_nlev
type(linsys_type), intent(inout) :: phaml_matrix
type(lapack_band_matrix), intent(inout) :: lapack_matrix
!----------------------------------------------------
! Local variables:

integer :: info,i,j,k,lev
!----------------------------------------------------
! Begin executable code

! if the number of equations is 0, there is nothing to solve

if (lapack_matrix%neq == 0) return

! Create the right hand side, by moving to the new numbering and eliminating
! Dirichlet boundary conditions

do lev=1,use_nlev
   do i=phaml_matrix%begin_level(lev),phaml_matrix%begin_level(lev+1)-1
      if (phaml_matrix%equation_type(i) == DIRICHLET) cycle
      j = lapack_matrix%renum(i)
      if (j == -1) cycle
      lapack_matrix%rhs(j,1) = phaml_matrix%rhs(i) + phaml_matrix%r_mine(i) + &
                               phaml_matrix%r_others(i)
      do k=phaml_matrix%begin_row(i),phaml_matrix%end_row(i)
         if (phaml_matrix%column_index(k) == NO_ENTRY) cycle
         if (phaml_matrix%equation_type(phaml_matrix%column_index(k)) == DIRICHLET) then
            lapack_matrix%rhs(j,1)=lapack_matrix%rhs(j,1) - &
               phaml_matrix%matrix_val(k)*phaml_matrix%solution(phaml_matrix%column_index(k))
         endif
      end do
   end do
end do

! Solve the system

if (my_real == kind(0.0e0)) then
   call sgbtrs("N",lapack_matrix%neq,lapack_matrix%halfbandwidth, &
               lapack_matrix%halfbandwidth,1,lapack_matrix%matrix, &
               3*lapack_matrix%halfbandwidth+1,lapack_matrix%ipiv, &
               lapack_matrix%rhs,lapack_matrix%neq,info)
elseif (my_real == kind(0.0d0)) then
   call dgbtrs("N",lapack_matrix%neq,lapack_matrix%halfbandwidth, &
               lapack_matrix%halfbandwidth,1,lapack_matrix%matrix, &
               3*lapack_matrix%halfbandwidth+1,lapack_matrix%ipiv, &
               lapack_matrix%rhs,lapack_matrix%neq,info)
else
   ierr = PHAML_INTERNAL_ERROR
   call fatal("my_real is neither default single nor double precision")
   return
endif
if (info /= 0) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("LAPACK _GBTRS failed.",intlist=(/info/))
   return
endif

! Copy the solution into the solution vector

do i=1,lapack_matrix%neq
   phaml_matrix%solution(lapack_matrix%inv_renum(i)) = lapack_matrix%rhs(i,1)
end do

end subroutine lapack_indef

!          -------------
subroutine lapack_precon(invec,outvec,choice,matrix,procs, &
                         solver_cntl,still_sequential)
!          -------------

!----------------------------------------------------
! This routine applies a preconditioner based on LAPACK, either using
! a solve on the coarse grid or a domain decomposition.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: invec(:)
real(my_real), intent(out) :: outvec(:)
integer, intent(in) :: choice
type(linsys_type), intent(inout) :: matrix
type(proc_info), intent(in) :: procs
type(solver_options), intent(in) :: solver_cntl
logical, intent(in) :: still_sequential
!----------------------------------------------------
! Local variables:

real(my_real) :: holdrhs(matrix%neq), &
                 holdsoln(0:matrix%neq)
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

where (matrix%equation_type == DIRICHLET) &
   matrix%solution(1:) = matrix%rhs

! Call the selected precondtioner

select case(choice)

case (FUDOP_DD_PRECONDITION)
   call fudop_dd_precon(matrix,procs,solver_cntl,still_sequential)
case (COARSE_GRID_PRECONDITION)
   call coarse_precon(matrix,solver_cntl)
   
end select

! Copy solution (which now contains the preconditioner times invec) to outvec

outvec = matrix%solution(1:)

! Restore rhs and solution

matrix%rhs = holdrhs
matrix%solution = holdsoln

end subroutine lapack_precon

!          ---------------
subroutine fudop_dd_precon(matrix,procs,solver_cntl, &
                           still_sequential)
!          ---------------

!----------------------------------------------------
! This routine performs a fudop domain decomposition preconditioner.
! Each processor solves its local problem exactly with LAPACK indefinite
! solver.  It then obtains the solution for unowned equations from the owner,
! makes the unowned points Dirichlet points, and solves again.  Obtaining the
! solution from other processors and solving it iterated some number of times.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(inout) :: matrix
type(proc_info), intent(in) :: procs
type(solver_options), intent(in) :: solver_cntl
logical, intent(in) :: still_sequential
!----------------------------------------------------
! Local variables:

logical :: gen_band_already_existed
integer :: i, ddit, toplev
real(my_real) :: resid(matrix%neq),hold_rhs(matrix%neq),hold_soln(matrix%neq)
!----------------------------------------------------
! Begin executable code

! set the top level of the hierarchy based on the number of equations

if (matrix%neq == matrix%neq_vert) then
   toplev = matrix%nlev
elseif (matrix%neq == matrix%neq_vert+matrix%neq_edge) then
   toplev = matrix%nlev+1
else
   toplev = matrix%nlev+2
endif

! make a LAPACK general band matrix if it doesn't already exist

if (.not. matrix%lapack_gen_band_exists) then
   gen_band_already_existed = .false.
   call make_lapack_gen_band(toplev,matrix,matrix%lapack_mat)
   matrix%lapack_gen_band_exists = .true.
else
   gen_band_already_existed = .true.
endif

! call lapack solver with the local matrix

call lapack_indef(toplev,matrix,matrix%lapack_mat)

! Domain decomposition iterations

do ddit = 1,solver_cntl%dd_iterations

! Compute the global residual and set to 0 at points I don't own

   call matrix_times_vector(matrix%solution(1:),resid,matrix,procs, &
                            still_sequential,1711,1712,1713,1714,1715,1716, &
                            nocomm2=.true.)
   do i=1,matrix%neq
      if (matrix%iown(i)) then
         resid(i) = matrix%rhs(i) - resid(i)
      else
         resid(i) = 0
      endif
   end do

! Set the rhs to the residual and Dirichlet points to 0, to set up an
! error correction problem

   hold_rhs = matrix%rhs
   hold_soln = matrix%solution(1:)
   matrix%rhs = resid
   where (matrix%equation_type == DIRICHLET) matrix%solution(1:) = 0.0_my_real

! convert to hierarchical basis, exchange rhs with other processors, and
! convert back to nodal basis

   if (.not. still_sequential) then
      call basis_change(toplev,TO_HIER,matrix)
      call exchange_fudop_vect(matrix%rhs,procs,matrix,1701,1702,1703)
      call basis_change(toplev,TO_NODAL,matrix)
   endif

! Solve the system

   call lapack_indef(toplev,matrix,matrix%lapack_mat)

! Add correction to solution, and reset rhs and Dirichlet b.c. values

   matrix%solution(1:) = hold_soln + matrix%solution(1:)
   matrix%rhs = hold_rhs

end do

! Get rid of the band matrix if this routine created it

if (.not. gen_band_already_existed) then
   call destroy_lapack_band(matrix%lapack_mat)
   matrix%lapack_gen_band_exists = .false.
endif

end subroutine fudop_dd_precon

!          -------------
subroutine coarse_precon(matrix,solver_cntl)
!          -------------

!----------------------------------------------------
! This routine uses LAPACK as a preconditioner by coarsening to a problem
! small enough to solve with LAPACK and interpolating back to the fine grid.
! Each processor uses LAPACK on the grid it sees with no communication.
! Extraction of the equations it owns and combination with other
! processor's results are the responsiility of the caller.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(inout) :: matrix
type(solver_options), intent(in) :: solver_cntl
!----------------------------------------------------
! Local variables:

integer, save :: maxeq, uselev, hold_neq, lev, eq, i, toplev
logical :: i_made_gen_band
!----------------------------------------------------
! Begin executable code

! set the top level of the hierarchy based on the number of equations

if (matrix%neq == matrix%neq_vert) then
   toplev = matrix%nlev
elseif (matrix%neq == matrix%neq_vert+matrix%neq_edge) then
   toplev = matrix%nlev+1
else
   toplev = matrix%nlev+2
endif

! determine how many levels make up a small enough grid

maxeq = solver_cntl%coarse_size
uselev = 1
do
   if (uselev+1 > toplev) exit
   if (matrix%begin_level(uselev+2) > maxeq) exit
   uselev = uselev + 1
end do

! coarsen to a small grid

do lev=toplev,uselev+1,-1
   call basis_change(lev,TO_HIER,matrix)
end do

! put the Dirichlet boundary values back to their nodal values

where (matrix%equation_type == DIRICHLET) &
   matrix%solution(1:) = matrix%rhs

! call lapack solver

hold_neq = matrix%neq
matrix%neq = matrix%begin_level(uselev+1)-1
i_made_gen_band = .false.
if (.not. matrix%lapack_gen_band_exists) then
   call make_lapack_gen_band(uselev,matrix,matrix%lapack_mat)
   matrix%lapack_gen_band_exists = .true.
   i_made_gen_band = .true.
endif
call lapack_indef(uselev,matrix,matrix%lapack_mat)
matrix%neq = hold_neq

! interpolate to the original grid

do lev=uselev+1,toplev
   call basis_change(lev,TO_NODAL,matrix)
end do

! perform a Gauss-Seidel iteration to get the solution out of the subspace

do eq=1,matrix%neq
   if (matrix%equation_type(eq) == DIRICHLET) cycle
   matrix%solution(eq) = matrix%rhs(eq) + matrix%r_mine(eq) + matrix%r_others(eq)
   do i=matrix%begin_row(eq)+1,matrix%end_row(eq)
      matrix%solution(eq) = matrix%solution(eq)-matrix%matrix_val(i)*matrix%solution(matrix%column_index(i))
   end do
   matrix%solution(eq) = matrix%solution(eq)/matrix%matrix_val(matrix%begin_row(eq))
end do

end subroutine coarse_precon

end module lapack_solve
