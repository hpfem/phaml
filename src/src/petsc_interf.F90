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

module petsc_interf

!----------------------------------------------------
! This module contains routines that interface to PETSc.  These are
! in a separate file, rather than being part of module linear_system,
! because PETSc requires the use of a C-type preprocessor (by using
! F90 instead of f90 as the suffix, with most compilers).
!
! communication tags in this module are of the form 16xx
!----------------------------------------------------

!----------------------------------------------------
! Other modules used are:

use global
use hash_mod
use hash_eq_mod
use message_passing
use petsc_type_mod
use gridtype_mod
use linsystype_mod
use hbmg
use lapack_solve
use linsys_util

!----------------------------------------------------

implicit none
private
public create_petsc_linear_system, create_petsc_linear_system_mf, &
       change_petsc_rhs, petsc_solve, destroy_petsc_linear_system, &
       petsc_lobpcg_solve_f

!----------------------------------------------------
! The PETSc include files.  Note the use of preprocessor #include instead of
! the Fortran include statement, because the include files contain
! preprocessor directives.

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscpc.h"

!----------------------------------------------------
! The following parameters are defined:

!----------------------------------------------------

!----------------------------------------------------
! The following types are defined:

type petsc_hold_data
   type(petsc_matrix_type), pointer :: petsc_matrix
   type(linsys_type), pointer :: phaml_matrix
   type(linsys_type), pointer :: phaml_full_matrix
   type(grid_type), pointer :: grid
   type(proc_info), pointer :: procs
   type(io_options), pointer :: io_cntl
   type(solver_options), pointer :: solver_cntl
   logical :: still_sequential
end type petsc_hold_data

!----------------------------------------------------

!----------------------------------------------------
! The following variables are defined:

! Pointers to maintain access to variables from callbacks

type(petsc_hold_data) :: petsc_hold

! A PETSc version of the uncondensed matrices, for BLOPEX
! TEMP would be better to have condensed mass and stiffness matrices

type(petsc_matrix_type), target :: petsc_M
type(petsc_matrix_type), target :: petsc_A

! The lapack factored coarse matrix for both stiffness and mass matrices,
! for BLOPEX with SHIFT_SQUARE and multigrid preconditioner

type(lapack_band_matrix) :: coarse_stiffness, coarse_mass

integer :: null_data ! no data to pass

! a big number by which to multiply the identity (Dirichlet) entries of A to
! avoid having a bunch of 1.0 eigenvalues, with BLOPEX

real(my_real) :: big = 1.0e3_my_real
!----------------------------------------------------

contains

!          --------------------------
subroutine create_petsc_linear_system(phaml_matrix,petsc_matrix, &
                                      still_sequential,procs)
!          --------------------------

!----------------------------------------------------
! This routine creates a PETSc linear system (matrix A and vector b) from
! a PHAML linear system.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), target :: phaml_matrix
type(petsc_matrix_type), intent(out), target :: petsc_matrix
logical, intent(in) :: still_sequential
type(proc_info), intent(in) :: procs
!----------------------------------------------------
! Local variables:

integer :: s, e, jerr, astat, i, j, k, pi
PetscScalar :: one
integer, allocatable :: d_nnz(:), o_nnz(:)
!----------------------------------------------------
! Begin executable code

one = 1.0_my_real
petsc_matrix%my_total_eq = phaml_matrix%neq

! Set the equation hash table

petsc_matrix%eq_hash => phaml_matrix%eq_hash

! Determine which equations I own

allocate(petsc_matrix%iown(petsc_matrix%my_total_eq), &
         d_nnz(petsc_matrix%my_total_eq), o_nnz(petsc_matrix%my_total_eq), &
         stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in create_petsc_linear_system")
   return
endif

petsc_matrix%iown = phaml_matrix%iown(1:petsc_matrix%my_total_eq)

! Determine how many nonzeroes in each row.

j = 0
do i=1,petsc_matrix%my_total_eq
   if (petsc_matrix%iown(i)) then
      j = j+1
      d_nnz(j) = 0
      o_nnz(j) = 0
      do k=phaml_matrix%begin_row(i),phaml_matrix%end_row(i)
         if (phaml_matrix%column_index(k) == NO_ENTRY) cycle
         if (petsc_matrix%iown(phaml_matrix%column_index(k))) then
            d_nnz(j) = d_nnz(j) + 1
         else
            o_nnz(j) = o_nnz(j) + 1
         endif
      end do
   endif
end do

! Create the matrix and rhs data structures

petsc_matrix%my_own_eq = count(petsc_matrix%iown)

if (still_sequential) then
   call VecCreateSeq(PETSC_COMM_SELF,petsc_matrix%my_own_eq,petsc_matrix%b,jerr)
else
   call VecCreateMPI(PETSC_COMM_WORLD,petsc_matrix%my_own_eq,PETSC_DECIDE, &
                     petsc_matrix%b,jerr)
endif

! Get the owned range for this processor and total number of equations.

call VecGetOwnershipRange(petsc_matrix%b,petsc_matrix%my_global_low, &
                          petsc_matrix%my_global_hi,jerr)
call VecGetSize(petsc_matrix%b,petsc_matrix%global_eq,jerr)

! Create the matrix data structure

if (still_sequential) then
   call MatCreateSeqAIJ(PETSC_COMM_SELF,petsc_matrix%my_own_eq, &
                        petsc_matrix%my_own_eq,0,d_nnz,petsc_matrix%A,jerr)
else
   call MatCreateMPIAIJ(PETSC_COMM_WORLD,petsc_matrix%my_own_eq, &
                        petsc_matrix%my_own_eq,PETSC_DECIDE,PETSC_DECIDE,0, &
                        d_nnz,0,o_nnz,petsc_matrix%A,jerr)
endif

deallocate(d_nnz,o_nnz,stat=astat)

! Determine the PETSc index for each equation on this processor by
! setting the ones that this processor owns and requesting the
! others from other processors.

allocate(petsc_matrix%petsc_index(petsc_matrix%my_total_eq),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in create_petsc_linear_system")
   return
endif

! Set the ones this processor owns

petsc_matrix%petsc_index = 0
pi = petsc_matrix%my_global_low
do i=1,petsc_matrix%my_total_eq
   if (petsc_matrix%iown(i)) then
      petsc_matrix%petsc_index(i) = pi
      pi = pi + 1
   endif
end do

if (pi /= petsc_matrix%my_global_hi) then
   call warning("count of owned equations not equal to size of petsc range", &
                intlist=(/pi,petsc_matrix%my_global_hi/))
endif

! Get the unowned ones from other processors

call update_shadows(petsc_matrix,phaml_matrix%gid,procs,1601,1602,1603, &
                    idata=petsc_matrix%petsc_index)

! copy the values from rhs to b

call VecSetValues(petsc_matrix%b,petsc_matrix%my_own_eq, &
                  (/(i+petsc_matrix%my_global_low,i=0,petsc_matrix%my_own_eq-1)/), &
                  pack(phaml_matrix%rhs(1:petsc_matrix%my_total_eq),mask=petsc_matrix%iown), &
                  INSERT_VALUES,jerr)
call VecAssemblyBegin(petsc_matrix%b,jerr)

! Copy the values in matrix_val into A

do i=1,petsc_matrix%my_total_eq
   if (.not. petsc_matrix%iown(i)) cycle
   s = phaml_matrix%begin_row(i); e = phaml_matrix%end_row(i)
   if (phaml_matrix%equation_type(i) == DIRICHLET) then
      call MatSetValues(petsc_matrix%A, 1, (/petsc_matrix%petsc_index(i)/), 1, &
                        (/petsc_matrix%petsc_index(i)/), &
                        (/one/), INSERT_VALUES, jerr)
   else
      call MatSetValues(petsc_matrix%A, 1, (/petsc_matrix%petsc_index(i)/), &
                        count(phaml_matrix%column_index(s:e)/=NO_ENTRY), &
                        petsc_matrix%petsc_index(pack(phaml_matrix%column_index(s:e), &
                              phaml_matrix%column_index(s:e)/=NO_ENTRY)), &
                        pack(phaml_matrix%matrix_val(s:e), &
                              phaml_matrix%column_index(s:e)/=NO_ENTRY), &
                        INSERT_VALUES, jerr)
   endif
end do
call MatAssemblyBegin(petsc_matrix%A,MAT_FINAL_ASSEMBLY,jerr)

! finish messages associated with assembling b and A

call VecAssemblyEnd(petsc_matrix%b,jerr)
call MatAssemblyEnd(petsc_matrix%A,MAT_FINAL_ASSEMBLY,jerr)

end subroutine create_petsc_linear_system

!          -----------------------------
subroutine create_petsc_linear_system_mf(phaml_matrix,petsc_matrix, &
                                         still_sequential)
!          -----------------------------

!----------------------------------------------------
! This routine creates a PETSc linear system shell for matrix free solvers.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), target :: phaml_matrix
type(petsc_matrix_type), intent(out) :: petsc_matrix
logical, intent(in) :: still_sequential
!----------------------------------------------------
! Local variables:

integer :: astat,i,jerr
!----------------------------------------------------
! Begin executable code

petsc_matrix%my_total_eq = phaml_matrix%neq

! Set the equation hash table

petsc_matrix%eq_hash => phaml_matrix%eq_hash

! Determine which equations I own

allocate(petsc_matrix%iown(petsc_matrix%my_total_eq), stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in create_petsc_linear_system")
   return
endif

petsc_matrix%iown = phaml_matrix%iown(1:petsc_matrix%my_total_eq)

! Create the matrix and rhs data structures

petsc_matrix%my_own_eq = count(petsc_matrix%iown)

if (still_sequential) then
   call MatCreateShell(PETSC_COMM_SELF,petsc_matrix%my_own_eq, &
                       petsc_matrix%my_own_eq,PETSC_DECIDE, &
                       PETSC_DECIDE,null_data,petsc_matrix%A,jerr)
   call VecCreateSeq(PETSC_COMM_SELF,petsc_matrix%my_own_eq,petsc_matrix%b,jerr)
else
   call MatCreateShell(PETSC_COMM_WORLD,petsc_matrix%my_own_eq, &
                       petsc_matrix%my_own_eq,PETSC_DECIDE, &
                       PETSC_DECIDE,null_data,petsc_matrix%A,jerr)
   call VecCreateMPI(PETSC_COMM_WORLD,petsc_matrix%my_own_eq,PETSC_DECIDE, &
                     petsc_matrix%b,jerr)
endif

! Get the owned range for this processor and total number of equations.

call VecGetOwnershipRange(petsc_matrix%b,petsc_matrix%my_global_low, &
                          petsc_matrix%my_global_hi,jerr)

! copy the values from rhs to b

call VecSetValues(petsc_matrix%b,petsc_matrix%my_own_eq, &
                  (/(i+petsc_matrix%my_global_low,i=0,petsc_matrix%my_own_eq-1)/), &
                  pack(phaml_matrix%rhs(1:petsc_matrix%my_total_eq),mask=petsc_matrix%iown), &
                  INSERT_VALUES,jerr)
call VecAssemblyBegin(petsc_matrix%b,jerr)
call VecAssemblyEnd(petsc_matrix%b,jerr)

! set the matrix multiply callback routine

call MatShellSetOperation(petsc_matrix%A,MATOP_MULT,matmult_mf,jerr)

nullify(petsc_matrix%petsc_index)

end subroutine create_petsc_linear_system_mf

!          ----------------
subroutine change_petsc_rhs(rhs,petsc_matrix)
!          ----------------

!----------------------------------------------------
! This routine copies a PHAML right hand side to a PETSc right hand side
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: rhs(:)
type(petsc_matrix_type), intent(inout) :: petsc_matrix
!----------------------------------------------------
! Local variables:

integer :: jerr, i
!----------------------------------------------------
! Begin executable code

! copy the values from rhs to b

call VecSetValues(petsc_matrix%b,petsc_matrix%my_own_eq, &
                  (/(i+petsc_matrix%my_global_low,i=0,petsc_matrix%my_own_eq-1)/), &
                  pack(rhs,mask=petsc_matrix%iown(:size(rhs))), &
                  INSERT_VALUES,jerr)
call VecAssemblyBegin(petsc_matrix%b,jerr)
call VecAssemblyEnd(petsc_matrix%b,jerr)

end subroutine change_petsc_rhs

!          -----------
subroutine petsc_solve(phaml_matrix,petsc_matrix,solver_cntl,io_cntl, &
                       still_sequential,grid,procs)
!          -----------

!----------------------------------------------------
! This routine solves the linear system by the specified method in PETSc
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(inout), target :: phaml_matrix
type(petsc_matrix_type), intent(inout), target :: petsc_matrix
type(solver_options), intent(in), target :: solver_cntl
type(io_options), intent(in), target :: io_cntl
logical, intent(in) :: still_sequential
type(grid_type), intent(in), target :: grid
type(proc_info), intent(in), target :: procs
!----------------------------------------------------
!----------------------------------------------------
! Local variables:

PC :: pc
KSP :: ksp
integer :: jerr
Vec :: x
real(kind(0.0d0)) :: temp1, temp2, temp3
integer :: temp4, temp5
!----------------------------------------------------
! Begin executable code

! For later access

petsc_hold%petsc_matrix => petsc_matrix
petsc_hold%phaml_matrix => phaml_matrix
nullify(petsc_hold%phaml_full_matrix)
petsc_hold%grid => grid
petsc_hold%procs => procs
petsc_hold%io_cntl => io_cntl
petsc_hold%solver_cntl => solver_cntl
petsc_hold%still_sequential = still_sequential

! set up the solver context

if (still_sequential) then
   call KSPCreate(PETSC_COMM_SELF,ksp,jerr)
else
   call KSPCreate(PETSC_COMM_WORLD,ksp,jerr)
endif
call KSPSetOperators(ksp,petsc_matrix%A,petsc_matrix%A, &
                      DIFFERENT_NONZERO_PATTERN,jerr)
call VecDuplicate(petsc_matrix%b,x,jerr)

! select the solver method

select case (solver_cntl%solver)
case (PETSC_RICHARDSON_SOLVER)
   call KSPSetType(ksp,KSPRICHARDSON,jerr)
case (PETSC_CHEBYCHEV_SOLVER)
   call KSPSetType(ksp,KSPCHEBYCHEV,jerr)
case (PETSC_CG_SOLVER)
   call KSPSetType(ksp,KSPCG,jerr)
case (PETSC_GMRES_SOLVER)
   call KSPSetType(ksp,KSPGMRES,jerr)
case (PETSC_TCQMR_SOLVER)
   call KSPSetType(ksp,KSPTCQMR,jerr)
case (PETSC_BCGS_SOLVER)
   call KSPSetType(ksp,KSPBCGS,jerr)
case (PETSC_CGS_SOLVER)
   call KSPSetType(ksp,KSPCGS,jerr)
case (PETSC_TFQMR_SOLVER)
   call KSPSetType(ksp,KSPTFQMR,jerr)
case (PETSC_CR_SOLVER)
   call KSPSetType(ksp,KSPCR,jerr)
case (PETSC_LSQR_SOLVER)
   call KSPSetType(ksp,KSPLSQR,jerr)
case (PETSC_BICG_SOLVER)
   call KSPSetType(ksp,KSPBICG,jerr)
case default
   ierr = USER_INPUT_ERROR
   call fatal("illegal solver choice in petsc_solve", &
              intlist=(/solver_cntl%solver/))
   return
end select

! select the preconditioning method

call KSPGetPC(ksp,pc,jerr)
select case (solver_cntl%preconditioner)
case (NO_PRECONDITION)
   call PCSetType(pc,PCSHELL,jerr)
   call PCShellSetApply(pc,precon_no,null_data,jerr)
case (MG_PRECONDITION)
   call PCSetType(pc,PCSHELL,jerr)
   call PCShellSetApply(pc,precon_mg,null_data,jerr)
case (FUDOP_DD_PRECONDITION)
   call PCSetType(pc,PCSHELL,jerr)
   call PCShellSetApply(pc,precon_fudop_dd,null_data,jerr)
case (COARSE_GRID_PRECONDITION)
   call PCSetType(pc,PCSHELL,jerr)
   call PCShellSetApply(pc,precon_coarse,null_data,jerr)
case (TEST_PRECONDITION)
   call PCSetType(pc,PCSHELL,jerr)
   call PCShellSetApply(pc,precon_test,null_data,jerr)
case (PETSC_JACOBI_PRECONDITION)
   call PCSetType(pc,PCJACOBI,jerr)
case (PETSC_BJACOBI_PRECONDITION)
   call PCSetType(pc,PCBJACOBI,jerr)
case (PETSC_SOR_PRECONDITION)
   call PCSetType(pc,PCSOR,jerr)
case (PETSC_EISENSTAT_PRECONDITION)
   call PCSetType(pc,PCEISENSTAT,jerr)
case (PETSC_ICC_PRECONDITION)
   call PCSetType(pc,PCICC,jerr)
case (PETSC_ILU_PRECONDITION)
   call PCSetType(pc,PCILU,jerr)
case (PETSC_ASM_PRECONDITION)
   call PCSetType(pc,PCASM,jerr)
case default
   ierr = USER_INPUT_ERROR
   call fatal("illegal preconditioner choice in petsc_solve", &
              intlist=(/solver_cntl%preconditioner/))
   return
end select

! set PETSc options

if (solver_cntl%petsc_cntl%petsc_richardson_damping_factor /= huge(0.0d0)) then
   call KSPRichardsonSetScale(ksp,solver_cntl%petsc_cntl%petsc_richardson_damping_factor,jerr)
endif

if (solver_cntl%petsc_cntl%petsc_chebychev_emin /= huge(0.0d0) .or. &
    solver_cntl%petsc_cntl%petsc_chebychev_emin /= huge(0.0d0)) then
   temp1 = solver_cntl%petsc_cntl%petsc_chebychev_emin
   if (temp1 == huge(0.0d0)) temp1 = PETSC_DEFAULT_DOUBLE_PRECISION
   temp2 = solver_cntl%petsc_cntl%petsc_chebychev_emax
   if (temp2 == huge(0.0d0)) temp2 = PETSC_DEFAULT_DOUBLE_PRECISION
   call KSPChebychevSetEigenvalues(ksp,temp1,temp2,jerr)
endif

if (solver_cntl%petsc_cntl%petsc_gmres_max_steps /= huge(0)) then
   call KSPGMRESSetRestart(ksp,solver_cntl%petsc_cntl%petsc_gmres_max_steps,jerr)
endif

if (solver_cntl%petsc_cntl%petsc_rtol /= huge(0.0d0) .or. &
    solver_cntl%petsc_cntl%petsc_atol /= huge(0.0d0) .or. &
    solver_cntl%petsc_cntl%petsc_dtol /= huge(0.0d0) .or. &
    solver_cntl%petsc_cntl%petsc_maxits /= huge(0)) then
   temp1 = solver_cntl%petsc_cntl%petsc_rtol
   if (temp1 == huge(0.0d0)) temp1 = PETSC_DEFAULT_DOUBLE_PRECISION
   temp2 = solver_cntl%petsc_cntl%petsc_atol
   if (temp2 == huge(0.0d0)) temp2 = PETSC_DEFAULT_DOUBLE_PRECISION
   temp3 = solver_cntl%petsc_cntl%petsc_dtol
   if (temp3 == huge(0.0d0)) temp3 = PETSC_DEFAULT_DOUBLE_PRECISION
   temp4 = solver_cntl%petsc_cntl%petsc_maxits
   if (temp4 == huge(0)) temp4 = PETSC_DEFAULT_INTEGER
   call KSPSetTolerances(ksp,temp1,temp2,temp3,temp4,jerr)
endif

if (solver_cntl%petsc_cntl%petsc_ilu_levels /= huge(0)) then

! For PETSc versions before 2.3.1
!   call PCILUSetLevels(pc,solver_cntl%petsc_cntl%petsc_ilu_levels,jerr)

! For PETSc version 2.3.1 and later
   call PCFactorSetLevels(pc,solver_cntl%petsc_cntl%petsc_ilu_levels,jerr)
endif

if (solver_cntl%petsc_cntl%petsc_icc_levels /= huge(0)) then

! For PETSc versions before 2.3.1
!   call PCICCSetLevels(pc,solver_cntl%petsc_cntl%petsc_icc_levels,jerr)

! For PETSc version 2.3.1 and later
   call PCFactorSetLevels(pc,solver_cntl%petsc_cntl%petsc_icc_levels,jerr)
endif

if (solver_cntl%petsc_cntl%petsc_ilu_dt /= huge(0.0d0) .or. &
    solver_cntl%petsc_cntl%petsc_ilu_dtcol /= huge(0.0d0) .or. &
    solver_cntl%petsc_cntl%petsc_ilu_maxrowcount /= huge(0)) then
   temp1 = solver_cntl%petsc_cntl%petsc_ilu_dt
   if (temp1 == huge(0.0d0)) temp1 = PETSC_DEFAULT_DOUBLE_PRECISION
   temp2 = solver_cntl%petsc_cntl%petsc_ilu_dtcol
   if (temp2 == huge(0.0d0)) temp2 = PETSC_DEFAULT_DOUBLE_PRECISION
   temp4 = solver_cntl%petsc_cntl%petsc_ilu_maxrowcount
   if (temp4 == huge(0)) temp4 = PETSC_DEFAULT_INTEGER

! For PETSc versions before 2.3.1
!   call PCILUSetUseDropTolerance(pc,temp1,temp2,temp4,jerr)

! For PETSc version 2.3.1 and later
   call PCFactorSetUseDropTolerance(pc,temp1,temp2,temp4,jerr)
endif

if (solver_cntl%petsc_cntl%petsc_sor_omega /= huge(0.0d0)) then
   call PCSORSetOmega(pc,solver_cntl%petsc_cntl%petsc_sor_omega,jerr)
endif

if (solver_cntl%petsc_cntl%petsc_sor_its /= huge(0) .or. &
    solver_cntl%petsc_cntl%petsc_sor_lits /= huge(0)) then
   temp4 = solver_cntl%petsc_cntl%petsc_sor_its
   if (temp4 == huge(0)) temp4 = PETSC_DEFAULT_INTEGER
   temp5 = solver_cntl%petsc_cntl%petsc_sor_lits
   if (temp5 == huge(0)) temp5 = PETSC_DEFAULT_INTEGER
   call PCSORSetIterations(pc,temp4,temp5,jerr)
endif

if (solver_cntl%petsc_cntl%petsc_eisenstat_nodiagscaling) then
   call PCEisenstatNoDiagonalScaling()
endif

if (solver_cntl%petsc_cntl%petsc_eisenstat_omega /= huge(0.0d0)) then
   call PCEisenstatSetOmega(pc,solver_cntl%petsc_cntl%petsc_eisenstat_omega,jerr)
endif

if (solver_cntl%petsc_cntl%petsc_asm_overlap /= huge(0)) then
   call PCASMSetOverlap(pc,solver_cntl%petsc_cntl%petsc_asm_overlap,jerr)
endif

! Print L2 norm of residual after each iteration

if (io_cntl%print_error_when == FREQUENTLY .or. &
    io_cntl%print_error_when == TOO_MUCH) then

! For PETSc versions before 2.3.3
!   call KSPSetMonitor(ksp,KSPTrueMonitor,PETSC_NULL,PETSC_NULL,jerr)

! For PETSc version 2.3.3 and later
  call KSPMonitorSet(ksp,KSPMonitorTrueResidualNorm,PETSC_NULL,PETSC_NULL,jerr)

endif

! solve the system

call KSPSolve(ksp,petsc_matrix%b,x,jerr)
if (jerr /= 0) call warning("PETSc solver returned error code", &
                            intlist=(/jerr/))

! extract the solution

call petscvec_to_phaml(x,phaml_matrix%solution(1:),petsc_matrix)

! free memory

call KSPDestroy(ksp,jerr)
call VecDestroy(x,jerr)

end subroutine petsc_solve

!          ---------------------------
subroutine destroy_petsc_linear_system(petsc_matrix)
!          ---------------------------

!----------------------------------------------------
! This routine frees the space for a PETSc matrix and right hand side
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(petsc_matrix_type), intent(inout) :: petsc_matrix
!----------------------------------------------------
! Local variables:

integer :: jerr, astat
!----------------------------------------------------
! Begin executable code

call MatDestroy(petsc_matrix%A,jerr)
call VecDestroy(petsc_matrix%b,jerr)
if (associated(petsc_matrix%iown)) deallocate(petsc_matrix%iown,stat=astat)
if (associated(petsc_matrix%petsc_index)) deallocate(petsc_matrix%petsc_index,stat=astat)

end subroutine destroy_petsc_linear_system

!          -----------------
subroutine petscvec_to_phaml(petscvec,phamlvec,petsc_matrix)
!          -----------------

!----------------------------------------------------
! This routine copies a PETSc Vec to a PHAML distributed vector
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real) :: phamlvec(:)
Vec :: petscvec
type(petsc_matrix_type), intent(in) :: petsc_matrix
!----------------------------------------------------
! Local variables:

PetscScalar :: x_array(1)
PetscOffset :: i_x
integer :: jerr, pi, i
!----------------------------------------------------
! Begin executable code

! Extract the local values from petscvec

call VecGetArray(petscvec,x_array,i_x,jerr)

! Copy the ones I own to the PHAML vector

pi = 0
do i=1,petsc_matrix%my_total_eq
   if (petsc_matrix%iown(i)) then
      pi = pi + 1
      phamlvec(i) = x_array(i_x + pi)
   endif
end do

! Request the unowned ones from other processors

if (.not. petsc_hold%still_sequential) then
   call update_shadows(petsc_matrix,petsc_hold%phaml_matrix%gid, &
                       petsc_hold%procs, 1604, 1605, 1606, rdata=phamlvec)
endif

! free memory

call VecRestoreArray(petscvec,x_array,i_x,jerr)

end subroutine petscvec_to_phaml

!          -----------------
subroutine phamlvec_to_petsc(phamlvec,petscvec,petsc_matrix)
!          -----------------

!----------------------------------------------------
! This routine copies a PHAML distributed vector to a PETSC vec
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real) :: phamlvec(:)
Vec :: petscvec
type(petsc_matrix_type), intent(in) :: petsc_matrix
!----------------------------------------------------
! Local variables:

integer :: jerr, i
!----------------------------------------------------
! Begin executable code

call VecSetValues(petscvec,petsc_matrix%my_own_eq, &
                  (/(i+petsc_matrix%my_global_low,i=0,petsc_matrix%my_own_eq-1)/), &
                  pack(phamlvec,mask=petsc_matrix%iown), &
                  INSERT_VALUES,jerr)
call VecAssemblyBegin(petscvec,jerr)
call VecAssemblyEnd(petscvec,jerr)

end subroutine phamlvec_to_petsc

!          --------------
subroutine update_shadows(petsc_matrix,gid,procs,tag1,tag2,tag3,rdata,idata)
!          --------------

!----------------------------------------------------
! This routine updates the values of rdata and/or idata at shadow points.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(petsc_matrix_type), intent(in) :: petsc_matrix
type(hash_key_eq), intent(in) :: gid(:)
type(proc_info), intent(in) :: procs
integer, intent(in) :: tag1, tag2, tag3
real(my_real), intent(inout), optional :: rdata(:)
integer, intent(inout), optional :: idata(:)
!----------------------------------------------------
! Local variables:

integer :: counter, i, p, lid, astat, nproc, my_processor, isub, rsub, &
           oisub, orsub, limit, KEY_SIZE_EQ
! newcomm
integer :: nisend
integer, allocatable :: isend(:), nisendv(:), nirecv(:), nrsend(:), nrrecv(:)
integer, pointer :: irecv(:)
real(my_real), allocatable :: rsend(:)
real(my_real), pointer :: rrecv(:)
!----------------------------------------------------
! Begin executable code

KEY_SIZE_EQ = KEY_SIZE+1

! allocate space for received messages

nproc = num_proc(procs)
my_processor = my_proc(procs)
allocate(nisendv(nproc),nirecv(nproc),nrsend(nproc),nrrecv(nproc),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in update_shadows")
   return
endif

nisend = (petsc_matrix%my_total_eq-petsc_matrix%my_own_eq)*KEY_SIZE_EQ
allocate(isend(nisend),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in update_shadows")
   return
endif

! Make a list of the ones this processor doesn't own

counter=1
do i=1,petsc_matrix%my_total_eq
   if (petsc_matrix%iown(i)) cycle
   call hash_pack_key(gid(i),isend,counter)
   counter = counter + KEY_SIZE_EQ
end do

! Send the request

call phaml_alltoall(procs,isend,nisend,irecv,nirecv,tag1)

! Reply with ones I own

! Count the number of responses

counter = 0
isub = 1
do p=1,nproc
   if (p == my_processor) then
      isub = isub + nirecv(p)
      cycle
   endif
   do i=1,nirecv(p)/KEY_SIZE_EQ
      lid = hash_decode_key(hash_unpack_key(irecv,isub,.true.), &
                            petsc_matrix%eq_hash)
      if (lid /= HASH_NOT_FOUND) then
         if (petsc_matrix%iown(lid)) then
            counter = counter + 1
         endif
      endif
      isub = isub + KEY_SIZE_EQ
   end do
end do

! allocate memory

deallocate(isend)
if (present(idata)) then
   allocate(isend(counter*(KEY_SIZE_EQ+1)), stat=astat)
else
   allocate(isend(counter*KEY_SIZE_EQ), stat=astat)
endif
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in update_shadows")
   return
endif
if (present(rdata)) then
   allocate(rsend(counter), stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in update_shadows")
      return
   endif
endif

! Make the arrays with the responses

counter = 1
isub = 1
oisub = isub
rsub = 1
orsub = rsub
do p=1,nproc
   if (p == my_processor) then
      nisendv(p) = 0
      nrsend(p) = 0
      counter = counter + nirecv(p)
      cycle
   endif
   do i=1,nirecv(p)/KEY_SIZE_EQ
      lid = hash_decode_key(hash_unpack_key(irecv,counter,.true.), &
                            petsc_matrix%eq_hash)
      if (lid /= HASH_NOT_FOUND) then
         if (petsc_matrix%iown(lid)) then
            isend(isub:isub+KEY_SIZE_EQ-1) = irecv(counter:counter+KEY_SIZE_EQ-1)
            isub = isub + KEY_SIZE_EQ
            if (present(rdata)) then
               rsend(rsub) = rdata(lid)
               rsub = rsub + 1
            endif
            if (present(idata)) then
               isend(isub) = idata(lid)
               isub = isub + 1
            endif
         endif
      endif
      counter = counter + KEY_SIZE_EQ
   end do
   nisendv(p) = isub - oisub
   oisub = isub
   nrsend(p) = rsub - orsub
   orsub = rsub
end do

if (associated(irecv)) deallocate(irecv,stat=astat)

! Send the replies

call phaml_alltoall(procs,isend,nisendv,irecv,nirecv,tag2)
if (present(rdata)) then
   call phaml_alltoall(procs,rsend,nrsend,rrecv,nrrecv,tag3)
else
   nullify(rrecv)
endif

deallocate(isend,rsend,stat=astat)

! Set the shadow values from the replies

isub = 1
rsub = 1
do p=1,nproc
   if (p == my_processor) cycle
   if (present(idata)) then
      limit = nirecv(p)/(KEY_SIZE_EQ+1)
   else
      limit = nirecv(p)/KEY_SIZE_EQ
   endif
   do i=1,limit
      lid = hash_decode_key(hash_unpack_key(irecv,isub,.true.), &
                            petsc_matrix%eq_hash)
      isub = isub+KEY_SIZE_EQ
      if (lid == HASH_NOT_FOUND) then
         call warning("received reply for an equation I don't have in update_shadows")
      else
         if (present(rdata)) then
            rdata(lid) = rrecv(rsub)
            rsub = rsub + 1
         endif
         if (present(idata)) then
            idata(lid) = irecv(isub)
            isub = isub + 1
         endif
      endif
   end do
end do

if (associated(irecv)) deallocate(irecv,stat=astat)
if (associated(rrecv)) deallocate(rrecv,stat=astat)
deallocate(nisendv,nirecv,nrsend,nrrecv,stat=astat)

end subroutine update_shadows

!          ---------
subroutine precon_no(nodata,avec,bvec,jerr)
!          ---------

!----------------------------------------------------
! This routine applies no preconditioning.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer :: nodata
Vec :: avec, bvec
integer :: jerr
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

call VecCopy(avec,bvec,jerr)

end subroutine precon_no

!          ---------
subroutine precon_mg(nodata,avec,bvec,jerr)
!          ---------

!----------------------------------------------------
! This routine provides the callback for V-cycle multigrid preconditioning.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer :: nodata
Vec :: avec, bvec
integer :: jerr
!----------------------------------------------------
! Local variables:

real(my_real), allocatable :: invec(:),outvec(:) 
integer :: isize

!----------------------------------------------------
! Begin executable code

! TEMP 04/04/08 the condensed number of equations should be sufficient (ifort
!      with array bounds checking ran fine), but gfortran gives an error during
!      free unless I use the full number of equations

isize = petsc_hold%phaml_matrix%neq_vert + &
        petsc_hold%phaml_matrix%neq_edge + &
        petsc_hold%phaml_matrix%neq_face

allocate(invec(isize),outvec(isize))

! convert PETSc avec to a normal vector, apply preconditioner, and convert back

call petscvec_to_phaml(avec,invec,petsc_hold%petsc_matrix)
call mg_precon(invec,outvec,petsc_hold%phaml_matrix, &
               petsc_hold%grid,petsc_hold%procs,petsc_hold%io_cntl, &
               petsc_hold%solver_cntl,petsc_hold%still_sequential)
call phamlvec_to_petsc(outvec,bvec,petsc_hold%petsc_matrix)

deallocate(invec,outvec)

end subroutine precon_mg

!          ---------------
subroutine precon_fudop_dd(nodata,avec,bvec,jerr)
!          ---------------

!----------------------------------------------------
! This routine provides the callback for FuDoP domain decomposition
! preconditioning.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer :: nodata
Vec :: avec, bvec
integer :: jerr
!----------------------------------------------------
! Local variables:

real(my_real) :: invec(petsc_hold%petsc_matrix%my_total_eq), &
                outvec(petsc_hold%petsc_matrix%my_total_eq)

!----------------------------------------------------
! Begin executable code

call petscvec_to_phaml(avec,invec,petsc_hold%petsc_matrix)
call lapack_precon(invec,outvec,FUDOP_DD_PRECONDITION,petsc_hold%phaml_matrix, &
                   petsc_hold%procs,petsc_hold%solver_cntl, &
                   petsc_hold%still_sequential)
call phamlvec_to_petsc(outvec,bvec,petsc_hold%petsc_matrix)

end subroutine precon_fudop_dd

!          -------------
subroutine precon_coarse(nodata,avec,bvec,jerr)
!          -------------

!----------------------------------------------------
! This routine provides the callback for coarse grid preconditioning.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer :: nodata
Vec :: avec, bvec
integer :: jerr
!----------------------------------------------------
! Local variables:

real(my_real) :: invec(petsc_hold%petsc_matrix%my_total_eq), &
                outvec(petsc_hold%petsc_matrix%my_total_eq)

!----------------------------------------------------
! Begin executable code

call petscvec_to_phaml(avec,invec,petsc_hold%petsc_matrix)
call lapack_precon(invec,outvec,COARSE_GRID_PRECONDITION, &
                   petsc_hold%phaml_matrix, &
                   petsc_hold%procs,petsc_hold%solver_cntl, &
                   petsc_hold%still_sequential)
call phamlvec_to_petsc(outvec,bvec,petsc_hold%petsc_matrix)

end subroutine precon_coarse

!          -----------
subroutine precon_test(nodata,avec,bvec,jerr)
!          -----------

!----------------------------------------------------
! This routine provides the callback for testing a new preconditioner.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer :: nodata
Vec :: avec, bvec
integer :: jerr
!----------------------------------------------------
! Local variables:

real(my_real) :: invec(petsc_hold%petsc_matrix%my_total_eq), &
                outvec(petsc_hold%petsc_matrix%my_total_eq)

!----------------------------------------------------
! Begin executable code

call petscvec_to_phaml(avec,invec,petsc_hold%petsc_matrix)
call test_precon(invec,outvec,petsc_hold%phaml_matrix,petsc_hold%grid, &
                 petsc_hold%procs,petsc_hold%solver_cntl, &
                 petsc_hold%still_sequential)
call phamlvec_to_petsc(outvec,bvec,petsc_hold%petsc_matrix)

end subroutine precon_test

!          -------------
subroutine matmult_mf(matrix,avec,bvec,jerr)
!          -------------

!----------------------------------------------------
! This routine provides the callback for matrix free matmult
! as an interface between PETSc and an external routine in linsys.f90.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

Mat :: matrix
Vec :: avec, bvec
integer :: jerr
!----------------------------------------------------
! Local variables:

real(my_real) :: invec(petsc_hold%petsc_matrix%my_total_eq), &
                outvec(petsc_hold%petsc_matrix%my_total_eq)

!----------------------------------------------------
! Begin executable code

call petscvec_to_phaml(avec,invec,petsc_hold%petsc_matrix)
call matrix_times_vector(invec,outvec,petsc_hold%phaml_matrix, &
                         petsc_hold%procs,petsc_hold%still_sequential, &
                         1621,1622,1623,1624,1625,1626)
call phamlvec_to_petsc(outvec,bvec,petsc_hold%petsc_matrix)

end subroutine matmult_mf

!--------------------------------------------------------------------
! BLOPEX interface
!--------------------------------------------------------------------

!          --------------------
subroutine petsc_lobpcg_solve_f(phaml_matrix,full_matrix,petsc_matrix, &
                                solver_cntl,io_cntl,still_sequential,grid,procs)
!          --------------------

!----------------------------------------------------
! This routine calls the C routine that calls the petsc version of lobpcg
! to solve the eigenproblem.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(inout), target :: phaml_matrix, full_matrix
type(petsc_matrix_type), intent(inout), target :: petsc_matrix
type(solver_options), intent(in), target :: solver_cntl
type(io_options), intent(in), target :: io_cntl
logical, intent(in) :: still_sequential
type(grid_type), intent(in), target :: grid
type(proc_info), intent(in), target :: procs
!----------------------------------------------------
! Local variables:

integer :: jerr
!----------------------------------------------------
! Begin executable code

! petsc_lobpcg_solve_c assumes eigenvalue is a double

if (my_real /= kind(0.0d0)) then
   ierr = UNCLASSIFIED_ERROR
   call fatal("lobpcg assumes double precision; change my_real in global.f90",&
              procs=procs)
   stop
endif

! For later access

petsc_hold%petsc_matrix => petsc_matrix
petsc_hold%phaml_matrix => phaml_matrix
petsc_hold%phaml_full_matrix => full_matrix
petsc_hold%grid => grid
petsc_hold%procs => procs
petsc_hold%io_cntl => io_cntl
petsc_hold%solver_cntl => solver_cntl
petsc_hold%still_sequential = still_sequential

! Make a PETSC matrix for the full mass matrix

full_matrix%matrix_val => full_matrix%mass
call create_petsc_linear_system(full_matrix,petsc_M, &
                                still_sequential,procs)
full_matrix%matrix_val => full_matrix%stiffness

! If doing shift and square, also need the full stiffness matrix

if (solver_cntl%transformation == SHIFT_SQUARE) then
   call create_petsc_linear_system(full_matrix,petsc_A,still_sequential,procs)
endif

! If doing shift and square and the preconditioner is multigrid, need
! coarse grid matrices for both the stiffness and mass matrix

if (solver_cntl%transformation == SHIFT_SQUARE .and. &
    solver_cntl%preconditioner == MG_PRECONDITION) then
   full_matrix%matrix_val => full_matrix%mass
   call make_lapack_symm_band(1,full_matrix,coarse_mass)
   full_matrix%matrix_val => full_matrix%stiffness
   call make_lapack_symm_band(1,full_matrix,coarse_stiffness)
   full_matrix%coarse_matrix = coarse_stiffness
   full_matrix%coarse_band_exists = .true.
endif

! call the C routine

call petsc_lobpcg_solve_c(petsc_M%b,petsc_M%my_own_eq, &
                          solver_cntl%num_eval, &
                          solver_cntl%blopex_cntl%maxit, &
                          solver_cntl%blopex_cntl%atol, &
                          solver_cntl%blopex_cntl%rtol, &
                          grid%eigenvalue, petsc_lobpcg_opA, &
                          petsc_lobpcg_opB, petsc_lobpcg_opT, &
                          petsc_lobpcg_return_evec, petsc_lobpcg_initial_guess,&
                          jerr)

if (jerr /= 0) then
   call phaml_send(petsc_hold%procs,MASTER,(/jerr,my_proc(petsc_hold%procs)/), &
                   2,(/0.0_my_real/),0,2201)
endif

! destroy objects created above

if (solver_cntl%transformation == SHIFT_SQUARE .and. &
    solver_cntl%preconditioner == MG_PRECONDITION) then
   call destroy_lapack_band(coarse_stiffness)
   call destroy_lapack_band(coarse_mass)
endif

if (solver_cntl%transformation == SHIFT_SQUARE) then
   call destroy_petsc_linear_system(petsc_A)
endif
call destroy_petsc_linear_system(petsc_M)

end subroutine petsc_lobpcg_solve_f

!          --------------------------
subroutine petsc_lobpcg_initial_guess(evec)
!          --------------------------

!----------------------------------------------------
! This routine sets the initial guess for the petsc lobpcg C routine
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

Vec :: evec
!----------------------------------------------------
! Local variables:

integer, save :: count=0
integer :: p, ni, nr
integer, pointer :: irecv(:)
real(my_real), pointer :: rrecv(:)
!----------------------------------------------------
! Begin executable code

if (count == 0) then
   if (maxval(abs(petsc_hold%phaml_full_matrix%solution(1:))) < 10*epsilon(0.0_my_real)) then
      if (.not. petsc_hold%still_sequential) then
         ierr = USER_INPUT_ERROR
         call fatal("First call to BLOPEX requires still_sequential")
         stop
      endif
      if (my_proc(petsc_hold%procs) == 1) then
         petsc_hold%phaml_full_matrix%solution(0) = 0
         call random_number(petsc_hold%phaml_full_matrix%solution(1:))
         where (petsc_hold%phaml_full_matrix%equation_type==DIRICHLET) &
            petsc_hold%phaml_full_matrix%solution(1:) = 0.0_my_real
         do p=2,num_proc(petsc_hold%procs)
            call phaml_send(petsc_hold%procs,p,(/0/),0, &
                            petsc_hold%phaml_full_matrix%solution(1:), &
                           size(petsc_hold%phaml_full_matrix%solution(1:)),1611)
         end do
      else
         call phaml_recv(petsc_hold%procs,p,irecv,ni,rrecv,nr,1611)
         petsc_hold%phaml_full_matrix%solution(1:) = rrecv
         deallocate(rrecv)
      endif
      call phamlvec_to_petsc(petsc_hold%phaml_full_matrix%solution(1:),evec, &
                             petsc_M)
   else
      call phamlvec_to_petsc(petsc_hold%phaml_full_matrix%solution(1:),evec, &
                             petsc_M)
   endif
else
   if (maxval(abs(petsc_hold%phaml_full_matrix%evecs(:,count))) < 10*epsilon(0.0_my_real)) then
      if (my_proc(petsc_hold%procs) == 1) then
         call random_number(petsc_hold%phaml_full_matrix%evecs(:,count))
         where (petsc_hold%phaml_full_matrix%equation_type==DIRICHLET) &
            petsc_hold%phaml_full_matrix%evecs(:,count) = 0.0_my_real
         do p=2,num_proc(petsc_hold%procs)
            call phaml_send(petsc_hold%procs,p,(/0/),0, &
                            petsc_hold%phaml_full_matrix%evecs(:,count), &
                            size(petsc_hold%phaml_full_matrix%evecs(:,count)), &
                            1612)
         end do
      else
         call phaml_recv(petsc_hold%procs,p,irecv,ni,rrecv,nr,1612)
         petsc_hold%phaml_full_matrix%evecs(:,count) = rrecv
         deallocate(rrecv)
      endif
      call phamlvec_to_petsc(petsc_hold%phaml_full_matrix%evecs(:,count),evec, &
                             petsc_M)
   else
      call phamlvec_to_petsc(petsc_hold%phaml_full_matrix%evecs(:,count),evec, &
                             petsc_M)
   endif
endif

count = count+1
if (count == petsc_hold%solver_cntl%num_eval) count = 0

end subroutine petsc_lobpcg_initial_guess

!          ------------------------
subroutine petsc_lobpcg_return_evec(evec)
!          ------------------------

!----------------------------------------------------
! This routine receives the eigenvectors from the petsc lobpcg C routine
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

Vec :: evec
!----------------------------------------------------
! Local variables:

integer, save :: count=0
!----------------------------------------------------
! Begin executable code

if (count == 0) then
   call petscvec_to_phaml(evec,petsc_hold%phaml_full_matrix%solution(1:), &
                          petsc_M)
else
   call petscvec_to_phaml(evec,petsc_hold%phaml_full_matrix%evecs(:,count), &
                          petsc_M)
endif

count = count+1
if (count == petsc_hold%solver_cntl%num_eval) count = 0

end subroutine petsc_lobpcg_return_evec

!          ----------------
subroutine petsc_lobpcg_opA(matrix,avec,bvec,jerr)
!          ----------------

!----------------------------------------------------
! This routine provides the callback for multiplication by the shift-invert
! matrix M (A - lambda0 M)^(-1) M for lobpcg, or by A if there is not shift.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

Mat :: matrix
Vec :: avec, bvec
integer :: jerr
!----------------------------------------------------
! Local variables:

real(my_real) :: invec(petsc_M%my_total_eq), &
                outvec(petsc_M%my_total_eq), &
                 cond_invec(petsc_M%my_total_eq), &
                cond_outvec(petsc_M%my_total_eq), &
                tempvec(petsc_M%my_total_eq)
type(petsc_hold_data) :: hold_petsc_hold
type(linsys_type), pointer :: full_matrix
real(my_real), pointer :: hold_full_matrix_val(:)
integer :: objtype,brank,srank,objlid,loc_neq,eq,info,i,j,full_neq,cond_neq
real(my_real), pointer :: block(:,:)
integer, pointer :: ipiv(:)
!----------------------------------------------------
! Begin executable code

! convert input PETSc vector to a normal vector

call petscvec_to_phaml(avec,invec,petsc_M)

! If lambda0 = -infinity, the operator is A

if (petsc_hold%solver_cntl%lambda0 == -huge(0.0_my_real)) then

   call matrix_times_vector(invec,outvec,petsc_hold%phaml_full_matrix, &
                            petsc_hold%procs,petsc_hold%still_sequential, &
                            1631,1632,1633,1634,1635,1636,nodirch=.true.)
! multiply Dirichlet rows by a big number so the corresponding eigenvalues
! are at the end of the spectrum instead of 1.0
   where (petsc_hold%phaml_full_matrix%equation_type == DIRICHLET) &
      outvec = big*outvec

! If lambda0 is not -infinity and the transformation is SHIFT_INVERT,
! the operator is M (A - lambda0 M)^-1 M

elseif (petsc_hold%solver_cntl%transformation == SHIFT_INVERT) then

! hold the current petsc_hold

   hold_petsc_hold = petsc_hold

! for convenience

   full_matrix => petsc_hold%phaml_full_matrix
   cond_neq = petsc_hold%phaml_matrix%neq
   full_neq = full_matrix%neq
   hold_full_matrix_val => full_matrix%matrix_val

! multiply the input vector by the mass matrix

   full_matrix%matrix_val => full_matrix%mass
   call matrix_times_vector(invec,outvec,full_matrix, &
                            petsc_hold%procs,petsc_hold%still_sequential, &
                            1641,1642,1643,1644,1645,1646,nodirch=.true.)

! perform static condensation on the result, putting the result in cond_invec
! b1 - A12^T A2^-1 b2 where 2 is face bases

   cond_invec = outvec

   full_matrix%matrix_val => full_matrix%shifted
   eq = cond_neq+1
   do
      if (eq > full_neq) exit
      call eq_to_grid(full_matrix,full_matrix%gid(eq),objtype,brank, &
                      srank,objlid,petsc_hold%grid)
      block => full_matrix%elem_block(objlid)%matrix
      ipiv => full_matrix%elem_block(objlid)%ipiv
      loc_neq = full_matrix%elem_block(objlid)%neq
      tempvec(1:loc_neq) = outvec(eq:eq+loc_neq-1)
      if (my_real == kind(0.0)) then
        call sgetrs("N",loc_neq,1,block,loc_neq,ipiv,tempvec,loc_neq,info)
      else
        call dgetrs("N",loc_neq,1,block,loc_neq,ipiv,tempvec,loc_neq,info)
      endif
      if (info /= 0) then
         ierr = PHAML_INTERNAL_ERROR
         call fatal("lapack sgetrs solution failed", &
                    intlist=(/info/),procs=petsc_hold%procs)
         stop
      endif
      do i=1,loc_neq
         do j=full_matrix%begin_row(eq+i-1),full_matrix%end_row(eq+i-1)
            if (full_matrix%column_index(j) == NO_ENTRY) cycle
            if (full_matrix%column_index(j) >= eq .and. &
                full_matrix%column_index(j) < eq+loc_neq) cycle
            cond_invec(full_matrix%column_index(j)) = &
                       cond_invec(full_matrix%column_index(j)) - &
                       full_matrix%matrix_val(j)*tempvec(i)
         end do
      end do
      eq = eq + loc_neq
   end do

! TEMP I don't know why the ARPACK version of this needs the sum, but
!      this one needs to not have it.  Also, I do it in opT.
!      if (cond_neq /= full_neq .and. .not. petsc_hold%still_sequential .and. &
!          num_proc(petsc_hold%procs) > 1) then
   if (.false. .and. &
       cond_neq /= full_neq .and. .not. petsc_hold%still_sequential .and. &
       num_proc(petsc_hold%procs) > 1) then
      call sum_fudop_vect(cond_invec(1:cond_neq),petsc_hold%procs, &
                          petsc_hold%phaml_matrix,1614,1615,1616)
   endif

! copy the result to PETSc matrix rhs

   call change_petsc_rhs(cond_invec(1:cond_neq),petsc_hold%petsc_matrix)

! solve the linear system

   call petsc_solve(petsc_hold%phaml_matrix,petsc_hold%petsc_matrix, &
                    petsc_hold%solver_cntl,petsc_hold%io_cntl, &
                   petsc_hold%still_sequential,petsc_hold%grid,petsc_hold%procs)
   petsc_hold = hold_petsc_hold

! set cond_outvec to (x1 b2)^T where x1 is the solution just computed and
! b2 is from the result of multiplying the input vector by M

   cond_outvec(1:cond_neq) = petsc_hold%phaml_matrix%solution(1:cond_neq)
   cond_outvec(cond_neq+1:full_neq) = outvec(cond_neq+1:full_neq)

! solve for the face bases; x2 = A2^-1 (b2 - A12 x1)

   eq = cond_neq+1
   do
      if (eq > full_neq) exit
      call eq_to_grid(full_matrix,full_matrix%gid(eq),objtype,brank, &
                      srank,objlid,petsc_hold%grid)
      block => full_matrix%elem_block(objlid)%matrix
      ipiv => full_matrix%elem_block(objlid)%ipiv
      loc_neq = full_matrix%elem_block(objlid)%neq
      do i=1,loc_neq
         do j=full_matrix%begin_row(eq+i-1),full_matrix%end_row(eq+i-1)
            if (full_matrix%column_index(j) == NO_ENTRY) cycle
            if (full_matrix%column_index(j) >= eq .and. &
                full_matrix%column_index(j) < eq+loc_neq) cycle
            cond_outvec(eq+i-1) = cond_outvec(eq+i-1)- &
              full_matrix%matrix_val(j)*cond_outvec(full_matrix%column_index(j))
         end do
      end do
      if (my_real == kind(0.0)) then
         call sgetrs("N",loc_neq,1,block,loc_neq,ipiv, &
                     cond_outvec(eq:eq+loc_neq-1),loc_neq,info)
      else
         call dgetrs("N",loc_neq,1,block,loc_neq,ipiv, &
                     cond_outvec(eq:eq+loc_neq-1),loc_neq,info)
      endif
      if (info /= 0) then
         ierr = PHAML_INTERNAL_ERROR
         call fatal("lapack spotrs solution failed for solving for face bases",&
                    intlist=(/info/),procs=petsc_hold%procs)
         stop
      endif
      eq = eq + loc_neq
   end do

! second multiplication by the mass matrix

   full_matrix%matrix_val => full_matrix%mass
   call matrix_times_vector(cond_outvec,outvec,full_matrix, &
                            petsc_hold%procs,petsc_hold%still_sequential, &
                            1651,1652,1653,1654,1655,1656,nodirch=.true.)

! restore the phaml_full_matrix%matrix_val to its original pointer

   full_matrix%matrix_val => hold_full_matrix_val

! If lambda0 is not -infinity and the transformation is SHIFT_SQUARE,
! the operator is (A - lambda0 M) M^-1 (A - lambda0 M)

else

! multiply by A - lambda0 M.  Since the shifted matrix has been condensed,
! we have to use the stiffness and mass matrices for this

   hold_full_matrix_val => petsc_hold%phaml_full_matrix%matrix_val
   petsc_hold%phaml_full_matrix%matrix_val => petsc_hold%phaml_full_matrix%stiffness
   call matrix_times_vector(invec,outvec,petsc_hold%phaml_full_matrix, &
                            petsc_hold%procs,petsc_hold%still_sequential, &
                            1631,1632,1633,1634,1635,1636,nodirch=.true.)
   petsc_hold%phaml_full_matrix%matrix_val => petsc_hold%phaml_full_matrix%mass
   call matrix_times_vector(invec,tempvec,petsc_hold%phaml_full_matrix, &
                            petsc_hold%procs,petsc_hold%still_sequential, &
                            1641,1642,1643,1644,1645,1646,nodirch=.true.)
   outvec = outvec - petsc_hold%solver_cntl%lambda0 * tempvec

! multiply by M^-1

   hold_petsc_hold = petsc_hold
   call change_petsc_rhs(outvec,petsc_M)
   if (petsc_hold%solver_cntl%transformation == SHIFT_SQUARE .and. &
       petsc_hold%solver_cntl%preconditioner == MG_PRECONDITION) then
      petsc_hold%phaml_full_matrix%coarse_matrix = coarse_mass
   endif
   call petsc_solve(petsc_hold%phaml_full_matrix,petsc_M,petsc_hold%solver_cntl, &
                       petsc_hold%io_cntl,petsc_hold%still_sequential, &
                       petsc_hold%grid,petsc_hold%procs)
   petsc_hold = hold_petsc_hold
   if (petsc_hold%solver_cntl%transformation == SHIFT_SQUARE .and. &
       petsc_hold%solver_cntl%preconditioner == MG_PRECONDITION) then
      petsc_hold%phaml_full_matrix%coarse_matrix = coarse_stiffness
   endif
   invec = petsc_hold%phaml_full_matrix%solution(1:)

! multiply by A - lambda0 M again

   petsc_hold%phaml_full_matrix%matrix_val => petsc_hold%phaml_full_matrix%stiffness
   call matrix_times_vector(invec,outvec,petsc_hold%phaml_full_matrix, &
                            petsc_hold%procs,petsc_hold%still_sequential, &
                            1651,1652,1653,1654,1655,1656,nodirch=.true.)
   petsc_hold%phaml_full_matrix%matrix_val => petsc_hold%phaml_full_matrix%mass
   call matrix_times_vector(invec,tempvec,petsc_hold%phaml_full_matrix, &
                            petsc_hold%procs,petsc_hold%still_sequential, &
                            1661,1662,1663,1664,1665,1666,nodirch=.true.)
   outvec = outvec - petsc_hold%solver_cntl%lambda0 * tempvec

   petsc_hold%phaml_full_matrix%matrix_val => hold_full_matrix_val

endif ! lambda0 == -inf or which transformation

! copy outvec to PETSc vector

call phamlvec_to_petsc(outvec,bvec,petsc_M)

end subroutine petsc_lobpcg_opA

!          ----------------
subroutine petsc_lobpcg_opB(matrix,avec,bvec,jerr)
!          ----------------

!----------------------------------------------------
! This routine provides the callback for multiplication by the mass
! matrix (operator B) for lobpcg.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

Mat :: matrix
Vec :: avec, bvec
integer :: jerr
!----------------------------------------------------
! Local variables:

real(my_real) :: invec(petsc_M%my_total_eq), &
                outvec(petsc_M%my_total_eq)
real(my_real), pointer :: hold_matval(:)

!----------------------------------------------------
! Begin executable code

call petscvec_to_phaml(avec,invec,petsc_M)
hold_matval => petsc_hold%phaml_full_matrix%matrix_val
petsc_hold%phaml_full_matrix%matrix_val => petsc_hold%phaml_full_matrix%mass
call matrix_times_vector(invec,outvec,petsc_hold%phaml_full_matrix, &
                         petsc_hold%procs,petsc_hold%still_sequential, &
                         1661,1662,1663,1664,1665,1666,nodirch=.true.)
petsc_hold%phaml_full_matrix%matrix_val => hold_matval
call phamlvec_to_petsc(outvec,bvec,petsc_M)

end subroutine petsc_lobpcg_opB

!          ----------------
subroutine petsc_lobpcg_opT(matrix,avec,bvec,jerr)
!          ----------------

!----------------------------------------------------
! This routine provides the callback for multiplication by the approximate
! inverse of A, i.e. preconditioning (operator T), for lobpcg, or by the
! identity if there is a shift lambda0.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

Mat :: matrix
Vec :: avec, bvec
integer :: jerr
!----------------------------------------------------
! Local variables:

real(my_real) :: invec(petsc_M%my_total_eq), &
                outvec(petsc_M%my_total_eq), &
                 cond_invec(petsc_M%my_total_eq), &
                tempvec(petsc_M%my_total_eq)
type(petsc_hold_data) :: hold_petsc_hold
type(linsys_type), pointer :: full_matrix
integer :: objtype,brank,srank,objlid,loc_neq,eq,info,i,j,full_neq,cond_neq
real(my_real), pointer :: block(:,:)
integer, pointer :: ipiv(:)
real(my_real), pointer :: hold_matval(:)

!----------------------------------------------------
! Begin executable code

! convert input PETSc vector to a normal vector

call petscvec_to_phaml(avec,invec,petsc_M)

! hold the current petsc_hold

hold_petsc_hold = petsc_hold

! for convenience

full_matrix => petsc_hold%phaml_full_matrix
cond_neq = petsc_hold%phaml_matrix%neq
full_neq = full_matrix%neq

! If lambda0 = -infinity, precondition with A^-1

if (petsc_hold%solver_cntl%lambda0 == -huge(0.0_my_real)) then

! inverse of the adjustment to operatorA to avoid many 1.0 eigenvalues

   where (petsc_hold%phaml_full_matrix%equation_type == DIRICHLET) invec = invec/big

! perform static condensation on invec, putting the result in cond_invec
! b1 - A12^T A2^-1 b2 where 2 is face bases

   cond_invec = invec

   eq = cond_neq+1
   do
      if (eq > full_neq) exit
      call eq_to_grid(full_matrix,full_matrix%gid(eq),objtype,brank, &
                      srank,objlid,petsc_hold%grid)
      block => full_matrix%elem_block(objlid)%matrix
      ipiv => full_matrix%elem_block(objlid)%ipiv
      loc_neq = full_matrix%elem_block(objlid)%neq
      tempvec(1:loc_neq) = invec(eq:eq+loc_neq-1)
      if (my_real == kind(0.0)) then
        call sgetrs("N",loc_neq,1,block,loc_neq,ipiv,tempvec,loc_neq,info)
      else
        call dgetrs("N",loc_neq,1,block,loc_neq,ipiv,tempvec,loc_neq,info)
      endif
      if (info /= 0) then
         ierr = PHAML_INTERNAL_ERROR
         call fatal("lapack sgetrs solution failed", &
                    intlist=(/info/),procs=petsc_hold%procs)
         stop
      endif
      do i=1,loc_neq
         do j=full_matrix%begin_row(eq+i-1),full_matrix%end_row(eq+i-1)
            if (full_matrix%column_index(j) == NO_ENTRY) cycle
            if (full_matrix%column_index(j) >= eq .and. &
                full_matrix%column_index(j) < eq+loc_neq) cycle
            cond_invec(full_matrix%column_index(j)) = &
                       cond_invec(full_matrix%column_index(j)) - &
                       full_matrix%matrix_val(j)*tempvec(i)
         end do
      end do
      eq = eq + loc_neq
   end do

   if (.not. petsc_hold%still_sequential .and. &
       num_proc(petsc_hold%procs) > 1) then
      call sum_fudop_vect(cond_invec(1:cond_neq),petsc_hold%procs, &
                          petsc_hold%phaml_matrix,1617,1618,1619)
   endif

! copy condensed invec to PETSc matrix rhs

   call change_petsc_rhs(cond_invec(1:cond_neq),petsc_hold%petsc_matrix)

! solve the linear system

   call petsc_solve(petsc_hold%phaml_matrix,petsc_hold%petsc_matrix, &
                    petsc_hold%solver_cntl,petsc_hold%io_cntl, &
                    petsc_hold%still_sequential,petsc_hold%grid,petsc_hold%procs)
   petsc_hold = hold_petsc_hold

! set outvec to (x1 b2)^T where x1 is the solution just computed and
! b2 is from invec

   outvec(1:cond_neq) = petsc_hold%phaml_matrix%solution(1:cond_neq)
   outvec(cond_neq+1:full_neq) = invec(cond_neq+1:full_neq)

! solve for the face bases; x2 = A2^-1 (b2 - A12 x1)

   eq = cond_neq+1
   do
      if (eq > full_neq) exit
      call eq_to_grid(full_matrix,full_matrix%gid(eq),objtype,brank, &
                      srank,objlid,petsc_hold%grid)
      block => full_matrix%elem_block(objlid)%matrix
      ipiv => full_matrix%elem_block(objlid)%ipiv
      loc_neq = full_matrix%elem_block(objlid)%neq
      do i=1,loc_neq
         do j=full_matrix%begin_row(eq+i-1),full_matrix%end_row(eq+i-1)
            if (full_matrix%column_index(j) == NO_ENTRY) cycle
            if (full_matrix%column_index(j) >= eq .and. &
                full_matrix%column_index(j) < eq+loc_neq) cycle
            outvec(eq+i-1) = outvec(eq+i-1)- &
               full_matrix%matrix_val(j)*outvec(full_matrix%column_index(j))
         end do
      end do
      if (my_real == kind(0.0)) then
         call sgetrs("N",loc_neq,1,block,loc_neq,ipiv, &
                     outvec(eq:eq+loc_neq-1),loc_neq,info)
      else
         call dgetrs("N",loc_neq,1,block,loc_neq,ipiv, &
                     outvec(eq:eq+loc_neq-1),loc_neq,info)
      endif
      if (info /= 0) then
         ierr = PHAML_INTERNAL_ERROR
         call fatal("lapack spotrs solution failed for solving for face bases", &
                    intlist=(/info/),procs=petsc_hold%procs)
         stop
      endif
      eq = eq + loc_neq
   end do

! If lambda0 is not -infinity and the transformation is SHIFT_INVERT,
! precondition with M^-1

elseif (petsc_hold%solver_cntl%transformation == SHIFT_INVERT) then


   call change_petsc_rhs(invec,petsc_M)
   full_matrix%matrix_val => full_matrix%mass
   call petsc_solve(full_matrix,petsc_M,petsc_hold%solver_cntl, &
                    petsc_hold%io_cntl,petsc_hold%still_sequential, &
                    petsc_hold%grid,petsc_hold%procs)
   petsc_hold = hold_petsc_hold
   outvec = full_matrix%solution(1:)
   full_matrix%matrix_val => full_matrix%stiffness

! If lambda0 is not -infinity and the transformation is SHIFT_SQUARE,
! precondition with A^-1 M A^-1

else

   hold_petsc_hold = petsc_hold
   call change_petsc_rhs(invec,petsc_A)
   hold_matval => petsc_hold%phaml_full_matrix%matrix_val
   petsc_hold%phaml_full_matrix%matrix_val => petsc_hold%phaml_full_matrix%stiffness
   call petsc_solve(petsc_hold%phaml_full_matrix,petsc_A,petsc_hold%solver_cntl, &
                       petsc_hold%io_cntl,petsc_hold%still_sequential, &
                       petsc_hold%grid,petsc_hold%procs)
   petsc_hold = hold_petsc_hold
   invec = petsc_hold%phaml_full_matrix%solution(1:)
   petsc_hold%phaml_full_matrix%matrix_val => petsc_hold%phaml_full_matrix%mass
   call matrix_times_vector(invec,outvec,petsc_hold%phaml_full_matrix, &
                         petsc_hold%procs,petsc_hold%still_sequential, &
                         1661,1662,1663,1664,1665,1666,nodirch=.true.)
   invec = outvec
   call change_petsc_rhs(invec,petsc_A)
   petsc_hold%phaml_full_matrix%matrix_val => petsc_hold%phaml_full_matrix%stiffness
   call petsc_solve(petsc_hold%phaml_full_matrix,petsc_A,petsc_hold%solver_cntl, &
                       petsc_hold%io_cntl,petsc_hold%still_sequential, &
                       petsc_hold%grid,petsc_hold%procs)
   petsc_hold = hold_petsc_hold
   petsc_hold%phaml_full_matrix%matrix_val => hold_matval
   invec = petsc_hold%phaml_full_matrix%solution(1:)

   outvec = invec

endif ! lambda0 == -inf

! copy outvec to PETSc vector

call phamlvec_to_petsc(outvec,bvec,petsc_M)

end subroutine petsc_lobpcg_opT

end module petsc_interf
