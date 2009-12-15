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

module superlu_mod

!----------------------------------------------------
! This module contains Fortran-side wrappers for the SuperLU get/set
! functions, with optional arguments so the user doesn't have to provide
! the full set of components.
!----------------------------------------------------

use superlutype_mod
implicit none

contains

subroutine get_SuperMatrix(A, nrow, ncol)
integer(superlu_ptr) :: A
integer, optional :: nrow, ncol
integer :: l_nrow, l_ncol

call f_get_SuperMatrix(A, l_nrow, l_ncol)

if (present(nrow)) nrow = l_nrow
if (present(ncol)) ncol = l_ncol

end subroutine get_SuperMatrix

subroutine set_SuperMatrix(A, nrow, ncol)
integer(superlu_ptr) :: A
integer, optional :: nrow, ncol
integer :: l_nrow, l_ncol

call f_get_SuperMatrix(A, l_nrow, l_ncol)

if (present(nrow)) l_nrow = nrow
if (present(ncol)) l_ncol = ncol

call f_set_SuperMatrix(A, l_nrow, l_ncol)

end subroutine set_SuperMatrix

subroutine get_superlu_options(opt, Fact, Trans, Equil, RowPerm, &
                               ColPerm, ReplaceTinyPivot, IterRefine, &
                               SolveInitialized, RefineInitialized, PrintStat)
integer(superlu_ptr) :: opt
integer, optional :: Fact, Trans, Equil, RowPerm, ColPerm, &
                     ReplaceTinyPivot, IterRefine, SolveInitialized, &
                     RefineInitialized, PrintStat
integer :: l_Fact, l_Trans, l_Equil, l_RowPerm, l_ColPerm, &
           l_ReplaceTinyPivot, l_IterRefine, l_SolveInitialized, &
           l_RefineInitialized, l_PrintStat

call f_get_superlu_options(opt, l_Fact, l_Trans, l_Equil, l_RowPerm, &
                           l_ColPerm, l_ReplaceTinyPivot, l_IterRefine, &
                           l_SolveInitialized, l_RefineInitialized, l_PrintStat)

if (present(Fact)) Fact = l_Fact
if (present(Trans)) Trans = l_Trans
if (present(Equil)) Equil = l_Equil
if (present(RowPerm)) RowPerm = l_RowPerm
if (present(ColPerm)) ColPerm = l_ColPerm
if (present(ReplaceTinyPivot)) ReplaceTinyPivot = l_ReplaceTinyPivot
if (present(IterRefine)) IterRefine = l_IterRefine
if (present(SolveInitialized)) SolveInitialized = l_SolveInitialized
if (present(RefineInitialized)) RefineInitialized = l_RefineInitialized
if (present(PrintStat)) PrintStat = l_PrintStat

end subroutine get_superlu_options

subroutine set_superlu_options(opt, Fact, Trans, Equil, RowPerm, &
                               ColPerm, ReplaceTinyPivot, IterRefine, &
                               SolveInitialized, RefineInitialized, PrintStat)
integer(superlu_ptr) :: opt
integer, optional :: Fact, Trans, Equil, RowPerm, ColPerm, &
                     ReplaceTinyPivot, IterRefine, SolveInitialized, &
                     RefineInitialized, PrintStat
integer :: l_Fact, l_Trans, l_Equil, l_RowPerm, l_ColPerm, &
           l_ReplaceTinyPivot, l_IterRefine, l_SolveInitialized, &
           l_RefineInitialized, l_PrintStat

call f_get_superlu_options(opt, l_Fact, l_Trans, l_Equil, l_RowPerm, &
                           l_ColPerm, l_ReplaceTinyPivot, l_IterRefine, &
                           l_SolveInitialized, l_RefineInitialized, l_PrintStat)

if (present(Fact)) l_Fact = Fact
if (present(Trans)) l_Trans = Trans
if (present(Equil)) l_Equil = Equil
if (present(RowPerm)) l_RowPerm = RowPerm
if (present(ColPerm)) l_ColPerm = ColPerm
if (present(ReplaceTinyPivot)) l_ReplaceTinyPivot = ReplaceTinyPivot
if (present(IterRefine)) l_IterRefine = IterRefine
if (present(SolveInitialized)) l_SolveInitialized = SolveInitialized
if (present(RefineInitialized)) l_RefineInitialized = RefineInitialized
if (present(PrintStat)) l_PrintStat = PrintStat

call f_set_superlu_options(opt, l_Fact, l_Trans, l_Equil, l_RowPerm, &
                           l_ColPerm, l_ReplaceTinyPivot, l_IterRefine, &
                           l_SolveInitialized, l_RefineInitialized, l_PrintStat)

end subroutine set_superlu_options

end module superlu_mod
