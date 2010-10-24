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

module linsystype_mod

!----------------------------------------------------
! This module contains defined types for lineary_system and related modules.
!----------------------------------------------------

!----------------------------------------------------
! Other modules used are:

use global
use hash_mod
use hash_eq_mod
use message_passing
use petsc_type_mod
use mumps_struc_mod
use hypretype_mod
use superlutype_mod
!----------------------------------------------------

implicit none
private
public NO_ENTRY, lapack_band_matrix, linsys_type, solver_options, &
       arpack_options, arpack_dummy, blopex_options, level1_ancestors, &
       encased_matrix, lid_arrays

!----------------------------------------------------
! The following parameters are defined:

! column_index for a nonexistent entry.  This must be a valid subscript
! for solution, and solution(NO_ENTRY) should be 0.0
integer, parameter :: NO_ENTRY = 0

!----------------------------------------------------

!----------------------------------------------------
! The following types are defined:

type encased_matrix
   real(my_real), pointer :: matrix(:,:)
   integer, pointer :: ipiv(:)
   integer :: neq
end type

! Storage in LAPACK band form, either symmetric or nonsymmetric
! (symmetric does not use ipiv).

type lapack_band_matrix
   real(my_real), pointer :: matrix(:,:), rhs(:,:)
   integer, pointer :: renum(:), inv_renum(:), ipiv(:)
   integer :: neq, halfbandwidth
end type

type lid_arrays
   integer, pointer :: lid(:)
end type lid_arrays

type linsys_type

! The matrix.  This is a compressed row storage with the rows in order.
! begin_row(i) is the index in matrix_val to start row i; the row ends at
! end_row(i).  There may be a gap between end_row(i) and begin_row(i+1) after
! restriction to a coarse grid.  end_row_linear gives the end of the entries
! for linear-basis equations (these columns come before high-order-basis
! equations), end_row_edge gives the end of the entries for edge equations,
! and end_row_face gives the end of the entries for face equations, i.e., the
! end of the row.  When the matrix is in it's usual p-hierarchical form,
! end_row points to end_row_face.  It points to end_row_edge when working
! with the condensed matrix from static condensation, and to end_row_linear
! when converting to h-hierarchical bases.
! column_index(k) is the column index for matrix_val(k).
! There may be nonexistant entries in matrix_val; for these the column_index
! is set to NO_ENTRY, which indexes a 0.0 entry in the solution vector so it
! can safely be used in computations (the matrix_val will be defined, but will
! not necessarily be 0.0).  
! For multicomponent solutions (systems of equations), rows and columns can be
! considered to be in blocks of size equal to the number of components
! (syssize).  In the first row of a block, the first syssize entries are the
! diagonal block in increasing order of local id, and the remaining entries
! are in increasing order of local id.  Subsequent rows of the block are the
! same except that in the first block of columns of the i'th row, the first and
! i'th entries are swapped so that the diagonal is the first entry. The matrix
! is stored as nonsymmetric, even though it is symmetric, so the storage is
! available for the matrix after h-hierarchical basis changes.  Dirichlet
! boundary rows contain the matrix values as if they were a Neuman point for
! the sake of basis changes; a separate flag indicates they should be handled
! special.  For a Dirichlet equation, the solution value is set at the time
! the equation is set up, and never changed (except possibly basis changes).

   integer, pointer :: column_index(:), begin_row(:), end_row(:), &
                       end_row_linear(:), end_row_edge(:), end_row_face(:)

! stiffness is the usual stiffness matrix (contains energy inner products
! of the basis functions).  It is the matrix in Ax=b.  It may contain the
! statically condensed values of the matrix.
! mass is the mass matrix (contains L2 inner products of the basis functions).
! It is the matrix M in the generalized eigenproblem Ax=lambda*Mx.
! shifted is the shifted matrix A-lambda0*M where lambda0 is close to the
! eigenvalues to be computed.  It actually contains the statically condensed
! values of the matrix.
! condensed is the Schur complement matrix formed by static condensation of
! the face bases in stiffness or shifted.  For elliptic boundary value
! problems it will point to stiffness unless the stiffness matrix is needed
! to compute residuals.  For elliptic eigenvalue problems it will point to
! shifted.
! matrix_val is used to access matrix values most of the time.  Usually it
! will point to condensed, but at times it will point to other matrices, e.g.
! it will point to mass when matrix_times_vector is called to multiply a
! vector by the mass matrix.
! All matrices have the same nonzero structure so they can share
! the rest of the data structure.

   real(my_real), pointer :: matrix_val(:)
   real(my_real), pointer :: stiffness(:), mass(:), shifted(:), condensed(:)

! The equations are numbered such that all the equations from linear basis
! functions (associated with vertices) come first, and those from the same
! level are contiguous.  begin_level(i) is the first row for level i; the last
! row is begin_level(i+1)-1.  begin_level(nlev+1) is the first row for edge
! bases; begin_level(nlev+2) is the first row for face bases.
! begin_level(nlev+3) is neq+1.  The edge (face) rows have those associated
! with an edge (face) contiguous and in the order that the p-hierarchical basis
! returns them.

   integer, pointer :: begin_level(:)

! The right hand side of the system, with and without static condensation

   real(my_real), pointer :: rhs(:), rhs_nocond(:), rhs_cond(:)

! The solution vector, like the other module variables, only exists while
! subroutine solve is running.  At the end of solve it is copied to the
! solution components of the grid data structure.

   real(my_real), pointer :: solution(:)

! Additional solutions for eigenvalue problems when multiple eigenvalues
! are computed

   real(my_real), pointer :: evecs(:,:)

! equation_type indicates the kind of vertex or edge (interior, Neuman,
! Dirichlet, etc.) this equation came from, to indicate the need for special
! handling.

   integer, pointer :: equation_type(:)

! For communication, each equation needs a gid.

   type(hash_key_eq), pointer :: gid(:)

! The hash table that gives local equation numbers that correspond to gids.

   type(hash_table_eq) :: eq_hash

! which equations are owned by this processor

   logical(small_logical), pointer :: iown(:)

! TEMP090126
! lists used in communication

   integer, pointer :: nn_comm_end_of_send(:,:)
   type(lid_arrays), pointer :: nn_comm_send_lid(:), nn_comm_recv_lid(:)
   logical(small_logical), pointer :: nn_comm_remote_neigh(:)
   integer, pointer :: fudop_comm_proxy(:), fudop_comm_nsend(:), &
                       fudop_comm_mess_size(:,:), fudop_comm_nproc_recv(:), &
                       fudop_comm_nchunk(:,:)
   type(lid_arrays), pointer :: fudop_comm_send_lid(:), fudop_comm_recv_lid(:)
   logical(small_logical), pointer :: fudop_comm_remote_neigh(:)
   integer, pointer :: resid_comm_nsend(:), resid_comm_mess_size(:,:), &
                       resid_comm_nproc_recv(:), resid_comm_nchunk(:,:)
   type(lid_arrays), pointer :: resid_comm_send_lid(:), resid_comm_recv_lid(:)
   logical(small_logical), pointer :: resid_comm_remote_neigh(:)

! The basis change between nodal and h-hierarchical requires that the
! solution be changed, including at the Dirichlet boundary points.  But
! when converted back to nodal it must have the original value.  hold_dirich
! keeps those values.

   real(my_real), pointer :: hold_dirich(:)

! The right hand side is modified during the multigrid cycle via the basis
! changes, but does not pick up the "residuals" from subtracting
! A sub br sup T times x sub r (b and r meaning black and red, using the
! notation of ETNA, 6 (1998) pp. 224-233).  These residuals are contained
! in r_mine and r_others.  See the ETNA paper.

   real(my_real), pointer :: r_mine(:), r_others(:)
   logical(small_logical), pointer :: need_r_others(:)

! Blocks corresponding to all the equations associated with an
! element interior are stored in PLU factored form

   type(encased_matrix), pointer :: elem_block(:)

! The number of equations, refinement levels and coupled PDEs, and maximum
! polynomial degree, number of vertex, edge and face equations.
! When treating the matrix as statically condensed, neq = neq_vert+neq_edge,
! otherwise neq_vert+neq_edge+neq_face

   integer :: neq, nlev, system_size, maxdeg, neq_vert, neq_edge, neq_face

! The minimum value of PDE coefficient r(x,y), used as lambda0 when
! the smallest eigenvalues are desired.

   real(my_real) :: rmin

! The matrices used for conversion between nodal and h-hierarchical bases.
! Multiplication by S converts the coefficient vector of a h-hierarchical
! basis to a nodal basis.  Create matricies for S and S^-1 for interior
! and boundary vertices.

   real(my_real), pointer :: s_int(:,:), s_int_inv(:,:), &
                             s_bnd(:,:), s_bnd_inv(:,:)

! The matrix for the coarse grid is also stored in LAPACK band form
! for direct solution.  coarse_neq does not include DIRICHLET boundary
! equations, which are excluded from the matrix.  coarse_halfbandwidth
! does not include the main diagonal.

   type(lapack_band_matrix) :: coarse_matrix
   logical :: coarse_band_exists

! Storage for the whole matrix as a band matrix, used for the LAPACK
! solvers and preconditioners, either as a symmetric band matrix or a
! general band matrix

   type(lapack_band_matrix) :: lapack_mat
   logical :: lapack_gen_band_exists, lapack_symm_band_exists

! Storage of the matrix in PETSc format

   type(petsc_matrix_type) :: petsc_matrix
   logical :: petsc_matrix_exists

! Storage of the matrix in MUMPS format

   type(mumps_matrix_type) :: mumps_matrix
   logical :: mumps_matrix_exists

! Storage of the matrix in hypre format

   type(hypre_matrix_type) :: hypre_matrix
   logical :: hypre_matrix_exists

! Storage of the matrix in SuperLU format

   type(superlu_matrix_type) :: superlu_matrix
   logical :: superlu_matrix_exists

end type linsys_type

type arpack_options
   real(my_real) :: tol
   integer :: ncv, maxit
end type arpack_options

! Can be used for solver_options constructors
type(arpack_options), parameter :: arpack_dummy=arpack_options(0.0_my_real,0,0)

type blopex_options
   real(my_real) :: atol, rtol
   integer :: maxit
end type blopex_options

type solver_options
   logical :: ignore_quad_err, petsc_matrix_free
   integer :: solver, preconditioner, ncycle, prerelax, postrelax, &
              prerelax_ho, postrelax_ho, dd_iterations, eq_type, &
              num_eval, system_size, coarse_size, coarse_method, &
              inc_quad_order, scale_evec, krylov_iter, krylov_restart, &
              eigensolver, lambda0_side, transformation, mg_comm
   real(my_real) :: lambda0, mg_tol, krylov_tol
   type(hypre_options) :: hypre_cntl
   type(arpack_options) :: arpack_cntl
   type(blopex_options) :: blopex_cntl
   type(petsc_options) :: petsc_cntl
end type solver_options

!----------------------------------------------------
! The following variables are defined:

integer :: level1_ancestors

end module linsystype_mod
