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

module eigen

!----------------------------------------------------

!----------------------------------------------------
! This module contains routines for solve eigenproblems
!
! communication tags in this module are of the form 2xxx
!----------------------------------------------------

!----------------------------------------------------
! Other modules used:

use linsystype_mod
use global
use message_passing
use gridtype_mod
use petsc_interf
use mumps_interf
use hypre_interf
use superlu_interf
use linsys_util
use sort_mod
use krylov
use lapack_solve
use hbmg
use make_linsys
!----------------------------------------------------

implicit none
private
public eigen_arpack, eigen_blopex

contains

!          ------------
subroutine eigen_arpack(grid,procs,linsys,io_cntl,solver_cntl,still_sequential)
!          ------------

!----------------------------------------------------
! This routine solves the eigenvalue problem Ax = lambda*Mx where A is 
! stiffness and M is mass.  It uses ARPACK to solve for the largest 
! eigenvalue(s) and eigenvector(s) of
!
!                -1
! (A - lambda0*M)  M x = mu*x
!
! which gives lambda = 1/mu + lambda0 to give the eigenvalues closest to
! lambda0 (but they are computed by the Raleigh quotient rather than this
! formula).  If lambda0 == -huge(0), then rmin (the minimum value of the
! coefficient of u in the PDE, which should be smaller than the smallest
! eigenvalue) is used to get the smallest eigenvalues.
!
! This is the traditional "shift and invert" spectral transformation.  In
! ARPACK notation
!                      -1
!  OP = (A - lambda0*M)  M  and B = M
!
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(proc_info), intent(in) :: procs
type(linsys_type), intent(inout), target :: linsys
type(io_options), intent(in) :: io_cntl
type(solver_options), intent(in) :: solver_cntl
logical, intent(in) :: still_sequential

!----------------------------------------------------
! Local variables:

! choose between using the symmetric solvers and nonsymmetric solvers

logical, parameter :: symm = .false.

real(my_real), allocatable :: arnoldi_vec(:,:), arp_workd(:), arp_resid(:), &
                 arp_workl(:), dr(:), di(:), evec(:,:), workev(:)
integer :: nev, ncv, lworkl, kev
integer :: arp_param(11), arp_pntr(14)
integer, allocatable :: iperm(:)
real(my_real) :: tol, maxsoln, minsoln
real(my_real) :: l2resid, l2rhs, solnMnorm, soln2norm, sums(2)
integer :: ido, info, iter, nproc, comm, neqown, allocstat, neq_total, &
           ndirch, neq, neq_full, objtype,brank,srank,objlid,loc_neq, ii, jj, &
           nresid_comp, lev
character(len=1) :: bmat, howmany
character(len=2) :: which
logical :: rvec, select(solver_cntl%arpack_cntl%ncv)
integer :: i, j, ss, eq
real(my_real) :: tempvec(linsys%neq_vert+linsys%neq_edge+linsys%neq_face), &
                 Ax(linsys%neq_vert+linsys%neq_edge+linsys%neq_face), &
                 Mx(linsys%neq_vert+linsys%neq_edge+linsys%neq_face)
real(my_real), pointer :: tmpptr(:), block(:,:)
integer, pointer :: ipiv(:)
real(my_real) :: my_resid(linsys%neq)
real(my_real) :: numer, denom, lambda0
logical :: i_made_lapack, i_made_petsc, i_made_mumps, i_made_superlu, &
           i_made_hypre
type(linsys_type) :: full_matrix
integer :: ni, nr, proc
integer, pointer :: irecv(:)
real(my_real), pointer :: rrecv(:)

!----------------------------------------------------
! Begin executable code

! in case we abort early

grid%arpack_iter = 0
grid%arpack_nconv = 0
grid%arpack_numop = 0
grid%arpack_numopb = 0
grid%arpack_numreo = 0
grid%arpack_info = 0

! master hangs around just to see if it needs to print a warning

if (my_proc(procs) == MASTER) then
   do
      call phaml_recv(procs,proc,irecv,ni,rrecv,nr,2201)
      select case (irecv(1))
      case (1) ! all done
         deallocate(irecv,stat=allocstat)
         exit
      case (2) ! ARPACK eigenvalue does not agree with Raleigh quotient
         call warning("Eigenvalue returned by ARPACK is not close to Raleigh quotient; may be wrong", &
                      "Eigenvalue number, eigenvalue and Raleigh quotient are:", &
                      (/irecv(2)/),rrecv)
         deallocate(irecv,rrecv,stat=allocstat)
      end select
   end do
   return
endif

! use linsys for solving the (statically condensed) shifted linear systems,
! and full_matrix for operations that need the whole matrix.  full_matrix
! points to entries in linsys with the following pointers/values different:
!  linsys%end_row => end_row_edge
!  linsys%matrix_val => condensed => shifted
!  linsys%rhs => rhs_cond
!  linsys%neq = neq_vert+neq_edge

full_matrix = linsys

full_matrix%end_row => linsys%end_row_face
full_matrix%matrix_val => linsys%stiffness
full_matrix%rhs => linsys%rhs_nocond
full_matrix%neq = linsys%neq_vert + linsys%neq_edge + linsys%neq_face

! convenience variables

nproc = num_proc(procs)
comm = slaves_communicator(procs)
nev = solver_cntl%num_eval
neq = linsys%neq
neq_full = full_matrix%neq
ss = grid%system_size
if (solver_cntl%lambda0 == -huge(0.0_my_real)) then
   lambda0 = linsys%rmin
else
   lambda0 = solver_cntl%lambda0
endif

! Make sure there is space in grid for the requested number of eigenvalues

if (size(grid%eigenvalue) /= nev) then
   deallocate(grid%eigenvalue,grid%eigenprob_l2_resid, &
              grid%eigenprob_variance,stat=allocstat)
   allocate(grid%eigenvalue(nev),grid%eigenprob_l2_resid(nev), &
            grid%eigenprob_variance(nev),stat=allocstat)
   if (allocstat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in eigen_arpack",procs=procs)
      return
   endif
   grid%eigenvalue = 0.0_my_real
   grid%eigenprob_l2_resid = 0.0_my_real
   grid%eigenprob_variance = 0.0_my_real
endif

! copy the solutions from the grid to the linear system to return the
! old solution if we abort

do i=1,neq_full
   call eq_to_grid(full_matrix,full_matrix%gid(i),objtype,brank,srank, &
                   objlid,grid)
   select case (objtype)
   case (VERTEX_ID)
      full_matrix%solution(i) = grid%vertex_solution(objlid,srank,1)
      if (nev > 1) then
         full_matrix%evecs(i,:) = grid%vertex_solution(objlid,srank,2:nev)
      endif
   case (EDGE_ID)
      full_matrix%solution(i) = grid%edge(objlid)%solution(brank,srank,1)
      if (nev > 1) then
         do j=1,nev-1
            full_matrix%evecs(i,j) = grid%edge(objlid)%solution(brank,srank,j)
         end do
      endif
   case (ELEMENT_ID)
      full_matrix%solution(i) = grid%element(objlid)%solution(brank,srank,1)
      if (nev > 1) then
         do j=1,nev-1
            full_matrix%evecs(i,j) = grid%element(objlid)%solution(brank,srank,j)
         end do
      endif
   end select
end do

! all processors must own an equation for ARPACK to work

if (phaml_global_min(procs,grid%nvert_own,2000) == 0) then
   call warning("Some processor owns no equations; ARPACK would fail.", &
                "Returning old solution.")
   return
endif

! if there aren't enough equations, reduce number of Arnoldi vectors

ncv = solver_cntl%arpack_cntl%ncv

neqown = 0
ndirch = 0
do i=1,neq_full
   if (full_matrix%iown(i)) then
      neqown = neqown + 1
      if (full_matrix%equation_type(i) == DIRICHLET) ndirch = ndirch + 1
   endif
end do
if (still_sequential) then
   neq_total = neqown - ndirch
else
   neq_total = phaml_global_sum(procs,neqown-ndirch,2001)
endif
if (neq_total < ncv) ncv = max(1,neq_total-1)

! if ncv is too small, reduce number of eigenvalues

if (ncv - nev < 2) then
   nev = ncv - 2
   call warning( &
     "Not enough equations or Arnoldi vectors for requested number of eigenvalues.", &
     "Computing", intlist=(/nev/))
endif

! allocate memory for ARPACK

lworkl = 3*ncv*ncv + 6*ncv
allocate(arnoldi_vec(neqown,ncv), arp_workd(3*neqown), arp_resid(neqown), &
         arp_workl(lworkl), dr(nev+1), di(nev+1), evec(neqown,nev+1), &
         workev(3*ncv), iperm(nev), stat=allocstat)
if (allocstat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in eigen_arpack",procs=procs)
   return
endif

! initalize convergence information

grid%eigen_linsys_max_l2_resid = 0.0_my_real
grid%eigen_linsys_ave_l2_resid = 0.0_my_real
nresid_comp = 0
grid%eigenprob_l2_resid = 0.0_my_real
grid%eigenprob_variance = 0.0_my_real

! create any needed alternate storage formats for solution methods

i_made_lapack = .false.
i_made_petsc = .false.
i_made_mumps = .false.
i_made_superlu = .false.
i_made_hypre = .false.
   
select case (solver_cntl%solver)

case (LAPACK_INDEFINITE_SOLVER)
   if (linsys%lapack_symm_band_exists) then
      call destroy_lapack_band(linsys%lapack_mat)
      linsys%lapack_symm_band_exists = .false.
   endif
   if (.not. linsys%lapack_gen_band_exists) then
      call make_lapack_gen_band(linsys%nlev+1,linsys,linsys%lapack_mat)
      linsys%lapack_gen_band_exists = .true.
      i_made_lapack = .true.
   endif

case (LAPACK_SPD_SOLVER)
   if (linsys%lapack_gen_band_exists) then
      call destroy_lapack_band(linsys%lapack_mat)
      linsys%lapack_gen_band_exists = .false.
   endif
   if (.not. linsys%lapack_symm_band_exists) then
      call make_lapack_symm_band(linsys%nlev+1,linsys,linsys%lapack_mat)
      if (ierr /= NO_ERROR) return
      linsys%lapack_symm_band_exists = .true.
      i_made_lapack = .true.
   endif

case (MUMPS_SPD_SOLVER, MUMPS_GEN_SOLVER, MUMPS_NONSYM_SOLVER)
   if (.not. linsys%mumps_matrix_exists) then
      call create_mumps_linear_system(linsys%mumps_matrix,linsys,procs, &
                                      solver_cntl,still_sequential)
      linsys%mumps_matrix_exists = .true.
      i_made_mumps = .true.
   endif

case (SUPERLU_SOLVER)
   if (.not. linsys%superlu_matrix_exists) then
      call create_superlu_linear_system(linsys%superlu_matrix,linsys, &
                                        procs,still_sequential)
      linsys%superlu_matrix_exists = .true.
      i_made_superlu = .true.
   endif

case (PETSC_RICHARDSON_SOLVER, PETSC_CHEBYCHEV_SOLVER, PETSC_CG_SOLVER, &
      PETSC_GMRES_SOLVER,      PETSC_TCQMR_SOLVER,     PETSC_BCGS_SOLVER,&
      PETSC_CGS_SOLVER,        PETSC_TFQMR_SOLVER,     PETSC_CR_SOLVER, &
      PETSC_LSQR_SOLVER,       PETSC_BICG_SOLVER)
   if (.not. linsys%petsc_matrix_exists) then
      if (solver_cntl%petsc_matrix_free) then
         call create_petsc_linear_system_mf(linsys,linsys%petsc_matrix, &
                                         still_sequential)
      else
         call create_petsc_linear_system(linsys,linsys%petsc_matrix, &
                                         still_sequential,procs)
      endif
      linsys%petsc_matrix_exists = .true.
      i_made_petsc = .true.
   endif

case (HYPRE_BOOMERAMG_SOLVER, HYPRE_PCG_SOLVER, HYPRE_GMRES_SOLVER)
   if (.not. linsys%hypre_matrix_exists) then
      call create_hypre_linear_system(linsys%hypre_matrix,linsys,procs, &
                                      still_sequential)
      linsys%hypre_matrix_exists = .true.
      i_made_hypre = .true.
   endif

end select

! chose between closest to lambda0 (LM), larger than lambda0 (LR),
! and smaller than lambda0 (SR) (assumes symm=.false., different if true)

select case (solver_cntl%lambda0_side)
case (EIGEN_BOTH)
   which = "LM"                         ! compute largest magnitude eigenvalues
case (EIGEN_RIGHT)
   which = "LR"                         ! compute largest real part eigenvalues
case (EIGEN_LEFT)
   which = "SR"                         ! compute smallest real part e-values
case default
   ierr = USER_INPUT_ERROR
   call fatal("bad value for lambda0_side in eigen_arpack", &
              intlist=(/solver_cntl%lambda0_side/),procs=procs)
   stop
end select

! set parameters for ARPACK

ido = 0                                 ! indicates first call
bmat = "G"                              ! generalized eigenvalue problem
tol = solver_cntl%arpack_cntl%tol       ! stopping criteria
arp_param = 0                           ! set unused entries to 0
arp_param(1) = 1                        ! shifts are not provided by the user
arp_param(3) = solver_cntl%arpack_cntl%maxit ! maximum number of iterations
arp_param(7) = 3                        ! shift-invert mode

! set initial guess for eigenvector
! TEMP would it be better to use some combination of existing eigenvectors
! instead of the first one?

call all_to_mine(full_matrix,full_matrix%solution(1:),arp_resid)
if (maxval(abs(arp_resid)) < 10*epsilon(0.0_my_real)) then
   info = 0                                 ! random initial residual vector
else
   info = 1                                 ! nonrandom initial residual vector
endif

! main loop

iter = 0
do
   iter = iter + 1

! call ARPACK routine [p| ][s|d][s|n]aupd for
! [parallel|sequential][single|double] precision [symmetric|nonsymmetric]
! eigensolver

   if (my_real == kind(1.0e0)) then
      if (nproc > 1 .and. .not. still_sequential) then
       if (symm) then
         call pssaupd(comm,ido,bmat,neqown,which,nev,tol,arp_resid,ncv, &
                     arnoldi_vec,neqown,arp_param,arp_pntr,arp_workd, &
                     arp_workl,lworkl,info)
       else
         call psnaupd(comm,ido,bmat,neqown,which,nev,tol,arp_resid,ncv, &
                     arnoldi_vec,neqown,arp_param,arp_pntr,arp_workd, &
                     arp_workl,lworkl,info)
       endif
      else
       if (symm) then
         call ssaupd(ido,bmat,neqown,which,nev,tol,arp_resid,ncv, &
                     arnoldi_vec,neqown,arp_param,arp_pntr,arp_workd, &
                     arp_workl,lworkl,info)
       else
         call snaupd(ido,bmat,neqown,which,nev,tol,arp_resid,ncv, &
                     arnoldi_vec,neqown,arp_param,arp_pntr,arp_workd, &
                     arp_workl,lworkl,info)
       endif
      endif
   elseif (my_real == kind(1.0d0)) then
      if (nproc > 1 .and. .not. still_sequential) then
       if (symm) then
         call pdsaupd(comm,ido,bmat,neqown,which,nev,tol,arp_resid,ncv, &
                     arnoldi_vec,neqown,arp_param,arp_pntr,arp_workd, &
                     arp_workl,lworkl,info)
       else
         call pdnaupd(comm,ido,bmat,neqown,which,nev,tol,arp_resid,ncv, &
                     arnoldi_vec,neqown,arp_param,arp_pntr,arp_workd, &
                     arp_workl,lworkl,info)
       endif
      else
       if (symm) then
         call dsaupd(ido,bmat,neqown,which,nev,tol,arp_resid,ncv, &
                     arnoldi_vec,neqown,arp_param,arp_pntr,arp_workd, &
                     arp_workl,lworkl,info)
       else
         call dnaupd(ido,bmat,neqown,which,nev,tol,arp_resid,ncv, &
                     arnoldi_vec,neqown,arp_param,arp_pntr,arp_workd, &
                     arp_workl,lworkl,info)
       endif
      endif
   else
      ierr = PHAML_INTERNAL_ERROR
      call fatal("my_real is neither single nor double precision. Can't call ARPACK",procs=procs)
      return
   endif

   if (info /= 0) then
      call warning("info from ARPACK *aupd is nonzero",intlist=(/info/))
   endif

! Determine which task to perform

   select case (ido)

   case (-1,1) ! multiply by (A-lambda0*M)^-1*M
               ! If ido==1, then leave off the multiplication by M

! Copy the vector from ARPACK into rhs

      if (ido == -1) then
         call mine_to_all(full_matrix,arp_workd(arp_pntr(1):arp_pntr(1)+neqown-1), &
                          full_matrix%rhs,0.0_my_real)
      else
         call mine_to_all(full_matrix,arp_workd(arp_pntr(3):arp_pntr(3)+neqown-1), &
                          full_matrix%rhs,0.0_my_real)
      endif

! For FuDoP multigrid, get rhs entries from other processors, making sure they
! match my level of the corresponding vertex

      if (solver_cntl%solver == MG_SOLVER .and. .not. still_sequential .and. &
          nproc > 1) then
         do lev=full_matrix%nlev+1,2,-1
            call basis_change(lev,TO_HIER,full_matrix,skip_matrix=.true.)
         end do
         call sum_fudop_vect(full_matrix%rhs,procs,full_matrix, &
                             2310+iter+1,2320+iter+1,2325+iter+1)
         do lev=2,full_matrix%nlev+1
            call basis_change(lev,TO_NODAL,full_matrix,skip_matrix=.true.)
         end do
! second mine_to_all to fix any changes at my points due to basis changes
         if (ido == -1) then
            call mine_to_all(full_matrix,arp_workd(arp_pntr(1):arp_pntr(1)+neqown-1), &
                             full_matrix%rhs)
         else
            call mine_to_all(full_matrix,arp_workd(arp_pntr(3):arp_pntr(3)+neqown-1), &
                             full_matrix%rhs)
         endif
      endif

! for ido==-1, multiply by M putting the result in rhs

      if (ido == -1) then
         tempvec = full_matrix%rhs
         full_matrix%matrix_val => full_matrix%mass
         call matrix_times_vector(tempvec,full_matrix%rhs,full_matrix,procs, &
                                  still_sequential,2100+iter,2110+iter, &
                                  2115+iter,2120+iter,2130+iter,2135+iter, &
                                  nocomm2=.true.)
! For FuDoP multigrid, get the product for equations owned by other processors,
! making sure it is converted to my level for that equation
         if (solver_cntl%solver == MG_SOLVER .and. &
             .not. still_sequential .and. nproc > 1) then
            do i=1,neq
               if (.not. full_matrix%iown(i)) &
                  full_matrix%rhs(i) = 0.0_my_real
            end do
            do lev=full_matrix%nlev+1,2,-1
               call basis_change(lev,TO_HIER,full_matrix,skip_matrix=.true.)
            end do
            call sum_fudop_vect(full_matrix%rhs,procs,full_matrix, &
                                2330+2*iter,2340+2*iter,2345+2*iter)
            do lev=2,full_matrix%nlev+1
               call basis_change(lev,TO_NODAL,full_matrix,skip_matrix=.true.)
            end do
         endif
      endif

! multiply by (A-lambda0*M)^-1, i.e., call a solver

! perform static condensation on the rhs of full_matrix, putting the
! result in the rhs of linsys; b1 - A21 A2^-1 b2 where 2 is face bases

      if (full_matrix%begin_level(full_matrix%nlev+2) /= &
          full_matrix%begin_level(full_matrix%nlev+3)) then

         linsys%rhs(1:neq) = full_matrix%rhs(1:neq)
         full_matrix%matrix_val => full_matrix%shifted
         eq = full_matrix%begin_level(full_matrix%nlev+2)
         do
            if (eq >= full_matrix%begin_level(full_matrix%nlev+3)) exit
            call eq_to_grid(full_matrix,full_matrix%gid(eq),objtype,brank, &
                            srank,objlid,grid)
            block => full_matrix%elem_block(objlid)%matrix
            ipiv => full_matrix%elem_block(objlid)%ipiv
            loc_neq = full_matrix%elem_block(objlid)%neq
            tempvec(1:loc_neq) = full_matrix%rhs(eq:eq+loc_neq-1)
            if (my_real == kind(0.0)) then
              call sgetrs("N",loc_neq,1,block,loc_neq,ipiv,tempvec,loc_neq,info)
            else 
              call dgetrs("N",loc_neq,1,block,loc_neq,ipiv,tempvec,loc_neq,info)
            endif 
            if (info /= 0) then
               ierr = PHAML_INTERNAL_ERROR
               call fatal("lapack sgetrs solution failed", &
                          intlist=(/info/),procs=procs)
               stop
            endif
            do i=1,loc_neq
               do j=full_matrix%begin_row(eq+i-1),full_matrix%end_row(eq+i-1)
                  if (full_matrix%column_index(j) == NO_ENTRY) cycle
                  if (full_matrix%column_index(j) >= eq .and. &
                      full_matrix%column_index(j) < eq+loc_neq) cycle
                  linsys%rhs(full_matrix%column_index(j)) = &
                                linsys%rhs(full_matrix%column_index(j)) - &
                                get_matval(full_matrix,full_matrix%matrix_val, &
                                full_matrix%column_index(j),eq+i-1)*tempvec(i)
               end do
            end do
            eq = eq + loc_neq
         end do

         if (.not. still_sequential .and. nproc > 1) then
            call sum_fudop_vect(linsys%rhs(1:neq),procs,linsys,2510+iter, &
                                2520+iter,2530+iter)
         endif

      endif

! solve the linear system by the chosen solver

      select case (solver_cntl%solver)
      case (MG_SOLVER)
         call multigrid(grid,procs,linsys,io_cntl,solver_cntl,still_sequential,&
                        no_master=.true.)
      case (CG_SOLVER)
         call phaml_cg(linsys,procs,grid,io_cntl,solver_cntl,still_sequential)
      case (GMRES_SOLVER)
         call phaml_gmres(linsys,procs,grid,io_cntl,solver_cntl, &
                          still_sequential)
      case (LAPACK_INDEFINITE_SOLVER)
         call lapack_indef(linsys%nlev+1,linsys,linsys%lapack_mat)
      case (LAPACK_SPD_SOLVER)
         call lapack_spd(linsys%nlev+1,linsys,linsys%lapack_mat)
      case (MUMPS_SPD_SOLVER, MUMPS_GEN_SOLVER, MUMPS_NONSYM_SOLVER)
         call change_mumps_rhs(linsys%mumps_matrix,linsys,linsys%rhs, &
                               procs,still_sequential)
         call mumps_solve(linsys%mumps_matrix,linsys,procs, &
                          still_sequential,noshadow=.true.)
      case (SUPERLU_SOLVER)
         call change_superlu_rhs(linsys%superlu_matrix,linsys,linsys%rhs, &
                                 procs,still_sequential)
         call superlu_solve(linsys%superlu_matrix,linsys,procs, &
                            still_sequential)
      case (PETSC_RICHARDSON_SOLVER, PETSC_CHEBYCHEV_SOLVER, PETSC_CG_SOLVER, &
            PETSC_GMRES_SOLVER,      PETSC_TCQMR_SOLVER,     PETSC_BCGS_SOLVER,&
            PETSC_CGS_SOLVER,        PETSC_TFQMR_SOLVER,     PETSC_CR_SOLVER, &
            PETSC_LSQR_SOLVER,       PETSC_BICG_SOLVER)
         where (linsys%equation_type == DIRICHLET) linsys%rhs = linsys%solution(1:)
         call change_petsc_rhs(linsys%rhs,linsys%petsc_matrix)
         call petsc_solve(linsys,linsys%petsc_matrix,solver_cntl,io_cntl, &
                          still_sequential,grid,procs)
      case (HYPRE_BOOMERAMG_SOLVER, HYPRE_PCG_SOLVER, HYPRE_GMRES_SOLVER)
         call change_hypre_rhs(linsys%hypre_matrix,linsys,linsys%rhs,procs, &
                               still_sequential)
         call zero_hypre_solution(linsys%hypre_matrix,procs,still_sequential)
         call hypre_solve(linsys%hypre_matrix,linsys,procs,solver_cntl, &
                          still_sequential)
      case default
         ierr = USER_INPUT_ERROR
         call fatal("illegal selection for solver",procs=procs)
         return
      end select

! get solution for unowned equations that neigbor an owned equation

      if (.not. still_sequential) then
         call exchange_neigh_vect(linsys%solution(1:neq),procs,linsys, &
                                  2360+iter,2370+iter,2380+iter)
      endif

! solve for the face bases; x2 = A2^-1 (b2 - A12 x1)
! note full_matrix%solution and linsys%solution point to the same place

      full_matrix%matrix_val => full_matrix%shifted
      eq = full_matrix%begin_level(full_matrix%nlev+2)
      do
         if (eq >= full_matrix%begin_level(full_matrix%nlev+3)) exit
         call eq_to_grid(full_matrix,full_matrix%gid(eq),objtype,brank, &
                         srank,objlid,grid)
         block => full_matrix%elem_block(objlid)%matrix
         ipiv => full_matrix%elem_block(objlid)%ipiv
         loc_neq = full_matrix%elem_block(objlid)%neq
         do i=1,loc_neq
            full_matrix%solution(eq+i-1) = full_matrix%rhs(eq+i-1)
            do j=full_matrix%begin_row(eq+i-1),full_matrix%end_row(eq+i-1)
               if (full_matrix%column_index(j) == NO_ENTRY) cycle
               if (full_matrix%column_index(j) >= eq .and. &
                   full_matrix%column_index(j) < eq+loc_neq) cycle
               full_matrix%solution(eq+i-1) = full_matrix%solution(eq+i-1)- &
                  full_matrix%matrix_val(j)*full_matrix%solution(full_matrix%column_index(j))
            end do
         end do
         if (my_real == kind(0.0)) then
            call sgetrs("N",loc_neq,1,block,loc_neq,ipiv, &
                        full_matrix%solution(eq:eq+loc_neq-1),loc_neq,info)
         else
            call dgetrs("N",loc_neq,1,block,loc_neq,ipiv, &
                        full_matrix%solution(eq:eq+loc_neq-1),loc_neq,info)
         endif
         if (info /= 0) then
            ierr = PHAML_INTERNAL_ERROR
            call fatal("lapack spotrs solution failed for solving for face bases", &
                       intlist=(/info/),procs=procs)
            stop
         endif
         eq = eq + loc_neq
      end do

! compute the residual from the linear system solver (does not include
! face bases)

      my_resid = 0.0_my_real
      call matrix_times_vector(linsys%solution(1:neq),my_resid,linsys, &
                               procs,still_sequential,2410+iter,2420+iter, &
                               2425+iter,2430+iter,2440+iter,2445+iter, &
                               nocomm2=.true.)
      my_resid = my_resid-linsys%rhs(1:neq)
      do i=1,neq
         if (.not. linsys%iown(i)) my_resid(i) = 0.0_my_real
         if (linsys%equation_type(i)==DIRICHLET) my_resid(i) = 0.0_my_real
      end do
      l2resid = sum(my_resid**2)
      l2rhs = sum(linsys%rhs(1:neq)**2, &
                  mask=linsys%equation_type(1:neq)/=DIRICHLET)
      sums = phaml_global_sum(procs,(/l2resid,l2rhs/),2450+iter)
      l2resid = sqrt(sums(1))
      l2rhs = sqrt(sums(2))
      if (l2rhs /= 0.0_my_real) l2resid = l2resid/l2rhs
      grid%eigen_linsys_max_l2_resid = max(grid%eigen_linsys_max_l2_resid, &
                                           l2resid)
      grid%eigen_linsys_ave_l2_resid = grid%eigen_linsys_max_l2_resid + l2resid
      nresid_comp = nresid_comp + 1

! copy the result into ARPACK's workspace

      call all_to_mine(full_matrix,full_matrix%solution(1:neq_full), &
                       arp_workd(arp_pntr(2):arp_pntr(2)+neqown-1))

   case (2) ! multiply by M

! Copy the vector from ARPACK into rhs

      call mine_to_all(full_matrix,arp_workd(arp_pntr(1):arp_pntr(1)+neqown-1), &
                       full_matrix%rhs,0.0_my_real)

! For FuDoP multigrid, get rhs entries from other processors, making sure they
! match my level of the corresponding vertex

      if (solver_cntl%solver == MG_SOLVER .and. .not. still_sequential .and. &
          nproc > 1) then
         do lev=full_matrix%nlev+1,2,-1
            call basis_change(lev,TO_HIER,full_matrix,skip_matrix=.true.)
         end do
         call sum_fudop_vect(full_matrix%rhs,procs,full_matrix, &
                             2350+iter+1,2360+iter+1,2365+iter+1)
         do lev=2,full_matrix%nlev+1
            call basis_change(lev,TO_NODAL,full_matrix,skip_matrix=.true.)
         end do
! second mine_to_all to fix any changes at my points due to basis changes
         call mine_to_all(full_matrix,arp_workd(arp_pntr(1):arp_pntr(1)+neqown-1), &
                          full_matrix%rhs)
      endif

! multiply by M putting the result in rhs

      tempvec = full_matrix%rhs
      full_matrix%matrix_val => full_matrix%mass
      call matrix_times_vector(tempvec,full_matrix%rhs,full_matrix,procs, &
                               still_sequential,2140+iter,2150+iter, &
                               2155+iter,2160+iter,2170+iter,2175+iter, &
                               nocomm2=.true.)

! copy the result into ARPACK's workspace

      call all_to_mine(full_matrix,full_matrix%rhs, &
                       arp_workd(arp_pntr(2):arp_pntr(2)+neqown-1))

   case default ! ARPACK done
      exit

   end select

end do

! finish computing the average residual

if (nresid_comp /= 0) then
   grid%eigen_linsys_ave_l2_resid = grid%eigen_linsys_ave_l2_resid/nresid_comp
endif

! gather some statistics from ARPACK

grid%arpack_iter = arp_param(3)
grid%arpack_nconv = arp_param(5)
grid%arpack_numop = arp_param(9)
grid%arpack_numopb = arp_param(10)
grid%arpack_numreo = arp_param(11)
grid%arpack_info = 0

! check for failure

if (info /= 0) then
   call warning("ARPACK *aupd failed. info and nconv are",&
                intlist=(/info,grid%arpack_nconv/))
   grid%arpack_info = info
endif

! extract eigenvalue and compute eigenvector

rvec = .true. ! compute eigenvector
howmany = "A" ! compute all eigenvectors

kev = nev
if (my_real == kind(1.0e0)) then
   if (nproc > 1 .and. .not. still_sequential) then
    if (symm) then
      call psseupd(comm,rvec,howmany,select,dr,evec,neqown,lambda0, &
                  bmat,neqown,which,kev,tol,arp_resid, &
                  ncv,arnoldi_vec,neqown,arp_param,arp_pntr,arp_workd, &
                  arp_workl,lworkl,info)
    else
      call psneupd(comm,rvec,howmany,select,dr,di,evec,neqown,lambda0, &
                  0.0_my_real,workev,bmat,neqown,which,kev,tol,arp_resid, &
                  ncv,arnoldi_vec,neqown,arp_param,arp_pntr,arp_workd, &
                  arp_workl,lworkl,info)
    endif
   else
    if (symm) then
      call sseupd(rvec,howmany,select,dr,evec,neqown,lambda0, &
                  bmat,neqown,which,kev,tol,arp_resid, &
                  ncv,arnoldi_vec,neqown,arp_param,arp_pntr,arp_workd, &
                  arp_workl,lworkl,info)
    else
      call sneupd(rvec,howmany,select,dr,di,evec,neqown,lambda0, &
                  0.0_my_real,workev,bmat,neqown,which,kev,tol,arp_resid, &
                  ncv,arnoldi_vec,neqown,arp_param,arp_pntr,arp_workd, &
                  arp_workl,lworkl,info)
    endif
   endif
elseif (my_real == kind(1.0d0)) then
   if (nproc > 1 .and. .not. still_sequential) then
    if (symm) then
      call pdseupd(comm,rvec,howmany,select,dr,evec,neqown,lambda0, &
                  bmat,neqown,which,kev,tol,arp_resid, &
                  ncv,arnoldi_vec,neqown,arp_param,arp_pntr,arp_workd, &
                  arp_workl,lworkl,info)
    else
      call pdneupd(comm,rvec,howmany,select,dr,di,evec,neqown,lambda0, &
                  0.0_my_real,workev,bmat,neqown,which,kev,tol,arp_resid, &
                  ncv,arnoldi_vec,neqown,arp_param,arp_pntr,arp_workd, &
                  arp_workl,lworkl,info)
    endif
   else
    if (symm) then
      call dseupd(rvec,howmany,select,dr,evec,neqown,lambda0, &
                  bmat,neqown,which,kev,tol,arp_resid, &
                  ncv,arnoldi_vec,neqown,arp_param,arp_pntr,arp_workd, &
                  arp_workl,lworkl,info)
    else
      call dneupd(rvec,howmany,select,dr,di,evec,neqown,lambda0, &
                  0.0_my_real,workev,bmat,neqown,which,kev,tol,arp_resid, &
                  ncv,arnoldi_vec,neqown,arp_param,arp_pntr,arp_workd, &
                  arp_workl,lworkl,info)
    endif
   endif
else
   ierr = PHAML_INTERNAL_ERROR
   call fatal("my_real is neither single nor double precision. Can't call ARPACK",procs=procs)
   return
endif

! check for failure

if (info /= 0) then
   call warning("ARPACK *eupd failed.",&
                intlist=(/info/))
   grid%arpack_info = info
endif

! sort the eigenvalues

call sort(dr,nev,iperm,1,info)

! For each eigenvalue/eigenvector

do jj=1,nev
   ii = iperm(jj)

! use tmpptr to point to the solution vector in the linear system data struct

   if (jj == 1) then
      tmpptr => full_matrix%solution(1:neq_full)
   else
      tmpptr => full_matrix%evecs(:,jj-1)
   endif

! copy the eigenvectors into the solutions

   call mine_to_all(full_matrix,evec(:,ii),tmpptr,0.0_my_real)
   if (.not. still_sequential) then
      call exchange_fudop_vect(tmpptr,procs,full_matrix,2060+ii,2070+ii, &
                               2075+ii)
   endif

! find the min and max solution values

   maxsoln = -huge(0.0_my_real)
   minsoln = huge(0.0_my_real)
   do i=1,neq_full
      if (full_matrix%equation_type(i) == DIRICHLET) cycle
      if (.not. full_matrix%iown(i)) cycle
      maxsoln = max(maxsoln,tmpptr(i))
      minsoln = min(minsoln,tmpptr(i))
   end do
   if (.not. still_sequential) then
      maxsoln = phaml_global_max(procs,maxsoln,2920+ii)
      minsoln = phaml_global_min(procs,minsoln,2930+ii)
   endif

! Scale the eigenvectors so the selected norm is 1.0 and so the maximum
! absolution value is in the positive direction

   select case (solver_cntl%scale_evec)

   case (SCALE_LINF)

! Scale eigenvectors to have L infinity norm 1.0

      if (abs(minsoln) > abs(maxsoln)) then
         if (minsoln /= 0.0_my_real) tmpptr = tmpptr/minsoln
      else
         if (maxsoln /= 0.0_my_real) tmpptr = tmpptr/maxsoln
      endif

   case (SCALE_L2)

! Scale eigenvectors to have L2 norm 1.0

      soln2norm = 0.0_my_real
      do i=1,neq_full
         if (full_matrix%equation_type(i) == DIRICHLET) cycle
         if (.not. full_matrix%iown(i)) cycle
         soln2norm = soln2norm + tmpptr(i)**2
      end do
      if (.not. still_sequential) then
         soln2norm = phaml_global_sum(procs,soln2norm,2810+ii)
      endif
      soln2norm = sqrt(soln2norm)
      if (abs(minsoln) > abs(maxsoln)) then
         soln2norm = -soln2norm
      endif
      if (soln2norm /= 0.0_my_real) tmpptr = tmpptr/soln2norm

   case (SCALE_M)

! Scale eigenvectors to have M-norm 1.0

      full_matrix%matrix_val => full_matrix%mass
      call matrix_times_vector(tmpptr,Mx,full_matrix,procs,still_sequential, &
                               2850+ii,2860+ii,2865+ii,2870+ii,2880+ii, &
                               2885+ii,nocomm2=.true.)
      solnMnorm = 0.0_my_real
      do i=1,neq_full
         if (full_matrix%equation_type(i) == DIRICHLET) cycle
         if (.not. full_matrix%iown(i)) cycle
         solnMnorm = solnMnorm + tmpptr(i)*Mx(i)
      end do
      if (.not. still_sequential) then
         solnMnorm = phaml_global_sum(procs,solnMnorm,2810+ii)
      endif
      solnMnorm = sqrt(solnMnorm)
      if (abs(minsoln) > abs(maxsoln)) then
         solnMnorm = -solnMnorm
      endif
      if (solnMnorm /= 0.0_my_real) tmpptr = tmpptr/solnMnorm

   end select

! Compute Ax and Mx

   Ax = 0.0_my_real
   Mx = 0.0_my_real
   full_matrix%matrix_val => full_matrix%stiffness
   call matrix_times_vector(tmpptr,Ax,full_matrix,procs,still_sequential, &
                            2210+ii,2220+ii,2225+ii,2230+ii,2240+ii, &
                            2245+ii,nocomm2=.true.)
   full_matrix%matrix_val => full_matrix%mass
   call matrix_times_vector(tmpptr,Mx,full_matrix,procs,still_sequential, &
                            2250+ii,2260+ii,2265+ii,2270+ii,2280+ii,2285+ii, &
                            nocomm2=.true.)

! Compute the Raleigh quotient x^T A x / x^T M x to use as eigenvalues

   numer = 0.0_my_real
   denom = 0.0_my_real
   do i=1,neq_full
      if (full_matrix%equation_type(i) == DIRICHLET) cycle
      if (.not. full_matrix%iown(i)) cycle
      numer = tmpptr(i)*Ax(i) + numer
      denom = tmpptr(i)*Mx(i) + denom
   end do
   if (.not. still_sequential) then
      numer = phaml_global_sum(procs,numer,2020+ii)
      denom = phaml_global_sum(procs,denom,2030+ii)
   endif
   grid%eigenvalue(jj) = numer/denom

! If the Raleigh quotient is not close to the eigenvalue that ARPACK returned,
! it may be wrong.  Tell the master to print a warning.

   if (dr(ii) /= 0.0_my_real) then
      if (abs((dr(ii)-grid%eigenvalue(jj))/dr(ii)) > 1.0_my_real) then
         if (my_proc(procs) == 1) then
            call phaml_send(procs,MASTER,(/2,jj/),2, &
                            (/dr(ii),grid%eigenvalue(jj)/),2,2201)
         endif
      endif
   endif

! Compute the l2 norm of the residual || Ax - lambda Mx || / || lambda Mx ||

   l2resid = 0.0_my_real
   l2rhs = 0.0_my_real
   do i=1,neq_full
      if (full_matrix%equation_type(i) == DIRICHLET) cycle
      if (.not. full_matrix%iown(i)) cycle
      l2resid = l2resid + (Ax(i)-grid%eigenvalue(jj)*Mx(i))**2
      l2rhs = l2rhs + Mx(i)**2
   end do
   if (.not. still_sequential) then
      l2resid = phaml_global_sum(procs,l2resid,2040+ii)
      l2rhs = phaml_global_sum(procs,l2rhs,2050+ii)
   endif
   l2resid = sqrt(l2resid)
   l2rhs = abs(grid%eigenvalue(jj))*sqrt(l2rhs)
   if (l2rhs /= 0.0_my_real) l2resid = l2resid/l2rhs
   grid%eigenprob_l2_resid(jj) = l2resid

! the variance x^T A^T M^-1 A x / x^T M x

   numer = 0.0_my_real
   denom = 0.0_my_real
   do i=1,neq_full
      if (full_matrix%equation_type(i) == DIRICHLET) cycle
      if (.not. full_matrix%iown(i)) cycle
      numer = numer + Ax(i)**2
      denom = denom + Mx(i)**2
   end do
   if (.not. still_sequential) then
      numer = phaml_global_sum(procs,numer,2080+ii)
      denom = phaml_global_sum(procs,denom,2090+ii)
   endif
   if (denom /= 0.0_my_real) then
      grid%eigenprob_variance(jj) = numer/denom - grid%eigenvalue(jj)**2
   else
      grid%eigenprob_variance(jj) = 0.0_my_real
   endif

end do ! jj

! free workspace

deallocate(arnoldi_vec, arp_workd, arp_resid, arp_workl, dr, di, evec, &
           workev, iperm, stat=allocstat)
if (allocstat /= 0) then
   call warning("deallocation failed in eigen_arpack")
endif
if (i_made_lapack) then
   call destroy_lapack_band(linsys%lapack_mat)
   linsys%lapack_symm_band_exists = .false.
   linsys%lapack_gen_band_exists = .false.
endif
if (i_made_petsc) then
   call destroy_petsc_linear_system(linsys%petsc_matrix)
   linsys%petsc_matrix_exists = .false.
endif
if (i_made_mumps) then
   call destroy_mumps_linear_system(linsys%mumps_matrix,procs)
   linsys%mumps_matrix_exists = .false.
endif
if (i_made_superlu) then
   call destroy_superlu_linear_system(linsys%superlu_matrix,procs)
   linsys%superlu_matrix_exists = .false.
endif
if (i_made_hypre) then
   call destroy_hypre_linear_system(linsys%hypre_matrix,procs)
   linsys%hypre_matrix_exists = .false.
endif

! tell the master we're done

if (my_proc(procs) == 1) then
   call phaml_send(procs,MASTER,(/1/),1,(/0.0_my_real/),0,2201)
endif

end subroutine eigen_arpack

!        ----------
function get_matval(linear_system,matval,row,col)
!        ----------

!----------------------------------------------------
! This routine get the (row,col) matrix value in matval
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(in) :: linear_system
real(my_real), intent(in) :: matval(:)
integer, intent(in) :: row,col
real(my_real) :: get_matval
!----------------------------------------------------
! Local variables:

integer :: i

!----------------------------------------------------
! Begin executable code

! TEMP this is used to build the transpose block in static condensation.
!      If that's all it is used for, this can be made more efficient when
!      the matrix is symmetric, and maybe other times, too

! seach row for column col to get the value; return 0 if col is not found

get_matval = 0.0_my_real
do i=linear_system%begin_row(row), linear_system%end_row(row)
   if (linear_system%column_index(i) == col) then
      get_matval = matval(i)
      exit
   endif
end do

end function get_matval

!          ------------
subroutine eigen_blopex(grid,procs,linsys,io_cntl,solver_cntl,still_sequential)
!          ------------

!----------------------------------------------------
! This routine solves the eigenvalue problem Ax = lambda*Mx where A is 
! stiffness and M is mass.  It uses BLOPEX.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
type(proc_info), intent(in) :: procs
type(linsys_type), intent(inout), target :: linsys
type(io_options), intent(in) :: io_cntl
type(solver_options), intent(in) :: solver_cntl
logical, intent(in) :: still_sequential

!----------------------------------------------------
! Local variables:

logical :: i_made_hypre, i_made_petsc
integer :: i, j, ii, nproc, nev, neq, neq_full, ss, comm, objtype, brank, &
           srank, objlid, allocstat, info, indx, istrt, indx0
real(my_real) :: maxsoln, minsoln, lambda0, soln2norm, solnMnorm, l2resid, &
                 l2rhs, numer, denom, temp, temp2, temp3
type(linsys_type) :: full_matrix
integer :: ni, nr, proc
integer, pointer :: irecv(:)
real(my_real), pointer :: rrecv(:)
real(my_real) :: Ax(linsys%neq_vert+linsys%neq_edge+linsys%neq_face), &
                 Mx(linsys%neq_vert+linsys%neq_edge+linsys%neq_face)
real(my_real), pointer :: tmpptr(:)
integer :: iperm(solver_cntl%num_eval)
!----------------------------------------------------
! Begin executable code

! master hangs around just to see if it needs to print a warning

if (my_proc(procs) == MASTER) then
   do
      call phaml_recv(procs,proc,irecv,ni,rrecv,nr,2201)
      select case (irecv(1))
      case (1) ! all done
         deallocate(irecv,stat=allocstat)
         exit
      case (2) ! lobpcg_solve returned an error code
         call warning("BLOPEX returned error code",intlist=(/irecv(2),irecv(3)/))
         deallocate(irecv,stat=allocstat)
         if (associated(rrecv)) deallocate(rrecv,stat=allocstat)
      end select
   end do
   return
endif

! use linsys for solving the (statically condensed) shifted linear systems,
! and full_matrix for operations that need the whole matrix.  full_matrix
! points to entries in linsys with the following pointers/values different:
!  linsys%end_row => end_row_edge
!  linsys%matrix_val => condensed => shifted
!  linsys%rhs => rhs_cond 
!  linsys%neq = neq_vert+neq_edge

full_matrix = linsys

full_matrix%end_row => linsys%end_row_face
full_matrix%matrix_val => linsys%stiffness
full_matrix%rhs => linsys%rhs_nocond
full_matrix%neq = linsys%neq_vert + linsys%neq_edge + linsys%neq_face

! convenience variables

nproc = num_proc(procs)
comm = slaves_communicator(procs)
nev = solver_cntl%num_eval
neq = linsys%neq
neq_full = full_matrix%neq
ss = grid%system_size
if (solver_cntl%lambda0 == -huge(0.0_my_real)) then
   lambda0 = solver_cntl%lambda0
elseif (solver_cntl%lambda0_side == EIGEN_RIGHT) then
   lambda0 = -solver_cntl%lambda0
else
   lambda0 = solver_cntl%lambda0
endif

! Make sure there is space in grid for the requested number of eigenvalues

if (size(grid%eigenvalue) /= nev) then
   deallocate(grid%eigenvalue,grid%eigenprob_l2_resid, &
              grid%eigenprob_variance,stat=allocstat)
   allocate(grid%eigenvalue(nev),grid%eigenprob_l2_resid(nev), &
            grid%eigenprob_variance(nev),stat=allocstat)
   if (allocstat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("memory allocation failed in eigen_blopex",procs=procs)
      return
   endif
   grid%eigenvalue = 0.0_my_real
   grid%eigenprob_l2_resid = 0.0_my_real
   grid%eigenprob_variance = 0.0_my_real
endif

! copy the solutions from the grid to the linear system to return the
! old solution if we abort

do i=1,neq_full
   call eq_to_grid(full_matrix,full_matrix%gid(i),objtype,brank,srank, &
                   objlid,grid)
   select case (objtype)
   case (VERTEX_ID)
      full_matrix%solution(i) = grid%vertex_solution(objlid,srank,1)
      if (nev > 1) then
         full_matrix%evecs(i,:) = grid%vertex_solution(objlid,srank,2:nev)
      endif
   case (EDGE_ID)
      full_matrix%solution(i) = grid%edge(objlid)%solution(brank,srank,1)
      if (nev > 1) then
         do j=1,nev-1
            full_matrix%evecs(i,j) = grid%edge(objlid)%solution(brank,srank,j)
         end do
      endif
   case (ELEMENT_ID)
      full_matrix%solution(i) = grid%element(objlid)%solution(brank,srank,1)
      if (nev > 1) then
         do j=1,nev-1
            full_matrix%evecs(i,j) = grid%element(objlid)%solution(brank,srank,j)
         end do
      endif
   end select
end do

! create any needed alternate storage formats for solution methods

i_made_hypre = .false.
i_made_petsc = .false.

select case (solver_cntl%solver)

! TEMP for later expansion
!case (HYPRE_BOOMERAMG_SOLVER, HYPRE_PCG_SOLVER, HYPRE_GMRES_SOLVER)
!   if (.not. linsys%hypre_matrix_exists) then
!      call create_hypre_linear_system(linsys%hypre_matrix,linsys,procs, &
!                                      still_sequential)
!      linsys%hypre_matrix_exists = .true.
!      i_made_hypre = .true.
!   endif

case (PETSC_RICHARDSON_SOLVER, PETSC_CHEBYCHEV_SOLVER, PETSC_CG_SOLVER, &
      PETSC_GMRES_SOLVER,      PETSC_TCQMR_SOLVER,     PETSC_BCGS_SOLVER,&
      PETSC_CGS_SOLVER,        PETSC_TFQMR_SOLVER,     PETSC_CR_SOLVER, &
      PETSC_LSQR_SOLVER,       PETSC_BICG_SOLVER)
   if (my_proc(procs) /= MASTER) then
      where (linsys%equation_type == DIRICHLET) linsys%rhs = linsys%solution(1:)
      if (.not. linsys%petsc_matrix_exists) then
         if (solver_cntl%petsc_matrix_free) then
            call create_petsc_linear_system_mf(linsys,linsys%petsc_matrix, &
                                               still_sequential)
         else
            call create_petsc_linear_system(linsys,linsys%petsc_matrix, &
                                            still_sequential,procs)
         endif
         linsys%petsc_matrix_exists = .true.
         i_made_petsc = .true.
      endif
   endif

case default
   call fatal("BLOPEX currently requires a petsc solver",procs=procs)
   stop

end select

! solve the eigenvalue problem

select case (solver_cntl%solver)

! TEMP for later expansion
!case (HYPRE_BOOMERAMG_SOLVER, HYPRE_PCG_SOLVER, HYPRE_GMRES_SOLVER)
!
!   if (my_proc(procs) /= MASTER) then
!      call hypre_lobpcg_solve_f(linsys,full_matrix,linsys%hypre_matrix, &
!                                solver_cntl,io_cntl,still_sequential,grid,procs)
!   endif

case (PETSC_RICHARDSON_SOLVER, PETSC_CHEBYCHEV_SOLVER, PETSC_CG_SOLVER, &
      PETSC_GMRES_SOLVER,      PETSC_TCQMR_SOLVER,     PETSC_BCGS_SOLVER,&
      PETSC_CGS_SOLVER,        PETSC_TFQMR_SOLVER,     PETSC_CR_SOLVER, &
      PETSC_LSQR_SOLVER,       PETSC_BICG_SOLVER)
   if (my_proc(procs) /= MASTER) then
      call petsc_lobpcg_solve_f(linsys,full_matrix,linsys%petsc_matrix, &
                                solver_cntl,io_cntl,still_sequential,grid,procs)
   endif

case default
   call fatal("BLOPEX currently requires a petsc solver",procs=procs)
   stop

end select

! For each eigenvalue/eigenvector

do ii=1,nev

! use tmpptr to point to the solution vector in the linear system data struct

   if (ii == 1) then
      tmpptr => full_matrix%solution(1:neq_full)
   else
      tmpptr => full_matrix%evecs(:,ii-1)
   endif

! find the min and max solution values

   maxsoln = -huge(0.0_my_real)
   minsoln = huge(0.0_my_real)
   do i=1,neq_full
      if (full_matrix%equation_type(i) == DIRICHLET) cycle
      if (.not. full_matrix%iown(i)) cycle
      maxsoln = max(maxsoln,tmpptr(i))
      minsoln = min(minsoln,tmpptr(i))
   end do
   if (.not. still_sequential) then
      maxsoln = phaml_global_max(procs,maxsoln,2920+ii)
      minsoln = phaml_global_min(procs,minsoln,2930+ii)
   endif

! Scale the eigenvectors so the selected norm is 1.0 and so the maximum
! absolution value is in the positive direction

   select case (solver_cntl%scale_evec)

   case (SCALE_LINF)

! Scale eigenvectors to have L infinity norm 1.0

      if (abs(minsoln) > abs(maxsoln)) then
         if (minsoln /= 0.0_my_real) tmpptr = tmpptr/minsoln
      else
         if (maxsoln /= 0.0_my_real) tmpptr = tmpptr/maxsoln
      endif

   case (SCALE_L2)

! Scale eigenvectors to have L2 norm 1.0

      soln2norm = 0.0_my_real
      do i=1,neq_full
         if (full_matrix%equation_type(i) == DIRICHLET) cycle
         if (.not. full_matrix%iown(i)) cycle
         soln2norm = soln2norm + tmpptr(i)**2
      end do
      if (.not. still_sequential) then
         soln2norm = phaml_global_sum(procs,soln2norm,2810+ii)
      endif
      soln2norm = sqrt(soln2norm)
      if (abs(minsoln) > abs(maxsoln)) then
         soln2norm = -soln2norm
      endif
      if (soln2norm /= 0.0_my_real) tmpptr = tmpptr/soln2norm
   
   case (SCALE_M)
   
! Scale eigenvectors to have M-norm 1.0

      full_matrix%matrix_val => full_matrix%mass
      call matrix_times_vector(tmpptr,Mx,full_matrix,procs,still_sequential, &
                               2850+ii,2860+ii,2865+ii,2870+ii,2880+ii, &
                               2885+ii,nocomm2=.true.)
      solnMnorm = 0.0_my_real
      do i=1,neq_full
         if (full_matrix%equation_type(i) == DIRICHLET) cycle
         if (.not. full_matrix%iown(i)) cycle
         solnMnorm = solnMnorm + tmpptr(i)*Mx(i)
      end do 
      if (.not. still_sequential) then
         solnMnorm = phaml_global_sum(procs,solnMnorm,2810+ii)
      endif
      solnMnorm = sqrt(solnMnorm)
      if (abs(minsoln) > abs(maxsoln)) then
         solnMnorm = -solnMnorm
      endif
      if (solnMnorm /= 0.0_my_real) tmpptr = tmpptr/solnMnorm

   end select 

! Compute Ax and Mx

   Ax = 0.0_my_real
   Mx = 0.0_my_real
   full_matrix%matrix_val => full_matrix%stiffness
   call matrix_times_vector(tmpptr,Ax,full_matrix,procs,still_sequential, &
                            2210+ii,2220+ii,2225+ii,2230+ii,2240+ii, &
                            2245+ii,nodirch=.true.,nocomm2=.true.)
   if (lambda0 /= -huge(0.0_my_real) .and. &
       solver_cntl%lambda0_side == EIGEN_RIGHT) then
      Ax = -Ax
   endif
   full_matrix%matrix_val => full_matrix%mass
   call matrix_times_vector(tmpptr,Mx,full_matrix,procs,still_sequential, &
                            2250+ii,2260+ii,2265+ii,2270+ii,2280+ii,2285+ii, &
                            nodirch=.true.,nocomm2=.true.)

! Compute the Raleigh quotient x^T A x / x^T M x to use as eigenvalues
! This is necessary for SHIFT_SQUARE because we don't know which square root
! to use, but we could use the transformed eigenvalue for SHIFT_INVERT.

   numer = 0.0_my_real
   denom = 0.0_my_real
   do i=1,neq_full
      if (full_matrix%equation_type(i) == DIRICHLET) cycle
      if (.not. full_matrix%iown(i)) cycle
      numer = tmpptr(i)*Ax(i) + numer
      denom = tmpptr(i)*Mx(i) + denom
   end do
   if (.not. still_sequential) then
      numer = phaml_global_sum(procs,numer,2020+ii)
      denom = phaml_global_sum(procs,denom,2030+ii)
   endif
   grid%eigenvalue(ii) = numer/denom

! Compute the l2 norm of the residual || Ax - lambda Mx || / || lambda Mx ||

   l2resid = 0.0_my_real
   l2rhs = 0.0_my_real
   do i=1,neq_full
      if (full_matrix%equation_type(i) == DIRICHLET) cycle
      if (.not. full_matrix%iown(i)) cycle
      l2resid = l2resid + (Ax(i)-grid%eigenvalue(ii)*Mx(i))**2
      call eq_to_grid(full_matrix,full_matrix%gid(i),objtype,brank,srank, &
                      objlid,grid)
      l2rhs = l2rhs + Mx(i)**2
   end do
   if (.not. still_sequential) then
      l2resid = phaml_global_sum(procs,l2resid,2040+ii)
      l2rhs = phaml_global_sum(procs,l2rhs,2050+ii)
   endif
   l2resid = sqrt(l2resid)
   l2rhs = abs(grid%eigenvalue(ii))*sqrt(l2rhs)
   if (l2rhs /= 0.0_my_real) l2resid = l2resid/l2rhs
   grid%eigenprob_l2_resid(ii) = l2resid

! the variance x^T A^T M^-1 A x / x^T M x

   numer = 0.0_my_real
   denom = 0.0_my_real
   do i=1,neq_full
      if (full_matrix%equation_type(i) == DIRICHLET) cycle
      if (.not. full_matrix%iown(i)) cycle
      numer = numer + Ax(i)**2
      denom = denom + Mx(i)**2
   end do
   if (.not. still_sequential) then
      numer = phaml_global_sum(procs,numer,2080+ii)
      denom = phaml_global_sum(procs,denom,2090+ii)
   endif
   if (denom /= 0.0_my_real) then
      grid%eigenprob_variance(ii) = numer/denom - grid%eigenvalue(ii)**2
   else
      grid%eigenprob_variance(ii) = 0.0_my_real
   endif

end do ! ii

! sort the eigenvalues, if they are on the left or both sides of lambda0

if (nev > 1 .and. lambda0 /= -huge(0.0_my_real)) then
   allocate(tmpptr(neq_full))

   call sort(grid%eigenvalue,nev,iperm,1,info)

   do i=1,nev
      indx=abs(iperm(i))
      if (iperm(indx) > 0) iperm(indx) = -iperm(indx)
   end do
   do istrt = 1,nev
      if (iperm(istrt) > 0) cycle
      indx = istrt
      indx0 = indx
      temp = grid%eigenvalue(istrt)
      temp2 = grid%eigenprob_l2_resid(istrt)
      temp3 = grid%eigenprob_variance(istrt)
      if (istrt == 1) then
         tmpptr = full_matrix%solution(1:)
      else
         tmpptr = full_matrix%evecs(:,istrt-1)
      endif
      do while (iperm(indx) < 0)
         grid%eigenvalue(indx) = grid%eigenvalue(-iperm(indx))
         grid%eigenprob_l2_resid(indx) = grid%eigenprob_l2_resid(-iperm(indx))
         grid%eigenprob_variance(indx) = grid%eigenprob_variance(-iperm(indx))
         if (indx == 1) then
            full_matrix%solution(1:) = full_matrix%evecs(:,-iperm(indx)-1)
         elseif (-iperm(indx) == 1) then
            full_matrix%evecs(:,indx-1) = full_matrix%solution(1:)
         else
            full_matrix%evecs(:,indx-1) = full_matrix%evecs(:,-iperm(indx)-1)
         endif
         indx0 = indx
         iperm(indx) = -iperm(indx)
         indx = iperm(indx)
      end do
      grid%eigenvalue(indx0) = temp
      grid%eigenprob_l2_resid(indx0) = temp2
      grid%eigenprob_variance(indx0) = temp3
      if (indx0 == 1) then
         full_matrix%solution(1:) = tmpptr
      else
         full_matrix%evecs(:,indx0-1) = tmpptr
      endif
   end do

   deallocate(tmpptr)
endif

! get rid of any extra versions of the linear system

if (i_made_hypre) then
   call destroy_hypre_linear_system(linsys%hypre_matrix,procs)
   linsys%hypre_matrix_exists = .false.
endif
if (i_made_petsc) then
   call destroy_petsc_linear_system(linsys%petsc_matrix)
   linsys%petsc_matrix_exists = .false.
endif

! tell the master we're done

if (my_proc(procs) == 1) then
   call phaml_send(procs,MASTER,(/1/),1,(/0.0_my_real/),0,2201)
endif

end subroutine eigen_blopex

end module eigen
