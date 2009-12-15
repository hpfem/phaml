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

module phaml

!----------------------------------------------------
! This module contains the user callable routines of the
! phaml package.
! The routines are: phaml_solve_pde, phaml_create, phaml_destroy, 
!                   phaml_evaluate, phaml_evaluate_old, phaml_copy_soln_to_old,
!                   phaml_integrate, phaml_connect, phaml_store, phaml_restore,
!                   phaml_store_matrix, phaml_popen, phaml_pclose, phaml_query,
!                   phaml_scale, phaml_compress
! See doc/USER_GUIDE for a complete description of the interface
! to these routines.
!
! communication tags in this module are of the form 1xx, 2xx, 3xx and 4xxx
!----------------------------------------------------

use stopwatch
use global
use quadrature_rules
use message_passing
use hash_mod
use gridtype_mod
use evaluate
use hypretype_mod
use petsc_type_mod
use linsystype_mod
use zoltan_interf
use grid_io
use grid_mod
use linear_system
use load_balance
use error_estimators
use linsys_io
use phaml_type_mod

!----------------------------------------------------

implicit none
private
public evaluate_slave ! TEMP080428 for time dependent Schroedinger
public phaml_solution_type, my_real, phaml_solve_pde, phaml_create, &
       phaml_destroy, phaml_evaluate, phaml_evaluate_old, &
       phaml_copy_soln_to_old, phaml_integrate, &
       phaml_connect, phaml_store, phaml_restore, phaml_store_matrix, &
       phaml_popen, phaml_pclose, phaml_slave, phaml_query, phaml_scale, &
       phaml_compress, pde, &
       master_to_slaves, my_pde_id, &
       MASTER, MASTER_ALL, SLAVES, EVERYONE, NO_ONE, NEVER, FINAL, PHASES, &
       FREQUENTLY, TOO_MUCH, LAST, LAST_AND_FINAL, &
       ENERGY_ERR, LINF_ERR, L2_ERR, ENERGY_LINF_ERR, ENERGY_L2_ERR, &
       LINF_L2_ERR, ENERGY_LINF_L2_ERR, ENERGY_ERREST, ENERGY_LINF_ERREST, &
       ENERGY_L2_ERREST, ENERGY_LINF_L2_ERREST, LINF_ERREST, &
       LINF_L2_ERREST, L2_ERREST, ABSOLUTE_ERROR, RELATIVE_ERROR, &
       CLOCK_W, CLOCK_C, CLOCK_CW, CLOCK_WC, &
       INTERIOR, DIRICHLET, NATURAL, MIXED, PERIODIC, &
       BALANCE_NONE, BALANCE_ELEMENTS, BALANCE_VERTICES, BALANCE_EQUATIONS, &
       HIERARCHICAL_COEFFICIENT, TRUE_DIFF, LOCAL_PROBLEM_H, LOCAL_PROBLEM_P, &
       INITIAL_CONDITION, EXPLICIT_ERRIND, EQUILIBRATED_RESIDUAL, &
       RTK, ZOLTAN_RCB, ZOLTAN_OCT, ZOLTAN_METIS, ZOLTAN_REFTREE, ZOLTAN_RIB, &
       ZOLTAN_HSFC, ZOLTAN_FILE, BALANCE_REFINE_SOLVE, SET_INITIAL,  &
       BALANCE_ONLY, SOLVE_ONLY, REFINE_ONLY
public MG_NO_TOL, MG_ERREST_TOL, &
       MG_SOLVER, LAPACK_INDEFINITE_SOLVER, LAPACK_SPD_SOLVER, &
       CG_SOLVER, GMRES_SOLVER, &
       PETSC_RICHARDSON_SOLVER, PETSC_CHEBYCHEV_SOLVER, PETSC_CG_SOLVER, &
       PETSC_GMRES_SOLVER, PETSC_TCQMR_SOLVER, PETSC_BCGS_SOLVER, &
       PETSC_CGS_SOLVER, PETSC_TFQMR_SOLVER, PETSC_CR_SOLVER, &
       PETSC_LSQR_SOLVER, PETSC_BICG_SOLVER, MUMPS_SPD_SOLVER, &
       MUMPS_NONSYM_SOLVER, MUMPS_GEN_SOLVER, &
       SUPERLU_SOLVER, HYPRE_BOOMERAMG_SOLVER, &
       HYPRE_PCG_SOLVER, HYPRE_GMRES_SOLVER, &
       ARPACK_SOLVER, BLOPEX_SOLVER, TEST_PRECONDITION, &
       NO_PRECONDITION, MG_PRECONDITION, COARSE_GRID_PRECONDITION,&
       PETSC_JACOBI_PRECONDITION, PETSC_BJACOBI_PRECONDITION,&
       PETSC_SOR_PRECONDITION, PETSC_EISENSTAT_PRECONDITION, &
       PETSC_ICC_PRECONDITION, PETSC_ILU_PRECONDITION, PETSC_ASM_PRECONDITION, &
       FUDOP_DD_PRECONDITION, HYPRE_DS_PRECONDITION, &
       HYPRE_BOOMERAMG_PRECONDITION, HYPRE_PARASAILS_PRECONDITION, &
       MGCOMM_NONE, MGCOMM_FUDOP, MGCOMM_CONVENTIONAL, &
       ELLIPTIC, EIGENVALUE, NORMAL_SPAWN, DEBUG_SLAVE, DEBUG_GRAPHICS, &
       DEBUG_BOTH, DOUBLE_NVERT, DOUBLE_NVERT_SMOOTH, DOUBLE_NELEM, &
       DOUBLE_NELEM_SMOOTH, DOUBLE_NEQ, HALVE_ERREST, KEEP_NVERT, &
       KEEP_NVERT_SMOOTH, KEEP_NELEM, KEEP_NELEM_SMOOTH, KEEP_NEQ, KEEP_ERREST,&
       DOUBLE_NEQ_SMOOTH, KEEP_NEQ_SMOOTH, ONE_REF, ONE_REF_HALF_ERRIND, &
       H_UNIFORM, H_ADAPTIVE, P_UNIFORM, P_ADAPTIVE, HP_ADAPTIVE, &
       HP_BIGGER_ERRIND, HP_APRIORI, HP_PRIOR2P_E, HP_PRIOR2P_H1, HP_T3S, &
       HP_ALTERNATE, HP_TYPEPARAM, HP_COEF_DECAY, HP_COEF_ROOT, HP_SMOOTH_PRED,&
       HP_NEXT3P, HP_REFSOLN_EDGE, HP_REFSOLN_ELEM, HP_NLP, MINIMUM_RULE, &
       MAXIMUM_RULE, SCALE_LINF, SCALE_L2, SCALE_M, EIGEN_LEFT, EIGEN_RIGHT, &
       EIGEN_BOTH, SHIFT_INVERT, SHIFT_SQUARE

!----------------------------------------------------
! The following generic interfaces are defined:

!----------------------------------------------------
! The following variables are defined:

type(phaml_solution_type), allocatable :: pde(:)
!----------------------------------------------------
! Module procedures:

contains

!          ---------------
subroutine phaml_solve_pde(phaml_solution, iterm, max_elem, max_vert, max_eq, &
   max_lev, max_deg, stop_on_maxlev, stop_on_maxdeg, max_refsolveloop,        &
   term_energy_err, term_Linf_err, term_L2_err, task,                         &
   print_grid_when, print_grid_who, print_error_when, print_error_who,        &
   print_error_what, print_errest_what, print_linsys_when, print_linsys_who,  &
   print_time_when, print_time_who, print_eval_when, print_eval_who,          &
   print_header_who, print_trailer_who, print_warnings, clocks,draw_grid_when,&
   pause_after_draw, pause_after_phases, pause_at_start, pause_at_end,        &
   solve_init, sequential_vert, inc_factor, error_estimator, errtype,         &
   reftype, refterm, reftol, hp_strategy, t3s_gamma, t3s_eta, t3s_nunif,      &
   t3s_maxref, t3s_maxdeginc, tp_gamma, sp_gamma_h, sp_gamma_p,               &
   nlp_max_h_dec, nlp_max_h_inc, nlp_max_p_dec, nlp_max_p_inc,                &
   refsoln_pbias, derefine, partition_method, edge_rule,                      &
   zoltan_param_file, prebalance, postbalance, petsc_matrix_free, solver,     &
   preconditioner, mg_cycles, mg_tol, mg_prerelax, mg_postrelax,              &
   mg_prerelax_ho, mg_postrelax_ho, dd_iterations, krylov_iter,               &
   krylov_restart, krylov_tol, mg_comm, ignore_quad_err, eigensolver,         &
   num_eval, lambda0, lambda0_side, transformation, scale_evec, arpack_ncv,   &
   arpack_maxit, arpack_tol, blopex_maxit, blopex_atol, blopex_rtol,          &
   degree, inc_quad_order,                                                    &
   hypre_BoomerAMG_MaxLevels,hypre_BoomerAMG_MaxIter,hypre_BoomerAMG_Tol,     &
   hypre_BoomerAMG_StrongThreshold,hypre_BoomerAMG_MaxRowSum,                 &
   hypre_BoomerAMG_CoarsenType,hypre_BoomerAMG_MeasureType,                   &
   hypre_BoomerAMG_CycleType,hypre_BoomerAMG_NumGridSweeps,                   &
   hypre_BoomerAMG_GridRelaxType,hypre_BoomerAMG_GridRelaxPoints,             &
   hypre_BoomerAMG_RelaxWeight,                                               &
   hypre_BoomerAMG_DebugFlag,hypre_ParaSails_thresh,hypre_ParaSails_nlevels,  &
   hypre_ParaSails_filter,hypre_ParaSails_sym,hypre_ParaSails_loadbal,        &
   hypre_ParaSails_reuse,hypre_ParaSails_logging,hypre_PCG_Tol,               &
   hypre_PCG_MaxIter, hypre_PCG_TwoNorm,hypre_PCG_RelChange,hypre_PCG_Logging,&
   hypre_GMRES_KDim, hypre_GMRES_Tol,hypre_GMRES_MaxIter,hypre_GMRES_Logging, &
   petsc_richardson_damping_factor, petsc_chebychev_emin,                     &
   petsc_chebychev_emax, petsc_gmres_max_steps, petsc_rtol, petsc_atol,       &
   petsc_dtol, petsc_maxits, petsc_ilu_levels, petsc_icc_levels, petsc_ilu_dt,&
   petsc_ilu_dtcol, petsc_ilu_maxrowcount, petsc_sor_omega, petsc_sor_its,    &
   petsc_sor_lits, petsc_eisenstat_nodiagscaling, petsc_eisenstat_omega,      &
   petsc_asm_overlap,coarse_size, coarse_method,                              &
   hypre_has_NumGridSweeps, hypre_has_GridRelaxType,                          &
   hypre_has_GridRelaxPoints, hypre_has_RelaxWeight)
!          ---------------


!----------------------------------------------------
! This is the top level subroutine for PHAML.  It performs initialization
! and then sits in a refine/redistribute/solve loop until a termination
! condition is met.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type (phaml_solution_type), intent(inout), target :: phaml_solution
integer, optional, intent(out) :: iterm
logical, optional, intent(in) :: pause_after_draw, pause_at_start, &
                                 pause_at_end, pause_after_phases, &
                                 ignore_quad_err, petsc_matrix_free, &
                                 derefine, print_warnings, solve_init, &
                                 stop_on_maxlev, stop_on_maxdeg
integer, optional, intent(in) :: max_elem, max_vert, max_eq, max_lev, max_deg, &
           print_grid_when, print_error_when, print_time_when, print_eval_when,&
           print_linsys_when, print_linsys_who, print_grid_who, &
           print_error_who, print_time_who, print_eval_who, print_header_who, &
           print_trailer_who, print_error_what, print_errest_what, &
           clocks, draw_grid_when, sequential_vert, errtype, mg_comm, &
           error_estimator, mg_cycles, max_refsolveloop, partition_method, &
           prebalance, postbalance, degree, inc_quad_order, &
           task, solver, preconditioner, mg_prerelax, mg_postrelax, &
           mg_prerelax_ho, mg_postrelax_ho, dd_iterations, edge_rule, &
           reftype, refterm, hp_strategy, t3s_nunif, t3s_maxref, t3s_maxdeginc,&
           nlp_max_h_dec, nlp_max_h_inc, nlp_max_p_dec, nlp_max_p_inc, &
           eigensolver, num_eval, arpack_ncv, arpack_maxit, blopex_maxit, &
           coarse_size, coarse_method, lambda0_side, transformation, &
           scale_evec, krylov_iter, krylov_restart
real(my_real), optional, intent(in) :: term_energy_err, term_Linf_err, &
                                       term_L2_err, inc_factor, lambda0, &
                                       arpack_tol, blopex_atol, blopex_rtol, &
                                       mg_tol, reftol, krylov_tol, &
                                       t3s_gamma, t3s_eta, tp_gamma, &
                                       sp_gamma_h, sp_gamma_p, refsoln_pbias
character(len=*), optional, intent(in) :: zoltan_param_file

integer, optional, intent(in) :: hypre_BoomerAMG_MaxLevels, &
   hypre_BoomerAMG_MaxIter,hypre_BoomerAMG_CoarsenType, &
   hypre_BoomerAMG_MeasureType,hypre_BoomerAMG_CycleType, &
   hypre_BoomerAMG_DebugFlag,hypre_ParaSails_nlevels, &
   hypre_ParaSails_sym,hypre_ParaSails_reuse,hypre_ParaSails_logging, &
   hypre_PCG_MaxIter,hypre_PCG_TwoNorm,hypre_PCG_RelChange, &
   hypre_PCG_Logging,hypre_GMRES_KDim,hypre_GMRES_MaxIter,hypre_GMRES_Logging

real(my_real), optional, intent(in) :: hypre_BoomerAMG_Tol, &
   hypre_BoomerAMG_StrongThreshold,hypre_BoomerAMG_MaxRowSum, &
   hypre_ParaSails_thresh,hypre_ParaSails_filter,hypre_ParaSails_loadbal, &
   hypre_PCG_Tol,hypre_GMRES_Tol

integer, optional, intent(in) :: hypre_BoomerAMG_NumGridSweeps(:), &
   hypre_BoomerAMG_GridRelaxType(:),hypre_BoomerAMG_GridRelaxPoints(:,:)

real(my_real), optional, intent(in) :: hypre_BoomerAMG_RelaxWeight(:)

real(my_real), optional, intent(in) :: petsc_richardson_damping_factor, &
   petsc_chebychev_emin, petsc_chebychev_emax, petsc_rtol, petsc_atol,  &
   petsc_dtol, petsc_ilu_dt, petsc_ilu_dtcol, petsc_sor_omega,          &
   petsc_eisenstat_omega

integer, optional, intent(in) :: petsc_gmres_max_steps, petsc_maxits,    &
   petsc_ilu_levels, petsc_icc_levels, petsc_ilu_maxrowcount, petsc_sor_its, &
   petsc_sor_lits, petsc_asm_overlap

logical, optional, intent(in) :: petsc_eisenstat_nodiagscaling

! These are not for the user to provide.  They are set true by a slave
! process if the corresponding hypre option was present in the user call
! to the master process.
logical, optional, intent(in) :: hypre_has_NumGridSweeps, &
   hypre_has_GridRelaxType, hypre_has_GridRelaxPoints, hypre_has_RelaxWeight

!----------------------------------------------------
! Local variables

! local copies of optional dummy arguments
integer :: loc_clocks, loc_sequential_vert, loc_max_refsolveloop, &
           loc_partition_method, loc_task, loc_prebalance, &
           loc_postbalance, loc_print_header_who, &
           loc_print_trailer_who, loc_print_eval_when, loc_print_eval_who, &
           loc_degree
logical :: loc_solve_init, loc_stop_on_maxlev, loc_stop_on_maxdeg
type (io_options) :: io_control
type (solver_options) :: solver_control
type(refine_options) :: refine_control
integer :: send_int(4),ni,nr,proc,loop,loop_end_sequential
real (my_real) :: no_reals(1)
real (my_real), allocatable :: send_real(:)
integer, pointer :: recv_int(:)
real (my_real), pointer :: recv_real(:)
logical :: still_sequential, ltemp
integer :: total_elem,total_vert, total_nlev, total_dof, astat, nvert1, &
           nvert2, nelem1, nelem2, ndof1, ndof2, total_deg
integer :: i, init_nvert, init_nelem, init_dof, ierr1, ierr2, comp, evec, soln
real(my_real) :: errest_energy,errest_Linf,errest_L2,normsoln
character(len=5) :: which_time
type(grid_type), pointer :: grid
type(proc_info), pointer :: procs
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

outunit = phaml_solution%outunit
errunit = phaml_solution%errunit
my_pde_id = phaml_solution%pde_id
ierr = NO_ERROR

procs => phaml_solution%procs

! set default values based on the value of or absence of optional arguments,
! and print header

call defaults
call print_header

! if I am the master then send the parameters to the
! slaves and tell them to solve the pde

if (my_proc(procs) == MASTER) then
   call slaves_solve
endif

! initialize

call init(phaml_solution,procs,io_control,solver_control,refine_control, &
          loc_partition_method,loc_task,loc_degree)

! For the Texas 3 Step hp strategy, perform the uniform refinements

if (refine_control%reftype == HP_ADAPTIVE .and. &
    refine_control%hp_strategy == HP_T3S) then
   refine_control%t3s_reftype = H_UNIFORM
   i = refine_control%t3s_nunif
   do while (i > 0)
      call refine(phaml_solution%grid,procs,refine_control,solver_control, &
                  io_control,phaml_solution%still_sequential,0,0,0,0, &
                  loc_prebalance,loc_prebalance/=BALANCE_NONE)
      i = i-1
   end do
endif

! Solve on the initial grid and print the initial error

if (loc_solve_init .and. &
    (loc_task == SET_INITIAL .or. loc_task == BALANCE_REFINE_SOLVE)) then

   call solve(phaml_solution%grid,procs,io_control,solver_control, &
              phaml_solution%still_sequential,loc_task==SET_INITIAL)

! update the REFSOLN error estimates, if needed

   if (refine_control%reftype == HP_ADAPTIVE) then
      if (refine_control%hp_strategy == HP_REFSOLN_EDGE) then
         call ext_refsoln_errest_edge(phaml_solution%grid,procs, &
                                    refine_control,solver_control,io_control, &
                                    phaml_solution%still_sequential)
      elseif (refine_control%hp_strategy == HP_REFSOLN_ELEM) then
         call ext_refsoln_errest_elem(phaml_solution%grid,procs, &
                                    refine_control,solver_control,io_control, &
                                    phaml_solution%still_sequential)
      endif
   endif

   call print_error_info(phaml_solution%grid,procs,io_control, &
                         phaml_solution%still_sequential, &
                         (/PHASES,FREQUENTLY,TOO_MUCH/),240,reduction=.true., &
                         errest=refine_control%error_estimator)
   call draw_grid(phaml_solution%grid,procs,io_control,refine_control, &
                  phaml_solution%i_draw_grid,phaml_solution%master_draws_grid, &
                  phaml_solution%still_sequential, &
                  (/PHASES/),loc_partition_method,phaml_solution%lb)

endif

! For the Alternate and Texas 3 Step hp strategies, compute the
! initial error estimate, and set up for the first refinement

if (refine_control%reftype == HP_ADAPTIVE .and. &
    (refine_control%hp_strategy == HP_T3S .or. &
     refine_control%hp_strategy == HP_ALTERNATE)) then
   call error_estimate(phaml_solution%grid,procs, &
                       refine_control%error_estimator, &
                       errest_energy=errest_energy)
   errest_energy = sqrt(phaml_global_sum(procs,errest_energy**2,244))
   refine_control%refterm = ONE_REF
   refine_control%t3s_reftype = H_ADAPTIVE
   refine_control%t3s_p_target = max(refine_control%t3s_eta*errest_energy, &
                                     refine_control%term_energy_err)
   refine_control%t3s_h_target = &
                      refine_control%t3s_gamma*refine_control%t3s_p_target
   refine_control%reftol = refine_control%t3s_h_target
endif

call print_error_info(phaml_solution%grid,procs,io_control, &
                      phaml_solution%still_sequential, &
                      (/FREQUENTLY,TOO_MUCH/),241,reduction=.true., &
                      errest=refine_control%error_estimator)

! pause before starting

if (present(pause_at_start)) then
   call pause_until_enter(procs,pause_at_start,dont_pausewatch=.true.)
endif

! convenience variables

grid => phaml_solution%grid
still_sequential = phaml_solution%still_sequential
if (.not. still_sequential) then
   loop_end_sequential = 0
else
   loop_end_sequential = -1
endif

! save initial number of vertices, elements and equations

call get_grid_info(grid,procs,still_sequential,250, &
                   total_nvert=init_nvert,total_nelem_leaf=init_nelem, &
                   total_dof=init_dof)

! refine/partition/redistribute/solve until an error condition occurs

loop = 0
mainloop: do
   if (ierr /= NO_ERROR) exit
   call reset_watch(ptotal)
   loop = loop+1

! For the Texas 3 Step hp strategy, change reftype

   if (refine_control%reftype == HP_ADAPTIVE .and. &
       refine_control%hp_strategy == HP_T3S) then
      if (loop > 1) then
         call error_estimate(phaml_solution%grid,procs, &
                             refine_control%error_estimator, &
                             errest_energy=errest_energy)
         if (refine_control%t3s_reftype == P_ADAPTIVE) then
            refine_control%t3s_p_target = max(refine_control%t3s_eta*errest_energy, &
                                              refine_control%term_energy_err)
            refine_control%t3s_h_target = &
                            refine_control%t3s_gamma*refine_control%t3s_p_target
            refine_control%t3s_reftype = H_ADAPTIVE
         else
            refine_control%t3s_reftype = P_ADAPTIVE
         endif
      endif
   endif

! For the Alternate hp strategy, determine if it is time to change reftype

   if (refine_control%reftype == HP_ADAPTIVE .and. &
       refine_control%hp_strategy == HP_ALTERNATE) then
      call error_estimate(phaml_solution%grid,procs, &
                          refine_control%error_estimator, &
                          errest_energy=errest_energy)
      if (refine_control%t3s_reftype == P_ADAPTIVE .and. &
          errest_energy <= refine_control%t3s_p_target) then
         refine_control%t3s_p_target = max(refine_control%t3s_eta*errest_energy, &
                                           refine_control%term_energy_err)
         refine_control%t3s_h_target = &
                            refine_control%t3s_gamma*refine_control%t3s_p_target
         refine_control%t3s_reftype = H_ADAPTIVE
         refine_control%reftol = refine_control%t3s_h_target
      elseif (refine_control%t3s_reftype == H_ADAPTIVE .and. &
              errest_energy <= refine_control%t3s_h_target) then
         refine_control%t3s_reftype = P_ADAPTIVE
         refine_control%reftol = refine_control%t3s_p_target
      endif
   endif

! if still building the grid sequentially, see if it's big enough to distribute

   if (still_sequential .and. PARALLEL/=SEQUENTIAL) then
      call check_end_sequential(phaml_solution,grid,procs,refine_control, &
                                io_control,loc_sequential_vert,loc_prebalance, &
                                loc_postbalance,loc_partition_method, &
                                loop,loop_end_sequential,still_sequential)
   endif

! load balance before refinement

   if (loc_task == BALANCE_REFINE_SOLVE .or. loc_task == BALANCE_ONLY .or. &
       loc_task == SET_INITIAL) then

      if (loc_prebalance /= BALANCE_NONE) then

         call balance(phaml_solution, grid, procs, io_control, refine_control, &
                      loc_partition_method, loc_prebalance, still_sequential, &
                      loop, loop_end_sequential, .true.)

      endif
   endif

! refine

   call get_grid_info(grid,procs,still_sequential,215+loop, &
                      total_nvert=nvert1,total_nelem_leaf=nelem1, &
                      total_dof=ndof1)

   if (loc_task == BALANCE_REFINE_SOLVE .or. loc_task == REFINE_ONLY .or. &
       loc_task == SET_INITIAL) then

      if (refine_control%reftype == P_ADAPTIVE .or. &
          refine_control%reftype == HP_ADAPTIVE) then
         ltemp = check_stall_special(grid,1)
      endif
      call refine(grid,procs,refine_control,solver_control,io_control, &
                  still_sequential,init_nvert,init_nelem,init_dof,loop, &
                  loc_prebalance,loc_prebalance/=BALANCE_NONE)
      call reconcile(grid,procs,refine_control,still_sequential)

      call print_grid_info(grid,procs,io_control,still_sequential, &
                           (/FREQUENTLY/),220+loop)
      call draw_grid(grid,procs,io_control,refine_control, &
                     phaml_solution%i_draw_grid, &
                     phaml_solution%master_draws_grid,still_sequential, &
                     (/FREQUENTLY/),loc_partition_method,phaml_solution%lb)

   endif

! load balance after refinement

   if (loc_task == BALANCE_REFINE_SOLVE .or. loc_task == BALANCE_ONLY .or. &
       loc_task == SET_INITIAL) then

      if (loc_postbalance /= BALANCE_NONE) then

         call balance(phaml_solution, grid, procs, io_control, refine_control, &
                      loc_partition_method, loc_postbalance, still_sequential, &
                      loop, -1, .false.)

      endif
   endif

! print and draw grid after refinement and load balancing phases

   call print_grid_info(grid,procs,io_control,still_sequential, &
                        (/PHASES/),230+loop)
   call draw_grid(grid,procs,io_control,refine_control, &
                  phaml_solution%i_draw_grid, &
                  phaml_solution%master_draws_grid,still_sequential, &
                  (/PHASES/),loc_partition_method,phaml_solution%lb)

! solve

   if (loc_task == BALANCE_REFINE_SOLVE .or. loc_task == SOLVE_ONLY .or. &
       loc_task == SET_INITIAL) then

      call solve(grid,procs,io_control,solver_control,still_sequential, &
                 loc_task==SET_INITIAL)

   endif

! update the REFSOLN error estimates, if needed

   if (refine_control%reftype == HP_ADAPTIVE) then
      if (refine_control%hp_strategy == HP_REFSOLN_EDGE) then
         call ext_refsoln_errest_edge(phaml_solution%grid,procs, &
                                    refine_control,solver_control,io_control, &
                                    phaml_solution%still_sequential)
      elseif (refine_control%hp_strategy == HP_REFSOLN_ELEM) then
         call ext_refsoln_errest_elem(phaml_solution%grid,procs, &
                                    refine_control,solver_control,io_control, &
                                    phaml_solution%still_sequential)
      endif
   endif

! print error and eigenvalues and draw grid after solution phase

   call print_error_info(grid,procs,io_control,still_sequential, &
                         (/PHASES,FREQUENTLY,TOO_MUCH/),240+loop, &
                         reduction=.true., &
                         errest=refine_control%error_estimator)

! print eigenvalues

   if (phaml_solution%eq_type == EIGENVALUE .and. &
       loc_print_eval_when == PHASES) then

      allocate(send_real(2*grid%num_eval))
      if (my_proc(procs) == MASTER) then
         if (loc_print_eval_who == MASTER .or. &
             loc_print_eval_who == EVERYONE) then
            send_real = 0.0_my_real
            do i=1,num_proc(phaml_solution%procs)
               call phaml_recv(procs,proc,recv_int,ni,recv_real,nr,111)
               if (proc==1) then
                  send_real(1:2*grid%num_eval:2)=recv_real(1:2*grid%num_eval:2)
                  send_real(2:2*grid%num_eval:2)=send_real(2:2*grid%num_eval:2)&
                                      + recv_real(2:2*grid%num_eval:2)**2
               elseif (.not. still_sequential) then
                  send_real(2:2*grid%num_eval:2)=send_real(2:2*grid%num_eval:2)&
                                      + recv_real(2:2*grid%num_eval:2)**2
               endif
               if (associated(recv_real)) deallocate(recv_real,stat=astat)
            end do
            send_real(2:2*grid%num_eval:2)=sqrt(send_real(2:2*grid%num_eval:2))
            write(outunit,"(A)")
            write(outunit,"(A)") "Eigenvalues and error estimates:"
            do i=1,size(send_real),2
               write(outunit,"(SS,1P,A,2E18.10E2)") "   ",send_real(i), &
                                                          send_real(i+1)
            end do
         endif

      else ! slave

         if (loc_print_eval_who == MASTER .or. loc_print_eval_who == EVERYONE &
                                          .or. loc_print_eval_who == SLAVES)then
            do i=0,grid%num_eval-1
               send_real(2*i+1) = grid%eigenvalue(i+1)
               call error_estimate(grid,procs,refine_control%error_estimator, &
                                   i+1,errest_eigenvalue=send_real(2*i+2))
            end do
         endif
         if (loc_print_eval_who == MASTER .or. loc_print_eval_who == EVERYONE) then
            call phaml_send(procs,MASTER,send_int,0,send_real,&
                            size(send_real),111)
         endif
         if (loc_print_eval_who == SLAVES .or. loc_print_eval_who == EVERYONE) then
            write(outunit,"(A)")
            write(outunit,"(A)") "Eigenvalues and error estimates:"
            do i=1,size(send_real),2
               write(outunit,"(SS,1P,A,2E18.10E2)") "   ",send_real(i),send_real(i+1)
            end do
         endif
         if (loc_print_eval_who == MASTER .or. loc_print_eval_who == EVERYONE &
                                          .or. loc_print_eval_who == SLAVES)then
         endif
      endif
      deallocate(send_real)
   endif

   call draw_grid(grid,procs,io_control,refine_control, &
                  phaml_solution%i_draw_grid, &
                  phaml_solution%master_draws_grid,still_sequential, &
                  (/PHASES,FREQUENTLY/),loc_partition_method, &
                  phaml_solution%lb)

! print times for this step

   call print_time_info(procs,io_control,(/PHASES/),'both',280+loop)

! pause, if requested

   if (present(pause_after_phases)) then
      call pause_until_enter(procs,pause_after_phases)
   endif

! check for termination criteria

   if (my_proc(procs) /= MASTER) then
      if (ierr == NO_ERROR .and. loc_task == BALANCE_ONLY) then
         ierr = DONE_BALANCE_ONLY
      endif
      if (ierr == NO_ERROR .and. loc_task == REFINE_ONLY) then
         ierr = DONE_REFINE_ONLY
      endif
      if (ierr == NO_ERROR .and. loc_task == SOLVE_ONLY) then
         ierr = DONE_SOLVE_ONLY
      endif
      if (ierr == NO_ERROR .and. loop >= loc_max_refsolveloop) then
         ierr = MAX_LOOP_ACHIEVED
      endif
      call get_grid_info(grid,procs,still_sequential,250+loop, &
                         total_nvert=total_vert,total_nelem_leaf=total_elem, &
                         max_nlev=total_nlev,total_dof=total_dof, &
                         no_master=.true.)
      if (ierr == NO_ERROR .and. total_vert >= refine_control%max_vert) then
         ierr = MAX_VERTICES_ACHIEVED
      endif
      if (ierr == NO_ERROR .and. total_elem >= refine_control%max_elem) then
         ierr = MAX_ELEMENTS_ACHIEVED
      endif
      if (ierr == NO_ERROR .and. total_dof >= refine_control%max_dof) then
         ierr = MAX_EQUATIONS_ACHIEVED
      endif
      if (ierr == NO_ERROR .and. total_nlev >= refine_control%max_lev .and. &
          loc_stop_on_maxlev) then
         ierr = MAX_LEV_ACHIEVED
      endif
      if (loc_stop_on_maxdeg) then
         call get_grid_info(grid,procs,still_sequential,260+loop, &
                            maxdeg=total_deg)
         if (ierr == NO_ERROR .and. total_deg >= refine_control%max_deg) then
            ierr = MAX_DEG_ACHIEVED
         endif
      endif
! error estimate termination is satisfied if all components of all eigenvectors
! are estimated smaller than the tolerance
      if (ierr == NO_ERROR .and. refine_control%term_energy_err > 0.0_my_real) then
         soln = 0
         ierr = ENERGY_ERREST_ACHIEVED
         outer_energy: do comp=1,1
            do evec=1,max(1,grid%num_eval)
               soln = soln+1
               call error_estimate(grid,procs,refine_control%error_estimator, &
                                   soln,errest_energy=errest_energy)
               errest_energy = sqrt(phaml_global_sum(procs,errest_energy**2, &
                                    4100+soln+loop))
               if (grid%errtype == RELATIVE_ERROR .and. &
                   .not.(refine_control%reftype==HP_ADAPTIVE .and. &
                         (refine_control%hp_strategy==HP_REFSOLN_EDGE .or. &
                          refine_control%hp_strategy==HP_REFSOLN_ELEM))) then
                  call norm_solution(grid,procs,still_sequential,comp,evec, &
                                     energy=normsoln)
                  normsoln = sqrt(phaml_global_sum(procs,normsoln**2, &
                                  4600+soln+loop))
                  if (normsoln /= 0.0_my_real) then
                     errest_energy = errest_energy/normsoln
                  endif
               endif
               if (errest_energy > refine_control%term_energy_err) then
                  ierr = NO_ERROR
                  exit outer_energy
               endif
            end do
         end do outer_energy
      endif
      if (ierr==NO_ERROR .and. refine_control%term_Linf_err>0.0_my_real .and. &
          (refine_control%reftype/=HP_ADAPTIVE .or. &
           (refine_control%hp_strategy/=HP_REFSOLN_EDGE .and. &
            refine_control%hp_strategy/=HP_REFSOLN_ELEM))) then
         soln = 0
         ierr = LINF_ERREST_ACHIEVED
         outer_Linf: do comp=1,grid%system_size
            do evec=1,max(1,grid%num_eval)
               soln = soln+1
               call error_estimate(grid,procs,refine_control%error_estimator, &
                                   soln,errest_Linf=errest_Linf)
               errest_Linf = phaml_global_max(procs,errest_Linf,4200+soln+loop)
               if (grid%errtype == RELATIVE_ERROR) then
                  call norm_solution(grid,procs,still_sequential,comp,evec, &
                                     linf=normsoln)
                  normsoln = phaml_global_max(procs,normsoln,4700+soln+loop)
                  if (normsoln /= 0.0_my_real) then
                     errest_Linf = errest_Linf/normsoln
                  endif
               endif
               if (errest_Linf > refine_control%term_Linf_err) then
                  ierr = NO_ERROR
                  exit outer_Linf
               endif
            end do
         end do outer_Linf
      endif
      if (ierr==NO_ERROR .and. refine_control%term_L2_err>0.0_my_real .and. &
          (refine_control%reftype/=HP_ADAPTIVE .or. &
           (refine_control%hp_strategy/=HP_REFSOLN_EDGE .and. &
            refine_control%hp_strategy/=HP_REFSOLN_ELEM))) then
         soln = 0
         ierr = L2_ERREST_ACHIEVED
         outer_L2: do comp=1,grid%system_size
            do evec=1,max(1,grid%num_eval)
               soln = soln+1
               call error_estimate(grid,procs,refine_control%error_estimator, &
                                   soln,errest_L2=errest_L2)
               errest_L2 = sqrt(phaml_global_sum(procs,errest_L2**2, &
                                4300+soln+loop))
               if (grid%errtype == RELATIVE_ERROR) then
                  call norm_solution(grid,procs,still_sequential,comp,evec, &
                                     L2=normsoln)
                  normsoln = sqrt(phaml_global_sum(procs,normsoln**2, &
                                  4800+soln+loop))
                  if (normsoln /= 0.0_my_real) then
                     errest_L2 = errest_L2/normsoln
                  endif
               endif
               if (errest_L2 > refine_control%term_L2_err) then
                  ierr = NO_ERROR
                  exit outer_L2
               endif
            end do
         end do outer_L2
      endif
      call get_grid_info(grid,procs,still_sequential,265+loop, &
                         total_nvert=nvert2,total_nelem_leaf=nelem2, &
                         total_dof=ndof2,no_master=.true.)
      if (ierr == NO_ERROR .and. nvert1==nvert2 .and. nelem1==nelem2 .and. ndof1==ndof2) then
         if (refine_control%reftype == P_ADAPTIVE .or. &
             refine_control%reftype == HP_ADAPTIVE) then
            if (check_stall_special(grid,2)) ierr = REFINEMENT_STALLED
         endif
      endif
   endif

! coordinate ierr using global max; processor 1 sends it to the master.
! To avoid slaves running away from the master, first send a message from
! the master to processor 1 (the global_max will cause the others to wait)

   if (my_proc(procs) /= MASTER) then
      if (my_proc(procs) == 1) then
         call phaml_recv(procs,proc,recv_int,ni,recv_real,nr,370+loop)
         if (associated(recv_int)) then
            ierr1 = recv_int(1)
            deallocate(recv_int,stat=astat)
            ierr2 = ierr
            if (ierr1 == NO_ERROR .or. ierr2 == NO_ERROR .or. ierr1 == ierr2) then
               if (ierr1 == NO_ERROR) then
                  ierr = ierr2
               else
                  ierr = ierr1
               endif
            else
               ierr = MULTIPLE_ERRORS
            endif
         endif
      endif
      ierr1 = phaml_global_min(procs,ierr,260+loop)
      ierr2 = phaml_global_max(procs,ierr,365+loop)
      if (ierr1 == NO_ERROR .or. ierr2 == NO_ERROR .or. ierr1 == ierr2) then
         if (ierr1 == NO_ERROR) then
            ierr = ierr2
         else
            ierr = ierr1
         endif
      else
         ierr = MULTIPLE_ERRORS
      endif
      if (my_proc(procs) == 1) then
         send_int(1) = ierr
         call phaml_send(procs,MASTER,send_int,1,no_reals,0,270+loop)
      endif
   else
      send_int(1) = ierr
      call phaml_send(procs,1,send_int,1,no_reals,0,370+loop)
      call phaml_recv(procs,proc,recv_int,ni,recv_real,nr,270+loop)
      if (associated(recv_int)) then
         ierr = recv_int(1)
         deallocate(recv_int,stat=astat)
      endif
   endif

end do mainloop ! next refine/partition/solve step

call stop_watch(ptotal)

! print final grid, error, eigenvalues, time and communication information,
! and graphics

call print_grid_info(grid,procs,io_control,still_sequential, &
                     (/FREQUENTLY,PHASES,FINAL/),120)
call print_error_info(grid,procs,io_control,still_sequential, &
                      (/FREQUENTLY,PHASES,FINAL/),121, &
                      errest=refine_control%error_estimator)
if (phaml_solution%eq_type == EIGENVALUE .and. &
    loc_print_eval_when == FINAL) then

   allocate(send_real(2*grid%num_eval))
   if (my_proc(procs) == MASTER) then
      if (loc_print_eval_who == MASTER .or. &
          loc_print_eval_who == EVERYONE) then
         send_real = 0.0_my_real
         do i=1,num_proc(phaml_solution%procs)
            call phaml_recv(procs,proc,recv_int,ni,recv_real,nr,112)
            if (proc==1) then
               send_real(1:2*grid%num_eval:2)=recv_real(1:2*grid%num_eval:2)
               send_real(2:2*grid%num_eval:2)=send_real(2:2*grid%num_eval:2)&
                                   + recv_real(2:2*grid%num_eval:2)**2
            elseif (.not. still_sequential) then
               send_real(2:2*grid%num_eval:2)=send_real(2:2*grid%num_eval:2)&
                                   + recv_real(2:2*grid%num_eval:2)**2
            endif
            if (associated(recv_real)) deallocate(recv_real,stat=astat)
         end do
         send_real(2:2*grid%num_eval:2)=sqrt(send_real(2:2*grid%num_eval:2))
         write(outunit,"(A)")
         write(outunit,"(A)") "Eigenvalues and error estimates:"
         do i=1,size(recv_real),2
            write(outunit,"(SS,1P,A,2E18.10E2)") "   ",send_real(i), &
                                                       send_real(i+1)
         end do
      endif

   else ! slave

      if (loc_print_eval_who == MASTER .or. loc_print_eval_who == EVERYONE &
                                       .or. loc_print_eval_who == SLAVES)then
         do i=0,grid%num_eval-1
            send_real(2*i+1) = grid%eigenvalue(i+1)
            call error_estimate(grid,procs,refine_control%error_estimator, &
                                i+1,errest_eigenvalue=send_real(2*i+2))
         end do
      endif
      if (loc_print_eval_who == MASTER .or. loc_print_eval_who == EVERYONE) then
         call phaml_send(procs,MASTER,send_int,0,send_real,&
                         size(send_real),112)
      endif
      if (loc_print_eval_who == SLAVES .or. loc_print_eval_who == EVERYONE) then
         write(outunit,"(A)")
         write(outunit,"(A)") "Eigenvalues and error estimates:"
         do i=1,size(send_real),2
            write(outunit,"(SS,1P,A,2E18.10E2)") "   ",send_real(i),send_real(i+1)
         end do
      endif
      if (loc_print_eval_who == MASTER .or. loc_print_eval_who == EVERYONE &
                                       .or. loc_print_eval_who == SLAVES)then
      endif
   endif
   deallocate(send_real)
endif
call draw_grid(grid,procs,io_control,refine_control,phaml_solution%i_draw_grid,&
               phaml_solution%master_draws_grid,still_sequential,(/FINAL/), &
               loc_partition_method,phaml_solution%lb)
select case(io_control%print_time_when)
   case(PHASES,FINAL)
      which_time = 'total'
   case(LAST)
      which_time = 'cycle'
   case(LAST_AND_FINAL)
      which_time = 'both'
end select
call print_time_info(procs,io_control,(/PHASES,FINAL,LAST,LAST_AND_FINAL/), &
                     which_time,122)
call print_trailer
call terminate(phaml_solution,procs,solver_control,still_sequential, &
               pause_at_end)
if (present(iterm)) iterm = ierr

return

contains ! internal routines for phaml_solve_pde

!          --------
subroutine defaults
!          --------

!----------------------------------------------------
! This contains the assignment of default values to all the optional
! arguments to subroutine phaml_solve_pde
!----------------------------------------------------

use hash_mod
logical :: doit

! print warnings

if (present(print_warnings)) then
   warn = print_warnings
else
   warn = .true.
endif

! task

if (present(task)) then
   loc_task = task
else
   loc_task = BALANCE_REFINE_SOLVE
endif

if (present(solve_init)) then
   loc_solve_init = solve_init
else
   loc_solve_init = .true.
endif

! refinement control

if (present(reftype)) then
   refine_control%reftype = reftype
else
   refine_control%reftype = H_ADAPTIVE
endif

if (present(refterm)) then
   refine_control%refterm = refterm
else
   refine_control%refterm = DOUBLE_NEQ_SMOOTH
endif

if (present(hp_strategy)) then
   refine_control%hp_strategy = hp_strategy
else
   refine_control%hp_strategy = HP_PRIOR2P_H1
endif

if (present(inc_factor)) then
   refine_control%inc_factor = inc_factor
else
   refine_control%inc_factor = 2.0_my_real
endif

if (refine_control%reftype == HP_ADAPTIVE .and. &
    (refine_control%hp_strategy == HP_REFSOLN_EDGE .or. &
     refine_control%hp_strategy == HP_REFSOLN_ELEM)) then
   refine_control%edge_rule = MINIMUM_RULE
elseif (present(edge_rule)) then
   refine_control%edge_rule = edge_rule
else
   refine_control%edge_rule = MINIMUM_RULE
endif

if (present(t3s_gamma)) then
   refine_control%t3s_gamma = t3s_gamma
else
   refine_control%t3s_gamma = 6.0_my_real
endif

if (present(t3s_eta)) then
   refine_control%t3s_eta = t3s_eta
else
   refine_control%t3s_eta = 0.1_my_real
endif

if (present(t3s_nunif)) then
   refine_control%t3s_nunif = t3s_nunif
else
   refine_control%t3s_nunif = 0
endif

if (present(t3s_maxref)) then
   refine_control%t3s_maxref = t3s_maxref
else
   refine_control%t3s_maxref = 3
endif

if (present(t3s_maxdeginc)) then
   refine_control%t3s_maxdeginc = t3s_maxdeginc
else
   refine_control%t3s_maxdeginc = 3
endif

if (present(tp_gamma)) then
   refine_control%tp_gamma = tp_gamma
else
   refine_control%tp_gamma = 0.2_my_real
endif

if (present(sp_gamma_h)) then
   refine_control%sp_gamma_h = sp_gamma_h
else
   refine_control%sp_gamma_h = 4.0_my_real
endif

if (present(sp_gamma_p)) then
   refine_control%sp_gamma_p = sp_gamma_p
else
   refine_control%sp_gamma_p = 0.4_my_real
endif

if (present(nlp_max_h_dec)) then
   refine_control%nlp_max_h_dec = nlp_max_h_dec
else
   refine_control%nlp_max_h_dec = 1
endif

if (present(nlp_max_h_inc)) then
   refine_control%nlp_max_h_inc = nlp_max_h_inc
else
   refine_control%nlp_max_h_inc = 1
endif

if (present(nlp_max_p_dec)) then
   refine_control%nlp_max_p_dec = nlp_max_p_dec
else
   refine_control%nlp_max_p_dec = 5
endif

if (present(nlp_max_p_inc)) then
   refine_control%nlp_max_p_inc = nlp_max_p_inc
else
   refine_control%nlp_max_p_inc = 2
endif

if (present(refsoln_pbias)) then
   refine_control%refsoln_pbias = refsoln_pbias
else
   refine_control%refsoln_pbias = 4.0_my_real
endif

if (present(derefine)) then
   refine_control%derefine = derefine
else
   refine_control%derefine = .true.
endif

if (refine_control%reftype == HP_ADAPTIVE .and. &
    refine_control%hp_strategy == HP_T3S .and. &
    refine_control%derefine) then
   call warning("setting derefine to false because hp strategy is HP_T3S")
   refine_control%derefine = .false.
endif

if (refine_control%reftype == HP_ADAPTIVE .and. &
    (refine_control%hp_strategy == HP_REFSOLN_EDGE .or. &
     refine_control%hp_strategy == HP_REFSOLN_ELEM)) then
   refine_control%error_estimator = REFSOLN_ERREST
elseif (present(error_estimator)) then
   refine_control%error_estimator = error_estimator
else
   refine_control%error_estimator = EXPLICIT_ERRIND
endif

if (loc_task == SET_INITIAL .and. &
    refine_control%error_estimator /= INITIAL_CONDITION) then
   call warning("error_estimator changed to INITIAL_CONDITION for task==SET_INITIAL")
   refine_control%error_estimator = INITIAL_CONDITION
endif

! load balancing control

if (present(sequential_vert)) then
   loc_sequential_vert = sequential_vert
else
   loc_sequential_vert = 100
endif

if (present(partition_method)) then
   loc_partition_method = partition_method
else
   loc_partition_method = RTK
endif

if (present(zoltan_param_file)) then
   if (len(zoltan_param_file) > FN_LEN) then
      call warning("Name of Zoltan param file is too long and will be truncated.",&
                   "Maximum length (FN_LEN in global.f90) is", &
                   intlist=(/FN_LEN/))
   endif
   phaml_solution%zoltan_param_file = zoltan_param_file
else
   phaml_solution%zoltan_param_file = "zoltan.params"
endif

if (present(prebalance)) then
   loc_prebalance = prebalance
else
   loc_prebalance = BALANCE_ELEMENTS
endif

if (present(postbalance)) then
   loc_postbalance = postbalance
else
   loc_postbalance = BALANCE_NONE
endif

! polynomial degree

if (present(degree)) then
   loc_degree = degree
else
   loc_degree = 0
endif

! solver choices

if (present(solver)) then
   solver_control%solver = solver
else
   solver_control%solver = MG_SOLVER
endif

if (present(preconditioner)) then
   solver_control%preconditioner = preconditioner
else
   select case(solver_control%solver)
   case (MG_SOLVER,LAPACK_INDEFINITE_SOLVER,LAPACK_SPD_SOLVER, &
         MUMPS_SPD_SOLVER, MUMPS_GEN_SOLVER, MUMPS_NONSYM_SOLVER, &
         SUPERLU_SOLVER, HYPRE_BOOMERAMG_SOLVER)
      solver_control%preconditioner = NO_PRECONDITION
   case (CG_SOLVER, GMRES_SOLVER)
      solver_control%preconditioner = MG_PRECONDITION
   case(PETSC_RICHARDSON_SOLVER, PETSC_CHEBYCHEV_SOLVER, PETSC_CG_SOLVER, &
        PETSC_GMRES_SOLVER, PETSC_TCQMR_SOLVER, PETSC_BCGS_SOLVER, &
        PETSC_CGS_SOLVER, PETSC_TFQMR_SOLVER, PETSC_CR_SOLVER, &
        PETSC_LSQR_SOLVER, PETSC_BICG_SOLVER)
      solver_control%preconditioner = MG_PRECONDITION
   case(HYPRE_GMRES_SOLVER, HYPRE_PCG_SOLVER)
      solver_control%preconditioner = HYPRE_BOOMERAMG_PRECONDITION
   case default
      call warning("illegal value for choice of solver")
      solver_control%preconditioner = NO_PRECONDITION
   end select
endif

if (present(inc_quad_order)) then
   solver_control%inc_quad_order = inc_quad_order
else
   solver_control%inc_quad_order = 0
endif

if (present(mg_tol)) then
   solver_control%mg_tol = mg_tol
else
   if (solver_control%solver == MG_SOLVER) then
      solver_control%mg_tol = MG_ERREST_TOL
   else
      solver_control%mg_tol = MG_NO_TOL
   endif
endif

if (present(mg_cycles)) then
   solver_control%ncycle = mg_cycles
else
   if (solver_control%solver == MG_SOLVER) then
      if (solver_control%mg_tol == MG_NO_TOL) then
         solver_control%ncycle = 1
      else
         solver_control%ncycle = huge(0)
      endif
   else
      solver_control%ncycle = 2
   endif
endif

if (present(mg_prerelax)) then
   solver_control%prerelax = mg_prerelax
else
   solver_control%prerelax = 1 ! just red relaxation
endif

if (present(mg_postrelax)) then
   solver_control%postrelax = mg_postrelax
else
   solver_control%postrelax = 2 ! red relaxation followed by local black
endif

if (present(mg_prerelax_ho)) then
   solver_control%prerelax_ho = mg_prerelax_ho
else
   solver_control%prerelax_ho = 1
endif

if (present(mg_postrelax_ho)) then
   solver_control%postrelax_ho = mg_postrelax_ho
else
   solver_control%postrelax_ho = 1
endif

if (present(mg_comm)) then
   solver_control%mg_comm = mg_comm
else
   if (solver_control%solver == MG_SOLVER) then
      solver_control%mg_comm = MGCOMM_FUDOP
   else
      solver_control%mg_comm = MGCOMM_NONE
   endif
endif

if (present(krylov_iter)) then
   solver_control%krylov_iter = krylov_iter
else
   solver_control%krylov_iter = 100
endif

if (present(krylov_restart)) then
   solver_control%krylov_restart = krylov_restart
else
   solver_control%krylov_restart = 20
endif

if (present(krylov_tol)) then
   solver_control%krylov_tol = krylov_tol
else
   solver_control%krylov_tol = KRYLOV_ERREST_TOL
endif

if (present(dd_iterations)) then
   solver_control%dd_iterations = dd_iterations
else
   solver_control%dd_iterations = 1
endif

if (present(petsc_matrix_free)) then
   solver_control%petsc_matrix_free = petsc_matrix_free
else
   solver_control%petsc_matrix_free = .false.
endif

if (present(ignore_quad_err)) then
   solver_control%ignore_quad_err = ignore_quad_err
else
   if (solver_control%solver == MG_SOLVER) then
      solver_control%ignore_quad_err = .false.
   else
      solver_control%ignore_quad_err = .true.
   endif
endif

if (present(coarse_size)) then
   solver_control%coarse_size = coarse_size
else
   solver_control%coarse_size = 5000
endif

if (present(coarse_method)) then
   solver_control%coarse_method = coarse_method
else
   solver_control%coarse_method = LAPACK_INDEFINITE_SOLVER
endif

if (present(eigensolver)) then
   solver_control%eigensolver = eigensolver
else
   solver_control%eigensolver = ARPACK_SOLVER
endif

if (present(num_eval)) then
   solver_control%num_eval = num_eval
else
   solver_control%num_eval = 1
endif

if (present(lambda0)) then
   solver_control%lambda0 = lambda0
else
   solver_control%lambda0 = -huge(0.0_my_real)
endif

if (present(scale_evec)) then
   solver_control%scale_evec = scale_evec
else
   solver_control%scale_evec = SCALE_LINF
endif

if (present(transformation)) then
   solver_control%transformation = transformation
else
   solver_control%transformation = SHIFT_INVERT
endif

if (present(lambda0_side)) then
   solver_control%lambda0_side = lambda0_side
else
   if (solver_control%eigensolver == ARPACK_SOLVER) then
      solver_control%lambda0_side = EIGEN_BOTH
   else ! BLOPEX
      if (solver_control%transformation == SHIFT_INVERT) then
         solver_control%lambda0_side = EIGEN_RIGHT
      else ! SHIFT_SQUARE
         solver_control%lambda0_side = EIGEN_BOTH
      endif
   endif
endif

if (present(arpack_ncv)) then
   solver_control%arpack_cntl%ncv = arpack_ncv
else
   solver_control%arpack_cntl%ncv = 20
endif

if (present(arpack_maxit)) then
   solver_control%arpack_cntl%maxit = arpack_maxit
else
   solver_control%arpack_cntl%maxit = 100
endif

if (present(arpack_tol)) then
   solver_control%arpack_cntl%tol = arpack_tol
else
   solver_control%arpack_cntl%tol = 1.0e-10_my_real
endif

if (present(blopex_maxit)) then
   solver_control%blopex_cntl%maxit = blopex_maxit
else
   solver_control%blopex_cntl%maxit = 100
endif

if (present(blopex_atol)) then
   solver_control%blopex_cntl%atol = blopex_atol
else
   solver_control%blopex_cntl%atol = 1.0e-6_my_real
endif

if (present(blopex_rtol)) then
   solver_control%blopex_cntl%rtol = blopex_rtol
else
   solver_control%blopex_cntl%rtol = 1.0e-6_my_real
endif

! termination criteria

if (present(max_elem)) then
   refine_control%max_elem = max_elem
else
   refine_control%max_elem = huge(0)
endif

if (present(max_vert)) then
   refine_control%max_vert = max_vert
else
   refine_control%max_vert = huge(0)
endif

if (present(max_eq)) then
   refine_control%max_dof = max_eq
else
   refine_control%max_dof = huge(0)
endif

if (present(max_refsolveloop)) then
   loc_max_refsolveloop = max_refsolveloop
else
   loc_max_refsolveloop = huge(0)
endif

if (present(max_lev)) then
   refine_control%max_lev = max_lev
else
   refine_control%max_lev = huge(0)
endif

if (present(max_deg)) then
   refine_control%max_deg = max_deg
else
   if (refine_control%reftype == HP_ADAPTIVE .and. &
       refine_control%hp_strategy == HP_NEXT3P) then
      refine_control%max_deg = 19
   elseif (refine_control%error_estimator == LOCAL_PROBLEM_P .or. &
           refine_control%error_estimator == EQUILIBRATED_RESIDUAL .or. &
           (refine_control%reftype == HP_ADAPTIVE .and. &
            (refine_control%hp_strategy == HP_REFSOLN_EDGE .or. &
             refine_control%hp_strategy == HP_REFSOLN_ELEM))) then
      refine_control%max_deg = 21
   else
      refine_control%max_deg = 22
   endif
endif

if (present(stop_on_maxlev)) then
   loc_stop_on_maxlev = stop_on_maxlev
else
   loc_stop_on_maxlev = .false.
endif

if (present(stop_on_maxdeg)) then
   loc_stop_on_maxdeg = stop_on_maxdeg
else
   loc_stop_on_maxdeg = .false.
endif

if (present(term_energy_err)) then
   refine_control%term_energy_err = term_energy_err
else
   refine_control%term_energy_err = 0.0_my_real
endif

if (present(term_Linf_err)) then
   refine_control%term_Linf_err = term_Linf_err
else
   refine_control%term_Linf_err = 0.0_my_real
endif

if (present(term_L2_err)) then
   refine_control%term_L2_err = term_L2_err
else
   refine_control%term_L2_err = 0.0_my_real
endif

if (present(reftol)) then
   refine_control%reftol = reftol
else
   refine_control%reftol = refine_control%term_energy_err/2.0_my_real
endif

! printed output control

if (present(print_grid_when)) then
   io_control%print_grid_when = print_grid_when
else
   io_control%print_grid_when = NEVER
endif

if (present(print_grid_who)) then
   io_control%print_grid_who = print_grid_who
   if (PARALLEL==SEQUENTIAL .and. &
       (print_grid_who == MASTER .or. print_grid_who == MASTER_ALL)) then
      io_control%print_grid_who = SLAVES
   endif
else
   if (PARALLEL/=SEQUENTIAL) then
      io_control%print_grid_who = MASTER
   else
      io_control%print_grid_who = SLAVES
   endif
endif

if (present(print_linsys_when)) then
   io_control%print_linsys_when = print_linsys_when
else
   io_control%print_linsys_when = NEVER
endif

if (present(print_linsys_who)) then
   io_control%print_linsys_who = print_linsys_who
   if (PARALLEL==SEQUENTIAL .and. &
       (print_linsys_who == MASTER .or. print_linsys_who == MASTER_ALL)) then
      io_control%print_linsys_who = SLAVES
   endif
else
   if (PARALLEL/=SEQUENTIAL) then
      io_control%print_linsys_who = MASTER
   else
      io_control%print_linsys_who = SLAVES
   endif
endif

if (present(print_error_when)) then
   io_control%print_error_when = print_error_when
else
   io_control%print_error_when = NEVER
endif

if (present(print_error_who)) then
   io_control%print_error_who = print_error_who
   if (PARALLEL==SEQUENTIAL .and. &
       (print_error_who == MASTER .or. print_error_who == MASTER_ALL)) then
      io_control%print_error_who = SLAVES
   endif
else
   if (PARALLEL/=SEQUENTIAL) then
      io_control%print_error_who = MASTER
   else
      io_control%print_error_who = SLAVES
   endif
endif

if (present(print_error_what)) then
   io_control%print_error_what = print_error_what
else
   io_control%print_error_what = NEVER
endif

if (present(print_errest_what)) then
   io_control%print_errest_what = print_errest_what
else
   io_control%print_errest_what = NEVER
endif

if (present(print_time_when)) then
   io_control%print_time_when = print_time_when
else
   io_control%print_time_when = NEVER
endif

if (present(print_time_who)) then
   io_control%print_time_who = print_time_who
   if (PARALLEL==SEQUENTIAL .and. &
       (print_time_who == MASTER .or. print_time_who == MASTER_ALL)) then
      io_control%print_time_who = SLAVES
   endif
else
   if (PARALLEL/=SEQUENTIAL) then
      io_control%print_time_who = MASTER
   else
      io_control%print_time_who = SLAVES
   endif
endif

if (present(print_eval_when)) then
   loc_print_eval_when = print_eval_when
else
   loc_print_eval_when = NEVER
endif

if (present(print_eval_who)) then
   loc_print_eval_who = print_eval_who
   if (PARALLEL==SEQUENTIAL .and. &
       (print_eval_who == MASTER .or. print_eval_who == MASTER_ALL)) then
      loc_print_eval_who = SLAVES
   endif
else
   if (PARALLEL/=SEQUENTIAL) then
      loc_print_eval_who = MASTER
   else
      loc_print_eval_who = SLAVES
   endif
endif

if (present(print_header_who)) then
   loc_print_header_who = print_header_who
   if (print_header_who == MASTER_ALL) loc_print_header_who = MASTER
   if (PARALLEL==SEQUENTIAL .and. loc_print_header_who == MASTER) then
      loc_print_header_who = SLAVES
   endif
else
   if (PARALLEL/=SEQUENTIAL) then
      loc_print_header_who = MASTER
   else
      loc_print_header_who = SLAVES
   endif
endif

if (present(print_trailer_who)) then
   loc_print_trailer_who = print_trailer_who
   if (PARALLEL==SEQUENTIAL .and. loc_print_trailer_who == MASTER) then
      loc_print_trailer_who = SLAVES
   endif
else
   if (PARALLEL/=SEQUENTIAL) then
      loc_print_trailer_who = MASTER
   else
      loc_print_trailer_who = SLAVES
   endif
endif

if (present(errtype)) then
   phaml_solution%grid%errtype = errtype
else
   phaml_solution%grid%errtype = ABSOLUTE_ERROR
endif

if (present(clocks)) then
   loc_clocks = clocks
else
   loc_clocks = CLOCK_W
endif

select case (loc_clocks)
case(CLOCK_W); call option_stopwatch('wall')
case(CLOCK_C); call option_stopwatch('cpu')
case(CLOCK_CW); call option_stopwatch((/'cpu ','wall'/))
end select

! graphical output control

if (present(draw_grid_when)) then
   io_control%draw_grid_when = draw_grid_when
else
   io_control%draw_grid_when = NEVER
endif

if (present(pause_after_draw)) then
   io_control%pause_after_draw = pause_after_draw
else
   io_control%pause_after_draw = .false.
endif

! hypre options

if (present(hypre_BoomerAMG_MaxLevels)) then
   solver_control%hypre_cntl%BoomerAMG_MaxLevels = hypre_BoomerAMG_MaxLevels
else
   solver_control%hypre_cntl%BoomerAMG_MaxLevels = huge(0)
endif

if (present(hypre_BoomerAMG_MaxIter)) then
   solver_control%hypre_cntl%BoomerAMG_MaxIter = hypre_BoomerAMG_MaxIter
else
   solver_control%hypre_cntl%BoomerAMG_MaxIter = huge(0)
endif

if (present(hypre_BoomerAMG_Tol)) then
   solver_control%hypre_cntl%BoomerAMG_Tol = hypre_BoomerAMG_Tol
else
   solver_control%hypre_cntl%BoomerAMG_Tol = huge(0.0d0)
endif

if (present(hypre_BoomerAMG_StrongThreshold)) then
   solver_control%hypre_cntl%BoomerAMG_StrongThreshold = hypre_BoomerAMG_StrongThreshold
else
   solver_control%hypre_cntl%BoomerAMG_StrongThreshold = huge(0.0d0)
endif

if (present(hypre_BoomerAMG_MaxRowSum)) then
   solver_control%hypre_cntl%BoomerAMG_MaxRowSum = hypre_BoomerAMG_MaxRowSum
else
   solver_control%hypre_cntl%BoomerAMG_MaxRowSum = huge(0.0d0)
endif

if (present(hypre_BoomerAMG_CoarsenType)) then
   solver_control%hypre_cntl%BoomerAMG_CoarsenType = hypre_BoomerAMG_CoarsenType
else
   solver_control%hypre_cntl%BoomerAMG_CoarsenType = huge(0)
endif

if (present(hypre_BoomerAMG_MeasureType)) then
   solver_control%hypre_cntl%BoomerAMG_MeasureType = hypre_BoomerAMG_MeasureType
else
   solver_control%hypre_cntl%BoomerAMG_MeasureType = huge(0)
endif

if (present(hypre_BoomerAMG_CycleType)) then
   solver_control%hypre_cntl%BoomerAMG_CycleType = hypre_BoomerAMG_CycleType
else
   solver_control%hypre_cntl%BoomerAMG_CycleType = huge(0)
endif

if (present(hypre_BoomerAMG_NumGridSweeps)) then
   if (present(hypre_has_NumGridSweeps)) then
      doit = hypre_has_NumGridSweeps
   else
      doit = .true.
   endif
else
   doit = .false.
endif
if (doit) then
   allocate(solver_control%hypre_cntl%BoomerAMG_NumGridSweeps(size(hypre_BoomerAMG_NumGridSweeps)))
   solver_control%hypre_cntl%BoomerAMG_NumGridSweeps = hypre_BoomerAMG_NumGridSweeps
else
   nullify(solver_control%hypre_cntl%BoomerAMG_NumGridSweeps)
endif

if (present(hypre_BoomerAMG_GridRelaxType)) then
   if (present(hypre_has_GridRelaxType)) then
      doit = hypre_has_GridRelaxType
   else
      doit = .true.
   endif
else
   doit = .false.
endif
if (doit) then
   allocate(solver_control%hypre_cntl%BoomerAMG_GridRelaxType(size(hypre_BoomerAMG_GridRelaxType)))
   solver_control%hypre_cntl%BoomerAMG_GridRelaxType = hypre_BoomerAMG_GridRelaxType
else
   nullify(solver_control%hypre_cntl%BoomerAMG_GridRelaxType)
endif

if (present(hypre_BoomerAMG_GridRelaxPoints)) then
   if (present(hypre_has_GridRelaxPoints)) then
      doit = hypre_has_GridRelaxPoints
   else
      doit = .true.
   endif
else
   doit = .false.
endif
if (doit) then
   allocate(solver_control%hypre_cntl%BoomerAMG_GridRelaxPoints( &
            size(hypre_BoomerAMG_GridRelaxPoints,1), &
            size(hypre_BoomerAMG_GridRelaxPoints,2)))
   solver_control%hypre_cntl%BoomerAMG_GridRelaxPoints = hypre_BoomerAMG_GridRelaxPoints
else
   nullify(solver_control%hypre_cntl%BoomerAMG_GridRelaxPoints)
endif

if (present(hypre_BoomerAMG_RelaxWeight)) then
   if (present(hypre_has_RelaxWeight)) then
      doit = hypre_has_RelaxWeight
   else
      doit = .true.
   endif
else
   doit = .false.
endif
if (doit) then
   allocate(solver_control%hypre_cntl%BoomerAMG_RelaxWeight(size(hypre_BoomerAMG_RelaxWeight)))
   solver_control%hypre_cntl%BoomerAMG_RelaxWeight = hypre_BoomerAMG_RelaxWeight
else
   nullify(solver_control%hypre_cntl%BoomerAMG_RelaxWeight)
endif

if (present(hypre_BoomerAMG_DebugFlag)) then
   solver_control%hypre_cntl%BoomerAMG_DebugFlag = hypre_BoomerAMG_DebugFlag
else
   solver_control%hypre_cntl%BoomerAMG_DebugFlag = huge(0)
endif

if (present(hypre_ParaSails_thresh)) then
   solver_control%hypre_cntl%ParaSails_thresh = hypre_ParaSails_thresh
else
   solver_control%hypre_cntl%ParaSails_thresh = huge(0.0d0)
endif

if (present(hypre_ParaSails_nlevels)) then
   solver_control%hypre_cntl%ParaSails_nlevels = hypre_ParaSails_nlevels
else
   solver_control%hypre_cntl%ParaSails_nlevels = huge(0)
endif

if (present(hypre_ParaSails_filter)) then
   solver_control%hypre_cntl%ParaSails_filter = hypre_ParaSails_filter
else
   solver_control%hypre_cntl%ParaSails_filter = huge(0.0d0)
endif

if (present(hypre_ParaSails_sym)) then
   solver_control%hypre_cntl%ParaSails_sym = hypre_ParaSails_sym
else
   solver_control%hypre_cntl%ParaSails_sym = huge(0)
endif

if (present(hypre_ParaSails_loadbal)) then
   solver_control%hypre_cntl%ParaSails_loadbal = hypre_ParaSails_loadbal
else
   solver_control%hypre_cntl%ParaSails_loadbal = huge(0.0d0)
endif

if (present(hypre_ParaSails_reuse)) then
   solver_control%hypre_cntl%ParaSails_reuse = hypre_ParaSails_reuse
else
   solver_control%hypre_cntl%ParaSails_reuse = huge(0)
endif

if (present(hypre_ParaSails_logging)) then
   solver_control%hypre_cntl%ParaSails_logging = hypre_ParaSails_logging
else
   solver_control%hypre_cntl%ParaSails_logging = huge(0)
endif

if (present(hypre_PCG_Tol)) then
   solver_control%hypre_cntl%PCG_Tol = hypre_PCG_Tol
else
   solver_control%hypre_cntl%PCG_Tol = huge(0.0d0)
endif

if (present(hypre_PCG_MaxIter)) then
   solver_control%hypre_cntl%PCG_MaxIter = hypre_PCG_MaxIter
else
   solver_control%hypre_cntl%PCG_MaxIter = huge(0)
endif

if (present(hypre_PCG_TwoNorm)) then
   solver_control%hypre_cntl%PCG_TwoNorm = hypre_PCG_TwoNorm
else
   solver_control%hypre_cntl%PCG_TwoNorm = huge(0)
endif

if (present(hypre_PCG_RelChange)) then
   solver_control%hypre_cntl%PCG_RelChange = hypre_PCG_RelChange
else
   solver_control%hypre_cntl%PCG_RelChange = huge(0)
endif

if (present(hypre_PCG_Logging)) then
   solver_control%hypre_cntl%PCG_Logging = hypre_PCG_Logging
else
   solver_control%hypre_cntl%PCG_Logging = huge(0)
endif

if (present(hypre_GMRES_KDim)) then
   solver_control%hypre_cntl%GMRES_KDim = hypre_GMRES_KDim
else
   solver_control%hypre_cntl%GMRES_KDim = huge(0)
endif

if (present(hypre_GMRES_Tol)) then
   solver_control%hypre_cntl%GMRES_Tol = hypre_GMRES_Tol
else
   solver_control%hypre_cntl%GMRES_Tol = huge(0.0d0)
endif

if (present(hypre_GMRES_MaxIter)) then
   solver_control%hypre_cntl%GMRES_MaxIter = hypre_GMRES_MaxIter
else
   solver_control%hypre_cntl%GMRES_MaxIter = huge(0)
endif

if (present(hypre_GMRES_Logging)) then
   solver_control%hypre_cntl%GMRES_Logging = hypre_GMRES_Logging
else
   solver_control%hypre_cntl%GMRES_Logging = huge(0)
endif

! PETSc options

if (present(petsc_richardson_damping_factor)) then
   solver_control%petsc_cntl%petsc_richardson_damping_factor = petsc_richardson_damping_factor
else
   solver_control%petsc_cntl%petsc_richardson_damping_factor = huge(0.0d0)
endif

if (present(petsc_chebychev_emin)) then
   solver_control%petsc_cntl%petsc_chebychev_emin = petsc_chebychev_emin
else
   solver_control%petsc_cntl%petsc_chebychev_emin = huge(0.0d0)
endif

if (present(petsc_chebychev_emax)) then
   solver_control%petsc_cntl%petsc_chebychev_emax = petsc_chebychev_emax
else
   solver_control%petsc_cntl%petsc_chebychev_emax = huge(0.0d0)
endif

if (present(petsc_rtol)) then
   solver_control%petsc_cntl%petsc_rtol = petsc_rtol
else
   solver_control%petsc_cntl%petsc_rtol = huge(0.0d0)
endif

if (present(petsc_atol)) then
   solver_control%petsc_cntl%petsc_atol = petsc_atol
else
   solver_control%petsc_cntl%petsc_atol = huge(0.0d0)
endif

if (present(petsc_dtol)) then
   solver_control%petsc_cntl%petsc_dtol = petsc_dtol
else
   solver_control%petsc_cntl%petsc_dtol = huge(0.0d0)
endif

if (present(petsc_ilu_dt)) then
   solver_control%petsc_cntl%petsc_ilu_dt = petsc_ilu_dt
else
   solver_control%petsc_cntl%petsc_ilu_dt = huge(0.0d0)
endif

if (present(petsc_ilu_dtcol)) then
   solver_control%petsc_cntl%petsc_ilu_dtcol = petsc_ilu_dtcol
else
   solver_control%petsc_cntl%petsc_ilu_dtcol = huge(0.0d0)
endif

if (present(petsc_sor_omega)) then
   solver_control%petsc_cntl%petsc_sor_omega = petsc_sor_omega
else
   solver_control%petsc_cntl%petsc_sor_omega = huge(0.0d0)
endif

if (present(petsc_eisenstat_omega)) then
   solver_control%petsc_cntl%petsc_eisenstat_omega = petsc_eisenstat_omega
else
   solver_control%petsc_cntl%petsc_eisenstat_omega = huge(0.0d0)
endif

if (present(petsc_gmres_max_steps)) then
   solver_control%petsc_cntl%petsc_gmres_max_steps = petsc_gmres_max_steps
else
   solver_control%petsc_cntl%petsc_gmres_max_steps = huge(0)
endif

if (present(petsc_maxits)) then
   solver_control%petsc_cntl%petsc_maxits = petsc_maxits
else
   solver_control%petsc_cntl%petsc_maxits = huge(0)
endif

if (present(petsc_ilu_levels)) then
   solver_control%petsc_cntl%petsc_ilu_levels = petsc_ilu_levels
else
   solver_control%petsc_cntl%petsc_ilu_levels = huge(0)
endif

if (present(petsc_icc_levels)) then
   solver_control%petsc_cntl%petsc_icc_levels = petsc_icc_levels
else
   solver_control%petsc_cntl%petsc_icc_levels = huge(0)
endif

if (present(petsc_ilu_maxrowcount)) then
   solver_control%petsc_cntl%petsc_ilu_maxrowcount = petsc_ilu_maxrowcount
else
   solver_control%petsc_cntl%petsc_ilu_maxrowcount = huge(0)
endif

if (present(petsc_sor_its)) then
   solver_control%petsc_cntl%petsc_sor_its = petsc_sor_its
else
   solver_control%petsc_cntl%petsc_sor_its = huge(0)
endif

if (present(petsc_sor_lits)) then
   solver_control%petsc_cntl%petsc_sor_lits = petsc_sor_lits
else
   solver_control%petsc_cntl%petsc_sor_lits = huge(0)
endif

if (present(petsc_asm_overlap)) then
   solver_control%petsc_cntl%petsc_asm_overlap = petsc_asm_overlap
else
   solver_control%petsc_cntl%petsc_asm_overlap = huge(0)
endif

if (present(petsc_eisenstat_nodiagscaling)) then
   solver_control%petsc_cntl%petsc_eisenstat_nodiagscaling = petsc_eisenstat_nodiagscaling
else
   solver_control%petsc_cntl%petsc_eisenstat_nodiagscaling = .false.
endif

! Check for unsupported combinations of parameters, and correct when possible

if (loc_degree < 0) then
   call warning("degree must be positive.  Using existing degrees")
   loc_degree = 0
endif

if (refine_control%refterm == ONE_REF .and. &
    .not. (present(reftol) .or. present(term_energy_err))) then
   call fatal("refterm == ONE_REF requires either reftol or term_energy_err be present")
   stop
endif

if (PARALLEL/=MPI1 .and. PARALLEL/=MPI2 .and. &
    phaml_solution%eq_type==EIGENVALUE .and. &
    num_proc(phaml_solution%procs)>1) then
   call fatal("Parallel ARPACK requires MPI",procs=procs)
   stop
endif

if (PARALLEL /= MPI1 .and. PARALLEL /= MPI2 .and. &
    num_proc(phaml_solution%procs) > 1 .and. &
    (loc_partition_method == ZOLTAN_RCB .or. &
     loc_partition_method == ZOLTAN_OCT .or. &
     loc_partition_method == ZOLTAN_METIS .or. &
     loc_partition_method == ZOLTAN_REFTREE .or. &
     loc_partition_method == ZOLTAN_RIB .or. &
     loc_partition_method == ZOLTAN_HSFC .or. &
     loc_partition_method == ZOLTAN_FILE)) then
   call fatal("Zoltan requires MPI",procs=procs)
endif

if (PARALLEL /= MPI1 .and. PARALLEL /= MPI2 .and. &
    num_proc(phaml_solution%procs) > 1 .and. &
    (solver_control%solver == PETSC_RICHARDSON_SOLVER .or. &
     solver_control%solver == PETSC_CHEBYCHEV_SOLVER .or. &
     solver_control%solver == PETSC_CG_SOLVER .or. &
     solver_control%solver == PETSC_GMRES_SOLVER .or. &
     solver_control%solver == PETSC_TCQMR_SOLVER .or. &
     solver_control%solver == PETSC_BCGS_SOLVER .or. &
     solver_control%solver == PETSC_CGS_SOLVER .or. &
     solver_control%solver == PETSC_TFQMR_SOLVER .or. &
     solver_control%solver == PETSC_CR_SOLVER .or. &
     solver_control%solver == PETSC_LSQR_SOLVER .or. &
     solver_control%solver == PETSC_BICG_SOLVER)) then
   call fatal("PETSc requires MPI",procs=procs)
endif

if (solver_control%petsc_matrix_free .and. &
    (solver_control%solver == PETSC_JACOBI_PRECONDITION .or. &
     solver_control%solver == PETSC_BJACOBI_PRECONDITION .or. &
     solver_control%solver == PETSC_SOR_PRECONDITION .or. &
     solver_control%solver == PETSC_EISENSTAT_PRECONDITION .or. &
     solver_control%solver == PETSC_ICC_PRECONDITION .or. &
     solver_control%solver == PETSC_ILU_PRECONDITION .or. &
     solver_control%solver == PETSC_ASM_PRECONDITION)) then
   call warning("PETSc preconditioners require the matrix in a PETSc data structure.", &
                "Setting petsc_matrix_free to .false.")
   solver_control%petsc_matrix_free = .false.
endif

if ((solver_control%solver == LAPACK_INDEFINITE_SOLVER .or. &
     solver_control%solver == LAPACK_SPD_SOLVER) .and. &
     num_proc(phaml_solution%procs) > 1) then
   call fatal("LAPACK solvers do not run in parallel; use nproc = 1",procs=procs)
endif

if (solver_control%solver == HYPRE_PCG_SOLVER .and. .not.( &
    solver_control%preconditioner == NO_PRECONDITION .or. &
    solver_control%preconditioner == HYPRE_DS_PRECONDITION .or. &
    solver_control%preconditioner == HYPRE_BOOMERAMG_PRECONDITION)) then
   call fatal("hypre PCG can only use DS, BoomerAMG, or no preconditioner",procs=procs)
endif

if (solver_control%solver == HYPRE_GMRES_SOLVER .and. .not.( &
    solver_control%preconditioner == NO_PRECONDITION .or. &
    solver_control%preconditioner == HYPRE_DS_PRECONDITION .or. &
    solver_control%preconditioner == HYPRE_BOOMERAMG_PRECONDITION .or. &
    solver_control%preconditioner == HYPRE_PARASAILS_PRECONDITION)) then
   call fatal("hypre GMRES can only use DS, BoomerAMG, ParaSails or no preconditioner",procs=procs)
endif

if ((solver_control%preconditioner == HYPRE_BOOMERAMG_PRECONDITION .or. &
    solver_control%preconditioner == HYPRE_DS_PRECONDITION .or. &
    solver_control%preconditioner == HYPRE_PARASAILS_PRECONDITION) .and. .not.(&
    solver_control%solver == HYPRE_GMRES_SOLVER .or. &
    solver_control%solver == HYPRE_PCG_SOLVER)) then
   call fatal("hypre preconditioners can only be used with hypre solvers",procs=procs)
endif

if ((solver_control%solver == CG_SOLVER .or. &
     solver_control%solver == GMRES_SOLVER) .and. &
    (solver_control%preconditioner /= NO_PRECONDITION .and. &
     solver_control%preconditioner /= MG_PRECONDITION)) then
   call fatal("native Krylov solvers can only be preconditioned with NO_PRECONDITION and MG_PRECONDITION",procs=procs)
endif

if (solver_control%coarse_method /= LAPACK_INDEFINITE_SOLVER .and. &
    solver_control%coarse_method /= MUMPS_GEN_SOLVER .and. &
    solver_control%coarse_method /= SUPERLU_SOLVER) then
   call fatal("coarse_method must be one of LAPACK_INDEFINITE_SOLVER, MUMPS_GEN_SOLVER, or SUPERLU_SOLVER",procs=procs)
endif

if (solver_control%coarse_method == MUMPS_GEN_SOLVER .or. &
    solver_control%coarse_method == SUPERLU_SOLVER) then
   call warning("Use of MUMPS or SUPERLU as coarse_method has not yet been implemented.", &
                "Using LAPACK instead.")
   solver_control%coarse_method = LAPACK_INDEFINITE_SOLVER
endif

if (solver_control%ncycle < 1) then
   call warning("mg_cycle must be positive; setting to 1")
   solver_control%ncycle = 1
endif

if (solver_control%mg_tol <= 0.0_my_real .and. &
    solver_control%mg_tol /= MG_NO_TOL .and. &
    solver_control%mg_tol /= MG_ERREST_TOL) then
   call warning("mg_tol must positive, MG_NO_TOL or MG_ERREST_TOL.  Setting to MG_ERREST_TOL")
   solver_control%mg_tol = MG_ERREST_TOL
endif

if (solver_control%num_eval < 1) then
   call warning("num_eval must be at least 1.  Setting to 1")
   solver_control%num_eval = 1
endif

if (solver_control%eigensolver == BLOPEX_SOLVER .and. &
    solver_control%transformation == SHIFT_INVERT .and. &
    solver_control%lambda0_side == EIGEN_BOTH) then
   call warning("BLOPEX_SOLVER does not currently support EIGEN_BOTH with SHIFT_INVERT.  Setting to EIGEN_RIGHT")
   solver_control%lambda0_side = EIGEN_RIGHT
endif

if (solver_control%eigensolver == BLOPEX_SOLVER .and. &
    solver_control%transformation == SHIFT_SQUARE .and. &
    solver_control%lambda0_side /= EIGEN_BOTH) then
   call warning("BLOPEX_SOLVER does not currently support EIGEN_LEFT or EIGEN_RIGHT with SHIFT_SQUARE.  Setting to EIGEN_BOTH")
   solver_control%lambda0_side = EIGEN_BOTH
endif

if (solver_control%eigensolver == ARPACK_SOLVER .and. &
    solver_control%transformation == SHIFT_SQUARE) then
   call warning("ARPACK SOLVER does not currently support SHIFT_SQUARE.  Setting to SHIFT_INVERT")
   solver_control%transformation = SHIFT_INVERT
endif

if (solver_control%arpack_cntl%ncv < 1) then
   call warning("arpack_ncv must be positive.  Setting to the default 20")
   solver_control%arpack_cntl%ncv = 20
endif

if (solver_control%arpack_cntl%maxit < 1) then
   call warning("arpack_maxit must be positive.  Setting to the default 100")
   solver_control%arpack_cntl%maxit = 100
endif

if (solver_control%arpack_cntl%tol <= 0.0_my_real) then
   call warning("arpack_tol must be positive.  Setting to the default 1.0d-10")
   solver_control%arpack_cntl%tol = 1.0e-10_my_real
endif

if (solver_control%blopex_cntl%maxit < 0) then
   call warning("blopex_maxit must be nonnegative.  Setting to the default 100")
   solver_control%blopex_cntl%maxit = 100
endif

if (solver_control%blopex_cntl%atol <= 0.0_my_real) then
   call warning("blopex_atol must be positive.  Setting to the default 1.0d-6")
   solver_control%blopex_cntl%atol = 1.0e-6_my_real
endif

if (solver_control%blopex_cntl%rtol <= 0.0_my_real) then
   call warning("blopex_rtol must be positive.  Setting to the default 1.0d-6")
   solver_control%blopex_cntl%rtol = 1.0e-6_my_real
endif

if (refine_control%reftype == HP_ADAPTIVE .and. &
    refine_control%hp_strategy == HP_SMOOTH_PRED .and. &
    refine_control%derefine) then
   call warning("derefinement not implemented for HP_SMOOTH_PRED strategy.  Setting derefine=.false.")
   refine_control%derefine = .false.
endif

! TEMP until implementation of REFSOLN_EDGE is completed

if (refine_control%reftype == HP_ADAPTIVE .and. &
    refine_control%hp_strategy == HP_REFSOLN_EDGE) then
   call fatal("REFSOLN_EDGE strategy not yet implemented")
   stop
endif

end subroutine defaults

!          ------------
subroutine print_header
!          ------------

!----------------------------------------------------
! This routine prints a header message and values of the parameters
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

if (loc_print_header_who == NO_ONE) return
if ((my_proc(procs) == MASTER .and. loc_print_header_who == SLAVES) .or. &
    (my_proc(procs) /= MASTER .and. loc_print_header_who == MASTER)) return

write(outunit,"(A)") "---------------------------------------------------------"
write(outunit,"(2A)") "           PHAML   Version ",version_number
write(outunit,"(A)")
if (my_proc(procs) == MASTER) then
   write(outunit,"(A)") "Output from the master process"
else
   write(outunit,"(A,I5)") "Output from slave number ",my_proc(procs)
endif
write(outunit,"(A,I5)") "  PDE id is ",phaml_solution%pde_id
write(outunit,"(A,I5)") "  Number of slave processors is ",num_proc(procs)
write(outunit,"(A)")
write(outunit,"(A)") "Subroutine phaml_solve_pde parameters:"
if (refine_control%max_elem == huge(0)) then
   write(outunit,"(A)") "  max_elem           infinity"
else
   write(outunit,"(A,I11)") "  max_elem          ",refine_control%max_elem
endif
if (refine_control%max_vert == huge(0)) then
   write(outunit,"(A)") "  max_vert           infinity"
else
   write(outunit,"(A,I11)") "  max_vert          ",refine_control%max_vert
endif
if (refine_control%max_dof == huge(0)) then
   write(outunit,"(A)") "  max_eq             infinity"
else
   write(outunit,"(A,I11)") "  max_eq            ",refine_control%max_dof
endif
if (refine_control%max_lev == huge(0)) then
   write(outunit,"(A)") "  max_lev            infinity"
else
   write(outunit,"(A,I11)") "  max_lev           ",refine_control%max_lev
endif
if (loc_stop_on_maxlev) then
   write(outunit,"(A)") "  stop_on_maxlev     true"
else
   write(outunit,"(A)") "  stop_on_maxlev     false"
endif
write(outunit,"(A,I11)") "  max_deg           ",refine_control%max_deg
if (loc_stop_on_maxdeg) then
   write(outunit,"(A)") "  stop_on_maxdeg     true"
else
   write(outunit,"(A)") "  stop_on_maxdeg     false"
endif
if (loc_max_refsolveloop == huge(0)) then
   write(outunit,"(A)") "  max_refsolveloop   infinity"
else
   write(outunit,"(A,I11)") "  max_refsolveloop  ",loc_max_refsolveloop
endif
if (refine_control%term_energy_err /= 0.0_my_real) &
   write(outunit,"(SS,1P,A,E18.10E2)") "  term_energy_err   ",refine_control%term_energy_err
if (refine_control%term_Linf_err /= 0.0_my_real) &
   write(outunit,"(SS,1P,A,E18.10E2)") "  term_Linf_err     ",refine_control%term_Linf_err
if (refine_control%term_L2_err /= 0.0_my_real) &
   write(outunit,"(SS,1P,A,E18.10E2)") "  term_L2_err       ",refine_control%term_L2_err
select case(loc_task)
case (BALANCE_REFINE_SOLVE)
   write(outunit,"(A)") "  task               BALANCE_REFINE_SOLVE"
case (SET_INITIAL)
   write(outunit,"(A)") "  task               SET_INITIAL"
case (BALANCE_ONLY)
   write(outunit,"(A)") "  task               BALANCE_ONLY"
case (REFINE_ONLY)
   write(outunit,"(A)") "  task               REFINE_ONLY"
case (SOLVE_ONLY)
   write(outunit,"(A)") "  task               SOLVE_ONLY"
case default
   write(outunit,"(A,I11)") "  task               ",loc_task
end select
if (loc_solve_init) then
   write(outunit,"(A)") "  solve_init         true"
else
   write(outunit,"(A)") "  solve_init         false"
endif
write(outunit,"(A,I11)") "  system_size       ",phaml_solution%system_size
select case(io_control%print_grid_when)
case (NEVER)
   write(outunit,"(A)") "  print_grid_when    NEVER"
case (FINAL)
   write(outunit,"(A)") "  print_grid_when    FINAL"
case (PHASES)
   write(outunit,"(A)") "  print_grid_when    PHASES"
case (FREQUENTLY)
   write(outunit,"(A)") "  print_grid_when    FREQUENTLY"
case (TOO_MUCH)
   write(outunit,"(A)") "  print_grid_when    TOO_MUCH"
case default
   write(outunit,"(A,I11)") "  print_grid_when    ",io_control%print_grid_when
end select
if (io_control%print_grid_when /= NEVER) then
   select case(io_control%print_grid_who)
   case (MASTER)
      write(outunit,"(A)") "  print_grid_who     MASTER"
   case (SLAVES)
      write(outunit,"(A)") "  print_grid_who     SLAVES"
   case (EVERYONE)
      write(outunit,"(A)") "  print_grid_who     EVERYONE"
   case (MASTER_ALL)
      write(outunit,"(A)") "  print_grid_who     MASTER_ALL"
   case (NO_ONE)
      write(outunit,"(A)") "  print_grid_who     NO_ONE"
   case default
      write(outunit,"(A,I11)") "  print_grid_who     ",io_control%print_grid_who
   end select
endif
select case(io_control%print_linsys_when)
case (NEVER)
   write(outunit,"(A)") "  print_linsys_when  NEVER"
case (PHASES)
   write(outunit,"(A)") "  print_linsys_when  PHASES"
case (FREQUENTLY)
   write(outunit,"(A)") "  print_linsys_when  FREQUENTLY"
case (TOO_MUCH)
   write(outunit,"(A)") "  print_linsys_when  TOO_MUCH"
case default
   write(outunit,"(A,I11)") "  print_linsys_when  ",io_control%print_linsys_when
end select
if (io_control%print_linsys_when /= NEVER) then
   select case(io_control%print_linsys_who)
   case (MASTER)
      write(outunit,"(A)") "  print_linsys_who   MASTER"
   case (SLAVES)
      write(outunit,"(A)") "  print_linsys_who   SLAVES"
   case (EVERYONE)
      write(outunit,"(A)") "  print_linsys_who   EVERYONE"
   case (MASTER_ALL)
      write(outunit,"(A)") "  print_linsys_who   MASTER_ALL"
   case (NO_ONE)
      write(outunit,"(A)") "  print_linsys_who   NO_ONE"
   case default
      write(outunit,"(A,I11)") "  print_linsys_who   ",io_control%print_linsys_who
   end select
endif

select case(io_control%print_error_when)
case (NEVER)
   write(outunit,"(A)") "  print_error_when   NEVER"
case (FINAL)
   write(outunit,"(A)") "  print_error_when   FINAL"
case (PHASES)
   write(outunit,"(A)") "  print_error_when   PHASES"
case (FREQUENTLY)
   write(outunit,"(A)") "  print_error_when   FREQUENTLY"
case (TOO_MUCH)
   write(outunit,"(A)") "  print_error_when   TOO_MUCH"
case default
   write(outunit,"(A,I11)") "  print_error_when   ",io_control%print_error_when
end select
if (io_control%print_error_when /= NEVER) then
   select case(io_control%print_error_who)
   case (MASTER)
      write(outunit,"(A)") "  print_error_who    MASTER"
   case (SLAVES)
      write(outunit,"(A)") "  print_error_who    SLAVES"
   case (EVERYONE)
      write(outunit,"(A)") "  print_error_who    EVERYONE"
   case (MASTER_ALL)
      write(outunit,"(A)") "  print_error_who    MASTER_ALL"
   case (NO_ONE)
      write(outunit,"(A)") "  print_error_who    NO_ONE"
   case default
      write(outunit,"(A,I11)") "  print_error_who    ",io_control%print_error_who
   end select
endif
if (io_control%print_error_when /= NEVER) then
   select case(io_control%print_error_what)
   case (NEVER)
      write(outunit,"(A)") "  print_error_what   NEVER"
   case (ENERGY_ERR)
      write(outunit,"(A)") "  print_error_what   ENERGY_ERR"
   case (LINF_ERR)
      write(outunit,"(A)") "  print_error_what   LINF_ERR"
   case (L2_ERR)
      write(outunit,"(A)") "  print_error_what   L2_ERR"
   case (ENERGY_LINF_ERR)
      write(outunit,"(A)") "  print_error_what   ENERGY_LINF_ERR"
   case (ENERGY_L2_ERR)
      write(outunit,"(A)") "  print_error_what   ENERGY_L2_ERR"
   case (LINF_L2_ERR)
      write(outunit,"(A)") "  print_error_what   LINF_L2_ERR"
   case (ENERGY_LINF_L2_ERR)
      write(outunit,"(A)") "  print_error_what   ENERGY_LINF_L2_ERR"
   case default
      write(outunit,"(A,I11)") "  print_error_what   ",io_control%print_error_what
   end select
endif
if (io_control%print_error_when /= NEVER) then
   select case(io_control%print_errest_what)
   case (NEVER)
      write(outunit,"(A)") "  print_errest_what  NEVER"
   case (ENERGY_ERREST)
      write(outunit,"(A)") "  print_errest_what  ENERGY_ERREST"
   case (LINF_ERREST)
      write(outunit,"(A)") "  print_errest_what  LINF_ERREST"
   case (L2_ERREST)
      write(outunit,"(A)") "  print_errest_what  L2_ERREST"
   case (ENERGY_LINF_ERREST)
      write(outunit,"(A)") "  print_errest_what  ENERGY_LINF_ERREST"
   case (ENERGY_L2_ERREST)
      write(outunit,"(A)") "  print_errest_what  ENERGY_L2_ERREST"
   case (LINF_L2_ERREST)
      write(outunit,"(A)") "  print_errest_what   LINF_L2_ERREST"
   case (ENERGY_LINF_L2_ERREST)
      write(outunit,"(A)") "  print_errest_what  ENERGY_LINF_L2_ERREST"
   case default
      write(outunit,"(A,I11)") "  print_errest_what  ",io_control%print_errest_what
   end select
endif
select case(io_control%print_time_when)
case (NEVER)
   write(outunit,"(A)") "  print_time_when    NEVER"
case (FINAL)
   write(outunit,"(A)") "  print_time_when    FINAL"
case (PHASES)
   write(outunit,"(A)") "  print_time_when    PHASES"
case (FREQUENTLY)
   write(outunit,"(A)") "  print_time_when    FREQUENTLY"
case (TOO_MUCH)
   write(outunit,"(A)") "  print_time_when    TOO_MUCH"
case (LAST)
   write(outunit,"(A)") "  print_time_when    LAST"
case (LAST_AND_FINAL)
   write(outunit,"(A)") "  print_time_when    LAST_AND_FINAL"
case default
   write(outunit,"(A,I11)") "  print_time_when    ",io_control%print_time_when
end select
if (io_control%print_time_when /= NEVER) then
   select case(io_control%print_time_who)
   case (MASTER)
      write(outunit,"(A)") "  print_time_who     MASTER"
   case (SLAVES)
      write(outunit,"(A)") "  print_time_who     SLAVES"
   case (EVERYONE)
      write(outunit,"(A)") "  print_time_who     EVERYONE"
   case (MASTER_ALL)
      write(outunit,"(A)") "  print_time_who     MASTER_ALL"
   case (NO_ONE)
      write(outunit,"(A)") "  print_time_who     NO_ONE"
   case default
      write(outunit,"(A,I11)") "  print_time_who     ",io_control%print_time_who
   end select
endif
if (phaml_solution%eq_type == EIGENVALUE) then
   select case(loc_print_eval_when)
   case (NEVER)
      write(outunit,"(A)") "  print_eval_when    NEVER"
   case (FINAL)
      write(outunit,"(A)") "  print_eval_when    FINAL"
   case (PHASES)
      write(outunit,"(A)") "  print_eval_when    PHASES"
   case default
      write(outunit,"(A,I11)") "  print_eval_when    ",loc_print_eval_when
   end select
   if (loc_print_eval_when /= NEVER) then
      select case(loc_print_eval_who)
      case (MASTER)
         write(outunit,"(A)") "  print_eval_who     MASTER"
      case (SLAVES)
         write(outunit,"(A)") "  print_eval_who     SLAVES"
      case (EVERYONE)
         write(outunit,"(A)") "  print_eval_who     EVERYONE"
      case (NO_ONE)
         write(outunit,"(A)") "  print_eval_who     NO_ONE"
      case default
         write(outunit,"(A,I11)") "  print_eval_who     ",loc_print_eval_who
      end select
   endif
endif
if (io_control%print_time_when /= NEVER) then
   select case(loc_clocks)
   case (CLOCK_W)
      write(outunit,"(A)") "  clocks             wall"
   case (CLOCK_C)
      write(outunit,"(A)") "  clocks             cpu"
   case (CLOCK_CW)
      write(outunit,"(A)") "  clocks             cpu and wall"
   case default
      write(outunit,"(A,I11)") "  clocks             ",loc_clocks
   end select
endif
if ((my_proc(procs) == MASTER .and. phaml_solution%master_draws_grid) .or. &
    (my_proc(procs) /= MASTER .and. phaml_solution%i_draw_grid)) then
   select case(io_control%draw_grid_when)
   case (NEVER)
      write(outunit,"(A)") "  draw_grid_when     NEVER"
   case (FINAL)
      write(outunit,"(A)") "  draw_grid_when     FINAL"
   case (PHASES)
      write(outunit,"(A)") "  draw_grid_when     PHASES"
   case (FREQUENTLY)
      write(outunit,"(A)") "  draw_grid_when     FREQUENTLY"
   case (TOO_MUCH)
      write(outunit,"(A)") "  draw_grid_when     TOO_MUCH"
   case default
      write(outunit,"(A,I11)") "  draw_grid_when     ",io_control%draw_grid_when
   end select
endif
if (io_control%pause_after_draw) then
   write(outunit,"(A)") "  pause_after_draw   true"
else
   write(outunit,"(A)") "  pause_after_draw   false"
endif
if (present(pause_after_phases)) then
   if (pause_after_phases) then
      write(outunit,"(A)") "  pause_after_phases true"
   else
      write(outunit,"(A)") "  pause_after_phases false"
   endif
else
   write(outunit,"(A)") "  pause_after_phases false"
endif
if (present(pause_at_start)) then
   if (pause_at_start) then
      write(outunit,"(A)") "  pause_at_start     true"
   else
      write(outunit,"(A)") "  pause_at_start     false"
   endif
else
   write(outunit,"(A)") "  pause_at_start     false"
endif
if (present(pause_at_end)) then
   if (pause_at_end) then
      write(outunit,"(A)") "  pause_at_end       true"
   else
      write(outunit,"(A)") "  pause_at_end       false"
   endif
else
   write(outunit,"(A)") "  pause_at_end       false"
endif
if (loc_degree == 0) then
   write(outunit,"(A,I11)") "  degree             existing"
else
   write(outunit,"(A,I11)") "  degree            ",loc_degree
endif
write(outunit,"(A,I11)") "  sequential_vert   ",loc_sequential_vert
write(outunit,"(SS,1P,A,E18.10E2)") "  inc_factor        ",refine_control%inc_factor
select case(refine_control%error_estimator)
case (HIERARCHICAL_COEFFICIENT)
   write(outunit,"(A)") "  error_estimator    HIERARCHICAL_COEFFICIENT"
case (TRUE_DIFF)
   write(outunit,"(A)") "  error_estimator    TRUE_DIFF"
case (LOCAL_PROBLEM_H)
   write(outunit,"(A)") "  error_estimator    LOCAL_PROBLEM_H"
case (LOCAL_PROBLEM_P)
   write(outunit,"(A)") "  error_estimator    LOCAL_PROBLEM_P"
case (INITIAL_CONDITION)
   write(outunit,"(A)") "  error_estimator    INITIAL_CONDITION"
case (EXPLICIT_ERRIND)
   write(outunit,"(A)") "  error_estimator    EXPLICIT_ERRIND"
case (EQUILIBRATED_RESIDUAL)
   write(outunit,"(A)") "  error_estimator    EQUILIBRATED_RESIDUAL"
case (REFSOLN_ERREST)
   write(outunit,"(A)") "  error_estimator    REFSOLN_ERREST"
case default
   write(outunit,"(A,I11)") "  error_estimator    ",refine_control%error_estimator
end select
select case(phaml_solution%grid%errtype)
case (ABSOLUTE_ERROR)
   write(outunit,"(A)") "  errtype            ABSOLUTE_ERROR"
case (RELATIVE_ERROR)
   write(outunit,"(A)") "  errtype            RELATIVE ERROR"
case default
   write(outunit,"(A,I11)") "  errtype            ",phaml_solution%grid%errtype
end select
select case(refine_control%reftype)
case (H_UNIFORM)
   write(outunit,"(A)") "  reftype            H_UNIFORM"
case (H_ADAPTIVE)
   write(outunit,"(A)") "  reftype            H_ADAPTIVE"
case (P_UNIFORM)
   write(outunit,"(A)") "  reftype            P_UNIFORM"
case (P_ADAPTIVE)
   write(outunit,"(A)") "  reftype            P_ADAPTIVE"
case (HP_ADAPTIVE)
   write(outunit,"(A)") "  reftype            HP_ADAPTIVE"
case default
   write(outunit,"(A,I11)") "  reftype            ",refine_control%reftype
end select
select case (refine_control%edge_rule)
case (MINIMUM_RULE)
   write(outunit,"(A)") "  edge_rule          MINIMUM_RULE"
case (MAXIMUM_RULE)
   write(outunit,"(A)") "  edge_rule          MAXIMUM_RULE"
case default
   write(outunit,"(A,I11)") "  edge_rule          ",refine_control%edge_rule
end select
if (refine_control%reftype == HP_ADAPTIVE) then
   select case(refine_control%hp_strategy)
   case (HP_BIGGER_ERRIND)
      write(outunit,"(A)") "  hp_strategy        HP_BIGGER_ERRIND"
   case (HP_APRIORI)
      write(outunit,"(A)") "  hp_strategy        HP_APRIORI"
   case (HP_PRIOR2P_E)
      write(outunit,"(A)") "  hp_strategy        HP_PRIOR2P_E"
   case (HP_PRIOR2P_H1)
      write(outunit,"(A)") "  hp_strategy        HP_PRIOR2P_H1"
   case (HP_T3S)
      write(outunit,"(A)") "  hp_strategy        HP_T3S"
   case (HP_ALTERNATE)
      write(outunit,"(A)") "  hp_strategy        HP_ALTERNATE"
   case (HP_TYPEPARAM)
      write(outunit,"(A)") "  hp_strategy        HP_TYPEPARAM"
   case (HP_COEF_DECAY)
      write(outunit,"(A)") "  hp_strategy        HP_COEF_DECAY"
   case (HP_COEF_ROOT)
      write(outunit,"(A)") "  hp_strategy        HP_COEF_ROOT"
   case (HP_SMOOTH_PRED)
      write(outunit,"(A)") "  hp_strategy        HP_SMOOTH_PRED"
   case (HP_NEXT3P)
      write(outunit,"(A)") "  hp_strategy        HP_NEXT3P"
   case (HP_REFSOLN_EDGE)
      write(outunit,"(A)") "  hp_strategy        HP_REFSOLN_EDGE"
   case (HP_REFSOLN_ELEM)
      write(outunit,"(A)") "  hp_strategy        HP_REFSOLN_ELEM"
   case (HP_NLP)
      write(outunit,"(A)") "  hp_strategy        HP_NLP"
   case default
      write(outunit,"(A,I11)") "  hp_strategy        ",refine_control%hp_strategy
   end select
endif
if (refine_control%reftype == HP_ADAPTIVE .and. &
    (refine_control%hp_strategy == HP_T3S .or. &
     refine_control%hp_strategy == HP_ALTERNATE)) then
   write(outunit,"(SS,1P,A,E18.10E2)") "  t3s_gamma         ",refine_control%t3s_gamma
   write(outunit,"(SS,1P,A,E18.10E2)") "  t3s_eta           ",refine_control%t3s_eta
   write(outunit,"(A,I11)") "  t3s_nunif          ",refine_control%t3s_nunif
   write(outunit,"(A,I11)") "  t3s_maxref         ",refine_control%t3s_maxref
   write(outunit,"(A,I11)") "  t3s_maxdeginc      ",refine_control%t3s_maxdeginc
endif
if (refine_control%reftype == HP_ADAPTIVE .and. &
    refine_control%hp_strategy == HP_TYPEPARAM) then
   write(outunit,"(SS,1P,A,E18.10E2)") "  tp_gamma          ",refine_control%tp_gamma
endif
if (refine_control%reftype == HP_ADAPTIVE .and. &
    refine_control%hp_strategy == HP_SMOOTH_PRED) then
   write(outunit,"(SS,1P,A,E18.10E2)") "  sp_gamma_h        ",refine_control%sp_gamma_h
   write(outunit,"(SS,1P,A,E18.10E2)") "  sp_gamma_p        ",refine_control%sp_gamma_p
endif
if (refine_control%reftype == HP_ADAPTIVE .and. &
    refine_control%hp_strategy == HP_NLP) then
   write(outunit,"(A,I11)") "  nlp_max_h_dec     ",refine_control%nlp_max_h_dec
   write(outunit,"(A,I11)") "  nlp_max_h_inc     ",refine_control%nlp_max_h_inc
   write(outunit,"(A,I11)") "  nlp_max_p_dec     ",refine_control%nlp_max_p_dec
   write(outunit,"(A,I11)") "  nlp_max_p_inc     ",refine_control%nlp_max_p_inc
endif
if (refine_control%reftype == HP_ADAPTIVE .and. &
    refine_control%hp_strategy == HP_REFSOLN_ELEM) then
   write(outunit,"(SS,1P,A,E18.10E2)") "  refsoln_pbias     ",refine_control%refsoln_pbias
endif
select case(refine_control%refterm)
case (DOUBLE_NVERT)
   write(outunit,"(A)") "  refterm            DOUBLE_NVERT"
case (DOUBLE_NVERT_SMOOTH)
   write(outunit,"(A)") "  refterm            DOUBLE_NVERT_SMOOTH"
case (DOUBLE_NELEM)
   write(outunit,"(A)") "  refterm            DOUBLE_NELEM"
case (DOUBLE_NELEM_SMOOTH)
   write(outunit,"(A)") "  refterm            DOUBLE_NELEM_SMOOTH"
case (DOUBLE_NEQ)
   write(outunit,"(A)") "  refterm            DOUBLE_NEQ"
case (HALVE_ERREST)
   write(outunit,"(A)") "  refterm            HALVE_ERREST"
case (KEEP_NVERT)
   write(outunit,"(A)") "  refterm            KEEP_NVERT"
case (KEEP_NVERT_SMOOTH)
   write(outunit,"(A)") "  refterm            KEEP_NVERT_SMOOTH"
case (KEEP_NELEM)
   write(outunit,"(A)") "  refterm            KEEP_NELEM"
case (KEEP_NELEM_SMOOTH)
   write(outunit,"(A)") "  refterm            KEEP_NELEM_SMOOTH"
case (KEEP_NEQ)
   write(outunit,"(A)") "  refterm            KEEP_NEQ"
case (KEEP_NEQ_SMOOTH)
   write(outunit,"(A)") "  refterm            KEEP_NEQ_SMOOTH"
case (KEEP_ERREST)
   write(outunit,"(A)") "  refterm            KEEP_ERREST"
case (DOUBLE_NEQ_SMOOTH)
   write(outunit,"(A)") "  refterm            DOUBLE_NEQ_SMOOTH"
case (ONE_REF)
   write(outunit,"(A)") "  refterm            ONE_REF"
case (ONE_REF_HALF_ERRIND)
   write(outunit,"(A)") "  refterm            ONE_REF_HALF_ERRIND"
case default
   write(outunit,"(A,I11)") "  refterm            ",refine_control%refterm
end select
if (refine_control%refterm == ONE_REF) then
   write(outunit,"(SS,1P,A,E18.10E2)") "  reftol             ",refine_control%reftol
endif
if (refine_control%derefine) then
   write(outunit,"(A)") "  derefine           true"
else
   write(outunit,"(A)") "  derefine           false"
endif
select case(loc_partition_method)
case (RTK)
   write(outunit,"(A)") "  partition_method   RTK"
case (ZOLTAN_RCB)
   write(outunit,"(A)") "  partition_method   ZOLTAN_RCB"
case (ZOLTAN_OCT)
   write(outunit,"(A)") "  partition_method   ZOLTAN_OCT"
case (ZOLTAN_METIS)
   write(outunit,"(A)") "  partition_method   ZOLTAN_METIS"
case (ZOLTAN_REFTREE)
   write(outunit,"(A)") "  partition_method   ZOLTAN_REFTREE"
case (ZOLTAN_RIB)
   write(outunit,"(A)") "  partition_method   ZOLTAN_RIB"
case (ZOLTAN_HSFC)
   write(outunit,"(A)") "  partition_method   ZOLTAN_HSFC"
case (ZOLTAN_FILE)
   write(outunit,"(A)") "  partition_method   ZOLTAN_FILE"
   write(outunit,"(A,A)") "  zoltan_param_file  ",trim(phaml_solution%zoltan_param_file)
case default
   write(outunit,"(A,I11)") "  partition_method   ",loc_partition_method
end select
select case(loc_prebalance)
case (BALANCE_NONE)
case (BALANCE_ELEMENTS)
   write(outunit,"(A)") "  prebalance         BALANCE_ELEMENTS"
case (BALANCE_VERTICES)
   write(outunit,"(A)") "  prebalance         BALANCE_VERTICES"
case (BALANCE_EQUATIONS)
   write(outunit,"(A)") "  prebalance         BALANCE_EQUATIONS"
case default
   write(outunit,"(A,I11)") "  prebalance         ",loc_prebalance
end select
select case(loc_postbalance)
case (BALANCE_NONE)
case (BALANCE_ELEMENTS)
   write(outunit,"(A)") "  postbalance        BALANCE_ELEMENTS"
case (BALANCE_VERTICES)
   write(outunit,"(A)") "  postbalance        BALANCE_VERTICES"
case (BALANCE_EQUATIONS)
   write(outunit,"(A)") "  postbalance        BALANCE_EQUATIONS"
case default
   write(outunit,"(A,I11)") "  postbalance        ",loc_prebalance
end select
if (solver_control%ignore_quad_err) then
   write(outunit,"(A)") "  ignore_quad_err    true"
else
   write(outunit,"(A)") "  ignore_quad_err    false"
endif
select case(solver_control%solver)
case (MG_SOLVER)
   write(outunit,"(A)") "  solver             MG_SOLVER"
case (CG_SOLVER)
   write(outunit,"(A)") "  solver             CG_SOLVER"
case (GMRES_SOLVER)
   write(outunit,"(A)") "  solver             GMRES_SOLVER"
case (LAPACK_INDEFINITE_SOLVER)
   write(outunit,"(A)") "  solver             LAPACK_INDEFINITE_SOLVER"
case (LAPACK_SPD_SOLVER)
   write(outunit,"(A)") "  solver             LAPACK_SPD_SOLVER"
case (MUMPS_NONSYM_SOLVER)
   write(outunit,"(A)") "  solver             MUMPS_NONSYM_SOLVER"
case (MUMPS_SPD_SOLVER)
   write(outunit,"(A)") "  solver             MUMPS_SPD_SOLVER"
case (MUMPS_GEN_SOLVER)
   write(outunit,"(A)") "  solver             MUMPS_GEN_SOLVER"
case (SUPERLU_SOLVER)
   write(outunit,"(A)") "  solver             SUPERLU_SOLVER"
case (PETSC_RICHARDSON_SOLVER)
   write(outunit,"(A)") "  solver             PETSC_RICHARDSON_SOLVER"
case (PETSC_CHEBYCHEV_SOLVER)
   write(outunit,"(A)") "  solver             PETSC_CHEBYCHEV_SOLVER"
case (PETSC_CG_SOLVER)
   write(outunit,"(A)") "  solver             PETSC_CG_SOLVER"
case (PETSC_GMRES_SOLVER)
   write(outunit,"(A)") "  solver             PETSC_GMRES_SOLVER"
case (PETSC_TCQMR_SOLVER)
   write(outunit,"(A)") "  solver             PETSC_TCQMR_SOLVER"
case (PETSC_BCGS_SOLVER)
   write(outunit,"(A)") "  solver             PETSC_BCGS_SOLVER"
case (PETSC_CGS_SOLVER)
   write(outunit,"(A)") "  solver             PETSC_CGS_SOLVER"
case (PETSC_TFQMR_SOLVER)
   write(outunit,"(A)") "  solver             PETSC_TFQMR_SOLVER"
case (PETSC_CR_SOLVER)
   write(outunit,"(A)") "  solver             PETSC_CR_SOLVER"
case (PETSC_LSQR_SOLVER)
   write(outunit,"(A)") "  solver             PETSC_LSQR_SOLVER"
case (PETSC_BICG_SOLVER)
   write(outunit,"(A)") "  solver             PETSC_BICG_SOLVER"
case (HYPRE_BOOMERAMG_SOLVER)
   write(outunit,"(A)") "  solver             HYPRE_BOOMERAMG_SOLVER"
case (HYPRE_PCG_SOLVER)
   write(outunit,"(A)") "  solver             HYPRE_PCG_SOLVER"
case (HYPRE_GMRES_SOLVER)
   write(outunit,"(A)") "  solver             HYPRE_GMRES_SOLVER"
case default
   write(outunit,"(A,I11)") "  solver             ",solver_control%solver
end select
if (solver_control%solver == PETSC_RICHARDSON_SOLVER .or. &
    solver_control%solver == PETSC_CHEBYCHEV_SOLVER  .or. &
    solver_control%solver == PETSC_CG_SOLVER         .or. &
    solver_control%solver == PETSC_GMRES_SOLVER      .or. &
    solver_control%solver == PETSC_TCQMR_SOLVER      .or. &
    solver_control%solver == PETSC_BCGS_SOLVER       .or. &
    solver_control%solver == PETSC_CGS_SOLVER        .or. &
    solver_control%solver == PETSC_TFQMR_SOLVER      .or. &
    solver_control%solver == PETSC_CR_SOLVER         .or. &
    solver_control%solver == PETSC_LSQR_SOLVER       .or. &
    solver_control%solver == PETSC_BICG_SOLVER       .or. &
    solver_control%solver == CG_SOLVER               .or. &
    solver_control%solver == GMRES_SOLVER            .or. &
    solver_control%solver == HYPRE_PCG_SOLVER        .or. &
    solver_control%solver == HYPRE_GMRES_SOLVER) then
   select case(solver_control%preconditioner)
   case (NO_PRECONDITION)
      write(outunit,"(A)") "  preconditioner     NO_PRECONDITION"
   case (MG_PRECONDITION)
      write(outunit,"(A)") "  preconditioner     MG_PRECONDITION"
   case (FUDOP_DD_PRECONDITION)
      write(outunit,"(A)") "  preconditioner     FUDOP_DD_PRECONDITION"
   case (COARSE_GRID_PRECONDITION)
      write(outunit,"(A)") "  preconditioner     COARSE_GRID_PRECONDITION"
   case (PETSC_JACOBI_PRECONDITION)
      write(outunit,"(A)") "  preconditioner     PETSC_JACOBI_PRECONDITION"
   case (PETSC_BJACOBI_PRECONDITION)
      write(outunit,"(A)") "  preconditioner     PETSC_BJACOBI_PRECONDITION"
   case (PETSC_SOR_PRECONDITION)
      write(outunit,"(A)") "  preconditioner     PETSC_SOR_PRECONDITION"
   case (PETSC_EISENSTAT_PRECONDITION)
      write(outunit,"(A)") "  preconditioner     PETSC_EISENSTAT_PRECONDITION"
   case (PETSC_ICC_PRECONDITION)
      write(outunit,"(A)") "  preconditioner     PETSC_ICC_PRECONDITION"
   case (PETSC_ILU_PRECONDITION)
      write(outunit,"(A)") "  preconditioner     PETSC_ILU_PRECONDITION"
   case (PETSC_ASM_PRECONDITION)
      write(outunit,"(A)") "  preconditioner     PETSC_ASM_PRECONDITION"
   case (HYPRE_DS_PRECONDITION)
      write(outunit,"(A)") "  preconditioner     HYPRE_DS_PRECONDITION"
   case (HYPRE_BOOMERAMG_PRECONDITION)
      write(outunit,"(A)") "  preconditioner     HYPRE_BOOMERAMG_PRECONDITION"
   case (HYPRE_PARASAILS_PRECONDITION)
      write(outunit,"(A)") "  preconditioner     HYPRE_PARASAILS_PRECONDITION"
   case default
      write(outunit,"(A,I11)") "  preconditioner     ",solver_control%preconditioner
   end select
endif
if (solver_control%preconditioner == COARSE_GRID_PRECONDITION) then
   write(outunit,"(A,I11)") "  coarse_size       ",solver_control%coarse_size
   select case (solver_control%coarse_method)
   case (LAPACK_INDEFINITE_SOLVER)
      write(outunit,"(A)") "  coarse_method      LAPACK_INDEFINITE_SOLVER"
   case (MUMPS_GEN_SOLVER)
      write(outunit,"(A)") "  coarse_method      MUMPS_GEN_SOLVER"
   case (SUPERLU_SOLVER)
      write(outunit,"(A)") "  coarse_method      SUPERLU_SOLVER"
   case default
      write(outunit,"(A,I11)") "  coarse_method      ",solver_control%coarse_method
   end select
endif
if (solver_control%inc_quad_order /= 0) then
   write(outunit,"(A,I11)") "  inc_quad_order    ",solver_control%inc_quad_order
endif
if (solver_control%solver == MG_SOLVER .or. &
    solver_control%preconditioner == MG_PRECONDITION) then
   if (solver_control%ncycle == huge(0)) then
      write(outunit,"(A)") "  mg_cycles          infinity"
   else
      write(outunit,"(A,I11)") "  mg_cycles         ",solver_control%ncycle
   endif
   if (solver_control%mg_tol == MG_NO_TOL) then
      write(outunit,"(A)") "  mg_tol             MG_NO_TOL"
   elseif (solver_control%mg_tol == MG_ERREST_TOL) then
      write(outunit,"(A)") "  mg_tol             MG_ERREST_TOL"
   else
      write(outunit,"(SS,1P,A,E18.10E2)") "  mg_tol            ",solver_control%mg_tol
   endif
   write(outunit,"(A,I11)") "  mg_prerelax       ",solver_control%prerelax
   write(outunit,"(A,I11)") "  mg_postrelax      ",solver_control%postrelax
   write(outunit,"(A,I11)") "  mg_prerelax_ho    ",solver_control%prerelax_ho
   write(outunit,"(A,I11)") "  mg_postrelax_ho   ",solver_control%postrelax_ho
   select case (solver_control%mg_comm)
   case (MGCOMM_NONE)
      write(outunit,"(A)") "  mg_comm            MGCOMM_NONE"
   case (MGCOMM_FUDOP)
      write(outunit,"(A)") "  mg_comm            MGCOMM_FUDOP"
   case (MGCOMM_CONVENTIONAL)
      write(outunit,"(A)") "  mg_comm            MGCOMM_CONVENTIONAL"
   case default
      write(outunit,"(A,I11)") "  mg_comm           ",solver_control%mg_comm
   end select
endif
if (solver_control%solver == CG_SOLVER .or. &
    solver_control%solver == GMRES_SOLVER) then
   write(outunit,"(A,I11)") "  krylov_iter       ",solver_control%krylov_iter
   if (solver_control%solver == GMRES_SOLVER) then
      write(outunit,"(A,I11)") "  krylov_restart    ",solver_control%krylov_restart
   endif
   if (solver_control%krylov_tol == KRYLOV_ERREST_TOL) then
      write(outunit,"(A)") "  krylov_tol         KRYLOV_ERREST_TOL"
   else
      write(outunit,"(SS,1P,A,E18.10E2)") "  krylov_tol        ",solver_control%krylov_tol
   endif
endif
if (solver_control%preconditioner == FUDOP_DD_PRECONDITION) then
   write(outunit,"(A,I11)") "  dd_iterations     ",solver_control%dd_iterations
endif
if (solver_control%solver == PETSC_RICHARDSON_SOLVER .or. &
    solver_control%solver == PETSC_CHEBYCHEV_SOLVER  .or. &
    solver_control%solver == PETSC_CG_SOLVER         .or. &
    solver_control%solver == PETSC_GMRES_SOLVER      .or. &
    solver_control%solver == PETSC_TCQMR_SOLVER      .or. &
    solver_control%solver == PETSC_BCGS_SOLVER       .or. &
    solver_control%solver == PETSC_CGS_SOLVER        .or. &
    solver_control%solver == PETSC_TFQMR_SOLVER      .or. &
    solver_control%solver == PETSC_CR_SOLVER         .or. &
    solver_control%solver == PETSC_LSQR_SOLVER       .or. &
    solver_control%solver == PETSC_BICG_SOLVER) then
   if (solver_control%petsc_matrix_free) then
      write(outunit,"(A)") "  petsc_matrix_free  true"
   else
      write(outunit,"(A)") "  petsc_matrix_free  false"
   endif
   if (solver_control%petsc_cntl%petsc_rtol == huge(0.0d0)) then
      write(outunit,"(A)") "  petsc_rtol         default"
   else
      write(outunit,"(SS,1P,A,E18.10E2)") "  petsc_rtol         ",solver_control%petsc_cntl%petsc_rtol
   endif
   if (solver_control%petsc_cntl%petsc_atol == huge(0.0d0)) then
      write(outunit,"(A)") "  petsc_atol         default"
   else
      write(outunit,"(SS,1P,A,E18.10E2)") "  petsc_atol         ",solver_control%petsc_cntl%petsc_atol
   endif
   if (solver_control%petsc_cntl%petsc_dtol == huge(0.0d0)) then
      write(outunit,"(A)") "  petsc_dtol         default"
   else
      write(outunit,"(SS,1P,A,E18.10E2)") "  petsc_dtol         ",solver_control%petsc_cntl%petsc_dtol
   endif
   if (solver_control%petsc_cntl%petsc_maxits == huge(0)) then
      write(outunit,"(A)") "  petsc_maxits       default"
   else
      write(outunit,"(A,I11)") "  petsc_maxits       ",solver_control%petsc_cntl%petsc_maxits
   endif
   if (solver_control%solver == PETSC_RICHARDSON_SOLVER) then
      if (solver_control%petsc_cntl%petsc_richardson_damping_factor == huge(0.0d0)) then
         write(outunit,"(A)") "  petsc_richardson_damping_factor   default"
      else
         write(outunit,"(SS,1P,A,E18.10E2)") "  petsc_richardson_damping_factor   ", &
         solver_control%petsc_cntl%petsc_richardson_damping_factor
      endif
   endif
   if (solver_control%solver == PETSC_CHEBYCHEV_SOLVER) then
      if (solver_control%petsc_cntl%petsc_chebychev_emin == huge(0.0d0)) then
         write(outunit,"(A)") "  petsc_chebychev_emin  default"
      else
         write(outunit,"(SS,1P,A,E18.10E2)") "  petsc_chebychev_emin  ",solver_control%petsc_cntl%petsc_chebychev_emin
      endif
      if (solver_control%petsc_cntl%petsc_chebychev_emax == huge(0.0d0)) then
         write(outunit,"(A)") "  petsc_chebychev_emax  default"
      else
         write(outunit,"(SS,1P,A,E18.10E2)") "  petsc_chebychev_emax  ",solver_control%petsc_cntl%petsc_chebychev_emax
      endif
   endif
   if (solver_control%solver == PETSC_GMRES_SOLVER) then
      if (solver_control%petsc_cntl%petsc_gmres_max_steps == huge(0)) then
         write(outunit,"(A)") "  petsc_gmres_max_steps  default"
      else
         write(outunit,"(A,I11)") "  petsc_gmres_max_steps  ",solver_control%petsc_cntl%petsc_gmres_max_steps
      endif
   endif
   if (solver_control%preconditioner == PETSC_ILU_PRECONDITION) then
      if (solver_control%petsc_cntl%petsc_ilu_levels == huge(0)) then
         write(outunit,"(A)") "  petsc_ilu_levels   default"
      else
         write(outunit,"(A,I11)") "  petsc_ilu_levels   ",solver_control%petsc_cntl%petsc_ilu_levels
      endif
      if (solver_control%petsc_cntl%petsc_ilu_dt == huge(0.0d0)) then
         write(outunit,"(A)") "  petsc_ilu_dt       default"
      else
         write(outunit,"(SS,1P,A,E18.10E2)") "  petsc_ilu_dt                      ",solver_control%petsc_cntl%petsc_ilu_dt
      endif
      if (solver_control%petsc_cntl%petsc_ilu_dtcol == huge(0.0d0)) then
         write(outunit,"(A)") "  petsc_ilu_dtcol    default"
      else
         write(outunit,"(SS,1P,A,E18.10E2)") "  petsc_ilu_dtcol    ",solver_control%petsc_cntl%petsc_ilu_dtcol
      endif
      if (solver_control%petsc_cntl%petsc_ilu_maxrowcount == huge(0)) then
         write(outunit,"(A)") "  petsc_ilu_maxrowcount  default"
      else
         write(outunit,"(A,I11)") "  petsc_ilu_maxrowcount  ",solver_control%petsc_cntl%petsc_ilu_maxrowcount
      endif
   endif
   if (solver_control%preconditioner == PETSC_ILU_PRECONDITION) then
      if (solver_control%petsc_cntl%petsc_icc_levels == huge(0)) then
         write(outunit,"(A)") "  petsc_icc_levels   default"
      else
         write(outunit,"(A,I11)") "  petsc_icc_levels   ",solver_control%petsc_cntl%petsc_icc_levels
      endif
   endif
   if (solver_control%preconditioner == PETSC_SOR_PRECONDITION) then
      if (solver_control%petsc_cntl%petsc_sor_omega == huge(0.0d0)) then
         write(outunit,"(A)") "  petsc_sor_omega    default"
      else
         write(outunit,"(SS,1P,A,E18.10E2)") "  petsc_sor_omega    ",solver_control%petsc_cntl%petsc_sor_omega
      endif
      if (solver_control%petsc_cntl%petsc_sor_its == huge(0)) then
         write(outunit,"(A)") "  petsc_sor_its      default"
      else
         write(outunit,"(A,I11)") "  petsc_sor_its      ",solver_control%petsc_cntl%petsc_sor_its
      endif
      if (solver_control%petsc_cntl%petsc_sor_lits == huge(0)) then
         write(outunit,"(A)") "  petsc_sor_lits     default"
      else
         write(outunit,"(A,I11)") "  petsc_sor_lits     ",solver_control%petsc_cntl%petsc_sor_lits
      endif
   endif
   if (solver_control%preconditioner == PETSC_EISENSTAT_PRECONDITION) then
      write(outunit,"(A,L1)") "  petsc_eisenstat_nodiagscaling ",solver_control%petsc_cntl%petsc_eisenstat_nodiagscaling
      if (solver_control%petsc_cntl%petsc_eisenstat_omega == huge(0.0d0)) then
         write(outunit,"(A)") "  petsc_eisenstat_omega  default"
      else
         write(outunit,"(SS,1P,A,E18.10E2)") "  petsc_eisenstat_omega  ",solver_control%petsc_cntl%petsc_eisenstat_omega
      endif
   endif
   if (solver_control%preconditioner == PETSC_ASM_PRECONDITION) then
      if (solver_control%petsc_cntl%petsc_asm_overlap == huge(0)) then
         write(outunit,"(A)") "  petsc_asm_overlap  default"
      else
         write(outunit,"(A,I11)") "  petsc_asm_overlap  ",solver_control%petsc_cntl%petsc_asm_overlap
      endif
   endif
endif
if (phaml_solution%eq_type == EIGENVALUE) then
   select case (solver_control%eigensolver)
   case (ARPACK_SOLVER)
      write(outunit,"(A)") "  eigensolver        ARPACK_SOLVER"
   case (BLOPEX_SOLVER)
      write(outunit,"(A)") "  eigensolver        BLOPEX_SOLVER"
   case default
      write(outunit,"(A,I11)") "  eigensolver       ",solver_control%eigensolver
   end select
   write(outunit,"(A,I11)") "  num_eval          ",solver_control%num_eval
   if (solver_control%lambda0 == -huge(0.0_my_real)) then
      write(outunit,"(A)") "  lambda0            -infinity"
   else
      write(outunit,"(SS,1P,A,E18.10E2)") "  lambda0           ",solver_control%lambda0
   endif
   select case (solver_control%scale_evec)
   case (SCALE_LINF)
      write(outunit,"(A)") "  scale_evec         SCALE_LINF"
   case (SCALE_L2)
      write(outunit,"(A)") "  scale_evec         SCALE_L2"
   case (SCALE_M)
      write(outunit,"(A)") "  scale_evec         SCALE_M"
   case default
      write(outunit,"(A,I11)") "  scale_evec        ",solver_control%scale_evec
   end select
   if (solver_control%lambda0 /= -huge(0.0_my_real)) then
      select case (solver_control%lambda0_side)
      case (EIGEN_LEFT)
         write(outunit,"(A)") "  lambda0_side       EIGEN_LEFT"
      case (EIGEN_RIGHT)
         write(outunit,"(A)") "  lambda0_side       EIGEN_RIGHT"
      case (EIGEN_BOTH)
         write(outunit,"(A)") "  lambda0_side       EIGEN_BOTH"
      case default
         write(outunit,"(A,I11)") "  lambda0_side      ",solver_control%lambda0_side
      end select
      select case (solver_control%transformation)
      case (SHIFT_INVERT)
         write(outunit,"(A)") "  transformation     SHIFT_INVERT"
      case (SHIFT_SQUARE)
         write(outunit,"(A)") "  transformation     SHIFT_SQUARE"
      case default
         write(outunit,"(A,I11)") "  transformation    ",solver_control%transformation
      end select
   endif
   select case (solver_control%eigensolver)
   case (ARPACK_SOLVER)
      write(outunit,"(A,I11)") "  arpack_ncv        ",solver_control%arpack_cntl%ncv
      write(outunit,"(A,I11)") "  arpack_maxit      ",solver_control%arpack_cntl%maxit
      write(outunit,"(SS,1P,A,E18.10E2)") "  arpack_tol        ",solver_control%arpack_cntl%tol
   case (BLOPEX_SOLVER)
      write(outunit,"(A,I11)") "  blopex_maxit      ",solver_control%blopex_cntl%maxit
      write(outunit,"(SS,1P,A,E18.10E2)") "  blopex_atol       ",solver_control%blopex_cntl%atol
      write(outunit,"(SS,1P,A,E18.10E2)") "  blopex_rtol       ",solver_control%blopex_cntl%rtol
   end select
endif
if (solver_control%solver == HYPRE_PCG_SOLVER) then
   if (solver_control%hypre_cntl%PCG_Tol == huge(0.0d0)) then
      write(outunit,"(A)") "  hypre_PCG_Tol                    default"
   else
      write(outunit,"(SS,1P,A,E18.10E2)") "  hypre_PCG_Tol                   ", solver_control%hypre_cntl%PCG_Tol
   endif
   if (solver_control%hypre_cntl%PCG_MaxIter == huge(0)) then
      write(outunit,"(A)") "  hypre_PCG_MaxIter                default"
   else
      write(outunit,"(A,I11)") "  hypre_PCG_MaxIter               ", solver_control%hypre_cntl%PCG_MaxIter
   endif
   if (solver_control%hypre_cntl%PCG_TwoNorm == huge(0)) then
      write(outunit,"(A)") "  hypre_PCG_TwoNorm                default"
   else
      write(outunit,"(A,I11)") "  hypre_PCG_TwoNorm               ", solver_control%hypre_cntl%PCG_TwoNorm
   endif
   if (solver_control%hypre_cntl%PCG_RelChange == huge(0)) then
      write(outunit,"(A)") "  hypre_PCG_RelChange              default"
   else
      write(outunit,"(A,I11)") "  hypre_PCG_RelChange             ", solver_control%hypre_cntl%PCG_RelChange
   endif
   if (solver_control%hypre_cntl%PCG_Logging == huge(0)) then
      write(outunit,"(A)") "  hypre_PCG_Logging                default"
   else
      write(outunit,"(A,I11)") "  hypre_PCG_Logging               ", solver_control%hypre_cntl%PCG_Logging
   endif
endif
if (solver_control%solver == HYPRE_GMRES_SOLVER) then
   if (solver_control%hypre_cntl%GMRES_KDim == huge(0)) then
      write(outunit,"(A)") "  hypre_GMRES_KDim                 default"
   else
      write(outunit,"(A,I11)") "  hypre_GMRES_KDim                ", solver_control%hypre_cntl%GMRES_KDim
   endif
   if (solver_control%hypre_cntl%GMRES_Tol == huge(0.0d0)) then
      write(outunit,"(A)") "  hypre_GMRES_Tol                  default"
   else
      write(outunit,"(SS,1P,A,E18.10E2)") "  hypre_GMRES_Tol                 ", solver_control%hypre_cntl%GMRES_Tol
   endif
   if (solver_control%hypre_cntl%GMRES_MaxIter == huge(0)) then
      write(outunit,"(A)") "  hypre_GMRES_MaxIter              default"
   else
      write(outunit,"(A,I11)") "  hypre_GMRES_MaxIter             ", solver_control%hypre_cntl%GMRES_MaxIter
   endif
   if (solver_control%hypre_cntl%GMRES_Logging == huge(0)) then
      write(outunit,"(A)") "  hypre_GMRES_Logging              default"
   else
      write(outunit,"(A,I11)") "  hypre_GMRES_Logging             ", solver_control%hypre_cntl%GMRES_Logging
   endif
endif
if (solver_control%solver == HYPRE_BOOMERAMG_SOLVER .or. &
    solver_control%preconditioner == HYPRE_BOOMERAMG_PRECONDITION) then
   if (solver_control%hypre_cntl%BoomerAMG_MaxLevels == huge(0)) then
      write(outunit,"(A)") "  hypre_BoomerAMG_MaxLevels        default"
   else
      write(outunit,"(A,I11)") "  hypre_BoomerAMG_MaxLevels       ", solver_control%hypre_cntl%BoomerAMG_MaxLevels
   endif
   if (solver_control%hypre_cntl%BoomerAMG_MaxIter == huge(0)) then
      write(outunit,"(A)") "  hypre_BoomerAMG_MaxIter          default"
   else
      write(outunit,"(A,I11)") "  hypre_BoomerAMG_MaxIter         ", solver_control%hypre_cntl%BoomerAMG_MaxIter
   endif
   if (solver_control%hypre_cntl%BoomerAMG_Tol == huge(0.0d0)) then
      write(outunit,"(A)") "  hypre_BoomerAMG_Tol              default"
   else
      write(outunit,"(SS,1P,A,E18.10E2)") "  hypre_BoomerAMG_Tol             ", solver_control%hypre_cntl%BoomerAMG_Tol
   endif
   if (solver_control%hypre_cntl%BoomerAMG_StrongThreshold == huge(0.0d0)) then
      write(outunit,"(A)") "  hypre_BoomerAMG_StrongThreshold  default"
   else
      write(outunit,"(SS,1P,A,E18.10E2)") "  hypre_BoomerAMG_StrongThreshold ", solver_control%hypre_cntl%BoomerAMG_StrongThreshold
   endif
   if (solver_control%hypre_cntl%BoomerAMG_MaxRowSum == huge(0.0d0)) then
      write(outunit,"(A)") "  hypre_BoomerAMG_MaxRowSum        default"
   else
      write(outunit,"(SS,1P,A,E18.10E2)") "  hypre_BoomerAMG_MaxRowSum       ", solver_control%hypre_cntl%BoomerAMG_MaxRowSum
   endif
   if (solver_control%hypre_cntl%BoomerAMG_CoarsenType == huge(0)) then
      write(outunit,"(A)") "  hypre_BoomerAMG_CoarsenType      default"
   else
      write(outunit,"(A,I11)") "  hypre_BoomerAMG_CoarsenType     ", solver_control%hypre_cntl%BoomerAMG_CoarsenType
   endif
   if (solver_control%hypre_cntl%BoomerAMG_MeasureType == huge(0)) then
      write(outunit,"(A)") "  hypre_BoomerAMG_MeasureType      default"
   else
      write(outunit,"(A,I11)") "  hypre_BoomerAMG_MeasureType     ", solver_control%hypre_cntl%BoomerAMG_MeasureType
   endif
   if (solver_control%hypre_cntl%BoomerAMG_CycleType == huge(0)) then
      write(outunit,"(A)") "  hypre_BoomerAMG_CycleType        default"
   else
      write(outunit,"(A,I11)") "  hypre_BoomerAMG_CycleType       ", solver_control%hypre_cntl%BoomerAMG_CycleType
   endif
   if (.not. associated(solver_control%hypre_cntl%BoomerAMG_NumGridSweeps)) then
      write(outunit,"(A)") "  hypre_BoomerAMG_NumGridSweeps    default"
   else
      write(outunit,"(A,100I11)") "  hypre_BoomerAMG_NumGridSweeps   ", &
          solver_control%hypre_cntl%BoomerAMG_NumGridSweeps(1: &
          min(100,size(solver_control%hypre_cntl%BoomerAMG_NumGridSweeps)))
   endif
   if (.not. associated(solver_control%hypre_cntl%BoomerAMG_GridRelaxType)) then
      write(outunit,"(A)") "  hypre_BoomerAMG_GridRelaxType    default"
   else
      write(outunit,"(A,100I11)") "  hypre_BoomerAMG_GridRelaxType   ", &
         solver_control%hypre_cntl%BoomerAMG_GridRelaxType(1: &
         min(100,size(solver_control%hypre_cntl%BoomerAMG_GridRelaxType)))
   endif
   if (.not. associated(solver_control%hypre_cntl%BoomerAMG_GridRelaxPoints)) then
      write(outunit,"(A)") "  hypre_BoomerAMG_GridRelaxPoints  default"
   else
      write(outunit,"(A,SS,1P,1000E18.10E2)") "  hypre_BoomerAMG_GridRelaxPoints ", &
          solver_control%hypre_cntl%BoomerAMG_GridRelaxPoints
   endif
   if (.not. associated(solver_control%hypre_cntl%BoomerAMG_RelaxWeight)) then
      write(outunit,"(A)") "  hypre_BoomerAMG_RelaxWeight      default"
   else
      write(outunit,"(SS,1P,A,100E18.10E2)") "  hypre_BoomerAMG_RelaxWeight     ", &
         solver_control%hypre_cntl%BoomerAMG_RelaxWeight(1: &
         min(100,size(solver_control%hypre_cntl%BoomerAMG_RelaxWeight)))
   endif
   if (solver_control%hypre_cntl%BoomerAMG_DebugFlag == huge(0)) then
      write(outunit,"(A)") "  hypre_BoomerAMG_DebugFlag        default"
   else
      write(outunit,"(A,I11)") "  hypre_BoomerAMG_DebugFlag        ",solver_control%hypre_cntl%BoomerAMG_DebugFlag
   endif
endif
if (solver_control%preconditioner == HYPRE_PARASAILS_PRECONDITION) then
   if (solver_control%hypre_cntl%ParaSails_thresh == huge(0.0d0)) then
      write(outunit,"(A)") "  hypre_ParaSails_thresh           default"
   else
      write(outunit,"(SS,1P,A,E18.10E2)") "  hypre_ParaSails_thresh          ", solver_control%hypre_cntl%ParaSails_thresh
   endif
   if (solver_control%hypre_cntl%ParaSails_nlevels == huge(0)) then
      write(outunit,"(A)") "  hypre_ParaSails_nlevels          default"
   else
      write(outunit,"(A,I11)") "  hypre_ParaSails_nlevels         ", solver_control%hypre_cntl%ParaSails_nlevels
   endif
   if (solver_control%hypre_cntl%ParaSails_filter == huge(0.0d0)) then
      write(outunit,"(A)") "  hypre_ParaSails_filter           default"
   else
      write(outunit,"(SS,1P,A,E18.10E2)") "  hypre_ParaSails_filter          ", solver_control%hypre_cntl%ParaSails_filter
   endif
   if (solver_control%hypre_cntl%ParaSails_sym == huge(0)) then
      write(outunit,"(A)") "  hypre_ParaSails_sym              default"
   else
      write(outunit,"(A,I11)") "  hypre_ParaSails_sym             ", solver_control%hypre_cntl%ParaSails_sym
   endif
   if (solver_control%hypre_cntl%ParaSails_loadbal == huge(0.0d0)) then
      write(outunit,"(A)") "  hypre_ParaSails_loadbal          default"
   else
      write(outunit,"(SS,1P,A,E18.10E2)") "  hypre_ParaSails_loadbal         ", solver_control%hypre_cntl%ParaSails_loadbal
   endif
   if (solver_control%hypre_cntl%ParaSails_reuse == huge(0)) then
      write(outunit,"(A)") "  hypre_ParaSails_reuse            default"
   else
      write(outunit,"(A,I11)") "  hypre_ParaSails_reuse           ", solver_control%hypre_cntl%ParaSails_reuse
   endif
   if (solver_control%hypre_cntl%ParaSails_logging == huge(0)) then
      write(outunit,"(A)") "  hypre_ParaSails_logging          default"
   else
      write(outunit,"(A,I11)") "  hypre_ParaSails_logging         ", solver_control%hypre_cntl%ParaSails_logging
   endif
endif

write(outunit,"(A)") "---------------------------------------------------------"

end subroutine print_header

!          -------------
subroutine print_trailer
!          -------------

!----------------------------------------------------
! This routine prints a trailer message
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

if (loc_print_trailer_who == NO_ONE) return
if ((my_proc(procs) == MASTER .and. loc_print_trailer_who == SLAVES) .or. &
    (my_proc(procs) /= MASTER .and. loc_print_trailer_who == MASTER)) return

write(outunit,"(A)")
write(outunit,"(A)") "---------------------------------------------------------"
write(outunit,"(A)") "phaml_solve_pde complete"
write(outunit,"(A)")
write(outunit,"(A,I11)") 'number of loops through phases = ',loop
write(outunit,"(A)") 'termination condition:'
select case (ierr)
   case (MAX_VERTICES_ACHIEVED)
      write(outunit,"(A)") '  maximum vertices achieved'
   case (MAX_ELEMENTS_ACHIEVED)
      write(outunit,"(A)") '  maximum elements achieved'
   case (MAX_EQUATIONS_ACHIEVED)
      write(outunit,"(A)") '  maximum equations achieved'
   case (MAX_LEV_ACHIEVED)
      write(outunit,"(A)") '  maximum refinement levels achieved'
   case (MAX_DEG_ACHIEVED)
      write(outunit,"(A)") '  maximum degree achieved'
   case (TOLERANCE_ACHIEVED)
      write(outunit,"(A)") '  tolerance achieved'
   case (MAX_LOOP_ACHIEVED)
      write(outunit,"(A)") '  maximum loops through phases achieved'
   case (DONE_BALANCE_ONLY)
      write(outunit,"(A)") '  one trip load balance'
   case (DONE_REFINE_ONLY)
      write(outunit,"(A)") '  one trip refinement'
   case (DONE_SOLVE_ONLY)
      write(outunit,"(A)") '  one trip solver'
   case (ENERGY_ERREST_ACHIEVED)
      write(outunit,"(A)") '  energy error estimate achieved'
   case (LINF_ERREST_ACHIEVED)
      write(outunit,"(A)") '  L infinity error estimate achieved'
   case (L2_ERREST_ACHIEVED)
      write(outunit,"(A)") '  L2 error estimate achieved'
   case (ALLOC_FAILED)
      write(outunit,"(A)") '  allocation failed'
   case (USER_INPUT_ERROR)
      write(outunit,"(A)") '  user input error'
   case (PHAML_INTERNAL_ERROR)
      write(outunit,"(A)") '  PHAML internal error'
   case (REFINEMENT_STALLED)
      write(outunit,"(A)") '  refinement stalled'
   case (UNCLASSIFIED_ERROR)
      write(outunit,"(A)") '  unclassified error'
   case (MULTIPLE_ERRORS)
      write(outunit,"(A)") '  different error codes on different processors'
   case default
      write(outunit,"(A)") '  unknown'
end select
write(outunit,"(A)") "---------------------------------------------------------"

end subroutine print_trailer

!          ------------
subroutine slaves_solve
!          ------------

!----------------------------------------------------
! In this routine, the master process sends a message to the slaves
! to call phaml_solve_pde, and provides the paramter list
!----------------------------------------------------

!----------------------------------------------------
! Local variables

integer :: i, proc, ni, nr, ipos, rpos, dim1, dim2
integer, allocatable :: send_int(:)
real(my_real), allocatable :: send_real(:)
!----------------------------------------------------

ni = 103 + FN_LEN
nr = 35

if (associated(solver_control%hypre_cntl%BoomerAMG_NumGridSweeps)) then
   ni = ni + 1 + size(solver_control%hypre_cntl%BoomerAMG_NumGridSweeps)
else
   ni = ni + 1
endif

if (associated(solver_control%hypre_cntl%BoomerAMG_GridRelaxType)) then
   ni = ni + 1 + size(solver_control%hypre_cntl%BoomerAMG_GridRelaxType)
else
   ni = ni + 1
endif

if (associated(solver_control%hypre_cntl%BoomerAMG_GridRelaxPoints)) then
   ni = ni + 2 + size(solver_control%hypre_cntl%BoomerAMG_GridRelaxPoints)
else
   ni = ni + 2
endif

if (associated(solver_control%hypre_cntl%BoomerAMG_RelaxWeight)) then
   ni = ni + 1
   nr = nr + size(solver_control%hypre_cntl%BoomerAMG_RelaxWeight)
else
   ni = ni + 1
endif

allocate(send_int(ni), send_real(nr), stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in slaves_solve",procs=procs)
   stop
endif

send_int( 1) = 1 ! signal to call phaml_solve_pde
send_int( 2) = refine_control%max_elem
send_int( 3) = refine_control%max_vert
send_int( 4) = refine_control%max_lev
send_int( 5) = loc_max_refsolveloop
send_int( 6) = loc_task
if (warn) then
   send_int( 7) = 1
else
   send_int( 7) = 0
endif
send_int( 8) = phaml_solution%grid%errtype
send_int( 9) = refine_control%t3s_nunif
send_int(10) = io_control%print_grid_when
send_int(11) = io_control%print_grid_who
send_int(12) = io_control%print_error_when
send_int(13) = io_control%print_error_who
send_int(14) = io_control%print_time_when
send_int(15) = io_control%print_time_who
send_int(16) = loc_print_header_who
send_int(17) = loc_print_trailer_who
send_int(18) = loc_clocks
send_int(19) = io_control%draw_grid_when
send_int(20) = solver_control%eigensolver
if (io_control%pause_after_draw) then
   send_int(21) = 1
else
   send_int(21) = 0
endif
send_int(22) = 0
if (present(pause_after_phases)) then; if (pause_after_phases) then
   send_int(22) = 1
endif; endif
send_int(23) = 0
if (present(pause_at_start)) then; if (pause_at_start) then
   send_int(23) = 1
endif; endif
send_int(24) = 0
if (present(pause_at_end)) then; if (pause_at_end) then
   send_int(24) = 1
endif; endif
send_int(25) = refine_control%reftype
send_int(26) = refine_control%hp_strategy
send_int(27) = loc_sequential_vert
send_int(28) = refine_control%error_estimator
send_int(29) = refine_control%refterm
if (refine_control%derefine) then
   send_int(30) = 1
else
   send_int(30) = 0
endif
send_int(31) = loc_partition_method
send_int(32) = solver_control%mg_comm
send_int(33) = solver_control%solver
send_int(34) = solver_control%preconditioner
send_int(35) = solver_control%ncycle
send_int(36) = solver_control%prerelax
send_int(37) = solver_control%postrelax
send_int(38) = solver_control%dd_iterations
if (solver_control%ignore_quad_err) then
   send_int(39) = 1
else
   send_int(39) = 0
endif
send_int(40) = solver_control%krylov_iter
send_int(41) = solver_control%krylov_restart
send_int(42) = solver_control%num_eval
send_int(43) = solver_control%arpack_cntl%ncv
if (phaml_solution%master_draws_grid) then
   send_int(44) = 1
else
   send_int(44) = 0
endif
if (phaml_solution%i_draw_grid) then
   send_int(45) = 1
else
   send_int(45) = 0
endif
send_int(46) = refine_control%t3s_maxref
send_int(47) = refine_control%t3s_maxdeginc
send_int(48) = loc_print_eval_when
send_int(49) = loc_print_eval_who
if (solver_control%petsc_matrix_free) then
   send_int(50) = 1
else
   send_int(50) = 0
endif
send_real( 1) = refine_control%inc_factor
send_real( 2) = solver_control%lambda0

send_int(51) = solver_control%hypre_cntl%BoomerAMG_MaxLevels
send_int(52) = solver_control%hypre_cntl%BoomerAMG_MaxIter
send_int(53) = solver_control%hypre_cntl%BoomerAMG_CoarsenType
send_int(54) = solver_control%hypre_cntl%BoomerAMG_MeasureType
send_int(55) = solver_control%hypre_cntl%BoomerAMG_CycleType
send_int(56) = solver_control%blopex_cntl%maxit
send_int(57) = solver_control%hypre_cntl%BoomerAMG_DebugFlag
send_int(58) = solver_control%hypre_cntl%ParaSails_nlevels
send_int(59) = solver_control%hypre_cntl%ParaSails_sym
send_int(60) = solver_control%hypre_cntl%ParaSails_reuse
send_int(61) = solver_control%hypre_cntl%ParaSails_logging
send_int(62) = solver_control%hypre_cntl%PCG_MaxIter
send_int(63) = solver_control%hypre_cntl%PCG_TwoNorm
send_int(64) = solver_control%hypre_cntl%PCG_RelChange
send_int(65) = solver_control%hypre_cntl%PCG_Logging
send_int(66) = solver_control%hypre_cntl%GMRES_KDim
send_int(67) = solver_control%hypre_cntl%GMRES_MaxIter
send_int(68) = solver_control%hypre_cntl%GMRES_Logging
send_int(69) = solver_control%arpack_cntl%maxit
send_real(3) = solver_control%hypre_cntl%BoomerAMG_Tol
send_real(4) = solver_control%hypre_cntl%BoomerAMG_StrongThreshold
send_real(5) = solver_control%hypre_cntl%BoomerAMG_MaxRowSum
send_real(6) = solver_control%hypre_cntl%ParaSails_thresh
send_real(7) = solver_control%hypre_cntl%ParaSails_filter
send_real(8) = solver_control%hypre_cntl%ParaSails_loadbal
send_real(9) = solver_control%hypre_cntl%PCG_Tol
send_real(10) = solver_control%hypre_cntl%GMRES_Tol
send_real(11) = solver_control%arpack_cntl%tol

send_int(70) = solver_control%petsc_cntl%petsc_gmres_max_steps
send_int(71) = solver_control%petsc_cntl%petsc_maxits
send_int(72) = solver_control%petsc_cntl%petsc_ilu_levels
send_int(73) = solver_control%petsc_cntl%petsc_icc_levels
send_int(74) = solver_control%petsc_cntl%petsc_ilu_maxrowcount
send_int(75) = solver_control%petsc_cntl%petsc_sor_its
send_int(76) = solver_control%petsc_cntl%petsc_sor_lits
send_int(77) = solver_control%petsc_cntl%petsc_asm_overlap
if (solver_control%petsc_cntl%petsc_eisenstat_nodiagscaling) then
   send_int(78) = 1
else
   send_int(78) = 0
endif
send_real(12) = solver_control%petsc_cntl%petsc_richardson_damping_factor
send_real(13) = solver_control%petsc_cntl%petsc_chebychev_emin
send_real(14) = solver_control%petsc_cntl%petsc_chebychev_emax
send_real(15) = solver_control%petsc_cntl%petsc_rtol
send_real(16) = solver_control%petsc_cntl%petsc_atol
send_real(17) = solver_control%petsc_cntl%petsc_dtol
send_real(18) = solver_control%petsc_cntl%petsc_ilu_dt
send_real(19) = solver_control%petsc_cntl%petsc_sor_omega
send_real(20) = solver_control%petsc_cntl%petsc_eisenstat_omega
send_real(21) = solver_control%petsc_cntl%petsc_ilu_dtcol
send_real(22) = solver_control%mg_tol
send_real(23) = refine_control%term_energy_err
send_real(24) = refine_control%tp_gamma
send_real(25) = refine_control%term_Linf_err
send_real(26) = refine_control%term_L2_err
send_real(27) = refine_control%t3s_gamma
send_real(28) = refine_control%t3s_eta
send_real(29) = refine_control%reftol
send_real(30) = solver_control%krylov_tol
send_real(31) = solver_control%blopex_cntl%atol
send_real(32) = solver_control%blopex_cntl%rtol
send_real(33) = refine_control%sp_gamma_h
send_real(34) = refine_control%sp_gamma_p
send_real(35) = refine_control%refsoln_pbias

send_int(79) = solver_control%lambda0_side
send_int(80) = solver_control%coarse_size
send_int(81) = solver_control%coarse_method
send_int(82) = solver_control%scale_evec
send_int(83) = loc_prebalance
send_int(84) = loc_postbalance
send_int(85) = loc_degree
send_int(86) = io_control%print_linsys_when
send_int(87) = io_control%print_linsys_who
send_int(88) = solver_control%inc_quad_order
send_int(89) = io_control%print_error_what
send_int(90) = io_control%print_errest_what
send_int(91) = refine_control%max_dof
send_int(92) = refine_control%max_deg
send_int(93) = solver_control%prerelax_ho
send_int(94) = solver_control%postrelax_ho
send_int(95) = solver_control%transformation
if (loc_solve_init) then
   send_int(96) = 1
else
   send_int(96) = 0
endif
if (loc_stop_on_maxlev) then
   send_int(97) = 1
else
   send_int(97) = 0
endif
if (loc_stop_on_maxdeg) then
   send_int(98) = 1
else
   send_int(98) = 0
endif
send_int(99) = refine_control%edge_rule
send_int(100) = refine_control%nlp_max_h_dec
send_int(101) = refine_control%nlp_max_h_inc
send_int(102) = refine_control%nlp_max_p_dec
send_int(103) = refine_control%nlp_max_p_inc

do i=1,FN_LEN
   send_int(103+i:103+i) = ichar(phaml_solution%zoltan_param_file(i:i))
end do

ipos = 103 + FN_LEN
rpos = 35

if (associated(solver_control%hypre_cntl%BoomerAMG_NumGridSweeps)) then
   dim1 = size(solver_control%hypre_cntl%BoomerAMG_NumGridSweeps)
   send_int(ipos+1) = dim1
   send_int(ipos+2:ipos+1+dim1) = solver_control%hypre_cntl%BoomerAMG_NumGridSweeps
   ipos = ipos + 1 + dim1
else
   send_int(ipos+1) = 0
   ipos = ipos + 1
endif

if (associated(solver_control%hypre_cntl%BoomerAMG_GridRelaxType)) then
   dim1 = size(solver_control%hypre_cntl%BoomerAMG_GridRelaxType)
   send_int(ipos+1) = dim1
   send_int(ipos+2:ipos+1+dim1) = solver_control%hypre_cntl%BoomerAMG_GridRelaxType
   ipos = ipos + 1 + dim1
else
   send_int(ipos+1) = 0
   ipos = ipos + 1
endif

if (associated(solver_control%hypre_cntl%BoomerAMG_GridRelaxPoints)) then
   dim1 = size(solver_control%hypre_cntl%BoomerAMG_GridRelaxPoints,1)
   dim2 = size(solver_control%hypre_cntl%BoomerAMG_GridRelaxPoints,2)
   send_int(ipos+1) = dim1
   send_int(ipos+2) = dim2
   send_int(ipos+3:ipos+2+dim1*dim2) = &
      reshape(solver_control%hypre_cntl%BoomerAMG_GridRelaxType,(/dim1*dim2/))
   ipos = ipos + 2 + dim1*dim2
else
   send_int(ipos+1) = 0
   send_int(ipos+2) = 0
   ipos = ipos + 2
endif

if (associated(solver_control%hypre_cntl%BoomerAMG_RelaxWeight)) then
   dim1 = size(solver_control%hypre_cntl%BoomerAMG_RelaxWeight)
   send_int(ipos+1) = dim1
   send_real(rpos+1:rpos+dim1) = solver_control%hypre_cntl%BoomerAMG_RelaxWeight
   ipos = ipos + 1
   rpos = rpos + dim1
else
   send_int(ipos+1) = 0
   ipos = ipos + 1
endif

do proc=1,num_proc(procs)
   call phaml_send(procs,proc,send_int,ni,send_real,nr,101)
end do

deallocate(send_int, send_real, stat=astat)

return
end subroutine slaves_solve

end subroutine phaml_solve_pde

!          ---------
subroutine terminate(phaml_solution,procs,solver_control,still_sequential, &
                     pause_at_end)
!          ---------

!----------------------------------------------------
! Dummy arguments:

type(phaml_solution_type), intent(inout) :: phaml_solution
type(proc_info), intent(in) :: procs
type (solver_options), intent(inout) :: solver_control
logical, intent(in) :: still_sequential
logical, intent(in), optional :: pause_at_end

!----------------------------------------------------
! This routine closes up everything at the end
!----------------------------------------------------

phaml_solution%still_sequential = still_sequential

call destroy_watchgroup(all_watches)
call destroy_watch(ttotal)
call destroy_watch(trefine)
call destroy_watch(trecon)
call destroy_watch(tpartition)
call destroy_watch(tdistribute)
call destroy_watch(tassemble)
call destroy_watch(tsolve)
call destroy_watch(ptotal)
call destroy_watch(prefine)
call destroy_watch(precon)
call destroy_watch(ppartition)
call destroy_watch(pdistribute)
call destroy_watch(passemble)
call destroy_watch(psolve)
call destroy_watch(ctrecon)
call destroy_watch(ctpartition)
call destroy_watch(ctdistribute)
call destroy_watch(ctassemble)
call destroy_watch(ctsolve)
call destroy_watch(cprecon)
call destroy_watch(cppartition)
call destroy_watch(cpdistribute)
call destroy_watch(cpassemble)
call destroy_watch(cpsolve)

! allocated options

if (associated(solver_control%hypre_cntl%BoomerAMG_NumGridSweeps)) &
   deallocate(solver_control%hypre_cntl%BoomerAMG_NumGridSweeps)

if (associated(solver_control%hypre_cntl%BoomerAMG_GridRelaxType)) &
   deallocate(solver_control%hypre_cntl%BoomerAMG_GridRelaxType)

if (associated(solver_control%hypre_cntl%BoomerAMG_GridRelaxPoints)) &
   deallocate(solver_control%hypre_cntl%BoomerAMG_GridRelaxPoints)

if (associated(solver_control%hypre_cntl%BoomerAMG_RelaxWeight)) &
   deallocate(solver_control%hypre_cntl%BoomerAMG_RelaxWeight)

! pause, if requested

if (present(pause_at_end)) then
   call pause_until_enter(procs,pause_at_end,dont_pausewatch=.true.)
endif

end subroutine terminate

!          -----------
subroutine phaml_slave
!          -----------

!----------------------------------------------------
! This is the main program for a slave process.  It sits in a loop waiting
! for commands from the master process, until it is told to terminate.
!----------------------------------------------------
implicit none

!----------------------------------------------------
! Local variables:

type (phaml_solution_type) :: my_pde
integer :: iterm, max_elem, max_vert, max_eq, max_lev, max_deg, &
           print_grid_when, print_error_when, print_time_when, print_eval_when,&
           print_linsys_when, print_linsys_who, &
           print_grid_who, print_error_who, print_time_who, print_eval_who, &
           print_header_who, print_trailer_who, &
           print_error_what, print_errest_what, &
           clocks, draw_grid_when, sequential_vert, &
           job, ni, nr, proc, error_estimator, max_refsolveloop, &
           partition_method, task, first_int, first_real, pde_index, &
           prebalance, postbalance, &
           solver, preconditioner, mg_cycles, mg_prerelax, mg_postrelax, &
           mg_prerelax_ho, mg_postrelax_ho, dd_iterations, errtype, edge_rule, &
           reftype, refterm, hp_strategy, t3s_nunif, t3s_maxref, t3s_maxdeginc,&
           nlp_max_h_dec, nlp_max_h_inc, nlp_max_p_dec, nlp_max_p_inc, &
           eigensolver, num_eval, arpack_ncv, arpack_maxit, blopex_maxit, &
           astat, kernel, comp1, comp2, eigen1, eigen2, p, q, coarse_size, &
           coarse_method, degree, inc_quad_order, scale_evec, krylov_iter, &
           krylov_restart, mg_comm, lambda0_side, transformation
integer, pointer :: recv_int(:)
real (my_real), pointer :: recv_real(:)
logical :: pause_after_draw, pause_at_start, pause_at_end, &
           pause_after_phases, petsc_matrix_free, ignore_quad_err, derefine, &
           solve_init, stop_on_maxlev, stop_on_maxdeg
real (my_real) :: term_energy_err, term_Linf_err, term_L2_err, inc_factor, &
                  lambda0, arpack_tol, blopex_atol, blopex_rtol, dum, mg_tol, &
                  reftol, krylov_tol, t3s_gamma, t3s_eta, tp_gamma, &
                  sp_gamma_h, sp_gamma_p, refsoln_pbias
character(len=FN_LEN) :: zoltan_param_file
integer :: hypre_BoomerAMG_MaxLevels, &
   hypre_BoomerAMG_MaxIter,hypre_BoomerAMG_CoarsenType, &
   hypre_BoomerAMG_MeasureType,hypre_BoomerAMG_CycleType, &
   hypre_BoomerAMG_DebugFlag,hypre_ParaSails_nlevels, &
   hypre_ParaSails_sym,hypre_ParaSails_reuse,hypre_ParaSails_logging, &
   hypre_PCG_MaxIter,hypre_PCG_TwoNorm,hypre_PCG_RelChange, &
   hypre_PCG_Logging,hypre_GMRES_KDim,hypre_GMRES_MaxIter,hypre_GMRES_Logging
real(my_real) :: hypre_BoomerAMG_Tol, &
   hypre_BoomerAMG_StrongThreshold,hypre_BoomerAMG_MaxRowSum, &
   hypre_ParaSails_thresh,hypre_ParaSails_filter,hypre_ParaSails_loadbal, &
   hypre_PCG_Tol,hypre_GMRES_Tol
integer, allocatable :: hypre_BoomerAMG_NumGridSweeps(:), &
   hypre_BoomerAMG_GridRelaxType(:),hypre_BoomerAMG_GridRelaxPoints(:,:)
real(my_real), allocatable :: hypre_BoomerAMG_RelaxWeight(:)
real(my_real) :: petsc_richardson_damping_factor, &
   petsc_chebychev_emin, petsc_chebychev_emax, petsc_rtol, petsc_atol,  &
   petsc_dtol, petsc_ilu_dt, petsc_ilu_dtcol, petsc_sor_omega, &
   petsc_eisenstat_omega
integer :: petsc_gmres_max_steps, petsc_maxits,    &
   petsc_ilu_levels, petsc_icc_levels, petsc_ilu_maxrowcount, petsc_sor_its, &
   petsc_sor_lits, petsc_asm_overlap
logical :: petsc_eisenstat_nodiagscaling


logical :: hypre_has_NumGridSweeps, hypre_has_GridRelaxType, &
           hypre_has_GridRelaxPoints, hypre_has_RelaxWeight

integer :: dim1, dim2, ipos, rpos, i
type(proc_info) :: invoker_procs
character(len=11) :: form

!----------------------------------------------------
! Begin executable code

! create the solution

call phaml_create(my_pde)

! receive job assignments and loop until the assignment is to terminate

do

! receive the message

   call phaml_recv(my_pde%procs,proc,recv_int,ni,recv_real,nr,101)
   job = recv_int(1)

! select which job to do

   select case (job)

   case (1) ! solve the pde

! unpack the parameters

      max_elem          = recv_int( 2)
      max_vert          = recv_int( 3)
      max_lev           = recv_int( 4)
      max_refsolveloop  = recv_int( 5)
      task              = recv_int( 6)
      warn              = (recv_int( 7)==1)
      errtype           = recv_int( 8)
      t3s_nunif         = recv_int( 9)
      print_grid_when   = recv_int(10)
      print_grid_who    = recv_int(11)
      print_error_when  = recv_int(12)
      print_error_who   = recv_int(13)
      print_time_when   = recv_int(14)
      print_time_who    = recv_int(15)
      print_header_who  = recv_int(16)
      print_trailer_who = recv_int(17)
      clocks            = recv_int(18)
      draw_grid_when    = recv_int(19)
      eigensolver       = recv_int(20)
      pause_after_draw  =(recv_int(21)==1)
      pause_after_phases=(recv_int(22)==1)
      pause_at_start    =(recv_int(23)==1)
      pause_at_end      =(recv_int(24)==1)
      reftype           = recv_int(25)
      hp_strategy       = recv_int(26)
      sequential_vert   = recv_int(27)
      error_estimator   = recv_int(28)
      refterm           = recv_int(29)
      derefine          =(recv_int(30)==1)
      partition_method  = recv_int(31)
      mg_comm           = recv_int(32)
      solver            = recv_int(33)
      preconditioner    = recv_int(34)
      mg_cycles         = recv_int(35)
      mg_prerelax       = recv_int(36)
      mg_postrelax      = recv_int(37)
      dd_iterations     = recv_int(38)
      ignore_quad_err   =(recv_int(39)==1)
      krylov_iter       = recv_int(40)
      krylov_restart    = recv_int(41)
      num_eval          = recv_int(42)
      arpack_ncv        = recv_int(43)
      my_pde%master_draws_grid    = (recv_int(44)==1)
      my_pde%i_draw_grid          = (recv_int(45)==1)
      t3s_maxref        = recv_int(46)
      t3s_maxdeginc     = recv_int(47)
      print_eval_when   = recv_int(48)
      print_eval_who    = recv_int(49)
      petsc_matrix_free = (recv_int(50)==1)
      inc_factor = recv_real( 1)
      lambda0    = recv_real( 2)

      hypre_BoomerAMG_MaxLevels   = recv_int(51)
      hypre_BoomerAMG_MaxIter     = recv_int(52)
      hypre_BoomerAMG_CoarsenType = recv_int(53)
      hypre_BoomerAMG_MeasureType = recv_int(54)
      hypre_BoomerAMG_CycleType   = recv_int(55)
      blopex_maxit                = recv_int(56)
      hypre_BoomerAMG_DebugFlag   = recv_int(57)
      hypre_ParaSails_nlevels     = recv_int(58)
      hypre_ParaSails_sym         = recv_int(59)
      hypre_ParaSails_reuse       = recv_int(60)
      hypre_ParaSails_logging     = recv_int(61)
      hypre_PCG_MaxIter           = recv_int(62)
      hypre_PCG_TwoNorm           = recv_int(63)
      hypre_PCG_RelChange         = recv_int(64)
      hypre_PCG_Logging           = recv_int(65)
      hypre_GMRES_KDim            = recv_int(66)
      hypre_GMRES_MaxIter         = recv_int(67)
      hypre_GMRES_Logging         = recv_int(68)
      arpack_maxit                = recv_int(69)
      hypre_BoomerAMG_Tol             = recv_real(3)
      hypre_BoomerAMG_StrongThreshold = recv_real(4)
      hypre_BoomerAMG_MaxRowSum       = recv_real(5)
      hypre_ParaSails_thresh          = recv_real(6)
      hypre_ParaSails_filter          = recv_real(7)
      hypre_ParaSails_loadbal         = recv_real(8)
      hypre_PCG_Tol                   = recv_real(9)
      hypre_GMRES_Tol                 = recv_real(10)
      arpack_tol                      = recv_real(11)

      petsc_gmres_max_steps           = recv_int(70)
      petsc_maxits                    = recv_int(71)
      petsc_ilu_levels                = recv_int(72)
      petsc_icc_levels                = recv_int(73)
      petsc_ilu_maxrowcount           = recv_int(74)
      petsc_sor_its                   = recv_int(75)
      petsc_sor_lits                  = recv_int(76)
      petsc_asm_overlap               = recv_int(77)
      petsc_eisenstat_nodiagscaling   = recv_int(78)==1
      petsc_richardson_damping_factor = recv_real(12)
      petsc_chebychev_emin            = recv_real(13)
      petsc_chebychev_emax            = recv_real(14)
      petsc_rtol                      = recv_real(15)
      petsc_atol                      = recv_real(16)
      petsc_dtol                      = recv_real(17)
      petsc_ilu_dt                    = recv_real(18)
      petsc_sor_omega                 = recv_real(19)
      petsc_eisenstat_omega           = recv_real(20)
      petsc_ilu_dtcol                 = recv_real(21)
      mg_tol                          = recv_real(22)
      term_energy_err                 = recv_real(23)
      tp_gamma                        = recv_real(24)
      term_Linf_err                   = recv_real(25)
      term_L2_err                     = recv_real(26)
      t3s_gamma                       = recv_real(27)
      t3s_eta                         = recv_real(28)
      reftol                          = recv_real(29)
      krylov_tol                      = recv_real(30)
      blopex_atol                     = recv_real(31)
      blopex_rtol                     = recv_real(32)
      sp_gamma_h                      = recv_real(33)
      sp_gamma_p                      = recv_real(34)
      refsoln_pbias                   = recv_real(35)

      lambda0_side                    = recv_int(79)
      coarse_size                     = recv_int(80)
      coarse_method                   = recv_int(81)
      scale_evec                      = recv_int(82)
      prebalance                      = recv_int(83)
      postbalance                     = recv_int(84)
      degree                          = recv_int(85)
      print_linsys_when               = recv_int(86)
      print_linsys_who                = recv_int(87)
      inc_quad_order                  = recv_int(88)
      print_error_what                = recv_int(89)
      print_errest_what               = recv_int(90)
      max_eq                          = recv_int(91)
      max_deg                         = recv_int(92)
      mg_prerelax_ho                  = recv_int(93)
      mg_postrelax_ho                 = recv_int(94)
      transformation                  = recv_int(95)
      solve_init                      = (recv_int(96)==1)
      stop_on_maxlev                  = (recv_int(97)==1)
      stop_on_maxdeg                  = (recv_int(98)==1)
      edge_rule                       = recv_int(99)
      nlp_max_h_dec                   = recv_int(100)
      nlp_max_h_inc                   = recv_int(101)
      nlp_max_p_dec                   = recv_int(102)
      nlp_max_p_inc                   = recv_int(103)

      do i=1,FN_LEN
         zoltan_param_file(i:i) = char(recv_int(103+i))
      end do

      ipos = 103 + FN_LEN
      rpos = 35

      dim1 = recv_int(ipos+1)
      if (dim1 /= 0) then
         allocate(hypre_BoomerAMG_NumGridSweeps(dim1),stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("allocation failed in phaml_slave main routine",procs=my_pde%procs)
            stop
         endif
         hypre_BoomerAMG_NumGridSweeps = recv_int(ipos+2:ipos+1+dim1)
         ipos = ipos + 1 + dim1
         hypre_has_NumGridSweeps = .true.
      else
         allocate(hypre_BoomerAMG_NumGridSweeps(1))
         ipos = ipos + 1
         hypre_has_NumGridSweeps = .false.
      endif

      dim1 = recv_int(ipos+1)
      if (dim1 /= 0) then
         allocate(hypre_BoomerAMG_GridRelaxType(dim1),stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("allocation failed in phaml_slave main routine",procs=my_pde%procs)
            stop
         endif
         hypre_BoomerAMG_GridRelaxType = recv_int(ipos+2:ipos+1+dim1)
         ipos = ipos + 1 + dim1
         hypre_has_GridRelaxType = .true.
      else
         allocate(hypre_BoomerAMG_GridRelaxType(1))
         ipos = ipos + 1
         hypre_has_GridRelaxType = .false.
      endif

      dim1 = recv_int(ipos+1)
      if (dim1 /= 0) then
         dim2 = recv_int(ipos+2)
         allocate(hypre_BoomerAMG_GridRelaxPoints(dim1,dim2),stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("allocation failed in phaml_slave main routine",procs=my_pde%procs)
            stop
         endif
         hypre_BoomerAMG_GridRelaxPoints = &
            reshape(recv_int(ipos+3:ipos+2+dim1*dim2),(/dim1,dim2/))
         ipos = ipos + 2 + dim1*dim2
         hypre_has_GridRelaxPoints = .true.
      else
         allocate(hypre_BoomerAMG_GridRelaxPoints(1,1))
         ipos = ipos + 2
         hypre_has_GridRelaxPoints = .false.
      endif

      dim1 = recv_int(ipos+1)
      if (dim1 /= 0) then
         allocate(hypre_BoomerAMG_RelaxWeight(dim1),stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("allocation failed in phaml_slave main routine",procs=my_pde%procs)
            stop
         endif
         hypre_BoomerAMG_RelaxWeight = recv_real(rpos+1:rpos+dim1)
         ipos = ipos + 1
         rpos = rpos + dim1
         hypre_has_RelaxWeight = .true.
      else
         allocate(hypre_BoomerAMG_RelaxWeight(1))
         ipos = ipos + 1
         hypre_has_RelaxWeight = .false.
      endif

      if (ni > 0) deallocate(recv_int,stat=astat)
      if (nr > 0) deallocate(recv_real,stat=astat)

      call phaml_solve_pde(my_pde, iterm, max_elem, max_vert, max_eq,        &
                     max_lev, max_deg, stop_on_maxlev, stop_on_maxdeg,       &
                     max_refsolveloop, term_energy_err,                      &
                     term_Linf_err, term_L2_err, task,                       &
                     print_grid_when, print_grid_who, print_error_when,      &
                     print_error_who, print_error_what, print_errest_what,   &
                     print_linsys_when, print_linsys_who, print_time_when,   &
                     print_time_who, print_eval_when, print_eval_who,        &
                     print_header_who, print_trailer_who, warn, clocks,      &
                     draw_grid_when, pause_after_draw, pause_after_phases,   &
                     pause_at_start, pause_at_end, solve_init,               &
                     sequential_vert, inc_factor, error_estimator, errtype,  &
                     reftype, refterm, reftol, hp_strategy, t3s_gamma,       &
                     t3s_eta, t3s_nunif, t3s_maxref, t3s_maxdeginc, tp_gamma,&
                     sp_gamma_h, sp_gamma_p, nlp_max_h_dec, nlp_max_h_inc,   &
                     nlp_max_p_dec, nlp_max_p_inc, refsoln_pbias, derefine,  &
                     partition_method, edge_rule, zoltan_param_file,         &
                     prebalance, postbalance, petsc_matrix_free,             &
                     solver, preconditioner, mg_cycles, mg_tol, mg_prerelax, &
                     mg_postrelax, mg_prerelax_ho, mg_postrelax_ho,          &
                     dd_iterations, krylov_iter, krylov_restart, krylov_tol, &
                     mg_comm, ignore_quad_err, eigensolver, num_eval,        &
                     lambda0, lambda0_side, transformation, scale_evec,      &
                     arpack_ncv, arpack_maxit, arpack_tol,                   &
                     blopex_maxit, blopex_atol, blopex_rtol,                 &
                     degree, inc_quad_order,                                 &
   hypre_BoomerAMG_MaxLevels,hypre_BoomerAMG_MaxIter,hypre_BoomerAMG_Tol,     &
   hypre_BoomerAMG_StrongThreshold,hypre_BoomerAMG_MaxRowSum,                 &
   hypre_BoomerAMG_CoarsenType,hypre_BoomerAMG_MeasureType,                   &
   hypre_BoomerAMG_CycleType,hypre_BoomerAMG_NumGridSweeps,                   &
   hypre_BoomerAMG_GridRelaxType,hypre_BoomerAMG_GridRelaxPoints,             &
   hypre_BoomerAMG_RelaxWeight,                                               &
   hypre_BoomerAMG_DebugFlag,hypre_ParaSails_thresh,hypre_ParaSails_nlevels,  &
   hypre_ParaSails_filter,hypre_ParaSails_sym,hypre_ParaSails_loadbal,        &
   hypre_ParaSails_reuse,hypre_ParaSails_logging,hypre_PCG_Tol,               &
   hypre_PCG_MaxIter, hypre_PCG_TwoNorm,hypre_PCG_RelChange,hypre_PCG_Logging,&
   hypre_GMRES_KDim, hypre_GMRES_Tol,hypre_GMRES_MaxIter,hypre_GMRES_Logging, &
   petsc_richardson_damping_factor, petsc_chebychev_emin,                     &
   petsc_chebychev_emax, petsc_gmres_max_steps, petsc_rtol, petsc_atol,       &
   petsc_dtol, petsc_maxits, petsc_ilu_levels, petsc_icc_levels, petsc_ilu_dt,&
   petsc_ilu_dtcol, petsc_ilu_maxrowcount, petsc_sor_omega, petsc_sor_its,    &
   petsc_sor_lits, petsc_eisenstat_nodiagscaling, petsc_eisenstat_omega,      &
   petsc_asm_overlap,coarse_size,coarse_method,                               &
   hypre_has_NumGridSweeps,hypre_has_GridRelaxType,hypre_has_GridRelaxPoints, &
   hypre_has_RelaxWeight)

   deallocate(hypre_BoomerAMG_NumGridSweeps, hypre_BoomerAMG_GridRelaxType, &
              hypre_BoomerAMG_GridRelaxPoints, hypre_BoomerAMG_RelaxWeight, &
              stat=astat)

   case (2) ! terminate

      if (ni > 0) deallocate(recv_int,stat=astat)
      if (nr > 0) deallocate(recv_real,stat=astat)
      call phaml_destroy(my_pde)
      exit

   case (3) ! evaluate the solution

      first_int = 9
      first_real = 1
      call unpack_procs(invoker_procs,recv_int,first_int,recv_real,first_real)
      first_real = first_real-1
      nr = nr - first_real
      call evaluate_slave(my_pde,recv_real(first_real+1:first_real+nr/2),&
                          recv_real(first_real+1+nr/2:first_real+nr), &
                          invoker_procs,recv_int(2),recv_int(3),recv_int(4:8))
      if (ni > 0) deallocate(recv_int,stat=astat)
      if (nr+first_real > 0) deallocate(recv_real,stat=astat)
      call cleanup_unpack_procs(invoker_procs)

   case (4) ! connect to another pde in pde(:)

      pde_index = recv_int(2)
      first_int = 3
      first_real = 1
      if (nr == 0) allocate(recv_real(1),stat=astat) ! avoid passing unassociated pointer
      if (astat /= 0) then
         call warning("unnecessary allocation of recv_real failed in phaml_slave")
      endif
      call unpack_procs(pde(pde_index)%procs,recv_int,first_int,recv_real, &
                        first_real)
      if (ni > 0) deallocate(recv_int,stat=astat)
      deallocate(recv_real,stat=astat)

   case (5) ! save phaml_solution

      call phaml_store(my_pde,recv_int(2))
      if (ni > 0) deallocate(recv_int,stat=astat)

   case (6) ! restore phaml_solution

      call phaml_restore(my_pde,recv_int(2),recv_int(3)==1,recv_int(4)==1)
      if (ni > 0) deallocate(recv_int,stat=astat)

   case (7) ! open a file

      if (recv_int(3) == 0) then
         form = "FORMATTED"
      else
         form = "UNFORMATTED"
      endif
      call phaml_popen(my_pde,recv_int(2),"",form)
      if (ni > 0) deallocate(recv_int,stat=astat)

   case (8) ! close a file

      call phaml_pclose(my_pde,recv_int(2))
      if (ni > 0) deallocate(recv_int,stat=astat)

   case (9) ! update user module variables

      call update_usermod(my_pde)
      if (ni > 0) deallocate(recv_int,stat=astat)

   case (10) ! increase process universe

      call increase_universe(my_pde%procs,recv_int)
      if (ni > 0) deallocate(recv_int,stat=astat)

   case (11) ! decrease process universe

      call decrease_universe(my_pde%procs)
      if (ni > 0) deallocate(recv_int,stat=astat)

   case (12) ! compute an integral of the solution

      kernel = recv_int(2)
      comp1 = recv_int(3)
      eigen1 = recv_int(4)
      comp2 = recv_int(5)
      eigen2 = recv_int(6)
      p = recv_int(7)
      q = recv_int(8)
      dum = phaml_integrate(my_pde,kernel,comp1,eigen1,comp2,eigen2,p,q)
      if (ni > 0) deallocate(recv_int,stat=astat)

   case (13) ! scale the computed solution

      call phaml_scale(my_pde,factor=recv_real(1),comp=recv_int(2), &
                       eigen=recv_int(3))
      if (ni > 0) deallocate(recv_int,stat=astat)
      if (nr > 0) deallocate(recv_real,stat=astat)

   case (14) ! query

      call phaml_query(my_pde,nvert_proc=recv_int)
      if (ni > 0) deallocate(recv_int,stat=astat)
      if (nr > 0) deallocate(recv_real,stat=astat)

   case (15) ! compress

      call phaml_compress(my_pde)
      if (ni > 0) deallocate(recv_int,stat=astat)
      if (nr > 0) deallocate(recv_real,stat=astat)

   case (16) ! copy solution to oldsoln

      call phaml_copy_soln_to_old(my_pde)
      if (ni > 0) deallocate(recv_int,stat=astat)
      if (nr > 0) deallocate(recv_real,stat=astat)

   case (17) ! store matrix

      call store_matrix(my_pde%grid,my_pde%procs,my_pde%still_sequential, &
                        my_pde%system_size,my_pde%eq_type,recv_int(2))
      if (ni > 0) deallocate(recv_int,stat=astat)
      if (nr > 0) deallocate(recv_real,stat=astat)

   end select
end do

end subroutine phaml_slave

!          ----------------
subroutine print_time_info(procs,io_control,this_time,which,tag)
!          ----------------

!----------------------------------------------------
! This routine prints execution times
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type (proc_info), intent(in) :: procs
type(io_options), intent(in) :: io_control
integer, intent(in) :: this_time(:)
character(len=*), intent(in) :: which
integer, intent(in) :: tag
!----------------------------------------------------

!----------------------------------------------------
! Local variables:

character(len=4) :: default_clocks(4)
integer :: nclock,proc,i,who,when,astat
real (my_real), allocatable :: my_times(:,:),all_times(:,:,:)
integer :: no_ints(1),ni,nr
integer, pointer :: recv_int(:)
real (my_real), pointer :: recv_real(:)
character(len=13) :: fmt
real, pointer :: get_times(:,:)
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

! If this is not the right time to print, return

who = io_control%print_time_who
when = io_control%print_time_when

if (.not. any(this_time == when) .or. who == NO_ONE) return

! If I'm the master and only slaves print, return

if (my_proc(procs) == MASTER .and. who == SLAVES) return

! stop the clocks

call pause_watch(all_watches)

nullify(get_times)

! get the timing information

call inquiry_stopwatch(default_clocks)
nclock = count(default_clocks /= ' ')
if (my_proc(procs) /= MASTER) then
   call read_watch(get_times,(/ &
      ttotal, trefine, trecon, tpartition, tdistribute, tassemble, tsolve, &
      ptotal, prefine, precon, ppartition, pdistribute, passemble, psolve, &
      ctrecon, ctpartition, ctdistribute, ctassemble, ctsolve, &
      cprecon, cppartition, cpdistribute, cpassemble, cpsolve /))
endif

! If the master will print, pass it the info

if (who /= SLAVES) then
   if (my_proc(procs) /= MASTER) then
      allocate(my_times(30,nclock),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in print_time_info",procs=procs)
         return
      endif
      my_times(1:24,:) = get_times
      my_times(25,:) = get_times(2,:) + get_times(3,:) ! tref+trecon
      my_times(26,:) = get_times(4,:) + get_times(5,:) ! tpart+tdist
      my_times(27,:) = get_times(6,:) + get_times(7,:) ! tassem+tsolve
      my_times(28,:) = get_times( 9,:) + get_times(10,:) ! pref+precon
      my_times(29,:) = get_times(11,:) + get_times(12,:) ! ppart+pdist
      my_times(30,:) = get_times(13,:) + get_times(14,:) ! passem+psolve
      call phaml_send(procs,MASTER,no_ints,0, &
                    reshape(my_times,(/30*nclock/)),30*nclock,tag)
      deallocate(my_times,stat=astat)
   else
      allocate(all_times(num_proc(procs),30,nclock),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in print_time_info",procs=procs)
         return
      endif
      do i=1,num_proc(procs)
         call phaml_recv(procs,proc,recv_int,ni,recv_real,nr,tag)
         if (associated(recv_real)) then
            all_times(proc,:,:) = reshape(recv_real,(/30,nclock/))
            deallocate(recv_real,stat=astat)
         else
            call warning("Did not receive times from proc ",intlist=(/proc/))
            all_times(proc,:,:) = 0.0_my_real
         endif
      end do
   endif
endif

! print the info, if requested

if (my_proc(procs) /= MASTER) then
 if (who == SLAVES .or. who == EVERYONE) then
   fmt = '(3a,f8.2)'
   if (which == 'cycle' .or. which == 'both') then
      write(outunit,"(A)")
      write(outunit,"(A)") 'Times for this cycle:'
      do i=1,nclock
         write(outunit,fmt) 'refine,       ',default_clocks(i),': ',get_times(9,i)
         write(outunit,fmt) 'reconcile,    ',default_clocks(i),': ',get_times(10,i)
         write(outunit,fmt) 'partition,    ',default_clocks(i),': ',get_times(11,i)
         write(outunit,fmt) 'distribute,   ',default_clocks(i),': ',get_times(12,i)
         write(outunit,fmt) 'assemble,     ',default_clocks(i),': ',get_times(13,i)
         write(outunit,fmt) 'solve,        ',default_clocks(i),': ',get_times(14,i)
         write(outunit,fmt) 'total,        ',default_clocks(i),': ',get_times(8,i)
         write(outunit,fmt) 'comm recon,   ',default_clocks(i),': ',get_times(20,i)
         write(outunit,fmt) 'comm part,    ',default_clocks(i),': ',get_times(21,i)
         write(outunit,fmt) 'comm distrib, ',default_clocks(i),': ',get_times(22,i)
         write(outunit,fmt) 'comm assemble,',default_clocks(i),': ',get_times(23,i)
         write(outunit,fmt) 'comm solve,   ',default_clocks(i),': ',get_times(24,i)
      end do
   endif
   if (which == 'total' .or. which == 'both') then
      write(outunit,"(A)")
      write(outunit,"(A)") 'Total times so far:'
      do i=1,nclock
         write(outunit,fmt) 'refine,       ',default_clocks(i),': ',get_times(2,i)
         write(outunit,fmt) 'reconcile,    ',default_clocks(i),': ',get_times(3,i)
         write(outunit,fmt) 'partition,    ',default_clocks(i),': ',get_times(4,i)
         write(outunit,fmt) 'distribute,   ',default_clocks(i),': ',get_times(5,i)
         write(outunit,fmt) 'assemble,     ',default_clocks(i),': ',get_times(6,i)
         write(outunit,fmt) 'solve,        ',default_clocks(i),': ',get_times(7,i)
         write(outunit,fmt) 'total,        ',default_clocks(i),': ',get_times(1,i)
         write(outunit,fmt) 'comm recon,   ',default_clocks(i),': ',get_times(15,i)
         write(outunit,fmt) 'comm part,    ',default_clocks(i),': ',get_times(16,i)
         write(outunit,fmt) 'comm distrib, ',default_clocks(i),': ',get_times(17,i)
         write(outunit,fmt) 'comm assemble,',default_clocks(i),': ',get_times(18,i)
         write(outunit,fmt) 'comm solve,   ',default_clocks(i),': ',get_times(19,i)
      end do
   endif
 endif
else
 write(fmt,'(a,i4,a)') '(3a,',num_proc(procs),'f8.2)'
 if (who == MASTER_ALL) then
   if (which == 'cycle' .or. which == 'both') then
      write(outunit,"(A)")
      write(outunit,"(A)") 'Times for this cycle, individual processors:'
      do i=1,nclock
         write(outunit,fmt) 'refine,       ',default_clocks(i),': ',all_times(:,9,i)
         write(outunit,fmt) 'reconcile,    ',default_clocks(i),': ',all_times(:,10,i)
         write(outunit,fmt) 'ref+recon,    ',default_clocks(i),': ',all_times(:,28,i)
         write(outunit,fmt) 'partition,    ',default_clocks(i),': ',all_times(:,11,i)
         write(outunit,fmt) 'distribute,   ',default_clocks(i),': ',all_times(:,12,i)
         write(outunit,fmt) 'part+dist,    ',default_clocks(i),': ',all_times(:,29,i)
         write(outunit,fmt) 'assemble,     ',default_clocks(i),': ',all_times(:,13,i)
         write(outunit,fmt) 'solve,        ',default_clocks(i),': ',all_times(:,14,i)
         write(outunit,fmt) 'assem+solve,  ',default_clocks(i),': ',all_times(:,30,i)
         write(outunit,fmt) 'total,        ',default_clocks(i),': ',all_times(:,8,i)
         write(outunit,fmt) 'comm recon,   ',default_clocks(i),': ',all_times(:,20,i)
         write(outunit,fmt) 'comm part,    ',default_clocks(i),': ',all_times(:,21,i)
         write(outunit,fmt) 'comm distrib, ',default_clocks(i),': ',all_times(:,22,i)
         write(outunit,fmt) 'comm assemble,',default_clocks(i),': ',all_times(:,23,i)
         write(outunit,fmt) 'comm solve,   ',default_clocks(i),': ',all_times(:,24,i)
      end do
   endif
   if (which == 'total' .or. which == 'both') then
      write(outunit,"(A)")
      write(outunit,"(A)") 'Total times so far, individual processors:'
      do i=1,nclock
         write(outunit,fmt) 'refine,       ',default_clocks(i),': ',all_times(:,2,i)
         write(outunit,fmt) 'reconcile,    ',default_clocks(i),': ',all_times(:,3,i)
         write(outunit,fmt) 'ref+recon,    ',default_clocks(i),': ',all_times(:,25,i)
         write(outunit,fmt) 'partition,    ',default_clocks(i),': ',all_times(:,4,i)
         write(outunit,fmt) 'distribute,   ',default_clocks(i),': ',all_times(:,5,i)
         write(outunit,fmt) 'part+dist,    ',default_clocks(i),': ',all_times(:,26,i)
         write(outunit,fmt) 'assemble,     ',default_clocks(i),': ',all_times(:,6,i)
         write(outunit,fmt) 'solve,        ',default_clocks(i),': ',all_times(:,7,i)
         write(outunit,fmt) 'assem+solve,  ',default_clocks(i),': ',all_times(:,27,i)
         write(outunit,fmt) 'total,        ',default_clocks(i),': ',all_times(:,1,i)
         write(outunit,fmt) 'comm recon,   ',default_clocks(i),': ',all_times(:,15,i)
         write(outunit,fmt) 'comm part,    ',default_clocks(i),': ',all_times(:,16,i)
         write(outunit,fmt) 'comm distrib, ',default_clocks(i),': ',all_times(:,17,i)
         write(outunit,fmt) 'comm assemble,',default_clocks(i),': ',all_times(:,18,i)
         write(outunit,fmt) 'comm solve,   ',default_clocks(i),': ',all_times(:,19,i)
      end do
   endif
 endif
 if (who == MASTER .or. who == EVERYONE .or. who == MASTER_ALL) then
   if (which == 'cycle' .or. which == 'both') then
      write(outunit,"(A)")
      write(outunit,"(A)") 'Times for this cycle, max over all processors:'
      do i=1,nclock
         write(outunit,fmt) 'refine,       ',default_clocks(i),': ',maxval(all_times(:,9,i))
         write(outunit,fmt) 'reconcile,    ',default_clocks(i),': ',minval(all_times(:,10,i))
         write(outunit,fmt) 'ref+recon,    ',default_clocks(i),': ',maxval(all_times(:,28,i))
         write(outunit,fmt) 'partition,    ',default_clocks(i),': ',maxval(all_times(:,11,i))
         write(outunit,fmt) 'distribute,   ',default_clocks(i),': ',maxval(all_times(:,12,i))
         write(outunit,fmt) 'part+dist,    ',default_clocks(i),': ',maxval(all_times(:,29,i))
         write(outunit,fmt) 'assemble,     ',default_clocks(i),': ',maxval(all_times(:,13,i))
         write(outunit,fmt) 'solve,        ',default_clocks(i),': ',maxval(all_times(:,14,i))
         write(outunit,fmt) 'assem+solve,  ',default_clocks(i),': ',maxval(all_times(:,30,i))
         write(outunit,fmt) 'total,        ',default_clocks(i),': ',maxval(all_times(:,8,i))
         write(outunit,fmt) 'comm recon,   ',default_clocks(i),': ',minval(all_times(:,20,i))
         write(outunit,fmt) 'comm part,    ',default_clocks(i),': ',maxval(all_times(:,21,i))
         write(outunit,fmt) 'comm distrib, ',default_clocks(i),': ',maxval(all_times(:,22,i))
         write(outunit,fmt) 'comm assemble,',default_clocks(i),': ',maxval(all_times(:,23,i))
         write(outunit,fmt) 'comm solve,   ',default_clocks(i),': ',maxval(all_times(:,24,i))
      end do
   endif
   if (which == 'total' .or. which == 'both') then
      write(outunit,"(A)")
      write(outunit,"(A)") 'Total times so far, max over all processors:'
      do i=1,nclock
         write(outunit,fmt) 'refine,       ',default_clocks(i),': ',maxval(all_times(:,2,i))
         write(outunit,fmt) 'reconcile,    ',default_clocks(i),': ',minval(all_times(:,3,i))
         write(outunit,fmt) 'ref+recon,    ',default_clocks(i),': ',maxval(all_times(:,25,i))
         write(outunit,fmt) 'partition,    ',default_clocks(i),': ',maxval(all_times(:,4,i))
         write(outunit,fmt) 'distribute,   ',default_clocks(i),': ',maxval(all_times(:,5,i))
         write(outunit,fmt) 'part+dist,    ',default_clocks(i),': ',maxval(all_times(:,26,i))
         write(outunit,fmt) 'assemble,     ',default_clocks(i),': ',maxval(all_times(:,6,i))
         write(outunit,fmt) 'solve,        ',default_clocks(i),': ',maxval(all_times(:,7,i))
         write(outunit,fmt) 'assem+solve,  ',default_clocks(i),': ',maxval(all_times(:,27,i))
         write(outunit,fmt) 'total,        ',default_clocks(i),': ',maxval(all_times(:,1,i))
         write(outunit,fmt) 'comm recon,   ',default_clocks(i),': ',minval(all_times(:,15,i))
         write(outunit,fmt) 'comm part,    ',default_clocks(i),': ',maxval(all_times(:,16,i))
         write(outunit,fmt) 'comm distrib, ',default_clocks(i),': ',maxval(all_times(:,17,i))
         write(outunit,fmt) 'comm assemble,',default_clocks(i),': ',maxval(all_times(:,18,i))
         write(outunit,fmt) 'comm solve,   ',default_clocks(i),': ',maxval(all_times(:,19,i))
      end do
   endif
 endif
endif

if (associated(get_times)) deallocate(get_times,stat=astat)
if (allocated(all_times)) deallocate(all_times,stat=astat)

call end_pause_watch(all_watches)

return
end subroutine print_time_info

!          ----
subroutine init(phaml_solution,procs,io_control,solver_control, &
                refine_control,partition_method,task, &
                degree)
!          ----

!----------------------------------------------------
! This routine performs the initializations.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type (phaml_solution_type), intent(inout), target :: phaml_solution
type (proc_info), intent(inout) :: procs
type (io_options), intent(in) :: io_control
type(refine_options), intent(in) :: refine_control
type (solver_options), intent(inout) :: solver_control
integer, intent(in) :: partition_method, task, degree

!----------------------------------------------------
! Local variables
 
integer :: ni, astat
real (my_real), pointer :: temp(:)
type(grid_type), pointer :: grid

!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

! create the timers

call create_watch(ttotal,name='total')
call create_watch(trefine,name='refine')
call create_watch(trecon,name='reconcile')
call create_watch(tpartition,name='partition')
call create_watch(tdistribute,name='distribute')
call create_watch(tassemble,name='assemble')
call create_watch(tsolve,name='solve')
call create_watch(ptotal,name='total')
call create_watch(prefine,name='refine')
call create_watch(precon,name='reconcile')
call create_watch(ppartition,name='partition')
call create_watch(pdistribute,name='distribute')
call create_watch(passemble,name='assemble')
call create_watch(psolve,name='solve')
call create_watch(ctrecon,name='comm recon')
call create_watch(ctpartition,name='comm partition')
call create_watch(ctdistribute,name='comm distribute')
call create_watch(ctassemble,name='comm assemble')
call create_watch(ctsolve,name='comm solve')
call create_watch(cprecon,name='comm recon')
call create_watch(cppartition,name='comm partition')
call create_watch(cpdistribute,name='comm distribute')
call create_watch(cpassemble,name='comm assemble')
call create_watch(cpsolve,name='comm solve')
call create_watchgroup((/ &
   ttotal, trefine, trecon, tpartition, tdistribute, tassemble, tsolve, &
   ptotal, prefine, precon, ppartition, pdistribute, passemble, psolve, &
   ctrecon, ctpartition, ctdistribute, ctassemble, ctsolve, &
   cprecon, cppartition, cpdistribute, cpassemble, cpsolve /),all_watches)

! wait for all processors (except master ) to get to here before starting timing

if (my_proc(procs) /= MASTER) then
   call phaml_barrier(procs)
endif

call start_watch((/ptotal,ttotal/))

grid => phaml_solution%grid

! put eq_type and system size into solver control

solver_control%eq_type = phaml_solution%eq_type
solver_control%system_size = phaml_solution%system_size

! if using the equilibrated residual error estimator, set the edge mass matrix

if (refine_control%error_estimator == EQUILIBRATED_RESIDUAL .or. &
    (refine_control%reftype == HP_ADAPTIVE .and. &
     refine_control%hp_strategy == HP_NLP) .or. &
    (refine_control%reftype == HP_ADAPTIVE .and. &
     refine_control%hp_strategy == HP_NEXT3P)) then
   if (my_proc(procs) /= MASTER) then
      call set_edge_mass(refine_control%max_deg)
   endif
endif

! initialize Zoltan

call zoltan_init(phaml_solution%lb,partition_method, &
                 phaml_solution%zoltan_param_file,phaml_solution%procs)

! make sure the allocation for solutions is right

call realloc_solution(grid,degree,phaml_solution%system_size, &
                      solver_control%num_eval)

! make sure the allocation for eigenvalues, etc., is right

if (phaml_solution%eq_type == EIGENVALUE) then
   if (.not. associated(grid%eigenvalue)) then
      allocate(grid%eigenvalue(solver_control%num_eval), &
               grid%eigenprob_l2_resid(solver_control%num_eval), &
               grid%eigenprob_variance(solver_control%num_eval),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("memory allocation failed in init",procs=procs)
         return
      endif
      grid%eigenvalue = 0.0_my_real
   elseif (size(grid%eigenvalue) /= solver_control%num_eval) then
      deallocate(grid%eigenprob_l2_resid,grid%eigenprob_variance, &
                 stat=astat)
      allocate(temp(solver_control%num_eval), &
               grid%eigenprob_l2_resid(solver_control%num_eval), &
               grid%eigenprob_variance(solver_control%num_eval),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("memory allocation failed in init",procs=procs)
         return
      endif
      ni=min(size(temp),size(grid%eigenvalue))
      temp = 0.0_my_real
      temp(1:ni) = grid%eigenvalue
      deallocate(grid%eigenvalue,stat=astat)
      grid%eigenvalue => temp
   endif
   grid%num_eval = solver_control%num_eval
else
   grid%num_eval = 0
endif

! Dirichlet boundary conditions and the exact solution may have changed,
! for example in a time dependent problem, so reset them

call reset_dirich_exact(grid)

! print and draw initial grid

call print_grid_info(grid,procs,io_control,phaml_solution%still_sequential, &
                     (/FREQUENTLY,PHASES/),310)
grid_changed=.true.
call draw_grid(grid,procs,io_control,refine_control, &
               phaml_solution%i_draw_grid,phaml_solution%master_draws_grid, &
               phaml_solution%still_sequential,(/FREQUENTLY,PHASES/), &
               partition_method,phaml_solution%lb)

end subroutine init

!          --------------------
subroutine check_end_sequential(phaml_solution,grid,procs,refine_control, &
                                io_control,sequential_vert,prebalance, &
                                postbalance,partition_method,loop, &
                                loop_end_sequential,still_sequential)
!          --------------------

!----------------------------------------------------
! This routine checks to see if the grid is big enough to change from solving
! sequentially to solving in parallel, and distributes the grid if it is
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(phaml_solution_type), intent(in) :: phaml_solution
type(grid_type), intent(inout) :: grid
type(proc_info), intent(in) :: procs
type(refine_options), intent(in) :: refine_control
type (io_options), intent(in) :: io_control
integer, intent(in) :: sequential_vert
integer, intent(in) :: prebalance, postbalance
integer, intent(in) :: partition_method
integer, intent(in) :: loop
integer, intent(inout) :: loop_end_sequential
logical, intent(inout) :: still_sequential
!----------------------------------------------------
! Local variables:

integer :: nvert
type(hash_key), pointer :: export_gid(:)
integer, pointer :: export_proc(:)
!----------------------------------------------------
! Begin executable code

! check size of grid

call get_grid_info(grid,procs,still_sequential,130+loop,total_nvert=nvert)
if (nvert > sequential_vert) then

! big enough; partition and distribute the grid

   still_sequential = .false.
   loop_end_sequential = loop
   if (prebalance /= BALANCE_NONE) then
      call partition(grid,procs,phaml_solution%lb,refine_control,.true., &
                     partition_method,prebalance,num_proc(procs),export_gid, &
                     export_proc,first_call=.true.)
   else
      call partition(grid,procs,phaml_solution%lb,refine_control,.false., &
                     partition_method,postbalance,num_proc(procs),export_gid, &
                     export_proc,first_call=.true.)
   endif
   call redistribute(grid,procs,refine_control,export_gid,export_proc, &
                     first_call=.true.)
   call reconcile(grid,procs,refine_control,still_sequential)
   if (associated(export_gid)) then
      deallocate(export_gid,export_proc)
   endif
   grid_changed = .true.
   call draw_grid(grid,procs,io_control,refine_control, &
                  phaml_solution%i_draw_grid,phaml_solution%master_draws_grid, &
                  still_sequential,(/FREQUENTLY/),partition_method, &
                  phaml_solution%lb)
endif

end subroutine check_end_sequential

!          ------
subroutine balance(phaml_solution, grid, procs, io_control, refine_control, &
                   partition_method, balance_what, still_sequential, loop, &
                   loop_end_sequential, prebalancecall)
!          ------

!----------------------------------------------------
! This routine performs the steps for load balancing
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(phaml_solution_type), intent(in) :: phaml_solution
type(grid_type), intent(inout) :: grid
type(proc_info), intent(in) :: procs
type (io_options), intent(in) :: io_control
type(refine_options), intent(in) :: refine_control
integer, intent(in) :: partition_method, balance_what, loop, loop_end_sequential
logical, intent(in) :: still_sequential, prebalancecall
!----------------------------------------------------
! Local variables:

logical :: repart, redist
integer :: nleaf, minleaf, maxleaf, numexp, astat
type(hash_key), pointer :: export_gid(:)
integer, pointer :: export_proc(:)
!----------------------------------------------------
! Begin executable code

! see if repartitioning is really necessary

! some cases where partitioning is not needed:
!   master does not participate
!   number of processors is 1
!   program compiled for sequential execution
!   haven't changed to parallel execution yet
!   first time in parallel execution (prebalancing only)
if (my_proc(procs) == MASTER .or. num_proc(procs) == 1 .or. &
    PARALLEL==SEQUENTIAL .or. still_sequential .or. &
    (prebalancecall .and. loop == loop_end_sequential)) then
   repart = .false.

! some cases where partitioning must be done
!   predictive load balancing, except the first time
! TEMP why am I required to balance if it is done before refinement?
elseif (prebalancecall) then
   repart = .true.

! otherwise, partition if the load is sufficiently out of balance
else
   call get_grid_info(grid,procs,still_sequential,160+loop, &
                      nelem_leaf_own=nleaf,no_master=.true.)
   minleaf = phaml_global_min(procs,nleaf,170+loop)
   maxleaf = phaml_global_max(procs,nleaf,180+loop)
   repart = ( (maxleaf-minleaf)/float(maxleaf) > 0.05)
endif

! partition

if (repart) then
   call partition(grid,procs,phaml_solution%lb,refine_control, &
                  prebalancecall, partition_method, balance_what, &
                  num_proc(procs),export_gid,export_proc) 

! redistribute

! see if enough elements are moved to be worthwhile
   if (associated(export_gid)) then
      numexp = size(export_gid)
   else
      numexp = 0
   endif
   numexp = phaml_global_sum(procs,numexp,190+loop)
   call get_grid_info(grid,procs,still_sequential,200+loop, &
                      total_nelem_leaf=nleaf,no_master=.true.)
   redist = (numexp/float(nleaf) > .05)
   if (redist) then
      call redistribute(grid,procs,refine_control,export_gid,export_proc)
      call reconcile(grid,procs,refine_control,still_sequential)
   endif
   if (associated(export_gid)) then
      deallocate(export_gid,export_proc,stat=astat)
   endif

endif

call print_grid_info(grid,procs,io_control,still_sequential, &
                     (/FREQUENTLY/),210+loop)
grid_changed = .true. ! because MASTER doesn't know if redistributed
call draw_grid(grid,procs,io_control,refine_control, &
               phaml_solution%i_draw_grid, &
               phaml_solution%master_draws_grid,still_sequential, &
               (/FREQUENTLY/),partition_method,phaml_solution%lb)

end subroutine balance

!        -------------------
function check_stall_special(grid,phase)
!        -------------------

!----------------------------------------------------
! This routine performs a test for refinement stalling in the special case
! of a degree 1 element with three degree 2 sides being p refined to degree 2,
! which does not change the number of vertices, elements or equations.
! If phase=1, it flags degree 1 elements.  If phase=2, it sees if any
! of the flagged elements were p refined.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: phase
logical :: check_stall_special
!----------------------------------------------------
! Local variables:

logical(small_logical), save, allocatable :: marked(:)
integer :: lev, elem
!----------------------------------------------------
! Begin executable code

select case (phase)
case (1)
   if (allocated(marked)) then
      if (size(marked) /= size(grid%element)) then
         deallocate(marked)
      endif
   endif
   if (.not. allocated(marked)) then
      allocate(marked(size(grid%element)))
   endif
   marked = .false.
   do lev=1,grid%nlev
      elem = grid%head_level_elem(lev)
      do while (elem /= END_OF_LIST)
         if (grid%element(elem)%degree == 1) then
            marked(elem) = .true.
         endif
         elem = grid%element(elem)%next
      end do
   end do
   check_stall_special = .false.

case(2)
   check_stall_special = .true.
outer: do lev=1,grid%nlev
      elem = grid%head_level_elem(lev)
      do while (elem /= END_OF_LIST)
         if (marked(elem) .and. grid%element(elem)%degree /= 1) then
            check_stall_special = .false.
            exit outer
         endif
         elem = grid%element(elem)%next
      end do
   end do outer

end select

end function check_stall_special

!          ------------
subroutine phaml_create(phaml_solution,nproc,draw_grid_who, &
                        spawn_form,debug_command,display,graphics_host, &
                        output_unit,error_unit,output_now,id,system_size, &
                        eq_type,max_blen,triangle_files,update_umod)
!          ------------

!----------------------------------------------------
! This routine creates a variable to contain a phaml_solution
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type (phaml_solution_type) :: phaml_solution
integer, optional :: nproc,draw_grid_who,spawn_form,output_unit,error_unit, &
                     output_now,id,system_size,eq_type
character(len=*), optional, intent(in) :: debug_command,display
character(len=*), optional, intent(in) :: graphics_host
real(my_real), optional, intent(in) :: max_blen
character(len=*), optional, intent(in) :: triangle_files
logical, optional, intent(in) :: update_umod
!----------------------------------------------------

!----------------------------------------------------
! Local variables
integer :: local_nproc, loc_spawn_form, loc_pde_id, astat, proc, ni, nr, &
           nalloc
integer, allocatable :: imess(:)
real(my_real) :: rmess(4)
integer, pointer :: recv_int(:)
real(my_real), pointer :: recv_real(:)
integer, parameter :: DEBUGLEN = 64
character(len=DEBUGLEN) :: loc_debug_command, loc_display
logical :: loc_update_umod

!----------------------------------------------------
! Begin executable code

! set I/O units

if (present(output_now)) then
   if (output_now < 0) then
      write(6,"(A)")
      write(6,"(A)") "------------------------------------------------------"
      write(6,"(3A)") "          PHAML Version ",version_number," ERROR"
      write(6,"(A,I11)") "phaml_create: output_now must be a nonnegative integer, it is ",output_now
      write(6,"(A)") "------------------------------------------------------"
      write(6,"(A)")
      stop
   endif
   outunit = output_now
   errunit = output_now
else
   outunit = 6
   errunit = 6
endif

if (present(output_unit)) then
   if (output_unit < 0) then
      call fatal("phaml_create: output_unit must be a nonnegative integer.  It is",intlist=(/output_unit/))
      stop
   endif
   phaml_solution%outunit = output_unit
else
   phaml_solution%outunit = 6
endif

if (present(error_unit)) then
   if (error_unit < 0) then
      call fatal("phaml_create: error_unit must be a nonnegative integer.  It is",intlist=(/error_unit/))
      stop
   endif
   phaml_solution%errunit = error_unit
else
   phaml_solution%errunit = 0
endif

! check the size of nproc

if (present(nproc)) then
   if (nproc <= 0) then
      call warning("nproc must be a positive integer.  Resetting to 1")
      local_nproc = 1
   else
      local_nproc = nproc
   endif
else
   local_nproc = 1
endif

! set the graphics options

if (present(draw_grid_who)) then
   select case (draw_grid_who)
   case (MASTER)
      if (PARALLEL==SEQUENTIAL) then
         phaml_solution%i_draw_grid = .true.
         phaml_solution%master_draws_grid = .false.
      else
         phaml_solution%i_draw_grid = .false.
         phaml_solution%master_draws_grid = .true.
      endif
   case (SLAVES)
      phaml_solution%i_draw_grid = .true.
      phaml_solution%master_draws_grid = .false.
   case (EVERYONE)
      phaml_solution%i_draw_grid = .true.
      phaml_solution%master_draws_grid = .true.
   case (NO_ONE)
      phaml_solution%i_draw_grid = .false.
      phaml_solution%master_draws_grid = .false.
   case default
      call warning("illegal value for draw_grid_who.  Setting to NO_ONE.")
      phaml_solution%i_draw_grid = .false.
      phaml_solution%master_draws_grid = .false.
   end select
else
   phaml_solution%i_draw_grid = .false.
   phaml_solution%master_draws_grid = .false.
endif

if (present(spawn_form)) then
   loc_spawn_form = spawn_form
else
   loc_spawn_form = NORMAL_SPAWN
endif

if (present(debug_command)) then
   loc_debug_command = debug_command
else
   loc_debug_command = "gdb"
endif

if (present(display)) then
   loc_display = display
else
   loc_display = "default"
endif

if (present(graphics_host)) then
   phaml_solution%graphics_host = trim(graphics_host)
else
   phaml_solution%graphics_host = "anywhere"
endif

! set initial grid specifier

if (present(max_blen)) then
   if (max_blen <= 0.0_my_real) then
      call fatal("phaml_create: max_blen must be positive.  It is ",reallist=(/max_blen/))
      stop
   endif
   phaml_solution%grid%max_blen = max_blen
else
   phaml_solution%grid%max_blen = huge(0.0_my_real)
endif

if (present(triangle_files)) then
   if (len(triangle_files) > FN_LEN) then
      call warning("Name of triangle files is too long and will be truncated.",&
                   "Maximum length (FN_LEN in global.f90) is", &
                   intlist=(/FN_LEN/))
   endif
   phaml_solution%grid%triangle_files = triangle_files
else
   phaml_solution%grid%triangle_files = "domain"
endif

if (present(update_umod)) then
   loc_update_umod = update_umod
else
   loc_update_umod = .false.
endif

if (present(id)) then
   loc_pde_id = id
else
   loc_pde_id = 0
endif

if (present(system_size)) then
   if (system_size <= 0) then
      call fatal("phaml_create: system_size must be positive.  It is ",intlist=(/system_size/))
      stop
   endif
   phaml_solution%system_size = system_size
else
   phaml_solution%system_size = 1
endif

if (present(eq_type)) then
   phaml_solution%eq_type = eq_type
else
   phaml_solution%eq_type = ELLIPTIC
endif

! initialize the message passing communications package

call init_comm(phaml_solution%procs,loc_spawn_form,phaml_solution%i_draw_grid, &
               phaml_solution%master_draws_grid, &
               phaml_solution%graphics_host,phaml_solution%outunit, &
               phaml_solution%errunit,phaml_solution%system_size, &
               phaml_solution%eq_type,loc_pde_id,local_nproc, &
               phaml_solution%grid%max_blen,phaml_solution%grid%triangle_files,&
               loc_update_umod,loc_debug_command,loc_display)

! update the usermod variables

if (loc_update_umod) then
   if (PARALLEL /= SEQUENTIAL) then
      if (my_proc(phaml_solution%procs) /= MASTER) then
         call phaml_recv(phaml_solution%procs,proc,recv_int,ni,recv_real,nr,101)
         if (recv_int(1) /= 9) then
            ierr = PHAML_INTERNAL_ERROR
            call fatal("In phaml_create, slave received code that is not update_usermod.")
            stop
         endif
      endif
      call update_usermod(phaml_solution)
   endif
endif

! create zoltan load balancing object

call zoltan_create_lb(phaml_solution%lb,phaml_solution%procs)

phaml_solution%still_sequential = .true.
phaml_solution%pde_id = loc_pde_id
my_pde_id = loc_pde_id

! send size(pde) to the slaves

if (PARALLEL /= SEQUENTIAL) then
   if (my_proc(phaml_solution%procs) == MASTER) then
      allocate(imess(1))
      if (allocated(pde)) then
         imess(1) = size(pde)
      else
         imess(1) = 0
      endif
      do proc=1,num_proc(phaml_solution%procs)
         call phaml_send(phaml_solution%procs,proc,imess,1,(/0._my_real/),0,301)
      end do
      deallocate(imess)
   else
      call phaml_recv(phaml_solution%procs,proc,recv_int,ni,recv_real,nr,301)
      if (recv_int(1) /= 0) then
         allocate(pde(recv_int(1)),stat=astat)
         if (astat /= 0) then
            call warning("allocation of pde failed.")
         endif
      endif
      if (associated(recv_int)) deallocate(recv_int,stat=astat)
   endif
endif

! initialize the grid

call set_grid_for_old_soln(phaml_solution%grid)

phaml_solution%grid%system_size = phaml_solution%system_size
if (phaml_solution%eq_type == EIGENVALUE) then
   phaml_solution%grid%num_eval = 1
else
   phaml_solution%grid%num_eval = 0
endif

if (my_proc(phaml_solution%procs) == MASTER) then
   nalloc = 1
else
   nalloc = 25000
endif
call allocate_grid(phaml_solution%grid,nalloc,phaml_solution%grid%num_eval, &
                   phaml_solution%eq_type,1)

call init_grid(phaml_solution%grid,phaml_solution%procs,1, &
               my_proc(phaml_solution%procs),.true.)

if (my_proc(phaml_solution%procs) /= MASTER) then
   if (my_proc(phaml_solution%procs) == 1) then
      rmess(1) = phaml_solution%grid%boundbox_min%x
      rmess(2) = phaml_solution%grid%boundbox_max%x
      rmess(3) = phaml_solution%grid%boundbox_min%y
      rmess(4) = phaml_solution%grid%boundbox_max%y
      call phaml_send(phaml_solution%procs,MASTER,(/0/),0,rmess,4,300)
   endif
else
   call phaml_recv(phaml_solution%procs,proc,recv_int,ni,recv_real,nr,300)
   phaml_solution%grid%boundbox_min%x = recv_real(1)
   phaml_solution%grid%boundbox_max%x = recv_real(2)
   phaml_solution%grid%boundbox_min%y = recv_real(3)
   phaml_solution%grid%boundbox_max%y = recv_real(4)
   if (associated(recv_real)) deallocate(recv_real,stat=astat)
endif

! initialize graphics
 
grid_changed = .true.
if ((my_proc(phaml_solution%procs)==MASTER .and. phaml_solution%master_draws_grid) &
    .or. &
    (my_proc(phaml_solution%procs)/=MASTER .and. phaml_solution%i_draw_grid)) then
   allocate(imess(4))
   imess(1) = GRAPHICS_INIT
   imess(2) = 100000
   imess(3) = 100000
   imess(4) = 25000
   rmess(1) = phaml_solution%grid%boundbox_min%x
   rmess(2) = phaml_solution%grid%boundbox_max%x
   rmess(3) = phaml_solution%grid%boundbox_min%y
   rmess(4) = phaml_solution%grid%boundbox_max%y
   if (PARALLEL == SEQUENTIAL) then
      call sequential_send(imess,4,rmess,4)
   else
      call phaml_send(phaml_solution%procs,graphics_proc(phaml_solution%procs),&
                      imess,4,rmess,4,101)
   endif
   deallocate(imess)
endif

end subroutine phaml_create

!          -------------
subroutine phaml_destroy(phaml_solution,finalize_mpi)
!          -------------

!----------------------------------------------------
! This routine frees up space by destroying phaml_solution
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type (phaml_solution_type) :: phaml_solution
logical, intent(in), optional :: finalize_mpi
!----------------------------------------------------

!----------------------------------------------------
! Local variables

integer :: send_int(1),proc,my_processor,astat
real (my_real) :: no_reals(1)
logical :: loc_finalize

!----------------------------------------------------
! Begin executable code

outunit = phaml_solution%outunit
errunit = phaml_solution%errunit

my_processor = my_proc(phaml_solution%procs)

if (my_processor == MASTER) then
   if (present(finalize_mpi)) then
      loc_finalize = finalize_mpi
   else
      loc_finalize = .true.
   endif
else
   loc_finalize = .true.
endif

! notify the slaves to destroy phaml_solution

if (my_processor == MASTER) then
   send_int(1) = 2
   do proc=1,num_proc(phaml_solution%procs)
      call phaml_send(phaml_solution%procs,proc,send_int,1,no_reals,0,101)
   end do
end if

! destroy zoltan load balancing object

call zoltan_destroy_lb(phaml_solution%lb,phaml_solution%procs)

! destroy phaml_solution

call deallocate_grid(phaml_solution%grid)
if (my_processor /= MASTER) then
   if (allocated(pde)) deallocate(pde,stat=astat)
endif
call terminate_comm(phaml_solution%procs,loc_finalize)

end subroutine phaml_destroy

!        ---------------
function phaml_integrate(phaml_solution,kernel,comp1,eigen1,comp2,eigen2,p,q)
!        ---------------

!----------------------------------------------------
! This function returns an integral of a computed solution or product of two
! computed solutions, or powers of computed solution, weighted by the
! function kernel, i.e. it computes
! 
!        /\
!       |             p          q
!       \  kernel * u        * u
!        |           choice1    choice2
!      \/
! 
! kernel is an integer that allows you to select among different kernel
! functions.  It is passed to the user written phaml_integral_kernel where
! it can be used, for example, in a case statement to determine the kernel
! to use.
! 
! comp1, comp2, eigen1 and eigen2 determine which component(s) and which
! eigenfunction(s) to use.  If comp1 (eigen1) is omitted, then
! comp2 (eigen2) must also be omitted and the first component of the first
! eigenfunction is used.  If comp2 and eigen2 are both omitted, then u_choice2
! is omitted from the integral.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type (phaml_solution_type), intent(in) :: phaml_solution
integer, intent(in) :: kernel
integer, optional, intent(in) :: comp1, eigen1, comp2, eigen2, p, q
real(my_real) :: phaml_integrate

!----------------------------------------------------
! Local variables

integer :: loc_comp1, loc_comp2, loc_eigen1, loc_eigen2, my_processor, proc, &
           ni, nr, astat, elem, loc_p, loc_q, i
real(my_real) :: loc_integrate
integer :: send_int(8)
integer, pointer :: recv_int(:)
real(my_real), pointer :: recv_real(:)
!----------------------------------------------------
! Begin executable code

my_processor = my_proc(phaml_solution%procs)

! Set local copies of arguments

if (present(comp1)) then
   loc_comp1 = comp1
else
   loc_comp1 = 1
endif
if (present(eigen1)) then
   loc_eigen1 = eigen1
else
   loc_eigen1 = 1
endif
if (.not. present(comp2) .and. .not. present(eigen2)) then
   loc_comp2 = 0
   loc_eigen2 = 0
else
   if (present(comp2)) then
      loc_comp2 = comp2
   else
      loc_comp2 = 1
   endif
   if (present(eigen2)) then
      loc_eigen2 = eigen2
   else
      loc_eigen2 = 1
   endif
endif
if (present(p)) then
   loc_p = p
else
   loc_p = 1
endif
if (present(q)) then
   loc_q = q
else
   loc_q = 1
endif

! Check validity of arguments

if (loc_comp1 < 1 .or. loc_comp1 > phaml_solution%system_size .or. &
    loc_eigen1 < 1 .or. loc_eigen1 > max(1,phaml_solution%grid%num_eval)) then
   call warning("phaml_integrate: comp1 or eigen1 is not a valid solution number", &
                "Returning 1.0")
   phaml_integrate = 1.0_my_real
   return
endif
if (.not. present(comp1) .and. present(comp2)) then
   call warning("phaml_integrate: comp1 is not present and comp2 is.", &
                "Returning 1.0")
   phaml_integrate = 1.0_my_real
   return
endif
if (.not. present(eigen1) .and. present(eigen2)) then
   call warning("phaml_integrate: eigen1 is not present and eigen2 is.", &
                "Returning 1.0")
   phaml_integrate = 1.0_my_real
   return
endif
if (my_processor == MASTER .and. present(comp2)) then
   if (comp2 == 0) then
      call warning("phaml_integrate: comp2 is not a valid solution number", &
                   "Returning 1.0")
      phaml_integrate = 1.0_my_real
      return
   endif
endif
if (my_processor == MASTER .and. present(eigen2)) then
   if (eigen2 == 0) then
      call warning("phaml_integrate: eigen2 is not a valid solution number", &
                   "Returning 1.0")
      phaml_integrate = 1.0_my_real
      return
   endif
endif
if (loc_comp2 < 0 .or. loc_comp2 > phaml_solution%system_size .or. &
    loc_eigen2 < 0 .or. loc_eigen2 > max(1,phaml_solution%grid%num_eval)) then
   call warning("phaml_integrate: comp2 or eigen2 is not a valid solution number", &
                "Returning 1.0")
   phaml_integrate = 1.0_my_real
   return
endif

! Master tells the slaves to compute the integral over their subdomains
! and sums them

if (my_processor == MASTER) then

   send_int(1) = 12 ! code for integrate
   send_int(2) = kernel
   send_int(3) = loc_comp1
   send_int(4) = loc_eigen1
   send_int(5) = loc_comp2
   send_int(6) = loc_eigen2
   send_int(7) = loc_p
   send_int(8) = loc_q

   do proc=1,num_proc(phaml_solution%procs)
      call phaml_send(phaml_solution%procs,proc,send_int,8,(/0.0_my_real/), &
                      0,101)
   end do

   phaml_integrate = 0.0_my_real

   do i=1,num_proc(phaml_solution%procs)
      call phaml_recv(phaml_solution%procs,proc,recv_int,ni,recv_real,nr,321)
      phaml_integrate = phaml_integrate + recv_real(1)
      deallocate(recv_real,stat=astat)
   end do

! Slave goes through elements to compute integral

else

   loc_integrate = 0.0_my_real

   elem = phaml_solution%grid%head_level_elem(1)
   do while (elem /= END_OF_LIST)
      loc_integrate = loc_integrate + &
                      recur_integrate(elem,phaml_solution%grid,kernel, &
                      loc_comp1,loc_eigen1,loc_comp2,loc_eigen2,loc_p,loc_q)
      elem = phaml_solution%grid%element(elem)%next
   end do

   call phaml_send(phaml_solution%procs,MASTER,(/0/),0,(/loc_integrate/), &
                      1,321)
   phaml_integrate = loc_integrate

endif

end function phaml_integrate

!                  ---------------
recursive function recur_integrate(elem,grid,kernel,comp1,eigen1,comp2,eigen2, &
                                   p,q) result(integ)
!                  ---------------

!----------------------------------------------------
! This routine computes the integral by traversing the tree
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(in) :: grid
integer, intent(in) :: elem, kernel, comp1, eigen1, comp2, eigen2, p, q
real(my_real) :: integ
!----------------------------------------------------
! Local variables:

integer :: children(MAX_CHILD), i, allc(MAX_CHILD), nqpoints, jerr, astat, order
real(my_real) :: xvert(3), yvert(3)
real(my_real), pointer :: qweights(:),xquad(:),yquad(:)
real(my_real), allocatable :: u(:,:,:), kern(:)

interface
   function phaml_integral_kernel(kernel,x,y)
   use global
   integer, intent(in) :: kernel
   real(my_real), intent(in) :: x,y
   real(my_real) :: phaml_integral_kernel
   end function phaml_integral_kernel
end interface

!----------------------------------------------------
! Begin executable code

integ = 0.0_my_real

allc = ALL_CHILDREN
children = get_child_lid(grid%element(elem)%gid,allc,grid%elem_hash)

! if elem is not a leaf, then sum the integrals over each of the children

if (children(1) /= NO_CHILD) then
   do i=1,MAX_CHILD
      integ = integ + recur_integrate(children(i),grid,kernel,comp1,eigen1, &
                                      comp2,eigen2,p,q)
   end do

! otherwise, compute the integral over this element, if I own it

else

   if (grid%element(elem)%iown) then

! RESTRICTION triangles

! get the vertex coordinates

      xvert(1) = grid%vertex(grid%element(elem)%vertex(1))%coord%x
      xvert(2) = grid%vertex(grid%element(elem)%vertex(2))%coord%x
      xvert(3) = grid%vertex(grid%element(elem)%vertex(3))%coord%x
      yvert(1) = grid%vertex(grid%element(elem)%vertex(1))%coord%y
      yvert(2) = grid%vertex(grid%element(elem)%vertex(2))%coord%y
      yvert(3) = grid%vertex(grid%element(elem)%vertex(3))%coord%y

! get the quadrature rule

      if (comp2==0) then
         order = min(MAX_QUAD_ORDER_TRI,grid%element(elem)%degree*p)
      else
         order = min(MAX_QUAD_ORDER_TRI,grid%element(elem)%degree*(p+q))
      endif
      call quadrature_rule_tri(order,xvert,yvert,nqpoints,qweights, &
                               xquad,yquad,jerr,stay_in=.true.)
      if (jerr /= 0) then
         call fatal("Failed to get quadrature rule in phaml_integrate.", &
                    "Integral not computed")
         integ = 0.0_my_real
         return
      endif

      allocate(u(2,2,nqpoints),kern(nqpoints),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("allocation failed in phaml_integrate")
         stop
      endif

! evaluate solution(s)

      call evaluate_soln_local(grid,xquad,yquad,elem,(/comp1/),(/eigen1/), &
                               u(1:1,1:1,:))
      if (comp2/=0 .and. eigen2==0) then
         call evaluate_soln_local(grid,xquad,yquad,elem,(/comp2/),(/1/), &
                                  u(2:2,2:2,:))
      elseif (comp2==0 .and. eigen2/=0) then
         call evaluate_soln_local(grid,xquad,yquad,elem,(/1/),(/eigen2/), &
                                  u(2:2,2:2,:))
      elseif (comp2/=0 .and. eigen2/=0) then
         call evaluate_soln_local(grid,xquad,yquad,elem,(/comp2/),(/eigen2/), &
                                  u(2:2,2:2,:))
      else
         u(2,2,:) = 1.0_my_real
      endif

! evaluate the kernel

      do i=1,nqpoints
         kern(i) = phaml_integral_kernel(kernel,xquad(i),yquad(i))
      end do

! compute the integral

      integ = 0.0_my_real
      do i=1,nqpoints
         integ = integ + qweights(i)*(kern(i) * u(1,1,i)**p * u(2,2,i)**q)
      end do

! free memory

      deallocate(kern,u,qweights,xquad,yquad,stat=astat)

   endif

endif

end function recur_integrate

!          -----------
subroutine phaml_scale(phaml_solution,factor,comp,eigen)
!          -----------

!----------------------------------------------------
! This routine scales the solution by multiplying by factor.  For
! multicomponent solutions, comp tells which component to scale.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(phaml_solution_type), intent(inout) :: phaml_solution
real(my_real), intent(in) :: factor
integer, intent(in), optional :: comp,eigen

!----------------------------------------------------
! Local variables

integer :: loc_comp, loc_eigen, my_processor, proc, lev, elem, edge, vert, i, &
           j
logical :: visited(size(phaml_solution%grid%edge))

!----------------------------------------------------
! Begin executable code

if (present(comp)) then
   loc_comp = comp
else
   loc_comp = 1
endif
if (present(eigen)) then
   loc_eigen = eigen
else
   loc_eigen = 1
endif
if (loc_comp < 1 .or. loc_comp > phaml_solution%system_size .or. &
    loc_eigen < 1 .or. loc_eigen > max(1,phaml_solution%grid%num_eval)) then
   call warning("phaml_scale: comp or eigen is not a valid solution number", &
                "Not scaling.")
   return
endif

my_processor = my_proc(phaml_solution%procs)

! Master tells the slaves to scale solution

if (my_processor == MASTER) then

! 13 is code for scale
   do proc=1,num_proc(phaml_solution%procs)
      call phaml_send(phaml_solution%procs,proc,(/13,loc_comp,loc_eigen/),3, &
                      (/factor/),1,101)
   end do

else ! slave

! traverse the elements, edges and vertices on each level, scaling the solution

   visited = .false.
   do lev=1,phaml_solution%grid%nlev
      elem = phaml_solution%grid%head_level_elem(lev)
      do while (elem /= END_OF_LIST)
         if (associated(phaml_solution%grid%element(elem)%solution)) then
            do j=1,size(phaml_solution%grid%element(elem)%solution,dim=1)
               phaml_solution%grid%element(elem)%solution(j,loc_comp,loc_eigen)=&
                  factor*phaml_solution%grid%element(elem)%solution(j,loc_comp,loc_eigen)
            end do
         endif
         do i=1,3
            edge = phaml_solution%grid%element(elem)%edge(i)
            if (.not. visited(edge) .and. &
                associated(phaml_solution%grid%edge(edge)%solution)) then
               do j=1,size(phaml_solution%grid%edge(edge)%solution,dim=1)
                  phaml_solution%grid%edge(edge)%solution(j,loc_comp,loc_eigen)=&
                     factor*phaml_solution%grid%edge(edge)%solution(j,loc_comp,loc_eigen)
               end do
               visited(edge) = .true.
            endif
         end do
         elem = phaml_solution%grid%element(elem)%next
      end do
      vert = phaml_solution%grid%head_level_vert(lev)
      do while (vert /= END_OF_LIST)
         phaml_solution%grid%vertex_solution(vert,loc_comp,loc_eigen) = &
            factor*phaml_solution%grid%vertex_solution(vert,loc_comp,loc_eigen)
         vert = phaml_solution%grid%vertex(vert)%next
      end do
   end do

   phaml_solution%grid%errind_up2date = .false.

endif

end subroutine phaml_scale

!          --------------
subroutine phaml_evaluate(phaml_solution,x,y,u,ux,uy,uxx,uyy,comp,eigen)
!          --------------

!----------------------------------------------------
! This routine evaluates the solution and derivatives at the points in the
! arrays (x,y) and returns the results in the arrays u*.  The return value
! is 0.0 for points that are outside the domain.  Any subset of the u's can
! be given.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type (phaml_solution_type), intent(in) :: phaml_solution
real(my_real), intent(in) :: x(:),y(:)
real(my_real), optional, intent(out) :: u(:),ux(:),uy(:),uxx(:),uyy(:)
integer, optional, intent(in) :: comp,eigen

!----------------------------------------------------
! Local variables:

integer :: loc_comp, loc_eigen
!----------------------------------------------------
! Begin executable code

if (size(y) /= size(x)) then
   call fatal("phaml_evaluate: x and y must have the same size.  They are ", &
              intlist=(/size(x),size(y)/))
   stop
endif
if (present(u)) then
   if (size(u) /= size(x)) then
      call fatal("phaml_evaluate: x and u must have the same size.  They are ", &
              intlist=(/size(x),size(u)/))
      stop
   endif
endif
if (present(ux)) then
   if (size(ux) /= size(x)) then
      call fatal("phaml_evaluate: x and ux must have the same size.  They are ", &
              intlist=(/size(x),size(ux)/))
      stop
   endif
endif
if (present(uy)) then
   if (size(uy) /= size(x)) then
      call fatal("phaml_evaluate: x and uy must have the same size.  They are ", &
              intlist=(/size(x),size(uy)/))
      stop
   endif
endif
if (present(uxx)) then
   if (size(uxx) /= size(x)) then
      call fatal("phaml_evaluate: x and uxx must have the same size.  They are ", &
              intlist=(/size(x),size(uxx)/))
      stop
   endif
endif
if (present(uyy)) then
   if (size(uyy) /= size(x)) then
      call fatal("phaml_evaluate: x and uyy must have the same size.  They are ", &
              intlist=(/size(x),size(uyy)/))
      stop
   endif
endif
if (present(comp)) then
   if (comp <= 0) then
      call fatal("phaml_evaluate: comp must be a positive integer.  It is ", &
                  intlist=(/comp/))
      stop
   endif
   if (comp > phaml_solution%system_size) then
      call fatal("phaml_evaluate: comp must be no larger than system size.  They are ", &
                 intlist=(/comp,phaml_solution%system_size/))
      stop
   endif
   loc_comp = comp
else
   loc_comp = 1
endif
if (present(eigen)) then
   if (eigen <= 0) then
      call fatal("phaml_evaluate: eigen must be a positive integer.  It is ", &
                 intlist=(/eigen/))
      stop
   endif
   if (eigen > max(1,phaml_solution%grid%num_eval)) then
      call fatal("phaml_evaluate: eigen must be no larger than the number of eigenvalues computed.  They are ", &
                 intlist=(/eigen,phaml_solution%grid%num_eval/))
      stop
   endif
   loc_eigen = eigen
else
   loc_eigen = 1
endif
if (PARALLEL == SEQUENTIAL) then
   call evaluate_soln_slave(phaml_solution%grid,phaml_solution%procs,x,y, &
                            loc_comp,loc_eigen,u,ux,uy,uxx,uyy)
else
   call evaluate_soln(phaml_solution%procs,x,y,u,ux,uy,uxx,uyy,loc_comp, &
                      loc_eigen)
endif

end subroutine phaml_evaluate

!          --------------
subroutine evaluate_slave(phaml_solution,x,y,invoker_procs,comp,eigen,which)
!          --------------

!----------------------------------------------------
! This is the slaves version of evaluate.  It does not have the
! solution argument because it sends the solution to the process that
! called subroutine evaluate, which in turn invoked the call to this
! routine by a sleeping process.  The invoking process sends its procs so
! I know who to send the return message to.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type (phaml_solution_type), intent(in) :: phaml_solution
real(my_real), intent(in) :: x(:),y(:)
type(proc_info), intent(in) :: invoker_procs
integer, intent(in) :: comp, eigen
integer, intent(in) :: which(5)

!----------------------------------------------------
! Begin executable code

call evaluate_soln_slave(phaml_solution%grid,invoker_procs,x,y,comp, &
                         eigen,which=which)

end subroutine evaluate_slave

!          ------------------
subroutine phaml_evaluate_old(x,y,u,ux,uy,uxx,uyy,comp,eigen)
!          ------------------

!----------------------------------------------------
! This routine returns the "old" solution and/or its derivatives
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

real(my_real), intent(in) :: x,y
real(my_real), intent(out), optional :: u, ux, uy, uxx, uyy
integer, intent(in), optional :: comp, eigen
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

! comp and eigen get checked in evaluate_oldsoln_local because this routine
! does not know how big they can be

call evaluate_oldsoln_local(x,y,u,ux,uy,uxx,uyy,comp,eigen)

end subroutine phaml_evaluate_old

!          ----------------------
subroutine phaml_copy_soln_to_old(phaml_solution)
!          ----------------------

!----------------------------------------------------
! This routine is the user interface for the routine that copies the solution
! in phaml_solution%grid to oldsoln
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(phaml_solution_type), intent(inout) :: phaml_solution
!----------------------------------------------------
! Local variables:

integer :: my_processor, proc
!----------------------------------------------------
! Begin executable code

my_processor = my_proc(phaml_solution%procs)

! Master tells the slaves to scale solution

if (my_processor == MASTER) then

! 16 is code for copy to old
   do proc=1,num_proc(phaml_solution%procs)
      call phaml_send(phaml_solution%procs,proc,(/16/),1,(/0.0_my_real/),0,101)
   end do

else ! slave

   call copy_old(phaml_solution%grid)

endif

end subroutine phaml_copy_soln_to_old

!          -----------
subroutine phaml_query(phaml_solution,nvert,nvert_proc,nvert_own,nelem, &
                       nelem_proc,nelem_own,neq,neq_proc,neq_own,nlev, &
                       min_degree,max_degree, linf_error,energy_error, &
                       l2_error,max_error_indicator, linf_error_estimate, &
                       energy_error_estimate,l2_error_estimate, linf_solution,&
                       l2_solution, energy_solution,linf_u,l2_u,energy_u, &
                       linf_true,l2_true,energy_true, &
                       eigenvalues,eigenvalue_error_estimate,max_linsys_resid, &
                       ave_linsys_resid,eigen_l2_resid,arpack_iter, &
                       arpack_nconv,arpack_numop,arpack_numopb, arpack_numreo, &
                       arpack_info,comp,eigen,error_estimator,eigen_variance)
!          -----------

!----------------------------------------------------
! This routine returns information about phaml_solution
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(phaml_solution_type), intent(inout), target :: phaml_solution
integer, intent(inout), optional :: nvert_proc(:)
integer, intent(out), optional :: nvert,nvert_own(:),nelem, &
                                  nelem_proc(:),nelem_own(:),neq, &
                                  neq_proc(:),neq_own(:),nlev, &
                                  min_degree,max_degree, &
                                  arpack_iter, arpack_nconv, arpack_numop, &
                                  arpack_numopb, arpack_numreo, arpack_info
real(my_real),intent(out), optional :: linf_error, energy_error, l2_error, &
                                  max_error_indicator, &
                                  linf_error_estimate, &
                                  energy_error_estimate,l2_error_estimate, &
                                  eigenvalue_error_estimate(:), &
                                  linf_solution, l2_solution, energy_solution, &
                                  linf_u, l2_u, energy_u, eigenvalues(:), &
                                  linf_true, l2_true, energy_true, &
                                  max_linsys_resid, ave_linsys_resid, &
                                  eigen_l2_resid(:), &
                                  eigen_variance(:)
integer, intent(in), optional :: comp, eigen, error_estimator
!----------------------------------------------------
! Local variables:

integer :: query_list(45)
integer :: my_processor, nproc
integer :: ni, nr, proc, i, j, p, lev, elem, neig, vert, which, loc_comp, &
           loc_eigen, loc_error_estimator, maxdeg
integer, allocatable :: send_int(:)
real(my_real), allocatable :: send_real(:)
integer, pointer :: recv_int(:)
real(my_real), pointer :: recv_real(:)
type(grid_type), pointer :: grid
!----------------------------------------------------
! Begin executable code

my_processor = my_proc(phaml_solution%procs)
nproc = num_proc(phaml_solution%procs)
grid => phaml_solution%grid

if (my_processor == MASTER .and. present(nvert_proc)) then
   if (size(nvert_proc) < nproc) then
      call fatal("phaml_query: size of nvert_proc must be at least the number of processors.  They are", &
                 intlist=(/size(nvert_proc),nproc/))
      stop 
   endif
endif

if (present(nvert_own)) then
   if (size(nvert_own) < nproc) then
      call fatal("phaml_query: size of nvert_own must be at least the number of processors.  They are", &
                 intlist=(/size(nvert_own),nproc/))
      stop 
   endif
endif

if (present(nelem_proc)) then
   if (size(nelem_proc) < nproc) then
      call fatal("phaml_query: size of nelem_proc must be at least the number of processors.  They are", &
                 intlist=(/size(nelem_proc),nproc/))
      stop 
   endif
endif

if (present(nelem_own)) then
   if (size(nelem_own) < nproc) then
      call fatal("phaml_query: size of nelem_own must be at least the number of processors.  They are", &
                 intlist=(/size(nelem_own),nproc/))
      stop 
   endif
endif

if (present(neq_proc)) then
   if (size(neq_proc) < nproc) then
      call fatal("phaml_query: size of neq_proc must be at least the number of processors.  They are", &
                 intlist=(/size(neq_proc),nproc/))
      stop 
   endif
endif

if (present(neq_own)) then
   if (size(neq_own) < nproc) then
      call fatal("phaml_query: size of neq_own must be at least the number of processors.  They are", &
                 intlist=(/size(neq_own),nproc/))
      stop 
   endif
endif

if ((present(eigenvalues) .or. present(max_linsys_resid) .or. &
    present(eigenvalue_error_estimate) .or. &
    present(ave_linsys_resid) .or. present(eigen_l2_resid) .or. &
    present(arpack_iter) .or. present(arpack_nconv) .or. &
    present(arpack_numop) .or. present(arpack_numopb) .or. &
    present(arpack_numreo) .or. present(arpack_info) .or. &
    present(eigen_variance)) .and. .not. associated(grid%eigenvalue)) then
   call warning("phaml_query: a query concerning eigenvalues was made,", &
                "but an eigenvalue problem was not solved.")
   return
endif

if (present(eigenvalues)) then
   if (size(eigenvalues) < max(1,grid%num_eval)) then
      call fatal("phaml_query: size of eigenvalues must be at least the number of eigenvalues computed.  They are ", &
                 intlist=(/size(eigenvalues),grid%num_eval/))
      stop
   endif
endif

if (present(eigenvalue_error_estimate)) then
   if (size(eigenvalue_error_estimate) < max(1,grid%num_eval)) then
      call fatal("phaml_query: size of eigenvalue_error_estimate must be at least the number of eigenvalues computed.  They are ", &
                 intlist=(/size(eigenvalue_error_estimate),grid%num_eval/))
      stop
   endif
endif

if (present(eigen_l2_resid)) then
   if (size(eigen_l2_resid) < max(1,grid%num_eval)) then
      call fatal("phaml_query: size of eigen_l2_resid must be at least the number of eigenvalues computed.  They are ", &
                 intlist=(/size(eigen_l2_resid),grid%num_eval/))
      stop
   endif
endif

if (present(comp)) then
   if (comp <= 0) then
      call fatal("phaml_query: comp must be a positive integer.  It is ", &
                  intlist=(/comp/))
      stop
   endif
   if (comp > grid%system_size) then
      call fatal("phaml_query: comp must be no larger than system size.  They are ", &
                 intlist=(/comp,grid%system_size/))
      stop
   endif
endif

if (present(eigen)) then
   if (eigen <= 0) then
      call fatal("phaml_query: eigen must be a positive integer.  It is ", &
                 intlist=(/eigen/))
      stop
   endif
   if (eigen > max(1,grid%num_eval)) then
      call fatal("phaml_query: eigen must be no larger than the number of eigenvalues computed.  They are ", &
                 intlist=(/eigen,grid%num_eval/))
      stop
   endif
endif

! might use a different error indicator than was used to solve the problem
! or used in a previous call to phaml_query

grid%errind_up2date = .false.

if (my_processor == MASTER) then

! Master

! Construct a list indicating which entities were queried

   query_list = 0
   query_list(1) = 14 ! code for query
   if (present(nvert))                 query_list(2) = 1
   if (present(nvert_proc))            query_list(3) = 1
   if (present(nvert_own))             query_list(4) = 1
   if (present(nelem))                 query_list(5) = 1
   if (present(nelem_proc))            query_list(6) = 1
   if (present(nelem_own))             query_list(7) = 1
   if (present(nlev))                  query_list(8) = 1
   if (present(linf_error))            query_list(9) = 1
   if (present(energy_error))          query_list(10) = 1
   if (present(max_error_indicator))   query_list(11) = 1
   if (present(linf_error_estimate))   query_list(12) = 1
   if (present(energy_error_estimate)) query_list(13) = 1
   if (present(linf_solution))         query_list(14) = 1
   if (present(l2_solution))           query_list(15) = 1
   if (present(eigenvalues))           query_list(16) = 1
   if (present(max_linsys_resid))      query_list(17) = 1
   if (present(ave_linsys_resid))      query_list(18) = 1
   if (present(eigen_l2_resid))        query_list(19) = 1
   if (present(arpack_iter))           query_list(20) = 1
   if (present(arpack_nconv))          query_list(21) = 1
   if (present(arpack_numop))          query_list(22) = 1
   if (present(arpack_numopb))         query_list(23) = 1
   if (present(arpack_numreo))         query_list(24) = 1
   if (present(arpack_info))           query_list(25) = 1
   if (present(eigen_variance))        query_list(26) = 1
   if (present(comp)) then
      query_list(27) = comp
   else
      query_list(27) = 1
   endif
   if (present(neq))                   query_list(28) = 1
   if (present(neq_proc))              query_list(29) = 1
   if (present(neq_own))               query_list(30) = 1
   if (present(l2_error))              query_list(31) = 1
!                                                          32 is available
   if (present(l2_error_estimate))     query_list(33) = 1
   if (present(linf_true))             query_list(34) = 1
   if (present(l2_true))               query_list(35) = 1
   if (present(energy_true))           query_list(36) = 1
   if (present(eigen)) then
      query_list(37) = eigen
   else
      query_list(37) = 1
   endif
   if (present(energy_solution))       query_list(38) = 1
   if (present(linf_u))                query_list(39) = 1
   if (present(l2_u))                  query_list(40) = 1
   if (present(energy_u))              query_list(41) = 1
   if (present(min_degree))            query_list(42) = 1
   if (present(max_degree))            query_list(43) = 1
   if (present(error_estimator)) then
      query_list(44) = error_estimator
   else
      query_list(44) = EXPLICIT_ERRIND
   endif
   if (present(eigenvalue_error_estimate)) query_list(45) = 1

! send list to slaves

   do proc=1,nproc
      call phaml_send(phaml_solution%procs,proc,query_list,size(query_list), &
                      (/0.0_my_real/),0,101)
   end do

! initialize results

   if (present(nvert)) nvert = 0
   if (present(nvert_proc)) nvert_proc = 0
   if (present(nvert_own)) nvert_own = 0
   if (present(nelem)) nelem = 0
   if (present(nelem_proc)) nelem_proc = 0
   if (present(nelem_own)) nelem_own = 0
   if (present(neq)) neq = 0
   if (present(neq_proc)) neq_proc = 0
   if (present(neq_own)) neq_own = 0
   if (present(nlev)) nlev = 0
   if (present(linf_error)) linf_error = 0.0_my_real
   if (present(energy_error)) energy_error = 0.0_my_real
   if (present(l2_error)) l2_error = 0.0_my_real
   if (present(max_error_indicator)) max_error_indicator = 0.0_my_real
   if (present(linf_error_estimate)) linf_error_estimate = 0.0_my_real
   if (present(energy_error_estimate)) energy_error_estimate = 0.0_my_real
   if (present(l2_error_estimate)) l2_error_estimate = 0.0_my_real
   if (present(eigenvalue_error_estimate)) eigenvalue_error_estimate=0.0_my_real
   if (present(linf_solution)) linf_solution = 0.0_my_real
   if (present(l2_solution)) l2_solution = 0.0_my_real
   if (present(eigenvalues)) eigenvalues = 0.0_my_real
   if (present(max_linsys_resid)) max_linsys_resid = 0.0_my_real
   if (present(ave_linsys_resid)) ave_linsys_resid = 0.0_my_real
   if (present(eigen_l2_resid)) eigen_l2_resid = 0.0_my_real
   if (present(arpack_iter)) arpack_iter = 0
   if (present(arpack_nconv)) arpack_nconv = 0
   if (present(arpack_numop)) arpack_numop = 0
   if (present(arpack_numopb)) arpack_numopb = 0
   if (present(arpack_numreo)) arpack_numreo = 0
   if (present(arpack_info)) arpack_info = 0
   if (present(eigen_variance)) eigen_variance = 0.0_my_real
   if (present(linf_true)) linf_true = 0.0_my_real
   if (present(l2_true)) l2_true = 0.0_my_real
   if (present(energy_true)) energy_true = 0.0_my_real
   if (present(energy_solution)) energy_solution = 0.0_my_real
   if (present(linf_u)) linf_u = 0.0_my_real
   if (present(l2_u)) l2_u = 0.0_my_real
   if (present(energy_u)) energy_u = 0.0_my_real
   if (present(min_degree)) min_degree = huge(0)
   if (present(max_degree)) max_degree = 0

! receive partial results from slaves and build results

   do p=1,nproc
      call phaml_recv(phaml_solution%procs,proc,recv_int,ni,recv_real,nr,350)
      i = 0; j = 0
      if (present(nvert)) then
         i = i + 1
         if (.not. phaml_solution%still_sequential .or. p==1) then
            nvert = nvert + recv_int(i)
         endif
      endif
      if (present(nvert_proc)) then
         i = i + 1
         nvert_proc(proc) = recv_int(i)
      endif
      if (present(nvert_own)) then
         i = i + 1
         nvert_own(proc) = recv_int(i)
      endif
      if (present(nelem)) then
         i = i + 1
         if (.not. phaml_solution%still_sequential .or. p==1) then
            nelem = nelem + recv_int(i)
         endif
      endif
      if (present(nelem_proc)) then
         i = i + 1
         nelem_proc(proc) = recv_int(i)
      endif
      if (present(nelem_own)) then
         i = i + 1
         nelem_own(proc) = recv_int(i)
      endif
      if (present(neq)) then
         i = i + 1
         if (.not. phaml_solution%still_sequential .or. p==1) then
            neq = neq + recv_int(i)
         endif
      endif
      if (present(neq_proc)) then
         i = i + 1
         neq_proc(proc) = recv_int(i)
      endif
      if (present(neq_own)) then
         i = i + 1
         neq_own(proc) = recv_int(i)
      endif
      if (present(nlev)) then
         i = i + 1
         nlev = max(nlev,recv_int(i))
      endif
      if (present(linf_error)) then
         j = j + 1
         linf_error = max(linf_error,recv_real(j))
      endif
      if (present(energy_error)) then
         j = j + 1
         energy_error = recv_real(j)
      endif
      if (present(l2_error)) then
         j = j + 1
         if (.not. phaml_solution%still_sequential .or. p==1) then
            l2_error = l2_error + recv_real(j)**2
         endif
      endif
      if (present(max_error_indicator)) then
         j = j + 1
         max_error_indicator = max(max_error_indicator,recv_real(j))
      endif
      if (present(linf_error_estimate)) then
         j = j + 1
         linf_error_estimate = max(linf_error_estimate,recv_real(j))
      endif
      if (present(energy_error_estimate)) then
         j = j + 1
         if (.not. phaml_solution%still_sequential .or. p==1) then
            energy_error_estimate = energy_error_estimate + recv_real(j)**2
         endif
      endif
      if (present(l2_error_estimate)) then
         j = j + 1
         if (.not. phaml_solution%still_sequential .or. p==1) then
            l2_error_estimate = l2_error_estimate + recv_real(j)**2
         endif
      endif
      if (present(eigenvalue_error_estimate)) then
         if (.not. phaml_solution%still_sequential .or. p==1) then
            neig = min(size(eigenvalue_error_estimate),size(grid%eigenvalue))
            eigenvalue_error_estimate(1:neig) = &
               eigenvalue_error_estimate(1:neig) + recv_real(j+1:j+neig)**2
         endif
         j = j + size(grid%eigenvalue)
      endif
      if (present(linf_solution)) then
         j = j + 1
         linf_solution = max(linf_solution,recv_real(j))
      endif
      if (present(l2_solution)) then
         j = j + 1
         if (.not. phaml_solution%still_sequential .or. p==1) then
            l2_solution = l2_solution + recv_real(j)
         endif
      endif
      if (present(eigenvalues)) then
         neig = min(size(eigenvalues),size(grid%eigenvalue))
         eigenvalues(1:neig) = recv_real(j+1:j+neig)
         j = j + size(grid%eigenvalue)
      endif
      if (present(max_linsys_resid)) then
         j = j + 1
         max_linsys_resid = recv_real(j)
      endif
      if (present(ave_linsys_resid)) then
         j = j + 1
         ave_linsys_resid = recv_real(j)
      endif
      if (present(eigen_l2_resid)) then
         neig = min(size(eigen_l2_resid),size(grid%eigenvalue))
         eigen_l2_resid(1:neig) = recv_real(j+1:j+neig)
         j = j + size(grid%eigenvalue)
      endif
      if (present(arpack_iter)) then
         i = i + 1
         arpack_iter = recv_int(i)
      endif
      if (present(arpack_nconv)) then
         i = i + 1
         arpack_nconv = recv_int(i)
      endif
      if (present(arpack_numop)) then
         i = i + 1
         arpack_numop = recv_int(i)
      endif
      if (present(arpack_numopb)) then
         i = i + 1
         arpack_numopb = recv_int(i)
      endif
      if (present(arpack_numreo)) then
         i = i + 1
         arpack_numreo = recv_int(i)
      endif
      if (present(arpack_info)) then
         i = i + 1
         arpack_info = recv_int(i)
      endif
      if (present(eigen_variance)) then
         neig = min(size(eigen_variance),size(grid%eigenvalue))
         eigen_variance(1:neig) = recv_real(j+1:j+neig)
         j = j + size(grid%eigenvalue)
      endif
      if (present(linf_true)) then
         j = j + 1
         linf_true = max(linf_true,recv_real(j))
      endif
      if (present(l2_true)) then
         j = j + 1
         if (.not. phaml_solution%still_sequential .or. p==1) then
            l2_true = l2_true + recv_real(j)
         endif
      endif
      if (present(energy_true)) then
         j = j + 1
         if (.not. phaml_solution%still_sequential .or. p==1) then
            energy_true = energy_true + recv_real(j)
         endif
      endif
      if (present(energy_solution)) then
         j = j + 1
         if (.not. phaml_solution%still_sequential .or. p==1) then
            energy_solution = energy_solution + recv_real(j)
         endif
      endif
      if (present(linf_u)) then
         j = j + 1
         linf_u = max(linf_u,recv_real(j))
      endif
      if (present(l2_u)) then
         j = j + 1
         if (.not. phaml_solution%still_sequential .or. p==1) then
            l2_u = l2_u + recv_real(j)
         endif
      endif
      if (present(energy_u)) then
         j = j + 1
         if (.not. phaml_solution%still_sequential .or. p==1) then
            energy_u = energy_u + recv_real(j)
         endif
      endif
      if (present(min_degree)) then
         i = i + 1
         if (.not. phaml_solution%still_sequential .or. p==1) then
            min_degree = min(min_degree,recv_int(i))
         endif
      endif
      if (present(max_degree)) then
         i = i + 1
         if (.not. phaml_solution%still_sequential .or. p==1) then
            max_degree = max(max_degree,recv_int(i))
         endif
      endif

      if (associated(recv_int)) deallocate(recv_int)
      if (associated(recv_real)) deallocate(recv_real)

   end do

   if (present(energy_error_estimate)) &
      energy_error_estimate = sqrt(energy_error_estimate)
   if (present(l2_error_estimate)) &
      l2_error_estimate = sqrt(l2_error_estimate)
   if (present(eigenvalue_error_estimate)) &
      eigenvalue_error_estimate = sqrt(eigenvalue_error_estimate)
   if (present(l2_error)) l2_error = sqrt(l2_error)
   if (present(l2_solution)) l2_solution = sqrt(l2_solution)
   if (present(l2_true)) l2_true = sqrt(l2_true)
   if (present(energy_true)) energy_true = sqrt(energy_true)
   if (present(energy_solution)) energy_solution = sqrt(energy_solution)
   if (present(l2_u)) l2_u = sqrt(l2_u)
   if (present(energy_u)) energy_u = sqrt(energy_u)

else

! Slave

! The query list is passed in through nvert_proc

   if (PARALLEL /= SEQUENTIAL) then

      query_list = nvert_proc

   else

! If running as a sequential program, construct the query list that the
! master would have sent in a parallel program

      query_list = 0
      query_list(1) = 14 ! code for query
      if (present(nvert))                 query_list(2) = 1
      if (present(nvert_proc))            query_list(3) = 1
      if (present(nvert_own))             query_list(4) = 1
      if (present(nelem))                 query_list(5) = 1
      if (present(nelem_proc))            query_list(6) = 1
      if (present(nelem_own))             query_list(7) = 1
      if (present(nlev))                  query_list(8) = 1
      if (present(linf_error))            query_list(9) = 1
      if (present(energy_error))          query_list(10) = 1
      if (present(max_error_indicator))   query_list(11) = 1
      if (present(linf_error_estimate))   query_list(12) = 1
      if (present(energy_error_estimate)) query_list(13) = 1
      if (present(linf_solution))         query_list(14) = 1
      if (present(l2_solution))           query_list(15) = 1
      if (present(eigenvalues))           query_list(16) = 1
      if (present(max_linsys_resid))      query_list(17) = 1
      if (present(ave_linsys_resid))      query_list(18) = 1
      if (present(eigen_l2_resid))        query_list(19) = 1
      if (present(arpack_iter))           query_list(20) = 1
      if (present(arpack_nconv))          query_list(21) = 1
      if (present(arpack_numop))          query_list(22) = 1
      if (present(arpack_numopb))         query_list(23) = 1
      if (present(arpack_numreo))         query_list(24) = 1
      if (present(arpack_info))           query_list(25) = 1
      if (present(eigen_variance))        query_list(26) = 1
      if (present(comp)) then
         query_list(27) = comp
      else
         query_list(27) = 1
      endif
      if (present(neq))                   query_list(28) = 1
      if (present(neq_proc))              query_list(29) = 1
      if (present(neq_own))               query_list(30) = 1
      if (present(l2_error))              query_list(31) = 1
      if (present(l2_error_estimate))     query_list(33) = 1
      if (present(linf_true))             query_list(34) = 1
      if (present(l2_true))               query_list(35) = 1
      if (present(energy_true))           query_list(36) = 1
      if (present(eigen)) then
         query_list(37) = eigen
      else
         query_list(37) = 1
      endif
      if (present(energy_solution))       query_list(38) = 1
      if (present(linf_u))                query_list(39) = 1
      if (present(l2_u))                  query_list(40) = 1
      if (present(energy_u))              query_list(41) = 1
      if (present(min_degree))            query_list(42) = 1
      if (present(max_degree))            query_list(43) = 1
      if (present(error_estimator)) then
         query_list(44) = error_estimator
      else
         query_list(44) = EXPLICIT_ERRIND
      endif
      if (present(eigenvalue_error_estimate)) query_list(45) = 1

   endif
   
   loc_comp = query_list(27)
   loc_eigen = query_list(37)
   which = loc_comp + (loc_eigen-1)*grid%system_size
   loc_error_estimator = query_list(44)

! if the error estimator is equilibrated residual, make sure the edge mass
! matrices have been set

   if (loc_error_estimator == EQUILIBRATED_RESIDUAL) then
      call get_grid_info(grid,phaml_solution%procs, &
                         phaml_solution%still_sequential,1, &
                         maxdeg=maxdeg)
      maxdeg = phaml_global_max(phaml_solution%procs,maxdeg,351)
      call set_edge_mass(maxdeg+1)
   endif

! if master had max_error_indicator present
   if (query_list(11) == 1) then
      call all_error_indicators(grid,loc_error_estimator)
   endif

! Count the number of responses

   ni = 0; nr = 0
   if (query_list(2) == 1) ni = ni + 1
   if (query_list(3) == 1) ni = ni + 1
   if (query_list(4) == 1) ni = ni + 1
   if (query_list(5) == 1) ni = ni + 1
   if (query_list(6) == 1) ni = ni + 1
   if (query_list(7) == 1) ni = ni + 1
   if (query_list(8) == 1) ni = ni + 1
   if (query_list(9) == 1) nr = nr + 1
   if (query_list(10) == 1) nr = nr + 1
   if (query_list(11) == 1) nr = nr + 1
   if (query_list(12) == 1) nr = nr + 1
   if (query_list(13) == 1) nr = nr + 1
   if (query_list(14) == 1) nr = nr + 1
   if (query_list(15) == 1) nr = nr + 1
   if (query_list(16) == 1) nr = nr + size(grid%eigenvalue)
   if (query_list(17) == 1) nr = nr + 1
   if (query_list(18) == 1) nr = nr + 1
   if (query_list(19) == 1) nr = nr + size(grid%eigenvalue)
   if (query_list(20) == 1) ni = ni + 1
   if (query_list(21) == 1) ni = ni + 1
   if (query_list(22) == 1) ni = ni + 1
   if (query_list(23) == 1) ni = ni + 1
   if (query_list(24) == 1) ni = ni + 1
   if (query_list(25) == 1) ni = ni + 1
   if (query_list(26) == 1) nr = nr + size(grid%eigenvalue)
   if (query_list(28) == 1) ni = ni + 1
   if (query_list(29) == 1) ni = ni + 1
   if (query_list(30) == 1) ni = ni + 1
   if (query_list(31) == 1) nr = nr + 1
   if (query_list(33) == 1) nr = nr + 1
   if (query_list(34) == 1) nr = nr + 1
   if (query_list(35) == 1) nr = nr + 1
   if (query_list(36) == 1) nr = nr + 1
   if (query_list(38) == 1) nr = nr + 1
   if (query_list(39) == 1) nr = nr + 1
   if (query_list(40) == 1) nr = nr + 1
   if (query_list(41) == 1) nr = nr + 1
   if (query_list(42) == 1) ni = ni + 1
   if (query_list(43) == 1) ni = ni + 1
   if (query_list(45) == 1) nr = nr + size(grid%eigenvalue)

   if (ni > 0) then
      allocate(send_int(ni))
   else
      allocate(send_int(1))
   endif
   if (nr > 0) then
      allocate(send_real(nr))
   else
      allocate(send_real(1))
   endif

! Build the response

   i = 0; j = 0
   if (query_list(2) == 1) then
      i = i + 1
      send_int(i) = grid%nvert_own
   endif
   if (query_list(3) == 1) then
      i = i + 1
      send_int(i) = grid%nvert
   endif
   if (query_list(4) == 1) then
      i = i + 1
      send_int(i) = grid%nvert_own
   endif
   if (query_list(5) == 1) then
      i = i + 1
      send_int(i) = grid%nelem_leaf_own
   endif
   if (query_list(6) == 1) then
      i = i + 1
      send_int(i) = grid%nelem_leaf
   endif
   if (query_list(7) == 1) then
      i = i + 1
      send_int(i) = grid%nelem_leaf_own
   endif
   if (query_list(28) == 1) then
      i = i + 1
      send_int(i) = grid%dof_own
   endif
   if (query_list(29) == 1) then
      i = i + 1
      send_int(i) = grid%dof
   endif
   if (query_list(30) == 1) then
      i = i + 1
      send_int(i) = grid%dof_own
   endif
   if (query_list(8) == 1) then
      i = i + 1
      send_int(i) = grid%nlev
   endif
   if (query_list(9) == 1) then
      j = j + 1
      call norm_error(phaml_solution%grid,phaml_solution%procs, &
                      phaml_solution%still_sequential,loc_comp,loc_eigen, &
                      my_Linf_norm=send_real(j))
   endif
   if (query_list(10) == 1) then
      j = j + 1
      call norm_error(phaml_solution%grid,phaml_solution%procs, &
                      phaml_solution%still_sequential,loc_comp,loc_eigen,  &
                      energy_norm=send_real(j))
   endif
   if (query_list(31) == 1) then
      j = j + 1
      call norm_error(phaml_solution%grid,phaml_solution%procs, &
                      phaml_solution%still_sequential,loc_comp,loc_eigen,  &
                      my_L2_norm=send_real(j))
   endif
   if (query_list(11) == 1) then
      j = j + 1
      send_real(j) = 0.0_my_real
      do lev = 1,grid%nlev
         elem = grid%head_level_elem(lev)
         do while (elem /= END_OF_LIST)
            if (grid%element(elem)%iown) then
               if (grid%element(elem)%isleaf) then
                  send_real(j) = max(send_real(j), &
           maxval(grid%element_errind(elem,:))/grid%element(elem)%work)
               endif
            endif
            elem = grid%element(elem)%next
         end do
      end do
   endif
   if (query_list(12) == 1) then
      j = j + 1
      call error_estimate(grid,phaml_solution%procs,loc_error_estimator,which, &
                          errest_Linf=send_real(j))
   endif
   if (query_list(13) == 1) then
      j = j + 1
      call error_estimate(grid,phaml_solution%procs,loc_error_estimator,which, &
                          errest_energy=send_real(j))
   endif
   if (query_list(33) == 1) then
      j = j + 1
      call error_estimate(grid,phaml_solution%procs,loc_error_estimator,which, &
                          errest_L2=send_real(j))
   endif
   if (query_list(45) == 1) then
      do p=1,size(grid%eigenvalue)
         j = j + 1
         call error_estimate(grid,phaml_solution%procs,loc_error_estimator,p, &
                             errest_eigenvalue=send_real(j))
      end do
   endif
   if (query_list(14) == 1) then
      j = j + 1
      send_real(j) = 0.0_my_real
      do lev=1,grid%nlev
         vert = grid%head_level_vert(lev)
         do while (vert /= END_OF_LIST)
            if (grid%element(grid%vertex(vert)%assoc_elem)%iown) then
               send_real(j) = max(send_real(j),abs(grid%vertex_solution(vert,loc_comp,loc_eigen)))
            endif
            vert = grid%vertex(vert)%next
         end do
      end do
   endif
   if (query_list(15) == 1) then
      j = j + 1
      send_real(j) = 0.0_my_real
      do lev=1,grid%nlev
         vert = grid%head_level_vert(lev)
         do while (vert /= END_OF_LIST)
            if (grid%element(grid%vertex(vert)%assoc_elem)%iown) then
               send_real(j) = send_real(j) + grid%vertex_solution(vert,loc_comp,loc_eigen)**2
            endif
            vert = grid%vertex(vert)%next
         end do
      end do
   endif
   if (query_list(16) == 1) then
      send_real(j+1:j+size(grid%eigenvalue)) = grid%eigenvalue
      j = j + size(grid%eigenvalue)
   endif
   if (query_list(17) == 1) then
      j = j + 1
      send_real(j) = grid%eigen_linsys_max_l2_resid
   endif
   if (query_list(18) == 1) then
      j = j + 1
      send_real(j) = grid%eigen_linsys_ave_l2_resid
   endif
   if (query_list(19) == 1) then
      send_real(j+1:j+size(grid%eigenvalue)) = grid%eigenprob_l2_resid
      j = j + size(grid%eigenvalue)
   endif
   if (query_list(20) == 1) then
      i = i + 1
      send_int(i) = grid%arpack_iter
   endif
   if (query_list(21) == 1) then
      i = i + 1
      send_int(i) = grid%arpack_nconv
   endif
   if (query_list(22) == 1) then
      i = i + 1
      send_int(i) = grid%arpack_numop
   endif
   if (query_list(23) == 1) then
      i = i + 1
      send_int(i) = grid%arpack_numopb
   endif
   if (query_list(24) == 1) then
      i = i + 1
      send_int(i) = grid%arpack_numreo
   endif
   if (query_list(25) == 1) then
      i = i + 1
      send_int(i) = grid%arpack_info
   endif
   if (query_list(26) == 1) then
      send_real(j+1:j+size(grid%eigenvalue)) = grid%eigenprob_variance
      j = j + size(grid%eigenvalue)
   endif
   if (query_list(34) == 1) then
      j = j + 1
      call norm_true(grid,phaml_solution%procs, &
                     phaml_solution%still_sequential, &
                     query_list(27),query_list(37), &
                     linf=send_real(j))
   endif
   if (query_list(35) == 1) then
      j = j + 1
      call norm_true(grid,phaml_solution%procs, &
                     phaml_solution%still_sequential, &
                     query_list(27),query_list(37), &
                     l2=send_real(j))
      send_real(j) = send_real(j)**2
   endif
   if (query_list(36) == 1) then
      j = j + 1
      call norm_true(grid,phaml_solution%procs, &
                     phaml_solution%still_sequential, &
                     query_list(27),query_list(37), &
                     energy=send_real(j))
      send_real(j) = send_real(j)**2
   endif
   if (query_list(38) == 1) then
      j = j + 1
      call norm_solution(grid,phaml_solution%procs, &
                         phaml_solution%still_sequential, &
                         query_list(27),query_list(37), &
                         discrete_energy=send_real(j))
      send_real(j) = send_real(j)**2
   endif
   if (query_list(39) == 1) then
      j = j + 1
      call norm_solution(grid,phaml_solution%procs, &
                         phaml_solution%still_sequential, &
                         query_list(27),query_list(37), &
                         linf=send_real(j))
   endif
   if (query_list(40) == 1) then
      j = j + 1
      call norm_solution(grid,phaml_solution%procs, &
                         phaml_solution%still_sequential, &
                         query_list(27),query_list(37), &
                         l2=send_real(j))
      send_real(j) = send_real(j)**2
   endif
   if (query_list(41) == 1) then
      j = j + 1
      call norm_solution(grid,phaml_solution%procs, &
                         phaml_solution%still_sequential, &
                         query_list(27),query_list(37), &
                         energy=send_real(j))
      send_real(j) = send_real(j)**2
   endif
   if (query_list(42) == 1) then
      i = i + 1
      call get_grid_info(grid,phaml_solution%procs, &
                         phaml_solution%still_sequential,1, &
                         mindeg=send_int(i))
   endif
   if (query_list(43) == 1) then
      i = i + 1
      call get_grid_info(grid,phaml_solution%procs, &
                         phaml_solution%still_sequential,1, &
                         maxdeg=send_int(i))
   endif

! send the results to the master

   if (PARALLEL /= SEQUENTIAL) then

      call phaml_send(phaml_solution%procs,MASTER,send_int,ni,send_real,nr,350)

   else

! if running as a sequential program, copy return values from the message
! that would have been sent to the master in a parallel program

      i = 0; j = 0; proc = 1
      if (present(nvert)) then
         i = i + 1
         nvert = send_int(i)
      endif
      if (present(nvert_proc)) then
         i = i + 1
         nvert_proc(proc) = send_int(i)
      endif
      if (present(nvert_own)) then
         i = i + 1
         nvert_own(proc) = send_int(i)
      endif
      if (present(nelem)) then
         i = i + 1
         nelem = send_int(i)
      endif
      if (present(nelem_proc)) then
         i = i + 1
         nelem_proc(proc) = send_int(i)
      endif
      if (present(nelem_own)) then
         i = i + 1
         nelem_own(proc) = send_int(i)
      endif
      if (present(neq)) then
         i = i + 1
         neq = send_int(i)
      endif
      if (present(neq_proc)) then
         i = i + 1
         neq_proc(proc) = send_int(i)
      endif
      if (present(neq_own)) then
         i = i + 1
         neq_own(proc) = send_int(i)
      endif
      if (present(nlev)) then
         i = i + 1
         nlev = send_int(i)
      endif
      if (present(linf_error)) then
         j = j + 1
         linf_error = send_real(j)
      endif
      if (present(energy_error)) then
         j = j + 1
         energy_error = send_real(j)
      endif
      if (present(l2_error)) then
         j = j + 1
         l2_error = send_real(j)
      endif
      if (present(max_error_indicator)) then
         j = j + 1
         max_error_indicator = send_real(j)
      endif
      if (present(linf_error_estimate)) then
         j = j + 1
         linf_error_estimate = send_real(j)
      endif
      if (present(energy_error_estimate)) then
         j = j + 1
         energy_error_estimate = send_real(j)
      endif
      if (present(l2_error_estimate)) then
         j = j + 1
         l2_error_estimate = send_real(j)
      endif
      if (present(eigenvalue_error_estimate)) then
         neig = min(size(eigenvalue_error_estimate),size(grid%eigenvalue))
         eigenvalue_error_estimate(1:neig) = send_real(j+1:j+neig)**2
         j = j + size(grid%eigenvalue)
      endif
      if (present(linf_solution)) then
         j = j + 1
         linf_solution = send_real(j)
      endif
      if (present(l2_solution)) then
         j = j + 1
         l2_solution = sqrt(send_real(j))
      endif
      if (present(eigenvalues)) then
         neig = min(size(eigenvalues),size(grid%eigenvalue))
         eigenvalues(1:neig) = send_real(j+1:j+neig)
         j = j + size(grid%eigenvalue)
      endif
      if (present(max_linsys_resid)) then
         j = j + 1
         max_linsys_resid = send_real(j)
      endif
      if (present(ave_linsys_resid)) then
         j = j + 1
         ave_linsys_resid = send_real(j)
      endif
      if (present(eigen_l2_resid)) then
         neig = min(size(eigen_l2_resid),size(grid%eigenvalue))
         eigen_l2_resid(1:neig) = send_real(j+1:j+neig)
         j = j + size(grid%eigenvalue)
      endif
      if (present(arpack_iter)) then
         i = i + 1
         arpack_iter = send_int(i)
      endif
      if (present(arpack_nconv)) then
         i = i + 1
         arpack_nconv = send_int(i)
      endif
      if (present(arpack_numop)) then
         i = i + 1
         arpack_numop = send_int(i)
      endif
      if (present(arpack_numopb)) then
         i = i + 1
         arpack_numopb = send_int(i)
      endif
      if (present(arpack_numreo)) then
         i = i + 1
         arpack_numreo = send_int(i)
      endif
      if (present(arpack_info)) then
         i = i + 1
         arpack_info = send_int(i)
      endif
      if (present(eigen_variance)) then
         neig = min(size(eigen_variance),size(grid%eigenvalue))
         eigen_variance(1:neig) = send_real(j+1:j+neig)
         j = j + size(grid%eigenvalue)
      endif
      if (present(linf_true)) then
         j = j + 1
         linf_true = send_real(j)
      endif
      if (present(l2_true)) then
         j = j + 1
         l2_true = sqrt(send_real(j))
      endif
      if (present(energy_true)) then
         j = j + 1
         energy_true = sqrt(send_real(j))
      endif
      if (present(energy_solution)) then
         j = j + 1
         energy_solution = sqrt(send_real(j))
      endif
      if (present(linf_u)) then
         j = j + 1
         linf_u = send_real(j)
      endif
      if (present(l2_u)) then
         j = j + 1
         l2_u = sqrt(send_real(j))
      endif
      if (present(energy_u)) then
         j = j + 1
         energy_u = sqrt(send_real(j))
      endif
      if (present(min_degree)) then
         i = i + 1
         min_degree = send_int(i)
      endif
      if (present(max_degree)) then
         i = i + 1
         max_degree = send_int(i)
      endif

   endif

   deallocate(send_int,send_real)
endif

end subroutine phaml_query

!          -------------
subroutine phaml_connect(pde1,pde2)
!          -------------

!----------------------------------------------------
! This routine "connects" two PDEs by providing the slaves associated
! with each one with the proc_info for the other one.  pde1 and pde2
! are indices into the array pde.
! This is called by the master; the slaves must be waiting to receive a
! job from the master.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer, intent(in) :: pde1, pde2

!----------------------------------------------------
! Local variables

integer, allocatable :: send_int(:)
real(my_real), allocatable :: send_real(:)
integer :: ni,nr,first_int,first_real,proc,astat

!----------------------------------------------------
! Begin executable code

! check validity of arguments

if (pde1 < 1 .or. pde2 < 1 .or. pde1 > size(pde) .or. pde2 > size(pde)) then
   call fatal("indices sent to subroutine phaml_connect must be >=1 and <= size(pde)",&
              "pde1, pde2 and size(pde) are",intlist=(/pde1,pde2,size(pde)/),procs=pde(1)%procs)
   stop
endif

! check that I am the master; this may also catch some cases where the
! pde's have not been created

if (my_proc(pde(pde1)%procs) /= MASTER .or. &
    my_proc(pde(pde2)%procs) /= MASTER) then
   call fatal("phaml_connect was called by a process that is not the master or with an argument that is not a created pde", &
              procs=pde(pde1)%procs)
   stop
endif

! send pde1 info to slaves for pde2

call pack_procs_size(pde(pde1)%procs,ni,nr)
ni = ni+2
allocate(send_int(ni),send_real(nr),stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in phaml_connect",procs=pde(pde1)%procs)
   return
endif
send_int(1) = 4 ! job code for connect
send_int(2) = pde1
first_int = 3
first_real = 1
call pack_procs(pde(pde1)%procs,send_int,first_int,send_real,first_real)
do proc=1,num_proc(pde(pde2)%procs)
   call phaml_send(pde(pde2)%procs,proc,send_int,ni,send_real,nr,101)
end do
deallocate(send_int,send_real,stat=astat)

! send pde2 info to slaves for pde1 if pde2/=pde1

if (pde2 /= pde1) then
   call pack_procs_size(pde(pde2)%procs,ni,nr)
   ni = ni+2
   allocate(send_int(ni),send_real(nr),stat=astat)
   if (astat /= 0) then
      ierr = ALLOC_FAILED
      call fatal("allocation failed in phaml_connect",procs=pde(pde1)%procs)
      return
   endif
   send_int(1) = 4 ! job code for connect
   send_int(2) = pde2
   first_int = 3
   first_real = 1
   call pack_procs(pde(pde2)%procs,send_int,first_int,send_real,first_real)
   do proc=1,num_proc(pde(pde1)%procs)
      call phaml_send(pde(pde1)%procs,proc,send_int,ni,send_real,nr,101)
   end do
   deallocate(send_int,send_real,stat=astat)
endif

end subroutine phaml_connect

!          -----------
subroutine phaml_store(phaml_solution,unit)
!          -----------

!----------------------------------------------------
! This routine stores information from phaml_solution into files for later use.
! unit is the unit number to write to, which should have been opened with
! phaml_popen.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type (phaml_solution_type), intent(in), target :: phaml_solution
integer, intent(in) :: unit
!----------------------------------------------------
! Local variables:

integer :: send_int(2)
real(my_real) :: send_real(1)
integer :: i,j,proc,nproc,me,sub,lev,elem,vert,edge,elem_dim_save, &
           vert_dim_save,edge_dim_save
type(grid_type), pointer :: grid
character(len=11) :: form
!----------------------------------------------------
! Begin executable code

outunit = phaml_solution%outunit
errunit = phaml_solution%errunit

inquire(unit,form=form)
if (form == "UNDEFINED") then
   ierr = USER_INPUT_ERROR
   call fatal("The file for phaml_store must be opened as either FORMATTED or UNFORMATTED", &
              "before calling phaml_store",procs=phaml_solution%procs)
   return
endif

me = my_proc(phaml_solution%procs)
nproc = num_proc(phaml_solution%procs)

! if I am the master, tell the slaves to save the data

if (me == MASTER) then
   send_int(1) = 5 ! code for saving solution
   send_int(2) = unit
   do proc=1,nproc
      call phaml_send(phaml_solution%procs,proc,send_int,2,send_real,0,101)
   end do
endif

! select between UNFORMATTED and FORMATTED.  The two segments of code
! should be identical except for write(unit) vs. write(unit,*)

if (form == 'UNFORMATTED') then

   write(unit) 'PHAML data'
   write(unit) 12 ! version number

! data from phaml_solution_type

   write(unit) phaml_solution%still_sequential
   write(unit) phaml_solution%pde_id
   write(unit) phaml_solution%system_size

! data from proc_info

   write(unit) nproc,me

! data from the grid

   grid => phaml_solution%grid
   if (me == MASTER) allocate(grid%initial_neighbor(1,1))
   elem_dim_save = 0
   vert_dim_save = 0
   edge_dim_save = 0
   do lev=1,grid%nlev
      elem = grid%head_level_elem(lev)
      do while (elem /= END_OF_LIST)
         elem_dim_save = max(elem,elem_dim_save)
         do vert = 1,VERTICES_PER_ELEMENT
            vert_dim_save = max(grid%element(elem)%vertex(vert),vert_dim_save)
         end do
         do edge = 1,EDGES_PER_ELEMENT
            edge_dim_save = max(grid%element(elem)%edge(edge),edge_dim_save)
         end do
         elem = grid%element(elem)%next
      end do
   end do
   if (associated(grid%eigenvalue)) then
      write(unit) vert_dim_save,elem_dim_save,edge_dim_save, &
                  size(grid%initial_neighbor,dim=2), &
                  size(grid%head_level_elem),size(grid%head_level_vert), &
                  size(grid%eigenvalue),grid%nsoln
   else
      write(unit) vert_dim_save,elem_dim_save,edge_dim_save, &
                  size(grid%initial_neighbor,dim=2), &
                  size(grid%head_level_elem),size(grid%head_level_vert), &
                  -1,grid%nsoln
   endif
   write(unit) size(grid%vertex_solution,2),size(grid%vertex_solution,3)
   if (grid%have_true) then
      write(unit) 1
   else
      write(unit) 0
   endif
   write(unit) EDGES_PER_ELEMENT
   elem = grid%next_free_elem
   if (elem > elem_dim_save) grid%next_free_elem = END_OF_LIST
   vert = grid%next_free_vert
   if (vert > vert_dim_save) grid%next_free_vert = END_OF_LIST
   edge = grid%next_free_edge
   if (edge > edge_dim_save) grid%next_free_edge = END_OF_LIST
   write(unit) grid%next_free_elem,grid%next_free_vert,grid%next_free_edge, &
               grid%partition, grid%nelem,grid%nelem_leaf,grid%nelem_leaf_own, &
               grid%nedge, grid%nedge_own, grid%nvert, grid%nvert_own, &
               grid%nlev,grid%dof,grid%dof_own
   write(unit) grid%oldsoln_exists
   grid%next_free_elem = elem
   grid%next_free_vert = vert
   grid%next_free_edge = edge
   write(unit) grid%head_level_elem
   write(unit) grid%head_level_vert
   if (me /= MASTER) then
      write(unit) grid%initial_neighbor
   endif
   write(unit) grid%boundbox_min,grid%boundbox_max
   if (associated(grid%eigenvalue)) then
      if (size(grid%eigenvalue) /= 0) then
         write(unit) grid%eigenvalue
         write(unit) grid%eigenprob_l2_resid
      endif
   endif
   call hash_table_store(grid%elem_hash,unit)
   call hash_table_store(grid%vert_hash,unit)
   call hash_table_store(grid%edge_hash,unit)
   if (vert_dim_save > 0) then
      call hash_print_key(grid%vertex(1:vert_dim_save)%gid,unit)
      write(unit) grid%vertex(1:vert_dim_save)%coord%x, &
                  grid%vertex(1:vert_dim_save)%coord%y
      write(unit) (grid%vertex_solution(j,:,:),j=1,vert_dim_save)
      if (grid%have_true) then
         write(unit) (grid%vertex_exact(j,:,:),j=1,vert_dim_save)
      endif
      write(unit) (grid%vertex_type(j,:),j=1,vert_dim_save)
      write(unit) grid%vertex(1:vert_dim_save)%bmark
      write(unit) grid%vertex(1:vert_dim_save)%assoc_elem
      write(unit) grid%vertex(1:vert_dim_save)%next
      write(unit) grid%vertex(1:vert_dim_save)%previous
      if (grid%oldsoln_exists) then
         write(unit) (grid%vertex_oldsoln(j,:,:),j=1,vert_dim_save)
      endif
   endif
   write(unit) VERTICES_PER_ELEMENT
   if (elem_dim_save > 0) then
      call hash_print_key(grid%element(1:elem_dim_save)%gid,unit)
      do j=1,VERTICES_PER_ELEMENT
         write(unit) grid%element(1:elem_dim_save)%vertex(j)
      end do
      call hash_print_key(grid%element(1:elem_dim_save)%mate,unit)
      write(unit) grid%element(1:elem_dim_save)%level
      write(unit) grid%element(1:elem_dim_save)%next
      write(unit) grid%element(1:elem_dim_save)%previous
      write(unit) grid%element(1:elem_dim_save)%iown
      do j=1,MAX_CHILD
         write(unit) grid%element(1:elem_dim_save)%order(j)
      end do
      write(unit) grid%element(1:elem_dim_save)%in
      write(unit) grid%element(1:elem_dim_save)%out
      write(unit) grid%element(1:elem_dim_save)%isleaf
      write(unit) grid%element(1:elem_dim_save)%degree
      do j=1,EDGES_PER_ELEMENT
         write(unit) grid%element(1:elem_dim_save)%edge(j)
      end do
      do j=1,elem_dim_save
         if (associated(grid%element(j)%solution)) then
            write(unit) size(grid%element(j)%solution,dim=1)
            write(unit) grid%element(j)%solution
         else
            write(unit) 0
         endif
         if (associated(grid%element(j)%exact)) then
            write(unit) size(grid%element(j)%exact,dim=1)
            write(unit) grid%element(j)%exact
         else
            write(unit) 0
         endif
         if (grid%oldsoln_exists) then
            if (associated(grid%element(j)%oldsoln)) then
               write(unit) size(grid%element(j)%oldsoln,dim=1), &
                           size(grid%element(j)%oldsoln,dim=2)
               write(unit) grid%element(j)%oldsoln
            else
               write(unit) 0,0
            endif
         endif
      end do
   endif
   if (edge_dim_save > 0) then
      call hash_print_key(grid%edge(1:edge_dim_save)%gid,unit)
      write(unit) grid%edge(1:edge_dim_save)%vertex(1), &
                  grid%edge(1:edge_dim_save)%vertex(2)
      write(unit) grid%edge(1:edge_dim_save)%degree
      write(unit) grid%edge(1:edge_dim_save)%bmark
      write(unit) grid%edge(1:edge_dim_save)%assoc_elem
      write(unit) grid%edge(1:edge_dim_save)%next
      do j=1,edge_dim_save
         write(unit) grid%edge_type(j,:)
      end do
      do j=1,edge_dim_save
         if (associated(grid%edge(j)%solution)) then
            write(unit) size(grid%edge(j)%solution,dim=1)
            write(unit) grid%edge(j)%solution
         else
            write(unit) 0
         endif
         if (associated(grid%edge(j)%exact)) then
            write(unit) size(grid%edge(j)%exact,dim=1)
            write(unit) grid%edge(j)%exact
         else
            write(unit) 0
         endif
         if (grid%oldsoln_exists) then
            if (associated(grid%edge(j)%oldsoln)) then
               write(unit) size(grid%edge(j)%oldsoln,dim=1), &
                           size(grid%edge(j)%oldsoln,dim=2)
               write(unit) grid%edge(j)%oldsoln
            else
               write(unit) 0,0
            endif
         endif
      end do
   endif

elseif (form == 'FORMATTED') then

   write(unit,*) '"PHAML data"'
   write(unit,*) 12 ! version number

! data from phaml_solution_type

   write(unit,*) phaml_solution%still_sequential
   write(unit,*) phaml_solution%pde_id
   write(unit,*) phaml_solution%system_size

! data from proc_info

   write(unit,*) nproc,me

! data from the grid

   grid => phaml_solution%grid
   if (me == MASTER) allocate(grid%initial_neighbor(1,1))
   elem_dim_save = 0
   vert_dim_save = 0
   edge_dim_save = 0
   do lev=1,grid%nlev
      elem = grid%head_level_elem(lev)
      do while (elem /= END_OF_LIST)
         elem_dim_save = max(elem,elem_dim_save)
         do vert = 1,VERTICES_PER_ELEMENT
            vert_dim_save = max(grid%element(elem)%vertex(vert),vert_dim_save)
         end do
         do edge = 1,EDGES_PER_ELEMENT
            edge_dim_save = max(grid%element(elem)%edge(edge),edge_dim_save)
         end do
         elem = grid%element(elem)%next
      end do
   end do
   write(unit,*) vert_dim_save,elem_dim_save,edge_dim_save, &
                size(grid%initial_neighbor,dim=2), &
                size(grid%head_level_elem)
   if (associated(grid%eigenvalue)) then
      write(unit,*) size(grid%head_level_vert), &
                    size(grid%eigenvalue),grid%nsoln
   else
      write(unit,*) size(grid%head_level_vert), &
                    -1,grid%nsoln
   endif
   write(unit,*) size(grid%vertex_solution,2),size(grid%vertex_solution,3)
   if (grid%have_true) then
      write(unit,*) 1
   else
      write(unit,*) 0
   endif
   write(unit,*) EDGES_PER_ELEMENT
   elem = grid%next_free_elem
   if (elem > elem_dim_save) grid%next_free_elem = END_OF_LIST
   vert = grid%next_free_vert
   if (vert > vert_dim_save) grid%next_free_vert = END_OF_LIST
   edge = grid%next_free_edge
   if (edge > edge_dim_save) grid%next_free_edge = END_OF_LIST
   write(unit,*) grid%next_free_elem,grid%next_free_vert,grid%next_free_edge, &
                grid%partition, grid%nelem
   write(unit,*) grid%nelem_leaf,grid%nelem_leaf_own,grid%nedge, &
                 grid%nedge_own,grid%nvert,grid%nvert_own,grid%nlev,grid%dof, &
                 grid%dof_own
   write(unit,*) grid%oldsoln_exists
   grid%next_free_elem = elem
   grid%next_free_vert = vert
   grid%next_free_edge = edge
   do sub=lbound(grid%head_level_elem,dim=1),ubound(grid%head_level_elem,dim=1)
      write(unit,*) grid%head_level_elem(sub)
   end do
   do sub=lbound(grid%head_level_vert,dim=1),ubound(grid%head_level_vert,dim=1)
      write(unit,*) grid%head_level_vert(sub)
   end do
   if (me /= MASTER) then
      do sub=lbound(grid%initial_neighbor,dim=2),ubound(grid%initial_neighbor,dim=2)
         write(unit,*) grid%initial_neighbor(:,sub)
      end do
   endif
   write(unit,*) grid%boundbox_min
   write(unit,*) grid%boundbox_max
   if (associated(grid%eigenvalue)) then
      if (size(grid%eigenvalue) /= 0) then
         write(unit,*) grid%eigenvalue
         write(unit,*) grid%eigenprob_l2_resid
      endif
   endif
   call hash_table_store(grid%elem_hash,unit)
   call hash_table_store(grid%vert_hash,unit)
   call hash_table_store(grid%edge_hash,unit)
   do sub=lbound(grid%vertex,dim=1),vert_dim_save
      call hash_print_key(grid%vertex(sub)%gid,unit)
      write(unit,*) grid%vertex(sub)%coord%x,grid%vertex(sub)%coord%y
      do i=1,size(grid%vertex_solution,2)
         do j=1,size(grid%vertex_solution,3)
            write(unit,*) grid%vertex_solution(sub,i,j)
         end do
      end do
      if (grid%have_true) then
         do i=1,size(grid%vertex_solution,2)
            do j=1,size(grid%vertex_solution,3)
               write(unit,*) grid%vertex_exact(sub,i,j)
            end do
         end do
      endif
      if (grid%oldsoln_exists) then
         do i=1,size(grid%vertex_solution,2)
            do j=1,size(grid%vertex_solution,3)
               write(unit,*) grid%vertex_oldsoln(sub,i,j)
            end do
         end do
      endif
      do j=1,grid%system_size
         write(unit,*) grid%vertex_type(sub,j)
      end do
      write(unit,*) grid%vertex(sub)%bmark
      write(unit,*) grid%vertex(sub)%assoc_elem
      write(unit,*) grid%vertex(sub)%next
      write(unit,*) grid%vertex(sub)%previous
   end do
   write(unit,*) VERTICES_PER_ELEMENT
   do sub=lbound(grid%element,dim=1),elem_dim_save
      call hash_print_key(grid%element(sub)%gid,unit)
      write(unit,*) grid%element(sub)%vertex
      call hash_print_key(grid%element(sub)%mate,unit)
      write(unit,*) grid%element(sub)%level
      write(unit,*) grid%element(sub)%next
      write(unit,*) grid%element(sub)%previous
      write(unit,"(L1)") grid%element(sub)%iown
      write(unit,*) grid%element(sub)%order
      write(unit,*) grid%element(sub)%in
      write(unit,*) grid%element(sub)%out
      write(unit,"(L1)") grid%element(sub)%isleaf
      write(unit,*) grid%element(sub)%degree
      write(unit,*) grid%element(sub)%edge
      if (associated(grid%element(sub)%solution)) then
         write(unit,*) size(grid%element(sub)%solution,dim=1)
         write(unit,*) grid%element(sub)%solution
      else
         write(unit,*) 0
      endif
      if (associated(grid%element(sub)%exact)) then
         write(unit,*) size(grid%element(sub)%exact,dim=1)
         write(unit,*) grid%element(sub)%exact
      else
         write(unit,*) 0
      endif
      if (grid%oldsoln_exists) then
         if (associated(grid%element(sub)%oldsoln)) then
            write(unit,*) size(grid%element(sub)%oldsoln,dim=1), &
                          size(grid%element(sub)%oldsoln,dim=2)
            write(unit,*) grid%element(sub)%oldsoln
         else
            write(unit,*) 0,0
         endif
      endif
   end do
   do sub=lbound(grid%edge,dim=1),edge_dim_save
      call hash_print_key(grid%edge(sub)%gid,unit)
      write(unit,*) grid%edge(sub)%vertex
      write(unit,*) grid%edge(sub)%degree
      write(unit,*) grid%edge(sub)%bmark
      write(unit,*) grid%edge(sub)%assoc_elem
      write(unit,*) grid%edge(sub)%next
      write(unit,*) grid%edge_type(sub,:)
      if (associated(grid%edge(sub)%solution)) then
         write(unit,*) size(grid%edge(sub)%solution,dim=1)
         write(unit,*) grid%edge(sub)%solution
      else
         write(unit,*) 0
      endif
      if (associated(grid%edge(sub)%exact)) then
         write(unit,*) size(grid%edge(sub)%exact,dim=1)
         write(unit,*) grid%edge(sub)%exact
      else
         write(unit,*) 0
      endif
      if (grid%oldsoln_exists) then
         if (associated(grid%edge(sub)%oldsoln)) then
            write(unit,*) size(grid%edge(sub)%oldsoln,dim=1), &
                          size(grid%edge(sub)%oldsoln,dim=2)
            write(unit,*) grid%edge(sub)%oldsoln
         else
            write(unit,*) 0,0
         endif
      endif
   end do

else

   ierr = USER_INPUT_ERROR
   call fatal("The file for phaml_store must be opened as either FORMATTED or UNFORMATTED", &
              "before calling phaml_store",procs=phaml_solution%procs)
endif

end subroutine phaml_store

!          -------------
subroutine phaml_restore(phaml_solution,unit,do_draw_grid,pause)
!          -------------

!----------------------------------------------------
! This routine restores information for phaml_solution from files created
! by subroutine phaml_store.  unit is the unit number to read from which
! should have been opened with phaml_popen
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type (phaml_solution_type), intent(inout), target :: phaml_solution
integer, intent(in) :: unit
logical, intent(in), optional :: do_draw_grid, pause
!----------------------------------------------------
! Local variables:

character(len=10) :: magic
integer :: send_int(4)
real(my_real) :: send_real(4)
integer :: i,j,k,proc,nproc,me,astat,dstat,version,size_vertex,sub, &
           size_element,size_initneigh,size_head_elem,size_head_vert, &
           size_eval,size_soln,size_edge,itemp,d1,d2, &
           saved_EDGES_PER_ELEMENT, saved_VERTICES_PER_ELEMENT, &
           size_soln2,size_soln3
type(grid_type), pointer :: grid
type (io_options) :: io_control
type (refine_options) :: ref_control
logical :: loc_draw_grid, loc_pause
character(len=11) :: form
!----------------------------------------------------
! Begin executable code

outunit = phaml_solution%outunit
errunit = phaml_solution%errunit

if (present(do_draw_grid)) then
   loc_draw_grid = do_draw_grid
else
   loc_draw_grid = .false.
endif
if (present(pause)) then
   loc_pause = pause
else
   loc_pause = .false.
endif
me = my_proc(phaml_solution%procs)
nproc = num_proc(phaml_solution%procs)

! if I am the master, tell the slaves to read the data

if (me == MASTER) then
   send_int(1) = 6 ! code for restoring solution
   send_int(2) = unit
   if (loc_draw_grid) then
      send_int(3) = 1
   else
      send_int(3) = 0
   endif
   if (loc_pause) then
      send_int(4) = 1
   else
      send_int(4) = 0
   endif
   do proc=1,nproc
      call phaml_send(phaml_solution%procs,proc,send_int,4,send_real,0,101)
   end do
endif

inquire(unit,form=form)

select case (form)
case ("UNDEFINED")
   ierr = USER_INPUT_ERROR
   call fatal("The file for phaml_restore must be opened before calling phaml_restore",procs=phaml_solution%procs)
   return

case ("FORMATTED")

   read(unit,*) magic
   if (magic /= "PHAML data") then
      ierr = USER_INPUT_ERROR
      call fatal("data file for phaml_restore does not appear to be a FORMATTED PHAML data file",procs=phaml_solution%procs)
      stop
   endif
   read(unit,*) version

   select case(version)

   case(10,11,12) ! version 10 PHAML data file
                  ! version 11 adds oldsoln
                  ! version 12 adds have_true, size_soln2, size_soln3

! data from phaml_solution_type

      read(unit,*) phaml_solution%still_sequential
      read(unit,*) phaml_solution%pde_id
      read(unit,*) phaml_solution%system_size

! data from proc_info

      read(unit,*) nproc,me
      if (nproc /= num_proc(phaml_solution%procs) .or. me /= my_proc(phaml_solution%procs)) then
         call fatal("num proc or my proc in data file does not match current configuration", &
                    intlist=(/nproc,num_proc(phaml_solution%procs),me, &
                    my_proc(phaml_solution%procs)/),procs=phaml_solution%procs)
         stop
      endif

! data from the grid

      grid => phaml_solution%grid
      grid%system_size = phaml_solution%system_size
      read(unit,*) size_vertex,size_element,size_edge,size_initneigh, &
                   size_head_elem
      read(unit,*) size_head_vert,size_eval,size_soln
      if (version >= 12) then
         read(unit,*) size_soln2, size_soln3
         read(unit,*) itemp
         grid%have_true = itemp==1
      else
         size_soln2 = grid%system_size
         size_soln3 = max(1,size_eval)
         grid%have_true = .true.
      endif
      read(unit,*) saved_EDGES_PER_ELEMENT
      if (saved_EDGES_PER_ELEMENT /= EDGES_PER_ELEMENT) then
         call fatal("phaml_restore: mismatch in number of edges per element", &
                    intlist=(/saved_EDGES_PER_ELEMENT,EDGES_PER_ELEMENT/),procs=phaml_solution%procs)
         stop
      endif
      allocate(grid%element(size_element),grid%vertex(size_vertex), &
               grid%initial_neighbor(EDGES_PER_ELEMENT,size_initneigh), &
               grid%head_level_elem(size_head_elem), &
               grid%head_level_vert(size_head_vert), grid%edge(size_edge), &
               stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("memory allocation failed in phaml_restore",procs=phaml_solution%procs)
         deallocate(grid%element,grid%vertex,grid%initial_neighbor, &
                    grid%head_level_elem,grid%head_level_vert, stat=dstat)
         deallocate(grid%edge, stat=dstat)
         return
      endif
      if (size_eval == -1) then
         nullify(grid%eigenvalue, grid%eigenprob_l2_resid, &
                 grid%eigenprob_variance)
      else
         allocate(grid%eigenvalue(size_eval), &
                  grid%eigenprob_variance(size_eval), &
                  grid%eigenprob_l2_resid(size_eval), stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("memory allocation failed in phaml_restore",procs=phaml_solution%procs)
            deallocate(grid%element,grid%vertex,grid%initial_neighbor, &
                       grid%head_level_elem,grid%head_level_vert, &
                       grid%eigenvalue, grid%eigenprob_l2_resid, &
                       grid%eigenprob_variance,stat=dstat)
            return
         endif
         grid%eigenvalue = 0.0_my_real
         grid%eigenprob_l2_resid = 0.0_my_real
         grid%eigenprob_variance = 0.0_my_real
      endif
      allocate(grid%vertex_type(size_vertex,phaml_solution%system_size), &
               grid%vertex_solution(size_vertex,size_soln2,size_soln3), &
               stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("memory allocation failed in phaml_restore",procs=phaml_solution%procs)
         deallocate(grid%element,grid%vertex,grid%initial_neighbor, &
                    grid%head_level_elem,grid%head_level_vert, stat=dstat)
         if (associated(grid%eigenvalue)) then
            deallocate(grid%eigenvalue, grid%eigenprob_l2_resid, &
                       grid%eigenprob_variance,stat=dstat)
         endif
         return
      endif
      if (grid%have_true) then
        allocate(grid%vertex_exact(size_vertex,size_soln2,size_soln3),stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("memory allocation failed in phaml_restore",procs=phaml_solution%procs)
            deallocate(grid%element,grid%vertex,grid%initial_neighbor, &
                       grid%head_level_elem,grid%head_level_vert, &
                       grid%vertex_solution, stat=dstat)
            if (associated(grid%eigenvalue)) then
               deallocate(grid%eigenvalue, grid%eigenprob_l2_resid, &
                          grid%eigenprob_variance,stat=dstat)
            endif
            return
         endif
      else
        nullify(grid%vertex_exact)
      endif
      if (grid%oldsoln_exists) then
        allocate(grid%vertex_oldsoln(size_vertex,size_soln2,size_soln3),stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("memory allocation failed in phaml_restore",procs=phaml_solution%procs)
            deallocate(grid%element,grid%vertex,grid%initial_neighbor, &
                       grid%head_level_elem,grid%head_level_vert, &
                       grid%vertex_solution, stat=dstat)
            if (grid%have_true) deallocate(grid%vertex_exact, stat=dstat)
            if (associated(grid%eigenvalue)) then
               deallocate(grid%eigenvalue, grid%eigenprob_l2_resid, &
                          grid%eigenprob_variance,stat=dstat)
            endif
            return
         endif
      else
        nullify(grid%vertex_oldsoln)
      endif
      grid%nsoln = size_soln
      grid%num_eval = max(0,size_eval)
      read(unit,*) grid%next_free_elem,grid%next_free_vert, &
                   grid%next_free_edge,grid%partition, grid%nelem
      read(unit,*) grid%nelem_leaf,grid%nelem_leaf_own,grid%nedge, &
                   grid%nedge_own,grid%nvert,grid%nvert_own,grid%nlev, &
                   grid%dof,grid%dof_own
      if (version >= 11) then
         read(unit,*) grid%oldsoln_exists
      else
         grid%oldsoln_exists = .false.
      endif
      do sub=lbound(grid%head_level_elem,dim=1),ubound(grid%head_level_elem,dim=1)
         read(unit,*) grid%head_level_elem(sub)
      end do
      do sub=lbound(grid%head_level_vert,dim=1),ubound(grid%head_level_vert,dim=1)
         read(unit,*) grid%head_level_vert(sub)
      end do
      if (me /= MASTER) then
         do sub=lbound(grid%initial_neighbor,dim=2),ubound(grid%initial_neighbor,dim=2)
            read(unit,*) grid%initial_neighbor(:,sub)
         end do
      endif
      read(unit,*) grid%boundbox_min
      read(unit,*) grid%boundbox_max
      if (associated(grid%eigenvalue)) then
         if (size_eval /= 0) then
            read(unit,*) grid%eigenvalue
            read(unit,*) grid%eigenprob_l2_resid
         endif
      endif
      call hash_table_restore(grid%elem_hash,unit)
      call hash_table_restore(grid%vert_hash,unit)
      call hash_table_restore(grid%edge_hash,unit)
      do sub=lbound(grid%vertex,dim=1),ubound(grid%vertex,dim=1)
         call hash_read_key(grid%vertex(sub:sub)%gid,unit)
         read(unit,*) grid%vertex(sub)%coord%x,grid%vertex(sub)%coord%y
         do i=1,size_soln2
            do j=1,size_soln3
               read(unit,*) grid%vertex_solution(sub,i,j)
            end do
         end do
         if (grid%have_true) then
            do i=1,size_soln2
               do j=1,size_soln3
                  read(unit,*) grid%vertex_exact(sub,i,j)
               end do
            end do
         endif
         if (grid%oldsoln_exists) then
            do i=1,size_soln2
               do j=1,size_soln3
                  read(unit,*) grid%vertex_oldsoln(sub,i,j)
               end do
            end do
         endif
         do j=1,phaml_solution%system_size
            read(unit,*) grid%vertex_type(sub,j)
         end do
         read(unit,*) grid%vertex(sub)%bmark
         read(unit,*) grid%vertex(sub)%assoc_elem
         read(unit,*) grid%vertex(sub)%next
         read(unit,*) grid%vertex(sub)%previous
      end do
      read(unit,*) saved_VERTICES_PER_ELEMENT
      if (saved_VERTICES_PER_ELEMENT /= VERTICES_PER_ELEMENT) then
         call fatal("phaml_restore: mismatch in number of vertices per element", &
                    intlist=(/saved_VERTICES_PER_ELEMENT,VERTICES_PER_ELEMENT/),procs=phaml_solution%procs)
         stop
      endif
      do sub=lbound(grid%element,dim=1),ubound(grid%element,dim=1)
         call hash_read_key(grid%element(sub:sub)%gid,unit)
         read(unit,*) grid%element(sub)%vertex
         call hash_read_key(grid%element(sub:sub)%mate,unit)
         read(unit,*) grid%element(sub)%level
         read(unit,*) grid%element(sub)%next
         read(unit,*) grid%element(sub)%previous
         read(unit,*) grid%element(sub)%iown
         grid%element(sub)%hrefined_unowned = .false.
         grid%element(sub)%prefined_unowned = .false.
         read(unit,*) grid%element(sub)%order
         read(unit,*) grid%element(sub)%in
         read(unit,*) grid%element(sub)%out
         read(unit,*) grid%element(sub)%isleaf
         read(unit,*) grid%element(sub)%degree
         read(unit,*) grid%element(sub)%edge
         read(unit,*) itemp
         if (itemp == 0) then
            nullify(grid%element(sub)%solution)
         else
            allocate(grid%element(sub)%solution(itemp,grid%system_size, &
                     max(1,grid%num_eval)),stat=astat)
            if (astat /= 0) then
               ierr = ALLOC_FAILED
               call fatal("memory allocation failed in phaml_restore",procs=phaml_solution%procs)
               do k=1,j-1
                  deallocate(grid%element(k)%solution,stat=dstat)
               end do
               deallocate(grid%vertex_type,grid%vertex_solution)
               if (grid%have_true) deallocate(grid%vertex_exact)
               if (grid%oldsoln_exists) deallocate(grid%vertex_oldsoln)
               deallocate(grid%element,grid%vertex,grid%initial_neighbor, &
                          grid%head_level_elem,grid%head_level_vert, stat=dstat)
               if (associated(grid%eigenvalue)) then
                  deallocate(grid%eigenvalue, grid%eigenprob_l2_resid, &
                             grid%eigenprob_variance,stat=dstat)
               endif
               return
            endif
            read(unit,*) grid%element(sub)%solution
         endif
         read(unit,*) itemp
         if (itemp == 0) then
            nullify(grid%element(sub)%exact)
         else
            allocate(grid%element(sub)%exact(itemp,grid%system_size, &
                     max(1,grid%num_eval)),stat=astat)
            if (astat /= 0) then
               ierr = ALLOC_FAILED
               call fatal("memory allocation failed in phaml_restore",procs=phaml_solution%procs)
               stop
            endif
            read(unit,*) grid%element(sub)%exact
         endif
         if (grid%oldsoln_exists) then
            read(unit,*) d1,d2
            if (d1 == 0) then
               nullify(grid%element(sub)%oldsoln)
            else
               allocate(grid%element(sub)%oldsoln(d1,d2,max(1,grid%num_eval)))
               read(unit,*) grid%element(sub)%oldsoln
            endif
         else
            nullify(grid%element(sub)%oldsoln)
         endif
      end do
      allocate(grid%edge_type(size_edge,grid%system_size),stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("memory allocation failed in phaml_restore",procs=phaml_solution%procs)
         stop
      endif
      do sub=lbound(grid%edge,dim=1),ubound(grid%edge,dim=1)
         call hash_read_key(grid%edge(sub:sub)%gid,unit)
         read(unit,*) grid%edge(sub)%vertex
         read(unit,*) grid%edge(sub)%degree
         read(unit,*) grid%edge(sub)%bmark
         read(unit,*) grid%edge(sub)%assoc_elem
         read(unit,*) grid%edge(sub)%next
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("memory allocation failed in phaml_restore",procs=phaml_solution%procs)
            stop
         endif
         read(unit,*) grid%edge_type(sub,:)
         read(unit,*) itemp
         if (itemp == 0) then
            nullify(grid%edge(sub)%solution)
         else
            allocate(grid%edge(sub)%solution(itemp,grid%system_size, &
                     max(1,grid%num_eval)),stat=astat)
            if (astat /= 0) then
               ierr = ALLOC_FAILED
               call fatal("memory allocation failed in phaml_restore",procs=phaml_solution%procs)
               do k=1,j-1
                  deallocate(grid%edge(k)%solution,stat=dstat)
               end do
               do k=1,size_element
                  deallocate(grid%element(k)%solution,stat=dstat)
               end do
               deallocate(grid%vertex_type,grid%vertex_solution)
               if (grid%have_true) deallocate(grid%vertex_exact)
               if (grid%oldsoln_exists) deallocate(grid%vertex_oldsoln)
               deallocate(grid%element,grid%vertex,grid%initial_neighbor, &
                          grid%head_level_elem,grid%head_level_vert, stat=dstat)
               if (associated(grid%eigenvalue)) then
                  deallocate(grid%eigenvalue, grid%eigenprob_l2_resid, &
                             grid%eigenprob_variance,stat=dstat)
               endif
               return
            endif
            read(unit,*) grid%edge(sub)%solution
         endif
         read(unit,*) itemp
         if (itemp == 0) then
            nullify(grid%edge(sub)%exact)
         else
            allocate(grid%edge(sub)%exact(itemp,grid%system_size, &
                     max(1,grid%num_eval)),stat=astat)
            if (astat /= 0) then
               ierr = ALLOC_FAILED
               call fatal("memory allocation failed in phaml_restore",procs=phaml_solution%procs)
               stop
            endif
            read(unit,*) grid%edge(sub)%exact
         endif
         if (grid%oldsoln_exists) then
            read(unit,*) d1,d2
            if (d1 == 0) then
               nullify(grid%edge(sub)%oldsoln)
            else
               allocate(grid%edge(sub)%oldsoln(d1,d2,max(1,grid%num_eval)))
               read(unit,*) grid%edge(sub)%oldsoln
            endif
         else
            nullify(grid%edge(sub)%oldsoln)
         endif
      end do

   case default ! version number not supported

      call fatal("unsupported version number for FORMATTED PHAML data file in phaml_restore", &
                 intlist=(/version/),procs=phaml_solution%procs)
      stop

   end select ! version

case ("UNFORMATTED")

   read(unit) magic
   if (magic /= "PHAML data") then
      call fatal("data file for phaml_restore does not appear to be an UNFORMATTED PHAML data file",procs=phaml_solution%procs)
      stop
   endif
   read(unit) version

   select case(version)

   case(10,11,12) ! version 10 PHAML data file
                  ! version 11 adds oldsoln
                  ! version 12 adds have_true, size_soln2, size_soln3

! data from phaml_solution_type

      read(unit) phaml_solution%still_sequential
      read(unit) phaml_solution%pde_id
      read(unit) phaml_solution%system_size

! data from proc_info

      read(unit) nproc,me
      if (nproc /= num_proc(phaml_solution%procs) .or. me /= my_proc(phaml_solution%procs)) then
         call fatal("num proc or my proc in data file does not match current configuration", &
                    intlist=(/nproc,num_proc(phaml_solution%procs),me, &
                    my_proc(phaml_solution%procs)/),procs=phaml_solution%procs)
         stop
      endif

! data from the grid

      grid => phaml_solution%grid
      grid%system_size = phaml_solution%system_size
      read(unit) size_vertex,size_element,size_edge,size_initneigh, &
                 size_head_elem,size_head_vert,size_eval,size_soln
      if (version >= 12) then
         read(unit) size_soln2, size_soln3
         read(unit) itemp
         grid%have_true = itemp==1
      else
         size_soln2 = grid%system_size
         size_soln3 = max(1,size_eval)
         grid%have_true = .true.
      endif
      read(unit) saved_EDGES_PER_ELEMENT
      if (saved_EDGES_PER_ELEMENT /= EDGES_PER_ELEMENT) then
         call fatal("phaml_restore: mismatch in number of faces per element", &
                    intlist=(/saved_EDGES_PER_ELEMENT,EDGES_PER_ELEMENT/),procs=phaml_solution%procs)
         stop
      endif
      allocate(grid%element(size_element),grid%vertex(size_vertex), &
               grid%initial_neighbor(EDGES_PER_ELEMENT,size_initneigh), &
               grid%head_level_elem(size_head_elem), &
               grid%head_level_vert(size_head_vert), grid%edge(size_edge), &
               stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("memory allocation failed in phaml_restore",procs=phaml_solution%procs)
         deallocate(grid%element,grid%vertex,grid%edge,grid%initial_neighbor, &
                    grid%head_level_elem,grid%head_level_vert,stat=dstat)
         deallocate(grid%edge, stat=dstat)
         return
      endif
      if (size_eval == -1) then
         nullify(grid%eigenvalue, grid%eigenprob_l2_resid, &
                 grid%eigenprob_variance)
      else
         allocate(grid%eigenvalue(size_eval), &
                  grid%eigenprob_variance(size_eval), &
                  grid%eigenprob_l2_resid(size_eval), stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("memory allocation failed in phaml_restore",procs=phaml_solution%procs)
            deallocate(grid%element,grid%vertex,grid%initial_neighbor, &
                       grid%head_level_elem,grid%head_level_vert, &
                       grid%eigenvalue, grid%eigenprob_l2_resid, &
                       grid%eigenprob_variance,stat=dstat)
            return
         endif
         grid%eigenvalue = 0.0_my_real
         grid%eigenprob_l2_resid = 0.0_my_real
         grid%eigenprob_variance = 0.0_my_real
      endif
      allocate(grid%vertex_type(size_vertex,phaml_solution%system_size), &
               grid%vertex_solution(size_vertex,size_soln2,size_soln3), &
               stat=astat)
      if (astat /= 0) then
         ierr = ALLOC_FAILED
         call fatal("memory allocation failed in phaml_restore",procs=phaml_solution%procs)
         deallocate(grid%element,grid%vertex,grid%initial_neighbor, &
                    grid%head_level_elem,grid%head_level_vert, stat=dstat)
         if (associated(grid%eigenvalue)) then
            deallocate(grid%eigenvalue, grid%eigenprob_l2_resid, &
                       grid%eigenprob_variance,stat=dstat)
         endif
         return
      endif
      if (grid%have_true) then
         allocate(grid%vertex_exact(size_vertex,size_soln2,size_soln3), &
                  stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("memory allocation failed in phaml_restore",procs=phaml_solution%procs)
            deallocate(grid%element,grid%vertex,grid%initial_neighbor, &
                       grid%head_level_elem,grid%head_level_vert, &
                       grid%vertex_solution,stat=dstat)
            if (associated(grid%eigenvalue)) then
               deallocate(grid%eigenvalue, grid%eigenprob_l2_resid, &
                          grid%eigenprob_variance,stat=dstat)
            endif
            return
         endif
      else
         nullify(grid%vertex_exact)
      endif
      if (grid%oldsoln_exists) then
         allocate(grid%vertex_oldsoln(size_vertex,size_soln2,size_soln3), &
                  stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("memory allocation failed in phaml_restore",procs=phaml_solution%procs)
            deallocate(grid%element,grid%vertex,grid%initial_neighbor, &
                       grid%head_level_elem,grid%head_level_vert, &
                       grid%vertex_solution,stat=dstat)
            if (grid%have_true) deallocate(grid%vertex_exact,stat=dstat)
            if (associated(grid%eigenvalue)) then
               deallocate(grid%eigenvalue, grid%eigenprob_l2_resid, &
                          grid%eigenprob_variance,stat=dstat)
            endif
            return
         endif
      else
         nullify(grid%vertex_oldsoln)
      endif
      grid%nsoln = size_soln
      grid%num_eval = max(0,size_eval)
      read(unit) grid%next_free_elem,grid%next_free_vert, &
                 grid%next_free_edge,grid%partition, &
                 grid%nelem,grid%nelem_leaf,grid%nelem_leaf_own,grid%nedge, &
                 grid%nedge_own,grid%nvert,grid%nvert_own,grid%nlev, &
                 grid%dof,grid%dof_own
      if (version >= 11) then
         read(unit) grid%oldsoln_exists
      else
         grid%oldsoln_exists = .false.
      endif
      read(unit) grid%head_level_elem
      read(unit) grid%head_level_vert
      if (me /= MASTER) then
         read(unit) grid%initial_neighbor
      endif
      read(unit) grid%boundbox_min,grid%boundbox_max
      if (associated(grid%eigenvalue)) then
         if (size_eval /= 0) then
            read(unit) grid%eigenvalue
            read(unit) grid%eigenprob_l2_resid
         endif
      endif
      call hash_table_restore(grid%elem_hash,unit)
      call hash_table_restore(grid%vert_hash,unit)
      call hash_table_restore(grid%edge_hash,unit)
      if (size_vertex > 0) then
         call hash_read_key(grid%vertex%gid,unit)
         read(unit) grid%vertex%coord%x,grid%vertex%coord%y
         read(unit) (grid%vertex_solution(j,:,:),j=1,size(grid%vertex))
         if (grid%have_true) then
            read(unit) (grid%vertex_exact(j,:,:),j=1,size(grid%vertex))
         endif
         read(unit) (grid%vertex_type(j,:),j=1,size(grid%vertex))
         read(unit) grid%vertex%bmark
         read(unit) grid%vertex%assoc_elem
         read(unit) grid%vertex%next
         read(unit) grid%vertex%previous
         if (grid%oldsoln_exists) then
            read(unit) (grid%vertex_oldsoln(j,:,:),j=1,size(grid%vertex))
         endif
      endif
      read(unit) saved_VERTICES_PER_ELEMENT
      if (saved_VERTICES_PER_ELEMENT /= VERTICES_PER_ELEMENT) then
         call fatal("phaml_restore: mismatch in number of vertices per element", &
                    intlist=(/saved_VERTICES_PER_ELEMENT,VERTICES_PER_ELEMENT/),procs=phaml_solution%procs)
         stop
      endif
      if (size_element > 0) then
         call hash_read_key(grid%element%gid,unit)
         do j=1,VERTICES_PER_ELEMENT
            read(unit) grid%element%vertex(j)
         end do
         call hash_read_key(grid%element%mate,unit)
         read(unit) grid%element%level
         read(unit) grid%element%next
         read(unit) grid%element%previous
         read(unit) grid%element%iown
         grid%element%hrefined_unowned = .false.
         grid%element%prefined_unowned = .false.
         do j=1,MAX_CHILD
            read(unit) grid%element%order(j)
         end do
         read(unit) grid%element%in
         read(unit) grid%element%out
         read(unit) grid%element%isleaf
         read(unit) grid%element%degree
         do j=1,EDGES_PER_ELEMENT
            read(unit) grid%element%edge(j)
         end do
         do j=1,size(grid%element)
            read(unit) itemp
            if (itemp == 0) then
               nullify(grid%element(j)%solution)
            else
               allocate(grid%element(j)%solution(itemp,grid%system_size, &
                        max(1,grid%num_eval)), stat=astat)
               if (astat /= 0) then
                  ierr = ALLOC_FAILED
                  call fatal("memory allocation failed in phaml_restore",procs=phaml_solution%procs)
                  do k=1,j-1
                     deallocate(grid%element(k)%solution,stat=dstat)
                  end do
                  deallocate(grid%vertex_type,grid%vertex_solution,stat=dstat)
                  if (grid%have_true) deallocate(grid%vertex_exact,stat=dstat)
                  if (grid%oldsoln_exists) deallocate(grid%vertex_oldsoln,stat=dstat)
                  deallocate(grid%element,grid%vertex,grid%edge, &
                             grid%initial_neighbor,grid%head_level_elem, &
                             grid%head_level_vert, stat=dstat)
                  if (associated(grid%eigenvalue)) then
                     deallocate(grid%eigenvalue, grid%eigenprob_l2_resid, &
                                grid%eigenprob_variance,stat=dstat)
                  endif
                  return
               endif
               read(unit) grid%element(j)%solution
            endif
            read(unit) itemp
            if (itemp == 0) then
               nullify(grid%element(j)%exact)
            else
               allocate(grid%element(j)%exact(itemp,grid%system_size, &
                        max(1,grid%num_eval)), stat=astat)
               if (astat /= 0) then
                  ierr = ALLOC_FAILED
                  call fatal("memory allocation failed in phaml_restore",procs=phaml_solution%procs)
                  stop
               endif
               read(unit) grid%element(j)%exact
            endif
            if (grid%oldsoln_exists) then
               read(unit) d1,d2
               if (d1 == 0) then
                  nullify(grid%element(j)%oldsoln)
               else
                  allocate(grid%element(j)%oldsoln(d1,d2,max(1,grid%num_eval)))
                  read(unit) grid%element(j)%oldsoln
               endif
            else
               nullify(grid%element(j)%oldsoln)
            endif
         end do
      endif
      if (size_edge > 0) then
         call hash_read_key(grid%edge%gid,unit)
         read(unit) grid%edge%vertex(1),grid%edge%vertex(2)
         read(unit) grid%edge%degree
         read(unit) grid%edge%bmark
         read(unit) grid%edge%assoc_elem
         read(unit) grid%edge%next
         allocate(grid%edge_type(size_edge,grid%system_size),stat=astat)
         if (astat /= 0) then
            ierr = ALLOC_FAILED
            call fatal("memory allocation failed in phaml_restore",procs=phaml_solution%procs)
            stop
         endif
         do j=1,size(grid%edge)
            read(unit) grid%edge_type(j,:)
         end do
         do j=1,size(grid%edge)
            read(unit) itemp
            if (itemp == 0) then
               nullify(grid%edge(j)%solution)
            else
               allocate(grid%edge(j)%solution(itemp,grid%system_size, &
                        max(1,grid%num_eval)), stat=astat)
               if (astat /= 0) then
                  ierr = ALLOC_FAILED
                  call fatal("memory allocation failed in phaml_restore",procs=phaml_solution%procs)
                  do k=1,j-1
                     deallocate(grid%edge(k)%solution,stat=dstat)
                  end do
                  do k=1,size_element
                     deallocate(grid%element(k)%solution,stat=dstat)
                  end do
                  deallocate(grid%vertex_type,grid%vertex_solution,stat=dstat)
                  if (grid%have_true) deallocate(grid%vertex_exact,stat=dstat)
                  if (grid%oldsoln_exists) deallocate(grid%vertex_oldsoln,stat=dstat)
                  deallocate(grid%element,grid%vertex,grid%edge, &
                             grid%initial_neighbor,grid%head_level_elem, &
                             grid%head_level_vert, stat=dstat)
                  if (associated(grid%eigenvalue)) then
                     deallocate(grid%eigenvalue, grid%eigenprob_l2_resid, &
                                grid%eigenprob_variance,stat=dstat)
                  endif
                  return
               endif
               read(unit) grid%edge(j)%solution
            endif
            read(unit) itemp
            if (itemp == 0) then
               nullify(grid%edge(j)%exact)
            else
               allocate(grid%edge(j)%exact(itemp,grid%system_size, &
                        max(1,grid%num_eval)), stat=astat)
               if (astat /= 0) then
                  ierr = ALLOC_FAILED
                  call fatal("memory allocation failed in phaml_restore",procs=phaml_solution%procs)
               endif
               read(unit) grid%edge(j)%exact
            endif
            if (grid%oldsoln_exists) then
               read(unit) d1,d2
               if (d1 == 0) then
                  nullify(grid%edge(j)%oldsoln)
               else
                  allocate(grid%edge(j)%oldsoln(d1,d2,max(1,grid%num_eval)))
                  read(unit) grid%edge(j)%oldsoln
               endif
            else
               nullify(grid%edge(j)%oldsoln)
            endif
         end do
      endif

   case default ! version number not supported

      call fatal("unsupported version number for UNFORMATTED PHAML data file in phaml_restore", &
                 intlist=(/version/),procs=phaml_solution%procs)
      stop

   end select ! version
end select ! form

allocate(grid%element_errind(size(grid%element),max(1,grid%num_eval)), &
         grid%errest_energy(max(1,grid%num_eval)), &
         grid%errest_Linf(grid%nsoln), grid%errest_L2(grid%nsoln), &
         grid%errest_eigenvalue(max(1,grid%num_eval)), stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("memory allocation failed in phaml_restore",procs=phaml_solution%procs)
   stop
endif

j = grid%next_free_elem
do while (j /= END_OF_LIST)
   if (grid%element(j)%next > size(grid%element)) then
      grid%element(j)%next = END_OF_LIST
      exit
   endif
   j = grid%element(j)%next
end do

j = grid%next_free_edge
do while (j /= END_OF_LIST)
   if (grid%edge(j)%next > size(grid%edge)) then
      grid%edge(j)%next = END_OF_LIST
      exit
   endif
   j = grid%edge(j)%next
end do

j = grid%next_free_vert
do while (j /= END_OF_LIST)
   if (grid%vertex(j)%next > size(grid%vertex)) then
      grid%vertex(j)%next = END_OF_LIST
      exit
   endif
   j = grid%vertex(j)%next
end do

grid%errind_up2date = .false.
grid_changed = .true.

! draw grid

if (loc_draw_grid) then
   io_control = io_options(NEVER,NO_ONE,NEVER,NO_ONE,NEVER,NO_ONE,NEVER,NEVER, &
                           NEVER,NO_ONE,PHASES,loc_pause)
   ref_control = refine_options(0.0_my_real,0.0_my_real,0.0_my_real, &
                               0.0_my_real,0.0_my_real,0.0_my_real,0.0_my_real,&
                                0.0_my_real,0.0_my_real,0.0_my_real, &
                                0.0_my_real,0.0_my_real,0.0_my_real, &
                                EXPLICIT_ERRIND, &
                                HP_ADAPTIVE,0,0,HP_BIGGER_ERRIND,0,0,0,0,0,0, &
                                0,0,0,0,0,0,0,.false.)
   call draw_grid(phaml_solution%grid,phaml_solution%procs, &
                  io_control,ref_control,phaml_solution%i_draw_grid, &
                  phaml_solution%master_draws_grid, &
                  phaml_solution%still_sequential, (/PHASES/), &
                  RTK, phaml_solution%lb)
endif

end subroutine phaml_restore

!          ------------------
subroutine phaml_store_matrix(phaml_solution,stiffness_unit,rhs_unit,mass_unit,&
                              inc_quad_order)
!          ------------------

!----------------------------------------------------
! This routine writes the stiffness matrix, mass matrix and/or right hand
! side to files in Matrix Market format.  The matrix is the global matrix
! collected by the master from the slaves.  The rhs is an Nx1 matrix.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(phaml_solution_type), intent(in) :: phaml_solution
integer, optional, intent(in) :: stiffness_unit, rhs_unit, mass_unit, &
                                 inc_quad_order
!----------------------------------------------------
! Local variables:

integer :: proc, loc_inc
!----------------------------------------------------
! Begin executable code

! check for invalid input

if (phaml_solution%eq_type == ELLIPTIC .and. present(mass_unit)) then
   call warning("phaml_store_matrix: Cannot write mass matrix without an eigenvalue problem.", &
                "No matrices stored.")
   return
endif

if (phaml_solution%eq_type == EIGENVALUE .and. present(rhs_unit)) then
   call warning("phaml_store_matrix: Cannot write right hand side for an eigenvalue problem.", &
                "No matrices stored.")
   return
endif

if (present(inc_quad_order)) then
   loc_inc = inc_quad_order
else
   loc_inc = 0
endif

! tell the slaves to store matrix

do proc=1,num_proc(phaml_solution%procs)
   call phaml_send(phaml_solution%procs,proc,(/17,loc_inc/),2,(/1.0_my_real/), &
                   0,101)
end do

! store matrix

call store_matrix(phaml_solution%grid,phaml_solution%procs, &
                  phaml_solution%still_sequential,phaml_solution%system_size, &
                  phaml_solution%eq_type,loc_inc,stiffness_unit,rhs_unit, &
                  mass_unit)

end subroutine phaml_store_matrix

!          --------------
subroutine phaml_compress(phaml_solution)
!          --------------

!----------------------------------------------------
! This routine compresses the element, edge and vertex data by moving all the
! unused elements and vertices after the used ones
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(phaml_solution_type), intent(inout), target :: phaml_solution
!----------------------------------------------------
! Local variables:

integer :: astat, elem, num_init_elem, big_elem, lev, edge, num_init_edge, &
           big_edge, vert, num_init_vert, big_vert, count, i, send_int(1), proc
type(grid_type), pointer :: grid
integer, allocatable :: renum_elem(:), renum_edge(:), renum_vert(:)
logical, allocatable :: visited(:)
!----------------------------------------------------
! Begin executable code

! if I am the master, tell the slaves to compress; master does not have a grid

if (my_proc(phaml_solution%procs) == MASTER) then
   send_int(1) = 15 ! code for compress
   do proc=1,num_proc(phaml_solution%procs)
      call phaml_send(phaml_solution%procs,proc,send_int,1,(/1.0_my_real/),0,101)
   end do
   return
endif

grid => phaml_solution%grid

! allocate space to store the renumberings

allocate(renum_elem(size(grid%element)), renum_edge(size(grid%edge)), &
         renum_vert(size(grid%vertex)), visited(size(grid%edge)), stat=astat)
if (astat /= 0) then
   ierr = ALLOC_FAILED
   call fatal("allocation failed in phaml_compress", &
              procs=phaml_solution%procs)
   return
endif

! count the number of elements in the first level, i.e. initial grid

num_init_elem = 0
big_elem = 0
lev = 1
elem = grid%head_level_elem(lev)
do while (elem /= END_OF_LIST)
   num_init_elem = num_init_elem + 1
   big_elem = max(big_elem,elem)
   elem = grid%element(elem)%next
end do

! make sure the elements are dense at the beginning of grid%element

if (big_elem /= num_init_elem) then
   call warning("Cannot compress data structures unless the initial grid uses lowest numbers.", &
      "Number of initial elements and biggest index are ", &
      intlist=(/num_init_elem,big_elem/))
   return
endif

! find the biggest index for any element

do lev=2,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      big_elem = max(big_elem,elem)
      elem = grid%element(elem)%next
   end do
end do

! go through all the elements and mark them as in-use in the renum array

renum_elem = 0
do lev = 1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      renum_elem(elem) = 1
      elem = grid%element(elem)%next
   end do
end do

! pass through the renum array to determine the renumbering

count = 0
do i=1,big_elem
   if (renum_elem(i) == 1) then
      count = count + 1
      renum_elem(i) = count
   endif
end do

! make sure the count came out right

if (count /= grid%nelem) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("count of elements does not agree with nelem in phaml_compress")
   return
endif

! count the number of edges in the first level, i.e. initial grid

num_init_edge = 0
big_edge = 0
lev = 1
visited = .false.
elem = grid%head_level_elem(lev)
do while (elem /= END_OF_LIST)
   do i=1,EDGES_PER_ELEMENT
      edge = grid%element(elem)%edge(i)
      if (.not. visited(edge)) then
         num_init_edge = num_init_edge + 1
         big_edge = max(big_edge,edge)
         visited(edge) = .true.
      endif
   end do
   elem = grid%element(elem)%next
end do

! make sure the edges are dense at the beginning of grid%edge

if (big_edge /= num_init_edge) then
   call warning("Cannot compress data structures unless the initial grid uses lowest numbers.", &
      "Number of initial edges and biggest index are ", &
      intlist=(/num_init_edge,big_edge/))
   return
endif

! find the biggest index for any edge

do lev=2,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      do i=1,3
         big_edge = max(big_edge,grid%element(elem)%edge(i))
      end do
      elem = grid%element(elem)%next
   end do
end do

! go through all the edges and mark them as in-use in the renum array

renum_edge = 0
do lev = 1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      do i=1,EDGES_PER_ELEMENT
         edge = grid%element(elem)%edge(i)
         renum_edge(edge) = 1
      end do
      elem = grid%element(elem)%next
   end do
end do

! pass through the renum array to determine the renumbering

count = 0
do i=1,big_edge
   if (renum_edge(i) == 1) then
      count = count + 1
      renum_edge(i) = count
   endif
end do

! make sure the count came out right

if (count /= grid%nedge) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("count of edges does not agree with nedge in phaml_compress")
   return
endif

! count the number of vertices in the first level, i.e. initial grid

num_init_vert = 0
big_vert = 0
lev = 1
vert = grid%head_level_vert(lev)
do while (vert /= END_OF_LIST)
   num_init_vert = num_init_vert + 1
   big_vert = max(big_vert,vert)
   vert = grid%vertex(vert)%next
end do

! make sure the vertices are dense at the beginning of grid%vertex

if (big_vert /= num_init_vert) then
   call warning("Cannot compress data structures unless the initial grid uses lowest numbers.", &
      "Number of initial vertices and biggest index are ", &
      intlist=(/num_init_vert,big_vert/))
   return
endif

do lev=2,grid%nlev
   vert = grid%head_level_vert(lev)
   do while (vert /= END_OF_LIST)
      big_vert = max(big_vert,vert)
      vert = grid%vertex(vert)%next
   end do
end do

! go through all the vertices and mark them as in-use in the renum array

renum_vert = 0
do lev = 1,grid%nlev
   vert = grid%head_level_vert(lev)
   do while (vert /= END_OF_LIST)
      renum_vert(vert) = 1
      vert = grid%vertex(vert)%next
   end do
end do

! pass through the renum array to determine the renumbering

count = 0
do i=1,big_vert
   if (renum_vert(i) == 1) then
      count = count + 1
      renum_vert(i) = count
   endif
end do

! make sure the count came out right

if (count /= grid%nvert) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("count of vertices does not agree with nvert in phaml_compress")
   return
endif

! go through the elements, moving them to the new position and correcting
! local indices.

elem = 0
do i=1,big_elem
   if (renum_elem(i) == 0) cycle
   elem = elem + 1
   if (elem /= i) then
      if (associated(grid%element(elem)%solution)) then
         deallocate(grid%element(elem)%solution)
      endif
      if (associated(grid%element(elem)%exact)) then
         deallocate(grid%element(elem)%exact)
      endif
      if (associated(grid%element(elem)%oldsoln)) then
         deallocate(grid%element(elem)%oldsoln)
      endif
      grid%element(elem) = grid%element(i)
      grid%element_errind(elem,:) = grid%element_errind(i,:)
      nullify(grid%element(i)%solution, grid%element(i)%exact, &
              grid%element(i)%oldsoln)
      call hash_remove(grid%element(elem)%gid,grid%elem_hash)
      call hash_insert(grid%element(elem)%gid,elem,grid%elem_hash)
   endif
   do vert=1,VERTICES_PER_ELEMENT
      grid%element(elem)%vertex(vert) = renum_vert(grid%element(elem)%vertex(vert))
   end do
   do edge=1,EDGES_PER_ELEMENT
      grid%element(elem)%edge(edge) = renum_edge(grid%element(elem)%edge(edge))
   end do
   grid%element(elem)%in = renum_vert(grid%element(elem)%in)
   grid%element(elem)%out = renum_vert(grid%element(elem)%out)
   if (grid%element(elem)%next /= END_OF_LIST) then
      grid%element(elem)%next = renum_elem(grid%element(elem)%next)
   endif
   if (grid%element(elem)%previous /= END_OF_LIST) then
      grid%element(elem)%previous = renum_elem(grid%element(elem)%previous)
   endif
end do

! set the next/previous for unused elements, and head element for each level

do i=grid%nelem+1,size(grid%element)-1
   grid%element(i)%next = i+1
end do
grid%element(size(grid%element))%next = END_OF_LIST
do i=grid%nelem+2,size(grid%element)
   grid%element(i)%previous = i-1
end do
grid%element(grid%nelem+1)%previous = END_OF_LIST
grid%next_free_elem = grid%nelem+1
do lev=2,grid%nlev
   if (grid%head_level_elem(lev) /= END_OF_LIST) then
      grid%head_level_elem(lev) = renum_elem(grid%head_level_elem(lev))
   endif
end do

! go through the edges, moving them to the new position and correcting
! local indices.

edge = 0
do i=1,big_edge
   if (renum_edge(i) == 0) cycle
   edge = edge + 1
   if (edge /= i) then
      if (associated(grid%edge(edge)%solution)) then
         deallocate(grid%edge(edge)%solution)
      endif
      if (associated(grid%edge(edge)%exact)) then
         deallocate(grid%edge(edge)%exact)
      endif
      if (associated(grid%edge(edge)%oldsoln)) then
         deallocate(grid%edge(edge)%oldsoln)
      endif
      grid%edge(edge) = grid%edge(i)
      grid%edge_type(edge,:) = grid%edge_type(i,:)
      nullify(grid%edge(i)%solution, grid%edge(i)%exact, grid%edge(i)%oldsoln)
      call hash_remove(grid%edge(edge)%gid,grid%edge_hash)
      call hash_insert(grid%edge(edge)%gid,edge,grid%edge_hash)
   endif
   do vert=1,2
      grid%edge(edge)%vertex(vert) = renum_vert(grid%edge(edge)%vertex(vert))
   end do
   grid%edge(edge)%assoc_elem = renum_elem(grid%edge(edge)%assoc_elem)
   if (any(grid%edge_type(edge,:) == PERIODIC_SLAVE)) then
      grid%edge(edge)%next = renum_edge(grid%edge(edge)%next)
   endif
end do

grid%next_free_edge = grid%nedge+1

! go through the vertices, moving them to the new position and correcting
! local indices.

vert = 0
do i=1,big_vert
   if (renum_vert(i) == 0) cycle
   vert = vert + 1
   if (vert /= i) then
      grid%vertex(vert) = grid%vertex(i)
      grid%vertex_type(vert,:) = grid%vertex_type(i,:)
      grid%vertex_solution(vert,:,:) = grid%vertex_solution(i,:,:)
      grid%vertex_exact(vert,:,:) = grid%vertex_exact(i,:,:)
      grid%vertex_oldsoln(vert,:,:) = grid%vertex_oldsoln(i,:,:)
      call hash_remove(grid%vertex(vert)%gid,grid%vert_hash)
      call hash_insert(grid%vertex(vert)%gid,vert,grid%vert_hash)
   endif
   grid%vertex(vert)%assoc_elem = renum_elem(grid%vertex(vert)%assoc_elem)
   if (grid%vertex(vert)%next /= END_OF_LIST) then
      grid%vertex(vert)%next = renum_vert(grid%vertex(vert)%next)
   endif
   if (grid%vertex(vert)%previous /= END_OF_LIST) then
      grid%vertex(vert)%previous = renum_vert(grid%vertex(vert)%previous)
   endif
end do

! set the next/previous for unused vertices, and head vertex for each level

do i=grid%nvert+1,size(grid%vertex)-1
   grid%vertex(i)%next = i+1
end do
grid%vertex(size(grid%vertex))%next = END_OF_LIST
do i=grid%nvert+2,size(grid%vertex)
   grid%vertex(i)%previous = i-1
end do
grid%vertex(grid%nvert+1)%previous = END_OF_LIST
grid%next_free_vert = grid%nvert+1
do lev=2,grid%nlev
   if (grid%head_level_vert(lev) /= END_OF_LIST) then
      grid%head_level_vert(lev) = renum_vert(grid%head_level_vert(lev))
   endif
end do

deallocate(renum_elem,renum_edge,renum_vert,visited,stat=astat)

end subroutine phaml_compress

!          -----------
subroutine phaml_popen(phaml_solution,unit,file,form)
!          -----------

!----------------------------------------------------
! This routine opens unit number unit on all processors in phaml_solution.
! See subroutine make_filename for an explanation of the actual filenames.
! The number of processors is limited to 9999.
! form indicates whether it should be opened as FORMATTED or UNFORMATTED;
! default is FORMATTED
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type (phaml_solution_type), intent(in) :: phaml_solution
integer, intent(in) :: unit
character(len=*), intent(in) :: file
character(len=*), intent(in), optional :: form
!----------------------------------------------------
! Local variables:

integer, parameter :: NAMELEN = 128 ! RESTRICTION 128 character base name
character(len=NAMELEN+4) :: full_filename ! RESTRICTION 9999 processors
integer :: send_int(3)
real(my_real) :: send_real(1)
integer :: proc,iostat
character(len=11) :: loc_form
!----------------------------------------------------
! Begin executable code

outunit = phaml_solution%outunit
errunit = phaml_solution%errunit

if (present(form)) then
   if (form == "FORMATTED" .or. form == "UNFORMATTED") then
      loc_form = form
   else
      call fatal("phaml_popen: form must be FORMATTED or UNFORMATTED")
      stop
   endif
else
   loc_form = "FORMATTED"
endif

if (unit < 0) then
   call fatal("phaml_popen: unit must be a nonnegative integer.  It is ", &
              intlist=(/unit/))
   stop
endif

! if I am the master, tell the slaves to open the unit

if (my_proc(phaml_solution%procs) == MASTER) then
   send_int(1) = 7 ! code for saving solution
   send_int(2) = unit
   if (loc_form == "FORMATTED") then
      send_int(3) = 0
   else
      send_int(3) = 1
   endif
   do proc=1,num_proc(phaml_solution%procs)
      call phaml_send(phaml_solution%procs,proc,send_int,3,send_real,0,101)
   end do
endif

! create the filename

call make_filename(file,full_filename,phaml_solution)

! open it

open(unit=unit,file=trim(full_filename),form=loc_form,iostat=iostat)
if (iostat /= 0) then
   call fatal("failed to open i/o unit",intlist=(/unit, iostat/),procs=phaml_solution%procs)
endif

end subroutine phaml_popen

!          ------------
subroutine phaml_pclose(phaml_solution,unit)
!          ------------

!----------------------------------------------------
! This routine closes unit number unit on all processors in phaml_solution.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type (phaml_solution_type), intent(in) :: phaml_solution
integer, intent(in) :: unit
!----------------------------------------------------
! Local variables:

integer :: send_int(2)
real(my_real) :: send_real(1)
integer :: proc,iostat
!----------------------------------------------------
! Begin executable code

outunit = phaml_solution%outunit
errunit = phaml_solution%errunit

! if I am the master, tell the slaves to close the unit

if (my_proc(phaml_solution%procs) == MASTER) then
   send_int(1) = 8 ! code for saving solution
   send_int(2) = unit
   do proc=1,num_proc(phaml_solution%procs)
      call phaml_send(phaml_solution%procs,proc,send_int,2,send_real,0,101)
   end do
endif

! close it

close(unit=unit,iostat=iostat)
if (iostat /= 0) then
   call warning("phaml_pclose: failed to close i/o unit",intlist=(/unit/))
endif

end subroutine phaml_pclose

!          -------------
subroutine make_filename(inname,outname,phaml_solution)
!          -------------

!----------------------------------------------------
! This routine creates a filename for this processor based on inname.
! If inname is of the form root.suffix then outname is
! rootXXXX.suffix for processor number XXXX, where the number of digits
! in XXXX is the minimum needed for the number of processors in
! phaml_solution.  The master processor is number 0.  If there is
! no "." in inname, then there is no suffix and XXXX is appended to
! inname.  The filename is limited to 128 characters and the number
! of processors is limited to 9999.  Only the master process uses
! inname; other processors receive it from the master.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

character(len=*), intent(in) :: inname
character(len=*), intent(out) :: outname
type (phaml_solution_type), intent(in) :: phaml_solution
!----------------------------------------------------
! Local variables:

integer, parameter :: NAMELEN = 128 ! RESTRICTION 128 character base name
character(len=NAMELEN) :: filename, rootname, suffixname
character(len=NAMELEN+4) :: ext_filename ! RESTRICTION 9999 processors
integer :: send_int(NAMELEN)
real(my_real) :: send_real(1)
integer, pointer :: recv_int(:)
real(my_real), pointer :: recv_real(:)
integer :: i,proc,nproc,me,start,sep,ni,nr,astat
!----------------------------------------------------
! Begin executable code

me = my_proc(phaml_solution%procs)
nproc = num_proc(phaml_solution%procs)

! truncate inname to NAMELEN characters, or pad with blanks

if (me == MASTER .or. PARALLEL==SEQUENTIAL) then
   if (len(inname) > NAMELEN) then
      call warning("filename is too long.  Truncating to",intlist=(/NAMELEN/))
   endif
   filename = inname
endif

! send the base name from the master to the slaves

if (me == MASTER .or. PARALLEL==SEQUENTIAL) then
   do i=1,NAMELEN
      send_int(i) = ichar(filename(i:i))
   end do
   do proc=1,nproc
      call phaml_send(phaml_solution%procs,proc,send_int,NAMELEN,send_real,0, &
                      320)
   end do
else
   call phaml_recv(phaml_solution%procs,proc,recv_int,ni,recv_real,nr,320)
   do i=1,NAMELEN
      filename(i:i) = char(recv_int(i))
   end do
   deallocate(recv_int,stat=astat)
endif

! break the filename into the root and suffix

sep = index(filename,".",back=.true.)
if (sep == 0) then
   rootname = filename
   suffixname = ""
else
   rootname = filename(1:sep-1)
   suffixname = filename(sep+1:NAMELEN)
endif

! construct the actual filename

ext_filename = rootname
start = len_trim(ext_filename)+1
if (nproc < 10) then
   write(ext_filename(start:start),"(i1)") me
elseif (nproc < 100) then
   if (me < 10) then
      write(ext_filename(start:start),"(a1)") "0"
      write(ext_filename(start+1:start+1),"(i1)") me
   else
      write(ext_filename(start:start+1),"(i2)") me
   endif
elseif (nproc < 1000) then
   if (me < 10) then
      write(ext_filename(start:start+1),"(a2)") "00"
      write(ext_filename(start+2:start+2),"(i1)") me
   elseif (me < 100) then
      write(ext_filename(start:start),"(a1)") "0"
      write(ext_filename(start+1:start+2),"(i2)") me
   else
      write(ext_filename(start:start+2),"(i3)") me
   endif
elseif (nproc < 10000) then
   if (me < 10) then
      write(ext_filename(start:start+2),"(a3)") "000"
      write(ext_filename(start+3:start+3),"(i1)") me
   elseif (me < 100) then
      write(ext_filename(start:start+1),"(a2)") "00"
      write(ext_filename(start+2:start+3),"(i2)") me
   elseif (me < 1000) then
      write(ext_filename(start:start),"(a1)") "0"
      write(ext_filename(start+1:start+3),"(i3)") me
   else
      write(ext_filename(start:start+3),"(i4)") me
   endif
else
! RESTRICTION 9999 processors
   call fatal("cannot create filenames for more than 9999 processors", &
               procs=phaml_solution%procs)
   return
endif
if (sep == 0) then
   outname = ext_filename
else
   outname = trim(ext_filename)//"."//suffixname
endif

end subroutine make_filename

!          ----------------
subroutine master_to_slaves(phaml_solution,iparam,rparam)
!          ----------------

!----------------------------------------------------
! This routine sends the arrays iparam and rparam from the master to
! the slaves.  It is called by the user provided routine update_usermod.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(phaml_solution_type), intent(in) :: phaml_solution
integer, intent(inout) :: iparam(:)
real(my_real), intent(inout) :: rparam(:)
!----------------------------------------------------
! Local variables:

integer, pointer :: irecv(:)
real(my_real), pointer :: rrecv(:)
integer :: nint, nreal, proc, astat, gp
!----------------------------------------------------
! Begin executable code

if (PARALLEL == SEQUENTIAL) return

gp = graphics_proc(phaml_solution%procs)

! the master sends the code for update_usermod (9) and the arrays

if (my_proc(phaml_solution%procs) == MASTER) then
   do proc=1,num_proc(phaml_solution%procs)
      call phaml_send(phaml_solution%procs,proc,(/9/),1,(/0.0_my_real/),0,101)
   end do
   if (gp /= -1) then
      call phaml_send(phaml_solution%procs,gp,(/9/),1,(/0.0_my_real/),0,101)
   endif
   do proc=1,num_proc(phaml_solution%procs)
      call phaml_send(phaml_solution%procs,proc,iparam,size(iparam), &
                      rparam,size(rparam),330)
   end do
   if (gp /= -1) then
      call phaml_send(phaml_solution%procs,gp,iparam,size(iparam), &
                      rparam,size(rparam),330)
   endif

! the slaves receive the arrays and copy them into the parameter arrays, and
! send data to their graphics processes if they exist

elseif (my_proc(phaml_solution%procs) <= num_proc(phaml_solution%procs)) then

   call phaml_recv(phaml_solution%procs,proc,irecv,nint,rrecv,nreal,330)
   nint = min(nint,size(iparam))
   nreal = min(nreal,size(rparam))
   iparam(1:nint) = irecv(1:nint)
   rparam(1:nreal) = rrecv(1:nreal)
   if (associated(irecv)) deallocate(irecv,stat=astat)
   if (associated(rrecv)) deallocate(rrecv,stat=astat)

   if (gp /= -1) then
      call phaml_send(phaml_solution%procs,gp,(/9/),1,(/0.0_my_real/),0,101)
      call phaml_send(phaml_solution%procs,gp,iparam,nint,rparam,nreal,330)
   endif

! graphics receive the arrays and copy them into the parameter arrays

else
   call phaml_recv(phaml_solution%procs,proc,irecv,nint,rrecv,nreal,330)
   nint = min(nint,size(iparam))
   nreal = min(nreal,size(rparam))
   iparam(1:nint) = irecv(1:nint)
   rparam(1:nreal) = rrecv(1:nreal)
   if (associated(irecv)) deallocate(irecv,stat=astat)
   if (associated(rrecv)) deallocate(rrecv,stat=astat)
endif

end subroutine master_to_slaves

end module phaml
