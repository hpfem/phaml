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

!----------------------------------------------------
! This file contains the user supplied external subroutines that define
! the PDE(s) to be solved, and other required external subroutines.
!   pdecoefs bconds boundary_point boundary_npiece boundary_param iconds trues
!   truexs trueys update_usermod phaml_integral_kernel regularity
!
! This version allows the subroutines to be changed at run time.
! Module pde_intf provides this capability through
!   1) subroutine init_pde_intf which MUST be called first
!   2) subroutine set_<user routine>(<user routine>_sub) for each of
!      the user supplied routines
!
! It also provides a second interface to PHAML's public routines which
! accepts an integer instead of the type(phaml_solution_type) argument.
! The integer is an index into an array of type(phaml_solution_type) that
! is maintained in this module.  The size of the array is given by the
! symbolic constant MAX_PHAML_SOLUTIONS defined below.  In most cases,
! the array only needs to be of length 1 and the integer argument is 1.
!
! I thank Ruediger Kessel for his suggestions for this "dynamic binding".
!----------------------------------------------------

!---------------------------------------------------------------------!
module pde_intf

! Several compilers have a bug where my use of generic interfaces for the
! phaml routines causes confusion.  Renaming is a workaround for that bug.

use phaml, phaml_compress_orig=>phaml_compress, &
           phaml_copy_soln_to_old_orig=>phaml_copy_soln_to_old, &
           phaml_create_orig=>phaml_create, &
           phaml_destroy_orig=>phaml_destroy, &
           phaml_evaluate_orig=>phaml_evaluate, &
           phaml_integrate_orig=>phaml_integrate, &
           phaml_pclose_orig=>phaml_pclose, &
           phaml_popen_orig=>phaml_popen, &
           phaml_query_orig=>phaml_query, &
           phaml_restore_orig=>phaml_restore, &
           phaml_scale_orig=>phaml_scale, &
           phaml_solve_pde_orig=>phaml_solve_pde, &
           phaml_store_orig=>phaml_store, &
           phaml_store_matrix_orig=>phaml_store_matrix
implicit none

!---------------------------------------------------------------------!
! The maximum number of type(phaml_solution_type) variables allowed
! through the interfaces that replace type(phaml_solution_type) with
! an integer.  This can be increased as needed.

integer, parameter, private :: MAX_PHAML_SOLUTIONS = 1

!---------------------------------------------------------------------!
! Interface blocks for the pointers to user provided subroutines

interface

   subroutine pdecoefs_proc(x,y,cxx,cxy,cyy,cx,cy,c,rs)
   double precision, intent(in) :: x,y
   double precision, intent(out) :: cxx(:,:),cxy(:,:),cyy(:,:),cx(:,:), &
                                    cy(:,:),c(:,:),rs(:)
   end subroutine pdecoefs_proc

   subroutine bconds_proc(x,y,bmark,itype,c,rs)
   double precision, intent(in) :: x,y
   integer, intent(in) :: bmark
   integer, intent(out) :: itype(:)
   double precision, intent(out) :: c(:,:),rs(:)
   end subroutine bconds_proc

   function iconds_proc(x,y,comp,eigen)
   double precision, intent(in) :: x,y
   integer, intent(in) :: comp,eigen
   double precision :: iconds_proc
   end function iconds_proc

   function trues_proc(x,y,comp,eigen)
   double precision, intent(in) :: x,y
   integer, intent(in) :: comp,eigen
   double precision :: trues_proc
   end function trues_proc

   function truexs_proc(x,y,comp,eigen)
   double precision, intent(in) :: x,y
   integer, intent(in) :: comp,eigen
   double precision :: truexs_proc
   end function truexs_proc

   function trueys_proc(x,y,comp,eigen)
   double precision, intent(in) :: x,y
   integer, intent(in) :: comp,eigen
   double precision :: trueys_proc
   end function trueys_proc

   subroutine boundary_point_proc(ipiece,s,x,y)
   integer, intent(in) :: ipiece
   double precision, intent(in) :: s
   double precision, intent(out) :: x,y
   end subroutine boundary_point_proc

   function boundary_npiece_proc(hole)
   integer, intent(in) :: hole
   integer :: boundary_npiece_proc
   end function boundary_npiece_proc

   subroutine boundary_param_proc(start,finish)
   double precision, intent(out) :: start(:), finish(:)
   end subroutine boundary_param_proc

   subroutine update_usermod_proc(phaml_solution)
   use phaml
   type(phaml_solution_type), intent(in) :: phaml_solution
   end subroutine update_usermod_proc

   function phaml_integral_kernel_proc(kernel,x,y)
   integer, intent(in) :: kernel
   double precision, intent(in) :: x,y
   double precision :: phaml_integral_kernel_proc
   end function phaml_integral_kernel_proc

   function regularity_proc(x,y)
   double precision, intent(in) :: x(3),y(3)
   double precision :: regularity_proc
   end function regularity_proc

end interface

!---------------------------------------------------------------------!
! Pointers to the user provided routines

procedure(pdecoefs_proc), pointer, save :: pdecoefs_ptr
procedure(bconds_proc), pointer, save :: bconds_ptr
procedure(iconds_proc), pointer, save :: iconds_ptr
procedure(trues_proc), pointer, save :: trues_ptr
procedure(truexs_proc), pointer, save :: truexs_ptr
procedure(trueys_proc), pointer, save :: trueys_ptr
procedure(boundary_point_proc), pointer, save :: boundary_point_ptr
procedure(boundary_npiece_proc), pointer, save :: boundary_npiece_ptr
procedure(boundary_param_proc), pointer, save :: boundary_param_ptr
procedure(update_usermod_proc), pointer, save :: update_usermod_ptr
procedure(phaml_integral_kernel_proc), pointer, save :: phaml_integral_kernel_ptr
procedure(regularity_proc), pointer, save :: regularity_ptr

!---------------------------------------------------------------------!
! An array of type(phaml_solution_type) that can be used instead of
! having the main program declare the type(phaml_solution_type) variable(s).

type(phaml_solution_type), save, private :: phaml_soluts(MAX_PHAML_SOLUTIONS)

!---------------------------------------------------------------------!
! Generic interfaces to overload PHAML's public routines with routines
! that accept an integer instead of the type(phaml_solution_type)
! argument.  The integer is an index into phaml_soluts.

interface phaml_compress
   module procedure phaml_compress_orig, phaml_compress_int
end interface

! don't need phaml_connect because it already takes integers

interface phaml_copy_soln_to_old
   module procedure phaml_copy_soln_to_old_orig, phaml_copy_soln_to_old_int
end interface

interface phaml_create
   module procedure phaml_create_orig, phaml_create_int
end interface

interface phaml_destroy
   module procedure phaml_destroy_orig, phaml_destroy_int
end interface

interface phaml_evaluate
   module procedure phaml_evaluate_orig, phaml_evaluate_int
end interface

! don't need phaml_evaluate_old because it doesn't have phaml_solution_type

interface phaml_integrate
   module procedure phaml_integrate_orig, phaml_integrate_int
end interface

interface phaml_pclose
   module procedure phaml_pclose_orig, phaml_pclose_int
end interface

interface phaml_popen
   module procedure phaml_popen_orig, phaml_popen_int
end interface

interface phaml_query
   module procedure phaml_query_orig, phaml_query_int
end interface

interface phaml_restore
   module procedure phaml_restore_orig, phaml_restore_int
end interface

interface phaml_scale
   module procedure phaml_scale_orig, phaml_scale_int
end interface

interface phaml_solve_pde
   module procedure phaml_solve_pde_orig, phaml_solve_pde_int
end interface

interface phaml_store
   module procedure phaml_store_orig, phaml_store_int
end interface

interface phaml_store_matrix
   module procedure phaml_store_matrix_orig, phaml_store_matrix_int
end interface

contains

!---------------------------------------------------------------------!
! Initialization routine

subroutine init_pde_intf

! Interfaces for the default user provided routines

interface

   subroutine pdecoefs_def(x,y,cxx,cxy,cyy,cx,cy,c,rs)
   double precision, intent(in) :: x,y
   double precision, intent(out) :: cxx(:,:),cxy(:,:),cyy(:,:),cx(:,:), &
                                    cy(:,:),c(:,:),rs(:)
   end subroutine pdecoefs_def

   subroutine bconds_def(x,y,bmark,itype,c,rs)
   double precision, intent(in) :: x,y
   integer, intent(in) :: bmark
   integer, intent(out) :: itype(:)
   double precision, intent(out) :: c(:,:),rs(:)
   end subroutine bconds_def

   function iconds_def(x,y,comp,eigen)
   double precision, intent(in) :: x,y
   integer, intent(in) :: comp,eigen
   double precision :: iconds_def
   end function iconds_def

   function trues_def(x,y,comp,eigen)
   double precision, intent(in) :: x,y
   integer, intent(in) :: comp,eigen
   double precision :: trues_def
   end function trues_def

   function truexs_def(x,y,comp,eigen)
   double precision, intent(in) :: x,y
   integer, intent(in) :: comp,eigen
   double precision :: truexs_def
   end function truexs_def

   function trueys_def(x,y,comp,eigen)
   double precision, intent(in) :: x,y
   integer, intent(in) :: comp,eigen
   double precision :: trueys_def
   end function trueys_def

   subroutine boundary_point_def(ipiece,s,x,y)
   integer, intent(in) :: ipiece
   double precision, intent(in) :: s
   double precision, intent(out) :: x,y
   end subroutine boundary_point_def

   function boundary_npiece_def(hole)
   integer, intent(in) :: hole
   integer :: boundary_npiece_def
   end function boundary_npiece_def

   subroutine boundary_param_def(start,finish)
   double precision, intent(out) :: start(:), finish(:)
   end subroutine boundary_param_def

   subroutine update_usermod_def(phaml_solution)
   use phaml
   type(phaml_solution_type), intent(in) :: phaml_solution
   end subroutine update_usermod_def

   function phaml_integral_kernel_def(kernel,x,y)
   integer, intent(in) :: kernel
   double precision, intent(in) :: x,y
   double precision :: phaml_integral_kernel_def
   end function phaml_integral_kernel_def

   function regularity_def(x,y)
   double precision, intent(in) :: x(3),y(3)
   double precision :: regularity_def
   end function regularity_def

end interface

! Assign default user provided routines

pdecoefs_ptr => pdecoefs_def
bconds_ptr => bconds_def
iconds_ptr => iconds_def
trues_ptr => trues_def
truexs_ptr => truexs_def
trueys_ptr => trueys_def
boundary_point_ptr => boundary_point_def
boundary_npiece_ptr => boundary_npiece_def
boundary_param_ptr => boundary_param_def
update_usermod_ptr => update_usermod_def
phaml_integral_kernel_ptr => phaml_integral_kernel_def
regularity_ptr => regularity_def

end subroutine init_pde_intf

!---------------------------------------------------------------------!
! Subroutines to change the user provided routines

subroutine set_pdecoefs(pdecoefs_sub)
interface
   subroutine pdecoefs_sub(x,y,cxx,cxy,cyy,cx,cy,c,rs)
   double precision, intent(in) :: x,y
   double precision, intent(out) :: cxx(:,:),cxy(:,:),cyy(:,:),cx(:,:), &
                                    cy(:,:),c(:,:),rs(:)
   end subroutine pdecoefs_sub
end interface
pdecoefs_ptr => pdecoefs_sub
end subroutine set_pdecoefs

subroutine set_bconds(bconds_sub)
interface
   subroutine bconds_sub(x,y,bmark,itype,c,rs)
   double precision, intent(in) :: x,y
   integer, intent(in) :: bmark
   integer, intent(out) :: itype(:)
   double precision, intent(out) :: c(:,:),rs(:)
   end subroutine bconds_sub
end interface
bconds_ptr => bconds_sub
end subroutine set_bconds

subroutine set_iconds(iconds_sub)
interface
   function iconds_sub(x,y,comp,eigen)
   double precision, intent(in) :: x,y
   integer, intent(in) :: comp,eigen
   double precision :: iconds_sub
   end function iconds_sub
end interface
iconds_ptr => iconds_sub
end subroutine set_iconds

subroutine set_trues(trues_sub)
interface
   function trues_sub(x,y,comp,eigen)
   double precision, intent(in) :: x,y
   integer, intent(in) :: comp,eigen
   double precision :: trues_sub
   end function trues_sub
end interface
trues_ptr => trues_sub
end subroutine set_trues

subroutine set_truexs(truexs_sub)
interface
   function truexs_sub(x,y,comp,eigen)
   double precision, intent(in) :: x,y
   integer, intent(in) :: comp,eigen
   double precision :: truexs_sub
   end function truexs_sub
end interface
truexs_ptr => truexs_sub
end subroutine set_truexs

subroutine set_trueys(trueys_sub)
interface
   function trueys_sub(x,y,comp,eigen)
   double precision, intent(in) :: x,y
   integer, intent(in) :: comp,eigen
   double precision :: trueys_sub
   end function trueys_sub
end interface
trueys_ptr => trueys_sub
end subroutine set_trueys

subroutine set_boundary_point(boundary_point_sub)
interface
   subroutine boundary_point_sub(ipiece,s,x,y)
   integer, intent(in) :: ipiece
   double precision, intent(in) :: s
   double precision, intent(out) :: x,y
   end subroutine boundary_point_sub
end interface
boundary_point_ptr => boundary_point_sub
end subroutine set_boundary_point

subroutine set_boundary_npiece(boundary_npiece_sub)
interface
   function boundary_npiece_sub(hole)
   integer, intent(in) :: hole
   integer :: boundary_npiece_sub
   end function boundary_npiece_sub
end interface
boundary_npiece_ptr => boundary_npiece_sub
end subroutine set_boundary_npiece

subroutine set_boundary_param(boundary_param_sub)
interface
   subroutine boundary_param_sub(start,finish)
   double precision, intent(out) :: start(:), finish(:)
   end subroutine boundary_param_sub
end interface
boundary_param_ptr => boundary_param_sub
end subroutine set_boundary_param

subroutine set_update_usermod(update_usermod_sub)
interface
   subroutine update_usermod_sub(phaml_solution)
   use phaml
   type(phaml_solution_type), intent(in) :: phaml_solution
   end subroutine update_usermod_sub
end interface
update_usermod_ptr => update_usermod_sub
end subroutine set_update_usermod

subroutine set_phaml_integral_kernel(phaml_integral_kernel_sub)
interface
   function phaml_integral_kernel_sub(kernel,x,y)
   integer, intent(in) :: kernel
   double precision, intent(in) :: x,y
   double precision :: phaml_integral_kernel_sub
   end function phaml_integral_kernel_sub
end interface
phaml_integral_kernel_ptr => phaml_integral_kernel_sub
end subroutine set_phaml_integral_kernel

subroutine set_regularity(regularity_sub)
interface
   function regularity_sub(x,y)
   double precision, intent(in) :: x(3),y(3)
   double precision :: regularity_sub
   end function regularity_sub
end interface
regularity_ptr => regularity_sub
end subroutine set_regularity

!---------------------------------------------------------------------!
! PHAML's public routines that accept an integer instead of a
! type(phaml_solution_type) variable.

subroutine phaml_compress_int(index)
integer, intent(in) :: index
call phaml_compress(phaml_soluts(index))
end subroutine phaml_compress_int

! don't need phaml_connect because it already takes integers

subroutine phaml_copy_soln_to_old_int(index)
integer, intent(in) :: index
call phaml_copy_soln_to_old(phaml_soluts(index))
end subroutine phaml_copy_soln_to_old_int

subroutine phaml_create_int(index,nproc,draw_grid_who, &
                            spawn_form,debug_command,display,graphics_host, &
                            output_unit,error_unit,output_now,id,system_size, &
                            eq_type,max_blen,triangle_files,update_umod)
integer, intent(in) :: index
integer, optional :: nproc,draw_grid_who,spawn_form,output_unit,error_unit, &
                     output_now,id,system_size,eq_type
character(len=*), optional, intent(in) :: debug_command,display
character(len=*), optional, intent(in) :: graphics_host
real(my_real), optional, intent(in) :: max_blen
character(len=*), optional, intent(in) :: triangle_files
logical, optional, intent(in) :: update_umod
if (index > MAX_PHAML_SOLUTIONS) then
   print *
   print *,"-----------------------------------------------------------------"
   print *,"ERROR -- phaml_create input index = ",index," is larger"
   print *,"         than MAX_PHAML_SOLUTIONS.  Increase MAX_PHAML_SOLUTIONS"
   print *,"         in file pde.f90."
   print *,"-----------------------------------------------------------------"
   stop
endif
call phaml_create(phaml_soluts(index),nproc,draw_grid_who, &
                  spawn_form,debug_command,display,graphics_host, &
                  output_unit,error_unit,output_now,id,system_size, &
                  eq_type,max_blen,triangle_files,update_umod)
end subroutine phaml_create_int

subroutine phaml_destroy_int(index,finalize_mpi)
integer, intent(in) :: index
logical, intent(in), optional :: finalize_mpi
call phaml_destroy(phaml_soluts(index),finalize_mpi)
end subroutine phaml_destroy_int

subroutine phaml_evaluate_int(index,x,y,u,ux,uy,uxx,uyy,comp,eigen)
integer, intent(in) :: index
real(my_real), intent(in) :: x(:),y(:)
real(my_real), optional, intent(out) :: u(:),ux(:),uy(:),uxx(:),uyy(:)
integer, optional, intent(in) :: comp,eigen
call phaml_evaluate(phaml_soluts(index),x,y,u,ux,uy,uxx,uyy,comp,eigen)
end subroutine phaml_evaluate_int

! don't need phaml_evaluate_old because it doesn't have phaml_solution_type

function phaml_integrate_int(index,kernel,comp1,eigen1,comp2,eigen2,p,q)
integer, intent(in) :: index
integer, intent(in) :: kernel
integer, optional, intent(in) :: comp1, eigen1, comp2, eigen2, p, q
real(my_real) :: phaml_integrate_int
phaml_integrate_int = phaml_integrate(phaml_soluts(index),kernel,comp1,eigen1, &
                                      comp2,eigen2,p,q)
end function phaml_integrate_int

subroutine phaml_pclose_int(index,unit)
integer, intent(in) :: index
integer, intent(in) :: unit
call phaml_pclose(phaml_soluts(index),unit)
end subroutine phaml_pclose_int

subroutine phaml_popen_int(index,unit,file,form)
integer, intent(in) :: index
integer, intent(in) :: unit
character(len=*), intent(in) :: file
character(len=*), intent(in), optional :: form
call phaml_popen(phaml_soluts(index),unit,file,form)
end subroutine phaml_popen_int

subroutine phaml_query_int(index,nvert,nvert_proc,nvert_own,nelem, &
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
integer, intent(in) :: index
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
call phaml_query(phaml_soluts(index),nvert,nvert_proc,nvert_own,nelem, &
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
end subroutine phaml_query_int

subroutine phaml_restore_int(index,unit,do_draw_grid,pause)
integer, intent(in) :: index
integer, intent(in) :: unit
logical, intent(in), optional :: do_draw_grid, pause
call phaml_restore(phaml_soluts(index),unit,do_draw_grid,pause)
end subroutine phaml_restore_int

subroutine phaml_scale_int(index,factor,comp,eigen)
integer, intent(in) :: index
real(my_real), intent(in) :: factor
integer, intent(in), optional :: comp,eigen
call phaml_scale(phaml_soluts(index),factor,comp,eigen)
end subroutine phaml_scale_int

subroutine phaml_solve_pde_int(index, iterm, max_elem, max_vert, max_eq, &
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
integer, intent(in) :: index
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
call phaml_solve_pde(phaml_soluts(index), iterm, max_elem, max_vert, max_eq, &
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
end subroutine phaml_solve_pde_int

subroutine phaml_store_int(index,unit)
integer, intent(in) :: index
integer, intent(in) :: unit
call phaml_store(phaml_soluts(index),unit)
end subroutine phaml_store_int

subroutine phaml_store_matrix_int(index,stiffness_unit,rhs_unit,mass_unit,&
                                  inc_quad_order)
integer, intent(in) :: index
integer, optional, intent(in) :: stiffness_unit, rhs_unit, mass_unit, &
                                 inc_quad_order
call phaml_store_matrix(phaml_soluts(index),stiffness_unit,rhs_unit,mass_unit,&
                        inc_quad_order)
end subroutine phaml_store_matrix_int

end module pde_intf

!---------------------------------------------------------------------!
! External default user provided routines define a Laplace equation with
! homogeneous Dirichlet boundary conditions

subroutine pdecoefs_def(x,y,cxx,cxy,cyy,cx,cy,c,rs)
double precision, intent(in) :: x,y
double precision, intent(out) :: cxx(:,:),cxy(:,:),cyy(:,:),cx(:,:),cy(:,:), &
                              c(:,:),rs(:)
cxx = 1; cyy = 1
cxy = 0; cx = 0; cy = 0; c = 0
rs = 0
end subroutine pdecoefs_def

subroutine bconds_def(x,y,bmark,itype,c,rs)
use phaml
double precision, intent(in) :: x,y
integer, intent(in) :: bmark
integer, intent(out) :: itype(:)
double precision, intent(out) :: c(:,:),rs(:)
itype = DIRICHLET
c = 0
rs = 0
end subroutine bconds_def

function iconds_def(x,y,comp,eigen)
double precision, intent(in) :: x,y
integer, intent(in) :: comp,eigen
double precision :: iconds_def
iconds_def = 0
end function iconds_def

function trues_def(x,y,comp,eigen)
double precision, intent(in) :: x,y
integer, intent(in) :: comp,eigen
double precision :: trues_def
trues_def = 0
end function trues_def

function truexs_def(x,y,comp,eigen)
double precision, intent(in) :: x,y
integer, intent(in) :: comp,eigen
double precision :: truexs_def
truexs_def = 0
end function truexs_def

function trueys_def(x,y,comp,eigen)
double precision, intent(in) :: x,y
integer, intent(in) :: comp,eigen
double precision :: trueys_def
trueys_def = 0
end function trueys_def

subroutine boundary_point_def(ipiece,s,x,y)
integer, intent(in) :: ipiece
double precision, intent(in) :: s
double precision, intent(out) :: x,y
x = 0; y = 0
end subroutine boundary_point_def

function boundary_npiece_def(hole)
integer, intent(in) :: hole
integer :: boundary_npiece_def
boundary_npiece_def = 0
end function boundary_npiece_def

subroutine boundary_param_def(start,finish)
double precision, intent(out) :: start(:), finish(:)
start = 0; finish = 0
end subroutine boundary_param_def

subroutine update_usermod_def(phaml_solution)
use phaml
type(phaml_solution_type), intent(in) :: phaml_solution
end subroutine update_usermod_def

function phaml_integral_kernel_def(kernel,x,y)
integer, intent(in) :: kernel
double precision, intent(in) :: x,y
double precision :: phaml_integral_kernel_def
phaml_integral_kernel_def = 1
end function phaml_integral_kernel_def

function regularity_def(x,y)
double precision, intent(in) :: x(3),y(3)
double precision :: regularity_def
regularity_def = huge(0.0d0)
end function regularity_def

!---------------------------------------------------------------------!
! External routines that are called by PHAML, call the pointed to routine

subroutine pdecoefs(x,y,cxx,cxy,cyy,cx,cy,c,rs)
use pde_intf
implicit none
double precision, intent(in) :: x,y
double precision, intent(out) :: cxx(:,:),cxy(:,:),cyy(:,:),cx(:,:),cy(:,:), &
                              c(:,:),rs(:)
call pdecoefs_ptr(x,y,cxx,cxy,cyy,cx,cy,c,rs)
end subroutine pdecoefs

subroutine bconds(x,y,bmark,itype,c,rs)
use pde_intf
implicit none
double precision, intent(in) :: x,y
integer, intent(in) :: bmark
integer, intent(out) :: itype(:)
double precision, intent(out) :: c(:,:),rs(:)
call bconds_ptr(x,y,bmark,itype,c,rs)
end subroutine bconds

function iconds(x,y,comp,eigen)
use pde_intf
implicit none
double precision, intent(in) :: x,y
integer, intent(in) :: comp,eigen
double precision :: iconds
iconds = iconds_ptr(x,y,comp,eigen)
end function iconds

function trues(x,y,comp,eigen)
use pde_intf
implicit none
double precision, intent(in) :: x,y
integer, intent(in) :: comp,eigen
double precision :: trues
trues = trues_ptr(x,y,comp,eigen)
end function trues

function truexs(x,y,comp,eigen)
use pde_intf
implicit none
double precision, intent(in) :: x,y
integer, intent(in) :: comp,eigen
double precision :: truexs
truexs = truexs_ptr(x,y,comp,eigen)
end function truexs

function trueys(x,y,comp,eigen)
use pde_intf
implicit none
double precision, intent(in) :: x,y
integer, intent(in) :: comp,eigen
double precision :: trueys
trueys = trueys_ptr(x,y,comp,eigen)
end function trueys

subroutine boundary_point(ipiece,s,x,y)
use pde_intf
implicit none
integer, intent(in) :: ipiece
double precision, intent(in) :: s
double precision, intent(out) :: x,y
call boundary_point_ptr(ipiece,s,x,y)
end subroutine boundary_point

function boundary_npiece(hole)
use pde_intf
implicit none
integer, intent(in) :: hole
integer :: boundary_npiece
boundary_npiece = boundary_npiece_ptr(hole)
end function boundary_npiece

subroutine boundary_param(start,finish)
use pde_intf
implicit none
double precision, intent(out) :: start(:), finish(:)
call boundary_param_ptr(start,finish)
end subroutine boundary_param

subroutine update_usermod(phaml_solution)
use phaml
use pde_intf
type(phaml_solution_type), intent(in) :: phaml_solution
call update_usermod_ptr(phaml_solution)
end subroutine update_usermod

function phaml_integral_kernel(kernel,x,y)
use pde_intf
integer, intent(in) :: kernel
double precision, intent(in) :: x,y
double precision :: phaml_integral_kernel
phaml_integral_kernel = phaml_integral_kernel_ptr(kernel,x,y)
end function phaml_integral_kernel

function regularity(x,y)
use pde_intf
double precision, intent(in) :: x(3),y(3)
double precision :: regularity
regularity = regularity_ptr(x,y)
end function regularity
