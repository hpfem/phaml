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

module linear_system

!----------------------------------------------------
! This module contains data structures and routines for
! manipulating the linear system of equations
!
! communication tags in this module are of the form 10xx
!----------------------------------------------------

!----------------------------------------------------
! Other modules used are:

   use global
   use message_passing
   use stopwatch
   use linsystype_mod
   use gridtype_mod
   use petsc_interf
   use mumps_interf
   use superlu_interf
   use eigen
   use hbmg
   use krylov
   use lapack_solve
   use make_linsys
   use linsys_io
   use hypre_interf
   use linsys_util
   use hash_mod
   use grid_util

!----------------------------------------------------

implicit none
private
public solve

!----------------------------------------------------
! Non-module procedures used are:

interface

   function iconds(x,y,comp,eigen)
   use global
   real(my_real), intent(in) :: x,y
   integer, intent(in) :: comp, eigen
   real(my_real) :: iconds
   end function iconds

end interface

!----------------------------------------------------

contains

!----------------------------------------------------------------
! CODE FOR SOLVING THE EQUATIONS
!----------------------------------------------------------------

!          -----
subroutine solve(grid,procs,io_cntl,solver_cntl,still_sequential,set_init, &
                 no_time)
!          -----

!----------------------------------------------------
! This routine creates the linear system, solves it putting the solution
! into the vertices of the grid, and deletes the linear system
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout), target :: grid
type(proc_info), intent(in), target :: procs
type(io_options), intent(in), target :: io_cntl
type(solver_options), intent(in), target :: solver_cntl
logical, intent(in) :: still_sequential, set_init
logical, intent(in), optional :: no_time

!----------------------------------------------------
! Local variables:

integer :: i,j,lev,vert,elem,edge,eq,objtype,brank,srank,lid,loc_neq,info,ss
type(linsys_type) :: linear_system
real(my_real), pointer :: block(:,:)
integer, pointer :: ipiv(:)
logical :: timeit

!----------------------------------------------------
! Begin executable code

grid_changed = .true. ! says we need to resend graphics data

! if set_init is true, set the solution from iconds and return

if (set_init) then
   call set_init_cond(grid)
   return
endif

! create the linear system

call create_linear_system(linear_system,grid,procs,solver_cntl, &
                          io_cntl,still_sequential,no_time)
if (ierr /= 0) then
   call warning("create_linear_system returned nonzero error code", &
                intlist=(/ierr/))
    return
endif

!call count_memory(linear_system)

! if printing error TOO_MUCH, then we are probably monitoring the
! convergence of the solver.  Set the initial guess to 0.0 instead
! of the solution of the previous grid to give more iterations
! before the solution error is smaller than the discretization error.

if (io_cntl%print_error_when == TOO_MUCH) then
   if (my_proc(procs) /= MASTER) then
      where (linear_system%equation_type /= DIRICHLET) &
         linear_system%solution(1:) = 0.0_my_real
   else
      call warning("setting solution initial guess to 0 for print_error_when==TOO_MUCH")
   endif
endif

! start timing the solution process

if (present(no_time)) then
   timeit = .not. no_time
else
   timeit = .true.
endif
if (timeit) then
   call reset_watch((/psolve,cpsolve/))
   call start_watch((/psolve,tsolve/))
endif

! solve the linear system

select case (solver_cntl%eq_type)

case (ELLIPTIC)

   select case (solver_cntl%solver)
   case (MG_SOLVER)
      call multigrid(grid,procs,linear_system,io_cntl,solver_cntl, &
                     still_sequential)
   case (CG_SOLVER)
      call phaml_cg(linear_system,procs,grid,io_cntl,solver_cntl, &
                    still_sequential)
   case (GMRES_SOLVER)
      call phaml_gmres(linear_system,procs,grid,io_cntl,solver_cntl, &
                       still_sequential)
   case (LAPACK_INDEFINITE_SOLVER)
      if (my_proc(procs) /= MASTER) then
         call make_lapack_gen_band(linear_system%nlev+1,linear_system, &
                                   linear_system%lapack_mat)
         linear_system%lapack_gen_band_exists = .true.
         call lapack_indef(linear_system%nlev+1,linear_system,linear_system%lapack_mat)
      endif
   case (LAPACK_SPD_SOLVER)
      if (my_proc(procs) /= MASTER) then
         call make_lapack_symm_band(linear_system%nlev+1,linear_system, &
                                    linear_system%lapack_mat)
         if (ierr /= NO_ERROR) return
         linear_system%lapack_symm_band_exists = .true.
         call lapack_spd(linear_system%nlev+1,linear_system,linear_system%lapack_mat)
      endif
   case (PETSC_RICHARDSON_SOLVER, PETSC_CHEBYCHEV_SOLVER, PETSC_CG_SOLVER, &
         PETSC_GMRES_SOLVER,      PETSC_TCQMR_SOLVER,     PETSC_BCGS_SOLVER, &
         PETSC_CGS_SOLVER,        PETSC_TFQMR_SOLVER,     PETSC_CR_SOLVER, &
         PETSC_LSQR_SOLVER,       PETSC_BICG_SOLVER)
      if (my_proc(procs) /= MASTER) then
         where (linear_system%equation_type == DIRICHLET) &
            linear_system%rhs = linear_system%solution(1:)
         if (solver_cntl%petsc_matrix_free) then
            call create_petsc_linear_system_mf(linear_system, &
                                               linear_system%petsc_matrix, &
                                               still_sequential)
         else
            call create_petsc_linear_system(linear_system, &
                                            linear_system%petsc_matrix, &
                                            still_sequential,procs)
         endif
         call petsc_solve(linear_system,linear_system%petsc_matrix, &
                          solver_cntl,io_cntl,still_sequential,grid,procs)
         call destroy_petsc_linear_system(linear_system%petsc_matrix)
      endif
   case (MUMPS_SPD_SOLVER, MUMPS_GEN_SOLVER, MUMPS_NONSYM_SOLVER)
      if (my_proc(procs) /= MASTER) then
         call create_mumps_linear_system(linear_system%mumps_matrix, &
                                         linear_system,procs, &
                                         solver_cntl,still_sequential)
         call mumps_solve(linear_system%mumps_matrix,linear_system, &
                          procs,still_sequential)
         call destroy_mumps_linear_system(linear_system%mumps_matrix,procs)
      endif
   case (SUPERLU_SOLVER)
      if (my_proc(procs) /= MASTER) then
         call create_superlu_linear_system(linear_system%superlu_matrix, &
                                           linear_system,procs, &
                                           still_sequential)
         call superlu_solve(linear_system%superlu_matrix,linear_system,&
                            procs,still_sequential)
         call destroy_superlu_linear_system(linear_system%superlu_matrix,procs)
      endif
   case (HYPRE_BOOMERAMG_SOLVER, HYPRE_PCG_SOLVER, HYPRE_GMRES_SOLVER)
      if (my_proc(procs) /= MASTER) then
         call create_hypre_linear_system(linear_system%hypre_matrix, &
                                         linear_system,procs, &
                                         still_sequential)
         call hypre_solve(linear_system%hypre_matrix,linear_system, &
                          procs,solver_cntl,still_sequential)
         call destroy_hypre_linear_system(linear_system%hypre_matrix,procs)
      endif
   case default
      ierr = USER_INPUT_ERROR
      call fatal("illegal selection for solver",procs=procs)
   end select

   if (io_cntl%print_error_when == FREQUENTLY .or. &
       io_cntl%print_error_when == TOO_MUCH) then
      if (my_proc(procs) == MASTER) write(outunit,"(A)") "After solution, condensed matrix:"
      call linsys_residual(linear_system,procs,still_sequential,0,1002,.true., &
                           .true.)
   endif

! return linear system delimiters to the non-condensed matrix; however, the
! matrix values and rhs for rows associated with vertex and edge bases remain
! the statically condensed Schur complement values, unless the stiffness matrix
! was kept

   linear_system%neq = linear_system%neq_vert + linear_system%neq_edge + &
                       linear_system%neq_face
   linear_system%end_row => linear_system%end_row_face
   linear_system%matrix_val => linear_system%stiffness
   linear_system%rhs => linear_system%rhs_nocond

! solve for the face bases

   if (my_proc(procs) /= MASTER) then
      eq = linear_system%begin_level(linear_system%nlev+2)
      do
         if (eq >= linear_system%begin_level(linear_system%nlev+3)) exit
         call eq_to_grid(linear_system,linear_system%gid(eq),objtype,brank, &
                         srank,lid,grid)
         block => linear_system%elem_block(lid)%matrix
         ipiv => linear_system%elem_block(lid)%ipiv
         loc_neq = linear_system%elem_block(lid)%neq
         do i=1,loc_neq
            linear_system%solution(eq+i-1) = linear_system%rhs(eq+i-1)
            do j=linear_system%begin_row(eq+i-1)+1,linear_system%end_row_edge(eq+i-1)
               if (linear_system%column_index(j) == NO_ENTRY) cycle
               if (linear_system%column_index(j) >= eq .and. &
                   linear_system%column_index(j) < eq+loc_neq) cycle
               linear_system%solution(eq+i-1) = linear_system%solution(eq+i-1) - &
                  linear_system%matrix_val(j)*linear_system%solution(linear_system%column_index(j))
            end do
         end do
         if (my_real == kind(0.0)) then
            call sgetrs("N",loc_neq,1,block,loc_neq,ipiv, &
                        linear_system%solution(eq:eq+loc_neq-1),loc_neq,info)
         else
            call dgetrs("N",loc_neq,1,block,loc_neq,ipiv, &
                        linear_system%solution(eq:eq+loc_neq-1),loc_neq,info)
         endif
         if (info /= 0) then
            ierr = PHAML_INTERNAL_ERROR
            call fatal("lapack sgetrs solution failed for solving for face bases", &
                       intlist=(/info/),procs=procs)
            stop
         endif
         eq = eq + loc_neq
      end do
   endif

   if (io_cntl%print_error_when == FREQUENTLY .or. &
       io_cntl%print_error_when == TOO_MUCH) then
      if (my_proc(procs) == MASTER) write(outunit,"(A)") "After solution, full matrix:"
      call linsys_residual(linear_system,procs,still_sequential,0,1003,.true., &
                           .true.)
   endif

case (EIGENVALUE)

   select case (solver_cntl%eigensolver)

   case (ARPACK_SOLVER)

      call eigen_arpack(grid,procs,linear_system,io_cntl,solver_cntl, &
                        still_sequential)

   case (BLOPEX_SOLVER)

      call eigen_blopex(grid,procs,linear_system,io_cntl,solver_cntl, &
                        still_sequential)

   case default

      ierr = USER_INPUT_ERROR
      call fatal("illegal selection for eigensolver",procs=procs)
      stop

   end select

! return linear system delimiters to the non-condensed matrix; however, the
! matrix values in shifted for rows associated with vertex and edge bases remain
! the statically condensed Schur complement values

   linear_system%neq = linear_system%neq_vert + linear_system%neq_edge + &
                       linear_system%neq_face
   linear_system%end_row => linear_system%end_row_face
   linear_system%matrix_val => linear_system%stiffness
   linear_system%rhs => linear_system%rhs_nocond

case default

   ierr = USER_INPUT_ERROR
   call fatal("illegal selection for equation type",procs=procs)

end select

if (my_proc(procs) /= MASTER) then

! Print linear system summary

   call print_linsys_info(linear_system,procs,io_cntl,still_sequential, &
                          (/PHASES,FREQUENTLY,TOO_MUCH/),1001)

! copy the new solution to the grid only if no error has occured

   if (ierr == 0) then

! copy the solution to the grid

      ss = grid%system_size
      do i=1,linear_system%neq
         call eq_to_grid(linear_system,linear_system%gid(i),objtype,brank, &
                         srank,lid,grid)
         select case (objtype)
         case (VERTEX_ID)
            if (grid%vertex_type(lid,srank) /= DIRICHLET .and. &
                grid%vertex_type(lid,srank) /= PERIODIC_MASTER_DIR .and. &
                grid%vertex_type(lid,srank) /= PERIODIC_SLAVE_DIR) then
               grid%vertex_solution(lid,srank,1) = linear_system%solution(i)
               if (solver_cntl%num_eval>1) then
                  grid%vertex_solution(lid,srank,2:solver_cntl%num_eval) = &
                        linear_system%evecs(i,:)
               endif
            endif
         case (EDGE_ID)
            if (grid%edge_type(lid,srank) /= DIRICHLET) then
               grid%edge(lid)%solution(brank,srank,1) = linear_system%solution(i)
               if (solver_cntl%num_eval>1) then
                  grid%edge(lid)%solution(brank,srank,2:solver_cntl%num_eval) =&
                      linear_system%evecs(i,:)
               endif
            endif
         case (ELEMENT_ID)
            grid%element(lid)%solution(brank,srank,1)=linear_system%solution(i)
            if (solver_cntl%num_eval>1) then
               grid%element(lid)%solution(brank,srank,2:solver_cntl%num_eval) =&
                        linear_system%evecs(i,:)
            endif
         end select
      end do

! copy periodic boundary solutions to the slave periodic vertices and edges

      do lev=1,grid%nlev
         vert = grid%head_level_vert(lev)
         do while (vert /= END_OF_LIST)
            do i=1,grid%system_size
               if (grid%vertex_type(vert,i) == PERIODIC_SLAVE .or. &
                   grid%vertex_type(vert,i) == PERIODIC_SLAVE_NAT .or. &
                   grid%vertex_type(vert,i) == PERIODIC_SLAVE_MIX) then
                  grid%vertex_solution(vert,i,:) = &
                     grid%vertex_solution(grid%vertex(vert)%next,i,:)
               endif
            end do
            vert = grid%vertex(vert)%next
         end do
      end do

      do lev=1,grid%nlev
         elem = grid%head_level_elem(lev)
         do while (elem /= END_OF_LIST)
            if (grid%element(elem)%isleaf) then
               do i=1,EDGES_PER_ELEMENT
                  edge = grid%element(elem)%edge(i)
                  if (grid%edge(edge)%degree < 2) cycle
                  do j=1,grid%system_size
                     if (grid%edge_type(edge,j) == PERIODIC_SLAVE) then
                        grid%edge(edge)%solution(:,j,:) = &
                                grid%edge(grid%edge(edge)%next)%solution(:,j,:)
                     endif
                  end do
               end do
            endif
            elem = grid%element(elem)%next
         end do
      end do

      grid%errind_up2date = .false.

   endif

! destroy any extra forms of the linear system

   if (linear_system%coarse_band_exists) then
      call destroy_lapack_band(linear_system%coarse_matrix)
      linear_system%coarse_band_exists = .false.
   endif
   if (linear_system%lapack_symm_band_exists) then
      call destroy_lapack_band(linear_system%lapack_mat)
      linear_system%lapack_symm_band_exists = .false.
   endif
   if (linear_system%lapack_gen_band_exists) then
      call destroy_lapack_band(linear_system%lapack_mat)
      linear_system%lapack_gen_band_exists = .false.
   endif
   if (linear_system%petsc_matrix_exists) then
      call destroy_petsc_linear_system(linear_system%petsc_matrix)
      linear_system%petsc_matrix_exists = .false.
   endif
   if (linear_system%hypre_matrix_exists) then
      call destroy_hypre_linear_system(linear_system%hypre_matrix,procs)
      linear_system%hypre_matrix_exists = .false.
   endif
   if (linear_system%mumps_matrix_exists) then
      call destroy_mumps_linear_system(linear_system%mumps_matrix,procs)
      linear_system%mumps_matrix_exists = .false.
   endif
   if (linear_system%superlu_matrix_exists) then
      call destroy_superlu_linear_system(linear_system%superlu_matrix,procs)
      linear_system%superlu_matrix_exists = .false.
   endif

else ! MASTER

   call print_linsys_info(linear_system,procs,io_cntl,still_sequential, &
                          (/PHASES,FREQUENTLY,TOO_MUCH/),1001)

endif

! destroy the linear system

call destroy_linear_system(linear_system)

! stop timing the solution process

if (timeit) call stop_watch((/psolve,tsolve/))

end subroutine solve

!          -------------
subroutine set_init_cond(grid)
!          -------------

!----------------------------------------------------
! This routine sets the solution from the function iconds
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(grid_type), intent(inout) :: grid
!----------------------------------------------------
! Local variables:

integer :: lev, vert, elem, i, j
!----------------------------------------------------
! Begin executable code

! set the linear (vertex) basis function coefficients

do lev=1,grid%nlev
   vert = grid%head_level_vert(lev)
   do while (vert /= END_OF_LIST)
      do i=1,grid%system_size
         if (grid%vertex_type(vert,i) /= DIRICHLET .and. &
             grid%vertex_type(vert,i) /= PERIODIC_MASTER_DIR .and. &
             grid%vertex_type(vert,i) /= PERIODIC_SLAVE_DIR) then
            do j=1,max(1,grid%num_eval)
               grid%vertex_solution(vert,i,j) = &
                  iconds(grid%vertex(vert)%coord%x, grid%vertex(vert)%coord%y, &
                         i,j)
            end do
         endif
      end do
      vert = grid%vertex(vert)%next
   end do
end do

! go through the elements to set all the edge basis function coefficients,
! and face (element) basis function coefficients

do lev=1,grid%nlev
   elem = grid%head_level_elem(lev)
   do while (elem /= END_OF_LIST)
      do i=1,grid%nsoln
         do j=1,EDGES_PER_ELEMENT
            if (grid%edge_type(grid%element(elem)%edge(j),i) /= DIRICHLET) then
               call edge_exact(grid,grid%element(elem)%edge(j),i,"i")
            endif
         end do
         call elem_exact(grid,elem,i,"i")
      end do
      elem = grid%element(elem)%next
   end do
end do

grid%errind_up2date = .false.

end subroutine set_init_cond

!          ------------
subroutine count_memory(ls)
!          ------------

!----------------------------------------------------
! This routine counts the memory in linear system ls
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(linsys_type), intent(in) :: ls
!----------------------------------------------------
! Local variables:

integer :: isize, rsize, lsize, lssize, csize, p0size, p1size, p2size, &
           p3size, p4size, hsize
integer :: mem, i

!----------------------------------------------------
! Begin executable code

isize = bit_size(1)/8

! TEMP need a good portable way to determine the number of bytes in real.
! I think the following should always work, at least on the computers of 2006.
! The standard intrinsic function digits returns the number of binary digits in
! the real model for the type of real given.  The number of binary digits cannot
! possibly be bigger than the number of bits, so we look for the first bit
! size that can hold that number of binary digits.  I assume the number of
! bytes is a power of 2, it's at least 4 and no more than 64.

rsize = digits(1.0_my_real)
if (rsize <= 32) then
   rsize = 4
elseif (rsize <= 64) then
   rsize = 8
elseif (rsize <= 128) then
   rsize = 16
elseif (rsize <= 256) then
   rsize = 32
else
   rsize = 64
endif 

! TEMP assume a default logical is 4 bytes and a small logical is 1 byte

lsize = 4
lssize = 1 

! TEMP assume a default character is 1 byte

csize = 1

! The size of pointer descriptors varies among compilers.  There is also
! lost memory for small pointer arrays because of alignment requirements;
! this can be very large (lf95 can loose up to 192 bytes for every allocated
! pointer).  I do not include the memory lost to alignment in the memory
! count in this routine.  The following byte counts for pointer descriptors
! and alignment requirements were determined experimentally for some Fortran
! compilers on Linux.
!          Lahey   Intel   g95   Absoft   Nag
! Rank 0     8       4      4      24      4
! Rank 1    28      36     28      36     20
! Rank 2    44      48     40      48     32
! Rank 3    60      60     52      60     44
! Rank 4    76      72     64      72     56
! alignment 192     16     32      96      8

p0size = 8
p1size = 28
p2size = 44
p3size = 60
p4size = 76

! hash key size

hsize = KEY_SIZE*isize

mem = 0

! the matrix

mem = mem + p1size*11
if (associated(ls%column_index)) mem = mem + isize*size(ls%column_index)
if (associated(ls%begin_row)) mem = mem + isize*size(ls%begin_row)
if (associated(ls%end_row)) then
   if (.not.associated(ls%end_row,ls%end_row_linear) .and. &
       .not.associated(ls%end_row,ls%end_row_edge) .and. &
       .not.associated(ls%end_row,ls%end_row_face)) &
       mem = mem + isize*size(ls%end_row)
endif
if (associated(ls%end_row_linear)) mem = mem + isize*size(ls%end_row_linear)
if (associated(ls%end_row_edge)) mem = mem + isize*size(ls%end_row_edge)
if (associated(ls%end_row_face)) mem = mem + isize*size(ls%end_row_face)
if (associated(ls%matrix_val)) then
   if (.not. associated(ls%matrix_val,ls%stiffness) .and. &
       .not. associated(ls%matrix_val,ls%mass) .and. &
       .not. associated(ls%matrix_val,ls%shifted) .and. &
       .not. associated(ls%matrix_val,ls%condensed)) &
      mem = mem + rsize*size(ls%matrix_val)
endif
if (associated(ls%stiffness)) then
   if (.not. associated(ls%stiffness,ls%mass) .and. &
       .not. associated(ls%stiffness,ls%shifted) .and. &
       .not. associated(ls%stiffness,ls%condensed)) &
      mem = mem + rsize*size(ls%stiffness)
endif
if (associated(ls%mass)) then
   if (.not. associated(ls%mass,ls%shifted) .and. &
       .not. associated(ls%mass,ls%condensed)) &
      mem = mem + rsize*size(ls%mass)
endif
if (associated(ls%shifted)) then
   if (.not. associated(ls%shifted,ls%condensed)) &
      mem = mem + rsize*size(ls%shifted)
endif
if (associated(ls%condensed)) mem = mem + rsize*size(ls%condensed)

! begin_level

mem = mem + p1size
if (associated(ls%begin_level)) mem = mem + isize*size(ls%begin_level)

! rhs

mem = mem + p1size*3
if (associated(ls%rhs)) then
   if (.not. associated(ls%rhs,ls%rhs_nocond) .and. &
       .not. associated(ls%rhs,ls%rhs_cond)) mem = mem + rsize*size(ls%rhs)
endif
if (associated(ls%rhs_nocond)) then
   if (.not. associated(ls%rhs_nocond,ls%rhs_cond)) mem = mem + rsize*size(ls%rhs_nocond)
endif
if (associated(ls%rhs_cond)) mem = mem + rsize*size(ls%rhs_cond)

! solutions

mem = mem + p1size + p2size
if (associated(ls%solution)) mem = mem + rsize*size(ls%solution)
if (associated(ls%evecs)) mem = mem + rsize*size(ls%evecs)

! various things

mem = mem + p1size
if (associated(ls%equation_type)) mem = mem + isize*size(ls%equation_type)

mem = mem + p1size
if (associated(ls%gid)) mem = mem + hsize*size(ls%gid)

! TEMP don't know size of a hash table

mem = mem + p1size
if (associated(ls%iown)) mem = mem + lssize*size(ls%iown)

mem = mem + p1size
if (associated(ls%hold_dirich)) mem = mem + rsize*size(ls%hold_dirich)

mem = mem + p1size*3
if (associated(ls%r_mine)) mem = mem + rsize*size(ls%r_mine)
if (associated(ls%r_others)) mem = mem + rsize*size(ls%r_others)
if (associated(ls%need_r_others)) mem = mem + rsize*size(ls%need_r_others)

! blocks for condensed matrix

mem = mem + p1size*2
if (associated(ls%elem_block)) then
   mem = mem + (p1size + p2size + isize)*size(ls%elem_block)
   do i=1,size(ls%elem_block)
      if (associated(ls%elem_block(i)%matrix)) mem = mem + rsize*size(ls%elem_block(i)%matrix)
      if (associated(ls%elem_block(i)%ipiv)) mem = mem + isize*size(ls%elem_block(i)%ipiv)
   end do
endif

! some scalars

mem = mem + isize*7
mem = mem + rsize

! basis conversion matrix

mem = mem + p2size*4
if (associated(ls%s_int)) mem = mem + rsize*size(ls%s_int)
if (associated(ls%s_int_inv)) mem = mem + rsize*size(ls%s_int_inv)
if (associated(ls%s_bnd)) mem = mem + rsize*size(ls%s_bnd)
if (associated(ls%s_bnd_inv)) mem = mem + rsize*size(ls%s_bnd_inv)

! lapack matrices

mem = mem + lsize*3
mem = mem + p2size*2 + p1size*3 + isize*2
if (ls%coarse_band_exists) then
   if (associated(ls%coarse_matrix%matrix)) mem = mem + rsize*size(ls%coarse_matrix%matrix)
   if (associated(ls%coarse_matrix%rhs)) mem = mem + rsize*size(ls%coarse_matrix%rhs)
   if (associated(ls%coarse_matrix%renum)) mem = mem + isize*size(ls%coarse_matrix%renum)
   if (associated(ls%coarse_matrix%inv_renum)) mem = mem + isize*size(ls%coarse_matrix%inv_renum)
   if (associated(ls%coarse_matrix%ipiv)) mem = mem + isize*size(ls%coarse_matrix%ipiv)
endif
mem = mem + p2size*2 + p1size*3 + isize*2
if (ls%lapack_gen_band_exists .or. ls%lapack_symm_band_exists) then
   if (associated(ls%lapack_mat%matrix)) mem = mem + rsize*size(ls%lapack_mat%matrix)
   if (associated(ls%lapack_mat%rhs)) mem = mem + rsize*size(ls%lapack_mat%rhs)
   if (associated(ls%lapack_mat%renum)) mem = mem + isize*size(ls%lapack_mat%renum)
   if (associated(ls%lapack_mat%inv_renum)) mem = mem + isize*size(ls%lapack_mat%inv_renum)
   if (associated(ls%lapack_mat%ipiv)) mem = mem + isize*size(ls%lapack_mat%ipiv)
endif

! TEMP don't include PETSc, MUMPS, hypre and SuperLU matrices because I only
! have a handle for the bulk of the storage, which is managed inside the
! package, and could only count some incidental stuff in my structures

write(outunit,"(A,I12,A)") "Memory for linear system ",mem," bytes"

end subroutine count_memory

end module linear_system
