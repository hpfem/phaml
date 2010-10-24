! This stub library is based on Appendices B and D of
! the OpenMP Specifications 3.0

module omp_lib

integer, parameter :: openmp_version = 200805 ! OpenMP Fortran API v3.0
integer, parameter :: omp_integer_kind = kind(0)
integer, parameter :: omp_logical_kind = kind(.false.)
integer, parameter :: omp_lock_kind = kind(0) 
integer, parameter :: omp_nest_lock_kind = kind(0)
integer, parameter :: omp_sched_kind = kind(0)
integer(kind=omp_sched_kind), parameter :: omp_sched_static = 1
integer(kind=omp_sched_kind), parameter :: omp_sched_dynamic = 2
integer(kind=omp_sched_kind), parameter :: omp_sched_guided = 3
integer(kind=omp_sched_kind), parameter :: omp_sched_auto = 4 

contains

subroutine omp_set_num_threads (number_of_threads_expr)
integer (kind=omp_integer_kind), intent(in) :: number_of_threads_expr 
end subroutine omp_set_num_threads

function omp_get_num_threads ()
integer (kind=omp_integer_kind) :: omp_get_num_threads
omp_get_num_threads = 1
end function omp_get_num_threads

function omp_get_max_threads ()
integer (kind=omp_integer_kind) :: omp_get_max_threads
omp_get_max_threads = 1
end function omp_get_max_threads

function omp_get_thread_num () 
integer (kind=omp_integer_kind) :: omp_get_thread_num
omp_get_thread_num = 0
end function omp_get_thread_num 

function omp_get_num_procs ()
integer (kind=omp_integer_kind) :: omp_get_num_procs
omp_get_num_procs = 1
end function omp_get_num_procs

function omp_in_parallel ()
logical (kind=omp_logical_kind) :: omp_in_parallel
omp_in_parallel = .false.
end function omp_in_parallel

subroutine omp_set_dynamic (enable_expr)
logical (kind=omp_logical_kind), intent(in) :: enable_expr 
end subroutine omp_set_dynamic

function omp_get_dynamic ()
logical (kind=omp_logical_kind) :: omp_get_dynamic
omp_get_dynamic = .false.
end function omp_get_dynamic

subroutine omp_set_nested (enable_expr)
logical (kind=omp_logical_kind), intent(in) :: enable_expr 
end subroutine omp_set_nested

function omp_get_nested ()
logical (kind=omp_logical_kind) :: omp_get_nested
omp_get_nested = .false.
end function omp_get_nested

subroutine omp_set_schedule (kind, modifier)
integer(kind=omp_sched_kind), intent(in) :: kind
integer(kind=omp_integer_kind), intent(in) :: modifier
end subroutine omp_set_schedule

subroutine omp_get_schedule (kind, modifier)
integer(kind=omp_sched_kind), intent(out) :: kind
integer(kind=omp_integer_kind), intent(out)::modifier
kind = omp_sched_static
modifier = 0
end subroutine omp_get_schedule

function omp_get_thread_limit()
integer (kind=omp_integer_kind) :: omp_get_thread_limit
omp_get_thread_limit = 1
end function omp_get_thread_limit

subroutine omp_set_max_active_levels(var)
integer (kind=omp_integer_kind), intent(in) :: var
end subroutine omp_set_max_active_levels 

function omp_get_max_active_levels()
integer (kind=omp_integer_kind) :: omp_get_max_active_levels 
omp_get_max_active_levels = 0
end function omp_get_max_active_levels

function omp_get_level()
integer (kind=omp_integer_kind) :: omp_get_level
omp_get_level = 0
end function omp_get_level

function omp_get_ancestor_thread_num(level)
integer (kind=omp_integer_kind), intent(in) :: level
integer (kind=omp_integer_kind) :: omp_get_ancestor_thread_num
if ( level .eq. 0 ) then
   omp_get_ancestor_thread_num = 0
else
   omp_get_ancestor_thread_num = -1
end if
end function omp_get_ancestor_thread_num

function omp_get_team_size(level)
integer (kind=omp_integer_kind), intent(in) :: level
integer (kind=omp_integer_kind) :: omp_get_team_size
if ( level .eq. 0 ) then
   omp_get_team_size = 1
else
   omp_get_team_size = -1
end if
end function omp_get_team_size

function omp_get_active_level()
integer (kind=omp_integer_kind) :: omp_get_active_level
omp_get_active_level = 0
end function omp_get_active_level

subroutine omp_init_lock (lock)
integer (kind=omp_lock_kind), intent(out) :: lock
! lock is
!  0 if the simple lock is not initialized
! -1 if the simple lock is initialized but not set
!  1 if the simple lock is set
lock = -1
end subroutine omp_init_lock

subroutine omp_destroy_lock (lock)
integer (kind=omp_lock_kind), intent(inout) :: lock 
lock = 0
end subroutine omp_destroy_lock

subroutine omp_set_lock (lock)
integer (kind=omp_lock_kind), intent(inout) :: lock
if (lock .eq. -1) then
   lock = 1
elseif (lock .eq. 1) then
   print *, 'error: deadlock in using lock variable'
   stop
else
   print *, 'error: lock not initialized'
   stop
endif
end subroutine omp_set_lock

subroutine omp_unset_lock (lock)
integer (kind=omp_lock_kind), intent(inout) :: lock
if (lock .eq. 1) then
   lock = -1
elseif (lock .eq. -1) then
   print *, 'error: lock not set'
   stop
else
   print *, 'error: lock not initialized'
   stop
endif
end subroutine omp_unset_lock 

function omp_test_lock (lock)
logical (kind=omp_logical_kind) :: omp_test_lock
integer (kind=omp_lock_kind), intent(inout) :: lock
if (lock .eq. -1) then
   lock = 1
   omp_test_lock = .true.
elseif (lock .eq. 1) then
   omp_test_lock = .false.
else
   print *, 'error: lock not initialized'
stop
endif
end function omp_test_lock

subroutine omp_init_nest_lock (nlock)
! nlock is
!  0 if the nestable lock is not initialized
! -1 if the nestable lock is initialized but not set
!  1 if the nestable lock is set
! no use count is maintained
integer (kind=omp_nest_lock_kind), intent(out) :: nlock
nlock = -1
end subroutine omp_init_nest_lock

subroutine omp_destroy_nest_lock (nlock)
integer (kind=omp_nest_lock_kind), intent(inout) :: nlock
nlock = 0
end subroutine omp_destroy_nest_lock

subroutine omp_set_nest_lock (nlock)
integer (kind=omp_nest_lock_kind), intent(inout) :: nlock
if (nlock .eq. -1) then
   nlock = 1
elseif (nlock .eq. 0) then
   print *, 'error: nested lock not initialized'
   stop
else
   print *, 'error: deadlock using nested lock variable'
   stop
endif
end subroutine omp_set_nest_lock

subroutine omp_unset_nest_lock (nlock)
integer (kind=omp_nest_lock_kind), intent(inout) :: nlock
if (nlock .eq. -1) then
   nlock = 1
elseif (nlock .eq. 0) then
   print *, 'error: nested lock not initialized'
   stop
else
   print *, 'error: deadlock using nested lock variable'
   stop
endif
end subroutine omp_unset_nest_lock

function omp_test_nest_lock (nlock)
integer (kind=omp_integer_kind) :: omp_test_nest_lock
integer (kind=omp_nest_lock_kind), intent(inout) :: nlock
if (nlock .eq. -1) then
   nlock = 1
   omp_test_nest_lock = 1
elseif (nlock .eq. 1) then
   omp_test_nest_lock = 0
else
   print *, 'error: nested lock not initialized'
   stop
endif
end function omp_test_nest_lock

function omp_get_wtick ()
double precision :: omp_get_wtick
integer :: rate
call system_clock(count_rate=rate)
if (rate == 0) then
   omp_get_wtick = 1
else
   omp_get_wtick = 1.0d0/rate
endif
end function omp_get_wtick

function omp_get_wtime ()
double precision :: omp_get_wtime
integer :: count, rate
call system_clock(count=count,count_rate=rate)
if (rate == 0) then
   omp_get_wtime = 0.0d0
else
   omp_get_wtime = real(count,kind(0.0d0))/rate
endif
end function omp_get_wtime

end module omp_lib 
