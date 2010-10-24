! **********************************************
! **********************************************
! **********************************************
! **********************************************
! This contains a workaround for the PGI generic interface bug
! **********************************************
! **********************************************
! **********************************************
! **********************************************






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

module hash_mod

!----------------------------------------------------
! This module contains the hash function for global indices.
! This version uses 2 integers in the hash key.
!----------------------------------------------------

!----------------------------------------------------
! Other modules used are:

use global

!----------------------------------------------------

implicit none
private
public hash_key, hash_table, & ! types
       NULL_KEY, KEY_SIZE, HASH_NOT_FOUND, LOG_MAX_KEY, & ! parameters
       hash_table_init, hash_table_destroy, & ! procedures
       hash_table_store, hash_table_restore, hash_table_copy, &
       hash_read_key, hash_overflow, &
       hash_insert, hash_remove, hash_decode_key, hash_print_key, &
       hash_pack_key, hash_unpack_key, hash_pack_key_pvm, hash_unpack_key_pvm, &
       mod,           & ! mod(key,integer)
       assignment(=), & ! key = integer
       operator(+),   & ! key + integer
       operator(-),   & ! key - integer
       operator(*),   & ! integer * key
       operator(/),   & ! key / integer
       operator(==),  & ! key == key, key == integer
       operator(<)      ! key < key

!----------------------------------------------------
! The following types are defined:

type hash_key
   private
   integer :: key,key2
end type

type hash_entry
   type(hash_key) :: key
   integer :: index
   type(hash_entry), pointer :: next
end type hash_entry

type entry_ptr
   type(hash_entry), pointer :: ptr
end type entry_ptr

type hash_table
   private
   integer :: size
   type(entry_ptr), pointer :: table(:)
end type hash_table

!----------------------------------------------------
! The following parameters are defined:

type(hash_key), parameter :: NULL_KEY = hash_key(-1,0)
integer, parameter :: KEY_SIZE = 2 ! integer words
integer, parameter :: HASH_NOT_FOUND = -20
integer, parameter :: LOG_MAX_KEY = 2*(bit_size(0) - 2)

!----------------------------------------------------
! The following variables are defined:

!----------------------------------------------------
! The public module procedures are:

!   subroutine hash_table_init(table,table_size)
! This routine initializes a hash table with size table_size
!   type(hash_table), intent(out) :: table
!   type(table_size), intent(in) :: table_size
!   end subroutine hash_table_init

!   subroutine hash_table_destroy(table)
! This routine destroys a hash table and frees the associated memory
!   type(hash_table), intent(inout) :: table
!   end subroutine hash_table_destroy

!   subroutine hash_insert(key,ind,table)
! This routine inserts (key,ind) in the hash table
!   type(hash_key), intent(in) :: key
!   integer, intent(in) :: ind
!   type(hash_table), intent(inout) :: table
!   end function hash_insert

!   subroutine hash_remove(key,table)
! This routine removes key from the hash table
!   type(hash_key), intent(in) :: key
!   type(hash_table), intent(inout) :: table
!   end subroutine hash_remove

!   function hash_decode_key(key,table)
! This routine returns the index associated with key
!   integer :: hash_decode_key
!   type(hash_key), intent(in) :: key
!   type(hash_table), intent(in) :: table
!   end function hash_decode_key

!   function hash_overflow(key,m,a)
! This routine returns true if m*key+a will overflow the hash values
!   logical :: hash_overflow
!   type(hash_key), intent(in) :: key
!   integer, intent(in) :: m, a
!   end function hash_overflow

!   subroutine hash_print_key(key,iounit,str)
! This routine prints a key
!   type(hash_key), intent(in) :: key
!OR type(hash_key), intent(in) :: key(:)
!   integer, intent(in) :: iounit
!   character(len=*), optional, intent(in) :: str
!   end subroutine hash_print_key

!   subroutine hash_read_key(key,iounit)
! This routine reads keys printed by hash_print_key with the default string
!   type(hash_key), intent(out) :: key(:)
!   integer, intent(in) :: iounit
!   end subroutine hash_read_key

!   subroutine hash_table_store(table,iounit)
! This routine writes the hash table to unit iounit.
!   type(hash_table), intent(in) :: table
!   integer, intent(in) :: iounit
!   end subroutine hash_table_store

!   subroutine hash_table_restore(table,iounit)
! This routine reads the hash table from unit iounit and places it in
! table.  Since there is no good way in Fortran 90 to check if table contains
! data, this creates a memory leak if the table is not currently empty.
!   type(hash_table), intent(in) :: table
!   integer, intent(in) :: iounit
!   end subroutine hash_table_restore

!   subroutine hash_table_copy(old_table,new_table)
! This routine copies old_table to new_table
!   type(hash_table), intent(in) :: old_table
!   type(hash_table), intent(out) :: new_table
!   end subroutine hash_table_copy

!   subroutine hash_pack_key(key,array,ind)
! This routine packs key into array starting at ind
!   type(hash_key), intent(in) :: key
!OR type(hash_key), intent(in) :: key(:)
!   integer, intent(inout) :: array(:)
!   integer, intent(in) :: ind
!   end subroutine hash_pack_key

!   function hash_unpack_key(array,ind)
! This routine unpacks the key in array starting at ind
!   type(hash_key) :: hash_unpack_key
!   integer, intent(in) :: array(:), ind
!   end function hash_unpack_key

!   subroutine hash_pack_key_pvm(keys,itype_pvm)
! This routine packs an array of keys into a PVM message buffer
!   type(hash_key), intent(in) :: keys(:)
!   integer, intent(in) :: itype_pvm
!   end subroutine hash_pack_key_pvm

!   subroutine hash_unpack_key_pvm(keys,nkey,itype_pvm)
! This routine unpacks an array of keys from a PVM message buffer
!   type(hash_key), intent(inout) :: keys(:)
!   integer, intent(in) :: nkey, itype_pvm
!   end subroutine hash_unpack_key_pvm

!----------------------------------------------------
! The following interface operators are defined:

interface assignment(=)
   module procedure key_equals_integer
end interface

interface operator(*)
   module procedure integer_times_key
end interface

interface operator(+)
   module procedure key_plus_integer
end interface

interface operator(-)
   module procedure key_minus_integer
end interface

interface operator(/)
   module procedure key_div_integer
end interface

interface operator(==)
   module procedure key_compare, key_integer_compare
end interface

interface operator(<)
   module procedure key_less_than
end interface

!----------------------------------------------------
! The following generic interfaces are defined:

interface hash_table_init
   module procedure hash_table_init_loc
end interface

interface hash_table_destroy
   module procedure hash_table_destroy_loc
end interface

interface hash_insert
   module procedure hash_insert_loc
end interface

interface hash_remove
   module procedure hash_remove_loc
end interface

interface hash_decode_key
   module procedure hash_decode_key_loc
end interface

interface hash_overflow
   module procedure hash_overflow_loc
end interface

interface mod
   module procedure hash_modfunc
end interface

interface hash_pack_key
   module procedure hash_pack_key_scalar, hash_pack_key_array
end interface

interface hash_unpack_key
   module procedure hash_unpack_key_loc
end interface

interface hash_pack_key_pvm
   module procedure hash_pack_key_pvm_loc
end interface

interface hash_unpack_key_pvm
   module procedure hash_unpack_key_pvm_loc
end interface

interface hash_print_key
   module procedure hash_print_key_scalar, hash_print_key_array
end interface

!----------------------------------------------------

contains

!---------------------------------------------------------------------
!  INITIALIZATION AND FINALIZATION ROUTINES
!---------------------------------------------------------------------

!          ---------------
subroutine hash_table_init_loc(table,table_size)
!          ---------------

!----------------------------------------------------
! This routine initializes a hash table with size table_size
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_table), intent(out) :: table
integer, intent(in) :: table_size

!----------------------------------------------------
! Local variables:

integer :: i, astat

!----------------------------------------------------
! Begin executable code

! set the size, allocate the array and nullify the pointers

table%size = table_size
allocate(table%table(table_size),stat=astat)
if (astat /= 0) then
   call fatal("memory allocation failed in hash_table_init")
   nullify(table%table)
   table%size = 0
   ierr = ALLOC_FAILED
   return
endif
do i=1,table_size
   nullify(table%table(i)%ptr)
end do

return
end subroutine hash_table_init_loc

!          ------------------
subroutine hash_table_destroy_loc(table)
!          ------------------

!----------------------------------------------------
! This routine destroys a hash table and frees the associated memory
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_table), intent(inout) :: table

!----------------------------------------------------
! Local variables:

integer :: i, dstat
type(hash_entry), pointer :: current, next

!----------------------------------------------------
! Begin executable code

! deallocate all the table entries

do i=1,table%size
   next => table%table(i)%ptr
   do while (associated(next))
      current => next
      next => current%next
      deallocate(current,stat=dstat)
   end do
end do

! deallocate the table

if (associated(table%table)) deallocate(table%table,stat=dstat)
table%size = 0

return
end subroutine hash_table_destroy_loc

!---------------------------------------------------------------------
!  HASH TABLE MANIPULATIONS
!---------------------------------------------------------------------

!          -----------
subroutine hash_insert_loc(key,ind,table)
!          -----------

!----------------------------------------------------
! This routine inserts (key,ind) in the hash table
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_key), intent(in) :: key
integer, intent(in) :: ind
type(hash_table), intent(inout) :: table

!----------------------------------------------------
! Local variables:

integer :: table_sub, astat
type(hash_entry), pointer :: newentry, next_entry
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

! compute the table subscript corresponding to key

table_sub = hash_function(key,table%size)

! if the entry is already in the table, don't need to add it

next_entry => table%table(table_sub)%ptr
do
   if (.not. associated(next_entry)) exit
   if (next_entry%key == key) return
   next_entry => next_entry%next
end do

! insert in the hash table, at the head of linked list indexed by table_sub

allocate(newentry,stat=astat)
if (astat /= 0) then
   call fatal("memory allocation failed in hash_insert")
   ierr = ALLOC_FAILED
   return
endif
newentry%key = key
newentry%index = ind
newentry%next => table%table(table_sub)%ptr
table%table(table_sub)%ptr => newentry

return
end subroutine hash_insert_loc

!          -----------
subroutine hash_remove_loc(key,table)
!          -----------

!----------------------------------------------------
! This routine removes key from the hash table
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_key), intent(in) :: key
type(hash_table), intent(inout) :: table

!----------------------------------------------------
! Local variables:

integer :: table_sub, dstat
type(hash_entry), pointer :: next_entry, previous

!----------------------------------------------------
! Begin executable code

! compute the table subscript corresponding to key

table_sub = hash_function(key,table%size)

! find the key in the linked list at table_sub

next_entry => table%table(table_sub)%ptr
nullify(previous)
do
   if (.not.associated(next_entry)) exit
   if (next_entry%key == key) exit
   previous => next_entry
   next_entry => next_entry%next
end do
if (.not.associated(next_entry)) then
   call warning("hash key not matched in hash_remove")
   return
endif

! remove it from the linked list

if (.not.associated(previous)) then

! it is the head of the list

   table%table(table_sub)%ptr => next_entry%next

else

! it is not the head of the list

   previous%next => next_entry%next

endif

deallocate(next_entry,stat=dstat)

return
end subroutine hash_remove_loc

!        ---------------
function hash_decode_key_loc(key,table)
!        ---------------

!----------------------------------------------------
! This routine returns the index associated with key
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer :: hash_decode_key_loc
type(hash_key), intent(in) :: key
type(hash_table), intent(in) :: table

!----------------------------------------------------
! Local variables:

integer :: table_sub
type(hash_entry), pointer :: next_entry

!----------------------------------------------------
! Begin executable code

! compute the table subscript corresponding to key

table_sub = hash_function(key,table%size)

! find the key in the linked list at table_sub

next_entry => table%table(table_sub)%ptr
do
   if (.not.associated(next_entry)) exit
   if (next_entry%key == key) exit
   next_entry => next_entry%next
end do

! return the index

if (.not.associated(next_entry)) then
   hash_decode_key_loc = HASH_NOT_FOUND
else
   hash_decode_key_loc = next_entry%index
endif

return
end function hash_decode_key_loc

!---------------------------------------------------------------------
!  PROCEDURES FOR OPERATOR INTERFACES
!---------------------------------------------------------------------

!          ------------------
subroutine key_equals_integer(key,int)
!          ------------------

!----------------------------------------------------
! This routine assigns integer int to a key
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_key), intent(out) :: key
integer, intent(in) :: int

!----------------------------------------------------
! Begin executable code

key%key = int
key%key2 = 0

return
end subroutine key_equals_integer

!                  -----------------
recursive function integer_times_key(int,key) result(res)
!                  -----------------

!----------------------------------------------------
! This routine multiplies an integer and a key
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_key) :: res
integer, intent(in) :: int
type(hash_key), intent(in) :: key

!----------------------------------------------------
! Begin executable code

! if int is 1, nothing to do

if (int == 1) then
   res = key
   return
endif

! lazy way of multiplying by 4

if (int == 4) then
   res = 2*(2*key)
   return
endif

if (int/=2) then
   call fatal("can only multiply hash key by 1, 2 or 4; requested",intlist=(/int/))
   stop
endif

if (key%key2 > huge(key%key2)/int) then
   call fatal("overflow in hash key.  int and key are",intlist=(/int/))
   call hash_print_key(key,errunit)
endif

res%key = int*key%key
res%key2 = int*key%key2

if (res%key > huge(key%key)/2) then
   res%key = res%key - huge(key%key)/2 + 1
   res%key2 = res%key2 + 1
endif

return
end function integer_times_key

!        ----------------
function key_plus_integer(key,int)
!        ----------------

!----------------------------------------------------
! This routine adds an integer to a key
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_key) :: key_plus_integer
type(hash_key), intent(in) :: key
integer, intent(in) :: int

!----------------------------------------------------
! Begin executable code

key_plus_integer%key = int+key%key
key_plus_integer%key2 = key%key2

if (key_plus_integer%key > huge(key%key)/2) then
   key_plus_integer%key = key_plus_integer%key - huge(key%key)/2 + 1
   key_plus_integer%key2 = key_plus_integer%key2 + 1
endif
 
return
end function key_plus_integer

!        -----------------
function key_minus_integer(key,int)
!        -----------------

!----------------------------------------------------
! This routine subtracts an integer from a key
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_key) :: key_minus_integer
type(hash_key), intent(in) :: key
integer, intent(in) :: int

!----------------------------------------------------
! Begin executable code

key_minus_integer%key = key%key - int
key_minus_integer%key2 = key%key2
if (key_minus_integer%key < 0) then
   key_minus_integer%key = key_minus_integer%key + huge(key%key)/2 - 1
   key_minus_integer%key2 = key_minus_integer%key2 - 1
endif

return
end function key_minus_integer

!                  ---------------
recursive function key_div_integer(key,int) result(res)
!                  ---------------

!----------------------------------------------------
! This routine divides a key by an integer
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_key) :: res
type(hash_key), intent(in) :: key
integer, intent(in) :: int

!----------------------------------------------------
! Begin executable code

if (int==4) then
   res = (key/2)/2
   return
endif

if (int/=2) then
   call fatal("can only divide hash key by 2 or 4; requested",intlist=(/int/))
   stop
endif

res%key = key%key/int
res%key2 = key%key2/int

if (res%key2*int /= key%key2) then
   res%key = res%key + huge(key%key)/4
endif
 
return
end function key_div_integer

!        -----------
function key_compare(key1,key2)
!        -----------

!----------------------------------------------------
! This routine determines if two keys are equal
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

logical :: key_compare
type(hash_key), intent(in) :: key1, key2

!----------------------------------------------------
! Begin executable code

key_compare = (key1%key==key2%key .and. key1%key2==key2%key2)

return
end function key_compare

!        -------------------
function key_integer_compare(key,int)
!        -------------------

!----------------------------------------------------
! This routine determines if a key is equal to an integer
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

logical :: key_integer_compare
type(hash_key), intent(in) :: key
integer, intent(in) :: int

!----------------------------------------------------
! Begin executable code

key_integer_compare = (key%key==int .and. key%key2==0)

return
end function key_integer_compare

!        -------------
function key_less_than(key1,key2)
!        -------------

!----------------------------------------------------
! This routine determines if key1 is less than key2
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

logical :: key_less_than
type(hash_key), intent(in) :: key1, key2

!----------------------------------------------------
! Begin executable code

if (key1%key2 == key2%key2) then
   key_less_than = key1%key < key2%key
else
   key_less_than = key1%key2 < key2%key2
endif

return
end function key_less_than

!---------------------------------------------------------------------
!  I/O
!---------------------------------------------------------------------

!          ---------------------
subroutine hash_print_key_scalar(key,iounit,str)
!          ---------------------

!----------------------------------------------------
! This routine prints a key
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_key), intent(in) :: key
integer, intent(in) :: iounit
character(len=*), optional, intent(in) :: str
character(len=11) :: form

!----------------------------------------------------
! Begin executable code

inquire(iounit,form=form)

select case(form)

case("FORMATTED")
   if (present(str)) then
      write(iounit,"(A,2I11)") str,key
   else
      write(iounit,"(A,2I11)") "Hash key:  ",key
   endif

case("UNFORMATTED")
   write(iounit) key

end select

return
end subroutine hash_print_key_scalar

!          --------------------
subroutine hash_print_key_array(key,iounit,str)
!          --------------------

!----------------------------------------------------
! This routine prints an array of keys
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_key), intent(in) :: key(:)
integer, intent(in) :: iounit
character(len=*), optional, intent(in) :: str
character(len=11) :: form

!----------------------------------------------------
! Begin executable code

inquire(iounit,form=form)

select case(form)

case("FORMATTED")
   if (present(str)) then
      write(iounit,"(A,200I11)") str,key
   else
      write(iounit,"(A,200I11)") "Hash key:  ",key(:)
   endif

case("UNFORMATTED")
   write(iounit) key(:)

end select

return
end subroutine hash_print_key_array

!          -------------
subroutine hash_read_key(key,iounit)
!          -------------

!----------------------------------------------------
! This routine reads keys printed by hash_print_key
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_key), intent(out) :: key(:)
integer, intent(in) :: iounit

!----------------------------------------------------
! local variables

character(len=10) :: str1, str2
character(len=11) :: form
!----------------------------------------------------
! Begin executable code

inquire(iounit,form=form)

select case(form)

case ("FORMATTED")
   read(iounit,*) str1,str2,key(:)

case ("UNFORMATTED")
   read(iounit) key(:)

end select

return
end subroutine hash_read_key

!          ----------------
subroutine hash_table_store(table,iounit)
!          ----------------

!----------------------------------------------------
! This routine writes the hash table to unit iounit.  It can be read
! back with hash_table_restore
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_table), intent(in) :: table
integer, intent(in) :: iounit
!----------------------------------------------------
! Local variables:

integer :: i
type(hash_entry), pointer :: entry
character(len=11) :: form
!----------------------------------------------------
! Begin executable code

inquire(iounit,form=form)

select case(form)

case ("UNDEFINED")
   call fatal("File appears unopened in hash_table_restore",intlist=(/iounit/))
   ierr = USER_INPUT_ERROR
   return

case ("FORMATTED")

   write(iounit,*) 3 ! version number
   write(iounit,*) KEY_SIZE

   write(iounit,*) table%size

   do i=1,table%size
      entry => table%table(i)%ptr
      do while (associated(entry))
         write(iounit,*) entry%key, entry%index
         entry => entry%next
      end do
   end do
   write(iounit,*) NULL_KEY, -1

case ("UNFORMATTED")

   write(iounit) 3 ! version number
   write(iounit) KEY_SIZE

   write(iounit) table%size

   do i=1,table%size
      entry => table%table(i)%ptr
      do while (associated(entry))
         write(iounit) entry%key, entry%index
         entry => entry%next
      end do
   end do
   write(iounit) NULL_KEY, -1

end select

return
end subroutine hash_table_store

!          ------------------
subroutine hash_table_restore(table,iounit)
!          ------------------

!----------------------------------------------------
! This routine reads the hash table from unit iounit and places it in
! table.  Since there is no good way in Fortran 90 to check if table contains
! data, this creates a memory leak if the table is not currently empty.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_table), intent(out) :: table
integer, intent(in) :: iounit
!----------------------------------------------------
! Local variables:

integer :: version, saved_key_size, table_size, index
type(hash_key) :: key
character(len=11) :: form
!----------------------------------------------------
! Begin executable code

inquire(iounit,form=form)

select case(form)

case ("UNDEFINED")
   call fatal("File appears unopened in hash_table_restore",intlist=(/iounit/))
   ierr = USER_INPUT_ERROR
   return

case ("FORMATTED")

   read(iounit,*) version

   select case(version)

   case(3)

      read(iounit,*) saved_key_size
      if (saved_key_size /= KEY_SIZE) then
         call fatal("hash_table_restore: different key size", &
                      intlist=(/saved_key_size,KEY_SIZE/))
         stop
      endif

      read(iounit,*) table_size
      call hash_table_init(table,table_size)

      do
         read(iounit,*) key, index
         if (index < 0) exit
         call hash_insert(key,index,table)
      end do

   case default

      call fatal("Unknown version number in hash_table_restore")
      stop

   end select

case ("UNFORMATTED")

   read(iounit) version

   select case(version)

   case(3)

      read(iounit) saved_key_size
      if (saved_key_size /= KEY_SIZE) then
         call fatal("hash_table_restore: different key size", &
                      intlist=(/saved_key_size,KEY_SIZE/))
         stop
      endif

      read(iounit) table_size
      call hash_table_init(table,table_size)

      do
         read(iounit) key, index
         if (index < 0) exit
         call hash_insert(key,index,table)
      end do

   case default

      call fatal("Unknown version number in hash_table_restore")
      stop

   end select

end select

return
end subroutine hash_table_restore

!          ---------------
subroutine hash_table_copy(old_table,new_table)
!          ---------------

!----------------------------------------------------
! This routine copies old_table to new_table
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_table), intent(in) :: old_table
type(hash_table), intent(out) :: new_table
!----------------------------------------------------
! Local variables:

integer :: i
type(hash_entry), pointer :: old_entry
!----------------------------------------------------
! Begin executable code

call hash_table_init(new_table,old_table%size)
do i=1,old_table%size
   old_entry => old_table%table(i)%ptr
   do while (associated(old_entry))
      call hash_insert(old_entry%key,old_entry%index,new_table)
      old_entry => old_entry%next
   end do
end do

end subroutine hash_table_copy

!---------------------------------------------------------------------
!  OTHER
!---------------------------------------------------------------------

!        -------------
function hash_overflow_loc(key,m,a)
!        -------------

!----------------------------------------------------
! This routine returns true if m*key+a will overflow the hash values
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

logical :: hash_overflow_loc
type(hash_key), intent(in) :: key
integer, intent(in) :: m, a
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

hash_overflow_loc = key%key2 > (huge(key%key2)-a)/m

end function hash_overflow_loc

!        -------------
function hash_function(key,table_size)
!        -------------

!----------------------------------------------------
! This routine evaluates the hash function for the given key.
! Based on Robert Jenkin's 32-bit mixing function, followed by a
! mod to place it in the range 1:table_size.  It only hashes on the
! low order word of the key.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer :: hash_function
type(hash_key), intent(in) :: key
integer, intent(in) :: table_size

!----------------------------------------------------
! Local variables:

integer :: lkey

!----------------------------------------------------
! Begin executable code

lkey = key%key

lkey = lkey + ishft(lkey,12)
lkey = ieor(lkey,ishft(lkey,-22))
lkey = lkey + ishft(lkey,4)
lkey = ieor(lkey,ishft(lkey,-9))
lkey = lkey + ishft(lkey,10)
lkey = ieor(lkey,ishft(lkey,-2))
lkey = lkey + ishft(lkey,7)
lkey = ieor(lkey,ishft(lkey,-12))

hash_function = 1 + mod(abs(lkey),table_size)

return
end function hash_function

!          --------------------
subroutine hash_pack_key_scalar(key,array,ind)
!          --------------------

!----------------------------------------------------
! This routine packs a single key into array starting at ind
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_key), intent(in) :: key
integer, intent(inout) :: array(:)
integer, intent(in) :: ind

!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

array(ind) = key%key
array(ind+1) = key%key2

return
end subroutine hash_pack_key_scalar

!          -------------------
subroutine hash_pack_key_array(key,array,ind)
!          -------------------

!----------------------------------------------------
! This routine packs an array of keys into array starting at ind
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_key), intent(in) :: key(:)
integer, intent(inout) :: array(:)
integer, intent(in) :: ind

!----------------------------------------------------
! Local variables:

integer :: i
!----------------------------------------------------
! Begin executable code

do i=0,size(key)-1
   array(ind+2*i  ) = key(i+1)%key
   array(ind+2*i+1) = key(i+1)%key2
end do

return
end subroutine hash_pack_key_array

!        ---------------
function hash_unpack_key_loc(array,ind)
!        ---------------

!----------------------------------------------------
! This routine unpacks the key in array starting at ind
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_key) :: hash_unpack_key_loc
integer, intent(in) :: array(:), ind

!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

hash_unpack_key_loc%key = array(ind)
hash_unpack_key_loc%key2 = array(ind+1)

return
end function hash_unpack_key_loc

!          -----------------
subroutine hash_pack_key_pvm_loc(keys,itype_pvm)
!          -----------------

!----------------------------------------------------
! This routine packs an array of keys into a PVM message buffer.  This is
! used so we don't have to copy the keys into an integer array which then
! gets packed into the message.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_key), intent(in) :: keys(:)
integer, intent(in) :: itype_pvm

!----------------------------------------------------
! Local variables:

integer :: info, i
!----------------------------------------------------
! Begin executable code

do i=1,size(keys)
   call pvmfpack(itype_pvm,(/keys(i)%key/),1,1,info)
   call pvmfpack(itype_pvm,(/keys(i)%key2/),1,1,info)
end do

return
end subroutine hash_pack_key_pvm_loc

!          -------------------
subroutine hash_unpack_key_pvm_loc(keys,nkey,itype_pvm)
!          -------------------

!----------------------------------------------------
! This routine unpacks an array of keys from a PVM message buffer.  This is
! used so we don't have to copy the keys from an integer array which was
! unpacked from the message.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_key), intent(inout) :: keys(:)
integer, intent(in) :: nkey, itype_pvm

!----------------------------------------------------
! Local variables:

integer :: info, i, unpacked(2)
!----------------------------------------------------
! Begin executable code

do i=1,nkey
   call pvmfunpack(itype_pvm,unpacked,2,1,info)
   keys(i)%key = unpacked(1)
   keys(i)%key2= unpacked(2)
end do

end subroutine hash_unpack_key_pvm_loc

!        ------------
function hash_modfunc(key,int)
!        ------------

!----------------------------------------------------
! This routine extends the generic mod function to hash keys
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer :: hash_modfunc
type(hash_key), intent(in) :: key
integer, intent(in) :: int

!----------------------------------------------------
! Begin executable code

hash_modfunc = mod(key%key,int)

end function hash_modfunc

! hash has its own version of warning and fatal because they are in
! module messpass, which uses hash_mod

!          -------
subroutine warning(msg,msg2,intlist,reallist)
!          -------

!----------------------------------------------------
! This routine handles warning messages
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

character(len=*), intent(in) :: msg
character(len=*), intent(in), optional :: msg2
integer, intent(in), optional, dimension(:) :: intlist
real(my_real), intent(in), optional, dimension(:) :: reallist
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

if (warn) then
   write(errunit,"(A)")
   write(errunit,"(A)") "------------------------------------------------------"
   write(errunit,"(3A)") "          PHAML Version ",version_number," WARNING"
   write(errunit,"(A)") msg
   if (present(msg2)) write(errunit,"(A)") msg2
   if (present(intlist)) write(errunit,"(7I11)") intlist
   if (present(reallist)) write(errunit,"(SS,1P,4E18.10E2)") reallist
   write(errunit,"(A)") "------------------------------------------------------"
   write(errunit,"(A)")
endif

end subroutine warning

!          -----
subroutine fatal(msg,msg2,intlist,reallist)
!          -----

!----------------------------------------------------
! This routine handles fatal errors
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

character(len=*), intent(in) :: msg
character(len=*), intent(in), optional :: msg2
integer, intent(in), optional, dimension(:) :: intlist
real(my_real), intent(in), optional, dimension(:) :: reallist
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

write(errunit,"(A)")
write(errunit,"(A)") "------------------------------------------------------"
write(errunit,"(3A)") "          PHAML Version ",version_number," ERROR"
write(errunit,"(A)") msg
if (present(msg2)) write(errunit,"(A)") msg2
if (present(intlist)) write(errunit,"(7I11)") intlist
if (present(reallist)) write(errunit,"(SS,1P,4E18.10E2)") reallist
write(errunit,"(A)") "------------------------------------------------------"
write(errunit,"(A)")

end subroutine fatal

end module hash_mod
