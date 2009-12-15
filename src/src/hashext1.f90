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

module hash_eq_mod

!----------------------------------------------------
! This module contains the hash function for global equation indices.
! This version uses one integer in the hash key and one integer to indicate
! whether the equation comes from a vertex, edge or face, and which rank
! it is on a face or edge.
!----------------------------------------------------

!----------------------------------------------------
! Other modules used are:

use global
use hash_mod

!----------------------------------------------------

implicit none
private
public hash_key_eq, hash_table_eq, & ! types
       NULL_KEY_EQ, KEY_SIZE, HASH_NOT_FOUND, LOG_MAX_KEY, & ! parameters
       hash_table_init, hash_table_destroy, & ! procedures
       hash_overflow, &
       hash_insert, hash_remove, hash_decode_key, hash_print_key, &
       hash_pack_key, hash_unpack_key, hash_pack_key_pvm, hash_unpack_key_pvm, &
       mod,           & ! mod(key,integer)
       assignment(=), & ! key = integer and rank
       operator(+),   & ! key + integer
       operator(-),   & ! key - integer
       operator(*),   & ! integer * key
       operator(/),   & ! key / integer
       operator(==),  & ! key == key, key == integer
       operator(<)      ! key < key

!----------------------------------------------------
! The following types are defined:

type hash_key_eq
   private
   integer :: key,rank
end type

type hash_entry_eq
   type(hash_key_eq) :: key
   integer :: index
   type(hash_entry_eq), pointer :: next
end type hash_entry_eq

type entry_ptr_eq
   type(hash_entry_eq), pointer :: ptr
end type entry_ptr_eq

type hash_table_eq
   private
   integer :: size
   type(entry_ptr_eq), pointer :: table(:)
end type hash_table_eq

!----------------------------------------------------
! The following parameters are defined:

type(hash_key_eq), parameter :: NULL_KEY_EQ = hash_key_eq(-1,0)

!----------------------------------------------------
! The following variables are defined:

!----------------------------------------------------
! The public module procedures are:

!   subroutine hash_table_init(table,table_size)
! This routine initializes a hash table with size table_size
!   type(hash_table_eq), intent(out) :: table
!   type(table_size), intent(in) :: table_size
!   end subroutine hash_table_init

!   subroutine hash_table_destroy(table)
! This routine destroys a hash table and frees the associated memory
!   type(hash_table_eq), intent(inout) :: table
!   end subroutine hash_table_destroy

!   subroutine hash_insert(key,ind,table)
! This routine inserts (key,ind) in the hash table
!   type(hash_key_eq), intent(in) :: key
!   integer, intent(in) :: ind
!   type(hash_table_eq), intent(inout) :: table
!   end function hash_insert

!   subroutine hash_remove(key,table)
! This routine removes key from the hash table
!   type(hash_key_eq), intent(in) :: key
!   type(hash_table_eq), intent(inout) :: table
!   end subroutine hash_remove

!   function hash_decode_key(key,table)
! This routine returns the index associated with key
!   integer :: hash_decode_key
!   type(hash_key_eq), intent(in) :: key
!   type(hash_table_eq), intent(in) :: table
!   end function hash_decode_key

!   function hash_overflow(key,m,a)
! This routine returns true if m*key+a will overflow the hash values
!   logical :: hash_overflow
!   type(hash_key), intent(in) :: key
!   integer, intent(in) :: m, a
!   end function hash_overflow

!   subroutine hash_print_key(key,iounit,str)
! This routine prints a key
!   type(hash_key_eq), intent(in) :: key
!OR type(hash_key_eq), intent(in) :: key(:)
!   integer, intent(in) :: iounit
!   character(len=*), optional, intent(in) :: str
!   end subroutine hash_print_key

!   subroutine hash_pack_key(key,array,ind)
! This routine packs key into array starting at ind
!   type(hash_key_eq), intent(in) :: key
!OR type(hash_key_eq), intent(in) :: key(:)
!   integer, intent(inout) :: array(:)
!   integer, intent(in) :: ind
!   end subroutine hash_pack_key

!   function hash_unpack_key(array,ind,extended)
! This routine unpacks the key in array starting at ind
!   type(hash_key_eq) :: hash_unpack_key
!   integer, intent(in) :: array(:), ind
!   logical, intent(in) :: extended ! needs to be present to resolve generic
!   end function hash_unpack_key

!   subroutine hash_pack_key_pvm(keys,itype_pvm)
! This routine packs an array of keys into a PVM message buffer
!   type(hash_key_eq), intent(in) :: keys(:)
!   integer, intent(in) :: itype_pvm
!   end subroutine hash_pack_key_pvm

!   subroutine hash_unpack_key_pvm(keys,nkey,itype_pvm)
! This routine unpacks an array of keys from a PVM message buffer
!   type(hash_key_eq), intent(inout) :: keys(:)
!   integer, intent(in) :: nkey, itype_pvm
!   end subroutine hash_unpack_key_pvm

!----------------------------------------------------
! The following interface operators are defined:

interface assignment(=)
   module procedure key_equals_integer_eq
end interface

interface operator(*)
   module procedure integer_times_key_eq
end interface

interface operator(+)
   module procedure key_plus_integer_eq
end interface

interface operator(-)
   module procedure key_minus_integer_eq
end interface

interface operator(/)
   module procedure key_div_integer_eq
end interface

interface operator(==)
   module procedure key_compare_eq, key_integer_compare_eq
end interface

interface operator(<)
   module procedure key_less_than_eq
end interface

!----------------------------------------------------
! The following generic interfaces are defined:

interface hash_table_init
   module procedure hash_table_init_eq
end interface

interface hash_table_destroy
   module procedure hash_table_destroy_eq
end interface

interface hash_insert
   module procedure hash_insert_eq
end interface

interface hash_remove
   module procedure hash_remove_eq
end interface

interface hash_decode_key
   module procedure hash_decode_key_eq
end interface

interface hash_overflow
   module procedure hash_overflow_eq
end interface

interface mod
   module procedure hash_modfunc_eq
end interface

interface hash_pack_key
   module procedure hash_pack_key_scalar_eq, hash_pack_key_array_eq
end interface

interface hash_unpack_key
   module procedure hash_unpack_key_eq
end interface

interface hash_pack_key_pvm
   module procedure hash_pack_key_pvm_eq
end interface

interface hash_unpack_key_pvm
   module procedure hash_unpack_key_pvm_eq
end interface

interface hash_print_key
   module procedure hash_print_key_scalar_eq, hash_print_key_array_eq
end interface

!----------------------------------------------------

contains

!---------------------------------------------------------------------
!  INITIALIZATION AND FINALIZATION ROUTINES
!---------------------------------------------------------------------

!          ------------------
subroutine hash_table_init_eq(table,table_size)
!          ------------------

!----------------------------------------------------
! This routine initializes a hash table with size table_size
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_table_eq), intent(out) :: table
integer, intent(in) :: table_size

!----------------------------------------------------
! Local variables:

integer :: i, astat

!----------------------------------------------------
! Begin executable code

! check that hash_mod has the same size hash key

if (KEY_SIZE /= 1) then
   ierr = PHAML_INTERNAL_ERROR
   call fatal("hash_mod and hash_mod_eq do not have the same size hash key", &
              intlist=(/KEY_SIZE,1/))
   stop
endif

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
end subroutine hash_table_init_eq

!          ---------------------
subroutine hash_table_destroy_eq(table)
!          ---------------------

!----------------------------------------------------
! This routine destroys a hash table and frees the associated memory
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_table_eq), intent(inout) :: table

!----------------------------------------------------
! Local variables:

integer :: i, dstat
type(hash_entry_eq), pointer :: current, next

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
end subroutine hash_table_destroy_eq

!---------------------------------------------------------------------
!  HASH TABLE MANIPULATIONS
!---------------------------------------------------------------------

!          --------------
subroutine hash_insert_eq(key,ind,table)
!          --------------

!----------------------------------------------------
! This routine inserts (key,ind) in the hash table
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_key_eq), intent(in) :: key
integer, intent(in) :: ind
type(hash_table_eq), intent(inout) :: table

!----------------------------------------------------
! Local variables:

integer :: table_sub, astat
type(hash_entry_eq), pointer :: newentry, next_entry
!----------------------------------------------------

!----------------------------------------------------
! Begin executable code

! compute the table subscript corresponding to key

table_sub = hash_function_eq(key,table%size)

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
   call fatal("memory allocation failed in hash_insert_eq")
   ierr = ALLOC_FAILED
   return
endif
newentry%key = key
newentry%index = ind
newentry%next => table%table(table_sub)%ptr
table%table(table_sub)%ptr => newentry

return
end subroutine hash_insert_eq

!          --------------
subroutine hash_remove_eq(key,table)
!          --------------

!----------------------------------------------------
! This routine removes key from the hash table
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_key_eq), intent(in) :: key
type(hash_table_eq), intent(inout) :: table

!----------------------------------------------------
! Local variables:

integer :: table_sub, dstat
type(hash_entry_eq), pointer :: next_entry, previous

!----------------------------------------------------
! Begin executable code

! compute the table subscript corresponding to key

table_sub = hash_function_eq(key,table%size)

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
   call warning("hash key not matched in hash_remove_eq")
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
end subroutine hash_remove_eq

!        ------------------
function hash_decode_key_eq(key,table)
!        ------------------

!----------------------------------------------------
! This routine returns the index associated with key
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer :: hash_decode_key_eq
type(hash_key_eq), intent(in) :: key
type(hash_table_eq), intent(in) :: table

!----------------------------------------------------
! Local variables:

integer :: table_sub
type(hash_entry_eq), pointer :: next_entry

!----------------------------------------------------
! Begin executable code

! compute the table subscript corresponding to key

table_sub = hash_function_eq(key,table%size)

! find the key in the linked list at table_sub

next_entry => table%table(table_sub)%ptr
do
   if (.not.associated(next_entry)) exit
   if (next_entry%key == key) exit
   next_entry => next_entry%next
end do

! return the index

if (.not.associated(next_entry)) then
   hash_decode_key_eq = HASH_NOT_FOUND
else
   hash_decode_key_eq = next_entry%index
endif

return
end function hash_decode_key_eq

!---------------------------------------------------------------------
!  PROCEDURES FOR OPERATOR INTERFACES
!---------------------------------------------------------------------

!          ---------------------
subroutine key_equals_integer_eq(key,int)
!          ---------------------

!----------------------------------------------------
! This routine assigns integer int and rank to a key
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_key_eq), intent(out) :: key
integer, intent(in) :: int(3)

!----------------------------------------------------
! Begin executable code

key%key = int(1)
key%rank = int(2)

return
end subroutine key_equals_integer_eq

!                  --------------------
recursive function integer_times_key_eq(int,key) result(res)
!                  --------------------

!----------------------------------------------------
! This routine multiplies an integer and a key
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_key_eq) :: res
integer, intent(in) :: int
type(hash_key_eq), intent(in) :: key

!----------------------------------------------------
! Begin executable code

if (key%key > huge(key%key)/int) then
   call fatal("overflow in hash key.  int and key are",intlist=(/int/))
   call hash_print_key(key,errunit)
endif

res%key = int*key%key
res%rank = key%rank

return
end function integer_times_key_eq

!        -------------------
function key_plus_integer_eq(key,int)
!        -------------------

!----------------------------------------------------
! This routine adds an integer to a key
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_key_eq) :: key_plus_integer_eq
type(hash_key_eq), intent(in) :: key
integer, intent(in) :: int

!----------------------------------------------------
! Begin executable code

key_plus_integer_eq%key = int+key%key
key_plus_integer_eq%rank = key%rank

return
end function key_plus_integer_eq

!        --------------------
function key_minus_integer_eq(key,int)
!        --------------------

!----------------------------------------------------
! This routine subtracts an integer from a key
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_key_eq) :: key_minus_integer_eq
type(hash_key_eq), intent(in) :: key
integer, intent(in) :: int

!----------------------------------------------------
! Begin executable code

key_minus_integer_eq%key = key%key - int
key_minus_integer_eq%rank = key%rank

return
end function key_minus_integer_eq

!                  ------------------
recursive function key_div_integer_eq(key,int) result(res)
!                  ------------------

!----------------------------------------------------
! This routine divides a key by an integer
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_key_eq) :: res
type(hash_key_eq), intent(in) :: key
integer, intent(in) :: int

!----------------------------------------------------
! Begin executable code

res%key = key%key/int
res%rank = key%rank

return
end function key_div_integer_eq

!        --------------
function key_compare_eq(key1,key2)
!        --------------

!----------------------------------------------------
! This routine determines if two keys are equal
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

logical :: key_compare_eq
type(hash_key_eq), intent(in) :: key1, key2

!----------------------------------------------------
! Begin executable code

key_compare_eq = (key1%key==key2%key .and.  key1%rank==key2%rank)

return
end function key_compare_eq

!        ----------------------
function key_integer_compare_eq(key,int)
!        ----------------------

!----------------------------------------------------
! This routine determines if a key is equal to an integer
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

logical :: key_integer_compare_eq
type(hash_key_eq), intent(in) :: key
integer, intent(in) :: int

!----------------------------------------------------
! Begin executable code

key_integer_compare_eq = key%key==int

return
end function key_integer_compare_eq

!        ----------------
function key_less_than_eq(key1,key2)
!        ----------------

!----------------------------------------------------
! This routine determines if key1 is less than key2
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

logical :: key_less_than_eq
type(hash_key_eq), intent(in) :: key1, key2

!----------------------------------------------------
! Begin executable code

key_less_than_eq = key1%key < key2%key

return
end function key_less_than_eq

!---------------------------------------------------------------------
!  I/O
!---------------------------------------------------------------------

!          ------------------------
subroutine hash_print_key_scalar_eq(key,iounit,str)
!          ------------------------

!----------------------------------------------------
! This routine prints a key
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_key_eq), intent(in) :: key
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
end subroutine hash_print_key_scalar_eq

!          -----------------------
subroutine hash_print_key_array_eq(key,iounit,str)
!          -----------------------

!----------------------------------------------------
! This routine prints an array of keys
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_key_eq), intent(in) :: key(:)
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
end subroutine hash_print_key_array_eq

!---------------------------------------------------------------------
!  OTHER
!---------------------------------------------------------------------

!        ----------------
function hash_overflow_eq(key,m,a)
!        ----------------

!----------------------------------------------------
! This routine returns true if m*key+a will overflow the hash values
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

logical :: hash_overflow_eq
type(hash_key_eq), intent(in) :: key
integer, intent(in) :: m, a
!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

hash_overflow_eq = key%key > (huge(key%key)-a)/m

end function hash_overflow_eq

!        ----------------
function hash_function_eq(key,table_size)
!        ----------------

!----------------------------------------------------
! This routine evaluates the hash function for the given key.
! Based on Robert Jenkin's 32-bit mixing function, followed by a
! mod to place it in the range 1:table_size.  It only hashes on the
! low order word of the key.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer :: hash_function_eq
type(hash_key_eq), intent(in) :: key
integer, intent(in) :: table_size

!----------------------------------------------------
! Local variables:

integer :: lkey

!----------------------------------------------------
! Begin executable code

lkey = key%key + key%rank

lkey = lkey + ishft(lkey,12)
lkey = ieor(lkey,ishft(lkey,-22))
lkey = lkey + ishft(lkey,4)
lkey = ieor(lkey,ishft(lkey,-9))
lkey = lkey + ishft(lkey,10)
lkey = ieor(lkey,ishft(lkey,-2))
lkey = lkey + ishft(lkey,7)
lkey = ieor(lkey,ishft(lkey,-12))

hash_function_eq = 1 + mod(abs(lkey),table_size)

return
end function hash_function_eq

!          -----------------------
subroutine hash_pack_key_scalar_eq(key,array,ind)
!          -----------------------

!----------------------------------------------------
! This routine packs a single key into array starting at ind
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_key_eq), intent(in) :: key
integer, intent(inout) :: array(:)
integer, intent(in) :: ind

!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

array(ind) = key%key
array(ind+1) = key%rank

return
end subroutine hash_pack_key_scalar_eq

!          ----------------------
subroutine hash_pack_key_array_eq(key,array,ind)
!          ----------------------

!----------------------------------------------------
! This routine packs an array of keys into array starting at ind
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_key_eq), intent(in) :: key(:)
integer, intent(inout) :: array(:)
integer, intent(in) :: ind

!----------------------------------------------------
! Local variables:

integer :: i
!----------------------------------------------------
! Begin executable code

do i=0,size(key)-1
   array(ind+2*i  ) = key(i+1)%key
   array(ind+2*i+1) = key(i+1)%rank
end do

return
end subroutine hash_pack_key_array_eq

!        ------------------
function hash_unpack_key_eq(array,ind,extended)
!        ------------------

!----------------------------------------------------
! This routine unpacks the key in array starting at ind
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_key_eq) :: hash_unpack_key_eq
integer, intent(in) :: array(:), ind
logical, intent(in) :: extended

!----------------------------------------------------
! Local variables:

!----------------------------------------------------
! Begin executable code

hash_unpack_key_eq%key = array(ind)
hash_unpack_key_eq%rank = array(ind+1)

return
end function hash_unpack_key_eq

!          --------------------
subroutine hash_pack_key_pvm_eq(keys,itype_pvm)
!          --------------------

!----------------------------------------------------
! This routine packs an array of keys into a PVM message buffer.  This is
! used so we don't have to copy the keys into an integer array which then
! gets packed into the message.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_key_eq), intent(in) :: keys(:)
integer, intent(in) :: itype_pvm

!----------------------------------------------------
! Local variables:

integer :: info, i
!----------------------------------------------------
! Begin executable code

do i=1,size(keys)
   call pvmfpack(itype_pvm,(/keys(i)%key/),1,1,info)
   call pvmfpack(itype_pvm,(/keys(i)%rank/),1,1,info)
end do

return
end subroutine hash_pack_key_pvm_eq

!          ----------------------
subroutine hash_unpack_key_pvm_eq(keys,nkey,itype_pvm)
!          ----------------------

!----------------------------------------------------
! This routine unpacks an array of keys from a PVM message buffer.  This is
! used so we don't have to copy the keys from an integer array which was
! unpacked from the message.
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

type(hash_key_eq), intent(inout) :: keys(:)
integer, intent(in) :: nkey, itype_pvm

!----------------------------------------------------
! Local variables:

integer :: info, i, unpacked(2)
!----------------------------------------------------
! Begin executable code

do i=1,nkey
   call pvmfunpack(itype_pvm,unpacked,3,1,info)
   keys(i)%key = unpacked(1)
   keys(i)%rank = unpacked(2)
end do

end subroutine hash_unpack_key_pvm_eq

!        ---------------
function hash_modfunc_eq(key,int)
!        ---------------

!----------------------------------------------------
! This routine extends the generic mod function to hash keys
!----------------------------------------------------

!----------------------------------------------------
! Dummy arguments

integer :: hash_modfunc_eq
type(hash_key_eq), intent(in) :: key
integer, intent(in) :: int

!----------------------------------------------------
! Begin executable code

hash_modfunc_eq = mod(key%key,int)

end function hash_modfunc_eq

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

end module hash_eq_mod
