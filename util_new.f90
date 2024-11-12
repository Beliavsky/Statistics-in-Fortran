module util_mod
use kind_mod, only: dp
implicit none
private
public :: bad_real, moving_average_vec, weighted_moving_average_vec, sech, &
   good_pos, grid, exclude, prepend, assert, default, combine_mat_tensor, &
   remove_char, min_optional, set_optional, match_string, file_exists, &
   first_false, close_file, num_fields, strip, subset_vec_mat, &
   weighted_sums, c, default_real, set_alloc, read_words_line, exe_name, &
   above_diag, print_wall_time_elapsed, spline_power_wgt
interface assert
   module procedure assert_scalar,assert_vec
end interface assert
interface default
   module procedure default_character, default_logical, default_real, &
      default_integer
end interface default
interface good_pos
   module procedure good_pos_int,good_pos_char
end interface good_pos
interface grid
   module procedure grid_int,grid_real
end interface grid
interface exclude
   module procedure exclude_vec_int,exclude_vec_real
end interface exclude
interface match_string
   module procedure match_string_scalar,match_string_vec
end interface match_string
interface min_optional
   module procedure min_optional_int
end interface min_optional
interface prepend
   module procedure prepend_char
end interface prepend
interface combine_mat_tensor
   module procedure combine_mat_tensor_real
end interface combine_mat_tensor
interface set_optional
   module procedure set_optional_real, set_optional_integer,set_optional_logical, &
                    set_optional_character
end interface set_optional
interface set_alloc
   module procedure set_alloc_real_vec, set_alloc_real_matrix
end interface set_alloc
character (len=*), parameter :: mod_str = "in util_mod::"
real(kind=dp), parameter :: bad_real = -999.0_dp
integer, parameter :: istdout = 6
contains
pure function moving_average_vec(xx, nma) result(xma)
! simple moving average with nma terms
! uses growing average for xma(1:nma-1)
real(kind=dp), intent(in) :: xx(:)
integer, intent(in) :: nma
real(kind=dp) :: xma(size(xx))
integer :: i
! Initialize the moving average for the first nma-1 elements
do i = 1, min(nma-1, size(xx))
    xma(i) = sum(xx(1:i)) / i
end do
! Calculate the moving average for the remaining elements
do i = nma, size(xx)
   xma(i) = sum(xx(i-nma+1:i)) / nma
end do
end function moving_average_vec
!
pure function weighted_moving_average_vec(xx,wgt) result(xma)
! weighted moving average with nma terms
! uses growing average for xma(1:nma-1)
real(kind=dp), intent(in) :: xx(:)
real(kind=dp), intent(in) :: wgt(:)
real(kind=dp)             :: xma(size(xx))
integer                   :: i,m,n,nma
nma = size(wgt)
if (nma < 1) return
n = size(xx)
if (n < 1) return
if (nma > n) then
   xma = moving_average_vec(xx,n)
   return
end if
do i=1,nma-1
   m = min(i,nma)
   xma(i) = sum(xx(i-m+1:i))/m       
end do
do i=nma,n
   xma(i) = sum(xx(i-nma+1:i)*wgt(nma:1:-1))
end do
end function weighted_moving_average_vec
!
elemental function sech(xx) result(yy)
real(kind=dp), intent(in) :: xx
real(kind=dp)             :: yy
yy = 1/cos(xx)
end function sech
!
subroutine good_pos_int(ivec, ipos, values_include, values_exclude, imin, imax)
    ! return in ipos the positions in ivec that are not in values_exclude and that are in the range (imin, imax), inclusive
    integer, intent(in)               :: ivec(:)
    integer, intent(out), allocatable :: ipos(:)
    integer, intent(in) , optional    :: values_include(:) ! if present, allowed values in ivec
    integer, intent(in) , optional    :: values_exclude(:) ! values to exclude from ivec
    integer, intent(in) , optional    :: imin     ! min allowable value in ivec
    integer, intent(in) , optional    :: imax     ! max allowable value in ivec
    integer :: i, ncount, n
    logical :: valid

    ! Initialize ncount to zero to keep track of valid positions
    ncount = 0

    ! First, determine how many valid positions there are
    do i = 1, size(ivec)
        valid = .true.

        ! Check if value is within the optional range [imin, imax]
        if (present(imin)) then
            if (ivec(i) < imin) valid = .false.
        end if
        if (present(imax)) then
            if (ivec(i) > imax) valid = .false.
        end if

        ! Check if value is in the optional values_include list
        if (present(values_include)) then
            if (.not. any(ivec(i) == values_include)) valid = .false.
        end if

        ! Check if value is in the optional values_exclude list
        if (present(values_exclude)) then
            if (any(ivec(i) == values_exclude)) valid = .false.
        end if

        if (valid) ncount = ncount + 1
    end do

    ! Allocate the ipos array based on the number of valid positions
    allocate(ipos(ncount))

    ! Now, populate the ipos array with the indices of the valid positions
    n = 0
    do i = 1, size(ivec)
        valid = .true.

        ! Repeat the same checks as above
        if (present(imin)) then
            if (ivec(i) < imin) valid = .false.
        end if
        if (present(imax)) then
            if (ivec(i) > imax) valid = .false.
        end if
        if (present(values_include)) then
            if (.not. any(ivec(i) == values_include)) valid = .false.
        end if
        if (present(values_exclude)) then
            if (any(ivec(i) == values_exclude)) valid = .false.
        end if

        if (valid) then
            n = n + 1
            ipos(n) = i
        end if
    end do
end subroutine good_pos_int
!
subroutine good_pos_char(ivec, ipos, values_include, values_exclude)
! return in ipos the positions in ivec that are not in values_exclude and are in values_include (if provided)
character(len=*), intent(in)               :: ivec(:)
integer         , intent(out), allocatable :: ipos(:)
character(len=*), intent(in) , optional    :: values_include(:) ! if present, allowed values in ivec
character(len=*), intent(in) , optional    :: values_exclude(:) ! values to exclude from ivec

integer :: i, count, n
logical :: valid

count = 0

do i = 1, size(ivec)
    valid = .true.

    if (present(values_include)) then
        if (.not. any(ivec(i) == values_include)) valid = .false.
    end if

    if (present(values_exclude)) then
        if (any(ivec(i) == values_exclude)) valid = .false.
    end if

    if (valid) count = count + 1
end do

allocate(ipos(count))

n = 0
do i = 1, size(ivec)
    valid = .true.

    if (present(values_include)) then
        if (.not. any(ivec(i) == values_include)) valid = .false.
    end if

    if (present(values_exclude)) then
        if (any(ivec(i) == values_exclude)) valid = .false.
    end if

    if (valid) then
        n = n + 1
        ipos(n) = i
    end if
end do

end subroutine good_pos_char
!
function grid_int(n, imin, inc) result(ivec)
! return a grid of n values, starting at imin,
! with increments of inc
integer, intent(in) :: n
integer, intent(in) :: imin, inc
integer :: ivec(n)

integer :: i

do i = 1, n
    ivec(i) = imin + (i-1) * inc
end do

end function grid_int
!
function grid_real(n, xmin, xh) result(xx)
! return a grid of n values, starting at xmin,
! with increments of xh
integer, intent(in) :: n
real(kind=dp), intent(in), optional :: xmin, xh
real(kind=dp) :: xx(n)

integer :: i
real(kind=dp) :: x_start, x_inc

! Set default values if xmin or xh are not provided
x_start = 0.0_dp
x_inc = 1.0_dp

if (present(xmin)) x_start = xmin
if (present(xh)) x_inc = xh

do i = 1, n
    xx(i) = x_start + (i-1) * x_inc
end do

end function grid_real
!
function exclude_vec_int(ivec, iexcl) result(jvec)
! return in jvec the elements in ivec except those in positions iexcl(:)
integer, intent(in)            :: ivec(:)
integer, intent(in), optional  :: iexcl(:)
integer, allocatable           :: jvec(:)

integer :: i, ncount
logical :: exclude_flag(size(ivec))

! Initialize the exclude_flag array to .false.
exclude_flag = .false.

! Mark the positions to be excluded
if (present(iexcl)) then
    do i = 1, size(iexcl)
        if (iexcl(i) >= 1 .and. iexcl(i) <= size(ivec)) then
            exclude_flag(iexcl(i)) = .true.
        end if
    end do
end if

! Count the number of elements that are not excluded
ncount = count(.not. exclude_flag)

! Allocate the result array
allocate(jvec(ncount))

! Fill jvec with the elements that are not excluded
ncount = 0
do i = 1, size(ivec)
    if (.not. exclude_flag(i)) then
        ncount = ncount + 1
        jvec(ncount) = ivec(i)
    end if
end do

end function exclude_vec_int
!
function exclude_vec_real(ivec, iexcl) result(jvec)
! return in jvec the elements in ivec except those in positions iexcl(:)
real(kind=dp), intent(in)            :: ivec(:)
integer      , intent(in), optional  :: iexcl(:)
real(kind=dp), allocatable           :: jvec(:)

integer :: i, ncount
logical :: exclude_flag(size(ivec))

! Initialize the exclude_flag array to .false.
exclude_flag = .false.

! Mark the positions to be excluded
if (present(iexcl)) then
    do i = 1, size(iexcl)
        if (iexcl(i) >= 1 .and. iexcl(i) <= size(ivec)) then
            exclude_flag(iexcl(i)) = .true.
        end if
    end do
end if

! Count the number of elements that are not excluded
ncount = count(.not. exclude_flag)

! Allocate the result array
allocate(jvec(ncount))

! Fill jvec with the elements that are not excluded
ncount = 0
do i = 1, size(ivec)
    if (.not. exclude_flag(i)) then
        ncount = ncount + 1
        jvec(ncount) = ivec(i)
    end if
end do

end function exclude_vec_real
!
function prepend_char(xadd, xx) result(yy)
! set yy equal to [xadd,xx]
character(len=*), intent(in)                  :: xadd(:)
character(len=*), intent(in)                  :: xx(:)
character(len=len(xx(1)))                     :: yy(size(xx) + size(xadd))

integer :: i

! Copy elements from xadd to the beginning of yy
do i = 1, size(xadd)
    yy(i) = xadd(i)
end do

! Copy elements from xx to the remaining positions in yy
do i = 1, size(xx)
    yy(size(xadd) + i) = xx(i)
end do

end function prepend_char
!
function combine_mat_tensor_real(xmat,xtens,idim,position) result(ytens)
! augment tensor xtens(:,:) by matrix xmat(:) along dimension idim
real(kind=dp)    , intent(in)           ::  xmat(:,:)    ! (m1,m2)
real(kind=dp)    , intent(in)           :: xtens(:,:,:)  ! (n1,n2,n3)
integer          , intent(in)           :: idim
character (len=*), intent(in), optional :: position
real(kind=dp)    , allocatable          :: ytens(:,:,:)
integer                                 :: n1,n2,n3
character (len=*), parameter            :: msg = mod_str // "combine_vec_mat, "
character (len=20)                      :: position_
position_ = default("before",position)
n1  = size(xtens,1)
n2  = size(xtens,2)
n3  = size(xtens,3)
call assert(idim>0 .and. idim<=3,msg // "need 1 <= idim <= 3",stop_error=.true.)
call assert(position_=="before" .or. position_ == "after","need position = 'before' or 'after'",stop_error=.true.)
if (any(shape(xmat) /= exclude(shape(xtens),iexcl=[idim]))) then
   write (*,*) msg,"idim =",idim," shape(xmat)=",shape(xmat)," shape(xtens)=",shape(xtens), &
               " cannot combine matrix with tensor along specified dimension"
   stop
end if
if (idim == 1) then ! add another row
   allocate (ytens(n1+1,n2,n3))
   if (position_ == "before") then
      ytens(1 ,:,:) = xmat
      ytens(2:,:,:) = xtens
   else
      ytens(n1+1,:,:) = xmat
      ytens(:n1 ,:,:) = xtens
   end if
else if (idim == 2) then ! add another column
   allocate (ytens(n1,n2+1,n3))
   if (position_ == "before") then
      ytens(:,1   ,:) = xmat
      ytens(:,2:  ,:) = xtens
   else
      ytens(:,n2+1,:) = xmat
      ytens(:,:n2 ,:) = xtens
   end if
else
   allocate (ytens(n1,n2,n3+1))
   if (position_ == "before") then
      ytens(:,:,1)  = xmat
      ytens(:,:,2:) = xtens
   else
      ytens(:,:,n3+1) = xmat
      ytens(:,:,:n3)  = xtens
   end if
end if
end function combine_mat_tensor_real
!
subroutine assert_scalar(tf,msg,stop_error)
! test if tf is true
logical          , intent(in)            :: tf
character (len=*), intent(in)            :: msg
logical          , intent(in) , optional :: stop_error
if (.not. tf) then
   write (*,*) "in " // trim(msg)
   if (present(stop_error)) then
      if (stop_error) then
         write (*,*) "stopping in " // mod_str // "assert_scalar"
         stop
      end if
   end if
end if
end subroutine assert_scalar
!
subroutine assert_vec(tf,msg,stop_error,jerr)
! test if all elements of tf are true
logical          , intent(in)            :: tf(:)
character (len=*), intent(in)            :: msg
logical          , intent(in) , optional :: stop_error
integer          , intent(out), optional :: jerr
integer                                  :: ierr
ierr = first_false(tf)
if (ierr /= 0) then
   write (*,*) "in " // msg // ", ierr =",ierr
   if (present(stop_error)) then
      if (stop_error) stop
   end if
end if
if (present(jerr)) jerr = ierr
end subroutine assert_vec
!
pure function first_false(tf) result(i1)
! return the location of the first false element in tf(:), 0 if all .true.
logical, intent(in) :: tf(:)
integer             :: i1

integer :: i

i1 = 0  ! Default to 0, assuming all elements are .true.

do i = 1, size(tf)
    if (.not. tf(i)) then
        i1 = i
        return
    end if
end do

end function first_false
!
elemental function default_real(def,opt) result(yy)
! return opt if it is present, otherwise def (real arguments and result)
real(kind=dp), intent(in)           :: def
real(kind=dp), intent(in), optional :: opt
real(kind=dp)                       :: yy
if (present(opt)) then
   yy = opt
else
   yy = def
end if
end function default_real
!
elemental function default_integer(def,opt) result(yy)
! return opt if it is present, otherwise def (real arguments and result)
integer, intent(in)           :: def
integer, intent(in), optional :: opt
integer                       :: yy
if (present(opt)) then
   yy = opt
else
   yy = def
end if
end function default_integer
!
elemental function default_character(def, opt) result(yy)
! return opt if it is present, otherwise def (character arguments and result)
character(len=*), intent(in)           :: def
character(len=*), intent(in), optional :: opt
character(len=len(def))                :: yy

if (present(opt)) then
    yy = opt
else
    yy = def
end if

end function default_character
!
elemental function default_logical(def, opt) result(tf)
! return opt if it is present, otherwise def (logical arguments and result)
logical, intent(in)           :: def
logical, intent(in), optional :: opt
logical                       :: tf

if (present(opt)) then
    tf = opt
else
    tf = def
end if

end function default_logical
!
elemental function remove_char(xx, xrem) result(yy)
! remove instances of character xrem from xx
character(len=*), intent(in)           :: xx
character(len=*), intent(in), optional :: xrem
character(len=len(xx))                 :: yy

integer :: i, pos
character(len=1) :: char_rem

! Determine the character to remove
if (present(xrem)) then
    char_rem = xrem(1:1)  ! Use the first character of xrem
else
    char_rem = ' '  ! Default to removing spaces if xrem is not provided
end if

! Initialize the result string yy to all spaces
yy = ' '

! Initialize the position counter
pos = 1

! Loop over each character in xx
do i = 1, len(xx)
    if (xx(i:i) /= char_rem) then
        yy(pos:pos) = xx(i:i)
        pos = pos + 1
    end if
end do

! Trim any trailing spaces from the result
yy = adjustl(yy)

end function remove_char
!
elemental function min_optional_int(ii,jj) result(mmin)
! return the minimum of ii and jj if jj is present, otherwise ii
integer, intent(in)           :: ii
integer, intent(in), optional :: jj
integer                       :: mmin
if (present(jj)) then
   mmin = min(ii,jj)
else
   mmin = ii
end if
end function min_optional_int
!
elemental subroutine set_optional_real(xx,xdef,xopt)
! set xx to (xopt,xdef) if xopt (is,is not) PRESENT
real(kind=dp), intent(out)          :: xx
real(kind=dp), intent(in)           :: xdef
real(kind=dp), intent(in), optional :: xopt
if (present(xopt)) then
   xx = xopt
else
   xx = xdef
end if
end subroutine set_optional_real
!
elemental subroutine set_optional_integer(i,idef,iopt)
! set i to (iopt,idef) if iopt (is,is not) PRESENT
integer, intent(out)          :: i
integer, intent(in)           :: idef
integer, intent(in), optional :: iopt
if (present(iopt)) then
   i = iopt
else
   i = idef
end if
end subroutine set_optional_integer
!
elemental subroutine set_optional_logical(xx,xdef,xopt)
! set xx to (xopt,xdef) if xopt (is,is not) PRESENT
logical, intent(out)          :: xx
logical, intent(in)           :: xdef
logical, intent(in), optional :: xopt
if (present(xopt)) then
   xx = xopt
else
   xx = xdef
end if
end subroutine set_optional_logical
!
elemental subroutine set_optional_character(xx,xdef,xopt)
! set xx to (xopt,xdef) if xopt (is,is not) PRESENT
character(len=*), intent(out)          :: xx
character(len=*), intent(in)           :: xdef
character(len=*), intent(in), optional :: xopt
if (present(xopt)) then
   xx = xopt
else
   xx = xdef
end if
end subroutine set_optional_character
!
pure function match_string_scalar(ii, ivec) result(ipos)
! return the position of the first element in ivec that
! matches ii, 0 otherwise
character(len=*), intent(in) :: ii, ivec(:)
integer                      :: ipos

integer :: i

ipos = 0  ! Default to 0, assuming no match is found

do i = 1, size(ivec)
    if (ii == ivec(i)) then
        ipos = i
        return  ! Exit the function as soon as a match is found
    end if
end do

end function match_string_scalar
!
pure function match_string_vec(ii, ivec) result(ipos)
! return the position of the first element in ivec that
! matches each element of ii, 0 otherwise
character(len=*), intent(in) :: ii(:), ivec(:)
integer                      :: ipos(size(ii))

integer :: i, j

! Initialize ipos array to 0
ipos = 0

! Loop over each element in ii
do i = 1, size(ii)
    ! Loop over each element in ivec
    do j = 1, size(ivec)
        if (ii(i) == ivec(j)) then
            ipos(i) = j
            exit  ! Exit the inner loop as soon as a match is found
        end if
    end do
end do

end function match_string_vec
!
function file_exists(xfile) result(tf)
! return .true. if a file xfile exists, otherwise .false.
character (len=*), intent(in) :: xfile
logical                       :: tf
inquire (file=xfile,exist=tf)
end function file_exists
!
subroutine close_file(xfile)
! close a file if it is open
character (len=*), intent(in) :: xfile
integer                       :: iu
inquire(file=xfile,number=iu)
if (iu > -1) close(iu)
end subroutine close_file
!
elemental function num_fields(text,delim,xstrip) result(n)
! find the number of fields when splitting text using delimiter delim
! if xstrip is .true. and presen, strip a delimiter appearing at the end
character (len=*), intent(in) :: text
character (len=1), intent(in), optional :: delim
logical          , intent(in), optional :: xstrip
character (len=1)                       :: delim_
logical                                 :: xstrip_
integer                                 :: n
xstrip_ = default(.true.,xstrip)
delim_ = default(",",delim)
if (xstrip_) then
   n = 1 + num_matching_char(strip(text,delim_),delim_)
else
   n = 1 + num_matching_char(text,delim_)
end if
end function num_fields
!
elemental function strip(xx, xchar) result(yy)
! replace trailing characters of xx that equal xchar with a space
character(len=*), intent(in) :: xx
character(len=1), intent(in) :: xchar
character(len=len(xx))       :: yy

integer :: i

yy = xx  ! Start by copying xx to yy

! Loop from the last character to the first
do i = len(xx), 1, -1
    if (yy(i:i) == xchar) then
        yy(i:i) = ' '  ! Replace trailing xchar with a space
    else
        exit  ! Stop when a non-xchar is found
    end if
end do

end function strip
!
elemental function num_matching_char(text,char_match) result(nc)
! return the number of characters of text matching char_match
character (len=*), intent(in) :: text
character (len=*), intent(in) :: char_match
integer                       :: i,nc
if (len(char_match) /= 1) then
   nc = -1
   return
end if
nc = 0
do i=1,len(text)
   if (text(i:i) == char_match) nc = nc + 1
end do
end function num_matching_char
!
subroutine subset_vec_mat(sym_find, sym, xmat, sym_found, ymat)
! return in sym_found the elements of sym_find found in sym and in ymat the corresponding columns of xmat
character(len=*), intent(in)               :: sym_find(:)
character(len=*), intent(in)               :: sym(:)       ! (nsym)
real(kind=dp),    intent(in)               :: xmat(:,:)    ! (n, nsym)
character(len=*), intent(out), allocatable :: sym_found(:) ! (nfound)
real(kind=dp),    intent(out), allocatable :: ymat(:,:)    ! (n, nfound)

integer :: i, j, nfound
integer :: nsym, n
logical :: found

nsym = size(sym)
n = size(xmat, 1)

! Initialize nfound to zero
nfound = 0

! First, determine how many elements of sym_find are found in sym
do i = 1, size(sym_find)
    found = .false.
    do j = 1, nsym
        if (sym_find(i) == sym(j)) then
            nfound = nfound + 1
            found = .true.
            exit
        endif
    end do
end do

! Allocate the result arrays based on the number of found elements
allocate(sym_found(nfound))
allocate(ymat(n, nfound))

! Fill the result arrays
nfound = 0
do i = 1, size(sym_find)
    do j = 1, nsym
        if (sym_find(i) == sym(j)) then
            nfound = nfound + 1
            sym_found(nfound) = sym(j)
            ymat(:, nfound) = xmat(:, j)
            exit
        endif
    end do
end do

end subroutine subset_vec_mat
!
pure function weighted_sums(xx,wgt) result(xma)
! products of xx(i-nwgt+1:i) and wgt(nwgt:1:-1)
real(kind=dp), intent(in) :: xx(:)
real(kind=dp), intent(in) :: wgt(:)
real(kind=dp)             :: xma(size(xx))
integer                   :: i,n,nwgt
xma = 0.0_dp
n = size(xx)
nwgt = size(wgt)
if (n < 1 .or. nwgt < 1 .or. nwgt > n) return
do i=nwgt,n
   xma(i) = sum(wgt(nwgt:1:-1)*xx(i-nwgt+1:i))
end do
end function weighted_sums

function c(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10) result(vec)
! return character array containing present arguments
character (len=*)  , intent(in), optional    :: x1,x2,x3,x4,x5,x6,x7,x8,x9,x10
character (len=100)            , allocatable :: vec(:)
character (len=100)            , allocatable :: vec_(:)
integer                                      :: n
allocate (vec_(10))
if (present(x1))  vec_(1)  = x1
if (present(x2))  vec_(2)  = x2
if (present(x3))  vec_(3)  = x3
if (present(x4))  vec_(4)  = x4
if (present(x5))  vec_(5)  = x5
if (present(x6))  vec_(6)  = x6
if (present(x7))  vec_(7)  = x7
if (present(x8))  vec_(8)  = x8
if (present(x9))  vec_(9)  = x9
if (present(x10)) vec_(10) = x10
n = count([present(x1),present(x2),present(x3),present(x4),present(x5), &
           present(x6),present(x7),present(x8),present(x9),present(x10)])
allocate (vec(n))
if (n > 0) vec = vec_(:n)
end function c
!
subroutine set_alloc_real_vec(xx,yy,nsize)
! for real vectors xx(:) and yy(:), allocate yy and set it to xx
real(kind=dp), intent(in)               :: xx(:)
real(kind=dp), intent(out), allocatable :: yy(:)
integer      , intent(out), optional    :: nsize
integer                                 :: n
n = size(xx)
allocate (yy(n))
if (present(nsize)) nsize = n
end subroutine set_alloc_real_vec
!
subroutine set_alloc_real_matrix(xx,yy)
! for real matrices xx(:,:) and yy(:,:), allocate yy and set it to xx
real(kind=dp) , intent(in)               :: xx(:,:)
real(kind=dp) , intent(out), allocatable :: yy(:,:)
allocate (yy(size(xx,1),size(xx,2)))
yy = xx
end subroutine set_alloc_real_matrix
!
subroutine read_words_line(iu, words, nread, echo, label)
! read words from line, where the line has the # of words followed by the words
! n word_1 word_2 ... word_n
integer, intent(in)               :: iu
character(len=*), intent(out), allocatable :: words(:)
integer, intent(out), optional    :: nread
character(len=*), intent(in), optional :: label
logical, intent(in), optional     :: echo

integer :: n, i
character(len=256) :: line

! Read the number of words (n) from the input unit
read(iu, '(A)') line
read(line, *) n

! Allocate the array to hold the words
allocate(words(n))

! Read the words from the line
read(line, *) n, (words(i), i=1, n)

! Echo the words if echo is true
if (present(echo)) then
    if (echo) then
        if (present(label)) then
            print *, trim(label), ": ", trim(line)
        else
            print *, trim(line)
        end if
    end if
endif

! Set the number of words read (nread)
if (present(nread)) then
    nread = n
endif

end subroutine read_words_line
!
function exe_name() result(xname)
! return the program name, including directory
character (len=1000) :: xname
character (len=1000) :: dir_name
call getcwd(dir_name)
call get_command_argument(0,xname)
xname = trim(dir_name) // "\" // trim(xname)
end function exe_name
!
pure function above_diag(xx) result(yy)
! return in yy(:) the above-diagonal elements of xx(:,:)
real(kind=dp), intent(in) :: xx(:,:)
real(kind=dp)             :: yy((size(xx,1)*size(xx,1)-size(xx,1))/2)

integer :: n, k, i, j

n = size(xx, 1)
k = 0

do j = 2, n
    do i = 1, j-1
        k = k + 1
        yy(k) = xx(i, j)
    end do
end do

end function above_diag
!
subroutine print_wall_time_elapsed(old_time,fmt_trailer,outu,t1_cpu)
! print the time elapsed since old_time by calling system_clock
integer          , intent(in), optional :: old_time ! initialized with call system_clock(old_time)
character (len=*), intent(in), optional :: fmt_trailer
integer          , intent(in), optional :: outu
real(kind=dp)    , intent(in), optional :: t1_cpu
integer                                 :: itick,new_time,outu_,old_time_
real(kind=dp)                           :: new_time_cpu
character (len=100)                     :: fmt_time_
if (present(old_time)) old_time_ = old_time
outu_ = default(istdout,outu)
fmt_time_ = "(/,1x,'wall time elapsed(s) = ',f9.3)"
call system_clock(new_time,itick)
write (outu_,fmt_time_) (new_time-old_time_)/dble(itick)
if (present(t1_cpu)) then
   call cpu_time(new_time_cpu)
   write (outu_,"(1x,' cpu time elapsed(s) = ',f9.3)") new_time_cpu-t1_cpu
end if
if (present(fmt_trailer)) then
   if (fmt_trailer /= "") write (outu_,fmt_trailer)
end if
end subroutine print_wall_time_elapsed
!
pure function spline_power_wgt(n,power) result(wgt)
! return weights that decay smoothly to zero by lag n+1
integer      , intent(in) :: n
real(kind=dp), intent(in) :: power
real(kind=dp)             :: wgt(n)
integer                   :: i
if (n < 1) return
do i=1,n
   wgt(i) = (n - i + 1) ** power
end do
wgt = wgt/sum(wgt)
end function spline_power_wgt
end module util_mod
