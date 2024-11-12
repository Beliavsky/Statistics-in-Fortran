module arch_var_mod
use kind_mod, only: wp => dp
use util_mod, only: moving_average_vec, weighted_moving_average_vec
implicit none
private
public :: arch_var_equal_wgt, arch_var, asymm_arch_var, shift_twist_square
integer, parameter :: itwist = 3, ishift = 4
contains
pure function arch_var_equal_wgt(x, nlags, arch_coeff, arch_var_min) result(xvar)
! compute ARCH variance predictions where the ARCH weights are equal up to lag nlags
real(kind=wp), intent(in) :: x(:)
integer      , intent(in) :: nlags
real(kind=wp), intent(in) :: arch_coeff(:) ! (2) [variance_constant, arch_coeff]
real(kind=wp), intent(in), optional :: arch_var_min
real(kind=wp)             :: xvar(size(x))
real(kind=wp)             :: avg_sq(size(x))
integer                   :: i
avg_sq = moving_average_vec(x**2,nlags)
xvar(1) = arch_coeff(1)
do i=2,size(x)
   xvar(i) = arch_coeff(1) + arch_coeff(2)*avg_sq(i-1)
end do
if (present(arch_var_min)) xvar = max(arch_var_min,xvar)
end function arch_var_equal_wgt
!
pure function arch_var(x, arch_wgt, arch_coeff, arch_var_min) result(xvar)
! compute ARCH variance predictions
real(kind=wp), intent(in)           :: x(:)
real(kind=wp), intent(in)           :: arch_wgt(:)
real(kind=wp), intent(in)           :: arch_coeff(:) ! (2) [variance_constant, arch_coeff]
real(kind=wp), intent(in), optional :: arch_var_min
real(kind=wp)                       :: xvar(size(x))
real(kind=wp)                       :: avg_sq(size(x))
integer                             :: i
avg_sq = weighted_moving_average_vec(x**2,arch_wgt)
xvar(1) = arch_coeff(1)
do i=2,size(x)
   xvar(i) = arch_coeff(1) + arch_coeff(2)*avg_sq(i-1)
end do
if (present(arch_var_min)) xvar = max(arch_var_min,xvar)
end function arch_var
!
function asymm_arch_var(x, arch_wgt, arch_alpha, arch_beta, arch_twist, arch_shift, arch_var_min) result(xvar)
! compute asymmetric ARCH variance predictions
real(kind=wp), intent(in)           :: x(:)
real(kind=wp), intent(in)           :: arch_wgt(:)
real(kind=wp), intent(in)           :: arch_alpha, arch_beta, arch_twist, arch_shift
real(kind=wp), intent(in), optional :: arch_var_min
real(kind=wp)                       :: xvar(size(x))
real(kind=wp)                       :: avg_sq(size(x))
integer                             :: i
avg_sq = weighted_moving_average_vec(shift_twist_square(x,shift=arch_shift,twist=arch_twist),arch_wgt)
xvar(1) = arch_alpha
do i=2,size(x)
   xvar(i) = arch_alpha + arch_beta*avg_sq(i-1)
end do
if (present(arch_var_min)) xvar = max(arch_var_min,xvar)
end function asymm_arch_var
!
elemental function shift_twist_square(x,shift,twist) result(xsq)
real(kind=wp), intent(in)           :: x
real(kind=wp), intent(in), optional :: shift, twist
real(kind=wp)                       :: xsq
real(kind=wp)                       :: xshift
if (present(shift)) then
   xshift = x - shift
   if (present(twist)) then
      xsq =  merge(1+twist,1-twist,xshift<0.0_wp) * xshift**2
   else
      xsq = xshift**2
   end if
else
   if (present(twist)) then
      xsq = merge(1+twist,1-twist,x<0.0_wp) * x**2
   else
      xsq = x**2
   end if
end if
end function shift_twist_square
end module arch_var_mod