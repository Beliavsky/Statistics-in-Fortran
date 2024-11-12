module util_mod
use kind_mod, only: dp
implicit none
private
public :: moving_average_vec, weighted_moving_average_vec, sech, bad_real
real(kind=dp), parameter :: bad_real = -huge(1.0_dp)
contains
recursive pure function moving_average_vec(xx,nma) result(xma)
! simple moving average with nma terms
! uses growing average for xma(1:nma-1)
real(kind=dp), intent(in) :: xx(:)
integer      , intent(in) :: nma
real(kind=dp)             :: xma(size(xx))
integer                   :: i,m,n
if (nma < 1) return
n = size(xx)
if (n < 1) return
if (nma > n) then
   xma = moving_average_vec(xx,n)
   return
end if
do i=1,nma
   m = min(i,nma)
   xma(i) = sum(xx(i-m+1:i))/m       
end do
do i=nma+1,n
   xma(i) = xma(i-1) + (xx(i)-xx(i-nma))/nma
end do
end function moving_average_vec
!
pure function weighted_moving_average_vec(xx,wgt) result(xma)
! simple moving average with nma terms
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
elemental function sech(x) result(res)
  ! Computes the hyperbolic secant of x
  real(kind=dp), intent(in) :: x
  real(kind=dp) :: res
  res = 1.0_dp / cosh(x)
end function sech
end module util_mod