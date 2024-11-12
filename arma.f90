module arma_mod
use kind_mod, only: dp
use random_mod, only: random_normal
use linear_algebra_mod, only: solve_linear_system
implicit none
private
public :: arma_1_1_sim, ar_1_sim, ma_1_sim, ma_sim, arma_1_1_acf, ma_1_acf, &
   ma_2_acf, ar_sim, ar_2_acf, yule_walker, ar_pacf, pacf_from_acf, pacf_ma1
contains
function arma_1_1_sim(n, ar1, ma1, xsd) result(x)
! simulate a zero-mean ARMA(1,1) process with normally distributed noise
integer      , intent(in) :: n
real(kind=dp), intent(in) :: ar1  ! AR(1) coefficient
real(kind=dp), intent(in) :: ma1  ! MA(1) coefficient
real(kind=dp), intent(in) :: xsd  ! standard deviation of noise
real(kind=dp)             :: x(n) ! simulated data
real(kind=dp)             :: noise, old_noise
integer                   :: i
if (n < 1) return
noise = xsd*random_normal()
x(1) = noise
do i=2,n
   old_noise = noise
   noise = xsd*random_normal()
   x(i) = ar1*x(i-1) + ma1*old_noise + noise
end do
end function arma_1_1_sim
!
function ma_1_sim(n, ma1, xsd) result(x)
! simulate a zero-mean MA(1) process with normally distributed noise
integer      , intent(in) :: n
real(kind=dp), intent(in) :: ma1  ! MA(1) coefficient
real(kind=dp), intent(in) :: xsd  ! standard deviation of noise
real(kind=dp)             :: x(n) ! simulated data
real(kind=dp)             :: noise, old_noise
integer                   :: i
if (n < 1) return
noise = xsd*random_normal()
x(1) = noise
do i=2,n
   old_noise = noise
   noise = xsd*random_normal()
   x(i) = ma1*old_noise + noise
end do
end function ma_1_sim
!
function ma_sim(n, ma, xsd) result(x)
! simulate a zero-mean MA(q) (moving average) process with normally 
! distributed noise
integer      , intent(in) :: n
real(kind=dp), intent(in) :: ma(:) ! MA coefficients
real(kind=dp), intent(in) :: xsd   ! standard deviation of noise
real(kind=dp)             :: x(n)  ! simulated data
real(kind=dp)             :: noise(size(x))
integer                   :: i, k, kmax, nma
nma = size(ma)
if (n < 1) return
do i=1,n
   noise(i) = xsd*random_normal()
end do
do i=1,n
   x(i) = noise(i)
   kmax = min(nma, i-1)
   do k=1,kmax
      x(i) = x(i) + ma(k) * noise(i-k)
   end do
end do
end function ma_sim
!
function ar_sim(n, ar, xsd) result(x)
! simulate a zero-mean AR(p) (moving average) process with normally 
! distributed noise
integer      , intent(in) :: n
real(kind=dp), intent(in) :: ar(:) ! AR coefficients
real(kind=dp), intent(in) :: xsd   ! standard deviation of noise
real(kind=dp)             :: x(n)  ! simulated data
real(kind=dp)             :: noise(size(x))
integer                   :: i, k, kmax, nar
nar = size(ar)
if (n < 1) return
do i=1,n
   noise(i) = xsd*random_normal()
end do
do i=1,n
   x(i) = xsd*random_normal()
   kmax = min(nar, i-1)
   do k=1,kmax
      x(i) = x(i) + ar(k)*x(i-k)
   end do
end do
end function ar_sim
!
function ar_1_sim(n, ar1, xsd) result(x)
! simulate a zero-mean AR(1) process with normally distributed noise
integer      , intent(in) :: n
real(kind=dp), intent(in) :: ar1  ! AR(1) coefficient
real(kind=dp), intent(in) :: xsd  ! standard deviation of noise
real(kind=dp)             :: x(n) ! simulated data
integer                   :: i
if (n < 1) return
x(1) = xsd*random_normal()
do i=2,n
   x(i) = ar1*x(i-1) + xsd*random_normal()
end do
end function ar_1_sim
!
function ma_1_acf(ma1) result(acf1)
  ! Returns the first order theoretical autocorrelations of an MA(1) process
  real(kind=dp), intent(in) :: ma1   ! MA(1) parameter
  real(kind=dp) :: acf1              ! ACF(1)
  acf1 = ma1/(1 + ma1**2)
end function ma_1_acf
!
function ma_2_acf(ma1, ma2) result(acfs)
! Returns the first 2 theoretical autocorrelations of an MA(2) process
real(kind=dp), intent(in) :: ma1, ma2
real(kind=dp)             :: acfs(2)
acfs = [(ma1 + ma1*ma2), ma2] / (1 + ma1**2 + ma2**2)
end function ma_2_acf
!
function arma_1_1_acf(ar1, ma1, nacf) result(xacf)
  ! Returns the first nacf theoretical autocorrelations of an ARMA(1,1) process
  real(kind=dp), intent(in) :: ar1         ! AR(1) parameter
  real(kind=dp), intent(in) :: ma1         ! MA(1) parameter
  integer, intent(in) :: nacf              ! Number of autocorrelations to compute
  real(kind=dp) :: xacf(nacf)              ! Output array for theoretical autocorrelations
  integer :: k
  if (nacf < 1) return
  xacf(1) = (ar1 + ma1) * (1 + ar1*ma1) / (1 + 2*ar1*ma1 + ma1**2)
  do k = 2, nacf
    xacf(k) = ar1**(k-1) * xacf(1)
  end do
end function arma_1_1_acf
!
function ar_2_acf(ar1, ar2, nacf) result(xacf)
  ! Returns the first nacf theoretical autocorrelations of an AR(2) process
  real(kind=dp), intent(in) :: ar1         ! AR(1) parameter
  real(kind=dp), intent(in) :: ar2         ! AR(2) parameter
  integer, intent(in) :: nacf              ! Number of autocorrelations to compute
  real(kind=dp) :: xacf(nacf)              ! Output array for theoretical autocorrelations
  integer :: k
  if (nacf < 1) return
  xacf(1) = ar1 / (1 - ar2)
  if (nacf < 2) return
  xacf(2) = ar1*xacf(1) + ar2
  do k = 3, nacf
    xacf(k) = ar1*xacf(k-1) + ar2*xacf(k-2)
  end do
end function ar_2_acf

function yule_walker(acf) result(ar)
  ! Fits an AR(n) model by solving the Yule-Walker equations using the Levinson-Durbin 
  ! algorithm. acf should include autocorrelations from lags 1 through n
  real(kind=dp), intent(in) :: acf(:)   ! Input autocorrelations (size n)
  real(kind=dp) :: ar(size(acf))        ! Output AR parameters (size n)
  integer :: n                          ! Order of the AR model
  real(kind=dp), allocatable :: a(:), k(:)
  integer :: i
  real(kind=dp) :: error
  n = size(acf)
  if (n < 1) return
  allocate(a(n), k(n))
  ! Initialize
  k(1) = acf(1)
  ar(1) = k(1)
  error = 1.0_dp - k(1)**2
  do i = 2, n
     k(i) = (acf(i) - sum(ar(1:i-1) * acf(i-1:1:-1))) / error
     a(1:i-1) = ar(1:i-1) - k(i) * ar(i-1:1:-1)
     ar(:i) = [a(1:i-1), k(i)]
     error = error * (1.0_dp - k(i)**2)
  end do
end function yule_walker

function ar_pacf(ar, nacf) result(pacf_vals)
  ! Returns the first nacf partial autocorrelations of an AR(p) process given the AR coefficients.
  !
  ! Arguments:
  ! ar        : real(kind=dp), intent(in)  - AR coefficients (phi_1 to phi_p)
  ! nacf      : integer, intent(in)        - Number of partial autocorrelations to compute
  !
  ! Returns:
  ! pacf_vals : real(kind=dp)              - Partial autocorrelations at lags 1 to nacf
  real(kind=dp), intent(in) :: ar(:)
  integer, intent(in) :: nacf
  real(kind=dp) :: pacf_vals(nacf)
  integer :: imax
  imax = min(nacf, size(ar))
  pacf_vals = 0.0_dp
  pacf_vals(1:imax) = ar(1:imax)
end function ar_pacf

function pacf_from_acf(acf) result(pacf_vals)
  ! Computes the partial autocorrelation function (PACF) from input autocorrelations.
  ! Arguments:
  ! acf       : real(kind=dp), intent(in)  - Array of autocorrelations at lags 1 to nacf
  ! Returns:
  ! pacf_vals : real(kind=dp)              - Partial autocorrelations at lags 1 to nacf
  real(kind=dp), intent(in) :: acf(:)         ! Input autocorrelations (size nacf)
  real(kind=dp) :: pacf_vals(size(acf))       ! Output PACF values at lags 1 to nacf
  integer :: nacf
  real(kind=dp), allocatable :: a(:), a_temp(:), k(:)
  integer :: i
  nacf = size(acf)
  if (nacf < 1) return
  allocate(a(nacf), k(nacf))
  ! Initialize first partial autocorrelation
  k(1) = acf(1)
  pacf_vals(1) = k(1)
  a(1) = k(1)
  if (nacf >= 2) then
    do i = 2, nacf
      k(i) = (acf(i) - sum(a(1:i-1) * acf(i-1:1:-1))) / (1.0_dp - sum(a(1:i-1)**2))
      a_temp = a(1:i-1)
      a(1:i-1) = a_temp - k(i) * a_temp(i-1:1:-1)
      a(i) = k(i)
      pacf_vals(i) = k(i)
    end do
  end if
end function pacf_from_acf

pure function pacf_ma1(ma1, npacf) result(pacf_vals)
! return the partial autocorrelations of an MA(1) process
! from notes by Donald B. Percival https://faculty.washington.edu/dbp/s519/
! at https://faculty.washington.edu/dbp/s519/PDFs/12-overheads-2020.pdf
real(kind=dp), intent(in) :: ma1
integer, intent(in) :: npacf
real(kind=dp) :: pacf_vals(npacf)
integer :: i
real(kind=dp) :: numer, denom
denom = 1.0_dp
do i=1,npacf
   denom = denom + ma1**(2*i)
   numer = - ((-ma1) ** i)
   pacf_vals(i) = numer / denom
end do
end function pacf_ma1

end module arma_mod
