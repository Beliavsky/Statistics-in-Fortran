module stats_mod
  use kind_mod, only: dp
  implicit none
  private
  public :: mean, sd, skew, kurtosis, acf, pacf, pacf_from_acf, acf_pacf

contains

  function mean(x) result(mean_val)
  ! return the mean of x
    real(kind=dp), intent(in) :: x(:)
    real(kind=dp) :: mean_val
    mean_val = sum(x) / size(x)
  end function mean

  function sd(x) result(sd_val)
  ! return the standard deviation of x
    real(kind=dp), intent(in) :: x(:)
    real(kind=dp) :: sd_val
    real(kind=dp) :: mean_x
    mean_x = mean(x)
    sd_val = sqrt(sum((x - mean_x)**2) / (size(x) - 1))
  end function sd

  function skew(x) result(skew_val)
  ! return the skewness of x
    real(kind=dp), intent(in) :: x(:)
    real(kind=dp) :: skew_val
    real(kind=dp) :: mean_x, sd_x
    integer :: n
    n = size(x)
    mean_x = mean(x)
    sd_x = sd(x)
    skew_val = sum(((x - mean_x) / sd_x)**3) / n
  end function skew

  function kurtosis(x) result(kurtosis_val)
  ! return the kurtosis of x
    real(kind=dp), intent(in) :: x(:)
    real(kind=dp) :: kurtosis_val
    real(kind=dp) :: mean_x, sd_x
    integer :: n
    n = size(x)
    mean_x = mean(x)
    sd_x = sd(x)
    kurtosis_val = sum(((x - mean_x) / sd_x)**4) / n - 3.0_dp
  end function kurtosis

function acf(x, nacf) result(xacf)
  ! return the autocorrelations at lags 1 through nacf
  real(kind=dp), intent(in) :: x(:)         ! Input array
  integer, intent(in) :: nacf               ! Number of autocorrelations to compute
  real(kind=dp) :: xacf(nacf)               ! Output array for autocorrelations
  real(kind=dp) :: denom
  real(kind=dp), allocatable :: xdm(:)      ! Demeaned version of x
  integer :: n, lag
  n = size(x)
  xdm = x - mean(x)                          ! Compute demeaned x
  denom = sum(xdm**2)
  ! Compute autocorrelation for each lag from 1 to nacf
  do lag = 1, nacf
    xacf(lag) = sum(xdm(1:n-lag) * xdm(lag+1:n)) / denom
  end do
end function acf

function pacf(x, npacf) result(xpacf)
  ! Computes the partial autocorrelations at lags 1 through npacf
  real(kind=dp), intent(in) :: x(:)          ! Input array
  integer, intent(in) :: npacf               ! Number of partial autocorrelations to compute
  real(kind=dp) :: xpacf(npacf)              ! Output array for partial autocorrelations
  real(kind=dp) :: denom
  real(kind=dp), allocatable :: xdm(:), acf_vals(:), a(:), k(:)
  integer :: n, lag
  n = size(x)
  xdm = x - mean(x)                          ! Demean x
  denom = sum(xdm**2)
  ! Compute the autocorrelations for lags 1 to npacf
  allocate(acf_vals(0:npacf))
  acf_vals(0) = 1.0_dp                       ! acf at lag 0 is 1
  do lag = 1, npacf
    acf_vals(lag) = sum(xdm(1:n-lag) * xdm(lag+1:n)) / denom
  end do
  ! Initialize arrays for Durbin-Levinson recursion
  allocate(a(npacf), k(npacf))
  ! Compute partial autocorrelations using recursion
  k(1) = acf_vals(1)
  xpacf(1) = k(1)
  if (npacf >= 2) then
    do lag = 2, npacf
      k(lag) = (acf_vals(lag) - sum(a(1:lag-1) * acf_vals(lag-1:1:-1))) / (1.0_dp - sum(a(1:lag-1)**2))
      a(1:lag) = [a(1:lag-1) - k(lag)*a(lag-1:1:-1), k(lag)]
      xpacf(lag) = k(lag)
    end do
  end if
end function pacf

function pacf_from_acf(acf) result(pacf_vals)
  ! Computes the partial autocorrelations (PACF) at lags 1 through nacf using autocorrelations.
  real(kind=dp), intent(in) :: acf(:)         ! Input autocorrelations (size nacf)
  real(kind=dp) :: pacf_vals(size(acf))       ! Output array for partial autocorrelations (size nacf)
  integer :: nacf                              ! Number of autocorrelations (and PACF values)
  real(kind=dp), allocatable :: a(:), k(:)
  integer :: i
  nacf = size(acf)
  allocate(a(nacf), k(nacf))
  ! Initialize first partial autocorrelation
  k(1) = acf(1)
  pacf_vals(1) = k(1)
  ! Compute partial autocorrelations using the Durbin-Levinson recursion
  if (nacf >= 2) then
    do i = 2, nacf
      k(i) = (acf(i) - sum(a(1:i-1) * acf(i-1:1:-1))) / (1.0_dp - sum(a(1:i-1)**2))
      ! Update the AR coefficients
      a(1:i) = [a(1:i-1) - k(i)*a(i-1:1:-1), k(i)]
      pacf_vals(i) = k(i)
    end do
  end if
end function pacf_from_acf

subroutine acf_pacf(x, xacf, xpacf)
! Computes both the autocorrelation function (ACF) and the partial autocorrelation function (PACF)
! for the input array x. The ACF and PACF values are stored in xacf and xpacf, respectively.
real(kind=dp), intent(in)  :: x(:) ! input time series
real(kind=dp), intent(out) :: xacf(:) ! Output array for autocorrelations at lags 1 through nacf
real(kind=dp), intent(out) :: xpacf(:) ! Output array for partial autocorrelations at lags 1 through nacf
integer                    :: nacf
nacf = size(xacf)
if (size(xpacf) /= nacf) then
   print*,"size(xacf), size(xpacf) =", nacf, size(xpacf), "must be equal"
   error stop
end if
xacf = acf(x, nacf)
xpacf = pacf_from_acf(xacf)
end subroutine acf_pacf

end module stats_mod
