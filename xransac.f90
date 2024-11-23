program ransac
! Simulate from a multiple linear regression model with outliers and compare the 
! fits of RANSAC and OLS.
use kind_mod, only: dp
use loss_functions_mod, only: mean_square_error, square_error_loss
use ransac_mod, only: ransac_fit
use random_mod, only: random_normal
use linear_regressor_mod, only: LinearRegressor_fit
implicit none
real(dp), allocatable :: X(:,:), y(:), ygood(:), coeff(:), best_params(:), &
   ols_params(:)
integer, parameter :: nsim=5
integer, parameter :: nobs_print=0 ! # of observations to print for debugging
integer :: i, n, k, d, isim
real(dp) :: t, yconst, best_error, sd_bad, sd_noise
integer :: nobs, nobs_bad, ncol
character (len=*), parameter :: fmt_cr = "(a20,':',*(1x, f10.4))", &
   fmt_ci = "((a20,':',*(1x, i10)))"
nobs = 1000
nobs_bad = 10
ncol = 2
sd_noise = 1.0_dp ! SD of noise for all observations
sd_bad = 1000.0_dp ! SD of extra noise for outliers
! Set RANSAC parameters
n = 100     ! Minimum number of data points to estimate parameters
k = 100     ! Maximum iterations allowed
t = 0.05_dp ! Threshold value to determine if points are fit well
d = 10      ! Number of close data points required to assert model fits well
write (*,fmt_ci) "#obs", nobs
write (*,fmt_ci) "#bad_obs", nobs_bad
write (*,fmt_ci) "#var", ncol
write (*,fmt_cr) "sd_noise", sd_noise
write (*,fmt_cr) "sd_bad", sd_bad
do isim=1, nsim
   yconst = 10.0_dp + random_normal()
   x = random_normal(nobs, ncol)
   coeff = random_normal(ncol)
   write (*,*)
   write (*,fmt_cr) "true param", yconst, coeff
   ygood = yconst + matmul(x, coeff) + sd_noise*random_normal(nobs)
   y = ygood
   y(:nobs_bad) = y(:nobs_bad) + sd_bad*random_normal(nobs_bad)
   ! Fit the model using RANSAC
   call ransac_fit(X, y, n, k, t, d, square_error_loss, mean_square_error, &
       best_params, best_error)
   call LinearRegressor_fit(X, y, ols_params)
   write (*,fmt_cr) "RANSAC param", best_params
   write (*,fmt_cr) "OLS param", ols_params
   write (*,fmt_cr) "Best error", best_error
   if (nobs_print > 0) then
      write (*,"(/,'data',/,*(a12))") "i", "y(i)", "x(i,:)"
      do i=1,min(nobs_print, nobs)
         print "(i12, *(f12.4))", i, y(i), x(i,:)
      end do
   end if
end do
end program ransac
