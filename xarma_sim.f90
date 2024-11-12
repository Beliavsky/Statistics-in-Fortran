program xarma_sim
! simulate from AR(1), MA(1), and ARMA(1,1) (autoregressive moving average) models
! and compute the empirical and theoretical autocorrelations
use kind_mod , only: dp
use arma_mod , only: ar_1_sim, ma_1_sim, arma_1_1_sim, ma_sim, ar_sim, &
   arma_1_1_acf, ma_1_acf, ma_2_acf, ar_2_acf
use stats_mod, only: acf
implicit none
integer, parameter :: n = 10**6, nacf=5, niter=5
integer :: iter
real(kind=dp) :: ar1, ar2, ma1, ma2, xsd, x(n)
character (len=*), parameter :: fmt_cr="(a20, ':', *(1x, f8.4))", &
   fmt_ci="(a20, ':', *(1x, i8))"
ar1 = 0.4_dp
ar2 = -0.7_dp
ma1 = 0.5_dp
ma2 = 0.3_dp
xsd = 10.0_dp
print fmt_ci, "#obs", n
print fmt_cr, "AR(1)", ar1
print fmt_cr, "true ACF", arma_1_1_acf(ar1, 0.0_dp, nacf)
do iter=1,niter
   x = ar_1_sim(n, ar1, xsd)
   print fmt_cr, "est. ACF", acf(x, nacf)
end do
print*
print fmt_cr, "MA(1)", ma1
print fmt_cr, "true ACF", arma_1_1_acf(0.0_dp, ma1, nacf)
print fmt_cr, "true ACF", ma_1_acf(ma1)
do iter=1, niter
   x = ma_1_sim(n, ma1, xsd)
   print fmt_cr, "est. ACF", acf(x, nacf)
end do
print*
print fmt_cr, "AR(1)", ar1
print fmt_cr, "MA(1)", ma1
print fmt_cr, "true ACF", arma_1_1_acf(ar1, ma1, nacf)
do iter=1, niter
   x = arma_1_1_sim(n, ar1, ma1, xsd)
   print fmt_cr, "est. ACF", acf(x, nacf)
end do
print*
print fmt_cr, "MA(1)", ma1
print fmt_cr, "MA(2)", ma2
print fmt_cr, "true ACF", ma_2_acf(ma1, ma2)
do iter=1, niter
   x = ma_sim(n, [ma1, ma2], xsd)
   print fmt_cr, "est. ACF", acf(x, nacf)
end do
print*
print fmt_cr, "AR(1)", ar1
print fmt_cr, "AR(2)", ar2
print fmt_cr, "true ACF", ar_2_acf(ar1, ar2, nacf)
do iter=1, niter
   x = ar_sim(n, [ar1, ar2], xsd)
   print fmt_cr, "est. ACF", acf(x, nacf)
end do
end program xarma_sim
