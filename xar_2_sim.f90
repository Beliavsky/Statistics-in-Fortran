program xar_2_sim
! simulate from an AR(2) process and compare the theoretical and empirical autocorrelations
use kind_mod , only: dp
use arma_mod , only: ar_sim, ar_2_acf
use stats_mod, only: acf
implicit none
integer, parameter :: n = 10**6, nacf=5, nar=2, niter=3
real(kind=dp) :: xsd, ar(nar), x(n)
integer :: iter
character (len=*), parameter :: fmt_cr="(a20, ':', *(1x, f8.4))", &
   fmt_ci="(a20, ':', *(1x, i8))"
print fmt_ci, "#obs", n
xsd = 1.0_dp
ar = [0.5_dp, 0.3_dp]
print fmt_cr, "AR", ar
print fmt_cr, "true ACF", ar_2_acf(ar(1), ar(2), nacf)
do iter=1, niter
   x = ar_sim(n, ar, xsd)
   print fmt_cr, "est. ACF", acf(x, nacf)
end do
end program xar_2_sim