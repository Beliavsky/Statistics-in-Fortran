program xar_sim
! simulate from an AR process and compare the theoretical and empirical autocorrelations
use kind_mod , only: dp
use arma_mod , only: ar_sim, yule_walker
use stats_mod, only: acf, pacf, pacf_from_acf, acf_pacf
implicit none
integer, parameter :: n = 10**6, nacf=5, niter=3
real(kind=dp) :: xsd, ar(nacf), x(n), xacf(nacf), xpacf(nacf)
integer :: iacf, iter
character (len=*), parameter :: fmt_cr="(a20, ':', *(1x, f8.4))", &
   fmt_ci="(a20, ':', *(1x, i8))"
print fmt_ci, "#obs", n
xsd = 1.0_dp
ar = 0.0_dp
ar(:2) = [0.4_dp, 0.2_dp]
print fmt_cr, "AR", ar
do iter=1, niter
   x = ar_sim(n, ar, xsd)
   xacf = acf(x, nacf)
   print*
   print fmt_cr, "est. ACF", xacf
   print fmt_cr, "est. PACF", pacf(x, nacf)
   print fmt_cr, "est. PACF", pacf_from_acf(xacf)
   call acf_pacf(x, xacf, xpacf)
   print*,"called acf_pacf"
   print fmt_cr, "est. ACF", xacf
   print fmt_cr, "est. PACF", xpacf
   do iacf=0, nacf
      print fmt_cr, "est. AR", yule_walker(xacf(:iacf))
   end do
end do
end program xar_sim