program xar_pacf
! test ar_pacf
use kind_mod, only: dp
use arma_mod, only: ar_pacf
implicit none
integer, parameter :: nar = 4
real(kind=dp) :: ar(nar)
character (len=*), parameter :: fmt_r = "(*(f8.4))"
ar = [0.5_dp, 0.3_dp, 0.0_dp, -0.5_dp]
print fmt_r, ar
print fmt_r, ar_pacf(ar, 0)
print fmt_r, ar_pacf(ar, 1)
print fmt_r, ar_pacf(ar, nar)
print fmt_r, ar_pacf(ar, 6)
end program xar_pacf
