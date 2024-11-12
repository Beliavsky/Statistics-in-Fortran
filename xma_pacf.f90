program xma_pacf
! test pacf_ma1
use kind_mod, only: dp
use arma_mod, only: pacf_ma1
implicit none
real(kind=dp) :: ma1
ma1 = 0.9_dp
print "(*(f8.4))", pacf_ma1(ma1, 5)
end program xma_pacf
