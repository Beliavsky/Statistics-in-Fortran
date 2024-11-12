module constants_mod
use kind_mod, only: dp
implicit none
private
public :: pi, pi_over_2, sqrt_two, one_over_sqrt_two_pi, log_two_pi
real(kind=dp), parameter :: pi=3.1415926535897931d0, pi_over_2=1.5707963267948966d0, &
   sqrt_two=1.4142135623730951d0, one_over_sqrt_two_pi=0.3989422804014327d0, &
   log_two_pi=1.8378770664093453d0
end module constants_mod