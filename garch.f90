module garch_mod
  use kind_mod, only: dp
  use random_mod, only: random_normal
  implicit none
  private
  public :: garch_var, gjr_garch_var, qgarch_var, tgarch_var, simulate_tgarch

contains

  function garch_var(var_old, ret, var_const, alpha, beta) result(var_new)
  ! Symmetric GARCH(1,1) model function
    real(kind=dp), intent(in) :: var_old   ! Previous variance
    real(kind=dp), intent(in) :: ret       ! Return at previous time step
    real(kind=dp), intent(in) :: var_const ! Constant term
    real(kind=dp), intent(in) :: alpha     ! ARCH term coefficient
    real(kind=dp), intent(in) :: beta      ! GARCH term coefficient
    real(kind=dp) :: var_new               ! New variance
    var_new = var_const + alpha*ret**2 + beta*var_old
  end function garch_var

  function gjr_garch_var(var_old, ret, var_const, alpha, gamma, beta) result(var_new)
  ! GJR-GARCH model function with asymmetric term
    real(kind=dp), intent(in) :: var_old   ! Previous variance
    real(kind=dp), intent(in) :: ret       ! Return at previous time step
    real(kind=dp), intent(in) :: var_const ! Constant term
    real(kind=dp), intent(in) :: alpha     ! ARCH term coefficient
    real(kind=dp), intent(in) :: gamma     ! Asymmetric term coefficient
    real(kind=dp), intent(in) :: beta      ! GARCH term coefficient
    real(kind=dp) :: var_new               ! New variance

    real(kind=dp) :: ret_sq, coeff_ret_sq
    ret_sq = ret**2
    coeff_ret_sq = merge(alpha + gamma, alpha, ret < 0.0_dp)
    var_new = var_const + coeff_ret_sq*ret_sq + beta*var_old
  end function gjr_garch_var

  function qgarch_var(var_old, ret, var_const, alpha, gamma, beta) result(var_new)
  ! QGARCH model function
    real(kind=dp), intent(in) :: var_old   ! Previous variance
    real(kind=dp), intent(in) :: ret       ! Return at previous time step
    real(kind=dp), intent(in) :: var_const ! Constant term
    real(kind=dp), intent(in) :: alpha     ! ARCH term coefficient
    real(kind=dp), intent(in) :: gamma     ! Quadratic term coefficient
    real(kind=dp), intent(in) :: beta      ! GARCH term coefficient
    real(kind=dp) :: var_new               ! New variance

    ! Compute new variance with quadratic term
    var_new = var_const + alpha*ret**2 + beta*var_old + gamma*ret**3
  end function qgarch_var

  function tgarch_var(var_old, ret, var_const, alpha, gamma, beta, shift) result(var_new)
  ! TGARCH model function with shift parameter
    real(kind=dp), intent(in) :: var_old    ! Previous variance
    real(kind=dp), intent(in) :: ret        ! Return at previous time step
    real(kind=dp), intent(in) :: var_const  ! Constant term
    real(kind=dp), intent(in) :: alpha      ! ARCH term coefficient
    real(kind=dp), intent(in) :: gamma      ! Asymmetric term coefficient
    real(kind=dp), intent(in) :: beta       ! GARCH term coefficient
    real(kind=dp), intent(in) :: shift      ! Shift parameter
    real(kind=dp) :: var_new                ! New variance

    real(kind=dp) :: ret_shifted, ret_shifted_sq, coeff_ret_shifted_sq
    ret_shifted = ret - shift
    ret_shifted_sq = ret_shifted**2
    coeff_ret_shifted_sq = merge(alpha + gamma, alpha, ret_shifted < 0.0_dp)
    var_new = var_const + coeff_ret_shifted_sq*ret_shifted_sq + beta*var_old
  end function tgarch_var

subroutine simulate_tgarch(var_const, alpha, gamma, beta, shift, sigma, returns)
  real(kind=dp), intent(in) :: var_const    ! Constant term in the TGARCH model
  real(kind=dp), intent(in) :: alpha        ! ARCH term coefficient
  real(kind=dp), intent(in) :: gamma        ! Asymmetric term coefficient
  real(kind=dp), intent(in) :: beta         ! GARCH term coefficient
  real(kind=dp), intent(in) :: shift        ! Shift parameter
  real(kind=dp), intent(out) :: sigma(:)    ! Conditional standard deviations
  real(kind=dp), intent(out) :: returns(:)  ! Simulated returns
  integer :: i, n
  real(kind=dp) :: var_ret

  ! Check that sigma and returns have the same size
  n = size(sigma)
  if (size(returns) /= n) then
    print *, "Error: sigma and returns arrays must have the same size."
    stop
  end if

  ! Initialize the process with arbitrary values
  sigma(1) = sqrt(var_const / (1.0_dp - alpha - beta))  ! Initial standard deviation
  var_ret = sigma(1)**2                                 ! Initial variance
  returns(1) = sigma(1) * random_normal()               ! Initial return

  ! Loop to simulate TGARCH process
  do i = 2, n
    var_ret = tgarch_var(var_ret, returns(i-1), var_const, alpha, gamma, beta, shift)
    sigma(i) = sqrt(var_ret)
    returns(i) = sigma(i) * random_normal()
  end do
end subroutine simulate_tgarch

end module garch_mod
