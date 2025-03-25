program test_fit
  use kind_mod, only: dp
  use random_mod, only: random_normal
  use prob_dist_mod, only: fit_normal_lognormal
  implicit none

  integer, parameter :: n = 50, niter = 1000
  real(dp) :: x(n), ll_norm(niter), ll_log(niter)
  integer :: i, iter
  logical, parameter :: call_print_stats = .false.
  call random_seed()
  ! Lognormal test
  print *, "Testing with lognormal data"
  do iter=1,niter
  do i = 1, n
    x(i) = exp(random_normal() * 0.5_dp + 1.0_dp)
  end do
  if (call_print_stats) call print_stats(x)
  call fit_normal_lognormal(x, ll_norm(iter), ll_log(iter))
  end do
  print*,"fraction assigned normal:", count(ll_norm > ll_log)/real(niter, kind=dp)

  ! Normal test
  print *, "Testing with normal data"
  do iter=1,niter
  do i = 1, n
    x(i) = random_normal() * 10.0_dp + 100.0_dp
    if (x(i) <= 0.0_dp) x(i) = 1.0_dp
  end do
  if (call_print_stats) call print_stats(x)
  call fit_normal_lognormal(x, ll_norm(iter), ll_log(iter))
  if (call_print_stats) call print_stats(x)
  end do
  print*,"fraction assigned normal:", count(ll_norm > ll_log)/real(niter, kind=dp)

contains
  subroutine print_stats(data)
    real(dp), intent(in) :: data(:)
    integer :: n
    real(dp) :: mean_val, std_val, skew_val, min_val, max_val

    n = size(data)
    mean_val = sum(data) / real(n, dp)
    std_val = sqrt(sum((data - mean_val)**2) / real(n, dp))
    skew_val = sum((data - mean_val)**3) / (real(n, dp) * std_val**3)
    min_val = minval(data)
    max_val = maxval(data)
    print *, "Data statistics:"
    print *, "  Mean =", mean_val
    print *, "  Std Dev =", std_val
    print *, "  Skew =", skew_val
    print *, "  Min =", min_val
    print *, "  Max =", max_val
  end subroutine print_stats

end program test_fit
