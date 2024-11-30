program test_mean_var_sd
    use kind_mod, only: dp
    use basic_stats_mod, only: mean_var_sd
    implicit none
    real(kind=dp) :: x(3) = [17.0_dp, 19.0_dp, 24.0_dp]
    real(kind=dp) :: mean, variance, stddev
    call mean_var_sd(x, mean, variance, stddev)
    print *, 'Mean = ', mean
    print *, 'Variance = ', variance
    print *, 'Standard Deviation = ', stddev
end program test_mean_var_sd