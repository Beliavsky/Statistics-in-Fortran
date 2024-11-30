module basic_stats_mod
    use kind_mod, only: dp
    implicit none
    private
    public :: mean_var_sd
contains
    subroutine mean_var_sd(x, mean, variance, stddev)
        real(kind=dp), intent(in) :: x(:)
        real(kind=dp), intent(out) :: mean, variance, stddev
        integer :: n, i
        real(kind=dp) :: m_oldM, m_newM, m_oldS, m_newS, xi
        n = 0
        m_oldM = 0.0_dp
        m_newM = 0.0_dp
        m_oldS = 0.0_dp
        m_newS = 0.0_dp
        do i = 1, size(x)
            xi = x(i)
            n = n + 1
            if (n == 1) then
                m_oldM = xi
                m_newM = xi
                m_oldS = 0.0_dp
            else
                m_newM = m_oldM + (xi - m_oldM) / n
                m_newS = m_oldS + (xi - m_oldM)*(xi - m_newM)
                m_oldM = m_newM
                m_oldS = m_newS
            end if
        end do
        mean = m_newM
        if (n > 1) then
            variance = m_newS / (n - 1)
        else
            variance = 0.0_dp
        end if
        stddev = sqrt(variance)
    end subroutine mean_var_sd
end module basic_stats_mod

