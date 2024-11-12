module linear_algebra_mod
  use kind_mod, only: dp
  implicit none
  private
  public :: solve_linear_system

contains

  subroutine solve_linear_system(A, b, x, success)
    ! Solves the linear system A * x = b using LU decomposition (Doolittle's method).
    ! Inputs:
    !   A(n,n) - Coefficient matrix
    !   b(n)   - Right-hand side vector
    ! Outputs:
    !   x(n)   - Solution vector
    !   success - Logical flag indicating if the solution was successful
    ! Assumes that A is a square, non-singular matrix.

    real(kind=dp), intent(in) :: A(:,:), b(:)
    real(kind=dp), intent(out) :: x(:)
    logical, intent(out) :: success
    integer :: n, i, j, k
    real(kind=dp), allocatable :: L(:,:), U(:,:), y(:)
    real(kind=dp) :: sum

    n = size(b)
    allocate(L(n,n), U(n,n), y(n))
    L = 0.0_dp
    U = 0.0_dp
    y = 0.0_dp
    success = .true.

    ! LU Decomposition without pivoting
    do i = 1, n
      ! Compute U(i,k)
      do k = i, n
        sum = 0.0_dp
        do j = 1, i-1
          sum = sum + L(i,j) * U(j,k)
        end do
        U(i,k) = A(i,k) - sum
      end do

      ! Check for zero pivot
      if (U(i,i) == 0.0_dp) then
        success = .false.
        deallocate(L, U, y)
        return
      end if

      ! Compute L(k,i)
      L(i,i) = 1.0_dp
      do k = i+1, n
        sum = 0.0_dp
        do j = 1, i-1
          sum = sum + L(k,j) * U(j,i)
        end do
        L(k,i) = (A(k,i) - sum) / U(i,i)
      end do
    end do

    ! Forward substitution to solve L * y = b
    do i = 1, n
      sum = 0.0_dp
      do j = 1, i-1
        sum = sum + L(i,j) * y(j)
      end do
      y(i) = b(i) - sum
    end do

    ! Backward substitution to solve U * x = y
    do i = n, 1, -1
      sum = 0.0_dp
      do j = i+1, n
        sum = sum + U(i,j) * x(j)
      end do
      x(i) = (y(i) - sum) / U(i,i)
    end do

    deallocate(L, U, y)
  end subroutine solve_linear_system

end module linear_algebra_mod
