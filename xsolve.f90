program test_solve_linear_system
! test solve_linear_system
  use kind_mod, only: dp
  use linear_algebra_mod, only: solve_linear_system
  implicit none
  integer :: n, i
  real(kind=dp), allocatable :: A(:,:), b(:), x(:), Ax(:), errors(:)
  logical :: success
  logical, parameter :: print_sol = .false.
  n = 10**3 ! size of linear system
  print*,"n:", n
  allocate(A(n,n), b(n), x(n))
  call random_seed()
  call random_number(A)
  call random_number(b)
  do i = 1, n
    A(i,i) = A(i,i) + real(n, kind=dp)
  end do
  call solve_linear_system(A, b, x, success)
  if (success) then
    Ax = matmul(A, x)
    errors = Ax - b
    if (print_sol) then
       print *, "Solution x:"
       do i = 1, n
          print *, "x(", i, ") = ", x(i)
       end do
       print *, "Comparison of A * x and b:"
       do i = 1, n
          print *, "A*x(", i, ") = ", Ax(i), "  b(", i, ") = ", b(i), &
               "  Difference = ", errors(i)
       end do
    end if
    print*,"min, max b:", minval(b), maxval(b)
    print*,"min, max Ax-b:", minval(errors), maxval(errors)
    print*,"sum of absolute errors:", sum(abs(errors))
  else
    print *, "Failed to solve the linear system."
  end if
end program test_solve_linear_system
