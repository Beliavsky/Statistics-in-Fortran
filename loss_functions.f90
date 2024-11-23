module loss_functions_mod
use kind_mod, only: dp
implicit none
private
public :: square_error_loss, mean_square_error
contains

pure function square_error_loss(y_true, y_pred) result(loss)
! return the squared differences of y_true(:) and y_pred(:)
real(dp), intent(in) :: y_true(:)
real(dp), intent(in) :: y_pred(:)
real(dp), allocatable :: loss(:)
loss = (y_true - y_pred) ** 2
end function square_error_loss

pure function mean_square_error(y_true, y_pred) result(mse)
! return the mean-squared error between y_true(:) and y_pred(:)
real(dp), intent(in) :: y_true(:)
real(dp), intent(in) :: y_pred(:)
real(dp) :: mse
mse = sum((y_true - y_pred) ** 2) / size(y_true)
end function mean_square_error

end module loss_functions_mod
