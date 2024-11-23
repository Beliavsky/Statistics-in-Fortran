module linear_regressor_mod
use kind_mod
use linear_algebra_mod, only: solve_linear_system
implicit none
private
public :: LinearRegressor_fit, LinearRegressor_predict

contains

subroutine LinearRegressor_fit(X, y, params)
real(dp), intent(in) :: X(:,:) ! LHS matrix of predictors
real(dp), intent(in) :: y(:)   ! RMS matrix of responses
real(dp), allocatable, intent(out) :: params(:) ! computed regression coefficients
real(dp), allocatable :: X_aug(:,:), XtX(:,:), Xty(:)
integer :: n_samples, n_features
logical :: success
n_samples = size(X, 1)
n_features = size(X, 2) + 1  ! Adding bias term
allocate(X_aug(n_samples, n_features))
X_aug(:, 1) = 1.0_dp  ! Bias term
X_aug(:, 2:) = X
XtX = matmul(transpose(X_aug), X_aug)
Xty = matmul(transpose(X_aug), y)
allocate(params(n_features))
call solve_linear_system(XtX, Xty, params, success)
if (.not. success) then
  print *, "Error solving linear system in LinearRegressor_fit."
  stop
end if
end subroutine LinearRegressor_fit

subroutine LinearRegressor_predict(X, params, y_pred)
! given data for predictors and regression coefficients, computed predictions
real(dp), intent(in) :: X(:,:)    ! predictor data
real(dp), intent(in) :: params(:) ! regression coefficients, with the intercept first
real(dp), allocatable, intent(out) :: y_pred(:) ! regression predictions
real(dp), allocatable :: X_aug(:,:)
integer :: n_samples, n_features
n_samples = size(X, 1)
n_features = size(params)
allocate(X_aug(n_samples, n_features))
X_aug(:, 1) = 1.0_dp  ! Bias term
X_aug(:, 2:) = X
y_pred = matmul(X_aug, params)
end subroutine LinearRegressor_predict

end module linear_regressor_mod
