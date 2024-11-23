module ransac_mod
  use kind_mod, only: dp
  use linear_regressor_mod, only: LinearRegressor_fit, LinearRegressor_predict
  use loss_functions_mod, only: square_error_loss, mean_square_error
  use random_mod, only: random_permutation
  implicit none
  private
  public :: ransac_fit
contains

  subroutine ransac_fit(X, y, n, k, t, d, loss, metric, best_params, best_error)
! Implements the RANSAC algorithm to fit a model robustly in the presence of outliers.
  real(kind=dp), intent(in) :: X(:,:)        ! Input data matrix (n_samples x n_features)
  real(kind=dp), intent(in) :: y(:)          ! Target values corresponding to the input data
  integer, intent(in) :: n                   ! Minimum number of points required to estimate the model
  integer, intent(in) :: k                   ! Maximum number of iterations allowed for RANSAC
  integer, intent(in) :: d                   ! Minimum number of inliers needed for a valid model
  real(kind=dp), intent(in) :: t             ! Threshold to determine whether a point is an inlier
    interface
      function loss(y_true, y_pred) result(loss_vec)
        import :: dp
        real(dp), intent(in) :: y_true(:)
        real(dp), intent(in) :: y_pred(:)
        real(dp), allocatable :: loss_vec(:)
      end function loss

      function metric(y_true, y_pred) result(metric_val)
        import :: dp
        real(dp), intent(in) :: y_true(:)
        real(dp), intent(in) :: y_pred(:)
        real(dp) :: metric_val
      end function metric
    end interface
    real(dp), allocatable, intent(out) :: best_params(:)
    real(dp), intent(out) :: best_error

    integer :: i, n_samples
    integer, allocatable :: ids(:), maybe_inliers(:), inlier_ids(:), inlier_points(:)
    real(dp), allocatable :: errors(:), y_pred(:), params(:)
    real(dp) :: this_error

    n_samples = size(X, 1)
    if (n_samples < n) then
       print*,"in ransac_fit, n, n_samples =", n, n_samples
       stop "need n <= n_samples"
    end if
    best_error = huge(1.0_dp)
    allocate(ids(n_samples))

    do i = 1, k
      call random_permutation(ids)
      maybe_inliers = ids(1:n)
      ! Fit model on maybe_inliers
      call LinearRegressor_fit(X(maybe_inliers, :), y(maybe_inliers), params)

      ! Predict on the rest
      call LinearRegressor_predict(X(ids(n+1:), :), params, y_pred)
      errors = loss(y(ids(n+1:)), y_pred)
      inlier_ids = pack(ids(n+1:), errors < t)

      if (size(inlier_ids) > d) then
        inlier_points = [maybe_inliers, inlier_ids]
        ! Fit better model on inlier_points
        call LinearRegressor_fit(X(inlier_points, :), y(inlier_points), params)
        ! Predict on inlier_points
        call LinearRegressor_predict(X(inlier_points, :), params, y_pred)
        this_error = metric(y(inlier_points), y_pred)
        if (this_error < best_error) then
          best_error = this_error
          best_params = params
        end if
      end if
    end do
  end subroutine ransac_fit

end module ransac_mod
