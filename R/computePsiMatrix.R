#'computePsiMatrix - Correlation (Covariance) Matrix
#'
#'@description
#'Computes Ψ, the correlation matrix, using an exponential kernel.
#'
#'@param Xt Numeric vector of time points.
#'@param Xp Numeric vector of time points (often identical to Xt).
#'@param w Numeric scalar decay parameter, controlling smoothness.
#'
#'@return A covariance matrix of dimension length(Xt) x length(Xp)
#'@keywords internal

computePsiMatrix <-  function(Xt, Xp, w=1) {
  range_scale = max(Xt) - min(Xt)
  Psi <- exp(-w*abs(outer(Xt, Xp, "-"))) / range_scale
  return (Psi)
}
