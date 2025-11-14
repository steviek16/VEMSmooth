#'helpers - psi
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
  dist_mat <- abs(outer(Xt, Xp, "-"))
  Psi <- exp(-w*dist_mat/ range_scale)
  return (Psi)
}

#'computeCovMatrix - Gaussian Process Covariance Matrix
#'
#'@description
#'Generates a covariance matrix using the specified Gaussian process, scaled by sigma2.
#'
#'@param Xt Numeric vector of time points.
#'@param Xp Numeric vector of time points.
#'@param sigma Numeric variance parameter.
#'@param w Numeric decay parameter.
#'
#'@return A covariance matrix of dimension length(Xt) x length(Xp)
#'@keywords internal

computeCovMatrix <- function(Xt, Xp, sigma=1, w=1) {
  sigma^2*computePsiMatrix(Xt, Xp, w)
}

#' psiDerivatives - Derivatives of Correlation Matrix with respect to w
#'
#' @description
#' Computes the first and second derivatives of the Ψ correlation matrix with respect to the decay parameter w.
#'
#'@param Xt Numeric vector of time points.
#'@param Xp Numeric vector of time points.
#'@param w Numeric decay parameter.
#'
#' @return A list with:
#'   \item{dPsi}{First derivative matrix (∂Ψ/∂w)}
#'   \item{d2Psi}{Second derivative matrix (∂²Ψ/∂w²)}
#' @keywords internal

devPsi <- function(Xt, Xp, w = 1) {
  range_scale <- max(Xt) - min(Xt)
  dist_mat <- abs(outer(Xt, Xp, "-"))

  # First derivative ∂Ψ/∂w
  dPsi <- (-dist_mat/range_scale) * exp(-w*dist_mat/range_scale)

  # Second derivative ∂²Ψ/∂w²
  d2Psi <- (dist_mat/range_scale)^2 * exp(-w*dist_mat/range_scale)

  list(dPsi = dPsi, d2Psi = d2Psi)
}

