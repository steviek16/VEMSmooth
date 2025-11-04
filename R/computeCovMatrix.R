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
