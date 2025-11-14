#' expectedInvTerms - Expected Values for Variance Parameters (sigma and tau)
#'
#' @description
#' Computes the expected inverse and logarithmic terms under the Inverse-gamm distributions for or τ² and σ².
#'
#'@param lambda1,lambda2 Shape and scale parameters for τ².
#'@param delta1,delta2 Shape and scale parameters for σ².
#'
#'@return Numeric expected values.
#'@keywords internal

expectedInvTau2 <- function(lambda1, lambda2) lambda1/lambda2
expectedInvSigma2 <- function(delta1, delta2) delta1/delta2
expectedLogSigma2 <- function(delta1, delta2) log(delta2) - digamma(delta1)
expectedLogTau2 <- function(lambda1, lambda2) log(lambda2) - digamma(lambda1)
