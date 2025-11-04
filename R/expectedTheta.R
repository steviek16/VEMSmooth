#' expectedTheta - Expected Log-Theta Calculations
#'
#' @description
#' Computes E[log(theta)] and E[log(1-theta)] for Beta distributions.
#'
#' @param a1_ki,a2_ki Shape parameters of the Beta distribution.
#'
#' @return Numeric expected values.
#' @keywords internal

expectedLogTheta <- function(a1_ki, a2_ki) {
  digamma(a1_ki) - digamma(a1_ki+a2_ki)
}

expectedLogOneMinusTheta <- function(a1_ki, a2_ki) {
  digamma(a2_ki) - digamma(a1_ki+a2_ki)
}
