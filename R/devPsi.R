#' devPsi - Derivatives of Correlation Matrix with respect to w
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

  # First derivative ∂Ψ/∂w
  dPsi <- outer(Xt, Xp, function(a, b) {
    (-abs(a - b) / range_scale) * exp(-w * abs(a - b) / range_scale)
  })

  # Second derivative ∂²Ψ/∂w²
  d2Psi <- outer(Xt, Xp, function(a, b) {
    (abs(a - b) / range_scale)^2 * exp(-w * abs(a - b) / range_scale)
  })

  list(dPsi = dPsi, d2Psi = d2Psi)
}
