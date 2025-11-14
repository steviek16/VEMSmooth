#' checkConvergence - ELBO Convergence test
#'
#' @description
#' Checks whether the algorithm has converged based on ELBO changes.
#'
#' @param elbo_c Current ELBO value.
#' @param elbo_prev Previous ELBO value.
#' @param convergence_threshold Numeric convergence threshold.
#'
#' @return Logical; TRUE if converged.
#' @keywords internal

checkConvergence <- function(elbo_c, elbo_prev, convergence_threshold) {
  if (is.null(elbo_prev)|| is.na(elbo_prev) || is.null(elbo_c) || is.na(elbo_c)) {
    return(FALSE)
  }

  abs_diff =abs(elbo_c - elbo_prev)
  if (abs_diff <= convergence_threshold) return(TRUE)
  else return(FALSE)
}
