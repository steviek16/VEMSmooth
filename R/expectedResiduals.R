#' expectedResiduals - Excpected Squared Residuals and Beta Moments
#'
#' @description
#' Computes expected residuals and beta2 moments under the variational posterior
#'
#'@param B List of basis matrices.
#'@param i Curve index.
#'@param y List of observations.
#'@param mu Matrix of posterior beta means.
#'@param Sigma Array of posterior beta covariance matrices.
#'@param p Matrix of inclusion probabilities.
#'@param iter Current iteration index.
#'@param psi Correlation matrix.
#'
#'@return Numeric expected squared values.
#'@keywords internal
#'

expectedResidualSq <- function(B, i, y, mu, Sigma, p, iter, psi) {
  p_i <- p[iter - 1, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu))]
  mu_bi <- mu[iter, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu))]
  Sigmai <- Sigma[,,i]

  var_z <- diag(p_i*(1-p_i))

  Var_zi_betai <- var_z*Sigmai + var_z*(mu_bi%*%t(mu_bi)) + Sigmai*(p_i%*%t(p_i))

  res <- t(y[[i]] - B[[i]]%*%(mu_bi*p_i))%*%solve(psi)%*%(y[[i]] - B[[i]]%*%(mu_bi*p_i)) + sum(diag(solve(psi)%*%(B[[i]]%*%Var_zi_betai%*%t(B[[i]]))))
  return(res)
}


expectedSumBetaSq <- function(i, mu, Sigma, iter) {
  mu_bi <- mu[iter, grep(paste0("^beta_[0-9]+_", i, "$"), colnames(mu))]
  var_bi <- diag(Sigma[,,i])
  sum(var_bi + mu_bi^2)
}


expectedBetaSq <- function(z, i, p, mu, Sigma, B, y, k, K, iter, psi) {
  mu_bi <- mu[iter, grep(paste0("^beta_[0-9]+_", i, "$"), colnames(mu))]
  Sigmai <- Sigma[,,i]

  # Replace k-th element with z
  pi_notk <- p
  pi_notk[k] <- z

  var_zi_notk <- pi_notk * (1 - pi_notk)
  var_zi_notk[k] <- 0
  Var_zi_notk <- diag(var_zi_notk)

  # Use matrix multiplication
  Var_zi_betai_notk <- Var_zi_notk * Sigmai +
    Var_zi_notk * (mu_bi %*% t(mu_bi)) +
    Sigmai * (pi_notk %*% t(pi_notk))

  Esquare_i_beta <- t(y[[i]] - t(pi_notk * t(B[[i]])) %*% mu_bi) %*%
    solve(psi) %*% (y[[i]] - t(pi_notk * t(B[[i]])) %*% mu_bi) +
    sum(diag(solve(psi) %*% (B[[i]] %*% Var_zi_betai_notk %*% t(B[[i]]))))

  as.numeric(Esquare_i_beta)
}
