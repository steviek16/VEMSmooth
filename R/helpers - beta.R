#' helpers - beta: beta-related expectation helpers
#'
#' @description
#' Expected squared residuals and beta^2 moments under the variational posterior.
#' These are used in the updates of q(sigma^2), q(tau^2), and q(Z_ki).
#'
#' @param B List of basis matrices (B[[i]] is n_i x K).
#' @param i Curve index.
#' @param y List of observations (y[[i]] is length n_i vector).
#' @param mu Matrix of posterior beta means (rows = iterations, cols = beta_k_i).
#' @param Sigma Array of posterior beta covariance matrices (K x K x m).
#' @param prob Matrix of inclusion probabilities (rows = iterations, cols = beta_k_i).
#' @param iter Current iteration index.
#' @param psi Correlation matrix for the errors.
#'
#' @return Numeric expected squared (or summed squared) values.
#' @keywords internal


#e_square_gi in the original code
expectedResidualSq <- function(B, i, y, mu, Sigma, prob, iter, psi) {
  #inclusion probs and beta means for curve i
  prob_i <- prob[iter - 1, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu))]
  mu_bi <- mu[iter, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu))]
  Sigmai <- Sigma[,,i]

  #Var(zi) - diag mat of Bernoulli variances
  var_z <- diag(prob_i*(1-prob_i))

  Var_zi_betai <- var_z*Sigmai + var_z*(mu_bi%*%t(mu_bi)) + Sigmai*(prob_i%*%t(prob_i))

  Psi_inv <- solve(psi)

  #quadratic form in residual mean
  resid <- y[[i]]-B[[i]]%*%(mu_bi * prob_i)
  quad_part <- t(resid)%*%Psi_inv%*%resid

  #trace term for residual variance
  trace_part <- sum(diag(Psi_inv %*% (B[[i]] %*% Var_zi_betai %*% t(B[[i]]))))

  as.numeric(quad_part+trace_part)
}


expectedSumBetaSq <- function(i, mu, Sigma, iter) {
  mu_bi <- mu[iter, grep(paste0("^beta_[0-9]+_", i, "$"), colnames(mu))]
  var_bi <- diag(Sigma[,,i])

  sum(var_bi + mu_bi^2)
}


expectedBetaSq <- function(z, i, prob, mu, Sigma, B, y, k, K, iter, psi) {
  mu_bi <- mu[iter, grep(paste0("^beta_[0-9]+_", i, "$"), colnames(mu))]
  Sigmai <- Sigma[,,i]

  # vector of probabilities w/ k-th element/position is replaced with z
  prob_i_notk <- prob
  prob_i_notk[k] <- z

  var_zi_notk <- prob_i_notk * (1 - prob_i_notk)
  var_zi_notk[k] <- 0
  Var_zi_notk <- diag(var_zi_notk)

  # Use matrix multiplication
  Var_zi_betai_notk <- Var_zi_notk * Sigmai +
    Var_zi_notk * (mu_bi %*% t(mu_bi)) +
    Sigmai * (prob_i_notk %*% t(prob_i_notk))

  Psi_inv <- solve(psi)

  # mean residual given prob_i_notk
  fitted <- t(prob_i_notk*t(B[[i]])) %*% mu_bi
  resid <- y[[i]] - fitted

  quad_part <- t(resid)%*%Psi_inv%*%resid
  trace_part <- sum(diag(Psi_inv%*% (B[[i]] %*% Var_zi_betai_notk %*% t(B[[i]]))))

  as.numeric(quad_part+trace_part)
}
