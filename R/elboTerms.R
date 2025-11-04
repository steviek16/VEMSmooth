#' elboTerms - ELBO component functions for VEMSmooth
#'
#' @description
#' Internal helper functions for computing the Evidence Lower Bound (ELBO) components within the Variational EM algorithm for correlated functional data.
#'
#' Each function corresponds to a specific ELBO sub-term:
#'  - expectedLogLikelihood(): expected log-likelihood for each curve.
#'  - elboInclusionTerm(): contribution from inclusion indicators (Z terms).
#'  - elboBetaTerm(): contribution from β coefficients.
#'  - elboThetaTerm(): contribution from θ prior parameters.
#'  - elboSigmaTerm(): contribution from σ² prior.
#'  - elboTauTerm(): contribution from τ² prior.
#'
#' All functions are internal and are not meant to be called directly by users.
#' @keywords internal
#'
#' @title expectedLogLikelihood
#' @description
#' Computes the expected log-likelihood term for curve i given correlated errors.
#'
#' @param y List of observed data for each curve.
#' @param ni Vector of sample sizes per curve.
#' @param B List of basis matrices.
#' @param i Curve index.
#' @param iter Current iteration index.
#' @param delta1_q Shape parameter of sigma²’s variational posterior.
#' @param delta2_values Rate parameters of sigma²’s variational posterior.
#' @param mu_beta_values Posterior mean matrix of beta.
#' @param Sigma_beta Posterior covariance array of beta.
#' @param p_values Matrix of inclusion probabilities.
#' @param psi Correlation matrix (Ψ).
#'
#' @return Numeric expected log-likelihood value.
#' @keywords internal

#first main term of the ELBO
expectedLogLikelihood <- function(y, ni, B, i, iter,
                                  delta1_q, delta2_values,
                                  mu_beta_values, Sigma_beta,
                                  p_values, psi) {
  delta2_q <- delta2_values[iter]
  mu_i_q <- mu_beta_values[iter, grep(paste0("^beta_[0-9]+_", i, "$"), colnames(mu_beta_values))]
  sigma2_i_q <- diag(Sigma_beta[,,i])

  # Extract inclusion probabilities
  p_qi <- p_values[iter, grep(paste0("^beta_[0-9]+_", i, "$"), colnames(mu_beta_values))]

  #expected squared residual term
  res_term <-  expectedResidualSq(B, i, y, mu_beta_values, Sigma_beta, p_values, iter, psi)

  res <- (-ni[i] / 2) * (log(2 * pi) + expectedLogSigma2(delta1_q, delta2_q)) -
    0.5 * log(det(psi)) -
    0.5 * expectedInvSigma2(delta1_q, delta2_q) * res_term

  as.numeric(res)
}

#inclusion probability term
elboInclusionTerm <- function(i, iter, K, p_values, mu_beta_values, a1_values, a2_values) {
  p_qi <- p_values[iter, grep(paste0("^beta_[0-9]+_", i, "$"), colnames(mu_beta_values))]
  a1_qi <- a1_values[iter, grep(paste0("^beta_[0-9]+_", i, "$"), colnames(mu_beta_values))]
  a2_qi <- a2_values[iter, grep(paste0("^beta_[0-9]+_", i, "$"), colnames(mu_beta_values))]

  term_sum <- numeric(K)
  for (k in 1:K) {
    resa <- if (p_qi[k] > 0) p_qi[k] * log(p_qi[k]) else 0
    resb <- if (p_qi[k] < 1) (1 - p_qi[k]) * log(1 - p_qi[k]) else 0

    term_sum[k] <- p_qi[k] * expectedLogTheta(a1_qi[k], a2_qi[k]) - resa +
      (1 - p_qi[k]) * expectedLogOneMinusTheta(a1_qi[k], a2_qi[k]) - resb
  }

  sum(term_sum)
}

#beta coefficient term
elboBetaTerm <- function(i, iter, K,
                         delta1_q, delta2_values,
                         lambda1_q, lambda2_values,
                         mu_beta_values, Sigma_beta) {
  delta2_q <- delta2_values[iter]
  lambda2_q <- lambda2_values[iter]

  mu_i_q <- mu_beta_values[iter, grep(paste0("^beta_[0-9]+_", i, "$"), colnames(mu_beta_values))]
  Sigma_i_q <- Sigma_beta[,,i]
  sigma2_i_q <- diag(Sigma_i_q)

  res <- -(K / 2) * (expectedLogSigma2(delta1_q, delta2_q) + expectedLogTau2(lambda1_q, lambda2_q)) -
    sum((sigma2_i_q + mu_i_q^2) / 2) *
    expectedInvSigma2(delta1_q, delta2_q) * expectedInvTau2(lambda1_q, lambda2_q) -
    (-0.5 * log(det(Sigma_i_q)) - 0.5 * K)

  as.numeric(res)
}

#Theta (Beta prior) term
elboThetaTerm <- function(i, iter, K, a1_values, a2_values, mu_beta_values, mu_ki) {
  a1_qi <- a1_values[iter, grep(paste0("^beta_[0-9]+_", i, "$"), colnames(mu_beta_values))]
  a2_qi <- a2_values[iter, grep(paste0("^beta_[0-9]+_", i, "$"), colnames(mu_beta_values))]

  term_sum <- numeric(K)
  for (k in 1:K) {
    term_sum[k] <- (mu_ki - a1_qi[k]) * expectedLogTheta(a1_qi[k], a2_qi[k]) +
      (1 - mu_ki - a2_qi[k]) * expectedLogOneMinusTheta(a1_qi[k], a2_qi[k]) +
      log(gamma(a1_qi[k]) * gamma(a2_qi[k]))
  }

  sum(term_sum)
}

#variance term
elboSigmaTerm <- function(iter, delta_1, delta_2, delta1_q, delta2_values) {
  delta2_q <- delta2_values[iter]

  res <- -delta1_q * log(delta2_q) +
    (delta1_q - delta_1) * expectedLogSigma2(delta1_q, delta2_q) +
    (delta2_q - delta_2) * expectedInvSigma2(delta1_q, delta2_q)

  as.numeric(res)
}

#tau variance term
elboTauTerm <- function(iter, lambda_1, lambda_2, lambda1_q, lambda2_values) {
  lambda2_q <- lambda2_values[iter]

  res <- -lambda1_q * log(lambda2_q) +
    (lambda1_q - lambda_1) * expectedLogTau2(lambda1_q, lambda2_q) +
    (lambda2_q - lambda_2) * expectedInvTau2(lambda1_q, lambda2_q)

  as.numeric(res)
}
