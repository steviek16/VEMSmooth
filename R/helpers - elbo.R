#' helpers - elbo: ELBO component functions for VEMSmooth
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
                                  prob_values, psi) {

  #sigma2 variance scale parameter
  delta2_q <- delta2_values[iter]

  #CHANGE - log det psi only be calculated once
  logdet_psi <- as.numeric(determinant(psi, logarithm = TRUE)$modulus)

  # Extract inclusion probabilities
  p_qi <- prob_values[iter, grep(paste0("^beta_[0-9]+_", i, "$"), colnames(mu_beta_values))]

  #expected squared residual term
  e_res <- expectedResidualSq(
    B   = B,
    i   = i,
    y   = y,
    mu  = mu_beta_values,
    Sigma = Sigma_beta,
    prob  = prob_values,
    iter = iter,
    psi  = psi
  )

  res <- (-ni[i] / 2) * (log(2 * pi) + expectedLogSigma2(delta1_q, delta2_q)) -
    0.5 * logdet_psi -
    0.5 * expectedInvSigma2(delta1_q, delta2_q) * e_res

  as.numeric(res)
}

#inclusion probability term
elboInclusionTerm <- function(i, iter, K, prob_values, mu_beta_values, a1_values, a2_values) {

  idx <- grep(paste0("^beta_[0-9]+_", i, "$"), colnames(mu_beta_values))

  p_qi <- prob_values[iter, idx]
  a1_qi <- a1_values[iter, idx]
  a2_qi <- a2_values[iter, idx]

  K_i <- length(idx)
  out <- numeric(K_i)

  for (k in 1:K) {

    p_k <- p_qi[k]

    # entropy of Bernoulli q(z_ki)
    resa <- if (p_k > 0) p_k * log(p_k) else 0
    resb <- if (p_k < 1) (1 - p_k) * log(1 - p_k) else 0

    out[k] <- p_k  * expectedLogTheta(a1_qi[k], a2_qi[k]) -
      resa +
      (1 - p_k) * expectedLogOneMinusTheta(a1_qi[k], a2_qi[k]) -
      resb
  }

  sum(out)
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
    (-0.5 * as.numeric(determinant(Sigma_i_q, logarithm=TRUE)$modulus) - 0.5 * K)

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

#full elbo for corelated functional VEM
elbo_corr <-  function(
    y, B, ni, m, K, iter,
    delta_1, delta_2, lambda_1, lambda_2,
    delta1_q, delta2_values,
    mu_beta_values, lambda1_q, lambda2_values,
    a1_values, a2_values,
    Sigma_beta, prob_values,
    mu_ki, psi
) {

  # expected log likelihood term
  term_loglik <- sum(
    sapply(1:m, function(i)
      expectedLogLikelihood(
        y = y, ni = ni, B = B, i = i, iter = iter,
        delta1_q = delta1_q, delta2_values = delta2_values,
        mu_beta_values = mu_beta_values,
        Sigma_beta = Sigma_beta,
        prob_values = prob_values,
        psi = psi
      )
    )
  )

  # 2. inclusion (Z) indicator term
  term_z <- sum(
    sapply(1:m, function(i)
      elboInclusionTerm(
        i = i, iter = iter, K = K,
        prob_values = prob_values,
        mu_beta_values = mu_beta_values,
        a1_values = a1_values, a2_values = a2_values
      )
    )
  )

  # beta coefficient term
  term_beta <- sum(
    sapply(1:m, function(i)
      elboBetaTerm(
        i = i, iter = iter, K = K,
        delta1_q = delta1_q, delta2_values = delta2_values,
        lambda1_q = lambda1_q, lambda2_values = lambda2_values,
        mu_beta_values = mu_beta_values,
        Sigma_beta = Sigma_beta
      )
    )
  )

  # theta prior term
  term_theta <- sum(
    sapply(1:m, function(i)
      elboThetaTerm(
        i = i, iter = iter, K = K,
        a1_values = a1_values, a2_values = a2_values,
        mu_beta_values = mu_beta_values,
        mu_ki = mu_ki
      )
    )
  )

  # sigma term
  term_sigma <- elboSigmaTerm(
    iter = iter, delta_1 = delta_1, delta_2 = delta_2,
    delta1_q = delta1_q, delta2_values = delta2_values
  )

  # 6. Tau² term
  term_tau <- elboTauTerm(
    iter = iter, lambda_1 = lambda_1, lambda_2 = lambda_2,
    lambda1_q = lambda1_q, lambda2_values = lambda2_values
  )

  # summing up to arrive at total ELBO
  total <- term_loglik + term_z + term_beta + term_theta + term_sigma + term_tau
  as.numeric(total)
}

# elbo_omega - ELBO as function of w only (used during M-step)
elbo_omega <- function(
    w, y, Xt, B, ni, m, K, iter,
    delta_1, delta_2, lambda_1, lambda_2,
    delta1_q, delta2_values,
    mu_beta_values, lambda1_q, lambda2_values,
    a1_values, a2_values,
    Sigma_beta, prob_values,
    mu_ki
) {

  # compute corr mat
  psi <- computePsiMatrix(Xt, Xt, w)

  # compute full ELBO at this w
  elbo_value <- elbo_corr(
    y = y,
    B = B,
    ni = ni,
    m = m, K = K,
    iter = iter,
    delta_1 = delta_1, delta_2 = delta_2,
    lambda_1 = lambda_1, lambda_2 = lambda_2,
    delta1_q = delta1_q,
    delta2_values = delta2_values,
    mu_beta_values = mu_beta_values,
    lambda1_q = lambda1_q,
    lambda2_values = lambda2_values,
    a1_values = a1_values, a2_values = a2_values,
    Sigma_beta = Sigma_beta,
    prob_values = prob_values,
    mu_ki = mu_ki,
    psi = psi
  )

  return(as.numeric(elbo_value))
}

# dev_elbo - derivative of ELBO wrt w
dev_elbo <- function(
    w, y, Xt, B, ni, m, K, iter,
    delta_1, delta_2, lambda_1, lambda_2,
    delta1_q, delta2_values,
    mu_beta_values, lambda1_q, lambda2_values,
    a1_values, a2_values,
    Sigma_beta, prob_values,
    mu_ki
) {

  # curr sigma2 expectation
  delta2_q <- delta2_values[iter]

  #compute Psi(w) and inverse
  psi <- computePsiMatrix(Xt, Xt, w)
  Psi_inv <- solve(psi)

  #computing derivatives wrt w
  dPsi_list <- devPsi(Xt, Xt, w)
  dPsi  <- dPsi_list$dPsi
  d2Psi <- dPsi_list$d2Psi

  E_inv_sigma2 <- expectedInvSigma2(delta1_q, delta2_q)

  #building A_i for each i
  A_list <- vector("list", m)

  for (i in 1:m) {
    idx_i <- grep(paste0("^beta_[0-9]+_", i, "$"), colnames(mu_beta_values))

    p_i  <- prob_values[iter, idx_i]
    mu_i <- mu_beta_values[iter, idx_i]
    Sig_i <- Sigma_beta[,,i]

    #bernoulli variance
    var_z <- diag(p_i*(1-p_i))

    Var_zi_betai <-var_z * Sig_i +var_z * (mu_i %*% t(mu_i)) +Sig_i * (p_i %*% t(p_i))

    # mean residual
    fitted <- t(p_i * t(B[[i]])) %*% mu_i
    resid  <- y[[i]] - fitted

    A_list[[i]] <- resid %*% t(resid) + B[[i]] %*% Var_zi_betai %*% t(B[[i]])
  }

  A_sum <- Reduce("+", A_list)

  # derivative of ELBO wrt w
  term1 <- -(m / 2) * sum(diag(Psi_inv %*% dPsi))
  term2 <-  0.5 * E_inv_sigma2 * sum(diag(Psi_inv %*% dPsi %*% Psi_inv %*% A_sum))

  dev1 <- term1 + term2
  return(as.numeric(dev1))
}


