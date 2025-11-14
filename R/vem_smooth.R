#'vem_smooth - Variational EM Algorithm for Bayesian Basis Selection
#'
#'Runs the Variational EM (VEM) algorithm for smoothing functional curves via basis function selection, accounting for correlated errors.
#'
#'@param y List of length m, each element a numeric vector of observed curve values (ie. outcome variable).
#'@param B List of length m, each element of a n_i x K matrix of basis evaluations (eg. B-spline values).
#'@param m Number of curves.
#'@param Xt Vector of time points (assumed to be identical for all curves).
#'@param K Number of basis functions.
#'@param lambda_1,lambda_2 Hyperparameter for the inverse-gamma prior τ².
#'@param delta_1,delta_2 Hyperparameters for the inverse-gamma prior on σ².
#'@param mu_ki Hyperparameter for the Beta prior on θ (inclusion probability).
#'@param initial_values List of starting values for for w, p, λ₂, and δ₂.
#'@param maxIter Maximum number of VEM iterations.
#'@param convergence_threshold ELBO convergence tolerance.
#'@param lower_opt Lower bound for w optimization.
#'
#'@return A list containing fitted variational parameters, ELBO trajectory, convergence information and final w.
#'@export
#'
#'@details
#'The algorithm alternates between:
#'- **E-step:** updating the variartional distributions for β, σ², τ², Z, and θ.
#'- **M-step:** maximizing the ELBO with respect to the correlation parameter w.

vem_smooth <- function(
  y,
  B,
  Xt = seq(0, 1, length.out=max(length(y))),
  m = length(y),
  K = ncol(B[[1]]),
  mu_ki = 0.5,
  lambda_1 = 1e-10, lambda_2 = 1e-10,
  delta_1 = 1e-10, delta_2 = 1e-10,
  maxIter=1000,
  initial_values,
  convergence_threshold = 0.01,
  lower_opt = 0.1
) {

  ### Initialization setup
  converged <- FALSE
  iter <- 1
  ni <- vapply(y, length, integer(1)) # num obs per curve

    #intialize w and covariance matrix (chosen by user)
  w_decay <- initial_values$w
  Psi_mat <- computePsiMatrix(Xt, Xt, w=w_decay)

    #storage for the matrices for variational parameters / iterative updates
  mu_beta_iter <- matrix(NA, nrow=maxIter, ncol = m*K) #posterior beta means (Norm dist)
  colnames(mu_beta_iter) <- as.vector(outer(paste0("beta_", 1:K, "_"), 1:m, paste0))

  Sigma_beta <- array(NA, dim=c(K, K, m)) #posterior beta covariances (inv gamma) (for each curve)
  prob_inclusion_iter <- matrix(NA, nrow=maxIter, ncol = m*K) # inclusion probabilities
  a1_iter <- a2_iter <- matrix(NA, nrow=maxIter, ncol=m*K) # theta parameters (Beta dist)
  tau2_var_iter <- sigma2_var_iter <- numeric(maxIter) # (rate) parameters for tau, sigma (inv gamma)
  w_trace <- numeric(maxIter) #decay param trace

    #Initialize (chosen by user)
  prob_inclusion_iter[1, ] <- initial_values$p
  tau2_var_iter[1] <- initial_values$lambda2
  sigma2_var_iter[1] <- initial_values$delta2
  w_trace[1] <- w_decay

    # fixed shape parameters (ie. not getting updated in VEM --> are fixed)
  delta1_q <- (sum(ni)+ m*K + 2 * delta_1)/2
  lambda1_q <- (m*K+2*lambda_1)/2

    #initialize ELBO tracking
  elbo_prev <- 0
  ELBO_values <- -Inf

  ### MAIN VEM LOOP

  while (!converged && iter < maxIter) {
    iter <- iter +1

    # E-step
    #Step 1 - Update q(Bi)

    #CHANGE - precompute psi_mat since common to all curves
    Psi_inv <- solve(Psi_mat)

    for (i in 1:m) {
      #inclusion probabilities for curve i from prev iter
      idx_i <- grep(paste0("^beta_[0-9]+_", i, "$"), colnames(mu_beta_iter))
      prob_i   <- prob_inclusion_iter[iter - 1, idx_i]

      #Expected basis contributions under q(zi)

      #expected basis design
      Gi <- prob_i * t(B[[i]])
      GPsiGt <- Gi %*% Psi_inv %*% t(Gi)
      yPsiG <- t(y[[i]]) %*% Psi_inv %*% t(Gi)

      #Expected inverse variances: under variational posteriors
      E_tau_inv <- expectedInvTau2(lambda1_q, tau2_var_iter[iter-1])
      E_sigma_inv <- expectedInvSigma2(delta1_q, sigma2_var_iter[iter-1])

      #variational covar of Bi
      V_i <- solve(diag(E_tau_inv, K)+GPsiGt)

      #full covariance and mean of q(Bi)
      Sigma_bi <- (1/E_sigma_inv) * V_i
      mu_bi <- (E_sigma_inv * yPsiG) %*% Sigma_bi

      #storing updated variational parameters for curve i
      mu_beta_iter[iter, idx_i] <- mu_bi
      Sigma_beta[,,i] <- Sigma_bi
    }

    #Step 2 - update q(sigma) - scale parameter for delta2_q (stored in sigma2_var_iter)

    sigma2_var_iter[iter] <- (
      sum(sapply(1:m, function(i){
        expectedResidualSq(B=B, i=i, y=y, mu=mu_beta_iter, Sigma=Sigma_beta, prob=prob_inclusion_iter, iter=iter, psi=Psi_mat)
      })) +
      sum(sapply(1:m, function(i){
        expectedSumBetaSq(i, mu_beta_iter, Sigma_beta, iter)
      })) * expectedInvTau2(lambda1_q, tau2_var_iter[iter-1]) + 2 * delta_2
    ) / 2

    #Step 3 - update q(tau): scale parameter for lambda2_q (stored in tau2_var_iter)

    tau2_var_iter[iter] <- (
      sum(sapply(1:m, function(i){
        expectedSumBetaSq(i, mu_beta_iter, Sigma_beta, iter)
      })) * expectedInvSigma2(delta1_q, sigma2_var_iter[iter]) + 2 * lambda_2
    ) / 2

    #Step 4 - update q(Z_Ki) and q(theta_Ki)
    for (i in 1:m) {
      idx_i <- grep(paste0("^beta_[0-9]+_", i, "$"), colnames(mu_beta_iter))

      # probabilities from previous iteration
      prob_prev <- prob_inclusion_iter[iter - 1, idx_i]

      prob_curr <-  prob_inclusion_iter[iter, idx_i]

      #fill NAs in curr row with previous iteration values
      na_idx <- is.na(prob_curr)
      prob_curr[na_idx] <- prob_prev[na_idx]

      #storage for updated values at this iter
      prob_new <- numeric(K)
      a1_new <- numeric(K)
      a2_new <- numeric(K)


      for (k in 1:K) {
        #variational parameters for theta_ki
        a1_q <- prob_prev[k] + mu_ki
        a2_q <- 2 - prob_prev[k] - mu_ki

        # unnormalized log posterior for Zki
        log_rho <- sapply(0:1, function(z) {
          (-ni[i] / 2) * expectedLogSigma2(delta1_q, sigma2_var_iter[iter]) -
            0.5 * expectedInvSigma2(delta1_q, sigma2_var_iter[iter]) *
            expectedBetaSq(z, i, prob_prev, mu_beta_iter, Sigma_beta, B, y, k, K, iter, Psi_mat) +
            z * expectedLogTheta(a1_q, a2_q) +
            (1-z) * expectedLogOneMinusTheta(a1_q, a2_q)
        })

        # convert log_rho for z=0,1 into prob p(Zki)
        denom <- sum(exp(log_rho))
        #fallback if numerical issues
       prob_ki_1 <- ifelse(denom==0||is.infinite(denom),
                           c(0,1)[which.max(log_rho)],
                           exp(log_rho[2])/denom)

        #coordinate ascent update
       prob_curr[k] <- prob_ki_1

       #storing output
       prob_new[k] <- prob_ki_1
        a1_new[k] <- a1_q
        a2_new[k] <- a2_q
    }

      #assigning updated values
      prob_inclusion_iter[iter, idx_i] <-prob_curr
      a1_iter[iter, idx_i] <- a1_new
      a2_iter[iter, idx_i] <- a2_new
    }

    #Step 5 - M-Step (Optimize w)
    w_decay_opt <- tryCatch(
      optim(
        par = w_decay,
        fn = elbo_omega,
        gr = dev_elbo,
        y = y, Xt = Xt, B = B, ni = ni, m = m, K = K, iter = iter,
        delta_1 = delta_1, delta_2 = delta_2, lambda_1 = lambda_1, lambda_2 = lambda_2,
        delta1_q = delta1_q, delta2_values = sigma2_var_iter,
        mu_beta_values = mu_beta_iter, lambda1_q = lambda1_q,
        lambda2_values = tau2_var_iter, a1_values = a1_iter, a2_values = a2_iter,
        Sigma_beta = Sigma_beta, prob_values = prob_inclusion_iter, mu_ki = mu_ki,
        control = list(fnscale = -1), method = "L-BFGS-B",
        lower = lower_opt, upper = 1e10
      )$par,
      error = function(e) NA_real_)

    #CHANGE

    #if update fails --> skip the update
    if (!is.na(w_decay_opt)) {
      w_decay <- w_decay_opt
      Psi_mat <- computePsiMatrix(Xt, Xt, w=w_decay)
      w_trace[iter] <- w_decay
    }

    #Step 6 Compute ELBO
    elbo_curr <-  elbo_corr(
      y, B, ni, m, K, iter, delta_1, delta_2, lambda_1, lambda_2,
      delta1_q, sigma2_var_iter, mu_beta_iter, lambda1_q, tau2_var_iter,
      a1_iter, a2_iter, Sigma_beta, prob_inclusion_iter, mu_ki, Psi_mat
    )

    ELBO_values <- c(ELBO_values, elbo_curr)

    converged <- checkConvergence(elbo_curr, elbo_prev, convergence_threshold)
    elbo_prev <- elbo_curr
  }

  #OUTPUT
  result <- list(
    data = y, B = B,
    mu_beta = mu_beta_iter[iter, ], Sigma_beta=Sigma_beta,
    prob=prob_inclusion_iter[iter, ],
    lambda1=lambda1_q, lambda2 = tau2_var_iter[iter],
    delta1=delta1_q, delta2=sigma2_var_iter[iter],
    a1=a1_iter[iter, ], a2 = a2_iter[iter, ],
    elbo_values = ELBO_values, n_iterations = iter,
    converged =converged, w=w_decay, w_trace=w_trace[1:iter]
  )

  return(result)
}
