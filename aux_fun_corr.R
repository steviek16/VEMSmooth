#Auxiliary functions for VEM with correlated errors
# Expected values -----

E_inv_tau2 <- function(lambda1, lambda2){
  lambda1/lambda2
}


E_inv_sigma2 <- function(delta1, delta2){
  delta1/delta2
}

E_log_sigma2 <- function(delta1, delta2){
  log(delta2) - digamma(delta1)
}

E_log_tau2 <- function(lambda1, lambda2){
  log(lambda2) - digamma(lambda1)
}


E_square_gi <- function(B, i, y, mu, Sigma, p, iter, psi){
  p_i <- p[iter - 1, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu))]
  mu_bi <- mu[iter, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu))]
  Sigmai <- Sigma[,,i]

  var_z <- diag(p_i*(1-p_i))

  Var_zi_betai <- var_z*Sigmai + var_z*(mu_bi%*%t(mu_bi)) + Sigmai*(p_i%*%t(p_i))

  res <- t(y[[i]] - B[[i]]%*%(mu_bi*p_i))%*%solve(psi)%*%(y[[i]] - B[[i]]%*%(mu_bi*p_i)) + sum(diag(solve(psi)%*%(B[[i]]%*%Var_zi_betai%*%t(B[[i]]))))
  return(res)
}

E_sumk_beta_i_corr <- function(i, mu, Sigma, iter){
  mu_bi <- mu[iter, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu))]
  var_bi <- diag(Sigma[,,i])

  return(sum(var_bi + mu_bi^2))
}


E_square_beta_i <- function(z, i, p, mu, Sigma, B, y, k, K, iter, psi){

  ni <- unlist(lapply(y,length))

  mu_bi <- mu[iter, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu))]
  Sigmai <- Sigma[,,i]

  # vector with probabilities when k-th position is z
  pi_notk <- p
  pi_notk[k] <- z

  var_zi_notk <- pi_notk*(1-pi_notk)
  var_zi_notk[k] <- 0

  Var_zi_notk <- diag(var_zi_notk)

  Var_zi_betai_notk <- Var_zi_notk*Sigmai + Var_zi_notk*(mu_bi%*%t(mu_bi)) + Sigmai*(pi_notk%*%t(pi_notk))

  Esquare_i_beta <- t(y[[i]] - t(pi_notk*t(B[[i]]))%*%mu_bi)%*%solve(psi)%*%(y[[i]] - t(pi_notk*t(B[[i]]))%*%mu_bi) + sum(diag(solve(psi)%*%(B[[i]]%*%Var_zi_betai_notk%*%t(B[[i]]))))

  return(Esquare_i_beta)

}


E_log_theta_ki <- function(a1_ki, a2_ki){

  return(digamma(a1_ki) - digamma(a1_ki + a2_ki))
}


E_log_theta_ki_c <- function(a1_ki, a2_ki){

  return(digamma(a2_ki) - digamma(a1_ki + a2_ki))
}

# Elbo -----


elbo_corr <- function(y, ni, B, m, K, iter, delta_1, delta_2, lambda_1, lambda_2, delta1_q, delta2_values, mu_beta_values, lambda1_q, lambda2_values, a1_values, a2_values, Sigma_beta, p_values, mu_ki, psi){

  res <- sum(sapply(1:m, function(i){E_log_like_i_corr(y = y, B = B, ni = ni, i = i, iter = iter , delta1_q = delta1_q, delta2_values = delta2_values, mu_beta_values = mu_beta_values, Sigma_beta = Sigma_beta, p_values = p_values, psi = psi)})) + sum(sapply(1:m, function(i){diff_z_i(i = i, iter = iter , K = K, p_values = p_values, mu_beta_values = mu_beta_values, a1_values = a1_values, a2_values = a2_values)})) + sum(sapply(1:m, function(i){diff_beta_i(i = i, iter = iter, K = K, delta1_q = delta1_q, delta2_values = delta2_values, lambda1_q = lambda1_q, lambda2_values = lambda2_values, mu_beta_values = mu_beta_values, Sigma_beta = Sigma_beta)})) + sum(sapply(1:m, function(i){diff_theta_i(i = i, iter = iter, K = K, a1_values = a1_values, a2_values = a2_values, mu_beta_values = mu_beta_values, mu_ki = mu_ki)})) + diff_sigma2(iter = iter, delta_1 = delta_1, delta_2 = delta_2, delta1_q = delta1_q, delta2_values = delta2_values) + diff_tau2(iter = iter, lambda_1 = lambda_1, lambda_2 = lambda_2, lambda1_q = lambda1_q, lambda2_values = lambda2_values)

  return(res)
}

check_convergence <- function(elbo_c, elbo_prev, convergence_threshold) {
  if(is.null(elbo_prev) == TRUE) {
    return(FALSE)
  }
  else{
    dif <- elbo_c - elbo_prev
    if(abs(dif)  <= convergence_threshold) return(TRUE)
    else return(FALSE)
  }
}


calPsi <- function(tt, s, w=1) {
  Psi <- matrix(rep(0, length(tt)*length(s)), nrow=length(tt))
  for (i in 1:nrow(Psi)) {
    for (j in 1:ncol(Psi)) {
      Psi[i,j] <- exp(- w * abs(tt[i] - s[j]) / (max(tt) - min(tt)))
    }
  }
  return(Psi)
}

# Generating correlation matrix based on the specified Gaussian process
calCov <- function(tt, s, sigma=1, w=1) {
  Cov <- sigma^2 * calPsi(tt, s, w)
  return(Cov)
}


# EM ----

# 1st derivative of psi with respect to w
dev_psi <- function(tt, s, w=1) {
  dev_Psi <- matrix(rep(0, length(tt)*length(s)), nrow=length(tt))
  for (i in 1:nrow(dev_Psi)) {
    for (j in 1:ncol(dev_Psi)) {
      dev_Psi[i,j] <- (-abs(tt[i]-s[j]) / (max(tt) - min(tt)))*exp(-w*abs(tt[i]-s[j]) / (max(tt) - min(tt)))
    }
  }
  return(dev_Psi)
}

# 2nd derivative of psi with respect to w
dev2_psi <- function(tt, s, w=1) {
  dev2_Psi <- matrix(rep(0, length(tt)*length(s)), nrow=length(tt))
  for (i in 1:nrow(dev2_Psi)) {
    for (j in 1:ncol(dev2_Psi)) {
      dev2_Psi[i,j] <- (abs(tt[i]-s[j]) / (max(tt) - min(tt)))^2*exp(-w*abs(tt[i]-s[j]) / (max(tt) - min(tt)))
    }
  }
  return(dev2_Psi)
}

# Objective function to be optimized (ELBO as a function of w only)
elbo_omega <- function(w, y, Xt, B, ni, m, K, iter, delta_1, delta_2, lambda_1, lambda_2, delta1_q, delta2_values, mu_beta_values, lambda1_q, lambda2_values, a1_values, a2_values, Sigma_beta, p_values, mu_ki){
   psi <- calPsi(Xt, Xt, w)

   elbo_value <- elbo_corr(y = y, B = B, ni = ni, m = m, K = K, iter = iter, delta_1 = delta_1, delta_2 = delta_2, lambda_1 = lambda_1, lambda_2 = lambda_2, delta1_q = delta1_q, delta2_values = delta2_values, mu_beta_values = mu_beta_values, lambda1_q = lambda1_q, lambda2_values = lambda2_values, a1_values = a1_values, a2_values = a2_values, Sigma_beta = Sigma_beta, p_values = p_values, mu_ki = mu_ki, psi = psi)
   return(elbo_value)
}

# 1st derivative of objective function with respect to w
dev_elbo <- function(w, y, Xt, B, ni, m, K, iter, delta_1, delta_2, lambda_1, lambda_2, delta1_q, delta2_values, mu_beta_values, lambda1_q, lambda2_values, a1_values, a2_values, Sigma_beta, p_values, mu_ki){
  delta2_q = delta2_values[iter]

  d_psi <- dev_psi(Xt, Xt, w)
  d2_psi <- dev2_psi(Xt, Xt, w)
  psi <- calPsi(Xt, Xt, w)

  A_matrix <- lapply(1:m, function(i){
    p_i <- p_values[iter, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu_beta_values))]
    mu_bi <- mu_beta_values[iter, grep(paste0("^beta_[0-9]+_",i,"$"), colnames(mu_beta_values))]
    Sigmai <- Sigma_beta[,,i]

    var_z <- diag(p_i*(1-p_i))

    Var_zi_betai <- var_z*Sigmai + var_z*(mu_bi%*%t(mu_bi)) + Sigmai*(p_i%*%t(p_i))

    (y[[i]] - t(p_i*t(B[[i]]))%*%mu_bi)%*%t(y[[i]] - t(p_i*t(B[[i]]))%*%mu_bi) + B[[i]]%*%Var_zi_betai%*%t(B[[i]])
  })

  dev_E <- E_inv_sigma2(delta1 = delta1_q, delta2 = delta2_q)
  s_psi <- solve(psi)

  dev1 <- -(m / 2) * sum(diag(s_psi %*% d_psi)) +
    0.5 * dev_E * sum(diag(s_psi %*% d_psi %*% s_psi %*% Reduce('+', A_matrix)))


  return(dev1)
}


