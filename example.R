rm(list = ls())
set.seed(16)

library(fda)
library(dqrng)
library(mvnfast)
library(VEMSmooth)

source("/Users/stephenkinsey/THESIS/VEMSmooth/aux_fun_corr.R")
source("/Users/stephenkinsey/THESIS/VEMSmooth/elbo_formulas_corr.R")
source("/Users/stephenkinsey/THESIS/VEMSmooth/VB_corr_vem.R")
source("/Users/stephenkinsey/THESIS/VEMSmooth/sim_data_corr.R")

## TOGGLE for Original vs Package
use_package <- TRUE #set true for package

if (!use_package) {
  cat(">>> Using ORIGINAL implementation (vb_bs_corr)\n")
} else {
  cat(">>> Using PACKAGE implementation (vem_smooth from VEMSmooth)\n")
}

# Generating a dataset with 5 curves with 100 evaluation points each using 10 B-splines basis functions
m <- 5 # number of curves
K <- 10 # number of basis functions
times <- 100 # number of evaluation points for each curve (same for all curves)
coef <- c(-2,0,1.5,1.5,0,-1,-0.5,-1,0,0) # Basis coefficients used to generate the curves
sigma <- 0.05 # standard deviation
w_true <- 10
data <- sim_data_corr(m = m, times = times, ordem = 4, K = K, coef = coef, seed = 1285, sigma = sigma, w = w_true, basis_type = "B-splines")

# Run VEM algorithm for only one simulated dataset (in practice we would run for 100 simulations to assess the performance of the algorithm)

# Specifying initial values for the required parameters
p_initial <- rep(1, K) #initialize all prob of inclusion as one for all curves
initial_values <- list(p = rep(p_initial, m),  delta2 = 5, lambda2 = 10000, w = 10)

if (!use_package) {
  out <- vb_bs_corr(y = data$y, B = data$B, m = m, mu_ki = 0.5, lambda_1 = 1e-6, lambda_2 =  1e-6, delta_1 = 10, delta_2 = 9*0.0025, maxIter = 1000, K = K, initial_values = initial_values, convergence_threshold = 0.01, Xt = data$Xt)
} else {
  out <- vem_smooth(y = data$y, B = data$B, m = m, mu_ki = 0.5, lambda_1 = 1e-6, lambda_2 =  1e-6, delta_1 = 10, delta_2 = 9*0.0025, maxIter = 1000, K = K, initial_values = initial_values, convergence_threshold = 0.01, Xt = data$Xt)
}

str(out, max.level=1)

# plot results
# Generating a list to define which values in out$mu_beta and out$p refer to each curve
seq_values <- lapply(c(seq(1, m*K, K)), function(x){seq(x,x+K-1)})

for (i in 1:m) {
  # raw data
  plot(
    data$Xt, data$y[[i]],
    cex = 0.8, pch = 16, col = "grey",
    ylab = expression(g[t]), xlab = expression(t),
    main = paste("Curve", i)
  )

  # posterior inclusion probs: p or prob
  p_vec <- if (!is.null(out$prob)) out$prob else out$p

  # estimated curve using posterior means & hard threshold on p
  beta_hat_i <- out$mu_beta[seq_values[[i]]]
  z_hat_i    <- ifelse(p_vec[seq_values[[i]]] > 0.5, 1, 0)

  fitted_i <- as.numeric((beta_hat_i * z_hat_i) %*% t(data$B[[i]]))

  lines(data$Xt, fitted_i, lwd = 2, col = "red")
}

# selected basis functions for each curve
cat(ifelse(use_package, "\nPackage selection:\n", "\nOriginal selection:\n"))
for(i in 1:m){
  cat('curve',i,':',ifelse(out$p[seq_values[[i]]] > 0.50, 1, 0), '\n')
}
cat('truth: ', ifelse(abs(coef) > 0, 1, 0))

## --- 7. Parameter comparison: beta, w, sigma ------------------------

## 7a. Beta coefficients (truth vs estimates)
# For each basis k, we can look at the estimate from curve 1 and the mean across curves.
beta_mat <- sapply(1:m, function(i) out$mu_beta[seq_values[[i]]]) # K x m
colnames(beta_mat) <- paste0("curve_", 1:m)

beta_mean <- rowMeans(beta_mat)

beta_summary <- data.frame(
  k             = 1:K,
  beta_true     = coef,
  beta_hat_c1   = beta_mat[,1],
  beta_hat_mean = beta_mean
)

cat("---- Beta coefficient comparison ----\n")
print(round(beta_summary, 4))

## 7b. w (decay parameter)
w_hat <- out$w
cat("\n---- Decay parameter w ----\n")
cat("w_true =", w_true, "\n")
cat("w_hat  =", round(w_hat, 4), "\n")

## 7c. sigma (error SD)
# Prior: sigma^2 ~ IG(delta_1, delta_2)
# Variational posterior: IG(delta1_q, delta2_q) = (out$delta1, out$delta2)
# For IG(shape=a, scale=b): E[sigma^2] = b / (a - 1), a > 1
delta1_q <- out$delta1
delta2_q <- out$delta2

sigma2_hat <- delta2_q / (delta1_q - 1)
sigma_hat  <- sqrt(sigma2_hat)

cat("\n---- Error scale sigma ----\n")
cat("sigma_true     =", sigma, "\n")
cat("sigma2_true    =", sigma^2, "\n")
cat("E[sigma^2]_hat ≈", round(sigma2_hat, 6), "\n")
cat("sigma_hat      ≈", round(sigma_hat, 6), "\n")
