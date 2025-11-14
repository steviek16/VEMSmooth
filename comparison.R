rm(list = ls())
set.seed(1285)

library(fda)
library(dqrng)
library(mvnfast)

## --- 1. Load original code into its own environment -----------------

source("/Users/stephenkinsey/THESIS/VEMSmooth/aux_fun_corr.R")
source("/Users/stephenkinsey/THESIS/VEMSmooth/elbo_formulas_corr.R")
source("/Users/stephenkinsey/THESIS/VEMSmooth/VB_corr_vem.R")
source("/Users/stephenkinsey/THESIS/VEMSmooth/sim_data_corr.R")

## --- 2. Simulate data ONCE using original sim_data_corr ------------

m     <- 5
K     <- 10
times <- 100
coef  <- c(-2,0,1.5,1.5,0,-1,-0.5,-1,0,0)
sigma <- 0.05

data <- sim_data_corr(
  m      = m,
  times  = times,
  ordem  = 4,
  K      = K,
  coef   = coef,
  seed   = 1285,
  sigma  = sigma,
  w      = 10,
  basis_type = "B-splines"
)

## --- 3. Initial values (same for both original + package) ----------

p_initial <- rep(1, K)
initial_values <- list(
  p       = rep(p_initial, m),
  delta2  = 5,
  lambda2 = 10000,
  w       = 10
)

lambda_1 <- 1e-6
lambda_2 <- 1e-6
delta_1  <- 10
delta_2  <- 9 * 0.0025
maxIter  <- 500
conv_thr <- 0.01
lower_opt <- 0.1

## --- 4. Run ORIGINAL vb_bs_corr (original function name) -----------

out_orig <- vb_bs_corr(
  y        = data$y,
  B        = data$B,
  m        = m,
  mu_ki    = 0.5,
  lambda_1 = lambda_1,
  lambda_2 = lambda_2,
  delta_1  = delta_1,
  delta_2  = delta_2,
  maxIter  = maxIter,
  K        = K,
  initial_values        = initial_values,
  convergence_threshold = conv_thr,
  Xt       = data$Xt,
  lower_opt = lower_opt
)

str(out_orig, max.level = 1)

## --- 5. Load your package and run vem_smooth -----------------------

# from package root during development:
# devtools::load_all(".")

library(VEMSmooth)  # if already installed

out_pkg <- vem_smooth(
  y        = data$y,
  B        = data$B,
  Xt       = data$Xt,
  m        = m,
  K        = K,
  mu_ki    = 0.5,
  lambda_1 = lambda_1,
  lambda_2 = lambda_2,
  delta_1  = delta_1,
  delta_2  = delta_2,
  maxIter  = maxIter,
  initial_values        = initial_values,
  convergence_threshold = conv_thr,
  lower_opt             = lower_opt
)

str(out_pkg, max.level = 1)

## Map package output to original-style names for comparison

mu_beta_orig <- out_orig$mu_beta
p_orig       <- out_orig$p

mu_beta_pkg  <- out_pkg$mu_beta
p_pkg        <- if (!is.null(out_pkg$prob)) out_pkg$prob else out_pkg$p  # adjust if needed

plot(out_orig$elbo_values, type = "l", lwd = 2,
     main = "ELBO: original vb_bs_corr vs package vem_smooth",
     ylab = "ELBO", xlab = "iteration")
lines(out_pkg$elbo_values, col = 2, lwd = 2, lty = 2)
legend("bottomright",
       legend = c("original vb_bs_corr", "package vem_smooth"),
       col    = c(1, 2), lty = c(1, 2), lwd = 2)

all.equal(mu_beta_orig, mu_beta_pkg, tolerance = 1e-6)
all.equal(p_orig,       p_pkg,       tolerance = 1e-6)

max(abs(mu_beta_orig - mu_beta_pkg))
max(abs(p_orig       - p_pkg))

seq_values <- lapply(seq(1, m*K, by = K), function(x){ seq(x, x+K-1) })

cat("Original selection:\n")
for(i in 1:m){
  cat('curve', i, ':',
      ifelse(p_orig[seq_values[[i]]] > 0.50, 1, 0), '\n')
}
cat('truth :', ifelse(abs(coef) > 0, 1, 0), "\n\n")

cat("Package selection:\n")
for(i in 1:m){
  cat('curve', i, ':',
      ifelse(p_pkg[seq_values[[i]]] > 0.50, 1, 0), '\n')
}
cat('truth :', ifelse(abs(coef) > 0, 1, 0), "\n")
