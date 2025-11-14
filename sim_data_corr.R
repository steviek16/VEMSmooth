# Simulating data (same number of time points for all curves, same coefficients and the basis used for all curves)



#parameter values are provided here --> use coef (for mu_beta), sigma (for var_sigma), w (corr decay, stated in the model)
sim_data_corr <- function(m = 5, times = 100, ordem = 4, K = 10, coef = c(-2,0,1.5,1.5,0,-1,-0.5,-1,0,0), seed = 1234, sigma = 0.1, w = 0.7, basis_type = "B-splines"){

  if(basis_type == "B-splines"){
    Xt <- seq(0, 1, length = times)
    basis_data <- create.bspline.basis(range(Xt), norder = ordem, nbasis = K)
    B_simulated_data <- getbasismatrix(Xt, basis_data, nderiv = 0)
  }  else{ #Fourier
        Xt = seq(0, 2*pi, length = times)
        basis_data <- create.fourier.basis(range(Xt), nbasis = K, dropind=1)
        B_simulated_data <- getbasismatrix(Xt, basis_data, nderiv = 0)
  }

  # Calculate the covariance matrix
  Cov_i <- calCov(Xt, Xt, sigma = sigma, w=w)

  # Generate m error functions from GP
  n.samples <- m
  et <- matrix(rep(0,length(Xt)*n.samples), ncol=n.samples)
  set.seed(seed)
  for (i in 1:n.samples) {
    # Each column represents a sample from a multivariate normal distribution
    # with zero mean and covariance Cov_i
    et[,i] <- mvrnorm(n=1, mu=rep(0, length(Xt)), Sigma=Cov_i)
  }

  # Generating m curves (using B-splines or Fourier)
  argvals <- lapply(1:m, function(s) Xt)
  if(basis_type == "B-splines"){
    y <- lapply(1:m, function(s) as.numeric(B_simulated_data %*% coef) + et[,s])
    B <- lapply(argvals, function(s) B_simulated_data)
  } else{
    y <- lapply(1:m, function(s) cos(Xt) + sin(2*Xt) + et[,s]) # same as Pedro's paper
    B <- lapply(argvals, function(s) B_simulated_data)
  }

  return(list(y = y, B = B, Xt = Xt))

}


# test:
# sim_data_corr(m = 5, times = 100, ordem = 4, K = 10, coef = c(-2,0,1.5,1.5,0,-1,-0.5,-1,0,0), seed = 1285, sigma = 0.05, w = 10, basis_type = "B-splines")
