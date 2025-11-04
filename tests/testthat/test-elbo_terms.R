#tests that ELBO components return finite scalars, etc
test_that("ELBO component functions return finite numeric values", {
  set.seed(42)
  m <- 2; K <- 3
  y <- list(rnorm(5), rnorm(5))
  B <- list(matrix(rnorm(15), 5, 3), matrix(rnorm(15), 5, 3))
  mu <- matrix(rnorm(10 * (m * K)), nrow = 10)
  colnames(mu) <- as.vector(outer(paste0("beta_", 1:K, "_"), 1:m, paste0))
  Sigma <- array(diag(K), dim = c(K, K, m))
  p <- matrix(runif(10 * (m * K), 0.1, 0.9), nrow = 10)
  a1 <- a2 <- matrix(runif(10 * (m * K), 1, 3), nrow = 10)
  psi <- computePsiMatrix(1:5, 1:5, 1)

  val1 <- expectedLogLikelihood(y, c(5, 5), B, 1, 5, 2, rep(1.5, 10), mu, Sigma, p, psi)
  expect_true(is.finite(val1))

  incTerm <- elboInclusionTerm(1, 5, K, p, mu, a1, a2)
  expect_true(is.finite(incTerm))

  betaTerm <- elboBetaTerm(1, 5, K, 2, rep(1.5, 10), 2, rep(1.5, 10), mu, Sigma)
  expect_true(is.finite(betaTerm))

  thetaTerm <- elboThetaTerm(1, 5, K, a1, a2, mu, 0.5)
  expect_true(is.finite(thetaTerm))

  sigmaTerm <- elboSigmaTerm(5, 1e-3, 1e-3, 2, rep(1.5, 10))
  tauTerm <- elboTauTerm(5, 1e-3, 1e-3, 2, rep(1.5, 10))

  expect_true(is.numeric(sigmaTerm) && is.numeric(tauTerm))
})
