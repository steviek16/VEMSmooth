test_that("expected value helpers return valid numeric results", {
  expect_type(expectedInvTau2(2, 1.5), "double")
  expect_gt(expectedInvTau2(2, 1.5), 0)

  expect_type(expectedInvSigma2(3, 2), "double")
  expect_gt(expectedInvSigma2(3, 2), 0)

  expect_true(is.finite(expectedLogSigma2(2, 1.5)))
  expect_true(is.finite(expectedLogTau2(3, 2)))
})

test_that("expected residuals and beta moments behave as expected", {
  set.seed(123)
  B <- list(matrix(rnorm(20), 5, 4))
  y <- list(rnorm(5))
  mu <- matrix(rnorm(40), nrow = 10)
  Sigma <- array(diag(4), dim = c(4, 4, 1))
  p <- matrix(runif(40, 0.1, 0.9), nrow = 10)
  psi <- computePsiMatrix(1:5, 1:5, w = 1)

  colnames(mu) <- colnames(p) <- paste0("beta_", 1:4, "_1")

  val <- expectedResidualSq(B, 1, y, mu, Sigma, p, 5, psi)
  expect_true(is.numeric(val) && val > 0)

  sumBeta <- expectedSumBetaSq(1, mu, Sigma, 5)
  expect_true(sumBeta > 0)

  betaSq <- expectedBetaSq(1, 1, runif(4), mu, Sigma, B, y, 2, 4, 5, psi)
  expect_true(is.numeric(betaSq) && betaSq > 0)
})
