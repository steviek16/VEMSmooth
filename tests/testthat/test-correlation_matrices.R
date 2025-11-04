test_that("correlation and covariance matrices behave as expected", {
  Xt <- seq(0, 1, length.out = 5)
  Psi <- computePsiMatrix(Xt, Xt, w = 1)
  Cov <- computeCovMatrix(Xt, Xt, sigma = 2, w = 1)
  devs <- devPsi(Xt, Xt, w = 1)

  expect_equal(dim(Psi), c(5, 5))
  expect_equal(dim(Cov), c(5, 5))
  expect_true(all(Cov >= 0))
  expect_true(all(Psi <= 1))
  expect_true(is.list(devs) && all(c("dPsi", "d2Psi") %in% names(devs)))
  expect_equal(dim(devs$dPsi), c(5, 5))
  expect_equal(dim(devs$d2Psi), c(5, 5))
})
