test_that("checkConvergence behaves correctly under typical differences", {
  expect_false(checkConvergence(10, 1, 0.01))
  expect_true(checkConvergence(1, 1.001, 0.01))
  expect_false(checkConvergence(NA, 1, 0.01))
  expect_false(checkConvergence(NULL, 1, 0.01))
})
