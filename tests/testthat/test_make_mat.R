context("Model matrices")

test_that("each penalty gets its own log-determinant", {
  set.seed(1)
  data <- data.frame(x = stats::runif(200), z = stats::runif(200))

  # Two smooths in one linear predictor. These used to share a single
  # determinant, that of the block diagonal of both penalties, leaving
  # log_det_S shorter than the number of penalties the likelihood loops over.
  mats <- make_matrices(list(a = ~ s(x, k = 6, bs = "cs") + s(z, k = 5, bs = "cs")),
                        data = data)

  expect_equal(ncol(mats$ncol_re), 2)
  expect_length(mats$log_det_S, 2)

  q1 <- diff(mats$ncol_re[, 1]) + 1
  q2 <- diff(mats$ncol_re[, 2]) + 1
  S1 <- as.matrix(mats$S[1:q1, 1:q1])
  S2 <- as.matrix(mats$S[q1 + (1:q2), q1 + (1:q2)])
  expect_equal(mats$log_det_S[1], gdeterminant(S1))
  expect_equal(mats$log_det_S[2], gdeterminant(S2))
})

test_that("penalties are counted across formulas too", {
  set.seed(1)
  data <- data.frame(x = stats::runif(200), z = stats::runif(200))
  mats <- make_matrices(list(a = ~ s(x, k = 6, bs = "cs"),
                             b = ~ s(z, k = 5, bs = "cs")),
                        data = data)

  expect_equal(ncol(mats$ncol_re), 2)
  expect_length(mats$log_det_S, 2)
})

test_that("log_det_S stays numeric when only some formulas have smooths", {
  set.seed(1)
  data <- data.frame(x = stats::runif(200), z = stats::runif(200))
  # The second formula contributes no penalty. An empty penalty list must not
  # turn log_det_S into a list, or TMB cannot read it as a vector.
  mats <- make_matrices(list(a = ~ s(x, k = 6, bs = "cs"), b = ~ z),
                        data = data)

  expect_type(mats$log_det_S, "double")
  expect_length(mats$log_det_S, 1)
})

test_that("a model with no smooths has no determinants", {
  set.seed(1)
  data <- data.frame(x = stats::runif(50))
  mats <- make_matrices(list(a = ~ x), data = data)

  expect_null(mats$ncol_re)
  expect_type(mats$log_det_S, "double")
  expect_length(mats$log_det_S, 0)
})
