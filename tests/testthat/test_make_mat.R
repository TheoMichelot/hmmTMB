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

# gam_shell ---------------------------------------------------------------

# The prediction matrix used to come from a gam fitted to a dummy response.
# These check that dropping that fit changes nothing, over the kinds of term
# whose prediction matrix depends on something established at setup: a basis
# and its knots, a factor's levels and contrasts, a by-variable.
test_that("gam_shell gives the same lpmatrix as a fitted gam", {
  set.seed(1)
  data <- data.frame(x = stats::runif(200, 0, 10),
                     z = stats::rnorm(200),
                     g = factor(sample(letters[1:3], 200, replace = TRUE)))
  new_data <- data.frame(x = seq(1, 9, length.out = 25),
                         z = stats::rnorm(25),
                         g = factor(rep(letters[1:3], length.out = 25),
                                    levels = letters[1:3]))

  forms <- list(~ s(x, k = 6, bs = "cs"),
                ~ s(x, k = 8, bs = "tp") + z,
                ~ s(x, k = 5, bs = "ps") + g,
                ~ s(x, k = 6, bs = "cs", by = g),
                ~ x + z + g,
                ~ s(g, bs = "re"))

  for(form in forms) {
    args <- list(formula = stats::update(form, dummy_response ~ .),
                 data = cbind(dummy_response = 1, data))
    G <- do.call(mgcv::gam, c(args, list(fit = FALSE)))
    fitted_gam <- do.call(mgcv::gam, args)

    from_fit <- stats::predict(fitted_gam, newdata = new_data, type = "lpmatrix")
    from_shell <- stats::predict(gam_shell(G), newdata = new_data, type = "lpmatrix")

    expect_equal(as.matrix(from_shell), as.matrix(from_fit),
                 info = deparse(form))
  }
})

test_that("make_matrices builds the same new_data matrices as before", {
  set.seed(1)
  data <- data.frame(x = stats::runif(200, 0, 10),
                     g = factor(sample(letters[1:3], 200, replace = TRUE)))
  new_data <- data.frame(x = seq(1, 9, length.out = 20),
                         g = factor(rep(letters[1:3], length.out = 20),
                                    levels = letters[1:3]))
  forms <- list(a = ~ s(x, k = 6, bs = "cs") + g)

  mats <- make_matrices(forms, data = data, new_data = new_data)

  # The route this replaced, spelled out
  args <- list(formula = stats::update(forms$a, dummy_response ~ .),
               data = cbind(dummy_response = 1, data))
  old <- stats::predict(do.call(mgcv::gam, args), newdata = new_data,
                        type = "lpmatrix")
  G <- do.call(mgcv::gam, c(args, list(fit = FALSE)))

  bare <- function(m) {
    m <- as.matrix(m)
    dimnames(m) <- NULL
    m
  }
  expect_equal(bare(mats$X_fe), bare(old[, 1:G$nsdf, drop = FALSE]))
  expect_equal(bare(mats$X_re), bare(old[, -(1:G$nsdf), drop = FALSE]))
})
