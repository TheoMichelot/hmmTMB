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

# The L convention --------------------------------------------------------

test_that("L is the identity for ordinary smooths", {
  set.seed(1)
  data <- data.frame(x = stats::runif(200), z = stats::runif(200),
                     g = factor(rep(1:5, length.out = 200)))

  one <- make_matrices(list(a = ~ s(x, k = 6, bs = "cs")), data = data)
  expect_equal(one$L, matrix(1, 1, 1))
  expect_equal(one$theta_start, 0)
  expect_equal(one$sp_names, "a.s(x)")   # prefixed by the formula it came from

  two <- make_matrices(list(a = ~ s(x, k = 6, bs = "cs") + s(z, k = 5, bs = "cs")),
                       data = data)
  expect_equal(two$L, diag(2))
  expect_equal(two$theta_start, c(0, 0))

  re <- make_matrices(list(a = ~ s(g, bs = "re")), data = data)
  expect_equal(re$L, matrix(1, 1, 1))
  expect_equal(re$theta_start, 0)
})

test_that("sp_names still names one smoothing parameter per smooth", {
  set.seed(1)
  data <- data.frame(x = stats::runif(200), z = stats::runif(200))
  mats <- make_matrices(list(a = ~ s(x, k = 6, bs = "cs"),
                             b = ~ s(z, k = 5, bs = "cs")),
                        data = data)
  # An ordinary smooth has one penalty and one parameter, so the per-parameter
  # names and the per-penalty column names of ncol_re still coincide
  expect_equal(mats$sp_names, unname(colnames(mats$ncol_re)))
  expect_length(mats$theta_start, 2)
})

test_that("theta_start starts every ordinary smooth at lambda = 1", {
  set.seed(1)
  data <- data.frame(x = stats::runif(200), z = NA)
  s <- 1
  for(i in 1:200) {
    data$z[i] <- stats::rnorm(1, c(0, 4)[s], 1)
    s <- sample(c(s, 3 - s), size = 1, prob = c(0.9, 0.1))
  }
  obs <- Observation$new(data = data, dists = list(z = "norm"), n_states = 2,
                         par = list(z = list(mean = c(0, 4), sd = c(1, 1))),
                         formulas = list(z = list(mean = ~ s(x, k = 6, bs = "cs"),
                                                  sd = ~1)))
  # One smooth per state, both starting where they always did
  expect_equal(as.vector(obs$lambda()), c(1, 1))
})

test_that("the likelihood reads penalty weights through L", {
  set.seed(1)
  n <- 200
  data <- data.frame(x = seq(-3, 3, length.out = n), z = NA)
  s <- 1
  for(i in 1:n) {
    data$z[i] <- stats::rnorm(1, c(0, 4)[s] + sin(data$x[i]), 0.7)
    s <- sample(c(s, 3 - s), size = 1, prob = c(0.9, 0.1))
  }
  obs <- Observation$new(data = data, dists = list(z = "norm"), n_states = 2,
                         par = list(z = list(mean = c(0, 4), sd = c(1, 1))),
                         formulas = list(z = list(mean = ~ s(x, k = 6, bs = "cs"),
                                                  sd = ~1)))
  hid <- MarkovChain$new(n_states = 2, data = data)
  hmm <- HMM$new(obs = obs, hid = hid)
  hmm$setup(silent = TRUE)

  args <- hmm$.__enclos_env__$private$tmb_args_
  # One smooth per state, so two penalties and two smoothing parameters
  n_pen <- nrow(args$data$L_obs)
  expect_equal(n_pen, 2L)
  expect_equal(args$data$L_obs, diag(n_pen))

  build <- function(L) {
    a <- args
    a$data$L_obs <- L
    a$data$include_smooths <- 1L
    do.call(TMB::MakeADFun, c(a, list(silent = TRUE)))
  }
  o1 <- build(diag(n_pen))
  o2 <- build(2 * diag(n_pen))

  i_lam <- which(names(o1$par) == "log_lambda_obs")
  expect_length(i_lam, n_pen)

  theta <- rep_len(c(0.37, -0.2), n_pen)
  # log(lambda) = L * theta, so doubling L at theta must equal L = I at 2 theta
  p1 <- o1$par; p1[i_lam] <- 2 * theta
  p2 <- o2$par; p2[i_lam] <- theta
  expect_equal(o2$fn(p2), o1$fn(p1))

  # and that this is not vacuous
  p3 <- o1$par; p3[i_lam] <- theta
  expect_false(isTRUE(all.equal(o1$fn(p3), o1$fn(p1))))
})
