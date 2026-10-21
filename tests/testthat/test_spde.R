context("SPDE smooths and Gaussian fields")

# The bookkeeping these rely on -- one determinant per penalty, and the L
# convention -- is tested in test_make_mat.R.

# The smoother itself ------------------------------------------------------

test_that("the SPDE smoother builds three penalties and two parameters", {
  skip_if_not_installed("fmesher")
  set.seed(1)
  data <- data.frame(x = stats::runif(300), y = stats::runif(300))
  mesh <- fmesher::fm_mesh_2d(loc = cbind(data$x, data$y),
                              max.edge = c(0.2, 0.5), cutoff = 0.1,
                              offset = c(0.1, 0.3))
  sm <- mgcv::smoothCon(mgcv::s(x, y, bs = "spde", xt = list(mesh = mesh)),
                        data = data, absorb.cons = FALSE)[[1]]

  expect_s3_class(sm, "spde.smooth")
  # The flag make_matrices() reads, rather than the class, is what says the
  # log-determinant belongs to the likelihood
  expect_true(isTRUE(sm$gmrf))
  expect_length(sm$S, 3)
  expect_equal(dim(sm$L), c(3L, 2L))
  expect_equal(sm$theta.names, c("sd", "range"))
  expect_equal(sm$null.space.dim, 0L)     # a proper field needs no constraint
  expect_equal(nrow(sm$C), 0L)
  expect_true(sm$no.rescale)
  expect_equal(ncol(sm$X), mesh$n)
  expect_true(inherits(sm$X, "Matrix"))   # the design matrix stays sparse

  # The (sd, range) parameterisation is exactly the (tau, kappa) one: at
  # theta = (log sigma, log rho), sum_i exp((L theta)_i) S_i must reproduce
  # tau^2 (kappa^4 C + 2 kappa^2 G1 + G2)
  sigma <- 0.7
  rho <- 0.4
  kappa <- 2 * sqrt(2) / rho
  tau <- 1 / (sigma * kappa * sqrt(4 * pi))
  lambda <- exp(sm$L %*% c(log(sigma), log(rho)))
  Q1 <- lambda[1] * sm$S[[1]] + lambda[2] * sm$S[[2]] + lambda[3] * sm$S[[3]]
  fem <- fmesher::fm_fem(mesh)
  Q2 <- tau^2 * (kappa^4 * fem$c0 + 2 * kappa^2 * fem$g1 + fem$g2)
  expect_lt(max(abs(as.matrix(Q1 - Q2))) / max(abs(as.matrix(Q2))), 1e-12)
})

test_that("the one-dimensional smoother uses the Matern constants for d = 1", {
  skip_if_not_installed("fmesher")
  set.seed(1)
  data <- data.frame(x = stats::runif(200, 2, 8))
  sm <- mgcv::smoothCon(mgcv::s(x, k = 12, bs = "spde"), data = data,
                        absorb.cons = FALSE)[[1]]

  # alpha = 2 gives nu = 2 - d/2, so nu = 3/2 here and nu = 1 in two
  # dimensions, and both the constants and the powers of rho differ. As in the
  # two-dimensional case, sum_i exp((L theta)_i) S_i has to reproduce
  # tau^2 (kappa^4 C + 2 kappa^2 G1 + G2) at theta = (log sigma, log rho).
  sigma <- 0.7
  rho <- 1.3
  kappa <- sqrt(8 * 1.5) / rho
  tau <- sqrt(gamma(1.5) / (gamma(2) * sqrt(4 * pi) * kappa^3 * sigma^2))
  lambda <- exp(sm$L %*% c(log(sigma), log(rho)))
  Q1 <- lambda[1] * sm$S[[1]] + lambda[2] * sm$S[[2]] + lambda[3] * sm$S[[3]]
  fem <- fmesher::fm_fem(sm$mesh)
  Q2 <- tau^2 * (kappa^4 * fem$c0 + 2 * kappa^2 * fem$g1 + fem$g2)
  expect_lt(max(abs(as.matrix(Q1 - Q2))) / max(abs(as.matrix(Q2))), 1e-12)

  # And so sd and range mean what they say: on a mesh fine enough for the
  # discretisation not to bite, the implied marginal standard deviation is
  # sigma and the correlation at lag rho is the Matern one at kappa * rho
  fine <- mgcv::smoothCon(mgcv::s(x, k = 400, bs = "spde"), data = data,
                          absorb.cons = FALSE)[[1]]
  lambda <- exp(fine$L %*% c(log(sigma), log(rho)))
  Q <- lambda[1] * fine$S[[1]] + lambda[2] * fine$S[[2]] + lambda[3] * fine$S[[3]]
  g <- seq(3, 7, length.out = 41)
  A <- as.matrix(fmesher::fm_basis(fine$mesh, g))
  V <- A %*% solve(as.matrix(Q)) %*% t(A)
  expect_equal(sqrt(V[21, 21]), sigma, tolerance = 0.01)
  r <- V[21, ] / sqrt(V[21, 21] * diag(V))
  expect_equal(stats::approx(g[21:41] - g[21], r[21:41], rho)$y,
               (1 + kappa * rho) * exp(-kappa * rho), tolerance = 0.02)
})

test_that("a one-dimensional SPDE smooth builds its own mesh", {
  skip_if_not_installed("fmesher")
  set.seed(1)
  data <- data.frame(x = stats::runif(200, 2, 8))

  # k is a basis dimension, as everywhere else in mgcv, but defaults to 15
  # rather than 10: it is the resolution of the field, not a cap on wiggliness
  # that a penalty pulls back from. A degree-2 mesh has one basis function
  # fewer than it has knots, so the two can only agree if the knots are
  # counted backwards from k.
  sm <- mgcv::smoothCon(mgcv::s(x, bs = "spde"), data = data,
                        absorb.cons = FALSE)[[1]]
  expect_s3_class(sm$mesh, "fm_mesh_1d")
  expect_equal(ncol(sm$X), 15L)
  expect_equal(sm$mesh$degree, 2)
  expect_equal(ncol(mgcv::smoothCon(mgcv::s(x, k = 25, bs = "spde"), data = data,
                                    absorb.cons = FALSE)[[1]]$X), 25L)

  # Evenly spaced knots over the range of the covariate, widened by a fifth on
  # each side so that the boundary condition does not inflate the field's
  # variance where the data are
  pad <- diff(range(data$x)) / 5
  expect_equal(range(sm$mesh$loc), range(data$x) + c(-pad, pad))
  expect_equal(diff(sm$mesh$loc), rep(diff(range(sm$mesh$loc)) / 15, 15))

  mats <- hmmTMB:::make_matrices(list(a = ~ s(x, k = 15, bs = "spde")), data = data)
  expect_equal(unname(diff(mats$ncol_re[, 1]) + 1), 15)

  q <- diff(mats$ncol_re[, 1]) + 1
  expect_equal(ncol(mats$ncol_re), 3)          # three penalties
  expect_equal(mats$gmrf, rep(1L, 3))
  expect_equal(dim(mats$L), c(3L, 2L))         # combined by two parameters
  expect_equal(mats$sp_gmrf, c(1L, 1L))
  expect_equal(mats$sp_names, c("a.s(x).sd", "a.s(x).range"))
  expect_equal(mats$theta_start[1], 0)         # unit marginal sd
  # All three penalties apply to the same coefficients, which is what groups
  # them into one smooth, and their blocks of S are stacked in order
  expect_equal(unname(mats$ncol_re[1, ]), rep(1, 3))
  expect_equal(unname(mats$ncol_re[2, ]), rep(q, 3))
  expect_equal(ncol(mats$S), 3 * q)
})

test_that("SPDE smooths are rejected where they cannot work", {
  skip_if_not_installed("fmesher")
  set.seed(1)
  data <- data.frame(x = stats::runif(50), y = stats::runif(50),
                     z = stats::runif(50))
  expect_error(
    mgcv::smoothCon(mgcv::s(x, y, bs = "spde"), data = data),
    "needs a mesh")
  expect_error(
    mgcv::smoothCon(mgcv::s(x, y, z, bs = "spde"), data = data),
    "one- or two-dimensional")
  expect_error(
    mgcv::smoothCon(mgcv::s(x, bs = "spde", xt = list(mesh = "not a mesh")),
                    data = data),
    "must be an fmesher mesh")
  expect_error(
    mgcv::smoothCon(mgcv::s(x, k = 2, bs = "spde"), data = data),
    "k of at least 3")
  # A mesh of the wrong dimension would otherwise be given the Matern
  # constants of the other one, silently
  mesh <- fmesher::fm_mesh_2d(loc = cbind(data$x, data$y),
                              max.edge = c(0.2, 0.5), cutoff = 0.1)
  expect_error(
    mgcv::smoothCon(mgcv::s(x, bs = "spde", xt = list(mesh = mesh)), data = data),
    "2-dimensional mesh was given for a smooth of 1 covariate")
  # And a constant covariate has no range to spread knots over
  expect_error(
    mgcv::smoothCon(mgcv::s(x, bs = "spde"),
                    data = data.frame(x = rep(1, 50))),
    "has to vary and be finite")
})

# End to end ---------------------------------------------------------------

## A 2-state HMM whose transition 2 -> 1 is driven by a smooth spatial field,
## observed through a Gaussian state-dependent distribution. Small and coarse
## on purpose, so that the whole thing fits in a few seconds.
sim_field <- function(n = 2000, seed = 3) {
  set.seed(seed)
  x <- y <- numeric(n)
  x[1] <- y[1] <- 0.5
  for (t in 2:n) {
    x[t] <- min(max(x[t - 1] + stats::rnorm(1, 0, 0.1), 0), 1)
    y[t] <- min(max(y[t - 1] + stats::rnorm(1, 0, 0.1), 0), 1)
  }
  u <- 1.5 * sin(4 * x) * cos(4 * y)
  s <- numeric(n)
  s[1] <- 1
  for (t in 2:n) {
    p <- if (s[t - 1] == 1) stats::plogis(-1.2) else stats::plogis(-0.5 + u[t])
    s[t] <- if (stats::runif(1) < p) 3 - s[t - 1] else s[t - 1]
  }
  data.frame(ID = 1, x = x, y = y, z = stats::rnorm(n, c(0, 4)[s], 1))
}

test_that("an SPDE field in the transition probabilities can be fitted", {
  skip_if_not_installed("fmesher")
  skip_on_cran()
  data <- sim_field()
  mesh <- fmesher::fm_mesh_2d(loc = cbind(data$x, data$y),
                              max.edge = c(0.15, 0.5), cutoff = 0.08,
                              offset = c(0.1, 0.3))
  form <- matrix(c(".", "~1",
                   "~ s(x, y, bs = 'spde', xt = list(mesh = mesh))", "."),
                 2, byrow = TRUE)
  hid <- MarkovChain$new(data = data, n_states = 2, formula = form,
                         initial_state = 1)
  obs <- Observation$new(data = data, n_states = 2, dists = list(z = "norm"),
                         par = list(z = list(mean = c(0, 4), sd = c(1, 1))))
  hmm <- HMM$new(hid = hid, obs = obs)

  # A field switches the banded forward algorithm on by itself
  expect_equal(hmm$bw(), 15L)

  # The two smoothing parameters are named, and start where the smoother says
  expect_equal(rownames(hid$lambda()),
               c("S2>S1.s(x,y).sd", "S2>S1.s(x,y).range"))
  expect_equal(unname(hid$lambda()[1, 1]), 1)

  hmm$fit(silent = TRUE)

  sp <- hmm$lambda()$hid
  expect_true(all(is.finite(sp)))
  expect_true(all(sp > 0))
  # Neither parameter has run off: a range approaching the size of the domain
  # would mean the field is unidentified
  expect_lt(sp[2, 1], 2)

  # 1/sqrt(lambda) is not a standard deviation for these two
  expect_true(all(is.na(hid$sd_re())))

  # The intercepts recover the transition probabilities that carry no field
  expect_equal(unname(hmm$hid()$coeff_fe()[1, 1]), -1.2, tolerance = 0.25)

  # And the field itself is recovered on a grid the model never saw
  grid <- expand.grid(x = seq(0.05, 0.95, length = 20),
                      y = seq(0.05, 0.95, length = 20))
  tpm <- hmm$predict("tpm", newdata = grid)
  fitted <- stats::qlogis(tpm[2, 1, ])
  truth <- -0.5 + 1.5 * sin(4 * grid$x) * cos(4 * grid$y)
  expect_gt(stats::cor(fitted, truth), 0.8)
})

test_that("a one-dimensional field can be fitted from the default mesh", {
  skip_if_not_installed("fmesher")
  skip_on_cran()
  set.seed(2)
  n <- 2000
  x <- sort(stats::runif(n, 0, 10))
  u <- 1.5 * sin(x)
  s <- numeric(n)
  s[1] <- 1
  for (t in 2:n) {
    p <- if (s[t - 1] == 1) stats::plogis(-1.2) else stats::plogis(-0.5 + u[t])
    s[t] <- if (stats::runif(1) < p) 3 - s[t - 1] else s[t - 1]
  }
  data <- data.frame(ID = 1, x = x, z = stats::rnorm(n, c(0, 4)[s], 1))

  # No mesh, no k: the whole of the user-facing side of a one-dimensional field
  form <- matrix(c(".", "~1", "~ s(x, bs = 'spde')", "."), 2, byrow = TRUE)
  hid <- MarkovChain$new(data = data, n_states = 2, formula = form,
                         initial_state = 1)
  obs <- Observation$new(data = data, n_states = 2, dists = list(z = "norm"),
                         par = list(z = list(mean = c(0, 4), sd = c(1, 1))))
  hmm <- HMM$new(hid = hid, obs = obs)
  hmm$fit(silent = TRUE)

  g <- data.frame(x = seq(0.2, 9.8, length.out = 200))
  tpm <- hmm$predict("tpm", newdata = g)
  expect_gt(stats::cor(stats::qlogis(tpm[2, 1, ]), -0.5 + 1.5 * sin(g$x)), 0.9)
  expect_equal(unname(hmm$hid()$coeff_fe()[1, 1]), -1.2, tolerance = 0.25)

  # The default mesh has to be fine enough that both parameters stay
  # interpretable. The failure mode of too coarse a mesh is not a flat fit but
  # a range that collapses towards zero with the standard deviation growing to
  # compensate, which k = 10 does on this data set and k = 15 does not.
  sp <- hmm$lambda()$hid
  spacing <- 1.4 * diff(range(data$x)) / 15     # the default mesh, widened
  expect_gt(sp[2, 1], spacing)
  expect_lt(sp[1, 1], 10)
})

test_that("a model without a field is left on the exact algorithm", {
  data <- sim_field(n = 300)
  hid <- MarkovChain$new(data = data, n_states = 2, formula = ~ s(x, k = 5, bs = "cs"),
                         initial_state = 1)
  obs <- Observation$new(data = data, n_states = 2, dists = list(z = "norm"),
                         par = list(z = list(mean = c(0, 4), sd = c(1, 1))))
  hmm <- HMM$new(hid = hid, obs = obs)
  expect_equal(hmm$bw(), 0L)
  expect_false(is.na(hid$sd_re()[1, 1]))
})
