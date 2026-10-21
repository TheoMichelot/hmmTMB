context("Posterior sampling")

# rmvn_prec ---------------------------------------------------------------

test_that("rmvn_prec draws from N(mu, Q^-1)", {
  set.seed(1)
  # A tridiagonal precision, so that the sparse path is the one exercised
  q <- 6
  Q <- Matrix::bandSparse(q, k = c(-1, 0, 1),
                          diagonals = list(rep(-0.4, q - 1), rep(2, q),
                                           rep(-0.4, q - 1)))
  mu <- seq_len(q) / 2

  draws <- rmvn_prec(n = 2e5, mu = mu, prec_mat = Q)

  expect_equal(dim(draws), c(2e5, q))
  expect_equal(colMeans(draws), mu, tolerance = 0.02)
  expect_equal(unname(stats::cov(draws)), as.matrix(solve(Q)), tolerance = 0.05)
})

test_that("rmvn_prec handles a dense precision the same way", {
  set.seed(2)
  A <- matrix(stats::rnorm(16), 4, 4)
  Q <- crossprod(A) + diag(4)
  mu <- c(-1, 0, 1, 2)

  draws <- rmvn_prec(n = 2e5, mu = mu, prec_mat = Q)

  expect_equal(colMeans(draws), mu, tolerance = 0.03)
  expect_equal(unname(stats::cov(draws)), solve(Q), tolerance = 0.05)
})

test_that("rmvn_prec falls back when the precision is not positive definite", {
  set.seed(3)
  # Singular: the last coefficient is unidentified
  Q <- diag(c(2, 2, 0))
  mu <- c(0, 0, 0)

  expect_warning(draws <- rmvn_prec(n = 100, mu = mu, prec_mat = Q))
  expect_equal(dim(draws), c(100L, 3L))
  expect_true(all(is.finite(draws)))
})

# post_coeff --------------------------------------------------------------

# A genuine smooth effect on the mean, so that the smoothing parameter has an
# interior optimum and sdreport returns a usable covariance matrix
set.seed(4)
n <- 300
data <- data.frame(x = seq(-3, 3, length.out = n), z = NA)
s <- 1
for(i in 1:n) {
  data$z[i] <- stats::rnorm(1, c(0, 5)[s] + sin(data$x[i]), 0.6)
  s <- sample(c(s, 3 - s), size = 1, prob = c(0.9, 0.1))
}
par0 <- list(z = list(mean = c(0, 5), sd = c(1, 1)))

test_that("a model with random effects samples from the joint precision", {
  obs <- Observation$new(data = data, dists = list(z = "norm"), n_states = 2,
                         par = par0,
                         formulas = list(z = list(mean = ~ s(x, k = 5, bs = "cs"),
                                                  sd = ~1)))
  hid <- MarkovChain$new(n_states = 2, data = data)
  hmm <- HMM$new(obs = obs, hid = hid)
  hmm$fit(silent = TRUE)

  rep <- hmm$tmb_rep()
  expect_false(is.null(rep$jointPrecision))

  set.seed(10)
  post <- hmm$post_coeff(n_post = 50)
  expect_equal(nrow(post), 50L)
  expect_true(all(is.finite(post)))

  # The route this replaced: invert the joint precision, then sample from it.
  # Same target distribution, so the two agree on the covariance.
  npar <- length(c(rep$par.fixed, rep$par.random))
  V_old <- prec_to_cov(rep$jointPrecision)
  set.seed(11)
  new_draws <- rmvn_prec(n = 3e4, mu = rep(0, npar),
                         prec_mat = rep$jointPrecision)
  expect_equal(unname(stats::cov(new_draws)), unname(as.matrix(V_old)),
               tolerance = 0.1)
})

test_that("a model without random effects is left on the covariance path", {
  obs <- Observation$new(data = data, dists = list(z = "norm"), n_states = 2,
                         par = par0)
  hid <- MarkovChain$new(n_states = 2, data = data)
  hmm <- HMM$new(obs = obs, hid = hid)
  hmm$fit(silent = TRUE)

  expect_null(hmm$tmb_rep()$jointPrecision)

  post <- hmm$post_coeff(n_post = 40)
  expect_equal(nrow(post), 40L)
  expect_true(all(is.finite(post)))
})
