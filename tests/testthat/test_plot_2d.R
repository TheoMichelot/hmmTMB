context("Plotting over two covariates")

# A correlated random walk over the unit square, as a tracked animal moves,
# with the transition out of state 2 driven by a spatial surface
sim_2d <- function(n = 800, seed = 5, band = FALSE) {
  set.seed(seed)
  x <- y <- numeric(n)
  x[1] <- y[1] <- 0.3
  for(t in 2:n) {
    x[t] <- min(max(x[t-1] + stats::rnorm(1, 0, 0.08), 0), 1)
    y[t] <- min(max(y[t-1] + stats::rnorm(1, 0, 0.08), 0), 1)
    # Optionally confine the walk to a diagonal band, so that much of the
    # bounding rectangle is extrapolation
    if(band && abs(x[t] - y[t]) > 0.25) { x[t] <- x[t-1]; y[t] <- y[t-1] }
  }
  u <- 1.5 * sin(4 * x) * cos(4 * y)
  s <- numeric(n); s[1] <- 1
  for(t in 2:n) {
    p <- if(s[t-1] == 1) stats::plogis(-1.2) else stats::plogis(-0.5 + u[t])
    s[t] <- if(stats::runif(1) < p) 3 - s[t-1] else s[t-1]
  }
  data.frame(ID = 1, x = x, y = y, z = stats::rnorm(n, c(0, 4)[s], 1))
}

# 'mesh' is an argument so that it is in scope where the formula strings are
# turned into formulas; a smooth's xt object has to be findable from there
build <- function(data, term, mesh = NULL) {
  form <- matrix(c(".", "~1", paste0("~ ", term), "."), 2, byrow = TRUE)
  hid <- MarkovChain$new(data = data, n_states = 2, formula = form,
                         initial_state = 1)
  obs <- Observation$new(data = data, n_states = 2, dists = list(z = "norm"),
                         par = list(z = list(mean = c(0, 4), sd = c(1, 1))))
  HMM$new(hid = hid, obs = obs)
}

# cov_grid_2d --------------------------------------------------------------

test_that("cov_grid_2d builds a regular lattice and pins the other covariates", {
  set.seed(1)
  data <- data.frame(x = stats::runif(80, 0, 4), y = stats::runif(80, -2, 2),
                     w = stats::rnorm(80),
                     g = factor(rep(c("a", "b"), 40)))
  grid <- cov_grid_2d("x", "y", data = data, formulas = list(a = ~ s(x, y) + w + g),
                      n_grid = 6)

  expect_equal(nrow(grid), 36L)
  expect_equal(sort(unique(grid$x)), seq(min(data$x), max(data$x), length = 6))
  expect_equal(sort(unique(grid$y)), seq(min(data$y), max(data$y), length = 6))
  # 'x' varies fastest, which is what geom_raster needs
  expect_equal(grid$x[1:6], seq(min(data$x), max(data$x), length = 6))
  expect_equal(length(unique(grid$y[1:6])), 1L)
  # Everything else is held at one value
  expect_equal(length(unique(grid$w)), 1L)
  expect_equal(unique(grid$w), mean(data$w))
  expect_equal(length(unique(grid$g)), 1L)
  # and can be set by hand
  grid2 <- cov_grid_2d("x", "y", data = data, covs = list(w = 3),
                       formulas = list(a = ~ s(x, y) + w), n_grid = 4)
  expect_equal(unique(grid2$w), 3)
})

test_that("cov_grid_2d rejects axes it cannot make a surface from", {
  set.seed(1)
  data <- data.frame(x = stats::runif(40), y = stats::runif(40),
                     g = factor(rep(c("a", "b"), 20)))
  expect_error(cov_grid_2d("x", "g", data = data, formulas = list(a = ~ s(x, y))),
               "not numeric")
  expect_error(cov_grid_2d("x", "nope", data = data, formulas = list(a = ~ s(x, y))),
               "Not a covariate")
})

# plot_2d ------------------------------------------------------------------

test_that("plot_2d draws a surface for any bivariate smoother", {
  skip_if_not_installed("fmesher")
  data <- sim_2d(n = 400)
  mesh <- fmesher::fm_mesh_2d(loc = cbind(data$x, data$y),
                              max.edge = c(0.2, 0.5), cutoff = 0.1,
                              offset = c(0.1, 0.3))
  terms <- c("s(x, y, k = 20)",                                   # thin plate
             "s(x, y, bs = 'gp', k = 20)",                        # Gaussian process
             "s(x, y, bs = 'spde', xt = list(mesh = mesh))")      # SPDE field

  for(term in terms) {
    hmm <- build(data, term, mesh = mesh)
    # Give the smooth some non-zero weights. An unfitted model has all of them
    # at zero and so a flat surface; this makes the surface depend on both
    # covariates through that basis, which is what is being tested here.
    set.seed(9)
    hmm$hid()$update_coeff_re(stats::rnorm(nrow(hmm$hid()$coeff_re()), 0, 0.5))

    # n_post = 0 needs no fit, so this tests the plotting itself
    p <- hmm$plot_2d("tpm", var = "x", var2 = "y", i = 2, j = 1,
                     n_grid = 12, too_far = 0)
    expect_s3_class(p, "ggplot")
    expect_equal(nrow(p$data), 144L, info = term)
    expect_true(all(c("var", "var2", "val") %in% names(p$data)), info = term)
    expect_true(all(is.finite(p$data$val)), info = term)
    # A probability, and not a constant one
    expect_true(all(p$data$val >= 0 & p$data$val <= 1), info = term)
    expect_gt(stats::sd(p$data$val), 0)
  }
})

test_that("plot_2d handles every model component", {
  data <- sim_2d(n = 400)
  hmm <- build(data, "s(x, y, k = 20)")

  expect_s3_class(hmm$plot_2d("tpm", "x", "y", n_grid = 8), "ggplot")
  expect_s3_class(hmm$plot_2d("delta", "x", "y", n_grid = 8), "ggplot")
  expect_s3_class(hmm$plot_2d("obspar", "x", "y", n_grid = 8), "ggplot")

  # tpm has one panel per transition, and i/j select among them
  all_tpm <- hmm$plot_2d("tpm", "x", "y", n_grid = 8)
  expect_equal(nrow(all_tpm$data), 4 * 64)
  one <- hmm$plot_2d("tpm", "x", "y", i = 2, j = 1, n_grid = 8)
  expect_equal(nrow(one$data), 64L)
  expect_equal(as.character(unique(one$data$from)), "State 2")
  expect_equal(as.character(unique(one$data$to)), "State 1")
})

test_that("too_far blanks the cells that are extrapolation", {
  data <- sim_2d(n = 600, band = TRUE)
  hmm <- build(data, "s(x, y, k = 20)")

  full <- hmm$plot_2d("tpm", "x", "y", i = 2, j = 1, n_grid = 20, too_far = 0)
  cut <- hmm$plot_2d("tpm", "x", "y", i = 2, j = 1, n_grid = 20, too_far = 0.1)

  expect_false(any(is.na(full$data$val)))
  expect_true(any(is.na(cut$data$val)))
  # The data sit on the diagonal, so it is the far corners that go
  blank <- cut$data[is.na(cut$data$val), ]
  expect_true(all(abs(blank$var - blank$var2) > 0.2))
  # and the cells that remain are unchanged
  keep <- !is.na(cut$data$val)
  expect_equal(cut$data$val[keep], full$data$val[keep])
})

test_that("show = 'ci' gives the width of the confidence interval", {
  skip_on_cran()
  data <- sim_2d(n = 1000)
  hmm <- build(data, "s(x, y, k = 12)")
  hmm$fit(silent = TRUE)

  p <- hmm$plot_2d("tpm", "x", "y", i = 2, j = 1, n_grid = 10,
                   n_post = 50, show = "ci", too_far = 0)
  expect_true(all(p$data$val >= 0))
  expect_equal(p$data$val, p$data$ucl - p$data$lcl)
  # and it differs from the estimate itself
  q <- hmm$plot_2d("tpm", "x", "y", i = 2, j = 1, n_grid = 10, too_far = 0)
  expect_false(isTRUE(all.equal(p$data$val, q$data$val)))
})

test_that("plot_2d rejects what it cannot plot", {
  data <- sim_2d(n = 300)
  hmm <- build(data, "s(x, y, k = 20)")

  expect_error(hmm$plot_2d("tpm", var = "x", var2 = "x"), "two different")
  expect_error(hmm$plot_2d("tpm", var = "x", var2 = "y", show = "ci"),
               "needs n_post")
  expect_error(hmm$plot_2d("tpm", var = "x", var2 = "nope"), "Not a covariate")
})

test_that("a one-dimensional field still plots through HMM$plot", {
  skip_if_not_installed("fmesher")
  skip_on_cran()
  set.seed(2)
  n <- 600
  x <- sort(stats::runif(n, 0, 10))
  u <- 1.5 * sin(x)
  s <- numeric(n); s[1] <- 1
  for(t in 2:n) {
    p <- if(s[t-1] == 1) stats::plogis(-1.2) else stats::plogis(-0.5 + u[t])
    s[t] <- if(stats::runif(1) < p) 3 - s[t-1] else s[t-1]
  }
  data <- data.frame(ID = 1, x = x, z = stats::rnorm(n, c(0, 4)[s], 1))
  hmm <- build(data, "s(x, bs = 'spde')")
  hmm$fit(silent = TRUE)

  # A field over one covariate needs nothing new: it is an mgcv smooth, and
  # HMM$plot() draws it like any other
  p <- hmm$plot("tpm", var = "x", i = 2, j = 1, n_post = 50)
  expect_s3_class(p, "ggplot")
  expect_true(all(c("var", "prob", "lcl", "ucl") %in% names(p$data)))
  expect_true(all(p$data$lcl <= p$data$prob & p$data$prob <= p$data$ucl))
})
