
context("Observation class")

# Create Observation object for tests
n <- 100
dists <- list(z = "norm")
data <- data.frame(z = rnorm(n), x1 = rnorm(n), x2 = rnorm(n))
f <- ~ x1*x2 + s(x2, k = 5, bs = "cs")
formulas <- list(z = list(mean = f, sd = ~1))
par0 <- list(z = list(mean = c(-1, 1), sd = c(0.5, 2)))
obs <- Observation$new(data = data, dists = dists, n_states = 2, 
                formulas = formulas, par = par0)

test_that("Specification stored correctly", {
    expect_equal(obs$data()[-4], data)
    expect_equal(obs$dists(), list(z = dist_norm))
    expect_equal(obs$nstates(), 2)
})

test_that("Correct parameter values returned", {
    expect_equal(as.numeric(obs$par()), c(-1, 0.5, 1, 2))
    expect_equal(as.numeric(obs$par(t = 1:5)), rep(c(-1, 0.5, 1, 2), 5))
})

test_that("Parameter transformation functions cancel out", {
    # Transform to linear predictor scale
    wpar <- obs$n2w(par = par0)
    # Transform back to natural scale
    par2 <- unlist(obs$w2n(wpar)[[1]])
    expect_equal(unname(unlist(par0)), unname(par2))
})

test_that("Model matrices have the right format", {
    terms <- obs$terms()
    expect_equal(dim(terms$X_fe), c(4*n, 10))
    expect_equal(dim(terms$X_re), c(4*n, 8))
    expect_equal(dim(terms$S), c(8, 8))
})

test_that("Model terms are consistent", {
    terms <- obs$terms()
    expect_equal(sum(terms$ncol_fe), ncol(terms$X_fe))
    expect_equal(sum(terms$ncol_re), ncol(terms$X_re))
})

test_that("linpred has correct length", {
    expect_equal(length(obs$linpred()), 4*n)
})

test_that("par has correct dimensions", {
    expect_equal(dim(obs$par()), c(2, 2, 1))
    expect_equal(dim(obs$par(t = 5)), c(2, 2, 1))
    expect_equal(dim(obs$par(t = 6:10)), c(2, 2, 5))
    expect_equal(dim(obs$par(t = "all")), c(2, 2, n))
    expect_equal(dim(obs$par(t = "all", linpred = seq(-2, 2, length = 4))), c(2, 2, 1))
})

test_that("obs_var has correct dimensions", {
    expect_equal(dim(obs$obs_var()), c(n, 1))
})

test_that("obs_probs has correct dimensions", {
    expect_equal(dim(obs$obs_probs()), c(n, 2))
})

test_that("Update methods work", {
    # Test update of par
    new_par <- list(z = list(mean = c(2, 3), sd = c(0.1, 0.2)))
    obs$update_par(par = new_par)
    expect_equal(c(t(obs$par()[,,1])), unname(unlist(new_par)))
    expect_equal(obs$coeff_fe()[c(1, 5, 9, 10)], c(new_par$z$mean, log(new_par$z$sd)))

    # Re-initialise object
    obs <- Observation$new(data = data, dists = dists, n_states = 2, 
                           formulas = formulas, par = par0)
    
    # Test update of coeff_fe/re
    new_fe <- 1:10
    new_re <- seq(0, 1, length = 8)
    obs$update_coeff_fe(coeff_fe = new_fe)
    obs$update_coeff_re(coeff_re = new_re)
    expect_equal(c(obs$coeff_fe()), new_fe)
    expect_equal(c(obs$coeff_re()), new_re)
    
    # Test update of lambda
    new_lambda <- c(10, 5)
    obs$update_lambda(lambda = new_lambda)
    expect_equal(c(obs$lambda()), new_lambda)
    
    # Re-initialise object
    obs <- Observation$new(data = data, dists = dists, n_states = 2, 
                           formulas = formulas, par = par0)
})

