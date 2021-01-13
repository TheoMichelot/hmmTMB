
# Tests to add:
# check dimension of par
# check dimension of obs_var
# check dimensions of obs_probs

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

test_that("Parameter transformation functions cancel out", {
    # Transform to linear predictor scale
    wpar <- obs$n2w(par = par0)
    # Transform back to natural scale
    par2 <- unlist(obs$w2n(coeff_fe = wpar)[[1]])
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
