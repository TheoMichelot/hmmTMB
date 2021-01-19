
context("MarkovChain class")

# Create MarkovChain object for tests
n <- 100
data <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
f <- ~ x1 + s(x2, k = 5, bs = "cs")
mc <- MarkovChain$new(n_states = 2, structure = f, data = data)

test_that("Parameter transformation functions cancel out", {
    tpm <- matrix(c(0.9, 0.1,
                    0.3, 0.7),
                  ncol = 2, byrow = TRUE)
    # Transform to linear predictor scale
    par <- mc$tpm2par(tpm)
    # Transform back to tpm
    tpm2 <- mc$par2tpm(par)
    expect_equal(tpm, tpm2)
})

test_that("Simulation output is as expected", {
    nsim <- 50
    sim <- mc$simulate(n = nsim, data = data)
    # Check length of output
    expect_equal(length(sim), nsim)
    # Check that output is in state space
    expect_equal(sort(unique(sim)), c(1, 2))
})

test_that("Model matrices have the right format", {
    terms <- mc$terms()
    expect_equal(dim(terms$X_fe), c(2*n, 4))
    expect_equal(dim(terms$X_re), c(2*n, 8))
    expect_equal(dim(terms$S), c(8, 8))
})

test_that("Model terms are consistent", {
    terms <- mc$terms()
    expect_equal(sum(terms$ncol_fe), ncol(terms$X_fe))
    expect_equal(sum(terms$ncol_re), ncol(terms$X_re))
})

test_that("linpred has correct length", {
    expect_equal(length(mc$linpred()), 2*n)
})

test_that("TPM has correct dimensions", {
    expect_equal(dim(mc$tpm()), c(2, 2, 1))
    expect_equal(dim(mc$tpm(t = 5)), c(2, 2, 1))
    expect_equal(dim(mc$tpm(t = 6:10)), c(2, 2, 5))
    expect_equal(dim(mc$tpm(t = "all")), c(2, 2, n))
    expect_equal(dim(mc$tpm(t = "all", linpred = seq(-2, 2, length = 4))), c(2, 2, 2))
})

test_that("Stationary distribution has correct dimensions", {
    expect_equal(dim(mc$delta()), c(1, 2))
    expect_equal(dim(mc$delta(t = 5)), c(1, 2))
    expect_equal(dim(mc$delta(t = 6:10)), c(5, 2))
    expect_equal(dim(mc$delta(t = "all")), c(n, 2))
    expect_equal(dim(mc$delta(t = "all", linpred = seq(-2, 2, length = 4))), c(2, 2))
})

test_that("Update methods work", {
    # Test update of delta and tpm
    new_delta <- c(0.3, 0.7)
    mc$update_delta(delta = new_delta)
    expect_equal(unname(mc$delta(t = NULL)), new_delta)
    new_tpm <- matrix(c(0.8, 0.3, 0.2, 0.7), ncol = 2)
    mc$update_tpm(tpm = new_tpm)
    expect_equal(mc$tpm()[,,1], new_tpm)
    
    # Re-initialise object
    mc <- MarkovChain$new(n_states = 2, structure = f, data = data)
    
    # Test update of coeff_fe/re
    new_fe <- 1:4
    new_re <- seq(-3, 3, length = 8)
    mc$update_coeff_fe(coeff_fe = new_fe)
    mc$update_coeff_re(coeff_re = new_re)
    expect_equal(c(mc$coeff_fe()), new_fe)
    expect_equal(c(mc$coeff_re()), new_re)
    
    # Test update of lambda
    new_lambda <- c(2, 5)
    mc$update_lambda(lambda = new_lambda)
    expect_equal(c(mc$lambda()), new_lambda)
    
    # Re-initialise object
    mc <- MarkovChain$new(n_states = 2, structure = f, data = data)
})

test_that("delta is computed correctly", {
    # Check stationary distribution for a few special cases of TPM
    mc$update_tpm(tpm = matrix(c(0.5, 0.5, 0.5, 0.5), ncol = 2))
    expect_equal(mc$delta()[1,], c(0.5, 0.5))
    mc$update_tpm(tpm = matrix(c(0.9, 0.1, 0.1, 0.9), ncol = 2))
    expect_equal(mc$delta()[1,], c(0.5, 0.5))
    mc$update_tpm(tpm = matrix(c(0.9, 0.2, 0.1, 0.8), ncol = 2))
    expect_equal(mc$delta()[1,], c(2/3, 1/3))
    mc$update_tpm(tpm = matrix(c(0.9, 0.3, 0.1, 0.7), ncol = 2))
    expect_equal(mc$delta()[1,], c(3/4, 1/4))
    mc$update_tpm(tpm = matrix(c(0.9, 0.4, 0.1, 0.6), ncol = 2))
    expect_equal(mc$delta()[1,], c(4/5, 1/5))
    mc$update_tpm(tpm = matrix(c(0.9, 0.5, 0.1, 0.5), ncol = 2))
    expect_equal(mc$delta()[1,], c(5/6, 1/6))
    mc$update_tpm(tpm = matrix(c(0.1, 0.1, 0.9, 0.9), ncol = 2))
    expect_equal(mc$delta()[1,], c(0.1, 0.9))
    
    # Check that error is throw for singular system
    mc$update_tpm(tpm = matrix(c(1, 0, 0, 1), ncol = 2))
    expect_error(mc$delta(), "singular system")
    
    # Re-initialise object
    mc <- MarkovChain$new(n_states = 2, structure = f, data = data)
})
