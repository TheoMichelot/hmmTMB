
# Tests to add:
# check dimensions of X_re and X_fe and S
# check that elements of terms are compatible (e.g. same length...)
# check update functions
# check dimensions of tpm
# check delta for one special case?
# check that delta throws error when no solution
# check dimension of linpred

context("MarkovChain class")

test_that("Parameter transformation functions cancel out", {
    data <- data.frame(ID = rep(1, 100))
    mc <- MarkovChain$new(n_states = 2, data = data)
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
    data <- data.frame(ID = rep(1, 100))
    mc <- MarkovChain$new(n_states = 2, data = data)
    n <- 50
    sim <- mc$simulate(n = n)
    # Check length of output
    expect_equal(length(sim), n)
    # Check that output is in state space
    expect_equal(sort(unique(sim)), c(1, 2))
})
