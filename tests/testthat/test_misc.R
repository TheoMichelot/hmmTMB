
# Tests to add:
# update.HMM updates the correct formula in the correct way
# check format of output of make_formulas

context("Miscellaneous")

# Create dummy data
data <- data.frame(z = rnorm(100), x1 = rnorm(100), x2 = rnorm(100))
# MarkovChain object
hid <- MarkovChain$new(n_states = 2, structure = ~x1, data = data)
# Observation object
dists <- list(z = "norm")
par <- list(z = list(mean = c(0, 1), sd = c(0.5, 2)))
formulas <- list(z = list(mean = ~x1, sd = ~x2))
obs <- Observation$new(data = data, dists = dists, n_states = 2, 
                       par = par, formulas = formulas)
# HMM object
hmm <- HMM$new(obs = obs, hidden = hid)

test_that("update.HMM modifies formulas correctly", {
    # Update in hidden state model
    hmm2 <- update(hmm, type = "hidden", i = 2, j = 1, 
                   change = ~ . + I(x2^2), fit = FALSE)
    
    expect_equal(hmm2$hidden()$structure()[1,2], "~x1")
    expect_equal(hmm2$hidden()$structure()[2,1], "~x1 + I(x2^2)")
    expect_equal(hmm2$hidden()$formulas()[[1]], ~x1)
    expect_equal(hmm2$hidden()$formulas()[[2]], ~x1 + I(x2^2))
    
    # Update in observation model
    hmm3 <- update(hmm, type = "obs", i = "z", j = "mean", 
                   change = ~ . - x1 + state1(x2), fit = FALSE)
    hmm3 <- update(hmm3, type = "obs", i = "z", j = "sd", 
                   change = ~ . + s(x1, k = 5, bs = 'cs'), fit = FALSE)
    
    expect_equal(hmm3$obs()$formulas()$z$mean$state1, ~1 + x2)
    expect_equal(hmm3$obs()$formulas()$z$mean$state2, ~1)
    expect_equal(hmm3$obs()$formulas()$z$sd$state1, 
                 ~1 + x2 + s(x1, k = 5, bs = "cs"))
    expect_equal(hmm3$obs()$formulas()$z$sd$state2, 
                 ~1 + x2 + s(x1, k = 5, bs = "cs"))
})
