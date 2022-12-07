
context("HMM class")

# Create data with NAs in covariates
data <- data.frame(z = 1:100, x1 = rnorm(100), x2 = rnorm(100))
data$x1[1:10] <- NA
data$x2[91:100] <- NA

# Initial parameters
par <- list(z = list(mean = c(0, 0), sd = c(1, 1)))

test_that("NA covariates are filled in", {
    # Create model
    f <- ~x1 + s(x2, k = 5, bs = "cs")
    formulas <- list(z = list(mean = f, sd = ~1))
    obs <- Observation$new(data = data, dists = list(z = "norm"), n_states = 2, 
                           par = par, formulas = formulas)
    hid <- MarkovChain$new(n_states = 2, data = data)
    # NAs should be filled in at this stage
    hmm <- HMM$new(obs = obs, hid = hid)
    
    # Check that NAs were filled in with the right value
    expect_true(!any(is.na(obs$data()$x1)))
    expect_true(all(obs$data()$x1[1:10] == obs$data()$x1[11]))
    expect_true(!any(is.na(obs$data()$x1)))
    expect_true(all(obs$data()$x2[91:100] == obs$data()$x2[90]))
})

test_that("edf works in case with no smooth", {
    # Model with no covariates
    obs <- Observation$new(data = data, dists = list(z = "norm"), n_states = 2, par = par)
    hid <- MarkovChain$new(n_states = 2, data = data)
    hmm <- HMM$new(obs = obs, hid = hid)
    expect_equal(hmm$edf(), length(hid$coeff_fe()) + length(obs$coeff_fe()) + 
                   length(hmm$coeff_list()$log_delta0))
    
    # Model with fixed effects only
    f <- ~ x1*x2 + I(x2^2)
    obs <- Observation$new(data = data, dists = list(z = "norm"), n_states = 2, 
                           par = par, formulas = list(z = list(mean = ~1, sd = f)))
    hid <- MarkovChain$new(n_states = 2, data = data, formula = f)
    hmm <- HMM$new(obs = obs, hid = hid)
    expect_equal(hmm$edf(), length(hid$coeff_fe()) + length(obs$coeff_fe()) +
                   length(hmm$coeff_list()$log_delta0))
})

# Create dummy data
data <- data.frame(z = rnorm(100), x1 = rnorm(100), x2 = rnorm(100))
# MarkovChain object
hid <- MarkovChain$new(n_states = 2, formula = ~x1, data = data)
# Observation object
dists <- list(z = "norm")
par <- list(z = list(mean = c(0, 1), sd = c(0.5, 2)))
formulas <- list(z = list(mean = ~x1, sd = ~x2))
obs <- Observation$new(data = data, dists = dists, n_states = 2, 
                       par = par, formulas = formulas)
# HMM object
hmm <- HMM$new(obs = obs, hid = hid)

test_that("update.HMM modifies formulas correctly", {
    # Update in hidden state model
    hmm2 <- update(hmm, type = "hid", i = 2, j = 1, 
                   change = ~ . + I(x2^2), fit = FALSE)
    
    expect_equal(hmm2$hid()$formula()[1,2], "~x1")
    expect_equal(hmm2$hid()$formula()[2,1], "~x1 + I(x2^2)")
    expect_equal(hmm2$hid()$formulas()[[1]], ~x1)
    expect_equal(hmm2$hid()$formulas()[[2]], ~x1 + I(x2^2))
    
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
