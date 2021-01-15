
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
    hmm <- HMM$new(obs = obs, hidden = hid)
    
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
    hmm <- HMM$new(obs = obs, hidden = hid)
    expect_equal(hmm$edf(), length(hid$coeff_fe()) + length(obs$coeff_fe()))
    
    # Model with fixed effects only
    f <- ~ x1*x2 + I(x2^2)
    obs <- Observation$new(data = data, dists = list(z = "norm"), n_states = 2, 
                           par = par, formulas = list(z = list(mean = ~1, sd = f)))
    hid <- MarkovChain$new(n_states = 2, data = data, structure = f)
    hmm <- HMM$new(obs = obs, hidden = hid)
    expect_equal(hmm$edf(), length(hid$coeff_fe()) + length(obs$coeff_fe()))
})
