
context("Forecast parameter validation")
options(warn = -1)  # disable warnings

# -- 1. Setup a simple HMM without covariates -------------------------------
n <- 10
dat_simple <- data.frame(ID = rep(1, n), count = rep(NA, n))
obs_simple <- Observation$new(
  data    = dat_simple,
  dists   = list(count = "pois"),
  n_states = 2,
  par     = list(count = list(rate = c(5, 10)))
)
hid_simple <- MarkovChain$new(n_states = 2, data = dat_simple)
mod_simple <- HMM$new(obs = obs_simple, hid = hid_simple)

test_that("`hmm` must be provided", {
  expect_error(
    Forecast$new(hmm = NULL, n = 5),
    "`hmm` must be provided"
  )
})

test_that("Supply either `n` or `forecast_data`", {
  expect_error(
    Forecast$new(hmm = mod_simple),
    "Supply either `n` or `forecast_data`"
  )
})

test_that("`n` must be a positive integer", {
  expect_error(
    Forecast$new(hmm = mod_simple, n = 0),
    "`n` must be a positive integer"
  )
})

test_that("`forecast_data` must be a data.frame", {
  expect_error(
    Forecast$new(hmm = mod_simple, forecast_data = list(a = 1)),
    "`forecast_data` must be a data.frame"
  )
})

# -- 2. Setup an HMM with covariates ---------------------
mod_df <- data.frame(ID = rep(1, 5), count = rep(0, 5))
mod_df$x <- runif(5)
cov_hid <- MarkovChain$new(data = mod_df, n_states = 2, 
                            tpm = matrix(c(0.8, 0.1, 0.2, 0.9), 2, 2),
                            formula = matrix(c(".", "~1", "~x", "."), 2, 2))
cov_obs <- Observation$new(data = mod_df, n_states = 2, 
                            dists = list(count = "pois"), 
                            formulas = list(count = list(rate = ~ x)),
                            par = list(count = list(rate = c(5, 20))))
mod_cov <- HMM$new(obs = cov_obs, hid = cov_hid)

good_df <- data.frame(ID = rep(1, 5), count = rep(0, 5), x = runif(5))
bad_df <- data.frame(ID = rep(1, 5), count = rep(0, 5)) # missing 'x'

test_that("Provide `forecast_data` when covariates are in the model", {
  expect_error(
    Forecast$new(hmm = mod_cov, n = 5),
    "Provide `forecast_data` when covariates are in the model"
  )
})

test_that("`forecast_data` is missing required covariates", {
  expect_error(
    Forecast$new(hmm = mod_cov, forecast_data = bad_df),
    "`forecast_data` is missing covariates: x"
  )
})

test_that("`preset_eval_range` must be a *named* list", {
  expect_error(
    Forecast$new(
      hmm           = mod_cov,
      forecast_data = good_df,
      preset_eval_range = list(runif(5))
    ),
    "`preset_eval_range` must be a \\*named\\* list"
  )
})

test_that("`starting_state_distribution` must be provided", {
  expect_error(
    Forecast$new(
      hmm                          = mod_cov,
      forecast_data                = good_df,
      starting_state_distribution  = NULL
    ),
    "`starting_state_distribution` must be provided"
  )
})

test_that("`starting_state_distribution` must be 'last' or 'stationary'", {
  expect_error(
    Forecast$new(
      hmm                          = mod_cov,
      forecast_data                = good_df,
      starting_state_distribution  = "foo"
    ),
    "Character `starting_state_distribution` must be 'last' or 'stationary'"
  )
})

test_that("`starting_state_distribution` numeric vector must have length n_states()", {
  expect_error(
    Forecast$new(
      hmm                          = mod_cov,
      forecast_data                = good_df,
      starting_state_distribution  = c(0.5, 0.5, 0.5)
    ),
    "`starting_state_distribution` must have length nstates\\(\\)"
  )
})

# -- 3. A valid call succeeds ------------------------------------------------
test_that("Valid arguments create a forecast object", {
  fc <- Forecast$new(
    hmm           = mod_cov,
    forecast_data = good_df
  )
  expect_true(inherits(fc, "Forecast"))
  expect_equal(ncol(fc$forecast_dists()[[1]]), nrow(good_df))
  expect_equal(names(fc$forecast_dists()), c('count'))
})

# -- 4. Forecasting with unsupported distributions --------------------------
test_that("`tweedie` distribution is not supported for forecasting", {
  dat_tweedie <- data.frame(ID = rep(1, 5), y = rep(1, 5))
  obs_tweedie <- Observation$new(
    data    = dat_tweedie,
    dists   = list(y = "tweedie"),
    n_states = 2,
    par     = list(y = list(mean = c(1, 2), p = c(1.5, 1.5), phi = c(1, 1)))
  )
  tpm <- matrix(c(0.15, 0.9, 0.85, 0.1), 2, 2)
  starting_state <- c(0.9, 0.1)
  hid_tweedie <- MarkovChain$new(n_states = 2, data = dat_tweedie, tpm = tpm)
  hid_tweedie$update_delta0(starting_state)
  mod_tweedie <- HMM$new(obs = obs_tweedie, hid = hid_tweedie)
  expect_error(
    Forecast$new(hmm = mod_tweedie, n = 5, starting_state_distribution = starting_state %*% tpm),
    "Forecasting not currently supported for 'cat', 'dirc', or 'tweedie' distributions"
  )
})
