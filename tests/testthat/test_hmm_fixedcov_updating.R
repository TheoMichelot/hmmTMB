
context("Test fixed covariate HMM with updating")
set.seed(58320)
# number of time steps
n <- 1000
# create an empty dataset
dat <- data.frame(ID = rep(0, n), count = rep(0, n))
# add covariates
dat$x <- rnorm(n)
dat$z <- rnorm(n)
dat$f <- factor(sample(1:3, size = n, replace = TRUE))
# create true model 
true_hid <- MarkovChain$new(data = dat, n_states = 2, 
                            tpm = matrix(c(0.8, 0.1, 0.2, 0.9), 2, 2),
                            formula = matrix(c(".", "~1", "~x", "."), 2, 2))
true_obs <- Observation$new(data = dat, n_states = 2, 
                            dists = list(count = "pois"), 
                            formulas = list(count = list(rate = ~ x + I(z^2) + f)),
                            par = list(count = list(rate = c(5, 20))))
true_mod <- HMM$new(obs = true_obs, hid = true_hid)
# set covariate effect on observation distribution 
par <- true_mod$coeff_fe()$obs
par <- c(par[1], -0.5, 0.08, 0.1, -0.8, 
         par[2], 0.8, 0.03, 0.1, 0.3)
true_mod$obs()$update_coeff_fe(par)
# set covariate effect on tpm 
tpmpar <- true_mod$coeff_fe()$hid
tpmpar <- c(tpmpar[1], 0.8, tpmpar[3])
true_mod$hid()$update_coeff_fe(tpmpar)
# simulate from true model
dat <- true_mod$simulate(n, data = dat, silent = TRUE)

# create model to fit 
hid <- MarkovChain$new(data = dat, n_states = 2, 
                       tpm = matrix(c(0.95, 0.05, 0.05, 0.95), 2, 2))
obs <- Observation$new(data = dat, n_states = 2, 
                       dists = list(count = "pois"), 
                       par = list(count = list(rate = c(5, 10))))
mod <- HMM$new(obs = obs, hid = hid)
# fit model
mod$fit(silent = TRUE)

# update observation model 
mod_x <- update(mod, "obs", "count", "rate", ~.+ x, silent = TRUE)
mod_z <- update(mod_x, "obs", "count", "rate", ~.+z, silent = TRUE)
mod_f <- update(mod_z, "obs", "count", "rate", ~.+f, silent = TRUE)
# update transition model 
mod_tpm0 <- update(mod_f, "hid", 1, 2, ~.+x, fit = FALSE, silent = TRUE)
mod_tpm <- update(mod_tpm0, "hid", 2, 1, ~.+x, silent = TRUE)

test_that("Updating works", {
  expect_equal(mod_x$obs()$formulas()$count$rate$state1, ~ 1 + x)
  expect_equal(mod_z$obs()$formulas()$count$rate$state1, ~ 1 + x + z)
  expect_equal(mod_f$obs()$formulas()$count$rate$state1, ~ 1 + x + z + f)
  expect_equal(mod_tpm0$hid()$formulas()[[1]], ~ x) 
  expect_equal(mod_tpm0$hid()$formulas()[[2]], ~ 1) 
  expect_equal(mod_tpm$hid()$formulas()[[1]], ~ x) 
  expect_equal(mod_tpm$hid()$formulas()[[2]], ~ x) 
})


test_that("AIC values are as expected", {
  expect_equal(as.numeric(AIC(mod, mod_x, mod_z, mod_f, mod_tpm)[,2]), 
               c(4316.4,3896.6,3896.8,3745.3,3728.4), tolerance = 0.1)
})
