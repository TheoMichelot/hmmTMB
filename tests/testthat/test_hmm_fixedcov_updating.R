
context("Test fixed covariate HMM with updating")
set.seed(58320)
# number of time steps
n <- 1000
# create an empty dataset
dat <<- data.frame(ID = rep(0, n), count = rep(0, n))
# add covariates
dat$x <<- rnorm(n)
dat$z <<- rnorm(n)
dat$f <<- factor(sample(1:3, size = n, replace = TRUE))
# create true model 
true_mod <<- HMM$new(file = "../../inst/examples/fixed_covs/updating_models/true_fixedcovsmod.hmm")
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
dat <<- true_mod$simulate(n, data = dat, silent = TRUE)
# create model to fit 
mod <<- HMM$new(file = "../../inst/examples/fixed_covs/updating_models/fixedcovs_mod.hmm")
# fit model
mod$fit(silent = TRUE)
# update observation model 
mod_x <<- update(mod, "obs", "count", "lambda", ~.+ x, silent = TRUE)
mod_z <<- update(mod_x, "obs", "count", "lambda", ~.+z, silent = TRUE)
mod_f <<- update(mod_z, "obs", "count", "lambda", ~.+f, silent = TRUE)
# update transition model 
mod_tpm0 <<- update(mod_f, "hid", 1, 2, ~.+x, fit = FALSE, silent = TRUE)
mod_tpm <<- update(mod_tpm0, "hid", 2, 1, ~.+x, silent = TRUE)

test_that("Updating works", {
  expect_equal(mod_x$obs()$formulas()$count$lambda$state1, ~ 1 + x)
  expect_equal(mod_z$obs()$formulas()$count$lambda$state1, ~ 1 + x + z)
  expect_equal(mod_f$obs()$formulas()$count$lambda$state1, ~ 1 + x + z + f)
  expect_equal(mod_tpm0$hid()$formulas()[[1]], ~ x) 
  expect_equal(mod_tpm0$hid()$formulas()[[2]], ~ 1) 
  expect_equal(mod_tpm$hid()$formulas()[[1]], ~ x) 
  expect_equal(mod_tpm$hid()$formulas()[[2]], ~ x) 
})


test_that("AIC values are as expected", {
  expect_equal(as.numeric(AIC(mod, mod_x, mod_z, mod_f, mod_tpm)[,2]), 
               c(4316.4,3896.6,3896.8,3745.3,3728.4), tolerance = 0.1)
})



