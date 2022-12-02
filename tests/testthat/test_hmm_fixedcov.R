
context("Tests for HMM with fixed covariates")

set.seed(58320)
# number of time steps
n <- 1000
# create an empty dataset
dat <- data.frame(ID = rep(0, n), count = rep(0, n))
# add covariates
dat$x <- rnorm(n)
# create true model 
true_hid <- MarkovChain$new(data = dat, n_states = 2, 
                            tpm = matrix(c(0.8, 0.1, 0.2, 0.9), 2, 2),
                            formula = matrix(c(".", "~1", "~x", "."), 2, 2))
true_obs <- Observation$new(data = dat, n_states = 2, 
                            dists = list(count = "pois"), 
                            formulas = list(count = list(rate = ~ x)),
                            par = list(count = list(rate = c(5, 20))))
true_mod <- HMM$new(obs = true_obs, hid = true_hid)
# set covariate effect on observation distribution
par <- true_mod$coeff_fe()$obs
par <- c(par[1], -0.05, par[3], 0.08)
true_mod$obs()$update_coeff_fe(par)
# set covariate effect on transition probability 
tpmpar <- true_mod$coeff_fe()$hid
tpmpar <- c(tpmpar[1], 0.1, tpmpar[3])
true_mod$hid()$update_coeff_fe(tpmpar)
# simulate from true model
dat <- true_mod$simulate(n, data = dat, silent = TRUE)

# create model to fit 
hid <- MarkovChain$new(data = dat, n_states = 2, 
                       tpm = matrix(c(0.95, 0.05, 0.05, 0.95), 2, 2),
                       formula = matrix(c(".", "~1", "~x", "."), 2, 2))
obs <- Observation$new(data = dat, n_states = 2, 
                       dists = list(count = "pois"),
                       formulas = list(count = list(rate = ~ x)),
                       par = list(count = list(rate = c(5, 10))))
mod <- HMM$new(obs = obs, hid = hid)
# fit model
mod$fit(silent = TRUE)

test_that("Formulas are understood", {
  expect_equal(ncol(mod$obs()$make_mat()$X_fe), 4)
  expect_equal(ncol(mod$obs()$make_mat()$X_fe), 4)
})

test_that("Parameters are reasonable", {
  expect_equal(as.numeric(mod$coeff_fe()$obs[,1]), c(log(5), -0.05, log(20), 0.08), tolerance = 0.2)
  expect_equal(as.numeric(mod$coeff_fe()$hid[,1]), c(qlogis(0.2), 0.1, qlogis(0.1)), tolerance = 0.2)
})

test_that("Predictions can be made over time", {
  pred <- mod$predict("obspar", t = c(2,3,12), level = 0.95)
  expect_equal(as.numeric(pred$mean), c(4.9, 20.1, 5.8, 17.4, 5.9, 17.0), tolerance = 0.1)
  expect_equal(as.numeric(pred$lcl), c(4.7, 20.5, 5.3, 16.73, 5.3, 16.3), tolerance = 0.1)
  expect_equal(as.numeric(pred$ucl), c(5.2,21.4,6.3,18.1,6.5,17.8), tolerance = 0.1)
  expect_equal(all(as.numeric(pred$lcl) < as.numeric(pred$ucl)), TRUE)
  
  pred2 <- mod$predict("tpm", t = c(5, 1, 12), level = 0.95)
  expect_equal(as.numeric(pred2$mean), c(0.81,0.14,0.19,0.86,0.8,0.14,0.2,0.86,0.82,0.14,0.18,0.86), tolerance = 0.01)
  expect_equal(as.numeric(pred2$lcl), c(0.84,0.11,0.16,0.89,0.84,0.11,0.16,0.89,0.89,0.11,0.11,0.89), tolerance = 0.01)
  expect_equal(as.numeric(pred2$ucl), c(0.76,0.17,0.24,0.83,0.75,0.17,0.25,0.83,0.72,0.17,0.28,0.83), tolerance = 0.01)
})

