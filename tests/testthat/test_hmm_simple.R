context("Basic tests for simple HMM")

# number of time steps
n <- 100
# create an empty dataset
dat <- data.frame(ID = rep(0, n), count = rep(0, n))
# create true model 
obs <- Observation$new(data = dat, 
                       dists = list(count = "pois"), 
                       n_states = 2, 
                       par = list(count = list(lambda = c(5, 20))))
hid <- MarkovChain$new(n_states = 2, data = dat)
true_mod <- HMM$new(obs = obs, hid = hid)
# simulate from true model
set.seed(57320)
dat <- true_mod$simulate(n, silent = TRUE)
# create model to fit 
obs <- Observation$new(data = dat, 
                       dists = list(count = "pois"), 
                       n_states = 2, 
                       par = list(count = list(lambda = c(5, 10))))
hid <- MarkovChain$new(n_states = 2, data = dat)
mod <- HMM$new(obs = obs, hid = hid)
# suggest better starting parameters?
ini <- mod$suggest_initial()

test_that("Simulate HMM data", {
  expect_equal(names(dat), c("ID", "count"))
  expect_equal(nrow(dat), 100)
  expect_equal(is_whole_number(dat[,2]), TRUE)
})

test_that("Model is setup correctly", {
  
  # check par 
  expect_equal(as.numeric(mod$par()$obspar), c(5, 10))
  expect_equal(as.numeric(mod$par()$tpm), c(0.9, 0.1, 0.1, 0.9))
  
  # check coeffs on link scale
  expect_equal(as.numeric(mod$coeff_fe()$obs), log(c(5, 10)))
  expect_equal(as.numeric(mod$coeff_fe()$hid), qlogis(c(0.1, 0.1)))
})

test_that("Model fits correctly", {
  # fit model 
  mod$fit(silent = TRUE)
  
  # check estimates are reasonable
  expect_equal(as.numeric(mod$par()$obspar), c(5, 20), tolerance = 1)
  expect_equal(as.numeric(mod$par()$tpm), c(0.9, 0.1, 0.1, 0.9), tolerance = 0.2)
  
  # get uncertainty, e.g., for count means 
  var <- mod$predict("obspar", level = 0.95)
  expect_equal(all(as.numeric(var$lcl) < as.numeric(var$ucl)), TRUE)
  expect_equal(as.numeric(var$lcl), c(5, 20), tolerance = 1)
  expect_equal(as.numeric(var$ucl), c(5, 20), tolerance = 1)
  
  # or for transition probabilities
  var <- mod$predict("tpm", level = 0.95)
  expect_equal(all(var$lcl[!diag(2)] < var$ucl[!diag(2)]), TRUE)
  expect_equal(var$lcl[!diag(2)], c(0.1, 0.1), tolerance = 0.1)
  expect_equal(var$ucl[!diag(2)], c(0.3, 0.3), tolerance = 0.1)
  
  # or for stationary distribution 
  var <- mod$predict("delta", level = 0.95)
  expect_equal(var$lcl[2] < var$ucl[2], TRUE)
  expect_equal(var$lcl[1], 0.57, tolerance = 0.01)
  expect_equal(var$ucl[1], 0.56, tolerance = 0.01)
})

test_that("State inference is correct", {
  # decode states
  states <- mod$viterbi()

  expect_equal(length(states), nrow(dat))  
  expect_equal(as.numeric(table(states)), c(57, 43))
  
  # sample possible state sequences 
  set.seed(52195)
  sim_states <- mod$sample_states(nsamp = 100)
  expect_equal(dim(sim_states), c(100, 100))
  expect_equal(all(sim_states %in% c(1, 2)), TRUE)
  
  # state probabilities 
  state_probs <- mod$state_probs()
  sim_probs <- apply(sim_states, 1, FUN = function(x) {sum(x == 1) / 100})
  expect_equal(all(abs(state_probs[,1] - sim_probs) < 0.02), TRUE)  
  expect_equal(rowSums(state_probs), rep(1, 100))
})

test_that("Pseudo-residuals look correct", {
  # get pseudo-residuals, should be normally distributed 
  resids <- mod$pseudores()
  expect_gt(ks.test(resids, "pnorm")$p.value, 0.4)
})

test_that("Goodness-of-fit metrics are correct", {
  # simulated-biased testing
  set.seed(51010)
  gof_stat <- function(x) {
    quantile(x$count, prob = seq(0.2, 0.8, 0.2))
  }
  # do simulation 
  sims <- mod$gof(gof_stat, nsims = 100, silent = TRUE)
  expect_equal(as.numeric(sims$obs_stat), c(4, 6, 15, 22), tolerance = 0.1)
})


