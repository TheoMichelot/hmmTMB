library(hmmTMB)

# number of simulations
nsims <- 100
# list to save fitted models in 
mods <- vector(mode = "list", length = nsims)

## setup simulation
set.seed(295291)
# number of time steps
n <- 1000
# number of individual time series
n_ID <- 10
# create an empty dataset
dat <- data.frame(ID = factor(rep(1:n_ID, each = n)), count = rep(0, n * n_ID))
# create true model 
true_mod <- HMM$new(file = "individual_random_effects_truemod.hmm")
# set random effects  
re_sd <- c(0.1, 0.3)

pb <- txtProgressBar(min = 0, max = nsims, style = 3)
for (i in 1:nsims) {
  setTxtProgressBar(pb, i)
  # simulate random effects 
  par <- true_mod$coeff_re()$obs
  par[1:n_ID] <- rnorm(n_ID, 0, re_sd[1])
  par[(n_ID + 1):nrow(par)] <- rnorm(n_ID, 0, re_sd[2]) 
  true_mod$obs()$update_coeff_re(par)
  # simulate data
  dat <- true_mod$simulate(n, data = dat, silent = TRUE)
  # fit model 
  mod <- HMM$new(file = "individual_random_effects_mod.hmm")
  mod$fit(silent = TRUE)
  # save model
  mods[[i]] <- mod$clone()
}

# check bias on the estimates
true_mod <- HMM$new(file = "individual_random_effects_truemod.hmm")
obs_ests <- sapply(mods, FUN = function(x) {x$par()$obspar[,,1]})
e_obs_ests <- rowMeans(obs_ests)
true_obs_ests <- true_mod$par()$obspar[,,1]
bias_obs_ests <- 100 * (e_obs_ests - true_obs_ests) / true_obs_ests
# relative percentage bias in observation parameter estimates
bias_obs_ests

tpm_ests <- sapply(mods, FUN = function(x) {diag(x$par()$tpm[,,1])})
e_tpm_ests <- rowMeans(tpm_ests)
true_tpm_ests <- diag(true_mod$par()$tpm[,,1])
bias_tpm_ests <- 100 * (e_tpm_ests - true_tpm_ests) / true_tpm_ests
# relative percentage bias in transition probability parameter estimates
bias_tpm_ests
