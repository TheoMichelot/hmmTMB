# Simulation to test performance on simple example 

library(hmmTMB)

# number of simulations
nsims <- 100
# list to save fitted models in 
mods <- vector(mode = "list", length = nsims)

## setup simulation
set.seed(295291)
# number of time steps
n <- 1000
# create an empty dataset
dat <- data.frame(ID = rep(0, n), count = rep(0, n))
# create true model 
true_mod <- HMM$new(file = "true_pois_mod.hmm")

for (i in 1:nsims) {
  cat(i, " / ", nsims, "\r")
  # simulate data
  dat <- true_mod$simulate(n, silent = TRUE)
  # fit model 
  mod <- HMM$new(file = "pois_mod.hmm")
  mod$fit(silent = TRUE)
  # save model
  mods[[i]] <- mod$clone()
}

# check bias on the estimates
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

# check confidence interval coverage 
preds <- lapply(mods, FUN = function(x) {x$predict("obspar", level = 0.95)})
lcls <- sapply(preds, FUN = function(x) {x$lcl[,,1]})
ucls <- sapply(preds, FUN = function(x) {x$ucl[,,1]})
cov <- rep(0, nrow(lcls))
for (i in 1:nrow(lcls)) {
  cov[i] <- sum((lcls[i,] < true_obs_ests[i]) & (ucls[i,] > true_obs_ests[i]))
}
# confidence interval coverage for observation parameters 
cov

tpmpreds <- lapply(mods, FUN = function(x) {x$predict("tpm", level = 0.95)})
tpmlcls <- sapply(tpmpreds, FUN = function(x) {diag(x$lcl[,,1])})
tpmucls <- sapply(tpmpreds, FUN = function(x) {diag(x$ucl[,,1])})
tpmcov <- rep(0, nrow(tpmlcls))
for (i in 1:nrow(tpmlcls)) {
  tpmcov[i] <- sum((tpmlcls[i,] > true_tpm_ests[i]) & (tpmucls[i,] < true_tpm_ests[i]))
}
# confidence interval coverage for observation parameters 
tpmcov


