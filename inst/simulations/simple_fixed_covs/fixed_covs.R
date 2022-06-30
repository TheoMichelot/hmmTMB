# Simple Poisson HMM example of using hmmTMB and maximum likelihood

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
dat$x <- rnorm(n)
# create true model 
true_mod <- HMM$new(file = "true_fixedcovsmod.hmm")
# set covariate effect on observation distribution
par <- true_mod$coeff_fe()$obs
par <- c(par[1], -0.05, par[3], 0.08)
true_mod$obs()$update_coeff_fe(par)
# set covariate effect on transition probability 
tpmpar <- true_mod$coeff_fe()$hid
tpmpar <- c(tpmpar[1], 0.1, tpmpar[3])
true_mod$hid()$update_coeff_fe(tpmpar)

for (i in 1:nsims) {
  cat(i, " / ", nsims, "\r")
  # simulate data
  dat <- true_mod$simulate(n, data = dat, silent = TRUE)
  # fit model 
  mod <- HMM$new(file = "fixedcovs_mod.hmm")
  mod$fit(silent = TRUE)
  # save model
  mods[[i]] <- mod$clone()
}

# check bias on the estimates
obs_ests <- sapply(mods, FUN = function(x) {x$coeff_fe()$obs[,1]})
e_obs_ests <- rowMeans(obs_ests)
true_obs_ests <- true_mod$coeff_fe()$obs[,1]
bias_obs_ests <- 100 * (e_obs_ests - true_obs_ests) / true_obs_ests
# relative percentage bias in observation parameter estimates
bias_obs_ests

tpm_ests <- sapply(mods, FUN = function(x) {x$coeff_fe()$hid[,1]})
e_tpm_ests <- rowMeans(tpm_ests)
true_tpm_ests <- true_mod$coeff_fe()$hid[,1]
bias_tpm_ests <- 100 * (e_tpm_ests - true_tpm_ests) / true_tpm_ests
# relative percentage bias in transition probability parameter estimates
bias_tpm_ests

# check confidence interval coverage 
true_obs_ests <- true_mod$par()$obs[,,1]
preds <- lapply(mods, FUN = function(x) {x$predict("obspar", level = 0.95)})
lcls <- sapply(preds, FUN = function(x) {x$lcl[,,1]})
ucls <- sapply(preds, FUN = function(x) {x$ucl[,,1]})
cov <- rep(0, nrow(lcls))
for (i in 1:nrow(lcls)) {
  cov[i] <- sum((lcls[i,] < true_obs_ests[i]) & (ucls[i,] > true_obs_ests[i]))
}
# confidence interval coverage for observation parameters 
cov

true_tpm_ests <- diag(true_mod$par()$tpm[,,1])
tpmpreds <- lapply(mods, FUN = function(x) {x$predict("tpm", level = 0.95)})
tpmlcls <- sapply(tpmpreds, FUN = function(x) {diag(x$lcl[,,1])})
tpmucls <- sapply(tpmpreds, FUN = function(x) {diag(x$ucl[,,1])})
tpmcov <- rep(0, nrow(tpmlcls))
for (i in 1:nrow(tpmlcls)) {
  tpmcov[i] <- sum((tpmlcls[i,] > true_tpm_ests[i]) & (tpmucls[i,] < true_tpm_ests[i]))
}
# confidence interval coverage for observation parameters 
tpmcov
