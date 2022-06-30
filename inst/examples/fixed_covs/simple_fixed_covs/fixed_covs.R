# Simple Poisson HMM example of using hmmTMB and maximum likelihood

library(hmmTMB)

# Simulate data -----------------------------------------------------------

# number of time steps
n <- 1000

# create an empty dataset
dat <- data.frame(ID = rep(0, n), count = rep(0, n))

# add covariates
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

# simulate from true model
set.seed(58320)
dat <- true_mod$simulate(n, data = dat)

# plot data
plot(dat$count, pch = 20, xlab = "time", ylab = "count")

# plot against covariate
plot(dat$x, dat$count, pch = 20, xlab = "x", ylab = "count")

# Fit basic model ------------------------------------------------------------

# create model to fit 
mod <- HMM$new(file = "fixedcovs_mod.hmm")

# check model formulation
mod

# fit model
mod$fit()

# Parameter Inference -----------------------------------------------------

# look at numerical estimates
mod$par()

# look at covariate effects
mod$coeff_fe()

# get uncertainty, e.g., for count means across first 10 time points 
# as x changes, so does the mean and uncertainty 
mod$predict("obspar", t = 1:10, level = 0.95)

# plot estimated relationship between lambda and x 
mod$plot("obspar", "x")

# get uncertainty on tpms for first 3 time points
mod$predict("tpm", t = 1:3, level = 0.95)

# plot estimated relationship between tpm(1, 2) entry and x
mod$plot("tpm", "x", i = 1, j = 2)


