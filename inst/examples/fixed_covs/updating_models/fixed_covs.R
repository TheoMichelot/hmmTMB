# Simple Poisson HMM example of using hmmTMB and maximum likelihood

library(hmmTMB)

# Simulate data -----------------------------------------------------------

# number of time steps
n <- 1000

# create an empty dataset
dat <- data.frame(ID = rep(0, n), count = rep(0, n))

# add covariates
dat$x <- rnorm(n)
dat$z <- rnorm(n)
dat$f <- factor(sample(1:3, size = n, replace = TRUE))

# create true model 
true_mod <- HMM$new(file = "true_fixedcovsmod.hmm")

# set covariate effect on observation distribution 
par <- true_mod$coeff_fe()$obs
par <- c(par[1], -0.5, 0.08, 0.1, -0.8, 
         par[2], 0.8, 0.03, 0.1, 0.3)
true_mod$obs()$update_coeff_fe(par)

# set covarite effect on tpm 
tpmpar <- true_mod$coeff_fe()$hid
tpmpar <- c(tpmpar[1], 0.8, tpmpar[3])
true_mod$hid()$update_coeff_fe(tpmpar)

# simulate from true model
set.seed(58320)
dat <- true_mod$simulate(n, data = dat)

# plot data
plot(dat$count, pch = 20, xlab = "time", ylab = "count")

# plot against covariate
plot(dat$x, dat$count, pch = 20, xlab = "x", ylab = "count")
plot(dat$z, dat$count, pch = 20, xlab = "z", ylab = "count")
plot(dat$count ~ dat$f)

# Fit basic model ------------------------------------------------------------

# create model to fit 
mod <- HMM$new(file = "fixedcovs_mod.hmm")

# check model formulation
mod

# fit model
mod$fit()


# Update with covariates --------------------------------------------------

## update observation model 
mod_x <- update(mod, "obs", "count", "lambda", ~.+ x)
AIC(mod, mod_x)

mod_z <- update(mod_x, "obs", "count", "lambda", ~.+z)
AIC(mod_x, mod_z)

mod_f <- update(mod_z, "obs", "count", "lambda", ~.+f)
AIC(mod_x, mod_z, mod_f)

# can also use BIC 
BIC(mod_x, mod_z, mod_f)

## update transition probability model
mod_tpm <- update(mod_f, "hid", 1, 2, ~.+x, fit = FALSE)
mod_tpm <- update(mod_tpm, "hid", 2, 1, ~.+x)

AIC(mod_f, mod_tpm)

# Parameter Inference -----------------------------------------------------

# look at numerical estimates
mod_tpm$par()

# look at covariate effects
mod_tpm$coeff_fe()

# get uncertainty, e.g., for count means across first 10 time points 
# as x changes, so does the mean and uncertainty 
mod_tpm$predict("obspar", t = 1:10, level = 0.95)

# predict parameter at certain covariate values
new_data <- data.frame(x = c(-1, 0, 1), z = c(-1.5, 1.5, 3), f = c(1, 1, 3))
mod_tpm$predict("obspar", t = "all", newdata = new_data, level = 0.95)

# plot estimated continuous relationships
mod_tpm$plot("obspar", "x")
mod_tpm$plot("obspar", "z")

# plot relationship with certain fixed values for other covariates
mod_tpm$plot("obspar", "x", covs = data.frame(z = 0, f = 2))

# plot estimated tpm relationship
mod_tpm$plot("tpm", "x", i = 2, j = 2)
