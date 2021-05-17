
library(hmmTMB)

# Simulate data -----------------------------------------------------------

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
par <- true_mod$coeff_re()$obs
par[1:n_ID] <- rnorm(n_ID, 0, re_sd[1])
par[(n_ID + 1):nrow(par)] <- rnorm(n_ID, 0, re_sd[2]) 
true_mod$obs()$update_coeff_re(par)

# simulate from true model
set.seed(75145)
dat <- true_mod$simulate(n * n_ID, data = dat)

plot(dat$count ~ dat$ID)

# Fit fixed effects only model --------------------------------------------

# create model to fit 
fxmod <- HMM$new(file = "fixed_effects_only_mod.hmm")

# check model formulation
fxmod

# fit model
fxmod$fit()


# Fit random effects model ------------------------------------------------

# create model to fit 
mod <- HMM$new(file = "individual_random_effects_mod.hmm")

# check model formulation
mod

# fit model
mod$fit()


# Parameter Inference -----------------------------------------------------

# look at numerical estimates
fxmod$par()
mod$par()
 
# look at random effect standard deviations
mod$vcomp()

# look at BLUP estimates of random effects
mod$coeff_re()

# compare to true random effects: 
plot(mod$coeff_re()$obs[,1], par[,1])
abline(a = 0, b = 1, col = "steelblue")

