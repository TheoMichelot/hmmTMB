
library(hmmTMB)

# Simulate data -----------------------------------------------------------

# number of time steps
n <- 1000

# create an empty dataset
dat <- data.frame(ID = rep(0, n), y = rep(0, n))

# create true model 
true_mod <- HMM$new(file = "true_mod.hmm")

# simulate from true model
dat <- true_mod$simulate(n)

# Fit model ------------------------------------------------------------

# create model to fit 
mod <- HMM$new(file = "mod.hmm")

mod$suggest_initial()

# fit model
mod$fit()

# look at numerical estimates
mod$par()
