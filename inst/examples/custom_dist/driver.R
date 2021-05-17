
library(hmmTMB)

# Create custom distribution ----------------------------------------------

# I'm going to create a copy of the weibull distribution (which is already
# included in hmmTMB) and call it weibull2, I only have to do this once and then
# this custom distribution will be stored in hmmTMB unless I reinstall the package

weibull2 <- Dist$new(
  name = "weibull2",
  pdf = dweibull,
  rng = rweibull,
  link = list(shape = log, scale = log),
  invlink = list(shape = exp, scale = exp),
  npar = 2, 
  parnames = c("shape", "scale"), 
  cpp = "weibull.cpp"
)

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

# fit model
mod$fit()

# look at numerical estimates
mod$par()
