## Driver script

library(hmmTMB)

# simulate data 
n <- 10000
empty <- data.frame(ID = 1, count = rep(0, n))
true_mod <- HMM$new("inst/devel/true_mod.hmm")
dat <- true_mod$simulate(n)

# fit model 
mod <- HMM$new("inst/devel/mod.hmm")
mod$fit()
mod$par()
