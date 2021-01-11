# Simple Poisson HMM example of using hmmTMB, tmbstan, and stan 

library(hmmTMB)

# Simulate data -----------------------------------------------------------

# number of time steps
n <- 1000

# create an empty dataset
dat <- data.frame(ID = rep(0, n), count = rep(0, n))

# create true model 
true_mod <- HMM$new(file = "true_pois_mod.hmm")

# simulate from true model
set.seed(58320)
dat <- true_mod$simulate(n)

# plot data
plot(dat$count, pch = 20, xlab = "time", ylab = "count")

# Create model ------------------------------------------------------------

# create model to fit 
mod <- HMM$new(file = "pois_mod.hmm")

# check model formulation
mod


# Set priors --------------------------------------------------------------

# look at current parameters on the working scale (in particular, their names)
mod$coeff_fe()

# priors are constrained to be normal distributions for all parameters on the
# link scale
# so for lambda, the link is log, and we must specify the prior for log(lambda)
# I must specify the mean and standard deviation in the same format as already
# stored in mod$priors() 
priors <- list(coeff_fe_obs = matrix(c(log(25), 0.5, 
                                       log(5), 0.5), nr = 2, byrow = TRUE))

# for transition probability matrices, I won't set a prior and so by default
# no prior is used, equivalent to an improper flat prior on the link scale 

# set priors 
mod$set_priors(priors)

# check they have been set properly
mod$priors()

# Fit model ---------------------------------------------------------------

# using stan MCMC and default stan arguments / random initial values
# but setting chains to 1 and iters to 100 as this is just an example 
mod$mcmc(chains = 1, iter = 100)

# Inference ---------------------------------------------------------------

# the fitted stan object can be accessed: 
mod$stan()

# and so you can use all the functions available to stan objects
# if you look the rstan package or equivalent
library(rstan)
stan_trace(mod$stan())
stan_dens(mod$stan())

# Alternatively, you can use the mod object 
# which by default has set all parameters to their posterior means 
mod$par()

# You can also set the parameters inside the model to any iteration
# (where iterations are combined across chains) 
mod$update_par(iter = 50)
mod$par()

# The $predict method is only applicable to maximum likelihood methods 
# where the model is fit using $fit, but you can extract all posterior
# iterations (after warm-up) for the parameters on the response scale easily 
iters <- mod$iters()
str(iters)
# e.g. credible interval for lambda in state 1 
quantile(iters[1,], prob = c(0.025, 0.975))


# State inference ---------------------------------------------------------

# I can decode conditional on the posterior means
# First, reset parameters inside model to means
mod$update_par(iter = "mean")

# Now decode states as usual
states <- mod$viterbi()

# or sample possible state sequences conditional on parameters 
sim_states <- mod$sample_states(nsamp = 100)

# or sample state sequences and posterior samples jointly 
sim_full_states <- mod$sample_states(nsamp = 100, full = TRUE)

# Goodness-of-fit ---------------------------------------------------------

# simulated-biased testing, again you do this conditional on parameters
# currently set in model or specify full=TRUE to sample goodness-of-fit
# statistic jointly with posterior parameter values 

# define test statistic which takes a sample as input 
# Here I'll use the 20% to 80% quantiles of the data 
gof_stat <- function(x) {
  quantile(x$count, prob = seq(0.2, 0.8, 0.2))
}
# do simulation conditional on posterior means
mod$update_par(iter = "mean")
sims <- mod$gof(gof_stat, nsims = 100, silent = FALSE)

# do simulation jointly with posterior sampling
sims_full <- mod$gof(gof_stat, nsims = 100, silent = FALSE, full = TRUE)




