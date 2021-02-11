# Simple Poisson HMM example of using hmmTMB and maximum likelihood

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

# Fit model ------------------------------------------------------------

# create model to fit 
mod <- HMM$new(file = "pois_mod.hmm")

# check model formulation
mod

# suggest better starting parameters?
ini <- mod$suggest_initial()

# set to new starting parameters (or you could edit the specification file)
mod$obs()$update_par(ini)

# fit model
mod$fit()

# Parameter Inference -----------------------------------------------------

# look at numerical estimates
mod$par()

# get uncertainty, e.g., for count means 
mod$predict("obspar", level = 0.95)
# or for transition probabilities
mod$predict("tpm", level = 0.95)
# or for stationary distribution 
mod$predict("delta", level = 0.95)

# State inference ---------------------------------------------------------

# decode states
states <- mod$viterbi()

# plot decoded states with data
mod$plot_ts("count")

# sample possible state sequences 
sim_states <- mod$sample_states(nsamp = 100)
plot(sim_states[,1], pch = 20, col = 1 + (sim_states[,1] < 1.5))

# get state probabilities 
state_probs <- mod$state_probs()
summary(state_probs)

# Goodness-of-fit ---------------------------------------------------------

# get pseudo-residuals, should be normally distributed 
resids <- mod$pseudores()
hist(resids)
qqnorm(resids)
qqline(resids) # bumps because count is discrete random variable 

# simulated-biased testing
# define test statistic which takes a sample as input 
# Here I'll use the 20% to 80% quantiles of the data 
gof_stat <- function(x) {
  quantile(x$count, prob = seq(0.2, 0.8, 0.2))
}
# do simulation 
sims <- mod$gof(gof_stat, nsims = 100, silent = FALSE)
# Could look at other aspects, e.g. proportion of time above a threshold
thres_stat <- function(x) {
  sum(x$count > 10)
}
thres_sims <- mod$gof(thres_stat, nsims = 100, silent = FALSE)




