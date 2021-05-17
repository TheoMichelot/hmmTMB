
library(hmmTMB)
set.seed(24829)

# Create weibull distribution ---------------------------------------------
weibull <- Dist$new(
  name = "Weibull2",
  pdf = dweibull,
  rng = rweibull,
  link = list(shape = log, scale = log),
  invlink = list(shape = exp, scale = exp),
  npar = 2, 
  parnames = c("shape", "scale"), 
  cpp = "inst/weibull.cpp"
)

# Simulate data -----------------------------------------------------------
n_sim <- 1000
n_states <- 2 
dist <- list(wait = "Weibull2")
par <- list(wait = list(shape = c(1, 2), scale = c(1, 10)))
dummy_data <- data.frame(ID = rep(1, n_sim), wait = rep(0, n_sim))
obs <- Observation$new(data = dummy_data, dists = dist, n_states = n_states, par = par)
par_hid0 <- c(0.05, 0.03)
hid <- MarkovChain$new(n_states = n_states)
mod <- HMM$new(obs = obs, hidden = hid)
sim <- mod$simulate(n = n_sim, data = data.frame(ID = 1, wait = rep(0, n_sim)))

# plot it 
plot(sim$wait)

# Fit model ---------------------------------------------------------------
# Initial parameters for estimation
par0 <- list(wait = list(shape = c(5, 5),
                         scale = c(10, 10)))
obs2 <- Observation$new(data = sim, dists = dist, 
                        n_states = n_states, par = par0)

hid2 <- MarkovChain$new(n_states = n_states, data = sim)

mod2 <- HMM$new(obs = obs2, hidden = hid2)

mod2$fit(silent = FALSE)

# Estimated observation parameters
mod2$obs()$par()

# Estimated parameters of the state process
mod2$hidden()$coeff_fe()

# Compare estimated states and true states
s <- mod2$viterbi()
table(s == sim$state)/n_sim
