
library(hmmTMB)
set.seed(5738)

###################
## Simulate data ##
###################
# Generate two covariates (random walk)
n_sim <- 1e4
covs <- data.frame(ID = 1, 
                   step = NA,
                   count = NA,
                   x1 = cumsum(rnorm(n_sim, 0, 0.1)),
                   x2 = cumsum(rnorm(n_sim, 0, 0.1)))

# Create observation model
n_states <- 2
dists <- list(step = dist_gamma,
              count = dist_pois)
par <- list(step = list(shape = c(1, 1),
                        scale = c(1, 10)),
            count = list(lambda = c(1, 10)))
obs <- Observation$new(data = covs, dists = dists, 
                       n_states = n_states, par = par)

# Create state process model
struct <- matrix(c(".", "~ x1 + x2",
                   "~ x2", "."), 
                 ncol = 2, byrow = TRUE)
par_hid0 <- c(-1.5, 0.2, -0.3, -3, 0.5)
hid <- MarkovChain$new(n_states = 2, formula = struct, 
                       data = covs)
hid$update_coeff_fe(par_hid0)

# Create HMM object and simulate data
mod <- HMM$new(obs = obs, hidden = hid)
sim <- mod$simulate(n = n_sim, data = covs)

###############
## Fit model ##
###############
# Initial parameters for estimation
par0 <- list(step = list(shape = c(0.5, 2),
                         scale = c(2, 8)),
             count = list(lambda = c(3, 7)))
obs2 <- Observation$new(data = sim, dists = dists, 
                        n_states = n_states, par = par)

hid2 <- MarkovChain$new(n_states = 2, formula = struct, data = sim)

mod2 <- HMM$new(obs = obs2, hidden = hid2)

mod2$fit(silent = FALSE)

# Estimated observation parameters
mod2$obs()$par()

# Estimated parameters of the state process
mod2$hidden()$coeff_fe()

# Compare estimated states and true states
s <- mod2$viterbi()
table(s == attr(sim, "state"))/n_sim
