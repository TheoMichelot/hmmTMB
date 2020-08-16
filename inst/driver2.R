
library(hmmTMB)
set.seed(5738)

###################
## Simulate data ##
###################
# Generate two covariates (random walk)
n_sim <- 1e4
covs <- data.frame(ID = 1, 
                   x1 = cumsum(rnorm(n_sim, 0, 0.1)),
                   x2 = cumsum(rnorm(n_sim, 0, 0.1)))
hmm_data <- HmmData$new(data = covs)

# Create observation model
n_states <- 2
dists <- list(step = dist_gamma,
              count = dist_pois)
par <- list(step = list(shape = c(1, 1),
                        scale = c(1, 10)),
            count = list(lambda = c(1, 10)))
obs <- Observation$new(data = hmm_data, dists = dists, 
                       n_states = n_states, par = par)

# Create state process model
struct <- matrix(c(".", "~ x1 + x2",
                   "~ x2", "."), 
                 ncol = 2, byrow = TRUE)
par_hid0 <- c(-1.5, 0.2, -0.3, -3, 0.5)
hid <- MarkovChain$new(n_states = 2, structure = struct, 
                       coeff_fe0 = par_hid0, data = hmm_data)

# Create HMM object and simulate data
mod <- Hmm$new(obs = obs, hidden = hid)
sim <- mod$simulate(n = n_sim, data = covs)

###############
## Fit model ##
###############
hmm_data2 <- HmmData$new(data = sim)

# Initial parameters for estimation
par0 <- list(step = list(shape = c(0.5, 2),
                         scale = c(2, 8)),
             count = list(lambda = c(3, 7)))
obs2 <- Observation$new(data = hmm_data2, dists = dists, 
                        n_states = n_states, par = par)

hid2 <- MarkovChain$new(n_states = 2, structure = struct, data = hmm_data)

mod2 <- Hmm$new(obs = obs2, hidden = hid2)

mod2$fit(silent = FALSE)

# Estimated observation parameters
mod2$obs()$par()

# Estimated parameters of the state process
mod2$hidden()$coeff_fe()

# Compare estimated states and true states
s <- mod2$viterbi()
table(s == sim$state)/n_sim

############################################
## Experiment: uncertainty quantification ##
############################################
rep <- mod2$tmb_rep()

par <- rep$par.fixed
V <- rep$cov.fixed

n_post <- 1e3
post <- rmvn(n = n_post, mu = par, V = V)

par_fe_obs <- post[, which(colnames(post) == "wpar_fe_obs")]
par_re_obs <- post[, which(colnames(post) == "wpar_re_obs")]

m <- obs$make_mat(new_data = data.frame(ID = 1))

foo <- sapply(1:n_post, function(i) {
  obs$par_all(X_fe = m$X_fe, X_re = m$X_re, 
              par_fe = par_fe_obs[i,], par_re = par_re_obs[i,])
})
