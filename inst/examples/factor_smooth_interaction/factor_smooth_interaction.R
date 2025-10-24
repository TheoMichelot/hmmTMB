
#' This code is just to check that models with factor-smooth interactions
#' of the form s(var, by = factor) are set up correctly. The true effect
#' is zero in the simulated data, so this doesn't check that the true
#' shape of the relationship can be recovered.
#' This does illustrate how priors can be used to mitigate numerical instability
#' (here caused by the fact that the "true" smoothness parameters are infinite)

library(moveHMM)
library(hmmTMB)

# Simulate from moveHMM
n <- 1000
data <- simData(nbAnimals = 1, 
                nbStates = 2, 
                stepPar = c(1, 10, 1, 10), 
                anglePar = c(pi, 0, 0.5, 3), 
                obsPerAnimal = n)

# Make up covariate values
data$hour <- (1:n) %% 24
data$habitat_type <- factor(sample(c("hab1", "hab2", "hab3", "hab4"), 
                                   size = n,
                                   replace = TRUE))

# Hidden state model
f <- ~ habitat_type + s(hour, by = habitat_type, bs = "cc", k = 5)
hid <- MarkovChain$new(data = data,
                       formula = f, 
                       n_states = 2)

# Initialise smoothness parameter
hid$update_lambda(rep(10, 8))

# Observation model
dists <- list(step = "gamma2", angle = "vm")
par0 <- list(step = list(mean = c(1, 10), sd = c(1, 10)),
             angle = list(mu = c(3, 0), kappa = c(0.5, 3)))
obs <- Observation$new(data = data, dists = dists, par = par0)

# Full HMM
hmm <- HMM$new(obs = obs, hid = hid)

# Add prior on smoothness parameters to avoid convergence issues
prior_smoothness <- 100
prior_log_smoothness_sd <- 1
priors <- matrix(c(log(prior_smoothness), prior_log_smoothness_sd,
                   log(prior_smoothness), prior_log_smoothness_sd,
                   log(prior_smoothness), prior_log_smoothness_sd,
                   log(prior_smoothness), prior_log_smoothness_sd,
                   log(prior_smoothness), prior_log_smoothness_sd,
                   log(prior_smoothness), prior_log_smoothness_sd,
                   log(prior_smoothness), prior_log_smoothness_sd,
                   log(prior_smoothness), prior_log_smoothness_sd),
                 nrow = 8, byrow = TRUE)
hmm$set_priors(new_priors = list(log_lambda_hid = priors))

# Fit model
hmm$fit(silent = FALSE)

# Check that estimated relationship is flat
hmm$plot("delta", var = "hour")
