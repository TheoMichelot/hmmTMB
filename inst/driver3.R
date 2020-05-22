# Driver script (with random effects)

library(hmmTMB)
set.seed(6932)

#' Simulate HMM
#'
#' @param n number of samples
#' @param shape_par parameters for gamma shape (matrix with
#' one row for each column of X, and one column for each state)
#' @param scale_par parameters of gamma scale for each state
#' @param lambda_par parameters of Poisson rate for each state
#' @param tpm transition probability matrix
#' @param X design matrix (covariates)
#'
#' @return vector of observed counts
simHMM <- function(nind, n, shape_par, scale_par, lambda_par, lambda_re, tpm, X) {
  n.states <- ncol(lambda_par)
  
  # Initial distribution (stationary)
  delta <- solve(t(diag(n.states) - tpm  + 1), rep(1, n.states))
  
  # Random intercepts
  lambda_ind <- rnorm(nind, 0, lambda_re)
  
  # Initialise data frame
  obs <- NULL
  
  # Loop over individuals
  for(zoo in 1:nind) {
    # Simulate state process
    s <- numeric(n)
    state.space <- 1:n.states
    s[1] <- sample(state.space, 1, prob = delta)
    for (t in 2:n) 
      s[t] <- sample(state.space, 1, prob = tpm[s[t - 1], ])
    
    # Get observation parameters  
    shape <- exp(X %*% shape_par)
    scale <- exp(X %*% scale_par)
    lambda <- exp(X %*% lambda_par + lambda_ind[zoo])
    
    # Simulate observations
    step <- rep(NA, n)
    count <- rep(NA, n)
    for(t in 1:n) {
      step[t] <- rgamma(1, shape = shape[t, s[t]], scale = scale[t, s[t]])
      count[t] <- rpois(1, lambda = lambda[t, s[t]])
    }
    
    obs <- rbind(obs, data.frame(ID = factor(zoo),
                                 step = step,
                                 count = count))
  }
  
  return(list(data = cbind(obs, X), state = s))
}

# Simulation parameters
nind <- 50
n <- 200
shape_par <- matrix(c(log(0.5), log(6),
                      -0.1, -0.1,
                      0, -0.5),
                    ncol = 2, byrow = TRUE)
scale_par <- matrix(c(log(2), log(3),
                      0, 0,
                      0, 0),
                    ncol = 2, byrow = TRUE)
lambda_par <- matrix(c(log(2), log(20),
                       0, 0,
                       0, 0),
                     ncol = 2, byrow = TRUE)
lambda_re <- 0.3 # SD of random effect distribution
tpm <- matrix(c(0.95, 0.1, 0.05, 0.9), nc = 2)

# Simulate covariates
X <- data.frame(Intercept = 1,
                x1 = rnorm(n),
                x2 = rnorm(n))

# Simulate HMM data
simdat <- simHMM(nind = nind, n = n, shape_par = shape_par, scale_par = scale_par, 
                 lambda_par = lambda_par, lambda_re = lambda_re, 
                 tpm = tpm, X = as.matrix(X))
data <- simdat$data
states <- simdat$states

# Observation distributions
dists <- list(step = dist_gamma, count = dist_pois)

# Initial parameters (working scale)
par0_shape <- c(log(1), log(4), 0, 0, 0, 0)
par0_scale <- c(log(1), log(4))
par0_lambda <- c(log(3), log(7))
wpar_fe <- c(par0_shape, par0_scale, par0_lambda)
wpar_re <- rep(0, 2*nind)

# Formulas on observation parameters
formulas <- list(step = list(shape = ~ x1 + x2, scale = ~ 1),
                 count = list(lambda = ~ s(ID, bs = "re")))

# Create objects
dat <- HmmData$new(data)
obs <- Observation$new(dat, dists = dists, n_states = 2, wpar = wpar_fe, 
                       wpar_re = wpar_re, formulas = formulas)
hid <- MarkovChain$new(matrix(c(".", "~1", "~1", "."), nr = 2),
                       matrix(c(0.8, 0.2, 0.2, 0.8), nr = 2))
mod <- Hmm$new(obs, hid)

# Fit model
mod$fit()

# Unpack parameters
wpar <- obs$wpar()

shape_est <- matrix(wpar[1:6], ncol = 2)
scale_est <- matrix(wpar[7:8], ncol = 2)
lambda_est <- matrix(wpar[9:10], ncol = 2)

# Random effect parameter
1/sqrt(exp(mod$res()$par["log_lambda_obs"]))

s <- mod$viterbi()
length(which(states == s))/length(states)