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
  
  # Initialise state vector
  states <- NULL
  
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
    
    states <- c(states, s)
  }
  
  return(list(data = cbind(obs, X), states = states, lambda = lambda))
}

# Simulation parameters
nind <- 30
n <- 50
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
                x1 = cumsum(rnorm(n, 0, 0.05)),
                x2 = cumsum(rnorm(n, 0, 0.05)))

# Simulate HMM data
simdat <- simHMM(nind = nind, n = n, shape_par = shape_par, scale_par = scale_par, 
                 lambda_par = lambda_par, lambda_re = lambda_re, 
                 tpm = tpm, X = as.matrix(X))
data <- simdat$data
states <- simdat$state

# Observation distributions
dists <- list(step = dist_gamma, count = dist_pois)

# Initial coefficients for fixed effects
par0_shape <- c(log(1), log(4), 0, 0, 0, 0)
par0_scale <- c(log(1), log(4))
par0_lambda <- c(log(3), log(7))
coeff_fe <- c(par0_shape, par0_scale, par0_lambda)

# Formulas on observation parameters
formulas <- list(step = list(shape = ~ x1 + x2, scale = ~ 1),
                 count = list(lambda = ~ s(ID, bs = "re")))

# Create objects
dat <- HMMData$new(data)
obs <- Observation$new(dat, dists = dists, n_states = 2, coeff_fe = coeff_fe, 
                       formulas = formulas)
hid <- MarkovChain$new(n_states = 2)
mod <- HMM$new(obs, hid)

# Fit model
mod$fit(silent = FALSE)


var <- "x1"
m <- obs$make_mat_grid(var = var, n_grid = 10)
ci <- mod$CI_obspar(X_fe = m$X_fe, X_re = m$X_re)
par <- obs$par_all(X_fe = m$X_fe, X_re = m$X_re)


mod$plot_obspar(var = var)

# Unpack parameters
coeff_fe <- obs$coeff_fe()

shape_est <- matrix(coeff_fe[1:6], ncol = 2)
scale_est <- matrix(coeff_fe[7:8], ncol = 2)
lambda_est <- matrix(coeff_fe[9:10], ncol = 2)

s <- mod$viterbi()
table(s == states)/length(s)
