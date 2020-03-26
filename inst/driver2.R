## Driver script (with covariates)

library(HmmTmb)
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
simHMM <- function(n, shape_par, scale_par, lambda_par, tpm, X) {
  n.states <- ncol(lambda_par)
  
  # Initial distribution (stationary)
  delta <- solve(t(diag(n.states) - tpm  + 1), rep(1, n.states))
  
  # Simulate state process
  s <- numeric(n)
  state.space <- 1:n.states
  s[1] <- sample(state.space, 1, prob = delta)
  for (t in 2:n) 
    s[t] <- sample(state.space, 1, prob = tpm[s[t - 1], ])
  
  # Get observation parameters  
  shape <- exp(X %*% shape_par)
  scale <- exp(X %*% scale_par)
  lambda <- exp(X %*% lambda_par)
  
  # Simulate observations
  step <- rep(NA, n)
  count <- rep(NA, n)
  for(t in 1:n) {
    step[t] <- rgamma(1, shape = shape[t, s[t]], scale = scale[t, s[t]])
    count[t] <- rpois(1, lambda = lambda[t, s[t]])
  }
  
  obs <- data.frame(step = step, count = count)
  return(cbind(obs, X))
}

# Simulation parameters
n <- 1e3
shape_par <- matrix(c(log(0.5), log(6),
                      -0.1, -0.1,
                      0, 0),
                    ncol = 2, byrow = TRUE)
scale_par <- matrix(c(log(2), log(3),
                      0, 0,
                      0, 0),
                    ncol = 2, byrow = TRUE)
lambda_par <- matrix(c(log(5), log(10),
                       0.5, 0.3,
                       0, -0.5),
                     ncol = 2, byrow = TRUE)
tpm <- matrix(c(0.9, 0.3, 0.1, 0.7), nc = 2)

# Simulate covariates
X <- data.frame(Intercept = 1,
                x1 = rnorm(n),
                x2 = rnorm(n))

# Simulate HMM data
simdat <- simHMM(n = n, shape_par = shape_par, scale_par = scale_par, 
                 lambda_par = lambda_par, tpm = tpm, X = as.matrix(X))

# Observation distributions
dists <- list(step = dist_gamma, count = dist_pois)

# Initial parameters (working scale)
wpar <- rep(0, 12)

# Formulas on observation parameters
formulas <- list(step = list(shape = ~ x1, scale = ~ 1),
                 count = list(lambda = ~ x1 + x2))

# Create objects
dat <- HmmData$new(simdat)
obs <- Observation$new(dat, dists = dists, wpar = wpar, 
                       formulas = formulas)
hid <- MarkovChain$new(matrix(c(".", "~1", "~1", "."), nr = 2),
                       matrix(c(0.8, 0.2, 0.2, 0.8), nr = 2))
mod <- Hmm$new(obs, hid)

# Fit model
mod$fit()

# Unpack parameters
wpar <- obs$tpar()
shape_est <- matrix(wpar[1:4], ncol = 2)
scale_est <- matrix(wpar[5:6], ncol = 2)
lambda_est <- matrix(wpar[7:12], ncol = 2)
