## Driver script

library(HmmTmb)

#' Simulate gamma HMM
#'
#' @param n number of samples
#' @param lambda rate of Poisson for each state
#' @param shape shape of gamma for each state
#' @param scale scale of gamma for each state
#' @param tpm transition probability matrix
#' @param n.states number of states
#'
#' @return vector of observed counts
simHMM <- function(n, lambda, shape, scale, tpm) {
  n.states <- length(lambda)
  
  # Initial distribution (stationary)
  delta <- solve(t(diag(n.states) - tpm  + 1), rep(1, n.states))
  
  # Simulate state process
  s <- numeric(n)
  state.space <- 1:n.states
  s[1] <- sample(state.space, 1, prob = delta)
  for (t in 2:n) 
    s[t] <- sample(state.space, 1, prob = tpm[s[t - 1], ])
  
  # Simulate observations
  data <- data.frame(step = rgamma(n = n, shape = shape[s], scale = scale[s]),
                     count = rpois(n = n, lambda = lambda[s]))
  
  return(data)
}

n <- 10000
shape <- c(0.5, 6)
scale <- c(2, 3)
lambda <- c(5, 10)
tpm <- matrix(c(0.9, 0.3, 0.1, 0.7), nc = 2)
n.states <- 2

simdat <- simHMM(n = n, lambda = lambda, shape = shape, scale = scale, 
                 tpm = tpm)

plot.ts(simdat)

dist_pois <- Dist$new(name = "pois", pdf = dpois,
                      link = list(lambda = log),
                      invlink = list(lambda = exp))

dist_gamma <- Dist$new(name = "custom", pdf = dgamma,
                       link = list(shape = log, scale = log),
                       invlink = list(shape = exp, scale = exp))

# create objects
dat <- HmmData$new(simdat)
dists <- list(step = dist_gamma, count = dist_pois)
par <- list(step = list(shape = c(1, 3), scale = c(1, 2)),
            count = list(lambda = c(3, 6)))
obs <- Observation$new(dat, dists = dists, par = par)
hid <- MarkovChain$new(matrix(c(".", "~1", "~1", "."), nr = 2),
                       matrix(c(0.8, 0.2, 0.2, 0.8), nr = 2))
mod <- Hmm$new(obs, hid)

#fit model
mod$fit()

# get parameter estimates
mod$est()
