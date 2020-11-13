## Driver script

library(hmmTMB)

#' Simulate Poisson HMM
#'
#' @param n number of samples
#' @param lambda means of Poisson for each state
#' @param tpm transition probability matrix
#' @param n.states number of states
#' @param delta initial distribution
#'
#' @return vector of observed counts
SimulatePoHmm <- function(n, lambda, tpm, n.states) {
  delta <- solve(t(diag(n.states) - tpm  + 1), rep(1, n.states))
  s <- numeric(n)
  state.space <- 1:n.states
  s <- sample(state.space, 1, prob = delta)
  for (t in 2:n) s[t] <- sample(state.space, 1, prob = tpm[s[t - 1], ])
  data <- rpois(n, lambda = lambda[s])
  return(list(data = data, state = s))
}

n <- 10000
lambda <- c(2, 20)
tpm <- matrix(c(0.95, 0.1, 0.05, 0.9), nc = 2)
n.states <- 2

simdat <- SimulatePoHmm(n, lambda, tpm, n.states)
counts <- simdat$data
states <- simdat$state

# create objects
dat <- data.frame(ID = 1,
                  count = counts)
dists <- list(count = "pois")
par <- list(count = list(lambda = c(3, 6)))
obs <- Observation$new(data = dat, dists = dists, n_states = 2, par = par)
hid <- MarkovChain$new(n_states = 2)
mod <- HMM$new(obs, hid)

#fit model
mod$fit(silent = FALSE)
mod$par()

s <- mod$viterbi()
table(s == states)/length(s)
