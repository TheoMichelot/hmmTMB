## Driver script

library(HmmTmb)

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
  return(data)
}

n <- 10000
lambda <- c(5, 10)
tpm <- matrix(c(0.9, 0.3, 0.1, 0.7), nc = 2)
n.states <- 2

simdat <- SimulatePoHmm(n, lambda, tpm, n.states)

plot(simdat, type = "b", pch = 19)

dist_pois <- Dist$new(name = "pois", pdf = dpois,
                      link = list(lambda = log),
                      invlink = list(lambda = exp))

# create objects
dat <- HmmData$new(data.frame(count = simdat))
dists <- list(count = dist_pois)
par <- list(count = list(lambda = c(3, 6)))
obs <- Observation$new(dat, dists = dists, par = par)
hid <- MarkovChain$new(matrix(c(".", "~1", "~1", "."), nr = 2),
                       matrix(c(0.8, 0.2, 0.2, 0.8), nr = 2))
mod <- Hmm$new(obs, hid)

#fit model
mod$fit()

ltpm <- mod$res()$par[1:2]
llam <- mod$res()$par[3:4]

exp(llam)
exp(ltpm) / (1 + exp(ltpm))
