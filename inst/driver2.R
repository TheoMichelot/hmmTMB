## Driver script

library(HmmTmb)

#' Simulate gamma HMM
#'
#' @param n number of samples
#' @param shape shape of gamma for each state
#' @param scale scale of gamma for each state
#' @param tpm transition probability matrix
#' @param n.states number of states
#'
#' @return vector of observed counts
SimulateGamHmm <- function(n, shape, scale, tpm, n.states) {
  delta <- solve(t(diag(n.states) - tpm  + 1), rep(1, n.states))
  s <- numeric(n)
  state.space <- 1:n.states
  s <- sample(state.space, 1, prob = delta)
  for (t in 2:n) 
    s[t] <- sample(state.space, 1, prob = tpm[s[t - 1], ])
  data <- rgamma(n = n, shape = shape[s], scale = scale[s])
  return(data)
}

n <- 10000
shape <- c(0.5, 6)
scale <- c(2, 3)
tpm <- matrix(c(0.9, 0.3, 0.1, 0.7), nc = 2)
n.states <- 2

simdat <- SimulateGamHmm(n = n, shape = shape, scale = scale, 
                         tpm = tpm, n.states = n.states)

plot(simdat, type = "b", pch = 19)

dist_gamma <- Dist$new(name = "custom", pdf = dgamma,
                       link = list(shape = log, scale = log),
                       invlink = list(shape = exp, scale = exp))

# create objects
dat <- HmmData$new(data.frame(step = simdat))
dists <- list(step = dist_gamma)
par <- list(step = list(shape = c(1, 3), scale = c(1, 2)))
obs <- Observation$new(dat, dists = dists, par = par)
hid <- MarkovChain$new(matrix(c(".", "~1", "~1", "."), nr = 2),
                       matrix(c(0.8, 0.2, 0.2, 0.8), nr = 2))
mod <- Hmm$new(obs, hid)

#fit model
mod$fit()

ltpm <- mod$res()$par[1:2]
llam <- mod$res()$par[-(1:2)]

exp(llam)
exp(ltpm) / (1 + exp(ltpm))
