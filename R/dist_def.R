
##' This file contains the distributions currently included in hmmTMB, but
##' user-specified distributions can be added with the name "custom"

# Poisson
dist_pois <- Dist$new(
  name = "pois", 
  pdf = dpois,
  link = list(lambda = log),
  invlink = list(lambda = exp),
  npar = 1
)

# gamma
dist_gamma <- Dist$new(
  name = "gamma", 
  pdf = dgamma,
  link = list(shape = log, scale = log),
  invlink = list(shape = exp, scale = exp),
  npar = 2
)

# normal
dist_norm <- Dist$new(
  name = "norm", 
  pdf = dnorm,
  link = list(mean = identity, sd = log),
  invlink = list(mean = identity, sd = exp),
  npar = 2
)

# beta
dist_beta <- Dist$new(
  name = "beta",
  pdf = dbeta,
  link = list(shape1 = log, shape2 = log),
  invlink = list(shape1 = exp, shape2 = exp),
  npar = 2
)

# von Mises
dist_vm <- Dist$new(
  name = "vm",
  pdf = function(x, mu = 0, kappa = 1, log = FALSE) {
    b <- besselI(kappa, 0)
    if(!log)
      val <- 1/(2 * pi * b) * exp(kappa * cos(x - mu))
    else
      val <- - log(2 * pi * b) + kappa * cos(x - mu)
    return(val)
  },
  link = list(mu = function(x) qlogis((x + pi) / (2 * pi)), 
              kappa = log),
  invlink = list(mu = function(x) 2 * pi * plogis(x) - pi, 
                 kappa = exp),
  npar = 2
)
