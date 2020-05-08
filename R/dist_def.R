
##' This file contains the distributions currently included in hmmTMB, but
##' user-specified distributions can be added with the name "custom"

# Poisson
dist_pois <- Dist$new(
  name = "pois", 
  pdf = dpois,
  link = list(lambda = log),
  invlink = list(lambda = log),
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
