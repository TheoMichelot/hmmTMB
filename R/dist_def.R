
# This file contains the distributions currently included in hmmTMB, but
# user-specified distributions can be added with the name "custom"

# Discrete distributions --------------------------------------------------

# Poisson ======================================
dist_pois <- Dist$new(
  name = "pois", 
  pdf = dpois,
  rng = rpois,
  link = list(lambda = log),
  invlink = list(lambda = exp),
  npar = 1
)

# Zero-Inflated Poisson ========================
dzip <- function(x, z, lambda, log = FALSE) {
  zero <- x < 1e-10 
  l <- z * zero + (1 - z) * dpois(x, lambda)
  if (log) l <- log(l)
  return(l)
}
rzip <- function(n, z, lambda) {
  zero <- rbinom(n, 1, z)
  y <- rpois(n, lambda)
  y[zero == 1] <- 0
  return(y)
}
dist_zip <- Dist$new(
  name = "zip", 
  pdf = dzip, 
  rng = rzip, 
  link = list(z = qlogis, lambda = log), 
  invlink = list(z = plogis, lambda = exp), 
  npar = 2
)

# Zero-Truncated Poisson ========================
dztp <- function(x, lambda, log = FALSE) {
  l <- dpois(x, lambda) / (1 - dpois(0, lambda))
  if (log) l <- log(l)
  return(l)
}
rztp <- function(n, lambda) {
  y <- NULL 
  while (length(y) < n) {
    z <- rpois(n, lambda)
    y <- c(y, z[z > 1e-10])
  }
  return(y)
}
dist_ztp <- Dist$new(
  name = "ztp", 
  pdf = dztp, 
  rng = rztp, 
  link = list(lambda = log), 
  invlink = list(lambda = exp), 
  npar = 1
)


# Binomial ===================================== 
dist_binom <- Dist$new(
  name = "binom", 
  pdf = dbinom, 
  rng = rbinom, 
  link = list(size = identity, prob = qlogis), 
  invlink = list(size = identity, prob = plogis),
  npar = 2, 
  fixed = c(size = TRUE, prob = FALSE)
)

# Zero-inflated Binomial =======================
dzib <- function(x, z, size, p, log = FALSE) {
  zero <- x < 1e-10 
  l <- z * zero + (1 - z) * dbinom(x, size, p)
  if (log) l <- log(l)
  return(l)
}
rzib <- function(n, z, size, p) {
  zero <- rbinom(n, 1, z)
  y <- rbinom(n, size, p)
  y[zero == 1] <- 0
  return(y)
}
dist_zib <- Dist$new(
  name = "zib", 
  pdf = dzib, 
  rng = rzib, 
  link = list(z = qlogis, size = identity, p = qlogis), 
  invlink = list(z = plogis, size = identity, p = plogis), 
  npar = 3, 
  fixed = c(z = FALSE, size = TRUE, p = FALSE) 
)

# Negative-binomial ============================
dist_nbinom <- Dist$new(
  name = "nbinom", 
  pdf = dnbinom, 
  rng = rnbinom, 
  link = list(size = log, prob = qlogis), 
  invlink = list(size = exp, prob = plogis), 
  npar = 2, 
)

# Categorical ============================
dcat <- function(x, ..., log = TRUE) {
  # get class probabilities
  p <- c(...) 
  p <- c(1 - sum(p), p)
  if (abs(sum(p) - 1) > 1e-8) stop("class probabilities must sum to one")
  n <- round(x)
  if (n < 0 | n > length(p)) stop("invalid input")
  val <- p[n + 1]
  if (log) val <- log(val)
  return(val)
}
rcat <- function(n, ...) {
  # get class probabilities
  p <- c(...) 
  p <- c(1 - sum(p), p)
  if (abs(sum(p) - 1) > 1e-8) stop("class probabilities must sum to one")
  # sample classes
  samp <- sample(1:length(p), size = n, prob = p, replace = TRUE) - 1
  return(samp)
}
mlogit <- function(x) {
  s <- 1 - sum(x)
  return(log(x / s))
}
invmlogit <- function(x) {
  y <- exp(x)
  s <- 1/(1 + sum(y))
  y <- y * s
  return(y)
}
mlogit_bystates <- function(x, n_states) {
  xmat <- unlist(x)
  xmat <- matrix(xmat, nr = n_states, nc = length(xmat) / n_states)
  ymat <- t(apply(xmat, 1, mlogit))
  return(as.vector(ymat))
}
invmlogit_bystates <- function(x, n_states) {
  xmat <- matrix(x, nr = n_states, nc = length(x) / n_states)
  ymat <- t(apply(xmat, 1, invmlogit))
  return(as.vector(ymat))
}
dist_cat <- Dist$new(
  name = "cat", 
  pdf = dcat, 
  rng = rcat, 
  link = mlogit_bystates, 
  invlink = invmlogit_bystates, 
  npar = 1, 
)

# Zero-inflated Negative-Binomial ==============
dzinb <- function(x, z, size, p, log = FALSE) {
  zero <- x < 1e-10 
  l <- z * zero + (1 - z) * dnbinom(x, size, p)
  if (log) l <- log(l)
  return(l)
}
rzinb <- function(n, z, size, p) {
  zero <- rbinom(n, 1, z)
  y <- rnbinom(n, size, p)
  y[zero == 1] <- 0
  return(y)
}
dist_zinb <- Dist$new(
  name = "zinb", 
  pdf = dzinb, 
  rng = rzinb, 
  link = list(z = qlogis, size = log, p = qlogis), 
  invlink = list(z = plogis, size = exp, p = plogis), 
  npar = 3, 
)

# Zero-Truncated Negative-Binomial ========================
dztnb <- function(x, size, p, log = FALSE) {
  l <- dnbinom(x, size, p) / (1 - dnbinom(0, size, p))
  if (log) l <- log(l)
  return(l)
}
rztnb <- function(n, size, p) {
  y <- NULL 
  while (length(y) < n) {
    z <- rnbinom(n, size, p)
    y <- c(y, z[z > 1e-10])
  }
  return(y)
}
dist_ztnb <- Dist$new(
  name = "ztnb", 
  pdf = dztnb, 
  rng = rztnb, 
  link = list(size = log, p = qlogis), 
  invlink = list(size = exp, p = plogis), 
  npar = 2, 
)

# Continuous distributions ------------------------------------------------

# gamma
dist_gamma <- Dist$new(
  name = "gamma", 
  pdf = dgamma,
  rng = rgamma,
  link = list(shape = log, scale = log),
  invlink = list(shape = exp, scale = exp),
  npar = 2
)

# normal
dist_norm <- Dist$new(
  name = "norm", 
  pdf = dnorm,
  rng = rnorm,
  link = list(mean = identity, sd = log),
  invlink = list(mean = identity, sd = exp),
  npar = 2
)

# beta
dist_beta <- Dist$new(
  name = "beta",
  pdf = dbeta,
  rng = rbeta,
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
  rng = function(n, mu, kappa) {
    # rvm and dvm use different parameter names
    # (also, translate rvm output from [0, 2pi] to [-pi, pi])
    rvm(n = n, mean = mu, k = kappa) - pi
  },
  link = list(mu = function(x) qlogis((x + pi) / (2 * pi)), 
              kappa = log),
  invlink = list(mu = function(x) 2 * pi * plogis(x) - pi, 
                 kappa = exp),
  npar = 2
)


# List of all distributions -----------------------------------------------

dist_list <- list(pois = dist_pois,
                  binom = dist_binom,
                  nbinom = dist_nbinom, 
                  gamma = dist_gamma,
                  norm = dist_norm,
                  beta = dist_beta,
                  vm = dist_vm, 
                  zip = dist_zip, 
                  zib = dist_zib, 
                  zinb = dist_zinb,
                  ztp = dist_ztp, 
                  ztnb = dist_ztnb, 
                  cat = dist_cat)

