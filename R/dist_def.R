
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
dist_zip <- Dist$new(
  name = "zip", 
  pdf = function(x, z, lambda, log = FALSE) {
    zero <- x < 1e-10 
    l <- z * zero + (1 - z) * dpois(x, lambda)
    if (log) l <- log(l)
    return(l)
  }, 
  rng = function(n, z, lambda) {
    zero <- rbinom(n, 1, z)
    y <- rpois(n, lambda)
    y[zero == 1] <- 0
    return(y)}, 
  link = list(z = qlogis, lambda = log), 
  invlink = list(z = plogis, lambda = exp), 
  npar = 2
)

# Zero-Truncated Poisson ========================
dist_ztp <- Dist$new(
  name = "ztp", 
  pdf = function(x, lambda, log = FALSE) {
    l <- dpois(x, lambda) / (1 - dpois(0, lambda))
    if (log) l <- log(l)
    return(l)
  }, 
  rng = function(n, lambda) {
    y <- NULL 
    while (length(y) < n) {
      z <- rpois(n, lambda)
      y <- c(y, z[z > 1e-10])
    }
    return(y)
  }, 
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
dist_zib <- Dist$new(
  name = "zib", 
  pdf = function(x, z, size, p, log = FALSE) {
    zero <- x < 1e-10 
    l <- z * zero + (1 - z) * dbinom(x, size, p)
    if (log) l <- log(l)
    return(l)
  }, 
  rng = function(n, z, size, p) {
    zero <- rbinom(n, 1, z)
    y <- rbinom(n, size, p)
    y[zero == 1] <- 0
    return(y)
  }, 
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
  pdf = function(x, ..., log = TRUE) {
    # get class probabilities
    p <- c(...) 
    p <- c(1 - sum(p), p)
    if (abs(sum(p) - 1) > 1e-8) stop("class probabilities must sum to one")
    n <- round(x)
    if (n < 0 | n > length(p)) stop("invalid input")
    val <- p[n + 1]
    if (log) val <- log(val)
    return(val)
  }, 
  rng = function(n, ...) {
    # get class probabilities
    p <- c(...) 
    p <- c(1 - sum(p), p)
    if (abs(sum(p) - 1) > 1e-8) stop("class probabilities must sum to one")
    # sample classes
    samp <- sample(1:length(p), size = n, prob = p, replace = TRUE) - 1
    return(samp)
  }, 
  link = mlogit_bystates, 
  invlink = invmlogit_bystates, 
  npar = 1, 
)

# Zero-inflated Negative-Binomial ==============
dist_zinb <- Dist$new(
  name = "zinb", 
  pdf = function(x, z, size, p, log = FALSE) {
    zero <- x < 1e-10 
    l <- z * zero + (1 - z) * dnbinom(x, size, p)
    if (log) l <- log(l)
    return(l)
  }, 
  rng = function(n, z, size, p) {
    zero <- rbinom(n, 1, z)
    y <- rnbinom(n, size, p)
    y[zero == 1] <- 0
    return(y)
  }, 
  link = list(z = qlogis, size = log, p = qlogis), 
  invlink = list(z = plogis, size = exp, p = plogis), 
  npar = 3, 
)

# Zero-Truncated Negative-Binomial ========================
dist_ztnb <- Dist$new(
  name = "ztnb", 
  pdf = function(x, size, p, log = FALSE) {
    l <- dnbinom(x, size, p) / (1 - dnbinom(0, size, p))
    if (log) l <- log(l)
    return(l)
  }, 
  rng = function(n, size, p) {
    y <- NULL 
    while (length(y) < n) {
      z <- rnbinom(n, size, p)
      y <- c(y, z[z > 1e-10])
    }
    return(y)
  }, 
  link = list(size = log, p = qlogis), 
  invlink = list(size = exp, p = plogis), 
  npar = 2, 
)

# Continuous distributions ------------------------------------------------

# Normal =======================================
dist_norm <- Dist$new(
  name = "norm", 
  pdf = dnorm,
  rng = rnorm,
  link = list(mean = identity, sd = log),
  invlink = list(mean = identity, sd = exp),
  npar = 2
)


# Truncated Normal =============================
dist_truncnorm <- Dist$new(
  name = "truncnorm", 
  pdf = function(x, mean, sd, min = -Inf, max = Inf, log = FALSE) {
    left <- pnorm(min, mean, sd)
    right <- pnorm(max, mean, sd)
    p <- dnorm(x, mean, sd) / (right - left)
    if (log) p <- log(p)
    return(p)
  },
  rng = function(n, mean, sd, min = -Inf, max = Inf) {
    u <- runif(n)
    left <- pnorm(min, mean, sd)
    right <- pnorm(max, mean, sd) - pnorm(min, mean, sd)
    x <- qnorm(left + u * right) * sd + mean 
    return(x)
  },
  link = list(mean = identity, sd = log, min = identity, max = identity),
  invlink = list(mean = identity, sd = exp, min = identity, max = identity),
  npar = 4, 
  fixed = c(mean = FALSE, sd = FALSE, min = TRUE, max = TRUE)
)

# Folded Normal ================================
dist_foldednorm <- Dist$new(
  name = "foldednorm", 
  pdf = function(x, mean, sd, log = FALSE) {
    p <- dnorm(x, mean, sd) + dnorm(-x, mean, sd)
    if (log) p <- log(p)
    return(p)
  },
  rng = function(n, mean, sd) {
    x <- rnorm(n, mean, sd)
    y <- abs(x)
    return(y)
  },
  link = list(mean = identity, sd = log),
  invlink = list(mean = identity, sd = exp),
  npar = 2
)

# Student's t distribution =====================
dist_t <- Dist$new(
  name = "t", 
  pdf = function(x, mean, scale, df, log = FALSE) {
    y <- (x - mean) / scale
    return(dt(y, df = df, log = log) / scale)
  }, 
  rng = function(n, mean, scale, df) {
    y <- rt(n, df = df)
    x <- y * scale + mean 
    return(x)
  }, 
  link = list(mean = identity, scale = log, df = log), 
  invlink = list(mean = identity, scale = exp, df = exp), 
  npar = 3
)

# Log-Normal ===================================
dist_lnorm <- Dist$new(
  name = "lnorm", 
  pdf = dlnorm,
  rng = rlnorm,
  link = list(meanlog = identity, sdlog = log),
  invlink = list(meanlog = identity, sdlog = exp),
  npar = 2
)

# Gamma ========================================
dist_gamma <- Dist$new(
  name = "gamma", 
  pdf = dgamma,
  rng = rgamma,
  link = list(shape = log, scale = log),
  invlink = list(shape = exp, scale = exp),
  npar = 2
)

# Weibull ========================================
dist_weibull <- Dist$new(
  name = "weibull", 
  pdf = dweibull,
  rng = rweibull,
  link = list(shape = log, scale = log),
  invlink = list(shape = exp, scale = exp),
  npar = 2
)

# Exponential ==================================
dist_exp <- Dist$new(
  name = "exp", 
  pdf = dexp, 
  rng = rexp, 
  link = list(rate = log), 
  invlink = list(rate = exp), 
  npar = 1
)

# Beta =========================================
dist_beta <- Dist$new(
  name = "beta",
  pdf = dbeta,
  rng = rbeta,
  link = list(shape1 = log, shape2 = log),
  invlink = list(shape1 = exp, shape2 = exp),
  npar = 2
)


# Mixed distributions -----------------------------------------------------

# Tweedie ======================================
dist_tweedie <- Dist$new(
  name = "tweedie", 
  pdf = function(x, mean, p, phi, log = FALSE) {
    l <- mgcv::ldTweedie(x, mu = mean, p = p + 1, phi = phi)
    if (!log) l <- exp(l)
    return(l)
  }, 
  rng = function(n, mean, p, phi) {
    return(rTweedie(rep(mean, n), p + 1, phi))
  }, 
  link = list(mean = identity, p = qlogis, phi = log), 
  invlink = list(mean = identity, p = plogis, phi = exp), 
  npar = 3
)

# Angular distributions ---------------------------------------------------

# Von Mises ====================================
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
                  cat = dist_cat,
                  exp = dist_exp,
                  lnorm = dist_lnorm,
                  weibull = dist_weibull,
                  truncnorm = dist_truncnorm, 
                  foldednorm = dist_foldednorm, 
                  t = dist_t, 
                  tweedie = dist_tweedie
                  )

