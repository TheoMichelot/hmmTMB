
# This file contains the distributions currently included in hmmTMB

# Discrete distributions --------------------------------------------------

# Poisson ======================================
dist_pois <- Dist$new(
  name = "pois", 
  name_long = "Poisson", 
  pdf = function(x, rate, log = FALSE) {
    dpois(x = x, lambda = rate, log = log)
  },
  cdf = function(q, rate) {
    ppois(q = q, lambda = rate)
  },
  rng = function(n, rate) {
    rpois(n = n, lambda = rate)
  },
  link = list(rate = log),
  invlink = list(rate = exp),
  npar = 1, 
  parnames = c("rate"), 
  parapprox = function(x) {
    return(mean(x))
  }
)

# Zero-Inflated Poisson ========================
dist_zipois <- Dist$new(
  name = "zipois", 
  name_long = "zero-inflated Poisson",
  pdf = function(x, rate, z, log = FALSE) {
    zero <- x == 0 
    l <- z * zero + (1 - z) * dpois(x, rate)
    if (log) l <- log(l)
    return(l)
  }, 
  cdf = function(q, rate, z) {
    return(NA)
  },
  rng = function(n, rate, z) {
    zero <- rbinom(n, 1, z)
    y <- rpois(n, rate)
    y[zero == 1] <- 0
    return(y)}, 
  link = list(rate = log, z = qlogis), 
  invlink = list(rate = exp, z = plogis), 
  npar = 2, 
  parnames = c("rate", "z"), 
  parapprox = function(x) {
    # Beckett et al. (2014). Zero-inflated Poisson (ZIP) distribution: 
    # parameter estimation and applications to model data from natural 
    # calamities. Involve, a Journal of Mathematics, 7(6), 751-767.
    mu <- mean(x)
    s2 <- var(x)
    if (mu > s2) {
      z <- 1e-3
      rate <- mu 
    } else {
      rate <- mu + s2 / mu - 1 
      z <- (s2 - mu)/(mu^2 + s2 - mu)
    }
    return(c(rate, z))
  }
)

# Zero-Truncated Poisson ========================
dist_ztpois <- Dist$new(
  name = "ztpois", 
  name_long = "zero-truncated Poisson",
  pdf = function(x, rate, log = FALSE) {
    l <- dpois(x, rate) / (1 - dpois(0, rate))
    if (log) l <- log(l)
    return(l)
  }, 
  cdf = function(q, rate) {
    return(NA)
  },
  rng = function(n, rate) {
    y <- NULL 
    while (length(y) < n) {
      z <- rpois(n, rate)
      y <- c(y, z[z > 1e-10])
    }
    return(y)
  }, 
  link = list(rate = log), 
  invlink = list(rate = exp), 
  npar = 1, 
  parnames = c("rate"), 
  parapprox = function(x) {
    warning("approximating parameters ignores zero truncation")
    return(mean(x))
  }
)


# Binomial ===================================== 
dist_binom <- Dist$new(
  name = "binom", 
  name_long = "binomial",
  pdf = dbinom, 
  cdf = pbinom,
  rng = rbinom, 
  link = list(size = identity, prob = qlogis), 
  invlink = list(size = identity, prob = plogis),
  npar = 2, 
  parnames = c("size", "prob"), 
  fixed = c(size = TRUE, prob = FALSE), 
  parapprox = function(x, size = 1) {
    p <- sum(x) / (size * length(x))
    return(c(size, p))
  }
)

# Zero-inflated Binomial =======================
dist_zibinom <- Dist$new(
  name = "zibinom", 
  name_long = "zero-inflated binomial",
  pdf = function(x, size, prob, z, log = FALSE) {
    zero <- x == 0 
    l <- z * zero + (1 - z) * dbinom(x, size, prob)
    if (log) l <- log(l)
    return(l)
  }, 
  cdf = function(q, size, prob, z) {
    return(NA)
  },
  rng = function(n, size, prob, z) {
    zero <- rbinom(n, 1, z)
    y <- rbinom(n, size, prob)
    y[zero == 1] <- 0
    return(y)
  }, 
  link = list(size = identity, prob = qlogis, z = qlogis), 
  invlink = list(size = identity, prob = plogis, z = plogis), 
  npar = 3, 
  parnames = c("size", "prob", "z"), 
  fixed = c(size = TRUE, prob = FALSE, z = FALSE), 
  parapprox = function(x, size = 1) {
    # approximates Binomial with Poisson 
    mu <- mean(x)
    s2 <- var(x)
    if (mu > s2) {
      z <- 1e-3
      rate <- mu 
    } else {
      rate <- mu + s2 / mu - 1 
      z <- s2 / mu - 1 / rate 
    }
    prob <- rate / size 
    return(c(size, prob, z))
  }
)

# Negative-binomial ============================
dist_nbinom <- Dist$new(
  name = "nbinom", 
  name_long = "negative binomial",
  pdf = dnbinom, 
  cdf = pnbinom,
  rng = rnbinom, 
  link = list(size = log, prob = qlogis), 
  invlink = list(size = exp, prob = plogis), 
  npar = 2, 
  parnames = c("size", "prob"), 
  parapprox = function(x) {
    mean <- mean(x)
    var <- var(x)
    if(var <= mean) {
      # needs overdispersion
      var <- 1.01 * mean
    }
    size <- mean^2 / (var - mean)
    prob <- mean / var
    return(c(size, prob))
  }
)

dist_nbinom2 <- Dist$new(
  name = "nbinom2", 
  name_long = "negative binomial",
  pdf = function(x, mean, shape, log = FALSE) {
    size <- shape
    prob <- shape / (mean + shape)
    dnbinom(x = x, size = size, prob = prob, log = log)
  }, 
  cdf = function(q, mean, shape) {
    size <- shape
    prob <- shape / (mean + shape)
    pnbinom(q = q, size = size, prob = prob)
  },
  rng = function(n, mean, shape) {
    size <- shape
    prob <- shape / (mean + shape)
    rnbinom(n = n, size = size, prob = prob)
  }, 
  link = list(mean = log, shape = log), 
  invlink = list(mean = exp, shape = exp), 
  npar = 2, 
  parnames = c("mean", "shape"), 
  parapprox = function(x) {
    mean <- mean(x)
    var <- var(x)
    # needs overdispersion
    shape <- ifelse(mean < var, 
                    yes = mean^2/(var-mean), 
                    no = 10000)
    return(c(mean, shape))
  }
)

# Categorical ============================
dist_cat <- Dist$new(
  name = "cat", 
  name_long = "categorical",
  pdf = function(x, ..., log = TRUE) {
    # get class probabilities
    p <- c(...) 
    p <- c(1 - sum(p), p)
    if (p[1] < -1e-10) stop("class probabilities must sum to one")
    n <- round(x)
    if (n < 1 | n > length(p)) stop("invalid input")
    val <- p[n]
    if (log) val <- log(val)
    return(val)
  }, 
  cdf = function(q, ...) {
    return(NA)
  },
  rng = function(n, ...) {
    # get class probabilities
    p <- c(...) 
    p <- c(1 - sum(p), p)
    if (p[1] < -1e-10) stop("class probabilities must sum to one")
    # sample classes
    samp <- sample(1:length(p), size = n, prob = p, replace = TRUE)
    return(samp)
  }, 
  link = function(x, n_states) {
    xmat <- unlist(x)
    xmat <- matrix(xmat, nr = n_states, nc = length(xmat) / n_states)
    ymat <- t(apply(xmat, 1, mlogit))
    return(as.vector(ymat))
  }, 
  invlink = function(x, n_states) {
    xmat <- matrix(x, nr = n_states, nc = length(x) / n_states)
    ymat <- t(apply(xmat, 1, invmlogit))
    return(as.vector(ymat))
  }, 
  # npar and parnames are updated based on data values in 
  # Observation$setup_cat()
  npar = 3,
  parnames = c("p2", "p3", "..."), 
  parapprox = function(x) {
    stop("parapprox() not implemented yet for categorical distribution")
    # # This would require parapprox knowing how many categories there are
    # p <- rep(0, length = npar + 1)
    # tab <- prop.table(table(x))
    # p[as.numeric(names(tab))] <- as.numeric(tab)
    # return(p[-1])
  },
  par_alt = function(par) {
    p <- c(1 - sum(par), par)
    names(p) <- paste0("p", 1:length(p))
    return(p)
  }
)

# Zero-inflated Negative-Binomial ==============
dist_zinbinom <- Dist$new(
  name = "zinbinom", 
  name_long = "zero-inflated negative binomial",
  pdf = function(x, size, prob, z, log = FALSE) {
    zero <- x == 0 
    l <- z * zero + (1 - z) * dnbinom(x, size, prob)
    if (log) l <- log(l)
    return(l)
  }, 
  cdf = function(q, size, prob, z) {
    return(NA)
  },
  rng = function(n, size, prob, z) {
    zero <- rbinom(n, 1, z)
    y <- rnbinom(n, size, prob)
    y[zero == 1] <- 0
    return(y)
  }, 
  link = list(size = log, prob = qlogis, z = qlogis), 
  invlink = list(size = exp, prob = plogis, z = plogis), 
  npar = 3, 
  parnames = c("size", "prob", "z"), 
  parapprox = function(x, size = 1) {
    # approximates Negative binomial with Poisson 
    mu <- mean(x)
    s2 <- var(x)
    if (mu > s2) {
      z <- 1e-3
      rate <- mu 
    } else {
      rate <- mu + s2 / mu - 1 
      z <- s2 / mu - 1 / rate 
    }
    # use minimum variance unbiased estimate of p
    n <- length(x)
    k <- rate * n
    prob <- (n * size - 1) / (n * size + k - 1)
    return(c(size, prob, z))
  }
)

# Zero-Truncated Negative-Binomial ========================
dist_ztnbinom <- Dist$new(
  name = "ztnbinom",
  name_long = "zero-truncated negative binomial",
  pdf = function(x, size, prob, log = FALSE) {
    l <- dnbinom(x, size, prob) / (1 - dnbinom(0, size, prob))
    if (log) l <- log(l)
    return(l)
  }, 
  cdf = function(q, size, prob) {
    return(NA)
  },
  rng = function(n, size, prob) {
    y <- NULL 
    while (length(y) < n) {
      z <- rnbinom(n, size, prob)
      y <- c(y, z[z > 1e-10])
    }
    return(y)
  }, 
  link = list(size = log, prob = qlogis), 
  invlink = list(size = exp, prob = plogis), 
  npar = 2,
  parnames = c("size", "prob"), 
  parapprox = function(x, size = 1) {
    warning("parameter approximation ignores zero truncation")
    # use minimum variance unbiased estimate of p
    k <- sum(x)
    n <- length(x)
    prob <- (n * size - 1) / (n * size + k - 1)
    return(c(size, prob))
  }
)

# Continuous distributions ------------------------------------------------

# Normal =======================================
dist_norm <- Dist$new(
  name = "norm", 
  name_long = "normal",
  pdf = dnorm,
  cdf = pnorm,
  rng = rnorm,
  link = list(mean = identity, sd = log),
  invlink = list(mean = identity, sd = exp),
  npar = 2, 
  parnames = c("mean", "sd"), 
  parapprox = function(x) {
    return(c(mean(x), sd(x)))
  }
)


# Truncated Normal =============================
dist_truncnorm <- Dist$new(
  name = "truncnorm", 
  name_long = "truncated normal",
  pdf = function(x, mean, sd, min = -Inf, max = Inf, log = FALSE) {
    left <- pnorm((min - mean) / sd)
    right <- pnorm((max - mean) / sd)
    p <- ifelse(x >= min & x <= max, 
                dnorm((x - mean) / sd) / (right - left),
                0)
    p <- p / sd
    if (log) p <- log(p)
    return(p)
  },
  cdf = function(q, mean, sd, min, max) {
    return(NA)
  },
  rng = function(n, mean, sd, min = -Inf, max = Inf) {
    u <- runif(n)
    left <- pnorm((min - mean) / sd)
    right <- pnorm((max - mean) / sd) - pnorm((min - mean) / sd)
    x <- qnorm(left + u * right) * sd + mean 
    return(x)
  },
  link = list(mean = identity, sd = log, min = identity, max = identity),
  invlink = list(mean = identity, sd = exp, min = identity, max = identity),
  npar = 4, 
  parnames = c("mean", "sd", "min", "max"), 
  parapprox = function(x) {
    warning("parameter approximations ignore truncations in distributions")
    return(c(mean(x), sd(x)))
  }, 
  fixed = c(mean = FALSE, sd = FALSE, min = TRUE, max = TRUE)
)

# Folded Normal ================================
dist_foldednorm <- Dist$new(
  name = "foldednorm", 
  name_long = "folded normal",
  pdf = function(x, mean, sd, log = FALSE) {
    p <- dnorm(x, mean, sd) + dnorm(-x, mean, sd)
    if (log) p <- log(p)
    return(p)
  },
  cdf = function(q, mean, sd) {
    return(NA)
  },
  rng = function(n, mean, sd) {
    x <- rnorm(n, mean, sd)
    y <- abs(x)
    return(y)
  },
  link = list(mean = identity, sd = log),
  invlink = list(mean = identity, sd = exp),
  npar = 2, 
  parnames = c("mean", "sd"), 
  parapprox = function(x) {
    warning("parameter approximations ignore folded distributions")
    return(c(mean(x), sd(x)))
  }
)

# Student's t distribution =====================
dist_t <- Dist$new(
  # assumes that degrees of freedom > 2 (so mean and variance are finite)
  name = "t", 
  name_long = "Student's t",
  pdf = function(x, mean, scale, log = FALSE) {
    y <- x - mean
    df <- 2 * scale^2 / (scale^2 - 1)
    return(dt(y, df = df, log = log))
  }, 
  cdf = function(q, mean, scale) {
    pt(q = q - mean, df = 2 * scale^2 / (scale^2 - 1))
  },
  rng = function(n, mean, scale) {
    df <- 2 * scale^2 / (scale^2 - 1)
    y <- rt(n, df = df)
    x <- y + mean 
    return(x)
  }, 
  link = list(mean = identity, scale = log), 
  invlink = list(mean = identity, scale = exp), 
  npar = 2, 
  parnames = c("mean", "scale"), 
  parapprox = function(x) {
    return(c(mean(x), sd(x)))
  }
)

# Log-Normal ===================================
dist_lnorm <- Dist$new(
  name = "lnorm", 
  name_long = "log-normal",
  pdf = dlnorm,
  cdf = plnorm,
  rng = rlnorm,
  link = list(meanlog = identity, sdlog = log),
  invlink = list(meanlog = identity, sdlog = exp),
  npar = 2, 
  parnames = c("meanlog", "sdlog"), 
  parapprox = function(x) {
    return(c(mean(log(x)), sd(log(x))))
  }
)

# Gamma (shape, scale) ========================================
dist_gamma <- Dist$new(
  name = "gamma",
  pdf = dgamma,
  cdf = pgamma,
  rng = rgamma,
  link = list(shape = log, scale = log),
  invlink = list(shape = exp, scale = exp),
  npar = 2, 
  parnames = c("shape", "scale"), 
  parapprox = function(x) {
    # method of moments estimators
    mean <- mean(x)
    sd <- sd(x)
    scale <- sd^2 / mean 
    shape <- mean / scale 
    return(c(shape, scale))
  }
)

# Gamma (mean, sd) ========================================
dist_gamma2 <- Dist$new(
  name = "gamma2", 
  pdf = function(x, mean, sd, log = FALSE) {
    scale <- sd^2 / mean 
    shape <- mean / scale 
    l <- dgamma(x, shape = shape, scale = scale, log = log)
    return(l)
  },
  cdf = function(q, mean, sd) {
    scale <- sd^2 / mean 
    shape <- mean / scale 
    p <- pgamma(q = q, shape = shape, scale = scale)
    return(p)
  },
  rng = function(n, mean, sd) {
    scale <- sd^2 / mean 
    shape <- mean / scale 
    return(rgamma(n, shape = shape, scale = scale))
  },
  link = list(mean = log, sd = log),
  invlink = list(mean = exp, sd = exp),
  npar = 2, 
  parnames = c("mean", "sd"), 
  parapprox = function(x) {
    # method of moments estimators
    mean <- mean(x)
    sd <- sd(x)
    return(c(mean, sd))
  }
)

# Zero-inflated gamma (shape, scale) ===================================
dist_zigamma <- Dist$new(
  name = "zigamma", 
  name_long = "zero-inflated gamma",
  pdf = function(x, shape, scale, z, log = FALSE) {
    l <- ifelse(x == 0, 
                yes = z,
                no = (1 - z) * dgamma(x, shape = shape, scale = scale))
    if (log) l <- log(l)
    return(l)
  }, 
  cdf = function(q, shape, scale, z) {
    return(NA)
  },
  rng = function(n, shape, scale, z) {
    zero <- rbinom(n, 1, z)
    y <- rgamma(n, shape = shape, scale = scale)
    y[zero == 1] <- 0
    return(y)},
  link = list(shape = log, scale = log, z = qlogis), 
  invlink = list(shape = exp, scale = exp, z = plogis), 
  npar = 3, 
  parnames = c("shape", "scale", "z"), 
  parapprox = function(x) {
    z <- length(which(x < 1e-10))/length(x)
    y <- x[which(x > 1e-10)]
    # method of moments estimators
    mean <- mean(y)
    sd <- sd(y)
    scale <- sd^2 / mean 
    shape <- mean / scale 
    return(c(shape, scale, z))
  }
)

# Zero-inflated gamma (mean, sd) ===================================
dist_zigamma2 <- Dist$new(
  name = "zigamma2", 
  name_long = "zero-inflated gamma2",
  pdf = function(x, mean, sd, z, log = FALSE) {
    shape <- mean^2 / sd^2
    scale <- sd^2 / mean
    l <- ifelse(x == 0, 
                yes = z,
                no = (1 - z) * dgamma(x, shape = shape, scale = scale))
    if (log) l <- log(l)
    return(l)
  }, 
  cdf = function(q, mean, sd, z) {
    return(NA)
  },
  rng = function(n, mean, sd, z) {
    shape <- mean^2 / sd^2
    scale <- sd^2 / mean
    zero <- rbinom(n, 1, z)
    y <- rgamma(n, shape = shape, scale = scale)
    y[zero == 1] <- 0
    return(y)},
  link = list(mean = log, sd = log, z = qlogis), 
  invlink = list(mean = exp, sd = exp, z = plogis), 
  npar = 3, 
  parnames = c("mean", "sd", "z"), 
  parapprox = function(x) {
    z <- length(which(x < 1e-10))/length(x)
    y <- x[which(x > 1e-10)]
    # method of moments estimators
    mean <- mean(y)
    sd <- sd(y)
    return(c(mean, sd, z))
  }
)

# Weibull ========================================
dist_weibull <- Dist$new(
  name = "weibull", 
  name_long = "Weibull",
  pdf = dweibull,
  cdf = pweibull,
  rng = rweibull,
  link = list(shape = log, scale = log),
  invlink = list(shape = exp, scale = exp),
  npar = 2, 
  parnames = c("shape", "scale"), 
  parapprox = function(x) {
    # regress empirical survival function against sampled values
    tmpfn <- ecdf(x)
    lx <- log(x)
    y <- -log(1 - tmpfn(x))
    lx <- lx[y < 1]
    y <- y[y < 1]
    m <- lm(y ~ lx)
    shape <- coef(m)[2]
    scale <- exp(-coef(m)[1] / shape)
    return(c(shape, scale))
  }
)

# Exponential ==================================
dist_exp <- Dist$new(
  name = "exp", 
  name_long = "exponential",
  pdf = dexp, 
  cdf = pexp,
  rng = rexp, 
  link = list(rate = log), 
  invlink = list(rate = exp), 
  npar = 1,
  parnames = c("rate"), 
  parapprox = function(x) {
    return(1 / mean(x))
  }
)

# Beta =========================================
dist_beta <- Dist$new(
  name = "beta",
  pdf = dbeta,
  cdf = pbeta,
  rng = rbeta,
  link = list(shape1 = log, shape2 = log),
  invlink = list(shape1 = exp, shape2 = exp),
  npar = 2,
  parnames = c("shape1", "shape2"), 
  parapprox = function(x) {
    # method of moments 
    mu <- mean(x)
    s2 <- var(x)
    tmp <- mu * (1 - mu) / s2 - 1
    if (tmp < 1e-10) {
      shape1 <- 1e-3
      shape2 <- 1e-3
    } else {
      shape1 <- mu * tmp
      shape2 <- (1 - mu) * tmp      
    }
    return(c(shape1, shape2))
  }
)

# Zero-one-inflated beta ========================
dist_zoibeta <- Dist$new(
  name = "zoibeta",
  pdf = function(x, shape1, shape2, zeromass, onemass, log = FALSE) {
    l <- rep(NA, length(x))
    l[which(x == 0)] <- zeromass
    l[which(x == 1)] <- onemass
    l[which(x > 0 & x < 1)] <- (1 - zeromass - onemass) * 
      dbeta(x = x, shape1 = shape1, shape2 = shape2)
    if(log) l <- log(l)
    return(l)
  },
  cdf = function(q, shape1, shape2, zeromass, onemass) {
    return(NA)
  },
  rng = function(n, shape1, shape2, zeromass, onemass) {
    y <- sample(0:2, size = n, replace = TRUE, 
                prob = c(zeromass, onemass, 1 - zeromass - onemass))
    # Indices of non-0 and non-1 observations
    ind <- which(y == 2)
    y[ind] <- rbeta(n = length(ind), shape1 = shape1, shape2 = shape2)
    return(y)
  },
  link = list(shape1 = log, shape2 = log, zeromass = qlogis, onemass = qlogis),
  invlink = list(shape1 = exp, shape2 = exp, zeromass = plogis, onemass = plogis),
  npar = 4,
  parnames = c("shape1", "shape2", "zeromass", "onemass"), 
  parapprox = function(x) {
    zeromass <- length(which(x == 0))/length(x)
    onemass <- length(which(x == 1))/length(x)
    x <- x[which(x > 0 & x < 1)]
    # method of moments 
    mu <- mean(x)
    s2 <- var(x)
    tmp <- mu * (1 - mu) / s2 - 1
    if (tmp < 1e-10) {
      shape1 <- 1e-3
      shape2 <- 1e-3
    } else {
      shape1 <- mu * tmp
      shape2 <- (1 - mu) * tmp      
    }
    return(c(shape1, shape2, zeromass, onemass))
  }
)

# Mixed distributions -----------------------------------------------------

# Tweedie ======================================
#' @importFrom mgcv ldTweedie rTweedie
dist_tweedie <- Dist$new(
  name = "tweedie", 
  name_long = "Tweedie",
  pdf = function(x, mean, p, phi, log = FALSE) {
    l <- ldTweedie(x, mu = mean, p = p + 1, phi = phi)[1,1]
    if (!log) l <- exp(l)
    return(l)
  }, 
  cdf = function(q, mean, p, phi) {
    return(NA)
  },
  rng = function(n, mean, p, phi) {
    return(rTweedie(rep(mean, n), p + 1, phi))
  }, 
  link = list(mean = identity, p = qlogis, phi = log), 
  invlink = list(mean = identity, p = plogis, phi = exp), 
  npar = 3, 
  parnames = c("mean", "p", "phi"), 
  parapprox = function(x) {
    p <- 0.5
    mean <- mean(x)
    phi <- sqrt(var(x) / mean^p)
    return(c(mean, p, phi))
  }
)

# Angular distributions ---------------------------------------------------

# Von Mises ====================================
dist_vm <- Dist$new(
  name = "vm",
  name_long = "von Mises",
  pdf = function(x, mu = 0, kappa = 1, log = FALSE) {
    dvm(x = x, mu = mu, kappa = kappa, log = log)
  },
  cdf = function(q, mu, kappa) {
    return(NA)
  },
  rng = function(n, mu, kappa) {
    rvm(n = n, mu = mu, kappa = kappa)
  },
  link = list(mu = function(x) qlogis((x + pi) / (2 * pi)),
              kappa = log),
  invlink = list(mu = function(x) 2 * pi * plogis(x) - pi,
                 kappa = exp), 
  npar = 2, 
  parnames = c("mu", "kappa"), 
  parapprox = function(x) {
    # approximate Von-Mises with Wrapped Normal
    mcosx <- mean(cos(x))
    msinx <- mean(sin(x))
    mu <- atan2(msinx, mcosx)
    r <- mcosx^2 + msinx^2 
    n <- length(x)
    r <- n * (r - 1 / n) / (n - 1)
    s2 <- log(1/r)
    kappa <- 1 - exp(-s2/2)
    return(c(mu, kappa))
  }
)

# Wrapped Cauchy ====================================
dist_wrpcauchy <- Dist$new(
  name = "wrpcauchy",
  name_long = "wrapped Cauchy",
  pdf = function(x, mu, rho, log = FALSE) {
    dwrpcauchy(x = x, mu = mu, rho = rho, log = log)
  },
  cdf = function(q, mu, rho) {
    return(NA)
  },
  rng = function(n, mu, rho) {
    rwrpcauchy(n = n, mu = mu, rho = rho)
  },
  link = list(mu = function(x) qlogis((x + pi) / (2 * pi)),
              rho = qlogis),
  invlink = list(mu = function(x) 2 * pi * plogis(x) - pi,
                 rho = plogis), 
  npar = 2, 
  parnames = c("mu", "rho"), 
  parapprox = function(x) {
    mcosx <- mean(cos(x))
    msinx <- mean(sin(x))
    mu <- atan2(msinx, mcosx)
    # Estimate of concentration taken from Wikipedia article on wrapped Cauchy
    # (they use gamma parameter where rho = e^-gamma)
    r <- mcosx^2 + msinx^2 
    n <- length(x)
    r <- n/(n-1) * (r - 1/n)
    rho <- sqrt(r)
    return(c(mu, rho))
  }
)

# Multivariate distributions ----------------------------------------------

# Multivariate Normal ========================== 
dist_mvnorm <- Dist$new(
  name = "mvnorm", 
  name_long = "multivariate normal",
  pdf = function(x, ...,  log = FALSE) {
    par <- c(...)
    # Dimension
    m <- quad_pos_solve(1, 3, - 2 * length(par))
    y <- do.call(cbind, as.matrix(x))
    # Unpack parameters
    mu <- par[1:m]
    sds <- par[(m + 1) : (2 * m)]
    corr <- par[(2 * m + 1) : (2 * m + (m^2 - m) / 2)]
    # Create covariance matrix
    V <- make_cov(sds, corr)
    p <- dmvn(y, mu, V)
    if (!log) p <- exp(p)
    return(p)
  }, 
  cdf = function(q, ...) {
    return(NA)
  },
  rng = function(n, ...) {
    par <- c(...)
    # Dimension
    m <- quad_pos_solve(1, 3, - 2 * length(par))
    # Unpack parameters
    mu <- par[1:m]
    sds <- par[(m + 1) : (2 * m)]
    corr <- par[(2 * m + 1) : (2 * m + (m^2 - m) / 2)]
    V <- make_cov(sds, corr)
    sims <- rmvn(n, mu, V)
    sims <- split(sims, 1:n)
    return(sims)
  }, 
  link = function(x, n_states) {
    xmat <- unlist(x)
    xmat <- matrix(xmat, nr = n_states, nc = length(xmat) / n_states)
    ymat <- t(apply(xmat, 1, mvnorm_link))
    return(as.vector(ymat))
  }, 
  invlink = function(x, n_states) {
    xmat <- matrix(x, nr = n_states, nc = length(x) / n_states)
    ymat <- t(apply(xmat, 1, mvnorm_invlink))
    return(as.vector(ymat))
  }, 
  npar = 8, 
  parnames = c("mu1", "mu2", "...", "sd1", "sd2", "...", "corr12", "..."), 
  parapprox = function(x) {
    y <- do.call(cbind, as.matrix(x))
    mu <- rowMeans(y)
    sds <- apply(y, 1, sd)
    corr <- cor(t(y))
    return(c(mu, sds, corr[upper.tri(corr)]))
  },
  par_alt = function(par) {
    # Dimension
    m <- quad_pos_solve(1, 3, - 2 * length(par))
    # Unpack parameters
    mu <- par[1:m]
    sds <- par[(m + 1) : (2 * m)]
    corr <- par[(2 * m + 1) : (2 * m + (m^2 - m) / 2)]
    # Create covariance matrix
    V <- make_cov(sds, corr)
    names(mu) <- paste0("mu", 1:m)
    rownames(V) <- colnames(V) <- 1:m
    return(list(mu = mu, Sigma = V))
  }
)

# Dirichlet Distribution =======================
dist_dir <- Dist$new(
  name = "dir", 
  name_long = "Dirichlet",
  pdf = function(x, ...,  log = FALSE) {
    alpha <- c(...)
    y <- do.call(cbind, as.matrix(x))
    p <- gamma(sum(alpha)) * prod(y ^ (alpha - 1)) /  prod(gamma(alpha))
    if (log) p <- log(p)
    return(p)
  }, 
  cdf = function(q, ...) {
    return(NA)
  },
  rng = function(n, ...) {
    alpha <- c(...)
    y <- rgamma(n * length(alpha), shape = alpha, scale = 1)
    x <- matrix(y, nr = length(alpha), nc = n)
    x <- apply(x, 2, FUN = function(r) {r / sum(r)})
    x <- split(x, rep(1:ncol(x), each = nrow(x)))
    return(x)
  }, 
  link = function(x, n_states) {log(unlist(x))}, 
  invlink = function(x, n_states) {exp(x)}, 
  npar = 2, 
  parnames = c("alpha1", "alpha2"), 
  parapprox = function(x) {
    y <- do.call(cbind, as.matrix(x))
    mu <- rowMeans(y)
    sds <- apply(y, 1, sd)
    a0 <- mean(mu * (1 - mu) / sds^2 - 1)
    alpha <- mu * a0 
    return(alpha)
  }
)

# List of all distributions -----------------------------------------------

dist_list <- list(beta = dist_beta,
                  binom = dist_binom,
                  cat = dist_cat,
                  dir = dist_dir,
                  exp = dist_exp,
                  foldednorm = dist_foldednorm, 
                  gamma = dist_gamma,
                  gamma2 = dist_gamma2,
                  lnorm = dist_lnorm,
                  mvnorm = dist_mvnorm, 
                  nbinom = dist_nbinom,
                  nbinom2 = dist_nbinom2,
                  norm = dist_norm,
                  pois = dist_pois,
                  t = dist_t, 
                  truncnorm = dist_truncnorm, 
                  tweedie = dist_tweedie, 
                  vm = dist_vm, 
                  weibull = dist_weibull,
                  wrpcauchy = dist_wrpcauchy,
                  zibinom = dist_zibinom, 
                  zigamma = dist_zigamma,
                  zigamma2 = dist_zigamma2,
                  zinbinom = dist_zinbinom,
                  zipois = dist_zipois,
                  zoibeta = dist_zoibeta,
                  ztnbinom = dist_ztnbinom, 
                  ztpois = dist_ztpois)

# Define unique distribution code (must match C++ side)
lapply(seq_along(dist_list), function(i) {
  dist_list[[i]]$set_code(new_code = i - 1)
})
