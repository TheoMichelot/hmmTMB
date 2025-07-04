
#' hmmTMB colour palette
hmmTMB_cols <- c("#00798c", "#d1495b", "#edae49", "#66a182", "#2e4057", "#8d96a3")

#' Create block diagonal matrix (safe version)
#' 
#' This version of bdiag checks whether the matrices passed as
#' arguments are NULL. This avoids errors that would arise if
#' using bdiag directly.
#' 
#' @param ... Matrix or list of matrices (only the first argument is used)
#' 
#' @return Block diagonal matrix
#' 
#' @importFrom Matrix bdiag
bdiag_check <- function(...) {
  # Only use first argument (matrix or list of matrices)
  args <- list(...)[[1]]
  
  # Check which matrices are non-NULL
  # (Needs two conditions to keep non-empty vectors which don't have
  # dimensions, and also keep empty matrices with set dimensions)
  check <- sapply(args, function(arg) {
    !is.null(dim(arg)) | length(arg) > 0
  })
  
  if(length(check) == 0)
    return(NULL)
  else
    return(bdiag(args[check]))
}

#' Generalized matrix determinant
#' 
#' Generalized determinant = product of non-zero eigenvalues 
#' (see e.g., Wood 2017). Used for (log)determinant of penalty matrices,
#' required in log-likelihood function. 
#' 
#' @param x Numeric matrix
#' @param eps Threshold below which eigenvalues are ignored (default: 1e-10)
#' @param log Logical: should the log-determinant be returned?
#' 
#' @return Generalized determinant of input matrix
gdeterminant <- function(x, eps = 1e-10, log = TRUE) {
  if(is.null(x)) {
    return(NULL)
  } else {
    # Compute sum of log of non-zero eigenvalues
    # (i.e., log generalized determinant)
    eigenpairs <- eigen(x)
    eigenvalues <- eigenpairs$values
    logdet <- sum(log(eigenvalues[eigenvalues > eps]))
    if(!log) {
      return(exp(logdet))
    } else{
      return(logdet)
    }        
  }
}

#' Fill in NAs
#' 
#' Replace NA entries in a vector by the last non-NA value. If the first
#' entry of the vector is NA, it is replaced by the first non-NA value. If the
#' vector passed as input doesn't contain NAs, it is returned as is.
#' 
#' @param x Vector in which NAs should be removed
#' 
#' @return Copy of x in which NAs have been replaced by nearest available value.
na_fill <- function(x) {
  # If no missing value, return x as is
  if(!any(is.na(x))) {
    return(x)
  }
  
  message("Replacing NAs in covariates by last non-NA values.")
  
  # Copy x
  y <- x
  # Replace first value by first non-NA value, just in case is.na(y[1])
  y[1] <- x[which(!is.na(x))[1]]
  
  # Loop over NA entries, and replace by previous entry
  na_ind <- which(is.na(y))
  for(i in na_ind) {
    y[i] <- y[i-1]
  }
  
  return(y)
}

#' Grid of covariates
#' 
#' @param var Name of variable
#' @param data Data frame containing the covariates. If not provided, data
#' are extracted from obj
#' @param obj HMM model object containing data and formulas
#' @param covs Optional named list for values of covariates (other than 'var') 
#' that should be used in the plot (or dataframe with single row). If this is
#' not specified, the mean value is used for numeric variables, and the
#' first level for factor variables.
#' @param formulas List of formulas used in the model
#' @param n_grid Grid size (number of points). Default: 1000.
#' 
#' @return Data frame of covariates, with 'var' defined over a grid,
#' and other covariates fixed to their mean (numeric) or first level
#' (factor).
cov_grid <- function(var, data = NULL, obj = NULL, covs = NULL, formulas, n_grid = 1e3) {
  # Get data set
  if(is.null(data)) {
    data <- obj$obs()$data()    
  }
  
  # Get covariate names
  if(!is.null(obj)) {
    var_names <- unique(c(rapply(obj$obs()$formulas(), all.vars), 
                          rapply(obj$hid()$formulas(), all.vars)))
  } else {
    var_names <- unique(rapply(formulas, all.vars))
  }
  # If no covariates in the model, only take covariate 'var'
  if(length(var_names) == 0)
    var_names <- var
  
  # pi might appear in the formulas (e.g. used in periodic terms),
  # in which case it needs to be added to the data frame
  if(any(var_names == "pi")) {
    data$pi <- pi
  }
  
  # Get data frame of covariates
  all_vars <- data[, var_names, drop = FALSE]
  
  # Grid of covariate
  if(is.factor(all_vars[, var])) {
    n_grid <- length(unique(all_vars[, var]))
    grid <- unique(all_vars[, var])
  } else {
    grid <- seq(min(all_vars[, var]), max(all_vars[, var]), length = n_grid)
  }
  
  # New data frame for covariate grid
  new_data <- matrix(NA, nrow = n_grid, ncol = ncol(all_vars))
  colnames(new_data) <- colnames(all_vars)
  new_data <- as.data.frame(new_data)
  new_data[, var] <- grid
  
  # Select value for other covariates
  covs_list <- lapply(seq_along(all_vars), function(i) {
    if(!is.null(covs[[var_names[i]]])) {
      # Set to user-provided value if possible
      return(covs[[var_names[i]]])
    } else {
      # No user-provided value
      col <- all_vars[, i]
      if(is.numeric(col)) {
        # If numeric, use mean value
        return(mean(col, na.rm = TRUE)) 
      } else {
        # If factor, use first factor level
        return(unique(col)[1])
      }
    }
  })
  covs <- as.data.frame(covs_list)
  colnames(covs) <- colnames(all_vars)
  
  # Fill columns for other covariates
  for(var_name in colnames(new_data)) {
    if(var_name != var)
      new_data[, var_name] <- covs[, var_name]
  }
  
  return(new_data)
}

#' Check if number of whole number 
#'
#' @param x number to check or vector of numbers 
#' @param tol how far away from whole number is ok? 
#'
#' @return TRUE if it is a whole number within tolerance 
#' @export
is_whole_number <- function(x, tol = 1e-10) {
  y <- round(x)
  return(all(abs(x - y) < tol, na.rm = TRUE))
}

#' Strip comments marked with a hash from a character vector 
#' 
#' @param str the character vector 
#' 
#' @return character vector with comments removed (and lines with only comments
#' completely removed) 
#' @export 
#' 
strip_comments <- function(str) {
  res <- str_trim(str_split_fixed(str, "#", 2)[, 1])
  res <- res[res != ""]
  return(res)
}

#' Solve for positive root of quadratic ax^2 + bx + c = 0 when it exists 
#' @param a coefficient of x^2
#' @param b coefficient of x 
#' @param c scalar coefficient 
#' @return real positive root if it exists 
#' @export
quad_pos_solve <- function(a, b, c) {
  deter <- b^2 - 4 * a * c
  if (deter < 0) stop("quadratic has no real roots")
  roots <- (-b + c(-1, 1) * sqrt(deter)) / (2 * a)
  if (all(roots < 1e-10)) stop("no positive roots for this quadratic")
  roots <- roots[roots > 1e-10]
  return(as.numeric(roots))
}

#' Multivariate logit function 
#'
#' @param x Numeric vector
#' 
#' @export 
mlogit <- function(x) {
  s <- 1 - sum(x)
  return(log(x / s))
}

#' Multivarite inverse logit function
#' 
#' @param x Numeric vector
#' 
#' @export 
invmlogit <- function(x) {
  y <- exp(x)
  s <- 1/(1 + sum(y))
  y <- y * s
  return(y)
}

#' Make covariance matrix from standard deviations and correlations
#' 
#' @param sds Vector of standard deviations
#' @param corr Vector of correlations (must be of length m*(m-1)/2 if
#' sds is of length m)
#' 
#' @return An m by m covariance matrix
make_cov <- function(sds, corr) {
  m <- length(sds)
  V <- diag(m)
  V[lower.tri(V)] <- corr 
  V[upper.tri(V)] <- t(V)[upper.tri(V)]
  for (i in 1:ncol(V)) {
    V[i,] <- V[i,] * sds[i]
    V[,i] <- V[,i] * sds[i]
  }
  return(V)
}

#' Multivariate Normal link function 
#'
#' @param x Vector of parameters on natural scale (in the order: means,
#' SDs, correlations)
#'
#' @export 
#' 
#' @importFrom stats qlogis
mvnorm_link <- function(x) {
  # Get dimension 
  m <- quad_pos_solve(1, 3, - 2 * length(x))
  
  # Mean parameters (untransformed)
  mu <- x[1:m]
  
  # Standard deviations and correlations
  sds <- x[(m + 1) : (2 * m)]
  corr <- x[(2 * m + 1) : (2 * m + (m^2 - m) / 2)]
  
  # Cholesky decomposition of covariance matrix
  cov_mat <- make_cov(sds, corr)
  cov_chol <- t(chol(cov_mat))
  
  # Diagonal elements of Cholesky factor are positive, so log
  # transform them to get unconstrained parameters
  # See: https://mc-stan.org/docs/reference-manual/transforms.html#cholesky-factors-of-covariance-matrices
  cov_par_diag <- log(diag(cov_chol))
  cov_par_offdiag <- cov_chol[lower.tri(cov_chol)]
  
  return(c(mu, cov_par_diag, cov_par_offdiag))
}

#' Multivariate Normal inverse link function 
#' 
#' @param x Vector of parameters on linear predictor scale (in the order:
#' means, SDs, correlations)
#' 
#' @export 
#'
#' @importFrom stats plogis
mvnorm_invlink = function(x) {
  # Get dimension 
  m <- quad_pos_solve(1, 3, - 2 * length(x))
  
  # Mean parameters (untransformed)
  mu <- x[1:m]
  
  # Unpack parameters of Cholesky factor of covariance matrix
  # Diagonal elements of factor must be positive so we 
  # exponentiate them here
  # See: https://mc-stan.org/docs/reference-manual/transforms.html#cholesky-factors-of-covariance-matrices
  cov_par_diag <- exp(x[(m+1):(2*m)])
  cov_par_offdiag <- x[(2*m+1):(2*m+(m^2-m)/2)]
  cov_chol <- diag(cov_par_diag)
  cov_chol[lower.tri(cov_chol)] <- cov_par_offdiag
  
  # Get covariance matrix from Cholesky factor
  cov_mat <- cov_chol %*% t(cov_chol)
  
  # Get standard deviations and correlations from covariance matrix
  sds <- sqrt(diag(cov_mat))
  corr_mat <- cov_mat / rep(sds, each = m) / rep(sds, m)
  corr <- corr_mat[lower.tri(corr_mat)]
  
  return(c(mu, sds, corr))
}

#' Log of sum of exponentials 
#' 
#' @param x Numeric vector
#' 
#' @export 
logsumexp <- function(x) {
  xmax <- max(x)
  val <- xmax + log(sum(exp(x - xmax)))
  return(val)
}

#' Read formula with as.character without splitting
#' 
#' @param x R formula
#' @param ... Unused
#' 
#' @details Citation: this function was taken from the R package
#' formula.tools: 
#'   Christopher Brown (2018). formula.tools: Programmatic Utilities for Manipulating Formulas,
#'   Expressions, Calls, Assignments and Other R Objects. R package version 1.7.1.
#'   https://CRAN.R-project.org/package=formula.tools
#' @export 
as_character_formula <- function (x, ...) 
{
  form <- paste(deparse(x), collapse = " ")
  form <- gsub("\\s+", " ", form, perl = FALSE)
  return(form)
}

#' Get covariance matrix from precision matrix
#' 
#' The covariance matrix is the inverse of the precision matrix. By default,
#' the function \code{solve} is used for inversion. If it fails (e.g.,
#' singular system), then \code{MASS::ginv} is used instead, and returns the
#' Moore-Penrose generalised inverse of the precision matrix.
#' 
#' @param prec_mat Precision matrix (either of 'matrix' type
#' or sparse matrix on which as.matrix can be used)
#' 
#' @return Precision matrix
#' 
#' @importFrom MASS ginv
#' @export
prec_to_cov <- function(prec_mat)
{
  cov_mat <- try(as.matrix(solve(prec_mat)), silent = TRUE)
  if(inherits(cov_mat, "try-error")) {
    message <- attr(cov_mat, 'condition')$message
    cov_mat <- ginv(as.matrix(prec_mat))
    warning(paste0("Inversion of precision matrix using 'solve' failed: ", 
                   message, ". Using 'MASS::ginv' instead (uncertainty ",
                   "estimates may be unreliable)."))
  }
  return(cov_mat)
}

#' Find s(, bs = "re") terms in formula
#' 
#' This function is used to identify the variables "x" which are 
#' included as s(x, bs = "re") in the formula, in particular to
#' check that they are factors.
#' 
#' @param form Model formula
#' 
#' @return Vector of names of variables for which a random
#' effect term is included in the model.
#'
#' @importFrom stats terms
find_re <- function(form) {
  term_labs <- attr(terms(form), "term.labels")
  # Regex description:
  # ^: start of string
  # s\\(: s followed by opening bracket
  # (.*): this part can be captured by gsub
  # , bs = "re"\\): the rest of the string
  # $: end of string
  pattern <- '^s\\((.*), bs = "re"\\)$'
  var_names <- gsub(pattern = pattern, "\\1", term_labs)
  which_re <- grep(pattern = pattern, term_labs)
  var_re <- var_names[which_re]
  if(length(var_re) > 0) {
    # This is needed for models with "by", e.g. ~s(ID, x, bs = "re"),
    # to make sure only "ID" is kept, rather than "ID, x"
    var_re <- strsplit(var_re, split = ",")[[1]][1]
  }
  return(var_re)
}

#' Transforms matrix to dgTMatrix
#' 
#' @param x Matrix or vector. If this is a vector, it is formatted into
#' a single-column matrix.
#' 
#' @return Sparse matrix of class dgTMatrix
#' 
#' @importFrom methods as
as_sparse <- function(x) {
  if(length(dim(x)) < 2) {
    x <- matrix(x, ncol = 1)
  }
  # # This is the syntax recommended by Matrix > 1.5.0, but doesn't seem
  # # to be compatible with earlier versions of Matrix.
  # mat <- as(as(as(x, "dMatrix"), "generalMatrix"), "TsparseMatrix")
  mat <- suppressMessages(as(x, "dgTMatrix"))
  return(mat)
}

#' Check values in vector are contiguous
#' 
#' @param x Vector of values (can be numeric, character, factor)
#' 
#' @return Logical: are values contiguous?
check_contiguous <- function(x) {
  vals <- unique(x)
  is_contiguous <- sapply(vals, function(val) {
    ind <- which(x == val)
    (ind[length(ind)] - ind[1]) == (length(ind) - 1)
  })
  return(all(is_contiguous))
}

#' Density function of von Mises distribution
#' 
#' @param x Angle
#' @param mu Mean parameter
#' @param kappa Concentration parameter
#' @param log Should log-density be returned?
#' 
#' @return Von Mises density
dvm <- function(x, mu, kappa, log = FALSE) {
  # The "- kappa" term below cancels out the expon.scaled
  b <- besselI(kappa, 0, expon.scaled = TRUE)
  val <- - log(2 * pi * b) + kappa * cos(x - mu) - kappa
  if(!log) {
    val <- exp(val)        
  }
  return(val)
}

#' Sample from von Mises distribution
#' 
#' @param n Number of samples
#' @param mu Mean parameter
#' @param kappa Concentration parameter
#' 
#' @return Vector of n samples from vm(mu, kappa)
#' 
#' @details Uses basic rejection sampling, based on dvm(), which might
#' be inefficient for large kappa. Could be improved following Best & Fisher 
#' (1979), Efficient simulation of the von Mises distribution, JRSSC, 28(2), 
#' 152-157.
#' 
#' @importFrom stats runif
rvm <- function(n, mu, kappa) {
  x <- rep(NA, n)
  n_accept <- 0
  pdf_max <- dvm(x = mu, mu = mu, kappa = kappa)
  while(n_accept < n) {
    x_star <- runif(1, min = -pi, max = pi)
    pdf_star <- dvm(x = x_star, mu = mu, kappa = kappa)
    accept_prob <- pdf_star/pdf_max
    if(runif(1) < accept_prob) {
      n_accept <- n_accept + 1
      x[n_accept] <- x_star
    }
  }
  return(x)
}

#' Density function of wrapped Cauchy distribution
#' 
#' @param x Angle
#' @param mu Mean parameter
#' @param rho Concentration parameter
#' @param log Should log-density be returned?
#' 
#' @return Wrapped Cauchy density
dwrpcauchy <- function(x, mu, rho, log = FALSE) {
  val <- (1 - rho^2) / 
    (2 * pi * (1 + rho^2 - 2 * rho * cos(x - mu)))
  if(log) {
    val <- log(val)
  }
  return(val)
}

#' Sample from wrapped Cauchy distribution
#' 
#' @param n Number of samples
#' @param mu Mean parameter
#' @param rho Concentration parameter
#' 
#' @return Vector of n samples from wrpcauchy(mu, rho)
#' 
#' @details Uses basic rejection sampling, based on dwrpcauchy(), which might
#' be inefficient for large rho.
rwrpcauchy <- function(n, mu, rho) {
  x <- rep(NA, n)
  n_accept <- 0
  pdf_max <- dwrpcauchy(x = mu, mu = mu, rho = rho)
  while(n_accept < n) {
    x_star <- runif(1, min = -pi, max = pi)
    pdf_star <- dwrpcauchy(x = x_star, mu = mu, rho = rho)
    accept_prob <- pdf_star/pdf_max
    if(runif(1) < accept_prob) {
      n_accept <- n_accept + 1
      x[n_accept] <- x_star
    }
  }
  return(x)
}
