
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
#' @export 
mlogit <- function(x) {
  s <- 1 - sum(x)
  return(log(x / s))
}

#' Multivarite inverse logit function
#' @export 
invmlogit <- function(x) {
  y <- exp(x)
  s <- 1/(1 + sum(y))
  y <- y * s
  return(y)
}

#' Multivariate Normal link function 
#' @export 
mvnorm_link <- function(x) {
  # get dimension 
  m <- quad_pos_solve(1, 3, - 2 * length(x))
  mu <- x[1:m]
  sds <- log(x[(m + 1) : (2 * m)])
  corr <- qlogis((x[(2 * m + 1) : (2 * m + (m^2 - m) / 2)] + 1) / 2)
  return(c(mu, sds, corr))
}

#' Multivariate Normal inverse link function 
#' @export 
mvnorm_invlink = function(x) {
  # get dimension 
  m <- quad_pos_solve(1, 3, - 2 * length(x))
  mu <- x[1:m]
  sds <- exp(x[(m + 1) : (2 * m)])
  corr <- 2 * plogis(x[(2 * m + 1) : (2 * m + (m^2 - m) / 2)]) - 1
  return(c(mu, sds, corr))
}

#' Log of sum of exponentials 
#' @export 
logsumexp <- function(x) {
  xmax <- max(x)
  val <- xmax + log(sum(exp(x - xmax)))
  return(val)
}

#' Read formula with as.character without splitting
#' Citation: this function was taken from the R package
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
#' @export
prec_to_cov <- function(prec_mat)
{
  cov_mat <- try(as.matrix(solve(prec_mat)), silent = TRUE)
  if(inherits(cov_mat, "try-error")) {
    message <- attr(cov_mat, 'condition')$message
    cov_mat <- MASS::ginv(as.matrix(prec_mat))
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
  return(var_re)
}
