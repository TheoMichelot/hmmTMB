
#' Create model matrices
#' 
#' @param formulas List of formulas (possibly nested, e.g. for use within Observation)
#' @param data Data frame including covariates
#' @param new_data Optional new data set, including covariates for which
#' the design matrices should be created. This needs to be passed in addition
#' to the argument '\code{data}', for cases where smooth terms or factor
#' covariates are included, and the original data set is needed to determine
#' the full range of covariate values.
#' 
#' @return A list of
#' \itemize{
#'   \item X_fe Design matrix for fixed effects
#'   \item X_re Design matrix for random effects
#'   \item S Smoothness matrix
#'   \item log_det_S Vector of log-determinants of smoothness matrices
#'   \item ncol_fe Number of columns of X_fe for each parameter
#'   \item ncol_re Number of columns of X_re and S for each random effect
#' }
#' 
#' @importFrom stats update predict
make_matrices = function(formulas, data, new_data = NULL) {
  # Initialise lists of matrices
  X_list_fe <- list()
  X_list_re <- list()
  S_list <- list()
  ncol_fe <- NULL
  ncol_re <- NULL
  names_fe <- NULL
  names_re <- NULL
  names_ncol_re <- NULL
  start <- 1
  
  # Unlist formulas so that this function works both for Observation and MarkovChain
  forms <- unlist(formulas)
  names <- names(forms)
  
  # Loop over formulas
  for(k in seq_along(forms)) {
    form <- forms[[k]]
    
    # Check that random effect variables are factors
    var_re <- find_re(form)
    for(var in var_re) {
      if(!inherits(data[[var]], "factor")) {
        data[[var]] <- factor(data[[var]])
        warning(paste0("'", var, "' is included as a random effect but is ",
                       "not a factor - changing to factor."))
      }
    }
    
    # Create matrices based on this formula
    if(is.null(new_data)) {
      gam_setup <- gam(formula = update(form, dummy_response ~ .), 
                       data = cbind(dummy_response = 1, data), 
                       fit = FALSE)
      Xmat <- gam_setup$X
      # Extract column names for design matrices
      term_names <- gam_setup$term.names
    } else {
      # Get design matrix for new data set
      gam_setup0 <- gam(formula = update(form, dummy_response ~ .), 
                       data = cbind(dummy_response = 1, data))
      gam_setup <- gam(formula = update(form, dummy_response ~ .), 
                        data = cbind(dummy_response = 1, data),
                       fit = FALSE)
      Xmat <- predict(gam_setup0, newdata = new_data, type = "lpmatrix")
      # Extract column names for design matrices
      term_names <- gam_setup$term.names
    }
    
    # Fixed effects design matrix
    X_list_fe[[k]] <- Xmat[, 1:gam_setup$nsdf, drop = FALSE]
    subnames_fe <- paste0(names[k], ".", term_names[1:gam_setup$nsdf])
    names_fe <- c(names_fe, subnames_fe)
    
    # Random effects design matrix
    X_list_re[[k]] <- Xmat[, -(1:gam_setup$nsdf), drop = FALSE]
    if(ncol(X_list_re[[k]]) > 0) {
      subnames_re <- paste0(names[k], ".", term_names[-(1:gam_setup$nsdf)])
      names_re <- c(names_re, subnames_re)                    
    }
    
    # Smoothing matrix
    S_list[[k]] <- bdiag_check(gam_setup$S)
    
    # Number of columns for fixed effects
    ncol_fe <- c(ncol_fe, gam_setup$nsdf)
    
    if(length(gam_setup$smooth) > 0) {
      sub_ncol_re <- matrix(1, nrow = 2, ncol = length(gam_setup$S))
      colnames(sub_ncol_re) <- 1:ncol(sub_ncol_re)
      start_s <- 1
      for (s in 1:length(gam_setup$smooth)) {
        # how many penalties for this smooth?
        npen <- length(gam_setup$smooth[[s]]$S)
        # how many parameters for this smooth? 
        npar <- ncol(gam_setup$smooth[[s]]$S[[1]])
        # where does this smooth's parameters start and end?
        sub_ncol_re[, (start_s:(start_s + npen - 1))] <- c(start, start + npar - 1)
        colnames(sub_ncol_re)[start_s:(start_s + npen - 1)] <- rep(gam_setup$smooth[[s]]$label, npen)
        # get names of smooth terms
        # regex from datascience.stackexchange.com/questions/8922
        s_terms <- gsub("(.*)\\..*", "\\1", names_re[sub_ncol_re[1, s]:sub_ncol_re[2, s]])
        names_ncol_re <- c(names_ncol_re, rep(unique(s_terms), npen))
        start <- start + npar
        start_s <- start_s + npen
      }
      ncol_re <- cbind(ncol_re, sub_ncol_re)
    }    
  }
  colnames(ncol_re) <- names_ncol_re
  
  # Store as block diagonal matrices
  X_fe <- bdiag_check(X_list_fe)
  colnames(X_fe) <- names_fe
  X_re <- bdiag_check(X_list_re)
  colnames(X_re) <- names_re
  S <- bdiag_check(S_list)
  
  # Get (log-)determinants of penalty matrices
  log_det_S <- unlist(sapply(S_list, gdeterminant))
  
  return(list(X_fe = X_fe, 
              X_re = X_re, 
              S = S,
              log_det_S = log_det_S,
              X_list_fe = X_list_fe, 
              X_list_re = X_list_re, 
              S_list = S_list, 
              ncol_fe = ncol_fe, 
              ncol_re = ncol_re))
}
