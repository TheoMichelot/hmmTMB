
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
#'   \item ncol_fe Number of columns of X_fe for each parameter
#'   \item ncol_re Number of columns of X_re and S for each random effect
#' }
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
  
  # Unlist formulas so that this function works both for Observation and MarkovChain
  forms <- unlist(formulas)
  names <- names(forms)
  
  # Loop over formulas
  for(k in seq_along(forms)) {
    form <- forms[[k]]
    
    # Create matrices based on this formula
    if(is.null(new_data)) {
      gam_setup <- gam(formula = update(form, dummy ~ .), 
                       data = cbind(dummy = 1, data), 
                       fit = FALSE)
      Xmat <- gam_setup$X
      # Extract column names for design matrices
      term_names <- gam_setup$term.names
    } else {
      # Get design matrix for new data set
      gam_setup <- gam(formula = update(form, dummy ~ .), 
                       data = cbind(dummy = 1, data))
      Xmat <- predict(gam_setup, newdata = new_data, type = "lpmatrix")
      # Extract column names for design matrices
      term_names <- names(gam_setup$coefficients)
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
    
    if(length(gam_setup$S) > 0) {
      # Number of columns for each random effect
      sub_ncol_re <- sapply(gam_setup$S, ncol)
      ncol_re <- c(ncol_re, sub_ncol_re)
      # Hacky way to get the names of smooth terms 
      # (one for each column of ncol_re)
      # regex from datascience.stackexchange.com/questions/8922
      s_terms_i1 <- cumsum(c(1, sub_ncol_re[-length(sub_ncol_re)]))
      s_terms <- gsub("(.*)\\..*", "\\1", subnames_re[s_terms_i1])
      names_ncol_re <- c(names_ncol_re, s_terms)
    }    
  }
  
  # Store as block diagonal matrices
  X_fe <- bdiag_check(X_list_fe)
  colnames(X_fe) <- names_fe
  X_re <- bdiag_check(X_list_re)
  colnames(X_re) <- names_re
  S <- bdiag_check(S_list)
  
  # Name elements of ncol_re
  names(ncol_re) <- names_ncol_re
  
  
  return(list(X_fe = X_fe, X_re = X_re, S = S, 
              ncol_fe = ncol_fe, ncol_re = ncol_re))
}
