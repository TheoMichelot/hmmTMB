
#' Create model matrices (observation model)
#' 
#' @param formulas Nested list of formulas, with three levels: (1) variable,
#' (2) parameter, and (3) state.
#' @param data Data frame including covariates
#' 
#' @return A list of
#' \itemize{
#'   \item X_fe Design matrix for fixed effects
#'   \item X_re Design matrix for random effects
#'   \item S Smoothness matrix
#'   \item ncol_re Number of columns of X_re and S for each random effect
#' }
make_mat_obs = function(formulas, data, new_data = NULL) {
  # Initialise lists of matrices
  X_list_fe <- list()
  X_list_re <- list()
  S_list <- list()
  ncol_fe <- NULL
  ncol_re <- NULL
  k <- 1
  
  # Loop over variables
  for(var_forms in formulas) {
    
    # Loop over parameters
    for(par_forms in var_forms) {
      
      # Loop over states
      for(form in par_forms) {
        # Create matrices based on this formula
        if(is.null(new_data)) {
          gam_setup <- gam(formula = update(form, dummy ~ .), 
                           data = cbind(dummy = 1, data), 
                           fit = FALSE)
          Xmat <- gam_setup$X
        } else {
          # Get design matrix for new data set
          gam_setup <- gam(formula = update(form, dummy ~ .), 
                           data = cbind(dummy = 1, data))
          Xmat <- predict(gam_setup, newdata = new_data, type = "lpmatrix")
        }
        
        # Fixed effects design matrix
        X_list_fe[[k]] <- Xmat[, 1:gam_setup$nsdf, drop = FALSE]
        
        # Random effects design matrix
        X_list_re[[k]] <- Xmat[, -(1:gam_setup$nsdf), drop = FALSE]
        
        # Smoothing matrix
        S_list[[k]] <- bdiag_check(gam_setup$S)
        
        # Number of columns for fixed effects
        ncol_fe <- c(ncol_fe, gam_setup$nsdf)
        
        # Number of columns for each random effect
        if(length(gam_setup$S) > 0)
          ncol_re <- c(ncol_re, sapply(gam_setup$S, ncol))
        
        k <- k + 1
      }
    }
  }
  
  # Store as block diagonal matrices
  X_fe <- bdiag_check(X_list_fe)
  X_re <- bdiag_check(X_list_re)
  S <- bdiag_check(S_list)
  
  return(list(X_fe = X_fe, X_re = X_re, S = S, 
              ncol_fe = ncol_fe, ncol_re = ncol_re))
}

#' Create model matrices (transition probabilities)
#' 
#' @param formulas List of formulas
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
#'   \item ncol_re Number of columns of X_re and S for each random effect
#' }
make_mat_hid = function(formulas, data, new_data = NULL) {
  # Initialise lists of matrices
  X_list_fe <- list()
  X_list_re <- list()
  S_list <- list()
  ncol_fe <- NULL
  ncol_re <- NULL
  k <- 1
  
  # Loop over formulas
  for(form in formulas) {
    # Create matrices based on this formula
    if(is.null(new_data)) {
      gam_setup <- gam(formula = update(form, dummy ~ .), 
                       data = cbind(dummy = 1, data), 
                       fit = FALSE)
      Xmat <- gam_setup$X
    } else {
      # Get design matrix for new data set
      gam_setup <- gam(formula = update(form, dummy ~ .), 
                       data = cbind(dummy = 1, data))
      Xmat <- predict(gam_setup, newdata = new_data, type = "lpmatrix")
    }

    # Fixed effects design matrix
    X_list_fe[[k]] <- Xmat[, 1:gam_setup$nsdf, drop = FALSE]
    
    # Random effects design matrix
    X_list_re[[k]] <- Xmat[, -(1:gam_setup$nsdf), drop = FALSE]
    
    # Smoothing matrix
    S_list[[k]] <- bdiag_check(gam_setup$S)
    
    # Number of columns for fixed effects
    ncol_fe <- c(ncol_fe, gam_setup$nsdf)
    
    # Number of columns for each random effect
    if(length(gam_setup$S) > 0)
      ncol_re <- c(ncol_re, sapply(gam_setup$S, ncol))
    
    k <- k + 1
  }
  
  # Store as block diagonal matrices
  X_fe <- bdiag_check(X_list_fe)
  X_re <- bdiag_check(X_list_re)
  S <- bdiag_check(S_list)
  
  return(list(X_fe = X_fe, X_re = X_re, S = S, 
              ncol_fe = ncol_fe, ncol_re = ncol_re))
}
