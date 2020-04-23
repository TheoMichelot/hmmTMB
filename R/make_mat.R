
#' Create model matrices
#' 
#' Used in \code{Observation} and \code{Hidden} classes, to construct design
#' and smoothing matrices for observation model and hidden process model,
#' respectively.
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
make_mat = function(formulas, data) {
  # Initialise lists of matrices
  X_list_fe <- list()
  X_list_re <- list()
  S_list <- list()
  ncol_re <- NULL
  k <- 1
  
  # Loop over variables
  for(var_forms in formulas) {
    
    # Loop over parameters
    for(par_forms in var_forms) {
      
      # Loop over states
      for(form in par_forms) {
        
        # Create matrices based on this formula
        gam_setup <- gam(formula = update(form, dummy ~ .), 
                         data = cbind(dummy = 1, data), 
                         fit = FALSE)
        
        # Fixed effects design matrix
        X_list_fe[[k]] <- gam_setup$X[, 1:gam_setup$nsdf, drop = FALSE]
        
        # Random effects design matrix
        X_list_re[[k]] <- gam_setup$X[, -(1:gam_setup$nsdf), drop = FALSE]
        
        # Smoothing matrix
        S_list[[k]] <- bdiag_check(gam_setup$S)
        
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
  
  return(list(X_fe = X_fe, X_re = X_re, S = S, ncol_re = ncol_re))
}
