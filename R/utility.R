
#' Grid of covariates
#' 
#' @param var Name of variable
#' @param data Data frame containing the covariates
#' @param covs Optional data frame with a single row and one column
#' for each covariate, giving the values that should be used. If this is
#' not specified, the mean value is used for numeric variables, and the
#' first level for factor variables.
#' @param formulas List of formulas used in the model
#' 
#' @return Data frame of covariates, with 'var' defined over a grid,
#' and other covariates fixed to their mean (numeric) or first level
#' (factor).
cov_grid <- function(var, data, covs = NULL, formulas) {
  # Get covariate names
  var_names <- unique(rapply(formulas, all.vars))
  
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
    n_grid <- 1e3
    grid <- seq(min(all_vars[, var]), max(all_vars[, var]), length = n_grid)
  }
  
  # New data frame for covariate grid
  new_data <- matrix(NA, nrow = n_grid, ncol = ncol(all_vars))
  colnames(new_data) <- colnames(all_vars)
  new_data <- as.data.frame(new_data)
  new_data[, var] <- grid
  
  # Select value for other covariates
  if(is.null(covs)) {
    covs_list <- lapply(all_vars, function(col) {
      if(is.numeric(col)) {
        # If numeric, use mean value
        return(mean(col, na.rm = TRUE)) 
      } else {
        # If factor, use first factor level
        return(unique(col)[1])
      }
    })
    covs <- as.data.frame(covs_list)
  }
  # Fill columns for other covariates
  for(var_name in colnames(new_data)) {
    if(var_name != var)
      new_data[, var_name] <- covs[, var_name]
  }
  
  return(new_data)
}
