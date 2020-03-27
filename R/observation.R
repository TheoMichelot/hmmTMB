
#' Hidden Markov observation class
#'
#' @description Encapsulates the observation model from a hidden Markov model.
#' Object can be created using $new with arguments:
#' \itemize{
#'   \item data: a HmmData object
#'   \item dists: named list of distributions for each data stream
#'   \item par: list of observation parameters (for covariate-free model)
#'   \item wpar: vector of observation parameters on working scale (for
#'   model with covariates)
#'   \item formulas: list of formulas for observation parameters
#' }
#'
#' Methods include:
#' \itemize{
#'  \item plot_dist: Plot histogram of observations, overlaid with the pdf
#'  of the specified distribution for that data stream. Helpful to select
#'  initial parameter values for model fitting.
#'  \item update_par: Update parameters to par
#'  \item update_wpar: Update working parameters to wpar
#'  \item make_X: Create design matrix for observation model
#'  \item obs_var: Data frame of observed (response) variables
#'  \item n2w: Transform parameters from natural to working scale
#'  (for model with no covariates)
#'  \item w2n: Transform parameters from working to natural scale
#'  (for model with no covariates)
#' }

Observation <- R6Class(
  classname = "Observation",
  
  public = list(
    initialize = function(data, dists, par = NULL, wpar = NULL, 
                          wpar_re = NULL, formulas = NULL) {
      private$data_ <- data
      private$dists_ <- dists
      if(is.null(formulas)) {
        # Case with no covariates
        private$par_ <- par 
        private$tpar_ <- self$n2w(par)
        private$formulas_ <- lapply(par, function(varpar) {
          f <- lapply(varpar, function(...) {
            return(~1) # Set all formulas to ~1
          })
          return(f)
        })
      } else if(is.null(wpar)) {
        stop("'wpar' needs to be specified if covariates in observation parameters")
      } else {
        # Case with covariates
        private$tpar_ <- wpar
        private$tpar_re_ <- wpar_re
        private$formulas_ <- formulas        
      }
    },
    
    # Accessors
    data = function() {return(private$data_)},
    dists = function() {return(private$dists_)},
    par = function() {return(private$par_)},
    tpar = function() {return(private$tpar_)},
    tpar_re = function() {return(private$tpar_re_)},
    formulas = function() {return(private$formulas_)},
    
    # Mutators
    update_par = function(par) {
      private$par_ <- par
      private$tpar_ <- self$n2w(par)
    },
    update_wpar = function(wpar, n_state) {
      private$tpar_ <- wpar
      if(!all(rapply(private$formulas_, function(f) { f == ~1 }))) {
        # Only update natural parameters if no covariates
        private$par_ <- self$w2n(wpar, n_state)
      }
    },
    
    # Data frame of response variables
    obs_var = function() {
      obs_names <- names(self$dists())
      obs_var <- self$data()$data()[, obs_names]
      return(obs_var)
    },

    # Create model matrices (same for all states for now)
    make_mat = function(n_states) {
      # Initialise lists of matrices
      X_list_fe <- list()
      X_list_re <- list()
      S_list <- list()
      ncol_re <- NULL
      k <- 1
      
      # Loop over variables
      for(varforms in self$formulas()) {
        # Loop over parameters
        for(form in varforms) {
          # Create matrices based on formula for this parameter
          gam_setup <- gam(formula = update(form, dummy ~ .), 
                           data = cbind(dummy = 1, self$data()$data()), 
                           fit = FALSE)
          
          # Fixed effects design matrix
          X_list_fe[[k]] <- gam_setup$X[, 1:gam_setup$nsdf, drop = FALSE]
          
          # Random effects design matrix
          X_list_re[[k]] <- gam_setup$X[, -(1:gam_setup$nsdf), drop = FALSE]
          
          # Smoothing matrix
          S_list[[k]] <- bdiag_check(gam_setup$S)
          
          # Number of columns for each random effect (rep to duplicate for each state)
          if(length(gam_setup$S) > 0)
            ncol_re <- c(ncol_re, rep(sapply(gam_setup$S, ncol), n_states))
          
          k <- k + 1
        }
      }
      
      # Store as block diagonal matrices
      X_fe <- bdiag_check(rep(X_list_fe, each = n_states))
      X_re <- bdiag_check(rep(X_list_re, each = n_states))
      S <- bdiag_check(rep(S_list, each = n_states))
      
      return(list(X_fe = X_fe, X_re = X_re, S = S, ncol_re = ncol_re))
    },
    
    # Natural to working parameter transformation
    # (No covariates)
    n2w = function(par) {
      wpar <- lapply(1:length(self$dists()), 
                     function(i) dists[[i]]$n2w(par[[i]]))
      names(wpar) <- names(par)
      wpar <- unlist(wpar)
      return(wpar)
    },
    
    # Working to natural parameter transformation
    # (No covariates)
    w2n = function(wpar, n_state) {
      # Initialise list of natural parameters
      par <- list()
      
      # Number of observed variables
      nvar <- length(self$dists())
      
      # Counter to subset observation parameters
      par_count <- 1
      
      # Loop over observed variables
      for(var in 1:nvar) {
        # Number of parameters for this distribution
        npar <- length(self$dists()[[var]]$link())
        # Subset and transform working parameters
        sub_wpar <- wpar[par_count:(par_count + npar*n_state - 1)]
        par_count <- par_count + npar*n_state
        par[[var]] <- self$dists()[[var]]$w2n(sub_wpar)
      }
      
      names(par) <- names(self$dists())
      return(par)
    },
    
    # Histogram of observations with overlaid pdf
    plot_dist = function(name, par = NULL) {
      # Extract observed values for relevant variable
      obs <- self$data()$data()[[name]]
      
      # Histogram of observations
      hist(obs, col = "lightgrey", border = "white", prob = TRUE, 
           main = "", xlab = name)
      
      # Create list of arguments for pdf
      grid <- seq(min(obs, na.rm = TRUE), max(obs, na.rm = TRUE), length = 1e3)
      args <- list(grid)
      if(!is.null(par)) {
        # if parameter values provided by user
        for(i in 1:length(par))
          args[[i+1]] <- par[[i]]        
      } else {
        # else, use default parameter values
        for(i in 1:length(self$par()))
          args[[i+1]] <- self$par()[[name]][[i]]
      }
      
      # Add pdf to histogram plot
      points(grid, do.call(self$dists()[[name]]$pdf(), args), type = "l")
    }
  ),
  
  private = list(
    data_ = NULL,
    dists_ = NULL,
    par_ = NULL,
    tpar_ = NULL,
    tpar_re_ = NULL,
    formulas_ = NULL
  )
)
