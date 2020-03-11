
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
    initialize = function(data, dists, par = NULL, 
                          wpar = NULL, formulas = NULL) {
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
        private$formulas_ <- formulas        
      }
    },
    
    # Accessors
    data = function() {return(private$data_)},
    dists = function() {return(private$dists_)},
    par = function() {return(private$par_)},
    tpar = function() {return(private$tpar_)},
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
      obs_names <- names(private$dists_)
      obs_var <- private$data_$data()[, obs_names]
      return(obs_var)
    },

    # Create block-diagonal design matrix (same for all states for now)
    make_X = function(n_states) {
      # List of design matrices (one for each parameter of each variable)
      X_list <- unlist(lapply(private$formulas_, function(varforms) {
        lapply(varforms, function(form) {
          # Use mgcv to create model matrices for each parameter
          gam_setup <- gam(formula = update(form, dummy ~ .), 
                           data = cbind(dummy = 1, private$data_$data()), 
                           fit = FALSE)
          return(gam_setup$X)
        })
      }), recursive = FALSE)        
      
      # Create block diagonal matrix
      X <- bdiag(rep(X_list, each = n_states))
      return(X)
    },
    
    # Natural to working parameter transformation
    # (No covariates)
    n2w = function(par) {
      wpar <- lapply(1:length(private$dists_), 
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
      nvar <- length(private$dists_)
      
      # Counter to subset observation parameters
      par_count <- 1
      
      # Loop over observed variables
      for(var in 1:nvar) {
        # Number of parameters for this distribution
        npar <- length(private$dists_[[var]]$link())
        # Subset and transform working parameters
        sub_wpar <- wpar[par_count:(par_count + npar*n_state - 1)]
        par_count <- par_count + npar*n_state
        par[[var]] <- private$dists_[[var]]$w2n(sub_wpar)
      }
      
      names(par) <- names(private$dists_)
      return(par)
    },
    
    # Histogram of observations with overlaid pdf
    plot_dist = function(name, par = NULL) {
      # Extract observed values for relevant variable
      obs <- private$data_$data()[[name]]
      
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
        for(i in 1:length(private$par_))
          args[[i+1]] <- private$par_[[name]][[i]]
      }
      
      # Add pdf to histogram plot
      points(grid, do.call(private$dists_[[name]]$pdf(), args), type = "l")
    }
  ),
  
  private = list(
    data_ = NULL,
    dists_ = NULL,
    par_ = NULL,
    tpar_ = NULL,
    formulas_ = NULL
  )
)
