
#' Hidden Markov observation class
#'
#' @description Encapsulates the observation model from a hidden Markov model.
#' Object can be created using $new with arguments:
#' \itemize{
#'   \item data: a HmmData object
#'   \item dists: named list of distributions for each data stream
#'   \item n_states: number of states (needed to construct model formulas)
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
    initialize = function(data, dists, n_states, par = NULL, wpar = NULL, 
                          wpar_re = NULL, formulas = NULL) {
      private$data_ <- data
      private$dists_ <- dists
      private$nstates_ <- n_states
      if(is.null(formulas)) {
        # Case with no covariates
        private$par_ <- par 
        private$tpar_ <- self$n2w(par)
        private$formulas_ <- lapply(par, function(varpar) {
          f <- lapply(varpar, function(...) {
            g <- lapply(1:n_states, function(...) {
              return(~1) # Set all formulas to ~1              
            })
            names(g) <- paste0("state", 1:n_states)
            return(g)
          })
          return(f)
        })
      } else if(is.null(wpar)) {
        stop("'wpar' needs to be specified if covariates in observation parameters")
      } else {
        # Case with covariates
        private$tpar_ <- wpar
        private$tpar_re_ <- wpar_re
        private$formulas_ <- make_formulas(formulas, n_states = n_states)        
      }
    },
    
    # Accessors
    data = function() {return(private$data_)},
    dists = function() {return(private$dists_)},
    nstates = function() {return(private$nstates_)},
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
      if(all(rapply(private$formulas_, function(f) { f == ~1 }))) {
        # Only update natural parameters if no covariates
        private$par_ <- self$w2n(wpar, n_state)
      }
    },
    update_wpar_re = function(wpar_re) {
      private$tpar_re_ <- wpar_re
    },
    
    # Data frame of response variables
    obs_var = function() {
      obs_names <- names(self$dists())
      obs_var <- self$data()$data()[, obs_names, drop = FALSE]
      return(obs_var)
    },
    
    # Create model matrices
    make_mat = function() {
      make_mat_obs(formulas = self$formulas(),
                   data = self$data()$data())
    },
    
    # Compute observation probabilities
    obs_probs = function(X_fe, X_re) {
      # Data frame of observations
      data <- self$obs_var()
      
      # Number of observations
      n <- nrow(data)
      # Number of states
      n_states <- self$nstates()
      # Number of variables
      n_var <- ncol(data)
      
      # Matrix of observation parameters
      wpar <- X_fe %*% self$tpar() + X_re %*% self$tpar_re()
      par_mat <- matrix(wpar, nrow = n)
      
      # Initialise matrix of probabilities to 1
      prob <- matrix(1, nrow = n, ncol = n_states)
      
      # Counter to subset parameter vector
      par_count <- 1
      
      # Loop over observed variables
      for(var in 1:n_var) {
        obsdist <- self$dists()[[var]]
        
        # Loop over observations (rows)
        for (i in 1:n) {
          # Subset and transform observation parameters
          sub_wpar <- par_mat[i, par_count:(par_count + obsdist$npar() * n_states - 1)]
          par <- obsdist$invlink_apply(sub_wpar, n_states)
          
          # Loop over states (columns)
          for (s in 1:n_states) {
            # Vector of parameters for state s
            subpar <- par[s,]
            
            prob[i, s] <- prob[i, s] * obsdist$pdf_apply(x = data[i, var], par = subpar)
          }
        }
        
        par_count <- par_count + obsdist$npar() * n_states
      }
      
      return(prob)
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
      
      # Matrix of parameters
      if(is.null(par)) {
        par <- self$dists()[[name]]$invlink_apply(wpar = self$tpar(), 
                                                  n_states = self$nstates())
      }
      
      # Grid over range of observed variable
      grid <- seq(min(obs, na.rm = TRUE), max(obs, na.rm = TRUE), length = 1e3)
      
      # Loop over states
      for(state in 1:self$nstates()) {
        # Define list of arguments to pass to pdf
        args <- list(grid)
        args <- c(args, par[state,])

        # Add pdf to histogram plot
        points(grid, do.call(self$dists()[[name]]$pdf(), args), type = "l")        
      }

    }
  ),
  
  private = list(
    data_ = NULL,
    dists_ = NULL,
    nstates_ = NULL,
    par_ = NULL,
    tpar_ = NULL,
    tpar_re_ = NULL,
    formulas_ = NULL
  )
)
