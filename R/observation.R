
#' R6 class for HMM observation model
#'
#' Contains the data, distributions, parameters, and formulas for
#' the observation model from a hidden Markov model.
Observation <- R6Class(
  classname = "Observation",
  
  public = list(
    #################
    ## Constructor ##
    #################
    #' @description Create new Observation object
    #' 
    #' @param data HmmData object
    #' @param dists Named list of Distribution objects for each data stream
    #' @param n_states Number of states (needed to construct model formulas)
    #' @param par List of observation parameters (for covariate-free model)
    #' @param wpar Vector of fixed effect parameters on working scale
    #' @param wpar_re Vector of random effect parameters. Defaults to a
    #' vector of zeros if not provided.
    #' @param formulas List of formulas for observation parameters
    #' 
    #' @return A new Observation object
    initialize = function(data, dists, n_states, par = NULL, wpar = NULL, 
                          wpar_re = NULL, formulas = NULL) {
      private$data_ <- data
      private$dists_ <- dists
      private$nstates_ <- n_states
      
      if(is.null(formulas)) {
        # Case with no covariates
        private$par_ <- par 
        private$wpar_ <- self$n2w(par)
        
        # Set all formulas to ~1
        private$formulas_ <- lapply(par, function(varpar) {
          f <- lapply(varpar, function(...) {
            g <- lapply(1:n_states, function(...) {
              return(~1)              
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
        private$wpar_ <- wpar
        private$formulas_ <- make_formulas(formulas, n_states = n_states)        
      }
      
      # Initialise random effect parameters
      mats <- self$make_mat()
      if(ncol(mats$X_re) == 0) {
        # integer(0) rather than NULL so that X_re %*% wpar_re 
        # is valid when X_re has zero columns
        private$wpar_re_ <- integer(0) 
      } else if(is.null(wpar_re)) {
        # if no value provided, wpar_re initialised to vector of zeros
        private$wpar_re_ <- rep(0, ncol(mats$X_re))
      } else {
        private$wpar_re_ <- wpar_re
      }
    },
    
    ###############
    ## Accessors ##
    ###############
    #' @description HmmData object
    data = function() {return(private$data_)},
    
    #' @description List of distributions
    dists = function() {return(private$dists_)},
    
    #' @description Number of states
    nstates = function() {return(private$nstates_)},
    
    #' @description Parameters on natural scale
    par = function() {return(private$par_)},
    
    #' @description Fixed effect parameters on working scale
    wpar = function() {return(private$wpar_)},
    
    #' @description Random effect parameters
    wpar_re = function() {return(private$wpar_re_)},
    
    #' @description List of model formulas for observation model
    formulas = function() {return(private$formulas_)},
    
    #' @description  Data frame of response variables
    obs_var = function() {
      obs_names <- names(self$dists())
      obs_var <- self$data()$data()[, obs_names, drop = FALSE]
      return(obs_var)
    },
    
    ##############
    ## Mutators ##
    ##############
    #' @description Update parameters
    #' 
    #' @param par New list of parameters
    update_par = function(par) {
      private$par_ <- par
      private$wpar_ <- self$n2w(par)
    },
    
    #' @description Update fixed effect parameters on working scale
    #' 
    #' @param wpar New vector of fixed effect parameters on working scale
    update_wpar = function(wpar) {
      names(wpar) <- NULL
      private$wpar_ <- wpar
      if(all(rapply(self$formulas(), function(f) { f == ~1 }))) {
        # Only update natural parameters if no covariates
        private$par_ <- self$w2n(wpar)
      }
    },
    
    #' @description Update random effect parameters
    #' 
    #' @param wpar_re New vector og random effect parameters
    update_wpar_re = function(wpar_re) {
      private$wpar_re_ <- wpar_re
    },
    
    ###################
    ## Other methods ##
    ###################
    #' @description Make model matrices
    #' 
    #' @return A list with elements:
    #' \itemize{
    #'   \item{X_fe}{Design matrix for fixed effects}
    #'   \item{X_re}{Design matrix for random effects}
    #'   \item{S}{Smoothness matrix for random effects}
    #'   \item{ncol_re}{Number of columns of X_re and S for each random effect}
    #' }
    make_mat = function() {
      make_mat_obs(formulas = self$formulas(),
                   data = self$data()$data())
    },
    
    #' @description Observation likelihoods
    #' 
    #' @param X_fe Design matrix for fixed effects
    #' @param X_re Design matrix for random effects
    #' 
    #' @return Matrix of likelihoods of observations, with one row for each 
    #' time step, and one column for each state.
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
      wpar <- X_fe %*% self$wpar() + X_re %*% self$wpar_re()
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
          # Don't update likelihood is observation is missing
          if(!is.na(data[i, var])) {
            # Subset and transform observation parameters
            sub_wpar <- par_mat[i, par_count:(par_count + obsdist$npar() * n_states - 1)]
            par <- obsdist$w2n(sub_wpar, as_matrix = TRUE)
            
            # Loop over states (columns)
            for (s in 1:n_states) {
              # Vector of parameters for state s
              subpar <- par[s,]
              
              prob[i, s] <- prob[i, s] * obsdist$pdf_apply(x = data[i, var], par = subpar)
            }            
          }
        }
        
        par_count <- par_count + obsdist$npar() * n_states
      }
      
      return(prob)
    },
    
    #' @description Natural to working parameter transformation
    #' 
    #' @param par List of parameters on natural scale
    #' 
    #' @return Vector of parameters on working scale
    n2w = function(par) {
      wpar <- lapply(1:length(self$dists()), 
                     function(i) dists[[i]]$n2w(par[[i]]))
      names(wpar) <- names(par)
      wpar <- unlist(wpar)
      return(wpar)
    },
    
    #' @description  Working to natural parameter transformation
    #'
    #' @param wpar Vector of parameters on working scale
    #' 
    #' @return List of parameters on natural scale
    w2n = function(wpar) {
      # Initialise list of natural parameters
      par <- list()
      
      # Number of observed variables
      nvar <- length(self$dists())
      # Number of states
      n_states <- self$nstates()
      
      # Counter to subset observation parameters
      par_count <- 1
      
      # Loop over observed variables
      for(var in 1:nvar) {
        # Number of parameters for this distribution
        npar <- length(self$dists()[[var]]$link())
        # Subset and transform working parameters
        sub_wpar <- wpar[par_count:(par_count + npar*n_states - 1)]
        par_count <- par_count + npar*n_states
        par[[var]] <- self$dists()[[var]]$w2n(sub_wpar)
      }
      
      names(par) <- names(self$dists())
      return(par)
    },

    #' @description Plot histogram of data and pdfs
    #' 
    #' Plot histogram of observations for the variable specified by the argument name, 
    #' overlaid with the pdf of the specified distribution for that data stream. 
    #' Helpful to select initial parameter values for model fitting, or to visualise 
    #' fitted state-dependent distributions.
    #' 
    #' @param name Name of response variable for which the histogram
    #' and pdfs should be plotted.
    #' @param par Optional matrix of parameters, with one row for each state
    #' and one column for each parameter. The columns of the matrix should be
    #' named with the names of the parameters, e.g. "mean" and "sd" for normal
    #' distribution. If not provided, the parameters stored in the object are
    #' used (default).
    #' @param weights Optional vector of length the number of pdfs that are
    #' plotted. Useful to visualise a mixture of distributions weighted by the
    #' proportion of time spent in the different states.
    #' 
    #' @return A ggplot object
    plot_dist = function(name, par = NULL, weights = NULL) {
      # Colour palette
      pal <- c("#00798c", "#d1495b", "#edae49", "#66a182", "#2e4057", "#8d96a3")
      
      # Extract observed values for relevant variable
      obs <- data.frame(val = self$data()$data()[[name]])
      
      # Matrix of parameters
      if(is.null(par)) {
        par <- matrix(unlist(self$par()[[name]]), nrow = self$nstates())
        colnames(par) <- names(self$dists()[[name]]$link())
      }
      
      # Weights for each state-dependent distribution
      if(is.null(weights))
        weights <- rep(1, self$nstates())
      
      # Grid over range of observed variable
      n_grid <- 1e3
      grid <- seq(min(obs, na.rm = TRUE), max(obs, na.rm = TRUE), length = n_grid)
      
      # Loop over states
      vals <- matrix(NA, nrow = n_grid, ncol = self$nstates() + 1)
      for(state in 1:self$nstates()) {
        # Define list of arguments to pass to pdf
        args <- list(grid)
        args <- c(args, par[state,])
        
        # Compute state-dependent pdf
        vals[,state] <- weights[state] * do.call(self$dists()[[name]]$pdf(), args)
      }
      # Weighted sum of state-dependent pdfs
      vals[, self$nstates() + 1] <- rowSums(vals[, 1:self$nstates()])
      
      # Data frame of state-dependent densities
      df_dens <- data.frame(
        state = rep(c(paste("State", 1:self$nstates()), "Total"), each = n_grid),
        grid = grid,
        val = as.vector(vals))
      
      # Create hist object to extract ylim
      breaks <- seq(min(obs, na.rm = TRUE), max(obs, na.rm = TRUE), length = 20)
      h <- hist(obs$val, breaks = breaks, plot = FALSE)
      
      # Create ggplot histogram
      p <- ggplot(obs, aes(x = val)) + xlab(name) +
        geom_histogram(breaks = breaks, aes(y=..density..), 
                       col = "white", bg = "lightgrey", na.rm = TRUE) + 
        geom_line(aes(grid, val, col = state, linetype = state), data = df_dens, size = 0.7) +
        scale_color_manual("", values = c(pal[1:self$nstates()], "black")) +
        scale_linetype_manual("", values = c(rep(1, self$nstates()), 2)) +
        coord_cartesian(ylim = c(0, 1.1 * max(h$density))) +
        theme_light()
      
      return(p)
    }
  ),
  
  private = list(
    data_ = NULL,
    dists_ = NULL,
    nstates_ = NULL,
    par_ = NULL,
    wpar_ = NULL,
    wpar_re_ = NULL,
    formulas_ = NULL
  )
)
