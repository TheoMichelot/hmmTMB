
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
    #' @param data HMMData object
    #' @param dists Named list of Distribution objects for each data stream
    #' @param n_states Number of states (needed to construct model formulas)
    #' @param par List of observation parameters (for covariate-free model)
    #' @param coeff_fe Vector of fixed effect parameters on working scale
    #' @param coeff_re Vector of random effect parameters. Defaults to a
    #' vector of zeros if not provided.
    #' @param formulas List of formulas for observation parameters
    #' 
    #' @return A new Observation object
    initialize = function(data, dists, n_states, par = NULL, coeff_fe = NULL, 
                          coeff_re = NULL, formulas = NULL) {
      private$check_args(data = data, dists = dists, n_states = n_states, par = par,
                         coeff_fe = coeff_fe, coeff_re = coeff_re, formulas = formulas)
      
      private$data_ <- data
      private$nstates_ <- n_states
      
      # Define observation distributions
      if(all(sapply(dists, is.character))) {
        # If distributions passed as strings (i.e., names), get corresponding
        # Dist objects
        private$dists_ <- lapply(dists, function(name) dist_list[[name]])
      } else {
        private$dists_ <- dists        
      }
      
      # Set formulas
      if(is.null(formulas)) {
        # Default: set all formulas to ~1
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
      } else {
        private$formulas_ <- make_formulas(formulas, n_states = n_states)        
      }
      
      # Save terms of model formulas
      mats <- self$make_mat()
      ncol_fe <- mats$ncol_fe
      ncol_re <- mats$ncol_re
      private$terms_ <- list(ncol_fe = ncol_fe,
                             ncol_re = ncol_re,
                             names_fe = colnames(mats$X_fe),
                             names_re_all = colnames(mats$X_re),
                             names_re = names(ncol_re))
      
      # Initialise parameters      
      self$update_coeff_fe(rep(0, sum(ncol_fe)))
      self$update_coeff_re(rep(0, sum(ncol_re)))
      self$update_lambda(rep(1, length(ncol_re)))
      
      # Fixed effect parameters     
      if(!is.null(par)) {
        self$update_par(par)
      } else if(!is.null(coeff_fe)) {
        self$update_coeff_fe(coeff_fe)
      } else {
        stop("Either 'par' or 'coeff_fe' must be provided")
      }
      
      # Random effect parameters
      if(!is.null(coeff_re)) {
        self$update_re(coeff_re)
      }
    },
    
    ###############
    ## Accessors ##
    ###############
    #' @description HMMData object
    data = function() {return(private$data_)},
    
    #' @description List of distributions
    dists = function() {return(private$dists_)},
    
    #' @description Number of states
    nstates = function() {return(private$nstates_)},
    
    #' @description Parameters on natural scale
    par = function() {return(private$par_)},
    
    #' @description Fixed effect parameters on working scale
    coeff_fe = function() {return(private$coeff_fe_)},
    
    #' @description Random effect parameters
    coeff_re = function() {return(private$coeff_re_)},
    
    #' @description Smoothness parameters
    lambda = function() {return(private$lambda_)},
    
    #' @description Variance components of smooth terms
    #' 
    #' @details This function transforms the smoothness parameter of
    #' each smooth term into a standard deviation, given by 
    #' SD = 1/sqrt(lambda). It is particularly helpful to get the
    #' standard deviations of independent normal random effects.
    vcomp = function() {return(1/sqrt(private$lambda_))},
    
    #' @description List of model formulas for observation model
    formulas = function() {return(private$formulas_)},
    
    #' @description Terms of model formulas
    terms = function() {return(private$terms_)},
    
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
    #' Updates the 'par' attribute to the list passed as input,
    #' and updates the intercept elements of 'coeff_fe' using
    #' the list passed as input
    #' 
    #' @param par New list of parameters
    update_par = function(par) {
      private$par_ <- par
      
      # Get index of first column of X_fe for each parameter
      ncol_fe <- self$terms()$ncol_fe
      n_par <- length(ncol_fe)
      i0 <- c(1, cumsum(ncol_fe)[-n_par] + 1)
      
      # Apply link to get parameters on working scale
      private$coeff_fe_[i0] <- self$n2w(par)
    },
    
    #' @description Update coefficients for fixed effect parameters
    #' 
    #' @param coeff_fe New vector of coefficients for fixed effect 
    #' parameters
    update_coeff_fe = function(coeff_fe) {
      private$coeff_fe_ <- matrix(coeff_fe)
      rownames(private$coeff_fe_) <- self$terms()$names_fe
      if(all(rapply(self$formulas(), function(f) { f == ~1 }))) {
        # Only update par if no covariates
        private$par_ <- self$w2n(coeff_fe)
      }
    },
    
    #' @description Update random effect parameters
    #' 
    #' @param coeff_re New vector of coefficients for random effect 
    #' parameters
    update_coeff_re = function(coeff_re) {
      private$coeff_re_ <- matrix(coeff_re)
      rownames(private$coeff_re_) <- self$terms()$names_re_all
    },
    
    #' @description Update smoothness parameters
    #' 
    #' @param lambda New smoothness parameter vector
    update_lambda = function(lambda) {
      private$lambda_ <- matrix(lambda)
      rownames(private$lambda_) <- self$terms()$names_re
    },
    
    ###################
    ## Other methods ##
    ###################
    #' @description Make model matrices
    #' 
    #' @param new_data Optional new data set, including covariates for which
    #' the design matrices should be created. If this argument is not specified,
    #' the design matrices are based on the original data frame. 
    #' 
    #' @return A list with elements:
    #' \itemize{
    #'   \item{X_fe}{Design matrix for fixed effects}
    #'   \item{X_re}{Design matrix for random effects}
    #'   \item{S}{Smoothness matrix for random effects}
    #'   \item{ncol_fe}{Number of columns of X_fe for each parameter}
    #'   \item{ncol_re}{Number of columns of X_re and S for each random effect}
    #' }
    make_mat = function(new_data = NULL) {
      make_matrices(formulas = self$formulas(),
                    data = self$data()$data(),
                    new_data = new_data)
    },
    
    #' Design matrices for grid of covariates
    #' 
    #' @param var Name of variable
    #' @param covs Optional data frame with a single row and one column
    #' for each covariate, giving the values that should be used. If this is
    #' not specified, the mean value is used for numeric variables, and the
    #' first level for factor variables.
    #' @param n_grid Grid size (number of points). Default: 1000.
    #' 
    #' @return A list with the same elements as the output of make_mat, 
    #' plus a data frame of covariates values.
    make_mat_grid = function(var, covs = NULL, n_grid = 1e3) {
      # Data frame for covariate grid
      new_data <- cov_grid(var = var, data = self$data()$data(), 
                           covs = covs, formulas = self$formulas(),
                           n_grid = n_grid)
      
      # Create design matrices
      mats <- self$make_mat(new_data = new_data)
      
      # Save data frame of covariate values
      mats$new_data <- new_data
      
      return(mats)
    },
    
    #' @description Natural to working parameter transformation
    #' 
    #' This function applies the link functions of the distribution
    #' parameters, to transform parameters from their natural scale
    #' to the working scale (i.e., linear predictor scale)
    #' 
    #' @param par List of parameters on natural scale
    #' 
    #' @return Vector of parameters on working scale
    n2w = function(par) {
      coeff_fe <- lapply(seq_along(self$dists()), 
                         function(i) self$dists()[[i]]$n2w(par[[i]]))
      names(coeff_fe) <- names(par)
      coeff_fe <- unlist(coeff_fe)
      return(coeff_fe)
    },
    
    #' @description  Working to natural parameter transformation
    #'
    #' This function applies the inverse link functions of the
    #' distribution parameters, to transform parameters from the working
    #' scale (i.e., linear predictor scale) to their natural scale.
    #'
    #' @param coeff_fe Vector of parameters on working scale
    #' 
    #' @return List of parameters on natural scale
    w2n = function(coeff_fe) {
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
        sub_coeff_fe <- coeff_fe[par_count:(par_count + npar*n_states - 1)]
        par_count <- par_count + npar*n_states
        par[[var]] <- self$dists()[[var]]$w2n(sub_coeff_fe)
      }
      
      names(par) <- names(self$dists())
      return(par)
    },
    
    #' @description Get observation parameters from design matrices
    #' 
    #' @param X_fe Design matrix for fixed effects, as returned
    #' by \code{make_mat}
    #' @param X_re Design matrix for random effects, as returned
    #' by \code{make_mat}
    #' @param coeff_fe Optional vector of coefficients for fixed effect
    #' parameters. If this isn't provided, the model parameters
    #' are used.
    #' @param coeff_re Optional vector of coefficients for random 
    #' effect parameters. If this isn't provided, the model parameters 
    #' are used.
    #' @param full_names Logical. If TRUE, the rows of the output
    #' are named in the format "variable.parameter" (default). If
    #' FALSE, the rows are names in the format "parameter". The
    #' latter is used in various internal functions, when the parameters
    #' need to be passed on to an R function.
    #' 
    #' @return Array with one slice for each time step, one row 
    #' for each observation parameter, and one column for each state.
    par_all = function(X_fe, X_re, coeff_fe = NULL, coeff_re = NULL, 
                       full_names = TRUE) {
      # Number of states
      n_states <- self$nstates()
      
      # Number of parameters on natural scale (in each state)
      n_par <- sum(sapply(self$dists(), function(d) d$npar()))
      
      # Define parameters
      if(length(coeff_fe) == 0)
        coeff_fe <- self$coeff_fe()
      if(length(coeff_re) == 0)
        coeff_re <- self$coeff_re()
      
      # Get linear predictor
      lp <- X_fe %*% coeff_fe + X_re %*% coeff_re
      lp_mat <- matrix(lp, ncol = n_par * n_states)
      
      # Number of observations
      n <- nrow(lp_mat)
      
      # Matrix of natural parameters
      par_mat <- apply(lp_mat, 1, function(lp_vec) {
        par_ls <- self$w2n(lp_vec)
        par_vec <- unlist(par_ls, use.names = FALSE)
        return(par_vec)
      })
      
      # Array of natural parameters
      par_array <- array(par_mat, dim = c(n_states, n_par, n))
      
      # Hacky way to get parameter names
      if(full_names) {
        par_names <- names(unlist(rapply(self$w2n(lp_mat[1,]),
                                         function(v) v[1])))
      } else {
        par_names <- unlist(lapply(self$w2n(lp_mat[1,]), names))
      }
      
      # Set dimension names for rows and columns
      dimnames(par_array) <- list(paste("state", 1:n_states),
                                  par_names,
                                  NULL)
      
      # Transpose each slice
      par_array <- aperm(par_array, perm = c(2, 1, 3))
      
      return(par_array)
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
      
      # State-dependent parameters
      par <- self$par_all(X_fe = X_fe, X_re = X_re, full_names = FALSE)
      
      # Initialise matrix of probabilities to 1
      prob <- matrix(1, nrow = n, ncol = n_states)
      
      # Counter to subset parameter vector
      par_count <- 1
      
      # Loop over observed variables
      for(var in 1:n_var) {
        obsdist <- self$dists()[[var]]
        par_ind <- par_count:(par_count + obsdist$npar() - 1)
        
        # Loop over observations (rows of prob)
        for (i in 1:n) {
          # Don't update likelihood is observation is missing
          if(!is.na(data[i, var])) {
            # Loop over states (columns of prob)
            for (s in 1:n_states) {
              prob[i, s] <- prob[i, s] * 
                obsdist$pdf_apply(x = data[i, var], par = par[par_ind, s, i])
            }            
          }
        }
        
        par_count <- par_count + obsdist$npar()
      }
      
      return(prob)
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
        geom_line(aes(grid, val, col = state, linetype = state), 
                  data = df_dens, size = 0.7) +
        scale_color_manual("", values = c(hmmTMB_cols[1:self$nstates()], "black")) +
        scale_linetype_manual("", values = c(rep(1, self$nstates()), 2)) +
        coord_cartesian(ylim = c(0, 1.1 * max(h$density))) +
        theme_light()
      
      return(p)
    },
    
    #' @description Print model formulation
    formulation = function() {
      cat("#######################\n")
      cat("## Observation model ##\n")
      cat("#######################\n")
      # List of distribution names
      d_list <- lapply(self$dists(), function(d) d$name())
      # List of parameter names
      p_list <- lapply(self$dists(), function(d) names(d$link()))
      
      # List of parameter formulas
      s_list <- lapply(self$formulas(), function(f) {
        ff <- unlist(f)
        s <- NULL
        for(i in seq_along(ff))
          s <- paste0(s, "  * ", names(ff)[i], " ~ ", as.character(ff[[i]])[2], "\n")
        return(s)
      })
      
      # Loop over observed variables
      for(i in seq_along(d_list)) {
        # Print variable distribution
        cat(paste0("+ ", names(d_list)[i], " ~ ", d_list[[i]], "(", 
                   paste0(p_list[[i]], collapse = ", "), ")"), "\n")
        # Print parameter formulas
        cat(s_list[[i]], "\n")
      }
    }
  ),
  
  private = list(
    ################
    ## Attributes ##
    ################
    data_ = NULL,
    dists_ = NULL,
    nstates_ = NULL,
    par_ = NULL,
    coeff_fe_ = NULL,
    coeff_re_ = NULL,
    lambda_ = NULL,
    formulas_ = NULL,
    terms_ = NULL,
    
    #################################
    ## Check constructor arguments ##
    #################################
    # (For argument description, see constructor)
    check_args = function(data, dists, n_states, par, coeff_fe, coeff_re, formulas) {
      if(!inherits(data, "HMMData")) {
        stop("'data' should be a HMMData object")
      }
      
      if(!is.list(dists)) {
        stop("'dists' should be a list")
      }
      
      if(!all(sapply(dists, inherits, "Dist")) & !all(sapply(dists, is.character))) {
        stop(paste("Elements of 'dists' should all be either character strings",
                   "(i.e., distribution names), or Dist objects"))
      }
      
      if(!all(names(dists) %in% colnames(data$data()))) {
        stop("Variable name in 'dists' not found in data")
      }
      
      if(!is.numeric(n_states) | n_states < 1) {
        stop("'n_states' should be a numeric >= 1")
      }
      
      if(!is.null(par)) {
        if(!is.list(par) | length(par) != length(dists)) {
          stop("'par' should be a list of same length as 'dists'")
        }
        
        if(!all(rapply(par, length) == n_states) | !all(rapply(par, is.numeric))) {
          stop("Elements of 'par' should be numeric vectors of length 'n_states'")
        }
        
        if(!all(names(par) == names(dists))) {
          stop("'par' should have the same names as 'dists'")
        }
      }
      
      if(!is.null(coeff_fe)) {
        if(!is.numeric(coeff_fe) | !is.vector(coeff_fe)) {
          stop("'coeff_fe' should be a numeric vector")
        }
      }
      
      if(!is.null(coeff_re)) {
        if(!is.numeric(coeff_re) | !is.vector(coeff_re)) {
          stop("'coeff_re' should be a numeric vector")
        }
      }
      
      if(!is.null(formulas)) {
        if(!is.list(formulas) |
           !all(rapply(formulas, function(x) inherits(x, "formula")))) {
          stop("'formulas' should be a list of R formulas")
        }
        
        if(!all(names(formulas) == names(dists))) {
          stop("'formulas' should have the same names as 'dists'")
        }
      }
    }
  )
)
