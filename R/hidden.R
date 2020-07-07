
#' R6 class for HMM hidden process model
#'
#' Contains the parameters and model formulas for the hidden process model.
MarkovChain <- R6Class(
  classname = "MarkovChain",
  
  public = list(
    #################
    ## Constructor ##
    #################
    #' @description Create new MarkovChain object
    #' 
    #' @param n_states Number of states. If not specified, then \code{structure} 
    #' needs to be provided as a matrix, and n_states is deduced from its dimensions.
    #' @param structure Either (1) matrix with an entry of "." on diagonal, a "0" for 
    #' transitions that are not allowed (not implemented yet), and a formula "~1" 
    #' for covariates affecting transitions that are to be estimated, or (2) single
    #' formula, assumed for all transition probabilities. (Default: no covariate
    #' dependence.)
    #' @param tpm Initial transition probability matrix. (Default: 0.9 on diagonal,
    #' and 0.1/(n_states - 1) for all other entries.)
    #' @param par Initial parameters of the state process, on the working scale.
    #' 
    #' @return A new MarkovChain object
    initialize = function(n_states = NULL, structure = NULL, tpm = NULL, par = NULL) {
      
      if(is.null(structure)) {
        # No covariate effects
        structure <- matrix("~1", nrow = n_states, ncol = n_states)
        diag(structure) <- "."
        private$nstates_ <- n_states
        
        # Set initial tpm (default: 0.9 on diagonal)
        if(is.null(tpm) & is.null(par)) {
          tpm <- matrix(0.1/(n_states-1), nrow = n_states, ncol = n_states)
          diag(tpm) <- 0.9
          par <- private$tpm2par(tpm)
        } else if(is.null(par)) {
          par <- private$tpm2par(tpm)
        } else {
          tpm <- private$par2tpm(par)
        }
        private$tpm_ <- tpm
        private$par_ <- par
      
      } else {
        # Covariate effects
        if(length(structure) == 1 | inherits(structure, "formula")) {
          # Same formula for all transitions
          structure <- matrix(format(structure), 
                              nrow = n_states, 
                              ncol = n_states)
          diag(structure) <- "."
        } else {
          n_states <- nrow(structure)
        }
        
        private$nstates_ <- n_states
        
        # Set initial parameters 
        private$par_ <- par
      }
      
      # Create list of formulas
      ls_form_char <- as.list(structure[!diag(self$nstates())])
      ls_form <- lapply(ls_form_char, function(form_char) {
        if(form_char == ".")
          return(NULL)
        else
          return(as.formula(form_char))
      })
      
      # Set remaining private attributes
      private$structure_ <- structure
      private$formulas_ <- ls_form
      private$par_re_ <- integer(0)  # so that X_re %*% par_re is valid
    },
    
    ###############
    ## Accessors ##
    ###############
    #' @description Structure of MarkovChain model
    structure = function() {return(private$structure_)},
    
    #' @description List of formulas for MarkovChain model
    formulas = function() {return(private$formulas_)},
    
    #' @description Current transition probability matrix
    tpm = function() {return(private$tpm_)},
    
    #' @description Current parameter estimates (fixed effects)
    par = function() {return(private$par_)},
    
    #' @description Current parameter estimates (random effects)
    par_re = function() {return(private$par_re_)},
    
    #' @description Number of states
    nstates = function() {return(private$nstates_)},
    
    ##############
    ## Mutators ##
    ##############
    #' @description Update transition probability matrix
    #' 
    #' @param newtpm New transition probability matrix
    update_tpm = function(newtpm) {
      private$tpm_ <- newtpm
      private$par_ <- private$tpm2par(newtpm)
    },
    
    #' @description Update parameters (fixed effects)
    #' 
    #' @param newpar New parameters (fixed effects)
    update_par = function(newpar) {
      private$par_ <- newpar
      if(all(self$structure() %in% c(".", "~1"))) {
        # Only update tpm if no covariates
        private$tpm_ <- private$par2tpm(newpar)        
      }
    },
    
    #' @description Update parameters (random effects)
    #' 
    #' @param newpar New parameters (random effects)
    update_par_re = function(newpar) {
      private$par_re_ <- newpar
    },
    
    ###################
    ## Other methods ##
    ###################
    #' @description Make model matrices
    #' 
    #' @param data Data frame containing all needed covariates
    #' @param new_data Optional new data set, including covariates for which
    #' the design matrices should be created. This needs to be passed in addition
    #' to the argument '\code{data}', for cases where smooth terms or factor
    #' covariates are included, and the original data set is needed to determine
    #' the full range of covariate values.
    #' 
    #' @return A list with elements:
    #' \itemize{
    #'   \item{X_fe}{Design matrix for fixed effects}
    #'   \item{X_re}{Design matrix for random effects}
    #'   \item{S}{Smoothness matrix for random effects}
    #'   \item{ncol_re}{Number of columns of X_re and S for each random effect}
    #' }
    make_mat = function(data, new_data = NULL) {
      make_mat_hid(formulas = self$formulas(), data = data, new_data = new_data)
    },
    
    #' @description Get transition probability matrices from design matrices
    #' 
    #' @param X_fe Design matrix for fixed effects, as returned
    #' by \code{make_mat}
    #' @param X_re Design matrix for random effects, as returned
    #' by \code{make_mat}
    #' 
    #' @return Array with one slice for each transition probability matrix
    tpm_all = function(X_fe, X_re) {
      n_states <- self$nstates()
      ltpm <- X_fe %*% self$par() + X_re %*% self$par_re()
      ltpm_mat <- matrix(ltpm, ncol = n_states * (n_states - 1))
      tpm <- apply(ltpm_mat, 1, private$par2tpm)
      tpm <- array(tpm, dim = c(n_states, n_states, nrow(ltpm_mat)))
      return(tpm)
    },
    
    #' @description Stationary distributions
    #' 
    #' @param X_fe Design matrix for fixed effects, as returned
    #' by \code{make_mat}
    #' @param X_re Design matrix for random effects, as returned
    #' by \code{make_mat}
    #'
    #' @return Matrix of stationary distributions. Each row corresponds to
    #' a row of the design matrices, and each column corresponds to a state.
    stat_dists = function(X_fe, X_re) {
      # Number of states
      n_states <- self$nstates()
      
      # Derive transition probability matrices
      tpms <- self$tpm_all(X_fe = X_fe, X_re = X_re)
      
      tryCatch({
        # For each transition matrix, get corresponding stationary distribution
        stat_dists <- apply(tpms, 3, function(tpm)
          solve(t(diag(n_states) - tpm + 1), rep(1, n_states)))
        stat_dists <- t(stat_dists)
      },
      error = function(e) {
        stop(paste("The stationary distributions cannot be calculated",
                   "for these covariate values (singular system)."))
      })
      
      return(stat_dists)
    },
    
    #' @description Design matrices for grid of covariates
    #' 
    #' Used in plotting functions such as Hmm$plot_tpm and Hmm$plot_stat_dist
    #' 
    #' @param var Name of variable
    #' @param data Data frame containing the covariates
    #' @param covs Optional data frame with a single row and one column
    #' for each covariate, giving the values that should be used. If this is
    #' not specified, the mean value is used for numeric variables, and the
    #' first level for factor variables.
    #' 
    #' @return A list with the same elements as the output of make_mat, 
    #' plus a data frame of covariates values.
    make_mat_grid = function(var, data, covs = NULL) {
      # Data frame for covariate grid
      new_data <- cov_grid(var = var, data = data, covs = covs, 
                           formulas = self$formulas())
      
      # Create design matrices
      mats <- self$make_mat(data = data, new_data = new_data)
      
      # Save data frame of covariate values
      mats$new_data <- new_data
      
      return(mats)
    },
    
    #' @description Simulate from Markov chain
    #' 
    #' @param n Number of time steps to simulate
    #' @param data Optional data frame containing all needed covariates
    #' @param new_data Optional new data set, including covariates for which
    #' the design matrices should be created. This needs to be passed in addition
    #' to the argument '\code{data}', for cases where smooth terms or factor
    #' covariates are included, and the original data set is needed to determine
    #' the full range of covariate values.
    #' 
    #' @return Sequence of states of simulated chain
    simulate = function(n, data = NULL, new_data = NULL) {
      # Number of states
      n_states <- self$nstates()
      
      # Time series ID
      if(is.null(new_data)) {
        ID <- rep(1, n)
      } else {
        ID <- new_data$ID
      }
      
      # If no covariates, 'data' is optional
      if(is.null(data)) {
        data <- data.frame(ID = ID)
      }
      
      # Create transition probability matrices
      mats_hid <- self$make_mat(data = data, new_data = new_data)
      tpms <- self$tpm_all(X_fe = mats_hid$X_fe, X_re = mats_hid$X_re)
      
      # Uniform initial distribution for now
      delta <- rep(1/n_states, n_states) 
      
      # Simulate state process      
      S <- rep(NA, n)
      S[1] <- sample(1:n_states, size = 1, prob = delta)
      for(i in 2:n) {
        if(round(i/n*100)%%10 == 0) {
          cat("\rSimulating states... ", round(i/n*100), "%", sep = "")        
        }
        
        if(ID[i] != ID[i-1]) {
          S[i] <- sample(1:n_states, size = 1, prob = delta)
        } else {
          S[i] <- sample(1:n_states, size = 1, prob = tpms[S[i-1], , i-1])          
        }
      }
      cat("\n")
      
      return(S)
    }
  ),
  
  private = list(
    structure_ = NULL,
    formulas_ = NULL,
    par_ = NULL,
    par_re_ = NULL,
    tpm_ = NULL,
    nstates_ = NULL,
    
    check_structure = function() {
      if (!all(diag(self$structure()) == ".")) {
        stop("Diagonal of structure should be '.'")
      }
      return(TRUE)
    },
    
    tpm2par = function(tpm) {
      ltpm <- log(tpm / diag(tpm))
      ltpm <- t(ltpm) # transpose to fill by rows (like in C++)
      par <- ltpm[!diag(self$nstates())]
      return(par)
    },
    
    par2tpm = function(par) {
      tpm <- diag(self$nstates())
      tpm[!diag(self$nstates())] <- exp(par)
      tpm <- t(tpm) # transpose to fill by rows (like in C++)
      tpm <- tpm / rowSums(tpm)
      return(tpm)
    }
    
  )
)





