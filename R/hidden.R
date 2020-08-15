
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
    #' @param tpm0 Initial transition probability matrix. (Default: 0.9 on diagonal,
    #' and 0.1/(n_states - 1) for all other entries.) If the model has covariates,
    #' then \code{tpm0} is used to set the intercept parameters for the transition
    #' probabilities, and the other parameters are set to 0.
    #' @param coeff_fe0 Initial parameters for fixed effects. 
    #' @param data HmmData object, needed if the model includes covariates
    #' 
    #' @return A new MarkovChain object
    initialize = function(n_states = NULL, structure = NULL, 
                          tpm0 = NULL, coeff_fe0 = NULL, 
                          coeff_re0 = NULL, data = NULL) {
      if(is.null(structure)) {
        # No covariate effects
        structure <- matrix("~1", nrow = n_states, ncol = n_states)
        diag(structure) <- "."
        private$nstates_ <- n_states
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
      }
      
      # Create list of formulas
      ls_form_char <- as.list(structure[!diag(self$nstates())])
      ls_form <- lapply(ls_form_char, function(form_char) {
        if(form_char == ".")
          return(NULL)
        else
          return(as.formula(form_char))
      })
      
      # Set structure and formulas attributes
      private$structure_ <- structure
      private$formulas_ <- ls_form
      
      # Does the hidden state model include covariates?
      no_covs <- all(structure %in% c(".", "~1"))

      # Get structure of design matrices      
      if(no_covs) {
        # If no covariates, N*(N-1) fixed effects and 0 random effects
        ncol_fe <- rep(1, n_states * (n_states - 1))
        ncol_re <- 0
      } else {
        if(is.null(data)) {
          stop("'data' must be provided if the model includes covariates")
        }
        mats <- self$make_mat(data = data$data())
        ncol_fe <- mats$ncol_fe
        ncol_re <- mats$ncol_re        
      }
      private$terms_ <- list(ncol_fe = ncol_fe,
                             ncol_re = ncol_re)
      
      # Initialise coeff_fe and coeff_re to 0
      self$update_coeff_fe(rep(0, sum(ncol_fe)))
      self$update_coeff_re(rep(0, sum(ncol_re))) 
      
      # Set fixed effect parameters, using either coeff_fe0 or tpm0
      if(!is.null(coeff_fe0)) {
        self$update_coeff_fe(coeff_fe0)
      } else {
        if(is.null(tpm0)) {
          # Defaults to diagonal of 0.9 if no initial tpm provided
          tpm0 <- matrix(0.1/(n_states-1), nrow = n_states, ncol = n_states)
          diag(tpm0) <- 0.9
        }
        self$update_tpm(tpm0)
      }
      
      # Set random effect parameters
      if(!is.null(coeff_re0)) {
        self$update_coeff_re(coeff_re0)
      }
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
    coeff_fe = function() {return(private$coeff_fe_)},
    
    #' @description Current parameter estimates (random effects)
    coeff_re = function() {return(private$coeff_re_)},
    
    #' @description Number of states
    nstates = function() {return(private$nstates_)},
    
    #' @description Terms of model formulas
    terms = function() {return(private$terms_)},
    
    ##############
    ## Mutators ##
    ##############
    #' @description Update transition probability matrix
    #' 
    #' @param tpm New transition probability matrix
    update_tpm = function(tpm) {
      n_states <- self$nstates()
      
      if(!is.matrix(tpm)) {
        stop("'tpm' should be a matrix")
      } else if(nrow(tpm) != n_states | ncol(tpm) != n_states) {
        stop("'tpm' should have ", n_states, " rows and ", 
             n_states, " columns")
      } else if(any(rowSums(tpm) != 1)) {
        stop("The rows of 'tpm' should sum to 1")
      }
      
      # Indices of intercepts in coeff_fe
      ncol_fe <- self$terms()$ncol_fe
      n_par <- length(ncol_fe)
      i0 <- c(1, cumsum(ncol_fe)[-n_par] + 1)
      
      # Update coeff_fe and tpm attributes
      private$coeff_fe_[i0] <- private$tpm2par(tpm)
      private$tpm_ <- tpm
    },
    
    #' @description Update coefficients for fixed effect parameters
    #' 
    #' @param coeff_fe Vector of coefficients for fixed effect parameters
    update_coeff_fe = function(coeff_fe) {
      ncol_total <- sum(self$terms()$ncol_fe)
      if(length(coeff_fe) != ncol_total) {
        stop("'coeff_fe' should be of length ", ncol_total, " (one ",
             "parameter for each column of the design matrix)")
      }
      private$coeff_fe_ <- coeff_fe
      if(all(self$structure() %in% c(".", "~1"))) {
        # Only update tpm if no covariates
        private$tpm_ <- private$par2tpm(coeff_fe)        
      }
    },
    
    #' @description Update coefficients for random effect parameters
    #' 
    #' @param coeff_re Vector of coefficients for random effect parameters
    update_coeff_re = function(coeff_re) {
      private$coeff_re_ <- newpar
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
      ltpm <- X_fe %*% self$coeff_fe() + X_re %*% self$coeff_re()
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
    },
    
    #' @description Print model formulation
    formulation = function() {
      cat("## State process model:\n")
      # Data frame of formulas on transition probabilities
      hid_forms <- as.data.frame(self$structure())
      n_states <- self$nstates()
      rownames(hid_forms) <- paste0("state ", 1:n_states)
      colnames(hid_forms) <- paste0("state ", 1:n_states)
      print(hid_forms)
    }
  ),
  
  private = list(
    structure_ = NULL,
    formulas_ = NULL,
    coeff_fe_ = NULL,
    coeff_re_ = NULL,
    tpm_ = NULL,
    nstates_ = NULL,
    terms_ = NULL,
    
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





