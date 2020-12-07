
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
    #' @param coeff_fe0 Initial coefficients for fixed effects parameters. If
    #' not provided, the intercepts are set using \code{tpm0}, and other 
    #' coefficients are set to 0.
    #' @param coeff_re0 Initial coefficients for random effects parameters.
    #' Defaults to 0 if not provided. 
    #' @param data Data frame, needed if the model includes covariates
    #' 
    #' @return A new MarkovChain object
    initialize = function(n_states = NULL, structure = NULL, 
                          tpm0 = NULL, coeff_fe0 = NULL, 
                          coeff_re0 = NULL, data = NULL) {
      # Check arguments
      private$check_args(n_states = n_states, structure = structure,
                         tpm0 = tpm0, coeff_fe0 = coeff_fe0, 
                         coeff_re0 = coeff_re0, data = data)
      
      # Define 'structure' as matrix
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
      
      # Create list of formulas  ('structure' is transposed to get 
      # the formulas in the order 1>2, 1>3, ..., 2>1, 2>3, ...)
      ls_form_char <- as.list(t(structure)[!diag(self$nstates())])
      ls_form <- lapply(ls_form_char, function(form_char) {
        if(form_char == ".")
          return(NULL)
        else
          return(as.formula(form_char))
      })
      # Names for transition probabilities
      tr_names <- paste0("S", rep(1:n_states, each = n_states), 
                         ">S", rep(1:n_states, n_states))
      names(ls_form) <- tr_names[-which(diag(n_states) == 1)]
      
      # Set structure and formulas attributes
      private$structure_ <- structure
      private$formulas_ <- ls_form
      
      # Does the hidden state model include covariates?
      no_covs <- all(structure %in% c(".", "~1"))
      
      # Get structure of design matrices
      if(no_covs) {
        # Create temporary dummy data set to pass to make_mat
        data <- data.frame(dummy = rep(1, 2))
      } else if(is.null(data)) {
        stop("'data' must be provided if the model includes covariates")
      }
      
      mats <- self$make_mat(data = data)
      ncol_fe <- mats$ncol_fe
      ncol_re <- mats$ncol_re       
      private$terms_ <- list(ncol_fe = ncol_fe,
                             ncol_re = ncol_re,
                             names_fe = colnames(mats$X_fe),
                             names_re_all = colnames(mats$X_re),
                             names_re = names(ncol_re))
      
      # Initialise coeff_fe and coeff_re to 0
      self$update_coeff_fe(rep(0, sum(ncol_fe)))
      self$update_coeff_re(rep(0, sum(ncol_re))) 
      self$update_lambda(rep(1, length(ncol_re)))
      self$update_delta(rep(1 / n_states, n_states))
      
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
    
    #' @description Current delta parameter estimates 
    delta = function() {return(private$delta_)}, 
    
    #' @description Current parameter estimates (random effects)
    coeff_re = function() {return(private$coeff_re_)},
    
    #' @description Smoothness parameters
    lambda = function() {return(private$lambda_)},
    
    #' @description Variance components of smooth terms
    #' 
    #' This function transforms the smoothness parameter of
    #' each smooth term into a standard deviation, given by 
    #' SD = 1/sqrt(lambda). It is particularly helpful to get the
    #' standard deviations of independent normal random effects.
    vcomp = function() {return(1/sqrt(private$lambda_))},
    
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
      private$coeff_fe_ <- matrix(rep(0, sum(self$terms()$ncol_fe)))
      private$coeff_fe_[i0] <- self$tpm2par(tpm)
      rownames(private$coeff_fe_) <- self$terms()$names_fe
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
      private$coeff_fe_ <- matrix(coeff_fe)
      rownames(private$coeff_fe_) <- self$terms()$names_fe
      if(all(self$structure() %in% c(".", "~1"))) {
        # Only update tpm if no covariates
        private$tpm_ <- self$par2tpm(coeff_fe)        
      }
    },
    
    #' @description Update coefficients for random effect parameters
    #' 
    #' @param coeff_re Vector of coefficients for random effect parameters
    update_coeff_re = function(coeff_re) {
      private$coeff_re_ <- matrix(coeff_re)
      rownames(private$coeff_re_) <- self$terms()$names_re_all
    },
    
    #' @description Update delta coefficients 
    update_delta = function(new_delta) {
      private$delta_ <- new_delta
      names(private$delta_) <- paste0("state", 1:length(new_delta))
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
    #'   \item{ncol_fe}{Number of columns of X_fe for each parameter}
    #'   \item{ncol_re}{Number of columns of X_re and S for each random effect}
    #' }
    make_mat = function(data, new_data = NULL) {
      make_matrices(formulas = self$formulas(), 
                    data = data, 
                    new_data = new_data)
    },
    
    #' @description Design matrices for grid of covariates
    #' 
    #' Used in plotting functions such as HMM$plot_tpm and HMM$plot_stat_dist
    #' 
    #' @param var Name of variable
    #' @param data Data frame containing the covariates
    #' @param covs Optional data frame with a single row and one column
    #' for each covariate, giving the values that should be used. If this is
    #' not specified, the mean value is used for numeric variables, and the
    #' first level for factor variables.
    #' @param n_grid Grid size (number of points). Default: 1000.
    #' 
    #' @return A list with the same elements as the output of make_mat, 
    #' plus a data frame of covariates values.
    make_mat_grid = function(var, data, covs = NULL, n_grid = 1e3) {
      # Data frame for covariate grid
      new_data <- cov_grid(var = var, data = data, covs = covs, 
                           formulas = self$formulas(), n_grid = n_grid)
      
      # Create design matrices
      mats <- self$make_mat(data = data, new_data = new_data)
      
      # Save data frame of covariate values
      mats$new_data <- new_data
      
      return(mats)
    },
    
    #' @description Transform transition probabilities to working scale
    #' 
    #' Apply the multinomial logit link function to get the corresponding parameters on the
    #' working scale (i.e., linear predictor scale).
    #' 
    #' @param tpm Transition probability matrix
    #' 
    #' @return Vector of parameters on linear predictor scale
    tpm2par = function(tpm) {
      ltpm <- log(tpm / diag(tpm))
      ltpm <- t(ltpm) # transpose to fill by rows (like in C++)
      par <- ltpm[!diag(self$nstates())]
      return(par)
    },
    
    #' @description Transform working parameters to transition probabilities
    #' 
    #' Apply the inverse multinomial logit link function to transform the parameters on
    #' the working scale (i.e., linear predictor scale) into the transition
    #' probabilities.
    #' 
    #' @param par Vector of parameters on working scale
    #' 
    #' @return Transition probability matrix
    par2tpm = function(par) {
      tpm <- diag(self$nstates())
      tpm[!diag(self$nstates())] <- exp(par)
      tpm <- t(tpm) # transpose to fill by rows (like in C++)
      tpm <- tpm / rowSums(tpm)
      return(tpm)
    },
    
    #' @description Get transition probability matrices from design matrices
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
    #' 
    #' @return Array with one slice for each transition probability matrix
    tpm_all = function(X_fe, X_re, coeff_fe = NULL, coeff_re = NULL) {
      n_states <- self$nstates()
      
      # Define parameters
      if(length(coeff_fe) == 0)
        coeff_fe <- self$coeff_fe()
      if(length(coeff_re) == 0)
        coeff_re <- self$coeff_re()
      
      # Linear predictor
      ltpm <- X_fe %*% coeff_fe + X_re %*% coeff_re
      ltpm_mat <- matrix(ltpm, ncol = n_states * (n_states - 1))
      
      # Transition probability matrices
      tpm <- apply(ltpm_mat, 1, self$par2tpm)
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
      if(is.null(new_data$ID)) {
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
      cat("#########################\n")
      cat("## State process model ##\n")
      cat("#########################\n")
      # Data frame of formulas on transition probabilities
      hid_forms <- as.data.frame(self$structure())
      n_states <- self$nstates()
      rownames(hid_forms) <- paste0("state ", 1:n_states)
      colnames(hid_forms) <- paste0("state ", 1:n_states)
      print(hid_forms)
      cat("\n")
    }
  ),
  
  private = list(
    ################
    ## Attributes ##
    ################
    structure_ = NULL,
    formulas_ = NULL,
    coeff_fe_ = NULL,
    coeff_re_ = NULL,
    delta_ = NULL, 
    lambda_ = NULL,
    tpm_ = NULL,
    nstates_ = NULL,
    terms_ = NULL,
    
    #################################
    ## Check constructor arguments ##
    #################################
    # (For argument description, see constructor)
    check_args = function(n_states, structure, tpm0, coeff_fe0, coeff_re0, data) {
      if(!is.null(n_states)) {
        if(!is.numeric(n_states) | n_states < 1) {
          stop("'n_states' should be a numeric >= 1")
        }        
      }
      
      if(!is.null(structure)) {
        if(is.matrix(structure)) {
          if(nrow(structure) != ncol(structure)) {
            stop("'structure' should be a square matrix")
          }
          
          if(!all(is.character(as.vector(structure)))) {
            stop("'structure' should be a matrix of character strings")
          }
          
          if (!all(diag(structure) == ".")) {
            stop("Diagonal of structure should be '.'")
          }
          
        } else if(!inherits(structure, "formula")) {
          stop("'structure' should be either a matrix or a formula")
        }
      }
      
      if(!is.null(tpm0)) {
        if(!is.matrix(tpm0)) {
          stop("'tpm0' should be a matrix")
        }
        
        if(nrow(tpm0) != ncol(tpm0)) {
          stop("'tpm0' should be a square matrix")
        }
        
        if(any(rowSums(tpm0) != 1)) {
          stop("The rows of 'tpm0' should sum to 1")
        }
      }
      
      if(!is.null(coeff_fe0)) {
        if(!is.numeric(coeff_fe0) | !is.vector(coeff_fe0)) {
          stop("'coeff_fe0' should be a numeric vector")
        }
      }
      
      if(!is.null(coeff_re0)) {
        if(!is.numeric(coeff_re0) | !is.vector(coeff_re0)) {
          stop("'coeff_re0' should be a numeric vector")
        }
      }
      
      if(!is.null(data)) {
        if(!inherits(data, "data.frame")) {
          stop("'data' should be a data.frame")
        }
      }
    }
  )
)
