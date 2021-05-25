
#' R6 class for HMM hidden process model
#'
#' Contains the parameters and model formulas for the hidden process model.
MarkovChain <- R6Class(
  classname = "MarkovChain",
  
  public = list(

    # Constructor -------------------------------------------------------------
    
    #' @description Create new MarkovChain object
    #' 
    #' @param data Data frame, needed to create model matrices
    #' @param structure Either (1) matrix with an entry of "." on diagonal, a "0" for 
    #' transitions that are not allowed (not implemented yet), and a formula "~1" 
    #' for covariates affecting transitions that are to be estimated, or (2) single
    #' formula, assumed for all transition probabilities. (Default: no covariate
    #' dependence.)
    #' @param n_states Number of states. If not specified, then \code{structure} 
    #' needs to be provided as a matrix, and n_states is deduced from its dimensions.
    #' @param tpm Optional transition probability matrix, to initialise the model
    #' parameters (intercepts in model with covariates). If not provided, the default 
    #' is a matrix with 0.9 on the diagonal. 
    #' @param stationary if TRUE then stationary distribution with respect to tpm for 
    #' first time point is used as initial distribution, if FALSE then initial distribution
    #' is estimated 
    #' 
    #' @return A new MarkovChain object
    initialize = function(data,
                          structure = NULL, 
                          n_states = NULL,
                          tpm = NULL,
                          stationary = FALSE) {
      # Check arguments
      private$check_args(n_states = n_states, 
                         structure = structure, 
                         stationary = stationary, 
                         data = data)
      
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
      
      # set whether delta to be stationary or not
      private$stationary_ <- stationary 
      
      # Create list of formulas  ('structure' is transposed to get 
      # the formulas in the order 1>2, 1>3, ..., 2>1, 2>3, ...)
      ls_form_char <- as.list(t(structure)[!diag(self$nstates())])
      ls_form <- lapply(ls_form_char, function(form_char) {
        if(form_char == ".")
          return(as.formula("~1"))
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
      
      # Create terms and necessary matrices 
      mats <- self$make_mat(data = data)
      ncol_fe <- mats$ncol_fe
      ncol_re <- mats$ncol_re       
      private$terms_ <- c(mats, list(names_fe = colnames(mats$X_fe),
                                     names_re_all = colnames(mats$X_re),
                                     names_re = names(ncol_re)))
      
      # Initialise coeff_fe and coeff_re to 0
      self$update_coeff_fe(rep(0, sum(ncol_fe)))
      self$update_coeff_re(rep(0, ncol(mats$X_re)))
      self$update_lambda(rep(1, ifelse(is.null(ncol_re), 0, ncol(ncol_re))))
      self$update_delta(rep(1 / n_states, n_states))
      
      # Initialise tpm
      if(is.null(tpm)) {
        # Defaults to diagonal of 0.9 if no initial tpm provided
        tpm <- matrix(0.1/(n_states-1), nrow = n_states, ncol = n_states)
        diag(tpm) <- 0.9
      }
      self$update_tpm(tpm)
    },
    

    # Accessors ---------------------------------------------------------------

    #' @description Structure of MarkovChain model
    structure = function() {return(private$structure_)},
    
    #' @description List of formulas for MarkovChain model
    formulas = function() {return(private$formulas_)},
    
    #' @description Get transition probability matrices 
    #' 
    #' @param t time point, default 1; if t = "all" then 
    #' all tpms are returned otherwise tpms for time points in t are returned
    #' @param linpred custom linear predictor 
    #' 
    #' @return Array with one slice for each transition probability matrix
    tpm = function(t = 1, linpred = NULL) {
      n_states <- self$nstates()
      npar <- n_states * (n_states - 1)
      if (is.null(linpred)) linpred <- self$linpred() 
      T <- length(linpred) / npar
      if (length(t) == 1) if (t == "all") t <- 1:T
      ind <- as.vector(sapply(1:npar, function(i) {t + (i - 1) * T}))
      linpred <- matrix(linpred[ind], ncol = n_states * (n_states - 1))
      val <- apply(linpred, 1, self$par2tpm)
      val <- array(val, dim = c(n_states, n_states, ncol(val)))
      return(val)
    },
    
    #' @description Current parameter estimates (fixed effects)
    coeff_fe = function() {return(private$coeff_fe_)},
    
    #' @description Stationary distributions
    #' @param t time point, default is 1; if t = "all" then 
    #' all deltas are returned otherwise deltas for time points in t are returned 
    #' @param linpred custom linear predictor 
    #' @return Matrix of stationary distributions. Each row corresponds to
    #' a row of the design matrices, and each column corresponds to a state.
    delta = function(t = NULL, linpred = NULL) {
      if (!is.null(t)) {
        # Number of states
        n_states <- self$nstates()
        
        # Derive transition probability matrices
        tpms <- self$tpm(t = t, linpred = linpred)
        
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
      } else {
        stat_dists <- private$delta_
      }
      return(stat_dists)
    }, 
    
    #' @description Use stationary distribution as initial distribution? 
    stationary = function() {return(private$stationary_)}, 
    
    #' @description Current parameter estimates (random effects)
    coeff_re = function() {return(private$coeff_re_)},
    
    #' @description Fixed effect design matrix 
    X_fe = function() {return(private$terms_$X_fe)}, 
    
    #' @description Random effect design matrix 
    X_re = function() {return(private$terms_$X_re)}, 
    
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
    

    # Mutators ----------------------------------------------------------------

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
    },
    
    #' @description Update coefficients for random effect parameters
    #' 
    #' @param coeff_re Vector of coefficients for random effect parameters
    update_coeff_re = function(coeff_re) {
      private$coeff_re_ <- matrix(coeff_re)
      rownames(private$coeff_re_) <- self$terms()$names_re_all
    },
    
    #' @description Update design matrix for fixed effects 
    #' 
    #' @param X_fe new design matrix for fixed effects
    update_X_fe = function(X_fe) {
      private$terms_$X_fe <- X_fe
    }, 
    
    #' @description Update design matrix for random effects 
    #' 
    #' @param X_re new design matrix for random effects
    update_X_re = function(X_re) {
      private$terms_$X_re <- X_re
    }, 
    
    #' @description Update delta coefficients 
    #' 
    #' @param delta Vector of delta coefficients
    update_delta = function(delta) {
      private$delta_ <- delta
      names(private$delta_) <- paste0("state", 1:length(delta))
    }, 
    
    #' @description Update smoothness parameters
    #' 
    #' @param lambda New smoothness parameter vector
    update_lambda = function(lambda) {
      private$lambda_ <- matrix(lambda)
      rownames(private$lambda_) <- self$terms()$names_re
    },
    

    # Other methods -----------------------------------------------------------

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
    #' @param covs Optional named list for values of covariates (other than 'var') 
    #' that should be used in the plot (or dataframe with single row). If this is
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
    
    #' @description Linear predictor for transition probabilities
    linpred = function() {
      linpred <- self$X_fe() %*% self$coeff_fe() + self$X_re() %*% self$coeff_re()
      return(linpred[,1])
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
    #' @param silent if TRUE then no messages are printed 
    #' 
    #' @return Sequence of states of simulated chain
    simulate = function(n, data = NULL, new_data = NULL, silent = FALSE) {
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
      X_fe_old <- mats_hid$X_fe
      X_re_old <- mats_hid$X_re
      self$update_X_fe(mats_hid$X_fe)
      self$update_X_re(mats_hid$X_re)
      tpms <- self$tpm(t = "all")
      
      # Uniform initial distribution for now
      delta <- rep(1/n_states, n_states) 
      
      # Simulate state process      
      S <- rep(NA, n)
      S[1] <- sample(1:n_states, size = 1, prob = delta)
      for(i in 2:n) {
        if(!silent & round(i/n*100)%%10 == 0) {
          cat("\rSimulating states... ", round(i/n*100), "%", sep = "")        
        }
        
        if(ID[i] != ID[i-1]) {
          S[i] <- sample(1:n_states, size = 1, prob = delta)
        } else {
          S[i] <- sample(1:n_states, size = 1, prob = tpms[S[i-1], , i-1])          
        }
      }
      if(!silent) cat("\n")
      
      # reset design matrices
      self$update_X_fe(X_fe_old)
      self$update_X_re(X_re_old)
      
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
    },
    
    #' @description Print MarkovChain object
    print = function() {
      self$formulation()
    }
  ),
  
  private = list(

    # Private data members ----------------------------------------------------
    
    structure_ = NULL,
    stationary_ = NULL, 
    formulas_ = NULL,
    coeff_fe_ = NULL,
    coeff_re_ = NULL,
    delta_ = NULL, 
    lambda_ = NULL,
    nstates_ = NULL,
    terms_ = NULL,
    
    #' Check constructor arguments ##
    # (For argument description, see constructor)
    check_args = function(n_states, structure, stationary, data) {
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
      
      if(!is.null(data)) {
        if(!inherits(data, "data.frame")) {
          stop("'data' should be a data.frame")
        }
      }
    }
  )
)
