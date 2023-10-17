
#' @title R6 class for HMM hidden process model
#'
#' @description
#' Contains the parameters and model formulas for the hidden process model.
#' 
#' @export
MarkovChain <- R6Class(
  classname = "MarkovChain",
  
  public = list(
    
    # Constructor -------------------------------------------------------------
    
    #' @description Create new MarkovChain object
    #' 
    #' @param data Data frame, needed to create model matrices, and to identify the
    #' number of time series (which each have a separate initial distribution)
    #' @param formula Either (1) R formula, used for all transition probabilities, 
    #' or (2) matrix of character strings giving the formula for each transition
    #' probability, with "." along the diagonal (or for reference elements; see
    #' \code{ref} argument). (Default: no covariate dependence.)
    #' @param n_states Number of states. If not specified, then \code{formula} 
    #' needs to be provided as a matrix, and n_states is deduced from its dimensions.
    #' @param tpm Optional transition probability matrix, to initialise the model
    #' parameters (intercepts in model with covariates). If not provided, the default 
    #' is a matrix with 0.9 on the diagonal. 
    #' @param initial_state Specify model for initial state distribution. There
    #' are five different options: 
    #' \itemize{
    #'    \item "estimated": a separate initial distribution is estimated for 
    #'    each ID (default)
    #'    \item "stationary": the initial distribution is fixed to the
    #'    stationary distribution of the transition probability matrix for the 
    #'    first time point of each ID
    #'    \item "shared": a common initial distribution is estimated for all IDs
    #'    \item integer value between 1 and n_states: used as the known initial 
    #'    state for all IDs
    #'    \item vector of integers between 1 and n_states (of length the number 
    #'    of IDs): each element is used as the known initial state for the 
    #'    corresponding ID
    #' }
    #' @param fixpar List with optional elements "hid" (fixed parameters for
    #' transition probabilities), "lambda_hid" (fixed smoothness parameters), 
    #' and "delta0" (fixed parameters for initial distribution). Each element is 
    #' a named vector of coefficients that should either be fixed (if the 
    #' corresponding element is set to NA) or estimated to a common value (using 
    #' integers or factor levels). 
    #' @param ref Vector of indices for reference transition probabilities, 
    #' of length \code{n_states}. The i-th element is the index for the 
    #' reference in the i-th row of the transition probability matrix. For 
    #' example, ref = c(1, 1) means that the first element of the first row
    #' Pr(1>1) and the first element of the second row Pr(2>1) are used as 
    #' reference elements and are not estimated. If this is not provided,
    #' the diagonal transition probabilities are used as references. 
    #' 
    #' @return A new MarkovChain object
    #' 
    #' @examples
    #' # Load data set from MSwM package
    #' data(energy, package = "MSwM")
    #' 
    #' # Create 2-state covariate-free model and initialise transition 
    #' # probability matrix
    #' hid <- MarkovChain$new(data = energy, n_states = 2,
    #'                        tpm = matrix(c(0.8, 0.3, 0.2, 0.7), 2, 2))
    #' 
    #' # Create 2-state model with non-linear effect of Oil on all transition 
    #' # probabilities
    #' hid <- MarkovChain$new(data = energy, n_states = 2,
    #'                        formula = ~ s(Oil, k = 5, bs = "cs"))
    #' 
    #' # Create 2-state model with quadratic effect of Oil on Pr(1 > 2)
    #' structure <- matrix(c(".", "~poly(Oil, 2)",
    #'                       "~1", "."),
    #'                     ncol = 2, byrow = TRUE)
    #' hid <- MarkovChain$new(data = energy, n_states = 2,
    #'                        formula = structure)
    initialize = function(data,
                          formula = NULL, 
                          n_states,
                          tpm = NULL,
                          initial_state = "estimated",
                          fixpar = NULL,
                          ref = 1:n_states) {
      # Check arguments
      private$check_args(n_states = n_states, 
                         formula = formula, 
                         data = data)
      private$nstates_ <- n_states
      private$ref_ <- ref
      
      # Matrix with 1 for reference element and 0 elsewhere
      # (used later to select relevant entries of tpm)
      private$ref_mat_ <- matrix(0, n_states, n_states)
      private$ref_mat_[cbind(1:n_states, ref)] <- 1
      
      # Define 'formula' as matrix
      if(is.null(formula)) {
        # No covariate effects
        formula <- matrix("~1", nrow = n_states, ncol = n_states)
      } else if(length(formula) == 1 | inherits(formula, "formula")) {
        # Same formula for all transitions
        # Use deparse to transform to string without linebreaks
        formula <- matrix(deparse(formula, width.cutoff = 500), 
                          nrow = n_states, 
                          ncol = n_states)
      } else {
        if(!all(formula[cbind(1:n_states, ref)] == ".")) {
          stop("Formulas for reference transition probabilities should be '.'.")
        }
      }
      formula[cbind(1:n_states, ref)] <- "."
      
      # How many time series are in data (used to initialise delta0)
      if(! "ID" %in% colnames(data)) {
        # Default = 1 if ID not provided
        private$unique_ID_ <- 1
      } else {
        private$unique_ID_ <- unique(data$ID)
      }
      
      # Should initial distribution delta0 be stationary?
      private$stationary_ <- FALSE
      if(length(initial_state) == 1) {
        if(initial_state == "stationary") {
          private$stationary_ <- TRUE          
        }
      }
      
      # Should initial distribution be estimated? If there are many IDs,
      # ask user to check that this is really the best option.
      if(length(initial_state) == 1) {
        if(initial_state == "estimated") {
          if(length(self$unique_ID()) >= 8) {
            message(paste("Large number of time series in the data.",
                          "Consider changing the 'initial_state' argument",
                          "to avoid estimating a different initial state",
                          "distribution for each time series."))
          }
        }
      }
      
      # Save fixpar for use in HMM
      private$fixpar_user_ <- fixpar
      private$fixpar_ <- fixpar
      
      # Create list of formulas  ('formula' is transposed to get 
      # the formulas in the order 1>2, 1>3, ..., 2>1, 2>3, ...)
      ls_form_char <- as.list(t(formula)[!t(self$ref_mat())])
      ls_form <- lapply(ls_form_char, function(form_char) {
        if(form_char == ".")
          return(as.formula("~1"))
        else
          return(as.formula(form_char))
      })
      
      # Names for transition probabilities
      tr_names <- paste0("S", rep(1:n_states, each = n_states), 
                         ">S", rep(1:n_states, n_states))
      names(ls_form) <- tr_names[-which(t(self$ref_mat()) == 1)]
      
      # Set formula and formulas attributes
      private$formula_ <- formula
      private$formulas_ <- ls_form
      
      # Get names of all covariates
      var_names <- unique(rapply(self$formulas(), all.vars))
      # Remove pi from list of covariates if it is in the formulas
      var_names <- var_names[which(var_names!="pi")]
      if(length(var_names) > 0) {
        # Remove NAs in covariates (replace by last non-NA value)
        data[,var_names] <- lapply(data[,var_names, drop=FALSE], 
                                   function(col) na_fill(col))
      }
      
      # Create terms and necessary matrices 
      mats <- self$make_mat(data = data)
      ncol_fe <- mats$ncol_fe
      ncol_re <- mats$ncol_re       
      private$terms_ <- c(mats, list(names_fe = colnames(mats$X_fe),
                                     names_re_all = colnames(mats$X_re),
                                     names_re = colnames(ncol_re)))
      
      # Initialise coeff_fe and coeff_re to 0
      self$update_coeff_fe(rep(0, sum(ncol_fe)))
      self$update_coeff_re(rep(0, ncol(mats$X_re)))
      self$update_lambda(rep(1, ifelse(is.null(ncol_re), 0, ncol(ncol_re))))
      
      # Setup initial distribution
      private$initial_state_ <- initial_state
      delta0 <- matrix(1/self$nstates(), nrow = length(self$unique_ID()), 
                       ncol = self$nstates())
      private$ref_delta0_ <- rep(ncol(delta0), length = nrow(delta0))
      if(is.numeric(initial_state)) {
        # Initial state is known (possible one for each ID)
        delta0[] <- 0
        delta0[cbind(1:nrow(delta0), initial_state)] <- 1
        private$ref_delta0_ <- rep(initial_state, length = nrow(delta0))
      }
      self$update_delta0(delta0 = delta0)
      
      # Initialise tpm
      if(is.null(tpm)) {
        # Defaults to diagonal of 0.9 if no initial tpm provided
        tpm <- matrix(0.1/(n_states-1), nrow = n_states, ncol = n_states)
        diag(tpm) <- 0.9
      }
      self$update_tpm(tpm)
      
      # Store information about fixed parameters
      self$update_fixpar(fixpar = fixpar)
    },
    
    
    # Accessors ---------------------------------------------------------------
    
    #' @description Formula of MarkovChain model
    formula = function() {return(private$formula_)},
    
    #' @description List of formulas for MarkovChain model
    formulas = function() {return(private$formulas_)},
    
    #' @description Get transition probability matrices 
    #' 
    #' @param t Time index or vector of time indices; default = 1. If t = "all" 
    #' then all transition probability matrices are returned.
    #' @param linpred Optional custom linear predictor 
    #' 
    #' @return Array with one slice for each transition probability matrix
    tpm = function(t = 1, linpred = NULL) {
      n_states <- self$nstates()
      npar <- n_states * (n_states - 1)
      if (is.null(linpred)) {
        linpred <- self$linpred() 
      }
      n <- length(linpred) / npar
      if (length(t) == 1) {
        if (t == "all") {
          t <- 1:n
        }
      }
      
      # Get the right entries of linpred for choice of t
      ind <- as.vector(sapply(1:npar, function(i) {t + (i - 1) * n}))
      linpred <- matrix(linpred[ind], ncol = n_states * (n_states - 1))
      
      # Get transition probs from linear predictor
      val <- apply(linpred, 1, self$par2tpm)
      val <- array(val, dim = c(n_states, n_states, ncol(val)))
      rownames(val) <- paste0("state ", 1:n_states)
      colnames(val) <- paste0("state ", 1:n_states)
      return(val)
    },
    
    #' @description Indices of reference elements in transition probability
    #' matrix
    ref = function() {return(private$ref_)},
    
    #' @description Matrix of reference elements in transition probability
    #' matrix
    ref_mat = function() {return(private$ref_mat_)},
    
    #' @description Indices of reference elements in initial distribution
    ref_delta0 = function() {return(private$ref_delta0_)},
    
    #' @description Current parameter estimates (fixed effects)
    coeff_fe = function() {return(private$coeff_fe_)},
    
    #' @description Stationary distribution
    #' 
    #' @param t Time point(s) for which stationary distribution should be returned. 
    #' If t = "all", all deltas are returned; else this should be a vector of
    #' time indices. If NULL (default), the stationary distribution for the first
    #' time step is returned. 
    #' @param linpred Optional custom linear predictor 
    #' 
    #' @return Matrix of stationary distributions. Each row corresponds to
    #' a row of the design matrices, and each column corresponds to a state.
    delta = function(t = NULL, linpred = NULL) {
      if(is.null(t)) {
        if(is.null(linpred)) {
          t <- 1
        } else {
          t <- "all"
        }
      }
      
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
      
      return(stat_dists)
    }, 
    
    #' @description Initial distribution
    #' 
    #' @param log Logical indicating whether to return the log of the initial
    #' probabilities (default: FALSE). If TRUE, then the last element is
    #' excluded, as it is not estimated.
    #' @param as_matrix Logical indicating whether the output should be
    #' formatted as a matrix (default). If as_matrix is FALSE and log is
    #' TRUE, the result is formatted as a column vector.
    #' 
    #' @return Matrix with one row for each time series ID, and one column
    #' for each state. For each ID, the i-th element of the corresponding 
    #' row is the probability Pr(S[1] = i)
    delta0 = function(log = FALSE, as_matrix = TRUE) {
      d0 <- private$delta0_
      if(!log) {
        return(d0)
      } else {
        ref <- self$ref_delta0()
        num <- matrix(sapply(1:nrow(d0), function(i) d0[i, -ref[i]]), 
                      nrow = nrow(d0), byrow = TRUE)
        denom <- sapply(1:nrow(d0), function(i) d0[i, ref[i]])
        log_d0 <- log(num / denom)
        
        if(!as_matrix) {
          states <- matrix(1:ncol(d0), nrow = nrow(d0), 
                           ncol = ncol(d0), byrow = TRUE)
          states <- as.vector(
            sapply(1:nrow(states), function(i) states[i, -ref[i]]))
          d0_names <- paste0("ID:", rep(1:nrow(log_d0), each = ncol(log_d0)), 
                             ".state", states)
          log_d0 <- matrix(log_d0, ncol = 1)
          rownames(log_d0) <- d0_names
        }
        
        return(log_d0)
      }
    },
    
    #' @description Use stationary distribution as initial distribution? 
    stationary = function() {return(private$stationary_)}, 
    
    #' @description Fixed parameters
    #' 
    #' @param all Logical. If FALSE, only user-specified fixed
    #' parameters are returned, but not parameters that are fixed
    #' for some other reason (e.g., from '.' in formula)
    fixpar = function(all = FALSE) {
      if(all) {
        return(private$fixpar_)
      } else {
        return(private$fixpar_user_)        
      }
    },
    
    #' @description Current parameter estimates (random effects)
    coeff_re = function() {return(private$coeff_re_)},
    
    #' @description Fixed effect design matrix 
    X_fe = function() {return(private$terms_$X_fe)}, 
    
    #' @description Random effect design matrix 
    X_re = function() {return(private$terms_$X_re)}, 
    
    #' @description Smoothness parameters
    lambda = function() {return(private$lambda_)},
    
    #' @description Standard deviation of smooth terms
    #' 
    #' This function transforms the smoothness parameter of
    #' each smooth term into a standard deviation, given by 
    #' SD = 1/sqrt(lambda). It is particularly helpful to get the
    #' standard deviations of independent normal random effects.
    sd_re = function() {return(1/sqrt(private$lambda_))},
    
    #' @description Number of states
    nstates = function() {return(private$nstates_)},
    
    #' @description Terms of model formulas
    terms = function() {return(private$terms_)},
    
    #' @description Number of time series
    unique_ID = function() {return(private$unique_ID_)},
    
    #' @description Initial state (see constructor argument)
    initial_state = function() {return(private$initial_state_)},
    
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
    
    #' @description Update initial distribution 
    #' 
    #' @param delta0 Either a matrix where the i-th row is the initial
    #' distribution for the i-th time series in the data, or a vector which is
    #' then used for all time series. Entries of each row of delta0 should sum
    #' to one.
    update_delta0 = function(delta0) {
      
      # If input is vector, copy into matrix
      if(is.null(dim(delta0))) {
        if(length(delta0) != self$nstates()) {
          stop(paste0("'delta0' should have length the number of states (", 
                      self$nstates(), ")"))
        }
        delta0 <- matrix(delta0, nrow = nrow(self$delta0()), 
                         ncol = self$nstates(), byrow = TRUE)
      }
      
      # Check format of matrix input
      if(nrow(delta0) != length(self$unique_ID()) | ncol(delta0) != self$nstates()) {
        stop(paste0("'delta0' should have ", length(self$unique_ID()), " rows and ", 
                    self$nstates(), " columns"))
      }
      
      # Normalise rows if necessary
      if(any(abs(rowSums(delta0) - 1) > 1e-10)) {
        wng <- paste0("At least one row of delta0 doesn't sum to 1 (sum = ",
                      round(max(abs(rowSums(delta0))), 2), 
                      "). Normalising probabilities.")
        warning(wng)
        delta0 <- delta0/rowSums(delta0)
      }
      
      # Copy delta0 into object attribute
      private$delta0_ <- delta0
      rownames(private$delta0_) <- self$unique_ID()
      colnames(private$delta0_) <- paste0("state", 1:ncol(delta0))
    }, 
    
    #' @description Update smoothness parameters
    #' 
    #' @param lambda New smoothness parameter vector
    update_lambda = function(lambda) {
      private$lambda_ <- matrix(lambda)
      rownames(private$lambda_) <- self$terms()$names_re
    },
    
    #' @description Update information about fixed parameters
    #' 
    #' @param fixpar New list of fixed parameters, in the same format
    #' expected by MarkovChain$new()
    update_fixpar = function(fixpar) {
      private$fixpar_ <- fixpar
      private$fixpar_user_ <- fixpar
      private$setup_fixpar()
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
    #' Apply the multinomial logit link function to get the corresponding 
    #' parameters on the working scale (i.e., linear predictor scale).
    #' 
    #' @param tpm Transition probability matrix
    #' 
    #' @return Vector of parameters on linear predictor scale
    tpm2par = function(tpm) {
      ltpm <- log(tpm / tpm[which(self$ref_mat() == 1)])
      ltpm <- t(ltpm) # transpose to fill by rows (like in C++)
      par <- ltpm[!t(self$ref_mat())]
      return(par)
    },
    
    #' @description Transform working parameters to transition probabilities
    #' 
    #' Apply the inverse multinomial logit link function to transform the 
    #' parameters on the working scale (i.e., linear predictor scale) into 
    #' the transition probabilities.
    #' 
    #' @param par Vector of parameters on working scale
    #' 
    #' @return Transition probability matrix
    par2tpm = function(par) {
      tpm <- matrix(1, nrow = self$nstates(), ncol = self$nstates())
      tpm[!t(self$ref_mat())] <- exp(par)
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
      X_fe_old <- self$X_fe()
      X_re_old <- self$X_re()
      self$update_X_fe(mats_hid$X_fe)
      self$update_X_re(mats_hid$X_re)
      tpms <- self$tpm(t = "all")
      
      # Initial distribution
      delta0 <- self$delta0()
      
      # Simulate state process
      S <- rep(NA, n)
      S[1] <- sample(1:n_states, size = 1, prob = delta0[1,])
      id <- 1 # time series ID
      for(i in 2:n) {
        if(!silent & round(i/n*100)%%10 == 0) {
          cat("\rSimulating states... ", round(i/n*100), "%", sep = "")        
        }
        
        if(ID[i] != ID[i-1]) {
          id <- id + 1
          S[i] <- sample(1:n_states, size = 1, prob = delta0[id,])
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
      hid_forms <- as.data.frame(self$formula())
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
    
    formula_ = NULL,
    stationary_ = NULL, 
    formulas_ = NULL,
    ref_ = NULL,
    ref_mat_ = NULL,
    ref_delta0_ = NULL,
    coeff_fe_ = NULL,
    coeff_re_ = NULL,
    delta0_ = NULL,
    lambda_ = NULL,
    nstates_ = NULL,
    unique_ID_ = NULL,
    terms_ = NULL,
    fixpar_user_ = NULL,
    fixpar_ = NULL,
    initial_state_ = NULL,
    
    # Setup fixed parameters
    setup_fixpar = function() {
      ldelta0 <- self$delta0(log = TRUE, as_matrix = FALSE)
      if(is.numeric(self$initial_state())) {
        # Don't estimate delta0 if initial state is known
        private$fixpar_$delta0 <- rep(NA, length = length(ldelta0))
        names(private$fixpar_$delta0) <- rownames(ldelta0)
      } else if(self$initial_state() %in% c("stationary", "uniform")) {
        # Don't estimate delta0 if stationary or uniform
        private$fixpar_$delta0 <- rep(NA, length = length(ldelta0))
        names(private$fixpar_$delta0) <- rownames(ldelta0)
      } else if(self$initial_state() == "shared") {
        # delta0 is shared for all IDs
        private$fixpar_$delta0 <- rep(
          1:(self$nstates() - 1), each = length(ldelta0)/(self$nstates() - 1))
        names(private$fixpar_$delta0) <- rownames(ldelta0)
      }
      
      # Check for transitions that have fixed probabilities 
      ls_form_char <- as.list(t(self$formula())[!t(self$ref_mat())])
      which_fixed <- which(sapply(ls_form_char, function(x) {x == "."}))
      if(length(which_fixed) > 0) {
        getnms <- rownames(self$coeff_fe())[which_fixed]
        oldnms <- names(private$fixpar_$hid)
        private$fixpar_$hid <- c(private$fixpar_$hid, rep(NA, length(getnms)))
        names(private$fixpar_$hid) <- c(oldnms, getnms)        
      }
    },
    
    # Check constructor arguments
    # (For argument description, see constructor)
    check_args = function(n_states, formula, data) {
      if(!is.null(n_states)) {
        if(!is.numeric(n_states) | n_states < 1) {
          stop("'n_states' should be a numeric >= 1")
        }        
      }
      
      if(!is.null(formula)) {
        if(is.matrix(formula)) {
          if(nrow(formula) != ncol(formula)) {
            stop("'formula' should be a square matrix")
          }
          
          if(!all(is.character(as.vector(formula)))) {
            stop("'formula' should be a matrix of character strings")
          }
          
          if (!all(diag(formula) == ".")) {
            stop("Diagonal of formula should be '.'")
          }
          
        } else if(!inherits(formula, "formula")) {
          stop("'formula' should be either a matrix or a formula")
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
