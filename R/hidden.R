
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
        if(length(structure) == 1) {
          # Same formula for all transitions
          structure <- matrix(structure, nrow = n_states, ncol = n_states)
          diag(structure) <- "."
        } else {
          n_states <- nrow(structure)
        }
        
        private$nstates_ <- n_states
        
        # Set initial parameters 
        private$par_ <- par
      }
      
      # Set remaining private attributes
      private$structure_ <- structure
      private$par_re_ <- integer(0)  # so that X_re %*% par_re is valid
    },
    
    ###############
    ## Accessors ##
    ###############
    #' @description Formulas for MarkovChain model
    structure = function() {return(private$structure_)},
    
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
    #' 
    #' @return A list with elements:
    #' \itemize{
    #'   \item{X_fe}{Design matrix for fixed effects}
    #'   \item{X_re}{Design matrix for random effects}
    #'   \item{S}{Smoothness matrix for random effects}
    #'   \item{ncol_re}{Number of columns of X_re and S for each random effect}
    #' }
    make_mat = function(data) {
      struct <- self$structure()[!diag(self$nstates())]
      formulas <- lapply(as.list(struct), function(string) {
        if(string == ".")
          return(NULL)
        else
          return(as.formula(string))
      })
      
      make_mat_hid(formulas = formulas, data = data)
    },
    
    #' @description Transition probability matrices
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
      ltpm_mat <- matrix(ltpm, ncol = n_states)
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
    }
  ),
  
  private = list(
    structure_ = NULL,
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





