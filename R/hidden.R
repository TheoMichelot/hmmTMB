
#' Hidden Markov chain class
#'
#' @description Encapsulates the Markov chain for the hidden component of the HMM.
#' Object can be created using $new with arguments:
#' \itemize{
#'   \item structure: matrix with an entry of "." on diagonal, a "0" for 
#'   transitions that are not allowed (not implemented yet), and a formula "~1" 
#'   for covariates affecting transitions that are to be estimated.
#'   \item tpm: an initial transition probability matrix.
#' }
#'
#' Methods include:
#' \itemize{
#'  \item  structure: returns the specified structure of the Markov chain
#'  \item tpm: return current transition probability matrix
#'  \item par: return current parameter estimates for transitions
#'  \item nstates: return number of states in Markov chain
#'  \item update_par (newpar): set parameters to newpar
#'  \item update_tpm (newtpm): set transition probability matrix to newtpm
#' }

MarkovChain <- R6Class(
  classname = "MarkovChain",
  
  public = list(
    initialize = function(n_states = 2, structure = NULL, tpm = NULL) {
      if(is.null(structure)) {
        # Default structure: no covariate effects
        structure <- matrix("~1", nrow = n_states, ncol = n_states)
        diag(structure) <- "."
      } else if(length(structure) == 1) {
        # Same formula for all transitions
        structure <- matrix(structure, nrow = n_states, ncol = n_states)
        diag(structure) <- "."
      } else {
        nstates_ <- nrow(structure)
      }
      
      # Set default initial tpm (0.9 on diagonal)
      if(is.null(tpm)) {
        tpm <- matrix(0.1/(n_states-1), nrow = n_states, ncol = n_states)
        diag(tpm) <- 0.9
      }
      
      # Set private attributes
      private$nstates_ <- n_states
      private$structure_ <- structure
      private$tpm_ <- tpm
      private$par_ <- private$tpm2par(tpm)
      private$par_re_ <- integer(0)  # so that X_re %*% par_re is valid
    },
    
    # Accessors
    structure = function() {return(private$structure_)},
    tpm = function() {return(private$tpm_)},
    par = function() {return(private$par_)},
    par_re = function() {return(private$par_re_)},
    nstates = function() {return(private$nstates_)},
    
    # Mutators
    update_tpm = function(newtpm) {
      private$tpm_ <- newtpm
      private$par_ <- private$tpm2par(newtpm)
    },
    update_par = function(newpar) {
      private$par_ <- newpar
      if(all(self$structure() %in% c(".", "~1"))) {
        # Only update tpm if no covariates
        private$tpm_ <- private$par2tpm(newpar)        
      }
    },
    update_par_re = function(newpar) {
      private$par_re_ <- newpar
    },
    
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
    
    tpm_all = function(X_fe, X_re, n) {
      ltpm <- X_fe %*% self$par() + X_re %*% self$par_re()
      ltpm_mat <- matrix(ltpm, nrow = n)
      tpm <- apply(ltpm_mat, 1, private$par2tpm)
      tpm <- array(tpm, dim = c(self$nstates(), self$nstates(), n))
      return(tpm)
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





