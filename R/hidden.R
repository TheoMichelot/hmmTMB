
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
    initialize = function(structure, tpm) {
      private$structure_ <- structure
      private$nstates_ <- nrow(structure)
      private$tpm_ <- tpm
      private$par_ <- private$tpm2par(tpm)
    },
    
    # Accessors
    structure = function() {return(private$structure_)},
    tpm = function() {return(private$tpm_)},
    par = function() {return(private$par_)},
    nstates = function() {return(private$nstates_)},
    
    # Mutators
    update_par = function(newpar) {
      private$par_ <- newpar
      private$tpm_ <- private$par2tpm(newpar)
    },
    update_tpm = function(newtpm) {
      private$tpm_ <- newtpm
      private$par_ <- private$tpm2par(newtpm)
    },
    
    make_mat = function(data) {
      struct <- self$structure()[!diag(self$nstates())]
      formulas <- lapply(as.list(struct), function(string) {
        if(string == ".")
          return(NULL)
        else
          return(as.formula(string))
      })
      
      make_mat_2(formulas = formulas, data = data)
    }
  ),
  
  
  private = list(
    structure_ = NULL,
    par_ = NULL,
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





