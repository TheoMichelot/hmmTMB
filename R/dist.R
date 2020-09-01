
#' R6 class for probability distribution
#' 
#' Contains the probability density/mass function, and the link and inverse link 
#' functions for a probability distribution.
Dist <- R6Class(
  classname = "Dist",
  
  public = list(
    #################
    ## Constructor ##
    #################
    #' @description Create a Dist object
    #' 
    #' @param name Name of distribution
    #' @param pdf Probability density/mass function of the distribution
    #' (e.g. \code{dnorm} for normal distribution).
    #' @param rng Random generator function of the distribution (e.g.
    #' \code{rnorm} for normal distribution).
    #' @param link Named list of link functions for distribution parameters
    #' @param invlink Named list of inverse link functions for distribution
    #' parameters
    #' @param npar Number of parameters of the distribution
    #' 
    #' @return A new Dist object
    initialize = function(name, pdf, rng, link, invlink, npar) {
      # Check arguments
      private$check_args(name = name, pdf = pdf, rng = rng, link = link, 
                         invlink = invlink, npar = npar)
      
      # Define object attributes
      private$name_ <- name
      private$pdf_ <- pdf
      private$rng_ <- rng
      private$link_ <- link
      private$invlink_ <- invlink
      private$npar_ <- npar
      
      # List of distributions included in the package
      distnames <- c("pois", "norm", "gamma", "beta", "vm", "custom")
      if(!name %in% distnames) {
        stop(paste0("'name' must be one of '", 
                    paste(distnames, sep = "", collapse = "', '"), 
                    "'"))
      } else {
        private$code_ <- which(distnames == name) - 1 # Starts at 0 for C++        
      }
    },
    
    
    ###############
    ## Accessors ##
    ###############
    #' @description Return name of Dist object
    name = function() {return(private$name_)},
    
    #' @description Return pdf of Dist object
    pdf = function() {return(private$pdf_)},
    
    #' @description Return random generator function of Dist object
    rng = function() {return(private$rng_)},
    
    #' @description Return link function of Dist object
    link = function() {return(private$link_)},
    
    #' @description Return inverse link function 
    #' of Dist object
    invlink = function() {return(private$invlink_)},
    
    #' @description Return number of parameters of Dist object
    npar = function() {return(private$npar_)},
    
    #' @description Return code of Dist object
    code = function() {return(private$code_)},
    
    ###################
    ## Other methods ##
    ###################
    #' @description Evaluate probability density/mass function
    #' 
    #' (Used in Observation$obs_probs)
    #' 
    #' @param x Value at which the function should be evaluated
    #' @param par Vector of parameters. The entries should be named if
    #' they are not in the same order as expected by the R function. (E.g.
    #' shape/scale rather than shape/rate for gamma distribution.)
    #' @param log Logical. If TRUE, the log-density is returned. 
    #' Default: FALSE.
    #' 
    #' @return Probability density/mass function evaluated at x for
    #' parameters par
    pdf_apply = function(x, par, log = FALSE) {
      args <- list(x = x)
      args <- c(args, par, log = log)
      do.call(self$pdf(), args)
    },
    
    #' @description Random number generator
    #' 
    #' @param n Number of realisations to generate
    #' @param par Vector of parameters. The entries should be named if
    #' they are not in the same order as expected by the R function. (E.g.
    #' shape/scale rather than shape/rate for gamma distribution.)
    #' 
    #' @return Vector of \code{n} realisations of this distribution
    rng_apply = function(n, par) {
      args <- list(n = n)
      args <- c(args, par)
      do.call(self$rng(), args)
    },
    
    #' @description Natural to working parameter transformation
    #' 
    #' @param par List of parameters
    #' 
    #' @return Vector of parameters on the working scale
    n2w = function(par) {
      # Apply link functions to natural parameters
      wpar_list <- Map(function(fn, arg) {fn(arg)}, self$link(), par)
      wpar <- unlist(wpar_list)
      return(wpar)
    },
    
    #' @description Working to natural parameter transformation
    #' 
    #' @param wpar Vector of working parameters
    #' @param as_matrix Logical. If TRUE, the natural parameters are
    #' returned as a matrix with one row for each state and one column
    #' for each parameter. If FALSE, the natural parameters are returned
    #' as a list (default).
    #' 
    #' @return List or matrix of parameters on natural scale
    w2n = function(wpar, as_matrix = FALSE) {
      invlink <- self$invlink()
      # Number of parameters for this distribution
      n_par <- length(invlink)
      # Number of states
      n_states <- length(wpar)/n_par
      
      # Apply inverse link functions
      par_list <- lapply(seq_along(invlink), function(i) {
        ind <- ((i-1) * n_states + 1) : (i * n_states)
        invlink[[i]](wpar[ind])
      })
      
      # Format into matrix or list
      if(as_matrix) {
        par <- matrix(unlist(par_list), nrow = n_states)
        colnames(par) <- names(invlink)        
      } else {
        par <- par_list
        names(par) <- names(invlink)        
      }
      
      return(par)
    }
  ),
  
  private = list(
    name_ = NULL,
    pdf_ = NULL,
    rng_ = NULL,
    link_ = NULL,
    invlink_ = NULL,
    npar_ = NULL,
    code_ = NULL,
    
    # Check arguments passed to constructor
    # (For argument description, see constructor)
    check_args = function(name, pdf, rng, link, invlink, npar) {
      if(!is.character(name)) {
        stop("'name' must be a character string")
      }
      
      if(!is.function(pdf)) {
        stop(paste("'pdf' must be an R function for the probability density",
                   "function of the distribution (and not the name of the function",
                   "as a string)."))
      }
      
      if(!is.function(rng)) {
        stop(paste("'rng' must be an R function for the random number generator",
                   "of the distribution (and not the name of the function",
                   "as a string)."))
      }
      
      if(length(link) != length(invlink) | !all(names(link) == names(invlink))) {
        stop(paste("'link' and 'invlink' should be lists of the same length, and with", 
                   "the same names (corresponding to the names of the distribution",
                   "parameters"))
      }
      
      if(!all(sapply(link, is.function)) | !all(sapply(invlink, is.function))) {
        stop(paste("Elements of 'link' and 'invlink' should be R functions"))
      }
      
      if(!is.numeric(npar) | npar < 1) {
        stop("'npar' should be a numeric >= 1")
      }
    }
  )
)
