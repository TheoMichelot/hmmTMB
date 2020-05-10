
#' Distribution class
#'
#' @description Distribution for observed variables.
#' Object can be created using $new with arguments:
#' \itemize{
#'   \item name Name of distribution
#'   \item pdf Probability density/mass function of the distribution
#'   \item link Named list of link functions for distribution parameters
#'   \item invlink Named list of inverse link functions for distribution
#'   parameters
#' }
#'
#' Methods include:
#' \itemize{
#'  \item n2w Transforms parameters from natural to working scale
#'  \item w2n Transforms parameters from working to natural scale
#' }

Dist <- R6Class(
  classname = "Dist",
  
  public = list(
    initialize = function(name, pdf, link, invlink, npar) {
      private$name_ <- name
      private$pdf_ <- pdf
      private$link_ <- link
      private$invlink_ <- invlink
      private$npar_ <- npar
      distnames <- c("pois", "norm", "gamma", "beta", "custom")
      private$code_ <- which(distnames == name) - 1 # Starts at 0 for C++
    },
    
    # Accessors
    name = function() {return(private$name_)},
    pdf = function() {return(private$pdf_)},
    link = function() {return(private$link_)},
    invlink = function() {return(private$invlink_)},
    npar = function() {return(private$npar_)},
    code = function() {return(private$code_)},
    
    # Evaluate the pdf of x for a generic parameter vector par
    pdf_apply = function(x, par, log = FALSE) {
      args <- list(x = x)
      args <- c(args, par, log = log)
      do.call(self$pdf(), args)
    },
    
    # Transform parameters from natural to working scale
    n2w = function(par) {
      # Apply link functions to natural parameters
      wpar_list <- Map(function(fn, arg) {fn(arg)}, self$link(), par)
      wpar <- unlist(wpar_list)
      return(wpar)
    },
    
    # Transform parameters from working to natural scale
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
    link_ = NULL,
    invlink_ = NULL,
    npar_ = NULL,
    code_ = NULL
  )
)
