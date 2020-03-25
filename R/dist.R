
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
    initialize = function(name, pdf, link, invlink) {
      private$name_ <- name
      private$pdf_ <- pdf
      private$link_ <- link
      private$invlink_ <- invlink
      distnames <- c("pois", "norm", "gamma", "beta", "custom")
      private$code_ <- which(distnames == name) - 1 # Starts at 0 for C++
    },
    
    # Accessors
    name = function() {return(private$name_)},
    pdf = function() {return(private$pdf_)},
    link = function() {return(private$link_)},
    invlink = function() {return(private$invlink_)},
    code = function() {return(private$code_)},
    
    # Transform parameters from natural to working scale
    n2w = function(par) {
      # Apply link functions to natural parameters
      wpar_list <- Map(function(fn, arg) {fn(arg)}, self$link(), par)
      wpar <- unlist(wpar_list)
      return(wpar)
    },
    
    # Transform parameters from working to natural scale
    w2n = function(wpar) {
      par <- list()
      
      # Number of parameters for this distribution
      n_par <- length(self$link())
      # Number of states
      n_state <- length(wpar)/n_par
      
      # Loop over parameters and apply inverse link functions
      for(i in 1:n_par) {
        i0 <- (i-1)*n_state + 1
        i1 <- i*n_state
        sub_wpar <- wpar[i0:i1]
        par[[i]] <- self$invlink()[[i]](sub_wpar)
        names(par[[i]]) <- NULL
      }
      
      names(par) <- names(self$link())
      return(par)
    }
  ),
  
  private = list(
    name_ = NULL,
    pdf_ = NULL,
    link_ = NULL,
    invlink_ = NULL,
    code_ = NULL
  )
)
