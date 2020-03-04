
#' Hidden Markov model class
#'
#' @description Encapsulates the hidden Markov model.
#' Object can be created using $new with arguments:
#' \itemize{
#'   \item obs: a Observations object
#'   \item hidden: a MarkovChain object
#' }
#'
#' Methods include:
#' \itemize{
#'  \item  NA
#' }

Hmm <- R6Class(
  classname = "Hmm",
  
  public = list(
    initialize = function(obs, hidden) {
      private$obs_ <- obs
      private$hidden_ <- hidden
    },
    
    # Accessors
    obs = function() {return(private$obs_)},
    hidden = function() {return(private$hidden_)},
    res = function() {
      if (is.null(private$fit_)) stop("Fit model first")
      return(private$fit_)
    },
    
    # Fitting
    fit = function() {
      # Vector of codes of observation distributions
      distcode <- as.vector(sapply(self$obs()$dists(), function(d) d$code()))
      
      tmb_dat <- list(data = as.matrix(self$obs()$data()$data()),
                      n_states = self$hidden()$nstates(),
                      distcode = distcode)
      
      tmb_par <- list(ltpm = self$hidden()$par(),
                      wpar = self$obs()$tpar())
      
      obj <- MakeADFun(tmb_dat, tmb_par, dll = "HmmTmb")
      
      private$fit_ <- do.call(optim, obj)
      
    }
    
  ),
  
  private = list(
    obs_ = NULL,
    hidden_ = NULL,
    fit_ = NULL
    
  )
)





