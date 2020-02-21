
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

Hmm <- R6Class("Hmm",

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

      tmb_dat <- list(data = self$obs()$data()$data()[,1],
                            n_states = self$hidden()$nstates(),
                      distname = self$obs()$dists()[[1]]$name())

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





