
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
      if (is.null(private$fit_)) {
        stop("Fit model first")
      }
      
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
    },
    
    # Parameter estimates
    est = function() {
      est_par <- self$res()$par
      n_state <- self$hidden()$nstates()
      
      # Observation parameters
      ind_wpar <- which(names(est_par) == "wpar")
      wpar <- est_par[ind_wpar]
      par <- self$obs()$w2n(wpar = wpar, n_state = n_state)
      
      # Transition probabilities
      ind_ltpm <- which(names(est_par) == "ltpm")
      ltpm <- est_par[ind_ltpm]
      self$hidden()$update_par(ltpm)
      
      return(list(obspar = par, tpm = self$hidden()$tpm()))
    }
  ),
  
  private = list(
    obs_ = NULL,
    hidden_ = NULL,
    fit_ = NULL
  )
)





