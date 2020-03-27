
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
#'  \item fit: fit the model
#'  \item res: fitted model object, after optimisation
#'  \item est: parameter estimates
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
      
      mod_mat <- self$obs()$make_mat(n_states = self$hidden()$nstates())
      X_fe <- mod_mat$X_fe
      X_re <- mod_mat$X_re
      S <- mod_mat$S
      ncol_re <- mod_mat$ncol_re

      tmb_par <- list(ltpm = self$hidden()$par(),
                      wpar_fe = self$obs()$tpar(),
                      wpar_re = 0,
                      log_lambda = 0)

      map <- NULL
      random <- NULL
      if(is.null(S)) {
        map <- c(map, list(wpar_re = factor(NA),
                           log_lambda = factor(NA)))
        S <- as(matrix(0, 1, 1), "sparseMatrix")
        ncol_re <- 0
      } else {
        random <- c(random, "wpar_re")
        tmb_par$wpar_re <- rep(0, ncol(S))
        tmb_par$log_lambda <- rep(0, length(ncol_re))
      }
      
      tmb_dat <- list(ID = self$obs()$data()$ID(),
                      data = as.matrix(self$obs()$obs_var()),
                      X_fe = X_fe,
                      X_re = X_re,
                      S = S,
                      ncol_re = ncol_re,
                      n_states = self$hidden()$nstates(),
                      distcode = distcode)

      obj <- MakeADFun(tmb_dat, tmb_par, dll = "HmmTmb", 
                       random = random,
                       map = map)
      
      private$fit_ <- do.call(optim, obj)
      
      # Update model parameters
      est_par <- self$res()$par
      n_state <- self$hidden()$nstates()
      
      # Observation parameters
      ind_wpar <- which(names(est_par) == "wpar_fe")
      wpar <- est_par[ind_wpar]
      self$obs()$update_wpar(wpar = wpar, n_state = n_state)
      
      # Transition probabilities
      ind_ltpm <- which(names(est_par) == "ltpm")
      ltpm <- est_par[ind_ltpm]
      self$hidden()$update_par(ltpm)
    },
    
    # Parameter estimates
    est = function() {
      par <- self$obs()$par()
      tpm <- self$hidden()$tpm()
      return(list(par = par, tpm = tpm))
    }
  ),
  
  private = list(
    obs_ = NULL,
    hidden_ = NULL,
    fit_ = NULL
  )
)





