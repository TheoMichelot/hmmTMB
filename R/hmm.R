
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
      
      # Number of states
      n_states <- self$hidden()$nstates()
      
      # Create model matrices of observation process
      # (Design matrices for fixed and random effects, and smoothing matrix)
      mod_mat_obs <- self$obs()$make_mat()
      X_fe_obs <- mod_mat_obs$X_fe
      X_re_obs <- mod_mat_obs$X_re
      S_obs <- mod_mat_obs$S
      ncol_re_obs <- mod_mat_obs$ncol_re
      
      # Create model matrices of hidden state process
      # (Design matrices for fixed and random effects, and smoothing matrix)
      mod_mat_hid <- self$hidden()$make_mat(data = self$obs()$data()$data())
      X_fe_hid <- mod_mat_hid$X_fe
      X_re_hid <- mod_mat_hid$X_re
      S_hid <- mod_mat_hid$S
      ncol_re_hid <- mod_mat_hid$ncol_re
      
      # Setup TMB parameters
      tmb_par <- list(wpar_fe_obs = self$obs()$tpar(),
                      wpar_re_obs = 0,
                      log_lambda_obs = 0,
                      wpar_fe_hid = rep(0, ncol(X_fe_hid)),
                      wpar_re_hid = 0,
                      log_lambda_hid = 0,
                      log_delta = rep(0, n_states - 1))
      
      # Initialise vectors of fixed parameters and random effects
      map <- NULL
      random <- NULL
      
      # Setup random effects in observation model
      if(is.null(S_obs)) {
        # If there are no random effects, 
        # wpar_re and log_lambda are not estimated
        map <- c(map, list(wpar_re_obs = factor(NA),
                           log_lambda_obs = factor(NA)))
        S_obs <- as(matrix(0, 1, 1), "sparseMatrix")
        ncol_re_obs <- 0
        X_re_obs <- as(rep(0, nrow(X_fe_obs)), "sparseMatrix")
      } else {
        # If there are random effects, 
        # set initial values for wpar_re and log_lambda
        random <- c(random, "wpar_re_obs")
        tmb_par$wpar_re_obs <- rep(0, ncol(S_obs))
        tmb_par$log_lambda_obs <- rep(0, length(ncol_re_obs))
      }
      
      # Setup random effects in hidden state model
      if(is.null(S_hid)) {
        # If there are no random effects, 
        # wpar_re and log_lambda are not estimated
        map <- c(map, list(wpar_re_hid = factor(NA),
                           log_lambda_hid = factor(NA)))
        S_hid <- as(matrix(0, 1, 1), "sparseMatrix")
        ncol_re_hid <- 0
        X_re_hid <- as(rep(0, nrow(X_fe_hid)), "sparseMatrix")
      } else {
        # If there are random effects, 
        # set initial values for wpar_re and log_lambda
        random <- c(random, "wpar_re_hid")
        tmb_par$wpar_re_hid <- rep(0, ncol(S_hid))
        tmb_par$log_lambda_hid <- rep(0, length(ncol_re_hid))
      }
      
      # Data for TMB
      tmb_dat <- list(ID = self$obs()$data()$ID(),
                      data = as.matrix(self$obs()$obs_var()),
                      n_states = n_states,
                      distcode = distcode,
                      X_fe_obs = X_fe_obs,
                      X_re_obs = X_re_obs,
                      S_obs = S_obs,
                      ncol_re_obs = ncol_re_obs,
                      X_fe_hid = X_fe_hid,
                      X_re_hid = X_re_hid,
                      S_hid = S_hid,
                      ncol_re_hid = ncol_re_hid)
      
      # Create TMB model
      obj <- MakeADFun(tmb_dat, tmb_par, dll = "HmmTmb", 
                       random = random,
                       map = map)
      
      # Fit model
      private$fit_ <- do.call(optim, obj)
      
      # Update model parameters
      est_par <- self$res()$par

      # Observation parameters
      ind_wpar <- which(names(est_par) == "wpar_fe_obs")
      wpar <- est_par[ind_wpar]
      self$obs()$update_wpar(wpar = wpar, n_state = n_states)
      
      # Transition probabilities
      ind_ltpm <- which(names(est_par) == "wpar_fe_hid")
      if(length(ind_ltpm) == n_states * (n_states - 1)) {
        # Only update if no covariates
        ltpm <- est_par[ind_ltpm]
        self$hidden()$update_par(ltpm)        
      }
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





