
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
    tmb_obj = function() {
      if(is.null(private$tmb_obj_)) {
        stop("Setup model first")
      }
      
      return(private$tmb_obj_)
    },
    tmb_rep = function() {
      if(is.null(private$tmb_rep_)) {
        stop("Fit model first")
      }
      
      return(private$tmb_rep_)
    },
    states = function() {
      if(is.null(private$states_)) {
        stop("Run viterbi first")
      }
      
      return(private$states_)
    },
    
    # Objective function
    nllk = function(par) {
      self$tmb_obj()$fn(par)
    },
    
    # TMB setup
    setup = function() {
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
      tmb_par <- list(wpar_fe_obs = self$obs()$wpar(),
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
      obj <- MakeADFun(tmb_dat, tmb_par, dll = "hmmTMB", 
                       random = random,
                       map = map)
      
      # Negative log-likelihood function
      private$tmb_obj_ <- obj
    },

    # Fitting
    fit = function() {
      # Setup if necessary
      if(is.null(private$tmb_obj_)) {
        self$setup()
      }
      
      # Number of states
      n_states <- self$hidden()$nstates()
      
      # Fit model
      private$fit_ <- do.call(optim, private$tmb_obj_)
      
      # Get estimates and precision matrix for all parameters
      private$tmb_rep_ <- sdreport(private$tmb_obj_)
      par_list <- as.list(private$tmb_rep_, "Estimate")
      
      # Observation parameters
      self$obs()$update_wpar(wpar = par_list$wpar_fe_obs, n_states = n_states)
      mats_obs <- self$obs()$make_mat()
      if(!is.null(mats_obs$ncol_re)) { # Only update if there are random effects
        self$obs()$update_wpar_re(wpar = par_list$wpar_re_obs)
      }
      
      # Transition probabilities
      self$hidden()$update_par(newpar = par_list$wpar_fe_hid)
      mats_hid <- self$hidden()$make_mat(data = self$obs()$data()$data())
      if(!is.null(mats_hid$ncol_re)) { # Only update if there are random effects
        self$hidden()$update_par_re(newpar = par_list$wpar_re_hid)        
      }
    },
    
    # Wald confidence intervals for the parameters on working scale
    CI_wpar = function(level = 0.95) {
      if(is.null(private$tmb_rep_)) {
        stop("Fit model first")
      }
      
      par_list <- as.list(private$tmb_rep_, "Estimate")
      se_list <- as.list(private$tmb_rep_, "Std. Error")
      
      lower <- lapply(seq_along(par_list), function(i) {
        par_list[[i]] - qnorm(1 - (1 - level)/2) * se_list[[i]]
      })
      upper <- lapply(seq_along(par_list), function(i) {
        par_list[[i]] + qnorm(1 - (1 - level)/2) * se_list[[i]]
      })
      
      return(cbind(estimate = unlist(par_list),
                   lower = unlist(lower),
                   upper = unlist(upper)))
    },
    
    # Model parameters
    par = function() {
      obspar <- self$obs()$par()
      tpm <- self$hidden()$tpm()
      return(list(obspar = obspar, tpm = tpm))
    },
    
    # Viterbi algorithm
    viterbi = function() {
      data <- self$obs()$data()$data()
      ID <- self$obs()$data()$ID()
      
      # Number of observations
      n <- nrow(data)
      # Number of states
      n_states <- self$hidden()$nstates()
      # Number of variables
      n_var <- length(self$obs()$dists())
      
      # Observation probabilities
      mod_mat_obs <-  self$obs()$make_mat()
      obs_probs <- self$obs()$obs_probs(X_fe = mod_mat_obs$X_fe, 
                                        X_re = mod_mat_obs$X_re)
      
      # Transition probability matrices      
      mod_mat_hid <- self$hidden()$make_mat(data = self$obs()$data()$data())
      tpm_all <- self$hidden()$tpm_all(X_fe = mod_mat_hid$X_fe, 
                                       X_re = mod_mat_hid$X_re,
                                       n = n)
      
      # Number of unique IDs
      n_id <- length(unique(ID))
      # First index for each ID
      i0 <- c(1, which(ID[-1] != ID[-n]) + 1, n + 1)
      
      # Initialise state sequence
      all_states <- NULL
      
      # For now, uniform initial distribution
      delta <- rep(1/n_states, n_states)
      
      # Loop over IDs
      for(id in 1:n_id) {
        # Subset to this ID
        ind_this_id <- which(ID == unique(ID)[id])
        sub_obs_probs <- obs_probs[ind_this_id,]
        sub_tpm_all <- tpm_all[,,ind_this_id]
        
        # Number of observations for this ID
        n_this_id <- length(ind_this_id)
        
        # Forward iterations
        xi <- matrix(NA, n_this_id, n_states)
        v <- delta * sub_obs_probs[1,]
        xi[1,] <- v/sum(v)
        for(i in 2:n_this_id) {
          v <- apply(xi[i-1,] * sub_tpm_all[,,i], 2, max) * sub_obs_probs[i,]
          xi[i,] <- v/sum(v)
        }
        
        # Backward iterations
        states <- rep(NA, n_this_id)
        states[n_this_id] <- which.max(xi[n_this_id,])
        for(i in (n_this_id - 1):1) {
          states[i] <- which.max(sub_tpm_all[, states[i+1], i+1] * xi[i,])
        }
        
        # Append estimated states for this ID
        all_states <- c(all_states, states)
      }
      
      # Save state sequence
      private$states_ <- all_states
      
      return(all_states)
    }
  ),
  
  private = list(
    obs_ = NULL,
    hidden_ = NULL,
    fit_ = NULL,
    tmb_obj_ = NULL,
    tmb_rep_ = NULL,
    states_ = NULL
  )
)
