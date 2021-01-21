
#' R6 class for hidden Markov model
#'
#' Encapsulates the observation and hidden state models for a hidden
#' Markov model.
HMM <- R6Class(
  classname = "HMM",
  
  public = list(

    # Constructor -------------------------------------------------------------
    #' @description Create new HMM object
    #' 
    #' @param file path to specification file for HMM 
    #' @param obs Observation object
    #' @param hidden MarkovChain object
    #' @param init HMM object to initialize parameters with 
    #' @param fixpar a named list of parameters in coeff_fe that you want fixed 
    #' (set to NA) or pool into common values (using factor levels)
    #' 
    #' @return A new HMM object
    initialize = function(file = NULL, 
                          obs = NULL, 
                          hidden = NULL, 
                          init = NULL, 
                          fixpar = NULL) {
      # Decide how model has been specified 
      if (is.null(file) & is.null(obs)) {
        stop("either 'file' must a file name for a file specifying the model or 
            both obs/hidden model objects should be supplied.")
      }
      if (!is.null(file) & is.null(obs)) {
        spec <- private$read_file(file)
        # create obs
        obs <- Observation$new(data = spec$data, 
                               dists = spec$dists,
                               n_states = spec$nstates, 
                               par = spec$par, 
                               formulas = spec$forms)
        hidden <- MarkovChain$new(n_states = spec$nstates, 
                                  structure = spec$tpm, 
                                  data = spec$data)
        if (!is.null(spec$fixed)) fixpar <- spec$fixed 
        if (!is.null(spec$tpm0))  hidden$update_tpm(spec$tpm0)
      }
      
      # Check arguments
      private$check_args(obs = obs, hidden = hidden, init = init)
      
      # Get names of all covariates
      var_names <- unique(c(rapply(hidden$formulas(), all.vars), 
                          rapply(obs$formulas(), all.vars)))
      # Remove pi from list of covariates if it is in the formulas
      var_names <- var_names[which(var_names!="pi")]
      if(length(var_names) > 0) {
        data <- obs$data()
        # Remove NAs in covariates (replace by last non-NA value)
        data[,var_names] <- lapply(data[,var_names, drop=FALSE], 
                                   function(col) na_fill(col))
        # Update data frame in obs
        obs$update_data(data)
      }
      
      # store sub-model components 
      private$obs_ <- obs
      private$hidden_ <- hidden
      
      # store fixed parameter 
      private$fixpar_ <- fixpar 
      
      # initialize model parameters if init provided 
      if (!is.null(init)) {
        private$obs_$update_coeff_fe(private$initialize_submodel(private$obs_$coeff_fe(), 
                                                         init$obs()$coeff_fe()))
        private$obs_$update_coeff_re(private$initialize_submodel(private$obs_$coeff_re(), 
                                                         init$obs()$coeff_re()))
        private$hidden_$update_coeff_fe(private$initialize_submodel(private$hidden_$coeff_fe(), 
                                                         init$hidden()$coeff_fe()))
        private$hidden_$update_coeff_re(private$initialize_submodel(private$hidden_$coeff_re(), 
                                                            init$hidden()$coeff_re()))
      }
      
      # initialize priors 
      self$set_priors()
    },
    

    # Accessors ---------------------------------------------------------------

    #' @description Observation object for this model
    obs = function() {return(private$obs_)},
    
    #' @description MarkovChain object for this model
    hidden = function() {return(private$hidden_)},
    
    #' @description Output of optimiser after model fitting
    out = function() {
      if (is.null(private$out_)) {
        stop("fit model first using $fit()")
      }
      return(private$out_)
    },
    
    #' @description Model object created by TMB. This is the output of 
    #' the TMB function \code{MakeADFun}, and it is a list including elements
    #' \itemize{
    #'   \item{\code{fn}}{Objective function}
    #'   \item{\code{gr}}{Gradient function of fn}
    #'   \item{\code{par}}{Vector of initial parameters on working scale}
    #' }
    tmb_obj = function() {
      if(is.null(private$tmb_obj_)) {
        stop("setup or fit model first")
      }
      return(private$tmb_obj_)
    },
    
    #' @description Output of the TMB function \code{sdreport}, which includes 
    #' estimates and standard errors for all model parameters.
    tmb_rep = function() {
      if(is.null(private$tmb_rep_)) {
        stop("fit model first")
      }
      return(private$tmb_rep_)
    },
    
    #' @description Vector of estimated states, after \code{viterbi} has
    #' been run
    states = function() {
      if(is.null(private$states_)) {
        stop("run viterbi first")
      }
      return(private$states_)
    },
    
    #' @description Coefficients for fixed effect parameters
    coeff_fe = function() {
      return(list(obs = self$obs()$coeff_fe(),
                  hidden = self$hidden()$coeff_fe()))
    },
    
    #' @description Coefficients for random effect parameters
    coeff_re = function() {
      return(list(obs = self$obs()$coeff_re(),
                  hidden = self$hidden()$coeff_re()))
    },
    
    #' @description Smoothness parameters
    lambda = function() {
      return(list(obs = self$obs()$lambda(),
                  hidden = self$hidden()$lambda()))
    },
    
    #' @description Update parameters stored inside model object
    #' 
    #' @param par_list a list for coeff_f(r)e_obs, coeff_f(r)e_hid, log_delta, 
    #'                   log_lambda_hid log_lambda_obs
    #' @param iter Optional argument to update model parameters based on MCMC
    #' iterations (if using rstan). Either the index of the iteration to use,
    #' or "mean" if the posterior mean should be used.
    update_par = function(par_list = NULL, iter = NULL) {
      if (is.null(par_list) & is.null(iter)) {
        stop("No new parameter values to update to")
      }
      if (!is.null(iter) & !is.null(par_list)) {
        stop("Either specify iter or par_list in update_par, not both")
      }
      if (!is.null(iter)) {
        # update to MCMC iteration 
        if (is.null(private$iters_)) {
          stop("Must run mcmc() before using iterations")
        }
        if (is.numeric(iter)) {
          if (iter > dim(private$iters_)[1]) {
            stop("iter exceeds number of mcmc iterations available")
          }
          samp <- private$iters_[iter,]
        } else if (iter == "mean") {
          samp <- colMeans(private$iters_)
        } else {
          stop("invalid iter to update_par()")
        }
        nms <- names(samp)
        nms <- gsub("\\[[^][]*\\]", "", nms)
        names(samp) <- NULL
        par_list <- split(samp,nms)
      }
      # Update observation parameters
      self$obs()$update_coeff_fe(coeff_fe = par_list$coeff_fe_obs)
      if(!is.null(self$obs()$terms()$ncol_re)) { 
        # Only update if there are random effects
        self$obs()$update_coeff_re(coeff = par_list$coeff_re_obs)
        self$obs()$update_lambda(exp(par_list$log_lambda_obs))
      }
      
      # Update transition probabilities
      self$hidden()$update_coeff_fe(coeff_fe = par_list$coeff_fe_hid)
      if(!is.null(self$hidden()$terms()$ncol_re)) { 
        # Only update if there are random effects
        self$hidden()$update_coeff_re(coeff_re = par_list$coeff_re_hid)
        self$hidden()$update_lambda(exp(par_list$log_lambda_hid))
      }
      
      # Update delta parameters 
      if (self$hidden()$stationary()) {
        tpms <- self$hidden()$tpm(t = "all")
        nstates <- self$hidden()$nstates()
        delta <- solve(t(diag(nstates) - tpms[,,1] + 1), rep(1, nstates))
        self$hidden()$update_delta(delta)
      } else {
        ldelta <- par_list$log_delta 
        delta <- c(exp(ldelta), 1)
        delta <- delta / sum(delta)
        self$hidden()$update_delta(delta)
      }
    }, 
    
    #' @description Variance components of smooth terms
    #' 
    #' This function transforms the smoothness parameter of
    #' each smooth term into a standard deviation, given by 
    #' SD = 1/sqrt(lambda). It is particularly helpful to get the
    #' standard deviations of independent normal random effects.
    vcomp = function() {
      if(is.null(private$tmb_rep_)) {
        stop("fit model first")
      }
      return(list(obs = self$obs()$vcomp(),
                  hidden = self$hidden()$vcomp()))
    },
    
    #' @description Model parameters
    #' @param t returns parameters at time t, default is t = 1 
    #' 
    #' @return A list with elements:
    #' \itemize{
    #'   \item{\code{obspar}}{Parameters of observation model}
    #'   \item{\code{tpm}}{Transition probability matrix of hidden state model}
    #' }
    par = function(t = 1) {
      obspar <- self$obs()$par(t = t)
      tpm <- self$hidden()$tpm(t = t)
      return(list(obspar = obspar, tpm = tpm))
    },
    
    #' @description Suggest initial parameters for the model based on kmeans
    #' clustering 
    #' 
    #' @return list of initial parameters for each observation variable 
    suggest_initial = function() {
      # do clustering
      cluster <- kmeans(self$obs()$obs_var(expand = TRUE), 
                        centers = self$hidden()$nstates())
      states <- cluster$cluster
      # get current parameters 
      current_par <- self$obs()$par(t = 1, full_names = FALSE)[,,1]
      # handle one parameter case
      if (is.null(dim(current_par))) {
        current_par <- matrix(current_par, nc = length(current_par))
      }
      par_count <- 1 
      # initial observation parameters 
      par <- vector(mode = "list", length = ncol(self$obs()$obs_var()))
      names(par) <- colnames(self$obs()$obs_var())
      # loop over observed variables 
      for (i in 1:length(self$obs()$dists())) {
        var <- self$obs()$obs_var()[,i]
        # possibly pass fixed parameters to parapprox function within dist
        par_ind <- par_count:(par_count + self$obs()$dists()[[i]]$npar() - 1)
        sub_current_par <- current_par[par_ind,,drop=FALSE]
        sub_current_par <- sub_current_par[self$obs()$dists()[[i]]$fixed(),,
                                           drop=FALSE]
        npar <- self$obs()$dists()[[i]]$npar()
        subpar <- vector(mode = "list", length = npar)
        for (j in 1:self$hidden()$nstates()) {
          args <- c(list(x = var[states == j]), as.list(sub_current_par[,j]))
          approx <- do.call(self$obs()$dists()[[i]]$parapprox(), args)
          for (k in 1:npar) subpar[[k]] <- c(subpar[[k]], approx[k])
        }
        par_count <- par_count + self$obs()$dists()[[i]]$npar()
        names(subpar) <- self$obs()$dists()[[i]]$parnames()
        par[[i]] <- subpar
      }
      return(par)
    }, 
    
    #' @description Set priors for coefficients 
    #' 
    #' @param new_priors is a list of matrices for optionally 
    #' coeff_fe_obs, coeff_fe_hid, log_lambda_obs log_lambda_hid 
    #' each matrix has two rows (first row = mean, second row = sd) 
    #' specifying parameters for Normal priors 
    set_priors = function(new_priors = NULL) {
      fe <- self$coeff_fe()
      if (!is.null(new_priors$coeff_fe_obs)) {
        coeff_fe_obs_prior <- new_priors$coeff_fe_obs 
      } else {
        coeff_fe_obs_prior <- matrix(NA, nr = length(fe$obs), nc = 2)
      }
      if (!is.null(new_priors$coeff_fe_hid)) {
        coeff_fe_hid_prior <- new_priors$coeff_fe_hid 
      } else {
        coeff_fe_hid_prior <- matrix(NA, nr = length(fe$hidden), nc = 2)
      }
      lam <- self$lambda()
      if (!is.null(new_priors$log_lambda_obs)) {
        log_lambda_obs_prior <- new_priors$log_lambda_obs 
      } else {
        log_lambda_obs_prior <- matrix(NA, nr = length(lam$obs), nc = 2)
      }
      if (!is.null(new_priors$log_lambda_hid)) {
        log_lambda_hid_prior <- new_priors$log_lambda_hid
      } else {
        log_lambda_hid_prior <- matrix(NA, nr = length(lam$hidden), nc = 2)
      }
      private$priors_ <- list(coeff_fe_obs = coeff_fe_obs_prior, 
                     coeff_fe_hid = coeff_fe_hid_prior, 
                     log_lambda_obs = log_lambda_obs_prior, 
                     log_lambda_hid = log_lambda_hid_prior)
      # Setup if necessary
      if(!is.null(private$tmb_obj_)) {
        self$setup(silent = TRUE)
      }
    }, 
    
    #' @description Extract stored priors 
    priors = function() {
      return(private$priors_)
    }, 
    
    #' @description Iterations from stan MCMC fit 
    #' 
    #' @param type Either "response" for parameters on the response (natural)
    #' scale, or "raw" for parameters on the linear predictor scale.
    #' 
    #' @return see output of as.matrix in stan 
    iters = function(type = "response") {
      if (is.null(private$iters_)) {
        stop("must run mcmc before extracting iterations")
      }
      if (type == "response") {
        return(private$par_iters_)
      } else if (type == "raw") {
        return(private$iters_)
      } else {
        stop("unknown type argument given to iters()")
      }
    }, 
    
    #' @description fitted stan object from MCMC fit 
    #' 
    #' @return the stanfit object 
    stan = function() {
      if (is.null(private$iters_)) {
        stop("must run mcmc before extracting stan fitted object")
      }
      return(private$mcmc_)
    }, 
    
    #' @description Log-likelihood at current parameters
    #' 
    #' @return Log-likelihood
    llk = function() {
      # Setup if necessary
      if(is.null(private$tmb_obj_)) {
        self$setup(silent = TRUE)
      }
      # set parameter vector to current values 
      delta <- self$hidden()$delta() 
      ldelta <- log(delta[-length(delta)] / delta[length(delta)])
      par <- c(self$obs()$coeff_fe(), 
               self$obs()$lambda(), 
               self$hidden()$coeff_fe(), 
               self$hidden()$lambda(), 
               ldelta)
      # compute log-likelihood
      return(-self$tmb_obj()$fn(par))
    },
    
    #' @description Compute effective degrees of freedom 
    #' 
    #' The degrees of freedom of the fixed effects are obtained as the number of
    #' fixed effect parameters. The effective degrees of freedom of the random 
    #' effects are computed in the comp_edf private method, based on the formula 
    #' from Section 5.4.2 from Wood (2017).
    #' 
    #' @references Wood (2017). Generalized additive models: an introduction with R. 
    #' CRC press.
    edf = function() {
      # DF for fixed effects 
      df <- nrow(self$obs()$coeff_fe()) + nrow(self$hidden()$coeff_fe())
      # EDF for hidden sub-model 
      mod_mat_hid <- self$hidden()$terms() 
      S_hid_list <- mod_mat_hid$S_list
      X_hid_list <- mod_mat_hid$X_list_re
      smoopar_hid <- self$lambda()$hidden[,1]
      edf <- 0 
      k <- 1
      for (i in seq_along(S_hid_list)) {
        if(!is.null(S_hid_list[[i]])) {
          edf <- edf + private$comp_edf(X_hid_list[[i]], 
                                        S_hid_list[[i]], 
                                        smoopar_hid[k])
          k <- k + 1 
        }
      }
      # EDF for observation sub-model 
      mod_mat_obs <- self$obs()$terms()
      S_obs_list <- mod_mat_obs$S_list
      X_obs_list <- mod_mat_obs$X_list_re
      smoopar_obs <- self$lambda()$obs[,1]
      k <- 1 
      for (i in seq_along(S_obs_list)) {
        if(!is.null(S_obs_list[[i]])) {
          edf <- edf + private$comp_edf(X_obs_list[[i]], 
                                        S_obs_list[[i]], 
                                        smoopar_obs[k])
          k <- k + 1 
        }
      }
      return(df + edf)
    }, 
    

    # Model fitting -----------------------------------------------------------

    #' @description TMB setup
    #' 
    #' This creates an attribute \code{tmb_obj}, which can be used to 
    #' evaluate the negative log-likelihood function.
    #' 
    #' @param silent Logical. If TRUE, all tracing outputs are hidden (default).
    setup = function(silent = TRUE) {
      # Vector of codes of observation distributions
      distcode <- as.vector(sapply(self$obs()$dists(), function(d) d$code()))
      # Vector of number of parameters for observation distributions
      distpar <- as.vector(sapply(self$obs()$dists(), function(d) d$npar()))
      
      # check distcodes exist
      if (any(is.null(distcode))) {
        stop("not all observation distributions are properly created. If you added any custom distributions
             recently, make sure you re-create them using Dist$new after you have recompiled the package.")
      }
      
      # Number of states
      n_states <- self$hidden()$nstates()
      
      # Create model matrices of observation process
      # (Design matrices for fixed and random effects, and smoothing matrix)
      mod_mat_obs <- self$obs()$terms()
      X_fe_obs <- mod_mat_obs$X_fe
      X_re_obs <- mod_mat_obs$X_re
      S_obs <- mod_mat_obs$S
      ncol_re_obs <- mod_mat_obs$ncol_re
      
      # Create model matrices of hidden state process
      # (Design matrices for fixed and random effects, and smoothing matrix)
      mod_mat_hid <- self$hidden()$terms()
      X_fe_hid <- mod_mat_hid$X_fe
      X_re_hid <- mod_mat_hid$X_re
      S_hid <- mod_mat_hid$S
      ncol_re_hid <- mod_mat_hid$ncol_re
      
      # Prepare delta initial parameter
      delta <- self$hidden()$delta() 
      ldelta <- log(delta[-length(delta)] / delta[length(delta)])
      
      # Setup TMB parameters
      tmb_par <- list(coeff_fe_obs = self$obs()$coeff_fe(),
                      log_lambda_obs = 0,
                      coeff_fe_hid = self$hidden()$coeff_fe(),
                      log_lambda_hid = 0,
                      log_delta = ldelta,
                      coeff_re_obs = 0,
                      coeff_re_hid = 0)
      
      # Initialise vectors of fixed parameters and random effects
      map <- NULL
      random <- NULL
      
      # Setup random effects in observation model
      if(is.null(S_obs)) {
        # If there are no random effects, 
        # coeff_re and log_lambda are not estimated
        map <- c(map, list(coeff_re_obs = factor(NA),
                           log_lambda_obs = factor(NA)))
        S_obs <- as(matrix(0, 1, 1), "sparseMatrix")
        ncol_re_obs <- 0
        X_re_obs <- as(rep(0, nrow(X_fe_obs)), "sparseMatrix")
      } else {
        # If there are random effects, 
        # set initial values for coeff_re and log_lambda
        random <- c(random, "coeff_re_obs")
        tmb_par$coeff_re_obs <- rep(0, ncol(S_obs))
        tmb_par$log_lambda_obs <- rep(0, length(ncol_re_obs))
      }
      
      # check if delta to be stationary (overrides custom map specified)
      statdist <- 0 
      if (self$hidden()$stationary()) {
        private$fixpar_$log_delta <- rep(NA, self$hidden()$nstates())
        statdist <- 1 
      }
      
      # check for parameters that must always be fixed (identified by having
      # a fixed = TRUE in dist)
      fixed <- unlist(lapply(self$obs()$dists(), function(d) d$fixed()))
      if (any(fixed)) {
        nms <- names(fixed)[fixed == TRUE]
        obsnms <- rownames(self$obs()$coeff_fe())
        for (i in 1:length(nms)) {
          getnms <- obsnms[grep(nms[i], obsnms)]
          oldnms <- names(private$fixpar_$obs)
          private$fixpar_$obs <- c(private$fixpar_$obs, rep(NA, length(getnms)))
          names(private$fixpar_$obs) <- c(oldnms, getnms)
        }
      }
      
      # check for transitions that have fixed probabilities 
      # find transitions that are fixed 
      ls_form_char <- as.list(t(self$hidden()$structure())[!diag(self$hidden()$nstates())])
      which_fixed <- sapply(ls_form_char, function(x) {x == "."})
      getnms <- rownames(self$hidden()$coeff_fe())[which_fixed]
      oldnms <- names(private$fixpar_$hid)
      private$fixpar_$hid <- c(private$fixpar_$hid, rep(NA, length(getnms)))
      names(private$fixpar_$hid) <- c(oldnms, getnms)
      
      # add custom mapping effects
      usernms <- c("obs", "hid", "lambda", "delta")
      comps <- c("coeff_fe_obs", "coeff_fe_hid", "log_lambda", "log_delta")
      for (i in seq_along(usernms)) {
        if (!is.null(private$fixpar_[[usernms[i]]])) {
          v <- tmb_par[[comps[i]]]
          tmp <- 1:length(v)
          if (is.matrix(v)) {
            nms <- rownames(v)
          } else {
            nms <- names(v)
          }
          tmp[nms %in% names(private$fixpar_[[usernms[i]]])] <- 
            as.numeric(private$fixpar_[[usernms[i]]])
          tmp <- factor(as.vector(tmp))
          ls <- list(tmp)
          names(ls) <- comps[i]
          map <- c(map, ls)
        }
      }

      # Setup random effects in hidden state model
      if(is.null(S_hid)) {
        # If there are no random effects, 
        # coeff_re and log_lambda are not estimated
        map <- c(map, list(coeff_re_hid = factor(NA),
                           log_lambda_hid = factor(NA)))
        S_hid <- as(matrix(0, 1, 1), "sparseMatrix")
        ncol_re_hid <- 0
        X_re_hid <- as(rep(0, nrow(X_fe_hid)), "sparseMatrix")
      } else {
        # If there are random effects, 
        # set initial values for coeff_re and log_lambda
        random <- c(random, "coeff_re_hid")
        tmb_par$coeff_re_hid <- rep(0, ncol(S_hid))
        tmb_par$log_lambda_hid <- rep(0, length(ncol_re_hid))
      }
      
      # Get stored priors 
      priors <- self$priors() 
      
      # Get variables for observation distributions
      obsvar <- self$obs()$obs_var(expand = TRUE)
      datadim <- attr(obsvar, "datadim")
      
      # Data for TMB
      tmb_dat <- list(ID = self$obs()$data()$ID,
                      data = as.matrix(obsvar),
                      datadim = datadim, 
                      known_states = as.vector(self$obs()$known_states()) - 1, 
                      n_states = n_states,
                      statdist = statdist, 
                      distcode = distcode,
                      distpar = distpar, 
                      X_fe_obs = X_fe_obs,
                      X_re_obs = X_re_obs,
                      S_obs = S_obs,
                      ncol_re_obs = ncol_re_obs,
                      X_fe_hid = X_fe_hid,
                      X_re_hid = X_re_hid,
                      S_hid = S_hid,
                      ncol_re_hid = ncol_re_hid, 
                      coeff_fe_obs_prior = priors$coeff_fe_obs, 
                      coeff_fe_hid_prior = priors$coeff_fe_hid, 
                      log_lambda_obs_prior = priors$log_lambda_obs, 
                      log_lambda_hid_prior = priors$log_lambda_hid)
      
      # Create TMB model
      obj <- MakeADFun(tmb_dat, tmb_par, dll = "hmmTMB", 
                       random = random,
                       map = map, 
                       silent = silent)
      
      nllk0 <- obj$fn(obj$par)
      if(is.nan(nllk0) | is.infinite(nllk0)) {
        stop(paste("log-likelihood is NaN or infinite at starting parameters.",
                   "Check that the data are within the domain of definition of the",
                   "observation distributions, and/or try other starting parameters."))
      }
      
      # Negative log-likelihood function
      private$tmb_obj_ <- obj
    },
    
    #' @description Fit model using tmbstan
    #' 
    #' @param ... Arguments passed to tmbstan
    #' @param silent Logical. If FALSE, all tracing outputs are shown (default).
    mcmc = function(..., silent = FALSE) {
      self$formulation()
      
      if (!requireNamespace("rstan", quietly = TRUE)) {
        stop("you need to install the package rstan to do this")
      }
      
      if (!requireNamespace("tmbstan", quietly = TRUE)) {
        stop("you need to install the package tmbstan to do this")
      }
      
      # Setup if necessary
      if(is.null(private$tmb_obj_)) {
        self$setup(silent = silent)
      }
      
      # do stan MCMC 
      private$mcmc_ <- tmbstan::tmbstan(private$tmb_obj_, ...)
      
      # store iterations 
      private$iters_ <- as.matrix(private$mcmc_)
      
      # store iterations on response scale 
      npar <- length(unlist(self$par()))
      niter <- nrow(private$iters_)
      par_iters <- matrix(0, nr = niter, nc = npar)
      for (i in 1:niter) {
        self$update_par(iter = i)
        par_iters[i,] <- unlist(self$par())
      }
      private$par_iters_ <- par_iters 
      
      # set coefficients to posterior means 
      self$update_par(iter = "mean")
      
    }, 
    
    #' @description Model fitting
    #' 
    #' The negative log-likelihood of the model is minimised using the
    #' function \code{optim}. TMB uses the Laplace approximation to integrate 
    #' the random effects out of the likelihood.
    #' 
    #' After the model has been fitted, the output of \code{optim} can be
    #' accessed using the method \code{out}.
    #' 
    #' @param silent Logical. If FALSE, all tracing outputs are hidden (default).
    #' @param ... other arguments to optimx which is used to optimise likelihood, 
    #' see ?optimx
    fit = function(silent = FALSE, ...) {
      if(!silent) self$formulation()
      
      # Setup if necessary
      if(is.null(private$tmb_obj_)) {
        self$setup(silent = silent)
      }
      
      # Number of states
      n_states <- self$hidden()$nstates()
      
      # Fit model
      args <- list(...)
      if (any(c("par", "fn", "gr", "he", "hessian") %in% names(args))) {
        stop("cannot supply arguments to fit of name par, fn, gr, he, or hessian. These are handled by TMB")
      }
      if (any(c("lower", "upper") %in% names(args))) {
        warning("lower and upper arguments to optimx are ignored in hmmTMB")
      }
      # change default method to nlminb
      private$tmb_obj_$method <- "nlminb"
      if ("method" %in% names(args)) {
        private$tmb_obj_$method <- args$method 
        args <- args[which(names(args) != "method")]
      }
      # create temporary optimization function
      opt_fn <- function(par) {as.vector(args$fn(par))}
      opt_gr <- function(par) {as.vector(args$gr(par))}
      # fit model 
      args <- c(private$tmb_obj_, args)
      private$out_ <- optimx(par = args$par, 
                             fn = opt_fn, 
                             gr = opt_gr, 
                             method = args$method,
                             itnmax = args$itnmax, 
                             hessian = args$hessian, 
                             control = args$control)
      best <- which.min(private$out_$value)
      if (private$out_$convcode[best] != 0) {
        warning("convergence code was not zero, indicating that the optimizer may
                not have converged to the correct estimates. Please check by consulting
                the out() function which shows what optimx returned.")
      }
      
      # Get estimates and precision matrix for all parameters
      best_par <- as.vector(private$out_[best, 1:length(private$tmb_obj_$par)])
      rownames(best_par) <- NULL
      names(best_par) <- names(private$tmb_obj_$par)
      private$tmb_obj_$par <- best_par
      private$tmb_rep_ <- sdreport(private$tmb_obj_,
                                   getJointPrecision = TRUE, 
                                   skip.delta.method = FALSE)
      par_list <- as.list(private$tmb_rep_, "Estimate")
      
      # update parameters 
      self$update_par(par_list)
      
    },
    
    #' @description Get maximum likelihood estimates once model fitted
    #' 
    #' @return list of maximum likelihood estimates as described as
    #' input for the function update_par() 
    mle = function() {
      if (is.null(private$out_)) stop("must fit model with fit() function first")
      par_list <- as.list(self$tmb_rep(), "Estimate")
      return(par_list)
    }, 
    

    # Forward-backward probabilities ------------------------------------------
    
    #' @description Forward-backward algorithm 
    #' 
    #' @return log-forward and log-backward probabilities 
    forward_backward = function() {
      delta <- self$hidden()$delta() 
      n <- nrow(self$obs()$data())
      lforw <- lback <- matrix(0, nr = self$hidden()$nstates(), nc = n)
      # get observation probabilities 
      obsprobs <- self$obs()$obs_probs()
      # get tpms 
      tpms <- self$hidden()$tpm(t = "all")
      # forward algorithm 
      p <- delta * obsprobs[1,]
      psum <- sum(p)
      llk <- log(psum)
      p <- p / psum
      lforw[, 1] <- log(p) + llk
      for (i in 2:n) {
        p <- p %*% tpms[, , i] * obsprobs[i, ]
        psum <- sum(p)
        llk <- llk + log(psum)
        p <- p / psum
        lforw[, i] <- log(p) + llk 
      }
      # backward algorithm
      lback[, n] <- rep(0, self$hidden()$nstates())
      p <- rep(1 / self$hidden()$nstates(), self$hidden()$nstates())
      llk <- log(self$hidden()$nstates())
      for (i in (n - 1):1) {
        p <- tpms[, , i] %*% (obsprobs[i + 1, ] * p)
        lback[, i] <- log(p) + llk
        psum <- sum(p)
        p <- p / psum
        llk <- llk + log(psum)
      }
      return(list(logforward = lforw, logbackward = lback))
    }, 
  
    # Conditional distributions -----------------------------------------------

    #' @description Compute conditional cumulative distribution functions 
    #' 
    #' @param ngrid how many cells on the grid that CDF is computed on 
    #' @param silent if TRUE then no messages are printed 
    #' 
    #' @return cdfs on grid for each variable 
    cond = function(ngrid = 100, silent = FALSE) {
      delta <- t(self$hidden()$delta())
      vars <- self$obs()$obs_var()
      nvars <- ncol(vars)
      n <- nrow(self$obs()$data())
      range <- matrix(0, nr = 2, nc = nvars)
      range[1,] <- sapply(vars, min, na.rm = TRUE)  
      range[2,] <- sapply(vars, max, na.rm = TRUE)
      # get forward-backward probabilities
      fb <- self$forward_backward() 
      lforw <- fb$logforward
      lback <- fb$logbackward
      lforw <- cbind(log(delta), lforw)
      # get transition matrices
      tpms <- self$hidden()$tpm(t = "all")
      # scaling 
      forwscale <- apply(lforw, 2, max)
      backscale <- apply(lback, 2, max)
      # compute conditional state probabilities
      if (!silent) cat("Computing conditional state probabilities...")
      cond <- matrix(0, nr = n, nc = self$hidden()$nstates()) 
      for (i in 1:n) {
        p <- (exp(lforw[,i] - forwscale[i]) %*% tpms[, , i]) * exp(lback[, i] - backscale[i])
        cond[i, ] <- p / sum(p)
      }
      if(!silent) cat("done\n")
      # list to store cdfs 
      pdfs <- array(0, dim = c(nvars, n, ngrid))
      grids <- vector(mode = "list", length = nvars)
      # compute cdf for each variable 
      obsmats <- self$obs()$terms()
      varnms <- names(vars)
      for (i in 1:nvars) {
        if (!silent) cat("Computing CDF for", varnms[i], "...")
        grid <- seq(range[1, i], range[2, i], length = ngrid)
        if (is_whole_number(vars[,i])) {
          grid <- unique(floor(grid))
        }
        grids[[i]] <- grid 
        for (g in 1:length(grid)) {
          tmp <- data.frame(var = rep(grid[g], n))
          colnames(tmp) <- varnms[i]
          probs <- self$obs()$obs_probs(data = tmp)
          pdfs[i, , g] <- rowSums(probs * cond)
        }
        if (!silent) cat("done\n")
      }
      return(list(grids = grids, pdfs = pdfs))
    }, 
    
    # Pseudo-residuals --------------------------------------------------------

    #' @description Pseudo-residuals
    #' 
    #' @details Compute pseudo-residuals for the fitted model. If the fitted model
    #' is the "true" model, the pseudo-residuals follow a standard normal distribution.
    #' Deviations from normality suggest lack of fit.
    #' 
    #' @return Matrix of pseudo-residuals, with one row for each response variable
    #' and one column for each observation
    pseudores = function() {
      n <- nrow(self$obs()$data())  
      cond <- self$cond(silent = TRUE)
      pdfs <- cond$pdfs
      grids <- cond$grids
      vars <- self$obs()$obs_var()
      nvars <- ncol(vars)
      varnms <- names(vars)
      # sum CDFs cumulatively 
      cdfs <- array(0, dim = dim(pdfs))
      for (v in 1:nvars) {
        cdfs[v,,] <- t(apply(pdfs[v,,], 1, cumsum))
        # if continuous then approximate the integral using Riemann sum
        if (!is_whole_number(vars[,v])) {
          dgrid <- diff(grids[[v]])[1]
          cdfs[v,,] <- cdfs[v,,] * dgrid 
        }
      }
      # do residuals for each variable 
      r <- matrix(0, nr = nvars, nc = n)
      for (v in 1:nvars) {
        for (i in 1:length(vars[,v])) {
          wh <- which.min(abs(vars[i, v] - grids[[v]]))
          r[v, i] <- qnorm(cdfs[v, i, wh])
        }
        rownames(r)[v] <- varnms[v]
      }
      return(r)
    }, 
  
    # State estimation --------------------------------------------------------
    
    #' @description Viterbi algorithm
    #' 
    #' @return Most likely state sequence
    viterbi = function() {
      data <- self$obs()$data()
      ID <- data$ID
      
      # Number of observations
      n <- nrow(data)
      # Number of states
      n_states <- self$hidden()$nstates()
      # Number of variables
      n_var <- length(self$obs()$dists())
      
      # Observation probabilities
      obs_probs <- self$obs()$obs_probs()
      
      # Transition probability matrices      
      tpm_all <- self$hidden()$tpm(t = "all")
      
      # Number of unique IDs
      n_id <- length(unique(ID))
      # First index for each ID
      i0 <- c(1, which(ID[-1] != ID[-n]) + 1, n + 1)
      
      # Initialise state sequence
      all_states <- NULL
      
      # Initial distribution
      delta <- self$hidden()$delta() 
      
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
    },
    
    #' @description Sample posterior state sequences using forward-filtering
    #' backward-sampling 
    #' 
    #' @param nsamp number of samples to produce 
    #' @param full if TRUE and model fit by mcmc then parameter estimates are 
    #' sampled from the posterior samples before simulating each sequence 
    #' 
    #' @return matrix where each column is a different sample of state sequences 
    sample_states = function(nsamp = 1, full = FALSE) {
      if (full & is.null(private$mcmc_)) stop("must fit model with $mcmc first")
      if (full) {
        n_full_samp <- nsamp
        nsamp <- 1 
      } else {
        n_full_samp <- 1
      }
      nstates <- self$hidden()$nstates()
      actual_samps <- max(n_full_samp, nsamp)
      n <- nrow(self$obs()$data())
      states <- matrix(0, nr = n, nc = actual_samps)
      for (k in 1:n_full_samp) {
        # sample a parameter iteration at random, if FULL asked for 
        if (full) self$update_par(iter = sample(1:nrow(mod$iters()), size = 1))
        # get forward-backward probabilities 
        fb <- self$forward_backward()
        # sample last states
        L <- log(sum(exp(fb$logforward[,n] - max(fb$logforward[,n])))) + max(fb$logforward[,n])
        prob <- exp(fb$logforward[,n] - L)
        prob <- prob / sum(prob)
        if (full) {
          states[n,k] <- sample(1:nstates, prob = prob, size = nsamp, replace = TRUE)
        } else {
          states[n,] <- sample(1:nstates, prob = prob, size = nsamp, replace = TRUE)
        }
        # get tpms 
        tpms <- self$hidden()$tpm(t = "all")
        # get observation probabilties
        obsprobs <- self$obs()$obs_probs()
        # sample backward
        for (s in 1:nsamp) {
          for (i in (n - 1):1) {
            prob <- exp(fb$logforward[, i] + log(tpms[, states[i + 1, s], i + 1]) + 
              log(obsprobs[i + 1, states[i + 1, s]]) + fb$logbackward[states[i + 1, s], i + 1] - L) 
            prob <- prob / sum(prob)
            if (full) {
              ind <- k
            } else {
              ind <- s 
            }
            states[i, ind] <- sample(1:nstates, prob = prob, size = 1)
          }
        }
      }
      return(states)
    }, 
    
    #' @description Compute posterior probability of being in each state 
    #' 
    #' @return matrix with a row for each observation and a column for each state 
    state_probs = function() {
      delta <- self$hidden()$delta
      n <- nrow(self$obs()$data())
      nstates <- self$hidden()$nstates()
      fb <- self$forward_backward()
      llk <- log(sum(exp(fb$logforward[,n] - max(fb$logforward[,n])))) + max(fb$logforward[,n])
      pr_state <- matrix(0, nr = n, nc = nstates)
      for (i in 1:n) {
        pr_state[i,] <- exp(fb$logforward[,i] + fb$logbackward[,i] - llk)
      }
      return(pr_state)
    }, 
    
    # Prediction and Uncertainty ----------------------------------------------
    
    #' @description Posterior sampling for model coefficients
    #' 
    #' @param n_post Number of posterior samples
    #' 
    #' @return Matrix with one column for each coefficient and one row
    #' for each posterior draw
    post_coeff = function(n_post) {
      # Get parameter estimates and covariance matrix
      rep <- self$tmb_rep()
      if(is.null(rep$jointPrecision)) {
        par <- rep$par.fixed
        V <- rep$cov.fixed
      } else {
        par <- c(rep$par.fixed, rep$par.random)
        V <- solve(rep$jointPrecision)
      }
      
      # Generate samples from MVN estimator distribution
      post <- rmvn(n = n_post, mu = par, V = V)
      
      # ensure it is a matrix
      post <- matrix(post, nc = length(par), nr = n_post)
      
      # parameter names
      colnames(post) <- names(par)
      
      return(post)
    },
    
    #' @description Posterior sampling for linear predictor 
    #' 
    #' @param n_post Number of posterior samples
    #' 
    #' @return Matrix with one column for each predictor and one row
    #' for each posterior draw
    post_linpred = function(n_post) {
      # save current parameters 
      coeff_fe_old <- self$coeff_fe() 
      coeff_re_old <- self$coeff_re()
      pars <- self$post_coeff(n_post)
      # observation submodel 
      obspars_fe <- pars[,which(colnames(pars) == "coeff_fe_obs")]
      hidpars_fe <- pars[,which(colnames(pars) == "coeff_fe_hid")]
      obspars_re <- pars[,which(colnames(pars) == "coeff_re_obs")]
      hidpars_re <- pars[,which(colnames(pars) == "coeff_re_hid")]
    
      lp <- NULL
      lp$obs <- matrix(0, nr = n_post, nc = nrow(self$obs()$X_fe()))
      lp$hidden <- matrix(0, nr = n_post, nc = nrow(self$hidden()$X_fe()))
      for (i in 1:n_post) {
        self$obs()$update_coeff_fe(obspars_fe[i,])
        self$obs()$update_coeff_re(obspars_re[i,])
        self$hidden()$update_coeff_fe(hidpars_fe[i,])
        self$hidden()$update_coeff_re(hidpars_re[i,])
        lp$obs[i,] <- self$obs()$linpred()
        lp$hidden[i,] <- self$hidden()$linpred()
      }
      
      # reset design matrices and parameters
      self$obs()$update_coeff_fe(coeff_fe_old$obs)
      self$obs()$update_coeff_re(coeff_re_old$obs)
      self$hidden()$update_coeff_fe(coeff_fe_old$hidden)
      self$hidden()$update_coeff_re(coeff_re_old$hidden)
      
      return(lp)
    }, 
    
    #' @description Create posterior simulations of a function of a model component 
    #' 
    #' @param fn the function which takes a vector of linear predictors as input
    #'           and produces either a scalar or vector output 
    #' @param comp is "obs" for observation model linear predictor, "hidden" for
    #'             hidden model linear predictor 
    #' @param n_post number of posterior simulations 
    #' @param ... arguments passed to fn
    #' @param level confidence interval level, default is 95\% 
    #' 
    #' @return a vector (for scalar outputs of fn) or a matrix (for vector outputs)
    #' with a column for each simulation 
    post_fn = function(fn, n_post, comp = NULL, ..., level = 0) {
      # get linear predictors
      lp <- self$post_linpred(n_post)
      # output
      res <- NULL
      # compute function of linear predictor
      res$samp <- lapply(1:nrow(lp[[comp]]), FUN = function(i) {fn(linpred = lp[[comp]][i,], ...)})
      # compute means 
      res$mean <- Reduce("+", res$samp) / length(res$samp)
      # compute confidence interval
      if (level > 0) {
        alp <- (1 - level) / 2
        arr <- simplify2array(res$samp)
        ci <- apply(arr, 1:(length(dim(arr)) - 1), quantile, prob = c(0.025, 0.975))
        nci <- length(dim(ci))
        ci <- aperm(ci, c(2:nci, 1))
        block <- prod(dim(ci)[-nci])
        lcl <- ci[1:block]
        ucl <- ci[(1:block) + block]
        dim(lcl) <- dim(ucl) <- dim(ci)[-nci]
        res$lcl <- lcl
        res$ucl <- ucl
      }
      return(res)
    }, 
    
    #' @description Predict estimates from a fitted model
    #' 
    #' @param name which estimates to predict? Options include 
    #' transition probability matrices "tpm", 
    #' stationary distributions "delta", or 
    #' observation distribution parameters "obspar"
    #' @param t time points to predict at 
    #' @param ... other arguments to the respective functions for hidden$tpm, hidden$delta, obs$par
    #' @param newdata new dataframe to use for prediction
    #' @param level if greater than zero, then produce confidence intervals with this level, e.g. CI = 0.95
    #'           will produce 95\% confidence intervals 
    #' @param n_post if greater than zero then n_post posterior samples are produced 
    #' @return named array of predictions and confidence interval, if requested
    predict = function(name, t = 1, newdata = NULL, level = 0, n_post = 0) {
      if (is.null(private$out_) & (level > 0 | n_post > 0)) stop("must fit model with fit() function first")
      if (!is.null(newdata)) {
        old <- list(X_fe_obs = self$obs()$X_fe(), 
                    X_re_obs = self$obs()$X_re(), 
                    X_fe_hid = self$hidden()$X_fe(), 
                    X_re_hid = self$hidden()$X_re())
        obsmats <- self$obs()$make_mat(new_data = newdata)
        hidmats <- self$hidden()$make_mat(data = self$obs()$data(), 
                                          new_data = newdata)
        self$obs()$update_X_fe(obsmats$X_fe)
        self$obs()$update_X_re(obsmats$X_re)
        self$hidden()$update_X_fe(hidmats$X_fe) 
        self$hidden()$update_X_re(hidmats$X_re)
      }
      
      # get appropriate prediction function 
      fn <- switch(name, 
                   tpm = self$hidden()$tpm,
                   delta = self$hidden()$delta,
                   obspar = self$obs()$par)
      
      # get appropriate model component 
      comp <- switch(name, tpm = "hidden", delta = "hidden", obspar = "obs")
      
      # just return predicted means if no confidence intervals wanted 
      # or posterior simulations 
      if (level == 0 & n_post == 0) {
        val <- fn(linpred = self[[comp]]()$linpred(), t = t)
      }
      
      # do posterior sampling if asked for  
      if (n_post > 0) {
        
        sim <- self$post_fn(fn = fn, 
                            n_post = n_post, 
                            comp = comp, 
                            t = t,
                            level = level)
        val <- sim
      }
      
      # delta method for CI 
      if (level > 0 & n_post < 1e-10) {
        oldpar <- list(obs_coeff_fe = self$obs()$coeff_fe(), 
                       obs_coeff_re = self$obs()$coeff_re(), 
                       hid_coeff_fe = self$hidden()$coeff_fe(), 
                       hid_coeff_re = self$hidden()$coeff_re())
        # get variance-covariance matrix of linear predictor 
        lfn <- function(par, oldpar, ind_fe, ind_re, comp) {
          new_par <- oldpar[[paste0(comp, "_coeff_fe")]]
          new_par[!(rownames(new_par) %in% names(private$fixpar_[[comp]]))] <- par[ind_fe]
          self[[comp]]()$update_coeff_fe(new_par)
          if(!is.null(ind_re)) self[[comp]]()$update_coeff_re(par[ind_re])
          return(self[[comp]]()$linpred())
        }
        # get variance-covariance of par 
        if(is.null(self$tmb_rep()$jointPrecision)) {
          par <- self$tmb_rep()$par.fixed
          V <- self$tmb_rep()$cov.fixed
        } else {
          par <- c(self$tmb_rep()$par.fixed, self$tmb_rep()$par.random)
          V <- solve(self$tmb_rep()$jointPrecision)
        }
        ind_fe <- grepl(paste0("coeff_fe_", substr(comp, 1, 3)), names(par))
        ind_re <- grepl(paste0("coeff_re_", substr(comp, 1, 3)), names(par))
        # compute jacobian 
        J <- numDeriv::jacobian(lfn, 
                                par, 
                                oldpar = oldpar, 
                                ind_fe = ind_fe, 
                                ind_re = ind_re, 
                                comp = comp)
        # compute variance-covariance for linear predictor
        Vlinpred <- J %*% V %*% t(J)
        sds <- sqrt(diag(Vlinpred))
        # compute confidence bounds on linear predictor scale
        mu <- self[[comp]]()$linpred()
        z <- qnorm(1 - (1 - level) / 2)
        lcl <- mu - z * sds 
        ucl <- mu + z * sds 
        # get results on response scale
        mu <- fn(linpred = mu, t = t)
        lcl <- fn(linpred = lcl, t = t)
        ucl <- fn(linpred = ucl, t = t)
        # reset parameters 
        self$obs()$update_coeff_fe(oldpar$obs_coeff_fe)
        self$obs()$update_coeff_re(oldpar$obs_coeff_re)
        self$hidden()$update_coeff_fe(oldpar$hid_coeff_fe)
        self$hidden()$update_coeff_re(oldpar$hid_coeff_re)
        val <- list(mean = mu, lcl = lcl, ucl = ucl)
      }
      
      if (!is.null(newdata)) {
        self$obs()$update_X_fe(old$X_fe_obs)
        self$hidden()$update_X_fe(old$X_fe_hid)
        self$obs()$update_X_re(old$X_re_obs)
        self$hidden()$update_X_re(old$X_re_hid)
      }
      
      return(val)
      
    }, 
    
    # Simulation --------------------------------------------------------------
    
    #' @description Simulate from hidden Markov model
    #' 
    #' @param n Number of time steps to simulate
    #' @param data Optional data frame including covariates
    #' @param silent if TRUE then no messages are printed
    #' 
    #' @return Data frame including columns of data (if provided), and simulated
    #' data variables
    simulate = function(n, data = NULL, silent = FALSE) {
      if(is.null(data)) {
        data <- data.frame(ID = rep(factor(1), n))
      } else if(is.null(data$ID)) {
        data$ID <- rep(factor(1), n)
      }
      
      # Number of states
      n_states <- self$hidden()$nstates()
      
      # Simulate state process      
      S <- self$hidden()$simulate(n = n, data = self$obs()$data(), 
                                  new_data = data, 
                                  silent = silent)
      
      # Create observation parameters
      mats_obs <- self$obs()$make_mat(new_data = data)
      X_fe_old <- self$obs()$X_fe()
      X_re_old <- self$obs()$X_re()
      self$obs()$update_X_fe(mats_obs$X_fe)
      self$obs()$update_X_re(mats_obs$X_re)
      par_obs <- self$obs()$par(t = "all", full_names = FALSE)
      
      # Simulate observation process
      obs_dists <- self$obs()$dists()
      obs_all <- data
      attributes(obs_all)$state <- S
      par_count <- 1
      for(var in seq_along(obs_dists)) {
        # Distribution, name, and parameter indices for this variable
        obsdist <- obs_dists[[var]]
        var_name <- names(obs_dists)[var]
        par_ind <- par_count:(par_count + obsdist$npar() - 1)
        
        # Simulate n realisations for variable "var"
        obs <- rep(NA, n)
        for(i in 1:n) {
          if(round(i/n*100)%%10 == 0) {
            if(!silent) cat("\rSimulating ", var_name, "... ", round(i/n*100), "%", sep = "")        
          }
          
          # Generate realisation
          obs[i] <- obsdist$rng_apply(n = 1, par = par_obs[par_ind, S[i], i])
        }
        
        # Add variable to data frame
        obs_all[[var_name]] <- obs
        par_count <- par_count + obsdist$npar()
        if(!silent) cat("\n")
      }
      
      # reset design matrices
      self$obs()$update_X_fe(X_fe_old)
      self$obs()$update_X_re(X_re_old)
      
      return(obs_all)
    },
    
    #' @description Compute goodness-of-fit statistics using simulation
    #' 
    #' @param gof_fn goodness-of-fit function which accepts "data" as input
    #'               and returns a statistic: either a vector or a single number. 
    #' @param nsims number of simulations to perform 
    #' @param full if model fit with MCMC then full set to TRUE will sample from
    #' posterior for each simulation 
    #' @param silent Logical. If FALSE, simulation progress is shown. 
    #' (Default: TRUE)
    #' 
    #' @return List with elements:
    #' \itemize{
    #'   \item{obs_stat}{Vector of values of goodness-of-fit statistics for the
    #'   observed data}
    #'   \item{stats}{Matrix of values of goodness-of-fit statistics for the
    #'   simulated data sets (one row for each statistic, and one column for each
    #'   simulation)}
    #'   \item{plot}{ggplot object}
    #' }
    gof = function(gof_fn, nsims = 100, full = FALSE, silent = FALSE) {
      # Evaluate statistics for observed data
      obs_stat <- gof_fn(self$obs()$data())
      
      # Simulate from model and evaluate statistics for simulated data
      stats <- matrix(0, nc = nsims, nr = length(obs_stat))
      for (sim in 1:nsims) {
        if (!silent) cat("Simulating", sim, " / ", nsims, "\r")
        
        # if full and mcmc then sample parameter
        if (full & !is.null(private$mcmc_)) {
          self$update_par(iter = sample(1:nrow(self$iters()), size = 1))
        }
        
        # simulate new data
        newdat <- self$simulate(n = nrow(self$obs()$data()), 
                                data = self$obs()$data(),
                                silent = TRUE) 
        # compute statistics
        stats[,sim] <- gof_fn(newdat)
      }
      
      # Get names of statistics
      if(is.null(names(obs_stat))) {
        stat_names <- paste("statistic", 1:length(obs_stat))
        names(obs_stat) <- stat_names
      } else {
        stat_names <- names(obs_stat)
      }
      rownames(stats) <- stat_names
      
      # Data frames for ggplot (observed and simulated values)
      df_obs <- as.data.frame.table(as.matrix(obs_stat))
      colnames(df_obs) <- c("stat", "sim", "val")
      df <- as.data.frame.table(stats)
      colnames(df) <- c("stat", "sim", "val")
      
      # Create plot
      p <- ggplot(df, aes(val)) + 
        geom_histogram(bins = 20, aes(y=..density..), col = "white", 
                       bg = "lightgrey", na.rm = TRUE) +
        facet_wrap("stat", scales = "free") + 
        geom_vline(aes(xintercept = val), data = df_obs) +
        xlab("statistic") + ggtitle("Vertical line is observed value") +
        theme_light() +
        theme(strip.background = element_blank(),
              strip.text = element_text(colour = "black"))
      plot(p)
      
      return(list(obs_stat = obs_stat, stats = stats, plot = p))
    }, 
    

    # Plotting ----------------------------------------------------------------
    
    #' @description Time series plot coloured by states
    #' 
    #' Creates a plot of the data coloured by the most likely state sequence,
    #' as estimated by the Viterbi algorithm. If one variable name is passed
    #' as input, it is plotted against time. If two variables are passed, they
    #' are plotted against each other. 
    #' 
    #' @param var1 Name of the variable to plot.
    #' @param var2 Optional name of a second variable, for 2-d plot.
    #' 
    #' @return A ggplot object
    plot_ts = function(var1, var2 = NULL) {
      # Data frame for plot
      data <- self$obs()$data()
      # State sequence as factor
      state <- as.factor(self$states())
      # Time series ID
      ID <- self$obs()$data()$ID
      
      if(is.null(var2)) {
        # 1d time series plot
        df <- data.frame(ID = ID,
                         index = 1:nrow(data),
                         x = data[[var1]])
        
        p <- ggplot(data = df, mapping = aes(index, x, col = state, group = ID)) +
          geom_line() + xlab("time") + ylab(var1)
      } else {
        # 2d plot
        df <- data.frame(ID = ID,
                         x = data[[var1]],
                         y = data[[var2]])
        
        p <- ggplot(data = df, mapping = aes(x, y, col = state, group = ID)) +
          geom_path() + xlab(var1) + ylab(var2)
      }
      
      p <- p + 
        scale_color_manual(values = hmmTMB_cols) +
        theme_light()
      
      return(p)
    },
    
    #' @description Plot a model component 
    #' 
    #' @param name name of model component: tpm, delta, or obspar 
    #' @param var name of covariate to plot on x-axis 
    #' @param covs Optional named list for values of covariates (other than 'var') 
    #' that should be used in the plot (or dataframe with single row). If this is
    #' not specified, the mean value is used for numeric variables, and the
    #' first level for factor variables.
    #' @param i if plotting tpm then rows of tpm, if plotting delta then state, 
    #' if plotting obspar then indices of parameter to plot 
    #' @param j if plotting tpm then cols of tpm to plot, if plotting delta then
    #' ignored, if plotting obspar then indices of states to plot 
    #' @param n_grid coarseness of grid over x-axis to create 
    #' 
    #' @return A ggplot object 
    plot = function(name, var = NULL, covs = NULL, i = NULL, j = NULL, n_grid = 50) {
      # Get relevant model component 
      comp <- switch(name, tpm = "hidden", delta = "hidden", obspar = "obs")
      # Get x-axis 
      # get newdata over a grid 
      newdata <- cov_grid(var = var, 
                          obj = self, 
                          covs = covs, 
                          formulas = self[[comp]]()$formulas(), 
                          n_grid = n_grid)
      
      # Number of states
      n_states <- self$hidden()$nstates()
      
      # Get predictions and uncertainty 
      preds <- self$predict(name, t = "all", newdata = newdata, level = 0.95)
      
      # Data frame for plot
      df <- as.data.frame.table(preds$mean)
      df$low <- as.vector(preds$lcl)
      df$upp <- as.vector(preds$ucl)
      if (name == "tpm") {
        colnames(df) <- c("from", "to", "var", "prob", "low", "upp")
        levels(df$from) <- paste("State", 1:n_states)
        levels(df$to) <- paste("State", 1:n_states)
        df$var <- rep(newdata[, var], each = n_states * n_states)
        if (!is.null(i)) df <- df[df$from == paste0("State ", i),]
        if (!is.null(j)) df <- df[df$to == paste0("State ", j),]
      } else if (name == "delta") {
        colnames(df) <- c("var", "state", "prob", "low", "upp")
        levels(df$state) <- paste("State", 1:n_states)
        df$var <- rep(newdata[, var], n_states)
        if (!is.null(i)) df <- df[df$state == paste0("State ", i),]
      } else if (name == "obspar") {
        colnames(df) <- c("par", "state", "var", "val", "low", "upp")
        levels(df$state) <- paste("State", 1:n_states)
        df$var <- rep(newdata[, var], each = nrow(df)/nrow(newdata))
        if (!is.null(i)) df <- df[df$par == i,]
        if (!is.null(j)) df <- df[df$state == paste0("State ", j),]
      }
      
      # Create caption with values of other (fixed) covariates      
      plot_txt <- NULL
      if(ncol(newdata) > 1) {
        other_covs <- newdata[1, which(colnames(newdata) != var), 
                              drop = FALSE]
        
        # Round numeric values, and transform factors to strings
        num_ind <- sapply(other_covs, is.numeric)
        other_covs[num_ind] <- lapply(other_covs[num_ind], function(cov) 
          round(cov, 2))
        fac_ind <- sapply(other_covs, is.factor)
        other_covs[fac_ind] <- lapply(other_covs[fac_ind], as.character)
        
        plot_txt <- paste(colnames(other_covs), "=", other_covs, 
                          collapse = ", ")
      }
      
      if (name == "tpm") {
        p <- ggplot(df, aes(var, prob)) + 
          geom_ribbon(aes(ymin = low, ymax = upp), alpha = 0.3) +
          geom_line() + 
          facet_wrap(c("from", "to"), 
                     strip.position = "left",
                     labeller = label_bquote("Pr("*.(from)*" -> "*.(to)*")")) +
          xlab(var) + ylab(NULL) + ggtitle(plot_txt) +
          theme_light() +
          theme(strip.background = element_blank(),
                strip.placement = "outside", 
                strip.text = element_text(colour = "black")) + 
          coord_cartesian(ylim = c(0, 1))
        if(is.factor(df$var)) {
          p <- p + geom_point(size = 0.7) +
            geom_segment(aes(x = var, y = low, xend = var, yend = upp), alpha = 0.3)
        } else {
          p <- p + geom_line(size = 0.7) +
            geom_ribbon(aes(ymin = low, ymax = upp), alpha = 0.3)
        }
      } else {
        if (name == "delta") {
          p <- ggplot(df, aes(var, prob, group = state, col = state)) +
            scale_color_manual("", values = hmmTMB_cols) +
            scale_fill_manual(values = hmmTMB_cols, guide = FALSE) +
            xlab(var) + ylab("State probabilities") + ggtitle(plot_txt) +
            theme_light() + 
            coord_cartesian(ylim = c(0, 1))
        } else if (name == "obspar") {
          p <- ggplot(df, aes(var, val, col = state)) + theme_light() +
            scale_color_manual("", values = hmmTMB_cols) +
            scale_fill_manual(values = hmmTMB_cols, guide = FALSE) +
            facet_wrap(c("par"), scales = "free_y",
                       strip.position = "left",
                       labeller = label_bquote(.(par))) +
            xlab(var) + ylab(NULL) + ggtitle(plot_txt) +
            theme(strip.background = element_blank(),
                  strip.placement = "outside", 
                  strip.text = element_text(colour = "black"))
        }
        if(is.factor(df$var)) {
          p <- p + geom_point(size = 0.7) +
            geom_segment(aes(x = var, y = low, xend = var, yend = upp), alpha = 0.3)
        } else {
          p <- p + geom_line(size = 0.7) +
            geom_ribbon(aes(ymin = low, ymax = upp, fill = state), col = NA, alpha = 0.3)
        }
      }
      return(p)
    }, 
    
    

    # Print methods -----------------------------------------------------------
    #' @description Print model formulation
    formulation = function() {
      self$obs()$formulation()
      self$hidden()$formulation()
    }, 
    
    #' @description Print HMM object
    print = function() {
      self$formulation()
    }
  ),
  
  private = list(

    # Private data members ----------------------------------------------------
    obs_ = NULL,
    hidden_ = NULL,
    out_ = NULL,
    tmb_obj_ = NULL,
    tmb_rep_ = NULL,
    priors_ = NULL, 
    mcmc_ = NULL, 
    iters_= NULL,
    par_iters_ = NULL, 
    fixpar_ = NULL, 
    states_ = NULL,
    
    # Reading from spec file --------------------------------------------------
    
    ## @description Read model specification files
    ## 
    ## @param file file location 
    ## 
    ## @return List with elements:
    ## \itemize{
    ##   \item{\code{data}}{Data frame}
    ##   \item{\code{nstates}}{Number of states}
    ##   \item{\code{dists}}{List of observation distributions}
    ##   \item{\code{forms}}{Formulas for observation model}
    ##   \item{\code{tpm}}{Structure of hidden state model}
    ##   \item{\code{par}}{Initial parameters for observation model}
    ##   \item{\code{tpm0}}{Initial transition probability matrix}
    ##   \item{\code{fixed}}{Fixed parameters}
    ## }
    read_file = function(file) {
      if (!file.exists(file)) stop("model specification file does not exist in this location:", file)
      spec <- scan(file = file, 
                   character(), 
                   sep = "\n", 
                   strip.white = TRUE, 
                   quiet = TRUE)
      # remove comments
      spec <- strip_comments(spec)
      # block names 
      blocknms <- c("DATA", "DISTRIBUTION", "FORMULA", "TPM", "INITIAL", "FIXED")
      # find blocks
      wh_blocks <- sapply(blocknms, FUN = function(b) {grep(b, spec)})
      wh_blocks <- unlist(wh_blocks[sapply(wh_blocks, length) > 0])
      wh_blocks <- sort(wh_blocks)
      read_nms <- names(wh_blocks)
      
      # DATA
      if ("DATA" %in% read_nms) {
        data_block <- private$read_block("DATA", wh_blocks, spec)
        read_data_block <- private$read_equals(data_block)
        dataset <- get(read_data_block$rhs[which(read_data_block$lhs == "dataset")])
        nstates <- as.numeric(read_data_block$rhs[which(read_data_block$lhs == "nstates")])
      } else {
        stop("DATA block is missing from model specification")
      }
        
      # DISTS 
      if ("DISTRIBUTION" %in% read_nms) {
        dist_block <- private$read_block("DISTRIBUTION", wh_blocks, spec)
        dists <- private$read_dists(dist_block)
      } else {
        dists <- NULL
      }
      
      # FORMULA
      if ("FORMULA" %in% read_nms) {
        form_block <- private$read_block("FORMULA", wh_blocks, spec)
        forms <- private$read_forms(form_block)
      } else {
        forms <- NULL
      }
      
      # TPM
      if ("TPM" %in% read_nms) {
        tpm_block <- private$read_block("TPM", wh_blocks, spec)
        tpm <- matrix(str_trim(unlist(str_split(tpm_block, ":"))), 
                      nr = nstates, 
                      nc = nstates, 
                      byrow = TRUE)
      } else {
        tpm <- NULL
      }
        
      # INITIAL
      if ("INITIAL" %in% read_nms) {
        ini_block <- private$read_block("INITIAL", wh_blocks, spec)
        ini <- private$read_forms(ini_block, ini = TRUE, nstates = nstates)
      } else {
        stop("INITIAL block is missing from model specification")
      }
      
      # FIXED
      if ("FIXED" %in% read_nms) {
        fixed_block <- private$read_block("FIXED", wh_blocks, spec)
        fixed <- list(obs = rep(NA, length(fixed_block)))
        names(fixed$obs) <- fixed_block
      } else {
        fixed <- NULL
      }
      
      return(list(data = dataset, 
                  nstates = nstates, 
                  dists = dists, 
                  forms = forms, 
                  tpm = tpm, 
                  par = ini$par, 
                  tpm0 = ini$tpm0, 
                  fixed = fixed))
      
    }, 
    
    ## @description Read a specified block in a model specification file 
    read_block = function(name, wh_blocks, spec) {
      find <- which(names(wh_blocks) == name)
      start_block <- wh_blocks[find] + 1 
      end_block <- ifelse(start_block > max(wh_blocks), length(spec), wh_blocks[find + 1] - 1)
      block <- spec[start_block:end_block]
      return(block)
    }, 
    
    ## @description Separate left hand side and right hand side variables from a equation 
    ## @param x character vector of equations 
    ## @return list of left hand sides (lhs) and right hand sides (rhs)
    read_equals = function(x) {
      rhs <- str_trim(gsub(".*=", "", x))
      lhs <- str_trim(gsub("=.*", "", x))
      return(list(lhs = lhs, rhs = rhs))
    }, 
    
    ## @description Read distribution block 
    ## @param dist character vector of distribtion block 
    ## @return list of distributions 
    read_dists = function(dists) {
      ds <- strsplit(dists, "\n")
      terms <- sapply(ds, FUN = function(x) {all.vars(as.formula(x))})
      ls <- as.list(terms[2,])
      names(ls) <- terms[1,]
      return(ls)
    }, 
    
    ## @description Read both formula and initial blocks 
    ## @param forms character vector of block to read
    ## @param ini if TRUE then read as if it is initial block otherwise assume it
    ## is the formula block 
    read_forms = function(forms, ini = FALSE, nstates = NULL) {
      # find variables 
      wh_vars <- grep(":", forms)
      vars <- str_trim(gsub(":", "", forms[wh_vars]))
      par <- NULL
      tpm0 <- NULL
      for (i in 1:length(vars)){
        # find sub-block of formula/initial values for that variable 
        if (i < length(vars)) {
          end <- wh_vars[i + 1] - 1
        } else {
          end <- length(forms)
        }
        subforms <- forms[(wh_vars[i] + 1):end]
        if (str_trim(vars[i]) == "tpm") {
          nstates <- length(subforms)
          tpm0 <- matrix(0, nr = nstates, nstates)
          for (s in 1:nstates) {
            tpm0[s,] <- as.numeric(strsplit(subforms[s], ",")[[1]])
          }
        } else {
          subpar <- NULL
          subparnms <- NULL
          for (j in 1:length(subforms)) {
            # get parameter name 
            if (ini) {
              par_name <- str_trim(gsub("\\=.*", "", subforms[j]))
            } else{
              par_name <- str_trim(gsub("\\~.*", "", subforms[j]))
            }
            # get rhs of formula or initial values 
            if (ini){
              par_ini <- as.numeric(str_split(sub(".*\\=", "", subforms[j]), ",")[[1]])
              if (length(par_ini) == 1) par_ini <- rep(par_ini, nstates)
            } else {
              par_ini <- as.formula(paste0("~", gsub(".*\\~", "", subforms[j])))
            }
            subparnms <- c(subparnms, par_name)
            subpar <- c(subpar, list(par_ini))
          }
          names(subpar) <- subparnms
          par <- c(par, list(subpar))
        }
      }
      names(par) <- vars[vars != "tpm"]
      if (ini) {
        res <- list(par = par, tpm0 = tpm0)
      } else {
        res <- par 
      }
      return(res)
    },
    
    # Other private methods ---------------------------------------------------
    
    ## Compute effective degrees of freedom for a GAM 
    comp_edf = function(X, S, lambda){
      
      if(length(lambda)==0){
        return(0)
      }
      # duplicate lambda enough times
      lambda <- rep(lambda, nrow(S))
      # calculate lambda*S
      Sbig <- S * lambda
      
      # calculate the hat matrix
      XtX <- t(X) %*% X
      Fi <- solve(XtX + Sbig)
      F <- Fi %*% XtX
      
      # return the trace
      return(sum(diag(F)))
    }, 
    
    ## Initialize submodel component 
    initialize_submodel = function(par, initpar) {
      # find which parts of initializing model occur in current model
      wh <- match(rownames(initpar), rownames(par))
      # find which parts of initializing model do not occur in current model
      wh2 <- 1:nrow(initpar)
      wh2 <- wh2[!is.na(wh)]
      wh <- wh[!is.na(wh)]
      # copy over communal part 
      if (length(wh) > 0) {
        par[wh, 1] <- initpar[wh2, 1]
      }
      return(par)
    },
    
    ## Check constructor arguments 
    # (For argument description, see constructor)
    check_args = function(obs, hidden, init) {
      if(!inherits(obs, "Observation")) {
        stop("'obs' should be an Observation object")
      }
      if(!inherits(hidden, "MarkovChain")) {
        stop("'hidden' should be an MarkovChain object")
      }
      if(obs$nstates() != hidden$nstates()) {
        stop(paste0("The observation model and hidden state model should have the ",
                    "same number of states. (obs$nstates = ", obs$nstates(), 
                    ", hidden$nstates = ", hidden$nstates(), ")"))
      }
      if(!is.null(init)) if (!("HMM" %in% class(init))) stop("init must be a HMM object.")
    }
  )
)
