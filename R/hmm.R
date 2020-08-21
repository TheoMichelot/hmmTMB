
#' R6 class for hidden Markov model
#'
#' Encapsulates the observation and hidden state models for a hidden
#' Markov model.
Hmm <- R6Class(
  classname = "Hmm",
  
  public = list(
    #################
    ## Constructor ##
    #################
    #' @description Create new Hmm object
    #' 
    #' @param obs Observation object
    #' @param hidden MarkovChain object
    #' 
    #' @return A new Hmm object
    initialize = function(obs, hidden) {
      private$obs_ <- obs
      private$hidden_ <- hidden
    },
    
    ###############
    ## Accessors ##
    ###############
    #' @description Observation object for this model
    obs = function() {return(private$obs_)},
    
    #' @description MarkovChain object for this model
    hidden = function() {return(private$hidden_)},
    
    #' @description Output of optimiser after model fitting
    out = function() {
      if (is.null(private$out_)) {
        stop("Fit model first")
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
        stop("Setup model first")
      }
      
      return(private$tmb_obj_)
    },
    
    #' @description Output of the TMB function \code{sdreport}, which includes 
    #' estimates and standard errors for all model parameters.
    tmb_rep = function() {
      if(is.null(private$tmb_rep_)) {
        stop("Fit model first")
      }
      
      return(private$tmb_rep_)
    },
    
    #' @description Vector of estimated states, after \code{viterbi} has
    #' been run
    states = function() {
      if(is.null(private$states_)) {
        stop("Run viterbi first")
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
    
    #' @description Variance components of smooth terms
    #' 
    #' @details This function transforms the smoothness parameter of
    #' each smooth term into a standard deviation, given by 
    #' SD = 1/sqrt(lambda). It is particularly helpful to get the
    #' standard deviations of independent normal random effects.
    vcomp = function() {
      return(list(obs = self$obs()$vcomp(),
                  hidden = self$hidden()$vcomp()))
    },
    
    #' @description Model parameters
    #' 
    #' @return A list with elements:
    #' \itemize{
    #'   \item{\code{obspar}}{Parameters of observation model}
    #'   \item{\code{tpm}}{Transition probability matrix of hidden state model}
    #' }
    par = function() {
      obspar <- self$obs()$par()
      tpm <- self$hidden()$tpm()
      return(list(obspar = obspar, tpm = tpm))
    },
    
    #' @description Objective function
    #' 
    #' @param par Vector of parameters for which the function should be
    #' evaluated (on the working scale).
    #' 
    #' @return Negative log-likelihood
    nllk = function(par) {
      self$tmb_obj()$fn(par)
    },
    
    #################################
    ## Model fitting and inference ##
    #################################
    #' @description TMB setup
    #' 
    #' This creates an attribute \code{tmb_obj}, which can be used to 
    #' evaluate the negative log-likelihood function.
    #' 
    #' @param silent Logical. If TRUE, all tracing outputs are hidden (default).
    setup = function(silent = TRUE) {
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
      tmb_par <- list(coeff_fe_obs = self$obs()$coeff_fe(),
                      coeff_re_obs = 0,
                      log_lambda_obs = 0,
                      coeff_fe_hid = self$hidden()$coeff_fe(),
                      coeff_re_hid = 0,
                      log_lambda_hid = 0,
                      log_delta = rep(0, n_states - 1))
      
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
                       map = map, 
                       silent = silent)
      
      # Negative log-likelihood function
      private$tmb_obj_ <- obj
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
    #' @param silent Logical. If TRUE, all tracing outputs are hidden (default).
    fit = function(silent = TRUE) {
      # Setup if necessary
      if(is.null(private$tmb_obj_)) {
        self$setup(silent = silent)
      }
      
      # Number of states
      n_states <- self$hidden()$nstates()
      
      # Fit model
      private$out_ <- do.call(optim, private$tmb_obj_)
      
      # Get estimates and precision matrix for all parameters
      private$tmb_rep_ <- sdreport(private$tmb_obj_, getJointPrecision = TRUE, 
                                   skip.delta.method = FALSE)
      par_list <- as.list(private$tmb_rep_, "Estimate")
      
      # Observation parameters
      self$obs()$update_coeff_fe(coeff_fe = par_list$coeff_fe_obs)
      mats_obs <- self$obs()$make_mat()
      if(!is.null(mats_obs$ncol_re)) { # Only update if there are random effects
        self$obs()$update_coeff_re(coeff = par_list$coeff_re_obs)
        self$obs()$update_lambda(exp(par_list$log_lambda_obs))
      }
      
      # Transition probabilities
      self$hidden()$update_coeff_fe(coeff_fe = par_list$coeff_fe_hid)
      mats_hid <- self$hidden()$make_mat(data = self$obs()$data()$data())
      if(!is.null(mats_hid$ncol_re)) { # Only update if there are random effects
        self$hidden()$update_coeff_re(coeff_re = par_list$coeff_re_hid)
        self$hidden()$update_lambda(exp(par_list$log_lambda_hid))
      }
    },
    
    #' @description Viterbi algorithm
    #' 
    #' @return Most likely state sequence
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
                                       X_re = mod_mat_hid$X_re)
      
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
    },
    
    #' @description Confidence intervals for model parameters on the working scale 
    #' 
    #' These are Wald confidence intervals, obtained from the standard errors 
    #' returned by the TMB function \code{sdreport}. See the TMB documentation 
    #' for more details.
    #' 
    #' @param level Confidence level (default: 0.95 for 95\% confidence intervals)
    #' 
    #' @return Matrix with three columns: (1) estimates, (2) lower bounds of
    #' confidence intervals, (3) upper bounds of confidence intervals.
    CI_coeff = function(level = 0.95) {
      # Extract parameter estimates and standard errors from TMB output
      par_list <- as.list(self$tmb_rep(), "Estimate")
      se_list <- as.list(self$tmb_rep(), "Std. Error")
      
      # Lower bounds
      lower <- lapply(seq_along(par_list), function(i) {
        par_list[[i]] - qnorm(1 - (1 - level)/2) * se_list[[i]]
      })
      # Upper bounds
      upper <- lapply(seq_along(par_list), function(i) {
        par_list[[i]] + qnorm(1 - (1 - level)/2) * se_list[[i]]
      })
      
      return(cbind(estimate = unlist(par_list),
                   lower = unlist(lower),
                   upper = unlist(upper)))
    },
    
    #' @description Confidence intervals for transition probabilities
    #'
    #' @param X_fe Design matrix for fixed effects, as returned
    #' by \code{make_mat}
    #' @param X_re Design matrix for random effects, as returned
    #' by \code{make_mat}
    #' @param level Confidence level (default: 0.95 for 95\% confidence 
    #' intervals)
    #' @param n_post Number of posterior samples from which the confidence
    #' intervals are calculated. Larger values will reduce approximation
    #' error, but increase computation time. Defaults to 1000.
    #' 
    #' @details This method generates confidence intervals by simulation.
    #' That is, it generates \code{n_post} posterior samples of 
    #' the estimated parameters from a multivariate normal distribution,
    #' where the mean is the vector of estimates and the covariance matrix 
    #' is provided by TMB. Then, transition probabilities are derived for 
    #' each set of posterior parameter values, and confidence intervals
    #' are obtained as quantiles of the posterior simulated transition
    #' probabilities.
    #' 
    #' @return List with elements:
    #' \itemize{
    #'   \item{\code{low}}{Array of lower bounds of confidence intervals.
    #'   Each slice of the array is a transiton probability matrix.}
    #'   \item{\code{upp}}{Array of upper bounds of confidence intervals.
    #'   Each slice of the array is a transiton probability matrix.}
    #' }
    CI_tpm = function(X_fe, X_re, level = 0.95, n_post = 1e3) {
      # Number of states
      n_states <- self$hidden()$nstates()
      # Number of time steps
      n_grid <- nrow(X_fe)/(n_states*(n_states-1))
      
      # Get parameter estimates and covariance matrix
      rep <- self$tmb_rep()
      if(is.null(rep$jointPrecision)) {
        par <- rep$par.fixed
        V <- rep$cov.fixed    
      } else {
        par <- c(rep$par.fixed, rep$par.random)
        V <- ginv(rep$jointPrecision)
      }
      
      # Generate samples from MVN estimator distribution
      post <- rmvn(n = n_post, mu = par, V = V)
      
      # Extract coefficients for transition probabilities
      post_coeff_fe <- post[, which(colnames(post) == "coeff_fe_hid")]
      post_coeff_re <- post[, which(colnames(post) == "coeff_re_hid")]
      
      # Get transition probabilities over rows of X_fe and X_re, for each 
      # posterior sample of coeff_fe and coeff_re
      post_tpm <- sapply(1:n_post, function(i) {
        hid$tpm_all(X_fe = X_fe, X_re = X_re, 
                    coeff_fe = post_coeff_fe[i,], 
                    coeff_re = post_coeff_re[i,])
      })
      
      # Get confidence intervals as quantiles of posterior tpms
      alpha <- (1 - level)/2
      CI <- t(apply(post_tpm, 1, quantile, probs = c(alpha, 1 - alpha)))
      
      # Format arrays for lower and upper bounds, where each layer is
      # a transition probability matrix
      low <- array(CI[,1], dim = c(n_states, n_states, n_grid))
      upp <- array(CI[,2], dim = c(n_states, n_states, n_grid))
      
      return(list(low = low, upp = upp))
    },
    
    #' @description Confidence intervals for stationary distribution
    #'
    #' @param X_fe Design matrix for fixed effects, as returned
    #' by \code{make_mat}
    #' @param X_re Design matrix for random effects, as returned
    #' by \code{make_mat}
    #' @param level Confidence level (default: 0.95 for 95\% confidence 
    #' intervals)
    #' @param n_post Number of posterior samples from which the confidence
    #' intervals are calculated. Larger values will reduce approximation
    #' error, but increase computation time. Defaults to 1000.
    #' 
    #' @details This method generates confidence intervals by simulation.
    #' That is, it generates \code{n_post} posterior samples of 
    #' the estimated parameters from a multivariate normal distribution,
    #' where the mean is the vector of estimates and the covariance matrix 
    #' is provided by TMB. Then, stationary state probabilities are derived 
    #' for each set of posterior parameter values, and confidence intervals
    #' are obtained as quantiles of the posterior simulated stationary state
    #' probabilities.
    #' 
    #' @return List with elements:
    #' \itemize{
    #'   \item{\code{low}}{Matrix of lower bounds of confidence intervals,
    #'   with one row for each time step, and one column for each state.}
    #'   \item{\code{upp}}{Matrix of upper bounds of confidence intervals,
    #'   with one row for each time step, and one column for each state.}
    #' }
    CI_statdist = function(X_fe, X_re, level = 0.95, n_post = 1e3) {
      # Number of states
      n_states <- self$hidden()$nstates()
      # Number of time steps
      n_grid <- nrow(X_fe)/(n_states*(n_states-1))
      
      # Get parameter estimates and covariance matrix
      rep <- self$tmb_rep()
      if(is.null(rep$jointPrecision)) {
        par <- rep$par.fixed
        V <- rep$cov.fixed    
      } else {
        par <- c(rep$par.fixed, rep$par.random)
        V <- ginv(rep$jointPrecision)
      }
      
      # Generate samples from MVN estimator distribution
      post <- rmvn(n = n_post, mu = par, V = V)
      
      # Extract coefficients for transition probabilities
      post_coeff_fe <- post[, which(colnames(post) == "coeff_fe_hid")]
      post_coeff_re <- post[, which(colnames(post) == "coeff_re_hid")]
      
      # Get stationary distributions over rows of X_fe and X_re, for each 
      # posterior sample of coeff_fe and coeff_re
      post_statdist <- sapply(1:n_post, function(i) {
        tpms <- hid$tpm_all(X_fe = X_fe, X_re = X_re, 
                            coeff_fe = post_coeff_fe[i,], 
                            coeff_re = post_coeff_re[i,])
        
        # For each transition matrix, get corresponding stationary distribution
        stat_dists <- apply(tpms, 3, function(tpm)
          solve(t(diag(n_states) - tpm + 1), rep(1, n_states)))
        stat_dists <- t(stat_dists)
        return(stat_dists)
      })
      
      # Get confidence intervals as quantiles of posterior tpms
      alpha <- (1 - level)/2
      CI <- t(apply(post_statdist, 1, quantile, probs = c(alpha, 1 - alpha)))
      
      # Format arrays for lower and upper bounds, where each layer is
      # a transition probability matrix
      low <- matrix(CI[,1], nrow = n_grid, ncol = n_states)
      upp <- matrix(CI[,2], nrow = n_grid, ncol = n_states)
      
      return(list(low = low, upp = upp))
    },
    
    #' @description Confidence intervals for observation parameters
    #'
    #' @param X_fe Design matrix for fixed effects, as returned
    #' by \code{make_mat}
    #' @param X_re Design matrix for random effects, as returned
    #' by \code{make_mat}
    #' @param level Confidence level (default: 0.95 for 95\% confidence 
    #' intervals)
    #' @param n_post Number of posterior samples from which the confidence
    #' intervals are calculated. Larger values will reduce approximation
    #' error, but increase computation time. Defaults to 1000.
    #' 
    #' @details This method generates confidence intervals by simulation.
    #' That is, it generates \code{n_post} posterior samples of 
    #' the estimated parameters from a multivariate normal distribution,
    #' where the mean is the vector of estimates and the covariance matrix 
    #' is provided by TMB. Then, observation parameters are derived for 
    #' each set of posterior parameter values, and confidence intervals
    #' are obtained as quantiles of the posterior simulated observation
    #' parameters.
    #' 
    #' @return List with elements:
    #' \itemize{
    #'   \item{\code{low}}{Array of lower bounds of confidence intervals,
    #'   with one slice for each time step, one row for each observation 
    #'   parameter, and one column for each state.}
    #'   \item{\code{upp}}{Array of upper bounds of confidence intervals,
    #'   with one slice for each time step, one row for each observation 
    #'   parameter, and one column for each state.}
    #' }
    CI_obspar = function(X_fe, X_re, level = 0.95, n_post = 1e3) {
      # Number of states
      n_states <- self$hidden()$nstates()
      # Number of observation parameters (in each state)
      n_par <- sum(sapply(self$obs()$dists(), function(d) d$npar()))
      # Number of time steps
      n_grid <- nrow(X_fe)/(n_par * n_states)
      
      # Get parameter estimates and covariance matrix
      rep <- self$tmb_rep()
      if(is.null(rep$jointPrecision)) {
        par <- rep$par.fixed
        V <- rep$cov.fixed    
      } else {
        par <- c(rep$par.fixed, rep$par.random)
        V <- ginv(rep$jointPrecision)
      }
      
      # Generate samples from MVN estimator distribution
      post <- rmvn(n = n_post, mu = par, V = V)
      
      # Extract coefficients for observation parameters
      post_coeff_fe <- post[, which(colnames(post) == "coeff_fe_obs")]
      post_coeff_re <- post[, which(colnames(post) == "coeff_re_obs")]
      
      # Get observation parameters over rows of X_fe and X_re, for each 
      # posterior sample of coeff_fe and coeff_re
      post_obspar <- sapply(1:n_post, function(i) {
        obs$par_all(X_fe = X_fe, X_re = X_re, 
                    coeff_fe = post_coeff_fe[i,], 
                    coeff_re = post_coeff_re[i,])
      })
      
      # Get confidence intervals as quantiles of posterior parameters
      alpha <- (1 - level)/2
      CI <- t(apply(post_obspar, 1, quantile, probs = c(alpha, 1 - alpha)))
      
      # Format arrays for lower and upper bounds, with one slice for each 
      # time step, one row for each observation parameter, and one column 
      # for each state.
      low <- array(CI[,1], dim = c(n_par, n_states, n_grid))
      upp <- array(CI[,2], dim = c(n_par, n_states, n_grid))
      
      # Hacky way to get parameter names
      par_names <- names(unlist(rapply(
        self$obs()$w2n(post_obspar[1,]), function(v) v[1])))
      
      # Set dimension names for rows and columns
      dimnames(low) <- list(par_names, paste("state", 1:n_states), NULL)
      dimnames(upp) <- list(par_names, paste("state", 1:n_states), NULL)
      
      return(list(low = low, upp = upp))
    },
    
    #' @description Simulate from hidden Markov model
    #' 
    #' @param n Number of time steps to simulate
    #' @param data Optional data frame including covariates
    #' 
    #' @return Data frame including columns of data (if provided), and simulated
    #' data variables
    simulate = function(n, data = NULL) {
      if(is.null(data)) {
        data <- data.frame(ID = rep(factor(1), n))
      } else if(is.null(data$ID)) {
        data$ID <- rep(factor(1), n)
      }
      
      # Number of states
      n_states <- self$hidden()$nstates()
      
      # Simulate state process      
      S <- self$hidden()$simulate(n = n, data = self$obs()$data()$data(), 
                                  new_data = data)
      
      # Create observation parameters
      mats_obs <- self$obs()$make_mat(new_data = data)
      par_obs <- self$obs()$par_all(X_fe = mats_obs$X_fe, X_re = mats_obs$X_re,
                                    full_names = FALSE)
      
      # Simulate observation process
      obs_dists <- self$obs()$dists()
      obs_all <- data.frame(state = S)
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
            cat("\rSimulating ", var_name, "... ", round(i/n*100), "%", sep = "")        
          }
          
          # Generate realisation
          obs[i] <- obsdist$rng_apply(n = 1, par = par_obs[par_ind, S[i], i])
        }
        
        # Add variable to data frame
        obs_all[[var_name]] <- obs
        par_count <- par_count + obsdist$npar()
        cat("\n")
      }
      
      # Combine original data set and simulated variables
      obs_all <- cbind(data, obs_all)
      
      return(obs_all)
    },
    
    ######################
    ## Plotting methods ##
    ######################
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
      data <- self$obs()$data()$data()
      # State sequence as factor
      state <- as.factor(self$states())
      # Time series ID
      ID <- self$obs()$data()$ID()
      
      if(is.null(var2)) {
        # 1d time series plot
        df <- data.frame(index = 1:nrow(data),
                         x = data[[var1]])
        
        p <- ggplot(data = df, mapping = aes(index, x, col = state, group = ID)) +
          geom_line() + xlab("time") + ylab(var1)
      } else {
        # 2d plot
        df <- data.frame(x = data[[var1]],
                         y = data[[var2]])
        
        p <- ggplot(data = df, mapping = aes(x, y, col = state, group = ID)) +
          geom_path() + xlab(var1) + ylab(var2)
      }
      
      p <- p + 
        scale_color_manual(values = hmmTMB_cols) +
        theme_light()
      
      return(p)
    },
    
    #' @description Plot transition probability matrix
    #' 
    #' @param var Name of covariate as a function of which the transition
    #' probabilities should be plotted
    #' @param covs Optional data frame with a single row and one column
    #' for each covariate, giving the values that should be used. If this is
    #' not specified, the mean value is used for numeric variables, and the
    #' first level for factor variables.
    #' 
    #' @return A ggplot object
    plot_tpm = function(var, covs = NULL) {
      # Number of states
      n_states <- self$hidden()$nstates()
      
      # Create design matrices for grid of "var" values
      mats <- self$hidden()$make_mat_grid(
        var = var, data = self$obs()$data()$data(), 
        covs = covs, n_grid = 50)
      # Get transition probability matrices
      tpms <- self$hidden()$tpm_all(X_fe = mats$X_fe, X_re = mats$X_re)
      # Get confidence intervals
      CIs <- self$CI_tpm(X_fe = mats$X_fe, X_re = mats$X_re)
      
      # Data frame for plot
      df <- as.data.frame.table(tpms)
      colnames(df) <- c("from", "to", "var", "prob")
      levels(df$from) <- paste("State", 1:n_states)
      levels(df$to) <- paste("State", 1:n_states)
      df$var <- rep(mats$new_data[, var], each = n_states * n_states)
      df$low <- as.vector(CIs$low)
      df$upp <- as.vector(CIs$upp)
      
      # Create caption with values of other (fixed) covariates      
      plot_txt <- NULL
      if(ncol(mats$new_data) > 1) {
        other_covs <- mats$new_data[1, which(colnames(mats$new_data) != var),
                                    drop = FALSE]
        plot_txt <- paste(colnames(other_covs), "=", round(other_covs, 2), 
                          collapse = ", ")
      }
      
      # Create plot using facets
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
      
      return(p)
    },
    
    #' @description Plot stationary state probabilities
    #' 
    #' @param var Name of covariate as a function of which the state
    #' probabilities should be plotted
    #' @param covs Optional data frame with a single row and one column
    #' for each covariate, giving the values that should be used. If this is
    #' not specified, the mean value is used for numeric variables, and the
    #' first level for factor variables.
    #' 
    #' @return A ggplot object
    plot_stat_dist = function(var, covs = NULL) {
      # Number of states
      n_states <- self$hidden()$nstates()
      
      # Create design matrices for grid of "var" values
      mats <- self$hidden()$make_mat_grid(
        var = var, data = self$obs()$data()$data(), 
        covs = covs, n_grid = 50)
      # Get stationary distributions
      stat_dists <- self$hidden()$stat_dists(X_fe = mats$X_fe, X_re = mats$X_re)
      # Get confidence intervals
      CIs <- self$CI_statdist(X_fe = mats$X_fe, X_re = mats$X_re)
      
      # Data frame for plot
      df <- as.data.frame.table(stat_dists)
      colnames(df) <- c("var", "state", "prob")
      levels(df$state) <- paste("State", 1:n_states)
      df$var <- rep(mats$new_data[, var], n_states)
      df$low <- as.vector(CIs$low)
      df$upp <- as.vector(CIs$upp)
      
      # Create caption with values of other (fixed) covariates      
      plot_txt <- NULL
      if(ncol(mats$new_data) > 1) {
        other_covs <- mats$new_data[1, which(colnames(mats$new_data) != var),
                                    drop = FALSE]
        plot_txt <- paste(colnames(other_covs), "=", round(other_covs, 2), 
                          collapse = ", ")
      }
      
      # Create plot
      p <- ggplot(df, aes(var, prob, group = state, col = state)) +
        geom_ribbon(aes(ymin = low, ymax = upp), alpha = 0.3) +
        geom_line(size = 0.7) + scale_color_manual("", values = hmmTMB_cols) +
        xlab(var) + ylab("State probabilities") + ggtitle(plot_txt) +
        theme_light() + 
        coord_cartesian(ylim = c(0, 1))
      
      return(p)
    },
    
    #' @description Plot observation parameters
    #' 
    #' @param var Name of covariate as a function of which the parameters
    #' should be plotted
    #' @param covs Optional data frame with a single row and one column
    #' for each covariate, giving the values that should be used. If this is
    #' not specified, the mean value is used for numeric variables, and the
    #' first level for factor variables.
    #' 
    #' @return A ggplot object
    plot_obspar = function(var, covs = NULL) {
      # Number of states
      n_states <- self$hidden()$nstates()
      
      # Create design matrices for grid of "var" values
      mats <- self$obs()$make_mat_grid(var = var, covs = covs, n_grid = 50)
      # Get observation parameters over covariate grid
      obs_par <- self$obs()$par_all(X_fe = mats$X_fe, X_re = mats$X_re)
      # Get confidence intervals over covariate grid
      CIs <- self$CI_obspar(X_fe = mats$X_fe, X_re = mats$X_re)
      
      # Data frame for plot
      df <- as.data.frame.table(obs_par)
      colnames(df) <- c("par", "state", "var", "val")
      levels(df$state) <- paste("State", 1:n_states)
      df$var <- rep(mats$new_data[, var], each = nrow(df)/nrow(mats$new_data))
      df$low <- as.vector(CIs$low)
      df$upp <- as.vector(CIs$upp)
      
      # Create caption with values of other (fixed) covariates      
      plot_txt <- NULL
      if(ncol(mats$new_data) > 1) {
        other_covs <- mats$new_data[1, which(colnames(mats$new_data) != var), 
                                    drop = FALSE]
        plot_txt <- paste(colnames(other_covs), "=", round(other_covs, 2), 
                          collapse = ", ")
      }
      
      # Create plot
      p <- ggplot(df, aes(var, val, col = state)) + theme_light() +
        geom_ribbon(aes(ymin = low, ymax = upp), alpha = 0.3) +
        geom_line(size = 0.7) + scale_color_manual("", values = hmmTMB_cols) +
        facet_wrap(c("par"), scales = "free_y",
                   strip.position = "left",
                   labeller = label_bquote(.(as.character(par)))) +
        xlab(var) + ylab(NULL) + ggtitle(plot_txt) +
        theme(strip.background = element_blank(),
              strip.placement = "outside", 
              strip.text = element_text(colour = "black"))
      
      return(p)
    },
    
    ###################
    ## Other methods ##
    ###################
    #' @description Print model formulation
    formulation = function() {
      self$obs()$formulation()
      self$hidden()$formulation()
    }
  ),
  
  private = list(
    obs_ = NULL,
    hidden_ = NULL,
    out_ = NULL,
    tmb_obj_ = NULL,
    tmb_rep_ = NULL,
    states_ = NULL
  )
)
