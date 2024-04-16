
#' @title R6 class for hidden Markov model
#'
#' @description
#' Encapsulates the observation and hidden state models for a hidden
#' Markov model.
#' 
#' @importFrom R6 R6Class
#' @importFrom mgcv gam rmvn dmvn
#' @importFrom ggplot2 ggplot aes theme_light geom_line theme scale_colour_manual
#' facet_wrap label_bquote xlab ylab ggtitle element_blank element_text geom_point
#' geom_ribbon scale_size_manual geom_histogram geom_vline geom_errorbar after_stat
#' @importFrom TMB MakeADFun sdreport
#' @importFrom stringr str_trim str_split str_split_fixed
#' @importFrom optimx optimx
#' @importFrom tmbstan tmbstan
#' 
#' @useDynLib hmmTMB, .registration = TRUE
#' 
#' @export
HMM <- R6Class(
  classname = "HMM",
  
  public = list(
    
    # Constructor -------------------------------------------------------------
    #' @description Create new HMM object
    #' 
    #' @param obs Observation object, created with \code{Observation$new()}.
    #' This contains the formulation for the observation model.
    #' @param hid MarkovChain object, created with \code{MarkovChain$new()}.
    #' This contains the formulation for the state process model.
    #' @param file Path to specification file for HMM. If this argument is
    #' used, then \code{obs} and \code{hid} are unnecessary.
    #' @param init HMM object, used to initialise the parameters for this model.
    #' If \code{init} is passed, then all parameters that are included in init
    #' and in the present model are copied. This may be useful when fitting
    #' increasingly complex models: start from a simple model, then pass it as
    #' init to create a more complex model, and so on. 
    #' @param fixpar Named list, with optional elements: 'hid', 'obs', 'delta0',
    #' 'lambda_obs', and 'lambda_hid'. Each element is a named vector of 
    #' parameters in coeff_fe that should either be fixed (if the corresponding
    #' element is set to NA) or estimated to a common value (using integers or
    #' factor levels).don See examples in the vignettes, and check the TMB
    #' documentation to understand the inner workings (argument \code{map}
    #' of \code{TMB::MakeADFun()}).
    #' 
    #' @return A new HMM object
    #' 
    #' @examples
    #' # Load data set (included with R)
    #' data(nottem)
    #' data <- data.frame(temp = as.vector(t(nottem)))
    #' 
    #' # Create hidden state and observation models
    #' hid <- MarkovChain$new(data = data, n_states = 2)
    #' par0 <- list(temp = list(mean = c(40, 60), sd = c(5, 5)))
    #' obs <- Observation$new(data = data, n_states = 2, 
    #'                        dists = list(temp = "norm"),
    #'                        par = par0)
    #' 
    #' # Create HMM
    #' hmm <- HMM$new(hid = hid, obs = obs)
    initialize = function(obs = NULL, 
                          hid = NULL,
                          file = NULL,
                          init = NULL, 
                          fixpar = NULL) {
      # Decide how model has been specified 
      if (is.null(file) & is.null(obs)) {
        stop(paste0("Either 'file' should be the name of a file specifying ",
                    "the model or both obs/hid model objects should be ",
                    "supplied."))
      }
      if (!is.null(file) & is.null(obs)) {
        spec <- private$read_file(file)
        # create obs
        obs <- Observation$new(data = spec$data, 
                               dists = spec$dists,
                               n_states = spec$nstates, 
                               par = spec$par, 
                               formulas = spec$forms)
        hid <- MarkovChain$new(n_states = spec$nstates, 
                               formula = spec$tpm, 
                               data = spec$data)
        if(!is.null(spec$delta0)) hid$update_delta0(spec$delta0)
        if (!is.null(spec$fixed)) fixpar <- spec$fixed 
        if (!is.null(spec$tpm0))  hid$update_tpm(spec$tpm0)
      }
      
      # Check arguments
      private$check_args(obs = obs, hid = hid, init = init)
      
      # Get names of all covariates
      var_names <- unique(c(rapply(hid$formulas(), all.vars), 
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
      
      # Store sub-model components 
      private$obs_ <- obs
      private$hid_ <- hid
      
      # Update/overwrite fixpar of each model component if necessary
      if(!is.null(fixpar)) {
        self$hid()$update_fixpar(fixpar = fixpar[
          intersect(names(fixpar), c("hid", "lambda_hid", "delta0"))])
        self$obs()$update_fixpar(fixpar = fixpar[
          intersect(names(fixpar), c("obs", "lambda_obs"))])
      }
      
      if (!is.null(init)) {
        # Copy parameters with matching names from 'init' model
        private$obs_$update_coeff_fe(
          private$initialize_submodel(private$obs_$coeff_fe(), init$obs()$coeff_fe()))
        private$obs_$update_coeff_re(
          private$initialize_submodel(private$obs_$coeff_re(), init$obs()$coeff_re()))
        private$hid_$update_coeff_fe(
          private$initialize_submodel(private$hid_$coeff_fe(), init$hid()$coeff_fe()))
        private$hid_$update_coeff_re(
          private$initialize_submodel(private$hid_$coeff_re(), init$hid()$coeff_re()))
      }
      
      # initialize priors 
      self$set_priors()
    },
    
    
    # Accessors ---------------------------------------------------------------
    
    #' @description Observation object for this model
    obs = function() {return(private$obs_)},
    
    #' @description MarkovChain object for this model
    hid = function() {return(private$hid_)},
    
    #' @description Output of optimiser after model fitting
    out = function() {
      if (is.null(private$out_)) {
        stop("Fit model first using $fit()")
      }
      return(private$out_)
    },
    
    #' @description Model object created by TMB. This is the output of 
    #' the TMB function \code{MakeADFun}, and it is a list including elements
    #' \describe{
    #'   \item{\code{fn}}{Objective function}
    #'   \item{\code{gr}}{Gradient function of fn}
    #'   \item{\code{par}}{Vector of initial parameters on working scale}
    #' }
    tmb_obj = function() {
      if(is.null(private$tmb_obj_)) {
        stop("Setup or fit model first")
      }
      return(private$tmb_obj_)
    },
    
    #' @description Model object created by TMB for the joint likelihood of
    #' the fixed and random effects. This is the output of the TMB function 
    #' \code{MakeADFun}, and it is a list including elements
    #' \describe{
    #'   \item{\code{fn}}{Objective function}
    #'   \item{\code{gr}}{Gradient function of fn}
    #'   \item{\code{par}}{Vector of initial parameters on working scale}
    #' }
    tmb_obj_joint = function() {
      if(is.null(private$tmb_obj_joint_)) {
        stop("Setup model first")
      }
      
      return(private$tmb_obj_joint_)
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
        return(self$viterbi())
      } else {
        return(private$states_)        
      }
    },
    
    #' @description Coefficients for fixed effect parameters
    coeff_fe = function() {
      return(list(obs = self$obs()$coeff_fe(),
                  hid = self$hid()$coeff_fe()))
    },
    
    #' @description Coefficients for random effect parameters
    coeff_re = function() {
      return(list(obs = self$obs()$coeff_re(),
                  hid = self$hid()$coeff_re()))
    },
    
    #' @description List of all model coefficients
    #' 
    #' These are the parameters estimated by the model, including
    #' fixed and random effect parameters for the observation parameters
    #' and the transition probabilities, (transformed) initial
    #' probabilities, and smoothness parameters.
    coeff_list = function() {
      list(coeff_fe_obs = self$obs()$coeff_fe(),
           log_lambda_obs = self$obs()$lambda(), 
           coeff_fe_hid = self$hid()$coeff_fe(), 
           log_lambda_hid = self$hid()$lambda(), 
           log_delta0 = self$hid()$delta0(log = TRUE, as_matrix = FALSE), 
           coeff_re_obs = self$obs()$coeff_re(), 
           coeff_re_hid = self$hid()$coeff_re())
    },
    
    #' @description Fixed parameters
    #' 
    #' @param all Logical. If FALSE, only user-specified fixed
    #' parameters are returned, but not parameters that are fixed
    #' by definition (e.g., size of binomial distribution).
    fixpar = function(all = FALSE) {
      c(self$hid()$fixpar(all = all), self$obs()$fixpar(all = all))
    }, 
    
    #' @description Array of working parameters
    coeff_array = function() {
      return(private$coeff_array_)
    },
    
    #' @description Smoothness parameters
    lambda = function() {
      return(list(obs = self$obs()$lambda(),
                  hid = self$hid()$lambda()))
    },
    
    #' @description Update parameters stored inside model object
    #' 
    #' @param par_list List with elements for coeff_fe_obs, 
    #' coeff_fe_hid, coeff_re_obs, coeff_re_hid, log_delta0, 
    #' log_lambda_hid, and log_lambda_obs
    #' @param iter Optional argument to update model parameters based on MCMC
    #' iterations (if using rstan). Either the index of the iteration to use,
    #' or "mean" if the posterior mean should be used.
    update_par = function(par_list = NULL, iter = NULL) {
      if (is.null(par_list) & is.null(iter)) {
        stop("No new parameter values to update to")
      }
      if (!is.null(iter) & !is.null(par_list)) {
        stop("Either specify iter or par_list, not both")
      }
      if (!is.null(iter)) {
        # update to MCMC iteration 
        if (is.null(private$iters_)) {
          stop("Must run fit_stan() before using iterations")
        }
        if (is.numeric(iter)) {
          if (iter > dim(private$iters_)[1]) {
            stop("'iter' exceeds number of MCMC iterations available")
          }
          samp <- private$iters_[iter,]
        } else if (iter == "mean") {
          samp <- colMeans(private$iters_)
        } else {
          stop("Invalid iter to update_par()")
        }
        par_list <- split(samp, names(samp))
      }
      # Update observation parameters
      self$obs()$update_coeff_fe(coeff_fe = par_list$coeff_fe_obs)
      if(!is.null(self$obs()$terms()$ncol_re)) { 
        # Only update if there are random effects
        self$obs()$update_coeff_re(coeff = par_list$coeff_re_obs)
        self$obs()$update_lambda(exp(par_list$log_lambda_obs))
      }
      
      # Update transition probabilities
      self$hid()$update_coeff_fe(coeff_fe = par_list$coeff_fe_hid)
      if(!is.null(self$hid()$terms()$ncol_re)) { 
        # Only update if there are random effects
        self$hid()$update_coeff_re(coeff_re = par_list$coeff_re_hid)
        self$hid()$update_lambda(exp(par_list$log_lambda_hid))
      }
      
      # Update initial distribution delta0
      ID <- self$obs()$data()$ID
      n_ID <- length(unique(ID))
      n_states <- self$hid()$nstates()
      if (self$hid()$stationary()) {
        # Find stationary distribution for each ID
        first_indices <- c(1, which(ID[-1] != ID[-length(ID)]) + 1)
        delta0 <- self$hid()$delta(t = first_indices)
      } else {
        # Fill delta0 except reference elements
        delta0 <- t(sapply(1:n_ID, function(i) {
          ld0 <- rep(0, length = n_states)
          ind <- ((n_states - 1)*(i - 1) + 1):((n_states - 1) * i)
          ld0[-self$hid()$ref_delta0()[i]] <- par_list$log_delta0[ind]
          d0 <- exp(ld0)/sum(exp(ld0))
          return(d0)
        }))
      }
      self$hid()$update_delta0(delta0)
      
      # Update coeff_array
      private$coeff_array_[,"value"] <- unlist(self$coeff_list(), use.names = FALSE)
    }, 
    
    #' @description Standard deviation of smooth terms (or random effects)
    #' 
    #' This function transforms the smoothness parameter of
    #' each smooth term into a standard deviation, given by 
    #' SD = 1/sqrt(lambda). It is particularly helpful to get the
    #' standard deviations of independent normal random effects.
    #' 
    #' @return List of standard deviations for observation model and
    #' hidden state model.
    sd_re = function() {
      if(is.null(private$tmb_rep_)) {
        stop("Fit model first")
      }
      return(list(obs = self$obs()$sd_re(),
                  hid = self$hid()$sd_re()))
    },
    
    #' @description Model parameters
    #' @param t returns parameters at time t, default is t = 1 
    #' 
    #' @return A list with elements:
    #' \describe{
    #'   \item{\code{obspar}}{Parameters of observation model}
    #'   \item{\code{tpm}}{Transition probability matrix of hidden state model}
    #' }
    par = function(t = 1) {
      obspar <- self$obs()$par(t = t)
      tpm <- self$hid()$tpm(t = t)
      return(list(obspar = obspar, tpm = tpm))
    },
    
    #' @description Set priors for coefficients 
    #' 
    #' @param new_priors is a named list of matrices with optional elements 
    #' coeff_fe_obs, coeff_fe_hid, log_lambda_obs, andlog_lambda_hid.
    #' Each matrix has two columns (first col = mean, second col = sd) 
    #' specifying parameters for normal priors. 
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
        coeff_fe_hid_prior <- matrix(NA, nr = length(fe$hid), nc = 2)
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
        log_lambda_hid_prior <- matrix(NA, nr = length(lam$hid), nc = 2)
      }
      
      # Name rows and columns for readability
      rownames(coeff_fe_obs_prior) <- rownames(fe$obs)
      rownames(coeff_fe_hid_prior) <- rownames(fe$hid)
      rownames(log_lambda_obs_prior) <- rownames(lam$obs)
      rownames(log_lambda_hid_prior) <- rownames(lam$hid)
      colnames(coeff_fe_obs_prior) <- c("mean", "sd")
      colnames(coeff_fe_hid_prior) <- c("mean", "sd")
      colnames(log_lambda_obs_prior) <- c("mean", "sd")
      colnames(log_lambda_hid_prior) <- c("mean", "sd")
      
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
        stop("Run fit_stan() before extracting iterations")
      }
      if (type == "response") {
        return(private$par_iters_)
      } else if (type == "raw") {
        return(private$iters_)
      } else {
        stop("Unknown type argument given to iters()")
      }
    }, 
    
    #' @description fitted stan object from MCMC fit 
    #' 
    #' @return the stanfit object 
    out_stan = function() {
      if (is.null(private$iters_)) {
        stop("Run fit_stan() first")
      }
      return(private$out_stan_)
    }, 
    
    #' @description Log-likelihood at current parameters
    #' 
    #' @return Log-likelihood
    llk = function() {
      # Setup if necessary
      if(is.null(private$tmb_obj_)) {
        self$setup(silent = TRUE)
      }
      # compute log-likelihood
      return(-self$tmb_obj()$fn(self$tmb_obj()$par))
    },
    
    #' @description Effective degrees of freedom
    #'
    #' @return Number of effective degrees of freedom
    #' (accounting for flexibility in non-parametric 
    #' terms implied by smoothing)
    edf = function() {
      # Degrees of freedom for fixed effects
      # (don't count smoothing parameters)
      edf <- nrow(self$obs()$coeff_fe()) + 
        nrow(self$hid()$coeff_fe())
      if(!self$hid()$stationary()) {
        edf <- edf + length(self$coeff_list()$log_delta0)
      }
      
      if(!is.null(private$tmb_rep_)) {
        if(!is.null(self$tmb_rep()$jointPrecision)) {
          # Joint covariance matrix
          Q <- self$tmb_rep()$jointPrecision
          V <- prec_to_cov(Q)
          
          # get Hessian
          par_all <- c(self$tmb_rep()$par.fixed, self$tmb_rep()$par.random)
          H <- self$tmb_obj_joint()$he(par_all)
          
          # Extract covariance for random effect components
          ind_re_hid <- which(colnames(Q) == "coeff_re_hid")
          ind_re_obs <- which(colnames(Q) == "coeff_re_obs")
          V_re_hid <- V[ind_re_hid, ind_re_hid]
          V_re_obs <- V[ind_re_obs, ind_re_obs]
          H_re_hid <- H[ind_re_hid, ind_re_hid]
          H_re_obs <- H[ind_re_obs, ind_re_obs]
          
          edf_hid <- sum(diag(H_re_hid %*% V_re_hid))
          edf_obs <- sum(diag(H_re_obs %*% V_re_obs))
          
          # Random effect EDF
          edf <- edf + edf_hid + edf_obs    
        }        
      }
      
      return(edf)
    },
    
    # Model fitting -----------------------------------------------------------
    
    #' @description Suggest initial parameter values
    #' 
    #' Uses K-means clustering to split the data into naive "states", and
    #' estimates observation parameters within each of these states. This is
    #' meant to help with selecting initial parameter values before model
    #' fitting, but users should still think about the values carefully,
    #' and try multiple set of initial parameter values to ensure
    #' convergence to the global maximum of the likelihood function.
    #' 
    #' @return List of initial parameters
    suggest_initial = function() {
      self$obs()$suggest_initial()
    },
    
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
        stop(paste0("Not all observation distributions are properly created. ",
                    "If you added any custom distributions recently, make sure ",
                    "you re-create them using Dist$new after you have ",
                    "recompiled the package."))
      }
      
      # Number of states
      n_states <- self$hid()$nstates()
      
      # Create model matrices of observation process
      # (Design matrices for fixed and random effects, and smoothing matrix)
      mod_mat_obs <- self$obs()$terms()
      X_fe_obs <- mod_mat_obs$X_fe
      X_re_obs <- mod_mat_obs$X_re
      S_obs <- mod_mat_obs$S
      ncol_re_obs <- mod_mat_obs$ncol_re
      
      # Create model matrices of hidden state process
      # (Design matrices for fixed and random effects, and smoothing matrix)
      mod_mat_hid <- self$hid()$terms()
      X_fe_hid <- mod_mat_hid$X_fe
      X_re_hid <- mod_mat_hid$X_re
      S_hid <- mod_mat_hid$S
      ncol_re_hid <- mod_mat_hid$ncol_re
      
      # Prepare initial distribution delta0
      ldelta0 <- self$hid()$delta0(log = TRUE, as_matrix = FALSE)
      
      # Setup TMB parameters
      tmb_par <- list(coeff_fe_obs = self$obs()$coeff_fe(),
                      log_lambda_obs = 0,
                      coeff_fe_hid = self$hid()$coeff_fe(),
                      log_lambda_hid = 0,
                      log_delta0 = ldelta0,
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
        S_obs <- as_sparse(matrix(0, 1, 1))
        ncol_re_obs <- matrix(-1, nr = 1, nc = 1)
        X_re_obs <- as_sparse(rep(0, nrow(X_fe_obs)))
      } else {
        # If there are random effects, 
        # set initial values for coeff_re and log_lambda
        random <- c(random, "coeff_re_obs")
        tmb_par$coeff_re_obs <- rep(0, ncol(X_re_obs))
        tmb_par$log_lambda_obs <- log(self$lambda()$obs)
      }
      
      # Setup random effects in hidden state model
      if(is.null(S_hid)) {
        # If there are no random effects, 
        # coeff_re and log_lambda are not estimated
        map <- c(map, list(coeff_re_hid = factor(NA),
                           log_lambda_hid = factor(NA)))
        S_hid <- as_sparse(matrix(0, 1, 1))
        ncol_re_hid <- matrix(-1, nr = 1, nc = 1)
        X_re_hid <- as_sparse(rep(0, nrow(X_fe_hid)))
      } else {
        # If there are random effects, 
        # set initial values for coeff_re and log_lambda
        random <- c(random, "coeff_re_hid")
        tmb_par$coeff_re_hid <- rep(0, ncol(X_re_hid))
        tmb_par$log_lambda_hid <- log(self$lambda()$hid)
      }
      
      # Add custom mapping effects, i.e., add them to the map argument
      # for TMB, and create a vector of 0/1 to record which parameters
      # are estimated and which are not (used e.g. in post_coeff)
      fixpar <- c(self$hid()$fixpar(all = TRUE), self$obs()$fixpar(all = TRUE))
      par_list <- self$coeff_list()
      usernms <- c("obs", "lambda_obs", "hid", "lambda_hid", "delta0", NA, NA)
      par_names <- names(par_list)
      fixpar_vec <- NULL
      # Loop over model components
      for (i in seq_along(par_list)) {
        # Vector of parameters
        v <- par_list[[par_names[i]]]
        
        # Map vector for TMB
        tmp <- seq_along(v) 
        # Check if user-specified constraint
        fixed <- fixpar[[usernms[i]]]
        mode(fixed) <- "integer"
        if (length(fixed) > 0) {
          # Increase fixed to make sure it's not between 1:length(v)
          fixed <- fixed + length(v)
          if(is.matrix(v)) {
            nms <- rownames(v)
          } else {
            nms <- names(v)
          }
          # Set map to user input
          tmp[nms %in% names(fixed)] <- fixed
          tmp <- factor(as.vector(tmp), levels = unique(as.vector(tmp)))
          ls <- list(tmp)
          names(ls) <- par_names[i]
          map <- c(map, ls)
          # # Which parameters are fixed
          # if(any(is.na(tmp))) {
          #   fix_or_not[which(is.na(tmp))] <- 1
          # }
        }
        fix_or_not <- as.numeric(tmp) + 
          ifelse(test = is.null(fixpar_vec) | all(is.na(fixpar_vec)), 
                 yes = 0, 
                 no = max(fixpar_vec, na.rm = TRUE))
        names(fix_or_not) <- rep(par_names[i], length(v))
        fixpar_vec <- c(fixpar_vec, fix_or_not)
      }
      coeff_array <- cbind(fixpar_vec, unlist(par_list, use.names = FALSE))
      colnames(coeff_array) <- c("fixed", "value")
      private$coeff_array_ <- coeff_array
      
      # Check if delta0 should be stationary
      statdist <- ifelse(self$hid()$stationary(), yes = 1, no = 0)
      
      # Get stored priors 
      priors <- self$priors() 
      
      # Get variables for observation distributions
      obsvar <- self$obs()$obs_var(expand = TRUE)
      datadim <- attr(obsvar, "datadim")
      
      # Data for TMB
      tmb_dat <- list(ID = self$obs()$data()$ID,
                      data = as.matrix(obsvar),
                      datadim = datadim, 
                      known_states = self$obs()$known_states(), 
                      n_states = n_states,
                      statdist = statdist, 
                      distcode = distcode,
                      distpar = distpar, 
                      X_fe_obs = as_sparse(X_fe_obs),
                      X_re_obs = as_sparse(X_re_obs),
                      S_obs = as_sparse(S_obs),
                      ncol_re_obs = ncol_re_obs,
                      X_fe_hid = as_sparse(X_fe_hid),
                      X_re_hid = as_sparse(X_re_hid),
                      S_hid = as_sparse(S_hid),
                      ncol_re_hid = ncol_re_hid,
                      include_smooths = 1, 
                      ref_tpm = self$hid()$ref(),
                      coeff_fe_obs_prior = priors$coeff_fe_obs, 
                      coeff_fe_hid_prior = priors$coeff_fe_hid, 
                      log_lambda_obs_prior = priors$log_lambda_obs, 
                      log_lambda_hid_prior = priors$log_lambda_hid)
      
      # Create TMB model
      obj <- MakeADFun(tmb_dat, tmb_par, DLL = "hmmTMB", 
                       random = random,
                       map = map, 
                       silent = silent)
      
      nllk0 <- obj$fn(obj$par)
      if(is.nan(nllk0) | is.infinite(nllk0)) {
        stop(paste("Log-likelihood is NaN or infinite at starting parameters.",
                   "Check that the data are within the domain of definition of the",
                   "observation distributions, or try other starting parameters."))
      }
      
      # Negative log-likelihood function
      private$tmb_obj_ <- obj
      
      # Joint negative log-likelihood function
      tmb_dat$include_smooths <- -1
      private$tmb_obj_joint_ <- MakeADFun(tmb_dat, tmb_par, DLL = "hmmTMB", 
                                          map = map, silent = silent)
    },
    
    #' @description Fit model using tmbstan
    #' 
    #' Consult documentation of the tmbstan package for more information.
    #' After this method has been called, the Stan output can be accessed
    #' using the method \code{out_stan()}. This Stan output can for example 
    #' be visualised using functions from the rstan package. The parameters
    #' stored in this HMM object are automatically updated to the mean 
    #' posterior estimate, although this can be changed using 
    #' \code{update_par()}.
    #'
    #' @param ... Arguments passed to tmbstan
    #' @param silent Logical. If FALSE, all tracing outputs are shown (default).
    fit_stan = function(..., silent = FALSE) {
      self$formulation()
      
      if (!requireNamespace("rstan", quietly = TRUE)) {
        stop("You need to install the package rstan to do this")
      }
      
      if (!requireNamespace("tmbstan", quietly = TRUE)) {
        stop("You need to install the package tmbstan to do this")
      }
      
      # Setup if necessary
      if(is.null(private$tmb_obj_)) {
        self$setup(silent = silent)
      }
      
      # Run Stan iterations 
      private$out_stan_ <- tmbstan(obj = private$tmb_obj_, init = "par", ...)
      post <- as.matrix(private$out_stan_)
      # Remove "lp__" column
      post <- post[,-ncol(post)]
      
      # Extract posterior samples 
      n_coeff <- nrow(self$coeff_array())
      n_post <- nrow(post)
      iters <- matrix(rep(self$coeff_array()[,"value"], each = n_post), 
                      nrow = n_post, ncol = n_coeff)
      colnames(iters) <- rownames(self$coeff_array())
      # Fill non-fixed columns with posterior samples
      # Do I need to account for shared parameter values? See post_coeff
      iters[,which(!is.na(self$coeff_array()[,"fixed"]))] <- post
      private$iters_ <- iters
      
      # Get iterations on response scale 
      # (only obs parameters and transition probs)
      npar <- length(unlist(self$par()))
      par_iters <- matrix(0, nrow = n_post, ncol = npar)
      for (i in 1:n_post) {
        self$update_par(iter = i)
        par <- self$par()
        # Parameter matrices are transposed for best column order in par_iters
        par_iters[i,] <- unlist(lapply(par, function(pararray) 
          t(pararray[,,1])))
      }
      
      # Name columns for readability
      obspar_names <- names(unlist(self$obs()$formulas()))
      n_states <- self$hid()$nstates()
      tr_names <- paste0(paste0("S", rep(1:n_states, each = n_states), 
                                ">S", rep(1:n_states, n_states)))
      colnames(par_iters) <- c(obspar_names, tr_names)
      private$par_iters_ <- par_iters
      
      # set coefficients to posterior means 
      self$update_par(iter = "mean")
      
    }, 
    
    #' @description Model fitting
    #' 
    #' The negative log-likelihood of the model is minimised using the
    #' function \code{optimx()}. TMB uses the Laplace approximation to integrate 
    #' the random effects out of the likelihood.
    #' 
    #' After the model has been fitted, the output of \code{optimx()} can be
    #' accessed using the method \code{out()}. The estimated parameters can
    #' be accessed using the methods \code{par()} (for the HMM parameters, 
    #' possibly dependent on covariates), \code{predict()} (for uncertainty
    #' quantification and prediction of the HMM parameters for new covariate 
    #' values), \code{coeff_fe()} (for estimated fixed effect coefficients on
    #' the linear predictor scale), and \code{coeff_re()} (for estimated random
    #' effect coefficients on the linear predictor scale).
    #' 
    #' @param silent Logical. If FALSE, all tracing outputs are shown (default).
    #' @param ... Other arguments to optimx which is used to optimise likelihood, 
    #' see ?optimx
    #'
    #' @examples
    #' # Load data set (included with R)
    #' data(nottem)
    #' data <- data.frame(temp = as.vector(t(nottem)))
    #' 
    #' # Create hidden state and observation models
    #' hid <- MarkovChain$new(data = data, n_states = 2)
    #' par0 <- list(temp = list(mean = c(40, 60), sd = c(5, 5)))
    #' obs <- Observation$new(data = data, n_states = 2, 
    #'                        dists = list(temp = "norm"),
    #'                        par = par0)
    #' 
    #' # Create HMM
    #' hmm <- HMM$new(hid = hid, obs = obs)
    #'
    #' # Fit HMM
    #' hmm$fit(silent = TRUE)
    fit = function(silent = FALSE, ...) {
      if(!silent) self$formulation()
      
      # Setup if necessary
      if(is.null(private$tmb_obj_)) {
        self$setup(silent = silent)
      }
      
      # Number of states
      n_states <- self$hid()$nstates()
      
      # Fit model
      args <- list(...)
      if (any(c("par", "fn", "gr", "he", "hessian") %in% names(args))) {
        stop(paste("Cannot supply arguments to HMM$fit() with name par,",
                   "fn, gr, he, or hessian. These are reserved by TMB."))
      }
      if (any(c("lower", "upper") %in% names(args))) {
        warning("'lower' and 'upper' arguments to optimx are ignored in hmmTMB")
      }
      # change default method to nlminb
      private$tmb_obj_$method <- "nlminb"
      if ("method" %in% names(args)) {
        private$tmb_obj_$method <- args$method 
        args <- args[which(names(args) != "method")]
      }
      if (!("control" %in% names(args))) {
        args$control <- list(kkt = FALSE, 
                             starttests = FALSE,
                             dowarn = FALSE)
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
        warning(paste("Convergence code was not zero, indicating that the",
                      "optimizer may not have converged to the correct",
                      "estimates. Please check by consulting the out() function",
                      "which shows what optimx returned."))
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
      if (is.null(private$out_)) stop("Fit model first")
      par_list <- as.list(self$tmb_rep(), "Estimate")
      return(par_list)
    }, 
    
    
    # Forward-backward probabilities ------------------------------------------
    
    #' @description Forward-backward algorithm 
    #' 
    #' The forward probability for time step t and state j 
    #' is the joint pdf/pmf of observations up to time t and of being in
    #' state j at time t, p(Z[1], Z[2], ..., Z[t], S[t] = j).
    #' The backward probability for time t and state j is the
    #' conditional pdf/pmf of observations between time t + 1 and n,
    #' given state j at time t, p(Z[t+1], Z[t+2], ..., Z[n] | S[t] = j).
    #' This function returns their logarithm, for use in other methods
    #' \code{state_probs}, and \code{sample_states}.
    #' 
    #' @return Log-forward and log-backward probabilities 
    forward_backward = function() {
      delta0 <- self$hid()$delta0() 
      n <- nrow(self$obs()$data())
      n_by_ID <- as.numeric(table(self$obs()$data()$ID))
      
      # State-dependent probabilities 
      obsprobs <- self$obs()$obs_probs()
      # Transition probability matrices 
      tpms <- self$hid()$tpm(t = "all")
      
      # Initialise log-forward/backward probabilities
      lforw <- matrix(0, nr = self$hid()$nstates(), nc = n)
      lback <- matrix(0, nr = self$hid()$nstates(), nc = n)
      
      # Loop over ID (time series)
      k <- 1 
      for (ind in 1:length(n_by_ID)) {
        # Forward algorithm 
        p <- delta0[ind,] * obsprobs[k,]
        psum <- sum(p)
        llk <- log(psum)
        p <- p / psum
        lforw[, k] <- log(p) + llk
        for (i in 2:n_by_ID[ind]) {
          p <- p %*% tpms[, , k + i - 2] * obsprobs[k + i - 1, ]
          psum <- sum(p)
          llk <- llk + log(psum)
          p <- p / psum
          lforw[, k + i - 1] <- log(p) + llk 
        }
        
        # Backward algorithm
        lback[, k + n_by_ID[ind] - 1] <- rep(0, self$hid()$nstates())
        p <- rep(1 / self$hid()$nstates(), self$hid()$nstates())
        llk <- log(self$hid()$nstates())
        for (i in (n_by_ID[ind] - 1):1) {
          p <- tpms[, , k + i - 1] %*% (obsprobs[k + i, ] * p)
          lback[, k + i - 1] <- log(p) + llk
          psum <- sum(p)
          p <- p / psum
          llk <- llk + log(psum)
        }
        
        k <- k + n_by_ID[ind]
      }
      
      return(list(logforward = lforw, logbackward = lback))
    }, 
    
    # Pseudo-residuals --------------------------------------------------------
    #' @description Pseudo-residuals
    #'
    #' Compute pseudo-residuals for the fitted model. If the fitted model
    #' is the "true" model, the pseudo-residuals follow a standard normal 
    #' distribution. Deviations from normality suggest lack of fit.
    #' 
    #' @return List (of length the number of variables), where each element is
    #' a vector of pseudo-residuals (of length the number of data points)
    pseudores = function() {
      ID <- self$obs()$data()$ID
      n <- length(ID)
      # Matrix of CDFs at observations
      cat("Computing CDFs... ")
      cdfs <- self$obs()$cdf()
      cat("done\n")
      n_var <- length(cdfs)
      # Transition probability matrices
      tpms <- self$hid()$tpm(t = "all")
      # Log-forward probabilities
      logforw <- self$forward_backward()$logforward
      
      # Loop over observed variables
      res_ls <- list()
      for(var in 1:n_var) {
        # Vector of residuals for this variable
        res <- rep(NA, n)
        
        if(all(is.na(cdfs[[var]]))) {
          message(paste0("Pseudo-residuals not implemented for '", 
                         self$obs()$dists()[[2]]$name(), "' distribution. ",
                         "Returning NA."))
        } else {
          cat("Computing residuals for", names(cdfs)[var], "... ")
          
          # Loop over time steps
          res[1] <- qnorm(t(self$hid()$delta0()[1,]) %*% cdfs[[var]][1,])
          for(i in 2:n) {
            if(ID[i] != ID[i-1]) {
              res[i] <- qnorm(t(self$hid()$delta0()[ID[i],]) %*% cdfs[[var]][i,])
            } else {
              # c cancels out below (to avoid numerical issues)
              c <- max(logforw[, i - 1])
              a <- exp(logforw[, i - 1] - c)
              res[i] <- qnorm(t(a) %*% (tpms[,,i]/sum(a)) %*% cdfs[[var]][i,])
            }
          }
          cat("done\n")
        }
        
        res_ls[[var]] <- res
      }
      names(res_ls) <- names(cdfs)
      
      return(res_ls)
    },
    
    # State estimation --------------------------------------------------------
    
    #' @description Viterbi algorithm
    #' 
    #' @return Most likely state sequence
    viterbi = function() {
      if(!is.null(private$states_)) {
        return(private$states_)
      }
      
      data <- self$obs()$data()
      ID <- data$ID
      
      # Number of observations
      n <- nrow(data)
      # Number of states
      n_states <- self$hid()$nstates()
      # Number of variables
      n_var <- length(self$obs()$dists())
      
      # Observation probabilities
      obs_probs <- self$obs()$obs_probs()
      
      # Transition probability matrices      
      tpm_all <- self$hid()$tpm(t = "all")
      
      # Number of unique IDs
      n_id <- length(unique(ID))
      # First index for each ID
      i0 <- c(1, which(ID[-1] != ID[-n]) + 1, n + 1)
      
      # Initialise state sequence
      all_states <- rep(NA, length = n)
      
      # Initial distribution
      delta0 <- self$hid()$delta0() 
      
      # Loop over IDs
      for(id in 1:n_id) {
        # Subset to this ID
        ind_this_id <- i0[id]:(i0[id + 1] - 1)
        sub_obs_probs <- obs_probs[ind_this_id,]
        sub_tpm_all <- tpm_all[,,ind_this_id]
        
        # Number of observations for this ID
        n_this_id <- length(ind_this_id)
        
        # Forward iterations
        xi <- matrix(NA, n_this_id, n_states)
        v <- delta0[id,] * sub_obs_probs[1,]
        xi[1,] <- v/sum(v)
        for(i in 2:n_this_id) {
          v <- apply(xi[i-1,] * sub_tpm_all[,,i-1], 2, max) * sub_obs_probs[i,]
          xi[i,] <- v/sum(v)
        }
        
        # Backward iterations
        states <- rep(NA, n_this_id)
        states[n_this_id] <- which.max(xi[n_this_id,])
        for(i in (n_this_id - 1):1) {
          states[i] <- which.max(sub_tpm_all[, states[i+1], i] * xi[i,])
        }
        
        # Append estimated states for this ID
        all_states[ind_this_id] <- states
      }
      
      # Save state sequence
      private$states_ <- all_states
      
      return(all_states)
    },
    
    #' @description Sample posterior state sequences using forward-filtering
    #' backward-sampling 
    #' 
    #' The forward-filtering backward-sampling algorithm returns a
    #' sequence of states, similarly to the Viterbi algorithm, but it generates
    #' it from the posterior distribution of state sequences, i.e., accounting
    #' for uncertainty in the state classification. Multiple generated sequences
    #' will therefore generally not be the same.
    #' 
    #' @param nsamp Number of samples to produce 
    #' @param full If TRUE and model fit by fit_stan then parameter estimates are 
    #' sampled from the posterior samples before simulating each sequence 
    #' 
    #' @return Matrix where each column is a different sample of state sequences,
    #' and each row is a time of observation
    sample_states = function(nsamp = 1, full = FALSE) {
      if (full & is.null(private$out_stan_)) {
        stop("Fit model with fit_stan() first")
      }
      if (full) {
        n_full_samp <- nsamp
        nsamp <- 1 
      } else {
        n_full_samp <- 1
      }
      nstates <- self$hid()$nstates()
      actual_samps <- max(n_full_samp, nsamp)
      n <- nrow(self$obs()$data())
      states <- matrix(0, nr = n, nc = actual_samps)
      
      # Loop over posterior draws if full = TRUE
      for (k in 1:n_full_samp) {
        if (full) {
          # Sample a parameter iteration at random
          self$update_par(iter = sample(1:nrow(mod$iters()), size = 1))
        }
        
        # Get forward-backward probabilities 
        fb <- self$forward_backward()
        # Get tpms 
        tpms <- self$hid()$tpm(t = "all")
        # Get observation probabilities
        obsprobs <- self$obs()$obs_probs()
        # Get number of observations per individual 
        n_by_ID <- as.numeric(table(self$obs()$data()$ID))
        cumn <- cumsum(n_by_ID)
        
        # Loop over time series / individuals
        for (ind in 1:length(n_by_ID)) {
          m <- cumn[ind]
          
          # Sample last state
          L <- logsumexp(fb$logforward[,m])
          prob <- exp(fb$logforward[,m] - L)
          prob <- prob / sum(prob)
          if (full) {
            states[m,k] <- sample(1:nstates, prob = prob, 
                                  size = nsamp, replace = TRUE)
          } else {
            states[m,] <- sample(1:nstates, prob = prob, 
                                 size = nsamp, replace = TRUE)
          }
          
          # Sample backward
          for (s in 1:nsamp) {
            for (i in 1:(n_by_ID[ind] - 1)) {
              lprob <- fb$logforward[, m - i] + 
                log(tpms[, states[m - i + 1, s], m - i]) + 
                log(obsprobs[m - i + 1, states[m - i + 1, s]]) + 
                fb$logbackward[states[m - i + 1, s], m - i + 1] - L 
              lprob <- lprob - logsumexp(lprob)
              prob <- exp(lprob)
              index <- ifelse(full, k, s)
              states[m - i, index] <- sample(1:nstates, prob = prob, size = 1)
            }
          }
        }
      }
      return(states)
    }, 
    
    #' @description Compute posterior probability of being in each state 
    #' 
    #' @return matrix with a row for each observation and a column for each state 
    state_probs = function() {
      n <- nrow(self$obs()$data())
      nstates <- self$hid()$nstates()
      n_by_ID <- as.numeric(table(self$obs()$data()$ID))
      cn <- cumsum(n_by_ID)
      fb <- self$forward_backward()
      pr_state <- matrix(0, nr = n, nc = nstates)
      k <- 0 
      for (ind in 1:length(n_by_ID)) {
        llk <- logsumexp(fb$logforward[,cn[ind]])
        for (i in 1:n_by_ID[ind]) {
          pr_state[k + i,] <- exp(fb$logforward[, k + i] + fb$logbackward[,k + i] - llk)
        }
        k <- k + n_by_ID[ind]
      }
      colnames(pr_state) <- paste0("state", 1:nstates)
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
        V <- prec_to_cov(rep$jointPrecision)
      }
      
      # Generate samples from MVN estimator distribution
      post <- rmvn(n = n_post, mu = par, V = V)
      
      # Matrix filled with estimates
      npar <- nrow(self$coeff_array())
      post_all <- matrix(rep(self$coeff_array()[,"value"], each = n_post), 
                         nrow = n_post, ncol = npar)
      colnames(post_all) <- rownames(self$coeff_array())
      
      # Fill non-fixed columns with posterior samples
      post_all[,which(!is.na(self$coeff_array()[,"fixed"]))] <- 
        post[,na.omit(self$coeff_array()[,"fixed"])]
      
      return(post_all)
    },
    
    #' @description Posterior sampling for linear predictor 
    #' 
    #' @param n_post Number of posterior samples
    #' 
    #' @return List with elements obs and hid, where each is a matrix 
    #' with one column for each predictor and one row for each posterior draw
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
      lp$hid <- matrix(0, nr = n_post, nc = nrow(self$hid()$X_fe()))
      for (i in 1:n_post) {
        self$obs()$update_coeff_fe(obspars_fe[i,])
        self$obs()$update_coeff_re(obspars_re[i,])
        self$hid()$update_coeff_fe(hidpars_fe[i,])
        self$hid()$update_coeff_re(hidpars_re[i,])
        lp$obs[i,] <- self$obs()$linpred()
        lp$hid[i,] <- self$hid()$linpred()
      }
      
      # reset design matrices and parameters
      self$obs()$update_coeff_fe(coeff_fe_old$obs)
      self$obs()$update_coeff_re(coeff_re_old$obs)
      self$hid()$update_coeff_fe(coeff_fe_old$hid)
      self$hid()$update_coeff_re(coeff_re_old$hid)
      
      return(lp)
    }, 
    
    #' @description Create posterior simulations of a function of a model component 
    #' 
    #' @param fn Function which takes a vector of linear predictors as input
    #' and produces either a scalar or vector output 
    #' @param n_post Number of posterior simulations 
    #' @param comp Either "obs" for observation model linear predictor, or
    #' "hid" for hidden model linear predictor 
    #' @param ... Arguments passed to fn
    #' @param level Confidence interval level if required (e.g., 0.95 for 95%
    #' confidence intervals). Default is 0, i.e., confidence intervals are not
    #' returned. 
    #' @param return_post Logical indicating whether to return the posterior
    #' samples. If FALSE (default), only mean estimates and confidence intervals
    #' are returned
    #' 
    #' @return A list with elements:
    #' \describe{
    #'   \item{post}{If return_post = TRUE, this is a vector (for scalar 
    #'   outputs of fn) or matrix (for vector outputs) with a column for 
    #'   each simulation}
    #'   \item{mean}{Mean over posterior samples}
    #'   \item{lcl}{Lower confidence interval bound (if level !=0)}
    #'   \item{ucl}{Upper confidence interval bound (if level !=0)}
    #' }
    #' 
    post_fn = function(fn, n_post, comp = NULL, ..., level = 0, return_post = FALSE) {
      # Get linear predictors
      lp <- self$post_linpred(n_post)
      # Output list
      res <- NULL
      
      # Compute function of linear predictor
      res$post <- lapply(1:nrow(lp[[comp]]), FUN = function(i) {
        fn(linpred = lp[[comp]][i,], ...)
      })
      
      # Compute means 
      res$mean <- Reduce("+", res$post) / length(res$post)
      
      # Compute confidence interval
      if (level > 0) {
        alpha <- 1 - level
        arr <- simplify2array(res$post)
        ci <- apply(X = arr, 
                    MARGIN = 1:(length(dim(arr)) - 1), 
                    FUN = quantile, 
                    prob = c(alpha/2, 1 - alpha/2))
        nci <- length(dim(ci))
        ci <- aperm(ci, c(2:nci, 1))
        block <- prod(dim(ci)[-nci])
        lcl <- ci[1:block]
        ucl <- ci[(1:block) + block]
        dim(lcl) <- dim(ucl) <- dim(ci)[-nci]
        res$lcl <- lcl
        res$ucl <- ucl
        dimnames(res$lcl) <- dimnames(res$mean)
        dimnames(res$ucl) <- dimnames(res$mean)
      }
      
      if(!return_post) {
        res$post <- NULL
      }
      
      return(res)
    }, 
    
    #' @description Predict estimates from a fitted model
    #' 
    #' By default, this returns point estimates of the HMM parameters
    #' for a new data frame of covariates. See the argument `n_post`
    #' to also get confidence intervals.
    #' 
    #' @param what Which estimates to predict? Options include 
    #' transition probability matrices "tpm", 
    #' stationary distributions "delta", or 
    #' observation distribution parameters "obspar"
    #' @param t Time points to predict at 
    #' @param ... Other arguments to the respective functions 
    #' for hid$tpm, hid$delta, obs$par
    #' @param newdata New dataframe to use for prediction
    #' @param n_post If greater than zero then n_post posterior 
    #' samples are produced, and used to create confidence intervals.
    #' @param level Level of the confidence intervals, e.g. CI = 0.95
    #' will produce 95\% confidence intervals (default) 
    #' @param return_post Logical. If TRUE, a list of posterior samples
    #' is returned.
    #' @param as_list Logical. If confidence intervals are required for the
    #' transition probabilities or observation parameters, this
    #' argument determines whether the MLE, lower confidence limit and upper
    #' confidence limit are returned as separate elements in a list (if
    #' TRUE; default), or whether they are combined into a single array (if
    #' FALSE). Ignored if \code{what = "delta"} or if \code{n_post = 0}.
    #' 
    #' @return Maximum likelihood estimates (\code{mle}) of predictions,
    #' and confidence limits (\code{lcl} and \code{ucl}) if requested. The
    #' format of the output depends on whether confidence intervals are
    #' required (specified through \code{n_post}), and on the argument
    #' \code{as_list}.
    #' 
    #' @examples
    #' # Load data set (included with R)
    #' data(nottem)
    #' data <- data.frame(temp = as.vector(t(nottem)))
    #' 
    #' # Create hidden state and observation models
    #' hid <- MarkovChain$new(data = data, n_states = 2)
    #' par0 <- list(temp = list(mean = c(40, 60), sd = c(5, 5)))
    #' obs <- Observation$new(data = data, n_states = 2, 
    #'                        dists = list(temp = "norm"),
    #'                        par = par0)
    #' 
    #' # Create HMM
    #' hmm <- HMM$new(hid = hid, obs = obs)
    #' 
    #' # Fit HMM
    #' hmm$fit(silent = TRUE)
    #' 
    #' # Get transition probability matrix with confidence intervals
    #' hmm$predict(what = "tpm", n_post = 1000)
    predict = function(what, t = 1, newdata = NULL, n_post = 0, level = 0.95, 
                       return_post = FALSE, as_list = TRUE) {
      if (is.null(private$out_) & n_post > 0) {
        stop("Fit model first")
      }
      if (!is.null(newdata)) {
        # Return predictions for all rows of newdata if provided
        t <- "all"
        
        # Save model matrices, then replace by matrices based on newdata
        old <- list(X_fe_obs = self$obs()$X_fe(), 
                    X_re_obs = self$obs()$X_re(), 
                    X_fe_hid = self$hid()$X_fe(), 
                    X_re_hid = self$hid()$X_re())
        obsmats <- self$obs()$make_mat(new_data = newdata)
        hidmats <- self$hid()$make_mat(data = self$obs()$data(), 
                                       new_data = newdata)
        self$obs()$update_X_fe(obsmats$X_fe)
        self$obs()$update_X_re(obsmats$X_re)
        self$hid()$update_X_fe(hidmats$X_fe) 
        self$hid()$update_X_re(hidmats$X_re)
      }
      
      # get appropriate prediction function 
      fn <- switch(what, 
                   tpm = self$hid()$tpm,
                   delta = self$hid()$delta,
                   obspar = self$obs()$par)
      
      # get appropriate model component 
      comp <- switch(what, tpm = "hid", delta = "hid", obspar = "obs")
      
      if (n_post == 0) {
        # just return predicted means if no confidence intervals wanted 
        val <- fn(linpred = self[[comp]]()$linpred(), t = t)
      } else {
        # return means and confidence intervals
        val <- self$post_fn(fn = fn, 
                            n_post = n_post, 
                            comp = comp, 
                            t = t,
                            level = level,
                            return_post = return_post)

        # Replace posterior mean from post_fn() by MLE        
        val$mean <- fn(linpred = self[[comp]]()$linpred(), t = t)
        names(val)[1] <- "mle"
        
        # Format as array for nicer output
        if(!as_list & what %in% c("tpm", "obspar")) {
          mle <- val$mle
          names <- paste0(rep(rownames(mle), each = ncol(mle)), 
                          " - ", rep(colnames(mle), nrow(mle)))
          a <- array(NA, dim = c(length(names), 3, dim(mle)[3]), 
                     dimnames = list(names, c("mle", "lcl", "ucl"), NULL))
          a[,1,] <- as.vector(aperm(mle, c(2, 1, 3)))
          a[,2,] <- as.vector(aperm(val$lcl, c(2, 1, 3)))
          a[,3,] <- as.vector(aperm(val$ucl , c(2, 1, 3)))
          val <- a
        }
      }
      
      # Reset to original model matrices
      if (!is.null(newdata)) {
        self$obs()$update_X_fe(old$X_fe_obs)
        self$hid()$update_X_fe(old$X_fe_hid)
        self$obs()$update_X_re(old$X_re_obs)
        self$hid()$update_X_re(old$X_re_hid)
      }
      
      return(val)
    }, 
    
    #' @description Confidence intervals for working parameters
    #' 
    #' This function computes standard errors for all fixed effect model
    #' parameters based on the diagonal of the inverse of the Hessian matrix,
    #' and then derives Wald-type confidence intervals.
    #' 
    #' @param level Level of confidence intervals. Defaults to 0.95, i.e., 95\%
    #' confidence intervals.
    #' 
    #' @return List of matrices with three columns: mle (maximum likelihood 
    #' estimate), lcl (lower confidence limit), ucl (upper confidence
    #' limit), and se (standard error). One such matrix is produced for 
    #' the working parameters of the observation model, the working parameters
    #' of the hidden state model, the smoothness parameters of the observation 
    #' model, and the smoothness parameters of the hidden state model.
    confint = function(level = 0.95) {
      # Get standard errors from covariance matrix
      rep <- self$tmb_rep()
      par <- rep$par.fixed
      if(is.null(rep$jointPrecision)) {
        V <- rep$cov.fixed
      } else {
        V <- prec_to_cov(rep$jointPrecision)
      }
      se <- sqrt(diag(V)[1:length(par)])
      
      # Unpack model components
      obspar_se <- se[which(names(par) == "coeff_fe_obs")]
      hidpar_se <- se[which(names(par) == "coeff_fe_hid")]
      obslam_se <- se[which(names(par) == "log_lambda_obs")]
      hidlam_se <- se[which(names(par) == "log_lambda_hid")]
      
      # Get Wald-type confidence intervals
      quant <- qnorm(1 - (1 - level)/2)
      obspar_ci <- cbind(mle = self$coeff_fe()$obs,
                         lcl = self$coeff_fe()$obs - quant * obspar_se,
                         ucl = self$coeff_fe()$obs + quant * obspar_se,
                         se = obspar_se)
      hidpar_ci <- cbind(mle = self$coeff_fe()$hid,
                         lcl = self$coeff_fe()$hid - quant * hidpar_se,
                         ucl = self$coeff_fe()$hid + quant * hidpar_se,
                         se = hidpar_se)
      obslam_ci <- cbind(mle = self$lambda()$obs,
                         lcl = self$lambda()$obs - quant * obslam_se,
                         ucl = self$lambda()$obs + quant * obslam_se,
                         se = obslam_se)
      hidlam_ci <- cbind(mle = self$lambda()$hid,
                         lcl = self$lambda()$hid - quant * hidlam_se,
                         ucl = self$lambda()$hid + quant * hidlam_se,
                         se = hidlam_se)
      
      out <- list(coeff_fe = list(obs = obspar_ci, hid = hidpar_ci),
                  lambda = list(obs = obslam_ci, hid = hidlam_ci)) 
      return(out)
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
      n_states <- self$hid()$nstates()
      
      # Simulate state process      
      S <- self$hid()$simulate(n = n, data = self$obs()$data(), 
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
            if(!silent) {
              cat("\rSimulating ", var_name, "... ", 
                  round(i/n*100), "%", sep = "")
            }
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
    #' Many time series are simulated from the fitted model, and the
    #' statistic(s) of interest are calculated for each. A histogram of 
    #' those values can for example be used to compare to the observed
    #' value of the statistic. An observation far in the tails of the
    #' distribution of simulated statistics suggests lack of fit.
    #' 
    #' @param check_fn Goodness-of-fit function which accepts "data" as input
    #' and returns a statistic (either a vector or a single number) to be
    #' compared between observed data and simulations. 
    #' @param nsims Number of simulations to perform 
    #' @param full If model fitted with `fit_stan`, then full = TRUE
    #'  will sample from posterior for each simulation 
    #' @param silent Logical. If FALSE, simulation progress is shown. 
    #' (Default: TRUE)
    #' 
    #' @return List with elements:
    #' \describe{
    #'   \item{obs_stat: }{Vector of values of goodness-of-fit statistics for the
    #'   observed data}
    #'   \item{stats: }{Matrix of values of goodness-of-fit statistics for the
    #'   simulated data sets (one row for each statistic, and one column for each
    #'   simulation)}
    #'   \item{plot: }{ggplot object}
    #' }
    check = function(check_fn, nsims = 100, full = FALSE, silent = FALSE) {
      # Evaluate statistics for observed data
      obs_stat <- check_fn(self$obs()$data())
      
      # Simulate from model and evaluate statistics for simulated data
      stats <- matrix(0, nc = nsims, nr = length(obs_stat))
      for (sim in 1:nsims) {
        if (!silent) cat("Simulating", sim, " / ", nsims, "\r")
        
        # if full and mcmc then sample parameter
        if (full & !is.null(private$out_stan_)) {
          self$update_par(iter = sample(1:nrow(self$iters()), size = 1))
        }
        
        # simulate new data
        newdat <- self$simulate(n = nrow(self$obs()$data()), 
                                data = self$obs()$data(),
                                silent = TRUE) 
        # compute statistics
        stats[,sim] <- check_fn(newdat)
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
      if (!silent) plot(p)
      
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
    #' @param var Name of the variable to plot.
    #' @param var2 Optional name of a second variable, for 2-d plot.
    #' @param line Logical. If TRUE (default), lines are drawn between
    #' successive data points. Can be set to FALSE if another geom is
    #' needed (e.g., geom_point).
    #' 
    #' @return A ggplot object
    plot_ts = function(var, var2 = NULL, line = TRUE) {
      if(is.null(private$states_)) {
        self$viterbi()
      }
      
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
                         x = data[[var]])
        
        p <- ggplot(data = df, mapping = aes(index, x, col = state, group = ID)) +
          xlab("time") + ylab(var)
        if(line) {
          p <- p + geom_line()          
        }
      } else {
        # 2d plot
        df <- data.frame(ID = ID,
                         x = data[[var]],
                         y = data[[var2]])
        
        p <- ggplot(data = df, mapping = aes(x, y, col = state, group = ID)) +
          xlab(var) + ylab(var2)
        if(line) {
          p <- p + geom_path()          
        }
      }
      
      p <- p + 
        scale_color_manual(values = hmmTMB_cols) +
        theme_light()
      
      return(p)
    },
    
    #' @description Plot observation distributions weighted by frequency in Viterbi 
    #' 
    #' This is a wrapper around Observation$plot_dist, where the
    #' distribution for each state is weighted by the proportion of time
    #' spent in that state (according to the Viterbi state sequence).
    #'
    #' @param var Name of data variable
    #' 
    #' @return Plot of distribution with data histogram 
    plot_dist = function(var) {
      if(is.null(private$states_)) {
        self$viterbi()
      }
      # Proportion of time spent in each state
      weights <- sapply(1:self$hid()$nstates(), function(s) 
        length(which(self$states() == s))/length(self$states()))
      return(self$obs()$plot_dist(var, weights = weights))
    },
    
    #' @description Plot a model component 
    #' 
    #' @param what Name of model component to plot: should be one of "tpm"
    #' (transition probabilities), "delta" (stationary state probabilities), 
    #' or "obspar" (state-dependent observation parameters) 
    #' @param var Name of covariate to plot on x-axis 
    #' @param covs Optional named list for values of covariates (other than 'var') 
    #' that should be used in the plot (or dataframe with single row). If this is
    #' not specified, the mean value is used for numeric variables, and the
    #' first level for factor variables.
    #' @param i If plotting tpm then rows of tpm; if plotting delta then indices
    #' of states to plot; if plotting obspar then full names of parameters 
    #' to plot (e.g., obsvar.mean) 
    #' @param j If plotting tpm then columnss of tpm to plot; if plotting delta 
    #' then this is ignored,; if plotting obspar then indices of states to plot 
    #' @param n_grid Number of points in grid over x-axis (default: 50) 
    #' @param n_post Number of posterior simulations to use when computing
    #' confidence intervals; default: 1000. See \code{predict} function for 
    #' more detail.
    #' 
    #' @return A ggplot object
    plot = function(what, var = NULL, covs = NULL, i = NULL, j = NULL, 
                    n_grid = 50, n_post = 1000) {
      # Get relevant model component 
      comp <- switch(what, tpm = "hid", delta = "hid", obspar = "obs")
      # Get x-axis 
      newdata <- cov_grid(var = var, 
                          obj = self, 
                          covs = covs, 
                          formulas = self[[comp]]()$formulas(), 
                          n_grid = n_grid)
      
      # Number of states
      n_states <- self$hid()$nstates()
      
      # Get predictions and uncertainty 
      preds <- self$predict(what = what, t = "all", newdata = newdata, 
                            level = 0.95, n_post = n_post)
      
      # Data frame for plot
      df <- as.data.frame.table(preds$mle)
      df$lcl <- as.vector(preds$lcl)
      df$ucl <- as.vector(preds$ucl)
      if (what == "tpm") {
        colnames(df) <- c("from", "to", "var", "prob", "lcl", "ucl")
        levels(df$from) <- paste("State", 1:n_states)
        levels(df$to) <- paste("State", 1:n_states)
        df$var <- rep(newdata[, var], each = n_states * n_states)
        if (!is.null(i)) df <- df[df$from == paste0("State ", i),]
        if (!is.null(j)) df <- df[df$to == paste0("State ", j),]
      } else if (what == "delta") {
        colnames(df) <- c("var", "state", "prob", "lcl", "ucl")
        levels(df$state) <- paste("State", 1:n_states)
        df$var <- rep(newdata[, var], n_states)
        if (!is.null(i)) df <- df[df$state == paste0("State ", i),]
      } else if (what == "obspar") {
        colnames(df) <- c("par", "state", "var", "val", "lcl", "ucl")
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
      
      if (what == "tpm") {
        p <- ggplot(df, aes(var, prob)) + 
          geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.3) +
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
            geom_errorbar(aes(x = var, ymin = lcl, ymax = ucl), 
                          alpha = 0.5, width = 0.2) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
        } else {
          p <- p + geom_line(linewidth = 0.7) +
            geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.3)
        }
      } else {
        if (what == "delta") {
          p <- ggplot(df, aes(var, prob, group = state, col = state)) +
            scale_color_manual("", values = hmmTMB_cols) +
            scale_fill_manual(values = hmmTMB_cols, guide = "none") +
            xlab(var) + ylab("Stationary state probabilities") + ggtitle(plot_txt) +
            theme_light() + 
            coord_cartesian(ylim = c(0, 1))
        } else if (what == "obspar") {
          p <- ggplot(df, aes(var, val, col = state)) + theme_light() +
            scale_color_manual("", values = hmmTMB_cols) +
            scale_fill_manual(values = hmmTMB_cols, guide = "none") +
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
            geom_errorbar(aes(x = var, ymin = lcl, ymax = ucl), 
                          alpha = 0.5, width = 0.2) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
        } else {
          p <- p + geom_line(linewidth = 0.7) +
            geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = state), 
                        col = NA, alpha = 0.3)
        }
      }
      return(p)
    }, 
    
    # AIC methods (for simulation experiments) --------------------------------
    #' @description Marginal Akaike Information Criterion
    #' 
    #' The marginal AIC is for example defined by 
    #' Wood (2017), as AIC = - 2L + 2k where L is the
    #' maximum marginal log-likelihood (of fixed 
    #' effects), and k is the number of degrees
    #' of freedom of the fixed effect component of
    #' the model
    #' 
    #' @return Marginal AIC
    AIC_marginal = function() {
      llk <- -self$out()$value
      npar <- nrow(self$obs()$coeff_fe()) + 
        nrow(self$hid()$coeff_fe()) +
        length(self$obs()$lambda()) +
        length(self$hid()$lambda())
      if(!self$hid()$stationary()) {
        npar <- npar + length(self$coeff_list()$log_delta0)
      }
      
      if(nrow(self$obs()$coeff_re()) + nrow(self$hid()$coeff_re()) > 0) {
        warning("AIC functions are experimental for models with random effects",
                " or splines. Use at your own risk.")
      }
      
      aic <- - 2 * llk + 2 * npar
      
      return(aic)
    },
    
    #' @description Conditional Akaike Information Criterion
    #' 
    #' The conditional AIC is for example defined by 
    #' Wood (2017), as AIC = - 2L + 2k where L is the
    #' maximum joint log-likelihood (of fixed and random
    #' effects), and k is the number of effective degrees
    #' of freedom of the model (accounting for flexibility
    #' in non-parametric terms implied by smoothing)
    #' 
    #' @return Conditional AIC
    AIC_conditional = function() {
      # Get all estimated parameters (fixed and random)
      par_all <- c(self$tmb_rep()$par.fixed, self$tmb_rep()$par.random)
      
      llk <- - self$tmb_obj_joint()$fn(par_all)
      npar <- self$edf()
      
      if(nrow(self$obs()$coeff_re()) + nrow(self$hid()$coeff_re()) > 0) {
        warning("AIC functions are experimental for models with random effects",
                " or splines. Use at your own risk.")
      }
      
      aic <- - 2 * llk + 2 * npar
      
      return(aic)
    },
    
    # Print methods -----------------------------------------------------------
    #' @description Print observation parameters at t = 1
    print_obspar = function() {
      if(is.null(private$out_) & is.null(private$out_stan_)) {
        cat("> Initial observation parameters (t = 1):\n")
      } else {
        cat("> Estimated observation parameters (t = 1):\n")
      }
      par <- self$par(t = 1)
      print(round(par$obspar[,,1], 3))
      cat("\n")
    },
    
    #' @description Print observation parameters at t = 1
    print_tpm = function() {
      if(is.null(private$out_) & is.null(private$out_stan_)) {
        cat("> Initial transition probabilities (t = 1):\n")
      } else {
        cat("> Estimated transition probabilities (t = 1):\n")
      }
      par <- self$par(t = 1)
      print(round(par$tpm[,,1], 3))
      cat("\n")
    },
    
    #' @description Print model formulation    
    formulation = function() {
      self$obs()$formulation()
      self$hid()$formulation()
    },
    
    #' @description Print HMM object
    print = function() {
      self$obs()$formulation()
      self$print_obspar()
      self$hid()$formulation()
      self$print_tpm()
    }
  ),
  
  private = list(
    
    # Private data members ----------------------------------------------------
    obs_ = NULL,
    hid_ = NULL,
    out_ = NULL,
    tmb_obj_ = NULL,
    tmb_obj_joint_ = NULL,
    tmb_rep_ = NULL,
    priors_ = NULL, 
    out_stan_ = NULL, 
    iters_= NULL,
    par_iters_ = NULL, 
    coeff_array_ = NULL,
    states_ = NULL,
    
    # Reading from spec file --------------------------------------------------
    
    ## @description Read model specification files
    ## 
    ## @param file File location 
    ## 
    ## @return List with elements:
    ## \describe{
    ##   \item{\code{data}}{Data frame}
    ##   \item{\code{nstates}}{Number of states}
    ##   \item{\code{dists}}{List of observation distributions}
    ##   \item{\code{forms}}{Formulas for observation model}
    ##   \item{\code{tpm}}{Formulas of hidden state model}
    ##   \item{\code{par}}{Initial parameters for observation model}
    ##   \item{\code{tpm0}}{Initial transition probability matrix}
    ##   \item{\code{fixed}}{Fixed parameters}
    ## }
    read_file = function(file) {
      if (!file.exists(file)) {
        stop("Model specification file does not exist:", file)
      }
      # Vector with one element for each line of 'file'
      spec <- scan(file = file,
                   character(), 
                   sep = "\n", 
                   strip.white = TRUE, 
                   quiet = TRUE)
      # Remove comments
      spec <- strip_comments(spec)
      # Block names 
      blocknms <- c("DATA", "DISTRIBUTION", "FORMULA", "TPM", "INITIAL", "FIXED")
      # Find blocks
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
        tpm <- matrix(str_trim(unlist(str_split(tpm_block, ";"))), 
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
        len <- length(fixed_block)
        fixed <- NULL
        if ("delta0" %in% fixed_block) {
          fixed$delta0 <- matrix(NA, nrow = length(unique(self$obs()$data()$ID)), 
                                 ncol = nstates - 1)
          colnames(fixed$delta0) <- paste0("state", 1:(nstates - 1))
          len <- len - 1
        }
        fixed$obs <- rep(NA, len)
        names(fixed$obs) <- fixed_block[fixed_block != "delta0"]
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
                  delta0 = ini$delta0, 
                  fixed = fixed))
      
    }, 
    
    ## @description Read a specified block in a model specification file 
    ## 
    ## Subset specification file, passed as vector with one string
    ## for each line, to only keep lines corresponding to specified block.
    read_block = function(name, wh_blocks, spec) {
      find <- which(names(wh_blocks) == name)
      start_block <- wh_blocks[find] + 1 
      end_block <- ifelse(start_block > max(wh_blocks), 
                          length(spec), 
                          wh_blocks[find + 1] - 1)
      block <- spec[start_block:end_block]
      return(block)
    }, 
    
    ## @description Separate left hand side and right hand side variables 
    ## from a equation 
    ## @param x Character vector of equations 
    ## @return List of left hand sides (lhs) and right hand sides (rhs)
    read_equals = function(x) {
      rhs <- str_trim(gsub(".*=", "", x))
      lhs <- str_trim(gsub("=.*", "", x))
      return(list(lhs = lhs, rhs = rhs))
    }, 
    
    ## @description Read distribution block 
    ## @param dist Character vector of distribution block 
    ## @return List of distributions 
    read_dists = function(dists) {
      ds <- strsplit(dists, "\n")
      terms <- sapply(ds, FUN = function(x) {all.vars(as.formula(x))})
      ls <- as.list(terms[2,])
      names(ls) <- terms[1,]
      return(ls)
    }, 
    
    ## @description Read both formula and initial blocks 
    ## @param forms Character vector of block to read
    ## @param ini If TRUE then read as if it is initial block otherwise assume it
    ## is the formula block 
    read_forms = function(forms, ini = FALSE, nstates = NULL) {
      par <- NULL
      tpm0 <- NULL
      delta0 <- NULL
      
      # Find variables 
      wh_vars <- grep(":", forms)
      vars <- str_trim(gsub(":", "", forms[wh_vars]))
      
      # Loop over variables
      for (i in 1:length(vars)){
        # Find sub-block of formula/initial values for that variable 
        if (i < length(vars)) {
          end <- wh_vars[i + 1] - 1
        } else {
          end <- length(forms)
        }
        subforms <- forms[(wh_vars[i] + 1):end]
        
        if (str_trim(vars[i]) == "tpm") {
          # Read transition probability matrix
          nstates <- length(subforms)
          tpm0 <- matrix(0, nr = nstates, nstates)
          for (s in 1:nstates) {
            tpm0[s,] <- as.numeric(strsplit(subforms[s], ",")[[1]])
          }
        } else if (str_trim(vars[i]) == "delta0") {
          # Read initial distribution
          delta0 <- as.numeric(strsplit(subforms, ",")[[1]])
        } else {
          # Read observation parameters
          subpar <- NULL
          subparnms <- NULL
          
          # Loop over parameters
          for (j in 1:length(subforms)) {
            # Get parameter name
            pattern <- ifelse(ini, "\\=.*", "\\~.*")
            par_name <- str_trim(gsub(pattern, "", subforms[j]))
            
            # Get rhs of formula or initial values 
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
      
      names(par) <- vars[!(vars %in% c("tpm", "delta0"))]
      if (ini) {
        res <- list(par = par, tpm0 = tpm0, delta0 = delta0)
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
    ## This is used to initialise parameters based on a previous model
    initialize_submodel = function(par, initpar) {
      # Find which parameters of initializing model occur in current model
      wh <- match(rownames(initpar), rownames(par))
      # Find which parameters of initializing model do not occur in current model
      wh2 <- 1:nrow(initpar)
      wh2 <- wh2[!is.na(wh)]
      wh <- wh[!is.na(wh)]
      # Copy over shared part 
      if (length(wh) > 0) {
        par[wh, 1] <- initpar[wh2, 1]
      }
      return(par)
    },
    
    ## Check constructor arguments 
    # (For argument description, see constructor)
    check_args = function(obs, hid, init) {
      if(!inherits(obs, "Observation")) {
        stop("'obs' should be an Observation object")
      }
      if(!inherits(hid, "MarkovChain")) {
        stop("'hid' should be an MarkovChain object")
      }
      if(obs$nstates() != hid$nstates()) {
        stop(paste0("The observation model and hidden state model should have the ",
                    "same number of states. (obs$nstates = ", obs$nstates(), 
                    ", hid$nstates = ", hid$nstates(), ")"))
      }
      if(!is.null(init)) {
        if (!("HMM" %in% class(init))) {
          stop("init must be a HMM object.")
        }
      }
    }
  )
)
