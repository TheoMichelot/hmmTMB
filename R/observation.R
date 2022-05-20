
#' R6 class for HMM observation model
#'
#' Contains the data, distributions, parameters, and formulas for
#' the observation model from a hidden Markov model.
Observation <- R6Class(
  classname = "Observation",
  
  public = list(
    
    # Constructor -------------------------------------------------------------
    #' @description Create new Observation object
    #' 
    #' @param data Data frame containing response variables and covariates
    #' @param formulas List of formulas for observation parameters
    #' @param n_states Number of states (needed to construct model formulas)
    #' @param dists Named list of Distribution objects for each data stream
    #' @param par List of observation parameters (for covariate-free model)
    #' 
    #' @return A new Observation object
    initialize = function(data, 
                          dists, 
                          formulas = NULL, 
                          n_states, 
                          par = NULL) {
      private$check_args(data = data, 
                         dists = dists, 
                         n_states = n_states, 
                         par = par, 
                         formulas = formulas)
      
      # Make sure there is an ID column in the data and it's a factor
      if(is.null(data$ID)) {
        data$ID <- factor(1)
      } else {
        data$ID <- factor(data$ID)
      }
      
      # Set data and nstates attributes
      private$data_ <- data
      private$nstates_ <- n_states
      private$inipar_ <- par
      
      # Check for observed states 
      if ("state" %in% colnames(data)) {
        kn <- lapply(strsplit(as.character(data$state), ","), FUN = as.numeric)
        private$known_states_data_ <- data$state
        known_states <- matrix(1, nr = nrow(data), nc = n_states)
        for (i in 1:nrow(data)) {
          known_states[i, -kn[[i]]] <- 0
        }
        private$known_states_ <- known_states
        wh <- which(colnames(data) == "state")
        private$data_ <- data[,-wh]
      } else {
        private$known_states_ <- matrix(NA, nrow = nrow(data), nc = n_states)
      }
      
      # Define observation distributions
      if(all(sapply(dists, is.character))) {
        # If distributions passed as strings (i.e., names), get corresponding
        # Dist objects
        private$dists_ <- lapply(dists, function(name) private$dist_maker(name))
      } else {
        private$dists_ <- dists        
      }
      
      # Set formulas
      var_names <- colnames(self$obs_var())
      par_names <- lapply(self$dists(), FUN = function(x) {x$parnames()})
      private$formulas_ <- make_formulas(formulas, 
                                         var_names = var_names,
                                         par_names = par_names, 
                                         n_states = n_states)
      if (is.null(formulas)) {
        # Default: set all formulas to ~1
        forms <- lapply(par, function(varpar) {
          f <- lapply(varpar, function(...) {
            return(~1)            
          })
          return(f)
        })
        private$raw_formulas_ <- forms
      } else {
        private$raw_formulas_ <- formulas 
      }
      
      # Get names of all covariates
      var_names <- unique(rapply(self$formulas(), all.vars))
      # Remove pi from list of covariates if it is in the formulas
      var_names <- var_names[which(var_names!="pi")]
      if(length(var_names) > 0) {
        # Remove NAs in covariates (replace by last non-NA value)
        data[,var_names] <- lapply(data[,var_names, drop=FALSE], 
                                   function(col) na_fill(col))
        self$update_data(data)
      }
      
      # Save terms of model formulas
      mats <- self$make_mat()
      ncol_fe <- mats$ncol_fe
      ncol_re <- mats$ncol_re
      private$terms_ <- c(mats, list(names_fe = colnames(mats$X_fe),
                                     names_re_all = colnames(mats$X_re),
                                     names_re = colnames(ncol_re)))
      
      # Initialise parameters      
      self$update_coeff_fe(rep(0, sum(ncol_fe)))
      self$update_coeff_re(rep(0, ncol(mats$X_re)))
      self$update_lambda(rep(1, ifelse(is.null(ncol_re), 0, ncol(ncol_re))))
      
      # Fixed effect parameters     
      if(!is.null(par)) {
        # make sure par is in right order 
        n_var <- ncol(self$obs_var())
        varnm <- colnames(self$obs_var())
        corrected_par <- vector(mode = "list", length = n_var)
        for (i in 1:n_var) {
          subvars <- self$dists()[[i]]$parnames()
          if (!all(subvars %in% names(par[[varnm[i]]]))) {
            stop("parameters for variable ", varnm[i], " are missing or have wrong name")
          }
          subpar <- vector(mode = "list", length = length(subvars))
          for (j in 1:length(subvars)) {
            subpar[[j]] <- par[[varnm[i]]][[subvars[j]]]
          }
          names(subpar) <- subvars
          corrected_par[[i]] <- subpar
        }
        names(corrected_par) <- varnm
        self$update_par(corrected_par)
      } else {
        stop("'par' must be provided")
      }
      
    },
    
    # Accessors ---------------------------------------------------------------
    
    #' @description Data frame
    data = function() {return(private$data_)},
    
    #' @description List of distributions
    dists = function() {return(private$dists_)},
    
    #' @description Number of states
    nstates = function() {return(private$nstates_)},
    
    #' @description Parameters on natural scale
    #' 
    #' @param t time point, default t = 1; 
    #' if "all" then return for all time points, otherwise return at time points given in t 
    #' @param full_names Logical. If TRUE, the rows of the output
    #' are named in the format "variable.parameter" (default). If
    #' FALSE, the rows are names in the format "parameter". The
    #' latter is used in various internal functions, when the parameters
    #' need to be passed on to an R function.
    #' @param linpred custom linear predictor 
    #' @param as_list Logical. If TRUE, the output is a nested list with three levels:
    #' (1) time step, (2) observed variable, (3) observation parameter. If FALSE (default),
    #' the output is an array with one row for each observation parameter, one column for
    #' each state, and one slice for each time step.
    #' 
    #' @return Array of parameters with one row for each observation parameter, 
    #' one column for each state, and one slice for each time step. (See as_list
    #' argument for alternative output format.)
    par = function(t = 1, full_names = TRUE, linpred = NULL, as_list = FALSE) {
      # Number of states
      n_states <- self$nstates()
      
      # Number of parameters on natural scale (in each state)
      n_par <- sum(sapply(self$dists(), function(d) d$npar()))
      
      # Get linear predictor
      if (is.null(linpred)) linpred <- self$linpred() 
      
      # Number of observations
      n <- length(linpred) / (n_par * n_states)
      
      # Subset by time
      if (length(t) == 1) if (t == "all") t <- 1:n
      ind <- as.vector(sapply(1:(n_states * n_par), function(i) {t + (i - 1) * n}))
      linpred <- linpred[ind]
      
      # Matrix of linear predictor
      lp_mat <- matrix(linpred, ncol = n_par * n_states)
      
      if(as_list) {
        # List with three levels: (1) time step, (2) observed variable, 
        # (3) observation parameter
        par <- apply(lp_mat, 1, self$w2n)
      } else {
        # Matrix of natural parameters
        par_mat <- apply(lp_mat, 1, function(lp_vec) {
          par_ls <- self$w2n(lp_vec)
          par_vec <- unlist(par_ls, use.names = FALSE)
          return(par_vec)
        })
        
        # Array of natural parameters
        par_array <- array(par_mat, dim = c(n_states, n_par, length(t)))
        
        # Get parameter names 
        par_names <- unlist(lapply(self$dists(), FUN = function(x) {x$parnames()}), use.names = FALSE)
        if(full_names) {
          var_names <- unlist(sapply(1:length(self$dists()), FUN = function(d) 
          {rep(names(self$dists())[d], self$dists()[[d]]$npar())}), use.names = FALSE)
          par_names <- paste0(var_names, ".", par_names)
        } 
        names(par_names) <- NULL 
        
        # Set dimension names for rows and columns
        dimnames(par_array) <- list(paste("state", 1:n_states),
                                    par_names,
                                    NULL)
        
        # Transpose each slice
        par <- aperm(par_array, perm = c(2, 1, 3))
      }
      
      return(par)
    },
    
    #' @description Return initial parameter values supplied 
    inipar = function() {return(private$inipar_)}, 
    
    #' @description Fixed effect parameters on working scale
    coeff_fe = function() {return(private$coeff_fe_)},
    
    #' @description Random effect parameters
    coeff_re = function() {return(private$coeff_re_)},
    
    #' @description Fixed effect design matrix 
    X_fe = function() {return(private$terms_$X_fe)}, 
    
    #' @description Random effect design matrix 
    X_re = function() {return(private$terms_$X_re)}, 
    
    #' @description Smoothness parameters
    lambda = function() {return(private$lambda_)},
    
    #' @description Variance components of smooth terms
    #' 
    #' This function transforms the smoothness parameter of
    #' each smooth term into a standard deviation, given by 
    #' SD = 1/sqrt(lambda). It is particularly helpful to get the
    #' standard deviations of independent normal random effects.
    vcomp = function() {return(1/sqrt(private$lambda_))},
    
    #' @description List of model formulas for observation model
    #' 
    #' @param raw Logical. If FALSE, returns the nested list created by
    #' make_formulas (default). If TRUE, returns formulas passed as input.
    formulas = function(raw = FALSE) {
      if (raw) {
        return(private$raw_formulas_)   
      } else {
        return(private$formulas_)
      }
    },
    
    #' @description Terms of model formulas
    terms = function() {return(private$terms_)},
    
    #' @description  Data frame of response variables
    #' @param expand if TRUE then multivariate variables in observations are expanded 
    #' to be univariate, creating extra columns 
    obs_var = function(expand = FALSE) {
      obs_names <- names(self$dists())
      obs_var <- self$data()[, obs_names, drop = FALSE]
      datadim <- rep(1, ncol(obs_var))
      if (expand) {
        multivar <- sapply(obs_var, is.list)
        if (any(multivar)) {
          wh <- which(multivar)
          for (i in 1:length(wh)) {
            v <- do.call(rbind, obs_var[[wh[i]]])
            datadim[wh[i]] <- ncol(v)
            tmp <- NULL
            tmpnms <- NULL
            if (wh > 1) {
              tmp <- cbind(tmp, obs_var[,1:(wh[i] - 1)])
              tmpnms <- c(tmpnms, colnames(obs_var)[1:(wh[i] - 1)])
            }
            tmp <- cbind(tmp, v)
            tmpnms <- c(tmpnms, rep(names(wh[i]), ncol(v)))
            if (wh < ncol(obs_var)) {
              tmp <- cbind(tmp, obs_var[,(wh[i] + 1):ncol(obs_var)])
              tmpnms <- c(tmpnms, colnames(obs_var)[(wh[i] + 1):ncol(obs_var)]) 
            }
            obs_var <- tmp
            colnames(obs_var) <- tmpnms
          }
        }
      }
      attributes(obs_var)$datadim <- datadim
      return(obs_var)
    },
    
    #' @description Vector of known states 
    known_states = function(mat = TRUE) {
      if (mat) {
        return(private$known_states_)
      } else {
        return(private$known_states_data_)
      }
    }, 
    
    # Mutators ----------------------------------------------------------------
    
    #' @description Update parameters
    #' 
    #' Updates the 'par' attribute to the list passed as input,
    #' and updates the intercept elements of 'coeff_fe' using
    #' the list passed as input
    #' 
    #' @param par New list of parameters
    update_par = function(par) {
      
      # Get index of first column of X_fe for each parameter
      ncol_fe <- self$terms()$ncol_fe
      n_par <- length(ncol_fe)
      i0 <- c(1, cumsum(ncol_fe)[-n_par] + 1)
      
      # Apply link to get parameters on working scale
      private$coeff_fe_[i0] <- self$n2w(par)
    },
    
    #' @description Update coefficients for fixed effect parameters
    #' 
    #' @param coeff_fe New vector of coefficients for fixed effect 
    #' parameters
    update_coeff_fe = function(coeff_fe) {
      private$coeff_fe_ <- matrix(coeff_fe)
      rownames(private$coeff_fe_) <- self$terms()$names_fe
    },
    
    #' @description Update random effect parameters
    #' 
    #' @param coeff_re New vector of coefficients for random effect 
    #' parameters
    update_coeff_re = function(coeff_re) {
      private$coeff_re_ <- matrix(coeff_re)
      rownames(private$coeff_re_) <- self$terms()$names_re_all
    },
    
    #' @description Update fixed effect design matrix
    #' 
    #' @param X_fe New fixed effect design matrix 
    update_X_fe = function(X_fe) {
      private$terms_$X_fe <- X_fe
    }, 
    
    #' @description Update random effect design matrix
    #' 
    #' @param X_re New random effect design matrix 
    update_X_re = function(X_re) {
      private$terms_$X_re <- X_re
    }, 
    
    #' @description Update smoothness parameters
    #' 
    #' @param lambda New smoothness parameter vector
    update_lambda = function(lambda) {
      private$lambda_ <- matrix(lambda)
      rownames(private$lambda_) <- self$terms()$names_re
    },
    
    #' @description Update data
    #' 
    #' @param data New data frame
    update_data = function(data) {
      private$data_ <- data
    },
    
    
    # Other methods -----------------------------------------------------------
    
    #' @description Make model matrices
    #' 
    #' @param new_data Optional new data set, including covariates for which
    #' the design matrices should be created. If this argument is not specified,
    #' the design matrices are based on the original data frame. 
    #' 
    #' @return A list with elements:
    #' \itemize{
    #'   \item{X_fe}{Design matrix for fixed effects}
    #'   \item{X_re}{Design matrix for random effects}
    #'   \item{S}{Smoothness matrix for random effects}
    #'   \item{ncol_fe}{Number of columns of X_fe for each parameter}
    #'   \item{ncol_re}{Number of columns of X_re and S for each random effect}
    #' }
    make_mat = function(new_data = NULL) {
      make_matrices(formulas = self$formulas(),
                    data = self$data(),
                    new_data = new_data)
    },
    
    #' Design matrices for grid of covariates
    #' 
    #' @param var Name of variable
    #' @param covs Optional named list for values of covariates (other than 'var') 
    #' that should be used in the plot (or dataframe with single row). If this is
    #' not specified, the mean value is used for numeric variables, and the
    #' first level for factor variables.
    #' @param n_grid Grid size (number of points). Default: 1000.
    #' 
    #' @return A list with the same elements as the output of make_mat, 
    #' plus a data frame of covariates values.
    make_newdata_grid = function(var, covs = NULL, n_grid = 1e3) {
      # Data frame for covariate grid
      new_data <- cov_grid(var = var, data = self$data(), 
                           covs = covs, formulas = self$formulas(),
                           n_grid = n_grid)
      
      return(new_data)
    },
    
    #' @description Natural to working parameter transformation
    #' 
    #' This function applies the link functions of the distribution
    #' parameters, to transform parameters from their natural scale
    #' to the working scale (i.e., linear predictor scale)
    #' 
    #' @param par List of parameters on natural scale
    #' 
    #' @return Vector of parameters on working scale
    n2w = function(par) {
      wpar <- lapply(seq_along(self$dists()), 
                     function(i) self$dists()[[i]]$n2w(par[[i]]))
      names(wpar) <- names(par)
      wpar <- unlist(wpar)
      return(wpar)
    },
    
    #' @description  Working to natural parameter transformation
    #'
    #' This function applies the inverse link functions of the
    #' distribution parameters, to transform parameters from the working
    #' scale (i.e., linear predictor scale) to their natural scale.
    #'
    #' @param wpar Vector of parameters on working scale
    #' 
    #' @return List of parameters on natural scale
    w2n = function(wpar) {
      # Initialise list of natural parameters
      par <- list()
      
      # Number of observed variables
      nvar <- length(self$dists())
      # Number of states
      n_states <- self$nstates()
      
      # Counter to subset observation parameters
      par_count <- 1
      
      # Loop over observed variables
      for(var in 1:nvar) {
        # Number of parameters for this distribution
        npar <- self$dists()[[var]]$npar()
        # Subset and transform working parameters
        sub_wpar <- wpar[par_count:(par_count + npar*n_states - 1)]
        par_count <- par_count + npar*n_states
        par[[var]] <- self$dists()[[var]]$w2n(sub_wpar)
      }
      
      names(par) <- names(self$dists())
      return(par)
    },
    
    #' @description Compute linear predictor 
    linpred = function() {
      linpred <- self$X_fe() %*% self$coeff_fe() + self$X_re() %*% self$coeff_re()
      return(linpred[,1])
    }, 
    
    #' @description Observation likelihoods
    #' 
    #' @param X_fe Design matrix for fixed effects
    #' @param X_re Design matrix for random effects
    #' @param data optional dataframe to include in form of obs_var() output 
    #' 
    #' @return Matrix of likelihoods of observations, with one row for each 
    #' time step, and one column for each state.
    obs_probs = function(data = NULL) {
      # Data frame of observations
      if (is.null(data)) {
        data <- self$obs_var()
        X_fe_old <- NULL
      } else {
        X_fe_old <- self$X_fe() 
        X_re_old <- self$X_re() 
        mats <- self$make_mat(data)
        self$update_X_fe(mats$X_fe)
        self$update_X_re(mats$X_re)
      }
      
      # Number of observations
      n <- nrow(data)
      # Number of states
      n_states <- self$nstates()
      # Number of variables
      n_var <- ncol(self$obs_var())
      
      # State-dependent parameters
      par <- self$par(t = "all", full_names = FALSE)
      
      # Initialise matrix of probabilities
      prob <- matrix(1, nrow = n, ncol = n_states)
      for(i in which(!is.na(self$known_states(mat = FALSE)))) {
        # Set other probabilities to zero if state is known
        prob[i,self$known_states()[i,] == 0] <- 0
      }
      
      # Counter to subset parameter vector
      par_count <- 1
      
      # Get variable names 
      givenvarnms <- colnames(data)
      varnms <- names(self$obs_var())
      
      # Loop over observed variables
      for(var in 1:n_var) {
        obsdist <- self$dists()[[var]]
        par_ind <- par_count:(par_count + obsdist$npar() - 1)
        
        if (varnms[var] %in% givenvarnms) {
          # Loop over observations (rows of prob)
          for (i in 1:n) {
            # Don't update likelihood is observation is missing
            if(!is.na(data[i, varnms[var]])) {
              # Loop over states (columns of prob)
              for (s in 1:n_states) {
                prob[i, s] <- prob[i, s] * 
                  obsdist$pdf_apply(x = data[i, varnms[var]], par = par[par_ind, s, i])
              }            
            }
          }
        }
        par_count <- par_count + obsdist$npar()
      }
      
      # reset design matrices
      if (!is.null(X_fe_old)) {
        self$update_X_fe(X_fe_old)
        self$update_X_re(X_re_old)
      }
      
      return(prob)
    },
    
    #' @description Plot histogram of data and pdfs
    #' 
    #' Plot histogram of observations for the variable specified by the argument name, 
    #' overlaid with the pdf of the specified distribution for that data stream. 
    #' Helpful to select initial parameter values for model fitting, or to visualise 
    #' fitted state-dependent distributions.
    #' 
    #' @param name Name of response variable for which the histogram
    #' and pdfs should be plotted.
    #' @param weights Optional vector of length the number of pdfs that are
    #' plotted. Useful to visualise a mixture of distributions weighted by the
    #' proportion of time spent in the different states.
    #' @param t Index of time step to use for covariates (default: 1).
    #' 
    #' @return A ggplot object
    plot_dist = function(name, weights = NULL, t = 1) {
      # Extract observed values for relevant variable
      obs <- data.frame(val = self$data()[[name]])
      
      # Matrix of parameters
      par <- matrix(unlist(self$par(t = t, as_list = TRUE)[[1]][[name]]), 
                    nrow = self$nstates())
      colnames(par) <- self$dists()[[name]]$parnames()
      
      # Weights for each state-dependent distribution
      if(is.null(weights)) weights <- rep(1, self$nstates())
      
      # Grid over range of observed variable
      n_grid <- 1e3
      grid <- seq(min(obs, na.rm = TRUE), max(obs, na.rm = TRUE), length = n_grid)
      
      # Check if variable is integer 
      if (is_whole_number(obs)) {
        grid <- unique(floor(grid))
        n_grid <- length(grid)
      }
      
      # Loop over states
      vals <- matrix(NA, nrow = n_grid, ncol = self$nstates() + 1)
      for(state in 1:self$nstates()) {
        # Define list of arguments to pass to pdf
        args <- list(grid)
        args <- c(args, par[state,])
        
        # Compute state-dependent pdf
        vals[,state] <- weights[state] * do.call(self$dists()[[name]]$pdf(), args)
      }
      # Weighted sum of state-dependent pdfs
      vals[, self$nstates() + 1] <- rowSums(vals[, 1:self$nstates()])
      
      # Data frame of state-dependent densities
      df_dens <- data.frame(
        state = rep(c(paste("State", 1:self$nstates()), "Total"), each = n_grid),
        grid = grid,
        val = as.vector(vals))
      
      # Create hist object to extract ylim
      breaks <- seq(min(obs, na.rm = TRUE), max(obs, na.rm = TRUE), length = 20)
      h <- hist(obs$val, breaks = breaks, plot = FALSE)
      
      # Create ggplot histogram
      p <- ggplot(obs, aes(x = val)) + xlab(name) +
        geom_histogram(breaks = breaks, aes(y=..density..), 
                       col = "white", bg = "lightgrey", na.rm = TRUE) + 
        geom_line(aes(grid, val, col = state, linetype = state), 
                  data = df_dens, size = 0.7) +
        scale_color_manual("", values = c(hmmTMB_cols[1:self$nstates()], "black")) +
        scale_linetype_manual("", values = c(rep(1, self$nstates()), 2)) +
        coord_cartesian(ylim = c(0, 1.1 * max(h$density))) +
        theme_light()
      
      return(p)
    },
    
    #' @description Print model formulation
    formulation = function() {
      cat("#######################\n")
      cat("## Observation model ##\n")
      cat("#######################\n")
      # List of distribution names
      d_list <- lapply(self$dists(), function(d) d$name())
      # List of parameter names
      p_list <- lapply(self$dists(), function(d) names(d$link()))
      # List of fixed parameters
      fix_list <- lapply(self$dists(), function(d) d$fixed())
      
      # List of parameter formulas
      s_list <- lapply(self$formulas(), function(f) {
        ff <- unlist(f)
        s <- NULL
        for(i in seq_along(ff))
          s <- paste0(s, "  * ", names(ff)[i], " ~ ", as.character.default(ff[[i]])[2], "\n")
        return(s)
      })
      
      # Remove formulas for fixed parameters
      nms <- names(s_list)
      for (i in 1:length(fix_list)) {
        wh <- which(nms == names(fix_list)[i])
        fix <- names(fix_list[[i]])[fix_list[[i]] == TRUE]
        if (length(fix) == 0) next 
        for (j in 1:length(fix)) {
          del <- paste0(".*?", fix[j], ".*?\n")
          s_list[[wh]] <- gsub(del, "", s_list[[wh]])
        }
      }
      
      # Loop over observed variables
      for(i in seq_along(d_list)) {
        # Print variable distribution
        cat(paste0("+ ", names(d_list)[i], " ~ ", d_list[[i]], "(", 
                   paste0(p_list[[i]], collapse = ", "), ")"), "\n")
        # Print parameter formulas
        cat(s_list[[i]], "\n")
      }
    },
    
    #' @description Print Observation object
    print = function() {
      self$formulation()
    }
    
  ),
  
  private = list(
    
    # Private data members ----------------------------------------------------
    
    data_ = NULL,
    known_states_ = NULL,
    known_states_data_ = NULL, 
    dists_ = NULL,
    nstates_ = NULL,
    coeff_fe_ = NULL,
    coeff_re_ = NULL,
    lambda_ = NULL,
    formulas_ = NULL,
    raw_formulas_ = NULL, 
    inipar_ = NULL, 
    terms_ = NULL,
    mats_ = NULL, 
    
    #' Check constructor arguments 
    # (For argument description, see constructor)
    check_args = function(data, dists, n_states, par, formulas) {
      if(!inherits(data, "data.frame")) {
        stop("'data' should be a data.frame")
      }
      
      # Check that time intervals are regular if 'time' is provided
      if(!is.null(data$time)) {
        # Get indices of start and end of time series
        if(!is.null(data$ID)) {
          i0 <- which(data$ID[-1] != data$ID[-nrow(data)])
          start <- c(1, i0 + 1)
          end <- c(i0, nrow(data))
        } else {
          start <- 1
          end <- nrow(data)
        }
        
        # Time intervals between data rows
        dt <- data$time[-start] - data$time[-end]
        
        # Length of range of time intervals
        dt_range <- as.numeric(diff(range(dt, na.rm = TRUE)))
        # Median time interval
        dt_median <- as.numeric(median(dt, na.rm = TRUE))
        
        # If (max - min) is longer than 0.5*median, send warning (arbitrary threshold)
        if(dt_range/dt_median > 0.5) {
          warning(paste("'data$time' seems to be irregular. Data rows should be at",
                        "regular time intervals."))
        }
      }
      
      if(!is.list(dists)) {
        stop("'dists' should be a list")
      }
      
      if(!all(sapply(dists, inherits, "Dist")) & !all(sapply(dists, is.character))) {
        stop(paste("Elements of 'dists' should all be either character strings",
                   "(i.e., distribution names), or Dist objects"))
      }
      
      if(!all(names(dists) %in% colnames(data))) {
        stop("Variable name in 'dists' not found in data")
      }
      
      if(!is.numeric(n_states) | n_states < 1) {
        stop("'n_states' should be a numeric >= 1")
      }
      
      if(!is.null(par)) {
        if(!is.list(par) | length(par) != length(dists)) {
          stop("'par' should be a list of same length as 'dists'")
        }
        
        if(!all(rapply(par, length) == n_states) | !all(rapply(par, is.numeric))) {
          stop("Elements of 'par' should be numeric vectors of length 'n_states'")
        }
        
        if(!all(names(par) == names(dists))) {
          stop("'par' should have the same names as 'dists'")
        }
      }
      
      if(!is.null(formulas)) {
        if(!is.list(formulas) |
           !all(rapply(formulas, function(x) inherits(x, "formula")))) {
          stop("'formulas' should be a list of R formulas")
        }
        
        if(!all(names(formulas) %in% names(dists))) {
          stop("'formulas' should have the same names as 'dists'")
        }
      }
    }, 
    
    # Create a distribution 
    # @param name name of distribution to create 
    dist_maker = function(name) {
      if (name %in% names(dist_list)) {
        # distribution with fixed parameter dimension
        return(dist_list[[name]]$clone())
      } else {
        # distribution with a variable dimension
        subname <- gsub("[0-9]+", "", name)
        if (!(subname %in% names(dist_list))) stop("distribution unknown")
        tmp <- dist_list[[subname]]$clone()
        getdim <-as.numeric(gsub("[^0-9.]", "", name))
        if (subname == "cat") {
          tmp$set_npar(getdim - 1)
          tmp$set_parnames(paste0("p", 1:(getdim - 1)))
        } else if (subname == "mvnorm") {
          tmp$set_npar(2 * getdim + (getdim^2 - getdim) / 2)
          V <- matrix(1:getdim, nr = getdim, nc = getdim)
          tV <- t(V)
          tmp$set_parnames(c(paste0("mu", 1:getdim), 
                             paste0("sd", 1:getdim), 
                             paste0("corr", V[upper.tri(V)], tV[upper.tri(tV)])))
        } else if (subname == "dir") {
          tmp$set_npar(getdim)
          tmp$set_parnames(paste0("alpha", 1:getdim))
        }
        return(tmp)
      }
    }
  )
)
