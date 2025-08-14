#' @title R6 Class for Forecast in Hidden Markov Models
#'
#' @description
#' Creates a forecast object from a fitted HMM model, projecting future state probabilities and unconditional observation densities for a user-defined horizon using a covariate data-frame if applicable. Forecasting occurs in the constructor; results are accessible via public methods.
#' 
#' @export
Forecast <- R6Class(
  classname = "Forecast",
  
  public = list(
    
    # Constructor -------------------------------------------------------------
    
    #' @description Construct a forecast object.
    #' 
    #' @details
    #' Algorithm overview
    #' \enumerate{
    #'   \item Validate inputs (private$validate_params()).
    #'   \item Build or accept future design matrix (forecast_data).
    #'   \item Obtain time-varying transition matrices (tpm_forecast) and
    #'         observation parameters (obs_par_forecast) via hmm$predict().
    #'   \item Forward-propagate hidden-state distribution (hidden_state_forecast).
    #'   \item Marginalise over states for unconditional pdfs on user-defined grid (forecast_dists).
    #' }
    #' 
    #' @param hmm A fitted HMM object.
    #' @param n Integer. Forecast horizon (ignored if forecast_data provided).
    #' @param forecast_data data.frame of future covariates / IDs.
    #'        Must include all covariates in fitted hmm.
    #' @param preset_eval_range Named list: elements are numeric vectors (univariate) or matrices
    #'        (multivariate, rows as dimension ex. mu1, mu2, columns as points).
    #'        Missing entries default to 90–110% range of training data.
    #' @param starting_state_distribution Numeric vector of length nstates(), or
    #'        "last" (Fitted model at final training row, propagated one step) or
    #'        "stationary" (model’s stationary distribution).
    #' 
    #' @return An R6 object with notable public accessor methods (invoke with ()):
    #' \describe{
    #'   \item{hidden_state_forecast}{nstates × n forecast steps.}
    #'   \item{obs_par_forecast}{Array of time-varying observation parameters.}
    #'   \item{tpm_forecast}{Array of time-varying transition matrices.}
    #'   \item{forecast_dists}{List of unconditional pdf matrices, one per observation variable.}
    #'   \item{eval_range}{List of grids on which each pdf was evaluated.}
    #' }
    #' 
    #' @examples
    #' \dontrun{
    #' data(faithful)
    #' 
    #' # Create training and forecast data frames
    #' training_fraction <- 0.8
    #' n_training <- floor(nrow(faithful) * training_fraction)
    #' training_df <- data.frame(
    #'   ID = rep(1, n_training),
    #'   waiting = faithful$waiting[1:n_training],
    #'   eruptions = faithful$eruptions[1:n_training]
    #' )
    #' forecast_df <- data.frame(
    #'   ID = rep(1, nrow(faithful) - n_training),
    #'   waiting = faithful$waiting[(n_training + 1):nrow(faithful)],
    #'   eruptions = faithful$eruptions[(n_training + 1):nrow(faithful)]
    #' )
    #' # Create 2-state model with non-linear effect of waiting on all transitions
    #' hid_model <- MarkovChain$new(
    #'   data = training_df,
    #'   n_states = 2,
    #'   formula = ~ s(waiting, k = 10, bs = "cs")
    #' )
    #' 
    #' # Create observation model with normal distribution for duration
    #' obs_model <- Observation$new(
    #'   data = training_df,
    #'   n_states = 2,
    #'   dists = list(eruptions = "norm"),
    #'   par = list(eruptions = list(
    #'     mean = c(0, 0),
    #'     sd = c(1, 1)
    #'   ))
    #' )
    #' # Update model parameters to suggested
    #' obs_model$update_par(par = obs_model$suggest_initial())
    #' # Create HMM model
    #' hmm <- HMM$new(
    #'   hid = hid_model,
    #'   obs = obs_model
    #' )
    #' hmm$fit(silent = TRUE)
    #' # Create forecast object
    #' forecast <- Forecast$new(
    #'   hmm = hmm,
    #'   forecast_data = forecast_df,
    #'   starting_state_distribution = "last"
    #'   )
    #' }
    initialize = function(hmm               = NULL,
                          n                 = NULL,
                          forecast_data     = NULL,
                          preset_eval_range     = NULL,
                          starting_state_distribution = "last") {

      ## -- 1  Input checking --------------------------------------------------
      private$validate_params(
        hmm                       = hmm,
        n                         = n,
        forecast_data             = forecast_data,
        preset_eval_range             = preset_eval_range,
        starting_state_distribution = starting_state_distribution
      )

      ## -- 2  Store the fitted model & basic metadata -------------------------
      private$hmm_             <- hmm
      private$observation_vars_ <- colnames(private$hmm_$obs()$obs_var())
      private$starting_state_distribution_ <- starting_state_distribution

      ## -- 3  Build or accept the future design matrix ------------------------
      private$forecast_data_ <- if (!is.null(forecast_data)) {
        forecast_data
      } else {
        as.data.frame(
          c(
            list(ID = rep(1, n)),
            stats::setNames(
              replicate(length(private$observation_vars_), rep(NA, n),
                        simplify = FALSE),
              private$observation_vars_
            )
          )
        )
      }

      ## -- 4  Predict forward-looking parameters ------------------------------
      private$obs_par_forecast_ <- private$hmm_$predict("obspar",
                                                newdata = private$forecast_data_)
      private$tpm_forecast_     <- private$hmm_$predict("tpm",
                                                newdata = private$forecast_data_)

      ## -- 5  Establish evaluation grids for each response variable -----------
      private$eval_range_ <- private$configure_eval_range(
        eval_range   = preset_eval_range,
        obs_vars = private$observation_vars_,
        data     = private$hmm_$obs()$data()
      )

      ## -- 6  Choose the initial hidden-state distribution --------------------
      if (is.character(starting_state_distribution)) {

        if (starting_state_distribution == "last") {
          # one-step-ahead distribution conditional on final training point
          last_sp  <- utils::tail(private$hmm_$state_probs(), 1)
          last_tpm <- private$hmm_$hid()$tpm(nrow(private$hmm_$obs()$data()))[, , 1]
          dist_0   <- last_sp %*% last_tpm
        } else if (starting_state_distribution == "stationary") {
          dist_0 <- private$hmm_$hid()$stationary()
        } else {
          stop("Character `starting_state_distribution` must be ",
               "'last' or 'stationary'", call. = FALSE)
        }
      } else {
        dist_0 <- starting_state_distribution
      }

      ## -- 7  Forward-propagate hidden states ---------------------------------
      n_steps <- nrow(private$forecast_data_)
      hidden_state_forecast <-
        array(NA_real_,
              dim = c(private$hmm_$hid()$nstates(), n_steps))

      hidden_state_forecast[, 1] <- dist_0

      if (n_steps > 1) {
        for (t in 2:n_steps) {
          hidden_state_forecast[, t] <-
            hidden_state_forecast[, t - 1] %*% private$tpm_forecast_[, , t - 1]
        }
      }
      private$hidden_state_forecast_ <- hidden_state_forecast
      
      ## -- 8  Build unconditional predictive pdfs -----------------------------
      private$forecast_dists_ <- vector("list", length(private$observation_vars_))
      names(private$forecast_dists_) <- private$observation_vars_

      for (obs_var in private$observation_vars_) {

        # Get distribution and parameters for the current observation variable
        obs_dists      <- private$hmm_$obs()$dists()[[obs_var]]

        model_params   <- names(private$obs_par_forecast_[, 1, 1])
        if (is.null(model_params)) {
          # Edge case: when obs_par_forecast has only 1 parameter it is unnamed
          # Use current_params = 1 to unpack the matrix correctly
          current_params <- 1
        } else {
          # Otherwise, obtain the parameters for the current dimension
          current_params <- grep(
            paste0("^", obs_var, "(\\.)"), model_params, value = TRUE
          )
        }

        obs_eval_range <- private$eval_range_[[obs_var]]

        # In the case where distribution is multivariate, (dirichlet, mvnorm) we
        # need a list where each element is a vector of x-values.
        if (is.matrix(obs_eval_range)) {
          multi_variate <- TRUE
          obs_eval_range <- as.list(
            as.data.frame(obs_eval_range, stringsAsFactors = FALSE)
          )
        } else {
          multi_variate <- FALSE
        }

        private$forecast_dists_[[obs_var]] <-
          array(NA_real_, dim = c(length(obs_eval_range), n_steps))

        for (i in seq_len(n_steps)) {

          # Evaluate each state's pdf at eval_range for time step i
          # pdf_matrix: matrix with dimensions |eval_range| x n_states
          pdf_matrix <- vapply(
            seq_len(private$hmm_$hid()$nstates()),   # loop over hidden states
            function(s) {                       # compute pdf for state s
              obs_dists$pdf_apply(
                x   = obs_eval_range,
                par = stats::setNames(
                  private$obs_par_forecast_[current_params, s, i],
                  obs_dists$parnames()
                )
              )
            },
            # vapply template: numeric vector of length |eval_range|
            numeric(length(obs_eval_range))
          )

          # Compute forecasted pdf
          private$forecast_dists_[[obs_var]][, i] <- pdf_matrix %*% private$hidden_state_forecast_[, i]

        }
      }
    },

    # Accessors ---------------------------------------------------------------
    #' @description Get predicted observation parameters.
    obs_par_forecast = function() {
      return(private$obs_par_forecast_)
    },

    #' @description Get predicted transition matrices.
    tpm_forecast = function() {
      return(private$tpm_forecast_)
    },

    #' @description Get the evaluation grid for each response variable.
    eval_range = function() {
      return(private$eval_range_)
    },

    #' @description Get the forecast data used for predictions.
    forecast_data = function() {
      return(private$forecast_data_)
    },

    #' @description Update the evaluation grid for forecast pdfs.
    #' This method reconfigures the evaluation range used to compute unconditional predictive pdfs.
    #' It validates the new eval_range, updates internal state, reinitializes forecast densities,
    #' and returns the object invisibly for chaining.
    #'
    #' @param eval_range A named list specifying the new evaluation grid for each observation variable.
    #' @return The Forecast object (invisibly), allowing method chaining.
    update_eval_range = function(eval_range) {
      
      # Validate the new evaluation range along with forecast data and model
      private$validate_params(
        hmm               = private$hmm_,
        forecast_data     = private$forecast_data_,
        preset_eval_range = eval_range
      )

      # Configure and update the internal evaluation grid using new input
      private$eval_range_ <- private$configure_eval_range(
        eval_range = eval_range,
        obs_vars   = private$observation_vars_,
        data       = private$hmm_$obs()$data()
      )

      # Recompute forecast distributions with the updated evaluation grid.
      # Note: This call to initialize() recalculates dependent quantities such as forecast_dists.
      self$initialize(
        hmm                         = private$hmm_,
        forecast_data               = private$forecast_data_,
        preset_eval_range           = eval_range,
        starting_state_distribution = private$hidden_state_forecast_[, 1]
      )

      # Return the forecast object invisibly for potential method chaining.
      invisible(self)
    },

    #' @description Get the forward state probabilities.
    hidden_state_forecast = function() {
      return(private$hidden_state_forecast_)
    },

    #' @description Get the unconditional predictive pdfs.
    forecast_dists = function() {
      return(private$forecast_dists_)
    }
  ),

  ## ---------------------------------------------------------------------------
  ## Private helpers -----------------------------------------------------------
  ## ---------------------------------------------------------------------------
  private = list(

    # Private data members
    hmm_                   = NULL,  # fitted HMM model
    observation_vars_      = NULL,  # names of response variables
    forecast_data_         = NULL,  # future covariate data-frame
    obs_par_forecast_      = NULL,  # predicted observation parameters
    tpm_forecast_          = NULL,  # predicted transition matrices
    eval_range_            = NULL,  # grid for each response variable
    hidden_state_forecast_ = NULL,  # forward state probabilities
    forecast_dists_        = NULL,  # unconditional predictive pdfs
    starting_state_distribution_ = NULL,  # initial state distribution

    ## -- 1  Comprehensive argument checks -------------------------------------
    validate_params = function(hmm = NULL, n = NULL, forecast_data = NULL,
                               preset_eval_range = NULL,
                               starting_state_distribution = NULL) {

      if (is.null(hmm))
        stop("`hmm` must be provided", call. = FALSE)
      if (!inherits(hmm, "HMM"))
        stop("`hmm` must inherit from class 'HMM'", call. = FALSE)

      dist_names <- lapply(hmm$obs()$dists(), function(x) x$name())
      if (any(unlist(dist_names) %in% c("cat", "dirc", "tweedie"))) {
        stop("Forecasting not currently supported for 'cat', 'dirc', or 'tweedie' distributions",
             call. = FALSE)
      }

      # n  vs  forecast_data ---------------------------------------------------
      if (is.null(n) && is.null(forecast_data))
        stop("Supply either `n` or `forecast_data`", call. = FALSE)

      if (!is.null(n) && (!is.numeric(n) || n <= 0))
        stop("`n` must be a positive integer", call. = FALSE)

      if (!is.null(forecast_data) && !is.data.frame(forecast_data))
        stop("`forecast_data` must be a data.frame", call. = FALSE)

      # Covariate coverage -----------------------------------------------------
      covariates <- unique(Filter(
        function(x) !is.null(x) && x != "",     # remove NULLs and empty strings
        c(
          sapply(hmm$hid()$formulas(), all.vars),
          rapply(hmm$obs()$formulas(), all.vars)
        )
      ))

      if (length(covariates)) {
        if (is.null(forecast_data))
          stop("Provide `forecast_data` when covariates are in the model",
               call. = FALSE)

        missing_covs <- setdiff(covariates, colnames(forecast_data))
        if (length(missing_covs))
          stop("`forecast_data` is missing covariates: ",
               paste(missing_covs, collapse = ", "), call. = FALSE)
      }

      # x-grid list ------------------------------------------------------------
      if (!is.null(preset_eval_range) &&
          (!is.list(preset_eval_range) ||
             is.null(names(preset_eval_range)) ||
             any(names(preset_eval_range) == ""))) {
        stop("`preset_eval_range` must be a *named* list", call. = FALSE)
      }

      # starting_state_distribution -------------------------------------------
      if (is.null(starting_state_distribution))
        stop("`starting_state_distribution` must be provided", call. = FALSE)

      if (is.character(starting_state_distribution) &&
          !starting_state_distribution %in% c("last", "stationary"))
        stop("Character `starting_state_distribution` must be ",
             "'last' or 'stationary'", call. = FALSE)

      if (!is.character(starting_state_distribution) &&
          (!is.numeric(starting_state_distribution) ||
             length(starting_state_distribution) != hmm$hid()$nstates()))
        stop("`starting_state_distribution` must have length nstates()",
             call. = FALSE)
    },

    ## -- 2  Build default x-grids if needed -----------------------------------
    configure_eval_range = function(eval_range   = NULL,
                                obs_vars = NULL,
                                data     = NULL) {

      if (is.null(eval_range)) eval_range <- list()

      for (obs_var in obs_vars) {

        if (is.null(eval_range[[obs_var]])) {
          warning(sprintf(
            "Using 100 evaluation points in range ±10%% of training data for observation variable '%s'. Setting eval_range explicitly is recommended.",
            obs_var
          ), call. = FALSE)
        }

        if (is.null(eval_range[[obs_var]])) {
            max_range <- max(data[[obs_var]], na.rm = TRUE)
            min_range <- min(data[[obs_var]], na.rm = TRUE)
            max_range <- if (max_range > 0) max_range * 1.1 else max_range * 0.9
            min_range <- if (min_range > 0) min_range * 0.9 else min_range * 1.1
          eval_range[[obs_var]] <- seq(min_range, max_range, length.out = 100)
        }
      }
      eval_range
    }
  )
)
