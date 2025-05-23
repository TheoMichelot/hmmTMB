#' @title R6 class for forecasting HMM models
#'
#' @description
#' Forecast Object which requires passing a fitted HMM model.
#' 
#' @export
forecast <- R6Class(
  "forecast",

  # Public methods
  public = list(
    hmm = NULL,
    observation_vars = NULL,
    forecast_data = NULL,
    obs_par_forecast = NULL,
    tpm_forecast = NULL,
    x_vals = NULL,
    hidden_state_forecast = NULL,
    forecasted_pdfs = NULL,

    initialize = function(
      hmm = NULL,
      n = NULL,
      forecast_data = NULL,
      preset_x_vals = NULL,
      starting_state_distribution = "last"
    ) {
      # Check Arguments
      private$validate_params(
        hmm = hmm,
        n = n,
        forecast_data = forecast_data,
        preset_x_vals = preset_x_vals,
        starting_state_distribution = starting_state_distribution
      )

      # Set public parameters
      self$hmm <- hmm
      self$observation_vars <- colnames(self$hmm$obs()$obs_var())
      self$forecast_data <- if (!is.null(forecast_data)) {
        forecast_data
      } else {
        as.data.frame(
          c(
            list(ID = rep(1, n)),
            setNames(
              replicate(length(self$observation_vars), rep(NA, n), simplify = FALSE),
              self$observation_vars
            )
          )
        )
      }
      self$obs_par_forecast <- self$hmm$predict("obspar", newdata = self$forecast_data)
      self$tpm_forecast <- self$hmm$predict("tpm", newdata = self$forecast_data)
      self$x_vals <- private$configure_x_vals(
        x_vals = preset_x_vals,
        obs_vars = self$observation_vars,
        data = self$hmm$obs()$data()
      )

      # Handle starting_state_distribution
      if (is.character(starting_state_distribution)) {
        # if 'last' then use last state distribution of the training data
        if (starting_state_distribution == "last") {
          # Get starting distribution of hidden states
          last_state_distribution <- tail(self$hmm$state_probs(), 1)
          last_training_tpm <- self$hmm$hid()$tpm(nrow(self$hmm$obs()$data()))[, , 1]
          dist_0 <- last_state_distribution %*% last_training_tpm
        }
        # if 'stationary' then use stationary distribution of the model
        if (starting_state_distribution == "stationary") {
          dist_0 <- self$hmm$hid()$stationary()
        }
      } else {
        dist_0 <- starting_state_distribution
      }

      # Initialize the forecast matrix to store results
      hidden_state_forecast <- array(NA, dim = c(self$hmm$hid()$nstates(), nrow(self$forecast_data)))

      # Set the initial distribution
      hidden_state_forecast[, 1] <- dist_0

      # Loop through the remaining time steps
      for (t in 2:nrow(forecast_data)) {
        hidden_state_forecast[, t] <- hidden_state_forecast[, t - 1] %*% self$tpm_forecast[, , t-1]
      }

      self$hidden_state_forecast <- hidden_state_forecast

      # Loop through the observation variables
      self$forecasted_pdfs <- list()
      for (dimension in self$observation_vars) {

        # Step 1 - Get observation distribution function
        obs_dists <- hmm$obs()$dists()[[dimension]]

        # Step 1.5 - Get the parameters of the observation distribution
        pdf_params <- paste0(dimension, ".", names(formals(obs_dists$pdf())))
        model_params <- names(self$obs_par_forecast[ , 1, 1])
        current_params <- intersect(pdf_params, model_params)

        # Step 2 - Loop through the forecasted parameters and
        # calculate the weighted pdf for each hidden state
        x_vals <- self$x_vals[[dimension]]
        forecasted_pdfs[[dimension]] <- array(NA, dim = c(length(x_vals), nrow(self$forecast_data)))
        for (i in seq_len(nrow(self$forecast_data))) {
          pdf_matrix <- sapply(seq_len(self$hmm$hid()$nstates()), function(s) {
            obs_dists$pdf_apply(
              x = x_vals,
              par = setNames(
                self$obs_par_forecast[current_params, s, i], 
                obs_dists$parnames()
              )
            )
          })
          # Normalize the probabilities
          forecasted_pdfs[[dimension]][, i] <- as.vector(pdf_matrix %*% self$hidden_state_forecast[, i])
        }
      }
    }
  ),

  # Private methods
  private = list(
    validate_params = function(
      hmm = NULL,
      n = NULL,
      forecast_data = NULL,
      preset_x_vals = NULL,
      starting_state_distribution = NULL
    ) {
      # Check hmm
      if (is.null(hmm)) {
        stop("hmm must be provided")
      }
      if (!inherits(hmm, "HMM")) {
        stop("hmm must be an object of class 'HMM'")
      }

      # Check that n or forecast_data is provided
      if (is.null(n) && is.null(forecast_data)) {
        stop("One of either n or forecast_data must be provided")
      }
      if (!is.null(n)) {
        if (!is.numeric(n) || n <= 0) {
          stop("n must be a positive integer")
        }
      }
      if (!is.null(forecast_data)) {
        if (!is.data.frame(forecast_data)) {
          stop("forecast_data must be a data.frame")
        }
      }

      # Check hmm model for covariates
      covariates <- unique(c(
        sapply(hmm$hid()$formulas(), all.vars),
        rapply(hmm$obs()$formulas(), all.vars)
      ))
      if (length(covariates) > 0) {
        if (is.null(forecast_data)) {
          stop("forecast_data must be provided when covariates are present")
        }
        missing_covs <- setdiff(covariates, colnames(forecast_data))
        if (length(missing_covs) > 0) {
          stop(
            sprintf(
              "forecast_data is missing columns for the following covariates: %s",
              paste(missing_covs, collapse = ", ")
            )
          )
        }
      }

      # Check for values of where to forecast the distribution
      if (!is.null(preset_x_vals)) {
        if (!is.list(preset_x_vals) || is.null(names(preset_x_vals)) || any(names(preset_x_vals) == "")) {
          stop("preset_x_vals must be a named list")
        }
      }

      # Check starting state distribution
      if (is.null(starting_state_distribution)) {
        stop("starting_state_distribution must be provided")
      } else if (
        is.character(starting_state_distribution)
        && !(starting_state_distribution %in% c("last", "stationary")) ) {
        stop("If starting_state_distribution is a string, it must be either 'last' or 'stationary'")
      } else if (!is.array(starting_state_distribution)) {
        stop("starting_state_distribution must be a numeric array")
      } else if (length(starting_state_distribution) != hmm$hid()$n_states()) {
        stop("starting_state_distribution must have the same length as the number of hidden states")
      }
    },

    configure_x_vals = function(
      x_vals = NULL, 
      obs_vars = NULL, 
      data = NULL
    ) {

      # If x_vals is NULL, create an empty list
      if (is.null(x_vals)) {
        x_vals <- list()
      }

      for (obs_var in obs_vars) {
        # If the observation variable is not in x_vals, set as pm 10%
        if (is.null(x_vals[[obs_var]])) {
          max_range <- max(data[[obs_var]], na.rm = TRUE) * 1.1
          min_range <- min(data[[obs_var]], na.rm = TRUE) * 0.9
          x_vals[[obs_var]] <- seq(min_range, max_range, length.out = 100)
        }
      }

      return(x_vals)
    }
  )
)