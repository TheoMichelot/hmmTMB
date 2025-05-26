#' @title Forecast R6 Class for Hidden Markov Models
#'
#' @description
#' Creates a *forecast* object from a fitted HMM model, projecting
#' future state probabilities and unconditional observation densities for a
#' user-defined horizon or covariate data-frame.
#' The forecasting is executed in the constructor; afterwards the results live
#' in public fields that can be queried or plotted directly.
#'
#' @details
#' **Algorithm overview**
#' \enumerate{
#'   \item Validate inputs (`private$validate_params()`).
#'   \item Build or accept a future design matrix (`forecast_data`).
#'   \item Obtain time-varying transition matrices (`tpm_forecast`) and
#'         observation parameters (`obs_par_forecast`) via `hmm$predict()`.
#'   \item Forward-propagate the hidden-state distribution
#'         (`hidden_state_forecast`).
#'   \item Marginalise over states to derive unconditional pdfs of each
#'         observation variable at a user-defined grid (`forecasted_pdfs`).
#' }
#'
#' @param hmm A fitted HMM object.
#' @param n Integer. Forecast horizon (ignored when `forecast_data` is supplied).
#' @param forecast_data `data.frame` of future covariates / IDs.
#'        Must contain every covariate used in the fitted `hmm`.
#' @param preset_x_vals *Named* list whose elements are numeric vectors giving
#'        the evaluation grid for each observation variable.  
#'        Missing names fall back to a 90–110 % range of the training data.
#' @param starting_state_distribution Numeric vector of length `nstates()`, or
#'        the string `"last"` (use the filtered distribution of the final
#'        training row) or `"stationary"` (use the model’s stationary
#'        distribution).
#'
#' @return An R6 object with the following **notable public fields**:  
#' \describe{
#'   \item{hidden_state_forecast}{`nstates × T` matrix of forward state
#'         probabilities.}
#'   \item{obs_par_forecast}{Array of time-varying observation parameters.}
#'   \item{tpm_forecast}{Array of time-varying transition matrices.}
#'   \item{forecasted_pdfs}{List of unconditional pdf matrices, one per
#'         observation variable.}
#'   \item{x_vals}{List of grids on which each pdf was evaluated.}
#' }
#'
#' @examples
#' \dontrun{
#'   mod <- HMM$new(file = "pois_mod.hmm")
#'   mod$fit()
#'
#'   fc  <- forecast$new(hmm = mod, n = 12)
#'
#'   # Plot the nth forecasted pdf for the first observation variable
#'   step <- 1
#'   ggplot() +
#'   geom_ridgeline(
#'     aes(
#'       x = fc$x_vals[[1]],
#'       y = 1,
#'       height = fc$forecasted_pdfs[[1]][, step]
#'     ),
#'     scale = 0.5
#'   ) +
#'   labs(
#'     title = sprintf("PDF for the %d Time Step", step),
#'     x = sprintf("Observation Value: %s", fc$observation_vars[[1]]),
#'     y = "Probability Density"
#'   )
#' }
#'
#' @export
forecast <- R6::R6Class(
  "forecast",

  ## Public fields -------------------------------------------------------------
  public = list(
    hmm                   = NULL,  # fitted HMM model
    observation_vars      = NULL,  # names of response variables
    forecast_data         = NULL,  # future covariate data-frame
    obs_par_forecast      = NULL,  # predicted observation parameters
    tpm_forecast          = NULL,  # predicted transition matrices
    x_vals                = NULL,  # grid for each response variable
    hidden_state_forecast = NULL,  # forward state probabilities
    forecasted_pdfs       = NULL,  # unconditional predictive pdfs

    #' @description
    #' Construct a forecast object.
    initialize = function(hmm               = NULL,
                          n                 = NULL,
                          forecast_data     = NULL,
                          preset_x_vals     = NULL,
                          starting_state_distribution = "last") {

      ## -- 1  Input checking ---------------------------------------------------
      private$validate_params(
        hmm                       = hmm,
        n                         = n,
        forecast_data             = forecast_data,
        preset_x_vals             = preset_x_vals,
        starting_state_distribution = starting_state_distribution
      )

      ## -- 2  Store the fitted model & basic metadata -------------------------
      self$hmm             <- hmm
      self$observation_vars <- colnames(self$hmm$obs()$obs_var())

      ## -- 3  Build or accept the future design matrix ------------------------
      self$forecast_data <- if (!is.null(forecast_data)) {
        forecast_data
      } else {
        as.data.frame(
          c(
            list(ID = rep(1, n)),
            stats::setNames(
              replicate(length(self$observation_vars), rep(NA, n),
                        simplify = FALSE),
              self$observation_vars
            )
          )
        )
      }

      ## -- 4  Predict forward-looking parameters ------------------------------
      self$obs_par_forecast <- self$hmm$predict("obspar",
                                                newdata = self$forecast_data)
      self$tpm_forecast     <- self$hmm$predict("tpm",
                                                newdata = self$forecast_data)

      ## -- 5  Establish evaluation grids for each response variable -----------
      self$x_vals <- private$configure_x_vals(
        x_vals   = preset_x_vals,
        obs_vars = self$observation_vars,
        data     = self$hmm$obs()$data()
      )

      ## -- 6  Choose the initial hidden-state distribution --------------------
      if (is.character(starting_state_distribution)) {

        if (starting_state_distribution == "last") {
          # one-step-ahead distribution conditional on final training point
          last_sp  <- utils::tail(self$hmm$state_probs(), 1)
          last_tpm <- self$hmm$hid()$tpm(nrow(self$hmm$obs()$data()))[, , 1]
          dist_0   <- last_sp %*% last_tpm
        }

        if (starting_state_distribution == "stationary") {
          dist_0 <- self$hmm$hid()$stationary()
        }

      } else dist_0 <- starting_state_distribution

      ## -- 7  Forward-propagate hidden states ---------------------------------
      n_steps <- nrow(self$forecast_data)
      hidden_state_forecast <-
        array(NA_real_,
              dim = c(self$hmm$hid()$nstates(), n_steps))

      hidden_state_forecast[, 1] <- dist_0

      if (n_steps > 1) {
        for (t in 2:n_steps) {
          hidden_state_forecast[, t] <-
            hidden_state_forecast[, t - 1] %*% self$tpm_forecast[, , t - 1]
        }
      }
      self$hidden_state_forecast <- hidden_state_forecast

      ## -- 8  Build unconditional predictive pdfs -----------------------------
      self$forecasted_pdfs <- vector("list", length(self$observation_vars))
      names(self$forecasted_pdfs) <- self$observation_vars

      for (dimension in self$observation_vars) {

        obs_dists      <- hmm$obs()$dists()[[dimension]]

        pdf_params     <- paste0(dimension, ".", names(formals(obs_dists$pdf())))
        model_params   <- names(self$obs_par_forecast[, 1, 1])
        current_params <- intersect(pdf_params, model_params)

        if (is.null(current_params)) {
          # Edge case: when obs_par_forecast has only 1 parameter it is unnamed
          # Use current_params = 1 to unpack the matrix correctly
          current_params <- 1
        }

        dim_x_vals <- self$x_vals[[dimension]]
        self$forecasted_pdfs[[dimension]] <-
          array(NA_real_, dim = c(length(dim_x_vals), n_steps))

        for (i in seq_len(n_steps)) {

          # Evaluate each state's pdf at x_vals for time step i
          # pdf_matrix: matrix with dimensions |x_vals| x n_states
          pdf_matrix <- vapply(
            seq_len(self$hmm$hid()$nstates()),   # loop over hidden states
            function(s) {                       # compute pdf for state s
              obs_dists$pdf_apply(
                x   = dim_x_vals,
                par = stats::setNames(
                  self$obs_par_forecast[current_params, s, i],
                  obs_dists$parnames()
                )
              )
            },
            numeric(length(dim_x_vals))        # vapply template: numeric vector of length |x_vals|
          )

          # Compute the weighted sum over states to get the unconditional forecasted pdf
          un_normalised_pdf <- pdf_matrix %*% self$hidden_state_forecast[, i]
          self$forecasted_pdfs[[dimension]][, i] <- un_normalised_pdf / sum(un_normalised_pdf)
        }
      }
    }
  ),

  ## --------------------------------------------------------------------------
  ## Private helpers -----------------------------------------------------------
  ## --------------------------------------------------------------------------
  private = list(
    ## -- 1  Comprehensive argument checks -------------------------------------
    validate_params = function(hmm = NULL, n = NULL, forecast_data = NULL,
                               preset_x_vals = NULL,
                               starting_state_distribution = NULL) {

      if (is.null(hmm))
        stop("`hmm` must be provided", call. = FALSE)
      if (!inherits(hmm, "HMM"))
        stop("`hmm` must inherit from class 'HMM'", call. = FALSE)

      # n  vs  forecast_data ----------------------------------------------------
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
      if (!is.null(preset_x_vals) &&
          (!is.list(preset_x_vals) ||
             is.null(names(preset_x_vals)) ||
             any(names(preset_x_vals) == ""))) {
        stop("`preset_x_vals` must be a *named* list", call. = FALSE)
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
    configure_x_vals = function(x_vals   = NULL,
                                obs_vars = NULL,
                                data     = NULL) {

      if (is.null(x_vals)) x_vals <- list()

      for (obs_var in obs_vars) {
        if (is.null(x_vals[[obs_var]])) {
          max_range <- max(data[[obs_var]], na.rm = TRUE) * 1.1
          min_range <- min(data[[obs_var]], na.rm = TRUE) * 0.9
          x_vals[[obs_var]] <- seq(min_range, max_range, length.out = 100)
        }
      }
      x_vals
    }
  )
)
