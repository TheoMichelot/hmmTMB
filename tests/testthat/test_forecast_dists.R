library(testthat)
library(devtools)

load_all("../../../hmmTMB")

# Load the cached test data
------------------------------------------
# See inst/devel/create_test_forecast_dists.ipynb for details
# on how this test data was created.
full_params <- readRDS("data/full_params.rds")

dists <- full_params[["dists"]]
params <- full_params[["params"]]
evaluation_vals <- full_params[["evaluation_vals"]]

# Cached simulated PDFs from running 1000 simulations
simulated_pdfs <- readRDS("data/simulated_pdfs.rds")

# Create the model
-----------------------------------------
n_training <- 2
training_df <- data.frame(
  ID = rep(1, n_training),
  matrix(
    rep(NA, n_training * length(dists)),
    nrow = n_training, dimnames = list(NULL, names(dists))
  )
)
training_df$mvnorm_obs <- I(replicate(n_training, c(NA, NA), simplify = FALSE))

tpm <- matrix(c(0.8, 0.3, 0.2, 0.7), 2, 2)
starting_state <- c(0.8, 0.2)
hid_mod <- MarkovChain$new(
  data = training_df,
  n_states = 2,
  tpm = tpm
)
hid_mod$update_delta0(starting_state)
obs_mod <- Observation$new(
  data = training_df,
  n_states = 2,
  dists = dists,
  par = params
)
true_mod <- HMM$new(obs = obs_mod, hid = hid_mod)

# Create the forecast object
-----------------------------------------
n_forecast <- 4
forecast_df <- data.frame(
  ID = rep(1, n_forecast),
  matrix(rep(NA, n_forecast * length(dists)), nrow = n_forecast, dimnames = list(NULL, names(dists)))
)
forecast_df$mvnorm_obs <- I(replicate(n_forecast, c(NA, NA), simplify = FALSE))

forecast = Forecast$new(
  hmm = true_mod,
  forecast_data = forecast_df,
  preset_x_vals = evaluation_vals,
  starting_state_distribution = starting_state,
)

# For multivariate distributions, only keep first dimension
------------------------------------------
forecasted_pdfs <- list()
# loop through each dimension and time step
for (obs in names(dists)) {

  if (obs %in% c("mvnorm_obs")) {

    multivar_pdf <- matrix(NA, nrow = length(unique(forecast$x_vals[[obs]][1, ])), ncol = n_forecast)
    for (i in seq_along(unique(forecast$x_vals[[obs]][1, ]))) {
      val <- unique(forecast$x_vals[[obs]][1, ])[i]
      mask <- forecast$x_vals[[obs]][1, ] == val
      multivar_pdf[i, ] <- apply(forecast$forecasted_pdfs[[obs]][mask, ], 2, sum)
    }
    forecasted_pdfs[[obs]] <- multivar_pdf
  } else {
    forecasted_pdfs[[obs]] <- forecast$forecasted_pdfs[[obs]]
  }
}

for (obs in names(dists)) {
  test_that(paste("Forecasted PDFs match simulated PDFs for ", dists[[obs]]), {
    expect_equal(
      simulated_pdfs[[obs]],
      forecasted_pdfs[[obs]],
      tolerance = 0.1
    )
  })
}

