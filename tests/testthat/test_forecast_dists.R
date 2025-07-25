library(testthat)
library(devtools)

load_all("../../../hmmTMB")

context("Testing forecast distributions")

full_params <- readRDS("data/full_params.rds")

dists <- full_params[["dists"]]
params <- full_params[["params"]]
evaluation_vals <- full_params[["evaluation_vals"]]

### Training Data Preparation
n_training <- 2
training_df <- data.frame(
  ID = rep(1, n_training),
  matrix(rep(NA, n_training * length(dists)), nrow = n_training, dimnames = list(NULL, names(dists)))
)
training_df$mvnorm_obs <- I(replicate(n_training, c(NA, NA), simplify = FALSE))

tpm <- matrix(c(0.15, 0.9, 0.85, 0.1), 2, 2)
starting_state <- c(0.9, 0.1)
# crafted such that hidden state at n=1 = (0.225 0.775)
# and hidden state at n=2 = (0.73125 0.26875)

## Create True Model ----------------------------------------------------
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

## Create Forecast Object ---------------------------------------------
# Use the true model to create a forecast object
n_forecast <- 2
forecast_df <- data.frame(
  ID = rep(1, n_forecast),
  matrix(rep(NA, n_forecast * length(dists)), nrow = n_forecast, dimnames = list(NULL, names(dists)))
)
forecast_df$mvnorm_obs <- I(replicate(n_forecast, c(NA, NA), simplify = FALSE))

forecast <- Forecast$new(
  hmm = true_mod,
  forecast_data = forecast_df,
  preset_eval_range = evaluation_vals,
  starting_state_distribution = starting_state %*% tpm,
)
forecast_dists <- list()
for (obs in names(dists)) {
  forecast_dists[[obs]] <- forecast$forecast_dists()[[obs]]
}
simulated_pdfs <- readRDS("data/simulated_pdfs.rds")

for (obs in names(dists)) {
  if (obs %in% c("dir_obs", "mvnorm_obs", "tweedie_obs")) { next }
  test_that(paste("Forecasted PDFs match simulated PDFs for ", dists[[obs]]), {
    for (i in seq_len(n_forecast)) {
      sim_pdf <- simulated_pdfs[[obs]][, i]
      forecast_pdf <- forecast_dists[[obs]][, i]
      
      # Ignore first and last elements which are biased by edge effects
      if (obs %in% c("exp_obs", "truncnorm_obs")) {
        sim_pdf <- sim_pdf[-c(1, length(sim_pdf))]
        forecast_pdf <- forecast_pdf[-c(1, length(forecast_pdf))]
      }
      
      ks <- ks.test(sim_pdf, forecast_pdf, simulate.p.value = TRUE, B = 1e4)
      expect_true(
        ks$p.value > 0.3,
        info = paste("KS test failed for", obs, "at forecast time step", i)
      )
    }
  })
}