library(expm)

context("Forecast Dists")

dists <- list(
  beta_obs       = "beta",
  binom_obs      = "binom",
  exp_obs        = "exp",
  gamma_obs      = "gamma",
  lnorm_obs      = "lnorm",
  norm_obs       = "norm",
  pois_obs       = "pois"
)

params <- list(
  beta_obs       = list(shape1   = c( 2,  5),   shape2   = c( 5, 2)),
  binom_obs      = list(size     = c(10, 20),   prob     = c(0.50, 0.30)),
  exp_obs        = list(rate     = c( 1,  3)),
  gamma_obs      = list(shape    = c( 2,  5),   scale    = c(1, 2)),
  lnorm_obs      = list(meanlog  = c( 0,  1),   sdlog    = c(1, 0.5)),
  norm_obs       = list(mean     = c( 0,  3),   sd       = c(1, 2)),
  pois_obs       = list(rate     = c( 2,  7))
)

evaluation_vals <- list(
  beta_obs       = seq(0.05, 0.995, by = 0.01),
  binom_obs      = 0:20,
  exp_obs        = seq(0.05, 10, by = 0.1),
  gamma_obs      = seq(0.05, 30, by = 0.1),
  lnorm_obs      = seq(0.01, 30, by = 0.1),
  norm_obs       = seq(-10, 10, by = 0.1),
  pois_obs       = 0:20
)

n_training <- 2
n_forecast <- 2
tpm <- matrix(c(0.15, 0.9, 0.85, 0.1), 2, 2)
starting_state <- c(0.9, 0.1)


training_df <- data.frame(
  ID = rep(1, n_training),
  beta_obs = c(NA, NA),
  binom_obs = c(NA, NA),
  exp_obs = c(NA, NA),
  gamma_obs = c(NA, NA),
  lnorm_obs = c(NA, NA),
  norm_obs = c(NA, NA),
  pois_obs = c(NA, NA)
)

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

forecast <- Forecast$new(
  hmm = true_mod,
  forecast_data = training_df,
  preset_eval_range = evaluation_vals,
  starting_state_distribution = starting_state %*% tpm
)

for (dist in names(dists)) {
  for (i in seq_len(n_forecast)) {
    hid_state <- starting_state %*% (tpm %^% i)
    
    # Compute theoretical mean
    par <- params[[dist]]
    if (dists[[dist]] == "beta") {
      means <- par$shape1 / (par$shape1 + par$shape2)
    } else if (dists[[dist]] == "binom") {
      means <- par$size * par$prob
    } else if (dists[[dist]] == "exp") {
      means <- 1 / par$rate
    } else if (dists[[dist]] == "gamma") {
      means <- par$shape * par$scale
    } else if (dists[[dist]] == "lnorm") {
      means <- exp(par$meanlog + par$sdlog^2 / 2)
    } else if (dists[[dist]] == "norm") {
      means <- par$mean
    } else if (dists[[dist]] == "pois") {
      means <- par$rate
    }
    
    theory_mean <- means %*% t(hid_state)
    
    dens <- forecast$forecast_dists()[[dist]][, i]
    vals <- evaluation_vals[[dist]]
    is_discrete <- dists[[dist]] %in% c("binom", "pois")
    
    if (is_discrete) {
      forecast_mean <- sum(dens * vals)
    } else {
      delta <- vals[2] - vals[1]
      forecast_mean <- sum(dens * vals * delta)
    }
    expect_equal(forecast_mean[1], theory_mean[1], tolerance = 0.01)
  }
}