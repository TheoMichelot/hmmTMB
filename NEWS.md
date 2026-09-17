
# hmmTMB 1.1.3

- Generalise the smoothing penalty to mgcv's `L` convention, so that a smooth
  may combine several penalty matrices through fewer smoothing parameters,
  with `log(lambda) = L * theta`. A smooth may also supply its own starting
  values and parameter names. `L` is the identity and the starting values are
  `lambda = 1` for every basis mgcv ships, so nothing changes for an existing
  model.
- Add the banded forward algorithm of Fischer (2026), through the new `bw`
  argument of `HMM$new()`, with `HMM$update_bw()` and `HMM$check_bw()`. It
  truncates the conditioning of each log-likelihood contribution at a fixed
  lag, which makes the Hessian with respect to latent variables spread over
  time banded rather than dense, at about twice the cost per evaluation and
  with an error that decays geometrically in the bandwidth. Off by default,
  so the exact forward algorithm is unchanged.
- `HMM$post_coeff()` now samples a model with random effects through a sparse
  Cholesky factorisation of the joint precision, rather than inverting it
  first. The inverse is dense even when the precision is not, so this is much
  cheaper and much smaller in memory for a model with many random effects.
  Models without random effects are unaffected.
- Build prediction matrices without fitting a throwaway `mgcv::gam()` to a
  dummy response. `predict.gam(type = "lpmatrix")` uses nothing that fitting
  produces, so the unfitted setup is enough.
- Fix the log-determinants of the penalty matrices when one linear predictor
  contains several smooth terms. `make_matrices()` returned one determinant per
  formula rather than one per penalty, so the likelihood read past the end of
  `log_det_S`, and the smoothing parameters of such a model were biased.
- New vignette on (semi-)supervised learning
- Fix parameter counts for models with constraints
- Use safe Hessian inversion even for models without random effects in `HMM$post_coeff()` and `HMM$confint()`
- Forecasting through the `Forecast` class

# hmmTMB 1.1.2

- Add hurdle negative binomial distribution
- Fix bug in `HMM$fit_stan()`
- Fix bug in `HMM$pseudores()`

# hmmTMB 1.1.1

- Add error message for `laplace = TRUE` in `HMM$fit_stan()`
- Add reference to JSS paper

# hmmTMB 1.1.0

- Add standard error to output of `HMM$confint()`
- Use ID-specific `delta0` when `initial_state = "stationary"`
- Add `par_alt()` functions for pretty display of cat and mvn parameters.
- Improve mvn distribution (automatically detect dimension), and fix mvnorm for dim > 2, using a Cholesky decomposition to get unconstrained parameters of covariance matrix.
- Add zero-one-inflated beta distribution
- Replace optimx by nlminb for model fitting
- Allow empty models for simulation
- Use generalized determinant for penalty matrices 
- Add nbinom(mean, shape) distribution
- Allow for user-specified arguments to be passed directly to mgcv::gam(); e.g., knots.
- Add built-in functions dvm, rvm, dwrpcauchy, rwrpcauchy, to remove dependence on CircStats (as requested by CRAN).

# hmmTMB 1.0.2

- Fix bug for categorical distribution
- Fix bug when using tibble input data
- Fix bugs for models with shared parameters
- Add `HMM$suggest_initial()`
