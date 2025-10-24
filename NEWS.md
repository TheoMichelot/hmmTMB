
# hmmTMB 1.1.2

- Fix bug in HMM$fit_stan()

# hmmTMB 1.1.1

- Add error message for laplace = TRUE in HMM$fit_stan()
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
