# Example analysis from 
#   McClintock, B.T., Langrock, R., Gimenez, O., Cam, E., Borchers, D.L., 
#   Glennie, R. and Patterson, T.A., 2020. Uncovering ecological state dynamics 
#   with hidden Markov models. Ecology letters, 23(12), pp.1878-1903.

# Data is on vertical sum of dynamic body acceleration (VDBA) of a striated 
# caracara (Phalcoboenus australis) measured every second over one hour. 

library(hmmTMB)

# Data --------------------------------------------------------------------

# Data set can be downloaded from 
# https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fele.13610&file=ele13610-sup-0003-Example.RData
download.file("https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fele.13610&file=ele13610-sup-0003-Example.RData", 
              destfile = "caracara.RData")

# load data
load("caracara.RData")
str(VDBA)

# put data in format for hmmTMB 
dat <- data.frame(ID = 1, y = VDBA)

# plot data 
plot(dat$y, type = "l")


# Starting values ---------------------------------------------------------

# get rough state means and sds 
cluster <- kmeans(dat$y, centers = 4)
means <- cluster$centers
sds <- tapply(dat$y, cluster$cluster, sd)

# transform to shape/scale for gamma distribution
scale <- sds^2 / means[,1] 
shape <- means[,1] / scale

# Fit Model ---------------------------------------------------------------

# specify model 
mod_spec <- cat(" 
DATA
dataset = dat
nstates = 4

DISTRIBUTION
y ~ gamma 

INITIAL
y: 
  shape = 33, 31, 1, 14
  scale = 0.01, 0.02, 0.01, 0.01
tpm:
  0.95, 0.02, 0.02, 0.01
  0.02, 0.95, 0.02, 0.01
  0.01, 0.02, 0.95, 0.02
  0.01, 0.02, 0.02, 0.95
", 
file = "caracara.hmm")

# create model object
mod <- HMM$new(file = "caracara.hmm")

# fit model
mod$fit()

# Alternative models ------------------------------------------------------

# Try models with 3 or 2 states 

mod3_spec <- # specify model 
  mod_spec <- cat(" 
DATA
dataset = dat
nstates = 3

DISTRIBUTION
y ~ gamma 

INITIAL
y: 
  shape = 33, 1, 14
  scale = 0.01, 0.01, 0.01
tpm:
  0.95, 0.02, 0.03
  0.02, 0.95, 0.03
  0.02, 0.03, 0.95
", 
file = "caracara3.hmm")

mod3 <- HMM$new(file = "caracara3.hmm")
mod3$fit()

mod2_spec <- # specify model 
  mod_spec <- cat(" 
DATA
dataset = dat
nstates = 2

DISTRIBUTION
y ~ gamma 

INITIAL
y: 
  shape = 33, 1
  scale = 0.01, 0.01
tpm:
  0.95, 0.05,
  0.05, 0.95
", 
                  file = "caracara2.hmm")

mod2 <- HMM$new(file = "caracara2.hmm")
mod2$fit()

# compare the models
AIC(mod2, mod3, mod)
BIC(mod2, mod3, mod)
# prefers model with 4 states 

# Model checking ----------------------------------------------------------

# residuals
resids <- mod$pseudores()
qqnorm(resids); qqline(resids)

# goodness-of-fit test by simulation
# quantiles
quantile_stat <- function(x) {
  quantile(x$y, prob = seq(0.1, 0.9, 0.1))
}
# compute for observed data to test fn 
quantile_stat(dat)
# run goodness-of-fit test 
quantile_sims <- mod$gof(quantile_stat, nsims = 100, silent = FALSE)

# autocorrelation 
autocor_stat <- function(x) {
  acf(x$y, plot = FALSE)$acf[,,1]
}
autocor_stat(dat)
autocor_sims <- mod$gof(autocor_stat, nsims = 100, silent = FALSE)

# Inference ---------------------------------------------------------------

# stationary distribution
round(mod$hidden()$delta(), 2)

# observation means and sds
means <- apply(mod$par()$obspar[,,1], 2, prod)
sds <- apply(mod$par()$obspar[,,1], 2, FUN = function(x) {sqrt(x[1]) * x[2]})
round(means, 3)
round(sds, 3)

# transition probability matrix 
round(mod$par()$tpm[,,1], 2)


# State Inference ---------------------------------------------------------

# global decoding 
states <- mod$viterbi() 

# local decoding 
state_pr <- mod$state_probs()
local_states <- apply(state_pr, 1, which.max) 

# plot states over data
mod$plot_ts("y")














