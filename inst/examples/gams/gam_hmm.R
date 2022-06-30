
library(hmmTMB)

# Simulate data -----------------------------------------------------------
set.seed(13763)

# number of time steps
n <- 2500

# create an empty dataset
dat <- data.frame(ID = rep(1, n), count = rep(0, n))

# create covariate
dat$x <- rnorm(n)

# create true model 
true_mod <- HMM$new(file = "true_gam_mod.hmm")

# set gam coefficients for observed 
par <- true_mod$coeff_re()$obs
S <- true_mod$obs()$make_mat()$S # precision matrix for smooth 
# exclude null space parameters 
samp <- rmvn(n = 1, rep(0, ncol(S) - 2), solve(S[-c(9, 18), -c(9, 18)]))
par <- c(samp[1:8], 0, samp[9:length(samp)], 0)
true_mod$obs()$update_coeff_re(par)

# set gam coefficients for tpm 
tpmpar <- true_mod$coeff_re()$hid
Shid <- true_mod$hid()$make_mat(data = dat)$S
samphid <- rmvn(n = 1, rep(0, ncol(Shid)- 1), solve(0.1 * Shid[-9, -9]))
tpmpar <- c(samphid, 0)
true_mod$hid()$update_coeff_re(tpmpar)

# simulate from true model
set.seed(52320)
dat <- true_mod$simulate(n, data = dat)

plot(dat$x, dat$count)

# Fit model --------------------------------------------

# create model to fit 
mod <- HMM$new(file = "gam_mod.hmm")

# check model formulation
mod

# fit model
mod$fit()


# Parameter Inference -----------------------------------------------------

# look at numerical estimates
mod$par()
 
# plot estimated relationship
pt <- mod$plot("obspar", "x")
pt

# plot true relationship with estimated 
xgr <- seq(-4, 4, length = 50)
pred <- true_mod$predict("obspar", t = 1:50, newdata = data.frame(x = xgr))
library(ggplot2)
library(reshape2)
preddat <- data.frame(t(pred[1,,]))
preddat <- melt(preddat)
preddat$x <- rep(xgr, 2)
pt + geom_line(aes(x = preddat$x, y = preddat$value, group = preddat$variable), size = 1.5, col = "black")

# plot tpm estimated relationship 
tpm_pt <- mod$plot("tpm", "x", i = 1, j = 2)
tpm_pt
tpmpred <- true_mod$predict("tpm", t = 1:50, newdata = data.frame(x = xgr))
tpmpreddat <- data.frame(p = tpmpred[1,2,])
tpmpreddat$x <- xgr
tpm_pt + geom_line(aes(x = tpmpreddat$x, y = tpmpreddat$p), size = 1.5, col = "firebrick")
