# Epileptic Seizures 
# This example is introduced on p135 of 
# Zucchini, W., MacDonald, I.L. and Langrock, R., 2017. 
# Hidden Markov models for time series: an introduction using R. CRC press.
# See that chapter for details. 

# Data --------------------------------------------------------------------

epil <- c(0, 3, 0, 0, 0, 0, 1, 1, 0, 2, 1, 1, 2, 0, 0, 1, 2, 1, 3, 1, 3,
          0, 4, 2, 0, 1, 1, 2, 1, 2, 1, 1, 1, 0, 1, 0, 2, 2, 1, 2, 1, 0,
          0, 0, 2, 1, 2, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
          0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0,
          0, 2, 1, 0, 1, 1, 0, 0, 0, 2, 2, 0, 1, 1, 3, 1, 1, 2, 1, 0, 3,
          6, 1, 3, 1, 2, 2, 1, 0, 1, 2, 1, 0, 1, 2, 0, 0, 2, 2, 1, 0, 1,
          0, 0, 2, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 3,
          0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 2, 1, 0,
          0, 0, 0, 0, 0, 1, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
          
# plot data over time 
plot(epil, type = "b", pch = 19, cex = 0.5)

# put into a dataframe
dat <- data.frame(epil = epil)

# One state model ---------------------------------------------------------

m1 <- glm(epil ~ 1, family = "poisson", data = dat)
# mean with uncertainty 
p <- predict(m1, se = TRUE, type = "response")
p$fit[1]
p$fit[1] + c(-1, 1) * qnorm(0.975) * p$se.fit[1]



# Higher state models -----------------------------------------------------

## Two states
cat("
DATA
dataset = dat
nstates = 2

DISTRIBUTION
epil ~ pois 

INITIAL
epil:
  lambda = 0.5, 1 
    
", file = "epil_2.hmm")

m2 <- HMM$new(file = "epil_2.hmm")
m2$fit()

## Three states
cat("
DATA
dataset = dat
nstates = 3

DISTRIBUTION
epil ~ pois 

INITIAL
epil:
  lambda = 0.1, 0.5, 1.5 
    
", file = "epil_3.hmm")

m3 <- HMM$new(file = "epil_3.hmm")
m3$fit()

## Four states
cat("
DATA
dataset = dat
nstates = 4

DISTRIBUTION
epil ~ pois 

INITIAL
epil:
  lambda = 0.5, 1, 2, 4 
    
", file = "epil_4.hmm")

m4 <- HMM$new(file = "epil_4.hmm")
m4$fit()

# Compare models ----------------------------------------------------------
# This reproduces Table 9.2
mods <- list(m1,m2,m3,m4)
# row of that table for each model
rows <- lapply(mods, FUN = function(m) {
  l <- logLik(m)
  c(attr(l, "df"), -l, AIC(m), BIC(m))
})
tab <- round(cbind(1:4, do.call("rbind", rows)), 2)
colnames(tab) <- c("no. of states", "k", "-l", "AIC", "BIC")
tab

# Look at 2-state model  --------------------------------------------------

m2$par()

# Model checking ----------------------------------------------------------

# can look at ACF 
acf_stat <- function(dat) {
  acf(dat, plot = FALSE)$acf[,1,1]
}
# observed
acf_stat(dat)
# do simulations from model
acf_sims <- m2$gof(acf_stat, nsims = 100)

# can look at observed v predicted days with certain number of 
# seizures 
pred_stat <- function(dat) {
  c(sum(dat$epil == 0), 
    sum(dat$epil == 1), 
    sum(dat$epil == 2), 
    sum(dat$epil == 3), 
    sum(dat$epil == 4), 
    sum(dat$epil == 5),
    sum(dat$epil > 5))
}
pred_stat(dat)
pred_sims <- m2$gof(pred_stat, nsims = 100)
