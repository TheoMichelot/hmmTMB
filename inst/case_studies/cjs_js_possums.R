# Cormack-Jolly-Seber and Jolly-Seber analysis of Brushtail Possums
# from openCR package, analyzed by Pledger et al. (2003) and available from
# Murray Efford 

library(hmmTMB)

# Data --------------------------------------------------------------------

# load data from package
library(openCR)
ch <- FebpossumCH
ch <- join(reduce(FebpossumCH, by = "all", verify = FALSE))
n <- dim(ch)[1]
nocc <- dim(ch)[2]

# get into format for hmmTMB 
dat <- data.frame(ID = rep(1:n, each = nocc), occ = rep(1:nocc, n))
dat$cap <- as.numeric(t(matrix(as.numeric(ch[,,1]), nr = 448, nc = 9)))
dat$sex <- covariates(ch)$sex[dat$ID]

# reduce to post-first detection only for cjs analysis
cjs_dat <- dat
cjs_dat$state <- NA
for (i in 1:n) {
  wh <- min(which(cjs_dat$cap[cjs_dat$ID == i] == 1)) - 1 
  start <- min(which(cjs_dat$ID == i))
  if (wh > nrow(ch)) {
    cjs_dat <- cjs_dat[cjs_dat$ID != i,]
  } else {
    cjs_dat$cap[start + wh] <- NA
    cjs_dat$state[start + wh] <- 1
    if (wh > 0) {
      cjs_dat <- cjs_dat[-(start:(start + wh - 1)),]
    }
  }
}

# CJS analysis ------------------------------------------------------------

cat("
DATA
dataset = cjs_dat
nstates = 2 # alive and dead 

DISTRIBUTION
cap ~ binom # seen or not seen 

INITIAL
cap: 
  size = 1, 1
  prob = 0.5, 0 
tpm:
  0.8, 0.2 
  0.0, 1.0 # if you are dead you stay dead 
delta:
  1.0, 0.0 
  
TPM
. ; ~ 1
. ; .  # dot means these probs are fixed at their initial value 

FIXED
delta
cap.prob.state2.(Intercept)
  
", file = "cjs.hmm")

cjs <- HMM$new(file = "cjs.hmm")
cjs$fit()
cjs$par()

# include time-varying detection
cjs_pt <- update(cjs, "obs", "cap", "prob", ~.+state1(s(occ, k = 9, bs = "cs")))
AIC(cjs_pt, cjs)

# include sex effect on detection
cjs_sex <- update(cjs, "obs", "cap", "prob", ~.+state1(sex))
AIC(cjs_sex, cjs)

# include time-varying survival
cjs_phit <- update(cjs_sex, "hid", 1, 2, ~.+s(occ, k = 9, bs = "cs"))
AIC(cjs_sex, cjs_phit)

# Try a finite mixture with 2 subpopulations 
# "alive" observation is now one of two possible states 
cjs_dat$state[cjs_dat$state == 1] <- "1,2"

cat("
DATA
dataset = cjs_dat
nstates = 3 # alive subpop1, alive supop2, and dead 

DISTRIBUTION
cap ~ binom # seen or not seen 

INITIAL
cap: 
  size = 1, 1, 1
  prob = 0.9, 0.1, 0 
tpm:
  0.8, 0.0, 0.2 
  0.0, 0.5, 0.5 
  0.0, 0.0, 1.0
delta:
  0.5, 0.4, 0.1 # initial equal chance of being in each subpop 
  
FORMULA
cap:
  prob ~ state1(sex) + state2(sex)
  
TPM
. ; . ; ~1
. ; .  ; ~1
. ; .  ; .  # dot means these probs are fixed at their initial value 

FIXED
cap.prob.state3.(Intercept)
  
", file = "cjs_h2.hmm")

cjs_h2 <- HMM$new(file = "cjs_h2.hmm")
cjs_h2$fit()
# sub-population probabilities
cjs_h2$par()
# sub-population proportions 
round(cjs_h2$hid()$delta(), 2)
AIC(cjs_h2, cjs_sex)



