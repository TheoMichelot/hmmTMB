# Dipper Cormack-Jolly-Seber model example 

library(hmmTMB)

# Data --------------------------------------------------------------------

# load data from package
data(dipper, package =  "marked")

# get data into hmmTMB format
ch <- sapply(strsplit(dipper$ch, ""), FUN = as.numeric)
dat <- data.frame(ID = rep(1:ncol(ch), each = nrow(ch)), 
                  cap = as.numeric(ch), 
                  occ =  rep(1:7, nrow(ch)), 
                  sex = rep(dipper$sex, each = nrow(ch)))

dat$state <- NA
# get first detection of each individual
# and take out data before first detection 
for (i in 1:ncol(ch)) {
  wh <- min(which(dat$cap[dat$ID == i] == 1)) - 1 
  start <- min(which(dat$ID == i))
  if (wh > nrow(ch)) {
    dat <- dat[dat$ID != i,]
  } else {
    dat$cap[start + wh] <- NA
    dat$state[start + wh] <- 1
    if (wh > 0) {
      dat <- dat[-(start:(start + wh - 1)),]
    }
  }
}

# Setup CJS model ---------------------------------------------------------

cat("
DATA
dataset = dat
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

# Fit Models --------------------------------------------------------------

# fit basic model
cjs$fit()

# fit time-varying detection and survival
cjs_pt <- update(cjs, "obs", "cap", "prob", ~ . + state1(s(occ, k = 7, bs = "cs")))
cjs_phit <- update(cjs, "hidden", 1, 2, ~.+s(occ, k = 6, bs = "cs"))
cjs_pt_phit <- update(cjs_pt, "hidden", 1, 2, ~.+s(occ, k = 6, bs = "cs"))

AIC(cjs, cjs_pt, cjs_phit, cjs_pt_phit)

# look at sex effect
cjs_psex <- update(cjs, "obs", "cap", "prob", ~. + sex)
AIC(cjs, cjs_psex)
cjs_phisex <- update(cjs_psex, "hidden", 1, 2, ~.+sex)
# won't fit, likely due sex effect being minimal 
