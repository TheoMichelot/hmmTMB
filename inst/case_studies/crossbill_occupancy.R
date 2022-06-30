# Crossbill data 
# Occupancy model 

library(hmmTMB)

# Data --------------------------------------------------------------------

# load data from package
data(crossbill, package = "unmarked")

# get into format for hmmTMB
nsites <- nrow(crossbill)
dat <- data.frame(ID = rep(1:nsites, each = 9),
                  year = rep(1:9, nsites), 
                  elev = rep(crossbill$ele, each = 9), 
                  forest = rep(crossbill$forest, each = 9), 
                  surveys = rep(crossbill$surveys, each = 9))
y <- apply(as.matrix(crossbill[,5:31]), 1, FUN = function(x) {tapply(as.numeric(x), rep(1:9, each = 3), FUN = function(r) {sum(r, na.rm = TRUE)})})
y <- as.numeric(y)
dat$y2 <- ifelse(dat$surveys == 2, y, NA)
dat$y3 <- ifelse(dat$surveys == 3, y, NA)
dat$forest <- as.numeric(dat$forest)
dat$elev <- as.numeric(dat$elev)

# Fit Models --------------------------------------------------------------

cat("
DATA
dataset = dat
nstates = 2 # not occupied and occupied 

DISTRIBUTION
y2 ~ binom
y3 ~ binom 

INITIAL
y2:
  size = 2, 2
  prob = 0, 0.5
y3:
  size = 3, 3
  prob = 0, 0.5
  
FIXED
y2.prob.state1.(Intercept)
y3.prob.state1.(Intercept)
    
", file = "cross.hmm")

m0 <- HMM$new(file = "cross.hmm")
m0$fit()
m0$par()

## I am going to assume detection is constant across time and 
## all other covariates. 

## Effects on occupancy 
# effect of year on colonisation
mt <- update(m0, "hid", i = 1, j = 2, ~.+s(year, bs = "cs", k = 9))
# effect of year on extinction
mt2 <- update(mt, "hid", i = 2, j = 1, ~.+s(year, bs = "cs", k = 9))
AIC(m0, mt, mt2)

# effect of forest on colonisation  
m_f <- update(mt2, "hid", i = 1, j = 2, ~.+s(forest, bs = "cs"))
# effect of forest on extinction
m_f2 <- update(m_f, "hid", i = 2, j = 1, ~.+s(forest, bs = "cs"))
AIC(m0, m_f, m_f2)

# effect of elevation on colonisation
m_e <- update(m_f2, "hid", i = 1, j = 2, ~.+s(elev, bs = "cs"))
# effect of elevation on extinction
m_e2 <- update(m_e, "hid", i = 2, j = 1, ~.+s(elev, bs = "cs"))
AIC(m_e2, m_e, m_f2)

# interaction between elevation and forest
m_int <- update(mt2, "hid", i = 1, j = 2, ~.+t2(forest, elev, bs = "cs"))
m_int$fit()
AIC(m_e, m_int)

m_e2$plot("tpm", "forest")
m_e2$plot("tpm", "elev")
m_e2$plot("delta", "forest", i = 2)
m_e2$plot("delta", "elev", i = 2)



