# Crossbill data 
# Occupancy model 

library(hmmTMB)

# Data --------------------------------------------------------------------

# load data from package
data(crossbill, package = "unmarked")

# get into format for hmmTMB
nsites <- nrow(crossbill)
dat <- data.frame(ID = rep(1:nsites, each = 9),
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
 
# effect of forest on turnover
m_f <- update(m0, "hidden", i = 1, j = 2, ~.+s(forest, bs = "cs", k = 5))
AIC(m0, m_f)


