# Lydia Pinkham Annual Advertising and Sales data 
# See Langrock, R., Kneib, T., Glennie, R. and Michelot, T., 2017. 
#     Markov-switching generalized additive models. Statistics and Computing, 
#     27(1), pp.259-270.

library(hmmTMB)

# Data --------------------------------------------------------------------

# load data from R package
data(pinkham, package = "mAr")

# create a lag variable and remove first data point
dat <- pinkham[-1,]
dat$sales1 <- pinkham$sales[-nrow(pinkham)]

# scale all variables
dat <- as.data.frame(scale(dat))

# plot it 
plot(dat$sales, type = "l")
plot(dat$advertising, dat$sales)
plot(dat$sales1, dat$sales)

# get rough clusters
k <- kmeans(dat$sales, centers = 2)
k$centers
tapply(dat$sales, k$cluster, sd)

# Linear model ------------------------------------------------------------

cat(" 
DATA 
dataset = dat
nstates = 2

DISTRIBUTION
sales ~ norm

INITIAL
sales:
  mean = 1, -0.5
  sd = 0.7, 0.5 
  
FORMULA
sales: 
  mean ~ advertising + sales1
", file = "pinkham.hmm")

m <- HMM$new(file = "pinkham.hmm")
m$fit()
m$coeff_fe()
states <- m$viterbi()
m$plot_ts("sales")
m$plot("obspar", "advertising", i = "sales.mean")


# GAM model ---------------------------------------------------------------

cat(" 
DATA 
dataset = dat
nstates = 2

DISTRIBUTION
sales ~ norm

INITIAL
sales:
  mean = 1, -0.5
  sd = 0.7, 0.5 
  
FORMULA
sales: 
  mean ~ s(advertising, bs = 'cs', k = 5) + sales1
", file = "pinkham_gam.hmm")

m_gam <- HMM$new(file = "pinkham_gam.hmm")
m_gam$fit()
m_gam$plot("obspar", "advertising", i = "sales.mean")


