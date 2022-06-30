# Old Faitful Geyser 

library(hmmTMB)

# Data --------------------------------------------------------------------

# load data
data(geyser, package = "MASS")

# look at data
hist(geyser$waiting)
hist(geyser$duration)
plot(geyser$waiting, geyser$duration)

# Fit Models --------------------------------------------------------------

## fit univariate normal distribution
## two states
cat("
DATA
dataset = geyser
nstates = 2

DISTRIBUTION
waiting ~ norm
duration ~ norm

INITIAL
waiting:
  mean = 55, 80
  sd = 2, 2
duration:
  mean = 2, 5
  sd = 1, 1
", file = "geyser_norm.hmm")

m2 <- HMM$new(file = "geyser_norm.hmm")
m2$fit()

m2$obs()$plot_dist("waiting", weights = m2$hid()$delta())
m2$obs()$plot_dist("duration", weights = m2$hid()$delta())

## try a bi-variate normal distribution of duration and waiting time 
# to do that need data in vector format 
y <- NULL
for (i in 1:nrow(geyser)) {
  y[[i]] <- c(geyser$waiting[i], geyser$duration[i])
}
geyser$y <- y 

cat("
DATA
dataset = geyser
nstates = 2

DISTRIBUTION
y ~ mvnorm2

INITIAL
y: 
  mu1 = 60, 82
  mu2 = 4.4, 6.3
  sd1 = 11, 6
  sd2 = 0.35, 1.0
  corr12 = -0.1, -0.5
", file = "mvnorm.hmm")

m_2d <- HMM$new(file = "mvnorm.hmm")
m_2d$fit()
m_2d$par()



