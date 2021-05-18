# Spanish energy prices example 
# See Langrock, R., Kneib, T., Glennie, R. and Michelot, T., 2017. 
#     Markov-switching generalized additive models. Statistics and Computing, 
#     27(1), pp.259-270.

library(hmmTMB)

# Data --------------------------------------------------------------------

# load data from a package 
data(energy, package = "MSwM")

# look at relationship
plot(energy$EurDol, energy$Price)

# scale data 
energy <- as.data.frame(scale(energy))

# get rough clusters
k <- kmeans(energy$Price, centers = 2)
k$centers
tapply(energy$Price, k$cluster, sd)

# Fit model ---------------------------------------------------------------

cat("
DATA
dataset = energy
nstates = 2

DISTRIBUTION
Price ~ norm

INITIAL
Price:
  mean = -0.5, 1
  sd = 0.5, 0.6

FORMULA
Price: 
  mean ~ s(EurDol, bs = 'cs')
    
", file = "energy.hmm")

m <- HMM$new(file = "energy.hmm")
m$fit()


states <- m$viterbi()
energy$state <- factor(paste0("State ", states))
m$plot("obspar", "EurDol", i = "Price.mean") + 
  geom_point(data = energy, aes(x = EurDol, y = Price, fill = state, col = state), alpha = 0.3)


