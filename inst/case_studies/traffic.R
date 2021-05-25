# Daily traffic causalties by car accidents in spain 
# See: Josep A. Sanchez-Espigares and Alberto Lopez-Moreno (2018). MSwM: Fitting Markov
#       Switching Models. R package version 1.4. https://CRAN.R-project.org/package=MSwM

library(hmmTMB)

# Data --------------------------------------------------------------------

data(traffic, package = "MSwM")

# plot data
plot(traffic$Temp, traffic$NDead)
plot(traffic$Prec, traffic$NDead)
hist(traffic$NDead)

# Fit Models --------------------------------------------------------------

## 2 states
cat("
DATA
dataset = traffic
nstates = 2

DISTRIBUTION
NDead ~ pois 

INITIAL
NDead:
  lambda = 2, 6

", file = "traffic2.hmm")

m2 <- HMM$new(file = "traffic2.hmm")
m2$fit()

## 3 states
cat("
DATA
dataset = traffic
nstates = 3

DISTRIBUTION
NDead ~ pois 

INITIAL
NDead:
  lambda = 1, 4, 6

", file = "traffic3.hmm")

m3 <- HMM$new(file = "traffic3.hmm")
m3$fit()
m3$fit()
m3$par()
AIC(m2, m3) # 2 states selected

## Temperature effect? 
m_temp <- update(m2, "obs", "NDead", "lambda", ~.+s(Temp, bs = "cs"))
AIC(m2, m_temp) # clear evidence of an effect 

## Precipitation effect?
m_pre <- update(m_temp, "obs", "NDead", "lambda", ~.+s(Prec, bs = "cs"))
AIC(m_temp, m_pre) # clear evidence of a preciptation effect

## Interaction? 
m_int <- update(m2, "obs", "NDead", "lambda", ~.+t2(Temp, Prec, bs = "cs"))
AIC(m_temp, m_int)

# Inference ---------------------------------------------------------------

## states
states <- m_int$viterbi()
m_int$plot_ts("NDead")

## distributions
m_int$plot_dist("NDead")

## covariate effects 
m_int$plot("obspar", "Temp")
m_int$plot("obspar", "Prec")






