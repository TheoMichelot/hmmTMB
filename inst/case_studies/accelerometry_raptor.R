# Soaring raptor accelerometry data from 
# Leos-Barajas, Vianey et al. (2017), Data from: Analysis of animal accelerometer data using hidden Markov models, Dryad, Dataset, https://doi.org/10.5061/dryad.6bm2c

library(hmmTMB)
library(lubridate)

# Data --------------------------------------------------------------------

# download zip file from Dryad
download.file("https://uc3-s3mrt1001-prd.s3.us-west-2.amazonaws.com/71275f33-f7fc-4731-b84e-f70cec36759d/data?response-content-disposition=attachment%3B%20filename%3Ddoi_10.5061_dryad.6bm2c__v1.zip&response-content-type=application%2Fzip&X-Amz-Security-Token=IQoJb3JpZ2luX2VjEAoaCXVzLXdlc3QtMiJHMEUCIBRio009i6tvgBgEQYMcWwpMyrMmmT5yC2YHkzahSvGhAiEAgu6coxiLAuUBZjfsS9LEU2H8nQsid2eudeEjl4lhkbcqgwQIo%2F%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FARAAGgw0NTE4MjY5MTQxNTciDB4eQKLswuFJlqsSFyrXAyInnezrfzjBYGdcn5gkqs0CRovmqHROXUCmtHPZ9yoIesq7tcgW2x804aBMq8MSRoDUHjUGASnkAjDfbt8ITPOptIBbAq%2BfKQpYORsY2SiKtv8AE8On4Z7IpCwQt8A1H0J0qZXXPX6%2FJuH7KlxKtKybwZDu1bOfcuorfZPCfwtq0puRxWkZ%2F2H9r2naOBuuIdm08lgOSMeXiZl0vLO6RVmR1jgb8dsmMZB2ouAAvMxZbD3iHh4iErXivNd7RKbzo5s2PUUBkJb0IWp7NAIyDVTrgECWsMeOPo389nkn4coEQdMPjSGZoC%2BNNSEK9O3NplZdd7qf%2B1Ah7Himqqlh8jJwa%2FxpGjtgtEss8Np2lzaGe%2Bm4qsPVEr848x%2F6Tyd5%2B%2BHf532GVxLKfy8BXz0Hj31DVCX1N%2Baw67GrirUp1PHpvn48C9Pbf6ciZy72ga6arCCTaHQ%2FjTkS3bgrD0JQmIezKzs1hu99LuX42gLF6tsL1mWKumLbWoUlb5y4hHrrJbdfeqvPNe6QHPFoUSglBWvLm5LIV84%2FT06fKPLlQ6Q62yu7Qry7vgO15hnSyszgdheOpjiO0fl2Mc%2FFoMv9cOy%2FFvsKC2RNoCT4sLz1%2FxZVbO4lkuX5%2BjC8jp6FBjqlATjp0%2FRFBasfKVz6yJnx3BNUug4qoh5oHd1SMvpIiT%2B2BEy2nEw0Pz3FIOWDyTNEEvLnnbaGc3Tj4gcr2zbAL%2Bgm%2FEB6TOdvTzEsxd3nEb7r0T80gDjiZ2ZH20YzRsDsLcuMxfucrocoVYh%2FA1ozaHJN3gV9CuLDglaaw3LieF%2FykN2dsZXi3MFclO4jZamcufNr8QVWlqMnxc6fm7LipfAu6jD%2F4g%3D%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20210521T105537Z&X-Amz-SignedHeaders=host&X-Amz-Expires=14400&X-Amz-Credential=ASIAWSMX3SNW6HI7GEMD%2F20210521%2Fus-west-2%2Fs3%2Faws4_request&X-Amz-Signature=c5fa20e1622a5b0b0c69cf072f67b087cbe60f6c9e19686a236ef7a2ce7aeb4f", 
              destfile = "raptor.zip")

# unzip file 
utils::unzip("raptor.zip")

# load datafile 
raw_dat <- read.table("Verreauxs.accel.txt", skip = 1, header = FALSE)
colnames(raw_dat) <- c("date", "time", "msa", "wind_speed", "temp", "ID")
str(raw_dat)
raw_dat$datetime <- ymd_hms(paste(raw_dat$date, raw_dat$time))

# put file into format for hmmTMB
dat <- data.frame(ID = raw_dat$ID, 
                  msa = raw_dat$msa, 
                  wind_speed = raw_dat$wind_speed, 
                  temp = raw_dat$temp, 
                  date = raw_dat$datetime)
dat$hour <- hour(dat$date)
dat$day <- day(dat$date)

# Exploratory -------------------------------------------------------------

# plot data
plot(dat$date, dat$msa)
hist(dat$msa)

# plot relationships
plot(dat$wind_speed, dat$msa)
plot(dat$temp, dat$msa)
plot(dat$hour, dat$msa)
plot(dat$day, dat$msa)

# do a simple 2-state cluster
k <- kmeans(dat$msa, centers = 2)
mu <- k$centers[,1]
sd <- tapply(dat$msa, k$cluster, sd)

# transform for gamma distribution parameters
scale <- sd^2 / mu
shape <- mu / scale 


# Fit Models --------------------------------------------------------------

# Going to restrict to 2-state models 
cat("
DATA
dataset = dat
nstates = 2

DISTRIBUTION
msa ~ gamma

INITIAL
msa:
  shape = 3.5, 0.5
  scale = 0.3, 0.2 
", file = "raptor_m0.hmm")

m0 <- HMM$new(file = "raptor_m0.hmm")
m0$fit()
m0$obs()$plot_dist("msa", weights = m0$hidden()$delta())

## look at effects on tpm
# temperature
m_temp <- update(m0, "hidden", 1, 2, ~.+s(temp, bs = "cs"), fit = FALSE)
m_temp <- update(m_temp, "hidden", 2, 1, ~.+s(temp, bs = "cs"))
AIC(m0, m_temp) # little evidence of an effect 

# try temp linear only 
m_temp2 <- update(m0, "hidden", 1, 2, ~.+temp, fit = FALSE)
m_temp2 <- update(m_temp2, "hidden", 2, 1, ~.+temp)
AIC(m0, m_temp2) # little evidence of an effect 

# wind speed 
m_wind <- update(m0, "hidden", 1, 2, ~.+s(wind_speed, bs = "cs"), fit = FALSE)
m_wind <- update(m_wind, "hidden", 2, 1, ~.+s(wind_speed, bs = "cs"))
AIC(m0, m_wind)

# try wind linear only 
m_wind2 <- update(m0, "hidden", 1, 2, ~.+wind_speed, fit = FALSE)
m_wind2 <- update(m_wind2, "hidden", 2, 1, ~.+wind_speed)
AIC(m0, m_wind2) # some evidence of a linear effect 

# look for hour of day effects
m_hr <- update(m_wind2, "hidden", 1, 2, ~.+s(hour, bs = "cc"), fit = FALSE)
m_hr <- update(m_hr, "hidden", 2, 1, ~.+s(hour, bs = "cc"))
AIC(m_wind2, m_hr) # nothing


# Inference ---------------------------------------------------------------

m_wind2$plot_dist("msa")
m_wind2$plot("tpm", "wind_speed")
m_wind2$plot("delta", "wind_speed")

states <- m_wind2$viterbi()
m_wind2$plot_ts("msa")

# Model checking ----------------------------------------------------------

resids <- m_wind2$pseudores()
qqnorm(resids); qqline(resids)
hist(resids)

