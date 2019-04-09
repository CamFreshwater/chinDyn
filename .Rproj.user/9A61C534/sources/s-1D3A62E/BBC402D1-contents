#*******************************************************************************
# estimateOU.R
# Date last revised 9 Jan 2019
# Inputs: CU-specific catch data and MU-specific forecast data
#*******************************************************************************

require(here); require(tidyverse); require(ggplot2); require(reshape2); 
require(samSim)

seasonDat <- read.csv(here("data", "sox", "catchAndForecastEstimates.csv"), 
                        stringsAsFactors=F) %>% 
  rename(finalCatch = postCatch_manTable)

## Calculate outcome uncertainty as deviation between mid-season TAC and 
# post-season catch
ouDat <- seasonDat %>%
  mutate(targetHR = midTAC / midForecast,
         realizedHR = finalCatch / postForecast) %>% 
  filter(!is.na(realizedHR))

# Replace target HRs from management tables with minimum ERs provided by AMH
# (used to account for incidental harvest on other MUs)
for (k in 1:nrow(ouDat)) {
  if (ouDat$targetHR[k] == 0) {
    ouDat$targetHR[k] <- ifelse(ouDat$manUnit[k] == "Lat", 0.2, 0.1)
  }
}

plot(targetHR ~ realizedHR, data = ouDat)

ouDat <- ouDat %>% 
  mutate(hrDev = realizedHR - targetHR)
hist(ouDat$hrDev)

trimOuDat <- ouDat %>%
  filter(!hrDev < -0.4)
# hist(trimOuDat$hrDev)

lateDat <- trimOuDat %>% 
  filter(manUnit == "Lat")

par(mfrow = c(2, 2))
hist(trimOuDat$realizedHR, breaks = 10)
hist(trimOuDat$finalCatch, breaks = 10)
hist(lateDat$realizedHR, breaks = 10)
hist(lateDat$finalCatch, breaks = 10)



# calculate mean and SD for additive deviations - we can use these if we decide
# to stick w/ log-normal
trimOuDat %>% 
  # group_by(manUnit) %>% 
  summarize(mean(hrDev),
            sd(hrDev))



### Simulations to a) check normal approximation works and 
# b) experiment with beta
N = 1000
xSeq <- seq(from = 0, to = 1, length.out = N)
meanHR = 0.4
obsSample <- sample(trimOuDat$hrDev, size = N, replace = TRUE)

# Normal
sigN <- 0.07
realHRsN <- rnorm(N, meanHR, sigN)
devsNorm <- data.frame(dev = realHRsN - meanHR,
                       dataType = "sim")
obsDat <- data.frame(dev = obsSample,
                     dataType = "obs")
devsN <- rbind(devsNorm, obsDat) %>% 
  mutate(shapeSD = "norm") %>% # add so it can be plotted with beta values
  select(shapeSD, dev, dataType)


# Beta
shapeSD <- seq(from = 0.05, to = 0.15, by = 0.01) #range of sigma values
realHRs <- matrix(NA, nrow = N, ncol = length(shapeSD))
hrDevs <- matrix(NA, nrow = N, ncol = length(shapeSD))

for(i in seq_along(shapeSD)) {
  location <-  meanHR^2 * ((1 - meanHR)/shapeSD[i]^2 - 1/(meanHR))
  scale <- location * (1/meanHR - 1)
  realHRs[ , i] <- dbeta(xSeq, location, scale)
  hrDevs[ , i] <- (rbeta(N, location, scale, ncp = 0) - meanHR)
}
colnames(realHRs) <- shapeSD
colnames(hrDevs) <- shapeSD

## Look at realized productivity curves for different sigma values
dum <- melt(realHRs) %>% 
  dplyr::rename("shapeSD" =  "Var2", "realHR" = "value") %>% 
  mutate(xVals = rep(xSeq, times = 11)) %>%
  dplyr::select(-Var1) 

ggplot(dum, aes(x = xVals, y = realHR, color = as.factor(shapeSD))) +
  geom_line() +
  theme_sleekX()

# This is the proportion of realized harvest rate deviations that exceed 20%
# the maximum observed deviation in the time series
sum(hrDevs[,3] > 0.2, hrDevs[,3] < -0.2) / nrow(hrDevs)


## look at deviations between target and realized
devs <- melt(hrDevs) %>% 
  dplyr::rename("shapeSD" =  "Var2", "dev" = "value") %>% 
  dplyr::select(-Var1) %>% 
  mutate(dataType = "sim")

#sample observed data N times and combine with simulated
obsDat <- data.frame(shapeSD = rep(unique(shapeSD), each = N),
                     dev = rep(obsSample, times = length(shapeSD)),
                     dataType = "obs")

## Combine all 3 dataframes
devsFull <- rbind(devs, obsDat, devsN)
# %>% 
  # filter(dataType == "sim") %>% 
  # filter(shapeSD < 0.08 | shapeSD == "norm")

ggplot(devsFull, aes(x = dev, color = dataType)) +
  geom_histogram(fill="white", alpha=0.5, position="identity")+
  theme_sleekX() +
  facet_wrap(~as.factor(shapeSD))

## Optimal OU SD for beta distribution looks to be 0.07; 0.15 as high treatment?

#Both log-normal and beta distributions do reasonable job recreating observed 
#deviations; for reasonable HRs and sigmas


### Code blurb by Sean Anderson showing impacts of different mu and sigma values
library(manipulate)
x <- seq(0.01, 0.99, length.out = 300)
manipulate({
  alpha <- mu^2 * (((1-mu)/sigma^2) - (1/mu))
  y <- dbeta(
    x = x,
    shape1 = alpha,
    shape2 = alpha * (1/mu - 1))
  plot(x = x, y = y, type = "l", ylim = c(0, max(y)))},
  mu = slider(0.02, 0.98, 0.8),
  sigma = slider(0.05, 1, 0.08))



### Repeat above analysis with stock specific realized HRs to help ID lower and 
## upper bounds
# Beta
shapeSD <- seq(from = 0.05, to = 0.15, by = 0.01) #range of sigma values
realHRs <- matrix(NA, nrow = N, ncol = length(shapeSD))
hrDevs <- matrix(NA, nrow = N, ncol = length(shapeSD))

for(i in seq_along(shapeSD)) {
  location <-  meanHR^2 * ((1 - meanHR)/shapeSD[i]^2 - 1/(meanHR))
  scale <- location * (1/meanHR - 1)
  hrDevs[ , i] <- (rbeta(N, location, scale, ncp = 0) - meanHR)
}
colnames(hrDevs) <- shapeSD

## look at deviations between target and realized
devs <- melt(hrDevs) %>% 
  dplyr::rename("shapeSD" =  "Var2", "dev" = "value") %>% 
  dplyr::select(-Var1) %>% 
  mutate(dataType = "sim")

#sample observed data N times and combine with simulated
muNames <- unique(trimOuDat$manUnit)
obsDat <- NULL
for(i in seq_along(muNames)) {
  dum <- trimOuDat %>% 
    filter(manUnit == muNames[i])
  sampleObs <- sample(dum$hrDev, size = N, replace = TRUE)
  devsByMU <- devs %>% 
    mutate(manUnit = muNames[i])
  dum2 <- data.frame(shapeSD = rep(unique(shapeSD), each = N),
                     dev = rep(sampleObs, times = length(shapeSD)),
                     dataType = "obs",
                     manUnit = muNames[i]) %>% 
    rbind(devsByMU)
  # obsDat <- rbind(obsDat, dum2)
  q <- ggplot(dum2, aes(x = dev, color = dataType)) +
    geom_histogram(fill="white", alpha=0.5, position="identity")+
    theme_sleekX() +
    labs(title = muNames[i]) +
    facet_wrap(~as.factor(shapeSD))
  plot(q)
}


## Doesn't really produce anything inteligible


### Look at variation in realized exploitation rates
muNames <- read.csv(here("data/sox/fraserCUpars.csv"), stringsAsFactors = F) %>% 
  select(stk, manUnit)
catchDat <- read.csv(here("data/sox/fraserCatchDatTrim.csv"), 
                     stringsAsFactors = F) %>% 
  filter(yr > 2008) %>% 
  left_join(., muNames, by = "stk")


ggplot(catchDat, aes(x = as.factor(manUnit), y = totalER)) +
  geom_boxplot() +
  theme_sleekX() +
  facet_wrap(~yr)
  
