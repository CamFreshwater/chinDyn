## linearModelTrends.R
# May 22, 2019
# Script to combine covariate data with modeling information. 
# -----

library(tidyverse)

cov <- read.csv(here::here("data/salmonData/survCovariateAnom.csv"), 
                stringsAsFactors = FALSE)
rotList <- readRDS(here::here("data", "dfaBayesFits", 
                          "coastWide_estTrends_topModels.rds"))
eyDat <- read.csv(here::here("data/salmonData/CLEANcwtInd_age2SR_OEY.csv"), 
                  stringsAsFactors = FALSE) 

yrSeq <- unique(eyDat$OEY)

## Looks at covariates
cov_long <- cov %>% 
  pivot_longer(cols = sealAnom:cciAnom, names_to = "metric", values_to = "anomaly")
cov_ts <- cov_long %>% 
  ggplot(.) +
  geom_line(aes(x = year, y = anomaly)) +
  ggsidekick::theme_sleek() +
  scale_x_continuous(breaks = c(1975, 1985, 1995, 2005, 2015)) +
  facet_wrap(~metric)

png(here::here("figs", "dfa", "covEffects", "cov_anomalies.png"), height = 5,
    width = 7, units = "in", res = 300)
cov_ts
dev.off()

## Prep trends data
pullTrends <- lapply(seq_along(rotList), function(h) {
  dum <- rotList[[h]]
  data.frame(group = names(rotList)[h],
             year = yrSeq,
             meanT1 = dum$trends_mean[1, ],
             meanT2 = dum$trends_mean[2, ]) %>% 
    full_join(., cov, by = "year") %>% 
    filter(!is.na(meanT1)) 
})

## Plot relationship
sapply(pullTrends, function(x) {
  plotD <- x %>% 
    gather(key = "var", value = "anomaly", -group, -year, -meanT1, -meanT2)
  fileName <- paste(abbreviate(unique(plotD$group), minlength = 5), 
                    "covScatter.pdf", sep = "_")
  p <- ggplot(plotD, aes(x = anomaly, y = meanT1)) +
    geom_point() +
    samSim::theme_sleekX() +
    labs(y = "DFA Trend 1",
         title = unique(plotD$group)) +
    facet_wrap(~var)
  p2 <- ggplot(plotD, aes(x = anomaly, y = meanT2)) +
    geom_point() +
    samSim::theme_sleekX() +
    labs(y = "DFA Trend 2",
         title = unique(plotD$group)) +
    facet_wrap(~var)
  
  pdf(here::here("figs", "dfa", "covEffects", fileName), height = 5,
      width = 7)
  print(p)
  print(p2)
  dev.off()
})

