## linearModelTrends.R
# May 22, 2019
# Script to combine covariate data with modeling information. 
# -----
cov <- read.csv(here("data/salmonData/survCovariateAnom.csv"), 
                stringsAsFactors = FALSE)
rotList <- readRDS(here::here("data", "dfaBayesFits", 
                          "coastWide_estTrends_topModels.rds"))
eyDat <- read.csv(here("data/salmonData/CLEANcwtInd_age2SR_OEY.csv"), 
                  stringsAsFactors = FALSE) 

yrSeq <- unique(eyDat$OEY)

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
  fileName <- paste(unique(plotD$group), "covScatter.png")
  p <- ggplot(plotD, aes(x = anomaly, y = meanT1)) +
    geom_point() +
    samSim::theme_sleekX() +
    labs(y = "DFA Trend 1",
         title = unique(plotD$group)) +
    facet_wrap(~var)
  
  png(here("figs", "dfa", "covEffects", fileName), height = 5,
      width = 5, units = "in", res = 300)
  print(p)
  dev.off()
})
