## threeTrendBayesDFA_allRegions
# April 24, 2019
# Script to fit three trend DFAs by region; similar to chinBayesDFA.Rmd
# -----

listOfPackages <- c("here", "bayesdfa", "tidyverse", "ggplot2", "parallel", 
                    "doParallel", "foreach", "tictoc")
lapply(listOfPackages, library, character.only = TRUE)

#helper functions to fit and post-process DFA
source(here("R/functions/dfaFunctions.R"))

eyDat <- read.csv(here("data/salmonData/CLEANcwtInd_age2SR_OEY.csv"), 
                  stringsAsFactors = FALSE)

eyDatFull <- eyDat %>% 
  mutate(lat = as.numeric(lat),
         long = as.numeric(long),
         aggReg = case_when(
           (is.na(lat)) ~ "north",
           (lat > 52 & !region == "UFR") ~ "north",
           (region %in% c("JFUCA", "LCOLR", "MCOLR", "ORCST", "UCOLR", "WACST",
                          "WCVI")) ~ "south",
           TRUE ~ "SS"
         )) %>% 
  mutate(grp = paste(smoltType, aggReg, sep = "_")) %>% 
  filter(!grp %in% c("oceantype_north", "streamtype_south")) %>% 
  arrange(desc(lat)) %>% 
  mutate(stock = factor(stock, unique(stock)),
         grp = factor(grp, unique(grp))) %>% 
  select(-stockName, -jurisdiction, -lat, -long)  %>% 
  mutate(grp = fct_recode(grp, "Yearling\nNorth" = "streamtype_north", 
                          "Yearling\nSalish Sea" = "streamtype_SS", 
                          "Subyearling\nSalish Sea" = "oceantype_SS",
                          "Subyearling\nSouth" = "oceantype_south"))

#split into lists 
convMat <- function(mat) {
  mat %>% 
    select(OEY, stock, surv) %>%
    spread(key = stock, value = surv) %>%
    select(-OEY) %>% 
    as.matrix() %>% 
    t()
}
survList <- split(eyDatFull, eyDatFull$grp) %>% 
  lapply(., convMat)

listStkNames <- lapply(survList, function(x) {
  name <- rownames(x)
  numSeq <- seq(1, length(name), by = 1)
  out <- vector(mode = "character", length = length(numSeq))
  setNames(name, numSeq)
})

options(mc.cores = parallel::detectCores())
Ncores <- detectCores()
cl <- makeCluster(Ncores - 4) #save four cores
registerDoParallel(cl)
clusterEvalQ(cl, c(library(bayesdfa), library(here), library(Rcpp),
                   library(RcppArmadillo), library(dplyr)))
clusterExport(cl, c("fit_dfa", "survList"), envir=environment())
tic("run in parallel")
dum <- parLapply(cl, survList, function(x) {
  fit_dfa(y = x, num_trends = 3, zscore = TRUE, iter = 4000, chains = 4, 
          thin = 1, control = list(adapt_delta = 0.97, max_treedepth = 20),
          estimate_nu = TRUE)
  # find_dfa_trends(y = x, kmin = 3, kmax = 3, zscore = TRUE, iter = 4000,
  #                 chains = 4, control = list(adapt_delta = 0.97,
  #                                            max_treedepth = 20),
  #                 compare_normal = TRUE, variance =  "unequal")
})
stopCluster(cl) #end cluster
toc()

saveRDS(dum, here::here("data", "dfaBayesFits", "coastWide_fitThreeTrends.rds"))


## Group trends
rotList <- lapply(dum, function(x) {
  # mod <- dum[[x]]$best_model
  rotate_trends(x)
})

names(rotList) <- names(survList)
# saveRDS(rotList, here::here("data", "dfaBayesFits",
#                             "coastWide_estTrends_topModels.rds"))

trendList <- lapply(seq_along(rotList), function(x) {
  p <- plot_trendsX(rotList[[x]], oneTrend = TRUE, startYr = 1972) +
    labs(title = names(survList)[x], x = "Ocean Entry Year", 
         y = "Survival Anomaly") +
    samSim::theme_sleekX()
  return(p)
})
png(here("figs", "dfa", "rangeWideBayes", "3Trends_EqualCov",
         "dominantCoastTrends.png"), height = 6, 
    width = 6, units = "in", res = 300)
ggpubr::ggarrange(trendList[[1]], trendList[[2]], trendList[[3]], 
                  trendList[[4]],  ncol = 2, nrow = 2)
dev.off()

fitList <- lapply(seq_along(dum), function(x) {
  mod <- dum[[x]]
  p <- plot_fittedX(mod, startYr = 1972) +
    ggtitle(names(survList)[x]) +
    facet_wrap("ID", labeller = as_labeller(listStkNames[[x]]), 
               scales = "free_y") +
    samSim::theme_sleekX()
  return(p)
})
for(x in seq_along(fitList)) {
  fileName <- paste(abbreviate(names(survList)[x], 6), "3TrendFits.png", 
                    sep = "_")
  png(here("figs", "dfa", "rangeWideBayes", "3Trends_EqualCov", fileName), 
      height = 5, 
      width = 6, units = "in", res = 300)
  plot(fitList[[x]])
  dev.off()
}
