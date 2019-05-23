## chinDFA.R
# May 22, 2019
# Script to fit two trend Bayesian DFAs to each stock group based on model 
# selection results in chinBayesDFA.Rmd
# -----

listOfPackages <- c("here", "bayesdfa", "tidyverse", "ggplot2", "parallel", 
                    "doParallel", "foreach", "tictoc")
lapply(listOfPackages, library, character.only = TRUE)

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
  arrange(grp) %>% 
  mutate(stock = factor(stock, unique(stock))) %>% 
  filter(!grp %in% c("oceantype_north", "streamtype_south")) %>% 
  select(-stockName, -jurisdiction, -lat, -long)

#Split groupe data into lists of matrices
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

#Surprisingly verbose chunk of code to create named character list for labelling
#subsequent plots
listStkNames <- lapply(survList, function(x) {
  name <- rownames(x)
  numSeq <- seq(1, length(name), by = 1)
  out <- vector(mode = "character", length = length(numSeq))
  setNames(name, numSeq)
})

#Fit models
#set up trend inputs to allow for parallel processing
# options(mc.cores = parallel::detectCores())
# Ncores <- detectCores()
# cl <- makeCluster(Ncores - 4) #save four cores
# registerDoParallel(cl)
# clusterEvalQ(cl, c(library(bayesdfa), library(here), library(Rcpp),
#                    library(RcppArmadillo), library(dplyr)))
# clusterExport(cl, c("fit_dfa", "survList"), envir=environment())
# tic("run in parallel")
# dum <- parLapply(cl, survList, function(x) {
#   find_dfa_trends(y = x, kmin = 2, kmax = 2, zscore = TRUE, iter = 4000,
#                   chains = 4, control = list(adapt_delta = 0.97, 
#                                              max_treedepth = 20),
#                   compare_normal = TRUE, variance =  c("equal", "unequal"))
# })
# stopCluster(cl) #end cluster
# toc()
# 
# saveRDS(dum, here::here("data", "dfaBayesFits", "coastWide_fitTwoTrends.rds"))


dum <- readRDS(here::here("data", "dfaBayesFits", "coastWide_fitTwoTrends.rds"))


# Check model fits
lapply(seq_along(survList), function(x) 
  dum[[x]]$summary %>% 
    arrange(looic) %>% 
    mutate(group = names(survList)[x])
) %>% 
  do.call(rbind, .)


# Check trends
trendsByGroups <- function(x) {
  mod <- x$best_model
  plot_trends(mod)
}
rotList <- lapply(seq_along(dum), function(x) {
  mod <- dum[[x]]$best_model
  rotate_trends(mod)
})

names(rotList) <- names(survList)
saveRDS(rotList, here::here("data", "dfaBayesFits",
                            "coastWide_estTrends_topModels_V2.rds"))

trendList <- lapply(seq_along(rotList), function(x) {
  p <- plot_trends(rotList[[x]]) +
    ggtitle(names(survList)[x])
  return(p)
})
png(here("figs", "dfa", "rangeWideBayes", "2Trends_MixedCov",
         "coastTrends.png"), height = 8, 
    width = 6, units = "in", res = 300)
ggpubr::ggarrange(trendList[[1]], trendList[[2]], trendList[[3]], 
                  trendList[[4]],  ncol = 1, nrow = 4)
dev.off()

loadList <- lapply(seq_along(rotList), function(x) {
  p <- plot_loadings(rotList[[x]]) +
    ggtitle(names(survList)[x]) +
    ylim(-2, 2)
  return(p)
})

png(here("figs", "dfa", "rangeWideBayes", "2Trends_MixedCov",
         "coastLoadings.png"), height = 9, 
    width = 6, units = "in", res = 300)
ggpubr::ggarrange(loadList[[1]], loadList[[2]], loadList[[3]], 
                  loadList[[4]],  ncol = 1, nrow = 4, common.legend = TRUE)
dev.off()

# Plot fits
fitList <- lapply(seq_along(dum), function(x) {
  mod <- dum[[x]]$best_model
  p <- plot_fitted(mod) +
    ggtitle(names(survList)[x]) +
    facet_wrap("ID", labeller = as_labeller(listStkNames[[x]]), 
               scales = "free_y")
  return(p)
})
for(x in seq_along(fitList)) {
  fileName <- paste(names(survList)[x], "2TrendFits.png", sep = "_")
  png(here("figs", "dfa", "rangeWideBayes", "2Trends_MixedCov", fileName), 
      height = 5, width = 6, units = "in", res = 300)
  plot(fitList[[x]])
  dev.off()
}
