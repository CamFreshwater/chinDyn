## chinDFA.R
# April 24, 2019
# Script to fit DFAs by region; output visualized in chinDFAOutput.Rmd
# -----

listOfPackages <- c("here", "MARSS", "tidyverse", "ggplot2", "parallel",
                    "doParallel", "foreach", "tictoc")
lapply(listOfPackages, require, character.only = TRUE)

# byDat <- read.csv(here("data/salmonData/CLEANcwtInd_age2SR_BY.csv"), 
#                   stringsAsFactors = FALSE)
eyDatF <- read.csv(here("data/salmonData/CLEANcwtInd_age2SR_OEY.csv"),
                  stringsAsFactors = FALSE)
# eyDat <- read.csv("C:/github/chinDyn/data/salmonData/CLEANcwtInd_age2SR_OEY.csv", 
#                   stringsAsFactors = FALSE)
source("C:/github/chinDyn/R/functions/dfaFunctions.R")

eyDat <- eyDatF %>% 
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
  group_by(stock) %>% 
  mutate(survZ = as.numeric(scale(surv))) %>%
  ungroup(stock) %>% 
  arrange(grp) %>% 
  mutate(stock = factor(stock, unique(stock))) %>% 
  select(-stockName, -jurisdiction, -lat, -long)


## Initial time series plots
plotSurvZ <- function(group) {
  dum <- eyDat %>% 
    filter(grp == group)
  q <- ggplot(dum, aes(x = OEY, y = survZ, colour = region)) +
    geom_line() +
    samSim::theme_sleekX() +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle(group) +
    facet_wrap(~stock)
  print(q)
}
grpSeq <- unique(eyDat$grp)
sapply(grpSeq, function(x) plotSurvZ(x))

### Fit DFA
#Convert to matrices (necessary for DFA)
eyMatZ <- eyDat %>%
  select(OEY, stock, survZ) %>% 
  spread(key = stock, value = survZ) %>% 
  select(-OEY) %>% 
  as.matrix() %>% 
  t()
entryYrs <- unique(eyDat$OEY)
colnames(eyMatZ) <- entryYrs

nStks <- nrow(eyMatZ)
nYrs <- ncol(eyMatZ)
stkID <- rownames(eyMatZ)

#preliminary model fit test
modelList <- list(m = 2, R = "diagonal and equal")
cntrList <- list(maxit = 200)
dfaTest <- MARSS(eyMatZ, model = modelList, z.score = TRUE, form = "dfa", 
                 control = cntrList, method = "kem")


## Fit models in parallel using multiple cores
inRSeq <- c("diagonal and equal", "diagonal and unequal", "equalvarcov")
inMList <- list(1, 2, 3, 4, 5)
subDir <- "full_OEY"
for (i in seq_along(inRSeq)) {
  Ncores <- detectCores()
  inRDum <- inRSeq[i]
  cl <- makeCluster(Ncores - 4) #save two cores
  registerDoParallel(cl)
  clusterEvalQ(cl, c(library(MARSS), library(here), library(Rcpp),
                     library(RcppArmadillo)))
  clusterExport(cl, c("eyMatZ", "inRDum", "inMList", "fitDFA", "subDir"),
                envir=environment())
  tic("run in parallel")
  parLapply(cl, inMList, function(x) {
    fitDFA(eyMatZ, inR = inRDum, inM = x, maxIteration = 2000, 
           subDirName = subDir)
  })
  stopCluster(cl) #end cluster
  toc()
}

summ <- getTopDFA(subDir)
summ[[2]]
### Difficulty getting even low trend models to converge but explore anyways

## Explore fit of top model
# mod1 <- summ[[1]]

mod1 <- readRDS(here::here("data", "dfaFits", subDir, 
                           "fitMod.2.equalvarcov.rds"))

estZ <- coef(mod1, type = "matrix")$Z
#retrieve rotated matrix
invH <- if (ncol(estZ) > 1) {
  varimax(estZ)$rotmat
} else if (ncol(estZ) == 1) {
  1
}

rotZ <- rotateLoadings(zIn = estZ, H = invH, stkNames = stkID, 
                       survDat = eyDat) %>% 
  inner_join(eyDat %>% select(stock, grp), 
             by = "stock") %>% 
  distinct() %>% 
  arrange(grp) %>% 
  mutate(stock = factor(stock, unique(stock)))

ggplot(rotZ, aes(x = stock, y = loading, fill = grp)) +
  geom_col() +
  samSim::theme_sleekX(axisSize = 9, legendSize = 0.7) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~trend)

rotTrends <- rotateTrends(modIn = mod1, H = invH)
ggplot(rotTrends, aes(x = year, y = est)) +
  geom_line() +
  samSim::theme_sleekX() +
  geom_hline(yintercept = 0, colour = "red") +
  facet_wrap(~trend)

modOutCI <- broom::augment(mod1, interval = "confidence") %>% 
  dplyr::rename(stock = .rownames) %>% 
  inner_join(eyDat %>% select(stock, region, grp), by = "stock") %>% 
  distinct() %>% 
  arrange(grp) %>% 
  mutate(year = t + 1984, #add so that year is correct 
         stock = factor(stock, unique(stock)))

ggplot(modOutCI) +
  geom_line(aes(x = year, y = .fitted)) +
  geom_point(aes(x = year, y = y, colour = grp)) +
  geom_ribbon(aes(x = year, ymin = .conf.low, ymax = .conf.up), linetype = 2, 
              alpha = 0.2) +
  facet_wrap(~stock)


