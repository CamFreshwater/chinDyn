## chinDFA.R
# April 24, 2019
# Script to fit DFAs by region; output visualized in chinDFAOutput.Rmd
# -----

listOfPackages <- c("here", "MARSS", "tidyverse", "ggplot2", "parallel", 
                    "doParallel", "foreach", "tictoc")
lapply(listOfPackages, library, character.only = TRUE)

eyDatFull <- read.csv(here("data/salmonData/CLEANcwtInd_age2SR_OEY.csv"), 
                  stringsAsFactors = FALSE)

source(here("R/functions/dfaFunctions.R"))

eyDat <- eyDatFull %>% 
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
eyMat <- eyDat %>%
  select(OEY, stock, surv) %>% 
  spread(key = stock, value = surv) %>% 
  select(-OEY) %>% 
  as.matrix() %>% 
  t()
entryYrs <- unique(eyDat$OEY)
colnames(eyMat) <- entryYrs

nStks <- nrow(eyMat)
nYrs <- ncol(eyMat)
stkID <- rownames(eyMat)

#preliminary model fit test
modelList <- list(m = 2, R = "diagonal and equal")
cntrList <- list(maxit = 200)
dfaTest <- MARSS(eyMat, model = modelList, z.score = TRUE, form = "dfa", 
                 control = cntrList, method = "kem")


## Fit models in parallel using multiple cores
# inRSeq <- c("diagonal and equal", "diagonal and unequal", "equalvarcov")
inMList <- list(2, 3, 4, 5)
subDir <- "full_OEY"
# for (i in seq_along(inRSeq)) {
  Ncores <- detectCores()
  inRDum <- "diagonal and unequal"#inRSeq[i]
  cl <- makeCluster(Ncores - 4) #save two cores
  registerDoParallel(cl)
  clusterEvalQ(cl, c(library(MARSS), library(here), library(Rcpp),
                     library(RcppArmadillo)))
  clusterExport(cl, c("eyMat", "inRDum", "inMList", "fitDFA", "subDir"),
                envir=environment())
  tic("run in parallel")
  parLapply(cl, inMList, function(x) {
    fitDFA(eyMat, inR = inRDum, inM = x, maxIteration = 3000, 
           subDirName = subDir)
  })
  stopCluster(cl) #end cluster
  toc()
}

summ <- getTopDFA(subDir)
summ[[2]]

## Explore fit of top model
mod1 <- summ[[1]]

mod1 <- readRDS(here::here("data", "dfaFits", subDir, 
                           "fitMod.4.diagonalAndUnequal.rds"))

estZ <- coef(mod1, type = "matrix")$Z
#retrieve rotated matrix
invH <- if (ncol(estZ) > 1) {
  varimax(estZ)$rotmat
} else if (ncol(estZ) == 1) {
  1
}

# Loadings
rotZ <- rotateLoadings(zIn = estZ, H = invH, stkNames = stkID, 
                       survDat = eyDat) %>% 
  inner_join(eyDat %>% select(stock, grp), 
             by = "stock") %>% 
  distinct() %>% 
  arrange(grp, stock) %>% 
  mutate(stock = factor(stock, unique(stock)),
         #invert T1 and T2 to make more intuitive
         loadingT = case_when(
           trend %in% c("1", "2", "4") ~ (loading * -1), 
           TRUE ~ loading
         ))

png(here("figs", "dfa", "globalMARSS", "4Trends_loadings.png"), height = 5.5, 
    width = 9.5, units = "in", res = 300)
ggplot(rotZ, aes(x = stock, y = loadingT, fill = grp)) +
  geom_col() +
  samSim::theme_sleekX(axisSize = 9, legendSize = 0.7) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~trend)
dev.off()

# Trends
rotTrends <- rotateTrends(modIn = mod1, H = invH) %>% 
  mutate(estT = case_when(
    trend %in% c("1", "2", "4") ~ (est * -1), 
    TRUE ~ est
  ))
  
png(here("figs", "dfa", "globalMARSS", "est4Trends.png"), height = 6, 
    width = 6, units = "in", res = 300)
ggplot(rotTrends, aes(x = year, y = estT)) +
  geom_line() +
  samSim::theme_sleekX() +
  geom_hline(yintercept = 0, colour = "red") +
  facet_wrap(~trend)
dev.off()

#Fits
modOutCI <- broom::augment(mod1, interval = "confidence") %>% 
  dplyr::rename(stock = .rownames) %>% 
  inner_join(eyDat %>% select(stock, region, grp), by = "stock") %>% 
  distinct() %>% 
  arrange(grp, stock) %>% 
  mutate(year = t + 1984, #add so that year is correct 
         stock = factor(stock, unique(stock)))

grps <- unique(modOutCI$grp)
sapply(seq_along(grps), function(x) {
  dum <- modOutCI %>% 
    filter(grp == grps[x])
  p <- ggplot(dum) +
    geom_line(aes(x = year, y = .fitted)) +
    geom_point(aes(x = year, y = y)) +
    geom_ribbon(aes(x = year, ymin = .conf.low, ymax = .conf.up), linetype = 2, 
                alpha = 0.2) +
    ggtitle(grps[x]) +
    samSim::theme_sleekX() +
    facet_wrap(~stock)
  
  fileName <- paste(grps[x], "4trends_fits.png", sep = "_")
  png(here("figs", "dfa", "globalMARSS", fileName), height = 5, 
      width = 6, units = "in", res = 300)
  print(p)
  dev.off()
})


### Repeat above but with seals as a covariate
sealDat <- read.csv(here("data/salmonData/survCovariateAnom.csv"), 
                stringsAsFactors = FALSE) %>% 
  select(OEY = year, sealAnom)

ssDat <- eyDat %>%
  filter(grp %in% c("oceantype_SS", "streamtype_SS")) %>% 
  full_join(sealDat, by = "OEY") %>% 
  filter(!is.na(sealAnom),
         !OEY < 1972) %>% 
  arrange(OEY)
ssMat <- ssDat %>% 
  select(OEY, stock, surv) %>% 
  spread(key = stock, value = surv) %>% 
  select(-OEY) %>% 
  as.matrix() %>% 
  t()
entryYrs <- unique(ssDat$OEY)
colnames(ssMat) <- entryYrs

nStks <- nrow(ssMat)
nYrs <- ncol(ssMat)
stkID <- rownames(ssMat)

sealMat <- ssDat %>%
  filter(stock == "SPS") %>% 
  select(sealAnom) %>% 
  as.matrix() %>%
  t()
colnames(sealMat) <- entryYrs
  

#preliminary model fit test
modelList <- list(m = 3, R = "diagonal and equal")
cntrList <- list(maxit = 2500)
fitMod <- MARSS(ssMat, model = modelList, z.score = TRUE, form = "dfa", 
                 control = cntrList, method = "kem", covariates = sealMat)
while(fitMod$convergence != 0) {
  fitMod <- MARSS(ssMat, model = modelList, z.score = TRUE,
                  control = list(maxit = 8000),
                  inits = as.matrix(coef(fitMod)[1]), form = "dfa", 
                  z.score  = TRUE, method = "BFGS-kf", covariates = sealMat)
}

saveRDS(fitMod, here::here("data", "dfaFits", "twoTrend_seal_fits.Rds"))
