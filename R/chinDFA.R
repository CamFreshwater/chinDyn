## chinDFA.R
# April 24, 2019
# Script to fit DFAs; output visualized in chinDFAOutput.Rmd
# -----

listOfPackages <- c("here", "MARSS", "tidyverse", "ggplot2", "parallel", 
                    "doParallel", "foreach", "tictoc")
lapply(listOfPackages, require, character.only = TRUE)

byDat <- read.csv(here("data/salmonData/CLEANcwtInd_age2SR_BY.csv"), 
                  stringsAsFactors = FALSE)
eyDat <- read.csv(here("data/salmonData/CLEANcwtInd_age2SR_OEY.csv"), 
                  stringsAsFactors = FALSE)

#focus on subset of BC pops for initial analyses
byDatTrim <- byDat %>% 
  filter(region %in% c("LFR", "MFR", "UFR", "ECVI", "SPGSD", "NPGSD"))%>% 
  group_by(stock) %>% 
  mutate(survZ = as.numeric(scale(surv))) %>%
  ungroup(stock) %>% 
  arrange(region) %>% 
  mutate(stock = factor(stock, unique(stock))) %>% 
  select(-stockName, -jurisdiction, -lat, -long)

## Initial time series plots
png(here("figs", "standardizedSurvival_salishSeaOnly.png"), height = 7,
    width = 7, units = "in", res = 300)
ggplot(byDatTrim, aes(x = BY, y = survZ, colour = region)) +
  geom_line() +
  theme_sleekX() +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_wrap(~stock)
dev.off()

byDatLH <- byDatTrim %>% 
  arrange(smoltType) %>% 
  mutate(stock = factor(stock, unique(stock)))
ggplot(byDatLH, aes(x = BY, y = survZ, colour = smoltType)) +
  geom_line() +
  samSim::theme_sleekX() +
  facet_wrap(~stock)
#no clear patterns by smolt type

byDatRT <- byDatTrim %>% 
  arrange(adultRunTiming) %>% 
  mutate(stock = factor(stock, unique(stock)))
ggplot(byDatRT, aes(x = BY, y = survZ, colour = adultRunTiming)) +
  geom_line() +
  samSim::theme_sleekX() +
  facet_wrap(~stock)
#no clear patterns by run timing


### Fit DFA
#Convert to matrices (necessary for DFA)
byMatZ <- byDatTrim %>%
  select(BY, stock, survZ) %>% 
  spread(key = stock, value = survZ) %>% 
  select(-BY) %>% 
  as.matrix() %>% 
  t()
broodYrs <- unique(byDatTrim$BY)
colnames(byMatZ) <- broodYrs

nStks <- nrow(byMatZ)
nYrs <- ncol(byMatZ)
stkID <- rownames(byMatZ)

#preliminary model fit test
modelList <- list(m = 2, R = "diagonal and equal")
cntrList <- list(maxit = 200)
dfaTest <- MARSS(byMatZ, model = modelList, z.score = TRUE, form = "dfa", 
                 control = cntrList, method = "BFGS-kf")


## Fit models in parallel using multiple cores
inRSeq <- c("diagonal and equal", "diagonal and unequal", "equalvarcov")
inMList <- list(2, 3, 4, 5, 6, 7)
subDir <- "salishSeaOnly"
for (i in seq_along(inRSeq)) {
  Ncores <- detectCores()
  inRDum <- inRSeq[i]
  cl <- makeCluster(Ncores - 2) #save two cores
  registerDoParallel(cl)
  clusterEvalQ(cl, c(library(MARSS), library(here), library(Rcpp),
                     library(RcppArmadillo)))
  clusterExport(cl, c("byMatZ", "inRDum", "inMList", "fitDFA", "subDir"),
                envir=environment())
  tic("run in parallel")
  parLapply(cl, inMList, function(x) {
    fitDFA(byMatZ, inR = inRDum, inM = x, maxIteration = 1000, 
           subDirname = subDir)
  })
  stopCluster(cl) #end cluster
  toc()
}


summ <- getTopDFA(subDir)
summ[[2]]
## Strong support for using diagonal and unequal covariance matrix and 
# intermediate number of trends

## Explore fit of top model
mod1 <- summ[[1]]

estZ <- coef(mod1, type = "matrix")$Z
#retrieve rotated matrix
invH <- if (ncol(estZ) > 1) {
  varimax(estZ)$rotmat
} else if (ncol(estZ) == 1) {
  1
}

## Factor loadings
#rotate factor loadings
rotZ <- (estZ %*% invH) %>% 
  as.data.frame() %>% 
  mutate(stock = stkID) %>% 
  inner_join(byDatTrim %>% select(stock, region), by = "stock") %>% 
  gather(key = "trend", value = "loading", -stock, -region) %>% 
  distinct() %>% 
  arrange(region) %>% 
  mutate(trend = fct_recode(as.factor(trend),
                            trend1 = "V1",
                            trend2 = "V2",
                            trend3 = "V3",
                            trend4 = "V4"),
         stock = factor(stock, unique(stock)))

png(here("figs", "dfa", paste(subDir, ncol(estZ), "Trend", "_loadings.png", 
                              sep = "")), height = 4, width = 5.5, units = "in",
    res = 300)
ggplot(rotZ, aes(x = stock, y = loading, fill = region)) +
  geom_col() +
  theme_sleekX(axisSize = 9, legendSize = 0.7) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~trend)
dev.off()


## Trends
#rotate trends
rotTrends <- (solve(invH) %*% mod1$states) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(year = broodYrs) %>% 
  gather(key = "trend", value = "est", -year) %>% 
  mutate(trend = fct_recode(as.factor(trend),
                            trend1 = "V1",
                            trend2 = "V2",
                            trend3 = "V3",
                            trend4 = "V4"))

png(here("figs", "dfa", paste(subDir, ncol(estZ), "Trend", "_trends.png", 
                              sep = "")), height = 4, width = 4, units = "in",
    res = 300)
ggplot(rotTrends, aes(x = year, y = est)) +
  geom_line() +
  theme_sleekX() +
  geom_hline(yintercept = 0, colour = "red") +
  facet_wrap(~trend)
dev.off()

## Plot fitted estimates
modOutCI <- broom::augment(mod1, interval = "confidence") %>% 
  dplyr::rename(stock = .rownames) %>% 
  inner_join(byDatTrim %>% select(stock, region), by = "stock") %>% 
  distinct() %>% 
  mutate(year = t + 1970, 
         stock = as.factor(stock)) %>% 
  arrange(region) 

ggplot(modOutCI) +
  geom_line(aes(x = year, y = .fitted)) +
  geom_point(aes(x = year, y = y, colour = region)) +
  geom_ribbon(aes(x = year, ymin = .conf.low, ymax = .conf.up), linetype = 2, 
              alpha = 0.2) +
  facet_wrap(~stock)



dum <- get_DFA_fits(mod1)
namesIn <- names(dum)
dum2 <- lapply(seq_along(dum), function (x) gatherList(dum[[x]], yrs = broodYrs, 
                                                       stks = stkID, 
                                                       names = namesIn[x])) 
modOutCI2 <- do.call(rbind, dum2) %>% 
  spread(key = estimate, value = value) %>% 
  inner_join(byDatTrim %>% select(stock, region), by = "stock") %>% 
  distinct() %>% 
  mutate(stock = as.factor(stock)) %>% 
  arrange(region) %>% 
  filter(year > 1985)
  
ggplot(modOutCI2) +
  geom_line(aes(x = year, y = ex)) +
  geom_point(aes(x = modOutCI$year, y = modOutCI$y, colour = modOutCI$region)) +
  geom_ribbon(aes(x = year, ymin = lo, ymax = up), linetype = 2, 
              alpha = 0.2) +
  facet_wrap(~stock)


