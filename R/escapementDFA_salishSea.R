## chinDFA.R
# May 27, 2019
# Script to fit DFAs of ESCAPEMENT data by by streams; initially constrained 
# only to Salish Sea; includes both MARSS and Bayesian DFA versions
# -----

escDat <- read.csv(here("data", "salmonData", "CLEANsalishSea_escData.csv"), 
                   stringsAsFactors = FALSE)

listOfPackages <- c("here", "MARSS", "tidyverse", "ggplot2", "parallel", 
                    "doParallel", "foreach", "tictoc", "bayesdfa")
lapply(listOfPackages, library, character.only = TRUE)

source(here("R/functions/dfaFunctions.R"))


## Time series of escapement
ggplot(escDat, aes(x = year, y = esc)) +
  geom_line() +
  samSim::theme_sleekX() +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_wrap(~stock, scales = "free_y")


## Fit DFA
### Fit DFA
#Convert to matrices (necessary for DFA)
escMat <- escDat %>%
  spread(key = stock, value = esc) %>% 
  select(-year) %>% 
  as.matrix() %>% 
  t()
retYrs <- unique(escDat$year)
colnames(escMat) <- retYrs

nStks <- nrow(escMat)
nYrs <- ncol(escMat)
stkID <- rownames(escMat)

#preliminary model fit test
modelList <- list(m = 2, R = "diagonal and equal")
cntrList <- list(maxit = 200)
dfaTest <- MARSS(escMat, model = modelList, z.score = TRUE, form = "dfa", 
                 control = cntrList, method = "kem")


## Fit models in parallel using up to five trends
inMList <- list(1, 2, 3, 4, 5)
subDir <- "escapementSalishSea"
Ncores <- detectCores()
cl <- makeCluster(Ncores - 4) #save two cores
registerDoParallel(cl)
clusterEvalQ(cl, c(library(MARSS), library(here), library(Rcpp),
                   library(RcppArmadillo)))
clusterExport(cl, c("escMat", "inMList", "fitDFA", "subDir"),
              envir=environment())
tic("run in parallel")
parLapply(cl, inMList, function(x) {
  fitDFA(escMat, inR = "diagonal and unequal", inM = x, maxIteration = 3000, 
         subDirName = subDir)
})
stopCluster(cl) #end cluster
toc()


#Model rankings
summ <- getTopDFA(subDir)
summ[[2]]


#Explore fit of top model
mod1 <- summ[[1]]

estZ <- coef(mod1, type = "matrix")$Z
#retrieve rotated matrix
invH <- if (ncol(estZ) > 1) {
  varimax(estZ)$rotmat
} else if (ncol(estZ) == 1) {
  1
}

# Loadings
rotZ <- (estZ %*% invH) %>% 
  as.data.frame() %>% 
  mutate(stock = stkID) %>% 
  inner_join(escDat %>% select(stock), by = "stock") %>% 
  gather(key = "trend", value = "loading", -stock) %>% 
  distinct() %>% 
  mutate(trend = as.numeric(as.factor(trend)),
         stock = factor(stock, unique(stock))) %>% 
  distinct() 
# %>% 
#   mutate(stock = factor(stock, unique(stock)),
#          #invert T1 and T2 to make more intuitive
#          # loadingT = case_when(
#          #   trend %in% c("1", "2", "4") ~ (loading * -1), 
#          #   TRUE ~ loading
#          # )
#          )

ggplot(rotZ, aes(x = stock, y = loading)) +
  geom_col() +
  samSim::theme_sleekX(axisSize = 9, legendSize = 0.7) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~trend)
#much less coherence in abundance than survival rates

rotTrends <- (solve(invH) %*% mod1$states) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(year = retYrs) %>% 
  gather(key = "trend", value = "est", -year) %>% 
  mutate(trend = as.numeric(as.factor(trend)))

ggplot(rotTrends, aes(x = year, y = est)) +
  geom_line() +
  samSim::theme_sleekX() +
  geom_hline(yintercept = 0, colour = "red") +
  facet_wrap(~trend)


modOutCI <- broom::augment(mod1, interval = "confidence") %>% 
  dplyr::rename(stock = .rownames) %>% 
  inner_join(escDat %>% select(stock), by = "stock") %>% 
  distinct() %>% 
  arrange(stock) %>% 
  mutate(year = t + min(escDat$year), #add so that year is correct 
         stock = factor(stock, unique(stock)))
ggplot(modOutCI) +
  geom_line(aes(x = year, y = .fitted)) +
  geom_point(aes(x = year, y = y)) +
  geom_ribbon(aes(x = year, ymin = .conf.low, ymax = .conf.up), linetype = 2, 
              alpha = 0.2) +
  # ggtitle(grps[x]) +
  samSim::theme_sleekX() +
  facet_wrap(~stock)


## Fit Bayesian DFA to better account for uncertainty
varIn <- list("unequal", "equal")

Ncores <- detectCores()
cl <- makeCluster(Ncores) #save two cores
registerDoParallel(cl)
clusterEvalQ(cl, c(library(bayesdfa), library(here), library(Rcpp),
                   library(RcppArmadillo)))
clusterExport(cl, c("escMat", "find_dfa_trends", "varIn"),
              envir=environment())
tic("run in parallel")
dum <- parLapply(cl, varIn, function(x) {
  find_dfa_trends(y = escMat, iter = 4000, kmin = 1, kmax = 4, chains = 1, 
                  compare_normal = FALSE, zscore = TRUE, variance = x,
                  control = list(adapt_delta = 0.96, max_treedepth = 20))
})
stopCluster(cl) #end cluster
toc()

saveRDS(dum, here::here("data", "dfaBayesFits",
                        "salishSeaEscapement_diffCorStructures_findTrends.rds"))

lapply(dum, function(x) x$summary) %>% 
  do.call(rbind, .) %>% 
  arrange(looic)

modOut <- dum[[1]]$best_model

r1 <- rotate_trends(modOut)
plot_trends(r1)
plot_fitted(modOut)
plot_loadings(r1) +
  ylim(-3, 3)


## Hacky version to alter plot_loadings
# plot_fitted <- function (modelfit, names = NULL) 
# {
stkData <- data.frame(ID = c("Cowichan", "Fr. Spr. 1.2", "Fr. Spr. 1.3", 
                "Fr. Summ. 0.3", "Fr. Summ. 1.3.", "Green", 
                "Harrison", "L. Shuswap", "Lake Wash.", "Nanaimo",
                "Nicola Spr. 1.2", "Nooksack Spr.", "Skagit Spr.", 
                "Skagit Summ.", "Snohomish", "Stillaguamish", "Upp. SoG")) %>% 
  mutate(lifeHistory = case_when(
    ID %in% c("Fr. Spr. 1.2", "Fr. Spr. 1.3", "Fr. Summ. 1.3.",
                   "Nicola Spr. 1.2", "Nooksack Spr.") ~ "yearling",
    TRUE ~ "subyearling/mixed"
  ))
  
n_ts <- dim(modOut$data)[1]
# n_years <- dim(modOut$data)[2]
pred <- predicted(modOut)
df <- data.frame(ID = rep(stkData$ID, n_years), 
                 Time = sort(rep(seq(1975, 2018, by = 1), n_ts)), 
                 mean = c(t(apply(pred, c(3, 4), mean))), 
                 lo = c(t(apply(pred, c(3, 4), quantile, 0.025))), 
                 hi = c(t(apply(pred, c(3, 4), quantile, 0.975))), 
                 y = c(modOut$data)) %>% 
  full_join(stkData, by = "ID") %>% 
  arrange(lifeHistory) %>% 
  # mutate(ID = as.factor(ID)) %>% 
  mutate(ID = factor(ID, unique(ID)))

png(here("figs", "dfa", "salishSeaEsc.png"), height = 5.25, 
    width = 8.5, units = "in", res = 300)
ggplot(df, aes_string(x = "Time", y = "mean")) + 
  geom_ribbon(aes_string(ymin = "lo", ymax = "hi"), alpha = 0.4) + 
  geom_line() + 
  geom_point(aes_string(x = "Time", y = "y", col = "lifeHistory"), 
             size = 1, alpha = 0.6) + 
  facet_wrap("ID", scales = "free_y") +
  labs(x = "Return Year", y = "Survival Anomaly", col = "Life History") +
  samSim::theme_sleekX(axisSize = 8, legendSize = 0.85)
dev.off()
