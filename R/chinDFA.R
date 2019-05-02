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
  filter(region %in% c("LFR", "MFR", "UFR", "ECVI", "SPGSD", "NPGSD"))
  
#Convert to matrices (necessary for DFA)
byMat <- byDatTrim %>%
  select(BY, stock, surv) %>% 
  spread(key = stock, value = surv) %>% 
  select(-BY) %>% 
  as.matrix() %>% 
  t()
broodYrs <- unique(byDatTrim$BY)
colnames(byMat) <- broodYrs

nStks <- nrow(byMat)
nYrs <- ncol(byMat)
stkID <- rownames(byMat)

#standardize
byMatZ <- MARSS::zscore(byMat)

#preliminary model fit test
modelList <- list(m = 2, R = "diagonal and equal")
cntrList <- list(maxit = 200)
dfaTest <- MARSS(byMatZ, model = modelList, z.score = TRUE, form = "dfa", 
                 control = cntrList, method = "BFGS-kf")


## Function to fit DFA and save output data
fitDFA <- function(mat, inR, inM, maxIteration = 500, subDirName) {
  cntrList <- list(minit = 100, maxit = maxIteration)
  dfaModel <- list(A = "zero", R = inR, m = inM)
  fitMod <- MARSS(byMatZ, model = dfaModel, control = cntrList,
                 form = "dfa", z.score  = TRUE)
  #If faster algorithm fails to converge use alternative
  while(fitMod$convergence != 0) {
    fitMod <- MARSS(byMatZ, model = dfaModel, control = list(maxit = 4000),
                    inits = as.matrix(coef(fitMod)[1]), form = "dfa", 
                    z.score  = TRUE, method = "BFGS-kf")
  }
  outName <- paste("fitMod", inM, rapportools::tocamel(inR), "rds", sep=".")
  
  #save data 
  dir.create(here::here("data", "dfaFits", subDirName),
             recursive = TRUE, showWarnings = FALSE)
  saveRDS(fitMod, here::here("data", "dfaFits", subDirName, outName))
}


## Fit models in parallel using multiple cores
inRSeq <- c("diagonal and equal", "diagonal and unequal", "equalvarcov",
            "unconstrained")
inMList <- list(2, 3, 4, 5, 6, 7)
subDir <- "salishSeaOnly"
# for (i in seq_along(inRSeq)) {
#   Ncores <- detectCores()
#   inRDum <- inRSeq[i]
#   cl <- makeCluster(Ncores - 2) #save two cores
#   registerDoParallel(cl)
#   clusterEvalQ(cl, c(library(MARSS), library(here), library(Rcpp), 
#                      library(RcppArmadillo)))
#   clusterExport(cl, c("byMatZ", "inRDum", "inMList", "fitDFA"), 
#                 envir=environment())
#   tic("run in parallel")
#   parLapply(cl, inMList, function(x) {
    # fitDFA(byMatZ, inR = inRDum, inM = x, maxIteration = 1000, 
    #        subDirname = subDir)
#   })
#   stopCluster(cl) #end cluster
#   toc()
# }


## Function to screen directory and pull files
getTopDFA <- function(subDirName) {
  dirPath <- here::here("data", "dfaFits", subDirName)
  modOutNames <- list.files(dirPath, pattern="\\.rds$")
  
  modOut <- data.frame(model = rep(NA, times = length(modOutNames)),
                       AICc = NA
                       )
  for(i in 1:length(modOutNames)){ #make list of lists!
    dum <- readRDS(paste(dirPath, modOutNames[i], sep="/"))
    modOut[i, "model"] <- modOutNames[i]
    modOut[i, "AICc"] <- ifelse(is.null(dum$AICc), NA, dum$AICc)
  }
  aicTable <- modOut %>% 
    arrange(AICc)
  topModelName <- aicTable %>% 
    filter(AICc == min(AICc, na.rm = TRUE)) %>% 
    select(model) %>% 
    as.character()
  topModel <- readRDS(paste(dirPath, topModelName, sep = "/"))
  return(list(topModel, aicTable))
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


### Generate some figures 
require(samSim)
ggplot(byDatTrim, aes(x = BY, y = surv, colour = stock)) +
  geom_line() +
  theme_sleekX() +
  facet_wrap(~region)

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
  arrange(region)
  
ggplot(modOutCI2) +
  geom_line(aes(x = year, y = ex)) +
  geom_point(aes(x = modOutCI$year, y = modOutCI$y, colour = modOutCI$region)) +
  geom_ribbon(aes(x = year, ymin = lo, ymax = up), linetype = 2, 
              alpha = 0.2) +
  facet_wrap(~stock)


#Function to scrape data from get_DFA_fits
gatherList <- function(x, yrs, stks, names) {
  rownames(x) <- stks
  tt <- x %>% 
    t() %>% 
    as.data.frame %>% 
    mutate(year = yrs, 
           estimate = names) %>% 
    gather(key = "stock", value = "value", -year, -estimate)
}


get_DFA_fits <- function(MLEobj,alpha=0.05) {
  ## empty list for results
  fits <- list()
  ## extra stuff for var() calcs
  Ey <- MARSS:::MARSShatyt(MLEobj)
  ## model params
  mod_par <- coef(MLEobj, type="matrix")
  ZZ <- mod_par$Z
  ## number of obs ts
  nn <- dim(Ey$ytT)[1]
  ## number of time steps
  TT <- dim(Ey$ytT)[2]
  ## get the inverse of the rotation matrix
  H_inv <- varimax(ZZ)$rotmat
  ## model expectation
  fits$ex <- ZZ %*% H_inv %*% MLEobj$states + matrix(mod_par$A,nn,TT)
  ## Var in model fits
  VtT <- MARSSkfss(MLEobj)$VtT
  VV <- NULL
  for(tt in 1:TT) {
    RZVZ <- mod_par$R - ZZ%*%VtT[,,tt]%*%t(ZZ)
    SS <- Ey$yxtT[,,tt] - Ey$ytT[,tt,drop=FALSE] %*% t(MLEobj$states[,tt,drop=FALSE])
    VV <- cbind(VV,diag(RZVZ + SS%*%t(ZZ) + ZZ%*%t(SS)))    
  }
  SE <- sqrt(VV)
  ## upper (1-alpha)% CI
  fits$up <- qnorm(1-alpha/2)*SE + fits$ex    
  ## lower (1-alpha)% CI
  fits$lo <- qnorm(alpha/2)*SE + fits$ex
  return(fits)
}
