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
getDfaAic <- function(subDirName) {
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
  modOut %>% 
    arrange(AICc)
}

getDfaAic(subDir)
## Strong support for using diagonal and unequal covariance matrix and 
# intermediate number of trends