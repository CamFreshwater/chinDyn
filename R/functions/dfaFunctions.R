## dfaFunctions.R
# May 3, 2019
# Companion functions to chinDFA.R
# -----

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


## Function to screen directory containing DFA outputs and pull files
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


##Alternative function to broom::augment() to estimate uncertainty 
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


##Function to scrape data from get_DFA_fits
gatherList <- function(x, yrs, stks, names) {
  rownames(x) <- stks
  tt <- x %>% 
    t() %>% 
    as.data.frame %>% 
    mutate(year = yrs, 
           estimate = names) %>% 
    gather(key = "stock", value = "value", -year, -estimate)
}
