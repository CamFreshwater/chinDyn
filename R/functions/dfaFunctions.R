## dfaFunctions.R
# May 3, 2019
# Companion functions to chinDFA.R
# -----

## Function to fit DFA and save output data
fitDFA <- function(mat, inR, inM, maxIteration = 500, subDirName) {
  cntrList <- list(minit = 100, maxit = maxIteration)
  dfaModel <- list(A = "zero", R = inR, m = inM)
  fitMod <- MARSS(mat, model = dfaModel, control = cntrList,
                  form = "dfa", z.score  = TRUE)
  #If faster algorithm fails to converge use alternative
  while(fitMod$convergence != 0) {
    fitMod <- MARSS(mat, model = dfaModel, control = list(maxit = 8000),
                    inits = as.matrix(coef(fitMod)[1]), form = "dfa", 
                    z.score  = TRUE, method = "BFGS-kf")
  }
  outName <- paste("fitMod", inM, rapportools::tocamel(inR), "rds", sep=".")
  
  #save data 
  # dirPath <- paste("C:/github/chinDyn/data/dfaFits/", subDirName, sep = "")
  # dir.create(dirPath, recursive = TRUE, showWarnings = FALSE)
  # saveRDS(fitMod, paste(dirPath, outName, sep = "/"))
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


##Function to rotate factor loadings
rotateLoadings <- function(zIn, H, stkNames = stkID, survDat) {
  (zIn %*% H) %>% 
    as.data.frame() %>% 
    mutate(stock = stkNames) %>% 
    inner_join(survDat %>% select(stock, region), by = "stock") %>% 
    gather(key = "trend", value = "loading", -stock, -region) %>% 
    distinct() %>% 
    arrange(region) %>% 
    mutate(trend = as.numeric(as.factor(trend)),
           stock = factor(stock, unique(stock)))
}


##Function to rotate estimated trends
rotateTrends <- function(modIn, H){
  (solve(H) %*% modIn$states) %>% 
    t() %>% 
    as.data.frame() %>% 
    mutate(year = entryYrs) %>% 
    gather(key = "trend", value = "est", -year) %>% 
    mutate(trend = as.numeric(as.factor(trend)))
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


##Modifications to plot trends
plot_trendsX <- function(rotated_modelfit, years = NULL, 
                          highlight_outliers = FALSE, threshold = 0.01,
                          oneTrend = FALSE, startYr = 1972) {
  rotated <- rotated_modelfit
  n_ts <- dim(rotated$Z_rot)[2]
  n_trends <- dim(rotated$Z_rot)[3]
  n_years <- dim(rotated$trends_mean)[2]
  if (is.null(years)) 
    years <- seq_len(n_years)
  df <- data.frame(x = c(t(rotated$trends_mean)), 
                   lo = c(t(rotated$trends_lower)), 
                   hi = c(t(rotated$trends_upper)), 
                   trend = paste0("Trend ", 
                                  sort(rep(seq_len(n_trends), n_years))), 
                   time = rep(years, n_trends))
  if (oneTrend == TRUE) {
    dum <- df %>% 
      dplyr::filter(trend == "Trend 1") %>% 
      dplyr::mutate(timeT = (startYr - 1 + time))
    p1 <- ggplot(dum, 
                 aes_string(x = "timeT", y = "x")) + 
      geom_ribbon(aes_string(ymin = "lo", ymax = "hi"), alpha = 0.4) + 
      geom_line() + 
      xlab("Year") + ylab("")
  } else {
    p1 <- ggplot(df, aes_string(x = "time", y = "x")) + 
      geom_ribbon(aes_string(ymin = "lo", ymax = "hi"), 
                  alpha = 0.4) + geom_line() + facet_wrap("trend") + 
      xlab("Time") + ylab("")
  }                                                                                                                   
  if (highlight_outliers) {
    swans <- find_swans(rotated, threshold = threshold)
    df$outliers <- swans$below_threshold
    p1 <- p1 + geom_point(data = df[which(df$outliers), ], 
                          color = "red")
  }
  p1
}


plot_fittedX <- function (modelfit, names = NULL, startYr = 1972) 
{
  n_ts <- dim(modelfit$data)[1]
  n_years <- dim(modelfit$data)[2]
  pred <- predicted(modelfit)
  df <- data.frame(ID = rep(seq_len(n_ts), n_years), 
                   Time = sort(rep(seq_len(n_years), n_ts)), 
                   mean = c(t(apply(pred, c(3, 4), mean))), 
                   lo = c(t(apply(pred, c(3, 4), quantile, 0.025))), 
                   hi = c(t(apply(pred, c(3,  4), quantile, 0.975))), 
                   y = c(modelfit$data)) %>% 
    dplyr::mutate(timeT = (startYr - 1 + Time))
  if (!is.null(names)) {
    df$ID <- names[df$ID]
  }
  p1 <- ggplot(df, aes_string(x = "timeT", y = "mean")) + 
    geom_ribbon(aes_string(ymin = "lo",  ymax = "hi"), alpha = 0.4) + 
    geom_line() + geom_point(aes_string(x = "timeT", y = "y"), col = "red", 
                             size = 0.5, alpha = 0.4) + 
    facet_wrap("ID",  scales = "free_y") + xlab("Ocean Entry Year") + 
    ylab("Survival Anomaly")
  p1
}
