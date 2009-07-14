pcropt1 <- function (object, fact = 3, opt = FALSE, plot = TRUE, ...) 
{
  
  window.l <- NULL
  window.u <- NULL    
  aic <- NULL
  aicc <- NULL
  resvar <- NULL
  eff <- NULL
  init.exp <- NULL
  init.sig <- NULL  
  
  START <- try(efficiency(object, plot = FALSE))
  if (inherits(START, "try-error")) stop("Could not initialize optimization. Try different 'fact'!")
  
  cpD1 <- round(START$cpD1)
  cpD2 <- round(START$cpD2)
  lower <- cpD1 - fact * (cpD1 - cpD2)
  upper <- cpD1 + fact * (cpD1 - cpD2)    
  lowerseq <- 1:lower
  upperseq <- nrow(object$DATA):upper    
    
  for (i in lowerseq) {
    for (j in upperseq) {
      newData <- object$DATA[i:j, 1:2]
     
      newCurve <- try(pcrfit(newData, 1, 2, model = object$MODEL), silent = TRUE)
      if (inherits(newCurve, "try-error")) break      
      
      if (opt) newCurve <- mselect(newCurve, verbose = FALSE, ...)       
            
      ITER <- try(efficiency(newCurve), silent = TRUE)   
      if (inherits(ITER, "try-error")) break
      
      window.l <- c(window.l, i)
      window.u <- c(window.u, j)        
      aic <- c(aic, ITER$AIC)
      aicc <- c(aicc, ITER$AICc)
      resvar <- c(resvar, ITER$resVar)
      eff <- c(eff, ITER$eff)
      init.sig <- c(init.sig, ITER$init1)
      init.exp <- c(init.exp, ITER$init2)      
      }
    }
    
    resMat <- cbind(as.numeric(window.l), as.numeric(window.u),  
                    as.numeric(aic), as.numeric(aicc), as.numeric(resvar), 
                    as.numeric(eff), as.numeric(init.exp), as.numeric(init.sig))
                    
    if (plot) {
      par(mfrow = c(2, 2))
      boxplot(resMat[, 4], main = "AICc")
      boxplot(resMat[, 6], main = "Eff")                 
      boxplot(resMat[, 7], main = "init.exp")       
      boxplot(resMat[, 8], main = "init.sig")      
    }         
    
    colnames(resMat) <- c("lower", "upper", "AIC", "AICc", "resVar", "Eff", "init.exp", "init.sig")      
    return(as.data.frame(resMat))
}
