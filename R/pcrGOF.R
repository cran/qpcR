pcrGOF <- function(object, PRESS = FALSE)
{
  retList <- list(Rsq = Rsq(object), Rsq.ad = Rsq.ad(object), AIC = AIC(object), 
            AICc = AICc(object), BIC = BIC(object), HQIC = HQIC(object), 
            resVar = resVar(object), RMSE = RMSE(object))    
  
  p.neill <- try(neill.test(object), silent = TRUE)
  if (!inherits(p.neill, "try-error") & !is.na(p.neill)) retList <- c(retList, p.neill = p.neill) 
  
  fcsq <- try(fitchisq(object), silent = TRUE)
  if (!inherits(fcsq, "try-error") & !is.na(fcsq)) retList <- c(retList, chi2.red = fcsq$chi2.red)
            
  if (PRESS) {
    P.square <- tryCatch(PRESS(object, verbose = TRUE)$P.square, error = function(e) NA)
    retList <- c(retList, P.square = P.square)
  }
  
  return(retList) 
}
