pcrGOF <- function(object, PRESS = FALSE, ...)
{
  retList <- list(Rsq = Rsq(object), Rsq.ad = Rsq.ad(object), AIC = AIC(object), 
            AICc = AICc(object), BIC = BIC(object), resVar = resVar(object), 
            RMSE = RMSE(object))    
  
  p.neill <- try(neill.test(object), silent = TRUE)
  if (!inherits(p.neill, "try-error")) retList <- c(retList, p.neill = p.neill)           
            
  if (PRESS) {
    P.square <- PRESS(object, verbose = FALSE)$P.square
    retList <- c(retList, P.square = P.square)
  }
  return(retList) 
}
