pcrGOF <- function(object, PRESS = FALSE)
{
  retList <- list(Rsq = Rsq(object), Rsq.ad = Rsq.ad(object), AIC = AIC(object), 
            AICc = AICc(object), BIC = BIC(object), resVar = resVar(object), 
            RMSE = RMSE(object), prob = fitprob(object))
  if (PRESS) {
    P.square <- PRESS(object, verbose = FALSE)$P.square
    retList <- c(retList, P.square = P.square)
  }
  return(retList) 
}
