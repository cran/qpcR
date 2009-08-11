pcrGOF <- function(object)
{
return(list(Rsq = Rsq(object), Rsq.ad = Rsq.ad(object), AIC = AIC(object), AICc = AICc(object), BIC = BIC(object), resVar = resVar(object), RMSE = RMSE(object)))
}
