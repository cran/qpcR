pcrGOF <- function(object)
{
return(list(Rsq = Rsq(object), AICc = AICc(object), AIC = AIC(object), resVar = resVar(object), RMSE = RMSE(object)))
}
