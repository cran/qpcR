expfit2 <- function (object, fitcyc = 5, plot = TRUE, pval = 0.01, ...) 
{
    cyc <- object$data[, 1]
    F <- object$data[, 2]
    which <- outlier(cyc, F, pval = pval)
    cyc.all <- which:(which + fitcyc - 1)
    if (plot) 
        pcrplot(object, xlab = "Cycles", ylab = "raw fluorescence")
    expMod <- multdrc(F[cyc.all] ~ cyc[cyc.all], fct = expGrowth(), 
        ...)
    rv <- resVar(expMod)
    aic <- AIC(expMod)
    rmse <- RMSE(expMod)
    EFF.curve1 <- object$data[cyc.all[-1], 2]
    EFF.curve2 <- object$data[cyc.all[-length(cyc.all)], 2]
    EFF.curve <- EFF.curve1/EFF.curve2
    if (plot) 
        pcrplot(expMod, add = TRUE, col = 3, lwd = 2)
    return(list(outlier = which, Eff.fit = as.numeric(exp(coef(expMod)[2])), 
        Eff.curve = EFF.curve, resVar = rv, AIC = aic, 
        RMSE = rmse, init = as.numeric(coef(expMod)[1]), mod = expMod))
}
