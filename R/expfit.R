expfit <- function (object, fitcyc = 5, plot = TRUE, crit = "resVar", start = 5, maxeff = 2, mineff = 1.5) 
{
    crit <- match.arg(crit, c("resVar", "AIC"))
    rv <- vector()
    Eff <- vector()
    aic <- vector()
    rmse <- vector()
    reseff <- efficiency(object, plot = FALSE)
    cpD1 <- reseff$cpD1
    cpD2 <- reseff$cpD2
    effgood <- vector()
    
    for (i in start:cpD1) {
        expMod <- update(object, data = object$data[i:(i + fitcyc - 
            1), ], fct = expGrowth())
        coeftemp <- as.numeric(coef(expMod)[2])
   	   rv[i] <- round(resVar(expMod), 8)
        aic[i] <- round(AIC(expMod), 2)
        rmse[i] <- round(RMSE(expMod), 5)
        Eff[i] <- round(exp(coeftemp), 3)
        
        if (Eff[i] < mineff || Eff[i] > maxeff) {
		  COL <- 2 
		  effgood[i] <- F
	   }
	   else {
		  COL <- 3
		  effgood[i] <- T
        }	
          
	   mtext <- paste("resVar:", rv[i], "\n", "AIC:", aic[i], 
            "\n", "Eff:", Eff[i], "\n")

        if (plot) 
            pcrplot(object, main = mtext, cex.main = 0.9, xlab = "Cycles", 
                ylab = "raw fluorescence")
        if (plot) 
            pcrplot(expMod, add = TRUE, col = COL, lwd = 2)
    }
            
    if (crit == "resVar") {
	   rv[effgood != TRUE] <- NA
        cyc.best <- which.min(rv)
        cyc.all <- cyc.best:(cyc.best + fitcyc - 1)
    }
    if (crit == "AIC") {
        aic[effgood != TRUE] <- NA
        cyc.best <- which.min(aic)
        cyc.all <- cyc.best:(cyc.best + fitcyc - 1)
    }
    
    EFF.curve1 <- object$data[cyc.all[-1], 2]
    EFF.curve2 <- object$data[cyc.all[-length(cyc.all)], 2]
    EFF.curve <- EFF.curve1/EFF.curve2
    expMod <- update(object, data = object$data[cyc.best:(cyc.best + 
        fitcyc - 1), ], fct = expGrowth())

    if (plot) {
        pcrplot(object)
        pcrplot(expMod, add = TRUE, col = 3, lwd = 2)
    }
    return(list(cyc.best = cyc.all, Eff.fit = exp(as.numeric(coef(expMod)[2])), 
        Eff.curve = EFF.curve, resVar = rv[cyc.best], AIC = aic[cyc.best], 
	   RMSE = rmse[cyc.best], init = as.numeric(coef(expMod)[1]),
        mod = expMod))
}
