expcomp <- function(object, ...)
{ 	 		
  fList <- list(b3, b4, b5, l3, l4, l5, w3, w4, baro5)
	fnList <- lapply(fList, function(x) x$name)
  print("Fitting all sigmoidal models...")
  flush.console()
	modList <- lapply(fList, function(x) pcrfit(object$DATA, 1, 2, x))   
	
	EXP <- expfit(object, plot = FALSE, ...)
	expMod <- EXP$mod
	expReg <- EXP$cycles

	RMSEs <- sapply(modList, function(x) RMSE(x, which = expReg))
	RMSEs <- c(RMSEs, EXP$RMSE)  	
	
	modList <- c(modList, list(expMod))
	fnList <- c(fnList, "expGrowth")
  cols <- rk <- rank(RMSEs)
  cols[cols != 1] <- "grey"   
	cols[cols == 1] <- "red"   
	lwds <- rep(1, length(RMSEs))
	lwds[cols != "grey"] <- 3 
  
  print(cols)  	
				
	for (i in 1:length(modList)) {
		if (i == 1) plot(modList[[i]], col = cols[i], subset = c(min(expReg) - 3, max(expReg) + 3),
						            lwd = lwds[i], xlab = "Cycles", ylab = "Raw fluorescence", 
						            main = "Fitting within the exponential region")
		else plot(modList[[i]], add = TRUE, col = cols[i], lwd = lwds[i])
	}
	return(cbind(model = fnList[order(rk)], RMSE = RMSEs[order(rk)]))
}	