expcomp <- function(object, fitcyc = 5)
{
	obj <- object
	modList <- list()
		
	fList <- list(b3(), b4(), b5(), l3(), l4(), l5(), 
				LL2.3(), LL2.4(), LL2.5(), baro5(),
				w3(), w4(), W2.3(), W2.4())
	fnList <- c("b3", "b4", "b5", "l3", "l4", "l5",
				"l3.2", "l4.2", "l5.2", "baro5", 
				"w3", "w4", "w3.2", "w4.2", "exp")
	
	for (i in 1:length(fList)) {
		newmod <- update(obj, fct = fList[[i]])
		modList[[i]] <- newmod
	}

	EXP <- expfit2(object, fitcyc = fitcyc)
	expMod <- EXP$mod
	expReg <- EXP$outlier:(EXP$outlier + fitcyc - 1)

	rmses <- sapply(modList, function(x) RMSE(x, which = expReg))
	rmses <- c(rmses, EXP$RMSE)
	modList[[length(modList) + 1]] <- expMod
	cols <- rk <- rank(rmses)
	cols[cols == 1] <- "red"
	cols[cols == 2] <- "orange"
	cols[cols == 3] <- "yellow"
	cols[cols != "red" & cols != "orange" & cols != "yellow"] <- "grey"
	lwds <- rep(1, length(rmses))
	lwds[cols != "grey"] <- 3
				
	for (i in 1:length(modList)) {
		if (i == 1) pcrplot(modList[[i]], col = cols[i], xlim = c(min(expReg) - 1, max(expReg) + 1),
						lwd = lwds[i], xlab = "Cycles", ylab = "Raw fluorescence", 
						main = "Fitting within the exponential region")
		else pcrplot(modList[[i]], type = "none", add = TRUE, col = cols[i], lwd = lwds[i])
	}
	return(cbind(model = fnList[order(rk)], RMSE = rmses[order(rk)]))
}	