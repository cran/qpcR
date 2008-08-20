mchoice <- function (object, fctList = NULL, sig.level = 0.05, verbose = TRUE, crit = c("ftest", "ratio", "weights"))
{
      crit <- match.arg(crit)
      mtype <- paste(qpcR:::typeid(object), "()", sep = "")
      
      if (crit == "ftest" && !is.null(fctList)) stop("Use Akaike weights for a list of non-nested models!")

      if (is.null(fctList)) {
            if (mtype == "b3()" || mtype == "b4()" || mtype == "b5()") {
		    fctList <- list(b3(), b4(), b5())
            }
      	if (mtype == "l3()" || mtype == "l4()" || mtype == "l5()") {
		    fctList <- list(l3(), l4(), l5())
		}
		if (mtype == "w4()" || mtype == "w3()") {
		    fctList <- list(w3(), w4())
		}
      }
      
      lenFL <- length(fctList)
	retMat <- matrix(nrow = lenFL, ncol = 6)
	rn <- NULL	

	for (i in 1:lenFL) {
		newObj <- update(object, fct = fctList[[i]])
		typeObj <- paste(qpcR:::typeid(newObj), "()", sep = "")  
		rn[i] <- typeObj
		retMat[i, 1] <- round(logLik(newObj), 2)
    		retMat[i, 2] <- round(AIC(newObj), 2)
    		retMat[i, 3] <- round(AICc(newObj), 2)
    		retMat[i, 4] <- round(resVar(newObj), 5)
    		if (i > 1 && crit == "ftest") retMat[i, 5] <- anova(newObj, prevObj, details = FALSE)[2, 5]
			else retMat[i, 5] <- NA
            if (i > 1) retMat[i, 6] <- LR(newObj, prevObj)$p.value else retMat[i, 6] <- NA
		prevObj <- newObj
	}
	
      aic.w <- round(akaike.weights(retMat[, 2])$weights, 3)
      aicc.w <- round(akaike.weights(retMat[, 3])$weights, 3)
      
      retMat <- cbind(retMat, aic.w, aicc.w)
	
	colnames(retMat) <- c("logLik", "AIC", "AICc", "resVar", "ftest", "LR.pval", "AIC.weights", "AICc.weights")
	rownames(retMat) <- rn

	if (verbose) print(retMat)
	
      if (crit == "ftest") {
            sigTRUE <- retMat[, 5] < sig.level
            sigTRUE[is.na(sigTRUE)] <- FALSE
            WHICH <- which(sigTRUE == TRUE)
            if (length(WHICH) > 0) SELECT <- max(WHICH) else SELECT <- 1
      }
      
      if (crit == "ratio") {
            sigTRUE <- retMat[, 6] < sig.level
            sigTRUE[is.na(sigTRUE)] <- FALSE
            WHICH <- which(sigTRUE == TRUE)
            if (length(WHICH) > 0) SELECT <- max(WHICH) else SELECT <- 1
      }

      if (crit == "weights") {
            SELECT <- which(retMat[, 7] == max(retMat[, 7], na.rm = TRUE))
      }
      
      fmodel <- eval(as.call(list(update, object, fct = fctList[[SELECT]])))

	fmodel$retMat <- retMat
    	return(fmodel)
}
