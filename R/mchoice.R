mchoice <- function (object, sig.level = 0.05, verbose = TRUE) 
{
	mtype <- paste(qpcR:::typeid(object), "()", sep = "")

	if (mtype == "b3()" || mtype == "b4()" || mtype == "b5()") {
		fctList <- list(b3(), b4(), b5())
	}

	if (mtype == "l3()" || mtype == "l4()" || mtype == "l5()") {
		fctList <- list(l3(), l4(), l5())
	}
	
	lenFL <- length(fctList)
	retMat <- matrix(rep(NA, 5), nrow = lenFL, ncol = 5, byrow = TRUE)
	rn <- NULL	

	for (i in 1:lenFL) {
		newObj <- update(object, fct = fctList[[i]])
		typeObj <- paste(qpcR:::typeid(newObj), "()", sep = "")  
		rn[i] <- typeObj
		retMat[i, 1] <- logLik(newObj)
    		retMat[i, 2] <- AIC(newObj)
    		retMat[i, 3] <- AICc(newObj)
    		retMat[i, 4] <- summary(newObj)$resVar
    		if (i > 1) retMat[i, 5] <- anova(newObj, prevObj, details = FALSE)[2, 5] 
			else retMat[i, 5] <- NA
		prevObj <- newObj
	}
	
	colnames(retMat) <- c("logLik", "AIC", "AICc", "Res var", "nested F-test")
	rownames(retMat) <- rn

	if (verbose) print(retMat)
	
	if (mtype == "l3()" || mtype == "l4()" || mtype == "l5()") {
		if (retMat[3, 5] < sig.level) {
          		fmodel <- update(object, fct = l5())
        	}
        	else {
            	if (retMat[2, 5] < sig.level) {
          			fmodel <- update(object, fct = l4())
            	}
            	else {
                	fmodel <- update(object, fct = l3())
            	}
        	}
    	}

	if (mtype == "b3()" || mtype == "b4()" || mtype == "b5()") {
        	if (retMat[3, 5] < sig.level) {
     			fmodel <- update(object, fct = b5())
        	}
        	else {
            	if (retMat[2, 5] < sig.level) {
                	fmodel <- update(object, fct = b4())
            	}
            	else {
                	fmodel <- update(object, fct = b3())
            	}
        	}
    	}

	fmodel$retMat <- retMat
    	return(fmodel)
}
