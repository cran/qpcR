mselect <- function(
object,
fctList = NULL,
sig.level = 0.05,
verbose = TRUE,
crit = c("ftest", "ratio", "weights", "chisq"),
do.all = FALSE, 
...
)
{
  crit <- match.arg(crit)
  if (any(class(object) == "replist")) object <- object[[1]]
  else if (any(class(object) == "pcrfit")) object <- object 
  else stop("'object' must be either of class 'pcrfit' or 'replist'!")
  
  mtype <- object$MODEL$name

  if (is.null(fctList)) {
    if (mtype %in% c("b3", "b4", "b5", "b6")) fctList <- list(b3, b4, b5, b6)
    if (mtype %in% c("l3", "l4", "l5", "l6")) fctList <- list(l3, l4, l5, l6)
    if (mtype == "w4" || mtype == "w3") fctList <- list(w4, w3)
  }
  
  if (do.all) {
      fctList <- list(l6, l5, l4, l3, b6, b5, b4, b3, w4, w3, baro5)
      crit <- "weights"
  }
  
  retMat <- matrix(nrow = length(fctList), ncol = 7)
  rn <- NULL	

  ml <- lapply(fctList, function(x) pcrfit(object$DATA, 1, 2, x, opt.method = object$opt.method))

  for (i in 1:length(ml)) {
    rn[i] <- ml[[i]]$MODEL$name
    retMat[i, 1] <- round(logLik(ml[[i]]), 2)
		retMat[i, 2] <- round(AIC(ml[[i]]), 2)
		retMat[i, 3] <- round(AICc(ml[[i]]), 2)
		retMat[i, 4] <- round(resVar(ml[[i]]), 5)   		
    if (i < length(ml)) retMat[i + 1 , 5] <- as.matrix(anova(ml[[i]], ml[[i + 1]]))[2, 6]
    if (i < length(ml)) retMat[i + 1, 6] <- LR(ml[[i]], ml[[i + 1]])$p.value  
    retMat[i, 7] <- fitchisq(ml[[i]], ...)$chi2.red            
  }           
	
  aic.w <- round(akaike.weights(retMat[, 2])$weights, 3)
  aicc.w <- round(akaike.weights(retMat[, 3])$weights, 3)       
  retMat <- cbind(retMat, aic.w, aicc.w)      
  	
  colnames(retMat) <- c("logLik", "AIC", "AICc", "resVar", "ftest", "LR", "Chisq", "AIC.weights", "AICc.weights")
  rownames(retMat) <- rn
	
  if (verbose) print(retMat)
	
  if (crit == "ftest") {
    modTRUE <- retMat[, 5] < sig.level     
    if(all(is.na(modTRUE))) stop("nested f-test was unsuccessful! Probably not nested (df = 0)?")   
    modTRUE[is.na(modTRUE)] <- FALSE     
    WHICH <- which(modTRUE)                
    SELECT <- max(WHICH)
    if (any(modTRUE == TRUE)) optModel <- fctList[[SELECT]] else optModel <- object$MODEL      
  }
      
  if (crit == "ratio") {     
    if (any(retMat[, 6] == 0, na.rm = TRUE)) stop("likelihood ratio p-value is 0! Probably not nested (df = 0)?")     
    modTRUE <- retMat[, 6] < sig.level   
    modTRUE[is.na(modTRUE)] <- FALSE      
    WHICH <- which(modTRUE)     
    SELECT <- max(WHICH) 
    optModel <- fctList[[SELECT]]   
  }

  if (crit == "weights") {
    SELECT <- which.max(retMat[, 9])
    optModel <- fctList[[SELECT]]
  }
  
  if (crit == "chisq") {
    SELECT <- which.min(retMat[, 7])
    optModel <- fctList[[SELECT]]
  }
      
  optMod <- pcrfit(object$DATA, 1, 2, optModel, opt.method = object$opt.method)    

  optMod$retMat <- retMat
  return(optMod)
}
