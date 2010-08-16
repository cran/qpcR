baseline <- function(
object,
refit = TRUE,
refit.model = NULL, 
verbose = TRUE,
...
)
{
  CLASS <- class(object)
  if (!(CLASS[1] %in% c("pcrfit", "modlist"))) stop("Object is not of class 'pcrfit'!")
  if (CLASS[1] == "pcrfit") object <- list(object) 
  NAMES <- lapply(object, function(x) x$names) 
  outLIST <- list()    
  
  for (i in 1:length(object)) {
    tempObject <- object[[i]]
    if (!(tempObject$MODEL$name %in% c("l6", "b6"))) stop("All objects must be fitted with 'l6' or 'b6' model")
    if (verbose) cat("Baselining ", if (refit) "and refitting " else "", NAMES[[i]], " ...\n", sep = "")
    flush.console()
    COEF <- coef(tempObject)
    X <- tempObject$DATA[, 1]
    Y <- tempObject$DATA[, 2]
    if (tempObject$MODEL$name == "l6") modY <- Y - COEF[2] - COEF[6] * log(X)
    if (tempObject$MODEL$name == "b6") modY <- Y - COEF[2] - COEF[6] * X
    newDATA <- cbind(Cycles = X, Fluo = modY)
    if (refit) {
      parList <- list(data = newDATA, cyc = 1, fluo = 2, model = if (is.null(refit.model)) tempObject$MODEL else refit.model,
                      start = COEF, ...)
      newObject <- do.call(pcrfit, parList)
    }
    newObject$DATA.base <- newDATA
    outLIST[[i]] <- newObject
  }
    
  cycLIST <- lapply(outLIST, function(x) x$DATA[, 1])
  CYC <- unique(matrix(do.call(cbind.na, cycLIST), ncol = 1))
  CYC <- CYC[!is.na(CYC)]
  fluoLIST <- lapply(outLIST, function(x) x$DATA[, 2])
  FLUO <- do.call(data.frame.na, fluoLIST)
  colnames(FLUO) <- NAMES 
  baseDATA <- cbind.na(Cycles = CYC, FLUO)  
  
  if (CLASS == "pcrfit") OUT <- outLIST[[1]] else OUT <- outLIST  
  class(OUT) <- CLASS
  attr(OUT, "baselined") <- baseDATA  
  return(OUT)
}       
