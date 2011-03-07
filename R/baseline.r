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
  outLIST <- vector("list", length = length(object))   
  
  for (i in 1:length(object)) {
    tempObject <- object[[i]]
    if (!(tempObject$MODEL$name %in% c("l6", "b6"))) stop("All objects must be fitted with 'l6' or 'b6' model")
    if (verbose) cat("Baselining ", if (refit) "and refitting " else "", NAMES[[i]], " ...\n", sep = "")
    flush.console()
    COEF <- coef(tempObject)
    X <- tempObject$DATA[, 1]
    Y <- tempObject$DATA[, 2]
    modY <- Y - COEF[2] - COEF[6] * X
    newDATA <- cbind(Cycles = X, Fluo = modY)
    if (refit) {
      parList <- list(data = newDATA, cyc = 1, fluo = 2, model = if (is.null(refit.model)) tempObject$MODEL else refit.model, 
                      nls.method = "default", start = COEF, ...)
      newObject <- do.call(pcrfit, parList)
    }
    newObject$DATA.base <- newDATA
    outLIST[[i]] <- newObject
  }
    
  if (CLASS == "pcrfit") OUT <- outLIST[[1]] else OUT <- outLIST  
  class(OUT) <- CLASS     
  return(OUT)
}       
