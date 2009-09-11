plot.pcrfit <- function(
x, 
fitted = TRUE,
subset = NULL,
confband = c("none", "confidence", "prediction"),
errbar = c("none", "sd", "se", "conf"),
add = FALSE, 
colvec = NULL,
level = 0.95,    
...) 
{
  object <- x
  confband <- match.arg(confband)
  errbar <- match.arg(errbar)
  if (!is.null(subset) && length(subset) != 2) stop("'subset' must be a 2-element vector!") 
  
  if (class(object) != "modlist") modList <- list(object) else modList <- object
  if (class(object) == "modlist") {
    minVal <- min(sapply(modList, function(x) x$DATA[, 2]), na.rm = TRUE)
    maxVal <- max(sapply(modList, function(x) x$DATA[, 2]), na.rm = TRUE)
  } else {
    minVal <- NULL
    maxVal <- NULL    
  }  
  
  if (is.null(colvec)) {
    colvec <- rep(1, length(modList))     
    if (!is.na(class(object)[2]) && class(object)[2] == "replist") colvec <- gl(attr(object, "nlevels"), 1)
  }      
  
  statfun <- function(object, errbar, level) {
    CYC <- object$DATA[, 1]     
    DATA <- object$DATA[, 2]
    fact <- qnorm(1 - ((1 - level)/2))       
    nobs <- round(length(DATA)/length(unique(CYC)))    
    statfun <- switch(errbar, sd = function(x) sd(x, na.rm = TRUE), 
                              se = function(x) sd(x, na.rm = TRUE)/sqrt(nobs), 
                              conf = function(x) sd(x, na.rm = TRUE)/sqrt(nobs) * fact)     
    STAT <- tapply(DATA, as.factor(CYC), function(x) statfun(x))         
  }
  
  for (i in 1:length(modList)) {   
    if (class(modList[[i]])[1] != "pcrfit") stop("object must be of class 'pcrfit'")       
   
    if (i > 1) add <- TRUE
    
    if (!add) plot(modList[[i]]$DATA, ylim = c(minVal, maxVal), xlim = c(subset[1], subset[2]), col = colvec[i], ...) 
        else points(modList[[i]]$DATA, col = colvec[i], ...)
    
    if (is.null(modList[[i]]$isReps)) {
      if (fitted) {
        CYCS <- modList[[i]]$DATA[, 1] 
        FITTED <- fitted(modList[[i]])
        lines(CYCS, FITTED, col = colvec[i])
      }
    } else {       
      LEN <- length(unique(modList[[i]]$DATA[, 1])) 
      CYCS <- modList[[i]]$DATA[1:LEN, 1]
      FITTED <- fitted(modList[[i]])[1:LEN]                 
      if (fitted) lines(CYCS, FITTED, col = colvec[i])
    }
  
    if (confband != "none") {
      CYC <- unique(modList[[i]]$DATA[, 1])          
      CONFINT <- predict(modList[[i]], interval = confband, level = level, ...)[CYC, ]             
      lines(CYC, CONFINT$Lower, col = 2, ...)
      lines(CYC, CONFINT$Upper, col = 2, ...)  
    }
    
    if (errbar != "none") {
      STAT <- statfun(modList[[i]], errbar, level)
      arrows(CYCS, FITTED + STAT, CYCS, FITTED - STAT, angle = 90, code = 3, col = colvec[i], length = 0.05)      
    }   
      
  }                
}  