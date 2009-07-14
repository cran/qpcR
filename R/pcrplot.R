pcrplot <- function(
object, 
fitted = TRUE,
confband = c("none", "confidence", "prediction"),
add = FALSE, 
colvec = NULL,    
...) 
{
  confband <- match.arg(confband)
  if (class(object) != "modlist") modList <- list(object) else modList <- object
  if (class(object) == "modlist") {
    minVal <- min(sapply(modList, function(x) x$DATA[, 2]), na.rm = TRUE)
    maxVal <- max(sapply(modList, function(x) x$DATA[, 2]), na.rm = TRUE)
  } else {
    minVal <- NULL
    maxVal <- NULL
  }  
  
  if (is.null(colvec)) colvec <- rep(1, length(modList))     
  
  for (i in 1:length(modList)) {   
    if (class(modList[[i]])[1] != "pcrfit") stop("object must be of class 'pcrfit'")       
   
    if (i > 1) add <- TRUE
    
    if (!add) plot(modList[[i]]$DATA, ylim = c(minVal, maxVal), col = colvec[i], ...) else points(modList[[i]]$DATA, col = colvec[i], ...)
    
    if (fitted) lines(modList[[i]]$DATA[, 1], fitted(modList[[i]]), col = colvec[i], ...)
  
    if (confband != "none") {
      CONFINT <- pcrpred(modList[[i]], interval = confband, ...)
      lines(modList[[i]]$DATA[, 1], CONFINT$Lower, col = 2, ...)
      lines(modList[[i]]$DATA[, 1], CONFINT$Upper, col = 2, ...)  
    }  
  }                
} 

 

