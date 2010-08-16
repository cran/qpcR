KOD <- function(
object,    
method = c("bar", "pam"),
efftype = c("sliwin", "sigfit", "expfit"),
train = TRUE,
remove = FALSE,
alpha = 0.05,
verbose = TRUE, 
...
)
{
  require(cluster, quietly = TRUE)
  method <- match.arg(method)
  efftype <- match.arg(efftype)
  CLASS <- class(object)
  tempList <- GROUP <- NULL
  if (class(object)[1] != "modlist") stop("Please supply either a 'modlist' or 'replist'!")   
   
  ### KOD as defined in Bar et al.
  baretal <- function(x) {
    x <- x[complete.cases(x)]
    if (!train) stat <- sapply(1:length(x), function(y) (x[y] - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE))
    else stat <- sapply(1:length(x), function(y) (x[y] - mean(x[-y], na.rm = TRUE))/sd(x[-y], na.rm = TRUE))
    pval <- 2 * (1 - pnorm(abs(stat)))
    which(pval < alpha)
  }
  
  ### PAM outlier detection
  pamout <- function(x, ...) {
    PAM <- pam(x, k = 2, stand = TRUE, ...)
    TABLE <- table(PAM$clustering)
    ID <- which(TABLE == 1)
    if (length(ID > 0)) OUTL <- which(PAM$clustering == ID) else OUTL <- NULL     
  }                       
    
  if (class(object)[2] == "replist") ITER <- length(object) else ITER <- 1       

  for (i in 1:ITER) {
    if (class(object)[2] == "replist") tempObj <- object[[i]]$modlist else tempObj <- object
    flush.console()
    NAMES <- sapply(tempObj, function(x) x$names) 
    if (verbose) cat("Calculating efficiencies...\n")
    eff <- sapply(tempObj, function(x) switch(efftype, sliwin = sliwin(x, plot = FALSE, ...)$eff,
                                                    sigfit = efficiency(x, plot = FALSE, ...)$eff,
                                                    expfit = expfit(x, plot = FALSE, ...)$eff))    

    if (verbose) cat("Calculating outlier(s)...\n")
    FAIL <- switch(method, bar = baretal(eff), pam = pamout(eff))    
    if (length(FAIL) == 0) PASS <- 1:length(tempObj) else PASS <- (1:length(tempObj))[-FAIL]
    
    for (i in FAIL) tempObj[[i]]$outlier <- TRUE
    for (j in PASS) tempObj[[j]]$outlier <- FALSE 
    
    if (length(FAIL) > 0) {
      if (verbose) cat(" Found kinetic outlier for", NAMES[FAIL], "\n")  
      flush.console() 
      if (remove) {
        if (verbose) cat(" Removing", NAMES[FAIL], "...\n")
        flush.console()
        tempObj <- tempObj[-FAIL]         
      } else {
        if (verbose) cat(" Tagging name of", NAMES[FAIL], "...\n")
        flush.console()
        for (i in FAIL) tempObj[[i]]$names <- paste("**", tempObj[[i]]$names, "**", sep = "") 
        flush.console()      
      }      
    }    
    cat("\n")
    if (class(object)[2] == "replist") {
      tempList <- c(tempList, tempObj)
      GROUP <- c(GROUP, length(tempObj))    
    }  
  }      
                   
  if (class(object)[2] == "replist") {
    class(tempList) <- c("modlist", "pcrfit")
    GROUP <- rep(1:length(GROUP), GROUP) 
    if (verbose) cat("Updating 'replist':\n")      
    tempObj <- replist(tempList, GROUP, verbose = verbose, ...)
  }
  
  class(tempObj) <- CLASS    
  return(tempObj)
}  
   