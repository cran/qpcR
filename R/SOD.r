SOD <- function(
object,    
remove = FALSE,
verbose = TRUE, 
...
)
{
  CLASS <- class(object)
  if (!(class(object)[1] %in% c("pcrfit", "modlist", "replist"))) stop("object must be of class 'pcrfit', 'modlist' or 'replist'!")
  if (class(object)[1] == "pcrfit") tempObj <- list(object) else tempObj <- object
  if (class(object)[2] == "replist") {
    tempObj <- rep2mod(object)
    GROUP <- attr(tempObj, "group")      
  }      
         
  effList <- lapply(tempObj, function(x) try(efficiency(x, plot = FALSE, ...), silent = TRUE)) 
  if (verbose) cat("\nCalculating second derivative maximum...\n")
  flush.console()
  cpD2 <- sapply(effList, function(x) if(inherits(x, "try-error")) 0 else x$cpD2)  
  if (verbose) cat("Calculating first derivative maximum...\n")   
  flush.console()
  cpD1 <- sapply(effList, function(x) if(inherits(x, "try-error")) 0 else x$cpD1) 
  if (verbose) cat("Calculating R-square...\n")  
  flush.console()
  RSQ <- sapply(tempObj, function(x) if(inherits(x, "try-error")) 0 else Rsq.ad(x)) 
  NAMES <- sapply(tempObj, function(x) x$names)  
  if (verbose) cat("Checking for sigmoidal consistency...\n")
  FAIL <- which(cpD2 > cpD1 |  cpD1 - cpD2 > 10 | RSQ < 0.9)
  if (length(FAIL) == 0) PASS <- 1:length(tempObj) else PASS <- (1:length(tempObj))[-FAIL]
    
  for (i in FAIL) tempObj[[i]]$outlier <- TRUE
  for (j in PASS) tempObj[[j]]$outlier <- FALSE  
       
  if (length(FAIL) > 0) {
    if (verbose) cat(" Found non-sigmoidal structure for", NAMES[FAIL], "...\n", sep = " ")  
    flush.console() 
    if (remove) {
     if (verbose) cat(" Removing", NAMES[FAIL], "...\n", sep = " ")
     flush.console()
     tempObj <- tempObj[-FAIL]
     if (class(object)[2] == "replist") GROUP <- GROUP[-FAIL]
    } else {
      if (verbose) cat(" Tagging name of", NAMES[FAIL], "...\n", sep = " ")
      flush.console()
      for (i in FAIL) tempObj[[i]]$names <- paste("**", tempObj[[i]]$names, "**", sep = "") 
      flush.console()      
    }      
  } 
  
  cat("\n")
    
  if (class(object)[1] == "pcrfit") tempObj <- tempObj[[1]]
  if (class(object)[2] == "replist") {
    class(tempObj) <- c("modlist", "pcrfit")  
    cat("Updating 'replist':\n")    
    tempObj <- replist(tempObj, GROUP, verbose = verbose, ...)
  }
  
  class(tempObj) <- CLASS      
  return(tempObj)
} 