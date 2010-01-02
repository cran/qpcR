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
  if (class(object)[1] != "modlist") stop("Please supply either a 'modlist' or 'replist'!")
  tempList <- list()
  tempGroup <- NULL   
  
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
    if (verbose) cat("Calculating efficiencies...\n")
    eff <- sapply(tempObj, function(x) switch(efftype, sliwin = sliwin(x, plot = FALSE, ...)$eff,
                                                    sigfit = efficiency(x, plot = FALSE, ...)$eff,
                                                    expfit = expfit(x, plot = FALSE, ...)$eff))    

    if (verbose) cat("Calculating outlier(s)...\n")
    outl <- switch(method, bar = baretal(eff), pam = pamout(eff))    
    if (length(outl) == 0) n.outl <- 1:length(tempObj) else n.outl <- (1:length(tempObj))[-outl]           
        
    ### in case of 'replist'          
    if (class(object)[2] == "replist") {           
      for (j in n.outl) { 
        flush.console()
        if (verbose) cat("Tagging", object[[i]]$names, "/", object[[i]]$modlist[[j]]$names, "as non-outlier...\n")               
        object[[i]]$modlist[[j]]$outlier <- FALSE         
      }           
      for (k in outl) {
        flush.console()
        if (verbose) cat("Tagging", object[[i]]$names, "/", object[[k]]$modlist[[j]]$names, "as outlier...\n")             
        object[[i]]$modlist[[k]]$outlier <- TRUE
      }
      if (remove) {          
        tempMod <- object[[i]]$modlist
        tempMod <- tempMod[n.outl]
        tempList <- c(tempList, tempMod)
        tempGroup <- c(tempGroup, rep(i, length(tempMod)))        
      }      
      object[[i]]$outlier <- outl                     
    }                                       
    ########################
    
    ### in case of 'modlist'
    if (class(object)[2] != "replist") {       
      for (j in n.outl) {
        flush.console()
        if (verbose) cat("Tagging", object[[j]]$names, "as non-outlier...\n")        
        object[[j]]$outlier <- FALSE
      }
      for (k in outl) {
        flush.console()
        if (verbose) cat("Tagging", object[[k]]$names, "as outlier...\n")
        object[[k]]$outlier <- TRUE               
      }     
      if (remove) {             
          object <- object[-outl] 
          cat("removed!\n")   
          class(object) <- c("modlist", "pcrfit")
      } else cat("\n")                         
    }
    #######################
    
    cat("\n")
  }
  
  if (class(object)[2] == "replist" && remove) {
    class(tempList) <- c("modlist", "pcrfit")
    if (verbose) cat("Updating 'replist'...\n")
    object <- replist(tempList, tempGroup)
  }    
  return(object)
}     