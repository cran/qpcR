replist <- function(object, group = NULL, opt = FALSE, ...)
{
  if (class(object) != "modlist") stop("Please supply an object of class 'modlist'!")
  if (is.null(group)) stop("Please define replicate groups!")
  if (length(group) != length(object)) stop("length of 'group' and 'object' must match!")   
  group <- as.factor(group)
  
  splitVec <- split(1:length(object), group)      
  nameVec <- sapply(object, function(x) x$names)
  nameVec <- split(nameVec, group)  
  
  repMod <- list()
  meanDATA <- list()
  tempModList <- list()
  tempDATA <- NULL
  
  for (i in 1:length(splitVec)) {   
   for (j in splitVec[[i]]) {      
    tempDATA <- rbind(tempDATA, object[[j]]$DATA)      
   }  
   tempModList[[i]] <- object[splitVec[[i]]]               
   repMod[[i]] <- tempDATA 
   meanFLUO <- tapply(tempDATA[, 2], as.factor(tempDATA[, 1]), function(x) mean(x, na.rm = TRUE))    
   uniqueCYC <- unique(tempDATA[, 1])
   tempDATA <- cbind(uniqueCYC, meanFLUO)
   colnames(tempDATA) <- c("Cycles", "Fluo")
   meanDATA[[i]] <- tempDATA  
   tempDATA <- NULL            
  }    
  
  finMod <- list()
  
  for (i in 1:length(meanDATA)) {           
    flush.console()
    cat("Making model for replicates:", nameVec[[i]]) 
    meanMod <- try(pcrfit(meanDATA[[i]], 1, 2, model = object[[splitVec[[i]][1]]]$MODEL), silent = TRUE)
    if (inherits(meanMod, "try-error")) stop("There was an error for the starting values!")
    fitObj <- try(pcrfit(repMod[[i]], 1, 2, model = object[[splitVec[[i]][1]]]$MODEL, do.optim = FALSE, start = coef(meanMod)), silent = TRUE)
    if (inherits(fitObj, "try-error")) cat(" => gave a fitting error!\n", sep = "")
    
    if (opt) {
 	    fitObj2 <- try(mselect(fitObj, verbose = FALSE, ...))             
      if (inherits(fitObj2, "try-error")) {
        fitObj <- fitObj
        cat(" => gave a model selection error!", sep = "")
      } else {
        fitObj <- fitObj2           
      }
    }     
    cat(" => ", fitObj$MODEL$name, "\n", sep = "")
    
    finMod[[i]] <- fitObj
    finMod[[i]]$isReps <- TRUE
    finMod[[i]]$names <- paste("group_", i, sep = "") 
    finMod[[i]]$DATA <- repMod[[i]]
    finMod[[i]]$modlist <- tempModList[[i]]     
  }      
  class(finMod) <- c("modlist", "replist", "pcrfit")
  attr(finMod, "nlevels") <- nlevels(group)
  attr(finMod, "nitems") <- length(group)/nlevels(group)     
  return(finMod)
}