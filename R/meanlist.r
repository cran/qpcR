meanlist <- function(object, group = NULL, type = c("mean", "median"))
{
  if (class(object) != "modlist") stop("Please supply an object of class 'modlist'!")
  if (is.null(group)) stop("Please define replicate groups!")
  if (length(group) != length(object)) stop("length of 'group' and 'object' must match!")   
  group <- as.factor(group)
  type <- match.arg(type)
  
  splitVec <- split(1:length(object), group)
  repMod <- list()
  meanDATA <- list()
  
  for (i in 1:length(splitVec)) {
    tempCYC <- sapply(splitVec[[i]], function(x) object[[x]]$DATA[, 1])
    tempCYC <- matrix(tempCYC, ncol = 1)
    tempFLUO <- sapply(splitVec[[i]], function(x) object[[x]]$DATA[, 2])
    tempFLUO <- matrix(tempFLUO, ncol = 1)
    tempDATA <- cbind(tempCYC, tempFLUO)
    colnames(tempDATA) <- c("Cycles", "Fluo")
    repMod[[i]] <- tempDATA
    meanFLUO <- tapply(tempDATA[, 2], as.factor(tempDATA[, 1]), function(x) switch(type, mean = mean(x, na.rm = TRUE), median = median(x, na.rm = TRUE)))
    uniqueCYC <- unique(tempDATA[, 1])
    tempDATA <- cbind(uniqueCYC, meanFLUO)
    colnames(tempDATA) <- c("Cycles", "Fluo")
    meanDATA[[i]] <- tempDATA
  }
  
  finMod <- list()
  for (i in 1:length(meanDATA)) {     
    finMod[[i]] <- pcrfit(meanDATA[[i]], 1, 2, model = object[[1]]$MODEL)
    finMod[[i]]$isReps <- FALSE
  }
  class(finMod) <- c("modlist", "meanlist", "pcrfit")
  attr(finMod, "nlevels") <- nlevels(group)
  attr(finMod, "nitems") <- length(group)/nlevels(group)
  return(finMod)
}