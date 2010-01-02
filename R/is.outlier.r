is.outlier <- function(object)
{
  if (class(object)[1] != "modlist") stop("Please supply a 'modlist' or 'replist'!")
  if (class(object)[2] != "replist") {
    OUTL <- sapply(object, function(x) x$outlier)
    NAMES <- sapply(object, function(x) x$names)
    names(OUTL) <- NAMES
  } else {
    OUTL <- NULL
    NAMES <- NULL
    for (i in 1:length(object)) {
      tempMod <- object[[i]]$modlist
      OUTL <- c(OUTL, sapply(tempMod, function(x) x$outlier))
      NAMES <- c(NAMES, sapply(tempMod, function(x) x$names))
    }
    names(OUTL) <- NAMES
  }
  return(OUTL)
}