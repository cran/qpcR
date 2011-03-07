is.outlier <- function(object)
{
  if (!(class(object) %in% c("modlist", "replist"))) stop("object must be of class 'modlist' or 'replist'!") 
  if (class(object)[2] == "replist") object <- rep2mod(object)  
      
  OUTL <- sapply(object, function(x) x$isOutlier)
  NAMES <- sapply(object, function(x) x$names)
  names(OUTL) <- NAMES
  return(OUTL)
}