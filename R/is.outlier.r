is.outlier <- function(object)
{
  if (!(class(object) %in% c("modlist", "replist", "pcrfit"))) stop("object must be of class 'pcrfit', 'modlist' or 'replist'!") 
  if (class(object)[1] == "pcrfit") object <- list(object) 
  else if (class(object)[2] == "replist") object <- rep2mod(object)
  else object <- object
      
  OUTL <- sapply(object, function(x) x$outlier)
  NAMES <- sapply(object, function(x) x$names)
  names(OUTL) <- NAMES
  return(OUTL)
}