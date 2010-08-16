rep2mod <- function(rl)
{
  if (class(rl)[2] != "replist") stop("object must be a 'replist'!")
  modList <- NULL
  for (i in 1:length(rl)) {
    temp <- rl[[i]]$modlist
    modList <- c(modList, temp)
  }
  nlev <- attr(rl, "nlevels")
  nitem <- attr(rl, "nitems")
  CLASS <- rep(1:nlev, nitem)
  class(modList) <- c("modlist", "pcrfit")
  attr(modList, "group") <- CLASS
  return(modList)
}


  

