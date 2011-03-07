replist <- function(
object, 
group = NULL, 
remove = TRUE,
opt = FALSE, 
verbose = TRUE, 
...)
{
  if (class(object) != "modlist") stop("Please supply an object of class 'modlist'!")
  if (is.null(group)) stop("Please define replicate groups!")
  if (length(group) != length(object)) stop("length of 'group' and 'object' must match!")    
  
  if (remove) {
    CLASS <- sapply(object, function(x) class(x)[2]) 
    NAMES <- sapply(object, function(x) x$names)  
    FAILS <- which(is.na(CLASS))     
    if (length(FAILS) > 0) {
      cat("Removing tagged runs and updating 'group':", NAMES[FAILS], "\n\n")
      object <- object[-FAILS]
      group <- group[-FAILS]
    }                     
  } 
    
  group <- as.factor(group)  
  
  splitVEC <- split(1:length(object), group)      
  nameVEC <- sapply(object, function(x) x$names)
  nameVEC <- split(nameVEC, group)  
  
  repMOD <- list()
  
  for (i in 1:length(splitVEC)) {
    CYCS <- NULL
    FLUO <- NULL  
     
    for (j in splitVEC[[i]]) {
      CYCS <- c(CYCS, object[[j]]$DATA[, 1])     
      FLUO <- c(FLUO, object[[j]]$DATA[, 2])     
      DATA <- cbind(Cycles = CYCS, Fluo = FLUO)    
    }  
    
    if (verbose) cat("Making model for replicates:", nameVEC[[i]], "\n") 
    flush.console()
    fitObj <- try(pcrfit(DATA, 1, 2, model = object[[splitVEC[[i]][1]]]$MODEL), silent = TRUE)
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
    
    if (verbose) cat(" => ", fitObj$MODEL$name, "\n\n", sep = "")
    flush.console()
    
    repMOD[[i]] <- fitObj
    repMOD[[i]]$isReps <- TRUE
    repMOD[[i]]$names <- paste("group_", i, sep = "") 
    repMOD[[i]]$DATA <- DATA     
    repMOD[[i]]$modlist <- object[unlist(splitVEC[i])]     
  }      
  
  class(repMOD) <- c("modlist", "replist", "pcrfit")
  attr(repMOD, "nlevels") <- nlevels(group)
  attr(repMOD, "nitems") <- as.numeric(table(group))
  attr(repMOD, "group") <- group
  return(repMOD)
}