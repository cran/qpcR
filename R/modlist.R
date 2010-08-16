modlist <- function(
x, 
cyc = 1, 
fluo = NULL, 
model = l4, 
remove = FALSE,
opt = FALSE, 
norm = FALSE,
backsub = NULL, 
opt.method =  "all",
nls.method = "all",
sig.level = 0.05, 
crit = "ftest",
verbose = TRUE,
...
)
{
  modList <- NULL
  counter <- 1   
  
  if (class(x)[1] == "pcrfit") {
    model <- x$MODEL
    x <- x$DATA            
  }                         
  
  if (is.null(fluo)) fluo <- 2:ncol(x)        
    
  if (!is.null(backsub) && !is.numeric(backsub)) 
      stop("'backsub' must be either NULL or a numeric sequence!")   

  for (i in fluo) {
    Cycles <- x[, cyc]      
    Fluo  <- x[, i]
    compl <- complete.cases(Fluo)
    Cycles <- Cycles[compl]
    Fluo <- Fluo[compl]
    NAME <- colnames(x)[i]
  
    if (norm) {
      Fluo <- Fluo - min(Fluo, na.rm = TRUE)
      Fluo <- Fluo/max(Fluo, na.rm = TRUE)
    }
        
    if (!is.null(backsub)) {
      back <- mean(Fluo[backsub], na.rm = TRUE)
      Fluo <- Fluo - back
    }                    
        
    DATA <- data.frame(Cycles = Cycles, Fluo = Fluo)     
    flush.console()
    if (verbose) cat("Making model for ", NAME, " (", model$name, ")\n", sep= "")                
    fitObj <- try(pcrfit(DATA, 1, 2, model, opt.method = opt.method, nls.method = nls.method, verbose = verbose, ...), silent = TRUE)
    
    if (inherits(fitObj, "try-error")) {  
      fitObj <- list()     
      if (verbose) cat(" => gave a fitting error!", sep = "")  
      if(remove) {
        if (verbose) cat(" => Removing ", NAME, "...\n\n", sep = "")
        next
      } 
      if (verbose) cat(" => Tagging name of ", NAME, "...\n", sep = "")
      NAME <- paste("*", NAME, "*", sep = "")                
      fitObj$DATA <- DATA
      class(fitObj) <- "pcrfit"      
    }
      
    if (opt) {
 	    fitObj2 <- try(mselect(fitObj, verbose = FALSE, sig.level = sig.level, crit = crit, ...), silent = TRUE)             
      if (inherits(fitObj2, "try-error")) {
        fitObj <- fitObj
        if (verbose) cat(" => ", colnames(x[i]), " gave a model selection error!\n", sep = "")
      } else {
        fitObj <- fitObj2
        cat(" => ", fitObj$MODEL$name, sep = "")
      }
    }   
    
    if (verbose) cat("\n")
    
    fitObj$call2$model <- fitObj$MODEL
    fitObj$call2$opt.method <- opt.method
    fitObj$call2$nls.method <- nls.method       
           
    modList[[counter]] <- fitObj
    modList[[counter]]$names <- NAME   
    counter <- counter + 1       
  }
  
  class(modList) <- c("modlist", "pcrfit")
  invisible(modList)
}

