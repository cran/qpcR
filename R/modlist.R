modlist <- function(
x, 
cyc = 1, 
fluo = NULL, 
model = l4, 
remove = FALSE,
opt = FALSE, 
norm = FALSE,
backsub = NULL,
smooth = c("none", "tukey", "lowess", "supsmu", "spline"), 
span = 0.1, 
factor = 1, 
opt.method =  "all",
nls.method = "all",
sig.level = 0.05, 
crit = "ftest",
verbose = TRUE,
...
)
{
  options(expressions = 50000)  
  smooth <- match.arg(smooth)    
    
  if (class(x)[1] == "pcrfit") {
    model <- x$MODEL
    x <- x$DATA            
  }                         
  
  if (is.null(fluo)) fluo <- 2:ncol(x) 
  modLIST <- list()       
    
  if (!is.null(backsub) && !is.numeric(backsub)) 
      stop("'backsub' must be either NULL or a numeric sequence!")   
  
  COUNTER <- 1
  
  for (i in fluo) {
    CYCLES <- x[, cyc]      
    FLUO  <- x[, i]      
    NAME <- colnames(x)[i]
  
    if (norm) {
      FLUO <- FLUO - min(FLUO, na.rm = TRUE)
      FLUO <- FLUO/max(FLUO, na.rm = TRUE)
    }
        
    if (!is.null(backsub)) {
      BACK <- mean(FLUO[backsub], na.rm = TRUE)
      FLUO <- FLUO - BACK
    }    
    
   if (smooth != "none") {
      if (smooth == "tukey") FLUO <- as.numeric(smooth(FLUO))
      if (smooth == "lowess") FLUO <- lowess(FLUO, f = span)$y
      if (smooth == "supsmu") FLUO <- supsmu(1:length(FLUO), FLUO, span = span)$y
      if (smooth == "spline") FLUO <- spline(1:length(FLUO), FLUO, n = length(FLUO))$y         
    }
    
    if (factor != 1) FLUO <- FLUO * factor                
        
    DATA <- data.frame(Cycles = CYCLES, Fluo = FLUO)     
    flush.console()
    if (verbose) cat("Making model for ", NAME, " (", model$name, ")\n", sep= "")                
    fitObj <- try(pcrfit(DATA, 1, 2, model, opt.method = opt.method, nls.method = nls.method, verbose = verbose, ...), silent = TRUE)
    
    if (inherits(fitObj, "try-error")) {  
      fitObj <- list()     
      if (verbose) cat(" => gave a fitting error!\n => Tagging name of ", NAME, "...\n", sep = "")  
      NAME <- paste("*", NAME, "*", sep = "")                
      fitObj$DATA <- DATA
      fitObj$isFitted <- FALSE
      class(fitObj) <- "pcrfit"
      if(remove) {
        if (verbose) cat(" => Removing ", NAME, "... \n\n", sep = "")
        next           
      }                
    } else fitObj$isFitted <- TRUE
      
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
           
    modLIST[[COUNTER]] <- fitObj
    modLIST[[COUNTER]]$names <- NAME    
    
    COUNTER <- COUNTER + 1      
  }    
     
  class(modLIST) <- c("modlist", "pcrfit")
  invisible(modLIST)
}

