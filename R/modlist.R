modlist <- function(
x, 
cyc = 1, 
fluo = NULL, 
model = l4, 
check = "uni2",
checkPAR = parKOD(),
remove = c("none", "fit", "KOD"),
labels = NULL, 
norm = FALSE,
backsub = NULL,
smooth = c("none", "smooth", "lowess", "supsmu", "spline"), 
smoothPAR = list(span = 0.1), 
factor = 1,
opt = FALSE,
optPAR = list(sig.level = 0.05, crit = "ftest"),
verbose = TRUE,
...
)
{
  options(expressions = 50000)  
  remove <- match.arg(remove)
  smooth <- match.arg(smooth) 
        
  ## convert from single fit 
  if (class(x)[1] == "pcrfit") {
    model <- x$MODEL
    x <- x$DATA            
  }                         
  
  if (is.null(fluo)) fluo <- 2:ncol(x) 
        
  if (!is.null(backsub) & !is.numeric(backsub)) 
      stop("'backsub' must be either 'NULL'' or a numeric sequence!")   
  
  ## from 1.3-5: define label vector
  if (!is.null(labels)) {
    LABNAME <- deparse(substitute(labels))     
    LABELS <- labels
    if (length(LABELS) != length(fluo)) stop("Number of labels and runs do not match!")
  } else {
    LABELS <- 1:length(fluo)
    LABNAME <- "lab"
  }
    
  modLIST <- vector("list", length = length(fluo))   
   
  CYCLES <- x[, cyc]
  allFLUO <- x[, fluo, drop = FALSE]
  NAMES <- colnames(allFLUO)    
  
  for (i in 1:ncol(allFLUO)) {
    FLUO  <- allFLUO[, i]      
    NAME <- NAMES[i]
  
    ## normalization
    if (norm) FLUO <- rescale(FLUO, 0, 1)    
    
    ## background subtraction
    if (!is.null(backsub)) {
      BACK <- mean(FLUO[backsub], na.rm = TRUE)
      FLUO <- FLUO - BACK
    }  
      
    ## smoothing
    if (smooth != "none") {
      if (smooth == "smooth") FLUO <- as.numeric(smooth(FLUO))
      if (smooth == "lowess") FLUO <- lowess(FLUO, f = smoothPAR$span)$y
      if (smooth == "supsmu") FLUO <- supsmu(1:length(FLUO), FLUO, span = smoothPAR$span)$y
      if (smooth == "spline") FLUO <- spline(1:length(FLUO), FLUO, n = length(FLUO))$y         
    }
      
    ## changing magnitude
    if (factor != 1) FLUO <- FLUO * factor                
        
    DATA <- data.frame(Cycles = CYCLES, Fluo = FLUO)    
      
    if (verbose) cat("Making model for ", NAME, " (", model$name, ")\n", sep= "")  
    flush.console()
    
    fitOBJ <- try(pcrfit(DATA, 1, 2, model, verbose = FALSE, ...), silent = TRUE)
       
    ## tag names if fit failed
    if (inherits(fitOBJ, "try-error")) {  
      fitOBJ <- list()     
      if (verbose) cat(" => Fitting failed. Tagging name of ", NAME, "...\n", sep = "")  
      flush.console()
      NAME <- paste("*", NAME, "*", sep = "")                
      fitOBJ$DATA <- DATA
      fitOBJ$isFitted <- FALSE
      fitOBJ$isOutlier <- FALSE
      class(fitOBJ) <- "pcrfit"        
    } else {
      if (verbose) cat(" => Fitting passed...\n", sep = "")
      flush.console()
      fitOBJ$isFitted <- TRUE
      fitOBJ$isOutlier <- FALSE
    }
    
    ## optional model selection  
    if (opt) {
 	    fitOBJ2 <- try(mselect(fitOBJ, verbose = FALSE, sig.level = optPAR$sig.level, crit = optPAR$crit), silent = TRUE)             
      if (inherits(fitOBJ2, "try-error")) {
        if (verbose) cat(" => Model selection failed! Using original model...\n", sep = "")
        fitOBJ$isFitted <- TRUE
        fitOBJ$isOutlier <- FALSE
        flush.console()              
      } else {
        if (verbose) cat(" => Model selection passed...", sep = "")
        flush.console()
        fitOBJ <- fitOBJ2
        fitOBJ$isFitted <- TRUE
        fitOBJ$isOutlier <- FALSE
        if (verbose) cat(" => ", fitOBJ$MODEL$name, "\n", sep = "")
        flush.console()
      }
    }    
    
    if (verbose) cat("\n")    
    
    fitOBJ$call2$model <- fitOBJ$MODEL
    fitOBJ$call2$opt.method <- "all"
    fitOBJ$call2$nls.method <- "all"       
           
    modLIST[[i]] <- fitOBJ
    modLIST[[i]]$names <- NAME    
  }  
         
  ## from 1.3-5: sigmoidal outlier detection by KOD
  if (!is.null(check)) {
    class(modLIST) <- c("modlist", "pcrfit")
    OUTL <- KOD(modLIST, method = check, par = checkPAR, plot = FALSE)   
    modLIST <- OUTL
  }
  
  ## from 1.3-5: remove failed fits, update label vector, assign new vector to global environment
  if (remove != "none") {
    logVEC <- vector("numeric", length = length(modLIST))
    ## set failed fits to 1
    if (remove %in% c("fit", "KOD")) {
      SEL <- sapply(modLIST, function(x) x$isFitted)
      logVEC[SEL == FALSE] <- 1      
    }
    ## set KOD's to 1
    if (remove == "KOD") {      
      SEL <- sapply(modLIST, function(x) x$isOutlier)
      logVEC[SEL == TRUE] <- 1      
    }
    ## remove and update LABELS vector
    SEL <- which(logVEC == 1)
    
    if (length(SEL) > 0) {
      if (verbose) cat(" => Removing from fit:", NAMES[SEL], "... \n", sep = " ")
      flush.console()
      modLIST <- modLIST[-SEL]      
      if (verbose) cat(" => Updating", LABNAME, "and Writing", paste(LABNAME, "_mod", sep = ""), "to global environment...\n\n", sep = " ")
      flush.console()    
      LABELS <- LABELS[-SEL]    
      assign(paste(LABNAME, "_mod", sep = ""), LABELS, envir = .GlobalEnv)        
    } 
  }
 
  class(modLIST) <- c("modlist", "pcrfit")
  invisible(modLIST)
}

