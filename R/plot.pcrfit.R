plot.pcrfit <- function(
x, 
fitted = TRUE, 
confband = c("none", "confidence", "prediction"),
errbar = c("none", "sd", "se", "conf"),
add = FALSE, 
col = NULL,
level = 0.95,
xlim = NULL,  
which = c("all", "single", "3D" ),
...) 
{
  require(rgl, quietly = TRUE)
  object <- x
  confband <- match.arg(confband)
  errbar <- match.arg(errbar)
  which <- match.arg(which)
  if (!is.null(xlim) && length(xlim) != 2) stop("'xlim' must be a 2-element vector!") 
  
  if (class(object) != "modlist") modList <- list(object) else modList <- object

  minVal <- min(unlist(sapply(modList, function(x) x$DATA[, 2])), na.rm = TRUE)
  maxVal <- max(unlist(sapply(modList, function(x) x$DATA[, 2])), na.rm = TRUE)
  CYC <- unique(as.vector(unlist(sapply(modList, function(x) x$DATA[, 1]))))

  LEN <- length(modList)
    
  if (is.null(col)) {
    colvec <- rep(1, LEN)     
    if (class(object)[2] == "replist") colvec <- gl(attr(object, "nlevels"), 1)
  } else colvec <- rep(col, length.out = LEN)  
  
  if (!add && which == "3D") {     
    plot3d(x = CYC, y = 1:LEN, z = maxVal, type = "n", axes = FALSE, box = FALSE,
           xlab = "", ylab = "", zlab = "", zlim = c(0, 1.1 * maxVal))       
    
    for (i in 1:LEN) {
     DATA <- modList[[i]]$DATA        
     CC <- complete.cases(DATA[, 2])       
     DATA <- DATA[CC, ] 
     FITTED <- fitted(modList[[i]])     
     if (nrow(DATA) > 0) points3d(x = DATA[, 1], y = i, z = DATA[, 2], col = colvec[i]) else next
     if (is.null(FITTED)) next
     if (class(modList)[2] != "replist") lines3d(x = DATA[, 1], y = i, z = FITTED[CC], col = colvec[i])
     else lines3d(x = DATA[CYC, 1], y = i, z = FITTED[CYC], col = colvec[i])
    }     
    
    axis3d('x', at = pretty(CYC), cex = 0.5)
    mtext3d("Cycles", 'x', line = 2) 
    NAMES <- sapply(modList, function(x) x$names)    
    axis3d('y', at = 1:LEN, label = NAMES, cex = 0.5)
    mtext3d("Run", 'y', line = 2)    
    axis3d('z', cex = 0.5)
    mtext3d("Fluo", 'z', line = 2)    
  return("Finished...")  
  }      
   
  if (!add && which == "all") plot(CYC, rep(maxVal, length(CYC)), ylim = c(minVal, maxVal), xlim = c(xlim[1], xlim[2]), type = "n")
  else if (!add && which == "single") {
    DIM <- ceiling(sqrt(LEN))
    par(mfrow = c(DIM, DIM))
    par(mar = c(0.2, 0.2, 1, 0.2))
  }     
    
  for (i in 1:LEN) { 
    if (class(modList[[i]]) != "pcrfit") stop("object must be of class 'pcrfit'")
    DATA <- modList[[i]]$DATA 
    CC <- complete.cases(DATA[, 2])   
    DATA <- DATA[CC, ]    
    FITTED <- fitted(modList[[i]])[CC] 
    
    if (!add && which == "single") 
      TRY <- try(plot(DATA, xlim = c(xlim[1], xlim[2]), main = modList[[i]]$names, cex.main = 0.7, type = "n", xlab = FALSE, ylab = FALSE, xaxt = "n", yaxt = "n"), silent = TRUE)
    else 
      TRY <- NA  
    if (inherits(TRY, "try-error")) next 
   
    points(DATA, col = colvec[i], ...)
       
    if (is.null(modList[[i]]$isReps)) {             
      if (!is.null(FITTED)) lines(DATA[, 1], FITTED[CC], col = colvec[i])      
    } else {          
      if (!is.null(FITTED)) lines(DATA[CYC, 1], FITTED[CYC], col = colvec[i])       
    }      
      
    if (confband != "none" || errbar == "conf") {
      if (errbar == "conf") confband <- "confidence"                 
      CONFINT <- predict(modList[[i]], interval = confband, level = level, ...)[unique(DATA[, 1]), ]           
      lines(CYC, CONFINT$Lower, col = 2, ...)
      lines(CYC, CONFINT$Upper, col = 2, ...)  
    }
    
    if (errbar != "none") {     
     DATA <- as.data.frame(DATA)       
     colnames(DATA) <- c("CYC", "FLUO")        
     DATA2 <- t(unstack(DATA, FLUO ~ CYC))      
     
     STAT <- switch(errbar, sd = apply(DATA2, 1, function(x) sd(x, na.rm = TRUE)), 
                            se = apply(DATA2, 1, function(x) sd(x, na.rm = TRUE)/sqrt(length(x[complete.cases(x)]))))
     if (errbar == "conf") arrows(DATA[, 1], CONFINT$Lower, DATA[, 1], CONFINT$Upper, angle = 90, code = 3, col = colvec[i], length = 0.05)
     else arrows(DATA[, 1], FITTED + STAT, DATA[, 1], FITTED - STAT, angle = 90, code = 3, col = colvec[i], length = 0.05)            
    }      
  }                
}  