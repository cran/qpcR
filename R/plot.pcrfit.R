plot.pcrfit <- function(
x, 
which = c("all", "single", "3D" ),
fitted = TRUE, 
add = FALSE,
col = NULL, 
confband = c("none", "confidence", "prediction"),
errbar = c("none", "sd", "se", "conf"),
par3D = list(),
par2D = list(),
parCI = list(),
parSD = list(), 
...) 
{
  require(rgl, quietly = TRUE)
  object <- x
  confband <- match.arg(confband)
  errbar <- match.arg(errbar)
  which <- match.arg(which)    
  if (class(object) != "modlist") modList <- list(object) else modList <- object      
  cycList <- lapply(modList, function(x) x$DATA[, 1])  
  fluoList <- lapply(modList, function(x) x$DATA[, 2])
  minVal <- min(sapply(fluoList, function(x) min(x, na.rm = TRUE)), na.rm = TRUE)   
  maxVal <- max(sapply(fluoList, function(x) max(x, na.rm = TRUE)), na.rm = TRUE)   
  CYC <- unique(matrix(do.call(cbind.na, cycList), ncol = 1))
  CYC <- CYC[!is.na(CYC)]
  LEN <- length(modList)
  NAMES <- sapply(modList, function(x) x$names)     
  
  if (is.null(col)) {
    colvec <- rep(1, LEN)    
    if (class(object)[2] == "replist") colvec <- rainbow(attr(object, "nlevels"))     
  } else colvec <- rep(col, length.out = LEN)   
    
  if (which == "3D") {
    do.call(plot3d, modifyList(list(x = CYC, y = 1:LEN, z = maxVal, type = "n", axes = FALSE, box = FALSE, xlab = "", 
           ylab = "", zlab = "", zlim = c(0, 1.1 * maxVal)), par3D))
    do.call(axis3d, modifyList(list('x', at = pretty(CYC), cex = 0.5), par3D))
    do.call(mtext3d, modifyList(list("Cycles", 'x', line = 2), par3D))     
    do.call(axis3d, modifyList(list('y', at = 1:LEN, label = NAMES, cex = 0.5), par3D))
    do.call(mtext3d, modifyList(list("Run", 'y', line = 2), par3D))
    do.call(axis3d, modifyList(list('z', cex = 0.5), par3D))
    do.call(mtext3d, modifyList(list("Fluo", 'z', line = 2), par3D))
  }   
  
  if (which == "all" && !add)   
    do.call(plot, modifyList(list(CYC, rep(maxVal, length(CYC)), ylim = c(minVal, maxVal), 
         xlab = "Cycles", ylab = "Raw fluorescence", type = "n"), par2D)) 
  
  if (which == "single") {
    DIM <- ceiling(sqrt(LEN))
    par(mfrow = c(DIM, DIM))
    par(mar = c(0.2, 0.2, 1, 0.2))
  } 
  
  for (i in 1:LEN) {
    DATA <- modList[[i]]$DATA
    CC <- complete.cases(DATA[, 2])
    DATA <- DATA[CC, ] 
    FITTED <- fitted(modList[[i]])[CC]     
         
    if (which == "3D") {
      do.call(points3d, modifyList(list(x = DATA[, 1], y = i, z = DATA[, 2], color = colvec[i]), par3D))
      if (!is.null(FITTED) && fitted) do.call(lines3d, modifyList(list(x = DATA[CYC, 1], y = i, z = FITTED[CYC], color = colvec[i]), par3D))      
    }
    
    if (which == "all") {
      do.call(points, modifyList(list(DATA[, 1], DATA[, 2], col = colvec[i]), par2D))
      if (!is.null(FITTED) && fitted) do.call(lines, modifyList(list(DATA[CYC, 1], FITTED[CYC], col = colvec[i]), par2D)) 
    } 
    
    if (which == "single") {
      NAME <- NAMES[i]
      if (grepl("\\*\\*[[:alnum:]]*", NAME)) colMain <- "blue" 
        else if (grepl("\\*[[:alnum:]]*", NAME)) colMain <- "red"
          else colMain <- "black"
      TRY <- try(do.call(plot, modifyList(list(DATA[, 1], DATA[, 2], main = NAME, cex.main = 0.7, col.main = colMain, type = "p", 
                         xlab = FALSE, ylab = FALSE, xaxt = "n", yaxt = "n", col = colvec[i]), par2D)), silent = TRUE)
      if (inherits(TRY, "try-error")) next      
      if (!is.null(FITTED) && fitted) do.call(lines, modifyList(list(DATA[CYC, 1], FITTED[CYC], col = colvec[i]), par2D))      
    } 
    
    if (confband != "none") {      
      CONFINT <- predict(modList[[i]], interval = confband, ...)
      do.call(lines, modifyList(list(CYC, CONFINT$Lower[CYC], col = 2), parCI))
      do.call(lines, modifyList(list(CYC, CONFINT$Upper[CYC], col = 2), parCI))
    }
    
    if (errbar != "none") {      
      if (class(object)[2] != "replist") stop("Error bars only possible on a 'replist'!")
      DATA <- as.data.frame(DATA)
      colnames(DATA) <- c("CYC", "FLUO")
      DATA2 <- t(unstack(DATA, FLUO ~ CYC))
      STAT <- switch(errbar, sd = apply(DATA2, 1, function(x) sd(x, na.rm = TRUE)),
                     se = apply(DATA2, 1, function(x) sd(x, na.rm = TRUE)/sqrt(length(x[complete.cases(x)]))))
      if (errbar == "conf") {
        CONFINT <- predict(modList[[i]], interval = "conf", ...)          
        do.call(arrows, modifyList(list(DATA[, 1], CONFINT$Lower, DATA[, 1], CONFINT$Upper, angle = 90, code = 3, col = colvec[i], length = 0.05), parSD))
      } else do.call(arrows, modifyList(list(DATA[, 1], FITTED + STAT, DATA[, 1], FITTED - STAT, angle = 90, code = 3, col = colvec[i], length = 0.05), parSD))
    } 
  }     
}  