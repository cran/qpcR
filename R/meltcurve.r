meltcurve <- function(data, temps = NULL, fluos = NULL, window = NULL, span = 11, ...)
{
  if (is.null(temps)) temps <- seq(from = 1, to = ncol(data), by = 2) 
  if (is.null(fluos)) fluos <- seq(from = 2, to = ncol(data), by = 2)
  if (length(temps) != length(fluos)) stop("Numbers of temperature columns and fluorescence columns do not match!")
  if (is.null(window)) window <- range(matrix(data[, temps], ncol = 1))   
  NAMES <- colnames(data[, fluos])    
  
  ### create dataframes with temp and fluo values
  TEMPS <- data[, temps, drop = FALSE]
  FLUOS <- data[, fluos, drop = FALSE]      
    
  meltDATA <- NULL
  tempDATA  <- NULL
  derivDATA <- NULL     
  tmDATA <- list()      
       
  for (i in 1:ncol(TEMPS)) {
    TEMP <- TEMPS[, i]
    FLUO <- FLUOS[, i]      
    
    ### cut off unimportant temperature regions
    SEL <- which(TEMP <= window[1] | TEMP > window[2])
    TEMP <- TEMP[-SEL]
    FLUO <- FLUO[-SEL]       
   
    ### cubic spline fitting and Friedman's Supersmoother
    ### on the first derivative curve
    SPLFN <- splinefun(TEMP, FLUO)
    meltDATA <- cbind.na(meltDATA, SPLFN(TEMP))
    tempDATA <- cbind.na(tempDATA, TEMP)    
    derivVEC <- SPLFN(TEMP, deriv = 1)
    SMOOTH <-  supsmu(TEMP, derivVEC)        
    derivDATA <- cbind.na(derivDATA, -SMOOTH$y)           
    
    ### find peaks in first derivative data
    PEAKS <- qpcR:::peaks(-SMOOTH$y, span = span, ...)$x
    TM <- TEMP[PEAKS] 
    tmDATA[[i]] <- TM                                      
  }     
  tempDATA <- tempDATA[, -1, drop = FALSE]
  meltDATA <- meltDATA[, -1, drop = FALSE]
  derivDATA <- derivDATA[, -1, drop = FALSE]     
 
  ### plot raw melt data and first derivatives
  ### including identified melting points   
  COL <- rainbow(ncol(meltDATA))
  DIM <- ceiling(sqrt(ncol(TEMPS)))   
  par(mfrow = c(DIM, DIM))   
   
  for (i in 1:ncol(tempDATA)) {   
    qpcR::: xyy.plot(tempDATA[, i], meltDATA[, i], derivDATA[, i], main = NAMES[i], 
                     y1.par = list(xlab = "", ylab = "", type = "l", lwd = 2),
                     y2.par = list(xlab = "", ylab = "", type = "l", lwd = 2, text = ""),
                     first = par(mar = c(3, 2, 2, 3)))
    abline(v = tmDATA[[i]], lwd = 1, lty = 2)          
  }  
  
  return(tmDATA)           
}     