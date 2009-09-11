pcropt2 <- function(object, plot = TRUE, which.par = "all", quan = 0.1, delete = c("low", "high"), ...)
{
  delete <- match.arg(delete)
  DATA <- eval(object$data)
  modelLIST <- list()
  parNAMES <- names(coef(object))
  if (which.par != "all" && !all(which.par %in% names(coef(object)))) stop("Parameter(s) '", which.par, "' not a coefficient of 'object'!")
  if (which.par == "all") parSEL <- parNAMES else parSEL <- which.par

  for (i in 1:nrow(DATA)) {
    newDATA <- DATA[-i, ]
    newMODEL <- update(object, data = newDATA, fluo = 2)
    modelLIST[[i]] <- newMODEL
    if (plot) plot(newMODEL)
  }
  
  coefMAT <- t(sapply(modelLIST, function(z) coef(z)))      
  dfbMAT <- t(sapply(modelLIST, function(z) abs(coef(z) - coef(object))/(summary(object)$parameters[, 2])))
  colSEL <- match(parSEL, colnames(dfbMAT))         
      
  if (plot) {
    par(mfrow = c(2, 1))
    par(mar = c(1, 3, 1, 1))
    plot(object, ...)
  }
  
  allSELECT <- NULL
       
  for (i in colSEL) {
    statVEC <- dfbMAT[, i]  
    if (delete == "low") SELECT <- which(statVEC <= quantile(statVEC, quan))
     else SELECT <- which(statVEC >= quantile(statVEC, 1 - quan))
    allSELECT <- cbind(allSELECT, SELECT)    
  }
  
  allREMOVE <- unique(matrix(allSELECT, ncol = 1))  
  
  if (plot) {
    plot(DATA[, 1], rep(1, length(DATA[, 1])), type = "n", ylim = c(0, 1), xaxt = "n", yaxt = "n", bty = "n", ...)  
    for (i in 1:ncol(allSELECT)) {
     points(allSELECT[, i], rep(1 - 0.05 * i, length(allSELECT[, i])), col = i, pch = 16, cex = 1, ...)                   
    }
    legend(round(max(DATA[, 1], na.rm = TRUE)/2), 0.5, names(coef(object))[colSEL], pch = 16, 
           col = 1:ncol(allSELECT), xjust = 0.5, horiz = TRUE, title = "coefficients", bty = "n", ...)
  }
  
  newOBJ <- update(object, data = DATA[-allREMOVE, ], fluo = 2)     
  return(newOBJ)      
} 