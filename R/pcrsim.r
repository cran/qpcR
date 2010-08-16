pcrsim <- function(
cyc = 1:30,
model = l4,
par = NULL,
nsim = 10,        
error = 0.02,
errfun = function(y) 1,
plot = TRUE,
fitmodel = NULL,
select = FALSE,
statfun = function(y) mean(y, na.rm = TRUE), 
...
)
{
    if (is.null(par)) stop("Please supply parameter estimates!")
    if (length(par) != length(model$parnames)) stop("Length of 'par' does not match number of parameters in ", model$name)

    initMod <- model$fct(cyc, par)
    fluoMat <- matrix(nrow = nsim, ncol = length(initMod))
    lenFM <- length(initMod)             
    
    for (i in 1:nsim) {
     ranVec <- sapply(initMod, function(x) rnorm(1, mean = x, sd = error * errfun(x)))       
     fluoMat[i, ] <- ranVec
    }

    if (plot) {
      plot(cyc, initMod, type = "n", ylim = c(min(fluoMat, na.rm = TRUE), max(fluoMat, na.rm = TRUE)),
           lwd = 2, col = 1, xlab = "Cycles", ylab = "Fluo", ...)
      apply(fluoMat, 1, function(x) points(cyc, x, cex = 0.5, ...))
    }        
    
    if (is.null(fitmodel)) fitmodel <- list(model) else fitmodel <- as.list(fitmodel)       
       
    fluoMatT <- t(fluoMat)
    fluoMatT <- cbind(cyc, fluoMatT)
    colnames(fluoMatT) <- c("Cycles", rep("Fluo", ncol(fluoMatT) - 1))
      
    coefList <- list()
    coefMat <- NULL
            
    gofList <- list()
    gofMat <- NULL            
    
    colvec <- rainbow(length(fitmodel))
    
    for (k in 1:length(fitmodel)) {
      cat(fitmodel[[k]]$name, "\n")         
      for (i in 2:ncol(fluoMatT)) {
        FIT <- pcrfit(fluoMatT, 1, i, fitmodel[[k]], verbose = FALSE, ...)          
        j <- i - 1 
        qpcR:::counter(j)        
        if (plot) lines(cyc, fitted(FIT), col = colvec[k], ...)  
        coefMat <- cbind(coefMat, coef(FIT))
        gofMat <- cbind(gofMat, unlist(pcrGOF(FIT, error = error)))                     
      }
      cat("\n\n")
      coefList[[k]] <- coefMat
      gofList[[k]] <- gofMat                  
      coefMat <- NULL
      gofMat <- NULL
    }          
    
    modMat <- NULL
    
    if (select) {  
     gofSel <- NULL
     selMat <- NULL 
          
     for (i in 1:nrow(gofList[[1]])) {
      for (j in 1:length(gofList)) {
        gofSel <- rbind(gofSel, gofList[[j]][i, ])
      }   
         
      if (i == 1) SEL <- apply(gofSel, 2, function(x) which.max(x))
      if (i == 2) SEL <- apply(gofSel, 2, function(x) which.max(x))
      if (i == 3) SEL <- apply(gofSel, 2, function(x) which.min(x))
      if (i == 4) SEL <- apply(gofSel, 2, function(x) which.min(x))
      if (i == 5) SEL <- apply(gofSel, 2, function(x) which.min(x))
      if (i == 6) SEL <- apply(gofSel, 2, function(x) which.min(x))
      if (i == 7) SEL <- apply(gofSel, 2, function(x) which.min(x))
      if (i == 8) SEL <- apply(gofSel, 2, function(x) which.min(abs(1 - x)))       
      
      modSel <- sapply(SEL, function(x) fitmodel[[x]]$name) 
      modMat <- rbind(modMat, modSel)     
      gofSel <- NULL                
     }       
     rownames(modMat) <- rownames(gofList[[1]])                            
    }  
    
    STAT <- lapply(gofList, function(x) apply(x, 1, statfun)) 
        
    invisible(list(cyc = cyc, fluoMat = t(fluoMat), coefList = coefList, 
                   gofList = gofList, statList = STAT, modelMat = modMat))
}     