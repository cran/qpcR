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
PRESS = FALSE,
...
)
{
    if (is.null(par)) stop("Please supply parameter estimates!")
    if (length(par) != length(model$parnames)) stop("Length of 'par' does not match number of parameters in ", model$name)

    ### make initial model
    initMOD <- model$fct(cyc, par)
    fluoMAT <- matrix(nrow = nsim, ncol = length(initMOD))
    lenFM <- length(initMOD)
    
    ### add random noise
    for (i in 1:nsim) {
     ranVEC <- sapply(initMOD, function(x) rnorm(1, mean = x, sd = error * errfun(x)))
     fluoMAT[i, ] <- ranVEC
    }

    ### make empty plot
    if (plot) {
      plot(cyc, initMOD, type = "n", ylim = c(min(fluoMAT, na.rm = TRUE), max(fluoMAT, na.rm = TRUE)),
           lwd = 2, col = 1, xlab = "Cycles", ylab = "Fluo", ...)
      apply(fluoMAT, 1, function(x) points(cyc, x, cex = 0.5, ...))
    }        
    
    ### select models to fit
    if (is.null(fitmodel)) fitMODEL <- list(model) else fitMODEL <- as.list(fitmodel)
       
    ### create data matrix to be fitted
    fluoMAT <- t(fluoMAT)
    fluoMAT <- cbind(cyc, fluoMAT)
    colnames(fluoMAT) <- c("Cycles", paste("Fluo.", 1:nsim, sep = ""))
    
    ### preallocate lists for increased speed
    coefLIST <- vector("list", length = length(fitMODEL))
    names(coefLIST) <- sapply(fitMODEL, function(x) x$name)
    gofLIST <- vector("list", length = length(fitMODEL))
    names(gofLIST) <- sapply(fitMODEL, function(x) x$name)

    colVEC <- rainbow(length(fitMODEL))
    
    ### for all models do...
    for (k in 1:length(fitMODEL)) {
      cat(fitMODEL[[k]]$name, "\n")

      ### preallocate matrix dimensions for increased speed
      coefMAT <- matrix(nrow = length(fitMODEL[[k]]$parnames), ncol = nsim)
      gofMAT <- matrix(nrow = ifelse(PRESS, 10, 9), ncol = nsim)
      
      ### for all simulated data do...
      for (i in 1:nsim) {
        FIT <- pcrfit(fluoMAT, 1, i + 1, fitMODEL[[k]], verbose = FALSE, ...)
        qpcR:::counter(i)
        
        ### plot fitted curve
        lines(cyc, fitted(FIT), col = colVEC[k], ...)

        ### put coef's in matrix
        if (i == 1) rownames(coefMAT) <- names(coef(FIT))
        coefMAT[, i] <- coef(FIT)

        ### obtain GOF measures
        GOFs <- unlist(pcrGOF(FIT))

        ### obtain reduced chi-square with error taken from
        ### simulation
        CHISQ <- fitchisq(FIT, error = error)$chi2.red
        GOFs <- c(GOFs, chi2.red = CHISQ)
        
        ### do PRESS, if selected and obtain P-square
        if (PRESS) {
           DATA <- fluoMAT[, c(1, i + 1)]
           MODEL <- fitMODEL[[k]]
           assign("MODEL", MODEL, envir = .GlobalEnv)
           FIT2 <- pcrfit(DATA, 1, 2, MODEL, verbose = FALSE)
           Psq <- PRESS(FIT2, verbose = FALSE)$P.square
           GOFs <- c(GOFs, P.square = Psq)
        }

        ### put GOF's in matrix
        if (i == 1) rownames(gofMAT) <- names(GOFs)
        gofMAT[, i] <- GOFs
      }
      
      cat("\n\n")
      coefLIST[[k]] <- coefMAT
      gofLIST[[k]] <- gofMAT
    }
    
    modMAT <- NULL
    colMAT <- NULL
    
    ### in case of model selection for all GOF measures...
    if (select) {
     gofSEL <- NULL
     selMAT <- NULL
          
     for (i in 1:nrow(gofLIST[[1]])) {
      for (j in 1:length(gofLIST)) {
        gofSEL <- rbind(gofSEL, gofLIST[[j]][i, ])
      }   
         
      if (i == 1) SEL <- apply(gofSEL, 2, function(x) which.max(x))
      if (i == 2) SEL <- apply(gofSEL, 2, function(x) which.max(x))
      if (i == 3) SEL <- apply(gofSEL, 2, function(x) which.min(x))
      if (i == 4) SEL <- apply(gofSEL, 2, function(x) which.min(x))
      if (i == 5) SEL <- apply(gofSEL, 2, function(x) which.min(x))
      if (i == 6) SEL <- apply(gofSEL, 2, function(x) which.min(x))
      if (i == 7) SEL <- apply(gofSEL, 2, function(x) which.min(x))
      if (i == 8) SEL <- apply(gofSEL, 2, function(x) which.min(x))
      if (i == 9) SEL <- apply(gofSEL, 2, function(x) which.min(abs(1 - x)))
      if (i == 10) SEL <- apply(gofSEL, 2, function(x) which.max(x))
      
      modSEL <- sapply(SEL, function(x) fitMODEL[[x]]$name)
      modMAT <- rbind(modMAT, modSEL)
      colMAT <- rbind(colMAT, SEL)
      gofSEL <- NULL
     }       
     rownames(modMAT) <- rownames(gofLIST[[1]])
    }  
    
    STAT <- lapply(gofLIST, function(x) apply(x, 1, statfun))
        
    OUT <-  list(cyc = cyc, fluoMat = t(fluoMAT), coefList = coefLIST,
                 gofList = gofLIST, statList = STAT, modelMat = modMAT)

    class(OUT) <- "pcrsim"
    return(OUT)
}     
