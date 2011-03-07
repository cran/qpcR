propagate <- function(
expr, 
data, 
type = c("raw", "stat"),  
do.sim = FALSE, 
use.cov = FALSE, 
nsim = 10000,
do.perm = FALSE, 
perm.crit = NULL, 
ties = NULL,
nperm = 2000,
alpha = 0.05,
plot = TRUE,    
logx = FALSE,
verbose = FALSE,      
...
)
{            
      require(MASS, quietly = TRUE)      
      type <- match.arg(type)
      if (!is.expression(expr) && !is.call(expr)) stop("'expr' must be an expression")
      crit.all <- c("perm > init", "perm == init", "perm < init")
      if (is.null(perm.crit)) perm.crit <- crit.all
        else if (!(perm.crit %in% crit.all)) stop("'perm.crit' must be one of 'perm > init', 'perm == init' or 'perm < init'!") 
      
      DATA <- as.matrix(data)
      EXPR <- expr 
      if (nrow(DATA) == 1) plot <- FALSE      
                 
      if (nsim > 0 && nsim < 5000) stop("Do at least 5000 simulations...")
      if (nsim < 0) stop("'nsim' must be >= 5000!")
      
      m <- match(all.vars(expr), colnames(data))
      if (any(is.na(m))) stop("Variable names of input dataframe and expression do not match!")
      if (length(unique(m)) != length(m)) stop("Some variable names are repetitive!")

      if (!is.logical(use.cov)) {
        if (!is.matrix(use.cov)) stop("'cov' must be a square covariance matrix!")
        if (is.matrix(use.cov)) {
          if (dim(use.cov)[1] != max(m) || dim(use.cov)[2] != max(m)) stop(paste("'use.cov' is not a ", max(m), "x", max(m), " matrix!", sep=""))
        }
      }
      
      if (!is.null(ties) && length(ties) != ncol(data)) stop("'ties' must have the same length as number of colums in 'data'!")

      if (type == "raw") {
        meanvals <- apply(DATA, 2, function(x) mean(x, na.rm = TRUE))
        sdvals <- apply(DATA, 2, function(x) sd(x, na.rm = TRUE))
      } else {
        meanvals <- DATA[1, ]
        sdvals <- DATA[2, ] 
      }                
      
      if (type == "raw" && use.cov == TRUE) SIGMA <- cov(DATA, use = "complete.obs")
      if (type == "raw" && use.cov == FALSE) SIGMA <- diag(diag(cov(DATA, use = "complete.obs")))
      if (type == "stat") SIGMA <- diag(DATA[2, ]^2) 
      
      if (all(!is.na(diag(SIGMA)))) {
        isCov <- TRUE
      } else {
        isCov <- FALSE
        SIGMA[is.na(SIGMA)] <- 0
      }        
      
      if (is.matrix(use.cov)) {
        m <- match(colnames(use.cov), colnames(DATA))            
        if (any(is.na(m))) stop("Names of input dataframe and var-cov matrix do not match!")             
        if (length(unique(m)) != length(m)) stop("Some names of the var-cov matrix are repetitive!")             
        if (is.unsorted(m)) stop("Names of input dataframe and var-cov matrix not in the same order!")             
        SIGMA <- use.cov
      }
      
      colnames(SIGMA) <- colnames(DATA)
      rownames(SIGMA) <- colnames(DATA)         
            
      ### Monte-Carlo simulation
      if (do.sim && isCov) {                    
        datSIM <- mvrnorm(nsim, mu = meanvals, Sigma = SIGMA, empirical = TRUE)
        colnames(datSIM) <- colnames(DATA)
        resSIM <- apply(datSIM, 1, function(x) eval(EXPR, envir = as.list(x)))
        confSIM <- quantile(resSIM, c(alpha/2, 1 - (alpha/2)), na.rm = TRUE)  
        if(do.sim && length(unique(resSIM)) == 1) print("Monte Carlo simulation gave unique repetitive error value! Are all derivatives constants?")   
        allSIM <- cbind(datSIM, resSIM)
      } else {
        resSIM <- datSIM <- confSIM <- allSIM <- NA
      }      
      
      ### permutation statistics with ties      
      if (do.perm) {
        if (is.null(ties)) ties <- 1:ncol(data)
        LEVELS <- unique(ties[!is.na(ties)])
        datPERM <- matrix(nrow = nperm, ncol = ncol(data))       
                
        for (i in 1:nrow(datPERM)) {                 
          SAMPLE <- sample(1:nrow(data), length(LEVELS), replace = TRUE)  
          for (j in 1:length(LEVELS)) {
            WHICH <- which(ties == LEVELS[j])              
            datPERM[i, WHICH] <- data[SAMPLE[j], WHICH]                       
           }          
        }  
                                   
        colnames(datPERM) <- colnames(data)                          
        resPERM <- apply(datPERM, 1, function(x) eval(EXPR, envir = as.list(x)))
        confPERM <- quantile(unique(resPERM), c(alpha/2, 1 - (alpha/2)), na.rm = TRUE)             
        
        ### permutation hypothesis testing
        datPERM2 <- datPERM                 
      
        for (i in 1:nrow(datPERM2)) {
          permLEVELS <- sample(LEVELS)
          for (j in 1:length(LEVELS)) {
            origPOS <- which(ties == LEVELS[j])            
            newPOS <- which(ties == permLEVELS[j])
            datPERM2[i, newPOS] <- datPERM[i, origPOS]
          }
        }            
        
        colnames(datPERM2) <- colnames(datPERM)
        resPERM2 <- apply(datPERM2, 1, function(x) eval(EXPR, envir = as.list(x)))   
        init <- resPERM 
        perm <- resPERM2           
        LOGIC <- lapply(perm.crit, function(x) eval(parse(text = x))) 
        names(LOGIC) <- perm.crit          
        pvalPERM <- lapply(LOGIC, function(x) sum(x == TRUE, na.rm = TRUE)/(length(x[!is.na(x)])))  
        names(pvalPERM) <- perm.crit                   
        datLOGIC <- as.data.frame(LOGIC)           
        allPERM <- cbind(datPERM, resPERM, datPERM2, resPERM2, datLOGIC)                                                 
      } else {
        datPERM <- resPERM <- resPERM2 <- confPERM <- pvalPERM <- allPERM <- NA
      }        
                
      ### error propagation        
      derivs <- try(lapply(colnames(DATA), D, expr = expr), silent = TRUE)
      if (inherits(derivs, "try-error")) stop(paste("Error within derivs:", derivs))      
      meanPROP <- eval(EXPR, envir = as.list(meanvals))
      NDERIVS <- sapply(derivs, eval, envir = as.list(meanvals))      
      errorPROP <- as.numeric(NDERIVS %*% SIGMA %*% matrix(NDERIVS))        
      confNORM <- abs(qnorm(alpha/2)) * sqrt(errorPROP)       
      confPROP <- c(meanPROP - confNORM, meanPROP + confNORM)
      distPROP <- rnorm(nsim, meanPROP, sqrt(errorPROP))       
                                       
      ### plotting simulations, permutations & propagations
      if (plot) {           
        par(mfrow = c(3, 1))
        par(mar = c(3, 1, 2, 1))   
        MAIN <- c("Monte-Carlo", "Permutation", "Error propagation")
                
        for (i in 1:3) {
          plotDATA <- switch(i, resSIM, resPERM, distPROP)
          confDATA <- switch(i, confSIM, confPERM, confPROP)  
          
          if (logx) {
            plotDATA <- suppressWarnings(log(plotDATA))  
            confDATA <- suppressWarnings(log(confDATA))          
          }     
                    
          FILTER <- quantile(plotDATA, c(0.01, 0.99), na.rm = TRUE)
          plotDATA <- plotDATA[plotDATA > FILTER[1] & plotDATA < FILTER[2]]
          
          if (length(plotDATA) <= 1) next() 
          if (!exists("XLIM")) XLIM <- range(plotDATA, na.rm = TRUE)         
                   
          HIST <- hist(plotDATA, xlab = "", ylab = "", col = "gray", yaxt = "n", breaks = 100, 
                       main = MAIN[i], xlim = XLIM, xaxt = "n", ...) 
          aT = axTicks(side = 1)          
          axis(side = 1, at = aT, labels = if(logx) round(exp(aT), 2) else aT)    
          suppressWarnings(rug(plotDATA))
          boxplot(plotDATA, horizontal = TRUE, add = TRUE, at = max(HIST$counts)/2, 
                  boxwex = diff(range(HIST$counts))/5, axes = FALSE, medcol = 2, boxfill = "gray", outline = FALSE, ...)
          rug(median(plotDATA, na.rm = TRUE), col = 2, lwd = 4, quiet = TRUE)
          abline(v = confDATA, lwd = 2, col = 4, ...)                              
        }
      }
            
      outDAT <- data.frame.na(Sim = mean(resSIM, na.rm = TRUE), Perm = mean(resPERM, na.rm = TRUE), Prop = meanPROP)          
      outDAT <- rbind.na(outDAT, c(sd(resSIM, na.rm = TRUE), sd(resPERM, na.rm = TRUE), sqrt(errorPROP)))
      outDAT <- rbind.na(outDAT, c(median(resSIM, na.rm = TRUE), median(resPERM, na.rm = TRUE)))
      outDAT <- rbind.na(outDAT, c(mad(resSIM, na.rm = TRUE), mad(resPERM, na.rm = TRUE)))    
      outDAT <- rbind.na(outDAT, c(confSIM[1], confPERM[1], confPROP[1]))
      outDAT <- rbind.na(outDAT, c(confSIM[2], confPERM[2], confPROP[2]))     
      row.names(outDAT) <- c("Mean", "Std.dev.", "Median", "MAD", "Conf.lower", "Conf.upper")       
           
      if (!all(is.na(pvalPERM))) {
        pVALS <- as.numeric(pvalPERM)         
        pMAT <- matrix(nrow = length(pVALS), ncol = ncol(outDAT)) 
        WHICH <- which(colnames(outDAT) == "Perm") 
        pMAT[, WHICH] <- pVALS          
        rownames(pMAT) <- names(pvalPERM)
        colnames(pMAT) <- colnames(outDAT)              
        outDAT <- rbind.na(outDAT, pMAT)
      }            
           
      if (verbose) OUT <- list(data.Sim = allSIM, data.Perm = allPERM, data.Prop = distPROP, derivs = derivs, covMat = SIGMA, summary = outDAT)
      else OUT <- list(summary = outDAT)
                   
      return(OUT)                                     
}