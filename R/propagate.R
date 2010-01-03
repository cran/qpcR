propagate <- function(
expr, 
data, 
type = c("raw", "stat"),  
do.sim = FALSE, 
use.cov = FALSE, 
nsim = 10000,
do.perm = FALSE, 
perm.crit = "perm > init", 
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
        if(do.sim && length(unique(resSIM)) == 1) print("Monte Carlo simulation gave unique repetitive error value! Are all derivations constants?")   
        checkSIM <- cbind(datSIM, resSIM)
      } else {
        resSIM <- datSIM <- confSIM <- checkSIM <- NA
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
        evalPERMdata <- apply(datPERM, 1, function(x) eval(EXPR, envir = as.list(x)))
        confPERM <- quantile(unique(evalPERMdata), c(alpha/2, 1 - (alpha/2)), na.rm = TRUE)             
        
        ### permutation hypothesis testing
        datPERM2 <- datPERM3 <- datPERM                 
      
        for (i in 1:nrow(datPERM2)) {
          permLEVELS <- sample(LEVELS)
          for (j in 1:length(LEVELS)) {
            origPOS <- which(ties == LEVELS[j])            
            newPOS <- which(ties == permLEVELS[j])
            datPERM2[i, newPOS] <- datPERM3[i, origPOS]
          }
        } 
        colnames(datPERM2) <- colnames(datPERM)
        evalPERMsamp <- apply(datPERM2, 1, function(x) eval(EXPR, envir = as.list(x)))                             
        init <- evalPERMdata 
        perm <- evalPERMsamp     
        LOGIC <- lapply(perm.crit, function(x) eval(parse(text = x)))
        names(LOGIC) <- perm.crit        
        pvalPERM <- lapply(LOGIC, function(x) sum(x == TRUE, na.rm = TRUE)/length(x[!is.na(x)]))
        names(pvalPERM) <- perm.crit            
        datLOGIC <- as.data.frame(LOGIC)       
        checkPERM <- cbind(datPERM, evalPERMdata, datPERM2, evalPERMsamp, datLOGIC)                                              
      } else {
        evalPERMdata <- evalPERMsamp <- confPERM <- pvalPERM <- checkPERM <- NA
      }          
                
      ### error propagation        
      derivs <- try(lapply(colnames(DATA), D, expr = expr), silent = TRUE)
      if (inherits(derivs, "try-error")) stop(paste("Error within derivs:", derivs))      
      propMEAN <- eval(EXPR, envir = as.list(meanvals))
      NDERIVS <- sapply(derivs, eval, envir = as.list(meanvals))      
      propERROR <- as.numeric(NDERIVS %*% SIGMA %*% matrix(NDERIVS))        
      confNORM <- abs(qnorm(alpha/2)) * sqrt(propERROR)       
      confPROP <- c(propMEAN - confNORM, propMEAN + confNORM)       
                                       
      if (plot) {           
        op <- par(no.readonly = TRUE)
        on.exit(par(op))
        layout(matrix(c(1, 2, 3), 3, 1), heights = c(1, 1, 1))
        par(mai = c(0.5, 0.5, 0.25, 0.5))       
       
        ### plot simulations
        if (length(resSIM) > 1) {
          if (logx) resSIM2 <- suppressWarnings(log(resSIM)) else resSIM2 <- resSIM  
          if (logx) confSIM2 <- suppressWarnings(log(confSIM)) else confSIM2 <- confSIM          
          HIST <- hist(resSIM2, xlab = "", ylab = "", col = "gray", yaxt = "n", breaks = 100, 
                       main = "Monte-Carlo", axes = FALSE, ...)           
          if (logx) LABELS <- round(exp(HIST$breaks), 2) else LABELS <- round(HIST$breaks, 2)         
          axis(side = 1, at = HIST$breaks, labels = LABELS)                            
          rug(resSIM2)
          boxplot(resSIM2, horizontal = TRUE, add = TRUE, at = max(HIST$counts)/2, 
                  boxwex = diff(range(HIST$counts))/5, axes = FALSE, medcol = 2, boxfill = "gray", outline = FALSE, ...)
          rug(median(resSIM2, na.rm = TRUE), col = 2, lwd = 3, quiet = TRUE)
          abline(v = confSIM2, lwd = 2, col = 4, ...)
          aT <- HIST$breaks  
        } else aT <- NULL
        
        ### plot permutations
        if (length(evalPERMdata) > 1) {
          if (logx) evalPERMdata2 <- suppressWarnings(log(evalPERMdata)) else evalPERMdata2 <- evalPERMdata  
          if (logx) confPERM2 <- suppressWarnings(log(confPERM)) else confPERM2 <- confPERM           
          if (is.null(aT)) aT <- pretty(evalPERMdata2)                            
          HIST2 <- hist(evalPERMdata2, xlab = "", ylab = "", col = "gray", yaxt = "n", breaks = 100, 
                        main = "Permutation", xlim = if (!is.null(aT)) c(min(aT), max(aT)), axes = FALSE, ...)              
          if (logx) LABELS <- round(exp(aT), 2) else LABELS <- round(aT, 2)
          axis(side = 1, at = aT, labels = LABELS)
          rug(evalPERMdata2, quiet = TRUE)
          boxplot(evalPERMdata2, horizontal = TRUE, add = TRUE, at = max(HIST2$counts)/2, 
                  boxwex = diff(range(HIST2$counts))/5, axes = FALSE, medcol = 2, boxfill = "gray", outline = FALSE, ...)
          rug(median(evalPERMdata2, na.rm = TRUE), col = 2, lwd = 3, quiet = TRUE)          
          abline(v = confPERM2, lwd = 2, col = 4, ...)
        }
        
        ### plot error propagation
        DISTR <- rnorm(nsim, propMEAN, sqrt(propERROR))
        if (logx) DISTR <- suppressWarnings(log(DISTR)) else  DISTR <- DISTR  
        if (logx) confPROP2 <- suppressWarnings(log(confPROP)) else confPROP2 <- confPROP         
        if (is.null(aT)) aT <- pretty(DISTR)       
        HIST3 <- hist(DISTR, xlab = "", ylab = "", col = "gray", yaxt = "n", breaks = 100, 
                      main = "Error propagation", xlim = if (!is.null(aT)) c(min(aT), max(aT)), axes = FALSE, ...) 
        if (logx) LABELS <- round(exp(aT), 2) else LABELS <- round(aT, 2)
        axis(side = 1, at = aT, labels = LABELS)
        rug(DISTR, quiet = TRUE)
        boxplot(DISTR, horizontal = TRUE, add = TRUE, at = max(HIST3$counts)/2, 
                boxwex = diff(range(HIST3$counts))/5, axes = FALSE, medcol = 2, boxfill = "gray", outline = FALSE, ...)
        rug(median(DISTR, na.rm = TRUE), col = 2, lwd = 3, quiet = TRUE)               
        abline(v = confPROP2, lwd = 2, col = 4, ...)
      }       
      
      outDat <- data.frame(Sim = mean(resSIM, na.rm = TRUE), Perm = mean(evalPERMdata, na.rm = TRUE), Prop = propMEAN)          
      outDat <- rbind(outDat, c(sd(resSIM, na.rm = TRUE), sd(evalPERMdata, na.rm = TRUE), sqrt(propERROR)))
      outDat <- rbind(outDat, c(median(resSIM, na.rm = TRUE), median(evalPERMdata, na.rm = TRUE), NA))
      outDat <- rbind(outDat, c(mad(resSIM, na.rm = TRUE), mad(evalPERMdata, na.rm = TRUE), NA))    
      outDat <- rbind(outDat, c(confSIM[1], confPERM[1], confPROP[1]))
      outDat <- rbind(outDat, c(confSIM[2], confPERM[2], confPROP[2]))     
      row.names(outDat) <- c("Mean", "s.d.", "Median", "MAD", "Conf.lower", "Conf.upper")       
      
      if (!all(is.na(pvalPERM))) {
        pvals <- as.numeric(pvalPERM)       
        pMat <- matrix(nrow = length(pvals), ncol = ncol(outDat)) 
        WHICH <- which(names(outDat) == "Perm")
        pMat[, WHICH] <- pvals       
        rownames(pMat) <- names(pvalPERM) 
        colnames(pMat) <- colnames(outDat)       
        outDat <- rbind(outDat, pMat)
      }
      
      if (verbose) outDat <- list(summary = outDat, data.Sim = checkSIM, data.Perm = checkPERM, derivs = derivs, covMat = SIGMA)
      else outDat <- list(summary = outDat)
                   
      invisible(outDat)                                     
}