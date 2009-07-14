propagate <- function(
expr, 
data, 
type = c("raw", "stat"),  
do.sim = TRUE, 
use.cov = FALSE, 
nsim = 10000,
do.perm = FALSE, 
perm.crit = "m1 < m2", 
ties = NULL,
nperm = 2000,
alpha = 0.05,
plot = TRUE, 
xlim = NULL,      
...
)
{
      require(MASS, quietly = TRUE)      
      type <- match.arg(type)
      if (!is.expression(expr) && !is.call(expr)) stop("'expr' must be an expression")
      DATA <- as.matrix(data)
      EXPR <- expr             
                 
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
      if (do.sim) {                    
        simDat <- mvrnorm(nsim, mu = meanvals, Sigma = SIGMA, empirical = TRUE)
        colnames(simDat) <- colnames(DATA)
        simRes <- apply(simDat, 1, function(x) eval(EXPR, envir = as.list(x)))
        CONFsim <- quantile(simRes, c(alpha/2, 1 - (alpha/2)), na.rm = TRUE)  
        if(do.sim && length(unique(simRes)) == 1) print("Monte Carlo simulation gave unique repetitive error value! Are all derivations constants?")   
      } else {
        simRes <- NA 
        simDat <- NA
        CONFsim <- NA
      }        

      ### statistics of permutation with replacement      
      if (do.perm) {
        if (is.null(ties)) ties <- 1:ncol(data)
        nlev <- unique(ties)
        p.count <- 0              
        all.count <- 0
        EVALperms <- NULL           
        for (j in 1:nperm) {
          DATAperm1 <- data
          DATAperm2 <- data
          DATAperm1[!is.na(DATAperm1)] <- NA  
          DATAperm2[!is.na(DATAperm2)] <- NA                            
          
          for (i in nlev) {              
            if (!is.na(i)) WHICH <- which(ties == i) else WHICH <- which(is.na(ties)) 
            DATAtemp <- as.data.frame(data[, WHICH])                      
            DATAtemp <- apply(DATAtemp, 2, function(x) sample(x, replace = TRUE))  ### resampled within columns
            DATAperm1[, WHICH] <- DATAtemp               
            if (!is.na(i)) DATAperm2[, WHICH] <- sample(data[, WHICH], replace = TRUE)    ### resamples between columns   
            else DATAperm2[, WHICH] <- DATAtemp ### but not if NA's (sample within)
          }                 
          
          EVALperm1 <- apply(DATAperm1, 1, function(x) eval(EXPR, envir = as.list(x)))           
          EVALperm2 <- apply(DATAperm2, 1, function(x) eval(EXPR, envir = as.list(x)))         
          m1 <- mean(EVALperm1, na.rm = TRUE)
          m2 <- mean(EVALperm2, na.rm = TRUE)             
          if  (eval(parse(text = perm.crit))) p.count <- p.count + 1    
          all.count <- all.count + 1           
          EVALperms <- c(EVALperms, EVALperm1)            
        }
        CONFperm <- quantile(EVALperms, c(alpha/2, 1 - (alpha/2)), na.rm = TRUE)         
        PVALperm <- 1-(p.count/all.count)                                 
      } else {
        EVALperms <- NA 
        CONFperm <- NA
        PVALperm <- NA
      }       
      
      ### error propagation
      derivs <- try(lapply(colnames(DATA), D, expr = expr), silent = TRUE)
      if (inherits(derivs, "try-error")) stop(paste("Error within derivs:", derivs))      
      MEANprop <- eval(EXPR, envir = as.list(meanvals))
      NDERIVS <- sapply(derivs, eval, envir = as.list(meanvals))      
      ERRORprop <- as.numeric(NDERIVS %*% SIGMA %*% matrix(NDERIVS))      
      CONFnorm <- abs(qnorm(alpha/2)) * sqrt(ERRORprop)
      CONFprop <- c(MEANprop - CONFnorm, MEANprop + CONFnorm)      
                                       
      if (plot) {   
        op <- par(no.readonly = TRUE)
        on.exit(par(op))
        layout(matrix(c(1, 2, 3), 3, 1), heights = c(1, 1, 1))
        par(mai = c(0.5, 0.5, 0.25, 0.5))  
        
        if (missing(xlim)) RANGE <- quantile(c(simRes, EVALperms), c(0.01, 0.99), na.rm = TRUE) else RANGE <- xlim
                    
        if (length(simRes) > 1) {
          HIST <- hist(simRes, xlab = "", ylab = "", col = "gray", yaxt = "n", breaks = 100, main = "Monte-Carlo", xlim = RANGE, ...)
          rug(simRes)
          boxplot(simRes, horizontal = TRUE, add = TRUE, at = max(HIST$counts)/2, 
                  boxwex = diff(range(HIST$counts))/5, axes = FALSE, medcol = 2, boxfill = "gray", ...)
          rug(median(simRes, na.rm = TRUE), col = 2, lwd = 3, quiet = TRUE)
          abline(v = CONFsim, lwd = 2, col = 4, ...)
        }
        
        if (length(EVALperms) > 1) {
          HIST2 <- hist(EVALperms, xlab = "", ylab = "", col = "gray", yaxt = "n", breaks = 100, main = "Permutation", xlim = RANGE, ...)
          rug(EVALperms, quiet = TRUE)
          boxplot(EVALperms, horizontal = TRUE, add = TRUE, at = max(HIST2$counts)/2, 
                  boxwex = diff(range(HIST2$counts))/5, axes = FALSE, medcol = 2, boxfill = "gray", ...)
          rug(median(EVALperms, na.rm = TRUE), col = 2, lwd = 3, quiet = TRUE)         
          abline(v = CONFperm, lwd = 2, col = 4, ...)
        }
        
        DISTR <- rnorm(10000, MEANprop, sqrt(ERRORprop))
        if (is.null(RANGE)) RANGE <- c(min(DISTR, na.rm = TRUE), max(DISTR, na.rm = TRUE))
        HIST3 <- hist(DISTR, xlab = "", ylab = "", col = "gray", yaxt = "n", breaks = 100, main = "Error propagation", xlim = RANGE, ...) 
        rug(DISTR, quiet = TRUE)
        boxplot(DISTR, horizontal = TRUE, add = TRUE, at = max(HIST3$counts)/2, 
                boxwex = diff(range(HIST3$counts))/5, axes = FALSE, medcol = 2, boxfill = "gray", ...)
        rug(median(DISTR, na.rm = TRUE), col = 2, lwd = 3, quiet = TRUE)       
        abline(v = CONFprop, lwd = 2, col = 4, ...)
      }       
      
      invisible(list(mean.Sim = mean(simRes, na.rm = TRUE),
                     sd.Sim = sd(simRes, na.rm = TRUE), 
                     med.Sim = median(simRes, na.rm = TRUE), 
                     mad.Sim = mad(simRes, na.rm = TRUE),                       
                     data.Sim = simDat,
                     eval.Sim = simRes,                      
                     conf.Sim = CONFsim,
                     mean.Perm = mean(EVALperms, na.rm = TRUE),
                     sd.Perm = sd(EVALperms, na.rm = TRUE),
                     med.Perm = median(EVALperms, na.rm = TRUE),
                     mad.Perm = mad(EVALperms, na.rm = TRUE),
                     conf.Perm = CONFperm,
                     pval.Perm = PVALperm,
                     eval.Prop = MEANprop,
                     error.Prop = sqrt(ERRORprop),
                     deriv.Prop = derivs,
                     conf.Prop = CONFprop,
                     covMat = SIGMA))                      
}