propagate2 <- function(
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
            
      ### error propagation (matrix notation)
      derivs <- try(lapply(colnames(DATA), D, expr = expr), silent = TRUE)
      if (inherits(derivs, "try-error")) stop(paste("Error within derivs:", derivs))      
      meanPROP <- eval(EXPR, envir = as.list(meanvals))
      NDERIVS <- sapply(derivs, eval, envir = as.list(meanvals))      
      errorPROP <- as.numeric(NDERIVS %*% SIGMA %*% matrix(NDERIVS))        
      confNORM <- abs(qnorm(alpha/2)) * sqrt(errorPROP)       
      confPROP <- c(meanPROP - confNORM, meanPROP + confNORM)
      distPROP <- rnorm(nsim, meanPROP, sqrt(errorPROP))
      
      ### error propagation (classical), first & second order Taylor expansion
      derivs <- try(lapply(colnames(DATA), D, expr = expr), silent = TRUE)
      if (inherits(derivs, "try-error")) stop(paste("Error within derivs:", derivs))      
      meanPROP <- eval(EXPR, envir = as.list(meanvals))
      NDERIVS <- sapply(derivs, eval, envir = as.list(meanvals))  
      
      VARNAMES <- all.vars(expr)
      names(NDERIVS) <- VARNAMES
      
      VAR <- diag(SIGMA)
      res <- sum(NDERIVS^2 * VAR)     
         
      COMBS <- expand.grid(VARNAMES, VARNAMES)
      SEL <- apply(COMBS, 1, function(x) x[1] != x[2])
      COMBS <- COMBS[SEL, ]     
      
      
      p1 <- as.numeric(NDERIVS[COMBS[, 1]])
      p2 <- as.numeric(NDERIVS[COMBS[, 2]])
      p3 <- as.numeric(diag(SIGMA[COMBS[, 1], COMBS[, 2]]))
      
      res2 <- res +  sum(p1 * p2 * p3)
          
              
      print(sqrt(res2))
               
              
}
