pcrfit <- function(
data, 
cyc = 1, 
fluo, 
model = l4, 
do.optim = TRUE,
opt.method = "all", 
nls.method = "all", 
start = NULL,
robust = FALSE, 
weights = NULL,
verbose = TRUE,  
...)
{            
  require(minpack.lm, quietly = TRUE)       
  options(warn = -1)   

  ## version 1.3-4: create (stacked) data, depending on replicates
  if (length(fluo) == 1) {
    CYC <- data[, cyc]
    FLUO <- data[, fluo]
  } else {
    CYC <- rep(data[, cyc], length(fluo))
    FLUO <- stack(data[, fluo])[, 1]
    ssDATA <- rowMeans(data[, fluo], na.rm = TRUE)    
  }   
  
  ## define weights for nls and rnls 
  if (is.null(weights)) WEIGHTS <- rep(1, length(FLUO)) 
  else if (is.function(weights)) {  
    WFCT <- weights
    WEIGHTS <- WFCT(CYC, FLUO)      
  } else WEIGHTS <- weights   
  if (length(WEIGHTS) != length(FLUO)) stop("'weights' and 'fluo' have unequal length!")  
  
  ## eliminate NAs
  allDAT <- cbind(CYC, FLUO, WEIGHTS)
  tempDAT <- na.omit(allDAT)
  CYC <- tempDAT[, 1]
  FLUO <- tempDAT[, 2]
  WEIGHTS <- tempDAT[, 3]

  ## version 1.3-4: get selfStart values
  if (is.null(start)) {
    if (length(fluo) == 1) {
      ssVal <- model$ssFct(CYC, FLUO)
    } else {
      ssVal <- model$ssFct(data[, cyc], ssDATA)
    }
  } else ssVal <- start
  
  ## get attribute 'subset' transferred from ssFct
  ## (as is the case in mak2/mak3/mak3n model, when only curve
  ## up till second derivative max is taken)
  ## or mak3n/chag model with rescaling within [0, 1]
  SCALE <- attr(ssVal, "scale")
  if (!is.null(SCALE)) FLUO <- qpcR:::rescale(FLUO, SCALE[1], SCALE[2])
  
  SUB <- attr(ssVal, "subset")
  if (!is.null(SUB)) {  
    m <-which(CYC %in% SUB)
    CYC <- CYC[m]    
    FLUO <- FLUO[m]
    WEIGHTS <- WEIGHTS[m]
  }
    
  ## initialize parameter matrix
  ssValMat <- NULL
  ssValMat <- rbind(ssValMat, c("start", ssVal)) 
    
  ## objective function for 'optim', returns residual sum-of-squares
  FCT <- function(x) {     
    SSR <- sum(sqrt(WEIGHTS) * (FLUO - model$fct(CYC, x))^2)
    SSR       
  }    
  
  ## objective function for 'nls.lm', returns residuals
  FCT2 <- function(x) {     
    RESID <- FLUO - model$fct(CYC, x)       
    RESID <- sqrt(WEIGHTS) * RESID 
    RESID       
  }        
  
  if (do.optim) {
    ## define all methods
    if (opt.method == "all") opt.method <- c("LM", "BFGS", "Nelder-Mead", "SANN") 
    
    ## initial fit for refinement of starting values either by 'nls.lm' (default) or 'optim)
    for (i in opt.method) {
      if (!(i %in% c("LM", "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"))) 
        stop("Not an available 'opt.method'! Try one of 'LM', 'Nelder-Mead', 'BFGS', 'CG', 'L-BFGS-B', 'SANN'...")      
      if (i == "LM") {        
        OPTIM <- try(nls.lm(ssVal, FCT2, control = nls.lm.control(maxiter = 1000), ...), silent = TRUE)    
      }        
      if (i %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN")) {
        OPTIM <- try(optim(ssVal, FCT, method = i, hessian = TRUE, control = list(parscale = abs(ssVal), maxit = 1000), ...), silent = TRUE)
      }               
                              
      ## if no convergence is reached, try next method...
      if (inherits(OPTIM, "try-error")) {
        if (verbose) cat("Method '", i, "' did not converge. Trying next method...\n", sep = "")
        ## if last method was unsuccessful, terminate with error
        if (inherits(OPTIM, "try-error") && i == "SANN") stop("Used all available optimization methods, but none converged...")
        next
      }
      
      ## attach parameter values to matrix
      ssValMat <- rbind(ssValMat, c(i, OPTIM$par))     
      ssVal <- OPTIM$par       
          
      ## break out of loop if any method converged
      if (verbose) cat("Using method '", i, "' converged.\n", sep = "")      
      break      
    }
  }              
  
  names(ssVal) <- model$parnames    
        
  ## define all 'nls'' methods
  if (nls.method == "all") nls.method <- c("port", "default", "plinear") 
  ## added 1.3-4 only use 'plinear' for models of type mak2/mak3/mak3n/chag
  if (model$name %in% c("mak2", "mak3", "chag")) nls.method <- "plinear"
  
  ## coerce to dataframe
  DATA <- as.data.frame(cbind(Cycles = CYC, Fluo = FLUO))    
    
  ## try all methods successively
  for (j in nls.method) {
    if (!robust) NLS <- try(nls(as.formula(model$expr), data = DATA, start = as.list(ssVal), model = TRUE, 
                          algorithm = j, control = nls.control(maxiter = 1000, warnOnly = TRUE), weights = WEIGHTS, ...), silent = TRUE)
    else NLS <- try(qpcR:::rnls(as.formula(model$expr), data = DATA, start = as.list(ssVal), nls.method = j, 
                    control = nls.control(maxiter = 1000, warnOnly = TRUE), weights = WEIGHTS, ...), silent = TRUE)     
        
    ## if no convergence is reached, try next method...
    if (inherits(NLS, "try-error")) {
      if (verbose) cat("Method '", j, "' did not converge. Trying next method...\n", sep = "")
      if (inherits(NLS, "try-error") && j == "plinear") stop("Used all available 'nls' methods, but none converged...")
      next
    }
    
    ## break out of loop if any method converged
    if (verbose) cat("Using method '", j, "' converged.\n", sep = "")
    break
  }  
     
  ## attach parameter values to matrix
  ssValMat <- rbind(ssValMat, c(class(NLS), coef(NLS)))      
                
  ## modify 'object'
  NLS$DATA <- DATA     
  NLS$MODEL <- model  
  NLS$parMat <- ssValMat
  NLS$opt.method <- opt.method
  
  CALL <- as.list(NLS$call)     
  CALL$formula <- as.formula(model$expr)
  CALL$start <- ssVal
  NLS$call <- as.call(CALL)
  NLS$call2 <- match.call()
  assign("DATA", DATA, envir = globalenv())
  NLS$names <- names(data)[fluo]
  
  class(NLS) <- c("pcrfit", "nls")    
  return(NLS)      
}
