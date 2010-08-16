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
control = nls.control(),
weights = NULL,
verbose = TRUE,
...)
{            
  require(minpack.lm, quietly = TRUE) 
  require(minqa, quietly = TRUE)        
  options(warn = -1)

  Cycles <- data[, cyc]
  Fluo <- data[, fluo]
  DATA <- as.data.frame(cbind(Cycles, Fluo))  
  DATA <- DATA[complete.cases(DATA), ]
  Cycles <- DATA[, 1]
  Fluo <- DATA[, 2] 

  if (is.null(start)) ssVal <- model$ssFct(Cycles, Fluo)
   else ssVal <- start
    
  ssValMat <- NULL
  ssValMat <- rbind(ssValMat, c("start", ssVal))
  
  if (is.null(weights)) weights <- rep(1, length(Fluo)) else weights <- abs(weights)
  if (length(weights) != length(Fluo)) stop("'weights' and 'fluo' have unequal length!")
    
  FCT <- function(x) {     
    SSR <- sum(sqrt(weights) * (Fluo - model$fct(Cycles, x))^2)
    SSR       
  }    
  
  FCT2 <- function(x) {     
    RESID <- Fluo - model$fct(Cycles, x)       
    RESID <- sqrt(weights) * RESID 
    RESID       
  }                         
  
  if (do.optim) {
    if (opt.method == "all") opt.method <- c("LM", "BFGS", "Nelder-Mead", "minqa", "SANN") 
    for (i in opt.method) {
      if (!(i %in% c("LM", "minqa", "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"))) 
        stop("Not an available 'opt.method'! Try one of 'LM', 'minqa', 'Nelder-Mead', 'BFGS', 'CG', 'L-BFGS-B', 'SANN'...")      
      if (i == "LM") {
        OPTIM <- try(nls.lm(ssVal, FCT2, ...), silent = TRUE)
      }
      if (i == "minqa") {
        OPTIM <- try(newuoa(ssVal, FCT, ...), silent = TRUE)  
      }
      if (i %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN")) {
        OPTIM <- try(optim(ssVal, FCT, method = i, hessian = TRUE, control = list(parscale = abs(ssVal), maxit = 50000), ...), silent = TRUE)
      }                           
      if (inherits(OPTIM, "try-error")) {
        if (verbose) cat("Method '", i, "' did not converge. Trying next method...\n", sep = "")
        if (inherits(OPTIM, "try-error") && i == "SANN") stop("Used all available optimization methods, but none converged...")
        next
      }
      ssValMat <- rbind(ssValMat, c(i, OPTIM$par))     
      ssVal <- OPTIM$par
      if (!is.null(OPTIM$hessian)) EIGEN <- eigen(OPTIM$hessian)$values else EIGEN <- 0        
      if (any(EIGEN < 0)) {
        if (verbose) cat("Negative hessian eigenvalues! Trying next method...\n")
        next
      }  
      if (verbose) cat("Using method '", i, "' converged.\n", sep = "")
      break      
    }
  }              
  
  names(ssVal) <- model$parnames    
  CONTROL <- nls.control()
  CONTROL$maxiter <-  50000
  CONTROL$warnOnly <-  TRUE
  CONTROL$tol <- 1e-5
  
  if (nls.method == "all") nls.method <- c("port", "default", "plinear") 
  
  for (j in nls.method) {
    if (!robust) NLS <- try(nls(as.formula(model$expr), data = DATA, start = as.list(ssVal), model = TRUE,
                          algorithm = j, control = CONTROL, weights = weights, ...), silent = TRUE)
    else NLS <- try(qpcR:::rnls(as.formula(model$expr), data = DATA, start = as.list(ssVal), nls.method = j, 
                    control = CONTROL, weights = weights, ...), silent = TRUE)
    if (inherits(NLS, "try-error")) {
      if (verbose) cat("Method '", j, "' did not converge. Trying next method...\n", sep = "")
      if (inherits(NLS, "try-error") && j == "plinear") stop("Used all available 'nls' methods, but none converged...")
      next
    }
    if (verbose) cat("Using method '", j, "' converged.\n", sep = "")
    break
  }  
   
  ssValMat <- rbind(ssValMat, c(class(NLS), coef(NLS)))      
                
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