pcrfit <- function(
data, 
cyc = 1, 
fluo, 
model = l4, 
do.optim = TRUE,
opt.method = "LM", 
nls.method = "port", 
start = NULL,
robust = FALSE,
control = nls.control(),
weights = NULL,
...)
{            
  require(minpack.lm, quietly = TRUE) 
  require(minqa, quietly = TRUE)        
  options(warn = -1)

  Cycles <- data[, cyc]
  Fluo <- data[, fluo]
  DATA <- as.data.frame(cbind(Cycles, Fluo))
  compl <- complete.cases(Fluo)
  Cycles <- Cycles[compl]
  Fluo <- Fluo[compl]

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
    for (i in opt.method) {
      if (!(i %in% c("LM", "minqa", "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"))) 
        stop("Not an available 'opt.method'! Try one of 'LM', 'minqa', 'Nelder-Mead', 'BFGS', 'CG', 'L-BFGS-B', 'SANN'...")      
      if (i == "LM") {
        OPTIM <- nls.lm(ssVal, FCT2, ...)
      }
      if (i == "minqa") {
        OPTIM <- try(newuoa(ssVal, FCT, ...), silent = TRUE)  
      }
      if (i %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN")) {
        OPTIM <- try(optim(ssVal, FCT, method = i, hessian = TRUE, control = list(parscale = abs(ssVal), maxit = 50000), ...), silent = TRUE)
      }                           
      if (inherits(OPTIM, "try-error")) stop("Try to use other 'opt.method'")
      ssValMat <- rbind(ssValMat, c(i, OPTIM$par))     
      ssVal <- OPTIM$par       
    }                 
    if (!is.null(OPTIM$hessian)) EIGEN <- eigen(OPTIM$hessian)$values else EIGEN <- 0        
    if (any(EIGEN < 0)) cat(" Negative hessian eigenvalues! Consider a different 'opt.method'...")
  }
  
  names(ssVal) <- model$parnames    
  control$maxiter <-  50000
  control$warnOnly <-  TRUE

  if (!robust) NLS <- try(nls(as.formula(model$expr), data = DATA, start = as.list(ssVal), model = TRUE,
                          algorithm = nls.method, control = control, weights = weights, ...), silent = TRUE)
   else NLS <- try(qpcR:::rnls(as.formula(model$expr), data = DATA, start = as.list(ssVal), control = control, weights = weights, ...), silent = TRUE)
    
  if (inherits(NLS, "try-error")) stop("There was a problem during 'nls'. Try other method...")                             
    
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





