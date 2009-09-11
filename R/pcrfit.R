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
  require(rgenoud, quietly = TRUE)
  require(minpack.lm, quietly = TRUE)     
  options(warn = -1)    
  if (is.null(weights)) wts <- rep(1, nrow(data)) else wts <- weights/max(weights, na.rm = TRUE)      
  
  Cycles <- data[, cyc]          
  Fluo <- data[, fluo]        
    
  DATA <- as.data.frame(cbind(Cycles, Fluo))     
  if (is.null(start)) ssVal <- model$ssFct(Cycles, Fluo)
   else ssVal <- start      
    
  ssValMat <- NULL
  ssValMat <- rbind(ssValMat, c("start", ssVal))
  
  FCT <- function(x) {     
    SSR <- sum(wts * ((Fluo - model$fct(Cycles, x))^2))    
    SSR       
  }    
  
  FCT2 <- function(x) {     
    RESID <- wts * (Fluo - model$fct(Cycles, x))      
    RESID       
  }                         
  
  if (do.optim) {
    for (i in opt.method) {      
      if (i == "LM") {
        OPTIM <- nls.lm(ssVal, FCT2, ...)              
      }
      if (i != "GA" && i != "LM") {
        OPTIM <- try(optim(ssVal, FCT, method = i, hessian = TRUE, control = list(parscale = abs(ssVal), maxit = 50000), ...), silent = TRUE)
      } 
      if (i == "GA") {
        OPTIM <- try(genoud(fn = FCT, nvars = length(model$parnames), starting.values = ssVal, 
                      hessian = TRUE, print.level = 0, pop.size = 100, ...), silent = TRUE)         
      }                       
      if (inherits(OPTIM, "try-error")) stop("Try to use other 'opt.method'")
      ssValMat <- rbind(ssValMat, c(i, OPTIM$par))     
      ssVal <- OPTIM$par       
    }                 
    EIGEN <- eigen(OPTIM$hessian)$values  
    if (any(EIGEN < 0)) print("One of the hessian matrix eigenvalues is negative! Consider a different 'opt.method'...")
  }
  
  names(ssVal) <- model$parnames    
  control$maxiter <-  50000
  control$warnOnly <-  TRUE
  

  if (!robust) NLS <- try(nls(as.formula(model$expr), data = DATA, start = as.list(ssVal), model = TRUE,
                          algorithm = nls.method, control = control, weights = wts, ...), silent = TRUE)
   else NLS <- try(qpcR:::rnls(as.formula(model$expr), data = DATA, start = as.list(ssVal), control = control, weights = wts, ...), silent = TRUE)
    
  if (inherits(NLS, "try-error")) stop("There was a problem during 'nls'. Try other method...")                             
    
  ssValMat <- rbind(ssValMat, c(class(NLS), coef(NLS)))      
                
  NLS$DATA <- DATA     
  NLS$MODEL <- model  
  NLS$parMat <- ssValMat
  NLS$opt.method <- opt.method
  NLS$weights <- wts   
  
  CALL <- as.list(NLS$call)     
  CALL$formula <- as.formula(model$expr)
  CALL$start <- ssVal
  NLS$call <- as.call(CALL)
  NLS$call2 <- match.call()
  assign("DATA", DATA, envir = globalenv())

  class(NLS) <- c("pcrfit", "nls")
  return(NLS)      
}





