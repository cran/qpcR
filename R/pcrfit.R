pcrfit <- function(
data, 
cyc = 1, 
fluo, 
model = l4, 
opt.method = rep("Nelder", 5), 
nls.method = "port", 
...)
{
  require(rgenoud, quietly = TRUE)     
  options(warn = -1)             
  
  Cycles <- data[, cyc]          
  Fluo <- data[, fluo]        
    
  DATA <- as.data.frame(cbind(Cycles, Fluo))     
  ssVal <- model$ssFct(Cycles, Fluo)   
    
  ssValMat <- NULL
  ssValMat <- rbind(ssValMat, c("start", ssVal))
  
  FCT <- function(x) {
    SSR <- sum((Fluo - model$fct(Cycles, x))^2)
    SSR
  }         
  
  FCT2 <- function(x) {
    RES <- Fluo - model$fct(Cycles, x)
    RES
  }        
  
  for (i in opt.method) {
    if (i != "GA") {
      OPTIM <- try(optim(ssVal, FCT, method = i, hessian = TRUE, control = list(parscale = abs(ssVal), maxit = 50000), ...), silent = TRUE)
    } else {
      OPTIM <- try(genoud(fn = FCT, nvars = length(model$parnames), starting.values = ssVal, 
                      hessian = TRUE, print.level = 0, pop.size = 100, ...), silent = TRUE)         
    }    
    if (inherits(OPTIM, "try-error")) stop("Try to use other 'opt.method'")
    ssValMat <- rbind(ssValMat, c(i, OPTIM$par))     
    ssVal <- OPTIM$par       
  }

  EIGEN <- eigen(OPTIM$hessian)$values  
  if (any(EIGEN < 0)) print("One of the hessian matrix eigenvalues is negative! Consider a different 'opt.method'...")
  
  names(ssVal) <- model$parnames

  NLS <- try(nls(as.formula(model$expr), data = DATA, start = as.list(ssVal), model = TRUE,
                    algorithm = nls.method, control = list(maxiter = 50000, warnOnly = TRUE), ...), silent = TRUE)
    
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
  assign("DATA", DATA, envir = globalenv())

  class(NLS) <- c("pcrfit", "nls")
  return(NLS)      
} 


