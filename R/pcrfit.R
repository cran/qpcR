pcrfit <- function(
data, 
cyc = 1, 
fluo, 
model = l4, 
start = NULL,
offset = 0, 
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
  
  ## define weights for nls
  if (is.null(weights)) WEIGHTS <- rep(1, length(FLUO))
  ## version 1.3-7: use an expression for weights
  if (is.character(weights)) WEIGHTS <- wfct(weights, CYC, FLUO, model = model,
                                             start = start, offset = offset, verbose = TRUE)
  if (is.numeric(weights)) WEIGHTS <- weights
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
  ## (as is the case in mak2/mak3 model, when only curve
  ## up till second derivative max is taken)  
  ## version 1.3-4: offset parameter from SDM of curve
  SUB <- attr(ssVal, "subset")
  if (!is.null(SUB)) {  
    if (offset < 0) SUB <- head(SUB, offset)
    else if (offset > 0) SUB <- c(SUB, max(SUB) + (1:offset))
    m <-which(CYC %in% SUB)
    CYC <- CYC[m]    
    FLUO <- FLUO[m]
    WEIGHTS <- WEIGHTS[m]
  }
  
  ## initialize parameter matrix
  ssValMat <- NULL
  ssValMat <- rbind(ssValMat, c("start", ssVal)) 
    
  names(ssVal) <- model$parnames      
        
  ## coerce to dataframe
  DATA <- as.data.frame(cbind(Cycles = CYC, Fluo = FLUO))  
      
  ## make nlsModel using 'nlsLM' from package 'minpack.lm'
  NLS <- nlsLM(as.formula(model$expr), data = DATA, start = as.list(ssVal), model = TRUE, 
              algorithm = "LM", control = nls.lm.control(maxiter = 1000, maxfev = 10000), weights = WEIGHTS, ...)
     
  ## attach parameter values to matrix
  ssValMat <- rbind(ssValMat, c(class(NLS), coef(NLS)))      
                
  ## modify 'object'
  NLS$DATA <- DATA     
  NLS$MODEL <- model  
  NLS$parMat <- ssValMat
  NLS$opt.method <- "LM"
  
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
