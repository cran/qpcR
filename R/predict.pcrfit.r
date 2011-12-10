predict.pcrfit <- function(
object, 
newdata,
which = c("y", "x"),
interval = c("none", "confidence", "prediction"),
level = 0.95,
...
)
{
  which <- match.arg(which)
  interval <- match.arg(interval)
    
  if (missing(newdata)) newDATA <- object$DATA
  else {
    if (which == "x") newDATA <- cbind(rep(1, nrow(newdata)), newdata)
    else newDATA <- newdata
  } 
  
  ### get predicted values
  if (which == "y") PRED <- object$MODEL$fct(newDATA[, 1], coef(object))
  else PRED <- object$MODEL$inv(newDATA[, 2], coef(object))
             
  ### make list with gradients     
  if (which == "y") DERIVS <- lapply(object$MODEL$parnames, function(x) D(object$MODEL$expr.grad, x))
  else DERIVS <- lapply(object$MODEL$parnames, function(x) D(object$MODEL$inv.grad, x))
   
  GRAD <- NULL
  resMAT <- NULL
  
  ### create t-statistic for confidence/prediction
  if (!identical(interval, "none")) {
    TQUAN <- qt(1 - (1 - level)/2, df.residual(object))     
  }       
      
  for (i in 1:nrow(newDATA)) {     
    ### create dataframe for gradient calculation
    tempDATA <- data.frame(newDATA[i, , drop = FALSE], t(coef(object))) 
           
    ### calculate gradients  
    dfEVAL <- as.numeric(lapply(DERIVS, function(x) eval(x, envir = tempDATA)))
    GRAD <- rbind(GRAD, as.numeric(dfEVAL))    
    
    ### calculate variance
    VAR <- dfEVAL %*% vcov(object) %*% dfEVAL           
  
    if (interval == "confidence") {     
      UPPER <- PRED[i] + TQUAN * sqrt(VAR)
      LOWER <- PRED[i] - TQUAN * sqrt(VAR)
      COLNAMES <- c("Prediction", "SE", "Lower", "Upper")         
    }  
  
    if (interval == "prediction") {
      UPPER <- PRED[i] + TQUAN * sqrt(VAR + resVar(object))        
      LOWER <- PRED[i] - TQUAN * sqrt(VAR + resVar(object))
      COLNAMES <- c("Prediction", "SE", "Lower", "Upper")  
    }
    
    if (interval == "none") {
      UPPER <- NULL
      LOWER <- NULL
      VAR <- NULL
      COLNAMES <- c("Prediction")  
    } 
    
    resMAT <-  rbind(resMAT, c(PRED[i], VAR, LOWER, UPPER))   
  }
  
  resMAT <- as.data.frame(resMAT)
  colnames(resMAT) <- COLNAMES      
   
  attr(resMAT, "gradient") <- as.matrix(GRAD)
  return(resMAT)   
}     