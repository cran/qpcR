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
  
  if (missing(newdata)) newdata <- object$DATA
  else {
    if (which == "x") newdata <- cbind(rep(1, nrow(newdata)), newdata)
  }     
  
  ### get predicted values
  if (which == "y") PRED <- object$MODEL$fct(newdata[, 1], coef(object))
  else PRED <- object$MODEL$inv(newdata[, 2], coef(object))
           
  ### make list with gradients     
  if (which == "y") DERIVS <- lapply(object$MODEL$parnames, function(x) D(object$MODEL$expr.grad, x))
       else DERIVS <- lapply(object$MODEL$parnames, function(x) D(object$MODEL$inv.grad, x))
   
  ### create dataframe for gradient calculation
  DATA <- cbind(newdata, t(coef(object)))   
    
  GRAD <- NULL
  retMat <- NULL
  
  ### create t-statistic for confidence/prediction
  if (!identical(interval, "none")) {
    tquan <- qt(1 - (1 - level)/2, df.residual(object))     
  }              
    
  for (i in 1:nrow(DATA)) {        
    ### calculate gradients     
    dfEval <- as.numeric(lapply(DERIVS, function(x) eval(x, envir = DATA[i, ])))      
    GRAD <- rbind(GRAD, as.numeric(dfEval))    
    ### calculate variance
    VAR <- dfEval %*% vcov(object) %*% dfEval           
  
    if (interval == "confidence") {     
      Upper <- PRED[i] + tquan * sqrt(VAR)
      Lower <- PRED[i] - tquan * sqrt(VAR)
      COLNAMES <- c("Prediction", "SE", "Lower", "Upper")         
    }  
  
    if (interval == "prediction") {
      Upper <- PRED[i] + tquan * sqrt(VAR + resVar(object))        
      Lower <- PRED[i] - tquan * sqrt(VAR + resVar(object))
      COLNAMES <- c("Prediction", "SE", "Lower", "Upper")  
    }
    
    if (interval == "none") {
      Upper <- NULL
      Lower <- NULL
      VAR <- NULL
      COLNAMES <- c("Prediction")  
    } 
    
    retMat <-  rbind(retMat, c(PRED[i], VAR, Lower, Upper))   
  }
  
  retMat <- as.data.frame(retMat)
  colnames(retMat) <- COLNAMES      
   
  attr(retMat, "gradient") <- as.matrix(GRAD)
  return(retMat)   
}     