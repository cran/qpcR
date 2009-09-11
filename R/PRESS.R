PRESS <- function(object, verbose = TRUE)
{
  if (!is.null(object$call$data)) DATA <- eval(object$call$data)
    else DATA <- as.data.frame(sapply(all.vars(object$call$formula), function(a) get(a, envir = .GlobalEnv)))
 
  VARS <- all.vars(object$call$formula)
  LHS <- VARS[1]
  RHS <- VARS[-1]   
  matchPRED <- which(!is.na(match(RHS, colnames(DATA)))) 
  PREDname <- RHS[matchPRED]
  PRED.pos <- which(colnames(DATA) == PREDname)
  RESP.pos <- which(colnames(DATA) == LHS)
  PRESS.res <- NULL
 
  for (i in 1:nrow(DATA)) {
    if (verbose) {
      if (i %% 10 == 0) cat(i) else cat(".")
      if (i %% 50 == 0) cat("\n")
      flush.console()
    }
  
    newDATA <- DATA[-i, ]      
    newMOD <- update(object, data = newDATA)        
    newPRED <- as.data.frame(DATA[i, PRED.pos])
    colnames(newPRED) <- PREDname
    y.hat <- as.numeric(predict(newMOD, newdata = newPRED))
    PRESS.res[i] <- DATA[i, RESP.pos] - y.hat
  }
  cat("\n")
  Yi <- residuals(object) - fitted(object)
  TSS <- sum((Yi - mean(Yi))^2)
  RSS <- sum(PRESS.res^2)
  P.square <- 1 - (RSS/TSS)    

  return(list(stat = sum(PRESS.res^2), residuals = PRESS.res, P.square = P.square))
}
  
  