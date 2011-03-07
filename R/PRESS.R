PRESS <- function(object, verbose = TRUE)
{
  fetchDATA <- qpcR:::fetchData(object)
  DATA <- fetchDATA$data
  PRED.pos <- fetchDATA$pred.pos
  RESP.pos <- fetchDATA$resp.pos
  PRED.name <- fetchDATA$pred.name
  PRESS.res <- vector("numeric", nrow(DATA))
 
  for (i in 1:nrow(DATA)) {
    if (verbose) {
      qpcR:::counter(i)  
      flush.console()
    }
     
    newDATA <- DATA[-i, ]    
    
    if (class(object) == "pcrfit") newMOD <- pcrfit(newDATA, cyc = 1, fluo = 2, model = object$MODEL, verbose = FALSE) 
    else newMOD <- update(object, data = newDATA)
    
    newPRED <- as.data.frame(DATA[i, PRED.pos])
    colnames(newPRED) <- PRED.name
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
  
  
