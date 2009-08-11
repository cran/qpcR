PRESS <- function(object)
{
  if (!is.null(object$call$data)) DATA <- eval(object$call$data)
    else DATA <- as.data.frame(sapply(all.vars(object$call$formula), function(a) get(a, envir = .GlobalEnv)))
  if (!is.na(class(object)[2]) && class(object)[2] == "pcrfit") DATA <- object$DATA
 
  CALL <- as.list(object$call)
  VARS <- all.vars(object$call$formula)
  LHS <- VARS[1]
  RHS <- VARS[-1]
  matchPRED <- which(!is.na(match(RHS, colnames(DATA))))
  PREDname <- RHS[matchPRED]
  PRED.pos <- which(colnames(DATA) == PREDname)
  RESP.pos <- which(colnames(DATA) == LHS)
  PRESS.res <- NULL
 
  for (i in 1:nrow(DATA)) {
    if (i %% 10 == 0) cat(i) else cat(".")
    if (i %% 50 == 0) cat("\n")
    flush.console()
  
    NEWCALL <- CALL
    NEWCALL$data <- DATA[-i, ]
    
    if (!is.na(class(object)[2]) && class(object)[2] == "pcrfit") {
      NEWMOD <- pcrfit(DATA, PRED.pos, RESP.pos, model = object$MODEL, opt.method = object$opt.method)        
    } else NEWMOD <- eval(as.call(NEWCALL))
    
    NEWPRED <- as.data.frame(DATA[i, PRED.pos])
    colnames(NEWPRED) <- PREDname
    y.hat <- as.numeric(predict(NEWMOD, newdata = NEWPRED))
    PRESS.res[i] <- DATA[i, RESP.pos] - y.hat
  }
  cat("\n")

  return(list(stat = sum(PRESS.res^2), residuals = PRESS.res))
}
  
  