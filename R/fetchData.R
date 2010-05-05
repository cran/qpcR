fetchData <- function(object)
{
  ### 'pcrfit' object
  if (class(object) == "pcrfit") DATA <- object$DATA
  
  ### any other object
  if (class(object$call$data) == "name") DATA <- eval(object$call$data)
  else if (class(object$call$data) == "data.frame"  || class(object$call$data) == "matrix") DATA <- object$call$data
  else if (is.null(object$call$data)) DATA <- as.data.frame(sapply(all.vars(object$call$formula), function(a) get(a, envir = .GlobalEnv)))

  ### get variables from formula
  VARS <- all.vars(object$call$formula)
  LHS <- VARS[1]
  RHS <- VARS[-1]
  matchPRED <- which(!is.na(match(RHS, colnames(DATA))))
  PREDname <- RHS[matchPRED]
  PRED.pos <- which(colnames(DATA) == PREDname)
  RESP.pos <- which(colnames(DATA) == LHS)
  
  return(list(data = DATA, pred.pos = PRED.pos, resp.pos = RESP.pos, pred.name = PREDname))
}
  
  
  
  