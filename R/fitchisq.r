fitchisq <- function(object, error = NULL)
{
  if (any(class(object) == "replist")) object <- object[[1]] else object <- object
    
  fetchDATA <- fetchData(object)
  DATA <- fetchDATA$data
  PRED.pos <- fetchDATA$pred.pos
  RESP.pos <- fetchDATA$resp.pos
  PRED.name <- fetchDATA$pred.name     

  if (length(DATA[, PRED.pos]) > length(unique(DATA[, PRED.pos]))) {
    error <- tapply(DATA[, RESP.pos], DATA[, PRED.pos], function(x) sd(x, na.rm = TRUE))
  } 
  else if (is.null(error)) return(list(chi2 = NA, chi2.red = NA, p.value = NA))
  else {
    error <- rep(error, length.out = length(DATA[, PRED.pos]))
    names(error) <- DATA[, PRED.pos]
  }      
    
  res <- residuals(object)
  n <- length(res)
  p <- length(coef(object))
  df <- n - p
  m <- match(DATA[, PRED.pos], names(error))
  CHISQ <- sum(res^2/error[m]^2)
  CHISQ.red <- CHISQ/df
  p.value <- 1 - pchisq(CHISQ, df)
  return(list(chi2 = CHISQ, chi2.red = CHISQ.red, p.value = p.value))
}