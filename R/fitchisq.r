fitchisq <- function(object, error = NULL)
{
  if (any(class(object) == "replist")) object <- object[[1]] else object <- object
  if (any(class(object) == "pcrfit")) DATA <- object$DATA else DATA <- eval(object$call$data)
  if (is.null(DATA)) stop("Please create model from a dataframe of values!")

  VARS <- all.vars(object$call$formula)
  LHS <- VARS[1]
  RHS <- VARS[-1]
  matchPRED <- which(!is.na(match(RHS, colnames(DATA))))
  PREDname <- RHS[matchPRED]
  
  X <- which(colnames(DATA) == PREDname)
  Y <- which(colnames(DATA) == LHS)

  if (length(DATA[, X]) > length(unique(DATA[, X]))) {
    error <- tapply(DATA[, Y], DATA[, X], function(x) sd(x, na.rm = TRUE))
  } else if (is.null(error)) return(list(chi2 = NA, chi2.red = NA, p.value = NA))
  else {
    error <- rep(error, length.out = length(DATA[, X]))
    names(error) <- DATA[, X]
  }
  
  res <- residuals(object)
  n <- length(res)
  p <- length(coef(object))
  df <- n - p
  m <- match(DATA[, X], names(error))
  CHISQ <- sum(res^2/error[m]^2)
  CHISQ.red <- CHISQ/df
  p.value <- 1 - pchisq(CHISQ, df)
  return(list(chi2 = CHISQ, chi2.red = CHISQ.red, p.value = p.value))
}