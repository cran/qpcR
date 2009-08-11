BIC <- function(object)
{
  npar <- length(coef(object))
  nobs <- length(residuals(object))
  as.numeric(-2 * logLik(object)) + (npar * log(nobs))
}
