HQIC <- function (object) 
{
    LL <- logLik(object)
    k <- length(coef(object))
    n <- length(residuals(object))
    -2 * as.numeric(LL) + 2 * k * log(log(n))
}