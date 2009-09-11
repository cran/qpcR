fitprob <- function(object)
{
    rss <- RSS(object)
    n <- length(residuals(object))
    p <- length(coef(object))
    pchisq(rss, n-p)
}







