Rsq <- function(object)
{
    response <- object$data[, 2]
    curve <- object$data[, 4]
    uniCurve <- unique(curve)
    lenUC <- length(uniCurve)
    numerator <- tapply(residuals(object)^2, curve, sum)
    denominator <- tapply((response - mean(response))^2, curve, 
        sum)
    totnum <- sum(residuals(object)^2)
    totden <- sum((response - mean(response))^2)
    if (lenUC == 1) {
        rsq <- matrix(c(1 - numerator/denominator), 1, 1)
    }
    else {
        rsq <- matrix(c(1 - numerator/denominator, 1 - totnum/totden), 
            lenUC + 1, 1)
    }
    return(as.vector(rsq))	
} 
