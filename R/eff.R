eff <- function (object, sequence = NULL, plot = FALSE) 
{
    if (is.null(sequence)) {
        sequence <- seq(min(object$data[, 1], na.rm = TRUE), 
            max(object$data[, 1], na.rm = TRUE), by = 0.01)
    }
    coefVec <- coef(object)
    obj <- object
    obj.fun <- summary(object)$fctName
    obj.fct <- object$fct$fct
    pMat <- matrix(coefVec, length(sequence), length(coefVec), 
        byrow = TRUE)
    E1 <- function(x) {
        FE1 <- obj.fct(x, pMat)
        FE2 <- obj.fct(x - 1, pMat)
        FE3 <- FE1/FE2
        return(FE3)
    }
    E1res <- E1(sequence)
    if (plot) plot(sequence, E1res, xlab = "Cycles", ylab = "Efficiency")
    return(list(eff.x = sequence, eff.y = E1res, effmax = max(E1res, na.rm = TRUE)))
}
