expGrowth <- function (fixed = c(NA, NA, NA), names = c("span", "K", "plateau")) 
{
    numParm <- 3
    if (!is.character(names) | !(length(names) == numParm)) {
        stop("Not correct 'names' argument")
    }
    if (!(length(fixed) == numParm)) {
        stop("Not correct 'fixed' argument")
    }
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]
    fct <- function(x, parm) {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
        parmMat[, 1] * exp(parmMat[, 2] * x) + parmMat[, 3]
    }
    ssfct <- function(data) {
        x <- data[, 1]
        y <- data[, 2]
        if (is.na(fixed[3])) {
            plateau <- 0.95 * min(y)
        }
        else {
            plateau <- fixed[3]
        }
        span <- max(y) - plateau
        tempY <- log((y - plateau))
        coefVec <- coef(lm(tempY ~ x))
        span <- exp(coefVec[1])
        K <- coefVec[2]
        return(c(span, K, plateau)[notFixed])
    }
    pnames <- names[notFixed]
    deriv1 <- function(x, parm) {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
        helper1 <- exp(parmMat[, 2] * x)
        cbind(helper1, parmMat[, 1] * helper1 * x, 1)
    }
    deriv2 <- NULL
    returnList <- list(fct = fct, ssfct = ssfct, names = pnames, 
        deriv1 = deriv1, deriv2 = deriv2)
    class(returnList) <- "Exponential growth"
    invisible(returnList)
}
