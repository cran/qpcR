pcropt1 <- function (object, fact = 3, opt = FALSE, ...) 
{
    resMat <- NULL
    window.l <- NULL
    window.u <- NULL
    aicc <- NULL
    resvar <- NULL
    eff <- NULL
    init.exp <- NULL
    init.sig <- NULL
    compl <- try(efficiency(object, plot = FALSE, ...))
    if (inherits(compl, "try-error")) 
        stop("Could not initialize optimization. Try different 'fact'!")
    cpD1 <- round(compl$cpD1)
    cpD2 <- round(compl$cpD2)
    lower <- cpD1 - fact * (cpD1 - cpD2)
    upper <- cpD1 + fact * (cpD1 - cpD2)
    lowerseq <- 1:(lower - 1)
    upperseq <- nrow(object$data):(upper + 1)
    for (i in lowerseq) {
        for (j in upperseq) {
            newData <- object$data[i:j, ]
            assign("newData", newData, envir = .GlobalEnv)
            newCurve <- update(object, data = newData)
 		  F0 <- pcrpred(newCurve, which = "y", newdata = 0)
            if (opt) 
                newCurve <- mchoice(newCurve)
            compl <- try(efficiency(newCurve, ...))
            window.l <- c(window.l, i - lower)
            window.u <- c(window.u, j - upper)
            aicc <- c(aicc, compl$AICc)
            resvar <- c(resvar, compl$resVar)
            eff <- c(eff, compl$eff)
            init.exp <- c(init.exp, compl$init)
            init.sig <- c(init.sig, F0)
        }
    }
    resMat <- cbind(window.l, window.u, aicc, resvar, eff, init.exp, init.sig)
    return(as.data.frame(resMat))
}
