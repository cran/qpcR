ratiocalc <- function (
data, 
group = NULL, 
which.eff = c("sig", "sli", "exp"),
type.eff = c("individual", "mean.single", "median.single",
              "mean.pair", "median.pair"), 
which.cp = c("cpD2", "cpD1", "cpE", "cpR", "cpT", "Cy0"),
...)
{      
    if (class(data) != "pcrbatch")
        stop("data is not of class 'pcrbatch'!")
    if (is.null(group))
        stop("Please define 'group'!")
    if (length(group) != ncol(data) - 1)
        stop("Length of 'group' and 'data' do not match!")
    if (is.numeric(which.eff))
        which.eff <- which.eff
    else which.eff <- match.arg(which.eff)
    which.cp <- match.arg(which.cp)
    type.eff <- match.arg(type.eff)   
    DATA <- data[, -1]
    cpDat <- ts()
    effDat <- ts()
    cpNames <- NULL
    effNames <- NULL
    PATTERN <- unique(group)
    if (all(regexpr("rs", group, perl = TRUE) == -1))
        refNo <- TRUE
    else refNo <- FALSE
    for (i in 1:length(PATTERN)) {
        WHICH <- which(group == PATTERN[i])
        cpSEL <- which(data[, 1] == paste("sig.", which.cp, sep = ""))
        effSEL <- which(data[, 1] == paste(which.eff, ".eff", sep = ""))
        tempCP <- as.numeric(DATA[cpSEL, WHICH])
        tempEff <- as.numeric(DATA[effSEL, WHICH])
        if (is.numeric(which.eff))
            tempEff <- rep(which.eff, length(WHICH))
        cpDat <- cbind(cpDat, ts(tempCP))
        effDat <- cbind(effDat, ts(tempEff))
        cpNames <- c(cpNames, paste("cp.", PATTERN[i], sep = ""))
        effNames <- c(effNames, paste("eff.", PATTERN[i], sep = ""))
    }

    cpDat <- cpDat[, -1]
    effDat <- effDat[, -1]

    if (is.numeric(which.eff))
        type.eff <- "individual"
    if (type.eff == "mean.single")
        effDat <- t(replicate(nrow(effDat), apply(effDat, 2,
            function(x) mean(x, na.rm = TRUE))))
    if (type.eff == "median.single")
        effDat <- t(replicate(nrow(effDat), apply(effDat, 2,
            function(x) median(x, na.rm = TRUE))))
    if (type.eff == "mean.pair") {
        effDat[, 1:2] <- mean(effDat[, 1:2], na.rm = TRUE)
        if (!refNo)
            effDat[, 3:4] <- mean(effDat[, 3:4], na.rm = TRUE)
    }
    if (type.eff == "median.pair") {
        effDat[, 1:2] <- median(effDat[, 1:2], na.rm = TRUE)
        if (!refNo)
            effDat[, 3:4] <- median(effDat[, 3:4], na.rm = TRUE)
    }
    
    cpDat <- matrix(cpDat, ncol = length(cpNames))
    effDat <- matrix(effDat, ncol = length(effNames))
    
    allDat <- cbind(cpDat, effDat)
    colnames(allDat) <- c(cpNames, effNames)        
    
    if (refNo) {
      EXPR <- expression(eff.gc^cp.gc/eff.gs^cp.gs)
      TIES <- c(1, 2, 1, 2)
    }    
    else {
      EXPR <- expression((eff.gc^cp.gc/eff.gs^cp.gs)/(eff.rc^cp.rc/eff.rs^cp.rs))
      TIES <- c(1, 2, 1, 2, 1, 2, 1, 2)
    }
        
    CRIT <- "perm > init"
    PROP <- propagate(EXPR, allDat, do.perm = TRUE, ties = TIES, perm.crit = CRIT, ...)
    PROP <- c(PROP, list(data = allDat))     
    class(PROP) <- "ratiocalc"
    return(PROP)
}
