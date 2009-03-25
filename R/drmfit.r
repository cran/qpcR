drmfit <- function(formula, curveid, weights, data = NULL, fct, na.action = na.omit, robust = "mean") 
{
    type <- "continuous"
    require(MASS, quietly = TRUE)
    options(na.action = deparse(substitute(na.action)))
    control <- drmc()
    useD <- control$useD
    constrained <- control$constr
    maxIt <- control$maxIt
    optMethod <- control$method
    relTol <- control$relTol
    warnVal <- control$warnVal
    rmNA <- control$rmNA
    errorMessage <- control$errorm
    noMessage <- control$noMessage
    options(warn = warnVal)
    selfStart <- TRUE
    bcVal <- NULL
    logDose <- NULL
    bcAdd <- NULL
    
    
    if ((!is.list(fct)) && (!is.function(fct))) {
        stop("No function or list given in argument 'fct'")
    }
    if (is.function(fct)) {
        fct <- drc:::fct2list(fct, 2)
    }
    if (is.null(names(fct))) {
        fct$fct <- fct[[1]]
        fct$ssfct <- fct[[2]]
        fct$names <- fct[[3]]
    }
    if (!is.function(fct$fct)) {
        stop("First entry in list to 'fct' NOT a function")
    }
    else {
        drcFct <- fct$fct
    }
    if (is.null(fct$ssfct)) {
        noSSfct <- TRUE
    }
    else {
        noSSfct <- FALSE
    }
    if ((!is.function(fct$ssfct)) && selfStart) {
        stop("Neither self starter function nor starting values provided")
    }
    else {
        ssfct <- fct$ssfct
    }
    if (is.null(fct$names) || (!is.character(fct$names))) {
        stop("Parameter names (as vector a strings) are NOT supplied")
    }
    else {
        parNames <- fct$names
        numNames <- length(parNames)
    }
    isDF <- is.function(fct$deriv1)
    if ((useD) && (isDF)) {
        dfct1 <- fct$deriv1
    }
    else {
        dfct1 <- NULL
    }
    if ((useD) && (is.function(fct$deriv2))) {
        dfct2 <- fct$deriv2
    }
    else {
        dfct2 <- NULL
    }
    anName <- deparse(substitute(curveid))
    if (length(anName) > 1) {
        anName <- anName[1]
    }
    mf <- match.call(expand.dots = FALSE)
    nmf <- names(mf)
    mnmf <- match(c("formula", "curveid", "data", "subset", "na.action", 
        "weights"), nmf, 0)
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf[c(1, mnmf)], parent.frame())
    mt <- attr(mf, "terms")
    dose <- model.matrix(mt, mf)[, -c(1)]
    resp <- model.response(mf, "numeric")
    origResp <- resp
    lenData <- length(resp)
    numObs <- length(resp)
    xDim <- ncol(as.matrix(dose))
    if (xDim > 1) {
        stop("drm() is only designed for 1-dim. dose vectors")
    }
    dimData <- xDim + 1
    varNames <- names(mf)
    varNames <- varNames[c(2:dimData, 1)]
    weights <- model.weights(mf)
    if (is.null(weights)) {
        weights <- rep(1, numObs)
    }
    missingIndices <- attr(mf, "na.action")
    if (is.null(missingIndices)) {
        removeMI <- function(x) {
            x
        }
    }
    else {
        removeMI <- function(x) {
            x[-missingIndices, ]
        }
    }
    assayNo <- model.extract(mf, "curveid")
    if (is.null(assayNo)) {
        assayNo <- rep(1, numObs)
    }
    uniqueNames <- unique(assayNo)
    colOrder <- order(uniqueNames)
    uniqueNames <- as.character(uniqueNames)
    assayNoOld <- assayNo
    colConvert <- function(vec) {
        len <- length(vec)
        assLev <- unique(vec)
        retVec <- rep(0, len)
        j <- 1
        for (i in 1:length(assLev)) {
            retVec[vec == assLev[i]] <- j
            j <- j + 1
        }
        return(retVec)
    }
    assayNo <- colConvert(assayNoOld)
    assayNames <- as.character(unique(assayNoOld))
    numAss <- length(assayNames)
    uniqueDose <- lapply(tapply(dose, assayNoOld, unique), length)
    udNames <- names(uniqueDose[uniqueDose == 1])
    if (length(udNames) > 0) {
        cm <- udNames
        if (!noMessage) {
            cat(paste("Control measurements detected for level: ", 
                udNames, "\n", sep = ""))
        }
        conInd <- assayNoOld %in% udNames
        assayNo[conInd] <- (assayNo[!conInd])[1]
        ciOrigIndex <- unique(assayNo)
        ciOrigLength <- numAss
        assayNames <- as.character(unique(assayNoOld[!conInd]))
        numAss <- length(assayNames)
        assayNo <- colConvert(assayNo)
        cm <- NULL
    }
    else {
        cm <- NULL
        ciOrigIndex <- unique(assayNo)
        ciOrigLength <- numAss
    }
    pmodelsList <- list()
    if (length(unique(assayNo)) == 1) {
            for (i in 1:numNames) {
                pmodelsList[[i]] <- matrix(1, numObs, 1)
            }
        }
        else {
            modelMat <- model.matrix(~factor(assayNo) - 1, level = unique(assayNo))
            for (i in 1:numNames) {
                pmodelsList[[i]] <- modelMat
            }
        }
        
    options(na.action = "na.omit")        
    pmodelsList2 <- list()
    for (i in 1:numNames) {
        colNames <- colnames(pmodelsList[[i]])
        if ((!is.null(cm)) && (!is.null(colNames))) {
            accm <- as.character(cm)
            pos <- grep(accm, colNames)
            if (length(pos) == 0) {
                candCol <- pmodelsList[[i]][, 1]
                if (!(length(assayNoOld[candCol == 1]) == 0) && 
                  (all(assayNoOld[candCol == 1] == accm))) {
                  pos <- 1
                }
            }
        }
        else {
            pos <- numeric(0)
        }
        if ((length(pos) > 0) && !(upperPos == i)) {
            pmodelsList2[[i]] <- as.matrix(pmodelsList[[i]][, 
                -pos])
        }
        else {
            pmodelsList2[[i]] <- as.matrix(pmodelsList[[i]])
        }
    }
    for (i in 1:numNames) {
        if (ncol(pmodelsList[[i]]) > numAss) {
            pmodelsList2[[i]] <- model.matrix(~factor(assayNo) - 
                1)
            colnames(pmodelsList2[[i]]) <- assayNames
        }
        else {
            pmodelsList2[[i]] <- as.matrix(pmodelsList[[i]])
        }
    }
    ncclVec <- rep(0, numNames)
    for (i in 1:numNames) {
        ncclVec[i] <- ncol(pmodelsList2[[i]])
    }
    parmPos <- c(0, cumsum(ncclVec)[-numNames])
    pnList <- drc:::drmParNames(numNames, parNames, pmodelsList2)
    parmVec <- pnList[[1]]
    parmVecA <- pnList[[2]]
    parmVecB <- pnList[[3]]
    parmIndex <- list()
    for (i in 1:numNames) {
        parmIndex[[i]] <- parmPos[i] + 1:ncclVec[i]
    }
    scaleFct <- fct$scaleFct
    if (!is.null(scaleFct)) {
        doseScaling <- 10^(floor(log10(median(dose))))
        if ((is.na(doseScaling)) || (doseScaling < 1e-10)) {
            doseScaling <- 1
        }
        respScaling <- 10^(floor(log10(median(resp))))
        if ((is.na(respScaling)) || (respScaling < 1e-10) || 
            (!identical(type, "continuous")) || (!is.null(bcVal))) {
            respScaling <- 1
        }
        longScaleVec <- rep(scaleFct(doseScaling, respScaling), 
            as.vector(unlist(lapply(parmIndex, length))))
    }
    else {
        doseScaling <- 1
        respScaling <- 1
        longScaleVec <- 1
    }
    startVecList <- list()
    if (!noSSfct) {
        startMat <- matrix(0, numAss, numNames)
        lenASS <- length(formals(ssfct))
        if (lenASS > 1) {
            doseresp <- data.frame(x = dose, y = origResp)
            ssFct <- function(dframe) {
                ssfct(dframe, doseScaling, respScaling)
            }
        }
        else {
            doseresp <- data.frame(x = dose/doseScaling, y = origResp/respScaling)
            ssFct <- ssfct
        }
        isfi <- is.finite(dose)
        for (i in 1:numAss) {
            indexT1 <- (assayNo == i)
            if (any(indexT1)) {
                logVec <- indexT1 & isfi
                startMat[i, ] <- ssFct(doseresp[logVec, ])
            }
            else {
                startMat[i, ] <- rep(NA, numNames)
            }
            if (sum(!is.na(startMat[i, ])) == 1) {
                upperPos <- (1:numNames)[!is.na(startMat[i, ])]
            }
        }                        
        nrsm <- nrow(startMat)
        for (i in 1:numNames) {
            sv <- rep(0, max(nrsm, ncclVec[i]))
            indVec <- 1:ncclVec[i]
            sv[1:nrsm] <- startMat[, i]
            sv <- sv[!is.na(sv)]
            isZero <- (sv == 0)
            sv[isZero] <- mean(sv)
            startVecList[[i]] <- sv[indVec]
        }
        startVec <- unlist(startVecList)
    }
    else {
        startVec <- start
    }
    if (!selfStart && !noSSfct) {
        lenReq <- length(startVec)
        if (length(start) == lenReq) {
            startVec <- start/longScaleVec
        }
        else {
            stop(paste("Wrong number of initial parameter values. ", 
                lenReq, " values should be supplied", sep = ""))
        }
    }
    if (selfStart) {
        startVec <- drc:::drmConvertParm(startVec, startMat, assayNo, 
            pmodelsList2)
    }
    startVecSc <- startVec
    parmMatrix <- matrix(0, numObs, numNames)
    parm2mat <- function(parm) {
        for (i in 1:numNames) {
            parmMatrix[, i] <- pmodelsList2[[i]] %*% parm[parmIndex[[i]]]
        }
        return(parmMatrix)
    }
    drcFct1 <- function(dose, parm) {
        drcFct(dose, parm2mat(parm))
    }
    if (!is.null(fct$retFct)) {
        drcFct <- fct$retFct(doseScaling, respScaling)
        drcFct1 <- function(dose, parm) {
            drcFct(dose, parm2mat(parm))
        }
    }
    if (is.null(cm)) {
        multCurves <- function(dose, parm) {
            drcFct1(dose, parm)
        }
    }
    else {
        iv <- assayNoOld == cm
        niv <- !iv
        fctEval <- rep(0, numObs)
        multCurves <- function(dose, parm) {
            parmVal <- parm2mat(parm)
            fctEval[iv] <- parmVal[iv, upperPos, drop = FALSE]
            fctEval[niv] <- drcFct(dose[niv], parmVal[niv, , 
                drop = FALSE])
            fctEval
        }
    }
    robustFct <- drc:::drmRobust(robust, match.call(), numObs, length(startVec))
    if (!is.null(dfct1)) {
        dmatfct <- function(dose, parm) {
            dfct1(dose, parm2mat(parm))
        }
    }
    else {
        dmatfct <- NULL
    }
    if (!is.null(bcVal)) {
        bcfct <- function(x, lambda, bctol, add = bcAdd) {
            if (abs(lambda) > bctol) {
                return(((x + add)^lambda - 1)/lambda)
            }
            else {
                return(log(x + add))
            }
        }
        bcTol <- 0.02
        resp <- bcfct(resp, bcVal, bcTol)
        multCurves2 <- function(dose, parm) {
            bcfct(multCurves(dose, parm), bcVal, bcTol)
        }
    }
    else {
        multCurves2 <- multCurves
    }
    if (type == "continuous") {
        estMethod <- drc:::drmEMls(dose, resp, multCurves2, startVecSc, 
            robustFct, weights, rmNA, dmf = dmatfct, doseScaling = doseScaling, 
            respScaling = respScaling)
    }
    if (identical(type, "binomial")) {
        estMethod <- drc:::drmEMbinomial(dose, resp, multCurves2, startVecSc, 
            robustFct, weights, rmNA, doseScaling = doseScaling)
    }
    if (identical(type, "Poisson")) {
        estMethod <- drc:::drmEMPoisson(dose, resp, multCurves2, startVecSc, 
            doseScaling = doseScaling)
    }
    opfct <- estMethod$opfct
    lowerLimits <- rep(-Inf, length(startVec))
    upperLimits <- rep(Inf, length(startVec))
    opdfctTemp <- estMethod$opdfct1
    appFct <- function(x, y) {
        tapply(x, y, sum)
    }
    if (!is.null(opdfctTemp)) {
        opdfct1 <- function(parm) {
            as.vector(apply(opdfctTemp(parm), 2, appFct, assayNo))
        }
    }
    else {
        opdfct1 <- NULL
    }
    startVecSc <- as.vector(startVecSc)
    nlsFit <- drc:::drmOpt(opfct, opdfct1, startVecSc, optMethod, constrained, 
        warnVal, upperLimits, lowerLimits, errorMessage, maxIt, 
        relTol, parmVec = parmVec, trace = control$trace)
    if (!nlsFit$convergence) {
        return(nlsFit)
    }
    if (!is.null(scaleFct)) {
        nlsFit$value <- nlsFit$value * (respScaling^2)
        nlsFit$par <- nlsFit$par * longScaleVec
        nlsFit$hessian <- nlsFit$hessian * (1/outer(longScaleVec/respScaling, 
            longScaleVec/respScaling))
    }
    if (!is.null(fct$retFct)) {
        drcFct <- fct$retFct(1, 1)
        drcFct1 <- function(dose, parm) {
            drcFct(dose, parm2mat(parm))
        }
    }
    nlsSS <- nlsFit$value
    nlsDF <- numObs - length(startVec)
    if (!is.null(cm)) {
        iVec <- (1:numAss)[!(uniqueNames == cm)]
    }
    else {
        iVec <- 1:numAss
    }
    pickCurve <- rep(0, length(iVec))
    for (i in iVec) {
        pickCurve[i] <- (1:numObs)[assayNo == i][1]
    }
    parmMat <- matrix(NA, numAss, numNames)
    fixedParm <- (estMethod$parmfct)(nlsFit)
    parmMat[iVec, ] <- parm2mat(fixedParm)[pickCurve, ]
    if (!is.null(cm)) {
        parmMat[-iVec, upperPos] <- parm2mat(fixedParm)[assayNoOld == 
            cm, , drop = FALSE][1, upperPos]
    }
    rownames(parmMat) <- assayNames
    pmFct <- function(fixedParm) {
        if (!is.null(cm)) {
            iVec <- (1:numAss)[!(uniqueNames == cm)]
        }
        else {
            iVec <- 1:numAss
        }
        if (!is.null(cm)) {
            parmMat[-iVec, upperPos] <- parm2mat(fixedParm)[assayNoOld == 
                cm, , drop = FALSE][1, upperPos]
        }
        rownames(parmMat) <- assayNames
        return(parmMat)
    }
    parmMat <- pmFct((estMethod$parmfct)(nlsFit))
    pfFct <- function(parmMat) {
        plotFct <- function(dose) {
            if (xDim == 1) {
                lenPts <- length(dose)
            }
            else {
                lenPts <- nrow(dose)
            }
            curvePts <- matrix(NA, lenPts, ciOrigLength)
            for (i in 1:numAss) {
                if (i %in% iVec) {
                  parmChosen <- parmMat[i, complete.cases(parmMat[i, 
                    ])]
                  parmMat2 <- matrix(parmChosen, lenPts, numNames, 
                    byrow = TRUE)
                  curvePts[, ciOrigIndex[i]] <- drcFct(dose, 
                    parmMat2)
                }
                else {
                  curvePts[, i] <- rep(NA, lenPts)
                }
            }
            return(curvePts)
        }
        return(plotFct)
    }
    plotFct <- pfFct(parmMat)
    predVec <- multCurves2(dose, fixedParm)
    resVec <- resp - predVec
    resVec[is.nan(predVec)] <- 0
    diagMat <- matrix(c(predVec, resVec), numObs, 2)
    colnames(diagMat) <- c("Predicted values", "Residuals")
    if (robust %in% c("median", "trimmed", "tukey", "winsor")) {
        nlsFit$value <- (mad(resVec, 0)^2) * nlsDF
    }
    if (robust %in% c("lms", "lts")) {
        scaleEst <- 1.4826 * (1 + 5/(numObs - length(nlsFit$par))) * 
            sqrt(median(resVec^2))
        w <- (resVec/scaleEst < 2.5)
        nlsFit$value <- sum(w * resVec^2)/(sum(w) - length(nlsFit$par))
    }
    robust <- switch(robust, median = "median", trimmed = "metric trimming", 
        tukey = "Tukey's biweight", winsor = "metric Winsorizing", 
        lms = "least median of squares", lts = "least trimmed squares")
    sumVec <- c(bcVal, NA, NA, NA, nlsSS, nlsDF, numObs)
    sumList <- list(lenData = numObs, alternative = NULL, df.residual = numObs - 
        length(startVec))
    callDetail <- match.call()
    if (is.null(callDetail$fct)) {
        callDetail$fct <- substitute(l4())
    }
    
    dataSet <- data.frame(dose, origResp, assayNo, assayNoOld, 
        weights)
    names(dataSet) <- c(varNames, anName, anName, "weights")
    hfct1 <- function(x) {
        uniVec <- unique(x[!is.na(x)])
        rv <- rep(NA, length(x))
        for (i in 1:length(uniVec)) {
            rv[abs(x - uniVec[i]) < 1e-12] <- i
        }
        rv
    }
    hfct2 <- function(x) {
        length(unique(x))
    }
    mat1 <- t(apply(t(parmMat), 1, hfct1))
    cnccl <- head(cumsum(ncclVec), -1)
    if (nrow(mat1) == 1) {
        mat1 <- t(mat1)
    }
    mat1[-1, ] <- mat1[-1, ] + cnccl
    if (isDF) {
        deriv1Mat <- fct$deriv1(dose, parmMat[assayNo, , drop = FALSE])
    }
    else {
        deriv1Mat <- NULL
    }
    returnList <- list(NULL, nlsFit, list(plotFct, logDose), 
        sumVec, startVecSc * longScaleVec, list(parmVec, parmVecA, 
            parmVecB), diagMat, callDetail, dataSet, t(parmMat), 
        fct, robust, estMethod, numObs - length(startVec), sumList, 
        function(x) {
            x
        }, pmFct, pfFct, type, mat1, logDose, cm, deriv1Mat, 
        anName, data, weights, list(dose = dose, resp = resp, 
            curveid = assayNoOld, origResp = origResp))
    names(returnList) <- c("varParm", "fit", "curve", "summary", 
        "start", "parNames", "predres", "call", "data", "parmMat", 
        "fct", "robust", "estMethod", "df.residual", "sumList", 
        "scaleFct", "pmFct", "pfFct", "type", "indexMat",  
        "cm", "deriv1", "curveVarNam", "origData", "weights", 
        "dataList")
    class(returnList) <- c("drc", class(fct))
    return(returnList)
}
