pcrplot <- function (x, ..., level = NULL, col = FALSE, 
    conLevel, conName, grid = 100, legend, legendText, legendPos, 
    type = "average", lty, log = "", pch, xlab, ylab, xlim, 
    ylim, xt = NULL, xtlab = NULL, yt = NULL, 
    ytlab = NULL, add = FALSE, confband = NULL) 
{
    getRange <- function (x, y, xlim) {
                logVec <- (x >= xlim[1] & x <= xlim[2])
                return(range(y[logVec]))
    }
    colour <- col
    object <- x
    origData <- object$data
    doseDim <- ncol(origData) - 4
    if (doseDim > 1) {
        stop("No plot features for plots in more than two dimensions")
    }
    dose <- origData[, 1:doseDim]
    resp <- origData[, doseDim + 1]
    assayNo <- origData[, 3]
    assayNoOld <- origData[, 4]
    numAss <- length(unique(assayNo))
    doPlot <- is.null(level) || any(unique(assayNoOld) %in% level)
    if (!doPlot) {
        stop("Nothing to plot")
    }
    plotFct <- (object$curve)[[1]]
    logDose <- (object$curve)[[2]]
    if (missing(conLevel)) {
        conLevel <- ifelse(is.null(logDose), 0.01, log(0.01))
    }
    if (missing(conName)) {
        if (is.null(logDose)) {
            conName <- expression(0)
        }
        else {
            conName <- expression(-infinity)
        }
    }
    varNames <- colnames(origData)[1:(doseDim + 1)]
    if (missing(xlab)) {
        if (varNames[1] == "") {
            xlab <- "Dose"
        }
        else {
            xlab <- varNames[1]
        }
    }
    if (missing(ylab)) {
        if (varNames[2] == "") {
            ylab <- "Response"
        }
        else {
            ylab <- varNames[2]
        }
    }
    logX <- TRUE
    if ((log == "") || (log == "y")) {
        logX <- FALSE
    }
    if (missing(xlim)) {
        xLimits <- c(min(dose), max(dose))
    }
    else {
        xLimits <- xlim
    }
    conNameYes <- FALSE
    if ((xLimits[1] < conLevel) && logX) {
        xLimits[1] <- conLevel
        smallDoses <- dose < conLevel
        dose[smallDoses] <- conLevel
        conNameYes <- TRUE
    }
    if (xLimits[1] >= xLimits[2]) {
        stop("Argument 'conLevel' is set too high")
    }
    if (doseDim == 1) {
        if ((is.null(logDose)) && (logX)) {
            dosePts <- exp(seq(log(xLimits[1]), log(xLimits[2]), 
                length = grid))
        }
        else {
            dosePts <- seq(xLimits[1], xLimits[2], length = grid)
        }
    }
    else {
    }
    if (is.null(logDose)) {
        plotMat <- plotFct(dosePts)
    }
    else {
        plotMat <- plotFct(logDose^(dosePts))
    }
    numCol <- ncol(plotMat)
    maxR <- max(resp)
    options(warn = -1)
    maxPM <- apply(plotMat, 2, max, na.rm = TRUE)
    if (max(maxPM) > maxR) {
        maxPM <- maxPM[which.max(maxPM)]
    }
    else {
        maxPM <- maxR
    }
    options(warn = 0)
    if (missing(ylim)) {
        if (missing(xlim)) {
            yLimits <- c(min(resp), maxPM)
        }
        else {
            yLimits <- getRange(dose, resp, xLimits)
        }
    }
    else {
        yLimits <- ylim
    }
    naPlot <- FALSE
    for (i in 1:numCol) {
        logVec <- !(plotMat[, i] >= yLimits[1] & plotMat[, i] <= 
            yLimits[2])
        if (any(!is.na(logVec)) && any(logVec)) {
            plotMat[logVec, i] <- NA
            naPlot <- TRUE
        }
    }
    if (missing(legend)) {
        if (numCol == 1) {
            legend <- FALSE
        }
        else {
            legend <- TRUE
        }
    }
    colourVec <- rep(1, numCol)
    if (is.logical(colour) && colour) {
        colourVec <- 1:numCol
    }
    if (!is.logical(colour) && length(colour) == numCol) {
        colourVec <- colour
    }
    if (!is.logical(colour) && (!(length(colour) == numCol))) {
        colourVec <- rep(colour, numCol)
    }
    if (!missing(lty)) {
        if (length(lty) == 1) {
            ltyVec <- rep(lty, numCol)
        }
        else {
            ltyVec <- lty
        }
    }
    else {
        ltyVec <- 1:numCol
    }
    if (!missing(pch)) {
        if (length(pch) == 1) {
            pchVec <- rep(pch, numCol)
        }
        else {
            pchVec <- pch
        }
    }
    else {
        pchVec <- 1:numCol
    }
    ivMid <- rep(TRUE, grid)
    par(las = 1)
    if (!is.null(logDose)) {
        if (log == "x") {
            log <- ""
        }
        if ((log == "xy") || (log == "yx")) {
            log <- "y"
        }
    }
    eps1 <- 1e-08
    logVec <- !((dose < xLimits[1] - eps1) | (dose > xLimits[2] + 
        eps1))
    dose <- dose[logVec]
    resp <- resp[logVec]
    uniAss <- unique(assayNoOld)
    for (i in 1:numCol) {
        ###added error bar function by ANS###
        plotPoints <- switch(type, average = cbind(as.numeric(names(tapVec <- tapply(resp[assayNo == i], 
                                           dose[assayNo == i], mean))), tapVec), 
                      none = c(max(dosePts) +  10, max(c(maxPM, max(resp))) + 10), 
                             points = cbind(dose[assayNo == i], resp[assayNo == i]), 
                      all = cbind(dose[assayNo == i], resp[assayNo == i]), 
                            obs = cbind(dose[assayNo == i], resp[assayNo == i]),
                      errbar = cbind(as.numeric(names(tapVec <- tapply(resp[assayNo == i], 
                               dose[assayNo == i], mean))), tapVec))

        ###added error bar function by ANS###
        if (type == "errbar") {stdev <- cbind(as.numeric(names(tapVec <- tapply(resp[assayNo == i], 
                               dose[assayNo == i], sd))), tapVec)              
        }
        #####################################
        if (!add) {
            if (i == 1) {
                if (is.null(level) || uniAss[i] %in% level) {
                  plot(plotPoints, xlab = xlab, ylab = ylab, 
                    log = log, xlim = xLimits, ylim = yLimits, 
                    axes = FALSE, frame.plot = TRUE, col = colourVec[i], 
                    pch = pchVec[i], ...)
                  ###added error bar function by ANS###
                  if (type == "errbar") {arrows(plotPoints[, 1],plotPoints[, 2],
                                        plotPoints[, 1],plotPoints[, 2] + stdev[, 2],
                                        angle = 90, length = 0.025,col = colourVec[i])
                  }
                  #####################################
                  yaxisTicks <- axTicks(2)
                  yLabels <- TRUE
                  if (!is.null(yt)) {
                    yaxisTicks <- yt
                    yLabels <- yt
                  }
                  if (!is.null(ytlab)) {
                    yLabels <- ytlab
                  }
                  axis(2, at = yaxisTicks, labels = yLabels)
                  xaxisTicks <- axTicks(1)
                  xLabels <- as.expression(xaxisTicks)
                  if (conNameYes) {
                    xLabels[1] <- conName
                  }
                  lenXT <- length(xaxisTicks)
                  if (lenXT > 6) {
                    halfLXT <- floor(lenXT/2) - 1
                    chosenInd <- 1 + 2 * (0:halfLXT)
                    xaxisTicks <- xaxisTicks[chosenInd]
                    xLabels <- xLabels[chosenInd]
                  }
                  if (!is.null(xt)) {
                    if (as.character(xt[1]) == as.character(eval(conName))) {
                      xaxisTicks <- c(xaxisTicks[1], xt[-1])
                      xLabels <- c(conName, xt[-1])
                    }
                    else {
                      xaxisTicks <- xt
                      xLabels <- xt
                    }
                  }
                  if (!is.null(xtlab)) {
                    xLabels <- xtlab
                  }
                  axis(1, at = xaxisTicks, labels = xLabels)                  
                }
            }
            else {
                matchLevel <- match(unique(assayNoOld)[i], level)
                if ((!is.null(level)) && (!is.na(matchLevel)) && 
                  (matchLevel == 1)) {
                  plot(plotPoints, xlab = xlab, ylab = ylab, 
                    log = log, xlim = xLimits, ylim = yLimits, 
                    axes = FALSE, frame.plot = TRUE, col = colourVec[i], 
                    pch = pchVec[i], ...)
                  yaxisTicks <- axTicks(2)
                  yLabels <- TRUE
                  if (!is.null(yt)) {
                    yaxisTicks <- yt
                    yLabels <- yt
                  }
                  if (!is.null(ytlab)) {
                    yLabels <- ytlab
                  }
                  axis(2, at = yaxisTicks)
                  xaxisTicks <- axTicks(1)
                  xLabels <- as.expression(xaxisTicks)
                  if (conNameYes) {
                    xLabels[1] <- conName
                  }
                  if (!is.null(xt)) {
                    if (as.character(xt[1]) == as.character(eval(conName))) {
                      xaxisTicks <- c(xaxisTicks[1], xt[-1])
                      xLabels <- c(conName, xt[-1])
                    }
                    else {
                      xaxisTicks <- xt
                      xLabels <- xt
                    }
                  }
                  if (!is.null(xtlab)) {
                    xLabels <- xtlab
                  }
                  axis(1, at = xaxisTicks, labels = xLabels)
                }
                if (is.null(level) || ((!is.na(matchLevel)) && 
                  (matchLevel > 1))) {
                  points(plotPoints, col = colourVec[i], pch = pchVec[i])
                  ###added error bar function by ANS###
             	   if (type == "errbar") {arrows(plotPoints[, 1], plotPoints[, 2], plotPoints[, 1],
                                         plotPoints[, 2]+stdev[, 2], angle = 90,
                                         length = 0.025, col = colourVec[i])
                  }
                  ##################################### 
                }
            }
        }
        else {
            if (is.null(level) || uniAss[i] %in% level) {
                points(plotPoints, pch = i, col = colourVec[i], 
                  ...)
            }
        }
    }
    noPlot <- rep(FALSE, numCol)
    if (!(type == "obs")) {
        for (i in 1:numCol) {
            if (any(is.na(plotMat[, i])) & (!naPlot)) {
                noPlot[i] <- TRUE
                next
            }
            if (is.null(level) || unique(assayNoOld)[i] %in% 
                level) {
                lines(dosePts[ivMid], plotMat[ivMid, i], lty = ltyVec[i], 
                  col = colourVec[i], ...)
            }
        }
    }
    legendLevels <- as.character(unique(assayNoOld))
    if (!missing(legendText)) {
        lenLT <- length(legendText)
        if (lenLT == numAss) {
            legendLevels <- legendText
        }
        if (lenLT == 1) {
            legendLevels <- rep(legendText, numAss)
        }
    }
    levInd <- 1:numAss
    if (!is.null(level)) {
        legendLevels <- level
        levInd <- (1:numAss)[unique(assayNoOld) %in% level]
    }
    ltyVec[noPlot] <- 0
    if (type == "obs") {
        ltyVec[levInd] <- 0
    }
    if (type == "none") {
        pchVec[levInd] <- NA
    }
    if (legend) {
        if (!missing(legendPos)) {
            if ((is.numeric(legendPos)) && (length(legendPos) == 
                2)) 
                xlPos <- legendPos[1]
            ylPos <- legendPos[2]
        }
        else {
            xlPos <- xLimits[2]
            ylPos <- yLimits[2]
        }
        legend(xlPos, ylPos, legendLevels, lty = ltyVec[levInd], 
            pch = pchVec[levInd], col = colourVec[levInd], bty = "n", 
            xjust = 1, yjust = 1)
    }
    ###added confidence band by ANS###
    if (is.numeric(confband)) {
        cb <- confband(object, level = confband)
        lines(cb$x, cb$clo, col=2)
        lines(cb$x, cb$cup, col=2)
    }   
    ##################################
    par(las = 0)
    retData <- data.frame(dosePts, as.data.frame(plotMat))
    colnames(retData) <- c(colnames(origData)[1:doseDim], as.character(unique(assayNoOld)))
    invisible(retData)
}
