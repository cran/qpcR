efficiency <- function (object, plot = TRUE, type = "cpD2", shift = 0, amount = NULL) 
{
    if (type != "cpD2" & type != "cpD1" & type != "maxE" & type != 
        "expR" & !is.numeric(type)) {
        stop("invalid estimation type")
    }
    if (!is.null(amount) && !is.numeric(amount)) 
        stop("amount must be numeric!")
    par(mar = c(5, 4, 4, 4))
    rv <- summary(object)$resVar
    aicc <- AICc(object)
    cyc.col <- which(colnames(object$data) == all.vars(object$call)[2])
    cyc <- object$data[, cyc.col]
    coefVec <- coef(object)
    obj.fct <- object$fct$fct
    effcurve <- eff(object)
    sequence <- effcurve$eff.x
    pMat <- matrix(coefVec, length(sequence), length(coefVec), 
        byrow = TRUE)
    res.grad <- object$fct$derivx(sequence, pMat)
    obj.fun <- qpcR:::typeid(object)
    if (obj.fun == "l5" || obj.fun == "l4" || obj.fun == "l3") {
        res.hess <- deriv2.l(sequence, pMat, obj.fun)
    }
    if (obj.fun == "b5" || obj.fun == "b4" || obj.fun == "b3") {
        res.hess <- deriv2.b(sequence, pMat, obj.fun)
    }
    E1res <- effcurve$eff.y
    maxD2 <- which.max(res.hess)
    cycmaxD2 <- sequence[maxD2]
    maxD1 <- which.max(res.grad)
    cycmaxD1 <- sequence[maxD1]
    cycmaxE1 <- NA
    expR <- NA
    shiftCyc <- NA
    exp.cyc <- NA

    if (type == "cpD2") {
        maxE1 <- E1res[maxD2 + (100 * shift)]
        CYC <- cycmaxD2 + shift 
	   if (shift != 0) shiftCyc <- cycmaxD2 + shift       
    }
    if (type == "cpD1") {
        maxE1 <- E1res[maxD1 + (100 * shift)]
  	   CYC <- cycmaxD1 + shift
        if (shift != 0) shiftCyc <- cycmaxD1 + shift        
    }
    if (type == "maxE") {
        maxE <- which.max(E1res)
        maxE1 <- E1res[maxE + (100 * shift)]
        cycmaxE1 <- sequence[maxE]
        CYC <- cycmaxE1 + shift
        if (shift != 0) shiftCyc <- cycmaxE1 + shift             
    }
    if (is.numeric(type)) {
        type.cyc <- which(sequence == type)
        maxE1 <- E1res[type.cyc]
        if (shift != 0) shiftCyc <- type + shift
        CYC <- type + shift
    }
    if (type == "expR") {
        expR <- maxD2 - (maxD1 - maxD2)
        exp.cyc <- sequence[expR]
        maxE1 <- E1res[expR + (100 * shift)]
        if (shift != 0) shiftCyc <- exp.cyc + shift
        CYC <- exp.cyc + shift
    }
    fluo <- pcrpred(object, "y", newdata = CYC)
    init <- fluo/(maxE1^CYC)
    if (is.numeric(amount)) 
        CF <- amount/init
    else CF <- NA
    if (plot) {
        pcrplot(object, lwd = 1.5, type = "all", main = NA, cex.main = 0.9)
        aT <- axTicks(side = 4)
        cF <- max(aT/1)
        axis(side = 4, at = aT, labels = round(aT/cF + 1, 2), 
            col = 4, col.axis = 4, las = 1, hadj = 0.3)
        lines(sequence, (E1res * cF) - cF, col = 4, lwd = 1.5)
        points(CYC, (maxE1 * cF) - cF, col = 4, pch = 16)
        points(CYC, fluo, col = 1, pch = 16)
        mtext(side = 4, "Efficiency", line = 2, col = 4)
        lines(sequence, res.grad, col = 2, lwd = 1.5)
        lines(sequence, res.hess, col = 3, lwd = 1.5)
        abline(h = (maxE1 * cF) - cF, lwd = 1.5, col = 4)
        abline(v = cycmaxD1, lwd = 1.5, col = 2)
        abline(h = fluo, col = 1)
        axis(side = 2, at = fluo, labels = round(fluo, 3), col = 1, 
            col.axis = 1, cex.axis = 0.7, las = 1)
        if (type == "maxE") 
            abline(v = cycmaxE1, lwd = 1.5, col = 4)
        if (type == "expR") 
            abline(v = exp.cyc, lwd = 1.5, col = 6)
        if (is.numeric(type)) 
            abline(v = type, lwd = 1.5, col = 6)
        abline(v = cycmaxD2, lwd = 1.5, col = 3)
        if (!is.null(shiftCyc)) 
            abline(v = shiftCyc, lwd = 1.5, col = 7)
        mtext(paste("cpD2:", cycmaxD2), line = 0, col = 3, adj = 0.65, 
            cex = 0.9)
        mtext(paste("cpD1:", cycmaxD1), line = 0, col = 2, adj = 0.35, 
            cex = 0.9)
        if (type == "maxE") 
            mtext(paste("cpE:", cycmaxE1), line = 1, col = 4, 
                adj = 0.65, cex = 0.9)
        if (type == "expR") 
            mtext(paste("cpR:", exp.cyc), line = 1, col = 6, adj = 0.65, 
                cex = 0.9)
        if (is.numeric(type)) 
            mtext(paste("ct:", type), line = 1, col = 6, adj = 0.65, 
                cex = 0.9)
        mtext(paste("Eff:", round(maxE1, digits = 3)), line = 1, 
            col = 4, adj = 0.35, cex = 0.9)
        mtext(paste("resVar:", round(rv, digits = 5)), line = 2, 
            col = 1, adj = 0.35, cex = 0.9)
        mtext(paste("AICc:", round(aicc, digits = 2)), line = 2, 
            col = 1, adj = 0.65, cex = 0.9)
        mtext(paste("Model:", obj.fun), line = 3, col = 1, adj = 0.5, 
            cex = 0.9)
    }
    return(list(eff = maxE1, resVar = round(summary(object)$resVar, 
        digits = 8), AICc = AICc(object), AIC = AIC(object), 
        Rsq = Rsq(object), cpD1 = cycmaxD1, cpD2 = cycmaxD2, 
        cpE = cycmaxE1, cpR = exp.cyc, fluo = fluo, init = init, 
        cf = CF))
}
