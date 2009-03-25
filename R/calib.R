calib <- function (data, cyc.col = 1, f.col = 2:ncol(data), groups = NULL, 
    thresh = 1, method = "var", dil = NULL, pred = NULL, fct = l5(), 
    plot = TRUE, conf.lev = 0.95, opt.min = NULL, opt.plot = FALSE) 
{
    if (is.null(dil)) 
        stop("Please define dilutions!")
    if (is.null(groups)) 
        DATA <- data[, f.col]
    else {
        DATA <- apply(data[, f.col], 1, function(x) tapply(x, 
            groups, function(x) mean(x, na.rm = TRUE)))
        DATA <- t(DATA)
    }
    rl <- list()
    for (i in 1:ncol(DATA)) {
        y <- DATA[, i]
        x <- data[, cyc.col]
        frame <- as.data.frame(cbind(x, y))
        m <- drm(y ~ x, data = frame, fct = fct)
        rl[[i]] <- m
    }
    cp <- efficiency(rl[[1]], plot = FALSE)$cpD2
    if (is.numeric(opt.min)) 
        cp <- efficiency(rl[[1]], plot = FALSE)$cpD1
    if (method == "cpD2" || is.numeric(opt.min)) 
        thresh <- pcrpred(rl[[1]], newdata = cp)
    ct <- NULL
    aic.opt <- NULL
    thresh.opt <- NULL
    if (plot) {
        par(mfrow = c(2, 1))
        par(mar = c(5, 4, 2, 2))
    }
    if (is.numeric(opt.min)) 
        endloop <- opt.min
    else endloop <- thresh
    for (j in seq(thresh, endloop, by = -signif(thresh/1000, 
        1))) {
        for (i in 1:length(rl)) {
            if (plot) {
                if (i == 1) {
                  pcrplot(rl[[i]], xlab = "Cycles", ylab = "raw fluorescence", 
                    ylim = c(0, 3 * j), type = "none", lwd = 2)
                  abline(h = 0, lwd = 2)
                  abline(h = j, lty = 2, lwd = 2)
                }
                else pcrplot(rl[[i]], col = i, add = TRUE, type = "none", 
                  lwd = 2)
            }
            m <- rl[[i]]
            predx <- pcrpred(m, j, "x")
            ct <- c(ct, predx)
            if (plot) 
                segments(predx, 0, predx, j, col = i, lwd = 2)
        }
        rsq <- NULL
        adj.rsq <- NULL
        aic <- NULL
        predres <- NULL
        if (!is.null(dil)) {
            dilut <- log10(dil)
            if (plot) 
                plot(dilut, ct, col = 1:length(rl), cex = 1.5, 
                  pch = 16, xlab = "log10(dilution)", ylab = "threshold cycle")
            FIT <- lm(ct ~ dilut)
            if (!is.null(conf.lev)) {
                conf.FIT <- predict(FIT, interval = "confidence", level = conf.lev) 
		lines(dilut, conf.FIT[, 2], col = 2)
                lines(dilut, conf.FIT[, 3], col = 2)
            }
            FIT.s <- summary(FIT)
            rsq <- FIT.s$r.squared
            adj.rsq <- FIT.s$adj.r.squared
            aic <- AIC(FIT)
            aic.opt <- c(aic.opt, aic)
            thresh.opt <- c(thresh.opt, j)
            if (plot) 
                abline(FIT)
            title(paste("AIC:", signif(aic, 4)), cex.main = 0.8)
        }
        if (!is.null(pred)) {
            predfun <- function(y) {
                (y - coef(FIT)[1])/coef(FIT)[2]
            }
            predres <- vector()
            for (i in 1:length(pred)) {
                val <- predfun(pred[i])
                val.ul <- 10^val
                predres <- as.numeric(c(predres, val.ul))
                if (plot) 
                  points(val, pred[i], pch = 1, cex = 1.5)
            }
        }
        Eff <- as.numeric(10^(-1/coef(FIT)[2]))
        if (is.numeric(opt.min)) {
            ct <- NULL
            cat("Thresh:", j, " aic:", aic, "Rsq:", rsq, "Eff:", 
                Eff, "\n")
        }
    }
    min.opt <- which.min(aic.opt)
    if (is.numeric(opt.min) && opt.plot) {
        par(mfrow = c(1, 1))
        plot(thresh.opt, aic.opt, xlab = "raw fluorescence", 
            ylab = "AIC of linear fit")
        abline(v = thresh.opt[min.opt], col = 2, lwd = 2)
    }
    par(mfrow = c(1, 1))
    par(mar = c(5.1, 4.1, 4.1, 2.1))
    return(list(ct = ct, r.squared = rsq, adj.r.squared = adj.rsq, 
        AIC = aic, Eff = Eff, pred = predres, thresh.opt = thresh.opt[min.opt], 
        AIC.opt = aic.opt[min.opt]))
}
