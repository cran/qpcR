calib2 <- function(refcurve, predcurve = NULL, thresh = "cpD2", term = NULL, dil = NULL, fct = l5(),
                  plot = TRUE, plot.map = TRUE, conf = 0.95, opt = c("none", "inter", "slope"),
                  opt.step = c(50, 50), quan = 0.5, slope = NULL, count = 1)
{
      require(gplots, quietly = TRUE)
      if (class(refcurve) != "modlist") stop("'refcurve' is not a 'modlist'!")
      if (!is.null(predcurve) & class(predcurve) != "modlist") stop("'predcurve' is not a 'modlist'!")
      if (thresh != "cpD2" && !is.numeric(thresh)) stop("'thresh' must be either 'cpD2' or numeric!")
      if (is.null(dil)) stop("Please define dilutions!")
      if (length(dil) != length(refcurve)) stop("Supply as many dilutions as number of PCRs in 'refcurve'!")
      opt <- match.arg(opt)

      lmod <- length(refcurve)
      lmod2 <- length(predcurve)        

      if (thresh == "cpD2") {
            CPD2 <- efficiency(refcurve[[1]], plot = FALSE)$cpD2
            start.y <- pcrpred(refcurve[[1]], CPD2, which = "y")
      } else start.y <- thresh

      if (opt != "none") {
            end.x <- qpcR:::outlier(refcurve[[lmod]])$outl
            end.y <- qpcR:::outlier(refcurve[[lmod]])$f.outl
      } else {
            end.y <- start.y
            opt.step[1] <- 1
      }
      
      if (is.numeric(term)) end.y <- term
      
      seqINTER <- seq(start.y, end.y, length.out = opt.step[1])

      if (opt == "slope") {
            delta.y <- start.y - end.y
            delta.x <- end.x
            if (opt.step[2] != 0) seqSLOPE <- seq(0, -(delta.y/delta.x), length.out = opt.step[2])
                  else seqSLOPE <- 0
      } else {
            seqSLOPE <- 0
      }

      seqALPHA <- atan(abs(seqSLOPE))*180/pi
      
      lsi <- length(seqINTER)
      lss <- length(seqSLOPE)
      
      resMat <- matrix(nrow = lsi, ncol = lss)      
      maxY <- max(sapply(refcurve, function(x) x$data[, 2]), na.rm = TRUE)
      maxX <- max(sapply(refcurve, function(x) x$data[, 1]), na.rm = TRUE)
      
      if (plot) {
            par(mfrow = c(2, 1))
            par(mar = c(5, 4, 2, 2))
      }
      
      dil.old <- dil
      dil <- log10(dil)

      ROOTFCT <- function(model, maxX, m, z) {
                        FCT <- function(x, parm, m, z) model$fct$fct(x, parm) - (m * x + z)
                        PARM <- t(coef(model))
                        pred.x <- try(uniroot(FCT, interval = c(1, maxX), parm = PARM, m = m, z = z), silent = TRUE)
                        if (inherits(pred.x, "try-error")) return()
                        pred.y <- predict(model, as.data.frame(pred.x$root))[1]
                        return(list(pred.x = pred.x$root, pred.y = pred.y))
      }
      
      LMFCT <- function(dil, ref, pred = NULL, conf) {
            linModY <- lm(ref ~ dil)
            conf.Y <- predict(linModY, interval = "confidence", level = conf)
            eff <- as.numeric(10^(-1/coef(linModY)[2]))
            FOM <- AIC(linModY)
            if (!is.null(pred)) {
                  linModX <- lm(dil ~ ref)
                  pred.conc <- sapply(as.numeric(pred), function(x) predict(linModX, data.frame(ref = x), interval = "confidence", level = conf))
            } else pred.conc <- NULL
            return(list(linModY = linModY, conf.Y = conf.Y, eff = eff, FOM = FOM, pred.conc = pred.conc[1, ], pred.conf = pred.conc[2:3, ]))
      }

      for (k in 1:lsi) {
            for (j in 1:lss) {
                  m <- seqSLOPE[j]
                  z <- seqINTER[k]
                  a <- seqALPHA[j]
                  
                  PRED <- lapply(refcurve, function(y) ROOTFCT(y, maxX, m, z))
                  PRED.X <- sapply(PRED, function(x) x$pred.x)
                  PRED.Y <- sapply(PRED, function(x) x$pred.y)

                  if (!is.null(predcurve) && count == 2) {
                        PRED2 <- lapply(predcurve, function(y) ROOTFCT(y, maxX, m, z))
                        PRED.X2 <- sapply(PRED2, function(x) x$pred.x)
                        PRED.Y2 <- sapply(PRED2, function(x) x$pred.y)
                  } else {
                        PRED.X2 <- NULL
                        PRED.Y2 <- NULL
                  }

            if (is.numeric(PRED.X)) LINMOD <- LMFCT(dil = dil, ref = PRED.X, pred = PRED.X2, conf = conf)
                  else break

            if (plot) {
                  MAIN1 <- paste("Intercept:", round(z, 2), "\nAlpha:", round(a, 2))
                  MAIN2 <- paste("AIC:", round(AIC(LINMOD$linModY), 2), "\nRsq:", round(Rsq(LINMOD$linModY), 5), "\nEff:", round(LINMOD$eff, 2))
                  pcrplot(refcurve[[1]], main = NULL, lwd = 1.3, xlab = "Cycles", ylab = "raw fluorescence", col = 1, ylim = c(0, maxY))
                  title(main = MAIN1, cex.main = 0.6)
                  sapply(2:lmod, function(x) pcrplot(refcurve[[x]], add = TRUE, lwd = 1.3, col = x))
                  abline(a = z, b = m, col = 1, lwd = 1.3)
                  abline(v = PRED.X, col = 1:lmod, lwd = 1.3)
                  points(PRED.X, PRED.Y, col = 1, pch = 16, cex = 0.5)
            
                  if (!is.null(predcurve) && count == 2) {
                        sapply(1:lmod2, function(x) pcrplot(predcurve[[x]], add = TRUE, lwd = 2, col = x))
                        points(PRED.X2, PRED.Y2, col = 1, pch = 15, cex = 0.5)
                        abline(v = PRED.X2, col = 1:lmod2, lwd = 2)
                  }

                  plot(dil, PRED.X, col = 1:lmod, pch = 16, cex = 1.3, xlab = "log(Dilution or copy number)", ylab = "threshold cycle")
                  title(main = MAIN2, cex.main = 0.6)
                  abline(LINMOD$linModY, lwd = 1.3)
                  lines(dil, LINMOD$conf.Y[, 2], col = 2)
                  lines(dil, LINMOD$conf.Y[, 3], col = 2)                  
                             
                  if (!is.null(predcurve) && count == 2) {
                        points(LINMOD$pred.conc, PRED.X2, pch = 15, col = 1:lmod2, cex = 1.5)
                        if (!all(is.na(LINMOD$pred.conc))) {                                                                                               
                              arrows(LINMOD$pred.conf[1, ], PRED.X2, LINMOD$pred.conf[2, ], PRED.X2, code = 3, angle = 90, length = 0.1)                              
                        }
                  }
                  
            }
            resMat[k, j] <- LINMOD$FOM
            if (!plot) cat("Intercept:", round(z, 2), "   Alpha:", round(a, 2), "   AIC:", round(AIC(LINMOD$linModY), 2), "   R-square:", Rsq(LINMOD$linModY), "\n", sep = "")
            if (PRED.Y[lmod] <= end.y) break
            }

      }

      if (count == 1) resMat.old <- resMat  
      
      dimnames(LINMOD$conf.Y) <- NULL
      if (count == 2) return(list(ref.cyc = PRED.X, pred.cyc = PRED.X2, pred.logconc = LINMOD$pred.conc, 
                                    pred.conc = 10^as.numeric(LINMOD$pred.conc), ref.conf = t(LINMOD$conf.Y[, 2:3]), pred.logconf = LINMOD$pred.conf, 
                                    pred.conf = 10^LINMOD$pred.conf, eff = LINMOD$eff, 
                                    aic = AIC(LINMOD$linModY), rsq = Rsq(LINMOD$linModY)))

      quan <- quantile(resMat, quan, na.rm = TRUE)
      resMat[resMat > quan] <- NA
      x.NA <- apply(resMat, 1, function(x) all(is.na(x)))
      y.NA <- apply(resMat, 2, function(x) all(is.na(x)))
      seqINTER <- seqINTER[!x.NA]
      seqALPHA <- seqALPHA[!y.NA]
      resMat <- as.matrix(resMat[!x.NA, !y.NA])
      BEST <- which(resMat == min(resMat, na.rm = TRUE), arr.ind = TRUE)
      bestINTER <- seqINTER[BEST[1]]
      bestSLOPE <- seqSLOPE[BEST[2]]
      bestFIT <- calib2(refcurve = refcurve, predcurve = predcurve, thresh = bestINTER, dil = dil.old, fct = fct,
                       plot = TRUE, conf = conf, opt = "none", slope = bestSLOPE, count = 2)
      
      bestFIT$aicMat <- resMat.old             

      if (opt != "none" && plot.map) {
            par(mfrow = c(1, 1))
            par(ask = TRUE)
            image(x = rev(seqINTER), y = seqALPHA, resMat, col = redblue(256), xlab = "Intercept", ylab = "Alpha", axes = FALSE)
            CEX.X <- 1 - (0.03 * sqrt(length(seqINTER)))
            CEX.Y <- 1 - (0.05 * sqrt(length(seqALPHA)))
            axis(1, at = seqINTER, labels = rev(round(seqINTER, 2)), cex.axis = CEX.X, las = 2)
            axis(2, at = seqALPHA, labels = round(seqALPHA, 2), las = 1, cex.axis = CEX.Y)
            INTspace <- 0.5 * (seqINTER[1] -  seqINTER[2])
            if (length(seqALPHA) > 1) ALPHAspace <- 0.5 * (seqALPHA[2] - seqALPHA[1]) else ALPHAspace <- 1
            rect(rev(seqINTER)[BEST[1]] - INTspace, seqALPHA[BEST[2]] - ALPHAspace,
                 rev(seqINTER)[BEST[1]] + INTspace, seqALPHA[BEST[2]] + ALPHAspace, lwd = 3)

            for (i in 1:length(seqINTER)) {
                  for (j in 1:length(seqALPHA)) {
                        text(rev(seqINTER)[i], seqALPHA[j], round(resMat[i, j], 2), cex = 0.3)
                  }
            }
      }
      return(bestFIT)
}           