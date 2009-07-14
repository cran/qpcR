calib <- function(
refcurve, 
predcurve = NULL, 
thresh = "cpD2", 
term = NULL, 
dil = NULL, 
plot = TRUE,
plot.map = TRUE, 
conf = 0.95, 
opt = c("none", "inter", "slope"),
opt.step = c(50, 50), 
quan = 0.5, 
slope = NULL, 
count = 1)
{
  require(gplots, quietly = TRUE)

  if (class(refcurve) != "modlist") stop("'refcurve' is not a 'modlist'!")
  if (!is.null(predcurve) & class(predcurve) != "modlist") stop("'predcurve' is not a 'modlist'!")
  if (thresh != "cpD2" && !is.numeric(thresh)) stop("'thresh' must be either 'cpD2' or numeric!")
  if (is.null(dil)) stop("Please define dilutions!")
  if (length(dil) != length(refcurve)) stop("Supply as many dilutions as number of PCRs in 'refcurve'!")
  opt <- match.arg(opt)

  lref <- length(refcurve)
  lpred <- length(predcurve)        

  if (thresh == "cpD2") {
    cpD2 <- efficiency(refcurve[[1]], plot = FALSE)$cpD2
    start.y <- as.numeric(pcrpred(refcurve[[1]], newdata = data.frame(Cycles = cpD2)))
  } else start.y <- thresh   

  if (opt != "none") {
    end.x <- outlier(refcurve[[lref]])$outl
    end.y <- outlier(refcurve[[lref]])$f.outl
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
  maxY <- max(sapply(refcurve, function(x) x$DATA[, 2]), na.rm = TRUE)
  maxX <- max(sapply(refcurve, function(x) x$DATA[, 1]), na.rm = TRUE)  
      
  if (plot) {
    par(mfrow = c(2, 1))
    par(mar = c(5, 4, 2, 2))
  }              
      
  dil.old <- dil
  dil <- log10(dil)

  ROOTFCT <- function(fit, maxX, m, z) {
    FCT <- function(x, parm, m, z) fit$MODEL$fct(x, parm) - (m * x + z)
    PARM <- coef(fit)       
    pred.x <- try(uniroot(FCT, interval = c(1, maxX), parm = PARM, m = m, z = z), silent = TRUE)      
    if (inherits(pred.x, "try-error")) return()
    pred.y <- as.numeric(pcrpred(fit, newdata = data.frame(Cycles = pred.x$root))[1])       
    return(list(pred.x = pred.x$root, pred.y = pred.y))
  }
      
  LMFCT <- function(dil, ref, pred = NULL, conf) {    
    linModY <- lm(ref ~ dil)
    conf.Y <- predict(linModY, interval = "confidence", level = conf)[, 2:3]        
    eff <- as.numeric(10^(-1/coef(linModY)[2]))
    FOM1 <- AIC(linModY)
    FOM2 <- AICc(linModY)
    FOM3 <- Rsq(linModY)
    FOM4 <- Rsq.ad(linModY)    
    
    if (!is.null(pred)) {
      linModX <- lm(dil ~ ref)
      pred.conc <- sapply(as.numeric(pred), function(x) predict(linModX, data.frame(ref = x), interval = "confidence", level = conf))      
    } else pred.conc <- NULL
        
    return(list(linModY = linModY, conf.Y = conf.Y, eff = eff, FOM1 = FOM1, FOM2 = FOM2,
                FOM3 = FOM3, FOM4 = FOM4, pred.conc = pred.conc[1, ], pred.conf = pred.conc[2:3, ]))
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
        MAIN2 <- paste("AIC:", round(LINMOD$FOM1, 2), "\nRsq:", round(LINMOD$FOM3, 5), "\nEff:", round(LINMOD$eff, 2))
        pcrplot(refcurve, main = NULL, lwd = 1.3, colvec = 1:lref)
        title(main = MAIN1, cex.main = 0.6)         
        abline(a = z, b = m, col = 1, lwd = 1.3)
        abline(v = PRED.X, col = 1:lref, lwd = 1.3)
        points(PRED.X, PRED.Y, col = 1, pch = 16, cex = 0.5)
            
        if (!is.null(predcurve) && count == 2) {
          pcrplot(predcurve, add = TRUE, lwd = 2, colvec = 1:lpred)
          points(PRED.X2, PRED.Y2, col = 1, pch = 15, cex = 0.5)
          abline(v = PRED.X2, col = 1:lpred, lwd = 2)
        }

        plot(dil, PRED.X, col = 1:lref, pch = 16, cex = 1.3, xlab = "log(Dilution or copy number)", ylab = "threshold cycle")
        title(main = MAIN2, cex.main = 0.6)
        abline(LINMOD$linModY, lwd = 1.3)
        lines(dil, LINMOD$conf.Y[, 1], col = 2)
        lines(dil, LINMOD$conf.Y[, 2], col = 2)                  
                             
        if (!is.null(predcurve) && count == 2) {
          points(LINMOD$pred.conc, PRED.X2, pch = 15, col = 1:lpred, cex = 1.5)
          if (!all(is.na(LINMOD$pred.conc))) {                                                                                               
            arrows(LINMOD$pred.conf[1, ], PRED.X2, LINMOD$pred.conf[2, ], PRED.X2, code = 3, angle = 90, length = 0.1)                              
          }
        }
      }
            
      resMat[k, j] <- LINMOD$FOM1
      if (!plot) cat("Intercept:", round(z, 2), "   Alpha:", round(a, 2), "   AIC:", round(LINMOD$FOM1, 2), "   R-square:", LINMOD$FOM, "\n", sep = "")
      if (PRED.Y[lref] <= end.y) break
    }
  }
  
  if (count == 1) resMat.old <- resMat  
  dimnames(LINMOD$conf.Y) <- NULL
      
  if (count == 2) return(list(refcyc = PRED.X, 
                              refcyc.conf = t(LINMOD$conf.Y), 
                              predcyc = PRED.X2,                                                        
                              predconc = 10^as.numeric(LINMOD$pred.conc),                              
                              predconc.conf = 10^LINMOD$pred.conf,                              
                              eff = LINMOD$eff, 
                              aic = LINMOD$FOM1, aicc = LINMOD$FOM2, rsq = LINMOD$FOM3, rsq.ad = LINMOD$FOM4))   

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
  bestFIT <- calib(refcurve = refcurve, predcurve = predcurve, thresh = bestINTER, dil = dil.old,
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