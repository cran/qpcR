sliwin <- function (object, wsize = 5, border = 7, plot = TRUE) 
{
  cpD2 <- efficiency(object, plot = FALSE)$cpD2
  xmin <- round(cpD2) - border
  xmax <- round(cpD2) + border
  x <- object$DATA[xmin:xmax, 1]
  y <- log(abs(object$DATA[xmin:xmax, 2]))
  rtable <- vector()
  slopetable <- vector()    

  for (i in 1:(length(x) - wsize)) {
    xcount <- x[i:(i + wsize - 1)]
    ycount <- y[i:(i + wsize - 1)]
    DATA = data.frame(xcount, ycount)
    DATA <- DATA[!is.infinite(DATA[, 2]), ]     
    
    fit <- try(lm(ycount ~ xcount, data = DATA))
        
    if (inherits(fit, "try-error")) next
    fit.sum <- summary(fit)
    r.squared <- fit.sum$r.squared
    rtable[i] <- r.squared
    slope <- fit.sum$coefficients[2, 1]
    slopetable[i] <- slope
  }

  efftable <- exp(slopetable)
  rmax <- max(rtable)
  crmax <- which.max(rtable)
  effmax <- efftable[crmax]      
  fit.best <- lm(y[crmax:(crmax + wsize - 1)] ~ x[crmax:(crmax + wsize - 1)])
  init <- exp(as.numeric(coef(fit.best)[1]))

  if (plot) {
    par(mar = c(5.1, 4.5, 3.0, 2.1))
    par(mfrow = c(3, 1))
    plot(x[1:(length(x) - wsize)], rtable, ylim = c(0.9, 
            1), xlab = "Cycles", ylab = "R-squared", cex.axis = 1.3,
            cex.lab = 1.5)
    points(xmin + crmax - 1, rmax, col = 4, cex = 2)
    plot(x[1:(length(x) - wsize)], efftable, xlab = "Cycles", 
            ylab = "Efficiency", cex.axis = 1.3, cex.lab = 1.5)
    points(xmin + crmax - 1, effmax, col = 4, cex = 2)
    plot(x, y, xlab = "Cycles", ylab = "log(RFU)", cex.axis = 1.3, 
            cex.lab = 1.5)
    points(x[crmax:(crmax + wsize - 1)], y[crmax:(crmax + 
            wsize - 1)], col = 4, cex = 2)
    abline(fit.best, col = 2)
  }
  return(list(eff = effmax, rmax = rmax, init = init))
}
