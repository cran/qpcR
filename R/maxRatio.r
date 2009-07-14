maxRatio <- function(ml, plot = TRUE, ...)
{
  COLS <- rainbow(length(ml))
  RES <- lapply(ml, function(x) efficiency(x, type = "maxE", plot = FALSE))
  RATIO <- round(sapply(RES, function(x) x$eff), 2)
  FCN <- sapply(RES, function(x) x$cpE)
  FCNA <- round(FCN - logb(RATIO, 2), 2)
  EFFS <- lapply(ml, function(x) eff(x))
  NAMES <- sapply(ml, function(x) x$names)

  if (plot) {
    par(mfrow = c(3, 1)) #split.screen(c(3, 1))    
    par(mar = c(5, 5, 1, 1))
    pcrplot(ml[[1]], col = COLS[1], xlab = "Cycles", ylab = "RFU", ...)
    sapply(2:length(ml), function(x) pcrplot(ml[[x]], add = TRUE, col = COLS[x], ...))
    #screen(2)
    par(mar = c(5, 5, 1, 1))
    plot(EFFS[[1]]$eff.x, EFFS[[1]]$eff.y - 1, type = "l", col = COLS[1],
         xlab = "Cycles", ylab = "Ratio", ylim = c(0, 2), ...)
    sapply(2:length(EFFS), function(x) lines(EFFS[[x]]$eff.x, EFFS[[x]]$eff.y - 1, col = COLS[x], ...))
    abline(v = FCN, col = COLS)
    #screen(3)
    par(mar = c(5, 5, 1, 1))
    plot(FCN, RATIO, pch = 16, col = COLS, ylab = "MR", ...)
  }
  
  #close.screen(all = TRUE)
  return(list(mr = RATIO, fcn = FCN, fcna = FCNA, names = NAMES))
}

