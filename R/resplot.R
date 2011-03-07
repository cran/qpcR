resplot <- function(object, overlay = TRUE, ...)
{
  if (!is.numeric(object)) RESID <- as.numeric(residuals(object)) else RESID <- object
  ORD <- order(abs(RESID))
  COL <- vector()
  LEN <- 1:length(RESID)

  COL[ORD] <- LEN
  
  if (!overlay) {
    barplot(RESID, col = rev(heat.colors(length(RESID)))[COL],
             ylab = "residual value", ylim = c(1.2 * min(RESID), 1.2 * max(RESID)), ...)
  } else {
    par(mar = c(5.1, 4.1, 4.1, 4.1))
    plot(object)
    par(new = TRUE)
    BP <- barplot(RESID, space = 0.8, axes = FALSE, plot = FALSE)
    barplot(RESID, space = 0.8, axes = FALSE, xlim = c(min(BP), max(BP)),
    col = rev(heat.colors(length(RESID)))[COL], ylim = c(2 * min(RESID), 2 * max(RESID)))
    axis(side = 4)
    mtext("Residual value", side = 4, line = 2.5)
  }
}
