resplot <- function(object, ...)
{
      if (!is.numeric(object)) resid <- as.numeric(residuals(object)) else resid <- object
      ORD <- order(abs(resid))
      COL <- vector()

      for (i in 1:length(resid)) {
            COL[ORD[i]] <- i
      }

      barplot(resid, col = rev(heat.colors(length(resid)))[COL],
                  ylab = "residual value", ylim = c(1.2 * min(resid), 1.2 * max(resid)), ...)
}