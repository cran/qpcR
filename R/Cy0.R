Cy0 <- function(object, plot = FALSE, add = FALSE, ...)
{
      cpD1 <- efficiency(object, plot = FALSE, type = "cpD1")$cpD1
      Fluo <- pcrpred(object, cpD1, which = "y")
      slope <- object$fct$derivx(cpD1, t(coef(object)))
      Cy0 <- cpD1 - (Fluo/slope)
      if (plot) {
            pcrplot(object, ...) 
            add = TRUE
      }
      if (add) {
            points(cpD1, Fluo, col = "blueviolet", pch = 16, ...)
            abline((-cpD1 * slope) + Fluo, slope)
            abline(h = 0)
            points(Cy0, 0, col = "blueviolet", pch = 16, ...)
      }
      return(round(Cy0, 2))
}