peakArea <- function(x, y)
{
  ### define background slope fit
  x.first <- head(x, 1)
  x.last <- tail(x, 1)
  y.first <- head(y, 1)
  y.last <- tail(y, 1)
  x.pair <- c(x.first, x.last)
  y.pair <- c(y.first, y.last)
  LM <- lm(y.pair ~ x.pair)

  ### calculate background values for all x
  BASELINE <- predict(LM, newdata = data.frame(x.pair = x))
  
  ### baseline data
  y.base <- y - BASELINE

  ### calculate peak area
  SPLFN <- splinefun(x, y.base)
  AREA <- integrate(SPLFN, min(x, na.rm = TRUE), max(x, na.rm = TRUE))$value
  
  return(list(area = AREA, baseline = BASELINE))
}