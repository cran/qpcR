TmFind <- function(
TEMP = NULL,
FLUO = NULL,
span.smooth = NULL,
span.peaks = NULL,
is.deriv = FALSE,
Tm.opt = NULL)
{
  ### set dataframes to zero
  meltDATA <- NULL
  tempDATA  <- NULL
  derivDATA <- NULL
  
  ### cubic spline fitting and Friedman's Supersmoother
  ### on the first derivative curve
  SPLFN <- try(splinefun(TEMP, FLUO), silent = TRUE)
  if (inherits(SPLFN, "try-error")) return()
  seqTEMP <- seq(min(TEMP, na.rm = TRUE), max(TEMP, na.rm = TRUE), length.out = 10 * length(TEMP))
  meltDATA <- cbind(meltDATA, Fluo = SPLFN(seqTEMP))
  tempDATA <- cbind(tempDATA, Temp = seqTEMP)
  if (!is.deriv) derivVEC <- SPLFN(seqTEMP, deriv = 1) else derivVEC <- -SPLFN(seqTEMP, deriv = 0)
  SMOOTH <-  try(supsmu(seqTEMP, derivVEC, span = span.smooth), silent = TRUE)
  if (inherits(SMOOTH, "try-error")) return()
  derivDATA <- cbind(derivDATA, df.dT = -SMOOTH$y)
  
  ### find peaks in first derivative data
  PEAKS <- try(peaks(-SMOOTH$y, span = span.peaks)$x, silent = TRUE)
  if (inherits(PEAKS, "try-error")) return()
  TMs <- seqTEMP[PEAKS]
  TMs <- TMs[!is.na(TMs)]
  
  ### calculate difference to Tm.opt if given
  ### by residual sum-of-squares
  if (!is.null(Tm.opt)) {
    length(TMs) <- length(Tm.opt)
    RSS <- sum((Tm.opt - TMs)^2)
  } else RSS <- NA
  
  ### return data
  outDATA <- data.frame.na(tempDATA, meltDATA, derivDATA, Pars = c(span.smooth, span.peaks), RSS, Tm = TMs)
  return(outDATA)
}
