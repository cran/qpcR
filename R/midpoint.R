midpoint <- function(object, first.cyc = 10)
{
      fluo <- object$data[, 2]
      Fmax <- max(fluo, na.rm = TRUE)
      Fnoise <- sd(fluo[1:first.cyc], na.rm = TRUE)
      mp <- Fnoise * sqrt(Fmax/Fnoise)
      cyc.mp <- pcrpred(object, mp, which = "x")
      return(list(f.mp = mp, mp = cyc.mp))
}