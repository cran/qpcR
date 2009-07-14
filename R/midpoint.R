midpoint <- function(object, noise.cyc = 1:5)
{
      Fluo <- object$DATA[, 2]
      Fmax <- max(Fluo, na.rm = TRUE)
      Fnoise <- sd(Fluo[noise.cyc], na.rm = TRUE)
      mp <- Fnoise * sqrt(Fmax/Fnoise)
      cyc.mp <- as.numeric(pcrpred(object, newdata = data.frame(Fluo = mp), which = "x"))
      return(list(f.mp = mp, mp = cyc.mp))
}