ERBCP <- function(object, corfact = 1)
{
      mod <- efficiency(object, plot = FALSE)
      expReg <- mod$cpD2 - (corfact*(mod$cpD1 - mod$cpD2))
      f.expReg <- pcrpred(object, expReg, which = "y")
      return(list(expReg = expReg, f.expReg = f.expReg))
}