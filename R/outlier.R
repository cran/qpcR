outlier <- function(object, pval = 0.05, nsig = 3)
{
      cyc <- object$data[, 1]
      fluo <- object$data[, 2]
      res <- vector()

      for (i in 5:length(cyc)) {
            mod <- lm(fluo[1:i] ~ cyc[1:i], na.action = na.exclude)
            st <- qpcR:::studtest(mod)
            st1 <- tail(st, 1)
            pst1 <- 1 - pt(st1, df = mod$df.residual)
            res <- c(res, pst1)
      }
      resl <- sapply(res, function(x) x < pval)
      resl[is.na(resl)] <- FALSE
      for (i in 1:length(resl - nsig - 1)) {
            window <- (0:(nsig - 1)) + i
            if (sum(resl[window] == TRUE) == nsig) {
                  outl <- as.numeric(names(resl[i]))
                  break
            }
      }
      fluo.outlier <- pcrpred(object, outl, which = "y")
      return(list(outl = outl, f.outl = fluo.outlier))
}