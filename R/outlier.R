outlier <- function(object, pval = 0.05, nsig = 3)
{
      require(MASS, quietly = TRUE)
      
      CYC <- object$DATA[, 1]
      FLUO <- object$DATA[, 2]
      res <- vector()

      for (i in 5:length(CYC)) {
            mod <- lm(FLUO[1:i] ~ CYC[1:i], na.action = na.exclude)
            st <- studres(mod)              
            st1 <- tail(st, 1)
            pst1 <- 2 * (1 - pt(st1, df = mod$df.residual))              
            res <- c(res, pst1)
      }
      resl <- sapply(res, function(x) x < pval)       
      resl[is.na(resl)] <- FALSE
           
      which.outl <- sapply(1:length(resl), function(x) all(resl[x:(x + nsig)]))
      min.outl <- min(which(which.outl == TRUE), na.rm = TRUE)
      outl <- as.numeric(names(resl[min.outl]))                    
      
      fluo.outlier <- as.numeric(pcrpred(object, newdata = data.frame(Cycles = outl)))
      return(list(outl = outl, f.outl = fluo.outlier))
}