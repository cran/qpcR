takeoff <- function(object, pval = 0.05, nsig = 3)
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
           
      which.top <- sapply(1:length(resl), function(x) all(resl[x:(x + nsig)]))
      min.top <- min(which(which.top == TRUE), na.rm = TRUE)
      top <- as.numeric(names(resl[min.top]))                    
      
      fluo.top <- as.numeric(predict(object, newdata = data.frame(Cycles = top)))
      return(list(top = top, f.top = fluo.top))
}