curvemean <- function(modlist, mean = c("amean", "gmean", "hmean"), fct = l5(), which = 1, plot = TRUE)
{
      mean <- match.arg(mean)
      lmod <- length(modlist)
      MODLIST <- modlist
      slmod <- 1:lmod
      nrdata <- nrow(MODLIST[[which]]$data)
      DF <- NULL
      SEQ <- seq(MODLIST[[which]]$data[1, 1], MODLIST[[which]]$data[nrdata, 1], by = 1)
      
      pred.y <- pcrpred(MODLIST[[which]], SEQ)
      
      for (i in slmod[-which]) {
             pred.x <- pcrpred(MODLIST[[i]], newdata = pred.y, which = "x")
             DF <- cbind(DF, pred.x)
      }

      DF <- cbind(SEQ, DF)

      gmean <- function(x) prod(x[!is.na(x)])^(1/length(x[!is.na(x)]))
      hmean <- function(x) length(x[!is.na(x)])/sum(1/x[!is.na(x)])

      mean.x <- switch(mean, amean = rowMeans(DF),
                               gmean = apply(DF, 1, function(x) gmean(x)),
                               hmean = apply(DF, 1, function(x) hmean(x)))
      sd.x <- apply(DF, 1, function(x) sd(x, na.rm = TRUE))
      nans <- which(is.nan(mean.x))

      DAT <- cbind(mean.x, pred.y)
      DAT <- DAT[complete.cases(DAT), ]
                           
      COL <- rainbow(lmod)

      newMod <- pcrfit(DAT, 1, 2, fct)
      if (plot) pcrplot(MODLIST[[which]], col = COL[1], lty = 2)
      
      if (plot) {
        for (i in slmod[-which]) {
          pcrplot(MODLIST[[i]], add = TRUE, col = COL[i], lty = 2)
        }
        pcrplot(newMod, add = TRUE, lwd = 2)          
      }
      
      newMod$mean.x <- mean.x
      newMod$sd.x <- sd.x
      return(newMod)
}