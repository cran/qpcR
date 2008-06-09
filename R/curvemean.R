curvemean <- function(modlist, mean = c("amean", "gmean", "hmean"), fct = l5(), which = 1)
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

      nans <- which(is.nan(mean.x))

      DAT <- cbind(mean.x, pred.y)
      DAT <- DAT[complete.cases(DAT), ]

      newMod <- pcrfit(DAT, 1, 2, fct)
      pcrplot(MODLIST[[which]])
      
      addypoints <- pcrpred(newMod, nans, which = "y")

      for (i in slmod[-which]) {
             pcrplot(MODLIST[[i]], add = TRUE, col = i)
      }

      pcrplot(newMod, add = TRUE, lwd = 1.5)

      return(newMod)
}