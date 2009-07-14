pcropt2 <- function(object, nsample = 5, B = 100, plot = TRUE, alpha = 0.05, ...)
{
  DATA <- object$DATA
  resMat <- NULL
  ndata <- nrow(DATA)

  for (i in 1:B) {
    SELECT <- sample(1:ndata, nsample)
    newData <- DATA[-SELECT, ]
    newModel <- pcrfit(newData, 1, 2, model = object$MODEL, ...)
    res <- efficiency(newModel, plot = FALSE, ...)
    resMat <- rbind(resMat, unlist(res))
    if (i %% 10 == 0) cat(i) else cat(".")
    if (i %% 50 == 0) cat("\n")
    flush.console()
}
  MEANS <- apply(resMat, 2, function(x) round(mean(x, na.rm = TRUE), 8))
  SDS <- apply(resMat, 2, function(x) round(sd(x, na.rm = TRUE), 8))
  CONFINT <- function(x, alpha) quantile(x, c(alpha/2, 1 - (alpha/2)), na.rm = TRUE)
  CONF <- apply(resMat, 2, function(x) round(CONFINT(x, alpha = alpha), 8))
  
  if (plot) {
    par(mfrow = c(4, 4))
    par(mar = c(1, 2, 2, 1))

    for (i in 1:ncol(resMat)) {
      if (!all(is.na(resMat[, i]))) boxplot(resMat[, i], main = colnames(resMat)[i], cex = 0.2, ...)
      abline(h = CONF[1, i], col = 2, lwd = 2)
      abline(h = CONF[2, i], col = 2, lwd = 2)
    }

  }
  return(list(raw = resMat, mean = MEANS, sd = SDS, upper = CONF[2, ], lower = CONF[1, ]))
}