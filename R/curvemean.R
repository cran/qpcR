curvemean <- function(ml, mean = c("amean", "gmean", "hmean"), model = l4, which = 1, plot = TRUE)
{
  mean <- match.arg(mean)
  REF.fluo <- fitted(ml[[which]])    
  DF <- NULL
  ml2 <- ml[-which]
  PRED.x <- as.data.frame(sapply(ml2, function(x) pcrpred(x, newdata = data.frame(Fluo = REF.fluo), which = "x")))
  DATA <- cbind(ml[[which]]$DATA[, 1], PRED.x)     

  gmean <- function(x) prod(x[!is.na(x)])^(1/length(x[!is.na(x)]))
  hmean <- function(x) length(x[!is.na(x)])/sum(1/x[!is.na(x)])

  MEAN <- switch(mean, amean = rowMeans(DATA, na.rm = TRUE),
                 gmean = apply(DATA, 1, function(x) gmean(x)),
                 hmean = apply(DATA, 1, function(x) hmean(x)))
                 
  SD <- apply(DATA, 1, function(x) sd(x, na.rm = TRUE))
   
  RES <- cbind(MEAN, REF.fluo)
  RES <- RES[complete.cases(RES), ]
  
  COL <- rainbow(length(ml))  

  newMod <- pcrfit(RES, 1, 2, model)  
  
  if (plot) {
    pcrplot(ml[[which]], col = COL[1], lty = 2)
    
    for (i in 1:length(ml2)) {
      pcrplot(ml2[[i]], add = TRUE, col = COL[i], lty = 2)
    }
    pcrplot(newMod, add = TRUE, lwd = 2)
  }  
       
  newMod$mean.x <- MEAN
  newMod$sd.x <- SD
  return(newMod)
}