curvemean <- function(
ml, 
type = c("fluo", "expval"), 
mean = c("amean", "gmean", "hmean"),
which = 1, 
plot = TRUE)
{
  mean <- match.arg(mean)
  type <- match.arg(type)
  REF.fluo <- fitted(ml[[which]])   
  DATA <- ml[[which]]$DATA[, 1] 
  ml2 <- ml[-which]       
      
  for (i in 1:length(ml2)) {
    PRED.x <- predict(ml2[[i]], newdata = data.frame(Fluo = REF.fluo), which = "x")[, 1]     
    DATA <- cbind(DATA, round(PRED.x, 2))  
  }       
   
  COMPL <- complete.cases(DATA)
  DATA <- DATA[COMPL, ]    
  
  gmean <- function(x) prod(x, na.rm = TRUE)^(1/length(x[!is.na(x)]))
  hmean <- function(x) length(x[!is.na(x)])/sum(1/x, na.rm = TRUE)    
  
  MEAN <- vector(length = nrow(DATA))
  CYCS <- DATA[, 1]
  EFFS <- sapply(ml, function(x) eff(x)$eff.y[CYCS * 100])  
  
  for (i in 1:nrow(DATA)) { 
    X <- DATA[i, ]
    if (type == "fluo") {
      MEAN[i] <- switch(mean, amean = mean(X, na.rm = TRUE),
                              gmean = gmean(X),
                              hmean = hmean(X))
    } else {
      EFF <- mean(EFFS[i, ], na.rm = TRUE)
      max.X <- max(X, na.rm = TRUE)      
      temp1 <- abs(max.X - X)
      temp2 <- switch(mean, amean = mean(EFF^temp1, na.rm = TRUE),
                            gmean = gmean(EFF^temp1),
                            hmean = hmean(EFF^temp1))
      temp3 <- logb(temp2, EFF)
      MEAN[i] <- max.X - temp3 
     }                                
  }    
 
  MEAN[!is.finite(MEAN)] <- NA
  RES <- cbind(MEAN, REF.fluo[COMPL])    
  RES <- RES[complete.cases(RES), ]    
  COL <- rainbow(length(ml))                     
  newMod <- pcrfit(RES, 1, 2, ml[[1]]$MODEL) 
    
  if (plot) {
    plot(ml[[which]], col = COL[1], lty = 2)     
    
    for (i in 1:length(ml2)) {
      plot(ml2[[i]], add = TRUE, col = COL[i], lty = 2)     
    }       
    
    plot(newMod, add = TRUE, lwd = 2)
  }  
       
  return(newMod)
}