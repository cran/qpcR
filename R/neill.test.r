neill.test <- function(object, grouping)
{
  fetchDATA <- fetchData(object)
  DATA <- fetchDATA$data
  PRED.pos <- fetchDATA$pred.pos
  RESP.pos <- fetchDATA$resp.pos
  PRED.name <- fetchDATA$pred.name
   
  X <- DATA[, PRED.pos]    
    
  if (length(X) > length(unique(X))) stop("Found replicate response values! Please use some other test (i.e. Lack-of-fit).")
  
  noCluster <- floor(length(X)/2)

  if (missing(grouping)) {
    for (i in noCluster:(length(coef(object)) + 1)) {
      grouping <- cutree(hclust(dist(X)), k = i)
      if (all(tapply(X, grouping, length) > 1)) break
    }
  }

  M <- length(unique(grouping))    
  N <- nrow(DATA)
  denDF <- N - M

  if (denDF <= 0) stop("Too many groups in 'grouping'")
    
  p <- N - df.residual(object)
  numDF <- M - p
    
  if (numDF <= 0) stop("Too few groups in 'grouping'")
    
  resVec <- residuals(object)
  resAver0 <- tapply(resVec, grouping, mean)
  resAver <- rep(resAver0, tapply(grouping, grouping, length))
  resDiff <- resVec - resAver
  F <- (denDF/numDF) * (sum(resAver * resAver)/(sum(resDiff * resDiff)))
  p <- pf(F, numDF, denDF, lower.tail = FALSE)
    
  return(p)
}


