REST <- function(eff.goi, eff.ref, data, B = 10000, alpha = 0.05, replace = FALSE, bootstrap = NULL, perm = c("pairwise", "single"))
{
  perm <- match.arg(perm)     
  SAMPLES <- matrix(nrow = B, ncol  = nrow(data))
  RATIOS1 <- NULL 
  RATIOS2 <- NULL    
  COUNTER <- 0    
  pVec <- NULL
  NROW <- 1:nrow(data) 
  repMat <- NULL    
  if (is.null(bootstrap)) sampleN <- nrow(data) else sampleN <- round(bootstrap * nrow(data)) 
                                      
  while (COUNTER < B) {
    RN1 <- sample(NROW, sampleN, replace = replace)     
    RN2 <- sample(NROW, sampleN, replace = replace) 
    if (perm == "single") {
      RN3 <- sample(NROW, sampleN, replace = replace)
      RN4 <- sample(NROW, sampleN, replace = replace) 
    }  else {
      RN3 <- NA
      RN4 <- NA
    }
      
    if (perm == "pairwise") {
      ref.c1 <- data[RN1, 1] 
      goi.c1 <- data[RN1, 2]
      ref.s1 <- data[RN2, 3]
      goi.s1 <- data[RN2, 4]
    } else {
      ref.c1 <- data[RN1, 1] 
      goi.c1 <- data[RN2, 2]
      ref.s1 <- data[RN3, 3]
      goi.s1 <- data[RN4, 4]
    }  
    
    nom1 <- eff.goi^(goi.c1 - goi.s1)
    denom1 <- eff.ref^(ref.c1 - ref.s1)
    RATIO1 <- nom1/denom1     
    RATIOS1 <- c(RATIOS1, RATIO1)     
    
    PERM <- cbind(RN1, RN2, RN3, RN4)
    repMat <- rbind(repMat, PERM)
         
    COUNTER <- COUNTER + length(RATIO1)                 
  } 
  
  uniques <- unique(repMat) 
  CONF <- quantile(RATIOS1, c(alpha/2, 1 - (alpha/2)), na.rm = TRUE)            
  invisible(list(ratios = sort(RATIOS1), conf = CONF, uniques = uniques)) 
}

