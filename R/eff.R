eff <- function(object, sequence = NULL, plot = FALSE) 
{
    if (is.null(sequence)) {
      MIN <- min(object$DATA[, 1], na.rm = TRUE)
      MAX <- max(object$DATA[, 1], na.rm = TRUE)
      DIVS <- 0.01        
    } else {
      MIN <- sequence[1]
      MAX <- sequence[2]
      DIVS <- sequence[3]      
    }    
    SEQ <- seq(MIN, MAX, by = DIVS)
    
    coefVec <- coef(object)       
    FCT <- object$MODEL$fct      
    EFF <- function(x) {
      F1 <- FCT(x, coefVec)     
      F2 <- FCT(x - 1, coefVec)     
      F1/F2            
    }
    EFFres <- EFF(SEQ)        
    maxCYC <-  SEQ[which.max(EFFres)]
            
    if (plot) {
      plot(SEQ, EFFres, xlab = "Cycles", ylab = "Efficiency")
      abline(v = maxCYC, col = 2, lwd = 2)
      mtext(maxCYC, side = 1, at = maxCYC, col = 2)
    }
    return(list(eff.x = SEQ, eff.y = EFFres, effmax.x = maxCYC, effmax.y = max(EFFres, na.rm = TRUE)))
}
