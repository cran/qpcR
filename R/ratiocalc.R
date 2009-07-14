ratiocalc <- function(
data, 
group = NULL, 
which.eff = c("sig", "sli", "exp"),
type.eff = c("individual", "mean.single", "median.single", "mean.pair", "median.pair"),
which.cp = c("cpD2", "cpD1", "cpE", "cpR", "cpT", "Cy0"),
perm = c("cp", "eff", "both", NULL),
pval = c("up", "down"), 
... 
)
{
      require(gtools, quietly = TRUE) 
          
      if (class(data) != "pcrbatch") stop("data is not of class 'pcrbatch'!")
      if (is.null(group)) stop("Please define 'group'!")  
      if (length(group) != ncol(data) - 1) stop("Length of 'group' and 'data' do not match!")     
      
      ### match selection from 'pcrbatch'
      if (is.numeric(which.eff)) which.eff <- which.eff else which.eff <- match.arg(which.eff)
      which.cp <- match.arg(which.cp)
      type.eff <- match.arg(type.eff)
      perm <- match.arg(perm)   
      pval <- match.arg(pval)    
      DATA <- data[, -1]
      cpDat <- NULL
      effDat <- NULL
      cpNames <- NULL
      effNames <- NULL                     
            
      ### pattern matching for 'group'
      PATTERN <- unique(group)        
      if (all(regexpr("rs", group, perl = TRUE) == -1)) refNo <- TRUE else refNo <- FALSE     ## is reference tag to be found?
          
      for (i in 1:length(PATTERN)) {          ## for each pattern item
        WHICH <- which(group == PATTERN[i])           ## get match in 'group'               
        cpSEL <- which(data[, 1] == paste("sig.", which.cp, sep = ""))
        effSEL <- which(data[, 1] == paste(which.eff, ".eff", sep = ""))
        tempCP <- as.numeric(DATA[cpSEL, WHICH])             
        tempEff <- as.numeric(DATA[effSEL, WHICH])
        if (is.numeric(which.eff)) tempEff <- rep(which.eff, length(WHICH))  
        cpDat <- cbind(cpDat, tempCP) 
        effDat <- cbind(effDat, tempEff)
        cpNames <- c(cpNames, paste("cp.", PATTERN[i], sep = ""))
        effNames <- c(effNames, paste("eff.", PATTERN[i], sep = ""))                      
      } 
           
      colnames(cpDat) <- cpNames
      colnames(effDat) <- effNames       
                      
      if (is.numeric(which.eff)) type.eff <- "individual"
      if (type.eff == "mean.single") effDat <- t(replicate(nrow(effDat), apply(effDat, 2, function(x) mean(x, na.rm = TRUE))))
      if (type.eff == "median.single") effDat <- t(replicate(nrow(effDat), apply(effDat, 2, function(x) median(x, na.rm = TRUE))))
      if (type.eff == "mean.pair") {
        effDat[, 1:2] <- mean(effDat[, 1:2], na.rm = TRUE)
        if (!refNo) effDat[, 3:4] <- mean(effDat[, 3:4], na.rm = TRUE)
      }
      if (type.eff == "median.pair") {
        effDat[, 1:2] <- median(effDat[, 1:2], na.rm = TRUE)
        if (!refNo) effDat[, 3:4] <- median(effDat[, 3:4], na.rm = TRUE)
      }
      
      allDat <- cbind(cpDat, effDat)      
      
      if (refNo) EXPR <- expression(eff.gc^cp.gc/eff.gs^cp.gs)        ## in case of no reference genes
      else EXPR <- expression((eff.gc^cp.gc/eff.gs^cp.gs)/(eff.rc^cp.rc/eff.rs^cp.rs))   ## in case of reference genes       
       
      ## define ties for the 'NULL' permutation dataset
      ## option 'cp' will swap crossing points between groups as in REST software
      if (!refNo) TIES <- switch(perm, cp = c(1, 1, 2, 2, NA, NA, NA, NA),
                                       eff = c(NA, NA, NA, NA, 1, 1, 2, 2),
                                       both = c(1, 1, 2, 2, 3, 3, 4, 4),
                                       NULL = NULL) 
      else TIES <- switch(perm, cp = c(1, 1, NA, NA),
                                eff = c(NA, NA, 1, 1),
                                both = c(1, 1, 2, 2),
                                NULL = NULL)
                                
      PERMCRIT <- switch(pval, up = "m1 >= m2", down = "m1 <= m2")
                                                               
      PROP <- propagate(EXPR, allDat, do.perm = TRUE, ties = TIES, perm.crit = PERMCRIT, ...)
      PROP <- c(PROP, list(data = allDat)) 
      class(PROP) <- "ratiocalc"            
      return(PROP)     
} 