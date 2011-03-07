ratiocalc <- function(
data, 
group = NULL, 
which.eff = c("sig", "sli", "exp", "mak"),
type.eff = c("individual", "mean.single", "median.single",
              "mean.pair", "median.pair"), 
which.cp = c("cpD2", "cpD1", "cpE", "cpR", "cpT", "Cy0"),
...)
{      
    if (class(data)[2] != "pcrbatch")
        stop("data is not of class 'pcrbatch'!")

    if (is.null(group))
        stop("Please define 'group'!")    

    if (length(group) != ncol(data) - 1)
        stop("Length of 'group' and 'data' do not match!")

    if (!is.numeric(which.eff)) which.eff <- match.arg(which.eff)

    if (!is.numeric(which.cp)) which.cp <- match.arg(which.cp)
    type.eff <- match.arg(type.eff)
    
    ## from 1.3-4: added option of external efficiencies or threshold cycles,
    ## either single value (recycled) or a vector of values.
    ## Check if all efficiencies are in [1, 2]
    ## if all is o.k., add to 'pcrbatch' data
    if (is.numeric(which.eff)) {
      if (length(which.eff) == 1) which.eff <- rep(which.eff, ncol(data) - 1)
      else {
        if (length(which.eff) != ncol(data) - 1) stop("Length of input efficiencies does not match number of runs!")
      }
      if (!all(which.eff >= 1 & which.eff <= 2)) stop("All efficiencies must be in [1, 2]. Consider adding 1 to each value, i.e. 'EFF <- EFF + 1'")          
      effDAT <- matrix(c("ext.eff", which.eff), nrow = 1)
      colnames(effDAT) <- colnames(data)
      data <- rbind(data, effDAT) 
      which.eff <- "ext"     
    }
      
    if (is.numeric(which.cp)) {
      if (length(which.cp) != ncol(data) - 1 ) stop("Length of input threshold cycles does not match number of runs!")
      cpDAT <- matrix(c("sig.ext", which.cp), nrow = 1)
      colnames(cpDAT) <- colnames(data)
      data <- rbind(data, cpDAT) 
      which.cp <- "ext"      
    }      
    
    ## from 1.3-4: added mak3 parameter
    if (which.eff == "mak") {
      if (length(grep("mak3", data[, 1])) == 0) stop("'data' has no mak3 model included! Please use 'pcrbatch' with 'do.mak = TRUE'!")
      WHICH <- which(data[, 1] == "mak3.D0")
      data[, 1] <- sub("mak3.D0", "mak.eff", data[, 1])
    }    
    
    DATA <- data[, -1]

    ## added removal of failed runs (either failed fits
    ## or SOD outlier) from DATA and 'group' by identification
    ## of *...* or **...** in sample name in version 1.3-5
    sampNAMES <- names(DATA)
    hasTag <- grep("\\*[[:print:]]*\\*", sampNAMES)
    if (length(hasTag) > 0) {
       DATA <- DATA[, -hasTag]
       group <- group[-hasTag]
    }
    
    cpNAMES <- effNAMES <- NULL
    PATTERN <- unique(group)
    
    ## test for presence of reference genes
    if (all(regexpr("rs", group, perl = TRUE) == -1)) refNo <- TRUE else refNo <- FALSE   
    
    ## check for replicate data, if not present set type.eff = "individual"
    REPS <- lapply(PATTERN, function(x) which(x == group))
    NREPS <- sapply(REPS, function(x) length(x))           
    if (!all(NREPS > 1)) type.eff <- "individual"
    
    ## initialize data as time-series
    cpDAT <- effDAT <- ts()     
    
    ## for all entries 'gs', 'gc', 'rs', 'rc' do...
    for (i in 1:length(PATTERN)) {
      WHICH <- which(group == PATTERN[i]) 
      cpSEL <- which(data[, 1] == paste("sig.", which.cp, sep = ""))
      effSEL <- which(data[, 1] == paste(which.eff, ".eff", sep = ""))
      if (length(WHICH) != 1) tempCP <- as.numeric(t(DATA[cpSEL, WHICH])) else tempCP <- as.numeric(as.vector(DATA[cpSEL, WHICH]))     
      if (length(WHICH) != 1) tempEff <- as.numeric(t(DATA[effSEL, WHICH])) else tempEff <- as.numeric(as.vector(DATA[effSEL, WHICH]))
      if (is.numeric(which.eff)) tempEff <- rep(which.eff, length(WHICH))
      cpDAT <- cbind(cpDAT, ts(tempCP))            
      effDAT <- cbind(effDAT, ts(tempEff))
      cpNAMES <- c(cpNAMES, paste("cp.", PATTERN[i], sep = ""))
      effNAMES <- c(effNAMES, paste("eff.", PATTERN[i], sep = ""))         
    }
        
    cpDAT <- cpDAT[, -1]      
    effDAT <- effDAT[, -1]         
    
    if (is.numeric(which.eff))
        type.eff <- "individual"
    if (type.eff == "mean.single")
        effDAT <- t(replicate(nrow(effDAT), apply(effDAT, 2, function(x) mean(x, na.rm = TRUE))))
    if (type.eff == "median.single")
        effDAT <- t(replicate(nrow(effDAT), apply(effDAT, 2, function(x) median(x, na.rm = TRUE))))
    if (type.eff == "mean.pair") {
        effDAT[, 1:2] <- mean(effDAT[, 1:2], na.rm = TRUE)
        if (!refNo)
            effDAT[, 3:4] <- mean(effDAT[, 3:4], na.rm = TRUE)
    }
    if (type.eff == "median.pair") {
        effDAT[, 1:2] <- median(effDAT[, 1:2], na.rm = TRUE)
        if (!refNo)
            effDAT[, 3:4] <- median(effDAT[, 3:4], na.rm = TRUE)
    }
    
    cpDAT <- matrix(cpDAT, ncol = length(cpNAMES))
    effDAT <- matrix(effDAT, ncol = length(effNAMES))
    
    allDAT <- cbind(cpDAT, effDAT)     
    colnames(allDAT) <- c(cpNAMES, effNAMES)     
            
    if (refNo) {
      EXPR <- expression(eff.gc^cp.gc/eff.gs^cp.gs)
      TIES <- c(1, 2, 1, 2)  
      
      ## added mak3 option (1.3-4) => we only need D0 for ratio calculation
      if (which.eff == "mak") {
        EXPR <- expression(eff.gs/eff.gc)
        TIES <- NULL
      }      
    }    
    else {
      EXPR <- expression((eff.gc^cp.gc/eff.gs^cp.gs)/(eff.rc^cp.rc/eff.rs^cp.rs))
      TIES <- c(1, 2, 1, 2, 1, 2, 1, 2)
      
      ## added mak3 option (1.3-4) => we only need D0 for ratio calculation
      if (which.eff == "mak") {
        EXPR <- expression((eff.gs/eff.gc)/(eff.rs/eff.rc))
        TIES <- NULL
      }      
    }       
            
    CRIT <- c("perm > init", "perm == init", "perm < init")
              
    PROP <- try(propagate(EXPR, allDAT, do.sim = TRUE, do.perm = TRUE, ties = TIES, perm.crit = CRIT, 
                      verbose = TRUE, logx = TRUE, ...))

    if (inherits(PROP, "try-error")) stop("'propagate' failed to calculate ratios! Try other 'which.eff', 'type.eff' or 'which.cp'.")
    PROP <- c(list(data = allDAT), PROP)     
    class(PROP) <- "ratiocalc"
    return(PROP)
}
