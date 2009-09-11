pcrbatch <- function(
x, 
cols = NULL, 
model = l4, 
group = NULL, 
type = "cpD2", 
opt = FALSE, 
smooth = c("none", "tukey", "lowess"), 
norm = FALSE, 
fact = 1, 
ave = c("mean", "median"), 
backsub = NULL,   
retPar = FALSE,
crit, 
...) 
{
  if (class(x) != "modlist" && is.null(cols)) cols <- 2:ncol(x)
  if (class(x) != "modlist" && names(x)[1] != "Cycles") stop("Column 1 should be 'Cycles'!")
  smooth <- match.arg(smooth)             
  ave <- match.arg(ave)
  outList <- list()         
    
  if (!is.null(backsub) && !is.numeric(backsub)) stop("'backsub' must be either NULL or a numeric sequence!")
  if (!is.null(cols) && min(cols) == 1) stop("'cols' must be > 1 because Column 1 must be 'Cycles'!")   
    
  if (class(x) == "modlist") {
    wdata <- sapply(x, function(x) x$DATA[, 2])
    namevec <- sapply(x, function(x) x$names)  
    models <- lapply(x, function(x) x$MODEL)    
    colnames(wdata) <- namevec
    Cycles <- x[[1]]$DATA[, 1]
  } else {
    wdata <- as.data.frame(x[, cols])
    namevec <- colnames(x)[cols] 
    models <- rep(list(model), ncol(wdata))        
    Cycles <- x[, 1]
  }
           
  if (!is.null(group)) {
    group <- as.factor(group)
    if (length(group) != length(cols)) stop("replicates and column numbers do not match!")
    data.pre <- wdata
    
    if (ave == "mean") centre <- function(x) mean(x, na.rm = TRUE)
    if (ave == "median") centre <- function(x) median(x, na.rm = TRUE)
    data.post <- apply(data.pre, 1, function(x) tapply(x, group, function(x) centre(x)))
    
    if (nlevels(group) > 1) data.post <- t(data.post)
    wdata <- as.data.frame(data.post)
    namevec <- paste("group", 1:length(levels(group)), sep = "")
  }
      
  for (i in 1:ncol(wdata)) {
    data <- wdata[, i] * fact
    
    if (smooth == "tukey") data <- smooth(data)
    if (smooth == "lowess") data <- lowess(data, f = 0.1)$y
    if (norm == TRUE) {
      data <- data - min(data, na.rm = TRUE)
      data <- data/max(data, na.rm = TRUE)
    }
    if (!is.null(backsub)) {
      back <- mean(data[backsub], na.rm = TRUE)
      data <- data - back
    }    
    
    l1 <- length(Cycles)
    l2 <- length(data)
    maxl <- max(c(l1, l2), na.rm = TRUE) 
    
    mat <- data.frame(cbind(Cycles[1:maxl], data[1:maxl]))
    FIT <- try(pcrfit(mat, 1, 2, model = models[[i]], ...), silent = TRUE)     
    
    if (missing(crit)) CRIT = "ftest" else CRIT <- crit
    if (opt) FIT <- mselect(FIT, verbose = FALSE, crit = CRIT)
       
    fctName <-  FIT$MODEL$name
        
    flush.console()       
    cat("Processing ", namevec[i], "...\n", sep = "") 
    cat("   Building sigmoidal model (", fctName, ")...\n", sep = "")
    EFF <- try(efficiency(FIT, plot = FALSE, type = type, ...), silent = TRUE) 
    if (retPar) EFF <- c(EFF, coef(FIT))            
    names(EFF) <- paste("sig.", names(EFF), sep = "")     
    
    cat("   Using window-of-linearity...\n")
    SLI <- try(sliwin(FIT, plot = FALSE, ...), silent = TRUE)
    names(SLI) <- paste("sli.", names(SLI), sep = "")       
            
    cat("   Fitting exponential model...\n")
    EXP <- try(expfit(FIT, plot = FALSE, ...)[-c(2, 4, 9)], silent = TRUE)
    names(EXP) <- paste("exp.", names(EXP), sep = "")
   
    out.all <- c(EFF, list(sig.model = NA), SLI, EXP)
    outList[[i]] <- out.all   
  }         
    
    allNAMES <- c("sig.eff", "sig.resVar", "sig.AICc", "sig.AIC", "sig.Rsq", "sig.Rsq.ad", "sig.cpD1", "sig.cpD2", "sig.cpE",
                  "sig.cpR", "sig.cpT", "sig.Cy0", "sig.fluo", "sig.init1", "sig.init2", "sig.cf", "sig.model", "sli.eff",
                  "sli.rmax", "sli.init", "exp.point", "exp.eff", "exp.AIC", "exp.resVar", "exp.RMSE", "exp.init")
                  
    resMat <- matrix(nrow = length(allNAMES), ncol = length(outList) + 1)
    resMat[, 1] <- allNAMES
    
    for (i in 1:length(outList)) {
      outVALS <- as.numeric(outList[[i]])
      outNAMES <- names(outList[[i]])
      m <- match(outNAMES, allNAMES)
      isNA <- which(is.na(m))
      m[isNA] <- isNA
      resMat[m, i + 1] <- outVALS
      resMat[which(resMat[, 1] == "sig.model"), i + 1] <- fctName
    }

    colnames(resMat) <- c("Vars", namevec)
    
    cat("Writing to clipboard...\n\n")
    write.table(resMat, file = "clipboard-64000", sep = "\t", row.names = FALSE)
    class(resMat) <- "pcrbatch"
    return(resMat)
}
