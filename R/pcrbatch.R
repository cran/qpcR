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
  smooth <- match.arg(smooth)
  ave <- match.arg(ave)
  outList <- list()         
  
  if (class(x) != "modlist" && is.null(cols)) cols <- 2:ncol(x)
  if (class(x) != "modlist" && names(x)[1] != "Cycles") stop("Column 1 should be 'Cycles'!")
  if (!is.null(backsub) && !is.numeric(backsub)) stop("'backsub' must be either NULL or a numeric sequence!")
  if (!is.null(cols) && min(cols) == 1) stop("'cols' must be > 1 because Column 1 must be 'Cycles'!")

  if (!is.null(group)) {
    if (class(x) == "modlist") stop("This is a 'modlist'. Use 'replist' for averaging...")
    group <- as.factor(group)
    if (length(group) != length(cols)) stop("replicates and column numbers do not match!")
    if (ave == "mean") centre <- function(x) mean(x, na.rm = TRUE)
    if (ave == "median") centre <- function(x) median(x, na.rm = TRUE)
    Cycles <- x[, 1]
    DATA <- x[, cols]     
    DATA2 <- apply(DATA, 1, function(x) tapply(x, group, function(x) centre(x)))
    if (nlevels(group) > 1) DATA2 <- t(DATA2)
    DATA2 <- as.data.frame(DATA2)
    namevec <- paste("group_", 1:length(levels(group)), sep = "")
    colnames(DATA2) <- namevec
    x <- cbind(Cycles, DATA2)
    cols <- 2:ncol(x)
  }
    
  if (class(x) == "modlist") {
      dataList <- lapply(x, function(x) x$DATA)
      modelList <- lapply(x, function(x) x$MODEL)
      nameList <- lapply(x, function(x) x$names)
  }
  else {
      dataList <- list()
      modelList <- list()
      nameList <- list()
      for (i in 1:length(cols)) {
            CC <- complete.cases(x[, cols[i]])
            dataList[[i]] <- cbind(x[CC, 1], x[CC, cols[i]])
            modelList[[i]] <- model
            nameList[[i]] <- colnames(x)[cols[i]]
      }
  }

  for (i in 1:length(dataList)) {
    Cycles <- dataList[[i]][, 1]
    Fluo <- dataList[[i]][, 2]
    
    if (fact != 1) Fluo <- Fluo * fact
    if (smooth == "tukey") Fluo <- smooth(Fluo)
    if (smooth == "lowess") Fluo <- lowess(Fluo, f = 0.1)$y
    if (norm == TRUE) {
      Fluo <- Fluo - min(Fluo, na.rm = TRUE)
      Fluo <- Fluo/max(Fluo, na.rm = TRUE)
    }
    if (!is.null(backsub)) {
      back <- mean(Fluo[backsub], na.rm = TRUE)
      Fluo <- Fluo - back
    }

    DATA <- data.frame(Cycles = Cycles, Fluo = Fluo)
    
    flush.console()
    cat("Processing ", nameList[[i]], "...\n", sep = "")
    cat("   Building sigmoidal model (", modelList[[i]]$name, ")...", sep = "")
    fitObj <- try(pcrfit(DATA, 1, 2, model = modelList[[i]], ...), silent = TRUE)
    if (inherits(fitObj, "try-error")) {
      cat(" => There was an error in sigmoidal fitting. Skipping...\n") 
      next
    }
    
    if (missing(crit)) crit = "ftest" else crit <- crit
    if (opt) {
 	    fitObj2 <- try(mselect(fitObj, verbose = FALSE, crit = crit, ...))
      if (inherits(fitObj2, "try-error")) {
        fitObj <- fitObj
        cat(" => ", nameList[[i]], " gave a model selection error!", sep = "")
      } else {
        fitObj <- fitObj2
        cat(" => ", fitObj$MODEL$name, sep = "")
      }
    }
    
    cat("\n")
    
    EFF <- try(efficiency(fitObj, plot = FALSE, type = type, ...), silent = TRUE)
    if (retPar) EFF <- c(EFF, coef(fitObj))
    names(EFF) <- paste("sig.", names(EFF), sep = "")

    cat("   Using window-of-linearity...\n")
    SLI <- try(sliwin(fitObj, plot = FALSE, ...), silent = TRUE)
    names(SLI) <- paste("sli.", names(SLI), sep = "")       
            
    cat("   Fitting exponential model...\n")
    EXP <- try(expfit(fitObj, plot = FALSE, ...)[-c(2, 4, 9)], silent = TRUE)
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
      resMat[which(resMat[, 1] == "sig.model"), i + 1] <- fitObj$MODEL$name
    }

    namevec <- unlist(nameList)
    colnames(resMat) <- c("Vars", namevec)
    
    cat("Writing to clipboard...\n\n")
    write.table(resMat, file = "clipboard-64000", sep = "\t", row.names = FALSE)
    class(resMat) <- "pcrbatch"
    return(resMat)
}
