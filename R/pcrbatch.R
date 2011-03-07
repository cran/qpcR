pcrbatch <- function(
x, 
cyc = 1, 
fluo = NULL, 
model = l4, 
check = "uni2",
checkPAR = parKOD(),
remove = c("none", "fit", "KOD"),
type = "cpD2",
labels = NULL, 
norm = FALSE,
backsub = NULL,
smooth = c("none", "smooth", "lowess", "supsmu", "spline"), 
smoothPAR = list(span = 0.1), 
factor = 1,
opt = FALSE,
optPAR = list(sig.level = 0.05, crit = "ftest"),
do.mak = FALSE,
group = NULL,
names = c("group", "first"),
plot = TRUE,
verbose = TRUE,
...) 
{
  remove <- match.arg(remove)
  smooth <- match.arg(smooth)  
  names <- match.arg(names)
        
  ## make initial 'modlist'
  if (class(x) != "modlist") {  
    if (names(x)[cyc] != "Cycles") stop("Column 1 should be named 'Cycles'!")
    modLIST <- try(modlist(x = x, cyc = cyc, fluo = fluo, model = model, check = check, checkPAR = checkPAR,
                       remove = remove, labels = labels, norm = norm, backsub = backsub, smooth = smooth,
                       smoothPAR = smoothPAR, factor = factor, opt = opt, optPAR = optPAR, 
                       verbose = verbose, ...), silent = TRUE)   
    if (inherits(modLIST, "try-error")) stop("There was an error during 'modlist' creation.")
  } else modLIST <- x 
                           
  ## if 'group' is defined, make a 'replist'
  if (!is.null(group)) {     
    repLIST <- try(replist(modLIST, group = group, check = check, checkPAR = checkPAR, remove = remove,
                           names = names, opt = opt, optPAR = optPAR, verbose = TRUE, ...), silent = TRUE) 
    if (inherits(repLIST, "try-error")) cat("There was an error during 'replist' creation. Continuing with original 'modlist'...\n")
    else modLIST <- repLIST
  }  
                              
  ## plot diagnostics, if selected
  if (plot) plot(modLIST, which = "single")
  outLIST <- vector("list", length = length(modLIST))
  
  ## for all single models in the 'modlist' do...
  for (i in 1:length(modLIST)) {    
    NAME <- modLIST[[i]]$name    
    fitOBJ <- modLIST[[i]]   
    
    cat("Analyzing", NAME, "...\n")
    flush.console()
    
    ## sigmoidal model
    cat("  Calculating 'eff' and 'ct' from sigmoidal model...\n")
    flush.console()
    EFF <- try(efficiency(fitOBJ, plot = FALSE, type = type, ...), silent = TRUE)
    if (!inherits(EFF, "try-error")) EFF <- c(EFF, coef(fitOBJ), model = fitOBJ$MODEL$name) else EFF <- list(eff = NA)
    names(EFF) <- paste("sig.", names(EFF), sep = "") 
        
    ## sliding window method
    cat("  Using window-of-linearity...\n")
    SLI <- try(sliwin(fitOBJ, plot = FALSE, ...)[1:3], silent = TRUE)
    if (inherits(SLI, "try-error")) SLI <- list(eff = NA) 
    names(SLI) <- paste("sli.", names(SLI), sep = "")       
            
    ## exponential model
    cat("  Fitting exponential model...\n")
    EXP <- try(expfit(fitOBJ, plot = FALSE, ...)[-c(2, 4, 9)], silent = TRUE)
    if (inherits(EXP, "try-error")) EXP <- list(eff = NA)
    names(EXP) <- paste("exp.", names(EXP), sep = "")
    
    ## from 1.3-4: if MAK model selected, make model to attach, added
    if (do.mak) {
      cat("  Fitting mak3 model...\n")
      MAK <- try(pcrfit(fitOBJ$DATA, 1, 2, mak3, verbose = FALSE), silent = TRUE)
      if (inherits(MAK, "try-error")) {
        cat("There was an error in buidling the mak3 model. Continuing without...\n")
        flush.console()
        MAK <- list(D0 = NA)
      } else {
        MAK <- coef(MAK)
        names(MAK) <- paste("mak3.", names(MAK), sep = "")
      }  
    } else MAK <- NULL
           
    outALL <- c(EFF, SLI, EXP, MAK)     
    outLIST[[i]] <- outALL
  }
  
  allNAMES <- unique(unlist(lapply(outLIST, function(x) names(x))))
  resMAT <- matrix(nrow = length(allNAMES), ncol = length(outLIST) + 1)
  resMAT <- as.data.frame(resMAT)
  resMAT[, 1] <- allNAMES
     
  ## aggregate all results into a dataframe by 'merge'
  for (i in 1:length(outLIST)) {
    tempDAT <- t(as.data.frame(outLIST[[i]]))
    m <- match(resMAT[, 1], rownames(tempDAT))
    resMAT[, i + 1] <- tempDAT[m, ]
  }
  
  colnames(resMAT)[1] <- "Vars"
  names(resMAT)[-1] <-  sapply(modLIST, function(x) x$name)       
  cat("Writing to clipboard...\n\n")
  write.table(resMAT, file = "clipboard-64000", sep = "\t", row.names = FALSE)
  class(resMAT) <- c("data.frame", "pcrbatch")   
  return(resMAT)
}
