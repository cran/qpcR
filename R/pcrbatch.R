pcrbatch <- function(
x,
cyc = 1,
fluo = NULL, 
model = l4, 
remove = FALSE, 
type = "cpD2", 
opt = FALSE, 
norm = FALSE,
backsub = NULL,
smooth = c("none", "tukey", "lowess", "supsmu", "spline"), 
span = 0.1, 
factor = 1, 
do.mak = FALSE,
opt.method =  "all",
nls.method = "all",
sig.level = 0.05, 
crit = "ftest", 
group = NULL,
plot = TRUE,
...) 
{
  ### make initial 'modlist'
  if (class(x) != "modlist") {
    if (is.null(fluo)) fluo <- 2:ncol(x)
    if (names(x)[cyc] != "Cycles") stop("Column 1 should be named 'Cycles'!")
    cat("Creating modlist...\n\n")
    modLIST <- modlist(x, cyc = cyc, fluo = fluo, model = model, remove = remove, opt = opt, norm = norm, backsub = backsub, 
                       smooth = smooth, span = span, factor = factor, opt.method = opt.method, nls.method = nls.method, sig.level = sig.level, 
                       crit = crit, ...)
  } else modLIST <- x
  
  ### if 'group' is defined, make a 'replist'
  if (!is.null(group)) {     
    if (length(group) != length(modLIST)) stop("Length of 'group' and sample number do not match!")
    repLIST <- try(replist(modLIST, group = group, opt = opt, verbose = TRUE, ...), silent = TRUE)
    if (inherits(repLIST, "try-error")) cat("There was an error during 'replist'. Continuing with original 'modlist'...")
      else modLIST <- repLIST
  }  
  
  flush.console()  
    
  ### plot diagnostics, if selected
  if (plot) plot(modLIST, which = "single")
  outLIST <- list() 
  
  ### for all single models in the 'modlist' do...
  for (i in 1:length(modLIST)) {    
    NAME <- modLIST[[i]]$name
    
    fitOBJ <- modLIST[[i]]
    
    cat("Analyzing", NAME, "...\n")
    
    ### sigmoidal model 
    cat("  Calculating 'eff' and 'ct' from sigmoidal model...\n")
    flush.console()
    EFF <- try(efficiency(fitOBJ, plot = FALSE, type = type, ...), silent = TRUE)
    if (!inherits(EFF, "try-error")) EFF <- c(EFF, coef(fitOBJ), model = fitOBJ$MODEL$name) else EFF <- list(eff = NA)
    names(EFF) <- paste("sig.", names(EFF), sep = "")       
    
    ### sliding window method
    cat("  Using window-of-linearity...\n")
    SLI <- try(sliwin(fitOBJ, plot = FALSE, ...), silent = TRUE)
    if (inherits(SLI, "try-error")) SLI <- list(eff = NA) 
    names(SLI) <- paste("sli.", names(SLI), sep = "")       
            
    ### exponential model
    cat("  Fitting exponential model...\n")
    EXP <- try(expfit(fitOBJ, plot = FALSE, ...)[-c(2, 4, 9)], silent = TRUE)
    if (inherits(EXP, "try-error")) EXP <- list(eff = NA)
    names(EXP) <- paste("exp.", names(EXP), sep = "")
    
    ### if MAK model selected, make model to attach
    if (do.mak) {
      cat("  Fitting mak3 model...\n")
      MAK <- try(pcrfit(fitOBJ$DATA, 1, 2, mak3), silent = TRUE)
      if (inherits(MAK, "try-error")) {
        cat("There was an error in buidling the mak3 model. Continuing without...\n")
        flush.console()
        MAK <- list(D0 = NA)
      } else {
        MAK <- coef(MAK)
        names(MAK) <- paste("mak3.", names(MAK), sep = "")
      }  
    } else MAK <- NULL
    
    cat("\n")  
   
    outALL <- c(EFF, SLI, EXP, MAK)     
    outLIST[[i]] <- outALL
  }     
    
  allNAMES <- unique(unlist(lapply(outLIST, function(x) names(x))))  
  resMAT <- data.frame(ROWNAMES = allNAMES) 
     
  ### aggregate all results into a dataframe by 'merge'
  for (i in 1:length(outLIST)) {
    tempDAT <- t(as.data.frame(outLIST[[i]]))       
    tempDAT <- cbind(ROWNAMES = rownames(tempDAT), tempDAT)       
    resMAT <- merge(resMAT, tempDAT, by.x = "ROWNAMES", by.y = "ROWNAMES", sort = FALSE, all = TRUE)         
  }   

  names(resMAT)[1] <- "Vars"  
  names(resMAT)[-1] <-  sapply(modLIST, function(x) x$name)       
  cat("Writing to clipboard...\n\n")
  write.table(resMAT, file = "clipboard-64000", sep = "\t", row.names = FALSE)
  class(resMAT) <- c("data.frame", "pcrbatch")   
  return(resMAT)
}
