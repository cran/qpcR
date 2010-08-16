ratiobatch <- function(
data, 
group = NULL, 
plot = TRUE, 
combs = c("same", "across", "all"),
type.eff = "mean.single",
which.cp = "cpD2",
which.eff = "sli",
dataout = "clip", 
...)
{
  combs <- match.arg(combs)
  if (all(class(data) != "pcrbatch")) stop("data must be of class 'pcrbatch'!")
  if (is.null(group)) {
    group <- colnames(data)[-1]
    CLASS <- sapply(group, function(x) class(x))
    if (!all(CLASS == "character")) stop("'group' definition must be of class 'character' (i.e. r1g1)")
  } 
  ANNO <- data[, 1]
  DATA <- data[, -1]    
  if (length(group) != ncol(DATA)) stop("'group' vector and 'data' columns are not of same length!")
  group <- sub("\\.\\d*", "", group, perl = TRUE)        
    
  ## detect r*c*, g*c*, r*s* and g*s* in 'group'
  RCs <- sort(unique(grep("r\\d*c\\d*", group, perl = TRUE, value = TRUE)))
  GCs <- sort(unique(grep("g\\d*c\\d*", group, perl = TRUE, value = TRUE)))
  RSs <- sort(unique(grep("r\\d*s\\d*", group, perl = TRUE, value = TRUE)))
  GSs <- sort(unique(grep("g\\d*s\\d*", group, perl = TRUE, value = TRUE))) 
  
  ## detect absence of reference runs
  if (length(RCs) > 0 && length(RSs) > 0) hasRef <- TRUE else hasRef <- FALSE
      
  ## check if reference and control genes are used equally in controls and treatment samples
  RefInCon <- unique(as.numeric(sub("r(\\d*)c\\d*", "\\1", RCs, perl = TRUE)))
  RefInSamp <- unique(as.numeric(sub("r(\\d*)s\\d*", "\\1", RSs, perl = TRUE)))
  GoiInCon <- unique(as.numeric(sub("g(\\d*)c\\d*", "\\1", GCs, perl = TRUE)))
  GoiInSamp <- unique(as.numeric(sub("g(\\d*)s\\d*", "\\1", GSs, perl = TRUE)))    
  if (!all(RefInCon == RefInSamp)) stop("Unequal number of reference genes in controls and treatment samples!")
  if (!all(GoiInCon == GoiInSamp)) stop("Unequal number of genes-of-interest in controls and treatment samples!")
    
  ## do combinations in presence/absence of reference genes
  if (hasRef) COMBS <- expand.grid(RCs, GCs, RSs, GSs, stringsAsFactors = FALSE)   
  else COMBS <- expand.grid(GCs, GSs, stringsAsFactors = FALSE)   
  
  ## remove 'nonpair' combinations, i.e. r1s2:r1s1    
  Cnum <- t(apply(COMBS, 1, function(x) gsub("[rgs]\\d*", "", x, perl = TRUE)))      
  Snum <- t(apply(COMBS, 1, function(x) gsub("[rgc]\\d*", "", x, perl = TRUE)))
  Gnum <- t(apply(COMBS, 1, function(x) gsub("[rsc]\\d*", "", x, perl = TRUE)))
  Rnum <- t(apply(COMBS, 1, function(x) gsub("[gsc]\\d*", "", x, perl = TRUE)))              
      
  Cnum <- t(apply(Cnum, 1, function(x) x[x != ""]))
  Snum <- t(apply(Snum, 1, function(x) x[x != ""]))
  Gnum <- t(apply(Gnum, 1, function(x) x[x != ""]))
  Rnum <- t(apply(Rnum, 1, function(x) x[x != ""]))  
      
  if (nrow(Cnum) == 1) Cnum <- t(Cnum)
  if (nrow(Snum) == 1) Snum <- t(Snum)
  if (nrow(Gnum) == 1) Gnum <- t(Gnum)
  if (nrow(Rnum) == 1) Rnum <- t(Rnum)      
      
  if (hasRef) {     
    if (combs == "across") {
      SELECT <- which(Cnum[, 1] == Cnum[, 2] & Snum[, 1] == Snum[, 2])
      COMBS <- COMBS[SELECT, ]
    } else 
    if (combs == "same") {
      SELECT <- which(Cnum[, 1] == Cnum[, 2] & Snum[, 1] == Snum[, 2] & Gnum[, 1] == Gnum[, 2] & Rnum[, 1] == Rnum[, 2])  
      COMBS <- COMBS[SELECT, ]
    } 
  } else {
    if (combs == "same") {
      SELECT <- which(Gnum[, 1] == Gnum[, 2]) 
      COMBS <- COMBS[SELECT, ]
    }        
  }         
       
  ncomb <- nrow(COMBS)
  outLIST <- list()
  nameLIST <- list()             
  counter <- 1      
  
  ## take combinations into dataframe  
  for (i in 1:nrow(COMBS)) {
    if (hasRef) {
      RCmatch <- grep(COMBS[i, 1], group, perl = TRUE)         
      RCdat <- as.data.frame(DATA[, RCmatch]) 
      GCmatch <- grep(COMBS[i, 2], group, perl = TRUE)     
      GCdat <- as.data.frame(DATA[, GCmatch])
      RSmatch <- grep(COMBS[i, 3], group, perl = TRUE)     
      RSdat <- as.data.frame(DATA[, RSmatch])
      GSmatch <- grep(COMBS[i, 4], group, perl = TRUE)     
      GSdat <- as.data.frame(DATA[, GSmatch])          
      finalDATA <- cbind(ANNO, RCdat, GCdat, RSdat, GSdat)      
    } else {          
      GCmatch <- grep(COMBS[i, 1], group, perl = TRUE)
      GCdat <- as.data.frame(DATA[, GCmatch])
      GSmatch <- grep(COMBS[i, 2], group, perl = TRUE)     
      GSdat <- as.data.frame(DATA[, GSmatch])          
      finalDATA <- cbind(ANNO, GCdat, GSdat) 
    }              
    finalNAME <-  as.vector(unlist(COMBS[i, ]))
    finalNAME <- paste(finalNAME, collapse = ":")      
    cat("Calculating ", finalNAME, " (", counter, " of ", ncomb, ")...\n", sep = "")
    flush.console()
    class(finalDATA) <- c("data.frame", "pcrbatch")
    if (hasRef) finalGROUP <- c(rep("rc", ncol(RCdat)), rep("gc", ncol(GCdat)), rep("rs", ncol(RSdat)), rep("gs", ncol(GSdat)))
    else finalGROUP <- c(rep("gc", ncol(GCdat)), rep("gs", ncol(GSdat)))
    
    ## do ratio calculation for all combinations
    outLIST[[counter]] <- ratiocalc(finalDATA, finalGROUP, plot = plot, type.eff = type.eff, 
                                    which.cp = which.cp, which.eff = which.eff, ...)$summary
    nameLIST[[counter]] <- finalNAME
    counter <- counter + 1 
  }  
  
  names(outLIST) <- nameLIST
  outFRAME <- sapply(outLIST, function(x) as.matrix(x, ncol = 1))
  rowNAMES <- c(paste(rownames(outLIST[[1]]), "Sim", sep = "."), 
                paste(rownames(outLIST[[1]]), "Perm", sep = "."), 
                paste(rownames(outLIST[[1]]), "Prop", sep = ".")) 
  rownames(outFRAME) <- rowNAMES 
  outFRAME <- outFRAME[complete.cases(outFRAME), ] 
  
  if (plot && length(outLIST) < 50) {
    DIM <- ceiling(sqrt(length(outLIST)))
    par(mfrow = c(DIM, DIM + 1)) 
    par(mar = c(1, 3, 2, 2))
    for (i in 1:length(outLIST)) {
      barplot(as.numeric(outLIST[[i]][1, ]), col = c("darkblue", "darkred", "darkgreen"), 
              main = nameLIST[[i]], cex.main = 0.8, log = "y", cex.axis = 0.8, las = 2)
    }
    par(mar = c(0.5, 0.5, 0.5, 0.5))
    plot(1, 1, type = "n", axes = FALSE)
    legend(x = 0.7, y = 1.4, legend = c("Mean of \nMonte-Carlo Sim.", "Mean of \nPermutation", "Mean of \nError Propagation"), 
           bty = "n", cex = 0.8, fill  = c("darkblue", "darkred", "darkgreen"), y.intersp = 2)
  } 
  outFRAME2 <- cbind(VALS = rownames(outFRAME), outFRAME)
  write.table(outFRAME2, file = ifelse(dataout == "clip", "clipboard-64000", dataout), row.names = FALSE, sep = "\t")
  return(list(resList = outLIST, resDat = outFRAME))     
}

