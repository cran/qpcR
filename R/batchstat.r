batchstat <- function(..., group = NULL,  do = c("cbind", "stat"), statfun = mean)
{
  do <- match.arg(do)
  DATA <- list(...)
  lapply(DATA, function(x) if (class(x) != "pcrbatch") stop("All data must be of class 'pcrbatch'!"))
  lapply(DATA, function(x) if (do == "stat" && ncol(x) - 1 != length(group)) stop("'group' length and number of runs must match!"))
  
  
  if (do == "cbind") {
    DATAout <- NULL
    anno <- DATA[[1]][, 1]
    for (i in 1:length(DATA)) {
     if (i == 1) DATAout <- cbind(DATAout, DATA[[i]])
      else DATAout <- cbind(DATAout, DATA[[i]][, -1])
     }
  }

  if (do == "stat") {
    if (is.null(group)) stop("Please define 'group'ing vector!")
    group <- as.factor(group)
    DATAout <- list()

    for (i in 1:length(DATA)) {
      anno <- DATA[[i]][ ,1]
      DATAtemp <- DATA[[i]][, -1]
      STAT <- apply(DATAtemp, 1, function(x) tapply(as.numeric(x), group, function(y) statfun(y, na.rm = TRUE)))
      if (nlevels(group) == 1) STAT <- t(t(STAT)) else STAT <- t(STAT)
      STAT <- cbind(anno, STAT)
      colnames(STAT) <- c("Vars", paste("group", 1:nlevels(group), sep = ""))
      DATAout[[i]] <- STAT
    }
  }
  return(DATAout)
}