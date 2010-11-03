batchstat <- function(..., group = NULL,  do = c("cbind", "stat"), statfun = mean)
{
  do <- match.arg(do)
  DATA <- list(...)
  lapply(DATA, function(x) if (class(x)[2] != "pcrbatch") stop("All data must be of class 'pcrbatch'!"))
  lapply(DATA, function(x) if (do == "stat" && ncol(x) - 1 != length(group)) stop("'group' length and number of runs must match!"))
    
  if (do == "cbind") {
    DATAout <- NULL
    aN <- unique(matrix(sapply(DATA, function(x) x[, 1]), ncol = 1))
    DATAout$Vars <- aN
    for (i in 1:length(DATA)) {
      DATAout <- merge(DATAout, DATA[[i]], by.x = "Vars", by.y = "Vars")
    }   
  }

  if (do == "stat") {
    if (is.null(group)) stop("Please define 'group'ing vector!")
    group <- as.factor(group)
    DATAout <- list()

    for (i in 1:length(DATA)) {
      anno <- DATA[[i]][ ,1, drop = FALSE]
      DATAtemp <- DATA[[i]][, -1, drop = FALSE]
      STAT <- apply(DATAtemp, 1, function(x) tapply(as.numeric(x), group, function(y) statfun(y, na.rm = TRUE)))
      STAT <- t(STAT)
      STAT <- cbind(anno, STAT)
      colnames(STAT) <- c("Vars", paste("group", 1:nlevels(group), sep = ""))
      DATAout[[i]] <- STAT
    }
  }
  return(DATAout)
}