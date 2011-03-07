getPar <- function(x, cp = "cpD2", eff = "sigfit", ...) {
  options(expressions = 20000)  
  if (class(x)[1] == "pcrfit") x <- modlist(x)
  if (class(x)[1] != "modlist") stop("'x' must be either of class 'pcrfit' or 'modlist'!")

  CP <- numeric(length(x))
  EFF <- numeric(length(x))
  NAMES <- sapply(x, function(a) a$names)

  for (i in 1:length(x)) {
    qpcR:::counter(i)
    flush.console()
    tempMod <- x[[i]]
    outName <- switch(cp, "cpD2" = "cpD2", "cpD1" = "cpD1", "maxE" = "cpE", "expR" = "cpR", "Cy0" = "Cy0", "CQ" = "cpCQ", "maxRatio" = "cpMR", stop())
    tempRes <- tryCatch(efficiency(tempMod, plot = FALSE, type = cp, ...), error = function(e) NA)
    if (!is.na(tempRes)) CP[i] <- tempRes[[outName]] else CP[i] <- NA
    EFF[i] <- switch(eff, "sigfit" =  if (!is.na(tempRes)) tempRes$eff else NA,
                          "expfit" = tryCatch(expfit(tempMod, plot = FALSE, ...)$eff, error = function(e) NA),
                          "sliwin" = tryCatch(sliwin(tempMod, plot = FALSE, ...)$eff, error = function(e) NA))
  }

  res <- rbind.na(CP, EFF)
  colnames(res) <- NAMES
  rownames(res) <- c(paste("sig.", outName, sep = ""), paste(substr(eff, 1, 3), ".eff", sep = ""))
  write.table(res, file = "clipboard", sep = "\t", row.names = FALSE)
  return(res)
}