repform <- function(x, group)
{
	x <- as.data.frame(x)
	Curve <- vector()
	if (class(group) != "factor") {
         group <- as.factor(group)
     }
	if (length(group) != ncol(x)) {
         stop("grouping and column number do not match")
     }
     for (i in 1:length(levels(group))) {
          lg <- levels(group)[i]
          if (lg == "0") {
              cyc.col <- which(group == lg)
              next
          }
          pos <- which(group == lg)
          rep <- length(pos)
          Curve <- c(Curve, rep(lg, rep * length(x[, cyc.col])))  
          Curve <- as.factor(Curve)   
     }
     ret <- data.frame(Curve = Curve, Cycles = rep(x[, cyc.col],
                       sum(group != 0)),stack(x[, which(group != 0)]))
     return(ret)
}
