outlier <- function(x, y, pval = 0.05)
{
	res <- vector()
	for (i in 5:length(x))
	    {mod <- lm(y[1:i] ~ x[1:i], na.action = na.exclude)
		st <- studtest(mod)
          st1 <- tail(st, 1)
		pst1 <- 1 - pt(st1, df = mod$df.residual)
		res <- c(res, pst1)					
	}
	resl <- sapply(res, function(x) x < pval)	
	resl[is.na(resl)] <- FALSE
	for (i in 1:length(resl - 2)) {
		if (resl[i] == TRUE && resl[i + 1] == TRUE &&  resl[i + 2] == TRUE)
		    return(which = as.numeric(names(resl[i])))
		}
}
