pcrpred <- function(object, newdata, which = c("y", "x"), ...)
{
	which <- match.arg(which)

	if (which == "y") {
		pred <- sapply(newdata, function(x) predict(object, newdata = data.frame(x), ...)[1])
	}

     	if (which == "x") {
		if (is.null(object$fct$inversion)) stop("no inverse function available for this model!")
		pred <- object$fct$inversion(newdata, object$parmMat)
	}
		
	pred <- as.numeric(pred)
	return(pred)
}
