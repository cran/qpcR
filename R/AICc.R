AICc <- function(object)
{
	if (class(object)[1] != "drc") stop("AICc works only for 'drc' models!")
	aic <- AIC(object)
	if (!is.numeric(aic)) stop("Cannot calculate AIC!")
	k <- length(object$fit$par) + 1
	n <- object$sumList$lenData
	aic + ((2 * k * (k + 1))/(n - k - 1))
} 