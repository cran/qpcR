AICc <- function(object)
{
	aic <- AIC(object)
    	if (!is.numeric(aic)) stop("Cannot calculate AIC!")
    	k <- length(coef(object)) + 1
   	n <- length(residuals(object))
    	aic + ((2 * k * (k + 1))/(n - k - 1))
} 