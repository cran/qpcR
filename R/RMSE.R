RMSE <- function(object, which = NULL)
{
	if (is.null(which)) which = 1:length(object$data[, 1]) 
	rmse <- sqrt(mean(residuals(object)[which]^2, na.rm = TRUE))
	return(rmse)
}