studtest <- function(model)
{	
	require(MASS, quietly = TRUE)
	sr <- studres(model)
	return(sr)	
}
