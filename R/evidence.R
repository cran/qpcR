evidence <- function(x, y, type = c("AICc", "AIC"))
{
	type <- match.arg(type)

	if (class(x) == "drc" && class(y) == "drc") {
		x1 <- switch(type, AIC = AIC(x), AICc = AICc(x))
		y1 <- switch(type, AIC = AIC(y), AICc = AICc(y))
	}
	else {
		if (is.numeric(x) && is.numeric(y)) {
			x1 <- x
			y1 <- y
		}
		else {
			stop("Input must (both) be either a fitted 'drc' model or numeric!")
		}
	}
	1/(exp(-0.5 * (y1 - x1)))	
}

