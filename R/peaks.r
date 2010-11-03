peaks <- function (series, span = 3, what = c("max", "min"), do.pad = TRUE, ...) 
{
    if ((span <- as.integer(span))%%2 != 1) 
        stop("'span' must be odd")
    if (!is.numeric(series)) 
        stop("`peaks' needs numeric input")
    what <- match.arg(what)
    if (is.null(dim(series)) || min(dim(series)) == 1) {
        series <- as.numeric(series)
        x <- seq(along = series)
        y <- series
    }
    else if (nrow(series) == 2) {
        x <- series[1, ]
        y <- series[2, ]
    }
    else if (ncol(series) == 2) {
        x <- series[, 1]
        y <- series[, 2]
    }
    if (span == 1) 
        return(list(x = x, y = y, pos = rep(TRUE, length(y))), 
            span = span, what = what, do.pad = do.pad)
    if (what == "min") 
        z <- embed(-y, span)
    else z <- embed(y, span)
    s <- span%/%2
    s1 <- s + 1
    v <- max.col(z, "first") == s1
    if (do.pad) {
        pad <- rep(FALSE, s)
        v <- c(pad, v, pad)
        idx <- v
    }
    else idx <- c(rep(FALSE, s), v)
    val <- list(x = x[idx], y = y[idx], pos = v, span = span, 
        what = what, do.pad = do.pad)
    val
}