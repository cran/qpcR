confband <- function (object, level = 0.95) 
{
    x <- unique(object$data[, 1])
    all <- predict(object, interval = "confidence", level = level)
    clo <- all[, 2]
    cup <- all[, 3]
    invisible(list(x = x, clo = clo, cup = cup))
}
