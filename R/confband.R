confband <- function (object, level = 0.95) 
{
    x <- unique(object$data[, 1])
    all <- predict(object, interval = "confidence", level = level)
    clo <- all[, 3]
    cup <- all[, 4]
    invisible(list(x = x, clo = clo, cup = cup))
}
