confband <- function (object, level = 0.95) 
{
    x <- unique(object$data[, 1])
    all <- as.data.frame(predict(object, interval = "confidence", level = level))
    clo <- all$Lower
    cup <- all$Upper
    invisible(list(x = x, clo = clo, cup = cup))
}

