typeid <- function (object) 
{
    cl <- class(object)[2]
    npar <- length(object$fct$names)
    mn <- tolower(substr(cl, 1, 1))
    type <- paste(mn, npar, sep = "")
    return(type)
}
