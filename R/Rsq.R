Rsq <- function(object) 
{
    rss <- sum(residuals(object)^2)
    mss <- sum((fitted(object) - mean(fitted(object)))^2)
    mss/(mss + rss)
}

