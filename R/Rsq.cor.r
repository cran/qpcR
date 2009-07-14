Rsq.cor <- function(object) 
{
  rss <- sum(residuals(object)^2)
  Yi <- residuals(object) - fitted(object)
  cor(Yi, fitted(object))^2
}



