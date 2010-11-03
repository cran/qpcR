xyy.plot <- function(x, y1, y2, y1.par = NULL, y2.par = NULL, first = NULL, ...)
{
  if (is.null(y1.par)) y1.par <- list()
  if (is.null(y2.par)) y2.par <- list()   
  par(mar = c(5, 4, 4, 5))   
  if (!is.null(first)) eval(first)
  suppressWarnings(do.call(plot, modifyList(list(x = x, y = y1, xlab = deparse(substitute(x)), 
         ylab = deparse(substitute(y1)), col = 1, ...), y1.par)))
  par(new = TRUE)
  suppressWarnings(do.call(plot, modifyList(list(x = x, y = y2, axes = FALSE, xlab = "", 
         ylab = "", col = 2), y2.par))) 
  suppressWarnings(do.call(axis, modifyList(list(side = 4, col = 2, col.ticks = 2, col.axis = 2), y2.par)))
  suppressWarnings(do.call(mtext, modifyList(list(text = deparse(substitute(y2)), 
                   side = 4, line = 3, col = 2), y2.par))) 
}