ratiocalc <- function(data, group = NULL, ratio.fun = c("E1-E2", "E1-E1", "2-2"), which.eff = c("sig", "sli", "exp"), ...)
{
      if (class(data) == "modlist") {
            ratiocalc.modlist(data = data, group = group, ratio.fun = ratio.fun, ...)
      } else {
            if (class(data) == "pcrbatch") {
                  ratiocalc.pcrbatch(data = data, group = group, ratio.fun = ratio.fun, which.eff = which.eff, ...)
            } else {
                  stop("Data must be either of class 'modlist' or 'pcrbatch'")
            }
      }
}
