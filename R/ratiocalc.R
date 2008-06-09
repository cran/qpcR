ratiocalc <- function(data, group = NULL, ratio = c("ind", "first"), which.eff = c("sig", "sli", "exp"),
                        iter = c("combs", "perms"), rep.all = TRUE, ttest = c("cp", "Ecp"), ...)
{
      if (class(data) == "modlist") {
            ratiocalc.modlist(data = data, group = group, ratio = ratio, iter = iter, rep.all = rep.all, ttest = ttest, ...)
      } else {
            if (class(data) == "pcrbatch") {
                  ratiocalc.pcrbatch(data = data, group = group, ratio = ratio, which.eff = which.eff, iter = iter, rep.all = rep.all, ttest = ttest, ...)
            } else {
                  stop("Data must be either of class 'modlist' or 'pcrbatch'")
            }
      }
}
