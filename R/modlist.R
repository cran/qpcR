modlist <- function(x, cyc = 1, fluo = 2:ncol(x), fct = l5())
{
      modList <- NULL
      counter <- 1

      mc <- match.call()
      cL <- as.list(mc)
      if (is.null(cL$fct)) cL$fct <- l5()

      for (i in fluo) {
            cat("Making model for", names(x)[i], "\n")
            Cycles <- x[, cyc]
            Fluo <- x[, i]
            modList[[counter]] <- eval(as.call(list(multdrc, Fluo ~ Cycles, fct = cL$fct)))
            counter <- counter + 1
      }
      invisible(modList)
}