modlist <- function (x, cyc = 1, fluo = 2:ncol(x), fct = l5(), opt = FALSE, norm = FALSE)
{
    modList <- NULL
    counter <- 1
    mc <- match.call()
    cL <- as.list(mc)
    if (is.null(cL$fct)) 
        cL$fct <- l5()
    for (i in fluo) {
        Cycles <- x[, cyc]
        Fluo <- x[, i]
        if (norm) Fluo <- Fluo/max(Fluo, na.rm = TRUE)
        Names <- colnames(x[i])
        m <- eval(as.call(list(multdrc, Fluo ~ Cycles, fct = cL$fct)))
        if (opt) {
        	  m <- try(mchoice(m, verbose = FALSE), silent = TRUE)
        }
        cat("Making model for ", names(x)[i], " (", qpcR:::typeid(m), ")\n", sep= "")
        modList[[counter]] <- m
        modList[[counter]]$names <- Names
        counter <- counter + 1
    }
    class(modList) <- "modlist"
    invisible(modList)
}
