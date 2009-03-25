pcrfit <- function (data, cyc = 1, fluo, fct = l5()) 
{
    mc <- match.call()
    cL <- as.list(mc)
    if (is.null(cL$fct)) cL$fct <- l5()
    Cycles <- data[, cyc]
    Fluo <- data[, fluo]
    DATA <- as.data.frame(cbind(Cycles, Fluo))
    mod <- eval(as.call(list(drmfit, Fluo ~ Cycles, data = DATA, fct = cL$fct)))
    invisible(mod)
}
