pcrbatch <- function (x, cols = 2:ncol(x), group = NULL, model = l4(), type = "cpD2", 
    opt = FALSE, smooth = c("none", "tukey", "lowess"), norm = FALSE, 
    fact = 1, ave = c("mean", "median"), plot = FALSE, retPar = FALSE) 
{
    if (min(cols) == 1) 
        stop("Column 1 should be cycle column!")
    smooth <- match.arg(smooth)
    ave <- match.arg(ave)
    resMat <- vector()
    wdata <- x[, cols]
    namevec <- colnames(x)[cols]
    colnames(wdata) <- namevec
    if (!is.null(group)) {
        group <- as.factor(group)
        if (length(group) != length(cols)) 
            stop("replicates and column numbers do not match!")
        data.pre <- wdata
        if (ave == "mean") 
            centre <- function(x) mean(x, na.rm = TRUE)
        if (ave == "median") 
            centre <- function(x) median(x, na.rm = TRUE)
        data.post <- apply(data.pre, 1, function(x) tapply(x, 
            group, function(x) centre(x)))
        if (nlevels(group) > 1) 
            data.post <- t(data.post)
        wdata <- data.post
        namevec <- paste("group", 1:length(levels(group)), sep = "")
    }
    for (i in 1:ncol(wdata)) {
        data <- wdata[, i] * fact
        if (smooth == "tukey") 
            data <- smooth(data)
        if (smooth == "lowess") 
            data <- lowess(data, f = 0.1)$y
        if (norm == TRUE) 
            data <- data/max(data, na.rm = TRUE)
        mat <- data.frame(cbind(x[, 1], data))
        m <- try(multdrc(data ~ V1, data = mat, fct = model), 
            silent = TRUE)
        if (opt) {
            assign("mat", mat, envir = .GlobalEnv)
            fct = paste(qpcR:::typeid(m), "()", sep = "")
            m <- try(mchoice(m, verbose = FALSE), silent = TRUE)
        }
        fct = paste(qpcR:::typeid(m), "()", sep = "")
        cat("Processing ", namevec[i], "...\n", sep = "")
        cat("   Building sigmoidal model ", fct, "...\n", sep = "")
        out.eff <- try(efficiency(m, plot = plot, type = type)[-12], 
            silent = TRUE)
        if (retPar) 
            out.eff <- c(out.eff, coef(m))
        if (is.null(out.eff$cpE)) 
            out.eff$cpE <- NA
        VAL.eff <- paste("sig.", names(out.eff), sep = "")
        cat("   Using window-of-linearity...\n")
        out.sli <- try(sliwin(m, plot = plot), silent = TRUE)
        VAL.sli <- paste("sli.", names(out.sli), sep = "")
        cat("   Fitting exponential model...\n")
        out.exp <- try(expfit(m, plot = plot)[-c(2, 4, 9)], silent = TRUE)
        VAL.exp <- paste("exp.", names(out.exp), sep = "")
        VALS <- c(VAL.eff, "sig.model", VAL.sli, VAL.exp)
        ul.eff <- as.vector(unlist(out.eff))
        ul.sli <- as.vector(unlist(out.sli))
        ul.exp <- as.vector(unlist(out.exp))
        if (!is.null(group)) {
            outall <- c(ul.eff, qpcR:::typeid(m), ul.sli, ul.exp)
        }
        else {
            outall <- c(ul.eff, qpcR:::typeid(m), ul.sli, ul.exp)
        }
        resMat <- cbind(resMat, outall)
    }
    colnames(resMat) <- namevec
    resMat <- cbind(VALS, resMat)
    cat("Writing to clipboard...\n\n")
    write.table(resMat, file = "clipboard-64000", sep = "\t", 
        row.names = FALSE)
    class(resMat) <- "pcrbatch"
    return(resMat)
}
