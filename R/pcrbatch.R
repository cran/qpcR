pcrbatch <- function (x, cols = NULL, fct = l4(), group = NULL, type = "cpD2", 
    opt = FALSE, smooth = c("none", "tukey", "lowess"), norm = FALSE, 
    fact = 1, ave = c("mean", "median"), backsub = NULL, plot = FALSE, retPar = FALSE, ...) 
{
    if (class(x) != "modlist" && is.null(cols)) cols <- 2:ncol(x)
    if (class(x) != "modlist" && names(x)[1] != "Cycles") stop("Column 1 should be 'Cycles'!")
    smooth <- match.arg(smooth)
    ave <- match.arg(ave)         
    
    if (!is.null(backsub) && !is.numeric(backsub)) 
        stop("'backsub' must be either NULL or a numeric sequence!")
    if (!is.null(cols) && min(cols) == 1) stop("'cols' must be > 1 because Column 1 must be 'Cycles'!")
    
    resMat <- vector()
    
    if (class(x) == "modlist") {
      wdata <- sapply(x, function(x) x$data[, 2])
      namevec <- sapply(x, function(x) x$names)
      colnames(wdata) <- namevec
      Cycles <- x[[1]]$data[, 1]
    } else {
      wdata <- x[, cols]
      namevec <- colnames(x)[cols]
      colnames(wdata) <- namevec
      Cycles <- x[, 1]
    }
           
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
        wdata <- as.data.frame(data.post)
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
        if (!is.null(backsub)) {
            back <- mean(data[backsub], na.rm = TRUE)
            data <- data - back
        }    
        
        mat <- data.frame(cbind(Cycles, data))
        m <- try(multdrc(data ~ Cycles, data = mat, fct = fct), 
            silent = TRUE)
        
        if (opt) {
            assign("mat", mat, envir = .GlobalEnv)   
            m <- try(mchoice(m, verbose = FALSE, ...), silent = TRUE)
        }       
        
        fctName <-  paste(qpcR:::typeid(m), "()", sep = "")
        
        flush.console()       
        cat("Processing ", namevec[i], "...\n", sep = "") 
                
        cat("   Building sigmoidal model ", fctName, "...\n", sep = "")
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
            outall <- c(ul.eff, fctName, ul.sli, ul.exp)
        }
        else {
            outall <- c(ul.eff, fctName, ul.sli, ul.exp)
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
