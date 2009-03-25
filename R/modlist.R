modlist <- function (x, cols = NULL, fct = l4(), opt = FALSE, norm = FALSE, backsub = NULL, ...)
{
    modList <- NULL
    counter <- 1
    
    if (is.null(cols)) cols <- 2:ncol(x)        
    
    if (!is.null(backsub) && !is.numeric(backsub)) 
        stop("'backsub' must be either NULL or a numeric sequence!")
    
    if (!is.null(cols) && min(cols) == 1) stop("'cols' must be > 1 because Column 1 must be 'Cycles'!")
    
    Cycles <- x[, 1]
    NAMES <- colnames(x[cols])
           
    for (i in cols) {
        data  <- x[, i]
               
        if (norm) data <- data/max(data, na.rm = TRUE)
        
        if (!is.null(backsub)) {
            back <- mean(data[backsub], na.rm = TRUE)
            data <- data - back
        }           
        
        m <- eval(as.call(list(drmfit, data ~ Cycles, fct = fct)))
        
        if (opt) {
        	  m <- try(mchoice(m, verbose = FALSE, ...), silent = TRUE)
        }
        
        flush.console()
        cat("Making model for ", NAMES[counter], " (", qpcR:::typeid(m), ")\n", sep= "")
        modList[[counter]] <- m
        modList[[counter]]$names <- NAMES[counter]
        counter <- counter + 1  
    }
    
    class(modList) <- "modlist"
    invisible(modList)
}
