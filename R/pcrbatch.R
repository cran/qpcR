pcrbatch <- function(
x, 
cols = 2:ncol(x), 
group = NULL,
model = l4(),
type = "cpD2",
opt = TRUE,
smooth = "tukey",
norm = FALSE,
fact = 1,
ave = "mean",
plot = FALSE
)
{
	if (min(cols) == 1) stop("Column 1 should be cycle column!")

	resMat <- vector()
	wdata <- x[, c(1, cols)]
    	namevec <- names(wdata)

    	if (!is.null(group)) {
        	group <- as.factor(group)
        	if (length(group) != length(cols)) stop("replicates and column numbers do not match!")
        	data.pre <- wdata[, cols]
        	if (ave == "mean") centre <- function(x) mean(x, na.rm = TRUE)
        	if (ave == "median") centre <- function(x) median(x, na.rm = TRUE)
        	data.post <- apply(data.pre, 1, function(x) tapply(x, 
            			group, function(x) centre(x)))
        	if (nlevels(group) > 1) data.post <- t(data.post)
        	data.post <- cbind(wdata[, 1], data.post)
        	wdata <- data.post
        	namevec <- c(NA, paste("group", 1:length(levels(group)), sep = ""))
    	}

	for (i in 2:ncol(wdata)) {
        	data <- x[, i] * fact
        	if (smooth == "tukey") data <- smooth(data)
        	if (smooth == "lowess") data <- lowess(data, f = 0.1)$y
        	if (norm == TRUE) data <- data/max(data, na.rm = TRUE)
        	mat <- data.frame(cbind(x[, 1], data))
        	m <- try(multdrc(data ~ V1, data = mat, fct = model), silent = TRUE)
        	if (opt) {
            		assign("mat", mat, envir = .GlobalEnv)
    		  	fct = paste(qpcR:::typeid(m), "()", sep = "")
            		m <- try(mchoice(m, verbose = FALSE), silent = TRUE)
        	}
	   	fct = paste(qpcR:::typeid(m), "()", sep = "")
        	cat("Processing ", namevec[i], "...\n", sep = "")
	   	cat("   Building sigmoidal model ", fct, "...\n", sep = "")
        	out.eff <- try(efficiency(m, plot = plot, type = type), silent = TRUE)
        	if (is.null(out.eff$cpE)) out.eff$cpE <- NA
        	VAL.eff <- paste("sig.", names(out.eff), sep = "")
        	cat("   Using window-of-linearity...\n")
        	out.sli <- try(sliwin(m, plot = plot), silent = TRUE)
        	VAL.sli <- paste("sli.", names(out.sli), sep = "")
        	cat("   Fitting exponential model with resVar...\n")
        	out.exp1 <- try(expfit(m, plot = plot)[-c(3, 8)], silent = TRUE)
        	out.exp1$cyc.best <- out.exp1$cyc.best[1]
        	VAL.exp1 <- paste("exp1.", names(out.exp1), sep = "")
        	cat("   Fitting exponential model with studentized residuals...\n")
        	out.exp2 <- try(expfit2(m, plot = plot)[-c(3, 8)], silent = TRUE)
        	VAL.exp2 <- paste("exp2.", names(out.exp2), sep = "")
        	VALS <- c(VAL.eff, "sig.model", VAL.sli, VAL.exp1, VAL.exp2)
        	ul.eff <- as.vector(unlist(out.eff))
       		ul.sli <- as.vector(unlist(out.sli))
        	ul.exp1 <- as.vector(unlist(out.exp1))
        	ul.exp2 <- as.vector(unlist(out.exp2))
        	if (!is.null(group)) {
			outall <- c(ul.eff, qpcR:::typeid(m), ul.sli, ul.exp1, ul.exp2)
        	}
        	else {
            		outall <- c(ul.eff, qpcR:::typeid(m), ul.sli, ul.exp1, ul.exp2)
        	}
        	resMat <- cbind(resMat, outall)
    	}

    	colnames(resMat) <- namevec[-1]
    	resMat <- cbind(VALS, resMat)
    	cat("Writing to clipboard...\n\n")
    	write.table(resMat, file = "clipboard-64000", sep = "\t", row.names = FALSE)
    	invisible(resMat)
}

