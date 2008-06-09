ratioplot <- function(data, errbar = c("prop", "propSim", "eval", "evalSim", "none"), sem = FALSE,
                        order = NA, y.fac = 1.5, normcol = NA, errwid = NULL, type = c("bar", "hbar"),
                        plot.t = c("none", "stars", "values"), offset = 1, ...)
{
      mc <- match.call(expand.dots = FALSE)$...

      if(class(data) != "ratiocalc") stop("Supplied data is not a result from 'ratiocalc'")
      
      errbar <- match.arg(errbar)
      type <- match.arg(type)
      plot.t <- match.arg(plot.t)
      
      if (!is.logical(sem)) stop("'sem' must be a logical")
      
      DATA <- data$ratio

      if (!all(is.na(order))) DATA <- DATA[, order]
      
      MEANDATA <- DATA[1, ]
      ERRDATA <- switch(errbar, none = rep(0, ncol(DATA)), prop = DATA[3, ], propSim = DATA[4, ], eval = DATA[5, ], evalSim = DATA[6, ])
      
      if (sem) ERRDATA <- ERRDATA/sqrt(DATA[7, ])
      
      if(!is.na(normcol)) {
            normfac <- 1/MEANDATA[normcol]
            MEANDATA <- normfac * MEANDATA
            ERRDATA <- normfac * ERRDATA
      }
      
      ERRDATA[is.na(ERRDATA)] <- 0
      
      if (type == "hbar") {
            m <- which(MEANDATA < 1)
            OLDMEANDATA <- MEANDATA[m]
            OLDERRDATA <- ERRDATA[m]
            MEANDATA[m] <- -1/MEANDATA[m]
            ERRDATA[m] <- MEANDATA[m] * (OLDERRDATA / OLDMEANDATA)
      }

      LIMDATA <- MEANDATA + ERRDATA

      if (type == "bar") MINLIMDATA <- 0 else MINLIMDATA <- 2 * y.fac * min(LIMDATA)
      MAXLIMDATA <- y.fac * max(LIMDATA)
      YLIM <- c(MINLIMDATA, MAXLIMDATA)

      if (any(names(mc) == "ylim")) barplot(LIMDATA, col = 0, border = 0, axisnames = FALSE, ...)
            else barplot(LIMDATA, col = 0, border = 0, axisnames = FALSE, ylim = YLIM, ...)
      
      ABSCISSA <- barplot(MEANDATA, names.arg = colnames(DATA), xlab = "Sample ratios",
                         ylab = "Ratio", main = "Ratios of all sample combinations", add = TRUE, ...)
                         
      if(is.null(errwid)) errwid <- 0.05 * (ABSCISSA[2] - ABSCISSA[1])
      if (length(ABSCISSA) == 1) errwid <- 0.5
      
      ERRDATA[ERRDATA == 0] <- NA
      
      if (errbar != "none") arrows(ABSCISSA, MEANDATA, ABSCISSA, MEANDATA + ERRDATA, angle = 90, lwd = 1.5, length = errwid, ...)
      
      if (plot.t == "stars") {
            PVALS <- DATA[8, ]
            pvals <- PVALS
            pvals[pvals < 0.001] <- "***"
            pvals[pvals < 0.01] <- "**"
            pvals[pvals < 0.05] <- "*"
            pvals[pvals > 0.05] <- ""
            OFFSET <- rep(offset, length(ERRDATA))
            if (type == "hbar") OFFSET[which(MEANDATA < 1)] <- -offset
            points(ABSCISSA, MEANDATA + ERRDATA + OFFSET, pch = pvals)
      }
      if (plot.t == "values") {
            PVALS = DATA[8, ]
            OFFSET <- rep(offset, length(ERRDATA))
            if (type == "hbar") OFFSET[which(MEANDATA < 1)] <- -offset
            text(ABSCISSA, MEANDATA + ERRDATA + OFFSET, round(PVALS, 5), cex = 0.5)
      }
      
      
}