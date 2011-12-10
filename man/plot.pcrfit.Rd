\name{plot.pcrfit}
\alias{plot.pcrfit}

\title{Plotting qPCR data with fitted curves/confidence bands/error bars}

\description{
A plotting function for data of class 'pcrfit' (single curves), 'modlist' (batch curves) or 'replist' (replicate curves) displaying the data points, the fitted curve and (if desired) confidence/prediction bands/error bars on replicates. Four different plot types are available, namely plotting all curves in a 2D graph, a 2D plot matrix, a 3D graph or a heatmap-like image plot.
}

\usage{
\method{plot}{pcrfit}(x, which = c("all", "single", "3D", "image"), 
        fitted = TRUE, add = FALSE, col = NULL, 
        confband = c("none", "confidence", "prediction"),
        errbar = c("none", "sd", "se", "conf"), par3D = list(), 
        par2D = list(), parCI = list(), parSD = list(),  ...)
}

\arguments{
  \item{x}{an object of class 'pcrfit', 'modlist' or 'replist'.}     
  \item{which}{plots all curves in 2D (\code{"all"}), a plot matrix with many curves (\code{"single"}), a 3D plot (\code{"3D"}) or a heatmap-like image plot (\code{image}).} 
  \item{fitted}{should the fitted lines be displayed?}      
  \item{add}{should the curve be added to an existing plot?}
  \item{col}{an optional color vector for the individual curves. Is recycled to the number of runs in \code{x}.} 
  \item{confband}{should confidence/prediction bands be displayed? See \code{\link{confint}}.}
  \item{errbar}{the type of error bar on the plot if replicates exist. See 'Examples'.}
  \item{par3D}{a list containing graphical parameters to change the 3D-plot: \code{\link{plot3d}}, \code{\link{points3d}}, \code{\link{lines3d}}, \code{\link{axis3d}} or \code{\link{mtext3d}}.}
  \item{par2D}{a list containing graphical parameters to change the 2D-plots: \code{\link{plot}}, \code{\link{points}} or \code{\link{lines}}.}
  \item{parCI}{a list containing graphical parameters to change the confidence band: \code{\link{lines}}.}
  \item{parSD}{a list containing graphical parameters to change the error bars: \code{\link{arrows}}.}
  \item{...}{other parameters to be passed to \code{\link{predict}}.}       
}

\value{
A 2D, multiple 2D, 3D or heatmap-like qPCR plot. 
If \code{object} was of class 'replist', colour coding of the curves is done automatically (but can be overridden).  
}

\details{
Uses the 'rgl' package for 3D plots. If the 'modlist' contains runs that failed to fit, these are displayed with RED asterisked names. Additionally, if an outlier method from \code{\link{KOD}} was applied on the 'modlist' with option \code{remove = FALSE}, outlier runs will be displayed in BLUE with double asterisked names. This approach makes the identification of failed runs easy and works only with \code{which = "single"}. See 'Examples'.\cr
For high-throughput data, the user of this function is encouraged to use the \code{"image"} kind of plot, as one can see quite nicely the differences in the amplification profiles of several hundred runs. Of course, this plot type does not display the fitted curve. See 'Examples'.
}
     
\author{
Andrej-Nikolai Spiess
}
        
\examples{
## single plot
m1 <- pcrfit(reps, 1, 2, l5)
plot(m1)

## add another plot in blue
## with 99\% confidence interval
m2 <- pcrfit(reps, 1, 12, l5)
plot(m2, add = TRUE, col = 4, confband = "confidence", level = 0.99)

## plot a 'modlist' batch with coloring of replicates
ml1 <- modlist(reps, 1, 2:13, model = l4)
plot(ml1, col = gl(3,4))   

## subset of data
plot(ml1, type = "n", col = rep(1:3, each = 4), 
     par2D = list(xlim = c(10, 30)))

## plot a 'replist'
rl1 <- replist(ml1, group = gl(3, 4))
plot(rl1)

## standard deviation instead of 
## replicate points; suppress plotting
## of point symbols
plot(rl1, type = "l", errbar = "sd",
     par2D = list(pch = ""))

## 95\% confidence values
plot(rl1, errbar = "conf", 
     par2D = list(pch = ""))

## plot single curves.
## good for diagnostics...
plot(ml1, which = "single", col = rep(1:3, each = 4))

## 3D plots of 'modlist's or 'replist's
plot(ml1, which = "3D", col = rep(1:3, each = 4))
rgl.close()
plot(rl1, which = "3D")
rgl.close()

\dontrun{
## example for "image" type when
## using large data
ml2 <- modlist(vermeulen2)
plot(ml2, which = "image")

## example for outlier identification:
## RED/*name* indicates failed fitting,
## BLUE/**name** indicates sigmoidal outlier
## using 'testdat' set
ml3 <- modlist(testdat, model = l5)
plot(ml3, which = "single") 
}
}

\keyword{models}
\keyword{nonlinear}
