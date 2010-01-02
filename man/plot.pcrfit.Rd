\name{plot.pcrfit}
\alias{plot.pcrfit}

\title{Plotting qPCR data with fitted curves/confidence bands/error bars}

\description{
A plotting function for data of class 'pcrfit' (single curves), 'modlist' (batch curves) or 'replist' (replicate curves) displaying the
 data points, the fitted curve and (if desired) confidence/prediction bands/error bars on replicates. Three different plot types are available,
  namely plotting all curves in 2D/3D graphs or a plot matrix with one graph per curve.
}

\usage{
\method{plot}{pcrfit}(x, fitted = TRUE,  
        confband = c("none", "confidence", "prediction"),
        errbar = c("none", "sd", "se", "conf"), add = FALSE, 
        col = NULL, level = 0.95, xlim = NULL, 
        which = c("all", "single", "3D"), ...)
}

\arguments{
  \item{x}{an object of class 'pcrfit', 'modlist' or 'replist'.}     
  \item{fitted}{should the fitted line be displayed?}      
  \item{confband}{should confidence/prediction bands be displayed? See \code{\link{confint}}.}
  \item{errbar}{the type of error bar on the plot if replicates exist. See 'Examples'.}
  \item{add}{should the curve be added to an existing plot?}
  \item{col}{an optional color vector for the individual curves. Is recycled to the number of runs in \code{x}.} 
  \item{level}{the confidence level used for \code{confband} or \code{errbar}.}
  \item{xlim}{a 2-element vector of the x-value (Cycles) region that should be plotted.}
  \item{which}{plots \code{all} curves in one 2D or \code{3D} graph or many graphs with \code{single} curves.} 
  \item{...}{other parameters to be passed to \code{\link{plot}} or \code{\link{predict}}.}       
}

\value{
A qPCR plot. 
If \code{object} was of class 'replist', colour coding of the curves is done automatically.  
}

\details{
Uses the 'rgl' package for 3D plots. 
If the 'modlist' contains data that failed to fit, these are missing in the plots.
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
ml <- modlist(reps, 1, 2:13, model = l4)
plot(ml, col = gl(3,4))   

## only the fitted curves
## and a subset of data
plot(ml, type = "n", col = gl(3,4), xlim = c(10, 30))

## plot a 'replist'
rl <- replist(ml, group = gl(3, 4))
plot(rl)

## standard deviation instead of 
## replicate points
plot(rl, type = "n", errbar = "sd")
## 95\% confidence values
plot(rl, type = "n", errbar = "conf")

## plot single curves.
## good for diagnostics...
plot(ml, which = "single", col = gl(3, 4))

## 3D plots of 'modlist's or 'replist's
plot(ml, which = "3D", col = gl(3, 4))
rgl.close()
plot(rl, which = "3D")
rgl.close()
}

\keyword{models}
\keyword{nonlinear}
