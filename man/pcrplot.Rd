\name{pcrplot}
\alias{pcrplot}

\title{Plotting PCR data with fitted curves}

\description{
  Displays fitted curves and observations in the same plot window, distinguishing between curves by different plot 
  symbols and line types. Various additional features can be plotted, i.e. error bars or confidence bands.
}

\usage{
 pcrplot(x, ..., level = NULL, col = FALSE, 
       conLevel, conName, grid = 100, legend, legendText, legendPos,
       type = "average", lty, log = "", pch, xlab, ylab, xlim, ylim, 
       xt = NULL, xtlab = NULL, yt = NULL, ytlab = NULL,
       add = FALSE, confband = NULL)
}

\arguments{
  \item{x}{an object of class 'drc'.}
  \item{...}{additional arguments.}
  \item{level}{vector of character strings. To plot only the curves
          specified by their names.}
  \item{col}{either logical or a vector of colours. If TRUE default
          colours are used. If FALSE (default) no colours are used.}
  \item{conLevel}{numeric. Dose level below which the dose is zero (the amount
          of stretching on the x-axis above zero). Default is 1e-2.}
  \item{conName}{character string. Name on x axis for dose zero. Default is '"0"'.}
  \item{grid}{ numeric. Number of points in the grid used for plotting the fitted curves.}
  \item{legend}{logical. If TRUE a legend is displayed.}
  \item{legendText}{a character string or vector of character strings specifying the legend text.}
  \item{legendPos}{numeric vector of length 2 giving the position of the legend.}
  \item{type}{a character string specifying how the originals observations
          should be plotted. There are 5 options:  "average" (default),
          "none" (only the fitted curve(s)), "obs" (only the data
          points), "all" or "points"  (all data points), "errbar" (errorbars +- stdev).}
  \item{lty}{a numeric vector specifying the line types.}
  \item{log}{a character string which contains '"x"' if the x axis is to
          be logarithmic, '"y"' if the y axis is to be logarithmic and
          '"xy"' or '"yx"' if both axes are to be logarithmic. The
          default is "", which yields the original axes.}
  \item{pch}{a vector of plotting characters or symbols.}
  \item{xlab}{an optional label for the x axis.}
  \item{ylab}{an optional label for the y axis.}
  \item{xlim}{a numeric vector of length two, containing the lower and
          upper limit for the x axis.}
  \item{ylim}{a numeric vector of length two, containing the lower and
          upper limit for the y axis.}
  \item{xt}{a numeric vector containing the positions of the tick marks
          on the x axis.}
       \item{xtlab}{a vector containing the tick marks on the x axis.}
  \item{yt}{a numeric vector containing the positions of the tick marks
          on the y axis.}
  \item{ytlab}{a vector containing the tick marks on the y axis.}
  \item{add}{logical. If TRUE then add to already existing plot.} 
  \item{confband}{numeric. The confidence level.} 
}

\value{
An invisible data frame with the values used for plotting the fitted curves. 
The first column contains the PCR cycles, and the following columns (one for each curve) contain
the fitted fluorescence data.
}

\details{
Is an extension of the \code{plot.drc} function from the 'drc' package.
}

\author{
  Christian Ritz & Andrej-Nikolai Spiess
}


\examples{
repData <- repform(reps[,1:5], c(0,1,1,1,1))
m <- multdrc(values ~ Cycles, curve = Curve, data = repData, fct = l5())
###plot the average
pcrplot(m)
###all observations
pcrplot(m, type = "all")
###mean with errorbars
pcrplot(m, type = "errbar")
###plot one curve and add 95 percent confidence band
m2 <- multdrc(F1.1 ~ Cycles, data = reps, fct = l5())
pcrplot(m2, confband = 0.95)
}

\keyword{models}
\keyword{nonlinear}
