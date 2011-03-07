\name{resplot}
\alias{resplot}

\title{A residuals bar-plot with colour-coded bars}

\description{
 A plotting function which displays a barplot of the residuals of any object from which \code{\link{residuals}}
 can be extracted. The bars are colour-coded with heat colours proportional to the residual value.
 As default, the residuals are displayed together with the points of the amplification curve.
}

\usage{
resplot(object, overlay = TRUE, ...)
}

\arguments{
 \item{object}{either residuals (numeric) or any fitted object, i.e. of class \code{nls}, \code{lm}, \code{drc} etc.} 
 \item{overlay}{logical. If \code{TRUE}, the residuals are plotted on top of the amplification curve, else alone.}
 \item{...}{any other parameters to be passed to \code{\link{barplot}}.} 	
}

\value{
A plot as described above.
}

\author{
Andrej-Nikolai Spiess
}

\examples{
## create l5 model and plot residuals
m1 <- pcrfit(reps, 1, 2, l5)
resplot(m1)
}

\keyword{models}
\keyword{nonlinear}
