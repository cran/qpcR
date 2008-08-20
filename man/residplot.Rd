\name{residplot}
\alias{residplot}

\title{A residuals bar-plot with colour-coded bars}

\description{
 A simple plotting function which displays a barplot of the residuals of any object from which \code{\link{residuals}}
 can be extracted. The bars are colour-coded with heat colours proportional to the residual value.  
}

\usage{
 residplot(object, ...)
}

\arguments{
 \item{object}{either residuals (numeric) or any fitted object, i.e. of class \code{lm}, \code{drc}, \code{nls} etc.} 
 \item{...}{any other parameters to be passed to \code{\link{barplot}}.} 	
}

\value{
 A plot as described above.  
}

\author{
  Andrej-Nikolai Spiess
}

\examples{
### create l5 model and plot residuals
m1 <- pcrfit(reps, 1, 2, l5())
residplot(m1)
### compare to l4 model (added smaller bars)
m2 <- pcrfit(reps, 1, 2, l4())
residplot(m2, add = TRUE, width = 0.5, space = c(1.4, 0.9))
}

\keyword{models}
\keyword{nonlinear}
