\name{RMSE}
\alias{RMSE}

\title{Root-mean-squared-error of a fitted model}

\description{
  Calculates the root-mean-squared-error (RMSE) for objects of class \code{drc}, \code{lm}, \code{glm}, \code{nls}
 or any other models where \code{\link{residuals}} can be extacted. 
}

\usage{
  RMSE(object, which = NULL)
}

\arguments{
  \item{object}{a fitted model.}
  \item{which}{the part of the curve to be used for RMSE calculation. If not defined, the complete curve is used.}
 }

\value{
The root-mean-squared-error from the fit or a part thereof.
}

\author{
  Andrej-Nikolai Spiess
}


\examples{
### for a part of the curve
m <- multdrc(F1.1 ~ Cycles, data = reps, fct = l5())
RMSE(m, 10:15)
}

\keyword{models}
\keyword{nonlinear}
