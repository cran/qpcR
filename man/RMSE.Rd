\name{RMSE}
\alias{RMSE}

\title{Root-mean-squared-error of a sigmoidal fit}

\description{
  Calculates the root-mean-squared-error (RMSE) of the fit or a region of the fit.
}

\usage{
  RMSE(object, which = NULL)
}

\arguments{
  \item{object}{an object of class 'drc'.}
  \item{which}{the part of the curve to be used for RMSE calculation. If not defined, the complete curve is used.}
 }

\value{
The root-mean-squared-error from the sigmoidal fit or a part thereof.
}

\author{
  Christian Ritz
}


\examples{
### for a part of the curve
m <- multdrc(F1.1 ~ Cycles, data = reps, fct = l5())
RMSE(m, 10:15)
}

\keyword{models}
\keyword{nonlinear}
