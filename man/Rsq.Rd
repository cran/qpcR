\name{Rsq}
\alias{Rsq}

\title{R-square value of a fitted model}

\description{
  Calculates the R-square value for objects of class \code{drc}, \code{lm}, \code{glm}, \code{nls}
 or any other models where \code{\link{fitted}} and \code{\link{residuals}} can be extacted.
}

\usage{
  Rsq(object)
}

\arguments{
  \item{object}{a fitted model.}
 }

\value{
The R-square value of the fit.
}

\author{
  Andrej-Nikolai Spiess
}


\examples{
m <- multdrc(F1.1 ~ Cycles, data = reps, fct = l5())
Rsq(m)
}

\keyword{models}
\keyword{nonlinear}
