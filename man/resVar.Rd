\name{resVar}
\alias{resVar}

\title{Residual variance of a fitted model}

\description{
  Calculates the residual variance for objects of class \code{drc}, \code{lm}, \code{glm}, \code{nls}
 or any other models where \code{\link{coef}} and \code{\link{residuals}} can be extacted.
}

\usage{
 resVar(object)
}

\arguments{
  \item{object}{a fitted model.}
 }

\value{
The residual variance of the fit.
}

\author{
  Andrej-Nikolai Spiess
}


\examples{
m <- multdrc(F1.1 ~ Cycles, data = reps, fct = l5())
resVar(m)
}

\keyword{models}
\keyword{nonlinear}
