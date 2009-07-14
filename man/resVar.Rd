\name{resVar}
\alias{resVar}

\title{Residual variance of a fitted model}

\description{
Calculates the residual variance for objects of class \code{nls}, \code{lm}, \code{glm}, \code{drc}
 or any other models from which \code{\link{coef}} and \code{\link{residuals}} can be extacted.
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
m <- pcrfit(reps, 1, 2, l5)
resVar(m)
}

\keyword{models}
\keyword{nonlinear}
