\name{resVar}
\alias{resVar}

\title{Residual variance of a sigmoidal fit}

\description{
  Calculates the residual variance of a sigmoidal fit.
}

\usage{
 resVar(object)
}

\arguments{
  \item{object}{an object of class 'drc'.}
 }

\value{
The residual variance from the sigmoidal fit.
}

\author{
  Christian Ritz
}


\examples{
m <- multdrc(F1.1 ~ Cycles, data = reps, fct = l5())
resVar(m)
}

\keyword{models}
\keyword{nonlinear}
