\name{Rsq}
\alias{Rsq}

\title{R-square of a sigmoidal fit}

\description{
  Calculates the R-square value of a sigmoidal fit.
}

\usage{
  Rsq(object)
}

\arguments{
  \item{object}{an object of class 'drc'.}
 }

\value{
The R-square value from the sigmoidal fit.
}

\author{
  Christian Ritz
}


\examples{
m <- multdrc(F1.1 ~ Cycles, data = reps, fct = l5())
Rsq(m)
}

\keyword{models}
\keyword{nonlinear}
