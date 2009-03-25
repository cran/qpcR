\name{RSS}
\alias{RSS}

\title{Residual sum-of-squares of a fitted model}

\description{
  Calculates the residual sum-of-squares for objects of class \code{drc}, \code{lm}, \code{glm}, \code{nls}
 or any other models where \code{\link{residuals}} can be extacted. 
}

\usage{
  RSS(object)
}

\arguments{
  \item{object}{a fitted model.}  
}

\value{
The residual sum-of-squares from the fit.
}

\author{
  Andrej-Nikolai Spiess
}


\examples{
m <- pcrfit(reps, 1, 2, l5())
RSS(m)
}

\keyword{models}
\keyword{nonlinear}
