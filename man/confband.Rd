\name{confband}
\alias{confband}

\title{Confidence bands for a sigmoidal fit}

\description{
  Calculates the confidence band for a fitted 'drc' model.
}

\usage{
  confband(object, level = 0.95)
}

\arguments{
  \item{object}{an object of class 'drc'.}
  \item{level}{the confidence level, defaults to 0.95.}
}

\value{
A list with the following components:
  \item{x}{the x values, i.e. cycle number.}
  \item{clo}{the response values of the lower confidence limit.}
  \item{cup}{the response values of the upper confidence limit.}  
}

\details{
Calls the \code{predict.drc} function from the 'drc' package.
The confidence intervals are based on asymptotic normality.
}

\author{
  Andrej-Nikolai Spiess & Christian Ritz
}

\seealso{
Also implemented in the \code{\link{pcrplot}} function for plotting the confidence bands.
}

\examples{
m <- multdrc(F1.1 ~ Cycles, data = reps, fct = l5())
###with 99 percent confidence interval
ci <- confband(m, level = 0.99)
print(ci)
###plotting the confidence bands
pcrplot(m, confband = 0.99)
}

\keyword{models}
\keyword{nonlinear}
