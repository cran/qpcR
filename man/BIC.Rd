\name{BIC}
\alias{BIC}      

\title{Bayesian Information Criterion}

\description{
 Calculates the Bayesian Information Criterion for objects of class \code{drc}, \code{lm}, \code{glm}, \code{nls}
 or any other models from which \code{\link{logLik}}, \code{\link{coef}} and \code{\link{residuals}} can be extacted.   
}

\usage{
BIC(object)
}

\arguments{
    \item{object}{a fitted model.}
}

\details{
\deqn{BIC = - 2 * logLik(object) + npar * log(nobs)}
with npar = number of parameters, nobs = number of observations.
}

\value{
The BIC value.
}

\author{
Andrej-Nikolai Spiess
}

\references{
Schwartz GE (1978).\cr
Estimating the dimension of a model.\cr
\emph{Annals of Statistics}, \bold{6}: 461-464.
}

\seealso{
\code{\link{AIC}}, \code{\link{logLik}}.
}

\examples{
m <- pcrfit(reps, 1, 2, l5)
BIC(m)
}

\keyword{models}
\keyword{nonlinear}
