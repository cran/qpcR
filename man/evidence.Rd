\name{evidence}
\alias{evidence}

\title{Evidence ratio for model comparison with AIC or AICc}

\description{
The evidence ratio \deqn{\frac{1}{exp(-0.5 * (AIC2 - AIC1))}} is calculated for either two fitted models or two numerical values.
Models can be compared that are not nested and where the f-test on residual-sum-of-squares is not applicable.
}

\usage{
evidence(x, y, type = c("AICc", "AIC"))
}

\arguments{
  \item{x}{a fitted object or numerical value.}
  \item{y}{a fitted object or numerical value.}
  \item{type}{Bias-corrected (\code{\link{AICc}}, default) or original (\code{\link{AIC}}) versions of the Akaike Information Criterion.} 
}

\details{
Small differences in AIC values can mean substantial more 'likelihood' of one model over the other.
For example, a model with AIC = -130 is nearly 150 times more likely than a model with AIC = -120.   
}

\value{
A value of the first model \code{x} being more likely than the second model \code{y}.
If large, first model is better. If small, second model is better.
}

\author{
Andrej-Nikolai Spiess
}

\examples{
## compare two four-parameter and five-parameter
## log-logistic models
m1 <- pcrfit(reps, 1, 2, l4)
m2 <- pcrfit(reps, 1, 2, l5)
evidence(m2, m1)

## ratio of two AIC's
evidence(-120, -123)
}

\keyword{models}
\keyword{nonlinear}
