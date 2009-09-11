\name{fitprob}
\alias{fitprob}      

\title{The (chi-square) fit probability}

\description{
Calculates the fit probability for objects of class \code{drc}, \code{lm}, \code{glm}, \code{nls}
 or any other models from which \code{\link{residuals}} and \code{\link{coef}} can be extracted.   
}

\usage{
fitprob(object)
}

\arguments{
\item{object}{a fitted model.}
}

\details{
Although the residual variance is a measure for the quality of a fit, it provides no information of how good a fit performs in
 relation to other fits in the (hypothetical) infinite population of assays performed under the same conditions. Under the assumption
 that the responses are approximately normally distributed and that the regression's expectation surface is approximately linear in the
 neighbourhood of the best fit, it can be shown that the residual sum-of-squares (RSS) obeys a \eqn{\chi^2} distribution with df = number of
 curve points - number of parameters. The p-value \eqn{\chi^2(RSS, df)} can be viewed as the fraction of an infinite number of assays that, under the same conditions,
 would be expected to have a better curve fit, i.e. a smaller RSS. 
}

\value{
The p-value.
}

\author{
Andrej-Nikolai Spiess
}

\references{
Draper NR, Smith H.\cr
\emph{Applied Regression Analysis, 3rd Ed}.\cr
Wiley, New York, 1998.     
}             

\examples{
## 'good model'
m1 <- pcrfit(reps, 1, 2, l4)
fitprob(m1)

## 'bad model'
m2 <- pcrfit(reps, 1, 2, w4)
fitprob(m1)       
}

\keyword{models}
\keyword{nonlinear}
