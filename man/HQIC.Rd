\name{HQIC}
\alias{HQIC}

\title{Hannan-Quinn Information Criterion}

\description{
Calculates the Hannan-Quinn Information Criterion for objects of class \code{pcrfit}, \code{nls}, \code{lm}, \code{glm} or any other models from which \code{\link{logLik}}, \code{\link{coef}} and \code{\link{residuals}} can be extracted. It is somewhat similar to \code{\link{BIC}}, but penalizes \eqn{n} even more by double logarithmation.
}

\usage{
HQIC(object)
}

\arguments{
\item{object}{a fitted model.}
}

\details{
\deqn{HQIC = -2 \cdot log(\mathcal{L}_{max}) + 2 \cdot k \cdot log(log(n))}
with \eqn{\mathcal{L}_{max}} = maximum likelihood, \eqn{k} = number of parameters and \eqn{n} = number of observations. 
}

\value{
The Hannan-Quinn Information Criterion.
}

\author{
Andrej-Nikolai Spiess
}

\references{
The Determination of the Order of an Autoregression.\cr
Hannan EJ & Quinn BG.\cr 
\emph{J Roy Stat Soc B} (1979), \bold{41}: 190-195.\cr
}

\note{
For \eqn{n > \sim 20} or so \code{\link{BIC}} is the strictest in penalizing loss of degree of freedom by having more parameters in the fitted model. For \eqn{n > \sim 40} \code{\link{AIC}} is the least strict of the three and \code{\link{HQIC}} holds the middle ground, or is the least penalizing for \eqn{n < \sim 20}.
}

\seealso{
\code{\link{AIC}}, \code{\link{BIC}}.
}

\examples{
m1 <- pcrfit(reps, 1, 2, l5)
HQIC(m1)
}

\keyword{models}
\keyword{nonlinear}
