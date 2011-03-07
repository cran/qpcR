\name{akaike.weights}
\alias{akaike.weights}

\encoding{latin1}

\title{Calculation of Akaike weights/relative likelihoods/delta-AICs}

\description{
Calculates Akaike weights from a vector of AIC values. 
}

\usage{
akaike.weights(x)
}

\arguments{
\item{x}{a vector containing the AIC values.}
}

\details{
Although Akaike's Information Criterion is recognized as a major measure for selecting models, it has one major drawback:
The AIC values lack intuitivity despite higher values meaning less goodness-of-fit. For this
purpose, Akaike weights come to hand for calculating the weights in a regime of several models. Additional measures can be derived,
such as delta-AIC's and relative likelihoods that demonstrate the probability of one model being in favor over the other. 
This is done by using the following formulae:


delta AICs: \deqn{\Delta_i(AIC) = AIC_i - min(AIC)}

relative likelihood: \deqn{L \propto exp\left(-\frac{1}{2}\Delta_i(AIC)\right)}

Akaike weights: \deqn{w_i(AIC) = \frac{exp\left(-\frac{1}{2}\Delta_i(AIC)\right)}{\sum_{k=1}^K exp\left(-\frac{1}{2}\Delta_k(AIC)\right)}}
}

\value{
A list containing the following items:
\item{deltaAIC}{the delta-AIC values.}
\item{rel.LL}{the relative likelihoods.}
\item{weights}{the Akaike weights.}
}

\author{
Andrej-Nikolai Spiess
}

\references{
Classical literature:\cr
Sakamoto Y, Ishiguro M, and Kitagawa G (1986). Akaike Information Criterion Statistics. D. Reidel Publishing Company.\cr
Burnham KP, Anderson DR. Model selection and inference: a practical information-theoretic approach (2002). Springer Verlag, New York, USA
 
A good summary:\cr
Wagenmakers EJ, Farrell Simon. AIC model selection using Akaike weights (2004). \emph{Psychonomic Bull Review}, \bold{11}: 192-196.
}

\seealso{
\code{\link{AIC}}, \code{\link{logLik}}.
}

\examples{
## apply a list of different sigmoidal models to data
## and analyze GOF statistics with Akaike weights
## on 6 different sigmoidal models 
modList <- list(l5, l4, l3, b5, b4, b3)
aics <- sapply(modList, function(x) AIC(pcrfit(reps, 1, 2, x))) 
akaike.weights(aics)$weights 
}

\keyword{models}
\keyword{nonlinear}
