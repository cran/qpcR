\name{AICc}
\alias{AICc}

\title{Akaike's second-order corrected Information Criterion for small sample sizes}

\description{
 Calculates the second-order corrected Akaike Information Criterion for objects of class \code{drc}, \code{lm}, \code{glm}, \code{nls}
 or any other models where \code{\link{coefficients}} and \code{\link{residuals}} can be extacted.
 This is a modified version of the original AIC which compensates for bias with small n.
 As qPCR data usually has n/par < 40 (see original reference), AICc was implemented to correct for this.
}

\usage{
  AICc(object)
}

\arguments{
  \item{object}{a fitted model.}
}

\details{
Extends the AIC such that \deqn{AICc = AIC+\frac{2k(k + 1)}{n - k - 1}} with k = number of parameters + 1, and n = number of observations.
For large n, AICc converges to AIC.
}

\value{
The second-order corrected AIC value.
}

\author{
  Andrej-Nikolai Spiess
}

\references{
Sakamoto Y, Ishiguro M, and Kitagawa G (1986). \emph{Akaike Information Criterion Statistics}. D. Reidel Publishing Company.
 
Hurvich CM & Tsai CL (1989). Regression and Time Series Model Selection in Small Samples. \emph{Biometrika} \bold{76}, 297-307.
}

\seealso{
\code{\link{AIC}}, \code{\link{logLik}}.
}

\examples{
m <- pcrfit(reps, 1, 2, l5())
AICc(m)
}

\keyword{models}
\keyword{nonlinear}
