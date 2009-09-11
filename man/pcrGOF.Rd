\name{pcrGOF}
\alias{pcrGOF}

\title{Summarize measures for the goodness-of-fit}

\description{
Calculates all implemented measures for the goodness-of-fit and returns them as a list.
Works for objects of class \code{pcrfit}, \code{drc}, \code{lm}, \code{glm}, \code{nls} and many others...
}

\usage{
pcrGOF(object, PRESS = FALSE)
}

\arguments{
  \item{object}{a fitted object.}
  \item{PRESS}{logical. If \code{TRUE}, the more calculation intensive P-square is also returned.}
}

\value{
A list with the following components:
  \item{Rsq}{the R-squared value.}
  \item{Rsq.ad}{the adjusted R-squared value.} 
  \item{AIC}{the Akaike Information Criterion.} 
  \item{AICc}{the bias-corrected Akaike Information Criterion.} 
  \item{BIC}{the Bayesian Information Criteron.}  
  \item{resVar}{the residual variance.}
  \item{RMSE}{the root-mean-squared-error.}
  \item{P.square}{the PRESS P-square, if \code{PRESS = TRUE}.}
}

\author{
Andrej-Nikolai Spiess
}


\examples{
m <- pcrfit(reps, 1, 2, l5)
pcrGOF(m, PRESS = TRUE)
}

\keyword{models}
\keyword{nonlinear}
